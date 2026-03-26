"""
5' End Aggregation Module for RECTIFY.

This module clusters reads by their 5' end positions (TSS) and
computes gene attribution based on the 3' ends of reads.

The output is a dataset suitable for:
- Transcription start site (TSS) analysis
- Promoter usage studies
- 5' capping / transcription initiation mapping

For splice-aware 5' end correction, use the corrected 5' positions
from the splice_aware_5prime module.

Author: Kevin R. Roy
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from ..analyze.gene_attribution import (
    get_read_3prime_position,
    get_read_5prime_position,
    get_read_body_interval,
    build_cds_interval_tree,
    compute_read_gene_attribution,
)

logger = logging.getLogger(__name__)


@dataclass
class FivePrimeCluster:
    """A cluster of reads sharing a 5' end region."""
    chrom: str
    strand: str
    position: int  # Representative position (mode or median)
    start: int  # Cluster start (leftmost 5' end)
    end: int  # Cluster end (rightmost 5' end)
    n_reads: int
    read_ids: List[str] = field(default_factory=list)

    # Attribution from 3' ends
    attribution_string: str = ""
    primary_gene: str = ""
    primary_fraction: float = 0.0
    n_intergenic: int = 0

    # Additional metrics
    mean_3prime_position: Optional[float] = None
    median_3prime_position: Optional[float] = None
    mean_read_length: Optional[float] = None


def cluster_5prime_ends(
    bam_path: str,
    chrom: Optional[str] = None,
    strand: Optional[str] = None,
    window: int = 10,
    min_reads: int = 1,
) -> Iterator[Tuple[Tuple[str, str, int], List[pysam.AlignedSegment]]]:
    """
    Cluster reads by their 5' end positions.

    Uses a sliding window approach to group reads with similar 5' ends.

    Args:
        bam_path: Path to sorted BAM file
        chrom: Optional chromosome filter
        strand: Optional strand filter ('+' or '-')
        window: Maximum distance between 5' ends in same cluster (bp)
        min_reads: Minimum reads per cluster

    Yields:
        Tuples of ((chrom, strand, position), [reads])
        where position is the modal 5' end position in the cluster
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')

    # Group reads by (chrom, strand, 5' position)
    reads_by_end: Dict[Tuple[str, str, int], List] = defaultdict(list)

    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        read_strand = '-' if read.is_reverse else '+'
        if strand and read_strand != strand:
            continue

        read_chrom = read.reference_name
        five_prime = get_read_5prime_position(read, read_strand)

        key = (read_chrom, read_strand, five_prime)
        reads_by_end[key].append(read)

    bam.close()

    # Group by (chrom, strand)
    chrom_strand_keys: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for (c, s, pos) in reads_by_end.keys():
        chrom_strand_keys[(c, s)].append(pos)

    # Process each (chrom, strand) group
    for (c, s), positions in chrom_strand_keys.items():
        positions.sort()

        if not positions:
            continue

        # Cluster adjacent positions
        cluster_start = positions[0]
        cluster_reads = []
        cluster_positions = []

        for pos in positions:
            if pos - cluster_start <= window or not cluster_positions:
                cluster_positions.append(pos)
                cluster_reads.extend(reads_by_end[(c, s, pos)])
            else:
                if len(cluster_reads) >= min_reads:
                    pos_counts = defaultdict(int)
                    for p in cluster_positions:
                        pos_counts[p] += len(reads_by_end[(c, s, p)])
                    modal_pos = max(pos_counts.keys(), key=lambda p: pos_counts[p])

                    yield ((c, s, modal_pos), cluster_reads)

                cluster_start = pos
                cluster_positions = [pos]
                cluster_reads = list(reads_by_end[(c, s, pos)])

        # Output final cluster
        if len(cluster_reads) >= min_reads:
            pos_counts = defaultdict(int)
            for p in cluster_positions:
                pos_counts[p] += len(reads_by_end[(c, s, p)])
            modal_pos = max(pos_counts.keys(), key=lambda p: pos_counts[p])

            yield ((c, s, modal_pos), cluster_reads)


def aggregate_attributions_for_5prime_end(
    attributions: List[Tuple[str, ...]],
    min_fraction: float = 0.01,
) -> str:
    """
    Aggregate gene attributions from multiple reads sharing a 5' end.

    This is similar to 3' end attribution but from the opposite perspective:
    for each 5' end cluster, we track what genes the reads' 3' ends fall into.

    Args:
        attributions: List of tuples, each containing gene IDs for one read.
        min_fraction: Minimum fraction to report (default 1%)

    Returns:
        Attribution string in format: "geneA,50%;geneA,geneB,25%;none,25%"
    """
    if not attributions:
        return "none,100%"

    counts = defaultdict(int)
    total = len(attributions)

    for genes in attributions:
        if not genes:
            key = ('none',)
        else:
            key = tuple(sorted(genes))
        counts[key] += 1

    parts = []
    for genes, count in sorted(counts.items(), key=lambda x: -x[1]):
        fraction = count / total
        if fraction >= min_fraction:
            pct = round(fraction * 100)
            gene_str = ','.join(genes)
            parts.append(f"{gene_str},{pct}%")

    if not parts:
        parts = ["none,100%"]

    return ';'.join(parts)


def aggregate_5prime_clusters(
    bam_path: str,
    annotation_df: pd.DataFrame,
    chrom: Optional[str] = None,
    strand: Optional[str] = None,
    window: int = 10,
    min_reads: int = 1,
    include_read_ids: bool = False,
) -> pd.DataFrame:
    """
    Aggregate reads by 5' end position with gene attribution from 3' ends.

    For each cluster of 5' ends, this computes:
    - Cluster position and boundaries
    - Gene attribution based on where 3' ends fall
    - 3' end distribution metrics

    Args:
        bam_path: Path to sorted BAM file
        annotation_df: Gene annotation DataFrame for attribution
        chrom: Optional chromosome filter
        strand: Optional strand filter
        window: Clustering window (bp)
        min_reads: Minimum reads per cluster
        include_read_ids: Include list of read IDs (increases memory)

    Returns:
        DataFrame with columns:
            - chrom, strand, position: Cluster identity
            - cluster_start, cluster_end: Cluster boundaries
            - n_reads: Number of reads
            - attribution: Gene attribution string (from 3' ends)
            - primary_gene: Most common gene (at 3' ends)
            - primary_fraction: Fraction with primary gene
            - mean_3prime, median_3prime: 3' end statistics
            - mean_length: Mean read length
    """
    # Build gene interval trees for 3' end attribution
    interval_trees = build_cds_interval_tree(annotation_df)

    clusters = []

    for (c, s, pos), reads in cluster_5prime_ends(
        bam_path, chrom, strand, window, min_reads
    ):
        three_prime_positions = []
        read_lengths = []
        attributions = []
        read_ids_list = []

        for read in reads:
            three_prime = get_read_3prime_position(read, s)
            three_prime_positions.append(three_prime)

            start, end = get_read_body_interval(read)
            read_lengths.append(end - start)

            # Gene attribution from read body (same as 3' end module)
            # This captures what gene the read belongs to
            genes = compute_read_gene_attribution(read, interval_trees, c)
            attributions.append(tuple(genes))

            if include_read_ids:
                read_ids_list.append(read.query_name)

        # Compute cluster boundaries
        five_prime_positions = [get_read_5prime_position(r, s) for r in reads]
        cluster_start = min(five_prime_positions)
        cluster_end = max(five_prime_positions)

        # Aggregate attribution
        attribution_string = aggregate_attributions_for_5prime_end(attributions)

        # Find primary gene
        gene_counts = defaultdict(int)
        n_intergenic = 0
        for genes in attributions:
            if genes:
                for g in genes:
                    gene_counts[g] += 1
            else:
                n_intergenic += 1

        if gene_counts:
            primary_gene = max(gene_counts.keys(), key=lambda g: gene_counts[g])
            primary_fraction = gene_counts[primary_gene] / len(attributions)
        else:
            primary_gene = 'none'
            primary_fraction = 1.0

        cluster_data = {
            'chrom': c,
            'strand': s,
            'position': pos,
            'cluster_start': cluster_start,
            'cluster_end': cluster_end,
            'n_reads': len(reads),
            'attribution': attribution_string,
            'primary_gene': primary_gene,
            'primary_fraction': round(primary_fraction, 3),
            'n_intergenic': n_intergenic,
            'mean_3prime': round(np.mean(three_prime_positions), 1),
            'median_3prime': int(np.median(three_prime_positions)),
            'mean_length': round(np.mean(read_lengths), 1),
        }

        if include_read_ids:
            cluster_data['read_ids'] = ','.join(read_ids_list)

        clusters.append(cluster_data)

    df = pd.DataFrame(clusters)

    if len(df) > 0:
        df = df.sort_values(['chrom', 'strand', 'position'])
        df = df.reset_index(drop=True)

    logger.info(f"Aggregated {len(df)} 5' end clusters from {bam_path}")

    return df


def export_5prime_clusters(
    clusters_df: pd.DataFrame,
    output_path: str,
    format: str = 'tsv',
) -> str:
    """
    Export 5' end clusters to file.

    Args:
        clusters_df: DataFrame from aggregate_5prime_clusters()
        output_path: Output file path
        format: 'tsv', 'csv', or 'parquet'

    Returns:
        Path to output file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format == 'tsv':
        clusters_df.to_csv(output_path, sep='\t', index=False)
    elif format == 'csv':
        clusters_df.to_csv(output_path, index=False)
    elif format == 'parquet':
        clusters_df.to_parquet(output_path, index=False)
    else:
        raise ValueError(f"Unknown format: {format}")

    logger.info(f"Exported {len(clusters_df)} clusters to {output_path}")
    return str(output_path)
