"""
Junction Aggregation Module for RECTIFY.

This module aggregates splice junction data from:
1. Fully spliced reads (CIGAR N operations)
2. Partial junction evidence (soft-clips at splice sites)

The output merges with nanocompass multi-aligner consensus and provides:
- Junction-level read counts
- Canonical vs non-canonical classification
- Annotation status (exact, near, novel)
- Confidence scores from multi-aligner agreement

Author: Kevin R. Roy
"""

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from ..terminal_exon_refiner import (
    PartialJunctionEvidence,
    detect_partial_junction_crossings,
    compute_junction_counts_with_ambiguity,
)

logger = logging.getLogger(__name__)


# Canonical splice site dinucleotides
CANONICAL_5SS = {'GT', 'GC'}  # 5' splice site (donor)
CANONICAL_3SS = {'AG'}  # 3' splice site (acceptor)


@dataclass
class JunctionStats:
    """Statistics for a single splice junction."""
    chrom: str
    strand: str
    intron_start: int  # 0-based, first base of intron
    intron_end: int  # 0-based exclusive, one past last base of intron
    intron_id: str = ""  # From annotation or auto-generated

    # Read counts
    full_junction_reads: int = 0  # Reads with CIGAR N operation
    partial_junction_reads: int = 0  # Reads with soft-clip evidence
    partial_confident: int = 0  # High-confidence partial evidence
    partial_ambiguous: int = 0  # Ambiguous partial evidence
    total_reads: float = 0.0  # After ambiguity resolution

    # Splice site motifs
    five_ss_dinuc: str = ""  # e.g., "GT"
    three_ss_dinuc: str = ""  # e.g., "AG"
    is_canonical: bool = False

    # Annotation status
    annotation_status: str = ""  # 'exact', 'near', 'novel'
    distance_from_annotated: int = 0  # bp from nearest annotated junction

    # Multi-aligner consensus (if available)
    n_aligners: int = 0
    confidence: str = ""  # 'high', 'medium', 'low'


def extract_junctions_from_cigar(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """
    Extract splice junctions from CIGAR N operations.

    Args:
        read: pysam AlignedSegment

    Returns:
        List of (intron_start, intron_end) tuples in 0-based coords.
        intron_end is exclusive (one past the last intron base).
    """
    if not read.cigartuples:
        return []

    junctions = []
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op == 3:  # N = skipped region (intron)
            intron_start = ref_pos
            intron_end = ref_pos + length
            junctions.append((intron_start, intron_end))
            ref_pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X consume reference
            ref_pos += length
        # I, S, H don't consume reference

    return junctions


def aggregate_junctions(
    bam_path: str,
    genome: Optional[Dict[str, str]] = None,
    annotation_df: Optional[pd.DataFrame] = None,
    min_reads: int = 1,
) -> pd.DataFrame:
    """
    Aggregate splice junctions from CIGAR N operations.

    Args:
        bam_path: Path to sorted BAM file
        genome: Optional genome dict for motif extraction
        annotation_df: Optional annotation for status classification
        min_reads: Minimum reads per junction

    Returns:
        DataFrame with junction statistics
    """
    # Count junctions
    junction_counts: Dict[Tuple[str, int, int, str], int] = defaultdict(int)
    junction_reads: Dict[Tuple[str, int, int, str], List[str]] = defaultdict(list)

    bam = pysam.AlignmentFile(bam_path, 'rb')

    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        strand = '-' if read.is_reverse else '+'
        chrom = read.reference_name

        junctions = extract_junctions_from_cigar(read)

        for intron_start, intron_end in junctions:
            key = (chrom, intron_start, intron_end, strand)
            junction_counts[key] += 1
            junction_reads[key].append(read.query_name)

    bam.close()

    # Build result DataFrame
    results = []

    for (chrom, intron_start, intron_end, strand), count in junction_counts.items():
        if count < min_reads:
            continue

        intron_length = intron_end - intron_start

        # Extract splice site motifs if genome provided
        five_ss_dinuc = ""
        three_ss_dinuc = ""
        is_canonical = False

        if genome and chrom in genome:
            genome_seq = genome[chrom]
            # 5'SS: first 2 bases of intron
            five_ss_dinuc = genome_seq[intron_start:intron_start + 2].upper()
            # 3'SS: last 2 bases of intron
            three_ss_dinuc = genome_seq[intron_end - 2:intron_end].upper()

            is_canonical = (five_ss_dinuc in CANONICAL_5SS and
                           three_ss_dinuc in CANONICAL_3SS)

        results.append({
            'chrom': chrom,
            'strand': strand,
            'intron_start': intron_start,
            'intron_end': intron_end,
            'intron_length': intron_length,
            'full_junction_reads': count,
            'five_ss_dinuc': five_ss_dinuc,
            'three_ss_dinuc': three_ss_dinuc,
            'is_canonical': is_canonical,
        })

    df = pd.DataFrame(results)

    if len(df) > 0:
        df = df.sort_values(['chrom', 'strand', 'intron_start'])
        df = df.reset_index(drop=True)

    logger.info(f"Aggregated {len(df)} junctions from {bam_path}")

    return df


def merge_with_partial_evidence(
    junction_df: pd.DataFrame,
    partial_results: Dict,
    ambiguous_mode: str = 'proportional',
) -> pd.DataFrame:
    """
    Merge junction counts with partial junction evidence.

    This adds counts from reads that are truncated at splice sites
    but have soft-clip evidence of splicing.

    Args:
        junction_df: DataFrame from aggregate_junctions()
        partial_results: Dict from detect_partial_junction_crossings()
        ambiguous_mode: How to handle ambiguous reads

    Returns:
        Enhanced DataFrame with partial evidence columns
    """
    if len(junction_df) == 0:
        return junction_df

    # Compute counts from partial evidence
    partial_counts = compute_junction_counts_with_ambiguity(
        partial_results, mode=ambiguous_mode
    )

    # Create lookup by junction coordinates
    partial_by_coords: Dict[Tuple[str, str, int, int], Dict] = {}

    for evidence in partial_results['evidence']:
        intron_start, intron_end = evidence.inferred_junction
        key = (evidence.chrom, evidence.strand, intron_start, intron_end)

        if key not in partial_by_coords:
            partial_by_coords[key] = {
                'partial_total': 0,
                'partial_confident': 0,
                'partial_ambiguous': 0,
            }

        partial_by_coords[key]['partial_total'] += 1
        if evidence.confidence == 'high':
            partial_by_coords[key]['partial_confident'] += 1
        if evidence.ambiguous:
            partial_by_coords[key]['partial_ambiguous'] += 1

    # Merge with junction DataFrame
    junction_df = junction_df.copy()

    partial_totals = []
    partial_confidents = []
    partial_ambiguous_list = []
    total_reads = []

    for _, row in junction_df.iterrows():
        key = (row['chrom'], row['strand'], row['intron_start'], row['intron_end'])
        partial_data = partial_by_coords.get(key, {})

        pt = partial_data.get('partial_total', 0)
        pc = partial_data.get('partial_confident', 0)
        pa = partial_data.get('partial_ambiguous', 0)

        partial_totals.append(pt)
        partial_confidents.append(pc)
        partial_ambiguous_list.append(pa)

        # Total = full junction + partial (after ambiguity handling)
        if ambiguous_mode == 'holdout':
            total = row['full_junction_reads'] + (pt - pa)
        else:
            total = row['full_junction_reads'] + pt

        total_reads.append(total)

    junction_df['partial_junction_reads'] = partial_totals
    junction_df['partial_confident'] = partial_confidents
    junction_df['partial_ambiguous'] = partial_ambiguous_list
    junction_df['total_reads'] = total_reads

    logger.info(f"Merged partial evidence for {len(junction_df)} junctions")

    return junction_df


def classify_annotation_status(
    junction_df: pd.DataFrame,
    annotation_introns: pd.DataFrame,
    tolerance: int = 5,
) -> pd.DataFrame:
    """
    Classify junctions as exact, near, or novel relative to annotation.

    Args:
        junction_df: DataFrame with junction coordinates
        annotation_introns: DataFrame with annotated intron coordinates
            Required columns: chrom, strand, start, end, intron_id
        tolerance: Maximum bp difference for 'near' classification

    Returns:
        DataFrame with annotation_status column added
    """
    if len(junction_df) == 0:
        return junction_df

    junction_df = junction_df.copy()

    # Build lookup for annotated introns
    annotated: Dict[Tuple[str, str, int, int], str] = {}
    annotated_by_chrom: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)

    for _, row in annotation_introns.iterrows():
        chrom = row['chrom']
        strand = row.get('strand', '+')
        start = int(row['start'])
        end = int(row['end'])
        intron_id = row.get('intron_id', f"{chrom}:{start}-{end}")

        annotated[(chrom, strand, start, end)] = intron_id
        annotated_by_chrom[(chrom, strand)].append((start, end, intron_id))

    statuses = []
    distances = []
    intron_ids = []

    for _, row in junction_df.iterrows():
        chrom = row['chrom']
        strand = row['strand']
        start = row['intron_start']
        end = row['intron_end']

        key = (chrom, strand, start, end)

        if key in annotated:
            statuses.append('exact')
            distances.append(0)
            intron_ids.append(annotated[key])
        else:
            # Check for near matches
            best_distance = float('inf')
            best_intron_id = ""

            for ann_start, ann_end, ann_id in annotated_by_chrom.get((chrom, strand), []):
                dist = abs(start - ann_start) + abs(end - ann_end)
                if dist < best_distance:
                    best_distance = dist
                    best_intron_id = ann_id

            if best_distance <= tolerance * 2:  # Both ends within tolerance
                statuses.append('near')
                distances.append(best_distance)
                intron_ids.append(best_intron_id)
            else:
                statuses.append('novel')
                distances.append(best_distance if best_distance < float('inf') else -1)
                intron_ids.append("")

    junction_df['annotation_status'] = statuses
    junction_df['distance_from_annotated'] = distances
    junction_df['intron_id'] = intron_ids

    # Log summary
    status_counts = junction_df['annotation_status'].value_counts()
    logger.info(f"Junction annotation status: {dict(status_counts)}")

    return junction_df


def export_junctions(
    junction_df: pd.DataFrame,
    output_path: str,
    format: str = 'tsv',
) -> str:
    """
    Export junction data to file.

    Args:
        junction_df: DataFrame with junction statistics
        output_path: Output file path
        format: 'tsv', 'csv', or 'parquet'

    Returns:
        Path to output file
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if format == 'tsv':
        junction_df.to_csv(output_path, sep='\t', index=False)
    elif format == 'csv':
        junction_df.to_csv(output_path, index=False)
    elif format == 'parquet':
        junction_df.to_parquet(output_path, index=False)
    else:
        raise ValueError(f"Unknown format: {format}")

    logger.info(f"Exported {len(junction_df)} junctions to {output_path}")
    return str(output_path)
