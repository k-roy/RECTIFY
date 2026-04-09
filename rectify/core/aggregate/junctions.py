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


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


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
            if length > 0:  # Skip malformed zero-length N-ops (invalid CIGAR)
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

        if genome and chrom in genome and intron_length >= 4:
            genome_seq = genome[chrom]

            if strand == '-':
                # On the minus strand the canonical donor is at the 3' end of
                # the intron (high coordinate) and the acceptor is at the 5'
                # end (low coordinate); both must be reverse-complemented.
                five_ss_dinuc = _reverse_complement(
                    genome_seq[intron_end - 2:intron_end].upper()
                )
                three_ss_dinuc = _reverse_complement(
                    genome_seq[intron_start:intron_start + 2].upper()
                )
            else:
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


def filter_junctions(
    junction_df: pd.DataFrame,
    min_reads: int = 5,
    min_aligners: int = 1,
    require_canonical: bool = False,
    require_annotated_or_near: bool = False,
) -> pd.DataFrame:
    """Apply multi-tier quality filters to a junction DataFrame.

    Filtering tiers:
      1. Minimum read support (min_reads)
      2. Multi-aligner consensus (min_aligners; based on n_aligners column if present)
      3. Canonical splice sites only (require_canonical; if is_canonical column present)
      4. Annotation status filter (require_annotated_or_near)

    Args:
        junction_df: DataFrame as produced by aggregate_junctions / classify_annotation_status.
        min_reads: Minimum full_junction_reads (or total_reads if former absent).
        min_aligners: Minimum n_aligners value. Ignored if column absent.
        require_canonical: If True, keep only canonical (GT-AG / GC-AG / AT-AC) junctions.
        require_annotated_or_near: If True, keep only 'exact' or 'near' annotation_status.

    Returns:
        Filtered copy of the DataFrame.
    """
    if junction_df.empty:
        return junction_df

    mask = pd.Series(True, index=junction_df.index)

    # Tier 1: read support
    read_col = 'full_junction_reads' if 'full_junction_reads' in junction_df.columns else 'total_reads'
    if read_col in junction_df.columns:
        mask &= junction_df[read_col] >= min_reads

    # Tier 2: multi-aligner consensus
    if min_aligners > 1 and 'n_aligners' in junction_df.columns:
        mask &= junction_df['n_aligners'] >= min_aligners

    # Tier 3: canonical splice sites
    if require_canonical and 'is_canonical' in junction_df.columns:
        mask &= junction_df['is_canonical']

    # Tier 4: annotation status
    if require_annotated_or_near and 'annotation_status' in junction_df.columns:
        mask &= junction_df['annotation_status'].isin(('exact', 'near'))

    filtered = junction_df[mask].copy()
    logger.info(
        f"filter_junctions: {len(junction_df)} → {len(filtered)} junctions "
        f"(min_reads={min_reads}, min_aligners={min_aligners}, "
        f"canonical={require_canonical}, annot={require_annotated_or_near})"
    )
    return filtered


def resolve_homopolymer_ambiguity(
    junction_df: pd.DataFrame,
    genome: Dict[str, str],
    motif_scorer=None,
    max_shift: int = 10,
) -> pd.DataFrame:
    """Collapse homopolymer-ambiguous junctions to their best representative.

    When identical bases flank an intron boundary (e.g., ...AAAA|AAAA...), the
    aligner may place the junction at several equivalent positions. This function
    clusters such junctions and keeps only the one with the best splice motif score
    (or, if no scorer is provided, the one with the most reads).

    Requires columns: chrom, intron_start, intron_end, strand.
    The best representative inherits the sum of all reads in the cluster.

    Args:
        junction_df: Junction DataFrame.
        genome: Dict[chrom, sequence] for sequence access.
        motif_scorer: SpliceMotifScorer instance (from rectify.utils.splice_motif).
            If None, selects the position with the most reads.
        max_shift: Maximum shift window for homopolymer detection.

    Returns:
        Deduplicated DataFrame with ambiguous clusters merged.
    """
    from rectify.utils.splice_motif import get_splice_site_sequences

    if junction_df.empty:
        return junction_df

    required = {'chrom', 'intron_start', 'intron_end', 'strand'}
    if not required.issubset(junction_df.columns):
        logger.warning("resolve_homopolymer_ambiguity: missing required columns, skipping")
        return junction_df

    read_col = 'full_junction_reads' if 'full_junction_reads' in junction_df.columns else 'total_reads'

    processed: set = set()
    keep_rows: List[int] = []
    drop_rows: set = set()
    merged_counts: Dict[int, float] = {}

    for idx, row in junction_df.iterrows():
        if idx in processed:
            continue

        chrom = row['chrom']
        start = int(row['intron_start'])
        end = int(row['intron_end'])
        strand = row['strand']
        seq = genome.get(chrom, '')

        if not seq:
            keep_rows.append(idx)
            processed.add(idx)
            continue

        # Detect equivalent positions by sliding the boundary through identical bases
        cluster_coords: List[Tuple[int, int]] = [(start, end)]

        # Slide left
        for shift in range(1, max_shift + 1):
            if start - shift < 0:
                break
            span = seq[start - shift:end - shift]
            if span and len(set(span)) == 1:
                cluster_coords.append((start - shift, end - shift))
            else:
                break

        # Slide right
        for shift in range(1, max_shift + 1):
            if end + shift > len(seq):
                break
            span = seq[start + shift:end + shift]
            if span and len(set(span)) == 1:
                cluster_coords.append((start + shift, end + shift))
            else:
                break

        if len(cluster_coords) == 1:
            # Not ambiguous
            keep_rows.append(idx)
            processed.add(idx)
            continue

        # Find all rows in this cluster
        cluster_mask = (
            (junction_df['chrom'] == chrom)
            & (junction_df['strand'] == strand)
            & (junction_df['intron_start'].isin([c[0] for c in cluster_coords]))
            & (junction_df['intron_end'].isin([c[1] for c in cluster_coords]))
        )
        cluster_idxs = junction_df.index[cluster_mask].tolist()
        cluster_rows = junction_df.loc[cluster_idxs]
        total_reads = cluster_rows[read_col].sum() if read_col in cluster_rows.columns else len(cluster_rows)

        # Choose best representative
        if motif_scorer is not None:
            best_score = float('inf')
            best_idx = cluster_idxs[0]
            for cidx in cluster_idxs:
                crow = junction_df.loc[cidx]
                five_seq, three_seq = get_splice_site_sequences(
                    genome, chrom, int(crow['intron_start']), int(crow['intron_end']), strand
                )
                score = motif_scorer.score_five_ss(five_seq) + motif_scorer.score_three_ss(three_seq)
                if score < best_score:
                    best_score = score
                    best_idx = cidx
        else:
            # Most reads wins
            best_idx = cluster_rows[read_col].idxmax() if read_col in cluster_rows.columns else cluster_idxs[0]

        keep_rows.append(best_idx)
        merged_counts[best_idx] = total_reads
        processed.update(cluster_idxs)

    result = junction_df.loc[keep_rows].copy()
    for idx, count in merged_counts.items():
        if idx in result.index and read_col in result.columns:
            result.loc[idx, read_col] = count

    n_merged = len(junction_df) - len(result)
    logger.info(
        f"resolve_homopolymer_ambiguity: {len(junction_df)} → {len(result)} junctions "
        f"({n_merged} homopolymer-ambiguous positions merged)"
    )
    return result.reset_index(drop=True)


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
