#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Splice-Aware 5' End Correction for RECTIFY

This module corrects the 5' end position of reads based on splice junction information.
Many introns are close to the TSS, and aligners may place the 5' end within an intron
rather than at the true transcription start site.

Problem:
    When the first exon is very short, the read's leftmost aligned position may be:
    - Within the first exon (correct 5' end)
    - Within an intron (needs correction to exon boundary)

Solution:
    1. Check if the read's 5' end falls within an annotated intron
    2. If so, shift the 5' end to the nearest upstream exon boundary
    3. Use consensus junctions to validate/refine the correction

Coordinate System:
    - 0-based, half-open coordinates (consistent with pysam/BED)
    - Plus strand: 5' end = leftmost position
    - Minus strand: 5' end = rightmost position

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-24
"""

from typing import List, Tuple, Dict, Optional, Set
from dataclasses import dataclass
import pysam

from ..utils.genome import standardize_chrom_name
from ..config import CHROM_TO_GENOME

try:
    from intervaltree import IntervalTree
    HAS_INTERVALTREE = True
except ImportError:
    HAS_INTERVALTREE = False


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class FivePrimeCorrection:
    """Result of 5' end correction."""
    five_prime_raw: int             # Original alignment position
    five_prime_corrected: int       # After splice-aware correction
    first_exon_start: Optional[int] # Start of first exon (may be None)
    starts_in_intron: bool          # True if correction was applied
    correction_bp: int              # Bases shifted (0 if no correction)
    correction_reason: str          # Explanation of correction


# =============================================================================
# Helper Functions
# =============================================================================

def get_read_5prime_position(read: pysam.AlignedSegment, strand: Optional[str] = None) -> int:
    """
    Get the 5' end genomic position from a read.

    The 5' end is the transcription start site (TSS) end of the RNA molecule.

    Args:
        read: pysam AlignedSegment
        strand: Optional strand override ('+' or '-')

    Returns:
        5' end position (0-based)
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    if strand == '+':
        return read.reference_start
    else:
        ref_end = read.reference_end
        if ref_end is None:
            return None  # Unmapped read — no valid position
        return ref_end - 1


def build_intron_interval_tree(
    annotation_df: 'pd.DataFrame',
    feature_types: Optional[List[str]] = None,
) -> Dict[Tuple[str, str], 'IntervalTree']:
    """
    Build IntervalTree index for introns from annotation.

    This identifies intron regions (gaps between CDS/exon features).

    Args:
        annotation_df: DataFrame with gene annotations
            Required columns: chrom, start, end, strand
            Optional columns: gene_id, feature_type
        feature_types: Feature types to use for exon boundaries (default: ['CDS', 'exon'])

    Returns:
        Dict mapping (chrom, strand) to IntervalTree of introns
    """
    if not HAS_INTERVALTREE:
        raise ImportError("intervaltree is required. Install with: pip install intervaltree")

    import pandas as pd
    from collections import defaultdict

    if feature_types is None:
        feature_types = ['CDS', 'exon', 'mRNA']

    # Filter to relevant feature types
    if 'feature_type' in annotation_df.columns:
        df = annotation_df[annotation_df['feature_type'].isin(feature_types)].copy()
    else:
        df = annotation_df.copy()

    # Build exon intervals per gene
    intron_trees = defaultdict(IntervalTree)

    # Group by gene and strand
    if 'gene_id' in df.columns:
        group_cols = ['chrom', 'strand', 'gene_id']
    else:
        group_cols = ['chrom', 'strand']

    for group_key, group_df in df.groupby(group_cols):
        if len(group_cols) == 3:
            chrom, strand, gene_id = group_key
        else:
            chrom, strand = group_key
            gene_id = "unknown"

        # Sort exons by position
        exons = group_df.sort_values('start')

        # Find introns (gaps between consecutive exons)
        exon_list = list(zip(exons['start'], exons['end']))

        for i in range(len(exon_list) - 1):
            _, exon1_end = exon_list[i]
            exon2_start, _ = exon_list[i + 1]

            # Intron is the gap between exons
            if exon2_start > exon1_end:
                intron_start = exon1_end
                intron_end = exon2_start

                key = (chrom, strand)
                intron_trees[key][intron_start:intron_end] = {
                    'gene_id': gene_id,
                    'upstream_exon_end': exon1_end,
                    'downstream_exon_start': exon2_start,
                    'intron_index': i,
                }

    return dict(intron_trees)


def find_overlapping_introns(
    chrom: str,
    position: int,
    strand: str,
    intron_trees: Dict[Tuple[str, str], 'IntervalTree'],
) -> List[Dict]:
    """
    Find introns that overlap a position.

    Args:
        chrom: Chromosome name
        position: Genomic position (0-based)
        strand: Strand ('+' or '-')
        intron_trees: Dict from build_intron_interval_tree()

    Returns:
        List of intron info dicts
    """
    key = (chrom, strand)

    if key not in intron_trees:
        return []

    tree = intron_trees[key]
    overlaps = tree[position]

    results = []
    for interval in overlaps:
        intron_info = interval.data.copy()
        intron_info['intron_start'] = interval.begin
        intron_info['intron_end'] = interval.end
        results.append(intron_info)

    return results


# =============================================================================
# Core Correction Functions
# =============================================================================

def correct_5prime_for_splicing(
    read: pysam.AlignedSegment,
    intron_trees: Optional[Dict[Tuple[str, str], 'IntervalTree']] = None,
    read_junctions: Optional[List[Tuple[int, int]]] = None,
    strand: Optional[str] = None,
) -> FivePrimeCorrection:
    """
    Correct 5' end position based on splice junction information.

    If the read starts within an intron (near TSS), the true 5' end
    should be at the upstream exon boundary, not the raw alignment start.

    Args:
        read: pysam AlignedSegment
        intron_trees: Dict from build_intron_interval_tree() (optional)
        read_junctions: List of (start, end) tuples for junctions in the read
        strand: Strand override, or None to infer from read

    Returns:
        FivePrimeCorrection with correction details
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    chrom = read.reference_name
    five_prime_raw = get_read_5prime_position(read, strand)

    # Default: no correction
    result = FivePrimeCorrection(
        five_prime_raw=five_prime_raw,
        five_prime_corrected=five_prime_raw,
        first_exon_start=None,
        starts_in_intron=False,
        correction_bp=0,
        correction_reason="",
    )

    # Strategy 1: Use read junctions to infer 5' correction
    # If the first junction in the read is very close to the 5' end,
    # the 5' end might be in an intron
    if read_junctions:
        if strand == '+':
            # For plus strand, 5' is leftmost
            # Check if first junction starts very close to 5' end
            first_junction_start = min(j[0] for j in read_junctions)
            first_junction_end = min(j[1] for j in read_junctions if j[0] == first_junction_start)

            distance_to_first_junction = first_junction_start - five_prime_raw

            # If 5' end is within the first junction, correct to upstream exon
            if five_prime_raw >= first_junction_start and five_prime_raw < first_junction_end:
                result.starts_in_intron = True
                result.five_prime_corrected = first_junction_start - 1
                result.first_exon_start = first_junction_start - 1
                result.correction_bp = five_prime_raw - (first_junction_start - 1)
                result.correction_reason = "5' end within read junction - shifted to exon boundary"
                return result

        else:
            # For minus strand, 5' is rightmost
            # Check if last junction ends very close to 5' end
            last_junction_end = max(j[1] for j in read_junctions)
            last_junction_start = max(j[0] for j in read_junctions if j[1] == last_junction_end)

            # If 5' end is within the last junction, correct to downstream exon
            if five_prime_raw >= last_junction_start and five_prime_raw < last_junction_end:
                result.starts_in_intron = True
                result.five_prime_corrected = last_junction_end
                result.first_exon_start = last_junction_end
                result.correction_bp = last_junction_end - five_prime_raw
                result.correction_reason = "5' end within read junction - shifted to exon boundary"
                return result

    # Strategy 2: Use annotation introns
    if intron_trees:
        overlapping_introns = find_overlapping_introns(chrom, five_prime_raw, strand, intron_trees)

        if overlapping_introns:
            # Find the most relevant intron (closest to 5' end)
            if strand == '+':
                # For plus strand, the 5' end (TSS) is at low coordinates.
                # If it lands inside an intron, correct to the last base of the
                # upstream exon: intron_start - 1.
                best_intron = min(overlapping_introns, key=lambda x: x['intron_start'])
                corrected_pos = best_intron['intron_start'] - 1
            else:
                # For minus strand, the 5' end (TSS) is at high coordinates.
                # If it lands inside an intron, correct to the first base of the
                # downstream exon: intron_end.
                best_intron = max(overlapping_introns, key=lambda x: x['intron_end'])
                corrected_pos = best_intron['intron_end']

            result.starts_in_intron = True
            result.five_prime_corrected = corrected_pos
            result.first_exon_start = corrected_pos if strand == '+' else corrected_pos + 1
            result.correction_bp = abs(corrected_pos - five_prime_raw)
            result.correction_reason = f"5' end in annotated intron ({best_intron.get('gene_id', 'unknown')}) - shifted to exon boundary"

    return result


def correct_5prime_batch(
    reads_with_junctions: List[Tuple[pysam.AlignedSegment, List[Tuple[int, int]]]],
    intron_trees: Optional[Dict[Tuple[str, str], 'IntervalTree']] = None,
) -> List[FivePrimeCorrection]:
    """
    Correct 5' end positions for a batch of reads.

    Args:
        reads_with_junctions: List of (read, junctions) tuples
        intron_trees: Dict from build_intron_interval_tree()

    Returns:
        List of FivePrimeCorrection objects
    """
    corrections = []

    for read, junctions in reads_with_junctions:
        correction = correct_5prime_for_splicing(
            read,
            intron_trees=intron_trees,
            read_junctions=junctions,
        )
        corrections.append(correction)

    return corrections


# =============================================================================
# Post-Consensus 3'SS Truncation Rescue
# =============================================================================

def _get_5prime_softclip_len(read: pysam.AlignedSegment) -> int:
    """Return the explicit 5' soft-clip length (S op adjacent to the 5' end)."""
    if not read.cigartuples:
        return 0
    if read.is_reverse:
        last_op, last_len = read.cigartuples[-1]
        return last_len if last_op == 4 else 0
    else:
        first_op, first_len = read.cigartuples[0]
        return first_len if first_op == 4 else 0


def _extract_5prime_rescue_seq(
    read: pysam.AlignedSegment,
    genome_seq: str = '',
    scan_ref_bp: int = 50,
) -> str:
    """Unified 5'-end rescue sequence extractor.

    Scans up to ``scan_ref_bp`` reference bases from the 5' alignment end and
    returns query bases from the 5' end up to and including the **last**
    imperfect alignment position in that window.

    **Trigger**: any S (soft-clip), X (mismatch), I (insertion), or D (deletion
    of any size) op detected in the CIGAR.  For reads that use M ops instead of
    =/X (e.g. mapPacBio), a reference-vs-query comparison is performed when
    ``genome_seq`` is provided; mismatching M positions are treated as imperfect.

    **Clipping to last error**: clean = ops beyond the last imperfect position
    are excluded, preventing correctly-aligned exon-2 tail bases from inflating
    edit distance in the downstream junction-matching step.

    Stops at N ops (existing splice junctions).
    Returns ``''`` if no imperfect op is found within the scan window.

    Replaces three earlier helpers:
        - ``_extract_5prime_terminal_error_seq``  (Case 2 mismatch scan)
        - ``_extract_5prime_deletion_bridged_seq`` (Case 2b large-D scan)
        and the explicit soft-clip slice (Case 1) in ``rescue_3ss_truncation``.
    """
    if not read.cigartuples or not read.query_sequence:
        return ''

    query_seq = read.query_sequence
    n_query = len(query_seq)
    _qc = frozenset([0, 1, 4, 7, 8])              # M, I, S, =, X (consume query)
    _rc = frozenset([0, 2, 7, 8])                  # M, D, =, X (consume reference)
    _explicit_imperfect = frozenset([1, 2, 4, 8])  # I, D, S, X (unambiguously bad)

    # --- Phase 1: CIGAR-based detection (I, D, S, X ops) ---
    # Iterate from the 5' end of the read (reversed CIGAR for minus strand).
    query_collected = 0   # query bases accumulated from 5' end
    ref_scanned = 0       # reference bases consumed (determines scan window)
    last_imp_q = 0        # query bases at the EXCLUSIVE end of last imperfect op
    found_explicit = False
    has_m_ops = False     # any M op seen (may need genome fallback)

    ops = reversed(read.cigartuples) if read.is_reverse else iter(read.cigartuples)
    for op, length in ops:
        if op == 5:   # H: not in query_sequence
            continue
        if op == 3:   # N: existing splice junction — stop
            break
        n_qb = length if op in _qc else 0
        n_rb = length if op in _rc else 0
        if op in _explicit_imperfect:
            found_explicit = True
            last_imp_q = query_collected + n_qb  # exclusive end (includes this op)
        if op == 0:
            has_m_ops = True
        query_collected += n_qb
        ref_scanned += n_rb
        if ref_scanned >= scan_ref_bp:
            break

    if found_explicit and last_imp_q > 0:
        if read.is_reverse:
            return query_seq[n_query - last_imp_q:]
        else:
            return query_seq[:last_imp_q]

    # --- Phase 2: M-op fallback (mapPacBio uses M instead of =/X) ---
    # Only runs when Phase 1 found no explicit imperfect ops but M ops are present.
    if not has_m_ops or not genome_seq:
        return ''

    try:
        pairs = read.get_aligned_pairs()
    except Exception:
        return ''

    gs = len(genome_seq)
    ref_count = 0
    last_imp_qp = -1  # query index of last mismatching position

    if read.is_reverse:
        for qp, rp in reversed(pairs):
            if rp is not None:
                ref_count += 1
            if qp is None:
                continue  # D op: no query base
            if rp is None:
                last_imp_qp = qp  # I/S: imperfect
            elif rp < gs:
                gb = genome_seq[rp].upper()
                rb = query_seq[qp].upper()
                if gb != 'N' and rb != gb:
                    last_imp_qp = qp
            if ref_count >= scan_ref_bp:
                break
    else:
        for qp, rp in pairs:
            if rp is not None:
                ref_count += 1
            if qp is None:
                continue
            if rp is None:
                last_imp_qp = qp
            elif rp < gs:
                gb = genome_seq[rp].upper()
                rb = query_seq[qp].upper()
                if gb != 'N' and rb != gb:
                    last_imp_qp = qp
            if ref_count >= scan_ref_bp:
                break

    if last_imp_qp < 0:
        return ''

    if read.is_reverse:
        # last_imp_qp is the query index; for minus strand the 5' end is at the
        # right of query_seq.  We want all bases from last_imp_qp to the right end.
        return query_seq[last_imp_qp:]
    else:
        return query_seq[:last_imp_qp + 1]


def _get_intronic_query_bases(
    read: pysam.AlignedSegment,
    clip_boundary: int,
    strand: str,
) -> str:
    """Return query bases that map to the intron-side of ``clip_boundary``.

    Iterates the CIGAR from the 5' end (reversed for minus strand) and
    accumulates query bases for every op whose reference span lies entirely
    at or beyond ``clip_boundary`` (i.e., inside the intron or at the exact
    boundary).  For an op that partially crosses the boundary only the
    proportional intronic slice is included.  Stops at N ops.

    For minus strand ``clip_boundary`` is ``intron_start`` (low ref coord);
    for plus strand it is ``intron_end`` (high ref coord).

    These bases are used by :func:`rescue_3ss_truncation` to drive the
    local alignment for the exon-1 CIGAR, ensuring the exon CIGAR has
    exactly as many query-consuming bases as the CIGAR trimming step in
    :func:`~rectify.core.bam_writer.reroute_intronic_tail_5prime_via_junction`
    will remove.
    """
    if not read.cigartuples or not read.query_sequence:
        return ''

    query_seq = read.query_sequence
    n_query = len(query_seq)
    _qc = frozenset([0, 1, 4, 7, 8])  # M, I, S, =, X — consume query
    _rc = frozenset([0, 2, 7, 8])      # M, D, =, X — consume reference

    query_bases = 0

    if strand == '-':
        ref_pos = read.reference_end
        if ref_pos is None:
            return ''
        for op, length in reversed(read.cigartuples):
            if op == 5:   # H: skip
                continue
            if op == 3:   # N: stop
                break
            n_qb = length if op in _qc else 0
            n_rb = length if op in _rc else 0
            if n_rb > 0:
                new_ref = ref_pos - n_rb
                if new_ref >= clip_boundary:
                    # Entirely inside intron
                    query_bases += n_qb
                    ref_pos = new_ref
                elif ref_pos > clip_boundary:
                    # Partially crosses boundary: include only the intronic slice
                    overlap_rb = ref_pos - clip_boundary
                    overlap_qb = round(n_qb * overlap_rb / n_rb) if n_qb else 0
                    query_bases += overlap_qb
                    break
                else:
                    break  # Entirely in exon 2
            else:
                # Query-only op (I/S): attach to the current ref position.
                # Use strict > so that insertions exactly AT clip_boundary are
                # excluded — the reroute trimmer also excludes them (its loop
                # breaks when cur_end <= clip_boundary, before consuming the
                # boundary insertion).
                if ref_pos > clip_boundary:
                    query_bases += n_qb
    else:
        ref_pos = read.reference_start
        for op, length in read.cigartuples:
            if op == 5:
                continue
            if op == 3:
                break
            n_qb = length if op in _qc else 0
            n_rb = length if op in _rc else 0
            if n_rb > 0:
                new_ref = ref_pos + n_rb
                if new_ref <= clip_boundary:
                    query_bases += n_qb
                    ref_pos = new_ref
                elif ref_pos < clip_boundary:
                    overlap_rb = clip_boundary - ref_pos
                    overlap_qb = round(n_qb * overlap_rb / n_rb) if n_qb else 0
                    query_bases += overlap_qb
                    break
                else:
                    break
            else:
                if ref_pos <= clip_boundary:
                    query_bases += n_qb

    if query_bases == 0:
        return ''
    if strand == '-':
        return query_seq[n_query - query_bases:]
    else:
        return query_seq[:query_bases]


def _get_n_op_intervals(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """Return (start, end) genomic intervals for every N-op (intron skip) in the CIGAR."""
    intervals: List[Tuple[int, int]] = []
    if not read.cigartuples or read.reference_start is None:
        return intervals
    pos = read.reference_start
    for op, length in read.cigartuples:
        if op == 3:   # N — intron skip
            intervals.append((pos, pos + length))
            pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X — consume reference
            pos += length
        # I (1), S (4), H (5), P (6) do not consume reference
    return intervals


def _edit_distance(s1: str, s2: str) -> int:
    """Simple edit distance (Levenshtein) for short sequences."""
    n, m = len(s1), len(s2)
    if n == 0:
        return m
    if m == 0:
        return n
    dp = list(range(m + 1))
    for i in range(1, n + 1):
        prev, dp[0] = dp[0], i
        for j in range(1, m + 1):
            temp = dp[j]
            dp[j] = min(dp[j] + 1, dp[j - 1] + 1,
                        prev + (0 if s1[i - 1] == s2[j - 1] else 1))
            prev = temp
    return dp[m]


def _hp_edit_distance(s1: str, s2: str) -> float:
    """Edit distance with 0.5 penalty for indels within homopolymer runs.

    Nanopore sequencers under/over-call homopolymer run lengths.  A deletion
    or insertion of a base that is part of a run of identical bases (i.e. the
    indel base matches its immediate neighbour in the same sequence) is given
    half the normal gap penalty (0.5 instead of 1.0).  Substitutions always
    cost 1.0 regardless of context.

    A base is considered part of a homopolymer run if it equals either the
    preceding or the following character in the *same* string.
    """
    n, m = len(s1), len(s2)
    if n == 0:
        return float(m)
    if m == 0:
        return float(n)

    def _del_cost(i: int) -> float:
        """Cost to delete s1[i-1] (gap in s2)."""
        c = s1[i - 1]
        if (i >= 2 and s1[i - 2] == c) or (i < n and s1[i] == c):
            return 0.5
        return 1.0

    def _ins_cost(j: int) -> float:
        """Cost to insert s2[j-1] (gap in s1)."""
        c = s2[j - 1]
        if (j >= 2 and s2[j - 2] == c) or (j < m and s2[j] == c):
            return 0.5
        return 1.0

    # 2-D DP (sequences are short, so O(n*m) space is fine)
    dp = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + _del_cost(i)
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + _ins_cost(j)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = min(
                    dp[i - 1][j - 1] + 1.0,   # substitution
                    dp[i - 1][j] + _del_cost(i),
                    dp[i][j - 1] + _ins_cost(j),
                )
    return dp[n][m]


# 3'SS acceptor dinucleotide priority (lower = more canonical).
# Plus strand:  last 2 bases of intron before exon 2 (genome[intron_end-2:intron_end])
# Minus strand: first 2 bases of intron (genome[intron_start:intron_start+2]),
#               which is the RC of the RNA-level 3'SS motif.
_ACCEPTOR_PRIORITY_PLUS  = {'AG': 0, 'CG': 1, 'TG': 2, 'AT': 3}
_ACCEPTOR_PRIORITY_MINUS = {'CT': 0, 'CG': 1, 'CA': 2, 'AT': 3}


def rescue_3ss_truncation(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int]],
    strand: Optional[str] = None,
    max_edit_frac: float = 0.2,
    junction_proximity_bp: int = 10,
    scan_bp: int = 50,
) -> Dict:
    """Rescue reads truncated or mis-aligned at the exon 2 / 3' splice site boundary.

    **General approach**: for any read whose 5' alignment end is near (within
    ``junction_proximity_bp`` bp of) a known 3'SS, or whose 5' end falls inside
    an annotated intron, extract the terminal query bases that show any alignment
    imperfection (S, X, I, D of any size) and attempt re-alignment against the
    upstream exon 1 sequence.

    Even a single mismatch or 1 bp intronic overlap is sufficient to trigger
    realignment — the downstream edit-distance check against both exon and intron
    sequence filters false positives.

    Cases handled (in priority order):

    1. **Any terminal imperfection** (S / X / I / D) within ``scan_bp`` ref bases
       of the 5' end: query bases from the 5' end up to and including the last
       imperfect position are extracted and aligned to each candidate junction's
       exon-1 sequence.  Covers soft-clip (Case 1), mapPacBio forced mismatches
       (Case 2), deletion-bridged alignments (Case 2b), single-bp boundary
       mismatches, and small indels.

    2. **Intronic snap** (Case 4): if sequence-based rescue produced no match but
       the 5' alignment end is strictly inside an annotated intron (and no N-op
       already covers it), the corrected position is snapped to the exon-1-side
       boundary.  Fires as a fallback after failed sequence rescue, e.g. when the
       terminal region has only 1 mismatched base that does not align to any
       known exon-1 sequence.

    3. **Proximity-only** (Case 3): alignment ends within ``junction_proximity_bp``
       of a 3'SS but has no imperfect op and is not inside the intron.  Records
       the junction hit without changing the 5' position.

    Args:
        read: pysam AlignedSegment (from rectified BAM)
        genome: chromosome → sequence dict
        candidate_junctions: Set of (chrom, intron_start, intron_end) to test
        strand: Strand override; inferred from read.is_reverse if None
        max_edit_frac: Edit distance / rescue_seq_len threshold for sequence match
        junction_proximity_bp: Max bp between alignment 5' end and intron edge
            to attempt rescue (both sequence-based and proximity-only)
        scan_bp: Reference bases to scan from the 5' end for imperfect-op detection

    Returns:
        Dict with keys:
            'rescued'              bool  — True for sequence-confirmed rescue
            'rescue_type'          str   — 'softclip' | 'mpb_mismatch' | 'intronic_snap'
                                           | 'proximity' | 'none'
            'five_prime_corrected' int   — Updated 5' genomic position
            'rescued_junction'     tuple — (chrom, intron_start, intron_end) or None
            'edit_distance'        float — HP-weighted edit distance; -1 for non-seq rescues
            'query_bp'             int   — Length of rescue sequence used; 0 for proximity
    """
    if strand is None:
        strand = '-' if read.is_reverse else '+'

    chrom = standardize_chrom_name(read.reference_name) if read.reference_name else read.reference_name
    if not chrom:
        return _no_rescue(read, strand)

    # Genome may be keyed by NCBI format (NC_001133.3) or canonical (chrI).
    # Try canonical first, then fall back to NCBI key via CHROM_TO_GENOME.
    genome_seq = genome.get(chrom) or genome.get(CHROM_TO_GENOME.get(chrom, ''))
    if not genome_seq:
        return _no_rescue(read, strand)

    five_clip = _get_5prime_softclip_len(read)

    # --- Determine 5' alignment boundary (first aligned base, after any soft-clip) ---
    # + strand: reference_start is already the first aligned base
    # − strand: reference_end - 1 is the first aligned base in transcript orientation
    if strand == '+':
        align_5prime = read.reference_start
    else:
        if read.reference_end is None:
            return _no_rescue(read, strand)
        align_5prime = read.reference_end - 1

    # --- Collect rescue sequence ---
    # General approach: scan scan_bp reference bases from the 5' end; any imperfect
    # alignment op (S, X, I, D of any size, or M-mismatch if genome is available)
    # triggers extraction of all query bases up to and including the last imperfect
    # position.  Clean exon-2 bases at the tail of the window are excluded to
    # prevent inflating edit distance in the junction-matching step.
    rescue_seq = _extract_5prime_rescue_seq(read, genome_seq=genome_seq, scan_ref_bp=scan_bp)
    rescue_type_candidate = (
        "softclip" if (five_clip > 0 and rescue_seq)
        else ("mpb_mismatch" if rescue_seq else "none")
    )

    # When a soft clip is present, the rescue sequence IS the soft-clipped bases
    # (exon sequence not aligned by the aligner).  The scan_bp extension can pull in
    # aligned bases beyond the soft clip that map INSIDE the intron; including those
    # intron-internal bases contaminates the exon-matching step and causes wrong
    # shifts to score better than the correct exon boundary.  Truncate to the
    # soft-clip length to prevent this.
    if five_clip > 0 and rescue_seq and len(rescue_seq) > five_clip:
        rescue_seq = rescue_seq[:five_clip]

    # --- Try sequence-based rescue against each candidate junction ---
    best_ed: float = -1.0
    best_junction = None
    best_five_prime_corrected = align_5prime
    best_is_canonical = False    # tiebreaker 1: canonical GT/GC donor
    best_in_amb = False          # tiebreaker 2: shift within ambiguity window
    best_shift_abs = 999         # tiebreaker 3: smallest |shift|
    best_acceptor_priority = 4   # tiebreaker 4: 3'SS quality (AG=0..AT=3..other=4)

    # Splice-site boundary ambiguity: when bases flanking a donor/acceptor are
    # repeated (homopolymer runs, tandem dinucleotides, etc.), the local aligner
    # cannot distinguish the correct intron start/end from nearby positions.
    # For each candidate junction we sample a range of shifts derived from the
    # *actual* run-length of matching bases on each side of the boundary:
    #
    #   right-of-boundary run  → how far right the junction can slide without
    #                            changing the sequence context (intron side)
    #   left-of-boundary run   → how far left the junction can slide (exon side)
    #
    # Constraining the search to this data-driven range prevents spurious matches
    # far from the junction (e.g. when exon1-end + intron-start looks identical
    # to intron-end + exon2-start — both windows could give edit_distance=0).
    # A minimum of ±1 is always included to handle coordinate-system off-by-one
    # errors (e.g. 1-based GFF vs 0-based half-open).  Canonical GT-AG is
    # preferred over non-canonical among tied candidates; smaller |shift| is the
    # final tiebreaker to preserve the annotated position when all else is equal.
    _MAX_SS_SHIFT = 15  # hard cap regardless of run length

    if rescue_seq:
        rescue_len = len(rescue_seq)
        _gs = len(genome_seq)
        for j_entry in candidate_junctions:
            j_chrom, intron_start, intron_end = j_entry[0], j_entry[1], j_entry[2]
            if j_chrom != chrom:
                continue

            if strand == '+':
                # Upstream intron: intron_end must be at or just before align_5prime.
                # When mapPacBio extends into the intron, align_5prime < intron_end
                # (dist < 0) — allow rescue if alignment starts inside the intron
                # (intron_start < align_5prime < intron_end), not completely before it.
                dist = align_5prime - intron_end
                if dist < 0:
                    if align_5prime <= intron_start:
                        continue  # alignment is upstream of intron entirely — skip
                    dist = 0  # treat as touching the boundary
                elif dist > junction_proximity_bp:
                    continue
                # Dynamic shift range from run-length at the annotated donor:
                #   r_amb = consecutive intron bases (going right) equal to last exon base
                #   l_amb = consecutive exon bases (going left) equal to first intron base
                # These tell us how many positions the junction can genuinely slide in
                # each direction while keeping the same alignment score.
                if 0 < intron_start < _gs:
                    _leb = genome_seq[intron_start - 1].upper()  # last exon base
                    _fib = genome_seq[intron_start].upper()      # first intron base
                    _r_amb = 0
                    while (_r_amb < _MAX_SS_SHIFT and intron_start + _r_amb < _gs
                           and genome_seq[intron_start + _r_amb].upper() == _leb):
                        _r_amb += 1
                    _l_amb = 0
                    while (_l_amb < _MAX_SS_SHIFT and intron_start - 1 - _l_amb >= 0
                           and genome_seq[intron_start - 1 - _l_amb].upper() == _fib):
                        _l_amb += 1
                else:
                    _r_amb = _l_amb = 0
                # Wide range for discovery (catches imprecise annotations / aligners),
                # at least ±5 bp regardless of local ambiguity.
                _shift_lo = -max(5, _l_amb)
                _shift_hi =  max(5, _r_amb)

                # 3'SS acceptor quality for this junction (fixed; used as outer tiebreaker).
                # Plus strand: last 2 bases of intron = genome[intron_end-2:intron_end]
                _acc_di = genome_seq[intron_end - 2:intron_end].upper() if intron_end >= 2 else 'NN'
                _acceptor_priority = _ACCEPTOR_PRIORITY_PLUS.get(_acc_di, 4)

                _best_local_ed: float = rescue_len + 1
                _best_local_canonical = False
                _best_in_amb = False
                _best_local_shift_abs = max(abs(_shift_lo), _shift_hi) + 1
                exon_seq = ""
                _eff_intron_start = intron_start

                for _shift in range(_shift_lo, _shift_hi + 1):
                    _eff_start = intron_start + _shift
                    if _eff_start <= 0 or _eff_start + 2 > _gs:
                        continue
                    # Canonical 5'SS donors: GT (major spliceosome) or GC (minor)
                    _donor_ok = genome_seq[_eff_start:_eff_start + 2].upper() in ('GT', 'GC')
                    # Whether this shift is within the natural sequence-ambiguity window
                    _in_amb = (-_l_amb <= _shift <= _r_amb)
                    # Cap _off when dist > 0: alignment is already past the intron_end
                    # (in exon-2); sliding the window further left than dist bp would
                    # reach exon-1 sequences that belong to a different, closer junction.
                    # When dist == 0 (alignment inside the intron) the full range is
                    # needed to back the window up to where the read body actually starts.
                    _off_limit = min(junction_proximity_bp, dist) if dist > 0 else junction_proximity_bp
                    for _off in range(_off_limit + 1):
                        _es = _eff_start - rescue_len - _off
                        if _es < 0:
                            continue
                        _cand = genome_seq[_es:_eff_start - _off].upper()
                        if len(_cand) < rescue_len:
                            continue
                        _ed = _hp_edit_distance(rescue_seq.upper(), _cand)
                        _shift_abs = abs(_shift)
                        # Two-phase scoring (lower tuple = better):
                        #   Phase 1 — discovery: minimise edit distance across wide range
                        #   Phase 2 — refinement: prefer canonical donor, then within
                        #             the ambiguity window, then smallest |shift|
                        _cur  = (not _donor_ok, not _in_amb, _shift_abs)
                        _best = (not _best_local_canonical, not _best_in_amb,
                                 _best_local_shift_abs)
                        if _ed < _best_local_ed or (
                                _ed == _best_local_ed and _cur < _best):
                            _best_local_ed = _ed
                            _best_local_canonical = _donor_ok
                            _best_in_amb = _in_amb
                            _best_local_shift_abs = _shift_abs
                            exon_seq = _cand
                            _eff_intron_start = _eff_start

                if not exon_seq:
                    continue
            else:
                # Minus strand: upstream intron (in transcript) has intron_start ≥ align_5prime.
                # When mapPacBio extends into the intron, align_5prime > intron_start
                # (dist < 0) — allow rescue if alignment ends inside the intron
                # (intron_start < align_5prime < intron_end), not completely past it.
                dist = intron_start - align_5prime
                if dist < 0:
                    if align_5prime >= intron_end:
                        continue  # alignment is downstream of intron entirely — skip
                    dist = 0  # treat as touching the boundary
                elif dist > junction_proximity_bp:
                    continue
                # Dynamic shift range for the minus-strand 5'SS boundary at intron_end:
                #   r_amb = consecutive exon bases (right of intron_end) equal to last intron base
                #   l_amb = consecutive intron bases (left of intron_end) equal to first exon base
                # Canonical 5'SS dinucleotide in genomic (minus) orientation: AC (= RC of GT).
                if 0 < intron_end <= _gs:
                    _lib = genome_seq[intron_end - 1].upper()                        # last intron base (genomic right)
                    _feb = genome_seq[intron_end].upper() if intron_end < _gs else 'N'  # first exon base
                    _r_amb = 0
                    while (_r_amb < _MAX_SS_SHIFT and intron_end + _r_amb < _gs
                           and genome_seq[intron_end + _r_amb].upper() == _lib):
                        _r_amb += 1
                    _l_amb = 0
                    while (_l_amb < _MAX_SS_SHIFT and intron_end - 1 - _l_amb >= 0
                           and genome_seq[intron_end - 1 - _l_amb].upper() == _feb):
                        _l_amb += 1
                else:
                    _r_amb = _l_amb = 0
                # Wide range for discovery (catches imprecise annotations / aligners),
                # at least ±5 bp regardless of local ambiguity.
                _shift_lo = -max(5, _l_amb)
                _shift_hi = max(5, _r_amb)

                # 3'SS acceptor quality for this junction (fixed; used as outer tiebreaker).
                # Minus strand: first 2 bases of intron (genomic) = RC of RNA-level 3'SS motif.
                _acc_di = genome_seq[intron_start:intron_start + 2].upper() if intron_start + 2 <= _gs else 'NN'
                _acceptor_priority = _ACCEPTOR_PRIORITY_MINUS.get(_acc_di, 4)

                _best_local_ed: float = rescue_len + 1
                _best_local_canonical = False
                _best_in_amb = False
                _best_local_shift_abs = max(abs(_shift_lo), _shift_hi) + 1
                exon_seq = ""
                _eff_intron_end = intron_end

                for _shift in range(_shift_lo, _shift_hi + 1):
                    _eff_end = intron_end + _shift
                    if _eff_end - 2 < 0 or _eff_end > _gs:
                        continue
                    # Canonical 5'SS on minus strand in genomic orientation:
                    # AC (RC of GT, major spliceosome) or GC (RC of GC, minor)
                    _donor_ok = genome_seq[_eff_end - 2:_eff_end].upper() in ('AC', 'GC')
                    # Whether this shift is within the natural sequence-ambiguity window
                    _in_amb = (-_l_amb <= _shift <= _r_amb)
                    _off_limit = min(junction_proximity_bp, dist) if dist > 0 else junction_proximity_bp
                    for _off in range(_off_limit + 1):
                        _cs = _eff_end + _off
                        _cand = genome_seq[_cs:_cs + rescue_len].upper()
                        if len(_cand) < rescue_len:
                            continue
                        _ed = _hp_edit_distance(rescue_seq.upper(), _cand)
                        _shift_abs = abs(_shift)
                        _cur  = (not _donor_ok, not _in_amb, _shift_abs)
                        _best = (not _best_local_canonical, not _best_in_amb,
                                 _best_local_shift_abs)
                        if _ed < _best_local_ed or (
                                _ed == _best_local_ed and _cur < _best):
                            _best_local_ed = _ed
                            _best_local_canonical = _donor_ok
                            _best_in_amb = _in_amb
                            _best_local_shift_abs = _shift_abs
                            exon_seq = _cand
                            _eff_intron_end = _eff_end

                if not exon_seq:
                    continue

            if len(exon_seq) < rescue_len:
                continue

            ed_exon = _hp_edit_distance(rescue_seq.upper(), exon_seq)
            # Compare against intronic sequence to avoid rescuing reads that match
            # the intron equally well.  Nanopore homopolymer undercalling means a
            # fixed edit-distance threshold is too strict; instead we rescue when
            # the exon match is ≥30% better than the intron match.
            if strand == '+':
                # Intronic sequence: the last rescue_len bases before the 3'SS
                # (i.e. just inside the intron from the splice-acceptor site)
                _ic_start = intron_end - rescue_len
                intron_cmp_seq = genome_seq[max(0, _ic_start):intron_end].upper()
            else:
                # Intronic sequence: the first rescue_len bases after the 3'SS
                intron_cmp_seq = genome_seq[intron_start:intron_start + rescue_len].upper()
            if len(intron_cmp_seq) == rescue_len:
                ed_intron = _hp_edit_distance(rescue_seq.upper(), intron_cmp_seq)
                # Rescue if: perfect exon match, OR exon edit-dist is ≥30% better
                # than intron edit-dist (comparative threshold, not absolute)
                rescue_ok = (ed_exon == 0) or (ed_intron > 0 and ed_exon < ed_intron * 0.70)
            else:
                # Near chromosome boundary — fall back to absolute threshold
                rescue_ok = (ed_exon / rescue_len <= max_edit_frac)
            if rescue_ok:
                # Full tiebreaking tuple (lower = better):
                #   1. edit distance (hp-aware, float)
                #   2. canonical 5'SS donor (GT/GC plus, AC/GC minus)
                #   3. shift within sequence-ambiguity window
                #   4. smallest |shift| from annotated position
                #   5. 3'SS acceptor quality: AG=0, CG=1, TG=2, AT=3, other=4
                _cur_outer  = (ed_exon, not _best_local_canonical, not _best_in_amb,
                               _best_local_shift_abs, _acceptor_priority)
                _best_outer = (best_ed, not best_is_canonical, not best_in_amb,
                               best_shift_abs, best_acceptor_priority)
                _overall_update = (best_ed < 0 or _cur_outer < _best_outer)
                if _overall_update:
                    best_ed = ed_exon
                    best_is_canonical = _best_local_canonical
                    best_in_amb = _best_in_amb
                    best_shift_abs = _best_local_shift_abs
                    best_acceptor_priority = _acceptor_priority
                    # Update 5' end using the effective donor/acceptor position
                    if strand == '+':
                        best_junction = (j_chrom, _eff_intron_start, intron_end)
                        best_five_prime_corrected = _eff_intron_start - 1
                    else:
                        best_junction = (j_chrom, intron_start, _eff_intron_end)
                        best_five_prime_corrected = _eff_intron_end

        if best_junction is not None:
            # Compute local alignment CIGAR for the exon portion so that
            # bam_writer can emit M/I/D ops instead of a flat nM block.
            #
            # IMPORTANT: use the ACTUAL intronic query bases (bases that map
            # to positions inside the intron) rather than the full rescue_seq.
            # rescue_seq may include exon-2 bases beyond the intron boundary
            # when a downstream imperfect op (e.g. an I/D inside exon 2)
            # pushed the last_imp_q cursor further than the true intron edge.
            # Passing too many bases to align_clip_to_exon produces a CIGAR
            # with more query-consuming ops than the BAM surgery step will
            # actually trim, causing the reroute sanity check to fail.
            _j_chrom, _intron_start, _intron_end = best_junction
            _clip_bd = _intron_start if strand == '-' else _intron_end
            _intronic_seq = _get_intronic_query_bases(read, _clip_bd, strand)
            _align_seq = _intronic_seq if _intronic_seq else rescue_seq
            _exon_cigar_str = ''
            try:
                from .local_aligner import align_clip_to_exon, cigar_ops_to_str
                _cigar_ops, _ = align_clip_to_exon(
                    _align_seq, genome_seq,
                    _intron_start, _intron_end, strand,
                )
                _exon_cigar_str = cigar_ops_to_str(_cigar_ops)
            except Exception as _e:
                logger.debug("Local alignment failed for read %s: %s", read.query_name, _e)
            return {
                'rescued': True,
                'rescue_type': rescue_type_candidate,
                'five_prime_corrected': best_five_prime_corrected,
                'rescued_junction': best_junction,
                'edit_distance': best_ed,
                'query_bp': len(rescue_seq),
                'five_prime_exon_cigar': _exon_cigar_str,
            }

    # --- Case 4: 5' end is strictly inside an annotated intron, no N-op for it ---
    # Fires as a fallback when:
    #   (a) no rescue_seq was extracted (all = ops, completely clean intronic match), OR
    #   (b) rescue_seq was extracted but no candidate junction produced a sequence match.
    # Even a single mismatch (1 bp rescue_seq) that fails the edit-distance check
    # should still be corrected via snap — the position is demonstrably wrong.
    # Only fires if no existing N-op in the CIGAR already covers this intron
    # (prevents double-rescue for reads that have a correct but off-by-a-few-bp N).
    _n_intervals = _get_n_op_intervals(read)
    for j_entry in candidate_junctions:
        j_chrom, intron_start, intron_end = j_entry[0], j_entry[1], j_entry[2]
        if j_chrom != chrom:
            continue
        # Condition: align_5prime inside [intron_start, intron_end).
        # intron_start is inclusive: a read whose rightmost base IS intron_start
        # (reference_end = intron_start + 1 for minus strand) is mapping into
        # the intron and must be snapped.
        if not (intron_start <= align_5prime < intron_end):
            continue
        # Skip if an existing N-op already approximates this intron
        already_has_n = any(
            abs(ns - intron_start) <= junction_proximity_bp and
            abs(ne - intron_end) <= junction_proximity_bp
            for ns, ne in _n_intervals
        )
        if already_has_n:
            continue
        # Snap to exon-1-side boundary
        snap_pos = intron_end if strand == '-' else intron_start
        intronic_depth = (align_5prime - intron_start) if strand == '-' else (intron_end - align_5prime)
        # Compute exon CIGAR so bam_writer can reroute the intronic tail.
        _clip_bd4 = intron_start if strand == '-' else intron_end
        _intronic_seq4 = _get_intronic_query_bases(read, _clip_bd4, strand)
        _exon_cigar_str4 = ''
        if _intronic_seq4:
            try:
                from .local_aligner import align_clip_to_exon, cigar_ops_to_str
                _cigar_ops4, _ = align_clip_to_exon(
                    _intronic_seq4, genome_seq, intron_start, intron_end, strand,
                )
                _exon_cigar_str4 = cigar_ops_to_str(_cigar_ops4)
            except Exception as _e4:
                logger.debug("Case 4 local alignment failed for read %s: %s",
                             read.query_name, _e4)
        return {
            'rescued': True,
            'rescue_type': 'intronic_snap',
            'five_prime_corrected': snap_pos,
            'rescued_junction': (j_chrom, intron_start, intron_end),
            'edit_distance': -1,
            'query_bp': intronic_depth,
            'five_prime_exon_cigar': _exon_cigar_str4,
        }

    # --- Case 3: proximity-only (no sequence to match, but start is at a 3'SS) ---
    for j_entry in candidate_junctions:
        j_chrom, intron_start, intron_end = j_entry[0], j_entry[1], j_entry[2]
        if j_chrom != chrom:
            continue
        if strand == '+':
            dist = align_5prime - intron_end
        else:
            dist = intron_start - align_5prime
        if 0 <= dist <= junction_proximity_bp:
            return {
                'rescued': False,
                'rescue_type': 'proximity',
                'five_prime_corrected': align_5prime,  # no change; no evidence
                'rescued_junction': (j_chrom, intron_start, intron_end),
                'edit_distance': -1,
                'query_bp': 0,
                'five_prime_exon_cigar': '',
            }

    return _no_rescue(read, strand)


def _no_rescue(read: pysam.AlignedSegment, strand: str) -> Dict:
    """Return a no-rescue result dict."""
    if strand == '-':
        _ref_end = read.reference_end
        five_prime = (_ref_end - 1) if _ref_end is not None else None
    else:
        five_prime = read.reference_start
    return {
        'rescued': False,
        'rescue_type': 'none',
        'five_prime_corrected': five_prime,
        'rescued_junction': None,
        'edit_distance': -1,
        'query_bp': 0,
        'five_prime_exon_cigar': '',
    }


# =============================================================================
# Validation
# =============================================================================

def validate_5prime_correction():
    """
    Print validation examples for 5' end correction.
    """
    print("=" * 60)
    print("SPLICE-AWARE 5' END CORRECTION VALIDATION")
    print("=" * 60)

    # Create mock read for testing
    class MockRead:
        def __init__(self, start, end, is_reverse, reference_name='chrI'):
            self.reference_start = start
            self.reference_end = end
            self.is_reverse = is_reverse
            self.reference_name = reference_name

    # Test case 1: Plus strand, 5' end in junction
    print("\nTest 1: Plus strand read with 5' end in junction")
    read1 = MockRead(start=1050, end=2000, is_reverse=False)
    # Junction from 1000-1100, so 5' end (1050) is in the junction
    junctions1 = [(1000, 1100)]
    result1 = correct_5prime_for_splicing(read1, read_junctions=junctions1)
    print(f"  Raw 5' end: {result1.five_prime_raw}")
    print(f"  Corrected 5' end: {result1.five_prime_corrected}")
    print(f"  Starts in intron: {result1.starts_in_intron}")
    print(f"  Correction: {result1.correction_bp} bp")
    print(f"  Reason: {result1.correction_reason}")
    if result1.starts_in_intron and result1.five_prime_corrected == 999:
        print("  PASSED")
    else:
        print("  FAILED")

    # Test case 2: Plus strand, 5' end NOT in junction
    print("\nTest 2: Plus strand read with 5' end NOT in junction")
    read2 = MockRead(start=900, end=2000, is_reverse=False)
    # Junction from 1000-1100, so 5' end (900) is before the junction
    junctions2 = [(1000, 1100)]
    result2 = correct_5prime_for_splicing(read2, read_junctions=junctions2)
    print(f"  Raw 5' end: {result2.five_prime_raw}")
    print(f"  Corrected 5' end: {result2.five_prime_corrected}")
    print(f"  Starts in intron: {result2.starts_in_intron}")
    if not result2.starts_in_intron and result2.five_prime_corrected == result2.five_prime_raw:
        print("  PASSED")
    else:
        print("  FAILED")

    # Test case 3: Minus strand, 5' end in junction
    print("\nTest 3: Minus strand read with 5' end in junction")
    read3 = MockRead(start=1000, end=2050, is_reverse=True)
    # Junction from 2000-2100, so 5' end (2050-1=2049) is in the junction
    junctions3 = [(2000, 2100)]
    result3 = correct_5prime_for_splicing(read3, read_junctions=junctions3)
    print(f"  Raw 5' end: {result3.five_prime_raw}")
    print(f"  Corrected 5' end: {result3.five_prime_corrected}")
    print(f"  Starts in intron: {result3.starts_in_intron}")
    print(f"  Correction: {result3.correction_bp} bp")
    if result3.starts_in_intron and result3.five_prime_corrected == 2100:
        print("  PASSED")
    else:
        print("  FAILED")

    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    validate_5prime_correction()
