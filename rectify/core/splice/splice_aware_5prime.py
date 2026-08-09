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
import array as _array_mod
import os
import pysam

from ...utils.genome import standardize_chrom_name
from ...config import CHROM_TO_GENOME
from .overhang_informativeness import (
    COUNTERS as _OI_COUNTERS,
    assess_overhang as _oi_assess,
    gate_alpha as _oi_gate_alpha,
    gate_enabled as _oi_gate_enabled,
)

try:
    from intervaltree import IntervalTree
    HAS_INTERVALTREE = True
except ImportError:
    HAS_INTERVALTREE = False


# =============================================================================
# Rescue tuning constants
# =============================================================================
# Module-level so the bam_processor pool-fetch window can be derived from the
# SAME values the rescue actually uses. The old _POOL_SEARCH_RADIUS was a magic
# 10000 sized for a junction_proximity_bp default of 5000 that no longer exists
# (real default 10) — exactly the vestigial-constant drift dev/PERF_AUDIT.md
# warns about. Keep these as the single source of truth.
MAX_SS_SHIFT = 15                  # max splice-site slide explored per candidate
DEFAULT_JUNCTION_PROXIMITY_BP = 10  # max bp between 5' align end and intron edge
DEFAULT_PEEL_MAX_BP = 100          # deepest terminal peel attempted


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
    :func:`~rectify.core.bam.read_edits.reroute_intronic_tail_5prime_via_junction`
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
            if op == 3:   # N op
                # A false junction entirely above clip_boundary (both ends inside
                # the intron) should be passed through rather than stopping — the
                # aligner placed a spurious small junction inside the target intron.
                # Only stop if the N spans below clip_boundary (real exon boundary).
                new_ref = ref_pos - length
                if new_ref >= clip_boundary:
                    ref_pos = new_ref   # pass through: N inside intron
                    continue
                else:
                    break               # N crosses exon/intron boundary — stop
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
                # Pass through false junctions entirely below clip_boundary
                new_ref = ref_pos + length
                if new_ref <= clip_boundary:
                    ref_pos = new_ref
                    continue
                else:
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


# ---------------------------------------------------------------------------
# Pre-allocated scratch buffers for _hp_edit_distance.
#
# Python's arena allocator never returns claimed arenas to the OS.  Each call
# to _hp_edit_distance previously allocated a 2D list-of-lists DP table
# (~1 MB for 192×192 sequences) plus two cost lists — all GC'd between calls
# but leaving arenas partially occupied.  With ~20% of reads triggering
# rescue_3ss_truncation × 6k calls/read × 37k floats/call, the arena pool
# grows to tens of GB and never shrinks.
#
# Fix: single module-level flat array.array buffers reused on every call.
# No Python heap allocation inside _hp_edit_distance at all.
# ---------------------------------------------------------------------------
_HP_ED_MAX_LEN = 200          # sequences longer than this are truncated
_HP_ED_STRIDE  = _HP_ED_MAX_LEN + 1          # row stride (m+1 dimension)
# Perf cap for the 3'SS-rescue DP (PERF_AUDIT.md rec #2, large_5prime_clip class):
# only the bases ADJACENT to the splice donor determine the boundary; distal
# exon-1 bases add _hp_edit_distance cost (up to 200x200 per call) without
# improving boundary precision. The rescue scoring sequence is capped to this
# many donor-adjacent bases. Surgery (align_clip_to_exon) still uses the full
# clip, so CIGAR geometry is unaffected — only junction *selection* is capped.
#
# Value 100 (not the 60 the audit sketched): an empirical pre-cap vs post-cap diff
# on real upf1Δ DRS (chrII+chrIV minimap2; 217 rescued >60bp clips) found cap=60
# produced one 1-bp donor slide where distal bases broke an ED-tie toward the
# canonical AC donor; cap>=80 reproduced the uncapped result exactly (0 diffs).
# 100 keeps full rescue-outcome equivalence with margin while still shrinking the
# DP ~4x (100x100 vs 200x200) on the median 244bp archetype clip. Lower toward 60
# if perf demands it and the 1-bp slide is acceptable.
_RESCUE_DP_CAP = 100
_hp_ed_dp  = _array_mod.array('d', [0.0] * _HP_ED_STRIDE * _HP_ED_STRIDE)
_hp_ed_del = _array_mod.array('d', [0.0] * _HP_ED_MAX_LEN)
_hp_ed_ins = _array_mod.array('d', [0.0] * _HP_ED_MAX_LEN)


def _hp_edit_distance(s1: str, s2: str) -> float:
    """Edit distance with 0.5 penalty for indels within homopolymer runs.

    Nanopore sequencers under/over-call homopolymer run lengths.  A deletion
    or insertion of a base that is part of a run of identical bases (i.e. the
    indel base matches its immediate neighbour in the same sequence) is given
    half the normal gap penalty (0.5 instead of 1.0).  Substitutions always
    cost 1.0 regardless of context.

    A base is considered part of a homopolymer run if it equals either the
    preceding or the following character in the *same* string.

    Empirical calibration note
    --------------------------
    The 0.5 HP indel penalty is a heuristic step-function (HP=1→1.0,
    HP≥2→0.5).  Empirical data from ``empirical_cigar_error_profiler.py``
    shows base-class- and length-specific costs that differ substantially:

        HP=8 A/T del cost ≈ 0.13   (this function gives 0.5)
        HP=8 C/G del cost ≈ 0.77   (this function gives 0.5)
        HP=1 non-STR del cost ≈ 0.58  (this function gives 1.0)

    To use empirical penalties here, pass a ``junction_refiner.HpPenaltyTable``
    and replace ``_del_cost``/``_ins_cost`` with
    ``table.del_cost(hp_run_length(s1, i-1), s1[i-1])`` and similarly for
    insertions.  The 3'SS rescue context (exon sequence vs intron boundary)
    makes the local HP run length in ``s1`` the relevant quantity, not the
    genome-wide run length.  This upgrade has not been implemented because
    the 3'SS rescue path (Module Cat3) is invoked with short sequences where
    the heuristic is adequate; revisit if Cat3 false-positives appear at
    long A/T HP exon-ends.
    """
    # Truncate to the pre-allocated buffer size (splice-site rescue does not
    # benefit from more than ~200 bp of context).
    n = len(s1)
    m = len(s2)
    if n > _HP_ED_MAX_LEN:
        n = _HP_ED_MAX_LEN
        s1 = s1[:n]
    if m > _HP_ED_MAX_LEN:
        m = _HP_ED_MAX_LEN
        s2 = s2[:m]
    if n == 0:
        return float(m)
    if m == 0:
        return float(n)

    # Fill HP-cost arrays in-place (no allocation).
    # _hp_ed_del[i-1] = cost to delete s1[i-1] (0.5 if in HP run, else 1.0).
    # _hp_ed_ins[j-1] = cost to insert s2[j-1].
    for i in range(1, n + 1):
        _hp_ed_del[i - 1] = (
            0.5 if (i >= 2 and s1[i - 2] == s1[i - 1]) or (i < n and s1[i - 1] == s1[i])
            else 1.0
        )
    for j in range(1, m + 1):
        _hp_ed_ins[j - 1] = (
            0.5 if (j >= 2 and s2[j - 2] == s2[j - 1]) or (j < m and s2[j - 1] == s2[j])
            else 1.0
        )

    # 2-D DP using flat pre-allocated array.  Index: dp[i][j] = _hp_ed_dp[i*S + j].
    _S = _HP_ED_STRIDE
    _dp = _hp_ed_dp
    _dp[0] = 0.0
    for j in range(1, m + 1):
        _dp[j] = _dp[j - 1] + _hp_ed_ins[j - 1]
    for i in range(1, n + 1):
        _row  = i * _S
        _prev = _row - _S
        dc    = _hp_ed_del[i - 1]
        s1c   = s1[i - 1]
        _dp[_row] = _dp[_prev] + dc
        for j in range(1, m + 1):
            if s1c == s2[j - 1]:
                _dp[_row + j] = _dp[_prev + j - 1]
            else:
                _dp[_row + j] = min(
                    _dp[_prev + j - 1] + 1.0,        # substitution
                    _dp[_prev + j]     + dc,           # deletion
                    _dp[_row  + j - 1] + _hp_ed_ins[j - 1],  # insertion
                )
    return _dp[n * _S + m]


# 3'SS acceptor dinucleotide priority (lower = more canonical).
# Plus strand:  last 2 bases of intron before exon 2 (genome[intron_end-2:intron_end])
# Minus strand: first 2 bases of intron (genome[intron_start:intron_start+2]),
#               which is the RC of the RNA-level 3'SS motif.
_ACCEPTOR_PRIORITY_PLUS  = {'AG': 0, 'CG': 1, 'TG': 2, 'AT': 3}
_ACCEPTOR_PRIORITY_MINUS = {'CT': 0, 'CG': 1, 'CA': 2, 'AT': 3}


def _log_peel_screen(path, read, strand, chrom, baseline, peel, depth) -> None:
    """Append one screening row (case-1 / new-rescue candidate) to ``path``.

    Best-effort and exception-swallowing: screening must NEVER affect correctness
    or crash a production run. Columns let candidates be categorised offline:
    a row with baseline_rescued=1 and a differing peel junction is a case-1
    override candidate (the deferred prior-aware override); baseline_rescued=0 is
    a new-rescue candidate (what fill-only adds).
    """
    try:
        def _j(d):
            rj = d.get('rescued_junction')
            return f"{rj[1]}-{rj[2]}" if rj else "NA"
        new = (not os.path.exists(path)) or os.path.getsize(path) == 0
        with open(path, 'a') as fh:
            if new:
                fh.write("read_id\tstrand\tchrom\tbaseline_rescued\tbaseline_junction\t"
                         "baseline_5p\tbaseline_ed\tpeel_junction\tpeel_5p\tpeel_ed\t"
                         "peel_depth\tpeel_rescue_type\n")
            fh.write("\t".join(str(x) for x in [
                read.query_name, strand, chrom, int(bool(baseline.get('rescued'))),
                _j(baseline), baseline.get('five_prime_corrected', 'NA'),
                baseline.get('edit_distance', 'NA'),
                _j(peel), peel.get('five_prime_corrected', 'NA'),
                peel.get('edit_distance', 'NA'), depth, peel.get('rescue_type', 'NA'),
            ]) + "\n")
    except Exception:
        pass


def _peel_candidate_depths(
    read: pysam.AlignedSegment,
    genome_seq: str,
    strand: str,
    max_peel: int,
    clean_anchor: int,
) -> List[int]:
    """Generate candidate terminal-peel depths (counts of 5'-terminal query bases).

    Walk the alignment inward from the read's 5' end. Each query-consuming
    position yields a candidate depth (number of query bases peeled from the 5'
    end so far). Stop generating at the first of:
      - an existing N op (a real/strong splice gap; detected as a ref-position
        jump > 1 between consecutive aligned reference bases),
      - ``clean_anchor`` consecutive clean matches (read base == genome base) —
        beyond this the read is well-anchored and deeper peels are pointless,
      - ``max_peel`` query bases.

    Returns depths in ascending order. Clean-match inference works for M-only
    CIGARs (mapPacBio) because it compares query vs genome directly rather than
    relying on =/X ops.
    """
    if not read.cigartuples or not read.query_sequence:
        return []
    try:
        pairs = read.get_aligned_pairs()
    except Exception:
        return []
    if read.is_reverse:
        pairs = list(reversed(pairs))
    # N ops (existing splice gaps) detected from the CIGAR, not from aligned-pair
    # ref jumps: pysam emits per-base (None, ref) pairs for N just like D, so a
    # ref-jump heuristic cannot distinguish them.
    n_intervals = _get_n_op_intervals(read)
    query_seq = read.query_sequence
    gs = len(genome_seq)
    depths: List[int] = []
    q_from_5p = 0
    consec_clean = 0
    for qp, rp in pairs:
        # Stop at the first existing N op — do not peel across an annotated /
        # aligner-called intron.
        if rp is not None and any(s <= rp < e for s, e in n_intervals):
            break
        if qp is None:
            continue  # D op: consumes reference only
        q_from_5p += 1
        if rp is not None and 0 <= rp < gs:
            gb = genome_seq[rp].upper()
            rb = query_seq[qp].upper()
            if gb != 'N' and rb == gb:
                consec_clean += 1
            else:
                consec_clean = 0
        else:
            consec_clean = 0  # I/S position: not a clean match
        depths.append(q_from_5p)
        if consec_clean >= clean_anchor or q_from_5p >= max_peel:
            break
    return depths


def _terminal_peel_rescue(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int]],
    strand: str,
    max_edit_frac: float,
    junction_proximity_bp: int,
    scan_bp: int,
    chrom: str,
    genome_seq: str,
    baseline: Dict,
    peel_max_bp: int,
    peel_clean_anchor: int,
    accept_margin: float,
) -> Optional[Dict]:
    """Multi-hypothesis terminal peel (Module 2F).

    Score several 5'-terminal peel depths against the candidate junctions and
    return a *strictly better* hypothesis than ``baseline`` (the default
    single-depth rescue result), or ``None`` to keep the baseline.

    Monotonic by construction: the baseline is always retained unless a deeper
    peel produces a confident sequence rescue with a strictly lower normalised
    HP-edit-distance (beyond ``accept_margin``). Ties favour the baseline /
    no-rescue, per the spec. Default-on behaviour therefore cannot regress below
    the pre-peel rescue by the edit-distance metric.

    Returns a dict with an added ``'terminal_peel'`` marker when it overrides.
    """
    # Screening (dry-run) mode: when RECTIFY_PEEL_SCREEN names a log path, run the
    # comparison even for already-rescued reads to surface case-1 override
    # candidates, but DO NOT change the (fill-only) production output.
    screen_path = os.environ.get('RECTIFY_PEEL_SCREEN')
    base_rescued = bool(baseline.get('rescued'))
    # Production fill-only fast path: never override an existing rescue.
    if base_rescued and not screen_path:
        return None

    # Gate A — terminal boundary evidence near the 5' end (S/X/I/D). No evidence
    # => clean terminal alignment => do not peel.
    if not _extract_5prime_rescue_seq(read, genome_seq=genome_seq, scan_ref_bp=scan_bp):
        return None

    # 5' alignment boundary (first aligned base in transcript orientation).
    if strand == '+':
        align_5prime = read.reference_start
    else:
        if read.reference_end is None:
            return None
        align_5prime = read.reference_end - 1

    # Gate B — collect candidate 3'SS within reach of a peel of <= peel_max_bp.
    # We narrow to this nearby subset and pass ONLY it to the per-depth body calls
    # below (the body scans its candidate set 3x per call with no internal
    # narrowing). The full junction pool can be ~17k entries in run-all mode, so
    # scanning it once per peel-depth per read is what blew up to a multi-hour OOM
    # on high-coverage loci (2026-05-24). A 5'-terminal peel re-aligns a segment
    # of <= peel_max_bp to an upstream exon, so a 3'SS beyond `reach` can never
    # win — excluding distant junctions is safe (does not drop a real winner).
    reach = junction_proximity_bp + peel_max_bp
    _nearby_with_dist: List[Tuple[int, Tuple]] = []
    for j in candidate_junctions:
        if j[0] != chrom:
            continue
        if len(j) >= 4 and j[3] not in (strand, '.', ''):
            continue
        edge = j[2] if strand == '+' else j[1]   # intron_end (+) / intron_start (-)
        dist = abs(align_5prime - edge)
        if dist <= reach:
            _nearby_with_dist.append((dist, j))
    if not _nearby_with_dist:
        return None

    # Cap to K closest junctions before the per-depth body calls.
    # Each depth calls _rescue_3ss_truncation_body which runs _hp_edit_distance
    # O(depths × N_junctions × N_shifts) times.  On junction-dense loci (e.g.
    # human chr5 SMN1/SMN2) the uncapped set can contain 50-200+ entries within
    # `reach`, making this O(100 × 200 × 40) = 800k DP calls per read.
    # Keeping only the K closest is provably non-regressing: a peel aligns the
    # 5' boundary to a 3'SS edge; among all nearby edges the closest ones are
    # the most likely geometric match.  The body's own proximity window further
    # narrows the set anyway, so distant junctions almost never win.
    _MAX_NEARBY_JUNCTIONS = 20
    if len(_nearby_with_dist) > _MAX_NEARBY_JUNCTIONS:
        # Coordinate is part of the key so equidistant candidates are cut
        # reproducibly: the input order comes from iterating a SET, so a
        # distance-only sort leaves the surviving slice PYTHONHASHSEED-dependent
        # (planning/649).
        _nearby_with_dist.sort(key=lambda x: (x[0], x[1][0], x[1][1], x[1][2]))
        _nearby_with_dist = _nearby_with_dist[:_MAX_NEARBY_JUNCTIONS]
    nearby_junctions: Set[Tuple] = {jd[1] for jd in _nearby_with_dist}
    max_edge_dist = max(jd[0] for jd in _nearby_with_dist)

    # Acceptance baseline. A peel must beat the baseline's normalised HP-edit
    # distance to win. When the baseline rescued via sequence (ed >= 0), that ed
    # is the bar (so screening surfaces only genuine case-1 candidates); when the
    # baseline did not rescue (or snapped, ed < 0), accept any confident peel.
    base_ed = baseline.get('edit_distance', -1.0)
    base_qbp = baseline.get('query_bp', 0) or 0
    base_norm = (base_ed / max(1, base_qbp)) if (base_rescued and base_ed is not None and base_ed >= 0) else float('inf')

    query_seq = read.query_sequence or ''
    n_query = len(query_seq)
    if n_query == 0:
        return None

    # Cap the peel depth at the farthest reachable 3'SS edge (+ the body's
    # boundary-ambiguity slack). A winning peel must align the 5' boundary onto a
    # nearby junction edge, so no winner lands deeper than `max_edge_dist`; deeper
    # depths are pure waste. On poorly-5'-anchored reads (soft-clip / mapPacBio
    # forced-mismatch) the depth generator would otherwise run all the way to
    # peel_max_bp=100, costing ~100 DP alignments/read — the per-read half of the
    # 2026-05-24 hang. The cap shortens the depth list to a strict prefix, so the
    # same winner is still found.
    _PEEL_EDGE_SLACK = MAX_SS_SHIFT  # matches the body's junction-slide window
    effective_max_peel = min(peel_max_bp, max_edge_dist + _PEEL_EDGE_SLACK)
    depths = _peel_candidate_depths(read, genome_seq, strand, effective_max_peel, peel_clean_anchor)
    if not depths:
        return None

    best_peel: Optional[Dict] = None
    best_depth: Optional[int] = None
    best_peel_norm = base_norm
    for d in depths:
        if d <= 0 or d > n_query:
            continue
        peeled = query_seq[:d] if strand == '+' else query_seq[n_query - d:]
        if not peeled:
            continue
        res = _rescue_3ss_truncation_body(
            read, genome, nearby_junctions, strand,
            max_edit_frac, junction_proximity_bp, scan_bp,
            chrom, genome_seq, rescue_seq_override=peeled,
        )
        if not res.get('rescued'):
            continue
        ed = res.get('edit_distance', -1.0)
        if ed is None or ed < 0:
            continue  # snap-type result from override path — ignore
        qbp = res.get('query_bp', 0) or 0
        norm = ed / max(1, qbp)
        # Strictly better than the current best (baseline or a prior peel),
        # beyond the acceptance margin. Ties favour the incumbent (baseline).
        if norm < best_peel_norm - accept_margin:
            best_peel_norm = norm
            best_peel = res
            best_depth = d

    if best_peel is None:
        return None
    # Screening: log the candidate (case-1 if baseline rescued, else new-rescue)
    # without changing output.
    if screen_path:
        _log_peel_screen(screen_path, read, strand, chrom, baseline, best_peel, best_depth)
    # Production fill-only: never override an existing baseline rescue.
    if base_rescued:
        return None
    best_peel = dict(best_peel)
    best_peel['terminal_peel'] = True
    return best_peel


def rescue_3ss_truncation(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int]],
    strand: Optional[str] = None,
    max_edit_frac: float = 0.2,
    junction_proximity_bp: int = DEFAULT_JUNCTION_PROXIMITY_BP,
    scan_bp: int = 50,
    terminal_peel: bool = True,
    peel_max_bp: int = DEFAULT_PEEL_MAX_BP,
    peel_clean_anchor: int = 10,
    peel_accept_margin: float = 0.0,
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

    4. **Intronic snap** (Case 4): if sequence-based rescue produced no match but
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

    # Short-circuit: a 5' soft-clip that is itself a long low-period repeat
    # expansion ((AAG)n / (CTT)n etc.) is a Nanopore DRS basecaller artifact
    # (motor slippage at RNA secondary structure), NOT a missed upstream intron.
    # Rescuing it is wasted work — the low-complexity clip is compared against
    # every candidate junction in the pool (~10k+ on human) and matches many
    # spuriously — and any "rescue" produced is wrong. On human RNA004 DRS ~90% of
    # reads entering rescue with a >=30 bp 5' soft-clip are these artifacts (clips
    # up to ~200 kb). Real Cat3 rescues have a diverse exon-1 5' clip that does NOT
    # trigger this, so they are unaffected. See core/splice/repeat_expansion.py.
    from .repeat_expansion import is_repeat_expansion
    _cigar = read.cigartuples
    _seq = read.query_sequence
    if _cigar and _seq:
        if strand == '+':
            _clip5 = _seq[:_cigar[0][1]] if _cigar[0][0] == 4 else ""
        else:
            _clip5 = _seq[-_cigar[-1][1]:] if _cigar[-1][0] == 4 else ""
        if is_repeat_expansion(_clip5):
            _res = _no_rescue(read, strand)
            _res['repeat_expansion'] = True
            return _res

    # 5'-edge reanchor pre-pass (for mapPacBio-style reads whose 5' edge has a
    # tight X/I/D cluster blocking soft-clip-based rescue). reanchor_5prime_for_rescue
    # walks from the 5' edge to find the first sustained match run (≥10 bp) and
    # collapses everything upstream into a leading soft-clip — making the read
    # look like a normal cat3 soft-clipped read to the rest of this function.
    #
    # We mutate the read in-place to keep the rest of the function unchanged,
    # then restore via try/finally so callers do not see the mutation. The same
    # reanchor will be reapplied deterministically in bam_writer (gated on the
    # reanchor_clip_len value we emit in the result dict).
    from ..bam.read_edits import reanchor_5prime_for_rescue as _reanchor_5p
    _saved_cigartuples = read.cigartuples
    _saved_reference_start = read.reference_start
    _cigar_before = list(_saved_cigartuples) if _saved_cigartuples else []
    _reanchor_clip_len = 0
    try:
        if _reanchor_5p(read, genome, anchor_min_run=10):
            _cigar_after = list(read.cigartuples) if read.cigartuples else []
            # Only treat as a "material" reanchor when the CIGAR actually changed
            # shape (a no-mutation reanchor — e.g. read already had a leading S
            # matching the first ≥10 match-run — produces an identical cigartuples
            # list and must not propagate a phantom reanchor_clip_len).
            if _cigar_after != _cigar_before:
                if strand == '+':
                    if _cigar_after and _cigar_after[0][0] == 4:
                        _reanchor_clip_len = _cigar_after[0][1]
                else:
                    if _cigar_after and _cigar_after[-1][0] == 4:
                        _reanchor_clip_len = _cigar_after[-1][1]
        # If reanchor did not materially change the CIGAR, restore now so the
        # rest of the function (and the finally restore) is a no-op on read state.
        if _reanchor_clip_len == 0:
            read.cigartuples = _saved_cigartuples
            read.reference_start = _saved_reference_start

        _result = _rescue_3ss_truncation_body(
            read, genome, candidate_junctions, strand,
            max_edit_frac, junction_proximity_bp, scan_bp,
            chrom, genome_seq,
        )
        # Module 2F: multi-hypothesis terminal peel. Monotonic — only overrides
        # the baseline when a deeper peel is a strictly-better sequence rescue.
        if terminal_peel:
            _peeled = _terminal_peel_rescue(
                read, genome, candidate_junctions, strand,
                max_edit_frac, junction_proximity_bp, scan_bp,
                chrom, genome_seq, baseline=_result,
                peel_max_bp=peel_max_bp, peel_clean_anchor=peel_clean_anchor,
                accept_margin=peel_accept_margin,
            )
            if _peeled is not None:
                _result = _peeled
    finally:
        # Restore original read state so callers see no mutation.
        read.cigartuples = _saved_cigartuples
        read.reference_start = _saved_reference_start

    if _result.get('rescued'):
        _result['reanchor_clip_len'] = _reanchor_clip_len
    return _result


def _rescue_3ss_truncation_body(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int]],
    strand: str,
    max_edit_frac: float,
    junction_proximity_bp: int,
    scan_bp: int,
    chrom: str,
    genome_seq: str,
    rescue_seq_override: Optional[str] = None,
) -> Dict:
    """Inner body of rescue_3ss_truncation — see that function's docstring.

    ``rescue_seq_override``: when not None, use this sequence as the 5' rescue
    sequence instead of extracting it via ``_extract_5prime_rescue_seq``. Used by
    the multi-hypothesis terminal-peel driver to score a specific peel depth.
    When None (the default) behaviour is byte-identical to the pre-peel code.

    Split out so the outer function can wrap the body in a try/finally that
    restores the read's CIGAR/reference_start after an optional reanchor pre-pass.
    All early returns in this function operate on a potentially-reanchored read;
    the caller is responsible for restore + emitting reanchor_clip_len.
    """
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
    if rescue_seq_override is not None:
        rescue_seq = rescue_seq_override
    else:
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
    #
    # IMPORTANT: for minus-strand reads, _extract_5prime_rescue_seq returns
    # query_seq[n - last_imp_q:] — the RIGHTMOST last_imp_q bases of query_seq.
    # The soft-clip occupies the RIGHTMOST five_clip of those bases (it is the
    # 5'-end of the RNA = trailing end of the BAM query sequence for is_reverse).
    # Truncating with [:five_clip] would take the LEFTMOST bases (aligned exon-2
    # body) instead.  Use [-five_clip:] for minus strand so the slice selects the
    # actual soft-clipped sequence.
    # Skip this truncation when an explicit override is supplied: the
    # terminal-peel driver passes a deliberately-deeper peel sequence and must
    # not have it clipped back to the soft-clip length.
    if rescue_seq_override is None and five_clip > 0 and rescue_seq and len(rescue_seq) > five_clip:
        if strand == '-':
            rescue_seq = rescue_seq[-five_clip:]
        else:
            rescue_seq = rescue_seq[:five_clip]

    # --- Informativeness gate (planning/641 §2; DARK by default) -----------
    # RECTIFY_OVERHANG_INFO_GATE=1 enables. A rescue sequence whose effective
    # information content cannot distinguish its true placement from chance in
    # any window (poly(A), homopolymer, periodic repeat) has its SEQUENCE
    # search refused — a first-class outcome, not a failure. The structural
    # evidence paths below (N-op match, intronic snap, proximity-only) do not
    # depend on the sequence and stay live. When the sequence is informative
    # but weakly so, W_max additionally tightens the candidate-narrowing
    # window (never the N-op-matched branch).
    _oi_w_max: Optional[int] = None
    if rescue_seq and _oi_gate_enabled():
        _oi_assessment = _oi_assess(rescue_seq, alpha=_oi_gate_alpha())
        if _oi_assessment.refused:
            rescue_seq = ''
            rescue_type_candidate = 'none'
        else:
            _oi_w_max = _oi_assessment.w_max_bp

    # --- Try sequence-based rescue against each candidate junction ---
    best_ed: float = -1.0
    best_junction = None
    best_five_prime_corrected = align_5prime
    best_is_canonical = False    # tiebreaker 1: canonical GT/GC donor
    best_in_amb = False          # tiebreaker 2: shift within ambiguity window
    best_shift_abs = 999         # tiebreaker 3: smallest |shift|
    best_acceptor_priority = 4   # tiebreaker 4: 3'SS quality (AG=0..AT=3..other=4)

    # Forced-snap fallback for the mapPacBio terminal-D overshoot pattern.
    # When the distance gate detects _n_match AND _leading_del (N-op proves the
    # read already spans the intron, but terminal D ops push align_5prime past
    # intron_end), we record the candidate junction here.  If all sequence-based
    # rescue attempts fail (garbled post-N query bases mismatch exon-1), we snap
    # directly to intron_end/-start based solely on the N-op evidence.
    _forced_snap_junction: Optional[Tuple[str, int, int]] = None
    _forced_snap_n_err: float = float('inf')  # N-op boundary error for best candidate

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
    _MAX_SS_SHIFT = MAX_SS_SHIFT  # hard cap regardless of run length

    # Pre-compute N-op intervals from the CIGAR.  Used in two places:
    # 1. Distance gate (below): when a read already has an N-op that matches a
    #    candidate junction, terminal D/I ops on the exon-1 side can push
    #    align_5prime past intron_end and trigger the "downstream of intron"
    #    skip.  We suppress that skip for reads whose CIGAR proves they span
    #    the intron (mapPacBio 10I1M15D pattern, Bug 2 fix).
    # 2. Case 4 intronic-snap (below): already used here.
    _n_intervals = _get_n_op_intervals(read)

    # --- Candidate narrowing (perf, 2026-05-24) ---------------------------------
    # In run-all mode candidate_junctions can hold ~17k pool junctions within
    # +/-10kb of the read; on dense alt-splice loci, running the per-candidate
    # HP-edit-distance DP in the rescue loops below over all of them was ~87% of
    # CPU (py-spy) and drove the inline-correct stall. Pre-filter ONCE here to the
    # junctions any loop could actually act on; the three loops below iterate this
    # narrowed list instead of the full pool.
    #
    # Provably non-regressing: a junction is kept iff it (a) overlaps a window
    # around align_5prime wide enough to cover loop 1's distance gate
    # (dist <= junction_proximity_bp + five_clip) plus the +/-_MAX_SS_SHIFT slide,
    # Case 4's containment, and Case 3's proximity; OR (b) matches one of the
    # read's existing N-ops within junction_proximity_bp (the mapPacBio
    # _leading_del N-op-match path, where the spanned intron can sit farther from
    # align_5prime). Anything dropped here would be skipped by the cheap gate at
    # the top of every loop anyway, so no winning candidate is excluded.
    _nb_W_full = junction_proximity_bp + five_clip + _MAX_SS_SHIFT + 5
    _nb_W = _nb_W_full
    if _oi_w_max is not None:
        # Informativeness gate: the rescue sequence cannot distinguish its
        # placement beyond W_max, so candidates farther out are chance hits
        # by construction. N-op-matched junctions (the OR-branch below) are
        # CIGAR evidence, not a sequence search, and keep the full window.
        _nb_W = min(_nb_W_full, _oi_w_max + _MAX_SS_SHIFT + 5)
    _nb_W2 = junction_proximity_bp + 5
    _nearby_junctions = []
    for _j in candidate_junctions:
        if _j[0] != chrom:
            continue
        if len(_j) >= 4 and _j[3] not in (strand, '.', ''):
            continue
        if (_j[1] <= align_5prime + _nb_W and _j[2] >= align_5prime - _nb_W) or any(
            abs(_j[1] - _ns) <= _nb_W2 and abs(_j[2] - _ne) <= _nb_W2
            for _ns, _ne in _n_intervals
        ):
            _nearby_junctions.append(_j)
        elif _nb_W < _nb_W_full and (
            _j[1] <= align_5prime + _nb_W_full and _j[2] >= align_5prime - _nb_W_full
        ):
            _OI_COUNTERS['candidates_skipped'] += 1

    # Deterministic candidate order (planning/649).  ``candidate_junctions`` is a
    # SET, so the iteration above inherits PYTHONHASHSEED — and equal-edit-distance
    # candidates are resolved by ENCOUNTER order once the 4-level tiebreaker at the
    # bottom of the loop runs out.  Measured: two runs of the same tree over the
    # same BAM and pool differed in 3/474 corrected rows (all five_prime_rescued).
    # Sorting by coordinate makes the winner reproducible and makes byte-identity
    # gates (planning/596) meaningful; it also fixes the distance-cap slice below,
    # whose ``sorted``/``sort`` are stable and therefore only as deterministic as
    # their input order.
    _nearby_junctions.sort(key=lambda _j: (_j[0], _j[1], _j[2]))

    # Cap to K closest junctions when the set is large.  On junction-dense loci
    # (human chr5 SMN1/SMN2) the proximity filter can still admit 50-200+ entries.
    # The rescue loops below call _hp_edit_distance O(N × N_shifts) per read; with
    # N=200+ this is ~10k DP calls per read. Keeping the K closest by edge distance
    # reduces N to ≤25 in the pathological case.
    #
    # Partitioned cap: junctions admitted via N-op match (mapPacBio _leading_del path,
    # line 1281 OR-branch) can sit far from align_5prime — sorting by edge distance
    # would drop them.  Preserve all N-op-matched junctions; apply the distance-based
    # cap only to the remainder.  When _n_intervals is empty (non-mapPacBio, typical
    # for MSP/WM2), the else-branch is the simple O(1) sort-and-slice.
    _MAX_RESCUE_JUNCTIONS = 25
    if len(_nearby_junctions) > _MAX_RESCUE_JUNCTIONS:
        _edge_key = lambda _j: abs(align_5prime - (_j[2] if strand == '+' else _j[1]))
        if _n_intervals:
            _n_matched = [
                _j for _j in _nearby_junctions
                if any(abs(_j[1] - _ns) <= _nb_W2 and abs(_j[2] - _ne) <= _nb_W2
                       for _ns, _ne in _n_intervals)
            ]
            _n_match_set = set(_n_matched)
            _prox_only = sorted(
                (_j for _j in _nearby_junctions if _j not in _n_match_set),
                key=_edge_key,
            )
            _budget = max(0, _MAX_RESCUE_JUNCTIONS - len(_n_matched))
            _nearby_junctions = _n_matched + _prox_only[:_budget]
        else:
            _nearby_junctions.sort(key=_edge_key)
            _nearby_junctions = _nearby_junctions[:_MAX_RESCUE_JUNCTIONS]

    # Pre-compute the "leading-D" pattern flags ONCE per read.  These are read
    # properties (CIGAR-only), independent of the candidate junction, so hoist
    # them out of the per-junction loop.  Without this hoist, the inner loop
    # walks the CIGAR for every (read, junction) pair where dist < 0 — a 4×
    # perf regression on full chunks (Bug 3 patch performance fix, 2026-05-10).
    #
    # _leading_del_plus  applies on plus  strand  (5' end = LEFT of CIGAR).
    # _leading_del_minus applies on minus strand  (5' end = RIGHT of CIGAR).
    # Both default to False — most reads (non-mapPacBio) never trigger them.
    _leading_del_plus = False
    _leading_del_minus = False
    if _n_intervals and read.cigartuples:
        # Collect first ≤3 non-H ops from the LEFT end (plus strand 5').
        _5p_ops_p: List[int] = []
        for _op, _ in read.cigartuples:
            if _op == 5:
                continue   # H: not in query
            _5p_ops_p.append(_op)
            if len(_5p_ops_p) >= 3:
                break
        # Pattern 1 (terminal D): CIGAR starts D …
        # Pattern 2 (D between match and N): CIGAR starts =/M/X D N
        if _5p_ops_p and _5p_ops_p[0] == 2:
            _leading_del_plus = True
        elif (len(_5p_ops_p) >= 3
              and _5p_ops_p[0] in (0, 7, 8)   # M / = / X
              and _5p_ops_p[1] == 2            # D
              and _5p_ops_p[2] == 3):          # N
            _leading_del_plus = True

        # Collect first ≤3 non-H ops from the RIGHT end (minus strand 5').
        _5p_ops_m: List[int] = []
        for _op, _ in reversed(read.cigartuples):
            if _op == 5:
                continue
            _5p_ops_m.append(_op)
            if len(_5p_ops_m) >= 3:
                break
        if _5p_ops_m and _5p_ops_m[0] == 2:
            _leading_del_minus = True
        elif (len(_5p_ops_m) >= 3
              and _5p_ops_m[0] in (0, 7, 8)
              and _5p_ops_m[1] == 2
              and _5p_ops_m[2] == 3):
            _leading_del_minus = True

    if rescue_seq:
        rescue_len = len(rescue_seq)
        _gs = len(genome_seq)
        for j_entry in _nearby_junctions:
            j_chrom, intron_start, intron_end = j_entry[0], j_entry[1], j_entry[2]
            if j_chrom != chrom:
                continue
            if len(j_entry) >= 4 and j_entry[3] not in (strand, '.', ''):
                continue

            if strand == '+':
                # Upstream intron: intron_end must be at or just before align_5prime.
                # When mapPacBio extends into the intron, align_5prime < intron_end
                # (dist < 0) — allow rescue if alignment starts inside the intron
                # (intron_start < align_5prime < intron_end), not completely before it.
                dist = align_5prime - intron_end
                if dist < 0:
                    # If the read already has an N-op for this junction AND has
                    # terminal D ops at the 5' end (mapPacBio artifact: D ops push
                    # align_5prime left of intron_start), suppress the skip so we
                    # can rescue the mangled exon-1 boundary.
                    # We do NOT suppress for correctly-aligned reads (Cat6) whose
                    # exon-1 alignment is clean (no leading D at the 5' end).
                    #
                    # Fast path: _leading_del_plus is precomputed once per read.
                    # For reads without that pattern (the vast majority), we
                    # collapse to the original `if align_5prime <= intron_start: continue`
                    # — restoring May-7 perf for the common case.
                    if _leading_del_plus:
                        _n_match = any(
                            abs(ns - intron_start) <= junction_proximity_bp
                            and abs(ne - intron_end) <= junction_proximity_bp
                            for ns, ne in _n_intervals
                        )
                    else:
                        _n_match = False
                    if not (_n_match and _leading_del_plus) and align_5prime <= intron_start:
                        continue  # alignment is upstream of intron entirely — skip
                    dist = 0  # treat as touching the boundary
                    # Track this junction as a forced-snap candidate (mapPacBio
                    # terminal-D overshoot: N-op proves the read spans the intron,
                    # but garbled post-N query bases will likely fail sequence rescue).
                    if _n_match and _leading_del_plus:
                        _n_err = min(
                            abs(ns - intron_start) + abs(ne - intron_end)
                            for ns, ne in _n_intervals
                            if abs(ns - intron_start) <= junction_proximity_bp
                            and abs(ne - intron_end) <= junction_proximity_bp
                        )
                        if _n_err < _forced_snap_n_err:
                            _forced_snap_n_err = _n_err
                            _forced_snap_junction = (j_chrom, intron_start, intron_end)
                elif dist > junction_proximity_bp + five_clip:
                    # Soft-clip bases extend the read toward the junction; only
                    # filter out junctions that are farther than proximity_bp
                    # even accounting for the unaligned clip length.
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

                # When the alignment's 5' end is inside the intron (dist==0 and
                # align_5prime < intron_end), only query bases mapped to intron
                # positions are candidates for exon-1 sequence.  Bases mapped to
                # exon-2 positions (>= intron_end) are NOT exon-1 sequence and must
                # not participate in candidate scoring — they inflate edit distance
                # against wrong windows and can cause false rescues to shifted donors.
                # Upper bound: intron_end - align_5prime reference positions (1 query
                # base per non-deletion ref base, so this slightly over-counts in the
                # presence of deletions, which is safe).
                if dist == 0 and align_5prime < intron_end:
                    _n_intr = intron_end - align_5prime
                    _rseq = rescue_seq[:_n_intr] if _n_intr < rescue_len else rescue_seq
                    _rlen = len(_rseq)
                else:
                    _rseq = rescue_seq
                    _rlen = rescue_len
                    # Perf cap: donor-adjacent bases (plus strand) sit at the END of
                    # _rseq (rescue_seq[-1] abuts intron_start). Keep the last
                    # _RESCUE_DP_CAP; the candidate window auto-shrinks (built from
                    # _rlen). Skipped for the dist==0 mapPacBio-overshoot branch above
                    # (self-limits via _n_intr) and for terminal-peel overrides, which
                    # pass a deliberately-deeper sequence (mirrors the L~1231 gate).
                    if rescue_seq_override is None and _rlen > _RESCUE_DP_CAP:
                        _rseq = _rseq[-_RESCUE_DP_CAP:]
                        _rlen = _RESCUE_DP_CAP
                if not _rseq:
                    continue

                _best_local_ed: float = _rlen + 1
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
                        _es = _eff_start - _rlen - _off
                        if _es < 0:
                            continue
                        _cand = genome_seq[_es:_eff_start - _off].upper()
                        if len(_cand) < _rlen:
                            continue
                        _ed = _hp_edit_distance(_rseq.upper(), _cand)
                        _shift_abs = abs(_shift)
                        # Two-step scoring (lower tuple = better):
                        #   Step 1 — match quality: minimise HP-edit-distance.
                        #            Only this ranks candidates as primary; all
                        #            other criteria are tiebreakers among ED-ties.
                        #   Step 2 — among ED-ties: prefer placements INSIDE the
                        #            sequence-ambiguity window. The ambig window
                        #            is a match-quality property — a slide inside
                        #            it is genuinely a no-op, while a slide outside
                        #            it changes the read's alignment geometry.
                        #   Step 3 — among ED+in_amb-ties: prefer canonical donor
                        #            (GT/GC). This is signal quality and must rank
                        #            below in_amb so canonical signal is never
                        #            purchased at the cost of a worse-anchored slide.
                        #   Step 4 — among ED+in_amb+donor-ties: smallest |shift|
                        #            from the annotated position.
                        _cur  = (not _in_amb, not _donor_ok, _shift_abs)
                        _best = (not _best_in_amb, not _best_local_canonical,
                                 _best_local_shift_abs)
                        if _ed < _best_local_ed or (
                                _ed == _best_local_ed and _cur < _best):
                            _best_local_ed = _ed
                            _best_local_canonical = _donor_ok
                            _best_in_amb = _in_amb
                            _best_local_shift_abs = _shift_abs
                            exon_seq = _cand
                            _eff_intron_start = _eff_start
                        if _ed == 0:
                            # Perfect match at THIS (shift, offset): later offsets of the
                            # same shift can't beat ED=0 and share its shift-fixed
                            # tiebreakers, so skip them. Must gate on _ed, NOT the
                            # cumulative _best_local_ed (initialized once before the shift
                            # loop): keying on the cumulative value makes every shift AFTER
                            # the first ED=0 shift break after only its first offset,
                            # skipping the offset a later (e.g. canonical-donor) shift needs.
                            # That regressed cat3_minus_2 (see AGENT_FIXES 2026-05-27).
                            break

                if not exon_seq:
                    continue
            else:
                # Minus strand: upstream intron (in transcript) has intron_start ≥ align_5prime.
                # When mapPacBio extends into the intron, align_5prime > intron_start
                # (dist < 0) — allow rescue if alignment ends inside the intron
                # (intron_start < align_5prime < intron_end), not completely past it.
                dist = intron_start - align_5prime
                if dist < 0:
                    # If the read already has an N-op for this junction AND has
                    # terminal D ops at the 5' end (mapPacBio artifact: the D ops
                    # push reference_end — and thus align_5prime — past intron_end),
                    # suppress the skip so we can rescue the mangled exon-1 boundary.
                    # We do NOT suppress for correctly-aligned reads (Cat6) whose
                    # exon-1 alignment is clean (no trailing D at the 5' end).
                    #
                    # Fast path: _leading_del_minus is precomputed once per read
                    # (see hoist near function entry). For reads without that
                    # pattern (the vast majority), this collapses to the original
                    # `if align_5prime >= intron_end: continue`.
                    if _leading_del_minus:
                        _n_match = any(
                            abs(ns - intron_start) <= junction_proximity_bp
                            and abs(ne - intron_end) <= junction_proximity_bp
                            for ns, ne in _n_intervals
                        )
                    else:
                        _n_match = False
                    if not (_n_match and _leading_del_minus) and align_5prime >= intron_end:
                        continue  # alignment is downstream of intron entirely — skip
                    dist = 0  # treat as touching the boundary
                    # Track this junction as a forced-snap candidate (mapPacBio
                    # terminal-D overshoot: N-op proves the read spans the intron,
                    # but garbled post-N query bases will likely fail sequence rescue).
                    if _n_match and _leading_del_minus:
                        _n_err = min(
                            abs(ns - intron_start) + abs(ne - intron_end)
                            for ns, ne in _n_intervals
                            if abs(ns - intron_start) <= junction_proximity_bp
                            and abs(ne - intron_end) <= junction_proximity_bp
                        )
                        if _n_err < _forced_snap_n_err:
                            _forced_snap_n_err = _n_err
                            _forced_snap_junction = (j_chrom, intron_start, intron_end)
                elif dist > junction_proximity_bp + five_clip:
                    # Soft-clip bases extend the read toward the junction; only
                    # filter out junctions that are farther than proximity_bp
                    # even accounting for the unaligned clip length.
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

                # Symmetric guard for minus strand: when align_5prime is inside the
                # intron (dist==0 and align_5prime > intron_start), only query bases
                # mapped to intron positions (< align_5prime in genomic coords, i.e.
                # > intron_start) are candidates for exon-1 sequence.  The right end
                # of the rescue_seq (minus strand 5' = rightmost query bases) may
                # extend into exon-2 territory (positions <= intron_start); those
                # bases must not participate in candidate scoring.
                if dist == 0 and align_5prime > intron_start:
                    _n_intr = align_5prime - intron_start
                    _rseq = rescue_seq[-_n_intr:] if _n_intr < rescue_len else rescue_seq
                    _rlen = len(_rseq)
                else:
                    _rseq = rescue_seq
                    _rlen = rescue_len
                    # Perf cap: donor-adjacent bases (minus strand) sit at the START
                    # of _rseq (rescue_seq[0] abuts intron_end). Keep the first
                    # _RESCUE_DP_CAP; the candidate window auto-shrinks (built from
                    # _rlen). Skipped for the dist==0 mapPacBio-overshoot branch above
                    # and for terminal-peel overrides (deliberately-deeper sequence).
                    if rescue_seq_override is None and _rlen > _RESCUE_DP_CAP:
                        _rseq = _rseq[:_RESCUE_DP_CAP]
                        _rlen = _RESCUE_DP_CAP
                if not _rseq:
                    continue

                _best_local_ed: float = _rlen + 1
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
                        _cand = genome_seq[_cs:_cs + _rlen].upper()
                        if len(_cand) < _rlen:
                            continue
                        _ed = _hp_edit_distance(_rseq.upper(), _cand)
                        _shift_abs = abs(_shift)
                        # Two-step scoring — see + strand block (lines ~1141-1156)
                        # for the full rationale. Tuple ordering matches the +
                        # branch so the two scoring loops are structural mirrors.
                        _cur  = (not _in_amb, not _donor_ok, _shift_abs)
                        _best = (not _best_in_amb, not _best_local_canonical,
                                 _best_local_shift_abs)
                        if _ed < _best_local_ed or (
                                _ed == _best_local_ed and _cur < _best):
                            _best_local_ed = _ed
                            _best_local_canonical = _donor_ok
                            _best_in_amb = _in_amb
                            _best_local_shift_abs = _shift_abs
                            exon_seq = _cand
                            _eff_intron_end = _eff_end
                        if _ed == 0:
                            # Gate on _ed, NOT cumulative _best_local_ed — see the
                            # plus-strand mirror above (regressed cat3_minus_2 when keyed
                            # on the cumulative value; AGENT_FIXES 2026-05-27).
                            break

                if not exon_seq:
                    continue

            if len(exon_seq) < _rlen:
                continue

            ed_exon = _hp_edit_distance(_rseq.upper(), exon_seq)
            # Compare against intronic sequence to avoid rescuing reads that match
            # the intron equally well.  Nanopore homopolymer undercalling means a
            # fixed edit-distance threshold is too strict; instead we rescue when
            # the exon match is ≥30% better than the intron match.
            if strand == '+':
                # Intronic sequence: the last _rlen bases before the 3'SS
                # (i.e. just inside the intron from the splice-acceptor site)
                _ic_start = intron_end - _rlen
                intron_cmp_seq = genome_seq[max(0, _ic_start):intron_end].upper()
            else:
                # Intronic sequence: the first _rlen bases after the 3'SS
                intron_cmp_seq = genome_seq[intron_start:intron_start + _rlen].upper()
            if len(intron_cmp_seq) == _rlen:
                ed_intron = _hp_edit_distance(_rseq.upper(), intron_cmp_seq)
                # Rescue if: perfect exon match, OR exon edit-dist is ≥30% better
                # than intron edit-dist (comparative threshold, not absolute)
                rescue_ok = (ed_exon == 0) or (ed_intron > 0 and ed_exon < ed_intron * 0.70)
            else:
                # Near chromosome boundary — fall back to absolute threshold
                rescue_ok = (ed_exon / _rlen <= max_edit_frac)
            if rescue_ok:
                # Cross-junction tiebreaking tuple (lower = better). Matches the
                # within-junction two-step ordering so the priority is consistent
                # at both scopes: match quality (ED) → match-anchor quality
                # (in_amb) → signal quality (canonical donor, then acceptor).
                #   1. edit distance (hp-aware, float)             — match quality
                #   2. shift within sequence-ambiguity window      — match anchor
                #   3. canonical 5'SS donor (GT/GC plus, AC/GC minus) — donor signal
                #   4. smallest |shift| from annotated position
                #   5. 3'SS acceptor quality: AG=0, CG=1, TG=2, AT=3, other=4
                _cur_outer  = (ed_exon, not _best_in_amb, not _best_local_canonical,
                               _best_local_shift_abs, _acceptor_priority)
                _best_outer = (best_ed, not best_in_amb, not best_is_canonical,
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
            # For soft-clip rescues, rescue_seq has already been truncated to
            # exactly the soft-clip bases (lines above).  _intronic_seq can
            # include one extra body base when the alignment ends exactly at
            # intron_start (the boundary straddle: last M op's final ref base
            # == intron_start), producing a CIGAR with one more query-consuming
            # op than five_prime_soft_clip_length.  The bam_writer guard would
            # then fall back to a flat M block (wrong geometry).  Use rescue_seq
            # directly for soft-clip rescues to keep query span == soft-clip len.
            # For mpb_mismatch rescues _intronic_seq correctly excludes exon-2
            # body bases beyond the intron, so keep the original logic there.
            if rescue_type_candidate == 'softclip':
                _align_seq = rescue_seq
            else:
                _align_seq = _intronic_seq if _intronic_seq else rescue_seq

            # Equivalence-extension: when the alignment overshoots the
            # annotated intron boundary on the body side, the k query bases
            # at the boundary-adjacent end of body M can be re-aligned to
            # the start of the downstream exon (- strand) or the end of the
            # upstream exon (+ strand). The intron length changes accordingly
            # and what would otherwise be a trailing/leading kD op gets
            # absorbed into the rescued M op.
            #
            # The k-sweep below tries k = 1..min(overshoot, _MAX_K) and picks
            # the largest k where the equivalence holds. This enables partial
            # absorption when the full overshoot doesn't qualify.
            #
            # Geometric criteria (overshoot only — undershoot is a separate
            # transformation; deferred):
            #
            #   - strand: ref_end > intron_start by `overshoot`. The k bases
            #     to absorb live at [ref_end - k, ref_end) = [intron_start,
            #     intron_start + k). After absorption they align at
            #     [intron_end, intron_end + k).
            #     Required: genome[intron_start : intron_start + k]
            #            == genome[intron_end : intron_end + k]
            #     AND the read bases at those query positions equal both
            #     (so the original body M was a match, not a mismatch).
            #
            #   + strand: ref_start > intron_end by `overshoot`. The k bases
            #     to absorb live at [ref_start, ref_start + k) = [intron_end
            #     + k, intron_end + 2k). After absorption they move into the
            #     upstream exon at [intron_start - k, intron_start).
            #     Required: genome[ref_start : ref_start + k]
            #            == genome[intron_start - k : intron_start]
            #     AND the read bases at those query positions equal both.
            _upstream_trim = 0
            _MAX_K = 10
            if (
                rescue_type_candidate == 'softclip'
                and strand == '-'
                and read.reference_end is not None
            ):
                _overshoot = read.reference_end - _intron_start
                if 0 < _overshoot:
                    _q = read.query_sequence or ''
                    _scl = _get_5prime_softclip_len(read)
                    _max_k = min(_overshoot, _MAX_K)
                    # Iterate from largest k downward so the first hit is
                    # the maximum-absorption case.
                    for _k_try in range(_max_k, 0, -1):
                        if _scl <= 0 or len(_q) < _scl + _k_try:
                            continue
                        # Borrowed read bases are the last k_try query bases
                        # of body M (just before the trailing soft-clip in BAM).
                        _borrowed = _q[len(_q) - _scl - _k_try : len(_q) - _scl].upper()
                        # Their OLD ref position (where body M aligned them).
                        # For overshoot, this is genome[ref_end - k, ref_end)
                        # = genome[intron_start, intron_start + k) when the
                        # k-window stays within the overshoot region.
                        _ref_old = genome_seq[
                            read.reference_end - _k_try : read.reference_end
                        ].upper()
                        _ref_new = genome_seq[
                            _intron_end : _intron_end + _k_try
                        ].upper()
                        if (
                            len(_ref_old) == _k_try
                            and len(_ref_new) == _k_try
                            and _borrowed == _ref_old        # body M was a real match here
                            and _borrowed == _ref_new        # bases also match at new position
                        ):
                            _upstream_trim = _k_try
                            _align_seq = _borrowed + _align_seq
                            break

            # + strand mirror: body M starts INSIDE the intron region from
            # the high (downstream) edge — i.e., `ref_start < intron_end`.
            # This is the structural mirror of - strand's `ref_end > intron_start`
            # overshoot. Body M's first k query bases sit at
            # `[ref_start, ref_start + k)` = `[intron_end - k, intron_end)`,
            # inside the intron region. Moving them to the upstream-exon-tail
            # `[intron_start - k, intron_start)` produces a canonical-coords
            # BAM (intron_len collapses to `intron_end - intron_start`).
            elif (
                rescue_type_candidate == 'softclip'
                and strand == '+'
                and read.reference_start is not None
            ):
                _undershoot = _intron_end - read.reference_start
                if 0 < _undershoot:
                    _q = read.query_sequence or ''
                    _scl = _get_5prime_softclip_len(read)
                    _max_k = min(_undershoot, _MAX_K)
                    for _k_try in range(_max_k, 0, -1):
                        if _scl <= 0 or len(_q) < _scl + _k_try:
                            continue
                        # Borrowed read bases are the first k_try body-M
                        # query bases (just after the leading 5' soft-clip).
                        _borrowed = _q[_scl : _scl + _k_try].upper()
                        # OLD position: body M's first k bases align at
                        # genome[ref_start : ref_start + k), inside intron.
                        _ref_old = genome_seq[
                            read.reference_start : read.reference_start + _k_try
                        ].upper()
                        # NEW position: upstream-exon-tail
                        # genome[intron_start - k : intron_start).
                        _ref_new = genome_seq[
                            _intron_start - _k_try : _intron_start
                        ].upper()
                        if (
                            len(_ref_old) == _k_try
                            and len(_ref_new) == _k_try
                            and _borrowed == _ref_old        # body M was a real match here
                            and _borrowed == _ref_new        # bases also match at upstream-exon tail
                        ):
                            _upstream_trim = _k_try
                            _align_seq = _align_seq + _borrowed  # + strand: append (clip is right-anchored at intron_start)
                            break

            _exon_cigar_str = ''
            try:
                from ..align.local_aligner import align_clip_to_exon, cigar_ops_to_str
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
                'five_prime_upstream_trim': _upstream_trim,
            }

    # --- Case 4.5: forced N-op snap for mapPacBio terminal-D overshoot ---
    # Fires when the distance gate confirmed _n_match AND _leading_del (the N-op
    # proves the read spans the intron, but terminal D ops push align_5prime past
    # intron_end on minus strand / before intron_start on plus strand), AND the
    # sequence rescue loop above failed to find a matching junction (garbled
    # post-N query bases from mapPacBio's intronic insertion distribution).
    # We trust the N-op evidence and snap directly to the intron boundary.
    # rescue_type='n_op_snap'; edit_distance=-1 signals no sequence comparison.
    if _forced_snap_junction is not None and best_junction is None:
        _fj_chrom, _fj_intron_start, _fj_intron_end = _forced_snap_junction
        _snap_pos = _fj_intron_end if strand == '-' else _fj_intron_start - 1
        return {
            'rescued': True,
            'rescue_type': 'n_op_snap',
            'five_prime_corrected': _snap_pos,
            'rescued_junction': _forced_snap_junction,
            'edit_distance': -1,
            'query_bp': 0,
            'five_prime_exon_cigar': '',
            'five_prime_upstream_trim': 0,
        }

    # --- Case 4: 5' end is strictly inside an annotated intron, no N-op for it ---
    # Fires as a fallback when:
    #   (a) no rescue_seq was extracted (all = ops, completely clean intronic match), OR
    #   (b) rescue_seq was extracted but no candidate junction produced a sequence match.
    # Even a single mismatch (1 bp rescue_seq) that fails the edit-distance check
    # should still be corrected via snap — the position is demonstrably wrong.
    # Only fires if no existing N-op in the CIGAR already covers this intron
    # (prevents double-rescue for reads that have a correct but off-by-a-few-bp N).
    # _n_intervals already computed above (before the sequence-rescue loop).
    for j_entry in _nearby_junctions:
        j_chrom, intron_start, intron_end = j_entry[0], j_entry[1], j_entry[2]
        if j_chrom != chrom:
            continue
        if len(j_entry) >= 4 and j_entry[3] not in (strand, '.', ''):
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
        # Snap to exon-1-side boundary.  five_prime_position is an inclusive
        # aligned-base coordinate, so the plus-strand upstream exon base is
        # intron_start - 1; intron_start itself is the first skipped base.
        snap_pos = intron_end if strand == '-' else intron_start - 1
        intronic_depth = (align_5prime - intron_start) if strand == '-' else (intron_end - align_5prime)
        # Compute exon CIGAR so bam_writer can reroute the intronic tail.
        _clip_bd4 = intron_start if strand == '-' else intron_end
        _intronic_seq4 = _get_intronic_query_bases(read, _clip_bd4, strand)

        # Guard: compare the intronic query bases against both the intron
        # reference (at align_5prime) and the exon-1 reference (at intron_end)
        # using HP-aware edit distance.  For a spliced mRNA mis-aligned into the
        # intron, the exon-1 sequence is the better match → snap is correct.
        # For an unspliced pre-mRNA, the intron sequence is the better match →
        # snap would be wrong, so skip this junction.
        # Also skip when rescue_seq is empty (no imperfect ops = the alignment
        # is clean inside the intron = the bases genuinely belong there).
        if _intronic_seq4:
            _n4 = len(_intronic_seq4)
            if strand == '-':
                # Intron reference: the intronic query bases span [intron_start,
                # align_5prime] in genomic coords; use intron_start as anchor.
                _intron_ref4 = genome_seq[intron_start : intron_start + _n4].upper()
                # Exon-1 reference: just past the 5'SS boundary.
                _exon_ref4 = genome_seq[intron_end : intron_end + _n4].upper()
            else:
                # Plus strand: intronic bases span [align_5prime, intron_end).
                _intron_ref4 = genome_seq[align_5prime : align_5prime + _n4].upper()
                # Exon-1 reference: just before the 5'SS boundary.
                _exon_ref4 = genome_seq[
                    max(0, intron_start - _n4) : intron_start
                ].upper()
            _seq4_upper = _intronic_seq4.upper()
            if len(_intron_ref4) == _n4 and len(_exon_ref4) == _n4:
                _ed_intron4 = _hp_edit_distance(_seq4_upper, _intron_ref4)
                _ed_exon4 = _hp_edit_distance(_seq4_upper, _exon_ref4)
                # Snap only when exon-1 is strictly a better match; ties favour
                # keeping the read in the intron (unspliced interpretation).
                if _ed_intron4 <= _ed_exon4:
                    continue  # Sequence belongs in intron → unspliced, skip snap
            elif not rescue_seq:
                # Cannot compare (boundary conditions) but alignment is clean →
                # conservatively treat as unspliced and skip the snap.
                continue
        elif not rescue_seq:
            # _intronic_seq4 is empty AND no terminal errors: completely clean
            # alignment inside the intron → sequence belongs there → skip snap.
            continue

        _exon_cigar_str4 = ''
        if _intronic_seq4:
            try:
                from ..align.local_aligner import align_clip_to_exon, cigar_ops_to_str
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
            'five_prime_upstream_trim': 0,
        }

    # --- Case 3: proximity-only (no sequence to match, but start is at a 3'SS) ---
    for j_entry in _nearby_junctions:
        j_chrom, intron_start, intron_end = j_entry[0], j_entry[1], j_entry[2]
        if j_chrom != chrom:
            continue
        if len(j_entry) >= 4 and j_entry[3] not in (strand, '.', ''):
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
                'five_prime_upstream_trim': 0,
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
        'five_prime_upstream_trim': 0,
        'repeat_expansion': False,
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
