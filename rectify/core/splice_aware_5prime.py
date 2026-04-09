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
                # For plus strand, we want the intron that starts closest to the 5' end
                # and shift to the downstream exon (higher coordinate)
                best_intron = min(overlapping_introns, key=lambda x: x['intron_start'])
                corrected_pos = best_intron['downstream_exon_start']
            else:
                # For minus strand, we want the intron that ends closest to the 5' end
                # and shift to the upstream exon (lower coordinate)
                best_intron = max(overlapping_introns, key=lambda x: x['intron_end'])
                corrected_pos = best_intron['upstream_exon_end'] - 1

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


def _extract_5prime_terminal_error_seq(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    five_clip: int,
    scan_bp: int = 25,
    window: int = 8,
    error_threshold: float = 0.40,
    min_errors: int = 2,
) -> str:
    """Return query bytes from the terminal mismatch/indel region at the 5' end.

    mapPacBio forces mismatches/indels at splice junction boundaries instead of
    producing a soft-clip. This function detects that region using the same
    greedy sliding-window scan used in consensus._get_effective_5prime_clip,
    then returns the actual query bytes so they can be matched against the
    upstream exon sequence.

    Returns empty string when no terminal error region is found or when the
    read is a zero-clip/zero-error truncation (no sequence to rescue with).
    """
    if not read.query_sequence or not read.cigartuples:
        return ""

    chrom = read.reference_name
    if not chrom or chrom not in genome:
        return ""

    ref_seq = genome[chrom]
    query_seq = read.query_sequence
    query_len = len(query_seq)

    try:
        pairs = read.get_aligned_pairs()
    except Exception:
        return ""

    errors: List[int] = []
    if read.is_reverse:
        cutoff_qp = query_len - 1 - five_clip
        for qp, rp in reversed(pairs):
            if qp is None:
                continue
            if qp > cutoff_qp:
                continue
            if len(errors) >= scan_bp:
                break
            if rp is None:
                errors.append(1)
            elif rp < len(ref_seq):
                ref_b = ref_seq[rp].upper()
                read_b = query_seq[qp].upper()
                errors.append(1 if (ref_b != 'N' and read_b != ref_b) else 0)
            else:
                errors.append(0)
    else:
        cutoff_qp = five_clip
        for qp, rp in pairs:
            if qp is None:
                continue
            if qp < cutoff_qp:
                continue
            if len(errors) >= scan_bp:
                break
            if rp is None:
                errors.append(1)
            elif rp < len(ref_seq):
                ref_b = ref_seq[rp].upper()
                read_b = query_seq[qp].upper()
                errors.append(1 if (ref_b != 'N' and read_b != ref_b) else 0)
            else:
                errors.append(0)

    if len(errors) < window:
        return ""

    # Greedy scan: extend terminal boundary while leading windows are high-error
    terminal_end = 0
    for i in range(len(errors) - window + 1):
        w = errors[i:i + window]
        if sum(w) / window >= error_threshold:
            terminal_end = i + window
        else:
            break

    if terminal_end == 0 or sum(errors[:terminal_end]) < min_errors:
        return ""

    # Clip back to the last actual error position so that clean bases at the
    # tail of the sliding window (overshoot) are excluded.  Without this, the
    # extracted sequence would contain correctly-aligned bases that inflate the
    # edit distance and prevent a valid rescue match.
    last_err_idx = max((i for i, e in enumerate(errors[:terminal_end]) if e), default=-1)
    if last_err_idx < 0:
        return ""
    terminal_end = last_err_idx + 1

    # Extract the corresponding query bytes
    if not read.is_reverse:
        region_start = five_clip
        region_end = five_clip + terminal_end
        if region_end <= query_len:
            return query_seq[region_start:region_end]
    else:
        region_end = query_len - five_clip
        region_start = region_end - terminal_end
        if region_start >= 0:
            return query_seq[region_start:region_end]

    return ""


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


def rescue_3ss_truncation(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
    candidate_junctions: Set[Tuple[str, int, int]],
    strand: Optional[str] = None,
    max_edit_frac: float = 0.2,
    junction_proximity_bp: int = 10,
    scan_bp: int = 25,
    window: int = 8,
    error_threshold: float = 0.40,
    min_errors: int = 2,
) -> Dict:
    """Rescue reads truncated or soft-clipped at the exon 2 / 3' splice site boundary.

    Handles three cases in priority order:

    1. **Soft-clip** — explicit S operation at the 5' end (minimap2, gapmm2, uLTRA):
       Clip sequence is compared against the upstream exon end across each candidate
       junction. Identical to the consensus-time _rescue_5prime_softclip logic.

    2. **mapPacBio forced-mismatch / indel** — no explicit soft-clip but terminal
       alignment errors at the 5' end (mapPacBio encodes junction boundary errors
       as mismatches/indels rather than S operations):
       The terminal error bytes are extracted with the same greedy sliding-window
       scan as consensus._get_effective_5prime_clip and compared to the upstream exon.

    3. **Zero-clip / zero-error truncation** — alignment starts at or within
       `junction_proximity_bp` of a known 3'SS with no clipped or mismatched bases:
       No sequence evidence is available; the rescue is recorded as a proximity hit
       (`rescue_type = 'proximity'`) without changing the attributed 5' position.

    Candidate junction pool should contain (chrom, intron_start, intron_end) tuples
    from annotated junctions UNION any novel junctions observed in the first-pass
    alignment (rectified BAM CIGAR junctions + per-aligner BAMs if available).

    Args:
        read: pysam AlignedSegment (from rectified BAM)
        genome: chromosome → sequence dict
        candidate_junctions: Set of (chrom, intron_start, intron_end) to test
        strand: Strand override; inferred from read.is_reverse if None
        max_edit_frac: Edit distance / rescue_seq_len threshold for sequence match
        junction_proximity_bp: Max bp between alignment 5' end and intron edge
            to attempt rescue (both sequence-based and proximity-only)
        scan_bp / window / error_threshold / min_errors: Parameters for MPB
            terminal error detection (same semantics as _get_effective_5prime_clip)

    Returns:
        Dict with keys:
            'rescued'            bool   — True for sequence-confirmed rescue (cases 1 & 2)
            'rescue_type'        str    — 'softclip' | 'mpb_mismatch' | 'proximity' | 'none'
            'five_prime_corrected' int  — Updated 5' genomic position (upstream exon end)
                                          for 'softclip' and 'mpb_mismatch'; raw 5' end
                                          for 'proximity' and 'none'
            'rescued_junction'   tuple  — (chrom, intron_start, intron_end) or None
            'edit_distance'      int    — Edit distance for sequence-based rescues; -1 otherwise
            'query_bp'           int    — Length of rescue sequence used; 0 for proximity
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

    # --- Collect rescue sequence (priority: soft-clip > MPB terminal errors) ---
    rescue_seq = ""
    rescue_type_candidate = "none"

    # Case 1: explicit soft-clip
    if five_clip > 0 and read.query_sequence:
        if not read.is_reverse:
            rescue_seq = read.query_sequence[:five_clip]
        else:
            rescue_seq = read.query_sequence[-five_clip:]
        rescue_type_candidate = "softclip"

    # Case 2: MPB forced-mismatch terminal region (only if no explicit soft-clip)
    if not rescue_seq:
        terminal_seq = _extract_5prime_terminal_error_seq(
            read, genome, five_clip,
            scan_bp=scan_bp, window=window,
            error_threshold=error_threshold, min_errors=min_errors,
        )
        if terminal_seq:
            rescue_seq = terminal_seq
            rescue_type_candidate = "mpb_mismatch"

    # --- Try sequence-based rescue against each candidate junction ---
    best_ed = -1
    best_junction = None
    best_five_prime_corrected = align_5prime

    if rescue_seq:
        rescue_len = len(rescue_seq)
        for j_entry in candidate_junctions:
            j_chrom, intron_start, intron_end = j_entry[0], j_entry[1], j_entry[2]
            if j_chrom != chrom:
                continue

            if strand == '+':
                # Upstream intron: intron_end must be at or just before align_5prime
                dist = align_5prime - intron_end
                if dist < 0 or dist > junction_proximity_bp:
                    continue
                # Exon1 end: last rescue_len bases before the intron donor
                exon_end_start = intron_start - rescue_len
                if exon_end_start < 0:
                    continue
                exon_seq = genome_seq[exon_end_start:intron_start].upper()
            else:
                # Minus strand: upstream intron (in transcript) has intron_start ≥ reference_end
                dist = intron_start - align_5prime
                if dist < 0 or dist > junction_proximity_bp:
                    continue
                # Exon1 end (upstream in transcript = right of intron_end in genome):
                # The BAM stores minus-strand sequence as RC of RNA, so the right soft-clip
                # equals genome_seq[intron_end:intron_end+rescue_len] directly.
                exon_seq = genome_seq[intron_end:intron_end + rescue_len].upper()
                if len(exon_seq) < rescue_len:
                    continue

            if len(exon_seq) < rescue_len:
                continue

            ed_exon = _edit_distance(rescue_seq.upper(), exon_seq)
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
                ed_intron = _edit_distance(rescue_seq.upper(), intron_cmp_seq)
                # Rescue if: perfect exon match, OR exon edit-dist is ≥30% better
                # than intron edit-dist (comparative threshold, not absolute)
                rescue_ok = (ed_exon == 0) or (ed_intron > 0 and ed_exon < ed_intron * 0.70)
            else:
                # Near chromosome boundary — fall back to absolute threshold
                rescue_ok = (ed_exon / rescue_len <= max_edit_frac)
            if rescue_ok:
                if best_ed == -1 or ed_exon < best_ed:
                    best_ed = ed_exon
                    best_junction = (j_chrom, intron_start, intron_end)
                    # Update 5' end: upstream boundary of the intron (= end of exon1)
                    if strand == '+':
                        best_five_prime_corrected = intron_start - 1
                    else:
                        best_five_prime_corrected = intron_end

        if best_junction is not None:
            return {
                'rescued': True,
                'rescue_type': rescue_type_candidate,
                'five_prime_corrected': best_five_prime_corrected,
                'rescued_junction': best_junction,
                'edit_distance': best_ed,
                'query_bp': len(rescue_seq),
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
    if result1.starts_in_intron and result1.five_prime_corrected == 1000:
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
    if result3.starts_in_intron and result3.five_prime_corrected == 2099:
        print("  PASSED")
    else:
        print("  FAILED")

    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    validate_5prime_correction()
