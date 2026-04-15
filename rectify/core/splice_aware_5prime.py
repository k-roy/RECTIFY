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


def _extract_5prime_deletion_bridged_seq(
    read: pysam.AlignedSegment,
    scan_query_bp: int = 50,
    min_deletion_ref_bp: int = 15,
) -> str:
    """Return query bases stranded past a large deletion near the 5' end.

    When the aligner bridges the exon-2/intron boundary using a large D-op
    (instead of an N-op or soft-clip), the bases between the 5' end and the
    deletion map to intronic reference but originate from exon 1.  The
    sliding-window mismatch scan in :func:`_extract_5prime_terminal_error_seq`
    misses these reads because the = (exact-match) ops against intronic sequence
    have zero error rate, even though the alignment is wrong.

    This function recognises the pattern:
        [5' end] …= (=|X|I)* D(≥min_deletion_ref_bp) [rest of read]
    and returns the query bases from the 5' end up to (but not including) the
    large D-op.  Those bases become ``rescue_seq`` for exon-1 local alignment.

    Parameters
    ----------
    scan_query_bp:
        Maximum query bases to scan from the 5' end before giving up.
    min_deletion_ref_bp:
        Minimum D-op length (reference bases deleted) to consider a "bridging"
        deletion.  Small D ops (1–2 bp) are routine nanopore alignment noise.

    Returns
    -------
    str
        Query sequence bases from the 5' end to the bridging deletion, or ``""``
        if no such pattern is found.
    """
    if not read.cigartuples or not read.query_sequence:
        return ''

    ct = read.cigartuples
    query_seq = read.query_sequence
    n_query = len(query_seq)

    _qc = frozenset([0, 1, 4, 7, 8])   # ops that consume query: M, I, S, =, X

    if read.is_reverse:
        # 5' end is at the right (high query index); iterate CIGAR from right.
        query_bases_scanned = 0
        for op, length in reversed(ct):
            if op == 5:          # H: skip (not in query_sequence)
                continue
            if op == 3:          # N: existing splice junction — stop
                break
            if op == 2:          # D: check size
                if length >= min_deletion_ref_bp and query_bases_scanned > 0:
                    return query_seq[n_query - query_bases_scanned:]
                # Small D: not a bridge candidate; keep scanning
            elif op in _qc:
                query_bases_scanned += length
                if query_bases_scanned >= scan_query_bp:
                    break
    else:
        # 5' end is at the left (low query index); iterate CIGAR from left.
        query_bases_scanned = 0
        for op, length in ct:
            if op == 5:
                continue
            if op == 3:
                break
            if op == 2:
                if length >= min_deletion_ref_bp and query_bases_scanned > 0:
                    return query_seq[:query_bases_scanned]
            elif op in _qc:
                query_bases_scanned += length
                if query_bases_scanned >= scan_query_bp:
                    break

    return ''


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

    # Case 2b: large-deletion bridge (only if no explicit soft-clip and no mismatch region)
    # Fires when the aligner used a large D-op to bridge from exon 2 to a region
    # inside the intron, leaving clean = ops near the 5' end.  The mismatch-rate
    # scanner misses this because the = ops have zero error rate against intron
    # reference, but the bases actually originate from exon 1.
    if not rescue_seq:
        deletion_bridged_seq = _extract_5prime_deletion_bridged_seq(read)
        if deletion_bridged_seq:
            rescue_seq = deletion_bridged_seq
            rescue_type_candidate = "deletion_bridge"

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
                    for _off in range(junction_proximity_bp + 1):
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
                    for _off in range(junction_proximity_bp + 1):
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
            _j_chrom, _intron_start, _intron_end = best_junction
            _exon_cigar_str = ''
            try:
                from .local_aligner import align_clip_to_exon, cigar_ops_to_str
                _cigar_ops, _ = align_clip_to_exon(
                    rescue_seq, genome_seq,
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
    # Fires when sequence-based rescue (Cases 1–2) produced no rescue_seq but
    # align_5prime landed inside [intron_start, intron_end).  This happens when an
    # aligner (minimap2, mapPacBio, gapmm2, …) uses D-ops, mismatches, or plain
    # truncation instead of an N-op, leaving the 5' end deep inside the intron.
    # We snap five_prime_corrected to the exon-1-side boundary (intron_end for minus,
    # intron_start for plus) — the nearest annotated splice donor site.
    # Only fires if no existing N-op in the CIGAR already covers this intron
    # (prevents double-rescue for reads that have a correct but off-by-a-few-bp N).
    if not rescue_seq:
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
            return {
                'rescued': True,
                'rescue_type': 'intronic_snap',
                'five_prime_corrected': snap_pos,
                'rescued_junction': (j_chrom, intron_start, intron_end),
                'edit_distance': -1,
                'query_bp': intronic_depth,
                'five_prime_exon_cigar': '',
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
