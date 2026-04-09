#!/usr/bin/env python3
"""
Poly(A) tail trimming for RECTIFY.

This module implements poly(A) tail detection and trimming for technologies that
sequence the poly(A) tail directly (nanopore, Helicos, QuantSeq, etc.).

Problem: Poly(A) tails can be partially sequenced and aligned to the genome,
causing the apparent 3' end to be shifted downstream. We need to identify and
trim these tails to recover the true cleavage and polyadenylation site (CPA).

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, Optional, Tuple
import pysam

from ..config import (
    POLYA_RICHNESS_THRESHOLD,
    POLYA_WINDOW_SIZE,
    MIN_POLYA_LENGTH,
    ADAPTER_POLY_T_MIN,
    ADAPTER_TC_MOTIFS,
)
from ..utils.alignment import extract_soft_clips, get_cigar_stats
from ..utils.genome import reverse_complement


def calculate_a_richness(sequence: str, window: int = POLYA_WINDOW_SIZE) -> float:
    """
    Calculate A-richness in sequence using sliding window.

    Args:
        sequence: DNA sequence (case-insensitive)
        window: Window size for A-richness calculation

    Returns:
        Maximum A-richness across all windows (0.0-1.0)
    """
    if len(sequence) == 0:
        return 0.0

    sequence = sequence.upper()

    # If shorter than window, use full sequence
    if len(sequence) < window:
        return sequence.count('A') / len(sequence)

    # Sliding window
    max_richness = 0.0
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        richness = window_seq.count('A') / window
        max_richness = max(max_richness, richness)

    return max_richness


def detect_adapter_pattern(sequence: str) -> bool:
    """
    Detect RTA adapter patterns in sequence.

    Looks for poly(T) stretches and TC motifs that indicate sequencing adapter.

    Args:
        sequence: DNA sequence (case-insensitive)

    Returns:
        True if adapter pattern detected
    """
    if len(sequence) == 0:
        return False

    sequence = sequence.upper()

    # Check for poly(T) stretch
    if 'T' * ADAPTER_POLY_T_MIN in sequence:
        return True

    # Check for TC motifs
    for motif in ADAPTER_TC_MOTIFS:
        if motif in sequence:
            return True

    return False


def score_polya_tail(
    sequence: str,
    threshold: float = POLYA_RICHNESS_THRESHOLD
) -> Dict:
    """
    Score sequence for poly(A) tail likelihood.

    Args:
        sequence: DNA sequence (typically soft-clip)
        threshold: Minimum A-richness to classify as poly(A)

    Returns:
        Dict with:
            - a_richness: Maximum A-richness across windows
            - is_polya: True if A-richness >= threshold
            - has_adapter: True if adapter pattern detected
            - length: Length of sequence
            - score: Combined score (0.0-1.0)
    """
    a_richness = calculate_a_richness(sequence)
    has_adapter = detect_adapter_pattern(sequence)

    # Classify as poly(A) if A-rich enough
    is_polya = a_richness >= threshold

    # Combined score (higher is more confident)
    score = a_richness

    return {
        'a_richness': a_richness,
        'is_polya': is_polya,
        'has_adapter': has_adapter,
        'length': len(sequence),
        'score': score,
    }


def find_polya_boundary(
    sequence: str,
    strand: str = '+',
    threshold: float = POLYA_RICHNESS_THRESHOLD
) -> int:
    """
    Find boundary between genomic sequence and poly(A) tail.

    Scans from right (+ strand) or left (- strand) to find where A-richness
    drops below threshold.

    Uses O(n) sliding window algorithm for efficiency.

    Args:
        sequence: Aligned read sequence
        strand: Gene strand
        threshold: A-richness threshold

    Returns:
        Position where poly(A) tail begins (0-based index into sequence)
        Returns len(sequence) if no poly(A) found
    """
    if len(sequence) == 0:
        return 0

    seq_len = len(sequence)
    sequence = sequence.upper()
    target_base = 'A' if strand == '+' else 'T'

    if strand == '+':
        # Scan from right (3' end) leftward with sliding window
        # Initialize window at rightmost position
        window_size = min(POLYA_WINDOW_SIZE, seq_len)
        if window_size < 3:
            return seq_len

        # Initial window count
        window_start = seq_len - window_size
        a_count = sum(1 for c in sequence[window_start:seq_len] if c == target_base)

        # Check initial window
        if a_count / window_size < threshold:
            return seq_len

        # Slide window leftward (i is the right edge of window)
        for i in range(seq_len - 1, -1, -1):
            # Determine window boundaries
            new_start = max(0, i - POLYA_WINDOW_SIZE + 1)
            new_end = i + 1
            curr_window_size = new_end - new_start

            if curr_window_size < 3:
                continue

            # Update count incrementally when window moves
            if i < seq_len - 1:
                # Character leaving window (right side)
                if sequence[i + 1] == target_base:
                    a_count -= 1
                # When window shrinks at left boundary, no char enters
                # When window is full size, a char enters at left
                old_start = max(0, i + 1 - POLYA_WINDOW_SIZE + 1)
                if new_start < old_start:
                    if sequence[new_start] == target_base:
                        a_count += 1

            a_richness = a_count / curr_window_size
            if a_richness < threshold:
                return i + 1

        # Entire sequence is A-rich
        return 0
    else:
        # Minus strand: poly-T is at the LEFT (beginning) of the read.
        # Scan left-to-right with an early-exit: advance as long as the
        # current window starting at position i is T-rich, and stop (return i)
        # the moment T-richness drops below threshold.  This prevents the
        # sliding-window from bridging across non-T-rich gaps and returning
        # an incorrectly large boundary.
        window_size = min(POLYA_WINDOW_SIZE, seq_len)
        if window_size < 3:
            return 0

        # Check whether the leftmost window is T-rich at all.
        t_count = sum(1 for c in sequence[0:window_size] if c == target_base)
        if t_count / window_size < threshold:
            return 0

        # Scan rightward: i is the left edge of the current window.
        # We already checked i=0 above; update incrementally from i=1 onward.
        for i in range(1, seq_len):
            window_end = min(seq_len, i + window_size)
            curr_window_size = window_end - i

            if curr_window_size < 3:
                # Window too small to be reliable — treat as end of poly-T.
                return i

            # Character leaving window on the left.
            if sequence[i - 1] == target_base:
                t_count -= 1
            # Character entering window on the right (only when window is full-size).
            old_end = min(seq_len, (i - 1) + window_size)
            if window_end > old_end and window_end <= seq_len:
                if sequence[window_end - 1] == target_base:
                    t_count += 1

            t_richness = t_count / curr_window_size
            if t_richness < threshold:
                # T-richness dropped — poly-T region ends here.
                return i

        # Entire sequence is T-rich.
        return seq_len


def calculate_full_polya_length(
    read: pysam.AlignedSegment,
    strand: str,
    atract_result: Optional[Dict] = None,
    threshold: float = POLYA_RICHNESS_THRESHOLD
) -> Dict:
    """
    Calculate the full observed poly(A) tail length from a read.

    The poly(A) length includes:
    1. Aligned A's: A bases in the genomic A-tract that the poly(A) tail aligned to
    2. Soft-clipped A's: A bases that extend beyond the genomic A-tract

    This represents the total observed poly(A) tail length, from the first non-A
    base (the true CPA site) through all A's until the sequencing adapter.

    NOTE: This function does NOT correct positions. Position correction is handled
    by the A-tract detector. This function only measures poly(A) tail length for
    QC and reporting purposes.

    Args:
        read: pysam AlignedSegment
        strand: Gene strand ('+' or '-')
        atract_result: Optional result from atract_detector.calculate_atract_ambiguity()
                      containing tract_length for aligned A's
        threshold: A-richness threshold for soft-clip poly(A) detection

    Returns:
        Dict with:
            - polya_length: Total poly(A) length (aligned A's + soft-clipped A's)
            - aligned_a_length: A's in aligned region (from A-tract)
            - soft_clip_a_length: A's in soft-clipped region
            - soft_clip_length: Total length of soft-clip at 3' end
            - has_polya: True if poly(A) tail detected (length > 0)
    """
    # Extract soft-clips
    soft_clips = extract_soft_clips(read)

    # Determine which clip is at 3' end
    if strand == '+':
        three_prime_clip = next((c for c in soft_clips if c['side'] == 'right'), None)
    else:
        three_prime_clip = next((c for c in soft_clips if c['side'] == 'left'), None)

    # Initialize result
    result = {
        'polya_length': 0,
        'aligned_a_length': 0,
        'soft_clip_a_length': 0,
        'soft_clip_length': 0,
        'has_polya': False,
    }

    # Get aligned A's from A-tract result (these are genomic A's that poly(A) aligned to)
    if atract_result and atract_result.get('tract_length', 0) > 0:
        result['aligned_a_length'] = atract_result['tract_length']

    # Score soft-clip for poly(A) content
    if three_prime_clip:
        clip_seq = three_prime_clip['seq']
        result['soft_clip_length'] = three_prime_clip['length']

        if clip_seq is None:
            # query_sequence absent from BAM record — skip soft-clip scoring
            return result

        # For - strand, reverse complement to get RNA orientation
        if strand == '-':
            clip_seq = reverse_complement(clip_seq)

        score = score_polya_tail(clip_seq, threshold)

        if score['is_polya'] and three_prime_clip['length'] >= MIN_POLYA_LENGTH:
            # Count actual A's in soft-clip (may be slightly less than length due to errors)
            result['soft_clip_a_length'] = clip_seq.upper().count('A')

    # Total poly(A) length = aligned A's + soft-clipped A's
    result['polya_length'] = result['aligned_a_length'] + result['soft_clip_a_length']
    result['has_polya'] = result['polya_length'] > 0

    return result


def trim_polya_from_read(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Optional[Dict[str, str]] = None,
    threshold: float = POLYA_RICHNESS_THRESHOLD,
    atract_result: Optional[Dict] = None
) -> Dict:
    """
    Detect poly(A) tail and calculate full poly(A) length from read.

    IMPORTANT: This function DETECTS poly(A) tails but does NOT shift positions.
    The soft-clipped bases are already excluded from the alignment by the aligner.
    Genomic A-tracts (which create positional ambiguity) are handled separately
    by the A-tract detector module.

    This function:
    1. Extracts soft-clips at 3' end
    2. Scores for poly(A) likelihood
    3. Calculates full poly(A) tail length (aligned A's + soft-clipped A's)

    Args:
        read: pysam AlignedSegment
        strand: Gene strand ('+' or '-')
        genome: Optional genome dict (unused, kept for API compatibility)
        threshold: A-richness threshold
        atract_result: Optional result from atract_detector with tract_length

    Returns:
        Dict with:
            - original_3prime: 3' end position from alignment (0-based)
            - corrected_3prime: Same as original (no shift from poly-A)
            - shift: Always 0 (soft-clips don't affect genomic position)
            - polya_length: Full poly(A) length (aligned A's + soft-clipped A's)
            - aligned_a_length: A's in aligned region (from A-tract)
            - soft_clip_a_length: A's in soft-clipped region
            - soft_clip_length: Length of soft-clip at 3' end
            - has_polya: True if poly(A) tail detected
    """
    # Get 3' position
    if strand == '+':
        original_3prime = read.reference_end - 1  # 0-based inclusive
    else:
        original_3prime = read.reference_start  # 0-based

    # Calculate full poly(A) length
    polya_result = calculate_full_polya_length(read, strand, atract_result, threshold)

    # Build result dict
    result = {
        'original_3prime': original_3prime,
        'corrected_3prime': original_3prime,  # No shift from poly-A detection
        'shift': 0,  # Soft-clips don't change genomic position
        'polya_length': polya_result['polya_length'],
        'aligned_a_length': polya_result['aligned_a_length'],
        'soft_clip_a_length': polya_result['soft_clip_a_length'],
        'soft_clip_length': polya_result['soft_clip_length'],
        'has_polya': polya_result['has_polya'],
    }

    return result


def calculate_polya_statistics(trimming_results: list) -> Dict:
    """
    Calculate summary statistics for poly(A) trimming.

    Args:
        trimming_results: List of trimming result dicts

    Returns:
        Dict with summary statistics
    """
    if not trimming_results:
        return {
            'total': 0,
            'with_polya': 0,
            'polya_rate': 0.0,
            'mean_polya_length': 0.0,
            'mean_shift': 0.0,
        }

    import numpy as np

    with_polya = [r for r in trimming_results if r['has_polya']]
    polya_lengths = [r['polya_length'] for r in with_polya]
    shifts = [abs(r['shift']) for r in with_polya]

    return {
        'total': len(trimming_results),
        'with_polya': len(with_polya),
        'polya_rate': len(with_polya) / len(trimming_results),
        'mean_polya_length': np.mean(polya_lengths) if polya_lengths else 0.0,
        'median_polya_length': np.median(polya_lengths) if polya_lengths else 0.0,
        'max_polya_length': max(polya_lengths) if polya_lengths else 0,
        'mean_shift': np.mean(shifts) if shifts else 0.0,
        'median_shift': np.median(shifts) if shifts else 0.0,
    }


def format_polya_report(stats: Dict) -> str:
    """
    Format poly(A) trimming statistics as human-readable report.

    Args:
        stats: Statistics dict from calculate_polya_statistics()

    Returns:
        Formatted report string
    """
    report = []
    report.append("=" * 60)
    report.append("Poly(A) Tail Trimming Summary")
    report.append("=" * 60)
    report.append("")

    report.append("Overall:")
    report.append(f"  Total reads:            {stats['total']:,}")
    report.append(f"  With poly(A):           {stats['with_polya']:,} ({stats['polya_rate']:.1%})")
    report.append("")

    if stats['with_polya'] > 0:
        report.append("Poly(A) Tail Lengths:")
        report.append(f"  Mean:                   {stats['mean_polya_length']:.1f} bp")
        report.append(f"  Median:                 {stats['median_polya_length']:.1f} bp")
        report.append(f"  Maximum:                {stats['max_polya_length']} bp")
        report.append("")

        report.append("Position Shifts:")
        report.append(f"  Mean:                   {stats['mean_shift']:.1f} bp")
        report.append(f"  Median:                 {stats['median_shift']:.1f} bp")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)


# =============================================================================
# BAM-level poly(A) tail removal
# =============================================================================

def trim_polya_from_bam_read(
    read: pysam.AlignedSegment,
    strand: str,
    threshold: float = POLYA_RICHNESS_THRESHOLD,
) -> Tuple[pysam.AlignedSegment, int]:
    """
    Remove the 3' poly(A) soft-clip from a BAM read in-place.

    Strips the poly(A)-rich soft-clip at the RNA 3' end from the read's
    CIGAR string, query sequence, and base-quality array.  The alignment
    coordinates are unchanged because the soft-clip is already excluded
    from the aligned region; only the stored query bases are affected.

    The function is conservative: it only removes bases that pass the
    A-richness threshold.  Any non-poly(A) soft-clip bases are preserved.

    Strand convention (consistent with bam_processor.py):
        + strand: 3' soft-clip is the RIGHT clip  (cigartuples[-1] with op=4)
        - strand: 3' soft-clip is the LEFT  clip  (cigartuples[0]  with op=4)

    Args:
        read:      pysam AlignedSegment (mutated in-place)
        strand:    '+' or '-'
        threshold: Minimum fraction of A's (plus) / T's (minus) required
                   for a soft-clip to be considered a poly(A) tail

    Returns:
        (read, n_trimmed) — the modified read and the number of bases removed.
        Returns (read, 0) unchanged if no trimming was applied.

    # TODO — Dorado poly(A) tail length integration:
    # Dorado (≥0.9) estimates the poly(A) tail length for each read and stores
    # it as the BAM tag  pt:i:<length>  in the unaligned BAM.  Once this tag
    # is propagated through alignment and into the consensus BAM, we can use it
    # to replace the A-richness heuristic here:
    #
    #   try:
    #       dorado_polya_len = read.get_tag('pt')
    #   except KeyError:
    #       dorado_polya_len = None
    #
    # When dorado_polya_len is present we should:
    #   1. Accept it as the authoritative tail length instead of measuring the
    #      soft-clip directly (it may be longer than the observable soft-clip
    #      because part of the tail aligned to the genomic A-tract).
    #   2. Pad the soft-clip with artificial A's up to dorado_polya_len before
    #      removing them, so the stored sequence reflects the full tail.
    #   3. Record the Dorado-derived length in the output TSV alongside the
    #      observed soft_clip_a_length for comparison.
    #
    # This feature is gated behind --use-dorado-polya (default: off) and
    # should be developed and tested once Dorado pt-tag support stabilises.
    """
    cigar = read.cigartuples
    seq   = read.query_sequence
    quals = read.query_qualities

    if not cigar or not seq:
        return read, 0

    if strand == '+':
        # 3' end is the right (last) soft-clip
        if cigar[-1][0] != 4:
            return read, 0
        clip_len  = cigar[-1][1]
        if clip_len == 0:
            return read, 0
        clip_seq  = seq[-clip_len:].upper()
        a_frac    = clip_seq.count('A') / clip_len
        if a_frac < threshold:
            return read, 0
        # Trim: remove last clip_len bases
        new_seq   = seq[:-clip_len]
        new_quals = quals[:-clip_len] if quals is not None else None
        new_cigar = cigar[:-1]          # drop the trailing S op
    else:
        # 3' end is the left (first) soft-clip
        if cigar[0][0] != 4:
            return read, 0
        clip_len  = cigar[0][1]
        if clip_len == 0:
            return read, 0
        clip_seq  = seq[:clip_len].upper()
        t_frac    = clip_seq.count('T') / clip_len
        if t_frac < threshold:
            return read, 0
        # Trim: remove first clip_len bases
        new_seq   = seq[clip_len:]
        new_quals = quals[clip_len:] if quals is not None else None
        new_cigar = cigar[1:]           # drop the leading S op

    if not new_cigar:
        # Degenerate: the entire read was a soft-clip — leave untouched
        return read, 0

    read.query_sequence  = new_seq
    read.query_qualities = new_quals
    read.cigartuples     = new_cigar

    return read, clip_len
