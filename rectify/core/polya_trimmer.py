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

    sequence = sequence.upper()
    target_base = 'A' if strand == '+' else 'T'

    if strand == '+':
        # Scan from right (3' end) leftward
        for i in range(len(sequence) - 1, -1, -1):
            # Check if we're still in A-rich region
            window_start = max(0, i - POLYA_WINDOW_SIZE + 1)
            window = sequence[window_start:i + 1]

            if len(window) < 3:  # Too short to evaluate
                continue

            a_count = window.count(target_base)
            a_richness = a_count / len(window)

            if a_richness < threshold:
                # Found boundary
                return i + 1

        # Entire sequence is A-rich
        return 0
    else:
        # Scan from left (3' end in - strand) rightward
        for i in range(len(sequence)):
            # Check if we're still in T-rich region
            window_end = min(len(sequence), i + POLYA_WINDOW_SIZE)
            window = sequence[i:window_end]

            if len(window) < 3:
                continue

            t_count = window.count(target_base)
            t_richness = t_count / len(window)

            if t_richness < threshold:
                # Found boundary
                return i

        # Entire sequence is T-rich
        return len(sequence)


def trim_polya_from_read(
    read: pysam.AlignedSegment,
    strand: str,
    genome: Optional[Dict[str, str]] = None,
    threshold: float = POLYA_RICHNESS_THRESHOLD
) -> Dict:
    """
    Detect poly(A) tail in soft-clipped 3' end of read.

    IMPORTANT: This function DETECTS poly(A) tails but does NOT shift positions.
    The soft-clipped bases are already excluded from the alignment by the aligner.
    Genomic A-tracts (which create positional ambiguity) are handled separately
    by the A-tract detector module.

    This function:
    1. Extracts soft-clips at 3' end
    2. Scores for poly(A) likelihood
    3. Reports poly(A) tail length for QC

    Args:
        read: pysam AlignedSegment
        strand: Gene strand ('+' or '-')
        genome: Optional genome dict (unused, kept for API compatibility)
        threshold: A-richness threshold

    Returns:
        Dict with:
            - original_3prime: 3' end position from alignment (0-based)
            - corrected_3prime: Same as original (no shift from poly-A)
            - shift: Always 0 (soft-clips don't affect genomic position)
            - polya_length: Length of poly(A) tail in soft-clip
            - soft_clip_length: Length of soft-clip at 3' end
            - soft_clip_score: Poly(A) score for soft-clip
            - has_polya: True if poly(A) tail detected in soft-clip
    """
    # Extract soft-clips
    soft_clips = extract_soft_clips(read)

    # Determine which clip is at 3' end
    if strand == '+':
        # 3' end is at right
        three_prime_clip = next((c for c in soft_clips if c['side'] == 'right'), None)
        original_3prime = read.reference_end - 1  # 0-based inclusive
    else:
        # 3' end is at left
        three_prime_clip = next((c for c in soft_clips if c['side'] == 'left'), None)
        original_3prime = read.reference_start  # 0-based

    # Initialize result
    # NOTE: corrected_3prime = original_3prime because soft-clips do NOT
    # affect the genomic position. The aligner already excluded them.
    # Genomic A-tract ambiguity is handled by atract_detector.py
    result = {
        'original_3prime': original_3prime,
        'corrected_3prime': original_3prime,  # No shift from poly-A detection
        'shift': 0,  # Soft-clips don't change genomic position
        'polya_length': 0,
        'soft_clip_length': 0,
        'soft_clip_score': None,
        'has_polya': False,
    }

    # Score soft-clip if present
    if three_prime_clip:
        clip_seq = three_prime_clip['seq']
        result['soft_clip_length'] = three_prime_clip['length']

        # For - strand, reverse complement to get RNA orientation
        if strand == '-':
            clip_seq = reverse_complement(clip_seq)

        score = score_polya_tail(clip_seq, threshold)
        result['soft_clip_score'] = score

        if score['is_polya']:
            result['has_polya'] = True
            result['polya_length'] = three_prime_clip['length']

    # NOTE: We intentionally do NOT look for A's in the aligned sequence here.
    # Aligned A's are GENOMIC sequence, not poly(A) tail. The ambiguity they
    # create (poly-A aligning to genomic A-tract) is handled by atract_detector.py
    # which looks at the genomic sequence downstream of the 3' position.

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
