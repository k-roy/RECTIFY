#!/usr/bin/env python3
"""
A-tract ambiguity detection for RECTIFY.

This module implements universal A-tract ambiguity detection that applies to ALL
poly(A)-tailed RNA-seq technologies. Genomic A-tracts near true 3' ends create
positional uncertainty due to poly(A) tail alignment.

Key insight: Poly(A) tails can align to genomic A-tracts, shifting the apparent
3' end downstream. The ambiguity window spans from the first non-A base upstream
to the observed position.

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, Optional
from pathlib import Path

from ..config import DOWNSTREAM_WINDOW_SIZE, CHROM_TO_GENOME, CHROM_SIZES, get_shift_from_acount
from ..utils.genome import clamp_position


def _get_sequence(genome: Dict[str, str], chrom: str) -> Optional[str]:
    """
    Get chromosome sequence, handling standard to NCBI name conversion.

    Args:
        genome: Genome dict with either standard (chrI) or NCBI (ref|NC_...|) keys
        chrom: Chromosome name (can be either format)

    Returns:
        Sequence string or None if not found
    """
    # Try direct lookup first
    if chrom in genome:
        return genome[chrom]

    # Try converting standard name to NCBI format
    if chrom in CHROM_TO_GENOME:
        ncbi_name = CHROM_TO_GENOME[chrom]
        if ncbi_name in genome:
            return genome[ncbi_name]

    return None


def _count_downstream_as(
    genome: Dict[str, str],
    chrom: str,
    position: int,
    strand: str,
    window_size: int = DOWNSTREAM_WINDOW_SIZE
) -> Optional[int]:
    """
    Count A's (or T's for - strand) in downstream window.

    For + strand: Count A's going rightward from position
    For - strand: Count T's going leftward from position (represents A's in RNA orientation)

    Args:
        genome: Genome dict
        chrom: Chromosome name
        position: Position to start from (0-based)
        strand: '+' or '-'
        window_size: Size of window to count in

    Returns:
        Count of A's (or T's) in window, or None if position is invalid
    """
    seq = _get_sequence(genome, chrom)
    if seq is None:
        return None

    chrom_len = len(seq)
    if position < 0 or position >= chrom_len:
        return None

    target_base = 'A' if strand == '+' else 'T'

    if strand == '+':
        # Downstream = rightward for + strand
        window_end = min(position + window_size, chrom_len)
        if window_end <= position:
            return None
        window_seq = seq[position:window_end]
    else:
        # Downstream = leftward for - strand
        window_start = max(position - window_size + 1, 0)
        if window_start >= position + 1:
            return None
        window_seq = seq[window_start:position + 1]

    return window_seq.count(target_base)


def find_atract_boundaries(
    genome: Dict[str, str],
    chrom: str,
    position: int,
    strand: str,
    max_search: int = 20
) -> Dict:
    """
    Find the boundaries of the A-tract containing/adjacent to the position.

    For + strand: Look for A's starting at position and extending RIGHT
    For - strand: Look for T's starting at position and extending LEFT

    Args:
        genome: Genome dict {chrom: sequence}
        chrom: Chromosome name (standard format: chrI, chrII, etc.)
        position: Position to start search (0-based)
        strand: Gene strand ('+' or '-')
        max_search: Maximum bp to search for tract boundaries

    Returns:
        Dict with:
            - tract_start: First position of A/T-tract (0-based)
            - tract_end: Last position of A/T-tract (0-based, inclusive)
            - tract_length: Length of tract (0 if position is not in a tract)
            - first_non_a: Position of first non-A base upstream of tract
    """
    # Get chromosome sequence (with name conversion)
    seq = _get_sequence(genome, chrom)
    if seq is None:
        return {
            'tract_start': position,
            'tract_end': position,
            'tract_length': 0,
            'first_non_a': position,
        }

    chrom_len = len(seq)

    # Determine target base
    target_base = 'A' if strand == '+' else 'T'

    # IMPORTANT: First check if the position itself is the target base
    # If not, there's no tract at this position
    if position < 0 or position >= chrom_len or seq[position] != target_base:
        return {
            'tract_start': position,
            'tract_end': position,
            'tract_length': 0,
            'first_non_a': position,
        }

    if strand == '+':
        # For + strand: A-tract extends rightward (downstream)
        # Poly-A tail would align to genomic A's at position and to the right

        # Find rightmost extent of A-tract (downstream boundary)
        tract_end = position
        while tract_end + 1 < chrom_len and tract_end - position < max_search:
            if seq[tract_end + 1] == target_base:
                tract_end += 1
            else:
                break

        # Find leftmost extent of A-tract (upstream boundary)
        tract_start = position
        while tract_start > 0 and position - tract_start < max_search:
            if seq[tract_start - 1] == target_base:
                tract_start -= 1
            else:
                break

        # First non-A is one position before tract_start
        first_non_a = tract_start - 1 if tract_start > 0 else 0

    else:
        # For - strand: T-tract extends leftward (downstream in gene orientation)
        # Poly-A tail (as T in genomic coords) would align to genomic T's

        # Find leftmost extent of T-tract (downstream boundary in gene coords)
        tract_start = position
        while tract_start > 0 and position - tract_start < max_search:
            if seq[tract_start - 1] == target_base:
                tract_start -= 1
            else:
                break

        # Find rightmost extent of T-tract (upstream boundary in gene coords)
        tract_end = position
        while tract_end + 1 < chrom_len and tract_end - position < max_search:
            if seq[tract_end + 1] == target_base:
                tract_end += 1
            else:
                break

        # First non-T is one position after tract_end (upstream in gene coords)
        first_non_a = tract_end + 1 if tract_end + 1 < chrom_len else tract_end

    tract_length = tract_end - tract_start + 1

    return {
        'tract_start': tract_start,
        'tract_end': tract_end,
        'tract_length': tract_length,
        'first_non_a': first_non_a,
    }


def calculate_atract_ambiguity(
    genome: Dict[str, str],
    chrom: str,
    position: int,
    strand: str,
    downstream_bp: int = DOWNSTREAM_WINDOW_SIZE
) -> Dict:
    """
    Calculate A-tract ambiguity window for a 3' end position.

    This is a UNIVERSAL correction that applies to all poly(A)-tailed RNA-seq.
    The ambiguity arises because poly(A) tails can align to genomic A-tracts,
    making it impossible to determine the exact CPA site without external validation.

    Algorithm:
    1. Count A's (or T's for - strand) in downstream window
    2. Look up expected shift from calibration table
    3. Calculate ambiguity window based on expected shift

    Args:
        genome: Genome dict from load_genome()
        chrom: Chromosome name (standard format: chrI, chrII, etc.)
        position: Observed 3' end position from BAM (0-based)
        strand: Gene strand ('+' or '-')
        downstream_bp: Window size for counting downstream A's

    Returns:
        Dict with:
            - ambiguity_min: Leftmost possible true CPA (0-based)
            - ambiguity_max: Rightmost possible true CPA (0-based)
            - ambiguity_range: Size of uncertainty window (bp)
            - tract_start: Start of contiguous A/T-tract
            - tract_end: End of contiguous A/T-tract
            - tract_length: Length of contiguous A/T-tract
            - has_ambiguity: True if ambiguity_range > 0
            - downstream_a_count: Number of A's in downstream window
            - expected_shift: Expected shift from calibration table

    Example:
        + strand read at position 1010, downstream 10bp = 'TCGATCGATC' has 3 A's:
        - downstream_a_count = 3
        - expected_shift = 0.4 (from config table)
        - ambiguity_min = int(1010 - 0.4) = 1009
        - ambiguity_max = 1010
        - ambiguity_range = 1
    """
    # Count A's in downstream window
    downstream_a_count = _count_downstream_as(genome, chrom, position, strand, downstream_bp)

    # Handle invalid positions
    if downstream_a_count is None:
        return {
            'ambiguity_min': position,
            'ambiguity_max': position,
            'ambiguity_range': 0,
            'tract_start': position,
            'tract_end': position,
            'tract_length': 0,
            'has_ambiguity': False,
            'downstream_a_count': None,
            'expected_shift': 0.0,
        }

    # Get expected shift from calibration table
    expected_shift = get_shift_from_acount(downstream_a_count)

    # Find contiguous A-tract boundaries for tract_length calculation
    tract_info = find_atract_boundaries(genome, chrom, position, strand, downstream_bp)
    tract_start = tract_info['tract_start']
    tract_end = tract_info['tract_end']
    tract_length = tract_info['tract_length']

    # Calculate ambiguity window based on expected shift
    if strand == '+':
        # + strand: poly-A shifts position RIGHTWARD (downstream)
        # True CPA is LEFTWARD (upstream) by expected_shift
        ambiguity_min = int(position - expected_shift)
        ambiguity_max = position
    else:
        # - strand: poly-A shifts position LEFTWARD (downstream in gene coords)
        # True CPA is RIGHTWARD (upstream in gene coords) by expected_shift
        ambiguity_min = position
        ambiguity_max = int(position + expected_shift) + 1  # +1 because shift can be fractional

    # Ensure min <= max
    if ambiguity_min > ambiguity_max:
        ambiguity_min, ambiguity_max = ambiguity_max, ambiguity_min

    # Clamp to valid chromosome range
    ambiguity_min = clamp_position(chrom, ambiguity_min)
    ambiguity_max = clamp_position(chrom, ambiguity_max)

    ambiguity_range = ambiguity_max - ambiguity_min

    return {
        'ambiguity_min': ambiguity_min,
        'ambiguity_max': ambiguity_max,
        'ambiguity_range': ambiguity_range,
        'tract_start': tract_start,
        'tract_end': tract_end,
        'tract_length': tract_length,
        'has_ambiguity': ambiguity_range > 0,
        'downstream_a_count': downstream_a_count,
        'expected_shift': expected_shift,
    }


def calculate_atract_ambiguity_batch(
    genome: Dict[str, str],
    positions: list,
    downstream_bp: int = DOWNSTREAM_WINDOW_SIZE
) -> list:
    """
    Calculate A-tract ambiguity for batch of positions.

    Args:
        genome: Genome dict
        positions: List of dicts with keys: chrom, position, strand
        downstream_bp: Window size for downstream A-count

    Returns:
        List of ambiguity dicts (same order as input)
    """
    results = []

    for pos_data in positions:
        ambiguity = calculate_atract_ambiguity(
            genome,
            pos_data['chrom'],
            pos_data['position'],
            pos_data['strand'],
            downstream_bp
        )
        results.append(ambiguity)

    return results


def get_ambiguity_category(ambiguity_range: int) -> str:
    """
    Categorize ambiguity range for reporting.

    Args:
        ambiguity_range: Size of ambiguity window (bp)

    Returns:
        Category: 'none', 'low', 'medium', or 'high'
    """
    if ambiguity_range == 0:
        return 'none'
    elif ambiguity_range <= 2:
        return 'low'
    elif ambiguity_range <= 4:
        return 'medium'
    else:
        return 'high'


def summarize_ambiguity_distribution(ambiguities: list) -> Dict:
    """
    Summarize distribution of ambiguity ranges.

    Args:
        ambiguities: List of ambiguity dicts from calculate_atract_ambiguity()

    Returns:
        Dict with summary statistics
    """
    if not ambiguities:
        return {
            'total': 0,
            'mean_range': 0.0,
            'median_range': 0.0,
            'by_category': {},
            'by_acount': {},
        }

    import numpy as np
    from collections import Counter

    ranges = [a['ambiguity_range'] for a in ambiguities]
    acounts = [a['downstream_a_count'] for a in ambiguities if a['downstream_a_count'] is not None]
    categories = [get_ambiguity_category(r) for r in ranges]

    return {
        'total': len(ambiguities),
        'mean_range': np.mean(ranges),
        'median_range': np.median(ranges),
        'max_range': max(ranges),
        'with_ambiguity': sum(1 for r in ranges if r > 0),
        'by_category': dict(Counter(categories)),
        'by_acount': dict(Counter(acounts)),
    }


def format_ambiguity_report(summary: Dict) -> str:
    """
    Format ambiguity summary as human-readable report.

    Args:
        summary: Summary dict from summarize_ambiguity_distribution()

    Returns:
        Formatted report string
    """
    report = []
    report.append("=" * 60)
    report.append("A-tract Ambiguity Summary")
    report.append("=" * 60)
    report.append("")

    report.append("Overall:")
    report.append(f"  Total positions:        {summary['total']:,}")
    report.append(f"  With ambiguity:         {summary['with_ambiguity']:,} ({100.0 * summary['with_ambiguity'] / summary['total']:.1f}%)")
    report.append(f"  Mean range:             {summary['mean_range']:.2f} bp")
    report.append(f"  Median range:           {summary['median_range']:.1f} bp")
    report.append(f"  Max range:              {summary['max_range']} bp")
    report.append("")

    report.append("By Category:")
    for cat in ['none', 'low', 'medium', 'high']:
        count = summary['by_category'].get(cat, 0)
        pct = 100.0 * count / summary['total'] if summary['total'] > 0 else 0.0
        report.append(f"  {cat:10s}          {count:7,} ({pct:5.1f}%)")
    report.append("")

    report.append("By Downstream A-count:")
    for a_count in sorted(summary['by_acount'].keys()):
        count = summary['by_acount'][a_count]
        pct = 100.0 * count / summary['total'] if summary['total'] > 0 else 0.0
        report.append(f"  {a_count:2d}A                {count:7,} ({pct:5.1f}%)")
    report.append("")

    report.append("=" * 60)

    return "\n".join(report)
