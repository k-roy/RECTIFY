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

from typing import Dict
from pathlib import Path

from ..config import DOWNSTREAM_WINDOW_SIZE, CHROM_TO_GENOME, CHROM_SIZES
from ..utils.genome import clamp_position


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
    # Get chromosome sequence
    seq = genome.get(chrom)
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
    1. Find the A-tract (or T-tract for - strand) containing/adjacent to position
    2. The true CPA is somewhere upstream of the observed position
    3. Ambiguity window spans from first non-A/T base to the observed position

    Args:
        genome: Genome dict from load_genome()
        chrom: Chromosome name (standard format: chrI, chrII, etc.)
        position: Observed 3' end position from BAM (0-based)
        strand: Gene strand ('+' or '-')
        downstream_bp: Maximum distance to search for tract boundaries

    Returns:
        Dict with:
            - ambiguity_min: Leftmost possible true CPA (0-based)
            - ambiguity_max: Rightmost possible true CPA (0-based)
            - ambiguity_range: Size of uncertainty window (bp)
            - tract_start: Start of A/T-tract
            - tract_end: End of A/T-tract
            - tract_length: Length of A/T-tract
            - has_ambiguity: True if ambiguity_range > 0
            - downstream_a_count: Number of A's in the tract (same as tract_length)
            - expected_shift: Expected shift distance due to poly(A) alignment (same as ambiguity_range)

    Example:
        + strand read at position 1005, genomic sequence GGGAAAAAATTT:
                                                            ^1005
        - A-tract spans 1003-1008 (AAAAAA)
        - First non-A is at 1002 (G)
        - Ambiguity window: [1002, 1005] (true CPA could be 1002, 1003, 1004, or 1005)

        - strand read at position 1000, genomic sequence AAATTTTTTTGGG:
                                                           ^1000
        - T-tract spans 1000-1006 (TTTTTTT)
        - First non-T is at 1007 (G)
        - Ambiguity window: [1000, 1007] (true CPA could be anywhere in T-tract)
    """
    # Find A-tract boundaries
    tract_info = find_atract_boundaries(genome, chrom, position, strand, downstream_bp)

    tract_start = tract_info['tract_start']
    tract_end = tract_info['tract_end']
    tract_length = tract_info['tract_length']
    first_non_a = tract_info['first_non_a']

    # Calculate ambiguity window
    # The true CPA is UPSTREAM of the observed position (poly-A shifts it downstream)
    # The ambiguity window includes all positions where the CPA could be:
    # - The observed position (if no poly-A aligned to genomic A-tract)
    # - Any position upstream within the A-tract (poly-A partially aligned)
    # - The position just before the A-tract (poly-A fully aligned)

    if strand == '+':
        # + strand: upstream = lower coords
        # If position is in an A-tract, CPA could be anywhere from first_non_a to position
        # first_non_a is the position RIGHT BEFORE the A-tract starts
        if tract_length > 0 and tract_start <= position <= tract_end:
            ambiguity_min = first_non_a  # Position just before A-tract (CPA if poly-A fully aligned)
            ambiguity_max = position  # Observed position (CPA if no alignment to A-tract)
        else:
            # Position is not in an A-tract, no ambiguity
            ambiguity_min = position
            ambiguity_max = position

    else:
        # - strand: upstream = higher coords (in gene orientation)
        # If position is in a T-tract, CPA could be anywhere from position to first_non_a
        # first_non_a is the position RIGHT AFTER the T-tract ends
        if tract_length > 0 and tract_start <= position <= tract_end:
            ambiguity_min = position  # Observed position
            ambiguity_max = first_non_a  # Position just after T-tract (CPA if poly-A fully aligned)
        else:
            # Position is not in a T-tract, no ambiguity
            ambiguity_min = position
            ambiguity_max = position

    # Ensure min <= max
    if ambiguity_min > ambiguity_max:
        ambiguity_min, ambiguity_max = ambiguity_max, ambiguity_min

    # Clamp to valid chromosome range
    ambiguity_min = clamp_position(chrom, ambiguity_min)
    ambiguity_max = clamp_position(chrom, ambiguity_max)

    ambiguity_range = ambiguity_max - ambiguity_min

    # Calculate downstream A-count (the number of A's in the tract)
    # This corresponds to how far the poly(A) tail could have shifted the position
    downstream_a_count = tract_length if tract_length > 0 else 0

    # Expected shift is the ambiguity range itself (how far the observed position
    # might be shifted from the true CPA due to poly(A) alignment to genomic A-tract)
    expected_shift = ambiguity_range

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
