#!/usr/bin/env python3
"""
AG mispriming detection for RECTIFY.

This module implements AG-richness screening to detect false poly(A) sites
caused by internal priming on A/G-rich genomic regions. This is the original
RECTIFY correction from Roy & Chanfreau 2019.

Problem: Reverse transcriptase with oligo-dT primers can prime internally on
stretches of genomic A's or G's (which pair with T's), creating false pA sites
that can be nucleotides to kilobases away from true termination sites.

Solution: Screen for downstream AG-richness. True pA sites typically have
lower AG content downstream, while misprimed sites show high AG-richness.

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, Optional
from collections import Counter

from ..config import AG_RICHNESS_WINDOW, AG_RICHNESS_THRESHOLD, AG_RICHNESS_MIN_WINDOW
from ..utils.genome import get_downstream_sequence, clamp_position


def calculate_ag_content(sequence: str) -> float:
    """
    Calculate AG-richness of sequence.

    Args:
        sequence: DNA sequence (case-insensitive)

    Returns:
        Fraction of A+G bases (0.0-1.0)
    """
    if len(sequence) == 0:
        return 0.0

    sequence = sequence.upper()
    ag_count = sequence.count('A') + sequence.count('G')

    return ag_count / len(sequence)


def screen_ag_mispriming(
    genome: Dict[str, str],
    chrom: str,
    position: int,
    strand: str,
    window: int = AG_RICHNESS_WINDOW,
    threshold: float = AG_RICHNESS_THRESHOLD,
    min_window: int = AG_RICHNESS_MIN_WINDOW
) -> Dict:
    """
    Screen for AG mispriming at a putative poly(A) site.

    Checks the downstream region for elevated AG-richness, which indicates
    likely internal priming on A/G-rich genomic sequence rather than true
    polyadenylation.

    Args:
        genome: Genome dict from load_genome()
        chrom: Chromosome name (standard format)
        position: Putative pA site position (0-based)
        strand: Gene strand ('+' or '-')
        window: Size of downstream window to check (default: 50bp)
        threshold: AG-richness threshold for flagging (default: 0.65)
        min_window: Minimum window size required (default: 20bp)

    Returns:
        Dict with:
            - ag_content: AG-richness in downstream window (0.0-1.0)
            - is_likely_misprimed: True if AG content >= threshold
            - confidence: 'high' (>=0.75), 'medium' (>=0.65), or 'low' (<0.65)
            - window_size: Actual window size used (may be less than requested near chrom end)
            - base_composition: Dict with counts of A, T, G, C

    Example:
        For a position with 70% AG content downstream:
        - ag_content = 0.70
        - is_likely_misprimed = True (if threshold is 0.65)
        - confidence = 'medium'
    """
    # Get downstream sequence
    downstream_seq = get_downstream_sequence(genome, chrom, position, strand, window)

    # Check if we got sufficient sequence
    actual_window = len(downstream_seq)

    if actual_window < min_window:
        # Insufficient sequence (near chromosome end)
        return {
            'ag_content': None,
            'is_likely_misprimed': False,
            'confidence': 'low',
            'window_size': actual_window,
            'base_composition': {},
            'insufficient_data': True,
        }

    # Calculate AG-richness
    ag_content = calculate_ag_content(downstream_seq)

    # Determine if likely misprimed
    is_misprimed = ag_content >= threshold

    # Assign confidence
    if ag_content >= 0.75:
        confidence = 'high'
    elif ag_content >= threshold:
        confidence = 'medium'
    else:
        confidence = 'low'

    # Get base composition
    seq_upper = downstream_seq.upper()
    base_counts = Counter(seq_upper)

    return {
        'ag_content': ag_content,
        'is_likely_misprimed': is_misprimed,
        'confidence': confidence,
        'window_size': actual_window,
        'base_composition': dict(base_counts),
        'insufficient_data': False,
    }


def screen_ag_mispriming_batch(
    genome: Dict[str, str],
    positions: list,
    window: int = AG_RICHNESS_WINDOW,
    threshold: float = AG_RICHNESS_THRESHOLD
) -> list:
    """
    Screen batch of positions for AG mispriming.

    Args:
        genome: Genome dict
        positions: List of dicts with keys: chrom, position, strand
        window: Downstream window size
        threshold: AG-richness threshold

    Returns:
        List of AG screening results (same order as input)
    """
    results = []

    for pos_data in positions:
        result = screen_ag_mispriming(
            genome,
            pos_data['chrom'],
            pos_data['position'],
            pos_data['strand'],
            window,
            threshold
        )
        results.append(result)

    return results


def calculate_ag_statistics(ag_results: list) -> Dict:
    """
    Calculate summary statistics for AG mispriming screening.

    Args:
        ag_results: List of AG screening dicts

    Returns:
        Dict with summary statistics
    """
    if not ag_results:
        return {
            'total': 0,
            'likely_misprimed': 0,
            'mispriming_rate': 0.0,
            'mean_ag_content': 0.0,
            'by_confidence': {},
        }

    import numpy as np

    # Filter out insufficient data
    valid_results = [r for r in ag_results if not r.get('insufficient_data', False)]

    if not valid_results:
        return {
            'total': len(ag_results),
            'likely_misprimed': 0,
            'mispriming_rate': 0.0,
            'mean_ag_content': None,
            'by_confidence': {},
            'insufficient_data': len(ag_results),
        }

    ag_contents = [r['ag_content'] for r in valid_results]
    misprimed_count = sum(1 for r in valid_results if r['is_likely_misprimed'])
    confidences = [r['confidence'] for r in valid_results]

    return {
        'total': len(ag_results),
        'valid': len(valid_results),
        'likely_misprimed': misprimed_count,
        'mispriming_rate': misprimed_count / len(valid_results),
        'mean_ag_content': np.mean(ag_contents),
        'median_ag_content': np.median(ag_contents),
        'max_ag_content': max(ag_contents),
        'by_confidence': dict(Counter(confidences)),
        'insufficient_data': len(ag_results) - len(valid_results),
    }


def format_ag_report(stats: Dict) -> str:
    """
    Format AG mispriming statistics as human-readable report.

    Args:
        stats: Statistics dict from calculate_ag_statistics()

    Returns:
        Formatted report string
    """
    report = []
    report.append("=" * 60)
    report.append("AG Mispriming Screening Summary")
    report.append("=" * 60)
    report.append("")

    report.append("Overall:")
    report.append(f"  Total positions:        {stats['total']:,}")

    if stats.get('valid'):
        report.append(f"  Valid positions:        {stats['valid']:,}")
        report.append(f"  Likely misprimed:       {stats['likely_misprimed']:,} ({stats['mispriming_rate']:.1%})")
        report.append("")

        report.append("AG Content:")
        if stats.get('mean_ag_content') is not None:
            report.append(f"  Mean:                   {stats['mean_ag_content']:.1%}")
            report.append(f"  Median:                 {stats['median_ag_content']:.1%}")
            report.append(f"  Maximum:                {stats['max_ag_content']:.1%}")
        report.append("")

        report.append("Confidence Levels:")
        for conf in ['high', 'medium', 'low']:
            count = stats['by_confidence'].get(conf, 0)
            pct = 100.0 * count / stats['valid'] if stats['valid'] > 0 else 0.0
            report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")

    if stats.get('insufficient_data', 0) > 0:
        report.append("")
        report.append(f"Insufficient data:        {stats['insufficient_data']:,}")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)


def get_ag_qc_flag(ag_result: Dict) -> str:
    """
    Get QC flag for AG mispriming result.

    Args:
        ag_result: Result dict from screen_ag_mispriming()

    Returns:
        QC flag string
    """
    if ag_result.get('insufficient_data', False):
        return 'INSUFFICIENT_DATA'
    elif ag_result['is_likely_misprimed']:
        if ag_result['confidence'] == 'high':
            return 'AG_RICH_HIGH'
        else:
            return 'AG_RICH_MEDIUM'
    else:
        return 'PASS'
