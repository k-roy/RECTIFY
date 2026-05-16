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

from ...config import (
    AG_RICHNESS_WINDOW, AG_RICHNESS_SCORE_THRESHOLD,
    AG_RICHNESS_MIN_WINDOW, AG_RICHNESS_THRESHOLD,
)
from ...utils.genome import get_downstream_sequence, clamp_position


def calculate_ag_content(sequence: str) -> float:
    """
    Calculate unweighted AG-richness of sequence.

    Returns:
        Fraction of A+G bases (0.0-1.0)
    """
    if len(sequence) == 0:
        return 0.0
    sequence = sequence.upper()
    return (sequence.count('A') + sequence.count('G')) / len(sequence)


def calculate_weighted_ag_score(sequence: str) -> float:
    """
    Compute the weighted A/G internal-priming score from Roy & Chanfreau 2019.

    The oligo-dT19 primer is 19 nt long.  A's proximal to the priming site
    contribute more to false priming than distal A's; G's also promote
    oligo-dT hybridisation but with lower efficiency.

    Weights (0-based position within the 19 bp window):
        A at positions 0-5  (first 6 bp):  2 points each
        A at positions 6-18 (next 13 bp):  1 point each
        G at any position   (all 19 bp):   0.5 points each

    Maximum possible score = 6×2 + 13×1 + 19×0.5 = 34.5 (pure A+G 19-mer)
    Threshold > 15 → likely internal priming artifact.

    Args:
        sequence: Up to 19 bp downstream sequence (case-insensitive).
                  Shorter sequences are scored over their actual length.

    Returns:
        Weighted score (float ≥ 0.0)
    """
    seq = sequence.upper()
    score = 0.0
    for i, base in enumerate(seq):
        if base == 'A':
            score += 2.0 if i < 6 else 1.0
        elif base == 'G':
            score += 0.5
    return score


def screen_ag_mispriming(
    genome: Dict[str, str],
    chrom: str,
    position: int,
    strand: str,
    window: int = AG_RICHNESS_WINDOW,
    threshold: float = AG_RICHNESS_SCORE_THRESHOLD,
    min_window: int = AG_RICHNESS_MIN_WINDOW
) -> Dict:
    """
    Screen for AG mispriming at a putative poly(A) site.

    Implements the weighted A/G scoring scheme from Roy & Chanfreau 2019
    (PMID 31128237) using a 19 bp downstream window matching the oligo-dT19
    primer.  A score > 15 indicates likely internal priming.

    Args:
        genome: Genome dict from load_genome()
        chrom: Chromosome name (standard format)
        position: Putative pA site position (0-based)
        strand: Gene strand ('+' or '-')
        window: Downstream window size (default: 19 bp, matching oligo-dT19)
        threshold: Weighted score threshold for flagging (default: 15.0)
        min_window: Minimum window size required near chromosome ends (default: 10 bp)

    Returns:
        Dict with:
            - ag_score: Weighted A/G score in downstream window
            - ag_content: Unweighted AG fraction (for backward compatibility)
            - is_likely_misprimed: True if ag_score > threshold
            - window_size: Actual window size used
            - base_composition: Dict with counts of A, T, G, C
    """
    # Get downstream sequence
    downstream_seq = get_downstream_sequence(genome, chrom, position, strand, window)

    # Check if we got sufficient sequence
    actual_window = len(downstream_seq)

    if actual_window < min_window:
        return {
            'ag_score': None,
            'ag_content': None,
            'is_likely_misprimed': False,
            'window_size': actual_window,
            'base_composition': {},
            'insufficient_data': True,
        }

    # Weighted score (Roy & Chanfreau 2019)
    ag_score = calculate_weighted_ag_score(downstream_seq)

    # Unweighted fraction (backward compatibility)
    ag_content = calculate_ag_content(downstream_seq)

    # Determine if likely misprimed
    is_misprimed = ag_score > threshold

    # Get base composition
    seq_upper = downstream_seq.upper()
    base_counts = Counter(seq_upper)

    return {
        'ag_score': ag_score,
        'ag_content': ag_content,
        'is_likely_misprimed': is_misprimed,
        'window_size': actual_window,
        'base_composition': dict(base_counts),
        'insufficient_data': False,
    }


def screen_ag_mispriming_batch(
    genome: Dict[str, str],
    positions: list,
    window: int = AG_RICHNESS_WINDOW,
    threshold: float = AG_RICHNESS_SCORE_THRESHOLD
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
        }

    import numpy as np

    # Filter out insufficient data
    valid_results = [r for r in ag_results if not r.get('insufficient_data', False)]

    if not valid_results:
        return {
            'total': len(ag_results),
            'likely_misprimed': 0,
            'mispriming_rate': 0.0,
            'mean_ag_score': None,
            'mean_ag_content': None,
            'insufficient_data': len(ag_results),
        }

    ag_scores = [r['ag_score'] for r in valid_results if r.get('ag_score') is not None]
    ag_contents = [r['ag_content'] for r in valid_results if r.get('ag_content') is not None]
    misprimed_count = sum(1 for r in valid_results if r['is_likely_misprimed'])

    return {
        'total': len(ag_results),
        'valid': len(valid_results),
        'likely_misprimed': misprimed_count,
        'mispriming_rate': misprimed_count / len(valid_results),
        'mean_ag_score': np.mean(ag_scores) if ag_scores else None,
        'median_ag_score': np.median(ag_scores) if ag_scores else None,
        'max_ag_score': max(ag_scores) if ag_scores else None,
        'mean_ag_content': np.mean(ag_contents) if ag_contents else None,
        'median_ag_content': np.median(ag_contents) if ag_contents else None,
        'max_ag_content': max(ag_contents) if ag_contents else None,
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

        report.append("AG Score (Roy & Chanfreau 2019 weighted, max=34.5):")
        if stats.get('mean_ag_score') is not None:
            report.append(f"  Mean:                   {stats['mean_ag_score']:.2f}")
            report.append(f"  Median:                 {stats['median_ag_score']:.2f}")
            report.append(f"  Maximum:                {stats['max_ag_score']:.2f}")
        if stats.get('mean_ag_content') is not None:
            report.append(f"  Mean AG fraction:       {stats['mean_ag_content']:.1%}")
        report.append("")

    if stats.get('insufficient_data', 0) > 0:
        report.append("")
        report.append(f"Insufficient data:        {stats['insufficient_data']:,}")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)


def get_ag_qc_flag(ag_result: Dict) -> str:
    """
    Get QC flag describing the *genomic* context at the read's called 3' end.

    The flag labels the 19 bp window of genome immediately downstream (in mRNA
    orientation) of the called 3' end, not the read itself. Reads carry an
    AG_RICH flag when that downstream window is A/G-enriched — i.e. it looks
    like the kind of substrate oligo-dT(V) primers will hybridise to, so the
    called 3' end is plausibly an internal-priming artefact rather than a
    bona-fide cleavage/polyadenylation site.

    Single-band classification using the weighted A/G score from Roy & Chanfreau 2019:
      AG_RICH           — score > AG_RICHNESS_SCORE_THRESHOLD (Youden-optimal cutoff, validated by DRS)
      PASS              — score ≤ threshold (likely true poly(A) site)
      INSUFFICIENT_DATA — downstream window shorter than min_window (e.g. near chrom end)

    Args:
        ag_result: Result dict from screen_ag_mispriming()

    Returns:
        QC flag string
    """
    if ag_result.get('insufficient_data', False):
        return 'INSUFFICIENT_DATA'
    if ag_result['is_likely_misprimed']:
        return 'AG_RICH'
    return 'PASS'
