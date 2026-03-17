#!/usr/bin/env python3
"""
Statistics and QC utilities for RECTIFY.

This module provides functions for:
- QC metrics calculation
- Confidence scoring
- Summary statistics

Author: Kevin R. Roy
Date: 2026-03-09
"""

from typing import Dict, List, Optional
import numpy as np
from collections import Counter

# =============================================================================
# Confidence Scoring
# =============================================================================

def calculate_confidence(correction_data: Dict) -> str:
    """
    Calculate confidence level for 3' end correction.

    Confidence is based on:
    - NET-seq support (if available)
    - Ambiguity range size
    - Correction type
    - QC flags

    Args:
        correction_data: Dict with correction metadata

    Returns:
        Confidence level: 'high', 'medium', or 'low'
    """
    # NET-seq support
    if correction_data.get('netseq_peak_signal', 0) > 1.0:
        return 'high'
    elif correction_data.get('netseq_peak_signal', 0) > 0.5:
        return 'medium'

    # Ambiguity range
    ambiguity_range = correction_data.get('ambiguity_range', 0)
    if ambiguity_range == 0:
        return 'high'
    elif ambiguity_range <= 2:
        return 'medium'
    else:
        return 'low'


def assign_qc_flag(correction_data: Dict) -> str:
    """
    Assign QC flag based on correction data.

    Flags:
    - PASS: No issues
    - AG_RICH: Likely AG mispriming
    - INDEL_ARTIFACT: Indel artifact detected
    - MULTI_PEAK: Multiple NET-seq peaks
    - AMBIGUOUS: Large ambiguity range without NET-seq support

    Args:
        correction_data: Dict with correction metadata

    Returns:
        QC flag string
    """
    flags = []

    # AG mispriming
    if correction_data.get('is_likely_misprimed', False):
        flags.append('AG_RICH')

    # Indel artifact
    if correction_data.get('indel_artifacts', 0) > 0:
        flags.append('INDEL_ARTIFACT')

    # Multiple NET-seq peaks
    if correction_data.get('n_netseq_peaks', 0) > 1:
        flags.append('MULTI_PEAK')

    # Large ambiguity without NET-seq support
    if (correction_data.get('ambiguity_range', 0) > 3 and
        correction_data.get('netseq_peak_signal', 0) < 0.5):
        flags.append('AMBIGUOUS')

    return ','.join(flags) if flags else 'PASS'


# =============================================================================
# Summary Statistics
# =============================================================================

def calculate_summary_stats(corrections: List[Dict]) -> Dict:
    """
    Calculate summary statistics for batch of corrections.

    Args:
        corrections: List of correction dicts

    Returns:
        Dict with summary statistics
    """
    if not corrections:
        return {
            'total_reads': 0,
            'corrected_reads': 0,
            'mean_shift': 0.0,
            'median_shift': 0.0,
            'mean_ambiguity_range': 0.0,
        }

    shifts = []
    ambiguity_ranges = []
    correction_types = []
    confidences = []
    qc_flags = []

    for corr in corrections:
        # Shift
        raw_pos = corr.get('raw_position', 0)
        corrected_pos = corr.get('corrected_position', raw_pos)
        shift = abs(corrected_pos - raw_pos)
        shifts.append(shift)

        # Ambiguity range
        ambiguity_ranges.append(corr.get('ambiguity_range', 0))

        # Correction type
        correction_types.append(corr.get('correction_type', 'none'))

        # Confidence
        confidences.append(corr.get('confidence', 'low'))

        # QC flags
        qc_flags.append(corr.get('qc_flags', 'PASS'))

    # Calculate stats
    corrected_reads = sum(1 for s in shifts if s > 0)

    stats = {
        'total_reads': len(corrections),
        'corrected_reads': corrected_reads,
        'correction_rate': corrected_reads / len(corrections) if len(corrections) > 0 else 0.0,
        'mean_shift': np.mean(shifts),
        'median_shift': np.median(shifts),
        'mean_ambiguity_range': np.mean(ambiguity_ranges),
        'median_ambiguity_range': np.median(ambiguity_ranges),
        'correction_types': Counter(correction_types),
        'confidences': Counter(confidences),
        'qc_flags': Counter(qc_flags),
    }

    return stats


def format_summary_report(stats: Dict) -> str:
    """
    Format summary statistics as human-readable report.

    Args:
        stats: Summary statistics dict

    Returns:
        Formatted report string
    """
    report = []
    report.append("=" * 60)
    report.append("RECTIFY Summary Report")
    report.append("=" * 60)
    report.append("")

    # Overall stats
    report.append("Overall Statistics:")
    report.append(f"  Total reads:            {stats['total_reads']:,}")
    report.append(f"  Corrected reads:        {stats['corrected_reads']:,}")
    report.append(f"  Correction rate:        {stats['correction_rate']:.1%}")
    report.append("")

    # Position shifts
    report.append("Position Shifts:")
    report.append(f"  Mean shift:             {stats['mean_shift']:.2f} bp")
    report.append(f"  Median shift:           {stats['median_shift']:.2f} bp")
    report.append("")

    # Ambiguity ranges
    report.append("Ambiguity Ranges:")
    report.append(f"  Mean range:             {stats['mean_ambiguity_range']:.2f} bp")
    report.append(f"  Median range:           {stats['median_ambiguity_range']:.2f} bp")
    report.append("")

    # Correction types
    report.append("Correction Types:")
    for corr_type, count in stats['correction_types'].most_common():
        pct = 100.0 * count / stats['total_reads']
        report.append(f"  {corr_type:20s}  {count:7,} ({pct:5.1f}%)")
    report.append("")

    # Confidences
    report.append("Confidence Levels:")
    for conf, count in stats['confidences'].most_common():
        pct = 100.0 * count / stats['total_reads']
        report.append(f"  {conf:20s}  {count:7,} ({pct:5.1f}%)")
    report.append("")

    # QC flags
    report.append("QC Flags:")
    for flag, count in stats['qc_flags'].most_common():
        pct = 100.0 * count / stats['total_reads']
        report.append(f"  {flag:20s}  {count:7,} ({pct:5.1f}%)")
    report.append("")

    report.append("=" * 60)

    return "\n".join(report)


# =============================================================================
# Distribution Analysis
# =============================================================================

def calculate_acount_distribution(corrections: List[Dict]) -> Dict[int, int]:
    """
    Calculate distribution of downstream A-counts.

    Args:
        corrections: List of correction dicts

    Returns:
        Dict mapping A-count to frequency
    """
    acounts = [corr.get('downstream_a_count', 0) for corr in corrections]
    return dict(Counter(acounts))


def calculate_shift_by_acount(corrections: List[Dict]) -> Dict[int, List[float]]:
    """
    Calculate position shifts grouped by A-count.

    Args:
        corrections: List of correction dicts

    Returns:
        Dict mapping A-count to list of shifts
    """
    shifts_by_acount = {}

    for corr in corrections:
        a_count = corr.get('downstream_a_count', 0)
        raw_pos = corr.get('raw_position', 0)
        corrected_pos = corr.get('corrected_position', raw_pos)
        shift = corrected_pos - raw_pos  # Signed shift

        if a_count not in shifts_by_acount:
            shifts_by_acount[a_count] = []

        shifts_by_acount[a_count].append(shift)

    return shifts_by_acount


def calculate_mean_shift_by_acount(corrections: List[Dict]) -> Dict[int, float]:
    """
    Calculate mean position shift for each A-count.

    Args:
        corrections: List of correction dicts

    Returns:
        Dict mapping A-count to mean shift
    """
    shifts_by_acount = calculate_shift_by_acount(corrections)

    mean_shifts = {}
    for a_count, shifts in shifts_by_acount.items():
        mean_shifts[a_count] = np.mean(shifts)

    return mean_shifts
