#!/usr/bin/env python3
"""
TSV output and summary reporting for RECTIFY correction results.

Author: Kevin R. Roy
"""

from typing import Dict, List

from .processing_stats import ProcessingStats, generate_stats_report


CORRECTION_TSV_HEADER = [
    'read_id', 'chrom', 'strand',
    'original_3prime', 'corrected_3prime',
    'five_prime_position',
    'five_prime_rescued',
    'five_prime_exon_cigar',
    'alignment_start', 'alignment_end',
    'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
    'polya_length',
    'aligned_a_length', 'soft_clip_a_length',
    'junctions', 'n_junctions',
    'five_prime_soft_clip_length', 'three_prime_soft_clip_length',
    'mapq',
    'correction_applied', 'confidence', 'qc_flags', 'fraction',
    'gene_id',
    'pt_tag',
    'polya_score',
    'polya_source',
    'sc_homopolymer_extension',
    'sc_rescued_seq',
    'sc_original_softclip_len',
    'five_prime_intron_clip_pos',
    'oc_homopolymer_extension',
    'oc_overcall_count',
    'oc_terminal_base',
    'five_prime_upstream_trim',
    'reanchor_clip_len',
    # How `strand` was determined. Populated by ONT PCR-cDNA only (other
    # protocols carry a fixed protocol-level strand rule); one of
    # polyA_3p / polyT_5p / gene_overlap / unassigned. Appended last so
    # existing positional consumers are unaffected.
    'strand_evidence',
]


def correction_result_to_tsv_row(result: Dict) -> List[str]:
    """Serialize one correction result using ``CORRECTION_TSV_HEADER`` order."""
    _pt = result.get('pt_tag')
    _ps = result.get('polya_score')
    return [
        result['read_id'],
        result['chrom'],
        result['strand'],
        str(result['original_3prime']),
        str(result['corrected_3prime']),
        str(result.get('five_prime_position', '')),
        '1' if result.get('five_prime_rescued') else '0',
        result.get('five_prime_exon_cigar') or '',
        str(result.get('alignment_start', '')),
        str(result.get('alignment_end', '')),
        str(result['ambiguity_min']),
        str(result['ambiguity_max']),
        str(result['ambiguity_range']),
        str(result.get('polya_length', 0)),
        str(result.get('aligned_a_length', 0)),
        str(result.get('soft_clip_a_length', 0)),
        result.get('junctions_str', ''),
        str(result.get('n_junctions', 0)),
        str(result.get('five_prime_soft_clip_length', 0)),
        str(result.get('three_prime_soft_clip_length', 0)),
        str(result.get('mapq', 0)),
        ','.join(result['correction_applied']) if result['correction_applied'] else 'none',
        result['confidence'],
        ','.join(result['qc_flags']),
        f"{result.get('fraction', 1.0):.6f}",
        result.get('gene_id') or '',
        str(_pt) if _pt is not None else '',
        f'{_ps:.4f}' if _ps is not None else '',
        result.get('polya_source', 'none'),
        str(result.get('sc_homopolymer_extension', 0)),
        result.get('sc_rescued_seq', ''),
        str(result.get('sc_original_softclip_len', 0)),
        str(result.get('five_prime_intron_clip_pos', -1)),
        str(result.get('oc_homopolymer_extension', 0)),
        str(result.get('oc_overcall_count', 0)),
        result.get('oc_terminal_base', ''),
        str(result.get('five_prime_upstream_trim', 0)),
        str(result.get('reanchor_clip_len', 0)),
        result.get('strand_evidence', '') or '',
    ]


def write_output_tsv(results: List[Dict], output_path: str):
    """
    Write correction results to TSV file.

    Args:
        results: List of correction result dicts
        output_path: Path to output TSV file
    """
    print(f"Writing output to {output_path}...")

    with open(output_path, 'w') as f:
        # Write header
        f.write('\t'.join(CORRECTION_TSV_HEADER) + '\n')

        # Write results
        for result in results:
            f.write('\t'.join(correction_result_to_tsv_row(result)) + '\n')

    print(f"  Wrote {len(results):,} corrected positions")


def generate_summary_report(results: List[Dict]) -> str:
    """
    Generate summary report from correction results.

    Uses single-pass accumulation for efficiency with large datasets.

    Args:
        results: List of correction result dicts

    Returns:
        Formatted report string
    """
    import numpy as np

    n_total = len(results)
    if n_total == 0:
        return "No reads processed."

    # Single-pass accumulation of all metrics
    n_with_ambiguity = 0
    n_corrected = 0
    ambiguity_ranges = []
    corrected_shifts = []
    by_confidence = {'high': 0, 'medium': 0, 'low': 0}

    for r in results:
        # Ambiguity check
        ambig_range = r['ambiguity_range']
        ambiguity_ranges.append(ambig_range)
        if ambig_range > 0:
            n_with_ambiguity += 1

        # Correction check
        shift = abs(r['corrected_3prime'] - r['original_3prime'])
        if shift > 0:
            n_corrected += 1
            corrected_shifts.append(shift)

        # Confidence
        conf = r['confidence']
        if conf in by_confidence:
            by_confidence[conf] += 1

    # Build report
    report = []
    report.append("=" * 60)
    report.append("RECTIFY Correction Summary")
    report.append("=" * 60)
    report.append("")
    report.append("Overall:")
    report.append(f"  Total reads:            {n_total:,}")
    report.append(f"  With ambiguity:         {n_with_ambiguity:,} ({100*n_with_ambiguity/n_total:.1f}%)")
    report.append(f"  Position corrected:     {n_corrected:,} ({100*n_corrected/n_total:.1f}%)")
    report.append("")

    report.append("Ambiguity Ranges:")
    report.append(f"  Mean:                   {np.mean(ambiguity_ranges):.2f} bp")
    report.append(f"  Median:                 {np.median(ambiguity_ranges):.1f} bp")
    report.append(f"  Maximum:                {max(ambiguity_ranges)} bp")
    report.append("")

    if n_corrected > 0:
        report.append("Position Corrections:")
        report.append(f"  Mean shift:             {np.mean(corrected_shifts):.2f} bp")
        report.append(f"  Median shift:           {np.median(corrected_shifts):.1f} bp")
        report.append(f"  Maximum shift:          {max(corrected_shifts)} bp")
        report.append("")

    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'low']:
        count = by_confidence.get(conf, 0)
        pct = 100.0 * count / n_total
        report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")

    report.append("")
    report.append("=" * 60)

    return "\n".join(report)


def generate_summary_from_stats(stats) -> str:
    """
    Generate summary report from ProcessingStats or legacy dict.

    Args:
        stats: ProcessingStats object or legacy Dict from process_bam_streaming

    Returns:
        Formatted report string
    """
    # Handle both ProcessingStats and legacy dict format
    if isinstance(stats, ProcessingStats):
        return generate_stats_report(stats)

    # Legacy dict format (backwards compatibility)
    n_total = stats.get('total_reads', 0)
    if n_total == 0:
        return "No reads processed."

    report = []
    report.append("=" * 60)
    report.append("RECTIFY Correction Summary (Streaming Mode)")
    report.append("=" * 60)
    report.append("")
    report.append("Overall:")
    report.append(f"  Total reads:            {n_total:,}")
    report.append(f"  With ambiguity:         {stats['with_ambiguity']:,} ({100*stats['with_ambiguity']/n_total:.1f}%)")
    report.append(f"  Position corrected:     {stats['position_corrected']:,} ({100*stats['position_corrected']/n_total:.1f}%)")
    report.append("")
    report.append("Confidence Levels:")
    for conf in ['high', 'medium', 'low']:
        count = stats['by_confidence'].get(conf, 0)
        pct = 100.0 * count / n_total
        report.append(f"  {conf:10s}          {count:7,} ({pct:5.1f}%)")
    report.append("")
    report.append("=" * 60)

    return "\n".join(report)
