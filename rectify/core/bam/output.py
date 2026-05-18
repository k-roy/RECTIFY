#!/usr/bin/env python3
"""
TSV output and summary reporting for RECTIFY correction results.

Author: Kevin R. Roy
"""

from typing import Dict, List

from .processing_stats import ProcessingStats, generate_stats_report


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
        header = [
            'read_id', 'chrom', 'strand',
            'original_3prime', 'corrected_3prime',
            'five_prime_position',  # TSS end of the read (v2.6.0)
            'five_prime_rescued',   # 1 if 5' end was corrected by junction rescue (v2.7.9)
            'five_prime_exon_cigar',  # SAM CIGAR for exon segment of Cat3 rescue (v2.8.0)
            'alignment_start', 'alignment_end',  # Full read body interval (v2.6.0)
            'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
            'polya_length',  # Total observed poly(A) tail length
            'aligned_a_length', 'soft_clip_a_length',  # Breakdown of poly(A)
            'junctions', 'n_junctions',  # Splice junctions (v2.7.0)
            'five_prime_soft_clip_length', 'three_prime_soft_clip_length',  # Soft clips (v2.7.0)
            'mapq',  # Mapping quality (v2.7.0)
            'correction_applied', 'confidence', 'qc_flags', 'fraction',
            'gene_id',  # Per-read gene attribution (optional)
            'pt_tag',      # dorado pt:i signal-level poly(A) length (blank if absent) (v2.9.0)
            'polya_score', # poly(A) model confidence 0-1 (blank if no model) (v2.9.0)
            'polya_source',  # 'pt_tag' | 'model' | 'none' (v2.9.0)
            'sc_homopolymer_extension',  # Cat2: under-called homopolymer bases (v2.9.1)
            'sc_rescued_seq',            # Cat2: non-poly-A bases matched to ref (v2.9.1)
            'sc_original_softclip_len',  # Cat2: original 3' soft-clip length (v2.9.1)
            'five_prime_intron_clip_pos',  # Case 4: exon-side intron boundary for BAM hard-clip (-1 if N/A)
            'oc_homopolymer_extension',  # over-call rescue: genomic HP extension past raw_pos (D op)
            'oc_overcall_count',         # over-call rescue: # of over-call stop-base bases (I op)
            'oc_terminal_base',          # over-call rescue: the terminal non-stop-base char (= op)
            'five_prime_upstream_trim',  # Cat3 equivalence-extension: k bases trimmed from upstream M (v2.9.9)
        ]
        f.write('\t'.join(header) + '\n')

        # Write results
        for result in results:
            _pt = result.get('pt_tag')
            _ps = result.get('polya_score')
            row = [
                result['read_id'],
                result['chrom'],
                result['strand'],
                str(result['original_3prime']),
                str(result['corrected_3prime']),
                str(result.get('five_prime_position', '')),  # 5' end (TSS)
                '1' if result.get('five_prime_rescued') else '0',  # 5' rescue flag
                result.get('five_prime_exon_cigar') or '',  # exon CIGAR for Cat3
                str(result.get('alignment_start', '')),  # Read body start
                str(result.get('alignment_end', '')),  # Read body end (exclusive)
                str(result['ambiguity_min']),
                str(result['ambiguity_max']),
                str(result['ambiguity_range']),
                str(result.get('polya_length', 0)),  # poly(A) length, default 0 if not computed
                str(result.get('aligned_a_length', 0)),  # Aligned A's
                str(result.get('soft_clip_a_length', 0)),  # Soft-clipped A's
                result.get('junctions_str', ''),  # Junctions as semicolon-separated string
                str(result.get('n_junctions', 0)),  # Number of junctions
                str(result.get('five_prime_soft_clip_length', 0)),  # 5' soft clip
                str(result.get('three_prime_soft_clip_length', 0)),  # 3' soft clip
                str(result.get('mapq', 0)),  # Mapping quality
                ','.join(result['correction_applied']) if result['correction_applied'] else 'none',
                result['confidence'],
                ','.join(result['qc_flags']),
                f"{result.get('fraction', 1.0):.6f}",
                result.get('gene_id') or '',  # Per-read gene attribution (empty if not computed)
                str(_pt) if _pt is not None else '',
                f'{_ps:.4f}' if _ps is not None else '',
                result.get('polya_source', 'none'),
                str(result.get('sc_homopolymer_extension', 0)),  # Cat2 CIGAR surgery
                result.get('sc_rescued_seq', ''),
                str(result.get('sc_original_softclip_len', 0)),
                str(result.get('five_prime_intron_clip_pos', -1)),  # Case 4 BAM clip
                str(result.get('oc_homopolymer_extension', 0)),  # over-call rescue D
                str(result.get('oc_overcall_count', 0)),          # over-call rescue I
                result.get('oc_terminal_base', ''),               # over-call rescue terminal base
                str(result.get('five_prime_upstream_trim', 0)),   # Cat3 equivalence-extension trim
            ]
            f.write('\t'.join(row) + '\n')

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
