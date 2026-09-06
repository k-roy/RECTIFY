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
    # polyA_3p / polyT_5p / gene_overlap / unassigned. Appended after the
    # pre-2026-08 columns so their positions are unaffected.
    'strand_evidence',
    # Multi-aligner consensus tags, carried through verbatim from the input
    # BAM record — see CONSENSUS_TAG_COLUMNS below. Appended last so every
    # column above keeps its index.
    'consensus_aligner',
    'consensus_confidence',
    'consensus_n_agree',
    'consensus_tied',
    # Why the 5' rescue this row found was NOT drawn into the corrected BAM.
    # '' = drawn (or no rescue). One of the bam_writer REFUSAL_* tokens:
    # extend_refused / reroute_refused / noncanonical_destination. When it is
    # set, five_prime_rescued is 0, five_prime_exon_cigar is '',
    # five_prime_intron_clip_pos is -1 and the rescued junction has been dropped
    # from `junctions` — so no consumer can draw a junction the BAM lacks. The
    # token keeps the decision auditable instead of silently dropping it.
    # Appended last (the codebase convention: newest column last, every existing column keeps its absolute index).
    'five_prime_rescue_refused',
    # Provenance of the 5' rescue's landing site (2026-09-05, ISSUE-017):
    # 1 = the rescued junction is in the annotation, 0 = a novel candidate
    # (pool junction, the read's own N-op), '' = not rescued. Lets a rescue be
    # partitioned by provenance OFFLINE from the TSV — re-deriving it from
    # coordinates is the leftmost-vs-motif trap that produced ISSUE-016.
    'five_prime_landing_annotated',
    # The novel-site evidence verdict for THIS rescue's placed segment
    # (splice_aware_5prime.NOVEL_EXON_REFUSALS), '' = passed or annotated site.
    # In RECTIFY_2F_NOVEL_GATE=report mode the rescue is still drawn and the
    # token here says what refuse mode would have done; in refuse mode the
    # same token also appears in five_prime_rescue_refused.
    'five_prime_novel_evidence',
    # ISSUE-026 invariant D (2026-09-05): junction-side 5' soft-clip bases that
    # lie over exon-2 positions (the alignment starts that many bases into exon
    # 2). The writer draws them as M between the N-op and the body, so the N-op
    # ends at the reported acceptor instead of the read's live edge. 0 = none.
    'five_prime_exon2_prefix',
    # ISSUE-028 invariant E (2026-09-06): the placed 5' block's shape, so every
    # threshold question is an offline join on the TSV (RULING 1 R4-3).
    # five_prime_exon_identity = matched / (matched + mismatches) over the
    # block, 2 decimals, a mismatch inside a genome homopolymer run >= 5
    # counting half; five_prime_exon_bits = the block's evidence score in bits
    # (a match +2, a mismatch / affine gap at the anchored aligner's constants
    # over 2, homopolymer-run errors half; local_aligner.evidence_shape), 1
    # decimal. Filled for every read whose 5' block was placed and judged —
    # drawn OR refused (the refusal token sits in five_prime_rescue_refused /
    # five_prime_novel_evidence); '' when no block was placed. Appended last.
    'five_prime_exon_identity',
    'five_prime_exon_bits',
]


#: Per-read multi-aligner consensus tags, stamped on the winning record by
#: ``consensus/consensus.py:971-975`` (``Xa``/``Xc`` also by
#: ``consensus/chimeric_consensus.py:989-992``), exposed as four SEPARATE
#: columns — never collapsed into a derived concordance class, because
#: downstream browsers consume the raw values live.
#:
#: Empty string when the tag is absent.  A per-aligner BAM carries none, so a
#: correct-first run (align → correct per aligner → merge → consensus) leaves
#: all four blank; they populate when ``correct`` runs on the align-stage
#: ``<prefix>.multialigned.bam``, which ``align_command.py:1212`` writes via
#: ``run_consensus_selection`` and ``commands/run/stages.py:79-89`` hands to
#: the correct stage.
#:
#: Order matters: these are the last four columns of ``CORRECTION_TSV_HEADER``
#: (pinned by ``tests/test_consensus_tag_columns.py``).
CONSENSUS_TAG_COLUMNS = (
    'consensus_aligner',      # Xa — winning aligner; comma list when chimeric
    'consensus_confidence',   # Xc — consensus confidence level
    'consensus_n_agree',      # Xn — number of aligners agreeing
    'consensus_tied',         # Xt — tied aligners, comma list; only when tied
)

#: ``CONSENSUS_TAG_COLUMNS`` entry → BAM tag it is read from.
_CONSENSUS_COLUMN_TAGS = (
    ('consensus_aligner', 'Xa'),
    ('consensus_confidence', 'Xc'),
    ('consensus_n_agree', 'Xn'),
    ('consensus_tied', 'Xt'),
)


def consensus_tag_fields(read) -> Dict[str, str]:
    """Read the consensus tags off *read* as strings, ``''`` when absent.

    Args:
        read: A ``pysam.AlignedSegment`` (duck-typed: only ``has_tag`` /
            ``get_tag`` are used).

    Returns:
        Dict with exactly the ``CONSENSUS_TAG_COLUMNS`` keys.

    This MUST stay exception-free.  ``bam_processor.correct_read_3prime``
    builds its chimeric (``Xz=1``) row inside a ``try/except KeyError``; an
    escaping ``KeyError`` would be swallowed there and silently reroute the
    read down the non-chimeric path — a wrong row, not a missing column.
    """
    fields: Dict[str, str] = {}
    for column, tag in _CONSENSUS_COLUMN_TAGS:
        value = ''
        try:
            if read.has_tag(tag):
                raw = read.get_tag(tag)
                if raw is not None:
                    value = str(raw)
        except (AttributeError, KeyError, TypeError, ValueError):
            value = ''
        fields[column] = value
    return fields


def _consensus_cell(result: Dict, column: str) -> str:
    """One consensus-tag TSV cell — ``''`` for a missing or ``None`` value."""
    value = result.get(column)
    return '' if value is None else str(value)


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
        # Multi-aligner consensus tags — populated by
        # ``consensus_tag_fields`` at intake; '' when the input BAM had none.
        _consensus_cell(result, 'consensus_aligner'),
        _consensus_cell(result, 'consensus_confidence'),
        _consensus_cell(result, 'consensus_n_agree'),
        _consensus_cell(result, 'consensus_tied'),
        result.get('five_prime_rescue_refused', '') or '',
        _consensus_cell(result, 'five_prime_landing_annotated'),
        result.get('five_prime_novel_evidence', '') or '',
        str(result.get('five_prime_exon2_prefix', 0) or 0),
        _shape_cell(result.get('five_prime_exon_identity'), '{:.2f}'),
        _shape_cell(result.get('five_prime_exon_bits'), '{:.1f}'),
    ]


def _shape_cell(value, fmt: str) -> str:
    """One invariant-E shape cell — ``''`` when no block was placed (None)."""
    if value is None or value == '':
        return ''
    try:
        return fmt.format(float(value))
    except (TypeError, ValueError):
        return str(value)


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
