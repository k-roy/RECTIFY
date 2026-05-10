#!/usr/bin/env python3
"""
Processing statistics tracking for RECTIFY.

Provides comprehensive tracking of read flow through each filtering
and correction stage, with TSV output.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from dataclasses import dataclass
from typing import Dict
import logging

logger = logging.getLogger(__name__)


@dataclass
class ProcessingStats:
    """Comprehensive statistics tracking for RECTIFY processing pipeline.

    Tracks read flow through each filtering and correction stage.

    Mode notes
    ----------
    Several fields are only meaningful for one sequencing protocol.
    Fields that are always populated regardless of mode are unmarked.

    DRS (direct RNA-seq, default)
    ─────────────────────────────
    • reads_with_polya / polya_length_* : poly(A) tail length measured from
      the pre-trimmed sequence; meaningful only for DRS.  Will be 0 for
      oligo-dT-primed cDNA (the poly-A tail is not sequenced in that protocol).
    • ends_flagged_ag_mispriming / ends_ag_* : ALWAYS 0 in DRS mode.
      AG mispriming is a reverse-transcriptase artefact that does not occur in
      native RNA sequencing.  These counters are populated only when
      --dT-primed-cDNA is passed.

    --dT-primed-cDNA (oligo-dT-primed cDNA, e.g. QuantSeq)
    ────────────────────────────────────────────────────────
    • reads_with_polya / polya_length_* : will be 0 because the poly-A tail is
      not in the read (the oligo-dT primer hybridises at its start).
    • ends_flagged_ag_mispriming / ends_ag_* : populated; report how many reads
      show evidence of internal AG-mispriming artefacts (Module AG).

    Both modes
    ──────────
    All other fields (A-tract ambiguity, indel correction, 5' rescue, NET-seq
    refinement, confidence) are populated in both modes.
    """
    # Input counts
    total_reads_in_bam: int = 0
    reads_unmapped: int = 0
    reads_secondary: int = 0
    reads_supplementary: int = 0

    # Spike-in filtering
    spikein_reads_filtered: int = 0

    # Reads processed (after filtering)
    reads_processed: int = 0

    # A-tract ambiguity (downstream A detection) — both modes
    # ends_with_downstream_A: count of reads where the genomic 10 bp immediately
    #   downstream of the corrected 3' end contains >=1 A (plus strand) or T
    #   (minus strand).  This is a raw A-content check, not a contiguous-run
    #   test — a single A anywhere in the window suffices.
    # ends_ambiguous_atract: subset of the above where the A-content is high
    #   enough to create real positional ambiguity (ambiguity_range > 0).
    #   This is the more meaningful diagnostic; ends_with_downstream_A is the
    #   broader "potentially affected" count.
    ends_with_downstream_A: int = 0  # >=1 A in downstream 10 bp window
    ends_ambiguous_atract: int = 0   # ambiguity_range > 0 (position uncertain)

    # Correction breakdown — both modes
    ends_corrected_indel: int = 0
    ends_corrected_false_junction: int = 0  # Poly(A)-artifact N-ops removed
    ends_shifted_atract_walking: int = 0    # Position changed due to A-tract
    ends_refined_netseq: int = 0
    ends_five_prime_rescued: int = 0        # 5' end rescued at 3'SS boundary

    # Poly(A) tail length statistics — DRS ONLY
    # Will be 0 for --dT-primed-cDNA (poly-A not sequenced in that protocol).
    reads_with_polya: int = 0        # Reads with detected poly(A) tail
    polya_length_sum: int = 0        # Sum of poly(A) lengths (for mean)
    polya_length_max: int = 0        # Maximum observed poly(A) length

    # AG mispriming flags — --dT-primed-cDNA ONLY (Module AG)
    # Always 0 in DRS mode; AG mispriming is a cDNA reverse-transcriptase
    # artefact that does not occur in native RNA sequencing.
    ends_flagged_ag_mispriming: int = 0
    ends_ag_high_confidence: int = 0
    ends_ag_medium_confidence: int = 0

    # Final confidence — both modes
    confidence_high: int = 0
    confidence_medium: int = 0
    confidence_low: int = 0

    # Position changes — both modes
    total_position_shifts: int = 0

    def to_dict(self) -> Dict[str, int]:
        """Convert to dictionary."""
        return {
            'total_reads_in_bam': self.total_reads_in_bam,
            'reads_unmapped': self.reads_unmapped,
            'reads_secondary': self.reads_secondary,
            'reads_supplementary': self.reads_supplementary,
            'spikein_reads_filtered': self.spikein_reads_filtered,
            'reads_processed': self.reads_processed,
            'ends_with_downstream_A': self.ends_with_downstream_A,
            'ends_ambiguous_atract': self.ends_ambiguous_atract,
            'ends_corrected_indel': self.ends_corrected_indel,
            'ends_corrected_false_junction': self.ends_corrected_false_junction,
            'ends_shifted_atract_walking': self.ends_shifted_atract_walking,
            'ends_refined_netseq': self.ends_refined_netseq,
            'ends_five_prime_rescued': self.ends_five_prime_rescued,
            'reads_with_polya': self.reads_with_polya,
            'polya_length_sum': self.polya_length_sum,
            'polya_length_max': self.polya_length_max,
            'ends_flagged_ag_mispriming': self.ends_flagged_ag_mispriming,
            'ends_ag_high_confidence': self.ends_ag_high_confidence,
            'ends_ag_medium_confidence': self.ends_ag_medium_confidence,
            'confidence_high': self.confidence_high,
            'confidence_medium': self.confidence_medium,
            'confidence_low': self.confidence_low,
            'total_position_shifts': self.total_position_shifts,
        }

    def merge(self, other: 'ProcessingStats') -> None:
        """Merge another stats object into this one (for parallel processing)."""
        self.total_reads_in_bam += other.total_reads_in_bam
        self.reads_unmapped += other.reads_unmapped
        self.reads_secondary += other.reads_secondary
        self.reads_supplementary += other.reads_supplementary
        self.spikein_reads_filtered += other.spikein_reads_filtered
        self.reads_processed += other.reads_processed
        self.ends_with_downstream_A += other.ends_with_downstream_A
        self.ends_ambiguous_atract += other.ends_ambiguous_atract
        self.ends_corrected_indel += other.ends_corrected_indel
        self.ends_corrected_false_junction += other.ends_corrected_false_junction
        self.ends_shifted_atract_walking += other.ends_shifted_atract_walking
        self.ends_refined_netseq += other.ends_refined_netseq
        self.ends_five_prime_rescued += other.ends_five_prime_rescued
        self.reads_with_polya += other.reads_with_polya
        self.polya_length_sum += other.polya_length_sum
        self.polya_length_max = max(self.polya_length_max, other.polya_length_max)
        self.ends_flagged_ag_mispriming += other.ends_flagged_ag_mispriming
        self.ends_ag_high_confidence += other.ends_ag_high_confidence
        self.ends_ag_medium_confidence += other.ends_ag_medium_confidence
        self.confidence_high += other.confidence_high
        self.confidence_medium += other.confidence_medium
        self.confidence_low += other.confidence_low
        self.total_position_shifts += other.total_position_shifts

    def update_from_result(self, result: Dict) -> None:
        """Update stats from a single read correction result.

        When NET-seq proportional splitting produces multiple output rows for
        one read, only the first row has is_primary_result=True.  Per-read
        stats (reads_processed, A-tract, corrections, poly(A), etc.) are gated
        on is_primary_result so they count distinct reads, not output rows.
        Per-row stats (confidence) are counted for every row.
        """
        is_primary = result.get('is_primary_result', True)

        if is_primary:
            self.reads_processed += 1

            # A-tract stats
            downstream_a = result.get('downstream_a_count', 0)
            if downstream_a and downstream_a > 0:
                self.ends_with_downstream_A += 1

            if result.get('ambiguity_range', 0) > 0:
                self.ends_ambiguous_atract += 1

            # Correction breakdown
            corrections = result.get('correction_applied', [])
            if 'indel_correction' in corrections:
                self.ends_corrected_indel += 1
            if 'false_junction' in corrections:
                self.ends_corrected_false_junction += 1
            if result.get('five_prime_rescued'):
                self.ends_five_prime_rescued += 1
            if 'atract_ambiguity' in corrections:
                # Check if position actually shifted
                if result.get('corrected_3prime') != result.get('original_3prime'):
                    self.ends_shifted_atract_walking += 1
            if 'netseq_refinement' in corrections:
                self.ends_refined_netseq += 1

            # Poly(A) length statistics
            polya_length = result.get('polya_length') or 0
            if polya_length > 0:
                self.reads_with_polya += 1
                self.polya_length_sum += polya_length
                if polya_length > self.polya_length_max:
                    self.polya_length_max = polya_length

            # Position shift (skip reads already counted via atract_ambiguity branch)
            if (result.get('corrected_3prime') != result.get('original_3prime')
                    and 'atract_ambiguity' not in corrections):
                self.total_position_shifts += 1

            # AG mispriming
            qc_flags = result.get('qc_flags', [])
            for flag in qc_flags:
                if 'AG_RICH' in flag:
                    self.ends_flagged_ag_mispriming += 1
                    ag_conf = result.get('ag_confidence', '')
                    if ag_conf == 'high':
                        self.ends_ag_high_confidence += 1
                    elif ag_conf == 'medium':
                        self.ends_ag_medium_confidence += 1
                    break

        # Final confidence is per output row (each NET-seq assignment has its own)
        conf = result.get('confidence', 'high')
        if conf == 'high':
            self.confidence_high += 1
        elif conf == 'medium':
            self.confidence_medium += 1
        else:
            self.confidence_low += 1


def write_stats_tsv(
    stats: ProcessingStats,
    output_path: str,
    protocol: str = 'drs',
) -> None:
    """Write processing statistics to TSV file.

    Args:
        stats:       ProcessingStats object.
        output_path: Path to output TSV file.
        protocol:    'drs' (direct RNA-seq, default) or 'dt_cdna'
                     (oligo-dT-primed cDNA, i.e. --dT-primed-cDNA).
                     Controls which mode-specific sections are included and
                     what the comment header says.
    """
    total = stats.total_reads_in_bam
    processed = stats.reads_processed
    is_drs = protocol not in ('dt_cdna', 'ont_cdna')
    is_ont_cdna = protocol == 'ont_cdna'

    with open(output_path, 'w') as f:
        # File-level comment header explaining mode-specific fields.
        # Lines starting with '#' are ignored by pandas read_csv / read_table
        # but are preserved in plain-text viewers.
        f.write("# RECTIFY correction statistics\n")
        if protocol == 'dt_cdna':
            f.write(
                "# Protocol: oligo-dT-primed cDNA (--dT-primed-cDNA, e.g. QuantSeq)\n"
                "#   reads_with_polya / polya_length_*  — always 0; poly-A tail is not sequenced in this protocol\n"
                "#   ends_flagged_ag_mispriming          — meaningful; reports internal AG-mispriming artefacts (Module AG)\n"
            )
        elif is_ont_cdna:
            f.write(
                "# Protocol: ONT PCR-cDNA (--ONT-cDNA, e.g. SQK-PCB114)\n"
                "#   reads_with_polya / polya_length_*  — poly(A) tail from 3' soft-clip (meaningful)\n"
                "#   ends_flagged_ag_mispriming          — always 0; no oligo-dT priming step in PCR-cDNA\n"
            )
        else:
            f.write(
                "# Protocol: DRS (direct RNA-seq)\n"
                "#   reads_with_polya / polya_length_*  — poly(A) tail length from pre-trimmed sequence (meaningful)\n"
                "#   ends_flagged_ag_mispriming          — always 0; AG mispriming does not occur in native RNA sequencing\n"
            )
        f.write(
            "#   ends_with_downstream_A              — reads with >=1 A in the 10 bp downstream window (both protocols)\n"
            "#   ends_ambiguous_atract               — subset: positional ambiguity range > 0 (the more diagnostic metric)\n"
            "# All percentages are relative to reads_processed unless noted otherwise.\n"
            "#\n"
        )
        f.write("metric\tcount\tpercent\tdescription\n")

        # Input section
        f.write(f"total_reads_in_bam\t{stats.total_reads_in_bam}\t100.00\tTotal reads in BAM file\n")

        if total > 0:
            f.write(f"reads_unmapped\t{stats.reads_unmapped}\t{100*stats.reads_unmapped/total:.2f}\tUnmapped reads (skipped)\n")
            f.write(f"reads_secondary\t{stats.reads_secondary}\t{100*stats.reads_secondary/total:.2f}\tSecondary alignments (skipped)\n")
            f.write(f"reads_supplementary\t{stats.reads_supplementary}\t{100*stats.reads_supplementary/total:.2f}\tSupplementary alignments (skipped)\n")
            f.write(f"spikein_reads_filtered\t{stats.spikein_reads_filtered}\t{100*stats.spikein_reads_filtered/total:.2f}\tSpike-in reads filtered\n")
            f.write(f"reads_processed\t{stats.reads_processed}\t{100*stats.reads_processed/total:.2f}\tReads with 3' ends corrected\n")

        # A-tract ambiguity (both protocols)
        if processed > 0:
            f.write(
                f"ends_with_downstream_A\t{stats.ends_with_downstream_A}\t"
                f"{100*stats.ends_with_downstream_A/processed:.2f}\t"
                f"3' ends with >=1 A in downstream 10 bp window (raw proximity; see ends_ambiguous_atract for actionable subset)\n"
            )
            f.write(
                f"ends_ambiguous_atract\t{stats.ends_ambiguous_atract}\t"
                f"{100*stats.ends_ambiguous_atract/processed:.2f}\t"
                f"3' ends where A-content creates real positional uncertainty (ambiguity_range > 0)\n"
            )

        # Correction breakdown (both protocols)
        if processed > 0:
            f.write(f"ends_corrected_indel\t{stats.ends_corrected_indel}\t{100*stats.ends_corrected_indel/processed:.2f}\t3' ends corrected for indel artifacts\n")
            f.write(f"ends_shifted_atract_walking\t{stats.ends_shifted_atract_walking}\t{100*stats.ends_shifted_atract_walking/processed:.2f}\t3' ends shifted by A-tract walking\n")
            f.write(f"ends_refined_netseq\t{stats.ends_refined_netseq}\t{100*stats.ends_refined_netseq/processed:.2f}\t3' ends refined by NET-seq\n")
            f.write(f"total_position_shifts\t{stats.total_position_shifts}\t{100*stats.total_position_shifts/processed:.2f}\tTotal 3' ends with position changed\n")

        # Poly(A) tail statistics — DRS and ONT PCR-cDNA; always 0 for oligo-dT cDNA
        if processed > 0:
            polya_note = " [always 0 for --dT-primed-cDNA]" if protocol == 'dt_cdna' else ""
            f.write(f"reads_with_polya\t{stats.reads_with_polya}\t{100*stats.reads_with_polya/processed:.2f}\tReads with detected poly(A) tail{polya_note}\n")
            mean_polya = stats.polya_length_sum / stats.reads_with_polya if stats.reads_with_polya > 0 else 0
            f.write(f"polya_length_mean\t{mean_polya:.1f}\t-\tMean poly(A) tail length (bp){polya_note}\n")
            f.write(f"polya_length_max\t{stats.polya_length_max}\t-\tMaximum poly(A) tail length (bp){polya_note}\n")

        # AG mispriming — --dT-primed-cDNA only; always 0 for DRS and ONT PCR-cDNA
        if processed > 0:
            ag_note = " [always 0 for DRS/ONT PCR-cDNA; AG mispriming is an oligo-dT artefact]" if protocol != 'dt_cdna' else ""
            f.write(
                f"ends_flagged_ag_mispriming\t{stats.ends_flagged_ag_mispriming}\t"
                f"{100*stats.ends_flagged_ag_mispriming/processed:.2f}\t"
                f"3' ends flagged for AG mispriming (Module AG){ag_note}\n"
            )

        # Confidence (both protocols; denominator = output rows, not reads_processed,
        # because NET-seq can split one read into multiple assignment rows)
        conf_total = stats.confidence_high + stats.confidence_medium + stats.confidence_low
        if conf_total > 0:
            f.write(f"confidence_high\t{stats.confidence_high}\t{100*stats.confidence_high/conf_total:.2f}\tHigh confidence assignments (% of output rows)\n")
            f.write(f"confidence_medium\t{stats.confidence_medium}\t{100*stats.confidence_medium/conf_total:.2f}\tMedium confidence assignments (% of output rows)\n")
            f.write(f"confidence_low\t{stats.confidence_low}\t{100*stats.confidence_low/conf_total:.2f}\tLow confidence assignments (% of output rows)\n")

    logger.info(f"Wrote processing statistics to {output_path}")


def write_consensus_stats_tsv(stats: dict, output_path: str) -> None:
    """Write consensus aligner selection statistics to TSV file.

    Args:
        stats: Stats dict returned by run_consensus_selection().
               Expected keys: total_reads, consensus_high, consensus_medium,
               consensus_low, 5prime_rescued, tied_score, by_aligner,
               chimeric_reads (optional).
        output_path: Destination TSV path.
    """
    total = stats.get('total_reads', 0)
    by_aligner = dict(stats.get('by_aligner', {}))

    with open(output_path, 'w') as f:
        f.write("metric\tcount\tpercent\tdescription\n")

        # Overall totals
        f.write(f"total_reads\t{total}\t100.00\tTotal reads written to consensus BAM\n")

        # Per-aligner wins
        for aligner, count in sorted(by_aligner.items(), key=lambda x: -x[1]):
            pct = 100.0 * count / total if total else 0.0
            f.write(f"aligner_wins_{aligner}\t{count}\t{pct:.2f}\tReads where {aligner} was selected\n")

        # Confidence
        if total:
            high   = stats.get('consensus_high', 0)
            medium = stats.get('consensus_medium', 0)
            low    = stats.get('consensus_low', 0)
            f.write(f"confidence_high\t{high}\t{100.0*high/total:.2f}\tHigh-confidence selections\n")
            f.write(f"confidence_medium\t{medium}\t{100.0*medium/total:.2f}\tMedium-confidence selections\n")
            f.write(f"confidence_low\t{low}\t{100.0*low/total:.2f}\tLow-confidence selections\n")

        # Other metrics
        rescued = stats.get('5prime_rescued', 0)
        tied    = stats.get('tied_score', 0)
        chimeric = stats.get('chimeric_reads', 0)
        if total:
            f.write(f"5prime_rescued\t{rescued}\t{100.0*rescued/total:.2f}\tReads rescued via 5' soft-clip splice detection\n")
            f.write(f"tied_score\t{tied}\t{100.0*tied/total:.2f}\tReads requiring tiebreaker (equal junction score)\n")
            if chimeric:
                f.write(f"chimeric_reads\t{chimeric}\t{100.0*chimeric/total:.2f}\tChimeric reads (multi-aligner segments)\n")

    logger.info(f"Wrote consensus aligner stats to {output_path}")


def generate_stats_report(stats: ProcessingStats, protocol: str = 'drs') -> str:
    """Generate formatted text report from stats.

    Args:
        stats:    ProcessingStats object.
        protocol: 'drs' (direct RNA-seq, default) or 'dt_cdna'
                  (--dT-primed-cDNA).  Adds a note next to sections that are
                  inactive for the current protocol.

    Returns:
        Formatted report string.
    """
    total = stats.total_reads_in_bam
    processed = stats.reads_processed
    is_drs = protocol not in ('dt_cdna', 'ont_cdna')
    is_ont_cdna = protocol == 'ont_cdna'

    lines = []
    lines.append("=" * 70)
    lines.append("RECTIFY Processing Statistics")
    if protocol == 'dt_cdna':
        proto_label = "oligo-dT-primed cDNA (--dT-primed-cDNA)"
    elif is_ont_cdna:
        proto_label = "ONT PCR-cDNA (--ONT-cDNA)"
    else:
        proto_label = "DRS (direct RNA-seq)"
    lines.append(f"Protocol: {proto_label}")
    lines.append("=" * 70)
    lines.append("")

    # Input
    lines.append("Input:")
    lines.append(f"  Total reads in BAM:       {total:>12,}")
    if total > 0:
        lines.append(f"  Unmapped (skipped):       {stats.reads_unmapped:>12,} ({100*stats.reads_unmapped/total:>5.1f}%)")
        lines.append(f"  Secondary (skipped):      {stats.reads_secondary:>12,} ({100*stats.reads_secondary/total:>5.1f}%)")
        lines.append(f"  Spike-in (filtered):      {stats.spikein_reads_filtered:>12,} ({100*stats.spikein_reads_filtered/total:>5.1f}%)")
        lines.append(f"  Reads processed:          {processed:>12,} ({100*processed/total:>5.1f}%)")
    lines.append("")

    # A-tract detection (both protocols)
    lines.append("A-tract Ambiguity Detection:")
    if processed > 0:
        lines.append(f"  With downstream A:        {stats.ends_with_downstream_A:>12,} ({100*stats.ends_with_downstream_A/processed:>5.1f}%)")
        lines.append( "    (>=1 A in 10 bp downstream window; raw proximity measure)")
        lines.append(f"  Ambiguous (range > 0):    {stats.ends_ambiguous_atract:>12,} ({100*stats.ends_ambiguous_atract/processed:>5.1f}%)")
        lines.append( "    (positional uncertainty > 0; the more diagnostic metric)")
    lines.append("")

    # Corrections (both protocols)
    lines.append("3' End Corrections Applied:")
    if processed > 0:
        lines.append(f"  Indel correction:         {stats.ends_corrected_indel:>12,} ({100*stats.ends_corrected_indel/processed:>5.1f}%)")
        lines.append(f"  False junction removed:   {stats.ends_corrected_false_junction:>12,} ({100*stats.ends_corrected_false_junction/processed:>5.1f}%)")
        lines.append(f"  A-tract walking:          {stats.ends_shifted_atract_walking:>12,} ({100*stats.ends_shifted_atract_walking/processed:>5.1f}%)")
        lines.append(f"  NET-seq refinement:       {stats.ends_refined_netseq:>12,} ({100*stats.ends_refined_netseq/processed:>5.1f}%)")
        lines.append(f"  5' soft-clip rescued:     {stats.ends_five_prime_rescued:>12,} ({100*stats.ends_five_prime_rescued/processed:>5.1f}%)")
        lines.append(f"  Total position shifts:    {stats.total_position_shifts:>12,} ({100*stats.total_position_shifts/processed:>5.1f}%)")
    lines.append("")

    # Poly(A) tail statistics — DRS and ONT PCR-cDNA; N/A for oligo-dT cDNA
    polya_note = " [N/A for --dT-primed-cDNA]" if protocol == 'dt_cdna' else ""
    lines.append(f"Poly(A) Tail Statistics:{polya_note}")
    if processed > 0:
        if protocol == 'dt_cdna':
            lines.append("  (poly-A tail not sequenced in oligo-dT-primed cDNA — all values are 0)")
        else:
            lines.append(f"  Reads with poly(A):       {stats.reads_with_polya:>12,} ({100*stats.reads_with_polya/processed:>5.1f}%)")
            if stats.reads_with_polya > 0:
                mean_polya = stats.polya_length_sum / stats.reads_with_polya
                lines.append(f"  Mean poly(A) length:      {mean_polya:>12.1f} bp")
                lines.append(f"  Max poly(A) length:       {stats.polya_length_max:>12,} bp")
    lines.append("")

    # QC — AG mispriming is --dT-primed-cDNA only; always 0 for DRS and ONT PCR-cDNA
    ag_note = " [always 0 — AG mispriming is an oligo-dT artefact]" if protocol != 'dt_cdna' else ""
    lines.append(f"QC Flags (Module AG):{ag_note}")
    if processed > 0:
        if protocol == 'dt_cdna':
            lines.append(f"  AG mispriming flagged:    {stats.ends_flagged_ag_mispriming:>12,} ({100*stats.ends_flagged_ag_mispriming/processed:>5.1f}%)")
        else:
            lines.append("  AG mispriming flagged:               0 (not applicable)")
    lines.append("")

    # Confidence (denominator = total output rows, which can exceed reads_processed
    # when NET-seq splits one read into multiple assignments)
    conf_total = stats.confidence_high + stats.confidence_medium + stats.confidence_low
    lines.append("Final Confidence (% of output rows):")
    if conf_total > 0:
        lines.append(f"  High:                     {stats.confidence_high:>12,} ({100*stats.confidence_high/conf_total:>5.1f}%)")
        lines.append(f"  Medium:                   {stats.confidence_medium:>12,} ({100*stats.confidence_medium/conf_total:>5.1f}%)")
        lines.append(f"  Low:                      {stats.confidence_low:>12,} ({100*stats.confidence_low/conf_total:>5.1f}%)")
    lines.append("")
    lines.append("=" * 70)

    return "\n".join(lines)
