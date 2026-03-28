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
    """
    Comprehensive statistics tracking for RECTIFY processing pipeline.

    Tracks read flow through each filtering and correction stage.
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

    # A-tract ambiguity (downstream A detection)
    ends_with_downstream_A: int = 0  # Has >=1 A in downstream 10bp window
    ends_ambiguous_atract: int = 0   # Has ambiguity range > 0

    # Correction breakdown
    ends_corrected_indel: int = 0
    ends_shifted_atract_walking: int = 0  # Position changed due to A-tract
    ends_refined_netseq: int = 0

    # Poly(A) tail length statistics
    reads_with_polya: int = 0        # Reads with detected poly(A) tail
    polya_length_sum: int = 0        # Sum of poly(A) lengths (for mean calculation)
    polya_length_max: int = 0        # Maximum observed poly(A) length

    # AG mispriming flags
    ends_flagged_ag_mispriming: int = 0
    ends_ag_high_confidence: int = 0
    ends_ag_medium_confidence: int = 0

    # Final confidence
    confidence_high: int = 0
    confidence_medium: int = 0
    confidence_low: int = 0

    # Position changes
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
            'ends_shifted_atract_walking': self.ends_shifted_atract_walking,
            'ends_refined_netseq': self.ends_refined_netseq,
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
        self.ends_shifted_atract_walking += other.ends_shifted_atract_walking
        self.ends_refined_netseq += other.ends_refined_netseq
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
        """Update stats from a single read correction result."""
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
        if 'atract_ambiguity' in corrections:
            # Check if position actually shifted
            if result.get('corrected_3prime') != result.get('original_3prime'):
                self.ends_shifted_atract_walking += 1
        if 'netseq_refinement' in corrections:
            self.ends_refined_netseq += 1

        # Poly(A) length statistics
        polya_length = result.get('polya_length', 0)
        if polya_length > 0:
            self.reads_with_polya += 1
            self.polya_length_sum += polya_length
            if polya_length > self.polya_length_max:
                self.polya_length_max = polya_length

        # Position shift
        if result.get('corrected_3prime') != result.get('original_3prime'):
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

        # Final confidence
        conf = result.get('confidence', 'high')
        if conf == 'high':
            self.confidence_high += 1
        elif conf == 'medium':
            self.confidence_medium += 1
        else:
            self.confidence_low += 1


def write_stats_tsv(stats: ProcessingStats, output_path: str) -> None:
    """
    Write processing statistics to TSV file.

    Args:
        stats: ProcessingStats object
        output_path: Path to output TSV file
    """
    total = stats.total_reads_in_bam
    processed = stats.reads_processed

    with open(output_path, 'w') as f:
        f.write("metric\tcount\tpercent\tdescription\n")

        # Input section
        f.write(f"total_reads_in_bam\t{stats.total_reads_in_bam}\t100.00\tTotal reads in BAM file\n")

        if total > 0:
            f.write(f"reads_unmapped\t{stats.reads_unmapped}\t{100*stats.reads_unmapped/total:.2f}\tUnmapped reads (skipped)\n")
            f.write(f"reads_secondary\t{stats.reads_secondary}\t{100*stats.reads_secondary/total:.2f}\tSecondary alignments (skipped)\n")
            f.write(f"reads_supplementary\t{stats.reads_supplementary}\t{100*stats.reads_supplementary/total:.2f}\tSupplementary alignments (skipped)\n")
            f.write(f"spikein_reads_filtered\t{stats.spikein_reads_filtered}\t{100*stats.spikein_reads_filtered/total:.2f}\tSpike-in reads filtered\n")
            f.write(f"reads_processed\t{stats.reads_processed}\t{100*stats.reads_processed/total:.2f}\tReads with 3' ends corrected\n")

        # Correction section (percentages relative to processed reads)
        if processed > 0:
            f.write(f"ends_with_downstream_A\t{stats.ends_with_downstream_A}\t{100*stats.ends_with_downstream_A/processed:.2f}\t3' ends with >=1 A in downstream 10bp window\n")
            f.write(f"ends_ambiguous_atract\t{stats.ends_ambiguous_atract}\t{100*stats.ends_ambiguous_atract/processed:.2f}\t3' ends in A-tract (ambiguity range > 0)\n")
            f.write(f"ends_corrected_indel\t{stats.ends_corrected_indel}\t{100*stats.ends_corrected_indel/processed:.2f}\t3' ends corrected for indel artifacts\n")
            f.write(f"ends_shifted_atract_walking\t{stats.ends_shifted_atract_walking}\t{100*stats.ends_shifted_atract_walking/processed:.2f}\t3' ends shifted by A-tract walking\n")
            f.write(f"ends_refined_netseq\t{stats.ends_refined_netseq}\t{100*stats.ends_refined_netseq/processed:.2f}\t3' ends refined by NET-seq\n")
            f.write(f"total_position_shifts\t{stats.total_position_shifts}\t{100*stats.total_position_shifts/processed:.2f}\tTotal 3' ends with position changed\n")

        # Poly(A) tail statistics
        if processed > 0:
            f.write(f"reads_with_polya\t{stats.reads_with_polya}\t{100*stats.reads_with_polya/processed:.2f}\tReads with detected poly(A) tail\n")
            mean_polya = stats.polya_length_sum / stats.reads_with_polya if stats.reads_with_polya > 0 else 0
            f.write(f"polya_length_mean\t{mean_polya:.1f}\t-\tMean poly(A) tail length (bp)\n")
            f.write(f"polya_length_max\t{stats.polya_length_max}\t-\tMaximum poly(A) tail length (bp)\n")

        # QC flags section
        if processed > 0:
            f.write(f"ends_flagged_ag_mispriming\t{stats.ends_flagged_ag_mispriming}\t{100*stats.ends_flagged_ag_mispriming/processed:.2f}\t3' ends flagged for AG mispriming\n")

        # Confidence section
        if processed > 0:
            f.write(f"confidence_high\t{stats.confidence_high}\t{100*stats.confidence_high/processed:.2f}\tHigh confidence assignments\n")
            f.write(f"confidence_medium\t{stats.confidence_medium}\t{100*stats.confidence_medium/processed:.2f}\tMedium confidence assignments\n")
            f.write(f"confidence_low\t{stats.confidence_low}\t{100*stats.confidence_low/processed:.2f}\tLow confidence assignments\n")

    logger.info(f"Wrote processing statistics to {output_path}")


def generate_stats_report(stats: ProcessingStats) -> str:
    """
    Generate formatted text report from stats.

    Args:
        stats: ProcessingStats object

    Returns:
        Formatted report string
    """
    total = stats.total_reads_in_bam
    processed = stats.reads_processed

    lines = []
    lines.append("=" * 70)
    lines.append("RECTIFY Processing Statistics")
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

    # A-tract detection
    lines.append("A-tract Ambiguity Detection:")
    if processed > 0:
        lines.append(f"  With downstream A:        {stats.ends_with_downstream_A:>12,} ({100*stats.ends_with_downstream_A/processed:>5.1f}%)")
        lines.append(f"  Ambiguous (range > 0):    {stats.ends_ambiguous_atract:>12,} ({100*stats.ends_ambiguous_atract/processed:>5.1f}%)")
    lines.append("")

    # Corrections
    lines.append("3' End Corrections Applied:")
    if processed > 0:
        lines.append(f"  Indel correction:         {stats.ends_corrected_indel:>12,} ({100*stats.ends_corrected_indel/processed:>5.1f}%)")
        lines.append(f"  A-tract walking:          {stats.ends_shifted_atract_walking:>12,} ({100*stats.ends_shifted_atract_walking/processed:>5.1f}%)")
        lines.append(f"  NET-seq refinement:       {stats.ends_refined_netseq:>12,} ({100*stats.ends_refined_netseq/processed:>5.1f}%)")
        lines.append(f"  Total position shifts:    {stats.total_position_shifts:>12,} ({100*stats.total_position_shifts/processed:>5.1f}%)")
    lines.append("")

    # Poly(A) tail statistics
    lines.append("Poly(A) Tail Statistics:")
    if processed > 0:
        lines.append(f"  Reads with poly(A):       {stats.reads_with_polya:>12,} ({100*stats.reads_with_polya/processed:>5.1f}%)")
        if stats.reads_with_polya > 0:
            mean_polya = stats.polya_length_sum / stats.reads_with_polya
            lines.append(f"  Mean poly(A) length:      {mean_polya:>12.1f} bp")
            lines.append(f"  Max poly(A) length:       {stats.polya_length_max:>12,} bp")
    lines.append("")

    # QC
    lines.append("QC Flags:")
    if processed > 0:
        lines.append(f"  AG mispriming flagged:    {stats.ends_flagged_ag_mispriming:>12,} ({100*stats.ends_flagged_ag_mispriming/processed:>5.1f}%)")
    lines.append("")

    # Confidence
    lines.append("Final Confidence:")
    if processed > 0:
        lines.append(f"  High:                     {stats.confidence_high:>12,} ({100*stats.confidence_high/processed:>5.1f}%)")
        lines.append(f"  Medium:                   {stats.confidence_medium:>12,} ({100*stats.confidence_medium/processed:>5.1f}%)")
        lines.append(f"  Low:                      {stats.confidence_low:>12,} ({100*stats.confidence_low/processed:>5.1f}%)")
    lines.append("")
    lines.append("=" * 70)

    return "\n".join(lines)
