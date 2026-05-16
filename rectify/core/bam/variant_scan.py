#!/usr/bin/env python3
"""
Variant-aware homopolymer rescue pre-scan.

Single-pass scan that builds a position-level variant frequency map so the
second pass can distinguish basecaller homopolymer errors from real variants.

Author: Kevin R. Roy
"""

import logging
from typing import Dict, Optional

import pysam

from ..correct.indel_corrector import VariantAwareHomopolymerRescue

logger = logging.getLogger(__name__)


def run_variant_aware_scan(
    bam_path: str,
    genome: Dict[str, str],
    min_variant_fraction: float = 0.8,
    min_reads_for_variant_call: int = 5,
    output_variants_path: Optional[str] = None,
) -> VariantAwareHomopolymerRescue:
    """
    Run first pass to scan reads and build variant frequency map.

    This enables variant-aware homopolymer rescue by identifying positions
    where high mismatch frequency indicates a true variant (not basecalling error).

    Args:
        bam_path: Path to BAM file
        genome: Pre-loaded genome dict
        min_variant_fraction: Threshold for variant call (default 0.8 = 80%)
        min_reads_for_variant_call: Minimum reads at position (default 5)
        output_variants_path: Optional path to write potential variants TSV

    Returns:
        VariantAwareHomopolymerRescue object ready for second pass
    """
    logger.info("Running variant-aware scan (first pass)...")

    rescue = VariantAwareHomopolymerRescue(
        min_variant_fraction=min_variant_fraction,
        min_reads_for_variant_call=min_reads_for_variant_call,
        min_homopolymer_len=4,
        max_rescue_bases=3,
    )

    n_scanned = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            strand = '-' if read.is_reverse else '+'
            rescue.scan_read(read, strand, genome, end='3prime')
            n_scanned += 1

            if n_scanned % 500000 == 0:
                logger.info(f"  Scanned {n_scanned:,} reads...")

    logger.info(f"  Total scanned: {n_scanned:,} reads")

    # Finalize scan
    rescue.finalize_scan()

    # Get statistics
    stats = rescue.get_statistics()
    logger.info(f"  Positions with mismatches: {stats['total_positions_scanned']:,}")
    logger.info(f"  With sufficient coverage: {stats['positions_with_sufficient_coverage']:,}")
    logger.info(f"  Potential variants: {stats['potential_variants_detected']:,}")

    # Write potential variants if requested
    if output_variants_path:
        variants = rescue.get_potential_variants()
        if variants:
            logger.info(f"  Writing {len(variants)} potential variants to {output_variants_path}")
            with open(output_variants_path, 'w') as f:
                header = ['chrom', 'position', 'ref_base', 'homopolymer_base',
                          'total_reads', 'mismatch_fraction',
                          'dominant_mismatch_base', 'dominant_mismatch_fraction']
                f.write('\t'.join(header) + '\n')
                for v in variants:
                    row = [
                        v['chrom'],
                        str(v['position']),
                        v['ref_base'],
                        v['homopolymer_base'],
                        str(v['total_reads']),
                        f"{v['mismatch_fraction']:.3f}",
                        v['dominant_mismatch_base'],
                        f"{v['dominant_mismatch_fraction']:.3f}",
                    ]
                    f.write('\t'.join(row) + '\n')

    return rescue
