#!/usr/bin/env python3
"""
Test variant-aware homopolymer rescue on wt_rep1.

This script demonstrates the two-pass approach:
1. First pass: Scan all reads to build mismatch frequency map
2. Second pass: Apply rescue with variant filtering

Author: Kevin R. Roy
Date: 2026-03-26
"""

import sys
from pathlib import Path
import pysam

# Add rectify to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.core.indel_corrector import (
    VariantAwareHomopolymerRescue,
    process_bam_with_variant_aware_rescue,
    correct_3prime_position,
)
from rectify.utils.genome import load_genome


def main():
    # Paths
    bam_path = "/path/to/wt_rep1.sorted.bam"                      # set to your BAM
    genome_path = "/path/to/S288C_reference_sequence.chrnames.fsa"  # set to your genome
    output_dir = Path("/path/to/variant_aware_test")               # set to your output dir
    output_dir.mkdir(parents=True, exist_ok=True)

    variants_path = output_dir / "potential_variants.tsv"

    print("=" * 70)
    print("Testing Variant-Aware Homopolymer Rescue")
    print("=" * 70)
    print(f"\nBAM: {bam_path}")
    print(f"Genome: {genome_path}")
    print(f"Output: {output_dir}")

    # Load genome
    print("\n[1/4] Loading genome...")
    genome = load_genome(genome_path)
    print(f"  Loaded {len(genome)} chromosomes")

    # Initialize rescue object
    print("\n[2/4] Initializing variant-aware rescue...")
    rescue = VariantAwareHomopolymerRescue(
        min_variant_fraction=0.8,       # 80% threshold for variant call
        min_reads_for_variant_call=5,   # Need at least 5 reads
        min_homopolymer_len=4,          # Require 4bp homopolymer
        max_rescue_bases=3,             # Max 3bp rescue
    )

    # First pass: Scan reads
    print("\n[3/4] First pass: Scanning reads to build mismatch frequencies...")
    n_scanned = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            strand = '-' if read.is_reverse else '+'
            rescue.scan_read(read, strand, genome, end='3prime')
            n_scanned += 1

            if n_scanned % 100000 == 0:
                print(f"  Scanned {n_scanned:,} reads...")

    print(f"  Total scanned: {n_scanned:,} reads")

    # Finalize scan
    rescue.finalize_scan()

    # Get statistics
    stats = rescue.get_statistics()
    print(f"\n  Scan Statistics:")
    print(f"    Positions with mismatches inside homopolymers: {stats['total_positions_scanned']:,}")
    print(f"    Positions with sufficient coverage (≥5 reads): {stats['positions_with_sufficient_coverage']:,}")
    print(f"    Potential variants detected (≥80% mismatch): {stats['potential_variants_detected']:,}")

    # Get and save potential variants
    variants = rescue.get_potential_variants()
    if variants:
        print(f"\n  Writing {len(variants)} potential variants to {variants_path}")
        with open(variants_path, 'w') as f:
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

        # Show first few variants
        print("\n  Sample potential variants:")
        print("  " + "-" * 80)
        for v in variants[:10]:
            print(f"    {v['chrom']}:{v['position']} "
                  f"ref={v['ref_base']} homopoly={v['homopolymer_base']} "
                  f"mismatch={v['mismatch_fraction']:.1%} ({v['total_reads']} reads) "
                  f"dominant={v['dominant_mismatch_base']} ({v['dominant_mismatch_fraction']:.1%})")
        if len(variants) > 10:
            print(f"    ... and {len(variants) - 10} more")
    else:
        print("\n  No potential variants detected (all mismatches look like basecalling errors)")

    # Second pass: Apply rescue with variant filtering (sample)
    print("\n[4/4] Second pass: Testing rescue with variant filtering (first 10,000 reads)...")

    n_tested = 0
    n_rescued = 0
    n_skipped_variant = 0
    n_skipped_too_large = 0
    n_no_rescue_needed = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            strand = '-' if read.is_reverse else '+'
            result = rescue.rescue_with_variant_filter(read, strand, genome, end='3prime')

            if result is None:
                n_no_rescue_needed += 1
            elif result['variant_check'] == 'SKIPPED_LIKELY_VARIANT':
                n_skipped_variant += 1
            elif result['variant_check'] == 'RESCUE_TOO_LARGE':
                n_skipped_too_large += 1
            elif result['variant_check'] == 'RESCUED':
                n_rescued += 1

            n_tested += 1
            if n_tested >= 10000:
                break

    print(f"\n  Rescue Results (first 10,000 reads):")
    print(f"    No rescue needed:       {n_no_rescue_needed:,} ({100*n_no_rescue_needed/n_tested:.1f}%)")
    print(f"    Rescued:                {n_rescued:,} ({100*n_rescued/n_tested:.1f}%)")
    print(f"    Skipped (likely SNP):   {n_skipped_variant:,} ({100*n_skipped_variant/n_tested:.1f}%)")
    print(f"    Skipped (too large):    {n_skipped_too_large:,} ({100*n_skipped_too_large/n_tested:.1f}%)")

    print("\n" + "=" * 70)
    print("Test complete!")
    print("=" * 70)

    if variants:
        print(f"\nPotential variants saved to: {variants_path}")
        print("These positions should be reviewed - they may be:")
        print("  - True SNPs in your strain vs reference")
        print("  - Reference assembly errors")
        print("  - Systematic sequencing artifacts")


if __name__ == "__main__":
    main()
