#!/usr/bin/env python3
"""
RECTIFY run command - all-in-one pipeline.

Runs both 'correct' and 'analyze' in sequence.

Author: Kevin R. Roy
Date: 2026-03-18
"""

import argparse
import sys
from pathlib import Path
from typing import Optional


def run(args: argparse.Namespace) -> None:
    """
    Run complete RECTIFY pipeline: correct + analyze.

    Args:
        args: Parsed command-line arguments
    """
    from . import correct_command
    from .analyze_command import run_analyze

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Output file for corrected positions
    corrected_tsv = output_dir / 'corrected_3ends.tsv'

    print("=" * 70)
    print("RECTIFY: Complete Pipeline")
    print("=" * 70)

    # =========================================================================
    # Step 1: Correct 3' end positions
    # =========================================================================
    print("\n[Step 1/2] Correcting 3' end positions...")
    print("-" * 50)

    # Build correct_args namespace
    correct_args = argparse.Namespace(
        bam=args.bam,
        genome=args.genome,
        annotation=args.annotation,
        output=corrected_tsv,
        netseq_dir=getattr(args, 'netseq_dir', None),
        aligner=getattr(args, 'aligner', 'minimap2'),
        polya_sequenced=getattr(args, 'polya_sequenced', False),
        threads=getattr(args, 'threads', 4),
        # Default values for other correct options
        min_mapq=10,
        skip_secondary=True,
        skip_supplementary=True,
        skip_ag_check=False,
        skip_atract=False,
        max_downstream_a=20,
        chunk_size=10000,
        debug=False,
    )

    try:
        correct_command.run(correct_args)
    except Exception as e:
        print(f"\nError in correction step: {e}", file=sys.stderr)
        sys.exit(1)

    if not corrected_tsv.exists():
        print(f"\nError: Corrected output file not created: {corrected_tsv}", file=sys.stderr)
        sys.exit(1)

    print(f"\nCorrection complete: {corrected_tsv}")

    # =========================================================================
    # Step 2: Analyze results
    # =========================================================================
    print("\n[Step 2/2] Analyzing results...")
    print("-" * 50)

    # Build analyze_args namespace
    analyze_args = argparse.Namespace(
        input=corrected_tsv,
        annotation=args.annotation,
        output_dir=output_dir,
        genome=args.genome,
        reference=getattr(args, 'reference', None),
        manifest=getattr(args, 'manifest', None),
        go_annotations=getattr(args, 'go_annotations', None),
        threads=getattr(args, 'threads', 4),
        # Default values for analyze options
        sample_column='sample',
        count_column=None,
        cluster_distance=25,
        min_reads=5,
        run_deseq2=True,
        run_motif=True,
        sample_sets=None,
        control=None,
        skip_browser_plots=False,
        max_browser_genes=10,
        top_genes=50,
    )

    try:
        exit_code = run_analyze(analyze_args)
        if exit_code != 0:
            print(f"\nAnalysis completed with warnings (exit code: {exit_code})")
    except Exception as e:
        print(f"\nError in analysis step: {e}", file=sys.stderr)
        sys.exit(1)

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print("\nKey outputs:")
    print(f"  - Corrected 3' ends: {corrected_tsv}")
    print(f"  - Cluster counts:    {output_dir / 'tables' / 'cluster_counts.tsv'}")
    print(f"  - DESeq2 results:    {output_dir / 'tables' / 'deseq2_genes_*.tsv'}")
    print(f"  - PCA plot:          {output_dir / 'plots' / 'pca.png'}")
    print(f"  - Motif results:     {output_dir / 'motifs' / ''}")
