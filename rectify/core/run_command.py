#!/usr/bin/env python3
"""
RECTIFY run command - all-in-one pipeline.

Runs 'correct', 'analyze', and 'junction aggregation' in sequence:
1. Correct 3' end positions (poly(A) trimming, indel correction, NET-seq refinement)
2. Analyze results (clustering, DESeq2, motifs)
3. Junction aggregation with partial rescue (requires GFF annotation)

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
    from ..data import ensure_netseq_data

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Output file for corrected positions
    corrected_tsv = output_dir / 'corrected_3ends.tsv'

    print("=" * 70)
    print("RECTIFY: Complete Pipeline")
    print("=" * 70)

    # Resolve NET-seq data (bundled or custom)
    from ..data import detect_organism
    custom_netseq = getattr(args, 'netseq_dir', None)
    resolved_organism = getattr(args, 'organism', None)

    if custom_netseq:
        print(f"\nUsing custom NET-seq data: {custom_netseq}")
        resolved_netseq_dir = custom_netseq
    else:
        # Auto-detect organism from genome/annotation if not provided
        if not resolved_organism:
            resolved_organism = detect_organism(args.genome, args.annotation)

        if resolved_organism:
            print(f"\nAuto-detected organism: {resolved_organism}")
            # Check if bundled data is available (just for info message)
            netseq_result = ensure_netseq_data(
                resolved_organism,
                auto_download=True,
                verbose=True
            )
            # 'bundled' marker means bundled data exists - let correct_command handle it
            # Only set resolved_netseq_dir if it's an actual path (custom data)
            if netseq_result and netseq_result != 'bundled':
                resolved_netseq_dir = netseq_result
            else:
                resolved_netseq_dir = None  # Let correct_command resolve bundled data
        else:
            print("\nCould not auto-detect organism. Running without NET-seq refinement.")
            print("(Provide --organism or --netseq-dir for NET-seq refinement)")
            resolved_netseq_dir = None

    # =========================================================================
    # Step 1: Correct 3' end positions
    # =========================================================================
    print("\n[Step 1/3] Correcting 3' end positions...")
    print("-" * 50)

    # Build correct_args namespace
    # Ensure all paths are Path objects (in case they're strings)
    input_path = Path(args.bam) if args.bam else None
    genome_path = Path(args.genome) if args.genome else None
    annotation_path = Path(args.annotation) if args.annotation else None
    netseq_path = Path(resolved_netseq_dir) if resolved_netseq_dir else None

    # Check for junction rescue requirements (GFF has intron features)
    has_gff = annotation_path and annotation_path.suffix.lower() in ('.gff', '.gff3')

    correct_args = argparse.Namespace(
        input=input_path,  # correct command uses 'input', not 'bam'
        genome=genome_path,
        annotation=annotation_path,
        output=corrected_tsv,
        netseq_dir=netseq_path,
        organism=resolved_organism,  # Pass organism so correct_command can resolve bundled data
        aligner=getattr(args, 'aligner', 'minimap2'),
        polya_sequenced=getattr(args, 'polya_sequenced', False),
        threads=getattr(args, 'threads', 4),
        # Default values for other correct options
        min_mapq=10,
        skip_secondary=True,
        skip_supplementary=True,
        skip_ag_check=False,
        skip_atract_check=False,
        skip_polya_trim=False,
        skip_indel_correction=False,
        polya_model=None,
        report=None,
        max_downstream_a=20,
        chunk_size=10000,
        debug=False,
        verbose=False,
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
    print("\n[Step 2/3] Analyzing results...")
    print("-" * 50)

    # Build analyze_args namespace (reuse already-converted Path objects)
    manifest_path = Path(args.manifest) if getattr(args, 'manifest', None) else None
    go_path = Path(args.go_annotations) if getattr(args, 'go_annotations', None) else None

    analyze_args = argparse.Namespace(
        input=corrected_tsv,
        annotation=annotation_path,
        output_dir=output_dir,
        genome=genome_path,
        reference=getattr(args, 'reference', None),
        manifest=manifest_path,
        go_annotations=go_path,
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
    # Step 3: Junction aggregation with partial rescue
    # =========================================================================
    junctions_dir = output_dir / 'junctions'
    junctions_tsv = junctions_dir / 'junctions.tsv'

    if genome_path and has_gff:
        print("\n[Step 3/3] Aggregating splice junctions with partial rescue...")
        print("-" * 50)

        from .aggregate.junctions import aggregate_junctions, merge_with_partial_evidence, export_junctions
        from .terminal_exon_refiner import load_splice_sites_from_gff, detect_partial_junction_crossings
        import pysam

        junctions_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Load genome for motif extraction
            print("Loading genome...")
            genome = {}
            fasta = pysam.FastaFile(str(genome_path))
            for chrom in fasta.references:
                genome[chrom] = fasta.fetch(chrom)
            fasta.close()

            # Basic junction aggregation from CIGAR N operations
            print("Aggregating junctions from alignments...")
            junction_df = aggregate_junctions(
                bam_path=str(input_path),
                genome=genome,
                min_reads=1,
            )
            print(f"  Found {len(junction_df)} junctions from CIGAR")

            # Rescue partial junction evidence using soft-clips
            print("Loading splice sites from annotation...")
            splice_index = load_splice_sites_from_gff(str(annotation_path))

            print("Detecting partial junction crossings...")
            partial_results = detect_partial_junction_crossings(
                bam_path=str(input_path),
                genome=genome,
                splice_index=splice_index,
                min_clip_length=1,
                ambiguous_mode='proportional',
            )

            n_rescued = partial_results['summary']['total_rescued']
            n_ambiguous = partial_results['summary']['total_ambiguous']
            print(f"  Rescued {n_rescued} partial crossings ({n_ambiguous} ambiguous)")

            # Merge partial evidence with junction counts
            junction_df = merge_with_partial_evidence(
                junction_df,
                partial_results,
                ambiguous_mode='proportional',
            )

            # Export
            export_junctions(junction_df, str(junctions_tsv), format='tsv')
            print(f"\nJunction aggregation complete: {junctions_tsv}")

        except Exception as e:
            print(f"\nWarning: Junction aggregation failed: {e}", file=sys.stderr)
            print("Continuing without junction output...")
    else:
        print("\n[Step 3/3] Skipping junction aggregation...")
        print("-" * 50)
        if not genome_path:
            print("  (requires --genome for splice site motif extraction)")
        if not has_gff:
            print("  (requires GFF/GFF3 annotation for intron features)")

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
    if junctions_tsv.exists():
        print(f"  - Junctions:         {junctions_tsv}")
