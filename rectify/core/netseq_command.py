"""
NET-seq command for RECTIFY CLI.

Usage:
    rectify netseq input.bam --genome sacCer3.fa --gff genes.gff -o output/

    # With deconvolution disabled (raw 3' ends only)
    rectify netseq input.bam --genome sacCer3.fa --no-deconvolution -o output/

    # Process multiple samples
    rectify netseq sample1.bam sample2.bam --genome sacCer3.fa -o output/

Author: Kevin R. Roy
"""

import sys
from pathlib import Path
from typing import List, Optional

from .netseq_bam_processor import process_netseq_bam, aggregate_positions
from .netseq_deconvolution import deconvolve_all_regions
from .netseq_output import export_netseq_results, write_exclusion_stats
from .exclusion_regions import ExclusionRegionDetector
from ..utils.genome import load_genome


def run_netseq(args) -> int:
    """
    Execute NET-seq processing pipeline.

    Args:
        args: Parsed command line arguments

    Returns:
        Exit code (0 for success)
    """
    print("=" * 60)
    print("RECTIFY NET-seq Processing Pipeline")
    print("=" * 60)

    # Validate inputs
    input_paths = args.input
    output_dir = Path(args.output_dir)
    genome_path = Path(args.genome)

    if not genome_path.exists():
        print(f"Error: Genome file not found: {genome_path}")
        return 1

    for input_path in input_paths:
        if not Path(input_path).exists():
            print(f"Error: Input file not found: {input_path}")
            return 1

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load genome (needed for deconvolution)
    genome = None
    if not args.no_deconvolution:
        print(f"\nLoading genome: {genome_path}")
        genome = load_genome(str(genome_path))
        print(f"  Loaded {len(genome)} chromosomes")

    # Set up exclusion regions
    exclusion_detector = None
    exclude_rdna = not getattr(args, 'include_rdna', False)
    exclude_pol3 = not getattr(args, 'include_pol3', False)

    if exclude_rdna or exclude_pol3:
        print("\nSetting up exclusion regions...")
        exclusion_detector = ExclusionRegionDetector(flanking_bp=args.pol3_flanking)

        # Load from GFF if provided
        gff_path = getattr(args, 'gff', None)
        if gff_path and Path(gff_path).exists():
            n_regions = exclusion_detector.load_from_gff(
                Path(gff_path),
                exclude_tRNA=exclude_pol3,
                exclude_snRNA=exclude_pol3,
                exclude_rDNA=exclude_rdna,
                exclude_mito=args.exclude_mito,
            )
            print(f"  Loaded {n_regions} exclusion regions from GFF")
        else:
            # Use default yeast rDNA coordinates
            if exclude_rdna:
                exclusion_detector.add_rdna_region()
                print("  Added default yeast rDNA region (chrXII:450,000-490,000)")

        # Show summary
        stats = exclusion_detector.get_stats_by_reason()
        for reason, count in stats.items():
            print(f"  {reason}: {count} regions")

    # Process each input BAM
    for input_path in input_paths:
        input_path = Path(input_path)
        sample_name = input_path.stem.replace('.bam', '').replace('_bbmap_correct', '')

        print(f"\n{'='*60}")
        print(f"Processing: {input_path.name}")
        print(f"Sample name: {sample_name}")
        print(f"{'='*60}")

        # Process BAM
        print("\n1. Processing BAM file...")
        records = list(process_netseq_bam(
            input_path,
            exclusion_detector=exclusion_detector,
            min_mapq=args.min_mapq if hasattr(args, 'min_mapq') else 0,
            min_a_fraction=0.8,
            max_reads=args.max_reads if hasattr(args, 'max_reads') else None,
            show_progress=True,
        ))

        if not records:
            print("  Warning: No records generated!")
            continue

        # Aggregate positions
        print("\n2. Aggregating positions...")
        raw_counts = aggregate_positions(iter(records))
        print(f"  Unique positions: {len(raw_counts):,}")

        # Deconvolution
        deconv_counts = None
        if not args.no_deconvolution and genome:
            print("\n3. Applying deconvolution...")
            deconv_counts, deconv_results = deconvolve_all_regions(
                raw_counts,
                genome,
                min_downstream_a=args.min_atract_length if hasattr(args, 'min_atract_length') else 3,
                show_progress=True,
            )
        else:
            print("\n3. Skipping deconvolution (--no-deconvolution)")

        # Export outputs
        print("\n4. Exporting results...")
        output_formats = args.output_format if hasattr(args, 'output_format') else ['parquet', 'bedgraph']
        normalize_rpm = not getattr(args, 'no_rpm_normalize', False)

        outputs = export_netseq_results(
            records=records,
            raw_counts=raw_counts,
            deconv_counts=deconv_counts,
            output_dir=output_dir,
            sample_name=sample_name,
            export_parquet='parquet' in output_formats,
            export_bedgraph='bedgraph' in output_formats,
            export_bigwig='bigwig' in output_formats,
            normalize_rpm=normalize_rpm,
        )

        print(f"\n  Generated {len(outputs)} output files")
        for name, path in outputs.items():
            print(f"    {name}: {path.name}")

    print(f"\n{'='*60}")
    print("NET-seq processing complete!")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")

    return 0


def add_netseq_parser(subparsers) -> None:
    """Add the netseq subcommand parser."""
    parser = subparsers.add_parser(
        'netseq',
        help='Process NET-seq BAM files to extract 3\' end positions',
        description=(
            'Extract and deconvolve 3\' end positions from NET-seq data. '
            'Outputs parquet files with per-read records and RPM-normalized '
            'bedgraph files with position-level signal.'
        ),
    )

    # Required arguments
    parser.add_argument(
        'input',
        nargs='+',
        type=Path,
        help='Input BAM file(s)',
    )

    parser.add_argument(
        '--genome', '-g',
        type=Path,
        required=True,
        help='Reference genome FASTA (required for deconvolution)',
    )

    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory',
    )

    # Annotation for exclusion
    parser.add_argument(
        '--gff',
        type=Path,
        help='GFF annotation file for exclusion region detection',
    )

    # Exclusion options
    excl_group = parser.add_argument_group('Exclusion regions')
    excl_group.add_argument(
        '--include-rdna',
        action='store_true',
        help='Include rDNA locus (default: exclude)',
    )
    excl_group.add_argument(
        '--include-pol3',
        action='store_true',
        help='Include Pol III genes (tRNAs, SNR6, etc.) (default: exclude)',
    )
    excl_group.add_argument(
        '--pol3-flanking',
        type=int,
        default=100,
        help='Flanking bp around Pol III genes to exclude (default: 100)',
    )
    excl_group.add_argument(
        '--exclude-mito',
        action='store_true',
        default=False,
        help='Exclude mitochondrial genome (default: False)',
    )

    # Deconvolution options
    deconv_group = parser.add_argument_group('Deconvolution')
    deconv_group.add_argument(
        '--no-deconvolution',
        action='store_true',
        help='Disable NNLS deconvolution (output raw 3\' ends only)',
    )
    deconv_group.add_argument(
        '--min-atract-length',
        type=int,
        default=3,
        help='Minimum downstream A\'s to trigger deconvolution (default: 3)',
    )

    # Output options
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument(
        '--output-format',
        nargs='+',
        choices=['parquet', 'bedgraph', 'bigwig', 'tsv'],
        default=['parquet', 'bedgraph'],
        help='Output formats to generate (default: parquet bedgraph)',
    )
    output_group.add_argument(
        '--no-rpm-normalize',
        action='store_true',
        help='Disable RPM normalization for bedgraph/bigwig',
    )

    # Processing options
    parser.add_argument(
        '--min-mapq',
        type=int,
        default=0,
        help='Minimum mapping quality (default: 0)',
    )

    parser.add_argument(
        '--max-reads',
        type=int,
        help='Maximum reads to process (for testing)',
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output',
    )

    parser.set_defaults(func=run_netseq)
