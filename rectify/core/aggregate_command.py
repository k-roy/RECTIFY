"""
Aggregate Command for RECTIFY.

Aggregates BAM reads into datasets:
- 3' ends: CPA clusters with 5' gene attribution
- 5' ends: TSS clusters with 3' gene attribution
- junctions: Splice junction counts with partial evidence

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

logger = logging.getLogger(__name__)


def create_aggregate_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Create aggregate subcommand parser."""
    parser = subparsers.add_parser(
        'aggregate',
        help='Aggregate reads into 3\' end, 5\' end, and junction datasets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument(
        'bam',
        type=Path,
        help='Input BAM file (sorted and indexed)'
    )

    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for aggregated datasets'
    )

    # Aggregation mode
    mode_group = parser.add_argument_group('Aggregation mode')
    mode_group.add_argument(
        '--mode',
        choices=['3prime', '5prime', 'junctions', 'all'],
        default='all',
        help='What to aggregate'
    )

    # Reference files
    ref_group = parser.add_argument_group('Reference files')
    ref_group.add_argument(
        '--genome',
        type=Path,
        help='Reference genome FASTA (for splice site motifs)'
    )

    ref_group.add_argument(
        '--annotation',
        type=Path,
        help='Gene annotation GFF/GTF (for gene attribution)'
    )

    ref_group.add_argument(
        '--gff',
        type=Path,
        help='GFF file with intron features (for partial junction rescue)'
    )

    # Clustering parameters
    cluster_group = parser.add_argument_group('Clustering parameters')
    cluster_group.add_argument(
        '--window',
        type=int,
        default=10,
        help='Clustering window size (bp)'
    )

    cluster_group.add_argument(
        '--min-reads',
        type=int,
        default=1,
        help='Minimum reads per cluster'
    )

    # Junction-specific options
    junction_group = parser.add_argument_group('Junction options')
    junction_group.add_argument(
        '--rescue-partial',
        action='store_true',
        help='Rescue reads truncated at splice sites using soft-clip evidence'
    )

    junction_group.add_argument(
        '--ambiguous-mode',
        choices=['proportional', 'holdout'],
        default='proportional',
        help='How to handle ambiguous partial evidence'
    )

    # Output options
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument(
        '--format',
        choices=['tsv', 'parquet'],
        default='tsv',
        help='Output format'
    )

    output_group.add_argument(
        '--include-read-ids',
        action='store_true',
        help='Include read IDs in output (increases file size)'
    )

    output_group.add_argument(
        '--prefix',
        default='',
        help='Output file prefix'
    )

    # Filter options
    filter_group = parser.add_argument_group('Filter options')
    filter_group.add_argument(
        '--chrom',
        help='Process only this chromosome'
    )

    filter_group.add_argument(
        '--strand',
        choices=['+', '-'],
        help='Process only this strand'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    return parser


def run_aggregate(args: argparse.Namespace) -> int:
    """Run aggregate command."""
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Validate input
    if not args.bam.exists():
        logger.error(f"BAM file not found: {args.bam}")
        return 1

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load genome if provided
    genome = None
    if args.genome and args.genome.exists():
        logger.info(f"Loading genome from {args.genome}")
        genome = {}
        fasta = pysam.FastaFile(str(args.genome))
        for chrom in fasta.references:
            genome[chrom] = fasta.fetch(chrom)
        fasta.close()
        logger.info(f"Loaded {len(genome)} chromosomes")

    # Load annotation if provided
    annotation_df = None
    if args.annotation and args.annotation.exists():
        logger.info(f"Loading annotation from {args.annotation}")
        annotation_df = load_annotation(args.annotation)
        logger.info(f"Loaded {len(annotation_df)} features")

    # Determine output prefix
    prefix = args.prefix if args.prefix else args.bam.stem

    # Run appropriate aggregation
    modes = ['3prime', '5prime', 'junctions'] if args.mode == 'all' else [args.mode]

    for mode in modes:
        logger.info(f"Running {mode} aggregation...")

        if mode == '3prime':
            run_3prime_aggregate(args, annotation_df, prefix)
        elif mode == '5prime':
            run_5prime_aggregate(args, annotation_df, prefix)
        elif mode == 'junctions':
            run_junction_aggregate(args, genome, prefix)

    logger.info("Aggregation complete!")
    return 0


def load_annotation(path: Path) -> pd.DataFrame:
    """Load annotation from GFF/GTF file."""
    records = []

    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attrs = fields

            # Parse attributes
            attr_dict = {}
            for attr in attrs.split(';'):
                attr = attr.strip()
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
                elif ' ' in attr:
                    parts = attr.split(' ', 1)
                    if len(parts) == 2:
                        attr_dict[parts[0]] = parts[1].strip('"')

            records.append({
                'chrom': chrom,
                'feature_type': feature,
                'start': int(start) - 1,  # Convert to 0-based
                'end': int(end),  # Half-open
                'strand': strand,
                'gene_id': attr_dict.get('ID', attr_dict.get('gene_id', '')),
                'gene_name': attr_dict.get('Name', attr_dict.get('gene_name', '')),
            })

    return pd.DataFrame(records)


def run_3prime_aggregate(args: argparse.Namespace, annotation_df: Optional[pd.DataFrame], prefix: str):
    """Run 3' end aggregation."""
    from .aggregate.three_prime import aggregate_3prime_clusters, export_3prime_clusters

    if annotation_df is None:
        logger.warning("No annotation provided - gene attribution will be empty")
        annotation_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'feature_type'])

    clusters_df = aggregate_3prime_clusters(
        bam_path=str(args.bam),
        annotation_df=annotation_df,
        chrom=args.chrom,
        strand=args.strand,
        window=args.window,
        min_reads=args.min_reads,
        include_read_ids=args.include_read_ids,
    )

    # Export
    suffix = 'parquet' if args.format == 'parquet' else 'tsv'
    output_path = args.output_dir / f"{prefix}_3prime_clusters.{suffix}"
    export_3prime_clusters(clusters_df, str(output_path), format=args.format)

    logger.info(f"Exported {len(clusters_df)} 3' end clusters to {output_path}")


def run_5prime_aggregate(args: argparse.Namespace, annotation_df: Optional[pd.DataFrame], prefix: str):
    """Run 5' end aggregation."""
    from .aggregate.five_prime import aggregate_5prime_clusters, export_5prime_clusters

    if annotation_df is None:
        logger.warning("No annotation provided - gene attribution will be empty")
        annotation_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'feature_type'])

    clusters_df = aggregate_5prime_clusters(
        bam_path=str(args.bam),
        annotation_df=annotation_df,
        chrom=args.chrom,
        strand=args.strand,
        window=args.window,
        min_reads=args.min_reads,
        include_read_ids=args.include_read_ids,
    )

    # Export
    suffix = 'parquet' if args.format == 'parquet' else 'tsv'
    output_path = args.output_dir / f"{prefix}_5prime_clusters.{suffix}"
    export_5prime_clusters(clusters_df, str(output_path), format=args.format)

    logger.info(f"Exported {len(clusters_df)} 5' end clusters to {output_path}")


def run_junction_aggregate(args: argparse.Namespace, genome: Optional[dict], prefix: str):
    """Run junction aggregation."""
    from .aggregate.junctions import (
        aggregate_junctions,
        merge_with_partial_evidence,
        export_junctions,
    )

    # Basic junction aggregation
    junction_df = aggregate_junctions(
        bam_path=str(args.bam),
        genome=genome,
        min_reads=args.min_reads,
    )

    # Optionally rescue partial evidence
    if args.rescue_partial and args.gff and args.gff.exists() and genome:
        logger.info("Rescuing partial junction evidence...")

        from .terminal_exon_refiner import (
            load_splice_sites_from_gff,
            detect_partial_junction_crossings,
        )

        splice_index = load_splice_sites_from_gff(str(args.gff))

        partial_results = detect_partial_junction_crossings(
            bam_path=str(args.bam),
            genome=genome,
            splice_index=splice_index,
            min_clip_length=1,
            ambiguous_mode=args.ambiguous_mode,
        )

        junction_df = merge_with_partial_evidence(
            junction_df,
            partial_results,
            ambiguous_mode=args.ambiguous_mode,
        )

        logger.info(f"Rescued partial evidence for {(junction_df['partial_junction_reads'] > 0).sum()} junctions")

    # Export
    suffix = 'parquet' if args.format == 'parquet' else 'tsv'
    output_path = args.output_dir / f"{prefix}_junctions.{suffix}"
    export_junctions(junction_df, str(output_path), format=args.format)

    logger.info(f"Exported {len(junction_df)} junctions to {output_path}")


def run(args: argparse.Namespace):
    """Entry point for CLI."""
    sys.exit(run_aggregate(args))
