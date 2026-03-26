#!/usr/bin/env python3
"""
RECTIFY command-line interface.

Provides commands for:
- correct: Correct 3' end positions in BAM files
- train-polya: Train poly(A) tail model from control data
- validate: Validate corrections against NET-seq or other ground truth
- analyze: Downstream analysis (clustering, DESeq2, PCA, motifs, etc.)

Author: Kevin R. Roy
Date: 2026-03-09
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from . import __version__


def create_parser() -> argparse.ArgumentParser:
    """Create main argument parser."""
    parser = argparse.ArgumentParser(
        prog='rectify',
        description='RECTIFY: Unified RNA 3\' End Correction Framework',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # QuantSeq (oligo-dT short-read)
  rectify correct quantseq.bam --genome sacCer3.fa --annotation genes.gtf --polya-sequenced -o corrected.tsv

  # Nanopore direct RNA-seq with NET-seq refinement
  rectify correct nanopore.bam --genome sacCer3.fa --annotation genes.gtf --polya-sequenced \\
          --aligner minimap2 --netseq-dir bigwigs/ -o corrected.tsv

  # Train poly(A) model
  rectify train-polya nanopore.bam --genome sacCer3.fa --control-sites cpa_clusters.tsv -o model.json

  # Analyze corrected positions
  rectify analyze corrected.tsv --annotation genes.gtf --run-deseq2 --genome sacCer3.fa \\
          --run-motif -o analysis_output/

Citation:
  Roy, K. R., & Chanfreau, G. F. (2019). RECTIFY: Identification and correction of mRNA
  mis-termination caused by oligo(dT)-primed internal priming. Nucleic Acids Research, 47(16), e96.
        """
    )

    parser.add_argument(
        '--version',
        action='version',
        version=f'RECTIFY {__version__}'
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # =========================================================================
    # correct command
    # =========================================================================
    correct_parser = subparsers.add_parser(
        'correct',
        help='Correct 3\' end positions in BAM or FASTQ file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    correct_parser.add_argument(
        'input',
        type=Path,
        help='Input file (BAM or FASTQ/FASTQ.GZ). FASTQ files will be aligned with minimap2.'
    )

    correct_parser.add_argument(
        '--genome',
        type=Path,
        help='Reference genome FASTA file. Required for FASTQ input unless using bundled genome.'
    )

    correct_parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output TSV file with corrected 3\' end positions'
    )

    # Optional arguments
    correct_parser.add_argument(
        '--annotation',
        type=Path,
        help='Gene annotation file (GTF/GFF) for AG mispriming context'
    )

    # Technology flags
    tech_group = correct_parser.add_argument_group('Technology settings')
    tech_group.add_argument(
        '--polya-sequenced',
        action='store_true',
        help='Poly(A) tail IS sequenced (enables poly(A) trimming and indel correction). '
             'Use for: nanopore direct RNA, Helicos, QuantSeq, etc.'
    )

    tech_group.add_argument(
        '--aligner',
        choices=['minimap2', 'bwa', 'star', 'auto'],
        default='auto',
        help='Aligner used (affects indel artifact detection)'
    )

    # Module flags
    module_group = correct_parser.add_argument_group('Module selection')
    module_group.add_argument(
        '--skip-atract-check',
        action='store_true',
        help='Skip A-tract ambiguity detection (for organisms without A-tracts)'
    )

    module_group.add_argument(
        '--skip-ag-check',
        action='store_true',
        help='Skip AG mispriming screening'
    )

    module_group.add_argument(
        '--skip-polya-trim',
        action='store_true',
        help='Skip poly(A) tail trimming (even if --polya-sequenced)'
    )

    module_group.add_argument(
        '--skip-indel-correction',
        action='store_true',
        help='Skip indel artifact correction (even if --polya-sequenced)'
    )

    module_group.add_argument(
        '--skip-variant-aware',
        action='store_true',
        help='Skip variant-aware homopolymer rescue (enabled by default). '
             'The two-pass variant-aware approach first scans all reads to '
             'identify positions where high mismatch frequency suggests true '
             'variants (not basecalling errors), then only rescues at '
             'low-frequency positions. Potential variants are written to '
             '*_potential_variants.tsv for review.'
    )

    # Poly(A) model
    polya_group = correct_parser.add_argument_group('Poly(A) tail model')
    polya_group.add_argument(
        '--polya-model',
        type=Path,
        help='Pre-trained poly(A) tail model (JSON). If not provided, uses built-in model.'
    )

    # NET-seq refinement
    netseq_group = correct_parser.add_argument_group('NET-seq refinement')
    netseq_group.add_argument(
        '--organism',
        help='Organism name for auto NET-seq (e.g., yeast, saccharomyces_cerevisiae). '
             'Bundled WT NET-seq data will be downloaded automatically if available.'
    )

    netseq_group.add_argument(
        '--netseq-dir',
        type=Path,
        help='Custom NET-seq BigWig directory (overrides bundled data). '
             'Use this for mutant-specific NET-seq data.'
    )

    netseq_group.add_argument(
        '--netseq-samples',
        nargs='+',
        help='NET-seq sample names to use (e.g., wt_2022_rep1). If not provided, auto-detect.'
    )

    # Thresholds and parameters
    param_group = correct_parser.add_argument_group('Parameters')
    param_group.add_argument(
        '--ag-threshold',
        type=float,
        default=0.65,
        help='AG-richness threshold for mispriming flagging (0.0-1.0)'
    )

    param_group.add_argument(
        '--polya-richness',
        type=float,
        default=0.8,
        help='A-richness threshold for poly(A) tail detection (0.0-1.0)'
    )

    param_group.add_argument(
        '--min-polya-length',
        type=int,
        default=15,
        help='Minimum poly(A) tail length for nanopore oligo-dT priming'
    )

    # Output options
    output_group = correct_parser.add_argument_group('Output options')
    output_group.add_argument(
        '--report',
        type=Path,
        help='Output QC report file (HTML or PDF)'
    )

    output_group.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    # Performance options
    perf_group = correct_parser.add_argument_group('Performance options')
    perf_group.add_argument(
        '-j', '--threads',
        type=int,
        default=0,
        help='Number of threads for parallel processing. '
             '0 = auto-detect from SLURM_CPUS_PER_TASK or system (default: 0)'
    )

    perf_group.add_argument(
        '--streaming',
        action='store_true',
        help='Use streaming output mode to minimize memory usage. '
             'Recommended for BAM files > 10GB. Writes directly to output file.'
    )

    perf_group.add_argument(
        '--chunk-size',
        type=int,
        default=10000,
        help='Number of reads per output chunk in streaming mode (default: 10000)'
    )

    # =========================================================================
    # train-polya command
    # =========================================================================
    train_parser = subparsers.add_parser(
        'train-polya',
        help='Train poly(A) tail model from control data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    train_parser.add_argument(
        'bam',
        type=Path,
        help='Input BAM file with soft-clipped poly(A) tails'
    )

    train_parser.add_argument(
        '--genome',
        type=Path,
        required=True,
        help='Reference genome FASTA file'
    )

    train_parser.add_argument(
        '--control-sites',
        type=Path,
        required=True,
        help='TSV file with control CPA sites (0A downstream A-count)'
    )

    train_parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output model file (JSON)'
    )

    train_parser.add_argument(
        '--min-reads',
        type=int,
        default=10,
        help='Minimum reads per control site for training'
    )

    train_parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    # =========================================================================
    # validate command
    # =========================================================================
    validate_parser = subparsers.add_parser(
        'validate',
        help='Validate corrections against ground truth',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    validate_parser.add_argument(
        'corrected',
        type=Path,
        help='Corrected 3\' ends TSV from RECTIFY'
    )

    validate_parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output validation results TSV'
    )

    # Ground truth sources (at least one required)
    truth_group = validate_parser.add_argument_group(
        'Ground truth sources (at least one required)'
    )

    truth_group.add_argument(
        '--netseq-dir',
        type=Path,
        help='NET-seq BigWig directory (optional - requires pyBigWig)'
    )

    truth_group.add_argument(
        '--netseq-samples',
        nargs='+',
        help='Specific NET-seq samples to use'
    )

    truth_group.add_argument(
        '--annotation',
        type=Path,
        help='Gene annotation GTF/GFF with known 3\' ends'
    )

    truth_group.add_argument(
        '--ground-truth',
        type=Path,
        help='TSV file with known true 3\' end positions'
    )

    # Validation parameters
    param_group = validate_parser.add_argument_group('Validation parameters')

    param_group.add_argument(
        '--tolerance',
        type=int,
        default=1,
        help='Position tolerance in bp for "correct" classification'
    )

    param_group.add_argument(
        '--min-signal',
        type=float,
        default=0.5,
        help='Minimum NET-seq signal for ground truth'
    )

    param_group.add_argument(
        '--search-window',
        type=int,
        default=10,
        help='Window size for finding nearest ground truth'
    )

    # Output options
    validate_parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose logging'
    )

    # =========================================================================
    # analyze command
    # =========================================================================
    from .core.analyze_command import create_analyze_parser
    create_analyze_parser(subparsers)

    # =========================================================================
    # export command (bedGraph/bigWig generation)
    # =========================================================================
    export_parser = subparsers.add_parser(
        'export',
        help='Export corrected 3\' ends to bedGraph/bigWig format',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    export_parser.add_argument(
        'input',
        type=Path,
        help='Corrected 3\' end TSV file from RECTIFY'
    )

    export_parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for bedGraph/bigWig files'
    )

    export_parser.add_argument(
        '--format',
        choices=['bigwig', 'bedgraph'],
        default='bigwig',
        help='Output format'
    )

    export_parser.add_argument(
        '--genome',
        type=Path,
        help='Reference genome FASTA (for chromosome sizes)'
    )

    export_parser.add_argument(
        '--chrom-sizes',
        type=Path,
        help='Chromosome sizes file (tab-separated: chrom, size)'
    )

    export_parser.add_argument(
        '--position-col',
        default='position',
        help='Column with corrected position'
    )

    export_parser.add_argument(
        '--per-replicate',
        action='store_true',
        help='Generate per-replicate files'
    )

    export_parser.add_argument(
        '--per-condition',
        action='store_true',
        help='Generate per-condition summed files'
    )

    # =========================================================================
    # batch command (multi-sample with SLURM support)
    # =========================================================================
    from .core.batch_command import create_batch_parser
    create_batch_parser(subparsers)

    # =========================================================================
    # align command (multi-aligner pipeline)
    # =========================================================================
    from .core.align_command import create_align_parser
    create_align_parser(subparsers)

    # =========================================================================
    # extract command (per-read BAM to TSV)
    # =========================================================================
    from .core.extract_command import create_extract_parser
    create_extract_parser(subparsers)

    # =========================================================================
    # aggregate command (3' ends, 5' ends, junctions)
    # =========================================================================
    from .core.aggregate_command import create_aggregate_parser
    create_aggregate_parser(subparsers)

    # =========================================================================
    # netseq command (NET-seq BAM processing)
    # =========================================================================
    from .core.netseq_command import add_netseq_parser
    add_netseq_parser(subparsers)

    # =========================================================================
    # run command (all-in-one: correct + analyze)
    # =========================================================================
    run_parser = subparsers.add_parser(
        'run',
        help='Run complete pipeline: correct 3\' ends + analyze results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    run_parser.add_argument(
        'bam',
        type=Path,
        help='Input BAM file (aligned RNA-seq reads)'
    )

    run_parser.add_argument(
        '--genome',
        type=Path,
        required=True,
        help='Reference genome FASTA file'
    )

    run_parser.add_argument(
        '--annotation',
        type=Path,
        required=True,
        help='Gene annotation file (GTF/GFF)'
    )

    run_parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory for all results'
    )

    # Optional correction arguments
    run_parser.add_argument(
        '--organism',
        default='yeast',
        help='Organism name (default: yeast). Bundled WT NET-seq used for refinement.'
    )

    run_parser.add_argument(
        '--netseq-dir',
        type=Path,
        help='Custom NET-seq directory (overrides bundled data for mutant-specific analysis)'
    )

    run_parser.add_argument(
        '--aligner',
        choices=['minimap2', 'star', 'bowtie2', 'bwa'],
        default='minimap2',
        help='Aligner used for BAM file'
    )

    run_parser.add_argument(
        '--polya-sequenced',
        action='store_true',
        help='Poly(A) tail was sequenced (not just primed)'
    )

    # Optional analysis arguments
    run_parser.add_argument(
        '--reference',
        help='Reference condition for DESeq2 (auto-detected if not specified)'
    )

    run_parser.add_argument(
        '--manifest',
        type=Path,
        help='Sample metadata TSV (sample, condition columns)'
    )

    run_parser.add_argument(
        '--go-annotations',
        type=Path,
        help='GO annotation file for enrichment analysis'
    )

    run_parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads'
    )

    return parser


def main(argv: Optional[list] = None):
    """Main entry point for RECTIFY CLI."""
    parser = create_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Import commands only when needed
    if args.command == 'correct':
        from .core import correct_command
        correct_command.run(args)
    elif args.command == 'train-polya':
        from .core import train_polya_command
        train_polya_command.run(args)
    elif args.command == 'validate':
        from .core import validate_command
        validate_command.run(args)
    elif args.command == 'analyze':
        from .core.analyze_command import run_analyze
        sys.exit(run_analyze(args))
    elif args.command == 'export':
        from .core import export_command
        sys.exit(export_command.run(args))
    elif args.command == 'run':
        from .core import run_command
        run_command.run(args)
    elif args.command == 'batch':
        from .core import batch_command
        sys.exit(batch_command.run(args))
    elif args.command == 'aggregate':
        from .core import aggregate_command
        aggregate_command.run(args)
    elif args.command == 'align':
        from .core import align_command
        align_command.run(args)
    elif args.command == 'extract':
        from .core import extract_command
        extract_command.run(args)
    elif args.command == 'netseq':
        from .core import netseq_command
        sys.exit(netseq_command.run_netseq(args))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
