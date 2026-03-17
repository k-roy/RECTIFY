#!/usr/bin/env python3
"""
RECTIFY command-line interface.

Provides commands for:
- correct: Correct 3' end positions in BAM files
- train-polya: Train poly(A) tail model from control data
- validate: Validate corrections against NET-seq or other ground truth

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
        help='Correct 3\' end positions in BAM file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    correct_parser.add_argument(
        'bam',
        type=Path,
        help='Input BAM file (aligned RNA-seq reads)'
    )

    correct_parser.add_argument(
        '--genome',
        type=Path,
        required=True,
        help='Reference genome FASTA file (indexed with .fai)'
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
        '--netseq-dir',
        type=Path,
        help='Directory containing NET-seq BigWig files (.bw) for refinement'
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
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
