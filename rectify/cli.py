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
  # QuantSeq or dT-primed cDNA (poly-A is NOT sequenced — oligo-dT primes at poly-A start)
  rectify correct quantseq.bam --genome sacCer3.fa --annotation genes.gtf --dT-primed-cDNA -o corrected.tsv

  # Nanopore direct RNA-seq with NET-seq refinement (poly-A IS sequenced)
  rectify correct nanopore.bam --genome sacCer3.fa --annotation genes.gtf \\
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

    # correct command (parser wired in rectify.core.commands.correct_command)
    from .core.commands.correct_command import create_correct_parser
    create_correct_parser(subparsers)

    # train-polya command (parser wired in rectify.core.commands.train_polya_command)
    from .core.commands.train_polya_command import create_train_polya_parser
    create_train_polya_parser(subparsers)

    # tag-polya command (parser wired in rectify.core.commands.tag_polya_command)
    from .core.commands.tag_polya_command import create_tag_polya_parser
    create_tag_polya_parser(subparsers)

    # validate command (parser wired in rectify.core.commands.validate_command)
    from .core.commands.validate_command import create_validate_parser
    create_validate_parser(subparsers)

    # =========================================================================
    # analyze command
    # =========================================================================
    from .core.commands.analyze_command import create_analyze_parser
    create_analyze_parser(subparsers)

    # export command (parser wired in rectify.core.commands.export_command)
    from .core.commands.export_command import create_export_parser
    create_export_parser(subparsers)

    # =========================================================================
    # batch command (multi-sample with SLURM support)
    # =========================================================================
    from .core.commands.batch_command import create_batch_parser
    create_batch_parser(subparsers)

    # =========================================================================
    # test command (installation smoke-test)
    # =========================================================================
    from .core.commands.test_command import create_test_parser
    create_test_parser(subparsers)

    # =========================================================================
    # align command (multi-aligner pipeline)
    # =========================================================================
    from .core.commands.align_command import create_align_parser
    create_align_parser(subparsers)

    # =========================================================================
    # split command (FASTQ chunker for parallel array alignment)
    # =========================================================================
    from .core.commands.split_command import create_split_parser
    create_split_parser(subparsers)

    # =========================================================================
    # install-aligners command (download/compile external aligners)
    # =========================================================================
    from .core.commands.install_aligners_command import create_install_aligners_parser
    create_install_aligners_parser(subparsers)

    # =========================================================================
    # consensus command (aligner selection on pre-built BAMs)
    # =========================================================================
    from .core.commands.consensus_command import create_consensus_parser
    create_consensus_parser(subparsers)

    # =========================================================================
    # extract command (per-read BAM to TSV)
    # =========================================================================
    from .core.commands.extract_command import create_extract_parser
    create_extract_parser(subparsers)

    # =========================================================================
    # prescan command (pre-compute variant scan + junction pool for chunked correction)
    # =========================================================================
    from .core.commands.prescan_command import create_prescan_parser
    create_prescan_parser(subparsers)

    # =========================================================================
    # aggregate command (3' ends, 5' ends, junctions)
    # =========================================================================
    from .core.commands.aggregate_command import create_aggregate_parser
    create_aggregate_parser(subparsers)

    # =========================================================================
    # netseq command (NET-seq BAM processing)
    # =========================================================================
    from .core.commands.netseq_command import add_netseq_parser
    add_netseq_parser(subparsers)

    # trim-polya command (parser wired in rectify.core.commands.drs_trim_command)
    from .core.commands.drs_trim_command import create_trim_polya_parser
    create_trim_polya_parser(subparsers)

    # trim-cdna-polya command (parser wired in rectify.core.commands.cdna_trim_command)
    from .core.commands.cdna_trim_command import create_trim_cdna_polya_parser
    create_trim_cdna_polya_parser(subparsers)

    # correct-cdna command (parser wired in rectify.core.commands.cdna_correct_command)
    from .core.commands.cdna_correct_command import create_correct_cdna_parser
    create_correct_cdna_parser(subparsers)

    # cdna-analyze command (parser wired in rectify.core.commands.cdna_analyze_command)
    from .core.commands.cdna_analyze_command import create_cdna_analyze_parser
    create_cdna_analyze_parser(subparsers)

    # restore-softclip command (parser wired in rectify.core.commands.restore_polya_command)
    from .core.commands.restore_polya_command import create_restore_softclip_parser
    create_restore_softclip_parser(subparsers)

    # run-all command (parser wired in rectify.core.commands.run_command)
    from .core.commands.run_command import create_run_parser
    create_run_parser(subparsers)

    return parser


def main(argv: Optional[list] = None):
    """Main entry point for RECTIFY CLI."""
    parser = create_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Global --organism / --Scer hook: fill args.genome / args.annotation /
    # args.go_annotations from bundled data for any subcommand that opted into
    # rectify.data.add_organism_args. Skipped for commands that perform their
    # own bespoke resolution (correct, run-all, batch — each decompresses the
    # bundled FASTA / merges with NET-seq path overrides before alignment).
    _BESPOKE_RESOLVERS = {'correct', 'run-all', 'batch'}
    if args.command not in _BESPOKE_RESOLVERS:
        from .data import resolve_reference_paths
        resolve_reference_paths(args, require_genome=False, verbose=False)

    # Import commands only when needed
    if args.command == 'correct':
        from .core.commands import correct_command
        correct_command.run(args)
    elif args.command == 'train-polya':
        from .core.commands import train_polya_command
        train_polya_command.run(args)
    elif args.command == 'tag-polya':
        from .core.commands import tag_polya_command
        sys.exit(tag_polya_command.run(args))
    elif args.command == 'validate':
        from .core.commands import validate_command
        validate_command.run(args)
    elif args.command == 'analyze':
        from .core.commands.analyze_command import run_analyze
        sys.exit(run_analyze(args))
    elif args.command == 'export':
        from .core.commands import export_command
        sys.exit(export_command.run(args))
    elif args.command == 'run-all':
        from .core.commands import run_command
        run_command.run(args)
    elif args.command == 'batch':
        from .core.commands import batch_command
        sys.exit(batch_command.run(args))
    elif args.command == 'aggregate':
        from .core.commands import aggregate_command
        aggregate_command.run(args)
    elif args.command == 'prescan':
        from .core.commands import prescan_command
        sys.exit(prescan_command.run(args))
    elif args.command == 'test':
        from .core.commands import test_command
        sys.exit(test_command.run(args))
    elif args.command == 'align':
        from .core.commands import align_command
        align_command.run(args)
    elif args.command == 'split':
        from .core.commands import split_command
        split_command.run(args)
    elif args.command == 'consensus':
        from .core.commands import consensus_command
        consensus_command.run(args)
    elif args.command == 'install-aligners':
        from .core.commands import install_aligners_command
        install_aligners_command.run(args)
    elif args.command == 'extract':
        from .core.commands import extract_command
        extract_command.run(args)
    elif args.command == 'netseq':
        from .core.commands import netseq_command
        sys.exit(netseq_command.run_netseq(args))
    elif args.command == 'trim-polya':
        from .core.commands import drs_trim_command
        sys.exit(drs_trim_command.run(args))
    elif args.command == 'trim-cdna-polya':
        from .core.commands import cdna_trim_command
        sys.exit(cdna_trim_command.run(args))
    elif args.command == 'correct-cdna':
        from .core.commands import cdna_correct_command
        sys.exit(cdna_correct_command.run(args))
    elif args.command == 'cdna-analyze':
        from .core.commands import cdna_analyze_command
        sys.exit(cdna_analyze_command.run(args))
    elif args.command == 'restore-softclip':
        from .core.commands import restore_polya_command
        sys.exit(restore_polya_command.run(args))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
