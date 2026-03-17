#!/usr/bin/env python3
"""
RECTIFY 'correct' command implementation.

Integrates all correction modules through the BAM processor.

Author: Kevin R. Roy
Date: 2026-03-09
"""

import sys
import logging
from pathlib import Path
from typing import Optional

# CRITICAL: Set thread limits BEFORE importing numpy/pandas
# This must happen before bam_processor imports numpy
from ..slurm import set_thread_limits, get_available_cpus, get_slurm_info

from . import bam_processor
from .processing_stats import write_stats_tsv, generate_stats_report
from ..utils import genome as genome_utils


def setup_logging(verbose: bool = False):
    """Configure logging for command execution."""
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def validate_inputs(args) -> dict:
    """
    Validate input files and arguments.

    Returns:
        Dict with validated paths and settings
    """
    errors = []

    # Check BAM file
    if not args.bam.exists():
        errors.append(f"BAM file not found: {args.bam}")

    # Check genome
    if not args.genome.exists():
        errors.append(f"Genome FASTA not found: {args.genome}")

    # Check annotation (optional)
    if args.annotation and not args.annotation.exists():
        errors.append(f"Annotation file not found: {args.annotation}")

    # Check NET-seq directory (optional)
    if args.netseq_dir:
        if not args.netseq_dir.exists():
            errors.append(f"NET-seq directory not found: {args.netseq_dir}")
        elif not args.netseq_dir.is_dir():
            errors.append(f"NET-seq path is not a directory: {args.netseq_dir}")

    # Check poly(A) model (optional)
    if args.polya_model:
        if not args.polya_model.exists():
            errors.append(f"Poly(A) model file not found: {args.polya_model}")

    # Check output directory exists
    if args.output:
        output_dir = args.output.parent
        if not output_dir.exists():
            errors.append(f"Output directory does not exist: {output_dir}")

    if errors:
        for error in errors:
            logging.error(error)
        sys.exit(1)

    # Determine which modules to apply
    config = {
        'bam_path': args.bam,
        'genome_path': args.genome,
        'annotation_path': args.annotation,
        'output_path': args.output,
        'apply_atract': not args.skip_atract_check,
        'apply_ag_mispriming': not args.skip_ag_check,
        'apply_polya_trim': False,  # Default False
        'apply_indel_correction': False,  # Default False
        'netseq_dir': args.netseq_dir,
        'netseq_samples': args.netseq_samples,
        'polya_model_path': args.polya_model,
        'threads': args.threads,
        'verbose': args.verbose,
    }

    # Enable poly(A) corrections if --polya-sequenced flag set
    if args.polya_sequenced:
        config['apply_polya_trim'] = not args.skip_polya_trim
        config['apply_indel_correction'] = not args.skip_indel_correction

    return config


def run(args):
    """
    Execute the 'correct' command.

    Args:
        args: Parsed command-line arguments from argparse
    """
    # Determine thread count and set limits BEFORE numpy import
    n_threads = args.threads if args.threads > 0 else get_available_cpus()
    set_thread_limits(n_threads)

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("=" * 70)
    logger.info("RECTIFY - RNA 3' End Correction Framework")
    logger.info("=" * 70)
    logger.info("")

    # Log SLURM info if available
    slurm_info = get_slurm_info()
    if slurm_info:
        logger.info(f"SLURM job ID: {slurm_info['job_id']}")
        if slurm_info.get('array_task_id'):
            logger.info(f"  Array task: {slurm_info['array_task_id']}")
        logger.info(f"  Allocated CPUs: {slurm_info.get('cpus', 'unknown')}")
        logger.info("")

    # Validate inputs
    logger.info("Validating inputs...")
    config = validate_inputs(args)

    # Log configuration
    logger.info("Configuration:")
    logger.info(f"  Input BAM:             {config['bam_path']}")
    logger.info(f"  Reference genome:      {config['genome_path']}")
    logger.info(f"  Output TSV:            {config['output_path']}")
    logger.info(f"  Threads:               {n_threads}")
    streaming_mode = getattr(args, 'streaming', False)
    logger.info(f"  Streaming mode:        {'ENABLED' if streaming_mode else 'DISABLED'}")
    logger.info("")

    logger.info("Correction modules:")
    logger.info(f"  A-tract ambiguity:     {'ENABLED' if config['apply_atract'] else 'DISABLED'}")
    logger.info(f"  AG mispriming:         {'ENABLED' if config['apply_ag_mispriming'] else 'DISABLED'}")
    logger.info(f"  Poly(A) trimming:      {'ENABLED' if config['apply_polya_trim'] else 'DISABLED'}")
    logger.info(f"  Indel correction:      {'ENABLED' if config['apply_indel_correction'] else 'DISABLED'}")
    logger.info(f"  NET-seq refinement:    {'ENABLED' if config['netseq_dir'] else 'DISABLED'}")
    logger.info("")

    # Process BAM file
    try:
        logger.info("Processing BAM file...")

        # Choose processing mode
        if streaming_mode:
            # Streaming mode - memory efficient for large BAMs
            chunk_size = getattr(args, 'chunk_size', 10000)
            stats = bam_processor.process_bam_streaming(
                bam_path=str(config['bam_path']),
                genome_path=str(config['genome_path']),
                output_path=str(config['output_path']),
                chunk_size=chunk_size,
                apply_atract=config['apply_atract'],
                apply_ag_mispriming=config['apply_ag_mispriming'],
                apply_polya_trim=config['apply_polya_trim'],
                apply_indel_correction=config['apply_indel_correction'],
                netseq_dir=str(config['netseq_dir']) if config['netseq_dir'] else None,
            )
            report = generate_stats_report(stats)
        else:
            # Standard parallel processing
            results, stats = bam_processor.process_bam_file_parallel(
                bam_path=str(config['bam_path']),
                genome_path=str(config['genome_path']),
                n_threads=n_threads,
                apply_atract=config['apply_atract'],
                apply_ag_mispriming=config['apply_ag_mispriming'],
                apply_polya_trim=config['apply_polya_trim'],
                apply_indel_correction=config['apply_indel_correction'],
                netseq_dir=str(config['netseq_dir']) if config['netseq_dir'] else None,
                output_path=str(config['output_path']) if config['output_path'] else None,
                show_progress=not args.verbose,  # Use verbose logging instead of progress bar
                return_stats=True,
            )
            report = generate_stats_report(stats)

        # Write processing statistics TSV
        if config['output_path']:
            stats_path = str(config['output_path']).replace('.tsv', '_stats.tsv')
            if stats_path == str(config['output_path']):
                stats_path = str(config['output_path']) + '_stats.tsv'
            write_stats_tsv(stats, stats_path)
            logger.info(f"Wrote processing statistics to {stats_path}")

        # Generate summary report
        logger.info("")
        logger.info("=" * 70)
        logger.info("CORRECTION SUMMARY")
        logger.info("=" * 70)

        print(report)

        # Save report if requested
        if args.report:
            logger.info(f"Saving detailed report to {args.report}...")
            with open(args.report, 'w') as f:
                f.write(report)

        logger.info("")
        logger.info("=" * 70)
        logger.info("RECTIFY completed successfully!")
        logger.info("=" * 70)

    except Exception as e:
        logger.error(f"Error during processing: {e}")
        if args.verbose:
            logger.exception("Full traceback:")
        sys.exit(1)


if __name__ == '__main__':
    # For testing
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument('bam', type=Path)
    parser.add_argument('--genome', type=Path, required=True)
    parser.add_argument('--annotation', type=Path)
    parser.add_argument('-o', '--output', type=Path)
    parser.add_argument('--polya-sequenced', action='store_true')
    parser.add_argument('--skip-atract-check', action='store_true')
    parser.add_argument('--skip-ag-check', action='store_true')
    parser.add_argument('--skip-polya-trim', action='store_true')
    parser.add_argument('--skip-indel-correction', action='store_true')
    parser.add_argument('--netseq-dir', type=Path)
    parser.add_argument('--netseq-samples', nargs='+')
    parser.add_argument('--polya-model', type=Path)
    parser.add_argument('--report', type=Path)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()
    run(args)
