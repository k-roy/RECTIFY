#!/usr/bin/env python3
"""
RECTIFY 'correct' command implementation.

Integrates all correction modules through the BAM processor.

Author: Kevin R. Roy
Date: 2026-03-09
"""

import sys
import dataclasses
import logging
from pathlib import Path
from typing import Optional

import pysam

# CRITICAL: Set thread limits BEFORE importing numpy/pandas
# This must happen before bam_processor imports numpy
from ..slurm import set_thread_limits, get_available_cpus, get_slurm_info

from . import bam_processor
from .processing_stats import write_stats_tsv, generate_stats_report
from .spikein_filter import filter_spikein_reads
from ..utils import genome as genome_utils
from ..utils.provenance import init_provenance


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

    Handles preprocessing of FASTQ files and bundled genomes.

    Returns:
        Dict with validated paths and settings
    """
    from .preprocess import detect_input_type, prepare_input, prepare_bundled_genome
    from ..data import (
        normalize_organism, is_bundled_genome_available,
        detect_organism, ensure_netseq_data
    )

    errors = []

    # Check input file
    if not args.input.exists():
        errors.append(f"Input file not found: {args.input}")
        # Can't proceed without input
        if errors:
            print("Input validation errors:")
            for err in errors:
                print(f"  - {err}")
            sys.exit(1)

    # Detect input type
    input_type = detect_input_type(args.input)
    logging.info(f"Input type detected: {input_type}")

    # Handle organism and bundled data
    organism = getattr(args, 'organism', None)

    # Try to get genome path
    genome_path = getattr(args, 'genome', None)
    annotation_path = getattr(args, 'annotation', None)

    # If no genome provided, check for bundled genome
    if genome_path is None and organism:
        org = normalize_organism(organism)
        if is_bundled_genome_available(org):
            logging.info(f"Using bundled genome for {org}")
            bundled_genome, bundled_ann = prepare_bundled_genome(
                organism=org,
                output_dir=args.output.parent if args.output else None,
                verbose=True
            )
            genome_path = bundled_genome
            if annotation_path is None:
                annotation_path = bundled_ann

    # For FASTQ input, genome is required
    if input_type in ('fastq', 'fastq.gz') and genome_path is None:
        errors.append(
            "Genome required for FASTQ input. Provide --genome or --organism with bundled genome."
        )

    # Check genome exists
    if genome_path and not genome_path.exists():
        errors.append(f"Genome FASTA not found: {genome_path}")

    # Preprocess FASTQ input (align with minimap2)
    if input_type in ('fastq', 'fastq.gz') and genome_path and not errors:
        try:
            output_dir = args.output.parent if args.output else args.input.parent
            threads = getattr(args, 'threads', 4)

            bam_path, _ = prepare_input(
                input_path=args.input,
                genome_path=genome_path,
                output_dir=output_dir,
                threads=threads,
                verbose=True
            )
            # Store the BAM path for later use
            args._bam_path = bam_path
        except Exception as e:
            errors.append(f"Failed to align FASTQ: {e}")
    else:
        # Input is already BAM
        args._bam_path = args.input

    # Store resolved paths
    args._genome_path = genome_path
    args._annotation_path = annotation_path

    # Check annotation (optional) - only if custom annotation provided
    if args.annotation and not args.annotation.exists():
        errors.append(f"Annotation file not found: {args.annotation}")
    elif annotation_path:
        # Use resolved annotation path (from bundled or custom)
        pass

    # Resolve NET-seq directory (organism bundled data or custom)
    resolved_netseq_dir = None
    if getattr(args, 'netseq_dir', None):
        # Custom dir provided - validate it exists
        if not args.netseq_dir.exists():
            errors.append(f"NET-seq directory not found: {args.netseq_dir}")
        elif not args.netseq_dir.is_dir():
            errors.append(f"NET-seq path is not a directory: {args.netseq_dir}")
        else:
            resolved_netseq_dir = args.netseq_dir
    else:
        # Try to use organism for bundled NET-seq data
        if not organism and genome_path:
            # Auto-detect from files
            organism = detect_organism(genome_path, annotation_path)
            if organism:
                logging.info(f"Auto-detected organism: {organism}")

        if organism:
            try:
                result = ensure_netseq_data(
                    organism,
                    auto_download=True,
                    verbose=True
                )
                # ensure_netseq_data returns 'bundled' for organisms with
                # pre-deconvolved TSV data included in the package.
                # Encode organism so bam_processor workers can self-load.
                if result == 'bundled':
                    resolved_netseq_dir = f'bundled:{organism}'
                else:
                    resolved_netseq_dir = result
            except Exception as e:
                logging.warning(f"Could not load bundled NET-seq data: {e}")
                logging.warning("Continuing without NET-seq refinement.")

    # Store resolved path for later use
    args._resolved_netseq_dir = resolved_netseq_dir

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
    # Use resolved paths from preprocessing
    bam_path = getattr(args, '_bam_path', args.input)
    genome_path = getattr(args, '_genome_path', getattr(args, 'genome', None))
    annotation_path = getattr(args, '_annotation_path', getattr(args, 'annotation', None))

    config = {
        'bam_path': bam_path,
        'genome_path': genome_path,
        'annotation_path': annotation_path,
        'output_path': args.output,
        'apply_atract': not args.skip_atract_check,
        'apply_ag_mispriming': not args.skip_ag_check,
        'ag_threshold': getattr(args, 'ag_threshold', 0.65),
        'apply_polya_trim': False,  # Default False
        'apply_indel_correction': False,  # Default False
        'netseq_dir': getattr(args, '_resolved_netseq_dir', getattr(args, 'netseq_dir', None)),
        'netseq_samples': getattr(args, 'netseq_samples', None),
        'polya_model_path': getattr(args, 'polya_model', None),
        'threads': getattr(args, 'threads', 4),
        'verbose': getattr(args, 'verbose', False),
        'variant_aware': not getattr(args, 'skip_variant_aware', False),
        'filter_spikein': getattr(args, 'filter_spikein', None),
        'output_bam': getattr(args, 'output_bam', None),
        'corrected_bam': getattr(args, 'corrected_bam', None),
        'softclipped_bam': getattr(args, 'softclipped_bam', None),
        'bedgraph_prefix': getattr(args, 'bedgraph_prefix', None),
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
    logger.info(f"  Variant-aware rescue:  {'ENABLED' if config['variant_aware'] else 'DISABLED'}")
    logger.info(f"  NET-seq refinement:    {'ENABLED' if config['netseq_dir'] else 'DISABLED'}")
    logger.info(f"  Spike-in filter:       {'ENABLED (' + ','.join(config['filter_spikein']) + ')' if config['filter_spikein'] else 'DISABLED'}")
    logger.info("")

    # Initialize provenance tracking
    provenance = None
    if config['output_path']:
        output_dir = Path(config['output_path']).parent
        provenance = init_provenance(
            output_dir,
            description="RECTIFY corrected 3' end positions",
            config=config
        )
        logger.info(f"Provenance tracking initialized for {output_dir}")

    # Process BAM file
    import time as _time
    _t_correct_total = _time.perf_counter()
    try:
        logger.info("Processing BAM file...")

        # Spike-in filtering (pre-processing step)
        bam_to_process = str(config['bam_path'])
        spikein_stats = {}
        if config['filter_spikein']:
            filtered_bam = str(config['output_path']).replace('.tsv', '_spikein_filtered.bam')
            spikein_report = str(config['output_path']).replace('.tsv', '_spikein_report.txt')
            logger.info(f"Filtering spike-in reads ({', '.join(config['filter_spikein'])})...")
            _t_spikein = _time.perf_counter()
            spikein_stats = filter_spikein_reads(
                input_bam=bam_to_process,
                output_bam=filtered_bam,
                known_genes=config['filter_spikein'],
                report_path=spikein_report,
            )
            bam_to_process = filtered_bam
            logger.info(f"  Spike-in reads removed: {spikein_stats.get('spikein_reads', 0):,}")
            logger.info(f"[TIMING] Spike-in filter: {_time.perf_counter() - _t_spikein:.1f}s")

        # Load annotated junctions for Module 2F (3'SS truncation rescue).
        # Without these, only reads whose own CIGAR contains an N operation near
        # the alignment's 5' end can trigger rescue — truncated reads that end
        # exactly at the 3'SS boundary (no N in CIGAR, just a soft-clip) are missed.
        _t_junc = _time.perf_counter()
        annotated_junctions = None
        if config.get('annotation_path'):
            from .consensus import load_annotated_junctions
            annotated_junctions = load_annotated_junctions(str(config['annotation_path']))
            logger.info(
                f"Loaded {len(annotated_junctions):,} annotated junctions for 3'SS rescue "
                f"({_time.perf_counter() - _t_junc:.1f}s)"
            )

        # Build per-read gene attribution interval trees (optional, requires annotation)
        gene_interval_trees = None
        if config.get('annotation_path'):
            try:
                from .analyze.gene_attribution import build_cds_interval_tree
                from .analyze_command import load_annotation as _load_annotation_for_trees
                _ann_df = _load_annotation_for_trees(str(config['annotation_path']),
                                                      normalize_chroms=False)
                gene_interval_trees = build_cds_interval_tree(_ann_df)
                logger.info(
                    "Built gene interval trees for per-read attribution (%d (chrom, strand) pairs)",
                    len(gene_interval_trees),
                )
            except Exception as e:
                logger.warning("Per-read gene attribution unavailable: %s", e)

        # Choose processing mode
        _t_proc = _time.perf_counter()
        if streaming_mode:
            # Streaming mode - memory efficient for large BAMs
            chunk_size = getattr(args, 'chunk_size', 10000)
            stats = bam_processor.process_bam_streaming(
                bam_path=bam_to_process,
                genome_path=str(config['genome_path']),
                output_path=str(config['output_path']),
                chunk_size=chunk_size,
                apply_atract=config['apply_atract'],
                apply_ag_mispriming=config['apply_ag_mispriming'],
                ag_threshold=config['ag_threshold'],
                apply_polya_trim=config['apply_polya_trim'],
                apply_indel_correction=config['apply_indel_correction'],
                netseq_dir=str(config['netseq_dir']) if config['netseq_dir'] else None,
                annotated_junctions=annotated_junctions,
                gene_interval_trees=gene_interval_trees,
            )
            report = generate_stats_report(stats)
        else:
            # Standard parallel processing
            # Determine variant output path
            variant_output_path = None
            if config['variant_aware'] and config['output_path']:
                variant_output_path = str(config['output_path']).replace('.tsv', '_potential_variants.tsv')

            results, stats = bam_processor.process_bam_file_parallel(
                bam_path=bam_to_process,
                genome_path=str(config['genome_path']),
                n_threads=n_threads,
                apply_atract=config['apply_atract'],
                apply_ag_mispriming=config['apply_ag_mispriming'],
                ag_threshold=config['ag_threshold'],
                apply_polya_trim=config['apply_polya_trim'],
                apply_indel_correction=config['apply_indel_correction'],
                netseq_dir=str(config['netseq_dir']) if config['netseq_dir'] else None,
                output_path=str(config['output_path']) if config['output_path'] else None,
                show_progress=not args.verbose,  # Use verbose logging instead of progress bar
                return_stats=True,
                variant_aware=config['variant_aware'],
                variant_output_path=variant_output_path,
                annotated_junctions=annotated_junctions,
                gene_interval_trees=gene_interval_trees,
            )
            report = generate_stats_report(stats)

        logger.info(f"[TIMING] BAM processing: {_time.perf_counter() - _t_proc:.1f}s")

        # Propagate spike-in count into stats
        if spikein_stats:
            stats.spikein_reads_filtered = spikein_stats.get('spikein_reads', 0)

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

        # Write poly(A)-trimmed BAM if requested
        if config.get('output_bam'):
            _t_bam = _time.perf_counter()
            logger.info(f"Writing poly(A)-trimmed BAM to {config['output_bam']}...")
            _unsorted_polya = str(config['output_bam']) + '.unsorted.bam'
            bam_stats = bam_processor.write_polya_trimmed_bam(
                str(bam_to_process),
                _unsorted_polya,
            )
            logger.info(
                f"  Reads written: {bam_stats['total']:,}  "
                f"trimmed: {bam_stats['trimmed']:,}  "
                f"bases removed: {bam_stats['bases_trimmed']:,}"
            )
            pysam.sort('-o', str(config['output_bam']), _unsorted_polya)
            pysam.index(str(config['output_bam']))
            import os as _os_polya
            _os_polya.unlink(_unsorted_polya)
            logger.info(f"  Sorted and indexed {config['output_bam']}")
            logger.info(f"[TIMING] Poly(A) BAM trim: {_time.perf_counter() - _t_bam:.1f}s")

        # Write fully-corrected BAM if requested (hard-clips aligned poly-A regions)
        if config.get('corrected_bam') and config.get('output_path'):
            import subprocess as _subprocess
            _t_cbam = _time.perf_counter()
            _unsorted_bam = str(config['corrected_bam']) + '.unsorted.bam'
            logger.info(f"Writing corrected BAM to {config['corrected_bam']}...")
            cbam_stats = bam_processor.write_corrected_bam(
                str(bam_to_process),
                str(config['output_path']),
                _unsorted_bam,
            )
            logger.info(
                f"  Reads written: {cbam_stats['total']:,}  "
                f"clipped: {cbam_stats['clipped']:,}  "
                f"unchanged: {cbam_stats['unchanged']:,}"
            )
            # Sort and index so the BAM is immediately usable in IGV
            pysam.sort('-o', str(config['corrected_bam']), _unsorted_bam)
            pysam.index(str(config['corrected_bam']))
            import os as _os
            _os.unlink(_unsorted_bam)
            logger.info(
                f"  Sorted and indexed {config['corrected_bam']}"
            )
            logger.info(f"[TIMING] Corrected BAM write: {_time.perf_counter() - _t_cbam:.1f}s")

        # Write soft-clipped BAM if requested (poly-A bases retained but soft-clipped)
        if config.get('softclipped_bam') and config.get('output_path'):
            _t_sbam = _time.perf_counter()
            _unsorted_sbam = str(config['softclipped_bam']) + '.unsorted.bam'
            logger.info(f"Writing soft-clipped BAM to {config['softclipped_bam']}...")
            sbam_stats = bam_processor.write_softclipped_bam(
                str(bam_to_process),
                str(config['output_path']),
                _unsorted_sbam,
            )
            logger.info(
                f"  Reads written: {sbam_stats['total']:,}  "
                f"clipped: {sbam_stats['clipped']:,}  "
                f"unchanged: {sbam_stats['unchanged']:,}"
            )
            pysam.sort('-o', str(config['softclipped_bam']), _unsorted_sbam)
            pysam.index(str(config['softclipped_bam']))
            import os as _os
            _os.unlink(_unsorted_sbam)
            logger.info(f"  Sorted and indexed {config['softclipped_bam']}")
            logger.info(f"[TIMING] Soft-clipped BAM write: {_time.perf_counter() - _t_sbam:.1f}s")

        # Write NET-seq bedgraph files if requested
        if config.get('bedgraph_prefix') and config.get('output_path'):
            import pandas as _pd
            from .netseq_output import write_bedgraph as _write_bedgraph
            _t_bg = _time.perf_counter()
            _tsv_path = str(config['output_path'])
            _prefix = str(config['bedgraph_prefix'])
            logger.info(f"Writing NET-seq bedgraph files (prefix: {_prefix})...")
            try:
                _df = _pd.read_csv(
                    _tsv_path, sep='\t',
                    usecols=lambda c: c in ('chrom', 'strand', 'corrected_3prime', 'fraction'),
                )
                if 'fraction' not in _df.columns:
                    logger.warning(
                        "--write-bedgraph: 'fraction' column absent from TSV; "
                        "bedgraph output skipped (NET-seq input required)."
                    )
                else:
                    _counts_series = (
                        _df.groupby(['chrom', 'strand', 'corrected_3prime'])['fraction'].sum()
                    )
                    _counts = {
                        (chrom, strand, int(pos)): float(val)
                        for (chrom, strand, pos), val in _counts_series.items()
                    }
                    for _strand, _name in [('+', 'plus'), ('-', 'minus')]:
                        _bg_path = Path(f"{_prefix}.{_name}.bedgraph")
                        _write_bedgraph(
                            _counts, _bg_path, _strand,
                            total_reads=1,
                            normalize_rpm=False,
                            track_name=f"{_bg_path.stem}",
                        )
                        _n = sum(1 for (_, s, _) in _counts if s == _strand)
                        logger.info(f"  Wrote {_bg_path} ({_n:,} positions)")
            except Exception as _bg_exc:
                logger.warning(f"Failed to write bedgraph: {_bg_exc}")
            logger.info(f"[TIMING] Bedgraph write: {_time.perf_counter() - _t_bg:.1f}s")

        # Save report if requested
        if args.report:
            logger.info(f"Saving detailed report to {args.report}...")
            with open(args.report, 'w') as f:
                f.write(report)

        # Save provenance
        if provenance:
            # Record output files
            output_path = Path(config['output_path'])
            if output_path.exists():
                provenance.add_output_file(
                    output_path,
                    source_files=[config['bam_path']],
                    metadata={'stats': dataclasses.asdict(stats)}  # serialize dataclass
                )
            stats_path = Path(str(config['output_path']).replace('.tsv', '_stats.tsv'))
            if stats_path.exists():
                provenance.add_output_file(stats_path)
            provenance.save()
            logger.info(f"Provenance saved to {provenance.output_dir}")

        logger.info("")
        logger.info("=" * 70)
        logger.info("RECTIFY completed successfully!")
        logger.info(f"[TIMING] Correction total: {_time.perf_counter() - _t_correct_total:.1f}s")
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
    parser.add_argument('--skip-variant-aware', action='store_true',
                        help='Skip variant-aware homopolymer rescue (enabled by default)')
    parser.add_argument('--netseq-dir', type=Path)
    parser.add_argument('--netseq-samples', nargs='+')
    parser.add_argument('--polya-model', type=Path)
    parser.add_argument('--report', type=Path)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()
    run(args)
