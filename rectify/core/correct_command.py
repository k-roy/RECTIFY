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


def _strip_aligner_prefix(bam_entry: str) -> str:
    """Strip optional ``aligner:`` prefix from a BAM path entry.

    ``--aligner-bams`` accepts both plain paths and ``aligner:path`` pairs
    (the same format used by ``rectify consensus``).  When the prefix is
    present (e.g. ``"minimap2:/path/to/file.bam"``), the aligner-name
    portion is not needed for junction pool construction and must be
    stripped so pysam receives a valid file path.
    """
    if ':' in bam_entry:
        prefix, _, rest = bam_entry.partition(':')
        if '/' not in prefix:   # aligner names never contain path separators
            return rest
    return bam_entry


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

    # Resolve NET-seq directory.
    # NET-seq refinement in 'rectify correct' is DISABLED by default.
    # It is only enabled when the user explicitly passes --netseq-dir.
    #
    # Rationale: auto-loading bundled NET-seq during BAM correction
    # caused widespread CPA mis-assignment.  The bundled yeast signal
    # is noisy in many loci (Pol II contamination, multi-mapper signal)
    # and pushes CPA positions far from where the DRS coverage boundary
    # indicates they should be.  NET-seq-based refinement is reserved for
    # the 3' end analysis step where signal quality can be assessed
    # independently.
    resolved_netseq_dir = None
    if getattr(args, 'netseq_dir', None):
        # Custom dir provided — validate it exists
        if not args.netseq_dir.exists():
            errors.append(f"NET-seq directory not found: {args.netseq_dir}")
        elif not args.netseq_dir.is_dir():
            errors.append(f"NET-seq path is not a directory: {args.netseq_dir}")
        else:
            resolved_netseq_dir = args.netseq_dir
            logging.info(f"NET-seq refinement ENABLED (explicit --netseq-dir)")

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

    # dT-primed cDNA flag enables AG mispriming detection (a cDNA-synthesis artifact).
    # Poly(A) trimming and indel correction are always enabled — the poly-A tail is
    # present in both DRS and dT-primed cDNA reads and always requires correction.
    is_dt_primed = getattr(args, 'dT_primed_cDNA', False) or getattr(args, 'polya_sequenced', False)

    config = {
        'bam_path': bam_path,
        'genome_path': genome_path,
        'annotation_path': annotation_path,
        'output_path': args.output,
        'apply_atract': not args.skip_atract_check,
        # AG mispriming: reverse-transcriptase slippage at AG runs — cDNA-specific artifact.
        # Disabled for DRS (default); enabled only when --dT-primed-cDNA is set.
        'apply_ag_mispriming': is_dt_primed and not args.skip_ag_check,
        'ag_threshold': getattr(args, 'ag_threshold', 0.65),
        # Poly(A) trimming and indel correction: always enabled — poly-A is present
        # in the read sequence for both DRS and dT-primed cDNA protocols.
        'apply_polya_trim': not args.skip_polya_trim,
        'apply_indel_correction': not args.skip_indel_correction,
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
        # Junction refinement (Module 2H)
        'aligner_bams':              [_strip_aligner_prefix(str(p)) for p in getattr(args, 'aligner_bams', []) or []],
        'junction_hp_pen':           getattr(args, 'junction_hp_pen', 0.25),
        'junction_search_radius':    getattr(args, 'junction_search_radius', 5000),
        'junction_window':           getattr(args, 'junction_window', 15),
        'junction_max_slide':        getattr(args, 'junction_max_slide', 10),
        'junction_max_boundary_shift': getattr(args, 'junction_max_boundary_shift', 50),
        'junction_penalty_table':    getattr(args, 'junction_penalty_table', None),
        'str_penalty_table':         getattr(args, 'str_penalty_table', None),
        # Pre-computed cache paths from `rectify prescan` (chunked correction pipelines).
        # When set, skip the corresponding first-pass computation and load from disk.
        'variant_scan_cache':        getattr(args, 'variant_scan_cache', None),
        'junction_pool_cache':       getattr(args, 'junction_pool_cache', None),
    }

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
    if config.get('variant_scan_cache'):
        logger.info(f"    Variant scan cache:  {config['variant_scan_cache']}")
    logger.info(f"  NET-seq refinement:    {'ENABLED' if config['netseq_dir'] else 'DISABLED'}")
    logger.info(f"  Spike-in filter:       {'ENABLED (' + ','.join(config['filter_spikein']) + ')' if config['filter_spikein'] else 'DISABLED'}")
    _n_abams = len(config.get('aligner_bams', []))
    _jpool_cache = config.get('junction_pool_cache')
    if _jpool_cache:
        logger.info(f"  Junction refinement:   ENABLED (pre-built pool cache: {_jpool_cache})")
    else:
        logger.info(f"  Junction refinement:   {'ENABLED (' + str(_n_abams) + ' aligner BAMs)' if _n_abams else 'DISABLED (pass --aligner-bams to enable)'}")
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
            # Derive output names from the input BAM stem, not the TSV path.
            # e.g. wt_by4742_rep1.consensus.sorted.bam → wt_by4742_rep1.spikein_filtered.bam
            _bam_stem = Path(str(config['bam_path'])).stem
            for _sfx in ('.consensus.sorted', '.consensus', '.rectified', '.sorted'):
                if _bam_stem.endswith(_sfx):
                    _bam_stem = _bam_stem[:-len(_sfx)]
                    break
            _out_dir = Path(str(config['output_path'])).parent
            filtered_bam = str(_out_dir / f"{_bam_stem}.spikein_filtered.bam")
            spikein_report = str(_out_dir / f"{_bam_stem}.spikein_report.txt")
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

        # Module 2H: Junction N-op boundary refinement (optional pre-processing step).
        # When --aligner-bams are provided (or a --junction-pool-cache pkl), replace
        # imprecise N-op boundaries in the consensus BAM with the best-supported
        # candidate junction using split-alignment scoring with homopolymer-aware
        # edit distance.
        _junction_pool_cache = config.get('junction_pool_cache')
        _run_2h = (config.get('aligner_bams') or _junction_pool_cache) and config.get('annotation_path')
        if _run_2h:
            _t_refine = _time.perf_counter()
            logger.info("Module 2H: Junction N-op boundary refinement...")
            try:
                from .junction_refiner import refine_bam_junctions
                from .consensus import load_annotated_junctions as _load_annot_j

                _annot_j = _load_annot_j(str(config['annotation_path']))
                _bam_stem = Path(bam_to_process).stem
                for _sfx in ('.consensus.sorted', '.consensus', '.rectified', '.sorted'):
                    if _bam_stem.endswith(_sfx):
                        _bam_stem = _bam_stem[:-len(_sfx)]
                        break
                _out_dir = Path(str(config['output_path'])).parent
                # Include a short unique ID so concurrent runs on the same sample
                # don't collide on the intermediate file.
                import uuid as _uuid
                _run_id = _uuid.uuid4().hex[:8]
                _refined_bam = str(_out_dir / f"{_bam_stem}.junction_refined_{_run_id}.bam")

                from ..utils.genome import load_genome as _load_genome_for_refine
                _refine_genome = _load_genome_for_refine(str(config['genome_path']))

                # Load pre-built junction pool from cache if available.
                _prebuilt_pool = None
                _prebuilt_annot_set = None
                if _junction_pool_cache and Path(_junction_pool_cache).exists():
                    import pickle as _pkl
                    logger.info(f"  Loading pre-built junction pool from cache: {_junction_pool_cache}")
                    with open(_junction_pool_cache, 'rb') as _pf:
                        _pool_data = _pkl.load(_pf)
                    _prebuilt_pool = _pool_data['all_junctions']
                    _prebuilt_annot_set = _pool_data['annotated_set']
                    logger.info(
                        "  Pre-built pool: %d junctions (%d annotated)",
                        len(_prebuilt_pool), len(_prebuilt_annot_set),
                    )

                _refine_stats = refine_bam_junctions(
                    input_bam=bam_to_process,
                    output_bam=_refined_bam,
                    aligner_bams=config['aligner_bams'],
                    annotated_junctions=_annot_j,
                    genome=_refine_genome,
                    hp_pen=config['junction_hp_pen'],
                    W=config['junction_window'],
                    max_slide=config['junction_max_slide'],
                    search_radius=config['junction_search_radius'],
                    max_boundary_shift=config['junction_max_boundary_shift'],
                    boundary_error_window=config.get('junction_boundary_error_window', 10),
                    sort_and_index=True,
                    penalty_table_path=config.get('junction_penalty_table'),
                    str_penalty_table_path=config.get('str_penalty_table'),
                    prebuilt_junction_pool=_prebuilt_pool,
                    prebuilt_annotated_set=_prebuilt_annot_set,
                )
                bam_to_process = _refined_bam
                logger.info(
                    "  Junction refinement: %d reads with N-ops, %d refined, %d unchanged",
                    _refine_stats['n_op_reads'], _refine_stats['refined'],
                    _refine_stats['unchanged'],
                )
                logger.info(f"[TIMING] Junction refinement: {_time.perf_counter() - _t_refine:.1f}s")
            except Exception as _exc:
                logger.warning(
                    "Module 2H junction refinement failed (non-fatal, continuing): %s", _exc
                )
            finally:
                # Explicitly release the genome dict to free ~500 MB RAM
                # before the main correction loads its own copy.
                try:
                    del _refine_genome
                except NameError:
                    pass
                import gc as _gc
                _gc.collect()

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
        _polya_model_path = str(config['polya_model_path']) if config.get('polya_model_path') else None

        if streaming_mode:
            if n_threads > 1:
                # Parallel streaming: region workers + stream output to disk
                # Best for large BAMs (> ~5M reads) — combines parallel throughput
                # with constant memory usage.
                variant_output_path = None
                if config['variant_aware'] and config['output_path']:
                    variant_output_path = str(config['output_path']).replace('.tsv', '_potential_variants.tsv')
                stats = bam_processor.process_bam_streaming_parallel(
                    bam_path=bam_to_process,
                    genome_path=str(config['genome_path']),
                    output_path=str(config['output_path']),
                    n_threads=n_threads,
                    apply_atract=config['apply_atract'],
                    apply_ag_mispriming=config['apply_ag_mispriming'],
                    ag_threshold=config['ag_threshold'],
                    apply_polya_trim=config['apply_polya_trim'],
                    apply_indel_correction=config['apply_indel_correction'],
                    netseq_dir=str(config['netseq_dir']) if config['netseq_dir'] else None,
                    variant_aware=config['variant_aware'],
                    variant_output_path=variant_output_path,
                    annotated_junctions=annotated_junctions,
                    gene_interval_trees=gene_interval_trees,
                    polya_model_path=_polya_model_path,
                    checkpoint_dir=getattr(args, 'checkpoint_dir', None),
                    variant_scan_cache=config.get('variant_scan_cache'),
                )
            else:
                # Single-threaded streaming for single-core or debugging
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
                    polya_model_path=_polya_model_path,
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
                polya_model_path=_polya_model_path,
                variant_scan_cache=config.get('variant_scan_cache'),
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

        # Write corrected BAMs (hard-clip and/or soft-clip).
        # When both are requested, use write_dual_bam for a single-pass read of the input BAM.
        _want_hc = bool(config.get('corrected_bam') and config.get('output_path'))
        _want_sc = bool(config.get('softclipped_bam') and config.get('output_path'))

        # Load genome once for homopolymer CIGAR surgery in BAM writers.
        _genome_for_bam = None
        if config.get('genome_path'):
            try:
                _genome_for_bam = genome_utils.load_genome(str(config['genome_path']))
            except Exception as _e:
                logger.warning(f"Could not load genome for homopolymer BAM surgery: {_e}")

        if _want_hc and _want_sc:
            _t_cbam = _time.perf_counter()
            _unsorted_bam  = str(config['corrected_bam'])  + '.unsorted.bam'
            _unsorted_sbam = str(config['softclipped_bam']) + '.unsorted.bam'
            logger.info(
                f"Writing hardclip + softclip BAMs in single pass "
                f"({config['corrected_bam']}, {config['softclipped_bam']})..."
            )
            cbam_stats, sbam_stats = bam_processor.write_dual_bam(
                str(bam_to_process),
                str(config['output_path']),
                _unsorted_bam,
                _unsorted_sbam,
                genome=_genome_for_bam,
            )
            logger.info(
                f"  Hardclip: total={cbam_stats['total']:,}  "
                f"clipped={cbam_stats['clipped']:,}  unchanged={cbam_stats['unchanged']:,}"
            )
            logger.info(
                f"  Softclip: total={sbam_stats['total']:,}  "
                f"clipped={sbam_stats['clipped']:,}  unchanged={sbam_stats['unchanged']:,}"
            )
            import os as _os
            for _unsorted, _final in [
                (_unsorted_bam,  config['corrected_bam']),
                (_unsorted_sbam, config['softclipped_bam']),
            ]:
                pysam.sort('-o', str(_final), _unsorted)
                pysam.index(str(_final))
                _os.unlink(_unsorted)
                logger.info(f"  Sorted and indexed {_final}")
            logger.info(f"[TIMING] Dual BAM write: {_time.perf_counter() - _t_cbam:.1f}s")

        elif _want_hc:
            import subprocess as _subprocess
            _t_cbam = _time.perf_counter()
            _unsorted_bam = str(config['corrected_bam']) + '.unsorted.bam'
            logger.info(f"Writing corrected BAM to {config['corrected_bam']}...")
            cbam_stats = bam_processor.write_corrected_bam(
                str(bam_to_process),
                str(config['output_path']),
                _unsorted_bam,
                genome=_genome_for_bam,
            )
            logger.info(
                f"  Reads written: {cbam_stats['total']:,}  "
                f"clipped: {cbam_stats['clipped']:,}  "
                f"unchanged: {cbam_stats['unchanged']:,}"
            )
            pysam.sort('-o', str(config['corrected_bam']), _unsorted_bam)
            pysam.index(str(config['corrected_bam']))
            import os as _os
            _os.unlink(_unsorted_bam)
            logger.info(f"  Sorted and indexed {config['corrected_bam']}")
            logger.info(f"[TIMING] Corrected BAM write: {_time.perf_counter() - _t_cbam:.1f}s")

        elif _want_sc:
            _t_sbam = _time.perf_counter()
            _unsorted_sbam = str(config['softclipped_bam']) + '.unsorted.bam'
            logger.info(f"Writing soft-clipped BAM to {config['softclipped_bam']}...")
            sbam_stats = bam_processor.write_softclipped_bam(
                str(bam_to_process),
                str(config['output_path']),
                _unsorted_sbam,
                genome=_genome_for_bam,
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

        # Write NET-seq-assigned bedgraph for Cat6 reads (rectified_ambiguous_pA_netseq_assigned)
        # Automatically produced whenever a corrected or soft-clipped BAM is written.
        if (config.get('corrected_bam') or config.get('softclipped_bam')) and config.get('output_path'):
            _t_nbg = _time.perf_counter()
            # Derive prefix from the corrected BAM path (or softclipped if no corrected)
            _bam_for_prefix = config.get('corrected_bam') or config.get('softclipped_bam')
            _bam_stem = Path(str(_bam_for_prefix)).stem  # e.g. "rectified"
            _bam_dir  = Path(str(_bam_for_prefix)).parent
            _nbg_prefix = str(_bam_dir / f"{_bam_stem}_ambiguous_pA_netseq_assigned")
            logger.info(f"Writing NET-seq assigned bedgraph (prefix: {_nbg_prefix})...")
            try:
                _nbg_counts = bam_processor.write_netseq_assigned_bedgraph(
                    str(config['output_path']),
                    _nbg_prefix,
                )
                logger.info(
                    f"  NET-seq bedgraph: plus={_nbg_counts.get('plus', 0):,} pos, "
                    f"minus={_nbg_counts.get('minus', 0):,} pos"
                )
            except Exception as _nbg_exc:
                logger.warning(f"Failed to write NET-seq assigned bedgraph: {_nbg_exc}")
            logger.info(f"[TIMING] NET-seq assigned bedgraph: {_time.perf_counter() - _t_nbg:.1f}s")

            # Write general corrected-3'-end bedgraph for ALL reads (Cat1–6).
            # Uses the 'fraction' column so Cat6 multi-peak rows contribute fractional
            # values; all other reads contribute 1.0.
            _t_c3bg = _time.perf_counter()
            # Use TSV stem ("corrected_3ends") so the name doesn't depend on
            # whatever the BAM happens to be called (avoids redundant names
            # like "rectified_corrected_3end_corrected_3ends" when the BAM
            # already contains "_corrected_3end" in its stem).
            _tsv_stem = Path(str(config['output_path'])).stem
            _c3bg_prefix = str(_bam_dir / _tsv_stem)
            logger.info(f"Writing corrected-3'-end bedgraph (prefix: {_c3bg_prefix})...")
            try:
                _c3bg_counts = bam_processor.write_corrected_3ends_bedgraph(
                    str(config['output_path']),
                    _c3bg_prefix,
                )
                logger.info(
                    f"  Corrected-3'-end bedgraph: plus={_c3bg_counts.get('plus', 0):,} pos, "
                    f"minus={_c3bg_counts.get('minus', 0):,} pos"
                )
            except Exception as _c3bg_exc:
                logger.warning(f"Failed to write corrected-3'-end bedgraph: {_c3bg_exc}")
            logger.info(f"[TIMING] Corrected-3'-end bedgraph: {_time.perf_counter() - _t_c3bg:.1f}s")

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

        # Clean up the intermediate junction_refined BAM (temporary file).
        # It's no longer needed once the main correction has completed.
        if config.get('aligner_bams'):
            try:
                import os as _os
                for _ext in ('', '.bai'):
                    _tmp = bam_to_process + _ext
                    if _os.path.exists(_tmp) and 'junction_refined_' in _tmp:
                        _os.remove(_tmp)
            except Exception:
                pass  # non-fatal if cleanup fails

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
    parser.add_argument('--dT-primed-cDNA', dest='dT_primed_cDNA', action='store_true')
    parser.add_argument('--polya-sequenced', dest='dT_primed_cDNA', action='store_true')  # deprecated
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
