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
from ...slurm import set_thread_limits, get_available_cpus, get_slurm_info

from ..bam import bam_processor
from ..bam import parallel as bam_parallel
from ..bam.processing_stats import write_stats_tsv, generate_stats_report
from ..spikein_filter import filter_spikein_reads
from ...utils import genome as genome_utils
from ...utils.provenance import init_provenance


def _derive_sample_id(args) -> str:
    """Derive a stable sample ID from args for sidecar naming.

    Priority:
    1. ``args.corrected_bam`` stem (most specific — the final corrected BAM path)
    2. ``args.input`` stem (input BAM / FASTQ stem)
    3. Fallback: ``"rectify_sample"``
    """
    if getattr(args, 'corrected_bam', None):
        return Path(str(args.corrected_bam)).stem
    if getattr(args, 'input', None):
        stem = Path(str(args.input)).stem
        # Strip common suffixes (e.g., .sorted, .consensus) for a cleaner ID.
        for sfx in ('.sorted', '.consensus', '.rectified', '.bam'):
            if stem.endswith(sfx):
                stem = stem[:-len(sfx)]
        return stem or "rectify_sample"
    return "rectify_sample"


def _estimate_n_reads(bam_path: str) -> int:
    """Estimate total read count from BAM index statistics.

    Uses ``pysam.idxstats`` which requires the BAM to be indexed. If the
    index is absent or the call fails, returns 0 (falls through to sequential
    path as a safe fallback).
    """
    try:
        lines = pysam.idxstats(str(bam_path)).strip().split('\n')
        total = 0
        for line in lines:
            parts = line.split('\t')
            if len(parts) >= 3:
                total += int(parts[2])  # mapped reads
        return total
    except Exception:
        return 0


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

    Handles preprocessing of FASTQ files and bundled genomes, resolves
    protocol flags, and assembles the correction config dict consumed by
    ``bam_parallel.process_bam_file_parallel`` /
    ``bam_parallel.process_bam_streaming_parallel``.

    Protocol flags:
        ``--dT-primed-cDNA`` (or deprecated alias ``--polya-sequenced``):
            Enables AG-mispriming detection (Module AG).  Poly(A) trimming
            and indel correction are always enabled regardless of this flag.
        ``--short-read``:
            Disables poly(A)-tail modules that require a sequenced poly(A) tail
            (``apply_polya_trim``, ``apply_indel_correction``).
            A-tract walk-back (``apply_atract``) is NOT disabled — CSP internal
            priming causes the poly-A complement to align into downstream genomic
            A-runs, shifting the 3' end rightward; the genomic walk-back corrects
            this regardless of read length.
            AG mispriming (``apply_ag_mispriming``) is also NOT disabled when
            ``--dT-primed-cDNA`` is set — QuantSeq is the primary use case.
            Use for Illumina/Aviti ≤150 bp data.

    Returns:
        Dict with the following keys:

        Paths:
            ``bam_path``, ``genome_path``, ``annotation_path``,
            ``output_path``, ``output_bam``, ``corrected_bam``,
            ``softclipped_bam``, ``bedgraph_prefix``,
            ``polya_model_path``, ``netseq_dir``,
            ``junction_penalty_table``, ``str_penalty_table``,
            ``variant_scan_cache``, ``junction_pool_cache``

        Module flags (bool):
            ``apply_atract``, ``apply_ag_mispriming``,
            ``apply_polya_trim``, ``apply_indel_correction``,
            ``variant_aware``, ``filter_spikein``

        Parameters:
            ``ag_threshold``, ``netseq_samples``, ``threads``,
            ``verbose``, ``aligner_bams``,
            ``junction_hp_pen``, ``junction_search_radius``,
            ``junction_window``, ``junction_max_slide``,
            ``junction_max_boundary_shift``
    """
    from ..align.preprocess import detect_input_type, prepare_input, prepare_bundled_genome
    from ...data import (
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

    # Guard: --dT-primed-cDNA (QuantSeq REV) with FASTQ input is unsafe
    # because the built-in FASTQ aligner (preprocess.prepare_input) uses
    # minimap2 only, which is the wrong panel for short-read antisense data.
    # Force the user to align separately with `rectify align --short-read`.
    is_dt_primed_arg = (getattr(args, 'dT_primed_cDNA', False)
                        or getattr(args, 'polya_sequenced', False))
    if input_type in ('fastq', 'fastq.gz') and is_dt_primed_arg:
        errors.append(
            "FASTQ input is not supported with --dT-primed-cDNA: the built-in "
            "aligner uses minimap2 only, which mis-aligns short antisense "
            "reads. Pre-align with `rectify align --short-read INPUT.fastq "
            "--genome GENOME -o ALIGNED.bam` first, then pass ALIGNED.bam to "
            "`rectify correct --dT-primed-cDNA --short-read`."
        )

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
    # ONT PCR-cDNA (e.g. SQK-PCB114): same strand convention as DRS (no flip);
    # poly-A tail IS present as a 3' soft-clip from minimap2 alignment.
    # AG mispriming disabled (no oligo-dT priming step).
    is_ont_cdna = getattr(args, 'ONT_cDNA', False)
    # Short-read mode (Illumina/Aviti): no poly(A) tail in reads; disable all
    # modules that assume a poly(A) soft-clip or A-tract walk-back.
    is_short_read = getattr(args, 'short_read', False)

    config = {
        'bam_path': bam_path,
        'genome_path': genome_path,
        'annotation_path': annotation_path,
        'output_path': args.output,
        # A-tract: detects reads whose 3' end has slid into a genomic A/T-run, and walks
        # the position back to the first non-A (non-T for minus strand) base upstream.
        # Uses genomic reference sequence — does NOT require a sequenced poly(A) tail.
        # Enabled for short reads: CSP internal priming causes poly-A complement T's to
        # reverse-complement into A's that align into downstream genomic A-runs, shifting
        # the reported 3' end rightward.  The same genomic walk-back corrects this.
        'apply_atract': not args.skip_atract_check,
        # AG mispriming: oligo-dT mispriming onto genomic A/G-rich regions creates false 3' ends.
        # Disabled for DRS (default) — direct RNA sequencing has no priming step, so no
        # mispriming artifact is possible.
        # Enabled when --dT-primed-cDNA is set (QuantSeq / oligo-dT cDNA), regardless of
        # --short-read: the module checks DOWNSTREAM GENOMIC SEQUENCE at the reported 3' end
        # via get_downstream_sequence() — it does not inspect read sequence — so it works
        # correctly for both short and long reads.  In fact, short-read (QuantSeq) data is
        # the primary use case: the reads stop at the cleavage site and the AG-rich context
        # is only detectable from the genome, not from the reads.
        'apply_ag_mispriming': is_dt_primed and not args.skip_ag_check,
        'ag_threshold': getattr(args, 'ag_threshold', 17.0),
        # Poly(A) trimming and indel correction: always enabled for long reads — poly-A
        # is present in the read sequence for both DRS and dT-primed cDNA protocols.
        # Disabled for short reads (no poly-A tail).
        'apply_polya_trim': not args.skip_polya_trim and not is_short_read,
        'apply_indel_correction': not args.skip_indel_correction and not is_short_read,
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
        # Antisense protocol flag: when True, gene strand = opposite of read strand.
        # Passed through to correct_read_3prime to flip 3' end assignment.
        'dt_primed_cDNA': is_dt_primed,
        # ONT PCR-cDNA flag: same strand convention as DRS (no strand flip).
        # Poly-A IS in the read as a right soft-clip; AG mispriming is disabled.
        # Used for stats protocol label only — no bam_processor behaviour changes needed.
        'ont_cDNA': is_ont_cdna,
        # Organism name (normalized, e.g. 'saccharomyces_cerevisiae') for bundled
        # data resolution. May be None when only a genome path was supplied.
        'organism': normalize_organism(organism) if organism else None,
        # Minimum MAPQ filter: reads with mapping_quality < min_mapq are skipped.
        # Default 0 = no filter. Use 1 to exclude MAPQ=0 multi-mapping reads.
        'min_mapq': getattr(args, 'min_mapq', 0),
        # Minimum aligned length filter: reads where reference_length < min_aligned_length
        # are skipped. Removes reads with very short alignments to low-complexity sequence
        # (e.g. T-tract alignments with only 10-20 bp of unique context).
        # Default 0 = no filter. Use 30 for QuantSeq REV short-read data.
        'min_aligned_length': getattr(args, 'min_aligned_length', 0),
    }

    return config


def run(args):
    """
    Execute the 'correct' command.

    Args:
        args: Parsed command-line arguments from argparse
    """
    import os as _os_run
    import sys as _sys_run
    import time as _time_run
    from datetime import datetime, timezone

    # Capture stage start time for provenance sidecar (must be before any work).
    _stage_started_at = datetime.now(timezone.utc).isoformat()

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
    is_dt_primed = getattr(args, 'dT_primed_cDNA', False) or getattr(args, 'polya_sequenced', False)
    is_ont_cdna = getattr(args, 'ONT_cDNA', False)
    is_short_read = getattr(args, 'short_read', False)

    # Log configuration
    logger.info("Configuration:")
    logger.info(f"  Input BAM:             {config['bam_path']}")
    logger.info(f"  Reference genome:      {config['genome_path']}")
    logger.info(f"  Output TSV:            {config['output_path']}")
    logger.info(f"  Threads:               {n_threads}")
    streaming_mode = getattr(args, 'streaming', False)
    logger.info(f"  Streaming mode:        {'ENABLED' if streaming_mode else 'DISABLED'}")
    if is_short_read:
        logger.info(f"  Protocol:              short-read (poly(A) modules disabled)")
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

    # ---------- Resume check (Commit B, 2026-05-20) ----------
    _sample_id = _derive_sample_id(args)
    _output_dir = Path(config['output_path']).parent if config.get('output_path') else None
    _emit_merged_tsv = getattr(args, 'emit_merged_tsv', False)
    _legacy_single_threaded = getattr(args, 'legacy_single_threaded', False)
    _tmp_dir = getattr(args, 'tmp_dir', None)

    if _output_dir is not None:
        from rectify.core.provenance import can_skip_stage
        from rectify.utils.version import get_rectify_git_sha as _get_sha
        _rectify_sha = _get_sha()
        _current_inputs = {
            'bam': str(config['bam_path']),
            'genome': str(config['genome_path']) if config.get('genome_path') else None,
        }
        if config.get('annotation_path'):
            _current_inputs['annotation'] = str(config['annotation_path'])
        # Remove None values — can_skip_stage expects existing paths.
        _current_inputs = {k: v for k, v in _current_inputs.items() if v}

        _skip_decision = can_skip_stage(
            stage='correct',
            sample_output=_output_dir,
            sample_id=_sample_id,
            current_inputs=_current_inputs,
            current_argv=_sys_run.argv,
            rectify_git_sha=_rectify_sha,
            force=getattr(args, 'force_all', False),
            force_stages=set(getattr(args, 'force_stage', None).split(','))
                         if getattr(args, 'force_stage', None) else None,
            accept_prior_provenance=getattr(args, 'accept_prior_provenance', False),
        )
        logger.info(
            "[RESUME] sample=%s stage=correct decision=%s reason=%s",
            _sample_id,
            'SKIP' if _skip_decision.skip else 'RUN',
            _skip_decision.reason,
        )
        if _skip_decision.skip:
            logger.info("[RESUME] Skipping correct stage — prior sidecar valid.")
            return 0

        # Dry-run mode: print decision and exit.
        if getattr(args, 'dry_run_resume', False):
            print(f"[dry-run-resume] stage=correct decision=RUN reason={_skip_decision.reason}")
            return 0
    else:
        _rectify_sha = 'unknown'
    # -------------------------------------------------------

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
        # Guard: detect empty/headerless BAM (e.g. aligner SIGSEGV fallback) before
        # spending time on junction refinement or genome loading.
        try:
            import pysam as _pysam
            with _pysam.AlignmentFile(str(config['bam_path']), 'rb', check_sq=False) as _bam_probe:
                _n_refs = _bam_probe.nreferences
        except Exception:
            _n_refs = None
        if _n_refs is not None and _n_refs == 0:
            logger.warning(
                "Input BAM has no reference sequences in header — likely an aligner "
                "fallback empty BAM. Writing empty outputs and skipping correction."
            )
            _EMPTY_TSV_COLS = [
                'read_id', 'chrom', 'strand',
                'original_3prime', 'corrected_3prime',
                'five_prime_position', 'five_prime_rescued', 'five_prime_exon_cigar',
                'alignment_start', 'alignment_end',
                'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
                'polya_length', 'aligned_a_length', 'soft_clip_a_length',
                'junctions', 'n_junctions',
                'five_prime_soft_clip_length', 'three_prime_soft_clip_length',
                'mapq', 'correction_applied', 'confidence', 'qc_flags', 'fraction',
                'gene_id', 'pt_tag', 'polya_score', 'polya_source',
                'sc_homopolymer_extension', 'sc_rescued_seq', 'sc_original_softclip_len',
            ]
            import pandas as _pd
            _pd.DataFrame(columns=_EMPTY_TSV_COLS).to_csv(
                str(config['output_path']), sep='\t', index=False
            )
            if config.get('write_corrected_bam'):
                import shutil as _shutil
                _shutil.copy(str(config['bam_path']), str(config['write_corrected_bam']))
            return
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
        # Per-chromosome sorted junction index for Module 2F pool lookup.
        # Built from the prescan pool if available; otherwise None (GFF-only mode).
        _pool_chrom_index = None
        if _run_2h:
            _t_refine = _time.perf_counter()
            logger.info("Module 2H: Junction N-op boundary refinement...")
            try:
                from ..splice.junction_refiner import refine_bam_junctions
                from ..consensus.consensus import load_annotated_junctions as _load_annot_j

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

                from ...utils.genome import load_genome as _load_genome_for_refine
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

                # Build bisect index for Module 2F pool lookup (novel junction rescue).
                # This is built once here and passed to all process_bam_* calls below.
                _pool_chrom_index = bam_processor._build_pool_chrom_index(_prebuilt_pool)

                # For ONT PCR-cDNA, try to build a per-UMI-bin PenaltyTableSet from
                # bundled tables.  This routes each read to the table calibrated at
                # its UMI cluster depth (umi1/umi2/umi3plus), which reduces false
                # corrections from over-confident singleton-read error rates.
                # Only activated when no explicit --junction-penalty-table was given
                # (a user-supplied path always uses single-table mode).
                _penalty_table_set = None
                if config.get('ont_cDNA') and not config.get('junction_penalty_table'):
                    try:
                        from ..splice.junction_refiner import PenaltyTableSet
                        from ..splice.hp_penalty import HpPenaltyTable as _HpPT
                        from ...data import (
                            get_bundled_junction_penalty_table,
                            get_bundled_str_penalty_table as _get_str_pt,
                        )
                        _org = config.get('organism')
                        if _org is None:
                            from ...data import detect_organism_from_genome
                            _gp = config.get('genome_path')
                            if _gp:
                                _org = detect_organism_from_genome(Path(str(_gp)))
                        if _org:
                            _pts_bins = {}
                            _str_p = _get_str_pt(_org, protocol='cdna')
                            _str_path = str(_str_p) if _str_p is not None else None
                            for _bin_name in ('umi1', 'umi2', 'umi3plus', 'pooled'):
                                _p = get_bundled_junction_penalty_table(
                                    _org, protocol='cdna', bin=_bin_name
                                )
                                if _p is not None:
                                    _pts_bins[_bin_name] = _HpPT.from_tsv(
                                        str(_p), str_path=_str_path
                                    )
                            if all(b in _pts_bins for b in ('umi1', 'umi2', 'umi3plus')):
                                _penalty_table_set = PenaltyTableSet(_pts_bins, umi_tag='XC')
                                logger.info(
                                    "  Module 2H: per-UMI-bin penalty tables loaded "
                                    "(%d bins: %s)",
                                    len(_pts_bins), ', '.join(sorted(_pts_bins)),
                                )
                    except Exception as _pts_exc:
                        logger.debug(
                            "  Module 2H: could not build PenaltyTableSet "
                            "(will use single-table mode): %s", _pts_exc
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
                    sort_threads=config.get('threads', 1),
                    n_workers=config.get('threads', 1),
                    penalty_table_set=_penalty_table_set,
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
            from ..consensus.consensus import load_annotated_junctions
            annotated_junctions = load_annotated_junctions(str(config['annotation_path']))
            logger.info(
                f"Loaded {len(annotated_junctions):,} annotated junctions for 3'SS rescue "
                f"({_time.perf_counter() - _t_junc:.1f}s)"
            )

        # Build per-read gene attribution interval trees (optional, requires annotation)
        gene_interval_trees = None
        if config.get('annotation_path'):
            try:
                from ..analyze.gene_attribution import build_cds_interval_tree
                from ..analyze.loaders import load_annotation as _load_annotation_for_trees
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
        _protocol = 'dt_cdna' if is_dt_primed else ('ont_cdna' if is_ont_cdna else 'drs')

        if streaming_mode:
            if n_threads > 1:
                # Parallel streaming: region workers + stream output to disk
                # Best for large BAMs (> ~5M reads) — combines parallel throughput
                # with constant memory usage.
                variant_output_path = None
                if config['variant_aware'] and config['output_path']:
                    _op = str(config['output_path'])
                    if _op.endswith('.tsv'):
                        variant_output_path = _op.replace('.tsv', '_potential_variants.tsv')
                    else:
                        # output_path is a directory — place variant TSV inside it
                        variant_output_path = str(Path(_op) / 'corrected_3ends_potential_variants.tsv')
                stats = bam_parallel.process_bam_streaming_parallel(
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
                    pool_chrom_index=_pool_chrom_index,
                    gene_interval_trees=gene_interval_trees,
                    polya_model_path=_polya_model_path,
                    checkpoint_dir=getattr(args, 'checkpoint_dir', None),
                    variant_scan_cache=config.get('variant_scan_cache'),
                    dt_primed_cDNA=config.get('dt_primed_cDNA', False),
                    min_mapq=config.get('min_mapq', 0),
                    min_aligned_length=config.get('min_aligned_length', 0),
                )
            else:
                # Single-threaded streaming for single-core or debugging
                chunk_size = getattr(args, 'chunk_size', 10000)
                stats = bam_parallel.process_bam_streaming(
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
                    pool_chrom_index=_pool_chrom_index,
                    gene_interval_trees=gene_interval_trees,
                    polya_model_path=_polya_model_path,
                    dt_primed_cDNA=config.get('dt_primed_cDNA', False),
                    min_mapq=config.get('min_mapq', 0),
                    min_aligned_length=config.get('min_aligned_length', 0),
                )
            report = generate_stats_report(stats, protocol=_protocol)
        else:
            # Standard parallel processing
            # Determine variant output path
            variant_output_path = None
            if config['variant_aware'] and config['output_path']:
                _op = str(config['output_path'])
                if _op.endswith('.tsv'):
                    variant_output_path = _op.replace('.tsv', '_potential_variants.tsv')
                else:
                    variant_output_path = str(Path(_op) / 'corrected_3ends_potential_variants.tsv')

            results, stats = bam_parallel.process_bam_file_parallel(
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
                pool_chrom_index=_pool_chrom_index,
                gene_interval_trees=gene_interval_trees,
                polya_model_path=_polya_model_path,
                variant_scan_cache=config.get('variant_scan_cache'),
                dt_primed_cDNA=config.get('dt_primed_cDNA', False),
                min_mapq=config.get('min_mapq', 0),
                min_aligned_length=config.get('min_aligned_length', 0),
            )
            report = generate_stats_report(stats, protocol=_protocol)

        logger.info(f"[TIMING] BAM processing: {_time.perf_counter() - _t_proc:.1f}s")

        # ---------- Manifest TSV emission (Commit B, 2026-05-20) ----------
        # The canonical artifact going forward is corrected_reads.manifest.tsv.
        # The merged corrected_reads.tsv is emitted only when --emit-merged-tsv
        # is explicitly passed.
        _manifest_path = None
        if config.get('output_path'):
            _tsv_path = Path(str(config['output_path']))
            _manifest_path = _tsv_path.parent / (
                _tsv_path.stem.replace('corrected_reads', 'corrected_reads') + '.manifest.tsv'
            )
            # Build manifest path next to the TSV
            _manifest_path = _tsv_path.parent / (_tsv_path.stem + '.manifest.tsv')
            if _tsv_path.exists():
                import hashlib as _hashlib
                # Create a single-entry manifest wrapping the existing TSV.
                # (Commit D will emit per-region TSVs from the parallel workers;
                # for now we wrap the monolithic TSV in a one-entry manifest.)
                with _tsv_path.open('rb') as _mf:
                    _sha = _hashlib.sha256(_mf.read()).hexdigest()
                _n_data_rows = sum(1 for _ in _tsv_path.open()) - 1  # subtract header
                _rel_tsv = _os_run.path.relpath(str(_tsv_path), str(_manifest_path.parent))
                with open(_manifest_path, 'w') as _mfh:
                    _mfh.write('region_id\tchrom\tstart\tend\ttsv_path\tn_rows\tsha256\n')
                    _mfh.write(f'region_000\t*\t0\t0\t{_rel_tsv}\t{_n_data_rows}\t{_sha}\n')
                logger.info(
                    "[CHANGED] correct now writes corrected_reads.manifest.tsv by default; "
                    "pass --emit-merged-tsv for the legacy concatenated form or use "
                    "'rectify export-merged-tsv <manifest>' on demand."
                )
                logger.info(f"  Manifest: {_manifest_path}")
                if not _emit_merged_tsv:
                    # Default: manifest-only. Rename the merged TSV so it's
                    # not at the canonical path (available via manifest).
                    _region_tsv_path = _tsv_path.parent / (_tsv_path.stem + '.region_000.tsv')
                    _tsv_path.rename(_region_tsv_path)
                    # Update manifest to point to renamed file.
                    _rel_region = _os_run.path.relpath(str(_region_tsv_path), str(_manifest_path.parent))
                    with open(_manifest_path, 'w') as _mfh:
                        _mfh.write('region_id\tchrom\tstart\tend\ttsv_path\tn_rows\tsha256\n')
                        _mfh.write(f'region_000\t*\t0\t0\t{_rel_region}\t{_n_data_rows}\t{_sha}\n')
                    # Update config['output_path'] to point to the region TSV so
                    # downstream BAM writers (which load corrections from TSV) can
                    # still find it via the manifest-aware loader.
                    config['output_path'] = str(_region_tsv_path)
                    logger.info(f"  Region TSV (manifest-only mode): {_region_tsv_path}")
                else:
                    logger.info(f"  Merged TSV (--emit-merged-tsv): {_tsv_path}")
        # -------------------------------------------------------------------

        # Propagate spike-in count into stats
        if spikein_stats:
            stats.spikein_reads_filtered = spikein_stats.get('spikein_reads', 0)

        # Write processing statistics TSV
        if config['output_path']:
            stats_path = str(config['output_path']).replace('.tsv', '_stats.tsv')
            if stats_path == str(config['output_path']):
                stats_path = str(config['output_path']) + '_stats.tsv'
            write_stats_tsv(stats, stats_path, protocol=_protocol)
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
            _t_cbam = _time.perf_counter()
            logger.info(f"Writing corrected BAM to {config['corrected_bam']}...")

            # Commit B: parallel BAM write gate.
            _input_size = _os_run.path.getsize(str(bam_to_process))
            _n_reads_est = _estimate_n_reads(str(bam_to_process))
            _use_parallel = (
                not _legacy_single_threaded
                and _input_size >= 100 * 1024 ** 2   # 100 MB
                and _n_reads_est >= 500_000
            )
            if _use_parallel:
                logger.info(
                    "[bam_writer] parallel path: size=%d MB, est_reads=%d, n_threads=%d",
                    _input_size // (1024 ** 2), _n_reads_est, n_threads,
                )
                from ..bam.bam_writer_parallel import write_corrected_bam_parallel
                cbam_parallel_stats = write_corrected_bam_parallel(
                    input_bam_path=str(bam_to_process),
                    corrected_tsv_path=str(config['output_path']),
                    output_bam_path=str(config['corrected_bam']),
                    n_threads=n_threads,
                    genome=_genome_for_bam,
                    tmp_dir=_tmp_dir,
                )
                cbam_stats = {
                    'total': cbam_parallel_stats.get('n_reads_out_total', 0),
                    'clipped': 0,  # parallel path doesn't track per-category
                    'unchanged': 0,
                }
                logger.info(
                    f"  Reads written: {cbam_stats['total']:,}  (parallel, {cbam_parallel_stats.get('n_regions', '?')} regions)"
                )
            elif _legacy_single_threaded:
                logger.info(
                    "[bam_writer] using legacy single-threaded path (--legacy-single-threaded)"
                )
                _unsorted_bam = str(config['corrected_bam']) + '.unsorted.bam'
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
            else:
                logger.info(
                    "[bam_writer] small-input gate fired (size=%d MB, est_reads=%d) — using sequential path",
                    _input_size // (1024 ** 2), _n_reads_est,
                )
                _unsorted_bam = str(config['corrected_bam']) + '.unsorted.bam'
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
            # Use TSV stem ("corrected_reads") so the name doesn't depend on
            # whatever the BAM happens to be called (avoids redundant names
            # like "rectified_corrected_3end_corrected_reads" when the BAM
            # already contains "_corrected_3end" in its stem).
            _tsv_stem = Path(str(config['output_path'])).stem
            _c3bg_prefix = str(_bam_dir / _tsv_stem)
            logger.info(f"Writing corrected-3'-end bedgraph (prefix: {_c3bg_prefix})...")
            try:
                _c3bg_counts = bam_processor.write_corrected_reads_bedgraph(
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
            from ..netseq.netseq_output import write_bedgraph as _write_bedgraph
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

        _wall_total_secs = _time.perf_counter() - _t_correct_total
        logger.info("")
        logger.info("=" * 70)
        logger.info("RECTIFY completed successfully!")
        logger.info(f"[TIMING] Correction total: {_wall_total_secs:.1f}s")
        logger.info("=" * 70)

        # ---------- Stage-level sidecar emission (Commit B, 2026-05-20) ----------
        # CRITICAL: sidecar MUST be written AFTER all durable outputs exist.
        # A crash before this point leaves no sidecar → next run correctly reruns.
        if _output_dir is not None:
            try:
                from datetime import datetime as _dt_sidecar, timezone as _tz_sidecar_cls
                from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar

                _sidecar_inputs = {}
                if config.get('bam_path'):
                    _sidecar_inputs['bam'] = str(config['bam_path'])
                if config.get('genome_path'):
                    _sidecar_inputs['genome'] = str(config['genome_path'])
                if config.get('annotation_path'):
                    _sidecar_inputs['annotation'] = str(config['annotation_path'])

                _sidecar_outputs = {}
                if config.get('output_path') and Path(str(config['output_path'])).exists():
                    _sidecar_outputs['corrected_tsv'] = str(config['output_path'])
                if _manifest_path and _manifest_path.exists():
                    _sidecar_outputs['corrected_tsv_manifest'] = str(_manifest_path)
                if config.get('corrected_bam') and Path(str(config['corrected_bam'])).exists():
                    _sidecar_outputs['corrected_bam'] = str(config['corrected_bam'])
                if config.get('corrected_bam') and Path(str(config['corrected_bam']) + '.bai').exists():
                    _sidecar_outputs['corrected_bam_index'] = (str(config['corrected_bam']) + '.bai', True)

                # Determine protocol subtype
                _subtype = 'drs'
                if is_dt_primed:
                    _subtype = 'cdna'
                elif is_ont_cdna:
                    _subtype = 'ont_cdna'
                elif is_short_read:
                    _subtype = 'short_read'

                _record = ProvenanceRecord.from_components(
                    stage='correct',
                    stage_subtype=_subtype,
                    sample_id=_sample_id,
                    sample_output_dir=_output_dir,
                    started_at=_stage_started_at,
                    completed_at=_dt_sidecar.now(_tz_sidecar_cls.utc).isoformat(),
                    exit_status=0,
                    inputs=_sidecar_inputs,
                    outputs=_sidecar_outputs,
                    stats={
                        'wall_seconds': _wall_total_secs,
                        'n_threads': n_threads,
                    },
                    skip_check_config={
                        'ignore_argv': [
                            '--n-threads', '--threads', '-j',
                            '--tmp-dir', '--verbose', '--keep-tmp',
                            '--legacy-single-threaded',
                            '--checkpoint-dir', '--chunk-size',
                        ],
                    },
                    argv=_sys_run.argv,
                    rectify_git_sha=_rectify_sha,
                )
                write_stage_sidecar(_record, sample_output=_output_dir)
                logger.info(
                    "[PROVENANCE] Wrote stage sidecar: %s.correct.provenance.json",
                    _sample_id,
                )
            except Exception as _sc_exc:
                logger.warning(
                    "[PROVENANCE] Failed to write correct stage sidecar (non-fatal): %s",
                    _sc_exc,
                )
        # -------------------------------------------------------------------------

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


def create_correct_parser(subparsers):
    """Wire the `correct` subcommand into the given subparsers group."""
    import argparse
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
        '--dT-primed-cDNA',
        dest='dT_primed_cDNA',
        action='store_true',
        help='Input was generated with oligo-dT priming (QuantSeq, dT-primed cDNA). '
             'The poly(A) tail is NOT in the read. Enables indel artifact correction '
             'and poly(A) trimming modules designed for this protocol. '
             'For nanopore direct RNA-seq the poly(A) tail IS sequenced — do NOT use this flag.'
    )
    tech_group.add_argument(
        '--polya-sequenced',
        dest='dT_primed_cDNA',
        action='store_true',
        help=argparse.SUPPRESS,  # Deprecated alias for --dT-primed-cDNA
    )
    tech_group.add_argument(
        '--ONT-cDNA',
        dest='ONT_cDNA',
        action='store_true',
        default=False,
        help='Input is Oxford Nanopore PCR-cDNA (e.g. SQK-PCB114). The poly(A) tail IS '
             'present in the read as a 3\' soft-clip (minimap2 alignment). Strand '
             'convention is the same as DRS: is_reverse=True → minus-strand gene. '
             'Enables poly(A) trimming and indel correction; disables AG mispriming '
             '(no oligo-dT priming step). Do NOT use --dT-primed-cDNA for this protocol.'
    )
    tech_group.add_argument(
        '--short-read',
        dest='short_read',
        action='store_true',
        default=False,
        help=(
            'Input is short-read data (Illumina/Aviti ≤150 bp). Disables '
            'poly(A)-tail trimming, A-tract correction, and indel modules '
            '(which assume a poly(A) soft-clip). Use with `rectify align '
            '--short-read` to select bbmap + bwa instead of the long-read '
            'aligner panel.'
        ),
    )

    tech_group.add_argument(
        '--aligner',
        choices=['minimap2', 'bwa', 'star', 'auto'],
        default='auto',
        help='Aligner used (affects indel artifact detection)'
    )

    # Spike-in filtering
    spikein_group = correct_parser.add_argument_group('Spike-in filtering')
    spikein_group.add_argument(
        '--filter-spikein',
        nargs='+',
        metavar='GENE',
        help='Remove spike-in reads by gene name (e.g., --filter-spikein ENO2). '
             'Pre-defined signatures: ENO2. Multiple genes accepted.'
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
        help='Skip poly(A) tail trimming (even if --dT-primed-cDNA)'
    )

    module_group.add_argument(
        '--skip-indel-correction',
        action='store_true',
        help='Skip indel artifact correction (even if --dT-primed-cDNA)'
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

    module_group.add_argument(
        '--min-mapq',
        type=int,
        default=0,
        metavar='MAPQ',
        help='Minimum mapping quality (MAPQ) to include a read. Reads with '
             'MAPQ < MIN_MAPQ are skipped before any correction. Default: 0 '
             '(no filter). Use --min-mapq 1 to exclude multi-mapping reads '
             '(MAPQ=0), which is recommended for dT-primed-cDNA data where '
             'internally-primed reads from repetitive T-rich regions are common.'
    )
    module_group.add_argument(
        '--min-aligned-length',
        type=int,
        default=0,
        metavar='BP',
        help='Minimum number of reference bases consumed by the alignment '
             '(reference_length) to include a read. Reads with fewer aligned '
             'bases are skipped before correction. Default: 0 (no filter). '
             'Use --min-aligned-length 30 for dT-primed-cDNA short-read data '
             'to remove reads that align only to low-complexity T-tracts with '
             'insufficient unique sequence context.'
    )

    # Poly(A) model
    polya_group = correct_parser.add_argument_group('Poly(A) tail model')
    polya_group.add_argument(
        '--polya-model',
        type=Path,
        help='Pre-trained poly(A) tail model (JSON). If not provided, uses built-in model.'
    )
    polya_group.add_argument(
        '--min-polya-score',
        type=float,
        default=None,
        metavar='SCORE',
        help='Minimum poly(A) model confidence score (0-1) to write to polya_pass column. '
             'Reads below this threshold are flagged polya_pass=0 but still written to output. '
             'Requires --polya-model or --dT-primed-cDNA to compute scores.'
    )

    # NET-seq refinement
    netseq_group = correct_parser.add_argument_group('NET-seq refinement')
    from rectify.data import add_organism_args
    add_organism_args(netseq_group)

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
        default=17.0,
        help='AG-richness weighted score threshold for mispriming flagging (0.0-34.5; default 17.0)'
    )

    # Output options
    output_group = correct_parser.add_argument_group('Output options')
    output_group.add_argument(
        '--output-bam',
        type=Path,
        metavar='BAM',
        help='Write a poly(A)-trimmed BAM alongside the corrected TSV. '
             'The 3\' poly(A) soft-clip is removed from each read; all '
             'other alignment fields and BAM tags are preserved unchanged. '
             'Requires --dT-primed-cDNA.'
    )

    output_group.add_argument(
        '--write-corrected-bam',
        dest='corrected_bam',
        type=Path,
        metavar='BAM',
        help='Write a corrected BAM where every read is hard-clipped at its '
             'corrected 3\' end. Unlike --output-bam, this reflects all '
             'correction categories (Cat1 poly-A walkback, Cat4 false-N '
             'absorption, Cat2 soft-clip rescue, etc.) — not just soft-clip '
             'trimming. The output is sorted and indexed. Best used with '
             '--annotation for full correction coverage.'
    )

    output_group.add_argument(
        '--write-softclipped-bam',
        dest='softclipped_bam',
        type=Path,
        metavar='BAM',
        help='Write a corrected BAM where poly(A) tail bases are soft-clipped '
             'rather than hard-clipped at the corrected 3\' end. The poly(A) '
             'bases remain in the query sequence and are visible in IGV when '
             '"Show soft-clipped bases" is enabled. Useful for visualising '
             'where the poly(A) tail was without affecting pileup or coverage. '
             'The output is sorted and indexed.'
    )

    output_group.add_argument(
        '--write-bedgraph',
        dest='bedgraph_prefix',
        metavar='PREFIX',
        help='Write strand-specific bedGraph files for NET-seq reads. '
             'Signal at each position is the sum of fractional contributions '
             'across all reads (the fraction column in the corrected TSV), so '
             'each read contributes exactly 1.0 of total signal distributed '
             'across its candidate positions. Output files: '
             'PREFIX.plus.bedgraph and PREFIX.minus.bedgraph. '
             'Requires NET-seq input (fraction column must be present).'
    )

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

    perf_group.add_argument(
        '--checkpoint-dir',
        type=str,
        default=None,
        metavar='DIR',
        help='Directory to store per-region checkpoint sentinels and the variant-scan pickle. '
             'When set, a failed run can be resumed: completed regions are skipped, '
             'the variant scan is reloaded from disk, and the partial output TSV is '
             'appended to rather than overwritten. Has no effect without --streaming.'
    )

    perf_group.add_argument(
        '--tmp-dir',
        type=str,
        default=None,
        metavar='DIR',
        help='Scratch directory for intermediate per-region BAM files used by the '
             'parallel BAM writer (write_corrected_bam_parallel). Defaults to a new '
             'directory under $TMPDIR / tempfile.gettempdir(). On H2/Sherlock, pass '
             '$L_SCRATCH/rectify_regions for fast local disk.'
    )

    perf_group.add_argument(
        '--emit-merged-tsv',
        dest='emit_merged_tsv',
        action='store_true',
        default=False,
        help='Also emit the legacy concatenated corrected_reads.tsv alongside the '
             'corrected_reads.manifest.tsv (default: manifest-only). Useful for '
             'downstream scripts that have not yet been migrated to accept manifests. '
             'Equivalent to running rectify export-merged-tsv <manifest> afterwards.'
    )

    perf_group.add_argument(
        '--legacy-single-threaded',
        dest='legacy_single_threaded',
        action='store_true',
        default=False,
        help='Force the pre-Commit-A single-threaded BAM writer regardless of input '
             'size or thread count. For debugging correctness regressions and '
             'memory-constrained environments. '
             'DEPRECATED: will be removed in a future release.'
    )

    # Resume / sidecar flags
    from rectify.core.provenance.cli import add_resume_args
    add_resume_args(correct_parser)
    perf_group.add_argument(
        '--variant-scan-cache',
        type=str,
        default=None,
        metavar='PKL',
        help='Path to a pre-computed variant scan pickle produced by `rectify prescan`. '
             'Skips the Pass 1 all-reads scan (Module 2D) and loads the '
             'VariantAwareHomopolymerRescue object directly. Use when running per-chunk '
             'correction with a shared scan built from the merged BAM.'
    )
    perf_group.add_argument(
        '--junction-pool-cache',
        type=str,
        default=None,
        metavar='PKL',
        help='Path to a pre-computed junction pool pickle produced by `rectify prescan`. '
             'Skips the aligner-BAM junction collection step (Module 2H) and loads the '
             'pool directly. The pickle must contain keys "all_junctions" and '
             '"annotated_set". Use when running per-chunk correction with a shared pool '
             'built from all aligner BAMs.'
    )

    # Junction refinement options
    junc_group = correct_parser.add_argument_group('Junction refinement')
    junc_group.add_argument(
        '--aligner-bams',
        action='append',
        type=Path,
        metavar='BAM',
        default=[],
        help='Per-aligner BAM file for junction candidate pool construction '
             '(repeat the flag for each aligner: --aligner-bams minimap2.bam '
             '--aligner-bams gapmm2.bam ...). Accepts plain paths or '
             '\'aligner:path\' pairs. When provided, every N-op in the '
             'consensus BAM is scored against all candidate junctions within '
             '--junction-search-radius bp and replaced with the best-supported '
             'junction (canonical GT-AG > annotated > highest split-alignment '
             'score). Improves junction accuracy for reads where the chosen '
             'aligner placed the intron boundary a few bp off from the true '
             'splice site.'
    )
    junc_group.add_argument(
        '--junction-hp-pen',
        type=float,
        default=0.25,
        metavar='FLOAT',
        help='Homopolymer indel penalty for junction scoring (0 < value ≤ 1.0). '
             'Lower values are more tolerant of poly-T/A undercalling errors near '
             'splice sites. Default 0.25 is recommended for Nanopore DRS data.'
    )
    junc_group.add_argument(
        '--junction-search-radius',
        type=int,
        default=5000,
        metavar='BP',
        help='Search radius (bp) around each N-op for candidate junction lookup '
             '(default: 5000). Covers the longest known yeast intron (~1 kb) '
             'with ample margin.'
    )
    junc_group.add_argument(
        '--junction-window',
        type=int,
        default=15,
        metavar='BP',
        help='Edit-distance window half-width for split-alignment scoring '
             '(default: 15 bp on each side of the junction).'
    )
    junc_group.add_argument(
        '--junction-max-slide',
        type=int,
        default=10,
        metavar='BP',
        help='Maximum query split displacement for split-alignment scoring '
             '(default: ±10 bp).'
    )
    junc_group.add_argument(
        '--junction-max-boundary-shift',
        dest='junction_max_boundary_shift',
        type=int,
        default=50,
        metavar='BP',
        help='Maximum allowed shift of either intron boundary when applying a '
             'junction replacement (default: 50 bp).  Candidates whose '
             'intron_start or intron_end differs from the current N-op by more '
             'than this threshold are excluded, preventing false-positive '
             'matches from junctions in neighbouring genes.'
    )
    junc_group.add_argument(
        '--junction-penalty-table',
        dest='junction_penalty_table',
        default=None,
        metavar='PATH',
        help='Path to empirical penalty table (penalty_scores.tsv) produced by '
             'empirical_cigar_error_profiler.py.  When provided, overrides the '
             'heuristic del/ins cost step-function with per-HP-length values '
             'derived from multi-aligner agreement on this dataset.  Recommended '
             'for fine-tuning junction scoring to the specific sequencing run.'
    )
    junc_group.add_argument(
        '--str-penalty-table',
        dest='str_penalty_table',
        default=None,
        metavar='PATH',
        help='Path to STR penalty table (str_penalty_scores.tsv) produced by '
             'empirical_cigar_error_profiler.py --str-repeat.  When provided '
             'alongside --junction-penalty-table, dinucleotide/trinucleotide '
             'repeat contexts (e.g. TATATA→TATA slippage) are scored with '
             'empirical STR-specific penalties instead of the HP=1 baseline.'
    )
