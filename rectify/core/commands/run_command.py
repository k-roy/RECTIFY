#!/usr/bin/env python3
"""
RECTIFY run-all command — complete end-to-end pipeline.

Handles everything from FASTQ (or BAM) to final results:

  Step 0  Alignment    triple-aligner consensus (minimap2 + mapPacBio + gapmm2)
                       Skipped automatically if {sample}.consensus.bam already exists,
                       or if input is already a BAM, or if --skip-alignment is set.

  Step 1  Correction   3' end correction (poly(A) trimming, indel correction,
                       NET-seq refinement, spike-in filtering)

  Step 2  Analysis     Clustering, DESeq2, GO enrichment, motif discovery
                       For multi-sample runs: runs after ALL samples are corrected
                       so DESeq2 has full statistical power across conditions.

Usage:
    # Single sample from FASTQ
    rectify run-all sample.fastq.gz --genome genome.fa --annotation genes.gff -o results/

    # Single sample from BAM (alignment skipped)
    rectify run-all sample.bam --genome genome.fa --annotation genes.gff -o results/

    # Multi-sample (manifest) — parallel correction + combined DESeq2
    rectify run-all --manifest manifest.tsv --genome genome.fa --annotation genes.gff -o results/

    # Multi-sample with SLURM
    rectify run-all --manifest manifest.tsv --genome genome.fa --annotation genes.gff -o results/ \\
        --profile my_cluster.yaml

This module is a thin orchestrator: ``run()`` dispatches to one of the runners
in the ``run`` sub-subpackage:

- ``run.helpers``       — reference resolution, BAM checks, canonical paths.
- ``run.stages``        — per-stage runners (align / correct / combine /
                          analyse / junctions).
- ``run.single_sample`` — single-sample pipeline + per-sample worker.
- ``run.multi_sample``  — manifest-driven multi-sample pipeline.
- ``run.chunked_batch`` — ``--chunked-alignment`` shell-script generator.

The names that external callers and the test suite reach in for
(``_run_analysis``, ``_combine_corrected_tsvs``, etc.) are re-exported below
so the move is invisible to ``batch_command.py`` and
``tests/test_run_command_wiring.py``.

Author: Kevin R. Roy
Date: 2026-03-28
"""

import argparse
import sys
from pathlib import Path

# Re-exports for backwards compatibility with batch_command.py and
# tests/test_run_command_wiring.py, which reach in via
# ``rectify.core.commands.run_command.<name>``.
from .run.helpers import (  # noqa: F401
    _bam_has_md_tags,
    _collect_per_aligner_bams,
    _rectified_bam_path,
    _resolve_reference_paths,
    _validate_bam_integrity,
)
from .run.stages import (  # noqa: F401
    _combine_corrected_tsvs,
    _run_alignment,
    _run_analysis,
    _run_correction,
    _run_correction_per_aligner,
    _run_junction_aggregation,
)
from .run.single_sample import (  # noqa: F401
    _process_one_sample,
    _run_single_sample,
)
from .run.multi_sample import (  # noqa: F401
    _run_analysis_manifest,
    _run_multi_sample,
)
from .run.chunked_batch import _generate_chunked_pipeline  # noqa: F401


__all__ = [
    'run',
    'create_run_parser',
    # Backwards-compat re-exports
    '_bam_has_md_tags',
    '_collect_per_aligner_bams',
    '_combine_corrected_tsvs',
    '_generate_chunked_pipeline',
    '_process_one_sample',
    '_rectified_bam_path',
    '_resolve_reference_paths',
    '_run_alignment',
    '_run_analysis',
    '_run_analysis_manifest',
    '_run_correction',
    '_run_correction_per_aligner',
    '_run_junction_aggregation',
    '_run_multi_sample',
    '_run_single_sample',
    '_validate_bam_integrity',
]


def run(args: argparse.Namespace) -> None:
    """Dispatch to single-sample or multi-sample pipeline."""
    from rectify.slurm import set_thread_limits
    set_thread_limits(getattr(args, 'threads', None))

    if getattr(args, 'manifest', None) and getattr(args, 'input', None):
        print(
            "ERROR: Cannot combine a positional input file with --manifest. "
            "Provide one or the other.",
            file=sys.stderr,
        )
        sys.exit(1)

    # --chunked-alignment: generate scripts and exit (don't run inline),
    # UNLESS the input is a BAM (chunked alignment not needed → run inline).
    if getattr(args, 'chunked_alignment', False):
        _resolve_reference_paths(args)
        # Check if it's a BAM — if so _generate_chunked_pipeline prints a warning
        # and returns 0; we then fall through to inline correction.
        input_path_str = getattr(args, 'input', None)
        is_bam_input = False
        if input_path_str:
            from ..align.preprocess import detect_input_type
            itype = detect_input_type(Path(str(input_path_str)))
            is_bam_input = itype not in ('fastq', 'fastq.gz')

        if is_bam_input:
            print(
                f"Warning: --chunked-alignment is only needed for FASTQ inputs. "
                f"Input is a BAM — proceeding with inline BAM correction.",
                file=sys.stderr,
            )
            # Fall through to normal single-sample processing below
        else:
            rc = _generate_chunked_pipeline(args)
            sys.exit(rc)

    if getattr(args, 'manifest', None):
        rc = _run_multi_sample(args)
    elif getattr(args, 'input', None):
        rc = _run_single_sample(args)
    else:
        print(
            "ERROR: Provide either an input file (FASTQ or BAM) or --manifest.",
            file=sys.stderr,
        )
        rc = 1

    sys.exit(rc)


def create_run_parser(subparsers):
    """Wire the `run-all` subcommand into the given subparsers group."""
    import argparse
    # =========================================================================
    # run-all command (all-in-one: correct + analyze)
    # =========================================================================
    run_parser = subparsers.add_parser(
        'run-all',
        help='Complete end-to-end pipeline: align → correct → analyze',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
    Examples:
      # S. cerevisiae — no need to specify genome/annotation (bundled)
      rectify run-all sample.fastq.gz --Scer -o results/sample/

      # Single sample with explicit references
      rectify run-all sample.fastq.gz --genome genome.fa --annotation genes.gff -o results/sample/

      # DRS: Dorado-aligned BAM — poly(A) trim (Step 0) + restore soft-clips (Step 4)
      rectify run-all sample_dorado.bam --drs --Scer -o results/sample/

      # Single sample from pre-aligned BAM (skips alignment)
      rectify run-all sample.bam --Scer -o results/sample/

      # Multi-sample from manifest (parallel correction + combined DESeq2)
      rectify run-all --manifest manifest.tsv --Scer -o results/

      # Multi-sample with SLURM
      rectify run-all --manifest manifest.tsv --Scer -o results/ \\
      --profile my_cluster.yaml

    Manifest format (TSV):
      sample_id    path                    condition
      wt_rep1      /path/wt_rep1.fastq.gz  WT
      ko_rep1      /path/ko_rep1.fastq.gz  KO
        """
    )

    # Input: single file (positional, optional) OR --manifest
    # Validated at runtime — exactly one must be provided.
    run_parser.add_argument(
        'input',
        type=Path,
        nargs='?',
        default=None,
        help='Input FASTQ.gz or BAM file (single sample). '
             'Omit when using --manifest.'
    )
    run_parser.add_argument(
        '--manifest', '-m',
        type=Path,
        help='Sample manifest TSV (sample_id, path, condition). '
             'Enables multi-sample parallel correction and combined DESeq2. '
             'Cannot be combined with a positional input file.'
    )

    run_parser.add_argument(
        '--genome',
        type=Path,
        default=None,
        help='Reference genome FASTA file. Not required when --Scer or --organism is set.'
    )

    run_parser.add_argument(
        '--annotation',
        type=Path,
        default=None,
        help='Gene annotation file (GTF/GFF). Not required when --Scer or --organism is set.'
    )

    run_parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        required=True,
        help='Output directory. Single sample: all outputs here. '
             'Multi-sample: per-sample subdirs + combined/ for DESeq2.'
    )

    run_parser.add_argument(
        '--profile',
        type=Path,
        metavar='PROFILE_YAML',
        help='SLURM profile YAML for cluster submission (see rectify/slurm_profiles/)'
    )

    run_parser.add_argument(
        '--drs',
        action='store_true',
        default=False,
        dest='drs',
        help=(
            'Input is a Dorado direct RNA-seq BAM. Automatically runs '
            'poly(A)+adapter pre-trimming (Step 0) before alignment and '
            'restores trimmed bases as soft-clips (Step 4) after correction. '
            'Has no effect on FASTQ inputs (assumed already trimmed).'
        ),
    )
    run_parser.add_argument(
        '--short-read',
        dest='short_read',
        action='store_true',
        default=False,
        help=(
            'Input is short-read data (Illumina/Aviti ≤150 bp). Uses bbmap + bwa '
            'instead of the long-read aligner panel (minimap2/mapPacBio/gapmm2) and '
            'disables poly(A)-tail trimming, A-tract correction, and indel modules.'
        ),
    )

    run_parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip triple-aligner alignment even for FASTQ input '
             '(use if you already have a rectified.bam or consensus.bam)'
    )

    run_parser.add_argument(
        '--chunked-alignment',
        dest='chunked_alignment',
        action='store_true',
        default=False,
        help=(
            'Generate scheduler array scripts for chunked parallel alignment instead of '
            'running alignment inline. Auto-sizes chunks from --target-reads-per-chunk. '
            'Run bash submit_pipeline.sh to launch the generated dependency chain.'
        )
    )
    run_parser.add_argument(
        '--force-no-chunking',
        dest='force_no_chunking',
        action='store_true',
        default=False,
        help=(
            'DANGER: Bypass the mandatory chunking requirement for FASTQ inputs. '
            'Running whole-FASTQ alignment inline causes severe NFS I/O contention on '
            'HPC clusters and is 10-100x slower than chunked array jobs. '
            'Only use this on a local workstation (not an HPC cluster).'
        )
    )
    run_parser.add_argument(
        '--target-reads-per-chunk',
        dest='target_reads_per_chunk',
        type=int,
        default=500_000,
        metavar='READS',
        help='Target reads per alignment chunk when --chunked-alignment is set (default: 500000)'
    )
    run_parser.add_argument(
        '--scheduler',
        choices=['slurm', 'uge', 'pbs'],
        default='slurm',
        help='Target scheduler for generated script headers (default: slurm)'
    )
    run_parser.add_argument(
        '--slurm-partition',
        default=None,
        dest='slurm_partition',
        help='SLURM partition(s) for generated scripts (used with --scheduler slurm)'
    )
    run_parser.add_argument(
        '--slurm-account',
        default=None,
        dest='slurm_account',
        help='SLURM account for generated scripts (used with --scheduler slurm)'
    )
    run_parser.add_argument(
        '--uge-queue',
        default='long.q',
        help='UGE/SGE queue name (used with --scheduler uge)'
    )
    run_parser.add_argument(
        '--uge-pe',
        default='smp',
        help='UGE/SGE parallel environment (used with --scheduler uge)'
    )
    run_parser.add_argument(
        '--pbs-queue',
        default='workq',
        help='PBS/Torque queue name (used with --scheduler pbs)'
    )
    run_parser.add_argument(
        '--python-path',
        default=None,
        help='Explicit path to conda Python for generated scripts (default: sys.executable)'
    )
    run_parser.add_argument(
        '--rectify-src',
        default=None,
        help='Path to RECTIFY source checkout for generated scripts (default: auto-detected)'
    )

    # Organism / bundled reference
    from rectify.data import add_organism_args
    add_organism_args(run_parser)

    # Optional correction arguments

    run_parser.add_argument(
        '--filter-spikein',
        nargs='+',
        metavar='GENE',
        help='Remove spike-in reads by gene name (e.g., --filter-spikein ENO2). '
             'Uses sequence-based classification to distinguish spike-in from endogenous reads.'
    )

    run_parser.add_argument(
        '--netseq-dir',
        type=Path,
        help='Custom NET-seq directory (overrides bundled data for mutant-specific analysis)'
    )

    run_parser.add_argument(
        '--aligner',
        choices=['minimap2', 'bwa', 'star', 'auto'],
        default='auto',
        help='Aligner used for BAM file'
    )

    run_parser.add_argument(
        '--dT-primed-cDNA',
        dest='dT_primed_cDNA',
        action='store_true',
        default=False,
        help='Input is dT-primed cDNA (QuantSeq, etc.) — poly(A) NOT in read. '
             'Enables indel artifact correction and poly(A) trimming modules. '
             'By default run-all assumes nanopore direct RNA (poly-A IS in read).'
    )
    run_parser.add_argument(
        '--no-polya-sequenced',
        dest='dT_primed_cDNA',
        action='store_true',
        default=False,
        help=argparse.SUPPRESS,  # Deprecated alias for --dT-primed-cDNA
    )

    # Optional analysis arguments
    run_parser.add_argument(
        '--reference',
        help='Reference condition for DESeq2 (auto-detected if not specified)'
    )

    run_parser.add_argument(
        '--go-annotations',
        type=Path,
        help='GO annotation file for enrichment analysis'
    )

    run_parser.add_argument(
        '--threads',
        type=int,
        default=0,
        help='Number of threads (0 = auto-detect from SLURM_CPUS_PER_TASK or CPU count)'
    )

    run_parser.add_argument(
        '--streaming',
        action='store_true',
        help='Use streaming output mode to minimise memory usage during correction. '
             'Recommended for BAM files > 10 GB.'
    )

    run_parser.add_argument(
        '--chunk-size',
        type=int,
        default=10000,
        help='Reads per output chunk in streaming mode (default: 10000)'
    )

    run_parser.add_argument(
        '--parallel-aligners',
        action='store_true',
        default=False,
        help=(
            'Run minimap2, mapPacBio, and gapmm2 in parallel during alignment. '
            'Threads are divided evenly across aligners (e.g., --threads 16 → 5 per aligner). '
            'Reduces wall-clock time for the alignment step at the cost of higher peak memory.'
        )
    )

    run_parser.add_argument(
        '--mapPacBio-chunks',
        type=int,
        default=1,
        metavar='N',
        help=(
            'Number of mapPacBio FASTQ chunks. When N > 1, run-all merges '
            'existing chunk BAMs ({sample}.mapPacBio.chunk_*_of_N.bam) instead '
            'of re-running mapPacBio. Submit the chunk array job first with '
            'rectify align --mapPacBio-chunks N --mapPacBio-chunk-idx K, then '
            'run-all with --dependency and this flag.'
        )
    )

    run_parser.add_argument(
        '--base-aligners',
        nargs='+',
        choices=['minimap2', 'mapPacBio', 'gapmm2', 'bbmap', 'bwa'],
        default=None,
        metavar='ALIGNER',
        dest='base_aligners',
        help=(
            'Base aligners for the consensus pool. When omitted, the default '
            'depends on --short-read: bbmap + bwa for short-read, '
            'minimap2 + mapPacBio + gapmm2 for long-read. '
            'Passing this flag overrides the default in both modes. '
            'Example (force short-read panel even without --short-read): '
            '--base-aligners bbmap bwa'
        )
    )

    run_parser.add_argument(
        '--junction-aligners',
        nargs='+',
        choices=['uLTRA', 'deSALT'],
        default=None,
        metavar='ALIGNER',
        help=(
            'Junction-aware aligners for the consensus pool '
            '(choices: uLTRA, deSALT). Requires --annotation. '
            'When omitted, defaults to [] (disabled) under --short-read and '
            '[uLTRA, deSALT] under long-read. '
            'Pass --no-junction-aligners to explicitly disable.'
        )
    )

    run_parser.add_argument(
        '--no-junction-aligners',
        dest='junction_aligners',
        action='store_const',
        const=[],
        help='Disable uLTRA and deSALT (use only the --base-aligners set).'
    )

    run_parser.add_argument(
        '--chimeric-consensus',
        action='store_true',
        default=True,
        help=(
            'Use chimeric consensus selection: independently pick the best aligner '
            'for each read segment. Enabled by default. Use --no-chimeric-consensus to disable.'
        )
    )

    run_parser.add_argument(
        '--no-chimeric-consensus',
        dest='chimeric_consensus',
        action='store_false',
        help='Disable chimeric consensus selection (use single-best-aligner mode).'
    )

    run_parser.add_argument(
        '--ultra-path',
        default='uLTRA',
        help='Path to uLTRA executable (used with --junction-aligners uLTRA)'
    )

    run_parser.add_argument(
        '--desalt-path',
        default='deSALT',
        help='Path to deSALT executable (used with --junction-aligners deSALT)'
    )

    # BAM output options
    bam_group = run_parser.add_argument_group('BAM output options')
    bam_group.add_argument(
        '--bam-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Directory to write alignment BAMs (per-aligner and rectified). '
             'Defaults to the sample output directory. '
             'Useful for inspecting per-aligner BAMs separately from corrected outputs.'
    )
    bam_group.add_argument(
        '--keep-aligner-bams',
        action='store_true',
        default=False,
        help='Retain per-aligner BAMs (minimap2, mapPacBio, gapmm2) after '
             'consensus selection. By default, per-aligner BAMs are excluded '
             'from the Oak sync to save disk space; only the rectified BAM is kept.'
    )

    run_parser.add_argument(
        '--continue-on-error',
        action='store_true',
        default=False,
        help='Continue processing remaining samples if one fails'
    )

    run_parser.add_argument(
        '--use-scratch',
        action='store_true',
        default=False,
        help='Stage I/O through $SCRATCH for better performance'
    )

    run_parser.add_argument(
        '--junction-penalty-table',
        dest='junction_penalty_table',
        default=None,
        metavar='PATH',
        help='Path to empirical HP-context penalty table (penalty_scores.tsv) produced by '
             'empirical_cigar_error_profiler.py. Passed through to rectify correct. '
             'Overrides heuristic del/ins costs with per-HP-length values derived from '
             'multi-aligner agreement on this dataset.'
    )

    run_parser.add_argument(
        '--str-penalty-table',
        dest='str_penalty_table',
        default=None,
        metavar='PATH',
        help='Path to STR penalty table (str_penalty_scores.tsv) produced by '
             'empirical_cigar_error_profiler.py --str-repeat. Passed through to '
             'rectify correct alongside --junction-penalty-table.'
    )

    run_parser.add_argument(
        '--junction-overhang-table',
        dest='junction_overhang_table',
        default=None,
        metavar='PATH',
        help='Path to empirical junction overhang table (overhang_table.tsv) produced by '
             'rectify/core/calibrate_junction_overhang.py. When provided, aligner×read '
             'pairs whose intron junctions lack sufficient flanking overhang for their '
             'intron size are penalised (sorted last) during winner selection. Short introns '
             '(< 500 bp) with high cross-read support are exempt from strict filtering. '
             'If absent, all junctions are treated as plausible regardless of overhang. '
             'Use OverhangTable.default() behaviour by passing the path to a pre-calibrated '
             'table, or generate one from a completed reference run with '
             'calibrate_junction_overhang.py.'
    )

    run_parser.add_argument(
        '--write-softclip-bam',
        action='store_true',
        default=False,
        dest='write_softclip_bam',
        help='Write a soft-clipped corrected BAM alongside the primary hard-clipped BAM. '
             'Cat2 soft-clip rescue bases are visible in IGV with "Show soft-clipped bases" '
             'enabled. Off by default — useful for QC/debugging Cat2 rescues only.'
    )

    run_parser.add_argument(
        '--write-polya-bam',
        action='store_true',
        default=False,
        dest='write_polya_bam',
        help='Write a BAM with the original poly(A) tail restored from the DRS trim '
             'parquet metadata as a 3\' soft clip. Off by default — useful for visually '
             'validating poly(A) tail handling in IGV. Has no effect without --drs.'
    )
