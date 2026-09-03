#!/usr/bin/env python3
"""
RECTIFY run-all command — complete end-to-end pipeline.

Handles everything from FASTQ (or BAM) to final results:

  Step 0  Alignment    multi-aligner panel (minimap2 + uLTRA + deSALT +
                       overhang resolver on the minimap2 arm)
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

from .run.netseq_pipeline import CHURCHMAN_3P_LINKER
from typing import Optional

# Re-exports for backwards compatibility with batch_command.py and
# tests/test_run_command_wiring.py, which reach in via
# ``rectify.core.commands.run_command.<name>``.
from .run.helpers import (  # noqa: F401
    _bam_has_md_tags,
    _collect_per_aligner_bams,
    _emit_legacy_consensus_symlink,
    _final_rectified_bam_path,
    _multialigned_bam_path,
    _rectified_bam_path,  # back-compat alias for _multialigned_bam_path
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
    '_validate_paired_end_args',
    '_validate_protocol_flags',
    # Backwards-compat re-exports
    '_bam_has_md_tags',
    '_collect_per_aligner_bams',
    '_combine_corrected_tsvs',
    '_generate_chunked_pipeline',
    '_process_one_sample',
    '_emit_legacy_consensus_symlink',
    '_final_rectified_bam_path',
    '_multialigned_bam_path',
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


#: run-all protocol flags, in help order. Each selects a different strand
#: convention and module set, so exactly one may be active.
#:
#: 🔴 ``--short-read`` is NOT one of them. It is a READ-LENGTH MODALITY that composes with a
#: protocol -- ``--short-read --dT-primed-cDNA`` (QuantSeq) is the documented QuantSeq invocation,
#: and treating it as a fourth mutually exclusive protocol made that command exit 1, i.e. made the
#: whole run-all QuantSeq path unreachable (planning 832 G-0, verified by executing the validator).
#: The chunked-alignment guard below already had to special-case it out for the same reason.
_PROTOCOL_FLAGS = (
    ('--drs', 'drs'),
    ('--dT-primed-cDNA', 'dT_primed_cDNA'),
    ('--ONT-cDNA', 'ONT_cDNA'),
    ('--netseq', 'netseq'),
)

#: Read-length modality flags -- orthogonal to the protocol, never mutually exclusive with one.
_MODALITY_FLAGS = (
    ('--short-read', 'short_read'),
)

#: Every run-all flag that must reach ``rectify correct``'s Namespace. ``stages._run_correction``
#: hand-builds that Namespace, so anything missing here is silently dropped: ``--short-read`` was,
#: which left poly(A) trimming and the position-shifting indel module ON for 150-bp Illumina reads
#: (832 G-1 / 833 C-1), and ``--netseq`` would have been, which would have run a NET-seq library
#: under the SENSE strand rule (834 §6.4).
_CORRECT_FORWARDED_FLAGS = ('drs', 'dT_primed_cDNA', 'ONT_cDNA', 'short_read', 'netseq')


def _validate_protocol_flags(args: argparse.Namespace) -> Optional[str]:
    """Reject more than one protocol flag.

    The protocols select mutually exclusive strand rules. Silently letting one
    win is how a library gets corrected under the wrong convention — e.g. ONT
    PCR-cDNA under ``--dT-primed-cDNA``'s fixed antisense rule, which puts the
    3' end at the wrong terminus for roughly half the reads.

    Returns an error message, or None when at most one protocol is active.
    Pure (no I/O, no sys.exit) so it is directly unit-testable.
    """
    active = [flag for flag, dest in _PROTOCOL_FLAGS if getattr(args, dest, False)]
    if len(active) <= 1:
        return None
    return (
        f"{' and '.join(active)} are mutually exclusive protocols — pass exactly "
        "one. For Oxford Nanopore PCR-cDNA (SQK-PCB114) use --ONT-cDNA; "
        "--dT-primed-cDNA is for oligo-dT protocols (QuantSeq) and applies a "
        "fixed antisense strand rule that is wrong for PCR-cDNA."
    )


def _validate_paired_end_args(args: argparse.Namespace) -> Optional[str]:
    """Validate the paired-end short-read (--read2) argument combination.

    Returns an error message describing the first problem found, or None when
    the combination is valid (including the common case of no --read2 at all).

    Kept pure (no I/O, no sys.exit) so it is directly unit-testable: input-type
    detection is extension-based, so no file needs to exist.
    """
    read2 = getattr(args, 'read2', None)
    if read2 is None:
        return None

    if not getattr(args, 'short_read', False):
        return (
            "--read2 requires --short-read. Paired-end alignment is only "
            "supported for short-read data (the COMPASS panel); long-read "
            "protocols (DRS / PCR-cDNA) are single-end."
        )

    if getattr(args, 'manifest', None):
        return (
            "--read2 cannot be combined with --manifest: the manifest schema "
            "(sample_id, path, condition) has no mate-2 column, so a single "
            "--read2 cannot describe multiple samples. Run each paired sample "
            "as its own `rectify run-all --short-read --read2 R2.fastq.gz "
            "R1.fastq.gz` invocation."
        )

    input_path = getattr(args, 'input', None)
    if input_path is not None:
        from ..align.preprocess import detect_input_type
        itype = detect_input_type(Path(str(input_path)))
        if itype in ('bam', 'sam'):
            return (
                f"--read2 is not valid with an aligned {itype.upper()} input "
                f"({input_path}): mates are already paired inside the "
                f"{itype.upper()}. Pass the R1 FASTQ as the positional input "
                f"to align from reads, or drop --read2 to correct the "
                f"{itype.upper()} as-is."
            )

    if getattr(args, 'chunked_alignment', False):
        return (
            "--read2 is not supported with --chunked-alignment: the chunked "
            "generator's alignment arrays are single-end/long-read only and "
            "would silently drop R2. For a chunked paired-end COMPASS run use "
            "the dedicated generator instead:\n"
            "    rectify split R1.fastq.gz --read2 R2.fastq.gz --generate-slurm "
            "--short-read \\\n"
            "        -n <n_chunks> -o <chunks_dir> --genome <genome.fa> "
            "--annotation <genes.gff>\n"
            "which emits per-chunk COMPASS align+consensus array scripts."
        )

    return None


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

    _pe_error = _validate_paired_end_args(args)
    if _pe_error:
        print(f"ERROR: {_pe_error}", file=sys.stderr)
        sys.exit(1)

    _proto_error = _validate_protocol_flags(args)
    if _proto_error:
        print(f"ERROR: {_proto_error}", file=sys.stderr)
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
            # 🔴 The chunked generator emits its own `rectify correct` command line
            # and does NOT forward protocol flags (verified: zero occurrences of
            # --ONT-cDNA/--drs/--dT-primed-cDNA in run/chunked_batch.py). A cDNA
            # library would therefore be corrected under the DRS strand rule --
            # wrong terminus for ~half the reads -- with nothing erroring. Refuse
            # rather than emit a script that produces a confident wrong answer.
            _active_proto = [
                flag for flag, dest in _PROTOCOL_FLAGS if getattr(args, dest, False)
            ]
            if _active_proto:
                print(
                    "ERROR: --chunked-alignment does not propagate protocol flags "
                    f"({', '.join(_active_proto)}) into the generated scripts, so the "
                    "data would be corrected under the default (DRS) strand rule.\n"
                    "       Run without --chunked-alignment, or generate the scripts "
                    "and add the flag to the emitted `rectify correct` command by hand.",
                    file=sys.stderr,
                )
                sys.exit(1)
            rc = _generate_chunked_pipeline(args)
            sys.exit(rc)

    # ── run-all --netseq: its own three-stage pipeline, not the align/correct chain ─────────
    # trim -> STAR (absolute floors) -> `rectify netseq`. Dispatched here rather than inside
    # _run_single_sample's Step 0/1/2 chain because it shares NONE of those stages: there is no
    # aligner panel, no consensus and no `correct` arm (planning 835 §2, Kevin 17:40).
    if getattr(args, 'netseq', False):
        if getattr(args, 'manifest', None):
            print("ERROR: --netseq does not support --manifest yet; run one library per "
                  "invocation (the mode writes per-library tracks + spike-in accounting).",
                  file=sys.stderr)
            sys.exit(1)
        _resolve_reference_paths(args)
        from .run.netseq_pipeline import run_netseq_pipeline
        from .run.single_sample import _canonical_sample_id
        from ..align.preprocess import detect_input_type
        _in = Path(str(args.input))
        _out = Path(args.output_dir)
        _work = Path(getattr(args, 'scratch_dir', None) or _out)
        try:
            rc = run_netseq_pipeline(
                args,
                input_path=_in,
                input_type=detect_input_type(_in),
                output_dir=_out,
                work_dir=_work,
                sample_id=_canonical_sample_id(_in),
                genome_path=Path(args.genome) if getattr(args, 'genome', None) else None,
                annotation_path=(Path(args.annotation)
                                 if getattr(args, 'annotation', None) else None),
                threads=getattr(args, 'threads', 4),
            )
        except Exception as exc:
            print(f"ERROR: run-all --netseq failed: {exc}", file=sys.stderr)
            sys.exit(1)
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

    def _aligner_concurrency_arg(value: str) -> str:
        if value == 'auto':
            return value
        try:
            parsed = int(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                "--aligner-concurrency must be 'auto' or a positive integer"
            ) from exc
        if parsed < 1:
            raise argparse.ArgumentTypeError(
                "--aligner-concurrency must be 'auto' or a positive integer"
            )
        return value

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

      # Paired-end short read (Illumina) — full COMPASS panel
      rectify run-all R1.fastq.gz --short-read --read2 R2.fastq.gz \\
      --genome genome.fa --annotation genes.gff -o results/sample/

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
            'Input is short-read data (Illumina/Aviti ≤150 bp). TruSeq-style '
            'RNA-seq uses the COMPASS splice-aware panel: bbmap + STAR×2 + '
            'HISAT2×2 single-end, plus magicblast + gsnap when paired '
            "(--read2). With --dT-primed-cDNA (QuantSeq-class 3'-end "
            'libraries): bbmap + bwa. Either way this replaces the long-read '
            'aligner panel and disables poly(A)-tail trimming, A-tract '
            'correction, and indel modules.'
        ),
    )
    run_parser.add_argument(
        '-2', '--read2',
        dest='read2',
        type=Path,
        default=None,
        help=(
            'Mate-2 (R2) FASTQ for paired-end short-read alignment. Requires '
            '--short-read and a FASTQ positional input (R1). Selects the full '
            'COMPASS panel (bbmap + STAR×2 + HISAT2×2 + magicblast + gsnap) run '
            'paired, with COMPASS-tiebreak consensus. Not supported with '
            '--manifest or --chunked-alignment.'
        ),
    )
    run_parser.add_argument(
        '--read-length',
        dest='read_length',
        type=int,
        default=150,
        help='Read length for STAR sjdbOverhang / index selection (default: 150). '
             'Only used by the paired-end COMPASS panel (--short-read --read2).',
    )

    run_parser.add_argument(
        '--skip-alignment',
        action='store_true',
        help='Skip triple-aligner alignment even for FASTQ input '
             '(use if you already have a multialigned.bam or consensus.bam)'
    )
    run_parser.add_argument(
        '--trust-existing-bams',
        dest='trust_existing_bams',
        action='store_true',
        default=False,
        help=(
            'Reuse pre-existing per-aligner / rectified BAMs without checking '
            'provenance sidecars. Use only when you have manually verified the '
            'BAMs were produced by a compatible rectify version.'
        )
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
    from rectify.data import add_organism_args, add_junction_anchor_args
    add_organism_args(run_parser)
    add_junction_anchor_args(run_parser)

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
             "With --short-read, also selects the 3'-end aligner panel "
             '(bbmap + bwa) instead of the TruSeq COMPASS panel. '
             'By default run-all assumes nanopore direct RNA (poly-A IS in read).'
    )
    run_parser.add_argument(
        '--no-polya-sequenced',
        dest='dT_primed_cDNA',
        action='store_true',
        default=False,
        help=argparse.SUPPRESS,  # Deprecated alias for --dT-primed-cDNA
    )
    run_parser.add_argument(
        '--ONT-cDNA',
        dest='ONT_cDNA',
        action='store_true',
        default=False,
        help='Input is Oxford Nanopore PCR-cDNA (e.g. SQK-PCB114). Reads arrive in '
             'BOTH orientations, so the RNA strand is resolved per read rather than '
             'from the BAM strand. Use this — NOT --dT-primed-cDNA — for ONT PCR-cDNA: '
             '--dT-primed-cDNA applies QuantSeq REV\'s fixed antisense rule and puts '
             'the 3\' end at the wrong terminus for roughly half the reads. '
             'DEFAULT behaviour is Path A: reads are UMI-collapsed to MOLECULES '
             '(trim -> pre-align -> correct-cdna) before the main pipeline, so '
             'counts are molecule counts. See --ont-cdna-path.'
    )
    run_parser.add_argument(
        '--ont-cdna-path',
        dest='ont_cdna_path',
        choices=['a', 'b'],
        default='a',
        help='ONT PCR-cDNA processing path. "a" (DEFAULT) = UMI-collapse to '
             'MOLECULES via correct-cdna, then align+correct; one output row per '
             'molecule. "b" = trim only; one output row per READ, NOT '
             'UMI-deduplicated. PCB114 is PCR-amplified, so a read is not a '
             'molecule and duplicates inflate counts unevenly between libraries: '
             'use "b" only for read-level diagnostics, never for abundance.'
    )
    # ═══════════════════════════════════════════════════════════════════════
    # NET-seq / mNET-seq  (planning 829 / 834 / 835 WP2 / 836)
    # ═══════════════════════════════════════════════════════════════════════
    netseq_group = run_parser.add_argument_group(
        'NET-seq (--netseq)',
        "Nascent Pol II 3'-end pipeline: cutadapt 3' linker trim -> STAR with ABSOLUTE match "
        "floors -> `rectify netseq`. \n"
        "🔴 There is deliberately NO `correct` / consensus stage. The junction-overhang resolver "
        "is a junction-DISCOVERY stage for long reads; the 1-10 nt soft-clip RE-PLACEMENT that "
        "NET-seq needs is `correct`'s logic, and it is implemented for the NET-seq geometry "
        "directly inside `rectify netseq` (donor-side junction-pool rescue against the exon-2 "
        "start, randomer-tolerant remainder, evidence-gated poly(A) walkback -- planning 836). "
        "Running the COMPASS panel plus a `correct` arm on 35-nt reads adds no placement "
        "information and costs the `correct` arm, which is 87x the whole align stage per BAM "
        "(planning 586 / 728).\n"
        "🔴 The 5' randomer is NEVER trimmed. planning 829 §4 measured these libraries as a "
        "MIXTURE (58-60 %% of reads carry no randomer), so `-u 6` would delete six genuine "
        "nucleotides from the majority class. The aligner soft-clips it and "
        "--netseq-umi-length accounts for it per read.",
    )
    netseq_group.add_argument(
        '--netseq',
        dest='netseq',
        action='store_true',
        default=False,
        help="Input is NET-seq / mNET-seq (Pol II-associated nascent RNA). Antisense chemistry: "
             "the RNA 3' end is the read 5' terminus and the gene strand is the INVERSE of the "
             "BAM strand (Churchman convention, --rna3p-at read5p). Selects the three-stage "
             "pipeline described above. For paired-end mNET-seq, isolate read 2 first "
             "(see `rectify netseq-cpa`).",
    )

    trim_help = netseq_group.add_argument
    trim_help(
        '--netseq-adapter',
        dest='netseq_adapter',
        default=CHURCHMAN_3P_LINKER,
        help="3' linker literal for cutadapt. Default is the Churchman 2011 linker, measured in "
             "81.4 %% of PRJNA1521488 reads (planning 829 §2).",
    )
    trim_help(
        '--netseq-min-length',
        dest='netseq_min_length', type=int, default=18,
        help="cutadapt -m: discard reads shorter than this after trimming. Adapter dimers are "
             "<18 nt (8.4 %% of PRJNA1521488).",
    )
    trim_help(
        '--netseq-nextseq-trim',
        dest='netseq_nextseq_trim', type=int, default=20,
        help="cutadapt --nextseq-trim: two-colour-aware quality trim. On NextSeq/NovaSeq a "
             "no-signal base reads as a HIGH-QUALITY G, so plain -q does not remove the poly-G "
             "tail. This is what planning 829 ran on all 52 libraries.",
    )
    trim_help(
        '--netseq-quality-cutoff',
        dest='netseq_quality_cutoff', type=int, default=None,
        help="Use cutadapt -q N instead of --nextseq-trim (four-colour chemistry).",
    )
    trim_help(
        '--netseq-skip-trim',
        dest='netseq_skip_trim', action='store_true', default=False,
        help="Input FASTQ is already linker-trimmed.",
    )

    trim_help(
        '--spikein-fasta',
        dest='spikein_fasta', type=Path, default=None,
        help="Spike-in genome FASTA (e.g. S. pombe ASM294v2 for the 1:10 spike-in of Bryll & "
             "Peterson 2022). Its contigs are prefixed and appended to the primary genome for a "
             "single combined STAR index; `samtools idxstats` then splits the two organisms and "
             "the spike-in fraction is written to <sample>.spikein.tsv. The BAM handed to "
             "`netseq` carries the PRIMARY contigs only. RPM erases a spike-in-normalised claim "
             "by construction, which is why this mode writes counts by default.",
    )
    trim_help(
        '--spikein-prefix',
        dest='spikein_prefix', default='Sp_',
        help="Prefix added to every spike-in contig name in the combined index.",
    )
    trim_help(
        '--star-index',
        dest='star_index', type=Path, default=None,
        help="Prebuilt STAR index directory. Used as-is and never rebuilt.",
    )
    trim_help(
        '--star-index-dir',
        dest='star_index_dir', type=Path, default=None,
        help="Where to build/cache the STAR index when --star-index is not given "
             "(default: <output>/star_index). Reused across runs; the SAindex file is the "
             "completion sentinel.",
    )
    trim_help(
        '--star-path', dest='star_path', default=None, help=argparse.SUPPRESS)
    trim_help(
        '--cutadapt-path', dest='cutadapt_path', default=None, help=argparse.SUPPRESS)
    trim_help(
        '--star-sa-index-nbases', dest='star_sa_index_nbases', type=int, default=11,
        help="STAR --genomeSAindexNbases for the index build. 11 is the manual's value for a "
             "~12 Mb genome; the mammalian default 14 will not build.",
    )

    trim_help(
        '--netseq-umi-length',
        dest='netseq_umi_length', type=int, default=0,
        help="Randomer/UMI length, forwarded to `rectify netseq --umi-length`. It does NOT trim: "
             "it lets the donor-side rescue explain a non-templated clip remainder and strips the "
             "distal N nt before the 3'-tail A-run is measured. 0 (default) = do not assume a "
             "randomer. Use 6 only after CONFIRMING one on the aligned 5'-clip histogram.",
    )
    trim_help(
        '--netseq-dedup',
        dest='netseq_dedup', action='store_true', default=False,
        help="Emit the molecule (UMI-deduplicated) track alongside the read track. Requires "
             "--netseq-umi-length. 🔴 Only valid when the randomer is UNIVERSAL: the overshoot "
             "correction shifts every read with a short clip, so on a MIXED library it moves the "
             "randomer-free class by the full UMI length (planning 829 §4 measured PRJNA1521488 "
             "at 52 %% zero-length clips, so dedup is WRONG for it).",
    )
    trim_help(
        '--netseq-umi-source', dest='netseq_umi_source', default=None,
        choices=['read5p', 'name'], help=argparse.SUPPRESS)
    trim_help(
        '--netseq-exclude-pol3',
        dest='netseq_exclude_pol3', action='store_true', default=False,
        help="Exclude Pol III loci. OFF by default here: the exclusion also drops every `snRNA` "
             "GFF feature, and the U5/SNR7 3' end (planning 829 §4: chrVII:939,521 -, 3.1 %% of "
             "all reads) is the QC observable that proves the strand convention. rDNA stays "
             "excluded either way.",
    )
    trim_help(
        '--netseq-exclude-mito', dest='netseq_exclude_mito', action='store_true', default=False,
        help="Exclude the mitochondrial genome (Pol II does not transcribe it).")
    trim_help(
        '--netseq-rpm',
        dest='netseq_rpm', action='store_true', default=False,
        help="RPM-normalise the tracks. OFF by default: NNLS deconvolution does not conserve "
             "mass so counting must use raw counts, the --netseq-dir refiner is calibrated on "
             "raw signal, and a spike-in-normalised claim needs counts.",
    )
    trim_help(
        '--netseq-output-format', dest='netseq_output_format', nargs='+', default=None,
        choices=['parquet', 'bedgraph', 'bigwig', 'tsv'],
        help="Forwarded to `rectify netseq --output-format`. bedgraph/bigwig STREAM; parquet/tsv "
             "are the read-level table and materialise every record (~23 GB / 50 M reads).",
    )
    trim_help(
        '--netseq-track-position', dest='netseq_track_position', default='corrected',
        choices=['raw', 'corrected'],
        help="Which 3' end drives the primary track. 'corrected' (default) is the post-rescue, "
             "post-walkback end; 'raw' is the bare alignment terminus.",
    )
    trim_help(
        '--netseq-junction-pool', dest='netseq_junction_pool', type=Path, default=None,
        help="External junction table merged with the annotated introns (e.g. long-read-derived).")
    trim_help(
        '--netseq-rescue-min-k', dest='netseq_rescue_min_k', type=int, default=None,
        help="Minimum recovered exon-2 length for a rescue with NO non-templated remainder.")
    trim_help(
        '--netseq-rescue-min-k-with-remainder',
        dest='netseq_rescue_min_k_with_remainder', type=int, default=None,
        help="Minimum recovered exon-2 length when a randomer remainder is invoked (the chance "
             "channel). Default 4 in `rectify netseq`.")
    trim_help(
        '--netseq-rescue-max-intronic', dest='netseq_rescue_max_intronic', type=int, default=None,
        help="How far past the donor a mis-extended aligned 3' end may sit and still be rescued.")
    trim_help(
        '--netseq-pool-include-trna', dest='netseq_pool_include_trna',
        action='store_true', default=False,
        help="Keep tRNA introns in the rescue pool (they are excised by the tRNA endonuclease, "
             "not the spliceosome, and can only manufacture rescues).")
    trim_help(
        '--netseq-pool-include-organellar', dest='netseq_pool_include_organellar',
        action='store_true', default=False,
        help="Keep mitochondrial/organellar introns in the rescue pool. OFF by default: they are "
             "group I/II SELF-SPLICING introns on a genome Pol II does not transcribe, and their "
             "parent feature type is `mRNA` so the tRNA filter misses them. Measured on wt_rep3, "
             "94 of 580 rescues (16%%) were fabricated on chrMito before this filter existed.")
    trim_help(
        '--netseq-walkback-unconditional', dest='netseq_walkback_unconditional',
        action='store_true', default=False,
        help="Walk back through terminal A's over genomic A's even with no non-templated A in the "
             "clip. Correct for poly(A)-selected input, wrong for nascent RNA (it moved 22 %% of "
             "ends and erased the RPL32 exon-1 3'-end peak, 33 reads -> 1).")
    trim_help(
        '--netseq-no-tail-detection', dest='netseq_no_tail_detection',
        action='store_true', default=False, help=argparse.SUPPRESS)
    trim_help(
        '--netseq-no-deconvolution', dest='netseq_no_deconvolution',
        action='store_true', default=False,
        help="Disable NNLS deconvolution (raw + corrected tracks only).")
    trim_help(
        '--netseq-min-atract-length', dest='netseq_min_atract_length', type=int, default=None,
        help=argparse.SUPPRESS)
    trim_help(
        '--netseq-min-mapq', dest='netseq_min_mapq', type=int, default=None,
        help=argparse.SUPPRESS)
    trim_help(
        '--netseq-max-reads', dest='netseq_max_reads', type=int, default=None,
        help="Process at most N reads (smoke tests).")

    run_parser.add_argument(
        '--ont-cdna-poa',
        dest='ont_cdna_poa',
        action='store_true',
        default=False,
        help='Path A only: build an abPOA consensus SEQUENCE per UMI cluster. '
             'Deduplication does not need one and POA is ~44%% of stage-1 runtime, '
             'so it is off by default.'
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
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            'Run the enabled panel aligners in parallel during alignment (default: enabled). '
            'Threads are divided evenly across aligners. '
            'Reduces wall-clock time for the alignment step at the cost of higher peak memory. '
            'Use --no-parallel-aligners to force sequential alignment.'
        )
    )

    run_parser.add_argument(
        '--aligner-concurrency',
        default='auto',
        type=_aligner_concurrency_arg,
        metavar='{auto,1,N}',
        help=(
            'Correct multiple aligners with a shared worker budget. auto disables '
            'this second parallelism axis on small laptops and reserves two CPUs '
            'for merge/main-process overhead on cluster nodes. Use 1 to preserve '
            'sequential-aligner correction behavior.'
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
        choices=['minimap2', 'mapPacBio', 'gapmm2', 'bbmap', 'bwa',
                 'STAR_default', 'STAR_noncanonical', 'HISAT2_default',
                 'HISAT2_noncanonical', 'magicblast', 'gsnap'],
        default=None,
        metavar='ALIGNER',
        dest='base_aligners',
        help=(
            'Base aligners for the consensus pool. When omitted, the default '
            'depends on --short-read: bbmap + bwa for short-read single-end, '
            'the COMPASS panel for short-read paired-end (--read2), '
            'minimap2 for long-read (mapPacBio and gapmm2 are opt-in arms). '
            'Passing this flag overrides the default in every mode. '
            'Example (force short-read panel even without --short-read): '
            '--base-aligners bbmap bwa. The COMPASS members '
            '(STAR_*/HISAT2_*/magicblast/gsnap) require --read2.'
        )
    )

    run_parser.add_argument(
        '--junction-aligners',
        nargs='+',
        choices=['uLTRA', 'deSALT', 'gmap', 'overhang_resolver'],
        default=None,
        metavar='ALIGNER',
        help=(
            'Junction-aware aligners for the consensus pool '
            '(choices: uLTRA, deSALT, gmap, overhang_resolver). uLTRA '
            'requires --annotation; gmap requires a pre-built db (--gmap-db). '
            'overhang_resolver is the Station-A post-pass on the minimap2 '
            'arm; its output SUBSTITUTES that arm downstream (planning/720 '
            'ADOPT verdict). '
            'When omitted, defaults to [] (disabled) under --short-read and '
            '[uLTRA, deSALT, overhang_resolver] under long-read. '
            'Pass --no-junction-aligners to explicitly disable.'
        )
    )

    run_parser.add_argument(
        '--no-junction-aligners',
        dest='junction_aligners',
        action='store_const',
        const=[],
        help=(
            'Disable uLTRA, deSALT and the overhang_resolver post-pass '
            '(use only the --base-aligners set).'
        )
    )

    run_parser.add_argument(
        '--resolver-acceptor-classes',
        choices=['canonical', 'prp18'],
        default='canonical',
        help=(
            "Acceptor candidate classes for the overhang_resolver post-pass. "
            "'prp18' adds the Roy et al. 2023 NAR alternative-3'SS classes "
            "(BG: TG/CG/GG + non-G: AT) — opt-in for splicing missions "
            "(upf1D, prp18D). Default 'canonical'. See planning/722b."
        )
    )

    run_parser.add_argument(
        '--max-intron',
        type=int,
        default=None,
        metavar='BP',
        help=(
            'Maximum intron size for the whole aligner panel. Default: '
            'derived from --annotation as 2x the longest annotated intron '
            '(rounded up to 100, clamped to [1000, 500000]); for the bundled '
            'S. cerevisiae annotation this derives the historical 5000. '
            'Fallback without an annotation: 5000.'
        )
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

    run_parser.add_argument(
        '--gmap-path',
        default='gmap',
        help='Path to gmap executable (used with --junction-aligners gmap)'
    )

    run_parser.add_argument(
        '--gmap-db',
        default=None,
        metavar='DIR',
        help=(
            'Pre-built GMAP database directory (used with --junction-aligners gmap). '
            'Build once with: gmap_build -D <dir> -d <name> <genome.fa>'
        )
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
        help='Retain per-aligner BAMs (one per panel arm) after '
             'consensus selection. By default, per-aligner BAMs are excluded '
             'from the Oak sync to save disk space; only the rectified BAM is kept.'
    )
    bam_group.add_argument(
        '--scratch-dir',
        type=Path,
        default=None,
        dest='scratch_dir',
        metavar='DIR',
        help=(
            'Base directory for intermediate BAM I/O (alignment output, unsorted '
            'corrected BAM, pysam sort temp files). Auto-detected from $SCRATCH, '
            '$SLURM_TMPDIR, or $TMPDIR when running inside an HPC batch job; no '
            'staging overhead on local workstations when no HPC scratch is found. '
            'Durable outputs (corrected_reads.tsv, bedgraphs, report) always stay '
            'in --output-dir.'
        ),
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
        default=True,
        help=argparse.SUPPRESS,  # Deprecated; scratch staging is on by default. Use --scratch-dir to override path.
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
