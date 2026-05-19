"""
Single-sample pipeline runners for ``rectify run-all``.

- ``_run_single_sample`` — top-level entry point for a single FASTQ or BAM
  input. Handles scratch staging, provenance tracking, junction aggregation,
  and DRS Step 0 (poly(A) pre-trim) + Step 4 (poly(A) restore) when ``--drs``
  is set.
- ``_process_one_sample`` — the per-sample worker run inside the
  multi-sample ``ThreadPoolExecutor``. Same Step 0 / align / correct / Step 4
  flow, but parameterised by a manifest row and isolated to its own subdir.
  Imported by ``multi_sample.py``.

Resume-on-restart logic (skip-if-output-exists) is intentionally *kept here*
rather than pushed down into ``stages`` — the decision to reuse a previous
run's output is sample-level, not stage-level.
"""

import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

from .helpers import (
    _rectified_bam_path,
    _validate_bam_integrity,
    _resolve_reference_paths,
)
from .stages import (
    _run_alignment,
    _run_correction,
    _run_correction_per_aligner,
    _run_analysis,
    _run_junction_aggregation,
)


def _process_one_sample(
    sample: Dict[str, str],
    output_dir: Path,
    genome_path: Path,
    annotation_path: Path,
    args,
) -> Tuple[str, int]:
    """
    Process one sample: [DRS trim →] alignment (if needed) → correction [→ DRS restore].
    DRS steps 0 and 4 run only when ``args.drs`` is True and the input is a BAM.
    Returns (sample_id, returncode).  Called from threads.
    """
    import re
    import logging
    _logger = logging.getLogger(__name__)
    sample_id = sample['sample_id']
    safe_id = re.sub(r'[^a-zA-Z0-9_\-.]', '_', sample_id)
    if safe_id != sample_id:
        _logger.warning(f"sample_id '{sample_id}' contains invalid characters, using '{safe_id}'")
        sample_id = safe_id
    input_path = Path(sample['path'] if 'path' in sample else sample['bam_path'])
    sample_output = output_dir / sample_id
    sample_output.mkdir(parents=True, exist_ok=True)

    log_file = sample_output / 'rectify_run.log'

    try:
        with open(log_file, 'w') as log:
            # Determine input type
            from ...align.preprocess import detect_input_type
            input_type = detect_input_type(input_path)

            bam_to_correct = input_path
            drs_metadata_path: Optional[Path] = None

            # ── DRS Step 0: Poly(A)+adapter pre-trimming ────────────────────
            if getattr(args, 'drs', False) and input_type not in ('fastq', 'fastq.gz'):
                import subprocess as _subprocess
                from ..drs_trim_command import trim_drs_bam_polya

                _drs_dir = sample_output / 'drs_trim'
                _drs_dir.mkdir(parents=True, exist_ok=True)
                _trimmed_bam   = _drs_dir / f"{input_path.stem}.bam"
                drs_metadata_path = sample_output / f"{input_path.stem}_polya_trim_metadata.parquet"
                _trimmed_fastq = _drs_dir / f"{input_path.stem}_trimmed.fastq.gz"

                print(f"  [{sample_id}] DRS poly(A) trimming…", flush=True)
                try:
                    _trim_stats = trim_drs_bam_polya(
                        input_bam_path=str(input_path),
                        output_bam_path=str(_trimmed_bam),
                        metadata_path=str(drs_metadata_path),
                        threads=getattr(args, 'threads', 4),
                    )
                    _n_total   = _trim_stats.get('total', 0)
                    _n_trimmed = _trim_stats.get('trimmed', 0)
                    log.write(
                        f"DRS trim: {_n_trimmed:,}/{_n_total:,} reads trimmed\n"
                    )
                    _subprocess.run(
                        [
                            'samtools', 'fastq',
                            '-@', str(max(1, getattr(args, 'threads', 4) - 1)),
                            '-0', str(_trimmed_fastq),
                            str(_trimmed_bam),
                        ],
                        check=True,
                        stdout=_subprocess.DEVNULL,
                        stderr=_subprocess.DEVNULL,
                    )
                    input_path = _trimmed_fastq
                    input_type = 'fastq.gz'
                    log.write(f"DRS trim FASTQ: {_trimmed_fastq}\n")
                except Exception as e:
                    log.write(f"DRS trim failed: {e}\n")
                    return sample_id, 1

            if input_type in ('fastq', 'fastq.gz') and not getattr(args, 'skip_alignment', False):
                print(f"  [{sample_id}] Aligning…", flush=True)
                # --bam-dir: per-sample subdir within the requested bam_dir
                _bam_dir_arg = getattr(args, 'bam_dir', None)
                if _bam_dir_arg:
                    _align_out = Path(_bam_dir_arg) / sample_id
                    _align_out.mkdir(parents=True, exist_ok=True)
                else:
                    _align_out = sample_output

                # Create a per-sample scratch dir for the consensus sort step.
                # This routes pysam.sort I/O to high-bandwidth scratch storage,
                # preventing NFS write hangs when output is on shared Oak NFS.
                from ....slurm import make_job_scratch_dir as _mks_consensus
                import shutil as _shutil_consensus
                _use_scratch = getattr(args, 'use_scratch', True)
                _consensus_ckpt_dir: Optional[str] = None
                _consensus_ckpt_path = None
                if _use_scratch:
                    _consensus_ckpt_path = _mks_consensus(
                        f'rectify_consensus_{sample_id[:12]}'
                    )
                    if _consensus_ckpt_path:
                        _consensus_ckpt_dir = str(_consensus_ckpt_path)

                try:
                    _sample_per_aligner_bams, bam_to_correct = _run_alignment(
                        input_path=input_path,
                        sample_id=sample_id,
                        sample_output_dir=_align_out,
                        genome_path=genome_path,
                        annotation_path=annotation_path,
                        threads=getattr(args, 'threads', 4),
                        parallel_aligners=getattr(args, 'parallel_aligners', False),
                        base_aligners=getattr(args, 'base_aligners', None),
                        junction_aligners=getattr(args, 'junction_aligners', []),
                        chimeric_consensus=getattr(args, 'chimeric_consensus', True),
                        ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
                        desalt_path=getattr(args, 'desalt_path', 'deSALT'),
                        mapPacBio_chunks=getattr(args, 'mapPacBio_chunks', 1),
                        checkpoint_dir=_consensus_ckpt_dir,
                        short_read=getattr(args, 'short_read', False),
                        trust_existing_bams=getattr(args, 'trust_existing_bams', False),
                    )
                    log.write(f"Alignment complete: {bam_to_correct}\n")
                except Exception as e:
                    log.write(f"Alignment failed: {e}\n")
                    return sample_id, 1
                finally:
                    # Clean up consensus scratch dir (batch BAMs already removed
                    # inside run_consensus_selection after copy to Oak)
                    if _consensus_ckpt_path and _consensus_ckpt_path.exists():
                        _shutil_consensus.rmtree(_consensus_ckpt_path, ignore_errors=True)
            elif input_type in ('fastq', 'fastq.gz'):
                # FASTQ but alignment skipped — look for existing rectified.bam
                rectified_bam = _rectified_bam_path(sample_id, sample_output)
                _legacy_bam = sample_output / f"{sample_id}.consensus.bam"
                if not rectified_bam.exists() and _legacy_bam.exists():
                    rectified_bam = _legacy_bam
                if _validate_bam_integrity(rectified_bam):
                    bam_to_correct = rectified_bam
                    log.write(f"Using existing rectified.bam: {rectified_bam}\n")
                else:
                    log.write(
                        f"ERROR: --skip-alignment set but no valid rectified.bam found: {rectified_bam}\n"
                    )
                    return sample_id, 1

            print(f"  [{sample_id}] Correcting 3' ends…", flush=True)
            # Stage correction I/O through $SCRATCH when available to avoid
            # NFS write hangs (position index gzip, large TSV writes under load).
            from ....slurm import make_job_scratch_dir, sync_to_oak as _sync_to_oak
            import shutil as _shutil_corr
            _use_scratch = getattr(args, 'use_scratch', True)
            _corr_scratch = make_job_scratch_dir(f'rectify_{sample_id}') if _use_scratch else None
            _corr_work = _corr_scratch if _corr_scratch else sample_output
            try:
                _run_correction(
                    bam_path=bam_to_correct,
                    output_dir=_corr_work,
                    genome_path=genome_path,
                    annotation_path=annotation_path,
                    args=args,
                )
                if _corr_scratch:
                    _sync_to_oak(_corr_scratch, sample_output)
                log.write("Correction complete\n")
            except Exception as e:
                log.write(f"Correction failed: {e}\n")
                return sample_id, 1
            finally:
                # Always clean up scratch directory, whether sync succeeded or failed
                if _corr_scratch and _corr_scratch.exists():
                    _shutil_corr.rmtree(_corr_scratch, ignore_errors=True)

            # ── DRS Step 4: Restore poly(A)+adapter as soft-clips ───────────
            if (getattr(args, 'drs', False)
                    and getattr(args, 'write_polya_bam', False)
                    and drs_metadata_path and drs_metadata_path.exists()):
                from ..restore_polya_command import restore_polya_softclips

                _corrected_tsv = sample_output / 'corrected_reads.tsv'
                _restored_bam = sample_output / f"{sample_id}.corrected_polya.bam"
                _aligner_bams = {k: str(v) for k, v in _sample_per_aligner_bams.items()
                                 if Path(str(v)).exists()}
                if _corrected_tsv.exists() and _aligner_bams:
                    print(f"  [{sample_id}] DRS poly(A) soft-clip restore…", flush=True)
                    try:
                        _sc_stats = restore_polya_softclips(
                            corrected_tsv_path=str(_corrected_tsv),
                            aligner_bam_paths=_aligner_bams,
                            metadata_path=str(drs_metadata_path),
                            output_bam_path=str(_restored_bam),
                            threads=getattr(args, 'threads', 4),
                        )
                        import pysam as _pysam_sc
                        _restored_tmp = str(_restored_bam) + '.sort_tmp.bam'
                        _pysam_sc.sort('-o', _restored_tmp, str(_restored_bam))
                        _restored_bam.unlink(missing_ok=True)
                        Path(_restored_tmp).rename(_restored_bam)
                        _pysam_sc.index(str(_restored_bam))
                        log.write(f"DRS restore: {_sc_stats.get('restored', 0):,} reads restored\n")
                    except Exception as e:
                        log.write(f"DRS restore failed (non-fatal): {e}\n")
                else:
                    log.write(
                        f"DRS restore skipped — corrected TSV or aligner BAMs not found\n"
                    )

        return sample_id, 0

    except Exception as e:
        with open(log_file, 'a') as log:
            log.write(f"Unexpected error: {e}\n")
        return sample_id, 1
def _run_single_sample(args) -> int:
    """
    Single-sample pipeline:
      Step 0 (DRS BAM + --drs): poly(A)+adapter pre-trimming → trimmed FASTQ
      Step 1 (if FASTQ input): multi-aligner alignment → rectified.bam
      Step 2: correction → corrected_reads.tsv
      Step 3: analysis (no DESeq2 — single sample)
      Step 4 (DRS + --drs): restore trimmed poly(A) as soft-clips
    """
    _resolve_reference_paths(args)

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_path = args.genome
    annotation_path = args.annotation

    # Pre-flight: verify required input files exist before starting any work
    for path, label in [
        (input_path, 'input'),
        (genome_path, 'genome'),
        (annotation_path, 'annotation GFF'),
    ]:
        if path and not Path(path).exists():
            raise FileNotFoundError(f"{label} file not found: {path}")

    print("=" * 70)
    print("RECTIFY: Single-Sample Pipeline")
    print("=" * 70)
    print(f"\nInput:      {input_path}")
    print(f"Output dir: {output_dir}")

    # Determine input type
    from ...align.preprocess import detect_input_type
    input_type = detect_input_type(input_path)
    print(f"Input type: {input_type}")

    bam_to_correct = input_path
    per_aligner_bams: Dict[str, Path] = {}

    # ── Scratch staging setup ────────────────────────────────────────────────
    # If $SCRATCH is available, run all I/O there and rsync back.
    # This avoids NFS contention across concurrent array tasks.
    # The rectified BAM is always copied back to the output dir even on fresh
    # alignments so it survives $SCRATCH's auto-purge and enables job resumption.
    from ....slurm import make_job_scratch_dir, sync_to_oak
    import shutil as _shutil

    use_scratch = getattr(args, 'use_scratch', True)  # default on when scratch available
    scratch_dir = make_job_scratch_dir('rectify_single') if use_scratch else None
    work_dir = scratch_dir if scratch_dir else output_dir

    # --bam-dir: explicit persistent directory for alignment BAMs.
    # --keep-aligner-bams: sync per-aligner BAMs to Oak (default: discard them).
    _bam_dir_arg = getattr(args, 'bam_dir', None)
    bam_dir: Optional[Path] = Path(_bam_dir_arg) if _bam_dir_arg else None
    if bam_dir:
        bam_dir.mkdir(parents=True, exist_ok=True)
    keep_aligner_bams: bool = getattr(args, 'keep_aligner_bams', False)

    if scratch_dir:
        print(f"\nScratch staging enabled: {scratch_dir}")

    # ── Provenance tracker ───────────────────────────────────────────────────
    # scratch_root/oak_root mapping ensures all sidecar JSON keys use canonical
    # Oak paths — so step_is_current() works on re-runs regardless of $SCRATCH.
    from ....utils.provenance import ProvenanceTracker
    tracker = ProvenanceTracker(output_dir=output_dir)
    tracker.set_command(sys.argv)

    import time as _time
    _pipeline_start = _time.perf_counter()

    # ── Step 0 (DRS only): Poly(A)+adapter pre-trimming ─────────────────────
    # When --drs is set and the input is a Dorado-aligned BAM, strip the
    # poly(A) tail and adapter stub before re-alignment.  The trimmed reads
    # are written as an unaligned BAM (preserving the pt:i: tag) and then
    # converted to FASTQ so the alignment step can run normally.
    # The metadata parquet is kept on Oak so Step 4 can restore the soft-clips.
    drs_mode = getattr(args, 'drs', False)
    drs_metadata_path: Optional[Path] = None

    if drs_mode and input_type not in ('fastq', 'fastq.gz'):
        from ..drs_trim_command import trim_drs_bam_polya
        import subprocess as _subprocess

        _sample_stem = input_path.stem
        _drs_dir = work_dir / 'drs_trim'
        _drs_dir.mkdir(parents=True, exist_ok=True)

        _trimmed_bam    = _drs_dir / f"{_sample_stem}.bam"
        drs_metadata_path = output_dir / f"{_sample_stem}_polya_trim_metadata.parquet"
        _trimmed_fastq  = _drs_dir / f"{_sample_stem}_trimmed.fastq.gz"

        print(f"\n[Step 0] DRS poly(A)+adapter pre-trimming...")
        print("-" * 50)
        _t0_drs = _time.perf_counter()

        _trim_stats = trim_drs_bam_polya(
            input_bam_path=str(input_path),
            output_bam_path=str(_trimmed_bam),
            metadata_path=str(drs_metadata_path),
            threads=getattr(args, 'threads', 4),
        )
        _n_total   = _trim_stats.get('total', 0)
        _n_trimmed = _trim_stats.get('trimmed', 0)
        print(
            f"  Trimmed {_n_trimmed:,} / {_n_total:,} reads "
            f"({100 * _n_trimmed / max(_n_total, 1):.1f}%)"
        )

        # Convert trimmed unaligned BAM → FASTQ.
        # Do NOT use -T pt: embedding the pt:i:N tag in the read name causes
        # mapPacBio to write it as a suffix (UUID_pt:i:N), creating duplicate
        # read IDs and breaking parquet lookups. Poly-A lengths are stored in
        # the parquet metadata; the FASTQ read name must be the bare UUID only.
        _subprocess.run(
            [
                'samtools', 'fastq',
                '-@', str(max(1, getattr(args, 'threads', 4) - 1)),
                '-0', str(_trimmed_fastq),
                str(_trimmed_bam),
            ],
            check=True,
        )

        print(f"  FASTQ written: {_trimmed_fastq}")
        print(f"[TIMING] DRS trim: {_time.perf_counter() - _t0_drs:.1f}s")
        tracker.record_step(
            'drs_trim',
            input_files=[input_path],
            output_files=[_trimmed_bam, drs_metadata_path],
        )

        # Re-classify as FASTQ so the alignment step proceeds on trimmed reads.
        input_path = _trimmed_fastq
        input_type = 'fastq.gz'

    # ── Step 0/1: Alignment ───────────────────────────────────────────────────
    if input_type in ('fastq', 'fastq.gz') and not getattr(args, 'skip_alignment', False):
        sample_id = input_path.stem.replace('.fastq', '').replace('.gz', '')
        _align_mode = 'short-read (bbmap + bwa)' if getattr(args, 'short_read', False) else 'long-read consensus'
        print(f"\n[Step 1/3] Aligning ({_align_mode})...")
        print("-" * 50)
        _t0 = _time.perf_counter()
        # Align to bam_dir if specified; otherwise use work_dir (scratch if available)
        align_output_dir = bam_dir if bam_dir else work_dir
        # When --bam-dir forces output to a potentially NFS-backed path, route
        # the consensus sort through scratch to avoid pysam.sort NFS hang.
        _single_ckpt_dir: Optional[str] = None
        if bam_dir and scratch_dir:
            _single_ckpt_dir = str(scratch_dir / 'consensus_ckpt')
        per_aligner_bams, bam_to_correct = _run_alignment(
            input_path=input_path,
            sample_id=sample_id,
            sample_output_dir=align_output_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            threads=getattr(args, 'threads', 4),
            parallel_aligners=getattr(args, 'parallel_aligners', False),
            base_aligners=getattr(args, 'base_aligners', None),
            junction_aligners=getattr(args, 'junction_aligners', []),
            chimeric_consensus=getattr(args, 'chimeric_consensus', True),
            ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
            desalt_path=getattr(args, 'desalt_path', 'deSALT'),
            mapPacBio_chunks=getattr(args, 'mapPacBio_chunks', 1),
            checkpoint_dir=_single_ckpt_dir,
            short_read=getattr(args, 'short_read', False),
            trust_existing_bams=getattr(args, 'trust_existing_bams', False),
        )
        print(f"\nAlignment complete: {bam_to_correct}")
        print(f"[TIMING] Alignment: {_time.perf_counter() - _t0:.1f}s")
        tracker.record_step('align', input_files=[input_path], output_files=[bam_to_correct])
        step_correction = 2
    elif input_type in ('fastq', 'fastq.gz'):
        sample_id = input_path.stem.replace('.fastq', '').replace('.gz', '')
        rectified_bam = _rectified_bam_path(sample_id, output_dir)
        _legacy_bam = output_dir / f"{sample_id}.consensus.bam"
        if not rectified_bam.exists() and _legacy_bam.exists():
            rectified_bam = _legacy_bam
        if _validate_bam_integrity(rectified_bam):
            print(f"\n[Step 1/3] Alignment — using existing rectified.bam")
            if scratch_dir:
                # Stage existing BAM to scratch
                print(f"  Staging to scratch: {scratch_dir}")
                _shutil.copy2(rectified_bam, scratch_dir / rectified_bam.name)
                bai = Path(str(rectified_bam) + '.bai')
                if bai.exists():
                    _shutil.copy2(bai, scratch_dir / bai.name)
                bam_to_correct = scratch_dir / rectified_bam.name
                tracker.register_staged(bam_to_correct, rectified_bam)
            else:
                bam_to_correct = rectified_bam
        else:
            print(
                f"ERROR: --skip-alignment set but no rectified.bam found: {rectified_bam}",
                file=sys.stderr,
            )
            if scratch_dir:
                _shutil.rmtree(scratch_dir, ignore_errors=True)
            return 1
        step_correction = 2
    else:
        print(f"\nBAM input detected — skipping alignment")
        if scratch_dir:
            # Stage BAM to scratch
            print(f"  Staging to scratch: {scratch_dir}")
            _shutil.copy2(input_path, scratch_dir / input_path.name)
            bai = Path(str(input_path) + '.bai')
            if bai.exists():
                _shutil.copy2(bai, scratch_dir / bai.name)
            bam_to_correct = scratch_dir / input_path.name
            tracker.register_staged(bam_to_correct, input_path)
        step_correction = 1

    # ── Step 1/2: Correction ─────────────────────────────────────────────────
    print(f"\n[Step {step_correction}/3] Correcting 3' end positions...")
    print("-" * 50)
    _t0 = _time.perf_counter()
    try:
        if per_aligner_bams:
            # New workflow: correct each aligner's BAM independently, then merge.
            # The final corrected_reads.tsv is selected from post-correction features
            # (five_prime_rescued, confidence, 3' agreement) rather than raw alignment
            # features — which are not cross-comparable across aligners.
            from ...consensus.corrected_consensus import merge_corrected_tsvs, identify_cat5_candidates
            print(f"    Running per-aligner correction ({len(per_aligner_bams)} aligners)...")
            per_aligner_tsvs, per_aligner_corrected_bams = _run_correction_per_aligner(
                per_aligner_bams=per_aligner_bams,
                output_dir=work_dir,
                genome_path=genome_path,
                annotation_path=annotation_path,
                args=args,
            )
            if per_aligner_tsvs:
                _per_aligner_dir = work_dir / 'per_aligner_corrected'
                _summary_tsv = _per_aligner_dir / 'comparison_summary.tsv'
                print(f"    Merging {len(per_aligner_tsvs)} per-aligner TSVs...")
                # Load junction overhang table if provided
                _jot_path = getattr(args, 'junction_overhang_table', None)
                _overhang_table = None
                if _jot_path:
                    from rectify.core.splice.calibrate_junction_overhang import OverhangTable as _OT
                    try:
                        _overhang_table = _OT.from_tsv(_jot_path)
                        print(f"    Junction overhang filter: {_jot_path}")
                    except Exception as _e:
                        print(f"    WARNING: Could not load overhang table {_jot_path}: {_e}",
                              file=sys.stderr)
                # Pass per_aligner_corrected_bams when available to activate
                # HP-edit-distance winner selection (the validated correct-first
                # path per CLAUDE.md PIPELINE ORDER). Falls back to the legacy
                # 5-key sort transparently when the dict is empty.
                corrected_tsv = merge_corrected_tsvs(
                    per_aligner_tsvs=per_aligner_tsvs,
                    output_tsv=work_dir / 'corrected_reads.tsv',
                    summary_tsv=_summary_tsv,
                    per_aligner_corrected_bams={
                        a: str(p) for a, p in per_aligner_corrected_bams.items()
                    } if per_aligner_corrected_bams else None,
                    overhang_table=_overhang_table,
                )
                # Identify Cat5 candidates (reads where aligners contribute unique introns)
                _cat5_tsv = _per_aligner_dir / 'cat5_candidates.tsv'
                identify_cat5_candidates(per_aligner_tsvs, output_tsv=_cat5_tsv)
            else:
                # All per-aligner corrections failed — fall back to consensus BAM
                print(
                    "    WARNING: No per-aligner correction succeeded; "
                    "falling back to consensus BAM correction.",
                    file=sys.stderr,
                )
                corrected_tsv = _run_correction(
                    bam_path=bam_to_correct,
                    output_dir=work_dir,
                    genome_path=genome_path,
                    annotation_path=annotation_path,
                    args=args,
                )
        else:
            # No per-aligner BAMs available (BAM input, or aligner BAMs were discarded).
            # Use existing single-BAM correction path.
            corrected_tsv = _run_correction(
                bam_path=bam_to_correct,
                output_dir=work_dir,
                genome_path=genome_path,
                annotation_path=annotation_path,
                args=args,
            )
    except Exception as e:
        print(f"ERROR in correction step: {e}", file=sys.stderr)
        if scratch_dir:
            _shutil.rmtree(scratch_dir, ignore_errors=True)
        return 1
    print(f"\nCorrection complete: {corrected_tsv}")
    print(f"[TIMING] Correction: {_time.perf_counter() - _t0:.1f}s")
    tracker.record_step('correct', input_files=[bam_to_correct], output_files=[corrected_tsv])

    # ── Early partial sync: corrected TSV → Oak ──────────────────────────────
    # After correction, provenance moves back to Oak. Sync the corrected TSV
    # (and its sidecar) immediately so analysis reads from Oak, not ephemeral
    # scratch — matching the canonical paths written into the sidecar.
    if scratch_dir:
        _oak_corrected_tsv = output_dir / corrected_tsv.name
        _shutil.copy2(corrected_tsv, _oak_corrected_tsv)
        _sidecar = corrected_tsv.parent / f".{corrected_tsv.name}.provenance.json"
        if _sidecar.exists():
            _shutil.copy2(_sidecar, output_dir / _sidecar.name)
        corrected_tsv = _oak_corrected_tsv
        work_dir = output_dir  # analysis outputs go directly to Oak

    # Add sample column (needed by analyze even for single sample)
    try:
        import pandas as pd
        df = pd.read_csv(corrected_tsv, sep='\t')
        sample_id = input_path.stem.replace('.fastq', '').replace('.gz', '').replace('.bam', '')
        df['sample'] = sample_id
        df.to_csv(corrected_tsv, sep='\t', index=False)
    except Exception:
        pass  # Not fatal — analyze can work without sample column

    # ── Step 2/3: Analysis ───────────────────────────────────────────────────
    print(f"\n[Step {step_correction + 1}/3] Analyzing results (single-sample)...")
    print("-" * 50)
    print("Note: DESeq2 requires multiple samples — skipped for single-sample run.")
    _t0 = _time.perf_counter()
    try:
        _run_analysis(
            corrected_tsv=corrected_tsv,
            output_dir=work_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            args=args,
            n_samples=1,
        )
    except Exception as e:
        print(f"ERROR in analysis step: {e}", file=sys.stderr)
        # Not fatal — correction output is still usable

    # ── Junction aggregation ─────────────────────────────────────────────────
    if genome_path and annotation_path:
        _run_junction_aggregation(
            bam_path=bam_to_correct,
            genome_path=genome_path,
            annotation_path=annotation_path,
            output_dir=work_dir,
            config=vars(args),
        )

    print(f"[TIMING] Analysis:   {_time.perf_counter() - _t0:.1f}s")

    # ── Copy from scratch to Oak ─────────────────────────────────────────────
    if scratch_dir:
        print(f"\nCopying outputs to Oak: {output_dir}")
        _t0 = _time.perf_counter()
        # When bam_dir is set, BAMs went directly there (not scratch), so exclude_aligner_bams
        # only matters when BAMs are on scratch. When keep_aligner_bams is False (default),
        # skip per-aligner BAMs in the scratch→Oak sync to save disk space.
        _exclude_aligner = (not keep_aligner_bams) and (bam_dir is None)
        sync_to_oak(scratch_dir, output_dir, exclude_aligner_bams=_exclude_aligner)
        print(f"[TIMING] Sync to Oak: {_time.perf_counter() - _t0:.1f}s")
        _shutil.rmtree(scratch_dir, ignore_errors=True)

    # ── Step 4 (DRS only, opt-in): Restore poly(A) from winning aligner's raw BAM ─
    # Pulls each read from the winning aligner's raw (pre-correction) BAM and
    # re-attaches the original Dorado poly(A) tail from parquet as a 3' soft clip.
    # The result (corrected_polya.bam) shows the full read as Dorado called it;
    # compare against corrected.bam in IGV to see exactly what Rectify changed.
    if (drs_mode and getattr(args, 'write_polya_bam', False)
            and drs_metadata_path and drs_metadata_path.exists()):
        from ..restore_polya_command import restore_polya_softclips
        _orig_stem = Path(args.input).stem
        _restored_bam = output_dir / f"{_orig_stem}.corrected_polya.bam"
        _aligner_bams_for_restore = {k: str(v) for k, v in per_aligner_bams.items()
                                     if Path(str(v)).exists()}
        if corrected_tsv and Path(corrected_tsv).exists() and _aligner_bams_for_restore:
            print(f"\n[Step 4] DRS poly(A) restoration from winning aligner raw BAMs...")
            print("-" * 50)
            _t0_sc = _time.perf_counter()
            _sc_stats = restore_polya_softclips(
                corrected_tsv_path=str(corrected_tsv),
                aligner_bam_paths=_aligner_bams_for_restore,
                metadata_path=str(drs_metadata_path),
                output_bam_path=str(_restored_bam),
                threads=getattr(args, 'threads', 4),
            )
            import pysam as _pysam_sc
            _restored_tmp = str(_restored_bam) + '.sort_tmp.bam'
            _pysam_sc.sort('-o', _restored_tmp, str(_restored_bam))
            _restored_bam.unlink(missing_ok=True)
            Path(_restored_tmp).rename(_restored_bam)
            _pysam_sc.index(str(_restored_bam))
            print(f"  Restored: {_sc_stats.get('restored', 0):,} reads")
            print(f"[TIMING] DRS restore: {_time.perf_counter() - _t0_sc:.1f}s")
            tracker.record_step(
                'drs_restore',
                input_files=list(_aligner_bams_for_restore.values()) + [str(drs_metadata_path)],
                output_files=[_restored_bam],
            )
        else:
            print(
                f"\n[Step 4] DRS poly(A) restoration skipped "
                f"— corrected TSV or aligner BAMs not found",
                file=sys.stderr,
            )

    # ── Provenance manifest ───────────────────────────────────────────────────
    tracker.save()
    manifest_path = tracker.prov_file

    # ── Summary ──────────────────────────────────────────────────────────────
    _pipeline_elapsed = _time.perf_counter() - _pipeline_start
    print("\n" + "=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"\nOutput directory: {output_dir}")
    print(f"  Corrected 3' ends: {corrected_tsv}")
    print(f"  Clusters:          {output_dir / 'cpa_clusters.tsv'}")
    print(f"  HTML report:       {output_dir / 'report.html'}")
    print(f"  Provenance:        {manifest_path}")
    print(f"\n[TIMING] Total pipeline: {_pipeline_elapsed:.1f}s ({_pipeline_elapsed/60:.1f} min)")

    return 0
