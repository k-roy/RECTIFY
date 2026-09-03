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
    _emit_legacy_consensus_symlink,
    _final_rectified_bam_path,
    _multialigned_bam_path,
    _validate_bam_integrity,
    _resolve_reference_paths,
)
from .stages import (
    _run_alignment,
    _run_browser_pack,
    _run_correction,
    _run_correction_per_aligner,
    _run_analysis,
    _run_junction_aggregation,
    _run_station_bc,
)


def _run_ont_cdna_path_a(
    args,
    input_path: Path,
    input_type: str,
    work_dir: Path,
    sample_output: Path,
    sample_id: str,
    log=None,
) -> Tuple[Path, str]:
    """ONT PCR-cDNA **Path A**: UMI-collapse to MOLECULES before the main pipeline.

    Returns ``(new_input_path, new_input_type)`` — a FASTQ of one consensus record
    per molecule, which the caller then feeds to the normal align → correct flow.

    🔴 WHY THIS IS THE DEFAULT (Kevin, standing policy 2026-08-03).
    *"cDNA analyses should always use umi deduped reads. DRS needs no such
    treatment."* SQK-PCB114 is PCR-amplified and carries a 27-nt UMI, so a READ is
    NOT a MOLECULE: duplicates inflate counts and do so **unevenly between
    libraries**, which is why the inflation does not cancel in a ratio. Path B
    (trim → align → correct) emits pre-collapse reads and cannot satisfy that
    policy no matter how well the rest of the pipeline behaves.

    THE CHAIN (Kevin's spec, 2026-08-04) — **NO PRE-TRIM**:
      1. a pre-alignment — ``correct-cdna`` derives its clustering anchor from
         alignment coordinates and skips unmapped reads
         (``cdna/read_info.py::extract_read_info``), so it requires an ALIGNED BAM.
         The reads go in **INTACT**: adapter, UMI and poly-A still attached.
      2. ``correct-cdna`` — groups reads into UMI molecules on their 5'/3' ends with
         window tolerance plus UMI edit distance, collapses each group to one
         consensus molecule, records the poly-A tail length, then ``pretrim_consensus``
         (``cdna/io.py:168``) strips the 3' tail+adapter and 5' UMI+adapter from the
         CONSENSUS, keeping the strip lengths as ``XQ``/``XK``.
      3. the caller re-aligns the collapsed mRNA body and corrects it as usual.

    🔴 **The trim happens AFTER collapse, never before.** An earlier version of this
    function ran ``trim-cdna-polya`` first; that was the defect. ``correct-cdna``
    reads the read SEQUENCE (``extract_read_info`` has ZERO ``get_tag`` calls), so
    pre-trimming removes the very structure it needs — see the comment at step 1.

    ``-y`` is required at every alignment so any FASTQ-comment tags survive into the
    BAM; ``stages._run_alignment`` handles that for the caller's alignment, and the
    pre-alignment here passes it explicitly.

    Raises on failure — a silently-skipped collapse would emit read-level counts
    under a molecule-level contract, which is precisely the defect this prevents.
    """
    import subprocess as _subprocess

    def _log(msg: str) -> None:
        print(f"  [{sample_id}] {msg}", flush=True)
        if log is not None:
            log.write(msg + "\n")

    if input_type not in ('fastq', 'fastq.gz'):
        raise ValueError(
            "ONT-cDNA Path A needs FASTQ input; got %r. Convert the BAM to FASTQ "
            "first, or pass an already-collapsed BAM with --ont-cdna-path b."
            % input_type
        )

    cdna_dir = work_dir / 'cdna_path_a'
    cdna_dir.mkdir(parents=True, exist_ok=True)

    base = input_path.name
    for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break

    # ---- 1. pre-alignment (correct-cdna needs alignment anchors) ------------
    # 🔴 NO PRE-TRIM HERE. `correct-cdna` reads the read SEQUENCE, not tags --
    # `cdna/read_info.py::extract_read_info` contains ZERO `get_tag` calls -- and it
    # derives the UMI, the orientation, the XF full-length tier and the tail length
    # from the adapter/UMI/poly-A structure still attached to the read. Trimming
    # first removes exactly what it needs. Measured (planning/567):
    #   detect_full_length_tier 2 -> 0 on every molecule (XF destroyed, and XF is
    #     the full-length gate every 3'-end analysis depends on)
    #   ~9% of reads (Type-2, SSP-less) hit read_info.py `return None` and are
    #     silently DROPPED, exit 0
    #   the survivors are a BIASED sample -- the trim fires only on reads with a
    #     detectable tail, i.e. selectively on the highest-quality reads
    # Smoke on 9,976 untrimmed reads: Type-2 = 13.1% (expected ~9%, NOT ~0%) and
    # XF=2 = 73.5% -- both correct only when the reads arrive intact.
    # The trim belongs AFTER collapse, and already happens there:
    # `pretrim_consensus` in cdna/io.py:168 strips adapter/UMI/poly-A from the
    # consensus and records the strip lengths as XQ/XK.
    pre_bam = cdna_dir / f"{base}_pre.bam"
    if pre_bam.exists() and pre_bam.stat().st_size > 0:
        _log(f"ONT-cDNA Path A: reusing existing pre-alignment {pre_bam.name}")
    else:
        _log("ONT-cDNA Path A step 1/2: pre-alignment for UMI anchors…")
        threads = str(getattr(args, 'threads', 4))
        mm2 = _subprocess.Popen(
            ['minimap2', '-y', '-t', threads, '-ax', 'splice', '--secondary=no',
             str(args.genome), str(input_path)],
            stdout=_subprocess.PIPE, stderr=_subprocess.DEVNULL,
        )
        sort = _subprocess.Popen(
            ['samtools', 'sort', '-@', '2', '-o', str(pre_bam)],
            stdin=mm2.stdout, stderr=_subprocess.DEVNULL,
        )
        mm2.stdout.close()
        sort.communicate()
        if sort.returncode != 0 or mm2.wait() != 0:
            raise RuntimeError("ONT-cDNA Path A pre-alignment failed")
        _subprocess.run(['samtools', 'index', str(pre_bam)], check=True)

    # ---- 2. UMI collapse: reads -> molecules --------------------------------
    stage1_dir = cdna_dir / 'stage1'
    collapsed = stage1_dir / 'stage1_consensus.fastq.gz'
    if collapsed.exists() and collapsed.stat().st_size > 0:
        _log(f"ONT-cDNA Path A: reusing existing UMI collapse {collapsed.name}")
    else:
        _log("ONT-cDNA Path A step 2/2: UMI clustering (reads → molecules)…")
        stage1_dir.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, '-m', 'rectify', 'correct-cdna', str(pre_bam),
               '-o', str(stage1_dir), '--reference', str(args.genome)]
        if getattr(args, 'annotation', None):
            cmd += ['--gff', str(args.annotation)]
        # Dedup does not need a consensus SEQUENCE, and POA is ~44% of stage-1
        # runtime, so skip it unless the caller explicitly wants consensus bases.
        if not getattr(args, 'ont_cdna_poa', False):
            cmd.append('--no-poa')
        rc = _subprocess.run(cmd).returncode
        if rc != 0:
            raise RuntimeError(f"rectify correct-cdna failed (rc={rc})")
    if not collapsed.exists() or collapsed.stat().st_size == 0:
        raise RuntimeError(
            "ONT-cDNA Path A produced no collapsed FASTQ at %s" % collapsed)

    _log(f"ONT-cDNA Path A complete → {collapsed} (one record per MOLECULE)")
    return collapsed, 'fastq.gz'


def _run_ont_cdna_prepare(
    args,
    input_path: Path,
    input_type: str,
    work_dir: Path,
    sample_output: Path,
    sample_id: str,
    log=None,
) -> Tuple[Path, str]:
    """Dispatch ONT PCR-cDNA preparation: Path A (default) or Path B.

    🔴 CALLED FROM BOTH ``_run_single_sample`` AND ``_process_one_sample`` — that is
    the point of this function. The previous implementation lived inline in
    ``_process_one_sample`` only, so ``run-all --ONT-cDNA reads.fastq`` (positional
    input, the command a lab member would naturally type) silently skipped cDNA
    preparation entirely and resolved every read by annotated-gene overlap: ~19% of
    3' ends land >250nt inside the CDS that way, versus ~1.7-2.2% for tag-resolved
    reads. Nothing errored. **Do not re-inline this; add callers instead.**

    Path A (``--ont-cdna-path a``, DEFAULT) → one record per MOLECULE (UMI-collapsed).
    Path B (``--ont-cdna-path b``)         → one record per READ (trim only).
    """
    mode = str(getattr(args, 'ont_cdna_path', 'a') or 'a').lower()

    if mode == 'a':
        return _run_ont_cdna_path_a(
            args, input_path, input_type, work_dir, sample_output, sample_id, log=log)

    # ---- Path B: trim only; output rows are pre-collapse READS --------------
    from ..cdna_trim_command import trim_cdna_fastq_polya

    msg = ("ONT-cDNA Path B selected: output rows are pre-collapse READS, not "
           "molecules. Counts/abundance from this output are NOT UMI-deduplicated.")
    print(f"  [{sample_id}] WARNING: {msg}", flush=True)
    if log is not None:
        log.write(msg + "\n")

    cdna_dir = work_dir / 'cdna_trim'
    cdna_dir.mkdir(parents=True, exist_ok=True)
    base = input_path.name
    for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    cdna_fastq = cdna_dir / f"{base}_cdna_trimmed.fastq.gz"
    cdna_meta = sample_output / f"{base}_cdna_trim_metadata.tsv"

    print(f"  [{sample_id}] ONT cDNA poly(A)/poly(T) trimming…", flush=True)
    stats = trim_cdna_fastq_polya(
        input_fastq_path=str(input_path),
        output_fastq_path=str(cdna_fastq),
        metadata_path=str(cdna_meta),
        trim_5p_polyt=True,
    )
    if log is not None:
        log.write(f"ONT cDNA trim stats: {stats}\n")
        log.write(f"ONT cDNA trim FASTQ: {cdna_fastq}\n")
    return cdna_fastq, 'fastq.gz'


def _run_samtools_fastq(input_bam: Path, output_fastq: Path, *, threads: int) -> None:
    import subprocess as _subprocess

    result = _subprocess.run(
        [
            'samtools', 'fastq',
            '-@', str(max(1, threads - 1)),
            '-0', str(output_fastq),
            str(input_bam),
        ],
        check=False,
        stdout=_subprocess.PIPE,
        stderr=_subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        stderr = (result.stderr or '').strip()
        stdout = (result.stdout or '').strip()
        detail = stderr or stdout or 'no samtools output captured'
        raise RuntimeError(f"samtools fastq failed (exit {result.returncode}): {detail}")


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
    import shutil as _shutil
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

    from ....slurm import sync_to_oak as _sync_to_oak
    _scratch_arg = getattr(args, 'scratch_dir', None)
    # Full sample_id (no truncation) is required for concurrent worker isolation:
    # two samples sharing a long prefix (e.g. by4742_drs_wt_rep1 / *_rep2) would
    # otherwise hash to the same scratch dir and clobber each other's BAMs.
    if _scratch_arg is not None:
        from ....slurm import resolve_scratch_dir
        _work: Path = resolve_scratch_dir(f'rectify_{sample_id}', base_dir=_scratch_arg)
        _using_scratch = True
    else:
        from ....slurm import make_job_scratch_dir
        _scratch = make_job_scratch_dir(f'rectify_{sample_id}')
        if _scratch is not None:
            _work = _scratch
            _using_scratch = True
        else:
            _work = sample_output
            _using_scratch = False

    _sample_per_aligner_bams: Dict[str, Path] = {}

    try:
        with open(log_file, 'w') as log:
            # Determine input type
            from ...align.preprocess import detect_input_type
            input_type = detect_input_type(input_path)

            bam_to_correct = input_path
            drs_metadata_path: Optional[Path] = None

            # ── DRS Step 0: Poly(A)+adapter pre-trimming ────────────────────
            if getattr(args, 'drs', False) and input_type not in ('fastq', 'fastq.gz'):
                from ..drs_trim_command import trim_drs_bam_polya

                _drs_dir = _work / 'drs_trim'
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
                    _run_samtools_fastq(
                        _trimmed_bam,
                        _trimmed_fastq,
                        threads=getattr(args, 'threads', 4),
                    )
                    input_path = _trimmed_fastq
                    input_type = 'fastq.gz'
                    log.write(f"DRS trim FASTQ: {_trimmed_fastq}\n")
                except Exception as e:
                    log.write(f"DRS trim failed: {e}\n")
                    return sample_id, 1

            # ── ONT PCR-cDNA Step 0: poly(A)/poly(T) + adapter pre-trimming ──
            # This stage is what makes --ONT-cDNA *work*, not merely parse.
            # correct --ONT-cDNA resolves the RNA strand per read with the
            # precedence XO/XY tag -> `ro` tag -> annotated-gene overlap. The
            # `ro` label is written ONLY here, and `-y` (already passed by
            # core/align/multi_aligner.py) carries it through alignment into the
            # BAM. Without this stage every read falls through to the annotated
            # -gene channel, which is by far the weakest: ~19% of its 3' ends
            # land >250 nt inside the CDS versus ~1.7-2.2% for the two tailed
            # classes -- and nothing errors, so the degradation is silent.
            # trim_5p_polyt=True is required: it is what labels the ANTISENSE
            # reads (~44% of a PCB114 library), which is the whole point.
            if getattr(args, 'ONT_cDNA', False) and input_type in ('fastq', 'fastq.gz'):
                try:
                    input_path, input_type = _run_ont_cdna_prepare(
                        args, input_path, input_type, _work, sample_output,
                        sample_id, log=log,
                    )
                except Exception as e:
                    log.write(f"ONT cDNA preparation failed: {e}\n")
                    print(f"ERROR: [{sample_id}] ONT cDNA preparation failed: {e}",
                          file=sys.stderr)
                    return sample_id, 1

            if input_type in ('fastq', 'fastq.gz') and not getattr(args, 'skip_alignment', False):
                print(f"  [{sample_id}] Aligning…", flush=True)
                # --bam-dir: per-sample subdir within the requested bam_dir
                _bam_dir_arg = getattr(args, 'bam_dir', None)
                if _bam_dir_arg:
                    _align_out = Path(_bam_dir_arg) / sample_id
                    _align_out.mkdir(parents=True, exist_ok=True)
                else:
                    _align_out = _work

                _consensus_ckpt_dir: Optional[str] = str(_work / 'consensus_ckpt')

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
                        junction_aligners=getattr(args, 'junction_aligners', None),
                        chimeric_consensus=getattr(args, 'chimeric_consensus', True),
                        ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
                        desalt_path=getattr(args, 'desalt_path', 'deSALT'),
                        gmap_path=getattr(args, 'gmap_path', 'gmap'),
                        gmap_db=getattr(args, 'gmap_db', None),
                        mapPacBio_chunks=getattr(args, 'mapPacBio_chunks', 1),
                        checkpoint_dir=_consensus_ckpt_dir,
                        short_read=getattr(args, 'short_read', False),
                        trust_existing_bams=getattr(args, 'trust_existing_bams', False),
                        read2=getattr(args, 'read2', None),
                        dt_primed_cdna=getattr(args, 'dT_primed_cDNA', False),
                        read_length=getattr(args, 'read_length', 150),
                        max_intron=getattr(args, 'max_intron', None),
                        resolver_acceptor_classes=getattr(
                            args, 'resolver_acceptor_classes', 'canonical'),
                    )
                    log.write(f"Alignment complete: {bam_to_correct}\n")
                except Exception as e:
                    log.write(f"Alignment failed: {e}\n")
                    return sample_id, 1
            elif input_type in ('fastq', 'fastq.gz'):
                # FASTQ but alignment skipped — look for existing multialigned.bam
                # (the pre-correction merged alignment artifact).
                multialigned_bam = _multialigned_bam_path(sample_id, sample_output)
                _legacy_bam = sample_output / f"{sample_id}.consensus.bam"
                if not multialigned_bam.exists() and _legacy_bam.exists():
                    multialigned_bam = _legacy_bam
                if _validate_bam_integrity(multialigned_bam):
                    bam_to_correct = multialigned_bam
                    log.write(f"Using existing multialigned.bam: {multialigned_bam}\n")
                else:
                    log.write(
                        f"ERROR: --skip-alignment set but no valid multialigned.bam found: {multialigned_bam}\n"
                    )
                    return sample_id, 1

            print(f"  [{sample_id}] Correcting 3' ends…", flush=True)
            try:
                if _sample_per_aligner_bams:
                    from ...consensus.corrected_consensus import (
                        merge_corrected_tsvs,
                        identify_cat5_candidates,
                        write_corrected_consensus_bam,
                        _stage_raw_bams,
                    )
                    per_aligner_tsvs, per_aligner_corrected_bams = _run_correction_per_aligner(
                        per_aligner_bams=_sample_per_aligner_bams,
                        output_dir=_work,
                        genome_path=genome_path,
                        annotation_path=annotation_path,
                        args=args,
                    )
                    if per_aligner_tsvs:
                        _per_aligner_dir = _work / 'per_aligner_corrected'
                        _summary_tsv = _per_aligner_dir / 'comparison_summary.tsv'
                        _jot_path = getattr(args, 'junction_overhang_table', None)
                        _overhang_table = None
                        if _jot_path:
                            from rectify.core.splice.calibrate_junction_overhang import OverhangTable as _OT
                            try:
                                _overhang_table = _OT.from_tsv(_jot_path)
                            except Exception as _e:
                                print(
                                    f"  [{sample_id}] WARNING: Could not load overhang table"
                                    f" {_jot_path}: {_e}",
                                    file=sys.stderr,
                                )
                        from rectify.utils.genome import load_genome as _load_genome_for_merge
                        _merge_genome = _load_genome_for_merge(str(genome_path))
                        _raw_bams_s1 = {
                            a: str(p) for a, p in _sample_per_aligner_bams.items()
                        } if _sample_per_aligner_bams else {}
                        from rectify.data import resolve_min_junction_anchor_bp as _resolve_anchor
                        with _stage_raw_bams(_raw_bams_s1) as _staged_s1:
                            merge_corrected_tsvs(
                                per_aligner_tsvs=per_aligner_tsvs,
                                output_tsv=_work / 'corrected_reads.tsv',
                                summary_tsv=_summary_tsv,
                                per_aligner_corrected_bams={
                                    a: str(p) for a, p in per_aligner_corrected_bams.items()
                                } if per_aligner_corrected_bams else None,
                                per_aligner_raw_bams=_staged_s1 if not per_aligner_corrected_bams else None,
                                genome=_merge_genome,
                                overhang_table=_overhang_table,
                                lazy_scoring_workers=max(1, int(getattr(args, 'threads', 1) or 1)),
                                min_junction_anchor_bp=_resolve_anchor(args),
                            )
                            _cat5_tsv = _per_aligner_dir / 'cat5_candidates.tsv'
                            identify_cat5_candidates(per_aligner_tsvs, output_tsv=_cat5_tsv)
                            # Emit the final corrected (rectified) BAM — symmetric with
                            # the chunked path.  Named <sample>.rectified.bam (the actually-
                            # rectified product); a legacy corrected_consensus.bam symlink is
                            # dropped for back-compat.  Non-fatal: merged TSV is the primary
                            # output; missing BAM is recoverable via post-hoc
                            # write_corrected_consensus_bam.
                            _consensus_bam_out = _final_rectified_bam_path(sample_id, _work)
                            try:
                                _consensus_stats = write_corrected_consensus_bam(
                                    per_aligner_raw_bams=_staged_s1,
                                    per_aligner_tsvs=per_aligner_tsvs,
                                    merged_tsv=_work / 'corrected_reads.tsv',
                                    output_bam=_consensus_bam_out,
                                    genome=_merge_genome,
                                    threads=max(1, int(getattr(args, 'threads', 1) or 1)),
                                    strict=True,
                                )
                                _emit_legacy_consensus_symlink(_consensus_bam_out)
                                print(
                                    f"  [{sample_id}] {_consensus_bam_out.name} written"
                                    f" ({_consensus_stats.get('written', '?')} reads)",
                                    flush=True,
                                )
                            except Exception as _cbam_exc:
                                print(
                                    f"  [{sample_id}] WARNING: {_consensus_bam_out.name} write failed"
                                    f" (non-fatal): {_cbam_exc}",
                                    file=sys.stderr,
                                )
                    else:
                        print(
                            f"  [{sample_id}] WARNING: No per-aligner correction succeeded; "
                            "falling back to consensus BAM correction.",
                            file=sys.stderr,
                        )
                        _run_correction(
                            bam_path=bam_to_correct,
                            output_dir=_work,
                            genome_path=genome_path,
                            annotation_path=annotation_path,
                            args=args,
                        )
                else:
                    _run_correction(
                        bam_path=bam_to_correct,
                        output_dir=_work,
                        genome_path=genome_path,
                        annotation_path=annotation_path,
                        args=args,
                    )
                # Sync corrected TSV immediately so DRS Step 4 reads from NFS path.
                if _using_scratch:
                    _tsv_src = _work / 'corrected_reads.tsv'
                    if _tsv_src.exists():
                        _shutil.copy2(_tsv_src, sample_output / 'corrected_reads.tsv')
                        _sidecar = _work / '.corrected_reads.tsv.provenance.json'
                        if _sidecar.exists():
                            _shutil.copy2(_sidecar, sample_output / _sidecar.name)
                    _cbam_src = _final_rectified_bam_path(sample_id, _work)
                    if _cbam_src.exists():
                        _cbam_dst = _final_rectified_bam_path(sample_id, sample_output)
                        _shutil.copy2(_cbam_src, _cbam_dst)
                        _cbam_bai = Path(str(_cbam_src) + '.bai')
                        if _cbam_bai.exists():
                            _shutil.copy2(_cbam_bai, Path(str(_cbam_dst) + '.bai'))
                        _emit_legacy_consensus_symlink(_cbam_dst)
                log.write("Correction complete\n")
            except Exception as e:
                log.write(f"Correction failed: {e}\n")
                return sample_id, 1

            # ── DRS Step 4: Restore poly(A)+adapter as soft-clips ───────────
            # Per-aligner BAMs are still on scratch (_align_out or _work) here.
            if (getattr(args, 'drs', False)
                    and getattr(args, 'write_polya_bam', False)
                    and drs_metadata_path and drs_metadata_path.exists()):
                from ..restore_polya_command import restore_polya_softclips

                _corrected_tsv = sample_output / 'corrected_reads.tsv'
                # Step 4 BAM lands on scratch (or sample_output when no scratch);
                # the final sync_to_oak ships it to NFS along with everything else.
                # Writing + pysam-sort here on scratch avoids the NFS pysam-sort hang.
                _restored_bam = _work / f"{sample_id}.corrected_polya.bam"
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
                        Path(_restored_tmp).replace(_restored_bam)
                        _pysam_sc.index(str(_restored_bam))
                        log.write(f"DRS restore: {_sc_stats.get('restored', 0):,} reads restored\n")
                    except Exception as e:
                        log.write(f"DRS restore failed (non-fatal): {e}\n")
                else:
                    log.write(
                        f"DRS restore skipped — corrected TSV or aligner BAMs not found\n"
                    )

            # ── Final sync: scratch → NFS ────────────────────────────────────
            if _using_scratch:
                _keep_aligner = getattr(args, 'keep_aligner_bams', False)
                _bam_dir_final = getattr(args, 'bam_dir', None)
                _exclude_aligner = (not _keep_aligner) and (_bam_dir_final is None)
                _sync_to_oak(_work, sample_output, exclude_aligner_bams=_exclude_aligner)
                log.write("Sync to NFS complete\n")

        return sample_id, 0

    except Exception as e:
        with open(log_file, 'a') as log:
            log.write(f"Unexpected error: {e}\n")
        return sample_id, 1
    finally:
        if _using_scratch and _work.exists():
            _shutil.rmtree(_work, ignore_errors=True)
#: Read-file suffixes stripped to derive a sample id, longest first.
_READ_SUFFIXES = ('.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz',
                  '.fastq', '.fq', '.fasta', '.fa',
                  '.bam', '.sam', '.cram')


def _canonical_sample_id(input_path) -> str:
    """Sample id from the ORIGINAL run-all input, before any Step-0 rewrite.

    Pure, so the naming contract is unit-testable without running a pipeline. Equivalent to the
    old ``stem.replace('.fastq','').replace('.gz','')`` for every real input; the point is WHERE
    it is called, not what it computes -- see the comment at the call site (planning 831).
    """
    name = Path(input_path).name
    lowered = name.lower()
    for sfx in _READ_SUFFIXES:
        if lowered.endswith(sfx):
            return name[: -len(sfx)]
    return Path(input_path).stem


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

    # 🔴 Canonical sample id — derived from the ORIGINAL input, ONCE, BEFORE any Step-0 rewrite.
    # Step 0 REPLACES ``input_path`` with a derived file: ``<sample>_trimmed.fastq.gz`` under
    # ``--drs``, ``stage1/stage1_consensus.fastq.gz`` under ``--ONT-cDNA`` Path A. Every downstream
    # site used to re-derive the id from ``input_path.stem``, so a positional
    # ``run-all --ONT-cDNA reads.fastq`` renamed the sample to **stage1_consensus** — every output
    # file and the ``sample`` column of ``corrected_reads.tsv`` (planning 831), which silently
    # breaks any multi-sample join keyed on the sample name.
    sample_id = _canonical_sample_id(input_path)
    print(f"Sample id:  {sample_id}")

    bam_to_correct = input_path
    per_aligner_bams: Dict[str, Path] = {}

    # ── Scratch staging setup ────────────────────────────────────────────────
    # If $SCRATCH is available, run all I/O there and rsync back.
    # This avoids NFS contention across concurrent array tasks.
    # The rectified BAM is always copied back to the output dir even on fresh
    # alignments so it survives $SCRATCH's auto-purge and enables job resumption.
    from ....slurm import sync_to_oak
    import shutil as _shutil

    _scratch_arg = getattr(args, 'scratch_dir', None)
    if _scratch_arg is not None:
        from ....slurm import resolve_scratch_dir
        scratch_dir: Optional[Path] = resolve_scratch_dir('rectify_single', base_dir=_scratch_arg)
    else:
        from ....slurm import make_job_scratch_dir
        scratch_dir = make_job_scratch_dir('rectify_single')
    work_dir: Path = scratch_dir if scratch_dir else output_dir

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
        _run_samtools_fastq(
            _trimmed_bam,
            _trimmed_fastq,
            threads=getattr(args, 'threads', 4),
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

    # ── Step 0 (ONT PCR-cDNA): UMI collapse (Path A) or trim-only (Path B) ───
    # 🔴 THIS CALL WAS MISSING. The cDNA preparation existed only in
    # `_process_one_sample` (the --manifest worker), so a positional invocation
    # -- `run-all --ONT-cDNA reads.fastq`, the natural one -- parsed the flag,
    # ran, exited 0, and silently corrected every read through the weakest
    # strand channel. Both entry points now call the SAME helper.
    if getattr(args, 'ONT_cDNA', False) and input_type in ('fastq', 'fastq.gz'):
        _t0_cdna = _time.perf_counter()
        _cdna_in = input_path
        input_path, input_type = _run_ont_cdna_prepare(
            args, input_path, input_type, work_dir, output_dir, sample_id,
        )
        print(f"[TIMING] ONT-cDNA prep: {_time.perf_counter() - _t0_cdna:.1f}s")
        tracker.record_step(
            'ont_cdna_prepare',
            input_files=[_cdna_in],
            output_files=[input_path],
        )

    # ── Step 0/1: Alignment ───────────────────────────────────────────────────
    if input_type in ('fastq', 'fastq.gz') and not getattr(args, 'skip_alignment', False):
        # sample_id is the canonical one computed above — NOT re-derived from input_path,
        # which Step 0 may have replaced with a derived FASTQ (planning 831).
        if getattr(args, 'short_read', False):
            if getattr(args, 'read2', None):
                _align_mode = 'short-read paired-end (COMPASS panel)'
            elif getattr(args, 'dT_primed_cDNA', False):
                _align_mode = "short-read dT-primed 3'-end (bbmap + bwa)"
            else:
                _align_mode = 'short-read single-end (COMPASS SE panel)'
        else:
            _align_mode = 'long-read consensus'
        print(f"\n[Step 1/3] Aligning ({_align_mode})...")
        print("-" * 50)
        _t0 = _time.perf_counter()
        # Align to bam_dir if specified; otherwise use work_dir (scratch if available)
        align_output_dir = bam_dir if bam_dir else work_dir
        # When --bam-dir forces output to a potentially NFS-backed path, route
        # the consensus sort through scratch to avoid pysam.sort NFS hang.
        _single_ckpt_dir: Optional[str] = str(scratch_dir / 'consensus_ckpt') if scratch_dir else None
        per_aligner_bams, bam_to_correct = _run_alignment(
            input_path=input_path,
            sample_id=sample_id,
            sample_output_dir=align_output_dir,
            genome_path=genome_path,
            annotation_path=annotation_path,
            threads=getattr(args, 'threads', 4),
            parallel_aligners=getattr(args, 'parallel_aligners', False),
            base_aligners=getattr(args, 'base_aligners', None),
            junction_aligners=getattr(args, 'junction_aligners', None),
            chimeric_consensus=getattr(args, 'chimeric_consensus', True),
            ultra_path=getattr(args, 'ultra_path', 'uLTRA'),
            desalt_path=getattr(args, 'desalt_path', 'deSALT'),
            gmap_path=getattr(args, 'gmap_path', 'gmap'),
            gmap_db=getattr(args, 'gmap_db', None),
            mapPacBio_chunks=getattr(args, 'mapPacBio_chunks', 1),
            checkpoint_dir=_single_ckpt_dir,
            short_read=getattr(args, 'short_read', False),
            trust_existing_bams=getattr(args, 'trust_existing_bams', False),
            read2=getattr(args, 'read2', None),
            dt_primed_cdna=getattr(args, 'dT_primed_cDNA', False),
            read_length=getattr(args, 'read_length', 150),
            max_intron=getattr(args, 'max_intron', None),
            resolver_acceptor_classes=getattr(
                args, 'resolver_acceptor_classes', 'canonical'),
        )
        print(f"\nAlignment complete: {bam_to_correct}")
        print(f"[TIMING] Alignment: {_time.perf_counter() - _t0:.1f}s")
        tracker.record_step('align', input_files=[input_path], output_files=[bam_to_correct])
        step_correction = 2
    elif input_type in ('fastq', 'fastq.gz'):
        # canonical sample_id from above (planning 831)
        multialigned_bam = _multialigned_bam_path(sample_id, output_dir)
        _legacy_bam = output_dir / f"{sample_id}.consensus.bam"
        if not multialigned_bam.exists() and _legacy_bam.exists():
            multialigned_bam = _legacy_bam
        if _validate_bam_integrity(multialigned_bam):
            print(f"\n[Step 1/3] Alignment — using existing multialigned.bam")
            if scratch_dir:
                # Stage existing BAM to scratch
                print(f"  Staging to scratch: {scratch_dir}")
                _shutil.copy2(multialigned_bam, scratch_dir / multialigned_bam.name)
                bai = Path(str(multialigned_bam) + '.bai')
                if bai.exists():
                    _shutil.copy2(bai, scratch_dir / bai.name)
                bam_to_correct = scratch_dir / multialigned_bam.name
            else:
                bam_to_correct = multialigned_bam
        else:
            print(
                f"ERROR: --skip-alignment set but no multialigned.bam found: {multialigned_bam}",
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
            from ...consensus.corrected_consensus import merge_corrected_tsvs, identify_cat5_candidates, write_corrected_consensus_bam, _stage_raw_bams
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
                # Use HP-edit-distance winner selection. If per-aligner
                # corrected BAMs were not materialized, score lazily from raw
                # BAMs + corrected TSVs.
                from rectify.utils.genome import load_genome as _load_genome_for_merge
                _merge_genome = _load_genome_for_merge(str(genome_path))
                _raw_bams_s2 = {
                    a: str(p) for a, p in per_aligner_bams.items()
                } if per_aligner_bams else {}
                from rectify.data import resolve_min_junction_anchor_bp as _resolve_anchor
                with _stage_raw_bams(_raw_bams_s2) as _staged_s2:
                    corrected_tsv = merge_corrected_tsvs(
                        per_aligner_tsvs=per_aligner_tsvs,
                        output_tsv=work_dir / 'corrected_reads.tsv',
                        summary_tsv=_summary_tsv,
                        per_aligner_corrected_bams={
                            a: str(p) for a, p in per_aligner_corrected_bams.items()
                        } if per_aligner_corrected_bams else None,
                        per_aligner_raw_bams=_staged_s2 if not per_aligner_corrected_bams else None,
                        genome=_merge_genome,
                        overhang_table=_overhang_table,
                        lazy_scoring_workers=max(1, int(getattr(args, 'threads', 1) or 1)),
                        min_junction_anchor_bp=_resolve_anchor(args),
                    )
                    # Identify Cat5 candidates (reads where aligners contribute unique introns)
                    _cat5_tsv = _per_aligner_dir / 'cat5_candidates.tsv'
                    identify_cat5_candidates(per_aligner_tsvs, output_tsv=_cat5_tsv)
                    # Emit the final corrected (rectified) BAM — symmetric with the
                    # chunked path.  Named <sample>.rectified.bam; a legacy
                    # corrected_consensus.bam symlink is dropped for back-compat.
                    # Non-fatal: merged TSV is the primary output; missing BAM is
                    # recoverable via post-hoc write_corrected_consensus_bam.
                    _consensus_bam_out = _final_rectified_bam_path(sample_id, work_dir)
                    try:
                        _consensus_stats = write_corrected_consensus_bam(
                            per_aligner_raw_bams=_staged_s2,
                            per_aligner_tsvs=per_aligner_tsvs,
                            merged_tsv=corrected_tsv,
                            output_bam=_consensus_bam_out,
                            genome=_merge_genome,
                            threads=max(1, int(getattr(args, 'threads', 1) or 1)),
                            strict=True,
                        )
                        _emit_legacy_consensus_symlink(_consensus_bam_out)
                        print(
                            f"    {_consensus_bam_out.name} written"
                            f" ({_consensus_stats.get('written', '?')} reads)",
                            flush=True,
                        )
                    except Exception as _cbam_exc:
                        print(
                            f"    {_consensus_bam_out.name} write failed"
                            f" (non-fatal): {_cbam_exc}",
                            file=sys.stderr,
                        )
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
        _cbam_scratch = _final_rectified_bam_path(sample_id, corrected_tsv.parent)
        if _cbam_scratch.exists():
            _cbam_dst = _final_rectified_bam_path(sample_id, output_dir)
            _shutil.copy2(_cbam_scratch, _cbam_dst)
            _cbam_bai = Path(str(_cbam_scratch) + '.bai')
            if _cbam_bai.exists():
                _shutil.copy2(_cbam_bai, Path(str(_cbam_dst) + '.bai'))
            _emit_legacy_consensus_symlink(_cbam_dst)
        corrected_tsv = _oak_corrected_tsv
        work_dir = output_dir  # analysis outputs go directly to Oak

    # Add sample column (needed by analyze even for single sample)
    #
    # 🔴 `corrected_tsv` may be the per-region MANIFEST, not a reads table — since Commit B made
    # `corrected_reads.manifest.tsv` the default artifact, `_run_correction` returns the manifest
    # whenever it exists. Writing a `sample` column into the MANIFEST appends an 8th column to its
    # 7-column header, which breaks every `_is_manifest()` signature check downstream
    # (`analyze/loaders.py`, `bam/bam_writer.py`, `commands/qc_command.py`). Those then fall back to
    # reading the manifest as if it were a reads TSV — one bogus row — so **analyze silently
    # produced no clusters at all and run-all still exited 0**. Observed end-to-end 2026-08-08:
    # `cpa_clusters.tsv` absent, empty `tables/`, "Pipeline Complete!".
    #
    # The sample column belongs in the REGION TSVs, which is what
    # `_load_manifest_as_dataframe` actually reads. The manifest is an index and must stay a
    # 7-column index.
    # canonical sample_id from above (planning 831) — was re-derived from a Step-0-rewritten
    # input_path here too, so the TSV's `sample` column disagreed with the file names.
    try:
        import pandas as pd
        from ...bam.tsv_partition import load_manifest
        from ...analyze.loaders import _is_manifest as _loader_is_manifest

        targets = []
        if _loader_is_manifest(str(corrected_tsv)):
            targets = [e['tsv_path'] for e in load_manifest(Path(corrected_tsv))]
        else:
            targets = [Path(corrected_tsv)]

        for _t in targets:
            df = pd.read_csv(_t, sep='\t')
            df['sample'] = sample_id
            df.to_csv(_t, sep='\t', index=False)
    except Exception as _e:
        # Still not fatal — analyze can work without the sample column — but it must not be
        # SILENT. The bug above survived because this handler swallowed everything.
        print(f"  WARNING: could not add the 'sample' column ({type(_e).__name__}: {_e}); "
              f"analyze will fall back to a single unnamed sample", file=sys.stderr)

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

    # ── Browser pack (LAST step; fail-soft) ──────────────────────────────────
    # Per-library QC + analysis.json, folded from the tables analyze just wrote.
    # Placed here rather than at the very end of the function so a later
    # optional stage (junction aggregation, DRS poly(A) restore) cannot cost the
    # browser bundle; nothing downstream feeds the packer. work_dir == output_dir
    # by this point on the scratch path too, so analysis.json lands durably.
    # ``sample_id`` is now bound UNCONDITIONALLY at the top of this function (from the original
    # input, before Step 0 — planning 831), so this no longer needs a local re-derivation. It used
    # to, because the name was only bound on the FASTQ branches; re-deriving it here also meant the
    # browser bundle was named ``stage1_consensus`` on the ONT-cDNA path.
    _bp_sample_id = sample_id
    _run_browser_pack(
        analysis_dir=work_dir,
        samples=[{
            'sample_id': _bp_sample_id,
            'bam': bam_to_correct,
            'tsv': corrected_tsv,
        }],
        args=args,
        annotation_path=annotation_path,
        modality=('DRS' if drs_mode else
                  ('cDNA' if getattr(args, 'ONT_cDNA', False) else '')),
    )

    # ── Re-aligner Stations B + C (fail-soft; B opt-in via --triage,
    #    C = the pool-gate junction admission report, on by default) ─────────
    _run_station_bc(
        work_dir=work_dir,
        sample_id=_bp_sample_id,
        fallback_bam=bam_to_correct,
        genome_path=genome_path,
        annotation_path=annotation_path,
        args=args,
    )

    # ── Junction aggregation ─────────────────────────────────────────────────
    if genome_path and annotation_path:
        try:
            _run_junction_aggregation(
                bam_path=bam_to_correct,
                genome_path=genome_path,
                annotation_path=annotation_path,
                output_dir=work_dir,
                config=vars(args),
            )
        except Exception as e:
            print(f"WARNING: Junction aggregation failed (non-fatal): {e}", file=sys.stderr)

    print(f"[TIMING] Analysis:   {_time.perf_counter() - _t0:.1f}s")

    # ── Step 4 (DRS only, opt-in): Restore poly(A) from winning aligner's raw BAM ─
    # Pulls each read from the winning aligner's raw (pre-correction) BAM and
    # re-attaches the original Dorado poly(A) tail from parquet as a 3' soft clip.
    # The result (corrected_polya.bam) shows the full read as Dorado called it;
    # compare against corrected.bam in IGV to see exactly what Rectify changed.
    if (drs_mode and getattr(args, 'write_polya_bam', False)
            and drs_metadata_path and drs_metadata_path.exists()):
        from ..restore_polya_command import restore_polya_softclips
        _orig_stem = Path(args.input).stem
        # Write Step-4 BAM to scratch (when available) and pysam-sort there;
        # the final sync_to_oak below ships it to NFS. Sorting on NFS hangs
        # under concurrent load, which scratch staging exists to avoid.
        _restore_dir = scratch_dir if scratch_dir else output_dir
        _restored_bam = _restore_dir / f"{_orig_stem}.corrected_polya.bam"
        _aligner_bams_for_restore = {k: str(v) for k, v in per_aligner_bams.items()
                                     if Path(str(v)).exists()}
        if corrected_tsv and Path(corrected_tsv).exists() and _aligner_bams_for_restore:
            print(f"\n[Step 4] DRS poly(A) restoration from winning aligner raw BAMs...")
            print("-" * 50)
            _t0_sc = _time.perf_counter()
            try:
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
                Path(_restored_tmp).replace(_restored_bam)
                _pysam_sc.index(str(_restored_bam))
                print(f"  Restored: {_sc_stats.get('restored', 0):,} reads")
                print(f"[TIMING] DRS restore: {_time.perf_counter() - _t0_sc:.1f}s")
                # Record provenance against the canonical Oak path so the
                # manifest survives scratch teardown.
                _restored_bam_oak = output_dir / _restored_bam.name
                tracker.record_step(
                    'drs_restore',
                    input_files=list(_aligner_bams_for_restore.values()) + [str(drs_metadata_path)],
                    output_files=[_restored_bam_oak],
                )
            except Exception as e:
                print(f"WARNING: DRS poly(A) restore failed (non-fatal): {e}", file=sys.stderr)
        else:
            print(
                f"\n[Step 4] DRS poly(A) restoration skipped "
                f"— corrected TSV or aligner BAMs not found",
                file=sys.stderr,
            )

    # ── Copy from scratch to Oak ─────────────────────────────────────────────
    # rmtree runs in finally so a sync_to_oak failure (which now raises
    # RuntimeError on rsync non-zero) does not orphan the scratch directory.
    if scratch_dir:
        try:
            print(f"\nCopying outputs to Oak: {output_dir}")
            _t0_sync = _time.perf_counter()
            # When bam_dir is set, BAMs went directly there (not scratch), so
            # exclude_aligner_bams only matters when BAMs are on scratch. When
            # keep_aligner_bams is False (default), skip per-aligner BAMs in
            # the scratch→Oak sync to save disk space.
            _exclude_aligner = (not keep_aligner_bams) and (bam_dir is None)
            sync_to_oak(scratch_dir, output_dir, exclude_aligner_bams=_exclude_aligner)
            print(f"[TIMING] Sync to Oak: {_time.perf_counter() - _t0_sync:.1f}s")
        finally:
            _shutil.rmtree(scratch_dir, ignore_errors=True)

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
