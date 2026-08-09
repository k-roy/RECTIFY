"""
Per-pipeline-stage runners for ``rectify run-all``.

Each function is a thin wrapper around one stage of the pipeline:

- ``_run_alignment``               — multi-aligner consensus (or re-use of an
                                     existing rectified BAM).
- ``_run_correction``              — single ``rectify correct`` invocation.
- ``_run_correction_per_aligner``  — per-aligner correction (one call per BAM)
                                     into ``output_dir/per_aligner_corrected/``.
- ``_combine_corrected_tsvs``      — concat per-sample corrected TSVs.
- ``_run_analysis``                — ``rectify analyze`` on a (combined) TSV.
- ``_run_junction_aggregation``    — junction aggregation with partial rescue
                                     (no-op when annotation has no introns).

The resume-on-restart logic (skip-if-output-exists) is sample-level and stays
in ``single_sample`` / ``multi_sample``; these runners just execute.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .helpers import (
    _bam_has_md_tags,
    _collect_per_aligner_bams,
    _multialigned_bam_path,
    _validate_bam_integrity,
)


def _run_alignment(
    input_path: Path,
    sample_id: str,
    sample_output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    threads: int,
    parallel_aligners: bool = False,
    base_aligners: Optional[List[str]] = None,
    junction_aligners: Optional[List[str]] = None,
    chimeric_consensus: bool = True,
    ultra_path: str = 'uLTRA',
    desalt_path: str = 'deSALT',
    gmap_path: str = 'gmap',
    gmap_db: Optional[str] = None,
    mapPacBio_chunks: int = 1,
    checkpoint_dir: Optional[str] = None,
    short_read: bool = False,
    trust_existing_bams: bool = False,
    read2: Optional[Path] = None,
    read_length: int = 150,
) -> Tuple[Dict[str, Path], Path]:
    """
    Run multi-aligner alignment and selection, or return existing multialigned.bam.

    Default aligners: minimap2 + mapPacBio + gapmm2 (long-read Tier 1).
    Pass base_aligners to restrict or change the set (e.g. ['mapPacBio']).
    Pass short_read=True to use bbmap + bwa instead of the long-read panel.
    Pass short_read=True *and* read2=<R2 FASTQ> to use the paired-end COMPASS
    panel (bbmap + STAR×2 + HISAT2×2 + magicblast + gsnap) — the same set
    `rectify align --short-read --read2 --aligners all` selects.
    Pass junction_aligners=[] to disable uLTRA + deSALT.

    Skips automatically if the multialigned.bam already exists — safe to re-run.

    Returns
    -------
    Tuple of (per_aligner_bams, multialigned_bam) where per_aligner_bams maps
    aligner name → BAM path for each per-aligner BAM found on disk.
    """
    multialigned_bam = _multialigned_bam_path(sample_id, sample_output_dir)

    # Backward compatibility: accept the older ``.consensus.bam`` name from prior
    # runs.  We deliberately do NOT accept the pre-2026-06-25 ``.rectified.bam``
    # name as the alignment artifact here: after the rename ``.rectified.bam`` is
    # the FINAL corrected output, so reusing it as the pre-correction artifact
    # would feed the corrected file back into correction.  A pre-rename run simply
    # re-aligns (correct, just slower) — the swap hazard is avoided.
    _legacy_bam = sample_output_dir / f"{sample_id}.consensus.bam"
    if not multialigned_bam.exists() and _legacy_bam.exists():
        multialigned_bam = _legacy_bam

    # Compute run provenance once — used by both the multialigned.bam reuse gate
    # and the per-aligner BAM provenance filter.
    import sys as _sys
    from rectify.utils.bam_provenance import (
        compute_run_provenance,
        expected_provenance_for_aligner,
        matches_strict,
        read_sidecar,
    )
    _run_provenance = compute_run_provenance(command=_sys.argv)

    if _validate_bam_integrity(multialigned_bam):
        if trust_existing_bams:
            print(f"    Skipping alignment — multialigned.bam exists (--trust-existing-bams): {multialigned_bam}")
            per_aligner_bams = _collect_per_aligner_bams(
                sample_id, sample_output_dir,
                run_provenance=None, trust_existing_bams=True,
            )
            return per_aligner_bams, multialigned_bam
        _expected = expected_provenance_for_aligner(_run_provenance, "consensus")
        _stored = read_sidecar(multialigned_bam)
        _ok, _reason = matches_strict(_stored, _expected)
        if _ok:
            print(f"    Skipping alignment — multialigned.bam exists (provenance match): {multialigned_bam}")
            per_aligner_bams = _collect_per_aligner_bams(
                sample_id, sample_output_dir,
                run_provenance=_run_provenance, trust_existing_bams=False,
            )
            return per_aligner_bams, multialigned_bam
        print(
            f"    multialigned.bam provenance mismatch ({_reason}); re-running alignment. "
            "Pass --trust-existing-bams to override."
        )

    # Panel selection. NOTE: run-all passes an EXPLICIT aligner list to run_align
    # (see align_args below), so align_command's own "all"-expansion never fires
    # from this path — the paired-COMPASS choice has to be made here too, or
    # --read2 would silently run bbmap+bwa paired instead of the COMPASS panel.
    from ..align_command import COMPASS_PE_ALIGNERS
    if base_aligners is not None:
        _base_aligners = base_aligners
    elif short_read and read2 is not None:
        _base_aligners = list(COMPASS_PE_ALIGNERS)
    elif short_read:
        _base_aligners = ['bbmap', 'bwa']
    else:
        _base_aligners = ['minimap2', 'mapPacBio', 'gapmm2']
    # Junction-aligner default depends on --short-read: BBMap's intronlen=20
    # already covers splicing for short reads, so uLTRA/deSALT aren't useful
    # and shouldn't be the default. Long-read protocols still default to
    # [uLTRA, deSALT]. Explicit user list (including []) overrides.
    if junction_aligners is None:
        _junction_aligners = [] if short_read else ['uLTRA', 'deSALT']
    else:
        _junction_aligners = junction_aligners
    all_aligners = _base_aligners + _junction_aligners
    aligner_desc = ' + '.join(all_aligners)
    print(f"    Running {len(all_aligners)}-aligner consensus ({aligner_desc})...")
    from ..align_command import run_align

    align_args = argparse.Namespace(
        reads=input_path,
        genome=genome_path,
        output_dir=sample_output_dir,
        annotation=annotation_path,
        threads=threads,
        aligners=_base_aligners,
        short_read=short_read,
        # Paired-end short-read (COMPASS panel): mate-2 FASTQ + read length for
        # STAR's --sjdbOverhang. Both are None/default-safe for every other mode.
        read2=read2,
        read_length=read_length,
        junction_aligners=_junction_aligners,
        no_consensus=False,
        chimeric_consensus=chimeric_consensus,
        junc_bonus=9,
        junc_bed=None,
        parallel_aligners=parallel_aligners,
        minimap2_path='minimap2',
        mapPacBio_path='mapPacBio.sh',
        gapmm2_path='gapmm2',
        ultra_path=ultra_path,
        desalt_path=desalt_path,
        gmap_path=gmap_path,
        gmap_db=gmap_db,
        mapPacBio_chunks=mapPacBio_chunks,
        mapPacBio_chunk_idx=None,  # merge mode: look for existing chunk BAMs
        # Use the canonical sample_id as prefix so the merged BAM lands at
        # <sample_id>.multialigned.bam (matching _multialigned_bam_path). With an
        # empty prefix, align_command falls back to args.reads.stem — which
        # for DRS inputs is "<sample>_trimmed" after Step 0, breaking the
        # post-align multialigned.bam lookup below.
        prefix=sample_id,
        keep_sam=False,
        sort=True,
        index=True,
        verbose=False,
        checkpoint_dir=checkpoint_dir,
        trust_existing_bams=trust_existing_bams,
        _run_provenance=_run_provenance,
    )

    # Temporarily hide scheduler array env vars so consensus.py runs in
    # single-sample mode. In run-all, array indices are for sample-level
    # parallelism; we never want within-sample read-level partitioning here.
    import os as _os
    import time as _time
    from datetime import datetime as _dt, timezone as _tz
    _align_started_at = _dt.now(_tz.utc).isoformat()
    _align_t0 = _time.perf_counter()
    _array_vars = {k: _os.environ.pop(k) for k in (
        # SLURM
        'SLURM_ARRAY_TASK_ID', 'SLURM_ARRAY_TASK_COUNT',
        'SLURM_ARRAY_TASK_MAX', 'SLURM_ARRAY_TASK_MIN',
        'SLURM_ARRAY_TASK_STEP',
        # UGE/SGE
        'SGE_TASK_ID', 'SGE_TASK_FIRST', 'SGE_TASK_LAST', 'SGE_TASK_STEPSIZE',
        # PBS/Torque
        'PBS_ARRAY_INDEX', 'PBS_ARRAYID',
    ) if k in _os.environ}
    try:
        rc = run_align(align_args)
    finally:
        _os.environ.update(_array_vars)
    if rc != 0:
        raise RuntimeError(f"Triple-aligner failed for {input_path}")

    if not multialigned_bam.exists():
        raise RuntimeError(
            f"Alignment completed but multialigned.bam not found: {multialigned_bam}"
        )

    per_aligner_bams = _collect_per_aligner_bams(
        sample_id, sample_output_dir,
        run_provenance=_run_provenance, trust_existing_bams=trust_existing_bams,
    )

    try:
        from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar
        _align_inputs: dict = {'reads': input_path, 'genome': genome_path}
        if annotation_path is not None:
            _align_inputs['annotation'] = annotation_path
        _align_outputs: dict = {'multialigned_bam': multialigned_bam}
        _align_outputs.update({
            f'per_aligner_bam_{name}': bam
            for name, bam in per_aligner_bams.items()
        })
        _align_record = ProvenanceRecord.from_components(
            stage='align',
            sample_id=sample_id,
            sample_output_dir=sample_output_dir,
            started_at=_align_started_at,
            completed_at=_dt.now(_tz.utc).isoformat(),
            exit_status=0,
            inputs=_align_inputs,
            outputs=_align_outputs,
            stats={
                'wall_seconds': _time.perf_counter() - _align_t0,
                'n_aligners': len(all_aligners),
                'aligners': all_aligners,
                'parallel_aligners': parallel_aligners,
                'threads': threads,
            },
        )
        write_stage_sidecar(_align_record, sample_output=sample_output_dir)
    except Exception as _sc_exc:
        import logging as _logging
        _logging.getLogger(__name__).warning("Failed to write align sidecar: %s", _sc_exc)

    return per_aligner_bams, multialigned_bam


def _run_correction(
    bam_path: Path,
    output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    args,
    aligner_bams: Optional[List[Path]] = None,
    reuse_pool_container: Optional[list] = None,
) -> Path:
    """
    Run rectify correct on a BAM, writing corrected_reads.tsv into output_dir.
    Returns path to the corrected TSV.

    aligner_bams: Optional list of per-aligner BAM paths to feed Module 2H
        (junction refinement). Required to match the CLAUDE.md "validated"
        correct-first pipeline order. When provided, every per-aligner
        correction gets the full cross-aligner junction candidate pool.
    """
    from .. import correct_command

    corrected_tsv = output_dir / 'corrected_reads.tsv'

    # Indel correction and variant-aware rescue require MD tags in the BAM.
    # Disable them gracefully when MD tags are absent (e.g. rectified BAMs
    # generated by older runs without --MD).
    has_md = _bam_has_md_tags(bam_path)
    # --dT-primed-cDNA enables poly-A trimming / indel-correction modules.
    # Old --no-polya-sequenced is a deprecated alias for NOT passing --dT-primed-cDNA.
    # Default (no flag): direct RNA mode — poly-A IS in the read.
    polya_seq = getattr(args, 'dT_primed_cDNA', False)
    # --ONT-cDNA selects the ONT PCR-cDNA protocol wrapper, which resolves the RNA
    # strand PER READ (ro/XO tag -> gene overlap) instead of reading it off the BAM
    # strand. Required for SQK-PCB114-style libraries, where reads arrive in both
    # orientations; see rectify.core.correct.protocols.ont_cdna.
    ont_cdna = getattr(args, 'ONT_cDNA', False)
    skip_indel = not has_md  # Can't correct indels without MD tags
    skip_variant = not has_md

    if not has_md and polya_seq:
        print(
            "    Note: BAM lacks MD tags — indel correction disabled. "
            "Re-run alignment with current RECTIFY to enable it."
        )

    # Auto-generate report path alongside the corrected TSV
    report_path = corrected_tsv.parent / (corrected_tsv.stem + '_report.html')

    # Derive sample prefix from input BAM stem (strip .rectified / .consensus suffixes)
    _stem = bam_path.stem
    for _sfx in ('.rectified', '.consensus'):
        if _stem.endswith(_sfx):
            _stem = _stem[:-len(_sfx)]
            break
    corrected_bam_path = (
        output_dir / f"{_stem}.rectified_corrected_3end.bam"
        if getattr(args, 'write_corrected_bam', False) else None
    )
    softclipped_bam_path = (
        output_dir / f"{_stem}.rectified_pA_tail_trimmed.bam"
        if getattr(args, 'write_softclip_bam', False) else None
    )

    # Module 2H (junction refinement) requires the full per-aligner BAM pool to
    # build the candidate junction set. Falling back to args.aligner_bams keeps
    # backwards-compat for callers that set the attribute directly.
    _aligner_bams = aligner_bams or getattr(args, 'aligner_bams', None) or []

    correct_args = argparse.Namespace(
        input=bam_path,
        genome=genome_path,
        annotation=annotation_path,
        output=corrected_tsv,
        corrected_bam=corrected_bam_path,
        softclipped_bam=softclipped_bam_path,
        netseq_dir=getattr(args, 'netseq_dir', None),
        organism=getattr(args, 'organism', None),
        aligner=getattr(args, 'aligner', 'minimap2'),
        aligner_bams=[str(p) for p in _aligner_bams],
        dT_primed_cDNA=polya_seq,
        polya_sequenced=polya_seq,  # deprecated attr kept for compat
        ONT_cDNA=ont_cdna,
        threads=getattr(args, 'threads', 4),
        filter_spikein=getattr(args, 'filter_spikein', None),
        streaming=getattr(args, 'streaming', False),
        # Defaults
        min_mapq=10,
        skip_secondary=True,
        skip_supplementary=True,
        skip_ag_check=False,
        skip_atract_check=False,
        skip_polya_trim=False,
        skip_indel_correction=skip_indel,
        skip_variant_aware=skip_variant,
        polya_model=None,
        report=report_path,
        max_downstream_a=20,
        chunk_size=getattr(args, 'chunk_size', 10000),
        junction_penalty_table=getattr(args, 'junction_penalty_table', None),
        str_penalty_table=getattr(args, 'str_penalty_table', None),
        junction_overhang_table=getattr(args, 'junction_overhang_table', None),
        junction_pool_cache=getattr(args, 'junction_pool_cache', None),
        debug=False,
        verbose=False,
        # Pool reuse: injected by _run_correction_per_aligner for ac>1; None otherwise.
        reuse_pool_container=reuse_pool_container,
    )

    correct_command.run(correct_args)

    # Commit B made the per-region manifest the default artifact:
    # corrected_reads.tsv is renamed to corrected_reads.region_000.tsv and
    # corrected_reads.manifest.tsv becomes the canonical entry point. Treat
    # either as valid output; downstream loaders (_read_corrected_tsv_or_manifest)
    # accept both. Pass --emit-merged-tsv to force the legacy concatenated form.
    manifest_path = corrected_tsv.parent / (corrected_tsv.stem + '.manifest.tsv')
    if manifest_path.exists():
        return manifest_path
    if corrected_tsv.exists():
        return corrected_tsv
    raise RuntimeError(
        f"Correction did not produce output: neither {corrected_tsv} "
        f"nor {manifest_path} exists"
    )


def _per_aligner_canonical_output(aligner_output_dir: Path) -> Optional[Path]:
    """Return the canonical corrected-output path for one aligner, or None.

    Prefers ``corrected_reads.manifest.tsv`` (Commit B default); falls back to
    ``corrected_reads.tsv`` (legacy ``--emit-merged-tsv`` mode). Used by both
    the resume skip-check and the post-call success validation in
    ``_run_correction_per_aligner``.
    """
    manifest_path = aligner_output_dir / 'corrected_reads.manifest.tsv'
    if manifest_path.exists():
        return manifest_path
    tsv_path = aligner_output_dir / 'corrected_reads.tsv'
    if tsv_path.exists():
        return tsv_path
    return None


def _readable_corrected_tsv(sample_dir: Path) -> Optional[Path]:
    """Return a path that can be ``pd.read_csv``'d as a READS table, or None.

    🔴 Deliberately different from :func:`_per_aligner_canonical_output`, which prefers the
    per-region *manifest*. The multi-sample analyze manifest carries one ``path`` per sample and
    every manifest-mode pass in ``analyze/manifest.py`` does ``pd.read_csv(path)`` on it directly
    (lines 244, 437, 583, 808). Handing those a per-region manifest is wrong twice over: it is not
    a reads table, and the header check silently ``continue``s, so the sample vanishes from the
    analysis with no error.

    Before this existed, ``multi_sample.py`` hardcoded ``corrected_reads.tsv`` — which has not
    existed since Commit B made manifest-only the default (the merged TSV is renamed to
    ``corrected_reads.region_000.tsv``). Pass 1 survived on the position index, then Pass 1b died
    with ``FileNotFoundError``. Reproduced end-to-end 2026-08-08.
    """
    flat = sample_dir / 'corrected_reads.tsv'
    if flat.exists():
        return flat

    manifest_path = sample_dir / 'corrected_reads.manifest.tsv'
    if not manifest_path.exists():
        return None

    try:
        from ...bam.tsv_partition import load_manifest
        regions = [Path(str(e['tsv_path'])) for e in load_manifest(manifest_path)]
    except Exception as exc:
        print(f"  WARNING: could not read {manifest_path}: {exc}", file=sys.stderr)
        return None

    regions = [r for r in regions if r.exists()]
    if len(regions) == 1:
        return regions[0]
    if not regions:
        return None
    # A single-path sample manifest cannot express several shards. `correct` currently emits a
    # one-entry manifest, so this is unreachable today — but it must fail LOUDLY if that changes,
    # not quietly analyse one shard and drop the rest.
    print(f"  WARNING: {manifest_path} names {len(regions)} region TSVs; the analyze sample "
          f"manifest holds one path per sample, so this sample would be analysed from only part "
          f"of its reads. Re-run correct with --emit-merged-tsv, or use "
          f"'rectify export-merged-tsv' to produce a single TSV.", file=sys.stderr)
    return None


def _run_correction_per_aligner(
    per_aligner_bams: Dict[str, Path],
    output_dir: Path,
    genome_path: Path,
    annotation_path: Optional[Path],
    args,
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    """
    Run ``rectify correct`` on each per-aligner BAM independently.

    Writes per-aligner outputs to ``output_dir/per_aligner_corrected/{aligner}/``.
    Skips any aligner whose canonical output (``corrected_reads.manifest.tsv``
    or ``corrected_reads.tsv``) already exists and parses (safe to resume).

    Reads ``args.aligner_concurrency`` (default ``'auto'``) and resolves it via
    :func:`rectify.core.utils.resources.resolve_aligner_concurrency`. The
    resolved value is currently informational only — the loop runs aligners
    sequentially. The shared cross-aligner region queue (Axis B performance)
    is deferred until ``correct_command.run`` exposes a region-task API.
    Use ``--aligner-concurrency 1`` to assert serial behavior explicitly.

    Returns
    -------
    Tuple ``(per_aligner_tsvs, per_aligner_corrected_bams)``:
        - ``per_aligner_tsvs``: Dict mapping aligner name → canonical
          corrected-output path (manifest TSV by default, raw TSV under
          ``--emit-merged-tsv``) for each aligner whose correction succeeded.
        - ``per_aligner_corrected_bams``: Dict mapping aligner name →
          rectified_corrected_3end.bam path for each aligner whose correction
          succeeded AND produced a corrected BAM. Empty when no corrected BAMs
          were materialized (e.g. ``--write-corrected-bam`` omitted upstream).
          When non-empty, this dict is the input ``merge_corrected_tsvs`` needs
          to activate HP-edit-distance mode for winner selection.
    """
    per_aligner_dir = output_dir / 'per_aligner_corrected'
    per_aligner_dir.mkdir(parents=True, exist_ok=True)

    # Commit C preparatory plumb: resolve the shared per-aligner concurrency
    # budget. Default policy keeps the loop sequential (auto -> 1 on laptops;
    # cluster nodes also stay serial until correct_command.run exposes a
    # region-task API for cross-aligner sharing). The value is logged so
    # downstream callers and tests can observe the resolved budget.
    from ...utils.resources import (
        detect_machine_class as _detect_machine_class,
        resolve_aligner_concurrency as _resolve_aligner_concurrency,
    )
    _ac_value = getattr(args, 'aligner_concurrency', 'auto') or 'auto'
    _total_threads = int(getattr(args, 'threads', 1) or 1)
    try:
        _resolved_ac = _resolve_aligner_concurrency(
            str(_ac_value), _total_threads, _detect_machine_class()
        )
    except ValueError as _ac_exc:
        print(
            f"    WARNING: --aligner-concurrency={_ac_value!r} invalid "
            f"({_ac_exc}); falling back to 1",
            file=sys.stderr,
        )
        _resolved_ac = 1
    # Shared pool path (ac > 1): one multiprocessing.Pool across all aligners.
    # MD tags must be uniform — apply_indel_correction is baked into the pool
    # worker state once at creation time and can't vary per aligner.
    _use_shared_pool = False
    if _resolved_ac > 1:
        _has_md_per_aligner = {
            name: _bam_has_md_tags(bam_path)
            for name, bam_path in per_aligner_bams.items()
        }
        if len(set(_has_md_per_aligner.values())) > 1:
            print(
                f"    WARNING: per-aligner BAMs have inconsistent MD tags "
                f"({_has_md_per_aligner}); falling back to sequential correction.",
                file=sys.stderr,
            )
        else:
            _use_shared_pool = True

    _mode_label = "shared pool" if _use_shared_pool else "sequential"
    print(
        f"    Aligner concurrency: --aligner-concurrency={_ac_value} "
        f"-> {_resolved_ac} ({_mode_label})"
    )

    per_aligner_tsvs: Dict[str, Path] = {}
    per_aligner_corrected_bams: Dict[str, Path] = {}

    if not _use_shared_pool:
        # Sequential path: one pool per aligner, or no pool (n_threads=1).
        # This branch is byte-identical to the pre-C-II.1 behavior.
        for aligner_name, bam_path in per_aligner_bams.items():
            aligner_output_dir = per_aligner_dir / aligner_name
            aligner_output_dir.mkdir(exist_ok=True)

            # Resume skip-check: accept manifest (Commit B default) OR raw TSV
            # (--emit-merged-tsv legacy). The previous tsv-only check meant
            # manifest-only runs never skipped on resume and silently re-ran the
            # whole correction stage on every attempt.
            _canonical = _per_aligner_canonical_output(aligner_output_dir)
            _output_valid = False
            if _canonical is not None:
                try:
                    import pandas as pd
                    _df = pd.read_csv(_canonical, sep='\t', nrows=1)
                    _output_valid = len(_df.columns) > 0
                except Exception:
                    _output_valid = False

            if _output_valid:
                print(f"    [{aligner_name}] Skipping — output exists: {_canonical}")
                per_aligner_tsvs[aligner_name] = _canonical
                continue

            print(f"    [{aligner_name}] Correcting {bam_path.name}...")
            # rectify correct uses pysam.fetch() which requires a coordinate-sorted,
            # indexed BAM.  Aligners that sort internally (uLTRA, deSALT) have a .bai
            # placed by align_command.py; aligners that emit unsorted output
            # (mapPacBio, minimap2, gapmm2) do not.  Sort + index to a temp file when
            # the .bai is absent, then clean up after correction.
            sorted_bam: Optional[Path] = None
            correction_bam = bam_path
            bai_path = Path(str(bam_path) + '.bai')
            if not bai_path.exists():
                try:
                    import pysam as _pysam
                    sorted_bam = bam_path.with_suffix('.coord_sorted.bam')
                    print(f"    [{aligner_name}] Coordinate-sorting for correction...")
                    _pysam.sort('-m', '1G', '-o', str(sorted_bam), str(bam_path))
                    _pysam.index(str(sorted_bam))
                    correction_bam = sorted_bam
                except Exception as _sort_exc:
                    print(
                        f"    WARNING: [{aligner_name}] sort failed ({_sort_exc}); "
                        "skipping correction",
                        file=sys.stderr,
                    )
                    if sorted_bam is not None and sorted_bam.exists():
                        sorted_bam.unlink(missing_ok=True)
                    continue
            try:
                # Pass ALL per-aligner BAMs so Module 2H (junction refinement) can
                # build the cross-aligner candidate junction pool — the documented
                # "correct-first" pipeline order in CLAUDE.md. Without this every
                # per-aligner correction would run with Module 2H silently disabled.
                _all_aligner_bams = list(per_aligner_bams.values())
                _run_correction(
                    bam_path=correction_bam,
                    output_dir=aligner_output_dir,
                    genome_path=genome_path,
                    annotation_path=annotation_path,
                    args=args,
                    aligner_bams=_all_aligner_bams,
                )
            except Exception as exc:
                print(
                    f"    WARNING: [{aligner_name}] correction failed: {exc}",
                    file=sys.stderr,
                )
            finally:
                if sorted_bam is not None:
                    sorted_bam.unlink(missing_ok=True)
                    sorted_bai = Path(str(sorted_bam) + '.bai')
                    sorted_bai.unlink(missing_ok=True)

            _canonical_after = _per_aligner_canonical_output(aligner_output_dir)
            if _canonical_after is not None:
                per_aligner_tsvs[aligner_name] = _canonical_after
                # Find the corrected BAM emitted by _run_correction.  Naming
                # convention: ``{input_bam.stem stripped of .rectified/.consensus}
                # .rectified_corrected_3end.bam``.  When ``--write-corrected-bam``
                # was disabled upstream this file won't exist; that's fine —
                # absence triggers legacy 5-key sort in merge_corrected_tsvs.
                _stem = (sorted_bam.stem if sorted_bam is not None else bam_path.stem)
                for _sfx in ('.rectified', '.consensus', '.coord_sorted'):
                    if _stem.endswith(_sfx):
                        _stem = _stem[:-len(_sfx)]
                        break
                _corr_bam = aligner_output_dir / f"{_stem}.rectified_corrected_3end.bam"
                if _corr_bam.exists():
                    per_aligner_corrected_bams[aligner_name] = _corr_bam
            else:
                print(
                    f"    WARNING: [{aligner_name}] correction produced no output",
                    file=sys.stderr,
                )

        return per_aligner_tsvs, per_aligner_corrected_bams

    # Shared pool path: one multiprocessing.Pool is created by the first aligner's
    # process_bam_file_parallel call (via correct_command.run) and reused for all
    # subsequent aligners.  Pool lifecycle is managed here via try/finally.
    _pool_container: list = []
    try:
        for aligner_name, bam_path in per_aligner_bams.items():
            aligner_output_dir = per_aligner_dir / aligner_name
            aligner_output_dir.mkdir(exist_ok=True)

            _canonical = _per_aligner_canonical_output(aligner_output_dir)
            _output_valid = False
            if _canonical is not None:
                try:
                    import pandas as pd
                    _df = pd.read_csv(_canonical, sep='\t', nrows=1)
                    _output_valid = len(_df.columns) > 0
                except Exception:
                    _output_valid = False

            if _output_valid:
                print(f"    [{aligner_name}] Skipping — output exists: {_canonical}")
                per_aligner_tsvs[aligner_name] = _canonical
                continue

            print(f"    [{aligner_name}] Correcting {bam_path.name} (shared pool)...")
            sorted_bam = None
            correction_bam = bam_path
            bai_path = Path(str(bam_path) + '.bai')
            if not bai_path.exists():
                try:
                    import pysam as _pysam
                    sorted_bam = bam_path.with_suffix('.coord_sorted.bam')
                    print(f"    [{aligner_name}] Coordinate-sorting for correction...")
                    _pysam.sort('-m', '1G', '-o', str(sorted_bam), str(bam_path))
                    _pysam.index(str(sorted_bam))
                    correction_bam = sorted_bam
                except Exception as _sort_exc:
                    print(
                        f"    WARNING: [{aligner_name}] sort failed ({_sort_exc}); "
                        "skipping correction",
                        file=sys.stderr,
                    )
                    if sorted_bam is not None and sorted_bam.exists():
                        sorted_bam.unlink(missing_ok=True)
                    continue
            try:
                _all_aligner_bams = list(per_aligner_bams.values())
                _run_correction(
                    bam_path=correction_bam,
                    output_dir=aligner_output_dir,
                    genome_path=genome_path,
                    annotation_path=annotation_path,
                    args=args,
                    aligner_bams=_all_aligner_bams,
                    reuse_pool_container=_pool_container,
                )
            except Exception as exc:
                print(
                    f"    WARNING: [{aligner_name}] correction failed: {exc}",
                    file=sys.stderr,
                )
            finally:
                if sorted_bam is not None:
                    sorted_bam.unlink(missing_ok=True)
                    sorted_bai = Path(str(sorted_bam) + '.bai')
                    sorted_bai.unlink(missing_ok=True)

            _canonical_after = _per_aligner_canonical_output(aligner_output_dir)
            if _canonical_after is not None:
                per_aligner_tsvs[aligner_name] = _canonical_after
                _stem = (sorted_bam.stem if sorted_bam is not None else bam_path.stem)
                for _sfx in ('.rectified', '.consensus', '.coord_sorted'):
                    if _stem.endswith(_sfx):
                        _stem = _stem[:-len(_sfx)]
                        break
                _corr_bam = aligner_output_dir / f"{_stem}.rectified_corrected_3end.bam"
                if _corr_bam.exists():
                    per_aligner_corrected_bams[aligner_name] = _corr_bam
            else:
                print(
                    f"    WARNING: [{aligner_name}] correction produced no output",
                    file=sys.stderr,
                )
    finally:
        if _pool_container:
            _pool_container[0].terminate()
            _pool_container[0].join()

    return per_aligner_tsvs, per_aligner_corrected_bams
def _combine_corrected_tsvs(
    samples: List[Dict[str, str]],
    output_dir: Path,
) -> Path:
    """
    Concatenate per-sample corrected TSVs into a single multi-sample TSV.

    Adds a `sample` column equal to sample_id so the analyze command can
    distinguish samples for DESeq2.
    """
    import pandas as pd

    combined_dir = output_dir / 'combined'
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_tsv = combined_dir / 'corrected_reads_combined.tsv'

    dfs = []
    missing = []
    for sample in samples:
        tsv_path = output_dir / sample['sample_id'] / 'corrected_reads.tsv'
        if not tsv_path.exists():
            missing.append(sample['sample_id'])
            print(f"  WARNING: {tsv_path} not found, skipping", file=sys.stderr)
            continue
        df = pd.read_csv(tsv_path, sep='\t')
        df['sample'] = sample['sample_id']
        dfs.append(df)
        print(f"  {sample['sample_id']}: {len(df):,} reads")

    if not dfs:
        raise RuntimeError(
            "No corrected TSVs found to combine. "
            "Ensure the correction step completed successfully."
        )

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(combined_tsv, sep='\t', index=False)
    print(f"  Combined → {combined_tsv}  ({len(dfs)} samples, {len(combined):,} reads)")

    if missing:
        print(
            f"  WARNING: {len(missing)} samples missing corrected TSV: {', '.join(missing)}",
            file=sys.stderr,
        )

    return combined_tsv
def _run_analysis(
    corrected_tsv: Path,
    output_dir: Path,
    genome_path: Optional[Path],
    annotation_path: Optional[Path],
    args,
    n_samples: int = 1,
) -> None:
    """Run the analyze command on a (possibly combined) corrected TSV."""
    from ..analyze_command import run_analyze

    # Only run DESeq2 when there are multiple samples
    run_deseq2 = n_samples > 1

    analyze_args = argparse.Namespace(
        input=corrected_tsv,
        output=output_dir,
        annotation=annotation_path,
        genome=genome_path,
        reference=getattr(args, 'reference', None),
        go_annotations=getattr(args, 'go_annotations', None),
        threads=getattr(args, 'threads', 4),
        # Clustering
        sample_column='sample',
        count_column=None,
        cluster_distance=25,
        min_reads=5,
        # Analysis flags
        run_deseq2=run_deseq2,
        run_motif=run_deseq2,   # motif discovery only makes sense with DESeq2 results
        sample_sets=None,
        # Filtering
        exclude_mito=True,
        include_mito=False,
        exclude_rdna=True,
        include_rdna=False,
        # Output options
        no_bedgraph=False,
        bedgraph_dir=None,
        no_genomic_distribution=False,
        # TSS clustering window (75 bp for DRS — 5' ends noisier than 3' ends)
        tss_cluster_distance=75,
        # Manifest mode (not used in single-file path)
        manifest=None,
        # Motif windows
        motif_upstream=100,
        motif_downstream=50,
    )

    exit_code = run_analyze(analyze_args)
    if exit_code != 0:
        print(f"\nAnalysis completed with warnings (exit code: {exit_code})")

def _run_junction_aggregation(
    bam_path: Path,
    genome_path: Path,
    annotation_path: Path,
    output_dir: Path,
    config: Optional[dict] = None,
) -> Optional[Path]:
    """
    Run junction aggregation with partial rescue.
    Only runs if annotation is GFF/GFF3 (has intron features).
    Returns path to junctions TSV, or None if skipped/failed.
    """
    # Handle compressed files: .gff.gz → stem is .gff, suffixes[-1] is .gz
    suffixes = [s.lower() for s in annotation_path.suffixes]
    has_gff = any(s in ('.gff', '.gff3') for s in suffixes)
    if not has_gff:
        print(
            "\n[Junctions] Skipping — requires GFF/GFF3 annotation "
            "(current annotation has no intron features)"
        )
        return None

    print("\n[Junctions] Aggregating splice junctions with partial rescue...")

    from ...aggregate.junctions import aggregate_junctions, merge_with_partial_evidence, export_junctions
    from ...splice.terminal_exon_refiner import load_splice_sites_from_gff, detect_partial_junction_crossings
    import pysam

    junctions_dir = output_dir / 'junctions'
    junctions_dir.mkdir(parents=True, exist_ok=True)
    junctions_tsv = junctions_dir / 'junctions.tsv'

    try:
        genome = {}
        fasta = pysam.FastaFile(str(genome_path))
        for chrom in fasta.references:
            genome[chrom] = fasta.fetch(chrom)
        fasta.close()

        junction_df = aggregate_junctions(
            bam_path=str(bam_path), genome=genome, min_reads=1,
        )
        print(f"  Found {len(junction_df)} junctions from CIGAR")

        splice_index = load_splice_sites_from_gff(str(annotation_path))
        partial_results = detect_partial_junction_crossings(
            bam_path=str(bam_path),
            genome=genome,
            splice_index=splice_index,
            min_clip_length=1,
            ambiguous_mode='proportional',
        )

        n_rescued = partial_results['stats']['reads_rescued_as_spliced']
        n_ambiguous = partial_results['stats']['reads_ambiguous']
        print(f"  Rescued {n_rescued} partial crossings ({n_ambiguous} ambiguous)")

        junction_df = merge_with_partial_evidence(
            junction_df, partial_results, ambiguous_mode='proportional',
        )
        export_junctions(junction_df, str(junctions_tsv), format='tsv')
        print(f"  Junctions written to {junctions_tsv}")

        # Provenance
        if junctions_tsv.exists():
            from ....utils.provenance import init_provenance
            prov = init_provenance(output_dir, description="RECTIFY junction aggregation", config=config)
            prov.add_output_file(
                junctions_tsv,
                source_files=[bam_path],
                metadata={
                    'n_junctions': len(junction_df),
                    'n_rescued': n_rescued,
                    'n_ambiguous': n_ambiguous,
                },
            )
            prov.save()

        return junctions_tsv

    except Exception as e:
        print(f"\nWarning: Junction aggregation failed: {e}", file=sys.stderr)
        print("Continuing without junction output...")
        return None


# ──────────────────────────────────────────────────────────────────────────────
# Browser pack — the LAST step of run-all, and the only optional one
# ──────────────────────────────────────────────────────────────────────────────

# Reads sampled per library for the QC pass. The selection is a hash on read_id,
# so this is a target, not a cap on the file — see qc_command's SAMPLING note.
_BROWSER_PACK_TARGET_READS = 200_000


def add_browser_pack_args(parser) -> None:
    """Attach ``--no-browser-pack`` to the ``run-all`` parser.

    Called from ``rectify.cli.create_parser`` against the parser that
    ``create_run_parser`` registered (that function does not return it).  A
    ``None`` parser is tolerated so a rename of the subcommand degrades to a
    missing flag rather than a CLI-wide crash.
    """
    if parser is None:
        return
    parser.add_argument(
        '--no-browser-pack',
        action='store_true',
        default=False,
        dest='no_browser_pack',
        help='Skip the final browser-pack step (per-sample `rectify qc` + '
             'analysis.json). The step is fail-soft and runs by default; use this '
             'to save the QC BAM pass when the browser bundle is not wanted.',
    )


def _resolve_qc_bam(
    sample_dir: Path,
    sample_id: str,
    fallback: Optional[Path] = None,
) -> Optional[Path]:
    """Pick the BAM whose geometry the QC block describes.

    Preference order: the pre-correction multi-aligner BAM (its read lengths,
    qualities and CIGARs are the alignment being characterised), then the final
    rectified BAM, then whatever the caller had in hand.
    """
    for cand in (
        _multialigned_bam_path(sample_id, sample_dir),
        sample_dir / f"{sample_id}.rectified.bam",
    ):
        try:
            if cand.exists() and cand.stat().st_size > 0:
                return cand
        except OSError:
            continue
    if fallback is not None:
        try:
            p = Path(fallback)
            if p.exists() and p.suffix.lower() == '.bam' and p.stat().st_size > 0:
                return p
        except OSError:
            pass
    return None


def _run_browser_pack(
    analysis_dir: Path,
    samples: List[Dict[str, object]],
    args,
    annotation_path: Optional[Path] = None,
    modality: str = '',
) -> Optional[Path]:
    """Per-sample QC + pack ``analysis.json``.  Runs LAST; never fatal.

    ``samples`` is a list of ``{'sample_id': str, 'bam': path|None,
    'tsv': path|None}``.

    🔴 FAIL-SOFT IS THE POINT.  A QC or packing failure must never fail a
    12-hour ``run-all`` that already produced correct BAMs and tables, so every
    step here is wrapped, the failure is printed loudly, and the function returns
    ``None``.

    Idempotent: a ``qc_<sample>.json`` that already exists and is non-empty is
    reused, so a re-run regenerates ``analysis.json`` from what is on disk
    without repeating the BAM pass.

    Fail-open: a sample missing its BAM or TSV is skipped, and ``analysis.json``
    is still written with whatever sections are available — the packer renders a
    missing section as "not computed", never as an empty table.

    Returns the ``analysis.json`` path, or ``None`` on any failure/skip.
    """
    if getattr(args, 'no_browser_pack', False):
        print("\n[browser-pack] skipped (--no-browser-pack)")
        return None

    try:
        import json as _json

        analysis_dir = Path(analysis_dir)
        qc_dir = analysis_dir / 'qc'
        qc_dir.mkdir(parents=True, exist_ok=True)

        print("\n[browser-pack] Per-library QC + analysis.json ...")
        print("-" * 50)

        # Protein-coding gene set (genes with >=1 CDS), derived once and shared.
        coding = None
        if annotation_path:
            try:
                from ..qc_command import coding_genes_from_gff
                _sfx = [s.lower() for s in Path(annotation_path).suffixes]
                if any(s in ('.gff', '.gff3') for s in _sfx):
                    coding = coding_genes_from_gff(annotation_path)
                    print(f"  coding-gene set: {len(coding):,} genes with >=1 CDS")
            except Exception as e:
                print(f"  WARNING: coding-gene set unavailable ({e}); "
                      f"the coding fraction will render as not computed",
                      file=sys.stderr)
                coding = None

        for s in samples:
            sid = str(s.get('sample_id') or '')
            if not sid:
                continue
            out_json = qc_dir / f"qc_{sid}.json"
            try:
                if out_json.exists() and out_json.stat().st_size > 0:
                    print(f"  [{sid}] QC already present — reusing {out_json.name}")
                    continue
                bam = s.get('bam')
                tsv = s.get('tsv')
                if not bam or not tsv or not Path(str(bam)).exists() \
                        or not Path(str(tsv)).exists():
                    print(f"  [{sid}] QC skipped — needs both a BAM and a corrected "
                          f"TSV (bam={bam}, tsv={tsv})", file=sys.stderr)
                    continue
                from ..qc_command import compute_qc
                block = compute_qc(
                    bam_path=str(bam),
                    tsv_path=str(tsv),
                    sample=sid,
                    modality=modality or '',
                    coding=coding,
                    target_reads=_BROWSER_PACK_TARGET_READS,
                )
                with open(out_json, 'w') as fh:
                    _json.dump(block, fh, indent=1)
                print(f"  [{sid}] QC -> {out_json.name}  "
                      f"err={block['error']['rate_pct']}%  "
                      f"N50={block['read_len']['n50']}")
            except Exception as e:
                print(f"  WARNING: QC failed for {sid} (non-fatal): {e}",
                      file=sys.stderr)

        from ..browser_pack_command import run_browser_pack
        out_path = analysis_dir / 'analysis.json'
        pack_args = argparse.Namespace(
            qc_dir=str(qc_dir),
            analyze_dir=str(analysis_dir),
            gene_index=getattr(args, 'gene_index', None),
            library_depth=str(analysis_dir / 'library_depth.json'),
            samples=None,
            pressure=None,
            top_n=getattr(args, 'browser_pack_top_n', 400),
            dataset=None,
            modality=modality or '',
            rectify_version=None,
            out=str(out_path),
        )
        rc = run_browser_pack(pack_args)
        if rc != 0 or not out_path.exists():
            print("WARNING: browser pack produced no analysis.json (non-fatal)",
                  file=sys.stderr)
            return None
        print(f"  Browser bundle: {out_path}")
        return out_path

    except Exception as e:
        print(f"WARNING: browser pack failed (non-fatal, pipeline outputs are "
              f"unaffected): {e}", file=sys.stderr)
        return None
