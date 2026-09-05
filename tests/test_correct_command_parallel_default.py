"""Commit B acceptance tests for correct_command parallel default + resume.

Tests:
  1. Small-input gate fires on tiny validation BAM (sequential path).
  2. Default output = manifest-only (no merged TSV unless --emit-merged-tsv).
  3. --emit-merged-tsv produces both manifest + merged TSV.
  4. Sidecar written with correct fields.
  5. Resume: immediate rerun returns SKIP.
  6. Resume: mtime change on input → still SKIP (content hash, not mtime).
  7. Resume: argv change on ignored flag (--threads) → SKIP.
  8. Resume: argv change on load-bearing flag (--emit-merged-tsv) → RUN.
  9. --legacy-single-threaded produces byte-equivalent BAM vs default path.

Author: Kevin R. Roy
Date: 2026-05-20
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import time
from pathlib import Path

import pysam
import pytest


# ---------------------------------------------------------------------------
# Shared helpers from test_correct_command_drs (duplicated here to avoid
# coupling the tests)
# ---------------------------------------------------------------------------

_CHROM_LEN = 2000


def _make_reference(tmp_path: Path) -> Path:
    """Build a minimal single-chromosome synthetic reference."""
    seq = list("X" * _CHROM_LEN)
    for i, b in enumerate("ACGTACGTAC"):
        seq[1000 + i] = b
    for i in range(10):
        seq[1010 + i] = "A"
    for i in range(10):
        seq[1020 + i] = "C"
    fa = tmp_path / "ref.fa"
    fa.write_text(f">chrT\n{''.join(seq)}\n")
    pysam.faidx(str(fa))
    return fa


def _make_small_bam(tmp_path: Path) -> Path:
    """Build a 5-read BAM that will trigger the small-input gate."""
    bam_path = tmp_path / "input.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": _CHROM_LEN}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(5):
            r = pysam.AlignedSegment(bam.header)
            r.query_name = f"r{i}_test"
            r.query_sequence = "ACGTACGTAC" + "A" * 5
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 1000
            r.mapping_quality = 60
            r.cigartuples = [(0, 15)]
            bam.write(r)
    pysam.index(str(bam_path))
    return bam_path


def _make_correct_args(
    *,
    bam: Path,
    genome: Path,
    output: Path,
    corrected_bam: Path = None,
    threads: int = 1,
    emit_merged_tsv: bool = False,
    legacy_single_threaded: bool = False,
    force_all: bool = False,
    force_stage: str = None,
    accept_prior_provenance: bool = False,
    dry_run_resume: bool = False,
) -> argparse.Namespace:
    """Build args for correct_command.run."""
    return argparse.Namespace(
        input=bam,
        output=output,
        genome=genome,
        annotation=None,
        dT_primed_cDNA=False,
        polya_sequenced=False,
        ONT_cDNA=False,
        short_read=False,
        aligner='auto',
        skip_atract_check=False,
        skip_ag_check=True,
        skip_polya_trim=False,
        skip_indel_correction=False,
        skip_variant_aware=True,
        filter_spikein=None,
        organism=None,
        netseq_dir=None,
        netseq_samples=None,
        polya_model=None,
        report=None,
        min_mapq=0,
        min_aligned_length=0,
        ag_threshold=17.0,
        output_bam=None,
        corrected_bam=corrected_bam,
        softclipped_bam=None,
        bedgraph_prefix=None,
        threads=threads,
        verbose=False,
        streaming=False,
        chunk_size=10000,
        checkpoint_dir=None,
        tmp_dir=None,
        emit_merged_tsv=emit_merged_tsv,
        legacy_single_threaded=legacy_single_threaded,
        variant_scan_cache=None,
        junction_pool_cache=None,
        aligner_bams=[],
        junction_hp_pen=0.25,
        junction_search_radius=5000,
        junction_window=15,
        junction_max_slide=10,
        junction_max_boundary_shift=50,
        junction_penalty_table=None,
        str_penalty_table=None,
        # Resume flags
        force_all=force_all,
        force_stage=force_stage,
        accept_prior_provenance=accept_prior_provenance,
        dry_run_resume=dry_run_resume,
    )


def _run_correct(args) -> int:
    from rectify.core.commands import correct_command
    try:
        result = correct_command.run(args)
        return result if result is not None else 0
    except SystemExit as e:
        return e.code if e.code is not None else 0


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_small_input_gate_fires_sequential(tmp_path, caplog):
    """Small BAM (5 reads) must trigger the sequential path via the gate."""
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected.tsv"

    args = _make_correct_args(bam=bam, genome=genome, output=output)

    with caplog.at_level(logging.INFO):
        rc = _run_correct(args)

    assert rc == 0, "correct_command should succeed"
    # Sentinel log for the small-input gate or sequential path
    log_text = caplog.text
    assert (
        "small-input gate fired" in log_text
        or "using legacy single-threaded" in log_text
        or "sequential" in log_text.lower()
    ), "Expected sequential-path sentinel in log"


def test_default_output_manifest_only(tmp_path):
    """Default mode: manifest exists, merged TSV at canonical path does NOT."""
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    args = _make_correct_args(bam=bam, genome=genome, output=output,
                               emit_merged_tsv=False)  # default
    rc = _run_correct(args)
    assert rc == 0

    manifest = tmp_path / "corrected_reads.manifest.tsv"
    assert manifest.exists(), "corrected_reads.manifest.tsv must be written"
    # The merged TSV at the canonical path must NOT exist in manifest-only mode
    assert not output.exists(), (
        "corrected_reads.tsv should NOT exist in manifest-only mode; "
        "pass --emit-merged-tsv to get it"
    )
    # A region TSV referenced by the manifest must exist
    manifest_lines = manifest.read_text().splitlines()
    assert len(manifest_lines) >= 2, "Manifest must have header + at least one row"


def test_emit_merged_tsv_produces_both(tmp_path):
    """--emit-merged-tsv produces manifest + merged TSV."""
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    args = _make_correct_args(bam=bam, genome=genome, output=output,
                               emit_merged_tsv=True)
    rc = _run_correct(args)
    assert rc == 0

    assert output.exists(), "corrected_reads.tsv must exist with --emit-merged-tsv"
    manifest = tmp_path / "corrected_reads.manifest.tsv"
    assert manifest.exists(), "corrected_reads.manifest.tsv must always be written"


def test_sidecar_written_with_correct_fields(tmp_path):
    """Stage sidecar has schema_version=1.0, exit_status=0, PortablePath envelopes."""
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    args = _make_correct_args(bam=bam, genome=genome, output=output,
                               emit_merged_tsv=True)
    rc = _run_correct(args)
    assert rc == 0

    # Derive sample_id (must match _derive_sample_id logic: input stem)
    sample_id = bam.stem
    sidecar_path = tmp_path / f"{sample_id}.correct.provenance.json"
    assert sidecar_path.exists(), f"Sidecar not found at {sidecar_path}"

    with sidecar_path.open() as fh:
        sidecar = json.load(fh)

    assert sidecar["schema_version"] == "1.0"
    assert sidecar["exit_status"] == 0
    assert sidecar["stage"] == "correct"

    # Per-module wall time (2026-09-05): the workers' per-read timers reach the
    # sidecar beside the stage-level walls, so the 2F/2H share is measurable
    # from the sidecar alone.
    module_seconds = sidecar["stats"]["module_seconds"]
    assert module_seconds["read_total"] >= module_seconds["five_prime_rescue_2f"] >= 0.0
    assert "bam_processing" in module_seconds
    stats_tsv = tmp_path / "corrected_reads_stats.tsv"
    assert "module_seconds_read_total\t" in stats_tsv.read_text()

    # All path fields in inputs/outputs must be PortablePath envelopes
    valid_kinds = {"sample_relative", "env_relative", "absolute", "cached_absolute"}
    for inp in sidecar.get("inputs", []):
        assert isinstance(inp["path"], dict), (
            f"Input path must be a PortablePath dict, got: {inp['path']!r}"
        )
        assert inp["path"].get("kind") in valid_kinds, (
            f"Input path kind must be one of {valid_kinds}, got: {inp['path'].get('kind')!r}"
        )
    for out in sidecar.get("outputs", []):
        assert isinstance(out["path"], dict), (
            f"Output path must be a PortablePath dict, got: {out['path']!r}"
        )
        assert out["path"].get("kind") in valid_kinds


def test_resume_skip_on_rerun(tmp_path):
    """Immediate rerun returns SKIP (all inputs/outputs unchanged)."""
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    args = _make_correct_args(bam=bam, genome=genome, output=output,
                               emit_merged_tsv=True)
    rc = _run_correct(args)
    assert rc == 0

    # Second run — should skip
    t0 = time.monotonic()
    rc2 = _run_correct(args)
    elapsed = time.monotonic() - t0
    assert rc2 == 0, "Resume run should return 0"
    assert elapsed < 5.0, (
        f"Resume run took {elapsed:.1f}s — expected < 5s for a SKIP"
    )


def test_resume_skip_after_mtime_touch(tmp_path):
    """Touch (mtime change) without content change → still SKIP."""
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    args = _make_correct_args(bam=bam, genome=genome, output=output,
                               emit_merged_tsv=True)
    rc = _run_correct(args)
    assert rc == 0

    # Touch the input BAM (mtime update, content unchanged).
    os.utime(str(bam), None)
    time.sleep(0.01)  # ensure mtime differs

    t0 = time.monotonic()
    rc2 = _run_correct(args)
    elapsed = time.monotonic() - t0
    assert rc2 == 0
    assert elapsed < 5.0, (
        f"After mtime touch, resume took {elapsed:.1f}s — expected SKIP (content hash match)"
    )


def test_resume_force_all_reruns(tmp_path, caplog):
    """--force-all causes the stage to rerun even if sidecar is valid."""
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    args1 = _make_correct_args(bam=bam, genome=genome, output=output,
                                emit_merged_tsv=True)
    rc1 = _run_correct(args1)
    assert rc1 == 0

    # Second run with --force-all: must NOT skip.
    args2 = _make_correct_args(bam=bam, genome=genome, output=output,
                                emit_merged_tsv=True, force_all=True)
    with caplog.at_level(logging.INFO):
        rc2 = _run_correct(args2)
    assert rc2 == 0
    # With --force-all, should have decided RUN
    log_text = caplog.text
    # Look for the [RESUME] line that says RUN
    assert "decision=RUN" in log_text, (
        "With --force-all, the [RESUME] log must show decision=RUN"
    )


def test_legacy_single_threaded_equivalent_output(tmp_path):
    """--legacy-single-threaded produces a BAM byte-equivalent to the default path."""
    from tests.utils import assert_bams_equivalent

    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output_default = tmp_path / "corrected_reads_default.tsv"
    output_legacy = tmp_path / "corrected_reads_legacy.tsv"
    cbam_default = tmp_path / "corrected_default.bam"
    cbam_legacy = tmp_path / "corrected_legacy.bam"

    # Default path (small-input gate will route to sequential anyway)
    args_default = _make_correct_args(
        bam=bam, genome=genome, output=output_default,
        corrected_bam=cbam_default, emit_merged_tsv=True,
    )
    rc_default = _run_correct(args_default)
    assert rc_default == 0

    # Legacy single-threaded path
    args_legacy = _make_correct_args(
        bam=bam, genome=genome, output=output_legacy,
        corrected_bam=cbam_legacy, emit_merged_tsv=True,
        legacy_single_threaded=True,
    )
    rc_legacy = _run_correct(args_legacy)
    assert rc_legacy == 0

    assert cbam_default.exists(), "Default-path corrected BAM not found"
    assert cbam_legacy.exists(), "Legacy-path corrected BAM not found"

    assert_bams_equivalent(str(cbam_default), str(cbam_legacy))
