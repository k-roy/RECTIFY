"""Tests for cdna_stage1 sidecar wiring (Commit E.2).

Covers:
- skip-check fires on re-entry when sidecar matches inputs
- RUN when input BAM changes
- sidecar placed at <out_dir>/<sample_id>.cdna_stage1.provenance.json
- --force-all overrides skip
"""

import os
from datetime import datetime, timezone
from pathlib import Path
from unittest import mock

import pytest

from rectify.core.provenance.sidecar import ProvenanceRecord, write_stage_sidecar
from rectify.core.provenance.skip_check import can_skip_stage


def _write_cdna_stage1_sidecar(
    out_dir: Path,
    *,
    sample_id: str,
    inputs: dict,
    git_sha: str = "deadbeef1234",
) -> Path:
    now = datetime.now(timezone.utc).isoformat()
    env_clear = {k: "" for k in ("TMPDIR", "L_SCRATCH") if os.environ.get(k)}
    with mock.patch.dict(os.environ, env_clear):
        record = ProvenanceRecord.from_components(
            stage="cdna_stage1",
            sample_id=sample_id,
            sample_output_dir=out_dir,
            started_at=now,
            completed_at=now,
            exit_status=0,
            inputs=inputs,
            outputs={},
            stats={"wall_seconds": 20.0},
            argv=["rectify", "correct-cdna"],
            rectify_git_sha=git_sha,
        )
    return write_stage_sidecar(record, out_dir)


def test_sidecar_at_correct_path(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM")
    out_dir = tmp_path / "sample"
    out_dir.mkdir()
    sidecar = _write_cdna_stage1_sidecar(
        out_dir, sample_id="sample", inputs={"input_bam": str(bam)}
    )
    assert sidecar.name == "sample.cdna_stage1.provenance.json"
    assert sidecar.parent == out_dir


def test_skip_on_reentry_matching_bam(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content v1")
    out_dir = tmp_path / "sample"
    out_dir.mkdir()
    git_sha = "abc123"
    inputs = {"input_bam": str(bam)}
    _write_cdna_stage1_sidecar(out_dir, sample_id="sample", inputs=inputs, git_sha=git_sha)

    decision = can_skip_stage(
        stage="cdna_stage1",
        sample_output=out_dir,
        sample_id="sample",
        current_inputs=inputs,
        current_argv=["rectify", "correct-cdna"],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is True


def test_recompute_when_bam_changes(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content v1")
    out_dir = tmp_path / "sample"
    out_dir.mkdir()
    git_sha = "abc123"
    inputs = {"input_bam": str(bam)}
    _write_cdna_stage1_sidecar(out_dir, sample_id="sample", inputs=inputs, git_sha=git_sha)

    bam.write_bytes(b"BAM content v2 updated")

    decision = can_skip_stage(
        stage="cdna_stage1",
        sample_output=out_dir,
        sample_id="sample",
        current_inputs=inputs,
        current_argv=["rectify", "correct-cdna"],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False


def test_force_all_overrides_skip(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content")
    out_dir = tmp_path / "sample"
    out_dir.mkdir()
    git_sha = "abc123"
    inputs = {"input_bam": str(bam)}
    _write_cdna_stage1_sidecar(out_dir, sample_id="sample", inputs=inputs, git_sha=git_sha)

    decision = can_skip_stage(
        stage="cdna_stage1",
        sample_output=out_dir,
        sample_id="sample",
        current_inputs=inputs,
        current_argv=["rectify", "correct-cdna"],
        rectify_git_sha=git_sha,
        force=True,
    )
    assert decision.skip is False


def test_no_sidecar_returns_run(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content")
    out_dir = tmp_path / "sample"
    out_dir.mkdir()

    decision = can_skip_stage(
        stage="cdna_stage1",
        sample_output=out_dir,
        sample_id="sample",
        current_inputs={"input_bam": str(bam)},
        current_argv=["rectify", "correct-cdna"],
        rectify_git_sha="abc123",
    )
    assert decision.skip is False
    assert "no prior sidecar" in decision.reason
