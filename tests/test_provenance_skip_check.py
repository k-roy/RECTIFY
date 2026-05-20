"""Tests for rectify.core.provenance.skip_check.can_skip_stage.

Tests each of the 4 RUN reasons + the SKIP case using parametrized table rows.
Also tests --force-all, --force-stage, and --accept-prior-provenance paths.
"""

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from unittest import mock

import pytest

from rectify.core.provenance.hashing import normalized_config_hash
from rectify.core.provenance.path_resolver import tokenize_argv_paths
from rectify.core.provenance.sidecar import (
    ProvenanceRecord,
    write_stage_sidecar,
)
from rectify.core.provenance.skip_check import SkipDecision, can_skip_stage


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------


def _write_valid_sidecar(
    tmp_path: Path,
    sample_dir: Path,
    *,
    stage: str = "correct",
    sample_id: str = "S1",
    git_sha: str = "abc1234abc1234",
    exit_status: int = 0,
    argv: list = None,
) -> tuple:
    """Write a valid sidecar and return (record, sidecar_path).

    $TMPDIR/$L_SCRATCH are cleared so that pytest's tmp_path (which lives
    under $TMPDIR on macOS) doesn't trigger the ephemeral-env-var guard in
    PortablePath.from_path for declared non-transient inputs.
    """
    argv = argv or ["rectify", "correct", "--sample", sample_id]

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")
    output_bam = sample_dir / "corrected.bam"
    output_bam.write_bytes(b"output_data")

    now = datetime.now(timezone.utc).isoformat()
    env_clear = {k: "" for k in ("TMPDIR", "L_SCRATCH") if os.environ.get(k)}
    with mock.patch.dict(os.environ, env_clear):
        record = ProvenanceRecord.from_components(
            stage=stage,
            sample_id=sample_id,
            sample_output_dir=sample_dir,
            started_at=now,
            completed_at=now,
            exit_status=exit_status,
            inputs={"bam": input_bam},
            outputs={"corrected_bam": output_bam},
            stats={"wall_seconds": 5.0},
            argv=argv,
            rectify_git_sha=git_sha,
            rectify_version="0.9.0",
        )
    sidecar_path = write_stage_sidecar(record, sample_dir)
    return record, sidecar_path


# ---------------------------------------------------------------------------
# SKIP case
# ---------------------------------------------------------------------------


def test_skip_when_everything_matches(tmp_path):
    """All conditions met → SKIP."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    git_sha = "abc1234abc1234"
    argv = ["rectify", "correct", "--sample", "S1"]

    record, _ = _write_valid_sidecar(
        tmp_path, sample_dir,
        git_sha=git_sha,
        argv=argv,
    )

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
    )
    assert decision.skip is True
    assert decision.prior is not None
    assert "prior valid run" in decision.reason.lower() or "completed_at" in decision.reason


# ---------------------------------------------------------------------------
# RUN reasons
# ---------------------------------------------------------------------------


def test_run_no_sidecar(tmp_path):
    """No prior sidecar → RUN."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"x")

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=["rectify", "correct"],
        rectify_git_sha="abc123",
    )
    assert decision.skip is False
    assert "no prior sidecar" in decision.reason.lower()


def test_run_prior_exit_status_nonzero(tmp_path):
    """Prior run exited non-zero → RUN."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    git_sha = "abc1234"
    argv = ["rectify", "correct", "--sample", "S1"]

    _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=argv, exit_status=1)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False
    assert "exit" in decision.reason.lower()


def test_run_input_sha256_changed(tmp_path):
    """Input sha256 changed since prior run → RUN."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"original_data")

    git_sha = "abc1234"
    argv = ["rectify", "correct", "--sample", "S1"]
    _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=argv)

    # Now modify the input file.
    input_bam.write_bytes(b"MODIFIED_data_completely_different")

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False
    assert "sha256" in decision.reason.lower() or "input" in decision.reason.lower()


def test_run_output_missing(tmp_path):
    """Declared output missing from disk → RUN."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    git_sha = "abc1234"
    argv = ["rectify", "correct", "--sample", "S1"]
    _, sidecar_path = _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=argv)

    # Delete the output file that the sidecar declared.
    output_bam = sample_dir / "corrected.bam"
    output_bam.unlink()

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False
    # Reason should mention missing output or sha256 or cannot resolve
    assert any(kw in decision.reason.lower() for kw in ("output", "missing", "sha256", "resolve"))


def test_run_git_sha_mismatch(tmp_path):
    """Git sha changed → RUN (default conservative behavior)."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    prior_sha = "aaaaaaa"
    current_sha = "bbbbbbb"
    argv = ["rectify", "correct", "--sample", "S1"]

    _write_valid_sidecar(tmp_path, sample_dir, git_sha=prior_sha, argv=argv)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=current_sha,
    )
    assert decision.skip is False
    assert "git_sha" in decision.reason.lower() or "sha" in decision.reason.lower()


def test_run_argv_changed(tmp_path):
    """Argv changed → RUN."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    git_sha = "abc1234"
    prior_argv = ["rectify", "correct", "--sample", "S1", "--mode", "drs"]
    current_argv = ["rectify", "correct", "--sample", "S1", "--mode", "cdna"]

    _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=prior_argv)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=current_argv,
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False
    assert "argv" in decision.reason.lower() or "config" in decision.reason.lower()


# ---------------------------------------------------------------------------
# --force flags
# ---------------------------------------------------------------------------


def test_force_all_overrides_skip(tmp_path):
    """force=True → RUN even if sidecar is valid."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    git_sha = "abc1234"
    argv = ["rectify", "correct", "--sample", "S1"]
    _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=argv)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
        force=True,
    )
    assert decision.skip is False
    assert "force-all" in decision.reason.lower()


def test_force_stage_overrides_skip(tmp_path):
    """force_stages={"correct"} → RUN for correct, even if sidecar is valid."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    git_sha = "abc1234"
    argv = ["rectify", "correct", "--sample", "S1"]
    _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=argv)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
        force_stages={"correct"},
    )
    assert decision.skip is False
    assert "force-stage" in decision.reason.lower()


def test_force_stage_only_affects_named_stage(tmp_path):
    """force_stages={"analyze"} does not force "correct" stage."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    git_sha = "abc1234"
    argv = ["rectify", "correct", "--sample", "S1"]
    _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=argv)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
        force_stages={"analyze"},  # different stage
    )
    assert decision.skip is True  # Should still skip "correct"


def test_accept_prior_provenance_skips_sha_mismatch(tmp_path):
    """accept_prior_provenance=True → SKIP even when git_sha changed."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"input_data")

    prior_sha = "aaaaaaa"
    current_sha = "bbbbbbb"
    argv = ["rectify", "correct", "--sample", "S1"]

    _write_valid_sidecar(tmp_path, sample_dir, git_sha=prior_sha, argv=argv)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=current_sha,
        accept_prior_provenance=True,
    )
    assert decision.skip is True
    assert decision.prior is not None


# ---------------------------------------------------------------------------
# SkipDecision constructors and repr
# ---------------------------------------------------------------------------


def test_skip_decision_run_has_no_prior():
    d = SkipDecision.RUN(reason="test run")
    assert d.skip is False
    assert d.prior is None
    assert "RUN" in repr(d)


def test_skip_decision_skip_carries_prior(tmp_path):
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"x")

    git_sha = "aaa"
    argv = ["rectify", "correct"]
    record, _ = _write_valid_sidecar(tmp_path, sample_dir, git_sha=git_sha, argv=argv)

    decision = can_skip_stage(
        stage="correct",
        sample_output=sample_dir,
        sample_id="S1",
        current_inputs={"bam": input_bam},
        current_argv=argv,
        rectify_git_sha=git_sha,
    )
    assert decision.skip is True
    assert isinstance(decision.prior, ProvenanceRecord)
    assert "SKIP" in repr(decision)
