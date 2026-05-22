"""Tests for analyze-stage sidecar wiring (Commit D).

Covers:
- _collect_analyze_inputs builds the right key/value dict for manifest vs
  single-input mode
- A valid analyze sidecar triggers SKIP on re-entry
- Changed inputs (TSV path or manifest content) trigger RUN
"""

import os
from argparse import Namespace
from datetime import datetime, timezone
from pathlib import Path
from unittest import mock

import pandas as pd
import pytest

from rectify.core.provenance.sidecar import ProvenanceRecord, write_stage_sidecar
from rectify.core.provenance.skip_check import can_skip_stage
from rectify.core.commands.analyze_command import _collect_analyze_inputs


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_analyze_sidecar(
    tmp_path: Path,
    sample_dir: Path,
    *,
    stage: str = "analyze",
    sample_id: str = "run1",
    inputs: dict,
    git_sha: str = "deadbeef1234",
    stage_subtype=None,
) -> Path:
    """Write a completed analyze sidecar and return its path."""
    now = datetime.now(timezone.utc).isoformat()
    env_clear = {k: "" for k in ("TMPDIR", "L_SCRATCH") if os.environ.get(k)}
    with mock.patch.dict(os.environ, env_clear):
        record = ProvenanceRecord.from_components(
            stage=stage,
            stage_subtype=stage_subtype,
            sample_id=sample_id,
            sample_output_dir=sample_dir,
            started_at=now,
            completed_at=now,
            exit_status=0,
            inputs=inputs,
            outputs={},
            stats={"wall_seconds": 10.0},
            argv=["rectify", "analyze", "--output", str(sample_dir)],
            rectify_git_sha=git_sha,
            rectify_version="0.9.0",
        )
    return write_stage_sidecar(record, sample_dir)


# ---------------------------------------------------------------------------
# _collect_analyze_inputs
# ---------------------------------------------------------------------------


def test_collect_analyze_inputs_single_mode(tmp_path):
    """Non-manifest mode: 'input' key present, no tsv_* keys."""
    tsv = tmp_path / "sample.tsv"
    tsv.write_text("chrom\tpos\n")
    args = Namespace(input=str(tsv), manifest=None, annotation=None)
    result = _collect_analyze_inputs(args)
    assert result == {"input": str(tsv)}


def test_collect_analyze_inputs_manifest_mode(tmp_path):
    """Manifest mode: manifest key + one tsv_<id> key per sample."""
    tsv_s1 = tmp_path / "s1.tsv"
    tsv_s1.write_text("chrom\tpos\n")
    tsv_s2 = tmp_path / "s2.tsv"
    tsv_s2.write_text("chrom\tpos\n")

    manifest = tmp_path / "manifest.tsv"
    pd.DataFrame(
        [{"sample_id": "s1", "path": str(tsv_s1)}, {"sample_id": "s2", "path": str(tsv_s2)}]
    ).to_csv(manifest, sep="\t", index=False)

    args = Namespace(input=None, manifest=str(manifest), annotation=None)
    result = _collect_analyze_inputs(args)
    assert result["manifest"] == str(manifest)
    assert result["tsv_s1"] == str(tsv_s1)
    assert result["tsv_s2"] == str(tsv_s2)
    assert len(result) == 3


def test_collect_analyze_inputs_includes_annotation(tmp_path):
    """annotation key is included when present."""
    tsv = tmp_path / "sample.tsv"
    tsv.write_text("chrom\tpos\n")
    ann = tmp_path / "anno.gff"
    ann.write_text("")
    args = Namespace(input=str(tsv), manifest=None, annotation=str(ann))
    result = _collect_analyze_inputs(args)
    assert result["annotation"] == str(ann)
    assert "input" in result


def test_collect_analyze_inputs_bad_manifest_is_graceful(tmp_path):
    """Malformed manifest falls back to just the manifest path (no crash)."""
    manifest = tmp_path / "bad.tsv"
    manifest.write_text("not\ta\tvalid\tmanifest\n")
    args = Namespace(input=None, manifest=str(manifest), annotation=None)
    result = _collect_analyze_inputs(args)
    assert result["manifest"] == str(manifest)
    # no tsv_* keys but no exception either
    assert all(k == "manifest" for k in result)


# ---------------------------------------------------------------------------
# Skip-check behaviour for analyze stage
# ---------------------------------------------------------------------------


def test_skip_on_reentry_matching_inputs(tmp_path):
    """Existing analyze sidecar with matching inputs → SKIP."""
    sample_dir = tmp_path / "run1"
    sample_dir.mkdir()
    tsv = tmp_path / "sample.tsv"
    tsv.write_bytes(b"chrom\tpos\n")

    git_sha = "deadbeef1234"
    inputs = {"input": str(tsv)}
    _write_analyze_sidecar(tmp_path, sample_dir, inputs=inputs, git_sha=git_sha)

    decision = can_skip_stage(
        stage="analyze",
        sample_output=sample_dir,
        sample_id="run1",
        current_inputs=inputs,
        current_argv=["rectify", "analyze", "--output", str(sample_dir)],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is True


def test_recompute_when_input_file_changes(tmp_path):
    """After sidecar written, changing input file content → RUN."""
    sample_dir = tmp_path / "run1"
    sample_dir.mkdir()
    tsv = tmp_path / "sample.tsv"
    tsv.write_bytes(b"chrom\tpos\n")

    git_sha = "deadbeef1234"
    inputs = {"input": str(tsv)}
    _write_analyze_sidecar(tmp_path, sample_dir, inputs=inputs, git_sha=git_sha)

    # Simulate the input changing
    tsv.write_bytes(b"chrom\tpos\nnew\tdata\n")

    decision = can_skip_stage(
        stage="analyze",
        sample_output=sample_dir,
        sample_id="run1",
        current_inputs=inputs,
        current_argv=["rectify", "analyze", "--output", str(sample_dir)],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False


def test_recompute_when_manifest_sample_tsv_changes(tmp_path):
    """In manifest mode: new per-sample TSV path in inputs dict → RUN."""
    sample_dir = tmp_path / "run1"
    sample_dir.mkdir()
    tsv_old = tmp_path / "old.tsv"
    tsv_old.write_bytes(b"chrom\tpos\n")
    tsv_new = tmp_path / "new.tsv"
    tsv_new.write_bytes(b"chrom\tpos\n")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_bytes(b"sample_id\tpath\ns1\told.tsv\n")

    git_sha = "deadbeef1234"
    old_inputs = {"manifest": str(manifest), "tsv_s1": str(tsv_old)}
    _write_analyze_sidecar(tmp_path, sample_dir, inputs=old_inputs, git_sha=git_sha)

    # Now inputs reflect a different per-sample TSV
    new_inputs = {"manifest": str(manifest), "tsv_s1": str(tsv_new)}
    decision = can_skip_stage(
        stage="analyze",
        sample_output=sample_dir,
        sample_id="run1",
        current_inputs=new_inputs,
        current_argv=["rectify", "analyze", "--manifest", str(manifest)],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False


def test_sidecar_written_in_correct_location(tmp_path):
    """write_stage_sidecar places file at <sample_output>/<sample_id>.analyze.provenance.json."""
    sample_dir = tmp_path / "run1"
    sample_dir.mkdir()
    tsv = tmp_path / "sample.tsv"
    tsv.write_bytes(b"chrom\tpos\n")

    sidecar_path = _write_analyze_sidecar(
        tmp_path, sample_dir, inputs={"input": str(tsv)}
    )
    assert sidecar_path.name == "run1.analyze.provenance.json"
    assert sidecar_path.parent == sample_dir
    assert sidecar_path.exists()


def test_force_all_overrides_skip(tmp_path):
    """--force-all causes RUN even when sidecar would match."""
    sample_dir = tmp_path / "run1"
    sample_dir.mkdir()
    tsv = tmp_path / "sample.tsv"
    tsv.write_bytes(b"chrom\tpos\n")

    git_sha = "deadbeef1234"
    inputs = {"input": str(tsv)}
    _write_analyze_sidecar(tmp_path, sample_dir, inputs=inputs, git_sha=git_sha)

    decision = can_skip_stage(
        stage="analyze",
        sample_output=sample_dir,
        sample_id="run1",
        current_inputs=inputs,
        current_argv=["rectify", "analyze", "--output", str(sample_dir)],
        rectify_git_sha=git_sha,
        force=True,
    )
    assert decision.skip is False
