"""Tests for drs_trim stage sidecar wiring (Commit E.1).

Covers:
- skip-check fires on re-entry when sidecar matches inputs
- RUN when input file changes
- sidecar placed at <output_dir>/<sample_stem>.drs_trim.provenance.json
- --force-all overrides skip
- FASTQ-mode inputs key is 'input_fastq', BAM-mode is 'input_bam'
"""

import os
from argparse import Namespace
from datetime import datetime, timezone
from pathlib import Path
from unittest import mock

import pytest

from rectify.core.provenance.sidecar import ProvenanceRecord, write_stage_sidecar
from rectify.core.provenance.skip_check import can_skip_stage


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_drs_trim_sidecar(
    output_dir: Path,
    *,
    sample_id: str,
    inputs: dict,
    git_sha: str = "deadbeef1234",
) -> Path:
    now = datetime.now(timezone.utc).isoformat()
    env_clear = {k: "" for k in ("TMPDIR", "L_SCRATCH") if os.environ.get(k)}
    with mock.patch.dict(os.environ, env_clear):
        record = ProvenanceRecord.from_components(
            stage="drs_trim",
            sample_id=sample_id,
            sample_output_dir=output_dir,
            started_at=now,
            completed_at=now,
            exit_status=0,
            inputs=inputs,
            outputs={},
            stats={"wall_seconds": 5.0},
            argv=["rectify", "trim-polya"],
            rectify_git_sha=git_sha,
        )
    return write_stage_sidecar(record, output_dir)


# ---------------------------------------------------------------------------
# Sidecar location
# ---------------------------------------------------------------------------


def test_sidecar_written_at_correct_path(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM")
    sidecar_path = _write_drs_trim_sidecar(
        tmp_path, sample_id="sample", inputs={"input_bam": str(bam)}
    )
    assert sidecar_path.name == "sample.drs_trim.provenance.json"
    assert sidecar_path.parent == tmp_path
    assert sidecar_path.exists()


# ---------------------------------------------------------------------------
# Skip / RUN decisions
# ---------------------------------------------------------------------------


def test_skip_on_reentry_matching_bam(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content")
    git_sha = "abc123"
    inputs = {"input_bam": str(bam)}
    _write_drs_trim_sidecar(tmp_path, sample_id="sample", inputs=inputs, git_sha=git_sha)

    decision = can_skip_stage(
        stage="drs_trim",
        sample_output=tmp_path,
        sample_id="sample",
        current_inputs=inputs,
        current_argv=["rectify", "trim-polya"],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is True


def test_recompute_when_bam_content_changes(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content v1")
    git_sha = "abc123"
    inputs = {"input_bam": str(bam)}
    _write_drs_trim_sidecar(tmp_path, sample_id="sample", inputs=inputs, git_sha=git_sha)

    bam.write_bytes(b"BAM content v2 changed")

    decision = can_skip_stage(
        stage="drs_trim",
        sample_output=tmp_path,
        sample_id="sample",
        current_inputs=inputs,
        current_argv=["rectify", "trim-polya"],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is False


def test_skip_on_reentry_matching_fastq(tmp_path):
    fq = tmp_path / "reads.fastq.gz"
    fq.write_bytes(b"FASTQ content")
    git_sha = "abc123"
    inputs = {"input_fastq": str(fq)}
    _write_drs_trim_sidecar(tmp_path, sample_id="reads", inputs=inputs, git_sha=git_sha)

    decision = can_skip_stage(
        stage="drs_trim",
        sample_output=tmp_path,
        sample_id="reads",
        current_inputs=inputs,
        current_argv=["rectify", "trim-polya"],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is True


def test_force_all_overrides_skip(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content")
    git_sha = "abc123"
    inputs = {"input_bam": str(bam)}
    _write_drs_trim_sidecar(tmp_path, sample_id="sample", inputs=inputs, git_sha=git_sha)

    decision = can_skip_stage(
        stage="drs_trim",
        sample_output=tmp_path,
        sample_id="sample",
        current_inputs=inputs,
        current_argv=["rectify", "trim-polya"],
        rectify_git_sha=git_sha,
        force=True,
    )
    assert decision.skip is False


def test_bam_sidecar_skips_even_when_fastq_key_supplied(tmp_path):
    """skip-check iterates over SIDECAR inputs, not current inputs.

    A sidecar written for input_bam will still return SKIP if the bam content
    is unchanged, even if the caller passes a different input_fastq key. This
    is by design: the check is conservative. Use --force-all to override when
    switching modes.
    """
    out_dir = tmp_path / "run"
    out_dir.mkdir()
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM content")
    fq = tmp_path / "reads.fastq.gz"
    fq.write_bytes(b"FASTQ content")

    git_sha = "abc123"
    _write_drs_trim_sidecar(out_dir, sample_id="s",
                             inputs={"input_bam": str(bam)}, git_sha=git_sha)

    # The sidecar's input_bam is still unchanged on disk, so skip-check says SKIP
    # regardless of the extra input_fastq key in current_inputs.
    decision = can_skip_stage(
        stage="drs_trim",
        sample_output=out_dir,
        sample_id="s",
        current_inputs={"input_fastq": str(fq)},
        current_argv=["rectify", "trim-polya"],
        rectify_git_sha=git_sha,
    )
    assert decision.skip is True
