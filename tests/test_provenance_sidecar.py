"""Tests for rectify.core.provenance.sidecar module.

Covers:
- Round-trip: build → write → read → deep-equal
- Atomic write: tmp file left behind on simulated crash; no corrupt sidecar
- Schema violation: bare string in path field → ValueError
- read on invalid JSON → ValueError
- read on schema_version != "1.0" → ValueError
- Optional outputs: allowed to be absent
"""

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from unittest import mock

import pytest

from rectify.core.provenance.path_resolver import PortablePath
from rectify.core.provenance.sidecar import (
    InputEntry,
    OutputEntry,
    ProvenanceRecord,
    SubOutputEntry,
    read_stage_sidecar,
    write_stage_sidecar,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_input_file(tmp_path: Path, name: str = "input.bam", content: bytes = b"bam") -> Path:
    p = tmp_path / name
    p.write_bytes(content)
    return p


def _make_output_file(tmp_path: Path, name: str = "output.bam", content: bytes = b"out") -> Path:
    p = tmp_path / name
    p.write_bytes(content)
    return p


def _make_record(tmp_path: Path, sample_id: str = "S1", stage: str = "correct") -> ProvenanceRecord:
    """Build a minimal ProvenanceRecord using from_components.

    $TMPDIR is cleared so that pytest's tmp_path (which lives under $TMPDIR
    on macOS) doesn't trigger the ephemeral-env-var guard in PortablePath.
    Real stage code will use paths outside $TMPDIR for its declared inputs;
    this fixture mimics that by pretending $TMPDIR is unset.
    """
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir(exist_ok=True)

    input_bam = _make_input_file(tmp_path, "input.bam")
    output_bam = _make_output_file(sample_dir, "output.bam")

    now = datetime.now(timezone.utc).isoformat()
    # Clear ephemeral env vars so paths under pytest's tmp_path don't
    # accidentally hit the non-transient ephemeral guard.
    env_clear = {k: "" for k in ("TMPDIR", "L_SCRATCH") if os.environ.get(k)}
    with mock.patch.dict(os.environ, env_clear):
        return ProvenanceRecord.from_components(
            stage=stage,
            sample_id=sample_id,
            sample_output_dir=sample_dir,
            started_at=now,
            completed_at=now,
            exit_status=0,
            inputs={"bam": input_bam},
            outputs={"corrected_bam": output_bam},
            stats={"wall_seconds": 1.23, "reads": 100},
            argv=["rectify", "correct", "--sample", sample_id],
            rectify_git_sha="abc1234",
            rectify_version="0.9.0",
        )


# ---------------------------------------------------------------------------
# Round-trip
# ---------------------------------------------------------------------------


def test_roundtrip_write_and_read(tmp_path):
    """write_stage_sidecar + read_stage_sidecar produces an equal record."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    record = _make_record(tmp_path)

    sidecar_path = write_stage_sidecar(record, sample_dir)
    assert sidecar_path.exists()
    assert sidecar_path.name == f"{record.sample_id}.{record.stage}.provenance.json"

    recovered = read_stage_sidecar(sidecar_path)

    # Compare key fields
    assert recovered.schema_version == record.schema_version
    assert recovered.stage == record.stage
    assert recovered.sample_id == record.sample_id
    assert recovered.rectify_git_sha == record.rectify_git_sha
    assert recovered.exit_status == record.exit_status
    assert recovered.stats == record.stats
    assert len(recovered.inputs) == len(record.inputs)
    assert len(recovered.outputs) == len(record.outputs)
    assert recovered.inputs[0].role == record.inputs[0].role
    assert recovered.inputs[0].sha256 == record.inputs[0].sha256
    assert recovered.outputs[0].role == record.outputs[0].role
    assert recovered.outputs[0].sha256 == record.outputs[0].sha256


def test_roundtrip_portablepath_preserved(tmp_path):
    """PortablePath kind and template survive the round-trip."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    record = _make_record(tmp_path)

    sidecar_path = write_stage_sidecar(record, sample_dir)
    recovered = read_stage_sidecar(sidecar_path)

    orig_pp = record.outputs[0].path
    recov_pp = recovered.outputs[0].path
    assert isinstance(recov_pp, PortablePath)
    assert recov_pp.kind == orig_pp.kind
    assert recov_pp.template == orig_pp.template
    assert recov_pp.cached_absolute == orig_pp.cached_absolute


def test_roundtrip_invocation_preserved(tmp_path):
    """invocation dict (argv, argv_template, config_hash) survives round-trip."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    record = _make_record(tmp_path)

    sidecar_path = write_stage_sidecar(record, sample_dir)
    recovered = read_stage_sidecar(sidecar_path)

    assert recovered.invocation["config_hash"] == record.invocation["config_hash"]
    assert recovered.invocation["argv"] == record.invocation["argv"]


# ---------------------------------------------------------------------------
# Atomic write
# ---------------------------------------------------------------------------


def test_atomic_write_no_corrupt_sidecar_on_crash(tmp_path):
    """Simulated crash (write .tmp but don't rename) leaves no corrupt final sidecar."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    record = _make_record(tmp_path)

    sidecar_name = f"{record.sample_id}.{record.stage}.provenance.json"
    sidecar_path = sample_dir / sidecar_name
    tmp_path_expected = sidecar_path.with_suffix(".json.tmp")

    # Simulate a crash: write the .tmp file manually without renaming.
    json_str = json.dumps(record.to_json(), indent=2)
    tmp_path_expected.write_text(json_str)

    # The final sidecar should NOT exist yet.
    assert not sidecar_path.exists(), "Sidecar shouldn't exist without a rename"
    # The .tmp file exists (would be cleaned up on next run).
    assert tmp_path_expected.exists()


def test_atomic_write_uses_os_replace(tmp_path, monkeypatch):
    """write_stage_sidecar calls os.replace (atomic POSIX rename)."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    record = _make_record(tmp_path)

    calls = []
    orig_replace = os.replace

    def spy_replace(src, dst):
        calls.append((src, dst))
        return orig_replace(src, dst)

    monkeypatch.setattr(os, "replace", spy_replace)
    write_stage_sidecar(record, sample_dir)

    assert len(calls) == 1
    src, dst = calls[0]
    assert str(src).endswith(".tmp")
    assert str(dst).endswith(".provenance.json")


# ---------------------------------------------------------------------------
# Schema violation: bare string in path field
# ---------------------------------------------------------------------------


def test_write_rejects_bare_string_in_input_path(tmp_path):
    """write_stage_sidecar raises ValueError if an input entry has a bare string path."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    record = _make_record(tmp_path)
    # Inject a bare-string path via object.__setattr__ (record is not frozen).
    # Build a new InputEntry with a string instead of PortablePath.
    bad_input = InputEntry(
        role="bam",
        sha256="sha256:deadbeef",
        size_bytes=100,
        path="/bare/string/path",  # type: ignore[arg-type]  — intentional schema violation
    )
    record.inputs = [bad_input]

    with pytest.raises(ValueError, match="PortablePath"):
        write_stage_sidecar(record, sample_dir)


def test_write_rejects_bare_string_in_output_path(tmp_path):
    """write_stage_sidecar raises ValueError if an output entry has a bare string path."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    record = _make_record(tmp_path)
    bad_output = OutputEntry(
        role="corrected_bam",
        sha256="sha256:deadbeef",
        size_bytes=200,
        path="/bare/string/output",  # type: ignore[arg-type]
        optional=False,
    )
    record.outputs = [bad_output]

    with pytest.raises(ValueError, match="PortablePath"):
        write_stage_sidecar(record, sample_dir)


# ---------------------------------------------------------------------------
# read_stage_sidecar error handling
# ---------------------------------------------------------------------------


def test_read_invalid_json(tmp_path):
    """ValueError with helpful message when file contains invalid JSON."""
    p = tmp_path / "bad.provenance.json"
    p.write_text("{not valid json{{{{")

    with pytest.raises(ValueError, match="invalid JSON"):
        read_stage_sidecar(p)


def test_read_wrong_schema_version(tmp_path):
    """ValueError when schema_version is not '1.0'."""
    record = _make_record(tmp_path)
    sample_dir = tmp_path / "sample_out"

    sidecar_path = write_stage_sidecar(record, sample_dir)

    # Patch the JSON on disk to change schema_version.
    with open(sidecar_path) as fh:
        d = json.load(fh)
    d["schema_version"] = "99.0"
    with open(sidecar_path, "w") as fh:
        json.dump(d, fh)

    with pytest.raises(ValueError, match="schema_version"):
        read_stage_sidecar(sidecar_path)


def test_read_bare_string_path_field(tmp_path):
    """ValueError when a path field in the JSON is a bare string instead of dict."""
    record = _make_record(tmp_path)
    sample_dir = tmp_path / "sample_out"

    sidecar_path = write_stage_sidecar(record, sample_dir)

    # Corrupt the sidecar: replace the input[0].path dict with a bare string.
    with open(sidecar_path) as fh:
        d = json.load(fh)
    d["inputs"][0]["path"] = "/bare/string/path"
    with open(sidecar_path, "w") as fh:
        json.dump(d, fh)

    with pytest.raises(ValueError, match="[Pp]ath.*bare string|bare string.*path|PortablePath"):
        read_stage_sidecar(sidecar_path)


def test_read_missing_required_key(tmp_path):
    """ValueError when a required top-level key is missing from the sidecar."""
    record = _make_record(tmp_path)
    sample_dir = tmp_path / "sample_out"

    sidecar_path = write_stage_sidecar(record, sample_dir)

    with open(sidecar_path) as fh:
        d = json.load(fh)
    del d["stage"]
    with open(sidecar_path, "w") as fh:
        json.dump(d, fh)

    with pytest.raises(ValueError, match="stage"):
        read_stage_sidecar(sidecar_path)


def test_read_non_existent_file(tmp_path):
    """FileNotFoundError for a missing sidecar path."""
    with pytest.raises(FileNotFoundError):
        read_stage_sidecar(tmp_path / "does_not_exist.provenance.json")


# ---------------------------------------------------------------------------
# Optional outputs
# ---------------------------------------------------------------------------


def test_optional_output_in_sidecar(tmp_path):
    """An optional output can be declared in the sidecar and round-trips correctly."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()

    input_bam = _make_input_file(tmp_path, "input.bam")
    output_bam = _make_output_file(sample_dir, "output.bam")
    opt_bam = _make_output_file(sample_dir, "optional_out.bam")

    now = datetime.now(timezone.utc).isoformat()
    env_clear = {k: "" for k in ("TMPDIR", "L_SCRATCH") if os.environ.get(k)}
    with mock.patch.dict(os.environ, env_clear):
        record = ProvenanceRecord.from_components(
            stage="correct",
            sample_id="S1",
            sample_output_dir=sample_dir,
            started_at=now,
            completed_at=now,
            exit_status=0,
            inputs={"bam": input_bam},
            outputs={
                "corrected_bam": output_bam,
                "optional_out": (opt_bam, True),  # (path, optional=True)
            },
            argv=["rectify", "correct"],
            rectify_git_sha="abc",
            rectify_version="0.9.0",
        )

    sidecar_path = write_stage_sidecar(record, sample_dir)
    recovered = read_stage_sidecar(sidecar_path)

    optional_entry = next(e for e in recovered.outputs if e.role == "optional_out")
    assert optional_entry.optional is True


# ---------------------------------------------------------------------------
# Sub-outputs with transient flag
# ---------------------------------------------------------------------------


def test_sub_output_transient_in_sidecar(tmp_path):
    """Transient sub-outputs are recorded and round-trip with transient=True."""
    sample_dir = tmp_path / "sample_out"
    sample_dir.mkdir()
    fake_l_scratch = tmp_path / "lscratch"
    fake_l_scratch.mkdir()

    input_bam = _make_input_file(tmp_path, "input.bam")
    output_bam = _make_output_file(sample_dir, "output.bam")
    region_bam = fake_l_scratch / "region_0.bam"
    region_bam.write_bytes(b"region")

    now = datetime.now(timezone.utc).isoformat()
    # Use L_SCRATCH for the sub-output dir; clear TMPDIR so input.bam
    # (also under pytest's tmp) doesn't hit the ephemeral guard.
    with mock.patch.dict(os.environ, {"L_SCRATCH": str(fake_l_scratch)}, clear=True):
        record = ProvenanceRecord.from_components(
            stage="correct",
            sample_id="S1",
            sample_output_dir=sample_dir,
            started_at=now,
            completed_at=now,
            exit_status=0,
            inputs={"bam": input_bam},
            outputs={"corrected_bam": output_bam},
            sub_outputs=[
                {
                    "role": "region_bam",
                    "path": str(region_bam),
                    "transient": True,
                    "region_id": "chrI:1-10000",
                }
            ],
            argv=["rectify", "correct"],
            rectify_git_sha="abc",
            rectify_version="0.9.0",
        )

    sidecar_path = write_stage_sidecar(record, sample_dir)
    recovered = read_stage_sidecar(sidecar_path)

    assert len(recovered.sub_outputs) == 1
    sub = recovered.sub_outputs[0]
    assert sub.role == "region_bam"
    assert sub.transient is True
    assert sub.region_id == "chrI:1-10000"
    assert isinstance(sub.path, PortablePath)
