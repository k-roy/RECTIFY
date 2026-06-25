"""Unit tests for rectify.utils.bam_provenance."""

import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from rectify.utils.bam_provenance import (
    SIDECAR_SCHEMA_VERSION,
    SIDECAR_SUFFIX,
    _ALIGNER_VERSION_CACHE,
    _GIT_SHA_CACHE,
    compute_run_provenance,
    expected_provenance_for_aligner,
    get_aligner_version,
    get_rectify_git_sha,
    matches_strict,
    read_sidecar,
    sidecar_path,
    write_sidecar,
)
import rectify.utils.bam_provenance as _bp_mod


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clear_caches():
    """Reset module-level caches before each test that exercises them."""
    _ALIGNER_VERSION_CACHE.clear()
    # _GIT_SHA_CACHE is a module global; reset via the module reference.
    _bp_mod._GIT_SHA_CACHE = None


# ---------------------------------------------------------------------------
# sidecar_path
# ---------------------------------------------------------------------------

def test_sidecar_path_str(tmp_path):
    bam = tmp_path / "sample.minimap2.bam"
    assert sidecar_path(bam) == Path(str(bam) + SIDECAR_SUFFIX)


def test_sidecar_path_accepts_string(tmp_path):
    bam = str(tmp_path / "sample.minimap2.bam")
    result = sidecar_path(bam)
    assert str(result).endswith(SIDECAR_SUFFIX)


# ---------------------------------------------------------------------------
# write_sidecar / read_sidecar round-trip
# ---------------------------------------------------------------------------

def test_round_trip_per_aligner(tmp_path):
    bam = tmp_path / "sample.minimap2.bam"
    bam.write_bytes(b"")
    run_prov = compute_run_provenance(command=["rectify", "run-all"])
    sc_path = write_sidecar(bam, run_prov, aligner_name="minimap2", aligner_version="2.26-r1175")
    assert sc_path.exists()
    stored = read_sidecar(bam)
    assert stored is not None
    assert stored["aligner_name"] == "minimap2"
    assert stored["aligner_version"] == "2.26-r1175"
    assert stored["schema_version"] == SIDECAR_SCHEMA_VERSION


def test_round_trip_consensus(tmp_path):
    """consensus aligner_version should round-trip via the special-case probe."""
    _clear_caches()
    bam = tmp_path / "sample.multialigned.bam"
    bam.write_bytes(b"")
    run_prov = compute_run_provenance(command=["rectify", "run-all"])
    # No explicit aligner_version — let auto-probe handle it.
    write_sidecar(bam, run_prov, aligner_name="consensus")
    stored = read_sidecar(bam)
    assert stored is not None
    assert stored["aligner_name"] == "consensus"
    # Version should equal the package version (not None).
    from rectify import __version__
    assert stored["aligner_version"] == __version__


def test_round_trip_consensus_matches_strict(tmp_path):
    """matches_strict must return True when comparing stored vs expected consensus."""
    _clear_caches()
    bam = tmp_path / "sample.multialigned.bam"
    bam.write_bytes(b"")
    run_prov = compute_run_provenance(command=sys.argv)
    write_sidecar(bam, run_prov, aligner_name="consensus")
    stored = read_sidecar(bam)
    expected = expected_provenance_for_aligner(run_prov, "consensus")
    ok, reason = matches_strict(stored, expected)
    assert ok, f"Expected match but got: {reason}"


def test_read_sidecar_missing(tmp_path):
    bam = tmp_path / "nonexistent.bam"
    assert read_sidecar(bam) is None


def test_read_sidecar_corrupt(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    sidecar_path(bam).write_text("NOT JSON {{{{")
    assert read_sidecar(bam) is None


def test_write_sidecar_atomic(tmp_path):
    """Sidecar must not leave a .tmp file behind on success."""
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"")
    run_prov = compute_run_provenance()
    write_sidecar(bam, run_prov, aligner_name="minimap2", aligner_version="2.26")
    tmp_files = list(tmp_path.glob("*.tmp"))
    assert tmp_files == [], f"Leftover .tmp files: {tmp_files}"


# ---------------------------------------------------------------------------
# matches_strict
# ---------------------------------------------------------------------------

def test_matches_strict_none_stored():
    current = {"rectify_git_sha": "abc", "aligner_name": "minimap2", "aligner_version": "2.26"}
    ok, reason = matches_strict(None, current)
    assert not ok
    assert "no sidecar" in reason


def test_matches_strict_all_match():
    d = {"rectify_git_sha": "abc123", "aligner_name": "minimap2", "aligner_version": "2.26"}
    ok, reason = matches_strict(d, d)
    assert ok
    assert reason == "match"


def test_matches_strict_sha_mismatch():
    stored = {"rectify_git_sha": "old_sha", "aligner_name": "minimap2", "aligner_version": "2.26"}
    current = {"rectify_git_sha": "new_sha", "aligner_name": "minimap2", "aligner_version": "2.26"}
    ok, reason = matches_strict(stored, current)
    assert not ok
    assert "rectify_git_sha" in reason


def test_matches_strict_aligner_name_mismatch():
    stored = {"rectify_git_sha": "abc", "aligner_name": "minimap2", "aligner_version": "2.26"}
    current = {"rectify_git_sha": "abc", "aligner_name": "gapmm2", "aligner_version": "2.26"}
    ok, reason = matches_strict(stored, current)
    assert not ok
    assert "aligner_name" in reason


def test_matches_strict_aligner_version_mismatch():
    stored = {"rectify_git_sha": "abc", "aligner_name": "minimap2", "aligner_version": "2.26"}
    current = {"rectify_git_sha": "abc", "aligner_name": "minimap2", "aligner_version": "2.27"}
    ok, reason = matches_strict(stored, current)
    assert not ok
    assert "aligner_version" in reason


def test_matches_strict_version_none_vs_none():
    """Both stored and current have version=None (probe unavailable) → still match."""
    stored = {"rectify_git_sha": "abc", "aligner_name": "uLTRA", "aligner_version": None}
    current = {"rectify_git_sha": "abc", "aligner_name": "uLTRA", "aligner_version": None}
    ok, _ = matches_strict(stored, current)
    assert ok


# ---------------------------------------------------------------------------
# get_rectify_git_sha
# ---------------------------------------------------------------------------

def test_get_rectify_git_sha_does_not_crash():
    """Should return a string or None — never raise."""
    _clear_caches()
    result = get_rectify_git_sha()
    assert result is None or isinstance(result, str)


def test_get_rectify_git_sha_cached():
    """Second call must return same value without hitting subprocess."""
    _clear_caches()
    first = get_rectify_git_sha()
    with patch("subprocess.run", side_effect=AssertionError("should be cached")):
        second = get_rectify_git_sha()
    assert first == second


# ---------------------------------------------------------------------------
# get_aligner_version
# ---------------------------------------------------------------------------

def test_get_aligner_version_consensus():
    """'consensus' special-case must return the package version, not None."""
    _clear_caches()
    from rectify import __version__
    assert get_aligner_version("consensus") == __version__


def test_get_aligner_version_unknown():
    _clear_caches()
    assert get_aligner_version("unknown_aligner_xyz") is None


def test_get_aligner_version_not_on_path():
    _clear_caches()
    with patch("shutil.which", return_value=None):
        result = get_aligner_version("minimap2")
    assert result is None


def test_get_aligner_version_minimap2_mock():
    """Mock subprocess.run to verify regex parsing."""
    _clear_caches()
    fake_result = type("R", (), {"returncode": 0, "stdout": "2.26-r1175\n", "stderr": ""})()
    with patch("shutil.which", return_value="/usr/bin/minimap2"), \
         patch("subprocess.run", return_value=fake_result):
        ver = get_aligner_version("minimap2")
    assert ver == "2.26-r1175"


# ---------------------------------------------------------------------------
# compute_run_provenance
# ---------------------------------------------------------------------------

def test_compute_run_provenance_fields():
    prov = compute_run_provenance(command=["rectify", "run-all"], extra={"custom": 1})
    assert prov["schema_version"] == SIDECAR_SCHEMA_VERSION
    assert "rectify_pkg_version" in prov
    assert "rectify_git_sha" in prov
    assert "created_at" in prov
    assert prov["command"] == ["rectify", "run-all"]
    assert prov["custom"] == 1
