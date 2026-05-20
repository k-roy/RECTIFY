"""Tests for rectify.core.provenance.path_resolver.PortablePath.

Six scenarios per spec §6.5.7:
  (a) sample-relative tokenization with symlinks
  (b) env_relative longest-prefix-first ordering
  (c) cached_absolute fallback when env var unset on read side
  (d) cross-cluster resolution (write with $OAK set; read without)
  (e) ValueError on non-transient env_relative with ephemeral env var
  (f) argv_template tokenization round-trip preserves config_hash across
      different $L_SCRATCH values
"""

import logging
import os
from pathlib import Path
from unittest import mock

import pytest

from rectify.core.provenance.hashing import normalized_config_hash
from rectify.core.provenance.path_resolver import (
    PathResolutionError,
    PortablePath,
    tokenize_argv_paths,
)


# ---------------------------------------------------------------------------
# (a) sample-relative tokenization with symlinks
# ---------------------------------------------------------------------------


def test_sample_relative_basic(tmp_path):
    """Path under sample_output_dir tokenizes as sample_relative."""
    sample_dir = tmp_path / "sample1"
    sample_dir.mkdir()
    output_file = sample_dir / "corrected.bam"
    output_file.write_text("data")

    pp = PortablePath.from_path(output_file, sample_output_dir=sample_dir)
    assert pp.kind == "sample_relative"
    assert pp.template == "corrected.bam"
    assert pp.env_var is None


def test_sample_relative_with_symlink(tmp_path):
    """Sample-relative tokenization works even when the sample dir is a symlink target."""
    real_dir = tmp_path / "real_sample1"
    real_dir.mkdir()
    symlink_dir = tmp_path / "symlink_sample1"
    symlink_dir.symlink_to(real_dir)

    output_file = real_dir / "output.bam"
    output_file.write_text("bam")

    # Tokenize via the symlink path — should still resolve to sample_relative
    # because from_path() calls os.path.realpath() on both sides.
    pp = PortablePath.from_path(output_file, sample_output_dir=symlink_dir)
    assert pp.kind == "sample_relative"
    assert pp.template == "output.bam"


def test_sample_relative_resolve(tmp_path):
    """sample_relative resolve_now returns the correct path."""
    sample_dir = tmp_path / "s1"
    sample_dir.mkdir()
    f = sample_dir / "foo.tsv"
    f.write_text("x")

    pp = PortablePath.from_path(f, sample_output_dir=sample_dir)
    resolved = pp.resolve_now(sample_output_dir=sample_dir)
    assert resolved.exists()
    assert resolved.name == "foo.tsv"


def test_sample_relative_requires_sample_output_dir_on_resolve(tmp_path):
    """sample_relative requires sample_output_dir on resolve; missing it raises ValueError."""
    sample_dir = tmp_path / "s1"
    sample_dir.mkdir()
    f = sample_dir / "foo.tsv"
    f.write_text("x")

    pp = PortablePath.from_path(f, sample_output_dir=sample_dir)
    with pytest.raises(ValueError, match="sample_output_dir"):
        pp.resolve_now(sample_output_dir=None)


# ---------------------------------------------------------------------------
# (b) env_relative longest-prefix-first ordering
# ---------------------------------------------------------------------------


def test_env_relative_longest_prefix_wins(tmp_path):
    """When $L_SCRATCH is a sub-path of $HOME, L_SCRATCH wins (longest prefix)."""
    home_dir = tmp_path / "fake_home"
    scratch_dir = home_dir / "lscratch"  # L_SCRATCH is under HOME
    scratch_dir.mkdir(parents=True)

    target_file = scratch_dir / "region_bam.bam"
    target_file.write_text("bam")

    env = {
        "HOME": str(home_dir),
        "L_SCRATCH": str(scratch_dir),
    }
    with mock.patch.dict(os.environ, env, clear=True):
        pp = PortablePath.from_path(target_file, transient=True)

    # L_SCRATCH should win because it's the longest matching prefix.
    assert pp.kind == "env_relative"
    assert pp.env_var == "L_SCRATCH"
    assert pp.template.startswith("$L_SCRATCH")


def test_env_relative_scratch_fallback(tmp_path):
    """When only $SCRATCH covers the path, uses $SCRATCH."""
    scratch_dir = tmp_path / "fake_scratch"
    scratch_dir.mkdir()
    f = scratch_dir / "data.bam"
    f.write_text("data")

    env = {"SCRATCH": str(scratch_dir)}
    with mock.patch.dict(os.environ, env, clear=True):
        pp = PortablePath.from_path(f)

    assert pp.kind == "env_relative"
    assert pp.env_var == "SCRATCH"


def test_env_relative_absolute_fallback_when_no_env_matches(tmp_path):
    """When no env var covers the path, kind='absolute'."""
    f = tmp_path / "no_env_match.bam"
    f.write_text("x")

    with mock.patch.dict(os.environ, {}, clear=True):
        pp = PortablePath.from_path(f)

    assert pp.kind == "absolute"


# ---------------------------------------------------------------------------
# (c) cached_absolute fallback when env var unset on read side
# ---------------------------------------------------------------------------


def test_cached_absolute_fallback_on_unset_env(tmp_path, caplog):
    """When env var is unset on resolve, falls back to cached_absolute."""
    fake_scratch = tmp_path / "scratch"
    fake_scratch.mkdir()
    f = fake_scratch / "temp_bam.bam"
    f.write_text("bam")

    # Write: env var set → env_relative
    with mock.patch.dict(os.environ, {"SCRATCH": str(fake_scratch)}, clear=True):
        pp = PortablePath.from_path(f)

    assert pp.kind == "env_relative"
    assert pp.env_var == "SCRATCH"

    # Resolve: env var unset → must fall back to cached_absolute
    with mock.patch.dict(os.environ, {}, clear=True):
        with caplog.at_level(logging.WARNING):
            resolved = pp.resolve_now()

    assert resolved == f.resolve() or resolved == Path(pp.cached_absolute)
    assert resolved.exists()
    # Should have logged a warning about the fallback
    assert any("fallback" in rec.message.lower() or "unset" in rec.message.lower()
               for rec in caplog.records), (
        "Expected a WARNING about env var unset / fallback to cached_absolute"
    )


# ---------------------------------------------------------------------------
# (d) cross-cluster resolution ($OAK write, no $OAK on read)
# ---------------------------------------------------------------------------


def test_cross_cluster_fallback_to_cached_absolute(tmp_path, caplog):
    """Write sidecar on cluster with $OAK set; read on machine where $OAK is unset."""
    fake_oak = tmp_path / "oak"
    fake_oak.mkdir()
    f = fake_oak / "genome.fa"
    f.write_text("ACGT")

    # Simulate writing on Sherlock (OAK set)
    with mock.patch.dict(os.environ, {"OAK": str(fake_oak)}, clear=True):
        pp = PortablePath.from_path(f)

    assert pp.kind == "env_relative"
    assert pp.env_var == "OAK"

    # Simulate resolving on a different cluster (OAK unset)
    # File still exists at cached_absolute.
    with mock.patch.dict(os.environ, {}, clear=True):
        with caplog.at_level(logging.WARNING):
            resolved = pp.resolve_now()

    assert resolved.exists()
    # Verify fallback warning was emitted
    assert any("fallback" in r.message.lower() or "unset" in r.message.lower()
               for r in caplog.records)


def test_cross_cluster_raises_when_neither_exists(tmp_path):
    """If neither template-expanded nor cached_absolute exists, raises PathResolutionError."""
    fake_oak = tmp_path / "oak"
    fake_oak.mkdir()
    f = fake_oak / "transient.bam"
    f.write_text("data")

    with mock.patch.dict(os.environ, {"OAK": str(fake_oak)}, clear=True):
        pp = PortablePath.from_path(f)

    # Delete the file — both candidates are now gone
    f.unlink()

    with mock.patch.dict(os.environ, {}, clear=True):
        with pytest.raises(PathResolutionError):
            pp.resolve_now()


# ---------------------------------------------------------------------------
# (e) ValueError on non-transient env_relative with ephemeral env var
# ---------------------------------------------------------------------------


def test_reject_non_transient_l_scratch_path(tmp_path):
    """Non-transient path under $L_SCRATCH → ValueError."""
    fake_l_scratch = tmp_path / "lscratch_dir"
    fake_l_scratch.mkdir()
    f = fake_l_scratch / "output.bam"
    f.write_text("data")

    with mock.patch.dict(os.environ, {"L_SCRATCH": str(fake_l_scratch)}, clear=True):
        with pytest.raises(ValueError, match="L_SCRATCH"):
            PortablePath.from_path(f, transient=False)


def test_reject_non_transient_tmpdir_path(tmp_path):
    """Non-transient path under $TMPDIR → ValueError."""
    fake_tmpdir = tmp_path / "tmpdir"
    fake_tmpdir.mkdir()
    f = fake_tmpdir / "output.bam"
    f.write_text("data")

    with mock.patch.dict(os.environ, {"TMPDIR": str(fake_tmpdir)}, clear=True):
        with pytest.raises(ValueError, match="TMPDIR"):
            PortablePath.from_path(f, transient=False)


def test_allow_transient_l_scratch_path(tmp_path):
    """Transient path under $L_SCRATCH → allowed (kind=env_relative)."""
    fake_l_scratch = tmp_path / "lscratch_dir"
    fake_l_scratch.mkdir()
    f = fake_l_scratch / "region.bam"
    f.write_text("data")

    with mock.patch.dict(os.environ, {"L_SCRATCH": str(fake_l_scratch)}, clear=True):
        pp = PortablePath.from_path(f, transient=True)

    assert pp.kind == "env_relative"
    assert pp.env_var == "L_SCRATCH"


# ---------------------------------------------------------------------------
# (f) argv_template tokenization round-trip preserves config_hash
# ---------------------------------------------------------------------------


def test_argv_template_stable_across_l_scratch_values(tmp_path):
    """config_hash is the same across two runs with different $L_SCRATCH values.

    The argv mentions a path in $L_SCRATCH. tokenize_argv_paths replaces it
    with the $L_SCRATCH template. normalized_config_hash then ignores the
    specific numeric value, so config_hash is stable.
    """
    sample_dir = tmp_path / "sample1"
    sample_dir.mkdir()

    # Run 1: L_SCRATCH = /lscratch/111
    lscratch1 = tmp_path / "lscratch_111"
    lscratch1.mkdir()
    bam1 = lscratch1 / "region_0.bam"
    bam1.write_text("x")

    argv1 = ["rectify", "correct", "--bam", str(bam1), "--sample", "S1"]
    with mock.patch.dict(os.environ, {"L_SCRATCH": str(lscratch1)}, clear=True):
        template1 = tokenize_argv_paths(argv1, sample_dir)
    hash1 = normalized_config_hash(template1)

    # Run 2: L_SCRATCH = /lscratch/222 (different job id)
    lscratch2 = tmp_path / "lscratch_222"
    lscratch2.mkdir()
    bam2 = lscratch2 / "region_0.bam"
    bam2.write_text("x")

    argv2 = ["rectify", "correct", "--bam", str(bam2), "--sample", "S1"]
    with mock.patch.dict(os.environ, {"L_SCRATCH": str(lscratch2)}, clear=True):
        template2 = tokenize_argv_paths(argv2, sample_dir)
    hash2 = normalized_config_hash(template2)

    # Both tokenizations should produce $L_SCRATCH/region_0.bam
    assert any("$L_SCRATCH" in tok for tok in template1), (
        f"Expected $L_SCRATCH template in {template1}"
    )
    assert template1 == template2, (
        f"argv_templates should be equal across L_SCRATCH runs:\n"
        f"  run1: {template1}\n  run2: {template2}"
    )
    assert hash1 == hash2, "config_hash should be stable across different $L_SCRATCH job ids"


# ---------------------------------------------------------------------------
# JSON round-trip
# ---------------------------------------------------------------------------


def test_portable_path_json_roundtrip_sample_relative(tmp_path):
    """PortablePath.to_json() + from_json() round-trip for sample_relative."""
    sample_dir = tmp_path / "s1"
    sample_dir.mkdir()
    f = sample_dir / "out.bam"
    f.write_text("data")

    pp = PortablePath.from_path(f, sample_output_dir=sample_dir)
    d = pp.to_json()
    pp2 = PortablePath.from_json(d)
    assert pp == pp2


def test_portable_path_json_roundtrip_env_relative(tmp_path):
    """PortablePath.to_json() + from_json() round-trip for env_relative."""
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    f = scratch / "file.bam"
    f.write_text("x")

    with mock.patch.dict(os.environ, {"SCRATCH": str(scratch)}, clear=True):
        pp = PortablePath.from_path(f)

    d = pp.to_json()
    pp2 = PortablePath.from_json(d)
    assert pp == pp2


def test_portable_path_from_json_rejects_non_dict():
    """from_json raises ValueError if given a bare string instead of a dict."""
    with pytest.raises(ValueError, match="expected dict"):
        PortablePath.from_json("/bare/string/path")


def test_portable_path_from_json_rejects_missing_keys():
    """from_json raises ValueError if required keys are missing."""
    with pytest.raises(ValueError, match="missing required keys"):
        PortablePath.from_json({"template": "foo"})
