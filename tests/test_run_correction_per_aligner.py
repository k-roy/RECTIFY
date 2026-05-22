"""Tests for ``_run_correction_per_aligner`` (Preparatory Commit C).

Two behaviors validated here:

1. **Manifest-output acceptance** — Commit B made
   ``corrected_reads.manifest.tsv`` the canonical artifact (with the merged
   ``corrected_reads.tsv`` renamed to ``corrected_reads.region_000.tsv``).
   The per-aligner driver must accept either form for both the resume
   skip-check and the post-call success check.

2. **Aligner-concurrency plumbing** — ``args.aligner_concurrency`` is
   resolved via :func:`resolve_aligner_concurrency` and observed. The loop
   itself stays sequential (Axis B performance is deferred), but the
   resolved value is reported and an invalid value falls back to ``1``
   instead of crashing.

These tests monkeypatch ``_run_correction`` so no real aligner BAMs are
needed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional

import pytest


MANIFEST_HEADER = "region_id\tchrom\tstart\tend\ttsv_path\tn_rows\tsha256\n"
SINGLE_TSV_HEADER = "read_id\tchrom\tstrand\tcorrected_3prime\n"


def _make_args(threads: int = 4, aligner_concurrency: str = "auto") -> argparse.Namespace:
    return argparse.Namespace(
        threads=threads,
        aligner_concurrency=aligner_concurrency,
        dT_primed_cDNA=False,
        polya_sequenced=False,
        organism=None,
        aligner='minimap2',
        chunk_size=10000,
        netseq_dir=None,
        junction_penalty_table=None,
        str_penalty_table=None,
        junction_overhang_table=None,
        junction_pool_cache=None,
        filter_spikein=None,
        streaming=False,
        write_corrected_bam=False,
        write_softclip_bam=False,
    )


def _make_per_aligner_bams(tmp_path: Path, names) -> Dict[str, Path]:
    """Create dummy .bam + .bai files so the indexed-bam guard doesn't sort."""
    out: Dict[str, Path] = {}
    for name in names:
        bam = tmp_path / f"{name}.rectified.bam"
        bai = tmp_path / f"{name}.rectified.bam.bai"
        bam.touch()
        bai.touch()
        out[name] = bam
    return out


def _write_manifest(aligner_output_dir: Path) -> Path:
    """Materialize the Commit-B default layout for one aligner."""
    aligner_output_dir.mkdir(parents=True, exist_ok=True)
    region_tsv = aligner_output_dir / 'corrected_reads.region_000.tsv'
    region_tsv.write_text(SINGLE_TSV_HEADER)
    manifest = aligner_output_dir / 'corrected_reads.manifest.tsv'
    manifest.write_text(
        MANIFEST_HEADER
        + f"region_000\t*\t0\t0\t{region_tsv.name}\t0\tdeadbeef\n"
    )
    return manifest


def _write_legacy_tsv(aligner_output_dir: Path) -> Path:
    """Materialize the ``--emit-merged-tsv`` legacy layout for one aligner."""
    aligner_output_dir.mkdir(parents=True, exist_ok=True)
    tsv = aligner_output_dir / 'corrected_reads.tsv'
    tsv.write_text(SINGLE_TSV_HEADER)
    return tsv


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_per_aligner_canonical_prefers_manifest(tmp_path: Path):
    """The canonical-output helper returns manifest when both forms exist."""
    from rectify.core.commands.run.stages import _per_aligner_canonical_output

    _write_manifest(tmp_path)
    legacy = _write_legacy_tsv(tmp_path)
    canonical = _per_aligner_canonical_output(tmp_path)

    assert canonical is not None
    assert canonical.name == 'corrected_reads.manifest.tsv'
    assert legacy.exists()  # legacy file still present, just not chosen


def test_per_aligner_canonical_falls_back_to_legacy_tsv(tmp_path: Path):
    from rectify.core.commands.run.stages import _per_aligner_canonical_output

    _write_legacy_tsv(tmp_path)
    canonical = _per_aligner_canonical_output(tmp_path)

    assert canonical is not None
    assert canonical.name == 'corrected_reads.tsv'


def test_per_aligner_canonical_none_when_missing(tmp_path: Path):
    from rectify.core.commands.run.stages import _per_aligner_canonical_output

    assert _per_aligner_canonical_output(tmp_path) is None


def test_per_aligner_runner_skips_when_manifest_already_exists(
    tmp_path: Path, monkeypatch
):
    """Resume path: pre-existing manifest output should be detected and skipped.

    Before the Commit B regression fix the runner only checked
    ``corrected_reads.tsv`` and would have re-run correction on every resume.
    """
    from rectify.core.commands.run import stages

    per_aligner_bams = _make_per_aligner_bams(tmp_path, ['minimap2'])
    output_dir = tmp_path / 'out'
    # Pre-seed the manifest as if a prior run completed.
    _write_manifest(output_dir / 'per_aligner_corrected' / 'minimap2')

    called = {'count': 0}

    def _fake_run_correction(**kwargs):  # pragma: no cover - defensive
        called['count'] += 1
        return kwargs['output_dir'] / 'corrected_reads.tsv'

    monkeypatch.setattr(stages, '_run_correction', _fake_run_correction)

    tsvs, _ = stages._run_correction_per_aligner(
        per_aligner_bams=per_aligner_bams,
        output_dir=output_dir,
        genome_path=tmp_path / 'genome.fa',
        annotation_path=None,
        args=_make_args(aligner_concurrency='1'),
    )

    assert called['count'] == 0, "should skip when manifest already exists"
    assert tsvs['minimap2'].name == 'corrected_reads.manifest.tsv'


def test_per_aligner_runner_accepts_manifest_after_correction(
    tmp_path: Path, monkeypatch
):
    """Post-call success path: writer that emits manifest-only is recognized.

    Mimics what ``correct_command.run`` does under the Commit B default
    (rename merged TSV to region_000.tsv; write manifest.tsv).
    """
    from rectify.core.commands.run import stages

    per_aligner_bams = _make_per_aligner_bams(tmp_path, ['minimap2', 'gapmm2'])
    output_dir = tmp_path / 'out'

    def _fake_run_correction(*, output_dir, **_kwargs):
        _write_manifest(output_dir)
        return output_dir / 'corrected_reads.manifest.tsv'

    monkeypatch.setattr(stages, '_run_correction', _fake_run_correction)

    tsvs, _ = stages._run_correction_per_aligner(
        per_aligner_bams=per_aligner_bams,
        output_dir=output_dir,
        genome_path=tmp_path / 'genome.fa',
        annotation_path=None,
        args=_make_args(aligner_concurrency='1'),
    )

    assert set(tsvs) == {'minimap2', 'gapmm2'}
    for path in tsvs.values():
        assert path.name == 'corrected_reads.manifest.tsv'
        assert path.exists()


def test_per_aligner_runner_accepts_legacy_tsv(tmp_path: Path, monkeypatch):
    """Back-compat: ``--emit-merged-tsv`` writers still pass success check."""
    from rectify.core.commands.run import stages

    per_aligner_bams = _make_per_aligner_bams(tmp_path, ['minimap2'])
    output_dir = tmp_path / 'out'

    def _fake_run_correction(*, output_dir, **_kwargs):
        _write_legacy_tsv(output_dir)
        return output_dir / 'corrected_reads.tsv'

    monkeypatch.setattr(stages, '_run_correction', _fake_run_correction)

    tsvs, _ = stages._run_correction_per_aligner(
        per_aligner_bams=per_aligner_bams,
        output_dir=output_dir,
        genome_path=tmp_path / 'genome.fa',
        annotation_path=None,
        args=_make_args(aligner_concurrency='1'),
    )

    assert tsvs['minimap2'].name == 'corrected_reads.tsv'


def test_per_aligner_runner_serial_loop_under_concurrency_auto(
    tmp_path: Path, monkeypatch, capsys
):
    """``--aligner-concurrency auto`` keeps the loop sequential and logs budget.

    Sequential is the Preparatory Commit C invariant: the resolved value is
    informational only. The acceptance test for Axis B performance is
    explicitly deferred until the shared region-task API exists below
    ``correct_command.run``.
    """
    from rectify.core.commands.run import stages

    per_aligner_bams = _make_per_aligner_bams(
        tmp_path, ['minimap2', 'gapmm2', 'uLTRA']
    )
    output_dir = tmp_path / 'out'

    call_order = []

    def _fake_run_correction(*, output_dir, **_kwargs):
        call_order.append(output_dir.name)
        _write_manifest(output_dir)
        return output_dir / 'corrected_reads.manifest.tsv'

    monkeypatch.setattr(stages, '_run_correction', _fake_run_correction)

    stages._run_correction_per_aligner(
        per_aligner_bams=per_aligner_bams,
        output_dir=output_dir,
        genome_path=tmp_path / 'genome.fa',
        annotation_path=None,
        args=_make_args(threads=16, aligner_concurrency='auto'),
    )

    # Sequential execution preserves dict insertion order.
    assert call_order == ['minimap2', 'gapmm2', 'uLTRA']

    out = capsys.readouterr().out
    assert '--aligner-concurrency=auto' in out
    assert 'loop stays sequential' in out


def test_per_aligner_runner_invalid_concurrency_falls_back_to_one(
    tmp_path: Path, monkeypatch, capsys
):
    """An unparseable value warns and falls back to 1; loop still runs."""
    from rectify.core.commands.run import stages

    per_aligner_bams = _make_per_aligner_bams(tmp_path, ['minimap2'])
    output_dir = tmp_path / 'out'

    def _fake_run_correction(*, output_dir, **_kwargs):
        _write_manifest(output_dir)
        return output_dir / 'corrected_reads.manifest.tsv'

    monkeypatch.setattr(stages, '_run_correction', _fake_run_correction)

    tsvs, _ = stages._run_correction_per_aligner(
        per_aligner_bams=per_aligner_bams,
        output_dir=output_dir,
        genome_path=tmp_path / 'genome.fa',
        annotation_path=None,
        args=_make_args(aligner_concurrency='not-a-number'),
    )

    err = capsys.readouterr().err
    assert 'invalid' in err.lower() or 'WARNING' in err
    assert tsvs['minimap2'].exists()
