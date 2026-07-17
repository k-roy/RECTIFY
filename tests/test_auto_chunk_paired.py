"""Paired (--read2) split must honor --auto-chunk.

Regression: ``run_split`` dispatched to ``_run_split_paired`` and returned BEFORE
reaching the single-end auto-chunk calibration block, and ``_run_split_paired``
never referenced ``auto_chunk`` — so ``rectify split --read2 R2 R1 --auto-chunk``
silently used ``DEFAULT_TARGET_READS`` (50k) instead of a calibrated target.

The fix resolves auto-chunk in ``run_split`` (shared ``_resolve_auto_chunk``)
BEFORE the paired dispatch, so ``_run_split_paired`` reads the calibrated target.
These tests stub the probe/dispatch so nothing shells out to real aligners.
"""
import argparse
from pathlib import Path

from rectify.core.commands import split_command
from rectify.core.commands.split_command import (
    DEFAULT_TARGET_READS,
    _resolve_auto_chunk,
    run_split,
)


def _paired_args(reads, read2, auto_chunk):
    return argparse.Namespace(
        reads=reads,
        read2=read2,
        auto_chunk=auto_chunk,
        verbose=False,
        target_reads_per_chunk=DEFAULT_TARGET_READS,
        n_chunks=None,
    )


# ── unit tests of the shared resolver ────────────────────────────────────────

def test_resolve_auto_chunk_sets_calibrated_target(monkeypatch):
    monkeypatch.setattr(split_command, 'count_reads', lambda p: 1_000_000)
    monkeypatch.setattr(split_command, '_run_calibration_probe', lambda a, f, t: 12345)
    args = argparse.Namespace(auto_chunk=True, reads=Path('R1.fq.gz'),
                              target_reads_per_chunk=DEFAULT_TARGET_READS, n_chunks=7)
    _resolve_auto_chunk(args)
    assert args.target_reads_per_chunk == 12345      # calibrated
    assert args.n_chunks is None                     # forces target-reads path


def test_resolve_auto_chunk_noop_when_disabled(monkeypatch):
    called = []
    monkeypatch.setattr(split_command, '_run_calibration_probe',
                        lambda *a: called.append(1) or 999)
    args = argparse.Namespace(auto_chunk=False, reads=Path('R1.fq.gz'),
                              target_reads_per_chunk=DEFAULT_TARGET_READS, n_chunks=7)
    _resolve_auto_chunk(args)
    assert not called                                # probe never runs
    assert args.target_reads_per_chunk == DEFAULT_TARGET_READS
    assert args.n_chunks == 7                         # untouched


def test_resolve_auto_chunk_calibration_failure_keeps_default(monkeypatch):
    monkeypatch.setattr(split_command, 'count_reads', lambda p: 1_000_000)
    monkeypatch.setattr(split_command, '_run_calibration_probe', lambda a, f, t: None)
    args = argparse.Namespace(auto_chunk=True, reads=Path('R1.fq.gz'),
                              target_reads_per_chunk=DEFAULT_TARGET_READS, n_chunks=3)
    _resolve_auto_chunk(args)
    assert args.target_reads_per_chunk == DEFAULT_TARGET_READS  # unchanged
    assert args.n_chunks == 3                                   # unchanged on failure


# ── the regression: the paired path must be wired to the resolver ────────────

def test_paired_split_honors_auto_chunk_before_dispatch(monkeypatch, tmp_path):
    r1 = tmp_path / 'R1.fastq.gz'; r1.write_text('')
    r2 = tmp_path / 'R2.fastq.gz'; r2.write_text('')
    monkeypatch.setattr(split_command, 'count_reads', lambda p: 1_000_000)
    monkeypatch.setattr(split_command, '_run_calibration_probe', lambda a, f, t: 12345)
    captured = {}

    def _fake_paired(a):
        captured['target'] = a.target_reads_per_chunk
        captured['n_chunks'] = a.n_chunks
        return 0

    monkeypatch.setattr(split_command, '_run_split_paired', _fake_paired)
    rc = run_split(_paired_args(r1, r2, auto_chunk=True))
    assert rc == 0
    # The calibrated target reached _run_split_paired — NOT the silent 50k default.
    assert captured['target'] == 12345
    assert captured['target'] != DEFAULT_TARGET_READS
    assert captured['n_chunks'] is None


def test_paired_split_without_auto_chunk_keeps_default(monkeypatch, tmp_path):
    r1 = tmp_path / 'R1.fastq.gz'; r1.write_text('')
    r2 = tmp_path / 'R2.fastq.gz'; r2.write_text('')
    probe_called = []
    monkeypatch.setattr(split_command, '_run_calibration_probe',
                        lambda *a: probe_called.append(1) or 999)
    captured = {}
    monkeypatch.setattr(split_command, '_run_split_paired',
                        lambda a: captured.update(target=a.target_reads_per_chunk) or 0)
    rc = run_split(_paired_args(r1, r2, auto_chunk=False))
    assert rc == 0
    assert not probe_called                          # no calibration without --auto-chunk
    assert captured['target'] == DEFAULT_TARGET_READS
