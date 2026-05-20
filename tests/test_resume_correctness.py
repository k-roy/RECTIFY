"""Commit B — sidecar emission ordering correctness test.

§3.6 from the briefing: sidecar MUST be written AFTER all durable outputs.
A crash BEFORE sidecar write leaves no sidecar → next run reruns (correct).
A crash AFTER outputs but BEFORE sidecar → also correct (reruns the stage).
A crash AFTER sidecar → SKIP on next run (correct — outputs are durable).

This module specifically tests the invariant that crashing mid-sidecar-write
does NOT leave a valid sidecar that causes a false SKIP on next run.

Author: Kevin R. Roy
Date: 2026-05-20
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from unittest import mock

import pysam
import pytest


# ---------------------------------------------------------------------------
# Shared fixture helpers (minimal inline version)
# ---------------------------------------------------------------------------

_CHROM_LEN = 2000


def _make_reference(tmp_path: Path) -> Path:
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
    emit_merged_tsv: bool = True,
    force_all: bool = False,
) -> argparse.Namespace:
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
        corrected_bam=None,
        softclipped_bam=None,
        bedgraph_prefix=None,
        threads=1,
        verbose=False,
        streaming=False,
        chunk_size=10000,
        checkpoint_dir=None,
        tmp_dir=None,
        emit_merged_tsv=emit_merged_tsv,
        legacy_single_threaded=False,
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
        force_all=force_all,
        force_stage=None,
        accept_prior_provenance=False,
        dry_run_resume=False,
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

def test_sidecar_not_written_on_crash_before_sidecar(tmp_path):
    """If correct_command crashes AFTER outputs but BEFORE sidecar is written,
    the sidecar must NOT exist. The next run must RUN (not SKIP).

    This is the §3.6 invariant: sidecar emission ordering.
    """
    # Patch at the import location inside correct_command module, not in sidecar module.
    import rectify.core.commands.correct_command as cc_mod

    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    sample_id = bam.stem
    sidecar_path = tmp_path / f"{sample_id}.correct.provenance.json"

    # First run: succeed normally.
    args = _make_correct_args(bam=bam, genome=genome, output=output)
    rc = _run_correct(args)
    assert rc == 0
    assert sidecar_path.exists(), "First run must write sidecar"

    # Second run: patch write_stage_sidecar at the call site inside correct_command.
    # Delete the old sidecar first to simulate a fresh-run environment.
    sidecar_path.unlink()
    assert not sidecar_path.exists()

    def _crashing_write(record, sample_output):
        raise RuntimeError("Simulated crash during sidecar write")

    # The sidecar write is inside run() via a lazy import. We need to patch
    # rectify.core.provenance.write_stage_sidecar which is called inside the
    # try/except in run(). Patch the function in the provenance package.
    with mock.patch('rectify.core.provenance.write_stage_sidecar', side_effect=_crashing_write):
        # run() catches the sidecar error (non-fatal): run completes, sidecar absent.
        # Force-all to bypass the skip check (no sidecar exists anyway).
        rc2 = _run_correct(args)

    assert rc2 == 0, "run() must not fail if sidecar write throws (it's non-fatal)"
    assert not sidecar_path.exists(), (
        "Sidecar must NOT exist after a crash during write_stage_sidecar"
    )

    # With no sidecar, can_skip_stage must return RUN.
    from rectify.core.provenance import can_skip_stage
    from rectify.utils.version import get_rectify_git_sha
    import sys

    decision = can_skip_stage(
        stage='correct',
        sample_output=tmp_path,
        sample_id=sample_id,
        current_inputs={'bam': str(bam), 'genome': str(genome)},
        current_argv=sys.argv,  # doesn't matter — no sidecar exists
        rectify_git_sha=get_rectify_git_sha(),
    )
    assert not decision.skip, (
        f"With no sidecar, can_skip_stage must return RUN, got: {decision!r}"
    )


def test_sidecar_written_on_success_enables_skip(tmp_path):
    """A successful run writes the sidecar; a second identical run returns SKIP.

    We verify the skip-behavior end-to-end by running correct_command.run() twice
    and checking the timing (SKIP should be near-instantaneous).
    """
    genome = _make_reference(tmp_path)
    bam = _make_small_bam(tmp_path)
    output = tmp_path / "corrected_reads.tsv"

    args = _make_correct_args(bam=bam, genome=genome, output=output,
                               emit_merged_tsv=True)
    rc = _run_correct(args)
    assert rc == 0

    sample_id = bam.stem
    sidecar_path = tmp_path / f"{sample_id}.correct.provenance.json"
    assert sidecar_path.exists(), "Sidecar must be written on success"

    # Second run: same args → must be fast (SKIP path).
    t0 = time.monotonic()
    rc2 = _run_correct(args)
    elapsed = time.monotonic() - t0
    assert rc2 == 0
    # SKIP should be much faster than the actual correction (~3s for 5 reads).
    assert elapsed < 3.0, (
        f"Second identical run took {elapsed:.1f}s — expected < 3s for a SKIP"
    )
