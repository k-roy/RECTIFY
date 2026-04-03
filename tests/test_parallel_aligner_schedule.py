#!/usr/bin/env python3
"""
Tests for the two-phase parallel aligner scheduling in align_command.run_align().

Phase 1: mapPacBio runs alone with all available threads.
Phase 2: remaining aligners (minimap2, gapmm2, uLTRA) run in parallel with
         equal thread shares; deSALT runs sequentially after the pool.

Author: Kevin R. Roy
Date: 2026-04-01
"""

import argparse
import types
from pathlib import Path
from unittest.mock import MagicMock, call, patch

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_args(
    threads=16,
    parallel_aligners=True,
    aligners=None,
    output_dir=None,
    tmp_path=None,
):
    """Build a minimal args namespace for run_align()."""
    base = tmp_path or Path('/tmp')
    # run_align validates that reads and genome files exist on disk
    reads = base / 'reads.fastq'
    genome = base / 'genome.fa'
    anno = base / 'anno.gff'
    reads.touch()
    genome.touch()
    anno.touch()

    ns = argparse.Namespace()
    ns.threads = threads
    ns.parallel_aligners = parallel_aligners
    ns.reads = reads
    ns.genome = genome
    ns.annotation = anno
    ns.output_dir = output_dir or base
    ns.prefix = None
    ns.minimap2_path = 'minimap2'
    ns.mapPacBio_path = 'mapPacBio.sh'
    ns.gapmm2_path = 'gapmm2'
    ns.ultra_path = 'uLTRA'
    ns.desalt_path = 'deSALT'
    ns.junc_bonus = 9
    ns.no_consensus = True   # skip consensus — we only test the aligner schedule
    ns.junction_aligners = ['uLTRA', 'deSALT']
    ns.verbose = False
    ns.junc_bed = None
    ns.aligners = aligners or ['minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT']
    return ns


def _aligner_names_from_calls(mock_fn):
    """Extract the 'aligner' keyword or first positional arg from each call."""
    names = []
    for c in mock_fn.call_args_list:
        # Each call is run_minimap2(...) / run_map_pacbio(...) etc. — capture
        # which function was called by inspecting the mock name.
        names.append(mock_fn._mock_name or str(c))
    return names


# ---------------------------------------------------------------------------
# Core scheduling tests
# ---------------------------------------------------------------------------

class TestTwoPhaseAlignerSchedule:
    """mapPacBio runs first (all threads), then others with equal shares."""

    @pytest.fixture()
    def _patch_aligners(self):
        """Patch all five aligner runners to record thread counts."""
        thread_log = {}

        def _make_runner(name):
            def _runner(*args, threads=1, **kwargs):
                thread_log.setdefault(name, []).append(threads)
                return None   # No real BAM produced — that's fine
            return _runner

        # Aligner functions are imported locally inside run_align() via
        # "from .multi_aligner import ..."; patch them on the source module.
        patches = [
            patch(
                f'rectify.core.multi_aligner.{fn}',
                side_effect=_make_runner(name),
            )
            for fn, name in [
                ('run_minimap2', 'minimap2'),
                ('run_map_pacbio', 'mapPacBio'),
                ('run_gapmm2', 'gapmm2'),
                ('run_ultra', 'uLTRA'),
                ('run_desalt', 'deSALT'),
            ]
        ]
        # Also patch check_aligner_available → always True
        patches.append(
            patch('rectify.core.multi_aligner.check_aligner_available', return_value=True)
        )
        for p in patches:
            p.start()
        yield thread_log
        for p in patches:
            p.stop()

    def _run(self, tmp_path, threads=16):
        from rectify.core.align_command import run_align
        args = _make_args(threads=threads, tmp_path=tmp_path)
        run_align(args)

    def test_mappacbio_gets_all_threads_phase1(self, tmp_path, _patch_aligners):
        self._run(tmp_path, threads=16)
        assert _patch_aligners.get('mapPacBio') == [16], (
            "mapPacBio should be called with all 16 threads in phase 1"
        )

    def test_remaining_aligners_get_equal_shares(self, tmp_path, _patch_aligners):
        self._run(tmp_path, threads=16)
        # 4 remaining (minimap2, gapmm2, uLTRA, deSALT) → 16 // 4 = 4 each
        for name in ('minimap2', 'gapmm2', 'uLTRA', 'deSALT'):
            assert _patch_aligners.get(name) == [4], (
                f"{name} should receive 4 threads (16 // 4) in phase 2"
            )

    def test_thread_shares_scale_with_cpu_count(self, tmp_path, _patch_aligners):
        self._run(tmp_path, threads=8)
        assert _patch_aligners.get('mapPacBio') == [8]
        # 4 remaining → 8 // 4 = 2 each
        for name in ('minimap2', 'gapmm2', 'uLTRA', 'deSALT'):
            assert _patch_aligners.get(name) == [2], (
                f"{name} should receive 2 threads (8 // 4) with 8 total CPUs"
            )

    def test_no_parallel_flag_gives_all_threads_to_each(self, tmp_path, _patch_aligners):
        """Without --parallel-aligners every aligner runs sequentially with full threads."""
        from rectify.core.align_command import run_align
        args = _make_args(threads=16, parallel_aligners=False, tmp_path=tmp_path)
        run_align(args)
        for name in ('minimap2', 'mapPacBio', 'gapmm2', 'uLTRA', 'deSALT'):
            assert _patch_aligners.get(name) == [16], (
                f"{name} should get all 16 threads in sequential mode"
            )

    def test_without_mappacbio_all_aligners_run_in_parallel_pool(
        self, tmp_path, _patch_aligners
    ):
        """If mapPacBio is excluded phase 1 is skipped; others share threads equally."""
        from rectify.core.align_command import run_align
        args = _make_args(threads=16, tmp_path=tmp_path)
        # Remove mapPacBio from the aligner list by marking it unavailable
        with patch(
            'rectify.core.multi_aligner.check_aligner_available',
            side_effect=lambda p: p != args.mapPacBio_path,
        ):
            run_align(args)

        assert 'mapPacBio' not in _patch_aligners, "mapPacBio should not run if unavailable"
        # 4 remaining (minimap2, gapmm2, uLTRA, deSALT) each get 16 // 4 = 4
        for name in ('minimap2', 'gapmm2', 'uLTRA', 'deSALT'):
            assert _patch_aligners.get(name) == [4]


class TestPhase2MinimumOneThread:
    """Thread allocation never drops below 1 even with many aligners."""

    def test_minimum_thread_floor(self, tmp_path):
        """With 1 CPU and 4 remaining aligners each should still get 1 thread."""
        thread_log = {}

        def _runner(name):
            def _r(*args, threads=1, **kwargs):
                thread_log.setdefault(name, []).append(threads)
            return _r

        patches = [
            patch(f'rectify.core.multi_aligner.{fn}', side_effect=_runner(nm))
            for fn, nm in [
                ('run_minimap2', 'minimap2'), ('run_map_pacbio', 'mapPacBio'),
                ('run_gapmm2', 'gapmm2'), ('run_ultra', 'uLTRA'),
                ('run_desalt', 'deSALT'),
            ]
        ] + [patch('rectify.core.multi_aligner.check_aligner_available', return_value=True)]

        for p in patches:
            p.start()
        try:
            from rectify.core.align_command import run_align
            args = _make_args(threads=1, tmp_path=tmp_path)
            run_align(args)
        finally:
            for p in patches:
                p.stop()

        assert thread_log.get('mapPacBio') == [1]
        for name in ('minimap2', 'gapmm2', 'uLTRA', 'deSALT'):
            assert thread_log.get(name) == [1], (
                f"{name} should get at least 1 thread even with 1 total CPU"
            )
