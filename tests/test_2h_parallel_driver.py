"""ISSUE-025 — Module 2H parallel driver: complements ``test_2h_parallel_driver_spawn.py``.

The spawn file pins the start-method contract itself (spawn / fork / platform
default at 2 workers give the sequential answer, table-less scoring) and
unit-tests the ``module_2h_failed`` ProcessingStats field.  This file covers
what that leaves open:

* 1-vs-2-thread per-read POS/CIGAR identity WITH the bundled yeast DRS penalty
  tables — the tables ``correct`` actually selects for this genome, i.e. the
  scoring path production runs (ISSUE-005: table-less scoring is a materially
  different algorithm).
* The driver's loud-failure contract: a pool failure (the pre-fix
  ``KeyError('header')`` shape, via a fake context) and a real worker raising
  (fork pool) both surface as a ``RuntimeError`` naming the driver, with an
  ERROR log — never a half-empty result.
* ``correct`` end to end (in-process ``correct_command.run`` on the upf1Δ
  mapPacBio BAM): a 2H failure writes the ``module_2h_failed`` row into the
  stats TSV and logs at ERROR, stays non-fatal by default,
  ``RECTIFY_2H_FAILURE_FATAL=1`` turns it into the CLI's ``SystemExit(1)``, and
  a clean run writes no such row.

Refined count on this input: 7 at 1 thread, with or without the tables.
"""
from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import os
from pathlib import Path

import pysam
import pytest

from rectify.core.commands import correct_command
from rectify.core.consensus.consensus import load_annotated_junctions
from rectify.core.splice import junction_refiner as jr
from rectify.data import get_bundled_annotation_path, get_bundled_genome_path
from rectify.utils.genome import load_genome

ALIGNERS = ['minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA']
VAL_DIR = Path(jr.__file__).resolve().parents[2] / 'data' / 'validation' / 'aligners'
INPUT_BAM = VAL_DIR / 'validation_reads_upf1d.mapPacBio.bam'
ALIGNER_BAMS = [VAL_DIR / f'validation_reads_upf1d.{a}.bam' for a in ALIGNERS]

EXPECTED_N_OP_READS = 20
EXPECTED_REFINED = 7

# 20 N-op reads -> 5 batches, so both workers really receive work (the default
# batch_size=200 would hand everything to one worker).
BATCH_SIZE = 4

pytestmark = pytest.mark.skipif(
    not INPUT_BAM.exists() or not all(b.exists() for b in ALIGNER_BAMS),
    reason='bundled upf1Δ per-aligner validation BAMs not present',
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def genome_path() -> Path:
    p = get_bundled_genome_path('saccharomyces_cerevisiae')
    if p is None or not Path(p).exists():
        pytest.skip('bundled yeast genome not available')
    return Path(p)


@pytest.fixture(scope='module')
def annotation_path() -> Path:
    p = get_bundled_annotation_path('saccharomyces_cerevisiae')
    if p is None or not Path(p).exists():
        pytest.skip('bundled yeast annotation not available')
    return Path(p)


@pytest.fixture(scope='module')
def genome(genome_path):
    return load_genome(str(genome_path))


@pytest.fixture(scope='module')
def annotated(annotation_path):
    return load_annotated_junctions(str(annotation_path))


@pytest.fixture(scope='module')
def penalty_tables(genome_path):
    """Exactly what ``correct`` selects for this genome: organism detected from
    the bundled FASTA, protocol = DRS.  Must resolve to real tables, or this
    file silently degrades into the spawn file's table-less coverage."""
    jpt, spt = correct_command.select_penalty_tables({'genome_path': str(genome_path)})
    assert jpt and Path(jpt).exists(), 'bundled yeast DRS junction penalty table not selected'
    return jpt, spt


def _refine(out_bam: Path, n_threads: int, genome, annotated, penalty_tables) -> dict:
    jpt, spt = penalty_tables
    return jr.refine_bam_junctions(
        input_bam=str(INPUT_BAM),
        output_bam=str(out_bam),
        aligner_bams=[str(b) for b in ALIGNER_BAMS],
        annotated_junctions=annotated,
        genome=genome,
        sort_and_index=True,
        penalty_table_path=jpt,
        str_penalty_table_path=spt,
        n_workers=n_threads,
        batch_size=BATCH_SIZE,
    )


def _alignments(bam: Path) -> list:
    """Sorted (name, flag, chrom, POS, CIGAR) per record.  The input carries a
    multi-mapper (``fa816f03…`` flag 16 on chrVII AND chrXII), so records are
    compared as a sorted multiset, not keyed by (name, flag)."""
    with pysam.AlignmentFile(str(bam), 'rb') as fh:
        return sorted((r.query_name, r.flag, r.reference_name, r.reference_start, r.cigarstring)
                      for r in fh.fetch(until_eof=True))


@pytest.fixture(scope='module')
def sequential(tmp_path_factory, genome, annotated, penalty_tables):
    out = tmp_path_factory.mktemp('seq') / 'refined_t1.bam'
    stats = _refine(out, 1, genome, annotated, penalty_tables)
    assert stats['n_op_reads'] == EXPECTED_N_OP_READS
    assert stats['refined'] == EXPECTED_REFINED
    assert stats['errors'] == 0
    return out, stats


# ---------------------------------------------------------------------------
# 1 thread vs 2 threads, platform default start method, production tables
# ---------------------------------------------------------------------------

def test_two_threads_match_sequential_with_penalty_tables(tmp_path, sequential, genome,
                                                          annotated, penalty_tables):
    seq_bam, seq_stats = sequential
    out = tmp_path / 'refined_t2.bam'
    stats = _refine(out, 2, genome, annotated, penalty_tables)
    assert stats == seq_stats
    assert stats['refined'] == EXPECTED_REFINED
    a, b = _alignments(seq_bam), _alignments(out)
    assert len(a) == len(b) == 36
    diffs = [(x, y) for x, y in zip(a, b) if x != y]
    assert not diffs, 'POS/CIGAR differ between 1 and 2 threads:\n' + '\n'.join(
        f'  seq={x}\n  par={y}' for x, y in diffs[:10])


# ---------------------------------------------------------------------------
# Loud failure — driver level
# ---------------------------------------------------------------------------

class _BoomPool:
    def __init__(self, *a, **k):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap_unordered(self, fn, batches):
        raise KeyError('header')          # the exact pre-fix worker failure


class _BoomContext:
    Pool = _BoomPool

    @staticmethod
    def get_start_method():
        return 'boom'


def test_driver_pool_failure_raises_and_logs_error(tmp_path, monkeypatch, caplog,
                                                   genome, annotated, penalty_tables):
    monkeypatch.setattr(jr, '_get_refiner_mp_context', lambda: _BoomContext())
    caplog.set_level(logging.ERROR, logger=jr.logger.name)
    with pytest.raises(RuntimeError, match='Module 2H parallel driver failed') as ei:
        _refine(tmp_path / 'boom.bam', 2, genome, annotated, penalty_tables)
    assert isinstance(ei.value.__cause__, KeyError)
    errs = [r for r in caplog.records
            if r.levelno == logging.ERROR and r.name == jr.logger.name]
    assert errs and 'parallel driver failed' in errs[0].getMessage()


def _boom_batch(_sam_strings):
    raise ValueError('worker exploded (ISSUE-025 test)')


@pytest.mark.skipif('fork' not in mp.get_all_start_methods(),
                    reason='fork start method unavailable')
def test_driver_worker_exception_raises(tmp_path, monkeypatch, genome, annotated,
                                        penalty_tables):
    """A real pool whose worker raises: the driver re-raises instead of
    returning a partial result.  Fork, so the monkeypatched worker is what the
    child runs (a spawned child would re-import the original)."""
    monkeypatch.setattr(jr, '_get_refiner_mp_context', lambda: mp.get_context('fork'))
    monkeypatch.setattr(jr, '_refine_read_batch', _boom_batch)
    with pytest.raises(RuntimeError, match='Module 2H parallel driver failed') as ei:
        _refine(tmp_path / 'boom_fork.bam', 2, genome, annotated, penalty_tables)
    assert isinstance(ei.value.__cause__, ValueError)


# ---------------------------------------------------------------------------
# Loud failure — `correct` end to end
# ---------------------------------------------------------------------------

def _correct_args(bam: Path, genome_path: Path, annotation_path: Path,
                  out_tsv: Path) -> argparse.Namespace:
    """Parse a real `rectify correct` command line so no default can drift."""
    parser = argparse.ArgumentParser()
    correct_command.create_correct_parser(parser.add_subparsers(dest='command'))
    argv = ['correct', str(bam), '--genome', str(genome_path),
            '--annotation', str(annotation_path), '-o', str(out_tsv),
            '--emit-merged-tsv',       # as the upf1Δ fixture runs it; keeps <out>_stats.tsv beside <out>
            '--threads', '1']
    for a, b in zip(ALIGNERS, ALIGNER_BAMS):
        argv += ['--aligner-bams', f'{a}:{b}']
    return parser.parse_args(argv)


def _run_correct_restoring_env(args):
    """``correct_command.run`` pins thread env vars (LOKY_MAX_CPU_COUNT, ...)
    in-process; undo that so it cannot leak into later tests (ISSUE-022)."""
    before = dict(os.environ)
    try:
        correct_command.run(args)
    finally:
        for k in set(os.environ) - set(before):
            del os.environ[k]
        for k, v in before.items():
            if os.environ.get(k) != v:
                os.environ[k] = v


def _stats_rows(out_tsv: Path) -> dict:
    stats_tsv = out_tsv.with_name(out_tsv.name.replace('.tsv', '_stats.tsv'))
    assert stats_tsv.exists(), f'stats TSV missing: {stats_tsv}'
    rows = {}
    for line in stats_tsv.read_text().splitlines():
        if line.startswith('#') or line.startswith('metric\t') or not line.strip():
            continue
        parts = line.split('\t')
        rows[parts[0]] = parts
    return rows


def test_correct_records_module_2h_failure_and_logs_error(tmp_path, monkeypatch, caplog,
                                                          genome_path, annotation_path):
    def _boom(**kwargs):
        raise RuntimeError('simulated 2H driver failure (ISSUE-025 test)')

    monkeypatch.setattr(jr, 'refine_bam_junctions', _boom)
    monkeypatch.delenv('RECTIFY_2H_FAILURE_FATAL', raising=False)
    caplog.set_level(logging.ERROR, logger=correct_command.__name__)
    out_tsv = tmp_path / 'corrected.tsv'
    _run_correct_restoring_env(_correct_args(INPUT_BAM, genome_path, annotation_path, out_tsv))

    assert out_tsv.exists(), 'non-fatal by default: correct must still finish'
    rows = _stats_rows(out_tsv)
    assert 'module_2h_failed' in rows, f'stats TSV lacks module_2h_failed: {sorted(rows)[:10]}'
    assert 'RuntimeError: simulated 2H driver failure' in rows['module_2h_failed'][3]
    errs = [r for r in caplog.records
            if r.levelno == logging.ERROR and r.name == correct_command.__name__
            and 'Module 2H junction refinement FAILED' in r.getMessage()]
    assert errs, 'no ERROR-level Module 2H failure log'
    assert 'simulated 2H driver failure' in errs[0].getMessage()


def test_correct_module_2h_failure_fatal_opt_in(tmp_path, monkeypatch, caplog,
                                                genome_path, annotation_path):
    """With RECTIFY_2H_FAILURE_FATAL=1 the 2H exception is re-raised and the
    CLI's outer handler turns it into ``sys.exit(1)``."""
    def _boom(**kwargs):
        raise RuntimeError('simulated 2H driver failure (fatal test)')

    monkeypatch.setattr(jr, 'refine_bam_junctions', _boom)
    monkeypatch.setenv('RECTIFY_2H_FAILURE_FATAL', '1')
    caplog.set_level(logging.ERROR, logger=correct_command.__name__)
    out_tsv = tmp_path / 'corrected_fatal.tsv'
    with pytest.raises(SystemExit) as ei:
        _run_correct_restoring_env(_correct_args(INPUT_BAM, genome_path, annotation_path, out_tsv))
    assert ei.value.code == 1
    msgs = [r.getMessage() for r in caplog.records
            if r.levelno == logging.ERROR and r.name == correct_command.__name__]
    assert any('Module 2H junction refinement FAILED' in m for m in msgs)
    assert any('simulated 2H driver failure (fatal test)' in m for m in msgs)


def test_correct_success_leaves_module_2h_failed_empty(tmp_path, genome_path, annotation_path):
    """The stats key is a failure marker only: absent from the TSV when 2H ran."""
    out_tsv = tmp_path / 'corrected_ok.tsv'
    _run_correct_restoring_env(_correct_args(INPUT_BAM, genome_path, annotation_path, out_tsv))
    rows = _stats_rows(out_tsv)
    assert 'module_2h_failed' not in rows
    # T0 4f9102d F3: the refine_bam_junctions stats reach the TSV as module_2h_<key>
    # rows (percent '-'), so the 2H counters are readable offline.
    assert rows['module_2h_n_op_reads'][1:3] == [str(EXPECTED_N_OP_READS), '-']
    assert rows['module_2h_refined'][1:3] == [str(EXPECTED_REFINED), '-']
    assert 'module_2h_noncanon_destination_refused' in rows
