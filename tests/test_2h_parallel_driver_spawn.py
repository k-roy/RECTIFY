"""ISSUE-025 — Module 2H's parallel driver must not depend on the fork start method.

``_run_parallel`` used to fill the module global ``_WORKER_POOL_STATE`` in the
parent and open a bare ``mp.Pool``; workers read the global.  That only works
under ``fork``.  On a ``spawn`` platform (macOS) every worker saw an empty dict,
raised ``KeyError('header')``, and ``correct`` logged "Module 2H junction
refinement failed (non-fatal, continuing)" — so 2H did NOTHING at >= 2 threads,
with no marker in any output.  The state now travels through the pool
initializer (``ctx.Pool(initializer=_init_worker_pool_state, initargs=...)``),
and a failure is a RuntimeError in the driver, an ERROR log plus a
``module_2h_failed`` row in the stats TSV in ``correct``.

These tests run the real ``refine_bam_junctions`` on the bundled upf1Δ mapPacBio
BAM (pool from all five aligner BAMs, bundled S. cerevisiae genome + GFF) with
``n_workers=2`` under an explicit ``spawn`` context and assert the refined count
and the written records equal the sequential run's.  Records — not raw bytes —
are compared because ``samtools sort`` writes each output's own path into its
@PG header line.

Author: Kevin R. Roy (agent, ISSUE-025)
"""

import multiprocessing as mp
import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_refiner as jr  # noqa: E402
from rectify.core.bam.processing_stats import ProcessingStats, write_stats_tsv  # noqa: E402

ALIGNERS = ['minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA']
ALIGNER_DIR = RECTIFY_ROOT / 'rectify' / 'data' / 'validation' / 'aligners'
INPUT_BAM = ALIGNER_DIR / 'validation_reads_upf1d.mapPacBio.bam'
# The sequential run on this bundle: 36 records, 20 with N-ops, 7 refined
# (the number the ISSUE-025 repro recorded at 1 thread).
EXPECTED_REFINED = 7


def _bundle_or_skip():
    bams = [ALIGNER_DIR / f'validation_reads_upf1d.{a}.bam' for a in ALIGNERS]
    if not INPUT_BAM.exists() or not all(b.exists() for b in bams):
        pytest.skip('bundled upf1d per-aligner BAMs not found')
    from rectify.data import get_bundled_genome_path, get_bundled_annotation_path
    g, a = get_bundled_genome_path('saccharomyces_cerevisiae'), get_bundled_annotation_path('saccharomyces_cerevisiae')
    if g is None or a is None:
        pytest.skip('bundled S. cerevisiae genome/annotation not available')
    return [str(b) for b in bams], g, a


@pytest.fixture(scope='module')
def context():
    bams, genome_path, annotation_path = _bundle_or_skip()
    from rectify.utils.genome import load_genome
    from rectify.core.consensus.consensus import load_annotated_junctions
    return dict(bams=bams, genome=load_genome(genome_path),
                annotated=load_annotated_junctions(str(annotation_path)))


def _records(path):
    with pysam.AlignmentFile(str(path)) as fh:
        return sorted(r.to_string() for r in fh.fetch(until_eof=True))


def _refine(ctx, out, n_workers):
    return jr.refine_bam_junctions(
        str(INPUT_BAM), str(out), ctx['bams'], ctx['annotated'], ctx['genome'],
        n_workers=n_workers, batch_size=5,     # 20 N-op reads -> 4 batches -> both workers busy
    )


@pytest.fixture(scope='module')
def sequential(context, tmp_path_factory):
    out = tmp_path_factory.mktemp('seq') / 'seq.bam'
    stats = _refine(context, out, 1)
    return stats, _records(out)


def test_sequential_baseline(sequential):
    stats, records = sequential
    assert stats['total'] == 36 and stats['n_op_reads'] == 20
    assert stats['refined'] == EXPECTED_REFINED
    assert len(records) == 36


def test_two_spawn_workers_match_sequential(context, sequential, tmp_path, monkeypatch):
    """The ISSUE-025 shape: explicit spawn, 2 workers.  Before the fix this raised
    KeyError('header') out of the pool."""
    monkeypatch.setenv('RECTIFY_2H_MP_START_METHOD', 'spawn')
    assert jr._get_refiner_mp_context().get_start_method() == 'spawn'
    out = tmp_path / 'spawn2.bam'
    stats = _refine(context, out, 2)
    assert stats['refined'] == EXPECTED_REFINED == sequential[0]['refined']
    assert stats == sequential[0]
    assert _records(out) == sequential[1]


def test_platform_default_two_workers_match_sequential(context, sequential, tmp_path, monkeypatch):
    """Whatever the platform default is (fork on Linux, spawn on macOS), the
    driver must give the sequential answer."""
    monkeypatch.delenv('RECTIFY_2H_MP_START_METHOD', raising=False)
    assert jr._get_refiner_mp_context().get_start_method() == mp.get_start_method()
    out = tmp_path / 'default2.bam'
    stats = _refine(context, out, 2)
    assert stats == sequential[0]
    assert _records(out) == sequential[1]


@pytest.mark.skipif('fork' not in mp.get_all_start_methods(), reason='fork unavailable here')
def test_fork_still_supported_when_requested(context, sequential, tmp_path, monkeypatch):
    """Linux production's zero-copy path: the initializer receives the state
    without pickling; the answer is the same."""
    monkeypatch.setenv('RECTIFY_2H_MP_START_METHOD', 'fork')
    out = tmp_path / 'fork2.bam'
    stats = _refine(context, out, 2)
    assert stats == sequential[0]
    assert _records(out) == sequential[1]


def test_invalid_start_method_falls_back_to_platform_default(monkeypatch, caplog):
    monkeypatch.setenv('RECTIFY_2H_MP_START_METHOD', 'teleport')
    with caplog.at_level('WARNING'):
        ctx = jr._get_refiner_mp_context()
    assert ctx.get_start_method() == mp.get_start_method()
    assert 'RECTIFY_2H_MP_START_METHOD' in caplog.text


def test_uninitialised_worker_fails_loudly(monkeypatch):
    """A worker without state must say so, not die on KeyError('header')."""
    monkeypatch.setattr(jr, '_WORKER_POOL_STATE', {})
    with pytest.raises(RuntimeError, match='ISSUE-025'):
        jr._refine_read_batch(['ignored'])


def test_initializer_rebuilds_the_header(monkeypatch):
    monkeypatch.setattr(jr, '_WORKER_POOL_STATE', {})
    header_dict = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': 1000}]}
    jr._init_worker_pool_state({'header_dict': header_dict, 'junctions_idx': {}, 'annotated_set': set(),
                                'genome': {}, 'kwargs': {}})
    assert isinstance(jr._WORKER_POOL_STATE['header'], pysam.AlignmentHeader)
    assert jr._WORKER_POOL_STATE['header'].to_dict()['SQ'] == header_dict['SQ']


def test_stats_tsv_carries_the_failure_marker(tmp_path):
    """`correct` records a survived 2H failure as the FIRST metric row."""
    st = ProcessingStats(total_reads_in_bam=3, module_2h_failed="KeyError: 'header'")
    out = tmp_path / 'x_stats.tsv'
    write_stats_tsv(st, str(out))
    rows = [l for l in out.read_text().splitlines() if l and not l.startswith('#')]
    assert rows[0].startswith('metric\t')
    assert rows[1].startswith('module_2h_failed\t1\t-\t')
    assert "KeyError: 'header'" in rows[1] and 'UNREFINED' in rows[1]
    assert 'module_2h_failed' in st.to_dict()
    # a clean run writes no such row
    clean = tmp_path / 'clean_stats.tsv'
    write_stats_tsv(ProcessingStats(total_reads_in_bam=3), str(clean))
    assert 'module_2h_failed' not in clean.read_text()


def test_merge_keeps_a_failure_from_either_side():
    a, b = ProcessingStats(), ProcessingStats(module_2h_failed='RuntimeError: x')
    a.merge(b)
    assert a.module_2h_failed == 'RuntimeError: x'
