"""A12 / A13 (dev/BUGS_TO_FIX.md, 2026-09-10): the overhang resolver on a
process pool, and ceiling abandonment as a reported correctness metric.

A12: `threads=N` used to be accepted and ignored (one core of resolver per
job however many slots were requested). The parallel driver scores batches
of records on a pool while one writer keeps input order — so the contract
under test is BYTE-IDENTICAL output and an identical tally across
`threads=1` and `threads>1`, batches deliberately smaller than the read set
so several tasks are in flight and one is a short tail.

A13: every clip abandoned on the candidate ceiling is a clip the rescue
never assessed. The count, the fraction and the per-contig split must
reach `ResolverStats.as_dict()` (the provenance JSON `align` writes), and
the end-of-run WARNING must say so and point at the knob (whether the
clips would have placed is data-dependent — real on human chr5, nothing on
yeast cDNA — so the warning never claims a lost-junction count).
"""

import argparse
import inspect
import logging

import pysam
import pytest

from rectify.core.align import overhang_resolver as ohr
from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    run_overhang_resolver,
)

from tests.test_overhang_resolver import (
    AMBIG_X, A_ACC, GENOME_SEQ, M_ACC, M_DON, P_ACC, P_DON, _header, _read,
)
from tests.test_resolver_candidate_ceiling import _clipped_bam, _planted_genome


def _mixed_bam(tmp_path, n_copies=6):
    """A stream that exercises every branch: plus/minus resolvable clips,
    a poly(A) clip (refused low-info), an ambiguous pair (rejected), a clean
    read (no clip), a secondary (passthrough) — repeated so the batch loop
    runs several full batches plus a tail."""
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chrI\n' + GENOME_SEQ + '\n')
    header = _header()
    reads = []
    for i in range(n_copies):
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        reads.append(_read(f'plus_{i}', clip + GENOME_SEQ[P_ACC:P_ACC + 60],
                           [(4, 30), (0, 60)], P_ACC, header=header))
        # minus intron [M_ACC, M_DON): aligned block ends at M_ACC, the
        # clip on the RIGHT carries the exon beyond the donor
        reads.append(_read(f'minus_{i}',
                           GENOME_SEQ[M_ACC - 60:M_ACC] + GENOME_SEQ[M_DON:M_DON + 30],
                           [(0, 60), (4, 30)], M_ACC - 60, reverse=True,
                           header=header))
        reads.append(_read(f'polya_{i}', 'A' * 30 + GENOME_SEQ[P_ACC:P_ACC + 60],
                           [(4, 30), (0, 60)], P_ACC, header=header))
        reads.append(_read(f'ambig_{i}', AMBIG_X + GENOME_SEQ[A_ACC:A_ACC + 60],
                           [(4, 30), (0, 60)], A_ACC, header=header))
        reads.append(_read(f'clean_{i}', GENOME_SEQ[100:190], [(0, 90)], 100,
                           header=header))
        sec = _read(f'sec_{i}', clip + GENOME_SEQ[P_ACC:P_ACC + 60],
                    [(4, 30), (0, 60)], P_ACC, header=header)
        sec.flag |= 256
        reads.append(sec)
    bam = tmp_path / 'in.bam'
    with pysam.AlignmentFile(bam, 'wb', header=header) as fh:
        for r in reads:
            fh.write(r)
    return bam, fasta, len(reads)


def _records(path):
    with pysam.AlignmentFile(path, 'rb', check_sq=False) as fh:
        pg = [p for p in fh.header.to_dict().get('PG', [])
              if p.get('ID') == 'rectify-overhang-resolver']
        return [r.to_string() for r in fh.fetch(until_eof=True)], pg


class TestParallelDriver:
    def test_parallel_output_is_byte_identical_to_serial(self, tmp_path):
        bam, fasta, n = _mixed_bam(tmp_path)
        ohr._BLOWUP_WARNED.clear()
        run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 'serial.bam'),
                              threads=1)
        serial = run_overhang_resolver.last_stats
        run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 'par.bam'),
                              threads=3, batch_size=7)   # 36 reads -> 6 tasks
        par = run_overhang_resolver.last_stats

        recs_s, pg_s = _records(tmp_path / 'serial.bam')
        recs_p, pg_p = _records(tmp_path / 'par.bam')
        assert len(recs_s) == n
        assert recs_p == recs_s                      # same records, same order
        assert par.as_dict() == serial.as_dict()      # same tally, every field
        # the run actually did something on every branch
        assert serial.resolved == 12 and serial.resolved_left == 6
        assert serial.resolved_right == 6
        assert serial.refused_low_info == 6
        assert serial.rejected_ambiguous == 6
        assert serial.passthrough_nonprimary == 6
        assert serial.reads == n
        assert 'threads=3' in pg_p[0]['CL'] and 'threads=1' in pg_s[0]['CL']

    @pytest.mark.parametrize('method', ['spawn', 'fork'])
    def test_both_start_methods_are_byte_identical(self, tmp_path, monkeypatch, method):
        """fork (the Linux default: workers inherit the parent's genome/index
        copy-on-write) and spawn (fresh interpreters that load their own)
        must produce the same records and tally."""
        import multiprocessing as _mp
        if method not in _mp.get_all_start_methods():
            pytest.skip(f'{method} unavailable here')
        monkeypatch.setenv('RECTIFY_RESOLVER_MP_START_METHOD', method)
        bam, fasta, n = _mixed_bam(tmp_path, n_copies=2)
        run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 's.bam'), threads=1)
        run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 'p.bam'),
                              threads=2, batch_size=5)
        assert _records(tmp_path / 'p.bam')[0] == _records(tmp_path / 's.bam')[0]
        assert run_overhang_resolver.last_stats.reads == n
        assert ohr._resolver_mp_context().get_start_method() == method

    def test_linux_defaults_to_fork_and_the_bam_knob_is_honoured(self, monkeypatch):
        import sys
        monkeypatch.delenv('RECTIFY_RESOLVER_MP_START_METHOD', raising=False)
        monkeypatch.delenv('RECTIFY_BAM_MP_START_METHOD', raising=False)
        expected = 'fork' if sys.platform.startswith('linux') else 'spawn'
        assert ohr._resolver_mp_context().get_start_method() == expected
        monkeypatch.setenv('RECTIFY_BAM_MP_START_METHOD', 'spawn')
        assert ohr._resolver_mp_context().get_start_method() == 'spawn'
        monkeypatch.setenv('RECTIFY_RESOLVER_MP_START_METHOD', 'nonsense')
        assert ohr._resolver_mp_context().get_start_method() == 'spawn'

    def test_parallel_with_a_single_short_batch(self, tmp_path):
        # batch larger than the stream: exactly one task, still ordered
        bam, fasta, n = _mixed_bam(tmp_path, n_copies=1)
        run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 'a.bam'), threads=1)
        run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 'b.bam'),
                              threads=2, batch_size=1000)
        assert _records(tmp_path / 'a.bam')[0] == _records(tmp_path / 'b.bam')[0]

    def test_threads_below_two_take_the_serial_path(self, tmp_path, monkeypatch):
        bam, fasta, _ = _mixed_bam(tmp_path, n_copies=1)
        called = []
        monkeypatch.setattr(ohr, '_run_parallel',
                            lambda *a, **k: called.append(1))
        for t in (0, 1, None):
            run_overhang_resolver(str(bam), str(fasta), str(tmp_path / f't{t}.bam'),
                                  threads=t)
        assert called == []

    def test_no_more_single_threaded_warning(self, tmp_path, caplog):
        bam, fasta, _ = _mixed_bam(tmp_path, n_copies=1)
        with caplog.at_level(logging.WARNING, logger=ohr.logger.name):
            run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 'o.bam'),
                                  threads=2, batch_size=4)
        assert not [r for r in caplog.records if 'SINGLE-THREADED' in r.getMessage()]


class TestAbandonmentIsAMetric:
    def test_per_contig_split_and_fraction_reach_the_stats_dict(self, tmp_path):
        ohr._BLOWUP_WARNED.clear()
        seq = _planted_genome(tmp_path)
        bam = _clipped_bam(tmp_path, seq)
        run_overhang_resolver(str(bam), str(tmp_path / 'g.fa'), str(tmp_path / 'o.bam'),
                              config=ResolverConfig(max_candidates_per_clip=5))
        d = run_overhang_resolver.last_stats.as_dict()
        assert d['refused_candidate_blowup'] == 1
        assert d['clips_assessed'] == 1
        assert d['abandoned_frac'] == 1.0
        assert d['blowup_by_contig'] == {'chrI': 1}
        ceiling, n_cand, window = d['blowup_first']['chrI']
        assert ceiling == 5 and n_cand == 6 and window > 0

    def test_kwarg_overrides_the_ceiling_without_a_config(self, tmp_path):
        ohr._BLOWUP_WARNED.clear()
        seq = _planted_genome(tmp_path)
        bam = _clipped_bam(tmp_path, seq)
        run_overhang_resolver(str(bam), str(tmp_path / 'g.fa'), str(tmp_path / 'o.bam'),
                              max_candidates_per_clip=5)
        assert run_overhang_resolver.last_stats.refused_candidate_blowup == 1
        run_overhang_resolver(str(bam), str(tmp_path / 'g.fa'), str(tmp_path / 'o2.bam'))
        assert run_overhang_resolver.last_stats.refused_candidate_blowup == 0

    def test_summary_warning_names_the_consequence(self, tmp_path, caplog):
        ohr._BLOWUP_WARNED.clear()
        seq = _planted_genome(tmp_path)
        bam = _clipped_bam(tmp_path, seq)
        with caplog.at_level(logging.WARNING, logger=ohr.logger.name):
            run_overhang_resolver(str(bam), str(tmp_path / 'g.fa'),
                                  str(tmp_path / 'o.bam'), max_candidates_per_clip=5)
        msgs = [r.getMessage() for r in caplog.records]
        summary = [m for m in msgs if 'ABANDONED on the candidate ceiling' in m]
        assert len(summary) == 1
        assert 'NOT ASSESSED' in summary[0]
        assert '(100.0%)' in summary[0]
        assert 'chrI=1' in summary[0]
        assert '--resolver-candidate-ceiling' in summary[0]

    def test_parallel_parent_warns_once_per_contig_and_merges_the_split(
            self, tmp_path, caplog):
        ohr._BLOWUP_WARNED.clear()
        seq = _planted_genome(tmp_path)
        bam = _clipped_bam(tmp_path, seq)
        with caplog.at_level(logging.WARNING, logger=ohr.logger.name):
            run_overhang_resolver(str(bam), str(tmp_path / 'g.fa'),
                                  str(tmp_path / 'o.bam'), threads=2,
                                  max_candidates_per_clip=5)
        st = run_overhang_resolver.last_stats
        assert st.refused_candidate_blowup == 1
        assert st.blowup_by_contig == {'chrI': 1}
        per_contig = [r for r in caplog.records
                      if "candidate blow-up on contig 'chrI'" in r.getMessage()]
        assert len(per_contig) == 1

    def test_merge_sums_counters_and_keeps_the_first_example(self):
        a = ResolverStats(reads=2, resolved=1, extra={'x': 1},
                          blowup_by_contig={'chrI': 1}, blowup_first={'chrI': (1, 2, 3)})
        b = ResolverStats(reads=3, resolved=0, extra={'x': 2, 'y': 5},
                          blowup_by_contig={'chrI': 2, 'chrII': 1},
                          blowup_first={'chrI': (9, 9, 9), 'chrII': (4, 5, 6)})
        a.merge(b)
        assert (a.reads, a.resolved) == (5, 1)
        assert a.extra == {'x': 3, 'y': 5}
        assert a.blowup_by_contig == {'chrI': 3, 'chrII': 1}
        assert a.blowup_first == {'chrI': (1, 2, 3), 'chrII': (4, 5, 6)}


class TestMoveTagIsNotShared:
    """A10: the resolver's move-family tag was `XB`, which the ONT cDNA
    pipeline also writes (strand split) and the consensus sidecar restore
    puts back on every cDNA read — so the resolver's record vanished from
    the final BAM. The tag is now `XE`; pin it against every other writer."""

    RESOLVER_TAGS = ('XJ', 'XE')

    def test_resolver_tags_collide_with_no_cdna_or_cma_tag(self):
        from rectify.core.consensus.consensus import _CDNA_COMMENT_TAGS
        from rectify.core.multialign.cma_schema import READ_INTRINSIC_TAGS
        for t in self.RESOLVER_TAGS:
            assert t not in _CDNA_COMMENT_TAGS, t
            assert t not in READ_INTRINSIC_TAGS, t

    def test_no_other_module_writes_the_resolver_tags(self):
        import re
        from pathlib import Path
        import rectify
        root = Path(rectify.__file__).parent
        pat = re.compile(r"set_tag\(\s*['\"](X[JE])['\"]")
        hits = {}
        for py in root.rglob('*.py'):
            if py.name == 'overhang_resolver.py':
                continue
            for m in pat.finditer(py.read_text()):
                hits.setdefault(m.group(1), []).append(str(py.relative_to(root)))
        assert hits == {}, hits

    def test_resolver_writes_xe_not_xb(self, tmp_path):
        bam, fasta, _ = _mixed_bam(tmp_path, n_copies=1)
        run_overhang_resolver(str(bam), str(fasta), str(tmp_path / 'o.bam'))
        with pysam.AlignmentFile(tmp_path / 'o.bam', 'rb', check_sq=False) as fh:
            recs = list(fh.fetch(until_eof=True))
        assert not any(r.has_tag('XB') for r in recs)
        # the clip placements in this fixture carry XJ; XE is the arbiter's
        # tag and is written only when a junction MOVES, which none do here
        assert sum(r.has_tag('XJ') for r in recs) == 2


class TestKnobPlumbing:
    def test_align_and_run_all_parsers_accept_the_ceiling(self):
        from rectify.core.commands.align_command import (
            _RESOLVER_CEILING_DEFAULT, create_align_parser)
        from rectify.core.commands.run_command import create_run_parser
        assert _RESOLVER_CEILING_DEFAULT == ResolverConfig().max_candidates_per_clip
        root = argparse.ArgumentParser()
        sub = root.add_subparsers(dest='command')
        create_align_parser(sub)
        create_run_parser(sub)
        a = root.parse_args(['align', 'r.fq', '--genome', 'g.fa', '-o', 'out',
                             '--resolver-candidate-ceiling', '20000'])
        assert a.resolver_candidate_ceiling == 20000
        a = root.parse_args(['align', 'r.fq', '--genome', 'g.fa', '-o', 'out'])
        assert a.resolver_candidate_ceiling is None
        r = root.parse_args(['run-all', 'r.fq', '--genome', 'g.fa', '-o', 'out',
                             '--resolver-candidate-ceiling', '20000'])
        assert r.resolver_candidate_ceiling == 20000

    def test_run_all_threads_the_ceiling_into_run_align(self):
        from rectify.core.commands.run import stages
        assert 'resolver_candidate_ceiling' in inspect.signature(
            stages._run_alignment).parameters
        src = inspect.getsource(stages._run_alignment)
        assert 'resolver_candidate_ceiling=resolver_candidate_ceiling' in src

    def test_run_multi_aligner_forwards_every_resolver_knob(self):
        """A14: the multi_aligner branch used to call the resolver with
        `threads` only, dropping --max-intron / --no-resolver-atac."""
        from rectify.core.align import multi_aligner as ma
        params = inspect.signature(ma.run_multi_aligner).parameters
        for k in ('resolver_acceptor_classes', 'resolver_atac',
                  'resolver_candidate_ceiling', 'max_intron'):
            assert k in params
        src = inspect.getsource(ma.run_multi_aligner)
        call = src[src.index("results['overhang_resolver'] = run_overhang_resolver("):]
        call = call[:call.index('\n                )')]
        for k in ('max_intron=max_intron', 'atac=resolver_atac',
                  'acceptor_classes=resolver_acceptor_classes',
                  'max_candidates_per_clip=resolver_candidate_ceiling'):
            assert k in call, k
