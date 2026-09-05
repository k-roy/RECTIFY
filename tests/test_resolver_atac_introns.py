"""Opt-in paired AT-AC intron class (2026-09-04).

Yeast splices AT-AC junctions through its MAJOR spliceosome — Talkish et al.
2019 PLoS Genet 15:e1008249 ("Rapidly evolving protointrons in Saccharomyces
genomes revealed by a hungry spliceosome") report an AT-AC junction in
SUT635 — and in human AT-AC is the signature class of the U12-type minor
spliceosome. These tests pin:

* the index carries the four AT-AC arrays (format v3) and the canonical
  kinds do NOT see them (the default query path is byte-identical);
* ``ResolverConfig(atac=True)`` resolves a clip across an AT..AC intron that
  the default config refuses, on both strands;
* the class is PAIRED: an AT donor with an AG acceptor is not enumerable even
  under ``atac=True`` (that would be the motif-free extension planning/722
  refuted), and ``_donor_rank`` ranks a genuine AT-AC pair below GT and GC;
* ``canonical_in_class`` accepts AT..AC only when asked (Station C and every
  other default caller are unchanged);
* the default config is unchanged and the CLI knob is plumbed.
"""

import random

import pysam
import pytest

from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    _boundary_kinds,
    _donor_rank,
    _site_kinds,
    _site_kinds_atac,
    resolve_read,
)
from rectify.core.splice.overhang_informativeness import (
    canonical_in_class,
    is_canonical_junction,
    reset_counters,
)
from rectify.core.splice.splice_site_index import SpliceSiteIndex, _KINDS

# ---------------------------------------------------------------------------
# Synthetic genome with AT-AC introns:
#   plus  intron [1200, 1500): AT ........ AC
#   minus intron [2200, 2500): transcript AT..AC => forward genome
#                              GT (rc of AC) ........ AT (rc of AT)
# Random backbone (seed 909); default-mode refusal is asserted behaviourally
# (no rewrite across the planted junction), as in the prp18 tests.
# ---------------------------------------------------------------------------

GLEN = 4200
P_DON, P_ACC = 1200, 1500
M_ACC, M_DON = 2200, 2500
ATAC_KINDS = {'don_at_plus', 'acc_ac_plus', 'don_at_minus', 'acc_ac_minus'}


def _make_genome():
    rng = random.Random(909)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]

    def put(pos, s):
        seq[pos:pos + len(s)] = list(s)

    put(P_DON, 'AT')
    put(P_ACC - 2, 'AC')
    put(M_ACC, 'GT')          # rc('AC') — the AT-AC acceptor on the minus strand
    put(M_DON - 2, 'AT')      # rc('AT') — the AT-AC donor on the minus strand
    return ''.join(seq)


GENOME_SEQ = _make_genome()
GENOME = {'chrI': GENOME_SEQ}
# the planted class explodes near x far on a 4.2 kb contig whose W spans it
# whole (as the prp18 minus-strand test documents); real windows are far smaller
BIG = dict(max_candidates_per_clip=50_000)


@pytest.fixture(scope='module')
def index():
    return SpliceSiteIndex.build(GENOME)


@pytest.fixture(autouse=True)
def _fresh_counters():
    reset_counters()
    yield
    reset_counters()


def _header():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': GLEN}],
    })


def _read(name, query, cigar, ref_start, reverse=False):
    r = pysam.AlignedSegment(_header())
    r.query_name = name
    r.query_sequence = query
    r.flag = 16 if reverse else 0
    r.reference_id = 0
    r.reference_start = ref_start
    r.mapping_quality = 60
    r.cigartuples = cigar
    return r


def _plus_read():
    # left clip on the plus strand: near site = the AC acceptor at P_ACC
    clip = GENOME_SEQ[P_DON - 30:P_DON]
    return _read('atac_plus', clip + GENOME_SEQ[P_ACC:P_ACC + 60],
                 [(4, 30), (0, 60)], P_ACC)


def _minus_read():
    # right clip on the minus strand: near site = the (GT = rc AC) acceptor
    query = GENOME_SEQ[M_ACC - 60:M_ACC] + GENOME_SEQ[M_DON:M_DON + 30]
    return _read('atac_minus', query, [(0, 60), (4, 30)], M_ACC - 60,
                 reverse=True)


class TestIndex:
    def test_v3_carries_the_four_paired_arrays(self, index):
        kinds = {k.split('|', 1)[1] for k in index._arrays}
        assert ATAC_KINDS <= kinds
        assert ATAC_KINDS <= set(_KINDS)

    def test_planted_sites_land_in_the_right_arrays(self, index):
        def has(kind, pos):
            return pos in {int(x) for x in index.sites_in('chrI', kind, pos, pos + 1)}
        assert has('don_at_plus', P_DON) and has('acc_ac_plus', P_ACC)
        assert has('acc_ac_minus', M_ACC) and has('don_at_minus', M_DON)
        # ...and the canonical kinds do not see them
        assert not has('don_plus', P_DON) and not has('acc_plus', P_ACC)
        assert not has('acc_minus', M_ACC) and not has('don_minus', M_DON)

    def test_default_kind_selection_never_touches_atac_arrays(self):
        for side in ('L', 'R'):
            for strand in '+-':
                assert not set(_site_kinds(side, strand)) & ATAC_KINDS
                assert set(_site_kinds_atac(side, strand)) <= ATAC_KINDS
        for strand in '+-':
            assert not set(_boundary_kinds(strand)) & ATAC_KINDS


class TestClipResolution:
    def test_default_mode_does_not_place_the_atac_junction(self, index):
        for r in (_plus_read(), _minus_read()):
            before = (r.reference_start, list(r.cigartuples))
            changed = resolve_read(r, GENOME, index, ResolverConfig(**BIG),
                                   ResolverStats())
            if changed:
                assert not (r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')
                            or r.get_tag('XJ').startswith(f'{M_ACC}-{M_DON}:'))
            else:
                assert (r.reference_start, list(r.cigartuples)) == before

    def test_atac_mode_resolves_plus(self, index):
        r = _plus_read()
        assert resolve_read(r, GENOME, index, ResolverConfig(atac=True, **BIG),
                            ResolverStats())
        assert r.cigartuples == [(0, 30), (3, P_ACC - P_DON), (0, 60)]
        assert r.reference_start == P_DON - 30
        assert r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')

    def test_atac_mode_resolves_minus(self, index):
        r = _minus_read()
        assert resolve_read(r, GENOME, index, ResolverConfig(atac=True, **BIG),
                            ResolverStats())
        assert r.cigartuples == [(0, 60), (3, M_DON - M_ACC), (0, 30)]
        assert r.get_tag('XJ').startswith(f'{M_ACC}-{M_DON}:')

    def test_the_class_is_paired_an_at_donor_needs_an_ac_acceptor(self):
        # AT donor + AG acceptor: not a class, must stay non-enumerable
        seq = list(GENOME_SEQ)
        seq[P_ACC - 2:P_ACC] = 'AG'
        g = {'chrI': ''.join(seq)}
        idx = SpliceSiteIndex.build(g)
        clip = g['chrI'][P_DON - 30:P_DON]
        r = _read('at_ag', clip + g['chrI'][P_ACC:P_ACC + 60],
                  [(4, 30), (0, 60)], P_ACC)
        changed = resolve_read(r, g, idx, ResolverConfig(atac=True, **BIG),
                               ResolverStats())
        if changed:
            assert not r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')


class TestGrammar:
    def test_donor_rank_orders_gt_gc_atac(self):
        g = GENOME_SEQ
        assert _donor_rank(g, 'L', '+', P_DON, P_ACC) == 2
        assert _donor_rank(g, 'R', '-', M_ACC, M_DON) == 2
        s = list(g); s[P_DON:P_DON + 2] = 'GT'; s[P_ACC - 2:P_ACC] = 'AG'
        assert _donor_rank(''.join(s), 'L', '+', P_DON, P_ACC) == 0
        s[P_DON:P_DON + 2] = 'GC'
        assert _donor_rank(''.join(s), 'L', '+', P_DON, P_ACC) == 1

    def test_canonical_in_class_accepts_atac_only_when_asked(self):
        g = GENOME_SEQ
        for j in ((P_DON, P_ACC), (M_ACC, M_DON)):
            assert not is_canonical_junction(g, *j)
            assert is_canonical_junction(g, *j, atac=True)
            assert not canonical_in_class(g, *j)
            assert canonical_in_class(g, *j, atac=True)
        # a GT..AG junction is canonical either way
        s = list(g); s[P_DON:P_DON + 2] = 'GT'; s[P_ACC - 2:P_ACC] = 'AG'
        assert is_canonical_junction(''.join(s), P_DON, P_ACC)
        assert is_canonical_junction(''.join(s), P_DON, P_ACC, atac=True)


class TestDefaultsAndPlumbing:
    def test_default_is_off(self):
        assert ResolverConfig().atac is False

    def test_run_overhang_resolver_accepts_the_knob(self):
        import inspect
        from rectify.core.align.overhang_resolver import run_overhang_resolver
        assert 'atac' in inspect.signature(run_overhang_resolver).parameters

    def test_cli_flag_parses_on_align(self, tmp_path):
        import argparse
        from rectify.core.commands.align_command import create_align_parser
        root = argparse.ArgumentParser()
        sub = root.add_subparsers(dest='command')
        create_align_parser(sub)
        base = ['align', 'reads.fastq', '--genome', 'g.fa', '-o', str(tmp_path),
                '--aligners', 'minimap2']
        assert root.parse_args(base + ['--resolver-atac']).resolver_atac is True
        assert root.parse_args(base).resolver_atac is False
