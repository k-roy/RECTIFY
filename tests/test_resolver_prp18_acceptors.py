"""Opt-in Prp18-class acceptor extension (planning/722b).

Roy et al. 2023 NAR (gkad968) measured genome-wide activation of non-YAG
3'SS in prp18 mutants — classes BG (TG/CG/GG) and non-G HAU (AT) — and the
same junctions accumulate in plain upf1Δ. These tests pin:

* the index always carries the extended acceptor arrays, and the ``*_all``
  union kinds expose them while the canonical kinds do NOT (the default
  query path is byte-identical to the planning/720-measured configuration);
* ``ResolverConfig(acceptor_classes='prp18')`` resolves a clip across a
  GT..TG (BG-class) intron that canonical mode refuses, on both strands;
* donors are NOT extended — an AT donor stays non-enumerable even under
  'prp18' (4/1,833 published alt-3'SS junctions had non-canonical donors;
  the extension is deliberately acceptor-only);
* the default config is unchanged.
"""

import random

import pysam
import pytest

from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    resolve_read,
)
from rectify.core.splice.overhang_informativeness import reset_counters
from rectify.core.splice.splice_site_index import SpliceSiteIndex

# ---------------------------------------------------------------------------
# Synthetic genome with Prp18-class introns:
#   plus  intron [1200, 1500): GT ........ TG   (BG-class acceptor)
#   minus intron [2200, 2500): transcript GT..TG => forward genome
#                              CA (acc, rc of TG) ........ AC (donor)
# The backbone is random with seed 727; canonical-mode refusal is asserted
# behaviorally (no rewrite), not by absence of candidates — incidental AG
# sites may exist in the backbone, but none can match the planted exon
# sequence, so canonical mode cannot accept a placement.
# ---------------------------------------------------------------------------

GLEN = 4200
P_DON, P_ACC = 1200, 1500
M_ACC, M_DON = 2200, 2500


def _make_genome():
    rng = random.Random(727)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]

    def put(pos, s):
        seq[pos:pos + len(s)] = list(s)

    put(P_DON, 'GT')
    put(P_ACC - 2, 'TG')      # BG-class acceptor (plus strand)
    put(M_ACC, 'CA')          # rc('TG') — BG-class acceptor on minus strand
    put(M_DON - 2, 'AC')      # canonical minus donor (transcript GT)
    return ''.join(seq)


GENOME_SEQ = _make_genome()
GENOME = {'chrI': GENOME_SEQ}


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


class TestIndexExtension:
    def test_ext_arrays_are_built_and_gated_behind_all_kinds(self, index):
        # canonical kind must NOT see the TG acceptor at the plus intron end
        assert P_ACC not in set(
            int(x) for x in index.sites_in('chrI', 'acc_plus', P_ACC - 1, P_ACC + 1))
        # the union kind must
        assert P_ACC in set(
            int(x) for x in index.sites_in('chrI', 'acc_plus_all', P_ACC - 1, P_ACC + 1))

    def test_minus_strand_revcomp_mapping(self, index):
        # forward-genome CA at the minus intron START = transcript TG acceptor
        assert M_ACC not in set(
            int(x) for x in index.sites_in('chrI', 'acc_minus', M_ACC, M_ACC + 1))
        assert M_ACC in set(
            int(x) for x in index.sites_in('chrI', 'acc_minus_all', M_ACC, M_ACC + 1))

    def test_kind_vocabulary_is_exactly_the_prp18_set(self, index):
        kinds = {k.split('|', 1)[1] for k in index._arrays}
        # 2026-09-04: fired, as designed, when the paired AT-AC class landed
        # (format v3, ResolverConfig.atac — see test_resolver_atac_introns.py).
        assert kinds == {'don_gt_plus', 'don_gc_plus', 'acc_plus', 'don_minus',
                         'acc_minus', 'acc_plus_ext', 'acc_minus_ext',
                         'don_at_plus', 'acc_ac_plus', 'don_at_minus', 'acc_ac_minus'}, (
            'the splice-site kind vocabulary changed — planning/722b pricing '
            'and the alpha calibration need review')


class TestClipResolution:
    """The behavioral pin: same read, same index — the config decides."""

    def _plus_read(self):
        # left clip on the plus strand: near site = the (TG) acceptor
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        return _read('bg_plus', query, [(4, 30), (0, 60)], P_ACC)

    def test_canonical_mode_does_not_place_the_bg_junction(self):
        r = self._plus_read()
        before = (r.reference_start, list(r.cigartuples))
        stats = ResolverStats()
        changed = resolve_read(r, GENOME, SpliceSiteIndex.build(GENOME),
                               ResolverConfig(), stats)
        # Whatever else happened, the read must not have been placed across
        # the BG-class junction (canonical mode cannot enumerate it).
        if changed:
            assert not r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')
        else:
            assert (r.reference_start, list(r.cigartuples)) == before

    def test_prp18_mode_resolves_the_bg_junction_plus(self, index):
        r = self._plus_read()
        stats = ResolverStats()
        cfg = ResolverConfig(acceptor_classes='prp18')
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(0, 30), (3, P_ACC - P_DON), (0, 60)]
        assert r.reference_start == P_DON - 30
        assert r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')

    def test_prp18_mode_resolves_the_bg_junction_minus(self, index):
        # right clip on the minus strand: near site = the (CA=rc TG) acceptor
        query = GENOME_SEQ[M_ACC - 60:M_ACC] + GENOME_SEQ[M_DON:M_DON + 30]
        r = _read('bg_minus', query, [(0, 60), (4, 30)], M_ACC - 60,
                  reverse=True)
        stats = ResolverStats()
        # The default 2000-candidate ceiling fires on this 4.2 kb synthetic
        # contig: extended acceptors are ~27% of positions and W spans the
        # whole contig, so near x far explodes — a live demonstration of the
        # planning/722b x4.8 price on a dense mini-genome, not a production
        # concern (real W windows are information-bounded far below 4 kb).
        cfg = ResolverConfig(acceptor_classes='prp18',
                             max_candidates_per_clip=50_000)
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(0, 60), (3, M_DON - M_ACC), (0, 30)]

    def test_donors_are_not_extended(self):
        # AT donor + TG acceptor: even 'prp18' mode must refuse — the
        # extension is acceptor-only by design.
        seq = list(GENOME_SEQ)
        seq[P_DON:P_DON + 2] = 'AT'
        g = {'chrI': ''.join(seq)}
        idx = SpliceSiteIndex.build(g)
        clip = g['chrI'][P_DON - 30:P_DON]
        query = clip + g['chrI'][P_ACC:P_ACC + 60]
        r = _read('at_donor', query, [(4, 30), (0, 60)], P_ACC)
        stats = ResolverStats()
        changed = resolve_read(r, g, idx, ResolverConfig(
            acceptor_classes='prp18'), stats)
        if changed:
            assert not r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')


class TestDefaults:
    def test_default_acceptor_classes_is_canonical(self):
        assert ResolverConfig().acceptor_classes == 'canonical'

    def test_run_overhang_resolver_accepts_the_knob(self):
        import inspect
        from rectify.core.align.overhang_resolver import run_overhang_resolver
        assert 'acceptor_classes' in inspect.signature(
            run_overhang_resolver).parameters
