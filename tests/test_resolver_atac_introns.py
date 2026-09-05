"""Paired AT-AC intron class — ON BY DEFAULT since 2026-09-05.

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
* the arbiter's two-boundary DISCOVERY paths are wired too (2026-09-05): the
  Case B1 intron-length-D -> N snap and both mismatch-flagged linear rescues
  (B2 right, B3 left) each run AT-AC as a separate paired pass;
* the default is ON (Kevin's 2026-09-05 call, every organism) and
  ``--no-resolver-atac`` / ``atac=False`` reproduces the pre-flip space.

2026-09-05 NOTE ON THE CONTRAST ARMS: these tests contrast the class ON against
the class OFF. Until the flip, "off" was spelled ``ResolverConfig()`` and "on"
``ResolverConfig(atac=True)``; now "off" must be spelled ``atac=False``. Only
the spelling changed — every assertion below still pins the same behaviour, so
the pre-flip evidence is preserved rather than deleted.
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
        """Unchanged across the 2026-09-05 default flip, on purpose: AT-AC is a
        SECOND pass (`_site_kinds_atac`), never a union into the canonical
        selectors, so these must stay disjoint however the default is set."""
        for side in ('L', 'R'):
            for strand in '+-':
                assert not set(_site_kinds(side, strand)) & ATAC_KINDS
                assert set(_site_kinds_atac(side, strand)) <= ATAC_KINDS
        for strand in '+-':
            assert not set(_boundary_kinds(strand)) & ATAC_KINDS


class TestClipResolution:
    def test_disabled_mode_does_not_place_the_atac_junction(self, index):
        # was `test_default_mode_...` with a bare ResolverConfig() until the
        # 2026-09-05 default flip; the behaviour pinned is identical.
        for r in (_plus_read(), _minus_read()):
            before = (r.reference_start, list(r.cigartuples))
            changed = resolve_read(r, GENOME, index,
                                   ResolverConfig(atac=False, **BIG),
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


class TestArbiterDiscoveryPaths:
    """The arbiter's two-boundary DISCOVERY paths — the ones that choose BOTH
    intron boundaries from the index rather than shifting an existing junction.
    Until 2026-09-05 these enumerated GT/GC..AG only, so an AT-AC intron the
    aligner had encoded as a deletion (Case B1) or smeared over linearly
    (Case B2/B3) could not be recovered even with ``atac=True``.

    Mirrors the GT..AG fixtures in ``tests/test_overhang_resolver.py``:
    ``test_dop_converted_to_snapped_intron`` (L612), the B2 right-side storm
    (L623) and the B3 head-storm mirror (L665)."""

    INTRON = P_ACC - P_DON

    @staticmethod
    def _junction(read):
        ref = read.reference_start
        for op, ln in read.cigartuples:
            if op == 3:
                return ref, ref + ln
            if op in (0, 2, 3, 7, 8):
                ref += ln
        return None

    @staticmethod
    def _run(g, read, **cfg_kw):
        genome = {'chrI': g}
        stats = ResolverStats()
        changed = resolve_read(read, genome, SpliceSiteIndex.build(genome),
                               ResolverConfig(alpha=0.01, max_intron=5000,
                                              **BIG, **cfg_kw),
                               stats)
        qlen = sum(ln for op, ln in read.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(read.query_sequence), read.cigartuples
        return changed, stats

    # -- Case B1: the intron arrived as a deletion --------------------------

    def _dop_read(self, g=GENOME_SEQ):
        query = g[P_DON - 60:P_DON] + g[P_ACC:P_ACC + 60]
        return _read('atac_dop', query,
                     [(0, 60), (2, self.INTRON), (0, 60)], P_DON - 60)

    def test_dop_stays_a_deletion_when_atac_is_off(self):
        # the pre-2026-09-05 default arm, now spelled explicitly
        r = self._dop_read()
        changed, stats = self._run(GENOME_SEQ, r, atac=False)
        assert stats.extra.get('arb_dop_checked') == 1, (
            'the fixture never reached Case B1')
        assert stats.extra.get('arb_dop_spliced', 0) == 0
        assert not changed
        assert r.cigartuples == [(0, 60), (2, self.INTRON), (0, 60)]

    def test_dop_becomes_an_intron_under_atac(self):
        r = self._dop_read()
        changed, stats = self._run(GENOME_SEQ, r, atac=True)
        assert changed and stats.extra.get('arb_dop_spliced') == 1, stats.extra
        assert r.cigartuples == [(0, 60), (3, self.INTRON), (0, 60)]
        assert self._junction(r) == (P_DON, P_ACC)
        assert r.get_tag('XB').startswith(f'dop:{P_DON}-{P_ACC}>')

    def test_dop_minus_strand(self):
        span = M_DON - M_ACC
        query = GENOME_SEQ[M_ACC - 60:M_ACC] + GENOME_SEQ[M_DON:M_DON + 60]
        r = _read('atac_dop_minus', query, [(0, 60), (2, span), (0, 60)],
                  M_ACC - 60, reverse=True)
        assert not self._run(GENOME_SEQ, r, atac=False)[0]
        r = _read('atac_dop_minus', query, [(0, 60), (2, span), (0, 60)],
                  M_ACC - 60, reverse=True)
        changed, stats = self._run(GENOME_SEQ, r, atac=True)
        assert changed and stats.extra.get('arb_dop_spliced') == 1, stats.extra
        assert self._junction(r) == (M_ACC, M_DON)

    def test_dop_pairing_is_enforced(self):
        # AT donor + AG acceptor: the union that would admit it is exactly what
        # planning/722 refuted, so even atac=True must leave the D alone.
        seq = list(GENOME_SEQ)
        seq[P_ACC - 2:P_ACC] = 'AG'
        g = ''.join(seq)
        r = self._dop_read(g)
        changed, stats = self._run(g, r, atac=True)
        assert stats.extra.get('arb_dop_checked') == 1
        if changed:
            assert self._junction(r) != (P_DON, P_ACC)
        else:
            assert r.cigartuples == [(0, 60), (2, self.INTRON), (0, 60)]

    # -- Case B2 / B3: the intron was aligned through linearly ---------------
    #
    # NOTE the default arm here does not simply refuse: on this dense 4.2 kb
    # synthetic contig it splices to a CHANCE GT..AG pair a few bp away
    # (measured: (1189, 1484) for B2, (1184, 1484) for B3, both at a worse edit
    # distance than the truth). That is the planning/721 residual in miniature —
    # a real non-canonical junction absorbed by a canonical decoy — so the
    # assertion is "only atac=True reaches the TRUE junction", not "the default
    # does nothing".

    def test_linear_read_reaches_the_true_junction_only_under_atac(self):
        query = GENOME_SEQ[P_DON - 200:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 160]
        r = _read('atac_mm', query, [(0, 360)], P_DON - 200)
        _, stats = self._run(GENOME_SEQ, r, atac=False)
        assert stats.extra.get('arb_mm_flagged', 0) >= 1, (
            'the fixture never reached the mismatch-flagged rescue')
        assert self._junction(r) != (P_DON, P_ACC)

        r = _read('atac_mm', query, [(0, 360)], P_DON - 200)
        changed, stats = self._run(GENOME_SEQ, r, atac=True)
        assert changed and stats.extra.get('arb_mm_spliced') == 1, stats.extra
        assert r.cigartuples == [(0, 200), (3, self.INTRON), (0, 160)]
        assert self._junction(r) == (P_DON, P_ACC)
        assert r.get_tag('XB').startswith(f'mm:{P_DON}-{P_ACC}:')

    def test_head_storm_left_mirror_only_under_atac(self):
        query = GENOME_SEQ[P_DON - 160:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 200]
        r = _read('atac_mmL', query, [(0, 360)], P_ACC - 160)
        _, stats = self._run(GENOME_SEQ, r, atac=False)
        assert stats.extra.get('arb_mm_flagged', 0) >= 1
        assert self._junction(r) != (P_DON, P_ACC)

        r = _read('atac_mmL', query, [(0, 360)], P_ACC - 160)
        changed, stats = self._run(GENOME_SEQ, r, atac=True)
        assert changed and stats.extra.get('arb_mm_spliced') == 1, stats.extra
        assert r.cigartuples == [(0, 160), (3, self.INTRON), (0, 200)]
        assert self._junction(r) == (P_DON, P_ACC)
        assert r.reference_start == P_DON - 160
        assert r.get_tag('XB').startswith(f'mmL:{P_DON}-{P_ACC}:')


class TestGrammar:
    def test_donor_rank_orders_gt_gc_atac(self):
        g = GENOME_SEQ
        assert _donor_rank(g, 'L', '+', P_DON, P_ACC) == 2
        assert _donor_rank(g, 'R', '-', M_ACC, M_DON) == 2
        s = list(g); s[P_DON:P_DON + 2] = 'GT'; s[P_ACC - 2:P_ACC] = 'AG'
        assert _donor_rank(''.join(s), 'L', '+', P_DON, P_ACC) == 0
        s[P_DON:P_DON + 2] = 'GC'
        assert _donor_rank(''.join(s), 'L', '+', P_DON, P_ACC) == 1

    def test_canonical_in_class_accepts_atac_by_default_and_can_be_disabled(self):
        # 2026-09-05 (Kevin): AT-AC is a canonical class everywhere by default;
        # atac=False is the opt-out (kept for A/B runs).
        g = GENOME_SEQ
        for j in ((P_DON, P_ACC), (M_ACC, M_DON)):
            assert is_canonical_junction(g, *j)
            assert is_canonical_junction(g, *j, atac=True)
            assert not is_canonical_junction(g, *j, atac=False)
            assert canonical_in_class(g, *j)
            assert not canonical_in_class(g, *j, atac=False)
        # a GT..AG junction is canonical either way
        s = list(g); s[P_DON:P_DON + 2] = 'GT'; s[P_ACC - 2:P_ACC] = 'AG'
        assert is_canonical_junction(''.join(s), P_DON, P_ACC)
        assert is_canonical_junction(''.join(s), P_DON, P_ACC, atac=True)


class TestDefaultsAndPlumbing:
    def test_default_is_on(self):
        """2026-09-05: Kevin flipped this default to True for every organism.

        AT-AC is not an exotic class some genomes carry — yeast splices it
        through the MAJOR spliceosome and human AT-AC is the U12-type class —
        so defaulting it off meant routinely snapping a real AT-AC junction
        onto a chance GT..AG (see TestArbiterDiscoveryPaths, where the off arm
        lands on a decoy 11-16 bp away at a WORSE edit distance)."""
        assert ResolverConfig().atac is True

    def test_run_overhang_resolver_defaults_the_knob_on(self):
        import inspect
        from rectify.core.align.overhang_resolver import run_overhang_resolver
        p = inspect.signature(run_overhang_resolver).parameters
        assert 'atac' in p
        assert p['atac'].default is True

    def test_disabling_reproduces_the_pre_flip_candidate_space(self, index):
        """The escape hatch has to actually work: atac=False must refuse the
        planted AT-AC junction exactly as the pre-flip default did."""
        r = _plus_read()
        changed = resolve_read(r, GENOME, index,
                               ResolverConfig(atac=False, **BIG), ResolverStats())
        if changed:
            assert not r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')

    def _align_parser(self, tmp_path):
        import argparse
        from rectify.core.commands.align_command import create_align_parser
        root = argparse.ArgumentParser()
        sub = root.add_subparsers(dest='command')
        create_align_parser(sub)
        return root, ['align', 'reads.fastq', '--genome', 'g.fa',
                      '-o', str(tmp_path), '--aligners', 'minimap2']

    def test_cli_defaults_on_align(self, tmp_path):
        root, base = self._align_parser(tmp_path)
        assert root.parse_args(base).resolver_atac is True

    def test_cli_no_flag_disables_on_align(self, tmp_path):
        root, base = self._align_parser(tmp_path)
        assert root.parse_args(base + ['--no-resolver-atac']).resolver_atac is False

    def test_cli_old_flag_is_a_compat_noop_on_align(self, tmp_path):
        """--resolver-atac was the opt-in switch; scripts written before the
        flip must keep working, now as a no-op."""
        root, base = self._align_parser(tmp_path)
        assert root.parse_args(base + ['--resolver-atac']).resolver_atac is True

    def _run_parser(self, tmp_path):
        import argparse
        from rectify.core.commands.run_command import create_run_parser
        root = argparse.ArgumentParser()
        sub = root.add_subparsers(dest='command')
        create_run_parser(sub)
        return root, ['run-all', 'reads.fastq', '--genome', 'g.fa',
                      '-o', str(tmp_path)]

    def test_cli_defaults_on_run_all(self, tmp_path):
        root, base = self._run_parser(tmp_path)
        assert root.parse_args(base).resolver_atac is True

    def test_cli_no_flag_disables_on_run_all(self, tmp_path):
        root, base = self._run_parser(tmp_path)
        assert root.parse_args(base + ['--no-resolver-atac']).resolver_atac is False

    def test_run_all_plumbing_defaults_on(self):
        """run-all reaches the resolver through run/stages.py, so its default
        has to agree with ResolverConfig's or the flag silently reverts."""
        import inspect
        from rectify.core.commands.run.stages import _run_alignment
        p = inspect.signature(_run_alignment).parameters
        assert p['resolver_atac'].default is True
