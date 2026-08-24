"""Tests for the information-bounded overhang resolver (planning/641).

Synthetic-genome geometry, all four clip-side x strand cases, the T1
zero-candidates-for-poly(A) assertion (on the COUNTER, per the SPEC), the
ambiguity rejection, the end-to-end BAM driver, and the env-gated
rescue_3ss_truncation wiring (dark by default; kills the poly(A)-clip
false-positive rescue when enabled).
"""

import random

import pysam
import pytest

from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    resolve_clip,
    resolve_read,
    run_overhang_resolver,
)
from rectify.core.splice.overhang_informativeness import COUNTERS, reset_counters
from rectify.core.splice.splice_site_index import SpliceSiteIndex

# ---------------------------------------------------------------------------
# Synthetic genome: two planted introns + an ambiguous decoy pair.
#
#   plus intron   [1200, 1500): GT ............................ AG
#   minus intron  [2200, 2500): CT ............................ AC
#   decoy:        acceptor end 3500; identical 30-mer X before GT donors at
#                 3000 and 3200 -> two ED-0 placements, different junctions
# ---------------------------------------------------------------------------

GLEN = 4200
P_DON, P_ACC = 1200, 1500      # plus intron [1200, 1500)
M_ACC, M_DON = 2200, 2500      # minus intron [2200, 2500)
A_ACC = 3500                   # ambiguous-case acceptor end
A_DON1, A_DON2 = 3000, 3200


def _make_genome():
    rng = random.Random(4242)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]

    def put(pos, s):
        seq[pos:pos + len(s)] = list(s)

    put(P_DON, 'GT')
    put(P_ACC - 2, 'AG')
    put(M_ACC, 'CT')
    put(M_DON - 2, 'AC')

    x = ''.join(rng.choice('ACGT') for _ in range(30))
    put(A_ACC - 2, 'AG')
    put(A_DON1 - 30, x)
    put(A_DON1, 'GT')
    put(A_DON2 - 30, x)
    put(A_DON2, 'GT')
    return ''.join(seq), x


GENOME_SEQ, AMBIG_X = _make_genome()
GENOME = {'chrI': GENOME_SEQ}


@pytest.fixture(scope='module')
def index():
    return SpliceSiteIndex.build(GENOME)


@pytest.fixture()
def cfg():
    return ResolverConfig(alpha=0.01, max_intron=5000)


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


def _read(name, query, cigar, ref_start, reverse=False, header=None):
    r = pysam.AlignedSegment(header or _header())
    r.query_name = name
    r.query_sequence = query
    r.flag = 16 if reverse else 0
    r.reference_id = 0
    r.reference_start = ref_start
    r.mapping_quality = 60
    r.cigartuples = cigar
    return r


class TestT1PolyARefusal:
    """641 T1: zero candidates evaluated for the poly(A) overhang —
    asserted on the counter, not the runtime."""

    def test_polya_clip_refused_zero_candidates(self, index, cfg):
        stats = ResolverStats()
        placement = resolve_clip(
            GENOME_SEQ, index, 'chrI', 'L', '+', 'A' * 30,
            edge=P_ACC, cfg=cfg, stats=stats,
        )
        assert placement is None
        assert stats.refused_low_info == 1
        assert stats.candidates_evaluated == 0
        assert COUNTERS['candidates_evaluated'] == 0

    def test_polya_read_emitted_unchanged(self, index, cfg):
        query = 'A' * 30 + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('polyA', query, [(4, 30), (0, 60)], P_ACC)
        stats = ResolverStats()
        changed = resolve_read(r, GENOME, index, cfg, stats)
        assert not changed
        assert r.cigartuples == [(4, 30), (0, 60)]
        assert r.reference_start == P_ACC
        assert stats.refused_low_info == 1
        assert COUNTERS['candidates_evaluated'] == 0

    def test_repeat_expansion_clip_refused(self, index, cfg):
        stats = ResolverStats()
        placement = resolve_clip(
            GENOME_SEQ, index, 'chrI', 'L', '+', 'AAG' * 12,
            edge=P_ACC, cfg=cfg, stats=stats,
        )
        assert placement is None
        assert stats.refused_repeat == 1
        assert stats.candidates_evaluated == 0


class TestResolution:
    def test_left_clip_plus_strand(self, index, cfg):
        # clip = exon-1 tail; aligned block starts at the intron end
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('left_plus', query, [(4, 30), (0, 60)], P_ACC)
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(0, 30), (3, P_ACC - P_DON), (0, 60)]
        assert r.reference_start == P_DON - 30
        assert r.get_tag('XJ').startswith(f'{P_DON}-{P_ACC}:')
        assert stats.resolved_left == 1

    def test_right_clip_plus_strand(self, index, cfg):
        # aligned over exon 1, right clip = exon-2 head
        query = GENOME_SEQ[P_DON - 60:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 30]
        r = _read('right_plus', query, [(0, 60), (4, 30)], P_DON - 60)
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(0, 60), (3, P_ACC - P_DON), (0, 30)]
        assert r.reference_start == P_DON - 60
        assert stats.resolved_right == 1

    def test_left_clip_minus_strand(self, index, cfg):
        # minus-strand transcript: CT..AC intron on the forward genome
        clip = GENOME_SEQ[M_ACC - 30:M_ACC]
        query = clip + GENOME_SEQ[M_DON:M_DON + 60]
        r = _read('left_minus', query, [(4, 30), (0, 60)], M_DON, reverse=True)
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(0, 30), (3, M_DON - M_ACC), (0, 60)]
        assert r.reference_start == M_ACC - 30
        assert stats.resolved_left == 1

    def test_right_clip_minus_strand(self, index, cfg):
        query = GENOME_SEQ[M_ACC - 60:M_ACC] + GENOME_SEQ[M_DON:M_DON + 30]
        r = _read('right_minus', query, [(0, 60), (4, 30)], M_ACC - 60, reverse=True)
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(0, 60), (3, M_DON - M_ACC), (0, 30)]
        assert stats.resolved_right == 1

    def test_long_clip_remainder_stays_softclipped(self, index, cfg):
        cfg.max_clip_match = 30
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = ('C' * 10 + clip) + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('long_clip', query, [(4, 40), (0, 60)], P_ACC)
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(4, 10), (0, 30), (3, P_ACC - P_DON), (0, 60)]

    def test_md_dropped_and_nm_recomputed_on_rewrite(self, index, cfg):
        """MD stays dropped (stale after CIGAR surgery, and not recomputed);
        NM is now RE-EMITTED against the rewritten placement.

        planning/769 defect 1 / planning/771 A5: the resolver used to drop NM
        and never put one back, so a 70%-mismatched block was invisible to
        every consumer that trusts NM — the single largest reason the P02
        artifact reached consensus unnoticed. This test previously asserted the
        old (defective) behaviour."""
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('tags', query, [(4, 30), (0, 60)], P_ACC)
        r.set_tag('NM', 0)
        r.set_tag('MD', '60')
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert not r.has_tag('MD')
        assert r.has_tag('NM')
        # the query was built verbatim from the reference at both placements
        assert r.get_tag('NM') == 0

    def test_nm_reports_a_mismatched_rescue(self, index, cfg):
        """A rescued block that does NOT match the reference must show it in NM."""
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        # corrupt 4 bases of the clip: they still place (hp-ED 4 <= 0.2*30=6)
        bad = list(clip)
        for i in (3, 9, 15, 21):
            bad[i] = 'A' if bad[i] != 'A' else 'C'
        query = ''.join(bad) + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('tags_mm', query, [(4, 30), (0, 60)], P_ACC)
        stats = ResolverStats()
        if resolve_read(r, GENOME, index, cfg, stats):
            assert r.has_tag('NM')
            assert r.get_tag('NM') >= 1

    def test_xj_appends_rather_than_overwrites(self, index, cfg):
        """Two resolved clips on one read must both appear in XJ.

        The previous unconditional ``set_tag('XJ', ...)`` silently dropped the
        first entry; every consumer already splits XJ on ','."""
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('xj', query, [(4, 30), (0, 60)], P_ACC)
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.has_tag('XJ')
        # XW, not XQ: XQ:i is the cDNA 5' pre-trim strip length (cdna/io.py:248,
        # READ_INTRINSIC in cma_schema.py). Appending to it clobbers a value
        # correct-cdna depends on, AND changes its type from i to Z.
        assert r.has_tag('XW')
        assert 'arm=' in str(r.get_tag('XW'))
        # XJ's three-field grammar is load-bearing for existing consumers
        for fld in str(r.get_tag('XJ')).split(','):
            span, ed, side = fld.split(':')
            assert '-' in span and float(ed) >= 0 and side in ('L', 'R', 'A1')


    def test_resolver_does_not_clobber_the_cdna_XQ_tag(self, index, cfg):
        """XQ:i belongs to correct-cdna (5' pre-trim strip length) and must survive.

        The resolver's own metadata rode on XQ until 2026-08-24, when a real
        production record read `XQ:Z:140,arm=arm2;...` — every read of the 748
        fixture carries an XQ:i, so the append silently overwrote it and changed
        its type. No test caught it because the fixtures here carry no XQ.
        """
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('xq_guard', query, [(4, 30), (0, 60)], P_ACC)
        r.set_tag('XQ', 140, 'i')
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert r.get_tag('XQ') == 140      # untouched, still an int
        assert r.has_tag('XW')             # resolver metadata went elsewhere


class TestRejection:
    def test_ambiguous_two_junctions_rejected(self, index, cfg):
        # identical exon tail X before donors at A_DON1 and A_DON2:
        # two ED-0 placements that are NOT the same junction -> reject
        query = AMBIG_X + GENOME_SEQ[A_ACC:A_ACC + 60]
        r = _read('ambig', query, [(4, 30), (0, 60)], A_ACC)
        stats = ResolverStats()
        assert not resolve_read(r, GENOME, index, cfg, stats)
        assert stats.rejected_ambiguous == 1
        assert r.cigartuples == [(4, 30), (0, 60)]

    def test_garbage_clip_no_accept(self, index, cfg):
        # complex (not refused) clip matching nothing near the junction
        rng = random.Random(555)
        clip = ''.join(rng.choice('ACGT') for _ in range(30))
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('garbage', query, [(4, 30), (0, 60)], P_ACC)
        stats = ResolverStats()
        assert not resolve_read(r, GENOME, index, cfg, stats)
        assert r.cigartuples == [(4, 30), (0, 60)]
        # it was assessed and searched, then rejected on edit distance
        assert stats.refused_low_info == 0

    def test_short_clip_ignored(self, index, cfg):
        query = GENOME_SEQ[P_DON - 4:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('short', query, [(4, 4), (0, 60)], P_ACC)
        stats = ResolverStats()
        assert not resolve_read(r, GENOME, index, cfg, stats)
        assert stats.clips_seen == 0

    def test_eq_encoded_clip_hard_fails(self, index, cfg):
        stats = ResolverStats()
        with pytest.raises(ValueError, match='corrupt'):
            resolve_clip(
                GENOME_SEQ, index, 'chrI', 'L', '+', '=' * 10 + 'ACGTACGTAC',
                edge=P_ACC, cfg=cfg, stats=stats,
            )


class TestDriver:
    def test_end_to_end_bam(self, tmp_path, cfg):
        fasta = tmp_path / 'genome.fa'
        with open(fasta, 'w') as fh:
            fh.write('>chrI\n' + GENOME_SEQ + '\n')

        header = _header()
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        reads = [
            _read('r1_resolve', clip + GENOME_SEQ[P_ACC:P_ACC + 60],
                  [(4, 30), (0, 60)], P_ACC, header=header),
            _read('r2_polya', 'A' * 30 + GENOME_SEQ[P_ACC:P_ACC + 60],
                  [(4, 30), (0, 60)], P_ACC, header=header),
            _read('r3_clean', GENOME_SEQ[100:190], [(0, 90)], 100, header=header),
        ]
        in_bam = tmp_path / 'in.bam'
        with pysam.AlignmentFile(in_bam, 'wb', header=header) as fh:
            for r in reads:
                fh.write(r)

        out_bam = tmp_path / 'out.bam'
        run_overhang_resolver(str(in_bam), str(fasta), str(out_bam))
        stats = run_overhang_resolver.last_stats
        assert stats.reads == 3
        assert stats.resolved == 1
        assert stats.refused_low_info == 1

        with pysam.AlignmentFile(out_bam, 'rb', check_sq=False) as fh:
            got = list(fh.fetch(until_eof=True))
            pg_ids = [pg.get('ID') for pg in fh.header.to_dict().get('PG', [])]
        assert [r.query_name for r in got] == ['r1_resolve', 'r2_polya', 'r3_clean']
        assert got[0].cigartuples == [(0, 30), (3, P_ACC - P_DON), (0, 60)]
        assert got[1].cigartuples == [(4, 30), (0, 60)]
        assert got[2].cigartuples == [(0, 90)]
        assert 'rectify-overhang-resolver' in pg_ids


class TestInsideEdge:
    """planning/644 T4c geometry: minimap2 mis-anchors the aligned start past
    the true acceptor when the intron tail duplicates the exon-1 end. The
    true near site then sits INSIDE the aligned block; the inside-edge
    extension must let it compete (and win) with the aligned bases
    re-assigned across the junction."""

    D, E = 1200, 1500

    @classmethod
    def _genome(cls):
        rng = random.Random(777)
        seq = [rng.choice('ACGT') for _ in range(4000)]
        D, E = cls.D, cls.E

        def put(pos, s):
            seq[pos:pos + len(s)] = list(s)

        put(D, 'GT')
        put(E - 2, 'AG')
        # exon-1 end whose last 2 nt are AG, duplicated as the intron tail —
        # the homology that makes minimap2 anchor 10 bp into the intron.
        put(D - 10, ''.join(rng.choice('ACGT') for _ in range(8)) + 'AG')
        seq[E - 10:E] = seq[D - 10:D]
        return ''.join(seq)

    def _read(self, g, inside_encoded=False):
        D, E = self.D, self.E
        clip20 = g[D - 30:D - 10]
        inside = '=' * 10 if inside_encoded else g[E - 10:E]
        query = clip20 + inside + g[E:E + 50]
        header = pysam.AlignmentHeader.from_dict({
            'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': len(g)}]})
        r = pysam.AlignedSegment(header)
        r.query_name = 'misanchor'
        r.query_sequence = query
        r.flag = 0
        r.reference_id = 0
        r.reference_start = E - 10
        r.mapping_quality = 60
        r.cigartuples = [(4, 20), (0, 60)]
        return r

    def _run(self, inside_encoded):
        from rectify.core.splice.overhang_informativeness import same_junction
        g = self._genome()
        genome = {'chrI': g}
        index = SpliceSiteIndex.build(genome)
        cfg = ResolverConfig(alpha=0.01, max_intron=5000)
        r = self._read(g, inside_encoded=inside_encoded)
        stats = ResolverStats()
        assert resolve_read(r, genome, index, cfg, stats), \
            f"not resolved: {stats.as_dict()}"
        assert stats.extra.get('resolved_inside_edge', 0) == 1
        # junction identity up to the (deliberate) 10-bp ambiguity class
        m = re_len = None
        ref = r.reference_start
        for op, ln in r.cigartuples:
            if op == 3:
                m, re_len = ref, ref + ln
                break
            if op in (0, 2, 7, 8):
                ref += ln
        assert m is not None, r.cigartuples
        assert same_junction(g, (m, re_len), (self.D, self.E))
        # query length must be conserved by the rewrite
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(r.query_sequence)
        return r, stats

    def test_misanchored_acceptor_recovered(self):
        r, stats = self._run(inside_encoded=False)
        assert r.get_tag('XJ')

    def test_eq_encoded_inside_bases_decoded(self):
        # calmd '='-encoded aligned bases must be decoded from the reference
        # before entering the comparison (the clip itself still hard-fails).
        self._run(inside_encoded=True)


class TestRearbitration:
    """planning/644b — the peelback move generalized: junction ASSIGNMENTS
    and suspicious linear structure are hypotheses, re-scored against
    index-derived competitors under the information budget."""

    D, E_TRUE = 1200, 1500     # true intron [D, E_TRUE), GT..AG
    E_DEC = 1488               # decoy AG 12 bp inside the intron tail

    @classmethod
    def _genome(cls):
        rng = random.Random(20260809)
        seq = [rng.choice('ACGT') for _ in range(4000)]

        def put(pos, s):
            seq[pos:pos + len(s)] = list(s)

        put(cls.D, 'GT')
        put(cls.E_TRUE - 2, 'AG')
        put(cls.E_DEC - 2, 'AG')       # the chance acceptor
        return ''.join(seq)

    @staticmethod
    def _mk_read(g, query, cigar, rs):
        header = pysam.AlignmentHeader.from_dict({
            'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': len(g)}]})
        r = pysam.AlignedSegment(header)
        r.query_name = 'arb'
        r.query_sequence = query
        r.flag = 0
        r.reference_id = 0
        r.reference_start = rs
        r.mapping_quality = 60
        r.cigartuples = cigar
        return r

    def _run(self, g, read, **cfg_kw):
        genome = {'chrI': g}
        index = SpliceSiteIndex.build(genome)
        cfg = ResolverConfig(alpha=0.01, max_intron=5000, **cfg_kw)
        stats = ResolverStats()
        changed = resolve_read(read, genome, index, cfg, stats)
        return changed, stats

    def _junction(self, read):
        ref = read.reference_start
        for op, ln in read.cigartuples:
            if op == 3:
                return ref, ref + ln
            if op in (0, 2, 7, 8):
                ref += ln
        return None

    def test_acceptor_shift_recovers_true_junction(self):
        # read spliced at (D, E_TRUE) but CIGAR asserts the decoy acceptor:
        # M_right's query (true exon-2) is misaligned over the intron tail
        g = self._genome()
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        cigar = [(0, 60), (3, self.E_DEC - self.D), (0, 60)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(g, r)
        assert changed and stats.extra.get('arb_shifted') == 1
        assert self._junction(r) == (self.D, self.E_TRUE)
        assert r.get_tag('XB').startswith('shift:')
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_donor_shift_recovers_true_junction(self):
        # decoy donor 8 bp downstream of the true one (GT planted at both);
        # the claimed CIGAR over-extends exon-1 by 8 query bases
        g = list(self._genome())
        g[self.D + 8:self.D + 10] = list('GT')
        g = ''.join(g)
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        cigar = [(0, 68), (3, self.E_TRUE - (self.D + 8)), (0, 52)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(g, r)
        assert changed and stats.extra.get('arb_shifted', 0) >= 1
        assert self._junction(r) == (self.D, self.E_TRUE)
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_correct_junction_not_shifted(self):
        # control: the true assignment must survive re-arbitration untouched
        g = self._genome()
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        cigar = [(0, 60), (3, self.E_TRUE - self.D), (0, 60)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(g, r)
        assert stats.extra.get('arb_shifted', 0) == 0
        assert self._junction(r) == (self.D, self.E_TRUE)

    def test_pseudo_slide_grammar_tiebreak(self):
        # SRC1-class dispute: the claimed junction is a pseudo-slide of the
        # truth (+4 on BOTH boundaries, beyond the legal window because the
        # 4-bp homology is partial: 1 mismatch). The DP near-ties (margin
        # can't move it) but only the true placement is canonical-class —
        # splicing grammar must adjudicate.
        g = list(self._genome())
        D, E = self.D, self.E_TRUE
        # partial homology g[D+i]==g[E+i] for i in {0,1,3}, differing at i=2
        g[E], g[E + 1] = g[D], g[D + 1]          # == 'G','T'
        g[E + 3] = g[D + 3] = 'A'
        g[D + 2] = 'C'
        g[E + 2] = 'G'
        g[D + 4] = 'A'                            # claimed donor dinuc 'A?' -> non-canonical
        g[E + 4] = 'T' if g[D + 4] != 'T' else 'C'
        g = ''.join(g)
        query = g[D - 60:D] + g[E:E + 60]
        cigar = [(0, 64), (3, E - D), (0, 56)]    # claimed = (D+4, E+4)
        r = self._mk_read(g, query, cigar, D - 60)
        changed, stats = self._run(g, r)
        assert changed and stats.extra.get('arb_shifted') == 1, stats.extra
        assert self._junction(r) == (D, E)
        # a grammar-admitted move is marked ':g' and counted, so the triage
        # layer can treat it as TRIAGED rather than high-confidence
        assert r.get_tag('XB').endswith(':g')
        assert stats.extra.get('arb_grammar_tiebreak') == 1
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_grammar_off_preserves_noncanonical_claim(self):
        # arb_grammar=False: the same SRC1-class dispute must NOT be snapped
        # to the canonical placement — the DP near-ties and margin alone
        # cannot move it. This is the discovery-mission setting (realigner
        # noncanon control: grammar snap flattens genuine non-canonical
        # junctions when a canonical decoy is in range).
        g = list(self._genome())
        D, E = self.D, self.E_TRUE
        g[E], g[E + 1] = g[D], g[D + 1]
        g[E + 3] = g[D + 3] = 'A'
        g[D + 2] = 'C'
        g[E + 2] = 'G'
        g[D + 4] = 'A'
        g[E + 4] = 'T' if g[D + 4] != 'T' else 'C'
        g = ''.join(g)
        query = g[D - 60:D] + g[E:E + 60]
        cigar = [(0, 64), (3, E - D), (0, 56)]    # claimed = (D+4, E+4)
        r = self._mk_read(g, query, cigar, D - 60)
        changed, stats = self._run(g, r, arb_grammar=False)
        assert stats.extra.get('arb_shifted', 0) == 0, stats.extra
        assert stats.extra.get('arb_grammar_tiebreak', 0) == 0
        assert self._junction(r) == (D + 4, E + 4)

    def test_ref_consuming_tail_refuses_shift_conserves_query(self):
        # v5.1 regression (chrII 98b21cfd, 2026-08-10): junction followed by
        # 1I 4M 4D 54M 1S — a shift here changes net reference consumption
        # and (pre-fix) the delta<0 M-swap overwrote the flanking 1I with a
        # manufactured M, losing a query base (htslib-unreadable record).
        # Post-fix: the junction enters arbitration but is frame-unsafe —
        # no shift, CIGAR untouched, query length conserved.
        g = self._genome()
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        nd = self.E_DEC - self.D
        cigar = [(0, 60), (3, nd), (1, 1), (0, 4), (2, 4), (0, 54), (4, 1)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(g, r)
        assert stats.extra.get('arb_shifted', 0) == 0, stats.extra
        assert stats.extra.get('arb_frame_unsafe_skip', 0) >= 1, stats.extra
        assert self._junction(r) == (self.D, self.E_DEC)   # held, not corrupted
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_interior_junction_shift_refused_downstream_anchored(self):
        # v5.1: shifting an interior junction displaces every ref-consuming
        # op downstream (a -delta interior left shift moves the NEXT
        # junction off its index sites). Interior junctions must be refused
        # shift arbitration; both junctions stay put.
        g = self._genome()
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        nd = self.E_DEC - self.D
        cigar = [(0, 60), (3, nd), (0, 20), (3, 200), (0, 40)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(g, r)
        assert stats.extra.get('arb_shifted', 0) == 0, stats.extra
        assert stats.extra.get('arb_frame_unsafe_skip', 0) >= 1, stats.extra
        juncs = []
        ref = r.reference_start
        for op, ln in r.cigartuples:
            if op == 3:
                juncs.append((ref, ref + ln))
            if op in (0, 2, 3, 7, 8):
                ref += ln
        assert juncs == [(self.D, self.E_DEC),
                         (self.E_DEC + 20, self.E_DEC + 220)]
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_skip_region_bypasses_all_rescue(self):
        # rDNA-filter contract (planning/644b perf lever): a read overlapping
        # a skip region is written through UNTOUCHED — the acceptor-shift
        # scenario that normally fires arb_shifted must not fire, and only
        # the skipped_region counter moves.
        g = self._genome()
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        cigar = [(0, 60), (3, self.E_DEC - self.D), (0, 60)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(
            g, r, skip_regions={'chrI': [(0, 4000)]})
        assert not changed
        assert stats.extra.get('skipped_region') == 1
        assert stats.extra.get('arb_shifted', 0) == 0
        assert stats.clips_seen == 0
        assert self._junction(r) == (self.D, self.E_DEC)

    def test_skip_region_elsewhere_does_not_bypass(self):
        # same read, region on a coordinate range the read never touches —
        # rescue proceeds normally (the shift fires as in the base test)
        g = self._genome()
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        cigar = [(0, 60), (3, self.E_DEC - self.D), (0, 60)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(
            g, r, skip_regions={'chrI': [(3900, 4000)]})
        assert changed and stats.extra.get('arb_shifted') == 1
        assert stats.extra.get('skipped_region', 0) == 0
        assert self._junction(r) == (self.D, self.E_TRUE)

    def test_boundary_deletion_merged_into_intron(self):
        # the SRC1 smoking gun (planning/644c): M D4 N — the aligner encodes
        # the GC-donor evidence as a 4-bp deletion while asserting the GT
        # donor 4 bp downstream; merging the D into the intron is
        # alignment-identical and lands on the canonical GC junction
        g = list(self._genome())
        D, E = self.D, self.E_TRUE
        g[D:D + 2] = list('GC')       # true donor (GC class)
        g[D + 4:D + 6] = list('GT')   # the decoy donor the aligner asserted
        g = ''.join(g)
        query = g[D - 60:D] + g[E:E + 60]
        cigar = [(0, 60), (2, 4), (3, E - (D + 4)), (0, 60)]
        r = self._mk_read(g, query, cigar, D - 60)
        changed, stats = self._run(g, r)
        assert changed and stats.extra.get('arb_dmerge') == 1, stats.extra
        assert self._junction(r) == (D, E)
        assert not any(op == 2 for op, _ in r.cigartuples)
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_boundary_deletion_not_merged_when_noncanonical(self):
        # control: absorbing the D must be refused when the merged junction's
        # ambiguity class has no canonical member
        g = list(self._genome())
        D, E = self.D, self.E_TRUE
        g[D:D + 2] = list('AA')       # merged donor would be non-canonical
        g[D + 2] = 'C'                # break any slide chain at the boundary
        g[D + 4:D + 6] = list('GT')
        g = ''.join(g)
        from rectify.core.splice.overhang_informativeness import canonical_in_class
        assert not canonical_in_class(g, D, E)
        query = g[D - 60:D] + g[E:E + 60]
        cigar = [(0, 60), (2, 4), (3, E - (D + 4)), (0, 60)]
        r = self._mk_read(g, query, cigar, D - 60)
        changed, stats = self._run(g, r)
        assert stats.extra.get('arb_dmerge', 0) == 0
        assert (2, 4) in r.cigartuples

    def test_dop_converted_to_snapped_intron(self):
        # the intron expressed as a 300-bp deletion at the exact bounds
        g = self._genome()
        query = g[self.D - 60:self.D] + g[self.E_TRUE:self.E_TRUE + 60]
        cigar = [(0, 60), (2, self.E_TRUE - self.D), (0, 60)]
        r = self._mk_read(g, query, cigar, self.D - 60)
        changed, stats = self._run(g, r)
        assert changed and stats.extra.get('arb_dop_spliced') == 1
        assert self._junction(r) == (self.D, self.E_TRUE)
        assert not any(op == 2 and ln >= 20 for op, ln in r.cigartuples)

    def test_linear_block_spliced_on_mismatch_cluster(self):
        # spliced molecule aligned LINEARLY through the intron: exon-2 query
        # misaligned over intron sequence => mismatch storm past the donor
        g = self._genome()
        query = g[self.D - 200:self.D] + g[self.E_TRUE:self.E_TRUE + 160]
        cigar = [(0, 360)]
        r = self._mk_read(g, query, cigar, self.D - 200)
        changed, stats = self._run(g, r)
        assert stats.extra.get('arb_mm_flagged', 0) >= 1
        assert changed and stats.extra.get('arb_mm_spliced') == 1
        assert self._junction(r) == (self.D, self.E_TRUE)
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_interior_dmerge_preserves_downstream_junction(self):
        # stage-1 vicinity triggers: a mis-assigned FIRST junction (M D4 N,
        # SRC1 pattern) in a two-intron read must be merged while the second
        # junction stays byte-identical (frame safety)
        g = list(self._genome())
        D, E = self.D, self.E_TRUE
        D2, E2 = 1700, 1900
        g[D:D + 2] = list('GC')
        g[D + 4:D + 6] = list('GT')
        g[D2:D2 + 2] = list('GT')
        g[E2 - 2:E2] = list('AG')
        g = ''.join(g)
        query = g[D - 60:D] + g[E:D2] + g[E2:E2 + 60]
        cigar = [(0, 60), (2, 4), (3, E - (D + 4)), (0, D2 - E), (3, E2 - D2), (0, 60)]
        r = self._mk_read(g, query, cigar, D - 60)
        changed, stats = self._run(g, r)
        assert changed and stats.extra.get('arb_dmerge') == 1
        juncs = []
        ref = r.reference_start
        for op, ln in r.cigartuples:
            if op == 3:
                juncs.append((ref, ref + ln))
            if op in (0, 2, 3, 7, 8):
                ref += ln
        assert juncs == [(D, E), (D2, E2)]
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_head_storm_spliced_left_mirror(self):
        # B3: linear alignment THROUGH the intron anchored on exon-2 — the
        # exon-1 query bases smear over the intron tail (storm at the block
        # HEAD); splicing re-anchors reference_start leftward
        g = self._genome()
        D, E = self.D, self.E_TRUE
        query = g[D - 160:D] + g[E:E + 200]
        cigar = [(0, 360)]
        r = self._mk_read(g, query, cigar, E - 160)
        changed, stats = self._run(g, r)
        assert changed and stats.extra.get('arb_mm_spliced') == 1, stats.extra
        assert self._junction(r) == (D, E)
        assert r.reference_start == D - 160
        assert r.get_tag('XB').startswith('mmL:')
        qlen = sum(ln for op, ln in r.cigartuples if op in (0, 1, 4, 7, 8))
        assert qlen == len(query)

    def test_linear_matching_block_untouched(self):
        # control: a genuinely linear (e.g. intron-retained) read must NOT be spliced
        g = self._genome()
        query = g[self.D - 200:self.D + 160]
        cigar = [(0, 360)]
        r = self._mk_read(g, query, cigar, self.D - 200)
        changed, stats = self._run(g, r)
        assert stats.extra.get('arb_mm_spliced', 0) == 0
        assert r.cigartuples == [(0, 360)]


# ---------------------------------------------------------------------------
# rescue_3ss_truncation gate wiring (dark by default)
# ---------------------------------------------------------------------------

class TestRescueGate:
    """The gate kills the poly(A)-clip false-positive rescue when enabled,
    changes nothing when the clip is informative, and is OFF by default."""

    # A-run exon end: genome[FP_DON-30:FP_DON] == 'A'*30, so a poly(A) clip
    # FALSELY sequence-matches the exon-1 side of this junction.
    FP_DON, FP_ACC = 700, 900

    @classmethod
    def _fp_genome(cls):
        seq = list(GENOME_SEQ)
        seq[cls.FP_DON - 30:cls.FP_DON] = list('A' * 30)
        seq[cls.FP_DON:cls.FP_DON + 2] = list('GT')
        seq[cls.FP_ACC - 2:cls.FP_ACC] = list('AG')
        return {'chrI': ''.join(seq)}

    def _polya_read(self, genome):
        query = 'A' * 30 + genome['chrI'][self.FP_ACC:self.FP_ACC + 60]
        return _read('fp', query, [(4, 30), (0, 60)], self.FP_ACC)

    def test_gate_off_polya_clip_false_rescues(self, monkeypatch):
        from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
        monkeypatch.delenv('RECTIFY_OVERHANG_INFO_GATE', raising=False)
        genome = self._fp_genome()
        res = rescue_3ss_truncation(
            self._polya_read(genome), genome,
            {('chrI', self.FP_DON, self.FP_ACC)},
        )
        # the defect the gate exists to fix: a poly(A) clip "confirms" the
        # A-run exon end
        assert res['rescued'] is True

    def test_gate_on_polya_clip_refused(self, monkeypatch):
        from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
        monkeypatch.setenv('RECTIFY_OVERHANG_INFO_GATE', '1')
        genome = self._fp_genome()
        res = rescue_3ss_truncation(
            self._polya_read(genome), genome,
            {('chrI', self.FP_DON, self.FP_ACC)},
        )
        assert res['rescued'] is False
        assert COUNTERS['refused'] >= 1

    def test_gate_on_informative_clip_identical_result(self, monkeypatch):
        from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
        genome = GENOME
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        junctions = {('chrI', P_DON, P_ACC)}

        monkeypatch.delenv('RECTIFY_OVERHANG_INFO_GATE', raising=False)
        r1 = _read('cmp', query, [(4, 30), (0, 60)], P_ACC)
        res_off = rescue_3ss_truncation(r1, genome, junctions)

        monkeypatch.setenv('RECTIFY_OVERHANG_INFO_GATE', '1')
        r2 = _read('cmp', query, [(4, 30), (0, 60)], P_ACC)
        res_on = rescue_3ss_truncation(r2, genome, junctions)

        assert res_off == res_on
        assert res_on['rescued'] is True
