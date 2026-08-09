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

    def test_md_nm_dropped_on_rewrite(self, index, cfg):
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = clip + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('tags', query, [(4, 30), (0, 60)], P_ACC)
        r.set_tag('NM', 0)
        r.set_tag('MD', '60')
        stats = ResolverStats()
        assert resolve_read(r, GENOME, index, cfg, stats)
        assert not r.has_tag('NM') and not r.has_tag('MD')


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
