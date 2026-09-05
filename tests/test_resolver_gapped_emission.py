"""The clip resolver must EMIT the indels its score assumed (item A7).

``resolve_clip`` scores a candidate placement with ``hp_edit_distance_bounded``,
which allows indels, but the emitter used to write the matched block as ONE flat
``M`` op — an alignment that cannot express the indels the score was built on.
Real example: read r060_2601 (P06) was written ``...405N24M`` for a 24-bp block
whose true fit needs 1-bp shifts on both flanks, i.e. 13/24 apparent identity,
and ``MD``/``NM`` are dropped by the same function so nothing downstream ever
noticed the disagreement.

These tests pin the fix and its boundaries:

* a gap-free block still emits EXACTLY the old flat ``M`` (all four
  side x strand cases — the CIGARs pinned in ``test_overhang_resolver.py`` must
  not move);
* a 1-bp deletion emits ``xM1DyM``, a 1-bp insertion ``xM1IyM``;
* query length is conserved and the emitted reference span agrees with
  ``reference_end - reference_start`` on every case — including LEFT-side
  rewrites, where ``reference_start`` must come from the block's ACTUAL
  reference consumption and not from ``placement.m``;
* the accepted junction never moves: the ``N`` op still lands exactly on
  ``(intron_start, intron_end)`` and keeps the accepted intron length;
* a degenerate alignment falls back to the flat ``M`` and is counted
  (``emit_fallback_flat``) rather than writing a junction-moving ``D``;
* the P06 shape: a block whose UNGAPPED identity is poor but whose gapped
  identity is high now emits the indel and reads >= 90% identical per base.
"""

import random

import pysam
import pytest

from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    _LEFT,
    _RIGHT,
    _block_cigar,
    resolve_read,
)
from rectify.core.splice.overhang_informativeness import reset_counters
from rectify.core.splice.splice_site_index import SpliceSiteIndex

# ---------------------------------------------------------------------------
# Synthetic genome (same shape as tests/test_resolver_prp18_acceptors.py):
#   plus  intron [1200, 1500): GT ........ AG
#   minus intron [2200, 2500): forward-genome CT ........ AC
# ---------------------------------------------------------------------------

GLEN = 4200
P_DON, P_ACC = 1200, 1500
M_ACC, M_DON = 2200, 2500


def _make_genome():
    rng = random.Random(70601)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]

    def put(pos, s):
        seq[pos:pos + len(s)] = list(s)

    put(P_DON, 'GT')
    put(P_ACC - 2, 'AG')
    put(M_ACC, 'CT')
    put(M_DON - 2, 'AC')
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


# ---------------------------------------------------------------------------
# Invariant helpers — every case runs all of them.
# ---------------------------------------------------------------------------

_Q_OPS = (0, 1, 4, 7, 8)     # consume query
_R_OPS = (0, 2, 3, 7, 8)     # consume reference


def assert_query_conserved(read, original_query):
    """Hard constraint (1): the rewrite may not change query length."""
    q = sum(ln for op, ln in read.cigartuples if op in _Q_OPS)
    assert q == len(original_query), (
        f'query length changed: CIGAR consumes {q}, SEQ is '
        f'{len(original_query)} ({read.cigartuples})')


def assert_reference_span_consistent(read):
    """Hard constraint (4): D consumes reference, I consumes query, so the
    emitted span must still agree with pysam's own end coordinate."""
    span = sum(ln for op, ln in read.cigartuples if op in _R_OPS)
    assert span == read.reference_end - read.reference_start, (
        f'reference span {span} != {read.reference_end} - '
        f'{read.reference_start} ({read.cigartuples})')


def junction_of(read):
    """(intron_start, intron_end) of the single N op, walked from
    reference_start — the check that catches a bad reference_start."""
    ref = read.reference_start
    for op, ln in read.cigartuples:
        if op == 3:
            return ref, ref + ln
        if op in _R_OPS:
            ref += ln
    return None


def block_ops(read, side):
    """The M/I/D ops the resolver placed on the far side of the junction."""
    ct = list(read.cigartuples)
    i = next(j for j, (op, _) in enumerate(ct) if op == 3)
    ops = ct[i + 1:] if side == _RIGHT else ct[:i]
    return [(op, ln) for op, ln in ops if op != 4]


def block_identity(read, side, genome_seq):
    """(matches, aligned_bases) over the placed block's M ops."""
    ct = list(read.cigartuples)
    i = next(j for j, (op, _) in enumerate(ct) if op == 3)
    lo, hi = (i + 1, len(ct)) if side == _RIGHT else (0, i)
    q = read.query_sequence
    qi, ref = 0, read.reference_start
    match = total = 0
    for j, (op, ln) in enumerate(ct):
        if op in (0, 7, 8):
            if lo <= j < hi:
                for d in range(ln):
                    total += 1
                    if q[qi + d].upper() == genome_seq[ref + d].upper():
                        match += 1
            qi += ln
            ref += ln
        elif op in (1, 4):
            qi += ln
        elif op in (2, 3):
            ref += ln
    return match, total


def _resolve(read, index, **cfg_kw):
    stats = ResolverStats()
    original = read.query_sequence
    changed = resolve_read(read, GENOME, index,
                           ResolverConfig(alpha=0.01, max_intron=5000, **cfg_kw),
                           stats)
    assert_query_conserved(read, original)
    assert_reference_span_consistent(read)
    return changed, stats


# ---------------------------------------------------------------------------
# Constraint (2): a gap-free block emits EXACTLY the old flat M.
# ---------------------------------------------------------------------------

class TestGapFreeIsByteIdentical:
    """These four CIGARs are the ones pinned in test_overhang_resolver.py
    L147/158/169/178. If one of them moves, the emitter is wrong — do not
    re-pin the other file."""

    def test_left_clip_plus(self, index):
        query = GENOME_SEQ[P_DON - 30:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('gapfree_Lp', query, [(4, 30), (0, 60)], P_ACC)
        assert _resolve(r, index)[0]
        assert r.cigartuples == [(0, 30), (3, P_ACC - P_DON), (0, 60)]
        assert r.reference_start == P_DON - 30

    def test_right_clip_plus(self, index):
        query = GENOME_SEQ[P_DON - 60:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 30]
        r = _read('gapfree_Rp', query, [(0, 60), (4, 30)], P_DON - 60)
        assert _resolve(r, index)[0]
        assert r.cigartuples == [(0, 60), (3, P_ACC - P_DON), (0, 30)]

    def test_left_clip_minus(self, index):
        query = GENOME_SEQ[M_ACC - 30:M_ACC] + GENOME_SEQ[M_DON:M_DON + 60]
        r = _read('gapfree_Lm', query, [(4, 30), (0, 60)], M_DON, reverse=True)
        assert _resolve(r, index)[0]
        assert r.cigartuples == [(0, 30), (3, M_DON - M_ACC), (0, 60)]
        assert r.reference_start == M_ACC - 30

    def test_right_clip_minus(self, index):
        query = GENOME_SEQ[M_ACC - 60:M_ACC] + GENOME_SEQ[M_DON:M_DON + 30]
        r = _read('gapfree_Rm', query, [(0, 60), (4, 30)], M_ACC - 60,
                  reverse=True)
        assert _resolve(r, index)[0]
        assert r.cigartuples == [(0, 60), (3, M_DON - M_ACC), (0, 30)]

    def test_long_clip_remainder_still_softclipped(self, index):
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        query = ('C' * 10 + clip) + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('remainder', query, [(4, 40), (0, 60)], P_ACC)
        assert _resolve(r, index, max_clip_match=30)[0]
        assert r.cigartuples == [(4, 10), (0, 30), (3, P_ACC - P_DON), (0, 60)]


# ---------------------------------------------------------------------------
# The fix itself: indels the score assumed are now spelled out.
# ---------------------------------------------------------------------------

class TestGappedEmission:
    def test_right_block_one_base_deletion(self, index):
        # exon-2 head missing its 15th base relative to the genome
        clip = GENOME_SEQ[P_ACC:P_ACC + 30]
        q = clip[:15] + clip[16:]
        query = GENOME_SEQ[P_DON - 60:P_DON] + q
        r = _read('del_R', query, [(0, 60), (4, len(q))], P_DON - 60)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert junction_of(r) == (P_DON, P_ACC)
        assert block_ops(r, _RIGHT) == [(0, 15), (2, 1), (0, 14)]
        assert stats.extra.get('emit_fallback_flat', 0) == 0
        m, t = block_identity(r, _RIGHT, GENOME_SEQ)
        assert m == t == 29

    def test_right_block_one_base_insertion(self, index):
        clip = GENOME_SEQ[P_ACC:P_ACC + 30]
        q = clip[:15] + 'T' + clip[15:]
        query = GENOME_SEQ[P_DON - 60:P_DON] + q
        r = _read('ins_R', query, [(0, 60), (4, len(q))], P_DON - 60)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert junction_of(r) == (P_DON, P_ACC)
        ops = block_ops(r, _RIGHT)
        assert [op for op, _ in ops] == [0, 1, 0]
        assert dict(  # one inserted base, 30 matched
            (op, ln) for op, ln in ops)[1] == 1
        assert sum(ln for op, ln in ops if op == 0) == 30
        assert stats.extra.get('emit_fallback_flat', 0) == 0

    def test_left_block_one_base_deletion(self, index):
        # exon-1 tail missing its 15th base; the block ENDS at intron_start,
        # so reference_start must absorb the extra reference base
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        q = clip[:15] + clip[16:]
        query = q + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('del_L', query, [(4, len(q)), (0, 60)], P_ACC)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert junction_of(r) == (P_DON, P_ACC)
        assert block_ops(r, _LEFT) == [(0, 15), (2, 1), (0, 14)]
        # 29 query bases spanning 30 reference bases: the old
        # `intron_start - placement.m` arithmetic would have been 1 bp off
        assert r.reference_start == P_DON - 30
        assert stats.extra.get('emit_fallback_flat', 0) == 0
        m, t = block_identity(r, _LEFT, GENOME_SEQ)
        assert m == t == 29

    def test_left_block_one_base_insertion(self, index):
        clip = GENOME_SEQ[P_DON - 30:P_DON]
        q = clip[:15] + 'T' + clip[15:]
        query = q + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('ins_L', query, [(4, len(q)), (0, 60)], P_ACC)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert junction_of(r) == (P_DON, P_ACC)
        ops = block_ops(r, _LEFT)
        assert [op for op, _ in ops] == [0, 1, 0]
        assert sum(ln for op, ln in ops if op == 1) == 1
        # 31 query bases over 30 reference bases
        assert r.reference_start == P_DON - 30

    def test_minus_strand_gapped_left(self, index):
        clip = GENOME_SEQ[M_ACC - 30:M_ACC]
        q = clip[:15] + clip[16:]
        query = q + GENOME_SEQ[M_DON:M_DON + 60]
        r = _read('del_Lm', query, [(4, len(q)), (0, 60)], M_DON, reverse=True)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert junction_of(r) == (M_ACC, M_DON)
        ops = block_ops(r, _LEFT)
        assert [op for op, _ in ops] == [0, 2, 0]
        assert r.reference_start == M_ACC - 30

    def test_remainder_softclip_composes_with_a_gapped_block_left(self, index):
        # clip longer than max_clip_match AND the used portion carries an
        # indel: the remainder S and the block ops have to compose, and the
        # block must still begin with an M (an `S I M` opening would be a
        # nonsense record). The gap-free version of this case takes the
        # zero-mismatch fast path, so only this one exercises the composition.
        tail = GENOME_SEQ[P_DON - 31:P_DON]           # 31 reference bases
        used = tail[:15] + tail[16:]                  # 30 query bases, 1 deleted
        query = 'C' * 10 + used + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('rem_gap_L', query, [(4, 40), (0, 60)], P_ACC)
        changed, stats = _resolve(r, index, max_clip_match=30)
        assert changed, stats.as_dict()
        assert r.cigartuples[0] == (4, 10)
        assert r.cigartuples[1][0] == 0, (
            f'block does not open with an M: {r.cigartuples}')
        ops = block_ops(r, _LEFT)
        assert [op for op, _ in ops] == [0, 2, 0]
        assert sum(ln for op, ln in ops if op == 0) == 30
        assert junction_of(r) == (P_DON, P_ACC)
        assert r.reference_start == P_DON - 31       # 30 query over 31 ref
        m, t = block_identity(r, _LEFT, GENOME_SEQ)
        assert m == t == 30

    def test_remainder_softclip_composes_with_a_gapped_block_right(self, index):
        head = GENOME_SEQ[P_ACC:P_ACC + 31]
        used = head[:15] + head[16:]                  # 30 query bases, 1 deleted
        query = GENOME_SEQ[P_DON - 60:P_DON] + used + 'C' * 10
        r = _read('rem_gap_R', query, [(0, 60), (4, 40)], P_DON - 60)
        changed, stats = _resolve(r, index, max_clip_match=30)
        assert changed, stats.as_dict()
        assert r.cigartuples[-1] == (4, 10)
        assert r.cigartuples[-2][0] == 0, (
            f'block does not close with an M: {r.cigartuples}')
        ops = block_ops(r, _RIGHT)
        assert [op for op, _ in ops] == [0, 2, 0]
        assert sum(ln for op, ln in ops if op == 0) == 30
        assert junction_of(r) == (P_DON, P_ACC)

    def test_minus_strand_gapped_right(self, index):
        clip = GENOME_SEQ[M_DON:M_DON + 30]
        q = clip[:15] + 'T' + clip[15:]
        query = GENOME_SEQ[M_ACC - 60:M_ACC] + q
        r = _read('ins_Rm', query, [(0, 60), (4, len(q))], M_ACC - 60,
                  reverse=True)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert junction_of(r) == (M_ACC, M_DON)
        assert [op for op, _ in block_ops(r, _RIGHT)] == [0, 1, 0]


# ---------------------------------------------------------------------------
# The P06 shape, in spirit.
# ---------------------------------------------------------------------------

class TestP06Shape:
    """A block whose ungapped identity is poor and whose gapped identity is
    high: the flat M was the bug, the gapped CIGAR is the fix."""

    def test_poor_ungapped_high_gapped(self, index):
        clip = GENOME_SEQ[P_ACC:P_ACC + 30]
        q = clip[:8] + clip[9:]                       # 1-bp deletion at 8
        # every base past position 8 is off by one, so a flat M reads badly
        ungapped = sum(1 for a, b in zip(q, GENOME_SEQ[P_ACC:P_ACC + len(q)])
                       if a.upper() == b.upper())
        assert ungapped / len(q) < 0.7, (
            'the control is broken: this block is not hard for a flat M')

        query = GENOME_SEQ[P_DON - 60:P_DON] + q
        r = _read('p06', query, [(0, 60), (4, len(q))], P_DON - 60)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert junction_of(r) == (P_DON, P_ACC)
        ops = block_ops(r, _RIGHT)
        assert any(op in (1, 2) for op, _ in ops), (
            f'no indel emitted for a block that needs one: {ops}')
        m, t = block_identity(r, _RIGHT, GENOME_SEQ)
        assert m / t >= 0.90, f'emitted block identity {m}/{t}'


# ---------------------------------------------------------------------------
# Degeneracy: fall back to the flat M rather than write a worse alignment.
# ---------------------------------------------------------------------------

class TestFallback:
    def test_junction_adjacent_deletion_falls_back(self, index):
        # exon-2 head missing its FIRST base: the only gapped spelling starts
        # with a D, which would silently make the intron 1 bp longer than the
        # N this call is writing. Refuse it and keep the flat M.
        clip = GENOME_SEQ[P_ACC:P_ACC + 30]
        query = GENOME_SEQ[P_DON - 60:P_DON] + clip[1:]
        r = _read('adj_del', query, [(0, 60), (4, 29)], P_DON - 60)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert stats.extra.get('emit_fallback_flat') == 1
        assert r.cigartuples == [(0, 60), (3, P_ACC - P_DON), (0, 29)]
        assert junction_of(r) == (P_DON, P_ACC)

    def test_junction_adjacent_insertion_falls_back_and_stays_shiftable(self, index):
        # An extra base AT the exon start. A gap at the ANCHORED end is the
        # aligner saying the block's first bases do not belong at intron_end —
        # a register shift, not an exon insertion — so `N 1I 30M` would bake in
        # a boundary the read does not support. Worse, it would be permanent:
        # the arbiter's Case A shift families need M ops on BOTH flanks
        # (`flank_m`), and the CIGAR shape persists in the output BAM, so the
        # junction could never be repaired on any later pass. Falling back
        # reproduces the old flat M and leaves the junction arbitrable.
        clip = GENOME_SEQ[P_ACC:P_ACC + 30]
        query = GENOME_SEQ[P_DON - 60:P_DON] + 'T' + clip
        r = _read('adj_ins', query, [(0, 60), (4, 31)], P_DON - 60)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert stats.extra.get('emit_fallback_flat') == 1
        assert r.cigartuples == [(0, 60), (3, P_ACC - P_DON), (0, 31)]
        ct = list(r.cigartuples)
        i = next(j for j, (op, _) in enumerate(ct) if op == 3)
        assert ct[i - 1][0] in (0, 7, 8) and ct[i + 1][0] in (0, 7, 8), (
            'the junction lost its M flanks and is no longer arbitrable: '
            f'{ct}')

    def test_scaled_allowance_admits_what_a_flat_five_refused(self, index):
        # 3 bp deleted at 15 and 3 bp inserted at 35 of a 60-bp block: 6 bp of
        # indel. The old flat _EMIT_MAX_INDEL=5 refused this; the block-scaled
        # allowance (ceil(max_edit_frac * m) = 12 here) spells it out, which is
        # the point — acceptance had already admitted up to that many edits.
        clip = GENOME_SEQ[P_ACC:P_ACC + 60]
        q = clip[:15] + clip[18:35] + 'TTT' + clip[35:]
        query = GENOME_SEQ[P_DON - 60:P_DON] + q
        r = _read('scaled', query, [(0, 60), (4, len(q))], P_DON - 60)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert stats.extra.get('emit_fallback_flat', 0) == 0
        ops = block_ops(r, _RIGHT)
        assert any(op in (1, 2) for op, _ in ops), ops
        # ...and the flat-5 allowance would still have refused it
        assert _block_cigar(_RIGHT, q, GENOME_SEQ, P_DON, P_ACC,
                            max_edit_frac=0.0) is None

    def test_long_block_with_scattered_indels_emits_gapped(self, index):
        # The case the scaled allowance exists for: a 200-bp ONT block at a
        # few percent indel. Eight scattered 1-bp events are far above a flat
        # 5 and far below ceil(0.2 * 200) = 40.
        head = list(GENOME_SEQ[P_ACC:P_ACC + 200])
        for pos in (170, 140, 110, 80):          # descending: indices stay valid
            del head[pos]
        for pos in (30, 60, 95, 125):
            head.insert(pos, 'T')
        q = ''.join(head)
        assert len(q) == 200
        query = GENOME_SEQ[P_DON - 60:P_DON] + q
        r = _read('long_gapped', query, [(0, 60), (4, 200)], P_DON - 60)
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert stats.extra.get('emit_fallback_flat', 0) == 0
        assert junction_of(r) == (P_DON, P_ACC)
        ops = block_ops(r, _RIGHT)
        assert sum(ln for op, ln in ops if op in (1, 2)) == 8, ops
        assert sum(ln for op, ln in ops if op in (0, 1)) == 200
        m, t = block_identity(r, _RIGHT, GENOME_SEQ)
        assert m / t >= 0.95, f'emitted block identity {m}/{t}'

    def test_clean_reads_never_count_a_fallback(self, index):
        query = GENOME_SEQ[P_DON - 30:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('clean', query, [(4, 30), (0, 60)], P_ACC)
        _, stats = _resolve(r, index)
        assert 'emit_fallback_flat' not in stats.extra


# ---------------------------------------------------------------------------
# Unit-level pins on the emitter itself.
# ---------------------------------------------------------------------------

class TestTagContract:
    """The module docstring promises XJ/XB are written ONLY on records the
    resolver changed, with no sentinel on an untouched read — a census that
    read "tag absent" as "rewritten with no move" would be counting wrong."""

    def test_untouched_read_carries_neither_tag(self, index):
        r = _read('linear', GENOME_SEQ[100:220], [(0, 120)], 100)
        changed, _ = _resolve(r, index)
        assert not changed
        assert not r.has_tag('XJ') and not r.has_tag('XB')

    def test_refused_clip_carries_neither_tag(self, index):
        query = 'A' * 30 + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('polya', query, [(4, 30), (0, 60)], P_ACC)
        changed, _ = _resolve(r, index)
        assert not changed
        assert not r.has_tag('XJ') and not r.has_tag('XB')

    def test_resolved_clip_carries_xj(self, index):
        query = GENOME_SEQ[P_DON - 30:P_DON] + GENOME_SEQ[P_ACC:P_ACC + 60]
        r = _read('resolved', query, [(4, 30), (0, 60)], P_ACC)
        assert _resolve(r, index)[0]
        assert r.get_tag('XJ') == f'{P_DON}-{P_ACC}:0.0:{_LEFT}'


class TestBlockCigarUnit:
    def test_gap_free_block_is_a_single_m(self):
        q = GENOME_SEQ[P_ACC:P_ACC + 40]
        assert _block_cigar(_RIGHT, q, GENOME_SEQ, P_DON, P_ACC) == (
            [(0, 40)], P_ACC)
        q = GENOME_SEQ[P_DON - 40:P_DON]
        assert _block_cigar(_LEFT, q, GENOME_SEQ, P_DON, P_ACC) == (
            [(0, 40)], P_DON - 40)

    def test_block_ref_start_lands_on_the_junction(self):
        # LEFT: block_ref_start + reference consumption must be intron_start
        q = GENOME_SEQ[P_DON - 40:P_DON]
        q = q[:20] + q[21:]                      # 39 query bases, 40 ref bases
        ops, rs = _block_cigar(_LEFT, q, GENOME_SEQ, P_DON, P_ACC)
        span = sum(ln for op, ln in ops if op in (0, 2))
        assert rs + span == P_DON
        assert sum(ln for op, ln in ops if op in (0, 1)) == len(q)

    def test_any_junction_adjacent_gap_is_refused(self):
        # A gap at the ANCHORED end is a register shift, not an exon indel:
        # refuse both ops, both sides. (All four shapes are reachable — a
        # 4,000-sample random-mutation sweep of 20-45 bp blocks produced
        # RIGHT leading I/D in ~6% of blocks and LEFT trailing I/D in ~0.2%.)
        cases = [
            (_RIGHT, GENOME_SEQ[P_ACC + 1:P_ACC + 41]),      # leading D
            (_RIGHT, 'T' + GENOME_SEQ[P_ACC:P_ACC + 40]),    # leading I
            (_LEFT, GENOME_SEQ[P_DON - 40:P_DON - 3]),       # trailing D
        ]
        for side, q in cases:
            assert _block_cigar(side, q, GENOME_SEQ, P_DON, P_ACC) is None, (
                f'{side} block {q[:8]}... should have fallen back')

    def test_interior_gaps_are_still_emitted(self):
        # the control for the rule above: one base further in and it is fine
        q = GENOME_SEQ[P_ACC:P_ACC + 1] + GENOME_SEQ[P_ACC + 2:P_ACC + 41]
        ops, rs = _block_cigar(_RIGHT, q, GENOME_SEQ, P_DON, P_ACC)
        assert ops[0][0] == 0 and any(op == 2 for op, _ in ops) and rs == P_ACC

    def test_allowance_scales_with_the_block(self):
        # a 60-bp block with a 10-bp deletion: refused at a flat 5, emitted
        # once the allowance follows max_edit_frac (ceil(0.2*60) = 12)
        q = GENOME_SEQ[P_ACC:P_ACC + 30] + GENOME_SEQ[P_ACC + 40:P_ACC + 70]
        assert _block_cigar(_RIGHT, q, GENOME_SEQ, P_DON, P_ACC,
                            max_edit_frac=0.0) is None
        ops, _ = _block_cigar(_RIGHT, q, GENOME_SEQ, P_DON, P_ACC,
                              max_edit_frac=0.2)
        assert sum(ln for op, ln in ops if op in (1, 2)) == 10

    def test_one_huge_deletion_still_falls_back(self):
        # 100-bp block, single 40-bp D: allowance is ceil(0.2*100) = 20, and
        # the reference window (m + allowance) cannot host the gap either —
        # this is the "wildly wrong placement" the bound exists for.
        q = GENOME_SEQ[P_ACC:P_ACC + 50] + GENOME_SEQ[P_ACC + 90:P_ACC + 140]
        assert len(q) == 100
        assert _block_cigar(_RIGHT, q, GENOME_SEQ, P_DON, P_ACC,
                            max_edit_frac=0.2) is None

    def test_empty_and_out_of_range_blocks_are_refused(self):
        assert _block_cigar(_RIGHT, '', GENOME_SEQ, P_DON, P_ACC) is None
        # a block longer than the reference left of intron_start
        assert _block_cigar(_LEFT, 'A' * 40, GENOME_SEQ, 10, 300) is None
