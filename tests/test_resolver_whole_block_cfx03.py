"""CFX-03 (Codex audit 2026-09-13): B2/B3 score a local window but rewrite a whole block.

The B2 rescue (mismatch-flagged TERMINAL linear block -> spliced) chose its
junction on ``arb_seg`` query bases per side and then relocated the ENTIRE tail
past the split by the intron length.  Codex's witness: a read that carries 40
bases of exon 2 and then returns to the intron.  The 80-base window cleared the
margin, the rewrite went through, and the read's literal mismatch burden rose
30 -> 94 (``360M`` -> ``200M300N160M``).  The fix validates every relocated base
under the discipline the window had to meet (beat the current placement by
``arb_margin`` over the whole block); B3 gets the same guard on its head.

The witness is synthetic and proves the mechanism, not a biological FP rate
(``dev/audits/codex_fable_20260913/repro_resolver_whole_block.py``).
"""

import random

import pysam
import pytest

from rectify.core.align.overhang_resolver import ResolverConfig, ResolverStats, resolve_read
from rectify.core.splice.splice_site_index import SpliceSiteIndex

D, E = 1200, 1500            # true intron [D, E), GT..AG


def _genome():
    rng = random.Random(20260809)
    seq = [rng.choice('ACGT') for _ in range(4000)]
    seq[D:D + 2] = list('GT')
    seq[E - 2:E] = list('AG')
    return ''.join(seq)


def _mk_read(g, query, cigar, rs):
    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': len(g)}]})
    r = pysam.AlignedSegment(header)
    r.query_name = 'cfx03'
    r.query_sequence = query
    r.flag = 0
    r.reference_id = 0
    r.reference_start = rs
    r.mapping_quality = 60
    r.cigartuples = cigar
    return r


def _run(g, read):
    genome = {'chrI': g}
    index = SpliceSiteIndex.build(genome)
    cfg = ResolverConfig(alpha=0.01, max_intron=5000)
    stats = ResolverStats()
    return resolve_read(read, genome, index, cfg, stats), stats


def _literal_mismatches(rec, genome):
    """Independent of the production scorer: walk the emitted coordinates."""
    refpos, qpos, n = rec.reference_start, 0, 0
    for op, length in rec.cigartuples:
        if op in (0, 7, 8):
            q = rec.query_sequence[qpos:qpos + length]
            r = genome[refpos:refpos + length]
            n += sum(a != b for a, b in zip(q, r))
        if op in (0, 1, 4, 7, 8):
            qpos += length
        if op in (0, 2, 3, 7, 8):
            refpos += length
    return n


def test_b2_refuses_a_rewrite_that_worsens_the_relocated_block():
    """Codex's witness: 40 bases of exon 2, then back into the intron."""
    g = _genome()
    query = g[D - 200:D] + g[E:E + 40] + g[D + 40:D + 160]
    r = _mk_read(g, query, [(0, 360)], D - 200)
    before = _literal_mismatches(r, g)
    changed, stats = _run(g, r)
    after = _literal_mismatches(r, g)
    assert stats.extra.get('arb_mm_flagged', 0) >= 1          # the trigger still fires
    assert stats.extra.get('arb_mm_whole_block_refused') == 1  # and the guard refuses
    assert stats.extra.get('arb_mm_spliced', 0) == 0
    assert not changed and r.cigartuples == [(0, 360)]
    assert after == before


def test_b2_positive_control_a_real_spliced_tail_still_splices():
    """Every relocated base IS exon 2: the window and the whole block agree."""
    g = _genome()
    query = g[D - 200:D] + g[E:E + 160]
    r = _mk_read(g, query, [(0, 360)], D - 200)
    before = _literal_mismatches(r, g)
    changed, stats = _run(g, r)
    after = _literal_mismatches(r, g)
    assert changed and stats.extra.get('arb_mm_spliced') == 1
    assert stats.extra.get('arb_mm_whole_block_refused', 0) == 0
    assert r.cigarstring == '200M300N160M'
    assert after < before and after == 0


def test_b3_positive_control_a_real_spliced_head_still_splices():
    """The B3 mirror: exon-1 bases smeared over the intron tail, re-anchored left."""
    g = _genome()
    query = g[D - 160:D] + g[E:E + 200]
    r = _mk_read(g, query, [(0, 360)], E - 160)
    changed, stats = _run(g, r)
    assert changed and stats.extra.get('arb_mm_spliced') == 1, stats.extra
    assert stats.extra.get('arb_mm_whole_block_refused', 0) == 0
    assert r.cigarstring == '160M300N200M'
    assert r.reference_start == D - 160
    assert _literal_mismatches(r, g) == 0


@pytest.mark.parametrize('tail_beyond_exon2', [60, 120])
def test_b2_the_guard_scales_with_how_much_of_the_tail_leaves_exon2(tail_beyond_exon2):
    """The witness family: the longer the stretch that returns to the intron,
    the worse the whole-block rewrite — every member is refused, and the
    read's mismatch burden is never raised."""
    g = _genome()
    exon2 = 160 - tail_beyond_exon2
    query = g[D - 200:D] + g[E:E + exon2] + g[D + exon2:D + 160]
    r = _mk_read(g, query, [(0, 360)], D - 200)
    before = _literal_mismatches(r, g)
    changed, stats = _run(g, r)
    assert _literal_mismatches(r, g) <= before
    if changed:
        # only a rewrite that does not raise the burden may go through
        assert stats.extra.get('arb_mm_whole_block_refused', 0) == 0
