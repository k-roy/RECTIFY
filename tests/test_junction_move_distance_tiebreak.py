"""Module 2H tie-break: at a genuine tie, take the smaller move (ISSUE-005(b)).

The fifth key of the candidate ordering tuple used to be `abs(delta)` from
`_score_junction`, which is the applied SLIDE and is structurally 0 for every
candidate — so the slot was a constant and equal-scoring candidates were
separated by `(js, je)`, i.e. by the LOWER GENOMIC COORDINATE.  That is an
arbitrary leftward bias with nothing to do with the read.  The key is now
`move_dist = |js - ns| + |je - ne|`, the L1 distance from the aligner's own
placement.

It is the FIFTH key.  Score, canonical tier, is_alt and annotation all decide
first, and the tests below pin that ordering as hard as they pin the new
behavior — a tie-break that can override a policy key is a policy change.

Scores are stubbed so "a tie" is exact rather than an artifact of DP arithmetic.

Author: Kevin R. Roy (agent S1)
"""

import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_refiner as jr  # noqa: E402

CHROM = "chrT"
GLEN = 600
REF_START = 150
CIGAR = [(0, 50), (3, 100), (0, 50)]     # exon [150,200) intron [200,300) exon [300,350)

# Donors and acceptors are spaced >= 6 bp apart: `_canonical_tier` reads the 3'SS
# TRInucleotide, so an acceptor planted within 2 bp of another would change its
# neighbor's tier from CAG (tier 0) to GAG (tier 1) and silently break the "equal
# on every prior key" premise.  test_the_candidates_are_indistinguishable... is
# the guard that catches exactly that.
DONORS = (188, 200, 212)
ACCEPTORS = (282, 294, 300, 306, 318)

INCUMBENT = (200, 300)
NEAR = (200, 306)          # move_dist 6,  higher coordinate
FAR = (200, 318)           # move_dist 18, higher coordinate
NEAR_LEFT = (200, 294)     # move_dist 6,  lower coordinate
FAR_LEFT = (200, 282)      # move_dist 18, lower coordinate
LEFT = (188, 300)          # move_dist 12
RIGHT = (212, 300)         # move_dist 12, equidistant with LEFT
NEAR_NONCANON = (200, 303)  # move_dist 3, but the acceptor is not AG


def _genome():
    g = list("C" * GLEN)
    for s in DONORS:
        g[s:s + 2] = list("GT")
    for e in ACCEPTORS:
        g[e - 2:e] = list("AG")
    return "".join(g)


GENOME = _genome()


def _read():
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )
    r = pysam.AlignedSegment(header)
    r.query_name = "tiebreak_read"
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.cigartuples = CIGAR
    r.query_sequence = "ACGT" * 25
    r.query_qualities = pysam.qualitystring_to_array("I" * 100)
    return r


def _winner(monkeypatch, scores, annotated=frozenset()):
    """Refine with a stubbed scorer; return the winning (start, end) or None.

    The incumbent is deliberately NOT in `scores` (i.e. not in the pool) unless a
    test puts it there — that is the real situation in which the tie-break
    decides anything, and it keeps the annotated-canonical gate (which needs an
    incumbent score) out of the way.
    """
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0

    monkeypatch.setattr(jr, "_score_junction", fake_score)
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
    repl = jr.refine_read_junctions(
        _read(), idx, {(CHROM, s, e) for s, e in annotated}, GENOME, "+",
        boundary_error_window=0,
    )
    return (repl[0][3], repl[0][4]) if repl else None


# ---------------------------------------------------------------------------
# Fixture guard: the candidates really are equal on every key above move_dist
# ---------------------------------------------------------------------------

def test_the_candidates_are_indistinguishable_above_the_tiebreak():
    equal = (INCUMBENT, NEAR, FAR, NEAR_LEFT, FAR_LEFT, LEFT, RIGHT)
    tiers = {c: jr._canonical_tier(*c, GENOME, "+") for c in equal}
    assert len(set(tiers.values())) == 1, tiers
    assert all(t < 4 for t in tiers.values()), tiers


def test_move_distances_are_what_the_names_say():
    ns, ne = INCUMBENT
    dist = {c: abs(c[0] - ns) + abs(c[1] - ne)
            for c in (INCUMBENT, NEAR, FAR, NEAR_LEFT, FAR_LEFT, LEFT, RIGHT)}
    assert dist[INCUMBENT] == 0
    assert dist[NEAR] == dist[NEAR_LEFT] == 6
    assert dist[FAR] == dist[FAR_LEFT] == 18
    assert dist[LEFT] == dist[RIGHT] == 12


# ---------------------------------------------------------------------------
# The new behavior
# ---------------------------------------------------------------------------

def test_at_a_tie_the_nearer_candidate_wins(monkeypatch):
    assert _winner(monkeypatch, {NEAR: 0.5, FAR: 0.5}) == NEAR


def test_the_result_does_not_depend_on_genomic_order(monkeypatch):
    """The old key was "lowest coordinate wins".

    Above, NEAR and FAR are both to the RIGHT of the incumbent, so the lowest
    coordinate happens to be the nearer one — a leftward bias would pass by
    accident.  Mirrored: here the nearer candidate has the HIGHER coordinate, so
    only a real distance key gets this right.
    """
    assert _winner(monkeypatch, {NEAR_LEFT: 0.5, FAR_LEFT: 0.5}) == NEAR_LEFT


def test_equidistant_candidates_are_still_resolved_deterministically(monkeypatch):
    """move_dist cannot separate a symmetric pair; (js, je) still breaks it.

    Documented, not endorsed — the residual arbitrariness is now confined to
    exact ties in BOTH score and distance.
    """
    assert _winner(monkeypatch, {LEFT: 0.5, RIGHT: 0.5}) == LEFT


def test_the_incumbent_is_its_own_nearest(monkeypatch):
    """Consistency with is_alt: distance 0, so the two keys can never disagree."""
    assert _winner(monkeypatch, {INCUMBENT: 0.5, NEAR: 0.5}) is None


# ---------------------------------------------------------------------------
# It must never override a key above it
# ---------------------------------------------------------------------------

def test_a_better_score_beats_a_smaller_move(monkeypatch):
    assert _winner(monkeypatch, {NEAR: 0.9, FAR: 0.2}) == FAR


def test_annotation_beats_a_smaller_move(monkeypatch):
    assert _winner(monkeypatch, {NEAR: 0.5, FAR: 0.5}, annotated={FAR}) == FAR


def test_canonical_tier_beats_a_smaller_move(monkeypatch):
    """A nearer, WORSE-tier candidate must not beat a farther canonical one."""
    assert (jr._canonical_tier(*NEAR_NONCANON, GENOME, "+")
            > jr._canonical_tier(*FAR, GENOME, "+"))
    assert _winner(monkeypatch, {NEAR_NONCANON: 0.5, FAR: 0.5}) == FAR


def test_the_incumbent_still_wins_at_equal_score(monkeypatch):
    """is_alt outranks move_dist, and both agree here — belt and braces."""
    assert _winner(monkeypatch, {INCUMBENT: 0.5, NEAR: 0.5, FAR: 0.5}) is None


def test_move_dist_does_not_reopen_the_annotated_canonical_gate(monkeypatch):
    """250a608 is upstream of this key and must still hold.

    A NEARER (3 bp) novel worse-tier candidate that wins on raw score by less
    than 1.0 must still not move an annotated canonical incumbent.
    """
    assert _winner(
        monkeypatch,
        {INCUMBENT: 0.6, NEAR_NONCANON: 0.2},
        annotated={INCUMBENT},
    ) is None
