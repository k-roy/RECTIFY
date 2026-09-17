"""Where Module 2H draws the line between "canonical" and "not".

`_canonical_tier` returns `t3` (0..4) when the donor is GT/GC and `4 + t3` when
it is not, and `t3 <= 1` is exactly "the acceptor ends in AG" (0 = YAG,
1 = RAG).  So `tier <= _CANONICAL_TIER_MAX` means a proper GT/GC..AG junction.

That line was drawn at `4` in three places, which asks only whether the DONOR is
intact.  A GT-AT (tier 3) or GT-GG (tier 2) N-op — canonical donor, no AG
acceptor at all — was therefore treated as canonical three times over:

  1. the clean-boundary pre-filter never scored it (no indel near the junction
     and "tier < 4", so 2H did not look at the read at all);
  2. it kept `is_alt` priority over a real GT-AG alternative at equal score;
  3. it collected the canonical prior itself, so the GT-AG alternative could not
     out-score it even with the discount.

Measured on the Sumner panel: bc1283c7 (GT-AT) and 12b2bc34 (GT-GG) each sat
6-7 nt from an annotated GT-AG and would not move.  Panel FN 3 -> 1, TP 8 -> 10,
FP unchanged at 23.

Author: Kevin R. Roy (agent S1)
"""

import sys
from pathlib import Path

import pysam
import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice import junction_refiner as jr  # noqa: E402
from rectify.core.splice.junction_scoring import (  # noqa: E402
    _CANONICAL_HP_PRIOR,
    _CANONICAL_TIER_MAX,
)

CHROM = "chrT"
GLEN = 600
REF_START = 150
CLEAN_CIGAR = [(0, 50), (3, 100), (0, 50)]   # no indel anywhere near the junction

INCUMBENT = (200, 300)
ALTERNATIVE = (200, 306)     # annotated, CAG acceptor -> tier 0


def _genome(acceptor_trinuc, alt_trinuc="CAG"):
    """Genome whose incumbent acceptor context is *acceptor_trinuc* (RNA order).

    *alt_trinuc* is the acceptor context of ALTERNATIVE, so a test can make the
    alternative canonical (the default, YAG -> tier 0) or not.
    """
    g = list("C" * GLEN)
    g[200:202] = list("GT")                  # donor, shared by both junctions
    g[297:300] = list(acceptor_trinuc)       # incumbent acceptor context
    g[303:306] = list(alt_trinuc)            # the alternative's acceptor context
    return "".join(g)


# incumbent acceptor context -> expected tier
CONTEXTS = {
    "CAG": 0,    # YAG  — a proper junction
    "AAG": 1,    # RAG  — still a proper junction
    "CGG": 2,    # NBG  — no AG acceptor  (the 12b2bc34 shape)
    "CAT": 3,    # NAT  — no AG acceptor  (the bc1283c7 shape)
}


def _read(cigar=CLEAN_CIGAR):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )
    r = pysam.AlignedSegment(header)
    r.query_name = "tier_read"
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.cigartuples = cigar
    n = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    r.query_sequence = ("ACGT" * (n // 4 + 1))[:n]
    r.query_qualities = pysam.qualitystring_to_array("I" * n)
    return r


def _winner(monkeypatch, genome, scores, annotated=(ALTERNATIVE,), window=10, cigar=CLEAN_CIGAR):
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0

    monkeypatch.setattr(jr, "_score_junction", fake_score)
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
    repl = jr.refine_read_junctions(
        _read(cigar), idx, {(CHROM, s, e) for s, e in annotated}, genome, "+",
        boundary_error_window=window,
    )
    return (repl[0][3], repl[0][4]) if repl else None



@pytest.fixture
def policy_only_surgery(monkeypatch):
    """POLICY-ONLY SCOPE.  The tests that take this fixture assert a ranking /
    gate DECISION on a synthetic read over a filler genome, where the move under
    test is not writable: a boundary shift on such a read can only be realized
    with a compensating I/D beside the N, which ISSUE-031 refuses.  Since 2H
    ranks only candidates the surgery can write, the probe is stubbed to
    "writable" here so the assertion stays about the decision.  Whether a
    decision can be MATERIALIZED is tested separately: the writable controls at
    the bottom of this file (reads carrying the deletion the move absorbs, final
    CIGAR asserted) and tests/test_2h_realizable_ranking.py."""
    monkeypatch.setattr(jr, "_realizable", lambda *a, **k: True)


# ---------------------------------------------------------------------------
# What the constant means
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("context,expected_tier", sorted(CONTEXTS.items()))
def test_tier_of_each_acceptor_context(context, expected_tier):
    assert jr._canonical_tier(*INCUMBENT, _genome(context), "+") == expected_tier


@pytest.mark.parametrize("context,tier", sorted(CONTEXTS.items()))
def test_the_line_is_exactly_acceptor_ends_in_AG(context, tier):
    """`tier <= _CANONICAL_TIER_MAX` must coincide with "the acceptor is AG"."""
    assert (tier <= _CANONICAL_TIER_MAX) == context.endswith("AG")


def test_the_alternative_is_canonical():
    assert jr._canonical_tier(*ALTERNATIVE, _genome("CAG"), "+") == 0


# ---------------------------------------------------------------------------
# 1. the clean-boundary pre-filter
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("context", ["CGG", "CAT"])
def test_a_clean_non_AG_acceptor_is_still_scored(context):
    """No indel near the junction, so only the tier bypass can let 2H look."""
    g = _genome(context)
    read = _read()
    assert not jr._has_boundary_error(read.cigartuples, *INCUMBENT,
                                      read.reference_start, 10)
    assert jr._canonical_tier(*INCUMBENT, g, "+") > _CANONICAL_TIER_MAX


@pytest.mark.parametrize("context", ["CGG", "CAT"])
def test_a_clean_non_AG_acceptor_moves_to_the_annotated_GT_AG(monkeypatch, context, policy_only_surgery):
    """The bc1283c7 / 12b2bc34 case, end to end, with the real pre-filter on.

    The alternative even scores slightly WORSE on read evidence (as it did for
    bc1283c7: 1.000 vs 0.834); the canonical prior is what carries it.
    """
    assert _winner(monkeypatch, _genome(context),
                   {INCUMBENT: 0.834, ALTERNATIVE: 1.000}) == ALTERNATIVE


@pytest.mark.parametrize("context", ["CAG", "AAG"])
def test_a_clean_proper_junction_is_still_skipped(monkeypatch, context):
    """A GT/GC..AG N-op with a clean flank is not 2H's business — unchanged."""
    assert _winner(monkeypatch, _genome(context),
                   {INCUMBENT: 0.834, ALTERNATIVE: 0.000}) is None


# ---------------------------------------------------------------------------
# 2. is_alt priority
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("context", ["CGG", "CAT"])
def test_a_non_AG_incumbent_forfeits_is_alt_priority(monkeypatch, context, policy_only_surgery):
    """At an EXACT tie the canonical alternative wins (the 12b2bc34 shape)."""
    assert _winner(monkeypatch, _genome(context),
                   {INCUMBENT: 1.407, ALTERNATIVE: 1.407}, window=0) == ALTERNATIVE


@pytest.mark.parametrize("context", ["CAG", "AAG"])
def test_a_proper_incumbent_keeps_is_alt_priority(monkeypatch, context):
    """The FN guard for the other direction: YAG and RAG both hold their ground."""
    assert _winner(monkeypatch, _genome(context),
                   {INCUMBENT: 1.407, ALTERNATIVE: 1.407}, window=0) is None


# ---------------------------------------------------------------------------
# 3. who collects the canonical prior
# ---------------------------------------------------------------------------

def test_a_non_AG_candidate_does_not_collect_the_canonical_prior(monkeypatch):
    """A candidate whose OWN acceptor is not AG must not be discounted into winning.

    Both placements are non-canonical here, so the prior is in play (the
    incumbent is tier 3) but neither side may collect it.  The alternative
    scores 0.3 WORSE — less than the 0.5 prior — so it would win if it were
    still being discounted, and must not.
    """
    g = _genome("CAT", alt_trinuc="CGG")        # incumbent tier 3, alternative tier 2
    assert jr._canonical_tier(*INCUMBENT, g, "+") > _CANONICAL_TIER_MAX
    assert jr._canonical_tier(*ALTERNATIVE, g, "+") > _CANONICAL_TIER_MAX
    gap = 0.3
    assert gap < _CANONICAL_HP_PRIOR            # the discount would flip this
    assert _winner(monkeypatch, g,
                   {INCUMBENT: 0.834, ALTERNATIVE: 0.834 + gap}, window=0) is None


def test_the_prior_is_what_carries_a_worse_scoring_canonical_alternative(monkeypatch, policy_only_surgery):
    """Sanity: inside the prior the canonical alternative wins, outside it loses."""
    g = _genome("CAT")
    inside = {INCUMBENT: 0.834, ALTERNATIVE: 0.834 + _CANONICAL_HP_PRIOR - 0.01}
    outside = {INCUMBENT: 0.834, ALTERNATIVE: 0.834 + _CANONICAL_HP_PRIOR + 0.01}
    assert _winner(monkeypatch, g, inside) == ALTERNATIVE
    assert _winner(monkeypatch, g, outside) is None


# ---------------------------------------------------------------------------
# Writable controls: the same line, on a move the surgery CAN write
# ---------------------------------------------------------------------------

class TestWritableControls:
    """``50M100N6D44M`` -> acceptor +6 absorbs the deletion (``50M106N44M``), so
    the canonical-prior decision materializes for real and the CIGAR is
    asserted.  The 6D sits at the boundary, so the real pre-filter (window 10)
    sees a boundary error and scores the N-op on its own."""

    D6 = [(0, 50), (3, 100), (2, 6), (0, 44)]

    def _proposals(self, monkeypatch, genome, scores):
        def fake_score(query, q_split, js, je, genome_seq, **kw):
            return scores[(js, je)], 0
        monkeypatch.setattr(jr, "_score_junction", fake_score)
        idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
        r = _read(self.D6)
        repl = jr.refine_read_junctions(
            r, idx, {(CHROM, *ALTERNATIVE)}, genome, "+", boundary_error_window=10)
        return r, repl

    def test_the_move_is_writable_by_the_real_surgery(self):
        assert jr._realizable(_read(self.D6), 1, *INCUMBENT, *ALTERNATIVE,
                              _genome("CGG"), "+", 0.25, 15)

    @pytest.mark.parametrize("context", ["CGG", "CAT"])
    def test_a_non_AG_acceptor_moves_and_the_cigar_follows(self, monkeypatch, context):
        r, repl = self._proposals(monkeypatch, _genome(context),
                                  {INCUMBENT: 0.834, ALTERNATIVE: 1.000})
        assert [(x[3], x[4]) for x in repl] == [ALTERNATIVE]
        out, applied = jr._apply_replacements_to_read(r, repl, _genome(context), "+", 0.25, 15)
        assert applied and out.cigarstring == "50M106N44M"

    @pytest.mark.parametrize("context", ["CAG", "AAG"])
    def test_a_proper_incumbent_holds_even_when_the_move_is_writable(self, monkeypatch, context):
        r, repl = self._proposals(monkeypatch, _genome(context),
                                  {INCUMBENT: 1.407, ALTERNATIVE: 1.407})
        assert repl == []
