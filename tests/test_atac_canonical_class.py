"""AT-AC is a canonical intron class for every organism (Kevin, 2026-09-05).

`_canonical_tier` implemented the GT/GC..AG grammar only, so an annotated U12
AT-AC junction scored tier 8 — worse than a broken donor. Module 2H therefore
neither protected an AT-AC incumbent nor let an AT-AC candidate compete.

It is a PAIRED class and that is the whole point: AT..AC on the forward genome
is canonical only together (a minus-strand transcript reads it as forward
GT..AT). AT..AG and GT..AC are NOT members and must keep scoring as the broken
junctions they are — the pair carries the meaning, not either half.

Rank is `_ATAC_TIER = 1`: below GT-AG (tier 0), so an AT-AC candidate never
outranks a GT-AG one at equal read score, but inside `_CANONICAL_TIER_MAX`, so
the annotated-canonical evidence gate PROTECTS an annotated AT-AC incumbent and
the clean-boundary pre-filter treats it as a proper junction.

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
    _ATAC_FORWARD_BY_STRAND,
    _ATAC_TIER,
    _CANONICAL_HP_PRIOR,
    _CANONICAL_TIER_MAX,
)
from rectify.core.splice.overhang_informativeness import (  # noqa: E402
    _ATAC_DINUCS,
    canonical_in_class,
    is_canonical_junction,
)

CHROM = "chrT"
GLEN = 600
REF_START = 150
CIGAR = [(0, 50), (3, 100), (0, 50)]     # exon [150,200) intron [200,300) exon [300,350)

INTRON = (200, 300)          # the read's own N-op
ALT_A = (206, 306)           # a competing junction
ALT_B = (212, 312)           # a second competing junction


def _genome(*specs):
    """Build a genome from explicit (start, end, donor, acceptor) junctions.

    Ends are kept >= 3 bp apart so the filler 'C' before each acceptor gives a
    YAG context — otherwise a neighbouring acceptor would silently downgrade its
    neighbour from tier 0 to tier 1.
    """
    g = list("C" * GLEN)
    for start, end, donor, acceptor in specs:
        g[start:start + 2] = list(donor)
        g[end - 2:end] = list(acceptor)
    return "".join(g)


def _read(cigar=CIGAR):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )
    r = pysam.AlignedSegment(header)
    r.query_name = "atac_read"
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.cigartuples = cigar
    n = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    r.query_sequence = ("ACGT" * (n // 4 + 1))[:n]
    r.query_qualities = pysam.qualitystring_to_array("I" * n)
    return r


def _winner(monkeypatch, genome, scores, annotated=(), strand="+", window=0):
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0

    monkeypatch.setattr(jr, "_score_junction", fake_score)
    idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
    repl = jr.refine_read_junctions(
        _read(), idx, {(CHROM, s, e) for s, e in annotated}, genome, strand,
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
# The grammar itself
# ---------------------------------------------------------------------------

def test_the_strand_map_is_derived_from_the_single_source_of_truth():
    assert set(_ATAC_FORWARD_BY_STRAND.values()) == set(_ATAC_DINUCS)


@pytest.mark.parametrize("donor,acceptor,strand", [
    ("AT", "AC", "+"),      # a plus-strand transcript reads AT..AC directly
    ("GT", "AT", "-"),      # a minus-strand transcript: forward GT..AT
])
def test_atac_is_canonical_in_its_transcript_frame(donor, acceptor, strand):
    g = _genome((*INTRON, donor, acceptor))
    tier = jr._canonical_tier(*INTRON, g, strand)
    assert tier == _ATAC_TIER
    assert tier <= _CANONICAL_TIER_MAX


@pytest.mark.parametrize("donor,acceptor,strand", [
    ("AT", "AC", "-"),      # right pair, wrong frame
    ("GT", "AT", "+"),
])
def test_the_wrong_transcript_frame_is_not_canonical(donor, acceptor, strand):
    g = _genome((*INTRON, donor, acceptor))
    assert jr._canonical_tier(*INTRON, g, strand) > _CANONICAL_TIER_MAX


@pytest.mark.parametrize("donor,acceptor", [("AT", "AG"), ("GT", "AC")])
def test_half_an_atac_pair_is_still_broken(donor, acceptor):
    """The pair carries the meaning: neither half earns the class on its own."""
    g = _genome((*INTRON, donor, acceptor))
    tier = jr._canonical_tier(*INTRON, g, "+")
    assert tier > _CANONICAL_TIER_MAX
    assert tier >= 4


def test_atac_ranks_below_gt_ag():
    """Canonical, but never better than the major class at equal read score."""
    gt_ag = jr._canonical_tier(*INTRON, _genome((*INTRON, "GT", "AG")), "+")
    at_ac = jr._canonical_tier(*INTRON, _genome((*INTRON, "AT", "AC")), "+")
    assert gt_ag == 0
    assert gt_ag < at_ac <= _CANONICAL_TIER_MAX


def test_the_shared_grammar_helper_accepts_atac_by_default():
    g = _genome((*INTRON, "AT", "AC"))
    assert is_canonical_junction(g, *INTRON)
    assert canonical_in_class(g, *INTRON)
    assert not is_canonical_junction(g, *INTRON, atac=False)


# ---------------------------------------------------------------------------
# Consumer 1: the evidence gate PROTECTS an annotated AT-AC incumbent
# ---------------------------------------------------------------------------

def test_an_annotated_atac_incumbent_is_protected_by_the_evidence_gate(monkeypatch):
    """A sub-noise-floor win must not move it — the R1 guarantee, now for U12."""
    g = _genome((*INTRON, "AT", "AC"), (*ALT_A, "GT", "AG"))
    assert jr._canonical_tier(*INTRON, g, "+") <= _CANONICAL_TIER_MAX
    assert _winner(monkeypatch, g, {INTRON: 0.6, ALT_A: 0.2},
                   annotated=(INTRON,)) is None


def test_that_protection_still_yields_to_a_full_unit_of_evidence(monkeypatch, policy_only_surgery):
    g = _genome((*INTRON, "AT", "AC"), (*ALT_A, "GT", "AG"))
    assert _winner(monkeypatch, g, {INTRON: 1.5, ALT_A: 0.2},
                   annotated=(INTRON,)) == ALT_A


def test_a_clean_atac_junction_is_not_re_examined_by_the_pre_filter(monkeypatch):
    """A proper junction with a clean flank is not 2H's business, as for GT-AG."""
    g = _genome((*INTRON, "AT", "AC"), (*ALT_A, "GT", "AG"))
    assert _winner(monkeypatch, g, {INTRON: 0.834, ALT_A: 0.0}, window=10) is None


# ---------------------------------------------------------------------------
# Consumer 2: an AT-AC CANDIDATE competes at canonical rank
# ---------------------------------------------------------------------------

def test_an_atac_candidate_collects_the_canonical_prior(monkeypatch, policy_only_surgery):
    """Incumbent broken, so the prior is live; the AT-AC alternative scores
    WORSE on read evidence and is carried by the discount."""
    g = _genome((*INTRON, "AA", "CC"), (*ALT_A, "AT", "AC"))
    assert jr._canonical_tier(*INTRON, g, "+") > _CANONICAL_TIER_MAX
    assert jr._canonical_tier(*ALT_A, g, "+") == _ATAC_TIER
    assert _winner(monkeypatch, g,
                   {INTRON: 0.834,
                    ALT_A: 0.834 + _CANONICAL_HP_PRIOR - 0.01}) == ALT_A


def test_an_atac_candidate_outside_the_prior_still_loses(monkeypatch):
    g = _genome((*INTRON, "AA", "CC"), (*ALT_A, "AT", "AC"))
    assert _winner(monkeypatch, g,
                   {INTRON: 0.834,
                    ALT_A: 0.834 + _CANONICAL_HP_PRIOR + 0.01}) is None


def test_a_gt_ag_candidate_beats_an_atac_one_at_equal_score(monkeypatch, policy_only_surgery):
    """`_ATAC_TIER` is 1, so the major class takes the tie."""
    g = _genome((*INTRON, "AA", "CC"), (*ALT_A, "AT", "AC"), (*ALT_B, "GT", "AG"))
    assert jr._canonical_tier(*ALT_A, g, "+") == _ATAC_TIER
    assert jr._canonical_tier(*ALT_B, g, "+") == 0
    assert _winner(monkeypatch, g,
                   {INTRON: 0.834, ALT_A: 0.2, ALT_B: 0.2}) == ALT_B


def test_an_atac_candidate_beats_a_broken_one_at_equal_score(monkeypatch, policy_only_surgery):
    """And it outranks anything outside the grammar."""
    g = _genome((*INTRON, "AA", "CC"), (*ALT_A, "AT", "AC"), (*ALT_B, "AT", "AG"))
    assert jr._canonical_tier(*ALT_B, g, "+") >= 4
    assert _winner(monkeypatch, g,
                   {INTRON: 0.834, ALT_A: 0.2, ALT_B: 0.2}) == ALT_A


# ---------------------------------------------------------------------------
# Writable control: an AT-AC candidate carried by the prior, materialized
# ---------------------------------------------------------------------------

class TestWritableControls:
    """A PURE SLIDE the fast path writes without any indel: the 8 bases
    [200,208) and [300,308) are both ``ATCCCCAC``, so the incumbent [200,300)
    (AT..CC — HALF a pair, broken) slides +8 onto [208,308) (AT..AC, the pair,
    `_ATAC_TIER`).  The prior decides on a writable move and the CIGAR follows
    (``50M100N50M`` -> ``58M100N42M``)."""

    ALT = (208, 308)

    @staticmethod
    def _genome():
        g = list("C" * GLEN)
        g[200:208] = list("ATCCCCAC")
        g[300:308] = list("ATCCCCAC")
        g[208:210] = list("AT")
        g[298:300] = list("CC")
        return "".join(g)

    def _run(self, monkeypatch, scores):
        def fake_score(query, q_split, js, je, genome_seq, **kw):
            return scores[(js, je)], 0
        monkeypatch.setattr(jr, "_score_junction", fake_score)
        idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
        r = _read()
        repl = jr.refine_read_junctions(r, idx, set(), self._genome(), "+",
                                        boundary_error_window=0)
        out, applied = jr._apply_replacements_to_read(r, repl, self._genome(), "+", 0.25, 15)
        return repl, (out.cigarstring if applied else None)

    def test_frame(self):
        g = self._genome()
        assert jr._canonical_tier(*INTRON, g, "+") > _CANONICAL_TIER_MAX
        assert jr._canonical_tier(*self.ALT, g, "+") == _ATAC_TIER
        assert jr._realizable(_read(), 1, *INTRON, *self.ALT, g, "+", 0.25, 15)

    def test_an_atac_candidate_carried_by_the_prior_materializes(self, monkeypatch):
        repl, cig = self._run(monkeypatch,
                              {INTRON: 0.834, self.ALT: 0.834 + _CANONICAL_HP_PRIOR - 0.01})
        assert [(x[3], x[4]) for x in repl] == [self.ALT]
        assert cig == "58M100N42M"

    def test_outside_the_prior_it_still_loses_even_though_writable(self, monkeypatch):
        repl, cig = self._run(monkeypatch,
                              {INTRON: 0.834, self.ALT: 0.834 + _CANONICAL_HP_PRIOR + 0.01})
        assert repl == [] and cig is None
