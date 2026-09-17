"""Annotated-canonical evidence gate in Module 2H (the R1 class).

An annotated + canonical junction is the strongest prior RECTIFY has.  Module 2H
ranked candidates on the raw edit-distance score whenever the incumbent was
canonical, so a sub-noise-floor win (measured on human DRS: 0.031 / 0.434 /
0.463 edit-distance units) could take an annotated GT-AG junction onto a novel
non-canonical one — the panel's R1 class.  The scoring policy in
``refine_read_junctions`` already said the opposite ("within the sub-integer HP
noise floor canonical annotated junctions are preferred"); only the
``tier_beats_alt`` branch implemented it.

These tests drive ``refine_read_junctions`` with a stubbed scorer so the margin
is exact and the assertions are about the GATE, not about DP arithmetic.

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
GLEN = 500
REF_START = 150
# 50M 100N 50M: exon1 [150,200), intron [200,300), exon2 [300,350)
CIGAR = [(0, 50), (3, 100), (0, 50)]

INCUMBENT = (200, 300)          # annotated, canonical GT-AG
ALT_NONCANON = (200, 301)       # novel, non-canonical acceptor
ALT_ANNOT_CANON = (200, 320)    # annotated, canonical GT-AG (isoform swap)


def _genome():
    g = list("C" * GLEN)
    g[200:202] = list("GT")     # donor, shared by every candidate
    g[298:300] = list("AG")     # acceptor of the incumbent
    g[300] = "C"                # -> [200,301) acceptor is "GC"... see below
    g[299] = "G"
    g[318:320] = list("AG")     # acceptor of the isoform-swap candidate
    return "".join(g)


GENOME = _genome()


def _read(cigar=CIGAR):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )
    r = pysam.AlignedSegment(header)
    r.query_name = "gate_read"
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.cigartuples = cigar
    n = sum(l for op, l in cigar if op in (0, 1, 4, 7, 8))
    r.query_sequence = ("ACGT" * (n // 4 + 1))[:n]
    r.query_qualities = pysam.qualitystring_to_array("I" * n)
    return r


@pytest.fixture
def policy_only_surgery(monkeypatch):
    """POLICY-ONLY SCOPE.  The tests that take this fixture assert the gate's
    DECISION on a synthetic ``50M100N50M`` read over a filler genome, where the
    moves under test are not writable: a boundary shift on such a read can only
    be realized with a compensating I/D beside the N, which ISSUE-031 refuses.
    Since 2H ranks only candidates the surgery can write, the probe is stubbed
    to "writable" here so the assertion stays about the gate.  Whether a
    decision can be MATERIALIZED is tested separately: the writable controls in
    ``TestWritableControls`` below (reads carrying the deletion the move
    absorbs, final CIGAR asserted) and tests/test_2h_realizable_ranking.py."""
    monkeypatch.setattr(jr, "_realizable", lambda *a, **k: True)


def _run(monkeypatch, scores, annotated, *, motif_blind=False, hold=None, cigar=CIGAR):
    """Refine one read with a stubbed scorer; return the replacement list."""
    def fake_score(query, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0

    monkeypatch.setattr(jr, "_score_junction", fake_score)
    if hold is not None:
        monkeypatch.setattr(jr, "_ANNOTATED_CANONICAL_HOLD", hold)

    pool = {(CHROM, s, e) for s, e in scores}
    idx = jr._build_junction_index(pool)
    annotated_set = {(CHROM, s, e) for s, e in annotated}
    return jr.refine_read_junctions(
        _read(cigar), idx, annotated_set, GENOME, "+",
        boundary_error_window=0,     # score every N-op; the filter is not under test
        motif_blind=motif_blind,
    )


def test_canonical_tier_of_the_fixture_is_what_the_tests_assume():
    """Guard the fixture itself: incumbent canonical, +1 shift non-canonical."""
    assert jr._canonical_tier(*INCUMBENT, GENOME, "+") < 4
    assert jr._canonical_tier(*ALT_NONCANON, GENOME, "+") >= 4
    assert jr._canonical_tier(*ALT_ANNOT_CANON, GENOME, "+") < 4


def test_sub_noise_floor_win_cannot_move_an_annotated_canonical_junction(monkeypatch):
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_NONCANON: 0.2},   # margin 0.4 < 1.0
        annotated={INCUMBENT},
    )
    assert repl == []


def test_full_edit_distance_unit_of_evidence_still_moves_it(monkeypatch, policy_only_surgery):
    repl = _run(
        monkeypatch,
        {INCUMBENT: 1.5, ALT_NONCANON: 0.2},   # margin 1.3 >= 1.0
        annotated={INCUMBENT},
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]


def test_gate_boundary_is_exactly_one_unit(monkeypatch, policy_only_surgery):
    """margin == 1.0 is enough; anything less is not."""
    assert _run(monkeypatch, {INCUMBENT: 1.0, ALT_NONCANON: 0.0},
                annotated={INCUMBENT}) != []
    assert _run(monkeypatch, {INCUMBENT: 0.999, ALT_NONCANON: 0.0},
                annotated={INCUMBENT}) == []


def test_isoform_swap_to_another_annotated_canonical_junction_is_not_gated(monkeypatch, policy_only_surgery):
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_ANNOT_CANON: 0.2},   # margin 0.4, but target is both
        annotated={INCUMBENT, ALT_ANNOT_CANON},
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_ANNOT_CANON]


def test_novel_incumbent_is_not_gated(monkeypatch, policy_only_surgery):
    """The corrections 2H exists for (novel -> annotated) must stay untouched."""
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_NONCANON: 0.2},
        annotated=set(),                          # incumbent is NOT annotated
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]


def test_motif_blind_bypasses_the_gate(monkeypatch, policy_only_surgery):
    """Station B decides on read evidence alone by construction."""
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_NONCANON: 0.2},
        annotated={INCUMBENT},
        motif_blind=True,
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]


def test_hold_can_be_disabled(monkeypatch, policy_only_surgery):
    """RECTIFY_ANNOT_CANON_HOLD=0 restores the pre-2026-09 behaviour."""
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_NONCANON: 0.2},
        annotated={INCUMBENT},
        hold=0.0,
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]


# ---------------------------------------------------------------------------
# Writable controls: the same gate, on moves the surgery CAN write
# ---------------------------------------------------------------------------

class TestWritableControls:
    """The policy tests above stub the surgery.  Here the moves are PURE SLIDES
    the fast path writes without any indel — the k genome bases that swap role
    between intron edge and exon edge are identical at both placements — so
    the gate decides on a writable move and the final CIGAR is asserted.  (A
    read carrying the deletion the move absorbs would be the D/N-merge twin,
    which is annotation's call, not the gate's — tests/test_dn_merge_twin.py.)

    Genome: incumbent [200,300) GT..CAG (annotated, canonical); the 10 bases
    [200,210) and [300,310) are both ``GTCCCCCCAG``, so
      slide +2  -> [202,302): donor CC, acceptor GT   = novel, non-canonical
      slide +10 -> [210,310): donor GT, acceptor CAG  = annotated, canonical
    """

    SLIDE2 = (202, 302)
    SLIDE10 = (210, 310)

    @staticmethod
    def _genome():
        g = list("C" * GLEN)
        g[200:210] = list("GTCCCCCCAG")
        g[300:310] = list("GTCCCCCCAG")
        g[210:212] = list("GT")
        g[298:300] = list("AG")
        return "".join(g)

    def _run(self, monkeypatch, scores, annotated):
        def fake_score(query, q_split, js, je, genome_seq, **kw):
            return scores[(js, je)], 0
        monkeypatch.setattr(jr, "_score_junction", fake_score)
        idx = jr._build_junction_index({(CHROM, s, e) for s, e in scores})
        r = _read()
        repl = jr.refine_read_junctions(r, idx, {(CHROM, s, e) for s, e in annotated},
                                        self._genome(), "+", boundary_error_window=0)
        out, applied = jr._apply_replacements_to_read(r, repl, self._genome(), "+", 0.25, 15)
        return repl, (out.cigarstring if applied else None)

    def test_frame(self):
        g = self._genome()
        assert jr._canonical_tier(*INCUMBENT, g, "+") == 0
        assert jr._canonical_tier(*self.SLIDE2, g, "+") >= 4
        assert jr._canonical_tier(*self.SLIDE10, g, "+") == 0
        assert jr._realizable(_read(), 1, *INCUMBENT, *self.SLIDE2, g, "+", 0.25, 15)
        assert jr._realizable(_read(), 1, *INCUMBENT, *self.SLIDE10, g, "+", 0.25, 15)

    def test_sub_noise_floor_win_is_held_even_when_writable(self, monkeypatch):
        repl, cig = self._run(monkeypatch, {INCUMBENT: 0.6, self.SLIDE2: 0.2},
                              annotated={INCUMBENT})
        assert repl == [] and cig is None

    def test_a_full_unit_moves_it_and_the_cigar_follows(self, monkeypatch):
        repl, cig = self._run(monkeypatch, {INCUMBENT: 1.5, self.SLIDE2: 0.2},
                              annotated={INCUMBENT})
        assert [(r[3], r[4]) for r in repl] == [self.SLIDE2]
        assert cig == "52M100N48M"

    def test_an_isoform_swap_is_not_gated_and_the_cigar_follows(self, monkeypatch):
        repl, cig = self._run(monkeypatch, {INCUMBENT: 0.6, self.SLIDE10: 0.2},
                              annotated={INCUMBENT, self.SLIDE10})
        assert [(r[3], r[4]) for r in repl] == [self.SLIDE10]
        assert cig == "60M100N40M"
