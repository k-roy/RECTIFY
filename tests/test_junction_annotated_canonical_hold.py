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


def _read():
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]}
    )
    r = pysam.AlignedSegment(header)
    r.query_name = "gate_read"
    r.reference_id = 0
    r.reference_start = REF_START
    r.mapping_quality = 60
    r.cigartuples = CIGAR
    r.query_sequence = "ACGT" * 25
    r.query_qualities = pysam.qualitystring_to_array("I" * 100)
    return r


def _run(monkeypatch, scores, annotated, *, motif_blind=False, hold=None):
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
        _read(), idx, annotated_set, GENOME, "+",
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


def test_full_edit_distance_unit_of_evidence_still_moves_it(monkeypatch):
    repl = _run(
        monkeypatch,
        {INCUMBENT: 1.5, ALT_NONCANON: 0.2},   # margin 1.3 >= 1.0
        annotated={INCUMBENT},
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]


def test_gate_boundary_is_exactly_one_unit(monkeypatch):
    """margin == 1.0 is enough; anything less is not."""
    assert _run(monkeypatch, {INCUMBENT: 1.0, ALT_NONCANON: 0.0},
                annotated={INCUMBENT}) != []
    assert _run(monkeypatch, {INCUMBENT: 0.999, ALT_NONCANON: 0.0},
                annotated={INCUMBENT}) == []


def test_isoform_swap_to_another_annotated_canonical_junction_is_not_gated(monkeypatch):
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_ANNOT_CANON: 0.2},   # margin 0.4, but target is both
        annotated={INCUMBENT, ALT_ANNOT_CANON},
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_ANNOT_CANON]


def test_novel_incumbent_is_not_gated(monkeypatch):
    """The corrections 2H exists for (novel -> annotated) must stay untouched."""
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_NONCANON: 0.2},
        annotated=set(),                          # incumbent is NOT annotated
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]


def test_motif_blind_bypasses_the_gate(monkeypatch):
    """Station B decides on read evidence alone by construction."""
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_NONCANON: 0.2},
        annotated={INCUMBENT},
        motif_blind=True,
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]


def test_hold_can_be_disabled(monkeypatch):
    """RECTIFY_ANNOT_CANON_HOLD=0 restores the pre-2026-09 behaviour."""
    repl = _run(
        monkeypatch,
        {INCUMBENT: 0.6, ALT_NONCANON: 0.2},
        annotated={INCUMBENT},
        hold=0.0,
    )
    assert [(r[3], r[4]) for r in repl] == [ALT_NONCANON]
