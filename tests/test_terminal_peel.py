"""Tests for the Module 2F multi-hypothesis terminal junction peel.

Covers the peel-depth generator's stop conditions and the fill-in-only safety
guarantee: the peel must NEVER override an existing baseline rescue (that
behaviour regressed cat3_minus_1 / cat6_minus_1 when acceptance used edit
distance alone). The peel only ADDS a rescue where the baseline produced none.
"""

import os

import pysam

from rectify.core.splice.splice_aware_5prime import (
    _log_peel_screen,
    _peel_candidate_depths,
    _terminal_peel_rescue,
)


def _make_read(*, chrom="chrI", start=1000, seq="", cigar=((0, 0),), is_reverse=False):
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 100_000}]}
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "test"
    r.reference_name = chrom
    r.reference_start = start
    r.cigartuples = list(cigar)
    r.is_reverse = is_reverse
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


# --- _peel_candidate_depths -------------------------------------------------

def test_depths_stop_at_existing_n_op():
    """An existing N op halts depth generation (ref jump > 1)."""
    # M10 N50 M30; first 10 query bases mismatch genome (no clean anchor).
    read = _make_read(start=1000, seq="A" * 10 + "C" * 30,
                      cigar=((0, 10), (3, 50), (0, 30)))
    depths = _peel_candidate_depths(read, "C" * 2000, "+", max_peel=100, clean_anchor=10)
    assert depths == list(range(1, 11)), "should stop at the N after 10 query bases"


def test_depths_stop_at_clean_anchor():
    """`clean_anchor` consecutive genome-matching bases halt generation."""
    # M30, query == genome (all matches) -> stop once 10 consecutive matches seen.
    read = _make_read(start=1000, seq="C" * 30, cigar=((0, 30),))
    depths = _peel_candidate_depths(read, "C" * 2000, "+", max_peel=100, clean_anchor=10)
    assert depths == list(range(1, 11)), "should stop at the clean-anchor run of 10"


def test_depths_capped_at_max_peel():
    """With no N and no clean anchor, generation stops at max_peel."""
    # M30, query all mismatch genome -> consec_clean never reaches anchor.
    read = _make_read(start=1000, seq="A" * 30, cigar=((0, 30),))
    depths = _peel_candidate_depths(read, "C" * 2000, "+", max_peel=15, clean_anchor=10)
    assert depths == list(range(1, 16)), "should cap at max_peel=15"


# --- _terminal_peel_rescue (fill-in-only guarantee) -------------------------

def test_peel_never_overrides_a_rescued_baseline():
    """Fill-only: if the baseline already rescued, the peel returns None even
    when both gates (terminal evidence + nearby candidate 3'SS) pass."""
    # 10S 20M plus-strand read => terminal evidence present (Gate A).
    read = _make_read(start=1000, seq="A" * 10 + "C" * 20, cigar=((4, 10), (0, 20)))
    genome = {"chrI": "C" * 2000}
    cands = {("chrI", 980, 1000, "+")}  # intron_end=1000 == align_5prime => Gate B passes
    res = _terminal_peel_rescue(
        read, genome, cands, "+", 0.2, 10, 50, "chrI", genome["chrI"],
        baseline={"rescued": True, "edit_distance": 2.0, "query_bp": 10},
        peel_max_bp=100, peel_clean_anchor=10, accept_margin=0.0,
    )
    assert res is None, "peel must never override an existing baseline rescue"


def test_peel_skips_clean_terminal_alignment():
    """Gate A: a clean terminal alignment (no S/X/I/D near the 5' end) does not peel."""
    read = _make_read(start=1000, seq="C" * 30, cigar=((0, 30),))  # clean 30M, matches genome
    genome = {"chrI": "C" * 2000}
    cands = {("chrI", 980, 1000, "+")}
    res = _terminal_peel_rescue(
        read, genome, cands, "+", 0.2, 10, 50, "chrI", genome["chrI"],
        baseline={"rescued": False},
        peel_max_bp=100, peel_clean_anchor=10, accept_margin=0.0,
    )
    assert res is None, "clean terminal alignment must not trigger a peel"


# --- cost bounding (2026-05-24 production-scale hang regression) ------------

def test_peel_bounds_work_against_large_candidate_pool(monkeypatch):
    """A poorly-5'-anchored read against a ~5k-junction pool must NOT scan the
    full pool on every peel depth, nor run ~100 depths.

    Regression for the 2026-05-24 run-all hang: the peel looped up to
    peel_max_bp=100 depths, each scanning the full ~17k candidate set 3x in the
    body -> multi-hour stall + OOM on high-coverage loci. The fix (a) narrows the
    body's candidate set to the nearby subset, and (b) caps depth at the farthest
    reachable junction edge. Here only 2 junctions are near; the rest are 50kb
    away, so the body must see exactly the 2 nearby ones, and depths must be
    capped far below 100 (max_edge_dist=3 -> effective_max_peel=18).
    """
    import rectify.core.splice.splice_aware_5prime as m

    # 5' soft-clip + 200 mismatched body bases => Gate A fires, no clean anchor,
    # so the unbounded generator would otherwise run all the way to peel_max_bp.
    read = _make_read(start=1000, seq="G" * 5 + "A" * 200, cigar=((4, 5), (0, 200)))
    genome = {"chrI": "C" * 60000}

    near = {("chrI", 980, 1000, "+"), ("chrI", 970, 1003, "+")}  # edges 1000,1003
    far = {("chrI", 50000 + i, 50100 + i, "+") for i in range(5000)}
    cands = near | far

    calls = {"n": 0, "max_cand_seen": 0}
    real_body = m._rescue_3ss_truncation_body

    def _spy(read_, genome_, candidate_junctions, *a, **kw):
        calls["n"] += 1
        calls["max_cand_seen"] = max(calls["max_cand_seen"], len(candidate_junctions))
        return {"rescued": False}

    monkeypatch.setattr(m, "_rescue_3ss_truncation_body", _spy)

    m._terminal_peel_rescue(
        read, genome, cands, "+", 0.2, 10, 50, "chrI", genome["chrI"],
        baseline={"rescued": False},
        peel_max_bp=100, peel_clean_anchor=10, accept_margin=0.0,
    )

    assert calls["n"] <= 25, f"depth cap failed: {calls['n']} body calls (pre-fix ~100)"
    assert calls["max_cand_seen"] == len(near), (
        f"body saw {calls['max_cand_seen']} candidates; must be narrowed to the "
        f"{len(near)} nearby (pre-fix would be {len(cands)})"
    )


# --- screening (dry-run) mode ----------------------------------------------

def test_screening_does_not_change_production_output(tmp_path, monkeypatch):
    """With RECTIFY_PEEL_SCREEN set, a rescued baseline still yields None
    (fill-only production output is unchanged by screening)."""
    monkeypatch.setenv("RECTIFY_PEEL_SCREEN", str(tmp_path / "screen.tsv"))
    read = _make_read(start=1000, seq="A" * 10 + "C" * 20, cigar=((4, 10), (0, 20)))
    genome = {"chrI": "C" * 2000}
    cands = {("chrI", 980, 1000, "+")}
    res = _terminal_peel_rescue(
        read, genome, cands, "+", 0.2, 10, 50, "chrI", genome["chrI"],
        baseline={"rescued": True, "edit_distance": 2.0, "query_bp": 10},
        peel_max_bp=100, peel_clean_anchor=10, accept_margin=0.0,
    )
    assert res is None, "screening must not change the fill-only production output"


def test_log_peel_screen_writes_header_and_row(tmp_path):
    """The screening logger writes a header plus one TSV row with both junctions."""
    log = tmp_path / "screen.tsv"
    read = _make_read(start=900000, seq="C" * 30, cigar=((0, 30),), is_reverse=True)
    _log_peel_screen(
        str(log), read, "-", "chrI",
        baseline={"rescued": True, "rescued_junction": ("chrI", 900, 901193),
                  "five_prime_corrected": 901193, "edit_distance": 3.0},
        peel={"rescued": True, "rescued_junction": ("chrI", 900, 901189),
              "five_prime_corrected": 901189, "edit_distance": 1.0, "rescue_type": "softclip"},
        depth=20,
    )
    content = log.read_text()
    assert content.startswith("read_id\t"), "expected a header row"
    assert "900-901193" in content and "900-901189" in content
    assert content.count("\n") == 2, "header + exactly one data row"
