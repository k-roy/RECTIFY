"""Tests for the dual-site (interval-tree) junction-pool index.

The 3'SS rescue fetches pool junctions near a read's 5' boundary. The index was
previously a per-chrom list sorted by ``intron_start`` queried with a fixed
radius, which can only surface a junction whose ``intron_start`` is within the
radius. For a plus-strand read the relevant (near) splice site is the ACCEPTOR
(``intron_end``); for a human intron the donor can sit >100 kb away, outside any
sane radius — so the old shape was correct for yeast only by accident.

The interval-tree index keys on the whole intron, so a junction is surfaced
whenever its near site sits by the read OR the read's 5' end falls inside the
intron (Case-4 snap), independent of intron length. These tests pin that
property.
"""

from rectify.core.bam.bam_processor import (
    _build_pool_chrom_index,
    _POOL_FETCH_HALF_WINDOW,
)


def _fetch(tree, pos, five_clip=0):
    """Mirror the correct_read_3prime fetch: intervals within the window of pos."""
    w = _POOL_FETCH_HALF_WINDOW + five_clip
    return {iv.data for iv in tree.overlap(pos - w, pos + w + 1)}


def test_build_returns_interval_tree_per_chrom():
    from intervaltree import IntervalTree

    idx = _build_pool_chrom_index({("chrI", 100, 200), ("chrII", 300, 400)})
    assert set(idx) == {"chrI", "chrII"}
    assert all(isinstance(t, IntervalTree) for t in idx.values())


def test_empty_or_none_pool_returns_none():
    assert _build_pool_chrom_index(None) is None
    assert _build_pool_chrom_index(set()) is None


def test_degenerate_intervals_skipped():
    # intron_end <= intron_start cannot be an IntervalTree interval; drop it
    # rather than crash.
    idx = _build_pool_chrom_index({("chrI", 500, 500), ("chrI", 600, 400), ("chrI", 100, 200)})
    hits = {iv.data for iv in idx["chrI"]}
    assert hits == {(100, 200)}


def test_either_site_surfaces_far_site_human_intron():
    # Two introns share an acceptor near 5200; one is short (yeast-scale), the
    # other spans 120 kb (human-scale donor far upstream). A read whose exon-2 5'
    # boundary aligns at ~5202 must surface BOTH — the old intron_start radius
    # would have missed the 120 kb one.
    short = (5000, 5200)
    human = (5200 - 120_000, 5200)
    distant = (9000, 9300)
    tree = _build_pool_chrom_index({("chrI", *short), ("chrI", *human), ("chrI", *distant)})["chrI"]
    near = _fetch(tree, 5202)
    assert short in near
    assert human in near, "far-site (human-scale) intron missed — the bug this fixes"
    assert distant not in near, "distant intron must not be pulled in"


def test_containment_surfaces_case4_deep_intron():
    # A read whose 5' end mis-aligns DEEP inside a long intron (neither splice
    # site near it) must still surface that intron for the Case-4 intronic snap.
    human = (0, 120_000)
    tree = _build_pool_chrom_index({("chrI", *human)})["chrI"]
    deep = _fetch(tree, 60_000)
    assert human in deep


def test_five_clip_extends_reach():
    # The baseline sequence rescue reaches a junction edge up to
    # junction_proximity_bp + five_clip away, so the fetch window must grow with
    # the read's 5' soft-clip length.
    edge = 5200
    j = (5000, edge)
    tree = _build_pool_chrom_index({("chrI", *j)})["chrI"]
    far_pos = edge + _POOL_FETCH_HALF_WINDOW + 200  # beyond the no-clip window
    assert j not in _fetch(tree, far_pos, five_clip=0)
    assert j in _fetch(tree, far_pos, five_clip=300)
