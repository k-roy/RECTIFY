"""Test that build_junction_pool uses ProcessPoolExecutor for multiple BAMs
and falls back to the sequential path for 0 or 1 BAMs.

The result must be identical regardless of which code path executes.
"""

from pathlib import Path

import pytest

from rectify.core.splice.junction_scoring import (
    _build_junction_index,
    _candidates_near,
    build_junction_pool,
    collect_junctions_from_bam,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_ALIGNERS_DIR = (
    Path(__file__).parent.parent / "rectify" / "data" / "validation" / "aligners"
)


def _get_aligner_bam(name: str) -> Path:
    p = _ALIGNERS_DIR / f"validation_reads.{name}.bam"
    return p


@pytest.fixture
def bam_minimap2():
    p = _get_aligner_bam("minimap2")
    if not p.exists():
        pytest.skip(f"minimap2 fixture BAM not found: {p}")
    return str(p)


@pytest.fixture
def bam_gapmm2():
    p = _get_aligner_bam("gapmm2")
    if not p.exists():
        pytest.skip(f"gapmm2 fixture BAM not found: {p}")
    return str(p)


@pytest.fixture
def bam_desalt():
    p = _get_aligner_bam("deSALT")
    if not p.exists():
        pytest.skip(f"deSALT fixture BAM not found: {p}")
    return str(p)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_empty_bam_list_returns_annotated_only():
    """With no BAMs, the pool equals the annotated set."""
    annotated = {("chrI", 100, 200), ("chrII", 500, 600)}
    pool, annot_3 = build_junction_pool([], annotated)
    # Pool should contain the annotated junctions (normalised to 3-tuples)
    for j in annotated:
        assert j in pool, f"Annotated junction {j} missing from pool"
    assert len(pool) == len(annotated), (
        f"Pool size {len(pool)} != annotated size {len(annotated)}"
    )


def test_single_bam_sequential_path(bam_minimap2):
    """With 1 BAM the sequential path fires and returns the correct union."""
    annotated: set = set()
    # Disable the anchor filter here: this test exercises the sequential/union
    # machinery, not the pool-quality filter (covered by test_junction_anchor_filter).
    pool, annot_3 = build_junction_pool([bam_minimap2], annotated, min_anchor_overhang=0)
    # Must contain at least what collect_junctions_from_bam returns directly.
    direct = collect_junctions_from_bam(bam_minimap2)
    for j in direct:
        assert j in pool, f"Junction {j} from direct scan missing from pool"


def test_two_bams_parallel_path_matches_sequential(bam_minimap2, bam_gapmm2):
    """With 2 BAMs the parallel path fires.

    The result must be the set-union of both BAMs' junctions.
    We verify: parallel result == sequential result.
    """
    annotated: set = set()

    # Sequential reference: call build_junction_pool on each BAM individually
    # and union manually.
    j1 = collect_junctions_from_bam(bam_minimap2)
    j2 = collect_junctions_from_bam(bam_gapmm2)
    expected = j1 | j2

    # Parallel path: 2 BAMs triggers ProcessPoolExecutor. Anchor filter off —
    # this verifies the parallel union machinery, not pool quality.
    pool_parallel, _ = build_junction_pool(
        [bam_minimap2, bam_gapmm2], annotated, min_anchor_overhang=0
    )

    assert pool_parallel == expected, (
        f"Parallel path returned {len(pool_parallel)} junctions, "
        f"expected {len(expected)}.\n"
        f"Extra: {pool_parallel - expected}\n"
        f"Missing: {expected - pool_parallel}"
    )


def test_three_bams_parallel_path_matches_sequential(bam_minimap2, bam_gapmm2, bam_desalt):
    """With 3 BAMs the parallel path fires and the result is the union of all three."""
    annotated: set = set()

    j1 = collect_junctions_from_bam(bam_minimap2)
    j2 = collect_junctions_from_bam(bam_gapmm2)
    j3 = collect_junctions_from_bam(bam_desalt)
    expected = j1 | j2 | j3

    pool_parallel, _ = build_junction_pool(
        [bam_minimap2, bam_gapmm2, bam_desalt], annotated, min_anchor_overhang=0
    )

    assert pool_parallel == expected, (
        f"Parallel path (3 BAMs) returned {len(pool_parallel)} junctions, "
        f"expected {len(expected)}.\n"
        f"Extra: {pool_parallel - expected}\n"
        f"Missing: {expected - pool_parallel}"
    )


def test_annotated_junctions_included_with_parallel(bam_minimap2, bam_gapmm2):
    """Annotated junctions must appear in the pool even when parallel path fires."""
    annotated = {("chrI", 1000, 2000), ("chrV", 3000, 4000)}
    pool, annot_3 = build_junction_pool([bam_minimap2, bam_gapmm2], annotated)

    for j in annotated:
        assert j in pool, f"Annotated junction {j} missing from parallel pool"
    # annot_3 should be the normalised annotated subset
    for j in annotated:
        assert j in annot_3, f"Annotated junction {j} missing from annot_3"


def test_min_observed_support_filters_observed_junctions_but_keeps_annotated(
    bam_minimap2, bam_gapmm2,
):
    """Observed junctions below support threshold are filtered; annotations stay."""
    j1 = collect_junctions_from_bam(bam_minimap2)
    j2 = collect_junctions_from_bam(bam_gapmm2)
    shared = j1 & j2
    only_one = (j1 ^ j2)
    if not shared or not only_one:
        pytest.skip("fixture BAMs do not have both shared and aligner-unique junctions")

    annotated_unique = next(iter(only_one))
    pool, annot_3 = build_junction_pool(
        [bam_minimap2, bam_gapmm2],
        {annotated_unique},
        min_observed_support=2,
        min_anchor_overhang=0,  # isolate the support filter from the anchor filter
    )

    assert annotated_unique in pool
    assert annotated_unique in annot_3
    assert shared <= pool
    assert not ((only_one - {annotated_unique}) & pool)


def test_candidates_near_uses_bounded_start_window_and_max_size():
    """Candidate lookup remains endpoint-bounded and can cap intron length."""
    junctions = {
        ("chrI", 100, 300),
        ("chrI", 900, 10950),   # near start/end but too large for yeast cap
        ("chrI", 1000, 1200),   # expected hit
        ("chrI", 1025, 1220),   # expected hit
        ("chrI", 4000, 4200),
        ("chrII", 1000, 1200),
    }
    idx = _build_junction_index(junctions)

    assert _candidates_near(
        idx,
        "chrI",
        1000,
        1200,
        radius=5000,
        start_radius=50,
        end_radius=50,
        max_junction_size=10000,
    ) == [(1000, 1200), (1025, 1220)]

    assert _candidates_near(
        idx,
        "chrI",
        900,
        10950,
        radius=5000,
        start_radius=50,
        end_radius=50,
        max_junction_size=10000,
    ) == []


def test_build_junction_pool_max_size_filters_observed_but_keeps_annotated(
    bam_minimap2,
):
    """The max-size cap filters BAM-derived calls only; annotations stay."""
    observed = collect_junctions_from_bam(bam_minimap2)
    if not observed:
        pytest.skip("fixture BAM has no observed junctions")

    annotated_long = ("chrI", 10, 20000)
    pool, annot_3 = build_junction_pool(
        [bam_minimap2],
        {annotated_long},
        max_junction_size=1,
    )

    assert annotated_long in pool
    assert annotated_long in annot_3
    assert not (observed & (pool - {annotated_long}))
