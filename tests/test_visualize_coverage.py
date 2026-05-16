"""
Smoke tests for rectify.visualize.coverage.

Covers:
- Module import + public function presence
- extract_coverage_from_bam: synthetic in-memory BAM
- extract_coverage_from_array: positions/counts → region arrays
- draw_coverage_track: renders without raising; produces fill artists
- draw_strand_coverage: stranded mirrored coverage rendering
- compare_coverage_tracks: multi-track overlay rendering
- Empty-input handling

All matplotlib output is rendered with the Agg backend (no display required).

Author: Smoke-test suite
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pysam
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize import coverage as coverage_mod
from rectify.visualize.coverage import (
    compare_coverage_tracks,
    draw_coverage_track,
    draw_strand_coverage,
    extract_coverage_from_array,
    extract_coverage_from_bam,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _close_figs():
    """Close all matplotlib figures after each test to keep memory low."""
    yield
    plt.close("all")


@pytest.fixture
def synthetic_bam(tmp_path):
    """Create a tiny indexed BAM with 5 reads on chrI.

    Reads:
        - read1: 100-200 forward
        - read2: 150-250 forward
        - read3: 300-400 reverse
        - read4: 105-205 forward (overlaps read1)
        - read5: 600-700 reverse
    """
    bam_path = tmp_path / "synthetic.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrI", "LN": 10000}],
    }
    reads_spec = [
        ("read1", 100, 100, False),
        ("read2", 150, 100, False),
        ("read3", 300, 100, True),
        ("read4", 105, 100, False),
        ("read5", 600, 100, True),
    ]
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for name, start, length, is_reverse in reads_spec:
            a = pysam.AlignedSegment(bam.header)
            a.query_name = name
            a.query_sequence = "A" * length
            a.flag = 16 if is_reverse else 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigar = [(0, length)]  # M
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            bam.write(a)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


# ---------------------------------------------------------------------------
# Import / presence
# ---------------------------------------------------------------------------


def test_coverage_imports():
    """Module imports and exposes the main public surface."""
    for fn_name in (
        "extract_coverage_from_bam",
        "extract_coverage_from_array",
        "draw_coverage_track",
        "draw_strand_coverage",
        "compare_coverage_tracks",
    ):
        assert callable(getattr(coverage_mod, fn_name)), f"{fn_name} not callable"


# ---------------------------------------------------------------------------
# extract_coverage_from_bam
# ---------------------------------------------------------------------------


def test_extract_coverage_from_bam_smoke(synthetic_bam):
    positions, depths = extract_coverage_from_bam(
        synthetic_bam, "chrI", 50, 300
    )
    assert positions.shape == depths.shape
    assert len(positions) == 250
    assert depths.sum() > 0, "should record some coverage"
    # In overlap region (150-200) we expect depth >= 2 (read1 + read2 + read4)
    overlap_idx = 150 - 50
    assert depths[overlap_idx] >= 2


def test_extract_coverage_from_bam_3prime_end(synthetic_bam):
    """use_3prime_end counts single base per read, not the read body."""
    positions, depths = extract_coverage_from_bam(
        synthetic_bam, "chrI", 50, 700, use_3prime_end=True
    )
    # 5 reads total → 5 counts (one per read's 3' end)
    assert depths.sum() == 5


def test_extract_coverage_from_bam_strand_filter(synthetic_bam):
    """Strand filter restricts to one strand only."""
    _, plus = extract_coverage_from_bam(
        synthetic_bam, "chrI", 50, 700, strand="+"
    )
    _, minus = extract_coverage_from_bam(
        synthetic_bam, "chrI", 50, 700, strand="-"
    )
    # Both should be non-empty (3 forward, 2 reverse reads)
    assert plus.sum() > 0
    assert minus.sum() > 0
    # They should not be the same
    assert not np.array_equal(plus, minus)


# ---------------------------------------------------------------------------
# extract_coverage_from_array
# ---------------------------------------------------------------------------


def test_extract_coverage_from_array_smoke():
    positions_in = np.array([100, 150, 200, 250])
    counts_in = np.array([5.0, 10.0, 3.0, 7.0])
    positions, depths = extract_coverage_from_array(
        positions_in, counts_in, 100, 251
    )
    assert len(positions) == 151
    # Spot-check the three positions that fell in range
    assert depths[0] == 5.0      # at 100
    assert depths[150 - 100] == 10.0  # at 150
    assert depths[250 - 100] == 7.0   # at 250


def test_extract_coverage_from_array_out_of_range():
    positions_in = np.array([5, 10000])
    counts_in = np.array([99.0, 99.0])
    _, depths = extract_coverage_from_array(positions_in, counts_in, 100, 200)
    assert depths.sum() == 0  # nothing in range


# ---------------------------------------------------------------------------
# draw_coverage_track
# ---------------------------------------------------------------------------


def test_draw_coverage_track_smoke():
    fig, ax = plt.subplots()
    positions = np.arange(100, 200)
    depths = np.linspace(0, 10, 100)

    draw_coverage_track(ax, positions, depths, 100, 200, label="test")

    # Has fill_between + line → at least one collection or line artist
    assert len(ax.collections) + len(ax.lines) >= 1
    assert ax.get_xlim() == (100, 200)


def test_draw_coverage_track_smoothing_and_log():
    fig, ax = plt.subplots()
    positions = np.arange(100, 200)
    depths = np.random.default_rng(0).integers(0, 100, size=100).astype(float)

    # Should not raise with smoothing + log scale combined
    draw_coverage_track(
        ax,
        positions,
        depths,
        100,
        200,
        smooth_window=5,
        log_scale=True,
        normalize=False,
    )
    assert "log10" in ax.get_ylabel().lower() or "log" in ax.get_ylabel().lower()


def test_draw_coverage_track_empty():
    """Zero-length input — should not raise."""
    fig, ax = plt.subplots()
    positions = np.array([], dtype=int)
    depths = np.array([], dtype=float)

    draw_coverage_track(ax, positions, depths, 0, 100)
    # No exception is the win condition


# ---------------------------------------------------------------------------
# draw_strand_coverage
# ---------------------------------------------------------------------------


def test_draw_strand_coverage_smoke():
    fig, ax = plt.subplots()
    pos = np.arange(100, 200)
    plus = np.linspace(0, 5, 100)
    minus = np.linspace(0, 3, 100)

    draw_strand_coverage(ax, pos, plus, pos, minus, 100, 200)

    # Should have axhline + 2 fills + 2 plot lines
    assert len(ax.collections) >= 2
    assert len(ax.lines) >= 1


# ---------------------------------------------------------------------------
# compare_coverage_tracks
# ---------------------------------------------------------------------------


def test_compare_coverage_tracks_smoke():
    fig, ax = plt.subplots()
    pos = np.arange(100, 200)
    cov_dict = {
        "WT": (pos, np.linspace(0, 10, 100)),
        "mut": (pos, np.linspace(0, 5, 100)),
    }
    compare_coverage_tracks(ax, cov_dict, 100, 200)

    # Two plot lines (one per condition)
    assert len(ax.lines) >= 2
    # Legend should exist
    leg = ax.get_legend()
    assert leg is not None
