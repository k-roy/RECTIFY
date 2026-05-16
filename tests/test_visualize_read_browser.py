"""
Smoke tests for rectify.visualize.read_browser.

Covers:
- Module import + public function presence
- assign_rows: non-overlapping, overlapping, gap-aware behavior
- parse_junction_strings: standard format, NaN/empty handling, malformed
- draw_stacked_reads: smoke (LineCollection produced)
- plot_stacked_read_panel: smoke + empty DataFrame returns 0 + n_rows correctness
- Junction-aware drawing produces both exon and intron LineCollections

Author: Smoke-test suite
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import LineCollection

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize import read_browser as rb_mod
from rectify.visualize.read_browser import (
    assign_rows,
    draw_stacked_reads,
    parse_junction_strings,
    plot_stacked_read_panel,
)


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


def test_read_browser_imports():
    for fn in (
        "assign_rows",
        "parse_junction_strings",
        "draw_stacked_reads",
        "plot_stacked_read_panel",
    ):
        assert callable(getattr(rb_mod, fn))


# ---------------------------------------------------------------------------
# assign_rows
# ---------------------------------------------------------------------------


def test_assign_rows_non_overlapping():
    """Reads separated by gap → all fit in row 0."""
    starts = [0, 500, 1000]
    ends = [100, 600, 1100]
    rows = assign_rows(starts, ends, gap=50)
    assert (rows == 0).all()


def test_assign_rows_overlapping():
    """Three overlapping reads → 3 distinct rows."""
    starts = [0, 50, 100]
    ends = [200, 250, 300]
    rows = assign_rows(starts, ends, gap=50)
    assert set(rows.tolist()) == {0, 1, 2}


def test_assign_rows_gap_respected():
    """Reads exactly gap apart should still need different rows."""
    starts = [0, 110]
    ends = [100, 210]
    # 100 to 110 = 10 bp gap; need >50 bp gap → row 1
    rows = assign_rows(starts, ends, gap=50)
    assert rows[0] == 0
    assert rows[1] == 1


# ---------------------------------------------------------------------------
# parse_junction_strings
# ---------------------------------------------------------------------------


def test_parse_junctions_standard():
    inputs = ["100-200,300-400", "500-600"]
    out = parse_junction_strings(inputs)
    assert out == [[(100, 200), (300, 400)], [(500, 600)]]


def test_parse_junctions_empty_and_nan():
    """NaN, empty string, and None all produce empty lists."""
    inputs = ["", None, float("nan"), "valid-junction-1-100,200-300"]
    out = parse_junction_strings(inputs)
    assert out[0] == []
    assert out[1] == []
    assert out[2] == []
    # malformed parts gracefully skipped; "200-300" is the only valid pair-like
    assert (200, 300) in out[3]


def test_parse_junctions_malformed():
    """Bad ints are silently skipped, not raised."""
    out = parse_junction_strings(["abc-def,100-200"])
    assert out == [[(100, 200)]]


# ---------------------------------------------------------------------------
# draw_stacked_reads
# ---------------------------------------------------------------------------


def test_draw_stacked_reads_smoke():
    fig, ax = plt.subplots()
    starts = np.array([100, 500, 1000])
    ends = np.array([200, 600, 1100])
    rows = np.array([0, 0, 0])
    colors = ["red", "red", "blue"]  # 2 color groups

    draw_stacked_reads(ax, starts, ends, rows, colors)

    # 2 color groups → 2 LineCollection objects (exons only, no junctions)
    lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
    assert len(lcs) == 2


def test_draw_stacked_reads_empty():
    """Zero-length input: should return cleanly without adding artists."""
    fig, ax = plt.subplots()
    draw_stacked_reads(
        ax,
        np.array([], dtype=int),
        np.array([], dtype=int),
        np.array([], dtype=int),
        "red",
    )
    assert len(ax.collections) == 0


def test_draw_stacked_reads_with_junctions():
    """Junction list produces both exon and intron LineCollections."""
    fig, ax = plt.subplots()
    starts = np.array([100])
    ends = np.array([1000])
    rows = np.array([0])
    junctions = [[(300, 500)]]  # one intron

    draw_stacked_reads(ax, starts, ends, rows, "red", junction_lists=junctions)

    lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
    # one exon LC + one intron LC for color "red"
    assert len(lcs) == 2


# ---------------------------------------------------------------------------
# plot_stacked_read_panel
# ---------------------------------------------------------------------------


def test_plot_stacked_read_panel_smoke():
    df = pd.DataFrame(
        {
            "alignment_start": [100, 500, 1000],
            "alignment_end": [200, 600, 1100],
            "color": ["red", "red", "blue"],
            "junctions": ["", "", ""],
        }
    )
    fig, ax = plt.subplots()
    n_rows = plot_stacked_read_panel(ax, df)
    assert n_rows == 1
    # ylim was set
    assert ax.get_ylim() == (-0.5, 0.5)


def test_plot_stacked_read_panel_empty():
    df = pd.DataFrame(columns=["alignment_start", "alignment_end", "color"])
    fig, ax = plt.subplots()
    n_rows = plot_stacked_read_panel(ax, df)
    assert n_rows == 0


def test_plot_stacked_read_panel_overlapping_rows():
    """Three overlapping reads → n_rows == 3."""
    df = pd.DataFrame(
        {
            "alignment_start": [0, 50, 100],
            "alignment_end": [200, 250, 300],
            "color": ["red"] * 3,
        }
    )
    fig, ax = plt.subplots()
    n_rows = plot_stacked_read_panel(ax, df, junction_col=None, gap=10)
    assert n_rows == 3
