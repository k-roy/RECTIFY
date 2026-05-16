"""
Smoke tests for rectify.visualize.multi_track.

Covers:
- Module import + MultiTrackFigure class presence
- Empty figure render() raises clear error
- Single gene track render
- Coverage + gene track combo render
- Custom track render
- create_gene_browser convenience function
- save() round-trip via tmp_path

Author: Smoke-test suite
"""

import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize import multi_track as mt_mod
from rectify.visualize.multi_track import MultiTrackFigure, create_gene_browser


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


def _make_cds(gene_id, start, end, strand="+", common_name=None):
    return SimpleNamespace(
        gene_id=gene_id,
        feature_type="CDS",
        start=start,
        end=end,
        strand=strand,
        common_name=common_name,
    )


@pytest.fixture
def gff_features():
    return {
        "chrI": [
            _make_cds("YAL001C", 100, 1000, "-", "TFC3"),
            _make_cds("YAL002W", 2000, 3500, "+", "VPS8"),
        ],
        "_gene_name_lookup": {"YAL001C": "TFC3", "YAL002W": "VPS8"},
    }


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


def test_multi_track_imports():
    assert hasattr(mt_mod, "MultiTrackFigure")
    assert hasattr(mt_mod, "create_gene_browser")
    assert callable(mt_mod.create_gene_browser)


# ---------------------------------------------------------------------------
# MultiTrackFigure
# ---------------------------------------------------------------------------


def test_multi_track_no_tracks_raises():
    mtf = MultiTrackFigure()
    with pytest.raises(ValueError, match="No tracks"):
        mtf.render(0, 1000)


def test_multi_track_gene_only_render(gff_features):
    mtf = MultiTrackFigure().add_gene_track(gff_features, highlight_gene="VPS8")
    fig = mtf.render(0, 5000)
    assert isinstance(fig, plt.Figure)
    # At least one axes
    assert len(fig.axes) >= 1


def test_multi_track_gene_plus_coverage(gff_features):
    pos = np.arange(0, 5000)
    depths = np.linspace(0, 100, 5000)

    mtf = (
        MultiTrackFigure(figsize=(8, 4))
        .add_gene_track(gff_features, highlight_gene="VPS8")
        .add_coverage_track("WT", (pos, depths), color="#0072B2")
    )
    fig = mtf.render(0, 5000)
    assert len(fig.axes) == 2


def test_multi_track_coverage_dict_form(gff_features):
    """Coverage data can also be supplied as a dict with positions/depths keys."""
    pos = np.arange(0, 5000)
    depths = np.linspace(0, 50, 5000)
    mtf = (
        MultiTrackFigure()
        .add_gene_track(gff_features, highlight_gene="VPS8")
        .add_coverage_track("WT", {"positions": pos, "depths": depths})
    )
    fig = mtf.render(0, 5000)
    assert len(fig.axes) == 2


def test_multi_track_custom_track():
    """A custom track simply gets the axes and renders into it."""
    sentinel = []

    def my_drawer(ax, data, region_start, region_end, **kwargs):
        sentinel.append((region_start, region_end))
        ax.axhline(y=0.5, color="purple")

    mtf = MultiTrackFigure().add_custom_track(my_drawer, data={"foo": 1})
    fig = mtf.render(100, 1000)
    assert sentinel == [(100, 1000)]
    assert len(fig.axes) == 1


def test_multi_track_save(tmp_path, gff_features):
    pos = np.arange(0, 5000)
    depths = np.ones(5000) * 5

    mtf = (
        MultiTrackFigure(figsize=(6, 3))
        .add_gene_track(gff_features, highlight_gene="VPS8")
        .add_coverage_track("WT", (pos, depths))
    )
    out_path = tmp_path / "browser.png"
    mtf.save(out_path, 0, 5000)
    assert out_path.exists()
    assert out_path.stat().st_size > 100  # non-empty


# ---------------------------------------------------------------------------
# create_gene_browser
# ---------------------------------------------------------------------------


def test_create_gene_browser_smoke(gff_features):
    pos = np.arange(0, 5000)
    depths = np.linspace(0, 30, 5000)

    fig = create_gene_browser(
        gff_features,
        gene_name="VPS8",
        region_start=0,
        region_end=5000,
        coverage_tracks={"WT": (pos, depths)},
    )
    assert isinstance(fig, plt.Figure)
    # 1 gene track + 1 coverage track
    assert len(fig.axes) == 2
