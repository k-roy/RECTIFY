"""
Smoke tests for rectify.visualize.ridge.

Covers:
- Module import + public function presence
- aggregate_replicate_profiles: single + multiple replicates + smoothing
- plot_ridge_panel: smoke render with/without SEM band
- create_overlapping_ridge_figure: panel count, optional gene-track axes
- Empty / error paths

Author: Smoke-test suite
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize import ridge as ridge_mod
from rectify.visualize.ridge import (
    aggregate_replicate_profiles,
    create_overlapping_ridge_figure,
    plot_ridge_panel,
)


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


def test_ridge_imports():
    for fn in (
        "aggregate_replicate_profiles",
        "plot_ridge_panel",
        "create_overlapping_ridge_figure",
    ):
        assert callable(getattr(ridge_mod, fn))


# ---------------------------------------------------------------------------
# aggregate_replicate_profiles
# ---------------------------------------------------------------------------


def test_aggregate_replicate_profiles_smoke():
    rng = np.random.default_rng(0)
    profiles = [rng.normal(5, 1, size=101) for _ in range(3)]
    mean, sem = aggregate_replicate_profiles(profiles)
    assert mean.shape == (101,)
    assert sem.shape == (101,)
    # SEM should be > 0 for multiple replicates
    assert (sem > 0).all()


def test_aggregate_replicate_profiles_single():
    """Single replicate → SEM is zeros."""
    mean, sem = aggregate_replicate_profiles([np.arange(10, dtype=float)])
    assert np.allclose(sem, 0)
    assert np.allclose(mean, np.arange(10))


def test_aggregate_replicate_profiles_smoothing():
    """smooth_sigma applies gaussian filter before aggregation."""
    profiles = [np.zeros(101) for _ in range(3)]
    for p in profiles:
        p[50] = 100.0  # impulse at center
    mean, _ = aggregate_replicate_profiles(profiles, smooth_sigma=3.0)
    # Impulse should be smoothed → peak at 50 but lower amplitude, side bins > 0
    assert mean[50] < 100.0
    assert mean[48] > 0.0
    assert mean[52] > 0.0


def test_aggregate_replicate_profiles_empty_raises():
    with pytest.raises(ValueError, match="empty"):
        aggregate_replicate_profiles([])


# ---------------------------------------------------------------------------
# plot_ridge_panel
# ---------------------------------------------------------------------------


def test_plot_ridge_panel_smoke():
    fig, ax = plt.subplots()
    x = np.arange(-50, 51)
    mean = np.exp(-(x ** 2) / (2 * 10 ** 2))  # bell curve

    plot_ridge_panel(ax, x, mean, color="#0072B2", label="WT")
    # fill_between collection + plot line + axhline = at least 1 collection + 2 lines
    assert len(ax.collections) >= 1
    assert len(ax.lines) >= 1
    # Y-axis is hidden (no yticks)
    assert len(ax.get_yticks()) == 0
    # ylim bottom should be clipped at 0
    assert ax.get_ylim()[0] == 0


def test_plot_ridge_panel_with_sem():
    fig, ax = plt.subplots()
    x = np.arange(-50, 51)
    mean = np.exp(-(x ** 2) / (2 * 10 ** 2))
    sem = np.ones_like(mean) * 0.05
    plot_ridge_panel(ax, x, mean, sem=sem, color="red")
    # 2 fill_betweens (main + sem) → at least 2 collections
    assert len(ax.collections) >= 2


def test_plot_ridge_panel_label_drawn():
    fig, ax = plt.subplots()
    x = np.arange(-10, 11)
    mean = np.ones_like(x, dtype=float)
    plot_ridge_panel(ax, x, mean, label="my-condition")
    texts = [t.get_text() for t in ax.texts]
    assert "my-condition" in texts


# ---------------------------------------------------------------------------
# create_overlapping_ridge_figure
# ---------------------------------------------------------------------------


def test_create_overlapping_ridge_figure_smoke():
    fig, axes, ax_track = create_overlapping_ridge_figure(n_panels=3)
    assert len(axes) == 3
    assert ax_track is None
    # All axes belong to the figure
    for ax in axes:
        assert ax.figure is fig


def test_create_overlapping_ridge_figure_with_gene_track():
    fig, axes, ax_track = create_overlapping_ridge_figure(
        n_panels=2, add_gene_track=True
    )
    assert len(axes) == 2
    assert ax_track is not None
    assert ax_track.figure is fig


def test_create_overlapping_ridge_figure_single_panel():
    """n_panels=1 still returns a list, not a bare axes."""
    fig, axes, _ = create_overlapping_ridge_figure(n_panels=1)
    assert isinstance(axes, list)
    assert len(axes) == 1
