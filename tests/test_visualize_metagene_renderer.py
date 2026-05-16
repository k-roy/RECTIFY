"""
Smoke tests for the metagene renderer methods of rectify.visualize.metagene.

The existing tests/test_metagene.py focuses on signal-extraction correctness;
this file covers the *plotting* side (plot_profile, plot_multiple_profiles)
that have no coverage today.

Covers:
- Module-level import of MetagenePipeline.plot_profile / plot_multiple_profiles
- plot_profile renders without raising on a real compute_profile() result
- plot_multiple_profiles renders both conditions on shared axes
- Empty-result rendering (n_loci=0) does not raise

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

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize import metagene as mg_mod
from rectify.visualize.metagene import (
    MetageneConfig,
    MetagenePipeline,
    PositionIndex,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


@pytest.fixture
def synthetic_pipeline_result():
    """Build a real compute_profile() result from synthetic in-memory data."""
    # Loci: 5 features, each 200 bp long on chrI '+' strand
    loci_rows = []
    sig_rows = []
    rng = np.random.default_rng(0)
    for i in range(5):
        start = 1000 + i * 5000
        end = start + 200
        loci_rows.append({"chrom": "chrI", "strand": "+", "start": start, "end": end})
        # Sprinkle signal across body
        for pos in range(start, end, 5):
            sig_rows.append(
                {"chrom": "chrI", "strand": "+", "position": pos, "count": int(rng.integers(1, 5))}
            )

    loci = pd.DataFrame(loci_rows)
    index = PositionIndex(pd.DataFrame(sig_rows), position_col="position")

    cfg = MetageneConfig(window_upstream=50, window_downstream=50, scaled_body_width=100)
    pipeline = MetagenePipeline(cfg)
    result = pipeline.compute_profile(loci, index)
    return pipeline, result


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


def test_metagene_renderer_imports():
    # Methods on the pipeline class
    assert hasattr(MetagenePipeline, "plot_profile")
    assert hasattr(MetagenePipeline, "plot_multiple_profiles")
    assert callable(MetagenePipeline.plot_profile)


# ---------------------------------------------------------------------------
# plot_profile
# ---------------------------------------------------------------------------


def test_plot_profile_smoke(synthetic_pipeline_result):
    pipeline, result = synthetic_pipeline_result
    assert result["n_loci"] > 0, "fixture should produce a non-empty profile"

    fig, ax = plt.subplots()
    pipeline.plot_profile(ax, result, color="#0072B2", label="WT")

    # Should have at least one line drawn + SEM fill_between collection
    assert len(ax.lines) >= 1
    # Feature boundary axvlines also count as lines → at least 3 lines total
    # (profile line + 2 axvlines)
    assert len(ax.lines) >= 3
    assert ax.get_xlabel().strip() != ""
    assert ax.get_ylabel().strip() != ""


def test_plot_profile_without_sem_and_bounds(synthetic_pipeline_result):
    pipeline, result = synthetic_pipeline_result
    fig, ax = plt.subplots()
    pipeline.plot_profile(
        ax, result, color="#D55E00",
        show_sem=False,
        show_feature_bounds=False,
    )
    # Only the profile line remains
    assert len(ax.lines) == 1
    # No fill_between collection
    assert len(ax.collections) == 0


# ---------------------------------------------------------------------------
# plot_multiple_profiles
# ---------------------------------------------------------------------------


def test_plot_multiple_profiles_smoke(synthetic_pipeline_result):
    pipeline, result = synthetic_pipeline_result

    results = {"WT": result, "mutant": result}  # reuse same profile
    fig, ax = plt.subplots()
    pipeline.plot_multiple_profiles(ax, results)

    # Two profile lines + first call's 2 axvlines
    assert len(ax.lines) >= 2
    leg = ax.get_legend()
    assert leg is not None, "legend should be drawn by default"
    # Both labels should appear (with n=K suffix)
    labels = [t.get_text() for t in leg.get_texts()]
    assert any("WT" in l for l in labels)
    assert any("mutant" in l for l in labels)


# ---------------------------------------------------------------------------
# Empty input
# ---------------------------------------------------------------------------


def test_plot_profile_empty_result():
    """plot_profile must not raise on an n_loci=0 result."""
    cfg = MetageneConfig(window_upstream=10, window_downstream=10, scaled_body_width=20)
    pipeline = MetagenePipeline(cfg)
    # Hand-construct an "empty" result matching the function's return shape
    profile_length = cfg.window_upstream + cfg.scaled_body_width + cfg.window_downstream
    empty_result = {
        "profile": np.zeros(profile_length),
        "profile_matrix": np.zeros((0, profile_length)),
        "sem": np.zeros(profile_length),
        "n_loci": 0,
        "x_positions": np.arange(profile_length) - cfg.window_upstream,
        "body_width": cfg.scaled_body_width,
    }
    fig, ax = plt.subplots()
    pipeline.plot_profile(ax, empty_result, label="empty")
    # No exception raised is the win
    assert len(ax.lines) >= 1
