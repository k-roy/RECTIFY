"""
Smoke tests for rectify.visualize.vep_panels.

Note: vep_panels is currently commented out of the visualize package's
__init__.py top-level exports. We import it via the submodule path directly.

Covers:
- Module import + public function presence
- extract_position_from_variants: standard, malformed, NaN/empty
- draw_evo2_panel: smoke (renders scatter on synthetic DataFrame) + empty
- draw_esm1v_panel: smoke + empty
- draw_shorkie_pergene_panel: smoke + empty
- draw_yorzoi_pergene_panel: smoke (uses log2 transform)
- draw_vep_panel_styled: smoke (the generic variant)
- Missing-column branches surface the "No ... data" fallback text

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

from rectify.visualize import vep_panels as vp_mod
from rectify.visualize.vep_panels import (
    draw_esm1v_panel,
    draw_evo2_panel,
    draw_shorkie_pergene_panel,
    draw_vep_panel_styled,
    draw_yorzoi_pergene_panel,
    extract_position_from_variants,
)


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------


def test_vep_panels_imports():
    for fn in (
        "extract_position_from_variants",
        "draw_evo2_panel",
        "draw_esm1v_panel",
        "draw_shorkie_pergene_panel",
        "draw_yorzoi_pergene_panel",
        "draw_vep_panel_styled",
    ):
        assert callable(getattr(vp_mod, fn))


# ---------------------------------------------------------------------------
# extract_position_from_variants
# ---------------------------------------------------------------------------


def test_extract_position_standard():
    assert extract_position_from_variants("530123_A_G_missense") == 530123


def test_extract_position_first_of_many():
    assert (
        extract_position_from_variants(
            "530123_A_G_missense,530125_C_T_synonymous"
        )
        == 530123
    )


def test_extract_position_nan():
    assert extract_position_from_variants(float("nan")) is None


def test_extract_position_empty():
    assert extract_position_from_variants("") is None


def test_extract_position_malformed():
    """Non-numeric prefix → None, not exception."""
    assert extract_position_from_variants("foo_bar") is None


# ---------------------------------------------------------------------------
# draw_evo2_panel
# ---------------------------------------------------------------------------


def test_draw_evo2_panel_smoke():
    df = pd.DataFrame(
        {
            "variant_coord": [100, 200, 300],
            "variant_score": [-1.0, 0.0, 1.0],
            "codon_variant_type": ["missense", "synonymous", "nonsense"],
            "is_sub_variant": [False, False, False],
            "is_filtered": [False, False, False],
        }
    )
    fig, ax = plt.subplots()
    draw_evo2_panel(ax, df, 0, 500)

    # 3 scatter PathCollections produced
    assert len(ax.collections) >= 3
    assert ax.get_xlim() == (0, 500)


def test_draw_evo2_panel_empty_df():
    fig, ax = plt.subplots()
    draw_evo2_panel(ax, pd.DataFrame(), 0, 500)
    texts = [t.get_text() for t in ax.texts]
    assert any("No Evo2" in t for t in texts)


def test_draw_evo2_panel_none():
    fig, ax = plt.subplots()
    draw_evo2_panel(ax, None, 0, 500)
    texts = [t.get_text() for t in ax.texts]
    assert any("No Evo2" in t for t in texts)


def test_draw_evo2_panel_no_score_column():
    """Missing score column → fallback 'No score column' message."""
    df = pd.DataFrame({"variant_coord": [100], "codon_variant_type": ["missense"]})
    fig, ax = plt.subplots()
    draw_evo2_panel(ax, df, 0, 500)
    texts = [t.get_text() for t in ax.texts]
    assert any("No score" in t for t in texts)


def test_draw_evo2_panel_filters_sub_variants():
    """is_sub_variant=True rows are excluded."""
    df = pd.DataFrame(
        {
            "variant_coord": [100, 200],
            "variant_score": [-1.0, 1.0],
            "codon_variant_type": ["missense", "missense"],
            "is_sub_variant": [True, True],   # both excluded
            "is_filtered": [False, False],
        }
    )
    fig, ax = plt.subplots()
    draw_evo2_panel(ax, df, 0, 500)
    texts = [t.get_text() for t in ax.texts]
    assert any("No complete variants" in t for t in texts)


# ---------------------------------------------------------------------------
# draw_esm1v_panel
# ---------------------------------------------------------------------------


def test_draw_esm1v_panel_smoke():
    df = pd.DataFrame(
        {
            "POS": [100, 200, 300],
            "ESM1v_Score": [-0.5, 0.2, 1.1],
            "codon_variant_type": ["missense"] * 3,
        }
    )
    fig, ax = plt.subplots()
    draw_esm1v_panel(ax, df, 0, 500)
    assert len(ax.collections) >= 3


def test_draw_esm1v_panel_empty():
    fig, ax = plt.subplots()
    draw_esm1v_panel(ax, pd.DataFrame(), 0, 500)
    texts = [t.get_text() for t in ax.texts]
    assert any("No ESM1v" in t for t in texts)


# ---------------------------------------------------------------------------
# draw_shorkie_pergene_panel
# ---------------------------------------------------------------------------


def test_draw_shorkie_pergene_panel_smoke():
    df = pd.DataFrame(
        {
            "variant_description": ["100_A_G_missense", "200_C_T_missense"],
            "gene_systematic": ["YAL001C", "YAL001C"],
            "gene_common": ["GENE1", "GENE1"],
            "gene_type": ["target", "target"],
            "shorkie_logSED": [-1.0, 0.5],
        }
    )
    fig, ax = plt.subplots()
    draw_shorkie_pergene_panel(ax, df, 0, 500, gene_name="GENE1")
    assert len(ax.collections) >= 1


def test_draw_shorkie_pergene_panel_empty():
    fig, ax = plt.subplots()
    draw_shorkie_pergene_panel(ax, pd.DataFrame(), 0, 500, gene_name="X")
    texts = [t.get_text() for t in ax.texts]
    assert any("No Shorkie" in t for t in texts)


# ---------------------------------------------------------------------------
# draw_yorzoi_pergene_panel
# ---------------------------------------------------------------------------


def test_draw_yorzoi_pergene_panel_smoke():
    df = pd.DataFrame(
        {
            "variant_description": ["100_A_G_missense", "200_C_T_missense"],
            "gene_systematic": ["YAL001C", "YAL001C"],
            "gene_common": ["GENE1", "GENE1"],
            "gene_type": ["target", "target"],
            "yorzoi_ratio": [0.5, 2.0],
        }
    )
    fig, ax = plt.subplots()
    draw_yorzoi_pergene_panel(ax, df, 0, 500, gene_name="GENE1")
    assert len(ax.collections) >= 1


def test_draw_yorzoi_pergene_panel_missing_ratio_column():
    df = pd.DataFrame(
        {"variant_description": ["100_A_G"], "gene_type": ["target"]}
    )
    fig, ax = plt.subplots()
    draw_yorzoi_pergene_panel(ax, df, 0, 500, gene_name="X")
    texts = [t.get_text() for t in ax.texts]
    assert any("No Yorzoi ratio" in t for t in texts)


# ---------------------------------------------------------------------------
# draw_vep_panel_styled (generic)
# ---------------------------------------------------------------------------


def test_draw_vep_panel_styled_smoke():
    df = pd.DataFrame(
        {
            "POS": [100, 200, 300],
            "my_score": [-1.0, 0.0, 1.0],
            "codon_variant_type": ["missense", "synonymous", "missense"],
        }
    )
    fig, ax = plt.subplots()
    draw_vep_panel_styled(
        ax, df, score_col="my_score", ylabel="Custom",
        region_start=0, region_end=500,
    )
    assert len(ax.collections) >= 3


def test_draw_vep_panel_styled_missing_score_col():
    df = pd.DataFrame({"POS": [100], "codon_variant_type": ["missense"]})
    fig, ax = plt.subplots()
    draw_vep_panel_styled(
        ax, df, score_col="nonexistent", ylabel="X",
        region_start=0, region_end=500,
    )
    texts = [t.get_text() for t in ax.texts]
    # When score_col is missing entirely, the function hits the
    # "No X data" fallback (it short-circuits before _filter_complete_variants).
    assert any("No X data" in t for t in texts)


def test_draw_vep_panel_styled_log2():
    df = pd.DataFrame(
        {
            "POS": [100, 200],
            "ratio": [2.0, 4.0],
            "codon_variant_type": ["missense", "missense"],
        }
    )
    fig, ax = plt.subplots()
    draw_vep_panel_styled(
        ax, df, score_col="ratio", ylabel="log2",
        region_start=0, region_end=500, use_log2=True,
    )
    # 2 scatter points (both positive ratios survive log2)
    assert len(ax.collections) >= 2
