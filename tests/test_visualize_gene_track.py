"""
Smoke tests for rectify.visualize.gene_track.

Covers:
- Module import + public function presence
- assign_feature_levels: non-overlapping vs overlapping layouts
- draw_gene_arrow: standalone arrow draw (plus and minus strand)
- draw_gene_track: full track with synthetic GeneFeature-shaped objects
- get_genes_in_region: dict lookup
- Edge cases: no features in region, gene name absent

Synthetic GeneFeature: a SimpleNamespace with the attributes the function
accesses (gene_id, feature_type, start, end, strand, common_name).

Author: Smoke-test suite
"""

import sys
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize import gene_track as gt_mod
from rectify.visualize.gene_track import (
    assign_feature_levels,
    draw_gene_arrow,
    draw_gene_track,
    get_genes_in_region,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


def _make_cds(gene_id, start, end, strand="+", common_name=None):
    """Build a CDS-like feature object the gene_track functions can introspect."""
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
    """A small synthetic GFF-features dict mirroring the real schema."""
    return {
        "chrI": [
            _make_cds("YAL001C", 100, 1000, "-", "TFC3"),
            _make_cds("YAL002W", 2000, 3500, "+", "VPS8"),
            _make_cds("YAL003W", 4000, 5500, "+", "EFB1"),
        ],
        "_gene_name_lookup": {
            "YAL001C": "TFC3",
            "YAL002W": "VPS8",
            "YAL003W": "EFB1",
        },
    }


# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------


def test_gene_track_imports():
    for fn_name in (
        "assign_feature_levels",
        "draw_gene_track",
        "get_genes_in_region",
        "draw_gene_arrow",
    ):
        assert callable(getattr(gt_mod, fn_name)), f"{fn_name} not callable"


# ---------------------------------------------------------------------------
# assign_feature_levels
# ---------------------------------------------------------------------------


def test_assign_feature_levels_non_overlapping():
    features = [
        {"start": 100, "end": 200, "name": "a"},
        {"start": 1000, "end": 1100, "name": "b"},
        {"start": 2000, "end": 2100, "name": "c"},
    ]
    out = assign_feature_levels(features, min_gap=100)
    # All separated by > 100 bp → all at level 0
    assert all(f["level"] == 0 for f in out)


def test_assign_feature_levels_overlapping():
    features = [
        {"start": 100, "end": 1000, "name": "a"},
        {"start": 200, "end": 1100, "name": "b"},  # overlaps a
        {"start": 300, "end": 1200, "name": "c"},  # overlaps a + b
    ]
    out = assign_feature_levels(features, min_gap=100)
    levels = sorted(f["level"] for f in out)
    # 3 mutually overlapping features → 3 distinct levels
    assert levels == [0, 1, 2]


def test_assign_feature_levels_empty():
    assert assign_feature_levels([]) == []


# ---------------------------------------------------------------------------
# draw_gene_arrow
# ---------------------------------------------------------------------------


def test_draw_gene_arrow_plus_strand():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 2000)
    ax.set_ylim(-1, 1)
    draw_gene_arrow(ax, 100, 1000, "+", label="GENE1")
    # 1 polygon + 1 text label
    assert len(ax.patches) == 1
    assert any(t.get_text() == "GENE1" for t in ax.texts)


def test_draw_gene_arrow_minus_strand():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 2000)
    ax.set_ylim(-1, 1)
    draw_gene_arrow(ax, 100, 1000, "-")
    assert len(ax.patches) == 1


# ---------------------------------------------------------------------------
# draw_gene_track
# ---------------------------------------------------------------------------


def test_draw_gene_track_smoke(gff_features):
    fig, ax = plt.subplots()
    draw_gene_track(
        ax,
        gene_name="VPS8",
        gff_features=gff_features,
        region_start=0,
        region_end=6000,
    )
    # 3 CDS features → 3 polygons
    assert len(ax.patches) == 3
    # X-limits enforced
    assert ax.get_xlim() == (0, 6000)


def test_draw_gene_track_no_features_in_region(gff_features):
    """Region with no genes — must produce a 'No CDS' message, not raise."""
    fig, ax = plt.subplots()
    draw_gene_track(
        ax,
        gene_name="VPS8",
        gff_features=gff_features,
        region_start=20000,
        region_end=30000,
    )
    # Find the 'No CDS in region' message
    texts = [t.get_text() for t in ax.texts]
    assert any("No CDS" in t for t in texts), f"got {texts}"


def test_draw_gene_track_gene_absent(gff_features):
    """Unknown gene_name — fallback text rendered, no exception."""
    fig, ax = plt.subplots()
    draw_gene_track(
        ax,
        gene_name="NONEXISTENT_GENE",
        gff_features=gff_features,
        region_start=0,
        region_end=6000,
    )
    texts = [t.get_text() for t in ax.texts]
    assert any("no GFF" in t or "NONEXISTENT" in t for t in texts)


# ---------------------------------------------------------------------------
# get_genes_in_region
# ---------------------------------------------------------------------------


def test_get_genes_in_region_smoke(gff_features):
    genes = get_genes_in_region(gff_features, "chrI", 0, 6000)
    # Should contain all 3 systematic names
    assert "YAL001C" in genes
    assert "YAL002W" in genes
    assert "YAL003W" in genes
    # Maps to common names
    assert genes["YAL002W"] == "VPS8"


def test_get_genes_in_region_no_overlap(gff_features):
    genes = get_genes_in_region(gff_features, "chrI", 20000, 30000)
    assert genes == {}


def test_get_genes_in_region_unknown_chrom(gff_features):
    genes = get_genes_in_region(gff_features, "chrUnknown", 0, 6000)
    assert genes == {}
