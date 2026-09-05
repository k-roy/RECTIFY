"""Tests for rectify.visualize.tokens -- the bundled RNA figure design tokens.

Pins the contract, not the taste: the bundle is a complete two-layer token file, the loader
refuses a layer-B identity as a data role, and every default the drawing modules resolve
from it exists.
"""
import json
import os

import pytest

from rectify.visualize import tokens as TOK


def test_bundle_exists_and_is_two_layers():
    t = TOK.load()
    assert TOK.BUNDLED.exists()
    a, b = set(TOK.layer_a()), set(TOK.layer_b())
    sem = set(t["palette"]["semantic"])
    assert a and b and not (a & b), "layers must be disjoint"
    assert a | b == sem, "layers must partition the semantic palette"


def test_every_role_has_a_declared_role_block():
    t = TOK.load()
    declared = set()
    for spec in t["roles"].values():
        if isinstance(spec, dict):
            declared |= set(spec.get("tokens", []))
    assert set(TOK.palette()) <= declared


def test_role_refuses_layer_b_and_unknown():
    with pytest.raises(ValueError):
        TOK.role("polya")
    with pytest.raises(ValueError):
        TOK.role("chartreuse")
    assert TOK.role("focal").startswith("#")
    assert TOK.role("neutral").startswith("#")


def test_color_refuses_unknown_instead_of_defaulting():
    with pytest.raises(KeyError):
        TOK.color("not_a_token")


def test_defaults_used_by_drawing_modules_exist():
    for name in ("ink", "hairline", "subtle", "neutral", "mute", "wash", "paper",
                 "splice", "polya", "mod", "focal", "reference"):
        assert TOK.color(name).startswith("#")
    g = TOK.track_geometry()
    for key in ("exon_height", "utr_height_ratio", "chevron_spacing_in", "arrow_tip_fraction",
                "min_gap_bp", "coverage_fill_alpha", "arc_height_ratio", "mark_size_pt",
                "lollipop_stem_in", "tss_arrow_in"):
        assert key in g
    for key in ("in_figure", "axis_label", "tick_label", "annotation", "panel_letter"):
        assert key in TOK.typography()
    assert set(TOK.stroke()) >= {"primary", "secondary", "hairline"}


def test_role_cycle_is_layer_a_only():
    hexes = set(TOK.role_cycle())
    a = {TOK.color(n) for n in TOK.layer_a()}
    assert hexes <= a


def test_env_override_points_at_another_file(tmp_path, monkeypatch):
    alt = tmp_path / "alt_tokens.json"
    t = json.loads(TOK.BUNDLED.read_text())
    t["palette"]["semantic"]["focal"] = "#123456"
    alt.write_text(json.dumps(t))
    monkeypatch.setenv(TOK.ENV_VAR, str(alt))
    TOK.load(force=True)
    try:
        assert TOK.color("focal") == "#123456"
        assert TOK.source() == alt
    finally:
        monkeypatch.delenv(TOK.ENV_VAR)
        TOK.load(TOK.BUNDLED, force=True)


def test_apply_rc_sets_house_family_and_floored_sizes():
    import matplotlib
    rc = TOK.apply_rc(matplotlib)
    T = TOK.typography()
    assert rc["font.sans-serif"][0] == "Arial"
    assert rc["xtick.labelsize"] == T["tick_label"]
    assert rc["pdf.fonttype"] == 42
