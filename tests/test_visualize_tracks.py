"""Tests for rectify.visualize.tracks -- the token-driven layer-B glyphs.

What these pin down (each was a design rule, not a taste):
  * a gene model and a mark take NO colour argument, and the colours they use are the
    token hexes (ink / hairline / splice / polya / mod), never anything else;
  * a signal glyph takes a ROLE and refuses a layer-B token (the boundary rule);
  * strand is position: a minus-strand Region reverses the x-axis, strand_coverage draws
    the minus strand below zero in the SAME colour as the plus strand;
  * the three marks are three different SHAPES (the declared greyscale cue);
  * UTRs draw at the token ratio of the exon height; introns are derived correctly.
"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402

from rectify.visualize import tokens as TOK  # noqa: E402
from rectify.visualize import tracks as T  # noqa: E402


@pytest.fixture
def ax():
    fig, ax = plt.subplots(figsize=(4, 1.5))
    yield ax
    plt.close(fig)


def _hex(c):
    return to_hex(c).upper()


# ----------------------------------------------------------------- data carriers
def test_transcript_from_introns_derives_exons_and_introns():
    tx = T.Transcript.from_introns("PDH1", "chrXVI", "+", 558219, 560137, introns=[(558423, 559811)])
    assert tx.exons == [(558219, 558423), (559811, 560137)]
    assert tx.introns == [(558423, 559811)]
    assert tx.tss == 558219 and tx.tes == 560137


def test_transcript_minus_strand_tss_is_high_coordinate():
    tx = T.Transcript("SWC4", "chrVII", "-", exons=[(498168, 498448), (498700, 499924)])
    assert tx.tss == 499924 and tx.tes == 498168


def test_transcript_rejects_bad_input():
    with pytest.raises(ValueError):
        T.Transcript("x", "chrI", "*", exons=[(1, 10)])
    with pytest.raises(ValueError):
        T.Transcript("x", "chrI", "+", exons=[(10, 10)])
    with pytest.raises(ValueError):
        T.Transcript.from_introns("x", "chrI", "+", 100, 200, introns=[(50, 150)])


def test_utr_cds_parts_split_each_exon():
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 200), (300, 400)], cds=(150, 350))
    parts = tx.utr_cds_parts()
    assert parts == [((100, 150), False), ((150, 200), True), ((300, 350), True), ((350, 400), False)]


def test_region_reversed_only_for_minus_strand():
    assert not T.Region("chrI", 1, 100).reversed
    assert T.Region("chrI", 1, 100, "-").reversed
    assert T.Region("chrI", 1, 100, "-").xlim() == (100, 1)
    with pytest.raises(ValueError):
        T.Region("chrI", 100, 1)


# ----------------------------------------------------------------- gene model
def test_gene_model_uses_only_grey_ramp(ax):
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 200), (300, 400)], cds=(150, 350))
    out = T.gene_model(ax, tx, region=("chrI", 0, 500), tss=True, intron_style="chevron")
    ink, hairline = TOK.color("ink").upper(), TOK.color("hairline").upper()
    assert {_hex(p.get_facecolor()) for p in out["exons"] + out["utrs"]} == {ink}
    assert {_hex(l.get_color()) for l in out["introns"] + out["chevrons"]} == {hairline}
    assert len(out["exons"]) == 2 and len(out["utrs"]) == 2 and len(out["introns"]) == 1


def test_gene_model_utr_height_is_token_ratio(ax):
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 200)], cds=(150, 200))
    out = T.gene_model(ax, tx, height=1.0, region=("chrI", 0, 300), label=False)
    ratio = TOK.track_geometry()["utr_height_ratio"]
    assert out["exons"][0].get_height() == pytest.approx(1.0)
    assert out["utrs"][0].get_height() == pytest.approx(ratio)


def test_gene_model_takes_no_colour_argument(ax):
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 200)])
    with pytest.raises(TypeError):
        T.gene_model(ax, tx, color="#FF0000")


def test_gene_model_label_is_italic(ax):
    tx = T.Transcript("PDH1", "chrXVI", "+", exons=[(100, 200)])
    out = T.gene_model(ax, tx, region=("chrXVI", 0, 300))
    assert out["label"][0].get_fontstyle() == "italic"
    assert out["label"][0].get_text() == "PDH1"


def test_gene_track_priority_takes_top_row(ax):
    a = T.Transcript("A", "chrI", "+", exons=[(100, 900)])
    b = T.Transcript("B", "chrI", "-", exons=[(500, 1400)])
    ys = T.gene_track(ax, [a, b], ("chrI", 0, 1500), priority=["B"], label_pos="left")
    assert ys["B"] == 0.0 and ys["A"] == -1.0
    ys2 = T.gene_track(ax, [a, b], ("chrI", 0, 1500), label_pos="left")
    assert ys2["A"] == 0.0 and ys2["B"] == -1.0


def test_gene_track_rows_are_spaced_wider_when_labels_sit_above(ax):
    a = T.Transcript("A", "chrI", "+", exons=[(100, 900)])
    b = T.Transcript("B", "chrI", "-", exons=[(500, 1400)])
    ys = T.gene_track(ax, [a, b], ("chrI", 0, 1500), label_pos="above")
    assert ys["B"] < -1.0                      # a label above row 1 must not rise into row 0
    ys2 = T.gene_track(ax, [a, b], ("chrI", 0, 1500), label_pos="above", row_pitch=2.5)
    assert ys2["B"] == -2.5


def test_gene_track_skips_transcripts_outside_region(ax):
    a = T.Transcript("A", "chrI", "+", exons=[(100, 200)])
    far = T.Transcript("Z", "chrI", "+", exons=[(5000, 6000)])
    ys = T.gene_track(ax, [a, far], ("chrI", 0, 1000))
    assert set(ys) == {"A"}


# ----------------------------------------------------------------- marks
def test_marks_are_three_shapes_in_their_own_colours(ax):
    ax.set_xlim(0, 100)
    tick = T.mark(ax, 10, "ss5")
    lol = T.mark(ax, 40, "polya")
    mod = T.mark(ax, 70, "mod")
    # colours are the layer-B tokens
    assert _hex(tick[0].get_color()) == TOK.color("splice").upper()
    assert _hex(lol[1].get_markerfacecolor()) == TOK.color("polya").upper()
    assert _hex(mod[1].get_markerfacecolor()) == TOK.color("mod").upper()
    # shapes differ: tick has no head marker beyond the triangle, lollipops differ in head
    assert lol[1].get_marker() == "o" and mod[1].get_marker() == "D"
    assert tick[1].get_marker() in (">", "<")


def test_pas_is_open_and_cpa_is_filled(ax):
    ax.set_xlim(0, 100)
    pas = T.mark(ax, 10, "pas")
    cpa = T.mark(ax, 20, "cpa")
    assert _hex(pas[1].get_markerfacecolor()) == TOK.color("paper").upper()
    assert _hex(cpa[1].get_markerfacecolor()) == TOK.color("polya").upper()


def test_splice_triangle_points_into_the_intron(ax):
    ax.set_xlim(0, 100)
    assert T.mark(ax, 10, "ss5", strand="+")[1].get_marker() == ">"
    assert T.mark(ax, 90, "ss3", strand="+")[1].get_marker() == "<"
    ax.set_xlim(100, 0)   # reversed axis, minus strand: 5'SS is on the LEFT, intron to its right
    assert T.mark(ax, 90, "ss5", strand="-")[1].get_marker() == ">"
    assert T.mark(ax, 10, "ss3", strand="-")[1].get_marker() == "<"


def test_mark_rejects_unknown_kind_and_colour(ax):
    with pytest.raises(ValueError):
        T.mark(ax, 1, "banana")
    with pytest.raises(TypeError):
        T.mark(ax, 1, "polya", color="#FF0000")


# ----------------------------------------------------------------- signal glyphs take a ROLE
def test_coverage_default_is_subtle_and_role_is_layer_a(ax):
    pos = np.arange(0, 100)
    dep = np.ones(100)
    poly = T.coverage(ax, pos, dep, region=("chrI", 0, 100))[0]
    assert _hex(poly.get_facecolor()[0]) == TOK.color("subtle").upper()
    poly2 = T.coverage(ax, pos, dep, role="focal")[0]
    assert _hex(poly2.get_facecolor()[0]) == TOK.color("focal").upper()


def test_coverage_refuses_layer_b_token(ax):
    with pytest.raises(ValueError):
        T.coverage(ax, np.arange(10), np.ones(10), role="polya")
    with pytest.raises(ValueError):
        T.arc(ax, 1, 5, role="splice")


def test_strand_coverage_mirrors_in_one_colour(ax):
    pos = np.arange(0, 100)
    out = T.strand_coverage(ax, pos, np.ones(100), pos, np.ones(100) * 2, region=("chrI", 0, 100))
    plus, minus = out[0], out[1]
    assert _hex(plus.get_facecolor()[0]) == _hex(minus.get_facecolor()[0])
    assert minus.get_paths()[0].vertices[:, 1].min() < 0   # drawn below the axis


def test_arc_annotation_is_splice_and_signal_is_role(ax):
    ax.set_ylim(0, 10)
    a = T.arc(ax, 10, 50, y=0, height=4)[0]
    assert _hex(a.get_edgecolor()) == TOK.color("splice").upper()
    b = T.arc(ax, 10, 50, y=0, height=4, role="reference", label="12")[0]
    assert _hex(b.get_edgecolor()) == TOK.color("reference").upper()


def test_reads_stack_and_default_grey(ax):
    st = np.array([10, 20, 30, 15])
    en = np.array([50, 60, 70, 55])
    n = T.reads(ax, st, en, region=("chrI", 0, 100))
    assert n >= 2
    cols = {_hex(c) for coll in ax.collections for c in coll.get_colors()}
    assert TOK.color("subtle").upper() in cols


# ----------------------------------------------------------------- the axis
def test_region_axis_reverses_for_minus_strand(ax):
    T.region_axis(ax, T.Region("chrVII", 497900, 500100, "-"))
    lo, hi = ax.get_xlim()
    assert lo > hi
    assert "chrVII" in ax.get_xlabel() and "kb" in ax.get_xlabel()


def test_region_axis_explicit_unit(ax):
    T.region_axis(ax, ("chrI", 0, 5000), unit="bp")
    assert "(bp)" in ax.get_xlabel()


# ----------------------------------------------------------------- audit 879 additions
def test_mark_label_alignment(ax):
    ax.set_xlim(0, 100)
    lab = T.mark(ax, 50, "polya", label="x", ha="right")[-1]
    assert lab.get_ha() == "right"
    with pytest.raises(ValueError):
        T.mark(ax, 50, "polya", label="x", ha="middle")


def test_reads_and_strand_coverage_take_xform(ax):
    xf = lambda p: p / 10.0
    n = T.reads(ax, [100, 200], [400, 500], junctions=[[(200, 300)], []], xform=xf)
    assert n >= 1
    lo, hi = ax.get_xlim()
    assert hi <= 60      # drawn in transformed space, not genome space
    pos = np.arange(0, 100)
    out = T.strand_coverage(ax, pos, np.ones(100), pos, np.ones(100), xform=xf)
    assert out[0].get_paths()[0].vertices[:, 0].max() <= 10.5
