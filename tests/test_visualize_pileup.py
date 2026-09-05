"""Tests for rectify.visualize.pileup -- rbrowse's read-alignment vocabulary on the
publication path: the merged raster, the ribbon wedge, the union-exon axis, isoform chains
and junction arcs. Pins the semantics that were ported, not the pixels:

  * a merged column is a MAJORITY (a minority tail never repaints the body);
  * a ribbon column is a FRACTION (the 5' staircase is an ECDF ramp);
  * the tail is the `polya` identity and the clip remainder is `mute`;
  * the spliced axis keeps covered/exonic segments, collapses the rest to a fixed gap,
    compresses long non-exonic segments, and to_t/from_t invert each other;
  * chains: exact junction key; retention through a boundary; DRS-suffix assignment;
    shared when two maximal parents; unspliced never shared; the dominant chain is starred;
  * isoform identity is a letter and never a hue.
"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402

from rectify.visualize import pileup as P  # noqa: E402
from rectify.visualize import tokens as TOK  # noqa: E402
from rectify.visualize import tracks as T  # noqa: E402


@pytest.fixture
def ax():
    fig, ax = plt.subplots(figsize=(4, 1.5))
    yield ax
    plt.close(fig)


def _hex(c):
    return to_hex(c).upper()


def R(i, s, e, blocks=None, tail=0, clip3=0, strand="+"):
    return P.Read(f"r{i}", strand, s, e, blocks or [(s, e - s)], tail=tail, clip3=clip3)


# ----------------------------------------------------------------- Read
def test_read_derives_junctions_from_block_gaps():
    r = P.Read("a", "+", 100, 500, [(100, 100), (300, 200)])
    assert r.junctions == [(200, 100)]
    assert r.three_prime == 500 and r.five_prime == 100
    assert r.spans_point(150) and not r.spans_point(250)


def test_read_from_rbrowse_record():
    rec = {"id": "x", "st": "-", "s": 10, "e": 60, "bl": [[10, 20], [40, 20]], "n": [[30, 10]], "tp": 12,
           "cl": "AAAAAAAAAAAAAAGT", "cr": "T"}
    r = P.Read.from_rbrowse(rec)
    assert r.strand == "-" and r.tail == 12 and r.clip3 == 16 and r.clip5 == 1
    assert r.junctions == [(30, 10)]


# ----------------------------------------------------------------- merged raster
def test_merged_column_is_a_majority(ax):
    # 8 reads with bodies 0-100; 1 of them has a 20-nt tail from 100 -- a minority
    reads = [R(i, 0, 100) for i in range(8)]
    reads[0] = R(0, 0, 100, tail=20, clip3=20)
    # and 8 reads whose bodies END at 60 with tails after -- tails win where bodies ended
    cb, ct, cc = P._columns(reads, 0, 120, 120)
    assert cb[50] == 8 and ct[50] == 0
    assert ct[110] == 1 and cb[110] == 0            # past the bodies the lone tail is present
    im = P.merged_reads(ax, reads, region=("chrI", 0, 120), nbins=120)
    img = im.get_array()
    body_rgb = np.array(matplotlib.colors.to_rgb(TOK.color("subtle")))
    tail_rgb = np.array(matplotlib.colors.to_rgb(TOK.color("polya")))
    assert np.allclose(img[0, 50, :3], body_rgb)    # bodies rule the column
    assert np.allclose(img[0, 110, :3], tail_rgb)   # tail wins once bodies have ended
    assert img[0, 50, 3] == pytest.approx(1.0)


def test_merged_minority_tail_does_not_repaint_bodies(ax):
    reads = [R(i, 0, 100) for i in range(8)] + [R(9, 0, 60, tail=40, clip3=40)]
    im = P.merged_reads(ax, reads, region=("chrI", 0, 120), nbins=120)
    img = im.get_array()
    body_rgb = np.array(matplotlib.colors.to_rgb(TOK.color("subtle")))
    assert np.allclose(img[0, 80, :3], body_rgb)    # 1 tail vs 8 bodies: bodies keep the column


def test_merged_takes_role_and_refuses_layer_b(ax):
    reads = [R(i, 0, 100) for i in range(3)]
    im = P.merged_reads(ax, reads, role="focal", region=("chrI", 0, 120), nbins=60)
    assert np.allclose(im.get_array()[0, 30, :3], matplotlib.colors.to_rgb(TOK.color("focal")))
    with pytest.raises(ValueError):
        P.merged_reads(ax, reads, role="polya", region=("chrI", 0, 120), nbins=60)


# ----------------------------------------------------------------- ribbon
def test_ribbon_heights_are_fractions(ax):
    reads = [R(i, s, 100, tail=10, clip3=10) for i, s in enumerate([0, 20, 40, 60])]
    out = P.ribbon(ax, reads, y=0, h=1.0, anchor=100, region=("chrI", 0, 120), nbins=120, letter="A")
    body = out[0]
    verts = body.get_paths()[0].vertices
    ys = verts[:, 1]
    assert ys.max() == pytest.approx(1.0, abs=0.02)         # all four present near the anchor
    # at x=10 only one read (the one starting at 0) is present -> quarter height
    xs = verts[:, 0]
    near = ys[(xs > 8) & (xs < 12)]
    assert near.max() == pytest.approx(0.25, abs=0.05)
    tail = out[1]
    assert _hex(tail.get_facecolor()[0]) == TOK.color("polya").upper()
    texts = [a for a in out if hasattr(a, "get_text")]
    assert texts and texts[0].get_text() == "A"


# ----------------------------------------------------------------- spliced axis
def test_spliced_axis_keeps_covered_and_exonic_collapses_the_rest():
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 200), (500, 600), (900, 1000)])
    reads = [P.Read("a", "+", 100, 1000, [(100, 100), (500, 100), (900, 100)])]
    ax_ = P.SplicedAxis((0, 1100), reads, [tx])
    kept = [(sg.s, sg.e) for sg in ax_.segs if sg.keep]
    assert (100, 200) in kept and (500, 600) in kept and (900, 1000) in kept
    assert (200, 500) not in kept and (600, 900) not in kept       # introns collapse
    # each collapsed intron takes exactly one gap
    gaps = [sg for sg in ax_.segs if not sg.keep and sg.t1 > sg.t0]
    assert all(sg.t1 - sg.t0 == ax_.gap_bp for sg in gaps)
    # monotone, invertible
    for p in (100, 150, 500, 950):
        assert ax_.from_t(ax_.to_t(p)) == pytest.approx(p)
    assert ax_.to_t(500) - ax_.to_t(200) == pytest.approx(ax_.gap_bp)   # the intron is one gap wide


def test_spliced_axis_compresses_long_non_exonic_segments():
    tx = T.Transcript("g", "chrI", "+", exons=[(0, 100), (5000, 5100)])
    reads = [P.Read("a", "+", 0, 5100, [(0, 5100)])]          # retention through the whole intron
    ax_ = P.SplicedAxis((0, 5200), reads, [tx])
    intron = [sg for sg in ax_.segs if (sg.s, sg.e) == (100, 5000)][0]
    assert intron.keep and intron.compressed and intron.wbp == ax_.seg_cap_bp


def test_spliced_axis_strict_has_no_gaps():
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 200), (500, 600)])
    reads = [P.Read("a", "+", 100, 600, [(100, 100), (500, 100)])]
    ax_ = P.SplicedAxis((0, 700), reads, [tx], strict=True)
    assert ax_.t_end == pytest.approx(200)


# ----------------------------------------------------------------- chains
def _model():
    return T.Transcript("g", "chrI", "+", exons=[(100, 200), (300, 400), (500, 600), (700, 800)])


def test_chains_fl_skip_retention_and_dominant():
    m = _model()
    fl = [P.Read(f"fl{i}", "+", 100, 800, [(100, 100), (300, 100), (500, 100), (700, 100)]) for i in range(5)]
    skip = [P.Read(f"sk{i}", "+", 100, 800, [(100, 100), (300, 100), (700, 100)]) for i in range(3)]
    ret = [P.Read("rt", "+", 100, 800, [(100, 100), (300, 300), (700, 100)])]     # runs through intron 2
    # with NO annotated skip transcript the skipping junction is novel: 'Δe3 novel' (rbrowse's rule)
    chains, _ = P.chain_clusters(fl + skip + ret, (0, 900), m)
    assert any(c.label == "Δe3 novel" and c.cls == "skip" and c.glyph == "✱" for c in chains)
    # with the skip ANNOTATED it is a plain skip, named after its transcript
    skip_tx = T.Transcript("g-202", "chrI", "+", exons=[(100, 200), (300, 400), (700, 800)])
    chains, summary = P.chain_clusters(fl + skip + ret, (0, 900), m, annotated=[m, skip_tx])
    by = {c.label: c for c in chains}
    assert chains[0].dominant and chains[0].label == "FL" and chains[0].n == 5
    assert "Δe3" in by and by["Δe3"].cls == "skip" and by["Δe3"].glyph == "⊘" and by["Δe3"].tx_name == "g-202"
    assert any(c.cls == "ret" and c.label.startswith("ret-i2") for c in chains)
    assert summary["unambiguous"] == 8 and summary["retention"] == 1


def test_chains_suffix_assignment_and_shared():
    m = _model()
    fl = [P.Read(f"fl{i}", "+", 100, 800, [(100, 100), (300, 100), (500, 100), (700, 100)]) for i in range(6)]
    skip = [P.Read(f"sk{i}", "+", 100, 800, [(100, 100), (300, 100), (700, 100)]) for i in range(6)]
    # a chain with only the FIRST junction (e1-e2), which FL and De3 both carry, and no
    # bases where they differ: compatible with both maximal chains -> shared
    suf = [P.Read(f"su{i}", "+", 150, 400, [(150, 50), (300, 100)]) for i in range(2)]
    # a suffix with the last TWO junctions (e2-e3, e3-e4): only FL extends it -> assigned
    suf2 = [P.Read(f"sv{i}", "+", 320, 800, [(320, 80), (500, 100), (700, 100)]) for i in range(2)]
    chains, summary = P.chain_clusters(fl + skip + suf + suf2, (0, 900), m)
    shared = [c for c in chains if c.shared]
    assert shared and shared[0].glyph == "≈" and shared[0].label.startswith("shared")
    assigned = [c for c in chains if c.assigned_to]
    assert assigned and assigned[0].label == "FL-compatible from e2"   # DRS suffix: says where it starts
    assert summary["shared"] == 2 and summary["compatible"] == 2


def test_unspliced_is_never_shared():
    m = _model()
    reads = [P.Read("u", "+", 320, 380, [(320, 60)])]
    chains, summary = P.chain_clusters(reads, (0, 900), m)
    assert chains[0].cls == "unspliced" and not chains[0].shared and summary["unspliced"] == 1


def test_isoform_rows_use_letters_not_hues(ax):
    m = _model()
    fl = [P.Read(f"fl{i}", "+", 100, 800, [(100, 100), (300, 100), (500, 100), (700, 100)]) for i in range(4)]
    skip = [P.Read(f"sk{i}", "+", 100, 800, [(100, 100), (300, 100), (700, 100)]) for i in range(2)]
    chains, _ = P.chain_clusters(fl + skip, (0, 900), m)
    axis = P.SplicedAxis((0, 900), fl + skip, [m])
    ys = P.isoform_rows(ax, chains, axis, denominator=6)
    assert len(ys) == 2
    letters = [t.get_text() for t in ax.texts if t.get_text() in ("A", "B")]
    assert letters == ["A", "B"]
    fills = {_hex(p.get_facecolor()) for p in ax.patches}
    assert fills == {TOK.color("subtle").upper()}          # one colour, identity by letter


# ----------------------------------------------------------------- junction arcs
def test_junction_counts_and_arcs(ax):
    reads = [P.Read(f"a{i}", "+", 100, 400, [(100, 100), (300, 100)]) for i in range(3)] + \
            [P.Read("b", "+", 100, 400, [(100, 50), (350, 50)])]
    counts = P.junction_counts(reads, (0, 500))
    assert counts == {(200, 300): 3, (150, 350): 1}
    ax.set_ylim(0, 10)
    out = P.junction_arcs(ax, counts, y=0, height=4, role="focal", annotated=[(200, 300), (250, 320)])
    patches = [a for a in out if a.__class__.__name__ == "PathPatch"]
    assert len(patches) == 3      # two supported junctions + one unsupported annotated arc
    cols = {_hex(p.get_edgecolor()) for p in patches}
    assert TOK.color("focal").upper() in cols and TOK.color("splice").upper() in cols


# ----------------------------------------------------------------- audit 879: the rules the first port lacked
def test_parent_floor_matches_rbrowse():
    """spliced.js:24-25 -- the SDHA incident set these; the browser and the figure must agree."""
    assert (P.MIN_PARENT_N, P.MIN_PARENT_FRAC, P.NEAR_DUP_NT) == (2, 0.2, 3)


def test_parent_floor_regime_where_it_separates():
    # chain A (n=10, one junction) whose extension B has n=3: with 2/0.2 B is A's parent (3 >= max(2, 2))
    m = _model()
    a = [P.Read(f"a{i}", "+", 100, 600, [(100, 100), (300, 100), (500, 100)]) for i in range(10)]
    b = [P.Read(f"b{i}", "+", 100, 800, [(100, 100), (300, 100), (500, 100), (700, 100)]) for i in range(3)]
    chains, summ = P.chain_clusters(a + b, (0, 900), m)
    assert summ["unambiguous"] == 3 and summ["compatible"] == 10
    assert any(c.label.startswith("FL-compatible") and c.n == 10 for c in chains)


def test_fl_from_exon_suffix():
    m = _model()
    fl = [P.Read(f"fl{i}", "+", 100, 800, [(100, 100), (300, 100), (500, 100), (700, 100)]) for i in range(4)]
    from3 = [P.Read(f"s{i}", "+", 520, 800, [(520, 80), (700, 100)]) for i in range(2)]
    chains, _ = P.chain_clusters(fl + from3, (0, 900), m)
    labels = {c.label for c in chains}
    assert "FL" in labels and "FL-compatible from e3" in labels


def test_read_refuses_start_end_blocks():
    with pytest.raises(ValueError):
        P.Read("x", "+", 100, 400, [(100, 200), (300, 400)])     # (start, end) given as blocks
    P.Read("ok", "+", 100, 400, [(100, 100), (300, 100)])        # (start, length) is fine


def test_spliced_axis_keeps_covered_3prime_flank():
    """spliced.js DOWNSTREAM_KEEP_BP: a covered 3' flank past the last exon stays 1:1 (not compressed)."""
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 300), (500, 700)])
    # 20 reads run 1,200 bp past the annotated end; 1 read stops at the annotation
    reads = [P.Read(f"r{i}", "+", 100, 1900, [(100, 200), (500, 1400)]) for i in range(20)]
    reads.append(P.Read("short", "+", 100, 700, [(100, 200), (500, 200)]))
    ax_ = P.SplicedAxis((0, 2200), reads, [tx])
    # read ends are not cut points, so the covered flank is ONE segment from 700 to the scope edge
    flank = [sg for sg in ax_.segs if sg.s >= 700 and sg.cov >= 20]
    assert flank and all(sg.keep and not sg.compressed for sg in flank)
    assert ax_.to_t(1900) - ax_.to_t(700) == pytest.approx(1200)   # 1:1, the whole covered flank
    # and the same flank with no coverage compresses as before
    ax2 = P.SplicedAxis((0, 2200), [reads[-1]], [tx])
    assert all((not sg.keep) or sg.compressed for sg in ax2.segs if sg.s >= 700 and sg.e <= 1900 and not sg.exonic)


def test_exon_only_reports_hidden_reads():
    tx = T.Transcript("g", "chrI", "+", exons=[(100, 200), (500, 600)])
    spliced = [P.Read(f"s{i}", "+", 100, 600, [(100, 100), (500, 100)]) for i in range(6)]
    retained = [P.Read(f"r{i}", "+", 100, 600, [(100, 500)]) for i in range(2)]
    ax_ = P.SplicedAxis((0, 700), spliced + retained, [tx], exon_only=True)
    assert len(ax_.hidden) == 1
    sg = ax_.hidden[0]
    assert (sg.s, sg.e) == (200, 500) and sg.hidden_n == 2 and sg.across == 6
    assert sg.hidden_frac == pytest.approx(0.25)


def test_junction_arcs_validates_role_even_when_empty(ax):
    with pytest.raises(ValueError):
        P.junction_arcs(ax, {}, role="polya")
