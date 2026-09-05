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
# ============================================================================
# SPACE ECONOMY -- the mm-budget arithmetic (Chanfreau planning/880)
#
# These pin the three claims the layout makes on the page, so a refactor that
# silently overflows a panel fails here instead of on a proof:
#   1. a panel never exceeds its mm budget;
#   2. a kept band is never below the mm floor;
#   3. a stack's height is the sum of its budgets plus the declared chrome.
# `test_budget_checker_can_fail` plants a violation of each and asserts the
# checker rejects it -- a checker that cannot fail is not a checker.
# ============================================================================
MM = 25.4


def _reads_at(n, anchor, *, start=1000, tail=12, tag="x"):
    """n reads sharing a 3' end at ``anchor``, with staggered 5' ends."""
    out = []
    for i in range(n):
        s = start + (i % 7) * 40
        out.append(P.Read(f"{tag}{anchor}_{i}", "+", s, anchor, [(s, anchor - s)],
                          tail=tail, clip3=tail + 2))
    return out


def check_panel_budget(plan, *, floor_mm):
    """The invariants, as assertions. Raises AssertionError on a violation."""
    assert plan.used_mm <= plan.budget_mm + 1e-6, (
        f"{plan.name}: used {plan.used_mm:.4f} mm > budget {plan.budget_mm:.4f} mm")
    for b in plan.bands:
        assert b.h_mm >= floor_mm - 1e-9, (
            f"{plan.name}: band {b.letter or 'other'} is {b.h_mm:.4f} mm, below the "
            f"{floor_mm} mm floor")
    if plan.mode == "squished":
        pf = TOK.track_geometry()["pitch_floor_mm"]
        assert plan.pitch_mm >= pf - 1e-9, f"{plan.name}: pitch {plan.pitch_mm} < floor {pf}"
        assert plan.n_rows * plan.pitch_mm <= plan.bands_mm + 1e-6, (
            f"{plan.name}: {plan.n_rows} rows x {plan.pitch_mm:.4f} mm leaves the "
            f"{plan.bands_mm:.4f} mm strip")
        drawn_mm = plan.n_rows * plan.lw_pt / 72.0 * MM
        assert drawn_mm <= plan.bands_mm + 1e-6, (
            f"{plan.name}: {drawn_mm:.2f} mm of drawn ink in a {plan.bands_mm:.2f} mm strip")
    return True


@pytest.fixture
def stack_samples():
    """Four synthetic samples over five 3'-end clusters, seeded."""
    rng = np.random.default_rng(880)
    anchors = [4000, 3700, 3400, 3100, 2800]
    out = []
    for si, (name, w) in enumerate([("WT", [60, 30, 12, 6, 3]), ("mutA", [20, 70, 15, 5, 2]),
                                    ("mutB", [40, 40, 20, 10, 5]), ("mutC", [90, 8, 4, 2, 1])]):
        reads = []
        for a, n in zip(anchors, w):
            reads.extend(_reads_at(n, a + int(rng.integers(-4, 5)), tag=f"s{si}"))
        out.append({"name": name, "reads": reads, "role": None})
    return out


def test_cluster_by_three_prime_groups_within_window():
    reads = _reads_at(5, 4000) + _reads_at(3, 4010) + _reads_at(4, 3000)
    cl = P.cluster_by_three_prime(reads, win=32)
    assert [c.n for c in cl] == [8, 4]                 # 4000/4010 fuse; 3000 apart
    assert cl[0].anchor in (4000, 4010) and cl[1].anchor == 3000


def test_union_clusters_rank_and_letter_across_samples(stack_samples):
    keep = P.union_clusters(stack_samples, k=4)
    assert [c.letter for c in keep] == ["A", "B", "C", "D"]
    counts = [c.n for c in keep]
    assert counts == sorted(counts, reverse=True)      # A is the biggest POOLED cluster


def test_panel_never_exceeds_its_budget(stack_samples):
    floor = TOK.track_geometry()["band_floor_mm"]
    keep = P.union_clusters(stack_samples, k=4)
    for budget in (8.0, 10.0, 12.0, 18.0, 25.0):
        bands_mm = budget - (2.5 + 0.8)
        counts = [P.effective_counts(s["reads"], keep, bands_mm=bands_mm, floor_mm=floor,
                                     band_gap_mm=0.4) for s in stack_samples]
        scale = P.solve_band_scale(counts, bands_mm, floor_mm=floor, band_gap_mm=0.4)
        for s in stack_samples:
            plan = P.plan_read_panel(s["name"], s["reads"], keep, budget_mm=budget, scale=scale)
            check_panel_budget(plan, floor_mm=floor)


def test_band_floor_holds_and_k_is_reduced_when_it_cannot(stack_samples):
    """A budget too small for k+1 floors POOLS clusters; it never shrinks a band."""
    floor = TOK.track_geometry()["band_floor_mm"]
    keep = P.union_clusters(stack_samples, k=4)
    scale = P.BandScale(0.0, floor)                   # every band at the floor
    s = stack_samples[0]
    wide = P.plan_read_panel("wide", s["reads"], keep, budget_mm=20.0, scale=scale)
    tight = P.plan_read_panel("tight", s["reads"], keep, budget_mm=6.0, scale=scale)
    assert len(tight.bands) < len(wide.bands)
    assert tight.n_pooled > wide.n_pooled             # what was dropped is VISIBLE as `other`
    for plan in (wide, tight):
        check_panel_budget(plan, floor_mm=floor)
        assert sum(b.n for b in plan.bands) == plan.n_reads   # no read is lost


def test_shared_scale_makes_heights_comparable(stack_samples):
    """Equal counts in two different samples get equal band heights -- the point of one
    scale across the stack (rbrowse `bandScale: 'shared'`)."""
    floor = TOK.track_geometry()["band_floor_mm"]
    keep = P.union_clusters(stack_samples, k=4)
    counts = [P.effective_counts(s["reads"], keep, bands_mm=8.7, floor_mm=floor,
                                 band_gap_mm=0.4) for s in stack_samples]
    scale = P.solve_band_scale(counts, 8.7, floor_mm=floor, band_gap_mm=0.4)
    plans = [P.plan_read_panel(s["name"], s["reads"], keep, budget_mm=12.0, scale=scale)
             for s in stack_samples]
    by_n = {}
    for pl in plans:
        for b in pl.bands:
            by_n.setdefault(b.n, set()).add(round(b.h_mm, 9))
    for n, hs in by_n.items():
        assert len(hs) == 1, f"count {n} drew {hs} -- the scale is not shared"
    # and the scale is the LARGEST that fits: nudging it up overflows some panel
    bigger = P.BandScale(scale.mm_per_read * 1.25, floor)
    assert not all(P._fits(c, bigger, 8.7, 0.4) for c in counts)


def test_squished_only_under_the_pitch_floor():
    floor = TOK.track_geometry()["band_floor_mm"]
    pf = TOK.track_geometry()["pitch_floor_mm"]
    scale = P.BandScale(0.02, floor)
    budget = 25.4 + 2.5 + 0.8                          # a 1.0 in band strip
    keep = []
    few = P.plan_read_panel("few", _reads_at(40, 4000), keep, budget_mm=budget, scale=scale)
    many = P.plan_read_panel("many", _reads_at(400, 4000), keep, budget_mm=budget, scale=scale)
    assert few.mode == "squished" and many.mode == "merged"
    assert few.pitch_mm >= pf
    check_panel_budget(few, floor_mm=floor)
    check_panel_budget(many, floor_mm=floor)
    # the drawn line width is DERIVED from the pitch -- the fixed 4.0 pt of tracks.reads
    # would need 40 * 1.411 = 56.4 mm of ink in a 25.4 mm strip
    assert few.lw_pt < 4.0
    assert few.n_rows * few.lw_pt / 72.0 * MM <= few.bands_mm + 1e-6
    # exactly at the threshold: 25.4 / 0.423 = 60 reads
    edge = P.plan_read_panel("edge", _reads_at(60, 4000), keep, budget_mm=budget, scale=scale)
    assert edge.mode == "squished"
    over = P.plan_read_panel("over", _reads_at(61, 4000), keep, budget_mm=budget, scale=scale)
    assert over.mode == "merged"


def test_pinned_scale_shaves_and_declares_the_cut(stack_samples):
    """A caller may PIN a scale (one shared with another figure). When it does not fit, the
    plan shaves proportionally, keeps the floor, stays inside the budget, and flags the
    bands `capped` so `ribbon` draws the dashed cut mark (rbrowse.js:17577)."""
    floor = TOK.track_geometry()["band_floor_mm"]
    keep = P.union_clusters(stack_samples, k=4)
    plan = P.plan_read_panel("pinned", stack_samples[0]["reads"], keep, budget_mm=12.0,
                             scale=P.BandScale(0.5, floor))     # far too generous
    assert plan.shaved_mm > 0
    assert any(b.capped for b in plan.bands)
    check_panel_budget(plan, floor_mm=floor)


def test_stack_height_is_the_sum_of_the_budgets(stack_samples):
    h = P.stack_height_mm(4, panel_mm=12.0, panel_gap_mm=2.5, model_mm=7.0, model_gap_mm=2.0,
                          axis_mm=8.0, top_mm=3.0, reserved_mm=30.0)
    assert h == pytest.approx(3.0 + 7.0 + 2.0 + 4 * 12.0 + 3 * 2.5 + 8.0 + 30.0)
    # and it is LINEAR in the sample count: one more sample costs exactly one pitch
    h5 = P.stack_height_mm(5, panel_mm=12.0, panel_gap_mm=2.5, model_mm=7.0, model_gap_mm=2.0,
                           axis_mm=8.0, top_mm=3.0, reserved_mm=30.0)
    assert h5 - h == pytest.approx(12.0 + 2.5)


def test_read_stack_renders_at_the_budgeted_geometry(stack_samples):
    """The rendered artists, not the intent: every panel axes is exactly `panel_mm` tall and
    the figure is exactly `stack_height_mm` tall."""
    tx = T.Transcript("GENE1", "chrI", "+", exons=[(1000, 1300), (1800, 2000), (3400, 4100)])
    reg = T.Region("chrI", 900, 4200)
    fig = plt.figure(figsize=(7.2, 3.0), layout="none")
    try:
        res = P.read_stack(fig, stack_samples, tx, region=reg, panel_mm=12.0, panel_gap_mm=2.5,
                           model_mm=7.0, axis_mm=8.0, top_mm=3.0, reserved_mm=0.0)
        H_mm = fig.get_figheight() * MM
        assert H_mm == pytest.approx(res["height_mm"], abs=1e-6)
        for ax in res["panels"]:
            h = ax.get_position().height * H_mm
            assert h == pytest.approx(12.0, abs=1e-6)
        assert res["model_ax"].get_position().height * H_mm == pytest.approx(7.0, abs=1e-6)
        floor = TOK.track_geometry()["band_floor_mm"]
        for plan in res["plans"]:
            check_panel_budget(plan, floor_mm=floor)
        # one model, one axis, N panels -- and nothing else
        assert len(fig.axes) == len(stack_samples) + 2
    finally:
        plt.close(fig)


def test_budget_checker_can_fail(stack_samples):
    """PLANT each violation and prove the checker rejects it."""
    floor = TOK.track_geometry()["band_floor_mm"]
    keep = P.union_clusters(stack_samples, k=4)
    counts = [P.effective_counts(s["reads"], keep, bands_mm=8.7, floor_mm=floor,
                                 band_gap_mm=0.4) for s in stack_samples]
    scale = P.solve_band_scale(counts, 8.7, floor_mm=floor, band_gap_mm=0.4)
    base = P.plan_read_panel("plant", stack_samples[0]["reads"], keep, budget_mm=12.0, scale=scale)
    check_panel_budget(base, floor_mm=floor)                      # clean first

    over = P.plan_read_panel("over", stack_samples[0]["reads"], keep, budget_mm=12.0, scale=scale)
    over.bands[0].h_mm += 5.0                                     # 1. blow the budget
    with pytest.raises(AssertionError, match="> budget"):
        check_panel_budget(over, floor_mm=floor)

    under = P.plan_read_panel("under", stack_samples[0]["reads"], keep, budget_mm=12.0, scale=scale)
    under.bands[-1].h_mm = floor / 2                               # 2. break the floor
    with pytest.raises(AssertionError, match="below the"):
        check_panel_budget(under, floor_mm=floor)

    sq = P.plan_read_panel("sq", _reads_at(40, 4000), (), budget_mm=28.7,
                           scale=P.BandScale(0.02, floor))
    assert sq.mode == "squished"
    sq.lw_pt = 4.0                                                 # 3. the tracks.reads defect
    with pytest.raises(AssertionError, match="drawn ink"):
        check_panel_budget(sq, floor_mm=floor)
    sq.lw_pt = 0.6
    sq.pitch_mm = 0.1                                              # 4. below the pitch floor
    with pytest.raises(AssertionError, match="pitch"):
        check_panel_budget(sq, floor_mm=floor)


def check_stack_height(fig, res, *, n, panel_mm, panel_gap_mm, **chrome):
    """Invariant 3: the rendered page height IS the sum of the panel budgets plus the
    declared chrome -- and every panel axes is exactly its budget."""
    H_mm = fig.get_figheight() * MM
    want = P.stack_height_mm(n, panel_mm=panel_mm, panel_gap_mm=panel_gap_mm, **chrome)
    assert H_mm == pytest.approx(want, abs=1e-6), (
        f"page is {H_mm:.4f} mm, the budgets sum to {want:.4f} mm")
    # PER PANEL, not the sum: a sum is blind to two panels swapping heights, or to one
    # shrinking by exactly what another grows -- both of which are broken layouts.
    for i, ax in enumerate(res["panels"]):
        h = ax.get_position().height * H_mm
        assert h == pytest.approx(panel_mm, abs=1e-6), (
            f"panel {i} is {h:.4f} mm, its budget is {panel_mm:.4f} mm")
    return True


def test_stack_height_checker_can_fail(stack_samples):
    """PLANT a violation of the height invariant and prove the checker rejects it."""
    tx = T.Transcript("GENE1", "chrI", "+", exons=[(1000, 1300), (1800, 2000), (3400, 4100)])
    reg = T.Region("chrI", 900, 4200)
    chrome = dict(model_mm=7.0, model_gap_mm=2.0, axis_mm=8.0, top_mm=3.0, reserved_mm=0.0)
    fig = plt.figure(figsize=(7.2, 3.0), layout="none")
    try:
        res = P.read_stack(fig, stack_samples, tx, region=reg, panel_mm=12.0, panel_gap_mm=2.5,
                           **chrome)
        check_stack_height(fig, res, n=4, panel_mm=12.0, panel_gap_mm=2.5, **chrome)  # clean

        # 1. the page is no longer the sum of the budgets
        fig.set_size_inches(7.2, fig.get_figheight() + 0.5)
        with pytest.raises(AssertionError, match="the budgets sum to"):
            check_stack_height(fig, res, n=4, panel_mm=12.0, panel_gap_mm=2.5, **chrome)

        # 2. a COMPENSATING pair: one panel grows by exactly what another loses, so the
        #    total is untouched. A sum-of-heights check would pass this; the per-panel one
        #    must not.
        fig.set_size_inches(7.2, fig.get_figheight() - 0.5)
        H_mm = fig.get_figheight() * MM
        d = 2.0 / H_mm
        for k, sign in ((0, +1), (1, -1)):
            pos = res["panels"][k].get_position()
            res["panels"][k].set_position([pos.x0, pos.y0, pos.width, pos.height + sign * d])
        total = sum(a.get_position().height * H_mm for a in res["panels"])
        assert total == pytest.approx(4 * 12.0, abs=1e-6)      # the SUM is still clean
        with pytest.raises(AssertionError, match="its budget is"):
            check_stack_height(fig, res, n=4, panel_mm=12.0, panel_gap_mm=2.5, **chrome)
    finally:
        plt.close(fig)


def test_shave_cannot_overflow_when_the_floors_exactly_fill_the_budget(stack_samples):
    """The tightest case for the shave path: `_reduce_k` leaves the maximum number of bands,
    whose floors fill the band region EXACTLY, and the caller pins a scale far too generous.

    The shave distributes the excess over each band's headroom above the floor, so it can only
    fail if `excess > total_headroom`, i.e. `nb * floor > room` -- which `_reduce_k` rules out
    by construction (it reduces `nb` until `nb * floor <= room`). This test is the proof made
    executable, because the failure mode would be a silent overflow at render time.
    """
    floor, gap = TOK.track_geometry()["band_floor_mm"], 0.4
    bands_mm = 4 * floor + 3 * gap                 # exactly four bands at the floor
    budget = bands_mm + 2.5 + 0.8
    keep = P.union_clusters(stack_samples, k=4)
    plan = P.plan_read_panel("tightest", stack_samples[0]["reads"], keep, budget_mm=budget,
                             scale=P.BandScale(0.5, floor), band_gap_mm=gap)
    assert len(plan.bands) == 4
    assert all(b.h_mm == pytest.approx(floor) for b in plan.bands)
    assert all(b.capped for b in plan.bands)       # the cut is DECLARED, not silent
    assert plan.used_mm == pytest.approx(budget, abs=1e-9)
    check_panel_budget(plan, floor_mm=floor)
