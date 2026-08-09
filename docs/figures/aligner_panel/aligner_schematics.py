#!/usr/bin/env python3
"""RECTIFY aligner-panel explainer schematics (3 figures).

F1 mechanism panel (factual) | F2 'chimeric' disambiguation (factual) |
F3 mosaic-stitching PROPOSAL (hypothesis, clearly labelled).

Sources: memory project-aligner-orthogonality-panel (verified Opus panel survey),
corrected_consensus.py (_add_chimera_flag / winner selection). Schematic — no data.

README / publication-grade restyle: Okabe-Ito-derived colorblind-safe palette
(five distinguishable per-aligner colors), unified type scale, consistent
box / badge / arrow styling. Content is verbatim — aesthetics only.
"""
import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch

SCRIPT = "aligner_schematics.py"
DATE = datetime.date.today().isoformat()
PROV = ("Schematic (no data) · sources: memory project-aligner-orthogonality-panel "
        "(Opus panel survey, primary-source verified) + corrected_consensus.py "
        f"(_add_chimera_flag / winner selection) · {SCRIPT} · {DATE}")

# ── Typography ───────────────────────────────────────────────────────────────
# A single sans stack + a consistent size hierarchy keeps all three figures
# visually cohesive. Absolute (0–100) placement means sizes and geometry are
# tuned together; treat these as the canonical scale.
plt.rcParams.update({
    # DejaVu Sans is the matplotlib default and is the only stack member
    # guaranteed to carry every load-bearing glyph here (→ ↑ ′ − · ≥ × — “ ”);
    # Helvetica/Arial lack the arrows, so they would silently drop to a tofu box.
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica Neue", "Helvetica",
                        "Arial", "sans-serif"],
    "svg.fonttype": "path",          # vectorize glyphs → self-contained SVG that
                                     # renders the load-bearing → ↑ ′ − ≥ glyphs
                                     # identically everywhere (no installed-font dep)
    "pdf.fonttype": 42,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

FS = dict(title=15.5, subtitle=10, h=9.2, body=8.2, small=7.6,
          chip=10.0, tiny=7.0, foot=5.6)

# ── Palette (Okabe-Ito derived; colorblind-safe, 5 distinguishable aligners) ──
C = dict(
    # per-aligner identity colors — deliberately spread across hue/value so they
    # stay distinguishable in deuteranopia/protanopia and in greyscale.
    mm2="#0072B2",      # blue
    ultra="#009E73",    # bluish green
    desalt="#7B52AB",   # violet
    mpb="#E69F00",      # orange
    gmap="#D55E00",     # vermillion
    # semantic valence (paired ALWAYS with +/− glyphs & words, never color alone)
    good="#009E73", bad="#C0392B", seam="#138D75",
    accent="#5B4FA8",   # proposal / construction violet
    gold="#B8860B",     # callout text
    # neutrals
    grey="#5f6368", ink="#1f2329", mute="#55585c",
    box="#F7F7F5", boxalt="#FBFBF9", edge="#C9CCD1", read="#DBD6C9",
    line="#B7BBC2", rule="#9AA0A6",
    # soft tints for step / badge fills
    tint_blue="#E7EFF7", tint_green="#E6F3EC", tint_amber="#FBEFDC",
    tint_violet="#EFEBF8", tint_neutral="#ECEEF1",
)


def footer(fig):
    fig.text(0.5, 0.011, PROV, ha="center", va="bottom", fontsize=FS["foot"],
             color="0.45", wrap=True)


def rbox(ax, x, y, w, h, fc, ec=None, lw=1.0, r=0.02, z=1, alpha=1.0):
    if ec is None:
        ec = C["edge"]
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, fc=fc, ec=ec, lw=lw, alpha=alpha, zorder=z,
        boxstyle=f"round,pad=0.002,rounding_size={r}", mutation_aspect=1))


def header(ax, title, subtitle, ty=96.5, sy=92.6):
    ax.text(50, ty, title, ha="center", fontsize=FS["title"],
            fontweight="bold", color=C["ink"])
    ax.text(50, sy, subtitle, ha="center", fontsize=FS["subtitle"], color="0.38")


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — five aligners, four independent splice models
# ─────────────────────────────────────────────────────────────────────────────
def fig1(path, return_fig=False):
    fig, ax = plt.subplots(figsize=(13.4, 7.2))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    header(ax, "The five RECTIFY aligners: seeding lineage × splice model",
           "What each brings to the consensus — and where its risk lives")

    # column headers
    cols = [(10.5, "ALIGNER"), (31, "SEEDING / ALIGNMENT LINEAGE"),
            (57.5, "SPLICE MODEL"), (81, "BRINGS  ·  RISK")]
    for x, t in cols:
        ax.text(x, 87.6, t, ha="center", fontsize=FS["h"], fontweight="bold",
                color=C["mute"])
    ax.plot([4, 96], [85.9, 85.9], color=C["rule"], lw=1.1)

    rows = [
        ("minimap2", C["mm2"], "minimizer seed → chain → base-DP",
         "shared", "reference family; fast & robust", "splice kernel = the shared noncan=9"),
        ("uLTRA", C["ultra"], "maximal exact matches + annotation 2-pass",
         "independent", "small exons others miss", "needs a GTF; annotation-biased"),
        ("deSALT", C["desalt"], "de-Bruijn graph index → 2-pass refine",
         "independent", "high sensitivity; non-minimizer", "long-range false junctions"),
        ("mapPacBio", C["mpb"], "k-mer seed → affine-gap DP  (BBMap)",
         "independent", "clean HP-aware indels at 3′ end", "short-anchor over-call gaming"),
        ("GMAP", C["gmap"], "oligomer-hash → diagonalize → chain → own splice-DP",
         "independent", "novel lineage; #1 ONT exon-span", "short-anchor whole-read gaming"),
    ]
    y = 79.5; rh = 13.2
    for i, (name, col, lineage, smodel, brings, risk) in enumerate(rows):
        rowfc = C["box"] if i % 2 == 0 else C["boxalt"]
        rbox(ax, 4, y - rh + 1.2, 92, rh - 1.6, rowfc, ec=C["edge"], r=0.018)
        # colored accent rail at the left edge of the row (per-aligner identity)
        ax.add_patch(Rectangle((4, y - rh + 1.2), 0.7, rh - 1.6, fc=col, ec="none",
                     zorder=2))
        # name chip
        rbox(ax, 5.7, y - rh + 3.2, 10, rh - 5.4, col, ec=col, r=0.04, z=2)
        ax.text(10.7, y - rh/2 + 0.6, name, ha="center", va="center", color="white",
                fontsize=FS["chip"], fontweight="bold", zorder=3)
        # lineage
        ax.text(18.4, y - rh/2 + 0.6, lineage, ha="left", va="center",
                fontsize=FS["body"], color=C["ink"])
        # splice model badge
        if smodel == "independent":
            rbox(ax, 50.8, y - rh/2 - 1.7, 13.4, 4.4, C["tint_green"],
                 ec=C["good"], r=0.06, z=2)
            ax.text(57.5, y - rh/2 + 0.5, "INDEPENDENT", ha="center", va="center",
                    fontsize=FS["small"], fontweight="bold", color=C["good"], zorder=3)
        else:
            rbox(ax, 50.8, y - rh/2 - 1.7, 13.4, 4.4, C["tint_neutral"],
                 ec=C["rule"], r=0.06, z=2)
            ax.text(57.5, y - rh/2 + 0.5, "shared kernel", ha="center", va="center",
                    fontsize=FS["small"], fontweight="bold", color=C["mute"], zorder=3)
        # brings / risk
        ax.text(67, y - rh/2 + 2.5, "+ " + brings, ha="left", va="center",
                fontsize=FS["small"], color=C["good"])
        ax.text(67, y - rh/2 - 1.6, "− " + risk, ha="left", va="center",
                fontsize=FS["small"], color=C["bad"])
        y -= rh

    # bottom callout
    rbox(ax, 4, 4.4, 92, 7.8, "#FBF5E6", ec="#E3C97A", lw=1.3, r=0.02)
    ax.text(50, 9.9, "Splice-model diversity comes ONLY from the four INDEPENDENT models "
            "(uLTRA · deSALT · mapPacBio · GMAP).",
            ha="center", fontsize=FS["h"], fontweight="bold", color=C["gold"])
    ax.text(50, 6.4, "minimap2 and every derivative (winnowmap2, minisplice, graphmap2, gapmm2) "
            "share one splice kernel. GMAP's is the most novel lineage — the reason it is worth this much effort.",
            ha="center", fontsize=FS["body"], color="0.32")
    footer(fig)
    if return_fig:
        return fig
    _save(fig, path)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — two meanings of "chimeric"
# ─────────────────────────────────────────────────────────────────────────────
def fig2(path, return_fig=False):
    fig, ax = plt.subplots(figsize=(12.4, 6.6))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    header(ax, "“Chimeric” is being used in two opposite senses",
           "Same word, opposite valence — worth disambiguating before we design",
           ty=96, sy=91.6)

    # ---- Panel A : QC reject (current code) ----
    rbox(ax, 3, 30, 45, 56, "#FBEEEE", ec=C["bad"], lw=1.5, r=0.018)
    ax.text(25.5, 82.5, "A.  chimeric ALIGNMENT  =  QC reject", ha="center",
            fontsize=11.5, fontweight="bold", color=C["bad"])
    ax.text(25.5, 78.6, "current code: _add_chimera_flag()", ha="center",
            fontsize=FS["body"], color=C["mute"], style="italic")
    # read forced across two unrelated loci
    ax.text(8, 71.6, "read", ha="left", fontsize=FS["body"], color="0.3")
    rbox(ax, 8, 67.5, 34, 3.0, C["read"], ec="0.5", r=0.03, z=2)
    # genome region 1 and 2 (far apart, unrelated)
    rbox(ax, 8, 58, 13, 3.2, "#CDD7E5", ec="0.5", r=0.03, z=2)
    ax.text(14.5, 55.3, "locus X", ha="center", fontsize=FS["small"], color=C["mute"])
    rbox(ax, 29, 58, 13, 3.2, "#CDD7E5", ec="0.5", r=0.03, z=2)
    ax.text(35.5, 55.3, "locus Y (unrelated)", ha="center", fontsize=FS["small"], color=C["mute"])
    # the spurious bridging junction
    ax.add_patch(FancyArrowPatch((21, 67.5), (29, 61.2), connectionstyle="arc3,rad=-0.35",
                 arrowstyle="-", lw=1.7, color=C["bad"], linestyle=(0, (3, 2)), zorder=3))
    ax.text(25.5, 64.2, "spurious long\njunction", ha="center", fontsize=FS["tiny"],
            color=C["bad"])
    ax.text(25.5, 48.5, "An alignment that fuses UNRELATED genomic\nsequence (false N-op, "
            "short anchor, no support).", ha="center", fontsize=FS["body"], color="0.2")
    rbox(ax, 9, 35.5, 33, 8.5, "white", ec=C["bad"], r=0.03)
    ax.text(25.5, 41.3, "RECTIFY response:  FLAG  →  demote / drop", ha="center",
            fontsize=8.8, fontweight="bold", color=C["bad"])
    ax.text(25.5, 37.6, "standard genomics usage of “chimeric”", ha="center",
            fontsize=FS["small"], color=C["mute"], style="italic")

    # ---- Panel B : constructive stitching (user intent) ----
    rbox(ax, 52, 30, 45, 56, "#EAF6F0", ec=C["good"], lw=1.5, r=0.018)
    ax.text(74.5, 82.5, "B.  chimeric STITCHING  =  construction", ha="center",
            fontsize=11.5, fontweight="bold", color=C["good"])
    ax.text(74.5, 78.6, "your intent: build a better alignment from parts", ha="center",
            fontsize=FS["body"], color=C["mute"], style="italic")
    ax.text(57, 71.6, "read", ha="left", fontsize=FS["body"], color="0.3")
    # best segment from aligner A then aligner B, fused at a seam
    rbox(ax, 57, 67.5, 17, 3.0, C["mm2"], ec="0.35", r=0.03, z=2)
    rbox(ax, 74, 67.5, 17, 3.0, C["gmap"], ec="0.35", r=0.03, z=2)
    ax.text(65.5, 71.4, "aligner A best", ha="center", fontsize=FS["tiny"], color=C["mm2"])
    ax.text(82.5, 71.4, "aligner B best", ha="center", fontsize=FS["tiny"], color=C["gmap"])
    ax.plot([74, 74], [66, 72], color=C["seam"], lw=2.4, zorder=3)
    ax.text(74, 64.4, "seam\n(ref=query concordant)", ha="center", fontsize=FS["tiny"],
            color=C["seam"])
    # one aligned consensus below
    rbox(ax, 57, 58, 34, 3.2, "#CDE9DF", ec=C["good"], r=0.03, z=2)
    ax.text(74, 55.3, "one stitched RECTIFY alignment", ha="center", fontsize=FS["small"],
            color=C["good"])
    ax.text(74.5, 48.5, "Fuse each aligner’s BEST region into one\nmosaic alignment, "
            "joined at concordant seams.", ha="center", fontsize=FS["body"], color="0.2")
    rbox(ax, 58, 35.5, 33, 8.5, "white", ec=C["good"], r=0.03)
    ax.text(74.5, 41.3, "RECTIFY response:  BUILD a mosaic", ha="center", fontsize=8.8,
            fontweight="bold", color=C["good"])
    ax.text(74.5, 37.6, "“chimeric” in the literal made-of-parts sense", ha="center",
            fontsize=FS["small"], color=C["mute"], style="italic")

    # bottom proposal
    rbox(ax, 14, 6.3, 72, 16, C["tint_violet"], ec=C["accent"], lw=1.3, r=0.02)
    ax.text(50, 18.3, "Proposal to stop the collision:", ha="center", fontsize=9.8,
            fontweight="bold", color=C["accent"])
    ax.text(50, 13.4, "Keep “chimeric” for the QC reject (A) — it matches standard "
            "usage and the existing _add_chimera_flag.", ha="center", fontsize=FS["body"], color="0.2")
    ax.text(50, 9.6, "Name the construction (B)  “MOSAIC  /  STITCHED CONSENSUS”.  "
            "Recommended standard: mosaic consensus.", ha="center", fontsize=FS["body"],
            fontweight="bold", color=C["accent"])
    footer(fig)
    if return_fig:
        return fig
    _save(fig, path)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — mosaic stitching status + the three gaps
# ─────────────────────────────────────────────────────────────────────────────
def fig3(path, return_fig=False):
    fig, ax = plt.subplots(figsize=(13.2, 8.0))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    ax.text(50, 97.3, "Chimeric stitching is ALREADY BUILT — the pipeline, and the three gaps",
            ha="center", fontsize=14.8, fontweight="bold", color=C["ink"])
    rbox(ax, 24, 91.2, 52, 4.6, C["tint_green"], ec=C["good"], lw=1.3, r=0.04)
    ax.text(50, 93.5, "chimeric_consensus.py (1,211 lines, tested) · DEFAULT-ON in `rectify run` "
            "(--chimeric-consensus=True)", ha="center", fontsize=FS["body"],
            fontweight="bold", color="#15663A")

    # ---- the real pipeline (function-anchored) ----
    ax.text(7, 86, "PIPELINE", ha="left", fontsize=FS["h"], fontweight="bold", color=C["mute"])
    steps = [
        ("per-aligner\nBAMs", C["tint_neutral"]),
        ("find_sync_points\nALL aligners agree\n(query→ref)", C["tint_blue"]),
        ("identify_segments\n5′ / interior / 3′", C["tint_blue"]),
        ("score_segment\nposition-aware\n(motif + annot)", C["tint_amber"]),
        ("build_chimeric_cigar\nfuse + D/N bridge", C["tint_blue"]),
        ("validate /\nfallback", C["tint_blue"]),
        ("consensus BAM\nXz=1, Xq/Xw/Xy", C["tint_green"]),
    ]
    n = len(steps); x = 4.5; w = 11.8; gap = (92 - n * w) / (n - 1); yb = 76.5; hh = 7.8
    cx = []
    for i, (lab, fc) in enumerate(steps):
        rbox(ax, x, yb, w, hh, fc, ec="#9FB0C4", r=0.03, z=2)
        ax.text(x + w/2, yb + hh/2, lab, ha="center", va="center", fontsize=FS["tiny"],
                color=C["ink"], zorder=3)
        cx.append(x + w/2)
        if i < n - 1:
            ax.annotate("", xy=(x + w + gap*0.94, yb + hh/2), xytext=(x + w, yb + hh/2),
                        arrowprops=dict(arrowstyle="-|>", color=C["rule"], lw=1.3))
        x += w + gap
    ax.text(cx[3], yb - 1.6, "↑ the scoring step is where GMAP lives — and where the gaps bite",
            ha="center", fontsize=FS["tiny"], color=C["gold"])

    # ---- the three gaps ----
    ax.text(7, 66.8, "THE GAPS (what stands between the backbone and a working 5-panel mosaic)",
            ha="left", fontsize=FS["h"], fontweight="bold", color=C["mute"])
    gaps = [
        ("1 · SEAM is too strict",
         "find_sync_points needs ALL aligners to agree on (query→ref). With 5 aligners, one "
         "bad aligner breaks every seam → it falls back to single-best. FIX: relax to k-of-n "
         "(≥2) agreement so the full panel actually stitches.", C["tint_blue"], "#7D97B8"),
        ("2 · Scoring philosophy CLASH",
         "score_segment rewards canonical (+5) / annotated (+8), penalizes non-canonical (−3) — "
         "the OPPOSITE of the motif-agnostic winner path. It defends vs GMAP’s spurious splices "
         "but PENALIZES GMAP’s real non-canonical junctions; the annotated reward is DEAD in "
         "production (0 loaded — the intron-GTF loader bug). FIX: port the cross-read "
         "support-concordance gate (anchor-gate machinery) as a motif-agnostic defense.",
         C["tint_amber"], "#CF9A55"),
        ("3 · GMAP was never tested HERE",
         "The gmap6 / orthovalue rejection used merge_corrected_tsvs — the WINNER-TAKE-ALL path, "
         "not this one. GMAP has never been evaluated as a SEGMENT contributor. FIX: wire the "
         "dead ChimericStats; run the 5-panel chimeric path; measure GMAP’s per-segment wins "
         "(Xq/Xw/Xy) and whether they’re correct.", C["tint_violet"], C["accent"]),
    ]
    gy = 60.5; gh = 16.5
    for title, body, fc, ec in gaps:
        rbox(ax, 4, gy - gh, 92, gh - 1.4, fc, ec=ec, lw=1.2, r=0.016)
        ax.text(6.5, gy - 3.3, title, ha="left", fontsize=FS["h"], fontweight="bold", color=C["ink"])
        ax.text(6.5, gy - gh + 3.0, body, ha="left", va="bottom", fontsize=FS["small"],
                color="0.2", wrap=True)
        gy -= gh

    # bottom honest strip
    rbox(ax, 4, 1.8, 92, 6.2, C["box"], ec=C["rule"], lw=1.0, r=0.02)
    ax.text(50, 6.0, "SOLID: GMAP fails as a whole-read WINNER (3 tests, the winner-take-all path). "
            "OPEN: does its splice-DP help as a SEGMENT here?", ha="center", fontsize=FS["small"], color="0.2")
    ax.text(50, 3.2, "Different question — the rejection neither confirms nor refutes it. Only closing "
            "the gaps and measuring decides it.", ha="center", fontsize=FS["small"], color="0.35", style="italic")
    footer(fig)
    if return_fig:
        return fig
    _save(fig, path)


def _save(fig, path):
    """Save both a hi-dpi PNG and an SVG next to it (path is the .png)."""
    fig.savefig(path, dpi=220, bbox_inches="tight")
    if path.endswith(".png"):
        fig.savefig(path[:-4] + ".svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    import os
    outdir = os.path.dirname(os.path.abspath(__file__))
    fig1(os.path.join(outdir, "fig1_aligner_panel.png"))
    fig2(os.path.join(outdir, "fig2_chimeric_disambiguation.png"))
    fig3(os.path.join(outdir, "fig3_chimeric_stitching_built.png"))
    print("wrote PNG+SVG: fig1_aligner_panel, fig2_chimeric_disambiguation, "
          "fig3_chimeric_stitching_built to", outdir)
