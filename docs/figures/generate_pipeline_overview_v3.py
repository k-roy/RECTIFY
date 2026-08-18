#!/usr/bin/env python3
"""
Hero pipeline_overview figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, Inter font stack, restrained PAL palette
  - Minimum text — viewer should grok the structure in <5 seconds
  - 3 chemistries converge → shared multi-aligner → correct (5'/intron/3')
    → rectified BAM → analyze + DESeq2
  - No success banners, no marketing framing, descriptive headers only
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 756
FONT = "Inter, Helvetica Neue, Helvetica, Arial, sans-serif"

PAL = dict(
    title   = "#1e293b",
    heading = "#475569",
    label   = "#64748b",
    muted   = "#94a3b8",
    border  = "#cbd5e1",
    divider = "#e2e8f0",
    bg      = "#ffffff",
    blue    = "#2563eb",
    blue_l  = "#dbeafe",
    teal    = "#0d9488",
    teal_l  = "#ccfbf1",
    orange  = "#d97706",
    orange_l= "#fed7aa",
    green   = "#059669",
    green_l = "#d1fae5",
    indigo  = "#4f46e5",
    indigo_l= "#e0e7ff",
    slate   = "#334155",
    slate_l = "#f1f5f9",
)

# Per-chemistry track colors — two families x two chemistries
CHEM = {
    "drs":    dict(color=PAL["blue"],   fill=PAL["blue_l"],   name="ONT DRS",           sub="direct RNA"),
    "cdna":   dict(color=PAL["teal"],   fill=PAL["teal_l"],   name="ONT PCR-cDNA",      sub="UMI, full-length"),
    "truseq": dict(color=PAL["green"],  fill=PAL["green_l"],  name="TruSeq / CORALL v2", sub="UMI, whole transcriptome"),
    "qsrev":  dict(color=PAL["orange"], fill=PAL["orange_l"], name="QuantSeq REV",      sub="3'-end, antisense"),
}

# Column centers for the four input lanes (two per family)
COL = {"drs": 137, "cdna": 293, "truseq": 467, "qsrev": 623}
CHIP_W = 150
CHIP_H = 38


# ─────────────────────────────────────────────────────────────────────────────
# Drawing helpers
# ─────────────────────────────────────────────────────────────────────────────

def svg_open(h):
    return (
        f'<?xml version="1.0" encoding="utf-8" ?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" height="{h}" '
        f'viewBox="0 0 {FIG_W} {h}">\n'
        f'<defs>'
        f'<marker id="arr" markerWidth="6" markerHeight="6" refX="6" refY="3" '
        f'orient="auto" markerUnits="userSpaceOnUse">'
        f'<path d="M0,0 L0,6 L6,3 z" fill="{PAL["label"]}"/></marker>'
        f'</defs>\n'
        f'<rect fill="{PAL["bg"]}" width="{FIG_W}" height="{h}"/>\n'
        f'<g font-family="{FONT}">'
    )

def svg_close():
    return "</g>\n</svg>"

def text(x, y, s, size=10, color=None, weight="400", anchor="start", letter_spacing=None):
    c = color or PAL["label"]
    ls = f' letter-spacing="{letter_spacing}"' if letter_spacing else ""
    return (f'<text fill="{c}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}"{ls} x="{x}" y="{y}">{s}</text>')

def section_head(x, y, s):
    return text(x, y, s, size=11, color=PAL["heading"], weight="700",
                letter_spacing="0.8")

def chip(cx, cy, w, h, label, color, fill, sub=None):
    """Rounded chip with label centered, optional subtitle below."""
    x = cx - w / 2
    y = cy - h / 2
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" '
        f'fill="{fill}" stroke="{color}" stroke-width="1.2"/>',
        text(cx, cy + 4, label, size=12, color=color, weight="600", anchor="middle"),
    ]
    if sub:
        parts.append(text(cx, cy + h / 2 + 13, sub, size=11,
                          color=PAL["muted"], anchor="middle"))
    return "\n".join(parts)

def shared_box(cx, cy, w, h, title, sub=None, title_y=None):
    """Wider neutral box for shared stages, with optional subtitle.

    By default the title sits near vertical center. Pass `title_y` (absolute
    y-coordinate) to override — needed when the box also contains inner pills
    that need clearance below the title.
    """
    x = cx - w / 2
    y = cy - h / 2
    if title_y is None:
        title_y = cy + (2 if sub else 4)
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" '
        f'fill="{PAL["slate_l"]}" stroke="{PAL["slate"]}" stroke-width="1.2"/>',
        text(cx, title_y, title, size=12,
             color=PAL["slate"], weight="700", anchor="middle"),
    ]
    if sub:
        parts.append(text(cx, title_y + 14, sub, size=11,
                          color=PAL["muted"], anchor="middle"))
    return "\n".join(parts)

def pill(cx, cy, w, h, label, color, fill):
    """Smaller secondary pill for sub-stages inside a parent box."""
    x = cx - w / 2
    y = cy - h / 2
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" '
        f'fill="{fill}" stroke="{color}" stroke-width="1"/>\n'
        + text(cx, cy + 4, label, size=12, color=color, weight="600", anchor="middle")
    )

def arrow(x1, y1, x2, y2, color=None):
    c = color or PAL["label"]
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{c}" '
            f'stroke-width="1.4" marker-end="url(#arr)"/>')

def hdivider(y, x1=60, x2=700, color=None):
    c = color or PAL["divider"]
    return f'<line x1="{x1}" x2="{x2}" y1="{y}" y2="{y}" stroke="{c}" stroke-width="1"/>'


# ─────────────────────────────────────────────────────────────────────────────
# Build the figure
# ─────────────────────────────────────────────────────────────────────────────

def build():
    out = [svg_open(FIG_H)]

    # Title block
    out.append(text(FIG_W / 2, 28, "RECTIFY", size=18,
                    color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W / 2, 48,
                    "one pipeline · ONT long reads + NGS short reads",
                    size=11.5, color=PAL["muted"], anchor="middle"))

    # ── INPUTS (two families: ONT long read · NGS short read) ────────────────
    y_inputs_head = 74
    out.append(section_head(60, y_inputs_head, "INPUTS"))

    # family group labels, each spanning its two lanes, with a light rule
    y_group = 94
    for label, gx1, gx2 in (("ONT LONG READ", 62, 368), ("NGS SHORT READ", 392, 698)):
        gcx = (gx1 + gx2) / 2
        out.append(text(gcx, y_group, label, size=10, color=PAL["heading"],
                        weight="700", anchor="middle", letter_spacing="0.8"))
        out.append(hdivider(y_group + 8, x1=gx1, x2=gx2))

    y_chip = 124
    for key in ("drs", "cdna", "truseq", "qsrev"):
        c = CHEM[key]
        out.append(chip(COL[key], y_chip, CHIP_W, CHIP_H,
                        c["name"], c["color"], c["fill"], sub=c["sub"]))

    # ── PRE-PROCESS (per-protocol labels, just below each chip's subtitle) ───
    y_pre_label = 178
    out.append(text(COL["drs"], y_pre_label, "trim-polya",
                    size=11, color=PAL["blue"], weight="600", anchor="middle"))
    out.append(text(COL["cdna"], y_pre_label, "correct-cdna",
                    size=11, color=PAL["teal"], weight="600", anchor="middle"))
    out.append(text(COL["truseq"], y_pre_label, "split --umi",
                    size=11, color=PAL["green"], weight="600", anchor="middle"))
    out.append(text(COL["qsrev"], y_pre_label, "(no pre-process)",
                    size=11, color=PAL["muted"], weight="400", anchor="middle"))
    out.append(text(FIG_W / 2, y_pre_label + 14,
                    "per-protocol — strip adapters / UMIs / poly(A)",
                    size=11, color=PAL["muted"], anchor="middle"))

    # ── LAYOUT GEOMETRY ──────────────────────────────────────────────────────
    # Single ARROW_GAP applied between every adjacent box → uniform visual
    # rhythm. Arrowhead tips land exactly on box top edges (marker refX=6,
    # userSpaceOnUse — see svg_open()).
    cx_funnel = FIG_W / 2
    ARROW_GAP = 24

    # Multi-aligner box
    aligner_h = 44
    y_aligner_top = 232
    y_aligner_bot = y_aligner_top + aligner_h          # 262
    y_aligner     = y_aligner_top + aligner_h / 2      # center, 240

    # Overhang resolver box (Station A — consumes the minimap2 arm)
    resolver_h = 32
    y_resolver_top = y_aligner_bot + ARROW_GAP         # 286
    y_resolver_bot = y_resolver_top + resolver_h       # 318
    y_resolver     = y_resolver_top + resolver_h / 2   # 302

    # Junction pool box
    jpool_h = 32
    y_jpool_top = y_resolver_bot + ARROW_GAP           # 342
    y_jpool_bot = y_jpool_top + jpool_h                # 374
    y_jpool     = y_jpool_top + jpool_h / 2            # 358

    # Correct box
    correct_h = 88
    y_correct_top = y_jpool_bot + ARROW_GAP            # 398
    y_correct_bot = y_correct_top + correct_h          # 486
    y_correct     = y_correct_top + correct_h / 2      # 442

    # Consensus box
    consensus_h = 32
    y_consensus_top = y_correct_bot + ARROW_GAP        # 510
    y_consensus_bot = y_consensus_top + consensus_h    # 542
    y_consensus     = y_consensus_top + consensus_h / 2  # 526

    # Rectified BAM box
    bam_h = 28
    y_bam_top = y_consensus_bot + ARROW_GAP            # 566
    y_bam_bot = y_bam_top + bam_h                      # 594
    y_bam     = y_bam_top + bam_h / 2                  # 580

    # ── CONVERGING ARROWS (3 inputs → multi-aligner) ─────────────────────────
    # Each arrow drops straight down on the column center and lands on the
    # aligner box top edge — uniform ARROW_GAP length matches downstream arrows.
    y_arrow_from = y_aligner_top - ARROW_GAP
    for key in ("drs", "cdna", "truseq", "qsrev"):
        out.append(arrow(COL[key], y_arrow_from, COL[key], y_aligner_top))

    # ── MULTI-ALIGNER (parallel, NOT consensus — consensus happens later) ────
    out.append(shared_box(cx_funnel, y_aligner, 520, aligner_h,
                          "rectify align — run aligners in parallel",
                          sub="minimap2 + uLTRA + deSALT · short reads: COMPASS (TruSeq) or bbmap + bwa (QuantSeq)"))

    # arrow down to overhang resolver
    out.append(arrow(cx_funnel, y_aligner_bot, cx_funnel, y_resolver_top))

    # ── OVERHANG RESOLVER (Station A — substitutes the minimap2 arm) ─────────
    out.append(shared_box(cx_funnel, y_resolver, 520, resolver_h,
                          "overhang resolver — re-places soft-clipped splice overhangs on the minimap2 arm",
                          sub=None))

    # arrow down to junction pool
    out.append(arrow(cx_funnel, y_resolver_bot, cx_funnel, y_jpool_top))

    # ── JUNCTION POOL ────────────────────────────────────────────────────────
    out.append(shared_box(cx_funnel, y_jpool, 420, jpool_h,
                          "junction pool — observed across aligners + annotated DB",
                          sub=None))

    # arrow down to correct
    out.append(arrow(cx_funnel, y_jpool_bot, cx_funnel, y_correct_top))

    # ── CORRECT (5' / introns / 3' sub-pills, applied per aligner) ───────────
    out.append(shared_box(cx_funnel, y_correct, 560, correct_h,
                          "rectify correct — applied per aligner (standardized scoring)",
                          sub=None,
                          title_y=y_correct - 22))
    # inner sub-pills below the title
    y_subpills = y_correct + 12
    sub_w = 130
    sub_h = 26
    sub_specs = [
        (cx_funnel - 165, "5' ends",  PAL["indigo"], PAL["indigo_l"]),
        (cx_funnel,       "introns",  PAL["indigo"], PAL["indigo_l"]),
        (cx_funnel + 165, "3' ends",  PAL["indigo"], PAL["indigo_l"]),
    ]
    for sx, label, color, fill in sub_specs:
        out.append(pill(sx, y_subpills, sub_w, sub_h, label, color, fill))
    # tiny footnote inside box
    out.append(text(cx_funnel, y_correct + 38,
                    "junction pool informs 5′ rescue + intron refinement",
                    size=11, color=PAL["muted"], anchor="middle"))

    # arrow down to consensus
    out.append(arrow(cx_funnel, y_correct_bot, cx_funnel, y_consensus_top))

    # ── CONSENSUS (pick best CORRECTED alignment per read) ───────────────────
    out.append(shared_box(cx_funnel, y_consensus, 380, consensus_h,
                          "consensus — pick best corrected alignment per read",
                          sub=None))

    # arrow down to rectified BAM
    out.append(arrow(cx_funnel, y_consensus_bot, cx_funnel, y_bam_top))

    # ── RECTIFIED BAM ────────────────────────────────────────────────────────
    out.append(shared_box(cx_funnel, y_bam, 220, bam_h,
                          "rectified BAM", sub=None))

    # ── RE-ALIGNER STATIONS B + C (post-stages on the rectified BAM) ─────────
    y_st = y_bam_bot + 30
    out.append(pill(cx_funnel - 168, y_st, 300, 26,
                    "Station B · triage re-align (--triage)",
                    PAL["indigo"], PAL["indigo_l"]))
    out.append(pill(cx_funnel + 168, y_st, 300, 26,
                    "Station C · pool-gate report (default)",
                    PAL["indigo"], PAL["indigo_l"]))

    # divider before downstream
    out.append(hdivider(y_st + 24))

    # ── ANALYZE + DESEQ2 (two rows of chips) ─────────────────────────────────
    y_down_head = 689
    out.append(section_head(60, y_down_head, "DOWNSTREAM"))

    # row 1: analyze outputs
    y_an = 709
    out.append(text(60, y_an + 3, "analyze", size=11,
                    color=PAL["heading"], weight="600"))
    analyze_items = [
        (240, "5'/3' clusters"),
        (380, "splice tables"),
        (530, "isoforms (cDNA)"),
    ]
    for cx, label in analyze_items:
        out.append(pill(cx, y_an, 130, 24, label, PAL["slate"], PAL["slate_l"]))

    # row 2: DESeq2
    y_de = 739
    out.append(text(60, y_de + 3, "DESeq2", size=11,
                    color=PAL["heading"], weight="600"))
    de_items = [
        (170, "gene"),
        (290, "3' cluster"),
        (425, "5' cluster*"),
        (570, "splice junction*"),
    ]
    for cx, label in de_items:
        out.append(pill(cx, y_de, 110, 24, label, PAL["slate"], PAL["slate_l"]))

    out.append(text(700, y_de + 3, "* planned",
                    size=10, color=PAL["muted"], anchor="end"))

    out.append(svg_close())
    return "\n".join(out)


def render(svg_str, name):
    svg_path = os.path.join(OUTDIR, f"{name}.svg")
    png_path = os.path.join(OUTDIR, f"{name}.png")
    with open(svg_path, "w") as f:
        f.write(svg_str)
    # Render at 2x for retina
    subprocess.run(
        ["cairosvg", svg_path, "-o", png_path, "--output-width", str(FIG_W * 2)],
        check=True,
    )
    print(f"wrote {svg_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    render(build(), "pipeline_overview")
