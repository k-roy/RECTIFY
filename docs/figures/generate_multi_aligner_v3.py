#!/usr/bin/env python3
"""
Multi-aligner consensus figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, single panel
  - Show the FORK / JUNCTION-POOL / PER-ALIGNER-CORRECT / MERGE structure
  - The key insight: correction happens PER ALIGNER (using a shared junction
    pool) BEFORE consensus, so the consensus pick compares apples to apples
    — every aligner's output graded on the same standardized scoring.
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 540
FONT = "Inter, Helvetica Neue, Helvetica, Arial, sans-serif"

PAL = dict(
    title   = "#1e293b",
    heading = "#475569",
    label   = "#64748b",
    muted   = "#94a3b8",
    border  = "#cbd5e1",
    divider = "#e2e8f0",
    bg      = "#ffffff",
    blue    = "#2563eb",   blue_l   = "#dbeafe",
    teal    = "#0d9488",   teal_l   = "#ccfbf1",
    orange  = "#d97706",   orange_l = "#fed7aa",
    green   = "#059669",   green_l  = "#d1fae5",
    indigo  = "#4f46e5",   indigo_l = "#e0e7ff",
    slate   = "#334155",   slate_l  = "#f1f5f9",
    amber   = "#b45309",   amber_l  = "#fef3c7",
)

# Per-aligner colors (consistent with hero figure)
ALIGNERS = [
    ("minimap2",  PAL["blue"],   PAL["blue_l"]),
    ("mapPacBio", PAL["teal"],   PAL["teal_l"]),
    ("gapmm2",    PAL["orange"], PAL["orange_l"]),
]


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

def text(x, y, s, size=10, color=None, weight="400", anchor="start"):
    c = color or PAL["label"]
    return (f'<text fill="{c}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" x="{x}" y="{y}">{s}</text>')

def chip(cx, cy, w, h, label, color, fill, sub=None):
    x, y = cx - w/2, cy - h/2
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" '
        f'fill="{fill}" stroke="{color}" stroke-width="1.2"/>',
        text(cx, cy + (1 if sub else 4), label, size=12, color=color,
             weight="600", anchor="middle"),
    ]
    if sub:
        parts.append(text(cx, cy + 14, sub, size=11, color=color,
                          anchor="middle"))
    return "\n".join(parts)

def arrow(x1, y1, x2, y2, color=None):
    c = color or PAL["label"]
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{c}" '
            f'stroke-width="1.4" marker-end="url(#arr)"/>')

def shared_bar(cx, cy, w, h, title, sub=None, fill=None, stroke=None):
    """Horizontal bar spanning multiple lanes (e.g., junction pool)."""
    f = fill or PAL["amber_l"]
    s = stroke or PAL["amber"]
    x, y = cx - w/2, cy - h/2
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" '
        f'fill="{f}" stroke="{s}" stroke-width="1.2"/>',
        text(cx, cy + (1 if sub else 4), title, size=12, color=s,
             weight="700", anchor="middle"),
    ]
    if sub:
        parts.append(text(cx, cy + 14, sub, size=10, color=s,
                          anchor="middle"))
    return "\n".join(parts)


def build():
    out = [svg_open(FIG_H)]

    out.append(text(FIG_W/2, 26, "Multi-aligner pipeline",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "correct each aligner first (standardized scoring), then pick the winner",
                    size=11, color=PAL["muted"], anchor="middle"))

    # ── Layout ───────────────────────────────────────────────────────────────
    cx = FIG_W / 2
    col_x = [180, 380, 580]   # three lanes
    COL_BOX_W = 160           # unified box width across aligner / correct columns
    y_input    = 80
    y_aligner  = 170   # +40px gap from FASTQ (was 130)
    y_pool     = 228   # maintained gap from aligner
    y_correct  = 310   # +30px extra gap below pool (was 252)
    y_consen   = 420   # maintained relative gap from correct
    y_bam      = 478   # gives 30px bottom margin at FIG_H=540

    # ── FASTQ at top ─────────────────────────────────────────────────────────
    out.append(chip(cx, y_input, 130, 32, "FASTQ", PAL["slate"], PAL["slate_l"]))

    # Fan-out arrows from FASTQ to each aligner
    for x in col_x:
        out.append(arrow(cx, y_input + 16, x, y_aligner - 16))

    # ── 3 aligner chips (same width as correct boxes below for clean columns)
    for (name, color, fill), x in zip(ALIGNERS, col_x):
        out.append(chip(x, y_aligner, COL_BOX_W, 30, name, color, fill, sub="raw BAM"))
        # arrow down to junction pool
        out.append(arrow(x, y_aligner + 15, x, y_pool - 17))

    # ── Junction pool (spans exactly the 3 column boxes) ─────────────────────
    pool_w = (col_x[-1] - col_x[0]) + COL_BOX_W
    out.append(shared_bar(cx, y_pool, pool_w, 38,
                          "junction pool",
                          sub="union of aligner junctions  +  annotated splice DB"))

    # Arrows from junction pool DOWN to each correct box
    for x in col_x:
        out.append(arrow(x, y_pool + 19, x, y_correct - 32))

    # ── 3 "correct" boxes (per-aligner) ──────────────────────────────────────
    for (name, color, fill), x in zip(ALIGNERS, col_x):
        # outer correct box — same width as the aligner chip directly above
        box_w, box_h = COL_BOX_W, 64
        bx = x - box_w/2
        by = y_correct - box_h/2
        out.append(
            f'<rect x="{bx}" y="{by}" width="{box_w}" height="{box_h}" rx="6" '
            f'fill="{PAL["indigo_l"]}" stroke="{PAL["indigo"]}" '
            f'stroke-width="1.2" opacity="0.95"/>'
        )
        out.append(text(x, y_correct - 16, "rectify correct",
                        size=12, color=PAL["indigo"], weight="700",
                        anchor="middle"))
        # mini sub-pills (5'/int/3') inline
        sub_y = y_correct + 8
        sub_labels = ["5′", "int", "3′"]
        sub_w_each = 36
        sub_gap = 4
        total_pills_w = len(sub_labels) * sub_w_each + (len(sub_labels) - 1) * sub_gap
        for i, lab in enumerate(sub_labels):
            sx = x - total_pills_w / 2 + i * (sub_w_each + sub_gap)
            out.append(
                f'<rect x="{sx}" y="{sub_y - 9}" width="{sub_w_each}" '
                f'height="18" rx="3" fill="{PAL["bg"]}" '
                f'stroke="{PAL["indigo"]}" stroke-width="0.8"/>'
            )
            out.append(text(sx + sub_w_each/2, sub_y + 4, lab,
                            size=11, color=PAL["indigo"], weight="600",
                            anchor="middle"))

        # arrow down to consensus
        out.append(arrow(x, y_correct + box_h/2, x, y_consen - 17))

    # ── Consensus (merge) ────────────────────────────────────────────────────
    out.append(shared_bar(cx, y_consen, 460, 34,
                          "consensus",
                          sub="pick best CORRECTED alignment per read",
                          fill=PAL["slate_l"], stroke=PAL["slate"]))

    # arrow down to BAM
    out.append(arrow(cx, y_consen + 17, cx, y_bam - 16))

    # ── Rectified BAM ────────────────────────────────────────────────────────
    out.append(chip(cx, y_bam, 220, 30, "rectified BAM",
                    PAL["slate"], PAL["slate_l"]))

    out.append(svg_close())
    return "\n".join(out)


def render(svg_str, name):
    svg_path = os.path.join(OUTDIR, f"{name}.svg")
    png_path = os.path.join(OUTDIR, f"{name}.png")
    with open(svg_path, "w") as f:
        f.write(svg_str)
    subprocess.run(
        ["cairosvg", svg_path, "-o", png_path, "--output-width", str(FIG_W * 2)],
        check=True,
    )
    print(f"wrote {svg_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    render(build(), "multi_aligner_consensus")
