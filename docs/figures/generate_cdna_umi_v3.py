#!/usr/bin/env python3
"""
ONT PCR-cDNA (PCB114.24) read architecture figure — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, single panel, one annotated read showing all segments
  - Segments: SSP · UMI (27 nt) · GGG bridge · mRNA · poly(A)
  - Small annotation under UMI shows the 4×(TT-VVVV)+TTT pattern
  - No three-class type table — that's a downstream detail
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 320
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
    green   = "#059669",   green_l  = "#d1fae5",
    orange  = "#d97706",   orange_l = "#fed7aa",
    indigo  = "#4f46e5",   indigo_l = "#e0e7ff",
    purple  = "#9333ea",   purple_l = "#f3e8ff",
)


def svg_open(h):
    return (
        f'<?xml version="1.0" encoding="utf-8" ?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" height="{h}" '
        f'viewBox="0 0 {FIG_W} {h}">\n'
        f'<rect fill="{PAL["bg"]}" width="{FIG_W}" height="{h}"/>\n'
        f'<g font-family="{FONT}">'
    )

def svg_close():
    return "</g>\n</svg>"

def text(x, y, s, size=10, color=None, weight="400", anchor="start"):
    c = color or PAL["label"]
    return (f'<text fill="{c}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" x="{x}" y="{y}">{s}</text>')

def block(x, y, w, h, label, fill, stroke, label_color):
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>\n'
        + text(x + w/2, y + h/2 + 4, label, size=10, color=label_color,
               weight="600", anchor="middle")
    )

def brace_below(x1, x2, y, color, label):
    """Brace opening downward, label centered below."""
    mid = (x1 + x2) / 2
    return (
        f'<path d="M{x1},{y} L{x1},{y+5} L{x2},{y+5} L{x2},{y}" '
        f'fill="none" stroke="{color}" stroke-width="1" opacity="0.7"/>\n'
        + text(mid, y + 17, label, size=8.5, color=color, anchor="middle",
               weight="600")
    )


def build():
    out = [svg_open(FIG_H)]

    out.append(text(FIG_W/2, 26, "ONT PCR-cDNA read architecture (SQK-PCB114.24)",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "both strands of each amplicon are sequenced — Dorado does not reverse-complement",
                    size=10, color=PAL["muted"], anchor="middle"))

    label_x = 80
    row_h = 30
    # Total content width (consistent across both rows)
    seg_w = dict(ssp=80, umi=90, ggg=32, mrna=260, pa=100)
    total_w = sum(seg_w.values())
    x_start = (FIG_W - total_w) / 2 + 12  # slight right-bias to leave room for 5′ label

    def draw_row(y, segs):
        """Draw a row of segments; returns list of (label, x_start, x_end)."""
        bounds = []
        x_cursor = x_start
        for (w, label, fill, stroke, lc) in segs:
            out.append(block(x_cursor, y, w, row_h, label, fill, stroke, lc))
            bounds.append((label, x_cursor, x_cursor + w))
            x_cursor += w
        return bounds

    # ── Row 1: orient = fwd ──────────────────────────────────────────────────
    y1 = 92
    out.append(text(label_x, y1 - 6, "orient = fwd   (≈50% of reads)",
                    size=10, color=PAL["teal"], weight="700"))
    out.append(text(label_x, y1 + row_h/2 + 4, "5′",
                    size=10, color=PAL["muted"], anchor="end"))

    fwd_segs = [
        (seg_w["ssp"],  "SSP",         PAL["blue_l"],   PAL["blue"],   PAL["blue"]),
        (seg_w["umi"],  "UMI (27 nt)", PAL["purple_l"], PAL["purple"], PAL["purple"]),
        (seg_w["ggg"],  "GGG",         PAL["orange_l"], PAL["orange"], PAL["orange"]),
        (seg_w["mrna"], "mRNA",        PAL["teal_l"],   PAL["teal"],   PAL["teal"]),
        (seg_w["pa"],   "poly(A)",     PAL["green_l"],  PAL["green"],  PAL["green"]),
    ]
    fwd_bounds = draw_row(y1, fwd_segs)
    out.append(text(x_start + total_w + 6, y1 + row_h/2 + 4, "3′",
                    size=10, color=PAL["muted"]))

    # Brace under the UMI showing the structured pattern (one annotation is enough)
    umi_x1, umi_x2 = fwd_bounds[1][1], fwd_bounds[1][2]
    out.append(brace_below(umi_x1, umi_x2, y1 + row_h + 4,
                           PAL["purple"], "(TT-VVVV) x 4 + TTT   V = A / C / G"))

    # ── Row 2: orient = rev (mirror layout) ──────────────────────────────────
    y2 = 220
    out.append(text(label_x, y2 - 6, "orient = rev   (≈50% of reads)",
                    size=10, color=PAL["orange"], weight="700"))
    out.append(text(label_x, y2 + row_h/2 + 4, "5′",
                    size=10, color=PAL["muted"], anchor="end"))

    rev_segs = [
        (seg_w["pa"],   "poly(T)",       PAL["green_l"],  PAL["green"],  PAL["green"]),
        (seg_w["mrna"], "mRNA (RC)",     PAL["teal_l"],   PAL["teal"],   PAL["teal"]),
        (seg_w["ggg"],  "CCC",           PAL["orange_l"], PAL["orange"], PAL["orange"]),
        (seg_w["umi"],  "UMI_RC (27 nt)", PAL["purple_l"], PAL["purple"], PAL["purple"]),
        (seg_w["ssp"],  "SSP_RC",        PAL["blue_l"],   PAL["blue"],   PAL["blue"]),
    ]
    rev_bounds = draw_row(y2, rev_segs)
    out.append(text(x_start + total_w + 6, y2 + row_h/2 + 4, "3′",
                    size=10, color=PAL["muted"]))

    # ── Footer note: correct-cdna normalises both ────────────────────────────
    out.append(text(FIG_W/2, 290,
                    "rectify correct-cdna strips SSP/UMI/GGG and poly(A/T) from both — output is one canonical orientation",
                    size=9, color=PAL["heading"], weight="600", anchor="middle"))

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
    render(build(), "cdna_umi_architecture")
