#!/usr/bin/env python3
"""
HP scoring / empirical deletion-penalty figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, SVG-native (matches the rest of the README figure set)
  - Line chart: deletion penalty score vs homopolymer length (HP 1-12)
  - Two series: AT empirical (blue solid), CG empirical (orange dashed)
  - Shows penalty floor at HP≥9 (AT) and HP≥7 (CG)
  - CG HP≥9 rendered as hollow markers (low-count / floor region)
  - Values from H. sapiens GM12878 IVT RNA004 (match yeast R10.4.1 closely)
  - Endpoint value labels at HP=1 and HP=12 for each curve
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 420
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
    orange  = "#d97706",
)

# ── Data ─────────────────────────────────────────────────────────────────────
# Values from H. sapiens GM12878 IVT RNA004 (SQK-RNA004, dorado rna004_130bps_sup@v5.1.0)
# AT closely matches yeast R10.4.1; CG floor is at HP≥7.

HP = list(range(1, 13))  # 1..12

AT_VALS = {1: 0.46, 2: 0.38, 3: 0.30, 4: 0.18, 5: 0.13, 6: 0.075,
           7: 0.047, 8: 0.032, 9: 0.029, 10: 0.029, 11: 0.029, 12: 0.029}
CG_VALS = {1: 0.76, 2: 0.76, 3: 0.67, 4: 0.36, 5: 0.16, 6: 0.080,
           7: 0.053, 8: 0.053, 9: 0.053, 10: 0.053, 11: 0.053, 12: 0.053}
CG_EXTRAPOLATED = set(range(9, 13))  # hollow markers — low-count / floor region

# ── Plot geometry ─────────────────────────────────────────────────────────────

PLOT_L = 90
PLOT_R = 700
PLOT_T = 100
PLOT_B = 340
PLOT_W = PLOT_R - PLOT_L
PLOT_H_PX = PLOT_B - PLOT_T

Y_MAX = 1.0
Y_TICKS = [0, 0.25, 0.50, 0.75, 1.00]
X_TICKS = HP

def px_x(hp):
    return PLOT_L + (hp - 1) / 11.0 * PLOT_W

def px_y(v):
    return PLOT_B - (v / Y_MAX) * PLOT_H_PX


# ── SVG helpers ───────────────────────────────────────────────────────────────

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

def text(x, y, s, size=10, color=None, weight="400", anchor="start", transform=""):
    c = color or PAL["label"]
    t_attr = f' transform="{transform}"' if transform else ""
    return (f'<text fill="{c}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}" x="{x}" y="{y}"{t_attr}>{s}</text>')


# ── Line builders ─────────────────────────────────────────────────────────────

def line_path(vals_dict, color, stroke_width=2.0, dasharray=None):
    coords = [(px_x(hp), px_y(vals_dict[hp])) for hp in HP]
    d = "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in coords)
    dash = f' stroke-dasharray="{dasharray}"' if dasharray else ""
    return (f'<path d="{d}" fill="none" stroke="{color}" '
            f'stroke-width="{stroke_width}" stroke-linejoin="round"{dash}/>')

def circle_marker(hp, val, color, r=4.5, filled=True):
    x, y = px_x(hp), px_y(val)
    fill = color if filled else PAL["bg"]
    return (f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r}" '
            f'fill="{fill}" stroke="{color}" stroke-width="1.5"/>')

def square_marker(hp, val, color, size=7, filled=True):
    x, y = px_x(hp), px_y(val)
    half = size / 2.0
    fill = color if filled else PAL["bg"]
    return (f'<rect x="{x - half:.2f}" y="{y - half:.2f}" '
            f'width="{size}" height="{size}" '
            f'fill="{fill}" stroke="{color}" stroke-width="1.5"/>')


# ── Main build ────────────────────────────────────────────────────────────────

def build():
    out = [svg_open(FIG_H)]

    # ── Title + subtitle ─────────────────────────────────────────────────────
    out.append(text(FIG_W / 2, 28,
                    "Empirical Nanopore deletion penalty vs homopolymer length",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W / 2, 48,
                    "S. cerevisiae (R10.4.1) + H. sapiens (GM12878 RNA004) — "
                    "penalty floor at HP ≥ 9 (AT) and HP ≥ 7 (CG)",
                    size=11, color=PAL["muted"], anchor="middle"))

    # ── Horizontal gridlines ─────────────────────────────────────────────────
    for v in Y_TICKS:
        y = px_y(v)
        out.append(f'<line x1="{PLOT_L}" x2="{PLOT_R}" y1="{y:.2f}" y2="{y:.2f}" '
                   f'stroke="{PAL["divider"]}" stroke-width="0.5"/>')

    # ── Axes ─────────────────────────────────────────────────────────────────
    out.append(f'<line x1="{PLOT_L}" x2="{PLOT_R}" '
               f'y1="{PLOT_B}" y2="{PLOT_B}" '
               f'stroke="{PAL["border"]}" stroke-width="1"/>')
    out.append(f'<line x1="{PLOT_L}" x2="{PLOT_L}" '
               f'y1="{PLOT_T}" y2="{PLOT_B}" '
               f'stroke="{PAL["border"]}" stroke-width="1"/>')

    # X tick labels
    for hp in X_TICKS:
        x = px_x(hp)
        out.append(f'<line x1="{x:.2f}" x2="{x:.2f}" '
                   f'y1="{PLOT_B}" y2="{PLOT_B + 4}" '
                   f'stroke="{PAL["border"]}" stroke-width="1"/>')
        out.append(text(x, PLOT_B + 17, str(hp), size=10,
                        color=PAL["muted"], anchor="middle"))

    out.append(text((PLOT_L + PLOT_R) / 2, PLOT_B + 36,
                    "Homopolymer length (HP)",
                    size=12, color=PAL["heading"], weight="600", anchor="middle"))

    # Y tick labels
    for v in Y_TICKS:
        y = px_y(v)
        out.append(f'<line x1="{PLOT_L - 4}" x2="{PLOT_L}" '
                   f'y1="{y:.2f}" y2="{y:.2f}" '
                   f'stroke="{PAL["border"]}" stroke-width="1"/>')
        out.append(text(PLOT_L - 8, y + 3.5, f"{v:.2f}", size=10,
                        color=PAL["muted"], anchor="end"))

    # Y axis label (rotated)
    y_mid = (PLOT_T + PLOT_B) / 2.0
    out.append(text(22, y_mid, "Deletion penalty score",
                    size=12, color=PAL["heading"], weight="600", anchor="middle",
                    transform=f"rotate(-90 22 {y_mid:.2f})"))

    # ── Data lines ───────────────────────────────────────────────────────────
    # CG empirical (orange dashed, drawn first so AT sits on top)
    out.append(line_path(CG_VALS, PAL["orange"], stroke_width=2.0, dasharray="6,3"))
    # AT empirical (blue solid)
    out.append(line_path(AT_VALS, PAL["blue"], stroke_width=2.2))

    # ── Data point markers ───────────────────────────────────────────────────
    # AT: filled circles
    for hp in HP:
        out.append(circle_marker(hp, AT_VALS[hp], PAL["blue"], filled=True))

    # CG: filled squares for HP=1-8, hollow for HP=9-12 (low-count / floor region)
    for hp in HP:
        filled = hp not in CG_EXTRAPOLATED
        out.append(square_marker(hp, CG_VALS[hp], PAL["orange"], filled=filled))

    # ── Endpoint value labels ─────────────────────────────────────────────────
    # AT HP=1 label (above-right of the point)
    x1_at, y1_at = px_x(1), px_y(AT_VALS[1])
    out.append(text(x1_at + 8, y1_at - 7, "0.46",
                    size=10, color=PAL["blue"], weight="600"))

    # AT HP=12 floor label (above-left of the last point)
    x12_at, y12_at = px_x(12), px_y(AT_VALS[12])
    out.append(text(x12_at - 8, y12_at - 9, "0.029",
                    size=10, color=PAL["blue"], weight="600", anchor="end"))

    # CG HP=1 label (above-right, offset to avoid overlap with AT label)
    x1_cg, y1_cg = px_x(1), px_y(CG_VALS[1])
    out.append(text(x1_cg + 8, y1_cg - 7, "0.76",
                    size=10, color=PAL["orange"], weight="600"))

    # CG HP=12 floor label (above-right of the last point, slightly below AT label)
    x12_cg, y12_cg = px_x(12), px_y(CG_VALS[12])
    out.append(text(x12_cg - 8, y12_cg + 17, "0.053†",
                    size=10, color=PAL["orange"], weight="600", anchor="end"))

    # ── Floor annotation lines ────────────────────────────────────────────────
    # AT floor (HP ≥ 9)
    y_at_floor = px_y(AT_VALS[9])
    x_floor_at_start = px_x(9)
    out.append(f'<line x1="{x_floor_at_start:.2f}" x2="{PLOT_R}" '
               f'y1="{y_at_floor:.2f}" y2="{y_at_floor:.2f}" '
               f'stroke="{PAL["blue"]}" stroke-width="0.8" '
               f'stroke-dasharray="3,4" opacity="0.4"/>')

    # CG floor (HP ≥ 7)
    y_cg_floor = px_y(CG_VALS[7])
    x_floor_cg_start = px_x(7)
    out.append(f'<line x1="{x_floor_cg_start:.2f}" x2="{PLOT_R}" '
               f'y1="{y_cg_floor:.2f}" y2="{y_cg_floor:.2f}" '
               f'stroke="{PAL["orange"]}" stroke-width="0.8" '
               f'stroke-dasharray="3,4" opacity="0.4"/>')

    # Extrapolated footnote
    out.append(text((PLOT_L + PLOT_R) / 2, PLOT_B + 58,
                    "† floor / low-count region (HP ≥ 9 CG, n ≤ 3,500 positions)",
                    size=10, color=PAL["muted"], anchor="middle"))

    # ── Legend ───────────────────────────────────────────────────────────────
    LEG_X = PLOT_L + PLOT_W * 0.52
    LEG_Y = PLOT_T + 14
    row_h = 20
    swatch_w = 28

    # Row 1: AT empirical
    out.append(f'<line x1="{LEG_X:.0f}" x2="{LEG_X + swatch_w:.0f}" '
               f'y1="{LEG_Y:.0f}" y2="{LEG_Y:.0f}" '
               f'stroke="{PAL["blue"]}" stroke-width="2.2"/>')
    out.append(f'<circle cx="{LEG_X + swatch_w/2:.0f}" cy="{LEG_Y:.0f}" r="4" '
               f'fill="{PAL["blue"]}" stroke="{PAL["blue"]}" stroke-width="1.5"/>')
    out.append(text(LEG_X + swatch_w + 7, LEG_Y + 4,
                    "AT bases (A, T)", size=10, color=PAL["heading"]))

    # Row 2: CG empirical
    ry2 = LEG_Y + row_h
    out.append(f'<line x1="{LEG_X:.0f}" x2="{LEG_X + swatch_w:.0f}" '
               f'y1="{ry2:.0f}" y2="{ry2:.0f}" '
               f'stroke="{PAL["orange"]}" stroke-width="2" '
               f'stroke-dasharray="6,3"/>')
    half = 3.5
    out.append(f'<rect x="{LEG_X + swatch_w/2 - half:.0f}" y="{ry2 - half:.0f}" '
               f'width="7" height="7" '
               f'fill="{PAL["orange"]}" stroke="{PAL["orange"]}" stroke-width="1.5"/>')
    out.append(text(LEG_X + swatch_w + 7, ry2 + 4,
                    "CG bases (C, G)  ·  hollow = floor / low-count (HP ≥ 9)",
                    size=11, color=PAL["heading"]))

    out.append(svg_close())
    return "\n".join(out)


def render(svg_str, name):
    svg_path = os.path.join(OUTDIR, f"{name}.svg")
    png_path = os.path.join(OUTDIR, f"{name}.png")
    with open(svg_path, "w") as f:
        f.write(svg_str)
    subprocess.run(
        ["cairosvg", svg_path, "-o", png_path,
         "--output-width", str(FIG_W * 2)],
        check=True,
    )
    print(f"wrote {svg_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    render(build(), "hp_scoring")
