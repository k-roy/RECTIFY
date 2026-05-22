#!/usr/bin/env python3
"""
Adaptive valley clustering figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, SVG-native (matches the rest of the README figure set)
  - Histogram of 3' end counts, three clusters separated by valleys
  - Cluster boundaries as dashed verticals; peaks marked with small triangles
  - No yellow algorithm-explainer box; no colored background panels
"""

import os
import subprocess
import math

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 280
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
    red     = "#dc2626",
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


# Plot area
PLOT_X_LEFT  = 80
PLOT_X_RIGHT = 720
PLOT_Y_TOP   = 80
PLOT_Y_BOT   = 230
PLOT_W = PLOT_X_RIGHT - PLOT_X_LEFT
PLOT_H = PLOT_Y_BOT - PLOT_Y_TOP

X_MAX = 100   # bp
Y_MAX = 40    # counts

def x_to_px(bp):
    return PLOT_X_LEFT + (bp / X_MAX) * PLOT_W

def y_to_px(count):
    return PLOT_Y_BOT - (count / Y_MAX) * PLOT_H


def gauss(center, sigma, x):
    return math.exp(-((x - center) ** 2) / (2 * sigma ** 2))


def build():
    out = [svg_open(FIG_H)]

    # Title + subtitle
    out.append(text(FIG_W/2, 26, "Adaptive valley clustering",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "group nearby CPA sites by peaks-then-valleys; boundaries at valley positions",
                    size=11, color=PAL["muted"], anchor="middle"))

    # Three cluster definitions: (center_bp, sigma, peak_height, color, label)
    clusters = [
        (22, 3.0, 38, PAL["blue"],   "Cluster 1"),
        (50, 3.5, 26, PAL["orange"], "Cluster 2"),
        (80, 3.2, 16, PAL["green"],  "Cluster 3"),
    ]
    # Valley positions (between adjacent clusters) — picked as midpoints
    valleys = [38, 68]
    # Cluster boundaries (dashed verticals) — extends to plot edges
    boundaries = [12, 38, 68, 92]

    # ── Axes ────────────────────────────────────────────────────────────────
    # X axis baseline
    out.append(f'<line x1="{PLOT_X_LEFT}" x2="{PLOT_X_RIGHT}" '
               f'y1="{PLOT_Y_BOT}" y2="{PLOT_Y_BOT}" '
               f'stroke="{PAL["border"]}" stroke-width="1"/>')
    # X tick labels (every 20 bp)
    for bp in [0, 20, 40, 60, 80, 100]:
        x = x_to_px(bp)
        out.append(f'<line x1="{x}" x2="{x}" y1="{PLOT_Y_BOT}" y2="{PLOT_Y_BOT+4}" '
                   f'stroke="{PAL["border"]}" stroke-width="1"/>')
        out.append(text(x, PLOT_Y_BOT + 16, str(bp), size=10,
                        color=PAL["muted"], anchor="middle"))
    out.append(text((PLOT_X_LEFT + PLOT_X_RIGHT)/2, PLOT_Y_BOT + 32,
                    "Genomic position (bp)",
                    size=12, color=PAL["heading"], anchor="middle", weight="600"))

    # Y axis labels (just gridline marks, no axis line — keeps things minimal)
    for cnt in [10, 20, 30, 40]:
        y = y_to_px(cnt)
        out.append(f'<line x1="{PLOT_X_LEFT - 4}" x2="{PLOT_X_LEFT}" '
                   f'y1="{y}" y2="{y}" stroke="{PAL["border"]}" stroke-width="1"/>')
        out.append(text(PLOT_X_LEFT - 8, y + 3, str(cnt), size=10,
                        color=PAL["muted"], anchor="end"))

    # Y axis label (rotated)
    y_mid = (PLOT_Y_TOP + PLOT_Y_BOT) / 2.0
    out.append(
        f'<text fill="{PAL["heading"]}" font-size="12" font-weight="600" '
        f'text-anchor="middle" x="20" y="{y_mid:.2f}" '
        f'transform="rotate(-90 20 {y_mid:.2f})">Read count</text>'
    )

    # ── Cluster boundaries (dashed verticals) ───────────────────────────────
    for bp in boundaries:
        x = x_to_px(bp)
        out.append(f'<line x1="{x}" x2="{x}" y1="{PLOT_Y_TOP}" y2="{PLOT_Y_BOT}" '
                   f'stroke="{PAL["muted"]}" stroke-width="1" '
                   f'stroke-dasharray="3,3" opacity="0.55"/>')

    # ── Histogram bars ──────────────────────────────────────────────────────
    bar_w = 4
    for (center, sigma, peak, color, _) in clusters:
        # Generate bars from center - 4*sigma to center + 4*sigma at 1bp resolution
        bp_min = max(0, int(center - 4 * sigma))
        bp_max = min(X_MAX, int(center + 4 * sigma) + 1)
        for bp in range(bp_min, bp_max + 1):
            v = peak * gauss(center, sigma, bp)
            if v < 0.6:
                continue
            x = x_to_px(bp) - bar_w / 2
            h = (v / Y_MAX) * PLOT_H
            y = PLOT_Y_BOT - h
            out.append(f'<rect x="{x}" y="{y}" width="{bar_w}" height="{h}" '
                       f'fill="{color}" opacity="0.85"/>')

    # ── Peak markers (small inverted triangles above peaks) ─────────────────
    for (center, sigma, peak, color, label) in clusters:
        x = x_to_px(center)
        y = y_to_px(peak) - 14
        out.append(
            f'<polygon points="{x-5},{y-7} {x+5},{y-7} {x},{y}" '
            f'fill="{PAL["red"]}" stroke="{PAL["red"]}"/>'
        )
        # Cluster label above the triangle
        out.append(text(x, y - 14, label, size=11,
                        color=color, weight="700", anchor="middle"))

    # ── Valley markers (small upright triangles at the bottom) ──────────────
    for bp in valleys:
        x = x_to_px(bp)
        y = PLOT_Y_BOT - 2
        out.append(
            f'<polygon points="{x-4},{y} {x+4},{y} {x},{y-6}" '
            f'fill="{PAL["teal"]}" stroke="{PAL["teal"]}"/>'
        )

    # ── Small inline legend (bottom right) ──────────────────────────────────
    leg_x = PLOT_X_RIGHT - 230
    leg_y = 60
    # peak swatch
    out.append(f'<polygon points="{leg_x-5},{leg_y-3} {leg_x+5},{leg_y-3} {leg_x},{leg_y+4}" '
               f'fill="{PAL["red"]}"/>')
    out.append(text(leg_x + 9, leg_y + 4, "peak",
                    size=11, color=PAL["heading"]))
    out.append(f'<polygon points="{leg_x+57},{leg_y+4} {leg_x+67},{leg_y+4} '
               f'{leg_x+62},{leg_y-3}" fill="{PAL["teal"]}"/>')
    out.append(text(leg_x + 71, leg_y + 4, "valley",
                    size=11, color=PAL["heading"]))
    out.append(f'<line x1="{leg_x+118}" x2="{leg_x+138}" y1="{leg_y}" y2="{leg_y}" '
               f'stroke="{PAL["muted"]}" stroke-dasharray="3,3"/>')
    out.append(text(leg_x + 142, leg_y + 4, "cluster boundary",
                    size=11, color=PAL["heading"]))

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
    render(build(), "adaptive_clustering")
