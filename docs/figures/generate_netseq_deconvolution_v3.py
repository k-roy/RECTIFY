#!/usr/bin/env python3
"""
NET-seq oligo(A) deconvolution figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, matches the README schematic figure set
  - Two side-by-side panels: A. Observed signal · B. Deconvolved signal
  - Observed signal: bars over an A-tract with downstream "spread" tail
  - Deconvolved signal: sharp single-base peaks at the true CPA positions
  - Inter font, restrained palette, descriptive headers, no marketing framing
"""

import os
import subprocess
import math

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
    red     = "#dc2626",   red_l    = "#fee2e2",
    purple  = "#9333ea",
    amber_l = "#fef3c7",   # soft A-tract background
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


# ── Shared plot geometry ────────────────────────────────────────────────────
X_MIN, X_MAX = 36, 80
Y_MAX = 30
PLOT_Y_TOP = 110
PLOT_Y_BOT = 260

# Two panels side-by-side
PANEL_L_X = 70
PANEL_L_W = 305
PANEL_R_X = PANEL_L_X + PANEL_L_W + 50   # 425
PANEL_R_W = 305

# Shared CPA positions and A-tract bounds
CPA_POSITIONS = [
    (51, "CPA 1", PAL["green"]),
    (55, "CPA 2", PAL["blue"]),
    (59, "CPA 3", PAL["red"]),
]
ATRACT_START = 50
ATRACT_END   = 62


def x_to_px(bp, panel_x, panel_w):
    return panel_x + (bp - X_MIN) / (X_MAX - X_MIN) * panel_w

def y_to_px(count):
    return PLOT_Y_BOT - (count / Y_MAX) * (PLOT_Y_BOT - PLOT_Y_TOP)


def panel_axes(panel_x, panel_w, header):
    """Axes for one panel: bottom line + Y ticks + X label, panel header above."""
    parts = []
    # Panel header
    parts.append(text(panel_x, PLOT_Y_TOP - 22,
                      header, size=12, color=PAL["title"], weight="700"))
    # X axis baseline
    parts.append(f'<line x1="{panel_x}" x2="{panel_x + panel_w}" '
                 f'y1="{PLOT_Y_BOT}" y2="{PLOT_Y_BOT}" '
                 f'stroke="{PAL["border"]}" stroke-width="1"/>')
    # X ticks at 40, 50, 60, 70, 80
    for bp in [40, 50, 60, 70, 80]:
        x = x_to_px(bp, panel_x, panel_w)
        parts.append(f'<line x1="{x}" x2="{x}" y1="{PLOT_Y_BOT}" y2="{PLOT_Y_BOT+4}" '
                     f'stroke="{PAL["border"]}" stroke-width="1"/>')
        parts.append(text(x, PLOT_Y_BOT + 16, str(bp), size=10,
                          color=PAL["muted"], anchor="middle"))
    # X axis label
    parts.append(text(panel_x + panel_w / 2, PLOT_Y_BOT + 34,
                      "Genomic position (bp)",
                      size=12, color=PAL["heading"], weight="600", anchor="middle"))
    # Y ticks at 10, 20, 30
    for cnt in [10, 20, 30]:
        y = y_to_px(cnt)
        parts.append(f'<line x1="{panel_x - 4}" x2="{panel_x}" '
                     f'y1="{y}" y2="{y}" stroke="{PAL["border"]}" stroke-width="1"/>')
        parts.append(text(panel_x - 8, y + 3, str(cnt), size=10,
                          color=PAL["muted"], anchor="end"))
    # Y axis label (rotated, on left margin of the panel)
    y_mid = (PLOT_Y_TOP + PLOT_Y_BOT) / 2.0
    parts.append(
        f'<text fill="{PAL["heading"]}" font-size="12" font-weight="600" '
        f'text-anchor="middle" x="{panel_x - 38}" y="{y_mid:.2f}" '
        f'transform="rotate(-90 {panel_x - 38} {y_mid:.2f})">Read counts</text>'
    )
    return "\n".join(parts)


def atract_shade(panel_x, panel_w):
    """Soft amber background block marking the genomic A-tract region."""
    x1 = x_to_px(ATRACT_START, panel_x, panel_w)
    x2 = x_to_px(ATRACT_END,   panel_x, panel_w)
    return (f'<rect x="{x1}" y="{PLOT_Y_TOP}" width="{x2 - x1}" '
            f'height="{PLOT_Y_BOT - PLOT_Y_TOP}" fill="{PAL["amber_l"]}" '
            f'opacity="0.55"/>')


def cpa_dashed_marker(panel_x, panel_w, bp, color):
    x = x_to_px(bp, panel_x, panel_w)
    return (f'<line x1="{x}" x2="{x}" y1="{PLOT_Y_TOP}" y2="{PLOT_Y_BOT}" '
            f'stroke="{color}" stroke-width="1" stroke-dasharray="3,3" '
            f'opacity="0.55"/>')


def bar(x_center, height_count, color, width=4, opacity=0.9):
    h = (height_count / Y_MAX) * (PLOT_Y_BOT - PLOT_Y_TOP)
    if h <= 0:
        return ""
    y = PLOT_Y_BOT - h
    return (f'<rect x="{x_center - width/2}" y="{y}" width="{width}" '
            f'height="{h}" fill="{color}" opacity="{opacity}"/>')


def build():
    out = [svg_open(FIG_H)]

    # ── Title + subtitle ────────────────────────────────────────────────────
    out.append(text(FIG_W / 2, 28, "NET-seq oligo(A) deconvolution",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W / 2, 48,
                    "NNLS deconvolution recovers sharp CPA peaks from the spread NET-seq signal over A-tracts",
                    size=11, color=PAL["muted"], anchor="middle"))

    # ── Panel A — Observed signal ────────────────────────────────────────────
    out.append(atract_shade(PANEL_L_X, PANEL_L_W))
    # dashed CPA markers
    for bp, _, color in CPA_POSITIONS:
        out.append(cpa_dashed_marker(PANEL_L_X, PANEL_L_W, bp, color))

    # observed bars: three peaks with exponential right-tail spread.
    # Smooth right-tailed distribution: peak * 0.54 at center, then
    # exponential decay with scale 3.5 to the right only.
    observed_peaks = [(51, 25), (55, 27), (59, 19)]
    total = {}
    for bp in range(X_MIN, X_MAX + 1):
        total[bp] = 0.0
        for c, amp in observed_peaks:
            if bp == c:
                total[bp] += amp * 0.54
            elif bp > c:
                d = bp - c
                total[bp] += amp * 0.46 * math.exp(-d / 3.5)
    # bonus singleton (control non-A-tract peak)
    total[41] += 1.5

    for bp, v in total.items():
        if v < 0.4:
            continue
        # color by region: inside A-tract = orange, downstream of A-tract = red,
        # upstream singleton control = blue
        if bp < ATRACT_START:
            color = PAL["blue"]
        elif bp <= ATRACT_END:
            color = PAL["orange"]
        else:
            color = PAL["red"]
        x = x_to_px(bp, PANEL_L_X, PANEL_L_W)
        out.append(bar(x, v, color))

    out.append(panel_axes(PANEL_L_X, PANEL_L_W, "A.  Observed signal"))

    # callout: "signal spreads downstream" pointing to the red tail
    callout_x = x_to_px(72, PANEL_L_X, PANEL_L_W)
    out.append(text(callout_x, 140, "signal spreads",
                    size=10.5, color=PAL["red"], weight="600", anchor="middle"))
    out.append(text(callout_x, 154, "downstream",
                    size=10.5, color=PAL["red"], weight="600", anchor="middle"))

    # A-tract label sits above the shaded region (not over the bars)
    atract_label_x = (x_to_px(ATRACT_START, PANEL_L_X, PANEL_L_W) +
                     x_to_px(ATRACT_END,   PANEL_L_X, PANEL_L_W)) / 2
    out.append(text(atract_label_x, PLOT_Y_TOP - 4, "A-tract",
                    size=10, color=PAL["orange"], weight="600", anchor="middle"))

    # ── Panel B — Deconvolved signal ─────────────────────────────────────────
    out.append(atract_shade(PANEL_R_X, PANEL_R_W))
    for bp, _, color in CPA_POSITIONS:
        out.append(cpa_dashed_marker(PANEL_R_X, PANEL_R_W, bp, color))

    # deconvolved peaks: sharp bars at the CPA positions, small residuals
    deconv = {51: 25, 52: 3, 55: 20, 56: 2.5, 59: 12, 60: 1.5}
    for bp, v in deconv.items():
        x = x_to_px(bp, PANEL_R_X, PANEL_R_W)
        # main peaks fully colored, residuals lighter
        is_main = bp in {51, 55, 59}
        out.append(bar(x, v, PAL["green"],
                       opacity=0.95 if is_main else 0.45))

    out.append(panel_axes(PANEL_R_X, PANEL_R_W, "B.  Deconvolved signal"))

    # CPA labels above the deconvolved peaks — each anchored at its own bar
    # height so labels don't crowd each other at the top
    cpa_label_heights = {51: 27, 55: 23, 59: 16}
    for bp, label, color in CPA_POSITIONS:
        x = x_to_px(bp, PANEL_R_X, PANEL_R_W)
        out.append(text(x, y_to_px(cpa_label_heights[bp]) - 4, label,
                        size=11, color=color, weight="700", anchor="middle"))

    # A-tract label above the shaded region (matches Panel A)
    atract_label_x_b = (x_to_px(ATRACT_START, PANEL_R_X, PANEL_R_W) +
                       x_to_px(ATRACT_END,   PANEL_R_X, PANEL_R_W)) / 2
    out.append(text(atract_label_x_b, PLOT_Y_TOP - 4, "A-tract",
                    size=10, color=PAL["orange"], weight="600", anchor="middle"))

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
    render(build(), "netseq_deconvolution")
