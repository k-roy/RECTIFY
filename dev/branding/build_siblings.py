"""Sibling logo concepts B, C, D for RECTIFY.

All three use the same palette tokens as build_logo.py / dev/figures/STYLE_GUIDE.md
so they sit side-by-side with concept A. Deliverables are light-mode square
mark + banner per concept; dark + 4x PNG raster gets generated only for the
concept the user picks to advance.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

OUT = Path(__file__).parent

# ---- shared design tokens (light only for siblings) ------------------------

PAL = dict(
    bg="#ffffff",
    baseline="#cbd5e1",
    text="#1e293b",
    text_muted="#475569",
    label="#64748b",       # slate-500
    muted="#94a3b8",       # slate-400
    spike="#0d9488",       # teal-600 — corrected
    spike_dark="#0f766e",  # teal-700
    # nucleotide fills (from STYLE_GUIDE.md)
    A_fill="#dcfce7", A_text="#166534",
    T_fill="#ffedd5", T_text="#9a3412",
    G_fill="#dbeafe", G_text="#1e40af",
    C_fill="#f3e8ff", C_text="#6b21a8",
    box_border="#cbd5e1",
)

FONT = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif"

# =========================================================================
# Concept B: Walkback chevrons
# =========================================================================
# 10 nucleotide boxes (CAGT | AAAAAA) centered horizontally, with a teal
# vertical CPA boundary at the T/A interface, a leftward walkback arrow
# above the A run, and a teal anchor caret below the boundary.

B_BASES = ["C", "A", "G", "T", "A", "A", "A", "A", "A", "A"]
B_CPA_INDEX = 4               # boundary sits between bases[3] and bases[4]
B_BOX_W, B_BOX_H, B_GAP, B_R = 16, 22, 3, 3
B_ROW_Y_TOP = 132             # box top y
B_ROW_LEN = len(B_BASES) * B_BOX_W + (len(B_BASES) - 1) * B_GAP
B_ROW_X0 = (256 - B_ROW_LEN) // 2  # center the row in 256-wide canvas


def _box(x: float, y: float, base: str) -> str:
    fill = PAL[f"{base}_fill"]
    text = PAL[f"{base}_text"]
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{B_BOX_W}" height="{B_BOX_H}" '
        f'rx="{B_R}" ry="{B_R}" fill="{fill}" stroke="{PAL["box_border"]}" '
        f'stroke-width="1" />'
        f'<text x="{x + B_BOX_W/2:.1f}" y="{y + B_BOX_H/2 + 4:.1f}" '
        f'text-anchor="middle" font-family="{FONT}" font-size="12" '
        f'font-weight="600" fill="{text}">{base}</text>'
    )


def mark_B(*, include_bg: bool = True) -> str:
    boxes = []
    for i, base in enumerate(B_BASES):
        x = B_ROW_X0 + i * (B_BOX_W + B_GAP)
        boxes.append(_box(x, B_ROW_Y_TOP, base))
    # CPA boundary: vertical teal line in the gap before bases[B_CPA_INDEX]
    cpa_x = B_ROW_X0 + B_CPA_INDEX * (B_BOX_W + B_GAP) - B_GAP / 2
    cpa_line = (
        f'<line x1="{cpa_x:.1f}" y1="{B_ROW_Y_TOP - 8}" '
        f'x2="{cpa_x:.1f}" y2="{B_ROW_Y_TOP + B_BOX_H + 8}" '
        f'stroke="{PAL["spike"]}" stroke-width="2.6" stroke-linecap="round" />'
    )
    # Walkback arrow above the A run, pointing left to the CPA
    a_run_x0 = B_ROW_X0 + B_CPA_INDEX * (B_BOX_W + B_GAP)
    a_run_x1 = B_ROW_X0 + B_ROW_LEN
    arrow_y = B_ROW_Y_TOP - 22
    arrow_path = (
        f'<defs><marker id="arrowL" markerWidth="9" markerHeight="9" '
        f'refX="2" refY="4.5" orient="auto" markerUnits="strokeWidth">'
        f'<path d="M 9,1 L 1,4.5 L 9,8 Z" fill="{PAL["label"]}" /></marker></defs>'
        f'<line x1="{a_run_x1 - 4:.1f}" y1="{arrow_y}" '
        f'x2="{a_run_x0 + 2:.1f}" y2="{arrow_y}" '
        f'stroke="{PAL["label"]}" stroke-width="2.4" stroke-linecap="round" '
        f'marker-end="url(#arrowL)" />'
    )
    # Small "walkback" caption ghost above the arrow? Skip for cleanness.
    # Anchor caret below the boundary
    anchor = (
        f'<path d="M {cpa_x-5:.1f},{B_ROW_Y_TOP + B_BOX_H + 14} '
        f'L {cpa_x:.1f},{B_ROW_Y_TOP + B_BOX_H + 24} '
        f'L {cpa_x+5:.1f},{B_ROW_Y_TOP + B_BOX_H + 14} Z" '
        f'fill="{PAL["spike_dark"]}" />'
    )
    bg_rect = (
        f'<rect width="256" height="256" fill="{PAL["bg"]}" />'
        if include_bg else ""
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 256 256"
     role="img" aria-label="RECTIFY walkback">
  {bg_rect}
  {arrow_path}
  {"".join(boxes)}
  {cpa_line}
  {anchor}
</svg>
'''


# =========================================================================
# Concept C: Anchored strand
# =========================================================================
# A wavy RNA backbone on the left resolves into a poly-A tail on the right.
# A teal vertical pin (with a teal dot at the top, like a map pin) marks the
# CPA at the transition. Six small "A" letters annotate the poly-A tail.

def mark_C(*, include_bg: bool = True) -> str:
    # Wavy backbone: sinusoidal path from (24, 128) to (132, 128)
    # Using cubic Beziers for smoothness.
    backbone = (
        "M 24,128 "
        "C 36,108 48,148 60,128 "
        "C 72,108 84,148 96,128 "
        "C 108,108 120,148 132,128"
    )
    # Straight tail line from (132, 128) to (224, 128)
    tail_line = (
        f'<line x1="132" y1="128" x2="224" y2="128" '
        f'stroke="{PAL["label"]}" stroke-width="2.6" stroke-linecap="round" />'
    )
    # 6 A letters along the tail
    a_positions = [148, 162, 176, 190, 204, 218]
    a_labels = "".join(
        f'<text x="{x}" y="116" text-anchor="middle" '
        f'font-family="{FONT}" font-size="13" font-weight="600" '
        f'fill="{PAL["A_text"]}">A</text>'
        for x in a_positions
    )
    # Faint A backgrounds (small rounded rectangles behind each A letter)
    a_bgs = "".join(
        f'<rect x="{x-6}" y="106" width="12" height="14" rx="2.5" '
        f'fill="{PAL["A_fill"]}" opacity="0.85" />'
        for x in a_positions
    )
    # CPA pin: vertical teal line + filled circle at top + anchor below
    pin = f'''
  <line x1="132" y1="92" x2="132" y2="142" stroke="{PAL["spike"]}" stroke-width="3" stroke-linecap="round" />
  <circle cx="132" cy="88" r="5.5" fill="{PAL["spike"]}" />
  <path d="M 127,148 L 132,160 L 137,148 Z" fill="{PAL["spike_dark"]}" />
'''
    bg_rect = (
        f'<rect width="256" height="256" fill="{PAL["bg"]}" />'
        if include_bg else ""
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 256 256"
     role="img" aria-label="RECTIFY anchored strand">
  {bg_rect}
  <path d="{backbone}" fill="none"
        stroke="{PAL["label"]}" stroke-width="2.6"
        stroke-linejoin="round" stroke-linecap="round" />
  {tail_line}
  {a_bgs}
  {a_labels}
  {pin}
</svg>
'''


# =========================================================================
# Concept D: Typographic monogram (I-as-spike)
# =========================================================================
# Wordmark-led: RECTIFY with the "I" replaced by a teal triangular spike
# carrying the visual DNA of concept A. Banner-only (square is awkward
# for typography); for square use, prefer A/B/C.

def banner_D() -> str:
    # Spike geometry: centered at x=cx, base width ~22, height ~160
    cx = 420
    base_y = 200
    apex_y = 36
    base_half = 14
    spike = (
        f'<path d="M {cx-base_half},{base_y} L {cx},{apex_y} '
        f'L {cx+base_half},{base_y} Z" '
        f'fill="{PAL["spike"]}" stroke="{PAL["spike_dark"]}" stroke-width="2" '
        f'stroke-linejoin="round" />'
    )
    anchor = (
        f'<path d="M {cx-5},{base_y+8} L {cx},{base_y+20} '
        f'L {cx+5},{base_y+8} Z" fill="{PAL["spike_dark"]}" />'
    )
    # Left text "RECT" ends just left of the spike; right text "FY" starts
    # just to the right. Use text-anchor to position around the spike.
    left_x = cx - base_half - 18
    right_x = cx + base_half + 18
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 840 256"
     role="img" aria-label="RECTIFY">
  <rect width="840" height="256" fill="{PAL["bg"]}" />
  <g font-family="{FONT}" font-weight="700" letter-spacing="3">
    <text x="{left_x}" y="178" text-anchor="end" font-size="140"
          fill="{PAL["text"]}">RECT</text>
    <text x="{right_x}" y="178" text-anchor="start" font-size="140"
          fill="{PAL["text"]}">FY</text>
  </g>
  {spike}
  {anchor}
  <text x="420" y="232" text-anchor="middle"
        font-family="{FONT}" font-size="16" font-weight="500"
        letter-spacing="3" fill="{PAL["text_muted"]}">NANOPORE  RNA  ALIGNMENT  CORRECTION</text>
</svg>
'''


# =========================================================================
# Banner wrappers for B and C (reuse the mark)
# =========================================================================

def banner_from_mark(mark_fn, label: str = "RECTIFY") -> str:
    inner_svg = mark_fn(include_bg=False)
    inner = inner_svg.split(">", 1)[1].rsplit("</svg>", 1)[0]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 840 256"
     role="img" aria-label="{label}">
  <rect width="840" height="256" fill="{PAL["bg"]}" />
  <svg x="0" y="0" width="256" height="256" viewBox="0 0 256 256">{inner}</svg>
  <g font-family="{FONT}">
    <text x="284" y="158" font-size="92" font-weight="700" letter-spacing="4"
          fill="{PAL["text"]}">RECTIFY</text>
    <text x="288" y="192" font-size="18" font-weight="500" letter-spacing="2.4"
          fill="{PAL["text_muted"]}">NANOPORE  RNA  ALIGNMENT  CORRECTION</text>
  </g>
</svg>
'''


# =========================================================================
# Build + rasterize
# =========================================================================

def build() -> list[Path]:
    pairs = {
        "rectify_mark_B_walkback.svg":     mark_B(),
        "rectify_banner_B_walkback.svg":   banner_from_mark(mark_B, "RECTIFY walkback"),
        "rectify_mark_C_anchored.svg":     mark_C(),
        "rectify_banner_C_anchored.svg":   banner_from_mark(mark_C, "RECTIFY anchored strand"),
        "rectify_banner_D_typographic.svg": banner_D(),
    }
    paths: list[Path] = []
    for name, content in pairs.items():
        p = OUT / name
        p.write_text(content)
        paths.append(p)
    return paths


def rasterize(svgs: list[Path]) -> list[Path]:
    pngs: list[Path] = []
    for svg in svgs:
        if "banner" in svg.name:
            w, h = 3360, 1024
        else:
            w, h = 1024, 1024
        png = svg.with_suffix(".png")
        subprocess.run(
            ["rsvg-convert", "-w", str(w), "-h", str(h), "-o", str(png), str(svg)],
            check=True,
        )
        pngs.append(png)
    return pngs


if __name__ == "__main__":
    svgs = build()
    pngs = rasterize(svgs)
    for p in svgs + pngs:
        print(f"  {p.name}  ({p.stat().st_size} bytes)")
