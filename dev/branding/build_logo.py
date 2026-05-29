"""Generate the RECTIFY logo set.

Concept: Signal rectifier. A noisy oscillation on the left (the smeared 3'-end
density of uncorrected reads, scattered across genomic A-tracts) resolves into
a single tall, sharp teal spike on the right (the rectified CPA — one true
cleavage site, anchored). Lab figure-style tokens from
rectify/dev/figures/STYLE_GUIDE.md.

Outputs (in this directory):
  rectify_mark_light.svg     square mark, light background  (256x256)
  rectify_mark_dark.svg      square mark, dark background   (256x256)
  rectify_banner_light.svg   horizontal lockup, light       (840x256)
  rectify_banner_dark.svg    horizontal lockup, dark        (840x256)
  *.png                      4x raster exports of each
"""

from __future__ import annotations

import subprocess
from pathlib import Path

OUT = Path(__file__).parent

# ---- design tokens (mirror dev/figures STYLE_GUIDE.md) --------------------

PAL_LIGHT = dict(
    bg="#ffffff",
    baseline="#cbd5e1",     # border
    noise="#64748b",        # label / slate-500 — the muted, uncorrected signal
    noise_soft="#94a3b8",   # muted / slate-400 — secondary noise echo
    spike="#0d9488",        # teal-600 — rectified / corrected
    spike_dark="#0f766e",   # teal-700 — spike outline
    text="#1e293b",         # title / slate-900 — wordmark
    text_muted="#475569",   # heading / slate-600 — tagline
)

PAL_DARK = dict(
    bg="#0b1220",            # very dark navy — sits well on GitHub dark theme
    baseline="#334155",      # slate-700
    noise="#94a3b8",         # slate-400 — readable on dark
    noise_soft="#64748b",    # slate-500
    spike="#2dd4bf",         # teal-400 — brighter for dark mode contrast
    spike_dark="#5eead4",    # teal-300 — spike outline (lighter than fill on dark)
    text="#f1f5f9",          # slate-100 — wordmark
    text_muted="#94a3b8",    # slate-400 — tagline
)

FONT_STACK = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif"

# ---- geometry ---------------------------------------------------------------

# Square mark viewBox is 256x256. Baseline sits at y=180.
MARK_VB = (0, 0, 256, 256)
BASELINE_Y = 180
BASELINE_X0, BASELINE_X1 = 22, 234

# Noisy input waveform: smooth-quadratic path through hand-placed waypoints
# between x=24 and x=132. Heights drawn from a fixed "random" pattern so the
# logo is reproducible byte-for-byte.
NOISE_PATH = (
    "M 24,180 "
    "Q 30,180 34,170 "
    "Q 39,156 45,168 "
    "Q 51,180 56,172 "
    "Q 62,160 66,146 "
    "Q 72,160 76,168 "
    "Q 82,180 86,164 "
    "Q 90,148 96,134 "
    "Q 102,150 106,150 "
    "Q 110,150 114,158 "
    "Q 120,172 126,180 "
    "L 132,180"
)

# Spike: a sharp filled triangle. Base 168..192, apex (180, 58).
SPIKE_PATH = "M 168,180 L 180,58 L 192,180 Z"

# Small caret beneath the baseline at the spike's anchor, x=180.
ANCHOR_PATH = "M 176,184 L 180,194 L 184,184 Z"

# Faint genome tick marks on the baseline to add scale + suggest "positions".
# Three short ticks under the noise region, one big tick under the spike anchor.
NOISE_TICKS_X = (44, 76, 108)
SPIKE_TICK_X = 180


def mark_svg(pal: dict, *, include_bg: bool = True) -> str:
    bg_rect = (
        f'<rect width="256" height="256" fill="{pal["bg"]}" />'
        if include_bg
        else ""
    )
    noise_ticks = "".join(
        f'<line x1="{x}" y1="{BASELINE_Y - 2}" x2="{x}" y2="{BASELINE_Y + 2}" '
        f'stroke="{pal["baseline"]}" stroke-width="1.4" stroke-linecap="round" />'
        for x in NOISE_TICKS_X
    )
    spike_tick = (
        f'<line x1="{SPIKE_TICK_X}" y1="{BASELINE_Y - 3}" '
        f'x2="{SPIKE_TICK_X}" y2="{BASELINE_Y + 4}" '
        f'stroke="{pal["spike_dark"]}" stroke-width="1.8" stroke-linecap="round" />'
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 256 256"
     role="img" aria-label="RECTIFY">
  {bg_rect}
  <!-- baseline -->
  <line x1="{BASELINE_X0}" y1="{BASELINE_Y}" x2="{BASELINE_X1}" y2="{BASELINE_Y}"
        stroke="{pal["baseline"]}" stroke-width="1.6" stroke-linecap="round" />
  {noise_ticks}
  {spike_tick}
  <!-- echo of the noisy trace, lighter, to suggest stacked uncorrected reads -->
  <path d="{NOISE_PATH}" fill="none"
        stroke="{pal["noise_soft"]}" stroke-width="2"
        stroke-linejoin="round" stroke-linecap="round"
        transform="translate(0,6) scale(1,0.78)"
        transform-origin="78 180"
        opacity="0.55" />
  <!-- primary noisy trace -->
  <path d="{NOISE_PATH}" fill="none"
        stroke="{pal["noise"]}" stroke-width="2.6"
        stroke-linejoin="round" stroke-linecap="round" />
  <!-- rectified spike -->
  <path d="{SPIKE_PATH}"
        fill="{pal["spike"]}"
        stroke="{pal["spike_dark"]}" stroke-width="2"
        stroke-linejoin="round" />
  <!-- anchor mark just below baseline at the CPA position -->
  <path d="{ANCHOR_PATH}" fill="{pal["spike_dark"]}" />
</svg>
'''


def banner_svg(pal: dict) -> str:
    # 840x256 lockup: mark on the left (scaled into a 256x256 region),
    # wordmark on the right with optional muted subtitle.
    mark_inline = mark_svg(pal, include_bg=False)
    # strip the outer <svg> wrapper and re-embed as a <g>
    inner = mark_inline.split(">", 1)[1].rsplit("</svg>", 1)[0]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 840 256"
     role="img" aria-label="RECTIFY">
  <rect width="840" height="256" fill="{pal["bg"]}" />
  <g transform="translate(0,0)">
    <svg x="0" y="0" width="256" height="256" viewBox="0 0 256 256">
      {inner}
    </svg>
  </g>
  <g font-family="{FONT_STACK}">
    <text x="284" y="158"
          font-size="92" font-weight="700"
          letter-spacing="4"
          fill="{pal["text"]}">RECTIFY</text>
    <text x="288" y="192"
          font-size="18" font-weight="500"
          letter-spacing="2.4"
          fill="{pal["text_muted"]}">NANOPORE  RNA  ALIGNMENT  CORRECTION</text>
  </g>
</svg>
'''


def write_set() -> list[Path]:
    files: list[Path] = []
    pairs = {
        "rectify_mark_light.svg": mark_svg(PAL_LIGHT),
        "rectify_mark_dark.svg": mark_svg(PAL_DARK),
        "rectify_banner_light.svg": banner_svg(PAL_LIGHT),
        "rectify_banner_dark.svg": banner_svg(PAL_DARK),
    }
    for name, content in pairs.items():
        p = OUT / name
        p.write_text(content)
        files.append(p)
    return files


def rasterize(svg_paths: list[Path]) -> list[Path]:
    """Render PNGs at 4x the viewBox using rsvg-convert."""
    pngs: list[Path] = []
    for svg in svg_paths:
        # Determine target size from filename
        if "mark" in svg.name:
            w, h = 1024, 1024  # 4x of 256
        else:  # banner
            w, h = 3360, 1024  # 4x of 840x256
        png = svg.with_suffix(".png")
        subprocess.run(
            [
                "rsvg-convert",
                "-w", str(w),
                "-h", str(h),
                "-o", str(png),
                str(svg),
            ],
            check=True,
        )
        pngs.append(png)
    return pngs


if __name__ == "__main__":
    svgs = write_set()
    pngs = rasterize(svgs)
    print("Wrote:")
    for p in svgs + pngs:
        print(f"  {p.name}  ({p.stat().st_size} bytes)")
