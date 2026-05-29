"""Full deliverable set for concept C (anchored strand) + the GitHub | Docs
two-QR poster sheet.

Outputs in this directory:
  rectify_mark_C_light.svg/png
  rectify_mark_C_dark.svg/png
  rectify_banner_C_light.svg/png       (same content as rectify_banner_C_anchored.svg)
  rectify_banner_C_dark.svg/png
  rectify_qr_sheet_light.svg/png       (GitHub + ReadTheDocs side-by-side)
  rectify_qr_sheet_dark.svg/png

The QR sheet uses the anchored-strand mark in the header so the page reads as
"a RECTIFY artefact" rather than two anonymous QRs.
"""

from __future__ import annotations

import io
import re
import subprocess
from pathlib import Path

import segno

OUT = Path(__file__).parent

URL_GITHUB = "https://github.com/k-roy/RECTIFY"
URL_DOCS = "https://rectify-rna.readthedocs.io"

FONT = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif"

# --- palettes ---------------------------------------------------------------

PAL_LIGHT = dict(
    bg="#ffffff",
    text="#1e293b",
    text_muted="#475569",
    strand="#64748b",         # slate-500 — backbone
    A_fill="#dcfce7", A_text="#166534",
    spike="#0d9488",          # teal-600
    spike_dark="#0f766e",     # teal-700
    qr_dark="#1e293b",        # QR modules — navy reads as black at small sizes
    qr_light="#ffffff",
    rule="#e2e8f0",           # subtle dividers
)

PAL_DARK = dict(
    bg="#0b1220",
    text="#f1f5f9",
    text_muted="#94a3b8",
    strand="#94a3b8",         # slate-400 on dark
    A_fill="#134e3a", A_text="#a7f3d0",  # darker green box, light A text
    spike="#2dd4bf",          # teal-400
    spike_dark="#5eead4",     # teal-300
    qr_dark="#f1f5f9",        # invert modules for dark sheet
    qr_light="#0b1220",
    rule="#1f2937",
)


# --- concept C mark (parameterised by palette) ------------------------------

def mark_C(pal: dict, *, include_bg: bool = True) -> str:
    backbone = (
        "M 24,128 "
        "C 36,108 48,148 60,128 "
        "C 72,108 84,148 96,128 "
        "C 108,108 120,148 132,128"
    )
    a_positions = [148, 162, 176, 190, 204, 218]
    a_bgs = "".join(
        f'<rect x="{x-6}" y="106" width="12" height="14" rx="2.5" '
        f'fill="{pal["A_fill"]}" opacity="0.9" />'
        for x in a_positions
    )
    a_labels = "".join(
        f'<text x="{x}" y="116" text-anchor="middle" '
        f'font-family="{FONT}" font-size="13" font-weight="600" '
        f'fill="{pal["A_text"]}">A</text>'
        for x in a_positions
    )
    pin = f'''
  <line x1="132" y1="92" x2="132" y2="142" stroke="{pal["spike"]}" stroke-width="3" stroke-linecap="round" />
  <circle cx="132" cy="88" r="5.5" fill="{pal["spike"]}" />
  <path d="M 127,148 L 132,160 L 137,148 Z" fill="{pal["spike_dark"]}" />
'''
    bg_rect = (
        f'<rect width="256" height="256" fill="{pal["bg"]}" />'
        if include_bg else ""
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 256 256"
     role="img" aria-label="RECTIFY">
  {bg_rect}
  <path d="{backbone}" fill="none"
        stroke="{pal["strand"]}" stroke-width="2.6"
        stroke-linejoin="round" stroke-linecap="round" />
  <line x1="132" y1="128" x2="224" y2="128"
        stroke="{pal["strand"]}" stroke-width="2.6" stroke-linecap="round" />
  {a_bgs}
  {a_labels}
  {pin}
</svg>
'''


def banner_C(pal: dict) -> str:
    inner = mark_C(pal, include_bg=False).split(">", 1)[1].rsplit("</svg>", 1)[0]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 840 256"
     role="img" aria-label="RECTIFY">
  <rect width="840" height="256" fill="{pal["bg"]}" />
  <svg x="0" y="0" width="256" height="256" viewBox="0 0 256 256">{inner}</svg>
  <g font-family="{FONT}">
    <text x="284" y="158" font-size="92" font-weight="700" letter-spacing="4"
          fill="{pal["text"]}">RECTIFY</text>
    <text x="288" y="192" font-size="18" font-weight="500" letter-spacing="2.4"
          fill="{pal["text_muted"]}">NANOPORE  RNA  ALIGNMENT  CORRECTION</text>
  </g>
</svg>
'''


# --- QR generation -----------------------------------------------------------

def qr_svg_inner(url: str, pal: dict) -> tuple[str, int]:
    """Return the QR <svg>...</svg> content and its module-grid size."""
    qr = segno.make(url, error="h", micro=False)
    buf = io.BytesIO()
    qr.save(
        buf, kind="svg", scale=1, border=2,
        dark=pal["qr_dark"], light=pal["qr_light"],
        xmldecl=False, svgns=True,
        svgclass=None, lineclass=None, omitsize=False,
    )
    svg_text = buf.getvalue().decode("utf-8")
    n = qr.symbol_size(scale=1, border=2)[0]
    return svg_text, n


def branded_qr(url: str, pal: dict) -> tuple[str, int]:
    """QR with the teal-spike circular badge at the centre."""
    svg, n = qr_svg_inner(url, pal)
    cx = cy = n / 2
    badge_r = n * 0.13
    spike_h = badge_r * 1.4
    spike_w = badge_r * 0.55
    overlay = f'''
  <circle cx="{cx:.2f}" cy="{cy:.2f}" r="{badge_r:.2f}"
          fill="{pal["qr_light"]}" stroke="{pal["qr_dark"]}" stroke-width="{badge_r*0.08:.2f}" />
  <path d="M {cx - spike_w/2:.2f},{cy + spike_h*0.42:.2f}
           L {cx:.2f},{cy - spike_h*0.62:.2f}
           L {cx + spike_w/2:.2f},{cy + spike_h*0.42:.2f} Z"
        fill="{pal["spike"]}" stroke="{pal["spike_dark"]}" stroke-width="{badge_r*0.04:.2f}"
        stroke-linejoin="round" />
'''
    return re.sub(r"</svg>\s*$", overlay + "</svg>\n", svg), n


# --- QR sheet ----------------------------------------------------------------

def qr_sheet(pal: dict) -> str:
    # Canvas: 1400 x 900 — landscape strip for poster footer / handout
    # Header: anchored-strand mark + RECTIFY wordmark, centered, y=40..200
    # Divider rule at y=232
    # Section labels at y=270 ("CODE" / "DOCS")
    # QR codes 460x460 at (270, 310) and (970, 310), bottom y=770
    # URL captions centered under each QR at y=820

    # The two QR images are themselves SVG snippets; embed via <svg> nested element
    # with their own viewBox = "0 0 n n".
    gh_svg, gh_n = branded_qr(URL_GITHUB, pal)
    docs_svg, docs_n = branded_qr(URL_DOCS, pal)
    # Extract inner contents (between opening <svg ...> and </svg>)
    def inner_of(svg: str) -> str:
        return re.sub(r"^<svg[^>]*>", "", svg.strip()).rsplit("</svg>", 1)[0]
    gh_inner = inner_of(gh_svg)
    docs_inner = inner_of(docs_svg)

    # Header mark (anchored-strand) without background, sized 220x220 at (130, 30)
    header_mark_inner = mark_C(pal, include_bg=False).split(">", 1)[1].rsplit("</svg>", 1)[0]

    qr_size = 460
    qr_l_x, qr_r_x = 220, 920    # left edge of each QR
    qr_y = 320

    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1400 900"
     role="img" aria-label="RECTIFY — GitHub and Docs">
  <rect width="1400" height="900" fill="{pal["bg"]}" />

  <!-- header: mark + wordmark, horizontally centered -->
  <g transform="translate(380,18)">
    <svg x="0" y="0" width="200" height="200" viewBox="0 0 256 256">{header_mark_inner}</svg>
  </g>
  <g font-family="{FONT}">
    <text x="600" y="128" font-size="80" font-weight="700" letter-spacing="4"
          fill="{pal["text"]}">RECTIFY</text>
    <text x="604" y="160" font-size="18" font-weight="500" letter-spacing="2.4"
          fill="{pal["text_muted"]}">NANOPORE  RNA  ALIGNMENT  CORRECTION</text>
  </g>

  <!-- divider -->
  <line x1="120" y1="240" x2="1280" y2="240" stroke="{pal["rule"]}" stroke-width="1.6" />

  <!-- left section: CODE / GitHub -->
  <g font-family="{FONT}" text-anchor="middle" fill="{pal["text_muted"]}">
    <text x="{qr_l_x + qr_size/2}" y="296" font-size="22" font-weight="700" letter-spacing="6">CODE</text>
  </g>
  <svg x="{qr_l_x}" y="{qr_y}" width="{qr_size}" height="{qr_size}"
       viewBox="0 0 {gh_n} {gh_n}">{gh_inner}</svg>
  <g font-family="{FONT}" text-anchor="middle">
    <text x="{qr_l_x + qr_size/2}" y="{qr_y + qr_size + 50}" font-size="22"
          font-weight="600" fill="{pal["text"]}">github.com/k-roy/RECTIFY</text>
    <text x="{qr_l_x + qr_size/2}" y="{qr_y + qr_size + 80}" font-size="15"
          font-weight="500" fill="{pal["text_muted"]}" letter-spacing="0.8">Source code, issues, releases</text>
  </g>

  <!-- right section: DOCS / ReadTheDocs -->
  <g font-family="{FONT}" text-anchor="middle" fill="{pal["text_muted"]}">
    <text x="{qr_r_x + qr_size/2}" y="296" font-size="22" font-weight="700" letter-spacing="6">DOCS</text>
  </g>
  <svg x="{qr_r_x}" y="{qr_y}" width="{qr_size}" height="{qr_size}"
       viewBox="0 0 {docs_n} {docs_n}">{docs_inner}</svg>
  <g font-family="{FONT}" text-anchor="middle">
    <text x="{qr_r_x + qr_size/2}" y="{qr_y + qr_size + 50}" font-size="22"
          font-weight="600" fill="{pal["text"]}">rectify-rna.readthedocs.io</text>
    <text x="{qr_r_x + qr_size/2}" y="{qr_y + qr_size + 80}" font-size="15"
          font-weight="500" fill="{pal["text_muted"]}" letter-spacing="0.8">User guide, API reference, methods</text>
  </g>
</svg>
'''


# --- build + rasterize ------------------------------------------------------

def build() -> list[Path]:
    files: dict[str, str] = {
        "rectify_mark_C_light.svg":    mark_C(PAL_LIGHT),
        "rectify_mark_C_dark.svg":     mark_C(PAL_DARK),
        "rectify_banner_C_light.svg":  banner_C(PAL_LIGHT),
        "rectify_banner_C_dark.svg":   banner_C(PAL_DARK),
        "rectify_qr_sheet_light.svg":  qr_sheet(PAL_LIGHT),
        "rectify_qr_sheet_dark.svg":   qr_sheet(PAL_DARK),
    }
    paths: list[Path] = []
    for name, content in files.items():
        p = OUT / name
        p.write_text(content)
        paths.append(p)
    return paths


def rasterize(svgs: list[Path]) -> list[Path]:
    out: list[Path] = []
    for svg in svgs:
        name = svg.name
        if "mark" in name:
            w, h = 1024, 1024
        elif "banner" in name:
            w, h = 3360, 1024
        elif "qr_sheet" in name:
            w, h = 2800, 1800   # 2x of 1400x900 — print-friendly without bloat
        else:
            w, h = 1024, 1024
        png = svg.with_suffix(".png")
        subprocess.run(
            ["rsvg-convert", "-w", str(w), "-h", str(h), "-o", str(png), str(svg)],
            check=True,
        )
        out.append(png)
    return out


if __name__ == "__main__":
    svgs = build()
    pngs = rasterize(svgs)
    for p in svgs + pngs:
        print(f"  {p.name}  ({p.stat().st_size} bytes)")
