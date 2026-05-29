"""Generate QR codes for the RECTIFY GitHub page.

Two variants:
  rectify_qr_plain.{svg,png}     — black on white, no overlay (max scan reliability)
  rectify_qr_branded.{svg,png}   — navy on white with the teal spike in a white
                                   circular badge centered on the QR (H-level error
                                   correction recovers the masked modules)
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import segno

URL = "https://github.com/k-roy/RECTIFY"
OUT = Path(__file__).parent

# Palette matches build_logo.py / dev/figures/STYLE_GUIDE.md
NAVY = "#1e293b"
TEAL = "#0d9488"
TEAL_DARK = "#0f766e"
WHITE = "#ffffff"


def make_qr(color: str) -> tuple[str, int]:
    """Return (svg_text, qr_module_count) for the URL at H-level ECC."""
    qr = segno.make(URL, error="h", micro=False)
    # Render at unit scale so we can post-process / wrap the SVG
    import io
    buf = io.BytesIO()
    qr.save(
        buf,
        kind="svg",
        scale=1,
        border=2,
        dark=color,
        light=WHITE,
        xmldecl=False,
        svgns=True,
        svgclass=None,
        lineclass=None,
        omitsize=False,
    )
    svg_text = buf.getvalue().decode("utf-8")
    # Module count includes the 2-module border on each side
    n = qr.symbol_size(scale=1, border=2)[0]
    return svg_text, n


def write_plain() -> Path:
    svg, _ = make_qr(color="#000000")
    p = OUT / "rectify_qr_plain.svg"
    p.write_text(svg)
    return p


def write_branded() -> Path:
    svg, n = make_qr(color=NAVY)
    # segno writes <svg ... width="N" height="N" viewBox="0 0 N N">...</svg>.
    # Inject a centered white circular badge with the teal spike at the centre.
    cx = cy = n / 2
    badge_r = n * 0.13          # ~13% of QR side; safe under H ECC (~30% recovery)
    spike_h = badge_r * 1.4
    spike_w = badge_r * 0.55
    overlay = f'''
  <circle cx="{cx:.2f}" cy="{cy:.2f}" r="{badge_r:.2f}"
          fill="{WHITE}" stroke="{NAVY}" stroke-width="{badge_r*0.08:.2f}" />
  <path d="M {cx - spike_w/2:.2f},{cy + spike_h*0.42:.2f}
           L {cx:.2f},{cy - spike_h*0.62:.2f}
           L {cx + spike_w/2:.2f},{cy + spike_h*0.42:.2f} Z"
        fill="{TEAL}" stroke="{TEAL_DARK}" stroke-width="{badge_r*0.04:.2f}"
        stroke-linejoin="round" />
'''
    # Splice overlay just before the closing </svg>
    branded = re.sub(r"</svg>\s*$", overlay + "</svg>\n", svg)
    p = OUT / "rectify_qr_branded.svg"
    p.write_text(branded)
    return p


def rasterize(svgs: list[Path], px: int = 1024) -> list[Path]:
    pngs: list[Path] = []
    for svg in svgs:
        png = svg.with_suffix(".png")
        subprocess.run(
            ["rsvg-convert", "-w", str(px), "-h", str(px), "-o", str(png), str(svg)],
            check=True,
        )
        pngs.append(png)
    return pngs


if __name__ == "__main__":
    svgs = [write_plain(), write_branded()]
    pngs = rasterize(svgs)
    for p in svgs + pngs:
        print(f"  {p.name}  ({p.stat().st_size} bytes)")
    print(f"\nEncoded URL: {URL}")
