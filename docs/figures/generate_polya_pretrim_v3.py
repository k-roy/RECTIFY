#!/usr/bin/env python3
"""
DRS poly(A) pre-trim figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, single panel, two rows (raw vs trimmed)
  - No three-pass detection internals; that's an implementation detail
  - Show: raw read = body + poly(A) + adapter stub; trimmed = body only,
    with poly(A)/adapter saved to soft-clip metadata for later restoration
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 210
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
    teal    = "#0d9488",   teal_l   = "#ccfbf1",
    green   = "#059669",   green_l  = "#d1fae5",
    red     = "#dc2626",   red_l    = "#fee2e2",
    orange  = "#d97706",   orange_l = "#fed7aa",
    indigo  = "#4f46e5",   indigo_l = "#e0e7ff",
)


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

def block(x, y, w, h, label, fill, stroke, label_color, dashed=False):
    da = ' stroke-dasharray="4,3"' if dashed else ""
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="1.2"{da}/>\n'
        + text(x + w/2, y + h/2 + 4, label, size=12, color=label_color,
               weight="600", anchor="middle")
    )


def build():
    out = [svg_open(FIG_H)]

    out.append(text(FIG_W/2, 26, "DRS poly(A) pre-trim",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "strip the basecalled poly(A) tail + adapter stub before alignment",
                    size=11, color=PAL["muted"], anchor="middle"))

    # Layout — 3 segments: mRNA body | poly(A) tail | adapter stub
    label_x = 130
    body_x, body_w   = 145, 320
    polya_x, polya_w = body_x + body_w,        160     # 465 → 625
    adapter_x, adapter_w = polya_x + polya_w,  60      # 625 → 685
    row_h = 30

    # Row 1: raw read
    y1 = 78
    out.append(text(label_x, y1 + row_h/2 + 4, "Raw read",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(block(body_x, y1, body_w, row_h, "mRNA body",
                     PAL["teal_l"], PAL["teal"], PAL["teal"]))
    out.append(block(polya_x, y1, polya_w, row_h, "poly(A) tail",
                     PAL["green_l"], PAL["green"], PAL["green"]))
    out.append(block(adapter_x, y1, adapter_w, row_h, "adapter",
                     PAL["orange_l"], PAL["orange"], PAL["orange"]))

    # 3' end label
    out.append(text(adapter_x + adapter_w + 8, y1 + row_h/2 + 4,
                    "3′", size=11, color=PAL["muted"]))

    # Arrow down with annotation
    arrow_x = body_x + body_w/2
    out.append(f'<line x1="{arrow_x}" y1="{y1 + row_h + 6}" '
               f'x2="{arrow_x}" y2="{y1 + row_h + 36}" '
               f'stroke="{PAL["label"]}" stroke-width="1.4" marker-end="url(#arr)"/>')
    out.append(text(arrow_x + 12, y1 + row_h + 24, "rectify trim-polya",
                    size=11, color=PAL["heading"], weight="600"))

    # Row 2: trimmed read
    y2 = 154
    out.append(text(label_x, y2 + row_h/2 + 4, "Trimmed",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(block(body_x, y2, body_w, row_h, "mRNA body",
                     PAL["teal_l"], PAL["teal"], PAL["teal"]))
    # Cached metadata block — dashed to show it's no longer in the read
    out.append(block(polya_x, y2, polya_w + adapter_w, row_h,
                     "cached as soft-clip metadata",
                     "#ffffff", PAL["muted"], PAL["muted"], dashed=True))

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
    render(build(), "polya_pretrim")
