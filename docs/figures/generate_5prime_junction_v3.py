#!/usr/bin/env python3
"""
5' end junction rescue figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, single panel
  - 3 rows: genome / read aligned (before) / read after rescue
  - No per-base nucleotide grids — use exon/aligned blocks at this scale
  - No "Problem / Solution" framing; descriptive row labels only
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 230
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
    blue_l  = "#dbeafe",
    teal    = "#0d9488",
    teal_l  = "#ccfbf1",
    red     = "#dc2626",
    red_l   = "#fee2e2",
    indigo  = "#4f46e5",
    indigo_l= "#e0e7ff",
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

def exon(x, y, w, h, label, fill=None, stroke=None, label_color=None):
    f = fill or PAL["indigo"]
    s = stroke or fill or PAL["indigo"]
    lc = label_color or "#ffffff"
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" '
        f'fill="{f}" stroke="{s}" stroke-width="1"/>\n'
        + text(x + w/2, y + h/2 + 4, label, size=12, color=lc,
               weight="600", anchor="middle")
    )

def softclip(x, y, w, h, label):
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" '
        f'fill="{PAL["red_l"]}" stroke="{PAL["red"]}" stroke-width="1" '
        f'stroke-dasharray="4,3"/>\n'
        + text(x + w/2, y + h/2 + 4, label, size=12, color=PAL["red"],
               weight="600", anchor="middle")
    )

def aligned_block(x, y, w, h, label, color=None):
    """Solid teal block for rescued/aligned regions."""
    c = color or PAL["teal"]
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" '
        f'fill="{PAL["teal_l"]}" stroke="{c}" stroke-width="1.2"/>\n'
        + text(x + w/2, y + h/2 + 4, label, size=12, color=c,
               weight="600", anchor="middle")
    )

def intron_line(x1, x2, y, label="intron"):
    """Dashed line connecting two exons."""
    return (
        f'<line x1="{x1}" x2="{x2}" y1="{y}" y2="{y}" '
        f'stroke="{PAL["muted"]}" stroke-width="1.2" stroke-dasharray="5,4"/>\n'
        + text((x1+x2)/2, y - 6, label, size=10, color=PAL["muted"],
               anchor="middle")
    )

def n_op_line(x1, x2, y, label="N (splice)"):
    """Dashed line for an N cigar op in a corrected read."""
    return (
        f'<line x1="{x1}" x2="{x2}" y1="{y}" y2="{y}" '
        f'stroke="{PAL["teal"]}" stroke-width="1.5" stroke-dasharray="4,3"/>\n'
        + text((x1+x2)/2, y - 6, label, size=10, color=PAL["teal"],
               anchor="middle")
    )


def build():
    out = [svg_open(FIG_H)]

    out.append(text(FIG_W/2, 26, "5′ end junction rescue",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "extend soft-clipped bases across the splice junction to the canonical donor",
                    size=11, color=PAL["muted"], anchor="middle"))

    # Layout columns: Exon 1 | intron | Exon 2
    x_e1 = 130
    w_e1 = 160
    x_int_a = x_e1 + w_e1                # 290
    x_int_b = 470
    w_int = x_int_b - x_int_a            # 180
    x_e2 = x_int_b                       # 470
    w_e2 = 160
    x_end = x_e2 + w_e2                  # 630

    row_h = 28
    label_x = 120

    # Row 1: Genome
    y1 = 78
    out.append(text(label_x, y1 + row_h/2 + 4, "Genome",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(exon(x_e1, y1, w_e1, row_h, "Exon 1"))
    out.append(intron_line(x_int_a, x_int_b, y1 + row_h/2, "intron (GT…AG)"))
    out.append(exon(x_e2, y1, w_e2, row_h, "Exon 2"))

    # Row 2: Read aligned (before)
    y2 = 130
    out.append(text(label_x, y2 + row_h/2 + 4, "Read (before)",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(softclip(x_e1, y2, w_e1, row_h, "soft-clipped"))
    # No alignment across intron — note the gap
    out.append(aligned_block(x_e2, y2, w_e2, row_h, "aligned"))
    # Annotation above the soft-clip
    out.append(text(x_e1 + w_e1/2, y2 - 6,
                    "bases actually match Exon 1",
                    size=11, color=PAL["red"], anchor="middle", weight="600"))

    # Row 3: Read after rescue
    y3 = 182
    out.append(text(label_x, y3 + row_h/2 + 4, "Read (rescued)",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(aligned_block(x_e1, y3, w_e1, row_h, "aligned"))
    out.append(n_op_line(x_int_a, x_int_b, y3 + row_h/2, "N (splice)"))
    out.append(aligned_block(x_e2, y3, w_e2, row_h, "aligned"))

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
    render(build(), "5prime_junction_rescue")
