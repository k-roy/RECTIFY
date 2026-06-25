#!/usr/bin/env python3
"""
Splice classification figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, single panel
  - Show annotated intron once + four read examples (one per class)
  - Use vertical donor/acceptor markers on the genome to show where annotated
    splice sites fall; each read's N-op endpoints either align or don't
  - Class label on the left in a distinct color per class
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
    indigo  = "#4f46e5",
    indigo_l= "#e0e7ff",
    teal    = "#0d9488",
    teal_l  = "#ccfbf1",
    green   = "#059669",
    blue    = "#2563eb",
    orange  = "#d97706",
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

def exon(x, y, w, h, label="", color=None, label_color="#ffffff"):
    c = color or PAL["indigo"]
    parts = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" '
             f'fill="{c}"/>']
    if label:
        parts.append(text(x + w/2, y + h/2 + 4, label, size=12,
                          color=label_color, weight="600", anchor="middle"))
    return "\n".join(parts)

def aligned_block(x, y, w, h, color=None):
    c = color or PAL["teal"]
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" '
            f'fill="{PAL["teal_l"]}" stroke="{c}" stroke-width="1.2"/>')

def n_op(x1, x2, y, color=None):
    """Thin dashed line representing an N-cigar op (skipped region)."""
    c = color or PAL["teal"]
    return (f'<line x1="{x1}" x2="{x2}" y1="{y}" y2="{y}" '
            f'stroke="{c}" stroke-width="1.5" stroke-dasharray="4,3"/>')

def donor_acceptor_marker(x, y1, y2, color=None):
    """Vertical dashed line marking an annotated splice site."""
    c = color or PAL["muted"]
    return (f'<line x1="{x}" x2="{x}" y1="{y1}" y2="{y2}" '
            f'stroke="{c}" stroke-width="1" stroke-dasharray="2,2" opacity="0.7"/>')

def hdivider(y, x1=60, x2=720):
    return f'<line x1="{x1}" x2="{x2}" y1="{y}" y2="{y}" stroke="{PAL["divider"]}" stroke-width="1"/>'


def build():
    out = [svg_open(FIG_H)]

    # Title
    out.append(text(FIG_W/2, 26, "Splice classification",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "classify each N-cigar op vs the reference annotation · \"novel\" = not in the annotation",
                    size=10, color=PAL["muted"], anchor="middle"))

    # Annotated reference layout — two annotated splice sites at fixed x
    # Exon 1 ends at donor_x; Exon 2 starts at acceptor_x
    x_label = 145
    x_e1_a, x_e1_b = 160, 305      # Exon 1
    donor_x        = 305            # annotated donor
    acceptor_x     = 505            # annotated acceptor
    x_e2_a, x_e2_b = 505, 660      # Exon 2
    row_h = 22

    # ── Reference row ─────────────────────────────────────────────────────────
    y_ref = 80
    out.append(text(x_label, y_ref + row_h/2 + 4, "annotated",
                    size=11, color=PAL["label"], anchor="end", weight="600"))
    out.append(exon(x_e1_a, y_ref, x_e1_b - x_e1_a, row_h, "Exon 1"))
    out.append(exon(x_e2_a, y_ref, x_e2_b - x_e2_a, row_h, "Exon 2"))
    # Faint intron link
    out.append(f'<line x1="{donor_x}" x2="{acceptor_x}" y1="{y_ref + row_h/2}" '
               f'y2="{y_ref + row_h/2}" stroke="{PAL["muted"]}" '
               f'stroke-width="0.8" stroke-dasharray="3,3"/>')
    out.append(text((donor_x + acceptor_x)/2, y_ref - 4,
                    "intron", size=10, color=PAL["muted"], anchor="middle"))

    # ── Vertical dashed lines extending down across all read rows ────────────
    # so the reader can eyeball whether N-op endpoints align with annotation
    y_rows_top = 130
    y_rows_bot = 305
    out.append(donor_acceptor_marker(donor_x, y_rows_top, y_rows_bot))
    out.append(donor_acceptor_marker(acceptor_x, y_rows_top, y_rows_bot))

    # Helper to draw a labeled read row
    def read_row(y, label, label_color, n_ops, aligned_spans):
        parts = [text(x_label, y + row_h/2 + 4, label, size=11,
                      color=label_color, anchor="end", weight="700")]
        for (xa, xb) in aligned_spans:
            parts.append(aligned_block(xa, y, xb - xa, row_h))
        for (xa, xb) in n_ops:
            parts.append(n_op(xa, xb, y + row_h/2))
        return "\n".join(parts)

    # Class 1: unspliced (no N-op, one long aligned block spanning intron)
    y1 = 140
    out.append(read_row(y1, "unspliced", PAL["heading"],
                        n_ops=[],
                        aligned_spans=[(160, 660)]))

    # Class 2: annotated (N-op endpoints exactly at donor and acceptor)
    y2 = 180
    out.append(read_row(y2, "annotated", PAL["green"],
                        n_ops=[(donor_x, acceptor_x)],
                        aligned_spans=[(160, donor_x), (acceptor_x, 660)]))

    # Class 3: alternative (one end aligns to annotation, other does not)
    # donor matches annotation, acceptor is 30 px to the right
    y3 = 220
    alt_acceptor = acceptor_x + 40
    out.append(read_row(y3, "one-side novel", PAL["orange"],
                        n_ops=[(donor_x, alt_acceptor)],
                        aligned_spans=[(160, donor_x), (alt_acceptor, 660)]))

    # Class 4: novel (neither end matches)
    y4 = 260
    nov_donor = donor_x + 35
    nov_acceptor = acceptor_x + 35
    out.append(read_row(y4, "both-side novel", PAL["red"],
                        n_ops=[(nov_donor, nov_acceptor)],
                        aligned_spans=[(160, nov_donor), (nov_acceptor, 660)]))

    # Labels above the dashed annotation lines (small, only on top row)
    out.append(text(donor_x, y_rows_top - 4, "donor",
                    size=10, color=PAL["muted"], anchor="middle"))
    out.append(text(acceptor_x, y_rows_top - 4, "acceptor",
                    size=10, color=PAL["muted"], anchor="middle"))

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
    render(build(), "splice_classification")
