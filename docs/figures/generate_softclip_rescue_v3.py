#!/usr/bin/env python3
"""
HP-aware soft-clip rescue figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, three-row layout (genome / read / rescued), mirrors walkback_v3
  - Show why the aligner truncates inside a homopolymer (nanopore HP undercall)
    and how RECTIFY recovers the soft-clipped downstream bases by extending
    the alignment through the missing HP positions to the true CPA
  - No BEFORE / AFTER framing — use factual headers ("Read (raw)",
    "Read (rescued)") that describe what each row IS
  - Inter font, restrained palette, no banners, no marketing copy
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 310
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
    green   = "#059669",
    teal    = "#0d9488",
    teal_l  = "#ccfbf1",
    orange  = "#d97706",
    red     = "#dc2626",
    red_l   = "#fee2e2",
    gap_f   = "#f1f5f9",   # soft slate for deletion / HP-gap fill
    gap_t   = "#94a3b8",
    A_f = "#dcfce7",  A_t = "#166534",
    T_f = "#ffedd5",  T_t = "#9a3412",
    G_f = "#dbeafe",  G_t = "#1e40af",
    C_f = "#f3e8ff",  C_t = "#6b21a8",
)
BASE_FILL = {"A": PAL["A_f"], "T": PAL["T_f"], "G": PAL["G_f"], "C": PAL["C_f"]}
BASE_TEXT = {"A": PAL["A_t"], "T": PAL["T_t"], "G": PAL["G_t"], "C": PAL["C_t"]}

BW, BH, BR = 28, 26, 3
GAP = 2


def svg_open(h):
    return (
        f'<?xml version="1.0" encoding="utf-8" ?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" height="{h}" '
        f'viewBox="0 0 {FIG_W} {h}">\n'
        f'<defs>'
        f'<marker id="arr_t" markerWidth="7" markerHeight="7" refX="7" refY="3.5" '
        f'orient="auto" markerUnits="userSpaceOnUse">'
        f'<path d="M0,0 L0,7 L7,3.5 z" fill="{PAL["teal"]}"/></marker>'
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


def base_box(x, y, b, fill=None, border=None, sw="0.7",
             dashed=False, text_color=None):
    """One nucleotide box at (x, y). Overrides for fill/border/styling."""
    f = fill if fill is not None else BASE_FILL.get(b, PAL["gap_f"])
    bdr = border if border is not None else PAL["border"]
    tc = text_color if text_color is not None else BASE_TEXT.get(b, PAL["gap_t"])
    dash_attr = ' stroke-dasharray="3,2"' if dashed else ''
    return (
        f'<rect fill="{f}" height="{BH}" rx="{BR}" width="{BW}" '
        f'x="{x}" y="{y}" stroke="{bdr}" stroke-width="{sw}"{dash_attr}/>'
        f'<text fill="{tc}" font-size="12" font-weight="600" '
        f'text-anchor="middle" x="{x + BW/2}" y="{y + BH/2 + 4}">{b}</text>'
    )


def gap_box(x, y, dashed=True):
    """Empty placeholder box — represents an HP-undercall gap in the read."""
    dash_attr = ' stroke-dasharray="3,2"' if dashed else ''
    return (
        f'<rect fill="{PAL["gap_f"]}" height="{BH}" rx="{BR}" width="{BW}" '
        f'x="{x}" y="{y}" stroke="{PAL["gap_t"]}" stroke-width="0.8"{dash_attr}/>'
        f'<text fill="{PAL["gap_t"]}" font-size="14" font-weight="400" '
        f'text-anchor="middle" x="{x + BW/2}" y="{y + BH/2 + 5}">·</text>'
    )


def brace_above(x1, x2, y, color, label, label_dy=-9):
    mid = (x1 + x2) / 2
    return (
        f'<path d="M{x1},{y} L{x1},{y-5} L{x2},{y-5} L{x2},{y}" '
        f'fill="none" stroke="{color}" stroke-width="1" opacity="0.7"/>\n'
        + text(mid, y + label_dy, label, size=10, color=color,
               weight="600", anchor="middle")
    )


def brace_below(x1, x2, y, color, label, label_dy=18):
    mid = (x1 + x2) / 2
    return (
        f'<path d="M{x1},{y} L{x1},{y+5} L{x2},{y+5} L{x2},{y}" '
        f'fill="none" stroke="{color}" stroke-width="1" opacity="0.7"/>\n'
        + text(mid, y + label_dy, label, size=10, color=color,
               weight="600", anchor="middle")
    )


def build():
    out = [svg_open(FIG_H)]

    # Title + subtitle
    out.append(text(FIG_W/2, 26, "HP-aware soft-clip rescue (3′ end)",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "extend the alignment downstream through under-called homopolymer bases "
                    "to recover the soft-clipped CPA region",
                    size=11, color=PAL["muted"], anchor="middle"))

    # Sequences (13 positions). Genome has a 9-T homopolymer then G-A-C-T.
    # The post-HP region ends on a non-A so the "true 3′ end" placement is
    # consistent with the walkback gate (which would otherwise eat trailing A's).
    # Read has only 6 T's basecalled (HP undercall) then G-A-C-T soft-clipped.
    # Rescued: alignment extended through 3 HP-gap positions, G-A-C-T aligned.
    n = 13
    content_w = n * (BW + GAP) - GAP
    x0 = (FIG_W - content_w) // 2 + 30          # +30 to balance left row-labels
    x_pos = lambda i: x0 + i * (BW + GAP)

    y_genome  = 80
    y_read    = 132
    y_rescued = 222

    # Row labels
    label_x = x0 - 14
    out.append(text(label_x, y_genome  + BH/2 + 4, "Reference",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(text(label_x, y_read    + BH/2 + 4, "Read (raw)",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(text(label_x, y_rescued + BH/2 + 4, "Read (rescued)",
                    size=11, color=PAL["label"], anchor="end"))

    # ── Reference row: T×9 then G A C T ─────────────────────────────────────
    genome_seq = "TTTTTTTTT" + "GACT"
    for i, b in enumerate(genome_seq):
        out.append(base_box(x_pos(i), y_genome, b))

    # Brace above the 9-T HP marking
    out.append(brace_above(x_pos(0), x_pos(8) + BW, y_genome - 4,
                           PAL["orange"], "9-T homopolymer"))

    # ── Read (raw) row: 6 aligned T's, then 3 phantom (no read bases at the
    # under-called positions), then 4 soft-clipped C G A A drawn UNDER their
    # matching genomic positions ────────────────────────────────────────────
    # Aligned T's at positions 0..5
    for i in range(6):
        out.append(base_box(x_pos(i), y_read, "T"))
    # Soft-clipped bases G A C T drawn at genomic positions 9..12 with dashed
    # red border + faint red fill (style guide: red = error/artifact)
    softclip_seq = "GACT"
    for j, b in enumerate(softclip_seq):
        i = 9 + j
        out.append(base_box(
            x_pos(i), y_read, b,
            fill=PAL["red_l"], border=PAL["red"], sw="1.0", dashed=True,
        ))

    # Apparent 3' end of the raw alignment (right edge of the 6th T)
    raw_end_x = x_pos(5) + BW + 1
    out.append(f'<line stroke="{PAL["red"]}" stroke-width="1.8" '
               f'stroke-dasharray="3,2" '
               f'x1="{raw_end_x}" x2="{raw_end_x}" '
               f'y1="{y_read - 4}" y2="{y_read + BH + 4}"/>')
    out.append(text(raw_end_x, y_read - 7, "apparent 3′ end",
                    size=10, color=PAL["red"], weight="600", anchor="middle"))

    # Brace below the soft-clipped tetramer
    sc_x1 = x_pos(9)
    sc_x2 = x_pos(12) + BW
    out.append(brace_below(sc_x1, sc_x2, y_read + BH + 4,
                           PAL["red"], "soft-clipped (matches genome 3 nt downstream)"))

    # ── Read (rescued) row: 6 aligned T's + 3 HP-gap boxes + 4 aligned bases
    for i in range(6):
        out.append(base_box(x_pos(i), y_rescued, "T"))
    # HP-gap fill (deletion ops bridging the under-called positions 6, 7, 8)
    for i in range(6, 9):
        out.append(gap_box(x_pos(i), y_rescued))
    # Newly aligned C G A A — highlight as rescued (teal accent)
    for j, b in enumerate(softclip_seq):
        i = 9 + j
        out.append(base_box(
            x_pos(i), y_rescued, b,
            fill=PAL["teal_l"], border=PAL["teal"], sw="1.2",
        ))

    # Brace above the gap region — explains the rescue mechanism
    gap_x1 = x_pos(6)
    gap_x2 = x_pos(8) + BW
    out.append(brace_above(gap_x1, gap_x2, y_rescued - 4,
                           PAL["heading"],
                           "3 HP-undercall gap (deletion ops)",
                           label_dy=-9))

    # True 3' end marker on rescued row (right edge of the last A box)
    true_end_x = x_pos(12) + BW + 1
    out.append(f'<line stroke="{PAL["teal"]}" stroke-width="2" '
               f'x1="{true_end_x}" x2="{true_end_x}" '
               f'y1="{y_rescued - 4}" y2="{y_rescued + BH + 4}"/>')
    out.append(text(true_end_x, y_rescued + BH + 18, "true 3′ end",
                    size=10, color=PAL["teal"], weight="700", anchor="middle"))

    # Downstream rescue arrow connecting raw 3' end → true 3' end (teal)
    arrow_y = y_rescued + BH + 38
    out.append(f'<line stroke="{PAL["teal"]}" stroke-width="1.5" '
               f'x1="{raw_end_x}" x2="{true_end_x - 3}" '
               f'y1="{arrow_y}" y2="{arrow_y}" marker-end="url(#arr_t)"/>')
    out.append(text((raw_end_x + true_end_x) / 2, arrow_y - 6,
                    "extend alignment downstream",
                    size=10, color=PAL["teal"], weight="600", anchor="middle"))

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
    render(build(), "softclip_rescue")
