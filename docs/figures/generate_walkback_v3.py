#!/usr/bin/env python3
"""
Read-vs-reference walkback figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, single panel (no old-vs-new comparison)
  - Show genome + read + walkback arrow + stop marker; nothing else
  - Three terminal-gate cases shown as compact braces over the read,
    not in a sidebar legend
  - No "WRONG/RIGHT" framing, no protocol-wiring box
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 260
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
    green   = "#059669",   green_l  = "#d1fae5",
    teal    = "#0d9488",
    teal_l  = "#ccfbf1",
    orange  = "#d97706",   orange_l = "#fed7aa",
    red     = "#dc2626",
    red_l   = "#fee2e2",
    A_f = "#dcfce7",  A_t = "#166534",
    T_f = "#ffedd5",  T_t = "#9a3412",
    G_f = "#dbeafe",  G_t = "#1e40af",
    C_f = "#f3e8ff",  C_t = "#6b21a8",
)
BASE_FILL = {"A": PAL["A_f"], "T": PAL["T_f"], "G": PAL["G_f"], "C": PAL["C_f"]}
BASE_TEXT = {"A": PAL["A_t"], "T": PAL["T_t"], "G": PAL["G_t"], "C": PAL["C_t"]}

BW, BH, BR = 28, 26, 3  # box dimensions
GAP = 2                 # px between boxes


def svg_open(h):
    return (
        f'<?xml version="1.0" encoding="utf-8" ?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" height="{h}" '
        f'viewBox="0 0 {FIG_W} {h}">\n'
        f'<defs>'
        f'<marker id="arr" markerWidth="7" markerHeight="7" refX="7" refY="3.5" '
        f'orient="auto" markerUnits="userSpaceOnUse">'
        f'<path d="M0,0 L0,7 L7,3.5 z" fill="{PAL["blue"]}"/></marker>'
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

def base_box(x, y, b, highlight=None):
    """One nucleotide box. `highlight` = ('fill', 'border') override."""
    fill = highlight[0] if highlight else BASE_FILL[b]
    bdr  = highlight[1] if highlight else PAL["border"]
    sw   = "1.2" if highlight else "0.7"
    return (
        f'<rect fill="{fill}" height="{BH}" rx="{BR}" width="{BW}" '
        f'x="{x}" y="{y}" stroke="{bdr}" stroke-width="{sw}"/>'
        f'<text fill="{BASE_TEXT[b]}" font-size="12" font-weight="600" '
        f'text-anchor="middle" x="{x + BW/2}" y="{y + BH/2 + 4}">{b}</text>'
    )

def row(x0, y, seq, highlights=None):
    """Draw a row of nucleotide boxes starting at x0."""
    highlights = highlights or {}
    out = []
    for i, b in enumerate(seq):
        x = x0 + i * (BW + GAP)
        out.append(base_box(x, y, b, highlight=highlights.get(i)))
    return "\n".join(out)

def brace_above(x1, x2, y, color, label, label_dy=-9):
    """Brace opening upward, label centered above."""
    mid = (x1 + x2) / 2
    return (
        f'<path d="M{x1},{y} L{x1},{y-5} L{x2},{y-5} L{x2},{y}" '
        f'fill="none" stroke="{color}" stroke-width="1" opacity="0.7"/>\n'
        + text(mid, y + label_dy, label, size=10, color=color, anchor="middle",
               weight="600")
    )

def h_arrow_left(x_from, x_to, y, color, label=None):
    """Leftward arrow from x_from to x_to (x_to < x_from)."""
    parts = [
        f'<line stroke="{color}" stroke-width="1.5" '
        f'x1="{x_from}" x2="{x_to}" y1="{y}" y2="{y}"/>',
        f'<polygon fill="{color}" '
        f'points="{x_to},{y} {x_to+7},{y-3.5} {x_to+7},{y+3.5}"/>',
    ]
    if label:
        parts.append(text((x_from + x_to) / 2, y + 14, label, size=11,
                          color=color, anchor="middle"))
    return "\n".join(parts)

def vert_marker(x, y1, y2, color, label, label_above=True):
    parts = [f'<line stroke="{color}" stroke-width="2" '
             f'x1="{x}" x2="{x}" y1="{y1}" y2="{y2}"/>']
    if label_above:
        parts.append(text(x, y1 - 5, label, size=10, color=color,
                          weight="700", anchor="middle"))
    else:
        parts.append(text(x, y2 + 14, label, size=10, color=color,
                          weight="700", anchor="middle"))
    return "\n".join(parts)


# ─────────────────────────────────────────────────────────────────────────────
# Build
# ─────────────────────────────────────────────────────────────────────────────

def build():
    out = [svg_open(FIG_H)]

    # Title
    out.append(text(FIG_W/2, 26, "Read-vs-reference walkback",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "walk past A=A pairs AND tail seq errors to the true CPA",
                    size=11, color=PAL["muted"], anchor="middle"))

    # Sequences.
    # Positions 0-9 = trimmed read aligned to genome.
    #   Position 2 = T/T match, non-A → true CPA, walkback stops here
    #   Positions 3-8 = A/A → case 2 (A=A, walk in)
    #   Position 9 = T/A → case 3 (terminal seq error after pre-trim, walk in)
    # Positions 10-13 = pA tail that pre-trim removed (non-genomic, faded green)
    # Positions 14-15 = DRS adapter stub "TC" that pre-trim removed (orange)
    genome    = "CCTAAAAAAA"  # 10 genomic bases
    read_aln  = "CCTAAAAAAT"  # 10 trimmed-read bases (aligned portion)
    pa_tail   = "AAAA"        # 4 pA tail bases removed by pre-trim
    adapter   = "TC"          # DRS adapter stub removed by pre-trim

    n_genome = len(genome)           # 10
    n_pa     = len(pa_tail)          # 4
    n_ada    = len(adapter)          # 2
    n_total  = n_genome + n_pa + n_ada  # 16

    # Layout — center on all 16 positions
    content_w = n_total * (BW + GAP) - GAP   # 478
    x0 = (FIG_W - content_w) // 2 + 30       # +30 balances row labels
    y_genome = 82
    y_read   = 134
    x_pos = lambda i: x0 + i * (BW + GAP)

    # Row labels
    out.append(text(x0 - 12, y_genome + BH/2 + 4, "Genome",
                    size=11, color=PAL["label"], anchor="end"))
    out.append(text(x0 - 12, y_read + BH/2 + 4, "Read",
                    size=11, color=PAL["label"], anchor="end"))

    # Genome row (positions 0-9 only)
    out.append(row(x0, y_genome, genome))
    # Genome continuation indicator — dashed line showing genome goes on
    genome_end_x = x_pos(n_genome - 1) + BW
    out.append(f'<line x1="{genome_end_x + 4}" x2="{x_pos(n_total - 1) + BW + 4}" '
               f'y1="{y_genome + BH/2}" y2="{y_genome + BH/2}" '
               f'stroke="{PAL["muted"]}" stroke-width="1" '
               f'stroke-dasharray="3,4" opacity="0.5"/>')

    # Read row — aligned portion (positions 0-9) with walkback highlights
    read_highlights = {
        i: (PAL["teal_l"], PAL["teal"]) for i in range(3, 9)
    }
    read_highlights[9] = (PAL["red_l"], PAL["red"])  # terminal seq error
    out.append(row(x0, y_read, read_aln, highlights=read_highlights))

    # Pre-trim separator (dashed vertical between position 9 and 10)
    sep_x = x_pos(n_genome) - (GAP // 2) - 1
    out.append(f'<line stroke="{PAL["heading"]}" stroke-width="1.4" '
               f'stroke-dasharray="4,3" '
               f'x1="{sep_x}" x2="{sep_x}" '
               f'y1="{y_genome - 4}" y2="{y_read + BH + 4}"/>')

    # pA tail boxes (positions 10-13) — faded green, dashed border
    for i, b in enumerate(pa_tail):
        xi = x_pos(n_genome + i)
        out.append(f'<rect fill="{PAL["green_l"]}" height="{BH}" rx="{BR}" '
                   f'width="{BW}" x="{xi}" y="{y_read}" '
                   f'stroke="{PAL["green"]}" stroke-width="0.8" '
                   f'stroke-dasharray="3,2" opacity="0.7"/>')
        out.append(f'<text fill="{PAL["green"]}" font-size="12" font-weight="600" '
                   f'text-anchor="middle" x="{xi + BW/2}" '
                   f'y="{y_read + BH/2 + 4}" opacity="0.7">{b}</text>')

    # Adapter boxes (positions 14-15) — orange
    for i, b in enumerate(adapter):
        xi = x_pos(n_genome + n_pa + i)
        out.append(f'<rect fill="{PAL["orange_l"]}" height="{BH}" rx="{BR}" '
                   f'width="{BW}" x="{xi}" y="{y_read}" '
                   f'stroke="{PAL["orange"]}" stroke-width="1.2"/>')
        out.append(f'<text fill="{PAL["orange"]}" font-size="12" font-weight="600" '
                   f'text-anchor="middle" x="{xi + BW/2}" '
                   f'y="{y_read + BH/2 + 4}">{b}</text>')

    # 3' label at far right
    out.append(text(x_pos(n_total - 1) + BW + 8, y_read + BH/2 + 4,
                    "3′", size=11, color=PAL["muted"]))

    # Braces over the read row
    # A=A pairs (positions 3-8)
    out.append(brace_above(x_pos(3), x_pos(8) + BW, y_read - 4,
                           PAL["teal"], "A = A · walk in", label_dy=-9))
    # Terminal seq error (position 9)
    out.append(brace_above(x_pos(9), x_pos(9) + BW, y_read - 4,
                           PAL["red"], "seq err", label_dy=-9))
    # Pre-trimmed region (positions 10-15)
    out.append(brace_above(x_pos(n_genome), x_pos(n_total - 1) + BW, y_read - 4,
                           PAL["muted"], "removed by pre-trim", label_dy=-9))

    # Stop marker between positions 2 and 3
    stop_x = x_pos(3) - 1
    out.append(f'<line stroke="{PAL["blue"]}" stroke-width="2" '
               f'x1="{stop_x}" x2="{stop_x}" '
               f'y1="{y_read - 6}" y2="{y_read + BH + 6}"/>')

    # Walkback arrow under read row (from trimmed-read right end to stop)
    y_walk = y_read + BH + 28
    out.append(h_arrow_left(x_pos(n_genome - 1) + BW, stop_x, y_walk,
                            PAL["blue"], "walk back from 3′ end"))

    # CPA label below stop marker
    out.append(text(stop_x, y_walk + 32, "true CPA (non-A match)",
                    size=11, color=PAL["blue"], weight="700", anchor="middle"))

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
    render(build(), "walkback_readvsref")
