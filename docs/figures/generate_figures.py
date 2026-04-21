#!/usr/bin/env python3
"""
Generate polished RECTIFY README figures — unified design system.

Design principles:
  - Restrained palette: 2-3 accent colors per figure max
  - Generous whitespace between sections
  - Consistent typography (one font family, standardized sizes)
  - Subtle annotations — no garish banners or success bars
  - Clean nucleotide boxes with muted borders, accent only on key elements
  - Consistent dimensions across all figures
"""

import os
import re
import cairosvg

OUTDIR = "/sessions/eager-gifted-keller/RECTIFY/docs/figures"
os.makedirs(OUTDIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════════
# Design tokens
# ═══════════════════════════════════════════════════════════════════════════════

FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"
FIG_W = 760  # all figures same width

# Base box dimensions
BW = 26   # width
BH = 24   # height
BR = 3    # border radius

# Colors — restrained, professional palette
PAL = dict(
    # Text hierarchy
    title   = "#1e293b",
    heading = "#475569",
    label   = "#64748b",
    muted   = "#94a3b8",
    # Structural
    border  = "#cbd5e1",
    divider = "#e2e8f0",
    bg      = "#ffffff",
    # Accents (used sparingly)
    blue    = "#2563eb",
    blue_l  = "#dbeafe",
    green   = "#059669",
    green_l = "#d1fae5",
    red     = "#dc2626",
    red_l   = "#fee2e2",
    teal    = "#0d9488",
    teal_l  = "#ccfbf1",
    orange  = "#d97706",
    orange_l= "#fef3c7",
    # Exon blocks
    exon    = "#6366f1",   # indigo — distinctive, professional
    exon_t  = "#ffffff",
    exon_l  = "#e0e7ff",
    # Nucleotide fills (soft, muted)
    A_f = "#dcfce7",  A_t = "#166534",
    T_f = "#ffedd5",  T_t = "#9a3412",
    G_f = "#dbeafe",  G_t = "#1e40af",
    C_f = "#f3e8ff",  C_t = "#6b21a8",
)

BASE_FILL  = {"A": PAL["A_f"], "T": PAL["T_f"], "G": PAL["G_f"], "C": PAL["C_f"],
              "-": PAL["red_l"], "N": PAL["red_l"]}
BASE_TEXT  = {"A": PAL["A_t"], "T": PAL["T_t"], "G": PAL["G_t"], "C": PAL["C_t"],
              "-": PAL["red"], "N": PAL["red"]}


# ═══════════════════════════════════════════════════════════════════════════════
# Shared drawing helpers
# ═══════════════════════════════════════════════════════════════════════════════

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

def fig_title(text, y=22):
    return (f'<text fill="{PAL["title"]}" font-size="14" font-weight="bold" '
            f'text-anchor="middle" x="{FIG_W//2}" y="{y}">{text}</text>')

def section_head(text, y, color=None):
    """Subtle uppercase section label."""
    c = color or PAL["heading"]
    return (f'<text fill="{c}" font-size="10" font-weight="bold" '
            f'letter-spacing="0.8" x="60" y="{y}">{text}</text>')

def row_label(text, y, x=14):
    """Left-side row label, vertically centered with row."""
    return (f'<text fill="{PAL["label"]}" font-size="10" x="{x}" '
            f'y="{y}">{text}</text>')

def hdivider(y, x1=50, x2=710):
    return f'<line stroke="{PAL["divider"]}" stroke-width="1" x1="{x1}" x2="{x2}" y1="{y}" y2="{y}"/>'

def base(x, y, b, highlight_fill=None, highlight_border=None):
    """Single nucleotide box."""
    fill = highlight_fill or BASE_FILL.get(b, "#f1f5f9")
    tc   = BASE_TEXT.get(b, PAL["label"])
    bdr  = highlight_border or PAL["border"]
    sw   = "1.2" if highlight_border else "0.7"
    return (
        f'<rect fill="{fill}" height="{BH}" rx="{BR}" width="{BW}" '
        f'x="{x}" y="{y}" stroke="{bdr}" stroke-width="{sw}"/>'
        f'<text fill="{tc}" font-size="11" font-weight="600" text-anchor="middle" '
        f'x="{x + BW // 2}" y="{y + BH // 2 + 4}">{b}</text>'
    )

def brace_above(x1, x2, y, color, label):
    """Minimal bracket above a region."""
    mid = (x1 + x2) / 2
    return (
        f'<path d="M{x1},{y} L{x1},{y-5} L{x2},{y-5} L{x2},{y}" '
        f'fill="none" stroke="{color}" stroke-width="1" opacity="0.6"/>\n'
        f'<text fill="{color}" font-size="8.5" text-anchor="middle" '
        f'x="{mid}" y="{y - 9}">{label}</text>'
    )

def brace_below(x1, x2, y, color, label):
    """Minimal bracket below a region."""
    mid = (x1 + x2) / 2
    return (
        f'<path d="M{x1},{y} L{x1},{y+5} L{x2},{y+5} L{x2},{y}" '
        f'fill="none" stroke="{color}" stroke-width="1" opacity="0.6"/>\n'
        f'<text fill="{color}" font-size="8.5" text-anchor="middle" '
        f'x="{mid}" y="{y + 17}">{label}</text>'
    )

def h_arrow(x1, x2, y, color, label=None):
    """Horizontal arrow (x1->x2) with optional centered label above."""
    parts = [
        f'<line stroke="{color}" stroke-width="1.5" x1="{x1}" x2="{x2}" y1="{y}" y2="{y}"/>',
    ]
    if x2 < x1:  # leftward
        parts.append(f'<polygon fill="{color}" points="{x2},{y} {x2+7},{y-3.5} {x2+7},{y+3.5}"/>')
    else:  # rightward
        parts.append(f'<polygon fill="{color}" points="{x2},{y} {x2-7},{y-3.5} {x2-7},{y+3.5}"/>')
    if label:
        parts.append(
            f'<text fill="{color}" font-size="9" text-anchor="middle" '
            f'x="{(x1 + x2) / 2}" y="{y - 9}">{label}</text>'
        )
    return "\n".join(parts)

def vert_marker(x, y1, y2, color, label=None, label_side="right"):
    """Vertical line marker with optional label."""
    parts = [f'<line stroke="{color}" stroke-width="2" x1="{x}" x2="{x}" y1="{y1}" y2="{y2}"/>']
    if label:
        if label_side == "right":
            parts.append(
                f'<text fill="{color}" font-size="9" font-weight="600" '
                f'x="{x + 5}" y="{y1 - 4}">{label}</text>'
            )
        else:
            parts.append(
                f'<text fill="{color}" font-size="9" font-weight="600" '
                f'text-anchor="end" x="{x - 5}" y="{y1 - 4}">{label}</text>'
            )
    return "\n".join(parts)

def exon(x, y, w, h, label, color=None):
    """Exon block."""
    c = color or PAL["exon"]
    return (
        f'<rect fill="{c}" height="{h}" rx="4" width="{w}" x="{x}" y="{y}"/>\n'
        f'<text fill="{PAL["exon_t"]}" font-size="11" font-weight="600" '
        f'text-anchor="middle" x="{x + w / 2}" y="{y + h / 2 + 4}">{label}</text>'
    )

def intron(x1, x2, y, label="intron"):
    """Dashed intron line."""
    return (
        f'<line stroke="{PAL["muted"]}" stroke-dasharray="6,4" stroke-width="1.5" '
        f'x1="{x1}" x2="{x2}" y1="{y}" y2="{y}"/>\n'
        f'<text fill="{PAL["muted"]}" font-size="8.5" text-anchor="middle" '
        f'x="{(x1 + x2) / 2}" y="{y - 7}">{label}</text>'
    )

def softclip_block(x, y, w, h, label="soft-clipped"):
    """Dashed red-tinted block for soft-clipped region."""
    return (
        f'<rect fill="{PAL["red_l"]}" height="{h}" rx="4" '
        f'stroke="{PAL["red"]}" stroke-dasharray="4,3" stroke-width="1" '
        f'width="{w}" x="{x}" y="{y}"/>\n'
        f'<text fill="{PAL["red"]}" font-size="9" text-anchor="middle" '
        f'x="{x + w / 2}" y="{y + h / 2 + 3}">{label}</text>'
    )

def aligned_block(x, y, w, h, label, color=None):
    """Semi-transparent aligned region block."""
    c = color or PAL["exon"]
    return (
        f'<rect fill="{c}" height="{h}" opacity="0.18" rx="4" '
        f'stroke="{c}" stroke-width="1" width="{w}" x="{x}" y="{y}"/>\n'
        f'<text fill="{c}" font-size="9" font-weight="600" text-anchor="middle" '
        f'x="{x + w / 2}" y="{y + h / 2 + 3}">{label}</text>'
    )

def rescued_block(x, y, w, h, label):
    """Teal-tinted rescued region."""
    return (
        f'<rect fill="{PAL["teal_l"]}" height="{h}" rx="4" '
        f'stroke="{PAL["teal"]}" stroke-width="1.2" '
        f'width="{w}" x="{x}" y="{y}"/>\n'
        f'<text fill="{PAL["teal"]}" font-size="9" font-weight="600" text-anchor="middle" '
        f'x="{x + w / 2}" y="{y + h / 2 + 3}">{label}</text>'
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: 3' End Indel Correction
# ═══════════════════════════════════════════════════════════════════════════════

def fig_indel_correction():
    """
    3' End Walk-Back on a pre-trimmed read.  After poly(A) pre-trimming,
    the read ends near the CPA site but the aligner may still misplace the
    3' end within the genomic A-tract (indels, T basecalling error).
    Walk-back resolves the true CPA position using genome context.
    """
    H = 370
    X0 = 70  # base start x
    L = []
    L.append(svg_open(H))
    L.append(fig_title("3\u2032 End Walk-Back: Recovering the True CPA Site"))

    # -- BEFORE --
    y_sec = 42
    L.append(section_head("ALIGNED READ (pre-trimmed, poly(A) already removed)", y_sec, PAL["red"]))

    # Genome row
    # Genome: C A G T C | A A A A A G A A A A | G T C
    #         0 1 2 3 4   5 6 7 8 9 10 11 12 13 14   15 16 17
    # CPA site is at position 5 (first A of the A-tract)
    y_g = y_sec + 28
    L.append(row_label("Genome", y_g + BH // 2 + 4))
    genome = list("CAGTCAAAAAGAAAAGTC")
    for i, b in enumerate(genome):
        hf = None
        hb = None
        if 5 <= i <= 14 and b == "A":
            hb = PAL["green"]
        if i == 10 and b == "G":
            hb = PAL["red"]
            hf = PAL["red_l"]
        if i >= 15:
            hb = PAL["orange"]
        L.append(base(X0 + i * BW, y_g, b, hf, hb))

    # Brackets above genome
    L.append(brace_above(X0 + 5 * BW, X0 + 15 * BW, y_g - 1, PAL["green"], "genomic A-tract"))
    L.append(brace_above(X0 + 15 * BW, X0 + 18 * BW, y_g - 1, PAL["orange"], "non-A"))

    # CPA site marker
    cpa_x = X0 + 5 * BW
    L.append(vert_marker(cpa_x, y_g - 14, y_g, PAL["green"], "CPA site", "left"))

    # Read row — pre-trimmed read aligned to genome
    # After pre-trimming: read ends with ...CAGTCAAAAAT (T is ambiguous base)
    # Aligner places: CAGTC matches, then AAAA aligns into A-tract,
    # then deletion at G (pos 10), then AAA aligns to second A-tract,
    # then T mismatches against A (pos 14)
    y_r = y_g + 50
    L.append(row_label("Read", y_r + BH // 2 + 4))
    # Transcript body: CAGTC
    for i, b in enumerate(list("CAGTC")):
        L.append(base(X0 + i * BW, y_r, b))

    # A's aligned into first A-tract
    for i in range(5):
        L.append(base(X0 + (5 + i) * BW, y_r, "A"))

    # Deletion at position 10 (G in genome)
    L.append(base(X0 + 10 * BW, y_r, "-", PAL["red_l"], PAL["red"]))

    # A's aligned into second A-tract
    for i in range(3):
        L.append(base(X0 + (11 + i) * BW, y_r, "A"))

    # T basecalling error (was poly(A), pre-trimmer left it)
    L.append(base(X0 + 14 * BW, y_r, "T", PAL["orange_l"], PAL["orange"]))

    # Labels
    L.append(brace_below(X0 + 5 * BW, X0 + 15 * BW, y_r + BH + 1, PAL["muted"],
                          "read extends past CPA into A-tract (del + T error from pre-trim)"))

    # Apparent vs true 3' end markers
    app_x = X0 + 15 * BW
    L.append(vert_marker(app_x, y_r - 4, y_r + BH + 4, PAL["red"], "apparent 3\u2032 end", "right"))

    # -- DIVIDER --
    y_div = y_r + BH + 30
    L.append(hdivider(y_div))

    # -- AFTER --
    y_sec2 = y_div + 18
    L.append(section_head("AFTER WALK-BACK", y_sec2, PAL["green"]))

    # Walk-back arrow
    y_arrow = y_sec2 + 18
    arrow_r = X0 + 14 * BW
    arrow_l = cpa_x + 4
    L.append(h_arrow(arrow_r, arrow_l, y_arrow, PAL["blue"],
                     "walk back: skip A\u2019s, deletion, T error \u2192 stop at C"))

    # Corrected row
    y_c = y_arrow + 24
    L.append(row_label("Corrected", y_c + BH // 2 + 4))
    for i, b in enumerate(list("CAGTC")):
        L.append(base(X0 + i * BW, y_c, b))

    # Green marker at true 3' end
    L.append(vert_marker(cpa_x, y_c - 8, y_c + BH + 4, PAL["green"], "true 3\u2032 end (CPA)"))

    # Note: no poly(A) to show — it was already trimmed
    L.append(f'<text fill="{PAL["muted"]}" font-size="9" font-style="italic" '
             f'x="{cpa_x + 8}" y="{y_c + BH + 22}">'
             f'poly(A) was removed during pre-trimming; CPA position recorded in corrected TSV</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: Soft-Clip Rescue
# ═══════════════════════════════════════════════════════════════════════════════

def fig_softclip_rescue():
    H = 350
    X0 = 70
    L = []
    L.append(svg_open(H))
    L.append(fig_title("Soft-Clip Rescue at Homopolymer Boundaries"))

    # -- BEFORE --
    y_sec = 42
    L.append(section_head("BEFORE CORRECTION", y_sec, PAL["red"]))

    # Reference row
    y_ref = y_sec + 20
    L.append(row_label("Reference", y_ref + BH // 2 + 4))
    for i in range(11):
        L.append(base(X0 + i * BW, y_ref, "T", None, PAL["orange"]))
    L.append(base(X0 + 11 * BW, y_ref, "C"))
    L.append(base(X0 + 12 * BW, y_ref, "G"))
    L.append(base(X0 + 13 * BW, y_ref, "A"))
    L.append(brace_above(X0, X0 + 11 * BW, y_ref - 1, PAL["orange"], "T-tract (11 bp)"))

    # Read row
    y_r = y_ref + 48
    L.append(row_label("Read", y_r + BH // 2 + 4))
    for i in range(8):
        L.append(base(X0 + i * BW, y_r, "T"))
    L.append(base(X0 + 8 * BW, y_r, "C", PAL["red_l"], PAL["red"]))
    for i in range(9):
        L.append(base(X0 + (9 + i) * BW, y_r, "A", PAL["red_l"], PAL["red"]))

    # Soft-clip boundary
    sc_x = X0 + 8 * BW
    L.append(f'<line stroke="{PAL["red"]}" stroke-dasharray="4,3" stroke-width="1.5" '
             f'x1="{sc_x}" x2="{sc_x}" y1="{y_r - 6}" y2="{y_r + BH + 4}"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="9" x="{sc_x + 4}" '
             f'y="{y_r - 9}">apparent 3\u2032 end (wrong)</text>')

    L.append(f'<text fill="{PAL["orange"]}" font-size="8.5" text-anchor="middle" '
             f'x="{X0 + 4 * BW}" y="{y_r + BH + 14}">8 T\u2019s called (should be 11)</text>')

    # -- DIVIDER --
    y_div = y_r + BH + 26
    L.append(hdivider(y_div))

    # -- AFTER --
    y_sec2 = y_div + 18
    L.append(section_head("AFTER CORRECTION", y_sec2, PAL["green"]))

    # Reference row 2
    y_ref2 = y_sec2 + 20
    L.append(row_label("Reference", y_ref2 + BH // 2 + 4))
    for i in range(11):
        L.append(base(X0 + i * BW, y_ref2, "T", None, PAL["orange"]))
    L.append(base(X0 + 11 * BW, y_ref2, "C"))
    L.append(base(X0 + 12 * BW, y_ref2, "G"))
    L.append(base(X0 + 13 * BW, y_ref2, "A"))

    # Corrected row
    y_c = y_ref2 + 40
    L.append(row_label("Corrected", y_c + BH // 2 + 4))
    for i in range(8):
        L.append(base(X0 + i * BW, y_c, "T"))
    for i in range(3):
        L.append(base(X0 + (8 + i) * BW, y_c, "T", PAL["teal_l"], PAL["teal"]))
    L.append(base(X0 + 11 * BW, y_c, "C", PAL["teal_l"], PAL["teal"]))
    for i in range(9):
        L.append(base(X0 + (12 + i) * BW, y_c, "A", PAL["blue_l"], PAL["blue"]))

    # Green marker at corrected 3' end
    corr_x = X0 + 12 * BW
    L.append(vert_marker(corr_x, y_c - 8, y_c + BH + 4, PAL["green"], "corrected 3\u2032 end"))

    # Brackets
    L.append(brace_below(X0 + 8 * BW, X0 + 12 * BW, y_c + BH + 1, PAL["teal"],
                          "+4 bp rescued (3 T\u2019s + CPA base)"))
    L.append(brace_below(X0 + 12 * BW, X0 + 21 * BW, y_c + BH + 1, PAL["blue"],
                          "poly(A) tail"))

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: 5' End Correction (Junction Soft-Clips)
# ═══════════════════════════════════════════════════════════════════════════════

def fig_5prime_junction():
    """
    Two-panel figure for Section 2: Unified 5' End and Junction Correction.
    Part 1 — Cat3: 5' soft-clip junction rescue (before/after).
    Part 2 — Module 2H: Post-consensus N-op junction refinement (before/after).
    """
    H = 560
    BLK_H = 26
    L = []
    L.append(svg_open(H))
    L.append(fig_title("Unified 5\u2032 End and Junction Correction"))

    # ════════════════════════════════════════════════════════════════
    # PART 1 — Cat3: 5' Soft-Clip Junction Rescue
    # ════════════════════════════════════════════════════════════════
    y_sec = 40
    L.append(section_head("PART 1 \u2014 5\u2032 SOFT-CLIP JUNCTION RESCUE (Cat3)", y_sec, PAL["blue"]))

    # Genome
    y_g = y_sec + 22
    L.append(row_label("Genome", y_g + BLK_H // 2 + 4))
    L.append(exon(90, y_g, 160, BLK_H, "Exon 1"))
    L.append(intron(250, 370, y_g + BLK_H // 2))
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" x="253" y="{y_g + BLK_H - 1}">GT</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" text-anchor="end" x="367" y="{y_g + BLK_H - 1}">AG</text>')
    L.append(exon(370, y_g, 160, BLK_H, "Exon 2"))

    # Read with soft-clip (BEFORE)
    y_r = y_g + 38
    L.append(row_label("Read", y_r + BLK_H // 2 + 4))
    L.append(softclip_block(170, y_r, 200, BLK_H, "soft-clipped (matches Exon 1)"))
    L.append(aligned_block(370, y_r, 160, BLK_H, "aligned to Exon 2", PAL["exon"]))
    # Annotation: NW alignment
    L.append(f'<text fill="{PAL["blue"]}" font-size="8" font-style="italic" '
             f'x="{540}" y="{y_r + BLK_H//2 - 2}">semi-global NW alignment</text>')
    L.append(f'<text fill="{PAL["blue"]}" font-size="8" font-style="italic" '
             f'x="{540}" y="{y_r + BLK_H//2 + 10}">places bases in Exon 1</text>')

    # Arrow showing rescue
    y_arr1 = y_r + BLK_H + 10
    L.append(f'<line stroke="{PAL["teal"]}" stroke-width="1.5" stroke-dasharray="4,3" '
             f'x1="270" x2="170" y1="{y_arr1}" y2="{y_arr1}"/>')
    L.append(f'<polygon fill="{PAL["teal"]}" points="170,{y_arr1} 177,{y_arr1-3.5} 177,{y_arr1+3.5}"/>')
    L.append(f'<text fill="{PAL["teal"]}" font-size="8" font-weight="600" '
             f'x="280" y="{y_arr1 + 4}">rescue \u2192 Exon 1</text>')

    # Corrected read (AFTER)
    y_c = y_arr1 + 14
    L.append(row_label("Corrected", y_c + BLK_H // 2 + 4))
    L.append(rescued_block(100, y_c, 150, BLK_H, "rescued \u2192 Exon 1"))
    L.append(intron(250, 370, y_c + BLK_H // 2, "N (splice)"))
    L.append(aligned_block(370, y_c, 160, BLK_H, "Exon 2", PAL["exon"]))
    L.append(f'<text fill="{PAL["green"]}" font-size="9" font-weight="600" '
             f'x="{540}" y="{y_c + BLK_H//2 + 4}">\u2713 true 5\u2032 end recovered</text>')

    # ════════════════════════════════════════════════════════════════
    # DIVIDER between panels
    # ════════════════════════════════════════════════════════════════
    y_div = y_c + BLK_H + 18
    L.append(hdivider(y_div))

    # ════════════════════════════════════════════════════════════════
    # PART 2 — Module 2H: Post-Consensus N-Op Junction Refinement
    # ════════════════════════════════════════════════════════════════
    y_sec2 = y_div + 16
    L.append(section_head("PART 2 \u2014 N-OP JUNCTION REFINEMENT (Module 2H)", y_sec2, PAL["teal"]))

    # Genome with two candidate junctions
    y_g2 = y_sec2 + 22
    L.append(row_label("Genome", y_g2 + BLK_H // 2 + 4))
    L.append(exon(90, y_g2, 140, BLK_H, "Exon 1"))
    # Show two possible introns: annotated (dashed) and true (solid candidate)
    ann_int_x1 = 230  # annotated junction 5'SS
    ann_int_x2 = 380  # annotated junction 3'SS
    true_int_x1 = 240  # true junction 5'SS (shifted)
    true_int_x2 = 370  # true junction 3'SS (shifted)
    L.append(intron(ann_int_x1, ann_int_x2, y_g2 + BLK_H // 2))
    L.append(f'<text fill="{PAL["muted"]}" font-size="7" x="{ann_int_x1 + 2}" '
             f'y="{y_g2 + BLK_H - 1}">GT</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="7" text-anchor="end" '
             f'x="{ann_int_x2 - 2}" y="{y_g2 + BLK_H - 1}">AG</text>')
    L.append(exon(380, y_g2, 150, BLK_H, "Exon 2"))

    # Read BEFORE (aligned to annotated junction, but with mismatches)
    y_r2 = y_g2 + 38
    L.append(row_label("Read", y_r2 + BLK_H // 2 + 4))
    L.append(aligned_block(90, y_r2, 140, BLK_H, "Exon 1", PAL["exon"]))
    # N-op at annotated position
    L.append(f'<line stroke="{PAL["muted"]}" stroke-dasharray="5,3" stroke-width="1.2" '
             f'x1="{ann_int_x1}" x2="{ann_int_x2}" y1="{y_r2+BLK_H//2}" y2="{y_r2+BLK_H//2}"/>')
    L.append(aligned_block(380, y_r2, 150, BLK_H, "Exon 2", PAL["exon"]))
    # Mismatch indicators at junction boundaries
    for dx in [-8, -4, 4, 8]:
        mx = ann_int_x2 + dx
        L.append(f'<rect fill="{PAL["red"]}" width="3" height="3" rx="1" '
                 f'x="{mx}" y="{y_r2 + BLK_H + 2}" opacity="0.7"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="7" text-anchor="middle" '
             f'x="{ann_int_x2}" y="{y_r2 + BLK_H + 14}">junction-proximal mismatches</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'x="540" y="{y_r2 + BLK_H//2 + 4}">annotated N-op</text>')

    # Arrow showing re-scoring
    y_arr2 = y_r2 + BLK_H + 22
    L.append(f'<text fill="{PAL["teal"]}" font-size="8" font-weight="600" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_arr2 + 4}">'
             f'HP-aware re-scoring: sequence match > stability > GT-AG > annotation > shift</text>')

    # Read AFTER (junction shifted to true position)
    y_c2 = y_arr2 + 14
    L.append(row_label("Corrected", y_c2 + BLK_H // 2 + 4))
    L.append(aligned_block(90, y_c2, 150, BLK_H, "Exon 1", PAL["exon"]))
    # N-op at refined position (shifted)
    L.append(f'<line stroke="{PAL["teal"]}" stroke-dasharray="5,3" stroke-width="1.5" '
             f'x1="{true_int_x1}" x2="{true_int_x2}" y1="{y_c2+BLK_H//2}" y2="{y_c2+BLK_H//2}"/>')
    L.append(f'<text fill="{PAL["teal"]}" font-size="7" font-weight="600" '
             f'x="{true_int_x1 + 2}" y="{y_c2 + BLK_H - 1}">GT</text>')
    L.append(f'<text fill="{PAL["teal"]}" font-size="7" font-weight="600" '
             f'text-anchor="end" x="{true_int_x2 - 2}" y="{y_c2 + BLK_H - 1}">AG</text>')
    L.append(aligned_block(370, y_c2, 160, BLK_H, "Exon 2", PAL["exon"]))

    # Show boundary shift arrows
    L.append(f'<line stroke="{PAL["teal"]}" stroke-width="1" '
             f'x1="{ann_int_x1}" x2="{true_int_x1}" '
             f'y1="{y_c2 - 4}" y2="{y_c2 - 4}"/>')
    L.append(f'<polygon fill="{PAL["teal"]}" points="{true_int_x1},{y_c2-4} '
             f'{true_int_x1-4},{y_c2-7} {true_int_x1-4},{y_c2-1}"/>')
    L.append(f'<line stroke="{PAL["teal"]}" stroke-width="1" '
             f'x1="{ann_int_x2}" x2="{true_int_x2}" '
             f'y1="{y_c2 - 4}" y2="{y_c2 - 4}"/>')
    L.append(f'<polygon fill="{PAL["teal"]}" points="{true_int_x2},{y_c2-4} '
             f'{true_int_x2+4},{y_c2-7} {true_int_x2+4},{y_c2-1}"/>')

    L.append(f'<text fill="{PAL["green"]}" font-size="9" font-weight="600" '
             f'x="540" y="{y_c2 + BLK_H//2 + 4}">\u2713 sequence-optimal junction</text>')

    # Priority note
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_c2 + BLK_H + 18}">'
             f'sequence evidence always overrides annotation; '
             f'fast path skips reads already at annotated canonical junctions</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: Multi-Aligner Consensus
# ═══════════════════════════════════════════════════════════════════════════════

def fig_multi_aligner_consensus():
    """
    Three-stage multi-aligner rectification pipeline matching the README:
      Stage 1 — Per-aligner rectification (3 parallel lanes)
      Stage 2 — Consensus selection (priority chain)
      Stage 3 — Chimeric reconstruction (optional)
    """
    H = 680
    BLK_H = 22
    L = []
    L.append(svg_open(H))
    L.append(fig_title("Multi-Aligner Rectification Pipeline"))

    # Aligner colors
    AC = {"minimap2": "#3b82f6", "mapPacBio": "#f59e0b", "gapmm2": "#22c55e"}
    ALIGNERS = ["minimap2", "mapPacBio", "gapmm2"]

    # Layout constants
    LANE_X = 100          # left edge of lane blocks
    LANE_W = 180          # width of each lane
    LANE_GAP = 20         # gap between lanes
    CMD_H = 26            # height of command boxes
    TSV_H = 20            # height of TSV output boxes

    def lane_x(i):
        return LANE_X + i * (LANE_W + LANE_GAP)

    def cmd_box(x, y, w, h, label, color):
        """Rounded command/process box."""
        return (
            f'<rect fill="{color}" height="{h}" rx="5" width="{w}" '
            f'x="{x}" y="{y}" opacity="0.15"/>\n'
            f'<rect fill="none" height="{h}" rx="5" width="{w}" '
            f'x="{x}" y="{y}" stroke="{color}" stroke-width="1.2"/>\n'
            f'<text fill="{color}" font-size="9" font-weight="600" '
            f'text-anchor="middle" x="{x + w/2}" y="{y + h/2 + 3.5}">{label}</text>'
        )

    def tsv_box(x, y, w, h, label, color):
        """File output box (smaller, muted fill)."""
        return (
            f'<rect fill="{color}" height="{h}" rx="3" width="{w}" '
            f'x="{x}" y="{y}" opacity="0.10"/>\n'
            f'<rect fill="none" height="{h}" rx="3" width="{w}" '
            f'x="{x}" y="{y}" stroke="{color}" stroke-width="0.8"/>\n'
            f'<text fill="{color}" font-size="8" font-weight="600" '
            f'text-anchor="middle" x="{x + w/2}" y="{y + h/2 + 3}">{label}</text>'
        )

    def v_arrow(x, y1, y2, color):
        """Vertical downward arrow."""
        return (
            f'<line stroke="{color}" stroke-width="1.3" x1="{x}" x2="{x}" '
            f'y1="{y1}" y2="{y2}"/>\n'
            f'<polygon fill="{color}" points="{x},{y2} {x-3.5},{y2-6} {x+3.5},{y2-6}"/>'
        )

    # ════════════════════════════════════════════════════════════════
    # STAGE 1 — Per-aligner rectification
    # ════════════════════════════════════════════════════════════════
    y_sec1 = 40
    L.append(section_head("STAGE 1 \u2014 PER-ALIGNER RECTIFICATION", y_sec1, PAL["blue"]))

    # BAM inputs
    y_bam = y_sec1 + 20
    for i, aligner in enumerate(ALIGNERS):
        cx = lane_x(i)
        c = AC[aligner]
        # Aligner name
        L.append(f'<text fill="{c}" font-size="10" font-weight="700" '
                 f'text-anchor="middle" x="{cx + LANE_W/2}" y="{y_bam + 4}">{aligner}</text>')
        # BAM box
        y_box = y_bam + 10
        L.append(tsv_box(cx, y_box, LANE_W, TSV_H, f"{aligner}.bam", c))
        # Arrow down to rectify correct
        L.append(v_arrow(cx + LANE_W/2, y_box + TSV_H, y_box + TSV_H + 16, c))
        # rectify correct box
        y_cmd = y_box + TSV_H + 16
        L.append(cmd_box(cx, y_cmd, LANE_W, CMD_H, "rectify correct", c))
        # Arrow down to corrected TSV
        L.append(v_arrow(cx + LANE_W/2, y_cmd + CMD_H, y_cmd + CMD_H + 16, c))
        # corrected_3ends.tsv output
        y_tsv = y_cmd + CMD_H + 16
        L.append(tsv_box(cx, y_tsv, LANE_W, TSV_H, "corrected_3ends.tsv", c))

    # Side annotation — what rectify correct does (right-aligned to fit)
    ann_x = FIG_W - 12
    y_ann = y_bam + 10 + TSV_H + 16
    L.append(f'<text fill="{PAL["muted"]}" font-size="7.5" font-style="italic" '
             f'text-anchor="end" x="{ann_x}" y="{y_ann + 6}">3\u2032 walk-back</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="7.5" font-style="italic" '
             f'text-anchor="end" x="{ann_x}" y="{y_ann + 16}">5\u2032 junction rescue</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="7.5" font-style="italic" '
             f'text-anchor="end" x="{ann_x}" y="{y_ann + 26}">soft-clip rescue</text>')

    # ════════════════════════════════════════════════════════════════
    # STAGE 2 — Consensus selection
    # ════════════════════════════════════════════════════════════════
    y_tsv_bottom = y_bam + 10 + TSV_H + 16 + CMD_H + 16 + TSV_H
    y_div1 = y_tsv_bottom + 16
    L.append(hdivider(y_div1))
    y_sec2 = y_div1 + 14
    L.append(section_head("STAGE 2 \u2014 CONSENSUS SELECTION", y_sec2, PAL["green"]))

    # Converging arrows from 3 TSVs into central consensus box
    y_conv = y_sec2 + 18
    cons_w = 260
    cons_x = (FIG_W - cons_w) / 2
    cons_y = y_conv + 26

    # Draw converging lines from each lane's TSV to the consensus box
    for i in range(3):
        src_x = lane_x(i) + LANE_W / 2
        dst_x = cons_x + cons_w / 2
        # Diagonal line from TSV bottom to consensus top
        L.append(f'<line stroke="{AC[ALIGNERS[i]]}" stroke-width="1.2" '
                 f'x1="{src_x}" x2="{dst_x}" y1="{y_conv}" y2="{cons_y}" opacity="0.5"/>')

    # Arrow head at consensus box
    mid_x = cons_x + cons_w / 2
    L.append(f'<polygon fill="{PAL["green"]}" points="{mid_x},{cons_y} '
             f'{mid_x-4},{cons_y-7} {mid_x+4},{cons_y-7}"/>')

    # rectify consensus box (wider, more prominent)
    L.append(cmd_box(cons_x, cons_y, cons_w, CMD_H + 4, "rectify consensus", PAL["green"]))

    # Priority chain below the consensus box
    y_pri = cons_y + CMD_H + 4 + 14
    L.append(f'<text fill="{PAL["heading"]}" font-size="9" font-weight="600" '
             f'text-anchor="middle" x="{FIG_W/2}" y="{y_pri}">per-read winner selection (priority order):</text>')

    # Priority items as a horizontal chain with arrows
    y_chain = y_pri + 16
    priorities = ["5\u2032 rescued", "confidence", "3\u2032 agreement", "span", "n_junctions"]
    # Calculate total width
    pri_widths = [68, 70, 82, 42, 72]
    pri_gap = 8
    total_pri_w = sum(pri_widths) + (len(priorities) - 1) * pri_gap
    px = (FIG_W - total_pri_w) / 2
    pri_h = 20

    for idx, (pri, pw) in enumerate(zip(priorities, pri_widths)):
        # Numbered priority box
        is_first = (idx == 0)
        fill_op = "0.12" if not is_first else "0.20"
        bdr = PAL["green"] if is_first else PAL["heading"]
        tc = PAL["green"] if is_first else PAL["heading"]
        L.append(
            f'<rect fill="{bdr}" height="{pri_h}" rx="10" width="{pw}" '
            f'x="{px}" y="{y_chain}" opacity="{fill_op}"/>\n'
            f'<rect fill="none" height="{pri_h}" rx="10" width="{pw}" '
            f'x="{px}" y="{y_chain}" stroke="{bdr}" stroke-width="0.8"/>\n'
            f'<text fill="{tc}" font-size="8" font-weight="600" '
            f'text-anchor="middle" x="{px + pw/2}" y="{y_chain + pri_h/2 + 3}">{pri}</text>'
        )
        # Arrow between items
        if idx < len(priorities) - 1:
            ax1 = px + pw + 1
            ax2 = ax1 + pri_gap - 2
            ay = y_chain + pri_h / 2
            L.append(f'<line stroke="{PAL["muted"]}" stroke-width="1" '
                     f'x1="{ax1}" x2="{ax2}" y1="{ay}" y2="{ay}"/>')
            L.append(f'<polygon fill="{PAL["muted"]}" points="{ax2},{ay} '
                     f'{ax2-4},{ay-2.5} {ax2-4},{ay+2.5}"/>')
        px += pw + pri_gap

    # Arrow from consensus to output
    y_out_arr = y_chain + pri_h + 10
    L.append(v_arrow(FIG_W/2, y_out_arr, y_out_arr + 18, PAL["green"]))

    # Output: consensus corrected_3ends.tsv
    y_out = y_out_arr + 18
    out_w = 200
    out_x = (FIG_W - out_w) / 2
    L.append(tsv_box(out_x, y_out, out_w, TSV_H + 2, "corrected_3ends.tsv (consensus)", PAL["green"]))

    # XA/XC/XN tags annotation
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'x="{out_x + out_w + 10}" y="{y_out + (TSV_H+2)/2 + 3}">'
             f'tags: XA (winner), XC (confidence), XN (agreement)</text>')

    # ════════════════════════════════════════════════════════════════
    # STAGE 3 — Chimeric reconstruction (optional)
    # ════════════════════════════════════════════════════════════════
    y_div2 = y_out + TSV_H + 2 + 16
    L.append(hdivider(y_div2))
    y_sec3 = y_div2 + 14
    L.append(section_head("STAGE 3 \u2014 CHIMERIC RECONSTRUCTION (optional)", y_sec3, PAL["teal"]))

    # Show two aligners each contributing a unique junction
    # Genome: Exon A -- intron1 -- Exon B -- intron2 -- Exon C (3 exons, 2 junctions)
    EXA_X, EXA_W = 90, 110
    EXB_X, EXB_W = 250, 100
    EXC_X, EXC_W = 410, 110
    y_gen = y_sec3 + 20
    L.append(row_label("Genome", y_gen + BLK_H//2 + 4, x=14))
    L.append(exon(EXA_X, y_gen, EXA_W, BLK_H, "Exon A"))
    L.append(intron(EXA_X + EXA_W, EXB_X, y_gen + BLK_H//2, "intron 1"))
    L.append(exon(EXB_X, y_gen, EXB_W, BLK_H, "Exon B"))
    L.append(intron(EXB_X + EXB_W, EXC_X, y_gen + BLK_H//2, "intron 2"))
    L.append(exon(EXC_X, y_gen, EXC_W, BLK_H, "Exon C"))

    # Aligner 1 (minimap2): has junction 1 (A-B) but misses junction 2
    y_a1 = y_gen + 32
    c1 = AC["minimap2"]
    L.append(f'<text fill="{c1}" font-size="9" font-weight="600" x="14" '
             f'y="{y_a1 + BLK_H//2 + 3}">minimap2</text>')
    L.append(aligned_block(EXA_X, y_a1, EXA_W, BLK_H, "A", c1))
    L.append(f'<line stroke="{c1}" stroke-dasharray="4,3" stroke-width="1" '
             f'x1="{EXA_X+EXA_W}" x2="{EXB_X}" y1="{y_a1+BLK_H//2}" y2="{y_a1+BLK_H//2}"/>')
    # B + C merged (no junction 2 found)
    L.append(aligned_block(EXB_X, y_a1, EXB_W + (EXC_X - EXB_X - EXB_W) + EXC_W, BLK_H, "B + C (unspliced)", c1))
    L.append(f'<text fill="{PAL["green"]}" font-size="8" x="{EXC_X + EXC_W + 10}" '
             f'y="{y_a1 + BLK_H//2 - 2}">\u2713 junction 1</text>')
    L.append(f'<text fill="{PAL["red"]}" font-size="8" x="{EXC_X + EXC_W + 10}" '
             f'y="{y_a1 + BLK_H//2 + 10}">\u2717 junction 2</text>')

    # Aligner 2 (mapPacBio): misses junction 1 but has junction 2
    y_a2 = y_a1 + BLK_H + 8
    c2 = AC["mapPacBio"]
    L.append(f'<text fill="{c2}" font-size="9" font-weight="600" x="14" '
             f'y="{y_a2 + BLK_H//2 + 3}">mapPacBio</text>')
    # A + B merged (no junction 1 found)
    L.append(aligned_block(EXA_X, y_a2, EXA_W + (EXB_X - EXA_X - EXA_W) + EXB_W, BLK_H, "A + B (unspliced)", c2))
    L.append(f'<line stroke="{c2}" stroke-dasharray="4,3" stroke-width="1" '
             f'x1="{EXB_X+EXB_W}" x2="{EXC_X}" y1="{y_a2+BLK_H//2}" y2="{y_a2+BLK_H//2}"/>')
    L.append(aligned_block(EXC_X, y_a2, EXC_W, BLK_H, "C", c2))
    L.append(f'<text fill="{PAL["red"]}" font-size="8" x="{EXC_X + EXC_W + 10}" '
             f'y="{y_a2 + BLK_H//2 - 2}">\u2717 junction 1</text>')
    L.append(f'<text fill="{PAL["green"]}" font-size="8" x="{EXC_X + EXC_W + 10}" '
             f'y="{y_a2 + BLK_H//2 + 10}">\u2713 junction 2</text>')

    # Arrow down to chimeric result
    y_arrow_ch = y_a2 + BLK_H + 6
    L.append(v_arrow(FIG_W/2, y_arrow_ch, y_arrow_ch + 16, PAL["teal"]))
    L.append(f'<text fill="{PAL["teal"]}" font-size="8" font-weight="600" '
             f'x="{FIG_W/2 + 8}" y="{y_arrow_ch + 10}">stitch</text>')

    # Chimeric result: all three exons properly spliced
    y_ch = y_arrow_ch + 16
    L.append(f'<text fill="{PAL["teal"]}" font-size="9" font-weight="600" x="14" '
             f'y="{y_ch + BLK_H//2 + 3}">chimeric</text>')
    L.append(rescued_block(EXA_X, y_ch, EXA_W, BLK_H, "A"))
    L.append(f'<line stroke="{PAL["teal"]}" stroke-dasharray="4,3" stroke-width="1.2" '
             f'x1="{EXA_X+EXA_W}" x2="{EXB_X}" y1="{y_ch+BLK_H//2}" y2="{y_ch+BLK_H//2}"/>')
    L.append(rescued_block(EXB_X, y_ch, EXB_W, BLK_H, "B"))
    L.append(f'<line stroke="{PAL["teal"]}" stroke-dasharray="4,3" stroke-width="1.2" '
             f'x1="{EXB_X+EXB_W}" x2="{EXC_X}" y1="{y_ch+BLK_H//2}" y2="{y_ch+BLK_H//2}"/>')
    L.append(rescued_block(EXC_X, y_ch, EXC_W, BLK_H, "C"))
    L.append(f'<text fill="{PAL["green"]}" font-size="8" font-weight="600" x="{EXC_X + EXC_W + 10}" '
             f'y="{y_ch + BLK_H//2 + 3}">\u2713 both junctions recovered</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 6: DRS Poly(A) Pre-Trimming
# ═══════════════════════════════════════════════════════════════════════════════

def fig_polya_pretrim():
    """
    Illustrates the two-pass poly(A) + adapter trimming algorithm used by
    `rectify trim-polya` before re-alignment.  The trimmer is deliberately
    conservative: Pass 0 strips the adapter stub, Pass 1 strips only
    consecutive pure A's from the 3' end, stopping at the first non-A base.
    The ambiguous boundary (terminal T that could be genomic or a seq error)
    is left for the post-alignment walk-back step to resolve.
    """
    H = 390
    X0 = 70   # left margin for bases
    L = []
    L.append(svg_open(H))
    L.append(fig_title("DRS Poly(A) Pre-Trimming"))

    # ── Color shortcuts ──
    ADAPTER_F = PAL["orange_l"]
    ADAPTER_B = PAL["orange"]
    POLYA_F   = PAL["blue_l"]
    POLYA_B   = PAL["blue"]
    AMB_F     = "#fef9c3"       # soft yellow for ambiguous base
    AMB_B     = "#ca8a04"       # darker yellow border
    BODY_C    = PAL["exon"]     # indigo for transcript body block

    # ── Sequence layout ──
    # Read in RNA 5'→3': [transcript body]...C T G A A A A A T A A A A [adapter]
    # The T between the A-tracts could be a basecalling error in the tail,
    # or the last genomic base before the CPA site.
    BODY_BLK_W = 160           # wide block representing transcript body
    body_tail  = list("CTG")   # last 3 genomic body bases
    a_tract    = list("AAAAA") # genomic A-tract / remaining tail A's
    amb_t      = "T"           # ambiguous base (seq error or genomic?)
    pure_a     = list("AAAA")  # pure poly(A) tail (will be trimmed)
    adapter    = list("TCT")   # adapter stub

    n_body   = len(body_tail)
    n_atract = len(a_tract)
    n_pure   = len(pure_a)
    n_adapt  = len(adapter)

    bx_body  = X0 + BODY_BLK_W + 4
    bx_at    = bx_body + n_body * BW
    bx_t     = bx_at + n_atract * BW
    bx_pure  = bx_t + BW
    bx_adapt = bx_pure + n_pure * BW

    # ── Helper: draw bases with optional fading ──
    def bases_row(x, y, seq, fill, border, text_c, fade=False):
        parts = []
        op = ' opacity="0.25"' if fade else ""
        for i, b in enumerate(seq):
            parts.append(
                f'<rect fill="{fill}" height="{BH}" rx="{BR}" width="{BW}" '
                f'x="{x + i*BW}" y="{y}" stroke="{border}" stroke-width="0.7"{op}/>'
                f'<text fill="{text_c}" font-size="11" font-weight="600" text-anchor="middle" '
                f'x="{x + i*BW + BW//2}" y="{y + BH//2 + 4}"{op}>{b}</text>'
            )
        return "\n".join(parts)

    def draw_read(y, fade_adapter=False, fade_pure=False):
        parts = []
        # Transcript body block
        parts.append(aligned_block(X0, y, BODY_BLK_W, BH, "transcript body", BODY_C))
        # Last body bases (always visible)
        for i, b in enumerate(body_tail):
            parts.append(base(bx_body + i*BW, y, b))
        # Genomic A-tract (always visible)
        parts.append(bases_row(bx_at, y, a_tract, POLYA_F, POLYA_B, POLYA_B))
        # Ambiguous T (always visible — highlighted)
        parts.append(
            f'<rect fill="{AMB_F}" height="{BH}" rx="{BR}" width="{BW}" '
            f'x="{bx_t}" y="{y}" stroke="{AMB_B}" stroke-width="1.2"/>'
            f'<text fill="{AMB_B}" font-size="11" font-weight="700" text-anchor="middle" '
            f'x="{bx_t + BW//2}" y="{y + BH//2 + 4}">T</text>'
        )
        # Pure poly(A) (may be faded)
        parts.append(bases_row(bx_pure, y, pure_a, POLYA_F, POLYA_B, POLYA_B, fade=fade_pure))
        # Adapter stub (may be faded)
        parts.append(bases_row(bx_adapt, y, adapter, ADAPTER_F, ADAPTER_B, ADAPTER_B, fade=fade_adapter))
        return "\n".join(parts)

    # ════════════════════════════════════════════════════════════════
    # Section 1: RAW DRS READ
    # ════════════════════════════════════════════════════════════════
    y_sec = 40
    L.append(section_head("RAW DRS READ (RNA 5\u2032 \u2192 3\u2032)", y_sec))

    y_r = y_sec + 20
    L.append(row_label("Read", y_r + BH//2 + 4, x=14))
    L.append(draw_read(y_r))

    # Braces below: poly(A) tail spans from first A-tract through pure A's
    L.append(brace_below(bx_at, bx_pure + n_pure * BW, y_r + BH + 1, POLYA_B,
                          "poly(A) tail (includes genomic A-tract)"))
    L.append(brace_below(bx_adapt, bx_adapt + n_adapt * BW, y_r + BH + 1, ADAPTER_B,
                          "adapter stub"))

    # Mark the ambiguous T
    L.append(f'<text fill="{AMB_B}" font-size="7" font-weight="600" '
             f'text-anchor="middle" x="{bx_t + BW//2}" y="{y_r - 4}">seq error or genomic?</text>')

    # 5'/3' labels
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" x="{X0 - 18}" '
             f'y="{y_r + BH//2 + 3}">5\u2032</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" x="{bx_adapt + n_adapt*BW + 6}" '
             f'y="{y_r + BH//2 + 3}">3\u2032</text>')

    # ════════════════════════════════════════════════════════════════
    # Section 2: TWO-PASS TRIMMING
    # ════════════════════════════════════════════════════════════════
    y_div = y_r + BH + 40
    L.append(hdivider(y_div))
    y_sec2 = y_div + 16
    L.append(section_head("TWO-PASS TRIMMING", y_sec2))

    # ── Pass 0: Adapter stub removal ──
    y_p0 = y_sec2 + 22
    L.append(f'<text fill="{ADAPTER_B}" font-size="9" font-weight="700" '
             f'x="14" y="{y_p0 + BH//2 + 4}">Pass 0</text>')
    L.append(draw_read(y_p0, fade_adapter=True))
    # Strikethrough on adapter
    L.append(f'<line stroke="{PAL["red"]}" stroke-width="2" '
             f'x1="{bx_adapt - 2}" x2="{bx_adapt + n_adapt*BW + 2}" '
             f'y1="{y_p0 + BH//2}" y2="{y_p0 + BH//2}" opacity="0.6"/>')
    L.append(f'<text fill="{ADAPTER_B}" font-size="8" font-style="italic" '
             f'x="{bx_adapt + n_adapt*BW + 8}" y="{y_p0 + BH//2 + 3}">'
             f'strip T[CT]{{0,10}}$</text>')

    # ── Pass 1: Pure-A scan (stop at first non-A) ──
    y_p1 = y_p0 + BH + 22
    L.append(f'<text fill="{POLYA_B}" font-size="9" font-weight="700" '
             f'x="14" y="{y_p1 + BH//2 + 4}">Pass 1</text>')
    L.append(draw_read(y_p1, fade_adapter=True, fade_pure=True))
    # Strikethrough on adapter (already gone) + pure poly(A)
    L.append(f'<line stroke="{PAL["red"]}" stroke-width="2" '
             f'x1="{bx_adapt - 2}" x2="{bx_adapt + n_adapt*BW + 2}" '
             f'y1="{y_p1 + BH//2}" y2="{y_p1 + BH//2}" opacity="0.3"/>')
    L.append(f'<line stroke="{PAL["red"]}" stroke-width="2" '
             f'x1="{bx_pure - 2}" x2="{bx_pure + n_pure*BW + 2}" '
             f'y1="{y_p1 + BH//2}" y2="{y_p1 + BH//2}" opacity="0.6"/>')
    # Scan arrow (leftward from 3' end toward T)
    scan_y = y_p1 - 6
    L.append(h_arrow(bx_pure + n_pure*BW - 4, bx_t + BW + 4, scan_y, POLYA_B,
                      "scan \u2190 (pure A only)"))
    # Stop marker at T
    L.append(f'<line stroke="{AMB_B}" stroke-width="2" '
             f'x1="{bx_t + BW + 2}" x2="{bx_t + BW + 2}" '
             f'y1="{y_p1 - 2}" y2="{y_p1 + BH + 2}"/>')
    L.append(f'<text fill="{AMB_B}" font-size="8" font-weight="600" '
             f'x="{bx_pure + n_pure*BW + 8}" y="{y_p1 + BH//2 + 3}">'
             f'stops at T \u2014 left for walk-back</text>')

    # ════════════════════════════════════════════════════════════════
    # Section 3: TRIMMED OUTPUT
    # ════════════════════════════════════════════════════════════════
    y_div2 = y_p1 + BH + 18
    L.append(hdivider(y_div2))
    y_sec3 = y_div2 + 16
    L.append(section_head("TRIMMED READ (for re-alignment)", y_sec3, PAL["green"]))

    y_out = y_sec3 + 20
    L.append(row_label("Output", y_out + BH//2 + 4, x=14))
    # Transcript body block
    L.append(aligned_block(X0, y_out, BODY_BLK_W, BH, "transcript body", BODY_C))
    # Body bases (kept)
    for i, b in enumerate(body_tail):
        L.append(base(bx_body + i*BW, y_out, b))
    # A-tract (kept — could be genomic)
    L.append(bases_row(bx_at, y_out, a_tract, POLYA_F, POLYA_B, POLYA_B))
    # Ambiguous T (kept)
    L.append(
        f'<rect fill="{AMB_F}" height="{BH}" rx="{BR}" width="{BW}" '
        f'x="{bx_t}" y="{y_out}" stroke="{AMB_B}" stroke-width="1.2"/>'
        f'<text fill="{AMB_B}" font-size="11" font-weight="700" text-anchor="middle" '
        f'x="{bx_t + BW//2}" y="{y_out + BH//2 + 4}">T</text>'
    )
    # Checkmark + note
    L.append(f'<text fill="{PAL["green"]}" font-size="10" font-weight="700" '
             f'x="{bx_t + BW + 10}" y="{y_out + BH//2 + 0}">'
             f'\u2713 ready for alignment</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'x="{bx_t + BW + 10}" y="{y_out + BH//2 + 12}">'
             f'ambiguous T resolved by 3\u2032 walk-back after alignment</text>')
    # Metadata note
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_out + BH + 22}">'
             f'tail length, adapter sequence, and pass number saved to metadata parquet</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: False Junction Walk-Back
# ═══════════════════════════════════════════════════════════════════════════════

def fig_false_junction_walkback():
    """
    False Junction Walk-Back on a pre-trimmed read.  When a pre-trimmed read
    is re-aligned, the aligner may introduce a spurious N-op (splice junction)
    to jump between two separated A-tracts, shifting the apparent 3' end
    downstream.  Walk-back absorbs the false N and recovers the true CPA.
    """
    H = 360
    X0 = 70
    GAP_W = 50
    L = []
    L.append(svg_open(H))
    L.append(fig_title("False Junction Walk-Back"))

    # -- BEFORE --
    y_sec = 42
    L.append(section_head("ALIGNED READ (pre-trimmed, poly(A) already removed)", y_sec, PAL["red"]))

    # Genome row: GTC | AAAA | --gap-- | AAAA | GTC
    y_g = y_sec + 28
    L.append(row_label("Genome", y_g + BH // 2 + 4))

    left_seq = list("GTC")
    for i, b in enumerate(left_seq):
        L.append(base(X0 + i * BW, y_g, b))

    a1_start = X0 + 3 * BW
    for i in range(4):
        L.append(base(a1_start + i * BW, y_g, "A", None, PAL["green"]))

    # Gap (genomic region between two A-tracts)
    gap_x1 = a1_start + 4 * BW
    gap_x2 = gap_x1 + GAP_W
    L.append(f'<line stroke="{PAL["muted"]}" stroke-dasharray="5,3" stroke-width="1.5" '
             f'x1="{gap_x1}" x2="{gap_x2}" y1="{y_g + BH // 2}" y2="{y_g + BH // 2}"/>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" text-anchor="middle" '
             f'x="{(gap_x1 + gap_x2) // 2}" y="{y_g + BH // 2 - 8}">gap</text>')

    a2_start = gap_x2
    for i in range(4):
        L.append(base(a2_start + i * BW, y_g, "A", None, PAL["green"]))

    ds_start = a2_start + 4 * BW
    for i, b in enumerate(list("GTC")):
        L.append(base(ds_start + i * BW, y_g, b))

    # Brackets
    L.append(brace_above(a1_start, a1_start + 4 * BW, y_g - 1, PAL["green"], "A-tract"))
    L.append(brace_above(a2_start, a2_start + 4 * BW, y_g - 1, PAL["green"], "A-tract"))

    # CPA marker
    L.append(vert_marker(a1_start, y_g - 14, y_g, PAL["green"], "true CPA", "left"))

    # Read row — pre-trimmed read with false N-op
    # The pre-trimmed read has ...GTCAAAAAT
    # Aligner jumped from first A-tract to second via false N-op to accommodate the T
    y_r = y_g + 50
    L.append(row_label("Read", y_r + BH // 2 + 4))

    for i, b in enumerate(left_seq):
        L.append(base(X0 + i * BW, y_r, b))

    # A's in first A-tract
    for i in range(4):
        L.append(base(a1_start + i * BW, y_r, "A"))

    # FALSE JUNCTION (N op) — aligner skipped the gap
    fj_x = gap_x1
    L.append(f'<rect fill="{PAL["red_l"]}" height="{BH}" rx="{BR}" '
             f'stroke="{PAL["red"]}" stroke-width="1.5" width="{GAP_W}" '
             f'x="{fj_x}" y="{y_r}"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="10" font-weight="700" '
             f'text-anchor="middle" x="{fj_x + GAP_W // 2}" y="{y_r + BH // 2 + 4}">N</text>')

    # A's in second A-tract
    for i in range(3):
        L.append(base(a2_start + i * BW, y_r, "A"))

    # Terminal T (ambiguous base from pre-trimmer)
    L.append(base(a2_start + 3 * BW, y_r, "T", PAL["orange_l"], PAL["orange"]))

    L.append(f'<text fill="{PAL["red"]}" font-size="9" font-weight="600" '
             f'text-anchor="middle" x="{fj_x + GAP_W // 2}" y="{y_r + BH + 14}">false N op</text>')

    # Apparent 3' end marker
    app_x = a2_start + 4 * BW
    L.append(vert_marker(app_x, y_r - 4, y_r + BH + 4, PAL["red"], "apparent 3\u2032 end", "right"))

    # -- DIVIDER --
    y_div = y_r + BH + 26
    L.append(hdivider(y_div))

    # -- AFTER --
    y_sec2 = y_div + 18
    L.append(section_head("AFTER WALK-BACK", y_sec2, PAL["green"]))

    y_arrow = y_sec2 + 18
    arrow_r = a2_start + 3 * BW
    arrow_l = a1_start + 4
    L.append(h_arrow(arrow_r, arrow_l, y_arrow, PAL["blue"],
                     "walk back: skip A\u2019s + T error, discard false N"))

    y_c = y_arrow + 24
    L.append(row_label("Corrected", y_c + BH // 2 + 4))
    for i, b in enumerate(left_seq):
        L.append(base(X0 + i * BW, y_c, b))

    L.append(vert_marker(a1_start, y_c - 8, y_c + BH + 4, PAL["green"], "true 3\u2032 end (CPA)"))

    # Note
    L.append(f'<text fill="{PAL["muted"]}" font-size="9" font-style="italic" '
             f'x="{a1_start + 8}" y="{y_c + BH + 22}">'
             f'false N-op discarded; CPA position recorded in corrected TSV</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    figures = {
        "indel_correction":       fig_indel_correction,
        "softclip_rescue":        fig_softclip_rescue,
        "5prime_junction_rescue": fig_5prime_junction,
        "multi_aligner_consensus": fig_multi_aligner_consensus,
        "false_junction_walkback": fig_false_junction_walkback,
        "polya_pretrim":          fig_polya_pretrim,
    }
    for name, fn in figures.items():
        svg_path = os.path.join(OUTDIR, f"{name}.svg")
        png_path = os.path.join(OUTDIR, f"{name}.png")
        svg_text = fn()
        with open(svg_path, "w") as f:
            f.write(svg_text)
        cairosvg.svg2png(bytestring=svg_text.encode(), write_to=png_path,
                         output_width=FIG_W * 2)
        print(f"  {name}: SVG + PNG ({FIG_W * 2}px wide)")
    print("Done - all figures generated.")
