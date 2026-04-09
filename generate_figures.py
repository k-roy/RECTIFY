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
    H = 370
    X0 = 70  # base start x
    L = []
    L.append(svg_open(H))
    L.append(fig_title("3\u2032 End Indel Correction"))

    # -- BEFORE --
    y_sec = 42
    L.append(section_head("BEFORE CORRECTION", y_sec, PAL["red"]))

    # Genome row — extra room for CPA marker above
    y_g = y_sec + 28
    L.append(row_label("Genome", y_g + BH // 2 + 4))
    genome = list("CAGTCAAAAAGAAAAGTC")
    for i, b in enumerate(genome):
        hf = None
        hb = None
        if 5 <= i <= 15 and b == "A":
            hb = PAL["green"]
        if i == 10 and b == "G":
            hb = PAL["red"]
            hf = PAL["red_l"]
        if i >= 16:
            hb = PAL["orange"]
        L.append(base(X0 + i * BW, y_g, b, hf, hb))

    # Brackets above genome
    L.append(brace_above(X0 + 5 * BW, X0 + 16 * BW, y_g - 1, PAL["green"], "genomic A-tract"))
    L.append(brace_above(X0 + 16 * BW, X0 + 18 * BW, y_g - 1, PAL["orange"], "non-A"))

    # CPA site marker — positioned below brackets, above genome row
    cpa_x = X0 + 5 * BW
    L.append(vert_marker(cpa_x, y_g - 14, y_g, PAL["green"], "CPA site", "left"))

    # Read row
    y_r = y_g + 50
    L.append(row_label("Read", y_r + BH // 2 + 4))
    transcript = list("CAGTC")
    for i, b in enumerate(transcript):
        L.append(base(X0 + i * BW, y_r, b))

    # Aligned A's with blue tint (poly(A) extending into A-tract)
    for i in range(4):
        L.append(base(X0 + (5 + i) * BW, y_r, "A", PAL["blue_l"], PAL["blue"]))

    # Deletion
    L.append(base(X0 + 9 * BW, y_r, "-", PAL["red_l"], PAL["red"]))

    # More aligned A's
    for i in range(4):
        L.append(base(X0 + (10 + i) * BW, y_r, "A", PAL["blue_l"], PAL["blue"]))

    # T error
    L.append(base(X0 + 14 * BW, y_r, "T", PAL["orange_l"], PAL["orange"]))

    # Another A
    L.append(base(X0 + 15 * BW, y_r, "A", PAL["blue_l"], PAL["blue"]))

    # Soft-clipped poly(A) tail
    for i in range(6):
        L.append(base(X0 + (16 + i) * BW, y_r, "A", "#bfdbfe", "#1d4ed8"))

    # Labels
    L.append(brace_below(X0 + 5 * BW, X0 + 16 * BW, y_r + BH + 1, PAL["blue"],
                          "poly(A) extended into A-tract (with deletion + T error)"))
    L.append(brace_above(X0 + 16 * BW, X0 + 22 * BW, y_r - 1, "#1d4ed8", "soft-clipped tail"))

    # -- DIVIDER --
    y_div = y_r + BH + 30
    L.append(hdivider(y_div))

    # -- AFTER --
    y_sec2 = y_div + 18
    L.append(section_head("AFTER CORRECTION", y_sec2, PAL["green"]))

    # Walk-back arrow
    y_arrow = y_sec2 + 18
    arrow_r = X0 + 21 * BW
    arrow_l = cpa_x + 4
    L.append(h_arrow(arrow_r, arrow_l, y_arrow, PAL["blue"],
                     "walk back: skip A\u2019s, deletion, T error \u2192 stop at C"))

    # Corrected row
    y_c = y_arrow + 24
    L.append(row_label("Corrected", y_c + BH // 2 + 4))
    for i, b in enumerate(list("CAGTC")):
        L.append(base(X0 + i * BW, y_c, b))

    # Green marker at true 3' end
    L.append(vert_marker(cpa_x, y_c - 8, y_c + BH + 4, PAL["green"], "true 3\u2032 end"))

    # Poly(A) tail -- clean
    for i in range(17):
        L.append(base(cpa_x + i * BW, y_c, "A", PAL["blue_l"], PAL["blue"]))
    L.append(brace_below(cpa_x, cpa_x + 17 * BW, y_c + BH + 1, PAL["blue"],
                          "poly(A) tail (all A\u2019s + deletions + T errors + soft-clipped A\u2019s)"))

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
    H = 310
    BLK_H = 30
    L = []
    L.append(svg_open(H))
    L.append(fig_title("5\u2032 End Correction: Junction Soft-Clips"))

    # -- BEFORE --
    y_sec = 42
    L.append(section_head("BEFORE CORRECTION", y_sec, PAL["red"]))

    # Genome
    y_g = y_sec + 22
    L.append(row_label("Genome", y_g + BLK_H // 2 + 4))
    L.append(exon(90, y_g, 180, BLK_H, "Exon 1"))
    L.append(intron(270, 390, y_g + BLK_H // 2))
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" x="273" y="{y_g + BLK_H - 1}">GT</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" text-anchor="end" x="387" y="{y_g + BLK_H - 1}">AG</text>')
    L.append(exon(390, y_g, 180, BLK_H, "Exon 2"))

    # Read (with soft-clip)
    y_r = y_g + 48
    L.append(row_label("Read", y_r + BLK_H // 2 + 4))
    L.append(softclip_block(220, y_r, 170, BLK_H, "soft-clipped bases"))
    L.append(aligned_block(390, y_r, 180, BLK_H, "aligned to Exon 2"))

    # Annotation
    sc_mid = 220 + 85
    L.append(f'<text fill="{PAL["red"]}" font-size="9" text-anchor="middle" '
             f'x="{sc_mid}" y="{y_r + BLK_H + 14}">these bases match Exon 1</text>')
    arr_y = y_r + BLK_H + 20
    L.append(f'<line stroke="{PAL["red"]}" stroke-width="1.2" x1="228" x2="382" y1="{arr_y}" y2="{arr_y}"/>')
    L.append(f'<polygon fill="{PAL["red"]}" points="228,{arr_y} 235,{arr_y-3} 235,{arr_y+3}"/>')
    L.append(f'<polygon fill="{PAL["red"]}" points="382,{arr_y} 375,{arr_y-3} 375,{arr_y+3}"/>')

    # -- DIVIDER --
    y_div = arr_y + 18
    L.append(hdivider(y_div))

    # -- AFTER --
    y_sec2 = y_div + 18
    L.append(section_head("AFTER CORRECTION", y_sec2, PAL["green"]))

    # Genome 2
    y_g2 = y_sec2 + 22
    L.append(row_label("Genome", y_g2 + BLK_H // 2 + 4))
    L.append(exon(90, y_g2, 180, BLK_H, "Exon 1"))
    L.append(intron(270, 390, y_g2 + BLK_H // 2))
    L.append(exon(390, y_g2, 180, BLK_H, "Exon 2"))

    # Corrected read
    y_c = y_g2 + 44
    L.append(row_label("Corrected", y_c + BLK_H // 2 + 4))
    L.append(rescued_block(100, y_c, 170, BLK_H, "rescued \u2192 Exon 1"))
    L.append(intron(270, 390, y_c + BLK_H // 2, "N (splice)"))
    L.append(aligned_block(390, y_c, 180, BLK_H, "Exon 2"))

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 4: Multi-Aligner Consensus
# ═══════════════════════════════════════════════════════════════════════════════

def fig_multi_aligner_consensus():
    """
    Same read, three aligners, three different structural decisions:
      - mapPacBio: clean alignment to annotated junction (wins)
      - minimap2: soft-clips the entire 5' end of the read (no Exon 1)
      - gapmm2: forces annotated junction but with junction-proximal indels
      - mapPacBio: discovers a novel alternative junction (diff 5'SS & 3'SS)
    """
    H = 530
    BLK_H = 24
    ROW_SP = 38
    L = []
    L.append(svg_open(H))
    L.append(fig_title("Multi-Aligner Consensus: Selecting the Optimal Alignment"))

    # Genome reference coordinates
    # Exon 1: 120-280, Intron: 280-390, Exon 2: 390-540
    EX1_X, EX1_W = 120, 160
    EX2_X, EX2_W = 390, 150
    INT_X1, INT_X2 = 280, 390

    # Aligner colors
    AC = {"minimap2": "#3b82f6", "gapmm2": "#22c55e", "mapPacBio": "#f59e0b"}

    # Junction shift for mapPacBio's novel junction
    NOVEL_5SS = 10
    NOVEL_3SS = 8

    # ── SAME READ, THREE DIFFERENT STRUCTURAL DECISIONS ──
    y_sec = 40
    L.append(section_head("SAME READ, THREE ALIGNERS \u2192 THREE DIFFERENT STRUCTURES", y_sec))

    # Genome annotation
    y_g = y_sec + 24
    L.append(row_label("Genome", y_g + BLK_H // 2 + 4))
    L.append(exon(EX1_X, y_g, EX1_W, BLK_H, "Exon 1"))
    L.append(intron(INT_X1, INT_X2, y_g + BLK_H // 2))
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" x="{INT_X1+4}" y="{y_g + BLK_H - 2}">GT</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" text-anchor="end" x="{INT_X2-4}" y="{y_g + BLK_H - 2}">AG</text>')
    L.append(exon(EX2_X, y_g, EX2_W, BLK_H, "Exon 2"))

    # ── minimap2: Soft-clips the entire 5' end of the read ──
    # minimap2 fails to split the read across the junction for this read.
    # It aligns Exon 2 cleanly but soft-clips everything upstream (Exon 1 bases).
    # The soft-clip is at the 5' end of the read — a valid terminal position.
    y = y_g + 40
    c = AC["minimap2"]
    L.append(f'<text fill="{c}" font-size="10" font-weight="600" x="14" y="{y + BLK_H//2 + 4}">minimap2</text>')
    # Soft-clipped Exon 1 region (entire 5' end of read)
    L.append(f'<rect fill="{PAL["red_l"]}" height="{BLK_H}" rx="4" '
             f'stroke="{PAL["red"]}" stroke-dasharray="4,2" stroke-width="0.8" opacity="0.6" '
             f'width="{EX1_W}" x="{EX1_X}" y="{y}"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="8" font-style="italic" text-anchor="middle" '
             f'x="{EX1_X + EX1_W//2}" y="{y + BLK_H//2 + 3}">soft-clipped</text>')
    # No junction — minimap2 didn't find the splice
    # Exon 2 aligned directly
    L.append(aligned_block(EX2_X, y, EX2_W, BLK_H, "Exon 2", c))
    L.append(f'<text fill="{PAL["red"]}" font-size="8" x="{EX2_X + EX2_W + 8}" '
             f'y="{y + BLK_H//2 + 0}">\u2717 no junction found</text>')
    L.append(f'<text fill="{PAL["red"]}" font-size="8" x="{EX2_X + EX2_W + 8}" '
             f'y="{y + BLK_H//2 + 12}">\u2717 5\u2032 end entirely soft-clipped</text>')

    # ── gapmm2: Annotated junction, but junction-proximal indels ──
    # gapmm2 forces the read onto the annotated splice site coordinates,
    # but introduces small insertions/deletions near the junction to make it fit.
    y += ROW_SP + 4
    c = AC["gapmm2"]
    L.append(f'<text fill="{c}" font-size="10" font-weight="600" x="14" y="{y + BLK_H//2 + 4}">gapmm2</text>')
    L.append(aligned_block(EX1_X, y, EX1_W, BLK_H, "Exon 1", c))
    L.append(f'<line stroke="{c}" stroke-dasharray="5,3" stroke-width="1" '
             f'x1="{INT_X1}" x2="{INT_X2}" y1="{y+BLK_H//2}" y2="{y+BLK_H//2}"/>')
    L.append(aligned_block(EX2_X, y, EX2_W, BLK_H, "Exon 2", c))
    # Indel markers — small triangles below the junction boundaries
    tri_y = y + BLK_H
    ins_x = INT_X1 - 4
    L.append(f'<polygon fill="{PAL["red"]}" opacity="0.8" '
             f'points="{ins_x},{tri_y} {ins_x+4},{tri_y+6} {ins_x+8},{tri_y}"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="6" text-anchor="middle" '
             f'x="{ins_x+4}" y="{tri_y+13}">ins</text>')
    del_x = EX2_X + 2
    L.append(f'<polygon fill="{PAL["red"]}" opacity="0.8" '
             f'points="{del_x},{tri_y} {del_x+4},{tri_y+6} {del_x+8},{tri_y}"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="6" text-anchor="middle" '
             f'x="{del_x+4}" y="{tri_y+13}">del</text>')
    L.append(f'<text fill="{PAL["orange"]}" font-size="8" x="{EX2_X + EX2_W + 8}" '
             f'y="{y + BLK_H//2 + 0}">forces annotated GT-AG junction</text>')
    L.append(f'<text fill="{PAL["red"]}" font-size="8" x="{EX2_X + EX2_W + 8}" '
             f'y="{y + BLK_H//2 + 12}">\u2717 junction-proximal indels</text>')

    # ── mapPacBio: Novel alternative junction (different 5'SS and 3'SS) ──
    y += ROW_SP + 18  # extra space for gapmm2's indel markers
    c = AC["mapPacBio"]
    L.append(f'<text fill="{c}" font-size="10" font-weight="600" x="14" y="{y + BLK_H//2 + 4}">mapPacBio</text>')
    # Exon 1 extends further (5'SS shifted into intron)
    mpb_ex1_w = EX1_W + NOVEL_5SS
    L.append(aligned_block(EX1_X, y, mpb_ex1_w, BLK_H, "Exon 1", c))
    # Junction from shifted positions
    mpb_int_x1 = INT_X1 + NOVEL_5SS
    mpb_int_x2 = INT_X2 + NOVEL_3SS
    L.append(f'<line stroke="{c}" stroke-dasharray="5,3" stroke-width="1" '
             f'x1="{mpb_int_x1}" x2="{mpb_int_x2}" y1="{y+BLK_H//2}" y2="{y+BLK_H//2}"/>')
    mpb_ex2_x = EX2_X + NOVEL_3SS
    mpb_ex2_w = EX2_W - NOVEL_3SS
    L.append(aligned_block(mpb_ex2_x, y, mpb_ex2_w, BLK_H, "Exon 2", c))
    L.append(f'<text fill="{PAL["green"]}" font-size="8" x="{EX2_X + EX2_W + 8}" '
             f'y="{y + BLK_H//2 + 0}">\u2713 novel GC-AG junction (correct)</text>')
    L.append(f'<text fill="{PAL["green"]}" font-size="8" x="{EX2_X + EX2_W + 8}" '
             f'y="{y + BLK_H//2 + 12}">\u2713 clean alignment, no indels</text>')

    # ── Callout ──
    y_note = y + BLK_H + 14
    L.append(f'<text fill="{PAL["heading"]}" font-size="9" font-style="italic" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_note}">'
             f'Each aligner makes different structural decisions \u2014 '
             f'scoring selects the best overall alignment per read</text>')

    # ── DIVIDER ──
    y_div = y_note + 12
    L.append(hdivider(y_div))

    # ── SCORE AND SELECT ──
    y_sec2 = y_div + 16
    L.append(section_head("SCORE EACH ALIGNMENT \u2192 SELECT BEST PER READ", y_sec2))

    # Column headers
    y_hdr = y_sec2 + 20
    score_cols = [("junction", 555), ("annotated", 610), ("indels", 662), ("5\u2032 cov.", 712), ("score", 748)]
    for label, cx in score_cols:
        L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-weight="600" '
                 f'text-anchor="middle" x="{cx}" y="{y_hdr}">{label}</text>')

    # Score rows — gapmm2 wins (annotated junction trumps indel penalty),
    # minimap2 middle (soft-clip rescuable), mapPacBio lowest (novel junction
    # penalized heavily even though it's biologically correct for this read)
    scores = [
        ("gapmm2",     AC["gapmm2"],     "GT-AG", "\u2713", "2 prox.", "full",  "4", True),
        ("minimap2",   AC["minimap2"],   "\u2014",  "\u2014",  "none",   "none",  "2", False),
        ("mapPacBio",  AC["mapPacBio"],  "GC-AG", "\u2717", "none",    "full",  "1", False),
    ]

    for i, (name, color, jn, ann, indels, cov, score, is_best) in enumerate(scores):
        y = y_hdr + 8 + i * 30
        L.append(f'<text fill="{color}" font-size="10" font-weight="600" x="14" '
                 f'y="{y + BLK_H//2 + 4}">{name}</text>')

        # Compact alignment bars
        bar_x = 120
        if name == "minimap2":
            # Soft-clipped Exon 1 + Exon 2 only
            L.append(f'<rect fill="{PAL["red_l"]}" height="{BLK_H}" rx="3" '
                     f'stroke="{PAL["red"]}" stroke-dasharray="3,2" stroke-width="0.6" opacity="0.4" '
                     f'width="{EX1_W}" x="{bar_x}" y="{y}"/>')
            L.append(aligned_block(EX2_X, y, EX2_W, BLK_H, "Exon 2", color))
        elif name == "gapmm2":
            L.append(aligned_block(bar_x, y, EX1_W, BLK_H, "Exon 1", color))
            L.append(f'<line stroke="{color}" stroke-dasharray="5,3" stroke-width="1" '
                     f'x1="{INT_X1}" x2="{INT_X2}" y1="{y+BLK_H//2}" y2="{y+BLK_H//2}"/>')
            L.append(aligned_block(EX2_X, y, EX2_W, BLK_H, "Exon 2", color))
            # Tiny indel marks
            L.append(f'<polygon fill="{PAL["red"]}" opacity="0.6" '
                     f'points="{INT_X1-2},{y+BLK_H} {INT_X1+1},{y+BLK_H+3} {INT_X1+4},{y+BLK_H}"/>')
            L.append(f'<polygon fill="{PAL["red"]}" opacity="0.6" '
                     f'points="{EX2_X+1},{y+BLK_H} {EX2_X+4},{y+BLK_H+3} {EX2_X+7},{y+BLK_H}"/>')
        elif name == "mapPacBio":
            L.append(aligned_block(bar_x, y, EX1_W + NOVEL_5SS, BLK_H, "Exon 1", color))
            L.append(f'<line stroke="{color}" stroke-dasharray="5,3" stroke-width="1" '
                     f'x1="{INT_X1 + NOVEL_5SS}" x2="{INT_X2 + NOVEL_3SS}" '
                     f'y1="{y+BLK_H//2}" y2="{y+BLK_H//2}"/>')
            L.append(aligned_block(EX2_X + NOVEL_3SS, y, EX2_W - NOVEL_3SS, BLK_H, "Exon 2", color))

        # Score values
        jn_color = PAL["green"] if jn == "GT-AG" else (PAL["orange"] if jn == "GC-AG" else PAL["muted"])
        L.append(f'<text fill="{jn_color}" font-size="9" text-anchor="middle" '
                 f'x="555" y="{y + BLK_H//2 + 4}">{jn}</text>')
        ann_color = PAL["green"] if ann == "\u2713" else (PAL["red"] if ann == "\u2717" else PAL["muted"])
        L.append(f'<text fill="{ann_color}" font-size="10" font-weight="600" text-anchor="middle" '
                 f'x="610" y="{y + BLK_H//2 + 4}">{ann}</text>')
        indel_color = PAL["green"] if indels == "none" else PAL["red"]
        L.append(f'<text fill="{indel_color}" font-size="9" font-weight="600" text-anchor="middle" '
                 f'x="662" y="{y + BLK_H//2 + 4}">{indels}</text>')
        cov_color = PAL["green"] if cov == "full" else PAL["red"]
        L.append(f'<text fill="{cov_color}" font-size="9" font-weight="600" text-anchor="middle" '
                 f'x="712" y="{y + BLK_H//2 + 4}">{cov}</text>')

        # Total score
        sc_weight = "700" if is_best else "600"
        sc_color = PAL["green"] if is_best else PAL["heading"]
        L.append(f'<text fill="{sc_color}" font-size="12" font-weight="{sc_weight}" '
                 f'text-anchor="middle" x="748" y="{y + BLK_H//2 + 4}">{score}</text>')

        if is_best:
            L.append(f'<text fill="{PAL["green"]}" font-size="9" font-weight="700" '
                     f'x="762" y="{y + BLK_H//2 + 4}">\u2190 best</text>')

    # Final output label
    y_out = y_hdr + 8 + 3 * 30 + 4
    L.append(hdivider(y_out))
    L.append(f'<text fill="{PAL["heading"]}" font-size="10" font-weight="600" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_out + 18}">'
             f'Best alignment per read \u2192 consensus BAM</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 5: False Junction Walk-Back
# ═══════════════════════════════════════════════════════════════════════════════

def fig_false_junction_walkback():
    H = 360
    X0 = 70
    GAP_W = 50
    L = []
    L.append(svg_open(H))
    L.append(fig_title("False Junction Walk-Back"))

    # -- BEFORE --
    y_sec = 42
    L.append(section_head("BEFORE CORRECTION", y_sec, PAL["red"]))

    # Genome row — extra room for CPA marker + brackets above
    y_g = y_sec + 28
    L.append(row_label("Genome", y_g + BH // 2 + 4))

    left_seq = list("GTC")
    for i, b in enumerate(left_seq):
        L.append(base(X0 + i * BW, y_g, b))

    a1_start = X0 + 3 * BW
    for i in range(4):
        L.append(base(a1_start + i * BW, y_g, "A", None, PAL["green"]))

    # Gap
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

    # CPA marker — positioned below brackets
    L.append(vert_marker(a1_start, y_g - 14, y_g, PAL["green"], "true CPA", "left"))

    # Read row
    y_r = y_g + 50
    L.append(row_label("Read", y_r + BH // 2 + 4))

    for i, b in enumerate(left_seq):
        L.append(base(X0 + i * BW, y_r, b))

    for i in range(4):
        L.append(base(a1_start + i * BW, y_r, "A", PAL["blue_l"], PAL["blue"]))

    # FALSE JUNCTION -- red N box
    fj_x = gap_x1
    L.append(f'<rect fill="{PAL["red_l"]}" height="{BH}" rx="{BR}" '
             f'stroke="{PAL["red"]}" stroke-width="1.5" width="{GAP_W}" '
             f'x="{fj_x}" y="{y_r}"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="10" font-weight="700" '
             f'text-anchor="middle" x="{fj_x + GAP_W // 2}" y="{y_r + BH // 2 + 4}">N</text>')

    for i in range(4):
        L.append(base(a2_start + i * BW, y_r, "A", PAL["blue_l"], PAL["blue"]))

    sc_start = a2_start + 4 * BW
    for i in range(3):
        L.append(base(sc_start + i * BW, y_r, "A", "#bfdbfe", "#1d4ed8"))
    L.append(brace_above(sc_start, sc_start + 3 * BW, y_r - 1, "#1d4ed8", "soft-clip"))

    L.append(f'<text fill="{PAL["red"]}" font-size="9" font-weight="600" '
             f'text-anchor="middle" x="{fj_x + GAP_W // 2}" y="{y_r + BH + 14}">false N op</text>')

    # -- DIVIDER --
    y_div = y_r + BH + 26
    L.append(hdivider(y_div))

    # -- AFTER --
    y_sec2 = y_div + 18
    L.append(section_head("AFTER CORRECTION", y_sec2, PAL["green"]))

    y_arrow = y_sec2 + 18
    arrow_r = sc_start + 2 * BW
    arrow_l = a1_start + 4
    L.append(h_arrow(arrow_r, arrow_l, y_arrow, PAL["blue"],
                     "walk back: eat A\u2019s, discard N"))

    y_c = y_arrow + 24
    L.append(row_label("Corrected", y_c + BH // 2 + 4))
    for i, b in enumerate(left_seq):
        L.append(base(X0 + i * BW, y_c, b))

    L.append(vert_marker(a1_start, y_c - 8, y_c + BH + 4, PAL["green"], "correct 3\u2032 end"))

    tail_n = 4 + 4 + 3
    for i in range(tail_n):
        L.append(base(a1_start + i * BW, y_c, "A", PAL["blue_l"], PAL["blue"]))
    L.append(brace_below(a1_start, a1_start + tail_n * BW, y_c + BH + 1, PAL["blue"],
                          "poly(A) tail (false N discarded)"))

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
