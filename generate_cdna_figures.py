#!/usr/bin/env python3
"""
Generate SVG + PNG figures for RECTIFY cDNA documentation.

Figures generated:
  1. cdna_pipeline_overview — 3-stage workflow (correct-cdna → align → cdna-analyze)
  2. cdna_umi_architecture  — PCB114.24 read structure, pipeline, subtypes
  3. cdna_poa_consensus     — UMI consensus pipeline flow with stats
  4. cdna_isoform_clustering — isoform grouping by (pos5, pos3) + T1/T2 reconciliation
  5. walkback_readvsref     — read-vs-reference walkback algorithm (correct algorithm + bug note)
  6. splice_classification  — GT-AG / AT-AC / non-canonical + chimeric voting

Run from /Users/kevinroy/work/rectify/:
    python3 generate_cdna_figures.py
"""

import os
import pathlib
import cairosvg

OUTDIR = pathlib.Path("/Users/kevinroy/work/rectify/docs/figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════
# Design tokens (matches existing figures exactly)
# ═══════════════════════════════════════════════════════════════════════

FONT  = "Helvetica Neue, Helvetica, Arial, sans-serif"
FIG_W = 760

PAL = dict(
    title    = "#1e293b",
    heading  = "#475569",
    label    = "#64748b",
    muted    = "#94a3b8",
    border   = "#cbd5e1",
    divider  = "#e2e8f0",
    bg       = "#ffffff",
    blue     = "#2563eb",
    blue_l   = "#dbeafe",
    green    = "#059669",
    green_l  = "#d1fae5",
    red      = "#dc2626",
    red_l    = "#fee2e2",
    teal     = "#0d9488",
    teal_l   = "#ccfbf1",
    orange   = "#d97706",
    orange_l = "#fef3c7",
    exon     = "#6366f1",
    exon_t   = "#ffffff",
    exon_l   = "#e0e7ff",
    A_f = "#dcfce7",  A_t = "#166534",
    T_f = "#fef3c7",  T_t = "#92400e",
    G_f = "#d1fae5",  G_t = "#065f46",
    C_f = "#dbeafe",  C_t = "#1e40af",
)

BW = 26   # nucleotide box width
BH = 24   # nucleotide box height
BR = 3    # border radius

BASE_FILL = {"A": PAL["A_f"], "T": PAL["T_f"], "G": PAL["G_f"], "C": PAL["C_f"],
             "-": PAL["red_l"], "N": PAL["red_l"]}
BASE_TEXT = {"A": PAL["A_t"], "T": PAL["T_t"], "G": PAL["G_t"], "C": PAL["C_t"],
             "-": PAL["red"],  "N": PAL["red"]}

# ═══════════════════════════════════════════════════════════════════════
# Shared drawing helpers
# ═══════════════════════════════════════════════════════════════════════

def svg_open(h):
    return (
        f'<?xml version="1.0" encoding="utf-8" ?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" height="{h}" '
        f'viewBox="0 0 {FIG_W} {h}">\n'
        f'<rect fill="{PAL["bg"]}" width="{FIG_W}" height="{h}"/>\n'
        f'<g font-family="{FONT}">'
    )

def svg_close():
    return "\n</g>\n</svg>"

def fig_title(text, y=22):
    return (f'<text fill="{PAL["title"]}" font-size="14" font-weight="bold" '
            f'text-anchor="middle" x="{FIG_W // 2}" y="{y}">{text}</text>')

def section_head(text, y, color=None):
    c = color or PAL["heading"]
    return (f'<text fill="{c}" font-size="10" font-weight="bold" '
            f'letter-spacing="0.8" x="60" y="{y}">{text}</text>')

def hdivider(y, x1=50, x2=710):
    return f'<line stroke="{PAL["divider"]}" stroke-width="1" x1="{x1}" x2="{x2}" y1="{y}" y2="{y}"/>'

def base_box(x, y, b, hi_fill=None, hi_border=None):
    """Single nucleotide box."""
    fill = hi_fill or BASE_FILL.get(b, "#f1f5f9")
    tc   = BASE_TEXT.get(b, PAL["label"])
    bdr  = hi_border or PAL["border"]
    sw   = "1.2" if hi_border else "0.7"
    return (
        f'<rect fill="{fill}" height="{BH}" rx="{BR}" width="{BW}" '
        f'x="{x}" y="{y}" stroke="{bdr}" stroke-width="{sw}"/>'
        f'<text fill="{tc}" font-size="11" font-weight="600" text-anchor="middle" '
        f'x="{x + BW // 2}" y="{y + BH // 2 + 4}">{b}</text>'
    )

def filled_block(x, y, w, h, label, fill, tc="#ffffff", rx=4):
    mid_x = x + w / 2
    mid_y = y + h / 2 + 4
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"/>'
        f'<text fill="{tc}" font-size="10" font-weight="600" text-anchor="middle" '
        f'x="{mid_x}" y="{mid_y}">{label}</text>'
    )

def outline_block(x, y, w, h, label, color, opacity=0.18, rx=4):
    """Semi-transparent outlined block (like cDNA body)."""
    mid_x = x + w / 2
    mid_y = y + h / 2 + 4
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="{color}" opacity="{opacity}"/>'
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="none" stroke="{color}" stroke-width="1"/>'
        f'<text fill="{color}" font-size="10" font-weight="600" text-anchor="middle" '
        f'x="{mid_x}" y="{mid_y}">{label}</text>'
    )

def flow_box(x, y, w, h, label, color, opacity=0.15, rx=5):
    """Pipeline flow box with border."""
    mid_x = x + w / 2
    mid_y = y + h / 2 + 4
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="{color}" opacity="{opacity}"/>'
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="none" stroke="{color}" stroke-width="1.2"/>'
        f'<text fill="{color}" font-size="9" font-weight="600" text-anchor="middle" '
        f'x="{mid_x}" y="{mid_y}">{label}</text>'
    )

def h_arrow(x1, x2, y, color, label=None):
    """Horizontal arrow with optional label above."""
    parts = [f'<line stroke="{color}" stroke-width="1.5" x1="{x1}" x2="{x2}" y1="{y}" y2="{y}"/>']
    if x2 > x1:
        parts.append(f'<polygon fill="{color}" points="{x2},{y} {x2-7},{y-3.5} {x2-7},{y+3.5}"/>')
    else:
        parts.append(f'<polygon fill="{color}" points="{x2},{y} {x2+7},{y-3.5} {x2+7},{y+3.5}"/>')
    if label:
        mx = (x1 + x2) / 2
        parts.append(f'<text fill="{color}" font-size="8.5" text-anchor="middle" x="{mx}" y="{y - 6}">{label}</text>')
    return "\n".join(parts)

def v_arrow(x, y1, y2, color):
    """Vertical downward arrow."""
    return (
        f'<line stroke="{color}" stroke-width="1.3" x1="{x}" x2="{x}" y1="{y1}" y2="{y2}"/>\n'
        f'<polygon fill="{color}" points="{x},{y2} {x-3.5},{y2-6} {x+3.5},{y2-6}"/>'
    )

def vert_marker(x, y1, y2, color, label=None, side="right"):
    """Vertical line marker."""
    parts = [f'<line stroke="{color}" stroke-width="2" x1="{x}" x2="{x}" y1="{y1}" y2="{y2}"/>']
    if label:
        if side == "right":
            parts.append(f'<text fill="{color}" font-size="9" font-weight="600" x="{x+4}" y="{y1-3}">{label}</text>')
        else:
            parts.append(f'<text fill="{color}" font-size="9" font-weight="600" text-anchor="end" x="{x-4}" y="{y1-3}">{label}</text>')
    return "\n".join(parts)

def brace_below(x1, x2, y, color, label):
    mid = (x1 + x2) / 2
    return (
        f'<path d="M{x1},{y} L{x1},{y+5} L{x2},{y+5} L{x2},{y}" '
        f'fill="none" stroke="{color}" stroke-width="1" opacity="0.6"/>\n'
        f'<text fill="{color}" font-size="8" text-anchor="middle" x="{mid}" y="{y+17}">{label}</text>'
    )

def dashed_line(x1, x2, y, color, sw=1.3):
    return (f'<line stroke="{color}" stroke-dasharray="5,3" stroke-width="{sw}" '
            f'x1="{x1}" x2="{x2}" y1="{y}" y2="{y}"/>')

def nt_row(bases, x0, y, highlights=None):
    """Draw a row of nucleotide boxes; highlights = {index: (fill, border)}."""
    parts = []
    for i, b in enumerate(bases):
        hf, hbdr = None, None
        if highlights and i in highlights:
            hf, hbdr = highlights[i]
        parts.append(base_box(x0 + i * BW, y, b, hf, hbdr))
    return "\n".join(parts)


# ═══════════════════════════════════════════════════════════════════════
# Figure 1: cdna_umi_architecture (760 × 420)
# ═══════════════════════════════════════════════════════════════════════

def fig_cdna_umi_architecture():
    H = 420
    L = []
    L.append(svg_open(H))
    L.append(fig_title("PCB114.24 ONT cDNA Read Architecture"))

    # ── Section 1: Read structure ──────────────────────────────────────
    y_sec = 42
    L.append(section_head("READ STRUCTURE  (5′ → 3′ basecall orientation)", y_sec))

    y_rd = y_sec + 18
    BH_RD = 28  # slightly taller blocks for this section

    # 5' label
    L.append(f'<text fill="{PAL["muted"]}" font-size="9" x="46" y="{y_rd + BH_RD//2 + 4}">5′</text>')

    # SSP (purple/exon)
    L.append(filled_block(58, y_rd, 80, BH_RD, "SSP", PAL["exon"]))
    # UMI 27nt (blue)
    L.append(filled_block(138, y_rd, 100, BH_RD, "UMI 27nt", PAL["blue"]))
    # GGG bridge (teal)
    L.append(filled_block(238, y_rd, 44, BH_RD, "GGG", PAL["teal"]))
    # cDNA body (outlined indigo)
    L.append(outline_block(282, y_rd, 258, BH_RD, "cDNA body", PAL["exon"]))
    # poly(A) (outlined green)
    L.append(outline_block(540, y_rd, 80, BH_RD, "poly(A)", PAL["green"]))
    # adapter (muted grey)
    L.append(outline_block(620, y_rd, 80, BH_RD, "adapter", PAL["muted"], opacity=0.25))

    # 3' label
    L.append(f'<text fill="{PAL["muted"]}" font-size="9" x="706" y="{y_rd + BH_RD//2 + 4}">3′</text>')

    # Tag annotations below blocks
    y_ann = y_rd + BH_RD + 6
    # XU tag under UMI
    L.append(f'<line x1="188" x2="188" y1="{y_ann}" y2="{y_ann+10}" stroke="{PAL["blue"]}" stroke-width="1" stroke-dasharray="3,2"/>')
    L.append(f'<text fill="{PAL["blue"]}" font-size="8" font-weight="600" text-anchor="middle" x="188" y="{y_ann+21}">XU:Z → UMI</text>')
    # XQ tag under SSP/UMI/GGG region
    L.append(f'<text fill="{PAL["muted"]}" font-size="7.5" text-anchor="middle" x="150" y="{y_ann+34}">(TT-VVVV)×4+TTT</text>')
    # XA tag under poly(A)
    L.append(f'<line x1="580" x2="580" y1="{y_ann}" y2="{y_ann+10}" stroke="{PAL["green"]}" stroke-width="1" stroke-dasharray="3,2"/>')
    L.append(f'<text fill="{PAL["green"]}" font-size="8" font-weight="600" text-anchor="middle" x="580" y="{y_ann+21}">XA:i → poly(A)</text>')
    # XQ strip indicator
    L.append(f'<text fill="{PAL["teal"]}" font-size="7.5" text-anchor="middle" x="282" y="{y_ann+34}">XQ:i (5′ stripped)</text>')
    # XK strip indicator
    L.append(f'<text fill="{PAL["teal"]}" font-size="7.5" text-anchor="middle" x="540" y="{y_ann+34}">XK:i (3′ stripped)</text>')

    # Divider
    y_div1 = y_ann + 44
    L.append(hdivider(y_div1))

    # ── Section 2: After pretrim ──────────────────────────────────────
    y_sec2 = y_div1 + 14
    L.append(section_head("AFTER PRETRIM  (stage1_consensus.fastq.gz)", y_sec2, PAL["green"]))

    y_pt = y_sec2 + 16
    BH_PT = 24
    # Only cDNA body remains
    L.append(filled_block(200, y_pt, 280, BH_PT, "cDNA body  (clean mRNA sequence)", PAL["exon"]))
    L.append(f'<text fill="{PAL["label"]}" font-size="8" font-style="italic" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_pt + BH_PT + 14}">'
             f'pretrim_consensus() strips SSP+UMI+GGG (5′) and poly(A) (3′) before alignment</text>')

    # Divider
    y_div2 = y_pt + BH_PT + 28
    L.append(hdivider(y_div2))

    # ── Section 3: Read subtypes ──────────────────────────────────────
    y_sec3 = y_div2 + 14
    L.append(section_head("READ SUBTYPES", y_sec3))

    BH_ST = 22
    COL_W = 220
    GAP   = 12
    COL_X = [40, 40 + COL_W + GAP, 40 + 2*(COL_W + GAP)]
    y_st  = y_sec3 + 16

    # umi_captured_fwd: SSP at 5'
    x = COL_X[0]
    L.append(filled_block(x, y_st, 44, BH_ST, "SSP", PAL["exon"]))
    L.append(filled_block(x+44, y_st, 50, BH_ST, "UMI", PAL["blue"]))
    L.append(outline_block(x+94, y_st, 106, BH_ST, "cDNA", PAL["exon"]))
    L.append(f'<text fill="{PAL["heading"]}" font-size="9" font-weight="bold" x="{x}" y="{y_st + BH_ST + 14}">UMI captured (fwd)</text>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8" x="{x}" y="{y_st + BH_ST + 26}">SSP+UMI at basecalled-5′  ·  XT:i:1</text>')
    L.append(f'<text fill="{PAL["blue"]}" font-size="8" x="{x}" y="{y_st + BH_ST + 37}">XY:Z=umi_captured_fwd</text>')

    # umi_captured_rev: SSP at 3' (reversed read)
    x = COL_X[1]
    L.append(outline_block(x, y_st, 106, BH_ST, "cDNA", PAL["exon"]))
    L.append(filled_block(x+106, y_st, 50, BH_ST, "UMI", PAL["blue"]))
    L.append(filled_block(x+156, y_st, 44, BH_ST, "SSP", PAL["exon"]))
    L.append(f'<text fill="{PAL["heading"]}" font-size="9" font-weight="bold" x="{x}" y="{y_st + BH_ST + 14}">UMI captured (rev)</text>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8" x="{x}" y="{y_st + BH_ST + 26}">SSP_RC near basecalled-3′  ·  XT:i:1</text>')
    L.append(f'<text fill="{PAL["blue"]}" font-size="8" x="{x}" y="{y_st + BH_ST + 37}">XY:Z=umi_captured_rev</text>')

    # umi_not_captured: no SSP
    x = COL_X[2]
    L.append(f'<rect x="{x}" y="{y_st}" width="44" height="{BH_ST}" rx="4" '
             f'fill="{PAL["muted"]}" opacity="0.25"/>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'text-anchor="middle" x="{x+22}" y="{y_st + BH_ST//2 + 4}">5′ trunc</text>')
    L.append(outline_block(x+44, y_st, 156, BH_ST, "cDNA (partial)", PAL["exon"]))
    L.append(f'<text fill="{PAL["heading"]}" font-size="9" font-weight="bold" x="{x}" y="{y_st + BH_ST + 14}">UMI not captured</text>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8" x="{x}" y="{y_st + BH_ST + 26}">5′-truncated before SSP  ·  XT:i:2</text>')
    L.append(f'<text fill="{PAL["orange"]}" font-size="8" x="{x}" y="{y_st + BH_ST + 37}">XY:Z=umi_not_captured</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════
# Figure 2: cdna_poa_consensus (760 × 380)
# ═══════════════════════════════════════════════════════════════════════

def fig_cdna_poa_consensus():
    H = 380
    L = []
    L.append(svg_open(H))
    L.append(fig_title("UMI Consensus Pipeline"))

    # ── Pipeline flow ─────────────────────────────────────────────────
    y_sec = 42
    L.append(section_head("PIPELINE STAGES", y_sec))

    BX_H = 52
    BX_W = 110
    BX_Y  = y_sec + 16
    GAP_X = 24
    X0    = 44

    stages = [
        ("Input BAM",            PAL["muted"],  "minimap2 alignment",      "#94a3b8"),
        ("UMI Extraction",       PAL["blue"],   "SSP+27nt UMI\nextract_read_info()", "#2563eb"),
        ("Directional\nClustering", PAL["exon"], "Lev≤3 UMI\nanchor ±25bp",    "#6366f1"),
        ("abPOA\nConsensus",     PAL["teal"],   "pyabpoa\nmulti-read→1 seq",   "#0d9488"),
        ("stage1_consensus\n.fastq.gz", PAL["green"], "XU/XO/XC/XQ/XK\ntags in comment", "#059669"),
    ]

    boxes = []
    x = X0
    for name, color, sublabel, tc in stages:
        boxes.append((x, color, name, sublabel, tc))
        x += BX_W + GAP_X

    # Draw rDNA filter annotation above arrow between stage 1 and 2
    box0_cx = boxes[0][0] + BX_W / 2
    box1_cx = boxes[1][0] + BX_W / 2

    for i, (bx, color, name, sublabel, tc) in enumerate(boxes):
        # Main box
        lines = name.split("\n")
        # box fill
        L.append(f'<rect x="{bx}" y="{BX_Y}" width="{BX_W}" height="{BX_H}" rx="6" '
                 f'fill="{color}" opacity="0.15"/>')
        L.append(f'<rect x="{bx}" y="{BX_Y}" width="{BX_W}" height="{BX_H}" rx="6" '
                 f'fill="none" stroke="{color}" stroke-width="1.4"/>')
        # Label lines
        n = len(lines)
        for j, line in enumerate(lines):
            ly = BX_Y + BX_H//2 - (n-1)*7 + j*14
            L.append(f'<text fill="{tc}" font-size="9" font-weight="600" text-anchor="middle" '
                     f'x="{bx + BX_W//2}" y="{ly}">{line}</text>')
        # Sub-label below box
        sublines = sublabel.split("\n")
        for j, sl in enumerate(sublines):
            L.append(f'<text fill="{PAL["label"]}" font-size="7.5" text-anchor="middle" '
                     f'x="{bx + BX_W//2}" y="{BX_Y + BX_H + 12 + j*11}">{sl}</text>')
        # Arrow to next box
        if i < len(boxes) - 1:
            ax1 = bx + BX_W + 2
            ax2 = boxes[i+1][0] - 2
            ay  = BX_Y + BX_H // 2
            L.append(h_arrow(ax1, ax2, ay, PAL["muted"]))

    # rDNA annotation on arrow 0→1
    arr_y0 = BX_Y + BX_H // 2
    arr_mx = (boxes[0][0] + BX_W + boxes[1][0]) / 2
    L.append(f'<text fill="{PAL["red"]}" font-size="7.5" font-weight="600" '
             f'text-anchor="middle" x="{arr_mx}" y="{arr_y0 - 8}">rDNA</text>')
    L.append(f'<text fill="{PAL["red"]}" font-size="7.5" text-anchor="middle" '
             f'x="{arr_mx}" y="{arr_y0 - 0}">masked</text>')

    # ── Stats section ─────────────────────────────────────────────────
    y_div = BX_Y + BX_H + 52
    L.append(hdivider(y_div))
    y_stats = y_div + 14
    L.append(section_head("chrI TEST RUN STATS", y_stats))

    stats = [
        ("40,312 input reads → 26,364 clusters", PAL["blue"],  "(65.4% dedup rate)"),
        ("21,564 singletons  ·  4,800 POA consensus  ·  190s on chrI", PAL["label"], ""),
    ]
    for i, (text, color, suffix) in enumerate(stats):
        y_s = y_stats + 18 + i * 18
        L.append(f'<text fill="{color}" font-size="9" font-weight="600" '
                 f'text-anchor="middle" x="{FIG_W//2}" y="{y_s}">{text}</text>')
        if suffix:
            tw = len(text) * 5.4  # rough estimate
            L.append(f'<text fill="{PAL["muted"]}" font-size="8.5" '
                     f'text-anchor="middle" x="{FIG_W//2}" y="{y_s + 14}">{suffix}</text>')

    # ── Tag summary ───────────────────────────────────────────────────
    y_tags = y_stats + 72
    L.append(hdivider(y_tags))
    y_t = y_tags + 14
    L.append(section_head("OUTPUT TAGS", y_t))

    tag_defs = [
        ("XU:Z", "canonical UMI",            PAL["blue"]),
        ("XO:Z", "orient fwd/rev",           PAL["exon"]),
        ("XC:i", "cluster size",             PAL["teal"]),
        ("XT:i", "read type 1/2",            PAL["exon"]),
        ("XQ:i / XK:i", "5′ / 3′ pretrim",   PAL["orange"]),
        ("XA:i", "poly(A) length",           PAL["green"]),
    ]
    TW = 110
    tx0 = (FIG_W - len(tag_defs) * TW) / 2
    for i, (tag, desc, color) in enumerate(tag_defs):
        tx = tx0 + i * TW
        ty = y_t + 18
        L.append(f'<text fill="{color}" font-size="9" font-weight="700" '
                 f'text-anchor="middle" x="{tx + TW//2}" y="{ty}">{tag}</text>')
        L.append(f'<text fill="{PAL["label"]}" font-size="7.5" '
                 f'text-anchor="middle" x="{tx + TW//2}" y="{ty + 11}">{desc}</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════
# Figure 3: cdna_isoform_clustering (760 × 380)
# ═══════════════════════════════════════════════════════════════════════

def fig_cdna_isoform_clustering():
    H = 380
    L = []
    L.append(svg_open(H))
    L.append(fig_title("Isoform Clustering from UMI Clusters"))

    # ── Left panel: cluster dots ──────────────────────────────────────
    y_sec = 42
    L.append(section_head("CLUSTER (pos5, pos3) GROUPING", y_sec))

    # Coordinate axes
    AX_X = 70
    AX_Y = y_sec + 20
    AX_W = 280
    AX_H = 180

    L.append(f'<line x1="{AX_X}" x2="{AX_X}" y1="{AX_Y}" y2="{AX_Y+AX_H}" '
             f'stroke="{PAL["border"]}" stroke-width="1.2"/>')
    L.append(f'<line x1="{AX_X}" x2="{AX_X+AX_W}" y1="{AX_Y+AX_H}" y2="{AX_Y+AX_H}" '
             f'stroke="{PAL["border"]}" stroke-width="1.2"/>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8" text-anchor="middle" '
             f'x="{AX_X + AX_W//2}" y="{AX_Y+AX_H+14}">pos5 (TSS)</text>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8" text-anchor="middle" '
             f'transform="rotate(-90,{AX_X-12},{AX_Y+AX_H//2})" '
             f'x="{AX_X-12}" y="{AX_Y+AX_H//2}">pos3 (CPA)</text>')

    # Three groups of cluster dots
    # Map logical coords → SVG coords within axis box
    def ax_pt(p5, p3):
        fx = (p5 - 80) / 160  # 0..1 over range 80..240
        fy = 1 - (p3 - 450) / 400  # 0..1 over range 450..850 (inverted y)
        return AX_X + fx * AX_W, AX_Y + fy * AX_H

    groups = [
        # (center_p5, center_p3, scatter_offsets, color, label)
        (100, 510, [(-4,-3),(2,4),(6,-5)],   PAL["blue"],   "Group 1\n(3 clusters)"),
        (200, 760, [(-3,3),(5,-4)],           PAL["teal"],   "Group 2\n(2 clusters)"),
        (155, 615, [(-6,3),(3,-5),(1,6),(-4,-2)], PAL["orange"], "Group 3\n(4 clusters)"),
    ]

    for gp5, gp3, offsets, color, label in groups:
        cx, cy = ax_pt(gp5, gp3)
        # Tolerance ellipse
        L.append(f'<ellipse cx="{cx:.1f}" cy="{cy:.1f}" rx="22" ry="22" '
                 f'fill="{color}" opacity="0.10" stroke="{color}" stroke-width="0.8" stroke-dasharray="4,2"/>')
        # Dot cluster
        for dx, dy in offsets:
            px, py = cx + dx * 2.5, cy + dy * 2.5
            L.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5" '
                     f'fill="{color}" opacity="0.85"/>')
        # Label
        lines = label.split("\n")
        for j, line in enumerate(lines):
            ly = cy - 30 + j * 11
            L.append(f'<text fill="{color}" font-size="8" font-weight="600" '
                     f'text-anchor="middle" x="{cx}" y="{ly}">{line}</text>')
        # Tol annotation
        L.append(f'<text fill="{PAL["muted"]}" font-size="7" text-anchor="middle" '
                 f'x="{cx}" y="{cy + 32}">tol=5bp</text>')

    # ── Middle: arrow ─────────────────────────────────────────────────
    ARR_X1 = AX_X + AX_W + 12
    ARR_X2 = ARR_X1 + 44
    ARR_Y  = AX_Y + AX_H // 2
    L.append(h_arrow(ARR_X1, ARR_X2, ARR_Y, PAL["green"], "group"))

    # ── Right panel: isoform manifest table ───────────────────────────
    TB_X = ARR_X2 + 10
    TB_W = 760 - TB_X - 16
    TB_Y = y_sec + 16
    TB_ROW = 22
    COLS = ["isoform_id", "chrom", "n_clust", "n_reads", "pos5_med", "pos3_med", "gene"]
    COL_W = TB_W / len(COLS)

    # Header
    L.append(f'<rect x="{TB_X}" y="{TB_Y}" width="{TB_W}" height="{TB_ROW}" rx="4" '
             f'fill="{PAL["green"]}" opacity="0.15"/>')
    for j, col in enumerate(COLS):
        L.append(f'<text fill="{PAL["green"]}" font-size="7.5" font-weight="600" '
                 f'text-anchor="middle" x="{TB_X + (j + 0.5)*COL_W}" y="{TB_Y + 14}">{col}</text>')

    # Example rows
    rows = [
        ["ISO_001", "chrI", "3", "12", "100", "510", "YAL001C"],
        ["ISO_002", "chrI", "2", "7",  "200", "760", "YAL002W"],
        ["ISO_003", "chrI", "4", "18", "155", "615", "YAL001C"],
    ]
    row_colors = [PAL["blue"], PAL["teal"], PAL["orange"]]
    for i, row in enumerate(rows):
        ry = TB_Y + (i + 1) * TB_ROW
        stripe = "0.05" if i % 2 == 0 else "0"
        L.append(f'<rect x="{TB_X}" y="{ry}" width="{TB_W}" height="{TB_ROW}" '
                 f'fill="{row_colors[i]}" opacity="{stripe}"/>')
        for j, val in enumerate(row):
            L.append(f'<text fill="{PAL["label"]}" font-size="8" text-anchor="middle" '
                     f'x="{TB_X + (j + 0.5)*COL_W}" y="{ry + 14}">{val}</text>')

    L.append(f'<rect x="{TB_X}" y="{TB_Y}" width="{TB_W}" height="{(len(rows)+1)*TB_ROW}" '
             f'rx="4" fill="none" stroke="{PAL["border"]}" stroke-width="0.8"/>')

    # ── Section 2: T1/T2 reconciliation ──────────────────────────────
    y_div = AX_Y + AX_H + 32
    L.append(hdivider(y_div))
    y_sec2 = y_div + 14
    L.append(section_head("T1 ⇔ T2 RECONCILIATION", y_sec2, PAL["teal"]))

    y_t1t2 = y_sec2 + 18
    BH2 = 28
    # T1 box (orange)
    L.append(flow_box(90, y_t1t2, 200, BH2, "Type-1 (pos5=TSS)", PAL["orange"]))
    L.append(f'<text fill="{PAL["label"]}" font-size="7.5" text-anchor="middle" '
             f'x="190" y="{y_t1t2 + BH2 + 12}">UMI-anchored, full-length</text>')
    # T2 box (blue)
    L.append(flow_box(470, y_t1t2, 200, BH2, "Type-2 (pos3=CPA only)", PAL["blue"]))
    L.append(f'<text fill="{PAL["label"]}" font-size="7.5" text-anchor="middle" '
             f'x="570" y="{y_t1t2 + BH2 + 12}">pos5 unknown (truncated)</text>')

    # Bidirectional arrow
    MID_Y = y_t1t2 + BH2 // 2
    L.append(f'<line x1="290" x2="470" y1="{MID_Y}" y2="{MID_Y}" '
             f'stroke="{PAL["teal"]}" stroke-width="1.5"/>')
    L.append(f'<polygon fill="{PAL["teal"]}" points="470,{MID_Y} 464,{MID_Y-3.5} 464,{MID_Y+3.5}"/>')
    L.append(f'<polygon fill="{PAL["teal"]}" points="290,{MID_Y} 296,{MID_Y-3.5} 296,{MID_Y+3.5}"/>')
    L.append(f'<text fill="{PAL["teal"]}" font-size="8" font-weight="600" '
             f'text-anchor="middle" x="380" y="{MID_Y - 8}">same anchor ±5bp → XL partner tag</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════
# Figure 4: walkback_readvsref (760 × 380)
# ═══════════════════════════════════════════════════════════════════════

def fig_walkback_readvsref():
    H = 380
    L = []
    L.append(svg_open(H))
    L.append(fig_title("Read-vs-Reference Poly-A Walkback Algorithm"))

    # ── Main walkback diagram ─────────────────────────────────────────
    y_sec = 42
    L.append(section_head("WALKBACK: STOP AT FIRST NON-A MATCH", y_sec))

    X0   = 88
    y_ref = y_sec + 18
    y_rd  = y_ref + BH + 28

    # Reference: ...G C T | A A A A A | G T
    ref_bases  = list("GCTAAAAAGT")
    a_start = 3  # index of first genomic A
    a_end   = 7  # index of last genomic A (inclusive)

    ref_hi = {i: (PAL["A_f"], PAL["green"]) for i in range(a_start, a_end+1)}
    L.append(f'<text fill="{PAL["label"]}" font-size="9" x="14" y="{y_ref + BH//2 + 4}">Reference</text>')
    L.append(nt_row(ref_bases, X0, y_ref, ref_hi))
    # Label A-tract above
    L.append(f'<path d="M{X0+a_start*BW},{y_ref-2} L{X0+a_start*BW},{y_ref-8} '
             f'L{X0+(a_end+1)*BW},{y_ref-8} L{X0+(a_end+1)*BW},{y_ref-2}" '
             f'fill="none" stroke="{PAL["green"]}" stroke-width="1" opacity="0.7"/>')
    L.append(f'<text fill="{PAL["green"]}" font-size="8" text-anchor="middle" '
             f'x="{X0 + (a_start + a_end + 1) * BW / 2}" y="{y_ref - 11}">genomic A-tract</text>')

    # Read: G C T A A A A A A A A A (extends through A-tract into poly-A)
    rd_bases = list("GCTAAAAAAAAA")
    # Poly-A starts at position a_end+1 = 8
    poly_start = a_end + 1  # 8
    rd_hi = {}
    for i in range(a_start, len(rd_bases)):
        rd_hi[i] = (PAL["A_f"], PAL["red_l"] if i >= poly_start else PAL["border"])

    L.append(f'<text fill="{PAL["label"]}" font-size="9" x="14" y="{y_rd + BH//2 + 4}">Read</text>')
    L.append(nt_row(rd_bases, X0, y_rd, rd_hi))
    # Label poly-A tail below
    L.append(brace_below(X0 + poly_start*BW, X0 + len(rd_bases)*BW, y_rd + BH,
                         PAL["red"], "poly(A) extends past genomic A-tract"))

    # Walk-back arrow
    y_arrow = y_rd + BH + 38
    arr_r = X0 + (len(rd_bases) - 1) * BW
    arr_l = X0 + (a_start - 1) * BW
    L.append(h_arrow(arr_r, arr_l, y_arrow, PAL["blue"], "walk back past A’s"))

    # True CPA marker
    cpa_x = X0 + a_start * BW
    L.append(vert_marker(cpa_x, y_ref - 14, y_rd + BH + 5, PAL["green"], "true CPA", "left"))

    # Stop rule caption
    y_rule = y_arrow + 18
    L.append(f'<text fill="{PAL["heading"]}" font-size="8.5" font-style="italic" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_rule}">'
             f'Stop at first position where read_base == ref_base AND ref_base ≠ A</text>')

    # Divider
    y_div = y_rule + 16
    L.append(hdivider(y_div))

    # ── Bug note ──────────────────────────────────────────────────────
    y_sec2 = y_div + 14
    L.append(section_head("CRITICAL: MUST NOT SKIP WALKBACK WHEN ref=A", y_sec2, PAL["red"]))

    y_bug = y_sec2 + 18
    BH_BUG = 30

    # Red X — wrong behavior
    L.append(f'<rect x="90" y="{y_bug}" width="240" height="{BH_BUG}" rx="5" '
             f'fill="{PAL["red"]}" opacity="0.10"/>')
    L.append(f'<rect x="90" y="{y_bug}" width="240" height="{BH_BUG}" rx="5" '
             f'fill="none" stroke="{PAL["red"]}" stroke-width="1.2"/>')
    L.append(f'<text fill="{PAL["red"]}" font-size="9" font-weight="600" x="96" y="{y_bug + 13}">'
             f'✗ OLD BUG: skip walkback if ref=A</text>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8" x="96" y="{y_bug + 25}">'
             f'Leaves read in internal-priming position</text>')

    # Green checkmark — correct behavior
    L.append(f'<rect x="380" y="{y_bug}" width="290" height="{BH_BUG}" rx="5" '
             f'fill="{PAL["green"]}" opacity="0.10"/>')
    L.append(f'<rect x="380" y="{y_bug}" width="290" height="{BH_BUG}" rx="5" '
             f'fill="none" stroke="{PAL["green"]}" stroke-width="1.2"/>')
    L.append(f'<text fill="{PAL["green"]}" font-size="9" font-weight="600" x="386" y="{y_bug + 13}">'
             f'✓ CORRECT: always walk back over A’s</text>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8" x="386" y="{y_bug + 25}">'
             f'Genomic A matching read A = internal-priming scenario</text>')

    # Protocol note
    y_note = y_bug + BH_BUG + 14
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_note}">'
             f'Protocol-agnostic: applies to QuantSeq REV, DRS, and PCR-cDNA  —  '
             f'fixed in walkback.py 2026-05-13</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════
# Figure 5: splice_classification (760 × 380)
# ═══════════════════════════════════════════════════════════════════════

def fig_splice_classification():
    H = 380
    L = []
    L.append(svg_open(H))
    L.append(fig_title("Per-Junction Splice Classification"))

    # ── Three junction type columns ───────────────────────────────────
    y_sec = 42
    L.append(section_head("JUNCTION TYPES", y_sec))

    COL_W = 200
    COL_H = 130
    GAP   = 30
    COL_X = [46, 46 + COL_W + GAP, 46 + 2*(COL_W + GAP)]
    y_col = y_sec + 16

    junc_defs = [
        ("GT-AG", "+5", "canonical",       PAL["blue"],   PAL["blue_l"],
         "GT", "AG", "Majority of annotated\njunctions in yeast"),
        ("AT-AC", "+3", "minor canonical",  PAL["teal"],   PAL["teal_l"],
         "AT", "AC", "Minor spliceosome;\nU12-type introns"),
        ("non-canonical", "−3", "non-canonical", PAL["red"],    PAL["red_l"],
         "??", "??", "Novel / spurious;\nflag for review"),
    ]

    for i, (dinuc, score, cat_label, color, light, donor, acceptor, note) in enumerate(junc_defs):
        x = COL_X[i]
        # Column box
        L.append(f'<rect x="{x}" y="{y_col}" width="{COL_W}" height="{COL_H}" rx="6" '
                 f'fill="{color}" opacity="0.10"/>')
        L.append(f'<rect x="{x}" y="{y_col}" width="{COL_W}" height="{COL_H}" rx="6" '
                 f'fill="none" stroke="{color}" stroke-width="1.2"/>')
        # Dinucleotide title
        L.append(f'<text fill="{color}" font-size="13" font-weight="bold" '
                 f'text-anchor="middle" x="{x + COL_W//2}" y="{y_col + 22}">{dinuc}</text>')
        # Score badge
        L.append(f'<rect x="{x + COL_W//2 - 18}" y="{y_col + 28}" width="36" height="16" '
                 f'rx="8" fill="{color}" opacity="0.2"/>')
        L.append(f'<text fill="{color}" font-size="9" font-weight="600" '
                 f'text-anchor="middle" x="{x + COL_W//2}" y="{y_col + 40}">score {score}</text>')
        # Category label
        L.append(f'<text fill="{color}" font-size="8.5" font-weight="600" '
                 f'text-anchor="middle" x="{x + COL_W//2}" y="{y_col + 56}">{cat_label}</text>')

        # Donor/Acceptor nucleotide boxes
        nb_y = y_col + 64
        nb_x_donor = x + COL_W//2 - 40
        nb_x_acc   = x + COL_W//2 + 8
        # donor dinucleotide
        for k, b in enumerate(donor):
            fill_b = BASE_FILL.get(b, light)
            text_b = BASE_TEXT.get(b, color)
            L.append(f'<rect x="{nb_x_donor + k*BW}" y="{nb_y}" width="{BW}" height="{BH}" '
                     f'rx="3" fill="{fill_b}" stroke="{color}" stroke-width="0.8"/>')
            L.append(f'<text fill="{text_b}" font-size="11" font-weight="600" '
                     f'text-anchor="middle" x="{nb_x_donor + k*BW + BW//2}" '
                     f'y="{nb_y + BH//2 + 4}">{b}</text>')
        # dashes for intron
        L.append(dashed_line(nb_x_donor + 2*BW + 2, nb_x_acc - 2, nb_y + BH//2, PAL["muted"]))
        # acceptor dinucleotide
        for k, b in enumerate(acceptor):
            fill_b = BASE_FILL.get(b, light)
            text_b = BASE_TEXT.get(b, color)
            L.append(f'<rect x="{nb_x_acc + k*BW}" y="{nb_y}" width="{BW}" height="{BH}" '
                     f'rx="3" fill="{fill_b}" stroke="{color}" stroke-width="0.8"/>')
            L.append(f'<text fill="{text_b}" font-size="11" font-weight="600" '
                     f'text-anchor="middle" x="{nb_x_acc + k*BW + BW//2}" '
                     f'y="{nb_y + BH//2 + 4}">{b}</text>')
        # donor / acceptor labels
        L.append(f'<text fill="{PAL["muted"]}" font-size="7" text-anchor="middle" '
                 f'x="{nb_x_donor + BW}" y="{nb_y + BH + 11}">donor</text>')
        L.append(f'<text fill="{PAL["muted"]}" font-size="7" text-anchor="middle" '
                 f'x="{nb_x_acc + BW}" y="{nb_y + BH + 11}">acceptor</text>')

        # Note lines
        for j, line in enumerate(note.split("\n")):
            L.append(f'<text fill="{PAL["label"]}" font-size="7.5" text-anchor="middle" '
                     f'x="{x + COL_W//2}" y="{y_col + COL_H - 14 + j*11}">{line}</text>')

    # ── Chimeric consensus voting ─────────────────────────────────────
    y_div = y_col + COL_H + 22
    L.append(hdivider(y_div))
    y_sec2 = y_div + 14
    L.append(section_head("CHIMERIC CONSENSUS VOTING  (select_best_chimeric())", y_sec2, PAL["teal"]))

    y_vote = y_sec2 + 18
    BH_V = 26

    # mapPacBio wins
    L.append(flow_box(80, y_vote, 210, BH_V, "mapPacBio calls GT-AG", PAL["orange"]))
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" text-anchor="middle" '
             f'x="185" y="{y_vote + BH_V + 12}">score +5</text>')

    # vs
    L.append(f'<text fill="{PAL["heading"]}" font-size="9" font-weight="600" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_vote + BH_V//2 + 4}">vs</text>')

    # minimap2 loses
    L.append(flow_box(470, y_vote, 210, BH_V, "minimap2 calls GT-TT", PAL["muted"]))
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" text-anchor="middle" '
             f'x="575" y="{y_vote + BH_V + 12}">score −3</text>')

    # Arrow → winner
    y_win = y_vote + BH_V + 30
    L.append(v_arrow(FIG_W//2, y_vote + BH_V + 16, y_win, PAL["teal"]))
    L.append(f'<rect x="270" y="{y_win}" width="220" height="{BH_V}" rx="5" '
             f'fill="{PAL["green"]}" opacity="0.15"/>')
    L.append(f'<rect x="270" y="{y_win}" width="220" height="{BH_V}" rx="5" '
             f'fill="none" stroke="{PAL["green"]}" stroke-width="1.2"/>')
    L.append(f'<text fill="{PAL["green"]}" font-size="9" font-weight="600" '
             f'text-anchor="middle" x="380" y="{y_win + BH_V//2 + 4}">'
             f'✓ select_best_chimeric() picks mapPacBio</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════
# Figure: cdna_pipeline_overview (760 × 300)
#   Top-level 3-stage view of the ONT PCR-cDNA workflow:
#     correct-cdna  →  rectify align  →  cdna-analyze
# ═══════════════════════════════════════════════════════════════════════

def fig_cdna_pipeline_overview():
    H = 300
    L = []
    L.append(svg_open(H))
    L.append(fig_title("ONT PCR-cDNA Pipeline (PCB114.24)"))

    # Three stage boxes
    STAGES = [
        dict(
            name="rectify correct-cdna",
            color=PAL["blue"],
            sub=("UMI extraction · directional clustering",
                 "abPOA consensus · pre-trim"),
            input="pre-aligned BAM",
            output="stage1_consensus.fastq.gz",
            tags="XU · XO · XC · XR · XM · XF · XA · XT · XY · XQ · XK · XB",
        ),
        dict(
            name="rectify align",
            color=PAL["teal"],
            sub=("minimap2 + mapPacBio + gapmm2",
                 "chimeric consensus selection"),
            input="per-cluster FASTQ",
            output="{prefix}.rectified.bam",
            tags="X[upper] tags pass through via minimap2 -y",
        ),
        dict(
            name="rectify cdna-analyze",
            color=PAL["green"],
            sub=("walkback (XA) · TSS walk-forward",
                 "isoform clustering (XI) · T1↔T2 (XL)"),
            input="post-align BAM",
            output="clusters.tsv · isoforms.tsv · t1t2_pairs.tsv · consensus_tagged.bam",
            tags="+ XG (gene) · XS (sense/antisense) · XI · XL",
        ),
    ]

    BX_W = 210
    BX_H = 64
    GAP_X = 12
    X0 = (FIG_W - 3 * BX_W - 2 * GAP_X) // 2
    BX_Y = 52

    for i, st in enumerate(STAGES):
        bx = X0 + i * (BX_W + GAP_X)
        # Box
        L.append(f'<rect x="{bx}" y="{BX_Y}" width="{BX_W}" height="{BX_H}" rx="8" '
                 f'fill="{st["color"]}" opacity="0.13"/>')
        L.append(f'<rect x="{bx}" y="{BX_Y}" width="{BX_W}" height="{BX_H}" rx="8" '
                 f'fill="none" stroke="{st["color"]}" stroke-width="1.6"/>')
        # Title
        L.append(f'<text fill="{st["color"]}" font-size="11" font-weight="700" '
                 f'text-anchor="middle" x="{bx + BX_W//2}" y="{BX_Y + 18}">'
                 f'{st["name"]}</text>')
        # Sub-lines
        for j, sub in enumerate(st["sub"]):
            L.append(f'<text fill="{PAL["heading"]}" font-size="8.5" '
                     f'text-anchor="middle" x="{bx + BX_W//2}" '
                     f'y="{BX_Y + 36 + j * 12}">{sub}</text>')

        # Input label above
        L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
                 f'text-anchor="middle" x="{bx + BX_W//2}" y="{BX_Y - 6}">'
                 f'in: {st["input"]}</text>')

        # Output label below
        L.append(f'<text fill="{st["color"]}" font-size="8.5" font-weight="600" '
                 f'text-anchor="middle" x="{bx + BX_W//2}" '
                 f'y="{BX_Y + BX_H + 14}">{st["output"]}</text>')

        # Tag list (smaller, muted)
        L.append(f'<text fill="{PAL["label"]}" font-size="7" '
                 f'text-anchor="middle" x="{bx + BX_W//2}" '
                 f'y="{BX_Y + BX_H + 28}">{st["tags"]}</text>')

        # Arrow to next stage
        if i < len(STAGES) - 1:
            ax1 = bx + BX_W + 1
            ax2 = X0 + (i + 1) * (BX_W + GAP_X) - 1
            ay = BX_Y + BX_H // 2
            L.append(h_arrow(ax1, ax2, ay, PAL["muted"]))

    # Bottom: namespace note
    y_note = BX_Y + BX_H + 72
    L.append(hdivider(y_note))
    L.append(f'<text fill="{PAL["heading"]}" font-size="9" font-weight="700" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_note + 18}">'
             f'Tag namespace</text>')
    L.append(f'<text fill="{PAL["label"]}" font-size="8.5" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_note + 34}">'
             f'X[upper] = upstream / user-visible metadata  ·  '
             f'X[lower] = rectify align internal aligner-selection bookkeeping (Xa, Xc, Xn, Xj, Xv, Xz, Xg, Xm, Xq, Xw, Xy)</text>')
    L.append(f'<text fill="{PAL["muted"]}" font-size="8" font-style="italic" '
             f'text-anchor="middle" x="{FIG_W//2}" y="{y_note + 48}">'
             f'FASTQ comment tags are TAB-separated so minimap2 -y emits each as a separate BAM aux field.</text>')

    L.append(svg_close())
    return "\n".join(L)


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    figures = {
        "cdna_pipeline_overview":  fig_cdna_pipeline_overview,
        "cdna_umi_architecture":   fig_cdna_umi_architecture,
        "cdna_poa_consensus":      fig_cdna_poa_consensus,
        "cdna_isoform_clustering": fig_cdna_isoform_clustering,
        "walkback_readvsref":      fig_walkback_readvsref,
        "splice_classification":   fig_splice_classification,
    }

    for name, fn in figures.items():
        svg_path = OUTDIR / f"{name}.svg"
        png_path = OUTDIR / f"{name}.png"
        svg_text = fn()
        svg_path.write_text(svg_text, encoding="utf-8")
        cairosvg.svg2png(
            url=str(svg_path),
            write_to=str(png_path),
            output_width=760,
            scale=2.0,
        )
        svg_kb = svg_path.stat().st_size // 1024
        png_kb = png_path.stat().st_size // 1024
        print(f"  {name}: SVG {svg_kb}KB  PNG {png_kb}KB")

    print(f"\nDone — {len(figures)} SVG + {len(figures)} PNG written to {OUTDIR}")
