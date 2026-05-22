#!/usr/bin/env python3
"""
Hero pipeline_overview figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, Inter font stack, restrained PAL palette
  - Minimum text — viewer should grok the structure in <5 seconds
  - 3 chemistries converge → shared multi-aligner → correct (5'/intron/3')
    → rectified BAM → analyze + DESeq2
  - No success banners, no marketing framing, descriptive headers only
"""

import os
import subprocess

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 760
FIG_H = 640
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
    orange  = "#d97706",
    orange_l= "#fed7aa",
    indigo  = "#4f46e5",
    indigo_l= "#e0e7ff",
    slate   = "#334155",
    slate_l = "#f1f5f9",
)

# Per-chemistry track colors
CHEM = {
    "drs":   dict(color=PAL["blue"],   fill=PAL["blue_l"],   name="ONT DRS",        sub="direct RNA"),
    "cdna":  dict(color=PAL["teal"],   fill=PAL["teal_l"],   name="ONT PCR-cDNA",   sub="UMI, full-length"),
    "qsrev": dict(color=PAL["orange"], fill=PAL["orange_l"], name="QuantSeq REV",   sub="short read, antisense"),
}

# Column centers for the three input lanes
COL = {"drs": 135, "cdna": 380, "qsrev": 625}
CHIP_W = 180
CHIP_H = 38


# ─────────────────────────────────────────────────────────────────────────────
# Drawing helpers
# ─────────────────────────────────────────────────────────────────────────────

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

def text(x, y, s, size=10, color=None, weight="400", anchor="start", letter_spacing=None):
    c = color or PAL["label"]
    ls = f' letter-spacing="{letter_spacing}"' if letter_spacing else ""
    return (f'<text fill="{c}" font-size="{size}" font-weight="{weight}" '
            f'text-anchor="{anchor}"{ls} x="{x}" y="{y}">{s}</text>')

def section_head(x, y, s):
    return text(x, y, s, size=11, color=PAL["heading"], weight="700",
                letter_spacing="0.8")

def chip(cx, cy, w, h, label, color, fill, sub=None):
    """Rounded chip with label centered, optional subtitle below."""
    x = cx - w / 2
    y = cy - h / 2
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" '
        f'fill="{fill}" stroke="{color}" stroke-width="1.2"/>',
        text(cx, cy + 4, label, size=12, color=color, weight="600", anchor="middle"),
    ]
    if sub:
        parts.append(text(cx, cy + h / 2 + 13, sub, size=11,
                          color=PAL["muted"], anchor="middle"))
    return "\n".join(parts)

def shared_box(cx, cy, w, h, title, sub=None, title_y=None):
    """Wider neutral box for shared stages, with optional subtitle.

    By default the title sits near vertical center. Pass `title_y` (absolute
    y-coordinate) to override — needed when the box also contains inner pills
    that need clearance below the title.
    """
    x = cx - w / 2
    y = cy - h / 2
    if title_y is None:
        title_y = cy + (2 if sub else 4)
    parts = [
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" '
        f'fill="{PAL["slate_l"]}" stroke="{PAL["slate"]}" stroke-width="1.2"/>',
        text(cx, title_y, title, size=12,
             color=PAL["slate"], weight="700", anchor="middle"),
    ]
    if sub:
        parts.append(text(cx, title_y + 14, sub, size=11,
                          color=PAL["muted"], anchor="middle"))
    return "\n".join(parts)

def pill(cx, cy, w, h, label, color, fill):
    """Smaller secondary pill for sub-stages inside a parent box."""
    x = cx - w / 2
    y = cy - h / 2
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" '
        f'fill="{fill}" stroke="{color}" stroke-width="1"/>\n'
        + text(cx, cy + 4, label, size=12, color=color, weight="600", anchor="middle")
    )

def arrow(x1, y1, x2, y2, color=None):
    c = color or PAL["label"]
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{c}" '
            f'stroke-width="1.4" marker-end="url(#arr)"/>')

def hdivider(y, x1=60, x2=700, color=None):
    c = color or PAL["divider"]
    return f'<line x1="{x1}" x2="{x2}" y1="{y}" y2="{y}" stroke="{c}" stroke-width="1"/>'


# ─────────────────────────────────────────────────────────────────────────────
# Build the figure
# ─────────────────────────────────────────────────────────────────────────────

def build():
    out = [svg_open(FIG_H)]

    # Title block
    out.append(text(FIG_W / 2, 28, "RECTIFY", size=18,
                    color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W / 2, 48,
                    "one pipeline · three RNA chemistries",
                    size=11.5, color=PAL["muted"], anchor="middle"))

    # ── INPUTS ───────────────────────────────────────────────────────────────
    y_inputs_head = 78
    out.append(section_head(60, y_inputs_head, "INPUTS"))

    y_chip = 100
    for key in ("drs", "cdna", "qsrev"):
        c = CHEM[key]
        out.append(chip(COL[key], y_chip, CHIP_W, CHIP_H,
                        c["name"], c["color"], c["fill"], sub=c["sub"]))

    # ── PRE-PROCESS (per-protocol pills, just below each chip) ───────────────
    y_pre_label = 153
    out.append(text(COL["drs"], y_pre_label, "trim-polya",
                    size=11, color=PAL["blue"], weight="600", anchor="middle"))
    out.append(text(COL["cdna"], y_pre_label, "correct-cdna",
                    size=11, color=PAL["teal"], weight="600", anchor="middle"))
    out.append(text(COL["qsrev"], y_pre_label, "(no pre-process)",
                    size=11, color=PAL["muted"], weight="400", anchor="middle"))
    out.append(text(FIG_W / 2, y_pre_label + 14,
                    "per-protocol — strip adapters / UMIs / poly(A)",
                    size=11, color=PAL["muted"], anchor="middle"))

    # ── LAYOUT GEOMETRY ──────────────────────────────────────────────────────
    # Single ARROW_GAP applied between every adjacent box → uniform visual
    # rhythm. Arrowhead tips land exactly on box top edges (marker refX=6,
    # userSpaceOnUse — see svg_open()).
    cx_funnel = FIG_W / 2
    ARROW_GAP = 24

    # Multi-aligner box
    aligner_h = 44
    y_aligner_top = 218
    y_aligner_bot = y_aligner_top + aligner_h          # 262
    y_aligner     = y_aligner_top + aligner_h / 2      # center, 240

    # Junction pool box
    jpool_h = 32
    y_jpool_top = y_aligner_bot + ARROW_GAP            # 286
    y_jpool_bot = y_jpool_top + jpool_h                # 318
    y_jpool     = y_jpool_top + jpool_h / 2            # 302

    # Correct box
    correct_h = 88
    y_correct_top = y_jpool_bot + ARROW_GAP            # 342
    y_correct_bot = y_correct_top + correct_h          # 430
    y_correct     = y_correct_top + correct_h / 2      # 386

    # Consensus box
    consensus_h = 32
    y_consensus_top = y_correct_bot + ARROW_GAP        # 454
    y_consensus_bot = y_consensus_top + consensus_h    # 486
    y_consensus     = y_consensus_top + consensus_h / 2  # 470

    # Rectified BAM box
    bam_h = 28
    y_bam_top = y_consensus_bot + ARROW_GAP            # 510
    y_bam_bot = y_bam_top + bam_h                      # 538
    y_bam     = y_bam_top + bam_h / 2                  # 524

    # ── CONVERGING ARROWS (3 inputs → multi-aligner) ─────────────────────────
    # Each arrow drops straight down on the column center and lands on the
    # aligner box top edge — uniform ARROW_GAP length matches downstream arrows.
    y_arrow_from = y_aligner_top - ARROW_GAP           # 194
    for key in ("drs", "cdna", "qsrev"):
        out.append(arrow(COL[key], y_arrow_from, COL[key], y_aligner_top))

    # ── MULTI-ALIGNER (parallel, NOT consensus — consensus happens later) ────
    out.append(shared_box(cx_funnel, y_aligner, 520, aligner_h,
                          "rectify align — run aligners in parallel",
                          sub="minimap2 + mapPacBio + gapmm2 · or bbmap + bwa-mem (--short-read) · each produces its own BAM"))

    # arrow down to junction pool
    out.append(arrow(cx_funnel, y_aligner_bot, cx_funnel, y_jpool_top))

    # ── JUNCTION POOL ────────────────────────────────────────────────────────
    out.append(shared_box(cx_funnel, y_jpool, 420, jpool_h,
                          "junction pool — observed across aligners + annotated DB",
                          sub=None))

    # arrow down to correct
    out.append(arrow(cx_funnel, y_jpool_bot, cx_funnel, y_correct_top))

    # ── CORRECT (5' / introns / 3' sub-pills, applied per aligner) ───────────
    out.append(shared_box(cx_funnel, y_correct, 560, correct_h,
                          "rectify correct — applied per aligner (standardized scoring)",
                          sub=None,
                          title_y=y_correct - 22))
    # inner sub-pills below the title
    y_subpills = y_correct + 12
    sub_w = 130
    sub_h = 26
    sub_specs = [
        (cx_funnel - 165, "5' ends",  PAL["indigo"], PAL["indigo_l"]),
        (cx_funnel,       "introns",  PAL["indigo"], PAL["indigo_l"]),
        (cx_funnel + 165, "3' ends",  PAL["indigo"], PAL["indigo_l"]),
    ]
    for sx, label, color, fill in sub_specs:
        out.append(pill(sx, y_subpills, sub_w, sub_h, label, color, fill))
    # tiny footnote inside box
    out.append(text(cx_funnel, y_correct + 38,
                    "junction pool informs 5′ rescue + intron refinement",
                    size=11, color=PAL["muted"], anchor="middle"))

    # arrow down to consensus
    out.append(arrow(cx_funnel, y_correct_bot, cx_funnel, y_consensus_top))

    # ── CONSENSUS (pick best CORRECTED alignment per read) ───────────────────
    out.append(shared_box(cx_funnel, y_consensus, 380, consensus_h,
                          "consensus — pick best corrected alignment per read",
                          sub=None))

    # arrow down to rectified BAM
    out.append(arrow(cx_funnel, y_consensus_bot, cx_funnel, y_bam_top))

    # ── RECTIFIED BAM ────────────────────────────────────────────────────────
    out.append(shared_box(cx_funnel, y_bam, 220, bam_h,
                          "rectified BAM", sub=None))

    # divider before downstream
    out.append(hdivider(y_bam_bot + 16))

    # ── ANALYZE + DESEQ2 (two rows of chips) ─────────────────────────────────
    y_down_head = 575
    out.append(section_head(60, y_down_head, "DOWNSTREAM"))

    # row 1: analyze outputs
    y_an = 595
    out.append(text(60, y_an + 3, "analyze", size=11,
                    color=PAL["heading"], weight="600"))
    analyze_items = [
        (240, "5'/3' clusters"),
        (380, "splice tables"),
        (530, "isoforms (cDNA)"),
    ]
    for cx, label in analyze_items:
        out.append(pill(cx, y_an, 130, 24, label, PAL["slate"], PAL["slate_l"]))

    # row 2: DESeq2
    y_de = 625
    out.append(text(60, y_de + 3, "DESeq2", size=11,
                    color=PAL["heading"], weight="600"))
    de_items = [
        (170, "gene"),
        (290, "3' cluster"),
        (425, "5' cluster*"),
        (570, "splice junction*"),
    ]
    for cx, label in de_items:
        out.append(pill(cx, y_de, 110, 24, label, PAL["slate"], PAL["slate_l"]))

    out.append(text(700, y_de + 3, "* planned",
                    size=10, color=PAL["muted"], anchor="end"))

    out.append(svg_close())
    return "\n".join(out)


def render(svg_str, name):
    svg_path = os.path.join(OUTDIR, f"{name}.svg")
    png_path = os.path.join(OUTDIR, f"{name}.png")
    with open(svg_path, "w") as f:
        f.write(svg_str)
    # Render at 2x for retina
    subprocess.run(
        ["cairosvg", svg_path, "-o", png_path, "--output-width", str(FIG_W * 2)],
        check=True,
    )
    print(f"wrote {svg_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    render(build(), "pipeline_overview")
