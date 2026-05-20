#!/usr/bin/env python3
"""
cDNA isoform clustering figure for the RECTIFY README — v3.

Design intent (per dev/figures/STYLE_GUIDE.md):
  - 760px wide, single panel
  - Show Type-1 (full-length, UMI-anchored) and Type-2 (truncated, no UMI)
    reads on a shared genome track
  - T1 reads cluster by (TSS, CPA); T2 reads cluster by CPA only
  - Same-molecule T1↔T2 pairs are linked by dashed line at matching 3' end
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
    blue    = "#2563eb",   blue_l   = "#dbeafe",
    teal    = "#0d9488",   teal_l   = "#ccfbf1",
    green   = "#059669",   green_l  = "#d1fae5",
    indigo  = "#4f46e5",   indigo_l = "#e0e7ff",
    purple  = "#9333ea",   purple_l = "#f3e8ff",
    orange  = "#d97706",   orange_l = "#fed7aa",
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

def gene_block(x, y, w, h):
    """Filled gene reference block."""
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" '
            f'fill="{PAL["indigo"]}"/>\n'
            + text(x + w/2, y + h/2 + 4, "Gene", size=10,
                   color="#ffffff", weight="600", anchor="middle"))

def read_bar(x1, x2, y, h, color, fill):
    return (f'<rect x="{x1}" y="{y - h/2}" width="{x2-x1}" height="{h}" rx="2" '
            f'fill="{fill}" stroke="{color}" stroke-width="1"/>')

def end_marker(x, y, h, color):
    """Solid vertical tick marking a read end."""
    return (f'<line x1="{x}" x2="{x}" y1="{y - h/2 - 2}" y2="{y + h/2 + 2}" '
            f'stroke="{color}" stroke-width="2"/>')


def build():
    out = [svg_open(FIG_H)]

    out.append(text(FIG_W/2, 26, "cDNA isoform clustering",
                    size=14, color=PAL["title"], weight="700", anchor="middle"))
    out.append(text(FIG_W/2, 44,
                    "T1 cluster by (TSS, CPA)  ·  T2 cluster by CPA only  ·  pairs linked at the same 3′ end",
                    size=10, color=PAL["muted"], anchor="middle"))

    # ── Genome reference ─────────────────────────────────────────────────────
    label_x = 70
    y_gene = 78
    gene_x = 110
    gene_w = 560
    gene_h = 26
    out.append(text(label_x, y_gene + gene_h/2 + 4, "Genome",
                    size=10, color=PAL["label"]))
    out.append(gene_block(gene_x, y_gene, gene_w, gene_h))

    # Two CPA cluster positions (3' end) — shared by T1 and T2
    cpa_a = gene_x + 380   # CPA cluster A
    cpa_b = gene_x + 510   # CPA cluster B (alternate)
    # TSS positions (5' end) — only used by T1
    tss = gene_x + 30

    # Mark CPA positions on the gene reference
    for cx, label in [(cpa_a, "CPA A"), (cpa_b, "CPA B")]:
        out.append(f'<line x1="{cx}" x2="{cx}" y1="{y_gene - 4}" '
                   f'y2="{y_gene + gene_h + 4}" stroke="{PAL["muted"]}" '
                   f'stroke-width="1" stroke-dasharray="2,2" opacity="0.6"/>')
        out.append(text(cx, y_gene - 8, label, size=8.5,
                        color=PAL["muted"], anchor="middle"))

    # Mark TSS position
    out.append(f'<line x1="{tss}" x2="{tss}" y1="{y_gene - 4}" '
               f'y2="{y_gene + gene_h + 4}" stroke="{PAL["muted"]}" '
               f'stroke-width="1" stroke-dasharray="2,2" opacity="0.6"/>')
    out.append(text(tss, y_gene - 8, "TSS", size=8.5,
                    color=PAL["muted"], anchor="middle"))

    # ── Type-1 reads (full-length, UMI-anchored) ─────────────────────────────
    # T1 reads span TSS → CPA. Two clusters: T1@(TSS, CPA_A) and T1@(TSS, CPA_B).
    y_t1_label = 132
    out.append(text(label_x, y_t1_label, "Type 1",
                    size=10, color=PAL["teal"], weight="700"))
    out.append(text(label_x, y_t1_label + 12, "(full-length)",
                    size=8.5, color=PAL["muted"]))

    t1_reads = [
        (132, cpa_a),
        (143, cpa_a),
        (154, cpa_b),
        (165, cpa_b),
    ]
    for y, cpa in t1_reads:
        out.append(read_bar(tss, cpa, y, 6, PAL["teal"], PAL["teal_l"]))
        out.append(end_marker(tss, y, 6, PAL["teal"]))
        out.append(end_marker(cpa, y, 6, PAL["teal"]))

    # ── Type-2 reads (truncated, no UMI — CPA-anchored only) ─────────────────
    y_t2_label = 198
    out.append(text(label_x, y_t2_label, "Type 2",
                    size=10, color=PAL["orange"], weight="700"))
    out.append(text(label_x, y_t2_label + 12, "(truncated)",
                    size=8.5, color=PAL["muted"]))

    # T2 reads have random 5' start, fixed 3' end at CPA_A or CPA_B
    t2_reads = [
        (198, gene_x + 200, cpa_a),
        (209, gene_x + 280, cpa_a),
        (220, gene_x + 250, cpa_a),
        (231, gene_x + 320, cpa_b),
        (242, gene_x + 360, cpa_b),
    ]
    for y, x_start, cpa in t2_reads:
        out.append(read_bar(x_start, cpa, y, 6, PAL["orange"], PAL["orange_l"]))
        out.append(end_marker(cpa, y, 6, PAL["orange"]))
        # Dashed left end indicates truncation (no 5' anchor)
        out.append(f'<line x1="{x_start}" x2="{x_start}" y1="{y - 5}" y2="{y + 5}" '
                   f'stroke="{PAL["orange"]}" stroke-width="1" stroke-dasharray="2,2"/>')

    # ── T1↔T2 pair linkage at matching CPA ───────────────────────────────────
    # Draw a curved bracket on the right side showing "T1 + T2 pooled at CPA A/B"
    y_pair_a = (t1_reads[1][0] + t2_reads[2][0]) / 2
    y_pair_b = (t1_reads[3][0] + t2_reads[4][0]) / 2

    # Linkage indicators
    out.append(f'<path d="M{cpa_a + 22},{t1_reads[0][0] - 2} '
               f'C{cpa_a + 38},{y_pair_a} {cpa_a + 38},{y_pair_a} '
               f'{cpa_a + 22},{t2_reads[2][0] + 2}" '
               f'fill="none" stroke="{PAL["heading"]}" stroke-width="1.2" '
               f'stroke-dasharray="3,3"/>')
    out.append(text(cpa_a + 46, y_pair_a + 4, "T1 + T2  ·  CPA A",
                    size=8.5, color=PAL["heading"], weight="600"))

    out.append(f'<path d="M{cpa_b + 22},{t1_reads[2][0] - 2} '
               f'C{cpa_b + 38},{y_pair_b} {cpa_b + 38},{y_pair_b} '
               f'{cpa_b + 22},{t2_reads[4][0] + 2}" '
               f'fill="none" stroke="{PAL["heading"]}" stroke-width="1.2" '
               f'stroke-dasharray="3,3"/>')
    out.append(text(cpa_b + 46, y_pair_b + 4, "T1 + T2  ·  CPA B",
                    size=8.5, color=PAL["heading"], weight="600"))

    # ── Footer note ──────────────────────────────────────────────────────────
    out.append(text(FIG_W/2, 285,
                    "linkage rule: same gene + same orientation + |Δ3′| ≤ 5 bp",
                    size=9, color=PAL["muted"], anchor="middle"))

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
    render(build(), "cdna_isoform_clustering")
