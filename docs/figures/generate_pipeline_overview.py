#!/usr/bin/env python3
"""
Generate the top-of-README workflow overview figure for RECTIFY.

Three protocol columns (DRS, ONT cDNA, QuantSeq REV) converging at the
shared multi-aligner consensus, read-vs-reference walkback, and downstream
analysis stages.

v2 (2026-05-15): redesign for legibility + aesthetic polish.
  - Larger title and clearer typographic hierarchy
  - Left-gutter row labels anchor each pipeline stage
  - Thin colored track rail down each column unifies per-protocol identity
  - Shared-stage rows get a "SHARED" pill in the corner and stronger border
  - QSREV Step 0 reframed as a small "passthrough" chip rather than a full card
  - Step 3 → Step 4 funnel uses real bezier curves converging on a single point
  - Output labels lifted into their own chips below each card (no edge crowding)
  - Body text bumped from 8.5 → 9.5pt for README-width legibility
"""

import os
import cairosvg

OUTDIR = os.path.dirname(os.path.abspath(__file__))

FIG_W = 920
FIG_H = 840
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"

PAL = dict(
    title="#0f172a",
    heading="#334155",
    label="#475569",
    muted="#64748b",
    faint="#94a3b8",
    border="#cbd5e1",
    divider="#e2e8f0",
    bg="#ffffff",
    panel="#f8fafc",
    blue="#2563eb",
    blue_l="#dbeafe",
    blue_ll="#eff6ff",
    green="#059669",
    green_l="#d1fae5",
    teal="#0d9488",
    teal_l="#ccfbf1",
    teal_ll="#f0fdfa",
    orange="#d97706",
    orange_l="#fed7aa",
    orange_ll="#fff7ed",
    indigo="#4f46e5",
    indigo_l="#e0e7ff",
)

TRACK = {
    "drs":   dict(color=PAL["blue"],   fill=PAL["blue_ll"],   accent=PAL["blue_l"]),
    "cdna":  dict(color=PAL["teal"],   fill=PAL["teal_ll"],   accent=PAL["teal_l"]),
    "qsrev": dict(color=PAL["orange"], fill=PAL["orange_ll"], accent=PAL["orange_l"]),
}

# Layout: 80 px gutter on left for row labels, three columns of 240 px each.
GUTTER_W = 88
COL_W = 240
COL_GAP = 28
COL_X = {
    "drs":   GUTTER_W + 0 * (COL_W + COL_GAP),
    "cdna":  GUTTER_W + 1 * (COL_W + COL_GAP),
    "qsrev": GUTTER_W + 2 * (COL_W + COL_GAP),
}
SHARED_X = GUTTER_W
SHARED_W = 3 * COL_W + 2 * COL_GAP

# Row positions (top y of each row, in document order)
ROW = dict(
    title=0,
    inputs=70,
    preprocess=148,
    align=258,
    walkback=362,
    correct=478,
    analyze=608,
    footer=720,
)
ROW_H = dict(
    inputs=66,
    preprocess=86,
    align=82,
    walkback=92,
    correct=104,
    analyze=88,
)


def _xml_escape(s: str) -> str:
    return (s.replace("&", "&amp;")
             .replace("<", "&lt;")
             .replace(">", "&gt;"))


def text(x, y, s, *, size=10, weight="normal", anchor="middle",
        fill=None, italic=False, ls=None, family=None):
    extras = []
    if ls is not None:
        extras.append(f'letter-spacing="{ls}"')
    if italic:
        extras.append('font-style="italic"')
    if family:
        extras.append(f'font-family="{family}"')
    extra = " " + " ".join(extras) if extras else ""
    fill = fill or PAL["heading"]
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}" '
            f'fill="{fill}"{extra}>{_xml_escape(s)}</text>')


def rounded_rect(x, y, w, h, *, rx=8, fill="none", stroke=None,
                 stroke_w=1.2, opacity=None, dash=None):
    attrs = [
        f'x="{x}"', f'y="{y}"',
        f'width="{w}"', f'height="{h}"', f'rx="{rx}"',
        f'fill="{fill}"',
    ]
    if stroke:
        attrs.append(f'stroke="{stroke}"')
        attrs.append(f'stroke-width="{stroke_w}"')
    if dash:
        attrs.append(f'stroke-dasharray="{dash}"')
    if opacity is not None:
        attrs.append(f'opacity="{opacity}"')
    return f'<rect {" ".join(attrs)}/>'


def track_card(track_key, x, y, w, h, title, body_lines, *,
               step_chip=None, output_chip=None):
    """A per-track card with rounded fill, accent left rail, title, body, optional output chip."""
    t = TRACK[track_key]
    parts = []

    # Card background (very pale track tint)
    parts.append(rounded_rect(x, y, w, h, rx=10, fill=t["fill"]))
    # Card border
    parts.append(rounded_rect(x, y, w, h, rx=10,
                              stroke=t["color"], stroke_w=1.3))

    # Left rail (thin solid bar in track color)
    parts.append(
        f'<rect x="{x}" y="{y}" width="3" height="{h}" '
        f'fill="{t["color"]}"/>'
    )

    cx = x + w / 2
    cy_title = y + 24

    # Optional step chip in upper-left
    if step_chip:
        chip_w = 56
        chip_x = x + 14
        chip_y = y + 12
        parts.append(rounded_rect(
            chip_x, chip_y, chip_w, 16, rx=8,
            fill=t["accent"], stroke=t["color"], stroke_w=0.8))
        parts.append(text(chip_x + chip_w / 2, chip_y + 11.5,
                          step_chip, size=8.5, weight="700",
                          fill=t["color"], ls="0.4"))
        cy_title = y + 42

    # Title
    parts.append(text(cx, cy_title, title, size=11.5, weight="700",
                      fill=t["color"]))

    # Body lines (centered)
    cy = cy_title + 18
    for line in body_lines:
        parts.append(text(cx, cy, line, size=9.5, fill=PAL["heading"]))
        cy += 13

    # Output chip near bottom
    if output_chip:
        out_w = w - 32
        out_x = x + 16
        out_y = y + h - 26
        parts.append(rounded_rect(
            out_x, out_y, out_w, 18, rx=4,
            fill="#ffffff", stroke=t["accent"], stroke_w=1))
        parts.append(text(out_x + out_w / 2, out_y + 13,
                          output_chip, size=9, weight="600",
                          fill=t["color"], italic=True,
                          family="ui-monospace, Menlo, Consolas, monospace"))

    return "\n".join(parts)


def shared_row(x, y, w, h, accent, step_label, title, body_lines, *,
               pill_text="SHARED"):
    """Full-width shared-stage row with stronger visual identity."""
    parts = []
    # Soft panel background
    parts.append(rounded_rect(x, y, w, h, rx=12,
                              fill=PAL["panel"]))
    # Border in accent color (solid this time, not dashed)
    parts.append(rounded_rect(x, y, w, h, rx=12,
                              stroke=accent, stroke_w=1.5))

    # Left accent bar
    parts.append(
        f'<rect x="{x}" y="{y}" width="4" height="{h}" '
        f'fill="{accent}"/>'
    )

    # Step chip (upper-left, inside)
    chip_x = x + 18
    chip_y = y + 12
    chip_w = 64
    parts.append(rounded_rect(
        chip_x, chip_y, chip_w, 18, rx=9,
        fill=accent, stroke="none"))
    parts.append(text(chip_x + chip_w / 2, chip_y + 13,
                      step_label, size=9, weight="700",
                      fill="#ffffff", ls="0.6"))

    # SHARED pill (upper-right)
    pill_w = 60
    pill_x = x + w - pill_w - 14
    pill_y = y + 12
    parts.append(rounded_rect(
        pill_x, pill_y, pill_w, 18, rx=9,
        fill="#ffffff", stroke=accent, stroke_w=1))
    parts.append(text(pill_x + pill_w / 2, pill_y + 13,
                      pill_text, size=8.5, weight="700",
                      fill=accent, ls="0.6"))

    # Title
    parts.append(text(x + 96, y + 25, title, size=11.5, weight="700",
                      fill=accent, anchor="start"))

    # Body lines (left-aligned, indented under the title)
    cy = y + 50
    for line in body_lines:
        parts.append(text(x + 96, cy, line, size=10, fill=PAL["heading"],
                          anchor="start"))
        cy += 15

    return "\n".join(parts)


def passthrough_chip(track_key, x, y, w, h, label, subline):
    """Small dashed 'skip' chip for tracks that bypass a stage."""
    t = TRACK[track_key]
    parts = []
    # Center the chip in the allocated card slot
    chip_w = w - 60
    chip_h = 40
    chip_x = x + (w - chip_w) / 2
    chip_y = y + (h - chip_h) / 2
    parts.append(rounded_rect(chip_x, chip_y, chip_w, chip_h, rx=10,
                              fill="#ffffff",
                              stroke=t["color"], stroke_w=1.2,
                              dash="4,3"))
    parts.append(text(chip_x + chip_w / 2, chip_y + 17,
                      label, size=10, weight="700",
                      fill=t["color"], ls="0.4"))
    parts.append(text(chip_x + chip_w / 2, chip_y + 31,
                      subline, size=8.5, italic=True,
                      fill=PAL["muted"]))
    return "\n".join(parts)


def vert_arrow(x, y1, y2, color, *, w=1.5):
    return (
        f'<line x1="{x}" x2="{x}" y1="{y1}" y2="{y2 - 5}" '
        f'stroke="{color}" stroke-width="{w}"/>'
        f'<polygon points="{x},{y2} {x - 4},{y2 - 6} {x + 4},{y2 - 6}" '
        f'fill="{color}"/>'
    )


def bezier_funnel(x_from, y_from, x_to, y_to, color, *, w=1.5):
    """Bezier curve converging from x_from to x_to."""
    cy = (y_from + y_to) / 2
    d = (f"M {x_from} {y_from} "
         f"C {x_from} {cy}, {x_to} {cy}, {x_to} {y_to - 5}")
    return (
        f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{w}"/>'
        f'<polygon points="{x_to},{y_to} {x_to - 4},{y_to - 6} '
        f'{x_to + 4},{y_to - 6}" fill="{color}"/>'
    )


def gutter_label(y, text_top, text_bottom=None):
    """Stage label in the left gutter."""
    parts = []
    parts.append(text(
        GUTTER_W - 18, y, text_top, size=10.5, weight="700",
        anchor="end", fill=PAL["muted"], ls="0.5"))
    if text_bottom:
        parts.append(text(
            GUTTER_W - 18, y + 14, text_bottom, size=9,
            anchor="end", fill=PAL["faint"]))
    return "\n".join(parts)


def main():
    svg = []
    svg.append('<?xml version="1.0" encoding="utf-8" ?>')
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" '
        f'height="{FIG_H}" viewBox="0 0 {FIG_W} {FIG_H}">'
    )
    svg.append(f'<rect fill="{PAL["bg"]}" width="{FIG_W}" height="{FIG_H}"/>')
    svg.append(f'<g font-family="{FONT}">')

    # ── Title block ────────────────────────────────────────────────────────
    svg.append(text(FIG_W / 2, 32, "RECTIFY pipeline overview",
                    size=20, weight="700", fill=PAL["title"]))
    svg.append(text(FIG_W / 2, 52,
                    "Three protocol tracks share a multi-aligner / walkback / analyze core",
                    size=11, italic=True, fill=PAL["muted"]))

    # ── Row 1: input metadata cards ────────────────────────────────────────
    y = ROW["inputs"]; h = ROW_H["inputs"]
    svg.append(gutter_label(y + 26, "INPUTS"))
    svg.append(track_card(
        "drs", COL_X["drs"], y, COL_W, h,
        "ONT direct RNA (DRS)",
        ["poly(A) tail in read • 3′ = right side",
         "RNA strand = BAM strand"],
    ))
    svg.append(track_card(
        "cdna", COL_X["cdna"], y, COL_W, h,
        "ONT PCR-cDNA (PCB114.24)",
        ["SSP+UMI bridge • poly(A) in read",
         "RNA strand = BAM strand"],
    ))
    svg.append(track_card(
        "qsrev", COL_X["qsrev"], y, COL_W, h,
        "QuantSeq REV (short read)",
        ["dT-primed antisense • no poly(A) in read",
         "RNA strand opposite BAM strand"],
    ))

    # connectors row1 → row2
    for key in ("drs", "cdna", "qsrev"):
        svg.append(vert_arrow(COL_X[key] + COL_W / 2,
                              y + h + 2, ROW["preprocess"] - 3,
                              TRACK[key]["color"]))

    # ── Row 2: per-track preprocess (Step 0) ───────────────────────────────
    y = ROW["preprocess"]; h = ROW_H["preprocess"]
    svg.append(gutter_label(y + 18, "STEP 0", "preprocess"))
    svg.append(track_card(
        "drs", COL_X["drs"], y, COL_W, h,
        "rectify trim-polya",
        ["3-pass poly(A) + adapter strip",
         "tail length cached for restore"],
        step_chip="STEP 0",
        output_chip="trimmed.bam",
    ))
    svg.append(track_card(
        "cdna", COL_X["cdna"], y, COL_W, h,
        "rectify correct-cdna",
        ["UMI cluster · abPOA consensus",
         "strip SSP+UMI+GGG (5′) and poly(A) (3′)"],
        step_chip="STEP 0",
        output_chip="stage1_consensus.fastq.gz",
    ))
    svg.append(passthrough_chip(
        "qsrev", COL_X["qsrev"], y, COL_W, h,
        "skip — no preprocess needed",
        "QSREV reads carry no poly(A)"))

    # connectors row2 → row3 (DRS, cDNA from card bottom; QSREV from passthrough bottom)
    align_top = ROW["align"]
    for key in ("drs", "cdna"):
        svg.append(vert_arrow(COL_X[key] + COL_W / 2,
                              y + h + 2, align_top - 3,
                              TRACK[key]["color"]))
    # QSREV: arrow from middle of passthrough chip area
    svg.append(vert_arrow(COL_X["qsrev"] + COL_W / 2,
                          y + h + 2, align_top - 3,
                          TRACK["qsrev"]["color"]))

    # ── Row 3: shared multi-aligner ────────────────────────────────────────
    y = ROW["align"]; h = ROW_H["align"]
    svg.append(gutter_label(y + 18, "STEP 1", "align"))
    svg.append(shared_row(
        SHARED_X, y, SHARED_W, h,
        PAL["indigo"], "STEP 1",
        "rectify align — multi-aligner consensus",
        ["long-read: minimap2 + mapPacBio + gapmm2  (optional tier-2: deSALT, uLTRA)",
         "short-read: bbmap + bwa     |     chimeric_consensus.select_best per read or cluster"],
    ))

    # connectors row3 → row4
    for key in ("drs", "cdna", "qsrev"):
        svg.append(vert_arrow(COL_X[key] + COL_W / 2,
                              y + h + 2, ROW["walkback"] - 3,
                              PAL["indigo"], w=1.3))

    # ── Row 4: shared walkback ─────────────────────────────────────────────
    y = ROW["walkback"]; h = ROW_H["walkback"]
    svg.append(gutter_label(y + 18, "STEP 2", "walkback"))
    svg.append(shared_row(
        SHARED_X, y, SHARED_W, h,
        PAL["green"], "STEP 2",
        "Read-vs-reference walkback   (rectify/core/correct/walkback.py)",
        ["3-case terminal gate:  (1) non-A anchored → no correction   "
         "(2) read-A & genome-A → walk inward   (3) mismatch → walk to first non-A read=ref",
         "walkback_drs (3′=right, RNA=BAM)   ·   cdna-analyze recompute   ·   "
         "walkback_quantseq_rev (3′=left, RNA inverted)"],
        pill_text="CORE + 3",
    ))

    # connectors row4 → row5 (back to track colors)
    for key in ("drs", "cdna", "qsrev"):
        svg.append(vert_arrow(COL_X[key] + COL_W / 2,
                              y + h + 2, ROW["correct"] - 3,
                              TRACK[key]["color"]))

    # ── Row 5: per-track correct (Step 3) ──────────────────────────────────
    y = ROW["correct"]; h = ROW_H["correct"]
    svg.append(gutter_label(y + 18, "STEP 3", "correct / analyze"))
    svg.append(track_card(
        "drs", COL_X["drs"], y, COL_W, h,
        "rectify correct",
        ["per-read CPA correction",
         "optional restore-softclip step",
         "(poly(A) re-attached for IGV)"],
        step_chip="STEP 3",
        output_chip="corrected_reads.tsv",
    ))
    svg.append(track_card(
        "cdna", COL_X["cdna"], y, COL_W, h,
        "rectify cdna-analyze",
        ["walkback on post-align coords",
         "isoform clustering (XI)",
         "Type-1 ↔ Type-2 pairing (XL)"],
        step_chip="STEP 3",
        output_chip="clusters · isoforms · pairs.tsv",
    ))
    svg.append(track_card(
        "qsrev", COL_X["qsrev"], y, COL_W, h,
        "rectify correct --dT-primed-cDNA",
        ["per-read CPA correction",
         "strand inverted to RNA",
         "internal-priming caught by walkback"],
        step_chip="STEP 3",
        output_chip="corrected_reads.tsv",
    ))

    # ── Row 5 → Row 6 funnel ───────────────────────────────────────────────
    analyze_y = ROW["analyze"]
    analyze_target_x = SHARED_X + SHARED_W / 2
    for key in ("drs", "cdna", "qsrev"):
        col_mid = COL_X[key] + COL_W / 2
        # All three funnel into the same column on the top of the analyze row,
        # spaced by ±28 px so the arrowheads don't overlap.
        if key == "drs":
            target = analyze_target_x - 28
        elif key == "cdna":
            target = analyze_target_x
        else:
            target = analyze_target_x + 28
        svg.append(bezier_funnel(col_mid, y + h + 2,
                                 target, analyze_y - 3,
                                 TRACK[key]["color"], w=1.4))

    # ── Row 6: shared analyze ──────────────────────────────────────────────
    y = ROW["analyze"]; h = ROW_H["analyze"]
    svg.append(gutter_label(y + 18, "STEP 4", "downstream"))
    svg.append(shared_row(
        SHARED_X, y, SHARED_W, h,
        PAL["blue"], "STEP 4",
        "rectify analyze — shared downstream",
        ["adaptive CPA clustering   ·   stranded bedGraph export   ·   "
         "DESeq2 (gene + cluster resolution)",
         "APA shift detection   ·   GO enrichment   ·   "
         "de novo motif discovery   ·   HTML report"],
    ))

    # ── Footer ─────────────────────────────────────────────────────────────
    foot_y = ROW["footer"]
    svg.append(
        f'<line x1="{SHARED_X}" x2="{SHARED_X + SHARED_W}" '
        f'y1="{foot_y}" y2="{foot_y}" stroke="{PAL["divider"]}" '
        f'stroke-width="1"/>'
    )
    svg.append(text(FIG_W / 2, foot_y + 24,
                    "Shared stages reuse one codebase. Per-track stages encode the "
                    "chemistry-specific differences",
                    size=10, italic=True, fill=PAL["muted"]))
    svg.append(text(FIG_W / 2, foot_y + 40,
                    "(strand convention, poly(A) presence, UMI architecture).",
                    size=10, italic=True, fill=PAL["muted"]))

    # Track legend (track colors only — stage colors are obvious from context)
    leg_y = foot_y + 70
    leg_items = [
        ("ONT direct RNA (DRS)",     PAL["blue"]),
        ("ONT PCR-cDNA",             PAL["teal"]),
        ("QuantSeq REV",             PAL["orange"]),
    ]
    item_w = 160
    total_w = item_w * len(leg_items)
    leg_x = (FIG_W - total_w) / 2
    for i, (lbl, col) in enumerate(leg_items):
        cx = leg_x + i * item_w
        svg.append(f'<rect x="{cx}" y="{leg_y - 9}" width="14" height="12" '
                   f'rx="2" fill="{col}"/>')
        svg.append(text(cx + 22, leg_y, lbl, size=10,
                        fill=PAL["heading"], anchor="start"))

    svg.append('</g>')
    svg.append('</svg>')

    svg_path = os.path.join(OUTDIR, "pipeline_overview.svg")
    png_path = os.path.join(OUTDIR, "pipeline_overview.png")
    with open(svg_path, "w") as f:
        f.write("\n".join(svg))

    cairosvg.svg2png(
        url=svg_path,
        write_to=png_path,
        output_width=FIG_W * 2,
        output_height=FIG_H * 2,
    )
    print(f"wrote {svg_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
