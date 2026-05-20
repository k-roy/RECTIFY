#!/usr/bin/env python3
"""Option A v2 — alignment-faithful rendering with improved insertion display.

Each row uses its own BAM's CIGAR placement (bases at ref positions per
that BAM's alignment). Same design as the original Option A, but with:

  * Inline color-coded inserted bases for insertions ≤ 4 bp:
    e.g., an insertion of "AGT" between ref columns 9841 and 9842 renders
    as small "A G T" characters (each in its own base color) just above
    the row, anchored at the inter-column boundary. No purple triangle
    clutter; the bases ARE the marker.
  * Larger insertions get a "+N" callout with the leading bases shown
    inline as well.
  * Multiple insertions within a few bp are vertically staggered so the
    labels don't collide.
  * More vertical breathing room between rows so insertion text doesn't
    bleed into the row above.

One-off renderer for cat1_minus_1 inspection.
"""
from __future__ import annotations

from pathlib import Path
import sys

import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

REPO = Path("/Users/kevinroy/work/rectify")
VAL = REPO / "rectify" / "data" / "validation"

QNAME = "77b392d9-d80a-4f7c-be62-4dbb854857bc"
CHROM = "chrII"
WIN_START = 9784
WIN_END = 9884
ORIG_3P = 9826
CORR_3P = 9834

GENOME = REPO / "rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
MAPPED_BAM = VAL / "aligners/validation_reads.minimap2.bam"
CORR_BAM = VAL / "rectified/rectified_pA_tail_trimmed.bam"
PAREST_BAM = VAL / "rectified/rectified_pA_tail_soft_clipped.bam"

BASE_COLOR = {
    "A": "#2ECC40", "T": "#FF4136", "C": "#0074D9", "G": "#FF851B",
    "N": "#888888", "-": "#B26C00", "|": "#666666",
}
MISMATCH_BG = "#FFD0D0"
MATCH_BG = "#FAFAFA"
SOFTCLIP_BG = "#F0F0F0"
SOFTCLIP_FG = "#999999"
DEL_BG = "#FFE0B2"


def load_ref(chrom, start, end):
    with pysam.FastaFile(str(GENOME)) as fa:
        return fa.fetch(chrom, max(0, start), end).upper()


def fetch_read(bam_path, qname):
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.query_name == qname:
                return r
    return None


def walk_cigar(read, win_start, win_end):
    """Return per-position state plus insertion + deletion lists.

    Insertions are recorded with the ref position they SIT IMMEDIATELY BEFORE
    (i.e., between ref_pos-1 and ref_pos). For minus-strand soft-clip at the
    leading CIGAR position, bases are placed at ref positions to the LEFT of
    aln_start.
    """
    n = win_end - win_start
    out = dict(
        bases=[None] * n, state=[None] * n,
        mismatches=set(), insertions=[], deletions=[],
        aln_start=None, aln_end=None,
        five_sc=0, three_sc=0,
        full_aln_start=None, full_aln_end=None,
    )
    if read is None:
        return out
    seq = read.query_sequence or ""
    cigar = read.cigartuples or []
    out["aln_start"] = read.reference_start
    out["aln_end"] = read.reference_end

    qp = 0
    rp = read.reference_start
    if cigar and cigar[0][0] == 4:
        clen = cigar[0][1]
        out["five_sc"] = clen
        for i_clip in range(clen):
            ref_pos = rp - clen + i_clip
            if win_start <= ref_pos < win_end and i_clip < len(seq):
                idx = ref_pos - win_start
                out["bases"][idx] = seq[i_clip]
                out["state"][idx] = "softclip5"
        qp = clen
    cigar_iter = iter(cigar)
    if cigar and cigar[0][0] == 4:
        next(cigar_iter)
    for op, length in cigar_iter:
        if op == 5:
            continue
        if op == 4:
            out["three_sc"] = length
            for i_clip in range(length):
                ref_pos = read.reference_end + i_clip
                q_idx = qp + i_clip
                if win_start <= ref_pos < win_end and q_idx < len(seq):
                    idx = ref_pos - win_start
                    out["bases"][idx] = seq[q_idx]
                    out["state"][idx] = "softclip3"
            qp += length
            continue
        if op in (0, 7, 8):
            for _ in range(length):
                if win_start <= rp < win_end and qp < len(seq):
                    idx = rp - win_start
                    out["bases"][idx] = seq[qp]
                    out["state"][idx] = "aligned"
                qp += 1
                rp += 1
        elif op == 1:
            ins_bases = seq[qp:qp + length] if qp < len(seq) else ""
            out["insertions"].append((rp, length, ins_bases))
            qp += length
        elif op == 2:
            out["deletions"].append((rp, length))
            for _ in range(length):
                if win_start <= rp < win_end:
                    idx = rp - win_start
                    out["bases"][idx] = "-"
                    out["state"][idx] = "deletion"
                rp += 1
        elif op == 3:
            for _ in range(length):
                if win_start <= rp < win_end:
                    idx = rp - win_start
                    out["bases"][idx] = "|"
                    out["state"][idx] = "intron"
                rp += 1
    out["full_aln_start"] = read.reference_start - out["five_sc"]
    out["full_aln_end"] = read.reference_end + out["three_sc"]
    return out


def compute_mismatches(a, ref):
    ms = set()
    for i, (b, st) in enumerate(zip(a["bases"], a["state"])):
        if b is None or st != "aligned":
            continue
        if b.upper() != ref[i].upper():
            ms.add(i)
    a["mismatches"] = ms


def _color(b, st, mm):
    if st in ("softclip5", "softclip3"):
        return SOFTCLIP_BG, SOFTCLIP_FG
    if st == "deletion":
        return DEL_BG, BASE_COLOR["-"]
    if st == "intron":
        return "#E0E0E0", BASE_COLOR["|"]
    if mm:
        return MISMATCH_BG, BASE_COLOR.get(b.upper(), "#000")
    return MATCH_BG, BASE_COLOR.get(b.upper(), "#000")


def draw_row_with_insertions(ax, a, ws, we, label):
    """Draw a row of bases + inline color-coded insertion characters above."""
    n = we - ws
    ax.set_xlim(0, n)
    ax.set_ylim(-0.05, 2.00)  # extra headroom for staggered insertion pills
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=9)

    # Main row bases
    for i, (b, st) in enumerate(zip(a["bases"], a["state"])):
        if b is None:
            continue
        bg, fg = _color(b, st, i in a["mismatches"])
        ax.add_patch(Rectangle((i, 0.05), 1, 0.9, facecolor=bg, edgecolor="none"))
        weight = "normal" if st in ("softclip5", "softclip3") else "bold"
        ax.text(i + 0.5, 0.5, b.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight=weight)

    # Insertion overlay — render inserted bases as FULL-SIZE chars wedged
    # between the surrounding ref columns. The bases sit on the SAME y as
    # the row's bases (vertically centered on the row) but with their own
    # narrow inter-column "slot". A vivid purple background pill spanning
    # the insertion makes the wedge unambiguous: these are bases in the
    # read that have no ref position (insertion).
    #
    # Slot geometry: each inserted base gets ~0.55 col-widths of horizontal
    # space, centered on the inter-column boundary. The pill extends ~0.10
    # col-widths into each adjacent ref column to visually separate the
    # inserted bases from the surrounding row bases.
    insertions = sorted(a["insertions"], key=lambda t: t[0])
    last_x_end = -100
    stagger = 0
    for rp, ln, bases in insertions:
        ix = rp - ws  # column edge where insertion happens (between rp-1 and rp)
        if ix < 0 or ix > n:
            continue
        # Display geometry
        char_width = 0.55  # roughly half the row's column width per inserted base
        total_w = ln * char_width
        start_x = ix - total_w / 2
        end_x = ix + total_w / 2
        # Stagger if overlapping with the previous insertion's pill
        if start_x < last_x_end + 0.3:
            stagger = 1 - stagger
        else:
            stagger = 0
        # Pill sits ABOVE the row when staggered=0, half-overlapping the row
        # when staggered=1, so adjacent insertions step vertically.
        pill_y = 0.55 if stagger == 0 else 1.15
        pill_h = 0.55
        # Thin vertical guide from row top down to the pill bottom (or up)
        guide_y1 = pill_y - 0.02 if stagger == 1 else pill_y + pill_h + 0.02
        # Render bases first so background pill is BELOW them visually
        ax.add_patch(Rectangle((start_x - 0.10, pill_y),
                               total_w + 0.20, pill_h,
                               facecolor="#E1BEE7", edgecolor="#7B1FA2",
                               linewidth=1.2, clip_on=False, zorder=3))
        if ln <= 4 and bases:
            for j, b in enumerate(bases):
                x = start_x + (j + 0.5) * char_width
                fg = BASE_COLOR.get(b.upper(), "#7B1FA2")
                ax.text(x, pill_y + pill_h / 2, b.upper(),
                        ha="center", va="center",
                        fontsize=8.5, family="monospace",
                        color=fg, fontweight="bold", clip_on=False, zorder=4)
        else:
            label_text = f"+{ln}"
            if bases:
                label_text += f" {bases[:3]}…" if ln > 3 else f" {bases}"
            ax.text(ix, pill_y + pill_h / 2, label_text, ha="center", va="center",
                    fontsize=7.0, color="#7B1FA2", fontweight="bold",
                    clip_on=False, zorder=4)
        last_x_end = end_x


def draw_ref(ax, ref, ws, we, orig, corr):
    n = we - ws
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_ylabel("ref", rotation=0, ha="right", va="center",
                  fontsize=9, fontweight="bold")
    for i, b in enumerate(ref):
        ax.add_patch(Rectangle((i, 0.05), 1, 0.7, facecolor="#F5F5F5", edgecolor="none"))
        ax.text(i + 0.5, 0.4, b.upper(), ha="center", va="center",
                fontsize=8, family="monospace",
                color=BASE_COLOR.get(b.upper(), "#000"), fontweight="bold")
    if ws <= orig < we:
        x = orig - ws + 0.5
        ax.plot(x, 1.05, marker="v", markersize=10, color="#d32f2f", clip_on=False)
        ax.text(x, 1.30, "orig", ha="center", va="bottom",
                fontsize=8, color="#d32f2f", clip_on=False, fontweight="bold")
    if ws <= corr < we:
        x = corr - ws + 0.5
        ax.plot(x, 1.05, marker="v", markersize=10, color="#2e7d32", clip_on=False)
        ax.text(x, 1.30, "corr", ha="center", va="bottom",
                fontsize=8, color="#2e7d32", clip_on=False, fontweight="bold")


def draw_ticks(ax, ws, we):
    n = we - ws
    ax.set_xlim(0, n)
    ax.set_ylim(0, 0.4)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    step = max(10, n // 8)
    for x in range(0, n + 1, step):
        ax.plot([x, x], [0.3, 0.4], color="#666", linewidth=0.7)
        ax.text(x, 0.0, f"{ws + x:,}", ha="center", va="bottom",
                fontsize=7, color="#444")


def main():
    ref = load_ref(CHROM, WIN_START, WIN_END)
    mapped = walk_cigar(fetch_read(MAPPED_BAM, QNAME), WIN_START, WIN_END)
    corrected = walk_cigar(fetch_read(CORR_BAM, QNAME), WIN_START, WIN_END)
    parest = walk_cigar(fetch_read(PAREST_BAM, QNAME), WIN_START, WIN_END)
    for a in (mapped, corrected, parest):
        compute_mismatches(a, ref)

    fig = plt.figure(figsize=(14, 8))
    fig.suptitle(
        f"cat1_minus_1 · Option A v2 — chrII:{WIN_START}-{WIN_END}  "
        f"(orig {ORIG_3P} → corr {CORR_3P}, −strand)",
        fontsize=12, y=0.995, fontweight="bold",
    )

    left, right = 0.10, 0.99
    width = right - left
    # Tall row spacing — each read row has 2.0 units of vertical extent
    # (incl. staggered insertion-pill headroom) so events don't collide
    # with neighbors
    panel_h = 0.16
    pad = 0.03
    y = 0.91

    ax = fig.add_axes([left, y - 0.10, width, 0.10])
    draw_ref(ax, ref, WIN_START, WIN_END, ORIG_3P, CORR_3P)
    y -= 0.10 + pad

    for label, a in [("mapped", mapped), ("corrected", corrected), ("pA-rest", parest)]:
        ax = fig.add_axes([left, y - panel_h, width, panel_h])
        draw_row_with_insertions(ax, a, WIN_START, WIN_END, label)
        y -= panel_h + pad

    ax = fig.add_axes([left, y - 0.04, width, 0.04])
    draw_ticks(ax, WIN_START, WIN_END)

    fig.text(
        0.5, 0.01,
        "Inserted bases shown INLINE between ref columns above each row (colored by base identity). "
        "Deletions show as orange '-' in the row. Soft-clip bases faded grey.",
        ha="center", va="bottom", fontsize=8, color="#666", style="italic",
    )

    out = Path("/tmp/igv_snapshots/cat1_minus_1_option_a_v2.png")
    fig.savefig(str(out), dpi=150, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
