#!/usr/bin/env python3
"""Option C — hybrid renderer.

All three rows (mapped, corrected, pA-rest) show the SAME query bases at
the SAME ref positions (taken from the corrected BAM's CIGAR — which is
the canonical post-correction alignment). What differs per row is the
**event overlay**:

  * Insertions: purple ↓ above the row with "+N (XYZ)" annotation
  * Deletions:  orange bracket above the row with "ΔN" annotation
  * pA-rest only: faded soft-clip extension at the 5' end (minus-strand
    poly-A) or 3' end (plus-strand poly-A)

The bases in every row are visually identical → reviewer sees at a glance
that "the read didn't change." Per-row event overlays show "this is what
each CIGAR thought happened."

One-off renderer for cat1_minus_1; if the user approves the design, this
moves into render_read_alignment.py and replaces option A.
"""
from __future__ import annotations

from pathlib import Path
import sys

import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

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
INS_COLOR = "#9c27b0"
DEL_COLOR = "#ff6f00"


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
    """Returns positional state for a read in the window. Captures:
        - bases: read base at each ref position (None if not covered or D/N)
        - state: 'aligned' / 'deletion' / 'intron' / 'softclip5' / 'softclip3'
        - mismatches: set of indices where base != ref
        - insertions: list of (ref_pos_before_insertion, length, bases)
        - deletions: list of (ref_pos_start, length)  (deletions vs read)
        - five_sc, three_sc: soft-clip lengths
        - aln_start, aln_end
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
            if win_start <= ref_pos < win_end and 0 <= i_clip < len(seq):
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
                if win_start <= ref_pos < win_end and 0 <= q_idx < len(seq):
                    idx = ref_pos - win_start
                    out["bases"][idx] = seq[q_idx]
                    out["state"][idx] = "softclip3"
            qp += length
            continue
        if op in (0, 7, 8):
            for _ in range(length):
                if win_start <= rp < win_end and 0 <= qp < len(seq):
                    idx = rp - win_start
                    out["bases"][idx] = seq[qp]
                    out["state"][idx] = "aligned"
                qp += 1
                rp += 1
        elif op == 1:
            ins_bases = seq[qp:qp + length] if qp < len(seq) else ""
            # store ref position AT which the insertion happens (between rp-1 and rp)
            out["insertions"].append((rp, length, ins_bases))
            qp += length
        elif op == 2:
            out["deletions"].append((rp, length))
            for _ in range(length):
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
    if st == "intron":
        return "#E0E0E0", BASE_COLOR["|"]
    if mm:
        return MISMATCH_BG, BASE_COLOR.get(b.upper(), "#000")
    return MATCH_BG, BASE_COLOR.get(b.upper(), "#000")


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def draw_bases(ax, bases, states, mismatches, ws, we):
    n = we - ws
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    for i, (b, st) in enumerate(zip(bases, states)):
        if b is None:
            continue
        bg, fg = _color(b, st, i in mismatches)
        ax.add_patch(Rectangle((i, 0.05), 1, 0.9, facecolor=bg, edgecolor="none"))
        weight = "normal" if st in ("softclip5", "softclip3") else "bold"
        ax.text(i + 0.5, 0.5, b.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight=weight)


def draw_insertion_overlay(ax, insertions, ws, we):
    """Purple ↓ above the row at each insertion site, with "+N (XYZ)" label."""
    n = we - ws
    for rp, ln, bases in insertions:
        ix = rp - ws  # insertion sits BETWEEN positions rp-1 and rp; visual is left edge of column rp
        if 0 <= ix <= n:
            ax.plot(ix, 1.0, marker="v", markersize=6, color=INS_COLOR, clip_on=False)
            label = f"+{ln} ({bases})" if 0 < ln <= 4 and bases else f"+{ln}"
            ax.text(ix, 1.15, label, ha="center", va="bottom",
                    fontsize=6.5, color=INS_COLOR, clip_on=False)


def draw_deletion_overlay(ax, deletions, ws, we):
    """Orange bracket below the row spanning the deleted ref positions, with "ΔN" label."""
    n = we - ws
    for rp_start, ln in deletions:
        x0 = max(0, rp_start - ws)
        x1 = min(n, rp_start + ln - ws)
        if x1 > x0:
            # Bracket: horizontal bar at y=-0.1 to -0.2 from top of row, with end ticks
            y_bar = -0.05
            ax.plot([x0, x1], [y_bar, y_bar], color=DEL_COLOR, linewidth=2.0, clip_on=False)
            ax.plot([x0, x0], [y_bar - 0.05, y_bar + 0.05],
                    color=DEL_COLOR, linewidth=2.0, clip_on=False)
            ax.plot([x1, x1], [y_bar - 0.05, y_bar + 0.05],
                    color=DEL_COLOR, linewidth=2.0, clip_on=False)
            mid = (x0 + x1) / 2
            ax.text(mid, y_bar - 0.20, f"Δ{ln}", ha="center", va="top",
                    fontsize=7, color=DEL_COLOR, clip_on=False, fontweight="bold")


def draw_ref(ax, ref, ws, we, orig, corr):
    n = we - ws
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel("ref", rotation=0, ha="right", va="center",
                  fontsize=9, fontweight="bold")
    for s in ax.spines.values():
        s.set_visible(False)
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


def draw_row(ax, bases, states, mismatches, insertions, deletions, ws, we, label):
    draw_bases(ax, bases, states, mismatches, ws, we)
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=9)
    draw_insertion_overlay(ax, insertions, ws, we)
    draw_deletion_overlay(ax, deletions, ws, we)


# ---------------------------------------------------------------------------
# Compose Option C
# ---------------------------------------------------------------------------

def main():
    ref = load_ref(CHROM, WIN_START, WIN_END)
    mapped = walk_cigar(fetch_read(MAPPED_BAM, QNAME), WIN_START, WIN_END)
    corrected = walk_cigar(fetch_read(CORR_BAM, QNAME), WIN_START, WIN_END)
    parest = walk_cigar(fetch_read(PAREST_BAM, QNAME), WIN_START, WIN_END)
    for a in (mapped, corrected, parest):
        compute_mismatches(a, ref)

    # Canonical bases/states/mismatches = corrected's
    canon_bases = corrected["bases"]
    canon_states = corrected["state"]
    canon_mm = corrected["mismatches"]

    # pA-rest's softclip5 bases extend leftward; merge them onto corrected's bases
    parest_bases = list(canon_bases)
    parest_states = list(canon_states)
    for i, (b, st) in enumerate(zip(parest["bases"], parest["state"])):
        if st == "softclip5" and parest_bases[i] is None:
            parest_bases[i] = b
            parest_states[i] = "softclip5"
    # Same idea for 3' softclip (which doesn't apply for cat1_minus_1 but
    # is needed for plus-strand reads).
    for i, (b, st) in enumerate(zip(parest["bases"], parest["state"])):
        if st == "softclip3" and parest_bases[i] is None:
            parest_bases[i] = b
            parest_states[i] = "softclip3"

    # Figure with stacked rows; all sharing the same x-axis range
    fig = plt.figure(figsize=(14, 5.5))
    fig.suptitle(
        f"cat1_minus_1 · Option C (hybrid) — chrII:{WIN_START}-{WIN_END}  "
        f"(orig {ORIG_3P} → corr {CORR_3P}, −strand)",
        fontsize=12, y=0.99, fontweight="bold",
    )

    # Layout: ref + 3 read rows + ticks
    # All rows have the same width; left margin reserved for labels
    left, right = 0.10, 0.99
    width = right - left
    height_ratios = [1.2, 1.0, 1.0, 1.0, 0.5]  # ref / mapped / corrected / pA-rest / ticks
    panel_pad = 0.025  # vertical padding between rows
    total_h = sum(height_ratios)
    avail_v = 0.83 - 0.04
    h_unit = (avail_v - panel_pad * (len(height_ratios) - 1)) / total_h
    y_cursor = 0.91

    # Ref
    ax = fig.add_axes([left, y_cursor - height_ratios[0] * h_unit, width, height_ratios[0] * h_unit])
    draw_ref(ax, ref, WIN_START, WIN_END, ORIG_3P, CORR_3P)
    y_cursor -= height_ratios[0] * h_unit + panel_pad

    # Mapped row — bases from corrected, events from mapped
    ax = fig.add_axes([left, y_cursor - height_ratios[1] * h_unit, width, height_ratios[1] * h_unit])
    draw_row(ax, canon_bases, canon_states, canon_mm,
             mapped["insertions"], mapped["deletions"],
             WIN_START, WIN_END, "mapped")
    y_cursor -= height_ratios[1] * h_unit + panel_pad

    # Corrected row — bases + events from corrected
    ax = fig.add_axes([left, y_cursor - height_ratios[2] * h_unit, width, height_ratios[2] * h_unit])
    draw_row(ax, canon_bases, canon_states, canon_mm,
             corrected["insertions"], corrected["deletions"],
             WIN_START, WIN_END, "corrected")
    y_cursor -= height_ratios[2] * h_unit + panel_pad

    # pA-rest row — bases include restored soft-clip, events from pA-rest
    ax = fig.add_axes([left, y_cursor - height_ratios[3] * h_unit, width, height_ratios[3] * h_unit])
    draw_row(ax, parest_bases, parest_states, canon_mm,  # mismatches unchanged
             parest["insertions"], parest["deletions"],
             WIN_START, WIN_END, "pA-rest")
    y_cursor -= height_ratios[3] * h_unit + panel_pad

    # Ticks
    ax = fig.add_axes([left, y_cursor - height_ratios[4] * h_unit, width, height_ratios[4] * h_unit])
    draw_ticks(ax, WIN_START, WIN_END)
    y_cursor -= height_ratios[4] * h_unit

    # Footer note explaining the convention
    fig.text(
        0.5, 0.01,
        "Bases in all 3 rows are taken from the corrected CIGAR (read content is identical). "
        "Each row's overlay shows that BAM's own indel events: purple ↓ above = insertion (+N (bases)); "
        "orange bracket below = deletion (ΔN). pA-rest also extends the 5'/3' soft-clip from the restored poly-A tail.",
        ha="center", va="bottom", fontsize=8, color="#666", style="italic",
    )

    out = Path("/tmp/igv_snapshots/cat1_minus_1_option_c.png")
    fig.savefig(str(out), dpi=150, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
