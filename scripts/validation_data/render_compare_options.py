#!/usr/bin/env python3
"""One-off comparison: render cat1_minus_1 in two visualization modes and
stack them into a single figure for user inspection.

  Option A (top): "alignment-faithful" — each row shows its BAM's own
  CIGAR placement. Same bases can appear at different ref positions across
  rows when CIGARs diverge (e.g., mapped's simple M ops vs corrected's
  =/X resolution). This is what we have today.

  Option B (bottom): "consistency-faithful" — anchor all rows to the
  corrected BAM's query→ref mapping so identical query bases line up
  vertically across rows. Mapped's distinct CIGAR is shown only via
  annotation diffs (e.g., "mapped had +3I here") above the row.
  pA-rest is rendered as corrected + the prepended poly-A as a
  faded prefix.

Also bakes in the three small fixes:
  - Short insertions (length ≤ 4) show their bases: "+2 (AC)"
  - Label margin widened so "overview"/"mapped"/etc don't collide
  - "orig"/"corr" arrow text moved above the triangle and clipped properly
"""
from __future__ import annotations

import sys
from pathlib import Path

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
MAPPED = VAL / "aligners/validation_reads.minimap2.bam"
CORR = VAL / "rectified/rectified_pA_tail_trimmed.bam"
PAREST = VAL / "rectified/rectified_pA_tail_soft_clipped.bam"
BG_PLUS = VAL / "rectified/corrected_3ends.plus.bedgraph"
BG_MINUS = VAL / "rectified/corrected_3ends.minus.bedgraph"

BASE_COLOR = {
    "A": "#2ECC40", "T": "#FF4136", "C": "#0074D9", "G": "#FF851B",
    "N": "#888", "-": "#B26C00", "|": "#666",
}
MISMATCH_BG = "#FFD0D0"
MATCH_BG = "#FAFAFA"
SOFTCLIP_BG = "#F0F0F0"
SOFTCLIP_FG = "#999"
DEL_BG = "#FFE0B2"
INTRON_BG = "#E0E0E0"


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


def walk_cigar(read, win_start, win_end, *, include_5sc=True, include_3sc=True):
    """Walk CIGAR. Returns per-position state with insertion bases captured."""
    n = win_end - win_start
    bases = [None] * n
    state = [None] * n
    insertions = []  # list of (ref_pos_before_ins, length, inserted_bases)
    out = dict(
        bases=bases, state=state, insertions=insertions, mismatches=set(),
        aln_start=None, aln_end=None,
        full_aln_start=None, full_aln_end=None,
        five_sc=0, three_sc=0,
        cigar_summary="",
    )
    if read is None:
        return out
    seq = read.query_sequence or ""
    cigar = read.cigartuples or []
    out["cigar_summary"] = read.cigarstring or ""
    out["aln_start"] = read.reference_start
    out["aln_end"] = read.reference_end
    qp = 0
    rp = read.reference_start
    if cigar and cigar[0][0] == 4:
        clen = cigar[0][1]
        out["five_sc"] = clen
        if include_5sc:
            for i_clip in range(clen):
                ref_pos = rp - clen + i_clip
                q_idx = i_clip
                if win_start <= ref_pos < win_end and 0 <= q_idx < len(seq):
                    idx = ref_pos - win_start
                    bases[idx] = seq[q_idx]
                    state[idx] = "softclip5"
        qp = clen
    cigar_iter = iter(cigar)
    if cigar and cigar[0][0] == 4:
        next(cigar_iter)
    for op, length in cigar_iter:
        if op == 5:
            continue
        elif op == 4:
            out["three_sc"] = length
            if include_3sc:
                for i_clip in range(length):
                    ref_pos = read.reference_end + i_clip
                    q_idx = qp + i_clip
                    if win_start <= ref_pos < win_end and 0 <= q_idx < len(seq):
                        idx = ref_pos - win_start
                        bases[idx] = seq[q_idx]
                        state[idx] = "softclip3"
            qp += length
            continue
        elif op in (0, 7, 8):
            for _ in range(length):
                if win_start <= rp < win_end and 0 <= qp < len(seq):
                    idx = rp - win_start
                    bases[idx] = seq[qp]
                    state[idx] = "aligned"
                qp += 1
                rp += 1
        elif op == 1:
            ins_bases = seq[qp:qp + length] if 0 <= qp < len(seq) else ""
            insertions.append((rp - 1, length, ins_bases))
            qp += length
        elif op == 2:
            for _ in range(length):
                if win_start <= rp < win_end:
                    idx = rp - win_start
                    bases[idx] = "-"
                    state[idx] = "deletion"
                rp += 1
        elif op == 3:
            for _ in range(length):
                if win_start <= rp < win_end:
                    idx = rp - win_start
                    bases[idx] = "|"
                    state[idx] = "intron"
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
        return INTRON_BG, BASE_COLOR["|"]
    if mm:
        return MISMATCH_BG, BASE_COLOR.get(b.upper(), "#000")
    return MATCH_BG, BASE_COLOR.get(b.upper(), "#000")


def _draw_row(ax, a, ref, ws, we, label, *, show_ins=True):
    n = we - ws
    ax.set_xlim(0, n); ax.set_ylim(0, 1)
    ax.set_yticks([]); ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=9)
    for s in ax.spines.values(): s.set_visible(False)
    ms = a["mismatches"]
    for i, (b, st) in enumerate(zip(a["bases"], a["state"])):
        if b is None: continue
        bg, fg = _color(b, st, i in ms)
        ax.add_patch(Rectangle((i, 0.05), 1, 0.9, facecolor=bg, edgecolor="none"))
        weight = "normal" if st in ("softclip5", "softclip3") else "bold"
        ax.text(i + 0.5, 0.5, b.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight=weight)
    if show_ins:
        for rp, ln, bases in a["insertions"]:
            ix = rp - ws + 1
            if 0 <= ix <= n:
                ax.plot(ix, 1.0, marker="v", markersize=6, color="#9c27b0", clip_on=False)
                label_text = f"+{ln} ({bases})" if 0 < ln <= 4 and bases else f"+{ln}"
                ax.text(ix, 1.12, label_text, ha="center", va="bottom",
                        fontsize=6.0, color="#9c27b0", clip_on=False)


def _draw_ref(ax, ref, ws, we, orig, corr):
    n = we - ws
    ax.set_xlim(0, n); ax.set_ylim(0, 1.5); ax.set_yticks([]); ax.set_xticks([])
    ax.set_ylabel("ref", rotation=0, ha="right", va="center",
                  fontsize=9, fontweight="bold")
    for s in ax.spines.values(): s.set_visible(False)
    for i, b in enumerate(ref):
        bg = "#F5F5F5"; fg = BASE_COLOR.get(b.upper(), "#000")
        ax.add_patch(Rectangle((i, 0.05), 1, 0.7, facecolor=bg, edgecolor="none"))
        ax.text(i + 0.5, 0.4, b.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight="bold")
    if ws <= orig < we:
        x = (orig - ws) + 0.5
        ax.plot(x, 1.05, marker="v", markersize=9, color="#d32f2f", clip_on=False)
        ax.text(x, 1.30, "orig", ha="center", va="bottom",
                fontsize=8, color="#d32f2f", clip_on=False, fontweight="bold")
    if ws <= corr < we:
        x = (corr - ws) + 0.5
        ax.plot(x, 1.05, marker="v", markersize=9, color="#2e7d32", clip_on=False)
        ax.text(x, 1.30, "corr", ha="center", va="bottom",
                fontsize=8, color="#2e7d32", clip_on=False, fontweight="bold")


def _draw_ticks(ax, ws, we):
    n = we - ws
    ax.set_xlim(0, n); ax.set_ylim(0, 0.4); ax.set_yticks([]); ax.set_xticks([])
    for s in ax.spines.values(): s.set_visible(False)
    step = max(10, n // 8)
    for x in range(0, n + 1, step):
        ax.plot([x, x], [0.3, 0.4], color="#666", linewidth=0.7)
        ax.text(x, 0.0, f"{ws+x:,}", ha="center", va="bottom",
                fontsize=7, color="#444")


def option_a(fig, gridspec_start, ref, mapped, corrected, parest):
    """Option A: each row uses its own BAM's CIGAR placement."""
    title = "Option A — alignment-faithful (each row uses its own CIGAR)"
    rows = [("mapped", mapped), ("corrected", corrected), ("pA-rest", parest)]
    return _draw_block(fig, gridspec_start, ref, rows, title, show_ins=True)


def option_b(fig, gridspec_start, ref, mapped, corrected, parest):
    """Option B: anchor all rows to corrected's positioning.
    mapped row reuses corrected's bases+state but is labeled distinctly;
    we surface the CIGAR diff (insertions present in mapped but not
    corrected, or vice versa) above the mapped row.
    """
    title = "Option B — consistency-faithful (corrected alignment is canonical; mapped diffs annotated)"
    # mapped row uses corrected's per-position bases/states (so the sequence
    # is visually identical) but we tag insertions distinct from corrected.
    mapped_view = {
        "bases": corrected["bases"][:],
        "state": corrected["state"][:],
        "insertions": [],  # we'll synthesize "mapped-only insertion" markers
        "mismatches": corrected["mismatches"],
    }
    # Find insertions present in mapped that are NOT in corrected (or are different)
    corr_ins = {(rp, ln) for rp, ln, _ in corrected["insertions"]}
    for rp, ln, bases in mapped["insertions"]:
        if (rp, ln) not in corr_ins:
            mapped_view["insertions"].append((rp, ln, bases))
    rows = [("mapped*", mapped_view), ("corrected", corrected), ("pA-rest", parest)]
    return _draw_block(fig, gridspec_start, ref, rows, title, show_ins=True,
                       note="* = bases shown per corrected CIGAR. Purple ↓ above 'mapped*' "
                            "marks insertions in the mapped BAM that the correction absorbed.")


def _draw_block(fig, gs_top, ref, rows, title, show_ins=True, note=None):
    """Draw a comparison block: title → ref → rows → ticks. Returns top y."""
    panel_h_in = 0.32  # inches per row
    ref_h_in = 0.55
    title_h_in = 0.25
    tick_h_in = 0.25
    note_h_in = 0.20 if note else 0

    total_in = title_h_in + ref_h_in + len(rows) * panel_h_in + tick_h_in + note_h_in
    fig_height = fig.get_figheight()
    y_top = gs_top
    y_cursor = y_top

    # Title (as text via fig.text)
    fig.text(0.55, y_cursor - title_h_in / fig_height * 0.5,
             title, ha="center", va="center", fontsize=11, fontweight="bold")
    y_cursor -= title_h_in / fig_height

    # Ref row
    ax = fig.add_axes([0.13, y_cursor - ref_h_in / fig_height, 0.85, ref_h_in / fig_height])
    _draw_ref(ax, ref, WIN_START, WIN_END, ORIG_3P, CORR_3P)
    y_cursor -= ref_h_in / fig_height

    # Read rows
    for label, a in rows:
        ax = fig.add_axes([0.13, y_cursor - panel_h_in / fig_height, 0.85, panel_h_in / fig_height])
        _draw_row(ax, a, ref, WIN_START, WIN_END, label, show_ins=show_ins)
        y_cursor -= panel_h_in / fig_height

    # Tick row
    ax = fig.add_axes([0.13, y_cursor - tick_h_in / fig_height, 0.85, tick_h_in / fig_height])
    _draw_ticks(ax, WIN_START, WIN_END)
    y_cursor -= tick_h_in / fig_height

    # Note (if any)
    if note:
        fig.text(0.55, y_cursor - note_h_in / fig_height * 0.5,
                 note, ha="center", va="center", fontsize=8,
                 color="#666", style="italic")
        y_cursor -= note_h_in / fig_height

    return y_cursor


def main():
    ref = load_ref(CHROM, WIN_START, WIN_END)
    mapped = walk_cigar(fetch_read(MAPPED, QNAME), WIN_START, WIN_END)
    corrected = walk_cigar(fetch_read(CORR, QNAME), WIN_START, WIN_END)
    parest = walk_cigar(fetch_read(PAREST, QNAME), WIN_START, WIN_END)
    for a in (mapped, corrected, parest):
        compute_mismatches(a, ref)

    fig = plt.figure(figsize=(14, 9))
    fig.suptitle(f"cat1_minus_1 — chrII:{WIN_START}-{WIN_END}  (orig {ORIG_3P} → corr {CORR_3P}, −strand)",
                 fontsize=12, y=0.985, fontweight="bold")
    y_after_a = option_a(fig, 0.94, ref, mapped, corrected, parest)
    y_after_b = option_b(fig, y_after_a - 0.04, ref, mapped, corrected, parest)
    out = Path("/tmp/igv_snapshots/cat1_minus_1_compare.png")
    fig.savefig(str(out), dpi=150, bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
