#!/usr/bin/env python3
"""Per-read alignment plot for visual review of RECTIFY corrections.

Pure pysam + matplotlib — no IGV, ~0.5s per read. Renders one PNG per
validation read showing the genomic context and the three BAM versions of
the read (mapped / corrected / pA-restored), with mismatches, indels,
soft-clips, intron gaps, and 3'-end markers all visible at base resolution.

Layout (top to bottom):
  * Title bar
  * Optional **overview** panel (when the alignment spans an intron or large
    shift): a horizontal track summary showing the full aln span with N-op
    gaps drawn as thin connectors labeled by intron size.
  * **3'-end pileup** bar chart (corrected_3ends.plus/minus.bedgraph) with
    dashed lines at orig_3p (red) and corr_3p (green).
  * **Reference row** — colored genomic bases with position-aware ↓ markers
    above pointing to orig_3p and corr_3p.
  * **3 alignment rows** (mapped, corrected, pA-restored) with:
      - matches: light gray background, base character colored by identity
      - mismatches: pink background
      - deletions: orange '-' character
      - insertions: purple ↓ marker between positions, labeled "+N"
      - soft-clip bases: faded grey, rendered past aln boundaries
      - N-op (intron): grey '|' separator
  * **Position-tick row** at the bottom.

Multi-scale design:
  By default, the window is centered on corrected_3p ± pad (default 50 bp).
  When the alignment contains an N-op OR spans > 200 bp, an additional
  overview panel is rendered above showing the full alignment.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Optional

import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# IGV-style base colors
BASE_COLOR = {
    "A": "#2ECC40",
    "T": "#FF4136",
    "C": "#0074D9",
    "G": "#FF851B",
    "N": "#888888",
    "-": "#B26C00",
    "|": "#666666",
}
MISMATCH_BG = "#FFD0D0"
MATCH_BG = "#FAFAFA"
SOFTCLIP_BG = "#F0F0F0"
SOFTCLIP_FG = "#999999"
DELETION_BG = "#FFE0B2"
INTRON_BG = "#E0E0E0"


# ---------------------------------------------------------------------------
# Data extraction
# ---------------------------------------------------------------------------

def load_window_ref(genome_fa: Path, chrom: str, start: int, end: int) -> str:
    """0-based half-open window. Pads if outside chromosome bounds."""
    with pysam.FastaFile(str(genome_fa)) as fa:
        seq = fa.fetch(chrom, max(0, start), end)
    if start < 0:
        seq = ("N" * (-start)) + seq
    return seq.upper()


def fetch_read(bam_path: Path, qname: str) -> Optional[pysam.AlignedSegment]:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.query_name == qname:
                return r
    return None


def walk_cigar(read: pysam.AlignedSegment, win_start: int, win_end: int) -> dict:
    """Walk the CIGAR and return per-window structures.

    bases[i]              : read base at ref position (win_start + i), or None
    state[i]              : one of 'aligned', 'deletion', 'intron', 'softclip5',
                            'softclip3'  (only set if position has a read base)
    mismatches            : set of i where read != ref (computed later)
    insertions            : list of (ref_pos_just_before, length)
    n_op_intervals        : list of (n_start, n_end) for intron gaps overlapping window
    aln_start / aln_end   : read.reference_start / reference_end
    full_aln_start/end    : EXTENDED by soft-clip lengths (so the overview shows
                            soft-clips too)
    """
    n = win_end - win_start
    bases: list[Optional[str]] = [None] * n
    state: list[Optional[str]] = [None] * n
    insertions: list[tuple[int, int]] = []
    n_op_intervals: list[tuple[int, int]] = []

    out = {
        "bases": bases, "state": state,
        "insertions": insertions, "mismatches": set(),
        "n_op_intervals": n_op_intervals,
        "aln_start": None, "aln_end": None,
        "full_aln_start": None, "full_aln_end": None,
        "five_softclip_bp": 0, "three_softclip_bp": 0,
    }
    if read is None:
        return out

    seq = read.query_sequence or ""
    cigar = read.cigartuples or []
    out["aln_start"] = read.reference_start
    out["aln_end"] = read.reference_end

    # Handle 5' soft-clip: render bases at ref positions BEFORE aln_start.
    # qp starts at 0; rp starts at reference_start. If first op is S, the
    # clip bases are seq[0:length]; we place them at rp - length .. rp - 1.
    qp = 0
    rp = read.reference_start
    if cigar and cigar[0][0] == 4:
        clen = cigar[0][1]
        out["five_softclip_bp"] = clen
        # Place clip bases at ref positions rp-clen .. rp-1 (best-effort visual)
        for i_clip in range(clen):
            ref_pos = rp - clen + i_clip
            q_idx = i_clip
            if win_start <= ref_pos < win_end and 0 <= q_idx < len(seq):
                idx = ref_pos - win_start
                bases[idx] = seq[q_idx]
                state[idx] = "softclip5"
        qp = clen

    # Track 3' soft-clip for later placement
    three_softclip = 0
    if len(cigar) >= 2 and cigar[-1][0] == 4:
        three_softclip = cigar[-1][1]
    if len(cigar) >= 2 and cigar[-1][0] == 5 and len(cigar) >= 3 and cigar[-2][0] == 4:
        three_softclip = cigar[-2][1]
    out["three_softclip_bp"] = three_softclip

    # Walk ops, skipping leading soft-clip (already handled)
    cigar_iter = iter(cigar)
    if cigar and cigar[0][0] == 4:
        next(cigar_iter)
    for op, length in cigar_iter:
        if op == 5:  # hard-clip
            continue
        elif op == 4:  # trailing soft-clip — handled after loop
            break
        elif op in (0, 7, 8):  # M, =, X
            for _ in range(length):
                if win_start <= rp < win_end and 0 <= qp < len(seq):
                    idx = rp - win_start
                    bases[idx] = seq[qp]
                    state[idx] = "aligned"
                qp += 1
                rp += 1
        elif op == 1:  # insertion
            insertions.append((rp - 1, length))
            qp += length
        elif op == 2:  # deletion
            for _ in range(length):
                if win_start <= rp < win_end:
                    idx = rp - win_start
                    bases[idx] = "-"
                    state[idx] = "deletion"
                rp += 1
        elif op == 3:  # N intron
            n_start = rp
            for _ in range(length):
                if win_start <= rp < win_end:
                    idx = rp - win_start
                    bases[idx] = "|"
                    state[idx] = "intron"
                rp += 1
            n_op_intervals.append((n_start, rp))

    # 3' soft-clip: place bases at ref positions AFTER aln_end (visual).
    if three_softclip > 0:
        clip_q_start = qp
        for i_clip in range(three_softclip):
            ref_pos = read.reference_end + i_clip
            q_idx = clip_q_start + i_clip
            if win_start <= ref_pos < win_end and 0 <= q_idx < len(seq):
                idx = ref_pos - win_start
                bases[idx] = seq[q_idx]
                state[idx] = "softclip3"

    out["full_aln_start"] = read.reference_start - out["five_softclip_bp"]
    out["full_aln_end"] = read.reference_end + out["three_softclip_bp"]
    return out


def compute_mismatches(aligned: dict, ref_seq: str) -> None:
    ms: set[int] = set()
    for i, (b, st) in enumerate(zip(aligned["bases"], aligned["state"])):
        if b is None or st != "aligned":
            continue
        if b.upper() != ref_seq[i].upper():
            ms.add(i)
    aligned["mismatches"] = ms


def load_bedgraph_window(path: Path, chrom: str, win_start: int, win_end: int) -> list[tuple[int, float]]:
    out: list[tuple[int, float]] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 4:
                continue
            c, s, e, v = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            if c != chrom or e <= win_start or s >= win_end:
                continue
            for pos in range(max(s, win_start), min(e, win_end)):
                out.append((pos, v))
    return out


# ---------------------------------------------------------------------------
# Rendering primitives
# ---------------------------------------------------------------------------

def _color_for_base(b: str, state: str, in_mismatch: bool) -> tuple[str, str]:
    """Return (bg_color, fg_color) for one base."""
    if state in ("softclip5", "softclip3"):
        return SOFTCLIP_BG, SOFTCLIP_FG
    if state == "deletion":
        return DELETION_BG, BASE_COLOR["-"]
    if state == "intron":
        return INTRON_BG, BASE_COLOR["|"]
    if in_mismatch:
        return MISMATCH_BG, BASE_COLOR.get(b.upper(), "#000")
    return MATCH_BG, BASE_COLOR.get(b.upper(), "#000")


def render_alignment_row(ax, aligned: dict, ref_seq: str,
                         win_start: int, win_end: int, label: str):
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

    mismatches = aligned["mismatches"]
    for i, (b, st) in enumerate(zip(aligned["bases"], aligned["state"])):
        if b is None:
            continue
        bg, fg = _color_for_base(b, st, i in mismatches)
        ax.add_patch(Rectangle((i, 0.05), 1, 0.9, facecolor=bg, edgecolor="none"))
        weight = "normal" if st in ("softclip5", "softclip3") else "bold"
        char = b.upper()
        ax.text(i + 0.5, 0.5, char, ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight=weight)
    # Insertions: small purple ↓ + "+N" label ABOVE the row, anchored AT the
    # ref position immediately preceding the insertion. Ref coords remain rigid.
    for ref_pos, length in aligned["insertions"]:
        i = ref_pos - win_start + 1
        if 0 <= i <= n:
            ax.plot(i, 1.0, marker="v", markersize=6, color="#9c27b0", clip_on=False)
            ax.text(i, 1.12, f"+{length}", ha="center", va="bottom",
                    fontsize=6.5, color="#9c27b0", clip_on=False)


def render_ref_row(ax, ref_seq: str, win_start: int, win_end: int,
                   orig_3p: Optional[int], corr_3p: Optional[int]):
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel("ref", rotation=0, ha="right", va="center",
                  fontsize=9, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i, b in enumerate(ref_seq):
        bg = "#F5F5F5"
        fg = BASE_COLOR.get(b.upper(), "#000")
        ax.add_patch(Rectangle((i, 0.05), 1, 0.7, facecolor=bg, edgecolor="none"))
        ax.text(i + 0.5, 0.4, b.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight="bold")
    # Place ↓ markers ABOVE the reference row
    label_y = 1.05
    if orig_3p is not None and win_start <= orig_3p < win_end:
        x = (orig_3p - win_start) + 0.5
        ax.plot(x, label_y, marker="v", markersize=10, color="#d32f2f", clip_on=False)
        ax.text(x, label_y + 0.13, "orig", ha="center", va="bottom",
                fontsize=7, color="#d32f2f", clip_on=False)
    if corr_3p is not None and win_start <= corr_3p < win_end:
        x = (corr_3p - win_start) + 0.5
        ax.plot(x, label_y, marker="v", markersize=10, color="#2e7d32", clip_on=False)
        ax.text(x, label_y + 0.13, "corr", ha="center", va="bottom",
                fontsize=7, color="#2e7d32", clip_on=False)


def render_bedgraph(ax, plus_pts, minus_pts, win_start: int, win_end: int,
                    orig_3p: Optional[int], corr_3p: Optional[int]):
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylabel("3' end pileup", rotation=0, ha="right", va="center", fontsize=8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.set_xticks([])
    if plus_pts:
        xs = [p - win_start + 0.5 for p, _ in plus_pts]
        ys = [v for _, v in plus_pts]
        ax.bar(xs, ys, width=1.0, color="#2e7d32", alpha=0.85)
    if minus_pts:
        xs = [p - win_start + 0.5 for p, _ in minus_pts]
        ys = [-v for _, v in minus_pts]
        ax.bar(xs, ys, width=1.0, color="#c62828", alpha=0.85)
    if plus_pts or minus_pts:
        ax.axhline(0, color="#aaa", linewidth=0.5)
    if orig_3p is not None and win_start <= orig_3p < win_end:
        x = orig_3p - win_start + 0.5
        ax.axvline(x, color="#d32f2f", linestyle=":", linewidth=0.8, alpha=0.7)
    if corr_3p is not None and win_start <= corr_3p < win_end:
        x = corr_3p - win_start + 0.5
        ax.axvline(x, color="#2e7d32", linestyle="--", linewidth=1.0, alpha=0.8)


def render_axis_ticks(ax, win_start: int, win_end: int):
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylim(0, 0.4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_yticks([])
    step = max(10, n // 8)
    for x in range(0, n + 1, step):
        pos = win_start + x
        ax.plot([x, x], [0.3, 0.4], color="#666", linewidth=0.7)
        ax.text(x, 0.05, f"{pos:,}", ha="center", va="bottom",
                fontsize=7, color="#444")
    ax.set_xticks([])


def render_overview(ax, aligned_set, win_start: int, win_end: int,
                    orig_3p: Optional[int], corr_3p: Optional[int]):
    """Schematic overview of the full alignment span with N-op gaps drawn
    as thin connector lines labeled by intron size. One row per BAM version,
    stacked vertically. Useful when the read spans an intron or has a large
    shift that extends beyond the detail window."""
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylim(0, len(aligned_set) + 0.5)
    ax.set_yticks([])
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_ylabel("overview", rotation=0, ha="right", va="center",
                  fontsize=8.5, fontweight="bold")
    for k, (label, aligned) in enumerate(aligned_set):
        y_center = len(aligned_set) - k - 0.5
        if aligned["aln_start"] is None:
            continue
        a, b = aligned["full_aln_start"], aligned["full_aln_end"]
        if a is None or b is None:
            continue
        # 5' soft-clip
        sc5 = aligned["five_softclip_bp"]
        if sc5 > 0:
            x0 = max(0, a - win_start)
            x1 = max(0, aligned["aln_start"] - win_start)
            if x1 > x0:
                ax.add_patch(Rectangle((x0, y_center - 0.18), x1 - x0, 0.36,
                                       facecolor=SOFTCLIP_BG, edgecolor=SOFTCLIP_FG,
                                       linewidth=0.5))
        # aligned (incl. N gaps); break the bar at each N-op interval
        cur = aligned["aln_start"]
        for nstart, nend in aligned["n_op_intervals"]:
            if nstart > cur:
                x0 = max(0, cur - win_start)
                x1 = min(n, nstart - win_start)
                if x1 > x0:
                    ax.add_patch(Rectangle((x0, y_center - 0.20), x1 - x0, 0.40,
                                           facecolor="#1976d2", edgecolor="none",
                                           alpha=0.85))
            # gap connector
            gx0 = max(0, nstart - win_start)
            gx1 = min(n, nend - win_start)
            if gx1 > gx0:
                ax.plot([gx0, gx1], [y_center, y_center],
                        color="#aaa", linewidth=0.8)
                ax.text((gx0 + gx1) / 2, y_center + 0.20,
                        f"{nend - nstart} bp", ha="center", va="bottom",
                        fontsize=6.5, color="#666")
            cur = nend
        if cur < aligned["aln_end"]:
            x0 = max(0, cur - win_start)
            x1 = min(n, aligned["aln_end"] - win_start)
            if x1 > x0:
                ax.add_patch(Rectangle((x0, y_center - 0.20), x1 - x0, 0.40,
                                       facecolor="#1976d2", edgecolor="none",
                                       alpha=0.85))
        # 3' soft-clip
        sc3 = aligned["three_softclip_bp"]
        if sc3 > 0:
            x0 = max(0, aligned["aln_end"] - win_start)
            x1 = min(n, b - win_start)
            if x1 > x0:
                ax.add_patch(Rectangle((x0, y_center - 0.18), x1 - x0, 0.36,
                                       facecolor=SOFTCLIP_BG, edgecolor=SOFTCLIP_FG,
                                       linewidth=0.5))
        ax.text(-2, y_center, label, ha="right", va="center", fontsize=7.5)
    # markers
    if orig_3p is not None and win_start <= orig_3p < win_end:
        x = orig_3p - win_start + 0.5
        ax.axvline(x, color="#d32f2f", linestyle=":", linewidth=0.8, alpha=0.7)
    if corr_3p is not None and win_start <= corr_3p < win_end:
        x = corr_3p - win_start + 0.5
        ax.axvline(x, color="#2e7d32", linestyle="--", linewidth=1.0, alpha=0.8)


# ---------------------------------------------------------------------------
# Top-level render
# ---------------------------------------------------------------------------

def needs_overview(aligned_set, win_start: int, win_end: int) -> bool:
    """Decide if an overview panel is helpful."""
    for _, a in aligned_set:
        if any(True for _ in a["n_op_intervals"]):
            return True
        if a["aln_start"] is not None and a["aln_end"] is not None:
            if (a["aln_start"] < win_start - 5) or (a["aln_end"] > win_end + 5):
                return True
    return False


def overview_window(aligned_set) -> tuple[int, int]:
    """Compute the overview window: union of all full_aln_{start,end} ± 20 bp."""
    starts = [a["full_aln_start"] for _, a in aligned_set if a["full_aln_start"] is not None]
    ends   = [a["full_aln_end"]   for _, a in aligned_set if a["full_aln_end"]   is not None]
    if not starts or not ends:
        return (0, 1)
    return (max(0, min(starts) - 20), max(ends) + 20)


def render(
    *,
    qname: str,
    chrom: str,
    win_start: int,
    win_end: int,
    orig_3p: int,
    corr_3p: int,
    genome_fa: Path,
    mapped_bam: Path,
    corrected_bam: Path,
    parestore_bam: Path,
    bg_plus: Path,
    bg_minus: Path,
    out_path: Path,
    title: str = "",
) -> Path:
    ref_seq = load_window_ref(genome_fa, chrom, win_start, win_end)
    aligned_set: list[tuple[str, dict]] = []
    for label, bp in [("mapped", mapped_bam),
                      ("corrected", corrected_bam),
                      ("pA-rest", parestore_bam)]:
        read = fetch_read(bp, qname)
        a = walk_cigar(read, win_start, win_end)
        compute_mismatches(a, ref_seq)
        aligned_set.append((label, a))

    plus_pts = load_bedgraph_window(bg_plus, chrom, win_start, win_end)
    minus_pts = load_bedgraph_window(bg_minus, chrom, win_start, win_end)

    show_overview = needs_overview(aligned_set, win_start, win_end)
    ov_start, ov_end = (win_start, win_end)
    overview_aligned_set: list[tuple[str, dict]] = []
    if show_overview:
        ov_start, ov_end = overview_window(aligned_set)
        # Re-walk CIGARs in the overview window
        for label, bp in [("mapped", mapped_bam),
                          ("corrected", corrected_bam),
                          ("pA-rest", parestore_bam)]:
            read = fetch_read(bp, qname)
            a = walk_cigar(read, ov_start, ov_end)
            overview_aligned_set.append((label, a))

    n = win_end - win_start
    fig_w = max(7.5, n * 0.13)
    panels = ["bedgraph", "ref", "mapped", "corrected", "pA-rest", "ticks"]
    height_ratios = [1.0, 0.85, 0.75, 0.75, 0.75, 0.35]
    if show_overview:
        panels = ["overview", "ov_ticks"] + panels
        height_ratios = [0.9 + 0.18 * len(aligned_set), 0.30] + height_ratios
    fig_h = sum(height_ratios) * 0.45 + 0.5
    fig, axes = plt.subplots(
        len(panels), 1, figsize=(fig_w, fig_h),
        gridspec_kw={"height_ratios": height_ratios},
    )
    if title:
        fig.suptitle(title, fontsize=10, y=0.995)
    ax_iter = iter(axes)
    if show_overview:
        ax_ov = next(ax_iter)
        render_overview(ax_ov, overview_aligned_set, ov_start, ov_end, orig_3p, corr_3p)
        ax_ov_ticks = next(ax_iter)
        render_axis_ticks(ax_ov_ticks, ov_start, ov_end)
    ax_bg = next(ax_iter)
    render_bedgraph(ax_bg, plus_pts, minus_pts, win_start, win_end, orig_3p, corr_3p)
    ax_ref = next(ax_iter)
    render_ref_row(ax_ref, ref_seq, win_start, win_end, orig_3p, corr_3p)
    for label, a in aligned_set:
        ax = next(ax_iter)
        render_alignment_row(ax, a, ref_seq, win_start, win_end, label)
    ax_ticks = next(ax_iter)
    render_axis_ticks(ax_ticks, win_start, win_end)

    plt.subplots_adjust(left=0.11, right=0.99, top=0.95, bottom=0.04, hspace=0.10)
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--qname", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--orig-3p", type=int, required=True)
    p.add_argument("--corr-3p", type=int, required=True)
    p.add_argument("--mapped-bam", type=Path, required=True)
    p.add_argument("--corrected-bam", type=Path, required=True)
    p.add_argument("--parestore-bam", type=Path, required=True)
    p.add_argument("--genome-fa", type=Path, required=True)
    p.add_argument("--bg-plus", type=Path, required=True)
    p.add_argument("--bg-minus", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--title", default="")
    args = p.parse_args()
    render(
        qname=args.qname, chrom=args.chrom,
        win_start=args.start, win_end=args.end,
        orig_3p=args.orig_3p, corr_3p=args.corr_3p,
        genome_fa=args.genome_fa,
        mapped_bam=args.mapped_bam,
        corrected_bam=args.corrected_bam,
        parestore_bam=args.parestore_bam,
        bg_plus=args.bg_plus, bg_minus=args.bg_minus,
        out_path=args.out,
        title=args.title,
    )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
