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

# Base complement for RNA-orientation rendering of minus-strand reads.
# When the rendered read is minus-strand RNA, every base shown (ref + read +
# soft-clip + insertion) is complemented so the panel reads in RNA 5'→3'
# direction. The x-axis is also inverted at render-time so high-coord (RNA 5')
# sits on the left and low-coord (RNA 3') sits on the right. Together these
# two transforms convert the canonical genome+ display into a natural RNA
# view where poly-A tails appear as poly-A (not poly-T) at the 3' end.
_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
    "a": "t", "c": "g", "g": "c", "t": "a", "n": "n",
}


def _comp(b: str) -> str:
    """Complement a single base; non-ACGT (e.g. '-', '|') pass through."""
    return _COMPLEMENT.get(b, b)


def _maybe_comp(b: str, do_rc: bool) -> str:
    """Complement only when do_rc is True (minus-strand RNA orientation)."""
    return _comp(b) if (do_rc and b is not None) else b


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
            ins_bases = seq[qp:qp + length] if qp < len(seq) else ""
            insertions.append((rp, length, ins_bases))
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


def _decode_eq_chars(aligned: dict, ref_seq: str) -> None:
    """Decode `=`-encoded sequence bases in walk_cigar's output.

    Per SAM/BAM spec, a `=` character in the BAM SEQ field means the read base
    matches the reference at that position. Some aligners (notably mapPacBio
    and minimap2 with =/X CIGAR ops) emit `=`-encoded sequences for space
    savings. walk_cigar reads seq[read_pos] literally, so the bases array
    ends up with `=` chars where the read actually matches the reference.

    Without this step, every `=` base would (a) display as a literal `=`
    character in the alignment row and (b) trip compute_mismatches' comparison
    against the genomic ref (`=` != A/C/G/T), painting the whole row pink.

    This function rewrites every `aligned` `=` base to the corresponding
    ref_seq character so downstream rendering treats them as true matches.
    Also decodes insertion-pill base strings.
    """
    bases = aligned["bases"]
    for i, b in enumerate(bases):
        if b == '=':
            bases[i] = ref_seq[i]
    # Decode `=` in insertion pill bases too (insertions don't have a ref
    # position, but if the SEQ field uses `=` an insertion base could be `=`).
    # In practice insertions are by definition NOT matching ref, so `=` is
    # vanishingly rare in insertion content, but handle defensively.
    new_inserts = []
    for entry in aligned["insertions"]:
        if len(entry) >= 3:
            rp, length, ins_bases = entry[0], entry[1], entry[2]
            # `=` in an insertion has no defined ref base; leave as 'N' as a
            # safe fallback.
            ins_bases = ''.join('N' if c == '=' else c for c in ins_bases)
            new_inserts.append((rp, length, ins_bases))
        else:
            new_inserts.append(entry)
    aligned["insertions"] = new_inserts


def compute_mismatches(aligned: dict, ref_seq: str) -> None:
    # Decode any `=`-encoded SEQ bases before comparing — otherwise everything
    # in a =-encoded BAM (mapPacBio, some minimap2 outputs) gets flagged as
    # mismatch.
    _decode_eq_chars(aligned, ref_seq)
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
                         win_start: int, win_end: int, label: str,
                         is_reverse: bool = False):
    """Render a read row with prominent INLINE insertion bases.

    Each row shows:
      - Bases at ref positions per the CIGAR (M/=/X → row chars;
        D → orange '-'; N → grey '|'; S → faded chars past aln bounds).
      - Insertion overlays as full-size color-coded base characters wrapped
        in vivid purple-bordered boxes, wedged between ref columns at the
        insertion site. Stagger heights when nearby insertions would
        collide so they don't overlap.

    ``is_reverse=True`` complements every displayed base (matches, mismatches,
    soft-clips, insertion contents) for minus-strand RNA orientation. The
    caller is responsible for ``ax.invert_xaxis()`` so columns also flip.
    """
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylim(-0.05, 2.60)  # headroom for staggered insertion pills + arrows
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

    mismatches = aligned["mismatches"]
    for i, (b, st) in enumerate(zip(aligned["bases"], aligned["state"])):
        if b is None:
            continue
        # Complement the displayed character for minus-strand RNA orientation.
        # Mismatch detection still uses the original b vs ref (computed upstream
        # in compute_mismatches); the complement is a display-only transform.
        b_disp = _maybe_comp(b, is_reverse)
        bg, fg = _color_for_base(b_disp, st, i in mismatches)
        ax.add_patch(Rectangle((i, 0.05), 1, 0.9, facecolor=bg, edgecolor="none"))
        weight = "normal" if st in ("softclip5", "softclip3") else "bold"
        ax.text(i + 0.5, 0.5, b_disp.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight=weight)

    # Insertions: full-size colored bases wedged in a vivid purple-bordered
    # box, positioned ABOVE the row with a black ▼ pointing down to the exact
    # insertion site (between ref columns rp-1 and rp). Nearby insertions
    # stagger vertically so their boxes don't overlap.
    insertions = sorted(aligned["insertions"], key=lambda t: t[0])
    last_x_end = -100
    stagger = 0
    for entry in insertions:
        if len(entry) >= 3:
            rp, length, bases = entry[0], entry[1], entry[2]
        else:
            rp, length = entry[0], entry[1]
            bases = ""
        ix = rp - win_start
        if ix < 0 or ix > n:
            continue
        char_width = 0.55
        total_w = length * char_width
        start_x = ix - total_w / 2
        end_x = ix + total_w / 2
        if start_x < last_x_end + 0.3:
            stagger = 1 - stagger
        else:
            stagger = 0
        # Two-tier vertical positioning, both ABOVE the row (y=1.0 = row top)
        pill_y = 1.25 if stagger == 0 else 1.85
        pill_h = 0.55
        # Black ▼ pointing down from pill bottom to the exact insertion site
        ax.plot(ix, 1.02, marker="v", markersize=7, color="#000",
                clip_on=False, zorder=5)
        # Thin black guide line from row top up to the pill bottom
        ax.plot([ix, ix], [1.05, pill_y - 0.02],
                color="#000", linewidth=0.6, clip_on=False, zorder=4)
        # Purple-bordered pill
        ax.add_patch(Rectangle((start_x - 0.10, pill_y),
                               total_w + 0.20, pill_h,
                               facecolor="#E1BEE7", edgecolor="#7B1FA2",
                               linewidth=1.2, clip_on=False, zorder=3))
        # Complement the inserted bases when rendering minus-strand RNA orientation.
        bases_disp = "".join(_maybe_comp(c, is_reverse) for c in bases) if bases else bases
        if length <= 4 and bases_disp:
            for j, b in enumerate(bases_disp):
                x = start_x + (j + 0.5) * char_width
                fg = BASE_COLOR.get(b.upper(), "#7B1FA2")
                ax.text(x, pill_y + pill_h / 2, b.upper(),
                        ha="center", va="center",
                        fontsize=8.5, family="monospace",
                        color=fg, fontweight="bold", clip_on=False, zorder=4)
        else:
            label_text = f"+{length}"
            if bases_disp:
                label_text += f" {bases_disp[:3]}…" if length > 3 else f" {bases_disp}"
            ax.text(ix, pill_y + pill_h / 2, label_text, ha="center", va="center",
                    fontsize=7.0, color="#7B1FA2", fontweight="bold",
                    clip_on=False, zorder=4)
        last_x_end = end_x


def render_ref_row(ax, ref_seq: str, win_start: int, win_end: int,
                   orig_3p: Optional[int], corr_3p: Optional[int],
                   label: str = "ref", is_reverse: bool = False):
    n = win_end - win_start
    ax.set_xlim(0, n)
    # Reserve space ABOVE the bases for the orig/corr markers + labels.
    # Keeping everything INSIDE the panel bounds prevents text bleeding
    # into the bedgraph panel above.
    ax.set_ylim(0, 1.40)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center",
                  fontsize=9, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Ref bases sit at y=0.05-0.55 so there's room for ↓ markers above.
    # For minus-strand RNA orientation, complement each base before display
    # (caller is responsible for inverting the x-axis).
    for i, b in enumerate(ref_seq):
        b_disp = _maybe_comp(b, is_reverse)
        bg = "#F5F5F5"
        fg = BASE_COLOR.get(b_disp.upper(), "#000")
        ax.add_patch(Rectangle((i, 0.05), 1, 0.50, facecolor=bg, edgecolor="none"))
        ax.text(i + 0.5, 0.30, b_disp.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight="bold")
    # ↓ markers ABOVE the ref bases. Text label is placed JUST ABOVE the
    # triangle (y=1.20) so it stays well inside the panel bounds (ylim=1.40)
    # and doesn't bleed into the bedgraph above. With more vertical headroom
    # for the marker stack, the labels are also further from the overview
    # ruler tick labels.
    tri_y = 0.95
    txt_y = 1.20
    if orig_3p is not None and win_start <= orig_3p < win_end:
        x = (orig_3p - win_start) + 0.5
        ax.plot(x, tri_y, marker="v", markersize=11, color="#d32f2f")
        ax.text(x, txt_y, "orig", ha="center", va="bottom",
                fontsize=8, color="#d32f2f", fontweight="bold")
    if corr_3p is not None and win_start <= corr_3p < win_end:
        x = (corr_3p - win_start) + 0.5
        ax.plot(x, tri_y, marker="v", markersize=11, color="#2e7d32")
        ax.text(x, txt_y, "corr", ha="center", va="bottom",
                fontsize=8, color="#2e7d32", fontweight="bold")


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
    """Render position labels along the bottom of a panel.

    Tick + label are placed at the CENTER of the ref column for that
    position (x = pos - win_start + 0.5). In 0-based half-open coords,
    position N occupies the column [N, N+1), so its base character sits
    at column center N+0.5. Aligning the tick to the center makes it
    visually point at the labeled base.
    """
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1.0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([])
    # Target ~6 labels per panel; round step up to multiple of 5 for readability
    raw_step = max(10, n // 6)
    step = max(10, (raw_step // 5) * 5 or raw_step)
    for x in range(0, n + 1, step):
        pos = win_start + x
        cx = x + 0.5  # center of the ref column for `pos`
        ax.plot([cx, cx], [0.70, 0.95], color="#666", linewidth=0.8)
        ax.text(cx, 0.05, f"{pos:,}", ha="center", va="bottom",
                fontsize=7.5, color="#333")


def render_overview(ax, aligned_set, win_start: int, win_end: int,
                    orig_3p: Optional[int], corr_3p: Optional[int]):
    """Schematic overview of the full alignment span with N-op gaps drawn
    as thin connector lines labeled by intron size. One row per BAM version,
    stacked vertically. Useful when the read spans an intron or has a large
    shift that extends beyond the detail window."""
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylim(0, len(aligned_set) + 0.5)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Use y-ticks for the per-row labels so matplotlib auto-spaces them
    # (no more collision with the panel-level "overview" ylabel).
    yticks = [len(aligned_set) - k - 0.5 for k in range(len(aligned_set))]
    ylabels = [lbl for lbl, _ in aligned_set]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7.5)
    ax.tick_params(axis="y", length=0, pad=2)
    # Small "overview" tag at the top-left corner of the panel, OUTSIDE the
    # y-tick label column so it doesn't stack on top of the row labels.
    ax.text(-0.02, 1.02, "overview", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7,
            color="#888", fontweight="bold")
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
        # (row label now drawn as a y-tick — see set_yticklabels above)
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
    minimap2_bam: Path | None = None,
    winner_label: str = "",
) -> Path:
    """Render per-base alignment plot with up to 4 alignment rows.

    Tracks rendered top-to-bottom:
      1. ``minimap2``        — reference aligner row (always shown if
                                ``minimap2_bam`` is given).
      2. ``winner: <name>``  — the winning aligner from consensus selection,
                                only added when ``mapped_bam != minimap2_bam``.
      3. ``corrected``       — the rectified BAM (CIGAR-surgery output).
      4. ``pA-rest``         — Step-4 BAM with poly-A tail restored as soft-clip.

    Args:
        mapped_bam:    The winning aligner's BAM. If ``minimap2_bam`` is None
                       or identical to ``mapped_bam``, this is the only
                       upstream row, labeled simply ``minimap2``.
        minimap2_bam:  Optional path to the minimap2 aligner BAM. If different
                       from ``mapped_bam``, an extra ``minimap2`` row is
                       inserted above the winner row.
        winner_label:  Display name for the winning aligner (e.g. ``mapPacBio``).
                       Ignored when winner == minimap2.
    """
    ref_seq = load_window_ref(genome_fa, chrom, win_start, win_end)

    # ---- Build the list of upstream alignment tracks (in display order) ----
    # If a minimap2 BAM is supplied AND differs from the winner BAM, we render
    # TWO upstream rows: minimap2 + winner. Otherwise we render one row
    # labeled "minimap2" (which is the winner — minimap2 is its own winner).
    show_minimap2_row = (
        minimap2_bam is not None and Path(minimap2_bam) != Path(mapped_bam)
    )
    upstream_tracks: list[tuple[str, Path]] = []
    if show_minimap2_row:
        upstream_tracks.append(("minimap2", Path(minimap2_bam)))
        winner_disp = f"winner: {winner_label}" if winner_label else "winner"
        upstream_tracks.append((winner_disp, mapped_bam))
    else:
        # Single upstream row. If winner_label is set and == "minimap2", or
        # not set, just call it minimap2. Otherwise show the winner name.
        single_label = winner_label or "minimap2"
        # When user passed winner_label == "minimap2" (or omitted it), keep it
        # simple; if it's some other aligner and minimap2_bam was None, label
        # it as that winner (rare validation-set path).
        upstream_tracks.append((single_label, mapped_bam))

    track_sources: list[tuple[str, Path]] = list(upstream_tracks) + [
        ("corrected", corrected_bam),
        ("pA-rest",   parestore_bam),
    ]

    aligned_set: list[tuple[str, dict]] = []
    mapped_read_obj = None
    for label, bp in track_sources:
        read = fetch_read(bp, qname)
        # ``mapped_read_obj`` = the winner's record; used for orig_3p override
        # below. The winner is the LAST upstream row (either minimap2 if
        # single, or winner row when both are shown).
        if bp == mapped_bam and mapped_read_obj is None:
            mapped_read_obj = read
        a = walk_cigar(read, win_start, win_end)
        compute_mismatches(a, ref_seq)
        aligned_set.append((label, a))

    # Override orig_3p with the mapped BAM's actual aln_start (− strand) or
    # aln_end - 1 (+ strand). The corrected_reads.tsv's `original_3prime` is
    # computed against the post-merge / consensus input BAM, which can have a
    # different aln_start than the per-aligner BAM we're rendering (the input
    # BAM may have absorbed the per-aligner soft-clip into M ops). Showing
    # the orig arrow against the BAM the row actually displays is the honest
    # representation of "where minimap2/mapPacBio called the 3' end."
    if mapped_read_obj is not None:
        if mapped_read_obj.is_reverse:
            orig_3p = mapped_read_obj.reference_start
        else:
            orig_3p = mapped_read_obj.reference_end - 1

    # Override corr_3p with the rectified BAM's actual final aln boundary.
    # `corrected_reads.tsv`'s `corrected_3prime` records the intermediate
    # output from the correction modules (Module 2C / 2E / 2G / etc.), but
    # `bam_writer._hardclip_trailing_a_run` may subsequently clip additional
    # genomic A-rich positions to land at a clean body boundary. For minus-
    # strand cat1_minus_1 this is a 7 bp gap (TSV records 9827; BAM ends at
    # 9834 = the G of `...AGTG` before the poly-A tail). Show the corr arrow
    # at the BAM's actual final position so the visualization reflects what
    # IGV would show.
    corrected_read_obj = fetch_read(corrected_bam, qname)
    if corrected_read_obj is not None:
        if corrected_read_obj.is_reverse:
            corr_3p = corrected_read_obj.reference_start
        else:
            corr_3p = corrected_read_obj.reference_end - 1

    plus_pts = load_bedgraph_window(bg_plus, chrom, win_start, win_end)
    minus_pts = load_bedgraph_window(bg_minus, chrom, win_start, win_end)

    show_overview = needs_overview(aligned_set, win_start, win_end)
    ov_start, ov_end = (win_start, win_end)
    overview_aligned_set: list[tuple[str, dict]] = []
    if show_overview:
        ov_start, ov_end = overview_window(aligned_set)
        # Re-walk CIGARs in the overview window. Use the SAME track_sources
        # list so the overview matches the per-base panel ordering exactly.
        for label, bp in track_sources:
            read = fetch_read(bp, qname)
            a = walk_cigar(read, ov_start, ov_end)
            overview_aligned_set.append((label, a))

    n = win_end - win_start
    fig_w = max(7.5, n * 0.13)
    # Each entry in track_sources becomes one alignment panel. Common case:
    # 3 panels (minimap2-only + corrected + pA-rest) or 4 panels (minimap2 +
    # winner + corrected + pA-rest) when winner != minimap2.
    track_labels = [lbl for lbl, _ in track_sources]
    panels = ["bedgraph", "ref"] + track_labels + ["ticks"]
    # Ref row a hair taller (1.0 → 1.2) for the orig/corr label stack.
    # Read rows much taller (2.0) for the insertion-pill stagger ABOVE the row.
    height_ratios = [1.0, 1.2] + [2.0] * len(track_sources) + [0.35]
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
    # For minus-strand reads, render in RNA 5'→3' orientation: complement
    # every displayed base AND invert the x-axis so the higher genomic coord
    # (which is the RNA 5' end) sits on the LEFT and the lower genomic coord
    # (RNA 3' end, where the CPA + poly-A tail are) sits on the RIGHT. This
    # makes the poly-A tail appear as poly-A (not poly-T) at the natural 3'
    # end position and lines up with how the user reads RNA sequences.
    is_reverse = bool(mapped_read_obj is not None and mapped_read_obj.is_reverse)

    ax_bg = next(ax_iter)
    render_bedgraph(ax_bg, plus_pts, minus_pts, win_start, win_end, orig_3p, corr_3p)
    ax_ref = next(ax_iter)
    render_ref_row(ax_ref, ref_seq, win_start, win_end, orig_3p, corr_3p,
                   is_reverse=is_reverse)
    for label, a in aligned_set:
        ax = next(ax_iter)
        render_alignment_row(ax, a, ref_seq, win_start, win_end, label,
                             is_reverse=is_reverse)
    ax_ticks = next(ax_iter)
    render_axis_ticks(ax_ticks, win_start, win_end)

    # Flip the x-axis on every panel for minus-strand reads (RNA 5'→3' = left
    # to right means higher coord on the left). Done after all rendering so
    # the internal column coordinates stay genomic-oriented and only the
    # display flips.
    if is_reverse:
        for ax in axes:
            ax.invert_xaxis()

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
    p.add_argument("--mapped-bam", type=Path, required=True,
                   help="Winning aligner's BAM. If equal to --minimap2-bam (or "
                        "the latter is omitted), only one upstream alignment "
                        "row is rendered.")
    p.add_argument("--minimap2-bam", type=Path, default=None,
                   help="Optional minimap2 BAM. When supplied AND different "
                        "from --mapped-bam, adds a top 'minimap2' row above "
                        "the winner row.")
    p.add_argument("--winner-label", default="",
                   help="Display name for the winning aligner (e.g. "
                        "'mapPacBio'). Used only when minimap2 and winner "
                        "are different BAMs.")
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
        minimap2_bam=args.minimap2_bam,
        winner_label=args.winner_label,
        corrected_bam=args.corrected_bam,
        parestore_bam=args.parestore_bam,
        bg_plus=args.bg_plus, bg_minus=args.bg_minus,
        out_path=args.out,
        title=args.title,
    )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
