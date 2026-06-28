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
import datetime
import gzip
from pathlib import Path
from typing import Optional

import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Font stack: Helvetica Neue → Helvetica → DejaVu Sans for sans-serif;
# Menlo → Monaco → DejaVu Sans Mono for monospace. macOS-native pairing
# (consistent across the figure, falls back gracefully on Linux/Windows).
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica Neue', 'Helvetica',
                                           'Arial', 'DejaVu Sans']
matplotlib.rcParams['font.monospace'] = ['Menlo', 'Monaco', 'Consolas',
                                          'DejaVu Sans Mono']

# Consolidated text-size constants. Most text in the figure is at TEXT;
# tiny inline CIGAR-op labels need to be SMALL.
SIZE_TITLE = 11
SIZE_TEXT  = 9
SIZE_SMALL = 6.5

# Toggle the amber divergence-band overlay in the overview track. When False,
# the per-position disagreement highlighting is suppressed entirely. Flip to
# True to re-enable for ad-hoc reviews of where aligners diverge.
SHOW_DIVERGENCE_SHADING = False

# Trim-metadata TSV produced by `rectify trim-polya`. Used to distinguish
# the pA-tail portion of the trailing soft-clip (the bases that were
# removed by trim-polya and re-attached after alignment) from the
# aligner-emitted soft-clip of the trimmed body.
_REPO_ROOT = Path(__file__).resolve().parents[2]
POLYA_TRIM_TSV = (_REPO_ROOT / "scripts" / "validation_data"
                  / "rebuild_2026_05" / "trimmed"
                  / "validation_reads_polya_trim_metadata.tsv")


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
SOFTCLIP_BG = "#F0F0F0"      # aligner-emitted soft-clip on the trimmed body
SOFTCLIP_FG = "#999999"
# pA-tail portion of the trailing soft-clip — the bases that were removed
# by `rectify trim-polya` and re-attached after alignment. Distinguishes
# the biological tail from the aligner's technical soft-clip of the body.
SOFTCLIP_PA_BG = "#C8E6C9"
SOFTCLIP_PA_FG = "#2E7D32"
# Force-aligned-then-walkback-removed portion of the tail: bases the ALIGNER
# aligned (force-aligned past the true 3' end) that RECTIFY's walkback clipped.
# Distinct from a NATIVE aligner soft-clip (grey, SOFTCLIP_BG) — the aligner
# CHOSE to clip those; the slate-blue ones it force-aligned and RECTIFY removed.
# Deliberately a SLATE-BLUE hue (not red): the old light-red was too close to the
# pink mismatch fill (MISMATCH_BG) and the two were confusable.
SOFTCLIP_WB_BG = "#B0C4DE"
SOFTCLIP_WB_FG = "#2c3e50"
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


_POLYA_LEN_CACHE: dict[Path, dict[str, int]] = {}


def _load_polya_lengths(parquet_tsv: Path) -> dict[str, int]:
    """Load read_id → length of the chunk that ``rectify trim-polya``
    removed (and that ``restore_polya_from_parquet.py`` re-attached as
    soft-clip) from the trim-metadata TSV. This is
    ``len(trimmed_3prime_seq)`` = ``polya_len`` + ``len(adapter_seq)``,
    i.e. what physically appears in the trailing soft-clip after restore.
    """
    if parquet_tsv in _POLYA_LEN_CACHE:
        return _POLYA_LEN_CACHE[parquet_tsv]
    out: dict[str, int] = {}
    try:
        import csv
        with open(parquet_tsv) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rid = row.get("read_id")
                if not rid:
                    continue
                seq = (row.get("trimmed_3prime_seq") or "").strip()
                if seq:
                    out[rid] = len(seq)
                    continue
                # Fallback: polya_len if no seq column (older parquet).
                try:
                    out[rid] = int(row.get("polya_len") or 0)
                except ValueError:
                    out[rid] = 0
    except Exception:
        pass
    _POLYA_LEN_CACHE[parquet_tsv] = out
    return out


def _mark_polya_softclip(aligned: dict, polya_len: int, is_reverse: bool) -> None:
    """Relabel the trim-attached portion of the trailing soft-clip in the
    state array so it can render in a distinct color from the aligner-
    emitted soft-clip of the trimmed body.

    The boundary is taken VERBATIM from the parquet: the last
    ``polya_len`` BAM positions of the soft-clip are the chunk that
    ``rectify trim-polya`` removed and ``restore_polya_from_parquet.py``
    re-attached. No inward extension — if the aligner soft-clipped a few
    trailing A's of the trimmed body those stay grey and visible, which
    is the honest signal that trim-polya did not walk all the way to the
    first non-A.

    pA bases in BAM orientation:
      - plus strand: bases live at the BAM-3' end (3' soft-clip).
      - minus strand: BAM is reverse-complemented so the chunk lives at
        the BAM-5' end (5' soft-clip).
    """
    if polya_len <= 0:
        return
    state = aligned.get("state") or []
    if is_reverse:
        target = "softclip5"
        sc_total = aligned.get("five_softclip_bp", 0)
    else:
        target = "softclip3"
        sc_total = aligned.get("three_softclip_bp", 0)
    if sc_total <= 0:
        return
    pa_n = min(polya_len, sc_total)
    indices = [i for i, st in enumerate(state) if st == target]
    if not indices:
        return
    if is_reverse:
        # 5' soft-clip: BAM pos 0 = outer end = indices[0]. The restored
        # chunk lives at BAM positions [0, pa_n).
        pa_indices = indices[:pa_n]
    else:
        # 3' soft-clip: BAM pos sc-1 = outer end = indices[-1]. The
        # restored chunk lives at BAM positions [sc-pa_n, sc).
        pa_indices = indices[-pa_n:]
    for idx in pa_indices:
        state[idx] = "softclip_pa"
    aligned["polya_len_effective"] = pa_n


def _mark_walkback_removed(aligned: dict, wb_n: int, is_reverse: bool) -> None:
    """Relabel the INNERMOST ``wb_n`` cells of the 3' soft-clip (those nearest
    the corrected body) as ``softclip_wb`` — the bases the aligner FORCE-ALIGNED
    past the true 3' end that RECTIFY's walkback clipped.

    Call AFTER :func:`_mark_polya_softclip`: it only relabels cells still tagged
    as the plain (grey) aligner soft-clip, so the parquet pA cells (green, the
    OUTER end) and these force-aligned cells (crimson, the INNER end) never
    overlap. Any cells left grey between them are the aligner's genuine native
    soft-clip. ``wb_n`` is clamped to the available grey cells, so when the
    force-aligned extent exceeds the reconstructed tail the row simply reads as
    all-crimson + pA (mapPacBio) vs mostly-grey + a crimson sliver (deSALT)."""
    if wb_n <= 0:
        return
    state = aligned.get("state") or []
    target = "softclip5" if is_reverse else "softclip3"
    # Only cells still grey (parquet pA already relabeled softclip_pa, excluded).
    indices = [i for i, st in enumerate(state) if st == target]
    if not indices:
        return
    wb_n = min(wb_n, len(indices))
    # Innermost = nearest the corrected body. Minus: 3' tail is at LOW genomic,
    # so the inner cells are the HIGH columns (end of the list). Plus: inner
    # cells are the LOW columns (start of the list).
    wb_indices = indices[-wb_n:] if is_reverse else indices[:wb_n]
    for idx in wb_indices:
        state[idx] = "softclip_wb"


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
        # Raw CIGAR + ref start + query SEQ for the overview-painter to
        # walk again (and expand M ops into =/X based on per-base match).
        "cigartuples": [], "ref_start_raw": None, "query_seq_raw": "",
    }
    if read is None:
        return out

    seq = read.query_sequence or ""
    cigar = read.cigartuples or []
    out["cigartuples"] = list(cigar)
    out["ref_start_raw"] = read.reference_start
    out["query_seq_raw"] = seq
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
    if state == "softclip_pa":
        return SOFTCLIP_PA_BG, SOFTCLIP_PA_FG
    if state == "softclip_wb":
        return SOFTCLIP_WB_BG, SOFTCLIP_WB_FG
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
                         is_reverse: bool = False,
                         corr_3p_aligner: Optional[int] = None):
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
    ax.set_ylim(-0.28, 2.60)  # headroom for insertion pills above + corr-3' mark below
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=SIZE_TEXT)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # RECTIFY-inferred corrected 3' end for THIS row's aligner (from the
    # per-aligner summary walkback result). A green ▲ just below the row + a
    # thin green guide marks where walkback placed this aligner's 3' end, so
    # the non-winner rows show their own walkback outcome (e.g. an over-
    # extending aligner converging — or not — to the consensus boundary).
    if corr_3p_aligner is not None and win_start <= corr_3p_aligner < win_end:
        # Place the ▲ at the BOUNDARY between the called 3'-end base and the
        # first downstream (removed-tail) nt, not centered under the 3'-end
        # base. The 3'-end base occupies data-cell [c, c+1] with c =
        # corr_3p-win_start. Downstream is ref+1 on plus strand (right edge,
        # offset +1.0) and ref-1 on minus strand (left edge, offset +0.0 —
        # the later invert_xaxis() flips it to the visual right, matching the
        # removed tail). See the 0-based half-open 3'-end convention in
        # CLAUDE.md.
        _off = 1.0 if not is_reverse else 0.0
        _cx = (corr_3p_aligner - win_start) + _off
        # ▲ below the row only — no vertical line through the bases (it obscured
        # the base letter at the corrected-3' column).
        ax.plot(_cx, -0.16, marker="^", markersize=7, color="#2e7d32",
                clip_on=False, zorder=7)

    # (Walkback-removed bases are shown by recoloring the innermost soft-clip
    # cells crimson via _mark_walkback_removed — see SOFTCLIP_WB_BG — so the
    # tail reads as non-overlapping crimson/grey/green segments without an
    # overlay obscuring the bases.)

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
        char_width = 1.0     # match the normal-cell width so inserted letters
                             # render at the same visual size as body bases
        total_w = length * char_width
        start_x = ix - total_w / 2
        end_x = ix + total_w / 2
        if start_x < last_x_end + 0.3:
            stagger = 1 - stagger
        else:
            stagger = 0
        # Two-tier vertical positioning, both ABOVE the row (y=1.0 = row top)
        pill_y = 1.25 if stagger == 0 else 1.85
        pill_h = 0.70    # taller pill so the inserted letter has buffer above/below
        # Black ▼ pointing down from pill bottom to the exact insertion site
        ax.plot(ix, 1.02, marker="v", markersize=7, color="#000",
                clip_on=False, zorder=5)
        # Thin black guide line from row top up to the pill bottom
        ax.plot([ix, ix], [1.05, pill_y - 0.02],
                color="#000", linewidth=0.6, clip_on=False, zorder=4)
        # Light purple highlight (no border, generous padding around the
        # letter so it visually breathes).
        ax.add_patch(Rectangle((start_x - 0.20, pill_y),
                               total_w + 0.40, pill_h,
                               facecolor="#E1BEE7", edgecolor="none",
                               clip_on=False, zorder=3))
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
                    fontsize=SIZE_SMALL, color="#7B1FA2", fontweight="bold",
                    clip_on=False, zorder=4)
        last_x_end = end_x


def render_ref_row(ax, ref_seq: str, win_start: int, win_end: int,
                   orig_3p: Optional[int], corr_3p: Optional[int],
                   label: str = "ref", is_reverse: bool = False,
                   samtools_3p: Optional[int] = None,
                   orig_5p: Optional[int] = None,
                   corr_5p: Optional[int] = None):
    n = win_end - win_start
    ax.set_xlim(0, n)
    # The ref bases sit at y=0.05-0.55. The old orig/corr/samtools markers
    # (which needed a tall ylim) are gone, so keep the panel TIGHT (top just
    # above the bases) — this lets the overview→ref zoom lines connect right at
    # the ref-seq top. Only the rare 5' markers (zoom-out windows that include
    # the 5' end) need the taller headroom; reserve it only then.
    _has5 = ((orig_5p is not None and win_start <= orig_5p < win_end)
             or (corr_5p is not None and win_start <= corr_5p < win_end))
    _ref_top = 1.75 if _has5 else 0.72
    ax.set_ylim(0, _ref_top)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center",
                  fontsize=SIZE_TEXT, fontweight="bold")
    # Anchor the ylabel at the ref-bases level (bases centre ~y=0.30).
    ax.yaxis.set_label_coords(-0.018, 0.30 / _ref_top)
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
    # The ref row is kept CLEAN (2026-06-27): the old orig=corr AND samtools
    # markers were removed. The naive "without-RECTIFY" 3' end (= samtools view
    # of the unrectified read) is now shown directly as the green ▲ on the
    # minimap2 (unrectified) per-aligner row, and each rectified row carries its
    # own corrected-3' ▲ — so all 3'-end attention lives on the per-aligner
    # tracks. tri_y / txt_y_* below are retained only for the rare 5' markers.
    tri_y = 0.95
    txt_y_lo = 1.15
    txt_y_hi = 1.60

    # 5' end markers — orig_5p (faded blue) and corr_5p (deep blue). These
    # only appear when the rendered window happens to include the 5' end
    # (rare for the default 3'-end-centered window, but present in zoom-
    # out views and overview windows that span the full read).
    blue_orig = "#42a5f5"
    blue_corr = "#1565c0"
    orig5_in = orig_5p is not None and win_start <= orig_5p < win_end
    corr5_in = corr_5p is not None and win_start <= corr_5p < win_end
    same5 = orig5_in and corr5_in and orig_5p == corr_5p
    if same5:
        x = (orig_5p - win_start) + 0.5
        ax.plot(x, tri_y, marker="v", markersize=11, color=blue_corr)
        ax.text(x, txt_y_lo, "orig=corr 5'", ha="center", va="bottom",
                fontsize=SIZE_TEXT, color=blue_corr, fontweight="bold")
    else:
        if orig5_in:
            x = (orig_5p - win_start) + 0.5
            ax.plot(x, tri_y, marker="v", markersize=11, color=blue_orig)
            ax.text(x, txt_y_lo, "orig 5'", ha="center", va="bottom",
                    fontsize=SIZE_TEXT, color=blue_orig, fontweight="bold")
        if corr5_in:
            x = (corr_5p - win_start) + 0.5
            ax.plot(x, tri_y, marker="v", markersize=11, color=blue_corr)
            close_5 = (
                orig5_in and abs(corr_5p - orig_5p) <= 3
            )
            y5 = txt_y_hi if close_5 else txt_y_lo
            ax.text(x, y5, "corr 5'", ha="center", va="bottom",
                    fontsize=SIZE_TEXT, color=blue_corr, fontweight="bold")


def render_bedgraph(ax, plus_pts, minus_pts, win_start: int, win_end: int,
                    orig_3p: Optional[int], corr_3p: Optional[int],
                    samtools_3p: Optional[int] = None):
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylabel("3' end pileup", rotation=0, ha="right", va="center", fontsize=SIZE_TEXT)
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
    # samtools (naive 3') line removed — the naive 3' end is now shown by the
    # green ▲ on the minimap2 (unrectified) per-aligner row instead.


_ALIGNER_ABBREV = {
    "minimap2": "mm2", "gapmm2": "gap", "mapPacBio": "mPB",
    "deSALT": "deS", "uLTRA": "uLT",
}


def render_aligner_3p_agreement(ax, corr3p_by_aligner: dict, winner_name: str,
                                win_start: int, win_end: int) -> None:
    """Per-aligner corrected-3'-end agreement track.

    Replaces the old bundle-bedgraph "pileup" (which, on the 1-read-per-locus
    validation set, was always a single count=1 bar). Instead, collapse every
    aligner's RECTIFY-corrected 3' end onto the genomic axis and stack them:
    each aligner is one cell at its corrected position, so a single tall column
    = full agreement (consensus), and multiple columns = the aligners disagree
    (which read/aligner landed where is then obvious). Winner cell in green,
    the rest steel-blue; abbreviation inside each cell.
    """
    from collections import defaultdict
    n = win_end - win_start
    ax.set_xlim(0, n)
    ax.set_ylabel("aligner 3′\nagreement", rotation=0, ha="right",
                  va="center", fontsize=SIZE_SMALL)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    by_pos: dict = defaultdict(list)
    for al, p in corr3p_by_aligner.items():
        if p is not None and win_start <= p < win_end:
            by_pos[p].append(al)
    if not by_pos:
        ax.set_ylim(0, 1)
        return
    max_count = max(len(v) for v in by_pos.values())
    ax.set_ylim(0, max(2, max_count) + 0.4)
    # Cell wide enough to hold the 3-char abbreviation at this window scale
    # (1 ref column is too narrow on a ~100 bp window).
    cell_w = max(2.0, n * 0.022)
    for p, als in by_pos.items():
        x = p - win_start + 0.5
        # winner first (bottom of the stack), then the rest alphabetically
        ordered = sorted(als, key=lambda a: (a != winner_name, a))
        ax.plot([x, x], [0, len(ordered)], color="#cfd8dc", linewidth=0.8, zorder=1)
        for i, al in enumerate(ordered):
            is_win = (al == winner_name)
            col = "#2e7d32" if is_win else "#607d8b"
            ax.add_patch(Rectangle((x - cell_w / 2, i + 0.08), cell_w, 0.84,
                                   facecolor=col, edgecolor="white",
                                   linewidth=0.5, zorder=2, clip_on=False))
            ax.text(x, i + 0.5, _ALIGNER_ABBREV.get(al, al[:3]),
                    ha="center", va="center", fontsize=SIZE_SMALL,
                    color="white", fontweight="bold", zorder=3, clip_on=False)


def render_axis_ticks(ax, win_start: int, win_end: int,
                      right_pad_frac: float = 0.0,
                      suppress_last_tick: bool = False,
                      is_reverse: bool = False):
    """Render position labels along the bottom of a panel.

    Tick + label are placed at the CENTER of the ref column for that
    position (x = pos - win_start + 0.5). In 0-based half-open coords,
    position N occupies the column [N, N+1), so its base character sits
    at column center N+0.5. Aligning the tick to the center makes it
    visually point at the labeled base.

    When the parent panel extends xlim past the data range (e.g. the
    overview's right pad for ED columns), pass ``right_pad_frac`` so
    this tick row uses the same xlim and stays aligned. Pass
    ``suppress_last_tick=True`` to drop the rightmost coord label —
    useful when the right edge bleeds into the ED column area above.
    """
    n = win_end - win_start
    # On minus-strand panels the parent renders with xlim extended to the
    # LEFT of the bars (-n*pad, n) so the ED-column padding ends up on
    # display-RIGHT after the later invert_xaxis(). Mirror that here so
    # tick positions stay aligned with the overview above.
    if is_reverse:
        ax.set_xlim(-n * right_pad_frac, n)
    else:
        ax.set_xlim(0, n * (1.0 + right_pad_frac))
    ax.set_ylim(0, 1.0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_yticks([])
    ax.set_xticks([])
    # Target ~6 labels per panel; round step up to multiple of 5 for readability
    raw_step = max(10, n // 6)
    step = max(10, (raw_step // 5) * 5 or raw_step)
    tick_xs = list(range(0, n + 1, step))
    if suppress_last_tick and tick_xs:
        # suppress_last_tick drops the tick that lands closest to the
        # ED-column area. On plus strand, the ED column is on display-
        # RIGHT and the highest data tick is rightmost (drop tail). On
        # minus strand, the axis is later inverted so data 0 ends up
        # near the ED column on display-RIGHT (drop head instead).
        tick_xs = tick_xs[1:] if is_reverse else tick_xs[:-1]
    for x in tick_xs:
        pos = win_start + x
        cx = x + 0.5  # center of the ref column for `pos`
        # Tick mark hugs the top; text sits at the bottom with a small gap.
        # The widened tick panel (height_ratio 0.55 main / 0.50 overview)
        # gives the gap visible breathing room.
        ax.plot([cx, cx], [0.78, 0.98], color="#666", linewidth=0.8)
        ax.text(cx, 0.00, f"{pos:,}", ha="center", va="bottom",
                fontsize=SIZE_TEXT, color="#333")


# ---------------------------------------------------------------------------
# Junction-zoom rendering
# ---------------------------------------------------------------------------

ZOOM_EXON_FLANK = 50   # bp of exon shown on each side of the junction
ZOOM_INTRON_STUB = 25  # bp of intron shown on each side of the break
ZOOM_GAP_VISUAL = 8    # visual columns used for the "≈ intron NNN bp" label
ZOOM_TOTAL_COLS = (
    ZOOM_EXON_FLANK + ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL
    + ZOOM_INTRON_STUB + ZOOM_EXON_FLANK
)  # 50 + 25 + 8 + 25 + 50 = 158


def _zoom_visual_x(ref_pos: int, donor: int, acceptor: int) -> Optional[float]:
    """Reference position → visual column in the junction-zoom panel.

    Visual layout (158 cols):
       0..49     upstream exon end  (ref [donor-50, donor))
      50..74     intron stub donor  (ref [donor, donor+25))
      75..82     gap region — intron length label, no bases drawn
      83..107    intron stub acceptor (ref [acceptor-25, acceptor))
     108..157    downstream exon start (ref [acceptor, acceptor+50))

    Returns None for ref positions inside the crushed intron interior.
    """
    if donor - ZOOM_EXON_FLANK <= ref_pos < donor:
        return float(ref_pos - (donor - ZOOM_EXON_FLANK))
    if donor <= ref_pos < donor + ZOOM_INTRON_STUB:
        return float(ZOOM_EXON_FLANK + (ref_pos - donor))
    if acceptor - ZOOM_INTRON_STUB <= ref_pos < acceptor:
        base = ZOOM_EXON_FLANK + ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL
        return float(base + (ref_pos - (acceptor - ZOOM_INTRON_STUB)))
    if acceptor <= ref_pos < acceptor + ZOOM_EXON_FLANK:
        base = ZOOM_EXON_FLANK + 2 * ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL
        return float(base + (ref_pos - acceptor))
    return None


def _zoom_ref_positions(donor: int, acceptor: int) -> list[tuple[int, float]]:
    """All (ref_pos, visual_col) pairs that get drawn in the zoom panel."""
    out: list[tuple[int, float]] = []
    for rp in range(donor - ZOOM_EXON_FLANK, donor + ZOOM_INTRON_STUB):
        v = _zoom_visual_x(rp, donor, acceptor)
        if v is not None:
            out.append((rp, v))
    for rp in range(acceptor - ZOOM_INTRON_STUB, acceptor + ZOOM_EXON_FLANK):
        v = _zoom_visual_x(rp, donor, acceptor)
        if v is not None:
            out.append((rp, v))
    return out


def _collect_corrected_junctions(bam_path: Path, qname: str) -> list[tuple[int, int]]:
    """Return [(donor, acceptor), ...] for every N-op in the read's corrected
    CIGAR. donor = first intronic ref pos; acceptor = first post-intronic ref
    pos. The interval is half-open [donor, acceptor)."""
    read = fetch_read(bam_path, qname)
    if read is None or not read.cigartuples:
        return []
    junctions: list[tuple[int, int]] = []
    rp = read.reference_start
    for op, length in read.cigartuples:
        if op == 4 or op == 5:  # soft / hard clip — don't advance ref
            continue
        if op == 3:  # N — junction
            junctions.append((rp, rp + length))
            rp += length
        elif op in (0, 7, 8, 2):  # M, =, X, D advance ref
            rp += length
        # op 1 (I) doesn't advance ref
    return junctions


# ---------------------------------------------------------------------------
# End-correction + Winner-selection panel
# ---------------------------------------------------------------------------

DEFAULT_ALIGNERS = ['minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA']


def _coalesce_adjacent_same_op(
    cigartuples: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Merge adjacent CIGAR ops with the same op code into a single op.
    Normalizes encoding so e.g. (X,1)(X,1) becomes (X,2). Used to make
    grouping by CIGAR insensitive to per-aligner op-splitting style."""
    out: list[tuple[int, int]] = []
    for op, length in cigartuples:
        if out and out[-1][0] == op:
            out[-1] = (op, out[-1][1] + length)
        else:
            out.append((op, length))
    return out


def _split_m_into_eqx(
    cigartuples: list[tuple[int, int]],
    ref_start: int,
    chrom_seq: Optional[str],
    query_seq: Optional[str],
) -> list[tuple[int, int]]:
    """Expand `M` ops (CIGAR op==0, "match or mismatch") into runs of
    `=` (op==7) and `X` (op==8) by comparing read vs. reference base-by-base,
    then coalesce adjacent same-op runs so the result is encoding-canonical.

    Aligners like minimap2 default to emitting `M` for both matches and
    mismatches, which makes the painted-CIGAR overview show large blue
    blocks that hide per-base disagreement. Re-expanding here normalizes
    visual output across aligners (mapPacBio/uLTRA already emit =/X).

    The coalesce step merges things like (X,1)(X,1) → (X,2) so that a
    natively-split CIGAR like gapmm2's `1X1X` ends up identical to the
    expanded-from-M version `2X` — important for the identical-CIGAR
    grouping in the overview.

    Falls back to the original cigartuples when genome or query sequence
    is unavailable.
    """
    if not cigartuples or chrom_seq is None or not query_seq:
        return _coalesce_adjacent_same_op(list(cigartuples)) if cigartuples else []
    out: list[tuple[int, int]] = []
    rp = ref_start
    qp = 0
    for op, length in cigartuples:
        if op == 4:  # S
            qp += length
            out.append((op, length))
        elif op == 5:  # H
            out.append((op, length))
        elif op == 0:  # M — expand
            cur_op: Optional[int] = None
            cur_len = 0
            for i in range(length):
                if rp + i >= len(chrom_seq) or qp + i >= len(query_seq):
                    # Beyond ref or query — treat as = (no useful info).
                    this_op = 7
                else:
                    rb = chrom_seq[rp + i].upper()
                    qb = query_seq[qp + i].upper()
                    # BAM spec: '=' in SEQ means "matches reference at
                    # this position". Raw aligner BAMs (not processed by
                    # rectify's =-SEQ decoder) carry this shorthand.
                    if qb == '=':
                        qb = rb
                    this_op = 7 if rb == qb else 8
                if cur_op is None:
                    cur_op, cur_len = this_op, 1
                elif this_op == cur_op:
                    cur_len += 1
                else:
                    out.append((cur_op, cur_len))
                    cur_op, cur_len = this_op, 1
            if cur_op is not None:
                out.append((cur_op, cur_len))
            rp += length
            qp += length
        elif op in (7, 8):  # =, X already split
            out.append((op, length))
            rp += length
            qp += length
        elif op == 1:  # I
            qp += length
            out.append((op, length))
        elif op == 2:  # D
            rp += length
            out.append((op, length))
        elif op == 3:  # N
            rp += length
            out.append((op, length))
        else:
            out.append((op, length))
    return _coalesce_adjacent_same_op(out)


# Display-only ED weighting: a soft/hard clip is a principled abstention
# ("I don't know how to align these bases"), whereas a forced X/I/D is a
# confident wrong call. Charging the clip LESS than a mismatch keeps an
# over-aggressive aligner (e.g. mapPacBio force-aligning where others honestly
# soft-clip) from looking competitive in the review ED columns. This affects
# ONLY the rendered "raw ED" review numbers — production winner selection lives
# in rectify.core.consensus.scoring.score_alignment and is unchanged.
_SOFTCLIP_ED_COST = 0.5

def _cigar_raw_edit_distance(read, chrom_seq: Optional[str]) -> float:
    """Compute the raw (flat-penalty) edit distance for a corrected read.

    Op costs:
      =          : 0
      X          : 1 per base
      M          : 0 if read base == genome base (when both available), else 1
      D          : 1 per base
      I          : 1 per base
      N (intron) : 0
      S, H       : ``_SOFTCLIP_ED_COST`` (0.5) per base — an honest abstention is
                   charged less than a forced mismatch (display-only asymmetry;
                   see module note above ``_SOFTCLIP_ED_COST``).
    """
    if read is None or read.cigartuples is None:
        return 0.0
    qseq = read.query_sequence
    ref_pos = read.reference_start
    q_pos = 0
    total = 0.0
    for op, length in read.cigartuples:
        if op == 7:                # =
            ref_pos += length; q_pos += length
        elif op == 8:              # X
            total += length
            ref_pos += length; q_pos += length
        elif op == 0:              # M
            if chrom_seq and qseq:
                ref_chunk = chrom_seq[ref_pos:ref_pos + length].upper()
                q_chunk = qseq[q_pos:q_pos + length].upper()
                total += sum(rb != qb for rb, qb in zip(ref_chunk, q_chunk))
            ref_pos += length; q_pos += length
        elif op == 1:              # I
            total += length
            q_pos += length
        elif op == 2:              # D
            total += length
            ref_pos += length
        elif op == 3:              # N
            ref_pos += length
        elif op == 4:              # S
            total += length * _SOFTCLIP_ED_COST
            q_pos += length
        elif op == 5:              # H
            total += length * _SOFTCLIP_ED_COST
        # P (6): no advance
    return total


def _parse_junction_str(s: str) -> list[tuple[int, int]]:
    """Parse a ';'-separated 'donor-acceptor' string into (donor, acceptor)
    tuples. Empty / '0' / missing → []."""
    if not s or s == '0':
        return []
    out: list[tuple[int, int]] = []
    for part in s.split(';'):
        part = part.strip()
        if not part:
            continue
        try:
            d_str, a_str = part.split('-')
            out.append((int(d_str), int(a_str)))
        except ValueError:
            continue
    return out


def _read_uncorrected_meta(
    bam_path: Path, qname: str
) -> Optional[dict]:
    """Pull RNA-direction 5'/3' coords + N-op junctions from an uncorrected
    aligner BAM. Returns None if the read is absent."""
    if not bam_path.exists():
        return None
    read = fetch_read(bam_path, qname)
    if read is None or read.cigartuples is None:
        return None
    rs = read.reference_start
    re = read.reference_end
    if read.is_reverse:
        orig_5p, orig_3p = re - 1, rs
    else:
        orig_5p, orig_3p = rs, re - 1
    juncs: list[tuple[int, int]] = []
    rp = rs
    for op, length in read.cigartuples:
        if op == 4 or op == 5:
            continue
        if op == 3:
            juncs.append((rp, rp + length))
            rp += length
        elif op in (0, 7, 8, 2):
            rp += length
    return {
        'is_reverse': read.is_reverse,
        'ref_start': rs, 'ref_end': re,
        'orig_5p': orig_5p, 'orig_3p': orig_3p,
        'orig_juncs': juncs,
    }


def load_correction_panel_data(
    qname: str,
    summary_tsv: Path,
    aligner_bam_dir: Path,
    aligners: Optional[list[str]] = None,
    per_aligner_corrected_dir: Optional[Path] = None,
    chrom_seq: Optional[str] = None,
    genome_fa: Optional[Path] = None,
) -> Optional[dict]:
    """Build the per-aligner panel data for `qname` from the comparison
    summary TSV + per-aligner uncorrected BAMs. Returns None if the read
    isn't found in the summary.

    Returned dict shape:
      {
        'strand': '+' | '-',
        'rows':   list[ per-aligner dict ],
        'n_introns':            int,
        'winner_introns_corr':  list[(donor, acceptor)],
      }
    """
    if aligners is None:
        aligners = DEFAULT_ALIGNERS
    if not summary_tsv or not Path(summary_tsv).exists():
        return None
    import csv
    rows_for_read: list[dict] = []
    with open(summary_tsv) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('read_id') == qname:
                rows_for_read.append(row)
    if not rows_for_read:
        return None

    strand: Optional[str] = None
    rows: list[dict] = []
    winner_introns_corr: list[tuple[int, int]] = []

    def _is_truthy(v: Optional[str]) -> bool:
        return v in ('1', 'True', 'true')

    for aligner in aligners:
        srow = next((r for r in rows_for_read if r.get('_aligner') == aligner), None)
        if srow is None:
            rows.append({'aligner': aligner, 'missing': True})
            continue
        bam_path = aligner_bam_dir / f'validation_reads.{aligner}.bam'
        meta = _read_uncorrected_meta(bam_path, qname)
        if meta is None:
            rows.append({'aligner': aligner, 'missing': True})
            continue
        if strand is None:
            strand = '-' if meta['is_reverse'] else '+'

        corr_5p_raw = srow.get('five_prime_position', '').strip()
        corr_3p_raw = srow.get('corrected_3prime', '').strip()
        corr_5p = int(corr_5p_raw) if corr_5p_raw and corr_5p_raw.lstrip('-').isdigit() else None
        corr_3p = int(corr_3p_raw) if corr_3p_raw and corr_3p_raw.lstrip('-').isdigit() else None
        corr_juncs = _parse_junction_str(srow.get('junctions', ''))

        # Sort junctions in RNA-direction (5'→3'). For + strand that's
        # ascending genomic; for − strand it's descending genomic.
        rna_reverse = meta['is_reverse']
        orig_juncs_sorted = sorted(meta['orig_juncs'], reverse=rna_reverse)
        corr_juncs_sorted = sorted(corr_juncs, reverse=rna_reverse)

        is_winner = _is_truthy(srow.get('_is_winner'))
        if is_winner:
            winner_introns_corr = corr_juncs_sorted

        # Raw (flat-penalty) ED computed from this aligner's CORRECTED BAM.
        # Each aligner may land the read on a different chromosome
        # (paralogous loci) so we load each read's actual chrom on
        # demand rather than reusing a single passed-in chrom_seq —
        # otherwise cross-chrom aligners get scored against the wrong
        # reference, inflating their raw ED with fake mismatches.
        raw_ed: Optional[float] = None
        if per_aligner_corrected_dir is not None:
            corr_bam_path = Path(per_aligner_corrected_dir) / f"{aligner}.trimmed.bam"
            if corr_bam_path.exists():
                corrected_read = fetch_read(corr_bam_path, qname)
                if corrected_read is not None:
                    this_chrom = corrected_read.reference_name
                    seq_for_ed: Optional[str] = None
                    if genome_fa is not None:
                        try:
                            with pysam.FastaFile(str(genome_fa)) as _fa_ed:
                                seq_for_ed = _fa_ed.fetch(this_chrom).upper()
                        except Exception:
                            seq_for_ed = None
                    if seq_for_ed is None:
                        # Caller didn't pass genome_fa — fall back to the
                        # passed chrom_seq (may be wrong for cross-chrom).
                        seq_for_ed = chrom_seq
                    raw_ed = _cigar_raw_edit_distance(corrected_read, seq_for_ed)

        # _chimera_ok value semantics in the summary TSV: 0 = plausible (good),
        # 1 = suspicious. Display flag is "ok" = (value == 0).
        chim_raw = srow.get('_chimera_ok', '0').strip()
        chimera_ok = chim_raw in ('0', '', 'False', 'false')

        try:
            hp_ed = float(srow.get('hp_edit_distance', 'inf'))
        except (TypeError, ValueError):
            hp_ed = float('inf')
        try:
            aligned_bases = int(srow.get('aligned_bases', '0'))
        except (TypeError, ValueError):
            aligned_bases = 0

        rows.append({
            'aligner': aligner,
            'missing': False,
            'orig_5p': meta['orig_5p'], 'corr_5p': corr_5p,
            'orig_3p': meta['orig_3p'], 'corr_3p': corr_3p,
            'orig_juncs': orig_juncs_sorted,
            'corr_juncs': corr_juncs_sorted,
            'hp_ed': hp_ed,
            'raw_ed': raw_ed,
            'aligned_bases': aligned_bases,
            'chimera_ok': chimera_ok,
            'five_rescued': _is_truthy(srow.get('_five_rescued')),
            'is_winner': is_winner,
            'correction_applied': srow.get('correction_applied', ''),
            'span': meta['ref_end'] - meta['ref_start'],
        })

    n_introns = max(
        (len(r['corr_juncs']) for r in rows if not r.get('missing')),
        default=0,
    )

    # Fallback for intron-column structure: when the winner has empty
    # corr_juncs (e.g., when summary TSV serializes 5'-rescued aligners'
    # junctions as the count "1" instead of donor-acceptor coords), use
    # the longest non-empty corr_juncs from any aligner so the intron
    # columns render against real coordinates. Debugger ticket tracks
    # the TSV-writer fix; this keeps the plotter robust to the bug.
    if not winner_introns_corr and n_introns > 0:
        for r in rows:
            if not r.get('missing') and len(r.get('corr_juncs', [])) == n_introns:
                winner_introns_corr = r['corr_juncs']
                break

    # Effective-utility clustering: group aligners by (corr_5p, corr_3p,
    # frozenset(corr_juncs)). Aligners in the same cluster carry the same
    # biological answer (interchangeable for RECTIFY); their HP-ED differs
    # only in CIGAR-level mismatch/indel placement. Convention: the winner's
    # cluster is always labeled 'A'; other clusters get B, C, ... in order
    # of first appearance.
    def _cluster_key(r: dict) -> Optional[tuple]:
        if r.get('missing'):
            return None
        return (
            r.get('corr_5p'), r.get('corr_3p'),
            tuple(r.get('corr_juncs', [])),
        )

    winner_key = next(
        (_cluster_key(r) for r in rows if r.get('is_winner') and not r.get('missing')),
        None,
    )
    cluster_letter: dict[tuple, str] = {}
    if winner_key is not None:
        cluster_letter[winner_key] = 'A'
    next_letter = ord('B') if winner_key is not None else ord('A')
    for r in rows:
        k = _cluster_key(r)
        if k is None or k in cluster_letter:
            continue
        cluster_letter[k] = chr(next_letter)
        next_letter += 1
    for r in rows:
        k = _cluster_key(r)
        r['eff_cluster'] = cluster_letter.get(k, '—') if k is not None else '—'

    return {
        'strand': strand or '+',
        'rows': rows,
        'n_introns': n_introns,
        'winner_introns_corr': winner_introns_corr,
    }


def _match_intron(
    aligner_juncs: list[tuple[int, int]],
    target_donor: int, target_acceptor: int,
    tolerance: int = 25,
) -> Optional[tuple[int, int]]:
    """Find the aligner's junction whose donor and acceptor are each within
    `tolerance` bp of the target. Returns (donor, acceptor) or None."""
    for d, a in aligner_juncs:
        if abs(d - target_donor) <= tolerance and abs(a - target_acceptor) <= tolerance:
            return (d, a)
    return None


def render_correction_panel(
    ax,
    panel_data: dict,
    qname: str,
    gene_id: str = '',
):
    """Draw the End-correction + Winner-selection panel on `ax`.

    One header strip + one row per aligner. Columns (left→right):
      aligner | 5' orig→corr Δ | [intron-k orig→corr Δd Δa for each k] |
      3' orig→corr Δ | HP-ED | span | chim | pick

    Coords are shown as raw genomic; the header indicates RNA strand so the
    reader can mentally orient. Δ is RNA-direction-signed: positive when the
    corrected end moved TOWARD the RNA 3' end relative to the original.
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    strand = panel_data['strand']
    rows = panel_data['rows']
    n_introns = panel_data['n_introns']
    winner_introns = panel_data['winner_introns_corr']

    # Column widths (chars). Tuned for yeast 7-digit coords.
    W_ALIGNER = 10
    W_COORD   = 7
    W_DELTA   = 5
    W_SPAN    = 13   # "1234567-1234567"
    W_HPED    = 7
    W_RAWED   = 7    # raw (flat-penalty) ED
    W_SPANBP  = 5
    W_CHIM    = 4
    W_EFF     = 3    # 'A' / 'B' / '—'
    W_PICK    = 22

    def _delta_str(orig: Optional[int], corr: Optional[int]) -> str:
        if orig is None or corr is None:
            return '   —'
        # RNA-direction sign: + strand → genomic delta; − strand → flip.
        d = (corr - orig) if strand == '+' else (orig - corr)
        return f"{d:+d}" if d != 0 else "  0"

    def _coord_str(v: Optional[int]) -> str:
        return f"{v}" if v is not None else "  —"

    # ---- Build header text ----
    def _pad_center(s: str, w: int) -> str:
        if len(s) >= w:
            return s[:w]
        left = (w - len(s)) // 2
        right = w - len(s) - left
        return ' ' * left + s + ' ' * right

    # First header row: section labels (5' / intron k / 3' / selection)
    h1_parts = [' ' * W_ALIGNER]
    # 5' end spans 3 columns of widths W_COORD W_COORD W_DELTA + 2 separator spaces
    w_5 = W_COORD * 2 + W_DELTA + 2
    h1_parts.append(_pad_center("── 5′ end ──", w_5))
    for k in range(n_introns):
        w_i = W_SPAN * 2 + W_DELTA * 2 + 3  # 2 spans + 2 deltas + 3 seps
        h1_parts.append(_pad_center(f"── intron {k+1} ──", w_i))
    w_3 = W_COORD * 2 + W_DELTA + 2
    h1_parts.append(_pad_center("── 3′ end ──", w_3))
    w_sel = W_HPED + W_RAWED + W_SPANBP + W_CHIM + W_EFF + W_PICK + 5
    h1_parts.append(_pad_center("── selection ──", w_sel))
    header1 = '  '.join(h1_parts)

    # Sub-header: column names within each section
    h2_parts = [f"{'aligner':<{W_ALIGNER}}"]
    h2_parts.append(f"{'orig':>{W_COORD}} {'corr':>{W_COORD}} {'Δ':>{W_DELTA}}")
    for _ in range(n_introns):
        h2_parts.append(
            f"{'orig span':>{W_SPAN}} {'corr span':>{W_SPAN}} {'Δd':>{W_DELTA}} {'Δa':>{W_DELTA}}"
        )
    h2_parts.append(f"{'orig':>{W_COORD}} {'corr':>{W_COORD}} {'Δ':>{W_DELTA}}")
    h2_parts.append(
        f"{'HP-ED':>{W_HPED}} {'raw ED':>{W_RAWED}} "
        f"{'span':>{W_SPANBP}} {'chim':>{W_CHIM}} "
        f"{'eff':>{W_EFF}} {'pick':<{W_PICK}}"
    )
    header2 = '  '.join(h2_parts)

    # ---- Build per-aligner rows ----
    body_lines: list[str] = []
    for row in rows:
        aligner = row['aligner']
        if row.get('missing'):
            line_parts = [f"{aligner:<{W_ALIGNER}}"]
            # missing-row placeholder
            line_parts.append(f"{'—':>{W_COORD}} {'—':>{W_COORD}} {'—':>{W_DELTA}}")
            for _ in range(n_introns):
                line_parts.append(
                    f"{'—':>{W_SPAN}} {'—':>{W_SPAN}} {'—':>{W_DELTA}} {'—':>{W_DELTA}}"
                )
            line_parts.append(f"{'—':>{W_COORD}} {'—':>{W_COORD}} {'—':>{W_DELTA}}")
            line_parts.append(
                f"{'—':>{W_HPED}} {'—':>{W_RAWED}} "
                f"{'—':>{W_SPANBP}} {'—':>{W_CHIM}} "
                f"{'—':>{W_EFF}} {'(no alignment)':<{W_PICK}}"
            )
            body_lines.append('  '.join(line_parts))
            continue

        line_parts = [f"{aligner:<{W_ALIGNER}}"]
        # 5' end block
        line_parts.append(
            f"{_coord_str(row['orig_5p']):>{W_COORD}} "
            f"{_coord_str(row['corr_5p']):>{W_COORD}} "
            f"{_delta_str(row['orig_5p'], row['corr_5p']):>{W_DELTA}}"
        )
        # Intron blocks (matched to winner's intron order)
        any_2h_refined = False
        for k, (w_d, w_a) in enumerate(winner_introns):
            orig = _match_intron(row['orig_juncs'], w_d, w_a)
            corr = _match_intron(row['corr_juncs'], w_d, w_a)
            orig_str = f"{orig[0]}-{orig[1]}" if orig else "—"
            corr_str = f"{corr[0]}-{corr[1]}" if corr else "—"
            if orig and corr:
                dd = (corr[0] - orig[0]) if strand == '+' else (orig[0] - corr[0])
                da = (corr[1] - orig[1]) if strand == '+' else (orig[1] - corr[1])
                if dd != 0 or da != 0:
                    any_2h_refined = True
                dd_str = f"{dd:+d}" if dd != 0 else " 0"
                da_str = f"{da:+d}" if da != 0 else " 0"
            else:
                dd_str = da_str = "—"
            line_parts.append(
                f"{orig_str:>{W_SPAN}} {corr_str:>{W_SPAN}} "
                f"{dd_str:>{W_DELTA}} {da_str:>{W_DELTA}}"
            )
        # 3' end block
        line_parts.append(
            f"{_coord_str(row['orig_3p']):>{W_COORD}} "
            f"{_coord_str(row['corr_3p']):>{W_COORD}} "
            f"{_delta_str(row['orig_3p'], row['corr_3p']):>{W_DELTA}}"
        )
        # Selection block
        hp_ed = row['hp_ed']
        if hp_ed == float('inf'):
            hp_str = '   inf'
        else:
            hp_str = f"{hp_ed:.1f}"
        chim_str = '✓' if row['chimera_ok'] else '✗'
        # Pick column: winner first, then badges (refined / 5' rescued / chim filt)
        pick_bits: list[str] = []
        if row['is_winner']:
            pick_bits.append('← winner')
        if not row['chimera_ok']:
            pick_bits.append('chim. filt.')
        if row['five_rescued']:
            pick_bits.append("5′ rescued")
        if any_2h_refined:
            pick_bits.append('refined (2H)')
        pick = ', '.join(pick_bits) if pick_bits else ''

        eff_str = row.get('eff_cluster', '—')
        raw_ed_val = row.get('raw_ed')
        raw_ed_str = f"{raw_ed_val:.1f}" if raw_ed_val is not None else '—'
        line_parts.append(
            f"{hp_str:>{W_HPED}} "
            f"{raw_ed_str:>{W_RAWED}} "
            f"{row['span']:>{W_SPANBP}} "
            f"{chim_str:>{W_CHIM}} "
            f"{eff_str:>{W_EFF}} "
            f"{pick:<{W_PICK}}"
        )
        body_lines.append('  '.join(line_parts))

    # ---- Top meta line ----
    n_intron_label = (
        f"{n_introns} intron{'s' if n_introns != 1 else ''}"
        if n_introns > 0 else "no introns"
    )
    gene_chunk = f"  {gene_id}" if gene_id else ""
    meta_line = f"strand {strand}{gene_chunk}    {n_intron_label}"

    # ---- Composite text ----
    # Winner CIGAR is now painted onto the overview track (see render_overview)
    # so we no longer need the raw-string footer here.
    full_text = '\n'.join(
        [meta_line, '', header1, header2] + body_lines
    )

    # Render. Use top-anchored text with small top margin so the text fills
    # the box without ascender clipping.
    ax.text(
        0.005, 0.92, full_text,
        family='monospace', fontsize=7.5,
        va='top', ha='left',
        transform=ax.transAxes,
    )
    # Subtle box around the panel for visual separation
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch(
        (0.0, 0.0), 1.0, 1.0,
        boxstyle="round,pad=0.005,rounding_size=0.01",
        transform=ax.transAxes,
        linewidth=0.6, edgecolor='#888888', facecolor='#FAFAFA',
        zorder=-1,
    ))


def render_zoom_ref_row(ax, ref_provider, donor: int, acceptor: int,
                        is_reverse: bool = False, label: str = "ref"):
    """Render the reference row of a junction-zoom panel.

    ``ref_provider`` is a callable ``ref_pos -> base char`` (or None when out
    of range). Intron-stub bases get a slightly different background + italic
    glyph; exon bases use the same look as the existing 3'-end ref row.
    Donor and acceptor splice motifs are labeled with ↓ markers above.
    """
    ax.set_xlim(0, ZOOM_TOTAL_COLS)
    ax.set_ylim(0, 1.45)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center",
                  fontsize=SIZE_TEXT, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    for rp, vx in _zoom_ref_positions(donor, acceptor):
        b = ref_provider(rp)
        if b is None:
            continue
        in_intron = (donor <= rp < donor + ZOOM_INTRON_STUB) or \
                    (acceptor - ZOOM_INTRON_STUB <= rp < acceptor)
        bg = "#EEEEF6" if in_intron else "#F5F5F5"
        b_disp = _maybe_comp(b, is_reverse)
        fg = BASE_COLOR.get(b_disp.upper(), "#000")
        ax.add_patch(Rectangle((vx, 0.05), 1, 0.50, facecolor=bg, edgecolor="none"))
        ax.text(vx + 0.5, 0.30, b_disp.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg,
                fontstyle="italic" if in_intron else "normal",
                fontweight="bold")
    # Donor (5' SS) / acceptor (3' SS) markers in RNA terms.
    #
    # On the plus strand, the RNA donor sits at the genomic donor coordinate
    # (lower coord) and the RNA acceptor at the genomic acceptor coordinate
    # (higher coord). On the minus strand, RNA reads 3' → 5' along the plus
    # genome, so the RNA donor is at the higher genomic coord and the RNA
    # acceptor at the lower genomic coord. The visual ax.invert_xaxis() at the
    # end of render() flips the panel for minus strand, so we want the
    # genomic→visual mapping to place the RNA donor on the (RNA-)5' side =
    # visual left after inversion.
    genomic_donor_vx    = ZOOM_EXON_FLANK + 0.0
    genomic_acceptor_vx = ZOOM_EXON_FLANK + 2 * ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL
    if is_reverse:
        rna_donor_vx     = genomic_acceptor_vx
        rna_acceptor_vx  = genomic_donor_vx
        # RNA donor motif = reverse-complement of plus[acceptor-2:acceptor]
        rna_donor_motif    = _comp(ref_provider(acceptor - 1) or "?") + \
                             _comp(ref_provider(acceptor - 2) or "?")
        # RNA acceptor motif = reverse-complement of plus[donor:donor+2]
        rna_acceptor_motif = _comp(ref_provider(donor + 1) or "?") + \
                             _comp(ref_provider(donor) or "?")
    else:
        rna_donor_vx     = genomic_donor_vx
        rna_acceptor_vx  = genomic_acceptor_vx
        rna_donor_motif    = (ref_provider(donor) or "?") + \
                             (ref_provider(donor + 1) or "?")
        rna_acceptor_motif = (ref_provider(acceptor - 2) or "?") + \
                             (ref_provider(acceptor - 1) or "?")
    ax.plot(rna_donor_vx, 0.95, marker="v", markersize=11, color="#7b1fa2")
    ax.text(rna_donor_vx, 1.22, f"donor {rna_donor_motif.upper()}",
            ha="center", va="bottom", fontsize=SIZE_TEXT,
            color="#7b1fa2", fontweight="bold")
    ax.plot(rna_acceptor_vx, 0.95, marker="v", markersize=11, color="#7b1fa2")
    ax.text(rna_acceptor_vx, 1.22, f"acceptor {rna_acceptor_motif.upper()}",
            ha="center", va="bottom", fontsize=SIZE_TEXT,
            color="#7b1fa2", fontweight="bold")
    # Intron-length label centered in the gap
    gap_center = ZOOM_EXON_FLANK + ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL / 2.0
    ax.text(gap_center, 0.30, f"≈{acceptor - donor:,} bp",
            ha="center", va="center", fontsize=SIZE_TEXT, color="#555",
            fontstyle="italic",
            bbox=dict(facecolor="#fffce0", edgecolor="#d9b738",
                      boxstyle="round,pad=0.2", linewidth=0.5))
    # Dashed separators bracketing the gap
    ax.axvline(ZOOM_EXON_FLANK + ZOOM_INTRON_STUB,
               color="#aaa", linestyle=":", linewidth=0.7)
    ax.axvline(ZOOM_EXON_FLANK + ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL,
               color="#aaa", linestyle=":", linewidth=0.7)


def render_zoom_alignment_row(ax, aligned_window: dict, donor: int, acceptor: int,
                              label: str, window_start: int,
                              is_reverse: bool = False):
    """Render one alignment row in the junction-zoom panel.

    ``aligned_window`` is the ``walk_cigar`` output for window
    ``[window_start, window_start + (acceptor + ZOOM_EXON_FLANK - window_start))``.
    Bases inside the crushed intron interior are skipped; bases in the visible
    intron stubs (±25 bp from each splice site) are drawn so an N-op shows up
    as missing bases bracketed by intron-stub ref context.
    """
    ax.set_xlim(0, ZOOM_TOTAL_COLS)
    # ylim extended ABOVE the row so insertion pills (purple-bordered boxes
    # with the inserted base inside) can sit above the bar without
    # clipping. Matches the headroom used in render_alignment_row.
    ax.set_ylim(0, 2.60)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(label, rotation=0, ha="right", va="center",
                  fontsize=SIZE_TEXT, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    bases = aligned_window["bases"]
    state = aligned_window["state"]
    mismatches = aligned_window["mismatches"]
    for rp, vx in _zoom_ref_positions(donor, acceptor):
        idx = rp - window_start
        if not (0 <= idx < len(bases)):
            continue
        b = bases[idx]
        st = state[idx]
        if b is None or st is None:
            continue
        if st == "intron":
            # N-op interior bar — visible in the intron stub but suppressed in
            # the crushed center (already absent from _zoom_ref_positions).
            ax.add_patch(Rectangle((vx, 0.10), 1, 0.80,
                                   facecolor="#F0F0F8", edgecolor="none"))
            ax.text(vx + 0.5, 0.50, "|", ha="center", va="center",
                    fontsize=10, family="monospace", color="#666")
            continue
        b_disp = _maybe_comp(b, is_reverse) if b != "-" else b
        bg, fg = _color_for_base(b_disp, st, idx in mismatches)
        ax.add_patch(Rectangle((vx, 0.10), 1, 0.80, facecolor=bg, edgecolor="none"))
        ax.text(vx + 0.5, 0.50, b_disp.upper(), ha="center", va="center",
                fontsize=8, family="monospace", color=fg, fontweight="bold")
    ax.axvline(ZOOM_EXON_FLANK + ZOOM_INTRON_STUB,
               color="#aaa", linestyle=":", linewidth=0.7)
    ax.axvline(ZOOM_EXON_FLANK + ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL,
               color="#aaa", linestyle=":", linewidth=0.7)

    # Insertion pills — same style as render_alignment_row: full-size
    # colored bases wedged in a vivid purple-bordered box, positioned
    # ABOVE the row with a black ▼ pointing down to the exact insertion
    # site. Nearby insertions stagger vertically to avoid pill overlap.
    insertions_sorted = sorted(aligned_window.get("insertions", []),
                               key=lambda t: t[0])
    last_x_end = -1e9
    stagger = 0
    for entry in insertions_sorted:
        if len(entry) >= 3:
            rp, length, bases_str = entry[0], entry[1], entry[2]
        else:
            rp, length = entry[0], entry[1]
            bases_str = ""
        vx = _zoom_visual_x(rp, donor, acceptor)
        if vx is None:
            continue
        char_width = 1.0     # match the normal-cell width so inserted letters
                             # render at the same visual size as body bases
        total_w = length * char_width
        start_x = vx - total_w / 2
        end_x = vx + total_w / 2
        if start_x < last_x_end + 0.3:
            stagger = 1 - stagger
        else:
            stagger = 0
        # Pill sits ABOVE the row (y=1.0 = top of base rectangles). Two
        # vertical tiers for stagger when insertions cluster.
        pill_y = 1.25 if stagger == 0 else 1.85
        pill_h = 0.70    # taller pill so the inserted letter has buffer above/below
        # Black ▼ pointing down to the insertion site
        ax.plot(vx, 1.02, marker="v", markersize=7, color="#000",
                clip_on=False, zorder=5)
        # Thin black guide line from row top up to pill bottom
        ax.plot([vx, vx], [1.05, pill_y - 0.02],
                color="#000", linewidth=0.6, clip_on=False, zorder=4)
        # Light purple highlight (no border, generous padding around the
        # letter so it visually breathes).
        ax.add_patch(Rectangle((start_x - 0.20, pill_y),
                               total_w + 0.40, pill_h,
                               facecolor="#E1BEE7", edgecolor="none",
                               clip_on=False, zorder=3))
        # Render inserted bases when ≤4; otherwise show "+N…" label.
        bases_disp = ("".join(_maybe_comp(c, is_reverse) for c in bases_str)
                      if bases_str else bases_str)
        if length <= 4 and bases_disp:
            for j, b in enumerate(bases_disp):
                bx = start_x + (j + 0.5) * char_width
                fg = BASE_COLOR.get(b.upper(), "#7B1FA2")
                ax.text(bx, pill_y + pill_h / 2, b.upper(),
                        ha="center", va="center",
                        fontsize=8.5, family="monospace",
                        color=fg, fontweight="bold",
                        clip_on=False, zorder=4)
        else:
            label_text = f"+{length}"
            if bases_disp:
                label_text += (f" {bases_disp[:3]}…" if length > 3
                               else f" {bases_disp}")
            ax.text(vx, pill_y + pill_h / 2, label_text,
                    ha="center", va="center",
                    fontsize=SIZE_SMALL, family="monospace",
                    color="#7B1FA2", fontweight="bold",
                    clip_on=False, zorder=4)
        last_x_end = max(last_x_end, end_x)


def render_zoom_ticks(ax, donor: int, acceptor: int):
    """Position ticks at the six landmarks of the broken axis."""
    ax.set_xlim(0, ZOOM_TOTAL_COLS)
    ax.set_ylim(0, 1.0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    landmarks = [
        (0,                                           donor - ZOOM_EXON_FLANK),
        (ZOOM_EXON_FLANK,                             donor),
        (ZOOM_EXON_FLANK + ZOOM_INTRON_STUB,          donor + ZOOM_INTRON_STUB),
        (ZOOM_EXON_FLANK + ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL,
         acceptor - ZOOM_INTRON_STUB),
        (ZOOM_EXON_FLANK + 2 * ZOOM_INTRON_STUB + ZOOM_GAP_VISUAL,
         acceptor),
        (ZOOM_TOTAL_COLS,                             acceptor + ZOOM_EXON_FLANK),
    ]
    for vx, refp in landmarks:
        ax.plot([vx, vx], [0.70, 0.95], color="#666", linewidth=0.8)
        ax.text(vx, 0.05, f"{refp:,}", ha="center", va="bottom",
                fontsize=7, color="#333")


def render_overview(ax, aligned_set, win_start: int, win_end: int,
                    orig_3p: Optional[int], corr_3p: Optional[int],
                    aligner_ed_map: Optional[dict] = None,
                    show_column_headers: bool = True,
                    is_reverse: bool = False):
    """Schematic overview of the full alignment span with N-op gaps drawn
    as thin connector lines labeled by intron size. One row per BAM version,
    stacked vertically.

    When ``aligner_ed_map`` is supplied, EER-ED and raw-ED values render
    on the display-RIGHT side of the panel (just past the painted CIGAR)
    with column headers aligned with the "overview" tag at top-left.

    On minus-strand reads the caller flips the axis later via
    ``invert_xaxis()`` so RNA 5'→3' reads left-to-right. To keep the ED
    column padding visible on display-RIGHT after that flip, we extend
    xlim to the LEFT of the bars (``(-n*0.14, n)``) when ``is_reverse``;
    otherwise we extend to the right (``(0, n*1.14)``). Bars themselves
    remain at data ``[0, n]`` in both cases.
    """
    n = win_end - win_start
    ed_right_frac = 0.14 if aligner_ed_map else 0.0
    if is_reverse:
        ax.set_xlim(-n * ed_right_frac, n)
    else:
        ax.set_xlim(0, n * (1.0 + ed_right_frac))
    ax.set_ylim(0, len(aligned_set) + 0.5)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Suppress y-ticks — we'll draw the row labels manually via ax.text
    # in transAxes coords so we can add aligned EER-ED / raw-ED columns
    # alongside them, with column headers at the top.
    ax.set_yticks([])

    # ---- Pre-compute per-row ED values ----
    # Aligners in the same group share the same post-normalization
    # CIGAR by definition (that's what defines the group), so their
    # EER ED / raw ED are equal. Show a single value per row — the
    # group's first aligner — for a much tighter left-margin layout.
    def _eds_for_label(label: str) -> tuple[str, str]:
        if not aligner_ed_map:
            return "", ""
        stripped = (label.split(":", 1)[-1] if ":" in label else label).strip()
        # Strip cross-chrom "[chrX]" tag — present on cross-chrom rows
        # whose label looks like "deSALT  [chrVI]".
        import re as _re
        stripped = _re.sub(r"\s*\[[^\]]*\]\s*$", "", stripped).strip()
        # First aligner in the group only (others are CIGAR-identical).
        first_name = stripped.split(",", 1)[0].strip()
        hp, raw = aligner_ed_map.get(first_name, (None, None))
        eer_s = f"{hp:.1f}" if hp is not None else "—"
        raw_s = f"{raw:.1f}" if raw is not None else "—"
        return eer_s, raw_s

    # Layout: row labels just LEFT of the axis (transAxes coords); ED
    # columns to the RIGHT of the painted CIGAR, also in transAxes so
    # they stay anchored to the display-right regardless of the
    # minus-strand invert_xaxis() flip applied later. The right 14% of
    # the axes was reserved via xlim above for these columns.
    x_lbl_ax = -0.02
    x_eer_ax = 0.94          # right-aligned EER value position (transAxes)
    x_raw_ax = 0.995         # right-aligned raw value position (transAxes)
    have_ed = bool(aligner_ed_map)

    # Background banding removed per style review — cleaner white
    # background. Winner indicated by bold label only (left-edge green
    # accent was dropped per user request — the "winner: ..." text
    # already conveys the same information).
    for k, (label, aligned) in enumerate(aligned_set):
        y_center = len(aligned_set) - k - 0.5
        y_frac = y_center / (len(aligned_set) + 0.5)
        is_winner = isinstance(label, str) and label.startswith("winner")
        ax.text(x_lbl_ax, y_frac, label, transform=ax.transAxes,
                ha="right", va="center", fontsize=SIZE_TEXT,
                fontweight=("bold" if is_winner else "normal"))
        if have_ed:
            eer_s, raw_s = _eds_for_label(label)
            ax.text(x_eer_ax, y_frac, eer_s, transform=ax.transAxes,
                    ha="right", va="center", fontsize=SIZE_TEXT,
                    family="monospace", color="#333")
            ax.text(x_raw_ax, y_frac, raw_s, transform=ax.transAxes,
                    ha="right", va="center", fontsize=SIZE_TEXT,
                    family="monospace", color="#333")

    # ---- Column headers ----
    if show_column_headers:
        ax.text(x_lbl_ax, 1.02, "overview", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=9,
                color="#555", fontweight="bold")
        if have_ed:
            ax.text(x_eer_ax, 1.02, "EER ED", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=9,
                    color="#555", fontweight="bold")
            ax.text(x_raw_ax, 1.02, "raw ED", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=9,
                    color="#555", fontweight="bold")
    # Per-CIGAR-op color palette for the overview (CIGAR-painted body).
    # Subtle hues for the common ops so the eye can scan for X/D/I outliers.
    CIGAR_OV_COLOR = {
        7: "#cfd8dc",  # = (sequence match) — neutral grey
        0: "#FFEB3B",  # M — should not survive _split_m_into_eqx; loud yellow
                       # so any residual ambiguous op is visible as an
                       # anomaly (bundle staleness signal). Dropped from
                       # the legend per user request; the loud color +
                       # absent legend entry triggers a "what's that?"
                       # response if the bundle regresses.
        8: "#e91e63",  # X (mismatch) — pink
        2: "#ff9800",  # D (deletion in read) — orange
        1: "#9c27b0",  # I (insertion in read) — purple (tick, no ref width)
        4: "#bbbbbb",  # S (soft-clip) — faded grey
        5: "#888888",  # H (hard-clip) — darker grey
        # 3 (N) handled as gap connector below.
    }
    CIGAR_OV_LETTER = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 7: "=", 8: "X"}

    # Width threshold for in-segment labels. Skip labels for tiny segments to
    # avoid clutter.
    MIN_LABEL_W = 6

    def _place_labels(label_specs: list[tuple[float, str, str]], y_c: float):
        """Place labels below y_c, redistributing colliding clusters
        laterally so each label has its own footprint. Original op-tick
        positions are connected to displaced labels with a thin leader line
        (color-matched to the segment so the eye can re-associate).

        Width estimate uses a conservative char-width in data units; if
        labels still touch the heuristic can be widened by raising
        ``CHAR_W_FACTOR``.
        """
        if not label_specs:
            return
        CHAR_W_FACTOR = 0.0042   # data-units per char (tuneable)
        char_w = max(0.7, n * CHAR_W_FACTOR)
        pad = char_w * 0.5

        specs = sorted(label_specs, key=lambda s: s[0])
        widths = [len(t) * char_w for _, t, _ in specs]

        # Pass 1: detect overlap clusters using original x-positions.
        clusters: list[list[int]] = [[0]]
        for i in range(1, len(specs)):
            prev_i = clusters[-1][-1]
            prev_right = specs[prev_i][0] + widths[prev_i] / 2
            this_left  = specs[i][0]      - widths[i]      / 2
            if this_left < prev_right + pad:
                clusters[-1].append(i)
            else:
                clusters.append([i])

        # Pass 2: redistribute x positions inside each multi-label cluster
        # so labels sit edge-to-edge with `pad` gaps, centered on the
        # cluster's natural midpoint.
        final_x: list[float] = [0.0] * len(specs)
        for cluster in clusters:
            if len(cluster) == 1:
                i = cluster[0]
                final_x[i] = specs[i][0]
                continue
            total_w = sum(widths[i] for i in cluster) + (len(cluster) - 1) * pad
            leftmost  = specs[cluster[0]][0]  - widths[cluster[0]]  / 2
            rightmost = specs[cluster[-1]][0] + widths[cluster[-1]] / 2
            center = (leftmost + rightmost) / 2
            cur_x = center - total_w / 2
            for i in cluster:
                final_x[i] = cur_x + widths[i] / 2
                cur_x += widths[i] + pad

        # Pass 3: render. Leader line drawn when label is visibly displaced.
        for i, (orig_x, text, color) in enumerate(specs):
            new_x = final_x[i]
            ax.text(new_x, y_c - 0.36, text, ha='center', va='top',
                    fontsize=SIZE_SMALL, color=color)
            if abs(new_x - orig_x) > char_w * 0.3:
                ax.plot([orig_x, new_x],
                        [y_c - 0.21, y_c - 0.33],
                        color=color, linewidth=0.4, alpha=0.55,
                        solid_capstyle="butt", zorder=2)

    # ----------------------------------------------------------------
    # Per-ref-position divergence highlighting
    # ----------------------------------------------------------------
    # Build a per-row vector of (op, ins_len_after) per ref position in
    # [win_start, win_end). Insertions are attributed to the ref position
    # immediately AFTER them (so an insertion between ref_pos K and K+1 is
    # recorded as ins_len_after at K+1). Positions outside any row's
    # alignment get op = None (off-row).
    n_cols = win_end - win_start

    def _row_state_vector(aligned: dict) -> list[tuple[Optional[int], int]]:
        state: list[tuple[Optional[int], int]] = [(None, 0)] * n_cols
        cigar = aligned.get("cigartuples") or []
        rp_seed = aligned.get("ref_start_raw")
        if rp_seed is None:
            return state
        rp = rp_seed
        pending_ins_len = 0
        for op, length in cigar:
            if op == 5:                # H — no advance
                continue
            if op == 4:                # S — no ref advance
                continue
            if op == 1:                # I — accumulates to next ref position
                pending_ins_len += length
                continue
            if op in (0, 7, 8):        # M / = / X — ref + query
                for i in range(length):
                    idx = (rp + i) - win_start
                    if 0 <= idx < n_cols:
                        # Insertion (if any) attributed to first base of run.
                        ins_here = pending_ins_len if i == 0 else 0
                        state[idx] = (op, ins_here)
                pending_ins_len = 0
                rp += length
            elif op == 2:              # D — ref only
                for i in range(length):
                    idx = (rp + i) - win_start
                    if 0 <= idx < n_cols:
                        ins_here = pending_ins_len if i == 0 else 0
                        state[idx] = (op, ins_here)
                pending_ins_len = 0
                rp += length
            elif op == 3:              # N — ref only (intron)
                for i in range(length):
                    idx = (rp + i) - win_start
                    if 0 <= idx < n_cols:
                        state[idx] = (op, 0)
                pending_ins_len = 0
                rp += length
        return state

    state_vectors: list[list[tuple[Optional[int], int]]] = []
    for label, aligned in aligned_set:
        state_vectors.append(_row_state_vector(aligned))

    # Identify divergent ref positions: at least 2 rows have a defined
    # state and not all defined states agree. Off-row (None) doesn't
    # count as divergence by itself (the row simply doesn't reach this
    # position — already visible by the row's bar not extending here).
    divergent = [False] * n_cols
    for idx in range(n_cols):
        defined = [sv[idx] for sv in state_vectors if sv[idx][0] is not None]
        if len(defined) >= 2 and len(set(defined)) > 1:
            divergent[idx] = True

    # Coalesce divergent positions into [start, end) ranges.
    divergent_ranges: list[tuple[int, int]] = []
    i = 0
    while i < n_cols:
        if not divergent[i]:
            i += 1
            continue
        j = i + 1
        while j < n_cols and divergent[j]:
            j += 1
        divergent_ranges.append((i, j))
        i = j

    # Render the shaded bands. Two layers:
    #  - A wide background patch (zorder<bars) at low alpha so the eye can
    #    trace divergence extent through the painted bars.
    #  - A narrow top-edge ribbon at high zorder + full opacity so 1-bp
    #    divergences (a single ref column wide) are still visible at the
    #    overview scale. The ribbon sits just above the top row to avoid
    #    hiding the bars.
    if divergent_ranges and SHOW_DIVERGENCE_SHADING:
        y_lo = 0.0
        y_hi = len(aligned_set) + 0.0
        # Floor: expand a 1-bp range to ~0.4% of the panel width so it has
        # at least a perceptible thickness. The expansion is centered on
        # the original range so the marker stays where it should.
        n_total = win_end - win_start
        min_w = max(1.0, n_total * 0.004)
        for x0, x1 in divergent_ranges:
            w = x1 - x0
            if w < min_w:
                cx = (x0 + x1) / 2
                x0_d, x1_d = cx - min_w / 2, cx + min_w / 2
            else:
                x0_d, x1_d = x0, x1
            # Background band (full-height, low alpha, behind bars).
            ax.add_patch(Rectangle(
                (x0_d, y_lo), x1_d - x0_d, y_hi - y_lo,
                facecolor="#FFD54F", edgecolor="none",
                alpha=0.40, zorder=0,
            ))
            # Top ribbon (thin, full opacity, ABOVE bars) — guarantees that
            # even 1-bp divergences are spotted.
            ribbon_y = y_hi - 0.08
            ribbon_h = 0.10
            ax.add_patch(Rectangle(
                (x0_d, ribbon_y), x1_d - x0_d, ribbon_h,
                facecolor="#F57C00",        # darker amber for definition
                edgecolor="none",
                alpha=0.95, zorder=3,
            ))

    for k, (label, aligned) in enumerate(aligned_set):
        y_center = len(aligned_set) - k - 0.5
        if aligned["aln_start"] is None:
            continue
        a, b = aligned["full_aln_start"], aligned["full_aln_end"]
        if a is None or b is None:
            continue

        # Per-row label collection — placed at the end with anti-collision.
        row_labels: list[tuple[float, str, str]] = []

        # 5' soft-clip — split into aligner-soft-clip (closer to body) and
        # the pA-tail portion (outer end) when this is a minus-strand read
        # (the pA tail lives at the BAM-5' end of minus-strand records).
        sc5 = aligned["five_softclip_bp"]
        if sc5 > 0:
            x0 = max(0, a - win_start)
            x1 = max(0, aligned["aln_start"] - win_start)
            if x1 > x0:
                # Prefer the trim-extended effective length when the
                # per-base walk-back set one; otherwise fall back to the
                # parquet polya_len. Only applies to minus strand (pA
                # lives in the 5' BAM soft-clip there).
                if aligned.get("is_reverse"):
                    _pa = aligned.get("polya_len_effective",
                                       aligned.get("polya_len", 0))
                else:
                    _pa = 0
                _pa = min(_pa, sc5)
                # In ref-position layout the 5' soft-clip runs from
                # `a` (low ref) at i_clip=0 to `aln_start-1` at i_clip=sc5-1.
                # The pA tail is i_clip ∈ [0, _pa) → ref positions [a, a+_pa).
                pa_w = (x0 + _pa) - x0 if _pa > 0 else 0
                if pa_w > 0:
                    ax.add_patch(Rectangle((x0, y_center - 0.18), pa_w, 0.36,
                                           facecolor=SOFTCLIP_PA_BG,
                                           edgecolor=SOFTCLIP_PA_FG,
                                           linewidth=0.5))
                # Non-pA clip [x0+pa_w, x1]: grey (native) + slate-blue
                # (walkback-removed, INNERMOST = nearest the body at x1).
                _rest_x0 = x0 + pa_w
                _rest_w = x1 - _rest_x0
                if _rest_w > 0:
                    _wb = min(aligned.get("ov_wb_n", 0), _rest_w)
                    _grey_w = _rest_w - _wb
                    if _grey_w > 0:
                        ax.add_patch(Rectangle((_rest_x0, y_center - 0.18),
                                               _grey_w, 0.36,
                                               facecolor=SOFTCLIP_BG,
                                               edgecolor=SOFTCLIP_FG,
                                               linewidth=0.5))
                    if _wb > 0:
                        ax.add_patch(Rectangle((_rest_x0 + _grey_w, y_center - 0.18),
                                               _wb, 0.36,
                                               facecolor=SOFTCLIP_WB_BG,
                                               edgecolor=SOFTCLIP_WB_FG,
                                               linewidth=0.5))

        # CIGAR-painted body: per-op colored segments. Labels collected.
        cigar = aligned.get("cigartuples") or []
        ref_start = aligned.get("ref_start_raw")
        rp = ref_start if ref_start is not None else aligned["aln_start"]
        cigar_iter = iter(cigar)
        if cigar and cigar[0][0] == 4:
            next(cigar_iter)
        for op, length in cigar_iter:
            if op == 5:           # H — rendered as outside-body patch below.
                continue
            if op == 4:           # trailing S — handled after the loop.
                break
            if op == 3:           # N — intron gap connector + label (stays ABOVE).
                gx0 = max(0, rp - win_start)
                gx1 = min(n, rp + length - win_start)
                if gx1 > gx0:
                    ax.plot([gx0, gx1], [y_center, y_center],
                            color="#aaa", linewidth=0.8)
                    ax.text((gx0 + gx1) / 2, y_center + 0.22,
                            f"{length} bp", ha="center", va="bottom",
                            fontsize=SIZE_SMALL, color="#666")
                rp += length
                continue
            if op == 1:           # I — no ref advance; draw a thin tick.
                x = rp - win_start
                if 0 <= x < n:
                    ax.plot([x, x], [y_center - 0.20, y_center + 0.20],
                            color=CIGAR_OV_COLOR[1], linewidth=1.0,
                            solid_capstyle="butt", zorder=3)
                    # Op-letter dropped; "+" prefix retained to signal "no
                    # ref consumption" (i.e., insertion vs. deletion).
                    row_labels.append((x, f"+{length}", CIGAR_OV_COLOR[1]))
                continue
            # M, =, X, D — ref-consuming.
            color = CIGAR_OV_COLOR.get(op, "#999")
            x0 = max(0, rp - win_start)
            x1 = min(n, rp + length - win_start)
            if x1 > x0:
                ax.add_patch(Rectangle((x0, y_center - 0.20), x1 - x0, 0.40,
                                       facecolor=color, edgecolor="none",
                                       alpha=0.90))
            # CIGAR-op length labels removed — segment color carries the
            # op type and the per-base view below shows lengths exactly.
            # Legend at top of figure documents the color → op mapping.
            rp += length

        # 3' soft-clip — split into aligner-soft-clip (closer to body) and
        # the pA-tail portion (outer end) when this is a plus-strand read
        # (the pA tail lives at the BAM-3' end of plus-strand records).
        sc3 = aligned["three_softclip_bp"]
        if sc3 > 0:
            x0 = max(0, aligned["aln_end"] - win_start)
            x1 = min(n, b - win_start)
            if x1 > x0:
                # Plus strand: pA lives in the 3' BAM soft-clip. Use the
                # trim-extended effective length when set by the per-base
                # walk-back; fall back to the parquet polya_len.
                if not aligned.get("is_reverse"):
                    _pa = aligned.get("polya_len_effective",
                                       aligned.get("polya_len", 0))
                else:
                    _pa = 0
                _pa = min(_pa, sc3)
                # In ref-position layout the 3' soft-clip runs from
                # `aln_end` (i_clip=0) up to `aln_end+sc3-1` (i_clip=sc3-1).
                # The pA tail is i_clip ∈ [sc3-_pa, sc3) → ref positions
                # [aln_end+sc3-_pa, aln_end+sc3), i.e. the OUTER end (right).
                aligner_w = (x1 - x0) - _pa if _pa > 0 else (x1 - x0)
                aligner_w = max(0, aligner_w)
                # Split the non-pA clip into slate-blue (walkback-removed,
                # INNERMOST = nearest the body) + grey (native soft-clip) so the
                # overview matches the detail's red/grey/green decomposition.
                _wb = min(aligned.get("ov_wb_n", 0), aligner_w)
                if _wb > 0:
                    ax.add_patch(Rectangle((x0, y_center - 0.18), _wb, 0.36,
                                           facecolor=SOFTCLIP_WB_BG,
                                           edgecolor=SOFTCLIP_WB_FG,
                                           linewidth=0.5))
                _grey_w = aligner_w - _wb
                if _grey_w > 0:
                    ax.add_patch(Rectangle((x0 + _wb, y_center - 0.18), _grey_w, 0.36,
                                           facecolor=SOFTCLIP_BG,
                                           edgecolor=SOFTCLIP_FG,
                                           linewidth=0.5))
                if _pa > 0:
                    pa_x0 = x0 + aligner_w
                    pa_w = (x1 - x0) - aligner_w
                    if pa_w > 0:
                        ax.add_patch(Rectangle((pa_x0, y_center - 0.18), pa_w, 0.36,
                                               facecolor=SOFTCLIP_PA_BG,
                                               edgecolor=SOFTCLIP_PA_FG,
                                               linewidth=0.5))

        _place_labels(row_labels, y_center)

        # Per-row RECTIFY-corrected 3' marker — a SMALL green ▲ just below
        # each row's bar, mirroring the detail view's per-aligner ▲ (only
        # smaller for the compressed overview). Replaces the old single
        # figure-level vertical green line, which spanned all rows and could
        # not show per-aligner divergence in the corrected 3' end.
        _ov_c3 = aligned.get("ov_corr3p")
        if _ov_c3 is not None and win_start <= _ov_c3 < win_end:
            # Boundary between the called 3'-end base and the downstream nt
            # (see render_alignment_row): +1.0 plus / +0.0 minus (inversion
            # flips minus to the visual right).
            _off = 1.0 if not is_reverse else 0.0
            _cx = (_ov_c3 - win_start) + _off
            ax.plot(_cx, y_center - 0.30, marker="^", markersize=4.5,
                    color="#2e7d32", clip_on=False, zorder=7)


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
    show_junction_zoom: bool = True,
    show_agreement_track: bool = False,
    summary_tsv: Path | None = None,
    aligner_bam_dir: Path | None = None,
    gene_id: str = "",
    dpi: int = 150,
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

    # ------------------------------------------------------------------
    # Per-aligner track build: replace the legacy (minimap2 / corrected /
    # pA-rest) trio with one row per aligner, using the pA-rest BAM each
    # aligner produced (corrected + softclipped + parquet pA-tail
    # restored as soft-clip). The merged-winner bundle BAMs `corrected_bam` and
    # `parestore_bam` are kept for the orig/corr ↓ marker computation
    # (downstream of this block) but not rendered as their own rows.
    # When per_aligner_corrected_dir is unavailable we fall back to the
    # legacy trio.
    # ------------------------------------------------------------------
    per_aligner_corrected_dir = None
    if aligner_bam_dir is not None:
        candidate = Path(aligner_bam_dir).parent / "rectified" / "per_aligner"
        if candidate.exists():
            per_aligner_corrected_dir = candidate

    track_sources: list[tuple[str, Path]] = []
    # Pre-RECTIFY baseline row: minimap2's faithful alignment of the
    # UNTRIMMED read, sourced from the Dorado BAM. Dorado uses minimap2
    # under the hood for the alignment, so this BAM IS minimap2's natural
    # behavior on the full read (including poly-A tail). The body
    # typically extends 1-3 bases into the genomic A-tract at the
    # boundary before soft-clipping the rest of the tail — i.e., the
    # samtools-style 3' end here can sit slightly past the corrected
    # one, exactly the contrast that motivates RECTIFY.
    samtools_3p: Optional[int] = None
    winner_name = winner_label or "minimap2"
    # aligner -> RECTIFY corrected 3' end (walkback result), populated from the
    # per-aligner summary below; drives the per-row corrected-end markers AND
    # the per-aligner 3'-agreement track.
    _corr3p_by_aligner: dict[str, int] = {}
    # aligner -> (raw reference_start, raw reference_end) from the RAW aligner
    # BAM (pre-RECTIFY). Used to mark the force-aligned-then-walkback-removed
    # span (raw extent beyond the corrected body) distinctly from native clips.
    _raw_span_by_aligner: dict[str, tuple] = {}
    if aligner_bam_dir is not None:
        for _al in ('minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA'):
            _rawp = Path(aligner_bam_dir) / f"validation_reads.{_al}.bam"
            if not _rawp.exists():
                continue
            _rr = fetch_read(_rawp, qname)
            if _rr is not None and not _rr.is_unmapped:
                _raw_span_by_aligner[_al] = (_rr.reference_start, _rr.reference_end)
    if aligner_bam_dir is not None:
        unrect_bam: Optional[Path] = None
        bundle = Path(aligner_bam_dir).parent
        # Preferred source: Dorado BAM (true untrimmed minimap2 alignment;
        # shows the body extending 1-3 bp into the genomic A-tract before
        # soft-clipping the rest of the pA tail — the visual motivator for
        # RECTIFY). Fallback: the pA-restored BAM, which aligns the trimmed
        # body and re-attaches the pA as soft-clip — body doesn't extend
        # into the A-tract but we still see minimap2's per-base shape.
        dorado = bundle / "validation_reads_dorado_source.bam"
        restored = bundle / "rectified" / "per_aligner" / "minimap2_unrectified.bam"
        # File-existence is not enough: the dorado BAM occasionally lags
        # the bundle and is missing reads that the per-aligner BAMs do
        # cover. Probe for the read and fall through to the restored BAM
        # on a per-read basis so we never drop the unrectified row when
        # any source has the read.
        # The dorado-source archive stores SOME reads (the build-X re-sourced
        # cat1/cat2/cat9 set) with a placeholder single all-M CIGAR rather than
        # a real per-base alignment. Painting that per-base frameshifts on the
        # read's true indels into a spurious ~all-mismatch (pink) row — the read
        # is actually well-mapped. So prefer the dorado source ONLY when it
        # carries a real (gapped/eqx) CIGAR; otherwise fall through to the
        # minimap2_unrectified BAM, which has a genuine =/X alignment.
        def _is_real_alignment(rd) -> bool:
            ct = rd.cigartuples or []
            return bool(ct) and not (len(ct) == 1 and ct[0][0] == 0)

        _dorado_read = fetch_read(dorado, qname) if dorado.exists() else None
        _restored_ok = restored.exists() and fetch_read(restored, qname) is not None
        if _dorado_read is not None and _is_real_alignment(_dorado_read):
            unrect_bam = dorado
        elif _restored_ok:
            unrect_bam = restored
        elif _dorado_read is not None:
            unrect_bam = dorado  # last resort — no real-CIGAR source available
        if unrect_bam is not None:
            track_sources.append(("minimap2 (unrectified)", unrect_bam))
            unrect_read = fetch_read(unrect_bam, qname)
            if unrect_read is not None:
                if unrect_read.is_reverse:
                    samtools_3p = unrect_read.reference_start
                else:
                    samtools_3p = unrect_read.reference_end - 1

    if per_aligner_corrected_dir is not None:
        # Identify the actual HP-ED winner for this read from summary_tsv
        # when available. The caller's `winner_label` is per-category default
        # (CAT_TO_ALIGNER_DRS) and may not match the read-level winner.
        winner_name = winner_label or "minimap2"
        if summary_tsv is not None and Path(summary_tsv).exists():
            try:
                import csv
                with open(summary_tsv) as _f:
                    for _row in csv.DictReader(_f, delimiter='\t'):
                        if _row.get('read_id') != qname:
                            continue
                        _al = _row.get('_aligner')
                        if _row.get('_is_winner') in ('1', 'True', 'true'):
                            winner_name = _al or winner_name
                        # Per-aligner RECTIFY corrected 3' end (walkback result),
                        # used to draw a corrected-end marker in every per-aligner
                        # row (not just the winner). Collect all aligners.
                        _c3 = (_row.get('corrected_3prime') or '').strip()
                        if _al and _c3:
                            try:
                                _corr3p_by_aligner[_al] = int(_c3)
                            except ValueError:
                                pass
            except Exception:
                pass
        # Group aligners by identical post-M-expansion CIGAR. Aligners in the
        # same group carry the same alignment shape (interchangeable for the
        # overview/zoom view); their HP-ED differs only in CIGAR mismatch
        # placement. Render one row per group with comma-separated labels.
        try:
            with pysam.FastaFile(str(genome_fa)) as _fa_grp:
                _grp_chrom_seq = _fa_grp.fetch(chrom).upper()
        except Exception:
            _grp_chrom_seq = None
        all_aligners = ['minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA']
        grp_order: list[tuple] = []
        grp_aligners: dict[tuple, list[str]] = {}
        grp_bam: dict[tuple, Path] = {}
        for aligner in all_aligners:
            bp = per_aligner_corrected_dir / f"{aligner}.pA_rest.bam"
            if not bp.exists():
                continue
            read = fetch_read(bp, qname)
            if read is None or read.cigartuples is None:
                continue
            norm = _split_m_into_eqx(
                list(read.cigartuples), read.reference_start,
                _grp_chrom_seq, read.query_sequence or "",
            )
            key = (read.reference_start, tuple(norm))
            if key not in grp_aligners:
                grp_aligners[key] = []
                grp_order.append(key)
                grp_bam[key] = bp
            grp_aligners[key].append(aligner)
            # Prefer the winner's BAM as the group's representative so the
            # winner's bases are shown in any per-base panels.
            if aligner == winner_name:
                grp_bam[key] = bp
        # Order groups: the one containing the winner first, then the rest
        # in original (first-occurrence) order.
        winner_key = next(
            (k for k in grp_order if winner_name in grp_aligners[k]),
            None,
        )
        if winner_key is not None:
            grp_order = [winner_key] + [k for k in grp_order if k is not winner_key]

        for key in grp_order:
            names = grp_aligners[key]
            # Put winner first within its own label if present.
            if winner_name in names:
                names = [winner_name] + sorted(a for a in names if a != winner_name)
                tag = "winner: " + ", ".join(names)
            else:
                tag = ", ".join(sorted(names))
            track_sources.append((tag, grp_bam[key]))

    if not track_sources:
        # Fallback: legacy 3-row layout.
        track_sources = list(upstream_tracks) + [
            ("corrected", corrected_bam),
            ("pA-rest",   parestore_bam),
        ]

    # ED values per aligner (HP-ED + raw ED, for label annotation +
    # cross-chrom comparison text).
    aligner_ed_map: dict[str, tuple[Optional[float], Optional[float]]] = {}
    if summary_tsv is not None and aligner_bam_dir is not None:
        per_aligner_corrected_dir = (
            Path(aligner_bam_dir).parent / "rectified" / "per_aligner"
        )
        if not per_aligner_corrected_dir.exists():
            per_aligner_corrected_dir = None
        try:
            with pysam.FastaFile(str(genome_fa)) as fa_ed:
                full_chrom_seq_for_ed = fa_ed.fetch(chrom).upper()
        except Exception:
            full_chrom_seq_for_ed = None
        ed_panel_data = load_correction_panel_data(
            qname=qname,
            summary_tsv=Path(summary_tsv),
            aligner_bam_dir=Path(aligner_bam_dir),
            per_aligner_corrected_dir=per_aligner_corrected_dir,
            chrom_seq=full_chrom_seq_for_ed,
            genome_fa=Path(genome_fa),
        )
        if ed_panel_data is not None:
            for r in ed_panel_data.get('rows', []):
                if r.get('missing'):
                    continue
                aligner_ed_map[r['aligner']] = (r.get('hp_ed'), r.get('raw_ed'))

    def _fmt_ed_for_label(label: str) -> str:
        """Append HP-ED / raw-ED to a row label like 'gapmm2' or
        'winner: deSALT, gapmm2, …'. For grouped labels, list ED per
        aligner so the reader sees within-group variance."""
        # Extract the aligner names from the label.
        s = label.split(":", 1)[-1].strip() if ":" in label else label
        names = [n.strip() for n in s.split(",")]
        eds = []
        for n in names:
            if n in aligner_ed_map:
                hp, raw = aligner_ed_map[n]
                hp_s = f"{hp:.1f}" if hp is not None else "—"
                raw_s = f"{raw:.1f}" if raw is not None else "—"
                eds.append(f"{n} HP {hp_s} / raw {raw_s}")
        if not eds:
            return label
        return f"{label}\n  " + " · ".join(eds)

    aligned_set: list[tuple[str, dict]] = []
    mapped_read_obj = None
    # Read this read's pA tail length from the trim-metadata TSV so we
    # can color the pA portion of the trailing soft-clip distinctly from
    # the aligner-emitted soft-clip of the trimmed body.
    _polya_len = _load_polya_lengths(POLYA_TRIM_TSV).get(qname, 0)
    # Cross-chrom aligners (e.g. deSALT picking a paralogous locus on a
    # different chromosome) get filtered out of the per-base view —
    # rendering them against the rendered chrom's ref bases would
    # produce a wall of fake mismatches. We collect them here and
    # surface as a "(diff chrom: chrX:N-M)" annotation below the title.
    # (label, chrom, start, end, hp_ed, raw_ed)
    cross_chrom_rows: list[tuple[str, str, int, int, Optional[float], Optional[float]]] = []
    pruned_track_sources: list[tuple[str, Path]] = []
    # Each cross-chrom entry also retains the BAM path so the mini-panel
    # can later re-walk the CIGAR and paint it against the aligner's own
    # chromosome.
    cross_chrom_bams: list[tuple[str, Path]] = []
    for label, bp in track_sources:
        read = fetch_read(bp, qname)
        # Drop the track if the BAM has no record for this read (rare —
        # e.g., the Dorado source BAM occasionally doesn't include a
        # read that the aligner BAMs do). Rendering an empty row is
        # worse than omitting it.
        if read is None:
            continue
        if read is not None and read.reference_name != chrom:
            # Pull ED from the first matching aligner name in the label
            # (cross-chrom rows are typically single-aligner anyway).
            name_for_ed = (label.split(":", 1)[-1].strip() if ":" in label
                           else label).split(",", 1)[0].strip()
            hp, raw = aligner_ed_map.get(name_for_ed, (None, None))
            cross_chrom_rows.append((
                label, read.reference_name,
                read.reference_start, read.reference_end, hp, raw,
            ))
            cross_chrom_bams.append((label, bp))
            continue
        # ``mapped_read_obj`` = the winner's record; used for orig_3p
        # override below. The winner's row is tagged "winner: …" in the
        # per-aligner layout, or `bp == mapped_bam` in the legacy fallback.
        if mapped_read_obj is None and (
            label.startswith("winner") or bp == mapped_bam
        ):
            mapped_read_obj = read
        a = walk_cigar(read, win_start, win_end)
        compute_mismatches(a, ref_seq)
        # The "minimap2 (unrectified)" row aligns the read with its pA
        # tail still attached — its soft-clip is whatever minimap2 chose
        # to clip on the FULL untrimmed read, not the result of a
        # rectify trim + re-attach. Don't impose the green/grey split.
        is_unrect = isinstance(label, str) and "unrect" in label.lower()
        row_polya = 0 if is_unrect else _polya_len
        _mark_polya_softclip(a, row_polya, bool(read.is_reverse))
        a["polya_len"] = row_polya
        a["is_reverse"] = bool(read.is_reverse)
        aligned_set.append((label, a))
        pruned_track_sources.append((label, bp))
    track_sources = pruned_track_sources

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

    # Compute orig_5p (from the winner's UNCORRECTED aligner BAM) and
    # corr_5p (from the rectified BAM, i.e. post-5'-rescue if it fired).
    # These show only if they fall inside the rendered window — usually
    # they won't for the default 3'-end-centered window, but they're
    # invaluable in overview-window views and for 5'-rescue Cat3 cases.
    orig_5p: Optional[int] = None
    corr_5p: Optional[int] = None
    if mapped_read_obj is not None:
        # The "winner" track in track_sources points at that aligner's
        # rectified pA_rest BAM (post-5'-rescue). To get the ORIGINAL 5'
        # position pre-rescue, peek at the same aligner's UNCORRECTED BAM.
        if mapped_read_obj.is_reverse:
            corr_5p = mapped_read_obj.reference_end - 1
        else:
            corr_5p = mapped_read_obj.reference_start
    if aligner_bam_dir is not None and winner_label:
        raw_winner_bp = Path(aligner_bam_dir) / f"validation_reads.{winner_label}.bam"
        if raw_winner_bp.exists():
            raw_read = fetch_read(raw_winner_bp, qname)
            if raw_read is not None:
                if raw_read.is_reverse:
                    orig_5p = raw_read.reference_end - 1
                else:
                    orig_5p = raw_read.reference_start

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
    # Load the full chrom once so the overview can expand M-ops into =/X
    # for per-base match visibility (independent of the aligner's CIGAR style).
    _full_chrom_seq: Optional[str] = None
    if show_overview:
        try:
            with pysam.FastaFile(str(genome_fa)) as _fa:
                _full_chrom_seq = _fa.fetch(chrom).upper()
        except Exception:
            _full_chrom_seq = None
        ov_start, ov_end = overview_window(aligned_set)
        # Re-walk CIGARs in the overview window. Use the SAME track_sources
        # list so the overview matches the per-base panel ordering exactly.
        for label, bp in track_sources:
            read = fetch_read(bp, qname)
            a = walk_cigar(read, ov_start, ov_end)
            # Expand M → =/X so the painted bars show match/mismatch even
            # for aligners that emit M-style CIGAR (minimap2 family).
            if a.get("cigartuples") and a.get("ref_start_raw") is not None:
                a["cigartuples"] = _split_m_into_eqx(
                    a["cigartuples"], a["ref_start_raw"],
                    _full_chrom_seq, a.get("query_seq_raw", ""),
                )
            is_unrect_ov = isinstance(label, str) and "unrect" in label.lower()
            ov_polya = 0 if is_unrect_ov else _polya_len
            is_rev_ov = bool(read.is_reverse) if read is not None else False
            # Walk back through the BAM-orientation soft-clip bases to
            # absorb any trailing pA bases the trim missed (and to compute
            # `polya_len_effective` that the soft-clip patches read).
            _mark_polya_softclip(a, ov_polya, is_rev_ov)
            a["polya_len"] = ov_polya
            a["is_reverse"] = is_rev_ov
            # Per-row corrected-3' (for the green ▲) + walkback extent (for the
            # slate-blue tail segment) so the overview strips carry the SAME
            # red/grey/green styling + corrected-end marker as the detail rows.
            a["ov_corr3p"] = None
            a["ov_wb_n"] = 0
            if is_unrect_ov:
                a["ov_corr3p"] = samtools_3p   # baseline → naive 3'
            else:
                _al0 = None
                for _al in label.replace("winner:", " ").replace(",", " ").split():
                    if _al in _corr3p_by_aligner and a["ov_corr3p"] is None:
                        a["ov_corr3p"] = _corr3p_by_aligner[_al]
                    if _al in _raw_span_by_aligner and _al0 is None:
                        _al0 = _al
                if _al0 is not None and a.get("aln_start") is not None:
                    _rs, _re = _raw_span_by_aligner[_al0]
                    a["ov_wb_n"] = max(
                        0, (a["aln_start"] - _rs) if is_rev_ov else (_re - a["aln_end"]))
            overview_aligned_set.append((label, a))

    # ---- Collect junctions for the zoom panels ----
    # Junction = N-op in the corrected (rectified) BAM's CIGAR. This is the
    # canonical post-rectification view of the read's splice structure.
    # For each junction we'll render one stacked panel showing 50 bp of exon
    # on each side + 25 bp of intron stub on each side + a crushed gap with
    # the intron length label.
    zoom_junctions: list[tuple[int, int]] = []
    zoom_aligned_sets: list[list[tuple[str, dict]]] = []
    zoom_ref_seqs: list[str] = []   # ref string covering [donor-50, acceptor+50)
    zoom_window_starts: list[int] = []
    if show_junction_zoom:
        zoom_junctions = _collect_corrected_junctions(corrected_bam, qname)
        for donor, acceptor in zoom_junctions:
            zwin_start = donor - ZOOM_EXON_FLANK
            zwin_end   = acceptor + ZOOM_EXON_FLANK
            zref = load_window_ref(genome_fa, chrom, zwin_start, zwin_end)
            zaligned: list[tuple[str, dict]] = []
            for label, bp in track_sources:
                read = fetch_read(bp, qname)
                a = walk_cigar(read, zwin_start, zwin_end)
                compute_mismatches(a, zref)
                zaligned.append((label, a))
            zoom_aligned_sets.append(zaligned)
            zoom_ref_seqs.append(zref)
            zoom_window_starts.append(zwin_start)

    # End-correction summary panel removed per the per-aligner-row redesign
    # (Phase 1, 2026-05-18). Per-aligner data is now visible directly in the
    # alignment rows; the panel's HP-ED / raw-ED / eff / pick metadata will
    # be re-introduced inline next to row labels in Phase 2.
    panel_data = None

    n = win_end - win_start
    fig_w = max(7.5, n * 0.13)

    # Each entry in track_sources becomes one alignment panel. Common case:
    # 3 panels (minimap2-only + corrected + pA-rest) or 4 panels (minimap2 +
    # winner + corrected + pA-rest) when winner != minimap2.
    track_labels = [lbl for lbl, _ in track_sources]
    # The per-aligner 3'-agreement track (panel "bedgraph") is opt-in
    # (show_agreement_track) — on the validation bundle it's redundant with the
    # per-aligner-row corrected-3' ▲s; keep it available for other plotting
    # scenarios (e.g. real multi-read pileups).
    _agree = ["bedgraph"] if show_agreement_track else []
    _agree_h = [1.7] if show_agreement_track else []
    panels = _agree + ["ref"] + track_labels + ["ticks"]
    # Read rows much taller (2.0) for the insertion-pill stagger ABOVE the row.
    # Tick row bumped from 0.35 → 0.55 so coord numbers don't crowd the
    # alignment-track baseline above them. ref row trimmed 1.65 → 0.9 now that
    # its ylim is tight (markers removed) — keeps the ref bases ≈ the per-aligner
    # base size. The agreement track (when shown) needs ~1.7 for the 5-cell stack.
    height_ratios = _agree_h + [0.9] + [2.0] * len(track_sources) + [0.55]

    # Cross-chrom mini-panels prepended at the TOP (right above the main
    # bedgraph). Each contributes a painted-CIGAR row + tick row with its
    # own ref-coord x-axis so the reader can judge the aligner's actual
    # placement on its own chromosome. Ordering: ED table → cross-chrom
    # → bedgraph / ref / per-aligner / ticks.
    cc_panels: list[str] = []
    cc_heights: list[float] = []
    for cc_idx in range(len(cross_chrom_rows)):
        cc_panels.extend([f"cc_aln_{cc_idx}", f"cc_ticks_{cc_idx}"])
        cc_heights.extend([1.4, 0.55])
    panels = cc_panels + panels
    height_ratios = cc_heights + height_ratios

    # (ED values are now rendered inline in the overview's left margin —
    # see render_overview's aligner_ed_map handling — so no separate
    # ED-table panel is added here anymore.)
    if show_overview:
        panels = ["overview", "ov_ticks"] + panels
        # Overview tick row bumped from 0.30 → 0.50 (same reason as above).
        # Overview body bumped to fit per-segment CIGAR labels above
        # (intron-len) and below (segment-length-and-op like "55=" / "3D").
        height_ratios = [1.2 + 0.32 * len(aligned_set), 0.50] + height_ratios
    if panel_data is not None:
        # Panel: top of figure. Height scales with row count.
        n_aligner_rows = len(panel_data['rows'])
        n_text_rows = n_aligner_rows + 4  # meta + blank + h1 + h2
        panel_height = 0.4 + 0.27 * n_text_rows
        panels = ["correction_panel"] + panels
        height_ratios = [panel_height] + height_ratios
    # Insert one zoom block per junction AFTER the overview (and ov_ticks) but
    # BEFORE the bedgraph. Each block: 1 ref-zoom row + N track-zoom rows + 1
    # tick row. They share the existing figure-wide x-axis width (158 cols),
    # narrower than the typical 3'-end window so they appear inset visually
    # without their own axes — the unequal widths are absorbed by matplotlib.
    insert_at = (2 if show_overview else 0)  # right after overview block
    for j_idx, (donor, acceptor) in enumerate(zoom_junctions):
        block = [f"zref_{j_idx}"] + [f"zaln_{j_idx}_{k}" for k in range(len(track_sources))] + [f"zticks_{j_idx}"]
        # zoom alignment rows now use ylim (0, 2.60) to give insertion
        # pills headroom above the row (matching render_alignment_row).
        # Bump from 2.0 → 3.2 so the base bars stay visually comparable.
        block_h = [1.2] + [3.2] * len(track_sources) + [0.35]
        panels[insert_at:insert_at] = block
        height_ratios[insert_at:insert_at] = block_h
        insert_at += len(block)
    fig_h = sum(height_ratios) * 0.45 + 0.5
    # Width: the zoom panels are 158 visual cols; let them set a minimum so
    # the figure isn't stretched too thin when the 3'-end window is small.
    zoom_min_width = (ZOOM_TOTAL_COLS * 0.10) if zoom_junctions else 0.0
    fig_w = max(fig_w, zoom_min_width)
    fig, axes = plt.subplots(
        len(panels), 1, figsize=(fig_w, fig_h),
        gridspec_kw={"height_ratios": height_ratios},
    )
    if title:
        fig.suptitle(title, fontsize=SIZE_TITLE, y=0.995)
    # Painted-CIGAR color legend (top center, just below the title).
    # Single horizontal strip — each op letter + its color swatch in a
    # compact key. Renders once per figure (not per panel).
    # M is omitted: rectify-corrected BAMs are =-SEQ decoded so the
    # painted-CIGAR overview expands M into =/X via _split_m_into_eqx.
    # If a stale aligner BAM slips through, any residual M segments are
    # rendered with the = color (CIGAR_OV_COLOR[0] below) so they merge
    # visually with match runs.
    # SINGLE merged key (was two stacked rows — too crowded). Drops the
    # duplicates: the old "S soft-clip" op is folded into "native soft-clip",
    # and the detail "mismatch" swatch into "X mismatch". The first cluster is
    # the overview's painted-CIGAR colors; the rest are the per-base tail-cell
    # shadings (exact cell facecolor+edgecolor) + the corrected-3' marker.
    # M is omitted: =-SEQ decoded BAMs expand M into =/X via _split_m_into_eqx;
    # any residual M renders with the = color so it merges with match runs.
    from matplotlib.patches import Patch as _Patch
    from matplotlib.lines import Line2D as _Line2D
    handles = [
        _Patch(facecolor="#cfd8dc", edgecolor="none", label="= match"),
        _Patch(facecolor="#e91e63", edgecolor="none", label="X mismatch"),
        _Patch(facecolor="#ff9800", edgecolor="none", label="D deletion"),
        _Patch(facecolor="#9c27b0", edgecolor="none", label="I insertion"),
        _Patch(facecolor="#aaaaaa", edgecolor="none", label="N intron"),
        _Patch(facecolor=SOFTCLIP_BG, edgecolor=SOFTCLIP_FG, linewidth=0.6,
               label="native soft-clip"),
        _Patch(facecolor=SOFTCLIP_WB_BG, edgecolor=SOFTCLIP_WB_FG, linewidth=0.6,
               label="walkback-removed"),
        _Patch(facecolor=SOFTCLIP_PA_BG, edgecolor=SOFTCLIP_PA_FG, linewidth=0.6,
               label="poly(A) restored"),
        _Line2D([0], [0], marker="^", color="none", markerfacecolor="#2e7d32",
                markeredgecolor="#2e7d32", markersize=7, label="corrected 3'"),
    ]
    leg = fig.legend(
        handles=handles, loc="upper center", ncol=len(handles),
        bbox_to_anchor=(0.5, 0.972),
        frameon=False, fontsize=SIZE_SMALL + 0.5,
        handlelength=1.0, handleheight=0.9, handletextpad=0.3,
        columnspacing=1.2, borderpad=0.1,
    )
    # Note: when the per-aligner layout is active, every alignment row shows
    # that aligner's corrected output with the parquet-trimmed poly(A) tail
    # restored as 3' soft-clip. So mismatches/indels visible in the body
    # are post-correction; the trailing soft-clip is two-tone — grey for
    # the aligner's soft-clip of the trimmed body, green for the pA tail.
    # Bottom italic banners — pushed below the per-base alignment block
    # with breathing room between them. Font bumped 7 → 8.5 for
    # readability. subplots_adjust(bottom=0.07) gives them room.
    if track_sources and any(lbl.startswith("winner:") or lbl == "minimap2"
                             for lbl, _ in track_sources):
        fig.text(0.5, 0.038,
                 "per-aligner rows: corrected CIGAR + poly(A) tail restored as "
                 "soft-clip from the trim parquet (tail-cell colors keyed above)",
                 ha="center", fontsize=SIZE_TEXT, color="#666", style="italic")
    if cross_chrom_rows:
        chunks = []
        for (lbl, cn, st, en, hp, raw) in cross_chrom_rows:
            hp_s = f"{hp:.1f}" if hp is not None else "—"
            raw_s = f"{raw:.1f}" if raw is not None else "—"
            chunks.append(f"{lbl} on {cn}:{st:,}–{en:,}  "
                          f"(EER ED {hp_s} / raw {raw_s})")
        fig.text(0.5, 0.010,
                 "cross-chrom (excluded from per-base view): " + "; ".join(chunks),
                 ha="center", fontsize=SIZE_TEXT, color="#b71c1c", style="italic",
                 fontweight="bold")
    # Provenance footer (project figure rule): small grey, self-documenting —
    # data source + key processing + script filename + render date.
    fig.text(
        0.5, 0.002,
        "source: wt_by4742_rep1 ONT-DRS · S. cerevisiae S288C (R64-5-1) · "
        "RECTIFY HP-ED winner-selection + poly(A) restore · "
        "render_read_alignment.py · " + datetime.date.today().isoformat(),
        ha="center", va="bottom", fontsize=5.5, color="0.45",
    )

    ax_iter = iter(axes)
    if panel_data is not None:
        ax_panel = next(ax_iter)
        render_correction_panel(ax_panel, panel_data, qname, gene_id=gene_id)
    if show_overview:
        ax_ov = next(ax_iter)
        _is_rev_main = bool(mapped_read_obj is not None and mapped_read_obj.is_reverse)
        render_overview(ax_ov, overview_aligned_set, ov_start, ov_end,
                        orig_3p, corr_3p, aligner_ed_map=aligner_ed_map,
                        is_reverse=_is_rev_main)
        ax_ov_ticks = next(ax_iter)
        # When the overview has ED columns on the right, the tick row
        # needs to match the extended xlim and drop the rightmost label
        # so it doesn't bleed under the ED area.
        _ov_ed_pad = 0.14 if aligner_ed_map else 0.0
        render_axis_ticks(ax_ov_ticks, ov_start, ov_end,
                          right_pad_frac=_ov_ed_pad,
                          suppress_last_tick=bool(aligner_ed_map),
                          is_reverse=_is_rev_main)

    # ---- ED-comparison table (DISABLED — moved inline into overview) ----
    if False and aligner_ed_map:
        ax_ed = next(ax_iter)
        ax_ed.set_xlim(0, 1)
        ax_ed.set_ylim(0, 1)
        ax_ed.set_axis_off()
        # Map aligner name → cross-chrom tag (when applicable).
        cc_chrom_by_name: dict[str, str] = {}
        for (cc_label, cc_chrom_, *_rest) in cross_chrom_rows:
            stripped = (cc_label.split(":", 1)[-1] if ":" in cc_label
                        else cc_label)
            for name_ in stripped.split(","):
                cc_chrom_by_name[name_.strip()] = cc_chrom_

        # Build rows in OVERVIEW order: per-base aligned_set rows first,
        # then cross-chrom rows. Each row carries the group label
        # (comma-separated aligners), with EER-ED + raw-ED values
        # comma-separated in matching order.
        def _per_aligner_vals(label: str) -> tuple[list[str], list[str]]:
            stripped = (label.split(":", 1)[-1] if ":" in label
                        else label).strip()
            names = [n.strip() for n in stripped.split(",")]
            eer_vals, raw_vals = [], []
            for n in names:
                hp, raw = aligner_ed_map.get(n, (None, None))
                eer_vals.append(f"{hp:.1f}" if hp is not None else "—")
                raw_vals.append(f"{raw:.0f}" if raw is not None else "—")
            return eer_vals, raw_vals

        row_specs: list[tuple[str, list[str], list[str]]] = []
        for (lbl, _a) in aligned_set:
            # Skip the "minimap2 (unrectified)" baseline row — its ED
            # isn't in aligner_ed_map and isn't comparable anyway.
            if "unrect" in lbl:
                continue
            eer_v, raw_v = _per_aligner_vals(lbl)
            row_specs.append((lbl, eer_v, raw_v))
        for (cc_label, cc_chrom_, cc_st, cc_en, cc_hp, cc_raw) in cross_chrom_rows:
            hp_s = f"{cc_hp:.1f}" if cc_hp is not None else "—"
            raw_s = f"{cc_raw:.0f}" if cc_raw is not None else "—"
            tagged_label = f"{cc_label}  [{cc_chrom_}]"
            row_specs.append((tagged_label, [hp_s], [raw_s]))

        # Column widths tuned for typical yeast 5-aligner case.
        W_LABEL  = 36
        W_EER    = 24
        W_RAW    = 22
        header = (f"{'':<{W_LABEL}}  "
                  f"{'EER ED':>{W_EER}}  {'raw ED':>{W_RAW}}")
        sep = ("─" * W_LABEL + "  " + "─" * W_EER + "  " + "─" * W_RAW)
        lines = [header, sep]
        for (lbl, eer_v, raw_v) in row_specs:
            # Truncate label if too long (keeps the table compact).
            lbl_disp = lbl if len(lbl) <= W_LABEL else (lbl[:W_LABEL - 1] + "…")
            lines.append(
                f"{lbl_disp:<{W_LABEL}}  "
                f"{', '.join(eer_v):>{W_EER}}  "
                f"{', '.join(raw_v):>{W_RAW}}"
            )
        ax_ed.text(0.005, 0.95, "\n".join(lines),
                   transform=ax_ed.transAxes,
                   family='monospace', fontsize=7.5,
                   va='top', ha='left', color='#222')
        from matplotlib.patches import FancyBboxPatch as _FBB
        ax_ed.add_patch(_FBB(
            (0.0, 0.0), 1.0, 1.0,
            boxstyle="round,pad=0.005,rounding_size=0.01",
            transform=ax_ed.transAxes,
            linewidth=0.5, edgecolor='#bbbbbb', facecolor='#FAFAFA',
            zorder=-1,
        ))

    # ---- Junction-zoom panels (one block per junction in corrected BAM) ----
    is_reverse_pre = bool(mapped_read_obj is not None and mapped_read_obj.is_reverse)
    for j_idx, (donor, acceptor) in enumerate(zoom_junctions):
        zref = zoom_ref_seqs[j_idx]
        zwin_start = zoom_window_starts[j_idx]
        def _ref_provider(rp: int, _zref=zref, _zs=zwin_start) -> Optional[str]:
            idx = rp - _zs
            if 0 <= idx < len(_zref):
                return _zref[idx]
            return None
        ax_zref = next(ax_iter)
        render_zoom_ref_row(ax_zref, _ref_provider, donor, acceptor,
                            is_reverse=is_reverse_pre,
                            label=f"ref (junc {j_idx + 1})")
        for label, a in zoom_aligned_sets[j_idx]:
            ax_z = next(ax_iter)
            render_zoom_alignment_row(ax_z, a, donor, acceptor, label,
                                      window_start=zwin_start,
                                      is_reverse=is_reverse_pre)
        ax_zticks = next(ax_iter)
        render_zoom_ticks(ax_zticks, donor, acceptor)
    # For minus-strand reads, render in RNA 5'→3' orientation: complement
    # every displayed base AND invert the x-axis so the higher genomic coord
    # (which is the RNA 5' end) sits on the LEFT and the lower genomic coord
    # (RNA 3' end, where the CPA + poly-A tail are) sits on the RIGHT. This
    # makes the poly-A tail appear as poly-A (not poly-T) at the natural 3'
    # end position and lines up with how the user reads RNA sequences.
    is_reverse = bool(mapped_read_obj is not None and mapped_read_obj.is_reverse)

    # ---- Cross-chrom mini-panels (rendered FIRST so they appear at the
    # TOP of the figure, right under the overview). Each has its own
    # painted-CIGAR strip + tick row in ITS chromosome's coords.
    # Anchor cross-chrom panels by aln_start (= the body start) — the
    # body is the dominant visual element and aligning bodies across
    # panels gives the cleanest visual comparison. Leading-softclip
    # differences mean the soft-clip starts don't align across rows
    # either way; bodies are what the reader scans first.
    main_anchor_offset = 20
    if mapped_read_obj is not None:
        main_anchor_offset = max(0, mapped_read_obj.reference_start - ov_start)
    main_ov_length = ov_end - ov_start
    cc_axis_indices: list[int] = []
    for cc_idx, ((cc_label, cc_chrom, cc_start, cc_end, cc_hp, cc_raw),
                 (_, cc_bp)) in enumerate(zip(cross_chrom_rows,
                                              cross_chrom_bams)):
        # Anchor by body-start (aln_start). The CC panel's body starts
        # at the same axis x as the main winner's body. Window length
        # matches the main's so the per-bp scale is identical.
        cc_win_start = max(0, cc_start - main_anchor_offset)
        cc_win_end   = cc_win_start + main_ov_length
        cc_read = fetch_read(cc_bp, qname)
        cc_aligned: list[tuple[str, dict]] = []
        cc_full_chrom_seq: Optional[str] = None
        if cc_read is not None:
            try:
                with pysam.FastaFile(str(genome_fa)) as fa_cc:
                    cc_full_chrom_seq = fa_cc.fetch(cc_chrom).upper()
            except Exception:
                cc_full_chrom_seq = None
            a_cc = walk_cigar(cc_read, cc_win_start, cc_win_end)
            if a_cc.get("cigartuples") and a_cc.get("ref_start_raw") is not None:
                a_cc["cigartuples"] = _split_m_into_eqx(
                    a_cc["cigartuples"], a_cc["ref_start_raw"],
                    cc_full_chrom_seq, a_cc.get("query_seq_raw", ""),
                )
            a_cc["polya_len"] = _polya_len
            a_cc["is_reverse"] = bool(cc_read.is_reverse)
            _mark_polya_softclip(a_cc, _polya_len, bool(cc_read.is_reverse))
            cc_label_full = f"{cc_label}  [{cc_chrom}]"
            cc_aligned.append((cc_label_full, a_cc))
        ax_cc_aln = next(ax_iter)
        ax_cc_ticks = next(ax_iter)
        cc_axis_indices.extend([axes.tolist().index(ax_cc_aln),
                                axes.tolist().index(ax_cc_ticks)])
        if cc_aligned:
            render_overview(ax_cc_aln, cc_aligned, cc_win_start, cc_win_end,
                            None, None, aligner_ed_map=aligner_ed_map,
                            show_column_headers=False)
            # No extra header text — the top-of-figure italic already
            # spells out chrom + ED, and the row label inside
            # render_overview tags this row with "[<chrom>]".
        # Same ED-pad treatment as the overview tick row.
        render_axis_ticks(ax_cc_ticks, cc_win_start, cc_win_end,
                          right_pad_frac=(0.14 if aligner_ed_map else 0.0),
                          suppress_last_tick=bool(aligner_ed_map))

    if show_agreement_track:
        ax_bg = next(ax_iter)
        render_aligner_3p_agreement(ax_bg, _corr3p_by_aligner, winner_name,
                                    win_start, win_end)
    ax_ref = next(ax_iter)
    render_ref_row(ax_ref, ref_seq, win_start, win_end, orig_3p, corr_3p,
                   is_reverse=is_reverse, samtools_3p=samtools_3p,
                   orig_5p=orig_5p, corr_5p=corr_5p)
    for label, a in aligned_set:
        ax = next(ax_iter)
        # Resolve this row's aligner(s) → its RECTIFY corrected 3' end so the
        # row can mark where walkback landed for it (skip the unrectified row).
        _is_unrect = "unrect" in label.lower()
        _rc3 = None
        if _is_unrect:
            # The unrectified baseline's green ▲ marks the NAIVE 3' end — where
            # the pA site would be called WITHOUT RECTIFY (the samtools view of
            # the raw read, pA still soft-clipped). The contrast against the
            # rectified rows' ▲ is the whole point of the baseline row.
            _rc3 = samtools_3p
        else:
            _al0 = None
            for _al in label.replace("winner:", " ").replace(",", " ").split():
                if _al in _corr3p_by_aligner and _rc3 is None:
                    _rc3 = _corr3p_by_aligner[_al]
                if _al in _raw_span_by_aligner and _al0 is None:
                    _al0 = _al
            # Crimson walkback-removed tail cells = force-aligned bases RECTIFY
            # clipped — a corrected-row concept; the unrectified baseline (above)
            # shows no correction, so it gets none.
            if _al0 is not None and a.get("aln_start") is not None:
                _rs, _re = _raw_span_by_aligner[_al0]
                _wb_n = (a["aln_start"] - _rs) if is_reverse else (_re - a["aln_end"])
                _mark_walkback_removed(a, max(0, _wb_n), is_reverse)
        render_alignment_row(ax, a, ref_seq, win_start, win_end, label,
                             is_reverse=is_reverse, corr_3p_aligner=_rc3)
    ax_ticks = next(ax_iter)
    render_axis_ticks(ax_ticks, win_start, win_end)

    # Flip the x-axis on every panel for minus-strand reads (RNA 5'→3' = left
    # to right means higher coord on the left). Done after all rendering so
    # the internal column coordinates stay genomic-oriented and only the
    # display flips. The correction-panel axis is excluded — its text content
    # already uses genomic coords with strand-direction Δ semantics baked in,
    # and inverting would mirror the table left/right.
    if is_reverse:
        skip_first = 1 if panel_data is not None else 0
        cc_set = set(cc_axis_indices)
        for i, ax in enumerate(axes):
            if i < skip_first:
                continue
            if i in cc_set:
                continue
            ax.invert_xaxis()

    # Zoom indicator: connect the detailed (ref) window to its location in the
    # whole-read overview, so the reader instantly sees WHERE in the read the
    # per-base view sits (e.g. the 3' end for these cats). A faint box marks the
    # detailed sub-region across the overview strips AND its coordinate (tick)
    # row, then two dotted lines fan from JUST BELOW the overview coordinates
    # (ov-ticks bottom) down to JUST ABOVE the ref seq (ref panel top). Added
    # AFTER the minus-strand invert so ConnectionPatch + the boxes resolve the
    # (possibly flipped) transData at draw.
    #
    # NOTE: when the optional pileup/agreement track is re-enabled it sits
    # between the ov-ticks and the ref row — the fan would then cross it, so the
    # ref-side endpoint (axesB) should be retargeted to that track's top.
    if show_overview:
        from matplotlib.patches import ConnectionPatch
        # Grey, MONOCHROME zoom "wedge" (span bracket + two slanted dotted edges
        # + a faint region tint). Deliberately grey, NOT amber: the overview
        # already uses amber/orange for divergence shading + deletions, so a
        # yellow cue would blend in and read as data. Drawn as an OUTLINE rather
        # than a filled trapezoid because a cross-axes fill (BboxConnectorPatch)
        # twists under the minus-strand axis inversion (TransformedBbox of an
        # inverted axis has reversed x-bounds). ConnectionPatch + single-axis
        # primitives resolve the flip correctly via data coords at draw time.
        _grey = "#9e9e9e"
        _ref_hi = ax_ref.get_ylim()[1]
        _ov_xL, _ov_xR = win_start - ov_start, win_end - ov_start
        _n_det = win_end - win_start
        # Just BELOW the overview coordinate numbers (which sit at y=0..~0.4 of
        # the ov-ticks panel, va="bottom" anchored at y=0). Negative y drops the
        # bracket into the gap beneath the numbers so it no longer overlaps them.
        _yb = -0.18
        # faint region tint over the zoomed sub-region in the overview strips
        _lo, _hi = ax_ov.get_ylim()
        ax_ov.add_patch(Rectangle(
            (_ov_xL, _lo), _ov_xR - _ov_xL, _hi - _lo,
            facecolor=_grey, edgecolor="none", alpha=0.12, zorder=0,
        ))
        # horizontal span bracket joining the two edges — the zoomed locus span
        ax_ov_ticks.plot([_ov_xL, _ov_xR], [_yb, _yb], color=_grey,
                         linewidth=1.2, clip_on=False, zorder=3)
        # two dotted edges fanning from the bracket down to the ref-seq top
        for _ovx, _refx in ((_ov_xL, 0), (_ov_xR, _n_det)):
            con = ConnectionPatch(
                xyA=(_ovx, _yb), coordsA="data", axesA=ax_ov_ticks,
                xyB=(_refx, _ref_hi), coordsB="data", axesB=ax_ref,
                color=_grey, linestyle=(0, (2, 2)), linewidth=0.9,
                alpha=0.85, zorder=1,
            )
            con.set_clip_on(False)
            fig.add_artist(con)

    # Left margin sized for the aligner-name column + the beige bg
    # behind it. Bottom expanded to 0.07 to give room for the two
    # bottom italic banners (per-aligner note + cross-chrom).
    plt.subplots_adjust(left=0.16, right=0.99, top=0.95, bottom=0.07, hspace=0.10)
    fig.savefig(str(out_path), dpi=dpi, bbox_inches="tight")
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
    p.add_argument("--no-junction-zoom", action="store_true",
                   help="Suppress the per-junction zoom panel(s) (50 bp exon "
                        "flanks + 25 bp intron stubs + crushed intron). Useful "
                        "when the existing 3'-end view is all you want.")
    p.add_argument("--summary-tsv", type=Path, default=None,
                   help="Per-aligner comparison TSV from merge_corrected_tsvs "
                        "(rectify/data/validation/rectified/per_aligner_summary.tsv). "
                        "When supplied with --aligner-bam-dir, renders an "
                        "End-correction + Winner-selection panel at the top "
                        "of the figure.")
    p.add_argument("--aligner-bam-dir", type=Path, default=None,
                   help="Directory holding the per-aligner uncorrected BAMs "
                        "(validation_reads.<aligner>.bam). Required with "
                        "--summary-tsv to compute the orig 5'/3' coords and "
                        "uncorrected junctions for the panel.")
    p.add_argument("--gene-id", default="",
                   help="Optional gene name to display in the panel header.")
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
        show_junction_zoom=(not args.no_junction_zoom),
        summary_tsv=args.summary_tsv,
        aligner_bam_dir=args.aligner_bam_dir,
        gene_id=args.gene_id,
    )
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
