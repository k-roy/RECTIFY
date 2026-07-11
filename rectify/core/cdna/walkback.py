"""cDNA walkback wrappers around the protocol-agnostic walkback core.

Two anchor corrections built on top of
``rectify.core.correct.walkback.walkback_3prime_with_qpos``:

  * :func:`walk_back_anchor_and_tail` — v1.12 polyA-side cleavage-site
    canonicalization with sequence-level tail-length measurement.
  * :func:`walk_forward_tss` — v1.19 bridge-G TSS canonicalization (analog of
    the 3' polyA walk-back, applied at the SSP/UMI/GGG boundary).

Also hosts the adapter-anchor finder used by both walkback and consensus
pretrim (kept here to avoid a read_info ↔ walkback import cycle).
"""
from __future__ import annotations

from typing import Optional, Tuple

import pysam

from rectify.core.correct.walkback import (
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime_with_qpos,
)

from ._constants import (
    ANCHOR_FWD,
    ANCHOR_LEN,
    ANCHOR_MAX_EDIT,
    ANCHOR_RC,
    ANCHOR_SEARCH_WIN,
    BRIDGE_LEN,
)

# edlib is an optional dep (loaded lazily so the module imports without it).
try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False


# ---- v1.19 chemistry-aware 5' boundary detection ---------------------------
# The canonical PCB114 5' boundary in BAM SEQ for orient=fwd reads is the
# 12-base pattern `TT-VVVV-TTT-GGG`:
#   - TT       UMI's penultimate TT (start of last TT-VVVV-TTT block)
#   - VVVV     random V positions (V = A/C/G, not T)
#   - TTT      UMI's trailing TTT
#   - GGG      TSO bridge (rGrGrG → GGG in cDNA after RT/PCR)
# Position right AFTER the GGG = the mRNA TSS.
#
# For orient=rev the pattern's RC = `CCC-AAA-BBBB-AA` (where B = T/G/C, not A).
# Position right BEFORE the CCC = the mRNA 5' end (= aln_end - 1 in ref coords).
_BOUNDARY_PATTERN_FWD = (
    {'T'}, {'T'},                       # 0, 1: TT
    {'A','C','G'}, {'A','C','G'},       # 2, 3: VV
    {'A','C','G'}, {'A','C','G'},       # 4, 5: VV
    {'T'}, {'T'}, {'T'},                # 6, 7, 8: TTT
    {'G'}, {'G'}, {'G'},                # 9, 10, 11: GGG (bridge)
)
_BOUNDARY_PATTERN_REV = (
    {'C'}, {'C'}, {'C'},                # 0, 1, 2: CCC (bridge, RC of GGG)
    {'A'}, {'A'}, {'A'},                # 3, 4, 5: AAA (RC of TTT)
    {'C','G','T'}, {'C','G','T'},       # 6, 7: BB
    {'C','G','T'}, {'C','G','T'},       # 8, 9: BB
    {'A'}, {'A'},                       # 10, 11: AA (RC of TT)
)
_BOUNDARY_LEN = 12
_BOUNDARY_MIN_SCORE = 10  # allow 2 mismatches across the 12-base window


def _find_adapter_anchor_pos(seq: str, orient: str) -> Optional[int]:
    """Return BAM-SEQ start position of the adapter anchor (fuzzy, Lev≤ANCHOR_MAX_EDIT),
    or None. For fwd: rightmost hit in last ANCHOR_SEARCH_WIN bp. For rev: leftmost
    hit in first ANCHOR_SEARCH_WIN bp."""
    if not HAS_EDLIB:
        return None
    if orient == "fwd":
        off = max(0, len(seq) - ANCHOR_SEARCH_WIN)
        win = seq[off:]
        r = edlib.align(ANCHOR_FWD, win, mode="HW", task="locations", k=ANCHOR_MAX_EDIT)
        if r["editDistance"] == -1 or not r["locations"]: return None
        # edlib HW/locations can return a location whose START is None (end found but
        # start not localizable); treat as "adapter anchor not confidently placed" and
        # fall back to unanchored poly-A detection in the caller (pretrim_consensus).
        start = r["locations"][-1][0]
        return None if start is None else off + start
    win = seq[:ANCHOR_SEARCH_WIN]
    r = edlib.align(ANCHOR_RC, win, mode="HW", task="locations", k=ANCHOR_MAX_EDIT)
    if r["editDistance"] == -1 or not r["locations"]: return None
    start = r["locations"][0][0]
    return None if start is None else start


def walk_back_anchor_and_tail(read: pysam.AlignedSegment,
                               chrom_seq: str,
                               orient: str) -> Tuple[int, int]:
    """Return (canonical_cleavage_anchor_ref_pos, polyA_tail_length_in_basecalled_orient).

    Walk-back direction is determined by **orient** (not is_reverse), because the
    polyA tail's position in BAM SEQ depends on where the SSP/UMI was found, not
    on the alignment direction. Strand-A and Strand-B reads of the same molecule
    can both produce orient=fwd records (Strand A on +gene → is_reverse=False;
    Strand B on −gene → is_reverse=True); in both cases the BAM SEQ structure is
    SSP-UMI-cDNA-polyA-adapter (left-to-right), so polyA is at the RIGHT regardless.

      - orient='fwd': polyA at RIGHT of BAM SEQ → side=right, stop_base='A'.
      - orient='rev': polyA-as-T at LEFT of BAM SEQ (RC of basecalled polyA in
                      + ref coords) → side=left, stop_base='T'.

    The read-vs-reference scan delegates to
    :func:`rectify.core.correct.walkback.walkback_3prime_with_qpos`, which is the
    same protocol-agnostic core used by ``walkback_drs`` and
    ``walkback_quantseq_rev``. This function adds cDNA-specific tail-length
    counting using the SSP/adapter anchor on top of the shared scan.

    Tail length is the count of basecalled-A bases (= ref-T for orient=rev) in the
    BAM-SEQ region between the cleavage anchor and the adapter anchor (or 200 bp
    fallback if adapter not found).
    """
    seq = read.query_sequence
    # Fallback anchor: polyA-side end of BAM SEQ in ref coords (orient-based).
    # orient=fwd → polyA at RIGHT → aln_end-1. orient=rev → polyA at LEFT → aln_start.
    fallback_anchor = ((read.reference_end or 1) - 1) if orient == "fwd" else read.reference_start
    if seq is None:
        return (fallback_anchor, 0)

    if orient == "fwd":
        side = THREE_PRIME_SIDE_RIGHT
        tail_base = 'A'
    else:
        side = THREE_PRIME_SIDE_LEFT
        tail_base = 'T'

    _orig, corr, _applied, anchor_qp = walkback_3prime_with_qpos(
        read, chrom_seq, side, stop_base=tail_base
    )
    if anchor_qp < 0:
        # No usable scan (empty pairs, missing sequence). Preserve the
        # historical fallback behavior — anchor at the BAM SEQ poly-A-side end
        # with tail length 0.
        return (fallback_anchor, 0)

    # Compute tail length using the SSP/adapter anchor on the basecalled-tail
    # side of the cleavage anchor.
    if orient == "fwd":
        adp = _find_adapter_anchor_pos(seq, orient)
        if adp is not None and adp > anchor_qp:
            tail_seg = seq[anchor_qp + 1: adp]
        else:
            tail_seg = seq[anchor_qp + 1: anchor_qp + 1 + 200]
    else:
        adp = _find_adapter_anchor_pos(seq, orient)
        if adp is not None:
            adp_end = adp + ANCHOR_LEN
            if adp_end < anchor_qp:
                tail_seg = seq[adp_end: anchor_qp]
            else:
                tail_seg = seq[max(0, anchor_qp - 200): anchor_qp]
        else:
            tail_seg = seq[max(0, anchor_qp - 200): anchor_qp]
    tail_a = tail_seg.count(tail_base)
    return (corr, tail_a)


def _score_boundary_window(window: str, pattern) -> int:
    """Score a 12-base window against the canonical pattern. Returns 0-12."""
    if len(window) != _BOUNDARY_LEN:
        return 0
    score = 0
    for i, allowed in enumerate(pattern):
        if window[i] in allowed:
            score += 1
    return score


def _find_boundary_match(seq: str, search_start: int, search_end: int,
                          rc: bool = False) -> Optional[Tuple[int, int]]:
    """Find the best 12-base canonical-pattern match in `seq[search_start:search_end]`.

    Returns (match_qpos, score) of the highest-scoring window, or None if no
    window meets the minimum score threshold AND has the 3 bridge bases intact.
    """
    pattern = _BOUNDARY_PATTERN_REV if rc else _BOUNDARY_PATTERN_FWD
    # Indices of the bridge bases within the pattern (where the chemistry-fixed
    # G's / C's sit — these must all match for a valid boundary call)
    bridge_indices = (0, 1, 2) if rc else (9, 10, 11)
    bridge_char = 'C' if rc else 'G'
    best_score = 0
    best_qpos = -1
    for qpos in range(search_start, search_end - _BOUNDARY_LEN + 1):
        window = seq[qpos: qpos + _BOUNDARY_LEN].upper()
        # Require all 3 bridge bases to match
        if not all(window[i] == bridge_char for i in bridge_indices):
            continue
        score = _score_boundary_window(window, pattern)
        if score > best_score:
            best_score = score
            best_qpos = qpos
    if best_qpos < 0 or best_score < _BOUNDARY_MIN_SCORE:
        return None
    return (best_qpos, best_score)


def walk_forward_tss(read: pysam.AlignedSegment, orient: str, chrom_seq: str
                     ) -> int:
    """v1.19 (corrected): bridge-G walk-forward — analog of 3' polyA walk-back.

    Chemistry: PCB114 TSO has rGrGrG. After RT, the cDNA carries the RC at its 3'
    end. The BAM SEQ for orient=fwd reads has the structure
        ...SSP_FWD - UMI (ending in TTT) - GGG (bridge) - mRNA TSS - mRNA body...
    The bridge GGG (3 bp, fixed by chemistry) is followed by the mRNA TSS.

    When the genome immediately upstream of the TSS happens to contain G's, the
    aligner can extend the alignment LEFT through the bridge G's into those
    upstream genomic G's (each G/G match adds to the alignment score). The
    reported aln_start ends up 1-3 bp UPSTREAM of the true TSS.

    Fix (symmetric to polyA walk-back): walk forward through aligned G's that
    match genomic G's, capped at the bridge length (3). The corrected position
    is at the first non-G match — the "leftmost possible TSS" convention
    (analogous to "leftmost possible CPA site" for polyA cleavage).

    Sanity-check landmark: before applying the walk, verify the soft-clip
    immediately before the aligned region contains the UMI's TTT tail (for
    orient=fwd) or its RC = AAA (for orient=rev). If the chemistry pattern is
    broken (truncated UMI, Type-2-like), we skip the correction.

    Trade-off: genes whose mRNA legitimately starts with G(s) get their TSS
    reported up to 3 bp downstream of truth. This is a systematic bias that
    matches the analogous bias for polyA cleavage at genomic A-tracts.
    """
    seq = read.query_sequence
    if seq is None:
        return read.reference_start if orient == "fwd" else (read.reference_end or 1) - 1

    pairs = read.get_aligned_pairs(matches_only=True)
    if not pairs:
        return read.reference_start if orient == "fwd" else (read.reference_end or 1) - 1

    # Approach: find the most-parsimonious match of the 12-base canonical
    # boundary pattern (TT-VVVV-TTT-GGG for fwd; CCC-AAA-BBBB-AA for rev) in a
    # window around the alignment boundary. The pattern's bridge bases (GGG /
    # CCC) MUST all match — they're the chemistry-fixed signal. Once we know
    # the pattern's location in qpos, the boundary between bridge and mRNA TSS
    # is at: (for fwd) pattern_qpos + 12; (for rev) pattern_qpos.
    if orient == "fwd":
        aln_start_qpos = pairs[0][0]
        aln_start_rpos = pairs[0][1]
        # Search window: 6 bp before to 3 bp after the expected pattern start.
        # Expected pattern start = aln_start_qpos - BOUNDARY_LEN.
        expected_pat_start = aln_start_qpos - _BOUNDARY_LEN
        search_start = max(0, expected_pat_start - 6)
        search_end = min(len(seq), expected_pat_start + 3 + _BOUNDARY_LEN)
        match = _find_boundary_match(seq, search_start, search_end, rc=False)
        if match is None:
            return aln_start_rpos  # chemistry pattern not detectable — skip
        pat_qpos, _ = match
        tss_qpos = pat_qpos + _BOUNDARY_LEN  # qpos right after GGG = the TSS
        # If tss_qpos <= aln_start_qpos: alignment already starts at or past TSS
        # (over-soft-clip case — don't make it worse).
        if tss_qpos <= aln_start_qpos:
            return aln_start_rpos
        # tss_qpos > aln_start_qpos: aligner under-clipped by (tss_qpos - aln_start_qpos)
        # bridge bases. Walk forward in the alignment to that qpos.
        shift_needed = tss_qpos - aln_start_qpos
        # Cap shift at BRIDGE_LEN (3) — chemistry says max bridge is 3 bases.
        shift_needed = min(shift_needed, BRIDGE_LEN)
        target_qpos = aln_start_qpos + shift_needed
        for q, r in pairs:
            if q == target_qpos:
                return r
        return aln_start_rpos
    else:
        # orient=rev: pattern starts at the bridge CCC, located right after the
        # alignment in BAM SEQ.
        aln_end_qpos = pairs[-1][0] + 1
        aln_end_rpos = pairs[-1][1]
        expected_pat_start = aln_end_qpos
        search_start = max(0, expected_pat_start - 3)
        search_end = min(len(seq), expected_pat_start + 6 + _BOUNDARY_LEN)
        match = _find_boundary_match(seq, search_start, search_end, rc=True)
        if match is None:
            return aln_end_rpos
        pat_qpos, _ = match
        # tss_qpos = qpos of the position right BEFORE the CCC (= last aligned)
        tss_qpos = pat_qpos - 1
        if tss_qpos >= aln_end_qpos - 1:
            return aln_end_rpos  # alignment already ends at or before pattern start
        # Aligner extended past the boundary into bridge C's. Walk back.
        shift_needed = (aln_end_qpos - 1) - tss_qpos
        shift_needed = min(shift_needed, BRIDGE_LEN)
        target_qpos = aln_end_qpos - 1 - shift_needed
        for q, r in reversed(pairs):
            if q == target_qpos:
                return r
        return aln_end_rpos
