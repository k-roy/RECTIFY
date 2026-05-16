#!/usr/bin/env python3
"""rectify correct-cdna — UMI-aware per-molecule consensus for PCR-cDNA.

v1: UMI extraction + Stage-1 clustering + per-cluster representative-read
    selection + Stage-1 consensus BAM output. (POA error-correction is a
    v2 upgrade; v1 uses umi_tools-style "longest-read-as-representative"
    which still does correct dedup.)

v1.11: fuzzy anchor matching via edlib HW (Lev≤2). The exact-match anchor
    finder caught only ~73% of full-length reads (matches the Q17 expected
    rate for a 13-mer with R10.4.1 hac); the ~27% miss was almost entirely
    real full-length reads with one or two basecaller errors in the adapter
    seed region. Lev≤2 lifts anchored detection to ~95–97% with effectively
    zero FP (genome+polyA combined: 0/50k random 300-bp chrI windows).

v2: Type-1↔Type-2 same-orient reconciliation by position (replaces the legacy
    cross-orient pair check, which had no biological substrate). Both physical
    PCR strands of a single molecule produce reads with the same orient label
    after BAM normalization, so cross-orient pairs cannot represent the same
    molecule. The real same-molecule signal split across populations is:
    Strand A sequenced → Type-1 (SSP captured); Strand B sequenced with SSP
    truncated at basecalled-5' → Type-2 (no UMI). v2 matches Type-1↔Type-2
    pairs at the same gene + same orient with |Δ5'|≤5 AND |Δ3'|≤5. Each pair
    gets `XL:Z:<partner cluster_id>` on both records and is emitted in
    t1t2_pairs.tsv.

v1.15: Type-2 (SSP-less) read parser. Reads where sequencing didn't reach the
    SSP/UMI region are no longer dropped — they are tagged `XT:i:2` and carried
    forward as independent observations without deduplication (no UMI available;
    reads are grouped by CPA position for isoform counting only).
    Type-1 reads (SSP+UMI captured, the standard case) carry `XT:i:1`. The two
    populations serve
    different biology:
      Type-1: full-length molecules where SSP template-switched at the original
              mRNA 5' cap → captures actual TSS.
      Type-2: 5'-truncated molecules. Likely co-translational decay
              intermediates (XRN1 paused at ribosomes → 5' positions frame-
              periodic, peak near stop codon, short polyA tail). Some technical
              truncations also live here.
    Distinguishing decay intermediates from technical noise: joint analysis
    of (Δframe, tail_length) in v1.17 (later).

v1.14: overlap-based XS classifier replaces anchor-proximity classifier. For
    each cluster's consensus alignment, find genes whose body overlaps; filter
    to candidates with gene_frac ≥ 0.3 OR read_frac ≥ 0.8. Sense candidates
    (gene strand matching implied mRNA strand from orient) are preferred over
    antisense; among sense candidates, pick the one whose TSS comes FIRST in
    the transcription direction — this attributes dicistronic readthroughs to
    the upstream-promoter gene and tolerates polyA-cleavage drift past
    annotated 3' ends without flipping to "antisense of downstream gene." Fixed
    the chrI:73,460 CDC19 misclassification where polyA drift 16 bp past
    CDC19's annotated end was being called "antisense of YAL037W."

v1.13: umi_tools-style directional UMI clustering replaces connected-components.
    Connected-components on a Lev-≤-3 graph chain-merges thousands of independent
    molecules whenever a polyA hotspot collides their UMIs within Lev-≤-3 hops;
    the v1.12 mitigation was to drop any bucket exceeding `per_cluster_cap`,
    which silently threw out ~78% of input reads on chrI. Directional clustering
    breaks transitivity via the 2× rule (a UMI can absorb a lower-count neighbor
    only if its own count is at least 2×−1 the neighbor's), so multi-read buckets
    correctly resolve into individual molecule clusters instead of one mega-cluster.
    Per-bucket cap is now advisory in the directional path.

v1.12: walk-back anchor canonicalization + sequence-level polyA tail length.
    Reads ending in polyA at a genomic A-tract had aln_end positions drifting
    by 5–50 bp across reads of the same molecule (tail-length variance +
    aligner extension into genomic A's + occasional false N-op into an
    internal A-stretch). Walk-back: starting at the basecalled-3' anchor,
    iterate get_aligned_pairs skipping ref-A (or ref-T for reversed reads),
    stop at the first matching non-A/non-T. False N-ops are handled for free
    because get_aligned_pairs omits skipped (N-op) ref positions. Tail length
    = A-count in the BAM-SEQ region between the canonical cleavage anchor
    and the adapter anchor (or 200 bp fallback). Reported per-read and
    aggregated to per-cluster median in the XA:i tag.

See ``docs/algorithms/cdna_correct.md`` (port of ont_cdna/DESIGN.md) for full
algorithmic description and empirical validation of the PCB114.24 UMI
architecture: SSP_FWD + (TT-VVVV)×4 + TTT UMI + GGG TSO bridge + cDNA body
+ polyA + adapter.

Output files (in --out directory):
  - stage1_consensus.fastq.gz    one consensus record per molecule. Per-cluster
                                  SAM-format tags are appended to each read's
                                  comment line for `rectify align -y` to pass
                                  through to the final aligned BAM.

Downstream: run `rectify align` on the FASTQ, then `rectify cdna-analyze` on
the resulting BAM to produce clusters.tsv, isoforms.tsv, and t1t2_pairs.tsv
(on post-align coordinates).

Tag glossary:
  XU  canonical UMI                          XA  polyA tail length
  XO  orient (fwd / rev)                     XG  primary gene (XS join key)
  XC  cluster size                           XI  isoform_id
  XR  input read IDs                         XT  read type (1=SSP+UMI, 2=SSP-less)
  XM  consensus method                       XB  strand-split (n_top/n_bottom)
  XS  sense/antisense/unannotated            XL  partner cluster_id (T1↔T2)
  XF  full-length tier (0/1/2)              XY  read subtype (umi_captured_fwd/rev, umi_not_captured)
  XQ  5' pre-trim bases stripped             XK  3' pre-trim bases stripped
       (SSP+UMI+GGG for T1 / polyT for rev)       (polyA for fwd / SSP_RC suffix for rev)

Algorithm versions are documented per-block above. Ported from
``ont_cdna/src/cdna_correct.py`` v1.19 + v2 (M1 dev, 2026-05-13);
argparse setup lives in ``rectify/cli.py``; the standalone entry point
lives in ont_cdna for algorithm prototyping.

Usage (via rectify CLI):
    rectify correct-cdna INPUT.bam --out OUTDIR [options]
"""

from __future__ import annotations

import argparse
import functools as _functools
import logging
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import gzip

import pysam
from intervaltree import IntervalTree
from rapidfuzz.distance import Levenshtein

from rectify.core.correct.walkback import (
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime_with_qpos,
)

# Optional dependencies (loaded lazily so the module imports without them)
try:
    import pyabpoa
    HAS_POA = True
except ImportError:
    HAS_POA = False
try:
    import mappy
    HAS_MAPPY = True
except ImportError:
    HAS_MAPPY = False
try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False

# ---- PCB114.24 chemistry constants ----------------------------------------
SSP_FWD = "TTTCTGTTGGTGCTGATATTGCT"
SSP_RC  = "AGCAATATCAGCACCAACAGAAA"
UMI_LEN = 27

# End-concordance check (XF tag): the polyA/polyT homopolymer at the OTHER end
# of the read (opposite the SSP/UMI side). For PCB114 chemistry, the basecalled
# 3' end goes through the pore first and is reliable; the basecalled 5' end is
# the truncation-prone side. The SSP/UMI side is at basecalled 5' for sense
# reads, basecalled 3' for antisense reads. So in BAM SEQ:
#   orient=fwd reads: polyA expected at RIGHT end (basecalled 3' for sense,
#                     truncation-prone basecalled 5' for antisense).
#   orient=rev reads: polyT expected at LEFT end (basecalled 5' for sense
#                     reads of - strand genes, vice versa for antisense).
# A read with the OTHER end's homopolymer detected is FULL-LENGTH (XF=1).
END_WINDOW_BP = 200            # unanchored window (last/first N bp of BAM SEQ)
MIN_HOMOPOLYMER_UNANCHORED = 10  # unanchored A/T-run threshold (XF=1 tier)
MIN_HOMOPOLYMER_ANCHORED = 6     # anchored threshold (XF=2 tier) — anchor specificity allows lower
ANCHOR_UPSTREAM_WIN = 30         # bp window upstream/downstream of anchor to test for polyA/polyT

# Adapter-start anchor — found empirically at ~72 bp from read 3' end in fwd-orient reads,
# 0 occurrences in 230 kb of chrI yeast genome (essentially zero FP combined with polyA).
ANCHOR_FWD = "TTGCGGGCGGCGG"   # fwd-orient reads: in last ~300 bp of BAM SEQ, polyA in 30 bp UPSTREAM
ANCHOR_RC  = "CCGCCGCCCGCAA"   # rev-orient reads: in first ~300 bp of BAM SEQ, polyT in 30 bp DOWNSTREAM
ANCHOR_LEN = len(ANCHOR_FWD)
ANCHOR_SEARCH_WIN = 300        # window within which to search for the anchor
# v1.11: fuzzy anchor via edlib HW mode (handles ONT R10.4.1 sub/ins/del). Lev≤2
# on 13-bp anchor → ~95–97% sensitivity (vs ~73% exact) with effectively zero FP
# when combined with the polyA-window constraint (empirical: 0/50k random 300-bp
# chrI windows pass both). Lev=3 starts to admit incidental yeast-genome hits.
ANCHOR_MAX_EDIT = 2

import re as _re
POLY_A_UNANCH_RE = _re.compile(r"A{%d,}" % MIN_HOMOPOLYMER_UNANCHORED)
POLY_T_UNANCH_RE = _re.compile(r"T{%d,}" % MIN_HOMOPOLYMER_UNANCHORED)
POLY_A_ANCH_RE   = _re.compile(r"A{%d,}" % MIN_HOMOPOLYMER_ANCHORED)
POLY_T_ANCH_RE   = _re.compile(r"T{%d,}" % MIN_HOMOPOLYMER_ANCHORED)


def _find_anchor_fuzzy(window: str, anchor: str, rightmost: bool) -> int:
    """Return position (start index within `window`) of best anchor hit at edit
    distance ≤ ANCHOR_MAX_EDIT, or -1 if no hit. Uses edlib HW (infix) mode;
    falls back to exact `find`/`rfind` if edlib is unavailable.

    edlib returns ALL co-optimal locations at the minimum edit distance found
    within the bound. For fwd-orient reads we want the RIGHTMOST occurrence
    (closest to read 3' end, i.e. the true adapter), for rev-orient the
    LEFTMOST.
    """
    if HAS_EDLIB:
        res = edlib.align(anchor, window, mode="HW", task="locations",
                          k=ANCHOR_MAX_EDIT)
        if res["editDistance"] == -1:
            return -1
        locs = res["locations"]
        if not locs:
            return -1
        # locations are (start, end_inclusive) within `window`
        return locs[-1][0] if rightmost else locs[0][0]
    # Fallback: exact match
    return window.rfind(anchor) if rightmost else window.find(anchor)


def detect_full_length_tier(seq: str, orient: str) -> int:
    """Return XF tier: 0 = not detected, 1 = unanchored only, 2 = anchored polyA at adapter.

    Tier 2 (HIGH confidence): adapter-start anchor present at expected position
    (edit distance ≤ ANCHOR_MAX_EDIT, default 2 — tolerates ONT R10.4.1 sub/ins/del)
    AND ≥6 consecutive A/T's in the 30-bp window directly adjacent to the anchor.
    Combined FP rate ≈ 0% (anchor + polyA constraint excludes random genome).

    Tier 1 (MEDIUM confidence): unanchored ≥10 A/T's in last/first 200 bp.
    For reads where even fuzzy anchor detection fails (e.g. the adapter region
    is fully truncated) but a clear homopolymer is still present.

    Tier 0: no polyA/polyT signature found.
    """
    if orient == "fwd":
        # Anchored: find rightmost fuzzy occurrence of fwd anchor in last 300 bp,
        # then test 30 bp UPSTREAM for polyA
        tail_off = max(0, len(seq) - ANCHOR_SEARCH_WIN)
        tail = seq[tail_off:]
        p_rel = _find_anchor_fuzzy(tail, ANCHOR_FWD, rightmost=True)
        if p_rel >= 0:
            anchor_abs = tail_off + p_rel
            upstream = seq[max(0, anchor_abs - ANCHOR_UPSTREAM_WIN):anchor_abs]
            if POLY_A_ANCH_RE.search(upstream):
                return 2
        # Unanchored fallback
        if POLY_A_UNANCH_RE.search(seq[-END_WINDOW_BP:]):
            return 1
        return 0
    else:  # rev
        # Anchored: find leftmost fuzzy occurrence of rev anchor in first 300 bp,
        # then test 30 bp DOWNSTREAM (immediately after the anchor) for polyT
        head = seq[:ANCHOR_SEARCH_WIN]
        p = _find_anchor_fuzzy(head, ANCHOR_RC, rightmost=False)
        if p >= 0:
            ds_start = p + ANCHOR_LEN
            downstream = seq[ds_start:ds_start + ANCHOR_UPSTREAM_WIN]
            if POLY_T_ANCH_RE.search(downstream):
                return 2
        if POLY_T_UNANCH_RE.search(seq[:END_WINDOW_BP]):
            return 1
        return 0


def detect_other_end_adapter(seq: str, orient: str) -> bool:
    """Legacy boolean shim — returns True for XF tier ≥ 1."""
    return detect_full_length_tier(seq, orient) >= 1

# ---- Sequence helpers ------------------------------------------------------
_COMPL = str.maketrans("ACGTacgtN", "TGCAtgcaN")
def revcomp(s: str) -> str:
    return s.translate(_COMPL)[::-1]


# ---- Per-read state --------------------------------------------------------
@dataclass(frozen=True)
class ReadInfo:
    """Minimal per-read state needed for clustering."""
    read_id: str
    chrom: str
    anchor: int          # canonical cleavage anchor (v1.12 walk-back) or aln-3' end (legacy)
    orient: str          # 'fwd' (SSP at LEFT of BAM SEQ) or 'rev' (SSP_RC at RIGHT)
    umi: str             # 27-nt UMI, basecalled-orient
    is_reverse: bool
    xf_tier: int         # XF: 0=not detected, 1=unanchored, 2=anchored adapter+polyA
    tail_len: int        # v1.12: sequence-level polyA tail length (A-count between cleavage anchor and adapter)
    aln_start: int       # v1.14: aligned region start in ref coords (for overlap-based XS classification)
    aln_end: int         # v1.14: aligned region end (exclusive) in ref coords
    read_type: int       # v1.15: 1 = SSP+UMI captured, 2 = SSP-less (5'-truncated, e.g. decay intermediate)
    pos5_corrected: int  # v1.19: TSS-side position corrected for SSP-bridge G-tract ambiguity (analog of 3' polyA walk-back)
    read_subtype: str    # "umi_captured_fwd" (Type-1, SSP+UMI at 5') / "umi_captured_rev" (Type-1, SSP+UMI at 3' via pA-first traversal) / "umi_not_captured" (Type-2, pA-first truncated before UMI)


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
        return off + r["locations"][-1][0]
    win = seq[:ANCHOR_SEARCH_WIN]
    r = edlib.align(ANCHOR_RC, win, mode="HW", task="locations", k=ANCHOR_MAX_EDIT)
    if r["editDistance"] == -1 or not r["locations"]: return None
    return r["locations"][0][0]


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


BRIDGE_LEN = 3     # the rGrGrG of the TSO that becomes the GGG bridge in BAM SEQ.
# v1.19 chemistry-aware 5' boundary detection:
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
#
# Detection: at each candidate 12-bp window near the alignment boundary, score
# the match against the canonical pattern. Pick the highest-scoring window
# with score >= threshold and 3-of-3 G's (or C's for rev) intact — those G's
# ARE the bridge and we use them to anchor the boundary.

# Pattern score templates: position → {allowed bases}
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


def extract_read_info(read: pysam.AlignedSegment,
                       chrom_seq_cache: Optional[Dict[str, str]] = None
                       ) -> Optional[ReadInfo]:
    """Extract UMI + anchor + orient + tail_len from a primary alignment.

    If `chrom_seq_cache` is provided (recommended), uses v1.12 walk-back to
    produce a canonical cleavage anchor and per-read polyA tail length. Without
    a cache, falls back to the polyA-side aln position with tail_len=0.

    ANCHOR LOGIC (v1.12): the anchor lives at the polyA cleavage end of the
    molecule in ref coords. Walk-back (when ref cache available) canonicalizes
    this across tail-length variance, aligner extension into genomic A-tracts,
    and false-intron N-op artifacts.

    Direction of "polyA side" is determined by **orient**, not is_reverse:
      orient=fwd → BAM SEQ = SSP-UMI-cDNA-polyA-adapter → polyA at RIGHT → anchor near aln_end
      orient=rev → BAM SEQ = adapter_RC-polyT-cDNA_RC-UMI_RC-SSP_RC → polyT at LEFT → anchor near aln_start
    The same molecule sequenced as Strand A vs Strand B can produce reads with
    different is_reverse but identical orient (depending on gene strand); using
    orient guarantees both land on the same polyA cleavage genomic position.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return None
    seq = read.query_sequence
    if seq is None:
        return None

    # ---- Type 1 detection (SSP+UMI captured) ----
    read_type = 1
    umi_basecalled: Optional[str] = None
    orient: Optional[str] = None
    p = seq.find(SSP_FWD)
    if p >= 0:
        umi_basecalled = seq[p + len(SSP_FWD): p + len(SSP_FWD) + UMI_LEN]
        if len(umi_basecalled) == UMI_LEN:
            orient = "fwd"
        else:
            umi_basecalled = None
    if orient is None:
        p = seq.find(SSP_RC)
        if p >= UMI_LEN:
            umi_rc = seq[p - UMI_LEN: p]
            umi_basecalled = revcomp(umi_rc)
            if len(umi_basecalled) == UMI_LEN:
                orient = "rev"
            else:
                umi_basecalled = None

    # ---- Type 2 fallback (SSP truncated; infer orient from polyA/adapter) ----
    if orient is None:
        # Test both BAM-SEQ ends for adapter+polyA pattern. detect_full_length_tier
        # returns 2 (anchored) or 1 (unanchored polyA only) when the pattern is
        # present at the orient-appropriate end. We try both orients and accept
        # the one with a tier ≥ 1; ties → prefer orient with anchored (tier 2).
        tier_fwd = detect_full_length_tier(seq, "fwd")
        tier_rev = detect_full_length_tier(seq, "rev")
        if tier_fwd == 0 and tier_rev == 0:
            return None
        if tier_fwd >= tier_rev:
            orient = "fwd"
        else:
            orient = "rev"
        read_type = 2
        umi_basecalled = ""  # placeholder for Type-2 (no UMI available)

    if chrom_seq_cache is not None and read.reference_name in chrom_seq_cache:
        chrom_seq = chrom_seq_cache[read.reference_name]
        anchor, tail_len = walk_back_anchor_and_tail(read, chrom_seq, orient)
        # v1.19: 5' TSS walk-forward (Type-1 only — Type-2 5' is the truncation
        # point, not the SSP-bridge boundary).
        if read_type == 1:
            pos5_corrected = walk_forward_tss(read, orient, chrom_seq)
        else:
            pos5_corrected = read.reference_start if orient == "fwd" else (read.reference_end or 1) - 1
    else:
        # Fallback: polyA-side aln position (orient-based, not is_reverse-based).
        # orient=fwd → polyA at RIGHT of BAM SEQ → aln_end-1.
        # orient=rev → polyT at LEFT of BAM SEQ → aln_start.
        anchor = ((read.reference_end or 0) - 1) if orient == "fwd" else read.reference_start
        tail_len = 0
        pos5_corrected = read.reference_start if orient == "fwd" else (read.reference_end or 1) - 1

    if read_type == 1:
        read_subtype = "umi_captured_fwd" if orient == "fwd" else "umi_captured_rev"
    else:
        read_subtype = "umi_not_captured"

    return ReadInfo(
        read_id=read.query_name,
        chrom=read.reference_name,
        anchor=anchor,
        orient=orient,
        umi=umi_basecalled or "",
        is_reverse=read.is_reverse,
        xf_tier=detect_full_length_tier(seq, orient),
        tail_len=tail_len,
        aln_start=read.reference_start,
        aln_end=(read.reference_end or read.reference_start + 1),
        read_type=read_type,
        pos5_corrected=pos5_corrected,
        read_subtype=read_subtype,
    )


# ---- Clustering ------------------------------------------------------------
def umi_components(umis: List[str], max_edit: int) -> List[List[int]]:
    """Connected-components clustering of UMIs by Levenshtein ≤ max_edit.

    LEGACY (pre-v1.13). Chains independent molecules together when many UMIs
    in a bucket are within Lev≤max_edit (the "polyA-hotspot chain-merge"
    pathology). Retained for comparison / debugging via --umi-clustering=components.
    """
    n = len(umis)
    if n == 0:
        return []
    parent = list(range(n))
    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    for i in range(n):
        for j in range(i + 1, n):
            if Levenshtein.distance(umis[i], umis[j], score_cutoff=max_edit) <= max_edit:
                union(i, j)
    comps: defaultdict[int, list] = defaultdict(list)
    for i in range(n):
        comps[find(i)].append(i)
    return list(comps.values())


def umi_components_directional(umis: List[str], max_edit: int) -> List[List[int]]:
    """umi_tools-style directional clustering (Smith, Heger & Sudbery 2017, Genome Res).

    Build a directed graph over UNIQUE UMIs in the bucket:
      node a → b iff Lev(a, b) ≤ max_edit AND count(a) ≥ 2*count(b) − 1.

    Then walk from the highest-count unvisited node via outgoing edges (BFS);
    everything reached forms one molecule cluster. Repeat with next highest-count
    unvisited node.

    The 2× rule constrains transitive merging: a small-count UMI (likely an error
    variant of a deeper UMI nearby) gets absorbed by its parent, but cannot itself
    absorb other UMIs of comparable depth. This breaks the chain-merge that
    connected-components suffers at polyA hotspots where thousands of independent
    molecules' UMIs land within Lev ≤ 3 of each other.

    For singleton-dominated buckets (every UMI has count 1, the 2× rule reduces
    to "1 ≥ 1", allowing edges between any two Lev-≤-max_edit singletons), the
    behavior collapses back to connected-components — that's fine; directional
    helps where it's needed (multi-read buckets with mixed coverage levels).
    """
    n = len(umis)
    if n == 0: return []
    if n == 1: return [[0]]

    # Group read indices by exact UMI → unique UMI list + counts.
    umi_to_read_indices: Dict[str, List[int]] = defaultdict(list)
    for i, u in enumerate(umis):
        umi_to_read_indices[u].append(i)
    unique = list(umi_to_read_indices.keys())
    counts = [len(umi_to_read_indices[u]) for u in unique]
    n_uniq = len(unique)

    # Order unique UMIs by count desc (ties broken by UMI lex order for determinism).
    order_desc = sorted(range(n_uniq), key=lambda i: (-counts[i], unique[i]))

    # Build adjacency: for each parent a, find children b later in order_desc whose
    # count satisfies the 2× rule, then Lev-check.
    adj: Dict[int, List[int]] = defaultdict(list)
    for ai_pos, a in enumerate(order_desc):
        ca = counts[a]
        max_cb = (ca + 1) // 2  # need cb ≤ this for 2× rule
        # counts in order_desc are non-increasing; skip until first b with cb ≤ max_cb.
        bi_start = ai_pos + 1
        while bi_start < n_uniq and counts[order_desc[bi_start]] > max_cb:
            bi_start += 1
        ua = unique[a]
        for bi_pos in range(bi_start, n_uniq):
            b = order_desc[bi_pos]
            d = Levenshtein.distance(ua, unique[b], score_cutoff=max_edit)
            if d <= max_edit:
                adj[a].append(b)

    # BFS from each unvisited highest-count root.
    visited: set = set()
    clusters: List[List[int]] = []
    for root in order_desc:
        if root in visited: continue
        component_uniqs = []
        queue = [root]
        while queue:
            u = queue.pop()
            if u in visited: continue
            visited.add(u)
            component_uniqs.append(u)
            for v in adj[u]:
                if v not in visited:
                    queue.append(v)
        # Expand back to read indices.
        read_indices: List[int] = []
        for uid in component_uniqs:
            read_indices.extend(umi_to_read_indices[unique[uid]])
        clusters.append(read_indices)
    return clusters


def position_components_directional(positions: List[int], weights: List[int],
                                     max_dist: int) -> List[List[int]]:
    """v1.17: directional clustering on integer positions for isoform grouping.

    Parallel to umi_components_directional but with two differences:
      1. Distance = abs(pos[a] - pos[b]) (integer), not Levenshtein.
      2. STRICT 2× rule: edge a→b iff w(a) ≥ 2·w(b) AND |pos[a]-pos[b]| ≤ max_dist.
         (Not the umi_tools "n ≥ 2m − 1" — strict prevents bridge-chains where
         count=1 valley positions link two real peaks through a chain of edges.)

    `positions` and `weights` are parallel arrays — one entry per input cluster.
    Multiple input clusters at the same position get their weights summed and
    are treated as one "node" for the directional graph; this is exactly the
    behavior we want (a peak's strength = total read support, not count of
    sub-clusters).

    Returns list of components, each as a list of indices into the input arrays.
    """
    n = len(positions)
    if n == 0: return []
    if n == 1: return [[0]]

    # Group input indices by exact position; weight per node = SUM of input weights.
    pos_to_indices: Dict[int, List[int]] = defaultdict(list)
    for i, p in enumerate(positions):
        pos_to_indices[p].append(i)
    unique = sorted(pos_to_indices.keys())  # sorted by position (helps locality)
    node_weight = [sum(weights[i] for i in pos_to_indices[p]) for p in unique]
    n_nodes = len(unique)

    # Order nodes by weight desc; ties broken by position for determinism.
    order_desc = sorted(range(n_nodes), key=lambda i: (-node_weight[i], unique[i]))

    # Build adjacency. STRICT 2× rule + max_dist.
    adj: Dict[int, List[int]] = defaultdict(list)
    for ai_pos, a in enumerate(order_desc):
        wa = node_weight[a]
        pa = unique[a]
        # Iterate later positions in count-desc order. Need w(b) ≤ wa/2 (strict)
        # AND |pa - pb| ≤ max_dist.
        for bi_pos in range(ai_pos + 1, n_nodes):
            b = order_desc[bi_pos]
            if node_weight[b] * 2 > wa:
                continue  # strict 2× rule fails
            if abs(unique[b] - pa) > max_dist:
                continue
            adj[a].append(b)

    # BFS from each unvisited highest-weight root.
    visited: set = set()
    clusters: List[List[int]] = []
    for root in order_desc:
        if root in visited: continue
        component_nodes = []
        queue = [root]
        while queue:
            u = queue.pop()
            if u in visited: continue
            visited.add(u)
            component_nodes.append(u)
            for v in adj[u]:
                if v not in visited:
                    queue.append(v)
        # Expand back to input-cluster indices.
        cluster_indices: List[int] = []
        for nid in component_nodes:
            cluster_indices.extend(pos_to_indices[unique[nid]])
        clusters.append(cluster_indices)
    return clusters


def position_components(reads: List[ReadInfo]) -> List[List[int]]:
    """v1.15: exact (aln_start, aln_end) grouping for Type-2 (SSP-less) reads.

    No UMI is available for Type-2 reads, so we collapse exact-position duplicates
    only. This is conservative — molecules sharing the SAME aln_start AND aln_end
    are treated as PCR/sequencing duplicates of one decay/truncation event.
    Wobble of even ±1 bp produces a new cluster; that's fine because real distinct
    decay-intermediate species typically have well-defined termini.
    """
    pos_to_indices: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, r in enumerate(reads):
        pos_to_indices[(r.aln_start, r.aln_end)].append(i)
    return list(pos_to_indices.values())


def cluster_reads(reads: List[ReadInfo], anchor_window: int, max_edit: int,
                  per_cluster_cap: int,
                  clustering_method: str = "directional"
                  ) -> Tuple[List[List[ReadInfo]], dict]:
    """Stage-1 clustering: anchor-bucket reads, then UMI-cluster within each bucket.

    clustering_method:
      - "directional" (v1.13 default): umi_tools-style; breaks polyA-hotspot
        chain-merge. Per-bucket cap is advisory only (large buckets logged but
        still processed).
      - "components" (legacy): connected-components; needs `per_cluster_cap` to
        drop oversized buckets to avoid catastrophic chain-merge.
    """
    # Bucket key includes read_type (v1.15) — Type-1 and Type-2 reads never
    # mix within a bucket, so they get different clustering algorithms.
    anchor_buckets: dict = defaultdict(list)
    for r in reads:
        bucket = (r.chrom, r.anchor // anchor_window, r.orient, r.read_type)
        anchor_buckets[bucket].append(r)

    stats = dict(
        anchor_buckets=len(anchor_buckets),
        buckets_dropped_polyA_pileup=0,
        reads_in_dropped_buckets=0,
        biggest_bucket_size=0,
        molecule_clusters=0,
        molecule_clusters_size_1=0,
        molecule_clusters_size_2=0,
        molecule_clusters_size_ge_3=0,
        type1_clusters=0,
        type2_clusters=0,
        type1_reads=sum(1 for r in reads if r.read_type == 1),
        type2_reads=sum(1 for r in reads if r.read_type == 2),
    )

    umi_cluster_fn = umi_components_directional if clustering_method == "directional" else umi_components

    out_clusters: List[List[ReadInfo]] = []
    for bucket_key, bucket_reads in anchor_buckets.items():
        read_type = bucket_key[3]
        bsize = len(bucket_reads)
        if bsize > stats["biggest_bucket_size"]:
            stats["biggest_bucket_size"] = bsize
        # Legacy components path still respects the cap to avoid chain-merge (Type-1 only).
        if read_type == 1 and clustering_method == "components" and bsize > per_cluster_cap:
            stats["buckets_dropped_polyA_pileup"] += 1
            stats["reads_in_dropped_buckets"] += bsize
            continue
        if read_type == 1:
            umis = [r.umi for r in bucket_reads]
            comps = umi_cluster_fn(umis, max_edit)
        else:
            comps = position_components(bucket_reads)
        for comp in comps:
            cluster = [bucket_reads[i] for i in comp]
            out_clusters.append(cluster)
            sz = len(cluster)
            if sz == 1: stats["molecule_clusters_size_1"] += 1
            elif sz == 2: stats["molecule_clusters_size_2"] += 1
            else: stats["molecule_clusters_size_ge_3"] += 1
            if read_type == 1: stats["type1_clusters"] += 1
            else: stats["type2_clusters"] += 1
    stats["molecule_clusters"] = len(out_clusters)
    return out_clusters, stats


# ---- v1: per-cluster consensus ---------------------------------------------
def pick_representative(reads_in_cluster: List[pysam.AlignedSegment]) -> pysam.AlignedSegment:
    """Choose the cluster representative (umi_tools-style longest-and-best).

    Used for singletons (size 1) and as fallback when pileup consensus fails.
    """
    def key(r: pysam.AlignedSegment) -> tuple:
        ref_len = (r.reference_end or 0) - r.reference_start
        return (-ref_len, -r.mapping_quality, r.query_name)
    return min(reads_in_cluster, key=key)


def pileup_consensus(reads_in_cluster: List[pysam.AlignedSegment]) -> Optional[Tuple[str, List[Tuple[int, int]], int, str]]:
    """Build pileup-based consensus across a multi-read cluster (v1.5).

    For each reference position covered by ≥ ceil(N/2) reads in the cluster, the
    consensus base is the majority across reads; positions with low M-coverage
    (insufficient reads have a matched base there) are emitted as 'D' (deletion)
    in the consensus CIGAR.

    Returns (consensus_seq, cigartuples, ref_start, qual_string) or None if
    insufficient data. cigartuples uses pysam op codes (0=M, 2=D).

    Limitations (deferred to v2 with real POA):
    - Insertions in input reads are dropped (not represented in consensus).
    - N-op intron skips are emitted as D (functional but not technically correct
      for spliced reads). For yeast this is fine since <5% of genes have introns.
    - Per-base quality is uniform Q30 (could be improved with depth weighting).
    """
    if not reads_in_cluster:
        return None

    n_reads = len(reads_in_cluster)
    base_counts: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    coverage: Dict[int, int] = defaultdict(int)

    span_start = min(r.reference_start for r in reads_in_cluster)
    span_end = max((r.reference_end or 0) for r in reads_in_cluster)
    if span_end <= span_start:
        return None

    for read in reads_in_cluster:
        seq = read.query_sequence
        if seq is None:
            continue
        # Walk M-aligned pairs (excludes I, S, H, N, D)
        for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
            if r_pos is None or q_pos is None:
                continue
            base = seq[q_pos].upper()
            base_counts[r_pos][base] += 1
            coverage[r_pos] += 1

    # Build consensus sequence + CIGAR walking the full span
    seq_chars: List[str] = []
    cigartuples: List[Tuple[int, int]] = []  # (op_code, length)
    threshold = max(1, (n_reads + 1) // 2)  # majority floor

    for r_pos in range(span_start, span_end):
        cov = coverage.get(r_pos, 0)
        if cov >= threshold:
            # Consensus base = most-frequent base; ties broken by alphabetical (deterministic)
            bdict = base_counts[r_pos]
            consensus_base = max(bdict.items(), key=lambda kv: (kv[1], -ord(kv[0])))[0]
            seq_chars.append(consensus_base)
            op = 0  # M
        else:
            op = 2  # D — most reads have no base here (intron skip or genuine deletion)

        if cigartuples and cigartuples[-1][0] == op:
            cigartuples[-1] = (op, cigartuples[-1][1] + 1)
        else:
            cigartuples.append((op, 1))

    if not seq_chars:
        return None

    # Trim trailing D from CIGAR (no SAM convention forbids them but they're
    # uninformative; cleaner output)
    while cigartuples and cigartuples[-1][0] == 2:
        cigartuples.pop()
    while cigartuples and cigartuples[0][0] == 2:
        # Leading D shifts ref_start
        leading_len = cigartuples[0][1]
        span_start += leading_len
        cigartuples.pop(0)

    consensus_seq = "".join(seq_chars)
    qual_string = "?" * len(consensus_seq)  # Q30 across (uniform; v2 can depth-weight)
    return consensus_seq, cigartuples, span_start, qual_string


def poa_consensus_from_strings(seqs: List[str]) -> Optional[str]:
    """Run POA on raw sequence strings (already in matching orientation).

    Used both by `poa_consensus` (reads-from-BAM path) and by the v1.18
    strand-aware consensus path (which feeds two per-strand sub-consensuses
    in for a final merge).
    """
    if not HAS_POA or len(seqs) < 2:
        return None
    aligner = pyabpoa.msa_aligner()
    res = aligner.msa(seqs, out_cons=True, out_msa=False)
    if not res.cons_seq:
        return None
    return res.cons_seq[0]


def poa_consensus(reads_in_cluster: List[pysam.AlignedSegment]) -> Optional[str]:
    """Build POA consensus sequence across cluster reads (v2 upgrade over pileup).

    Operates on the BAM SEQ (in reference orientation) of each read. POA handles
    indels in homopolymer regions much better than pileup. Returns the consensus
    as a string, or None if abPOA isn't available or the cluster has too few reads.

    The consensus needs to be re-aligned to the reference afterward (use mappy)
    to obtain a valid CIGAR — the POA itself doesn't produce reference alignment.
    """
    seqs = [r.query_sequence for r in reads_in_cluster if r.query_sequence]
    return poa_consensus_from_strings(seqs)


def poa_consensus_strand_aware(reads_in_cluster: List[pysam.AlignedSegment]
                                ) -> Optional[str]:
    """v1.18: strand-split-then-merge POA consensus.

    Split cluster reads by is_reverse:
      - is_reverse=True reads (top-strand sequencing) have their reliable
        basecalled-3' mapped to the LEFT of BAM SEQ (SSP/UMI side).
      - is_reverse=False reads (bottom-strand sequencing) have their reliable
        basecalled-3' mapped to the RIGHT of BAM SEQ (polyA side).

    Each sub-pool has correlated strand-specific error modes. Building a
    per-strand sub-consensus first, then POA-merging the two sub-consensuses,
    cancels strand-specific systematic biases that POA-on-all-reads would
    leave on the consensus.

    Falls back to single-strand POA when one side has too few reads (<2).
    """
    if not HAS_POA or len(reads_in_cluster) < 2:
        return None
    top_seqs = [r.query_sequence for r in reads_in_cluster
                if r.is_reverse and r.query_sequence]
    bot_seqs = [r.query_sequence for r in reads_in_cluster
                if (not r.is_reverse) and r.query_sequence]

    # If one strand has only 0-1 reads, can't usefully build a sub-consensus.
    # Fall back to the all-reads POA on the other strand (or both pooled).
    if len(top_seqs) < 2 and len(bot_seqs) < 2:
        return poa_consensus_from_strings(top_seqs + bot_seqs)
    if len(top_seqs) < 2:
        return poa_consensus_from_strings(bot_seqs)
    if len(bot_seqs) < 2:
        return poa_consensus_from_strings(top_seqs)

    # Both strands have enough reads — build sub-consensus per strand and merge.
    top_cons = poa_consensus_from_strings(top_seqs)
    bot_cons = poa_consensus_from_strings(bot_seqs)
    if top_cons is None and bot_cons is None: return None
    if top_cons is None: return bot_cons
    if bot_cons is None: return top_cons
    # Both sub-consensuses are in BAM-SEQ orientation, so POA them directly.
    return poa_consensus_from_strings([top_cons, bot_cons])


def pretrim_consensus(consensus_seq: str, orient: str, read_type: int
                      ) -> Tuple[str, int, int]:
    """Strip SSP/UMI/GGG from 5' and poly-A/adapter from 3' of a consensus sequence.

    This ensures aligners receive clean mRNA sequence — no adapter contamination
    at either end. Downstream `rectify align` re-aligns the trimmed consensus
    and the aligner produces appropriate soft-clips for the stripped bases.

    5' trim (Type-1 only — Type-2 has no SSP/UMI):
      - orient=fwd: SSP_FWD (23nt) + UMI (27nt) + GGG bridge (3nt) = 53 nt
      - orient=rev: SSP_RC at the RIGHT of BAM SEQ + UMI_RC + CCC bridge
        For rev, the mRNA 5' end is at the RIGHT; the poly-A (= polyT in BAM SEQ)
        is at the LEFT. We strip the CCC+AAA+UMI_RC+SSP_RC suffix.

    3' trim (both types):
      - orient=fwd: poly-A is at the RIGHT of BAM SEQ. Strip from the start of
        the first homopolymer A run of length >= MIN_HOMOPOLYMER_ANCHORED.
      - orient=rev: poly-T is at the LEFT. Strip from position 0 up to the end
        of the first homopolymer T run.

    Returns:
      (trimmed_seq, q_trim_5, pretrim_pa_len) where:
        trimmed_seq    — the sequence passed to aligners
        q_trim_5       — bases stripped from 5' (for restoring soft-clip)
        pretrim_pa_len — A's stripped from 3' (for XA fallback if walkback fails)
    """
    q_trim_5 = 0
    pretrim_pa_len = 0

    # ---- 5' trim (SSP/UMI/GGG) ----
    if read_type == 1:
        if orient == "fwd":
            p = consensus_seq.find(SSP_FWD)
            if p >= 0:
                q_trim_5 = p + len(SSP_FWD) + UMI_LEN + BRIDGE_LEN
            else:
                # SSP not found in consensus (degenerate) — skip 5' trim
                q_trim_5 = 0
        else:  # orient == "rev"
            p = consensus_seq.find(SSP_RC)
            if p >= UMI_LEN + BRIDGE_LEN:
                # SSP_RC found; trim everything from (p - UMI_LEN - BRIDGE_LEN) onward
                # The mRNA body occupies consensus_seq[0 : p - UMI_LEN - BRIDGE_LEN]
                # We express this as a suffix trim (handled below as rev 3' trim equiv).
                # For rev: the "5'" of the mRNA is at the RIGHT, and "5' trim" means
                # removing from the right end. We encode as q_trim_5 on rev in terms of
                # the sequence orientation: suffix_len = len(consensus_seq) - (p - UMI_LEN - BRIDGE_LEN)
                # We use a separate variable to keep the logic clear:
                pass  # handled in the 3' trim block for rev

    # ---- 3' trim (poly-A/adapter) ----
    if orient == "fwd":
        # For fwd: poly-A is at the RIGHT of BAM SEQ. Strip from last A-run.
        # Search for the poly-A start using the adapter anchor first (Tier 2),
        # fall back to unanchored poly-A detection (Tier 1).
        seq_to_search = consensus_seq[q_trim_5:]
        adp_pos = _find_adapter_anchor_pos(seq_to_search, "fwd")
        if adp_pos is not None:
            # Poly-A is in the 30 bp UPSTREAM of the adapter anchor
            upstream_start = max(0, adp_pos - ANCHOR_UPSTREAM_WIN)
            upstream = seq_to_search[upstream_start:adp_pos]
            m = POLY_A_ANCH_RE.search(upstream)
            if m:
                pa_start = q_trim_5 + upstream_start + m.start()
                pretrim_pa_len = m.group(0).count("A")
                trimmed_seq = consensus_seq[q_trim_5:pa_start]
                return (trimmed_seq, q_trim_5, pretrim_pa_len)
        # Unanchored fallback: find rightmost A-run (>= MIN_HOMOPOLYMER_ANCHORED) in
        # last END_WINDOW_BP bases. Using the anchored threshold (6+) rather than
        # the conservative unanchored threshold (10+) because in a cDNA consensus
        # sequence we know the poly-A is at the 3' end. Taking the RIGHTMOST match
        # avoids stripping on internal A-runs (e.g. gene body A's).
        window = seq_to_search[-END_WINDOW_BP:] if len(seq_to_search) > END_WINDOW_BP else seq_to_search
        window_off = len(seq_to_search) - len(window)
        all_pa_matches = list(POLY_A_ANCH_RE.finditer(window))
        if all_pa_matches:
            m = all_pa_matches[-1]  # rightmost = closest to 3' end
            pa_start = q_trim_5 + window_off + m.start()
            pretrim_pa_len = m.group(0).count("A")
            trimmed_seq = consensus_seq[q_trim_5:pa_start]
            return (trimmed_seq, q_trim_5, pretrim_pa_len)
        # No poly-A found: return 5'-trimmed only
        trimmed_seq = consensus_seq[q_trim_5:]
        return (trimmed_seq, q_trim_5, 0)
    else:  # orient == "rev"
        # For rev: poly-T is at the LEFT of BAM SEQ (basecalled poly-A in +ref coords).
        # The mRNA body is to the RIGHT; the SSP_RC/UMI_RC/CCC are also to the RIGHT.
        # Strip poly-T from left, and SSP_RC suffix from right.

        # First: find poly-T at left end (= 3' trim in rev orientation).
        seq = consensus_seq
        adp_pos = _find_adapter_anchor_pos(seq, "rev")
        if adp_pos is not None:
            # poly-T is in 30 bp DOWNSTREAM of the anchor (adapter at LEFT, polyT just after)
            ds_end = adp_pos + ANCHOR_LEN + ANCHOR_UPSTREAM_WIN
            downstream = seq[adp_pos + ANCHOR_LEN: min(ds_end, len(seq))]
            m = POLY_T_ANCH_RE.search(downstream)
            if m:
                pt_end = adp_pos + ANCHOR_LEN + m.end()
                pretrim_pa_len = downstream[m.start():m.end()].count("T")
                seq = seq[pt_end:]
                q_trim_5 = pt_end
        else:
            # Unanchored: find first T-run in first END_WINDOW_BP bases
            window = seq[:END_WINDOW_BP]
            m = POLY_T_UNANCH_RE.search(window)
            if m:
                pt_end = m.end()
                pretrim_pa_len = window[m.start():m.end()].count("T")
                seq = seq[pt_end:]
                q_trim_5 = pt_end

        # Second: strip SSP_RC/UMI_RC/CCC suffix from the RIGHT (rev 5' trim)
        p = seq.find(SSP_RC)
        if p >= UMI_LEN + BRIDGE_LEN:
            ssp_trim_start = p - UMI_LEN - BRIDGE_LEN  # start of CCC in the seq
            seq = seq[:ssp_trim_start]
        # Note: q_trim_5 for rev tracks left-side trim; right-side trim is not
        # encoded in q_trim_5. The right suffix trim is implicit in the trimmed_seq
        # length: suffix_trim = (len(consensus_seq) - q_trim_5) - len(seq).

        return (seq, q_trim_5, pretrim_pa_len)


def write_stage1_fastq(input_bam: Path, output_fastq: Path,
                       clusters: List[List[ReadInfo]],
                       umi_canonical: Dict[int, str],
                       cluster_xf_tier: Dict[int, int],
                       cluster_tail_len: Dict[int, int],
                       use_poa: bool = False,
                       strand_aware_consensus: bool = False) -> dict:
    """Emit per-cluster consensus sequences as a FASTQ for downstream alignment.

    `rectify align` will run the multi-aligner on this FASTQ to produce the final
    aligned BAM. Per-cluster SAM-format tags are appended to each read's comment
    line so `minimap2 -y` (and equivalent flags on mapPacBio/gapmm2) can
    propagate them through to the output BAM.

    Only alignment-independent tags are emitted here. Gene/strand/isoform/pair
    tags (XG, XS, XI, XL) and any pre-align poly-A length are recomputed by
    `rectify cdna-analyze` on the post-align CIGAR, so emitting them here would
    be misleading by the time downstream consumers read them. Tail length is
    still emitted as XA for upstream debuggability; cdna-analyze overwrites it.

    Read naming: `cluster_<cid>`.
    """
    import gzip

    # Build read_id → cluster_id index
    rid_to_cluster: Dict[str, int] = {}
    cluster_read_ids: Dict[int, List[str]] = {}
    for cid, c in enumerate(clusters):
        cluster_read_ids[cid] = [r.read_id for r in c]
        for r in c:
            rid_to_cluster[r.read_id] = cid

    cluster_strand_split: Dict[int, Tuple[int, int]] = {}
    for cid, c in enumerate(clusters):
        n_top = sum(1 for r in c if r.is_reverse)
        n_bot = sum(1 for r in c if not r.is_reverse)
        cluster_strand_split[cid] = (n_top, n_bot)

    # Bucket reads per cluster from input BAM
    cluster_segments: Dict[int, List[pysam.AlignedSegment]] = defaultdict(list)
    n_in = 0
    with pysam.AlignmentFile(str(input_bam), "rb") as inbam:
        for read in inbam.fetch(until_eof=True):
            n_in += 1
            cid = rid_to_cluster.get(read.query_name)
            if cid is None:
                continue
            cluster_segments[cid].append(read)

    is_gz = str(output_fastq).endswith(".gz")
    opener = gzip.open if is_gz else open

    written = singleton = multi_pileup = multi_fallback = 0

    with opener(str(output_fastq), "wt") as fq:
        for cid, c in enumerate(clusters):
            segs = cluster_segments.get(cid, [])
            if not segs:
                continue

            method = "rep"
            seq: Optional[str] = None

            if len(segs) == 1:
                seg = segs[0]
                seq = seg.query_sequence or ""
                # Restore basecalled orientation: BAM SEQ is RC'd for minus-strand
                # alignments; we want the original read sequence for re-alignment.
                if seg.is_reverse:
                    seq = seq[::-1].translate(_COMPL)
                singleton += 1
            else:
                if use_poa:
                    poa_fn = poa_consensus_strand_aware if strand_aware_consensus else poa_consensus
                    seq = poa_fn(segs)
                    if seq:
                        method = "poa"
                if not seq:
                    pileup = pileup_consensus(segs)
                    if pileup is not None:
                        seq = pileup[0]
                        method = "pileup"
                if not seq:
                    rep = pick_representative(segs)
                    seq = rep.query_sequence or ""
                    if rep.is_reverse:
                        seq = seq[::-1].translate(_COMPL)
                    method = "rep_fallback"
                    multi_fallback += 1
                else:
                    multi_pileup += 1

            if not seq:
                continue

            # Strip adapter/UMI/GGG prefix and poly-A/T suffix so rectify align
            # receives clean mRNA sequence. Store strip lengths as XQ/XK tags so
            # the aligner can reconstruct the full-length CIGAR with soft-clips.
            orient = c[0].orient
            read_type = c[0].read_type
            trimmed_seq, q_trim_5, _ = pretrim_consensus(seq, orient, read_type)
            if not trimmed_seq:
                trimmed_seq = seq
                q_trim_5 = 0
            q_trim_3 = len(seq) - q_trim_5 - len(trimmed_seq)

            n_top, n_bot = cluster_strand_split[cid]
            tag_parts = [
                f"XU:Z:{umi_canonical[cid]}",
                f"XO:Z:{orient}",
                f"XC:i:{len(c)}",
                f"XR:Z:{','.join(cluster_read_ids[cid])}",
                f"XM:Z:{method}",
                f"XF:i:{cluster_xf_tier[cid]}",
                f"XA:i:{cluster_tail_len.get(cid, 0)}",
                f"XT:i:{read_type}",
                f"XY:Z:{c[0].read_subtype}",
                f"XQ:i:{q_trim_5}",
                f"XK:i:{q_trim_3}",
                f"XB:Z:{n_top}/{n_bot}",
            ]

            # Tab-separate the read name and tags so `minimap2 -y` (and the
            # equivalent flags in mapPacBio / gapmm2) parse each `XX:T:value`
            # as a separate SAM aux field. Space separators silently collapse
            # the comment into a single Z-string aux on the first tag.
            comment = "\t".join(tag_parts)
            qual = "?" * len(trimmed_seq)
            fq.write(f"@cluster_{cid}\t{comment}\n{trimmed_seq}\n+\n{qual}\n")
            written += 1

    return dict(input_reads=n_in, written=written, from_singletons=singleton,
                from_multi_pileup=multi_pileup, from_multi_fallback=multi_fallback)


def canonical_umi(cluster_umis: List[str]) -> str:
    """Pick the canonical UMI for a cluster (most frequent; ties → first lexicographically)."""
    counts = defaultdict(int)
    for u in cluster_umis:
        counts[u] += 1
    return min(counts, key=lambda u: (-counts[u], u))


# ---- v1.17: isoform clustering --------------------------------------------
def assign_isoforms(clusters: List[List[ReadInfo]],
                     cluster_xg: Dict[int, Optional[str]],
                     tol5: int, tol3: int
                     ) -> Tuple[Dict[int, str], List[dict]]:
    """v1.17: group Stage-1 clusters into isoform clusters.

    Within each (gene, orient, read_type) group:
      Type-1 (SSP+UMI captured, 5' = real mRNA 5' end):
        - directional clustering on 5' positions (tol5, strict 2× rule)
        - directional clustering on 3' positions (tol3, strict 2× rule)
        - isoform = (gene, orient, type=1, 5'-group, 3'-group)
      Type-2 (SSP-less, 5' is random truncation point):
        - skip 5' clustering (positions are noisy)
        - directional clustering on 3' positions (tol3, strict 2× rule)
        - isoform = (gene, orient, type=2, "*", 3'-group)

    Weights = cluster.n_reads (PCR-collapsed support).

    Returns:
      cluster_to_isoform_id: dict cluster_id → isoform_id string
      isoform_records: list of dicts (one per isoform) with aggregated stats
    """
    cluster_to_iso: Dict[int, str] = {}
    iso_records: List[dict] = []

    # Group cluster indices by (gene, orient, read_type)
    groups: Dict[Tuple[str, str, int], List[int]] = defaultdict(list)
    for cid, c in enumerate(clusters):
        gene = cluster_xg.get(cid)
        if gene is None: continue  # skip unannotated
        groups[(gene, c[0].orient, c[0].read_type)].append(cid)

    for (gene, orient, rt), cids in groups.items():
        # 5'/3' position per cluster (mRNA-relative coords).
        # orient=fwd: 5' = aln_start, 3' = aln_end-1.
        # orient=rev: 5' = aln_end-1, 3' = aln_start.
        pos5: List[int] = []
        pos3: List[int] = []
        weights: List[int] = []
        for cid in cids:
            r0 = clusters[cid][0]
            # v1.19: use pos5_corrected (bridge walk-forward applied for Type-1;
            # raw aln-side position for Type-2).
            pos5.append(r0.pos5_corrected)
            if orient == "fwd":
                pos3.append(r0.aln_end - 1)
            else:
                pos3.append(r0.aln_start)
            weights.append(len(clusters[cid]))

        # 3' clustering (used by both Type 1 and Type 2)
        comps3 = position_components_directional(pos3, weights, tol3)
        cid_to_3g = {}
        for g3_idx, comp in enumerate(comps3):
            for local_idx in comp:
                cid_to_3g[cids[local_idx]] = g3_idx

        if rt == 1:
            comps5 = position_components_directional(pos5, weights, tol5)
            cid_to_5g = {}
            for g5_idx, comp in enumerate(comps5):
                for local_idx in comp:
                    cid_to_5g[cids[local_idx]] = g5_idx
            # Build isoform IDs and group clusters by (5'-group, 3'-group)
            iso_buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
            for cid in cids:
                iso_buckets[(cid_to_5g[cid], cid_to_3g[cid])].append(cid)
            for (g5, g3), iso_cids in iso_buckets.items():
                iso_id = f"{gene}_t1_5g{g5}_3g{g3}"
                _emit_isoform(iso_id, gene, orient, 1, iso_cids, clusters,
                              cluster_to_iso, iso_records, has_5g=True)
        else:
            # Type-2: only 3' grouping
            iso_buckets3: Dict[int, List[int]] = defaultdict(list)
            for cid in cids:
                iso_buckets3[cid_to_3g[cid]].append(cid)
            for g3, iso_cids in iso_buckets3.items():
                iso_id = f"{gene}_t2_3g{g3}"
                _emit_isoform(iso_id, gene, orient, 2, iso_cids, clusters,
                              cluster_to_iso, iso_records, has_5g=False)

    return cluster_to_iso, iso_records


def reconcile_t1_t2_pairs(clusters: List[List[ReadInfo]],
                            cluster_xg: Dict[int, Optional[str]],
                            tol5: int, tol3: int
                            ) -> Tuple[Dict[int, int], List[dict]]:
    """v2: link Type-1 ↔ Type-2 clusters of the SAME molecule (same gene + same orient
    + both 5' and 3' termini within tol bp).

    Biology: a single dsDNA molecule produces some reads from Strand A (SSP captured
    at reliable basecalled-3' → Type-1, UMI extractable) and some from Strand B
    where SSP at basecalled-5' gets truncated → Type-2 (no UMI). The two read
    populations end up in different Stage-1 clusters (Type-1 in UMI cluster,
    Type-2 in position cluster). Reconciliation rescues the link.

    Cross-orient pairs (the v1.13 Stage-2 idea) don't exist in PCB114 chemistry —
    both physical PCR strands of the same molecule give the same orient label
    after BAM normalization. The Type-1↔Type-2 same-orient pair IS the real
    biological substrate.

    Matching is one-to-one greedy: each Type-2 cluster pairs with its nearest
    Type-1 candidate (closest sum |Δ5'|+|Δ3'|). Each Type-1 can be claimed by
    at most one Type-2 (first-match-by-distance wins).

    Returns:
      cluster_xl: dict cluster_id → partner_cluster_id (bidirectional links)
      pair_records: list of dicts (one per pair) for the TSV manifest
    """
    cluster_xl: Dict[int, int] = {}
    pair_records: List[dict] = []

    # Group cluster indices by (gene, orient)
    groups: Dict[Tuple[str, str], Dict[int, List[int]]] = defaultdict(
        lambda: {1: [], 2: []})
    for cid, c in enumerate(clusters):
        gene = cluster_xg.get(cid)
        if gene is None: continue
        groups[(gene, c[0].orient)][c[0].read_type].append(cid)

    for (gene, orient), by_type in groups.items():
        t1_cids = by_type[1]
        t2_cids = by_type[2]
        if not t1_cids or not t2_cids:
            continue
        # Build (5', 3') for each
        def termini(cid: int) -> Tuple[int, int]:
            r0 = clusters[cid][0]
            # v1.19 pos5_corrected (bridge walk-forward); 3' is raw polyA-side
            if orient == "fwd":
                return (r0.pos5_corrected, r0.aln_end - 1)
            return (r0.pos5_corrected, r0.aln_start)

        t1_termini = {cid: termini(cid) for cid in t1_cids}
        t2_termini = {cid: termini(cid) for cid in t2_cids}

        # Greedy match: for each Type-2, find closest available Type-1 within tol
        claimed_t1: set = set()
        # Sort Type-2 by n_reads desc (claim higher-coverage Type-2 first)
        t2_sorted = sorted(t2_cids, key=lambda c: -len(clusters[c]))
        for t2 in t2_sorted:
            p5_2, p3_2 = t2_termini[t2]
            best_t1 = None
            best_dist = None
            for t1 in t1_cids:
                if t1 in claimed_t1: continue
                p5_1, p3_1 = t1_termini[t1]
                d5 = abs(p5_1 - p5_2)
                d3 = abs(p3_1 - p3_2)
                if d5 > tol5 or d3 > tol3: continue
                dist = d5 + d3
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best_t1 = t1
            if best_t1 is None: continue
            claimed_t1.add(best_t1)
            cluster_xl[best_t1] = t2   # Type-1 → Type-2 partner cid
            cluster_xl[t2] = best_t1   # Type-2 → Type-1 partner cid
            p5_1, p3_1 = t1_termini[best_t1]
            pair_records.append({
                "t1_cid": best_t1, "t2_cid": t2,
                "gene": gene, "orient": orient,
                "t1_pos5": p5_1, "t1_pos3": p3_1,
                "t2_pos5": p5_2, "t2_pos3": p3_2,
                "d5": abs(p5_1 - p5_2), "d3": abs(p3_1 - p3_2),
                "t1_n_reads": len(clusters[best_t1]),
                "t2_n_reads": len(clusters[t2]),
                "t1_umi": umi_canon_placeholder(clusters[best_t1]),
            })
    return cluster_xl, pair_records


def umi_canon_placeholder(cluster: List[ReadInfo]) -> str:
    """Return the canonical UMI of a Type-1 cluster (most-frequent UMI string).
    For Type-2 clusters this returns empty string."""
    counts: defaultdict[str, int] = defaultdict(int)
    for r in cluster:
        counts[r.umi] += 1
    if not counts: return ""
    return min(counts, key=lambda u: (-counts[u], u))


def _emit_isoform(iso_id: str, gene: str, orient: str, read_type: int,
                   iso_cids: List[int], clusters: List[List[ReadInfo]],
                   cluster_to_iso: Dict[int, str], iso_records: List[dict],
                   has_5g: bool) -> None:
    pos5_vals = []; pos3_vals = []; n_reads_total = 0
    tail_lens: List[int] = []
    chrom = clusters[iso_cids[0]][0].chrom
    read_subtype = clusters[iso_cids[0]][0].read_subtype
    for cid in iso_cids:
        cluster_to_iso[cid] = iso_id
        r0 = clusters[cid][0]
        pos5_vals.append(r0.pos5_corrected)
        pos3_vals.append(r0.aln_end - 1 if orient == "fwd" else r0.aln_start)
        n_reads_total += len(clusters[cid])
        # Each read's tail_len weighted by appearing in the cluster
        for r in clusters[cid]:
            tail_lens.append(r.tail_len)
    def _modal(xs):
        if not xs: return -1
        c = defaultdict(int)
        for x in xs: c[x] += 1
        return max(c, key=lambda k: (c[k], -k))
    def _median(xs):
        if not xs: return 0
        s = sorted(xs); return s[len(s) // 2]
    iso_records.append({
        "isoform_id": iso_id,
        "gene": gene,
        "chrom": chrom,
        "orient": orient,
        "read_type": read_type,
        "read_subtype": read_subtype,
        "n_clusters": len(iso_cids),
        "n_reads_total": n_reads_total,
        "pos5_modal": _modal(pos5_vals) if has_5g else -1,
        "pos3_modal": _modal(pos3_vals),
        "tail_len_median": _median(tail_lens),
        "cluster_ids": ",".join(str(c) for c in iso_cids),
    })


# ---- v1.6: gene-strand–aware sense/antisense classification ----------------
def _parse_gff_gene_name(attrs: str) -> Optional[str]:
    """Extract a stable gene identifier from a GFF attributes column.

    Preference order (matches SGD GFF conventions + cross-project handoff needs):
      1. ID=        — systematic name (e.g. YAL038W); the canonical join key
      2. Name=      — same as ID for SGD
      3. gene=      — common name (e.g. CDC19); fallback only
    Returns None if no identifier found.
    """
    for kv in attrs.split(";"):
        if kv.startswith("ID="):
            return kv[3:].strip()
    for kv in attrs.split(";"):
        if kv.startswith("Name="):
            return kv[5:].strip()
    for kv in attrs.split(";"):
        if kv.startswith("gene="):
            return kv[5:].strip()
    return None


def load_rdna_intervals(gff_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """Parse rRNA_gene features from GFF → {chrom: [(start, end), ...]} (0-based half-open).

    Used for rDNA masking: reads overlapping these intervals are dropped in
    stream_reads() before UMI extraction and clustering, preventing the O(n²)
    Levenshtein bottleneck that occurs when thousands of rDNA reads share
    identical anchor sequences (observed: 261k chrXII reads → 4h42m clustering).
    """
    result: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "rRNA_gene":
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GFF is 1-based inclusive → 0-based half-open
            end = int(parts[4])
            result[chrom].append((start, end))
    return dict(result)


def load_gff_genes(gff_path: Path) -> Dict[str, IntervalTree]:
    """Load gene-like features → chrom → IntervalTree.

    Each interval's data is a tuple `(strand, gene_name)` for both XS
    classification and downstream join-key emission (XG tag, v1.14).
    """
    trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    # Gene-level types only (no mRNA/transcript — those expose transcript IDs
    # like "YAL038W_id001" which break the systematic-name join key used by
    # cross-project taxonomies. SGD GFF reliably emits a 'gene' feature at
    # every protein-coding locus.)
    GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene",
                   "tRNA_gene", "snoRNA_gene", "snRNA_gene", "rRNA_gene"}
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in GENE_TYPES:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            if start < end and strand in ("+", "-"):
                name = _parse_gff_gene_name(parts[8]) or f"{chrom}:{start}-{end}({strand})"
                trees[chrom][start:end] = (strand, name)
    return dict(trees)


def classify_sense_antisense(chrom: str, aln_start: int, aln_end: int, orient: str,
                              gene_trees: Dict[str, IntervalTree],
                              min_gene_frac: float = 0.3,
                              min_read_frac: float = 0.8
                              ) -> Tuple[str, Optional[str]]:
    """v1.14: Overlap-based gene attribution with sense-strand preference and first-TSS tiebreaker.

    For each gene whose body overlaps the cluster's alignment:
      - gene_frac = overlap / gene_length  (how much of the gene is covered)
      - read_frac = overlap / read_aln_length  (how much of the read is in the gene)

    Keep candidates passing `gene_frac >= min_gene_frac` OR `read_frac >= min_read_frac`.
    Among these, classify by sense-strand match: the implied mRNA strand is "+"
    for orient=fwd, "-" for orient=rev.

      - If any candidate gene's strand matches mRNA strand → SENSE. Tie-break
        by picking the gene whose TSS is FIRST in the transcription direction
        (lowest TSS for + reads, highest TSS for - reads). This naturally
        attributes dicistronic readthroughs to the upstream-promoter gene.
      - Else if any candidate gene exists on the opposite strand → ANTISENSE.
        Tie-break by maximum gene_frac.
      - Else → unannotated.

    Replaces the legacy anchor-proximity classifier (pre-v1.14), which mis-
    classified polyA-drift reads as antisense of downstream genes whose CAP
    start happened to lie just past the read's drift cleavage site.
    """
    tree = gene_trees.get(chrom)
    if tree is None:
        return ("unannotated", None)
    aln_len = aln_end - aln_start
    if aln_len <= 0:
        return ("unannotated", None)

    overlaps = list(tree.overlap(aln_start, aln_end))
    if not overlaps:
        return ("unannotated", None)

    mRNA_strand = "+" if orient == "fwd" else "-"

    # (gene_frac, read_frac, tss, gene_strand, gene_name)
    candidates: List[Tuple[float, float, int, str, str]] = []
    for iv in overlaps:
        ov_len = max(0, min(iv.end, aln_end) - max(iv.begin, aln_start))
        if ov_len <= 0:
            continue
        gene_len = iv.end - iv.begin
        gene_frac = ov_len / gene_len
        read_frac = ov_len / aln_len
        if gene_frac < min_gene_frac and read_frac < min_read_frac:
            continue
        gene_strand, gene_name = iv.data
        tss = iv.begin if gene_strand == "+" else iv.end
        candidates.append((gene_frac, read_frac, tss, gene_strand, gene_name))

    if not candidates:
        return ("unannotated", None)

    # Sense-preferred: pick first by transcription direction TSS.
    sense_candidates = [c for c in candidates if c[3] == mRNA_strand]
    if sense_candidates:
        # First TSS in transcription direction = lowest for +, highest for -.
        if mRNA_strand == "+":
            sense_candidates.sort(key=lambda c: c[2])
        else:
            sense_candidates.sort(key=lambda c: -c[2])
        return ("sense", sense_candidates[0][4])

    # Antisense fallback: pick the gene whose body the read overlaps most.
    anti_candidates = [c for c in candidates if c[3] != mRNA_strand]
    if anti_candidates:
        anti_candidates.sort(key=lambda c: -c[0])  # max gene_frac
        return ("antisense", anti_candidates[0][4])

    return ("unannotated", None)


# ---- I/O -------------------------------------------------------------------
def stream_reads(bam_path: Path, region: Optional[str],
                 reference: Optional[Path] = None,
                 rdna_intervals: Optional[Dict[str, List[Tuple[int, int]]]] = None,
                 ) -> Tuple[List[ReadInfo], int]:
    """Stream UMI-extractable reads; return (reads, n_rdna_masked).

    rdna_intervals: {chrom: [(start, end), ...]} from load_rdna_intervals().
    Reads whose alignment overlaps any rRNA_gene interval are skipped before
    UMI extraction, preventing the O(n²) clustering bottleneck on rDNA chroms.
    """
    chrom_cache: Optional[Dict[str, str]] = None
    if reference is not None:
        chrom_cache = {}
        with pysam.FastaFile(str(reference)) as fa:
            for c in fa.references:
                chrom_cache[c] = fa.fetch(c).upper()
    reads: List[ReadInfo] = []
    n_masked = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        it = bam.fetch(region=region) if region else bam.fetch(until_eof=True)
        for read in it:
            if rdna_intervals and read.reference_name in rdna_intervals:
                rs = read.reference_start
                re = read.reference_end or rs
                if any(s < re and rs < e for s, e in rdna_intervals[read.reference_name]):
                    n_masked += 1
                    continue
            info = extract_read_info(read, chrom_cache)
            if info is not None:
                reads.append(info)
    return reads, n_masked


# ---- Main ------------------------------------------------------------------
def run(args) -> int:
    """v3 RECTIFY integration entry point. The argparse Namespace is built by
    rectify/cli.py's `correct-cdna` subparser; we coerce string paths to Path
    objects here (the rectify CLI uses bare `default=None` strings, while the
    algorithm code below expects pathlib.Path)."""
    args.bam = Path(args.bam) if not isinstance(args.bam, Path) else args.bam
    args.out = Path(args.out) if not isinstance(args.out, Path) else args.out
    if getattr(args, "reference", None) is not None:
        args.reference = Path(args.reference) if not isinstance(args.reference, Path) else args.reference
    for attr, default in [
        ("umi_clustering", "directional"),
        ("strand_aware_consensus", False),
        ("no_mask_rdna", False),
    ]:
        if not hasattr(args, attr): setattr(args, attr, default)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("cdna_correct")

    if not args.bam.exists():
        log.error("BAM not found: %s", args.bam)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    log.info("Streaming reads from %s region=%s ...", args.bam, args.region or "all")
    if args.reference is None:
        log.warning("No --reference provided: walk-back anchor canonicalization + "
                    "polyA tail-length measurement DISABLED (legacy aln-end anchor used)")

    # rDNA masking (default ON): skip reads overlapping rRNA_gene intervals to
    # avoid O(n²) UMI clustering on chrXII rDNA tandem repeats.
    rdna_intervals: Optional[Dict[str, List[Tuple[int, int]]]] = None
    if args.gff and not args.no_mask_rdna:
        rdna_intervals = load_rdna_intervals(args.gff)
        n_rdna_loci = sum(len(v) for v in rdna_intervals.values())
        if n_rdna_loci:
            log.info("rDNA masking ON: %d rRNA_gene loci across %d chrom(s) "
                     "(--no-mask-rdna to disable)", n_rdna_loci, len(rdna_intervals))

    reads, n_rdna_masked = stream_reads(args.bam, args.region,
                                        reference=args.reference,
                                        rdna_intervals=rdna_intervals)
    log.info("  %d reads with extractable UMIs (%.1fs)", len(reads), time.time() - t0)
    if n_rdna_masked:
        log.info("  rDNA masked: %d reads skipped (rRNA_gene overlap)", n_rdna_masked)

    t1 = time.time()
    log.info("Stage 1: clustering (anchor window %d, UMI Lev≤%d, method=%s)",
             args.anchor_window_bp, args.umi_edit_distance, args.umi_clustering)
    clusters, stats = cluster_reads(reads, args.anchor_window_bp,
                                     args.umi_edit_distance, args.per_cluster_cap,
                                     clustering_method=args.umi_clustering)
    log.info("  %d molecule clusters (%.1fs)  biggest bucket=%d reads",
             stats["molecule_clusters"], time.time() - t1, stats["biggest_bucket_size"])

    # Canonical UMI per cluster
    umi_canon = {cid: canonical_umi([r.umi for r in c]) for cid, c in enumerate(clusters)}

    # Per-cluster XF tier: max of constituent reads (any one full-length is enough)
    cluster_xf_tier = {cid: max((r.xf_tier for r in c), default=0) for cid, c in enumerate(clusters)}

    # Per-cluster polyA tail length: median across reads (sequence-level, biased
    # short for long tails — see v1.12 docstring). Singletons report the lone
    # read's per-read tail length.
    def _median(xs: List[int]) -> int:
        if not xs: return 0
        s = sorted(xs)
        return s[len(s) // 2]
    cluster_tail_len = {cid: _median([r.tail_len for r in c])
                        for cid, c in enumerate(clusters)}

    # Gene assignment (XG/XS), isoform clustering (XI), and Type-1↔Type-2
    # pairing (XL) moved to `rectify cdna-analyze`. They run on post-align
    # coordinates from `rectify align`, which are more accurate than the
    # pre-align positions available here. correct-cdna now emits ONLY the
    # alignment-independent stage-1 outputs: the consensus FASTQ.

    # Stage 1 FASTQ — per-cluster consensus sequences for `rectify align` to
    # remap with the triple aligner. Per-cluster tags ride along on the FASTQ
    # comment line so the downstream BAM picks them up via `minimap2 -y`.
    t2 = time.time()
    out_fastq = args.out / "stage1_consensus.fastq.gz"
    log.info("Writing Stage-1 consensus FASTQ → %s", out_fastq)
    # POA still uses pyabpoa for sequence-level consensus; no aligner needed here.
    use_poa = HAS_POA and not args.no_poa

    if use_poa:
        log.info("Multi-read cluster consensus: POA (pyabpoa)")
    else:
        log.info("Multi-read cluster consensus: pileup-style (substitution-corrective only)")

    fastq_stats = write_stage1_fastq(args.bam, out_fastq, clusters, umi_canon,
                                      cluster_xf_tier, cluster_tail_len,
                                      use_poa=use_poa,
                                      strand_aware_consensus=args.strand_aware_consensus)
    log.info("  wrote %d records (%d singletons, %d pileup consensus, %d rep fallback) in %.1fs",
             fastq_stats["written"], fastq_stats["from_singletons"],
             fastq_stats["from_multi_pileup"], fastq_stats["from_multi_fallback"],
             time.time() - t2)

    # Final report
    print()
    print("=" * 70)
    print("cdna_correct v1 — Stage-1 dedup complete")
    print("=" * 70)
    print(f"input reads (UMI-extractable): {len(reads):>8d}")
    print(f"output records (one per molecule): {fastq_stats['written']:>8d}"
          f"  ({100 * fastq_stats['written'] / max(1, len(reads)):.1f}% of input)")
    print(f"  singletons (passed through):     {fastq_stats['from_singletons']:>8d}")
    print(f"  multi-read (pileup consensus):   {fastq_stats['from_multi_pileup']:>8d}")
    print(f"  multi-read (rep fallback):       {fastq_stats['from_multi_fallback']:>8d}")
    print(f"polyA-pileup buckets dropped:    {stats['buckets_dropped_polyA_pileup']:>8d}"
          f"  ({stats['reads_in_dropped_buckets']} reads — these need POA + position-aware handling)")
    n_t1 = stats.get('type1_reads', 0); n_t2 = stats.get('type2_reads', 0)
    n_total_typed = max(1, n_t1 + n_t2)
    print(f"Read-type breakdown:")
    print(f"  Type 1 (SSP+UMI captured):       {n_t1:>8d}  ({100*n_t1/n_total_typed:.1f}%)  → {stats.get('type1_clusters', 0)} clusters")
    print(f"  Type 2 (SSP-less, 5'-truncated): {n_t2:>8d}  ({100*n_t2/n_total_typed:.1f}%)  → {stats.get('type2_clusters', 0)} clusters")
    print()
    print("XF tier breakdown (full-length confidence):")
    tier_counts = defaultdict(int)
    for v in cluster_xf_tier.values():
        tier_counts[v] += 1
    n_clust = max(1, len(clusters))
    print(f"  XF=2 (anchored, HIGH confidence):     {tier_counts[2]:>8d} clusters"
          f"  ({100*tier_counts[2]/n_clust:.1f}%)")
    print(f"  XF=1 (unanchored, MEDIUM confidence): {tier_counts[1]:>8d} clusters"
          f"  ({100*tier_counts[1]/n_clust:.1f}%)")
    print(f"  XF=0 (not detected):                  {tier_counts[0]:>8d} clusters"
          f"  ({100*tier_counts[0]/n_clust:.1f}%)")
    n_full = tier_counts[1] + tier_counts[2]
    print(f"  XF≥1 (any full-length):               {n_full:>8d} clusters"
          f"  ({100*n_full/n_clust:.1f}%)")
    # PolyA tail length distribution (sequence-level; bias-uniform across conditions)
    print()
    print("PolyA tail-length distribution (per-cluster median, sequence-level, XF=2 only):")
    bins = [(0, 1), (1, 10), (10, 20), (20, 30), (30, 50), (50, 75),
            (75, 100), (100, 150), (150, 250), (250, 10_000)]
    bin_counts = [0] * len(bins)
    for cid in range(len(clusters)):
        if cluster_xf_tier[cid] < 2: continue
        tl = cluster_tail_len[cid]
        for i, (lo, hi) in enumerate(bins):
            if lo <= tl < hi:
                bin_counts[i] += 1; break
    n_with_tl = sum(bin_counts)
    print(f"  (restricted to XF=2 anchored clusters: N={n_with_tl})")
    print(f"  {'range':<12} {'count':>8} {'pct':>6}")
    for (lo, hi), n in zip(bins, bin_counts):
        label = f"{lo}-{hi-1}" if hi < 10_000 else f"≥{lo}"
        pct = 100 * n / max(1, n_with_tl)
        print(f"  {label:<12} {n:>8d} {pct:>5.1f}%")
    print()
    print(f"Output FASTQ: {out_fastq}")
    print("Next step: `rectify align` on this FASTQ → `rectify cdna-analyze` "
          "for clusters.tsv / isoforms.tsv / t1t2_pairs.tsv.")
    log.info("Total runtime: %.1fs", time.time() - t0)
    return 0


# Module is invoked via rectify CLI: `rectify correct-cdna ...` calls run(args).
# No __main__ block — the standalone entry point lives in the ont_cdna staging
# repo at /Users/kevinroy/work/ont_cdna/src/cdna_correct.py.


def create_correct_cdna_parser(subparsers):
    """Wire the `correct-cdna` subcommand into the given subparsers group."""
    import argparse
    # =========================================================================
    # correct-cdna command (UMI-aware per-molecule consensus from cDNA BAM)
    # =========================================================================
    correct_cdna_parser = subparsers.add_parser(
        'correct-cdna',
        help='UMI-aware per-molecule consensus for PCR-cDNA BAMs (PCB114.24 chemistry)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            'Cluster aligned PCR-cDNA reads (SQK-PCB114.24 SSP + 27-nt UMI architecture) '
            'by (locus, orientation, UMI) and emit one consensus record per starting RNA '
            'molecule.\n\n'
            'Operates AFTER alignment, complementing trim-cdna-polya which runs on the '
            'FASTQ before alignment. Emits a representative-read or pileup-based '
            'consensus per cluster (POA if pyabpoa is available).\n\n'
            'Output: stage1_consensus.fastq.gz — one consensus sequence per UMI cluster, '
            'with alignment-independent SAM-tag comments (XU/XO/XC/XR/XM/XF/XA/XT/XY/XQ/XK/XB) '
            'for `rectify align -y` to propagate into the post-align BAM. Gene assignment, '
            'isoform clustering, and Type-1↔Type-2 pairing run downstream in '
            '`rectify cdna-analyze`.'
        ),
    )
    correct_cdna_parser.add_argument(
        'bam',
        help='Input BAM (aligned PCR-cDNA reads, indexed)',
    )
    correct_cdna_parser.add_argument(
        '-o', '--out',
        dest='out',
        required=True,
        help='Output directory (will contain stage1_consensus.fastq.gz)',
    )
    correct_cdna_parser.add_argument(
        '--umi-edit-distance',
        dest='umi_edit_distance',
        type=int,
        default=3,
        help='Max Levenshtein distance between UMIs in the same cluster',
    )
    correct_cdna_parser.add_argument(
        '--anchor-window-bp',
        dest='anchor_window_bp',
        type=int,
        default=25,
        help='Window around alignment-end anchor for same-locus clustering (bp)',
    )
    correct_cdna_parser.add_argument(
        '--per-cluster-cap',
        dest='per_cluster_cap',
        type=int,
        default=200,
        help='Max reads per cluster (larger = potential PCR-jackpot signal, slower POA)',
    )
    correct_cdna_parser.add_argument(
        '--gff',
        default=None,
        help='Genome annotation GFF3 (gzip OK). Used for rDNA masking (reads '
             'overlapping rRNA_gene loci are excluded by default to prevent the '
             'O(n²) UMI clustering blow-up on chrXII tandem repeats). Sense/antisense '
             'XS classification, XG gene-name tagging, and isoform clustering now '
             'run downstream in `rectify cdna-analyze`.',
    )
    correct_cdna_parser.add_argument(
        '--reference',
        default=None,
        help='Reference FASTA (gzip OK). Required for walk-back anchor '
             'canonicalization + poly-A tail-length measurement during UMI '
             'extraction; without it, the legacy aln-end anchor is used.',
    )
    correct_cdna_parser.add_argument(
        '--no-poa',
        action='store_true',
        dest='no_poa',
        default=False,
        help='Force pileup-only consensus even if abPOA is available',
    )
    correct_cdna_parser.add_argument(
        '--region',
        default=None,
        help='Restrict to a single BAM region (e.g. chrI) — useful for testing',
    )
    correct_cdna_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        default=False,
        help='Verbose DEBUG-level logging',
    )
    # v1.13+ flags (added during v3 RECTIFY integration):
    correct_cdna_parser.add_argument(
        '--umi-clustering',
        dest='umi_clustering',
        choices=('directional', 'components'),
        default='directional',
        help='UMI clustering method (default: directional, umi_tools-style 2× rule)',
    )
    correct_cdna_parser.add_argument(
        '--strand-aware-consensus',
        dest='strand_aware_consensus',
        action='store_true',
        default=False,
        help='v1.18 strand-aware POA: split reads by is_reverse, build a per-strand '
             'sub-consensus, then merge. Cancels strand-specific systematic errors.',
    )
    correct_cdna_parser.add_argument(
        '--no-mask-rdna',
        dest='no_mask_rdna',
        action='store_true',
        default=False,
        help='Disable rDNA masking. By default, reads overlapping rRNA_gene loci in '
             '--gff are excluded before clustering to prevent the O(n²) UMI '
             'bottleneck on chrXII tandem repeats (observed: 261k reads → 4h42m).',
    )
