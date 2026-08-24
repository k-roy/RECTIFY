"""Per-read state for the cdna pipeline.

Defines :class:`ReadInfo`, the minimal per-read record carried through
clustering and consensus, plus the BAM-level extractor that builds one from a
primary alignment (UMI extraction, orient detection, walkback-canonical anchor,
TSS bridge walk-forward, full-length classification).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import pysam

from ._constants import (
    ANCHOR_FWD,
    ANCHOR_LEN,
    ANCHOR_MAX_EDIT,
    ANCHOR_RC,
    ANCHOR_SEARCH_WIN,
    ANCHOR_UPSTREAM_WIN,
    COMPLEMENT_TABLE,
    END_WINDOW_BP,
    POLY_A_ANCH_RE,
    POLY_A_UNANCH_RE,
    POLY_T_ANCH_RE,
    POLY_T_UNANCH_RE,
    SSP_FWD,
    SSP_RC,
    UMI_LEN,
)
from ..bam.bam_writer import _decode_eq_seq_inplace
from .walkback import walk_back_anchor_and_tail, walk_forward_tss

# edlib is an optional dep (loaded lazily so the module imports without it).
try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False


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
        start = locs[-1][0] if rightmost else locs[0][0]
        # edlib HW/locations can return a location whose START is None — the match END is
        # found but the start is not localizable. Guarding only on editDistance == -1 and an
        # empty `locations` list misses this, and the None then flows into the callers'
        # `if p >= 0` tests as `TypeError: '>=' not supported between 'NoneType' and 'int'`.
        # An un-localizable start carries no usable position, so report it as "no hit" using
        # this function's documented sentinel (-1), matching the find/rfind fallback below.
        # Same defect and same fix as `cdna/walkback.py::_find_adapter_anchor_pos`, where it
        # DID fire in production and crashed every ONT-cDNA run on Sherlock (planning/623).
        if start is None:
            return -1
        return start
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


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT_TABLE)[::-1]


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
        # Decode calmd '=' before the walkback reads query_sequence — the same
        # guard bam_processor.py:316 already applies on its own path. Path A
        # feeds this a bare `minimap2 | samtools sort` pre-alignment that was
        # never calmd'd, so this is a no-op there; it matters if correct-cdna is
        # ever pointed at a `rectify align` output, where SEQ is ~99.8 % '='-
        # encoded and the walkback would silently return the raw alignment end
        # with tail_len 0 (planning/578). _decode_eq_seq_inplace is idempotent
        # and returns immediately when SEQ carries no '='.
        _decode_eq_seq_inplace(read, chrom_seq_cache)
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
