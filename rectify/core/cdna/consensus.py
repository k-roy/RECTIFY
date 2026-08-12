"""Per-cluster consensus building + 5'/3' adapter pretrim for the cdna pipeline.

Three consensus paths:
  * :func:`pileup_consensus` (v1.5) — majority-vote per ref position; emits a
    pysam-style CIGAR. Drops insertions, treats N-op intron skips as deletions.
    Fine for yeast (<5% intron-containing genes) and is the always-available
    fallback when abPOA isn't installed.
  * :func:`poa_consensus` (v2) — partial-order alignment via abPOA. Handles
    homopolymer indels much better than pileup. Returns a sequence string only;
    requires re-alignment downstream for a valid CIGAR.
  * :func:`poa_consensus_strand_aware` (v1.18) — split reads by ``is_reverse``,
    build a per-strand sub-consensus, then POA-merge the two sub-consensuses to
    cancel strand-specific systematic errors.

:func:`pretrim_consensus` strips the SSP/UMI/GGG 5' prefix and the polyA/adapter
3' suffix from a finished consensus so downstream ``rectify align`` receives
clean mRNA sequence.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, NamedTuple, Optional, Tuple

import pysam

from ._constants import (
    ANCHOR_LEN,
    ANCHOR_UPSTREAM_WIN,
    BRIDGE_LEN,
    END_WINDOW_BP,
    POLY_A_ANCH_RE,
    POLY_T_ANCH_RE,
    POLY_T_UNANCH_RE,
    SSP_FWD,
    SSP_MAX_EDIT,
    SSP_RC,
    SSP_SEARCH_WIN,
    UMI_LEN,
)
from .walkback import _find_adapter_anchor_pos

# Optional dependencies (loaded lazily so the module imports without them).
try:
    import pyabpoa
    HAS_POA = True
except ImportError:
    HAS_POA = False

try:
    import edlib
    HAS_EDLIB = True
except ImportError:
    HAS_EDLIB = False


def restore_eq_seq(seg: pysam.AlignedSegment,
                   fasta: "Optional[pysam.FastaFile]") -> Optional[str]:
    """Return ``seg``'s query sequence with calmd ``=`` placeholders resolved.

    ``rectify align`` runs ``samtools calmd -e`` (see ``_apply_calmd_eq``), which
    rewrites every M base matching the reference to ``=`` in SEQ. pysam's
    ``query_sequence`` then returns ``=`` for those positions, so re-emitting it
    to the Stage-1 consensus FASTQ yields an unmappable read — this silently
    dropped ~99% of cDNA molecules at Stage-C align2. The original match base is
    NOT recoverable from the BAM (MD stores match *run lengths*, not bases), so
    it is restored from the *reference genome*: for each M/=/X position that
    reads ``=``, substitute the reference base. Soft-clips / insertions /
    mismatches keep their stored (real) base. Returns the BAM-oriented sequence
    (same orientation as ``query_sequence``); if ``fasta`` is None, SEQ has no
    ``=``, or the read is unmapped, the raw sequence is returned unchanged.
    """
    seq = seg.query_sequence
    if not seq or '=' not in seq or fasta is None or seg.is_unmapped:
        return seq
    ref = fasta.fetch(seg.reference_name, seg.reference_start, seg.reference_end).upper()
    chars = list(seq)
    qi = 0  # index into query (chars)
    ri = 0  # index into fetched reference span
    for op, ln in (seg.cigartuples or []):
        if op in (0, 7, 8):        # M / = / X : consume query + reference
            for k in range(ln):
                if chars[qi + k] == '=':
                    chars[qi + k] = ref[ri + k]
            qi += ln
            ri += ln
        elif op in (1, 4):         # I / S : consume query only
            qi += ln
        elif op in (2, 3):         # D / N : consume reference only
            ri += ln
        # H (5) / P (6): consume neither
    return ''.join(chars)


def pileup_consensus(reads_in_cluster: List[pysam.AlignedSegment], fasta: "Optional[pysam.FastaFile]" = None) -> Optional[Tuple[str, List[Tuple[int, int]], int, str]]:
    """Build pileup-based consensus across a multi-read cluster (v1.5; 653 quality-weighted).

    For each reference position covered by ≥ ceil(N/2) reads in the cluster, the
    consensus base is the phred-quality-weighted majority across reads (raw count
    when quality strings are absent); positions with low M-coverage (insufficient
    reads have a matched base there) are emitted as 'D' (deletion) in the
    consensus CIGAR.

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
        seq = restore_eq_seq(read, fasta)
        if seq is None:
            continue
        # 653: weight each vote by the base's phred quality instead of +1. PCR-cDNA read pools
        # mix the two sequencing-time regimes (planning/650: pore-entry bases run ~1.7x the
        # error of pore-exit bases, a gradient along each read), so on an even strand split a
        # raw-count vote is a coin flip while the quality-weighted vote resolves toward the
        # strand that read this REGION at pore exit. Reads without quality strings vote with
        # weight 1 (the pre-653 behaviour); `coverage` stays a read count, so the majority
        # coverage threshold below is unchanged.
        quals = read.query_qualities
        for q_pos, r_pos in read.get_aligned_pairs(matches_only=True):
            if r_pos is None or q_pos is None:
                continue
            base = seq[q_pos].upper()
            base_counts[r_pos][base] += quals[q_pos] if quals else 1
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


def poa_consensus(reads_in_cluster: List[pysam.AlignedSegment], fasta: "Optional[pysam.FastaFile]" = None) -> Optional[str]:
    """Build POA consensus sequence across cluster reads (v2 upgrade over pileup).

    Operates on the BAM SEQ (in reference orientation) of each read. POA handles
    indels in homopolymer regions much better than pileup. Returns the consensus
    as a string, or None if abPOA isn't available or the cluster has too few reads.

    The consensus needs to be re-aligned to the reference afterward (use mappy)
    to obtain a valid CIGAR — the POA itself doesn't produce reference alignment.
    """
    seqs = [s for r in reads_in_cluster if (s := restore_eq_seq(r, fasta))]
    return poa_consensus_from_strings(seqs)


def poa_consensus_strand_aware(reads_in_cluster: List[pysam.AlignedSegment],
                                fasta: "Optional[pysam.FastaFile]" = None
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
    top_seqs = [s for r in reads_in_cluster
                if r.is_reverse and (s := restore_eq_seq(r, fasta))]
    bot_seqs = [s for r in reads_in_cluster
                if (not r.is_reverse) and (s := restore_eq_seq(r, fasta))]

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


class PretrimResult(NamedTuple):
    """What :func:`pretrim_consensus` stripped, in both coordinate systems.

    The old 3-tuple ``(trimmed_seq, q_trim_5, pretrim_pa_len)`` conflated two different
    things: on the ``fwd`` branch ``q_trim_5`` was the mRNA 5' trim, but on the ``rev``
    branch it was the LEFT (poly-T = mRNA **3'**) trim, and the mRNA 5' trim was not
    returned at all — the function's own comment admitted "right-side trim is not encoded
    in q_trim_5". The caller then derived ``XK`` by subtraction, so ``XQ``/``XK`` were
    **swapped on every rev molecule**. Nothing consumes those tags yet (planning/679 CP9a),
    so this widens the contract while it is still free to do so.

    ``seq == input[left_trim : len(input) - right_trim]`` always holds.
    """
    seq: str            # the trimmed sequence handed to the aligner
    left_trim: int      # bases removed from the LEFT of the input (frame-relative)
    right_trim: int     # bases removed from the RIGHT of the input (frame-relative)
    trim_5p: int        # bases removed from the mRNA 5' end  → XQ
    trim_3p: int        # bases removed from the mRNA 3' end  → XK
    pa_len: int         # A's stripped with the poly-A tail (XA fallback)
    frame: str          # the frame actually trimmed in ('fwd' / 'rev')
    frame_mismatch: bool  # the sequence's frame disagreed with the `orient` label


def _find_ssp(seq: str, frame: str) -> int:
    """Return the start index of the SSP in `seq` for `frame`, or -1.

    ``fwd`` → ``SSP_FWD``, expected near the LEFT; ``rev`` → ``SSP_RC``, near the RIGHT
    (hence ``rfind``: the rightmost hit is the real adapter).

    Exact search runs first and is unwindowed, preserving the historical behaviour — an
    exact 23-mer arises by chance at ~4⁻²³. The **fuzzy** fallback is deliberately
    WINDOW-GATED to the end where the SSP belongs (planning/681): edlib ``HW`` over a full
    ~2 kb consensus finds a ≤3-edit 23-mer in ordinary genomic sequence and would trim real
    mRNA — worst on the pileup molecules, which carry no adapter at all and so have nothing
    but mRNA to hit. Mirrors ``_find_adapter_anchor_pos``, which is already fuzzy and
    already windowed on the 3' anchor.
    """
    pattern = SSP_FWD if frame == "fwd" else SSP_RC
    p = seq.find(pattern) if frame == "fwd" else seq.rfind(pattern)
    if p >= 0 or not HAS_EDLIB:
        return p
    if frame == "fwd":
        off, window = 0, seq[:SSP_SEARCH_WIN]
    else:
        off = max(0, len(seq) - SSP_SEARCH_WIN)
        window = seq[off:]
    r = edlib.align(pattern, window, mode="HW", task="locations", k=SSP_MAX_EDIT)
    if r["editDistance"] == -1 or not r["locations"]:
        return -1
    # edlib HW can return a location whose START is None (end found, start not
    # localizable). It carries no usable position — report the documented sentinel, the
    # same defect and fix as read_info._find_anchor_fuzzy and walkback._find_adapter_anchor_pos.
    start = r["locations"][-1][0] if frame == "rev" else r["locations"][0][0]
    return -1 if start is None else off + start


def _detect_frame(consensus_seq: str, orient: str, read_type: int) -> Tuple[str, bool]:
    """Return (frame_to_trim_in, disagreed_with_label).

    ``orient`` is a **BAM-SEQ-frame** label (``read_info.py:171`` calls it on
    ``read.query_sequence``; see the ``ReadInfo.orient`` docstring), but some consensus
    branches emit a basecalled-frame sequence. Rather than trust the label, key on the
    pattern actually present: exactly one of SSP_FWD / SSP_RC decides the frame. When both
    or neither is present (chimera, or an ONT-degraded SSP) fall back to the label.

    Type-2 molecules have no SSP by design, so the label is all there is — which is exactly
    why the caller must still propagate a correct frame (planning/681 CP0).
    """
    if read_type != 1:
        return orient, False
    has_fwd = SSP_FWD in consensus_seq
    has_rev = SSP_RC in consensus_seq
    if has_fwd == has_rev:            # both, or neither → no evidence; trust the label
        return orient, False
    detected = "fwd" if has_fwd else "rev"
    return detected, detected != orient


def pretrim_consensus(consensus_seq: str, orient: str, read_type: int) -> PretrimResult:
    """Strip SSP/UMI/GGG from the mRNA 5' end and poly-A/adapter from the 3' end.

    This ensures aligners receive clean mRNA sequence — no adapter contamination at either
    end. Downstream `rectify align` re-aligns the trimmed consensus and the aligner
    produces appropriate soft-clips for the stripped bases.

    `orient` is a **hint**, not a selector: the frame is taken from whichever SSP pattern
    the sequence actually carries (:func:`_detect_frame`), because the label is computed in
    BAM-SEQ frame while some consensus branches emit basecalled frame. See planning/681.

    5' trim (Type-1 only — Type-2 has no SSP/UMI):
      - fwd: SSP_FWD (23nt) + UMI (27nt) + GGG bridge (3nt) = 53 nt off the LEFT
      - rev: the mRNA 5' end is at the RIGHT — strip the CCC+UMI_RC+SSP_RC suffix.

    3' trim (both types):
      - fwd: poly-A at the RIGHT. Strip from the start of the first homopolymer A run of
        length >= MIN_HOMOPOLYMER_ANCHORED.
      - rev: poly-T at the LEFT. Strip position 0 through the end of the first T run.

    Returns a :class:`PretrimResult`.
    """
    frame, frame_mismatch = _detect_frame(consensus_seq, orient, read_type)
    n = len(consensus_seq)
    left_trim = right_trim = pa_len = 0

    if frame == "fwd":
        # ---- 5' (LEFT): SSP + UMI + GGG ----
        if read_type == 1:
            p = _find_ssp(consensus_seq, "fwd")
            if p >= 0:
                left_trim = p + len(SSP_FWD) + UMI_LEN + BRIDGE_LEN

        # ---- 3' (RIGHT): poly-A. Anchored (Tier 2) first, unanchored (Tier 1) second ----
        seq_to_search = consensus_seq[left_trim:]
        pa_start = None
        adp_pos = _find_adapter_anchor_pos(seq_to_search, "fwd")
        if adp_pos is not None:
            # Poly-A is in the 30 bp UPSTREAM of the adapter anchor
            upstream_start = max(0, adp_pos - ANCHOR_UPSTREAM_WIN)
            m = POLY_A_ANCH_RE.search(seq_to_search[upstream_start:adp_pos])
            if m:
                pa_start = left_trim + upstream_start + m.start()
                pa_len = m.group(0).count("A")
        if pa_start is None:
            # Unanchored fallback: rightmost A-run (>= MIN_HOMOPOLYMER_ANCHORED) in the last
            # END_WINDOW_BP bases. The anchored threshold (6+) rather than the conservative
            # unanchored one (10+) because in a cDNA consensus we know the poly-A is at the
            # 3' end; RIGHTMOST avoids stripping on internal A-runs (e.g. gene-body A's).
            window = seq_to_search[-END_WINDOW_BP:] if len(seq_to_search) > END_WINDOW_BP else seq_to_search
            window_off = len(seq_to_search) - len(window)
            matches = list(POLY_A_ANCH_RE.finditer(window))
            if matches:
                m = matches[-1]
                pa_start = left_trim + window_off + m.start()
                pa_len = m.group(0).count("A")
        if pa_start is not None:
            right_trim = n - pa_start

    else:  # frame == "rev" — mRNA runs right-to-left in this sequence
        # ---- 3' (LEFT): poly-T (the RC of the basecalled poly-A) ----
        adp_pos = _find_adapter_anchor_pos(consensus_seq, "rev")
        if adp_pos is not None:
            # poly-T is in the 30 bp DOWNSTREAM of the anchor (adapter at LEFT, polyT after)
            ds_start = adp_pos + ANCHOR_LEN
            downstream = consensus_seq[ds_start:min(ds_start + ANCHOR_UPSTREAM_WIN, n)]
            m = POLY_T_ANCH_RE.search(downstream)
            if m:
                left_trim = ds_start + m.end()
                pa_len = m.group(0).count("T")
        else:
            # Unanchored: first T-run in the first END_WINDOW_BP bases
            m = POLY_T_UNANCH_RE.search(consensus_seq[:END_WINDOW_BP])
            if m:
                left_trim = m.end()
                pa_len = m.group(0).count("T")

        # ---- 5' (RIGHT): CCC + UMI_RC + SSP_RC suffix ----
        if read_type == 1:
            rest = consensus_seq[left_trim:]
            p = _find_ssp(rest, "rev")
            if p >= UMI_LEN + BRIDGE_LEN:
                right_trim = len(rest) - (p - UMI_LEN - BRIDGE_LEN)

    # Degenerate guard: a spurious hit at each end can ask for more than the sequence has.
    # Trimming to empty makes the molecule unmappable, which is the failure mode
    # restore_eq_seq was written to prevent — drop the trim entirely instead.
    if left_trim < 0 or right_trim < 0 or left_trim + right_trim >= n:
        left_trim = right_trim = pa_len = 0

    trimmed = consensus_seq[left_trim:n - right_trim] if right_trim else consensus_seq[left_trim:]
    if frame == "fwd":
        trim_5p, trim_3p = left_trim, right_trim
    else:
        trim_5p, trim_3p = right_trim, left_trim
    return PretrimResult(trimmed, left_trim, right_trim, trim_5p, trim_3p,
                         pa_len, frame, frame_mismatch)
