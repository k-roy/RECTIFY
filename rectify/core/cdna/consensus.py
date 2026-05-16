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
from typing import Dict, List, Optional, Tuple

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
    SSP_RC,
    UMI_LEN,
)
from .walkback import _find_adapter_anchor_pos

# Optional dependency (loaded lazily so the module imports without it).
try:
    import pyabpoa
    HAS_POA = True
except ImportError:
    HAS_POA = False


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
