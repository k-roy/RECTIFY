"""rectify correct-cdna — UMI-aware per-molecule consensus for PCR-cDNA.

Implementation of the `rectify correct-cdna` subcommand. Operates on a BAM
of aligned PCR-cDNA reads (typically SQK-PCB114.24 chemistry); clusters reads
by (locus, orientation, UMI) and emits one consensus record per starting
RNA molecule. Complements the FASTQ-level `rectify trim-cdna-polya` subcommand
(see ``cdna_trim_command.py``) which runs *before* alignment.

Implementation notes
--------------------
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

See ``docs/algorithms/cdna_correct.md`` (ported from cdna_dev DESIGN.md) for
the full algorithmic description and empirical validation of the PCB114.24
UMI architecture.

Ported from ``cdna_dev/src/cdna_correct.py`` v1.11 (M1 dev branch) on
2026-05-11; argparse setup moved to ``rectify/cli.py`` to integrate with
the rectify CLI surface.
"""

from __future__ import annotations

import argparse
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
    anchor: int
    orient: str          # 'fwd' (SSP at LEFT of BAM SEQ) or 'rev' (SSP_RC at RIGHT)
    umi: str             # 27-nt UMI, basecalled-orient
    is_reverse: bool
    xf_tier: int         # XF: 0=not detected, 1=unanchored, 2=anchored adapter+polyA


def extract_read_info(read: pysam.AlignedSegment) -> Optional[ReadInfo]:
    """Extract UMI + anchor + orient from a primary alignment. Returns None if no SSP found.

    ANCHOR LOGIC (fixed 2026-05-10): basecalled-3' anchor = aln_end-1 if not is_reverse
    else aln_start. This is independent of orient — for any read, basecalled-3' aligns
    to BAM SEQ's RIGHT end when FLAG=0 (BAM SEQ matches basecalled) and to BAM SEQ's
    LEFT end (= aln_start in + ref coords) when FLAG=16 (BAM SEQ is RC of basecalled).
    Earlier versions used `(orient=='fwd') != is_reverse` which inverted the anchor
    for - strand gene reads and produced flipped sense/antisense labels.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return None
    seq = read.query_sequence
    if seq is None:
        return None

    p = seq.find(SSP_FWD)
    if p >= 0:
        umi_basecalled = seq[p + len(SSP_FWD): p + len(SSP_FWD) + UMI_LEN]
        if len(umi_basecalled) != UMI_LEN:
            return None
        orient = "fwd"
    else:
        p = seq.find(SSP_RC)
        if p < UMI_LEN:
            return None
        umi_rc = seq[p - UMI_LEN: p]
        umi_basecalled = revcomp(umi_rc)
        if len(umi_basecalled) != UMI_LEN:
            return None
        orient = "rev"

    # basecalled-3' anchor: where the cDNA's basecalled 3' end (= the reliable
    # translocation end, since 3' enters the pore first) sits in genomic coords.
    if not read.is_reverse:
        anchor = (read.reference_end or 0) - 1
    else:
        anchor = read.reference_start

    return ReadInfo(
        read_id=read.query_name,
        chrom=read.reference_name,
        anchor=anchor,
        orient=orient,
        umi=umi_basecalled,
        is_reverse=read.is_reverse,
        xf_tier=detect_full_length_tier(seq, orient),
    )


# ---- Clustering ------------------------------------------------------------
def umi_components(umis: List[str], max_edit: int) -> List[List[int]]:
    """Connected-components clustering of UMIs by Levenshtein ≤ max_edit."""
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


def cluster_reads(reads: List[ReadInfo], anchor_window: int, max_edit: int,
                  per_cluster_cap: int) -> Tuple[List[List[ReadInfo]], dict]:
    anchor_buckets: dict = defaultdict(list)
    for r in reads:
        bucket = (r.chrom, r.anchor // anchor_window, r.orient)
        anchor_buckets[bucket].append(r)

    stats = dict(
        anchor_buckets=len(anchor_buckets),
        buckets_dropped_polyA_pileup=0,
        reads_in_dropped_buckets=0,
        molecule_clusters=0,
        molecule_clusters_size_1=0,
        molecule_clusters_size_2=0,
        molecule_clusters_size_ge_3=0,
    )

    out_clusters: List[List[ReadInfo]] = []
    for bucket_reads in anchor_buckets.values():
        if len(bucket_reads) > per_cluster_cap:
            stats["buckets_dropped_polyA_pileup"] += 1
            stats["reads_in_dropped_buckets"] += len(bucket_reads)
            continue
        umis = [r.umi for r in bucket_reads]
        for comp in umi_components(umis, max_edit):
            cluster = [bucket_reads[i] for i in comp]
            out_clusters.append(cluster)
            sz = len(cluster)
            if sz == 1: stats["molecule_clusters_size_1"] += 1
            elif sz == 2: stats["molecule_clusters_size_2"] += 1
            else: stats["molecule_clusters_size_ge_3"] += 1
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


def poa_consensus(reads_in_cluster: List[pysam.AlignedSegment]) -> Optional[str]:
    """Build POA consensus sequence across cluster reads (v2 upgrade over pileup).

    Operates on the BAM SEQ (in reference orientation) of each read. POA handles
    indels in homopolymer regions much better than pileup. Returns the consensus
    as a string, or None if abPOA isn't available or the cluster has too few reads.

    The consensus needs to be re-aligned to the reference afterward (use mappy)
    to obtain a valid CIGAR — the POA itself doesn't produce reference alignment.
    """
    if not HAS_POA or len(reads_in_cluster) < 2:
        return None
    seqs = [r.query_sequence for r in reads_in_cluster if r.query_sequence]
    if len(seqs) < 2:
        return None
    aligner = pyabpoa.msa_aligner()
    res = aligner.msa(seqs, out_cons=True, out_msa=False)
    if not res.cons_seq:
        return None
    return res.cons_seq[0]


def realign_consensus(consensus_seq: str, mp_aligner: "mappy.Aligner") -> Optional[Tuple[str, int, List[Tuple[int, int]], int, int, str]]:
    """Re-align a consensus sequence to the reference via mappy/minimap2.

    Returns (ref_name, ref_start, cigartuples, mapq, flag, query_seq_for_record)
    of the best hit, or None.

    The returned cigartuples include soft-clips (op=4) for any query bases
    outside the aligned region. The returned query_seq_for_record matches the
    CIGAR's expected query length: if mappy aligned to the - strand, the
    consensus is reverse-complemented (so BAM SEQ matches reference orient).
    """
    if not HAS_MAPPY:
        return None
    best = None
    for hit in mp_aligner.map(consensus_seq):
        if hit.is_primary:
            best = hit
            break
    if best is None:
        return None

    # mappy hit.cigar: list of (length, op_int) — BAM op codes
    aln_cigar = [(op, length) for length, op in best.cigar]

    # Soft-clip prefix and suffix to make CIGAR sum equal full query length
    pre_clip = best.q_st
    post_clip = len(consensus_seq) - best.q_en
    cigartuples: List[Tuple[int, int]] = []
    if pre_clip > 0:
        cigartuples.append((4, pre_clip))
    cigartuples.extend(aln_cigar)
    if post_clip > 0:
        cigartuples.append((4, post_clip))

    flag = 16 if best.strand == -1 else 0
    # If aligned to - strand, BAM convention: SEQ is RC'd to match + ref
    if best.strand == -1:
        seq_for_record = consensus_seq[::-1].translate(_COMPL)
        # Soft-clip positions also flip
        cigartuples = []
        if post_clip > 0:
            cigartuples.append((4, post_clip))
        cigartuples.extend(reversed(aln_cigar))
        if pre_clip > 0:
            cigartuples.append((4, pre_clip))
    else:
        seq_for_record = consensus_seq

    return (best.ctg, best.r_st, cigartuples, best.mapq, flag, seq_for_record)


def make_consensus_record(template: pysam.AlignedSegment,
                           consensus_seq: str,
                           cigartuples: List[Tuple[int, int]],
                           ref_start: int,
                           qual_string: str,
                           umi: str,
                           orient: str,
                           cluster_size: int,
                           read_ids: List[str],
                           header: pysam.AlignmentHeader) -> pysam.AlignedSegment:
    """Build a synthetic consensus AlignedSegment from a multi-read cluster.

    The template provides reference_id and FLAG (orientation); everything
    sequence-related comes from the consensus.
    """
    rec = pysam.AlignedSegment(header)
    rec.query_name = f"consensus_{umi}_{orient}_{template.reference_name}_{ref_start}_{cluster_size}reads"
    rec.flag = template.flag & ~0x900  # primary; clear secondary/supplementary
    rec.reference_id = template.reference_id
    rec.reference_start = ref_start
    rec.mapping_quality = 60
    rec.cigartuples = cigartuples
    rec.query_sequence = consensus_seq
    rec.query_qualities = pysam.qualitystring_to_array(qual_string)
    rec.next_reference_id = -1
    rec.next_reference_start = -1
    rec.template_length = 0
    return rec


def write_stage1_bam(input_bam: Path, output_bam: Path,
                      clusters: List[List[ReadInfo]],
                      umi_canonical: Dict[int, str],
                      cluster_xs: Dict[int, str],
                      cluster_xf_tier: Dict[int, int],
                      paired_partner: Dict[int, str],
                      mp_aligner: Optional["mappy.Aligner"] = None,
                      use_poa: bool = False) -> dict:
    """Re-stream the input BAM and emit one record per cluster.

    For singletons: pass-through with added tags.
    For multi-read: pick representative read, copy with added tags.

    Tags written:
      XU:Z:<umi_canonical>           — canonical UMI of the cluster
      XO:Z:<orient>                  — 'fwd' or 'rev'
      XC:i:<cluster_size>            — number of input reads in this cluster
      XR:Z:<all_read_ids>            — comma-separated input read IDs
    """
    # Build read_id -> cluster_id index
    rid_to_cluster: Dict[str, int] = {}
    cluster_read_ids: Dict[int, List[str]] = {}
    for cid, c in enumerate(clusters):
        cluster_read_ids[cid] = [r.read_id for r in c]
        for r in c:
            rid_to_cluster[r.read_id] = cid

    # Bucket: per cluster, collect AlignedSegment as we stream the BAM
    cluster_segments: Dict[int, List[pysam.AlignedSegment]] = defaultdict(list)

    n_in = 0
    with pysam.AlignmentFile(str(input_bam), "rb") as inbam:
        header = inbam.header.to_dict()
        # Add @PG line
        pg = header.setdefault("PG", [])
        pg.append({"ID": "rectify-correct-cdna",
                    "PN": "rectify",
                    "VN": "0.1.0-v1",
                    "CL": " ".join(sys.argv)})
        for read in inbam.fetch(until_eof=True):
            n_in += 1
            cid = rid_to_cluster.get(read.query_name)
            if cid is None:
                continue
            cluster_segments[cid].append(read)

    written = singleton = multi_pileup = multi_fallback = 0
    bam_header = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(str(output_bam), "wb", header=bam_header) as out:
        for cid, c in enumerate(clusters):
            segs = cluster_segments.get(cid, [])
            if not segs:
                continue
            if len(segs) == 1:
                # Singleton — pass-through with tags
                rec = segs[0]
                rec.set_tag("XU", umi_canonical[cid], value_type="Z")
                rec.set_tag("XO", c[0].orient, value_type="Z")
                rec.set_tag("XC", len(c), value_type="i")
                rec.set_tag("XR", ",".join(cluster_read_ids[cid]), value_type="Z")
                rec.set_tag("XM", "rep", value_type="Z")  # consensus method
                rec.set_tag("XS", cluster_xs.get(cid, "unannotated"), value_type="Z")
                if cid in paired_partner:
                    rec.set_tag("XP", paired_partner[cid], value_type="Z")
                rec.set_tag("XF", cluster_xf_tier[cid], value_type="i")
                out.write(rec)
                singleton += 1
            else:
                rec = None
                # Try POA + re-align first (best quality)
                if use_poa and mp_aligner is not None:
                    poa_seq = poa_consensus(segs)
                    if poa_seq:
                        realigned = realign_consensus(poa_seq, mp_aligner)
                        if realigned:
                            ref_name, ref_start, cigartuples, mapq, new_flag, rec_seq = realigned
                            template = pick_representative(segs)
                            rec = make_consensus_record(template, rec_seq, cigartuples,
                                                         ref_start, "?" * len(rec_seq),
                                                         umi_canonical[cid], c[0].orient,
                                                         len(c), cluster_read_ids[cid], bam_header)
                            rec.flag = new_flag  # set strand from re-alignment
                            rec.mapping_quality = mapq
                            rec.set_tag("XM", "poa", value_type="Z")
                            multi_pileup += 1  # re-using counter; v2 will split
                # Fallback: pileup-based consensus
                if rec is None:
                    pileup = pileup_consensus(segs)
                    if pileup is None:
                        rec = pick_representative(segs)
                        rec.set_tag("XM", "rep_fallback", value_type="Z")
                        multi_fallback += 1
                    else:
                        consensus_seq, cigartuples, ref_start, qual_str = pileup
                        template = pick_representative(segs)
                        rec = make_consensus_record(template, consensus_seq, cigartuples,
                                                     ref_start, qual_str,
                                                     umi_canonical[cid], c[0].orient,
                                                     len(c), cluster_read_ids[cid], bam_header)
                        rec.set_tag("XM", "pileup", value_type="Z")
                        multi_pileup += 1
                rec.set_tag("XU", umi_canonical[cid], value_type="Z")
                rec.set_tag("XO", c[0].orient, value_type="Z")
                rec.set_tag("XC", len(c), value_type="i")
                rec.set_tag("XR", ",".join(cluster_read_ids[cid]), value_type="Z")
                rec.set_tag("XS", cluster_xs.get(cid, "unannotated"), value_type="Z")
                if cid in paired_partner:
                    rec.set_tag("XP", paired_partner[cid], value_type="Z")
                rec.set_tag("XF", cluster_xf_tier[cid], value_type="i")
                out.write(rec)
            written += 1

    return dict(input_reads=n_in, written=written, from_singletons=singleton,
                from_multi_pileup=multi_pileup, from_multi_fallback=multi_fallback)


def canonical_umi(cluster_umis: List[str]) -> str:
    """Pick the canonical UMI for a cluster (most frequent; ties → first lexicographically)."""
    counts = defaultdict(int)
    for u in cluster_umis:
        counts[u] += 1
    return min(counts, key=lambda u: (-counts[u], u))


# ---- v1.6: gene-strand–aware sense/antisense classification ----------------
def load_gff_genes(gff_path: Path) -> Dict[str, IntervalTree]:
    """Load gene-like features from a GFF file → chrom → IntervalTree(strand)."""
    trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    GENE_TYPES = {"gene", "mRNA", "transcript", "pseudogene", "ncRNA_gene",
                   "tRNA_gene", "snoRNA_gene", "snRNA_gene", "rRNA_gene"}
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            if parts[2] not in GENE_TYPES:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            if start < end and strand in ("+", "-"):
                trees[chrom][start:end] = strand
    return dict(trees)


def classify_sense_antisense(chrom: str, anchor: int,
                              gene_trees: Dict[str, IntervalTree],
                              terminal_window: int = 200) -> str:
    """Classify a cluster's basecalled-3' anchor as 'sense' / 'antisense' / 'unannotated'.

    Logic: For PCB114 cDNA, the basecalled-3' anchor lies at the polyA end of the
    cDNA molecule. If that anchor is closer to the gene's polyA-side end (3' end
    of mRNA: gene_end on + strand, gene_start on - strand), the cDNA matched the
    mRNA → SENSE. If closer to the cap-side end, the cDNA was RC of mRNA → ANTISENSE.

    `terminal_window` (default 200 bp) extends gene boundaries to catch reads
    just outside the annotation — common at 3' UTRs not fully annotated.
    """
    tree = gene_trees.get(chrom)
    if tree is None:
        return "unannotated"
    # Direct overlap first
    overlaps = list(tree[anchor])
    if not overlaps:
        # Allow small terminal window
        overlaps = list(tree.overlap(anchor - terminal_window, anchor + terminal_window))
    if not overlaps:
        return "unannotated"
    # Pick the gene whose polyA end is closest to the anchor
    best_call = None
    best_d = float("inf")
    for iv in overlaps:
        if iv.data == "+":
            polyA_pos = iv.end
            cap_pos = iv.begin
        else:
            polyA_pos = iv.begin
            cap_pos = iv.end
        d_polyA = abs(anchor - polyA_pos)
        d_cap = abs(anchor - cap_pos)
        if d_polyA < d_cap:
            call = "sense"
            d = d_polyA
        elif d_cap < d_polyA:
            call = "antisense"
            d = d_cap
        else:
            call = "ambiguous"
            d = d_polyA
        if d < best_d:
            best_d = d
            best_call = call
    return best_call or "unannotated"


# ---- I/O -------------------------------------------------------------------
def stream_reads(bam_path: Path, region: Optional[str]) -> Iterator[ReadInfo]:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        it = bam.fetch(region=region) if region else bam.fetch(until_eof=True)
        for read in it:
            info = extract_read_info(read)
            if info is not None:
                yield info


# ---- Entry point -----------------------------------------------------------
def run(args) -> int:
    """Entry point called from rectify/cli.py for ``rectify correct-cdna``.

    Expects the args namespace populated by the rectify subparser; argument
    registration lives in cli.py to keep the rectify CLI surface consistent.
    """
    # Normalise args: cli.py passes strings/None for path-like options.
    args.bam = Path(args.bam)
    args.out = Path(args.out)
    args.gff = Path(args.gff) if getattr(args, 'gff', None) else None
    args.reference = Path(args.reference) if getattr(args, 'reference', None) else None

    logging.basicConfig(
        level=logging.DEBUG if getattr(args, 'verbose', False) else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s")
    log = logging.getLogger("cdna_correct")

    if not args.bam.exists():
        log.error("BAM not found: %s", args.bam)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    log.info("Streaming reads from %s region=%s ...", args.bam, args.region or "all")
    reads = list(stream_reads(args.bam, args.region))
    log.info("  %d reads with extractable UMIs (%.1fs)", len(reads), time.time() - t0)

    t1 = time.time()
    log.info("Stage 1: clustering (anchor window %d, UMI Lev≤%d, cap %d)",
             args.anchor_window_bp, args.umi_edit_distance, args.per_cluster_cap)
    clusters, stats = cluster_reads(reads, args.anchor_window_bp,
                                     args.umi_edit_distance, args.per_cluster_cap)
    log.info("  %d molecule clusters (%.1fs)", stats["molecule_clusters"], time.time() - t1)

    # Canonical UMI per cluster
    umi_canon = {cid: canonical_umi([r.umi for r in c]) for cid, c in enumerate(clusters)}

    # Per-cluster XF tier: max of constituent reads (any one full-length is enough)
    cluster_xf_tier = {cid: max((r.xf_tier for r in c), default=0) for cid, c in enumerate(clusters)}

    # GFF-based sense/antisense classification
    cluster_xs: Dict[int, str] = {}
    if args.gff:
        log.info("Loading GFF: %s", args.gff)
        gene_trees = load_gff_genes(args.gff)
        log.info("  loaded gene trees for %d chromosomes (%d total intervals)",
                 len(gene_trees), sum(len(t) for t in gene_trees.values()))
        for cid, c in enumerate(clusters):
            cluster_xs[cid] = classify_sense_antisense(c[0].chrom, c[0].anchor, gene_trees)
    else:
        for cid in range(len(clusters)):
            cluster_xs[cid] = "unannotated"

    # Stage 2 cross-orient pair check (compute first, before BAM write, so we can override XS=paired)
    fwd = [(cid, umi_canon[cid], c[0].chrom, c[0].anchor, len(c))
            for cid, c in enumerate(clusters) if c[0].orient == "fwd"]
    rev = [(cid, umi_canon[cid], c[0].chrom, c[0].anchor, len(c))
            for cid, c in enumerate(clusters) if c[0].orient == "rev"]
    paired_pairs = []
    fwd_by_chrom = defaultdict(list)
    for e in fwd:
        fwd_by_chrom[e[2]].append(e)
    for r_cid, r_umi, r_ch, r_a, r_sz in rev:
        for f_cid, f_umi, f_ch, f_a, f_sz in fwd_by_chrom.get(r_ch, []):
            if abs(r_a - f_a) > args.max_cross_orient_span:
                continue
            d = Levenshtein.distance(f_umi, r_umi, score_cutoff=args.umi_edit_distance)
            if d <= args.umi_edit_distance:
                paired_pairs.append((f_cid, r_cid, f_ch, f_a, r_a, f_sz, r_sz, d, abs(r_a - f_a)))

    # Override XS for paired clusters; build paired_partner mapping
    paired_partner: Dict[int, str] = {}
    for f_cid, r_cid, *_ in paired_pairs:
        # Use the canonical UMI of the fwd cluster as the molecule ID
        molecule_id = umi_canon[f_cid]
        paired_partner[f_cid] = molecule_id
        paired_partner[r_cid] = molecule_id
        cluster_xs[f_cid] = "paired"
        cluster_xs[r_cid] = "paired"

    # Stage 1 BAM (write unsorted, then sort+index)
    t2 = time.time()
    unsorted_bam = args.out / "stage1_consensus.unsorted.bam"
    out_bam = args.out / "stage1_consensus.bam"
    log.info("Writing Stage-1 consensus BAM (representative-read mode) → %s", out_bam)
    # Set up POA aligner if requested
    use_poa = HAS_POA and not args.no_poa and args.reference is not None
    mp_aligner = None
    if use_poa:
        if not HAS_MAPPY:
            log.warning("mappy not available; falling back to pileup consensus")
            use_poa = False
        else:
            log.info("Building mappy aligner from %s ...", args.reference)
            mp_aligner = mappy.Aligner(str(args.reference), preset="splice")
            if not mp_aligner:
                log.warning("Failed to build mappy aligner; falling back to pileup")
                use_poa = False
    if use_poa:
        log.info("Multi-read cluster consensus: POA (pyabpoa) + re-align (mappy/minimap2)")
    else:
        log.info("Multi-read cluster consensus: pileup-style (substitution-corrective only)")

    bam_stats = write_stage1_bam(args.bam, unsorted_bam, clusters, umi_canon,
                                   cluster_xs, cluster_xf_tier, paired_partner,
                                   mp_aligner=mp_aligner, use_poa=use_poa)
    log.info("  wrote %d records (%d singletons, %d pileup consensus, %d rep fallback) in %.1fs",
             bam_stats["written"], bam_stats["from_singletons"],
             bam_stats["from_multi_pileup"], bam_stats["from_multi_fallback"],
             time.time() - t2)
    pysam.sort("-o", str(out_bam), str(unsorted_bam))
    unsorted_bam.unlink()
    pysam.index(str(out_bam))
    log.info("  sorted + indexed → %s(.bai)", out_bam)

    # Cluster manifest with XS classification
    manifest = args.out / "clusters.tsv"
    with manifest.open("w") as f:
        f.write("cluster_id\tchrom\tanchor\torient\tn_reads\tumi_canonical\txs\tread_ids\n")
        for cid, c in enumerate(clusters):
            f.write(f"{cid}\t{c[0].chrom}\t{c[0].anchor}\t{c[0].orient}"
                    f"\t{len(c)}\t{umi_canon[cid]}\t{cluster_xs[cid]}"
                    f"\t{','.join(r.read_id for r in c)}\n")
    log.info("Wrote cluster manifest → %s", manifest)

    pairs_tsv = args.out / "stage2_pairs.tsv"
    with pairs_tsv.open("w") as f:
        f.write("fwd_cid\trev_cid\tchrom\tfwd_anchor\trev_anchor\tfwd_size\trev_size\tumi_levdist\tspan_bp\n")
        for p in paired_pairs:
            f.write("\t".join(str(x) for x in p) + "\n")
    log.info("Stage 2: %d cross-orient pairs → %s", len(paired_pairs), pairs_tsv)

    # Final report
    print()
    print("=" * 70)
    print("cdna_correct v1 — Stage-1 dedup complete")
    print("=" * 70)
    print(f"input reads (UMI-extractable): {len(reads):>8d}")
    print(f"output records (one per molecule): {bam_stats['written']:>8d}"
          f"  ({100 * bam_stats['written'] / max(1, len(reads)):.1f}% of input)")
    print(f"  singletons (passed through):     {bam_stats['from_singletons']:>8d}")
    print(f"  multi-read (pileup consensus):   {bam_stats['from_multi_pileup']:>8d}")
    print(f"  multi-read (rep fallback):       {bam_stats['from_multi_fallback']:>8d}")
    print(f"polyA-pileup buckets dropped:    {stats['buckets_dropped_polyA_pileup']:>8d}"
          f"  ({stats['reads_in_dropped_buckets']} reads — these need POA + position-aware handling)")
    print(f"Stage 2 cross-orient pairs (span≤{args.max_cross_orient_span}, Lev≤{args.umi_edit_distance}):"
          f" {len(paired_pairs):>5d}")
    print()
    print("XS classification breakdown (per cluster):")
    xs_counts = defaultdict(int)
    for v in cluster_xs.values():
        xs_counts[v] += 1
    for cat in ("sense", "antisense", "paired", "unannotated", "ambiguous"):
        if xs_counts.get(cat, 0):
            print(f"  XS={cat:<12} {xs_counts[cat]:>8d} clusters"
                  f"  ({100*xs_counts[cat]/len(cluster_xs):.1f}%)")
    print()
    print("XF tier breakdown (full-length confidence):")
    tier_counts = defaultdict(int)
    for v in cluster_xf_tier.values():
        tier_counts[v] += 1
    n_clust = len(clusters)
    print(f"  XF=2 (anchored, HIGH confidence):  {tier_counts[2]:>8d} clusters"
          f"  ({100*tier_counts[2]/n_clust:.1f}%)")
    print(f"  XF=1 (unanchored, MEDIUM confidence): {tier_counts[1]:>8d} clusters"
          f"  ({100*tier_counts[1]/n_clust:.1f}%)")
    print(f"  XF=0 (not detected):                {tier_counts[0]:>8d} clusters"
          f"  ({100*tier_counts[0]/n_clust:.1f}%)")
    n_full = tier_counts[1] + tier_counts[2]
    print(f"  XF≥1 (any full-length):              {n_full:>8d} clusters"
          f"  ({100*n_full/n_clust:.1f}%)")
    print()
    print("Cross-table XS × XF tier:")
    print(f"  {'XS':<14} {'XF=2':>8} {'XF=1':>8} {'XF=0':>8}  {'%XF≥1':>8}  {'%XF=2':>8}")
    xs_xf = defaultdict(lambda: defaultdict(int))
    for cid in range(len(clusters)):
        xs_xf[cluster_xs.get(cid, "unannotated")][cluster_xf_tier[cid]] += 1
    for cat in ("sense", "antisense", "paired", "unannotated", "ambiguous"):
        row = xs_xf.get(cat)
        if not row:
            continue
        n2 = row[2]; n1 = row[1]; n0 = row[0]
        tot = n2 + n1 + n0
        if tot == 0: continue
        pct_full = 100 * (n2 + n1) / tot
        pct_high = 100 * n2 / tot
        print(f"  {cat:<14} {n2:>8d} {n1:>8d} {n0:>8d}  {pct_full:>7.1f}%  {pct_high:>7.1f}%")
    print()
    print(f"Output BAM: {out_bam}")
    print(f"Cluster manifest: {manifest}")
    print(f"Stage-2 pairs: {pairs_tsv}")
    log.info("Total runtime: %.1fs", time.time() - t0)
    return 0
