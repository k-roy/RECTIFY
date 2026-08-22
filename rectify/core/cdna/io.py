"""BAM streaming + Stage-1 consensus FASTQ emission for the cdna pipeline.

  * :func:`stream_reads` — fetch records from a BAM, mask rRNA_gene overlaps,
    extract :class:`ReadInfo` for every UMI-extractable primary alignment.
  * :func:`write_stage1_fastq` — for each cluster, build a consensus
    (POA → pileup → representative-read fallback), pretrim adapters, and emit
    a FASTQ record with the alignment-independent SAM tags packed into the
    comment line for ``rectify align -y`` to pass through.
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam

from ._constants import COMPLEMENT_TABLE
from .cluster import pick_representative
from .consensus import (
    pileup_consensus,
    poa_consensus,
    poa_consensus_strand_aware,
    pretrim_consensus,
    restore_eq_seq,
)
from .read_info import ReadInfo, extract_read_info

# Reverse-complementing a consensus swaps which end the SSP sits at, which is exactly the
# fwd/rev distinction (`ReadInfo.orient`: "fwd" = SSP at LEFT of BAM SEQ, "rev" = SSP_RC at
# RIGHT). So flipping the label alongside the sequence is definitionally exact.
_FLIP_ORIENT = {"fwd": "rev", "rev": "fwd"}


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


def write_stage1_fastq(input_bam: Path, output_fastq: Path,
                       clusters: List[List[ReadInfo]],
                       umi_canonical: Dict[int, str],
                       cluster_xf_tier: Dict[int, int],
                       cluster_tail_len: Dict[int, int],
                       use_poa: bool = False,
                       strand_aware_consensus: bool = False,
                       cluster_name_prefix: str = "cluster",
                       reference: Optional[Path] = None) -> dict:
    """Emit per-cluster consensus sequences as a FASTQ for downstream alignment.

    `rectify align` will run the multi-aligner on this FASTQ to produce the final
    aligned BAM. Per-cluster SAM-format tags are appended to each read's comment
    line so `minimap2 -y` (and equivalent flags on mapPacBio/gapmm2) can
    propagate them through to the output BAM.

    Only alignment-independent tags are emitted here. Gene/strand/isoform/pair
    tags (XG, XS, XI, XL) are recomputed by `rectify cdna-analyze` on the
    post-align CIGAR, so emitting them here would be misleading by the time
    downstream consumers read them. Tail length IS emitted as XA and is
    authoritative whenever the tail was stripped pre-alignment (b3a8c35+):
    cdna-analyze only overrides it with its own walkback when the tail is still
    present in the aligned SEQ.

    Read naming: `cluster_<cid>`.
    """
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
    n_frame_flipped = n_frame_mismatch = n_noop_5 = n_noop_3 = 0

    # Reference is needed to resolve calmd '=' (match-to-ref) placeholders in SEQ
    # back to real bases before re-emitting; without it the FASTQ is unmappable.
    fasta = pysam.FastaFile(str(reference)) if reference is not None else None

    with opener(str(output_fastq), "wt") as fq:
        for cid, c in enumerate(clusters):
            segs = cluster_segments.get(cid, [])
            if not segs:
                continue

            method = "rep"
            seq: Optional[str] = None
            # Frame bookkeeping (planning/681). `orient` is a BAM-SEQ-frame label — it is
            # computed by extract_read_info on `read.query_sequence`, which pysam returns
            # reference-forward. The singleton and rep_fallback branches then RC the
            # sequence into basecalled frame, so for a minus-strand source read the label
            # and the sequence end up in OPPOSITE frames and every trim searches the wrong
            # pattern at the wrong end. That silently no-op'd the trim on 39.7% of all
            # molecules and sent ~145 nt of adapter+UMI+poly-A into the aligner as a soft
            # clip. Track which frame the assembled sequence is actually in and flip the
            # label with it; by the codebase's own definition ("fwd" = SSP_FWD at the LEFT
            # of BAM SEQ) the RC of a fwd sequence IS a rev sequence, so the flip is exact.
            seq_in_bam_frame = True

            if len(segs) == 1:
                seg = segs[0]
                seq = restore_eq_seq(seg, fasta) or ""
                # Restore basecalled orientation: BAM SEQ is RC'd for minus-strand
                # alignments; we want the original read sequence for re-alignment.
                if seg.is_reverse:
                    seq = seq[::-1].translate(COMPLEMENT_TABLE)
                    seq_in_bam_frame = False
                singleton += 1
            else:
                if use_poa:
                    poa_fn = poa_consensus_strand_aware if strand_aware_consensus else poa_consensus
                    seq = poa_fn(segs, fasta)
                    if seq:
                        method = "poa"
                if not seq:
                    pileup = pileup_consensus(segs, fasta)
                    if pileup is not None:
                        seq = pileup[0]
                        method = "pileup"
                if not seq:
                    rep = pick_representative(segs)
                    seq = restore_eq_seq(rep, fasta) or ""
                    # Key the flip on the segment that was actually RC'd. pick_representative
                    # need not return c[0], so `rep.is_reverse` and c[0]'s strand can differ.
                    if rep.is_reverse:
                        seq = seq[::-1].translate(COMPLEMENT_TABLE)
                        seq_in_bam_frame = False
                    method = "rep_fallback"
                    multi_fallback += 1
                else:
                    multi_pileup += 1
                # poa_consensus / pileup_consensus both build in BAM (reference) frame and
                # apply no RC, so seq_in_bam_frame stays True for them — which is CORRECT
                # and must stay that way. A pileup consensus is built from
                # get_aligned_pairs(matches_only=True) and so contains no soft-clip, hence
                # no adapter to strip at all (planning/681 CP2: its molecules carry 13.7 nt
                # clips vs 145.2 nt for the genuinely untrimmed singletons). Adding the RC
                # here would break a correctly-framed branch.

            if not seq:
                continue

            # Strip adapter/UMI/GGG prefix and poly-A/T suffix so rectify align
            # receives clean mRNA sequence. Store strip lengths as XQ/XK tags so
            # the aligner can reconstruct the full-length CIGAR with soft-clips.
            orient = c[0].orient
            read_type = c[0].read_type
            trim_orient = orient if seq_in_bam_frame else _FLIP_ORIENT[orient]
            if not seq_in_bam_frame:
                n_frame_flipped += 1
            pre = pretrim_consensus(seq, trim_orient, read_type)
            trimmed_seq, q_trim_5, q_trim_3 = pre.seq, pre.trim_5p, pre.trim_3p
            if not trimmed_seq:
                trimmed_seq = seq
                q_trim_5 = q_trim_3 = 0
            # Fail-loud instrumentation. A silent trim no-op is what let this defect run
            # for months as a throughput mystery; surface it as a number instead.
            if pre.frame_mismatch:
                n_frame_mismatch += 1
            if read_type == 1 and q_trim_5 == 0:
                n_noop_5 += 1
            if q_trim_3 == 0:
                n_noop_3 += 1

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
            fq.write(f"@{cluster_name_prefix}_{cid}\t{comment}\n{trimmed_seq}\n+\n{qual}\n")
            written += 1

    if fasta is not None:
        fasta.close()

    return dict(input_reads=n_in, written=written, from_singletons=singleton,
                from_multi_pileup=multi_pileup, from_multi_fallback=multi_fallback,
                trim_frame_flipped=n_frame_flipped,
                trim_frame_mismatch=n_frame_mismatch,
                trim_noop_5p=n_noop_5, trim_noop_3p=n_noop_3)
