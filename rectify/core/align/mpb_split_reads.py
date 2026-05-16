"""
mapPacBio split-read handling for RECTIFY.

BBTools mapPacBio.sh crashes on reads >6 kb.  This module provides:

  1. split_long_reads()  – split reads >MAX_MPB_READ_LENGTH into ≤6 kb chunks
  2. stitch_split_bam()  – recombine chunk alignments into full-read BAM records,
                           inserting N (intron) operations for reference gaps

Designed to be called from run_map_pacbio() in multi_aligner.py.
"""

import gzip
import logging
import re
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)

MAX_MPB_READ_LENGTH = 6000
CHUNK_TAG = "__mpbsplit_"

# ── pysam CIGAR operation codes ──────────────────────────────────────
_M, _I, _D, _N, _S, _H, _P, _EQ, _X = range(9)
_REF_OPS  = {_M, _D, _N, _EQ, _X}      # consume reference
_QUERY_OPS = {_M, _I, _S, _EQ, _X}      # consume query (incl. soft-clip)
_ALIGN_OPS = {_M, _I, _EQ, _X}          # consume query alignment (no S)


# ── 1.  Splitting ────────────────────────────────────────────────────

def split_long_reads(input_fastq, output_fastq, max_length=MAX_MPB_READ_LENGTH):
    """Write *output_fastq* with reads >max_length split into even chunks.

    Short reads are copied verbatim.

    Returns
    -------
    chunk_map : dict
        chunk_name → dict with keys:
            orig_name, chunk_idx, n_chunks, q_start, q_end, total_len
    n_split : int
        Number of reads that were split.
    """
    chunk_map = {}
    n_split = 0
    n_total = 0

    _open_in = gzip.open if str(input_fastq).endswith('.gz') else open
    with _open_in(input_fastq, 'rt') as fin, open(output_fastq, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            header = header.strip()
            if not header.startswith("@"):
                continue
            seq  = fin.readline().strip()
            plus = fin.readline().strip()
            qual = fin.readline().strip()
            n_total += 1

            read_name = header[1:].split()[0]

            # Normalize RNA notation to DNA: BBMap does not handle U bases and
            # throws a Java AssertionError on reads from RNA004/Dorado basecallers
            # that output U instead of T.
            seq = seq.replace('U', 'T').replace('u', 't')

            if len(seq) <= max_length:
                fout.write(f"{header}\n{seq}\n+\n{qual}\n")
                continue

            # -- split --
            n_split += 1
            n_chunks = -(-len(seq) // max_length)          # ceil division
            chunk_size = -(-len(seq) // n_chunks)           # even size

            for i in range(n_chunks):
                s = i * chunk_size
                e = min(s + chunk_size, len(seq))
                cname = f"{read_name}{CHUNK_TAG}{i}_{n_chunks}"
                fout.write(f"@{cname}\n{seq[s:e]}\n+\n{qual[s:e]}\n")
                chunk_map[cname] = dict(
                    orig_name=read_name, chunk_idx=i,
                    n_chunks=n_chunks, q_start=s, q_end=e,
                    total_len=len(seq),
                )

    if n_split:
        logger.info(
            "mapPacBio: split %d / %d reads (>%d bp) into chunks",
            n_split, n_total, max_length,
        )
    return chunk_map, n_split


# ── 2.  Stitching ────────────────────────────────────────────────────

def stitch_split_bam(input_bam, output_bam, chunk_map, threads=1):
    """Read *input_bam*, stitch chunk alignments, write *output_bam*.

    Non-chunk records are passed through unchanged.
    """
    import pysam

    chunk_groups = defaultdict(list)   # orig_name → [(info, AlignedSegment), …]
    regular = []

    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        header = bam_in.header.to_dict()

        for seg in bam_in:
            qname = seg.query_name
            if qname in chunk_map:
                chunk_groups[chunk_map[qname]["orig_name"]].append(
                    (chunk_map[qname], seg)
                )
            else:
                regular.append(seg)

    n_stitched = 0
    n_single   = 0

    with pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
        for seg in regular:
            bam_out.write(seg)

        for orig_name, chunks in chunk_groups.items():
            result = _stitch_group(orig_name, chunks, header)
            if result is not None:
                bam_out.write(result)
                if any(CHUNK_TAG in seg.query_name for _, seg in chunks if not seg.is_unmapped):
                    n_stitched += 1
                else:
                    n_single += 1

    logger.info(
        "mapPacBio stitch: %d stitched, %d single-chunk, %d passthrough",
        n_stitched, n_single, len(regular),
    )


# ── internal helpers ─────────────────────────────────────────────────

def _ref_len(cigar):
    """Reference bases consumed by pysam cigartuples list."""
    return sum(l for op, l in cigar if op in _REF_OPS)


def _query_aln_len(cigar):
    """Query bases consumed (alignment only, no soft-clip)."""
    return sum(l for op, l in cigar if op in _ALIGN_OPS)


def _strip_edge_clips(cigar):
    """Remove leading/trailing S operations, return (lead_S, core, trail_S)."""
    lead = 0
    trail = 0
    core = list(cigar)
    if core and core[0][0] == _S:
        lead = core[0][1]
        core = core[1:]
    if core and core[-1][0] == _S:
        trail = core[-1][1]
        core = core[:-1]
    return lead, core, trail


def _stitch_group(orig_name, chunks, header):
    """Stitch one group of chunk alignments → single AlignedSegment."""
    import pysam

    mapped = [(info, seg) for info, seg in chunks if not seg.is_unmapped]
    if not mapped:
        logger.debug("mapPacBio stitch: %s — all chunks unmapped; read dropped", orig_name)
        return None

    # Sort by chunk_idx (query order)
    mapped.sort(key=lambda x: x[0]["chunk_idx"])
    total_len = mapped[0][0]["total_len"]

    if len(mapped) == 1:
        return _promote_single(orig_name, mapped[0], total_len, header)

    # All on same chrom + strand?
    refs    = {seg.reference_name for _, seg in mapped}
    strands = {seg.is_reverse     for _, seg in mapped}

    if len(refs) > 1 or len(strands) > 1:
        # Fall back to longest single chunk
        best = max(mapped, key=lambda x: x[1].query_alignment_length)
        logger.debug(
            "mapPacBio stitch: %s chunks span multiple chromosomes/strands (%s) → using longest chunk",
            orig_name, refs,
        )
        return _promote_single(orig_name, best, total_len, header)

    is_reverse = mapped[0][1].is_reverse

    if is_reverse:
        # Reverse-strand stitching is complex (RC flips chunk order);
        # fall back to longest chunk for now.  Most DRS reads are sense-strand.
        best = max(mapped, key=lambda x: x[1].query_alignment_length)
        logger.debug(
            "mapPacBio stitch: %s reverse-strand multi-chunk → using longest chunk",
            orig_name,
        )
        return _promote_single(orig_name, best, total_len, header)

    # Forward strand: sort by reference start (== query order for + strand)
    by_ref = sorted(mapped, key=lambda x: x[1].reference_start)

    return _build_stitched(orig_name, by_ref, total_len, header)


def _promote_single(orig_name, item, total_len, header):
    """Wrap a single mapped chunk as a full-read alignment."""
    import pysam

    info, seg = item
    q_start = info["q_start"]
    q_end   = info["q_end"]

    cigar = list(seg.cigartuples) if seg.cigartuples else []
    chunk_seq  = seg.query_sequence or ""
    chunk_qual = list(seg.query_qualities) if seg.query_qualities is not None else None

    # Pad sequence / quality with placeholder for soft-clipped flanks
    before = q_start
    after  = total_len - q_end

    if seg.is_reverse:
        # In SAM convention the CIGAR reads left-to-right on the reference.
        # For a reverse-strand read the original query is reversed, so the
        # chunk that was at the *end* of the original read now appears first
        # on the reference.  Swap the padding sides accordingly.
        before, after = after, before

    # Prepend / append soft-clips to existing CIGAR
    if before > 0:
        if cigar and cigar[0][0] == _S:
            cigar[0] = (_S, cigar[0][1] + before)
        else:
            cigar.insert(0, (_S, before))
    if after > 0:
        if cigar and cigar[-1][0] == _S:
            cigar[-1] = (_S, cigar[-1][1] + after)
        else:
            cigar.append((_S, after))

    # Build padded SEQ / QUAL
    full_seq = "N" * before + chunk_seq + "N" * after
    if chunk_qual is not None:
        full_qual = [0] * before + chunk_qual + [0] * after
    else:
        full_qual = None

    new = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
    new.query_name      = orig_name
    new.flag             = seg.flag
    new.reference_id     = seg.reference_id
    new.reference_start  = seg.reference_start
    new.mapping_quality  = seg.mapping_quality
    new.cigartuples      = cigar
    new.query_sequence   = full_seq
    if full_qual is not None:
        new.query_qualities = full_qual

    # Carry over tags (NM, cs, etc.)
    for tag, val in seg.get_tags():
        try:
            new.set_tag(tag, val)
        except Exception:
            pass

    return new


def _build_stitched(orig_name, by_ref, total_len, header):
    """Stitch multiple forward-strand chunks sorted by reference position.

    Each chunk's CIGAR may have leading/trailing soft-clips.  The core
    (non-S) part is the aligned region.  Between consecutive chunks:

    - prev chunk's trailing soft-clip bases → I (insertion)
    - next chunk's leading soft-clip bases  → I (insertion)
    - reference gap between cores           → N (intron)

    Leading S of the first chunk and trailing S of the last chunk become
    the overall soft-clips of the stitched read.
    """
    import pysam

    # -- decompose each chunk's CIGAR --
    entries = []
    for info, seg in by_ref:
        cig = list(seg.cigartuples) if seg.cigartuples else []
        lead_s, core, trail_s = _strip_edge_clips(cig)
        core_ref = _ref_len(core)
        entries.append({
            "info":      info,
            "seg":       seg,
            "lead_s":    lead_s,
            "core":      core,
            "trail_s":   trail_s,
            "ref_start": seg.reference_start,
            "core_ref":  core_ref,
        })

    # ── build combined CIGAR ──
    combined = []

    # Overall leading soft-clip = bases before first chunk + first chunk's lead_S
    lead_clip = entries[0]["info"]["q_start"] + entries[0]["lead_s"]
    if lead_clip > 0:
        combined.append((_S, lead_clip))

    # First chunk core
    combined.extend(entries[0]["core"])

    prev_core_ref_end = entries[0]["ref_start"] + entries[0]["core_ref"]

    for i in range(1, len(entries)):
        prev = entries[i - 1]
        curr = entries[i]

        # Query bases between aligned cores of adjacent chunks:
        # = trailing S of prev chunk + leading S of curr chunk
        internal_q_gap = prev["trail_s"] + curr["lead_s"]

        # Reference gap between aligned cores
        ref_gap = curr["ref_start"] - prev_core_ref_end

        # Emit the internal query gap as insertion
        if internal_q_gap > 0:
            combined.append((_I, internal_q_gap))

        # Emit the reference gap as intron
        if ref_gap > 0:
            combined.append((_N, ref_gap))
        elif ref_gap < 0:
            logger.debug(
                "mapPacBio stitch: %s has %d bp ref overlap between chunks",
                orig_name, -ref_gap,
            )

        # This chunk's core alignment
        combined.extend(curr["core"])
        prev_core_ref_end = curr["ref_start"] + curr["core_ref"]

    # Overall trailing soft-clip = last chunk's trail_S + bases after last chunk
    trail_clip = entries[-1]["trail_s"] + (total_len - entries[-1]["info"]["q_end"])
    if trail_clip > 0:
        combined.append((_S, trail_clip))

    merged = _merge_cigar(combined)

    # ── reconstruct full-read sequence (chunk order = query order for + strand) ──
    idx_to_seg = {e["info"]["chunk_idx"]: e["seg"] for e in entries}
    all_indices = sorted(idx_to_seg.keys())
    parts_seq  = [idx_to_seg[i].query_sequence or "" for i in all_indices]
    parts_qual = [
        list(idx_to_seg[i].query_qualities)
        if idx_to_seg[i].query_qualities is not None else None
        for i in all_indices
    ]

    # Chunks cover [q_start_first .. q_end_last]; pad with N for overall clips
    pre_pad  = entries[0]["info"]["q_start"]   # normally 0
    post_pad = total_len - entries[-1]["info"]["q_end"]  # normally 0
    chunk_seq = "".join(parts_seq)
    full_seq  = "N" * pre_pad + chunk_seq + "N" * post_pad

    if all(q is not None for q in parts_qual):
        flat_qual = []
        for q in parts_qual:
            flat_qual.extend(q)
        full_qual = [0] * pre_pad + flat_qual + [0] * post_pad
    else:
        full_qual = None

    # ── sanity check: CIGAR query length must equal SEQ length ──
    cigar_q = sum(l for op, l in merged if op in _QUERY_OPS)
    if cigar_q != len(full_seq):
        logger.warning(
            "mapPacBio stitch: %s CIGAR query len %d != SEQ len %d — falling back",
            orig_name, cigar_q, len(full_seq),
        )
        best = max(by_ref, key=lambda x: x[1].query_alignment_length)
        return _promote_single(
            orig_name, (best[0], best[1]),  # (info, seg) from by_ref
            total_len, header,
        )

    # ── build output segment ──
    template = by_ref[0][1]
    new = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
    new.query_name      = orig_name
    new.flag             = template.flag
    new.reference_id     = template.reference_id
    new.reference_start  = entries[0]["ref_start"]
    new.mapping_quality  = min(e["seg"].mapping_quality for e in entries)
    new.cigartuples      = merged
    new.query_sequence   = full_seq
    if full_qual is not None:
        new.query_qualities = full_qual

    # Copy tags from template
    for tag, val in template.get_tags():
        try:
            new.set_tag(tag, val)
        except Exception:
            pass

    return new


def _merge_cigar(ops):
    """Merge runs of the same CIGAR operation and drop zero-length ops."""
    if not ops:
        return ops
    # Filter zero-length first so adjacent same-ops can merge
    filtered = [(op, l) for op, l in ops if l > 0]
    if not filtered:
        return []
    merged = [filtered[0]]
    for op, length in filtered[1:]:
        if op == merged[-1][0]:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return merged
