"""CMA v1 reader — expand() reconstructs each aligner's pre-correct placement.

expand() rebuilds the per-aligner ``Dict[aligner -> record]`` **from the Za tags**
(planning/254 §3.1, HARD GUARD judge #5) — it MUST NOT re-feed the CMA through
consensus._iter_name_grouped_bams, whose per-key max(_n_ops) collapse would drop
variant placements.

SEQ/QUAL reconstruction frame (planning/254 §2.4 — the #1 silent-failure risk,
fails on reverse-strand reads only):

    same   = variant.is_reverse == payload.is_reverse
    base   = payload_SEQ  if same else revcomp(payload_SEQ)
    baseq  = payload_QUAL if same else payload_QUAL[::-1]     # reverse, not complement
    qs, qe = Zq (or leading_H / L-trailing_H from CIGAR)      # soft-clips ride inside
    SEQ, QUAL = base[qs:qe], baseq[qs:qe]
"""

from __future__ import annotations

import logging
from array import array
from typing import Dict, Iterator, List, Optional, Tuple

import pysam

from .cma_schema import (
    TAG_ALIGNERS,
    TAG_MAPQ,
    TAG_PAYLOAD,
    TAG_QSLICE,
    _STRUCTURAL_TAGS,
    hclip_bounds,
    revcomp,
)

logger = logging.getLogger(__name__)


def _read_key(rec: pysam.AlignedSegment):
    """Group key: the canonical (already-normalized) QNAME the writer stored on
    every record of a read. Contiguity is guaranteed by the writer, so streaming
    grouping is exact."""
    return rec.query_name


def _strip_structural_tags(rec: pysam.AlignedSegment) -> None:
    for t in _STRUCTURAL_TAGS:
        if rec.has_tag(t):
            rec.set_tag(t, None)


def _reconstruct_seq_qual(
    rec: pysam.AlignedSegment,
    payload_seq: str,
    payload_qual,
    payload_is_reverse: bool,
    payload_len: int,
):
    """Return (seq, qual_array) for a SEQ=* record from the payload."""
    same = rec.is_reverse == payload_is_reverse
    base_seq = payload_seq if same else revcomp(payload_seq)
    if payload_qual is None:
        base_qual = None
    else:
        base_qual = payload_qual if same else payload_qual[::-1]
    if rec.has_tag(TAG_QSLICE):
        qs, qe = (int(x) for x in rec.get_tag(TAG_QSLICE))
    else:
        qs, qe = hclip_bounds(rec.cigartuples, payload_len)
    seq = base_seq[qs:qe]
    qual = None if base_qual is None else array("B", base_qual[qs:qe])
    return seq, qual


def _expand_group(records: List[pysam.AlignedSegment], header) -> Dict[str, pysam.AlignedSegment]:
    # Locate the payload (SEQ-bearing, Zp=1). Fall back to the first SEQ-bearing
    # record for defensiveness.
    payload = None
    for r in records:
        if r.has_tag(TAG_PAYLOAD) and r.query_sequence is not None:
            payload = r
            break
    if payload is None:
        for r in records:
            if r.query_sequence is not None:
                payload = r
                break
    payload_seq = payload.query_sequence if payload is not None else None
    payload_qual = payload.query_qualities if payload is not None else None
    payload_rev = bool(payload.is_reverse) if payload is not None else False
    payload_len = len(payload_seq) if payload_seq is not None else 0

    out: Dict[str, pysam.AlignedSegment] = {}
    for rec in records:
        aligners = (rec.get_tag(TAG_ALIGNERS).split(",") if rec.has_tag(TAG_ALIGNERS)
                    else ["?"])
        if rec.has_tag(TAG_MAPQ):
            mapqs = [int(x) for x in rec.get_tag(TAG_MAPQ)]
        else:
            mapqs = [rec.mapping_quality] * len(aligners)

        if rec.query_sequence is None and rec.has_tag(TAG_QSLICE):
            # A compressed variant — slice SEQ/QUAL out of the payload.
            seq, qual = _reconstruct_seq_qual(
                rec, payload_seq, payload_qual, payload_rev, payload_len
            )
        else:
            # SEQ present (payload / uncompressed fallback record), OR a
            # genuinely SEQ-less record with no donor (Zq absent) — keep as-is,
            # which may legitimately be SEQ=* (planning/254 §8-Q2).
            seq, qual = rec.query_sequence, rec.query_qualities

        for i, aligner in enumerate(aligners):
            clone = pysam.AlignedSegment.fromstring(rec.to_string(), header)
            clone.flag = clone.flag & ~0x900  # restore primary
            # Restore SEQ/QUAL (order matters: setting SEQ clears QUAL). Never
            # assign SEQ=None onto a record that is already SEQ=* (a no-op that
            # would still wipe QUAL); skip when there is nothing to restore.
            if seq is not None and clone.query_sequence != seq:
                clone.query_sequence = seq
                clone.query_qualities = qual
            clone.mapping_quality = mapqs[i] if i < len(mapqs) else rec.mapping_quality
            _strip_structural_tags(clone)
            out[aligner] = clone
    return out


def expand(cma_path: str, genome: Optional[Dict[str, str]] = None) -> Iterator[Tuple[object, Dict[str, pysam.AlignedSegment]]]:
    """Stream ``(read_key, {aligner: record})`` from a CMA BAM.

    Each yielded dict is field-equivalent (planning/254 §3.3 load-bearing view)
    to what the K-way merge produced pre-CMA, so extract_alignment_info /
    select_best_alignment consume it unchanged. ``genome`` is accepted for API
    parity (MD/NM ride the records verbatim in v1, so it is currently unused).
    """
    with pysam.AlignmentFile(cma_path, "rb") as cma:
        header = cma.header
        cur_key = object()
        buf: List[pysam.AlignedSegment] = []
        for rec in cma:
            k = _read_key(rec)
            if buf and k != cur_key:
                yield cur_key, _expand_group(buf, header)
                buf = []
            cur_key = k
            buf.append(rec)
        if buf:
            yield cur_key, _expand_group(buf, header)
