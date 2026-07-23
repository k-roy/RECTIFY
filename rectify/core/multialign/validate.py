"""CMA v1 structural validator.

Checks the invariants expand() relies on (planning/254 §4 validate.py):
  * exactly one payload (Zp=1, SEQ present) per read,
  * each SEQ=* variant's Zq interval length == its CIGAR query-consuming length,
  * Za aligner-sets are disjoint across a read's records,
  * Zm (per-aligner MAPQ) length == len(Za) when present.

Returns a list of human-readable problem strings; empty == valid.
"""

from __future__ import annotations

from typing import Dict, List

import pysam

from .cma_schema import (
    TAG_ALIGNERS,
    TAG_MAPQ,
    TAG_PAYLOAD,
    TAG_QSLICE,
    cigar_query_len,
)


def _validate_read(read_key, records: List[pysam.AlignedSegment]) -> List[str]:
    problems: List[str] = []
    n_payload = 0
    seen_aligners: Dict[str, int] = {}
    for rec in records:
        za = rec.get_tag(TAG_ALIGNERS).split(",") if rec.has_tag(TAG_ALIGNERS) else []
        if not za:
            problems.append(f"{read_key}: record missing {TAG_ALIGNERS}")
        for a in za:
            seen_aligners[a] = seen_aligners.get(a, 0) + 1

        if rec.has_tag(TAG_MAPQ):
            zm = list(rec.get_tag(TAG_MAPQ))
            if len(zm) != len(za):
                problems.append(
                    f"{read_key}: {TAG_MAPQ} len {len(zm)} != {TAG_ALIGNERS} len {len(za)}"
                )

        is_payload = rec.has_tag(TAG_PAYLOAD) and int(rec.get_tag(TAG_PAYLOAD)) == 1
        if is_payload:
            n_payload += 1
            if rec.query_sequence is None:
                problems.append(f"{read_key}: payload record has SEQ=*")
        elif rec.query_sequence is None and rec.has_tag(TAG_QSLICE):
            # a compressed variant — its Zq slice must match the CIGAR query span
            qs, qe = (int(x) for x in rec.get_tag(TAG_QSLICE))
            span = qe - qs
            clen = cigar_query_len(rec.cigartuples)
            if span != clen:
                problems.append(
                    f"{read_key}: {TAG_QSLICE} span {span} != CIGAR query len {clen}"
                )
        # NB: a SEQ=* record WITHOUT Zq is a genuine unreconstructable placement
        # (no-donor fallback, planning/254 §8-Q2) — allowed, not a defect.

    if n_payload != 1:
        problems.append(f"{read_key}: expected exactly 1 payload, found {n_payload}")
    dupes = {a: c for a, c in seen_aligners.items() if c > 1}
    if dupes:
        problems.append(f"{read_key}: aligner(s) appear in >1 placement (Za not disjoint): {dupes}")
    return problems


def validate_cma(cma_path: str) -> List[str]:
    """Validate a CMA BAM; returns a list of problems (empty == valid)."""
    problems: List[str] = []
    with pysam.AlignmentFile(cma_path, "rb") as cma:
        cur_key = object()
        buf: List[pysam.AlignedSegment] = []
        for rec in cma:
            k = rec.query_name
            if buf and k != cur_key:
                problems.extend(_validate_read(cur_key, buf))
                buf = []
            cur_key = k
            buf.append(rec)
        if buf:
            problems.extend(_validate_read(cur_key, buf))
    return problems
