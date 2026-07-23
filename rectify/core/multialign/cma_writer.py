"""CMA v1 writer — build a compressed-multialign BAM from per-aligner placements.

Input is a stream of ``(read_key, {aligner: pysam.AlignedSegment})`` — the same
per-read grouping the K-way consensus merge produces
(consensus._iter_name_grouped_bams). For each read the writer:

  * groups the aligners' records by ``collapse_key`` (distinct placements),
  * chooses ONE payload donor (full-SEQ soft-clip record) that stores SEQ/QUAL,
  * emits one payload record (primary, SEQ/QUAL/MD + read-intrinsic-tag union),
  * emits one SEQ=* variant record per other distinct placement (Zq slice),
  * records per-aligner MAPQ (Zm) whenever it varies inside a shared placement.

Each read's records are emitted CONTIGUOUSLY (payload first), so expand() can
stream-group them without a global sort (planning/254 §2.1).
"""

from __future__ import annotations

import logging
from array import array
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import pysam

from .cma_schema import (
    SCHEMA_VERSION,
    TAG_ALIGNERS,
    TAG_K,
    TAG_MAPQ,
    TAG_PAYLOAD,
    TAG_QSLICE,
    READ_INTRINSIC_TAGS,
    choose_payload_donor,
    collapse_key,
    decode_eq_seq,
    hclip_bounds,
)

logger = logging.getLogger(__name__)


def _clone(rec: pysam.AlignedSegment, header) -> pysam.AlignedSegment:
    """Deep copy a record into ``header`` via SAM text (RNAME resolves by name)."""
    return pysam.AlignedSegment.fromstring(rec.to_string(), header)


def _panel_order(aligners: List[str], panel: Optional[Iterable[str]]) -> List[str]:
    ordered = [a for a in (panel or ()) if a in aligners]
    for a in sorted(aligners):
        if a not in ordered:
            ordered.append(a)
    return ordered


def _collect_intrinsic_union(aligner_records) -> List[Tuple[str, object, str]]:
    """Union of read-intrinsic tags across a read's records (first value wins)."""
    seen: Dict[str, Tuple[str, object, str]] = {}
    for r in aligner_records.values():
        for tag, val, typ in r.get_tags(with_value_type=True):
            if tag in READ_INTRINSIC_TAGS and tag not in seen:
                seen[tag] = (tag, val, typ)
    return list(seen.values())


def _apply_intrinsic_union(rec: pysam.AlignedSegment, intrinsic) -> None:
    existing = {t for t, _ in rec.get_tags()}
    for tag, val, typ in intrinsic:
        if tag not in existing:
            rec.set_tag(tag, val, value_type=typ)


def _augment_header(template_header, panel: Iterable[str]):
    hdr = template_header.to_dict()
    hdr.setdefault("CO", []).append(f"{SCHEMA_VERSION} panel={','.join(panel)}")
    hdr.setdefault("PG", []).append(
        {
            "ID": "rectify-cma",
            "PN": "rectify",
            "DS": f"CMA {SCHEMA_VERSION} compressed-multialign writer",
        }
    )
    return pysam.AlignmentHeader.from_dict(hdr)


def _emit_read_records(aligner_records, panel, header, genome=None) -> List[pysam.AlignedSegment]:
    if not aligner_records:
        return []

    from ..consensus.consensus import _normalize_bam_read_name

    order = _panel_order(list(aligner_records), panel)
    groups: Dict[tuple, List[Tuple[str, pysam.AlignedSegment]]] = {}
    for a in order:
        r = aligner_records[a]
        groups.setdefault(collapse_key(r), []).append((a, r))
    K = len(groups)

    payload_aligner, payload_rec = choose_payload_donor(aligner_records, panel)
    have_donor = payload_rec is not None
    payload_key = collapse_key(payload_rec) if have_donor else None
    payload_len = len(payload_rec.query_sequence) if have_donor else None

    canonical = _normalize_bam_read_name(next(iter(aligner_records.values())).query_name or "")
    intrinsic = _collect_intrinsic_union(aligner_records)

    out: List[pysam.AlignedSegment] = []
    payload_emitted = False
    for gkey, members in groups.items():
        aligners = [a for a, _ in members]
        mapqs = [rc.mapping_quality for _, rc in members]
        is_payload = have_donor and (not payload_emitted) and gkey == payload_key
        src = payload_rec if is_payload else members[0][1]
        rec = _clone(src, header)
        rec.query_name = canonical

        if is_payload:
            rec.flag = rec.flag & ~0x900  # force primary
            # Store fully-explicit nucleotides (resolve any '=' match-encoding),
            # so variants at OTHER placements reconstruct real bases (no-op on
            # explicit-SEQ production BAMs). Setting SEQ clears QUAL → save/restore.
            decoded = decode_eq_seq(payload_rec, genome)
            if decoded is not None and decoded != rec.query_sequence:
                q = rec.query_qualities
                rec.query_sequence = decoded
                rec.query_qualities = q
            rec.set_tag(TAG_PAYLOAD, 1, "i")
            rec.set_tag(TAG_K, K, "i")
            _apply_intrinsic_union(rec, intrinsic)
            payload_emitted = True
        else:
            rec.flag = rec.flag | 0x100  # structural secondary
            if have_donor:
                qs, qe = hclip_bounds(rec.cigartuples, payload_len)
                rec.query_sequence = None  # SEQ=* (also clears QUAL)
                rec.set_tag(TAG_QSLICE, array("i", [qs, qe]))  # B,i
            # else: no-donor fallback keeps this record's own SEQ (uncompressed)

        rec.set_tag(TAG_ALIGNERS, ",".join(aligners), "Z")
        if len(set(mapqs)) > 1:
            rec.set_tag(TAG_MAPQ, array("B", mapqs))  # B,C per-aligner MAPQ
        out.append(rec)

    if not payload_emitted:
        # No-donor fallback (planning/254 §8-Q2): no full-SEQ soft-clip record
        # exists. Every record kept its SEQ above; anoint the first SEQ-BEARING
        # record as payload so "exactly one payload per read" holds and the
        # payload is never SEQ=*. (If truly none is SEQ-bearing — a fully
        # SEQ-less read — the first record is used; expand yields it SEQ=*.)
        anchor = next((r for r in out if r.query_sequence is not None), out[0])
        anchor.flag = anchor.flag & ~0x900
        anchor.set_tag(TAG_PAYLOAD, 1, "i")
        anchor.set_tag(TAG_K, K, "i")
        _apply_intrinsic_union(anchor, intrinsic)
        logger.debug("CMA no-donor fallback for read %s (kept %d SEQ-bearing records)", canonical, len(out))

    return out


def build_cma(read_stream, template_header, out_path: str, panel: Iterable[str], genome=None) -> Dict[str, int]:
    """Write a CMA BAM from a ``(read_key, {aligner: record})`` stream.

    ``genome`` ({reference_name: str}) is only consulted to decode SAM ``=``
    match-encoded SEQ (calmd -e / the DRS fixture); pass it whenever inputs may
    be ``=``-encoded. Returns ``{'reads': n, 'records': m}``.
    """
    panel = list(panel)
    header = _augment_header(template_header, panel)
    n_reads = n_records = 0
    with pysam.AlignmentFile(out_path, "wb", header=header) as out:
        for _read_key, aligner_records in read_stream:
            for rec in _emit_read_records(aligner_records, panel, header, genome):
                out.write(rec)
                n_records += 1
            n_reads += 1
    logger.info("CMA written: %s (%d reads, %d records)", out_path, n_reads, n_records)
    return {"reads": n_reads, "records": n_records}


def build_cma_from_bams(aligner_bams: Dict[str, str], out_path: str, panel=None, genome=None) -> Dict[str, int]:
    """Build a CMA from a ``{aligner: bam_path}`` map, streaming and scale-safe.

    Name-sorts inputs as needed (the K-way merge requires queryname order), then
    streams them through the production grouper (consensus._iter_name_grouped_bams)
    into build_cma. Used by ``rectify cma build`` and ``align --emit-cma``.
    """
    import os
    import tempfile

    from ..consensus.consensus import _iter_name_grouped_bams

    panel = list(panel) if panel else list(aligner_bams)
    with tempfile.TemporaryDirectory(prefix="cma_build_") as td:
        sorted_paths = {}
        for aligner, path in aligner_bams.items():
            with pysam.AlignmentFile(path, "rb") as f:
                so = f.header.to_dict().get("HD", {}).get("SO")
            if so == "queryname":
                sorted_paths[aligner] = path
            else:
                sp = os.path.join(td, os.path.basename(path) + ".nsort.bam")
                pysam.sort("-n", "-o", sp, path)
                sorted_paths[aligner] = sp
        first = next(iter(sorted_paths.values()))
        with pysam.AlignmentFile(first, "rb") as f:
            header = f.header
        return build_cma(_iter_name_grouped_bams(sorted_paths), header, out_path, panel, genome=genome)


def load_aligner_records(bam_paths: Dict[str, str]) -> Iterator[Tuple[object, Dict[str, pysam.AlignedSegment]]]:
    """In-memory per-read grouper for small inputs (fixtures, migration ingest).

    Reads PRIMARY mapped records, groups by integer RN when every input carries
    it, else by normalized QNAME. Yields ``(read_key, {aligner: record})`` in
    sorted key order. Production DAG wiring uses
    consensus._iter_name_grouped_bams on name-sorted BAMs instead (planning/254 §4).
    """
    from ..consensus.consensus import _normalize_bam_read_name, _first_primary_has_rn

    use_rn = bool(bam_paths) and all(_first_primary_has_rn(p) for p in bam_paths.values())
    grouped: Dict[object, Dict[str, pysam.AlignedSegment]] = {}
    for aligner, path in bam_paths.items():
        with pysam.AlignmentFile(path, "rb") as bam:
            for r in bam:
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if use_rn:
                    try:
                        key = int(r.get_tag("RN"))
                    except KeyError:
                        key = _normalize_bam_read_name(r.query_name or "")
                else:
                    key = _normalize_bam_read_name(r.query_name or "")
                grouped.setdefault(key, {}).setdefault(aligner, r)
    for key in sorted(grouped, key=lambda k: (isinstance(k, str), k)):
        yield key, grouped[key]
