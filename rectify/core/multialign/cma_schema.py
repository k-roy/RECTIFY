"""CMA v1 — schema constants, collapse key, and reconstruction helpers.

Losslessness contract (planning/254 §3.3): expand() reproduces, per aligner, the
load-bearing FIELD VIEW — QNAME / FLAG / RNAME / POS / MAPQ / CIGAR / SEQ / QUAL
plus the consumed aux subset {NM, MD, RN, RX, pt, ...}. Aligner-internal scratch
tags (AS, ms, de, ...) are NOT contracted to round-trip (they are preserved
best-effort on the group representative but not per-aligner).
"""

from __future__ import annotations

from typing import Dict, Iterable, Optional, Tuple

import pysam

SCHEMA_VERSION = "CMA:v1"

# --- Tag schema: all in the verified-free Z* user namespace (planning/254 §2.3).
TAG_ALIGNERS = "Za"  # Z   : comma-list of aligners at THIS placement
TAG_PAYLOAD = "Zp"   # i   : 1 on the single SEQ/QUAL-bearing payload record
TAG_K = "Zk"         # i   : number of distinct placements for this read
TAG_QSLICE = "Zq"    # B,i : [qs, qe] slice of payload SEQ this record covers
TAG_MAPQ = "Zm"      # B,C : per-aligner MAPQ, order-matched to Za (only when it varies)

_STRUCTURAL_TAGS = (TAG_ALIGNERS, TAG_PAYLOAD, TAG_K, TAG_QSLICE, TAG_MAPQ)

# Read-intrinsic tags ride the payload once (identical across a read's aligners).
READ_INTRINSIC_TAGS = frozenset(
    {
        "RN",
        "RX",
        "pt",
        # cDNA comment-tag block (rectify/core/cdna/io.py)
        "XU", "XO", "XC", "XR", "XM", "XF", "XA", "XT", "XY", "XB", "XQ", "XK",
    }
)

# Aux tags that MUST round-trip through the field view (used by tests + the
# deletion gate). Xa/Xc/Xn/Xt are post-correct selection metadata, absent from
# the pre-correct per-aligner BAMs the CMA stores, so they are not listed here.
FIELD_VIEW_AUX = ("NM", "MD", "RN", "RX", "pt")

# Preferred payload donors: minimap2 / mapPacBio always soft-clip and retain the
# full read in SEQ (planning/254 §2.4). gapmm2 (PAF→BAM, SEQ=None) and deSALT
# (hard-clip truncation, in production) are last resorts.
PAYLOAD_DONOR_PREFERENCE = ("minimap2", "mapPacBio", "uLTRA", "deSALT", "gapmm2")

# CIGAR op codes
_OP_H = 5
_CONSUMES_QUERY = frozenset({0, 1, 4, 7, 8})  # M I S = X

_COMPLEMENT = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                            "TGCANtgcanYRSWMKVHDByrswmkvhdb")


def revcomp(seq: str) -> str:
    """Reverse-complement (IUPAC-aware)."""
    return seq.translate(_COMPLEMENT)[::-1]


def decode_eq_seq(read: pysam.AlignedSegment, genome=None):
    """Return the read SEQ with any SAM ``=`` (match-to-reference) bases resolved
    to explicit nucleotides, using ``genome`` (a ``{reference_name: str}`` map).

    Some BAMs (e.g. anything through ``samtools calmd -e``, and the bundled DRS
    fixture) store SEQ with ``=`` where the base equals the reference. That
    encoding is PLACEMENT-RELATIVE, so it is not read-intrinsic and cannot be
    transferred to a variant at a different locus. The CMA payload therefore
    stores fully-explicit nucleotides. On BAMs that already use explicit SEQ
    (the production default) this is a no-op fast path that never touches
    ``genome``.

    ``query_sequence`` is reference-oriented, so ``=`` at query position ``qpos``
    equals the forward-reference base at its aligned ``rpos`` — substitute
    directly (no complement, even for reverse reads). Returns SEQ unchanged if
    ``=`` is present but no reference is available (caller must supply one).
    """
    seq = read.query_sequence
    if seq is None or "=" not in seq:
        return seq
    if genome is None:
        return seq
    ref = genome.get(read.reference_name)
    if ref is None:
        return seq
    decoded = list(seq)
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if qpos is not None and rpos is not None and decoded[qpos] == "=":
            decoded[qpos] = ref[rpos].upper()
    return "".join(decoded)


def collapse_key(read: pysam.AlignedSegment) -> Tuple[str, bool, int, str]:
    """Distinct-placement key.

    Uses reference **name** (not numeric id) so it is independent of @SQ order
    across per-aligner BAM headers. The reference component is REQUIRED so two
    aligners placing a read at the same start+CIGAR+strand on DIFFERENT
    chromosomes do not wrongly collapse (planning/254 §2.2, fix 2026-07-17).
    """
    return (
        read.reference_name,
        bool(read.is_reverse),
        int(read.reference_start),
        read.cigarstring or "",
    )


def cigar_query_len(cigartuples) -> int:
    """Total query-consuming (M/I/S/=/X) bases in a CIGAR."""
    return sum(length for op, length in (cigartuples or []) if op in _CONSUMES_QUERY)


def hclip_bounds(cigartuples, payload_len: int) -> Tuple[int, int]:
    """(qs, qe) into the payload SEQ that a record's SEQ covers.

    Only hard-clips trim the slice; soft-clips ride INSIDE it (their bases are
    present in SEQ). planning/254 §2.4.
    """
    if not cigartuples:
        return 0, payload_len
    lead = cigartuples[0][1] if cigartuples[0][0] == _OP_H else 0
    trail = cigartuples[-1][1] if cigartuples[-1][0] == _OP_H else 0
    return lead, payload_len - trail


def is_full_seq_softclip_donor(read: pysam.AlignedSegment) -> bool:
    """A valid payload donor: SEQ-bearing and soft-clip-only (no hard-clip → the
    entire read is in SEQ)."""
    if read.query_sequence is None:
        return False
    if read.cigartuples and any(op == _OP_H for op, _ in read.cigartuples):
        return False
    return True


def choose_payload_donor(
    aligner_records: Dict[str, pysam.AlignedSegment],
    panel: Optional[Iterable[str]] = None,
):
    """Pick (aligner, record) to bear SEQ/QUAL for a read.

    Prefers full-SEQ soft-clip donors in PAYLOAD_DONOR_PREFERENCE order, then the
    given panel order, then any remaining aligner. Returns (None, None) when no
    record qualifies — the no-full-SEQ-donor case (planning/254 §8-Q2); the
    writer then keeps every record's SEQ verbatim (uncompressed but lossless).
    """
    order = []
    for a in PAYLOAD_DONOR_PREFERENCE:
        if a in aligner_records:
            order.append(a)
    for a in (panel or ()):
        if a in aligner_records and a not in order:
            order.append(a)
    for a in aligner_records:
        if a not in order:
            order.append(a)
    for a in order:
        r = aligner_records.get(a)
        if r is not None and is_full_seq_softclip_donor(r):
            return a, r
    return None, None
