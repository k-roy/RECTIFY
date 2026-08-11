"""
Tests for the simple-slide fast path in _apply_junction_replacement.

When a candidate junction is a pure k-bp slide of the current N-op and the
genome bases that flip role between intron and exon are identical at both
placements, the CIGAR transformation is trivial: extend the boundary =/M op
on one side of N by k, shrink the boundary =/M op on the other side by k.

This regression pins the cat3_plus_2 mapPacBio case (chrI:142252-142618 with
G==G at both endpoints, sliding +1 to canonical GT|AG at chrI:142253-142619).
The expected refiner output is exactly:
    14=1D9=366N50=...
and NOT the boundary-D/I form (14=1D8=1D366N1I50=...) that the general path
emits, which the downstream indel corrector mangles into degenerate output.

Author: Kevin R. Roy
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pytest
import pysam

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.splice.junction_refiner import _apply_junction_replacement
from rectify.core.splice.junction_scoring import (
    _D, _EQ, _H, _I, _M, _N, _S, _X,
    _QUERY_CONSUMING, _REF_CONSUMING,
)

GENOME_PATH = (
    RECTIFY_ROOT
    / "rectify/data/genomes/saccharomyces_cerevisiae"
    / "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
)


# ---------------------------------------------------------------------------
# Genome fixture (chrI only)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def chrI_seq() -> str:
    """Load chrI sequence from the bundled S. cerevisiae genome."""
    if not GENOME_PATH.exists():
        pytest.skip(f"Genome not found: {GENOME_PATH}")
    # Avoid the BioPython dependency: parse FASTA directly.
    with gzip.open(GENOME_PATH, "rt") as fh:
        cur_name = None
        chunks: List[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur_name == "chrI":
                    return "".join(chunks).upper()
                cur_name = line[1:].split()[0]
                chunks = []
            else:
                if cur_name == "chrI":
                    chunks.append(line.strip())
        if cur_name == "chrI":
            return "".join(chunks).upper()
    pytest.skip("chrI not found in bundled genome")


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

_CIGAR_OP_BY_CHAR = {"M": _M, "I": _I, "D": _D, "N": _N, "S": _S, "H": _H, "=": _EQ, "X": _X}
_CIGAR_CHAR_BY_OP = {v: k for k, v in _CIGAR_OP_BY_CHAR.items()}


def parse_cigar_str(cigar_str: str) -> List[Tuple[int, int]]:
    """Parse a CIGAR string like '14=1D8=' into [(op_code, length), ...]."""
    out: List[Tuple[int, int]] = []
    num = ""
    for ch in cigar_str:
        if ch.isdigit():
            num += ch
        else:
            if ch not in _CIGAR_OP_BY_CHAR:
                raise ValueError(f"Unknown CIGAR op char: {ch}")
            out.append((_CIGAR_OP_BY_CHAR[ch], int(num)))
            num = ""
    return out


def cigar_to_str(cigar: List[Tuple[int, int]]) -> str:
    """Format CIGAR tuples back to a string."""
    return "".join(f"{length}{_CIGAR_CHAR_BY_OP[op]}" for op, length in cigar)


def build_query_for_cigar(
    cigar: List[Tuple[int, int]],
    ref_seq: str,
    ref_start: int,
) -> str:
    """Build a query string consistent with the given CIGAR against ref_seq.

    For =/M ops: copy ref bases (so they are exact matches).
    For X ops:   substitute a non-matching base (cyclic A/C/G/T fallback).
    For I/S ops: use 'A'.
    For D/N ops: no query bases consumed.
    """
    q_parts: List[str] = []
    r_pos = ref_start
    for op, length in cigar:
        if op in (_EQ, _M):
            q_parts.append(ref_seq[r_pos : r_pos + length])
            r_pos += length
        elif op == _X:
            # Mismatch: pick the next base (A→C, C→G, G→T, T→A) so query != ref.
            sub_map = {"A": "C", "C": "G", "G": "T", "T": "A", "N": "A"}
            for i in range(length):
                q_parts.append(sub_map.get(ref_seq[r_pos + i], "A"))
            r_pos += length
        elif op == _D:
            r_pos += length
        elif op == _N:
            r_pos += length
        elif op in (_I, _S):
            q_parts.append("A" * length)
        elif op == _H:
            pass  # H consumes nothing
        else:
            raise ValueError(f"Unhandled op: {op}")
    return "".join(q_parts)


def make_aligned_segment(
    cigar_str: str,
    ref_start: int,
    ref_seq: str,
    qname: str = "test_read",
    is_reverse: bool = False,
) -> pysam.AlignedSegment:
    """Build a pysam.AlignedSegment with a chrI header and the given CIGAR."""
    header = pysam.AlignmentHeader.from_references(["chrI"], [len(ref_seq)])
    read = pysam.AlignedSegment(header)
    read.query_name = qname
    read.reference_id = 0
    read.reference_start = ref_start
    read.mapping_quality = 60
    read.is_unmapped = False
    read.is_reverse = is_reverse
    cigar = parse_cigar_str(cigar_str)
    seq = build_query_for_cigar(cigar, ref_seq, ref_start)
    read.query_sequence = seq
    read.cigartuples = cigar
    return read


# ---------------------------------------------------------------------------
# Pinned regression: cat3_plus_2 mapPacBio +1 slide
# ---------------------------------------------------------------------------

CAT3_INPUT_CIGAR = (
    "14=1D8=366N51=1I2=2D32=3X18=2X1=2X13=3D38=1X124=1I78=2I99=2X60="
    "1D27=3I55=1D49=1D13=2X9=2D30=3D6=1X32="
)
CAT3_EXPECTED_CIGAR = (
    "14=1D9=366N50=1I2=2D32=3X18=2X1=2X13=3D38=1X124=1I78=2I99=2X60="
    "1D27=3I55=1D49=1D13=2X9=2D30=3D6=1X32="
)
CAT3_REF_START = 142229
CAT3_OLD_NS, CAT3_OLD_NE = 142252, 142618
CAT3_NEW_NS, CAT3_NEW_NE = 142253, 142619
# Index of the N-op in the CIGAR: 0=14=, 1=1D, 2=8=, 3=366N
CAT3_N_CIGAR_IDX = 3


def test_simple_slide_genome_bases_match(chrI_seq):
    """Sanity: the slide window has identical genome bases (G == G)."""
    assert chrI_seq[CAT3_OLD_NS] == chrI_seq[CAT3_OLD_NE] == "G"


def test_cat3_plus_2_mapPacBio_slide_cigar(chrI_seq):
    """Trivial +1 slide for cat3_plus_2: 8= → 9=, 51= → 50=.

    Regression for the bug where the general path emitted 14=1D8=1D366N1I50=,
    which the downstream indel corrector then mangled into
    22D1M20I1M366N1I50=.  The fast path emits the clean transformation
    14=1D9=366N50=.
    """
    read = make_aligned_segment(
        CAT3_INPUT_CIGAR,
        ref_start=CAT3_REF_START,
        ref_seq=chrI_seq,
    )
    pre_query = read.query_sequence

    modified = _apply_junction_replacement(
        read,
        cigar_idx=CAT3_N_CIGAR_IDX,
        old_ns=CAT3_OLD_NS,
        old_ne=CAT3_OLD_NE,
        new_ns=CAT3_NEW_NS,
        new_ne=CAT3_NEW_NE,
        genome_seq=chrI_seq,
        strand="+",
        hp_pen=0.25,
        W=15,
    )
    assert modified, "_apply_junction_replacement should have modified the CIGAR"
    out_cigar_str = cigar_to_str(list(read.cigartuples))
    assert out_cigar_str == CAT3_EXPECTED_CIGAR, (
        f"Expected:\n  {CAT3_EXPECTED_CIGAR}\nGot:\n  {out_cigar_str}"
    )
    # reference_start is unchanged
    assert read.reference_start == CAT3_REF_START
    # query span unchanged (pure slide does not change which bases are aligned)
    assert read.query_sequence == pre_query

    # Verify the new N-op spans (142253, 142619)
    r_pos = read.reference_start
    n_ops_seen = []
    for op, length in read.cigartuples:
        if op == _N:
            n_ops_seen.append((r_pos, r_pos + length))
        if op in _REF_CONSUMING or op == _N:
            r_pos += length
    assert (CAT3_NEW_NS, CAT3_NEW_NE) in n_ops_seen, (
        f"Expected N-op at ({CAT3_NEW_NS}, {CAT3_NEW_NE}); got {n_ops_seen}"
    )


def test_minus_strand_slide_uses_same_fast_path(chrI_seq):
    """The fast path is strand-agnostic: a +1 slide on the same coordinates
    must produce the same CIGAR transformation regardless of is_reverse.

    CIGAR ops live in genomic order, so the structural transformation is
    identical between strands.
    """
    read = make_aligned_segment(
        CAT3_INPUT_CIGAR,
        ref_start=CAT3_REF_START,
        ref_seq=chrI_seq,
        is_reverse=True,
    )
    modified = _apply_junction_replacement(
        read,
        cigar_idx=CAT3_N_CIGAR_IDX,
        old_ns=CAT3_OLD_NS,
        old_ne=CAT3_OLD_NE,
        new_ns=CAT3_NEW_NS,
        new_ne=CAT3_NEW_NE,
        genome_seq=chrI_seq,
        strand="-",
        hp_pen=0.25,
        W=15,
    )
    assert modified
    assert cigar_to_str(list(read.cigartuples)) == CAT3_EXPECTED_CIGAR


def test_negative_slide_falls_through_if_bases_differ(chrI_seq):
    """When the slide window bases do NOT match (e.g. delta_start = -1 and
    genome[old_ns - 1] != genome[old_ne - 1]), the fast path must NOT fire.

    This pins the safety of the fast path: it must only fire when the slide
    is provably equivalent.

    UPDATED for the compensating-indel invariant (e40ca00, merged 2026-08-09):
    the general path used to REALIZE this both-boundary move by emitting
    compensating boundary I/D ops (modified=True, non-trivial CIGAR). That
    realization is exactly the fabrication signature the invariant refuses —
    the exon sequence never moves, only the N-op coordinate does (~91% of
    apparent fabrication on real DRS; dev/REPLACER_COMPENSATING_INDEL_BUG.md).
    The correct outcome is now stronger: the unsupported move is REFUSED
    outright and the read keeps its original CIGAR.
    """
    # Pick a coord pair where genome[ns-1] != genome[ne-1].
    # We need an N-op on chrI where genome[ns-1] != genome[ne-1].
    # The cat3 case has genome[142252]=G, genome[142618]=G, so genome at
    # ns-1=142251 and ne-1=142617 may differ — check it.
    if chrI_seq[CAT3_OLD_NS - 1] == chrI_seq[CAT3_OLD_NE - 1]:
        pytest.skip("test coordinates happen to match — pick different ones")

    # Try a -1 slide on the same N-op (new_ns = 142251, new_ne = 142617).
    new_ns = CAT3_OLD_NS - 1
    new_ne = CAT3_OLD_NE - 1

    read = make_aligned_segment(
        CAT3_INPUT_CIGAR,
        ref_start=CAT3_REF_START,
        ref_seq=chrI_seq,
    )
    modified = _apply_junction_replacement(
        read,
        cigar_idx=CAT3_N_CIGAR_IDX,
        old_ns=CAT3_OLD_NS,
        old_ne=CAT3_OLD_NE,
        new_ns=new_ns,
        new_ne=new_ne,
        genome_seq=chrI_seq,
        strand="+",
        hp_pen=0.25,
        W=15,
    )
    # The fast path must not fire (bases differ), and the general path must
    # REFUSE the compensating-indel realization: no move, read unchanged.
    assert modified is False
    out = cigar_to_str(list(read.cigartuples))
    trivial = (
        "14=1D7=366N52=1I2=2D32=3X18=2X1=2X13=3D38=1X124=1I78=2I99=2X60="
        "1D27=3I55=1D49=1D13=2X9=2D30=3D6=1X32="
    )
    assert out != trivial, (
        "Fast path fired when genome bases differ — this is unsafe."
    )
    assert out == CAT3_INPUT_CIGAR, (
        "Refused move must leave the read's original CIGAR untouched."
    )
