#!/usr/bin/env python3
"""extract_deletions must never hand pysam an MD tag that disagrees with the CIGAR (planning/746).

pysam's ``get_reference_sequence()`` reconstructs reference bases from MD over the read's full
reference span; when MD and CIGAR disagree its C code overruns the heap — SIGABRT
``free(): invalid next size`` that no try/except can catch, often thousands of reads later.
Long-read splice aligners emit such reads (spurious 100-190 kb "introns", MD shorter than the
CIGAR). The 2026-07-22 guard only skipped max(N) > 50 kb; reads with several N-ops each
<= 50 kb passed it and took down 4/48 DRS samples (2026-08-22). Both fixes are pinned here:
(1) MD-consistency gate, (2) genome-sourced deleted bases that bypass MD entirely.
"""
from __future__ import annotations

from unittest.mock import patch

import pysam

from rectify.utils import alignment as A


def _hdr():
    return pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": 200000}]})


def _read(cigar, md, start=100, seqlen=None):
    r = pysam.AlignedSegment(_hdr())
    r.query_name = "r"
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = cigar
    n = seqlen or sum(l for o, l in r.cigartuples if o in (0, 1, 4, 7, 8))
    r.query_sequence = "A" * n
    if md is not None:
        r.set_tag("MD", md, "Z")
    return r


def test_md_reference_length():
    assert A.md_reference_length("10") == 10
    assert A.md_reference_length("5^AC3") == 10
    assert A.md_reference_length("2G0T4") == 8


def test_consistent_md_is_detected():
    assert A.md_consistent_with_cigar(_read("10M2D5M", "10^GG5"))
    assert not A.md_consistent_with_cigar(_read("10M2D5M", "10^GG3"))
    assert not A.md_consistent_with_cigar(_read("10M", None))


def test_inconsistent_md_never_reaches_pysam():
    """The crash class: three 45 kb N-ops (each under the old 50 kb max-N guard), MD shorter
    than the CIGAR's M/D span. pysam must not be called; the deletion is still reported."""
    r = _read("20M45000N20M45000N10M3D20M45000N20M", "40^GGG20")   # MD claims 63, CIGAR says 93
    with patch.object(A, "_reference_sequence_from_md", side_effect=AssertionError("must not be called")) as m:
        dels = A.extract_deletions(r)
    assert m.call_count == 0
    assert len(dels) == 1 and dels[0]["length"] == 3 and dels[0]["ref_seq"] == ""


def test_consistent_md_uses_pysam_when_no_genome():
    r = _read("10M3D10M", "10^GTC10")
    with patch.object(A, "_reference_sequence_from_md", return_value="A" * 10 + "GTC" + "A" * 10) as m:
        dels = A.extract_deletions(r)
    assert m.call_count == 1
    assert dels[0]["ref_seq"] == "GTC"


def test_genome_path_bypasses_md_entirely():
    """With the genome dict (what the correct-stage worker holds) the deleted bases come from
    the genome and pysam is never consulted — even on an inconsistent MD."""
    genome = {"chrT": "N" * 100 + "A" * 10 + "GTC" + "A" * 10 + "N" * 100}
    r = _read("10M3D10M", "10^GT3")   # inconsistent on purpose
    with patch.object(A, "_reference_sequence_from_md", side_effect=AssertionError("must not be called")) as m:
        dels = A.extract_deletions(r, genome)
    assert m.call_count == 0
    assert dels[0]["ref_pos"] == 110 and dels[0]["ref_seq"] == "GTC"


def test_total_intron_guard_blocks_many_small_nops():
    """Two 30 kb N-ops: each passes a max-N guard, the total does not."""
    r = _read("10M30000N10M30000N10M", "30")    # MD consistent
    with patch.object(A, "_reference_sequence_from_md", side_effect=AssertionError("must not be called")) as m:
        A.extract_deletions(r)
    assert m.call_count == 0
