#!/usr/bin/env python3
"""Tests for the junction-local perfect-match anchor gate.

The gate flags an N-op alignment as suspicious (so it loses winner-selection to a
clean soft-clipping competitor regardless of the free-N-op HP-ED score) when the
exact-match read run immediately abutting the intron is shorter than K — UNLESS
every junction in the read is rescued by the cross-read support relaxation.

Validated empirically on A549 chr5 DRS (real introns ~13 bp min-flank anchor vs
~2 bp for spurious mapPacBio splits); see memory project-hped-anchor-gate.

Author: Kevin R. Roy
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pandas as pd
import pytest
from unittest.mock import MagicMock

from rectify.core.consensus.corrected_consensus import (
    _cigar_min_junction_anchor,
    _add_chimera_flag,
    _NO_JUNCTION_ANCHOR,
    _MIN_JUNCTION_ANCHOR,
)
from rectify.core.splice.calibrate_junction_overhang import OverhangTable

# Deterministic non-repetitive 200 bp genome.
_BASES = "ACGTAGCTGATCGTACGATCGTTAGCATCGATCGTACGTAGCATGCATCGA"
GENOME = {"chr1": (_BASES * 4)[:200]}
G = GENOME["chr1"]

# intron occupies ref [20, 120); flanking exons ref [10,20) and [120,130)
L_EXON = G[10:20]    # 10 bp, abuts the donor
R_EXON = G[120:130]  # 10 bp, abuts the acceptor


def make_read(cigar_tuples, ref_start, query_seq, reference_name="chr1"):
    read = MagicMock()
    read.cigartuples = cigar_tuples
    read.reference_start = ref_start
    read.reference_name = reference_name
    read.query_sequence = query_seq
    return read


# ── _cigar_min_junction_anchor ───────────────────────────────────────────────

def test_clean_junction_full_anchor():
    """Both flanks match the genome exactly -> anchor = flank length (10)."""
    read = make_read([(7, 10), (3, 100), (7, 10)], 10, L_EXON + R_EXON)
    assert _cigar_min_junction_anchor(read, GENOME) == 10


def test_zero_anchor_when_junction_adjacent_bases_mismatch():
    """A mismatch immediately abutting the intron on each flank -> anchor 0."""
    # last base of left exon and first base of right exon flipped to a non-matching base
    def flip(b):
        return "A" if b != "A" else "C"
    q = L_EXON[:-1] + flip(L_EXON[-1]) + flip(R_EXON[0]) + R_EXON[1:]
    read = make_read([(7, 10), (3, 100), (7, 10)], 10, q)
    assert _cigar_min_junction_anchor(read, GENOME) == 0


def test_anchor_is_min_of_two_flanks():
    """Left flank matches 10, right flank breaks after 3 -> min anchor = 3."""
    def flip(b):
        return "A" if b != "A" else "C"
    # break the right exon at offset 3 (4th base)
    q = L_EXON + R_EXON[:3] + flip(R_EXON[3]) + R_EXON[4:]
    read = make_read([(7, 10), (3, 100), (7, 10)], 10, q)
    assert _cigar_min_junction_anchor(read, GENOME) == 3


def test_eq_encoded_seq_counts_as_match():
    """BBMap-style '=' SEQ encoding is treated as a match on both flanks."""
    q = ("=" * 10) + ("=" * 10)
    read = make_read([(7, 10), (3, 100), (7, 10)], 10, q)
    assert _cigar_min_junction_anchor(read, GENOME) == 10


def test_no_njunction_returns_sentinel():
    read = make_read([(7, 20)], 10, G[10:30])
    assert _cigar_min_junction_anchor(read, GENOME) == _NO_JUNCTION_ANCHOR


def test_no_genome_returns_sentinel():
    read = make_read([(7, 10), (3, 100), (7, 10)], 10, L_EXON + R_EXON)
    assert _cigar_min_junction_anchor(read, None) == _NO_JUNCTION_ANCHOR


def test_m_op_literal_bases_supported():
    """Aligners that emit op M (0) with literal bases are handled like =/X."""
    read = make_read([(0, 10), (3, 100), (0, 10)], 10, L_EXON + R_EXON)
    assert _cigar_min_junction_anchor(read, GENOME) == 10


# ── _add_chimera_flag gate behaviour ─────────────────────────────────────────

def _rep_df(anchor):
    """One aligner×read row whose junction passes the existing overhang checks."""
    return pd.DataFrame([{
        "read_id": "r1",
        "_aligner": "mapPacBio",
        "junctions": "20-120",
        "alignment_start": 0,
        "alignment_end": 140,
        "hp_edit_distance": 1.0,
        "aligned_bases": 120,
        "min_junction_anchor": anchor,
    }])


OV = OverhangTable.default()


def test_gate_off_by_default_is_inert():
    """K defaults to 0 -> a short-anchor row is NOT flagged (byte-identical behavior)."""
    assert _MIN_JUNCTION_ANCHOR == 0
    flags = _add_chimera_flag(_rep_df(anchor=2), OV, junction_stats={})
    assert int(flags.iloc[0]) == 0


def test_gate_flags_unsupported_short_anchor():
    """K=10: sub-K anchor + no cross-read support -> flagged suspicious."""
    flags = _add_chimera_flag(_rep_df(anchor=2), OV, junction_stats={},
                              min_junction_anchor_bp=10)
    assert int(flags.iloc[0]) == 1


def test_gate_spares_well_supported_short_anchor():
    """K=10: sub-K anchor but the junction is well-supported across reads -> spared."""
    stats = {(20, 120): (5, 10)}  # cnt>=min_support(5), max_ov>=good_overhang(10)
    flags = _add_chimera_flag(_rep_df(anchor=2), OV, junction_stats=stats,
                              min_junction_anchor_bp=10)
    assert int(flags.iloc[0]) == 0


def test_gate_passes_long_anchor():
    """K=10: anchor >= K -> not flagged by the gate."""
    flags = _add_chimera_flag(_rep_df(anchor=15), OV, junction_stats={},
                              min_junction_anchor_bp=10)
    assert int(flags.iloc[0]) == 0


def test_gate_absent_column_is_inert():
    """If the min_junction_anchor column is missing, the gate never fires."""
    df = _rep_df(anchor=2).drop(columns=["min_junction_anchor"])
    flags = _add_chimera_flag(df, OV, junction_stats={}, min_junction_anchor_bp=10)
    assert int(flags.iloc[0]) == 0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
