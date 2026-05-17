#!/usr/bin/env python3
"""Tests for chimeric-read short-circuit in :func:`correct_read_3prime`.

Pins the contract that ``Xz=1`` reads, while skipped by the correction
pipeline, still report their N-op junctions in the returned dict so they
flow through to ``corrected_reads.tsv`` for the validation tests and the
junction-aggregation step.
"""
from __future__ import annotations

import os
import sys

import pytest
from unittest.mock import MagicMock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rectify.core.bam.bam_processor import correct_read_3prime


def _make_chimeric_read(
    cigar_tuples,
    ref_start: int,
    ref_name: str = 'chrI',
    is_reverse: bool = False,
    xz: int = 1,
    query_name: str = 'chim_read',
):
    """Minimal mock of a chimeric-stitched pysam.AlignedSegment."""
    read = MagicMock()
    read.query_name = query_name
    read.cigartuples = cigar_tuples
    read.reference_start = ref_start
    # reference_end = ref_start + sum of ref-consuming op lengths
    _ref_consuming = {0, 2, 3, 7, 8}  # M, D, N, =, X
    read.reference_end = ref_start + sum(
        length for op, length in cigar_tuples if op in _ref_consuming
    )
    read.reference_name = ref_name
    _query_consuming = {0, 1, 4, 7, 8}  # M, I, S, =, X
    read.query_length = sum(
        length for op, length in cigar_tuples if op in _query_consuming
    )
    read.query_sequence = 'A' * read.query_length
    read.query_qualities = None
    read.is_reverse = is_reverse
    read.is_unmapped = False
    read.is_secondary = False
    read.is_supplementary = False

    # Tag lookup: Xz returns the requested value; anything else raises KeyError.
    def _get_tag(name):
        if name == 'Xz':
            return xz
        raise KeyError(name)
    read.get_tag.side_effect = _get_tag
    read.has_tag.side_effect = lambda n: n == 'Xz'
    return read


def test_chimeric_short_circuit_populates_junctions():
    """A chimeric (Xz=1) read with a single N op must surface that junction
    in the returned dict — previously these fields were omitted entirely.
    """
    # 50M 100N 50M — a single splice junction at ref [1100, 1200)
    cigar = [(0, 50), (3, 100), (0, 50)]
    read = _make_chimeric_read(cigar, ref_start=1050)

    rows = correct_read_3prime(read, genome={})

    assert len(rows) == 1, 'chimeric short-circuit must return a single row'
    row = rows[0]
    # The three fields the previous early-return omitted:
    assert 'junctions' in row
    assert 'junctions_str' in row
    assert 'n_junctions' in row
    # 50M starts at 1050 → ends at 1100; 100N spans [1100, 1200); 50M resumes at 1200.
    # extract_junctions_simple returns 0-based intron coordinates.
    assert row['n_junctions'] == 1
    assert row['junctions'] == [(1100, 1200)]
    # junctions_str is a human-readable form; assert it carries the same coords.
    assert '1100' in row['junctions_str']
    assert '1200' in row['junctions_str']


def test_chimeric_short_circuit_multi_junction():
    """Two N ops produce two junctions; both must appear."""
    # 30M 50N 20M 80N 30M
    cigar = [(0, 30), (3, 50), (0, 20), (3, 80), (0, 30)]
    read = _make_chimeric_read(cigar, ref_start=1000)

    rows = correct_read_3prime(read, genome={})

    assert len(rows) == 1
    row = rows[0]
    assert row['n_junctions'] == 2
    assert row['junctions'] == [(1030, 1080), (1100, 1180)]


def test_chimeric_short_circuit_no_junction():
    """A chimeric stitch without an N op (rare but possible) reports zero
    junctions cleanly — empty list, empty string, not missing fields."""
    cigar = [(0, 100)]
    read = _make_chimeric_read(cigar, ref_start=2000)

    rows = correct_read_3prime(read, genome={})

    assert len(rows) == 1
    row = rows[0]
    assert row['n_junctions'] == 0
    assert row['junctions'] == []
    assert row['junctions_str'] == ''
