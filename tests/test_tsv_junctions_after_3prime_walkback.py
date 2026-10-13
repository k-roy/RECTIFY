#!/usr/bin/env python3
"""ISSUE-088: the TSV must list exactly the junctions the written record carries after the 3' walkback.

A walkback that moves the 3' end across an N (a terminal "exon" made of tail) lets the writer's clip remove that N,
but the row kept the aligner's junction list, so the TSV claimed a junction the BAM does not have (CF-1, read
5cac4ae6, 3 of 62,602 Sumner chr5 cohort reads; the writer audit tagged it `jx-:`). read-review §5: TSV junctions
equal the written N ops; the BAM is the record. Each case below clips a read with the writer's own function and
checks the row against the record.
"""

import pysam
import pytest

from rectify.core.bam.bam_processor import _drop_junctions_past_3prime
from rectify.core.bam.bam_writer import _n_op_intervals
from rectify.core.bam.read_edits import clip_read_to_corrected_3prime
from rectify.utils.alignment import format_junctions_string


def _read(start, cigar, flag):
    header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': 100000}]})
    r = pysam.AlignedSegment(header)
    r.query_name = 'r'
    r.flag = flag
    r.reference_id = 0
    r.reference_start = start
    r.cigarstring = cigar
    n = r.query_length
    r.query_sequence = 'A' * n
    r.query_qualities = pysam.qualitystring_to_array('I' * n)
    r.mapping_quality = 60
    return r


def _row(read, strand, corrected_3prime):
    js = list(_n_op_intervals(read))
    return {'strand': strand, 'corrected_3prime': corrected_3prime, 'junctions': js,
            'junctions_str': format_junctions_string(js), 'n_junctions': len(js)}


def _check(read, strand, c3):
    row = _row(read, strand, c3)
    _drop_junctions_past_3prime(row)
    clip_read_to_corrected_3prime(read, c3, strand)
    carried = list(_n_op_intervals(read))
    assert sorted(row['junctions']) == sorted(carried)
    assert row['n_junctions'] == len(carried)
    assert row['junctions_str'] == format_junctions_string(row['junctions'])
    return row


def test_plus_strand_walkback_into_exon_one_drops_the_terminal_junction():
    # 1,000: 40M (exon 1, 1000..1039) 100N (1040..1139) 12M (1140..1151); the walk lands on exon 1's last base
    r = _read(1000, '40M100N12M', 0)
    row = _check(r, '+', 1039)
    assert row['junctions'] == [] and row['n_junctions'] == 0


def test_minus_strand_walkback_to_the_first_base_after_the_intron_drops_it():
    # CF-1's shape: 12M (tail) 100N then the body; the corrected 3' end is the first base after the intron
    r = _read(1000, '12M100N40M', 16)
    row = _check(r, '-', 1112)
    assert row['junctions'] == []


def test_end_inside_the_last_exon_keeps_the_junction():
    r = _read(1000, '40M100N12M', 0)
    row = _check(r, '+', 1140)           # one base into exon 2
    assert row['junctions'] == [(1040, 1140)]


def test_minus_strand_end_left_of_the_intron_keeps_it():
    r = _read(1000, '12M100N40M', 16)
    row = _check(r, '-', 1011)           # last base of the leftmost block
    assert row['junctions'] == [(1012, 1112)]


def test_only_the_junctions_past_the_end_go():
    r = _read(1000, '30M50N30M60N10M', 0)  # introns 1030..1079 and 1110..1169
    row = _check(r, '+', 1100)
    assert row['junctions'] == [(1030, 1080)]


def test_uncorrected_rows_are_untouched():
    r = _read(1000, '40M100N12M', 0)
    row = _row(r, '+', r.reference_end - 1)
    before = dict(row)
    _drop_junctions_past_3prime(row)
    assert row == before
