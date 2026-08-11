"""668: trim-polya must emit exactly ONE record per molecule (the primary).

Two duplicate-QNAME sources, both live in production until 2026-08-11 (caught by the align
RN-collision guard on wt_rep1: 130,311 dup QNAMEs, 3.2%):
  (a) _trim_region_task fetched by OVERLAP with no start-in-region filter — boundary-spanning
      reads emitted by both neighbouring workers (the 624 defect class, third instance);
  (b) supplementary records (SEQ-bearing) share the primary's QNAME and the trim output is
      rebuilt un-flagged, so they were indistinguishable downstream.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.commands.drs_trim_command import _trim_region_task

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6"},
    "SQ": [{"SN": "chrI", "LN": 200000}],
})


def _rec(name, start, flag=0, seq_len=500):
    r = pysam.AlignedSegment(HEADER)
    r.query_name = name
    r.flag = flag
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigartuples = [(0, seq_len)]
    r.query_sequence = "A" * seq_len
    r.query_qualities = pysam.qualitystring_to_array("?" * seq_len)
    return r


@pytest.fixture()
def toy_bam(tmp_path):
    unsorted = tmp_path / "toy.unsorted.bam"
    p = tmp_path / "toy.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=HEADER) as bam:
        bam.write(_rec("in_region", 1_000))
        bam.write(_rec("boundary_spanner", 99_800, seq_len=500))   # spans the 100k boundary
        bam.write(_rec("in_second_region", 150_000))
        bam.write(_rec("supp_of_in_region", 5_000, flag=0x800))    # supplementary, SEQ-bearing
        bam.write(_rec("in_region", 120_000, flag=0x100))          # secondary dup of in_region
    pysam.sort("-o", str(p), str(unsorted))
    pysam.index(str(p))
    return p


def _names(out_path):
    with pysam.AlignmentFile(str(out_path), "rb", check_sq=False) as bam:
        return [r.query_name for r in bam.fetch(until_eof=True)]


def test_region_task_one_record_per_molecule(toy_bam, tmp_path):
    out1 = tmp_path / "r1.bam"
    out2 = tmp_path / "r2.bam"
    _trim_region_task(str(toy_bam), "chrI", 0, 100_000, str(out1), 0.1, 3, 100)
    _trim_region_task(str(toy_bam), "chrI", 100_000, 200_000, str(out2), 0.1, 3, 100)
    n1, n2 = _names(out1), _names(out2)
    # boundary spanner: starts in region 1 -> region 1 ONLY (was: both)
    assert "boundary_spanner" in n1 and "boundary_spanner" not in n2
    # supplementary + secondary excluded entirely
    assert "supp_of_in_region" not in n1 + n2
    combined = n1 + n2
    assert sorted(combined) == sorted(set(combined)), "duplicate QNAMEs across regions"
    assert set(combined) == {"in_region", "boundary_spanner", "in_second_region"}
