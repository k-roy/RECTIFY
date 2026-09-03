"""Path-A poly(A): `correct --ONT-cDNA` must fall back to the stage-1 `XA` tail length.

WHY THIS EXISTS (planning/860, measured, not predicted).  `run-all --ONT-cDNA` Path A -- the
DEFAULT for this datatype -- never runs `trim-cdna-polya`, so the `pl` tag the consumer used to
read is absent.  `correct-cdna` stage 1 instead pretrims the tail off the emitted CONSENSUS and
records its length as `XA:i` (`core/cdna/io.py`).  The aligned molecule therefore has no tail
left to measure, exactly as after the Path-B trim that `pl` was introduced to survive.

Measured on 49,831 real in-house PCB114 reads (wtaa_rep1, rectify fd2e2d2):

    stage-1 consensus FASTQ   XA:i present 100.0%, NON-ZERO 95.5%   (per-cluster median 20-30 nt)
    corrected_reads.tsv       polya_source='none' 100.0%
                              polya_length non-zero 27.0%, median 1 nt, only 40/4726 rows >= 20 nt
                              pt_tag EMPTY on 100% of rows

i.e. the tail was measured by stage 1 and then discarded.  This is the failure the module
docstring of `protocols/ont_cdna.py` marks "WHY THIS IS NOT OPTIONAL", reached by a different
route than the Path-B case it documents.

The tests are deliberately structural AND end-to-end: a real-data smoke can pass while silently
zeroing this column -- that is exactly what happened.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.correct.protocols.ont_cdna import (
    POLYA_SOURCE_STAGE1,
    POLYA_SOURCE_TRIM,
    STAGE1_TAIL_LEN_TAG,
    TAIL_LEN_TAG,
    stage1_tail_length,
    trim_stage_tail_length,
)

CHROM = "chrTest"
CHROM_LEN = 6000
# '+'-strand sense read; a non-A/T genome below keeps walkback from moving the terminus.
READ_SPAN = (1800, 2000)


class _FakeRead:
    """Minimal pysam.AlignedSegment stand-in for tag lookup."""

    def __init__(self, tags):
        self._tags = tags

    def get_tag(self, tag):
        if tag not in self._tags:
            raise KeyError(tag)
        return self._tags[tag]


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    )


def _read(pl=None, xa=None, ro="S", is_reverse=False):
    start, end = READ_SPAN
    a = pysam.AlignedSegment(_header())
    a.query_name = "r"
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.flag = 16 if is_reverse else 0
    n = end - start
    a.query_sequence = "C" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.cigartuples = [(0, n)]
    if ro is not None:
        a.set_tag("ro", ro, value_type="A")
    if pl is not None:
        a.set_tag(TAIL_LEN_TAG, int(pl), value_type="i")
    if xa is not None:
        a.set_tag(STAGE1_TAIL_LEN_TAG, int(xa), value_type="i")
    return a


# ---------------------------------------------------------------------------
# The reader
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value,expected", [(0, 0), (18, 18), (42, 42), ("25", 25)])
def test_stage1_tail_length_reads_the_tag(value, expected):
    assert stage1_tail_length(_FakeRead({STAGE1_TAIL_LEN_TAG: value})) == expected


def test_stage1_tail_length_absent_returns_none_not_zero():
    """Absent tag must be None (fall back), NOT 0 (assert a tail-less molecule).

    Conflating the two would make every non-Path-A BAM silently acquire
    zero-length tails -- the same trap `trim_stage_tail_length` guards against.
    """
    assert stage1_tail_length(_FakeRead({})) is None


def test_stage1_tail_length_zero_is_a_real_zero():
    """`XA:i:0` is a measured zero-length tail and must NOT read as 'absent'."""
    assert stage1_tail_length(_FakeRead({STAGE1_TAIL_LEN_TAG: 0})) == 0
    assert stage1_tail_length(_FakeRead({STAGE1_TAIL_LEN_TAG: 0})) is not None


def test_stage1_tail_length_rejects_garbage():
    assert stage1_tail_length(_FakeRead({STAGE1_TAIL_LEN_TAG: "not_an_int"})) is None
    assert stage1_tail_length(_FakeRead({STAGE1_TAIL_LEN_TAG: -3})) is None


def test_the_two_sources_are_distinguishable():
    """`XA` is a stage-1 CONSENSUS measurement, `pl` a trim-stage one.

    They must not share a `polya_source` value, or a downstream reader cannot
    tell which route produced the number.
    """
    assert POLYA_SOURCE_STAGE1 != POLYA_SOURCE_TRIM
    assert POLYA_SOURCE_STAGE1 == "cdna_stage1"


# ---------------------------------------------------------------------------
# End-to-end through correct_read_3prime -- proves the fallback is actually wired
# ---------------------------------------------------------------------------

class TestCorrectReadWiring:
    GENOME = {CHROM: "C" * CHROM_LEN}

    def _run(self, read, ont_cDNA=True):
        from rectify.core.bam.bam_processor import correct_read_3prime

        out = correct_read_3prime(read, self.GENOME, ont_cDNA=ont_cDNA,
                                  apply_polya_trim=True, apply_3ss_rescue=False)
        assert len(out) == 1
        return out[0]

    def test_xa_populates_polya_length_when_pl_is_absent(self):
        """THE REGRESSION.  Path A writes XA and no `pl`; the tail must survive."""
        row = self._run(_read(pl=None, xa=27))
        assert row["polya_length"] == 27
        assert row["polya_source"] == POLYA_SOURCE_STAGE1

    def test_pl_still_wins_when_both_are_present(self):
        """Path-B precedence is unchanged: `pl` is the trim-stage truth."""
        row = self._run(_read(pl=22, xa=27))
        assert row["polya_length"] == 22
        assert row["polya_source"] == POLYA_SOURCE_TRIM

    def test_xa_zero_is_recorded_as_a_measured_zero(self):
        row = self._run(_read(pl=None, xa=0))
        assert row["polya_length"] == 0
        assert row["polya_source"] == POLYA_SOURCE_STAGE1

    def test_neither_tag_falls_through_to_the_post_alignment_value(self):
        """No `pl`, no `XA` => the old behaviour, and `polya_source` stays 'none'."""
        row = self._run(_read(pl=None, xa=None))
        assert row["polya_source"] == "none"

    def test_the_fallback_is_reachable_only_via_ont_cDNA(self):
        """Regression guard for DRS and every other protocol.

        `XA` is a cDNA stage-1 tag; a DRS run must never pick it up, even if a
        BAM somehow carries one.
        """
        row = self._run(_read(pl=None, xa=27, ro=None), ont_cDNA=False)
        assert row["polya_source"] == "none"
        assert row["polya_length"] != 27 or row["polya_length"] == 0
