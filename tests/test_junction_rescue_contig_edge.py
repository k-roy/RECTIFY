"""Contig-edge refusal in the 5' junction rescue (planning/719).

The bug this pins: the Cat3 5' rescue converted a leading ``25S`` to ``25M``
on a read at POS 0 and walked 25 bases off the START of chrXV, writing a
record at POS -23. pysam writes negative POS without complaint (verified),
``samtools index`` then refuses the WHOLE file, and one read kills a
30-minute correct run — silently, a stage later. A second observed class
overshot the edge by 80,234 bp via a bogus five_prime_position; both funnel
through the same ``new_ref_start`` computation, so one guard catches both.

The fix is a REFUSAL, not a clamp: an extension that leaves the contig means
the proposed placement is impossible, so the read must stay unmodified.
These tests pin refusal on both strands (start and end edges), that the
refusal leaves the read byte-identical, anti-vacuity (the same rescue away
from the edge still fires), and the two-sided write-time invariant in
consensus's checkpoint validation (the original guard checked only
``reference_end > contig_len`` and never fired on a negative POS).
"""

import pysam
import pytest

from rectify.core.bam.read_edits import extend_read_5prime_for_junction_rescue
from rectify.core.consensus.consensus import _validate_bam_sample

CONTIG = 'chrT'
CONTIG_LEN = 1000
HEADER = pysam.AlignmentHeader.from_dict({
    'HD': {'VN': '1.6'},
    'SQ': [{'SN': CONTIG, 'LN': CONTIG_LEN}],
})


def _read(ref_start: int, cigartuples, name: str = 'r1') -> pysam.AlignedSegment:
    r = pysam.AlignedSegment(header=HEADER)
    r.query_name = name
    r.flag = 0
    r.reference_id = 0
    r.reference_start = ref_start
    r.cigartuples = cigartuples
    qlen = sum(l for op, l in cigartuples if op in (0, 1, 4, 7, 8))
    r.query_sequence = 'A' * qlen
    r.mapping_quality = 60
    return r


class TestPlusStrandStartEdge:
    def test_refuses_extension_past_contig_start(self):
        # The live failure shape: POS 0, 25 bp leading clip, rescue target
        # upstream of base 0. Flat fallback exon_ref_span = 25, so
        # new_ref_start = 10 - 25 + 1 = -14.
        read = _read(0, [(4, 25), (0, 50)])
        before = (read.reference_start, list(read.cigartuples))
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=10, soft_clip_len=25, strand='+') is False
        # Refusal must leave the read untouched — a half-applied edit is
        # exactly the malformed-record class the guard exists to prevent.
        assert (read.reference_start, list(read.cigartuples)) == before

    def test_refuses_wild_five_prime_position(self):
        # The 80 kb-overshoot class: a five_prime_position nowhere near the
        # read puts new_ref_start hugely negative. Same guard, same refusal.
        read = _read(0, [(4, 25), (0, 50)])
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=-80_000, soft_clip_len=25, strand='+') is False
        assert read.reference_start == 0

    def test_rescue_away_from_edge_still_fires(self):
        # Anti-vacuity: identical geometry shifted into the contig interior
        # must still rescue (True) with the expected exon/N/body CIGAR.
        read = _read(500, [(4, 25), (0, 50)])
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=400, soft_clip_len=25, strand='+') is True
        assert read.reference_start == 376          # 400 - 25 + 1
        assert read.cigartuples == [(0, 25), (3, 99), (0, 50)]

    def test_boundary_exact_zero_is_allowed(self):
        # new_ref_start == 0 is legal (first base of the contig): fpp=24,
        # exon_ref_span=25 -> new_ref_start = 0. The guard is strictly < 0.
        read = _read(500, [(4, 25), (0, 50)])
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=24, soft_clip_len=25, strand='+') is True
        assert read.reference_start == 0


class TestMinusStrandEndEdge:
    def test_refuses_extension_past_contig_end(self):
        # Body ends at 950; intron 30 + exon 25 would end at 1005 > 1000.
        read = _read(900, [(0, 50), (4, 25)])
        before = (read.reference_start, list(read.cigartuples))
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=980, soft_clip_len=25, strand='-') is False
        assert (read.reference_start, list(read.cigartuples)) == before

    def test_refuses_no_intron_append_past_contig_end(self):
        # intron_len <= 0 path: plain append of 25M past base 1000.
        read = _read(940, [(0, 50), (4, 25)])
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=985, soft_clip_len=25, strand='-') is False
        assert read.cigartuples == [(0, 50), (4, 25)]

    def test_rescue_away_from_edge_still_fires(self):
        read = _read(900, [(0, 50), (4, 25)])
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=960, soft_clip_len=25, strand='-') is True
        assert read.cigartuples == [(0, 50), (3, 10), (0, 25)]
        assert read.reference_end == 985

    def test_boundary_exact_contig_end_is_allowed(self):
        # fpp=975: intron 25, end = 900+50+25+25 = 1000 == CONTIG_LEN — legal
        # (reference_end is half-open). The guard is strictly > contig_len.
        read = _read(900, [(0, 50), (4, 25)])
        assert extend_read_5prime_for_junction_rescue(
            read, five_prime_position=975, soft_clip_len=25, strand='-') is True
        assert read.reference_end == CONTIG_LEN


class TestTwoSidedWriteInvariant:
    """consensus._validate_bam_sample must reject BOTH contig edges."""

    def _write_bam(self, tmp_path, ref_start, cigartuples):
        path = str(tmp_path / 'batch.bam')
        with pysam.AlignmentFile(path, 'wb', header=HEADER) as f:
            f.write(_read(ref_start, cigartuples, name='edge'))
        return path

    def test_negative_pos_is_fatal(self, tmp_path):
        # The planning/719 record class: pysam writes POS -23 without
        # complaint (verified live), so the validator is the last gate
        # before `samtools index` detonates a stage later.
        path = self._write_bam(tmp_path, -23, [(0, 30)])
        with pytest.raises(RuntimeError, match='OFF A CONTIG EDGE'):
            _validate_bam_sample(path)

    def test_past_end_is_still_fatal(self, tmp_path):
        # The original one-sided check must keep working.
        path = self._write_bam(tmp_path, 990, [(0, 30)])
        with pytest.raises(RuntimeError, match='OFF A CONTIG EDGE'):
            _validate_bam_sample(path)

    def test_clean_read_passes(self, tmp_path):
        path = self._write_bam(tmp_path, 100, [(0, 30)])
        _validate_bam_sample(path)  # must not raise
