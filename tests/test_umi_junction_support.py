"""Test the artifact-resistant junction support metric: n_distinct_umis.

The whole scientific point of the UMI feature: a genuine novel junction is backed
by MANY distinct UMIs; a PCR-jackpot artifact is ONE molecule amplified (many
reads, few distinct UMIs). aggregate_junctions(umi_tag='RX') exposes exactly that.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.aggregate.junctions import aggregate_junctions


def _header():
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrI", "LN": 5000}],
    })


def _spliced_read(header, qname, umi=None, start=100):
    """A read spanning one junction: 50M 200N 50M starting at `start`."""
    r = pysam.AlignedSegment(header)
    r.query_name = qname
    r.flag = 0
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = "50M200N50M"
    r.query_sequence = "A" * 100
    r.query_qualities = pysam.qualitystring_to_array("I" * 100)
    if umi is not None:
        r.set_tag("RX", umi, value_type="Z")
    return r


def _write_sorted(path, reads, header):
    unsorted = str(path) + ".u.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        for r in reads:
            bam.write(r)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    return path


def test_distinct_umis_separates_real_junction_from_jackpot(tmp_path):
    hdr = _header()
    # A REAL junction at chrI:150-350: 5 reads, 5 DISTINCT UMIs.
    real = [_spliced_read(hdr, f"real{i}", umi=u, start=100)
            for i, u in enumerate(
                ["AAAACCCCGGGG", "TTTTGGGGCCCC", "ACGTACGTACGT",
                 "CCCCAAAATTTT", "GGGGTTTTAAAA"])]
    # A JACKPOT junction at chrI:1150-1350: 5 reads but ONE UMI (amplified).
    jackpot = [_spliced_read(hdr, f"jack{i}", umi="ACACACACACAC", start=1100)
               for i in range(5)]
    bam = _write_sorted(tmp_path / "in.bam", real + jackpot, hdr)

    df = aggregate_junctions(str(bam), umi_tag="RX")
    assert "n_distinct_umis" in df.columns

    real_row = df[df.intron_start == 150].iloc[0]
    jack_row = df[df.intron_start == 1150].iloc[0]
    # Same read count...
    assert real_row.full_junction_reads == jack_row.full_junction_reads == 5
    # ...but distinct-UMI support tells them apart.
    assert real_row.n_distinct_umis == 5
    assert jack_row.n_distinct_umis == 1


def test_no_umi_tag_omits_column(tmp_path):
    hdr = _header()
    bam = _write_sorted(tmp_path / "in.bam",
                        [_spliced_read(hdr, "r", umi="ACGTACGTACGT")], hdr)
    df = aggregate_junctions(str(bam))
    assert "n_distinct_umis" not in df.columns
    assert df.iloc[0].full_junction_reads == 1


def test_reads_without_umi_contribute_zero_distinct(tmp_path):
    hdr = _header()
    bam = _write_sorted(tmp_path / "in.bam",
                        [_spliced_read(hdr, "r1"), _spliced_read(hdr, "r2")], hdr)
    df = aggregate_junctions(str(bam), umi_tag="RX")
    assert df.iloc[0].full_junction_reads == 2
    assert df.iloc[0].n_distinct_umis == 0   # no RX tags present
