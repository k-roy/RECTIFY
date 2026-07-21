"""Integration tests for the short-read UMI carry-through: FASTQ RX -> BAM RX:Z.

The COMPASS PE aligners (STAR/HISAT2/magicblast/gsnap) and bbmap/bwa do NOT copy
FASTQ comments into the BAM (only minimap2 -y does), so ``rectify split --umi``
writes the UMI as an ``RX:Z:`` comment token and the runners re-attach it
post-alignment via ``_load_fastq_umi_map`` + ``_inject_rn_into_bam``. These tests
exercise that path on a real (tiny, on-disk) BAM + FASTQ -- no aligner binaries.

Also covers the sidecar RX helpers and the format/parse round-trip.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pysam
import pytest

from rectify.core.align.multi_aligner import _load_fastq_umi_map, _inject_rn_into_bam
from rectify.core.chunking.sidecar import (
    format_fastq_header_with_rn,
    parse_rn_from_fastq_header,
    parse_rx_from_comment,
    parse_rx_from_fastq_header,
)

_REF_NAME = "chrI"


def _header():
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": _REF_NAME, "LN": 1000}],
    })


def _read(header, qname, pos=100, flag=0):
    r = pysam.AlignedSegment(header)
    r.query_name = qname
    r.flag = flag
    r.reference_id = 0
    r.reference_start = pos
    r.mapping_quality = 60
    r.cigarstring = "50M"
    r.query_sequence = "A" * 50
    r.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return r


def _write_bam(path, reads, header):
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for r in reads:
            bam.write(r)
    return path


def _read_tags(path):
    """Return {qname: {tag: value}} for every record in a BAM."""
    out = {}
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for r in bam:
            out.setdefault(r.query_name, []).append(dict(r.get_tags()))
    return out


# ---------------------------------------------------------------------------
# sidecar RX helpers + round-trip
# ---------------------------------------------------------------------------

def test_format_header_includes_rx_when_umi_given():
    h = format_fastq_header_with_rn("READ1", "", read_num=7, umi="ACGTACGTACGT")
    assert "RN:i:7" in h and "RX:Z:ACGTACGTACGT" in h
    assert h.startswith("@READ1\t")


def test_format_header_omits_rx_when_no_umi():
    h = format_fastq_header_with_rn("READ1", "", read_num=7)
    assert "RX:Z:" not in h and "RN:i:7" in h


def test_format_header_backward_compatible_signature():
    """Existing callers pass no umi kwarg -- must behave exactly as before."""
    assert format_fastq_header_with_rn("R", "", 3) == "@R\tRN:i:3\n"


def test_rn_and_rx_coexist_and_parse_back():
    h = format_fastq_header_with_rn("READ1", "CB:Z:cellx", read_num=2, umi="TTTTGGGGCCCC")
    qn_rn, rn = parse_rn_from_fastq_header(h)
    qn_rx, rx = parse_rx_from_fastq_header(h)
    assert qn_rn == qn_rx == "READ1"
    assert rn == 2
    assert rx == "TTTTGGGGCCCC"


def test_rx_survives_alongside_preserved_aux_comment():
    """A pre-existing valid aux token in the comment must not clobber RX."""
    h = format_fastq_header_with_rn("READ1", "CB:Z:cellx", read_num=1, umi="ACGTACGTACGT")
    assert "CB:Z:cellx" in h and "RX:Z:ACGTACGTACGT" in h


def test_parse_rx_from_comment_absent_returns_none():
    assert parse_rx_from_comment("RN:i:5\tCB:Z:x") is None


def test_parse_rx_from_comment_empty_value_returns_none():
    assert parse_rx_from_comment("RX:Z:") is None


# ---------------------------------------------------------------------------
# _load_fastq_umi_map
# ---------------------------------------------------------------------------

def _write_fastq(path, records, gz=False):
    """records: list of (header_line_without_newline, seq)."""
    opener = gzip.open if gz else open
    with opener(str(path), "wt") as fh:
        for header, seq in records:
            fh.write(header + "\n" + seq + "\n+\n" + "I" * len(seq) + "\n")


def test_load_umi_map_parses_rx(tmp_path):
    fq = tmp_path / "r1.fastq"
    _write_fastq(fq, [
        ("@READ1\tRN:i:0\tRX:Z:ACGTACGTACGT", "GGGGCCCC"),
        ("@READ2\tRN:i:1\tRX:Z:TTTTGGGGCCCC", "GGGGCCCC"),
    ])
    m = _load_fastq_umi_map(str(fq))
    assert m == {"READ1": "ACGTACGTACGT", "READ2": "TTTTGGGGCCCC"}


def test_load_umi_map_gzip(tmp_path):
    fq = tmp_path / "r1.fastq.gz"
    _write_fastq(fq, [("@READ1\tRN:i:0\tRX:Z:ACGTACGTACGT", "GGGG")], gz=True)
    assert _load_fastq_umi_map(str(fq)) == {"READ1": "ACGTACGTACGT"}


def test_load_umi_map_empty_when_no_rx(tmp_path):
    """No RX tags -> {} -> downstream injection is a no-op (non-UMI runs safe)."""
    fq = tmp_path / "r1.fastq"
    _write_fastq(fq, [("@READ1\tRN:i:0", "GGGGCCCC")])
    assert _load_fastq_umi_map(str(fq)) == {}


def test_load_umi_map_strips_at_and_uses_bare_qname(tmp_path):
    fq = tmp_path / "r1.fastq"
    _write_fastq(fq, [("@A00197:1:1:1:1\tRN:i:0\tRX:Z:ACGTACGTACGT", "GGGG")])
    assert _load_fastq_umi_map(str(fq)) == {"A00197:1:1:1:1": "ACGTACGTACGT"}


# ---------------------------------------------------------------------------
# _inject_rn_into_bam with UMI (the money path)
# ---------------------------------------------------------------------------

def test_inject_rx_tags_both_mates(tmp_path):
    """Both mates share one QNAME and must BOTH receive the fragment UMI."""
    hdr = _header()
    bam = tmp_path / "aln.bam"
    _write_bam(bam, [
        _read(hdr, "PAIR1", pos=100, flag=0x1 | 0x40),   # R1
        _read(hdr, "PAIR1", pos=300, flag=0x1 | 0x80),   # R2, same QNAME
    ], hdr)
    n = _inject_rn_into_bam(str(bam), {"PAIR1": 0}, {"PAIR1": "ACGTACGTACGT"})
    tags = _read_tags(bam)
    assert len(tags["PAIR1"]) == 2
    for rec in tags["PAIR1"]:
        assert rec.get("RX") == "ACGTACGTACGT"
        assert rec.get("RN") == 0
    assert n == 2  # RN tagged count


def test_inject_rn_only_when_umi_map_none(tmp_path):
    """Passing no UMI map preserves original RN-only behaviour exactly."""
    hdr = _header()
    bam = tmp_path / "aln.bam"
    _write_bam(bam, [_read(hdr, "R1")], hdr)
    _inject_rn_into_bam(str(bam), {"R1": 5})
    tags = _read_tags(bam)["R1"][0]
    assert tags.get("RN") == 5 and "RX" not in tags


def test_inject_rx_without_rn_map(tmp_path):
    """UMI-only injection works even with an empty RN map (RX is independent)."""
    hdr = _header()
    bam = tmp_path / "aln.bam"
    _write_bam(bam, [_read(hdr, "R1")], hdr)
    _inject_rn_into_bam(str(bam), {}, {"R1": "TTTTGGGGCCCC"})
    tags = _read_tags(bam)["R1"][0]
    assert tags.get("RX") == "TTTTGGGGCCCC" and "RN" not in tags


def test_inject_skips_records_absent_from_umi_map(tmp_path):
    hdr = _header()
    bam = tmp_path / "aln.bam"
    _write_bam(bam, [_read(hdr, "HAS"), _read(hdr, "MISSING")], hdr)
    _inject_rn_into_bam(str(bam), {"HAS": 0, "MISSING": 1}, {"HAS": "ACGTACGTACGT"})
    tags = _read_tags(bam)
    assert tags["HAS"][0].get("RX") == "ACGTACGTACGT"
    assert "RX" not in tags["MISSING"][0]
    assert tags["MISSING"][0].get("RN") == 1  # RN still applied


def test_inject_empty_maps_is_noop(tmp_path):
    hdr = _header()
    bam = tmp_path / "aln.bam"
    _write_bam(bam, [_read(hdr, "R1")], hdr)
    assert _inject_rn_into_bam(str(bam), {}, {}) == 0
    assert "RX" not in _read_tags(bam)["R1"][0]


def test_end_to_end_fastq_to_bam_rx(tmp_path):
    """The full carry-through: RX in FASTQ header -> loaded -> injected into BAM."""
    fq = tmp_path / "r1.fastq"
    _write_fastq(fq, [
        ("@FRAG1\tRN:i:0\tRX:Z:ACGTACGTACGT", "GGGGCCCCAAAA"),
        ("@FRAG2\tRN:i:1\tRX:Z:TTTTGGGGCCCC", "GGGGCCCCAAAA"),
    ])
    hdr = _header()
    bam = tmp_path / "aln.bam"
    _write_bam(bam, [_read(hdr, "FRAG1"), _read(hdr, "FRAG2")], hdr)

    umi_map = _load_fastq_umi_map(str(fq))
    _inject_rn_into_bam(str(bam), {"FRAG1": 0, "FRAG2": 1}, umi_map)

    tags = _read_tags(bam)
    assert tags["FRAG1"][0].get("RX") == "ACGTACGTACGT"
    assert tags["FRAG2"][0].get("RX") == "TTTTGGGGCCCC"
