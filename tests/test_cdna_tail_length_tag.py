"""Tests for the ONT PCR-cDNA tail-length tag (`pl`) and the anchored orientation floor.

Both behaviours were added after planning/550 found that the trim stage REMOVES the poly(A)
tail before alignment, so `rectify correct` measured a read that by construction has none:
96.0% of sense and 95.0% of antisense reads arrived with polya_length == 0 despite the trim
stage having measured a median 22 nt / 18 nt tail for the same reads.

The tests below are deliberately structural -- they assert the tag is EMITTED, is
ORIENTATION-AWARE, and is CONSUMED -- because a real-data smoke can pass while silently
dropping the tail (that is exactly the failure being fixed).
"""
import gzip
from pathlib import Path

import pytest

from rectify.core.commands.cdna_trim_command import (
    _CRTA,
    _CRTA_RC,
    _MIN_ORIENTATION_TAIL,
    _MIN_ORIENTATION_TAIL_ANCHORED,
    _TAIL_LEN_TAG,
    trim_cdna_fastq_polya,
)
from rectify.core.correct.protocols.ont_cdna import (
    POLYA_SOURCE_TRIM,
    TAIL_LEN_TAG,
    trim_stage_tail_length,
)


BODY = "GTTATGTCCTGTCTTTGGTTCAGTTATTGAACCAATGTCACAGGCCTTCCTCCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTT"


def _write_fastq(path, records):
    with gzip.open(path, "wt") as fh:
        for name, seq in records:
            fh.write("@%s\n%s\n+\n%s\n" % (name, seq, "I" * len(seq)))


def _read_fastq_comments(path):
    """Return {read_id: {tag: value}} parsed from the output FASTQ header comments."""
    out = {}
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            fh.readline(); fh.readline(); fh.readline()
            parts = h[1:].rstrip("\n").split("\t")
            tags = {}
            for p in parts[1:]:
                bits = p.split(":")
                if len(bits) == 3:
                    tags[bits[0]] = bits[2]
            out[parts[0]] = tags
    return out


class _FakeRead:
    """Minimal pysam.AlignedSegment stand-in for tag lookup."""

    def __init__(self, tags):
        self._tags = tags

    def get_tag(self, tag):
        if tag not in self._tags:
            raise KeyError(tag)
        return self._tags[tag]


# --------------------------------------------------------------------------
# The `pl` tag is emitted, and is orientation-aware
# --------------------------------------------------------------------------

def test_sense_read_emits_3p_polya_length(tmp_path):
    """A sense read (3' poly-A + CRTA) reports the poly-A length on `pl`."""
    tail = 25
    seq = BODY + "A" * tail + _CRTA
    fq = tmp_path / "in.fastq.gz"
    _write_fastq(fq, [("sense_read", seq)])
    out = tmp_path / "out.fastq.gz"
    trim_cdna_fastq_polya(str(fq), str(out), str(tmp_path / "meta.tsv"),
                          trim_5p_polyt=True)
    tags = _read_fastq_comments(out)["sense_read"]
    assert tags["ro"] == "S"
    assert int(tags[_TAIL_LEN_TAG]) == tail


def test_antisense_read_emits_5p_polyt_length(tmp_path):
    """An antisense read reports the 5' poly-T length -- NOT 0, and NOT the 3' value.

    This is the half of the library that planning/550 found being recorded as
    tail-less; the tail is present, it is simply poly-T on the other strand.
    """
    tail = 18
    seq = _CRTA_RC + "T" * tail + BODY
    fq = tmp_path / "in.fastq.gz"
    _write_fastq(fq, [("antisense_read", seq)])
    out = tmp_path / "out.fastq.gz"
    trim_cdna_fastq_polya(str(fq), str(out), str(tmp_path / "meta.tsv"),
                          trim_5p_polyt=True)
    tags = _read_fastq_comments(out)["antisense_read"]
    assert tags["ro"] == "A"
    assert int(tags[_TAIL_LEN_TAG]) == tail


def test_unresolved_read_emits_zero(tmp_path):
    """A read with no tail at either end reports 0 -- explicitly, not by omission."""
    fq = tmp_path / "in.fastq.gz"
    _write_fastq(fq, [("bare_read", BODY)])
    out = tmp_path / "out.fastq.gz"
    trim_cdna_fastq_polya(str(fq), str(out), str(tmp_path / "meta.tsv"),
                          trim_5p_polyt=True)
    tags = _read_fastq_comments(out)["bare_read"]
    assert tags["ro"] == "U"
    assert int(tags[_TAIL_LEN_TAG]) == 0


def test_tag_present_on_every_read(tmp_path):
    """Every emitted read carries `pl` -- a missing tag would silently fall back."""
    fq = tmp_path / "in.fastq.gz"
    _write_fastq(fq, [
        ("r_sense", BODY + "A" * 22 + _CRTA),
        ("r_anti", _CRTA_RC + "T" * 15 + BODY),
        ("r_bare", BODY),
    ])
    out = tmp_path / "out.fastq.gz"
    trim_cdna_fastq_polya(str(fq), str(out), str(tmp_path / "meta.tsv"),
                          trim_5p_polyt=True)
    comments = _read_fastq_comments(out)
    assert len(comments) == 3
    for rid, tags in comments.items():
        assert _TAIL_LEN_TAG in tags, "read %s lost the tail-length tag" % rid


# --------------------------------------------------------------------------
# The anchored orientation floor
# --------------------------------------------------------------------------

def test_short_tail_with_adapter_is_called_sense(tmp_path):
    """An adapter-anchored SHORT tail is sense, not 'U'.

    Measured justification (planning/550): among 2,088,487 adapter-positive reads,
    ZERO have polya_3p_len == 0 -- splint ligation cannot attach the adapter
    anywhere but a genuine 3'-terminal A-tail. 7.6% of them carry a 1-9 nt tail and
    were being discarded by the flat >=10 gate.
    """
    short = 6
    assert short < _MIN_ORIENTATION_TAIL
    assert short >= _MIN_ORIENTATION_TAIL_ANCHORED
    seq = BODY + "A" * short + _CRTA
    fq = tmp_path / "in.fastq.gz"
    _write_fastq(fq, [("short_anchored", seq)])
    out = tmp_path / "out.fastq.gz"
    trim_cdna_fastq_polya(str(fq), str(out), str(tmp_path / "meta.tsv"),
                          trim_5p_polyt=True)
    tags = _read_fastq_comments(out)["short_anchored"]
    assert tags["ro"] == "S", "adapter-anchored short tail must not fall through to U"
    assert int(tags[_TAIL_LEN_TAG]) == short


def test_short_tail_without_adapter_stays_unresolved(tmp_path):
    """A short A-run with NO adapter stays 'U' -- genomic A-tracts must not become sense.

    This is the guard that keeps the relaxed floor from admitting internal-priming
    lookalikes; it is the reason the floor is conditional rather than global.
    """
    seq = BODY + "A" * 4
    fq = tmp_path / "in.fastq.gz"
    _write_fastq(fq, [("short_naked", seq)])
    out = tmp_path / "out.fastq.gz"
    trim_cdna_fastq_polya(str(fq), str(out), str(tmp_path / "meta.tsv"),
                          trim_5p_polyt=True)
    assert _read_fastq_comments(out)["short_naked"]["ro"] == "U"


# --------------------------------------------------------------------------
# The consumer side
# --------------------------------------------------------------------------

@pytest.mark.parametrize("value,expected", [(0, 0), (18, 18), ("25", 25)])
def test_trim_stage_tail_length_reads_the_tag(value, expected):
    assert trim_stage_tail_length(_FakeRead({TAIL_LEN_TAG: value})) == expected


def test_trim_stage_tail_length_absent_returns_none():
    """Absent tag must be None (fall back), NOT 0 (assert a tail-less read).

    Conflating the two is precisely how a BAM that never went through the trim
    stage would silently acquire zero-length tails.
    """
    assert trim_stage_tail_length(_FakeRead({})) is None


def test_trim_stage_tail_length_rejects_garbage():
    assert trim_stage_tail_length(_FakeRead({TAIL_LEN_TAG: "not_an_int"})) is None
    assert trim_stage_tail_length(_FakeRead({TAIL_LEN_TAG: -3})) is None


def test_polya_source_constant_is_distinct():
    """The source label must be distinguishable from the existing values."""
    assert POLYA_SOURCE_TRIM not in ("pt_tag", "model", "none")
