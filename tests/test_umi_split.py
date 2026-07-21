"""Tests for UMI extraction in split_fastq_paired (rectify split --umi).

Covers both entry modes on tiny in-memory FASTQ pairs (no aligners):
  * read-id  -- pre-extracted data (umi_tools convention). THE critical path
                (Angel's RES CORALL libraries arrive this way).
  * r1-start -- raw CORALL (slice N12 off Read 1).

The invariant that matters most: both mates of a pair carry the SAME RX:Z UMI,
and the bare QNAME is identical across mates (so paired aligners still pair them).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from rectify.core.commands.split_command import split_fastq_paired


def _write_fastq(path, records):
    with open(str(path), "w") as fh:
        for header, seq in records:
            fh.write(header + "\n" + seq + "\n+\n" + "I" * len(seq) + "\n")


def _chunk_records(gz_path):
    """Yield (header, seq) for each record in a gzipped chunk."""
    with gzip.open(str(gz_path), "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline(); fh.readline(); fh.readline()
            yield h.rstrip("\n"), s.rstrip("\n")


def _bare_qname(header):
    return header.split()[0].lstrip("@").split("\t")[0]


def _rx(header):
    for tok in header.replace("\t", " ").split():
        if tok.startswith("RX:Z:"):
            return tok[5:]
    return None


# ---------------------------------------------------------------------------
# read-id mode (the critical path)
# ---------------------------------------------------------------------------

def test_split_umi_read_id_tags_both_mates(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    umis = ["ACGTACGTACGT", "TTTTGGGGCCCC"]
    _write_fastq(r1, [(f"@READ{i}_{u} 1:N:0:AT", "G" * 50) for i, u in enumerate(umis)])
    _write_fastq(r2, [(f"@READ{i}_{u} 2:N:0:AT", "C" * 50) for i, u in enumerate(umis)])

    r1p, r2p = split_fastq_paired(r1, r2, tmp_path / "out", n_chunks=1, prefix="s",
                                  umi=True, umi_location="read-id", umi_length=12)
    r1recs = list(_chunk_records(r1p[0]))
    r2recs = list(_chunk_records(r2p[0]))
    assert [_rx(h) for h, _ in r1recs] == umis
    assert [_rx(h) for h, _ in r2recs] == umis          # R2 gets the fragment UMI too
    # QNAME identical across mates and stripped of the UMI
    for (h1, _), (h2, _) in zip(r1recs, r2recs):
        assert _bare_qname(h1) == _bare_qname(h2)
        assert "ACGT" not in _bare_qname(h1) or _bare_qname(h1).startswith("READ")
        assert _bare_qname(h1).startswith("READ")


def test_split_umi_read_id_strips_umi_from_qname(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    _write_fastq(r1, [("@FRAG_ACGTACGTACGT 1:N:0:AT", "G" * 50)])
    _write_fastq(r2, [("@FRAG_ACGTACGTACGT 2:N:0:AT", "C" * 50)])
    r1p, _ = split_fastq_paired(r1, r2, tmp_path / "out", n_chunks=1, prefix="s",
                                umi=True, umi_location="read-id", umi_length=12)
    h, _ = next(_chunk_records(r1p[0]))
    assert _bare_qname(h) == "FRAG"          # the _UMI suffix is gone
    assert _rx(h) == "ACGTACGTACGT"


def test_split_umi_read_id_seq_untouched(tmp_path):
    """read-id mode must NOT alter the sequence (UMI was never in it)."""
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    _write_fastq(r1, [("@FRAG_ACGTACGTACGT 1:N:0:AT", "GATTACA" * 7)])
    _write_fastq(r2, [("@FRAG_ACGTACGTACGT 2:N:0:AT", "C" * 49)])
    r1p, _ = split_fastq_paired(r1, r2, tmp_path / "out", n_chunks=1, prefix="s",
                                umi=True, umi_location="read-id", umi_length=12)
    _, seq = next(_chunk_records(r1p[0]))
    assert seq == "GATTACA" * 7


# ---------------------------------------------------------------------------
# r1-start mode (raw CORALL)
# ---------------------------------------------------------------------------

def test_split_umi_r1_start_slices_umi_off_r1(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    umi, insert = "ACGTACGTACGT", "GGGGCCCCAAAATTTT" * 3
    _write_fastq(r1, [("@FRAG 1:N:0:AT", umi + insert)])
    _write_fastq(r2, [("@FRAG 2:N:0:AT", "C" * 48)])
    r1p, r2p = split_fastq_paired(r1, r2, tmp_path / "out", n_chunks=1, prefix="s",
                                  umi=True, umi_location="r1-start", umi_length=12)
    h1, seq1 = next(_chunk_records(r1p[0]))
    h2, seq2 = next(_chunk_records(r2p[0]))
    assert _rx(h1) == umi and _rx(h2) == umi
    assert seq1 == insert              # UMI sliced off R1
    assert seq2 == "C" * 48            # R2 untouched (CORALL has no R2 UMI)


def test_split_umi_r1_start_qual_sliced_with_seq(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    _write_fastq(r1, [("@FRAG 1:N:0:AT", "ACGTACGTACGT" + "G" * 40)])
    _write_fastq(r2, [("@FRAG 2:N:0:AT", "C" * 40)])
    r1p, _ = split_fastq_paired(r1, r2, tmp_path / "out", n_chunks=1, prefix="s",
                                umi=True, umi_location="r1-start", umi_length=12)
    with gzip.open(str(r1p[0]), "rt") as fh:
        fh.readline(); seq = fh.readline().rstrip(); fh.readline(); qual = fh.readline().rstrip()
    assert len(seq) == len(qual) == 40


def test_split_umi_r1_start_drops_too_short(tmp_path):
    """A read shorter than the UMI has no usable insert -> dropped, not emitted broken."""
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    _write_fastq(r1, [("@GOOD 1:N:0:AT", "ACGTACGTACGT" + "G" * 30),
                      ("@TOOSHORT 1:N:0:AT", "ACGT")])       # 4 < 12
    _write_fastq(r2, [("@GOOD 2:N:0:AT", "C" * 30),
                      ("@TOOSHORT 2:N:0:AT", "C" * 30)])
    r1p, _ = split_fastq_paired(r1, r2, tmp_path / "out", n_chunks=1, prefix="s",
                                umi=True, umi_location="r1-start", umi_length=12)
    names = [_bare_qname(h) for h, _ in _chunk_records(r1p[0])]
    assert names == ["GOOD"]           # TOOSHORT dropped


# ---------------------------------------------------------------------------
# umi disabled -> unchanged behaviour
# ---------------------------------------------------------------------------

def test_split_no_umi_writes_no_rx(tmp_path):
    r1, r2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    _write_fastq(r1, [("@READ0_ACGTACGTACGT 1:N:0:AT", "G" * 50)])
    _write_fastq(r2, [("@READ0_ACGTACGTACGT 2:N:0:AT", "C" * 50)])
    r1p, _ = split_fastq_paired(r1, r2, tmp_path / "out", n_chunks=1, prefix="s")  # umi=False
    h, _ = next(_chunk_records(r1p[0]))
    assert _rx(h) is None
    # Without --umi, the umi-tools suffix is preserved verbatim in the QNAME
    assert "ACGTACGTACGT" in h
