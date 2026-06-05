#!/usr/bin/env python3
"""Unit tests for the ``netseq-cpa`` per-read parquet sidecar
(:mod:`rectify.core.netseq_cpa.reads_parquet`).

Needs pyarrow to read the output back, so the whole module ``importorskip``s it
(absent on the M1 base env; present in the Sherlock rectify env, where this runs).
Covers: (1) the additive per-read emit preserves the Tier-1 TSV byte-for-byte;
(2) the mapped-row count equals ``reads_used``; (3) ``at_cpa`` is True/False with a
CPA set and **null** without one; (4) the writer no-ops when disabled.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pysam
import pytest

pa = pytest.importorskip("pyarrow")
import pyarrow.parquet as pq  # noqa: E402

from rectify.core.netseq_cpa import pileup, rescue  # noqa: E402
from rectify.core.netseq_cpa.reads_parquet import COLUMNS, ParquetReadWriter  # noqa: E402

CHROM = "chrI"
REF_PLUS = "X" * 1000 + "ACGTC" + "GGGGG" + "X" * 100


def _header(chrom=CHROM, ln=2000):
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": ln}]}
    )


def _read(hdr, *, name="r", start=1000, seq="", cigar=((0, 0),),
          is_reverse=False, mapq=60, qual=None):
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.is_unmapped = False
    r.reference_name = CHROM
    r.reference_start = start
    r.cigartuples = list(cigar)
    r.mapping_quality = mapq
    r.is_reverse = is_reverse
    r.is_secondary = False
    r.is_supplementary = False
    r.query_sequence = seq
    if qual is not None:
        r.query_qualities = pysam.qualitystring_to_array(qual)
    return r


def _write_genome(tmp_path: Path, seq: str) -> Path:
    fa = tmp_path / "genome.fa"
    fa.write_text(f">{CHROM}\n{seq}\n")
    pysam.faidx(str(fa))
    return fa


def _write_bam(tmp_path: Path, reads, name="in.bam") -> Path:
    bam = tmp_path / name
    with pysam.AlignmentFile(str(bam), "wb", header=_header()) as out:
        for r in reads:
            out.write(r)
    return bam


def _two_plus_reads(hdr):
    # two + reads: 5-base anchor + 3 non-templated clipped A's -> CPA 1004, oaNT 3
    return [
        _read(hdr, name="a", start=1000, seq="ACGTC" + "AAA",
              cigar=((0, 5), (4, 3)), is_reverse=True, mapq=44, qual="I" * 8),
        _read(hdr, name="b", start=1000, seq="ACGTC" + "AAA",
              cigar=((0, 5), (4, 3)), is_reverse=True, mapq=60, qual="I" * 8),
    ]


def test_parquet_emit_preserves_pileup_tsv(tmp_path):
    """The TSV is byte-identical with and without the parquet writer (parity)."""
    fa = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    bam = _write_bam(tmp_path, _two_plus_reads(hdr))

    out_no = tmp_path / "no.tsv.gz"
    pileup.walkback_pileup(bam, fa, out_no, "smoke")
    out_yes = tmp_path / "yes.tsv.gz"
    w = ParquetReadWriter(tmp_path / "reads.parquet")
    pileup.walkback_pileup(bam, fa, out_yes, "smoke", reads_parquet=w,
                           cpa_set={(CHROM, "+", 1004)})
    w.close()
    assert gzip.open(out_no, "rb").read() == gzip.open(out_yes, "rb").read()


def test_mapped_rows_equal_reads_used_and_at_cpa_true(tmp_path):
    fa = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    bam = _write_bam(tmp_path, _two_plus_reads(hdr))
    pqpath = tmp_path / "reads.parquet"
    w = ParquetReadWriter(pqpath)
    stats = pileup.walkback_pileup(bam, fa, tmp_path / "p.tsv.gz", "smoke",
                                   reads_parquet=w, cpa_set={(CHROM, "+", 1004)})
    n_rows = w.close()

    assert n_rows == stats["reads_used"] == 2
    t = pq.read_table(pqpath)
    assert t.column_names == COLUMNS
    d = t.to_pydict()
    assert sorted(d["read_id"]) == ["a", "b"]
    assert set(d["tier"]) == {"mapped"}
    assert d["chrom"] == [CHROM, CHROM]
    assert d["cpa_pos"] == [1004, 1004]
    assert d["gene_strand"] == ["+", "+"]
    assert d["oaNT_tail_len"] == [3, 3]
    assert d["at_cpa"] == [True, True]          # corr 1004 in cpa_set
    assert sorted(d["mapq"]) == [44, 60]


def test_at_cpa_false_off_site(tmp_path):
    fa = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    bam = _write_bam(tmp_path, _two_plus_reads(hdr))
    pqpath = tmp_path / "reads.parquet"
    w = ParquetReadWriter(pqpath)
    # CPA set at a DIFFERENT coordinate -> at_cpa False
    pileup.walkback_pileup(bam, fa, tmp_path / "p.tsv.gz", "smoke",
                           reads_parquet=w, cpa_set={(CHROM, "+", 500)})
    w.close()
    assert pq.read_table(pqpath).to_pydict()["at_cpa"] == [False, False]


def test_at_cpa_null_without_cpa_set(tmp_path):
    """No DRS map supplied -> at_cpa is null (unknowable != not-at-CPA)."""
    fa = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    bam = _write_bam(tmp_path, _two_plus_reads(hdr))
    pqpath = tmp_path / "reads.parquet"
    w = ParquetReadWriter(pqpath)
    pileup.walkback_pileup(bam, fa, tmp_path / "p.tsv.gz", "smoke",
                           reads_parquet=w, cpa_set=None)
    w.close()
    assert pq.read_table(pqpath).to_pydict()["at_cpa"] == [None, None]


def test_batch_flush_persists_all_rows(tmp_path):
    """batch_size smaller than the read count -> multiple flushes, all rows kept."""
    fa = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    reads = []
    for i in range(7):
        reads.append(_read(hdr, name=f"r{i}", start=1000, seq="ACGTC" + "AAA",
                           cigar=((0, 5), (4, 3)), is_reverse=True, qual="I" * 8))
    bam = _write_bam(tmp_path, reads)
    pqpath = tmp_path / "reads.parquet"
    w = ParquetReadWriter(pqpath, batch_size=2)   # forces 4 flushes
    stats = pileup.walkback_pileup(bam, fa, tmp_path / "p.tsv.gz", "smoke",
                                   reads_parquet=w, cpa_set={(CHROM, "+", 1004)})
    n_rows = w.close()
    assert n_rows == stats["reads_used"] == 7
    assert pq.read_table(pqpath).num_rows == 7


def test_disabled_writer_is_noop(tmp_path):
    """A disabled writer accepts adds, writes nothing, never raises."""
    pqpath = tmp_path / "reads.parquet"
    w = ParquetReadWriter(pqpath)
    w.disabled = True
    w.add(read_id="x", chrom=CHROM, cpa_pos=1, gene_strand="+",
          oaNT_tail_len=3, at_cpa=True, tier="mapped", mapq=60)
    assert w.close() == 0
    assert not pqpath.exists()
