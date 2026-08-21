#!/usr/bin/env python3
"""Unit tests for the PARE (5'P degradome) arm (:mod:`rectify.core.pare`).

All M1-runnable: synthetic pysam reads + a tiny indexed FASTA, no bbmap/bbduk.
Covers the SENSE geometry (the mirror of the NET-seq antisense arm), the
strand-aware 5'P extraction + mandatory clip gate, the Tier-1 single-pass
pileup outputs (CPA pileup, 5'P pileup, bedgraph convention, read-length
census), the unmapped 3' poly-A rescue trim, and the rescue quantify's
single-molecule (5'P, pA-junction) pairing. The bbmap-dependent alignment
steps are exercised only on the cluster integration smoke.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pysam
import pytest

from rectify.core.netseq_cpa import pileup as nq_pileup
from rectify.core.netseq_cpa.pileup import LEFT, RIGHT, PILEUP_COLUMNS
from rectify.core.pare import pileup, rescue

CHROM = "chrI"
# pos: 0..999 X | 1000-1004 ACGTC | 1005-1009 GGGGG | 1010.. X   (+ gene case)
REF_PLUS = "X" * 1000 + "ACGTC" + "GGGGG" + "X" * 100
# pos: 0..999 X | 1000-1004 GGGGG | 1005-1009 ACGTC | 1010.. X   (- gene case)
REF_MINUS = "X" * 1000 + "GGGGG" + "ACGTC" + "X" * 100


def _header(chrom=CHROM, ln=2000):
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": ln}]}
    )


def _read(hdr, *, name="r", start=1000, seq="", cigar=((0, 0),),
          is_reverse=False, unmapped=False, mapq=60, qual=None):
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.is_unmapped = unmapped
    if unmapped:
        r.reference_id = -1
        r.reference_start = -1
    else:
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


def _write_genome(tmp_path: Path, seq: str, chrom=CHROM) -> Path:
    fa = tmp_path / "genome.fa"
    fa.write_text(f">{chrom}\n{seq}\n")
    pysam.faidx(str(fa))
    return fa


def _write_bam(tmp_path: Path, reads, name="in.bam") -> Path:
    bam = tmp_path / name
    hdr = _header()
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as out:
        for r in reads:
            out.write(r)
    return bam


class _StubWriter:
    """Collects parquet-writer add() kwargs without needing pyarrow."""

    def __init__(self):
        self.rows = []

    def add(self, **kw):
        self.rows.append(kw)


# --------------------------------------------------------------------------- #
# Geometry: PARE is SENSE — the exact mirror of the NET-seq antisense mapping
# --------------------------------------------------------------------------- #
def test_geometry_is_sense_mirror_of_netseq():
    hdr = _header()
    fwd = _read(hdr, seq="ACGTC", cigar=((0, 5),))
    rev = _read(hdr, seq="ACGTC", cigar=((0, 5),), is_reverse=True)
    assert pileup.geometry_for_read(fwd) == ("+", RIGHT, "A")
    assert pileup.geometry_for_read(rev) == ("-", LEFT, "T")
    # mirror check: NET-seq assigns the same (side, stop) to the OPPOSITE flag
    assert nq_pileup.geometry_for_read(rev)[1:] == pileup.geometry_for_read(fwd)[1:]
    assert nq_pileup.geometry_for_read(fwd)[1:] == pileup.geometry_for_read(rev)[1:]


def test_five_prime_of_forward_and_reverse():
    hdr = _header()
    fwd = _read(hdr, start=1000, seq="ACGTC", cigar=((0, 5),))
    assert pileup.five_prime_of(fwd) == (1000, 0)
    fwd_clip = _read(hdr, start=1000, seq="GGACGTC", cigar=((4, 2), (0, 5)))
    assert pileup.five_prime_of(fwd_clip) == (1000, 2)
    rev = _read(hdr, start=1005, seq="ACGTC", cigar=((0, 5),), is_reverse=True)
    assert pileup.five_prime_of(rev) == (1009, 0)
    rev_clip = _read(hdr, start=1005, seq="ACGTCGG", cigar=((0, 5), (4, 2)),
                     is_reverse=True)
    assert pileup.five_prime_of(rev_clip) == (1009, 2)
    # a clip on the RNA-3' side must NOT count as a 5' clip
    fwd_3clip = _read(hdr, start=1000, seq="ACGTCAAA", cigar=((0, 5), (4, 3)))
    assert pileup.five_prime_of(fwd_3clip) == (1000, 0)


# --------------------------------------------------------------------------- #
# Tier-1 single pass: CPA pileup + 5'P pileup + bedgraphs + length census
# --------------------------------------------------------------------------- #
def _read_tsv_gz(path: Path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh]


def test_pare_pileup_plus_strand_tail(tmp_path):
    """+ gene: FORWARD read, RIGHT soft-clip of 3 A's over non-A genome."""
    genome = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    tailed = _read(hdr, name="tailed", start=1000, seq="ACGTC" + "AAA",
                   cigar=((0, 5), (4, 3)))
    plain = _read(hdr, name="plain", start=1000, seq="ACGTC", cigar=((0, 5),))
    bam = _write_bam(tmp_path, [tailed, plain])
    w = _StubWriter()
    stats = pileup.pare_pileup(bam, genome, tmp_path / "out", "s1",
                               reads_parquet=w)
    assert stats["reads_used"] == 2
    assert stats["fivep_used"] == 2 and stats["fivep_clip_excluded"] == 0
    assert stats["oaNT_ge2_reads"] == 1 and stats["sum_oaNT"] == 3

    cpa = _read_tsv_gz(tmp_path / "out" / "cpa_pileup.tsv.gz")
    assert [r for r in cpa if r["pos"] == "1004" and r["strand"] == "+"], cpa
    row = next(r for r in cpa if r["pos"] == "1004")
    assert row["n_reads"] == "2" and row["n_oaNT_ge2"] == "1"
    assert list(row) == PILEUP_COLUMNS

    fp = _read_tsv_gz(tmp_path / "out" / "fivep_pileup.tsv.gz")
    assert len(fp) == 1 and fp[0]["pos"] == "1000" and fp[0]["n_reads"] == "2"

    bg = (tmp_path / "out" / "s1.5p_plus.bedgraph").read_text().splitlines()
    assert bg == ["chrI\t1000\t1001\t2"]
    assert (tmp_path / "out" / "s1.5p_minus.bedgraph").read_text() == ""

    lengths = _read_tsv_gz(tmp_path / "out" / "read_lengths.tsv.gz")
    q = {r["value"]: r["count"] for r in lengths if r["metric"] == "query_len"}
    assert q == {"8": "1", "5": "1"}
    oa = {r["value"]: r["count"] for r in lengths if r["metric"] == "oaNT"}
    assert oa == {"0": "1", "3": "1"}

    assert len(w.rows) == 2
    t = next(r for r in w.rows if r["read_id"] == "tailed")
    assert t["cpa_pos"] == 1004 and t["gene_strand"] == "+"
    assert t["oaNT_tail_len"] == 3 and t["five_p_pos"] == 1000
    assert t["five_p_clip"] == 0 and t["tier"] == "mapped"
    assert t["at_cpa"] is None            # no DRS map -> unknowable, not False


def test_pare_pileup_minus_strand_tail(tmp_path):
    """- gene: REVERSE read, LEFT soft-clip of 3 T's (= RNA poly-A) over non-T genome."""
    genome = _write_genome(tmp_path, REF_MINUS)
    hdr = _header()
    tailed = _read(hdr, name="tailed", start=1005, seq="TTT" + "ACGTC",
                   cigar=((4, 3), (0, 5)), is_reverse=True)
    bam = _write_bam(tmp_path, [tailed])
    w = _StubWriter()
    stats = pileup.pare_pileup(bam, genome, tmp_path / "out", "s1",
                               reads_parquet=w)
    assert stats["reads_used"] == 1 and stats["sum_oaNT"] == 3
    cpa = _read_tsv_gz(tmp_path / "out" / "cpa_pileup.tsv.gz")
    assert cpa[0]["pos"] == "1005" and cpa[0]["strand"] == "-"
    fp = _read_tsv_gz(tmp_path / "out" / "fivep_pileup.tsv.gz")
    assert fp[0]["pos"] == "1009" and fp[0]["strand"] == "-"
    bg = (tmp_path / "out" / "s1.5p_minus.bedgraph").read_text().splitlines()
    assert bg == ["chrI\t1009\t1010\t1"]
    assert w.rows[0]["five_p_pos"] == 1009 and w.rows[0]["cpa_pos"] == 1005


def test_pare_pileup_templated_a_tract_absorbed(tmp_path):
    """Read A's aligned OVER a genomic A-tract are templated -> oaNT excludes them."""
    ref = "X" * 1000 + "ACGTC" + "AAAAA" + "X" * 100   # genomic A-tract at 1005-1009
    genome = _write_genome(tmp_path, ref)
    hdr = _header()
    # aligner absorbed 3 tail A's as matches over the genomic tract
    absorbed = _read(hdr, name="absorbed", start=1000, seq="ACGTC" + "AAA",
                     cigar=((0, 8),))
    bam = _write_bam(tmp_path, [absorbed])
    stats = pileup.pare_pileup(bam, genome, tmp_path / "out", "s1")
    assert stats["sum_oaNT"] == 0          # all-templated: no CPA evidence
    cpa = _read_tsv_gz(tmp_path / "out" / "cpa_pileup.tsv.gz")
    # walkback still walks the tract back to the last non-A anchored base
    assert cpa[0]["pos"] == "1004" and cpa[0]["n_oaNT_ge1"] == "0"


def test_pare_pileup_fivep_clip_gate(tmp_path):
    """A 5'-clipped read is excluded from the 5'P pileup but kept for CPA."""
    genome = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    clipped = _read(hdr, name="clipped", start=1000, seq="GG" + "ACGTC",
                    cigar=((4, 2), (0, 5)))
    clean = _read(hdr, name="clean", start=1000, seq="ACGTC", cigar=((0, 5),))
    bam = _write_bam(tmp_path, [clipped, clean])
    stats = pileup.pare_pileup(bam, genome, tmp_path / "out", "s1")
    assert stats["reads_used"] == 2
    assert stats["fivep_used"] == 1 and stats["fivep_clip_excluded"] == 1
    assert stats["fivep_clip_fraction"] == pytest.approx(0.5)
    fp = _read_tsv_gz(tmp_path / "out" / "fivep_pileup.tsv.gz")
    assert len(fp) == 1 and fp[0]["n_reads"] == "1"
    cpa = _read_tsv_gz(tmp_path / "out" / "cpa_pileup.tsv.gz")
    assert cpa[0]["n_reads"] == "2"        # CPA side keeps both


# --------------------------------------------------------------------------- #
# Rescue: unmapped 3' poly-A trim + (5'P, pA-junction) quantify
# --------------------------------------------------------------------------- #
def test_trim_unmapped_polya(tmp_path):
    hdr = _header()
    anchor = "ACGTACGTACGTACGTACG"                     # 19 nt
    reads = [
        _read(hdr, name="keep", unmapped=True, seq=anchor + "A" * 10,
              qual="I" * 29),
        _read(hdr, name="short_tail", unmapped=True, seq=anchor + "A" * 5),
        _read(hdr, name="short_anchor", unmapped=True, seq="ACGTACGTAC" + "A" * 10),
        _read(hdr, name="trailing_n", unmapped=True, seq=anchor + "A" * 8 + "NN"),
        _read(hdr, name="mapped", seq=anchor, cigar=((0, 19),)),
    ]
    bam = _write_bam(tmp_path, reads)
    out_fq = tmp_path / "anchors.fq.gz"
    stats = rescue.trim_unmapped_polya(bam, out_fq, min_pa=8, min_anchor=15)
    assert stats == {"unmapped": 4, "rescued": 2}
    with gzip.open(out_fq, "rt") as fh:
        recs = fh.read().splitlines()
    names = [ln[1:] for ln in recs[0::4]]
    seqs = recs[1::4]
    assert names == ["keep_pa10", "trailing_n_pa8"]
    assert seqs == [anchor, anchor]


def test_quantify_rescue_pairs_fivep_with_cpa(tmp_path):
    hdr = _header()
    fwd = _read(hdr, name="f_pa12", start=1000, seq="A" * 19, cigar=((0, 19),),
                mapq=30)
    rev = _read(hdr, name="r_pa9", start=1005, seq="A" * 5, cigar=((0, 5),),
                is_reverse=True, mapq=30)
    low = _read(hdr, name="low_pa9", start=1000, seq="A" * 5, cigar=((0, 5),),
                mapq=1)
    bam = _write_bam(tmp_path, [fwd, rev, low])
    drs = tmp_path / "clusters.tsv"
    drs.write_text("chrom\tstrand\tmodal_position\nchrI\t+\t1018\n")
    w = _StubWriter()
    stats = rescue.quantify_rescue(bam, drs, min_mapq=3, reads_parquet=w)
    assert stats["mapped"] == 2            # low-MAPQ anchor dropped
    assert stats["at_cpa"] == 1 and stats["frac_at_cpa"] == pytest.approx(0.5)
    rows = {r["read_id"]: r for r in w.rows}
    assert rows["f"]["cpa_pos"] == 1018 and rows["f"]["five_p_pos"] == 1000
    assert rows["f"]["gene_strand"] == "+" and rows["f"]["at_cpa"] is True
    assert rows["r"]["cpa_pos"] == 1005 and rows["r"]["five_p_pos"] == 1009
    assert rows["r"]["gene_strand"] == "-" and rows["r"]["at_cpa"] is False
    assert all(r["tier"] == "rescued" and r["five_p_clip"] == 0
               for r in w.rows)


def test_quantify_rescue_without_drs_map(tmp_path):
    """No DRS map -> at_cpa is None per row, counts stay 0 (unknowable)."""
    hdr = _header()
    fwd = _read(hdr, name="f_pa12", start=1000, seq="A" * 19, cigar=((0, 19),),
                mapq=30)
    bam = _write_bam(tmp_path, [fwd])
    w = _StubWriter()
    stats = rescue.quantify_rescue(bam, None, reads_parquet=w)
    assert stats["mapped"] == 1 and stats["at_cpa"] == 0
    assert w.rows[0]["at_cpa"] is None


# --------------------------------------------------------------------------- #
# Parquet sidecar (schema round-trip; skipped without pyarrow)
# --------------------------------------------------------------------------- #
def test_pare_parquet_roundtrip(tmp_path):
    pa = pytest.importorskip("pyarrow")
    pq = pytest.importorskip("pyarrow.parquet")
    from rectify.core.pare.reads_parquet import COLUMNS, PareReadWriter
    path = tmp_path / "reads.parquet"
    w = PareReadWriter(path)
    w.add(read_id="r1", chrom="chrI", cpa_pos=1004, gene_strand="+",
          oaNT_tail_len=3, at_cpa=None, tier="mapped", mapq=60,
          five_p_pos=1000, five_p_clip=0)
    w.add(read_id="r2", chrom="chrI", cpa_pos=1018, gene_strand="-",
          oaNT_tail_len=12, at_cpa=True, tier="rescued", mapq=30,
          five_p_pos=1036, five_p_clip=0)
    assert w.close() == 2
    t = pq.read_table(path)
    assert t.column_names == COLUMNS
    assert t.num_rows == 2
    assert t.column("five_p_pos").to_pylist() == [1000, 1036]
    assert t.column("at_cpa").to_pylist() == [None, True]
