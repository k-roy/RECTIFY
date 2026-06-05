#!/usr/bin/env python3
"""Unit tests for the NET-seq CPA-intermediate capture arm
(:mod:`rectify.core.netseq_cpa`).

All M1-runnable: synthetic pysam reads + a tiny indexed FASTA, no bbmap/bowtie2.
Covers the Tier-1 genome-aware poly-A metric (and its parity with the canonical
``walkback_netseq`` wrapper), the unmapped poly-T rescue trim, mapped CPA
concordance, paired-end read-2 isolation, and the Tier-2 window-reference
builder. The bowtie2-dependent Tier-2 alignment is exercised only on the cluster
integration smoke (no bowtie2 on M1).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pysam
import pytest

from rectify.core.netseq_cpa import (
    NCBI_TO_CHROM,
    concordance,
    cpa_window_rescue,
    pe_reduce,
    pileup,
    rescue,
)
from rectify.core.netseq_cpa.pileup import PILEUP_COLUMNS

CHROM = "chrI"
# pos: 0..999 X | 1000-1004 ACGTC | 1005-1009 GGGGG | 1010.. X   (+ strand case)
REF_PLUS = "X" * 1000 + "ACGTC" + "GGGGG" + "X" * 100
# pos: 0..999 X | 1000-1004 GGGGG | 1005-1009 ACGTC | 1010.. X   (- strand case)
REF_MINUS = "X" * 1000 + "GGGGG" + "ACGTC" + "X" * 100


def _header(chrom=CHROM, ln=2000):
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": ln}]}
    )


def _read(hdr, *, name="r", start=1000, seq="", cigar=((0, 0),),
          is_reverse=False, unmapped=False, qual=None):
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
        r.mapping_quality = 60
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


# --------------------------------------------------------------------------- #
# Tier-1 metric: walkback counts + genome-aware soft-clip
# --------------------------------------------------------------------------- #
def test_softclip_run_nt_plus_strand():
    """+ gene: RIGHT clip of 3 A's over non-A genome -> oaNT soft-clip = 3."""
    hdr = _header()
    r = _read(hdr, start=1000, seq="ACGTC" + "AAA", cigar=((0, 5), (4, 3)), is_reverse=True)
    assert pileup.softclip_run_nt(r, REF_PLUS, pileup.RIGHT, "A") == 3
    _orig, corr, n_tmpl, n_nt = pileup.walkback_and_count(r, REF_PLUS, pileup.RIGHT, "A")
    assert corr == 1004 and n_nt == 0   # terminal aligned C is anchored; tail is the soft-clip


def test_softclip_run_nt_minus_strand():
    """- gene: LEFT clip of 3 T's over non-T genome -> oaNT soft-clip = 3."""
    hdr = _header()
    r = _read(hdr, start=1005, seq="TTT" + "ACGTC", cigar=((4, 3), (0, 5)), is_reverse=False)
    assert pileup.softclip_run_nt(r, REF_MINUS, pileup.LEFT, "T") == 3
    _orig, corr, n_tmpl, n_nt = pileup.walkback_and_count(r, REF_MINUS, pileup.LEFT, "T")
    assert corr == 1005 and n_nt == 0


def test_softclip_run_nt_subtracts_genomic_run():
    """Genome-aware: a clipped A-run sitting over a genomic A-tract is NOT poly-A."""
    hdr = _header()
    ref = "X" * 1000 + "ACGTC" + "AAAAA" + "X" * 100   # genomic A-run downstream
    r = _read(hdr, start=1000, seq="ACGTC" + "AAA", cigar=((0, 5), (4, 3)), is_reverse=True)
    # clip run 3 minus genomic A-run (>=3) -> 0 non-templated
    assert pileup.softclip_run_nt(r, ref, pileup.RIGHT, "A") == 0


def test_parity_with_walkback_netseq_both_strands():
    """walkback_and_count corrected CPA == canonical walkback_netseq corrected CPA."""
    hdr = _header()
    r_plus = _read(hdr, start=1000, seq="ACGTC" + "AAA", cigar=((0, 5), (4, 3)), is_reverse=True)
    assert pileup.corrected_position_matches_walkback_netseq(r_plus, REF_PLUS)
    r_minus = _read(hdr, start=1005, seq="TTT" + "ACGTC", cigar=((4, 3), (0, 5)), is_reverse=False)
    assert pileup.corrected_position_matches_walkback_netseq(r_minus, REF_MINUS)


# --------------------------------------------------------------------------- #
# Tier-1 end-to-end pileup
# --------------------------------------------------------------------------- #
def test_walkback_pileup_end_to_end(tmp_path):
    fa = _write_genome(tmp_path, REF_PLUS)
    hdr = _header()
    # two + reads, each a 5-base anchor + 3 non-templated clipped A's at CPA 1004
    reads = [
        _read(hdr, name="a", start=1000, seq="ACGTC" + "AAA",
              cigar=((0, 5), (4, 3)), is_reverse=True, qual="I" * 8),
        _read(hdr, name="b", start=1000, seq="ACGTC" + "AAA",
              cigar=((0, 5), (4, 3)), is_reverse=True, qual="I" * 8),
    ]
    bam = _write_bam(tmp_path, reads)
    out = tmp_path / "pileup.tsv.gz"
    stats = pileup.walkback_pileup(bam, fa, out, "smoke")
    assert stats["reads_used"] == 2
    assert stats["oaNT_ge2_reads"] == 2
    rows = list(gzip.open(out, "rt"))
    assert rows[0].rstrip("\n").split("\t") == PILEUP_COLUMNS
    data = [r.rstrip("\n").split("\t") for r in rows[1:]]
    assert len(data) == 1
    chrom, pos, strand = data[0][0], int(data[0][1]), data[0][2]
    assert (chrom, pos, strand) == (CHROM, 1004, "+")
    n_oaNT_ge2 = int(data[0][PILEUP_COLUMNS.index("n_oaNT_ge2")])
    assert n_oaNT_ge2 == 2


# --------------------------------------------------------------------------- #
# Rescue: unmapped poly-T trim
# --------------------------------------------------------------------------- #
def test_trim_unmapped_polyt(tmp_path):
    hdr = _header()
    anchor = "ACGT" * 6                       # 24 bp clean anchor (>= MIN_ANCHOR 20)
    seq = "N" + "T" * 10 + anchor             # leading N + 10 poly-T + anchor
    r = _read(hdr, name="u1", seq=seq, unmapped=True, qual="I" * len(seq))
    # an unmapped read with too-short anchor is dropped
    r2 = _read(hdr, name="u2", seq="N" + "T" * 10 + "ACGT", unmapped=True, qual="I" * 15)
    bam = _write_bam(tmp_path, [r, r2])
    fq = tmp_path / "anchors.fq.gz"
    stats = rescue.trim_unmapped_polyt(bam, fq, min_pt=8, min_anchor=20)
    assert stats["unmapped"] == 2 and stats["rescued"] == 1
    lines = gzip.open(fq, "rt").read().splitlines()
    assert lines[0] == "@u1_pa10"             # poly-A length encoded
    assert lines[1] == anchor                 # leading N + poly-T stripped


# --------------------------------------------------------------------------- #
# Concordance
# --------------------------------------------------------------------------- #
def test_mapped_cpa_concordance(tmp_path):
    pil = tmp_path / "pileup.tsv.gz"
    with gzip.open(pil, "wt") as fh:
        fh.write("\t".join(PILEUP_COLUMNS) + "\n")
        # row AT a CPA: 10 reads, 6 oaNT>=2, sum_oaNT 18, oaNT_ge1 8
        fh.write(f"{CHROM}\t1004\t+\t10\t8\t6\t4\t5\t20\t8\t6\t4\t18\tsmoke\n")
        # row OFF a CPA: 5 reads, 2 oaNT>=2
        fh.write(f"{CHROM}\t500\t+\t5\t3\t2\t1\t2\t6\t3\t2\t1\t6\tsmoke\n")
    clusters = tmp_path / "drs.tsv"
    clusters.write_text("chrom\tstrand\tmodal_position\n" f"{CHROM}\t+\t1004\n")
    res = concordance.mapped_cpa_concordance(pil, clusters, label="smoke", nuclear_bp=1_000_000)
    assert res["oaNT_reads"] == 8          # 6 + 2
    assert res["at_cpa"] == 6              # only the 1004 row is within +-5 of CPA 1004
    assert res["frac_at_cpa"] == pytest.approx(6 / 8)
    assert res["enrichment"] > 1


# --------------------------------------------------------------------------- #
# Paired-end read-2 isolation
# --------------------------------------------------------------------------- #
def test_reduce_bam_to_read2(tmp_path):
    hdr = _header()
    r1 = _read(hdr, name="p", start=1000, seq="ACGTC", cigar=((0, 5),))
    r1.is_paired = True
    r1.is_read1 = True
    r2 = _read(hdr, name="p", start=1100, seq="ACGTC", cigar=((0, 5),))
    r2.is_paired = True
    r2.is_read2 = True
    bam = _write_bam(tmp_path, [r1, r2])
    out = tmp_path / "r2.bam"
    stats = pe_reduce.reduce_bam_to_read2(bam, out)
    assert stats == {"seen": 2, "written": 1}
    recs = list(pysam.AlignmentFile(str(out), "rb").fetch(until_eof=True))
    assert len(recs) == 1
    assert not recs[0].is_paired and not recs[0].is_read2
    assert recs[0].reference_start == 1100


# --------------------------------------------------------------------------- #
# Tier-2 window-reference builder
# --------------------------------------------------------------------------- #
def test_build_cpa_window_reference(tmp_path):
    fa = _write_genome(tmp_path, REF_PLUS)
    clusters = tmp_path / "drs.tsv"
    clusters.write_text("chrom\tstrand\tmodal_position\n" f"{CHROM}\t+\t1004\n")
    out_fa = tmp_path / "windows.fa"
    n = cpa_window_rescue.build_cpa_window_reference(
        clusters, fa, out_fa, gene_bodyward=4, downstream=5
    )
    assert n == 1
    txt = out_fa.read_text().splitlines()
    # + strand window = [modal-4, modal+5+1) = [1000, 1010): "ACGTC" + "GGGGG"
    assert txt[0] == f">{CHROM}|+|1004|1000"
    assert txt[1] == "ACGTCGGGGG"


def test_ncbi_chrom_map_shape():
    assert NCBI_TO_CHROM["ref|NC_001133|"] == "chrI"
    assert NCBI_TO_CHROM["ref|NC_001148|"] == "chrXVI"
    assert NCBI_TO_CHROM["ref|NC_001224|"] == "chrMito"
