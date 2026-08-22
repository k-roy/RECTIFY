#!/usr/bin/env python3
"""Tests for `rectify cdna-analyze` (Stage 3 of the cDNA pipeline).

This module recomputes walkback / TSS-walk-forward / gene assignment /
isoform clustering / Type-1↔Type-2 pairing on the post-alignment coordinates
of the per-cluster consensus BAM produced by ``rectify align`` after
``rectify correct-cdna``.

Tests come in three layers:

  1. Import smoke test — guards against `cdna_analyze_command` itself failing
     to import (e.g. broken pysam tag schema after a refactor).
  2. Synthetic minimal smoke — builds a 3-cluster BAM and an in-memory
     synthetic reference + GFF, runs the full ``run(args)`` entry point,
     and checks the three TSVs land on disk with the expected gene
     assignments and isoform groupings.
  3. On-disk chrI smoke — opt-in, only fires when the M1 test fixtures are
     present (the full three-stage chain canary lives in
     ``test_cdna_chain_canary.py``).
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pysam
import pytest


# ---------------------------------------------------------------------------
# 1) Import smoke
# ---------------------------------------------------------------------------
def test_cdna_analyze_imports():
    """The module imports and exposes a `run` callable."""
    from rectify.core.commands import cdna_analyze_command

    assert callable(cdna_analyze_command.run)
    # The helper that builds a ReadInfo from a BAM record is also part of
    # the public-ish surface — assert it's there so refactors don't silently
    # remove it.
    assert callable(cdna_analyze_command._read_info_from_bam_record)


# ---------------------------------------------------------------------------
# Helpers for synthetic test fixtures
# ---------------------------------------------------------------------------

# A tiny synthetic chromosome:
#   - bases   0–199  : "X" (background)
#   - 200–229          gene_A body (+ strand)  "ACGTACGTAC" × 3
#   - 230–249          poly-A tail genomic context for gene_A reads
#   - 250–299          gap
#   - 300–329          gene_B body (- strand)
#   - rest             "X" padding to length 1000
_TINY_CHROM_LEN = 1000


def _build_tiny_reference(tmp_path: Path) -> Path:
    """Write a single-chromosome FASTA + index."""
    seq = list("X" * _TINY_CHROM_LEN)
    # gene_A region: simple ACGT pattern, then 20 A's for a tail context
    for i, b in enumerate("ACGTACGTACACGTACGTAC"):
        seq[200 + i] = b
    for i in range(20):
        seq[220 + i] = "A"
    # gene_B region: same pattern but as a different stretch
    for i, b in enumerate("TGCATGCATGCATGCATGCA"):
        seq[300 + i] = b

    fa = tmp_path / "tiny.fa"
    fa.write_text(f">chrT\n{''.join(seq)}\n")
    pysam.faidx(str(fa))
    return fa


def _build_tiny_gff(tmp_path: Path) -> Path:
    """Write a tiny SGD-style GFF with two genes.

    Use bare ID values (no "gene:" namespace prefix) to match SGD's emission
    style, since `_parse_gff_gene_name` returns the ID verbatim.
    """
    gff = tmp_path / "tiny.gff"
    gff.write_text(
        "##gff-version 3\n"
        # 1-based, inclusive — gene_A on + strand spans 201..240
        "chrT\tsynthetic\tgene\t201\t240\t.\t+\t.\tID=gene_A;Name=gene_A\n"
        # gene_B on - strand spans 301..320
        "chrT\tsynthetic\tgene\t301\t320\t.\t-\t.\tID=gene_B;Name=gene_B\n"
    )
    return gff


def _build_synthetic_consensus_bam(tmp_path: Path) -> Path:
    """Build a 3-record consensus BAM resembling `rectify align` output.

    Records:
      cl0  gene_A sense   (orient=fwd, aln 200..230, ends in A-tract)
      cl1  gene_A sense   (orient=fwd, aln 200..230, ends in A-tract; same
                           gene + termini → same isoform as cl0)
      cl2  gene_B anti    (orient=fwd, aln 300..320 — overlaps gene_B on - strand;
                           orient=fwd → implied mRNA is + → classified antisense)
    """
    bam_path = tmp_path / "consensus.bam"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chrT", "LN": _TINY_CHROM_LEN}],
    }
    # Use a unique UMI per Type-1 record so the canonical-UMI choice is deterministic.
    umis = [
        "ACGTACGTACGTACGTACGTACGTACG",
        "TGCATGCATGCATGCATGCATGCATGC",
        "GGGGGGGGGGGGGGGGGGGGGGGGGGG",  # cl2 placeholder
    ]
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # cl0
        r = pysam.AlignedSegment(bam.header)
        r.query_name = "cluster_0"
        # 10 bp of mRNA body + 20 bp poly-A tail aligned over genomic A's
        r.query_sequence = "ACGTACGTAC" + "A" * 20
        r.flag = 0
        r.reference_id = 0
        r.reference_start = 200
        r.mapping_quality = 60
        r.cigartuples = [(0, 30)]
        r.set_tag("XU", umis[0], "Z")
        r.set_tag("XO", "fwd", "Z")
        r.set_tag("XT", 1, "i")
        r.set_tag("XY", "umi_captured_fwd", "Z")
        r.set_tag("XC", 5, "i")
        r.set_tag("XF", 2, "i")
        r.set_tag("XR", "rd_a1,rd_a2,rd_a3,rd_a4,rd_a5", "Z")
        bam.write(r)
        # cl1
        r = pysam.AlignedSegment(bam.header)
        r.query_name = "cluster_1"
        r.query_sequence = "ACGTACGTAC" + "A" * 20
        r.flag = 0
        r.reference_id = 0
        r.reference_start = 200
        r.mapping_quality = 60
        r.cigartuples = [(0, 30)]
        r.set_tag("XU", umis[1], "Z")
        r.set_tag("XO", "fwd", "Z")
        r.set_tag("XT", 1, "i")
        r.set_tag("XY", "umi_captured_fwd", "Z")
        r.set_tag("XC", 3, "i")
        r.set_tag("XF", 2, "i")
        r.set_tag("XR", "rd_b1,rd_b2,rd_b3", "Z")
        bam.write(r)
        # cl2 — gene_B (- strand) overlap; orient=fwd → +mRNA → antisense
        r = pysam.AlignedSegment(bam.header)
        r.query_name = "cluster_2"
        r.query_sequence = "TGCATGCATGCATGCATGCA"
        r.flag = 0
        r.reference_id = 0
        r.reference_start = 300
        r.mapping_quality = 60
        r.cigartuples = [(0, 20)]
        r.set_tag("XU", umis[2], "Z")
        r.set_tag("XO", "fwd", "Z")
        r.set_tag("XT", 1, "i")
        r.set_tag("XY", "umi_captured_fwd", "Z")
        r.set_tag("XC", 2, "i")
        r.set_tag("XF", 1, "i")
        r.set_tag("XR", "rd_c1,rd_c2", "Z")
        bam.write(r)
    pysam.index(str(bam_path))
    return bam_path


def _make_args(bam: Path, out: Path, gff: Path, ref: Path) -> argparse.Namespace:
    """Build the minimum argparse.Namespace that `cdna_analyze_command.run` reads."""
    return argparse.Namespace(
        bam=bam,
        out=out,
        gff=gff,
        reference=ref,
        # Defaults applied by run() if missing; we set them explicitly for
        # clarity and to make the test resilient to default-value changes.
        min_gene_frac=0.3,
        min_read_frac=0.8,
        isoform_tol_5=5,
        isoform_tol_3=5,
        t1t2_tol_5=5,
        t1t2_tol_3=5,
    )


# ---------------------------------------------------------------------------
# 2) Synthetic minimal smoke
# ---------------------------------------------------------------------------
def test_cdna_analyze_synthetic_minimal(tmp_path):
    """End-to-end ``run(args)`` against a 3-record synthetic BAM.

    Checks:
      * exit code 0,
      * all four output files exist (3 TSVs + tagged BAM with .bai),
      * cluster manifest has 3 rows,
      * gene assignments / strand call match expectations,
      * cl0 and cl1 share an isoform_id (same gene + same termini),
      * t1t2_pairs.tsv exists and is well-formed (header only — there are no
        Type-2 records in this synthetic, so 0 data rows is correct).
    """
    from rectify.core.commands import cdna_analyze_command

    ref = _build_tiny_reference(tmp_path)
    gff = _build_tiny_gff(tmp_path)
    bam = _build_synthetic_consensus_bam(tmp_path)
    out = tmp_path / "out"

    args = _make_args(bam=bam, out=out, gff=gff, ref=ref)
    rc = cdna_analyze_command.run(args)
    assert rc == 0, "run() should return 0 on success"

    for name in ("clusters.tsv", "corrected_reads.tsv", "isoforms.tsv",
                 "t1t2_pairs.tsv",
                 "consensus_tagged.bam", "consensus_tagged.bam.bai"):
        assert (out / name).exists(), f"missing output: {name}"

    # ── clusters.tsv ─────────────────────────────────────────────────────
    lines = (out / "clusters.tsv").read_text().rstrip("\n").split("\n")
    header, rows = lines[0], lines[1:]
    assert len(rows) == 3, f"expected 3 cluster rows, got {len(rows)}"

    cols = header.split("\t")
    idx_gene = cols.index("gene")
    idx_xs = cols.index("xs")
    idx_iso = cols.index("isoform_id")
    idx_orient = cols.index("orient")

    by_qname = {}
    for row in rows:
        parts = row.split("\t")
        # cluster_id column appears as a sequential integer — the read_ids
        # column at the end pins each row to a known synthetic record.
        read_ids = parts[cols.index("read_ids")]
        by_qname[read_ids] = parts

    # cl0 + cl1 — gene_A sense
    cl0 = by_qname["rd_a1,rd_a2,rd_a3,rd_a4,rd_a5"]
    cl1 = by_qname["rd_b1,rd_b2,rd_b3"]
    cl2 = by_qname["rd_c1,rd_c2"]

    assert cl0[idx_gene] == "gene_A", f"cl0 gene={cl0[idx_gene]!r}"
    assert cl0[idx_xs] == "sense"
    assert cl1[idx_gene] == "gene_A"
    assert cl1[idx_xs] == "sense"
    # Same gene + same termini + same orient → same isoform_id.
    assert cl0[idx_iso] == cl1[idx_iso] != "", (
        f"cl0 and cl1 should share isoform_id: {cl0[idx_iso]!r} vs {cl1[idx_iso]!r}"
    )

    # cl2 — gene_B antisense (orient=fwd overlaps - strand gene_B → antisense)
    assert cl2[idx_gene] == "gene_B", f"cl2 gene={cl2[idx_gene]!r}"
    assert cl2[idx_xs] == "antisense"
    assert cl2[idx_orient] == "fwd"

    # ── corrected_reads.tsv (DRS-equivalent per-molecule) ────────────────
    cr_lines = (out / "corrected_reads.tsv").read_text().rstrip("\n").split("\n")
    cr_header, cr_rows = cr_lines[0], cr_lines[1:]
    cr_cols = cr_header.split("\t")
    # DRS-identical column names so a single loader works across modalities.
    for required in ("chrom", "strand", "corrected_3prime", "original_3prime",
                     "alignment_start", "alignment_end", "gene_id"):
        assert required in cr_cols, f"corrected_reads.tsv missing column {required!r}"
    # One row per molecule (= one per cluster).
    assert len(cr_rows) == 3, f"expected 3 molecule rows, got {len(cr_rows)}"
    _ci = {c: i for i, c in enumerate(cr_cols)}
    # All synthetic records are orient=fwd → strand '+', and
    # corrected_3prime must equal the cluster anchor.
    cl_anchor_idx = cols.index("anchor")
    anchors = {r.split("\t")[cl_anchor_idx] for r in rows}
    for row in cr_rows:
        p = row.split("\t")
        assert p[_ci["strand"]] == "+", f"orient fwd should map to strand '+', got {p[_ci['strand']]!r}"
        assert p[_ci["corrected_3prime"]] in anchors, (
            f"corrected_3prime {p[_ci['corrected_3prime']]!r} not among cluster anchors {anchors}")
        # original_3prime (+ strand) = alignment_end - 1
        assert int(p[_ci["original_3prime"]]) == int(p[_ci["alignment_end"]]) - 1

    # ── isoforms.tsv ─────────────────────────────────────────────────────
    iso_lines = (out / "isoforms.tsv").read_text().rstrip("\n").split("\n")
    # Header + at least 1 data row (gene_A's shared isoform). gene_B
    # antisense doesn't get an isoform (only "sense" genes are isoform-grouped
    # in assign_isoforms because cluster_xg is set to None for "antisense"
    # in classify_sense_antisense? — actually classify returns the gene name
    # even for antisense, so it MAY be assigned. We check for presence of
    # gene_A as a minimum.)
    iso_header = iso_lines[0]
    iso_rows = iso_lines[1:]
    iso_cols = iso_header.split("\t")
    assert "isoform_id" in iso_cols
    assert "gene" in iso_cols
    iso_gene_idx = iso_cols.index("gene")
    assert any(r.split("\t")[iso_gene_idx] == "gene_A" for r in iso_rows), (
        f"gene_A missing from isoforms.tsv ({len(iso_rows)} rows)"
    )

    # ── t1t2_pairs.tsv ───────────────────────────────────────────────────
    pair_lines = (out / "t1t2_pairs.tsv").read_text().rstrip("\n").split("\n")
    pair_header = pair_lines[0]
    # All 3 records here are Type-1 — no Type-2 → 0 data rows expected.
    assert pair_header.startswith("t1_cid"), (
        f"unexpected t1t2_pairs.tsv header: {pair_header}"
    )

    # ── tagged BAM ───────────────────────────────────────────────────────
    with pysam.AlignmentFile(str(out / "consensus_tagged.bam"), "rb") as bam_in:
        recs = list(bam_in.fetch(until_eof=True))
    assert len(recs) == 3
    for rec in recs:
        tags = {t for t, _ in rec.tags}
        assert "XA" in tags, f"{rec.query_name} missing recomputed XA"
        assert "XS" in tags, f"{rec.query_name} missing recomputed XS"


# ---------------------------------------------------------------------------
# 3) On-disk chrI smoke (opt-in)
# ---------------------------------------------------------------------------
_TEST_DATA = Path.home() / "work/ont_cdna/test_data"  # opt-in local smoke data
_CHRI_BAM = _TEST_DATA / "wt_rep1.chrI.bam"
_GFF = _TEST_DATA / "saccharomyces_cerevisiae_R64-5-1_20240529.gff"
_FASTA = _TEST_DATA / "S288C_reference_sequence_R64-5-1_20240529.fsa"


@pytest.mark.slow
@pytest.mark.skipif(
    not (_CHRI_BAM.exists() and _GFF.exists() and _FASTA.exists()),
    reason="chrI test data not on disk (dev M1 only)",
)
def test_cdna_analyze_chri_smoke(tmp_path):
    """Run cdna-analyze against the chrI fixture BAM if available.

    NOTE: ``wt_rep1.chrI.bam`` is the *raw* input BAM — it has not yet been
    through correct-cdna → align, so it does not carry XU/XC/XR/XT/XY/XF/XB
    tags. The full chain canary lives in test_cdna_chain_canary.py. Here we
    only check that the entry point handles "no usable records" gracefully
    and produces empty manifests with headers intact.
    """
    from rectify.core.commands import cdna_analyze_command

    out = tmp_path / "out"
    args = _make_args(bam=_CHRI_BAM, out=out, gff=_GFF, ref=_FASTA)
    rc = cdna_analyze_command.run(args)
    assert rc == 0

    # The BAM has no records carrying the required XU/XC/XO/XT/XY/XF tags,
    # so every record falls into the "missing required tags" branch. The
    # output manifests must still have their header rows.
    for name in ("clusters.tsv", "isoforms.tsv", "t1t2_pairs.tsv"):
        lines = (out / name).read_text().splitlines()
        assert len(lines) >= 1, f"{name} missing header row"

    # The cluster row count should be 0 (no records carry the required tags).
    rows = (out / "clusters.tsv").read_text().splitlines()[1:]
    assert len(rows) == 0, (
        f"raw chrI BAM has no XU/XC/XR/XT/XY/XF/XB tags — clusters.tsv "
        f"should be empty, got {len(rows)} rows"
    )


# ---------------------------------------------------------------------------
# 4) Carried poly(A) tail length (XA) survives a pretrimmed consensus record
# ---------------------------------------------------------------------------
def _build_pretrimmed_consensus_bam(tmp_path: Path) -> Path:
    """Three records that exercise the three tail-length provenance paths.

      cluster_0  pretrimmed (no tail in SEQ) + carried XA:i:23  → tail_len 23 (carried)
      cluster_1  tail still in SEQ            + stale XA:i:7     → walkback value (>0, ≠7)
      cluster_2  pretrimmed, no XA tag at all                    → tail_len 0
    """
    bam_path = tmp_path / "consensus_pretrimmed.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": _TINY_CHROM_LEN}]}
    umis = [
        "ACGTACGTACGTACGTACGTACGTACG",
        "TGCATGCATGCATGCATGCATGCATGC",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCC",
    ]
    specs = [
        ("cluster_0", "ACGTACGTAC", 23, "rd_p1,rd_p2"),
        ("cluster_1", "ACGTACGTAC" + "A" * 20, 7, "rd_w1,rd_w2"),
        ("cluster_2", "ACGTACGTAC", None, "rd_z1"),
    ]
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for (qname, seq, xa, rids), umi in zip(specs, umis):
            r = pysam.AlignedSegment(bam.header)
            r.query_name = qname
            r.query_sequence = seq
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 200
            r.mapping_quality = 60
            r.cigartuples = [(0, len(seq))]
            r.set_tag("XU", umi, "Z")
            r.set_tag("XO", "fwd", "Z")
            r.set_tag("XT", 1, "i")
            r.set_tag("XY", "umi_captured_fwd", "Z")
            r.set_tag("XC", len(rids.split(",")), "i")
            r.set_tag("XF", 2, "i")
            r.set_tag("XR", rids, "Z")
            if xa is not None:
                r.set_tag("XA", xa, "i")
            bam.write(r)
    pysam.index(str(bam_path))
    return bam_path


def test_cdna_analyze_honours_carried_xa_on_pretrimmed_records(tmp_path):
    """Regression for AGENT_FIXES 2026-08-21: after b3a8c35 `correct-cdna` strips
    the poly(A) tail BEFORE alignment and carries the measured length as XA:i.
    cdna-analyze used to recompute by walkback on the tail-less record and
    overwrite XA with 0 on 100 % of molecules (exit 0, no warning).

    The carried value must survive into clusters.tsv (tail_len),
    corrected_reads.tsv (polya_length) and the tagged BAM (XA); a record whose
    tail IS still in SEQ must keep the walkback value (it is the post-align
    measurement, more accurate than the pre-align one).
    """
    from rectify.core.cdna.walkback import walk_back_anchor_and_tail
    from rectify.core.commands import cdna_analyze_command

    ref = _build_tiny_reference(tmp_path)
    gff = _build_tiny_gff(tmp_path)
    bam = _build_pretrimmed_consensus_bam(tmp_path)
    out = tmp_path / "out"

    # Independent expectation for the walkback path: same record, same genome.
    with pysam.FastaFile(str(ref)) as fa:
        chrom_seq = fa.fetch("chrT").upper()
    with pysam.AlignmentFile(str(bam), "rb") as b:
        rec1 = [r for r in b.fetch(until_eof=True) if r.query_name == "cluster_1"][0]
    _anchor, expected_walkback = walk_back_anchor_and_tail(rec1, chrom_seq, "fwd")
    assert expected_walkback > 0, "fixture must have a walkable tail in SEQ"
    assert expected_walkback != 7, "fixture's stale XA must differ from the walkback value"

    rc = cdna_analyze_command.run(_make_args(bam=bam, out=out, gff=gff, ref=ref))
    assert rc == 0

    lines = (out / "clusters.tsv").read_text().rstrip("\n").split("\n")
    cols = lines[0].split("\t")
    by_rids = {row.split("\t")[cols.index("read_ids")]: row.split("\t") for row in lines[1:]}
    tl = cols.index("tail_len")
    assert int(by_rids["rd_p1,rd_p2"][tl]) == 23, "carried XA must survive a pretrimmed record"
    assert int(by_rids["rd_w1,rd_w2"][tl]) == expected_walkback, "walkback wins when the tail is in SEQ"
    assert int(by_rids["rd_z1"][tl]) == 0, "no tail anywhere → 0"

    cr = (out / "corrected_reads.tsv").read_text().rstrip("\n").split("\n")
    cr_cols = cr[0].split("\t")
    pl = {row.split("\t")[0]: int(row.split("\t")[cr_cols.index("polya_length")]) for row in cr[1:]}
    assert pl["ACGTACGTACGTACGTACGTACGTACG"] == 23
    assert pl["TGCATGCATGCATGCATGCATGCATGC"] == expected_walkback
    assert pl["CCCCCCCCCCCCCCCCCCCCCCCCCCC"] == 0

    with pysam.AlignmentFile(str(out / "consensus_tagged.bam"), "rb") as b:
        xa = {r.query_name: r.get_tag("XA") for r in b.fetch(until_eof=True)}
    assert xa == {"cluster_0": 23, "cluster_1": expected_walkback, "cluster_2": 0}


def test_cdna_analyze_warns_when_tail_len_is_dead(tmp_path, caplog):
    """A library where neither SEQ nor XA carries a tail must warn loudly — a
    silently dead column is exactly the failure the fix is for."""
    import logging as _logging
    from rectify.core.commands import cdna_analyze_command

    ref = _build_tiny_reference(tmp_path)
    gff = _build_tiny_gff(tmp_path)
    bam_path = tmp_path / "dead_tail.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": _TINY_CHROM_LEN}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(3):
            r = pysam.AlignedSegment(bam.header)
            r.query_name = f"cluster_{i}"
            r.query_sequence = "ACGTACGTAC"
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 200
            r.mapping_quality = 60
            r.cigartuples = [(0, 10)]
            r.set_tag("XU", "ACGT" * 6 + chr(65 + i) * 3, "Z")
            r.set_tag("XO", "fwd", "Z")
            r.set_tag("XT", 1, "i")
            r.set_tag("XY", "umi_captured_fwd", "Z")
            r.set_tag("XC", 1, "i")
            r.set_tag("XF", 2, "i")
            r.set_tag("XR", f"rd_{i}", "Z")
            r.set_tag("XA", 0, "i")
            bam.write(r)
    pysam.index(str(bam_path))

    with caplog.at_level(_logging.WARNING, logger="cdna-analyze"):
        rc = cdna_analyze_command.run(_make_args(bam=bam_path, out=tmp_path / "out", gff=gff, ref=ref))
    assert rc == 0
    assert any("tail_len is 0 for" in m for m in caplog.messages), caplog.messages
