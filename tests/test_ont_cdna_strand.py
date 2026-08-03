#!/usr/bin/env python3
"""ONT PCR-cDNA per-read strand resolution (rectify.core.correct.protocols.ont_cdna).

The bug these tests lock down: ``rectify correct --ONT-cDNA`` used to apply the
DRS strand rule (gene strand := ``is_reverse``), which is correct only for reads
that happen to be the mRNA sense strand.  PCR-cDNA libraries are double-stranded,
so ~a third of reads are antisense and the fixed rule assigns them the wrong gene
strand AND the wrong terminus.  See ``planning/541_ont_cdna_strand_fix.md``.

The FOUR cases that matter (sense/antisense x +/- strand gene) are enumerated
explicitly, because a real-data smoke can silently contain zero antisense reads
and pass without ever exercising the flip.
"""
from __future__ import annotations

import pandas as pd
import pysam
import pytest

from rectify.core.analyze.gene_attribution import build_cds_interval_tree
from rectify.core.correct.protocols.ont_cdna import (
    EVIDENCE_GENE_OVERLAP,
    EVIDENCE_POLYA_3P,
    EVIDENCE_POLYT_5P,
    EVIDENCE_UNASSIGNED,
    EVIDENCE_XO_ORIENT,
    resolve_rna_strand,
    three_prime_position,
)

CHROM = "chrTest"
CHROM_LEN = 6000

# GENE_P: '+' strand, [1000, 2000) -> annotated 3' terminus at 1999
# GENE_M: '-' strand, [3000, 4000) -> annotated 3' terminus at 3000
GENE_P = ("GENE_P", 1000, 2000, "+")
GENE_M = ("GENE_M", 3000, 4000, "-")

# A read covering the 3'-most 200 nt of each gene.
READ_SPAN_P = (1800, 2000)   # RNA 3' end at 1999 (the '+' gene's high coord)
READ_SPAN_M = (3000, 3200)   # RNA 3' end at 3000 (the '-' gene's low coord)


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    )


def _read(start: int, end: int, is_reverse: bool, ro: str | None = None,
          xo: str | None = None, xy: str | None = None):
    """Minimal gapless aligned segment spanning [start, end)."""
    a = pysam.AlignedSegment(_header())
    a.query_name = "r"
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.flag = 16 if is_reverse else 0
    n = end - start
    a.query_sequence = "A" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.cigartuples = [(0, n)]  # nM
    if ro is not None:
        a.set_tag("ro", ro, value_type="A")
    if xo is not None:
        a.set_tag("XO", xo)
    if xy is not None:
        a.set_tag("XY", xy)
    return a


@pytest.fixture(scope="module")
def gene_trees():
    df = pd.DataFrame(
        [
            {"chrom": CHROM, "start": s, "end": e, "strand": st,
             "gene_id": g, "gene_name": g, "feature_type": "gene"}
            for (g, s, e, st) in (GENE_P, GENE_M)
        ]
    )
    return build_cds_interval_tree(df)


# ---------------------------------------------------------------------------
# The four cases, resolved from read-intrinsic tail evidence (the `ro` tag)
# ---------------------------------------------------------------------------
class TestTailEvidence:
    """gene strand is '+' iff (read_is_sense XOR is_reverse)."""

    @pytest.mark.parametrize(
        "span,is_reverse,ro,exp_strand,exp_3p,exp_evidence",
        [
            # --- '+' strand gene: RNA 3' end is the HIGH coordinate (1999) ---
            # sense read aligns forward
            (READ_SPAN_P, False, "S", "+", 1999, EVIDENCE_POLYA_3P),
            # antisense read aligns reverse -- the case the DRS rule got wrong
            (READ_SPAN_P, True, "A", "+", 1999, EVIDENCE_POLYT_5P),
            # --- '-' strand gene: RNA 3' end is the LOW coordinate (3000) ---
            # sense read aligns reverse
            (READ_SPAN_M, True, "S", "-", 3000, EVIDENCE_POLYA_3P),
            # antisense read aligns forward -- also wrong under the DRS rule
            (READ_SPAN_M, False, "A", "-", 3000, EVIDENCE_POLYT_5P),
        ],
        ids=["plus_sense", "plus_antisense", "minus_sense", "minus_antisense"],
    )
    def test_four_cases(self, span, is_reverse, ro, exp_strand, exp_3p, exp_evidence):
        read = _read(span[0], span[1], is_reverse, ro=ro)
        strand, evidence = resolve_rna_strand(read, None)
        assert strand == exp_strand
        assert evidence == exp_evidence
        assert three_prime_position(read, strand) == exp_3p

    def test_drs_rule_would_be_wrong_for_antisense(self):
        """Guard the premise: on antisense reads the old rule inverts both
        the strand and the terminus.  If this ever passes trivially the test
        above has stopped testing anything."""
        read = _read(*READ_SPAN_P, is_reverse=True, ro="A")
        drs_strand = "-" if read.is_reverse else "+"
        assert drs_strand == "-"                                  # old answer
        assert three_prime_position(read, drs_strand) == 1800     # old terminus
        assert resolve_rna_strand(read, None)[0] == "+"           # new answer
        assert three_prime_position(read, "+") == 1999            # new terminus

    def test_both_tails_resolve_as_sense(self):
        """`ro:A:B` (poly-A AND poly-T) behaves as sense: measured 98.3%
        agreement with the annotated gene strand, same as the pure sense class."""
        read = _read(*READ_SPAN_P, is_reverse=False, ro="B")
        assert resolve_rna_strand(read, None) == ("+", EVIDENCE_POLYA_3P)


# ---------------------------------------------------------------------------
# The canonical `correct-cdna` orientation tag (XO / XY)
# ---------------------------------------------------------------------------
class TestCanonicalXOTag:
    """`correct-cdna` stage 1 LABELS orientation as XO:Z:fwd|rev and `align -y`
    carries it into the BAM; `cdna-analyze` maps {fwd:'+', rev:'-'}.  A BAM from
    that path must resolve here identically, without needing `ro` or annotation.

    XO is defined on BAM SEQ, so it gives the gene strand directly and is
    INDEPENDENT of is_reverse — both is_reverse values must give the same answer
    for a given XO.
    """

    @pytest.mark.parametrize("is_reverse", [False, True], ids=["fwd_aln", "rev_aln"])
    @pytest.mark.parametrize("xo,exp", [("fwd", "+"), ("rev", "-")])
    def test_xo_determines_strand_regardless_of_is_reverse(self, xo, exp, is_reverse):
        read = _read(*READ_SPAN_P, is_reverse=is_reverse, xo=xo)
        assert resolve_rna_strand(read, None) == (exp, EVIDENCE_XO_ORIENT)

    @pytest.mark.parametrize("xy,exp", [("umi_captured_fwd", "+"),
                                        ("umi_captured_rev", "-")])
    def test_xy_subtype_used_when_xo_absent(self, xy, exp):
        read = _read(*READ_SPAN_P, is_reverse=False, xy=xy)
        assert resolve_rna_strand(read, None) == (exp, EVIDENCE_XO_ORIENT)

    def test_umi_not_captured_falls_through(self, gene_trees):
        read = _read(*READ_SPAN_P, is_reverse=False, xy="umi_not_captured")
        strand, ev = resolve_rna_strand(read, gene_trees, chrom=CHROM)
        assert ev == EVIDENCE_GENE_OVERLAP

    def test_xo_outranks_ro(self):
        """XO is per-MOLECULE consensus; `ro` is one read's tail. XO wins."""
        read = _read(*READ_SPAN_P, is_reverse=False, ro="A", xo="fwd")
        assert resolve_rna_strand(read, None) == ("+", EVIDENCE_XO_ORIENT)

    def test_agrees_with_cdna_analyze_mapping(self):
        """Lock the two XO consumers together — if cdna_analyze_command's
        {fwd:'+', rev:'-'} ever changes, this fails rather than silently
        diverging."""
        import re
        from pathlib import Path
        from rectify.core.correct.protocols.ont_cdna import ORIENT_TO_STRAND
        src = Path(__file__).resolve().parents[1] / (
            "rectify/core/commands/cdna_analyze_command.py")
        m = re.search(r'_orient_to_strand\s*=\s*(\{[^}]*\})', src.read_text())
        assert m, "cdna_analyze_command no longer defines _orient_to_strand"
        theirs = eval(m.group(1))  # noqa: S307 - fixed literal from our own repo
        assert {k: v for k, v in theirs.items() if k in ("fwd", "rev")} == ORIENT_TO_STRAND


# ---------------------------------------------------------------------------
# Fallback: maximally-overlapping annotated gene, queried on BOTH strands
# ---------------------------------------------------------------------------
class TestGeneOverlapFallback:
    @pytest.mark.parametrize(
        "span,is_reverse,exp_strand,exp_3p",
        [
            (READ_SPAN_P, False, "+", 1999),
            (READ_SPAN_P, True, "+", 1999),
            (READ_SPAN_M, True, "-", 3000),
            (READ_SPAN_M, False, "-", 3000),
        ],
        ids=["plus_fwd", "plus_rev", "minus_rev", "minus_fwd"],
    )
    def test_resolves_without_tail_evidence(self, gene_trees, span, is_reverse,
                                            exp_strand, exp_3p):
        read = _read(span[0], span[1], is_reverse)  # no `ro` tag
        strand, evidence = resolve_rna_strand(read, gene_trees, chrom=CHROM)
        assert (strand, evidence) == (exp_strand, EVIDENCE_GENE_OVERLAP)
        assert three_prime_position(read, strand) == exp_3p

    def test_tail_evidence_outranks_gene_overlap(self, gene_trees):
        """Precedence: read-intrinsic evidence wins, so a genuinely antisense
        transcript is not force-assigned to the sense gene."""
        read = _read(*READ_SPAN_P, is_reverse=False, ro="A")
        strand, evidence = resolve_rna_strand(read, gene_trees, chrom=CHROM)
        assert (strand, evidence) == ("-", EVIDENCE_POLYT_5P)

    def test_intergenic_read_is_unassigned_not_guessed(self, gene_trees):
        read = _read(5000, 5200, is_reverse=False)
        assert resolve_rna_strand(read, gene_trees, chrom=CHROM) == (
            None, EVIDENCE_UNASSIGNED)

    def test_no_annotation_and_no_tag_is_unassigned(self):
        read = _read(*READ_SPAN_P, is_reverse=False)
        assert resolve_rna_strand(read, None) == (None, EVIDENCE_UNASSIGNED)


# ---------------------------------------------------------------------------
# End-to-end through correct_read_3prime -- proves the flag is actually wired
# ---------------------------------------------------------------------------
class TestCorrectReadWiring:
    @pytest.mark.parametrize(
        "span,is_reverse,ro,exp_strand,exp_3p",
        [
            (READ_SPAN_P, False, "S", "+", 1999),
            (READ_SPAN_P, True, "A", "+", 1999),
            (READ_SPAN_M, True, "S", "-", 3000),
            (READ_SPAN_M, False, "A", "-", 3000),
        ],
        ids=["plus_sense", "plus_antisense", "minus_sense", "minus_antisense"],
    )
    def test_ont_cdna_flag_sets_strand_and_terminus(self, span, is_reverse, ro,
                                                    exp_strand, exp_3p):
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read(span[0], span[1], is_reverse, ro=ro)
        # Non-A/T genome so no walkback fires and the terminus stays put.
        genome = {CHROM: "C" * CHROM_LEN}
        out = correct_read_3prime(read, genome, ont_cDNA=True,
                                  apply_3ss_rescue=False)
        assert len(out) == 1
        assert out[0]["strand"] == exp_strand
        assert out[0]["original_3prime"] == exp_3p
        assert out[0]["strand_evidence"] in (EVIDENCE_POLYA_3P, EVIDENCE_POLYT_5P)

    def test_without_the_flag_the_old_drs_rule_still_applies(self):
        """Regression guard for DRS and every other protocol: the new code path
        must be reachable ONLY via ont_cDNA=True."""
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read(*READ_SPAN_P, is_reverse=True, ro="A")
        genome = {CHROM: "C" * CHROM_LEN}
        out = correct_read_3prime(read, genome, apply_3ss_rescue=False)
        assert out[0]["strand"] == "-"
        assert out[0]["original_3prime"] == 1800
        assert out[0]["strand_evidence"] == ""


class TestGeneAttributionUsesResolvedStrand:
    """`gene_id` must be looked up on the RESOLVED RNA strand, not `is_reverse`.

    Found 2026-08-02 on WT_BY4742_rep1: reads whose orientation resolved from a
    5' poly-T received a `gene_id` only 53.4% of the time (vs 99.1% for 3'
    poly-A reads) — and the ones that did resolve were matched to a gene on the
    OPPOSITE strand, i.e. mis-attributed rather than merely dropped. Every
    gene-based metric silently lost or corrupted the antisense half.
    """

    def test_antisense_read_is_attributed_to_the_correct_gene(self, gene_trees):
        from rectify.core.analyze.gene_attribution import compute_read_gene_attribution

        # Antisense read over the '+' gene: aligns reverse, so is_reverse-based
        # lookup queries the '-' tree and finds nothing.
        read = _read(*READ_SPAN_P, is_reverse=True, ro="A")
        assert resolve_rna_strand(read, None)[0] == "+"

        by_alignment = compute_read_gene_attribution(read, gene_trees, chrom=CHROM)
        by_rna = compute_read_gene_attribution(read, gene_trees, chrom=CHROM,
                                               rna_strand="+")
        assert by_alignment == [], "precondition: the old behaviour finds nothing"
        assert by_rna == ["GENE_P"], "resolved strand must find the real gene"

    def test_sense_read_is_unaffected(self, gene_trees):
        from rectify.core.analyze.gene_attribution import compute_read_gene_attribution

        read = _read(*READ_SPAN_P, is_reverse=False, ro="S")
        assert (compute_read_gene_attribution(read, gene_trees, chrom=CHROM)
                == compute_read_gene_attribution(read, gene_trees, chrom=CHROM,
                                                 rna_strand="+")
                == ["GENE_P"])

    def test_none_falls_back_to_alignment_strand(self, gene_trees):
        """Other protocols pass nothing and must keep the old behaviour."""
        from rectify.core.analyze.gene_attribution import compute_read_gene_attribution

        read = _read(*READ_SPAN_P, is_reverse=False)
        assert compute_read_gene_attribution(read, gene_trees, chrom=CHROM,
                                             rna_strand=None) == ["GENE_P"]
