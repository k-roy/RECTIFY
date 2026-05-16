#!/usr/bin/env python3
"""Unit tests for the QuantSeq REV walkback wrapper.

QuantSeq REV chemistry is *antisense*: the cDNA reads are reverse-complementary
to the mRNA. After pysam-driven reverse-complementation:

    is_reverse = False  ->  gene is '-' strand, 3' end at the LEFT
                          ->  polyT (RC of basecalled polyA) at the LEFT
    is_reverse = True   ->  gene is '+' strand, 3' end at the RIGHT
                          ->  polyA at the RIGHT

These tests mirror the structure of ``tests/test_walkback_readvsref.py``'s
``TestDrsWrapper`` but with the strand inversion that QuantSeq REV adds.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.correct.protocols.quantseq_rev import walkback_quantseq_rev
from rectify.core.correct.walkback import APPLIED_NONE, APPLIED_WALKBACK


def _make_read(
    *,
    chrom: str = "chrI",
    start: int = 1000,
    seq: str = "",
    cigar=((0, 0),),
    is_reverse: bool = False,
) -> pysam.AlignedSegment:
    """Build a minimal ``pysam.AlignedSegment`` for unit testing."""
    hdr = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6"},
            "SQ": [{"SN": chrom, "LN": 100_000}],
        }
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "test"
    r.reference_name = chrom
    r.reference_start = start
    r.cigartuples = list(cigar)
    r.is_reverse = is_reverse
    r.is_unmapped = False
    r.is_secondary = False
    r.is_supplementary = False
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


class TestQuantSeqRevStrandInversion:
    """Verify the BAM-strand -> gene-strand inversion."""

    def test_is_reverse_true_yields_plus_strand_and_right_side(self):
        """is_reverse=True -> gene is '+', 3' end at RIGHT side.

        Same V-primer-tip pattern as the core right-side test but routed
        through the QuantSeq REV wrapper.
        """
        # genome + 5 ref bases + 6 A's; read carries 5 ref + 5 A's + V-primer G
        ref = "X" * 1000 + "ACGTC" + "AAAAAA" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAAAA" + "G",
            cigar=((0, 11),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch
        assert orig == 1010  # rightmost ref position (reference_end - 1)

    def test_is_reverse_false_yields_minus_strand_and_left_side(self):
        """is_reverse=False -> gene is '-', 3' end at LEFT side.

        BAM SEQ for is_reverse=False is in basecall orientation. The
        cDNA's polyA-side carries T's (the dT primer's T-tract); after
        alignment those T's appear at the LEFT of the alignment over
        a genomic T-tract (= RC of the mRNA polyA-over-A-tract on the
        gene's - strand). Walkback walks past the T-T matches and
        anchors at the first non-T read=ref agreement inward.
        """
        ref = "X" * 1000 + "TTTTT" + "CGTAC" + "X" * 100
        read = _make_read(
            start=1000,
            seq="TTTTT" + "CGTAC",
            cigar=((0, 10),),
            is_reverse=False,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "-"
        assert applied == APPLIED_WALKBACK
        assert orig == 1000
        assert corr == 1005  # C immediately after the T-stretch

    def test_inversion_is_opposite_of_drs(self):
        """Sanity check: a read with is_reverse=True should map to gene '+'
        in QuantSeq REV but to gene '-' in DRS. Use the DRS wrapper for
        the contrast.
        """
        from rectify.core.correct.walkback import walkback_drs

        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        rd_qs = _make_read(
            start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),), is_reverse=True
        )
        rd_drs = _make_read(
            start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),), is_reverse=True
        )
        _, _, _, qs_strand = walkback_quantseq_rev(rd_qs, ref)
        _, _, _, drs_strand = walkback_drs(rd_drs, ref)
        assert qs_strand == "+"
        assert drs_strand == "-"
        assert qs_strand != drs_strand


class TestQuantSeqRevWalkbackFires:
    """Walkback fires on internal-priming over a genomic A-run."""

    def test_over_extension_into_genomic_a_tract_right_side(self):
        """is_reverse=True, gene '+', poly-A over-extension at the RIGHT.

        Genome:  ...ACGTC AAA TT X...  (3 genomic A's then T at pos 1008)
        Read:    ...ACGTC AAA AA       (5 A's; last A over a genomic T -> mismatch)
        """
        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAA" + "AA",
            cigar=((0, 10),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch
        assert orig == 1009

    def test_terminal_match_no_walkback(self):
        """Non-A terminal matching genome -> no walkback (Case 1 gate)."""
        ref = "X" * 1000 + "ACGTACGTAC" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTACGTAC",
            cigar=((0, 10),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_NONE
        assert orig == corr == 1009

    def test_internal_priming_terminal_a_on_genomic_a_fires_walkback(self):
        """Core RECTIFY case for QuantSeq REV: poly-A aligning to genomic A-tract.

        is_reverse=True → gene '+', 3' end at RIGHT.
        Genome:  ...ACGTC AAAAA X...
        Read:    ...ACGTC AAAAA      (all A's on genomic A; walkback must fire)
        The scan anchors at C (pos 1004), the last non-A read-genome match.
        """
        ref = "X" * 1000 + "ACGTC" + "AAAAA" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAAAA",
            cigar=((0, 10),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK
        assert orig == 1009
        assert corr == 1004  # C before the A-stretch

    def test_internal_priming_left_side_terminal_t_on_genomic_t(self):
        """Core RECTIFY case on LEFT side, antisense form: gene '-', polyT
        (RC of basecalled polyA) aligning over a genomic T-tract.

        is_reverse=False → gene '-', 3' end at LEFT, stop_base='T'.
        On the gene's - strand the mRNA polyA over-aligns onto a genomic
        A-run; on the + strand reference the same scenario appears as
        polyT over a genomic T-run. Walkback fires and anchors at the
        first non-T inward agreement.
        """
        ref = "X" * 1000 + "TTTTT" + "CGTAC" + "X" * 100
        read = _make_read(
            start=1000,
            seq="TTTTT" + "CGTAC",
            cigar=((0, 10),),
            is_reverse=False,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "-"
        assert applied == APPLIED_WALKBACK
        assert orig == 1000
        assert corr == 1005  # C immediately after the T-stretch


class TestQuantSeqRevBoundary:
    """Edge cases the wrapper should pass through to the core."""

    def test_unmapped_no_aligned_pairs(self):
        """Read with all soft-clip CIGAR -> wrapper returns APPLIED_NONE."""
        ref = "A" * 100
        read = _make_read(start=10, seq="AAAA", cigar=((4, 4),), is_reverse=False)
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert applied == APPLIED_NONE
        assert orig == corr
        assert gene_strand == "-"


class TestQuantSeqRevTerminalGateAndGappedReads:
    """Audit-identified edge cases (2026-05-15).

    The QSREV wrapper has solid coverage for the basic polyA-walkback cases,
    but three specific edge cases were uncovered:

    1. The Case 1 terminal gate — a read ending in a non-A base that matches
       the genome should *immediately* return APPLIED_NONE, with no
       walkback. Existing ``test_terminal_match_no_walkback`` uses a C-final
       sequence, but the gate's role in *skipping the entire scan* was not
       asserted at the gate level.
    2. A read composed entirely of A's at the polyA side, scanned far enough
       to reach the first non-A read-genome agreement upstream — covers the
       full traversal of the walk-back loop.
    3. A read with an N-op (splice junction) in the CIGAR, is_reverse=True
       → QSREV maps to gene '+', RIGHT-side 3' end. The walkback must scan
       the RIGHT (post-splice) span, NOT the left side. This is the most
       likely place for a strand-inversion bug to silently land on the
       wrong end of the read.
    """

    def test_terminal_gate_case1_nonA_match_no_walkback(self):
        """Case 1 gate: terminal non-A read base matching genome — no walkback.

        is_reverse=True → QSREV gene '+', 3' end at RIGHT. Read terminal is
        C and genomic position is also C; gate returns APPLIED_NONE without
        entering the walk-back loop at all. The corrected position equals
        the original (= reference_end - 1).
        """
        # Genome at 1000..1009: ACGTACGTAC (terminal pos 1009 = C)
        ref = "X" * 1000 + "ACGTACGTAC" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGTACGTAC",
            cigar=((0, 10),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_NONE, (
            f"non-A terminal matching genome must skip walkback "
            f"(Case 1 gate); got applied={applied!r}"
        )
        assert orig == corr == 1009, (
            f"orig and corr should both be 1009 (reference_end - 1); "
            f"got orig={orig} corr={corr}"
        )

    def test_terminal_gate_case2_all_A_full_walkback(self):
        """Read composed entirely of A's at the polyA side — walk-back through
        the whole A-tract, anchor at the first non-A read-genome agreement
        upstream.

        is_reverse=True → QSREV gene '+', 3' end at RIGHT.

        Genome at 1000..1007:  AC AAAAAA   (positions 1000=A, 1001=C, then 6 A's)
        Read same as genome:   ACAAAAAA    (terminal is A; matches genomic A)
        Scan from RIGHT:
          rp=1007..1002: read A == ref A, ref is stop_base → skip
          rp=1001:       read C == ref C, C != stop → anchor at 1001
        """
        ref = "X" * 1000 + "AC" + "A" * 6 + "X" * 100  # pos 1000=A, 1001=C, 1002-1007=A
        read = _make_read(
            start=1000,
            seq="AC" + "A" * 6,
            cigar=((0, 8),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK, (
            f"all-A polyA tail over genomic A-tract must trigger walkback; "
            f"got applied={applied!r}"
        )
        assert orig == 1007, f"orig should be 1007 (reference_end - 1); got {orig}"
        assert corr == 1001, (
            f"walkback should anchor at 1001 (first non-A read-genome "
            f"agreement upstream); got corr={corr}"
        )

    def test_strand_inverted_gapped_read_walkback(self):
        """Spliced read (N-op in CIGAR), is_reverse=True → QSREV gene '+',
        RIGHT side. The walkback MUST scan the RIGHT (post-splice) span,
        not the LEFT (pre-splice) one — this is the strand-inversion-bug
        canary.

        Genome layout (0-based):
          1000..1004  ACGTC      first exon body
          1005..1054  X*50       intron (skipped by N-op in CIGAR)
          1055..1059  ACGTC      second exon body
          1060..1067  AAAAAAAA   genomic A-tract (poly-A walkback target)
          1068+       X*100      padding

        Read (CIGAR 5M 50N 13M, is_reverse=True):
          ACGTC + skip + ACGTC + AAAAAAAA   (5 + 5 + 8 = 18 query bases)
        Aligned pairs cover ref positions 1000..1004 and 1055..1067.
        Scan from RIGHT (rp=1067) walks back through the 8 genomic A's and
        anchors at rp=1059 (the C terminating exon 2).
        """
        # Build the genome
        ref = (
            "X" * 1000
            + "ACGTC"        # 1000..1004
            + "X" * 50       # 1005..1054 (intron)
            + "ACGTC"        # 1055..1059
            + "A" * 8        # 1060..1067
            + "X" * 100
        )
        read = _make_read(
            start=1000,
            seq="ACGTC" + "ACGTC" + "A" * 8,
            cigar=((0, 5), (3, 50), (0, 13)),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+", (
            f"is_reverse=True must invert to gene '+' under QSREV; "
            f"got gene_strand={gene_strand!r}"
        )
        assert applied == APPLIED_WALKBACK, (
            f"polyA over genomic A-tract on the RIGHT (post-splice) span "
            f"must trigger walkback; got applied={applied!r}"
        )
        # Original 3' end (reference_end - 1) is at the END of the 8-A tract.
        assert orig == 1067, f"orig should be 1067 (reference_end - 1); got {orig}"
        # Walkback should anchor at the C terminating exon 2 = pos 1059.
        # If the wrapper had wrongly used the LEFT side, it would have
        # anchored at pos 1000 or scanned an entirely different span.
        assert corr == 1059, (
            f"walkback should anchor at 1059 (last non-A read-genome match "
            f"on the RIGHT span); got corr={corr}. If corr is near 1000, "
            f"the wrapper is mistakenly walking the LEFT side."
        )


class TestQuantSeqRevGuardedRouting:
    """The QSrev wrapper now delegates to walkback_3prime_guarded so the
    K-context tail false-stop guard and the right-side N-op bridging guard
    are both available to QSrev reads.
    """

    def test_k_context_false_stop_skipped_inside_polya(self):
        """K-context tail-context guard fires inside the poly-A run.

        Construct a read whose poly-A tail aligns over a region whose
        genome has a single non-A base (T) buried inside a longer A-tract.
        The read at that position is A (basecalled poly-A), so it
        MISMATCHES the genomic T. The candidate non-A anchor sits one
        position further inward at a C that *does* match the genome —
        but only one A separates this C from the poly-A tail noise.

        Layout (gene '+', is_reverse=True, RIGHT side, stop_base=A):
            ref pos 1000..1014:   ACGCAGGAGT  AAAAA   (last 5 are genomic A's)
            read positions:       ACGCAGGAGT  AAAAA   (read = ref content)

        With tail_context_k=4 lookahead, the candidate at the rightmost
        non-A match (T at refp 1009) is checked: 4 inward positions look
        like GAGGA — these are NOT all stop-base, so the candidate is NOT
        treated as a false stop. Anchor lands at 1009.

        This test verifies the wrapper exposes the guarded core's
        tail-context behavior; it is the QSrev-side analogue of the
        cat3/cat6 validation reads in tests/test_walkback_readvsref.py.
        """
        ref = "X" * 1000 + "ACGCAGGAGT" + "AAAAA" + "X" * 100
        read = _make_read(
            start=1000,
            seq="ACGCAGGAGT" + "AAAAA",
            cigar=((0, 15),),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK, (
            f"poly-A over genomic A-tract must trigger walkback through the "
            f"guarded core; got applied={applied!r}"
        )
        assert orig == 1014
        assert corr == 1009, (
            f"K-context lookahead at T (1009) sees mixed bases (GAGGA) "
            f"inward — not a false stop. Anchor must land at 1009; "
            f"got corr={corr}"
        )

    def test_nop_bridging_returns_no_correction(self):
        """Right-side N-op bridging guard: an aligner-inserted intron
        within poly_noise_window of the 3' end, whose entire post-intron
        span is poly-A over a genomic A-tract, must NOT let the walkback
        cross the intron and anchor in the pre-intron exon body.

        Without the guard the scan would walk through the post-intron
        A-tract, cross the N-op silently (the scratch array skips it),
        and anchor at the last non-A position of the pre-intron span —
        a corrected position that propagates the bridging artifact. The
        guard clips scan_lo to the first post-intron scratch index; with
        no non-A anchor in the (all-A) post-intron span, the function
        returns None (= APPLIED_NONE in the wrapper).

        Layout (gene '+', is_reverse=True, RIGHT side, stop_base=A):
            refp 1000..1004:   ACGTC          (pre-intron exon body)
            refp 1005..1054:   (intron, 50bp)
            refp 1055..1064:   AAAAAAAAAA     (post-intron, all A's)
            poly_noise_window default = 50; original_3prime - last_n_end
            = 1064 - 1055 = 9 ≤ 50 → guard fires.
        """
        ref = (
            "X" * 1000
            + "ACGTC"          # 1000..1004 pre-intron
            + "X" * 50         # 1005..1054 intron (X is non-A but skipped)
            + "A" * 10         # 1055..1064 post-intron A-tract
            + "X" * 100
        )
        read = _make_read(
            start=1000,
            seq="ACGTC" + "A" * 10,  # 15 read bases total
            cigar=((0, 5), (3, 50), (0, 10)),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_NONE, (
            f"N-op bridging guard must suppress the walkback when the "
            f"post-intron span is entirely poly-A; got applied={applied!r}. "
            f"If applied=APPLIED_WALKBACK with corr near 1004, the guard "
            f"is not firing and the scan crossed the intron into the "
            f"pre-intron exon body."
        )
        assert orig == 1064
        assert corr == 1064

    def test_nop_outside_noise_window_no_clipping(self):
        """N-ops far from the 3' end are legitimate introns, not bridging
        artifacts; the right-side N-op guard must NOT fire for them.

        Layout (gene '+', is_reverse=True, RIGHT side, stop_base=A):
            refp 1000..1004:    ACGTC           (first exon)
            refp 1005..1054:    (intron, 50bp)
            refp 1055..1149:    95 bases ending in CGTAC (last exon body)
            refp 1150..1154:    AAAAA           (terminal poly-A noise)

            original_3prime    = 1154
            last_n_end         = 1055
            distance           = 1154 - 1055 = 99 > 50 (poly_noise_window)
            → guard does NOT fire; scan must operate over both exons and
            anchor at the last non-A match (refp 1149, C).
        """
        # Build a 95-bp post-intron exon body: 95 chars with no A's,
        # ending in T at refp 1149. Then 5 trailing A's of poly-A noise.
        exon2_body = "CGCG" * 23 + "CGT"  # 92 + 3 = 95 chars (no A's)
        assert "A" not in exon2_body
        ref = (
            "X" * 1000
            + "ACGTC"          # 1000..1004 first exon
            + "X" * 50         # 1005..1054 intron
            + exon2_body       # 1055..1149 (95 bases)
            + "AAAAA"          # 1150..1154 trailing poly-A
            + "X" * 100
        )
        read = _make_read(
            start=1000,
            seq="ACGTC" + exon2_body + "AAAAA",  # 5 + 95 + 5 = 105 bases
            cigar=((0, 5), (3, 50), (0, 100)),
            is_reverse=True,
        )
        orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
        assert gene_strand == "+"
        assert applied == APPLIED_WALKBACK, (
            f"Intron 99 bp from 3' end is outside poly_noise_window=50; "
            f"the N-op guard must not fire and walkback must anchor "
            f"normally; got applied={applied!r}"
        )
        assert orig == 1154
        # last non-A read=ref match: T at refp 1149 (last char of exon2_body).
        assert corr == 1149, (
            f"Anchor must land at the last non-A match (refp 1149, T at "
            f"end of exon2_body). If corr is < 1100, the N-op guard fired "
            f"incorrectly; got corr={corr}"
        )
