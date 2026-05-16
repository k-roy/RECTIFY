#!/usr/bin/env python3
"""Unit tests for ``rectify.core.correct.walkback``.

The walk-back algorithm under test compares the read base to the reference
base at each aligned position and stops at the first non-stop-base match.
These tests use synthetic ``pysam.AlignedSegment`` objects so no real BAM
or FASTA on disk is required.
"""
from __future__ import annotations

import pysam
import pytest

from rectify.core.correct.walkback import (
    APPLIED_NONE,
    APPLIED_WALKBACK,
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime,
    walkback_3prime_guarded,
    walkback_drs,
    walkback_drs_full,
)


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


# ---------------------------------------------------------------------------
# Terminal-position gating (no walk-back needed)
# ---------------------------------------------------------------------------
class TestTerminalGate:
    def test_terminal_nonstop_match_no_walkback(self):
        """Terminal non-A read base matching genome → correctly anchored."""
        # 10 bp aligned, read ends in C, reference at that pos is also C.
        ref = "X" * 1000 + "ACGTACGTAC" + "X" * 100
        read = _make_read(start=1000, seq="ACGTACGTAC", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_NONE
        assert orig == corr == 1009

    def test_terminal_stop_base_matching_genome_walkback_fires(self):
        """Terminal A on genomic A → walkback fires to find true CPA.

        Genome:  ...CGTACGT AAA
        Read:    ...CGTACGT AAA
        The last 3 A's could be poly-A tail aligned over the genomic A-run.
        Walkback scans inward and anchors at T (pos 1006), the last non-A
        read-genome agreement.
        """
        ref = "X" * 1000 + "CGTACGT" + "AAA" + "X" * 100
        read = _make_read(start=1000, seq="CGTACGT" + "AAA", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert orig == 1009
        assert corr == 1006  # T immediately before the A-stretch


# ---------------------------------------------------------------------------
# Walk-back actually fires when poly-A over-extends into genomic A-stretch
# ---------------------------------------------------------------------------
class TestWalkbackFires:
    def test_polya_extended_into_genomic_a_stretch_right(self):
        """Core RECTIFY case: poly-A tail aligning into a genomic A-tract.

        Genome:  ...ACGTC AAAAA...
        Read:    ...ACGTC AAAAA   (last 5 A's are basecalled tail over genomic A)
        The terminal A matches the genomic A, but walkback MUST still fire —
        non-genomically encoded poly-A routinely lands on genomic A-runs.
        The scan walks back to C (the last non-A read-genome agreement).
        """
        ref = "X" * 1000 + "ACGTC" + "AAAAA" + "X" * 100
        read = _make_read(start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert orig == 1009
        assert corr == 1004  # C immediately before the A-stretch

    def test_polya_extends_past_genomic_tract_mismatching_terminal(self):
        """Unambiguous over-extension: read has more A's than the genomic A-run.

        Genome:  ...ACGTC AAA T  X...   (3 A's then T at position 1008)
        Read:    ...ACGTC AAAA A         (4 A's then A at position 1008 — mismatch!)
        The last A in the read is basecalled tail aligned over a genomic T;
        walkback finds the canonical anchor at the C (position 1004).
        """
        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        read = _make_read(start=1000, seq="ACGTC" + "AAA" + "AA", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch
        assert orig == 1009

    def test_polya_v_primer_tip_artifact(self):
        """QuantSeq V-primer artifact: read ends in ...AAAAAAAG (terminal G)
        over a genomic A-run. Terminal G mismatches genomic A → walkback fires,
        scan continues through the A's to the canonical anchor.
        """
        ref = "X" * 1000 + "ACGTC" + "AAAAAA" + "X" * 100  # 5 ref bases + 6 A's
        read = _make_read(
            start=1000,
            seq="ACGTC" + "AAAAA" + "G",  # 5 ref + 5 A's + 1 G (artifact)
            cigar=((0, 11),),
        )
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_WALKBACK
        assert corr == 1004  # the C before the A-stretch

    def test_polya_extension_left_side_v_primer_artifact(self):
        """Mirror of the V-primer-tip artifact on the LEFT side of the alignment.

        For a minus-strand gene under a left-3'-end protocol, the read's
        terminal base (at ``reference_start``) is the V-primer tip (a G or
        other non-A) sitting over a genomic A. That terminal mismatch fires
        Case 3 of the gate and the walkback proceeds.

        ref pos 1000+:  A A A A A C G T A C
        read         :  G A A A A A G T A C   (terminal G is the artifact)
        Scan from left:
          (0,1000) G/A mismatch → skip
          (1,1001) A/A stop_base → skip
          (2,1002) A/A stop_base → skip
          (3,1003) A/A stop_base → skip
          (4,1004) A/A stop_base → skip
          (5,1005) A/C mismatch → skip
          (6,1006) G/G MATCH, non-stop → anchor at 1006.
        """
        ref = "X" * 1000 + "AAAAAC" + "GTAC" + "X" * 100
        read = _make_read(start=1000, seq="GAAAAA" + "GTAC", cigar=((0, 10),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_LEFT)
        assert applied == APPLIED_WALKBACK
        assert orig == 1000
        assert corr == 1006


# ---------------------------------------------------------------------------
# Bad input handling
# ---------------------------------------------------------------------------
class TestInputValidation:
    def test_bad_three_prime_side_raises(self):
        ref = "A" * 100
        read = _make_read(start=10, seq="ACGT", cigar=((0, 4),))
        with pytest.raises(ValueError, match="three_prime_side"):
            walkback_3prime(read, ref, "middle")

    def test_bad_stop_base_raises(self):
        ref = "A" * 100
        read = _make_read(start=10, seq="ACGT", cigar=((0, 4),))
        with pytest.raises(ValueError, match="stop_base"):
            walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT, stop_base="X")

    def test_unmapped_read_returns_unchanged(self):
        """Read with no aligned pairs → no correction, no crash."""
        ref = "A" * 100
        # CIGAR all soft-clip → no aligned (qp, rp) pairs.
        read = _make_read(start=10, seq="AAAA", cigar=((4, 4),))
        orig, corr, applied = walkback_3prime(read, ref, THREE_PRIME_SIDE_RIGHT)
        assert applied == APPLIED_NONE
        assert orig == corr


# ---------------------------------------------------------------------------
# DRS wrapper: BAM strand passes through to gene strand
# ---------------------------------------------------------------------------
class TestDrsWrapper:
    def test_drs_plus_strand_uses_right_side(self):
        ref = "X" * 1000 + "ACGTC" + "AAA" + "TT" + "X" * 100
        read = _make_read(
            start=1000, seq="ACGTC" + "AAAAA", cigar=((0, 10),), is_reverse=False
        )
        orig, corr, applied, strand = walkback_drs(read, ref)
        assert strand == "+"
        assert applied == APPLIED_WALKBACK
        assert corr == 1004

    def test_drs_minus_strand_uses_left_side_and_T_stopbase(self):
        """Realistic DRS minus-strand pattern.

        pysam returns ``query_sequence`` in alignment orientation, so a
        minus-strand DRS read has the basecalled poly(A) reverse-complemented
        to **T's at the LEFT of query_sequence**. The genome at those + strand
        positions has T's (the - strand has the A's of the genomic A-run that
        the basecaller miscalled into the poly-A tail).

        Walk-back must therefore use ``stop_base='T'`` on the read side and
        stop at the first non-T match. Here the alignment starts with 5 T's
        of poly-A noise; the true CPA anchor is at the first non-T match
        (``GTAC`` at ref offset 1005).
        """
        ref = "X" * 1000 + "TTTTT" + "GTAC" + "X" * 100
        read = _make_read(
            start=1000, seq="TTTTT" + "GTAC", cigar=((0, 9),), is_reverse=True
        )
        orig, corr, applied, strand = walkback_drs(read, ref)
        assert strand == "-"
        assert applied == APPLIED_WALKBACK
        assert corr == 1005

    def test_drs_minus_strand_terminal_nonT_match_no_walkback(self):
        """A minus-strand DRS read whose terminal base is already a non-T
        match should be left alone — the terminal gate fires."""
        ref = "X" * 1000 + "GTAC" + "X" * 100
        read = _make_read(
            start=1000, seq="GTAC", cigar=((0, 4),), is_reverse=True
        )
        orig, corr, applied, strand = walkback_drs(read, ref)
        assert strand == "-"
        assert applied == APPLIED_NONE
        assert corr == orig


# ---------------------------------------------------------------------------
# Parity: walkback_drs_full == find_polya_boundary on the bundled validation
# BAM. find_polya_boundary is the legacy DRS production walkback;
# walkback_drs_full is the protocol-agnostic guarded core. This test is the
# gate for the bam_processor.py swap.
# ---------------------------------------------------------------------------
class TestGuardedParityWithFindPolyaBoundary:
    """walkback_drs_full must produce byte-identical output to
    find_polya_boundary on every primary read in validation_reads.bam.
    """

    @pytest.fixture(scope="class")
    def genome_dict(self):
        """Load the bundled S288C genome as a dict[chrom -> str]."""
        from pathlib import Path
        try:
            from rectify.utils.genome import load_genome
        except Exception:
            pytest.skip("rectify.utils.genome.load_genome not importable")
        data_dir = Path(__file__).parent.parent / "rectify" / "data" / "genomes" / "saccharomyces_cerevisiae"
        fasta = data_dir / "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
        if not fasta.exists():
            pytest.skip(f"Bundled genome not present: {fasta}")
        return load_genome(str(fasta))

    @pytest.fixture(scope="class")
    def validation_reads(self):
        """Load all primary reads from the bundled validation BAM."""
        from pathlib import Path
        bam_path = Path(__file__).parent.parent / "rectify" / "data" / "validation" / "validation_reads.bam"
        if not bam_path.exists():
            pytest.skip(f"Validation BAM not present: {bam_path}")
        reads = []
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for r in bam:
                if r.is_secondary or r.is_supplementary or r.is_unmapped:
                    continue
                reads.append(r)
        return reads

    def test_parity_all_reads(self, genome_dict, validation_reads):
        """For every primary read in validation_reads.bam, both functions
        must agree on whether to correct, where to anchor, and the corrected
        position. Any disagreement is a regression to investigate before
        wiring walkback_drs_full into the DRS production path.
        """
        from rectify.core.indel_corrector import find_polya_boundary
        from rectify.config import CHROM_TO_GENOME

        disagreements = []
        for read in validation_reads:
            strand = "-" if read.is_reverse else "+"
            chrom = read.reference_name
            chrom_seq = genome_dict.get(chrom) or genome_dict.get(
                CHROM_TO_GENOME.get(chrom, ""), ""
            )
            if not chrom_seq:
                continue

            legacy = find_polya_boundary(read, strand, genome_dict)
            unified = walkback_drs_full(read, chrom_seq)

            # Both None → agree (no correction).
            if legacy is None and unified is None:
                continue

            # One is None, other isn't → disagreement.
            if (legacy is None) != (unified is None):
                disagreements.append(
                    (read.query_name, strand,
                     f"legacy={legacy} unified={unified}")
                )
                continue

            # Both not None → compare fields. Use 'corrected_pos' and
            # 'correction_bp' as the load-bearing values; 'original_pos' and
            # 'polya_aligned_bp' are derivative.
            if legacy["corrected_pos"] != unified["corrected_pos"]:
                disagreements.append(
                    (read.query_name, strand,
                     f"legacy.corrected_pos={legacy['corrected_pos']} "
                     f"unified.corrected_pos={unified['corrected_pos']}")
                )
            if legacy["correction_bp"] != unified["correction_bp"]:
                disagreements.append(
                    (read.query_name, strand,
                     f"legacy.correction_bp={legacy['correction_bp']} "
                     f"unified.correction_bp={unified['correction_bp']}")
                )

        if disagreements:
            msg = "\n  ".join(f"{r}/{s}: {d}" for r, s, d in disagreements)
            pytest.fail(
                f"walkback_drs_full diverges from find_polya_boundary on "
                f"{len(disagreements)} read(s):\n  {msg}"
            )
