"""Unit tests for _normalize_bam_read_name and consensus invariant check.

Regression coverage for the v0.9.1 wt_R1 bug where BBmap retained FASTQ
comments (``"SRR22434624.1654499 1654499 length=76"``) while BWA truncated
to the bare accession (``"SRR22434624.869918"``), causing the K-way merge
in consensus.py to never match cross-aligner reads — silently emitting
both rows instead of one consensus winner per read.
"""

import pytest

from rectify.core.consensus.consensus import _normalize_bam_read_name


class TestNormalizeBamReadName:
    def test_bare_name_passthrough(self):
        assert _normalize_bam_read_name("SRR22434624.869918") == "SRR22434624.869918"

    def test_bbmap_fastq_comment_space_separated(self):
        # BBmap retains the FASTQ comment after the accession.
        assert _normalize_bam_read_name(
            "SRR22434624.1654499 1654499 length=76"
        ) == "SRR22434624.1654499"

    def test_bbmap_fastq_comment_after_samtools_sort(self):
        # samtools may convert internal spaces in QNAME to underscores after sort.
        assert _normalize_bam_read_name(
            "SRR22434624.1654499_1654499_length=76"
        ) == "SRR22434624.1654499"

    def test_mappacbio_pt_tag_space(self):
        # Pre-existing mapPacBio case (regression coverage).
        assert _normalize_bam_read_name("abcdef-uuid pt:i:25") == "abcdef-uuid"

    def test_mappacbio_pt_tag_underscore(self):
        assert _normalize_bam_read_name("abcdef-uuid_pt:i:25") == "abcdef-uuid"

    def test_bwa_already_clean(self):
        # BWA truncates to first whitespace token at BAM emission — already clean.
        assert _normalize_bam_read_name("SRR22434624.869918") == "SRR22434624.869918"

    def test_cross_aligner_match(self):
        # The bug: bbmap and bwa names must normalize to the same key.
        bbmap = "SRR22434624.1654499 1654499 length=76"
        bwa = "SRR22434624.1654499"
        assert _normalize_bam_read_name(bbmap) == _normalize_bam_read_name(bwa)

    def test_underscored_match(self):
        # After samtools sort, BBmap's space → underscore. Still must match BWA's bare.
        bbmap = "SRR22434624.1654499_1654499_length=76"
        bwa = "SRR22434624.1654499"
        assert _normalize_bam_read_name(bbmap) == _normalize_bam_read_name(bwa)

    def test_dot_in_accession_preserved(self):
        # The `.` between SRR and integer is intra-accession, not a separator.
        assert _normalize_bam_read_name("SRR22434624.1") == "SRR22434624.1"

    def test_empty_input(self):
        assert _normalize_bam_read_name("") == ""

    def test_only_whitespace_suffix(self):
        # Edge case: trailing whitespace only.
        assert _normalize_bam_read_name("read1 ") == "read1"

    def test_uuid_with_dash_no_strip(self):
        # ONT/PacBio UUIDs contain dashes; never strip those.
        assert _normalize_bam_read_name(
            "abc-def-123e4567-e89b-12d3-a456"
        ) == "abc-def-123e4567-e89b-12d3-a456"

    def test_ont_uuid_with_runid_comment(self):
        # ONT FASTQ headers commonly carry runid= sample= etc. after the UUID.
        assert _normalize_bam_read_name(
            "abc-def-123e4567 runid=xyz sampleid=foo"
        ) == "abc-def-123e4567"
