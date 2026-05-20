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

    # ── Extended-regex coverage (2026-05-20 survey gaps) ─────────────────

    def test_tab_aux_leak_xa_z(self):
        # mapPacBio (pre-sanitizer-fix) converted tab to underscore, leaving
        # the minimap2 -y / cDNA-pipeline aux tag in the QNAME.
        assert _normalize_bam_read_name("uuid_XA:Z:foo") == "uuid"

    def test_tab_aux_leak_xc_i(self):
        assert _normalize_bam_read_name("uuid_XC:i:42") == "uuid"

    def test_tab_aux_leak_ts_a(self):
        # gapmm2 strand tag in the BBmap-style leaked form.
        assert _normalize_bam_read_name("uuid_ts:A:+") == "uuid"

    def test_tab_aux_leak_multi(self):
        # Multiple aux tags concatenated by samtools sort.
        assert _normalize_bam_read_name(
            "uuid_XA:Z:foo_XC:i:42_XF:Z:bar"
        ) == "uuid"

    def test_dorado_runid_underscore_encoded(self):
        # Dorado metadata after samtools sort encodes spaces as underscores.
        assert _normalize_bam_read_name(
            "abc-def-123_runid=xyz_ch=42_start_time=2024"
        ) == "abc-def-123"

    def test_dorado_flow_cell_id_underscore_encoded(self):
        assert _normalize_bam_read_name(
            "uuid_flow_cell_id=FAL12345"
        ) == "uuid"

    def test_legitimate_underscore_qname_preserved(self):
        # Illumina-style flow-cell IDs contain underscores but no `=` or `:T:`.
        # These must NOT be mangled — only specific suffix forms are stripped.
        assert _normalize_bam_read_name(
            "M00123_45_000000000-ABCDE_1_1101_12345_67890"
        ) == "M00123_45_000000000-ABCDE_1_1101_12345_67890"

    def test_sra_accession_with_dot_underscore_preserved(self):
        # SRA-style `SRR.<id>` should survive even with trailing length=.
        assert _normalize_bam_read_name(
            "SRR123456.789_789_length=76"
        ) == "SRR123456.789"
