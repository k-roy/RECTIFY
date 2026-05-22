"""Tests for mapPacBio per-record QNAME sanitizer and post-alignment validator.

Regression coverage for the 2026-05-20 survey findings:

  1. mapPacBio sanitizer was gated on the first FASTQ record. If the first
     read had a bare QNAME, the whole file skipped sanitization and Dorado
     metadata leaked into the BAM QNAME via samtools sort space→underscore.
  2. Sanitizer checked only space, not tab — minimap2 -y / cDNA-pipeline
     tab-aux FASTQs bypassed it entirely.
  3. No post-alignment check caught these mutations until consensus time.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import pysam
import pytest

from rectify.core.align.multi_aligner import _sanitize_mpb_fastq
from rectify.core.align.qname_validator import (
    _strip_fastq_whitespace,
    validate_post_alignment_qnames,
)


# ─── _sanitize_mpb_fastq ──────────────────────────────────────────────────────


def _read_qnames(path: Path) -> list[str]:
    out = []
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            f.readline()
            f.readline()
            f.readline()
            out.append(h[1:].rstrip('\n'))
    return out


def _write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    with open(path, 'w') as f:
        for header, seq, qual in records:
            f.write(f'@{header}\n{seq}\n+\n{qual}\n')


class TestSanitizeMpbFastq:
    def test_bare_qname_passthrough(self, tmp_path: Path):
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        _write_fastq(src, [('uuid-1', 'ACGT', 'IIII')])
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1']

    def test_space_comment_stripped(self, tmp_path: Path):
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        _write_fastq(src, [('uuid-1 runid=abc ch=42', 'ACGT', 'IIII')])
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1']

    def test_tab_aux_stripped(self, tmp_path: Path):
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        _write_fastq(src, [('uuid-1\tXA:Z:foo\tXC:i:42', 'ACGT', 'IIII')])
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1']

    def test_mixed_first_bare_subsequent_dorado(self, tmp_path: Path):
        """The exact bug the survey caught: bare first read masked all
        subsequent comment-bearing reads from sanitization."""
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        _write_fastq(src, [
            ('uuid-1', 'ACGT', 'IIII'),
            ('uuid-2 runid=abc ch=42', 'ACGT', 'IIII'),
            ('uuid-3 runid=xyz ch=99', 'ACGT', 'IIII'),
        ])
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1', 'uuid-2', 'uuid-3']

    def test_mixed_first_bare_subsequent_tab_aux(self, tmp_path: Path):
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        _write_fastq(src, [
            ('uuid-1', 'ACGT', 'IIII'),
            ('uuid-2\tXA:Z:foo', 'ACGT', 'IIII'),
        ])
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1', 'uuid-2']

    def test_overlength_truncated_to_254(self, tmp_path: Path):
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        long_id = 'a' * 500
        _write_fastq(src, [(long_id, 'ACGT', 'IIII')])
        _sanitize_mpb_fastq(str(src), str(dst))
        out = _read_qnames(dst)
        assert len(out) == 1 and len(out[0]) == 254

    def test_gzipped_input(self, tmp_path: Path):
        src = tmp_path / 'in.fastq.gz'
        dst = tmp_path / 'out.fastq'
        with gzip.open(src, 'wt') as f:
            f.write('@uuid-1 comment\nACGT\n+\nIIII\n')
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1']

    def test_tab_takes_precedence_over_later_space(self, tmp_path: Path):
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        _write_fastq(src, [('uuid-1\ttag space content', 'ACGT', 'IIII')])
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1']

    def test_space_takes_precedence_over_later_tab(self, tmp_path: Path):
        src = tmp_path / 'in.fastq'
        dst = tmp_path / 'out.fastq'
        _write_fastq(src, [('uuid-1 space\ttab', 'ACGT', 'IIII')])
        _sanitize_mpb_fastq(str(src), str(dst))
        assert _read_qnames(dst) == ['uuid-1']


# ─── validate_post_alignment_qnames ───────────────────────────────────────────


def _minimal_header() -> dict:
    return {
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': 230218}],
    }


def _write_bam(path: Path, qnames: list[str]) -> None:
    with pysam.AlignmentFile(str(path), 'wb', header=_minimal_header()) as bam:
        for i, qn in enumerate(qnames):
            a = pysam.AlignedSegment()
            a.query_name = qn
            a.query_sequence = 'ACGT'
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 100 + i
            a.mapping_quality = 60
            a.cigartuples = [(0, 4)]
            a.query_qualities = pysam.qualitystring_to_array('IIII')
            bam.write(a)


class TestStripFastqWhitespace:
    def test_bare(self):
        assert _strip_fastq_whitespace('uuid') == 'uuid'

    def test_space(self):
        assert _strip_fastq_whitespace('uuid comment') == 'uuid'

    def test_tab(self):
        assert _strip_fastq_whitespace('uuid\tXA:Z:foo') == 'uuid'

    def test_tab_before_space(self):
        assert _strip_fastq_whitespace('uuid\ttag space') == 'uuid'

    def test_space_before_tab(self):
        assert _strip_fastq_whitespace('uuid space\ttag') == 'uuid'


class TestValidatePostAlignmentQnames:
    def test_clean_passes(self, tmp_path: Path):
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [('uuid-1', 'ACGT', 'IIII'),
                          ('uuid-2', 'ACGT', 'IIII')])
        _write_bam(bam, ['uuid-1', 'uuid-2'])
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')

    def test_fastq_comment_round_trips(self, tmp_path: Path):
        # The aligner saw FASTQ with comments but emitted bare QNAMEs (correct).
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [('uuid-1 runid=abc', 'ACGT', 'IIII'),
                          ('uuid-2\tXA:Z:foo', 'ACGT', 'IIII')])
        _write_bam(bam, ['uuid-1', 'uuid-2'])
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')

    def test_coordinate_sorted_bam_sample_can_match_late_fastq_records(self, tmp_path: Path):
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [(f'uuid-{i} runid=abc', 'ACGT', 'IIII')
                          for i in range(100)])
        _write_bam(bam, [f'uuid-{i}' for i in range(80, 85)])
        validate_post_alignment_qnames(
            str(bam), str(fq), 'test-aligner', sample_size=5,
        )

    def test_auto_sanitize_repairs_whitespace(self, tmp_path: Path):
        # Default behavior: detect whitespace, repair in place, re-validate.
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [('uuid-1', 'ACGT', 'IIII')])
        _write_bam(bam, ['uuid-1 leaked comment'])
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')
        # Post-repair: BAM QNAMEs should be normalized.
        with pysam.AlignmentFile(str(bam), 'rb') as out:
            qnames = [r.query_name for r in out]
        assert qnames == ['uuid-1']

    def test_auto_sanitize_repairs_tab(self, tmp_path: Path):
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [('uuid-1', 'ACGT', 'IIII')])
        _write_bam(bam, ['uuid-1\tXA:Z:foo'])
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')
        with pysam.AlignmentFile(str(bam), 'rb') as out:
            qnames = [r.query_name for r in out]
        assert qnames == ['uuid-1']

    def test_auto_sanitize_repairs_underscore_aux_leak(self, tmp_path: Path):
        # Historical mapPacBio bug: tab → underscore, aux tag leaked into QNAME.
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [(f'uuid-{i}\tXA:Z:foo', 'ACGT', 'IIII')
                          for i in range(5)])
        _write_bam(bam, [f'uuid-{i}_XA:Z:foo' for i in range(5)])
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')
        with pysam.AlignmentFile(str(bam), 'rb') as out:
            qnames = [r.query_name for r in out]
        assert qnames == [f'uuid-{i}' for i in range(5)]

    def test_auto_sanitize_disabled_via_env_raises(self, tmp_path: Path, monkeypatch):
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [('uuid-1', 'ACGT', 'IIII')])
        _write_bam(bam, ['uuid-1 leaked comment'])
        monkeypatch.setenv('RECTIFY_NO_AUTO_QNAME_SANITIZE', '1')
        with pytest.raises(RuntimeError, match=r'test-aligner.*whitespace'):
            validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')

    def test_auto_sanitize_disabled_via_param_raises(self, tmp_path: Path):
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [('uuid-1', 'ACGT', 'IIII')])
        _write_bam(bam, ['uuid-1\tXA:Z:foo'])
        with pytest.raises(RuntimeError, match=r'test-aligner.*whitespace'):
            validate_post_alignment_qnames(
                str(bam), str(fq), 'test-aligner', auto_sanitize=False,
            )

    def test_uncovered_mutation_raises_even_with_auto_sanitize(self, tmp_path: Path):
        # Mutation that neither whitespace-strip nor _normalize_bam_read_name
        # can recover — auto-sanitize runs, doesn't help, validator raises.
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [(f'uuid-{i}', 'ACGT', 'IIII') for i in range(10)])
        _write_bam(bam, [f'uuid-{i}-phantom-suffix' for i in range(10)])
        with pytest.raises(RuntimeError, match=r'test-aligner.*do not normalize'):
            validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')

    def test_skip_validation_env(self, tmp_path: Path, monkeypatch):
        # Full skip: validator is a no-op. BAM is left untouched.
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [('uuid-1', 'ACGT', 'IIII')])
        _write_bam(bam, ['uuid-1 leaked comment'])
        monkeypatch.setenv('RECTIFY_SKIP_POST_ALIGN_VALIDATION', '1')
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')
        # BAM untouched — original leaked QNAME still present.
        with pysam.AlignmentFile(str(bam), 'rb') as out:
            qnames = [r.query_name for r in out]
        assert qnames == ['uuid-1 leaked comment']

    def test_empty_fastq_no_crash(self, tmp_path: Path):
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        fq.write_text('')
        _write_bam(bam, [])
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')

    def test_dorado_underscore_leak_is_caught(self, tmp_path: Path):
        # Simulates the historical mapPacBio bug: BAM QNAMEs leaked Dorado
        # metadata as underscore-encoded suffixes. After fix (b), the
        # extended _normalize_bam_read_name regex strips these and the
        # validator passes — round-trip OK.
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [(f'uuid-{i} runid=abc ch={i}', 'ACGT', 'IIII')
                          for i in range(5)])
        # mapPacBio (pre-fix) would have produced these:
        _write_bam(bam, [f'uuid-{i}_runid=abc_ch={i}' for i in range(5)])
        # Extended regex strips `_runid=...` so round-trip succeeds.
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')

    def test_tab_aux_underscore_leak_is_caught(self, tmp_path: Path):
        # Same shape but for the tab-aux leak the survey identified.
        fq = tmp_path / 'reads.fastq'
        bam = tmp_path / 'aln.bam'
        _write_fastq(fq, [(f'uuid-{i}\tXA:Z:foo{i}', 'ACGT', 'IIII')
                          for i in range(5)])
        _write_bam(bam, [f'uuid-{i}_XA:Z:foo{i}' for i in range(5)])
        validate_post_alignment_qnames(str(bam), str(fq), 'test-aligner')


# ─── TSV-side normalizer (corrected_consensus._normalize_read_id) ────────────


class TestNormalizeReadIdVectorized:
    """The TSV-merge hot path uses a vectorized pandas form of
    _normalize_bam_read_name. Verify it matches the scalar version for every
    documented suffix shape."""

    def test_matches_scalar_for_all_documented_patterns(self):
        import pandas as pd
        from rectify.core.consensus.corrected_consensus import _normalize_read_id
        from rectify.core.consensus.consensus import _normalize_bam_read_name

        cases = [
            'bare-uuid',
            'uuid pt:i:25',
            'uuid_pt:i:25',
            'SRR123.4 4 length=76',
            'SRR123.4_4_length=76',
            'uuid_XA:Z:foo',
            'uuid_XC:i:42',
            'uuid_ts:A:+',
            'uuid_XA:Z:foo_XC:i:42_XF:Z:bar',
            'uuid_runid=abc_ch=42',
            'uuid_flow_cell_id=FAL12345',
            'M00123_45_000000000-ABCDE_1_1101_12345_67890',
            'uuid\tXA:Z:foo',
            'uuid space comment',
        ]
        series = pd.Series(cases)
        vectorized = _normalize_read_id(series).tolist()
        scalar = [_normalize_bam_read_name(c) for c in cases]
        assert vectorized == scalar
