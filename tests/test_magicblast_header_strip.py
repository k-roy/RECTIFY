"""Regression tests for two COMPASS short-read panel fixes in multi_aligner:

1. ``_write_bare_qname_fastq`` — strips multi-token FASTQ headers (rectify-split's
   ``RN:i:N`` + Casava comment) to the bare QNAME so Magic-BLAST 1.5.0 under
   ``-no_query_id_trim`` does not spill those tokens into SAM columns 2-3
   (which made ``samtools view`` reject the SAM → 0 records).
2. ``check_aligner_available(None)`` — must not raise. The COMPASS aligners pass
   ``exec_path=None`` to the availability pre-check (binary resolved in-wrapper);
   ``shutil.which(None)`` would TypeError and abort the whole align step.
"""
import gzip

from rectify.core.align.multi_aligner import (
    _write_bare_qname_fastq,
    check_aligner_available,
)

_MULTITOKEN = (
    "@K00151:464:H3HL2BBXY:2:1101:27620:1384 RN:i:0 1:N:0:ACAGTG\n"
    "ACGTACGTACGT\n"
    "+\n"
    "FFFFFFFFFFFF\n"
    "@K00151:464:H3HL2BBXY:2:1101:27620:9999 RN:i:1 1:N:0:ACAGTG\n"
    "TTTTGGGGCCCC\n"
    "+\n"
    "FFFFFFFFFFFF\n"
)


def _read_headers(path):
    with open(path) as fh:
        return [l.rstrip('\n') for i, l in enumerate(fh) if i % 4 == 0]


class TestBareQnameFastq:

    def test_strips_multitoken_header_to_bare_qname(self, tmp_path):
        src = tmp_path / 'in.fastq'
        src.write_text(_MULTITOKEN)
        dst = tmp_path / 'out.fastq'
        _write_bare_qname_fastq(str(src), str(dst))
        assert _read_headers(dst) == [
            '@K00151:464:H3HL2BBXY:2:1101:27620:1384',
            '@K00151:464:H3HL2BBXY:2:1101:27620:9999',
        ]

    def test_sequence_and_quality_preserved(self, tmp_path):
        src = tmp_path / 'in.fastq'
        src.write_text(_MULTITOKEN)
        dst = tmp_path / 'out.fastq'
        _write_bare_qname_fastq(str(src), str(dst))
        lines = dst.read_text().splitlines()
        assert lines[1] == 'ACGTACGTACGT'
        assert lines[3] == 'FFFFFFFFFFFF'
        assert lines[5] == 'TTTTGGGGCCCC'

    def test_decompresses_gz_input(self, tmp_path):
        src = tmp_path / 'in.fastq.gz'
        with gzip.open(src, 'wt') as fh:
            fh.write(_MULTITOKEN)
        dst = tmp_path / 'out.fastq'
        _write_bare_qname_fastq(str(src), str(dst))
        assert _read_headers(dst)[0] == '@K00151:464:H3HL2BBXY:2:1101:27620:1384'

    def test_already_bare_header_unchanged(self, tmp_path):
        src = tmp_path / 'in.fastq'
        src.write_text("@read1\nACGT\n+\nFFFF\n")
        dst = tmp_path / 'out.fastq'
        _write_bare_qname_fastq(str(src), str(dst))
        assert _read_headers(dst) == ['@read1']


class TestCheckAlignerAvailable:

    def test_none_returns_false_not_raises(self):
        # The bug: shutil.which(None) -> TypeError, aborting the align step.
        assert check_aligner_available(None) is False

    def test_missing_binary_returns_false(self):
        assert check_aligner_available('definitely_not_a_real_binary_xyz123') is False
