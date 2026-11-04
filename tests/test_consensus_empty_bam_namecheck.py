"""An aligner that wrote an EMPTY BAM (the deSALT bug #49 fallback) must not fail the read-name compatibility check (Chanfreau 1248)."""
import pysam
import pytest

from rectify.core.consensus.consensus import _check_read_name_compatibility


def _bam(path, names):
    hdr = {'HD': {'VN': '1.6', 'SO': 'queryname'}, 'SQ': [{'SN': 'I', 'LN': 1000}]}
    with pysam.AlignmentFile(str(path), 'wb', header=hdr) as out:
        for n in names:
            a = pysam.AlignedSegment(out.header)
            a.query_name = n; a.query_sequence = 'ACGT' * 5; a.flag = 0; a.reference_id = 0; a.reference_start = 10
            a.mapping_quality = 60; a.cigarstring = '20M'; a.query_qualities = pysam.qualitystring_to_array('I' * 20)
            out.write(a)
    return str(path)


def test_empty_bam_is_skipped(tmp_path):
    names = ['r000_%d' % i for i in range(20)]
    bams = {'minimap2': _bam(tmp_path / 'm.bam', names), 'uLTRA': _bam(tmp_path / 'u.bam', names),
            'deSALT': _bam(tmp_path / 'd.bam', [])}
    _check_read_name_compatibility(bams)          # must not raise


def test_real_mismatch_still_raises(tmp_path):
    bams = {'minimap2': _bam(tmp_path / 'm.bam', ['r000_%d' % i for i in range(20)]),
            'deSALT': _bam(tmp_path / 'd.bam', ['other_%d' % i for i in range(20)])}
    with pytest.raises(RuntimeError, match='Read-id format mismatch'):
        _check_read_name_compatibility(bams)
