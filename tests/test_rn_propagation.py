import pysam

from rectify.core.align.multi_aligner import _inject_rn_into_bam, _load_fastq_rn_map


def _header():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': 1000}],
    })


def _read(header, name):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = 10
    r.mapping_quality = 60
    r.cigarstring = '4M'
    r.query_sequence = 'ACGT'
    r.query_qualities = pysam.qualitystring_to_array('IIII')
    return r


def test_rn_map_and_injection(tmp_path):
    fq = tmp_path / 'reads.fastq'
    fq.write_text('@read1\tRN:i:7\tXA:Z:x\nACGT\n+\nIIII\n@read2\tRN:i:8\nACGT\n+\nIIII\n')
    assert _load_fastq_rn_map(str(fq)) == {'read1': 7, 'read2': 8}

    bam = tmp_path / 'out.bam'
    header = _header()
    with pysam.AlignmentFile(bam, 'wb', header=header) as out:
        out.write(_read(header, 'read1'))
        out.write(_read(header, 'read2'))
    n = _inject_rn_into_bam(str(bam), {'read1': 7, 'read2': 8})
    assert n == 2
    with pysam.AlignmentFile(bam, 'rb') as inp:
        assert [r.get_tag('RN') for r in inp] == [7, 8]
