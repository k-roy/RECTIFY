from collections import defaultdict

import pysam

from rectify.core.chunking.sidecar import ReadNumSidecar, ReadNumSidecarWriter
from rectify.core.consensus import consensus as c
from rectify.core.consensus.consensus import ConsensusResult, _iter_name_grouped_bams


def _header():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': 1000}],
    })


def _read(header, name, rn=None, pos=10):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = pos
    r.mapping_quality = 60
    r.cigarstring = '4M'
    r.query_sequence = 'ACGT'
    r.query_qualities = pysam.qualitystring_to_array('IIII')
    if rn is not None:
        r.set_tag('RN', rn, value_type='i')
    return r


def _write_bam(path, reads):
    header = _header()
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for r in reads:
            out.write(r)


def test_iter_name_grouped_bams_uses_rn_when_all_inputs_have_it(tmp_path):
    header = _header()
    a = tmp_path / 'a.bam'
    b = tmp_path / 'b.bam'
    _write_bam(a, [_read(header, 'mutated_a', rn=1)])
    _write_bam(b, [_read(header, 'different_b', rn=1)])
    groups = list(_iter_name_grouped_bams({'a': str(a), 'b': str(b)}))
    assert len(groups) == 1
    assert groups[0][0] == 1
    assert set(groups[0][1]) == {'a', 'b'}


def test_iter_name_grouped_bams_falls_back_when_any_input_lacks_rn(tmp_path):
    header = _header()
    a = tmp_path / 'a.bam'
    b = tmp_path / 'b.bam'
    _write_bam(a, [_read(header, 'read1', rn=1)])
    _write_bam(b, [_read(header, 'read1')])
    groups = list(_iter_name_grouped_bams({'a': str(a), 'b': str(b)}))
    assert len(groups) == 1
    assert groups[0][0] == 'read1'


def test_process_batch_sidecar_tags_overwrite_cdna_keep_guard_for_others(tmp_path, monkeypatch):
    """Since f20a8e6 the sidecar/FASTQ-comment is AUTHORITATIVE for the
    Stage-1 cDNA tag set (_CDNA_COMMENT_TAGS): a colliding aligner-emitted
    tag (uLTRA's own XC:Z:NO_SPLICE / XA:Z:) must be OVERWRITTEN, or
    cdna-analyze drops the record as "missing required tags" (the 7% biased
    uLTRA-won molecule drop, planning/659). Tags OUTSIDE that set keep the
    original don't-clobber guard."""
    sidecar_path = tmp_path / 'sample.read_num_sidecar.parquet'
    with ReadNumSidecarWriter(sidecar_path, sample_id='sample') as w:
        w.add(3, 'read3', 'XA:Z:tail\tXC:i:42\tXF:Z:full\tZZ:Z:sidecar',
              'chunk_000_of_001', 'ACGT', 'IIII')
    sc = ReadNumSidecar.open(sidecar_path)

    header = _header()
    read = _read(header, 'read3', rn=3)
    read.set_tag('XA', 'existing')      # cDNA tag: sidecar must overwrite
    read.set_tag('ZZ', 'aligner')       # non-cDNA tag: guard must hold
    out_path = tmp_path / 'out.bam'
    with pysam.AlignmentFile(out_path, 'wb', header=header) as out:
        monkeypatch.setattr(c, 'select_best_alignment', lambda *args, **kwargs: ConsensusResult(
            read_id='read3', best_aligner='mapPacBio', best_alignment=None,
            aligners_compared=['mapPacBio'], confidence='high', n_aligners_agree=1,
        ))
        stats = {
            'consensus_high': 0, 'consensus_medium': 0, 'consensus_low': 0,
            '5prime_rescued': 0, 'tied_score': 0, 'by_aligner': defaultdict(int),
            'by_aligner_combo': defaultdict(int),
        }
        c._process_and_write_batch(
            [('read3', {'mapPacBio': object()})], [('read3', {'mapPacBio': read})],
            {'chrI': 'A' * 1000}, None, out, stats, read_num_sidecar=sc,
        )
    with pysam.AlignmentFile(out_path, 'rb') as inp:
        rec = next(inp)
        assert rec.get_tag('XA') == 'tail'      # authoritative sidecar wins
        assert rec.get_tag('XC') == 42
        assert rec.get_tag('XF') == 'full'
        assert rec.get_tag('ZZ') == 'aligner'   # non-cDNA guard preserved
