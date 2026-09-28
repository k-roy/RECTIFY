"""The writer audits every record it emits against the TSV row it replayed (2026-09-28).

The corrected TSV is written a pipeline stage before any BAM writer runs. Two ways it has advertised
something the BAM did not carry, both silent until a census compared them:

* the first ISSUE-083 re-split census lost 63 junctions: the writer's row lacked a column, the 5'
  rescue fell to another branch, and the TSV kept both the junction and an empty refusal;
* a 3' walkback that clips away a terminal exon leaves the TSV listing a junction the record no longer
  has (ISSUE-085; 3 of 62,602 cohort reads, all ``polya_walkback``).

Now every written record whose N ops or 5' verdict contradict its row is tagged ``Xh:Z`` with tokens
``5p:<tsv>><bam>``, ``jx-:<s>-<e>`` (a TSV junction the record lacks), ``jx+:<s>-<e>`` (an N op the TSV
lacks), counted, and summarized as a WARNING; the writers' stats carry ``tsv_bam_disagree``.
"""
import logging
import random
from array import array

import pysam
import pytest

from rectify.core.bam import bam_writer as BW

random.seed(928)

_TSV_COLUMNS = ['read_id', 'corrected_3prime', 'strand', 'five_prime_position', 'five_prime_rescued',
                'five_prime_soft_clip_length', 'five_prime_exon_cigar', 'five_prime_exon2_prefix',
                'five_prime_exon2_cigar', 'five_prime_rescue_refused', 'junctions']


@pytest.fixture(autouse=True)
def _fresh_counts():
    BW.report_writer_audit('', log=False)
    yield
    BW.report_writer_audit('', log=False)


def _rand(n):
    return ''.join(random.choice('ACGT') for _ in range(n))


def _read(header, name, start, ops, seq, reverse):
    r = pysam.AlignedSegment(header)
    r.query_name = name
    r.reference_id = 0
    r.reference_start = start
    r.flag = 16 if reverse else 0
    r.cigartuples = ops
    r.query_sequence = seq
    r.query_qualities = array('B', [30] * len(seq))
    return r


def _row(**over):
    """A correction dict shaped like `_load_corrections_from_single_tsv`'s output."""
    row = dict(corrected_3prime=0, strand='+', five_prime_position=None, five_prime_rescued=False,
               five_prime_soft_clip=0, five_prime_exon_cigar='', five_prime_upstream_trim=0,
               reanchor_clip_len=0, five_prime_exon2_prefix=0, five_prime_exon2_cigar='',
               five_prime_intron_clip_pos=-1, five_prime_clip_origin='', tail_correction_enabled=True,
               sc_homopolymer_extension=0, sc_rescued_seq='', sc_original_softclip_len=0,
               oc_homopolymer_extension=0, oc_overcall_count=0, oc_terminal_base='',
               station_b_microexons='', station_b_alternatives='', station_b_applied=0,
               station_b_intron_start='', station_b_intron_end='',
               tsv_five_prime_rescue_refused='', tsv_junctions='')
    row.update(over)
    return row


def _spliced(header):
    # 20M 100N 30M at 1000: one N at (1020, 1120)
    return _read(header, 'r1', 1000, [(0, 20), (3, 100), (0, 30)], 'A' * 50, reverse=False)


# --------------------------------------------------------------------------- the token rules

def test_an_agreeing_record_is_untagged_and_a_stale_tag_is_removed():
    header = pysam.AlignmentHeader.from_references(['chr5'], [10_000])
    read = _spliced(header)
    read.set_tag(BW.WRITER_AUDIT_TAG, 'jx-:1-2')                     # left by an earlier pass
    assert BW.audit_written_record(read, _row(tsv_junctions='1020-1120'), '') == []
    assert not read.has_tag(BW.WRITER_AUDIT_TAG)
    assert BW.report_writer_audit('t') == {}


def test_each_kind_of_disagreement_is_named_counted_and_reported(caplog):
    header = pysam.AlignmentHeader.from_references(['chr5'], [10_000])
    read = _spliced(header)
    row = _row(five_prime_rescued=True, tsv_five_prime_rescue_refused='',
               tsv_junctions='900-950;1020-1120')
    tokens = BW.audit_written_record(read, row, 'softclipped_no_junction')
    assert tokens == ['5p:->softclipped_no_junction', 'jx-:900-950']
    assert read.get_tag(BW.WRITER_AUDIT_TAG) == '5p:->softclipped_no_junction;jx-:900-950'

    other = _spliced(header)
    other.query_name = 'r2'
    assert BW.audit_written_record(other, _row(tsv_junctions=''), '') == ['jx+:1020-1120']

    third = _spliced(header)
    third.query_name = 'r3'
    BW.audit_written_record(third, _row(tsv_junctions=''), '', count=False)   # tagged, not counted
    assert third.has_tag(BW.WRITER_AUDIT_TAG)

    with caplog.at_level(logging.WARNING, logger=BW.logger.name):
        counts = BW.report_writer_audit('unit')
    assert counts == {'records': 2, '5p': 1, 'jx-': 1, 'jx+': 1}
    text = caplog.text
    assert 'do NOT match their corrected-TSV row' in text and 'r1' in text and 'r2' in text
    assert BW.report_writer_audit('unit') == {}                          # reset after reporting


def test_the_5p_verdict_is_only_audited_on_rows_that_claim_a_rescue():
    header = pysam.AlignmentHeader.from_references(['chr5'], [10_000])
    read = _spliced(header)
    # A row the processor downgraded (refused, five_prime_rescued cleared): the writer skips the
    # surgery, so its empty verdict is the expected one, not a disagreement.
    row = _row(five_prime_rescued=False, tsv_five_prime_rescue_refused='extend_refused',
               tsv_junctions='1020-1120')
    assert BW.audit_written_record(read, row, '') == []


def test_a_row_from_an_older_tsv_schema_is_not_audited():
    header = pysam.AlignmentHeader.from_references(['chr5'], [10_000])
    read = _spliced(header)
    row = _row(five_prime_rescued=True, tsv_five_prime_rescue_refused=None, tsv_junctions=None)
    assert BW.audit_written_record(read, row, 'extend_refused') == []
    assert not read.has_tag(BW.WRITER_AUDIT_TAG)


# --------------------------------------------------------------------------- through the writers

def _write_tsv(path, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(_TSV_COLUMNS) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in _TSV_COLUMNS) + '\n')


def _write_bam(path, header, reads):
    with pysam.AlignmentFile(str(path), 'wb', header=header) as fh:
        for r in sorted(reads, key=lambda x: x.reference_start):
            fh.write(r)
    pysam.index(str(path))


def _written(path):
    with pysam.AlignmentFile(str(path), 'rb') as fh:
        return {r.query_name: r for r in fh}


@pytest.mark.parametrize('writer', ['hard', 'soft', 'dual'])
def test_walkback_that_clips_away_an_exon_is_flagged(tmp_path, writer, caplog):
    """ISSUE-085's geometry: `25S 6M 650I 2046N 45M 3I 5M 2S` (minus) clipped to 40832463 loses the
    whole 6M, so the 2046-bp N goes with it, while the row still lists 40832463-40834509."""
    header = pysam.AlignmentHeader.from_references(['chr5'], [200_000_000])
    ops = [(4, 25), (0, 6), (1, 650), (3, 2046), (0, 45), (1, 3), (0, 5), (4, 2)]
    # Mixed bases, not the ISSUE-085 test's homopolymers: the writer also hard-clips a trailing genomic
    # A-run (a T-run at a minus read's left end), which would eat an all-T body block here.
    seq = _rand(25) + 'GCAGCA' + _rand(650) + 'G' + _rand(44) + _rand(3) + 'GCAGC' + _rand(2)
    clipped = _read(header, 'clipped', 40832457, ops, seq, reverse=True)
    kept = _read(header, 'kept', 40832457, ops, seq, reverse=True)
    kept.query_name = 'kept'
    bam = tmp_path / 'in.bam'
    _write_bam(bam, header, [clipped, kept])
    tsv = tmp_path / 'c.tsv'
    _write_tsv(tsv, [
        dict(read_id='clipped', corrected_3prime=40832463, strand='-', five_prime_rescued=0,
             five_prime_rescue_refused='', junctions='40832463-40834509'),
        # the control keeps its 3' end where the aligner put it, so its N survives
        dict(read_id='kept', corrected_3prime=40832457, strand='-', five_prime_rescued=0,
             five_prime_rescue_refused='', junctions='40832463-40834509'),
    ])
    with caplog.at_level(logging.WARNING, logger=BW.logger.name):
        if writer == 'hard':
            stats = BW.write_corrected_bam(str(bam), str(tsv), str(tmp_path / 'o.bam'))
            outs = [tmp_path / 'o.bam']
        elif writer == 'soft':
            stats = BW.write_softclipped_bam(str(bam), str(tsv), str(tmp_path / 'o.bam'))
            outs = [tmp_path / 'o.bam']
        else:
            stats, _ = BW.write_dual_bam(str(bam), str(tsv), str(tmp_path / 'h.bam'),
                                         str(tmp_path / 's.bam'))
            outs = [tmp_path / 'h.bam', tmp_path / 's.bam']
    assert stats['tsv_bam_disagree'] == 1                                # dual counts once
    for out in outs:
        recs = _written(out)
        assert recs['clipped'].get_tag(BW.WRITER_AUDIT_TAG) == 'jx-:40832463-40834509'
        assert not recs['kept'].has_tag(BW.WRITER_AUDIT_TAG)
    assert 'clipped' in caplog.text and 'do NOT match' in caplog.text


def _rescue_fixture(tmp_path):
    """A plus-strand read whose aligner clipped 20 exon-1 bases + 2 exon-2 bases (the re-split tests'
    geometry): exon 1 [0,60) · GT…AG intron [60,260) · exon 2 [260,320)."""
    exon1 = _rand(58) + 'CA'
    intron = 'GT' + _rand(196) + 'AG'
    exon2 = 'GG' + _rand(58)
    g = exon1 + intron + exon2 + _rand(40)
    s, e = 60, 260
    header = pysam.AlignmentHeader.from_references(['chrT'], [len(g)])
    seq = g[s - 20:s] + g[e:e + 50]

    def read(name):
        return _read(header, name, e + 2, [(4, 22), (0, 48)], seq, reverse=False)

    return g, s, e, header, read


def test_a_rescue_the_writer_refuses_while_the_tsv_says_drawn_is_flagged(tmp_path, caplog):
    g, s, e, header, read = _rescue_fixture(tmp_path)
    bam = tmp_path / 'in.bam'
    _write_bam(bam, header, [read('drawn'), read('refused')])
    # five_prime_position is the exon-1 base beside the donor (s - 1), as the re-split tests pass it;
    # corrected_3prime is the read's own 3' end, so no 3' edit runs.
    rescue = dict(strand='+', corrected_3prime=e + 49, five_prime_position=s - 1, five_prime_rescued=1,
                  five_prime_soft_clip_length=22, five_prime_exon_cigar='20M', five_prime_exon2_prefix=2,
                  five_prime_rescue_refused='', junctions='%d-%d' % (s, e))
    _write_tsv(tmp_path / 'c.tsv', [
        dict(rescue, read_id='drawn'),
        # an exon-2 head claiming 62 query bases from a 48-base body cannot be cut: extend refuses,
        # yet the row says drawn (empty refusal) and lists the junction, as the first census's rows did
        dict(rescue, read_id='refused', five_prime_exon2_cigar='62M'),
    ])
    with caplog.at_level(logging.WARNING, logger=BW.logger.name):
        stats = BW.write_corrected_bam(str(bam), str(tmp_path / 'c.tsv'), str(tmp_path / 'o.bam'),
                                       genome={'chrT': g})
    recs = _written(tmp_path / 'o.bam')
    assert [(o, n) for o, n in recs['drawn'].cigartuples if o == 3] == [(3, e - s)]
    assert not recs['drawn'].has_tag(BW.WRITER_AUDIT_TAG)
    assert recs['refused'].get_tag(BW.WRITER_AUDIT_TAG) == '5p:->%s;jx-:%d-%d' % (BW.REFUSAL_EXTEND, s, e)
    assert stats['tsv_bam_disagree'] == 1
    assert 'refused' in caplog.text


def test_the_parallel_writer_reports_the_region_workers_total(tmp_path, caplog):
    pytest.importorskip('rectify.core.bam.bam_writer_parallel')
    from rectify.core.bam.bam_writer_parallel import write_corrected_bam_parallel
    g, s, e, header, read = _rescue_fixture(tmp_path)
    bam = tmp_path / 'in.bam'
    _write_bam(bam, header, [read('drawn'), read('refused')])
    # five_prime_position is the exon-1 base beside the donor (s - 1), as the re-split tests pass it;
    # corrected_3prime is the read's own 3' end, so no 3' edit runs.
    rescue = dict(strand='+', corrected_3prime=e + 49, five_prime_position=s - 1, five_prime_rescued=1,
                  five_prime_soft_clip_length=22, five_prime_exon_cigar='20M', five_prime_exon2_prefix=2,
                  five_prime_rescue_refused='', junctions='%d-%d' % (s, e))
    _write_tsv(tmp_path / 'c.tsv', [dict(rescue, read_id='drawn'),
                                    dict(rescue, read_id='refused', five_prime_exon2_cigar='62M')])
    with caplog.at_level(logging.WARNING):
        result = write_corrected_bam_parallel(
            str(bam), str(tmp_path / 'c.tsv'), str(tmp_path / 'o.bam'),
            n_threads=1, genome={'chrT': g}, tmp_dir=str(tmp_path / 'tmp'))
    assert result['tsv_bam_disagree'] == 1
    recs = _written(tmp_path / 'o.bam')
    assert recs['refused'].get_tag(BW.WRITER_AUDIT_TAG).startswith('5p:->')
    assert 'refused' in caplog.text
