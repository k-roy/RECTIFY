"""Successful micro-exon edits retain exact junction provenance through BAM output."""

import copy
import json

import pysam
import pytest

from rectify.core.aggregate.junctions import aggregate_junctions
from rectify.core.bam.bam_writer import apply_corrected_edits_to_read, apply_station_b_microexons
from rectify.core.consensus.chimeric_consensus import (
    ChimericResult, ChimericSegment, build_chimeric_read, cigar_to_events,
)


HEADER = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chrT', 'LN': 2000}]})
ORPHAN = '60M400N6I100M80N40M'
DRAWN = '60M140N6M254N100M80N40M'
CALL = {'intron': [60, 460], 'exons': [[200, 206]], 'alternatives': 'chrT:300-306'}


def _read(cigar=ORPHAN, reverse=False, name='r'):
    read = pysam.AlignedSegment(HEADER)
    read.query_name = name
    read.reference_id = 0
    read.reference_start = 0
    read.is_reverse = reverse
    read.cigarstring = cigar
    length = sum(n for op, n in read.cigartuples if op in (0, 1, 4, 7, 8))
    read.query_sequence = 'C' * length
    read.query_qualities = [10 + i % 25 for i in range(length)]
    return read


def _correction(reverse=False, **overrides):
    row = {'station_b_applied': 1, 'station_b_microexons': 'chrT:200-206',
           'station_b_intron_start': 60, 'station_b_intron_end': 460,
           'station_b_alternatives': 'chrT:300-306',
           'corrected_3prime': 0 if reverse else 679, 'strand': '-' if reverse else '+'}
    row.update(overrides)
    return row


def _literal_tag(read, calls=None, **overrides):
    payload = {'v': 1, 'chrom': 'chrT', 'strand': '-' if read.is_reverse else '+',
               'calls': [copy.deepcopy(CALL)] if calls is None else calls}
    payload.update(overrides)
    read.set_tag('Xb', json.dumps(payload))


def _aggregate(tmp_path, reads):
    path = tmp_path / 'reads.bam'
    with pysam.AlignmentFile(path, 'wb', header=HEADER) as bam:
        for read in reads:
            bam.write(read)
    df = aggregate_junctions(str(path))
    return {(r.intron_start, r.intron_end): r for r in df.itertuples()}


@pytest.mark.parametrize('reverse', [False, True])
def test_live_writer_preserves_cdna_and_validation_tags_and_localizes_credit(tmp_path, reverse):
    read = _read(reverse=reverse)
    read.set_tag('XB', '7/5')
    read.set_tag('XV', 'cat9_fixture')
    original_seq, original_qual = read.query_sequence, list(read.query_qualities)
    assert apply_corrected_edits_to_read(read, _correction(reverse))
    assert read.cigarstring == DRAWN
    assert read.get_tag('XB') == '7/5'
    assert read.get_tag('XV') == 'cat9_fixture'
    assert read.query_sequence == original_seq and list(read.query_qualities) == original_qual
    assert json.loads(read.get_tag('Xb')) == {
        'v': 1, 'chrom': 'chrT', 'strand': '-' if reverse else '+', 'calls': [CALL]}

    rows = _aggregate(tmp_path, [read])
    assert {j: r.station_b_reads for j, r in rows.items()} == {
        (60, 200): 1, (206, 460): 1, (560, 640): 0}
    assert rows[560, 640].station_b_alternatives == ''
    assert rows[60, 200].station_b_alternatives == 'chrT:300-306'


def test_applied_row_whose_surgery_refuses_cannot_claim_a_draw(tmp_path):
    read = _read('60M400N2I4I100M80N40M')  # no adjacent 6I exists
    read.set_tag('XB', '4/2')
    read.set_tag('XV', 'fixture')
    apply_corrected_edits_to_read(read, _correction())
    assert read.cigarstring == '60M400N2I4I100M80N40M'
    assert not read.has_tag('Xb')
    assert read.get_tag('XB') == '4/2' and read.get_tag('XV') == 'fixture'
    assert all(r.station_b_reads == 0 for r in _aggregate(tmp_path, [read]).values())


@pytest.mark.parametrize('failed', [0, 1, None])
def test_partial_multicall_success_keeps_the_correct_alternatives(tmp_path, failed):
    insertions = ['6I', '6I']
    if failed is not None:
        insertions[failed] = '2I4I'
    read = _read(f'60M400N{insertions[0]}160M400N{insertions[1]}100M')
    row = _correction(station_b_microexons='chrT:200-206|chrT:700-706',
                      station_b_intron_start='60,620', station_b_intron_end='460,1020',
                      station_b_alternatives='chrT:300-306|chrT:800-806')
    assert apply_station_b_microexons(read, row)
    expected = [CALL, {'intron': [620, 1020], 'exons': [[700, 706]],
                       'alternatives': 'chrT:800-806'}]
    assert json.loads(read.get_tag('Xb'))['calls'] == [c for i, c in enumerate(expected) if i != failed]
    rows = _aggregate(tmp_path, [read])
    for i, junction in enumerate([(60, 200), (620, 700)]):
        if i != failed:
            assert rows[junction].station_b_reads == 1
            assert rows[junction].station_b_alternatives == expected[i]['alternatives']

    previous = read.to_string()
    assert not apply_station_b_microexons(read, row)
    assert read.to_string() == previous  # replay preserves successful provenance exactly


@pytest.mark.parametrize('reverse', [False, True])
def test_legacy_coordinates_are_unverified_and_cdna_counts_are_not_draws(tmp_path, caplog, reverse):
    old = _read(DRAWN, reverse, 'legacy')
    old.set_tag('XB', 'chrT:200-206')
    old.set_tag('XV', 'chrT:300-306')
    cdna = _read(DRAWN, reverse, 'cdna')
    cdna.set_tag('XB', '4/2')
    rows = _aggregate(tmp_path, [old, cdna])
    assert all(r.full_junction_reads == 2 and r.station_b_reads == 0
               and r.station_b_unverified_reads == 1 for r in rows.values())
    assert 'Reprocess the original alignments' in caplog.text


@pytest.mark.parametrize('overrides', [
    {'v': 2}, {'v': True}, {'chrom': 'chrOther'}, {'strand': '-'}, {'calls': {}},
    {'calls': [dict(CALL, exons=[[200, 199]])]},
    {'calls': [dict(CALL, junctions=[[560, 640]])]},
])
def test_unknown_or_invalid_schema_never_credits_a_junction(tmp_path, overrides):
    read = _read(DRAWN)
    _literal_tag(read, **overrides)
    assert all(r.station_b_reads == 0 and r.station_b_unverified_reads == 1
               for r in _aggregate(tmp_path, [read]).values())


def test_later_cigar_changes_do_not_transfer_credit_to_new_coordinates(tmp_path):
    read = _read('60M141N6M253N100M80N40M')
    _literal_tag(read)
    rows = _aggregate(tmp_path, [read])
    assert all(r.station_b_reads == 0 for r in rows.values())


def test_fallback_provenance_comes_from_the_winner_not_the_sequence_donor(tmp_path):
    winner = _read(DRAWN)
    donor = _read(DRAWN)
    result = ChimericResult('r', False, [('interior', 'winner', 0, 206)],
                           winner.cigartuples, 0, 'high', anchor_aligner='winner', is_fallback=True)
    for tagged, want in [('donor', 0), ('winner', 1)]:
        for read in [winner, donor]:
            read.set_tag('Xb', None)
        _literal_tag(donor if tagged == 'donor' else winner)
        out = build_chimeric_read(donor, 0, winner.cigartuples, result, HEADER,
                                  anchor_read=winner, aligner_reads={'winner': winner, 'donor': donor})
        rows = _aggregate(tmp_path, [out])
        assert rows[60, 200].station_b_reads == want
        assert rows[560, 640].station_b_reads == 0


def test_chimeric_provenance_follows_each_selected_junction(tmp_path):
    # Same final coordinates, two sources: only the first N is selected from the
    # arm that drew the exon. A call-wide flag would wrongly credit the second N.
    drawn = _read(DRAWN)
    unaided = _read(DRAWN)
    _literal_tag(drawn)
    events = cigar_to_events(drawn.cigartuples, 0)
    left = ChimericSegment(0, 66, 'interior', winning_aligner='drawn',
                           cigar_events={'drawn': events[:3]})
    right = ChimericSegment(66, 206, 'interior', winning_aligner='unaided',
                            cigar_events={'unaided': events[3:]})
    result = ChimericResult('r', True, [('interior', 'drawn', 0, 66),
                                      ('interior', 'unaided', 66, 206)],
                           drawn.cigartuples, 0, 'high', all_segment_scores=[left, right],
                           anchor_aligner='drawn')
    out = build_chimeric_read(unaided, 0, drawn.cigartuples, result, HEADER,
                              anchor_read=drawn, aligner_reads={'drawn': drawn, 'unaided': unaided})
    rows = _aggregate(tmp_path, [out])
    assert {j: r.station_b_reads for j, r in rows.items()} == {
        (60, 200): 1, (206, 460): 0, (560, 640): 0}
    assert json.loads(out.get_tag('Xb'))['calls'][0]['junctions'] == [[60, 200]]
