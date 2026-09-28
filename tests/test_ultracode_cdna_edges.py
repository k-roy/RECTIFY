"""cDNA frame markers must survive real TSV and writer handoffs."""
from array import array
import copy

import pysam

from rectify.core.bam import bam_writer as bw
from rectify.core.bam.bam_processor import correct_read_3prime
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.commands.cdna_analyze_command import _read_info_from_bam_record
from rectify.core.correct.protocols.ont_cdna import has_rna_sense_frame, resolve_rna_strand


def frame_record(reverse, kind, hard=0, equals=False):
    h = pysam.AlignmentHeader.from_references(['chrT'], [1000])
    r = pysam.AlignedSegment(h)
    r.query_name = f'{kind}_R{int(reverse)}_H{hard}_eq{int(equals)}'
    r.reference_id, r.reference_start, r.mapping_quality = 0, 300, 50
    r.flag, r.cigarstring = (16 if reverse else 0), (f'{hard}H80M7H' if hard else '80M')
    r.query_sequence = ('=' if equals else 'C') * 80
    r.query_qualities = array('B', [24 + i % 13 for i in range(80)])
    for tag, value in [('XU', 'ACG'*9), ('XO', 'fwd' if reverse else 'rev'),
                       ('XT', 2 if hard else 1), ('XY', 'umi_not_captured' if hard else 'umi_captured_fwd'),
                       ('XC', 3), ('XF', 1), ('XA', 12), ('XR', 'source1,source2,source3'),
                       ('XB', '2/1'), ('XQ', 0 if hard else 53), ('XK', 36),
                       ('XD', 3), ('XP', 18.5), ('XW', 1.5)]:
        r.set_tag(tag, value)
    markers = {'int1': (1, 'i'), 'uint1': (1, 'C'), 'text1': ('1', 'Z'),
               'float1p5': (1.5, 'f'), 'floatinf': (float('inf'), 'f'),
               'floatnan': (float('nan'), 'f'), 'array1': (array('i', [1]), None)}
    if kind in markers:
        value, typ = markers[kind]
        r.set_tag('XN', value, value_type=typ)
    r.set_tag('ZZ', array('h', [-8, 2, 900]))
    return r


def run_frame_roundtrip(directory):
    directory.mkdir(exist_ok=True, parents=True)
    records, rows, results = [], [], []
    genome = {'chrT': 'C' * 1000}
    for reverse in (False, True):
        for marker in ('missing', 'int1', 'uint1', 'text1', 'float1p5', 'floatinf', 'floatnan', 'array1'):
            for hard in (0, 5):
                for equals in (False, True):
                    r = frame_record(reverse, marker, hard, equals)
                    valid_marker = marker in ('int1', 'uint1')
                    expected_strand = ('-' if reverse else '+') if valid_marker else ('+' if reverse else '-')
                    expected_orient = 'fwd' if expected_strand == '+' else 'rev'
                    expected_end = 379 if expected_strand == '+' else 300
                    assert has_rna_sense_frame(r) == valid_marker
                    strand, evidence = resolve_rna_strand(r)
                    assert strand == expected_strand
                    row = correct_read_3prime(copy.deepcopy(r), genome, ont_cDNA=True,
                                             apply_3ss_rescue=False, apply_atract=False)[0]
                    assert row['strand'] == strand
                    assert row['original_3prime'] == row['corrected_3prime'] == expected_end
                    records.append(r); rows.append(row)
                    results.append(dict(read_id=r.query_name, input_sam=r.to_string(), marker=marker,
                                        expected_orient=expected_orient, expected_end=expected_end,
                                        strand_evidence=evidence, decision=row, outputs={}))
    bam, tsv = directory / 'stock.bam', directory / 'corrections.tsv'
    with pysam.AlignmentFile(str(bam), 'wb', header=records[0].header) as fh:
        for r in records:
            fh.write(r)
    with tsv.open('w') as fh:
        fh.write('\t'.join(CORRECTION_TSV_HEADER) + '\n')
        for row in rows:
            fh.write('\t'.join(correction_result_to_tsv_row(row)) + '\n')
    bw.write_corrected_bam(str(bam), str(tsv), str(directory / 'hard.bam'), genome)
    bw.write_softclipped_bam(str(bam), str(tsv), str(directory / 'soft.bam'), genome)
    bw.write_dual_bam(str(bam), str(tsv), str(directory / 'dual_hard.bam'), str(directory / 'dual_soft.bam'), genome)
    by_name = {r.query_name: r for r in records}
    for arm in ('hard', 'soft', 'dual_hard', 'dual_soft'):
        with pysam.AlignmentFile(str(directory / (arm + '.bam')), 'rb') as fh:
            outputs = {r.query_name: r for r in fh}
        assert len(outputs) == len(records)
        for result in results:
            r = outputs[result['read_id']]
            stock = by_name[r.query_name]
            info, size = _read_info_from_bam_record(r, genome['chrT'])
            assert (info.orient, info.anchor, size) == (result['expected_orient'], result['expected_end'], 3)
            assert r.cigarstring == stock.cigarstring and r.reference_start == stock.reference_start
            assert r.query_sequence == 'C' * 80 and r.query_qualities == stock.query_qualities
            for tag in ('XO', 'XU', 'XT', 'XY', 'XR', 'XB', 'XQ', 'XK', 'XD', 'XP', 'XW', 'ZZ'):
                assert r.get_tag(tag, with_value_type=True) == stock.get_tag(tag, with_value_type=True)
            result['outputs'][arm] = dict(sam=r.to_string(), info=vars(info), cluster_size=size)
    for r in results:
        assert r['outputs']['hard']['sam'] == r['outputs']['dual_hard']['sam']
        assert r['outputs']['soft']['sam'] == r['outputs']['dual_soft']['sam']
    return results


def test_typed_frame_marker_all_writers(tmp_path):
    assert len(run_frame_roundtrip(tmp_path)) == 64


def test_all_bam_integer_storage_widths_establish_frame():
    for reverse in (False, True):
        for typ in ('c', 'C', 's', 'S', 'i', 'I'):
            r = frame_record(reverse, 'missing')
            r.set_tag('XN', 1, value_type=typ)
            assert has_rna_sense_frame(r)
            assert resolve_rna_strand(r)[0] == ('-' if reverse else '+')


def test_malformed_required_counts_refuse_without_crashing():
    for reverse in (False, True):
        for field in ('XT', 'XC', 'XF'):
            for value, typ in ((float('inf'), 'f'), (float('nan'), 'f'),
                               (array('i', [1]), None), ('not-an-integer', 'Z')):
                r = frame_record(reverse, 'int1')
                r.set_tag(field, value, value_type=typ)
                assert _read_info_from_bam_record(r, 'C' * 1000) is None
