"""Guarded walkback must compare the query bases named by the live CIGAR."""
import csv
import random

import pysam

from rectify.core.bam import bam_processor as bp, bam_writer as bw
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.correct import walkback as wb
from rectify.utils.genome import register_genome_contigs


def _read(header, name, start, ops, sequence, reverse=False):
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.reference_id = 0
    read.reference_start = start
    read.flag = 16 if reverse else 0
    read.mapping_quality = 60
    read.cigartuples = ops
    read.query_sequence = sequence
    read.query_qualities = pysam.qualitystring_to_array('I' * len(sequence))
    read.set_tag('ZZ', 'preserve')
    return read


def test_guarded_cursor_matches_live_query_with_h_s_i_d_n():
    rng = random.Random(19052026)
    genome = ''.join(rng.choices('ACGT', k=600))
    header = pysam.AlignmentHeader.from_references(['cursorT'], [600])
    shapes = [[(0, 20)], [(1, 3), (0, 20)],
              [(0, 10), (2, 2), (0, 10)], [(0, 8), (3, 30), (0, 12)]]
    for reverse in (False, True):
        for soft in (0, 6):
            for body in shapes:
                plain_ops = ([(4, soft)] if soft else []) + body + [(4, 4), (5, 5)]
                query, reference = '', 100
                for op, length in plain_ops:
                    if op == 4:
                        query += 'T' * length
                    elif op == 1:
                        query += 'ACG'
                    elif op == 0:
                        query += genome[reference:reference + length]
                    if op in (0, 2, 3):
                        reference += length
                results = []
                for hard in (0, 7):
                    ops = ([(5, hard)] if hard else []) + plain_ops
                    read = _read(header, 'cursor', 100, ops, query, reverse)
                    original = read.to_string()
                    result = wb.walkback_drs_full(read, genome)
                    # Observer-only inspection after actual production scoring.
                    count = sum(n for op, n in body if op in (0, 2))
                    actual = [(int(wb._g_rp[i]), int(wb._g_refp[i]))
                              for i in range(count) if wb._g_rp[i] >= 0]
                    assert actual == read.get_aligned_pairs(matches_only=True)
                    assert read.to_string() == original
                    results.append(result)
                assert results[0] == results[1]


def terminal_cursor_inputs():
    rng = random.Random(19190926)
    bases = list(''.join(rng.choices('ACGT', k=500)))
    bases[120:122], bases[188:190], bases[201], bases[202:206] = 'GT', 'AG', 'C', 'AAAA'
    forward = ''.join(bases)
    rc = str.maketrans('ACGT', 'TGCA')
    genomes = {'cursorP': forward, 'cursorM': forward.translate(rc)[::-1]}
    register_genome_contigs(genomes)
    annotation = {('cursorP', 120, 190), ('cursorM', 310, 380)}
    header = pysam.AlignmentHeader.from_references(list(genomes), [500, 500])
    reads = []
    for reverse in (False, True):
        for soft in (0, 6):
            for hard in (0, 2, 17):
                ops = [(5, 3), (0, 20), (3, 70), (0, 16)]
                ops += [(4, soft)] if soft else []
                ops += [(5, hard)] if hard else []
                sequence = forward[100:120] + forward[190:206] + 'A' * soft
                if reverse:
                    ops, sequence = ops[::-1], sequence.translate(rc)[::-1]
                read = _read(header, f'cursor_{reverse}_{soft}_{hard}',
                             294 if reverse else 100, ops, sequence, reverse)
                read.reference_name = 'cursorM' if reverse else 'cursorP'
                reads.append(read)
    return genomes, annotation, header, reads


def run_terminal_cursor_roundtrip(directory):
    """Also runnable under plain Python; uses real correction and fixed TSV."""
    genomes, annotation, header, reads = terminal_cursor_inputs()
    rows = []
    for read in reads:
        row = bp.correct_read_3prime(read.__copy__(), genomes,
                                    annotated_junctions=annotation)[0]
        assert row['corrected_3prime'] == (298 if read.is_reverse else 201)
        rows.append(row)
    with pysam.AlignmentFile(str(directory / 'stock.bam'), 'wb', header=header) as out:
        for read in reads:
            out.write(read)
    with (directory / 'corrections.tsv').open('w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(CORRECTION_TSV_HEADER)
        for row in rows:
            writer.writerow(correction_result_to_tsv_row(row))
    bw.write_dual_bam(str(directory / 'stock.bam'), str(directory / 'corrections.tsv'),
                      str(directory / 'hard.bam'), str(directory / 'soft.bam'), genomes)
    outputs = {}
    for arm in ('hard', 'soft'):
        with pysam.AlignmentFile(str(directory / f'{arm}.bam'), 'rb') as inp:
            outputs[arm] = {read.query_name: read for read in inp}
        assert len(outputs[arm]) == len(reads)
        for original in reads:
            result = outputs[arm][original.query_name]
            _, _, soft, hard = original.query_name.split('_')
            soft, hard = int(soft), int(hard)
            assert result.infer_read_length() == original.infer_read_length()
            expected = [(5, 3), (0, 20), (3, 70), (0, 12)]
            if arm == 'hard':
                expected += [(5, hard + soft + 4)]
                wanted = original.query_sequence[soft + 4:] if original.is_reverse else original.query_sequence[:-(soft + 4)]
                qualities = original.query_qualities[soft + 4:] if original.is_reverse else original.query_qualities[:-(soft + 4)]
            else:
                expected += [(4, soft + 4)]
                expected += [(5, hard)] if hard else []
                wanted, qualities = original.query_sequence, original.query_qualities
            assert result.cigartuples == (expected[::-1] if original.is_reverse else expected)
            assert result.query_sequence == wanted and result.query_qualities == qualities
            assert result.get_tag('ZZ') == 'preserve'
            assert result.get_tag('cp') == (298 if result.is_reverse else 201)
            offset = soft + 4 if arm == 'hard' and original.is_reverse else 0
            stock_map = dict(original.get_aligned_pairs(matches_only=True))
            for q, p in result.get_aligned_pairs(matches_only=True):
                assert stock_map[q + offset] == p
                assert result.query_sequence[q] == genomes[result.reference_name][p]
            assert set(bw._n_op_intervals(result)) == set(bw._n_op_intervals(original))
    return reads, rows, outputs


def test_correction_fixed_tsv_and_writers_keep_exact_terminal_anchor(tmp_path):
    run_terminal_cursor_roundtrip(tmp_path)
