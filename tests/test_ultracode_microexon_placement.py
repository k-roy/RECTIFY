"""Actual 2F decisions, fixed TSV and all writers must agree on microexons."""
from array import array
import csv
import random

import pysam

from rectify.core.bam import bam_processor as bp
from rectify.core.bam import bam_writer as bw
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.splice import microexon as mx
from rectify.core.splice.microexon_provenance import load_microexon_calls
from rectify.utils.genome import register_genome_contigs

RC = str.maketrans('ACGT', 'TGCA')
KINDS = ('overlap_success', 'overlap_refused', 'reanchor_refused',
         'disjoint', 'overlap_and_disjoint')


def placement_inputs():
    rng = random.Random(191926)
    bases = list(''.join(rng.choices('ACGT', k=1200)))
    for start, end in [(60, 460), (60, 200), (206, 460),
                       (500, 900), (500, 700), (706, 900)]:
        bases[start:start + 2] = 'GT'
        bases[end - 2:end] = 'AG'
    bases[559] = bases[999] = 'C'  # 3' tail handling is inert in this fixture.
    normal = ''.join(bases)
    bases[200:206] = bases[54:60]
    repeated = ''.join(bases)
    genomes = {}
    for key, genome in [('repeat', repeated), ('normal', normal)]:
        genomes[key + 'P'] = genome
        genomes[key + 'M'] = genome.translate(RC)[::-1]
    register_genome_contigs(genomes)
    index = {c: ([(200, 206), (700, 706)] if c.endswith('P')
                 else [(494, 500), (994, 1000)]) for c in genomes}
    annotated = {(c, s, e) for c in genomes
                 for s, e in ([(60, 460), (500, 900)] if c.endswith('P')
                              else [(740, 1140), (300, 700)])}
    header = pysam.AlignmentHeader.from_references(list(genomes), [1200] * 4)
    return genomes, index, annotated, header


def placement_read(kind, reverse, encoding, genomes, header):
    repeated = kind in ('overlap_success', 'overlap_and_disjoint')
    prefix = 'repeat' if repeated else 'normal'
    g = genomes[prefix + 'P']
    if kind == 'disjoint':
        start = 460
        ops = [(4, 12), (0, 40), (3, 400), (1, 6), (0, 100)]
        seq = g[48:60] + g[460:500] + g[700:706] + g[900:1000]
        expected_start = 48
        expected_ops = [(0, 12), (3, 400), (0, 40), (3, 200),
                        (0, 6), (3, 194), (0, 100)]
    elif kind == 'overlap_refused':
        start = expected_start = 36
        ops = [(2, 12), (0, 12), (3, 400), (1, 6), (0, 100)]
        seq = g[48:60] + g[200:206] + g[460:560]
        expected_ops = [(2, 12), (0, 12), (3, 140), (0, 6), (3, 254), (0, 100)]
    else:
        start = 36
        ops = [(0, 12), (2, 12), (3, 400), (1, 6)]
        seq = (g[42:54] if repeated else g[48:60]) + g[200:206]
        if kind == 'overlap_and_disjoint':
            ops += [(0, 40), (3, 400), (1, 6), (0, 100)]
            seq += g[460:500] + g[700:706] + g[900:1000]
            expected_start = 42
            expected_ops = [(0, 18), (3, 400), (0, 40), (3, 200),
                            (0, 6), (3, 194), (0, 100)]
        else:
            ops += [(0, 100)]
            seq += g[460:560]
            expected_start = 42 if repeated else 460
            expected_ops = [(0, 18), (3, 400), (0, 100)] if repeated else [(4, 18), (0, 100)]
    if reverse:
        start = 1200 - start - sum(n for op, n in ops if op in (0, 2, 3, 7, 8))
        ops = list(reversed(ops))
        seq = seq.translate(RC)[::-1]
        expected_start = 1200 - expected_start - sum(
            n for op, n in expected_ops if op in (0, 2, 3, 7, 8))
        expected_ops = list(reversed(expected_ops))
    read = pysam.AlignedSegment(header)
    read.query_name = f'{kind}_{"minus" if reverse else "plus"}_{encoding}'
    read.reference_name = prefix + ('M' if reverse else 'P')
    read.reference_start = start
    read.flag = 16 if reverse else 0
    read.mapping_quality = 60
    read.cigartuples = ops
    read.query_sequence = seq
    read.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    read.set_tag('ZA', array('h', [1, -2, 300]))
    read.set_tag('ZZ', 'preserved')
    read.set_tag('XB', '2/1')  # cDNA tag namespace is independent of Xb.
    expected = read.__copy__()
    expected.reference_start = expected_start
    expected.cigartuples = expected_ops
    if encoding == 'equals':
        encoded = list(seq)
        for q, p in read.get_aligned_pairs(matches_only=True):
            if encoded[q] == genomes[read.reference_name][p]:
                encoded[q] = '='
        qualities = read.query_qualities
        read.query_sequence = ''.join(encoded)
        read.query_qualities = qualities
    return read, expected


def write_rows(path, rows):
    with path.open('w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(CORRECTION_TSV_HEADER)
        for row in rows:
            writer.writerow(correction_result_to_tsv_row(row))


def run_placement_roundtrip(directory):
    """Also callable by a plain-Python audit witness; no fixture or score mocks."""
    genomes, index, annotated, header = placement_inputs()
    previous = mx.microexon_index()
    previous_species_tier = mx._SPECIES_MAX_3SS_TIER
    mx.set_species('homo_sapiens')
    mx.set_microexon_index(index)
    records, expected, rows = {}, {}, {}
    try:
        for kind in KINDS:
            for reverse in (False, True):
                for encoding in ('literal', 'equals'):
                    read, want = placement_read(kind, reverse, encoding, genomes, header)
                    # Every original read genuinely contains a Station B substrate.
                    decoded = read.__copy__()
                    bw._decode_eq_seq_inplace(decoded, genomes)
                    calls = mx.recover_all_microexons(decoded, genomes[read.reference_name],
                                                    '-' if reverse else '+', index)
                    assert len(calls) == (2 if kind == 'overlap_and_disjoint' else 1)
                    row = bp.correct_read_3prime(read.__copy__(), genomes,
                                               annotated_junctions=annotated)[0]
                    assert row['five_prime_rescued'] == (kind in ('overlap_success', 'disjoint', 'overlap_and_disjoint'))
                    assert bool(row['station_b_applied']) == (kind in ('overlap_refused', 'disjoint', 'overlap_and_disjoint'))
                    if kind in ('overlap_refused', 'reanchor_refused'):
                        assert row['five_prime_rescue_refused'] == 'extend_refused'
                    else:
                        assert row['five_prime_rescue_refused'] == ''
                    records[read.query_name], expected[read.query_name], rows[read.query_name] = read, want, row
        stock = directory / 'stock.bam'
        tsv = directory / 'corrections.tsv'
        with pysam.AlignmentFile(str(stock), 'wb', header=header) as f:
            for read in records.values():
                f.write(read)
        write_rows(tsv, rows.values())
        bw.write_corrected_bam(str(stock), str(tsv), str(directory / 'hard.bam'), genomes)
        bw.write_softclipped_bam(str(stock), str(tsv), str(directory / 'soft.bam'), genomes)
        bw.write_dual_bam(str(stock), str(tsv), str(directory / 'dual_hard.bam'),
                          str(directory / 'dual_soft.bam'), genomes)
        outputs = {}
        for arm in ('hard', 'soft', 'dual_hard', 'dual_soft'):
            with pysam.AlignmentFile(str(directory / (arm + '.bam')), 'rb') as f:
                outputs[arm] = {r.query_name: r for r in f}
            assert set(outputs[arm]) == set(records)
            for name, read in outputs[arm].items():
                want, row = expected[name], rows[name]
                assert (read.reference_start, read.cigartuples, read.query_sequence) == (
                    want.reference_start, want.cigartuples, want.query_sequence), (arm, name)
                assert read.query_qualities == want.query_qualities
                assert read.get_aligned_pairs() == want.get_aligned_pairs(), (arm, name)
                assert set(bw._n_op_intervals(read)) == {tuple(j) for j in row['junctions']}
                assert read.get_tag('ZA') == array('h', [1, -2, 300])
                assert read.get_tag('ZZ') == 'preserved'
                assert read.get_tag('XB') == '2/1'
                assert read.has_tag('Xb') == bool(row['station_b_applied'])
                if read.has_tag('Xb'):
                    calls = load_microexon_calls(read)
                    assert calls and len(calls) == 1
                assert read.to_string() == outputs['hard'][name].to_string(), (arm, name)
        return genomes, records, rows, outputs
    finally:
        mx.set_microexon_index(previous)
        mx._SPECIES_MAX_3SS_TIER = previous_species_tier


def test_live_2f_station_b_all_writer_modes(tmp_path, monkeypatch):
    monkeypatch.setenv('RECTIFY_STATION_B', 'apply')
    run_placement_roundtrip(tmp_path)


def test_report_mode_preserves_call_without_drawing(tmp_path, monkeypatch):
    monkeypatch.setenv('RECTIFY_STATION_B', 'report')
    genomes, index, annotated, header = placement_inputs()
    previous = mx.microexon_index()
    mx.set_microexon_index(index)
    try:
        read, _ = placement_read('disjoint', False, 'equals', genomes, header)
        row = bp.correct_read_3prime(read.__copy__(), genomes, annotated_junctions=annotated)[0]
        assert row['five_prime_rescued'] and row['station_b_microexons']
        assert not row['station_b_applied']
        correction = tmp_path / 'report.tsv'
        write_rows(correction, [row])
        loaded = bw._load_corrections_from_tsv(str(correction))[read.query_name]
        bw.apply_corrected_edits_to_read(read, loaded, genomes)
        assert not read.has_tag('Xb')
        assert set(bw._n_op_intervals(read)) == {tuple(j) for j in row['junctions']}
    finally:
        mx.set_microexon_index(previous)


def test_only_live_microexon_draw_invalidates_alignment_tags(tmp_path, monkeypatch):
    monkeypatch.setenv('RECTIFY_STATION_B', 'apply')
    genomes, index, annotated, header = placement_inputs()
    previous = mx.microexon_index()
    mx.set_microexon_index(index)
    tags = {'MD': '24', 'NM': 6, 'AS': 42, 'ms': 42,
            'cs': ':12+acggat:100', 'de': .1, 'dv': .1, 'UQ': 12}
    fasta = tmp_path / 'tag_reference.fa'
    fasta.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in genomes.items()))
    pysam.faidx(str(fasta))
    try:
        for reverse in (False, True):
            for encoding in ('literal', 'equals'):
                read, _ = placement_read('overlap_refused', reverse, encoding, genomes, header)
                row = bp.correct_read_3prime(read.__copy__(), genomes,
                                           annotated_junctions=annotated)[0]
                assert row['station_b_applied']
                tsv = tmp_path / 'live.tsv'
                write_rows(tsv, [row])
                correction = bw._load_corrections_from_tsv(str(tsv))[read.query_name]
                for tag, value in tags.items():
                    read.set_tag(tag, value)
                    # Duplicate obsolete fields must not leave a stale second
                    # occurrence behind when pysam deletes only the first.
                    read.set_tag(tag, value, replace=False)
                assert bw.apply_corrected_edits_to_read(read, correction, genomes)
                # ISSUE075: the outer writer now recomputes final NM/MD after
                # the B primitive invalidates obsolete fields. Every other
                # derived score remains absent; duplicate old tags cannot leak.
                assert not any(read.has_tag(tag) for tag in tags if tag not in ('NM', 'MD'))
                final_bam = tmp_path / 'tag_final.bam'
                with pysam.AlignmentFile(str(final_bam), 'wb', header=header) as out:
                    out.write(read)
                oracle_bam = tmp_path / 'tag_oracle.bam'
                oracle_bam.write_bytes(pysam.calmd('-b', str(final_bam), str(fasta)))
                with pysam.AlignmentFile(str(oracle_bam)) as bam:
                    oracle = next(bam)
                for tag in ('NM', 'MD'):
                    assert [v for key, v in read.get_tags() if key == tag] == [oracle.get_tag(tag)]
                assert read.cigarstring == oracle.cigarstring
                assert read.query_sequence == oracle.query_sequence and read.qual == oracle.qual
                assert read.has_tag('Xb') and read.get_tag('XB') == '2/1'
                assert read.get_tag('ZA') == array('h', [1, -2, 300])
                # A replay cannot redraw this consumed substrate. Newly supplied
                # tags remain untouched when the primitive makes no placement.
                for tag, value in tags.items():
                    read.set_tag(tag, value)
                before = read.to_string()
                assert not bw.apply_station_b_microexons(read, correction)
                assert read.to_string() == before
                correction['station_b_applied'] = 0
                assert not bw.apply_station_b_microexons(read, correction)
                assert read.to_string() == before
    finally:
        mx.set_microexon_index(previous)
