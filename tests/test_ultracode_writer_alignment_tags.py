"""Final writer NM/MD must describe final bases, verified by actual calmd."""
from array import array
import json
from pathlib import Path
import random

import pysam

from rectify.core.bam import bam_processor as bp, bam_writer as bw
from rectify.core.bam.alignment_tags import (
    calculate_nm_md, decode_original_sequence, finalize_alignment_tags, placement_state,
)
from rectify.core.bam.output import write_output_tsv
from rectify.utils.genome import register_genome_contigs


def _raises(call, text):
    try:
        call()
    except ValueError as exc:
        assert text in str(exc), str(exc)
    else:
        raise AssertionError('Expected ValueError: ' + text)


def _read(header, name, start, ops, sequence, reverse=False, contig=0):
    read = pysam.AlignedSegment(header)
    read.query_name, read.reference_id, read.reference_start = name, contig, start
    read.flag, read.mapping_quality = (16 if reverse else 0), 60
    read.cigartuples, read.query_sequence = ops, sequence
    if sequence is not None:
        read.query_qualities = [30 + i % 10 for i in range(len(sequence))]
    return read


def _fasta(directory, genome):
    path = directory / 'genome.fa'
    path.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in genome.items()))
    pysam.faidx(str(path))
    return path


def _bam(path, reads, header):
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for read in reads:
            out.write(read)


def _load(path):
    with pysam.AlignmentFile(str(path)) as bam:
        return {r.query_name: r for r in bam}


def _calmd(path, fasta, output, equal=False):
    # Never use pysam save_stdout: its input-file mutation was independently
    # witnessed during this audit. Direct binary capture preserves the input.
    original = path.read_bytes()
    args = (['-e'] if equal else []) + ['-b', str(path), str(fasta)]
    data = pysam.calmd(*args)
    assert isinstance(data, bytes)
    output.write_bytes(data)
    assert path.read_bytes() == original
    return _load(output)


def _geometry(read):
    return (read.reference_id, read.reference_start, read.cigarstring,
            read.query_sequence, read.qual)


def _modes(stock, tsv, directory, genome):
    bw.write_corrected_bam(str(stock), str(tsv), str(directory / 'hard.bam'), genome)
    bw.write_softclipped_bam(str(stock), str(tsv), str(directory / 'soft.bam'), genome)
    bw.write_dual_bam(str(stock), str(tsv), str(directory / 'dual_hard.bam'),
                      str(directory / 'dual_soft.bam'), genome)
    return ['hard', 'soft', 'dual_hard', 'dual_soft']


def terminal_tag_inputs():
    rng = random.Random(750921)
    text = list(''.join(rng.choice('ACGT') for _ in range(800)))
    text[100:181] = list(''.join(rng.choice('CGT') for _ in range(81)))
    text[180], text[181:193], text[193] = 'C', list('A' * 12), 'G'
    forward = ''.join(text)
    rc = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    genome = {'tagP': forward, 'tagM': rc(forward)}
    register_genome_contigs(genome)
    header = pysam.AlignmentHeader.from_references(list(genome), [800, 800])
    reads = []
    for reverse in (False, True):
        for kind in ('clean', 'genomic', 'mismatch', 'soft'):
            sequence, ops, end = forward[100:181], [(0, 81)], 181
            if kind in ('genomic', 'mismatch'):
                sequence += 'A' * 11 + ('C' if kind == 'mismatch' else 'A')
                ops, end = [(0, 93)], 193
            elif kind == 'soft':
                sequence += 'A' * 12
                ops += [(4, 12)]
            read = _read(header, f'{kind}_{reverse}', 800-end if reverse else 100,
                         ops[::-1] if reverse else ops, rc(sequence) if reverse else sequence,
                         reverse, int(reverse))
            for i, subtype in enumerate(('b', 'B', 'h', 'H', 'i', 'I', 'f')):
                read.set_tag('Z' + str(i), array(subtype, [1, 2, 3]))
            read.set_tag('XN', 1)
            read.set_tag('XO', 'fwd')
            reads.append(read)
    return genome, header, reads


def run_terminal_tag_roundtrip(directory):
    genome, header, reads = terminal_tag_inputs()
    fasta = _fasta(directory, genome)
    raw = directory / 'raw.bam'
    _bam(raw, reads, header)
    records = []
    for equal in (False, True):
        out = directory / ('equals' if equal else 'literal')
        out.mkdir()
        stock = out / 'stock.bam'
        inputs = _calmd(raw, fasta, stock, equal)
        rows = [bp.correct_read_3prime(r.__copy__(), genome, apply_atract=False,
                                      apply_3ss_rescue=False)[0] for r in inputs.values()]
        tsv = out / 'corrections.tsv'
        write_output_tsv(rows, str(tsv))
        (out / 'rows.json').write_text(json.dumps(rows, indent=2, default=str))
        for mode in _modes(stock, tsv, out, genome):
            emitted = _load(out / (mode + '.bam'))
            actual = _calmd(out / (mode + '.bam'), fasta, out / (mode + '.calmd.bam'))
            for name, read in emitted.items():
                oracle = actual[name]
                assert _geometry(read) == _geometry(oracle)
                assert (read.get_tag('NM'), read.get_tag('MD')) == (oracle.get_tag('NM'), oracle.get_tag('MD'))
                assert '=' not in read.query_sequence
                assert read.get_tag('XN') == 1 and read.get_tag('XO') == 'fwd'
                for i, subtype in enumerate(('b', 'B', 'h', 'H', 'i', 'I', 'f')):
                    assert read.get_tag('Z' + str(i)).typecode == subtype
                    assert list(read.get_tag('Z' + str(i))) == [1, 2, 3]
                if equal:
                    literal = _load(directory / 'literal' / (mode + '.bam'))[name]
                    assert read.to_string() == literal.to_string()
                records.append({'equal': equal, 'mode': mode, 'name': name,
                                'input_sam': inputs[name].to_string(),
                                'emitted_sam': read.to_string(), 'calmd_sam': oracle.to_string(),
                                'querymap': [(q, p) for q, p in read.get_aligned_pairs() if q is not None]})
    (directory / 'records.json').write_text(json.dumps(records, indent=2))
    return records


def test_terminal_all_modes_match_actual_calmd(tmp_path):
    assert len(run_terminal_tag_roundtrip(tmp_path)) == 64


def test_calmd_m_i_d_n_h_s_x_and_ambiguous_bases(tmp_path):
    rng = random.Random(750922)
    text = list(''.join(rng.choice('ACGT') for _ in range(1500)))
    text[100:117] = list('ACMGRSVTWYHKDBNNA')
    genome = {'calmd': ''.join(text)}
    header = pysam.AlignmentHeader.from_references(['calmd'], [1500])
    shapes = [[(0, 17)], [(5, 7), (4, 3), (0, 17), (4, 2), (5, 9)],
              [(0, 10), (1, 3), (0, 7)], [(0, 8), (2, 3), (0, 9)],
              [(0, 8), (3, 500), (0, 9)],
              [(0, 5), (2, 2), (2, 1), (3, 50), (7, 4), (1, 2), (8, 6)],
              [(0, 5), (3, 80), (0, 6), (3, 70), (0, 7)]]
    reads = []
    for reverse in (False, True):
        for i, ops in enumerate(shapes):
            sequence, position = '', 100
            for op, n in ops:
                if op in (0, 7, 8):
                    sequence += genome['calmd'][position:position+n]
                elif op in (1, 4):
                    sequence += 'C' * n
                if op in (0, 2, 3, 7, 8):
                    position += n
            # Actual bases govern tags, regardless of the supplied =/X op label.
            sequence = ('T' if sequence[0] != 'T' else 'A') + sequence[1:]
            reads.append(_read(header, f'cigar_{i}_{reverse}', 100, ops, sequence, reverse))
    fasta = _fasta(tmp_path, genome)
    raw = tmp_path / 'raw.bam'
    _bam(raw, reads, header)
    actual = _calmd(raw, fasta, tmp_path / 'calmd.bam')
    for read in reads:
        assert calculate_nm_md(read, genome) == (actual[read.query_name].get_tag('NM'), actual[read.query_name].get_tag('MD'))


def test_large_n_is_not_walked_or_counted():
    class Reference:
        def __init__(self): self.lookups = 0
        def __len__(self): return 100000100
        def __getitem__(self, index):
            assert isinstance(index, int)
            self.lookups += 1
            return 'A'
    reference = Reference()
    header = pysam.AlignmentHeader.from_references(['longN'], [100000100])
    read = _read(header, 'longN', 10, [(0, 10), (3, 50000000), (0, 12)], 'A' * 22)
    assert calculate_nm_md(read, {'longN': reference}) == (0, '22')
    assert reference.lookups == 22


def test_original_equals_decode_precedes_placement_change():
    header = pysam.AlignmentHeader.from_references(['frame'], [80])
    genome = {'frame': 'A' * 20 + 'CGTC' + 'T' * 56}
    read = _read(header, 'frame', 20, [(5, 7), (0, 4), (5, 9)], '====')
    quality = read.qual
    assert decode_original_sequence(read, genome)
    before = placement_state(read)
    read.reference_start = 30
    assert finalize_alignment_tags(read, before, genome)
    assert read.query_sequence == 'CGTC' and read.qual == quality
    assert read.get_tag('NM') == 3 and read.get_tag('MD') == '0T0T1T0'


def test_unresolved_equals_refuses_without_mutation():
    header = pysam.AlignmentHeader.from_references(['guard'], [80])
    for ops, sequence, genome, start in [
        ([(0, 4)], '====', None, 10), ([(0, 4)], '====', {'other': 'A' * 80}, 10),
        ([(0, 4)], '====', {'guard': 'A' * 12}, 10),
        ([(4, 1), (0, 3)], '=AAA', {'guard': 'A' * 80}, 10),
        ([(1, 1), (0, 3)], '=AAA', {'guard': 'A' * 80}, 10),
        ([(0, 4)], '====', {'guard': 'A' * 80}, -1),
        ([(0, 4)], '====', {'guard': '=' * 80}, 10),
    ]:
        read = _read(header, 'guard', start, ops, sequence)
        before = read.to_string()
        _raises(lambda: decode_original_sequence(read, genome), 'Cannot decode original')
        assert read.to_string() == before


def test_unavailable_reference_clears_only_obsolete_tags():
    header = pysam.AlignmentHeader.from_references(['unknown'], [100])
    for genome in (None, {}, {'unknown': 'A' * 11}):
        read = _read(header, 'unknown', 10, [(0, 5)], 'ACGTA')
        for tag in ('NM', 'MD', 'AS', 'ms', 'cs', 'de', 'dv', 'UQ'):
            for _ in range(2): read.set_tag(tag, 'historical marker')
        read.set_tag('ZB', array('I', [1, 4000000000]))
        read.set_tag('XO', 'rev'); read.set_tag('XN', 1)
        before = placement_state(read)
        read.cigartuples = [(4, 1), (0, 4)]
        read.reference_start += 1
        assert finalize_alignment_tags(read, before, genome)
        assert not any(read.has_tag(t) for t in ('NM', 'MD', 'AS', 'ms', 'cs', 'de', 'dv', 'UQ'))
        assert read.get_tag('ZB') == array('I', [1, 4000000000])
        assert read.get_tag('XO') == 'rev' and read.get_tag('XN') == 1


def test_missing_sequence_invalidates_without_fabrication():
    header = pysam.AlignmentHeader.from_references(['unknown'], [100])
    read = _read(header, 'noSEQ', 10, [(0, 5)], None)
    read.set_tag('NM', 3); read.set_tag('MD', '2A2')
    before = placement_state(read)
    read.reference_start += 1
    assert finalize_alignment_tags(read, before, {'unknown': 'A' * 100})
    assert read.query_sequence is None and not read.has_tag('NM') and not read.has_tag('MD')


def test_changed_duplicates_become_single_fresh_nm_md():
    header = pysam.AlignmentHeader.from_references(['dup'], [100])
    read = _read(header, 'dup', 10, [(0, 5)], 'ACGTA')
    for tag in ('NM', 'MD', 'AS', 'ms', 'cs', 'de', 'dv', 'UQ'):
        for _ in range(2): read.set_tag(tag, 'historical marker', replace=False)
    before = placement_state(read)
    read.cigartuples, read.reference_start = [(4, 1), (0, 4)], 11
    assert finalize_alignment_tags(read, before, {'dup': 'A' * 100})
    tags = read.get_tags(with_value_type=True)
    assert [(v, t) for k, v, t in tags if k == 'NM'] == [(3, 'i')]
    assert [(v, t) for k, v, t in tags if k == 'MD'] == [('0A0A0A1', 'Z')]
    assert not any(k in ('AS', 'ms', 'cs', 'de', 'dv', 'UQ') for k, _, _ in tags)


def test_unchanged_records_keep_complete_sam():
    header = pysam.AlignmentHeader.from_references(['same'], [100])
    read = _read(header, 'same', 10, [(0, 5)], 'ACGTA')
    read.set_tag('NM', 4); read.set_tag('MD', '0A0A0A0A1')
    read.set_tag('AS', 17); read.set_tag('ZZ', array('h', [-2, 300]))
    before = read.to_string()
    assert not finalize_alignment_tags(read, placement_state(read), {'same': 'A' * 100})
    assert not bw.apply_corrected_edits_to_read(read, None, None)
    assert read.to_string() == before


def test_atomic_all_modes_preserve_previous_outputs_on_late_invalid(tmp_path):
    genome, header, reads = terminal_tag_inputs()
    valid = reads[0]
    invalid = valid.__copy__(); invalid.query_name = 'late_invalid'
    invalid.query_sequence = '=' + invalid.query_sequence[1:]
    invalid.cigartuples = [(4, 1), (0, 80)]
    raw = tmp_path / 'raw.bam'; _bam(raw, [valid, invalid], header)
    tsv = tmp_path / 'rows.tsv'
    tsv.write_text('read_id\tcorrected_3prime\tstrand\ttail_correction_enabled\n'
                   f'{valid.query_name}\t180\t+\t1\nlate_invalid\t180\t+\t1\n')
    for mode in ('hard', 'soft', 'dual'):
        paths = [tmp_path / (mode + '.bam')]
        if mode == 'dual': paths.append(tmp_path / 'dual_second.bam')
        for i, path in enumerate(paths): path.write_bytes(f'previous-{i}'.encode())
        if mode == 'dual':
            call = lambda: bw.write_dual_bam(str(raw), str(tsv), *(str(p) for p in paths), genome)
        else:
            writer = bw.write_corrected_bam if mode == 'hard' else bw.write_softclipped_bam
            call = lambda: writer(str(raw), str(tsv), str(paths[0]), genome)
        _raises(call, 'equals in insertion or soft clip')
        assert [p.read_bytes() for p in paths] == [f'previous-{i}'.encode() for i in range(len(paths))]
        assert not list(tmp_path.glob('.*.tmp.bam'))


def test_live_2f_and_refused_reanchor_have_final_tags(tmp_path):
    from tests.test_ultracode_microexon_placement import placement_inputs, placement_read, write_rows
    from rectify.core.splice import microexon as mx
    genome, _, annotation, header = placement_inputs()
    previous = mx.microexon_index()
    mx.set_microexon_index({})
    try:
        fasta = _fasta(tmp_path, genome)
        raw = tmp_path / 'raw.bam'
        reads = [placement_read(kind, reverse, 'literal', genome, header)[0]
                 for kind in ('overlap_success', 'reanchor_refused') for reverse in (False, True)]
        _bam(raw, reads, header)
        for equal in (False, True):
            out = tmp_path / ('equals' if equal else 'literal'); out.mkdir()
            stock = out / 'stock.bam'
            inputs = _calmd(raw, fasta, stock, equal)
            rows = [bp.correct_read_3prime(r.__copy__(), genome, annotated_junctions=annotation)[0]
                    for r in inputs.values()]
            assert all(not r['station_b_applied'] for r in rows)
            assert sum(bool(r['five_prime_rescued']) for r in rows) == 2
            assert sum(r['five_prime_rescue_refused'] == 'extend_refused' for r in rows) == 2
            tsv = out / 'rows.tsv'; write_rows(tsv, rows)
            for mode in _modes(stock, tsv, out, genome):
                emitted = _load(out / (mode + '.bam'))
                actual = _calmd(out / (mode + '.bam'), fasta, out / (mode + '.calmd.bam'))
                for name, read in emitted.items():
                    oracle = actual[name]
                    assert _geometry(read) == _geometry(oracle)
                    assert read.get_tag('NM') == oracle.get_tag('NM') == 0
                    assert read.get_tag('MD') == oracle.get_tag('MD')
                    assert not read.has_tag('Xb')
                    if equal:
                        assert read.to_string() == _load(tmp_path / 'literal' / (mode + '.bam'))[name].to_string()
    finally:
        mx.set_microexon_index(previous)
