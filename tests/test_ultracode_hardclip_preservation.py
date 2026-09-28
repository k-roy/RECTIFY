"""Original molecule extent must survive2F H/S and reanchor geometry."""
from array import array
import importlib.util
from pathlib import Path

import pysam

from rectify.core.bam import read_edits as e, bam_writer as w, bam_processor as bp
from rectify.core.splice import microexon as mx
from rectify.core.splice.splice_aware_5prime import _get_5prime_softclip_len


def read_at(ops, start=100, reverse=False, query=None):
    header = pysam.AlignmentHeader.from_references(['clipT'], [1000])
    read = pysam.AlignedSegment(header)
    read.query_name = 'H_geometry'
    read.reference_id = 0
    read.flag = 16 if reverse else 0
    read.reference_start = start
    read.cigartuples = ops
    read.query_sequence = query or 'A' * sum(n for o, n in ops if o in (0, 1, 4, 7, 8))
    read.query_qualities = array('B', [35] * len(read.query_sequence))
    return read


def terminal_h(read):
    ops = read.cigartuples
    left = sum(n for o, n in ops[:next((i for i, x in enumerate(ops) if x[0] != 5), len(ops))])
    tail = list(reversed(ops))
    right = sum(n for o, n in tail[:next((i for i, x in enumerate(tail) if x[0] != 5), len(tail))])
    return left, right


def without_h(ops):
    return [(o, n) for o, n in ops if o != 5]


def test_slow_fast_reanchor_h_m_i_and_h_s_both_ends():
    for reverse in (False, True):
        for edge in ([(0, 3), (1, 2)], [(4, 5)]):
            ops = [(5, 7)] + edge + [(0, 20), (5, 3)]
            query = 'C' * 5 + 'A' * 20
            if reverse:
                ops, query = ops[::-1], query[::-1]
            read = read_at(ops, reverse=reverse, query=query)
            slow, fast = read.__copy__(), read.__copy__()
            assert e.reanchor_5prime_for_rescue(slow, {'clipT': 'A' * 1000})
            assert e._apply_reanchor_from_clip_len(fast, 5)
            assert slow.to_string() == fast.to_string()
            assert terminal_h(slow) == terminal_h(read)
            assert _get_5prime_softclip_len(slow, '-' if reverse else '+') == 5
            assert slow.query_sequence == read.query_sequence and slow.query_qualities == read.query_qualities
            no_h = read.__copy__()
            no_h.cigartuples = without_h(no_h.cigartuples)
            assert e.reanchor_5prime_for_rescue(no_h, {'clipT': 'A' * 1000})
            assert without_h(slow.cigartuples) == no_h.cigartuples
            assert slow.reference_start == no_h.reference_start
            assert slow.get_aligned_pairs() == no_h.get_aligned_pairs()


def test_h_s_projection_extend_and_refusal():
    for reverse in (False, True):
        strand = '-' if reverse else '+'
        ops = [(5, 7), (4, 12), (0, 30), (5, 3)]
        if reverse: ops = ops[::-1]
        original = read_at(ops, reverse=reverse)
        assert _get_5prime_softclip_len(original, strand) == 12
        expected_edge = original.reference_end if reverse else original.reference_start
        assert e.projected_5prime_rescue_intron_edge(original, 12, strand) == expected_edge
        target = original.reference_end + 100 if reverse else original.reference_start - 21
        read = original.__copy__()
        assert e.extend_read_5prime_for_junction_rescue(read, target, 12, strand, '12M')
        assert terminal_h(read) == terminal_h(original)
        assert read.query_sequence == original.query_sequence and read.query_qualities == original.query_qualities
        for cigar in ('11M', '13M'):
            refused = original.__copy__()
            before = refused.to_string()
            assert not e.extend_read_5prime_for_junction_rescue(refused, target, 12, strand, cigar)
            assert refused.to_string() == before


def test_h_m_i_reroute_preserves_molecule_extent():
    for reverse in (False, True):
        strand = '-' if reverse else '+'
        ops = [(5, 7), (0, 4), (1, 2), (0, 30), (5, 3)]
        if reverse: ops = ops[::-1]
        read = read_at(ops, reverse=reverse)
        original = read.__copy__()
        boundary = read.reference_end - 4 if reverse else read.reference_start + 4
        target = boundary + 100 if reverse else boundary - 101
        assert e.reroute_intronic_tail_5prime_via_junction(read, boundary, target, '6M', strand)
        assert terminal_h(read) == terminal_h(original)
        assert read.query_sequence == original.query_sequence and read.query_qualities == original.query_qualities
        plain = original.__copy__()
        plain.cigartuples = without_h(plain.cigartuples)
        assert e.reroute_intronic_tail_5prime_via_junction(plain, boundary, target, '6M', strand)
        assert without_h(read.cigartuples) == plain.cigartuples
        assert read.get_aligned_pairs() == plain.get_aligned_pairs()
        refused = original.__copy__()
        before = refused.to_string()
        assert not e.reroute_intronic_tail_5prime_via_junction(refused, boundary, target, '5M', strand)
        assert refused.to_string() == before


def run_hardclip_correction_matrix(directory):
    """Production correction -> fixed TSV -> all writers; no score mocks."""
    spec = importlib.util.spec_from_file_location('microexon_clip_inputs', Path(__file__).with_name('test_ultracode_microexon_placement.py'))
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    genomes, index, annotation, header = fixture.placement_inputs()
    previous = mx.microexon_index()
    mx.set_microexon_index(index)
    records, rows, expected = [], [], {}
    try:
        for kind in fixture.KINDS:
            for reverse in (False, True):
                for encoding in ('literal', 'equals'):
                    read, _ = fixture.placement_read(kind, reverse, encoding, genomes, header)
                    plain_row = bp.correct_read_3prime(read.__copy__(), genomes, annotated_junctions=annotation)[0]
                    plain = read.__copy__()
                    read.query_name = 'H_' + read.query_name
                    left, right = (3, 7) if reverse else (7, 3)
                    read.cigartuples = [(5, left)] + read.cigartuples + [(5, right)]
                    row = bp.correct_read_3prime(read.__copy__(), genomes, annotated_junctions=annotation)[0]
                    for key in ('five_prime_rescued', 'five_prime_position', 'five_prime_exon_cigar',
                                'five_prime_rescue_refused', 'reanchor_clip_len', 'station_b_applied'):
                        assert row.get(key) == plain_row.get(key), (read.query_name, key, row.get(key), plain_row.get(key))
                    records.append(read)
                    rows.append(row)
                    expected[read.query_name] = (plain, plain_row)
        source, tsv = directory / 'input.bam', directory / 'corrections.tsv'
        with pysam.AlignmentFile(str(source), 'wb', header=header) as out:
            for read in records: out.write(read)
        fixture.write_rows(tsv, rows)
        w.write_corrected_bam(str(source), str(tsv), str(directory / 'hard.bam'), genomes)
        w.write_softclipped_bam(str(source), str(tsv), str(directory / 'soft.bam'), genomes)
        w.write_dual_bam(str(source), str(tsv), str(directory / 'dual_hard.bam'), str(directory / 'dual_soft.bam'), genomes)
        # Paired no-H controls go through the same production TSV/writer path.
        plain_source, plain_tsv = directory / 'plain.bam', directory / 'plain.tsv'
        with pysam.AlignmentFile(str(plain_source), 'wb', header=header) as out:
            for plain, _ in expected.values(): out.write(plain)
        fixture.write_rows(plain_tsv, [row for _, row in expected.values()])
        w.write_dual_bam(str(plain_source), str(plain_tsv), str(directory / 'plain_hard.bam'),
                         str(directory / 'plain_soft.bam'), genomes)
        plain_outputs = {}
        for arm in ('hard', 'soft'):
            with pysam.AlignmentFile(str(directory / ('plain_' + arm + '.bam')), 'rb') as f:
                plain_outputs[arm] = {r.query_name:r for r in f}
        outputs = {}
        original = {r.query_name:r for r in records}
        for arm in ('hard', 'soft', 'dual_hard', 'dual_soft'):
            with pysam.AlignmentFile(str(directory / (arm + '.bam')), 'rb') as f:outputs[arm] = {r.query_name:r for r in f}
            for name, out in outputs[arm].items():
                before = original[name]
                assert terminal_h(out) == terminal_h(before), (name, arm, out.cigarstring)
                plain = plain_outputs['hard' if 'hard' in arm else 'soft'][expected[name][0].query_name]
                assert without_h(out.cigartuples) == plain.cigartuples
                assert out.reference_start == plain.reference_start
                assert out.get_aligned_pairs() == plain.get_aligned_pairs()
                left_h = terminal_h(out)[0]
                full_map = {left_h + q:p for q,p in out.get_aligned_pairs() if q is not None}
                expected_map = {terminal_h(before)[0] + q:p for q,p in plain.get_aligned_pairs() if q is not None}
                assert full_map == expected_map
                literal = before.__copy__()
                w._decode_eq_seq_inplace(literal, genomes)
                assert out.query_sequence == literal.query_sequence and out.query_qualities == literal.query_qualities
                assert len(out.query_sequence) + sum(terminal_h(out)) == len(before.query_sequence) + sum(terminal_h(before))
        return genomes, original, rows, outputs
    finally:
        mx.set_microexon_index(previous)


def test_production_2f_fixed_tsv_four_writers_keep_h(tmp_path, monkeypatch):
    monkeypatch.setenv('RECTIFY_STATION_B', 'apply')
    run_hardclip_correction_matrix(tmp_path)


def test_h_does_not_hide_repeat_refusal_or_clip_origin():
    from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation, clip_origin
    for reverse in (False, True):
        strand = '-' if reverse else '+'
        repeat = 'AAG' * 15
        ops, sequence = [(4, 45), (0, 30)], repeat + 'A' * 30
        if reverse:
            ops, sequence = ops[::-1], sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        plain = read_at(ops, reverse=reverse, query=sequence)
        hard = plain.__copy__()
        hard.cigartuples = [(5, 7)] + hard.cigartuples + [(5, 3)]
        for read in (plain, hard):
            before = read.to_string()
            result = rescue_3ss_truncation(read, {'clipT': 'A' * 1000}, set(), strand)
            assert result.get('repeat_expansion') and not result['rescued']
            assert read.to_string() == before
        ops, sequence = [(4, 12), (0, 30)], 'ACGTTCGAGTCC' + 'A' * 30
        if reverse:
            ops, sequence = ops[::-1], sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        plain = read_at(ops, reverse=reverse, query=sequence)
        hard = plain.__copy__()
        hard.cigartuples = [(5, 7)] + hard.cigartuples + [(5, 3)]
        a = clip_origin(plain, strand, 'A' * 1000, ('clipT', 30, 80), 20, True)
        b = clip_origin(hard, strand, 'A' * 1000, ('clipT', 30, 80), 20, True)
        assert a == b and a[0] != 'none'
