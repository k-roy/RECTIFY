"""Consensus evidence and placements cannot depend on SAM '=' spelling."""
from array import array
import copy
from pathlib import Path
import random

import pysam

from rectify.core.consensus.consensus import run_consensus_selection
from rectify.core.consensus.extract import extract_alignment_info
from rectify.core.consensus.scoring import (
    _get_effective_5prime_clip, _get_effective_3prime_clip,
    _count_junction_proximity_errors, score_alignment,
)
from rectify.core.consensus.chimeric_consensus import select_best_chimeric, build_chimeric_read
from rectify.core.consensus.sequence import decoded_alignment_copy
from rectify.utils.genome import register_genome_contigs


def make_family(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rng = random.Random(2009066)
    bases = list(''.join(rng.choices('ACGT', k=1600)))
    bases[300:302], bases[394:396], bases[474:476] = 'CG', 'GC', 'GC'
    comp = str.maketrans('ACGT', 'TGCA')
    bases[900:996] = ''.join(bases[300:396]).translate(comp)
    bases[1020:1076] = ''.join(bases[420:476]).translate(comp)
    for start, end in ((340, 420), (940, 1020)):
        bases[start:start + 2], bases[end - 2:end] = 'GT', 'AG'
    forward = ''.join(bases)
    genome = {'eqP': forward, 'eqM': forward.translate(comp)[::-1]}
    register_genome_contigs(genome)
    header = pysam.AlignmentHeader.from_references(list(genome), [1600, 1600])
    reference = directory / 'genome.fa'
    reference.write_text(''.join('>' + chrom + '\n' + seq + '\n' for chrom, seq in genome.items()))
    pysam.faidx(str(reference))
    paths = {'literal': {}, 'calmd': {}}
    expected, inputs = {}, {}
    for arm, start in (('minimap2', 900), ('uLTRA', 300)):
        original = directory / (arm + '.original.bam')
        with pysam.AlignmentFile(str(original), 'wb', header=header) as out:
            for reverse in (False, True):
                for kind in ('plain', 'HSIN'):
                    for quality in (False, True):
                        name = f'{kind}_R{int(reverse)}_Q{int(quality)}'
                        if kind == 'plain':
                            ops, query = [(0, 96)], forward[300:396]
                        else:
                            ops = [(5, 4), (4, 3), (0, 40), (1, 2), (3, 80), (0, 56), (4, 4), (5, 7)]
                            query = 'GCC' + forward[300:340] + 'CG' + forward[420:476] + 'CCGC'
                        span = sum(n for op, n in ops if op in (0, 2, 3, 7, 8))
                        r = pysam.AlignedSegment(header)
                        r.query_name, r.reference_name = name, 'eqM' if reverse else 'eqP'
                        r.reference_start = 1600 - start - span if reverse else start
                        r.flag, r.mapping_quality = (16 if reverse else 0), 55
                        r.cigartuples = ops[::-1] if reverse else ops
                        r.query_sequence = query.translate(comp)[::-1] if reverse else query
                        if quality:
                            r.query_qualities = array('B', [24 + i % 12 for i in range(len(query))])
                        r.set_tag('ZZ', array('h', [-2, 71, 800]))
                        r.set_tag('XO', 'fwd' if reverse else 'rev')
                        r.set_tag('XN', 1)
                        r.set_tag('RN', len(expected) if name not in expected else list(expected).index(name))
                        expected[name] = r.query_sequence
                        out.write(r)
        encoded = directory / (arm + '.calmd.bam')
        encoded.write_bytes(pysam.calmd('-e', '-b', str(original), str(reference)))
        literal = directory / (arm + '.literal.bam')
        paths['calmd'][arm], paths['literal'][arm] = str(encoded), str(literal)
        with pysam.AlignmentFile(str(encoded)) as inp, pysam.AlignmentFile(str(literal), 'wb', header=inp.header) as out:
            inputs[arm] = []
            for r in inp:
                # Independent decoding oracle uses only the tiny fixture's
                # aligned query/ref pairs; production must not allocate N spans.
                seq = list(r.query_sequence)
                for q, p in r.get_aligned_pairs(matches_only=True):
                    if seq[q] == '=':
                        seq[q] = genome[r.reference_name][p]
                literal_r = copy.deepcopy(r)
                quals = literal_r.query_qualities
                literal_r.query_sequence = ''.join(seq)
                literal_r.query_qualities = quals
                assert literal_r.query_sequence == expected[r.query_name]
                inputs[arm].append((r, literal_r))
                out.write(literal_r)
    return genome, paths, expected, inputs


def assert_refused(fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except ValueError as exc:
        assert "Cannot decode consensus SEQ '='" in str(exc)
    else:
        raise AssertionError('Unresolved reference-relative SEQ must refuse')


def test_direct_evidence_and_chimeric_selection_preserve_inputs(tmp_path):
    genome, _, expected, inputs = make_family(tmp_path)
    for arm, pairs in inputs.items():
        for encoded, literal in pairs:
            original = encoded.to_string()
            decoded = decoded_alignment_copy(encoded, genome)
            assert (decoded is not encoded) == ('=' in encoded.query_sequence)
            assert decoded.to_string() == literal.to_string()
            assert decoded_alignment_copy(literal, genome) is literal
            left = extract_alignment_info(encoded, arm, genome)
            right = extract_alignment_info(literal, arm, genome)
            assert left == right
            assert score_alignment(left, genome) == score_alignment(right, genome)
            for fn in (_get_effective_5prime_clip, _get_effective_3prime_clip, _count_junction_proximity_errors):
                assert fn(encoded, genome) == fn(literal, genome)
            assert encoded.to_string() == original
    for index in range(8):
        encoded = {arm: pairs[index][0] for arm, pairs in inputs.items()}
        literal = {arm: pairs[index][1] for arm, pairs in inputs.items()}
        before = {arm: read.to_string() for arm, read in encoded.items()}
        assert select_best_chimeric(encoded, genome) == select_best_chimeric(literal, genome)
        assert before == {arm: read.to_string() for arm, read in encoded.items()}
        result = select_best_chimeric(encoded, genome)
        template = encoded[result.anchor_aligner]
        original = template.to_string()
        if '=' in template.query_sequence:
            try:
                build_chimeric_read(template, result.chimeric_ref_start, result.chimeric_cigar, result, template.header, anchor_read=template, aligner_reads=encoded)
            except ValueError as exc:
                assert 'decode the template' in str(exc)
            else:
                raise AssertionError('Standalone builder must require explicit template SEQ')
        decoded_template = decoded_alignment_copy(template, genome)
        built = build_chimeric_read(decoded_template, result.chimeric_ref_start, result.chimeric_cigar, result, template.header, anchor_read=decoded_template, aligner_reads=literal)
        assert built.query_sequence == expected[template.query_name]
        assert built.query_qualities == template.query_qualities
        assert template.to_string() == original


def test_full_standard_and_chimeric_consensus_encoding_parity(tmp_path):
    genome, paths, expected, _ = make_family(tmp_path)
    inputs = {p: Path(p).read_bytes() for arms in paths.values() for p in arms.values()}
    for chimeric in (False, True):
        outputs = {}
        for mode in ('literal', 'calmd'):
            path = tmp_path / f'{mode}.{int(chimeric)}.selected.bam'
            stats = run_consensus_selection(paths[mode], genome, str(path), n_workers=1, use_chimeric=chimeric)
            assert stats['total_reads'] == 8
            with pysam.AlignmentFile(str(path)) as bam:
                records = list(bam)
            assert len(records) == 8 and len({r.query_name for r in records}) == 8
            outputs[mode] = {r.query_name: r.to_string() for r in records}
            for r in records:
                assert r.query_sequence == expected[r.query_name]
                assert r.get_tag('Xa') == 'uLTRA'
        assert outputs['literal'] == outputs['calmd']
    assert all(Path(p).read_bytes() == content for p, content in inputs.items())


def test_unresolved_equals_refuse_and_literal_no_reference_is_unchanged():
    header = pysam.AlignmentHeader.from_references(['eq'], [200])
    r = pysam.AlignedSegment(header)
    r.query_name, r.reference_id, r.reference_start = 'bad', 0, 10
    r.cigarstring, r.query_sequence = '4M', '===='
    r.query_qualities = [31] * 4
    initial = r.to_string()
    for genome in (None, {}, {'other': 'C' * 200}, {'eq': 'C' * 12}):
        assert_refused(decoded_alignment_copy, r, genome)
        assert r.to_string() == initial
    negative = copy.deepcopy(r); negative.reference_start = -1
    assert_refused(decoded_alignment_copy, negative, {'eq': 'C' * 200})
    for cigar, seq in (('1S3M', '=CCC'), ('1I3M', '=CCC'), ('3M1S', 'CCC=')):
        unaligned = copy.deepcopy(r); unaligned.cigarstring = cigar; unaligned.query_sequence = seq
        for fn in (decoded_alignment_copy, extract_alignment_info, _get_effective_5prime_clip, _get_effective_3prime_clip, _count_junction_proximity_errors):
            if fn is extract_alignment_info:
                assert_refused(fn, unaligned, 'minimap2', {'eq': 'C' * 200})
            else:
                assert_refused(fn, unaligned, {'eq': 'C' * 200})
    literal = copy.deepcopy(r); literal.query_sequence = 'CCCC'
    assert decoded_alignment_copy(literal, None) is literal
    literal.query_sequence = None
    assert decoded_alignment_copy(literal, None) is literal


def test_late_decode_failure_removes_partial_output(tmp_path):
    header = pysam.AlignmentHeader.from_references(['eq'], [200])
    source = tmp_path / 'input.bam'
    with pysam.AlignmentFile(str(source), 'wb', header=header) as out:
        for i, seq in enumerate(('CCCC', '====')):
            r = pysam.AlignedSegment(header)
            r.query_name, r.reference_id, r.reference_start = f'read{i}', 0, 10
            r.cigarstring, r.query_sequence = '4M', seq
            out.write(r)
    before = source.read_bytes()
    for chimeric in (False, True):
        output = tmp_path / f'failed{int(chimeric)}.bam'
        assert_refused(run_consensus_selection, {'minimap2': str(source)}, {}, str(output), n_workers=1, batch_size=1, use_chimeric=chimeric)
        assert not output.exists()
        assert source.read_bytes() == before


def test_large_intron_decode_uses_no_aligned_pair_allocation():
    class NoPairs(pysam.AlignedSegment):
        def get_aligned_pairs(self, *args, **kwargs):
            raise AssertionError('Production decoder must walk CIGAR directly')
    h = pysam.AlignmentHeader.from_references(['eq'], [100004])
    r = NoPairs(h)
    r.query_name, r.reference_id, r.reference_start = 'longN', 0, 0
    r.cigarstring, r.query_sequence = '2M100000N2M', '===='
    result = decoded_alignment_copy(r, {'eq': 'AC' + 'N' * 100000 + 'GT'})
    assert result.query_sequence == 'ACGT' and r.query_sequence == '===='
