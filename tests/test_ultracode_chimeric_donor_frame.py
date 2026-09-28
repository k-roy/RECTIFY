"""Native donor search preserves the original retained molecule interval."""
import hashlib
import random
from pathlib import Path

import pysam
import pytest

from rectify.core.consensus.chimeric_consensus import select_best_chimeric
from rectify.core.consensus.consensus import run_consensus_selection
from rectify.utils.genome import register_genome_contigs
from tests.test_ultracode_chimeric_hardclips import full_molecule_map


def donor_family(reverse=False, opposite=False, wrong_first=True):
    rng = random.Random(680921)
    reference = ''.join(rng.choice('ACGT') for _ in range(1500))
    rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    genome = {'donorP': reference, 'donorM': rc(reference)}
    register_genome_contigs(genome)
    header = pysam.AlignmentHeader.from_references(list(genome), [1500, 1500])
    start, stored = 177, 103
    query = reference[start:start + stored]
    full_molecule = 'CGTACGC' + query + 'AGCTGACGTCC'
    qualities = [25 + i % 12 for i in range(stored)]
    def read(name, seq, left, right):
        r = pysam.AlignedSegment(header)
        r.query_name = 'physical_donor_interval'
        r.flag = 16 if reverse else 0
        r.reference_id = 0
        r.reference_start = start if name == 'anchor' else start + 12
        r.mapping_quality = 60
        r.cigartuples = [(5, left)] + ([(0, stored)] if name == 'anchor' else
                                      [(4, 12), (0, stored - 12)]) + [(5, right)]
        r.query_sequence = seq
        if seq is not None:
            r.query_qualities = qualities
        if name == 'valid' and opposite:
            r.flag ^= 16
            r.reference_id = 1
            r.reference_start = 1500 - (start + stored)
            r.cigartuples = list(reversed(r.cigartuples))
            r.query_sequence = rc(seq)
            r.query_qualities = list(reversed(qualities))
        return r
    anchor = read('anchor', None, 7, 11)
    valid = read('valid', query, 7, 11)
    wrong = read('wrong', full_molecule[11:114], 11, 7)
    names = [('wrong', wrong), ('valid', valid)] if wrong_first else [
        ('valid', valid), ('wrong', wrong)]
    arms = dict([('anchor', anchor)] + names)
    expected = read('anchor', query, 7, 11)
    return arms, expected, genome


def write_and_run(arms, genome, tmp_path):
    paths = {}
    for name, read in arms.items():
        path = tmp_path / (name + '.bam')
        with pysam.AlignmentFile(path, 'wb', header=read.header) as bam:
            bam.write(read)
        paths[name] = str(path)
    before = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths.values()}
    output = tmp_path / 'selected.bam'
    try:
        run_consensus_selection(paths, genome, str(output), n_workers=1, use_chimeric=True)
    finally:
        assert before == {p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                          for p in paths.values()}
    with pysam.AlignmentFile(output) as bam:
        return list(bam)


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('opposite', [False, True])
@pytest.mark.parametrize('wrong_first', [False, True])
def test_native_search_reaches_compatible_donor(reverse, opposite, wrong_first, tmp_path):
    arms, expected, genome = donor_family(reverse, opposite, wrong_first)
    saved = {name: read.to_string() for name, read in arms.items()}
    assert select_best_chimeric(arms, genome).anchor_aligner == 'anchor'
    outputs = write_and_run(arms, genome, tmp_path)
    assert len(outputs) == 1
    out = outputs[0]
    assert out.reference_name == expected.reference_name
    assert out.reference_start == expected.reference_start
    assert out.flag == expected.flag
    assert out.cigartuples == expected.cigartuples
    assert out.query_sequence == expected.query_sequence
    assert list(out.query_qualities) == list(expected.query_qualities)
    assert full_molecule_map(out) == full_molecule_map(expected)
    assert saved == {name: read.to_string() for name, read in arms.items()}


@pytest.mark.parametrize('reverse', [False, True])
def test_native_no_compatible_donor_refuses_without_final_bam(reverse, tmp_path):
    arms, _, genome = donor_family(reverse)
    arms.pop('valid')
    saved = {name: read.to_string() for name, read in arms.items()}
    assert select_best_chimeric(arms, genome).anchor_aligner == 'anchor'
    with pytest.raises(ValueError, match='compatible original query frame'):
        write_and_run(arms, genome, tmp_path)
    assert not (tmp_path / 'selected.bam').exists()
    assert saved == {name: read.to_string() for name, read in arms.items()}


@pytest.mark.parametrize('reverse', [False, True])
def test_native_all_sequence_missing_keeps_existing_empty_output(reverse, tmp_path):
    arms, _, genome = donor_family(reverse)
    # One real placement without any available query donor takes the existing
    # documented warning/skip path; this test does not claim read retention.
    arms = {'anchor': arms['anchor']}
    assert write_and_run(arms, genome, tmp_path) == []
