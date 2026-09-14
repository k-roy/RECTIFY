"""Canonical credit must survive extraction and whole-read consensus selection."""

import pysam
import pytest

from rectify.core.consensus.chimeric_consensus import (
    cigar_to_events, score_segment, select_best_chimeric,
)
from rectify.core.consensus.extract import (
    check_canonical_splice_sites, extract_alignment_info,
)
from rectify.core.consensus.select import select_best_alignment


# Literal genomic pairs: these expectations do not use the production RC helper.
PAIRS = [
    ('+', 'GT', 'AG'), ('+', 'GC', 'AG'), ('+', 'AT', 'AC'),
    ('-', 'CT', 'AC'), ('-', 'CT', 'GC'), ('-', 'GT', 'AT'),
]


def _fixture(left, right, strand, wrong_first=True):
    seq = list('C' * 400)
    seq[100:102] = left
    seq[198:200] = right
    genome = {'chrT': ''.join(seq)}
    header = pysam.AlignmentHeader.from_dict({'SQ': [{'SN': 'chrT', 'LN': 400}]})
    reads = {}
    pairs = [('wrong', '44M100N56M'), ('canonical', '40M100N60M')]
    if not wrong_first:
        pairs.reverse()
    for name, cigar in pairs:
        read = pysam.AlignedSegment(header)
        read.query_name = 'canonical-credit'
        read.reference_id = 0
        read.reference_start = 60
        read.is_reverse = strand == '-'
        read.cigarstring = cigar
        read.query_sequence = 'C' * 100
        read.query_qualities = [30] * 100
        reads[name] = read
    return genome, reads


@pytest.mark.parametrize('strand,left,right', PAIRS)
def test_extract_counts_paired_motifs_on_the_read_strand(strand, left, right):
    genome, reads = _fixture(left, right, strand)
    info = extract_alignment_info(reads['canonical'], 'canonical', genome)
    assert (info.canonical_count, info.non_canonical_count) == (1, 0)
    assert check_canonical_splice_sites([(100, 200)], 'chrT', genome,
                                        '-' if strand == '+' else '+') == (0, 1)


@pytest.mark.parametrize('strand,left,right', [
    ('+', 'AT', 'AG'), ('+', 'GT', 'AC'), ('+', 'GC', 'AC'),
    ('-', 'CT', 'AT'), ('-', 'GT', 'AC'), ('-', 'GT', 'GC'),
])
def test_mixed_pairs_never_get_canonical_credit(strand, left, right):
    genome, reads = _fixture(left, right, strand)
    assert check_canonical_splice_sites([(100, 200)], 'chrT', genome, strand) == (0, 1)
    events = cigar_to_events(reads['canonical'].cigartuples, 60)
    score = score_segment(events, 'interior', 'chrT', genome, strand=strand)
    assert score.n_canonical_junctions == 0
    assert score.score == -3


@pytest.mark.parametrize('strand,left,right', PAIRS)
@pytest.mark.parametrize('shift', [0, 1])
def test_segment_credit_includes_sequence_equivalent_placements(strand, left, right, shift):
    genome, _ = _fixture(left, right, strand)
    seq = list(genome['chrT'])
    if shift:
        # The canonical junction can slide right by one without changing SEQ.
        seq[200] = seq[100]
    genome['chrT'] = ''.join(seq)
    events = cigar_to_events([(0, 40 + shift), (3, 100), (0, 60 - shift)], 60)
    score = score_segment(events, 'interior', 'chrT', genome, strand=strand)
    assert score.n_canonical_junctions == 1
    assert score.n_novel_canonical_junctions == 1
    assert score.score == 5


@pytest.mark.parametrize('strand,left,right', PAIRS)
@pytest.mark.parametrize('wrong_first', [True, False])
def test_whole_selection_keeps_the_canonical_junction(strand, left, right, wrong_first):
    genome, reads = _fixture(left, right, strand, wrong_first)
    infos = {name: extract_alignment_info(read, name, genome) for name, read in reads.items()}
    result = select_best_alignment(infos, genome)
    assert result.best_aligner == 'canonical'
    assert result.best_alignment.cigar_string == '40M100N60M'

    stitched = select_best_chimeric(reads, genome)
    assert stitched.chimeric_ref_start == 60
    assert stitched.chimeric_cigar == [(0, 40), (3, 100), (0, 60)]
    assert stitched.anchor_is_reverse == (strand == '-')
    # Exercise real segment scoring as well as final geometry; no mocked winners.
    if not stitched.is_fallback:
        scores = [s.scores['canonical'] for s in stitched.all_segment_scores
                  if 'canonical' in s.scores]
        assert sum(s.n_canonical_junctions for s in scores) == 1


@pytest.mark.parametrize('start,end', [(-1, 20), (10, 401), (10, 10), (20, 10), (10, 12)])
def test_invalid_or_overlapping_boundaries_cannot_be_canonical(start, end):
    genome = {'chrT': 'GT' * 200}
    canonical, _ = check_canonical_splice_sites([(start, end)], 'chrT', genome)
    assert canonical == 0
