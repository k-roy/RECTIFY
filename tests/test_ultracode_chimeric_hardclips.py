"""Original molecule coordinates survive real cross-arm CIGAR assembly."""
import copy

import pysam
import pytest

from rectify.core.consensus.chimeric_consensus import (
    _cigar_query_frame, build_chimeric_read, build_query_ref_map,
    select_best_chimeric,
)
from tests.test_ultracode_chimeric_boundary_ownership import make_family


def with_hardclips(read, left, right):
    read = copy.copy(read)
    read.cigartuples = ([(5, left)] if left else []) + list(read.cigartuples) + (
        [(5, right)] if right else [])
    return read


def full_molecule_map(read):
    left, stored, right = _cigar_query_frame(read.cigartuples)
    extent = left + stored + right
    return {(extent - 1 - (left + q) if read.is_reverse else left + q): p
            for q, p in build_query_ref_map(read).items()}


def build(arms, genome, annotation, template=None):
    result = select_best_chimeric(arms, genome, annotation)
    anchor = arms[result.anchor_aligner]
    out = build_chimeric_read(
        template if template is not None else anchor,
        result.chimeric_ref_start, result.chimeric_cigar, result,
        anchor.header, anchor_read=anchor, aligner_reads=arms,
    )
    return result, out


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('soft', [False, True])
@pytest.mark.parametrize('reorder', [False, True])
@pytest.mark.parametrize('ends', [(0, 0), (5, 0), (0, 9), (5, 9)])
def test_cross_arm_original_query_frame(reverse, soft, reorder, ends):
    arms, truth, genome, annotation = make_family(reverse, shared='I', soft=soft)
    arms = {name: with_hardclips(read, *ends) for name, read in arms.items()}
    truth = with_hardclips(truth, *ends)
    if reorder:
        arms = dict(reversed(list(arms.items())))
    before = {name: read.to_string() for name, read in arms.items()}
    result, out = build(arms, genome, annotation)
    assert result.is_chimeric and not result.is_fallback
    assert out.cigartuples == truth.cigartuples
    assert out.query_sequence == truth.query_sequence
    assert list(out.query_qualities) == list(truth.query_qualities)
    assert _cigar_query_frame(out.cigartuples) == _cigar_query_frame(truth.cigartuples)
    assert full_molecule_map(out) == full_molecule_map(truth)
    assert before == {name: read.to_string() for name, read in arms.items()}


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('reorder', [False, True])
def test_different_original_intervals_fall_back(reverse, reorder):
    arms, _, genome, annotation = make_family(reverse, soft=True)
    arms = {name: with_hardclips(read, 5 if name == 'minimap2' else 9,
                                 9 if name == 'minimap2' else 5)
            for name, read in arms.items()}
    if reorder:
        arms = dict(reversed(list(arms.items())))
    result, out = build(arms, genome, annotation)
    anchor = arms[result.anchor_aligner]
    assert result.is_fallback and not result.is_chimeric
    assert out.cigartuples == anchor.cigartuples
    assert out.query_sequence == anchor.query_sequence
    assert full_molecule_map(out) == full_molecule_map(anchor)


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('opposite_donor', [False, True])
@pytest.mark.parametrize('missing_anchor_sequence', [False, True])
def test_builder_donor_frame_and_strand(reverse, opposite_donor, missing_anchor_sequence):
    arms, _, genome, annotation = make_family(reverse, soft=True)
    anchor = with_hardclips(arms['minimap2'], 5, 9)
    result = select_best_chimeric({'minimap2': anchor}, genome, annotation)
    donor = copy.copy(anchor)
    if opposite_donor:
        seq, qual = donor.query_sequence, donor.query_qualities
        donor.flag ^= 16
        donor.cigartuples = list(reversed(donor.cigartuples))
        donor.query_sequence = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        donor.query_qualities = qual[::-1]
    before = donor.to_string()
    expected = copy.copy(anchor)
    if missing_anchor_sequence:
        anchor.query_sequence = None
    out = build_chimeric_read(donor, result.chimeric_ref_start, result.chimeric_cigar,
                              result, anchor.header, anchor_read=anchor)
    assert out.cigartuples == expected.cigartuples
    assert out.query_sequence == expected.query_sequence
    assert list(out.query_qualities) == list(expected.query_qualities)
    assert full_molecule_map(out) == full_molecule_map(expected)
    assert donor.to_string() == before


@pytest.mark.parametrize('reverse', [False, True])
def test_builder_refuses_wrong_origin_even_equal_length_and_extent(reverse):
    arms, _, genome, annotation = make_family(reverse, soft=True)
    anchor = with_hardclips(arms['minimap2'], 5, 9)
    donor = with_hardclips(arms['minimap2'], 9, 5)
    result = select_best_chimeric({'minimap2': anchor}, genome, annotation)
    before = (anchor.to_string(), donor.to_string())
    with pytest.raises(ValueError, match='original query frames'):
        build_chimeric_read(donor, result.chimeric_ref_start, result.chimeric_cigar,
                            result, anchor.header, anchor_read=anchor)
    assert before == (anchor.to_string(), donor.to_string())


def test_internal_hardclip_is_not_a_valid_query_frame():
    with pytest.raises(ValueError, match='internal hard clip'):
        _cigar_query_frame([(0, 10), (5, 3), (0, 10)])
