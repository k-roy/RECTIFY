"""ISSUE-085: a minus-strand 3' clip that strips a dangling D/N must move the record by that op.

`clip_read_to_corrected_3prime` / `softclip_read_to_corrected_3prime` walk reference bases off the left end of a
minus-strand record and then strip any D/N/I left at the new start (a record cannot begin with D or N). The
stripped D/N consumed reference, so the first surviving base sits `length` bp further right; the code kept
`reference_start = current_start + n_ref_removed` without it, and every surviving base was drawn upstream of
where the aligner put it. Found on the chr5 cohort census (job 44521649): 5cac4ae6, `25S6M650I2046N45M...`
clipped to 40832462, was written 2,046 bp off on master, on the branch and on every census variant; the TSV
`junctions` column disagreed with the BAM N ops on 3 of 62,602 reads, all `polya_walkback`.

The invariant pinned here: every base that survives the clip keeps the reference position the aligner gave it.
"""
from array import array

import pysam
import pytest

from rectify.core.bam.read_edits import (
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
)


def _read(start, ops, seq, reverse=True):
    header = pysam.AlignmentHeader.from_references(['chr5'], [200_000_000])
    r = pysam.AlignedSegment(header)
    r.query_name = 'issue085'
    r.reference_id = 0
    r.reference_start = start
    r.flag = 16 if reverse else 0
    r.cigartuples = ops
    r.query_sequence = seq
    r.query_qualities = array('B', [30] * len(seq))
    return r


def _surviving_pairs_unchanged(before, after, removed_query):
    """Every aligned base still in the record maps to the reference position it had before the clip
    (`removed_query` query bases were cut from the left in hard mode; 0 in soft mode)."""
    original = dict(before.get_aligned_pairs(matches_only=True))
    pairs = after.get_aligned_pairs(matches_only=True)
    assert pairs, 'nothing survived'
    for q, r in pairs:
        assert original[q + removed_query] == r, (q, r, original.get(q + removed_query))


@pytest.mark.parametrize('mode', ['hard', 'soft'])
def test_5cac4ae6_geometry_dangling_n_after_walk(mode):
    """`25S 6M 650I 2046N 45M 3I 5M 2S` at 40832457, clipped to 40832463 (the whole 6M goes): the 650I and the
    2046N dangle after the walk and are stripped; the 45M must stay at 40834509. (The census record carried
    `681H45M…` at 40832463, i.e. the 6M gone and the start moved by 6 only.)"""
    ops = [(4, 25), (0, 6), (1, 650), (3, 2046), (0, 45), (1, 3), (0, 5), (4, 2)]
    seq = 'A' * 25 + 'C' * 6 + 'G' * 650 + 'T' * 45 + 'C' * 3 + 'A' * 5 + 'G' * 2
    read = _read(40832457, ops, seq)
    before = read.__copy__()
    fn = clip_read_to_corrected_3prime if mode == 'hard' else softclip_read_to_corrected_3prime
    assert fn(read, 40832463, '-')
    assert read.reference_start == 40832457 + 6 + 2046 == 40834509
    body = [op for op in read.cigartuples if op[0] not in (4, 5)]
    assert body[0] == (0, 45)
    assert body[0][0] not in (2, 3)
    removed = 25 + 6 + 650 if mode == 'hard' else 0
    if mode == 'hard':
        assert read.cigartuples[0] == (5, 681)
        assert read.query_sequence == seq[681:]
    else:
        assert read.cigartuples[0] == (4, 681)
        assert read.query_sequence == seq
    assert read.infer_read_length() == before.infer_read_length()
    _surviving_pairs_unchanged(before, read, removed)


def test_5cac4ae6_writer_path_clip_then_a_run():
    """The writer's real sequence: `clip_read_to_corrected_3prime` to 40832462 leaves `30H 1M 650I 2046N …`,
    then `_hardclip_trailing_a_run` removes the one remaining T (an A in RNA) through the same function, which
    strips the dangling 650I and 2046N. The census record was `681H45M…` at 40832463; it must be at 40834509."""
    from rectify.core.bam.read_edits import _hardclip_trailing_a_run
    ops = [(4, 25), (0, 6), (1, 650), (3, 2046), (0, 45), (1, 3), (0, 5), (4, 2)]
    seq = 'A' * 25 + 'C' * 5 + 'T' + 'G' * 650 + 'T' * 45 + 'C' * 3 + 'A' * 5 + 'G' * 2
    read = _read(40832457, ops, seq)
    before = read.__copy__()
    assert clip_read_to_corrected_3prime(read, 40832462, '-')
    assert read.cigartuples[:3] == [(5, 30), (0, 1), (1, 650)] and read.reference_start == 40832462
    assert _hardclip_trailing_a_run(read, '-')
    assert read.cigartuples[0] == (5, 681)
    assert read.cigartuples[1] == (0, 45)
    assert read.reference_start == 40834509
    assert read.infer_read_length() == before.infer_read_length()
    _surviving_pairs_unchanged(before, read, 681)


@pytest.mark.parametrize('mode', ['hard', 'soft'])
def test_deletion_span_minus_keeps_surviving_bases_in_place(mode):
    """`6=4D6=` at 1000, corrected_3prime 1006 (inside the deletion): the walk removes the left 6=, the 4D
    dangles and is stripped, and the surviving 6= must stay at 1010, not 1006."""
    read = _read(1000, [(7, 6), (2, 4), (7, 6)], 'A' * 6 + 'G' * 6)
    before = read.__copy__()
    fn = clip_read_to_corrected_3prime if mode == 'hard' else softclip_read_to_corrected_3prime
    assert fn(read, 1006, '-')
    assert read.reference_start == 1010
    _surviving_pairs_unchanged(before, read, 6 if mode == 'hard' else 0)


def test_plus_strand_unchanged_by_the_fix():
    """The plus-strand mirror strips trailing D/N without moving the start; it is not touched."""
    read = _read(1000, [(7, 6), (2, 4), (7, 6)], 'A' * 6 + 'G' * 6, reverse=False)
    before = read.__copy__()
    assert clip_read_to_corrected_3prime(read, 1006, '+')
    assert read.reference_start == 1000
    assert read.cigartuples == [(7, 6), (5, 6)]
    _surviving_pairs_unchanged(before, read, 0)
