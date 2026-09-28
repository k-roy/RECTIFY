"""Preserve original molecule coordinates during an active RNA3 clip."""
from array import array

import pysam
import pytest

from rectify.core.bam.read_edits import (
    clip_read_to_corrected_3prime,
    softclip_read_to_corrected_3prime,
)


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('mode', ['hard', 'soft'])
@pytest.mark.parametrize('old_h', [0, 3, 9])
def test_existing_rna3_hardclip_keeps_original_molecule_frame(reverse, mode, old_h):
    header = pysam.AlignmentHeader.from_references(['clipT'], [1000])
    read = pysam.AlignedSegment(header)
    read.query_name = f'tail_{reverse}_{mode}_{old_h}'
    read.reference_id = 0
    read.reference_start = 100
    read.flag = 16 if reverse else 0
    # Include an interior junction, terminal insertion and soft clip. The
    # explicit clip target remains within the final M and does not cross N.
    ops = [(5, 2), (0, 10), (3, 50), (0, 20), (1, 2), (4, 6)]
    if old_h:
        ops.append((5, old_h))
    if reverse:
        ops = ops[::-1]
    read.cigartuples = ops
    seq = 'ACGT' * 8 + 'AAAAAA'
    read.query_sequence = seq[::-1] if reverse else seq
    read.query_qualities = array('B', range(20, 58))
    before = read.__copy__()
    fn = clip_read_to_corrected_3prime if mode == 'hard' else softclip_read_to_corrected_3prime
    target = 104 if reverse else 175
    assert fn(read, target, '-' if reverse else '+')
    assert read.infer_read_length() == before.infer_read_length()
    expected = [(5, 2), (0, 10), (3, 50), (0, 16)]
    if mode == 'hard':
        expected.append((5, old_h + 12))
        wanted = before.query_sequence[12:] if reverse else before.query_sequence[:-12]
        quals = before.query_qualities[12:] if reverse else before.query_qualities[:-12]
        assert read.query_sequence == wanted
        assert read.query_qualities == quals
    else:
        expected.append((4, 12))
        if old_h:
            expected.append((5, old_h))
        assert read.query_sequence == before.query_sequence
        assert read.query_qualities == before.query_qualities
    assert read.cigartuples == (expected[::-1] if reverse else expected)
    assert (read.reference_start if reverse else read.reference_end - 1) == target
    original_map = dict(before.get_aligned_pairs(matches_only=True))
    offset = 12 if reverse and mode == 'hard' else 0
    for q, r in read.get_aligned_pairs(matches_only=True):
        assert original_map[q + offset] == r
    # Applying the already satisfied target never touches any part of SAM.
    unchanged = read.to_string()
    assert not fn(read, target, '-' if reverse else '+')
    assert read.to_string() == unchanged
