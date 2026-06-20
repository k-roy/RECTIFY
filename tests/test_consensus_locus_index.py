"""Regression tests for the locus-scoped annotated-junction index in
``select.select_best_alignment``.

The 5' soft-clip rescue pool used to be built by scanning every annotated
junction on the read's whole chromosome (``for j in annotated_junctions: if
j[0] == chrom``). At genome-wide short-read scale (~50k reads × tens of
thousands of junctions/chrom) that is O(reads × junctions) and hung consensus
for >2h. It is now a memoized per-chrom position index queried by a locus
window that is a strict superset of the rescue's ``search_window_bp=300`` reach,
so scoring is byte-identical. These tests lock both the index correctness and
that equivalence.
"""
from typing import List, Tuple

from rectify.core.consensus.extract import AlignmentInfo
from rectify.core.consensus.select import (
    select_best_alignment,
    _query_annot_window,
    _annot_index,
)


def _ai(
    aligner: str,
    chrom: str = 'chrI',
    strand: str = '+',
    start: int = 1000,
    end: int = 1100,
    effective_five_prime_clip: int = 0,
    five_prime_softclip_seq: str = '',
    junctions: List[Tuple[int, int]] = None,
) -> AlignmentInfo:
    return AlignmentInfo(
        read_id='r1',
        aligner=aligner,
        chrom=chrom,
        strand=strand,
        reference_start=start,
        reference_end=end,
        cigar_string='100M',
        junctions=junctions or [],
        five_prime_softclip=effective_five_prime_clip,
        three_prime_softclip=0,
        effective_five_prime_clip=effective_five_prime_clip,
        effective_three_prime_clip=0,
        five_prime_softclip_seq=five_prime_softclip_seq,
        three_prime_atract_depth=0,
        corrected_3prime=end - 1,
        canonical_count=0,
    )


class TestQueryAnnotWindow:

    def test_start_in_window_included(self):
        annot = {('chrI', 1000, 1200, '+')}
        assert _query_annot_window(annot, 'chrI', 900, 1100) == annot

    def test_long_intron_captured_by_either_endpoint(self):
        # The + strand rescue filters on intron_end, the - strand on intron_start,
        # so a long intron with EITHER endpoint in the window must be captured.
        start_inside = {('chrI', 1_050, 200_000, '+')}   # start in window, end far
        assert _query_annot_window(start_inside, 'chrI', 1000, 1100) == start_inside
        end_inside = {('chrI', 50, 1_050, '+')}          # start far, end in window
        assert _query_annot_window(end_inside, 'chrI', 1000, 1100) == end_inside

    def test_outside_window_excluded(self):
        annot = {('chrI', 500_000, 500_200, '+'), ('chrI', 1000, 1100, '+')}
        got = _query_annot_window(annot, 'chrI', 900, 1200)
        assert got == {('chrI', 1000, 1100, '+')}

    def test_other_chrom_excluded(self):
        annot = {('chrII', 1000, 1100, '+')}
        assert _query_annot_window(annot, 'chrI', 900, 1200) == set()

    def test_empty_for_unknown_chrom(self):
        assert _query_annot_window(set(), 'chrI', 0, 100) == set()

    def test_memoized_index_is_reused(self):
        annot = frozenset({('chrI', 1000, 1100, '+')})
        idx1 = _annot_index(annot)
        idx2 = _annot_index(annot)
        assert idx1 is idx2  # same object -> cached, not rebuilt


class TestLocusWindowEquivalence:
    """select_best_alignment must be invariant to annotated junctions outside
    the read's locus window — the whole point of the index optimisation."""

    def _genome(self, size=1_000_000):
        return {'chrI': 'A' * size}

    def test_far_junctions_do_not_change_scoring(self):
        # Clipped alignment that 5'-rescues against a NEAR junction
        # (intron 100-200, + strand): clip 'AAAAA' matches genome[95:100].
        near = {('chrI', 100, 200, '+')}
        far = {('chrI', 500_000 + i * 100, 500_000 + i * 100 + 50, '+')
               for i in range(5000)}

        def fresh():
            return {
                'a_clipped': _ai('a_clipped', start=200, end=250,
                                 effective_five_prime_clip=5,
                                 five_prime_softclip_seq='AAAAA'),
                'b_clean': _ai('b_clean', start=200, end=250),
            }

        genome = self._genome()
        al_near = fresh()
        res_near = select_best_alignment(al_near, genome, annotated_junctions=near)
        al_both = fresh()
        res_both = select_best_alignment(al_both, genome,
                                         annotated_junctions=near | far)

        # The 5000 far junctions must not perturb the outcome or any score.
        assert res_near.best_aligner == res_both.best_aligner
        assert (al_near['a_clipped'].junction_score
                == al_both['a_clipped'].junction_score)
        # And the near junction must still rescue the clip (score 0, not -10).
        assert al_both['a_clipped'].junction_score == 0.0
