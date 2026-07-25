"""Junction-pool admission gate (rectify.core.consensus.pool_gate)."""

import pytest

from rectify.core.consensus.extract import AlignmentInfo
from rectify.core.consensus.pool_gate import (
    admit_junction,
    gated_alignment_junctions,
    junction_min_anchors,
)


def _aln(cigar, junctions, chrom='I', strand='+'):
    return AlignmentInfo(
        read_id='r', aligner='minimap2', chrom=chrom, strand=strand,
        reference_start=100, reference_end=500, cigar_string=cigar,
        junctions=list(junctions),
    )


class TestJunctionMinAnchors:
    def test_single_junction_takes_the_shorter_flank(self):
        assert junction_min_anchors('30M100N5M') == [5]
        assert junction_min_anchors('5M100N30M') == [5]

    def test_one_value_per_n_op_in_cigar_order(self):
        # 40M | N | 12M | N | 7M  ->  anchors min(40,12)=12 and min(12,7)=7
        assert junction_min_anchors('40M100N12M200N7M') == [12, 7]

    def test_soft_clip_does_not_count_as_anchor(self):
        # The 20S is clipped, not aligned: the anchor is the 3M only.
        assert junction_min_anchors('20S3M100N50M') == [3]

    def test_indel_terminates_the_contiguous_run(self):
        # 30M2I4M | N | 50M -> left run stops at the I, so the anchor is 4.
        assert junction_min_anchors('30M2I4M100N50M') == [4]

    def test_eq_and_x_ops_count_as_aligned(self):
        assert junction_min_anchors('10=5X100N40M') == [15]

    def test_adjacent_introns_give_zero_anchor(self):
        assert junction_min_anchors('30M100N200N40M') == [0, 0]

    def test_no_junction(self):
        assert junction_min_anchors('100M') == []

    @pytest.mark.parametrize('cigar', ['', None, 'garbage', '30M100'])
    def test_malformed_cigar_returns_empty(self, cigar):
        assert junction_min_anchors(cigar) == []


class TestAdmitJunction:
    def test_disabled_by_default(self):
        assert admit_junction(1, 500_000) is True

    def test_short_anchor_rejected(self):
        assert admit_junction(3, 200, min_anchor_bp=8) is False
        assert admit_junction(8, 200, min_anchor_bp=8) is True

    def test_mega_intron_rejected(self):
        assert admit_junction(30, 150_000, max_intron_len=3000) is False
        assert admit_junction(30, 2_999, max_intron_len=3000) is True

    def test_unknown_anchor_is_admitted_not_dropped(self):
        # A malformed CIGAR must not silently discard real junctions.
        assert admit_junction(None, 200, min_anchor_bp=12) is True


class TestGatedAlignmentJunctions:
    def test_off_by_default_admits_everything(self):
        aln = _aln('3M100N3M', [(200, 300)])
        assert gated_alignment_junctions(aln) == [('I', 200, 300, '+')]

    def test_short_anchor_junction_is_dropped(self):
        aln = _aln('3M100N3M', [(200, 300)])
        assert gated_alignment_junctions(aln, min_anchor_bp=8) == []

    def test_long_anchor_junction_is_kept(self):
        aln = _aln('30M100N40M', [(200, 300)])
        assert gated_alignment_junctions(aln, min_anchor_bp=8) == [('I', 200, 300, '+')]

    def test_annotated_junction_bypasses_the_gate(self):
        aln = _aln('3M100N3M', [(200, 300)])
        annotated = {('I', 200, 300, '+')}
        assert gated_alignment_junctions(
            aln, annotated_junctions=annotated, min_anchor_bp=12,
            max_intron_len=50) == [('I', 200, 300, '+')]

    def test_mixed_alignment_keeps_only_the_well_anchored_junction(self):
        # 40M | N(100) | 12M | N(200) | 3M  -> anchors 12 and 3
        aln = _aln('40M100N12M200N3M', [(200, 300), (312, 512)])
        assert gated_alignment_junctions(aln, min_anchor_bp=8) == [('I', 200, 300, '+')]

    def test_mega_intron_dropped_by_length(self):
        aln = _aln('40M150000N40M', [(200, 150_200)])
        assert gated_alignment_junctions(aln, max_intron_len=3000) == []

    def test_anchor_count_mismatch_falls_back_to_admitting(self):
        # CIGAR reports one N op but two junctions were recorded: rather than mis-pair
        # anchors with junctions, treat every anchor as unknown.
        aln = _aln('3M100N3M', [(200, 300), (400, 500)])
        assert len(gated_alignment_junctions(aln, min_anchor_bp=8)) == 2

    def test_minus_strand_key_uses_alignment_strand(self):
        aln = _aln('30M100N40M', [(200, 300)], strand='-')
        assert gated_alignment_junctions(aln, min_anchor_bp=8) == [('I', 200, 300, '-')]


class TestSelectIntegration:
    def test_select_best_alignment_accepts_the_new_kwargs(self):
        from rectify.core.consensus.select import select_best_alignment
        genome = {'I': 'A' * 1000}
        aln = _aln('30M100N40M', [(200, 300)])
        res = select_best_alignment({'minimap2': aln}, genome,
                                    pool_min_anchor_bp=8, pool_max_intron_len=3000)
        assert res.best_aligner == 'minimap2'
