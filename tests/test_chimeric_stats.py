#!/usr/bin/env python3
"""Tests for the ChimericStats instrumentation.

Covers the measurement layer wired in to separate genuine multi-aligner
segment wins from degenerate fallback pass-throughs, and to expose WHY each
contested segment was decided (the annotation-off confound check).  See the
`_classify_term` priority order and the `decisive_term` aggregate.

Author: Kevin R. Roy
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rectify.core.consensus.chimeric_consensus import (
    ChimericResult,
    ChimericSegment,
    ChimericStats,
    SegmentScore,
    _classify_term,
)


def _result(**kw):
    base = dict(
        read_id='r', is_chimeric=False, segment_winners=[],
        chimeric_cigar=[], chimeric_ref_start=0, confidence='high',
    )
    base.update(kw)
    return ChimericResult(**base)


# ── fallback / single-aligner / chimeric classification ─────────────────────

def test_fallback_counted_separately_from_single_aligner():
    st = ChimericStats()
    st.update(_result(is_fallback=True, n_segments=1,
                      segment_winners=[('whole', 'GMAP', 0, 100)]))
    st.update(_result(is_fallback=False, is_chimeric=False, n_segments=1,
                      segment_winners=[('whole', 'minimap2', 0, 100)]))
    st.update(_result(is_fallback=False, is_chimeric=True, n_segments=2,
                      segment_winners=[('five_prime', 'GMAP', 0, 50),
                                       ('three_prime', 'minimap2', 50, 100)]))
    assert st.fallback_reads == 1
    assert st.single_aligner_reads == 1
    assert st.chimeric_reads == 1
    assert st.total_reads == 3


# ── _classify_term priority order ───────────────────────────────────────────

def test_classify_annotated_beats_canonical():
    w = SegmentScore(aligner='mm2', score=10, n_annotated_junctions=1,
                     n_canonical_junctions=1)
    l = SegmentScore(aligner='GMAP', score=2, n_annotated_junctions=0,
                     n_canonical_junctions=1)
    assert _classify_term(w, l) == 'annotated'


def test_classify_canonical_when_annotation_off():
    """Annotation dead (both 0 annotated) → decision falls to canonical."""
    w = SegmentScore(aligner='mm2', score=5, n_annotated_junctions=0,
                     n_canonical_junctions=1)
    l = SegmentScore(aligner='GMAP', score=-3, n_annotated_junctions=0,
                     n_canonical_junctions=0)
    assert _classify_term(w, l) == 'canonical'


def test_classify_false_3prime():
    w = SegmentScore(aligner='mm2', score=0, has_false_3prime_junction=False)
    l = SegmentScore(aligner='GMAP', score=-15, has_false_3prime_junction=True)
    assert _classify_term(w, l) == 'false_3prime'


def test_classify_softclip_overhang():
    """The 'went for it' aggressive-overhang term."""
    w = SegmentScore(aligner='mm2', score=0, n_softclip=2)
    l = SegmentScore(aligner='GMAP', score=-20, n_softclip=20)
    assert _classify_term(w, l) == 'softclip'


def test_classify_residue_fallback():
    w = SegmentScore(aligner='mm2', score=5, n_matches=50)
    l = SegmentScore(aligner='GMAP', score=4, n_matches=49, n_mismatches=1)
    assert _classify_term(w, l) == 'residue'


# ── decisive_term aggregation on a contested segment ────────────────────────

def test_decisive_term_aggregated_for_contested_segment():
    seg = ChimericSegment(q_start=0, q_end=50, position='interior',
                          winning_aligner='minimap2')
    seg.scores = {
        'minimap2': SegmentScore(aligner='minimap2', score=5,
                                 n_canonical_junctions=1),
        'GMAP': SegmentScore(aligner='GMAP', score=-3,
                             n_canonical_junctions=0),
    }
    res = _result(is_chimeric=True, n_segments=1,
                  segment_winners=[('interior', 'minimap2', 0, 50)],
                  all_segment_scores=[seg])
    st = ChimericStats()
    st.update(res)
    assert st.decisive_term['canonical'] == 1
    assert st.segments_won_by_aligner['minimap2'] == 1
    assert st.segments_competed_by_aligner['GMAP'] == 1
    assert st.loss_reason_by_aligner['GMAP']['canonical'] == 1
    assert st.win_reason_by_aligner['minimap2']['canonical'] == 1


def test_uncontested_segment_not_in_loss_reasons():
    """Identical scores → not contested → no loss-reason attribution."""
    seg = ChimericSegment(q_start=0, q_end=50, position='interior',
                          winning_aligner='minimap2')
    seg.scores = {
        'minimap2': SegmentScore(aligner='minimap2', score=5),
        'GMAP': SegmentScore(aligner='GMAP', score=5),
    }
    res = _result(is_chimeric=False, n_segments=1,
                  segment_winners=[('interior', 'minimap2', 0, 50)],
                  all_segment_scores=[seg])
    st = ChimericStats()
    st.update(res)
    assert sum(st.decisive_term.values()) == 0
    # still credited as competed/won (just not a contested decision)
    assert st.segments_won_by_aligner['minimap2'] == 1


def test_summary_dict_is_json_serializable():
    import json
    seg = ChimericSegment(q_start=0, q_end=50, position='interior',
                          winning_aligner='GMAP')
    seg.scores = {
        'GMAP': SegmentScore(aligner='GMAP', score=8, n_annotated_junctions=1),
        'minimap2': SegmentScore(aligner='minimap2', score=2),
    }
    st = ChimericStats()
    st.update(_result(is_chimeric=True, all_segment_scores=[seg],
                      segment_winners=[('interior', 'GMAP', 0, 50)]))
    st.update(_result(is_fallback=True))
    blob = json.dumps(st.summary_dict())  # must not raise
    d = json.loads(blob)
    assert d['fallback_reads'] == 1
    assert d['decisive_term']['annotated'] == 1
    assert d['win_reason_by_aligner']['GMAP']['annotated'] == 1


if __name__ == '__main__':
    import pytest
    sys.exit(pytest.main([__file__, '-q']))
