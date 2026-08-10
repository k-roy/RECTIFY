"""653: quality-aware consensus construction for PCR-cDNA clusters.

PCR-cDNA clusters mix two sequencing-time regimes (UMI/SSP end at pore entry vs exit,
planning/650) with ~1.7x different error rates and a per-read positional gradient, so both
the representative pick and the pileup vote should listen to basecall quality.
"""
from __future__ import annotations

import pysam

from rectify.core.cdna.cluster import pick_representative
from rectify.core.cdna.consensus import pileup_consensus

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6"},
    "SQ": [{"SN": "chrT", "LN": 1000}],
})


def _seg(name, seq, qual_char=None, start=0, mapq=60):
    seg = pysam.AlignedSegment(HEADER)
    seg.query_name = name
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = start
    seg.mapping_quality = mapq
    seg.cigartuples = [(0, len(seq))]
    seg.query_sequence = seq
    if qual_char is not None:
        seg.query_qualities = pysam.qualitystring_to_array(qual_char * len(seq))
    return seg


class TestPickRepresentative:
    def test_prefers_higher_mean_quality_over_span(self):
        long_noisy = _seg("long_noisy", "A" * 100, qual_char="+")   # Q10
        short_clean = _seg("short_clean", "A" * 90, qual_char="?")  # Q30
        assert pick_representative([long_noisy, short_clean]) is short_clean

    def test_no_qualities_falls_back_to_span(self):
        long_r = _seg("long", "A" * 100)
        short_r = _seg("short", "A" * 90)
        assert pick_representative([long_r, short_r]) is long_r

    def test_deterministic_on_full_tie(self):
        a = _seg("a", "A" * 50, qual_char="?")
        b = _seg("b", "A" * 50, qual_char="?")
        assert pick_representative([b, a]) is a  # name tiebreak


class TestQualityWeightedPileup:
    def test_quality_outvotes_count_tie(self):
        # 2v2 strand-split disagreement at position 5: 'T' carried at Q40, 'G' at Q10.
        # A raw-count vote ties (old code broke the tie alphabetically toward 'G');
        # the quality-weighted vote must take 'T'.
        seq_t = "AAAAATAAAA"
        seq_g = "AAAAAGAAAA"
        cluster = [
            _seg("t1", seq_t, qual_char="I"),  # Q40
            _seg("t2", seq_t, qual_char="I"),
            _seg("g1", seq_g, qual_char="+"),  # Q10
            _seg("g2", seq_g, qual_char="+"),
        ]
        res = pileup_consensus(cluster)
        assert res is not None
        assert res[0] == seq_t

    def test_unweighted_when_no_qualities(self):
        # without quality strings every read votes 1: 2 'T' vs 1 'G' -> majority 'T'
        cluster = [
            _seg("t1", "AAAAATAAAA"),
            _seg("t2", "AAAAATAAAA"),
            _seg("g1", "AAAAAGAAAA"),
        ]
        res = pileup_consensus(cluster)
        assert res is not None
        assert res[0] == "AAAAATAAAA"
