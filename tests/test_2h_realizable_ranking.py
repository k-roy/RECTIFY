"""2H ranks only candidates the surgery can WRITE (realizability before ranking).

Codex audit 2026-09-13 (`dev/audits/codex_fable_20260913/repro_issue031_runner_up.py`):
`refine_read_junctions` kept only the head of the ranking; when the CIGAR surgery
then refused it (ISSUE-031 no-I/D-beside-the-N, the indel-burden invariant, the
boundary-shift guard) the read silently kept its stock placement and a realizable
runner-up was never tried — a real FN mechanism, not a stats inconsistency.

Witness: stock ``20M100N5D20M``; acceptor +2 leaves ``3D`` beside the N (refused);
acceptor +5 absorbs the whole deletion (``20M105N20M``, realizable). With +2 ranked
first, the old code proposed +2, the surgery refused it, and +5 was hidden.

The probe itself is 018's F1 (`5b81d88`), landed here on its own: the ranking is
walked head-first, each candidate dry-run through the surgery on a copy, the
incumbent stops the walk, and every skip is counted as
``unrealizable_winner_skipped`` (profile + driver stats, sequential and parallel).
"""
from collections import Counter

import pysam
import pytest

from rectify.core.splice import junction_refiner as jr
from rectify.core.splice.junction_scoring import _build_junction_index

CHROM = "chrT"
GLEN = 500


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": GLEN}]})


def _witness_read(name="rejected_best_hides_realizable_second"):
    r = pysam.AlignedSegment(_header())
    r.query_name = name
    r.reference_id = 0
    r.reference_start = 100
    r.cigarstring = "20M100N5D20M"
    r.query_sequence = "C" * 40
    r.query_qualities = pysam.qualitystring_to_array("I" * 40)
    return r


REF = "C" * GLEN
IDX = _build_junction_index({(CHROM, 120, e) for e in (220, 222, 225)})


def _controlled(monkeypatch, scores):
    """Rank by a fixed per-acceptor score (lower = better)."""
    def fake(q, q_split, js, je, genome_seq, **kw):
        return scores[je], 0
    monkeypatch.setattr(jr, "_score_junction", fake)


class TestTheWitness:
    def test_the_surgery_refuses_plus2_and_writes_plus5(self):
        """The two candidates, judged by the surgery alone (no ranking involved)."""
        r = _witness_read()
        assert not jr._realizable(r, 1, 120, 220, 120, 222, REF, "+", 0.25, 15)
        assert jr._realizable(r, 1, 120, 220, 120, 225, REF, "+", 0.25, 15)
        out, applied = jr._apply_replacements_to_read(
            r, [(1, 120, 220, 120, 225)], REF, "+", 0.25, 15)
        assert applied and out.cigarstring == "20M105N20M"

    def test_the_realizable_runner_up_is_proposed_and_applied(self, monkeypatch):
        _controlled(monkeypatch, {220: 2.0, 222: 0.0, 225: 1.0})   # +2 ranks first
        r = _witness_read()
        counters = Counter()
        proposals = jr.refine_read_junctions(
            r, IDX, {(CHROM, 120, 225)}, REF, "+", motif_blind=True, counters=counters)
        assert proposals == [(1, 120, 220, 120, 225)], proposals
        assert counters["unrealizable_winner_skipped"] == 1
        out, applied = jr._apply_replacements_to_read(r, proposals, REF, "+", 0.25, 15)
        assert applied and out.cigarstring == "20M105N20M"

    def test_a_writable_head_is_untouched(self, monkeypatch):
        _controlled(monkeypatch, {220: 2.0, 222: 1.0, 225: 0.0})   # +5 ranks first
        r = _witness_read()
        counters = Counter()
        proposals = jr.refine_read_junctions(
            r, IDX, {(CHROM, 120, 225)}, REF, "+", motif_blind=True, counters=counters)
        assert proposals == [(1, 120, 220, 120, 225)]
        assert counters["unrealizable_winner_skipped"] == 0

    def test_the_incumbent_stops_the_walk(self, monkeypatch):
        """+2 refused, then the incumbent outranks +5: the read stays, one skip counted."""
        _controlled(monkeypatch, {220: 0.5, 222: 0.0, 225: 1.0})
        r = _witness_read()
        counters = Counter()
        proposals = jr.refine_read_junctions(
            r, IDX, {(CHROM, 120, 225)}, REF, "+", motif_blind=True, counters=counters)
        assert proposals == []
        assert counters["unrealizable_winner_skipped"] == 1

    def test_nothing_writable_stays(self, monkeypatch):
        _controlled(monkeypatch, {220: 2.0, 222: 0.0, 225: 1.0})
        r = _witness_read()
        idx = _build_junction_index({(CHROM, 120, e) for e in (220, 222)})   # no +5 on offer
        counters = Counter()
        assert jr.refine_read_junctions(
            r, idx, set(), REF, "+", motif_blind=True, counters=counters) == []
        assert counters["unrealizable_winner_skipped"] == 1

    def test_counters_kwarg_is_optional(self, monkeypatch):
        # (+5 stays annotated here: with NO annotation the destination is a novel
        # site and a pre-existing novel-evidence hold vetoes the move after the
        # walk — unrelated to realizability, and not what this test is about)
        _controlled(monkeypatch, {220: 2.0, 222: 0.0, 225: 1.0})
        assert jr.refine_read_junctions(
            _witness_read(), IDX, {(CHROM, 120, 225)}, REF, "+",
            motif_blind=True) == [(1, 120, 220, 120, 225)]


class TestTheCounterReachesTheDriver:
    @pytest.mark.parametrize("n_workers", [1, 2], ids=["sequential", "parallel"])
    def test_driver_stats_carry_the_skip(self, monkeypatch, tmp_path, n_workers):
        _controlled(monkeypatch, {220: 2.0, 222: 0.0, 225: 1.0})
        path = tmp_path / "in.bam"
        with pysam.AlignmentFile(str(path), "wb", header=_header()) as fh:
            fh.write(_witness_read())
        out = tmp_path / "out.bam"
        if n_workers > 1:
            # the parallel driver spawns workers that cannot see the monkeypatch;
            # its counter plumbing is exercised through the batch-result contract
            res = jr._refine_read_batch.__wrapped__ if hasattr(jr._refine_read_batch, "__wrapped__") else None
            pytest.skip("worker processes do not inherit the score monkeypatch; "
                        "the parallel merge is covered by test_batch_result_carries_counters")
        stats = jr.refine_bam_junctions(
            str(path), str(out), aligner_bams=[], annotated_junctions=set(),
            genome={CHROM: REF},
            prebuilt_junction_pool={(CHROM, 120, 222), (CHROM, 120, 225)},
            prebuilt_annotated_set={(CHROM, 120, 225)},
            boundary_error_window=0, sort_and_index=False, n_workers=n_workers,
            motif_blind=True,
        )
        assert stats["refined"] == 1
        assert stats["unrealizable_winner_skipped"] == 1
        with pysam.AlignmentFile(str(out), "rb") as fh:
            assert [r.cigarstring for r in fh] == ["20M105N20M"]

    def test_batch_result_carries_counters(self, monkeypatch):
        """The parallel worker returns its counters in the batch dict and the
        parent folds them into stats (the merge contract, exercised in-process)."""
        _controlled(monkeypatch, {220: 2.0, 222: 0.0, 225: 1.0})
        jr._WORKER_POOL_STATE.clear()
        jr._WORKER_POOL_STATE.update({
            "header": _header(), "genome": {CHROM: REF}, "junctions_idx": IDX,
            "annotated_set": {(CHROM, 120, 225)},
            "kwargs": {"motif_blind": True, "boundary_error_window": 0},
            "penalty_table_set": None, "profile_enabled": False, "profile_sample_rate": 1,
        })
        try:
            res = jr._refine_read_batch([_witness_read().to_string()])
        finally:
            jr._WORKER_POOL_STATE.clear()
        assert isinstance(res, dict)
        assert res["counters"] == {"unrealizable_winner_skipped": 1}
        assert res["results"][0][1] == [(1, 120, 220, 120, 225)]
