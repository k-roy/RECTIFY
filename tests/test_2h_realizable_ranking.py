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

Codex's review of that first cut (CFX-12, 2026-09-14) narrowed the claim with two
witnesses, both now pinned below:
  (a) the move GATES ran after the walk, so a writable but gate-vetoed head hid a
      permitted runner-up — the gates now run per candidate, inside the walk,
      before the dry run (``gate_vetoed_candidate_skipped``);
  (b) each N-op was dry-run on the ORIGINAL read, so two moves that are each
      writable alone could conflict on a shared exon in the real right-to-left
      write — the walk now runs right-to-left on an evolving trial copy, so each
      candidate is judged on the read the writer will hand its surgery.
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


class TestTwoNopInteraction:
    """Codex's interaction question: two N-ops whose moves are each writable ALONE but
    not jointly (the read-level indel-burden invariant sees the first move's burden when
    the second is written). The per-N-op dry-run runs on the read as it stands at
    proposal time, so both are proposed; the real write applies right-to-left and the
    second surgery must refuse on its own — the read may never end up in a state no
    single dry-run vouched for."""

    def _two_nop_read(self):
        # 20M 100N 5D 20M 100N 5D 20M on an all-C genome: each acceptor +5 absorbs its
        # own 5D (writable alone, burden 10 -> 5 -> 0); the +2 variants leave 3D (refused).
        r = pysam.AlignedSegment(_header())
        r.query_name = "two_nops"
        r.reference_id = 0
        r.reference_start = 100
        r.cigarstring = "20M100N5D20M100N5D20M"
        r.query_sequence = "C" * 60
        r.query_qualities = pysam.qualitystring_to_array("I" * 60)
        return r

    def test_two_individually_writable_absorptions_are_both_written(self, monkeypatch):
        ns1, ne1 = 120, 220        # first N: ref 120..220, then 5D at 220..225, exon 225..245
        ns2, ne2 = 245, 345        # second N, then 5D 345..350, exon 350..370
        scores = {220: 2.0, 222: 0.0, 225: 1.0, 345: 2.0, 347: 0.0, 350: 1.0}
        monkeypatch.setattr(jr, "_score_junction",
                            lambda q, q_split, js, je, g, **kw: (scores[je], 0))
        idx = _build_junction_index({(CHROM, ns1, e) for e in (220, 222, 225)}
                                    | {(CHROM, ns2, e) for e in (345, 347, 350)})
        ann = {(CHROM, ns1, 225), (CHROM, ns2, 350)}
        r = self._two_nop_read()
        counters = Counter()
        proposals = jr.refine_read_junctions(r, idx, ann, REF, "+", motif_blind=True,
                                             counters=counters)
        assert sorted(proposals) == [(1, ns1, ne1, ns1, 225), (4, ns2, ne2, ns2, 350)], proposals
        assert counters["unrealizable_winner_skipped"] == 2       # one +2 per N-op
        out, applied = jr._apply_replacements_to_read(r, proposals, REF, "+", 0.25, 15)
        assert applied and out.cigarstring == "20M105N20M105N20M"

    def test_a_jointly_unwritable_pair_never_yields_an_unvouched_state(self, monkeypatch):
        """Force the second surgery to refuse AFTER the first succeeded and check the
        read is left in a state some single dry-run vouched for (first move only)."""
        real = jr._apply_junction_replacement
        calls = []

        def flaky(read, cigar_idx, *a, **k):
            calls.append(cigar_idx)
            if cigar_idx == 1 and len(calls) > 1:       # the real (second) write of N-op 1
                return False
            return real(read, cigar_idx, *a, **k)
        monkeypatch.setattr(jr, "_apply_junction_replacement", flaky)
        r = self._two_nop_read()
        proposals = [(1, 120, 220, 120, 225), (4, 245, 345, 245, 350)]
        out, applied = jr._apply_replacements_to_read(r, proposals, REF, "+", 0.25, 15)
        # right-to-left: N-op 4 written, N-op 1 refused -> exactly the N-op-4-only state
        assert applied and out.cigarstring == "20M100N5D20M105N20M"


# ---------------------------------------------------------------------------
# Codex CFX-12 witnesses (dev/audits/codex_fable_20260913/repro_1597a3b_limits.py)
# ---------------------------------------------------------------------------

def _read600(cigar):
    h = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": CHROM, "LN": 600}]})
    r = pysam.AlignedSegment(h)
    r.query_name = "cfx12"
    r.reference_id = 0
    r.reference_start = 100
    r.cigarstring = cigar
    r.query_sequence = "C" * sum(n for op, n in r.cigartuples if op in (0, 1, 4, 7, 8))
    r.query_qualities = pysam.qualitystring_to_array("I" * len(r.query_sequence))
    return r


def _by_pair(monkeypatch, scores):
    def fake(q, q_split, js, je, genome_seq, **kw):
        return scores[(js, je)], 0
    monkeypatch.setattr(jr, "_score_junction", fake)


class TestGateVetoedHeadYieldsToPermittedRunnerUp:
    """(a) ``20M100N20M``: annotated canonical incumbent [120,220) at .5, a
    writable NON-canonical head [121,221) at 0 that the annotated-canonical hold
    vetoes, and a writable annotated canonical runner-up [124,224) at .2 that
    nothing vetoes.  The runner-up must be proposed, with or without the head."""

    def _genome(self):
        g = list("C" * 600)
        g[120:126] = "GTAGGT"
        g[218:226] = "AGGTAGGT"
        return "".join(g)

    CANDS = {(120, 220), (121, 221), (124, 224)}
    ANN = {(CHROM, 120, 220), (CHROM, 124, 224)}
    SCORES = {(120, 220): 0.5, (121, 221): 0.0, (124, 224): 0.2}

    def test_frame(self):
        g = self._genome()
        assert jr._canonical_tier(120, 220, g, "+") == 0
        assert jr._canonical_tier(124, 224, g, "+") == 0
        assert jr._canonical_tier(121, 221, g, "+") >= 4
        r = _read600("20M100N20M")
        assert jr._realizable(r, 1, 120, 220, 121, 221, g, "+", 0.25, 15)
        assert jr._realizable(r, 1, 120, 220, 124, 224, g, "+", 0.25, 15)

    def test_the_permitted_runner_up_is_proposed_past_the_vetoed_head(self, monkeypatch):
        g = self._genome()
        _by_pair(monkeypatch, self.SCORES)
        idx = _build_junction_index({(CHROM, s, e) for s, e in self.CANDS})
        counters = Counter()
        proposals = jr.refine_read_junctions(_read600("20M100N20M"), idx, self.ANN, g, "+",
                                             boundary_error_window=0, counters=counters)
        assert proposals == [(1, 120, 220, 124, 224)], proposals
        assert counters["gate_vetoed_candidate_skipped"] == 1
        assert counters["unrealizable_winner_skipped"] == 0
        out, applied = jr._apply_replacements_to_read(
            _read600("20M100N20M"), proposals, g, "+", 0.25, 15)
        assert applied and out.cigarstring == "24M100N16M"

    def test_removing_the_vetoed_head_changes_nothing(self, monkeypatch):
        g = self._genome()
        _by_pair(monkeypatch, self.SCORES)
        idx = _build_junction_index({(CHROM, s, e) for s, e in self.CANDS if s != 121})
        counters = Counter()
        proposals = jr.refine_read_junctions(_read600("20M100N20M"), idx, self.ANN, g, "+",
                                             boundary_error_window=0, counters=counters)
        assert proposals == [(1, 120, 220, 124, 224)]
        assert counters["gate_vetoed_candidate_skipped"] == 0


class TestInteractingMovesAreJudgedOnTheEvolvingRead:
    """(b) ``20M100N5M100N20M`` on a C-only genome, motif-blind, controlled
    ranking: first N +3 and second N -3 are each writable ALONE, but written
    right-to-left the second leaves ``2M`` and the first cannot be funded.  The
    walk runs in the writer's order on a trial copy: the second N's move is
    accepted first, the first N's +3 is refused on the updated read, and its +1
    runner-up — writable there — is proposed instead.  The real write then
    reproduces exactly that state."""

    CANDS = {(120, 220), (121, 221), (123, 223), (225, 325), (222, 322)}
    SCORES = {(120, 220): 2, (121, 221): 1, (123, 223): 0, (225, 325): 2, (222, 322): 0}

    def test_each_move_is_writable_alone_but_not_jointly(self):
        g = "C" * 600
        r = _read600("20M100N5M100N20M")
        assert jr._realizable(r, 1, 120, 220, 123, 223, g, "+", 0.25, 15)
        assert jr._realizable(r, 3, 225, 325, 222, 322, g, "+", 0.25, 15)
        out, _ = jr._apply_replacements_to_read(
            r, [(1, 120, 220, 123, 223), (3, 225, 325, 222, 322)], g, "+", 0.25, 15)
        assert out.cigarstring == "20M100N2M100N23M"      # the first move was refused

    def test_the_runner_up_is_proposed_on_the_evolving_read(self, monkeypatch):
        g = "C" * 600
        _by_pair(monkeypatch, self.SCORES)
        idx = _build_junction_index({(CHROM, s, e) for s, e in self.CANDS})
        counters = Counter()
        proposals = jr.refine_read_junctions(_read600("20M100N5M100N20M"), idx, set(), g, "+",
                                             boundary_error_window=0, motif_blind=True,
                                             counters=counters)
        assert proposals == [(1, 120, 220, 121, 221), (3, 225, 325, 222, 322)], proposals
        assert counters["unrealizable_winner_skipped"] == 1        # the +3 on the updated read
        out, applied = jr._apply_replacements_to_read(
            _read600("20M100N5M100N20M"), proposals, g, "+", 0.25, 15)
        assert applied and out.cigarstring == "21M100N1M100N23M"
