"""Phase-0 harness unit tests — the spec-mandated cases.

Runnable on M1 without the human data (validates the verdict machinery itself):
    python dev/round2_phase0/test_phase0.py        # plain runner, prints PASS/FAIL
    pytest dev/round2_phase0/test_phase0.py         # or via pytest

Covers, from SPEC_round2_cdna_scoring_and_liftover_DESIGNER2 + the MASTER spec:
  - lift-over plus/minus strand + N-op insertion + invariants (§2.5, incl. the CRITICAL
    minus-strand orientation invariant that catches a missed revcomp);
  - §1.3 trivial-win defeat (clean spliced read: cDNA ties -> genome kept);
  - micro-exon WIN (genome soft-clips a micro-exon it can't seed; cDNA recovers -> wins,
    attributed 'placed_micro_exon');
  - LOSE #1 intron-retention (cDNA worse -> genome kept);
  - LOSE #3 no-shrink veto (cDNA lower abs HP-ED but fewer aligned bases -> kept);
  - BLOCKER-1 (§7): a lifted candidate with sub-K anchor AND five_prime_rescued=1 must be
    REJECTED -- the gate is on the RAW anchor, never on _effective_chimera_ok;
  - gate-off safety precondition (§3.2): Round-2 refuses to run with the gate at 0.
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover import (  # noqa: E402
    Block, CdnaModel, lift, revcomp,
    assert_query_length, reference_span, assert_no_n_merge, assert_orientation_invariant,
)
from score_phase0 import (  # noqa: E402
    Candidate, decide, anchor_gate_passes, run_verdict, _NO_JUNCTION_ANCHOR,
)
from rectify.core.consensus.corrected_consensus import (  # noqa: E402
    _cigar_hp_edit_distance, _cigar_aligned_bases, _cigar_min_junction_anchor,
)

M, I, D, N, S, EQ, X = 0, 1, 2, 3, 4, 7, 8


def _genome(seed=0, length=400):
    random.seed(seed)
    return {"T1": "".join(random.choice("ACGT") for _ in range(length))}


# ───────────────────────────── lift-over tests ─────────────────────────────

def test_plus_strand_single_junction():
    g = _genome()
    seq = g["T1"]
    # exon1 g[100:120], intron g[120:170] (50bp), exon2 g[170:190]
    cdna = CdnaModel(
        cdna_id="plus", chrom="T1", strand="+",
        blocks=[Block(0, 100, 120), Block(20, 170, 190)], length=40,
    )
    cdna.validate(g)
    read = seq[100:120] + seq[170:190]            # perfect-match read across the junction
    lifted = lift(cdna, [(M, 40)], 0, read)
    assert lifted.cigartuples == [(M, 20), (N, 50), (M, 20)], lifted.cigartuples
    assert lifted.reference_start == 100
    assert not lifted.is_reverse
    assert_query_length(lifted, len(read))
    assert reference_span(lifted) == 90           # 20 + 50 + 20
    assert_no_n_merge(lifted)
    f, l = assert_orientation_invariant(lifted, g)
    assert f and l
    # production scorers on the LiftedRead:
    assert _cigar_hp_edit_distance(lifted, g, None) == 0.0
    assert _cigar_aligned_bases(lifted) == 40
    assert _cigar_min_junction_anchor(lifted, g) == 20


def test_minus_strand_single_junction():
    g = _genome()
    seq = g["T1"]
    # minus gene: RNA-sense cDNA = revcomp; 5' exon is the higher genome coord (170:190)
    cdna = CdnaModel(
        cdna_id="minus", chrom="T1", strand="-",
        blocks=[Block(0, 170, 190), Block(20, 100, 120)], length=40,
    )
    cdna.validate(g)
    read = revcomp(seq[170:190]) + revcomp(seq[100:120])   # RNA-sense, as it aligns to cDNA
    lifted = lift(cdna, [(M, 40)], 0, read)
    assert lifted.is_reverse
    assert lifted.reference_start == 100
    assert lifted.cigartuples == [(M, 20), (N, 50), (M, 20)], lifted.cigartuples
    # the lifted query must be reference-forward = genome[100:120]+genome[170:190]
    assert lifted.query_sequence == seq[100:120] + seq[170:190]
    assert_query_length(lifted, len(read))
    f, l = assert_orientation_invariant(lifted, g)
    assert f and l, "minus-strand orientation invariant failed (missed revcomp?)"
    assert _cigar_hp_edit_distance(lifted, g, None) == 0.0
    assert _cigar_min_junction_anchor(lifted, g) == 20


def test_orientation_invariant_catches_missed_revcomp():
    """A deliberately un-revcomp'd minus-strand query must FAIL the orientation invariant
    while still passing the (revcomp-blind) query-length + ref-span invariants (§2.5)."""
    g = _genome()
    seq = g["T1"]
    cdna = CdnaModel("minus", "T1", "-", [Block(0, 170, 190), Block(20, 100, 120)], 40)
    read = revcomp(seq[170:190]) + revcomp(seq[100:120])
    good = lift(cdna, [(M, 40)], 0, read)
    # simulate the bug: keep the cigar/coords but forget to revcomp the query
    from liftover import LiftedRead
    bad = LiftedRead(good.cigartuples, good.reference_name, good.reference_start,
                     read, is_reverse=True)
    assert_query_length(bad, len(read))                 # revcomp-blind: still passes
    assert reference_span(bad) == reference_span(good)  # revcomp-blind: still passes
    f, l = assert_orientation_invariant(bad, g)
    assert not (f and l), "orientation invariant should have caught the missed revcomp"


def _4exon_plus(g):
    """4-exon plus-strand gene incl. a 9 bp micro-exon. Returns (CdnaModel, cdna_seq)."""
    seq = g["T1"]
    blocks = [Block(0, 100, 120), Block(20, 200, 209), Block(29, 300, 320), Block(49, 400, 420)]
    cdna = CdnaModel("plus4", "T1", "+", blocks, 69)
    cdna.validate(g)
    cseq = seq[100:120] + seq[200:209] + seq[300:320] + seq[400:420]
    return cdna, cseq


def _4exon_minus(g):
    """4-exon minus-strand gene (RNA-sense, descending genome coords)."""
    seq = g["T1"]
    blocks = [Block(0, 400, 420), Block(20, 300, 320), Block(40, 200, 209), Block(49, 100, 120)]
    cdna = CdnaModel("minus4", "T1", "-", blocks, 69)
    cdna.validate(g)
    cseq = revcomp(seq[400:420]) + revcomp(seq[300:320]) + revcomp(seq[200:209]) + revcomp(seq[100:120])
    return cdna, cseq


def test_boundary_break_emits_n():
    """ADVISOR REGRESSION: ops that break EXACTLY at a block boundary must still emit the
    intron. [(M,20),(X,1),(M,19)] across two non-contiguous blocks => N must appear."""
    g = _genome(length=500)
    cdna = CdnaModel("plus2", "T1", "+", [Block(0, 100, 120), Block(20, 170, 190)], 40)
    seq = g["T1"]
    read = list(seq[100:120] + seq[170:190])
    # force a mismatch at the first base of exon2 (the X)
    read[20] = "A" if seq[170] != "A" else "C"
    read = "".join(read)
    lifted = lift(cdna, [(M, 20), (X, 1), (M, 19)], 0, read)
    assert lifted.cigartuples == [(M, 20), (N, 50), (X, 1), (M, 19)], lifted.cigartuples
    assert any(op == N for op, _ in lifted.cigartuples), "intron dropped at op/block-boundary coincidence"
    assert_query_length(lifted, len(read))
    assert _cigar_hp_edit_distance(lifted, g, None) == 1.0   # the single X


def test_multi_junction_three_n():
    """>=3 N-ops: a clean read across a 4-exon gene (incl. a 9 bp micro-exon)."""
    g = _genome(length=500)
    cdna, cseq = _4exon_plus(g)
    lifted = lift(cdna, [(M, 69)], 0, cseq)
    assert lifted.cigartuples == [(M, 20), (N, 80), (M, 9), (N, 91), (M, 20), (N, 80), (M, 20)], lifted.cigartuples
    assert sum(1 for op, _ in lifted.cigartuples if op == N) == 3
    assert reference_span(lifted) == (20 + 80 + 9 + 91 + 20 + 80 + 20)
    assert _cigar_hp_edit_distance(lifted, g, None) == 0.0
    assert _cigar_min_junction_anchor(lifted, g) == 9   # weakest flank = the 9 bp micro-exon


def test_insertion_at_junction_e6():
    """E6: an insertion abutting a junction attaches to the donor (5') side: `=, I, N, =`."""
    g = _genome(length=500)
    cdna, cseq = _4exon_plus(g)
    read = cseq[:20] + "GGG" + cseq[20:]            # 3 inserted bases at the exon0/exon1 boundary
    lifted = lift(cdna, [(M, 20), (I, 3), (M, 49)], 0, read)
    # the I must sit BEFORE the first N (donor side), never inside the intron
    ops = lifted.cigartuples
    i_idx = next(k for k, (o, _) in enumerate(ops) if o == I)
    n_idx = next(k for k, (o, _) in enumerate(ops) if o == N)
    assert i_idx < n_idx and ops[i_idx] == (I, 3) and ops[i_idx + 1][0] == N
    assert sum(1 for o, _ in ops if o == N) == 3
    assert_query_length(lifted, len(read))


def test_deletion_at_junction_e5():
    """E5: a deletion at the exon end stays a separate D op, never merged into the N."""
    g = _genome(length=500)
    cdna, cseq = _4exon_plus(g)
    # read deletes the last 2 bases of exon0 then continues; cDNA-local: M18 D2 M49
    read = cseq[:18] + cseq[20:]
    lifted = lift(cdna, [(M, 18), (D, 2), (M, 49)], 0, read)
    ops = lifted.cigartuples
    assert (D, 2) in ops
    # the D is not absorbed into an N (no op of type N with length != a real intron gap)
    n_lens = sorted(l for o, l in ops if o == N)
    assert n_lens == [80, 80, 91], n_lens
    assert_query_length(lifted, len(read))


def test_minus_multijunction_e11():
    """E6xE11: minus-strand multi-junction lift — reversal + revcomp + N-insertion together;
    the orientation invariant (the one revcomp-blind check can't catch) must hold."""
    g = _genome(length=500)
    cdna, cseq = _4exon_minus(g)
    lifted = lift(cdna, [(M, 69)], 0, cseq)
    assert lifted.is_reverse
    assert sum(1 for op, _ in lifted.cigartuples if op == N) == 3
    # genome CIGAR reads low->high: M20 N80 M9 N91 M20 N80 M20 (symmetric layout here)
    assert lifted.cigartuples == [(M, 20), (N, 80), (M, 9), (N, 91), (M, 20), (N, 80), (M, 20)], lifted.cigartuples
    assert lifted.reference_start == 100
    f, l = assert_orientation_invariant(lifted, g)
    assert f and l, "minus-strand multi-junction orientation invariant failed"
    assert _cigar_hp_edit_distance(lifted, g, None) == 0.0
    assert _cigar_min_junction_anchor(lifted, g) == 9


# ───────────────────────────── win-guard tests ─────────────────────────────

def _cand(aligner, hp, ab, anchor, **kw):
    return Candidate(read_id="r", aligner=aligner, hp_ed=hp, aligned_bases=ab,
                     min_junction_anchor=anchor, **kw)


def test_trivial_win_defeat():
    """§1.3: a read the genome already splices correctly -> cDNA ties -> genome kept."""
    gwin = _cand("minimap2", 0.0, 215, 120, n_junctions=1)
    cdna = _cand("cdna_round2", 0.0, 215, 95, n_junctions=1)
    d = decide(cdna, gwin, min_junction_anchor_bp=10)
    assert not d.cdna_wins and d.reason == "not_strictly_better"


def test_micro_exon_win_attributed():
    """WIN #1: genome soft-clips a micro-exon (HP-ED 12); cDNA recovers it (HP-ED 0)."""
    # micro-exon skip: genome aligner emits mismatches over the 9 bp exon (HP-ED 12),
    # one fused junction; cDNA places the extra micro-exon junction cleanly.
    gwin = _cand("minimap2", 12.0, 200, 120, n_junctions=1, genome_softclip_bp=9)
    cdna = _cand("cdna_round2", 0.0, 209, 15, n_junctions=2, genome_softclip_bp=9)
    ultra = _cand("uLTRA", 10.0, 201, 120, n_junctions=1)
    d = decide(cdna, gwin, ultra=ultra, min_junction_anchor_bp=10)
    assert d.cdna_wins and d.cause == "placed_micro_exon"
    assert d.beats_ultra is True


def test_intron_retention_lose():
    """LOSE #1: retained-intron read; library lacks the intron -> cDNA worse -> genome kept."""
    gwin = _cand("minimap2", 5.0, 300, _NO_JUNCTION_ANCHOR, n_junctions=0)
    cdna = _cand("cdna_round2", 15.0, 250, 40, n_junctions=1)
    d = decide(cdna, gwin, min_junction_anchor_bp=10)
    assert not d.cdna_wins and d.reason == "not_strictly_better"


def test_no_shrink_veto():
    """LOSE #3: cDNA has lower ABSOLUTE HP-ED but aligned fewer bases -> vetoed."""
    gwin = _cand("minimap2", 3.0, 200, 120, n_junctions=1)
    cdna = _cand("cdna_round2", 2.0, 80, 40, n_junctions=1)   # 2 <= 3-1 passes (2)
    d = decide(cdna, gwin, min_junction_anchor_bp=10)
    assert not d.cdna_wins and d.reason == "no_shrink_veto"


def test_blocker1_rescued_does_not_bypass_gate():
    """BLOCKER-1 (§7): sub-K anchor + five_prime_rescued=1 must be REJECTED.
    The gate is on the RAW anchor; the Cat3 exemption must NOT switch it off."""
    gwin = _cand("minimap2", 12.0, 200, 120, n_junctions=1)
    cdna = _cand("cdna_round2", 0.0, 209, 3, n_junctions=1, five_prime_rescued=1)
    assert not anchor_gate_passes(cdna, 10)        # raw anchor 3 < 10 -> fails
    d = decide(cdna, gwin, min_junction_anchor_bp=10)
    assert not d.cdna_wins and d.reason == "cdna_failed_anchor_gate"


def test_gate_off_safety_precondition():
    """§3.2: Round-2 must refuse to run with the anchor gate off (yeast default)."""
    cdna = _cand("cdna_round2", 0.0, 209, 3, n_junctions=1)
    try:
        anchor_gate_passes(cdna, 0)
    except ValueError:
        return
    raise AssertionError("gate-off precondition did not raise")


def test_verdict_aggregation_go():
    """A mixed batch where the majority of wins are attributable + beat uLTRA -> GO."""
    triples = []
    # 3 genuine micro-exon wins (attributed, beat uLTRA)
    for _ in range(3):
        triples.append((
            _cand("cdna_round2", 0.0, 209, 15, n_junctions=2, genome_softclip_bp=25),
            _cand("minimap2", 12.0, 200, 120, n_junctions=1, genome_softclip_bp=25),
            _cand("uLTRA", 10.0, 201, 120, n_junctions=1),
        ))
    # 2 reads kept on genome (ties)
    for _ in range(2):
        triples.append((
            _cand("cdna_round2", 0.0, 215, 95, n_junctions=1),
            _cand("minimap2", 0.0, 215, 120, n_junctions=1),
            _cand("uLTRA", 0.0, 215, 120, n_junctions=1),
        ))
    v = run_verdict(triples, min_junction_anchor_bp=10)
    assert v.n_cdna_wins == 3
    assert v.trivial_win_leak == 0.0
    assert v.ultra_specific_win_frac == 1.0
    assert v.net_hp_ed_reduction == 36.0
    assert v.go() is True


def _run_all():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    fails = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except Exception as e:
            fails += 1
            print(f"FAIL  {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - fails}/{len(tests)} passed")
    return fails


if __name__ == "__main__":
    sys.exit(1 if _run_all() else 0)
