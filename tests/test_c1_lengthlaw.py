"""C1 HP-length-law gap-cost: mechanism tests (NOT the scientific ablation).

Asserts the wiring is correct and SAFE:
  1. penalty_table=None is byte-identical to the legacy DP (the Cat3 / junction-
     rescue regression guard);
  2. the gap-open delta is a baseline-anchored log-odds (zero at hp=1, monotone
     up with run length) computed from rate_mean, NOT from the mislabeled
     penalty_score column;
  3. a strong in-run discount pulls the deletion INTO the homopolymer run.

The scientific claim (does the discount improve placement vs a matched baseline,
without hallucinating indels) is the separate ablation, scripts/benchmark/c1_ablation.py.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from rectify.core.align.local_aligner import (  # noqa: E402
    align_exon_block_global, _homopolymer_run_len,
)
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa: E402


def _synthetic_table():
    """HpPenaltyTable with rate_mean RISING with hp (like real DRS deletions)."""
    del_rate = {"AT": {1: 0.005, 4: 0.02, 8: 0.06, 12: 0.12}, "CG": {1: 0.004}}
    ins_rate = {"AT": {1: 0.003, 12: 0.01}, "CG": {1: 0.003}}
    # penalty_score tables (reciprocal heuristic) — present but UNUSED by the delta
    pen = {"AT": {1: 0.44, 4: 0.17}, "CG": {1: 0.5}}
    return HpPenaltyTable(pen, pen, default_ins=1.25,
                          del_rate_tables=del_rate, ins_rate_tables=ins_rate)


def _indel_positions(cig):
    """Ref positions of D ops in a (op,len) CIGAR (ref_start=0)."""
    pos, rpos = [], 0
    for op, ln in cig:
        if op == 0:       # M
            rpos += ln
        elif op == 2:     # D
            pos.append((rpos, ln)); rpos += ln
        # I (op==1) consumes no ref
    return pos


def test_run_len_helper():
    seq = "CCCC" + "A" * 7 + "GG"
    assert _homopolymer_run_len(seq, 5) == (7, "A")
    assert _homopolymer_run_len(seq, 0)[1] == "C" and _homopolymer_run_len(seq, 0)[0] == 4
    assert _homopolymer_run_len(seq, -1) == (0, "")


def test_none_is_byte_identical():
    """penalty_table=None (and the default) must not change any CIGAR — Cat3 guard."""
    cases = [
        ("ACGTACGT", "ACGTACGT"),                       # clean, no HP
        ("ACGT", "ACGTACGT"),                           # deletion, no HP
        ("CCCCAAAAAAAGGGG", "CCCCAAAAAAAAAGGGG"),        # HP undercall
        ("CCCCAAAAAAATGGGG", "CCCCAAAAAAAAAGGGG"),       # HP + boundary sub
    ]
    for q, r in cases:
        base = align_exon_block_global(q, r, chrom_ref=r)
        none = align_exon_block_global(q, r, chrom_ref=r, penalty_table=None)
        assert base == none, (q, r, base, none)
        # and identical to the no-chrom_ref legacy call when chrom_ref omitted
        assert align_exon_block_global(q, r) == align_exon_block_global(q, r, penalty_table=None)


def test_delta_baseline_anchored_and_monotone():
    pt = _synthetic_table()
    assert pt.del_open_delta(1, "A") == 0.0                  # anchored at hp=1
    assert pt.del_open_delta(4, "A") == math.log(0.02 / 0.005)   # ln 4
    assert pt.del_open_delta(12, "A") == math.log(0.12 / 0.005)  # ln 24
    # monotone increasing with run length (rate rises -> cheaper gap)
    d = [pt.del_open_delta(hp, "A") for hp in (1, 4, 8, 12)]
    assert d == sorted(d) and d[0] == 0.0 and d[-1] > 0
    # lam scales it
    assert pt.del_open_delta(12, "A", lam=2.0) == 2.0 * pt.del_open_delta(12, "A")
    # no rate table -> zero delta (legacy file safety)
    empty = HpPenaltyTable({"AT": {1: 0.4}}, {"AT": {1: 0.4}})
    assert empty.del_open_delta(8, "A") == 0.0


def test_strong_discount_pulls_deletion_into_run():
    """With a large in-run discount, a 2-base undercall in a long A-run plus a
    boundary substitution is resolved as an IN-RUN deletion (not absorbed
    out-of-run). The legacy arm may place it elsewhere; the law arm must keep
    the deletion inside the run span [10, 22)."""
    left = "CGCGCGCGCG"          # 10 bp, ends in G
    right = "GCGCGCGCGC"         # 10 bp
    ref = left + "A" * 12 + right            # run at ref [10, 22)
    # query: run shortened to 10 A's + a boundary substitution (first 'right' G->T)
    query = left + "A" * 10 + "T" + right[1:]
    pt = _synthetic_table()
    law = align_exon_block_global(query, ref, chrom_ref=ref,
                                  penalty_table=pt, lam=5.0)
    dels = _indel_positions(law)
    assert dels, ("law produced no deletion", law)
    # every deleted base lies within the run span [10, 22)
    for pos, ln in dels:
        for p in range(pos, pos + ln):
            assert 10 <= p < 22, (f"deletion at ref {p} is OUT of the run", law)
    # consumes query and ref fully (global invariant)
    qcons = sum(l for op, l in law if op in (0, 1))
    rcons = sum(l for op, l in law if op in (0, 2))
    assert qcons == len(query) and rcons == len(ref), law


if __name__ == "__main__":
    import inspect
    g = dict(globals())
    tests = sorted((n, f) for n, f in g.items()
                   if n.startswith("test_") and callable(f)
                   and not inspect.signature(f).parameters)
    for n, f in tests:
        f(); print(f"  PASS {n}")
    print(f"all {len(tests)} C1 length-law tests passed")
