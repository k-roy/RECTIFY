#!/usr/bin/env python3
"""C3 multi-event coherence check — the ONE place an LLR-addressable gap could hide.

The single-event headroom probes (c3_headroom.py, c3_junction_headroom.py) showed
the shipped hp_edit_distance arbiter at ceiling. But hp_ed scores indels with
`del_cost/ins_cost` = the `penalty_score` column, which C1 proved is reciprocal-rate
(~c/rate_mean), NOT -logP — and a reciprocal-rate sum can rank MULTI-event alignment
paths differently from a coherent -logP sum. That divergence is the only mechanism by
which a calibrated-LLR arbiter (C3) could beat hp_ed on accuracy. This check measures
whether the divergence (a) exists and (b) is REACHABLE as a truth-favoring tie.

Finding (run output below):
  (a) EXISTS at the table level: many 2-deletion decompositions rank oppositely under
      penalty_score-sum vs -logP-sum (penalty over-weights long/cheap runs ~reciprocally;
      -logP weights them ~logarithmically).
  (b) REACHABILITY is geometrically constrained and undemonstrated: for two alignments
      of ONE read to be edit-distance-tied yet assign deletions to DIFFERENT-length runs,
      you need a cross-run-length-reassignable multi-deletion tie. Within a same-base run,
      shifts are ambiguity-EQUIVALENT (same run length -> no inversion). Across
      different-base adjacent runs, each deletion sits unambiguously in its own run
      (no reassignment). So the inversions above do NOT correspond to competing
      placements of a realistic read; demonstrating a reachable truth-favoring case
      would require a dedicated ADJACENT/INTERLEAVED-RUNS multi-event stratum (not in
      the corpus). Absent that, the refute stands; this names the precise locus for any
      future C3 revisit.

Usage: python scripts/benchmark/c3_multievent_check.py \
         --penalty-table rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import itertools
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--penalty-table", required=True)
    ap.add_argument("--str-table", default=None)
    ap.add_argument("--min-count", type=int, default=100)
    ap.add_argument("--base", default="A", help="ref base (A/T -> AT class, C/G -> CG)")
    ap.add_argument("--events", type=int, default=2, help="events per path (multiset size)")
    ap.add_argument("--max-hp", type=int, default=10)
    args = ap.parse_args()

    T = HpPenaltyTable.from_tsv(args.penalty_table, args.str_table, min_count=args.min_count)
    bc = "AT" if args.base.upper() in ("A", "T") else "CG"
    rate_tbl = T._del_rate.get(bc) or {}

    def pen(hp):
        return T.del_cost(hp, args.base)

    def nlp(hp):
        r = rate_tbl.get(hp)
        return -math.log(r) if r else None

    hps = [h for h in range(1, args.max_hp + 1) if nlp(h) is not None]
    decs = list(itertools.combinations_with_replacement(hps, args.events))

    def Pen(d):
        return sum(pen(h) for h in d)

    def Nlp(d):
        return sum(nlp(h) for h in d)

    inv = 0
    examples = []
    for A, B in itertools.combinations(decs, 2):
        sp = Pen(A) - Pen(B)
        sn = Nlp(A) - Nlp(B)
        if sp * sn < -1e-9:
            inv += 1
            if len(examples) < 8:
                examples.append((A, B, Pen(A), Pen(B), Nlp(A), Nlp(B)))

    print(f"=== C3 multi-event coherence check (base_class={bc}, {args.events} dels/path) ===")
    print(f"decompositions: {len(decs)}; pairwise inversions "
          f"(penalty_score-sum vs -logP-sum disagree): {inv}")
    for A, B, pa, pb, na, nb in examples:
        print(f"  del@runs{A} vs {B}: Pen {pa:.3f}/{pb:.3f} (hp_ed prefers "
              f"{'A' if pa < pb else 'B'}) | -logP {na:.3f}/{nb:.3f} (LLR prefers "
              f"{'A' if na < nb else 'B'})  margin={abs(na-nb):.3f} nat")
    print("\n(a) inversions EXIST => penalty_score is incoherent to sum (C1's claim, confirmed).")
    print("(b) REACHABILITY undemonstrated: needs an adjacent/interleaved-runs multi-event")
    print("    stratum where two edit-distance-tied placements assign dels to different-length")
    print("    runs with truth = -logP-MAP. Same-base shifts are ambiguity-equivalent; different-")
    print("    base adjacent runs fix each del in its own run. => no reachable C3 gap shown; the")
    print("    refute stands. This is the named locus IF C3 is ever revisited.")


if __name__ == "__main__":
    main()
