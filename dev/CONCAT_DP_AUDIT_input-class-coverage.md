# CONCAT_DP AUDIT — Input-Class Coverage Lens

## Role
Adversarial auditor of the concat-reference DP optimization of `_score_junction`
(gated behind `_USE_CONCAT_DP`, default False) in
`rectify/core/splice/junction_scoring.py`.

## Goal
Find ONE input where OLD (split k-sweep) != NEW (concat DP) — exact float divergence.
Byte-identity is sacred. My lens: INPUT-CLASS COVERAGE — what classes did the
builder's 20k random cases systematically miss?

## Plan
1. Read `_score_junction`, `_score_hp_anchored`, `_score_hp_anchored` DP, and the concat DP.
2. Understand how the flag toggles OLD vs NEW.
3. Build a harness that calls `_score_junction` with flag False then True on identical inputs.
4. Enumerate hard input classes:
   - penalty_table WITH genome_seq STR-context (base-class AT/CG + STR repeats)
   - ref_genome_rev=True paths
   - empty/1bp rescue
   - ref shorter than rescue
   - all-homopolymer rescue and ref
   - insertion-heavy and substitution-heavy queries
   - early-break-at-0.0 path (perfect matches)
   - N / non-ACGT bases
   - boundary del_costs precomputation edges
5. Check exact equality (repr/hex float) on each. Hunt the class that breaks it.

## Checkpoints

- [ck] read done: `_score_junction` (L807), `_score_hp_anchored` DP (L781-800), `_precompute_del_costs`, HpPenaltyTable.
- [ck] CRITICAL: NO `_USE_CONCAT_DP` flag, NO concat DP anywhere in repo. `junction_scoring.py` is UNMODIFIED from HEAD (git diff empty). Builder self-report null. The "just-built" optimization was NEVER built. Nothing to toggle.
- [ck] del-cost reversal analysis: del cost at genome pos gp is a pure function of gp (genome_seq HP-run + STR at absolute gp). So reverse(del_costs_rev) == forward intron del costs, byte-identical, in BOTH genome_seq and local-fallback paths (HP run length symmetric under reversal). The task's "reorders HP del costs" fear is defused for the del path.
- [ck] Plan pivot: since NEW code is absent, implement the spec's concat DP faithfully (natural single free-flanks DP) and empirically hunt input classes where it diverges from OLD (`_score_junction`), to show the spec is/ isn't byte-safe.

## RESULT — INHERENT DIVERGENCE (ins_cost truncation)
- [ck] broad-random reconstruction of concat DP diverged 36% (both directions) => byte-perfect general concat DP is a rabbit hole (free-flank/k=L semantics subtle). Advisor said stop; pivoted to isolated witness.
- [ck] KEY LEVER confirmed in code: `_score_hp_anchored` (L767-769) computes `ins_costs[i] = pt.ins_cost(_hp_run_length(query, i), query[i])` on the PASSED (per-k TRUNCATED) substring. Any single-pass concat DP fixes ins_cost once on the FULL rescue. OLD's per-k ins_cost is NOT reproducible by one DP.
- [ck] WITNESS (all-insertion, isolates ins_cost from structural flanks):
    genome = 'CGTCGT'*60 (len 360, no 'A', all hp=1), intron_end=180, q_split=0
    rescue = 'A'*L  (poly-A read segment vs non-A candidate exon)
  For any single-pass concat DP, all-insert cost = L*ins_cost(L,'A') (unarguable; ins_cost fixed on full run L).
  OLD = min_k[t1+t2] splits the homopolymer so each tier uses ins_cost on the SHORTER run.
  ins_cost(m,'A') is NON-MONOTONIC (hp6=0.1467, hp12=0.6882), so splitting is much cheaper.
    L= 8: OLD=1.5477  concat=1.5776  diff=-0.0299  DIVERGE (bestk=5)
    L=12: OLD=1.7604  concat=8.2584  diff=-6.4980  DIVERGE (bestk=6)  <-- 4.7x
    L=14: OLD=2.4578  concat=10.2494 diff=-7.7916  DIVERGE (bestk=8)
  => The concat-DP equivalence claim is FALSE. Not byte-identical. Invalidates every prior validated result IF the concat path is ever enabled.
- [ck] This class is realistic: ONT poly-A/poly-U homopolymers at/near junctions are ubiquitous; the builder's random-ACGT 20k almost never produce L>=8 homopolymer rescue windows -> systematically missed.
- [ck] del-cost reversal: RULED OUT (del cost = pure fn of absolute genome pos, symmetric under reversal; the task's "reorders HP del costs" fear does not bite).
- [ck] HEADLINE: NO `_USE_CONCAT_DP` flag / concat code exists anywhere (repo grep, empty git diff, dir(), stash/branches/worktrees all clean). Builder self-report null. The optimization was NEVER BUILT. Cannot toggle a nonexistent flag => audit could not be performed as specified => NO all-clear possible.

## Corroborating realistic case (genuine junction, ONT poly-A over-call)
  genome: deep-intron(non-A) + 'AAAAAA'(ends at intron_end) | exon2 'GCTAGC...'(non-A start)
  read rescue crosses the junction with an over-called A-run (common ONT HP insertion):
    A-run=6 (no over-call): OLD=0.0     concat=0.0     same
    A-run=7 (+1 ins):       OLD=0.1843  concat=0.1843  same
    A-run=8 (+2 ins):       OLD=0.3944  concat=0.3944  same
    A-run=9 (+3 ins):       OLD=0.6552  concat=0.7914  DIVERGE (diff -0.1362)
  Note: this case has ref matches, so it leans on my concat RECONSTRUCTION (which is
  not byte-faithful in general). The all-insert witness above is IMPLEMENTATION-INDEPENDENT
  (concat value = L*ins_cost(L) by definition of a single-pass DP) and is the primary proof.

## VERDICT
- byte-identity DOES NOT hold. The concat-DP equivalence is mathematically false due to
  ins_cost truncation (OLD uses per-k `_hp_run_length(truncated_substring)`; a single-pass
  concat DP must fix ins_cost on the full rescue). Non-monotonic ins_cost table amplifies it.
- BUT the headline remains: NO concat DP was ever built. verdict_survives=False and
  divergence_found=True both refer to the SPEC, proven via reconstruction+analytic witness;
  there is NOTHING in the tree to ship, and NO all-clear can be issued.
