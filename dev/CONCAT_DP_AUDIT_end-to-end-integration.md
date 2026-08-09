# CONCAT_DP_AUDIT — END-TO-END INTEGRATION LENS

## Task
Adversarial audit of concat-reference DP optimization of `_score_junction`
(gated behind `_USE_CONCAT_DP`, default False) in
`rectify/core/splice/junction_scoring.py`.

My lens: END-TO-END INTEGRATION. Scalar-score equality is necessary but not
sufficient. Run the ACTUAL refiner OLD vs NEW on a real workload and compare
refined outputs byte-for-byte (chosen (js,je), CIGARs, emitted BAM), the
(score, delta) tuple, early-exit interaction, and profile counters.

## Plan
1. Read junction_scoring.py — understand _score_junction, _score_hp_anchored,
   the _USE_CONCAT_DP flag, and how the flag toggles behavior.
2. Find the refiner entry point (_make_arm_e.py / refine_bam_junctions).
3. Locate sim panel BAM: scripts/benchmark/noncanon_sim/mix_fair_out/aligned.sorted.bam
4. Build a harness: refine flag=False -> outputs A; refine flag=True -> outputs B.
5. Compare A vs B byte-for-byte (BAM, chosen js/je, CIGARs, score/delta, counters).
6. Run full refiner+guard test suite with flag ON.
7. Persist any divergence input+values IMMEDIATELY.

## Checkpoints

- CHECKPOINT: orientation done. CRITICAL: `_USE_CONCAT_DP` flag + concat DP code
  DO NOT EXIST. grep of junction_scoring.py = 0 hits; git diff HEAD on that file = empty;
  git stash empty; all 11 worktrees = 0 concat-hits. Builder self-report null.
  => There is NO builder NEW implementation to toggle. Task #22 (refiner perf) is still
  `pending` in the task list. The optimization was never built in this tree.
- DECISION: since I hold the exact spec, I will (a) load OLD _score_junction as-is,
  (b) construct a faithful concat-DP NEW per the spec as a standalone fn, monkeypatch it
  in-process (NOT editing the file), and (c) run the END-TO-END refiner OLD vs NEW to test
  whether the DESIGN is byte-identical. Primary suspect per spec's CRITICAL SUBTLETY:
  reversed(A)-vs-forward(A) HP del-cost asymmetry + boundary-spanning deletion cost.

## CORRECTION (was reading MAIN checkout; now on worktree)
- Worktree HAS `dev/CONCAT_DP_BUILD.md` (builder's design doc) + sibling auditor files
  (dp-boundary-correctness, input-class-coverage). BUT: builder's checkpoint list is ALL
  UNCHECKED and "Checkpoints log" is EMPTY. `_USE_CONCAT_DP` flag ABSENT from worktree
  junction_scoring.py (md5 == main checkout: 05a6f25412c87f2d90eb3a23c4412d99). => concat
  DP was DESIGNED but NEVER CODED. Nothing to toggle. This is the PRIMARY integration reality.
- Builder's own doc "KNOWN RISK": penalty_table != None MAY diverge due to (1) slice-relative
  ins costs truncating a boundary-straddling HP run; (2) non-additive column-0 quirk
  `curr[0]=i*ins_costs[i-1]`. Builder predicted None=identical, table!=None=may diverge, but
  NEVER RAN THE TEST (empty log).
- `_make_arm_e.py` calls refine_bam_junctions with penalty_table_path=None => the cited "actual
  workload" is exactly the config builder predicts to be identical. To surface the divergence
  the optimization would introduce in PRODUCTION DRS (which uses a penalty table), must run
  refiner WITH a penalty table.
- Refiner uses score RAW/unbinned in candidate_tuple[0], compared with `<` (junction_refiner.py
  L664,L705/708,L709). Any float score diff can flip chosen (js,je) AND the move-margin guard
  (L738). So a scalar divergence propagates end-to-end.
- PLAN: (a) scalar harness OLD vs faithful-concat NEW with DRS penalty table -> find divergence,
  verify mechanism by hand. (b) end-to-end: monkeypatch faithful concat, run refiner WITH DRS
  table on smoke panel OLD vs NEW, diff BAMs.

## SCALAR RESULT (v1, DRS penalty table)
- DRS ins_cost(A) by hp: {1:1.25, 2:.414, 3:.218, 4:.194, 5:.179, 6:.147, 7:.184, 8:.197} (NON-monotone; dips at hp6).
- del-cost reversal claim (reverse(del_costs_rev)==forward Aref del costs): HOLDS.
- TARGETED (HP A-run straddling junction + 1 inserted A): OLD=0.1943 NEW=0.1467 (DIVERGE).
  Mechanism CLEAN: split cuts the A-run at crossing => ins_cost(hp=4)=0.1943; concat sees full
  run => ins_cost(hp=6)=0.1467. Exactly the builder's predicted slice-relative-ins divergence.
  read='GCTTGACTTGACTTAAAAAAGCTAGCT' q_split=14 intron_end=16 rescue='AAAAAAGCTAGCT'
- RANDOM (DRS table): 5412/40000 (13.5%) diverged. Both directions (some NEW>OLD by >1.0).
  => must run penalty_table=None CONTROL (v2) to confirm concat is faithful for constant-ins
  and isolate the mechanism. If None==0 divergences, table divergences are real (penalty-dependent).

## MECHANISM (why table diverges, None should not)
Given del costs are a PURE function of genome position (orientation-invariant), the ONLY
OLD-vs-concat difference is the INS cost vector:
 - OLD split: ins_costs computed per-SLICE: t1 uses hp_run_length(rescue[k:], .), t2 uses
   hp_run_length(reverse(rescue[:k]), .). A homopolymer run straddling crossing k is TRUNCATED
   in the slice.
 - concat NEW: ONE ins vector over the WHOLE rescue -> full (untruncated) run length.
 - DRS ins_cost(A) is NON-MONOTONE in hp -> truncation can move cost either way => NEW<OLD OR NEW>OLD.
 - Plus the non-additive col-0 quirk i*ins_costs[i-1] (slice vs whole ins).
For penalty_table=None: ins is constant 1.25 in BOTH; all costs are exact dyadic floats
(0.5,1.0,1.25) so no FP-reorder either => concat MUST == split if structure is faithful.
=> None-control divergences would mean MY concat has a structural bug; None==0 validates it.

## END-TO-END plan (production-like arm-C = penalty table ON)
- e2e_refine.py: subset first 400 N-op reads of mix_fair_out/aligned.sorted.bam; run
  refine_bam_junctions WITH penalty_table_path=DRS penalty_scores.tsv, n_workers=1, motif_blind=True.
  MODE=old (stock _score_junction) vs MODE=new (monkeypatch JR._score_junction=faithful concat).
  Diff the two refined BAMs (chosen js/je via N-op coords + CIGAR) read-by-read.

## *** DIVERGENCE FOUND — TWO independent mechanisms; byte-identity FAILS ***

### MECH 1 (penalty-INDEPENDENT: k=L boundary-column off-by-one) — refutes "None==identical"
None-control: 25/20000 diverged. First case debugged FULLY:
  read='GATGGCGGGG' q_split=2 intron_end=45  rescue='TGGCGGGG' (L=8), pt=None
  OLD _score_junction = 2.0  ;  concat NEW = 1.5   (NEW < OLD)
OLD split min_k[t1+t2] over k in range(L)=[0,8): best=2.0 at k=5 (rescue[5:]='GGG' matches exon2).
concat backtrace: ALL 8 rescue bases align within Aref (intron side), ending EXACTLY at column
R2=45 (=genome[ie-1], last intron base); cost = 1 sub + 1 HP-del = 1.0+0.5 = 1.5.
ROOT CAUSE: OLD loops `for k in range(L)` => suffix rescue[k:] is ALWAYS non-empty => OLD NEVER
scores "entire rescue aligns to intron end, 0 bases in exon2". The concat's free suffix over
columns [R2, n] INCLUDES column R2 (empty-exon2 / k=L), so it admits that alignment. Builder's
OWN decomposition ("LAST cell in col R2 at row k, Suffix=rescue[k:]") INCLUDES k=L (empty suffix)
but the OLD code EXCLUDES it. => concat search space strictly LARGER => concat <= split, strictly
less when whole-rescue-to-intron is cheapest. Penalty-INDEPENDENT. Builder's "penalty_table=None
should be byte-identical" is FALSE.

### MECH 2 (penalty-table-dependent: slice-relative ins costs) — irreducible for DRS
TARGETED clean case: read='GCTTGACTTGACTTAAAAAAGCTAGCT' q_split=14 ie=16 rescue='AAAAAAGCTAGCT'
  OLD=0.1943 NEW=0.1467. HP A-run straddles crossing; split truncates run -> ins_cost(hp=4)=.1943;
  concat sees full run -> ins_cost(hp=6)=.1467. DRS ins_cost(A) is NON-monotone -> both directions.
Random with DRS table: 3142/20000 (15.7%) diverge.

### PROPAGATION: score feeds RAW (unbinned) into candidate_tuple[0] (junction_refiner L705/708),
compared with `<` (L709) -> a diverging candidate score can flip chosen (js,je); also drives the
move-margin guard (L738). So both mechanisms propagate to the refined BAM.

## MECH SEPARATION CONFIRMED (v3, 20000 cases)
- None,  concat INCLUDES boundary col (k=L allowed) : 23/20000 diverged  (MECH 1)
- None,  concat EXCLUDES boundary col (k=L FIXED)   :  0/20000 diverged  => MECH 1 is EXACTLY the
  k=L/boundary-column off-by-one; FULLY isolated and FIXABLE (exclude col R2 from free suffix).
- DRS table, concat EXCLUDES boundary col (FIXED)   : 3142/20000 (15.7%) diverged  => MECH 2 (slice-
  relative ins costs) SURVIVES the k=L fix. IRREDUCIBLE for any single-DP concat under a penalty table.
CONCLUSION: byte-identity FAILS. MECH1 fixable; MECH2 inherent to the optimization approach whenever
penalty_table != None (production DRS config). Builder's own doc predicted MECH2.

## END-TO-END REFINE (production-like: penalty table ON, arm-C config)
- Confirmed TARGETED case is MECH2 (survives k=L fix): OLD=0.1943, concat-FIXED=0.1467.
- e2e uses the k=L-FIXED concat (only MECH2 remains) for clean attribution (advisor Guardrail B).
- Guardrail A: concat writes a per-pid marker file on first call to PROVE the monkeypatch fired
  in the forked pool workers (else an empty BAM diff could be a silent no-op patch = false pass).
- Subset: first 400 N-op reads of mix_fair_out/aligned.sorted.bam; refine_bam_junctions with
  penalty_table_path=DRS penalty_scores.tsv, n_workers=1, motif_blind=True. OLD then NEW, serial.
- Checkpoint: OLD refine launched.

## *** END-TO-END RESULT — refined BAM DIVERGES (production penalty-table config) ***
Subset: first 400 N-op reads, refine_bam_junctions WITH DRS penalty table, motif_blind, n_workers=1.
NEW = k=L-FIXED concat (MECH1 removed; ONLY the irreducible MECH2 ins-cost divergence remains).
- Guardrail A PASSED: patch-fired marker present (worker pid 87853 ran the concat) => NOT a silent
  no-op; the NEW codepath genuinely executed.
- Aggregate: OLD refined 309/400 (91 unchanged); NEW refined 252/400 (148 unchanged).
- Per-read: 97/400 (24.25%) have DIFFERENT chosen N-op (js,je) AND different CIGAR.
- SAM-body md5: OLD=f5ad9b51b25217e5a1b77db9c6540f0c  NEW=e881aa8e288b544936837a18edb18e32 (DIFFER).
- Examples: tid_0_r399 OLD junc=(180,275)->NEW (180,270); tid_0_r369 (180,276)->(180,270);
  tid_0_r393 (178,272)->(180,270); tid_0_r237 (180,275)->(178,272). Intron-end shifts 2-6 bp;
  N-op length + surrounding D/I ops rewritten.
CONCLUSION: Even the k=L-FIXED concat (best-case: MECH1 bug removed) produces a DIFFERENT refined
BAM on 24% of reads under the production DRS penalty-table config. MECH2 (single ins vector vs OLD's
per-slice hp_run_length) is inherent to ANY single-DP concat and propagates through the raw-score
argmin (junction_refiner L709) + move-margin guard (L738) to the emitted BAM. BYTE-IDENTITY FAILS
end-to-end. NOTE: arm-E (task-cited _make_arm_e.py) uses penalty_table=None -> would NOT surface this;
the divergence appears specifically in the penalty-table (arm-C / production DRS) config.

## FINAL VERDICT: divergence_found=TRUE, verdict_survives=FALSE.
Headline: the concat DP was NEVER BUILT (flag absent, junction_scoring.py md5==main
05a6f25412c87f2d90eb3a23c4412d99, CONCAT_DP_BUILD.md checkpoints all unchecked + empty log). All
divergences are vs a FAITHFUL RECONSTRUCTION of the design doc (validated: k=L-fixed concat is
None-byte-identical, 0/20000). Three independent reasons byte-identity cannot hold:
 (1) MECH2 ins-cost (irreducible, penalty table, builder-predicted) — scalar 15.7%, e2e 24% of reads;
 (2) MECH1 k=L boundary-column (penalty-independent, but FIXABLE) — refutes builder's None claim;
 (3) FP summation-order (non-dyadic penalty costs) can differ at 1 ULP on mathematically-equal paths.
