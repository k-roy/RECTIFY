# INSCOST reference-column ins_cost — build + 3-way A/B

## Task
Implement REFERENCE-COLUMN ins_cost behind a new flag `_USE_REFCOL_INS` (default False),
keeping `_USE_FULL_RUN_INS` selectable. Ref-column = charge insertion cost by the HP-run
length of the GENOME/REFERENCE context at the aligned position (mirror _precompute_del_costs:
genome_seq, ref_genome_start, ref_genome_rev). Confirm cut-independence. Then 3-way A/B on
mix_fair_out + mix_r3b_out: per-cut (both flags off) vs full-run vs ref-column.

KEY QUESTION: does ref-column keep full-run's ~35 anti-fab wins while ELIMINATING the
4-6/panel canonical-demotion (0->4) losses?

## Plan
1. Read junction_scoring.py: _score_junction, _score_hp_anchored, _precompute_del_costs,
   _USE_FULL_RUN_INS prototype. Understand how ins_costs enters the k-sweep DP.
2. Read dev/INSCOST_INVESTIGATION.md, INSCOST_SYNTHESIS.md, INSCOST_AUDIT_model-correctness.md.
3. Design ref-column: precompute ins_cost per absolute genome position from genome HP context.
4. Implement behind _USE_REFCOL_INS (default False). Verify cut-independence by construction.
5. Locate the A/B harness (_make_arm_e.py / refine_bam_junctions pattern), panels
   mix_fair_out / mix_r3b_out. Reuse aligned bams. Lean, one arm at a time (M1 kills heavy refines).
6. Run 3-way: per-cut / full-run / ref-column on both panels.
7. Compute: per-junction call-change %, tier 4->0 (anti-fab wins) vs 0->4 (canonical-demotion
   losses), recovery delta, R3 discovery.
8. A/B with hp_drift_margin=3.0 (guard interaction).
9. Recommend variant with evidence. NO commit, NO default flip.

## Checkpoints
- [init] plan written.

## DESIGN CONFIRMED (2026-07-08, advisor-folded)
- Ref-column = REF-COLUMN-INDEXED DP, not a query-indexed vector swap. Insertion aligns
  to NO ref base -> its cost is a property of the DP COLUMN (gap) where `above` fires.
- Calibration axis VERIFIED in scripts/calibration/empirical_cigar_error_profiler.py
  Phase 5 (lines 1319-1324): hp_left=hp_arr[pos-1], hp_right=hp_arr[pos],
  hp_len=max(hp_left,hp_right); hp_base = ref_seq[pos-1] if hp_left>=hp_right else ref_seq[pos].
  => insertion cost is a PER-GAP quantity: length R+1 vector (gaps 0..R), NOT length R.
    gap 0: only right neighbor (ref[0]); gap j (1..R-1): max(run(gp_{j-1}), run(gp_j));
    gap R: only left neighbor (ref[R-1]).
- CUT-INDEPENDENCE (stronger than full-run): the ins_col vector is a function of the FIXED
  ref window + absolute genome position ONLY -> computed ONCE outside the k-loop (like
  del_costs_fwd/rev), passed UNCHANGED to every k. NO [k:]/[:k][::-1] slicing. That IS the proof.
- DP wiring: curr[0] = i * ins_col[0]; above = prev[j] + ins_col[j] for j in 1..R.
- numba TRAP: existing _score_hp_dp_numba takes ic_arr length Q (query-indexed). Ref-col
  ic is length R+1 -> WRONG-indexed if fed to that kernel. FORCE pure-Python branch when refcol.
- Both-flags-on: refcol WINS (mutually exclusive precedence).
- WITNESS target (audit W-model, 8-A genome acceptor, read over-calls to 12-A): expect
  refcol ~= 0.7888 (= 4*ins_cost(8)); full-run 2.7528; per-cut 0.7772.
- CHECKPOINT: design confirmed, calibration axis verified. Implementing next.

## CODE WRITTEN (2026-07-08)
- hp_penalty.py: added _precompute_refcol_ins_costs (length R+1, per-gap,
  max(run_left,run_right) genome context, penalty_table.ins_cost lookup; flat fallback).
- junction_scoring.py:
  * flag _USE_REFCOL_INS (env RECTIFY_REFCOL_INS=1, default False), precedence over full-run.
  * import _precompute_refcol_ins_costs.
  * _score_hp_anchored: new param precomputed_refcol_ins (length R+1, column-indexed);
    R==0 returns Q*ins_col[0]; dedicated pure-Python DP branch (above=prev[j]+ins_col[j],
    curr[0]=i*ins_col[0]); FORCES pure-Python (skips numba Q-indexed kernel).
  * _score_junction: refcol_ins_fwd/rev precomputed ONCE outside k-loop, passed UNCHANGED
    to t1/t2 (no slicing) -> cut-independence by construction.
- CHECKPOINT: code written. Running witness next.

## WITNESS (3-way) — 2026-07-08
ins_cost: ic4=0.1943 ic6=0.1467 ic8=0.1972 ic12=0.6882
- Wm (audit W-model: genome exon2 starts 8-A ref run; read over-calls to 12-A):
    per-cut = 0.7772  |  full-run = 2.7528  |  REF-COL = 0.7888
    => ref-col == 4*ins_cost(8) = 0.7888 EXACTLY -> lands on the calibration-correct
       axis (audit's predicted number). CONFIRMS implementation is on the reference
       axis, NOT read-run (2.7528) and NOT per-cut-collapse (0.7772).
- W1 (all-A(12) rescue, A-free genome flanks):
    per-cut = 1.7604  |  full-run = 8.2584  |  REF-COL = 6.8748  [investigating W1]
- CHECKPOINT: witness Wm confirms ref axis. Verifying W1 mechanism next.

## W1 mechanism CONFIRMED (not a bug)
- W1 A-free genome: every insertion gap has genome run=1 -> ins_cost(1,base)=1.25.
  Ref-col DP optimum 6.8748 charges the 12 pure-over-call A's against the genome's
  (absent) A context. This is CORRECT anti-fabrication: an over-call with NO genome
  HP to explain it is EXPENSIVE (vs per-cut's 1.7604 run-splitting discount). W1 is an
  artificial pure-gaming control; the DECISIVE realistic case is Wm (0.7888, exact).

## BATCH PLAN (12 arms, one subprocess at a time — M1)
Per panel {mix_fair_out, mix_r3b_out} x variant {percut, fullrun, refcol} x guard {0.0, 3.0}
Labels: <variant>_ng (no guard), <variant>_g3 (guard 3.0).
Builder: _make_ins_arms_refcol.py (honors RECTIFY_FULL_RUN_INS / RECTIFY_REFCOL_INS).
Adjudicator: _ins_compare.py <wd> arm_off arm_on --rungs ...
Emits arm_mbF_<label>.bam. ~140s each -> ~28min. Driver: scratchpad/drive_refcol.sh (logs per arm).
- CHECKPOINT: inputs verified, batch launching.

## SMOKE + BYTE-IDENTITY (2026-07-08)
- Smoke arm refcol_ng on r3b: flag captured (_USE_REFCOL_INS=True), refined 1134/4800
  in 140.6s. Harness OK.
- Byte-identity (flag OFF default): tests/test_junction_refiner.py + test_hp_drift_guard.py
  => 46 passed, 17 skipped. Default path unchanged by the code additions.
- CHECKPOINT: launching full 12-arm batch (refcol_ng r3b already done -> 11 to go).

## RESUME (if this session dies mid-batch)
- Batch driver: scratchpad/drive_refcol.sh (idempotent: skips arms whose bam exists).
  Background id b17rh1hrv; per-arm logs scratchpad/refcol_logs/<panel>_<label>.log;
  driver stdout scratchpad/drive_refcol.out.
- To resume: re-run `bash scratchpad/drive_refcol.sh` (skips completed arms).
- Arms produced per panel: arm_mbF_{percut,fullrun,refcol}_{ng,g3}.bam (12 total).
- After ALL arms exist, adjudicate with _ins_compare.py:
    cd scripts/benchmark/noncanon_sim
    PY=/Users/kevinroy/miniconda3/bin/python
    # KEY 3-way no-guard: per-cut vs full-run vs ref-col
    $PY _ins_compare.py mix_fair_out mbF_percut_ng mbF_fullrun_ng --rungs R0flank,INTRONFREE
    $PY _ins_compare.py mix_fair_out mbF_percut_ng mbF_refcol_ng  --rungs R0flank,INTRONFREE
    $PY _ins_compare.py mix_fair_out mbF_fullrun_ng mbF_refcol_ng --rungs R0flank,INTRONFREE
    $PY _ins_compare.py mix_r3b_out  mbF_percut_ng mbF_fullrun_ng --rungs R3,R0,R1,WT
    $PY _ins_compare.py mix_r3b_out  mbF_percut_ng mbF_refcol_ng  --rungs R3,R0,R1,WT
    $PY _ins_compare.py mix_r3b_out  mbF_fullrun_ng mbF_refcol_ng --rungs R3,R0,R1,WT
    # guard interaction (hp_drift_margin=3.0): repeat with _g3 labels
  NOTE: _ins_compare takes bare labels; arms are named arm_mbF_<label>.bam ->
  pass off/on as "mbF_percut_ng" etc (it prepends arm_).
- KEY QUESTION: does refcol keep full-run's ~35 (fair) 4->0 anti-fab wins while
  DROPPING the 4-6/panel 0->4 canonical-demotion losses? Compare percut->refcol
  tier_shift vs percut->fullrun tier_shift.

## DP CORRECTNESS CHECK (2026-07-08)
- Rolling-row refcol DP vs a brute-force full-2D reference DP (column-indexed ins):
  small controlled case (genome 8-A run, query over-calls to 12-A) -> BOTH = 0.7888,
  match to 1e-9. len(ic)=R+1 confirmed. Column-indexed insertion transition correct.
- CUT-INDEPENDENCE: refcol_ins_fwd/rev are precomputed ONCE (functions of the fixed
  ref window + absolute genome pos only) and passed UNCHANGED to every k in the
  k-sweep (no [k:]/[:k][::-1] slicing) -> ins cost is a function of absolute genome
  position, NOT the cut k. Cut-independent BY CONSTRUCTION -> unlocks single-pass DP.
- CHECKPOINT: DP verified. Batch running (b17rh1hrv); awaiting arms for A/B.

## BATCH RELAUNCH (2026-07-08) — nohup-detached
- First run_in_background driver died at shell-state reset (only pre-existing arms present).
- Relaunched DETACHED: nohup bash scratchpad/drive_refcol.sh (pid in scratchpad/.driver_pid).
  Monitored via kill -0 <pid> + arm-count. Prior investigation's arm_mbF_{off,on}_ng exist
  (=per-cut/full-run) but I rebuild percut/fullrun/refcol FRESH so all 3 come from identical
  current code. Adjudication ready: scratchpad/run_ab.sh -> scratchpad/ab_results.txt.
- CHECKPOINT: detached driver alive, refining fair/percut_ng.

## CODE FINAL-REVIEW (2026-07-08)
- t1/t2 call sites thread precomputed_refcol_ins=refcol_ins_fwd/rev (unchanged per k).
- Precedence: when refcol active, full_ins_costs stays None (elif) -> _score_hp_anchored
  refcol branch (checked first) fires, query-indexed path bypassed. Both-flags-on ->
  refcol wins. Byte-identity OFF confirmed (46 tests). Ready for A/B once arms land.

## DUPLICATE-DRIVER CLEANUP (2026-07-08)
- Found TWO drivers racing (orig run_in_background b17rh1hrv pid 51141 + nohup 51448),
  both writing drive_refcol.out (interleaved) and both refining -> 20 workers = M1 risk.
- Killed the run_in_background one (explicit pid 51141 -> SIGTERM; its zsh wrapper exited).
  Kept nohup 51448. Now ONE driver, ONE arm's workers (10). Idempotent SKIP means no
  corruption (arms build to a temp then sort_and_index; no arm_mbF_refcol_ng.bam on fair yet).
- Waiter b7wc3c1n5 tracks pid 51448 (still valid). Progress so far: fair percut_ng(1090),
  fullrun_ng(1061) done; refcol_ng refining.
- CHECKPOINT: single driver restored; awaiting 12-arm completion.

## FAIR no-guard REFINE COUNTS (early signal, 2026-07-08)
- fair/percut_ng  refined 1090/5600
- fair/fullrun_ng refined 1061/5600
- fair/refcol_ng  refined  953/5600
  => refcol makes FEWER moves than per-cut (953 vs 1090) AND fewer than full-run (1061):
     consistent with refcol charging poly-A over-calls at genome context (not gaming
     splits, not read-run inflation) -> fewer spurious into-HP junction moves.
- CHECKPOINT: fair no-guard arms done; awaiting remaining 7 arms (fair g3 x3, r3b x6-1).

## STATUS (mid-batch, 2026-07-08)
- 4/12 arms done: fair {percut_ng 1090, fullrun_ng 1061, refcol_ng 953}, r3b {refcol_ng 1134}.
- Now: fair percut_g3 refining (slowed by 2 lingering orphan refcol_ng workers from the
  killed dup driver stealing CPU; they'll drain). Driver 51448 alive. Waiter b7wc3c1n5 armed.
- Remaining: fair {percut_g3, fullrun_g3, refcol_g3}, r3b {percut_ng, fullrun_ng, percut_g3,
  fullrun_g3, refcol_g3}. Resume = re-run scratchpad/drive_refcol.sh (idempotent).

## FAIR guard=3.0 REFINE COUNTS (2026-07-08)
- fair/percut_g3  refined 236/5600
- fair/fullrun_g3 refined 240/5600  (refcol_g3 pending)
  => guard heavily suppresses moves (236-240 vs ~1000 no-guard), matching prior
     investigation (234/238). FAIR panel 6/6 done. Driver -> r3b (6 arms, refcol_ng SKIP).

## STATUS (2026-07-08, mid-r3b): FAIR DONE 6/6; M1 SWAP-BOUND
- FAIR (5600) all 6 arms built: percut_ng 1090, fullrun_ng 1061, refcol_ng 953,
  percut_g3 236, fullrun_g3 240, refcol_g3 (in progress). r3b refcol_ng 1134 (smoke).
- M1 swap 13.5/14.3 GB used -> arms slow (~3-6 min each) but progressing (380% CPU, not stalled).
- Waiter b7wc3c1n5 armed on driver pid 51448. Remaining: refcol_g3 fair + r3b {percut_ng,
  fullrun_ng, percut_g3, fullrun_g3, refcol_g3}.
- ADJUDICATION READY: scratchpad/run_ab.sh -> scratchpad/ab_results.txt (all 3-way + guard).
  Parse tier_shift 4->0 vs 0->4 for percut->fullrun vs percut->refcol; R3/R0flank recovery.

## FAIR COMPLETE 6/6 (2026-07-08): refcol_g3 322/5600. Driver -> r3b panel.

## *** FAIR no-guard 3-way A/B — DECISIVE, SURPRISING (2026-07-08) ***
percut->fullrun (5600 reads):
  CALL-CHANGE 42/5600 = 0.75%
  tier shift: 4->0 = 35 (anti-fab wins) | 0->4 = 6 (canon-demotion losses) | 4->4 = 1
  R0flank recovery 0.778 -> 0.783 (+0.005)
  => matches prior investigation EXACTLY.

percut->refcol (5600 reads):
  CALL-CHANGE 694/5600 = 12.39%  (16x MORE churn than full-run!)
  tier shift: 4->0 = 358 | 0->2 = 20 | 0->4 = 196 | 2->0 = 2 | 4->4 = 118
  R0flank recovery 0.778 -> 0.789 (+0.011)  (higher than full-run's +0.005)

*** VERDICT ON THE AUDIT HYPOTHESIS: FALSIFIED. ***
The audit (INSCOST_AUDIT_model-correctness KEY FINDING #3) predicted ref-column would
KEEP full-run's ~35 anti-fab wins while ELIMINATING the 4-6/panel 0->4 canonical-demotion
losses. INSTEAD ref-column produces a MUCH LARGER change: 358 wins (4->0) BUT ALSO 196
losses (0->4) -- the demotion class BALLOONS (196 vs full-run's 6), it does NOT vanish.
Net recovery IS higher (+0.011 vs +0.005) but via massively more churn, both directions.
Ref-column does NOT cleanly dominate full-run; it is a far more AGGRESSIVE re-scorer.
- CHECKPOINT: fair no-guard A/B done, persisted. Investigating mechanism next.
