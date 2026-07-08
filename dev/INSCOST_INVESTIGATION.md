# ins_cost full-run vs per-cut truncation investigation

**Agent:** scoring-semantics investigator (worktree agent-a25a2c1e784ad37dc)
**Started:** 2026-07-08
**Status:** IN PROGRESS

## Question
`_score_hp_anchored` (junction_scoring.py:693) computes
`ins_costs = [ins_cost(_hp_run_length(query,i), query[i]) for i in Q]` on the
PER-K TRUNCATED substring (rescue[k:] or reverse(rescue[:k])). A homopolymer run
straddling the cut k is measured as SHORTER on each side. The DRS ins_cost table
is NON-MONOTONIC (hp6=0.1467, hp12=0.6882), so the k-sweep can SPLIT a run to get
a cheaper score. Witness: rescue='A'*12 -> old=1.7604 at k=6 vs full-run
12*ins_cost(12)=8.2584.

## Hypothesis (PI-approved to investigate)
Per-cut truncation is an IMPLEMENTATION ARTIFACT, not a model. Nanopore over-call
error is governed by the READ's FULL physical homopolymer run, so ins_cost should
be computed on the FULL rescue run (cut-independent), NOT the truncated substring.
Full-run ins_cost: (a) sounder error model; (b) removes run-splitting gaming that
can FABRICATE junctions; (c) makes ins_costs cut-INDEPENDENT -> single-pass DP
becomes valid -> unlocks ~30x speedup. It CHANGES results -> full re-validation.

## Subtlety to CHECK not assume
Interaction with the HP-DRIFT GUARD (hp_drift_margin, junction_refiner.py) which
already governs junction-slides-within-a-homopolymer.

## Plan
1. [ ] Read helpers: _hp_run_length, _precompute_del_costs, ins_cost table
2. [ ] Read junction_refiner.py hp_drift_margin guard
3. [ ] Locate mix_fair_out / mix_r3b_out refine harness + aligned bams
4. [ ] Implement _USE_FULL_RUN_INS flag (default False): precompute full-run
       ins_costs ONCE in _score_junction on the full rescue window, thread
       cut-independently into per-k _score_hp_anchored calls (like del_costs)
5. [ ] MEASURE (a) FABRICATION witness: all-A + realistic polyA over-call,
       old split-score vs full-run
6. [ ] MEASURE (b) CALL-CHANGE RATE on mix_fair_out + mix_r3b_out flag off/on
7. [ ] MEASURE (c) GUARD interaction with hp_drift_margin=3.0: R3 discovery
       (0.284) + canonical drift-fix preserved?
8. [ ] (d) Steelman per-cut truncation
9. [ ] VERDICT + SWITCH design + RE-VALIDATION plan
10. [ ] Do NOT commit, do NOT flip default

## CHECKPOINTS (append-only)
- 2026-07-08 start: plan written; read _score_junction + _score_hp_anchored.

## FINDINGS — orientation (2026-07-08)
- Real yeast DRS ins_cost table (I/AT): hp1=1.2500, hp2=0.4138, hp3=0.2184,
  hp4=0.1943, hp5=0.1785, **hp6=0.1467 (MIN)**, hp7=0.1843, hp8=0.1972,
  hp9=0.2638, hp10=0.3239, hp11=0.3472, **hp12=0.6882**. U-shaped -> splitting a
  long run into ~hp6 chunks minimizes per-base ins cost. CONFIRMS gaming.
- ins_costs computed on the READ (query) hp run; del_costs on the GENOME hp run.
  So ins=read-context, del=genome-context. Hypothesis (over-call governed by
  read's full physical run) aligns with using full READ run. SOUND direction.
- GAMING FIRES ONLY when penalty_table is not None (line 766). With
  penalty_table=None, ins_costs=[ins]*Q (flat, already cut-independent) -> flag
  is a no-op. IMPLICATION: _make_arm_ff.py / _make_r3b_arms.py use
  penalty_table_path=None -> the ins-cost flag has ZERO effect on arm_Ff/arm_E.
  To exercise ins-cost, must refine with penalty_table (arm-C style:
  motif_blind + penalty_table_path=DEFAULT_PENALTY_TABLE).
- Panels: mix_fair_out=5800 reads (R0flank 5600 = canonical flank-HP drift test),
  mix_r3b_out=5100 reads (R3=1600 cryptic discovery, R0/R1 canonical, WT).
  R3 discovery 0.284 = recovery on R3 rung; drift-fix = R0flank recovery.
- Indexing for cut-independent threading (analogous to query slicing):
  t1 query=rescue[k:] -> ins costs = full_ins_costs[k:]
  t2 query=rescue[:k][::-1] -> ins costs = full_ins_costs[:k][::-1]
- CHECKPOINT: orientation done; implementing _USE_FULL_RUN_INS flag next.

## ANCHOR PINNED: 0.284 (2026-07-08)
- Source: dev/HANDOFF_DEL_OPEN_STAIRCASE.md: "R3 recovery: B=0.608 E=0.608 ...
  (HP cell all 0.284)". Also HANDOFF_NATIVE_ALIGNER_VETTING.md: arm-B 0.284 ==
  arm-E 0.284 at every margin; guard touches 0 reads on R3-HP.
- => 0.284 is a **penalty_table=None** number (arm-B/E). My flag lives on the
  penalty_table path -> STRUCTURAL NO-OP on the arms that produced 0.284. So
  full-run "preserves" 0.284 trivially. The MEANINGFUL (c) test is a NEW arm:
  motif_blind + penalty_table + hp_drift_margin=3.0, flag off (new baseline) vs on.
- Advisor guidance adopted: (1) gate flag by env var RECTIFY_FULL_RUN_INS=1, one
  arm per subprocess (fresh forks + M1 memory cap); (2) byte-identity via test
  suite flag-off, NOT diff vs stale arm_C.bam; (3) build penalty_table arms for
  (b)/(c); (4) direction = classify changed calls by _canonical_tier, off<->on
  diff (not vs arm-B); (5) don't build concat DP, just argue validity.
- CHECKPOINT: anchor pinned; implementing flag.

## (a) FABRICATION WITNESS — DONE (2026-07-08)
Implementation validated: reproduces the stated witness numbers EXACTLY.
Script: scratchpad/witness.py (direct _score_junction calls, controlled genome).
- W1 all-A(12), A-free genome (pure gaming demo):
    flag OFF (per-cut split @k=6): **1.7604**  (= 6*ic6 + 6*ic6, ic6=0.1467)
    flag ON  (full-run):          **8.2584**  (= 12*ic12, ic12=0.6882)
    => run-splitting discount REMOVED; 4.69x more expensive. Full-run charges the
       12-A insertion at its true physical run length, not two gamed hp6 halves.
- W2 realistic poly-A over-call (genome exon2 starts 8-A; read rescue 12-A):
    flag OFF: 0.7772   flag ON: 2.7528   delta +1.9756
    => the 4 over-called bases are charged at the full run length under full-run.
- W3 control non-HP rescue (GTAC...): OFF=ON=0.0000 -> flag only touches HP runs.
CHECKPOINT: (a) done, numbers persisted.

## BYTE-IDENTITY (flag OFF) — CONFIRMED (2026-07-08)
`pytest -m "not slow"` on junction refiner/scoring/validator/false-junction +
test_validation_reads: **201 passed, 25 skipped** (167s). Flag OFF takes the exact
legacy branch (precomputed_ins_costs=None -> old per-substring ins_costs). Green.
CHECKPOINT: byte-identity verified.

## CONCAT-DP UNLOCK ANALYSIS (from CONCAT_DP_AUDIT_end-to-end-integration.md)
The prior audit proved single-pass concat DP canNOT be byte-identical to the
CURRENT (per-cut) split for THREE reasons:
- MECH2 = per-SLICE ins costs (hp_run_length on the truncated substring). Called
  "IRREDUCIBLE for ANY single-DP concat under a penalty table." scalar 15.7%
  divergence, e2e 24% of reads. **This is EXACTLY the artifact _USE_FULL_RUN_INS
  removes.**
- MECH1 = k=L boundary column off-by-one (penalty-INDEPENDENT). Audit showed it is
  "FULLY isolated and FIXABLE (exclude col from free suffix)"; k=L-fixed concat is
  None-byte-identical 0/20000.
- MECH3 = FP summation order, 1 ULP on non-dyadic costs (tolerance issue).

KEY UNLOCK ARGUMENT (durable):
- OLD split assigns ins cost by SLICE-relative position -> concat (one ins vector
  over whole rescue) cannot match -> MECH2 irreducible.
- FULL-RUN split assigns ins cost by ABSOLUTE rescue position (cut-independent
  vector full_ins_costs[j]) -> this is EXACTLY what a single-pass concat DP
  naturally computes. The audit's own control (penalty_table=None == constant ins
  = trivially cut-independent) is concat==split byte-identical (0/20000). Full-run
  generalizes that property to the penalty-table case: any per-position ins vector
  indexed by ABSOLUTE position makes concat==split. => MECH2 DISSOLVES.
- Therefore making full-run DEFAULT removes the exact property the audit cited as
  irreducible. Single-pass DP becomes VALID (still needs MECH1 boundary fix, which
  audit already proved isolatable, + MECH3 FP tolerance). ~30x speedup unlockable.
- SCOPE: I did NOT build the concat DP (out of scope). This is the validity
  argument only.

## (d) STEELMAN per-cut truncation
- Principled defense FAILS: ins_cost fires ONLY on OVER-CALLED bases with no
  reference A to match (true poly-A that matches genome A's are MATCHES, governed
  by del/match, not ins). An over-call is a property of the READ's ONE continuous
  physical homopolymer the pore traversed -> the full read run is unambiguously the
  right context; a k-cut splitting it has no physical meaning. There is NO scenario
  where truncating the run at an arbitrary alignment cut point is the correct error
  model. Full-run DOMINATES per-cut (and even full-run still truncates at the
  30bp rescue-window cap -> a maximally-correct model would use the full READ run).
- Only HONEST (practical, not principled) defense: the ins_cost TABLE and the
  downstream margins (hold_margin, hp_drift_margin) were calibrated/tuned WITH the
  per-cut score distribution in play. Full-run raises HP-insertion scores (W1:
  1.76->8.26), shifting the score scale -> margins may need re-tuning. That is the
  real switching cost, and the reason full re-validation (not just a flip) is
  required. It does NOT make per-cut the sounder MODEL.
CHECKPOINT: concat unlock + steelman persisted.

## REFINE TIMING + BATCH (2026-07-08)
- Single refine mix_r3b_out (arm-C style, penalty_table, n_workers=4): 142.4s;
  flag OFF refined 726 of 4800 N-op reads. ins_off_ng bam = 1338911 bytes.
- Batching remaining 7 arms sequentially (fresh subprocess each = clean env flag):
  mix_r3b_out {ins_on_ng, ins_off_g3, ins_on_g3}; mix_fair_out {ins_off_ng,
  ins_on_ng, ins_off_g3, ins_on_g3}. ~17 min. Driver: scratchpad/drive_ins.sh.
CHECKPOINT: batch launched.

## (b) CALL-CHANGE — mix_r3b_out (penalty_table, no guard), flag OFF vs ON
- ins_off_ng refined 726/4800; ins_on_ng refined 725/4800 (135s).
- CALL-CHANGE: **10/4800 = 0.21%** (small; this panel is cryptic-discovery, not
  poly-A-fabrication heavy).
- Direction of 10 changed calls (recovery x canonicality):
    3  recov_GAINED  x became_MORE_canonical (tier 4->0: full-run FIXED a wrong
       non-canonical, recovered truth — the anti-fabrication win)
    2  recov_GAINED  x tier_same
    4  recov_LOST    x became_LESS_canonical (tier 0->4: full-run moved a
       true/canonical call to non-canonical — a cost)
    1  recov_neutral x became_LESS_canonical
  => roughly balanced on THIS panel; NET recovery ~0.
- Per-rung recovery OFF->ON: R3 0.589->0.589 (+0.001, DISCOVERY PRESERVED),
  R0 0.979->0.978, R1 0.866->0.867, WT 0.774->0.774. (NB penalty-table R3 baseline
  is 0.589, NOT the 0.284 penalty_table=None HP-cell number.)
CHECKPOINT: (b) r3b done.

## (c) GUARD INTERACTION — mix_r3b_out (hp_drift_margin=3.0)
Guard cuts refinements 726->427 (off) / 725->429 (on) — guard vetoes into-HP moves.
- WITH guard, flag OFF vs ON: CALL-CHANGE **5/4800 = 0.10%** (guard SUPPRESSES the
  small full-run-induced changes, 0.21%->0.10%). R3 0.604->0.604 (DISCOVERY
  PRESERVED), R0/R1/WT within +-0.001. NO conflict.
- Full-run ON, guard OFF vs guard=3.0: CALL-CHANGE 296/4800 = **6.17%** -> the guard
  STILL FIRES robustly under full-run, and IMPROVES recovery: WT +0.166, R1 +0.048,
  R3 +0.014, R0 +0.000. => the guard's drift-fix mechanism is fully intact when
  full-run is active (guard acts on the MOVE decision / into-HP margin, orthogonal
  to the ins-cost SCORE). 
- ANSWER: full-run and the guard are COMPATIBLE and complementary. The 0.284 anchor
  (penalty_table=None) is untouched by the flag (structural no-op) -> trivially
  preserved; on the penalty_table config where the flag acts, R3 discovery
  (0.589 no-guard / 0.604 guard) is preserved to +-0.001 under the flag flip, both
  with and without the guard. NO conflict with the HP-drift guard.
CHECKPOINT: (c) r3b guard interaction done.

## (b) CALL-CHANGE — mix_fair_out (R0flank HP-drift), penalty_table, no guard
- ins_off_ng refined 1085/5600; ins_on_ng refined 1056/5600.
- CALL-CHANGE: **42/5600 = 0.75%**.
- DIRECTION (the anti-fabrication signal):
    **32 recov_GAINED x became_MORE_canonical** (tier 4->0)
     5 recov_LOST   x became_LESS_canonical (tier 0->4)
     3 recov_neutral x became_MORE_canonical
     2 recov_neutral x became_LESS_canonical
  Tier shift: **35 reads 4->0 (non-canonical -> canonical)** vs only 6 0->4, 1 4->6.
- R0flank recovery 0.777 -> 0.781 (**+0.005**, net positive).
=> On the HP-drift panel full-run is DIRECTIONAL and net-positive: it removes
   fabricated non-canonical calls (poly-A over-call was being cheaply absorbed by
   run-splitting at a wrong junction; full-run makes that expensive so the true
   canonical junction wins). FEWER fabricated non-canonical junctions. CONFIRMS the
   hypothesis' fabrication-removal claim on a realistic panel.
- Combined (b): r3b 0.21% (balanced), fair 0.75% (net anti-fabrication 35:6).
CHECKPOINT: (b) fair done — directional anti-fabrication confirmed.

## (c) GUARD INTERACTION — mix_fair_out (R0flank drift, hp_drift_margin=3.0)
- ins_off_g3 refined 234/5600; ins_on_g3 refined 238/5600 (guard cuts moves hard).
- WITH guard, flag OFF vs ON: CALL-CHANGE **4/5600 = 0.07%**; R0flank 0.899->0.898
  (-0.001). All 4 changes 0->4 (tiny recov_LOST). Under the guard the flag is
  essentially inert on drift.
- Full-run ON, guard OFF vs guard=3.0: CALL-CHANGE 818/5600 = **14.61%**; R0flank
  0.781 -> 0.898 (**+0.117**). Guard's drift-fix is FULLY INTACT under full-run.
- SYNTHESIS: guard and full-run target OVERLAPPING failure modes (both discourage
  HP-insertion junction moves). Guard is the DOMINANT drift-fix (+0.117); full-run
  alone is a weaker complementary version (+0.005 on the drift panel). guard+full-run
  (0.898) ~= guard+per-cut (0.899). NO conflict; full-run's marginal drift benefit
  is largely subsumed by the guard when both are on. This is WHY flag-off-vs-on
  call-change collapses to 0.07% under the guard.
CHECKPOINT: (c) fair guard interaction done. ALL MEASUREMENTS COMPLETE.

## SUMMARY TABLE (penalty_table = yeast DRS, motif_blind)
| panel | guard | refined off/on | call-change off->on | key recovery off->on |
|-------|-------|----------------|---------------------|----------------------|
| r3b   | none  | 726/725        | 10/4800 = 0.21%     | R3 0.589->0.589      |
| r3b   | 3.0   | 427/429        | 5/4800 = 0.10%      | R3 0.604->0.604      |
| fair  | none  | 1085/1056      | 42/5600 = 0.75%     | R0flank 0.777->0.781 (35 reads 4->0) |
| fair  | 3.0   | 234/238        | 4/5600 = 0.07%      | R0flank 0.899->0.898 |
Guard-fires-under-full-run: r3b 6.17% (WT +0.166), fair 14.61% (R0flank +0.117).

## FLAG-ON TEST SUITE (default-on readiness) — 2026-07-08
- "46 tests" reconciled: = tests/test_junction_refiner.py + tests/test_hp_drift_guard.py
  (the narrow suite in HANDOFF_DEL_OPEN_STAIRCASE.md / DEL_OPEN_DELTA_FINDING.md).
- 46-suite flag OFF: 46 passed, 17 skipped. flag ON: **46 passed, 17 skipped**. CLEAN.
- Larger subset (refiner/scoring/validator/false-junction + validation_reads),
  flag ON: [pending bs4rk8jvh].

## FLAG-ON TEST SUITE (cont.)
- Larger subset flag ON: **201 passed, 25 skipped** (222s) — IDENTICAL to flag OFF
  (201/25). => default-on is TEST-CLEAN on the junction/scoring/validation subset:
  none of the fixture reads fall in the 0.21-0.75% of production calls that change.
- CAVEAT: numba is UNAVAILABLE in this env (miniconda base) -> the flag's numba DP
  path (precomputed reversed-slice list -> np.array for t2) is structurally sound
  but UNEXERCISED on the cluster's numba-ON deployment path. FULL `not slow` suite
  (1603 tests) not yet run flag-ON.

## VERDICT (2026-07-08)
FULL-RUN IS THE SOUNDER MODEL. Primary case (NOT recovery, which is small/redundant):
1. PRINCIPLED ERROR MODEL: a Nanopore over-call (insertion with no ref base to
   match) is a property of the READ's ONE continuous physical homopolymer the pore
   traversed. Per-cut truncation charges that run at an arbitrary alignment cut k
   with no biological meaning; it only shrinks the cost via the U-shaped table
   minimum (hp6). The steelman shows per-cut is NEVER *more* correct than full-run.
2. DISSOLVES MECH2: full-run assigns ins cost by ABSOLUTE rescue position ->
   cut-independent -> exactly what a single-pass concat DP computes. It removes the
   precise property the prior audit called "irreducible," unlocking the ~30x DP.
CORROBORATING (not primary): on the HP-drift panel, full-run moves 35 reads
   non-canonical->canonical vs 6 the other way (anti-fabrication behaves as
   predicted); recovery +0.005 fair / ~0 r3b (small, largely redundant with guard).
COUNTER-SIGNAL (honest): a small consistent 0->4 population (4-6 reads/panel) where
   full-run moves a canonical call to non-canonical & loses recovery — mechanism:
   raising the poly-A insertion cost at the TRUE junction can let an alternative
   that absorbs the run as matches/dels win. Bounded, understood, net-positive.

## SWITCH DESIGN (make full-run default; do NOT do it here)
- Flip: `_USE_FULL_RUN_INS: bool = os.environ.get("RECTIFY_FULL_RUN_INS", "1") != "0"`
  (default True; env can force-off for A/B). KEEP the `precomputed_ins_costs is None`
  fallback branch in _score_hp_anchored ALIVE — the flat/no-table path and any
  non-precomputing caller still need it; only _score_junction precomputes.
- DP UNLOCK: with full-run default, OLD reference behavior == full-run == the single
  ins vector a concat computes. The single-pass DP then computes the SAME score (not
  free byte-identity): still needs (i) MECH1 boundary-column fix (exclude the k=L
  free-suffix column; audit proved isolatable, None-byte-identical 0/20000) and
  (ii) an FP summation tolerance for non-dyadic costs (MECH3, ~1 ULP). MECH2 gone.

## RE-VALIDATION PLAN (all 5, with expected direction)
1. FULL test suite flag-ON (`pytest -m "not slow"`, 1603): expect green (subset
   already 247/247 clean); this is the default-on gate.
2. Yeast DRS real-data HP-drift transfer (task #17, Sherlock): expect mirror of sim
   — drift-fix preserved, small anti-fabrication, guard still dominant.
3. Human ONT DRS transfer (task #18, GENCODE truth): UNEXERCISED direction — human
   HP-length distribution differs; must run (this is where the balance could shift).
4. del_open verdict (DEL_OPEN_DELTA_FINDING / arm-F): reached under PER-CUT ins;
   full-run shifts the score scale (W1 1.76->8.26) so the arm-F/del_open comparison
   must be RE-RUN, not assumed to carry.
5. Guard-margin re-tune (task #16): full-run does part of the drift-suppression
   itself -> predict the smallest hp_drift_margin with zero discovery cost may
   DECREASE from 3.0. Testable; re-sweep.
Plus: re-run cat1-cat9 validation fixtures flag-ON (part of #1) and re-confirm
del_open verdict on the same panels.

## LEAST-SURE RISK
The numba-ON DP path is unverified on the deployment target: all local validation
ran numba-OFF (pure-Python DP). The flag's numba branch feeds a precomputed list
(incl. the reversed slice full_ins_costs[:k][::-1]) into np.array->_score_hp_dp_numba;
structurally sound but UNEXERCISED where it actually ships. Must run the cluster
numba build flag-ON before trusting the numbers there.
CHECKPOINT: verdict + switch + revalidation + risk persisted. INVESTIGATION COMPLETE.
