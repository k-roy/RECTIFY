# INSCOST AUDIT — LENS: MODEL CORRECTNESS

Adversarial auditor. Claim under audit: FULL-RUN ins_cost is the sounder Nanopore
error model (vs per-cut truncation) AND is safe to switch (unlocks single-pass DP).

My lens: **MODEL CORRECTNESS**. Is full-run ins_cost REALLY the correct Nanopore
error model, or does it OVER-penalize genuine cases? Steelman per-cut HARD: find a
real biology case where full-run's single-large-run cost is WRONG and the split
cost is right. Does full-run break the legitimate junction-slide-within-a-homopolymer
that the hp-drift guard handles?

## Plan
1. Read dev/INSCOST_INVESTIGATION.md (durable record of the investigation).
2. Read junction_scoring.py: _score_junction (807), _score_hp_anchored (693),
   ins_cost, _hp_run_length. Understand exactly what is truncated and how full-run
   changes it.
3. Read junction_refiner.py hp-drift guard (hp_drift_margin) — the ONE subtlety.
4. Steelman per-cut: construct biology cases:
   (a) splice truly inside a homopolymer (junction slide within HP run).
   (b) a read run that is genuinely TWO separate genomic runs joined only by the read
       (exon2 starts with A-run, intron-end also A-run; read run straddles).
   (c) over-penalization: does full-run charge ins for bases that are actually
       genomic matches (poly-A that IS in the reference)?
5. Toggle flag in-process, compute scores on the constructed cases, PERSIST numbers.
6. Verdict: certify full-run sound, OR exhibit the case where per-cut is right.

## Checkpoints
- (init) file created, plan written.
- Read investigation record + junction_scoring.py (_score_junction 839, _score_hp_anchored 710),
  hp_penalty.py (_hp_run_length 87, ins_cost 293), junction_refiner.py guard (730-741),
  calibrator empirical_cigar_error_profiler.py.

## KEY FINDING #1 — full-run indexes ins_cost on the WRONG AXIS (calibration mismatch)
The penalty table ins_cost(hp, base) is CALIBRATED by the REFERENCE homopolymer run
length, not the read run length:
- profiler docstring line 15: "Look up the homopolymer run length at that REFERENCE position."
- Phase-5 insertion tally (lines 1320-1322): hp_left=hp_arr[pos-1], hp_right=hp_arr[pos],
  hp_len = max(hp_left, hp_right); hp_arr = _precompute_hp_array(ref_seq) = REFERENCE run.
  hp_base is a REFERENCE base (ref_seq[...]). So ins_cost(h) answers: "insertion adjacent
  to a REFERENCE h-run -> penalty."
BUT the scorer indexes by the READ run: _score_hp_anchored line 800 uses
`_hp_run_length(query, i)`, and _score_junction line 940 (full-run) uses
`_hp_run_length(rescue, j)` over the WHOLE rescue window.
- DEL side is calibration-consistent (indexed by genome_seq/reference run).
- INS side is NOT: it uses the read run. Full-run AMPLIFIES this mismatch by using the
  LARGEST read run (over the full window), which INCLUDES the over-called bases.
Physical truth: an over-call in a reference h-run makes read_run = h + (#overcalls) > h.
Calibrated index = h. Full-run uses h+overcalls, overshooting INTO the U-shape's steep
expensive tail (hp8=0.197, hp10=0.324, hp12=0.688). => full-run systematically OVER-charges
genuine over-calls relative to the calibrated model. The investigation's principled claim
(durable record (d), "the full read run is UNAMBIGUOUSLY the right context; per-cut is NEVER
more correct") is FALSE: the calibrated axis is the REFERENCE run (shorter), so per-cut's
truncation can land CLOSER to the calibrated cost than full-run does.
- CHECKPOINT: calibration-axis mismatch identified; quantifying next.

## KEY FINDING #2 — the insertion RATE FALLS with reference run length (physics)
From penalty_scores.tsv, I/AT rows (rate_mean = calibrated per-ref-column ins prob):
  hp1 .0012 | hp4 .0074 | hp5 .0081 | **hp6 .0098 (PEAK RATE)** | hp7 .0078 | hp8 .0073
  | hp9 .0055 | hp10 .0044 | hp11 .0041 | hp12 .0021 | ...
=> the per-column INSERTION rate PEAKS at hp6 and DECLINES for longer reference runs
(long HPs UNDER-call: the D/AT rate CLIMBS hp8 .073 -> hp12 .127). penalty = c/rate is
therefore U-shaped with a MINIMUM at hp6 and a rising, jittery tail. The investigation's
physical intuition ("charge the over-call at its true physical run length -> MORE
expensive is more correct") is BACKWARDS w.r.t. the data: longer runs have LOWER
per-column insertion rates. Full-run's higher cost comes NOT from a real higher error
probability but from (i) indexing by the read run (=ref+overcalls) and (ii) the
reciprocal-rate penalty amplifying the sparse tail.
count_total (calibration reliability): hp6=354206, hp8=48609, hp10=7052, hp12=**2165**,
tail jittery (12=.688, 13=.535, 14=.732). => full-run routes genuine over-calls (common,
at well-sampled short/mid ref runs) into a RARE, NOISY tail bin.

## WITNESS (worktree flag toggled in-process) — PERSISTED
Env note: the _USE_FULL_RUN_INS flag lives ONLY in the worktree copy
(/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc/rectify), NOT the
installed main checkout (/Users/kevinroy/work/rectify/rectify). `python -c` from the
worktree root resolves to the worktree copy (flag present); `python script.py` puts the
script dir on sys.path[0] and resolves to the main checkout (NO flag) -> flag looks inert.
witness_model.py now forces the worktree path (sys.path.insert). All numbers below are the
worktree/flag-active copy, cross-checked by a manual per-k DP loop.

W-model: genome exon2 starts with a REFERENCE 8-A run; read over-calls to 12 A (4 excess),
then an A-free tag matching exon2 downstream. Score the TRUE acceptor:
- flag OFF (per-cut) : **0.7772**  (k-sweep splits: t1=0 at k=4, t2 = 4 x ins_cost(4)=0.1943)
- flag ON  (full-run): **2.7528**  (k=0: t1 = 4 x ins_cost(12)=0.6882; run-splitting removed)
- CALIBRATION-CORRECT (ref-run=8) : **0.7888** (= 4 x ins_cost(8)=0.1972)
=> per-cut (0.7772) sits ~1.5% from the calibrated cost (0.7888); full-run (2.7528)
   OVER-charges by **3.49x/base (3.54x total)**. FALSIFIES the investigation's absolute
   claim ("per-cut is NEVER more correct than full-run; full read run is UNAMBIGUOUSLY the
   right context"): for a GENUINE over-call — the only regime ins_cost fires — per-cut is
   MEASURABLY closer to the calibrated cost than full-run. (Per advisor: per-cut is not
   "correct" either — it's off-axis via truncation; the near-match is mechanism, not
   principle. The narrow, solid claim is that "NEVER more correct" is false.)
W1 control (all-A(12), A-free genome, pure gaming): OFF=1.7604, ON=8.2584 (reproduces the
  investigation's fabrication witness EXACTLY -> toggle verified).
- CHECKPOINT: witness numbers persisted.

## KEY FINDING #3 (CENTERPIECE) — false dichotomy; a THIRD axis dominates full-run
The investigation frames the choice as per-cut vs full-run and defends full-run as a
NET-POSITIVE tradeoff (fair panel: 35 anti-fabrication wins 4->0 vs 6 losses 0->4). That
tradeoff is FALSE because a third option was never costed:
  INDEX ins_cost BY THE REFERENCE HP RUN AT THE INSERTION COLUMN — exactly what del_cost
  ALREADY does (_precompute_del_costs uses genome_seq HP run; ins uses _hp_run_length(query)).
- This is CUT-INDEPENDENT (it's a per-REFERENCE-COLUMN vector over the fixed ref window,
  like del_costs[j]), so it unlocks the SAME single-pass concat DP (~30x) that full-run
  unlocks. The DP-unlock is therefore NOT unique to full-run and cannot justify full-run's
  MODEL choice.
- It is CALIBRATION-CORRECT (the table's own axis), so it removes the run-splitting gaming
  (the anti-fabrication wins) WITHOUT the 3.5x over-penalization -> it does NOT create the
  0->4 losses. The investigation's conceded "6 losses / 4-6 reads/panel" are therefore NOT
  the intrinsic price of de-gaming; they are an AVOIDABLE artifact of choosing the wrong
  cut-independent axis.
=> Full-run is DOMINATED: same DP-unlock, strictly worse error model, and it manufactures a
   class of false canonical->non-canonical demotions that reference-run avoids. "Full-run is
   the sounder model" is WRONG; the honest statement is "full-run is cut-independent but on
   the OUTCOME axis (read run incl. over-calls), dominated by the untested per-reference-
   column axis (cause)."
Reinforcing (hp_penalty.py ~183/305, code's own comment): penalty_score is "a reciprocal-
rate heuristic, c/rate_mean, NOT -logP, and incoherent to sum in an additive DP." Neither
per-cut nor full-run is a coherent probabilistic error model; "sounder error model" is
ill-posed, and full-run merely sums LARGER incoherent tail values.
- CHECKPOINT: centerpiece (reference-run dominance) persisted.

## GUARD SUB-QUESTION (does full-run break junction-slide-within-HP?)
NO break: the investigation's own R3 0.589->0.589 and guard-fires-under-full-run
(fair +0.117) show full-run pushes the SAME direction as the guard, not against it. The
real, narrower point: hp_drift_margin=3.0 is an ABSOLUTE additive constant (junction_refiner
line 734: eff_margin += hp_drift_margin; veto if best > incumbent - eff_margin). Full-run
rescales HP-insertion scores ~4.7x (W1 1.76->8.26), so the score DELTAS the 3.0 margin
gates against change scale. "Flip the default to full-run and keep 3.0" leaves the guard at
an UNKNOWN effective stringency; the investigation's own revalidation #5 says the margin
must be re-swept, yet it still calls the switch "safe." That is internally inconsistent:
safe-to-switch cannot be asserted BEFORE the re-tune it also demands. Frame = "not safe
without re-tune," not "breaks."
- CHECKPOINT: guard interaction persisted.

## VERDICT (MODEL-CORRECTNESS LENS) — FAULT FOUND; "full-run is sounder" does NOT survive
fault_found = TRUE.

CERTIFICATION REQUESTED WAS: "full-run IS the sound Nanopore error model." DENIED.

The fault: the penalty table's ins_cost(hp, base) is CALIBRATED by the REFERENCE
homopolymer run length flanking the insertion column (profiler Phase 5, lines 1320-1322:
max(hp_arr[pos-1], hp_arr[pos]) on the REFERENCE sequence; docstring line 15). Full-run
indexes it by the READ run length measured over the whole rescue window — which, for a
genuine over-call, equals reference_run + (#over-called bases), overshooting the table INTO
its sparse, jittery, expensive tail. Empirically (worktree flag, in-process): a real 4-base
poly-A over-call at an 8-A reference acceptor is charged 2.7528 by full-run vs the
calibration-correct 0.7888 (ref-run=8) — a 3.5x OVER-penalization. This directly falsifies
the investigation's load-bearing principled claim that "per-cut is NEVER more correct than
full-run" and "the full read run is UNAMBIGUOUSLY the right context." In the over-call
regime (the ONLY regime ins_cost fires), per-cut (0.7772) is closer to the calibrated cost
than full-run.

Physics compounds it: the calibrated per-column insertion RATE PEAKS at hp6 and FALLS for
longer reference runs (long HPs under-call, they don't over-call); full-run's higher cost is
an artifact of the read-run index + reciprocal-rate tail noise (hp12 n=2165, jittery), not a
real higher error probability.

Centerpiece (why this DISQUALIFIES the switch, not just dings it): the investigation poses a
FALSE dichotomy. A third axis — index ins_cost by the REFERENCE run at the insertion column,
exactly as del_cost already does — is BOTH cut-independent (unlocks the identical ~30x
single-pass DP) AND calibration-correct (removes the run-splitting gaming WITHOUT the 3.5x
over-penalization, so it does NOT manufacture the 0->4 canonical->non-canonical demotions the
investigation concedes). Full-run is therefore DOMINATED: it shares the DP-unlock (so the
speedup cannot justify its model choice) but is strictly worse as an error model and creates
avoidable false demotions. The investigation's "net-positive tradeoff (35:6)" is illusory —
the 6 losses are the direct symptom of the wrong axis and are avoidable.

Safe-to-switch: also fails. (i) Model is not sounder (above). (ii) hp_drift_margin=3.0 is an
absolute constant tuned against per-cut score deltas; full-run rescales HP scores ~4.7x, so
flipping the default while keeping 3.0 leaves the guard at an unknown effective stringency —
the investigation itself demands a margin re-sweep (revalidation #5) yet calls the switch
safe, which is inconsistent. (iii) Code's own comment: penalty_score is a reciprocal-rate
heuristic, "NOT -logP, incoherent to sum in an additive DP" -> "sounder error model" is
ill-posed for BOTH variants; full-run just sums larger incoherent tail values.

What survives from the investigation: (1) the per-cut run-splitting IS a real gaming artifact
worth removing; (2) cut-independence IS the property that unlocks the DP; (3) the fabrication
witness is real. What does NOT survive: that FULL-RUN is the correct way to achieve (1)/(2).
The reference-column axis achieves both, correctly, and should be evaluated before any
default flip.

Recommendation: do NOT flip _USE_FULL_RUN_INS to default-True. Prototype reference-column
ins indexing (mirror _precompute_del_costs on the ins side) and A/B it against full-run on
the same mix_fair_out / mix_r3b_out panels; expect it to keep the 35 anti-fabrication wins
while eliminating the 4-6 reads/panel 0->4 losses. Re-tune hp_drift_margin under whichever
axis is chosen BEFORE claiming safe-to-switch.
- CHECKPOINT: verdict written. AUDIT COMPLETE.
