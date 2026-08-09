# MICROHOM AUDIT V5 — Task "byte-identity-architecture" — AUDITOR B

**Status:** COMPLETE — VERDICT: CLEAR on byte-identity; design-debt HOLD-note on architecture (non-blocking; no correctness fault)
**Auditor:** B (independent; no coordination with A)
**Mandate:** (a) Re-verify guard+cap+positional stack is byte-identical at defaults
(all of microhom_drift_margin / hp_drift_margin / drift_near_tie_cap / drift_positional_gate = 0.0).
Independent BAM diff vs pre-close parent (sequential raw-order SHA256 + parallel position-sorted).
pytest -m "not slow" green. Confirm `_positional_signal` / `_semiglobal_ed` never called at default.
(b) ARCHITECTURE critique: 4th stacked guard param vs fixing the scorer; second alignment redundancy;
gate/cap interaction; complexity/maintainability vs benefit on a default-off feature.

Working READ-ONLY in worktree agent-a25a2c1e784ad37dc.
Python: /Users/kevinroy/miniconda3/bin/python

---

## CHECKPOINT LOG (append-only)

### CP0 — record initialized
- Scratch: /private/tmp/claude-501/-Users-kevinroy-work-rectify/798e5203-2185-4c27-87b4-f4e50c56c3ba/scratchpad/audit_v5/byte-identity-architecture_B
- No prior record/harness. Fresh start.
- Next: read junction_refiner.py (the CLOSE + drift guards), git log the named commits.

### CP1 — STATIC byte-identity analysis DONE (before empirical run)
Commits: CLOSE=cd3de46, PARENT=4440c2c. CLOSE touched ONLY junction_refiner.py (executable)
+ 4 docs. No default was changed.

**Two INDEPENDENT reasons the positional/ed path is DEAD at default:**
1. Veto block ENTRY guard (line 915): `if moves and eff_margin > 0.0 and incumbent_score is not None`.
   At default hold_margin=0.0, hp_drift_margin=0.0, microhom_drift_margin=0.0 → eff_margin=0.0 →
   `eff_margin > 0.0` is FALSE → the ENTIRE veto/cap/gate block is skipped. So even the cap +
   `_effective_veto_margin` are unreached at default; `_positional_signal`/`_semiglobal_ed` never called.
2. GATE guard INSIDE the block (line 921): `if drift_positional_gate > 0.0:`. Default 0.0 → dead branch
   → `spared` stays False → `if not spared: moves=False; profile.inc('move_margin_vetoes')` — this is
   BYTE-IDENTICAL control flow to the parent's `moves=False; profile.inc('move_margin_vetoes')`.

**Refactor diff (the ONLY change to an existing executable path):** the parent's 2-line veto body
(`moves=False` + profile inc) was wrapped in `spared=False; if gate>0: ...; if not spared: <parent body>`.
At gate=0.0 this is provably equivalent. The new helpers `_semiglobal_ed`/`_positional_signal` are
pure additions, called ONLY from the gate>0 branch.

**Wiring:** real caller = correct_command.py:746 refine_bam_junctions(...) passes NO drift kwargs →
all default 0.0. Forwarded verbatim through refine_bam_junctions → _run_sequential/_run_parallel →
refine_read_junctions. UNWIRED confirmed (grep of commands/ + bam/ finds zero drift-kwarg passes).

Static verdict: byte-identical at default is SOUND by construction. Next: empirical confirmation
(pytest + real BAM diff parent-vs-CLOSE, sequential SHA + sorted).

### CP2 — EMPIRICAL BAM diff DONE (validation_reads.bam, n_workers=1, sort_and_index=False)
Harness: scratch/run_refiner.py injects a chosen junction_refiner.py version as the canonical
module, runs refine_bam_junctions at ALL-DEFAULT drift kwargs on real S288C validation reads
(bgzipped genome via pysam.FastaFile, GTF annotation). Instruments _semiglobal_ed/_positional_signal
with call counters.

Ran THREE versions: PARENT=4440c2c (pre-close), CLOSE=cd3de46 (the close), HEAD=a97ff6d (working
tree = CLOSE + docstring-only expansion; confirmed a97ff6d's .py delta is docstring-only).

STATS all three IDENTICAL: total=36, n_op_reads=16, refined=4, unchanged=32, errors=0.
HELPER_CALL_COUNTS on CLOSE & HEAD: _semiglobal_ed=0, _positional_signal=0  ← proven never called
at default (PARENT lacks the fns; harness only wraps if present).

BAM diff (36 records):
- (A) SEQUENTIAL raw-order SHA256 (order-sensitive): all three = 3d93cd22ed480b414748fbffdd5798c5ac0a9fbd7a54b15c72afb9aa1ca4c862
- (B) PARALLEL position-sorted SHA256 (order-insensitive): all three = e976505c5dab1dd9c102bd40895a5349276dfc8d59224f613ce4422812f18473

=> BYTE-IDENTICAL at default, empirically, sequential AND parallel-sorted. Next: parallel path
(n_workers>1), a 2nd BAM (upf1d), and a POSITIVE CONTROL (gate>0 must DIFFER — proves the harness
can detect a change).

### CP3 — 2nd BAM + POSITIVE CONTROL DONE (harness proven non-blind)
- upf1d BAM (validation_reads_upf1d.bam), default n_workers=1: PARENT==CLOSE byte-identical,
  sorted_sha256=c705abfe3be9b797eda6911d7937fa156a5d8f88b2313c4710c2500d162ea164 (n=32).
  Helpers 0 calls.
- PARALLEL n_workers>1 via module-injection HANGS (multiprocessing fork re-imports the REAL module
  in children, not my injected one → not a valid injection test; harness limitation, NOT a code
  fault). Parallel default-equivalence instead follows from: (i) sequential empirical identity,
  (ii) the CLOSE's _run_parallel diff being pure default-kwarg forwarding (verified statically),
  (iii) the in-tree parallel test coverage. NOT independently BAM-diffed by me for n_workers>1.
- POSITIVE CONTROL (scratch/positive_control.py, PG_* fixture, force veto+gate path on CLOSE & HEAD):
    PG_NE=170 (incumbent), PG_JE=176 (cryptic).
    gate=0.0 → acceptor=170 (VETOED), helper_calls {_semiglobal_ed:0, _positional_signal:0}
    gate=1.0 → acceptor=176 (SPARED),  helper_calls {_semiglobal_ed:2, _positional_signal:1}
  => The harness + instrumentation DO detect a difference when the gate engages. So the
  'identical-at-default / helpers=0' results are a TRUE NULL, not a dead harness. The behavioral
  delta is gated PURELY by drift_positional_gate>0.0.
- Note on my first naive positive control (hold_margin=5, gate=1 on the real BAM): did NOT fire —
  hold_margin alone opens eff_margin>0 but does NOT drift-FLAG a move, and no real read hit a
  drift-flagged near-tie move. This is why the fixture (microhom_drift_margin=8.0) is the correct
  positive control. Lesson: eff_margin>0 is necessary but the gate only runs inside the drift-veto.
Next: pytest -m "not slow".

### CP4 — pytest + ARCHITECTURE probes (part b) — in progress
PYTEST:
- Targeted refiner/guard suites (test_microhom_drift_guard + test_junction_refiner + test_hp_drift_guard
  + test_junction_scoring_parallel): 85 passed, 17 skipped in 5.5s. GREEN.
- Full `pytest -m "not slow"`: launched (PID 26283). Slow tail (test_validation_reads*_upf1d spawn real
  rectify.cli correct subprocesses; ~1700 tests). NOTE: I accidentally spawned duplicate full-suite runs
  via monitor loops → CPU contention on the 8GB M1; killed the duplicates (29196/29194) + monitor wrappers
  (27892/29141), kept original 26283. One `E` seen early: test_restore_polya_from_parquet.py — pending
  classification (likely a parquet/pyarrow env dependency, UNRELATED to the refiner). Will confirm final
  summary + that E when the run lands.

ARCHITECTURE (part b) — the redundancy question (scratch/redundancy_probe.py, escape_probe.py):
- On the PG showcase fixture (the CLOSE's own positive example): _positional_signal = 6 (ed_inc=6,
  ed_mov=0). The SCORER'S OWN delta_improve = s_inc - s_mov = 6.0 — IDENTICAL sign AND magnitude.
  t1(k=0) == _semiglobal_ed exactly on both arms (6.0/0.0). min_k COLLAPSE gap = 0.0 → the free-k
  soft-clip 'escape' did NOT fire on the design's showcase case. Here the 2nd alignment is REDUNDANT.
  Also: the fixture is vetoed only because microhom_drift_margin=8.0 (artificial); at the shipped
  operating point cap=2, delta=6 > cap so it would NOT be drift-vetoed at all → gate moot for it.
- 400-case true-cryptic sweep (spanning delta): scorer delta_improve and ed_signal are ~the SAME
  distribution (both median 5, min 1). The 'escape' (min_k soft-clips incumbent, gap>0.5) fires in only
  ~4.8% of true cryptics; median gap = 0.0. => For the DISCOVERY side, delta_improve already carries the
  signal; ed_signal is largely redundant there.
- STILL TO TEST (decisive): the FAB (error-driven drift) side. The CLOSE's value is separating fab
  (matches incumbent, esig<=0) from true cryptic (esig>0) IN the overlap band. Redundancy verdict hinges
  on whether delta_improve and esig DISAGREE for fab reads. Extending the probe now.

### CP5 — PYTEST FINAL + FAB redundancy probe
PYTEST (full `-m "not slow"`, PID 26283, ran 16m26s to completion, harness exit 0):
  **1670 passed, 39 skipped, 4 deselected, 1 xfailed, 1 ERROR**.
  The 1 ERROR = tests/test_restore_polya_from_parquet.py::test_restore_cat3_plus_2 — AssertionError
  at setup of module fixture `parquet_row`: missing validation-bundle data file
  scripts/validation_data/rebuild_2026_05/trimmed/validation_reads_polya_trim_metadata.tsv (NOT
  git-tracked in this worktree). VERIFIED INDEPENDENT of the CLOSE:
    * that test references junction_refiner/refine_ ZERO times;
    * it was NOT modified in 4440c2c..HEAD;
    * it errors on missing DATA, deterministically, on PARENT and HEAD alike.
  => byte-identity pytest gate SATISFIED: all refiner/guard/validation tests green; the lone error is a
  pre-existing environmental data-provisioning gap, orthogonal to the guard. (17 skips in the refiner
  subset are pre-existing conditional skips.)

FAB REDUNDANCY PROBE (scratch/fab_redundancy_probe.py — my INDEPENDENT reproduction of the panel
construction, CRY=truth-move vs FAB=truth-stay, mh>=0.6, ONT sub+indel errors, n=600 each):
  - Scorer delta_improve ALONE separates CRY (median 6, min 2) from FAB (median 0, max 1) at
    balanced-accuracy = 1.000 — IDENTICAL to ed_signal (1.000), WITH indels.
  - In my reproduction there is NO real delta overlap band (only 3-27 FAB in 0<delta<=2; those have
    esig<=0, same call as delta). => In MY construction the scorer already discriminates; _positional_signal
    is REDUNDANT.
  - TENSION with panel CP2 (claims real overlap: cry ~24% <3, fab ~97% <2). My model did NOT reproduce
    that overlap. Running the panel's OWN dev/discovery_loss_panel.py to reconcile: does the panel's own
    `delta` separate as well as `esig` on ITS construction? (in progress, PID 31297). NOTE: this is a
    SIZE-OF-BENEFIT question (architecture part b), NOT a byte-identity question — byte-identity is SETTLED.

### CP6 — RECONCILED with the panel's own run (I CREDIT the panel; my probe was too kind to the scorer)
Ran dev/discovery_loss_panel.py --n 40 --seed 1 (5760 reads). It DOES produce a real overlap band my
simpler probe missed:
  - DISCOVERY (cry) delta: median 4.0, but 24% <3, 10% <2 (genuine low-delta tail).
  - FAB drift delta: median 1.0, 97% <2.  => cry & fab OVERLAP in delta in [0.5,2].
  In the overlap band [0.5,1.5] where aggregate delta is ~at chance:
    naive positional : balanced-acc 78%
    ed_signal (CLOSE): balanced-acc 99% (cry 97% / fab 100%)   <- SEPARATES where delta cannot.
  cap-alone: disc-loss 9.8% / fab 2.9% ;  WIRED m3/c2/gate1: disc-loss 0.7% / fab 4.3% (== shipped, matches
  DISCOVERY_LOSS_PANEL_RESULT.md CP5 within rounding: doc says 0.4%, my seed-1 run 0.7% — same order).
  => The panel's core claim REPRODUCES. _positional_signal is NOT fully redundant: in the overlap band it
     carries real orthogonal info the scorer's min_k aggregate discarded. My earlier "redundant" reading
     held only because my construction lacked the low-delta cry tail. HONEST CORRECTION logged.

REFINED redundancy point (the ACCURATE architecture critique): ed_signal's MOVED arm
_semiglobal_ed(rescue, g[je:]) is structurally t1(k=0) of _score_junction for the moved candidate; the
INCUMBENT arm is t1(k=0) for the incumbent. The scorer COMPUTES t1(k=0) internally (its k=0 loop iteration)
then DISCARDS it by taking min_k (the free-k soft-clip 'escape'). So the CLOSE re-derives a k=0-anchored
score the scorer already forms but throws away. The architecturally-cleaner fix is to have the SCORER expose
an anchored (k=0 / no-free-prefix) score, not to bolt on a 2nd full O(W*W) alignment as a separate 4th gate
parameter. Caveats to verify: (a) _semiglobal_ed uses UNIT costs vs t1's HP-aware penalty_table costs — same
on the flat-cost path? (b) the anchored incumbent score must also be formed for the INCUMBENT placement.
Verifying (a) next.

### CP7 — cap/gate ordering + t1(k=0) vs _semiglobal_ed (architecture correctness)
CAP/GATE ORDERING (scratch/cap_gate_ordering.py, PG fixture, delta_improve=6, mh flag trips):
  margin=8, cap=0 : gate_off HELD(veto)   / gate_on MOVED(spared)
  margin=8, cap=2 : gate_off MOVED(spared)/ gate_on MOVED(spared)   <- CAP alone spares (delta 6 > cap 2)
  margin=8, cap=6 : gate_off MOVED        / gate_on MOVED           <- boundary (strict > )
  margin=8, cap=10: gate_off HELD(veto)   / gate_on MOVED(spared)   <- cap>=eff inactive; margin 8 vetoes
  => ORDERING IS CLEAN + MONOTONE. veto_margin = max(hold, min(eff, cap)); the move only enters the veto
     branch if within the CAPPED band; the gate can only ADD sparing inside that band, never remove a
     cap-spared move. Cap (delta-axis magnitude) and gate (positional evidence) are complementary and never
     contradict. spare-then-veto within-branch is correct: gate checked BEFORE moves=False.

t1(k=0) vs _semiglobal_ed (redundancy precision): 102/300 (34%) random cases DISAGREE. So _semiglobal_ed
(UNIT edit distance) is NOT numerically the scorer's t1(k=0) (which uses HP-length-dependent del costs from
_precompute_del_costs even on the 'flat' path). => The CLOSE is NOT literally recomputing the scorer's
internal number; it computes a DIFFERENT, simpler unit-cost anchored metric. The redundancy is STRUCTURAL,
not literal: both are 'hard-anchored align rescue->exon2 ref, free suffix', a computation the scorer already
performs (for both candidates) and DISCARDS via min_k. Architecturally-clean alternative: expose the
scorer's own anchored (k=0 / no-free-prefix) score, or add a unit-cost anchored output to the scorer, rather
than a bolt-on 4th gate parameter running a 2nd O(W^2) alignment. (Correcting my CP6 'recomputes' wording:
it re-derives the SHAPE, not the exact value.)

---

## CONSOLIDATED VERDICT (Auditor B)

### (a) BYTE-IDENTITY — CONFIRMED (no fault)
- STATIC: two independent reasons the positional/ed path is dead at default — (i) the veto-block entry
  guard `eff_margin > 0.0` is False when hold/hp/microhom margins all = 0.0; (ii) the inner `gate > 0.0`
  branch is dead, and the fall-through `if not spared: moves=False; profile.inc(...)` is byte-identical to
  the parent's veto body. The CLOSE's only edit to an existing executable path is a behavior-preserving wrap.
- EMPIRICAL: refine_bam_junctions at ALL-default drift kwargs, PARENT(4440c2c) vs CLOSE(cd3de46) vs
  HEAD(a97ff6d), on 2 real validation BAMs:
    validation_reads.bam: raw-order SHA256 3d93cd22.. AND sorted SHA256 e976505c.. IDENTICAL across all 3.
    validation_reads_upf1d.bam: sorted SHA256 c705abfe.. IDENTICAL PARENT==CLOSE.
  Helper call counts at default: _semiglobal_ed=0, _positional_signal=0 (instrumented) — NEVER CALLED.
- POSITIVE CONTROL proves the harness is NOT blind: with gate=1 on the veto-band fixture, output CHANGES
  (acc 170->176) and helpers FIRE (_positional_signal=1). So the default nulls are TRUE nulls.
- PYTEST -m "not slow": 1670 passed / 39 skip / 1 xfail / 1 ERROR. The single error
  (test_restore_polya_from_parquet::test_restore_cat3_plus_2) is a MISSING validation-bundle TSV, references
  the refiner zero times, unmodified in 4440c2c..HEAD → pre-existing + orthogonal. Gate SATISFIED.
- a97ff6d (HEAD) is docstring-only over cd3de46 → runtime-identical.
- PARALLEL PATH (CP8, gap now CLOSED): working tree (=CLOSE) at n_workers=2 (fork = the production Linux
  path), all drift kwargs default, validation_reads.bam → sorted SHA256 = e976505c5dab.. == the proven
  sequential/PARENT sorted hash. Chain: CLOSE-parallel(sorted) == CLOSE-sequential(sorted) ==
  PARENT-sequential(sorted). Requirement (a) parallel diff FULLY met. (The earlier hang was a harness bug —
  injection-Pool on macOS spawn — not the refiner.) No remaining byte-identity caveat.

### (b) ARCHITECTURE — the design is CORRECT but NOT the cleanest; a genuine (non-blocking) critique
1. GATE/CAP INTERACTION: CLEAN. veto_margin = max(hold, min(eff, cap)); a move only enters the veto branch
   if inside the capped band; the gate can only ADD sparing there; never contradicts the cap. spare-then-veto
   ordering verified correct. No fault.
2. SECOND-ALIGNMENT REDUNDANCY: STRUCTURAL, not literal. _positional_signal's two arms are hard-anchored
   'rescue -> exon2 ref' alignments — the SAME SHAPE _score_junction already computes as t1(k=0) for each
   candidate and then DISCARDS via min_k (the free-k soft-clip 'escape'). But numerically _semiglobal_ed
   (unit cost) != t1(k=0) (HP-aware del costs) in 34% of cases, so it is a DIFFERENT metric, not a literal
   recompute. The cleaner design is to make the SCORER expose an anchored (no-free-prefix / k=0) score so
   delta_improve ITSELF discriminates in the overlap band — rather than a bolt-on 4th gate parameter running
   a second O(W^2) DP. This is a maintainability/altitude critique, NOT a correctness fault.
3. THE BENEFIT IS REAL (I initially under-credited it; corrected): the panel's overlap band is genuine
   (cry 24% delta<3, fab 97% <2, overlapping [0.5,2]); in that band ed_signal separates at balanced-acc 99%
   where aggregate delta is at chance. So the positional signal is NOT snake-oil — it recovers discovery
   (9.8%->0.7% loss) the cap alone cannot, at a modest fab cost (2.9%->4.3%). My own simpler reproduction
   failed to make the overlap and mislead me toward 'fully redundant'; the panel's construction is the
   authority and it holds.
4. 4th-PARAMETER / DEFAULT-OFF COMPLEXITY: The guard stack is now 4 stacked knobs
   (hp_drift_margin, microhom_drift_margin, drift_near_tie_cap, drift_positional_gate) + microhom_threshold
   + hp_drift_min_run — all DEFAULT-OFF and UNWIRED (correct_command passes none). This is a lot of surface
   for a feature no production path exercises, validated on SYNTHETIC panels only, with COMPASS real-data as
   the still-unmet enablement gate. The byte-identity discipline is impeccable, but the design debt is real:
   the discrimination logically belongs in the scorer (one anchored score), and 3 of the 4 knobs
   (hp_drift_margin, microhom_drift_margin, near_tie_cap) are superseded in value by the gate per the panel
   (margin/cap alone lose 10-24% discovery; only gate closes it) — they persist mainly as structural scaffolding.
   RECOMMENDATION: before enabling, consolidate the drift-margin+cap+gate cascade into a scorer-level anchored
   discrimination (fewer knobs, discrimination where placement is decided), and gate enablement on COMPASS
   real-data — do NOT ship the 4-knob synthetic-only stack as the permanent interface.

VERDICT: CLEAR on byte-identity (the audited claim). HOLD-worthy DESIGN DEBT on architecture (4 default-off
knobs; discrimination bolted beside the scorer instead of inside it) — non-blocking for the byte-identity
close, but should be resolved before the guard is ever wired/enabled. No correctness fault found.

### CP8 — PARALLEL PATH CLOSED (advisor-flagged gap eliminated)
The n_workers>1 caveat is now RESOLVED. The earlier hang was a HARNESS bug (module-injection Pool launched
at module top level with no __main__ guard on macOS spawn → recursive re-spawn), NOT the refiner's parallel
path failing. Correct test: spawn-safe driver (scratch/run_parallel.py, __main__ guard), NO sys.modules
injection (run the WORKING TREE = shipped CLOSE directly), force fork (the refiner's _run_parallel is
fork-based — populates _WORKER_POOL_STATE in the parent for forked workers; this is the production Linux
path; macOS spawn breaks header inheritance -> KeyError 'header').

RESULT (validation_reads.bam, n_workers=2, ALL drift kwargs default):
  STATS identical: total=36, n_op_reads=16, refined=4, unchanged=32, errors=0.
  PARALLEL sorted SHA256 = e976505c5dab1dd9c102bd40895a5349276dfc8d59224f613ce4422812f18473
  == the PROVEN sequential/PARENT sorted hash.  MATCH: TRUE.
=> Chain closed: CLOSE-parallel(sorted) == CLOSE-sequential(sorted) == PARENT-sequential(sorted).
   The parallel path is byte-identical at default. Requirement (a) FULLY met (sequential raw-order +
   parallel position-sorted, both diffed). (Raw-order is not a meaningful invariant for a reordering
   worker path; sorted is the correct bar — which is why the task pairs 'parallel' with 'position-sorted'.)
