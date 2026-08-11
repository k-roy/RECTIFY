# INSCOST 3-way ins_cost synthesis — per-cut vs full-run vs reference-column

Adjudicator synthesis over 3 audits (model-correctness, A/B-measurement-validity,
revalidation-and-coupling). Builder returned null — no builder self-report; synthesis
rests on the 3 auditor records + their durable .md files.

## Plan
1. Read the 3 auditor .md files for the concrete numbers/witnesses. [in progress]
2. Read dev/INSCOST_SYNTHESIS.md for the 8-item re-validation checklist.
3. Decide WINNER (ref-column / full-run / keep per-cut) + whether auditors agree.
4. Decide safe-to-proceed vs blocking fault.
5. Ordered re-validation checklist for the winner.
6. The one result most likely to change.

## Checkpoints
- [init] File created, plan written.
- [read] All 3 refcol audits + BUILD.md + prior INSCOST_SYNTHESIS.md (8-item checklist) read.
  Verdict is UNANIMOUS: ref-column FALSIFIED. Writing synthesis.

---

# VERDICT

## (1) WINNER: FULL-RUN (adopt as reference model, gated flag). Reference-column is ELIMINATED.
Keep per-cut as the safe status-quo hold until full-run's checklist is green. Do NOT flip
default to ref-column under any circumstance.

### Do the auditors agree? YES — unanimous against ref-column, on 3 independent lenses.
| Audit lens | fault_found | verdict survives | winner from that lens |
| --- | --- | --- | --- |
| refcol model-correctness | TRUE | NO (falsified) | Neither full-run nor refcol; certifies refcol must NOT be adopted |
| A/B measurement-validity | TRUE | NO (falsified) | Full-run (or per-cut) on discovery axis; refcol regresses |
| revalidation + coupling | TRUE | NO (falsified) | Refcol eliminated; full-run credible candidate, per-cut safe hold |

All three converge: the refcol recommendation ("keep full-run's ~35 anti-fab wins AND
eliminate the 4-6/panel 0->4 canonical-demotion losses") is EMPIRICALLY FALSE. All three
independently re-confirm the builder's own 3-way A/B:
- percut->fullrun (fair, no-guard): 42/5600 = 0.75% churn; tier 4->0 = 35 wins, 0->4 = 6 losses; R0flank +0.005
- percut->refcol  (fair, no-guard): 694/5600 = 12.39% churn (16x); tier 4->0 = 358, 0->4 = 196, 4->4 = 118; R0flank +0.011

Ref-column does NOT eliminate the 0->4 demotion class — it BALLOONS it 33x (196 vs 6).

### The three convergent reasons ref-column loses:
1. **Structural (model-correctness):** refcol charges insertion by GENOME-COLUMN context ALONE
   with no term protecting the canonical dinucleotide. Placing the acceptor DEEPER in a genome
   poly-A run (cheap in-run insertions) is always <= placing it at the run EDGE (true canonical
   acceptor), so the minimizing DP slides the acceptor into the run (slide-delta +1..+10) and
   fabricates a non-canonical AA junction. 133/196 (68%) of demotions are base-MATCHED A-into-A
   slides — a fully base-AWARE ins_cost would STILL demote ~133 (still 22x full-run's 6). NOT a
   fixable base-identity bug; the demotion class is structural. Secondary (minority, 61/196=31%):
   base-identity-independence hard mis-charge (T=46,G=36,C=2 charged the poly-A/T discount;
   witness: 4 inserted C's next to a genome 8-A run score 0.7888 = 4*ins_cost(8,'A'), identical
   to 4 genuine A over-calls). Ref-column RELOCATES per-cut's run-splitting gaming into
   column-choice (drag-into-run) gaming rather than removing it.
2. **Discovery regression (A/B validity — DECISIVE):** the refcol recommendation was graded on
   the FAIR panel, where EVERY true junction is CANONICAL (0 non-canonical truths) — a panel that
   structurally CANNOT observe non-canonical discovery. On the r3b R3 rung (the real discovery
   target, 1600 reads), measured from IDENTICAL code: full-run is exactly NEUTRAL (R3 net 0,
   0.258->0.258, GAINED=0 LOST=0) while ref-column REGRESSES (R3 net -33 reads, 0.258->0.237,
   -0.021; GAINED=9 LOST=42). Per distinct junction: refcol net-worse on 3 of 4 R3 junctions,
   better on 0 — survives the per-read->per-junction unit correction. Full-run has no such regression.
3. **Coupling / re-validation burden:** refcol's scale is BIMODAL by design (≈per-cut for genuine
   over-calls: Wm refcol 0.7888 vs per-cut 0.7772; ≈full-run in the gaming/A-free regime: W1 refcol
   6.8748 vs per-cut 1.7604, ~4x). Its 12.39% churn makes it behaviorally the FARTHEST of the three
   from per-cut, demanding the MOST re-validation, not the least. The 196 demotions are also not
   curable by a prior re-tune (a subset fire with the prior OFF at tier<4 incumbents).

### Where ref-column IS right (conceded, does not rescue it):
For an ISOLATED pore over-call, refcol's per-column charge is biologically correct in BOTH task
sub-cases: 4 A's adjacent to a genome NON-HP column = 4.58 (expensive, correct — no run to explain
them); 4 A's inside a genome poly-A(8) = 0.79 (cheap, correct — a run's own-base stutter). The DP
correctness was verified against a 2D brute force (both = 0.7888, 1e-9). And it is genuinely
cut-independent (unlocks the single-pass DP). But locally-correct per-column cost becomes a
DRIFT-INDUCER in junction comparison. Axis-match is NECESSARY, not SUFFICIENT.

## (2) SAFE TO PROCEED? NO — there is a BLOCKING FAULT for the recommended-under-audit variant.
The variant under review (reference-column) is FALSIFIED and must not be adopted; adopting it
would introduce a real -0.021 non-canonical discovery regression on the project's primary
objective. That is a blocking fault against ref-column.

For the ACTUAL winner (full-run) the picture is: no BLOCKING fault, but NOT yet default-ship-safe.
Full-run remains a gated adoption exactly as the prior INSCOST_SYNTHESIS concluded — land the
flag, keep default OFF, complete the checklist below before flipping default-ON. Safe to PROCEED
to full-run's re-validation; NOT safe to ship default-ON, and NOT safe to adopt ref-column at all.

safe_to_proceed (schema) = TRUE — the schema boolean pairs with `winner`=full-run, and Q2 asks
"safe to proceed TO RE-VALIDATION." Full-run has NO blocking fault and is the gated-but-sound
path, so it is CLEARED to enter the re-validation checklist below (matching the prior
INSCOST_SYNTHESIS "PROCEED-TO-SWITCH, GATED" posture = proceed-with-gating). The BLOCKING FAULT
is against REF-COLUMN adoption specifically (it would introduce the -0.021 R3 discovery
regression) and lives in final_call — it does not gate full-run's re-validation track.

## (3) ORDERED RE-VALIDATION CHECKLIST FOR THE WINNER (full-run)
Carried from INSCOST_SYNTHESIS's 8-item list, adjusted for full-run (ref-column-specific items dropped):
1. Commit the full-run flag on / rebase to the ACTUAL merge target (master); diff the benchmark
   worktree (carries motif_blind + both uncommitted ins flags) vs master (has neither) for any
   OTHER scoring-relevant divergence. Re-establish the baseline on the branch the switch lands on.
   [All A/B ran on the benchmark worktree — cross-branch equivalence is UNVERIFIED.]
2. FULL suite flag-ON on the merge target (`pytest -m "not slow"`, ~1603) incl. cat1-cat9 fixtures.
   Local green gate. (Flag-ON is only PARTIALLY suite-tested today.)
3. numba-ON DP path, flag-ON, on the CLUSTER build: verify the reversed-slice list->float64 array
   feeds `_score_hp_dp_numba` correctly BEFORE trusting any cluster run. (Full-run has a numba path
   structurally present; refcol did NOT — its column-indexed R+1 ic vector is wrong for the
   Q-indexed kernel, another reason refcol's DP-unlock was never actually written.)
4. `_CANONICAL_HP_PRIOR` re-confirm on the full-run scale: verify 0.5 still = intended noise floor
   (anchor = one HP-DELETION equivalent; DEL cost is UNCHANGED by the flag, so the anchor HOLDS —
   only the prior's effective strength vs HP-insertion-driven gaps erodes ~3.5x). Goal is "still
   the right noise floor," NOT "beat per-cut" (per-cut is a degeneracy; tuning to beat it overfits).
   Confirm NO net discovery regression on fair + r3b, no-guard (prior always-on = shipping proxy).
5. Guard `hp_drift_margin` re-sweep (task #16): smallest zero-discovery-cost margin may DECREASE
   below 3.0 since full-run does part of the drift suppression itself. NOTE: guard-ON is NOT a
   safe substitute for re-tuning — the guard MASKS the ins-cost differences.
6. del_open / arm-F verdict RE-RUN under the full-run scale (DEL_OPEN_DELTA_FINDING). The per-cut-era
   verdict was reached on the 1.76 scale; scale shifted ~4.69x (1.76->8.26) — do NOT assume it carries.
7. Yeast DRS real-data HP-drift transfer (task #17, Sherlock).
8. Human ONT DRS transfer (task #18, GENCODE truth) — decisive external-validity test; direction
   unknown; watch the 0->4 canonical-demotion balance.
Flip default-ON to ship ONLY after 1/2/3/4 are green and 7/8 show no regression.

## (4) THE ONE RESULT MOST LIKELY TO CHANGE
The **HUMAN ONT DRS transfer (checklist item 8)** — the go/no-go bottom line most likely to move.
Every neutral/positive full-run number to date is SIM (fair + r3b), and both panels are
low-diversity for non-canonical junctions (fair = 0 non-canonical truths; r3b = only 4 distinct R3
junctions over 1600 reads). Human HP-length distribution differs from yeast/sim, so the
fabrication/recovery balance could shift; direction is genuinely UNKNOWN (not assumed worse). The
specific quantity to watch is the 0->4 canonical-demotion count under shipping always-on
_CANONICAL_HP_PRIOR on real human reads.

SEPARABLE (most likely to LITERALLY flip, but not the go/no-go): del_open / arm-F (item 6) — an
existing sub-verdict decided on the old 1.76 scale, now facing a ~4.69x rescale; it "must be re-run,
not assumed to carry." It is a separable sub-verdict, not the bottom line.

## NOTES ON THE BUILDER (null)
No builder self-report was returned, but the builder's durable INSCOST_REFCOL_BUILD.md IS complete
(the durability policy worked): it records the design, code-written checkpoints, DP-vs-brute-force
correctness check, byte-identity of the default-OFF path (46 tests), and the decisive fair no-guard
3-way A/B with the FALSIFIED verdict already called. All three auditors independently reproduced
the builder's headline numbers exactly. The synthesis rests on a complete, cross-verified record.

STATUS: COMPLETE.
