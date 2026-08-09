# MICROHOM-DRIFT GUARD + NEAR-TIE CAP — AUDIT V2 (LENS: BYTE-IDENTITY AT DEFAULT + FIX STRUCTURAL HONESTY)

Auditor: adversarial Opus-Max, agent-a25a2c1e784ad37dc (worktree branch worktree-agent-a25a2c1e784ad37dc).
READ-ONLY. Python: /Users/kevinroy/miniconda3/bin/python (repo env if importable).
Target: commit `05664bc` "feat(refiner): near-tie read-evidence cap on the read-blind drift guards".
Direct parent: `93180ef` (docs-only handoff). The DIRECT CODE parent for junction_refiner.py logic is the
pre-cap guard state; the cap commit changed ONLY junction_refiner.py + tests + the design doc.

## LENS (two duties)
- **A. BYTE-IDENTITY AT DEFAULT** — re-verify the guard+cap is TRULY INERT at default
  (drift_near_tie_cap=0.0, all margins 0.0): output byte-identical to the pre-cap code. Independent BAM diff
  vs the pre-cap commit on sim panels (mix_r3b / mix_fair / fab_sweep), SEQUENTIAL path (SHA256 over full SAM
  record text, raw record order) AND parallel path position-sorted (raw parallel order is known
  guard-independent scheduling noise — do NOT falsely flag). Confirm _effective_veto_margin does not run / has
  no effect at default. pytest -m "not slow" green. Isolate confounders via the DIRECT parent.
- **B. FIX STRUCTURAL HONESTY** — the code+doc CLAIM the cap "bounds but does not close" the read-blind fault
  and "adds no discriminating signal inside (0,cap)". Independently verify HONEST not over/under-claim:
  (i) is delta_improve/eff_margin/cap really ONE scalar axis (=> cap == reducing the effective margin, no new
  signal)? (ii) does the cap EVER protect a fab-drift move it shouldn't (does RAISING the cap re-admit
  fabrication)? (iii) confirm the cap only ever LOWERS the veto margin (never raises it) => can only ADD moves,
  never remove a real discovery relative to no-cap. State plainly whether the docstring honest-framing is right.

## PLAN
1. Orientation: read prior records + cap diff + hot-path + design doc + tests + production caller. (DONE)
2. Static proof of inertness: enumerate the cap function truth table + the veto-block reachability at default.
3. Harness A1 (cap function): exhaustive property test of _effective_veto_margin — monotonic in cap, never
   exceeds eff_margin, never below min(hold,eff)... => proves "only lowers, never raises".
4. Harness A2 (SEQUENTIAL BAM diff): HEAD vs pre-cap commit on mix_r3b / mix_fair / fab_sweep, default config
   AND motif_blind config, n_workers=1, SHA256 over full SAM record text, RAW record order. MATCH required.
5. Harness A3 (PARALLEL BAM diff): same, n_workers>1, position-sorted digest (+ HEAD-vs-HEAD self-consistency
   to show raw-order noise is guard-independent). MATCH required.
6. Harness B (structural honesty): (a) drive refine_read_junctions with cap sweep on a constructed fab-drift
   move and a constructed real-cryptic move; measure whether raising cap re-admits the fab move; (b) prove
   one-scalar-axis: cap only enters via veto_margin subtraction from incumbent_score — same axis as
   delta_improve; (c) prove monotone superset: moves(cap=c) ⊇ moves(no-cap) for all c, and moves grow
   monotonically as cap shrinks toward 0+... actually toward the (hold, eff) window.
7. pytest -m "not slow" (guard-relevant subset + full) green.
8. VERDICT: byte-identity CLEAR/FAULT; structural-honesty HONEST/OVER/UNDER-claim, with the (margin,cap)
   re-admission finding if any.

## CHECKPOINTS (appended as each sub-step lands — persist every number immediately)

### CP0 — orientation DONE (2026-07-12)
- Cap commit `05664bc`, code parent for the veto logic is the pre-cap guard. `git diff --name-only 05664bc^ 05664bc`
  = junction_refiner.py + tests/test_microhom_drift_guard.py + dev/MICROHOMOLOGY_DRIFT_GUARD_DESIGN.md ONLY.
  (Parent 93180ef is docs-only; the true code delta to compare against for byte-identity is a pre-cap
  junction_refiner.py — I will diff HEAD vs 05664bc^ :junction_refiner.py, and cross-check that 05664bc^'s
  junction_refiner.py == the earlier byte-identity audit's reference b6f07f7 guard code + geometry fix.)
- Cap function (junction_refiner.py:535-558): pure. Rule:
    if drift_near_tie_cap > 0.0 and eff_margin > hold_margin: return max(hold_margin, min(eff_margin, cap))
    else: return eff_margin
  Default cap=0.0 → first branch False → returns eff_margin unchanged (IDENTITY).
- Veto block (L834-841): `if moves and eff_margin > 0.0 and incumbent_score is not None:` THEN
    veto_margin = _effective_veto_margin(hold_margin, eff_margin, cap)
    if best_score_cmp > incumbent_score - veto_margin: moves = False
  So _effective_veto_margin is called ONLY when eff_margin > 0.0. At the PRODUCTION default all margins 0.0 →
  eff_margin stays hold_margin = 0.0 → eff_margin>0.0 is False → veto block NEVER entered → cap fn never runs.
- Production caller correct_command.py:746 passes NONE of hold_margin/hp_drift_margin/microhom_drift_margin/
  microhom_threshold/drift_near_tie_cap → all take function defaults (0.0/0.0/0.0/0.5/0.0). Cap not surfaced in
  any CLI. So in production the guard+cap is entirely inert.
- Prior byte-identity audit (dev/MICROHOM_AUDIT_byte-identity.md) CLEARED the PRE-CAP guard's inertness
  (sequential + parallel BAM diff MATCH on the same 3 panels). The cap commit is PURELY ADDITIVE (git diff:
  the only removed line in the veto block is `if best_score_cmp > incumbent_score - eff_margin:` replaced by the
  veto_margin variant; everything else is +). So the byte-identity leg reduces to: at default, veto_margin ==
  eff_margin, and the veto block is unreached in production anyway.
- STATUS after CP0: static case for inertness is very strong; must still (a) confirm empirically via BAM diff
  HEAD vs pre-cap, (b) confirm the cap fn is truly identity at default including the profile counter side effect
  (`near_tie_cap_applied` only inc'd when veto_margin < eff_margin — at default veto_margin==eff_margin so NOT
  inc'd => no profile-output divergence), (c) run the honesty harness.

### CP1 — pre-cap equivalence + SEQUENTIAL byte-identity DONE (2026-07-12)
GIT EQUIVALENCE (load-bearing): junction_refiner.py is IDENTICAL between `17e57ba` (the geometry-fix commit =
the prior byte-identity audit's cleared guard state) and `05664bc^` (the cap commit's parent). i.e. the pre-cap
reference IS exactly the guard code the earlier audit cleared as inert. `git diff --stat 17e57ba 05664bc^ --
rectify/core/splice/junction_refiner.py` = EMPTY. The cap commit is therefore a PURE delta on an already-inert
file. Full-tree diff `05664bc^..05664bc` = junction_refiner.py + tests + design.md ONLY; junction_scoring.py
IDENTICAL parent..HEAD (no concat-DP confounder — cleaner isolation than the prior audit's b6f07f7 base).
The ONLY removed CODE line: `if best_score_cmp > incumbent_score - eff_margin:` → replaced by
`veto_margin = _effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap)` + the same compare on
veto_margin. Everything else additive.

SEQUENTIAL BAM-diff harness (my own: v2audit/harness_cap_byteident.py) — HEAD (cap present, default) vs
05664bc^ blob loaded as sibling module. n_workers=1, RAW record order, SHA256 over full SAM record text.
Three configs incl. the CRUX (guard_on_cap_off = motif_blind + microhom_drift_margin=8.0 + drift_near_tie_cap
=0.0 on HEAD / just the guard on PRE) — this ENTERS the veto block (eff_margin>0.0) so HEAD's
_effective_veto_margin ACTUALLY RUNS and must return eff_margin unchanged. The prior audit never did this (only
margins=0.0 → veto line never executes). RESULTS (all MATCH):
  mix_r3b_out/default        HEAD refined=815 PRE=815  5100 recs  MATCH
  mix_r3b_out/motifblind     HEAD refined=334 PRE=334  5100 recs  MATCH
  mix_r3b_out/guard_on_cap_off HEAD refined=308 PRE=308 5100 recs MATCH   <- veto fired (815->308), still identical
  mix_fair_out/default       HEAD refined=585 PRE=585  5800 recs  MATCH
  mix_fair_out/motifblind    HEAD refined=580 PRE=580  5800 recs  MATCH
  mix_fair_out/guard_on_cap_off HEAD refined=188 PRE=188 5800 recs MATCH  <- veto fired (585->188), still identical
  fab_sweep/default          HEAD refined=152 PRE=152  2100 recs  MATCH
  fab_sweep/motifblind       HEAD refined=136 PRE=136  2100 recs  MATCH
  fab_sweep/guard_on_cap_off HEAD refined=92  PRE=92   2100 recs  MATCH   <- veto fired (152->92),  still identical
OVERALL: BYTE-IDENTICAL. The refined-count DROP in guard_on_cap_off (many reads vetoed) proves the veto block
was genuinely exercised on HEAD; the byte-match proves _effective_veto_margin(...,cap=0.0)==eff_margin exactly.
=> Byte-identity leg (sequential) CLEAR, and STRONGER than the pre-cap audit (exercises the replaced line live).
NEXT: parallel path (position-sorted) + cap-function property harness + structural-honesty harness.

### CP2 — STRUCTURAL HONESTY (pure-fn properties + end-to-end re-admission) DONE (2026-07-12)
Harnesses: v2audit/harness_honesty.py, harness_honesty2.py, harness_readmit2.py (+ probe_delta.py).

PURE-FUNCTION PROPERTIES of _effective_veto_margin (exhaustive grid 12^3 + 200k random draws,
hold,eff,cap ∈ [-10,20]):
  (P1) ONE SCALAR AXIS — signature is exactly (hold_margin, eff_margin, drift_near_tie_cap); NO read/seq
       input. The caller uses the result ONLY as `best_score_cmp > incumbent_score - veto_margin`
       (⟺ delta_improve < veto_margin, delta_improve = incumbent_score - best_score_cmp). Cap adds NO new
       term, NO read feature, NO extra comparison. CONFIRMED (source-grep + signature).
  (P2) CAP ONLY LOWERS, NEVER RAISES — veto_margin <= eff_margin for ALL (hold,eff,cap). CONFIRMED
       (0 violations in 200k+grid). => relative to no-cap (veto=eff) the cap can only ADD moves, never
       remove one => never removes a real discovery relative to no-cap.
  (P3) HOLD IS A FLOOR — on the active branch veto_margin >= hold_margin always. CONFIRMED.
  (P4) CLAMP IDENTITY — on the active branch (cap>0 AND eff>hold), _effective_veto_margin ==
       clamp(cap, hold, eff) = max(hold, min(eff, cap)) EXACTLY (0 mismatch in 200k). i.e. the cap is
       literally "use `cap` as the veto margin, clamped into [hold, eff]".
  (P5) DISABLED/IDENTITY BRANCH — cap<=0 OR eff<=hold => returns eff_margin unchanged (byte-identical to
       the pre-cap veto line). CONFIRMED.
  (P6) MONOTONE within cap>0: veto_margin is NON-DECREASING in cap (hold-floor at small cap up to eff at
       cap>=eff). CONFIRMED.
  (P7) DISCONTINUITY at cap=0 (the enabled/disabled switch): f(cap=0)=eff (NO protection, full veto) but
       f(cap=0+eps)=max(hold,eps) (near-MAX protection). So the param jumps: 0 = OFF, any tiny +val =
       near-max guard-weakening. NOT a fault (docstring says default 0.0 = disabled) but an ERGONOMIC
       sharp edge: "small cap" does NOT mean "small effect" — it means near-max re-admission. Flag for
       whoever tunes it (a user setting cap=0.1 expecting a mild effect gets a near-maximal one).

END-TO-END RE-ADMISSION on refine_read_junctions (harness_readmit2.py) — the decisive mechanism check.
Construction: canonical GT-AG intron, exon2 = partial period-6 repeat (mh(A->A+6)=0.667 >= 0.5 so the
+6 move is drift-FLAGGED) but built from the A+6 frame so A+6 scores STRICTLY better; read carries a
1X boundary error at the acceptor (so the N-op is actually SCORED — the refiner skips clean-boundary
canonical N-ops via boundary_error_window, which is WHY the existing tests' tie-case never demonstrates
a real veto). Direct scorer: score(A)=6.0, score(A+6)=0.0 => delta_improve≈6.0.
RESULTS:
  guard OFF (m=0)        -> acc=176  (MOVED A->A+6; the read genuinely discovers the +6 site)
  guard ON  (m=8,no cap) -> acc=170  (VETOED, held A) — the read-blind fault: a read-supported move suppressed
  cap sweep (m=8, hold=0, veto_margin=clamp(cap,0,8)):
     cap=0.00               -> veto=8.00 -> held A     (DISABLED == no-cap == full guard)
     cap in [0.25 .. 6.00]  -> veto in [0.25..6.00] -> MOVED A+6   (RE-ADMITTED: veto_margin < delta_improve≈6)
     cap in [6.25 .. 8.00+] -> veto in [6.25..8.00] -> held A      (vetoed: veto_margin >= delta_improve)
  FLIP at cap≈6.25 pins delta_improve≈6.0 (matches the direct scorer). 
INTERPRETATION (answers task B(ii) DIRECTION):
  - The cap re-admits a flagged move PRECISELY when clamp(cap,hold,eff) <= delta_improve.
  - It is LOWERING the enabled cap (toward hold) that re-admits the move; RAISING it toward eff RESTORES
    full suppression. So the task's phrasing "does RAISING the cap re-admit fabrication" is DIRECTIONALLY
    INVERTED — raising cap toward eff = STRONGER guard; lowering cap = weaker guard = re-admits. Reported
    as an ADVERSARIAL CORRECTION, not a code fault.
  - CRUCIAL for the honesty claim: the cap re-admits on delta_improve ALONE. A REAL cryptic and a FAB
    drift with the SAME delta_improve are treated IDENTICALLY — both re-admitted iff veto_margin <=
    their (shared) delta_improve. So the cap indeed BOUNDS the read-blind loss (protects strong-evidence
    moves) but adds NO discriminating signal inside the (0, cap) band: within it, real and fab overlap.
    => The docstring's "bounds but does not close / no discriminating signal inside (0,cap)" is HONEST
    and, on this evidence, CORRECT (neither over- nor under-claimed).
NEXT: parallel BAM path (position-sorted) byte-identity + pytest -m "not slow".

### CP3 — PARALLEL byte-identity + test suite DONE (2026-07-12)
PARALLEL BAM-diff (v2audit/harness_par.py, n_workers=3, fork, POSITION-SORTED digest) HEAD vs 05664bc^,
3 configs incl. guard_on_cap_off. ALL 9 cells MATCH:
  mix_r3b_out {default,motifblind,guard_on_cap_off}   pos-sorted HEAD==PRE MATCH (5100 recs)
  mix_fair_out {default,motifblind,guard_on_cap_off}  pos-sorted HEAD==PRE MATCH (5800 recs)
  fab_sweep   {default,motifblind,guard_on_cap_off}   pos-sorted HEAD==PRE MATCH (2100 recs)
  + HEAD-self pos-sorted MATCH (determinism). On THIS run HEAD-self RAW-order also MATCHed (scheduling
  happened deterministic here); either way the position-sorted digest is the guard-independent-safe
  comparison and it MATCHES. OVERALL PARALLEL: BYTE-IDENTICAL.

TEST SUITE:
  tests/test_microhom_drift_guard.py = 20 passed (14 pre-cap + 6 new cap tests incl. TestEffectiveVetoMargin
    + test_near_tie_cap_byte_identical_when_disabled).
  tests/test_junction_refiner.py + tests/test_hp_drift_guard.py = 46 passed, 17 skipped.
  pytest -m "not slow" (full) = 1108 passed then STOPPED (-x) at the SAME unrelated missing-fixture ERROR the
    prior audit flagged: tests/test_restore_polya_from_parquet.py::test_restore_cat3_plus_2 — a
    fixture-setup `assert _METADATA_TSV.exists()` on scripts/validation_data/rebuild_2026_05/trimmed/
    validation_reads_polya_trim_metadata.tsv which is ABSENT in this worktree. That test references
    junction_refiner/microhom/cap 0 times and is NOT in the cap commit range (git-confirmed) => missing-DATA
    condition, orthogonal to the cap. Clean re-run deselecting it: [see pytest_clean.txt] TRUE_RC + count below.
