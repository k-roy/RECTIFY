# MICROHOM-DRIFT GUARD + NEAR-TIE CAP — ADVERSARIAL AUDIT V4 (DISCOVERY-LOSS LENS, AUDITOR A)

Auditor: A (independent of B). Task: discovery-loss. Model: Opus 4.8.
Branch: worktree-agent-a25a2c1e784ad37dc. READ-ONLY (no repo edits / commits / default flips).
Independent redundant audit — my own harness, my own verdict. Not resuming B's harness.

## THE FIX UNDER AUDIT (commit 05664bc)
`_effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap)`
  = max(hold_margin, min(eff_margin, cap))   when cap>0 AND eff_margin>hold_margin
  = eff_margin                               otherwise (byte-identical to pre-cap)
Threaded through refine_read_junctions / _run_sequential / _run_parallel / refine_bam_junctions.
CLAIM: a move with delta_improve >= cap is NEVER drift-vetoed (cap BOUNDS discovery loss).
COUNTER (to test empirically): cap adds NO discriminating signal inside (0,cap) — real cryptics and
fab-drift OVERLAP there. If they share a delta_improve cell in (0,cap), cap is provably insufficient.

## GUARD MECHANICS (verified in junction_refiner.py — CP1 below)
- Veto (line 838): `if best_score_cmp > incumbent_score - veto_margin: moves=False`.
  score_cmp LOWER=better. delta_improve := incumbent_score - best_score_cmp (>=0 for chosen move).
  Move SURVIVES iff delta_improve >= veto_margin.
- eff_margin (802-820): hold_margin (+hp_drift if into_hp) (+microhom_drift if mh>=microhom_threshold).
- cap (_effective_veto_margin 535-558): with hold=0 -> veto_margin = min(eff_margin, cap).
  => any move with delta_improve>=cap survives. STRUCTURAL claim confirmed. Question is EMPIRICAL magnitude.
- SCORER (_score_junction): t1(k)=rescue[k:] left-anchored to genome[je:je+buf] with FREE-PREFIX
  soft-clip; t2(k)=rescue[:k] rev right-anchored to intron_end. score=min_k t1+t2. penalty_table=None
  => flat ins, delta_improve in ~edit-distance units matching cap/margin {1,2,3,8}.

## METHOD (advisor-endorsed direct-drive)
Drive refine_read_junctions DIRECTLY with synthetic pysam AlignedSegments — faster, deterministic,
per-read controllable; guard+veto+scorer all live inside that fn so it is faithful. NO pbsim/BAM.

DISCOVERY-LOSS defined OPERATIONALLY (not by construction):
  discoverable := refine_read_junctions(microhom_drift_margin=0) returns canonical->cryptic replacement.
  loss        := discoverable AND held when guard ON (margin>0, possibly capped).
This subsumes candidate-pool + boundary-filter constraints (a read guard-OFF doesn't move is out of
the guard's domain, not counted). boundary_error_window kept at production default (10).

PER-READ MEASUREMENT (advisor's anti-false-CLEAR mandate — do NOT assume from construction):
  * delta_improve: call _score_junction on incumbent (ns,ne) and cryptic (js,je); incumbent canonical
    => tier_beats_alt=False => raw scores, no _CANONICAL_HP_PRIOR confound => clean subtraction.
  * mh: call _move_microhomology(genome, ns, ne, js, je); confirm >=0.5 (else guard never fires, trivial).
A "no loss" cell is only believable if BOTH mechanisms verified to fire.

TWO TRAPS that produce a FALSE CLEAR (guard against):
  1. delta_improve actually 0 (soft-clip absorbed distinguishing bases). Distinguishing evidence must be
     in exon2 body BEYOND the microhomology repeat, surrounding genome non-repetitive (no alt soft-clip
     offset rescues incumbent). Use _no_run + random fill.
  2. mh never reaches 0.5 -> guard never fires -> trivial "no loss".

GRID: microhom_drift_margin in {3,8} x drift_near_tie_cap in {off,1,2,3} x microhomology_frac (sweep).
  Per (margin,cap): discovery-loss RATE + fab FDR. Populate (0,3) delta_improve band DENSELY.

DECISIVE OUTPUT: delta_improve DISTRIBUTION real-cryptics vs fab-drift on IDENTICAL harness/units;
  do they OVERLAP inside (0,cap)? Fab run through the SAME direct-drive builder (read from canonical,
  drift target in pool, does refiner wrongly move) — NOT the BAM pipeline, else distributions not comparable.

## CHECKPOINTS (appended as I go)

### CP1 — orientation + code verified (independent read)
junction_refiner.py verified as summarized above (lines 462-558 detectors + _effective_veto_margin;
730-841 scoring loop + eff_margin + veto). junction_scoring.py _score_junction (983-1168): confirmed
free-prefix soft-clip t1 + reversed t2, penalty_table=None fast path (_USE_CONCAT_DP) byte-identical.
_candidates_near (703-759): candidates come from a junction INDEX -> BOTH canonical incumbent AND cryptic
target must be seeded into the pool within max_boundary_shift(50). Fab harness pattern (spikein_fab.py)
studied for read-builder + microhomology installation (exon2[0:k] := near-copy of exon2[k:2k], n_mis flips).
NEXT: write harness (build reads, drive refiner, measure delta+mh per read), persist panel to disk.

### CP2 — harness written + FIRST decisive read confirmed
Harness: scratchpad/audit_v4/discovery-loss_A/harness.py (persisted). Direct-drive
refine_read_junctions with synthetic pysam reads; per-read MEASURED delta_improve + mh.
FIRST decisive read (k=6, mh_frac=0.8): delta_improve=6.0 (inc=6.0, cryptic=0.0 exact match),
mh_measured=0.667 (>=0.5 => guard fires), discoverable=TRUE (guard OFF => refiner moves
canonical->cryptic). => a REAL cryptic, genuinely discoverable, that TRIPS mh>=0.5. This is exactly
the delta>0-yet-mh>=0.5 case V2 named but never built. At delta=6: cap<=3 lets it survive (6>3);
margin=8 uncapped VETOES it (6<8). Need to sweep n_distinct/mh_frac/k to populate the (0,3) band
where cap decides + fab overlaps. NEXT: verify guard-ON veto fires as predicted; build full sweep.
