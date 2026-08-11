# MICROHOM-DRIFT GUARD + NEAR-TIE CAP — ADVERSARIAL AUDIT V2 (DISCOVERY-LOSS LENS)

Auditor: adversarial V2 (discovery-loss lens, the LOAD-BEARING leg). Model: Opus 4.8.
Branch: worktree-agent-a25a2c1e784ad37dc. READ-ONLY (no edits/commits/default flips).
Resumes the prior auditor who STALLED at CHECKPOINT 3 (`dev/MICROHOM_AUDIT_discovery-loss.md`).

## THE FIX UNDER AUDIT (Phase 4, commit 05664bc on this branch)
`_effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap)`:
  = max(hold_margin, min(eff_margin, drift_near_tie_cap))   when cap active (>0)
  = eff_margin                                              when cap==0 (disabled, byte-identical)
Threaded through refine_read_junctions / _run_sequential / _run_parallel / refine_bam_junctions.
CLAIM: a move with delta_improve >= cap is NEVER drift-vetoed (cap BOUNDS discovery loss).
COUNTER-CLAIM to test: cap adds NO discriminating signal inside (0,cap) — real cryptics and
fab-drift OVERLAP there, so the (0,cap) overlap is IRREDUCIBLE with an aggregate-delta signal.

## PRIOR AUDITOR'S ESTABLISHED FACTS (build on, do not re-derive)
- Pure mh tie (delta_improve=0): guard-IRRELEVANT. Scorer soft-clips the drift-distance PREFIX
  (t1 free-prefix soft-clip: rescue[k:] vs genome[acc:acc+buf]), so at incumbent C the first
  `drift` bases are clipped and rescue[drift:] matches genome[C:] perfectly -> delta=0 -> refiner
  tie-break holds incumbent REGARDLESS of guard. Confirmed by instrumentation.
- Genuine guard-CAUSED loss REQUIRES: cryptic scores STRICTLY BETTER (delta_improve>0, discovered
  guard-OFF) yet still trips mh>=0.5 and is vetoed because delta_improve < eff_margin.
- The distinguishing evidence must live where the incumbent soft-clip CANNOT absorb it.

## MY PLAN (V2)
1. Read guard + _effective_veto_margin + scorer + _move_microhomology CODE (verify summary).
2. Read the fab benchmark harness (spikein_fab.py, spikein_score.py, gen_reads.py, build_panel.py,
   mh_characterize.py) + regression tests to understand truth model + read builder polarity.
3. BUILD the decisive panel the prior auditor named: REAL cryptic (3'SS/5'SS) reads whose TRUE site
   sits next to a SHORT direct repeat (mh>=0.5), incumbent N-op at the WRONG canonical site, read
   query = cryptic exon2 -> read WANTS to move. Sweep delta_improve ~0.5..8+ by controlled in-window
   base differences the C-soft-clip cannot absorb.
4. Measure DISCOVERY-LOSS RATE = f(microhom_drift_margin in {3,8}) x (cap in {off,1,2,3}) x
   (microhomology_frac). Same sweep: re-measure fab suppression (spike-in fab). Per-read truth,
   seed excluded.
5. VERDICT: does cap@margin=3 BOUND discovery loss acceptably? Is there a (margin,cap) with fab~0
   AND ~zero discovery loss, or is the (0,cap) overlap IRREDUCIBLE (=> cap insufficient, a per-base
   positional-distinctiveness signal the soft-clip can't absorb is required)? Numbers, not a hunch.
   Recommend a concrete (margin,cap) IF defensible, else state the cap is insufficient.

## CHECKPOINTS (appended as I go)

---
### CHECKPOINT 1 (V2) — orientation complete, code + harness verified independently
Branch head 05664bc = "near-tie read-evidence cap". Prior byte-identity auditor CLEARED off=inert;
synthesis verdict = HOLD (audit gate unsatisfied; 2/3 auditors stalled). I resume discovery-loss.

CODE VERIFIED (junction_refiner.py):
- Veto (line 838): `if best_score_cmp > incumbent_score - veto_margin: moves=False`.
  score_cmp LOWER=better. delta_improve := incumbent_score - best_score_cmp (>=0 for a move that
  the sort picked as best). Move SURVIVES iff best_score_cmp <= incumbent_score - veto_margin
  i.e. delta_improve >= veto_margin.
- eff_margin (802-820): hold_margin (+hp_drift if into_hp) (+microhom_drift if mh>=threshold).
- cap (_effective_veto_margin, 535-558): veto_margin = max(hold_margin, min(eff_margin, cap)) when
  cap>0 and eff_margin>hold_margin, else eff_margin. => with hold=0, veto_margin=min(eff,cap).
  STRUCTURAL CLAIM CONFIRMED CORRECT: any move with delta_improve>=cap survives (veto_margin<=cap).
  The code's OWN docstring (549-554) concedes the cap "does not add a discriminating signal inside
  the (0,cap) band, where real cryptics and fabricated drift still overlap." <- THIS is exactly the
  claim I must test empirically: is that overlap real & irreducible?

SCORER (junction_scoring.py:983 _score_junction): t1(k)=rescue[k:] left-anchored to genome[je:je+buf]
with FREE-PREFIX soft-clip (first k rescue bases dropped free); t2(k)=rescue[:k] rev right-anchored
to intron_end. score=min_k t1+t2. _MAX_RESCUE=30. Fast path _USE_CONCAT_DP when penalty_table=None
(flat ins) — my panel will use penalty_table=None (default) so evidence units are ~edit-distance.
=> the SOFT-CLIP is why a pure-mh drift ties (delta=0): at incumbent, first `drift` rescue bases
soft-clipped, remainder matches. Distinguishing evidence must be OUTSIDE the soft-clipped prefix.

FAB HARNESS (spikein_fab.py): builds CANONICAL truth + non-canonical drift target @+k, with
microhom_mismatch knob (frac of the k drift-flank bases that differ, exon2[0:k] vs exon2[k:2k]).
Reads spliced at TRUE canonical site; a drifted call = fabrication FP. This is the MIRROR of my
discovery panel. Its microhom_mismatch is the SAME evidence axis as my delta_improve.

MH_CHARACTERIZE claim (the n=2 under-powering to attack): fab-drift mh~0.66 vs R3-real mh~0.33,
"a threshold separates them." BUT R3 panel exon2 = random.Random => real cryptic sits at a random
transition (mh~0.25). The audit hinges on: do REAL cryptics ever sit next to a short direct repeat
(mh>=0.5) yet stay genuinely discoverable (delta_improve>0)? If yes, and if fab-drift also lives in
the same (delta_improve, mh) cell, the cap cannot separate them.

NEXT: build the decisive READ-CRYPTIC panel (own read builder, query=cryptic exon2, incumbent N-op
at canonical site). Sweep delta_improve by controlled in-window bases OUTSIDE the soft-clip prefix.
