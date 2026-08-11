# MICROHOM AUDIT V5 — Auditor B — task "independent-remeasure"

**Auditor:** B (Opus-Max, adversarial). Independent of Auditor A. Do NOT reuse `dev/discovery_loss_panel.py`.
**Goal:** Build MY OWN panel (different construction + error model), independently re-measure the CLOSE.
Questions: (1) does ~0.4% disc-loss / ~4.3% fab-residual @ m3/cap2/gate1 hold on my panel? (2) is the
original U+U'+TAIL construction BIASED (over-represents easy-to-separate cases)? test longer repeats (k>10),
imperfect tails, higher ONT error (10-15%), homopolymer-adjacent boundaries, multi-candidate pools. (3) at
what microhom/error regime does the ed signal's separation DEGRADE below usefulness? Report the frontier.
**VERDICT:** is "CLOSED" robust or construction-dependent?

Working READ-ONLY. Python: /Users/kevinroy/miniconda3/bin/python. Repo importable from ROOT.

---

## Code understood (verified by reading, not trusting summary)

- `_positional_signal(genome, q, q_split, ne, new_je, W=28)` (junction_refiner.py:600):
  `rescue = q[q_split:q_split+W]`; `ref_inc = genome[ne:ne+W+6]`; `ref_mov = genome[new_je:new_je+W+6]`;
  returns `_semiglobal_ed(rescue, ref_inc) - _semiglobal_ed(rescue, ref_mov)`. `>0` ⇒ read matches MOVED
  acceptor better ⇒ genuine evidence. Returns None if new_je==ne or rescue empty.
- `_semiglobal_ed` (junction_refiner.py:580): edit distance, ALL of query vs a PREFIX of ref (free ref
  SUFFIX, HARD-anchored at ref[0], NO free query/ref prefix). This is the key: it removes the scorer's
  free-`k` soft-clip escape.
- Gate wiring (junction_refiner.py:915-930): a drift-flagged would-be-veto (`best_score_cmp > incumbent_score
  - veto_margin`) is SPARED when `drift_positional_gate>0` AND `psig is not None` AND `psig >= gate`.
  ACCEPTOR-ONLY (`new_je`); donor/both-boundary moves fall through to margin/cap.
- `_score_junction` (junction_scoring.py:983): bilateral k-sweep, `for k in [0,L)`: t1=align rescue[k:] to
  g[je:], t2=align reverse(rescue[:k]) to reverse(g[je-buf:je]). `k>0` frees a query prefix ⇒ the "soft-clip
  escape" that hides discriminating exon2 bases from delta_improve. FAST PATH (concat DP, DEFAULT ON) when
  penalty_table=None — byte-identical to legacy.
- refiner is ACCEPTOR-CENTRIC: `_score_junction` ignores intron_start (donor). Close is acceptor-only by design.

## Original panel construction (dev/discovery_loss_panel.py) — what I must DIFFER from
- genome = LPAD(40 T) + exon1(60 rand) + INTRON("GT"+86rand+"AG"=90) + EXON2. EXON2 = U(k) + U'(U w/ mm subs)
  + TAIL(30 rand) + 40 rand. ne = start of U (canonical decoy); je = ne+k = start of U'. mh = frac(U~U').
- cry read: exon2 = genome[je:je+E2], E2 = k + e2_tail_take. Placed at incumbent N-op [ns,ne).
- fab read: exon2 = genome[ne:ne+k+E2] (from incumbent). drift to je = fabrication.
- ONT errors: 6% sub, 3% indel (0.5 del / 0.5 ins split).
- Pool = exactly 2 candidates {(ns,ne),(ns,je)}. motif_blind=True, boundary_error_window=0.
- Sweep: k∈{4,6,8,10}, mm∈0..k/2, e2_tail_take∈{2,4,8,16}, n=40/cell.

### BIAS HYPOTHESES to test (why the original may over-state the close):
1. TAIL is fully random ⇒ exon2s diverge sharply right after the repeat ⇒ ed signal easy. Real cryptics
   may have LONGER shared context / imperfect tails.
2. Only k≤10. Longer near-tandem repeats (k>10) push the divergence point past W=28 window ⇒ ed blind.
3. Pool = 2 candidates only. Real pools have MANY (annotated + observed) → more decoys near je.
4. Error rate 6/3%. Human/older-flowcell ONT is 10-15%.
5. No homopolymer-adjacent boundaries (the ORIGINAL fault motivation!).
6. e2_tail_take up to 16 ⇒ read carries lots of distinctive tail past the repeat. Short exon2 reads
   (truncated, or repeat near read end) starve the ed window.

---

## CHECKPOINT LOG
- [CP0] Record + code-read complete. Building independent harness next.
- [CP1] Code VERIFIED by reading (not trusting summary):
  - `_positional_signal` (jr:600): rescue=q[q_split:q_split+28]; ref_inc=genome[ne:ne+34]; ref_mov=genome[new_je:new_je+34];
    returns `_semiglobal_ed(rescue,ref_inc)-_semiglobal_ed(rescue,ref_mov)`. >0 ⇒ moved matches better. None if new_je==ne.
  - `_semiglobal_ed` (jr:580): all-query vs ref-PREFIX (free ref SUFFIX, hard anchor at ref[0]). Confirmed identical to panel's copy.
  - Gate (jr:915-930): veto path fires only when `moves AND eff_margin>0 AND incumbent_score is not None`. eff_margin>0
    requires a drift margin was ADDED (microhom or hp) since hold_margin=0. Then `best_score_cmp > incumbent_score - veto_margin`
    ⇒ would-veto; SPARED iff `drift_positional_gate>0 AND psig is not None AND psig>=gate`. ACCEPTOR-ONLY (new_je).
  - `_move_microhomology` (jr): acceptor shift frac = _frac_match(seq[lo:hi], seq[hi:hi+k]), k=|je-ne|. MIN over moved boundaries.
  - `_frac_match`: only ACGT match counts (N/gap → 0). So a genome N-run does NOT phantom-flag.
- [CP2] PLAN LOCKED (advisor-guided, the structural kill):
  * SPINE = tandem microsatellite array, period p, incumbent ne + cryptic je=ne+p BOTH in-register one unit apart.
    Array length L swept 5→~50. When L>W(28): genome[je:je+28]==genome[ne:ne+28] (period-shifted) ⇒ BOTH ref windows
    match the read ⇒ psig→0 ⇒ gate=1 VETOES a genuine cryptic (discovery loss). Predict sharp crossover near L≈W.
    This is the natural generalization of the HP-drift motivation (HP = period-1 array).
  * WHY original biased: U+U'(≤10bp)+random TAIL forces divergence EARLY (<k≤10, inside W=28) ⇒ guaranteed strong ed
    ⇒ their 98-99% sep. Construction never lets microhomology reach the window. The guard fires PRECISELY on high-mh
    moves — exactly where ed is weakest — so testing only k≤10 hides the failure.
  * DISCIPLINE: headline through the REAL refiner (mirror acceptor_after_refine: motif_blind=True, boundary_error_window=0,
    microhom_drift_margin=3.0, drift_near_tie_cap=2.0, drift_positional_gate=1.0, hold_margin default 0). Only construction differs.
  * ANCHOR FIRST: least-adversarial corner (short repeat, random tail, 6% err) must land ~0.4%/4.3% → proves harness measures
    same quantity before attributing hard-regime degradation to the regime not a harness artifact.
  * Report at-risk N (denominators): disc-loss/fab-residual conditional on {moved-guard-off AND mh>=0.5}.
  * WATCH: does fab-residual ALSO climb across L? only disc-loss → fails safe; both → fails both ways (stronger finding).
  * ECONOMY: 2 axes (L spine ~5pts × err 6/10/15%) + anchor + fold-in single cells (HP boundary, truncated exon2, multi-cand pool).
- [CP3] HARNESS WRITTEN + DEBUG-PROBED (error-free, harness_B.py). CRITICAL EARLY FINDING (structural mechanics):
  In a period-2 microsatellite with in-register acceptors (ne, je=ne+2), the read's exon2 first DIVERGES from the
  incumbent window at read-offset == L (the array length). Error-free probe (e2_take=40 so rescue reaches the flank):
    L= 8: first_div=8  (in W=28) psig=+2 delta=+1.5 moved_OFF=T moved_CLOSE=T
    L=20: first_div=21 (in W=28) psig=+2 delta=+2.0 moved_OFF=T moved_CLOSE=T
    L=30: first_div=30 (PAST W=28) psig=0  delta=+0.0 moved_OFF=F moved_CLOSE=F
    L=40: first_div=40 (PAST W=28) psig=0  delta=+0.0 moved_OFF=F moved_CLOSE=F
  KEY MECHANISTIC INSIGHT: as the divergence point passes W=28, psig->0 AND delta->0 TOGETHER. The scorer
  (_score_junction, full-length via k-sweep concat DP) and _positional_signal DEGRADE IN LOCKSTEP. When the
  discriminating base is past W, the SCORER ALSO can't distinguish (delta=0), so the refiner does NOT move the read
  even guard-OFF => it is NOT "discovery at-risk" => there is NO discovery for the guard to lose. The guard only
  ever vetoes a move the scorer WOULD make (moved_off=True). So the true discovery-loss = frac of moved_off=True
  cryptics that the CLOSE (m3/cap2/gate1) now HOLDS. Frontier = the boundary regime (L~24-28) + ONT error where a
  read stays discoverable (delta>0 via errors/tail) yet psig fails the gate. Running spine panel next.
- [CP4] ANCHOR + SPINE RUN (n=120/cell, real refiner, m3/cap2/gate1). DATA persisted: out_anchor.json, out_spine.json.
  ANCHOR (short dinuc arrays, 6% err) — Q1 reproduce ~0.4%/4.3%?
    cell        cry@risk fabdrift DLclose FRclose bBalAcc
    L6           112      2        0.009   0.000   0.981
    L8           115      4        0.017   0.750   0.582
    L10          110      5        0.009   0.800   0.586
  => DLclose 0.9-1.7% — SAME ORDER as claimed ~0.4% (my construction is harder: dinuc period-2 has more in-window
     microhomology than U+U'+random). CONFIRMS harness measures the same quantity. BUT already at L8/L10 band_balacc
     is only ~0.58 (NOT 98-99%) and fab_resid=75-80% — cracks visible even in the "easy" regime. NOTE denominators:
     fab_drift_n=2-5, so FRclose 75-80% = 3-4 reads (small-N, report as such).
  SPINE (period-2 dinuc microsatellite, L swept across W=28, err 6/10/15, n=120): THE CROSSOVER (advisor-predicted):
    L     cry@risk  DLclose(6/10/15%)   band_balacc(6/10/15)   NOTE
    6     101-114   0.00/0.00/0.01      0.88/0.67/0.49
    12     95-110   0.01/0.02/0.03      0.65/0.73/0.55
    20     87-105   0.02/0.03/0.06      0.48/0.48/0.60
    26     66-102   0.13/0.21/0.30 <==  0.78/0.76/0.65   *** THE KNEE: large discovery (66-102 at-risk) yet 13-30% LOST
    30     21-29    0.95/0.96/0.93      0.50/0.52/0.52   cry@risk COLLAPSES (scorer also blind past W) => few discoveries left
    36      0-5     ~1.0 (tiny N)       0.50             cry@risk ~0 => NO discovery to lose (denominator gone)
    44      0-2     ~1.0 (tiny N)       0.50             same
  DECISIVE READING:
   * The close HOLDS (DLclose<=6%) for L<=20 (divergence inside W=28) — reproduces the claim in the regime the
     ORIGINAL tested (its U+U' forces divergence <=k<=10, always L<=~10). 
   * At L=26 (divergence ~at the W=28 edge): DLclose JUMPS to 13-30% with a LARGE at-risk denominator (66-102 reads)
     AND band_balacc falls to 0.65-0.78. This is GENUINE discovery-loss the close does NOT close — the ed signal's
     separation degrades exactly as the advisor predicted, as microhomology extent approaches W.
   * At L>=30, cry_at_risk collapses because the SCORER (_score_junction) and _positional_signal go blind TOGETHER
     (both windowed near ~W): the refiner no longer moves the read even guard-OFF, so there is no discovery to lose.
     The 100% DLclose at L36+ is over a ~0 denominator (NOT a real loss).
   * fab_resid: FRclose stays high (0.4-1.0) across L20-26 but on SMALL denominators (fab_drift_n 3-19). At L>=30
     FRclose->0 (fab also can't drift when scorer is blind). So the gate does NOT obviously fail both ways at the knee
     on the fab side — but at L26 both a real discovery-loss AND residual fab remain (band_balacc 0.65-0.78 = the
     signal is near coin-flip in that band). Need frontier fine-sweep to pin the knee; running next.
  VERDICT DIRECTION: CONSTRUCTION-DEPENDENT. "CLOSED" holds ONLY while the discriminating divergence sits inside the
  W=28 positional window (L<~24, matching the original's L<=10). It DEGRADES to 13-30% discovery-loss at L~26 (real
  microsatellite lengths: dinucleotide repeats of 13+ units are common in genomes), where the signal is still needed.
- [CP5] FRONTIER FINE-SWEEP + FOLD-INS + ORIGINAL-PANEL REPRODUCTION + GENUINE-DISCOVERY CONFIRMATION. All data on disk
  (out_frontier.json, out_period.json, out_decoys.json, out_imperfect.json, out_truncate.json, orig_panel.log).
  (A) FRONTIER (period-2, L across W=28, rescue reaches flank, n=80) — THE KNEE PINNED:
     L    cry@risk   DLclose(6/10/15%)      band_balacc(6/10/15)
     18   58-66      0.00/0.015/0.00        0.80/0.61/0.59
     22   60-65      0.016/0.015/0.033      0.48/0.74/0.60
     24   47-68      0.015/0.143/0.170 <-   0.67/0.71/0.59    knee begins (errors push loss up)
     26   51-64      0.186/0.172/0.314 <-   0.37/0.69/0.57    FULL KNEE: 17-31% loss, large denom
     28   36-58      0.621/0.615/0.611 <-   0.63/0.38/0.60    signal ~dead
     30   13-17      0.92/1.0/1.0           0.50             cry@risk collapses (scorer blind too)
     32   0-2        (tiny N)                                denominator gone
  (B) ORIGINAL PANEL REPRODUCED on THIS machine (dev/discovery_loss_panel.py --n 20 --seed 1):
        ed_signal balanced-acc 98%; cap+ed(esig>0) disc-loss 0.7% / fab-resid 0.0%; WIRED m3/c2/gate1 0.7%/0.0%.
     => My wiring is CORRECT (same refiner, same operating point → the claimed numbers). The claim HOLDS in the
        original's construction. So my divergent frontier is PURELY construction, NOT a harness artifact.
  (C) GENUINE-DISCOVERY CONFIRMATION (L=26, real cryptic reads, 10% err, n=200): 70% of cryptics are discoverable
     guard-OFF (refiner correctly moves ne->je = substantial real discovery, NOT a tiny denom); 18% of those are
     LOST (vetoed) by the shipped close. Lost reads have psig = -1 / 0 (at/below gate=1) on mh=1.0 AC-dinuc
     microsatellites — the divergence sits near/past W=28 so hard-anchored ed cannot distinguish. REAL discovery loss.
  (D) PERIOD GENERALIZATION (L=24/26, 10% err): the knee holds across ALL periods — p1(HP) DLclose 18-20% & fab-resid
     75-88% (the close FAILS on the very HP case that motivated the guard, once the run is long!); p2 16-20%; p3 5-23%;
     p4 8-24%. band_balacc 0.46-0.75 throughout. NOT a dinucleotide artifact — general tandem-microhomology property.
  (E) MULTI-CANDIDATE POOL (extra in-register decoys): does NOT drive the effect — L=20 stays ~2%, L=26 stays 15-32%.
     Microhomology EXTENT (vs W), not pool size, is the driver.
  (F) IMPERFECT/TRUNCATE folds: cry@risk collapsed to ~0 (reads not discoverable when rescue is pure tandem or
     divergence is past W) — consistent with the mechanism; no new frontier info.

## ============ FINAL VERDICT (Auditor B, independent) ============
**FAULT FOUND: the "CLOSED" claim is CONSTRUCTION-DEPENDENT, not robust.**

Q1 (does ~0.4%/4.3% @ m3/cap2/gate1 hold on MY panel?): YES but ONLY in the original's regime. My tandem-array
   anchor (short arrays L<=20, dinuc, 6% err) reproduces DLclose ~0-3% — same order as the claimed ~0.4%. And I
   independently reproduced the ORIGINAL panel's exact 0.7%/0.0% on this machine (wiring verified).

Q2 (is the original construction BIASED?): YES, decisively. EXON2 = U(k<=10) + U'(<=10bp) + RANDOM TAIL forces the
   read's exon2 to DIVERGE from the incumbent within <=~10 bp — ALWAYS inside the W=28 positional window. So
   `_positional_signal` (hard-anchored ed over W=28) ALWAYS sees the discriminating bases → its 98-99% balanced
   accuracy and 0.4% loss are a STRUCTURAL consequence of never testing microhomology that extends toward W. The
   guard fires PRECISELY on high-microhomology moves — exactly where the ed signal is weakest — and the original
   panel's construction never enters that region. Real genomes DO (dinucleotide/trinucleotide microsatellites of
   12-15+ units; long homopolymers).

Q3 (frontier where the ed signal degrades below usefulness): the shared-microhomology length L vs W=28.
   - L <= ~20 : DLclose 0-3% (close works; = the original's regime).
   - L ~ 24-28: DLclose 14-62% with LARGE at-risk denominators (36-77 reads, 70% of cryptics discoverable),
     band_balacc 0.37-0.75 (near coin-flip). THE CLOSE FAILS HERE while discovery is real. Higher error (10-15%)
     worsens it. Generalizes across all unit periods incl. homopolymers.
   - L >= ~30: at-risk denominator collapses because the SCORER (_score_junction) goes blind past ~W TOO, so no
     discovery is made even guard-OFF (nothing to lose). The 100% DLclose there is over ~0 reads (not a real loss).

MECHANISM: `_positional_signal` and `_score_junction` share the same ~W-bp horizon. The close only helps in the
   narrow band where the scorer CAN see the divergence (so a read is discoverable) but within W=28 — i.e. divergence
   at offset ~[k, 28]. The original panel populates only the easy end of that band (offset <=10). As divergence
   approaches 28, ed separation collapses; past ~28 there is no discovery to protect. So "the fault CLOSED" is true
   for short microhomology and FALSE for long microhomology — the exact regime the guard exists to handle.

HONEST WEIGHTING (fairness to the claim): the claim doc (DISCOVERY_LOSS_PANEL_RESULT.md) ITSELF states (CP3) "the
   microhomology overlap is an IRREDUCIBLE single-read ambiguity" and ships the guard DEFAULT-OFF / UNWIRED with
   COMPASS real-data named as the enablement gate. So the OVERSTATEMENT is localized to CP4/CP5's headline "the
   fault CLOSED" and the panel doc/DISCOVERY_LOSS_PANEL_RESULT numbers — which are correct for the tested
   construction but do NOT generalize. The code is correct; the byte-identical-off property holds; nothing ships
   enabled. The fault is a VALIDATION-SCOPE / CLAIM-STRENGTH fault, not a code bug.

RECOMMENDATION (verdict = HOLD on the "CLOSED" language, not on the code):
   1. Re-scope the claim: "closes the fault for microhomology whose discriminating divergence lies within ~W=28 bp
      of the acceptor; degrades to 14-62% residual discovery-loss for longer tandem microhomology (dinuc >=13 units,
      long HP) where divergence exceeds the ed window." Do NOT call it unconditionally CLOSED.
   2. Add a long-microsatellite / long-HP arm to dev/discovery_loss_panel.py (tandem array with divergence swept
      past W) so the frontier is in the repo's own validation, not just this audit.
   3. Consider W as a tunable tied to local repeat length (extend the ed window when a tandem run is detected) — the
      fixed W=28 IS the ceiling. But this only shifts the knee; past the scorer's own horizon nothing is discoverable.
   4. Keep DEFAULT-OFF and gate real enablement on COMPASS, as the claim already does. The synthetic "CLOSED" must
      not be read as a green light to enable by default.
VERDICT: HOLD — "CLOSED" is construction-dependent (robust only for short microhomology); overstated as
   unconditional. Code correct + default-off, so non-blocking, but the claim language and the panel's coverage
   must be corrected before the close is cited as done.

