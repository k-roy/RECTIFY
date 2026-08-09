# MICROHOM AUDIT V5 — Task: strategic-realdata — Auditor A

**Role:** Adversarial auditor (Opus-Max) vetting RECTIFY's POSITIONAL-DISTINCTIVENESS CLOSE for the
motif-blind junction re-placer's microhomology-drift guard.
**Task:** STRATEGIC / REAL-DATA VALIDITY — "should this ship" review.
**Auditor:** A (independent; no coordination with B).
**Started:** 2026-07-14
**Mode:** READ-ONLY on branch worktree-agent-a25a2c1e784ad37dc.

## VERDICT: HOLD (fault_found=TRUE on the ENABLING claim; code safe default-off)
The POSITIONAL CLOSE's VALIDITY as an enabling claim ("the read-blind microhomology fault is CLOSED")
is **NOT established** — it is SYNTHETIC-ONLY INSURANCE whose worth is UNPROVEN. Two separable things:
(a) is the code safe to carry default-OFF? YES — verified unwired/byte-identical; NOT a fault.
(b) is "the fault is CLOSED, safe to ENABLE" validated? NO. HOLD on (b).

### ***Q1 CORRECTED — END-TO-END through the real refiner (this supersedes the block below)***
SHIPPED close = m3/cap2/GATE1 (drift_positional_gate=1, panel claim 0.4%/4.3%). GATE2 = the panel's
ALTERNATIVE "balanced knob" (panel claim 6.2%/2.1%). E2E LOSS = end-to-end discovery loss through the
REAL refiner, among cryptics discoverable guard-off (directly comparable to the panel's 0.4%). n=400,
seeds 1&2 agree:

  d    mh    Dintr   SHIPPED gate1 loss    ALT gate2 loss    regime
  0   1.00    0      100% (n~12-15)        100%              IRREDUCIBLE (identical windows; no method wins)
  1   0.97    2      1.6-2.4%              44-47% ***        PARALOG PSV regime (SMN1/SMN2 ~1 PSV/1-2kb)
  2   0.93    3      0.8%                  7%
  3   0.89    5      0-0.5%               1-2%
  5   0.82    7      0%                    0-0.5%            <- PANEL at-risk cryptics live here (ndisc med=7)
  10  0.65   12      0%                    0%
  20  0.29   15      0%                    0%                PANEL RANDOM-TAIL regime (mh<0.5, guard OFF)

HONEST HEADLINE: (i) at the SHIPPED gate=1 the close degrades GRACEFULLY: ~2% loss at d=1, ~0% by d=3.
More robust than a first pass suggested, but still 4-6x the panel's 0.4% claim, concentrated in the
guard's TARGET low-d paralog regime the panel never tests. (ii) The REAL danger is the gate=2 "balanced
knob" the panel PRESENTS AS AN OPTION (6.2%/2.1% looks mild): at d=1 paralog loci it costs ~45% discovery
loss — a precision-prioritizing user picking gate=2 on the panel's numbers unknowingly sacrifices HALF the
discovery at exactly the interesting loci, because the panel's at-risk cryptics sit at d~5-7 not d=1-2.
(iii) CLEANEST EXHIBIT (no gate needed): the PANEL'S OWN at-risk cryptics (n=1418) have ndisc median=7,
only 3.2% at d<=1, 9.6% at d<=2 (esig median=4). A RANDOM TAIL always follows the microhomology inside W,
so d stays high => the panel undersamples its own target population (high-mh-AND-low-d paralogs).

### THE FAULT I FOUND (Q1 — the panel is NOT representative of the guard's target regime)
[CORRECTED — see "Q1 CORRECTED" block below; gate1 is the SHIPPED point, not gate2.]
The panel's 98-99% balanced-accuracy claim was measured in the WRONG REGIME. My d-sweep (real
`_positional_signal` from the codebase, periodic-K microhomology, ONT 6/3%, seeds 1&2 agree) shows the
signal per +/-d where d = distinguishing bases in W=28. [ISOLATED-signal spare rates, gate=2:]
  d=1 (paralog PSV regime): real-cryptic SPARE rate = 0.44  → the guard VETOES ~56% of real cryptics.
  d=0 (perfect microhom):   spare = 0.00 (irreducible — no method wins; reported separately, weaker).
  d>=5: clean (spare 1.00).  d=20 (PANEL regime): spare 1.00 BUT mh=0.30 → the guard doesn't even trip.
Measured on the PANEL'S OWN at-risk cryptics (n=1418): ndisc median=7, only 3.2% have d<=1, 9.6% <=2;
esig median=4. The panel sweeps mh up to 1.0 but a RANDOM TAIL always follows the microhomology inside
W, so its d stays high (median 7) — it structurally NEVER samples high-mh-AND-low-d, which is exactly
the real paralog regime (SMN1/SMN2 ~1 PSV / 1-2 kb; rDNA; tandem dup) the guard TARGETS. The
"close" (0.4% loss / 4.3% fab) is an artifact of testing where the signal is easy.

### Q2 — the COMPASS harness AS-RUN never exercises the microhom guard, and the fab number was never measured
- The deployed COMPASS `refined.bam` = **arm-E (`hp_drift_margin=3.0`) = the HP-drift guard, NOT the
  microhom guard** (RECALL doc line 58, 72; §3 caveat 175-176 "microhom guard NOT yet applied"). So the
  headline 32x recall (0.54%->17.46%) was measured with `microhom_drift_margin=0` — the guard under
  audit is ABSENT from every real-data number that exists.
- §4b — the ONLY analysis that would measure the real-data FABRICATION rate the guard is meant to
  suppress AND confirm the operating point — **NEVER COMPLETED** (OOM exit 137, then 30-min TIMEOUT;
  algorithmic O(142k × per-chrom-list)). The doc declares it "no longer blocking" for the HOLD, but that
  concedes: the real fabrication rate is UNMEASURED. So the gate is NOT measurable with the existing
  harness as-configured — it needs BOTH the §4b binned-index fix AND `microhom_drift_margin>0` enabled.
- WHAT COMPASS MUST SHOW TO JUSTIFY ENABLING: run §4b with the microhom guard ON at the proposed
  operating point (m3/cap2/gate1) vs OFF, on the SMA paralog loci (the spike-in doc's real drift set:
  SNRPN ~6bp, UBA1 ~26bp, PCBP2 ~28bp), and demonstrate (i) the real fabrication rate it removes and
  (ii) that it does NOT drop motif-agnostic-short-read-CONFIRMED real non-canonical recall. Until that
  exists, enabling is unjustified.
- The recall doc's line-176 hand-wave — "removed drift largely absent from SR truth, so guard-on recall
  preserved" — is EXACTLY what my low-d finding rebuts: if the guard vetoes real low-d cryptics that BBMap
  CAN see (and BBMap is motif-agnostic, so it can), guard-on recall COULD drop; that is precisely the thing
  §4b-with-guard-ON would measure and never did. The preservation claim is assumed, not shown.

### Q3 — net-positive for a DISCOVERY tool? NO, on the current evidence
The tool's proven value is DISCOVERY: the real-data 32x non-canonical recall win. The guard suppresses a
~1.3-2.4% (synthetic/spike-in) fabrication whose REAL rate is unmeasured. At the SHIPPED gate=1 the
discovery cost is modest (~2% at d=1 paralogs, ~0% elsewhere) — so shipping-as-code default-off is benign.
But the guard exposes a gate=2 "balanced knob" that, on the panel's own numbers (6.2%/2.1%), looks like a
mild precision/recall trade while ACTUALLY costing ~45% discovery at d=1 paralog loci. For a discovery tool,
adding 4 params (margin/threshold/cap/gate) and a foot-gun operating point to suppress an UNMEASURED real
fabrication rate is NOT net-positive until COMPASS shows the real fab-removed and confirms recall is held.
Marginal at best; the burden of proof is on enabling, and that proof (real-data, guard-on) does not exist.

### Q4 — steelman: the honest stop is "acceptor close done, default-OFF, gate on COMPASS"
STRONGLY endorsed. The acceptor-side positional close is a legitimate, well-engineered, byte-identical-off
piece of CODE. But its EFFICACY claim is synthetic-only and over-favored by construction. The honest,
defensible position: SHIP AS CODE, DEFAULT-OFF; STOP building (no donor term, no 5th param); make the SOLE
remaining gate a COMPASS real-data measurement that runs the microhom guard ON over the SMA/paralog set and
reports real fab-removed vs real-recall-lost. "Keep building" adds complexity to an unenabled guard whose
fundamental limit (single-read ambiguity at low d) my sweep shows no scalar signal can escape. Do not enable
on synthetic evidence.

### CONFIRMS prior audit rounds
Consistent with V1/V2/V4 HOLD and the byte-identity-cleared-but-declined-to-endorse pattern. My independent
contribution: a QUANTIFIED regime-representativeness fault (d-sweep) + the sharpened Q2 fact that the
COMPASS harness as-run never exercised the guard at all.

## Checkpoints
- [x] CP0: durable record created
- [x] CP1: code read (junction_refiner.py _positional_signal, wiring, defaults) — CONFIRMED unwired
- [x] CP2: discovery-loss claim read + spot-checked
- [x] CP3: COMPASS recall result read (is enablement gate measurable?) — §4b NEVER COMPLETED
- [x] CP4: guard design + handoff read
- [x] CP5: independent harness on synthetic-representativeness question — d-sweep DONE (fault found)
- [x] CP6: verdict — HOLD, fault_found=TRUE on enabling claim

## ARTIFACTS
- Harness: /private/tmp/claude-501/-Users-kevinroy-work-rectify/798e5203-2185-4c27-87b4-f4e50c56c3ba/scratchpad/audit_v5/strategic-realdata_A/dsweep.py
- Results: .../strategic-realdata_A/RESULTS.txt (d-sweep + panel-instrumentation, seeds 1&2)
- Reproduce: /Users/kevinroy/miniconda3/bin/python dsweep.py --n 400 --seed {1,2}

## CP1 — FACTUAL ANCHORS (verified by code read)
1. **Guard is TRULY UNWIRED in production.** `correct_command.py:746` calls `refine_bam_junctions`
   passing NO drift kwargs. grep for `drift_positional_gate|microhom_drift_margin|drift_near_tie_cap`
   across `rectify/` (excluding junction_refiner.py) returns EMPTY. All 4 params default 0.0 = OFF =
   byte-identical. So in the shipped tool, NONE of this code path ever executes. It is dead-at-default.
2. **_positional_signal (line 600-629): ACCEPTOR-ONLY by design.** Returns None if `new_je == ne`
   (acceptor didn't move). rescue = q[q_split:q_split+28]. Compares hard-anchored ed of rescue to
   incumbent exon2 window vs moved exon2 window. Signal>0 ⇒ read favors moved acceptor. Docstring
   explicitly: donor-side ed term would be net-harmful → acceptor-only is deliberate, not incomplete.
3. **_semiglobal_ed (line 580): standard DP, free ref SUFFIX, hard-anchored ref[0].** Correct as a
   "does the read's exon2 match this candidate's downstream genome better" measure. O(|q|·|ref|) but
   only on the rare would-be-veto path.
4. **Veto path (915-930):** only fires when `moves and eff_margin>0 and incumbent_score is not None`.
   eff_margin>0 requires a drift margin was ADDED (hp OR microhom flagged). So the whole positional
   machinery is gated behind a drift flag AND a would-be-veto. Narrow blast radius by construction.

## CP2 — DISCOVERY-LOSS CLAIM (the CLOSE's evidence) — read, structural read
- The claim (`DISCOVERY_LOSS_PANEL_RESULT.md`, CP4/CP5): WIRED m3/cap2/gate1 → 0.4% discovery loss /
  4.3% fab-residual on a SYNTHETIC panel. Balanced accuracy 98-99% in the delta overlap band.
- CONSTRUCTION (lines 11-18 of that doc + discovery_loss_panel.py): genome = LPAD + exon1 +
  canonical GT..AG intron + EXON2; EXON2 = U(k-mer) + U'(U with `mm` mismatches) + random TAIL.
  Microhomology = U~U'. Real cryptic read carries genome[je:] = U'+TAIL; placed at incumbent; refiner
  should MOVE ne->je. ONT error model 6% sub / 3% indel. fab = canonical read that ONT-error-drifts.
- **This is a SINGLE-INTRON, RANDOM-TAIL, RANDOM-SEQUENCE toy.** The discriminating power of the ed
  signal comes ENTIRELY from the random TAIL beyond the microhomology (U'+random). See CP5 for the
  adversarial break: real genomic loci where BOTH candidate exon2 windows are near-identical for the
  FULL 28bp window (paralog/tandem/rDNA) have NO distinguishing tail → ed signal → 0 → no separation.

## CP3 — COMPASS REAL-DATA GATE — the decisive strategic finding
- COMPASS_RECALL_RESULT.md §3: the motif-blind re-placer lifts non-canonical recall **~32x
  (0.54%->17.46%)** vs raw minimap2 against an INDEPENDENT motif-agnostic BBMap short-read truth,
  holding canonical flat (94.63%->94.41%). THIS IS THE TOOL'S CORE VALUE and it is REAL-DATA proven.
- **BUT §4b — the AGGREGATE REVERSE CHECK that classifies arm-B's revealed non-canonical junctions
  as recovery / drift(fabrication) / inconclusive — NEVER COMPLETED.** It OOM'd (exit 137), then
  TIMED OUT at 30min (algorithmic O(142k × per-chrom-list): billions of ops). The doc itself
  (§4b UPDATE, line 183-192) declares §4b "NO LONGER BLOCKING" for the HOLD verdict — but that means
  **the real-data FABRICATION RATE the guard is supposed to suppress HAS NEVER BEEN MEASURED.**
- The guard exists to suppress a fabrication rate. The synthetic estimate is ~1.3-2.4%. The spike-in
  estimate is 1.31%. **Neither is real-data. The one real-data measurement that would confirm both the
  fabrication rate AND the 0.5/margin operating point is the §4b that never ran to completion.**
- Note the recall table's own caveat (line 175-176): "microhom guard NOT yet applied (its removed
  drift is largely absent from SR truth anyway, so guard-on recall on this set ≈ preserved)." So even
  the 32x recall number was measured WITHOUT the microhom guard active — the guard's effect on real
  recall is inferred, not measured.

## CP4 — GUARD DESIGN + HANDOFF — the audit trail itself concedes the gaps
- HANDOFF_MICROHOM_GUARD.md is candid: guard COMMITTED, default OFF, "until the triple-audit is
  all-clear AND the COMPASS real-data threshold is confirmed." Both prior audit rounds (V1, V2, V4)
  returned HOLD; the load-bearing empirical legs repeatedly stalled.
- The design's own "REMAINING before enable" (line 90-92): (ii) RE-AUDIT the close; (iii) COMPASS
  real-data confirmation (independent). The authors themselves list COMPASS real-data as an unmet
  ship prerequisite. This is not a hidden flaw — it is a self-declared open gate.


---

## CP4b — ADVISOR-SHARPENED PLAN (2026-07-14)
- **Q1 decisive test = distinguishing-base-density (d) sweep**, NOT identical-windows (that's the
  degenerate null — d=0 = irreducible ambiguity, easy to rebut as "nothing lost"). `_positional_signal`
  = ed(rescue,inc) - ed(rescue,mov). For a true cryptic read: ed(mov)≈ONT_err, ed(inc)≈ONT_err + d
  where **d = # distinguishing bases between the two W=28 windows**. So signal ≈ d.
  - Panel: post-microhom tail is RANDOM → d≈20+ → signal huge → 98-99% acc. THE ONLY REGIME TESTED.
  - Real paralogs (SMN1/SMN2 ~1 PSV per 1-2kb), rDNA, tandem dup: **d=0 or 1 in a 28bp window is COMMON.**
  - At d=1, signal ≈ 1 ± noise ≈ gate=1 → real cryptics vetoed a large fraction. Evidence EXISTS but
    the signal can't resolve it under error. **THE GENUINE FAULT, on exactly the loci the guard targets.**
  - Guardrails: keep divergence high enough that refiner still attempts the move AND mh≥threshold trips;
    report d=0 (irreducible, weaker) separately from d=1-2 (the real fault). Expect collapse at d≤2.
- **Q2 sharpened:** recall doc lines 175-176 — the microhom guard was NOT applied when 32x recall
  measured (that's arm-E = HP-guard only). Combined with §4b never completing: **the COMPASS harness
  as-run NEVER EXERCISES the microhom guard at all.** Gate is NOT measurable with the existing harness
  as-configured — needs the §4b binned-index fix AND microhom_drift_margin>0 actually enabled.
- **Verdict framing — split two things:** (a) safe to carry default-off? YES (byte-identical/unwired,
  not a fault). (b) is "the fault is CLOSED" validated for ENABLING? NO (synthetic-only, panel
  over-favors signal, real fab rate never measured, gate harness never ran the guard). HOLD on (b).

## Working notes

