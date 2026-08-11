# MICROHOM AUDIT V5 — task "strategic-realdata" — Auditor B

**Role:** Adversarial auditor (Opus-Max), INDEPENDENT of Auditor A. Read-only.
**Task:** STRATEGIC / REAL-DATA VALIDITY — "should this ship" review of the POSITIONAL CLOSE
for the motif-blind junction re-placer's microhomology-drift guard.
**Verdict field:** CLEAR (validity established) vs HOLD (synthetic-only insurance, worth unproven).

Started: 2026-07-14 (session date). ATTEMPT 1/2.

---

## Checkpoint log (append-only)

### CP0 — record initialized
Skeleton written. Scratch dir created. No prior record/harness found — fresh start.

Plan:
1. Read the CLOSE code (junction_refiner.py: _positional_signal, _semiglobal_ed, drift_positional_gate wiring).
2. Read the claimed-result docs (DISCOVERY_LOSS_PANEL_RESULT.md, discovery_loss_panel.py).
3. Read the COMPASS harness/result (COMPASS_RECALL_RESULT.md) — the enablement gate.
4. Verify default-OFF / unwired claim (correct_command.py).
5. Adversarial questions:
   (1) synthetic panel representativeness vs real ONT cryptic-in-microhomology loci
   (2) what COMPASS must show + is it measurable with existing harness
   (3) net-positive for a discovery tool? (suppress ~2% fab vs discovery risk + 4 params)
   (4) steelman "honest stop" vs "keep building"
6. Verdict.

### CP1 — code + result docs read
**CLOSE code (junction_refiner.py):**
- `_semiglobal_ed(query, ref)` (L580): edit distance aligning ALL of query to a PREFIX of ref (free ref
  suffix, HARD-anchored at ref[0], NO free query/ref prefix). Standard DP, O(m·n). Verified correct.
- `_positional_signal(genome, q, q_split, ne, new_je, W=28)` (L600): rescue = q[q_split:q_split+W];
  ref_inc = genome[ne:ne+W+6], ref_mov = genome[new_je:new_je+W+6];
  returns ed(rescue, ref_inc) - ed(rescue, ref_mov). >0 ⇒ rescue matches MOVED acceptor better.
  Returns None if new_je==ne (acceptor didn't move) or rescue empty. ACCEPTOR-ONLY by design.
- Wiring (L916-923): veto_margin = _effective_veto_margin(...); if delta_improve < veto_margin
  (would-be-veto) AND drift_positional_gate>0: psig = _positional_signal(...); if psig>=gate → SPARE
  the veto. Default gate=0.0 → OFF → byte-identical (the >0.0 guard).

**Claimed result (DISCOVERY_LOSS_PANEL_RESULT.md):** CP4 indel-robust ed_signal 98% balanced-acc,
cap+ed(>0) → 0.4% loss/4.3% fab; CP5 WIRED end-to-end 0.4%/4.3% confirmed, STILL default OFF.
NOTE: pre-CP4 VERDICT block (L133-144) still says "does NOT cleanly close... IRREDUCIBLE... deep scorer
surgery buys a tradeoff-shift not a clean close, hard to justify". CP4/CP5 claim to SUPERSEDE it. Tension.

**COMPASS_RECALL_RESULT.md:** §1-3 DURABLE: 32× recall win (raw-mm2 0.54% → arm-B 17.46% vs independent
motif-agnostic BBMap short-read truth; canonical held flat 94.63→94.41%). §4b (the REAL-DATA FABRICATION
rate — recovery/drift/inconclusive of arm-B's recurrent revealed non-canonical junctions) **NEVER
COMPUTED**: OOM then TIMEOUT; binned-index fix specified but NOT run; declared "context, not a gate."
⚠ The guard's TARGET metric (real-data fabrication) was NEVER quantified on COMPASS. And the microhom
guard + positional close were NOT applied on COMPASS at all (RESULT.md L176).

NEXT: read discovery_loss_panel.py (error+microhomology model realism); verify default-OFF in
correct_command.py; assess synthetic representativeness vs real cryptic-in-microhomology loci.

### CP2 — DEFAULT-OFF + UNWIRED verified (stronger than the brief claimed)
Production caller `correct_command.py:746 refine_bam_junctions(...)` passes:
- NO `motif_blind` (defaults to **False**) → the ENTIRE motif-blind re-placer arm is OFF in production,
  not just the guard. Default (motif-aware) mode only slides to ANNOTATED junctions.
- NO drift kwargs → all default 0.0 → OFF.
Grep across all of rectify/ (non-test, non-junction_refiner.py): motif_blind set True = **0 hits**;
drift kwargs passed = **0 hits**; no CLI flag, no config key. ⇒ NO production enablement path exists.
The close guards a feature (motif_blind=True) that is itself DARK CODE reachable only from dev/test.
So "default-OFF byte-identical" is TRUE and STRONGER than stated: the whole arm is inert in shipped behavior.

**Synthetic panel construction (discovery_loss_panel.py) — key properties for representativeness:**
- Genome FULLY RANDOM i.i.d. uniform ACGT (exon1, intron interior, TAIL, pads all `_rand_seq`).
- EXACTLY ONE engineered microhomology: EXON2 = U(k)+U'=mutate(U,mm)+TAIL(random 30)+random(40).
  mh = 1 - mm/k, swept k∈{4,6,8,10}, mm∈0..k/2 → mh 1.0..0.5.
- ed_signal discriminates because rescue=U'+TAIL and PAST the k-bp microhomology the two candidate exon2s
  DIVERGE at ~0.75/base (random TAIL). The signal's power is MANUFACTURED by random-TAIL distinctiveness.
- ONT error = uniform 6% sub / 3% indel i.i.d. per base; NO homopolymer-context, NO systematic bias,
  NO position/quality correlation.
- The MECHANICS are real (_score_junction + refine_read_junctions); only the SEQUENCE CONTEXT is synthetic.

NEXT: reproduce panel headline; adversarial representativeness (random TAIL vs paralog/rDNA/SMA reality);
COMPASS-gate measurability; net-value for a discovery tool.

### CP3 — BASELINE REPRODUCED (their panel is real, close numbers hold)
`python dev/discovery_loss_panel.py --n 40 --seed {1,2}` reproduced INDEPENDENTLY:
- seed 1: ed_signal bal-acc 99% in delta[0.5,1.5]; WIRED m3/cap2/gate1 → disc-loss 0.7% / fab-residual 4.3%.
- seed 2: ed_signal bal-acc 100%; WIRED close → disc-loss 0.4% / fab-residual 0.0%.
Matches DISCOVERY_LOSS_PANEL_RESULT.md CP5 (0.4%/4.3%) within seed noise. ⇒ The close's synthetic result
is NOT fabricated — on THEIR random-TAIL panel the ed_signal genuinely separates the overlap band.
The question is entirely whether that panel represents real cryptic-in-microhomology loci.

**KEY DIAGNOSTIC (why the panel's signal is manufactured):** verified by hand — in their construction
the two candidate exon2 windows genome[ne:ne+28] vs genome[je:je+28] diverge at ~0.75/base (random TAIL),
so ed(rescue,inc) is LARGE and esig is strongly positive for a real cryptic. In a PERFECT period-k tandem
(paralog reality) the two windows are IDENTICAL → esig collapses to exactly 0 (verified: esig==0 for a
clean read; esig==2 with a single sparse paralog SNP). So the ed_signal's discriminating power is a direct
function of downstream distinctiveness — HIGH in their random panel, ~ABSENT in a tight paralog/rDNA array.

NEXT: run the paralog-tandem variant (perfect period-8 repeat + n_snp ground-truth SNPs, 6% ONT error) to
quantify how far the close degrades when the downstream is NOT random.

### CP4 — PARALOG-TANDEM VARIANT: the ed_signal is TECHNICALLY ROBUST (null result — strengthens HOLD)
Perfect period-8 tandem (near-identical windows) + n_snp ground-truth paralog SNPs in the cryptic copy,
6% sub / 3% indel, n=60/cell, seeds 1&2. Harness:
`/private/tmp/.../audit_v5/strategic-realdata_B/repeat_tail_variant.py`.
```
 seed1  n_snp=1..4:  esig>0 sep bal-acc = 100%/100%/100%/100%;  CLOSE disc-loss 0.0% all
 seed2  n_snp=1..4:  bal-acc = 83%/100%/50%/100%;               CLOSE disc-loss 0.0% all
```
**FINDING (honest, against my initial hypothesis):** the ed_signal does NOT collapse even in a tight
paralog tandem, as long as a recoverable SNP sits INSIDE the rescue window W=28. Reason: esig = ed_inc -
ed_cry is a DIFFERENCE — common-mode ONT error cancels; only the differential paralog base survives. So a
single in-window SNP holds esig ≈ +2 despite 6% error. ⇒ **The technical close is SOUND on the acceptor
side wherever recoverable truth lies within W.** This is a NULL for "the signal is fragile" — and it makes
the HOLD airtight: the fault is NOT a math defect; it is strategic (untested-where-it-counts).
CAVEAT (do not over-read): the fab-drift denominators are TINY (0–3 of 60) because canonical reads barely
drift in a tight tandem (fab-rate ~2–5%) → the "fab-residual 33%" at seed2/n_snp=1 is 1/3 reads = NOISE,
NOT signal. Informative direction: in exactly the paralog regime, the fabrication the guard exists to
suppress is ITSELF RARE (~2–5%), which further undercuts the guard's marginal value at paralog loci.

### CP4b — the ONE genuinely-different structural gap: OUT-OF-WINDOW SNP (reasoning + 1 check, not tuned)
The acceptor ed_signal only sees distinguishing bases within W=28bp of the junction. Verified:
```
 SNP at je+40 (BEYOND W=28):  esig(clean read) = 0   -> gate SPARE (psig>=1) FALSE -> cryptic VETOED by cap-alone
 SNP at je+5  (INSIDE  W):    esig(clean read) = 2   -> gate fires -> cryptic spared
```
Real paralog loci whose nearest discriminating position is >28bp from the cryptic acceptor are INVISIBLE
to the gate → the close silently reverts to cap-alone (~10% documented discovery loss) precisely there.
Their panel NEVER tests this (TAIL divergence always begins at offset k). But this is the CONSERVATIVE
failure direction (discovery loss, NOT fabrication), so it degrades the close's discovery-preservation
claim, not its safety. Still: a representativeness gap the synthetic evidence structurally omits.

### CP5 — THE STRUCTURAL REPRESENTATIVENESS GAP NO PANEL CAN CLOSE
In synthetic data, read ORIGIN is known by construction, so "fabrication vs discovery" is a clean label.
In real paralog / rDNA / SMA loci, BOTH the incumbent and the cryptic acceptor are REAL genomic positions;
the correct label depends on WHICH paralog copy the read physically came from — which single-read ONT
(6% error, no phasing) frequently CANNOT determine. No panel with baked-in ground truth can test the case
where the ground-truth LABEL ITSELF is unknowable. This is the irreducible representativeness gap, and it
is why real-data confirmation (COMPASS §4b) is not optional decoration — it is the ONLY evidence class that
can measure the guard's real fabrication/discovery trade. It was never produced.

---

## ★ VERDICT — task "strategic-realdata" — Auditor B: **HOLD** (fault_found = TRUE)

**The fault is a STRATEGIC VALIDITY GAP, explicitly NOT a technical defect in the ed_signal** (I verified
the acceptor signal holds even under an adversarial tight-paralog attack). "Validity established" is
STRUCTURALLY IMPOSSIBLE to claim from the evidence that exists, for four independently-sufficient reasons:

1. **The guard's TARGET metric was never measured on real data.** COMPASS §4b — the recovery/drift/
   inconclusive classification that IS the real-data fabrication rate — OOM'd, then TIMED OUT (algorithmic
   O(142k×list)), then was declared "context, not a gate." The single number the close exists to reduce has
   NO real-data value. (COMPASS_RECALL_RESULT.md §4b UPDATE.)
2. **The close was never RUN on COMPASS at all.** The 32× recall win (§1–3, raw-mm2 0.54% → arm-B 17.46%)
   is arm-E = `hp_drift_margin=3.0` ONLY — NOT the microhom margin + cap + positional gate. `refined.bam`
   does not exercise the close. So neither the recall side NOR the precision side of COMPASS touches the
   feature under audit. (RESULT.md L176: "microhom guard NOT yet applied".)
3. **The close guards a feature that cannot be turned on.** `correct_command.py` passes NO `motif_blind`
   (default False) and NO drift kwargs; there is NO CLI flag, NO config key, ZERO non-test caller in
   rectify/ that sets `motif_blind=True` or any drift param. The entire motif-blind re-placer + its 4 guard
   params are DARK CODE reachable only from dev/test scripts. 4 params of complexity guard an unshippable arm.
4. **All close evidence is one synthetic panel whose signal is manufactured by a property real target loci
   lack.** The ed_signal's discriminating power = downstream distinctiveness. Their panel's TAIL is RANDOM
   (~0.75/base divergence starting at offset k) → strong esig. Real paralog/rDNA/microsatellite loci (the
   SMA fabrication-enriched set) are near-identical tandems where distinctiveness is SPARSE and often
   OUT-OF-WINDOW (>28bp) → esig weak-to-zero at exactly the hardest loci (CP4b). And the label-unknowability
   of real paralog origin (CP5) is untestable by any ground-truth panel.

**Is the guard net-positive for a discovery tool? Not demonstrated, and arguably net-negative unguarded-by-
evidence.** It trades the tool's CORE VALUE (discovery — the 32× recall win) to suppress a fabrication class
(~1–2% synthetic; UNMEASURED on real data) that a discovery workflow ALREADY gates downstream (novel
non-canonical junctions are flagged for orthogonal validation regardless). The correct disposition for a
use-case-dependent precision/recall knob is DEFAULT-OFF, opt-in — which is EXACTLY what shipped. That
supports "honest stop," not "keep building."

**Q2 (measurability) — the sharpest point:** the existing COMPASS harness AS-RUN cannot answer the
enablement question, and finishing §4b would NOT fix it. Two gaps: (a) §4b never computed; (b) even a
completed §4b measures ARM-E's fabrication, not the close's (refined.bam has no microhom/cap/gate).
Justifying enablement requires a NEW COMPASS run with the close ENABLED on both arms — does it cut real-data
fabrication without eroding the 32× recall? — THEN the §4b classification. That new run is the honest gate;
it is feasible (binned-index fix specified) but has NOT been done.

**Q4 steelman — "acceptor close done, default-OFF, gate on COMPASS" WINS decisively over "keep building":**
the acceptor ed_signal is technically complete and adversarially robust (my CP4 null); the doc's own pre-CP4
verdict block (L133–144, never reconciled) already reached "deep scorer surgery buys a tradeoff-shift, hard
to justify on an unwired default-off guard." The rational next action is NOT more synthetic scorer work
(donor/both-boundary signal, out-of-window handling) on dark code — it is the ONE real-data measurement
(close-enabled COMPASS + §4b) that alone can decide whether the guard is worth its complexity. Until that
exists, the close is SYNTHETIC-ONLY INSURANCE whose worth is UNPROVEN.

**Bottom line:** the close is internally VALID and CORRECTLY SCOPED (default-OFF, byte-identical, acceptor-
only by design, signal technically sound) — but its VALIDITY AS A SHIPPING DECISION is NOT established. The
technical work is at a legitimate stopping point; the strategic claim ("the fault CLOSED") is premature
because the only evidence that could substantiate it (real-data, close-enabled) was never produced. HOLD:
mark the acceptor close DONE + default-OFF; DO NOT enable and DO NOT build further scorer complexity until a
close-enabled COMPASS run + the §4b fabrication classification exist. Independent of Auditor A.

