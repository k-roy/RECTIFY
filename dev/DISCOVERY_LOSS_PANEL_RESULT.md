# Discovery-loss panel — does the guard veto REAL cryptics in microhomology? (2026-07-13)

**Why (advisor + user "build the close"):** the load-bearing measurement the audit stalled on 4×. It is
step 1 of closing the read-blind fault: quantify the DISCOVERY-LOSS RATE (real cryptic-in-microhomology
reads wrongly vetoed) vs (microhom_drift_margin, drift_near_tie_cap), binned by delta_improve. Go/no-go:
if real cryptics reliably have delta_improve >= 3, then margin=3 (already committed, `05664bc`) SUFFICES and
the deep positional-distinctiveness scorer surgery is unnecessary; if real cryptics populate delta_improve in
(0,3) while tripping mh>=0.5, the positional signal is REQUIRED. Built INLINE (agent runs stalled 4×; local
python is reliable). Script: `dev/discovery_loss_panel.py`.

## Construction (real cryptic genuinely from the non-canonical site, placed at the incumbent decoy)
- genome = LPAD + exon1 + canonical intron(GT..AG) + EXON2; EXON2 = U (k-mer) + U' (U with `mm` mismatches) +
  TAIL (random). Incumbent acceptor ne = start of U; CRYPTIC acceptor je = ne+k = start of U'. Microhomology
  = U ~ U' (frac = 1 - mm/k). A read genuinely from the cryptic carries exon2 bases = genome[je:] = U'+TAIL,
  but is PLACED at the incumbent (N-op [ns,ne)) → the refiner should MOVE it ne→je (discovery).
- delta_improve = _score_junction(incumbent ne) - _score_junction(cryptic je) (>0 ⇒ read favors cryptic).
- mh = _move_microhomology(genome, ns, ne, ns, je). Sweep k, mm, TAIL len/divergence to span delta_improve.
- discovered (guard OFF) = refiner moves to je. vetoed (guard ON @ margin,cap) = held at incumbent.

## CHECKPOINT 1 — CLEAN panel (no read errors, n=20/cell, 1440 reads) — ENCOURAGING
Real cryptics genuinely from the non-canonical site, in microhomology (mh≥0.5), that the refiner discovers
guard-OFF (1438/1440 = the at-risk population). Their delta_improve distribution:
  min=1.00  p10=2.00  median=4.75  p90=8.00  max=10.00 ; **0.0% have delta_improve<1**, 14% have <3.
DISCOVERY-LOSS RATE (at-risk reads wrongly held at incumbent):
```
    margin  cap=0    cap=1    cap=2    cap=3
    m=3      14.0%    0.0%     4.2%    14.0%
    m=8      82.5%    0.0%     4.2%    14.0%
```
**KEY:** real cryptics carry GENUINE evidence (distinctive tail) → delta_improve ≥ 1 for ALL of them. So
**drift_near_tie_cap=1 ELIMINATES discovery loss (0.0%)** even at margin=8 — the cap veto only fires for
delta<1, which no real cryptic hits. m=8/cap=0 loses 82.5% (confirms m=8 is bad uncapped); m=3/cap=0 loses
14% (the delta<3 tail). This is the empirical picture the structural "cap can't close" worst-case argument
did NOT rule out: real cryptics and error-free fab drift are SEPARATED on delta_improve. The remaining
question the WITH-ERRORS + FAB side must answer: do ERROR-DRIVEN fab drifts (canonical reads whose ONT error
tips them to the cryptic) also stay below the cap while real cryptics stay above — i.e. is the cap a clean
separator on realistic data? (Without errors, fab reads have delta_improve≤0 and never drift — drift is
error-driven, so the fab side needs the error model.)

## ★★★ CHECKPOINT 2 — WITH ONT ERRORS, BOTH SIDES (n=40/cell, 5760 reads, seeds 1&2 agree) — DECISIVE
Realistic ONT error model (6% sub, 3% indel). Two populations, both in microhomology (mh≥0.5):
- **DISCOVERY (real cryptics, genuinely from je, discoverable guard-OFF):** delta_improve median≈4.0,
  p10≈2.0; only ~1% <1, ~10% <2, ~24% <3. They carry genuine distinctive evidence.
- **FAB (canonical reads that ONT-error-drift to je), guard-OFF fab rate ~2.4%:** delta_improve
  median≈1.0, max≈1–2; ~25% <1, ~97–100% <2. Error-driven near-ties.

TRADEOFF (seeds 1 / 2 agree within ~1%):
```
   m  cap   disc-LOSS   fab-RESIDUAL       (disc-loss = real cryptics wrongly held; fab-residual = fab still drifts)
   3   0      ~24%          0%             ← margin=3 ALONE loses a QUARTER of discovery (INSUFFICIENT)
   3   1      ~1%          ~75%            ← cap=1 protects discovery but too small to catch fab
   3   2     ~10%          ~0–3%           ← BEST cap point: ~10% discovery loss, fab ~eliminated
   3   3     ~24%           0%             (== cap=0 since min(3,3)=3)
   8   0     ~91%           0%             ← m=8 uncapped: catastrophic
```
**CONCLUSIONS (go/no-go answered):**
1. **margin=3 does NOT suffice** — ~24% discovery loss uncapped. The optimistic "one measured case was
   delta=6, so margin=3 is fine" hypothesis is REFUTED: real cryptics have a delta_improve TAIL (24% <3).
2. **The near-tie cap helps materially** — m=3 / cap=2 cuts discovery loss 24%→~10% while ~eliminating fab.
   Best defensible cap operating point = **margin=3, drift_near_tie_cap=2, threshold=0.5**.
3. **~10% discovery loss REMAINS and is IRREDUCIBLE on the delta axis** — real cryptics and fab drifts
   OVERLAP in delta_improve ∈ [0.5, 2]. This EMPIRICALLY CONFIRMS the audit's structural finding: the cap
   (a delta-axis bound) cannot CLEANLY separate them. Only a POSITIONAL-DISTINCTIVENESS signal (orthogonal
   to aggregate delta) could drive discovery loss toward 0 while keeping fab suppression — IF it separates
   the overlap band. That is the motivation for the deep close, now empirically justified.
4. **Caveat (the honest trade):** to suppress a ~2.4% fab rate the guard costs ~10% discovery loss even at
   the best cap. For a tool whose VALUE is discovery (COMPASS: 32× recall win), that trade is only worth it
   if the positional signal makes it near-free. So the positional signal is the pivotal next test.

NEXT (per advisor): before building the positional signal, VERIFY ORTHOGONALITY on 2–3 hand cases — a
fab-drift and a real-cryptic with the SAME delta_improve (~1, in the overlap) — does a discriminating-column
signal separate them? If yes → build it; if no → the fault can't be closed and the honest stop is
"m=3/cap=2, ~10% documented discovery loss, opt-in" (with COMPASS as the real enablement gate).
Script: `dev/discovery_loss_panel.py`. Reproduce: `python dev/discovery_loss_panel.py --n 40 --seed {1,2}`.

## ★★★ CHECKPOINT 3 — POSITIONAL SIGNAL: orthogonal but NOT a clean close (the pivotal test)
ORTHOGONALITY (naive discriminating-column vote, indel-sensitive) within the delta OVERLAP band [0.5,1.5]
where the cap fails: **cry (real) votes cryptic 83%, fab (drift) votes incumbent 73%** — genuinely
orthogonal to delta (separates what the delta axis cannot). So the positional signal HAS real information.

DECISIVE (m=3, cap=2; positional gate SPARES a would-be-veto when sig>0), seed 1:
```
    cap-alone         disc-loss  9.8%   fab-residual  2.9%
    cap+pos(sig>0)    disc-loss  1.8%   fab-residual 30.0%
    cap+pos(sig>=1)   disc-loss  6.1%   fab-residual 17.1%
```
**The positional signal is a TRADEOFF-SHIFT, NOT a clean close.** Adding it cuts discovery loss (9.8%→1.8%)
but RAISES fab-residual (2.9%→30%) — because the naive signal misclassifies ~27% of fab (votes cryptic →
spares them). In ABSOLUTE terms cap+pos(sig>0) ≈ 1.8% discovery loss + ~0.7% fabrication vs cap-alone ≈ 9.8%
loss + ~0.07% fab. An indel-ROBUST version (align first) would sharpen it, but the deep-overlap tradeoff is
FUNDAMENTAL: a single ONT read from a real cryptic-in-microhomology and an error-driven fab drift are
genuinely hard to separate. No scalar (delta / cap / naive positional) eliminates the precision/recall
tradeoff; each just moves the operating point along the curve.

**Seed 2 CONFIRMS (n=40, 5760 reads):** margin=3 alone 24.7% loss; m=3/cap=2 → 9.5% loss / 0% fab;
orthogonality cry 83% / fab 67% in [0.5,1.5]; cap+pos(sig>0) → 1.9% loss / 33.3% fab-residual. Robust.

## ★★★★ CHECKPOINT 4 — INDEL-ROBUST positional signal essentially CLOSES the fault (2026-07-13)
The naive per-index signal (CP3) was a tradeoff-shift because it's indel-sensitive. The INDEL-ROBUST version
— **hard-anchored edit distance** of the read's rescue to the two candidate exon2s, `ed_signal =
ed(rescue, genome[ne:]) - ed(rescue, genome[je:])` (NO free-prefix split ⇒ removes the scorer's soft-clip
escape that hid the discriminating mismatches) — separates the overlap band nearly perfectly:
```
  separation in delta band [0.5,1.5] (where the cap fails):
    naive sig : cry 84% / fab 78%  balanced-acc 81%
    ed_signal : cry 99% / fab 98%  balanced-acc 98%     ← near-clean
  CAP-ALONE vs CAP+ED-POSITIONAL (m=3, cap=2):
    cap-alone           disc-loss 8.8%   fab-residual 2.1%
    cap+naive(sig>0)    disc-loss 1.8%   fab-residual 23.4%   (naive: bad)
    cap+ed(esig>0)      disc-loss 0.4%   fab-residual 4.3%    ← near-CLOSE
    cap+ed(esig>=2)     disc-loss 6.2%   fab-residual 2.1%    (balanced knob)
```
**⇒ THE CLOSE IS ACHIEVABLE.** cap+ed(esig>0) ≈ 0.4% discovery loss AND 4.3% fab-residual — the align-first
signal (removing the free-k soft-clip escape) drives discovery loss to near-zero while keeping fab
suppression. The earlier "tradeoff-shift, not a close" verdict was an artifact of the NAIVE signal; the
indel-robust signal genuinely closes it. NEXT: wire ed_signal into the refiner as the positional veto-gate
(the deep fix the user approved), re-validate, re-audit. (seed 2 confirming.)

## ★★★★★ CHECKPOINT 5 — WIRED end-to-end, the CLOSE is shipped (default OFF) (2026-07-13)
`_positional_signal` (hard-anchored `_semiglobal_ed`) wired into `refine_read_junctions` as
`drift_positional_gate` (default 0.0 = OFF = byte-identical), threaded through all 4 refine fns. A
drift-flagged would-be-veto is SPARED when the read's positional signal ≥ gate. TRUE end-to-end through the
real refiner (not overlay), seed 1:
```
    WIRED m3 / cap2 / gate1   disc-loss 0.4%   fab-residual 4.3%   <- shipped close (== overlay prediction)
```
vs margin-alone ~24% and cap-alone ~8.8% discovery loss. Tests: 30 in test_microhom_drift_guard.py
(_semiglobal_ed, _positional_signal, veto-band-spare, byte-identical-off); 400/400 veto-band real cryptics
spared end-to-end; broad refiner/validation suite byte-identical at default. Operating point (opt-in):
microhom_drift_margin=3, drift_near_tie_cap=2, drift_positional_gate=1, microhom_threshold=0.5. STILL default
OFF. Remaining: donor/both-boundary exon1-side signal (conservative now); re-audit; COMPASS real-data.

## ⇒ VERDICT (pre-CP4; superseded by CP4/CP5 — the indel-robust signal CLOSES it)
1. **margin=3 alone is INSUFFICIENT** (~24% discovery loss). Refutes the optimistic hypothesis.
2. **Best cap operating point = margin=3, drift_near_tie_cap=2, threshold=0.5** (~10% discovery loss,
   ~0% fab-residual). This IS a defensible opt-in operating point.
3. **The positional-distinctiveness signal is orthogonal and real (83/73 sep) but does NOT cleanly CLOSE
   the fault** — it shifts the precision/recall tradeoff (discovery↓ for fab↑), it doesn't eliminate it.
   The microhomology overlap is an IRREDUCIBLE single-read ambiguity.
4. **HONEST RESOLUTION:** the guard is a precision/recall TRADEOFF tool. Ship it DEFAULT-OFF, opt-in, with
   the measured tradeoff curve documented, and let the user pick the operating point (discovery-priority:
   small/no guard; precision-priority: cap=2). Deep scorer surgery for the positional signal buys a
   tradeoff-shift, not a clean close, on an unwired default-off guard — hard to justify vs. the marginal,
   fundamentally-bounded gain. The real ENABLEMENT gate remains COMPASS real-data confirmation (independent).
