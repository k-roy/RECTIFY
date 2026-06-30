# Flat-Q / per-base-quality arbitration axis — gate

## VERDICT (2026-06-30) — Q adds NO gap over the error-type table; flat-Q caveat MOOT

**DECISIVE NULL.** `scripts/benchmark/flatq_headroom.py` (reps=120, 5628 indel-strata
reads, oracle Q, injected sub-only burst corruption at measured sub_rate=4.4%):

| stratum | n | ceil | hp_ed | Qcost | disagr | Qgain | Qbreak | Qgain\|dis |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HP | 2635 | 0.745 | 0.745 | 0.734 | 0.015 | 0.000 | 0.011 | 0.000 |
| HP_HARD | 2400 | 0.573 | 0.573 | 0.384 | 0.217 | 0.000 | 0.190 | 0.000 |
| STR | 593 | 0.728 | 0.728 | 0.728 | 0.000 | 0.000 | 0.000 | 0.000 |
| **TOTAL** | 5628 | 0.670 | 0.670 | 0.584 | 0.099 | **0.000** | 0.086 | **0.000** |

- **hp_ed == ceiling on EVERY stratum** — the incumbent error-type arbiter already
  picks a truth-correct member on 100% of recoverable reads, even on HP_HARD where
  members genuinely DISAGREE on 21.7% of reads. This extends C3's headroom=0 finding
  from clean reads to ERROR-CONCENTRATED reads (the advisor's "run it where it
  matters" requirement).
- **Qgain = 0.0000** (and Qgain|disagreement = 0.0000): even ORACLE Q — the
  perfectly-calibrated basecaller, the most generous possible case — recovers ZERO
  reads hp_ed misses. A saturated arbiter leaves no room for ANY orthogonal
  covariate (Panel C's information-theoretic point: I(correct;pick) is already
  maxed).
- **Q is strictly WORSE as an arbiter** (Qcost 0.584 < hp_ed 0.670; Qbreak=8.6%):
  the full-query Q-emission cost is dominated by per-base substitution likelihood
  and is BLIND to the indel-placement signal hp_ed captures, so swapping in Q
  actively breaks reads. Quality is not a free add-on here — it discards the very
  signal the error-type table exists to weigh.

**Therefore: the flat-Q caveat (SPEC:225) is MOOT. The error-type gates stand; the
quality-blind panel loses nothing by ignoring `query_qualities` on this corpus.**
No Q decoder should be built. The result is robust because oracle Q is an UPPER
BOUND (no noisier real Q can beat what an omniscient one already fails to win), so
it needs no SIRV magnitude calibration.

### Triple adversarial panel — CONVERGED (all 3 endorse the null)
- **A (circularity):** oracle Q makes a NULL decisive but a WIN non-credible (the
  injector hazard has no reference-context term). read_mult is inert for same-read
  arbitration; only burst localizes. → one-sided framing adopted; null is clean.
- **B (zero-evidence):** arbitration zero-evidence guard sound; Q must not be built
  if headroom=0 (it is 0). Discovery worst-case (clean-read + junction-local burst)
  is UNINSTANTIABLE on the current corpus → discovery stays SKETCH-only.
- **C (steelman):** arbitration airtight (saturated arbiter → no covariate adds);
  the only false-positive risk is DISCOVERY-by-construction — guard with
  flank-windowed-density RESIDUAL + realistic-Q + decoupled-hotness before any claim.

### DISCOVERY axis — SKETCH ONLY (not built; verdict: no actionable gap shown)
The exonic-error-DENSITY proxy (`novel_support_probe.py` axis-b) already tracks the
hot-read hazard at r=0.955. Q is a pre-alignment per-base covariate; the ONLY place
it could structurally beat density is LOCAL (junction-flank, where the suspect
alignment distorts density). But (Panels B+C): on this corpus spurious junctions are
planted in hot reads and the non-circular Q is drawn from the same hazard, so any
Q-beats-density result is construction, not orthogonal signal; AND the
clean-but-locally-bursty class needed to exercise the genuine local advantage does
not exist in the corpus. **Where it would plug in IF ever built:** a soft per-call
down-weight in `novel_support_probe.score_axis`, gated by `head_tail_autocorr` so
only GLOBALLY-elevated reads are down-weighted (preserving a real novel that merely
has a local junction dip). Non-circularity bar (pre-registered, Panel C): Q must
beat the density RESIDUAL (nested/partial-AUC) under realistic quantized+miscalibrated
Q with spurious-placement decoupled from hotness (purity→p_spectrum). Until a corpus
with those properties exists (SIRV-gated), the discovery quality-gap is UNPROVEN, not
demonstrated — consistent with the moot verdict.

---
# (design notes below)



The one untested orthogonality lever the C3 gate left open (SPEC:225). Every gate
so far (C1..C6) ran on FLAT-quality sim reads (`'I'*len`), so each tested only the
ERROR-TYPE axis. The panel's headline orthogonality claim is that RECTIFY is
**quality-BLIND**: `query_qualities` is propagated-but-unconsumed. This gate asks:
**does per-base quality carry arbitration/discovery signal the error-type table
(hp_edit_distance / the C1 penalty table) does NOT — a quality-addressable gap?**

Discipline (mirrors C1/C3): MEASUREMENT before any decoder, fitness = TRUTH,
pre-committed null, zero-evidence guard. NO LLR / Q-decoder is built unless the
LLR-free headroom probe shows a gap.

## Pre-committed verdict (BEFORE running)
- **Q adds no headroom over the error-type table** → flat-Q caveat MOOT; the
  error-type gates stand. (The expected outcome IF the incumbent arbiter is already
  at ceiling — C3 showed headroom=0.000 on indel strata; an orthogonal axis cannot
  fill a gap that does not exist.)
- **Q opens a real, non-circular gap** → report as the live quality-axis
  opportunity (SKETCH where it plugs in; do NOT build).

## Step 0 — PREMISE confirmed (query_qualities propagated-but-unconsumed)
- WRITE sites (passenger copy only): `chimeric_consensus.py:952`
  (`out.query_qualities = template_read.query_qualities`), `consensus.py:558`
  (`best_read.query_qualities = donor_read.query_qualities`), `mpb_split_reads.py`
  (chunk reassembly 234/275/371/414), `multi_aligner.py:1789/1792` (FASTQ inject).
- READ sites: `terminal_exon_refiner.py:714` slices `qual[:length]` ONLY to carry
  the soft-clip's quality string as a passenger — NOT consumed in any cost/score.
- grep across consensus/correct/splice/align: NO arbitration, scoring,
  correction, or selection function reads a base quality. **PREMISE CONFIRMED.**
- `cigar_records_to_bam` Q gotcha: TODO confirm (the sim BAM-izer path).

## The circularity problem (panel angle i — the headline risk)
A Q-win is FAKE if Q is injected as a deterministic readout of the realized error
positions and the scorer just reads them back. NON-CIRCULAR design rule:
- Q is drawn from the **latent per-position HAZARD** (`base_rate × read_mult ×
  burst_factor`), NOT from the realized `ErrorEvent` track. A low-Q base is MORE
  likely an error but not deterministically one; errors can occur at moderate Q.
- The error-type table is a function of **reference context only** (HP length,
  base class). The hazard's `read_mult` (layer 1) and `burst_factor` (layer 2) are
  per-read / per-position stochastic and **independent of reference context** — so
  Q carries an axis the error-type table structurally cannot see. That orthogonality
  is the non-circularity argument, and it is also why a Q arbiter cannot
  hill-climb the reference-context table.

## ARBITRATION axis (the c3_headroom analog, Q-aware)
For arbitration between two PLACEMENTS of the SAME read, Q is identical for both
placements (same query bases). Q can only re-rank if the two placements assign the
read's low-Q bases to different roles (error vs match). KEY: if the error-type
arbiter is already AT CEILING (C3: headroom 0.000), there is NO room for Q to add
arbitration accuracy. Measurement: inject informative Q, add a Q-aware arbiter,
ask whether it recovers ANY read hp_ed misses (expect 0) AND breaks none hp_ed
gets right (safety).

## DISCOVERY axis (vs novel_support_probe axis-b)
`novel_support_probe.py` already uses **exonic error DENSITY** (alignment-derived,
axis b) as the deployable reliability proxy. Q is a PRE-alignment, per-base
covariate. Genuine gap candidate: near a fabricated junction the alignment (hence
the density proxy) is itself distorted; Q reads reliability in the junction-flank
WITHOUT depending on the suspect alignment. Test: does Q add LIFT OVER the
exonic-density proxy in the overlap (low-support) zone? If Q is just a noisier
readout of the same hotness density already captures → no gap.

## PANEL A (circularity) — decisive refinements [integrated]
1. **read_mult is INERT for same-read arbitration.** It is constant within a read,
   scales every position equally → cannot re-rank two placements of the SAME read.
   Only the **burst factor** (layer 2) carries arbitration signal. (My orthogonality
   claim wrongly leaned on read_mult; corrected.) read_mult DOES matter for
   DISCOVERY (per-read reliability varies across reads).
2. **The injector hazard has NO reference-context term** (`base_rate` is a scalar;
   only read/position stochastics modulate it). So injected errors are HP-blind,
   while hp_ed is keyed to HP context. ⇒ oracle Q makes a **NULL decisive** (even
   omniscient Q can't beat hp_ed → caveat moot) but a **WIN non-credible** (oracle Q
   would be reading the generative latent at zero noise, graded against the same
   latent; real basecaller Q is strongly ref-context-correlated, unlike here).
3. **Pre-registered consequence:** frame arbitration as a ONE-SIDED test. If even
   oracle Q yields ZERO headroom over hp_ed on error-concentrated/disagreement
   reads → DECISIVE NULL, caveat moot. If oracle Q DOES fill headroom → report it as
   NON-CREDIBLE (requires (a) re-run under degraded/calibrated Q AND (b) a
   ref-context hazard term so hp_ed and Q genuinely compete) → still no actionable
   gap on this corpus.

## Implementation notes
- `cigar_records_to_bam` (scorer.py:280-289) never assigns `query_qualities` →
  member BAMs carry NO Q. Carry oracle-Q in a SIDE dict `{read_id: qual_array}`
  computed at injection; the Q-aware arbiter walks each member CIGAR and accumulates
  a Q-weighted cost from that dict. No BAM-izer change needed.
- Background injection restricted to SUBSTITUTIONS (advisor): injected background
  INDELs trip the scorer's `has_unexplained` gate (any CIGAR indel outside a truth
  span → not position-exact) and collapse concordance spuriously (SESSION-2 lesson).
- Arbiter consumes ONLY (candidate alignment, realized Q string) — never
  hazard/mult/burst/ErrorEvent. That barrier = deployment conditions (non-circular
  by construction regardless of how Q was generated).

## Zero-evidence discipline (panel angle ii)
A read clean in exons but bursty at the junction must stay RECOVERABLE — Q must be
a SOFT down-weight, never a hard gate. A read whose low-Q is global (hot) vs local
(junction burst) must be distinguishable (the head_tail_autocorr axis).

## Files
- `dev/FLATQ_QUALITY_AXIS.md` (this doc)
- `scripts/benchmark/flatq_headroom.py` (the LLR-free Q-headroom probe — TODO)
- `scripts/benchmark/flatq_headroom_result.txt` (result — TODO)
- Shared-file edits (flagged separately, default-OFF + byte-identical proof): TBD
