# C2 — 3'/poly-A CPA placement: gate result + design (DRAFT 2026-06-29)

Native placement of the transcript 3' end (cleavage-and-polyadenylation site).
The hard case: a true CPA abutting a GENOMIC A-tract — the poly-A tail (read)
blends into a downstream genomic A-run (reference), so a pure-alignment 3' end
drifts downstream into the tract.

Mirrors the C1 gate discipline (`dev/C1_DESIGN.md`): establish the stratum is
DISCRIMINATING (incumbent below ceiling vs TRUTH) AND C2-ADDRESSABLE (a
change-point decoder fixes it, not generic noise) BEFORE designing the decoder,
with a pre-committed null and an over-call safety control.

## The decisive finding: the incumbent is the WALKBACK, not raw minimap2

RECTIFY already ships a guarded heuristic walkback
(`rectify/core/correct/walkback.py::walkback_drs_full` →
`walkback_3prime_guarded`): from the 3' alignment end it walks INWARD through the
stop-base (A) run to the first non-stop read==ref match, with three artifact
guards (large-deletion pre-scan, real-intron clip, tail-context false-stop). It
also ships `PolyAModel` (`rectify/core/polya/polya_model.py`) — a learned
tail-length law (mean 28.5, std 15.2), a position profile, and non-A frequencies
(G 4% / C 2% / T 1%) — trained on CONTROL SITES (0 downstream A's) precisely
because the project already knows genomic-A is the confound.

So the C2 matched-arm trap (the analog of C1's `chrom_ref`): scoring `est =
reference_end - 1` measures RAW minimap2, which drifts the full tract and which
the walkback ALREADY closes — a phantom facet. The baseline arm MUST be the
shipped walkback. C2 only earns its keep if it beats `walkback_drs_full`.

## The gate (RUNNABLE): `scripts/benchmark/c2_gate.py`

Three matched arms on a FIXED deterministically-drifted alignment (M through the
genomic A-region, S the tail; no aligner invoked on toy contigs):
- **raw** — `reference_end-1` (context only; already-solved by walkback).
- **walkback** — `walkback_drs_full` (the REAL incumbent to beat).
- **decoder** — prototype 2-state templated-vs-tail change-point: pull the 3' end
  to the 5' start of the maximal **gap-tolerant** 3'-terminal A-run (a degenerate
  tail-emission model tolerating up to `max_gap` non-A, the ~7% PolyAModel allows).

Generator (`scripts/benchmark/sim/controlled.py::genomic_a_cpa_cells`, wired into
`generate_corpus` as the `GENOMIC_A_CPA` stratum; truth carries `true_cpa` +
`downstream_a_count`, which the scorer's `cpa_abs_errors` already consumes):

| cell | construct | true CPA | identifiable? | role |
| --- | --- | --- | --- | --- |
| tract_start | body + A×g + nonA | body end (whole tract non-transcribed) | yes | NON-DISCRIMINATING control |
| readthrough | body + A×g + nonA, cleave k into the run | body end + k | **NO** | pre-committed NULL (iron triangle) |
| interrupted | body + A×a1 + X + A×a2 + nonA | body end | **disputed** | candidate addressable cell |
| clean_g0 | body(nonA end) + nonA, no A-tract | body end | yes | over-call control |
| terminal_A | body + A×t (REAL templated) + nonA | last real A | NO | SHARED over-trim limitation |

## RESULT (2026-06-29, reps=100, median / mean (exact-rate))

| cell | ident | raw | walkback | decoder |
| --- | --- | --- | --- | --- |
| tract_start | yes | 8.0 / 8.50 (0.00) | **0.0 / 0.00 (1.00)** | 0.0 / 0.00 (1.00) |
| interrupted | NO | 9.0 / 9.50 (0.00) | 4.0 / 4.00 (0.00) | 0.0 (artifact) |
| readthrough | NO | 4.0 / 4.50 (0.00) | 4.0 / 4.00 (0.00) | 4.0 / 4.00 (0.00) |
| clean_g0 | yes | 0.0 / 0.00 (1.00) | 0.0 / 0.00 (1.00) | 0.0 / 0.00 (1.00) |
| overcall_Aend | yes | **0.0 / 0.00 (1.00)** | **0.0 / 0.00 (1.00)** | 2.0 / 2.72 (0.00) |
| terminal_A | NO | **0.0 / 0.00 (1.00)** | 4.5 (0.00) | 4.5 (0.00) |

raw drift grows EXACTLY linearly with the genomic-A length g (g=3→3, 6→6, 10→10,
15→15); walkback = 0 across all g. `interrupted` is labelled UNIDENTIFIABLE (its
"win" is the truth-definition artifact below). `overcall_Aend` is the STRONG
over-call control: the decoder over-pulls (2.72) where the walkback is exact (0.00).

### What the gate proves
1. **The canonical genomic-A-CPA drift IS real vs raw** (raw absorbs the full
   tract, error = g). DISCRIMINATING vs raw minimap2.
2. **The shipped walkback ALREADY solves the canonical case at ceiling**
   (tract_start err 0.00 across all g). So the canonical stratum is
   NON-DISCRIMINATING against the real incumbent — the HP-vertical-slice
   saturation lesson, one facet over. C2-as-placement earns nothing on the
   canonical drift; that gap is already shipped.
3. **readthrough is the pre-committed NULL and it holds**: decoder == walkback ==
   k. Both snap to the run start; the body/tail boundary inside a contiguous
   A-run carries ZERO sequence information and the tail-length law (std 15.2) is
   far too broad to pin ±k. A decoder "win" here would be the C1 Claim-B artifact.
   It correctly does NOT win.
4. **Over-call: the decoder FAILS the strong control.** The weak `clean_g0`
   passes (err 0) but is rigged (`_nonA_flank` forces 3 trailing non-A). The
   STRONG `overcall_Aend` control (body ends `...A Z`, truth = Z) shows the decoder
   over-pulls across the true CPA (2.72) while the walkback is exact (0.00) — the
   prototype is NET-HARMFUL on realistic A-ending 3'UTRs.
5. **terminal_A exposes a SHARED limitation**: raw is CORRECT (0.00) but the
   shipped walkback OVER-TRIMS real templated 3'-terminal A's (4.5); the decoder
   inherits the same error. From the read+alignment, terminal_A and tract_start
   are STRUCTURALLY IDENTICAL (M through some A's, then S of tail A's) — only
   whether the aligned A's were templated differs, which is unidentifiable. The
   walkback's "all near-3' aligned A's are tail" prior helps tract_start and HURTS
   terminal_A; no sequence-only decoder can separate them.

### `interrupted` — ADJUDICATED ARTIFACT (relabelled UNIDENTIFIABLE 2026-06-29)
The decoder appears to beat the walkback on `interrupted` (0.0 vs 4.0) ONLY
because the cell's truth is set to the body end Z. The 3-reviewer adversarial
panel (+ 2 advisor passes) refuted this on four independent grounds:

1. **Truth-flip inversion (model-free, decisive).** Relabel the cell's truth from
   Z to the likelihood-MAP position X (= `z+a1+1`, the last *provably*-templated
   base) and the two arms EXACTLY SWAP — walkback 4.0→0.0, decoder 0.0→4.0. No
   sequence signal moved; the "win" is 100% a function of which coordinate is
   declared truth. The gate measured a DEFINITION, not a capability.
2. **The likelihood favors read-through 26-113×.** Transcript = genome[:CPA+1] +
   tail. Truth-at-Z needs a non-A tail base landing exactly on genomic X
   (~p_X·pA^a1, a ~1-4%×positional coincidence); read-through makes X templated
   ((1-err)^(a1+1)). LR(read-through/Z) = 26-113× (≈37× at err=0.15, ≈20× at
   err=0.30; the tail-length prior moves it <2×). The walkback's anchor-AT-X is
   the maximum-likelihood estimate; the decoder's pull-past-X is the ~50× LESS
   likely interpretation. Same iron triangle that sank C1 Claim B, subtler.
3. **Construction-tuned (red-team, empirical).** The decoder wins iff `max_gap`
   exactly equals the injected interruption width: 1-base X + `max_gap=1` → win;
   2-base interruption → win vanishes (ties walkback). Tuned to the cell.
4. **One basecall error erases it.** The corpus is noise-free, which structurally
   DISABLES the walkback's tail-context guard (fires only on `ctx_has_mismatch`).
   Inject a single non-A read error in the tract zone and the guard skips X →
   walkback 0.0 on the realistic g=10,15 cells, tying the decoder. The surviving
   "win" is confined to noise-free g=3,6 — the least realistic cell.

### Over-call: the decoder is NET-HARMFUL (strong control, red-team Attack 5)
The original `clean_g0` over-call control was RIGGED — `_nonA_flank` forces 3
trailing non-A, the exact pattern that hides over-call. The added STRONG control
`overcall_Aend` (body ends `...A Z`, Z a non-A CPA with a genomic A immediately
5', no downstream tract; truth unambiguously = Z) exposes it: **decoder mean err
2.72 vs walkback 0.00** — the gap-tolerant decoder over-pulls across the true CPA
into the body A, while the walkback is exactly right. There is no single `max_gap`
that both absorbs realistic interruptions (needs ≥2) and avoids over-calling
A-ending 3'UTRs (needs ≤1). This is the C1 insertion-discount lesson exactly:
a permissive cost that rewrites correct calls — and real bodies routinely end in A.

## VERDICT — C2-as-PLACEMENT is REFUTED (panel-unanimous, advisor-confirmed)
The shipped walkback is AT CEILING (err 0.00 across all g) on the identifiable
canonical genomic-A-CPA drift; readthrough is sequence-unidentifiable (decoder
ties, the held NULL); terminal_A is a shared irreducible over-trim (≡ tract_start
input with opposite truth — direct proof the boundary is unidentifiable); the only
apparent decoder win (`interrupted`) is a truth-definition artifact; and the
prototype decoder is NET-HARMFUL (over-calls A-ending bodies). **No identifiable
cell has a defensible decoder win.** The gate's printed verdict now reads REFUTED.
This is a valid, valuable gate outcome — the discipline working exactly as in C1's
insertion-discount gate-off: prove-don't-assert killed a plausible facet.

## The pivot — C2-as-SOFT-POSTERIOR (not earned as new code yet; see caveats)
The change-point model's only non-refuted product is a CALIBRATED UNCERTAINTY over
the boundary (the genomic A-run length IS the uncertainty). TWO hard constraints
the panel imposed, both of which currently BLOCK shipping new code:
- **The uncertainty width ALREADY ships.** `bam_processor.py:851-859` already sets
  `ambiguity_min/ambiguity_max/ambiguity_range = wb_bp` (the walk-back span) and
  tags `correction_applied=polya_walkback`. A change-point posterior would add only
  a per-base shape INSIDE that already-emitted window — and there is **no consumer**
  (consensus uses a categorical `confidence∈{high,med,low}`, no continuous CPA-LR
  arbitration). Building a posterior with nothing to read it is gold-plating.
- **It can NEVER be a hard trim/don't-trim rule.** Because terminal_A ≡ tract_start
  as inputs, their posteriors are IDENTICAL; any hard guard either declines on both
  (re-breaks tract_start, where trimming is correct) or trims on both (re-breaks
  terminal_A). It is a prevalence-dependent tradeoff, not a strict improvement.
  Validating such a guard on a terminal_A-only stratum would be "C2-Claim-B-prime"
  — the same fiat-truth trap one level up. Posterior is legitimate ONLY as a soft
  signal feeding a consumer that can use uncertainty, never a placement decision.

## Decoder design (RECORDED for completeness; do NOT build per the verdict)
2-state change-point: state T (templated, emit = genome match) vs P (tail, emit
from PolyAModel: A 0.91/G 0.04/C 0.02/T 0.01 + position profile); one T→P
transition at the CPA; MAP boundary under the tail-length prior. Plug-in point =
a REFINEMENT of the inherited 3' end where `walkback_drs_full` sits in
`bam_processor` (not a re-decode). The prototype `changepoint_cpa_plus` is the
degenerate gap-tolerant form; both the prototype AND the full likelihood form
inherit the read-through unidentifiability (they cannot beat the walkback's MLE
anchor) and the A-ending over-call, so neither is shippable as a placement facet.

## Real-data check (NATIVE DRS only — cDNA has no real poly-A 3' end)
Real DRS has NO per-read CPA truth, so a real-data arm CANNOT measure placement
accuracy. It can only be a SAFETY/consistency check at 0-downstream-A control
sites (how `PolyAModel` is trained). Moot given the REFUTE — recorded so nobody
later expects a real-data placement number.

## Build plan (what actually ships)
1. **DONE — the runnable gate**: `scripts/benchmark/c2_gate.py` +
   `controlled.genomic_a_cpa_cells` (the `GENOMIC_A_CPA` stratum, wired into
   `generate_corpus`; truth carries `true_cpa`+`downstream_a_count`, consumed by
   the scorer's `cpa_abs_errors`). Prints REFUTED with the full per-cell table.
2. **Documentation only** (the C1-idiom outcome — no product-code facet):
   - this design doc records the gate result + the panel adjudication;
   - one sentence on `walkback_drs_full`/`walkback_3prime_guarded` noting the
     terminal-templated-A over-trim is the ACCEPTED cost of the DRS CPA+tail prior,
     cross-referencing that `ambiguity_range` already exposes the affected span.
3. **Deferred (behind a real consumer)**: a soft CPA posterior would only earn its
   keep once consensus has continuous-LR arbitration (concept #3) to consume it;
   gate it on that, mirroring C1's "re-enable only after the validator exists."
4. **NOT a `c2_ablation.py`** — there is no defensible identifiable win to ablate;
   the gate IS the deliverable. (Contrast C1, which had a real Claim-A ablation.)

## Open risks / honest residue
- The `GENOMIC_A_CPA` stratum is now in the corpus; the `terminal_A`/`readthrough`
  cells carry truths that are unidentifiable by construction — they are CONTROLS,
  not scored placement targets. Any future consumer must treat them as such.
- The gate is NOISE-FREE; the red-team showed that overstates the incumbent's gap
  (one basecall error closes `interrupted`). If C2 is ever revisited, inject
  realistic basecall noise (the `error_injector.py` path) before trusting any win.
- The decoder prototype is retained in `c2_gate.py` ONLY as the refuted arm; it is
  not exported to `rectify/` and must not be wired into the pipeline.
