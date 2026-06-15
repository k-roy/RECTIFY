# Round-2 cDNA "discovery → assignment" — independent senior review (Opus)

**Date:** 2026-06-15 · **Reviewer:** independent Opus pass, at Kevin's request, to pressure-test
the Sonnet Phase-0 NO-GO (`PHASE0_VERDICT.md`, `EXPLORED_ROUTES.md §1`).
**Scope:** re-examined the apparatus, re-derived the verdict logic, and — for the first time —
**looked at the individual winning reads** (the step the prior pass explicitly deferred: "attributed,
not proven-genuine, not individually IGV-vetted").

**Bottom line:** the prior NO-GO **stands on frequency**, but its framing was incomplete in both
directions. The mechanism is **real** (3 clean existence proofs found), AND the negative is
**sharper than reported** once the attribution leak is understood. The productive reframe:
*the rescuable population is real but structurally tiny on DRS, and the reason points to a
specific untested substrate (full-length PCR-cDNA), not to more DRS cell lines.*

---

## 1. What holds from the prior verdict (load-bearing negatives — confirmed)

- **F3 — tiny, stable structural-win count.** ~7 spliced structural wins on A549 (4,579 eligible)
  and ~7 on brain (1,065 eligible) — flat across a 4.3× substrate-size range, RNA002↔RNA004,
  SRRM4-low↔SRRM4-high. A signal that does not scale with substrate size or biology is a noise
  floor. **Solid.**
- **F4 — 0/844,262 genome-unmapped reads get a confident spliced cDNA rescue** (identity-ceiling
  robust at ≤0.83). The unmapped pool is exactly where genuine seed-failures would live; it is
  empty. **This is the single most decisive negative and I did not find a way around it.**

## 2. What the prior framing got *wrong in the method's favor* (the negative is sharper)

The headline leak (90% A549 / 94.9% brain) and the "single-exon floor 8.87%" are computed over
**all eligible reads**, but the `attribute_win` function leaks a non-structural cause:

- `recovered_clip` fires on `genome_softclip ≥ 20 ∧ cdna_aligned_bases > genome` with **no junction
  requirement**. The overwhelming majority of attributed wins are `recovered_clip` on **single-exon
  reads** (`cdna_njunc=0`, anchor sentinel `1e9`), winning by **exactly 1.0 HP-ED** (the ε floor),
  several with **negative `cpa_shift`** (−9, −26 bp).
- Mechanism: the padded cDNA (1–2 kb genomic flank) swallows a 3′ soft-clip that the genome
  walkback **correctly** hard-clipped (poly-A / non-genomic). This is the walkback/CPA asymmetry
  the HANDOFF flagged for the CPA *gate* — but it also contaminates the whole-read HP-ED *scorer*
  and the `recovered_clip` *attribution*. These "wins" are not rescues; a subset are **active CPA
  regressions** (harmful for a 3′-end pipeline).

**Consequence:** the real result is not "90% noise leak" — it is "**~3 genuine rescues out of
5,644 eligible reads (0.05%)**, plus a cloud of ε-floor jitter and walkback artifacts." Same NO-GO,
cleaner reasoning. *Do not quote the leak ratio as the decision metric* (the advisor's false-GO
trap: removing the jitter could drop the ratio below 5% while structural wins stay at ~3).

## 3. What the prior framing got *wrong against the method* (the mechanism IS real)

Filtering spliced wins to those with a real HP-ED gain (≫ the 1.0 jitter floor) **and** an
increased junction count (genome demonstrably missed a junction) yields **3 clean existence proofs**,
confirmed at the CIGAR level:

| Read | Locus | Genome CIGAR (corrected, all 3 aligners agree) | cDNA → lifted | Gain |
|---|---|---|---|---|
| `8a7102d4` (A549) | TRIO | `11S 38M 2I 30M 1D 24M 2912N 19M…` — 2 junc, **11 bp 5′ clip** | no 5′ clip, **3 junc** | **+13** |
| `e89261eb` (brain) | TRIO | `9S 125M 3938N 48M…` — 4 junc, **9 bp 5′ clip** | 2 bp clip, **5 junc** | **+8** |
| `f750fc18` (brain) | MPC1 | `57S 371M 1D 7M 1D 148M 11S` — **0 junc, 57 bp 5′ clip** | **1 junc**, clip recovered | **+11** |

All three: anchor gate passes (10–11), CPA stable (shift ≤2). In each, the genome aligner
**soft-clipped a read end rather than seed across a junction**, and the pre-spliced cDNA placed the
read fully. **This is precisely the hypothesized mechanism, working.** The prior pass never looked,
so it could neither confirm nor exclude this — it is the gap Kevin named ("in-depth examination of
read alignment edge cases").

The multi-junction "wins" (equal junction count, gain ≈1–3.7) are NOT existence proofs — they are
a few recovered indel bases / jitter. The multi-junction use case is genuinely dead, as the spec
predicted ("well-powered, tested-and-weak").

## 4. The diagnostic insight — *why* the population is tiny, and which substrate would enlarge it

**Every one of the 3 real rescues is a junction near the read's 5′ end with a short anchor.** The
genome aligner clips because the terminal exon is too short to seed across the junction; the cDNA
reference removes the need to seed.

This is the bottleneck: **the rescuable read must terminate near a junction.** ONT **direct RNA is
3′-biased and 5′-truncated** — most reads never reach a 5′-proximal short exon, so the population
that *can* be rescued is structurally suppressed on DRS *regardless of how many micro-exons the
biology contains.* The brain retest confirmed this from the other side: its micro-exons were
**constitutive** (TRIO 95.8%, DIAPH1 100%, ABLIM3 100%) → minimap2 sees them in every read → no
seed-failure to rescue. Constitutive inclusion is the *wrong* test for a seed-failure mechanism.

**Therefore the data-driven next substrate is NOT another DRS cell line, NOT SMN1/SMN2 (that is
multi-mapping/paralog hardness — normal-sized exons, minimap2 seeds them fine, the cDNA round
cannot disambiguate near-identical paralogs), and NOT yeast (anchor gate off → Round 2 refuses to
run). It is full-length PCR-cDNA**, where reads reach both ends and the
"read-terminates-near-a-short-anchored-junction" population is mechanically larger. Kevin's pipeline
already handles PCR-cDNA (`ont_cdna`); a `by4742_cdna_*` or human cDNA chr5/6 set is the right
Phase-0 retest. Secondary candidate: a divergent strain/species where k-mer seeding breaks
systematically (already in `EXPLORED_ROUTES` "conditions that would justify revisiting").

## 5. Recommendation

1. **Record the reframe, keep the NO-GO for DRS.** The mechanism is real but the DRS rescuable
   population is ~0.05% — not worth the production library/lift-over machinery for DRS substrates.
2. **One more Phase-0, on full-length PCR-cDNA** (the substrate the §4 diagnosis points to), with
   two apparatus fixes:
   - **Pre-register the decision rule as the absolute count of junction-increasing wins with
     gain ≥ ~5 HP-ED on the spliced/weak subset** — never the leak ratio (false-GO trap).
   - **Fix the attribution + scorer confound:** require a junction for `recovered_clip`; score
     HP-ED raw-vs-raw or over the spliced interior only, so the walkback/CPA asymmetry stops
     contaminating both the win count and the CPA gate.
   GO only if the junction-increasing-rescue rate is materially above the DRS 0.05% floor.
3. **If PCR-cDNA also floors out → close the route for good.** At that point the bottleneck is
   intrinsic (minimap2 splice-mode already places the overwhelming majority of junction-spanning
   reads at any read terminus), and no substrate will create headroom.

## 5b. RESULT of the PCR-cDNA retest (2026-06-15) — route CLOSED

Ran the §5 recommendation: SG-NEx A549 PCR-cDNA (4,058 reads, 3 reps) vs the SAME A549 cells/loci/
42-cDNA library, only the protocol changed. Apparatus fixes applied as specified (attribution
requires a junction; pre-registered strong-junction-win decision metric). F.0 re-validated (975
reads, 0 mismatches). Workdir: `sherlock:.../rectify_round2_phase0/cdna_a549/`.

- **Strong junction-increasing wins: 4 / 3,324 eligible = 0.120%** (~3 genuine after dropping one
  with a 1,112 bp CPA shift). Controlled vs A549 dRNA: 1 → ~3 genuine reads. DRS floor was 0.05%.
- **NO-GO.** The truncation hypothesis is directionally right but the magnitude is negligible
  (0.09%, same order as the floor). CPA-regression wins rose to 57 — net-harmful for 3′ ends.
- **Root cause confirmed across 3 substrates:** minimap2 splice-mode already seeds junction-spanning
  reads at any read terminus; read completeness doesn't create a genome-unseedable pool.

Decision: **close the route** (EXPLORED_ROUTES §1 updated). Do not build Phases 1–3. Only a
divergent-seeding substrate (non-model organism/strain, gate ON) with a prior existence-proof of a
large unseedable pool would justify re-opening — not more cell-line DRS/cDNA.

## 6. What stays valid (unchanged from prior pass)
The apparatus (`liftover.py`, `cdna_models.py`, `score_phase0.py`, F.0 round-trip) is correct —
F.0 returned 0 mismatches on 833+98 real reads. The corrected-genome baseline is the right baseline.
The kill-gate design caught the dead end cheaply. The two lift bugs F.0 caught are guarded by
`test_projection_local.py`. None of that changes; this review only re-reads the *verdict*, not the
machinery.
