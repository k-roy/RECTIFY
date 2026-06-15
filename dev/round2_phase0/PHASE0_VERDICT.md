# Round-2 cDNA "discovery → assignment" — Phase-0 empirical kill-gate: **NO-GO**

**Date:** 2026-06-14 · **Substrate:** SG-NEx A549 directRNA (rep1+4+5+6, RNA002) · 7 loci
**Spec:** `dev/specs/SPEC_round2_cdna_discovery_assignment_20260614.md` (§11: this was the single
riskiest assumption to test first). **Verdict gates the whole project: NO-GO ⇒ do not build the
production library/lift-over** without changing the substrate (Kevin's call — see §5).

---

## 1. The question

Does a non-trivial subset of Round-1-weak reads get placed **strictly better** by a contiguous
aligner on a pre-spliced cDNA than by the genome consensus (including uLTRA's own record), with
the k≥10 perfect-match anchor gate ON? (Recovers exonic bases the genome soft-clipped at a junction
it couldn't seed; places micro-exons seed-aligners skip.)

## 2. What was built and validated

- **Apparatus** (`dev/round2_phase0/`, harness + `cluster/`): lift-over (transcript→genome, N-op
  insertion, minus-strand), win-guard + BLOCKER-1 raw-anchor gate + uLTRA-specific metric, all using
  the **production** `_cigar_*` scorers and the real human RNA004 penalty table. 15/15 unit tests.
- **Step D** (Round-1): 4,665 reads at the 7 loci re-aligned with minimap2 + uLTRA + deSALT, each
  run through `rectify correct` (DRS, walkback fired: 34% position shifts). **Baseline is the
  CORRECTED genome winner + corrected uLTRA** (advisor: a raw baseline is a false-GO generator —
  the motif-agnostic 5'/soft-clip rescue competes with the cDNA round's exact mechanism).
- **Step C** (library): **42 cDNAs = 7 MANE ∪ 35 Round-1 read-supported k≥10 gate-passed de-novo
  chains** (advisor: MANE-only manufactures the novelty-suppression false-GO). Per-base 1:1 identity
  validated; novel chains admitted (e.g. `DIAPH1__denovo0_s107`).
- **Step F.0** (lift validation): **833 real spliced reads round-tripped, 660(+)/173(−), 0 mismatches**
  (semantic matched-(query,ref)-pair comparison; this gate caught and forced the fix of two real
  lift/projection bugs before any verdict was trusted).

## 3. Result — NO-GO, robust against three controls

### (a) Genome-mapped population (4,579 eligible reads)
- 417 total cDNA "wins", but **only 11 are spliced** and just **7 carry an attributable structural
  cause** (1 placed-micro-exon + 5 spanned-multi-junction + 1 clip-on-spliced). These 7 are *below
  the noise floor* (next bullet) — so even the best-case wins don't clear noise; "attributed", not
  proven-genuine (not individually IGV-vetted — that's a Kevin step if ever revisited).
- **Single-exon reads are a null control** (Round-2 ≡ Round-1 by construction → cDNA should never
  win): 406 of them won (417 total − 11 spliced). Single-exon win rate ≥ **406 / 4,579 = 8.87%**
  — a *hard lower bound* (the true single-exon-eligible denominator is < 4,579, so the real rate is
  higher), already **above the spec's 5% trivial-win-leak bar.** The win-guard's false-win floor from
  alignment jitter alone exceeds the GO threshold. Leak = 90% (all reads), 36% (spliced only).
- (uLTRA-specific 100%, net HP-ED reduction +1060 — both positive, but moot given the leak.)

### (b) Genome-UNMAPPED population (844,262 reads) — Kevin's selection-bias check
The region-extracted Round-1 set was biased toward genome-mappable reads. Tested the excluded
population directly: aligned all 844k unmapped reads to the cDNA library.
- 642 map to a cDNA (TRIO 451, DIAPH1 96, …) → 182 confident (≥80% aligned, ≥85% identity).
- **0 confident SPLICED rescues.** All 182 confident placements are single-exon. The 17 spliced
  cDNA-hits are **68–83%-identity force-fits** of low-quality reads (ANXA6 25-junction hit = 74%
  identity), several genome-mappable anyway. Genome-unmapped reads are unmapped because they are
  genuinely low-quality; the pre-spliced cDNA cannot make them clean.
- **This is threshold-robust:** the rescue signature (clean read, spans junctions, genome couldn't
  seed it) would be ident ≥0.85 ∧ njunc ≥1 ∧ genome-unmapped. The 17 spliced hits cap at **0.83
  identity** — it is the *identity ceiling*, not the frac cut, that excludes them, so the absence of
  spliced rescues holds **at any coverage threshold.**

**Both populations agree: the cDNA round provides no meaningful spliced-read rescue on this substrate.**

## 3b. Why NO-GO — mechanism (baseline ladder)

Re-scored the cDNA round against three baselines (same 4,579 reads):

| Baseline | cDNA wins | per-win HP-ED gain | leak |
|---|---|---|---|
| Vanilla minimap2 (raw, single) | **34 (0.7%)** | 10.5 | 79% |
| Raw 3-aligner consensus (uncorrected) | **33** | 3.5 | 82% |
| Corrected 3-aligner consensus (the verdict) | 417 | 2.5 | 90% |

- **The cDNA round beats even plain vanilla minimap2 on only ~34/4,579 reads (0.7%), 79% noise.**
  It is NOT the case that rectify's framework is so complete it leaves nothing — there is simply
  almost nothing to rescue on this substrate.
- **The multi-aligner framework is not what absorbs the cDNA's potential (34 ≈ 33):** going from a
  single raw aligner to the raw 3-aligner consensus barely changes the win count — so the cDNA has
  no latent gains the consensus is hiding.
- **The 417 (vs corrected) is inflated, not real superiority:** `rectify correct`'s 3′ walkback
  hard-clips genome read ends (H costs 1.0/bp by design → higher HP-ED, fewer aligned bases), so the
  un-walked cDNA candidate clears the win-guard against ~380 extra reads by a hair (gain 2.5, leak
  90%). That is the walkback/CPA asymmetry, not the junction mechanism. The apples-to-apples ceiling
  is ~33–34 noise-dominated wins.
- **Root cause = substrate has no headroom-creating reads:** A549 (SRRM4-low) has almost no
  included-micro-exon or genome-unseedable-multi-junction reads, so vanilla minimap2 already handles
  essentially everything. (This is exactly why a neuronal/SRRM4-high RNA004 retest is the next test.)

## 4. Why "single-exon reads" dominate

In ONT direct-RNA, ~88% of reads at these multi-exon loci cross no splice junction (heavy 3′-bias /
truncation → reads cover one exon). "Single-exon read" = the *read's alignment* spans no junction,
not a single-exon gene. Junction-spanning reads are correctly classified spliced (junctions counted
from the lifted read, where N-ops are re-inserted).

## 5. Caveat and the decision that is Kevin's

- **Honest substrate limitation:** A549 is **SRRM4-low** (micro-exon-poor) and **RNA002** (not
  RNA004). The **multi-junction** use case *was* well-powered (TRIO 56 junctions, XPO5 31) and is
  **tested-and-weak** — that part of the NO-GO is solid. The **micro-exon** sub-case is **under-tested**.
- **Recommendation:** record **NO-GO on this substrate** and **stop** (Phase 0 = cheapest test, then
  stop). Whether to test a **neuronal / SRRM4-high RNA004 substrate** before fully killing the
  feature is **Kevin's call**, not an auto-launch. Phase 1–3 (production library/lift-over) remain
  blocked until/unless a substrate change flips Phase 0.

## 6. Artifacts
- Verdict: `<cluster>/round2/verdict_primary/{verdict.txt,verdict.json,wins.tsv}`
- Bias check: `<cluster>/round2/bias/rescue/{rescue_summary.txt,rescues.tsv}`
- Provenance + table hashes: `dev/round2_phase0/PROVENANCE.json`
- Cluster workdir: Sherlock `/oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/round2/`
