# SMN1/SMN2 differential splicing — SCOPE **(REVISED 2026-07-17 after provenance dig)**

> ⚠️ **CORRECTION TO THE FIRST DRAFT OF THIS DOC.** The first version scoped this as new work. It is not:
> a **first pass was completed 2026-06-12** and it already independently converged on the same design
> (exon-8 PSV avoiding the c.840 trap). **Do not re-derive it.** Canonical prior artifacts:
> - `~/work/JHU/sumner_lab/smn_exon7/RESULTS_smn_exon7_20260612.{md,tsv}` — the results
> - `smn_exon7/smn_exon7.py` (PSI+PSV caller), `smn1_genotype.py` (2-PSV noise check),
>   `smn_compare_arms.py` (aligner-arm comparison)
> - `~/work/JHU/sumner_lab/SUMNER_SMN_STATUS.md` §6.A — status + its own NEXT list
> - `~/work/JHU/sumner_lab/SAMPLE_PROVENANCE_coriell.md`, `sample_sheet_genomewide.tsv` — genotypes
> - Sherlock: `/scratch/users/kevinroy/sumner_lab/{smn_exon7/,gate2_smn_anchor.py,smn_region_heatmap_*.tsv}`
> This doc is now **only** "what genuinely remains."

---

## 1. GATE #1 (SMN1 genotype) — **RESOLVED**, no need to ask Stephen
Two independent lines of evidence agree.

**(a) Documented (Coriell-verified, `SAMPLE_PROVENANCE_coriell.md`):**
- **SMA_7.12** (GM09677) — SMN1 **0, homozygous EX7-8 deletion**, SMN2=3
- **SMA_3.6** (GM03813) — SMN1 **0, homozygous EX7-8 deletion**, SMN2=3
- **WT_4.2** (GM03814) — SMN1 **1 (carrier)**, SMN2=5; the **carrier MOTHER of SMA_3.6** (related pair!)
- **WT_21.8** (GM02183) — genetically normal (HTT 33/18); no SMN record (HD subcollection)
- SMN2 copy number known for all 7 SMA: **Type I = 2–3, Type II = 3**

**(b) Empirical, from our own DRS via the exon-8 PSV** (`RESULTS_smn_exon7_20260612.tsv`) — SMN1-assigned
reads stratify into three clean tiers that MATCH documented copy number (a dose-response that both
validates the assignment method AND extends genotype to the undocumented Greenstone/Iperian lines):

| tier | SMN1 PSV reads | samples |
|---|---|---|
| **0 copies (SMN1-null)** | 0–3 | **all 7 SMA** (2 Coriell-confirmed) |
| **1 copy (carrier)** | 13–14 | WT_4.2 (confirmed carrier), **WT_3939** (corroborates its fingerprint carrier flag) |
| **2 copies (normal)** | 91–133 | WT_21.8, WT_HB53 |

⇒ **SMA = SMN1-null.** The paralog problem collapses for SMA samples (every SMN-locus read is SMN2, and
any SMN1-locus read is definitionally mismapped = free dev-validation truth). It does NOT collapse for
controls (they carry both) — the PSV assigner is still required there.

## 2. The biology headline — **already measured** (first pass, pooled)
From `RESULTS_smn_exon7_20260612`:
- **SMN1 PSI ≈ 98–100% vs SMN2 PSI ≈ 45–53%** — the c.840/C6T effect, cleanly reproduced.
- **SMN2 inclusion is IDENTICAL in SMA (48%) and control (49%)** ⇒ an **intrinsic SMN2 property, not
  dysregulated in SMA**. This is exactly the dosage-vs-regulation separation: the SMA deficit is **gene
  dosage** (SMN1 gone), *not* altered SMN2 exon-7 regulation.
- Method validation: PSV specificity ≤0.4% (SMA_8.2: 247 SMN2 vs 1 SMN1); PSI identical across
  rectified/uLTRA/deSALT arms; ambig=0.
- Caveats carried forward: **pooled only** (not powered per-sample / Type I vs II); **CNTL_HB53
  lower-confidence** (degraded, SMN1 PSI 76% = base-call leakage); observed ~48% inclusion vs textbook
  ~10–15% is expected (NMD depletes Δ7 from steady-state poly-A+; the metric undercounts inclusion →
  conservative).

## 3. What genuinely REMAINS (the actual work)
Ranked by value. Items 1–2 are new to this session; 3–5 are the prior pass's own NEXT list.

1. **★ Re-run on FULL-DEPTH BAMs with the FIXED refiner (the new, highest-value item).** The first pass
   used chr5 / 3-aligner-consensus BAMs — produced **before** the compensating-indel fix (e40ca00). Two
   questions only we can now answer: (a) does the fix change any SMN PSI number? (Prediction: **no** —
   PSI was already identical across aligner arms, and exon in/out is a large presence/absence call robust
   to boundary drift. Confirming this is a **clean, honest robustness check** of the fix at the hardest
   paralog locus.) (b) full depth may lift the "pooled only" limit → **per-sample and Type I vs II PSI**.
   Full-depth distinct primary SMN reads: SMA 19–563, WT 21–333 (pooled ~1,600 SMA / ~1,280 WT).
2. **Mismapping-rate metric (dev arm, now trivially available).** Since SMA is SMN1-null, ANY SMA read
   whose primary lands at the SMN1 locus is definitionally mismapped. Report that rate **raw vs refined**
   — a self-contained RECTIFY metric needing no external truth. (NB all SMN primaries are **MAPQ 0** —
   position never identifies paralog here; keep the multimappers, assign by PSV, never by MAPQ.)
3. Realign locus reads to a **single SMN2 target** to collapse coordinates and resolve the few
   cross-placed reads.
4. Add **c.840 corroboration** for the degraded **CNTL_HB53** only (`smn1_genotype.py` already does this
   2-PSV check).
5. **Absolute-quantify SMN2 copy dosage** with depth normalization (ties PSI to the severity modifier).

## 4. Design constraints that still bind (do not re-litigate)
- **Exon-7 in/out PSI must stay refiner-FREE** for the biology readout (the re-placer must not be able to
  manufacture the effect); the refiner appears only in the dev arm (item 2) and as the robustness check (1).
- **Paralog assignment = exon-8 PSV `chr5:70,077,254` (A=SMN2 / G=SMN1)**, NOT c.840 — c.840 is inside
  exon 7 and therefore ABSENT from the skipped (Δ7) reads we most need to count; using it would silently
  drop skipped reads and inflate inclusion. Use c.840 only as a corroborating second PSV on reads that
  include exon 7 (phasing both on one molecule is a DRS advantage).
- **Relatedness:** WT_4.2 is the mother of SMA_3.6 — encode as related; not independent samples.
- **Confounds:** age ↔ condition (controls adult, Type I toddlers); repository/batch (Coriell/Greenstone/
  Iperian/collaborator); WT_4.2 is a carrier with SMN2=5, so it is not a clean "normal" control.

## 5. Coordinates (verified, GRCh38 / GENCODE v44 MANE)
- SMN1 chr5:70,925,030–70,953,942 (`ENST00000380707.9`); SMN2 chr5:70,049,638–70,078,522
  (`ENST00000380743.9`); ~875 kb apart, + strand, ~99.9% identical.
- **Classic SMA "exon 7" = MANE exon 8 = the 54 bp exon**: SMN2 **70,076,521–70,076,574**;
  SMN1 **70,951,941–70,951,994**. Flanks (SMN2): classic exon 6 = 70,070,641–70,070,751;
  classic exon 8 = 70,077,019–70,077,595.
- Junction logic (validated in `gate2_smn_anchor.py`): for an intron ending at the terminal acceptor
  (SMN2 0-based **70,077,018**), donor == **70,070,751** ⇒ **SKIP**; donor > that ⇒ **INCLUSION**.
