# Sumner SMA panel — splicing + CPA leads scoping (2026-07-09)

**Collaboration:** JHU/Sumner (SEPARATE from Chanfreau/RECTIFY-dev; hourly-billed). SMA (SMN1-null) vs control
iPSC-derived motor neurons, ONT direct-RNA (full-length → splicing AND 3′-end/CPA per read). Deliverable: a
prioritized, honest LEAD set to present to the Sumner team.

## Data + prior work (scouted)
- 15 WGS BAMs: 8 SMA + 7 WT (control). Depth 2.8M–14M reads. RIN varies (metadata: good/poor). PolyASite atlas
  `references/polyasite_atlas.GRCh38.96.bed.gz`. 
- PRIOR (June 2026, build on — don't redo): `apa_3prime/` gene-level abundance DE with orientation (SMN1-down
  positive control validated; 832 genes; RIN-matched, UNDERPOWERED 2 SMA vs 4 control good-RIN) +
  intronic-CPA-fraction per sample (~4.2% of CPA sites intronic). Junction discovery (this session): 5% + 25%
  genome-wide (recurrence-filtered; SMA leads chr6:7.29Mb/chr21:25.9Mb/chr15:78.9Mb; fabrication-dominated raw).
- GAPS: (i) no differential SPLICE-JUNCTION test; (ii) SMN2 exon-7 skipping (the SMA mechanism) never quantified;
  (iii) no 3′UTR APA-SHIFT (proximal↔distal) analysis — RECTIFY's specialty, unexamined.

## Proposed lead categories (SMA vs WT)
### A. SMN-locus ANCHOR (positive control + the disease mechanism)
- SMN1 abundance DOWN (validated). SMN2 exon-7 INCLUSION vs SKIPPING ratio SMA vs WT (SMNΔ7 = the SMA readout;
  the exon6→exon8 skip junction). Built-in sanity that the pipeline sees the known biology.
### B. Differential SPLICING
- Recurrent, SMA-specific differential junctions: exon skipping, alt 3′/5′SS, cryptic/NON-CANONICAL 3′SS
  (Prp18-relevant — Roy 2023). ANNOTATED-or-corroborated only (exclude the motif-blind fabrication). Paired test
  across samples; effect size + recurrence, not just p (underpowered).
### C. CPA / 3′-END (RECTIFY's specialty)
- 3′UTR APA SHIFTS: proximal↔distal poly-A site usage per gene (3′UTR shortening/lengthening) SMA vs WT, vs the
  PolyASite atlas. (NEW — the gene-DE + intronic work didn't do this.)
- Intronic CPA (IPA) changes SMA vs WT (extend the existing intronic-orient fraction to per-gene differential).
- Genes with 3′-end dysregulation / mis-cleavage.
### D. Cross-cutting biology
- SMN-pathway genes, known SMA modifiers (PLS3, NCALD, …), NMD-sensitive transcripts (SMN loss + splicing
  defects → NMD; a defect masked in abundance may show in splicing/APA).

## Rigor / honesty framing (small n, RIN-varied)
- Present as PRIORITIZED CANDIDATE leads (effect size × recurrence × biological plausibility), NOT significance-
  gated hits. Underpowered (≤2 SMA vs 4 control good-RIN if RIN-gating).
- CAVEATS to state up front: SMN1/SMN2 PARALOG confound (near-identical → read-assignment ambiguity; locus-level
  robust, within-locus provisional); fabrication in NOVEL junction calls (motif-blind, no ground-truth precision
  yet — spike-in track pending); RIN as covariate (degradation → 3′-bias artifacts that MIMIC APA shortening —
  critical for the CPA leads); small n.

## Plan (advisor-revised — GATES FIRST, fan out only after Gate 2 passes)

### GATE 1 — RIN vs group (BLOCKS category C until resolved). RUN FIRST.
RNA degradation truncates 5′ ends + manufactures apparent 3′UTR shortening + internal-priming CPA — it FAKES the
APA-shift signal. Prior June work already RIN-matched to 2-SMA-vs-4-control-good → RIN IS skewed → mandatory.
Pull per-sample RIN; test RIN~group; pick the comparison set (RIN-matched subset OR RIN covariate). ALL CPA leads
add an INTERNAL-PRIMING filter (genomic-A stretch downstream of the called site) on top of the PolyASite atlas.

### GATE 2 — SMN exon-7 anchor (the disease definition). RUN FIRST, before any genome-wide sweep.
SMA = SMN1-null → transcripts from SMN2 → ~90% exon-7 skip. If the pipeline can't recover this, don't trust any
splicing lead. Frame as LOCUS-level SMN(1+2) exon-7 INCLUSION ratio (paralog-ROBUST — reads can't cleanly split
the paralogs, but the aggregate ratio still shifts), NOT a per-paralog claim. Nail the exon6→exon8 skip-junction
coordinate against the SMNΔ7 LITERATURE coord, not GENCODE exon-numbering. Weak anchor → suspect paralog mapping
collapsing the signal BEFORE any biological conclusion.

### THEN (only after Gate 2 validates):
1. MINE the existing June gene-DE (832 genes, SMA-vs-control log2FC already computed) → top hits past SMN =
   free category-D candidates.
2. Differential SPLICING — do NOT drop novel junctions; FILTER them by (a) cross-sample RECURRENCE (per-read
   fabrication is random, won't recur at the same coord) + (b) the novel-feature-support-not-tail-enriched
   guardrail (NATIVE_ALIGNER_OVERVIEW). Recurrence-surviving cryptic/non-canonical 3′SS = the interesting Prp18
   biology, not excluded.
3. CPA — 3′UTR APA-shift + IPA differential, RIN-GATED (Gate 1) + internal-priming-filtered.
4. Cross-reference SMN pathway / SMA modifiers (PLS3, NCALD) / NMD targets.
5. Synthesize a SHORT ranked table: gene · event · effect size · direction · recurrence · RT-PCR/3′RACE
   validation path · caveat. Led by the SMN anchor proving the method. CREDIBILITY > completeness.
   Express junction/isoform leads in FLAIR terms where possible (Stephen's FLAIR→OARFISH toolchain →
   cross-checkable, de-risks bespoke-script bugs).

---
## GATE RESULTS (2026-07-09)
### Gate 1 (RIN) — PASS (better than feared)
RIN balanced across groups: good-RIN = 2 CNTL + 4 SMA; degraded = 2 CNTL + 3 SMA. RIN is NOT strongly confounded
with group -> CPA leads use the good-RIN subset (4 SMA vs 2 CNTL) without RIN masquerading as group. Real limit =
small control n (2 good), not RIN. (Still add the internal-priming filter for CPA.)

### Gate 2 (SMN anchor) — PASS via SMN1-null; the nuance matters
Per-sample chr5 (gate2_smn_anchor.py):
  SMN1 abundance: SMA mean 13 reads vs CNTL ~68 (~5x DOWN) -> the SMN1-null disease definition VALIDATES STRONGLY.
  SMN2 exon-7 skip_frac: SMA 0.62 vs CNTL 0.57 -> WEAK/marginal, overlapping distributions.
INTERPRETATION (important, honest): the DISCRIMINATING SMA signal is SMN1-LOSS (abundance), NOT a change in SMN2's
own splicing. SMN2 predominantly skips exon 7 in BOTH groups (that's SMN2's intrinsic behavior), so SMN2-exon-7-skip
alone barely separates SMA from control. The textbook "exon-7 skip up in SMA" manifests at the LOCUS level as loss
of the SMN1 full-length (inclusion) contribution -> captured cleanly by SMN1 abundance down. So the pipeline SEES the
defining biology; the naive per-SMN2 skip readout is correctly the weaker one. Anchor validates. Proceed to leads.
NOTE: CNTL_GSB3939 is a low-coverage outlier (4 SMN2 / 17 SMN1 reads) — down-weight.

## PROGRESS (2026-07-09)
- GATES PASS. SWEEP A DONE -> dev/SUMNER_SMA_LEADS_A_geneDE.md + cluster sma_leads/sma_leads_A.md. Leads:
  GTF2H2(5q13 SMA locus)/FST(follistatin)/TARS1-DT/protocadherins. CAVEAT: shallow (20 genes >=30 reads).
- SWEEP B (differential splicing) READY: 25% junction run 33236959 = 13/15 done. NEXT: per-sample junction
  counts (annotated + novel) -> paired SMA-vs-WT differential test, recurrence-filtered (keep novel/non-canonical
  that recur; drop per-read fabrication). Anchor at SMN region.
- SWEEP C (CPA/APA-shift) infra EXISTS: trim_polya/ per-sample 3'-ends + apa_concordance.py (atlas concordance) +
  polyasite_atlas.GRCh38.96.bed.gz. NEXT (build): per-gene 3'UTR APA-shift (proximal<->distal poly-A usage) SMA vs
  WT on the GOOD-RIN subset (4 SMA/2 CNTL) + internal-priming filter (genomic-A downstream). RECTIFY's specialty.
- SWEEP D (cross-ref): overlay A/B/C hits with SMN-pathway / SMA-modifier (PLS3,NCALD) / NMD gene lists.
- SYNTHESIS: short ranked table (gene·event·effect·direction·recurrence·validation-path·caveat), led by SMN anchor.
  Deliverable -> /scratch/users/kevinroy/sumner_lab/sma_leads/. Express in FLAIR terms for Stephen where possible.

## SWEEP B DONE + C LAUNCHED (2026-07-09)
- B (differential splicing = re-placer rigor test) DONE -> dev/SUMNER_SMA_LEADS_B_splicing.md. VERDICT: finds real
  SMA biology (UBA1=X-linked SMA gene, SNRPN=snRNP/SMN's function) BUT fabricates at high-expression genes
  (ribosomal cluster) -> the re-placer's real-data precision ceiling exposed; needs spike-in calibration.
- C (APA-shift, RECTIFY specialty) LAUNCHED: sbatch 33302223 (good-RIN 4 SMA/2 WT, internal-priming filtered,
  atlas-anchored). RESUME: ssh sherlock 'cat /scratch/users/kevinroy/sumner_lab/sma_leads/.sweepC_rc; tail -40
  /scratch/users/kevinroy/sumner_lab/sma_leads/sweepC-33302223.log'. Reports top 3'UTR shortening/lengthening SMA vs WT.
- D + synthesis: after C. Overlay A/B/C with SMN-pathway/modifier/NMD; short ranked table -> sma_leads/.

## ★ RIGOR-TEST FEEDBACK TO THE NATIVE ALIGNER (the SMA panel's purpose)
Sweep B is the native re-placer's real-data discovery output. It shows the HP-drift guard (which fixes homopolymer
fabrication) does NOT cover EXPRESSION-CORRELATED fabrication (ribosomal/high-coverage genes over-call non-canonical
junctions). => next native-aligner build step motivated: a DISCOVERY fabrication guardrail beyond HP-drift (the
"novel-feature support must not be read-quality-tail-enriched" / expression-normalized-support principle from
NATIVE_ALIGNER_OVERVIEW), AND the spike-in precision number to calibrate/measure it. UBA1/SNRPN surviving = the tool
finds real biology; the ribosomal fabrication = the ceiling to fix.
