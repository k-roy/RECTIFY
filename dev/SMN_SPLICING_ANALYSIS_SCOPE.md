# SMN1/SMN2 differential splicing analysis — SCOPE (Sumner SMA ONT-DRS)

Date: 2026-07-17. Data: Sumner lab ONT direct-RNA (DRS) motor-neuron iPSC panel, full-genome BAMs
`/scratch/users/kevinroy/sumner_lab/full_genome_bams/*.bam` (GRCh38 / GENCODE v44). Dual purpose:
(1) the plain-face SMA biology (SMN2 exon-7 skipping), (2) a robustness check for RECTIFY at a hard
paralog locus. This is a PLAN, not results.

---

## 1. The biology in one paragraph
SMA is caused by loss of functional **SMN1**. A near-identical paralog **SMN2** differs by a critical
**c.840 C→T** in the classic **exon 7** (a translationally-silent SNP that disrupts an exonic splicing
enhancer / creates a silencer), so ~80–90% of SMN2 transcripts **skip exon 7** → unstable truncated
**SMNΔ7** protein. Disease severity tracks SMN2 copy number (more SMN2 → more full-length SMN → milder).
The canonical DRS readout is therefore **exon-7 inclusion vs skipping (PSI)**, split by paralog.

## 2. Exact coordinates (GRCh38, + strand, GENCODE v44 MANE-Select)
- **SMN1** gene chr5:70,925,030–70,953,942 (ENSG00000172062), MANE `ENST00000380707.9`
- **SMN2** gene chr5:70,049,638–70,078,522 (ENSG00000205571), MANE `ENST00000380743.9` — ~875 kb apart.
- **The critical exon** = classic SMA "exon 7" = **MANE exon 8**, the **54 bp** exon:
  - SMN2: chr5:**70,076,521–70,076,574** ; SMN1: chr5:**70,951,941–70,951,994**
  - c.840C>T lives INSIDE this exon.
- MANE flanking exons (SMN2): exon7 = 70,070,641–70,070,751 (=classic exon 6); exon9 = 70,077,019–70,077,595
  (=classic exon 8 / 3′UTR). SMN1 equivalents: exon7 70,944,658–70,944,753; exon9 70,952,439–70,953,015.
- **Diagnostic junctions** (0-based half-open intron [donor, acceptor); per paralog):
  - INCLUSION: upstream `exon7_end → crit_exon_start` + downstream `crit_exon_end → exon9_start`
  - SKIPPING : `exon7_end → exon9_start` (the intron that skips the 54 bp exon)
  - (compute exact donor/acceptor from the coords above at build time.)

## 3. Feasibility (measured, full-depth BAMs)
- **The apparent depth is an illusion of multi-mapping.** `samtools view -c` at SMN gave ~320 reads/locus,
  but primary-only (`-F 0x900`) is **~5–58 at SMN1, ~4–505 at SMN2**, and **every primary read is MAPQ 0**
  — minimap2 cannot assign the near-identical paralogs, so it emits one arbitrary primary + a secondary at
  the other locus. Position ≠ paralog here.
- **Distinct primary SMN reads per sample** (SMN1_prim+SMN2_prim): SMA 19–563 (pooled ~1,600);
  WT 21–333 (pooled ~1,280). ⇒ **Pooled group PSI is well-powered; per-sample is fine for deep samples
  (SMA_7.12_rep2=563, SMA_8.2=378, SMA_3.6=279, SMA_2713=237; WT_21.8_rep3=333, WT_4.2_rep2=235,
  WT_HB53=226) and marginal for the thin ones (SMA_2945=19, SMA_191=23, WT_3939=21).** Plan for pooled +
  per-deep-sample; report CIs.
- **Genotype signal already visible:** SMA samples skew hard to SMN2-primary (e.g. SMA_3.6 16 vs 263,
  SMA_2713 10 vs 227); WT is balanced (WT_21.8 98 vs 85). Consistent with **SMN1 deletion in SMA** — but
  MAPQ0 makes it soft; CONFIRM genotype (see gates).

## 4. Two axes — DIFFERENT treatment (hard wall between them)
### Axis A — exon-7 in/out PSI = the biology. **Refiner-FREE, from RAW full-depth alignments.**
- A read is INCLUSION if its CIGAR carries an aligned block over the 54 bp critical exon; SKIPPING if it
  has the single N-op spanning exon7_end→exon9_start. This is a large presence/absence call, robust to
  boundary drift AND to paralog mismapping (the exon block shows regardless of which locus the read landed).
- **Do NOT route PSI through the refiner** — we just proved the re-placer fabricates at paralog loci; letting
  it touch the readout risks manufacturing the effect. Use raw `sub`/full-depth CIGARs.
- Per read: classify include/skip/ambiguous; require the read to span exon7_end..exon9_start. PSI =
  include/(include+skip). Report pooled-by-group + per-deep-sample, with binomial CIs.

### Axis B — paralog assignment (SMN1 vs SMN2) = the hard axis. **Variant panel on the READ sequence.**
- Position/MAPQ is useless (all MAPQ0). Assign by paralog-distinguishing bases the read carries.
- **The c.840 trap:** c.840C>T is INSIDE the critical exon, so it is ABSENT from exactly the SKIPPED
  (SMNΔ7) transcripts we care about. Do NOT classify Δ7 reads by c.840 alone — build a **multi-marker
  panel** including **exon-8/3′UTR (MANE exon 9) paralog SNPs** and intronic paralog markers that survive
  skipping. Classify each read by a **likelihood over all markers it covers** (ONT ~5–10% error → never
  trust a single base). Track "unclassifiable" as its own bin; note that classifiability differing between
  include and skip can BIAS PSI if unhandled → compute paralog-split PSI only on confidently-assigned reads
  and report the assignment rate.

## 5. SMA-vs-WT: separate DOSAGE from REGULATION
- Total exon-7 PSI will be lower in SMA largely because SMN1 (all full-length) is gone — that is **gene
  dosage, not a splicing change.** Report it, but the **regulatory** readout is **SMN2-specific exon-7 PSI**,
  expected ~constant (~10–20%) unless something modifies it (a real differential-splicing finding).
- Stratify by **SMN2 copy number** if available. Primary comparison: SMN2-specific PSI, SMA vs WT, pooled +
  per-sample; secondary: total PSI (dosage).

## 6. Dev-validation arm (secondary; honest reframe)
- **The hard problem at SMN is paralog DISAMBIGUATION (a mapping/MAPQ0 problem), which RECTIFY's junction
  re-aligner does NOT solve** — it refines junction boundaries, not locus assignment. So SMN is a weak venue
  for demonstrating the re-aligner's headline value (that lives in the genome-wide recall win already shown).
- What SMN DOES test: (a) does the **compensating-indel fix** keep the re-placer well-behaved at a
  near-identical paralog locus (no fabricated junctions here)? (b) a self-contained **mismapping-rate**
  metric — of reads whose primary is the SMN1 locus, what fraction carry SMN2-specific bases (raw vs
  refined). If genotype = SMN1-deletion, mismapping rate = fraction of SMN1-locus reads outright.
- Re-refine ONLY the SMN locus (cheap) with the fixed refiner for this arm; keep it walled off from Axis A.

## 7. Gates & open questions (resolve before/while building)
1. **[BLOCKING for interpretation] SMN1 genotype per patient** (homozygous exon-7 deletion vs point-mutant)
   + **SMN2 copy number**. Ask Stephen Brown (sbrow208@jh.edu) / check line metadata. If SMN1-deleted →
   paralog problem collapses (all SMN2; any SMN1-locus read is definitionally mismapped = free dev truth).
2. **Cross-check** against the independent short-read **COMPASS** junctions at SMN, and against Stephen's
   **FLAIR→OARFISH** isoform calls (`sumner_lab/flair_oarfish/`) + the existing **`smn_region_heatmap_*.tsv`**
   paralog-expression work (reconcile, don't duplicate).
3. Confirm strand/orientation of DRS reads at SMN (DRS is stranded; + gene).
4. Decide ambiguous-read policy (reads spanning too little; MAPQ0 secondary handling).

## 8. Concrete first steps (in order)
1. Email Stephen for genotype + SMN2 copy number (gate #1) — do this first, it may halve the problem.
2. Build the exact junction coordinate set (§2) + the paralog-marker panel (§4B) from the two MANE
   transcripts + a SMN1/SMN2 sequence diff over the gene body.
3. Axis A prototype on the 4 deepest samples (raw full-depth): per-read include/skip caller → PSI + CI.
4. Axis B: marker-likelihood paralog classifier; report assignment rate + the c.840-trap handling.
5. Paralog-split, dosage-separated SMA-vs-WT PSI (pooled + per-deep-sample) → the headline table + figure
   (exon-7 inclusion, provenance footer).
6. Dev arm: mismapping rate raw vs refined at the re-refined SMN locus.

## Files / provenance
- BAMs: `sumner_lab/full_genome_bams/*.bam` (full depth; the downsampled `panel_deep/*.sub.bam` are TOO thin
  for SMN — do NOT use for Axis A).
- Annotation: `compass_a549/COMPASS/genome_references/GRCh38_gencode_v44.gtf` (MANE tx above).
- Prior SMN work: `sumner_lab/smn_region_heatmap_{long,matrix}.tsv`, `sumner_lab/flair_oarfish/`.
- Fixed refiner (dev arm): `/scratch/users/kevinroy/rectify_guard2` (has commit e40ca00 + `_positional_signal`).
