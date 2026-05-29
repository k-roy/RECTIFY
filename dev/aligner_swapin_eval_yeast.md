# Aligner swap-in eval: winnowmap2 / minisplice_mm2 vs mapPacBio on yeast upf1Δ DRS

**Date:** 2026-05-26
**Author:** aligner-eval agent (standalone `align`-stage run)
**Branch/HEAD:** `drs-validation-rebuild` @ `1190129` (= `c52966c` + branding-only commit; no code delta vs the HEAD Sherlock is synced to)
**Question:** Are **winnowmap2** and/or **minisplice_mm2** good swap-ins for **mapPacBio** in the *yeast* rectify aligner panel? (mapPacBio is being dropped from the *human*-DRS panel — see `docs/aligners/mapPacBio.md` — but is "kept for yeast" there. This eval tests whether the two new aligners would be better yeast panel members.)

---

## CRITICAL CAVEAT — read before the recommendation

`vi2-7k.kan` is the upstream **VERTEBRATE (7k) ship model**. Applying it to *yeast*
splice sites is an **organism mismatch**: yeast introns are short with a strong
`GTATGT` 5′SS / `TACTAAC` branchpoint / `YAG` 3′SS — structurally unlike vertebrate
introns. minisplice's deep-learning splice scores are therefore expected to be
**poorly calibrated on yeast**, so minisplice-on-yeast results must be read
cautiously. **winnowmap2** carries **no splice model** (it is organism-agnostic
weighted-minimizer minimap2), so it is the *a priori* cleaner candidate for a yeast
swap-in. This caveat governs the recommendation below.

**Mechanism (why the mismatch can actively hurt, not just fail to help):**
`run_minisplice_mm2` invokes `minimap2 --spsc=<scores> -ax splice -uf -k14 -G 5000
--splice-flank=no ...`. The `--splice-flank=no` flag **disables minimap2's built-in
GT-AG canonical-dimer prior** and hands junction biology entirely to the supplied
splice-score model. On vertebrate data the trained model is a fine substitute; on
yeast, if the vertebrate scores under-rate canonical yeast `GTATGT…YAG` introns,
minisplice is left with *neither* the canonical rule *nor* a calibrated model — so it
can produce *worse* junctions than plain minimap2 (which keeps its GT-AG prior).
winnowmap2 uses the standard `-ax splice -uf -k14` with NO `--spsc` and NO
`--splice-flank=no`, so it retains minimap2's canonical prior intact.

---

## Why this was run standalone (not via run-all)

winnowmap2 and minisplice_mm2 are wired into `rectify align` but **NOT** into
`run-all` (run-all's `--base-aligners` choices = `[minimap2, mapPacBio, gapmm2,
bbmap, bwa]` only, with no `--minisplice-model` passthrough). So the `align` stage
was driven directly. mapPacBio params in the in-tree `run_map_pacbio` are the
**FIXED** values (`intronlen=10`, `maxindel=max(200000, max_intron)`) verified at
this HEAD, so for yeast (`--max-intron 5000`) mapPacBio correctly relabels all
yeast introns (50–1000 bp) as N-ops — **the head-to-head is fair** (not the
near-zero-N misconfiguration that bit the human run before the 2026-05-25 fix).

## Inputs & method

- **Sample:** upf1Δ rep1, 5% subset (`$SCRATCH/wt_upf1_manifest_test_20260523/upf1d_rep1.subset.bam`).
  upf1Δ is NMD-deficient → enriched for novel/alternative intron usage → a strong splice-aligner stressor.
- **Sanity sample:** wt by4742 rep1 trimmed FASTQ (`$SCRATCH/rectify_wt_by4742_rep1_25846844_0/drs_trim/wt_by4742_rep1.subset_trimmed.fastq.gz`).
- **Prep (production-faithful):** `rectify trim-polya` BAM mode → `samtools fastq` → trimmed FASTQ (mirrors `single_sample.py` DRS Step 0).
- **Align panel:** `--aligners minimap2 gapmm2 mapPacBio winnowmap2 minisplice_mm2 --junction-aligners uLTRA deSALT --Scer --minisplice-model vi2-7k.kan --minisplice-model-cali vi2-7k.kan.cali --max-intron 5000`.
- **Genome:** bundled R64-5-1 (uncompressed `.fsa`, explicit `--genome` to keep meryl/minisplice happy); caches written to the run output dir.
- **Two align variants:** `--no-consensus` (per-aligner BAMs for intron analysis) and default (consensus → `consensus_aligner_stats.tsv` for HP-ED wins).
- **Metrics:** per-aligner N-op intron counts; annotated (±5 bp vs the 378 bundled yeast GFF `intron` rows) vs novel; per-read unique-to-aligner introns (the mapPacBio.md disagreement-rate methodology); pairwise distinct-junction Jaccard; terminal soft-clip rates; HP-ED consensus wins; runtime.

**Output dir (Sherlock):** `$SCRATCH/aligner_swapin_eval_yeast_20260526/`
**Submit job:** 26225965

---

## RESULTS

_(filled in after job 26225965 completes — see tables below)_

### Per-aligner BAM paths

- _pending_

### Intron / N-op comparison (upf1Δ)

- _pending_

### Unique-to-aligner disagreement rate (the mapPacBio.md signal)

- _pending_

### HP-ED consensus wins per aligner

- _pending_

### Runtime per aligner

- _pending_

---

## RECOMMENDATION

_(pending results)_
