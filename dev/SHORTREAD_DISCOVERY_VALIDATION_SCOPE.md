# Short-read discovery validation — SCOPE (2026-07-07)

**Goal.** Use an orthogonal short-read (Illumina) view of which junctions are *real* to answer the
two questions annotation-truth could NOT on the human chr5 DRS run:
- **Q1 (guard on novels) — tractable now.** The guard changed ~944 chr5 reads at *non-annotated* junctions
  that GENCODE couldn't score. Are those changes **fixes** (guard moves the junction to a short-read-supported
  position that arm-B had drifted off) or **harms**? This directly extends the do-no-harm / drift-fix result
  into the novel-junction domain.
- **Q2 (motif-blind discovery) — the make-or-break, gated on a sensitive truth.** Does the motif-blind
  re-aligner **recover real novel junctions that minimap2 flattens** (high recall), **without fabricating**
  (high precision)? This is the whole tool's reason to exist.

---

## What's on the cluster (inventory, 2026-07-07)
- **Re-aligner outputs to validate:** `/scratch/users/kevinroy/human_drs_out/real_arm_B.bam` (motif-blind),
  `real_arm_Bguard.bam` (motif-blind + guard). chr5, A549 direct-RNA rep1.
- **Short-read data:** `/scratch/users/kevinroy/sgnex_a549_illumina/replicate{1,3,5}/SGNex_A549_Illumina_replicateN_run1.bam`
  (genome BAMs, GRCh38, RF/dUTP stranded). Prior `illumina_validation.json` (the 14.5%-sensitivity STEP-0).
- **COMPASS A549 sensitive detector:** BUILT but **NOT finished** — smoke test `SMOKE_DONE_RC=1` (failed);
  only `compass_a549/COMPASS/COMPASS_junctions/A549_rep1_subsampled_COMPASS_splice_junctions.tsv` (100k subsample)
  exists. Open issues (see dev/COMPASS_HANDOFF.md): compare_splice OOM at human scale, BBMap params. Reviving
  it = a substantial separate track.
- **Target novel sets:** `/scratch/users/kevinroy/deliverable_b/gmap_corroboration.json` (the 111 GMAP-only
  recurrent novels + 609 corroborated). `regtools` NOT on PATH (source/install, or use a pysam split-read counter).

---

## The load-bearing decision: which short-read TRUTH?
Sensitivity vs precision vs effort. The prior 14.5% number is why COMPASS-A549 was started.

| Truth source | Sensitivity to novels | Effort | Good for |
|---|---|---|---|
| **regtools / pysam split-reads on the existing Illumina BAMs** | LOW (~1-pass, ~14.5% on known) | hours | **Q1** (support at a *specific* placement is a coverage question, not a detection-sensitivity one) |
| **STAR 2-pass on the Illumina FASTQ** | MODERATE (novel SJ discovery) | ~1 day (align 3 reps) | a first **Q2** pass |
| **COMPASS A549 6-aligner detector** (the PI's purpose-built sensitive tool) | HIGH | days (revive P2–P5) | the gold-standard **Q2** |

**Key insight that de-risks Q1:** judging whether the guard moved a junction *to* a short-read-supported
position needs only **coverage at that exact position**, not sensitive *de novo* novel detection — so the low
Illumina sensitivity does NOT block Q1. Q1 is tractable with the existing BAMs + a lightweight split-read
counter. Q2's recall estimate IS sensitivity-limited → needs STAR-2pass or COMPASS.

---

## Q1 pipeline (PRIMARY — tractable now)
1. **Short-read junction support:** from each Illumina rep BAM, count split (N-cigar) reads per (chrom, donor,
   acceptor) with adequate anchor (≥8 bp) on both flanks (pysam, ~1 pass). Keep junctions with ≥K reads in
   ≥2/3 reps → `sr_support` set (position-exact + a ±ambiguity-window variant).
2. **The guard's non-annotated changes:** from `real_arm_B.bam` vs `real_arm_Bguard.bam`, take the reads whose
   junction differs AND is NOT in GENCODE (the ~944). For each, get arm-B's placement and arm-Bguard's placement.
3. **Classify against `sr_support`:** guard placement supported & arm-B placement not → **FIX**; reverse →
   **HARM**; both or neither supported → **inconclusive** (low SR coverage / genuinely novel-unseen).
4. **Decision rule:** FIX >> HARM (mirroring yeast 17/0 and human's 0-harm on annotated) → the guard's
   drift-fix holds in the novel-junction domain. HARM > 0 → investigate (a real novel acceptor inside a run —
   the undecidable case the guard is designed to hold conservatively).

## Q2 pipeline (make-or-break — needs the sensitivity decision first)
1. Build a **sensitive novel-junction truth** (STAR-2pass or COMPASS) → short-read-supported junctions,
   split into ANNOTATED vs NOVEL (non-GENCODE).
2. **Recall:** of short-read-supported NOVEL junctions in expressed loci, what fraction does the motif-blind
   re-aligner place (vs minimap2 snapping them to canonical)? Compare arm-B (re-aligner) vs raw minimap2.
3. **Precision:** of the NOVEL junctions the re-aligner calls, what fraction have short-read support (vs
   fabrication)? Compare arm-B vs arm-Bguard (does the guard cut fabrication without cutting recall?).
4. **Calibrated controls (from COMPASS_HANDOFF):** positives = annotated junctions in the same expressed loci
   (must validate HIGH — the 14.5% floor is the yardstick); negatives = a shuffled-junction null of matched
   intron length (fabrication baseline). Frame every rate RELATIVE to the positive-control floor.

---

## Traps to honor (from prior COMPASS work + the guard's own logic)
- **Strand:** SG-NEx is RF/dUTP (0.999) — matches COMPASS default; no flip, but enforce strand-aware junctions.
- **Ambiguity:** score junctions ambiguity-aware (a call one base into a donor/acceptor repeat, or an indel
  inside a homopolymer, is NOT an error) — the same normalize_junction the DRS side uses.
- **The undecidable case:** a *real* non-canonical acceptor that genuinely sits inside a homopolymer is
  unrecoverable by anyone; if the guard "harms" one of these, that's correct conservative behavior, not a bug.
  Q1's HARM bucket must be read in that light.
- **Illumina sensitivity floor:** never report a bare Q2 recall — always relative to the annotated-positive
  floor in the same loci, or the claim is uninterpretable (this is exactly the 14.5% lesson).
- **Chrom names:** Illumina BAMs may be Ensembl `5` vs our `chr5` — harmonize (the DRS run's reheader trick).

---

## Recommended sequencing
1. **Q1 now** — lightweight split-read counter on the 3 Illumina reps + classify the guard's ~944 non-annotated
   changes. Days, existing data, directly extends the validated do-no-harm result into novels.
2. **Q2 truth decision (PI):** STAR-2pass (fast first pass) vs revive COMPASS-A549 (gold standard) — see the
   AskUserQuestion. Then run Q2 recall/precision with calibrated controls.
