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

---
# ★ VETTED (v2, 2026-07-07) — 3-critic adversarial pass; must-fixes folded in
Verdict: **direction sound, but the validity assessment REVERSES by question** (the draft's blanket "short
reads are orthogonal" hid it). Q1 is the VALID one; Q2 as drafted is RIGGED.

## The reversal (why)
- **Q1 (guard / HP-drift) = VALID.** The drift is caused by homopolymer *undercount* — a DRS-specific error
  Illumina does NOT share — so Illumina breaks the tie on the correct axis. (After the 5 must-fixes.)
- **Q2 (non-canonical discovery) = RIGGED as drafted.** STAR/HISAT2/regtools share minimap2's short-anchor-
  over-GT-AG *motif-snapping* (COMPASS handoff Step E admits it); STAR `SJ.out.tab` is post-filtered on
  canonical motif. So an unstratified Q2 scores the tool's TRUE non-canonical recoveries as *fabrications*,
  against a canonical floor that can't bound the non-canonical target — re-importing the exact real-data bias
  the simulation was built to escape. Q2's fixes are MANDATORY, and it's gated on a motif-blind truth.

## Q1 — refined, tractable NOW (chr5 rep1 DRS + all 3 Illumina reps' chr5)
0. **FIRST, the load-bearing positive control (the decisive omission all 3 critics named):** measure base-exact
   Illumina split-read concordance at the ~3,007 already-enumerated ANNOTATED HP-abutting chr5 junctions
   (>=2/3 reps). If Illumina lands base-exact there → Q1 is answerable and every Q1 rate is reported RELATIVE
   to this HP-abutting floor. If not → Q1 is unanswerable for the HP subset, learned cheaply before any effort.
   (This closes the sharpest risk: a null "no harm" is otherwise confounded with "Illumina can't see the move.")
1. **Re-partition the 944 against comprehensive GENCODE v44 + RefSeq FIRST** — some are annotated in the fuller
   reference → adjudicated for free; carry only the truly-novel residual into Illumina Q1.
2. **Per-JUNCTION, not per-read:** collapse the changed reads to distinct (chrom, donor, acceptor); tally DRS
   reads per junction; classify each junction once. Report #distinct-junctions (honest denominator) + #reads.
3. **BASE-EXACT (±0) Illumina support — do NOT ambiguity-normalize** (the Nanopore window is the artifact under
   test; normalizing collapses arm-B==arm-Bguard and destroys the contrast — must-fix #1). Restrict Q1 to guard
   moves that cross OUT of the ambiguity window; pre-report the resolvable fraction of the residual.
4. **Decidable denominator as a primary output:** junctions with >=K split reads (K set by calibration + a
   matched-null false-support rate) in >=2/3 reps. Report FIX/HARM as a fraction of the DECIDABLE subset; if it
   is trivially small → declare Q1 underpowered (exhaust free levers: all 3 Illumina reps, fuller annotation —
   NOT more DRS, which is gated on generalizing standardize_chrom_name + the 0.24s/read cost).
5. **FIX = support-RATIO margin** (guard-position reads >> arm-B-position reads), not mere presence (guards vs
   the Illumina HP-artifact false-HARM). **HARM is one-directional / non-clean** (a real in-run acceptor the
   guard rightly holds presents identically) → treat HARM as an INVESTIGATION SET, not an error rate.

## Q2 — mandatory fixes; GATED on a motif-blind truth (defer)
- Truth = **motif-BLIND split counter** (raw N-CIGAR from a permissive alignment, long unique anchors that
  force gap placement by identity) cross-checked against **COMPASS-A549** — NEVER STAR `SJ.out.tab`.
- Stratify every non-canonical novel by **canonical-decoy-within-snap-window** (no-decoy = clean; decoy-present
  = confounded/non-probative).
- Floor = **annotated NON-CANONICAL** junctions (GC-AG, AT-AC minor), not canonical — measures Illumina's OWN
  non-canonical sensitivity.
- Control = **arm-A** (motif-biased refinement, SAME engine): motif-blind contribution = arm-B − arm-A, NEVER
  arm-B − raw-minimap2 (which confounds motif-blindness with the re-placement engine).
- Confirm the SG-NEx Illumina reps are the SAME cell-line/sample as DRS rep1 (else biology masquerades as fab).
- **Q2 truth-set = revive COMPASS-A549** (the PI's motif-agnostic detector), NOT STAR-2pass. So Q2 is gated on
  COMPASS being finished (smoke rc=1 currently). STAR-2pass allowed ONLY as a canonical sanity cross-check.

## Sequencing (final)
1. **Q1 now** — positive control → re-partition → per-junction base-exact adjudication on the decidable subset.
   Unblocked by COMPASS; needs none of Q2's machinery.
2. **Q2 later** — after the motif-blind counter + COMPASS-A549 revival + non-canonical floor + arm-A control.
