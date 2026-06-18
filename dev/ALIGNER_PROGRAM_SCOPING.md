# RECTIFY native-aligner program — deliverable scoping (2026-06-18)

Cycle scope (user-approved): **validated simulation benchmark (the gate) + crafter design doc only.**
Native member DEFERRED to a later cycle, gated on the benchmark. Pressure-test verdict:
`build_groundtruth_first` (result: `tasks/wy10ip5ax.output`).

Build-vs-reuse derived from repo inventory (Explore `ac2c32ec`, branch `drs-validation-rebuild`).

---

## Deliverable A — simulation-based ground-truth benchmark (THE GATE)

**Why:** every "this primitive wins read-class X" claim is unfalsifiable until truth exists. This is
also the paper thesis ("no benchmark isolates aligner-algorithm-class on ONT direct-RNA"). Fitness
function for any future native member; replaces the internal score that was provably an artifact 6
commits ago.

### Reuse (exists)
- `scripts/benchmark/aligner_contribution_analysis.py` — per-aligner win-rate cross-tab (evaluates
  existing BAMs, NOT simulated truth). Scaffolding for the report layer.
- `rectify/core/commands/validate_command.py` — 3'-end (CPA) truth scoring (`CorrectedPosition`,
  `ValidationResult`). Pattern to mirror for junctions; does NOT do junctions today.
- `rectify/core/splice/junction_validator.py` (COMPASS 3-pass) + `junction_scoring.py` motif/anchor
  classification — reusable for the SCORER's junction-call side.

### Build (absent)
1. **Read-simulator wrapper** — call pbsim3 / badread / nanosim (DRS + cDNA error models) over an
   LRGASP-style transcript GFF; emit FASTQ. (No simulator referenced anywhere in repo.)
2. **Per-read junction-truth propagation** — carry each simulated read's KNOWN donor/acceptor
   coordinates + NIC/NNC class as SAM aux tags / a sidecar truth table. (Existing `XV` tag is a
   correction-CATEGORY label, not junction truth.)
3. **Per-junction accuracy scorer** — TP/FP/FN per aligner, stratified by junction class
   (canonical / annotated / NIC / NNC) AND by read-class, against the truth table. Must use the
   ambiguity-aware match (`normalize_junction` / `_canonical_within_window`) so a correct call one
   bp into the donor/acceptor repeat is not charged as FP (the exact trap that produced GMAP 0.09).

### `benchmark_requirements` — LANDED (panel `w03wt9tmh`) → see `ALIGNER_IDEATION_SYNTHESIS.md`
Concrete spec is now in `dev/ALIGNER_IDEATION_SYNTHESIS.md` §"Benchmark requirements". Headline:
**framing metric = EXACT INDEL POSITION CONCORDANCE WITH TRUTH (not edit distance — tied by
construction).** Contain: HP runs A/C/G/T×1-12, genomic-A-abutting CPA, NIC/NNC canonical+noncanonical
novels, edit-distance-tied STR positions, paralog loci (SMN1/SMN2-style), constructed panel-failure
reads (size the discovery-ceiling tail), coverage strata + calibrated/miscalibrated phred deciles.
Meta: held-out train/test split for any table; fitness = truth set NEVER the internal score.

Open decision already surfaced to user: simulation-first NOW (chosen) vs. SIRV/sequin spike-ins
(none in repo, needs data gen; GM12878-IVT is unspliced → can bound FP, cannot validate a novel TP).

---

## Deliverable B — read-level corroboration of the ~127 GMAP-only recurrent GT-AG novels (CHEAP WIN)

**Why:** turns "recurrence + GT-AG candidates" into validated novel junctions (or refutes them);
central to the paper thesis; cheap on existing A549 data. Recurrence + GT-AG is suggestive, NOT proof.

### Reuse (exists — this is mostly assembly)
- `collect_junction_counts_from_bam()` (`junction_scoring.py:473-485`) → per-junction
  `{raw_count, anchored_count, read_ids}` with a 10bp anchor gate. The core evidence extractor.
- `normalize_junction` / `junction_ambiguity_window` / `_canonical_within_window` /
  `_normalized_annotation_set` (`chimeric_consensus.py:59-155`) — ambiguity-aware match + GT-AG.
- `load_annotated_junctions` (`consensus.py:1222+`) — to confirm the 127 are genuinely NOVEL (not a
  catalogued junction shifted into the ambiguity window).

### Build (thin)
1. **Aligner-stratified evidence wrapper** — extend/wrap `collect_junction_counts_from_bam()` to
   return `{junction: {aligner: anchored_count}}` across the 5 A549 per-aligner BAMs, so we can show
   which INDEPENDENT aligners (if any) corroborate each GMAP-only candidate.
2. **Candidate loader** — read the 127 from the Sherlock jval outputs
   (`/scratch/users/kevinroy/rectify_chimeric_measure/jval/`); normalize each to leftmost coord;
   drop any that collapse onto an annotated junction within the ambiguity window.
3. **Corroboration report** — per candidate: total anchored read support, #independent aligners,
   recurrence (#reads ≥5 gate), GT-AG within window confirmed, annotated-overlap excluded. Bucket
   into {independently-corroborated / GMAP-only-but-high-support / GMAP-only-singleton-noise}.

### Resume specifics
- jval tooling that produced the candidates lives on Sherlock
  (`/scratch/users/kevinroy/rectify_chimeric_measure/measure/junction_validate.py`); the in-repo
  equivalents above are what a clean implementation should use. A549 per-aligner BAMs are on Sherlock
  (`/scratch/users/kevinroy/rectify_human_validation/` per workspace CLAUDE.md). Run on Sherlock
  (chunked, owners partition, AVX-512 constraint) — do NOT relay BAMs through the M1.

---

## Sequencing
1. **NOW (panel-independent):** Deliverable B — thin, high-value, no truth set needed.
2. **On `w03wt9tmh` synthesis:** report verdict + concepts + `benchmark_requirements` to user (gate the
   crafter design-doc fan-out), then finalize Deliverable A's sim-transcript set + metric list folding
   in `benchmark_requirements`.
3. **Design doc** = crafter fan-out, gated and reported-before-spend (user's standing instruction).
