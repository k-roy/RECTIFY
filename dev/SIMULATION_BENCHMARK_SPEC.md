# RECTIFY simulation benchmark — spec (Deliverable A, the GATE) — 2026-06-18

Turns the ideation panel's `benchmark_requirements` (`ALIGNER_IDEATION_SYNTHESIS.md`) into a buildable
plan. This is the **fitness function** for the native-aligner program: every concept is a
hypothesis-pending-ablation against THIS truth set, **never** the internal score (which was provably
artifact-prone — the 0.09→1.07 flip). Status: DRAFT (pending advisor check + crafter `benchmark_coupling`).

## Framing metric (the one that matters)
**EXACT INDEL-POSITION CONCORDANCE WITH TRUTH, not edit distance.** At every contested position edit
distance is tied by construction, so it cannot separate any of the 5 concepts — only "which tied
placement matches truth" does. The whole design is built around position-exact truth, scored with the
**ambiguity-aware match** (`normalize_junction` / `_canonical_within_window`, `chimeric_consensus.py:59-155`)
so a correct call one bp into a donor/acceptor repeat is not charged FP.

## Two-tier architecture
The requirements split cleanly into "discriminate the concepts" (needs controlled, by-construction truth)
and "external validity / size the tail" (needs realistic whole-transcriptome reads). Build both.

### Tier 1 — CONTROLLED micro-benchmark (the discriminating tier; in-repo, no external WGS sim)
Hand-construct reference mini-loci with KNOWN truth per failure-mode stratum, then inject ONT-realistic
errors with a calibrated error model (badread's ONT model, or RECTIFY's own `empirical_cigar_error_profiler`
table for self-consistency — but the model used to SCORE must be HELD OUT from the table used to CALIBRATE
the HP-length-law concept, else the win is overfitting). Truth is exact by construction. This tier is where
position-exact concordance is scored; it needs no genome-scale simulator.

Strata (one generator each, parameterized, with a per-read truth row):
| Stratum | Construct | Truth recorded | Scores which concept(s) |
| --- | --- | --- | --- |
| HP runs | A/C/G/T × len 1-12, flanked, del-dominant length-dependent miscall | true run length, boundary positions | HP-length-law DP |
| Genomic-A CPA | true cleavage site abutting a genomic A-tract of varied length | true CPA coord, genomic-A length | CPA decoder |
| NIC/NNC novels | junctions NOT in GTF, BOTH GT-AG and non-canonical | donor/acceptor, canonical?, NIC vs NNC | discovery + FDR guard, LR arbitration |
| Edit-distance-tied STR | tandem repeats where ≥2 indel placements are ED-tied | the true placement | tie-break vs left-normalization |
| Paralog loci | SMN1/SMN2-style near-duplicates, known true locus | true locus per read | window-selection, POA mis-clustering |
| Panel-failure | reads constructed to defeat seed-chain (high divergence / chimeric) | true origin | FracMinHash fallback + TAIL SIZE |
| Coverage × Q | singleton→deep pools; Q deciles, BOTH calibrated AND deliberately-miscalibrated phred | true bases, true Q | POA pooling, soft-decision (isolate recalibration gain) |

### Tier 2 — REALISTIC transcriptome simulation (external validity + global recall/FDR)
Whole-transcriptome ONT reads with per-read origin → junction truth, for global novel-junction recall+FDR
and to SIZE the panel-failure tail on realistic data (open question Q1).
- **Simulator: NanoSim transcriptome mode** (trained ONT model; emits reads from a reference transcriptome
  with known transcript-of-origin → derivable junction truth) is the lead. Alternative/cross-check:
  **badread** (explicit, tunable ONT error model; better for controlled error injection but no native
  transcriptome-truth mode — would need a truth-propagation wrapper). Evaluate both; do NOT hand-roll a
  read simulator (none in repo today — confirmed).
- Organisms: **yeast** (simple introns, the saturation control) + a **human locus panel** incl. SMN1/SMN2
  and a NIC/NNC-rich set (mirror the A549 chr5 in-hand data).
- DRS **and** PCR-cDNA error profiles (RECTIFY serves both).

## What to BUILD (none exist in repo — confirmed by inventory)
1. **`sim/` generators** — Tier-1 stratum constructors (deterministic, seeded) emitting FASTQ + a per-read
   truth table; Tier-2 NanoSim/badread wrappers + a truth-propagation step (transcript origin → per-read
   junction/CPA truth).
2. **Truth format** — a per-read truth table (parquet/TSV) keyed by `read_id`: `true_locus`,
   `true_junctions` (list of donor/acceptor + canonical? + annotated/NIC/NNC), `true_cpa`,
   `boundary_true_run_lengths`, `stratum`. Optionally also as BAM aux tags for in-pipeline carry.
3. **`benchmark/` scorer** — consumes an aligner/consensus BAM + the truth table → the position-exact
   metrics below. Mirror `validate_command.py`'s pattern (it does this for 3' ends already); reuse
   `collect_junction_counts_from_bam` (`junction_scoring.py:473`) + the ambiguity helpers for the
   junction-call side; reuse `junction_validator` motif classification for canonical/non-canonical.

## What to MEASURE (per aligner AND per consensus config)
- exact-indel-position accuracy at HP boundaries; corrected-run-length accuracy
- |est − true CPA| stratified by genomic-A abutment (vs heuristic walkback, vs each panel end)
- novel-junction **RECALL and FDR**, stratified annotated/de-novo **×** canonical/non-canonical
  (this is the CO-PRIMARY-discovery anti-minisplice guard — FDR must be reported, not hidden)
- tied-indel placement agreement on STR; window-selection accuracy on paralogs
- fallback recovery rate vs a random-window null **+ the size of the panel-failure tail**
- consensus accuracy under integer-score **vs** calibrated-LR (replicate 0.09→1.07; show LR holds)

## Two META-requirements every ablation inherits
1. **Held-out train/test split** for any calibrated table (HP-length-law, CPA) — else the win is
   overfitting, not calibration. Re-estimate per chemistry/basecaller version (silent maintenance dep).
2. **Fitness = THIS truth set, NEVER the internal score.**

## Where it runs
Tier-1 generation is light (M1-OK). Tier-2 NanoSim/badread + the 5-aligner panel over simulated reads is
heavy → **Sherlock, chunked/owners/AVX-512** (or H2). Do not relay simulated BAMs through the M1.

## Reconciliation with crafter `benchmark_coupling` (design-side view, folded 2026-06-18)
The design doc's per-concept ablations impose extra requirements my draft missed — fold these in:
- **min cell count = 100 per (run_length, base) and per (STR unit, n_copies):** below `min_count=100`
  (`hp_penalty.py:184`) the length-law silently nullifies to flat → a FALSE REFUTE. Tier-1 HP/STR
  generators MUST size every cell ≥100, including clean **no-error** runs (the FP side).
- **ONE shared genomic-region-disjoint TRAIN/TEST split tag**, used by EVERY calibration table (HP-length-law,
  CPA, Q-recalibration), so a win cannot leak across concepts. Make it a first-class column in the truth table.
- **NIC/NNC junction truth labels required from P0** — without them the co-primary discovery FDR is not
  measurable at all. canonical/non-canonical FDR are SEPARATE tracks, never pooled.
- **Oversample** the tail-over-genomic-A-tract class (C2) and the guaranteed-panel-unplaced population (C5,
  injected at elevated error/repeat/novel-junction contexts, origin windows retained — its measured size is
  the C5 dep-commit gate).
- **Per-read which-aligners-placed labels** + the 5 panel raw ends/walkback/atract estimates carried in the
  truth table so all comparators score identical reads; C3 needs a flat-affine-tie aligner family + the
  0.09→1.07 re-weighting replay metadata at fixed placements.
- **Precomputed ambiguity-equivalence SETS per truth position** so the position-exact metric is ambiguity-aware.

## Open decisions (carry to user / fold crafter `benchmark_coupling`)
- NanoSim vs badread for Tier-2 (recommend: NanoSim primary for transcriptome-truth, badread for controlled
  error cross-check).
- Whether Tier-1's error model is badread's ONT model or RECTIFY's own empirical table (held-out either way).
- Human locus-panel composition (SMN1/SMN2 + how many NIC/NNC-rich loci; mirror A549 chr5).
