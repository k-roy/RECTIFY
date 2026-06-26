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

## Two distinct validity claims (do NOT let one stand in for the other)
Simulation validates **placement MECHANICS** — did the DP put the indel/junction where truth says — but
it CANNOT validate that badread/NanoSim's error distribution matches real RNA004. Tuning the
HP-length-law / emission against simulated errors and trusting a green number = **hill-climbing into the
simulator's error model** (the same artifact trap as the internal score, one level up). Therefore:
- **Simulation = placement-mechanics truth** (this benchmark).
- **Error-model REALISM = real-data corroboration** (Deliverable B's read-level GMAP corroboration; future
  SIRV/sequin spike-ins). A simulation win is necessary, not sufficient — it must transfer to real data.
Keep the two claims separate in every report; never let a simulation number substitute for transfer.

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
| **Standing variants (C6)** | reads carrying KNOWN SNPs/indels at known positions: het vs hom, near vs far from junctions/CPA, at non-Mendelian VAFs (mimic aneuploid A549) | true variant set per locus + true junction/CPA truth | C6 variant/haplotype-aware alignment; measures variant-induced discovery FDR |

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
- **variant-induced junction/indel FDR with vs without C6** (variant-aware emission), stratified
  het/hom and variant-adjacent vs variant-distant junctions — the measurement that decides whether C6 is
  material or descoped (a variant near a splice site can fabricate a "novel" non-canonical junction →
  this is a first-class CO-PRIMARY-discovery FDR guard, complementary to the §8 abstain band)

## Two META-requirements every ablation inherits
1. **Held-out train/test split** for any calibrated table (HP-length-law, CPA) — else the win is
   overfitting, not calibration. Re-estimate per chemistry/basecaller version (silent maintenance dep).
2. **Fitness = THIS truth set, NEVER the internal score.**

## Where it runs
Tier-1 generation is light (M1-OK). Tier-2 NanoSim/badread + the 5-aligner panel over simulated reads is
heavy → **Sherlock, chunked/owners/AVX-512** (or H2). Do not relay simulated BAMs through the M1.

## ⚠ VERTICAL-SLICE FINDING (2026-06-18) — the HP stratum MUST be discriminating, v1 was not
The thin slice (`dev/bench/hp_vertical_slice.py`) ran end-to-end and exposed a real spec gap: an
**isolated, cleanly-flanked HP run is NON-DISCRIMINATING** — **BOTH** minimap2 AND the live flat-affine DP
(`align_exon_block_global`, the one C1 upgrades) scored **1.000** run-length accuracy AND position-exact
concordance on every length cell (L=1..12, n=480 each). Reason: when a read genuinely has L−k bases,
faithfully aligning it to the ref(L) yields a k-gap somewhere in the run, and any in-run placement is
ambiguity-equivalent → every aligner is trivially correct. Both arms saturating identically PROVES the
harness works AND that the stratum cannot separate flat-affine from anything — including C1. **A benchmark where the
incumbent scores 100% cannot separate the concepts.**
Implication — the HP stratum (and C1's claimed value) lives at the HARD boundary cases, NOT isolated runs:
- **indel-vs-substitution tradeoff at run boundaries** (where `homo_mismatch=−2` and the length-law DIVERGE
  from flat affine — a boundary substitution can be scored as a mismatch OR absorbed as an indel-shift);
- **run-bleeds-into-flank ambiguity** (flank shares the run base → the run LENGTH itself is ambiguous and
  the length-law prior must pick it);
- **adjacent runs / cross-boundary placement** (left-alignment lands in the wrong run);
- **combined background noise** (substitutions co-occurring with the length error).
v2 of the generator must construct these. A candidate ablation only counts if minimap2 is BELOW ceiling on
it. This may also re-prioritize C1: its isolated-run benefit looks small; the action is at boundaries/noise.

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

## SIMULATOR DECISION (2026-06-26 — RESOLVED, advisor-checked) → pbsim3 for Tier-2
The Tier-2 simulator is **pbsim3** (bioconda 3.0.5, installed into Sherlock env
`aligner_bench`). **Decisive reason: of pbsim3 / badread / nanosim, only pbsim3
emits a per-read MAF (read↔template ground-truth alignment).** The framing metric
is EXACT INDEL-POSITION CONCORDANCE, which needs a per-read ground-truth edit
script; badread and nanosim report transcript-of-origin + identity but NO per-read
edit, so they can validate junction truth (construction-derived) but NOT exact-indel
truth on realistic reads. pbsim3's MAF composed with the known transcript↔genome
exon structure (`TranscriptModel.transcript_pos_to_genome`) yields the exact
read↔genome alignment — indels AND junctions. DRS = `--method errhmm --errhmm
ERRHMM-ONT.model`; PCR-cDNA = `ERRHMM-ONT-HQ.model`.
- **Tier-1 stays self-injected controlled errors** (`scripts/benchmark/sim/controlled.py`,
  generalizing `dev/bench/hp_vertical_slice.py`) — truth by construction, no
  simulator, exact per-position edit. This is the DISCRIMINATING tier; pbsim3 is
  for realism / global recall+FDR / tail-sizing.
- Implementation: `scripts/benchmark/sim/pbsim3_wrapper.py` (MAF parse + projection,
  validated locally on a synthetic MAF); scorer + schema in `rectify/core/benchmark/`.
- badread is also installed (cross-check / error-model realism), not the primary.

## Open decisions (carry to user / fold crafter `benchmark_coupling`)
- ~~NanoSim vs badread for Tier-2~~ → **RESOLVED: pbsim3** (per-read MAF; see above).
- Whether Tier-1's error model is badread's ONT model or RECTIFY's own empirical table (held-out either way).
  (Current Tier-1 uses a parametric del-dominant K_DIST; badread cross-check is an option.)
- Human locus-panel composition (SMN1/SMN2 + how many NIC/NNC-rich loci; mirror A549 chr5).
