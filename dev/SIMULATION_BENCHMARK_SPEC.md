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

## VARIANT/C6 stratum — BUILT + VERIFIED DISCRIMINATING (2026-06-26)
The standing-variant (C6) Tier-1 generator is now in
`scripts/benchmark/sim/controlled.py` (`gen_variant_stratum`), wired into
`generate_corpus` and gated by smoke assertion **(E)**. It measures the
variant-induced discovery FDR the flat **haploid** reference panel suffers — the
first-class C6 discovery-FDR guard the design names (Addendum (b) / §8).

**Discriminating construct, verified vs minimap2 -ax splice -uf (empirically, not
assumed):** a **GT..AG-flanked DELETION variant ≥ ~40bp** is re-expressed by
minimap2 as a spurious **intron (N-op)** instead of a deletion (D) → a FABRICATED,
**variant-adjacent** junction FP (scorer `fp_variant_adjacent`) AND the deletion is
not scored position-exact. Empirical thresholds (5 reads/len): 20/30bp stay a
correct `D`; 40/60/100bp all become `N`. Smoke (reps=20): driver
`fp_variant_adjacent=60`, junction FDR=1.0 on the driver sub-case.

**Specificity controls (kept BECAUSE minimap2 is robust to them):** a SNP 3bp from
the donor, and a SNP ~100bp from any junction, are both correctly called as a
mismatch (`…1X…`) with the real junction preserved (0 FP); a 25bp splice-mimic
deletion stays a correct `D`. The smoke asserts controls produce **0**
variant-adjacent FP — so the stratum proves the FDR is SPECIFIC to the
splice-mimic-deletion context, NOT "any variant near a junction." This guards the
future C6 member against degenerating into a blunt abstain-near-every-variant rule.
Zygosity (HET/HOM) + non-Mendelian VAF (0.33, aneuploid-A549-style) are recorded
per `VariantTruth` for downstream stratification. **Whether C6 REDUCES this FDR is
the next-cycle ablation; the gate only has to MEASURE it, which it now does.**

**Anti-overcount (the artifact-class this benchmark exists to prevent — advisor-
caught 2026-06-26):** each VARIANT read gets its OWN freshly-randomized contig
carrying the same variant KIND (HP varies `k`, STR varies `drop`; VARIANT must NOT
replicate one construct `reps` times or a reported FDR has effective n = #sub-cases,
false confidence at scale). So at `reps=R` the driver carries `3×R` INDEPENDENT
constructs (not 3). `true_locus` is the GROUP label for stratification; the per-read
`chrom` is unique; the region-disjoint split is decided per-contig (each contig = one
read = one partition, never split train↔test). Empirical check (reps=20, 60 distinct
driver constructs): all 60 fabricate; 60 distinct control constructs: 0
variant-adjacent FP. (One control read showed a generic NON-variant-adjacent junction
misplacement — counted in junction FP but NOT in `fp_variant_adjacent`, which
confirms the adjacency tag isolates variant-INDUCED FPs from generic placement
misses.)

**Version-pin caveat:** the discrimination (40bp→N) AND specificity (25bp→D, SNP→X)
both rest on minimap2's intron-vs-deletion length threshold, pinned to **minimap2
2.28 + `-k 14` -ax splice -uf**. 25bp was not probed directly (monotonicity from
30bp→D makes it safe); a minimap2/version or `-k`/seed change could shift the
threshold and silently flip the control. Re-probe the threshold if the aligner
version or splice flags change.

**Triage of the remaining SPEC strata (advisor-reviewed 2026-06-26) — why VARIANT
went first and the others are deferred/scaffold:** a stratum only counts when the
incumbent is shown BELOW ceiling via a metric the scorer *already emits*.
- **VARIANT/C6 — DONE** (above): fully wired (`fp_variant_adjacent`), incumbent
  below ceiling, specificity proven.
- **PARALOG/C4 — DONE (2026-06-26).** Added the missing scorer readout
  (`AlignerScore.locus_accuracy` / `locus_correct` / `locus_incorrect` /
  `locus_mapq0`: mapped contig vs truth-origin `chrom`) FIRST, then
  `gen_paralog_stratum` + smoke **(F)**.
  - **Construct (advisor-corrected — the first cut was the vertical-slice trap):**
    a window-EXCLUDING fragment carries ZERO distinguishing bases → it is
    information-theoretically unrecoverable by ANY method INCLUDING C4 pooling, so
    its below-ceiling 0.5 is structural noise, NOT a closeable gap. The real
    C4-addressable slice is a WEAK-evidence fragment covering EXACTLY ONE
    distinguishing SNP (SMN-like; SNPs spread, frag centered on one): per-read the
    lone SNP is sometimes noise-flipped → minimap2 below ceiling, but POOLING
    denoises.
  - **Verified vs minimap2 (reps=20):** window-SPANNING reads (cover all spread
    SNPs) assign at ceiling locus_acc=1.000 (proves the metric isn't trivially
    failing); WEAK 1-SNP fragments drop to locus_acc≈0.94 (below ceiling, some
    confidently wrong); and the **pooling-majority base at the lone SNP recovers the
    true copy in 6/6 pools** — proving the gap is closeable by C4 pooling NOW, from
    truth, without building C4 (the (D)-equivalent). Whether C4 ACTUALLY closes it
    is the next-cycle ablation.
  - Each contig pair + all its reads sit in ONE split partition (no paralog leaks
    train↔test). **Effective construct diversity = `n_families` (default 3); `reps`
    is per-locus DEPTH (a feature for C4 pooling, not padding) — scale `n_families`
    (not `reps`) for diversity in the Sherlock run.**
- **COVERAGE×Q — cannot be validated this cycle.** Nothing consumes phred (the
  corpus hardcodes `'I'*len(seq)`; the scorer ignores Q). The consumer is C3 (next
  cycle). A Q stratum now is pure scaffold that cannot show discrimination — the
  vertical-slice trap again. Defer or label unvalidated.
- **PANEL_FAILURE/C5 — real validity is Tier-2.** The SPEC ties tail-SIZE to the
  real 5-aligner panel on Sherlock; a Tier-1 "defeat minimap2" read is artificial
  (one aligner ≠ the panel). Build as truth scaffold only; do not let a smoke
  assertion read as "the tail is sized."
- **GENOMIC_A_CPA/C2 — partially wired** (scorer emits `|est−true CPA|`) but its
  comparator (the walkback heuristic) is not wired as an arm, so same partial-
  wiring caveat as Q; a generator is worthwhile but it can't fully discriminate
  until the walkback arm exists.

## TIER-2 BRANCH-A RESULT (2026-06-26) — yeast saturation control, RUN on Sherlock
First scaled Tier-2 run (job 31628436, COMPLETED 2m17s): ~500 spliced yeast
transcripts × 20 copies = 10k reads, both error models, minimap2 -ax splice.
Scope: ANNOTATED-recall + spurious-FDR of ONE aligner (the labels travel in each
`tier2_results/*_summary.json` `_scope` block).

| model | reads placed | ANNOTATED recall | spurious-FDR | TP/FP/FN |
| --- | --- | --- | --- | --- |
| DRS (ERRHMM-ONT)     | 9810/10000 (98.1%) | 0.816 | 0.051 | 8418/451/1902 |
| cDNA (ERRHMM-ONT-HQ) | 9976/10000 (99.8%) | 0.843 | 0.077 | 8698/724/1622 |

**ERROR-FREE saturation baseline (the actual saturation control — run locally on
M1, 250 spliced transcripts, no pbsim, perfect reads): recall=0.985, FDR=0.008**
(TP=255 FP=2 FN=4; the 4 FN are edge soft-clips/1 in-span, NOT coordinate errors).
This is the load-bearing validation: it ISOLATES harness correctness from
minimap2-under-noise. Two things it proves:
1. **The harness/loader is correct** — perfect reads of annotation-derived
   transcripts round-trip to ~ceiling recall with near-zero FDR, so the GFF
   1-based→0-based conversion, projection, and ANNOTATED classification are sound
   (a coordinate off-by-one would crater error-free recall with FN+FP pairs; it did
   not). Combined with the projection fix (0 truth junctions outside a read's
   covered+anchored span), the saturation control PASSES.
2. **The drop to ~0.82–0.84 under simulated error is NOISE-DRIVEN, not artifact and
   not an intrinsic short-intron effect** — if minimap2 intrinsically called these
   short yeast introns as deletions, the ERROR-FREE reads would show it too; they
   do not (0.985). So DRS/cDNA error degrades junction placement, which is the real,
   reportable Tier-2 finding. (Earlier wording claimed a specific "short-intron→
   deletion" mechanism — REFUTED by the error-free baseline; corrected here.)

cDNA (lower error) places more reads + higher recall but ALSO higher spurious-FDR.
A projection bug (every transcript intron assigned to every fragment read) was
found+fixed during real-coord verification BEFORE this run — see the pbsim3_wrapper
commit. The error-free check is `dev/bench/`-style and rerunnable on the M1.

Infra note: the run landed on `sh03-07n10` (an AVX-512-trap node) DESPITE
`--exclude`; the rectify conda env did NOT SIGILL there (the AVX-512 trap is
build-specific — this env is safe), so the exclude not being honored was harmless
here. For Branch B, prefer `--partition=owners` or re-verify the exclude.

**NOT measured here (BRANCH B, gated):** the panel-failure TAIL (needs the wired
multi-aligner panel so `_panel_unplaced_fraction` = placed-by-NO-aligner, + an
injected hard sub-population — a clean run reports tail≈0, a false negative) and
novel-junction NIC/NNC recall (needs isoform injection: exon-skip→NIC, novel-site→
NNC). Aligner inventory for Branch B (rectify env): minimap2+uLTRA+deSALT+gapmm2
present, gmap in other envs, mapPacBio absent.

## ERROR-REALISM — in-silico injector VERDICT (measured 2026-06-27)
The empirical table gives MARGINAL error rates only, and its `--isolation-flank 10` recipe
EXCLUDES clustered/bursty errors by construction (counts only errors with ≥10 exact matches
both sides — confirmed in `empirical_cigar_error_profiler.py`). Measured the error CORRELATION
structure of pbsim3 vs a real haploid-yeast (BY4742) DRS BAM (29,796 reads;
`/scratch/users/kevinroy/rectify_wt_by4742_rep1_*/...minimap2.namesorted.bam`), pbsim
Bernoulli-thinned to the real 1.9% marginal rate (raw ERRHMM-ONT is ~7× too hot — old-R9 vs
modern dorado — so only SHAPE is comparable). **Verdict: pbsim3 is MORE iid/uniform than real
on every axis → an injector at the marginal rate alone is unrealistic where it matters.** Three
data-grounded layers the in-silico injector needs (direction unambiguous):
1. **Per-read OVER-DISPERSION** (hot-read mixture) — real length-adj dispersion index 9.95 vs
   rate-matched pbsim 3.55 (~2.8× under); real per-read-rate tail heavier (p90/median 1.98 vs 1.59).
2. **Within-read BURST/clustering** — real 2.83× excess sub-5bp gaps over its own geometric null
   vs pbsim ~1.2× (cleanest gap; not a rate-ceiling artifact).
3. **Longer multi-base INDEL runs** — pbsim-DRS 19% indels ≥2bp vs real 39% (pbsim-cDNA/ONT-HQ
   happens to match at 39%).
**Magnitude must be calibrated against ABSOLUTE truth, not read-vs-reference:** real-yeast 2.83×
is an UPPER BOUND — read-vs-reference conflates true error with RNA-modification basecall errors
(genuine DRS, motif-clustered) AND minimap2 pile-up near indels/junctions (an ALIGNMENT artifact
a PRE-alignment injector must NOT reproduce). Clean target = SIRV/LRGASP absolute truth
(read-vs-known-sequence removes the alignment-artifact confound). So: measurement = direction;
SIRV = magnitude; the LRGASP three-way also tests whether NanoSim shares the same blind spot.
**Design (deferred to the injector-build cycle):** add a per-read rate multiplier (over-dispersed
mixture) + a self-exciting/Markov burst process + a longer indel-run length model, on top of the
marginal table; calibrate the three parameters against SIRV absolute-truth error structure; keep
HP-only marginal as the backward-compat default. Caveats: shape-only (chemistry mismatch); thinning
makes matched-pbsim clustering a conservative (low) estimate; full analysis `/tmp/err_corr.py` +
`/scratch/users/kevinroy/err_corr_work/` on Sherlock.

### INJECTOR BUILD (2026-06-27, session-5) — M1-local code DONE; fitting/validation SHERLOCK-gated
The 3-layer error-realism injector is BUILT, composition-verified on the M1, and the Sherlock
validation is turnkey. **Files** (all benchmark-only): `scripts/benchmark/sim/error_injector.py`
(model + measurement + calibration + `self_check`), `tests/test_error_injector.py` (8 fast tests),
`scripts/benchmark/read_reliability_probe.py` (PI-#2 mechanism), `scripts/benchmark/error_realism_validate.py`
(Sherlock `measure-bam` / `inject-fastq`).

- **Model** — a generative process applied to ANY base read (clean controlled read OR a pbsim read;
  source-decoupled per the validation protocol): (1) per-read multiplier (Gamma E[m]=1 / 2-comp mixture)
  → over-dispersion; (2) 2-state cold/hot burst HMM modulating the LOCAL hazard (marginal-preserving:
  cold factor derived so the stationary-weighted mean factor = 1) → within-read clustering; (3) a
  fat-tailed geometric indel-run-length pmf. `InjectorParams.null()` reproduces the marginal/Poisson
  single-base baseline (backward-compat).
- **TRUTH RULE (the SESSION-2 lesson, advisor-reinforced):** injected background errors are a SEPARATE
  per-read ERROR TRACK (`List[ErrorEvent]`), **NEVER `IndelTruth`** — writing them as IndelTruth
  re-triggers the scorer's `has_unexplained` zeroing that broke (C)/(D) for pbsim. The track IS PI-#2's
  per-read hotness label.
- **Calibration** — coordinate descent against the COMBINED simulated output (NOT three independent MoM
  solves — the layers cross-contaminate: the burst HMM inflates the dispersion index, so layer-1's Gamma
  shape converges to ~1.26 not ~tiny, the loop accounting for the burst's dispersion contribution).
  Converges robustly across seeds to the PLACEHOLDER targets: **disp 9.45–10.24 / gap5x 2.71–2.95 /
  run≥2 0.40** (targets 9.95 / 2.83 / 0.39). Gap statistic defined on error EVENTS (a k-bp run = one
  event), so layer-3 runs do NOT contaminate layer-2's sub-5bp-gap count.
- **The M1 `self_check` is a COMPOSITION/INTERACTION check, NOT realism validation** (avoiding the
  SESSION-2 overclaim): NULL→Poisson baseline; each layer moves its own statistic; layers compose
  without one zeroing another; calibration reaches the (just-identified) moments. **Realism is
  established ONLY by SIRV magnitude calibration + a DISTRIBUTIONAL comparison — both Sherlock-gated.**
- **Magnitude is UNKNOWN, not "to refine"** — 2.83× is an UPPER BOUND (read-vs-ref conflates RNA-mod +
  minimap2 pile-up). M1 params are PLACEHOLDER; the SIRV absolute-truth re-fit sets the real magnitude.
- **PI-#2 MECHANISM probe (M1-safe qualitative half; advisor-gated):** exonic error density vs
  junction-window error across reads — null r≈−0.02, **over-dispersion-only r=0.955** (global hotness →
  exonic predicts junction), **burst-only r=0.011** (local; distant regions independent), combined
  r=0.841. So the reliability covariate has signal BUT is imperfect → **SOFT down-weight / posterior
  input, NOT a hard filter** (confirms PI refinement (a)). The FDR-LIFT NUMBER is magnitude-sensitive /
  SIRV-gated and was deliberately NOT computed.
- **SHERLOCK-gated next steps (auth down this session — do NOT re-open the master):** (a) RE-PROFILE the
  real BAM with `empirical_cigar_error_profiler --isolation-flank 0` (gate OFF, captures bursts; default
  is already 0); (b) `error_realism_validate.py measure-bam` on {real, pbsim, pbsim+injector} with
  `--thin 0.019` → confirm +injector closes the pbsim→real gap; (c) re-fit `calibrate_params` to SIRV
  absolute-truth targets (the real magnitude); (d) THEN the PI-#2 FDR-lift on stratum (G). Protocol in
  `error_realism_validate.py` docstring.
- **GUARDS to bake into the Sherlock validation (advisor, session-5) — the self-checks are internally
  circular (calibrate to `REAL_TARGETS`, confirm you hit them), so these are the only things that make
  the validation MEANINGFUL:**
  1. **RE-DERIVE the targets with THIS module's `measure_error_structure`** — the hardcoded 9.95/3.55/2.83
     came from the LOST `err_corr.py` (possibly a different statistic definition); first Sherlock action is
     to re-measure real AND pbsim, and `calibrate_params` to the FRESHLY-measured real target, not the stale
     constant.
  2. **LENGTH-CONTROL the comparison** — the dispersion index is ~linear in read length (≈ 1 + rate·L·v), so
     compare real vs pbsim+injector binned-by-length / at a common window, and calibrate to the real reads'
     actual length distribution (not the flat 600bp M1 default), else over-dispersion is confounded with
     length differences between the two sources.
  3. **Verify the two measurement paths AGREE** before trusting the cross-source comparison: inject into
     clean reads → align → `measure-bam` should reproduce the injected-truth `measure_error_structure`. (Known
     discrepancy to check: `events_from_alignment` can place an insertion and a substitution at the SAME ref
     pos → gap=0, which `inject` never emits (min gap 1) → can inflate the BAM's sub-5 fraction.)
  4. **Derive the error TYPE split from the empirical table** — `frac_sub/ins/del` (0.55/0.15/0.30) and
     `base_rate` are PLACEHOLDER; DRS is deletion-dominant, so re-fit them from the real profiler output.

### SHERLOCK VALIDATION RESULT (2026-06-27, session-6) — RAN; verdict empirically nailed
Sherlock auth restored; ran the validation with the NEW `measure_error_structure` (self-consistent code
across real/pbsim/injector). **All measurements on real BY4742 DRS (`...rep1_26167419_0/...minimap2.namesorted.bam`)
and the `err_corr_work` pbsim reads.**
- **Guard #1 (re-derive targets) — PASSED.** New code REPRODUCES the lost err_corr.py verdict: dispersion_index
  9.64 (prior 9.95), indel_run≥2 0.39 (prior 0.39), p90/median 2.03 (prior 1.98). Only `sub5_gap_excess` differs
  — **5.28 vs prior 2.83** — a definitional difference (this stat = gaps between error EVENTS in ref coords). Real
  marginal event rate 0.0153; length-invariant `overdisp_v` 0.70 (stable across the [600,1000] window).
- **pbsim is DECISIVELY shape-deficient (the SPEC verdict, empirically confirmed with self-consistent code).**
  pbsim DRS aligned to its templates, THINNED to the real 0.0153 rate: `overdisp_v` **0.054 vs real 0.70**,
  `sub5_gap_excess` **1.16 vs 5.28** (≈no clustering — matches the SPEC's pbsim ~1.2×), `indel_run≥2` **0.25 vs
  0.39**. (The "~13×" overdisp ratio is DIRECTIONAL, not precise: Bernoulli-thinning a clustered process breaks
  up runs, so thinned-pbsim UNDERSTATES pbsim's own clustering/over-dispersion — the deficiency verdict is robust
  but the multiplier is soft. The relative comparison uses the SAME code both sides, so it is robust regardless of
  the absolute-definition question on 5.28.) (Raw pbsim rate 0.119 = ~7× too hot, R9-era — only SHAPE is comparable; the literal
  "pbsim+injector" STACKING is rate-incoherent and dropped — pbsim's role was the shape-deficient baseline,
  injector-on-clean-reads is the real vehicle.)
- **Real within-read AUTOCORRELATION (PI-#2(i) on REAL data) = 0.34 (frac 0.3) / 0.28 (frac 0.2)** — a MODERATE
  positive head-vs-tail error-density correlation. So a GLOBAL hotness component exists (the PI-#2 reliability
  covariate has REAL signal) but is PARTIAL (r~0.3, not ~0.95) → **SOFT down-weight, not hard filter** — PI
  refinement (a) confirmed on real reads (stronger evidence than the synthetic r=0.955 over-dispersion-only probe).
  This real r is the THIRD constraint that identifies the global-vs-local split (overdisp_v + gap5x cannot).
  **ATTRIBUTION CAVEAT (control for it in the FDR-lift test):** head-vs-tail autocorrelation conflates per-MOLECULE
  hotness (the over-dispersion we want) with a per-TRANSCRIPT alignability effect — both read ends share transcript
  identity, and some transcripts are intrinsically harder to align (low-complexity / paralogy / GC), elevating both
  ends independent of any per-molecule over-dispersion. r≈0.30 is correct as a measurement; the attribution to
  per-molecule hotness is one inferential step beyond it. Harmless here ("partial signal → soft not hard" holds
  regardless of WHY a read is error-prone), but for the gated PI-#2 FDR-lift it is load-bearing: if r is partly
  per-transcript, down-weighting hot reads could suppress novel-junction support from hard-to-align transcripts (a
  DISCOVERY BIAS, not just lost power). So stratify by transcript / compare within- vs across-transcript autocorr to
  isolate per-molecule hotness before claiming the lift.
- **MODEL-EXPRESSIVENESS FINDING (key input to the SIRV-fit cycle):** the 2-state burst HMM CANNOT jointly match
  real's (overdisp_v 0.70, gap5x ~5, autocorr 0.30) — longer bursts (needed for gap5x) raise per-read hot-fraction
  variance → raise autocorr, so strong clustering and moderate global correlation are COUPLED in the model, whereas
  REAL DECOUPLES them (real over-dispersion is mostly LOCAL). Matching real likely needs a 3-state / self-exciting
  (Hawkes) burst. Locked a HAND-PICKED `placeholder_params` (gamma_shape 5, hot_factor 8, p_hot_to_cold 0.2,
  p_cold_to_hot 0.025) matching the autocorr SPLIT (0.30) so it does NOT overstate the PI-#2 covariate; achieved
  overdisp_v ~0.24 < real 0.70 is the documented model gap. PLACEHOLDER-PENDING-SIRV; auto-fit NOT used (advisor).
- **Guard #3 (alignment-only inflation) — GAP-FREE case only: inflation ~0; junction case UNRESOLVED.**
  Inject clean errors (truth gap5x 4.44) → re-align → aligned gap5x **3.99** (slightly LOWER, not higher). This
  establishes the WEAKER claim: minimap2 does NOT fabricate clustering out of scattered clean errors in GAP-FREE
  alignment. It does NOT speak to real's 5.28, because the test differs on BOTH axes that matter: it used
  `-ax map-ont` to the gap-free SAME templates (no introns), whereas real's 5.28 is `-ax splice` read-vs-GENOME
  crossing introns — and the SPEC's artifact concern is specifically minimap2 pile-up NEAR junctions/true indels,
  the one condition this test cannot exercise. So: **alignment inflation of gap5x is ~0 in the gap-free case; the
  junction-spanning case is unresolved → SIRV (junction-spanning, splice-aligned) DECIDES whether the true target
  is near 3 or near 5. Do NOT pre-bias the SIRV fit.**
- **STILL SIRV-gated (the genuine open work):** absolute magnitude (read-vs-ref is an upper bound, though guard #3
  shrinks the alignment-artifact worry); the error TYPE split (DRS deletion-dominant); and the 3-state/self-exciting
  burst upgrade to decouple clustering from autocorr. Then the PI-#2 FDR-LIFT on stratum (G).

### SIRV ABSOLUTE-TRUTH RESULT (2026-06-28, session-7) — RAN; clean magnitude leans LOW; confounds noted
First absolute-truth (read-vs-known-sequence, mod-free, junction-spanning) error-structure measurement.
Job 31823846 COMPLETED rc=0 (2m13s). Real WTC11 DRS (ENCODE `ENCFF155CFF`) aligned to SIRV-Set 4
(`minimap2 -ax splice -uf -k14 --MD`); 4424 primary-mapped spike-in reads; measured with the SAME
`measure_error_structure` as yeast/pbsim/injector. Outputs `/scratch/users/kevinroy/sirv_work/`.

**Numbers (rate+length-controlled — the valid comparison):** SIRV thinned to the yeast rate 0.0153 in
the [600,1000] window: **overdisp_v 0.17, sub5_gap_excess 2.27, indel_run≥2 0.36, p90/med 1.87, autocorr
r 0.20** (508 reads). vs real-yeast read-vs-genome **0.70 / 5.28 / 0.39 / 2.03 / 0.30**.

**READING (HUMBLE — the magnitude is NOT yet trustworthy; a measurement bug + confounds block it):**
The defensible facts: (a) the pipeline works end-to-end and absolute-truth SIRV reads are obtained;
(b) the [600,1000] window is **junction-spanning-dominated** (369 reads with an N-op / SIRV1–7 vs 139
monoexonic ERCCs — so NOT an ERCC artifact); (c) **`indel_run≥2` = 0.36 vs yeast 0.39** is the single most
transferable number (rate/architecture-robust); (d) **autocorr r ≈ 0.20 (positive)** corroborates
"partial global hotness → SOFT down-weight, not hard filter" on absolute-truth data (this stat uses
head/tail density correlation and survives the confounds below). The clustering/over-dispersion MAGNITUDE
(gap5x 2.27, overdisp_v 0.17) is **NOT yet usable** — do not conclude "leans LOW (~3 not ~5)" from it.

**WHY the magnitude is blocked — a measurement bug (advisor-caught) + confounds:**
1. **N-in-span bug in `measure_error_structure` (the blocker).** `events_from_alignment` does
   `ref_span += length` AND `rpos += length` for the N (intron) op. So for SPLICED reads the span is
   INFLATED by total intron length → rate DEFLATED (this is exactly the full-set 0.0076 vs windowed 0.059
   anomaly), AND `rpos` jumps the intron so error GAPS cross introns → the sub5_gap_excess / dispersion
   metrics are corrupted for any junction-spanning read. **FIX (next): measure error rate PER EXONIC BASE
   (exclude N from span) and compute gaps WITHIN exons (don't let a gap straddle an intron)** — a per-exon
   mode in `measure_error_structure`. The yeast 5.28 has the SAME bug but yeast introns are short → smaller
   effect, so even the yeast↔SIRV comparison is muddied until both are re-measured per-exonic-base.
2. **Chemistry — RNA002-era.** ENCODE WTC11 DRS native in-window rate **0.059 (~4× yeast)** = OLDER
   chemistry, not the RNA004 target nor the modern-dorado yeast. Shape-comparable at best (pbsim status).
   **Confirm on RNA004 SIRV (LongBench) before locking.**
3. **Thinning distortion.** Comparing SIRV (0.059) to yeast (0.0153) needs thinning, which BREAKS UP runs
   → understates clustering. No clean apples-to-apples across two chemistries.
**NEXT (in order):** (i) add the per-exonic-base mode to `measure_error_structure` + re-measure SIRV
spliced reads AND re-measure real yeast the same way (makes them comparable); (ii) confirm on RNA004 SIRV
(LongBench); (iii) THEN re-fit or confirm the placeholder. Auto-refit NOT done (advisor discipline). Do
NOT drop the 3-state/Hawkes upgrade on this evidence — the magnitude question is still open.

**CORRECTED (2026-06-28, session-7b) — the N-in-span fix landed; numbers re-measured EXONIC both sides.**
The fix (exclude intron N from span + exon-local positions, `error_injector.events_from_alignment(...,
exonic_coords=True)`) materially changed things — and confirmed the bug was inflating yeast clustering:
| metric (EXONIC) | SIRV abs-truth (thin 0.0153, spliced SIRV1-7) | yeast real (native 0.0187) |
| --- | --- | --- |
| overdisp_v | 0.10 | 0.54 |
| sub5_gap_excess | 1.71 | **4.28** (was 5.28 with the bug) |
| indel_run≥2 | 0.37 | 0.39 |
| autocorr r | 0.50 | 0.47 |
Reads-in-window jumped (SIRV 508→1109, yeast →114k) because exonic spans are shorter. **READING (still
needs the adversarial panel + RNA004):** (1) the bug inflated yeast gap5x 5.28→4.28 and overdisp 0.70→0.54
— real clustering is LOWER than previously stated. (2) `indel_run≥2` (0.37 vs 0.39) remains the rock-solid
transferable number. (3) **NEW AMBIGUITY:** SIRV shows LOW clustering/over-dispersion (gap5x 1.7,
overdisp_v 0.10) but HIGH autocorr (0.50) — internally suspicious (strong global correlation with weak
over-dispersion). Most likely the **per-TRANSCRIPT alignability confound** (SPEC attribution caveat),
sharpened because SIRV has only ~7 multi-exon transcripts → reads from a hard transcript are error-rich at
both ends → autocorr inflated without per-molecule hotness. This is load-bearing for PI-#2 (per-molecule
vs per-transcript) → stratify autocorr by transcript before trusting r. Still RNA002-confounded; confirm
on RNA004 (LongBench) + SG-NEx. Adversarial panel pending on (3).

### ORGANISM-SPECIFIC ERROR MODEL — design note (2026-06-28, advisor-checked; decisions DEFERRED, not pre-answered)
Question raised (PI): yeast vs human differ in (a) RNA modifications, (b) pA-tail length, (c) exon length /
multi-intron fraction — do we need DISTINCT yeast vs human error models? Decomposition + the measurement-design
confounds, recorded as OPEN questions (we do NOT pre-answer "one model vs two"). The three named differences
split across THREE different layers, not all of them the error model:

- **(1) Basecaller-INTRINSIC error (HP / k-mer-dependent miscalls) — the per-context error FUNCTION is
  organism-agnostic** (the pore+dorado don't know the species; the same k-mers miscall the same way). BUT the
  AGGREGATE statistics it produces are NOT organism-agnostic — they depend on sequence composition (human HP-length
  distribution, GC, low-complexity content differ from yeast). **PRECONDITION on any yeast→human transfer claim:**
  the injector must be driven by a CONTEXT-CONDITIONED error map (the empirical k-mer/HP table — currently
  PLACEHOLDER `base_rate`+`frac_*`) and APPLIED to human clean reads (run the same context function over human
  sequence); you must NOT port yeast AGGREGATE params. With a flat base_rate the transfer claim is false.
- **(2) Modification-driven error — the ONE genuinely organism-specific error mechanism.** Vegetative
  *S. cerevisiae* mRNA is m6A-POOR (Ime4/m6A is meiosis-induced); human mRNA carries pervasive DRACH m6A + Ψ.
  DRS basecallers are trained on canonical bases → a modified base gives a SYSTEMATIC, motif-localized miscall
  (the same signal m6A/Ψ detectors exploit). Mechanism if/when needed: a MOTIF-CONDITIONED hazard multiplier
  (elevated sub/del at DRACH / known-Ψ motifs) — deterministic-on-SEQUENCE (fixed positions), NOT a random burst;
  cheap to layer on. **Do NOT build speculatively — characterize first.**
- **(3) Transcript ARCHITECTURE (shorter exons, more multi-intron genes) — NOT an error-model property at all.**
  It changes TASK difficulty (junction-discovery FDR, NIC/NNC opportunity), handled by the transcript PANEL +
  strata (yeast=saturation control; human SMN1/SMN2 + A549-chr5 NIC/NNC-rich loci), not the injector.
- **pA tail is DOUBLE-COUNTED on purpose:** architecture (CPA truth / C2, a panel property) AND an out-of-
  calibration ERROR regime — we calibrated on ~50–70 nt yeast pA; 150–250 nt human pA is a known DRS
  stall/length-undercall regime, and HP error is already length-dependent (possibly super-linear), so the HP
  mechanism is EXTRAPOLATED OUT OF RANGE for human pA, not "more of the same." Flag both.

**MEASUREMENT-DESIGN CONFOUNDS (advisor — do NOT design the SIRV cycle around the naive subtractions):**
- **SIRV is the basecaller floor for SIRV's OWN composition — NOT a mod-free proxy for human sequence.** SIRV/ERCC
  are IVT constructs with their own designed k-mer/GC. So **`real-human − SIRV` does NOT isolate the mod term** —
  it conflates modification with the SIRV-vs-human composition difference (same guard-#3 / per-transcript shape:
  the control differs on a second axis). The CLEAN mod isolation is a **WITHIN-MOTIF modified-vs-unmodified
  contrast**: error at high-stoichiometry m6A sites vs low/zero-stoichiometry instances of the SAME DRACH motif,
  using an ORTHOGONAL map (miCLIP / GLORI). (DRACH-vs-genome-average is still composition-confounded — not all
  DRACH are methylated.)
- **Yeast-real is NOT the mod-free anchor — only SIRV is.** The real-yeast 5.28 already contains yeast's own mod +
  tRNA/rRNA-contamination error. Call yeast "LOWER-mod than human," never "mod-free reference."
- **Per-molecule vs per-transcript (ties to the PI-#2 attribution caveat above):** mod-driven error is
  per-TRANSCRIPT/per-POSITION (sequence-deterministic), NOT per-molecule hotness → a high-mod human benchmark
  inflates the per-TRANSCRIPT autocorrelation component, making the per-molecule reliability covariate HARDER to
  isolate. Same control (stratify by transcript) handles both.

**RECORDED PLAN (prove-don't-assert; nothing pre-decided):** one injector ENGINE, organism-specific PARAMETER
SETS (not two models). SIRV = basecaller floor (for SIRV composition); the pbsim3-vs-NanoSim-vs-real-SIRV three-way
answers the SIMULATOR-realism question (the original SPEC purpose); the yeast-vs-human MOD question needs the
within-motif contrast (orthogonal map). **"One model vs two" stays UNDECIDED pending those measurements** — do not
build a human error model, and do not assume one suffices, until the within-motif mod magnitude + the
context-conditioned-transfer check are in.

## EXTERNAL-VALIDITY DATA PLAN — real ONT with ground truth (2026-06-26, researched)
Simulation validates placement MECHANICS but cannot validate that pbsim's error
model matches real ONT (the §"Two validity claims" gap). The transfer check needs
REAL ONT reads with independent truth. Our scorer is read-source-agnostic (scores
any BAM vs a truth table), so spike-in reads drop in via the vendor GTF + the
`gff_panel` loader. The decisive distinction:
- **ABSOLUTE truth** = synthetic spike-ins (Lexogen SIRV-Set, Garvan Sequins, ERCC):
  exact junctions/3′/indels by construction, IMMUNE to ONT systematic error — the
  ONLY truth that can validate the homopolymer/CPA-sensitive parts (where C1/C2 live).
- **BIASED-but-real truth** = per-molecule consensus (UMI / R2C2): removes RANDOM but
  NOT SYSTEMATIC error — the methods' own authors confirm consensus still undercalls
  long homopolymers (>10nt) and KIV-2 HP errors survive in >50% of UMI consensuses.
  So NEVER use consensus truth to validate HP-indel correction; it is wrong in the
  SAME direction as the reads. Use only for non-HP junction/indel tolerance.
- Only **native DRS** gives real poly-A 3′ ends; all cDNA/R2C2/UMI 3′ ends are
  oligo-dT/template-switch confounded (not native-CPA truth). DRS-UMI: none public.

**Ranked datasets (bring in this order):**
1. **LongBench** (WEHI/Ritchie 2025) — START HERE. ONT DRS on **RNA004** (current
   chemistry) + cDNA + matched Illumina, 8 cell lines, **SIRV-Set 4 + Sequins**
   spike-ins = absolute truth on the chemistry the simulator must match. AWS Open
   Data `s3://longbench-data/` (ap-southeast-2, --no-sign-request); bioRxiv
   2025.09.11.675724; GitHub mritchielab/LongBench.io. NOTE: spike-in reference
   FASTA/GTF NOT in bucket (get Lexogen SIRVsuite + Garvan sequinstandards.com);
   S3-only (no GEO/ENA); mind Dorado RNA v5.2.0 3′-coverage deficit (fixed v5.3.0).
2. **SG-NEx** — broadest spike-in corpus, RNA002 DRS + cDNA, ships a ready-made
   combined GRCh38+Sequin+SIRV+ERCC GTF (easiest to encode at scale). `s3://sg-nex-data/`;
   ENA PRJEB44348. Caveat: DRS is RNA002 (tests mechanics on real reads, not RNA004
   transfer); the one RNA004 sample has no spike-in.
3. **LRGASP** — uniquely ships BOTH a simulated ONT track (read-origin truth) AND real
   SIRV reads on the SAME samples → turnkey simulation-vs-real comparison (their sim
   is NanoSim, informative vs our pbsim3). ENCODE internal_tags=LRGASP; SIRV-Set 4 on
   Synapse syn25683367/syn25683630. Caveat: DRS kit version not confirmed RNA002.
4. R2C2 (Vollmers; PRJNA971991 etc.) / UMI cDNA (Sicelore GSE130708): biased-but-real,
   junction/indel only, NEVER for HP or native-CPA. Supplementary.

### LRGASP acquisition plan (URL-verified 2026-06-27) — the turnkey THREE-WAY
LRGASP uniquely enables a same-sample, same-scorer triangulation: **our pbsim3 vs
LRGASP's NanoSim sim vs real SIRV reads** → two simulators bracketing reality (if
BOTH under-cluster vs real SIRV = a shared simulator blind spot to inject around).
**Role split:** LRGASP is **RNA002-era** (date-inferred 2021) → it answers the
error-MODEL / sim-realism question; **LongBench is the RNA004** source for the
current-chemistry transfer question. Use both for their distinct roles.

All open/no-auth (Synapse syn IDs are gated; UCSC CGL mirrors are not). Human WTC11:
- **Sim ONT (one COMBINED file — DRS+cDNA NOT split on the public mirror):**
  `https://cgl.gi.ucsc.edu/data/LRGASP/data/simulation/human_simulation/human.ONT.simulated.fq.gz`
  (~30.3 GiB). Truth join: `.../ground_truth/ONT.simulated.read_to_isoform.tsv.gz`
  (211 MB; **no header, 2 cols `read_id<TAB>transcript_id`**; read_id=`ONT_simulated_read_N`;
  transcript_id = versioned Ensembl, human ENST + **mouse ENSMUST decoys interleaved**
  = their built-in novel-discovery test). Annotation: `.../ground_truth/hs_GENCODE38.basic_annotation.gtf.gz`
  (28 MB, **exon-feature GTF**) + `...counts.tsv.gz` + `...novel_isoforms.tsv.gz`.
  CPA truth = annotated 3' terminus (polyA appended pre-sim).
- **Real ONT + SIRV (WTC11), ENCODE @@download, no auth:** DRS = **ENCSR392BGY**
  (3 reps ~3.8 GiB: ENCFF155CFF/ENCFF600LIU/ENCFF771DIX); cDNA = **ENCSR539ZXJ**
  (3 reps ~46 GiB; 1 rep ~9.5 GiB pilot). URL form
  `https://www.encodeproject.org/files/<ACC>/@@download/<ACC>.fastq.gz`.
- **SIRV-Set 4 reference (ABSOLUTE truth), open CGL mirror:**
  `.../references/lrgasp_sirv4.fasta.gz` (132 KB) + `.../references/lrgasp_sirvs4.gtf.gz`
  (4.6 KB; **176 transcripts, exon features**; multi-exon SIRV1–7 + monoexonic
  long-SIRVs + ~92 ERCCs). Combined genome+SIRV for aligning real reads:
  `.../references/lrgasp_grch38_sirvs.fasta.gz` (886 MB) + `lrgasp_gencode_v38_sirvs.gtf.gz` (50 MB).
- **Staging:** minimal pilot ~42 GiB, full ~81 GiB. Cluster only, never the M1.
- **Shared build dep — BUILT (2026-06-28, session-7):** both LRGASP-sim (GENCODE) and
  SIRV truth are **exon-feature GTF** → `gff_panel.build_panel_from_gtf` (exon rows grouped
  by `transcript_id`, introns derived from adjacent-exon gaps) is the drop-in loader, returning
  the same `(models, pairs, donors, acceptors)` as the GFF path. **Validated on real SIRV-Set 4**
  (176 transcripts, 288 junctions all ANNOTATED, 99.3% canonical at derived boundaries; tests in
  `tests/test_gff_panel_gtf.py`). Transcript key = verbatim `transcript_id` (version included) so
  the `read_to_isoform` join matches. The `read_to_isoform` join + GTF → per-read junctions/CPA is
  the (still-Sherlock-gated) sim-truth path.
- **CORRECTION (verified):** use the **`hs_GENCODE38.basic_annotation.with_mm_M27.gtf.gz`**
  (29.8 MB), NOT the human-only `basic` GTF — **~20% of the "human" sim reads are MOUSE
  transcripts** (`ENSMUST`, the artificial-novel decoys: 79,489 ENST vs 20,511 ENSMUST in
  the first 100k rows). Scoring against the human-only GTF falsely fails ~20% of reads as
  "no truth." (Or filter read_to_isoform to `ENST` only, dropping those reads.)
- **The three-way is DISTRIBUTIONAL, not locus-matched:** NanoSim reads are GENCODE-origin;
  the real ABSOLUTE truth is SIRV-only (endogenous real reads have annotation-level, not
  per-read, truth). So compare the SHAPES of error distributions (junction-placement error,
  indel-by-HP-length, CPA-offset) across pbsim3/NanoSim/real-SIRV — which is exactly what
  "do both simulators share an error-realism blind spot vs real" needs. A strictly
  locus-matched SIRV three-way would require re-running NanoSim on SIRV (not turnkey) —
  optional rigorous follow-up.
- **CPA caveat:** only real DRS + SIRV give native poly-A 3' ends; real cDNA and the NanoSim
  track 3' ends reflect priming → weight CPA conclusions to the DRS-real + SIRV subset.
- Caveats: RNA002 (not RNA004; date-inferred, "RNA002" not in metadata); sim is one
  combined FASTQ (per-chemistry split would need re-running lrgasp-simulation); DRS
  reads carry uracils + nonstandard ONT names (benign ENCODE flags). All URLs curl-verified
  (CGL mirror fully open; Synapse not needed). Pilot ~37.7 GB / full ~87.5 GB.

**Recommended first transfer experiment:** align LongBench RNA004 DRS (+cDNA) spike-in
reads to the SIRV/Sequin reference → score junction recall/FDR with our scorer →
compare to the pbsim3 Tier-2 numbers (DRS 0.816 / cDNA 0.843). A large gap = the
simulator's error model doesn't transfer (re-tune or caveat); a small gap = the
simulation gate's junction conclusions hold on real data. Heavy → run on the cluster,
never relay reads through the M1.

## Open decisions (carry to user / fold crafter `benchmark_coupling`)
- ~~NanoSim vs badread for Tier-2~~ → **RESOLVED: pbsim3** (per-read MAF; see above).
- Whether Tier-1's error model is badread's ONT model or RECTIFY's own empirical table (held-out either way).
  (Current Tier-1 uses a parametric del-dominant K_DIST; badread cross-check is an option.)
- Human locus-panel composition (SMN1/SMN2 + how many NIC/NNC-rich loci; mirror A549 chr5).
