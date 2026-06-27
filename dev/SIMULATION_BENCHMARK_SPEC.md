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
- **Shared build dep:** both LRGASP-sim (GENCODE) and SIRV truth are **exon-feature
  GTF** → need the exon-GTF loader variant of `gff_panel` (same one SIRV needs). The
  `read_to_isoform` join + GTF → per-read junctions/CPA is the sim-truth path.
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
