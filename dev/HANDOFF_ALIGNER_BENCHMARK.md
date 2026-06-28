# HANDOFF — Deliverable A: the simulation ground-truth benchmark (the GATE)

**Agent:** dedicated benchmark-builder (isolated worktree, branch
`worktree-agent-a25a2c1e784ad37dc`, based on `drs-validation-rebuild` so the reuse
primitives + design docs are present). **NEVER commit to `drs-validation-rebuild`.**
**Updated:** 2026-06-27 (session 4).

---

## SESSION-6 (2026-06-27) — SHERLOCK VALIDATION RAN; verdict empirically nailed; placeholder LOCKED

Sherlock auth restored. Ran the injector validation with the NEW `measure_error_structure`
(self-consistent code across real/pbsim/injector). Full numbers in SPEC §"SHERLOCK VALIDATION
RESULT". This is a CLEAN CYCLE BOUNDARY (advisor): the deliverable is validated; the remaining
magnitude fit is the SIRV gate working as designed, not incompleteness.

**RESULTS (all gates green: smoke, injector self_check, 8 injector tests):**
- **Guard #1 PASSED** — new code reproduces the lost-err_corr verdict (disp 9.64 vs 9.95, run≥2 0.39
  vs 0.39, p90/med 2.03 vs 1.98). `sub5_gap_excess` = 5.28 (vs prior 2.83 — definitional). Real rate
  0.0153, length-invariant `overdisp_v` 0.70.
- **pbsim DECISIVELY shape-deficient** (thinned to real rate, self-consistent code): overdisp_v 0.054
  vs 0.70 (~13×), gap5x 1.16 vs 5.28, run≥2 0.25 vs 0.39. SPEC verdict empirically confirmed.
- **Real within-read autocorrelation r ≈ 0.30** (PI-#2(i) on REAL data) — moderate global hotness →
  covariate has REAL but PARTIAL signal → SOFT down-weight not hard filter (confirms PI refinement (a)
  on real reads). This r is the third constraint pinning the global/local split.
- **Model-expressiveness finding:** the 2-state HMM CANNOT jointly match real (overdisp_v 0.70, gap5x ~5,
  autocorr 0.30) — clustering↔autocorr are coupled in-model, real decouples (mostly LOCAL over-dispersion).
  Locked a HAND-PICKED `placeholder_params` matching the autocorr split (0.30); achieved overdisp_v ~0.24
  is the documented gap. SIRV fit likely needs a 3-state/self-exciting burst. Auto-fit NOT used (advisor).
- **Guard #3 (alignment-only inflation) — GAP-FREE case ~0; junction case UNRESOLVED.** Inject clean
  errors (truth gap5x 4.44) → realign → 3.99 (LOWER). Establishes only: minimap2 doesn't fabricate
  clustering from scattered clean errors in gap-free alignment. Does NOT speak to real's 5.28 (that is
  `-ax splice` read-vs-genome crossing introns; this test was `-ax map-ont` to gap-free templates — the
  junction pile-up the SPEC worried about is untested). → SIRV (junction-spanning) decides whether the
  true target is near 3 or 5; do NOT pre-bias the SIRV fit toward 5.

**CODE CHANGES this session** (committed): calibration metric → length-invariant `overdisp_v` (not
length-dependent `dispersion_index`); `REAL_TARGETS` re-measured + `head_tail_autocorr` 0.30 added;
hand-picked `placeholder_params` (PLACEHOLDER-PENDING-SIRV); `self_check` (3) decoupled from the
contaminated magnitude (asserts reachable composition bands + autocorr split); `error_realism_validate.py`
gained `autocorr-bam`, length-window filter, truth-track report on `inject-fastq`; fixed the burst
monotonicity bug (lower p_cold_to_hot → more clustering). New stat `overdisp_v`/`head_tail_autocorr`.

**RESUME (session-6 boundary) — SIRV-gated work only; nothing in flight. When resuming:**
1. **SIRV/LRGASP absolute-truth fit** (SPEC §EXTERNAL-VALIDITY has URL-verified accessions): build the
   exon-GTF `gff_panel` loader variant (still unbuilt), align SIRV reads to the SIRV reference,
   `measure-bam` + `autocorr-bam` → the CLEAN magnitude targets (junction-spanning, so it also re-checks
   guard #3's alignment-inflation on junctions). Re-fit; derive the error TYPE split (DRS deletion-dominant).
2. **3-state / self-exciting burst upgrade** — to decouple clustering (gap5x) from autocorr, which the
   2-state HMM cannot (this session's finding). Only if the SIRV targets confirm the real decoupling.
3. **PI-#2 FDR-LIFT on stratum (G)** — inject hot reads (calibrated, post-SIRV) into JUNCTION_DISCOVERY,
   measure novel-junction FDR with vs without hot-read down-weighting. Prove the lift before shipping.
   **CONTROL (advisor):** the real r≈0.30 conflates per-MOLECULE hotness with per-TRANSCRIPT alignability;
   down-weighting hot reads could suppress novel-junction support from hard-to-align transcripts (a DISCOVERY
   BIAS). Stratify by transcript / compare within- vs across-transcript autocorr to isolate per-molecule hotness
   before claiming the lift. (SPEC §SHERLOCK VALIDATION "ATTRIBUTION CAVEAT".)
4. **ORGANISM-SPECIFIC error model (design note recorded, decisions DEFERRED)** — before any human panel:
   yeast/human differ in RNA mods (the one organism-specific ERROR mechanism: yeast m6A-poor, human DRACH-pervasive),
   pA length (out-of-calibration HP regime), exon/intron architecture (panel property, not error model). Plan: ONE
   engine + per-organism PARAM SETS; "one model vs two" stays UNDECIDED. CONFOUNDS to honor: `human−SIRV` does NOT
   isolate the mod term (SIRV is IVT, own composition) — use a WITHIN-MOTIF modified-vs-unmodified contrast
   (miCLIP/GLORI); yeast-real is lower-mod, NOT mod-free; transfer requires a CONTEXT-CONDITIONED error map applied
   to human sequence (not ported yeast aggregates). Full note: SPEC §"ORGANISM-SPECIFIC ERROR MODEL".
5. C1 member code remains the NEXT GATED CYCLE.
Sherlock paths: real BAM `/scratch/users/kevinroy/rectify_wt_by4742_rep1_26167419_0/...minimap2.namesorted.bam`;
pbsim reads + templates `/scratch/users/kevinroy/err_corr_work/`; bench code
`/home/groups/larsms/users/kevinroy/aligner_bench_live/` (rsync the worktree `scripts/benchmark/` first).
Gates: `pysam-python scripts/benchmark/smoke_roundtrip.py --out /tmp/x --reps 20` + `pysam-python
scripts/benchmark/sim/error_injector.py` → both exit 0.

---

## SESSION-5 (2026-06-27) — ERROR INJECTOR BUILT (M1-local DONE; fitting/validation SHERLOCK-gated)

The FIRST TASK (the error injector the measurement proved is needed) is BUILT and
composition-verified on the M1; the Sherlock fitting/validation is turnkey. **No member
code started** (correctly — that is the next gated cycle). Regression gate GREEN.

**DONE this session (all benchmark-only paths; gate stays green):**
- `scripts/benchmark/sim/error_injector.py` — the 3-layer generative injector: (1) per-read
  multiplier (Gamma/mixture) → over-dispersion, (2) 2-state cold/hot burst HMM (marginal-preserving)
  → within-read clustering, (3) fat-tailed geometric indel-run pmf. `InjectorParams.null()` =
  backward-compat marginal Poisson. Source-DECOUPLED `inject(seq,...)` works on clean OR pbsim reads.
  `measure_error_structure` + a pysam BAM adapter (`events_from_bam_read` / pure
  `events_from_alignment`) compute the SAME 3 stats on real/pbsim/pbsim+injector.
- **TRUTH RULE enforced:** injected errors = a SEPARATE `List[ErrorEvent]` track, NEVER `IndelTruth`
  (the SESSION-2 `has_unexplained`-zeroing lesson). The track IS PI-#2's hotness label.
- **Calibration = coordinate descent on the COMBINED output** (NOT independent MoM — layers
  cross-contaminate; burst inflates dispersion so Gamma shape converges ~1.26). Robust across seeds:
  disp 9.45–10.24 / gap5x 2.71–2.95 / run≥2 0.40 (targets 9.95/2.83/0.39).
- `tests/test_error_injector.py` (8 fast tests, all PASS run directly — no local pytest env) +
  `error_injector.self_check` (9s; the local composition gate, NOT realism — labeled as such).
- **PI-#2 MECHANISM probe** `scripts/benchmark/read_reliability_probe.py`: over-dispersion drives
  exonic→junction predictability (r 0.955 vs null −0.02), burst-only stays local (r 0.011) → the
  reliability covariate has signal but is imperfect → SOFT down-weight not hard filter (PI refinement
  (a) CONFIRMED). FDR-lift NUMBER deliberately NOT computed (SIRV-gated/magnitude-sensitive).
- `scripts/benchmark/error_realism_validate.py` — turnkey Sherlock `measure-bam`/`inject-fastq`
  (full protocol in its docstring). `inject-fastq` verified on a synthetic FASTQ; **`measure-bam`
  NOT VERIFIED on a real BAM** (Sherlock down) — only the pure `events_from_alignment` is unit-tested.

**KEY DISCIPLINE (advisor-reinforced, don't regress):** the M1 self-check proves COMPOSITION, not
realism. Realism = SIRV absolute-truth magnitude + DISTRIBUTIONAL comparison, both Sherlock-gated.
Magnitude is UNKNOWN (2.83× is an upper bound — read-vs-ref conflates RNA-mod + minimap2 pile-up).

**RESUME (session-5 boundary) — Sherlock auth was DOWN (Permission denied gssapi); ask Kevin to
re-auth (`ssh sherlock true`), do NOT re-open the master. Nothing in flight. Then, in order:**
1. **Re-sync code to Sherlock** (worktree changed): `rsync -az scripts/benchmark/ sherlock:$D/scripts/benchmark/`
   where `D=/home/groups/larsms/users/kevinroy/aligner_bench_live`.
2. **Re-profile real BAM with the burst gate OFF:** `empirical_cigar_error_profiler.py --isolation-flank 0`
   (default is already 0) on `/scratch/users/kevinroy/rectify_wt_by4742_rep1_*/...minimap2.namesorted.bam`.
3. **Validate pbsim+injector vs real** (the FIRST-TASK validation): `error_realism_validate.py measure-bam`
   on {real, pbsim round-trip, pbsim+injector (via `inject-fastq` then re-align)} all with `--thin 0.019`;
   success = +injector closes the pbsim→real gap on disp/gap5x/run≥2.
   **4 GUARDS (advisor, the self-checks are circular — these make validation MEANINGFUL; full text in
   SPEC §INJECTOR BUILD "GUARDS"):** (i) RE-DERIVE targets with THIS module's `measure_error_structure`
   (the hardcoded 9.95/2.83/0.39 are from the LOST err_corr.py — don't trust); (ii) LENGTH-CONTROL the
   comparison (dispersion ≈ linear in read len → calibrate to real reads' length dist, not flat 600bp);
   (iii) verify the two measurement paths AGREE (inject→align→`measure-bam` reproduces injected-truth
   stats — `events_from_alignment` can emit ins+sub at the same ref pos → gap=0, which `inject` never
   does); (iv) derive `frac_sub/ins/del`+`base_rate` from the real profiler (DRS is deletion-dominant,
   not the PLACEHOLDER 0.55/0.15/0.30).
4. **Re-fit `calibrate_params` to SIRV ABSOLUTE-truth targets** (the real magnitude; not the read-vs-ref
   upper bound) — needs the SIRV/LRGASP reference + the exon-GTF `gff_panel` loader variant (still unbuilt).
5. **THEN PI-#2 FDR-lift on stratum (G):** inject hot reads into JUNCTION_DISCOVERY, measure novel-junction
   FDR with vs without hot-read down-weighting. Prove the lift BEFORE shipping the member feature.
6. C1 member code remains the NEXT GATED CYCLE — not a benchmark session.
Regression gate (any session): `PATH=... pysam-python scripts/benchmark/smoke_roundtrip.py --out /tmp/x
--reps 20` → exit 0; plus `pysam-python scripts/benchmark/sim/error_injector.py` → exit 0 (composition).

---

## SESSION-4b (2026-06-27) — REDIRECT: in-silico error-REALISM first (real-data = calibrator, deferred)

User redirected: prioritize **in-silico** synthetic reads (errors injected into perfect
reads = absolute truth) and get the ERROR MODEL right, BEFORE real-data. Real/spike-in
data is reframed as the **error-process CALIBRATOR** (the only contamination-free way to
measure real error correlation structure), not the primary truth set. Two user priorities:
- **(4, most important) junction-discovery bias stratum** — canonical/non-canonical ×
  annotated/unannotated, to measure which aligners over-SNAP to GT-AG/annotated vs
  reliably discover de-novo introns. PROBED + CONFIRMED: minimap2 -ax splice -uf snaps
  non-canonical junctions to nearby canonical motifs (canonical GT-AG=EXACT; CA-TC→shift
  (197,397); AT-AC→(205,400); even with -C0), and snaps even ERROR-FREE → member-addressable
  motif-bias. NOT yet built (design in this session's transcript: gen_junction_discovery_stratum
  + scorer `fp_canonical_snap` + smoke (G)).
- **error-distribution realism** — we have MARGINAL rates only; missing (2) spatial
  clustering/burst + (3) per-read over-dispersion. CONFIRMED IN CODE: the empirical profiler's
  `--isolation-flank 10` recipe (its documented usage) counts only errors with ≥10 exact
  matches both sides → by construction EXCLUDES clusters/bursts. So an injector driven by it
  under-produces the bursty/hot-read regime where the convergent-error problem lives.

**ERROR-STRUCTURE MEASUREMENT: DONE** (was in-flight; result in SPEC §ERROR-REALISM).
VERDICT: pbsim3 is MORE iid/uniform than real ONT on every axis (rate-matched to the real
1.9%). Injector needs THREE layers: (1) per-read over-dispersion (real disp-index 9.95 vs
pbsim 3.55), (2) within-read burst/clustering (real 2.83x excess sub-5bp gaps vs pbsim ~1.2x),
(3) longer multi-base indel runs (pbsim 19% vs real 39% ≥2bp). Magnitude must be calibrated
vs SIRV ABSOLUTE truth (real-yeast 2.83x is an upper bound — read-vs-ref conflates RNA-mod +
alignment artifacts a pre-alignment injector must not reproduce). Injector BUILD = next cycle.

**POINT 4 (JUNCTION_DISCOVERY) — DONE** (commits 6c40f60 + 478ad88). gen_junction_discovery_stratum
(canonical/non-canonical × annotated/unannotated) + scorer `fp_canonical_snap` + smoke (G).
VERIFIED: minimap2 recovers canonical at ceiling (1.000), SNAPS non-canonical (recall 0.000,
25 snaps), and addressability is PROVEN not asserted (all 40/40 snaps cost NM>=1, so the true
non-canonical site strictly wins on a motif-blind score → evidence-weighing member recovers it;
the smoke fails if >10% snaps are NM==0 = the zero-evidence trap). Annotation axis flat for
minimap2 (annotation-agnostic); moves uLTRA/gmap in the panel run.

**DONE this session-4b:** junction-snapping probe (above); empirical-table clean-flank bias
confirmed in code; LRGASP URL-verified acquisition plan (committed to SPEC §EXTERNAL-VALIDITY,
the pbsim3-vs-NanoSim-vs-real-SIRV three-way; RNA002 so it answers sim-realism, LongBench=RNA004
for transfer); the 4-bug red-team fixes + error-realism framework.

**RESUME (session-4b — at the HANDOFF BOUNDARY; steps 1-2 DONE this session):**
The gate phase is complete + documented. Remaining work is fresh-session/next-cycle:
1. **Error injector build** (highest priority — the error-realism verdict above). Add to the
   in-silico path: (a) per-read rate multiplier (over-dispersed mixture), (b) self-exciting /
   Markov BURST process, (c) longer multi-base indel-run model — on top of the marginal table;
   HP-only marginal stays the backward-compat default. CALIBRATE the 3 params against SIRV
   ABSOLUTE-truth error structure (read-vs-known-sequence; NOT read-vs-reference). Measure-first
   (the err_corr.py approach on Sherlock) before/after to confirm pbsim+injector matches real.
2. **EXON-feature GTF loader** (`gff_panel` variant: group exon rows by transcript_id) — the
   shared dep unlocking ALL real-data integration (SIRV + LRGASP-GENCODE are exon-GTF). Small.
3. **Real-data acquisition** (cluster, never M1): LRGASP three-way (SPEC §EXTERNAL-VALIDITY has
   curl-verified URLs, pilot ~38GB; pbsim3-vs-NanoSim-vs-real-SIRV, DISTRIBUTIONAL comparison,
   use the with_mm_M27 GTF) for sim-realism; LongBench RNA004 (s3://longbench-data/) for the
   current-chemistry transfer check. Both gated on Sherlock auth being up.
4. **C1 member code** = the NEXT GATED CYCLE (do NOT start in a benchmark session): arm-LAW vs
   arm-flat ablation via `align_exon_block_global(..., penalty_table=...)`; plug-in point is
   `run_flat_affine_arm`/`cigar_records_to_bam` in the smoke.
5. **Open red-team hardening** (from the RED-TEAM section): non-canonical/minus-strand smoke
   coverage, (E) addressability, discriminating-scale smoke, Tier-2 indel-projection coords,
   verify the reverse-strand MAF risk on Sherlock.
Smoke regression gate (any session): `PATH=... pysam-python scripts/benchmark/smoke_roundtrip.py
--out /tmp/x --reps 20` → exit 0 = (A)(A2)(B)(B2)(C)(C')(D)(E)(F)(G) green.

---

## PI DESIGN NOTES (2026-06-27) — assessed + AGREED; fold into the error-injector + C3/§8 cycles

Two PI ideas off the error-distribution measurement. Both agreed (with refinements); both are
benchmark-TESTABLE (keep the prove-don't-assert discipline).

**(1) Model the error CORRELATION structure, not just marginals (the injector model).**
AGREE — but it is NOT more columns on the marginal table (clustering + per-read hotness are
sequential/joint, not per-position). Build a lightweight generative ERROR-PROCESS model:
  - per-read rate MIXTURE → over-dispersion (the hot-read tail; real max/median ~13x vs pbsim ~4x);
  - a 2-3-state BURST HMM / self-exciting process → within-read clustering (real 2.83x excess
    sub-5bp gaps vs pbsim ~1.2x);
  - a fat-tailed INDEL-RUN-LENGTH distribution (real 39% vs pbsim 19% indels ≥2bp).
  Fit by EM/HMM/method-of-moments — **dependency-light, NOT a neural net** (project no-heavy-deps
  rule; escalate only if the simple model fails). **CRITICAL DATA CAVEAT:** the current empirical
  table CANNOT see bursts — it was built with `--isolation-flank 10`, which excludes clustered
  errors by construction. Re-profile with the gate OFF, ideally on **SIRV absolute-truth** reads
  (so real variants / alignment artifacts don't masquerade as bursts). This model feeds the
  injector (benchmark realism) AND enables (2).

**(2) Per-read "hotness" as a NOVEL-JUNCTION-DISCOVERY FDR signal (a member/§8 facet).**
STRONGLY AGREE — likely the higher-payoff idea. Mechanism: a globally-hot (error-riddled) read's
junction region is also likely error-riddled → it can map "best" to a spurious UNANNOTATED junction
by chance → counting it as support for a novel junction INFLATES discovery FDR. The enabler: a
read's hotness is ESTIMABLE FROM ITS EXONIC PORTION (error density vs reference in well-anchored
exon blocks) WITHOUT junction truth → a self-calibrating per-read RELIABILITY covariate. Works
because of the measured per-read over-dispersion. Refinements: (a) SOFT down-weight / posterior
input, NOT a hard filter — the localized-burst component means a read can be clean in the exon yet
bursty at the junction, so exonic density predicts junction reliability only probabilistically
(fits §8.2 "soft prior, never hard gate"); (b) slots into §8 discovery-FDR + the C3 calibrated
posterior as a NEW ORTHOGONAL read-reliability signal, complementary to the abstain band + anchor
gate + canonical/non-canonical tracks. **BENCHMARK TEST (build it):** once the injector makes hot
reads, measure (i) does exonic error density predict junction-region error (within-read error
autocorrelation), and (ii) does down-weighting hot reads improve novel-junction FDR on the
JUNCTION_DISCOVERY stratum (G)? The over-dispersion injector layer + stratum (G) are exactly the
pieces to test it — do NOT ship the signal as a member feature until the benchmark shows the lift.

**Thread connection:** (1) makes hot reads exist in the benchmark → which lets us validate (2).
So sequence: injector over-dispersion layer first, then the read-reliability FDR test, then (if it
shows lift) the member facet.

---

## SESSION-4 (2026-06-27) — chosen next: real-data transfer check (LongBench); BLOCKED on Sherlock re-auth

User chose to run the real-ONT external-validity transfer check (does the pbsim3
Tier-2 junction recall/FDR transfer to real reads?). Started recon; then Sherlock's
ControlMaster died (daily Kerberos/2FA lapse at the date boundary) — **needs Kevin
to re-auth** (`ssh sherlock true`); do NOT re-open the master yourself.

**Recon DONE (from M1 via public-S3 HTTP listing — no aws cli on M1):**
- LongBench bucket `s3://longbench-data/` (ap-southeast-2, --no-sign-request).
  Reads in `raw/fastq/`: per cell line (H146,H1975,H211,H2228,H526,H69,HCC827,SHP77)
  `{S}_dRNA_ONT.fastq.gz` = **RNA004 DRS** (~12–28 GB ea), `{S}_bulk_ONT.fastq.gz` =
  cDNA (~25 GB), plus PacBio/Illumina. A `fastqs_topup/` dir has more dRNA.
- ⚠ **OPEN/UNVERIFIED:** the bucket README does NOT mention SIRV/Sequin spike-ins and
  no reference/annotation/SIRV/sequin keys exist in the bucket. MUST verify spike-ins
  are actually present in these reads BEFORE downloading 16+GB (align a small sample
  to a SIRV/Sequin reference and check for mapped reads). The SIRV-Set 4 + Sequin
  reference (FASTA+GTF) must come from Lexogen (SIRVsuite / lexogen.com) + Garvan
  (sequinstandards.com), or LRGASP Synapse syn25683367/syn25683630.

**RESUME (real-data transfer check):**
1. Confirm Sherlock is back: `ssh sherlock true` (Kevin re-auths if it 2FA-prompts).
2. Build the SIRV/Sequin loader FIRST (M1-local, no cluster): extend `gff_panel.py`
   with an exon-feature GTF loader (`build_panel_from_gtf`: group exon features by
   transcript_id, sort, derive introns) — SIRV/Sequin are standard exon-GTF, unlike
   the yeast mRNA+intron GFF. Test on a synthetic exon-GTF.
3. On Sherlock (NOT M1 for the reads): fetch the SIRV-Set 4 + Sequin reference;
   `aws s3 cp --no-sign-request s3://longbench-data/raw/fastq/H146_dRNA_ONT.fastq.gz`
   (+ a cDNA `H146_bulk_ONT.fastq.gz`) to scratch; align to the SIRV/Sequin reference
   with minimap2 -ax splice (endogenous human reads won't map → mapped = spike-in
   subset); build truth via the GTF loader; `score_bam` → junction recall/FDR.
4. COMPARE to pbsim3 Tier-2 (DRS 0.816 / cDNA 0.843). Small gap = sim transfers;
   large gap = error model needs re-tuning/caveat. Use `sherlock-sbatch`; heavy →
   cluster only, never relay reads through M1.
Full dataset plan + ranked alternatives (SG-NEx, LRGASP): `SIMULATION_BENCHMARK_SPEC.md`
§"EXTERNAL-VALIDITY DATA PLAN". Nothing is mid-flight (no running job); marker cleared.

---

## SESSION-3 ADDENDUM (2026-06-26) — VARIANT/C6 stratum BUILT + verified discriminating

Continued the next-increment strata. Per advisor triage, led with **VARIANT/C6**
(the only remaining stratum that is fully wired to a scorer metric the incumbent is
below ceiling on); the other three are deferred/scaffold with reasons recorded.

- **`gen_variant_stratum` in `scripts/benchmark/sim/controlled.py`** (wired into
  `generate_corpus`, `include_variant=True`) + smoke assertion **(E)**.
- **Discriminating construct, verified EMPIRICALLY vs minimap2 -ax splice -uf:** a
  GT..AG-flanked **deletion variant ≥ ~40bp** is re-expressed by minimap2 as a
  spurious **intron (N-op)** → a FABRICATED **variant-adjacent** junction FP
  (`fp_variant_adjacent`) + the deletion mis-scored. Probed thresholds: 20/30bp
  stay `D` (correct), 40/60/100bp all become `N`. Smoke reps=20: driver
  `fp_variant_adjacent=60`, driver junction FDR=1.0.
- **SPECIFICITY proven (not just sensitivity):** SNP 3bp from donor, SNP ~100bp
  away, and a 25bp deletion are ALL handled correctly by minimap2 (0 FP). Smoke (E)
  asserts controls = 0 variant-adjacent FP, so the stratum shows the FDR is
  SPECIFIC to the splice-mimic-deletion context — guarding the future C6 member
  from a blunt abstain-near-every-variant rule. het/hom + non-Mendelian VAF (0.33)
  recorded per `VariantTruth`.
- Whole-corpus (C) indel concordance moves 0.999→0.940 — EXPECTED + correct: the 60
  splice-mimic deletions minimap2 turned into introns are now counted as genuine
  incumbent indel failures the benchmark measures (not a regression; (C) only
  requires ≥1 correct). Schema variant round-trip is exercised end-to-end (the
  `fp_variant_adjacent` count depends on variants surviving the TSV reload).
- **Anti-overcount (advisor-caught):** each VARIANT read now gets its OWN
  freshly-randomized contig carrying the same variant KIND (HP/STR already vary
  per-read; VARIANT had been replicating ONE construct `reps` times → effective
  n=#sub-cases). Now `3×reps` INDEPENDENT driver constructs; region-disjoint split
  preserved (per-contig). Verified: 60 distinct drivers all fabricate, 60 distinct
  controls give 0 variant-adjacent FP. Caveat recorded: 25bp→D / 40bp→N threshold
  is pinned to minimap2 2.28 `-k 14` — re-probe on a version/flag change.
- **PARALOG/C4 also landed this session** — see DONE/OPEN below + SPEC §Triage.
- **Smoke GREEN** with assertions (A)(A2)(E)(F)(D)(C)(C')(B)(B2) (`--reps 20`).

### SESSION-3 TIER-2 ADDENDUM (2026-06-26) — transcriptome run drafted, verified, SUBMITTED
The Tier-2 external-validity tier (Branch A) is built, verified on REAL yeast
coords, and submitted to Sherlock.
- **Aligner inventory decided the branch** (rectify env, verified): minimap2 +
  uLTRA + deSALT (binaries) + gapmm2 + mappy present; **gmap absent** (lives in
  `aligner_bench`/`apanel`/`compass` envs), **mapPacBio absent**, spoa/abpoa absent.
  None wired into the benchmark driver. So a true 5-aligner panel TAIL is NOT
  runnable today → **Branch A = minimap2 baseline** (honest labels in the JSON);
  **Branch B (panel + hard-read injection) enumerated, gated.**
- **NEW code:** `scripts/benchmark/sim/gff_panel.py` (GFF→TranscriptModel loader;
  yeast has no `exon` feature → exons = mRNA span minus `intron` features),
  `scripts/benchmark/tier2_run.py` (driver), `pbsim3_wrapper` annotation-aware
  classification + the projection fix below. `run_tier2.sbatch` (DRS+cDNA).
- **Projection bug FOUND+FIXED during real-coord verification** (the reason to
  verify before cluster spend): `project_maf_record` assigned EVERY transcript
  intron to EVERY read; pbsim reads are fragments → false FN, recall deflated
  0.79. Fixed to per-read spanned+anchored (≥10bp flank) junctions only → 0
  out-of-span junctions; real minimap2 DRS ANNOTATED recall ≈0.80–0.85 (short
  yeast introns are genuinely hard — minimap2 calls some as deletions), FDR
  ≈0.03–0.07. **This is the saturation control (harness validation), NOT
  discrimination.**
- **Job 31628436 COMPLETED (2m17s, rc=0).** DRS: 9810/10000 placed, ANNOTATED
  recall=0.816 FDR=0.051. cDNA: 9976/10000 placed, recall=0.843 FDR=0.077.
  Summaries in `$D/tier2_results/{drs,cdna}_summary.json` (with `_scope` labels);
  numbers + reading recorded in `SIMULATION_BENCHMARK_SPEC.md` §"TIER-2 BRANCH-A
  RESULT". Saturation control PASSED (harness validated on real coords at scale);
  ~0.82–0.84 recall is real minimap2/short-yeast-intron behavior, not artifact.
  (Infra: landed on AVX-512 node sh03-07n10 despite --exclude but the rectify env
  did NOT SIGILL — exclude was harmless-but-not-honored; use owners for Branch B.)
- **Remaining strata deferred with reasons** (PARALOG needs a scorer locus-
  concordance metric first; COVERAGE×Q can't discriminate until C3 consumes phred;
  PANEL_FAILURE's real validity is the Tier-2 5-aligner panel) — see
  `SIMULATION_BENCHMARK_SPEC.md` §"VARIANT/C6 stratum … Triage".

---

## SESSION-2 ADDENDUM (2026-06-26) — gate was NOT green on a clean checkout; now is

The prior VERIFIED claim "GATE smoke is GREEN" **overclaimed**: it was green only
because that session's `.fai` happened to be fresh. On a clean re-run the smoke
**crashed** (`KeyError 'hph_A_04'` in `run_flat_affine_arm`). Two fixes landed +
the cell-size MUST-HAVE was audited:

- **`021ef97` stale-.fai guard** — `pysam.FastaFile` reuses a `.fai` older than its
  FASTA (no mtime check). The benchmark regenerates `ref.fa` in place; when the
  HP_HARD (`hph_*`) contigs were added, a leftover 61-contig index from a prior run
  exposed only the OLD contig set. `score_bam` tolerated it (`genome.get(chrom,'')`
  masked the miss) but `run_flat_affine_arm` did `genome[t.chrom]` → crash, so the
  gate's two-arm DP path **never actually ran clean**. New `scorer.open_fasta()`
  drops a stale `.fai` so pysam rebuilds; routed `load_genome` + `cigar_records_to_bam`
  (header `get_tid` would also silently truncate) + the smoke's shifted-junction
  opener through it. Verified by a targeted guard test (rewrite FASTA + new contig
  with a backdated `.fai` → guard rebuilds, plain pysam does not).
- **`c680906` missing-contig now LOUD** — `genome.get(chrom,'')` silently disabled
  `normalize_junction` (junctions compared RAW = the GMAP-0.09 FP the gate prevents).
  Now counted (`reads_missing_contig`), warned once, surfaced in `summary()`.
- **Cell-size MUST-HAVE AUDITED at `--reps 400`** (the documented gate-validity bar):
  HP indel cells span **lengths 2–12 × {A,C,G,T}** (min n=183), HP_HARD {4,6,8,10,12}
  (n=400), STR ≥155 — **all clear `min_count=100`**. **Length 1 has no deletion cell
  by construction** (`_draw_k(L=1)`→0; a length-1 run can't undercall its length) but
  is fully present as **clean reads** (400/base — the FP/false-indel control). So
  `reps=400` produces a VALID C1 corpus. (Default `--reps 120` is fine for the smoke
  but the C1 ablation corpus MUST use ≥400.)
- **Smoke re-verified GREEN** after both fixes (`--reps 20`, all GATE assertions incl.
  the (D) two-arm HP_HARD-noisy DP path that was crashing).
- **`9f71770` pbsim3 LIVE round-trip — now VERIFIED on Sherlock** (was the top OPEN
  item). The fresh `pbsim3` conda env built clean (rc=0; `ERRHMM-ONT.model` DRS +
  `ERRHMM-ONT-HQ.model` cDNA present). The live run exposed **3 real mechanical bugs
  in `pbsim3_wrapper.py`** that the synthetic-MAF unit test could not — all fixed:
  - this pbsim3 build emits **gzipped single-file** outputs (`sim.maf.gz`/`sim.fq.gz`),
    not `<prefix>*.maf`/`*.fastq` → MAF discovery + FASTQ concat found nothing → 0
    truth. `parse_maf` now gz-aware; `run_pbsim3` matches `.maf(.gz)` and **RAISES**
    (not silent 0) if none; concat handles `.fq/.fastq(.gz)`.
  - templ mode emits **ONE read per template RECORD** (`--depth` does NOT multiply
    count; `--pass-num` does). The new live driver replicates each gene into N
    records → N reads/gene.
  - New `scripts/benchmark/sim/live_roundtrip.py` (self-contained ±strand spliced
    genes → pbsim → minimap2 → scorer). **VERIFIED end-to-end:** pbsim ran, gzipped
    MAF parsed, **166/180 reads placed, junction recall=0.667 FDR=0.277 on real
    pbsim reads, 0 missing contigs.**
  - **Indel exact-concordance is REPORTED but NOT gated** — it is a TIER-1 metric.
    pbsim's MAF encodes every stochastic error as an indel, so the projected per-read
    "indel truth" is thousands of scattered edits; minimap2 redistributes them and the
    scorer's per-read `has_unexplained` gate (built for Tier-1's single known indel)
    zeroes the read. This is the SPEC's two-tier split working as designed (Tier-2 =
    junction recall/FDR; Tier-1 = exact-indel), NOT a defect.

---

## RED-TEAM (2026-06-26, 4 independent auditor agents) — gate is TRUSTWORTHY, no blocker
4 read-only auditors swept truth-construction (Tier-1 + Tier-2), the scorer, and
gate-validity/claims. Verdict: fitness genuinely scores vs truth (never the internal
score), ambiguity-aware match works, the flat-affine ablation path runs. **4 verified
bugs FIXED** (commit `2f80ee5`): VARIANT eq-span too narrow (193/480, ambiguity-span
class); scorer insertion right-boundary FN; `fp_variant_adjacent` donor-only (missed
acceptor); HP `true_cigar` missing `-k` (latent).

**OPEN red-team findings (NOT yet fixed) — for the next hardening pass:**
- **Coverage: non-canonical + minus-strand FDR tracks UNEXERCISED.** Whole corpus is
  `+`-strand/canonical, so the §8 non-canonical track and the minus-strand canonicity
  fix are never hit by the smoke (classifier verified correct, but a regression would
  pass green). Add a non-canonical-NNC + a minus-strand spliced locus to the smoke.
- **(E) VARIANT proves discriminating+specific but NOT addressable** (the bar (D)/(F)
  meet). A ≥40bp-deletion-as-intron may leave no pileup signal a C6-MVP emission could
  consume — add an addressability proof or downgrade the claim.
- **reps=20 smoke is not discriminating-scale:** (D)'s C1-addressability verdict is
  decided on ~3 reads (fragile); `min_cell_reads`=9 vs the 100 floor; "cell-audited"
  has no CI guard. Make (C') a hard assert at C1 scale; run (D) at a reps with a
  non-trivial failure denominator.
- **Tier-2 indel-projection coords (REPORT-ONLY, not gated):** insertion off-by-one +
  not coalesced; minus-strand multi-base UNIQUE deletion `eq` anchored at the wrong
  end (`pbsim3_wrapper.py:142-170`). Fix for truth.tsv honesty.
- **UNVERIFIED risk (verify on Sherlock):** Tier-2 reverse-strand MAF projection
  assumes pbsim keeps the template 's' line forward; if pbsim revcomps it, `-`-read
  JUNCTIONS (gated) would be wrong. Check a real `.maf` block: template strand for a
  `-` read line.
- MINOR: `reads_missing_contig` still folds raw-coord TP/FP/FN into totals (warns but
  doesn't exclude); `locus_accuracy` diluted by unique-contig strata; STR AT-only
  cells <100 at reps=120; SPEC strata table still calls STR discriminating.

## DONE
Three ABSENT components built (everything else is reuse), all under the
benchmark-only paths the brief allows (`rectify/core/benchmark/`, `scripts/benchmark/`):

1. **Truth-propagation schema (built FIRST — it constrains the other two)** —
   `rectify/core/benchmark/truth_schema.py`. Per-read `ReadTruth` with:
   - **exact-indel truth** as a per-position `IndelTruth` carrying an
     ambiguity-equivalence span `[eq_start, eq_end)` (HP run / STR period / unique)
     — the net-indel-in-run model, so the framing metric is ambiguity-aware;
   - **junction truth** `JunctionTruth` left-normalized via
     `chimeric_consensus.normalize_junction`, with canonicity + **NIC/NNC class**
     (ANNOTATED / NIC / NNC) and the ambiguity window;
   - **CPA truth** + downstream-A context (C2); **standing-variant truth** (C6);
   - the **shared genomic-region-disjoint TRAIN/TEST split tag** (region-decided);
   - C1 HP/STR cell metadata (`base_class`, `run_copies`) for the `min_count=100` audit.
   - Lossless TSV round-trip (`write_truth_table`/`read_truth_table`). NOT overloaded
     onto the `XV` tag (a sidecar table, per the brief).

2. **Read-simulator wrapper (pbsim3)** — `scripts/benchmark/sim/`:
   - `transcript_model.py` — the truth-propagation backbone: `TranscriptModel`
     (exon structure -> spliced FASTA, `transcript_pos_to_genome` map, introns,
     NIC/NNC `junction_truths`). Junction/CPA/NIC-NNC truth is derived from THIS
     construction, NOT from the simulator.
   - `pbsim3_wrapper.py` — Tier-2 realistic reads. **pbsim3 chosen because it is
     the ONLY one of pbsim3/badread/nanosim that emits a per-read MAF
     (read<->template alignment)**, which the EXACT-INDEL framing metric requires;
     badread/nanosim give origin+identity but no per-read edit script. Parses the
     MAF and composes (read->transcript) ∘ (transcript->genome) into exact
     `ReadTruth` (genome CIGAR + indels + junctions). `--method errhmm`;
     ERRHMM-ONT model = DRS, ERRHMM-ONT-HQ = PCR-cDNA.
   - `controlled.py` — **Tier-1 controlled-error generators (the DISCRIMINATING
     tier)**: HP (A/C/G/T x 1-12, del-dominant), STR (ED-tied placements),
     JUNCTION_AMB (a GT-AG intron whose donor sits in a repeat → ambiguity test),
     CLEAN (false-indel control). Emits ref.fa + reads.fastq + truth.tsv.

3. **Ambiguity-aware per-junction + per-indel scorer** —
   `rectify/core/benchmark/scorer.py`. Reuses
   `chimeric_consensus.normalize_junction / junction_ambiguity_window /
   _canonical_within_window` so a call 1bp into a donor/acceptor repeat is NOT
   charged FP. Junction TP/FP/FN stratified by class (ANNOTATED/NIC/NNC) AND
   canonicity (separate FDR tracks, §8); indel **position-exact concordance** +
   false-indel rate stratified by C1 cell; CPA `|est-true|`; variant-adjacent FP
   tagging (C6); `score_panel` fills which-aligners-placed + panel-unplaced fraction.

4. **End-to-end smoke** — `scripts/benchmark/smoke_roundtrip.py`.

## VERIFIED — the gate is BUILT, CONTAINS a C1-addressable stratum (validated vs
## TRUTH), and the internal-DP ablation harness is RUNNABLE. (Proving the length-law
## CLOSES the gap — arm-LAW > arm-flat — is the next-cycle ablation; see OPEN.)
- **GATE smoke is GREEN** (`scripts/benchmark/smoke_roundtrip.py`, minimap2 -ax
  splice on the Tier-1 corpus):
  - (A) known junction round-trips truth->reads->minimap2->scorer as **NNC TP=30**;
  - (A2) **NIC TP=30 and ANNOTATED TP=30** — the discovery-class classifier is
    verified end-to-end (not just present); separate FDR tracks exercised;
  - (B) a deliberately 1bp-shifted junction call (into the repeat) scores
    **TP=1 FP=0** (ambiguity-aware match, the GMAP-0.09 trap defended);
  - (B2) NEGATIVE control: a junction shifted 40bp (beyond the window) is correctly **FP**;
  - (C) indel **position-exact concordance = 0.997**, false-indel-rate=0.0;
  - (D) **the gate CONTAINS a C1-addressable stratum, validated against TRUTH:**
    the live flat-affine DP (`align_exon_block_global`, the arm C1 upgrades) scores
    **0.980 on HP_HARD-noisy (BELOW ceiling vs truth)**, and **24/24 (100%) of its
    failures are systematic out-of-run deletion misplacement OR a substitution
    "repaired" with a spurious indel** — exactly the indel-vs-substitution tradeoff
    C1 targets (a run/length-aware cost keeps the deletion in-run and stops the
    mismatch-repair). The internal-DP ablation path RUNS — the DP is BAM-ized
    (`scorer.cigar_records_to_bam`) and scored on 2400 reads — so the arm-LAW vs
    arm-flat comparison is runnable; isolated-HP control = 1.000 (the SPEC
    VERTICAL-SLICE FINDING), STR = 1.000.
- **KEY method note (advisor-driven):** minimap2 and the flat-affine DP are the
  SAME error family, so they AGREE by construction (boundary_sub = 1.000 == 1.000
  after the truth-corruption fix below). C1-discrimination is therefore decided vs
  **TRUTH** (incumbent below ceiling + failures C1-addressable, both shown above),
  NOT by a flat-vs-flat separation. Whether the length-law CLOSES the measured gap
  is the next cycle's ablation, not claimable here.
- **Two scorer correctness fixes:**
  1. a per-indel TP now requires net-in-span == truth net AND no unexplained indel
     OUTSIDE every truth span (the vertical-slice ``out_run == 0`` rule); insertion
     off-by-one made consistent with the half-open span test;
  2. **indel ambiguity span now uses the base-equality SLIDE** (`truth_schema.
     deletion_ambiguity_span`, the same primitive style as the junction slide), not
     just full-period run detection — this captures HALF-PERIOD bleeds (a 'T' before
     an (AT)n run) that were charging ambiguity-EQUIVALENT STR placements as false
     errors (STR 0.91 -> 1.000). This separated the real C1 signal (HP_HARD out-of-
     run misplacement) from the scorer artifact (STR bleed).
- **Truth-corruption fix:** the boundary substitution in HP_HARD is never allowed to
  become the run base (it would silently change the effective run length).
- Schema lossless TSV round-trip verified (junction normalization, canonicity,
  HP-cell accounting, variant/split round-trip).
- pbsim3 MAF→genome **projection** verified locally on a synthetic MAF AND now on a
  **LIVE pbsim3 run on Sherlock** (session 2, `9f71770`) — see the addendum above.

## OPEN
- **pbsim3 live run: DONE (session 2, `9f71770`).** Fresh `pbsim3` conda env on
  Sherlock (rc=0). Live round-trip VERIFIED via
  `scripts/benchmark/sim/live_roundtrip.py` (166/180 reads placed, junction
  recall=0.667 FDR=0.277, 0 missing contigs); 3 wrapper I/O bugs found + fixed. See
  the addendum for the full result + the Tier-1-vs-Tier-2 indel-metric caveat.
  Code is synced to Sherlock at
  `/home/groups/larsms/users/kevinroy/aligner_bench_live/` (run with the `rectify`
  env python + `PYTHONPATH` there + pbsim3-env binaries; see RESUME).
- **Tier-1 cell sizing: RESOLVED (audited session 2).** `--reps 400` clears
  `min_count=100` for every indel cell (HP lengths 2–12, HP_HARD {4,6,8,10,12}, STR)
  and every clean cell incl. L=1 (the FP control). Use `--reps 400` for the C1 corpus.
- **NEXT-CYCLE ABLATION (gated on this gate) — does the length-law CLOSE the gap?**
  The gate now CONTAINS a validated C1-addressable stratum (HP_HARD-noisy: flat-
  affine 0.980 vs truth, 24/24 failures C1-addressable). The remaining proof is the
  C1 cycle's job: build the length-law arm and show **arm-LAW > arm-flat** on
  HP_HARD-noisy (recovering those 24/24 out-of-run misplacements). The harness is
  ready: `run_flat_affine_arm` (smoke) + `scorer.cigar_records_to_bam` are the
  plug-in point — add a `run_lengthlaw_arm` that calls
  `align_exon_block_global(..., penalty_table=...)` once C1 Phase-1 exists.
- **STR is non-discriminating after the scorer fix** (1.000 — the earlier 0.57/0.91
  was the half-period-bleed scorer artifact, now corrected). If an
  edit-distance-TIED STR discriminator is wanted, construct a case where the true
  unit-deletion placement is NOT the leftmost slide (so left-alignment genuinely
  diverges from truth) — the current whole-unit deletion is slide-equivalent and
  correctly scores 1.0.
- **Tier-2 (Branch A, yeast minimap2 baseline): DONE — job 31628436 COMPLETED rc=0.**
  Results recorded in SPEC §"TIER-2 BRANCH-A RESULT" (DRS recall 0.816/FDR 0.051;
  cDNA 0.843/0.077). Re-run any time: `ssh sherlock 'cd <D> && sbatch run_tier2.sbatch'`.
  Branch B (the real tail + cross-aligner FDR) is the gated follow-on: NOT this
  job. It needs (1) the multi-aligner panel WIRED + index-prepped (minimap2 +
  gapmm2 + uLTRA + deSALT are in the env; gmap is in `aligner_bench`, mapPacBio
  absent) so `score_panel`'s `_panel_unplaced_fraction` is "placed by NO aligner",
  and (2) an INJECTED hard sub-population (elevated error / repeat) — a clean run
  reports tail≈0, which is a false negative, not "no C5 needed". Novel-junction
  (NIC/NNC) recall also needs isoform injection (exon-skip→NIC, novel-site→NNC);
  Branch A measures ANNOTATED-recall + spurious-FDR only.
- **standing-variant/C6: DONE (session 3)** — `gen_variant_stratum`, smoke (E),
  verified discriminating + specific (see session-3 addendum).
- **paralog/C4: DONE (session 3)** — added the scorer locus-concordance readout
  (`AlignerScore.locus_accuracy`/`locus_correct`/`locus_incorrect`/`locus_mapq0`)
  FIRST, then `gen_paralog_stratum` + smoke (F). **Advisor-corrected:** the first
  cut used a window-EXCLUDING (zero-evidence) fragment — below ceiling at 0.5 but
  UNRECOVERABLE by any method incl. C4 (the vertical-slice trap; "C4 could recover"
  was false). Fixed to a WEAK-evidence fragment (covers exactly ONE distinguishing
  SNP): per-read minimap2 below ceiling (locus_acc≈0.94) AND pooling-majority at the
  lone SNP recovers the true copy 6/6 pools — the gap is provably C4-addressable
  from truth (the (D)-equivalent). Spanning reads at ceiling 1.000 (control).
  Effective diversity = n_families (3); reps = per-locus depth (a C4-pooling
  feature). See SPEC §Triage PARALOG bullet.
- Remaining strata (panel-failure/C5, coverage×Q) DEFERRED with recorded reasons
  (Tier-2-only / next-cycle-consumer) — see `SIMULATION_BENCHMARK_SPEC.md`
  §"Triage of the remaining SPEC strata". Best next: GENOMIC_A_CPA/C2 (scorer
  already emits `|est−true CPA|`; needs a walkback comparator arm to fully
  discriminate) OR the real Tier-2 transcriptome run on Sherlock.
- No C1/member code built (correctly — that is the next, gated cycle).

## RESUME
Both gating items (pbsim3 live run + Tier-1 cell sizing) are DONE. The gate is
green, valid, cell-audited, and Tier-2-mechanically verified — a clean stopping
point. The next increment is the remaining strata + the real Tier-2 run.

0. **Local smoke = regression gate first (≈30s at low reps):**
   `cd <this worktree> && PATH="/Users/kevinroy/miniconda3/bin:/opt/homebrew/bin:$PATH" \
    /Users/kevinroy/miniconda3/envs/pysam/bin/python scripts/benchmark/smoke_roundtrip.py \
    --out /tmp/bench_smoke_chk --reps 20` (exit 0 = green). Full-scale is `--reps 120`
   (slow: ~2+ min — the two-arm DP runs the Python aligner on ~2400 reads).
1. **Re-confirm the live pbsim3 round-trip (regression, any time):**
   `ssh sherlock 'D=/home/groups/larsms/users/kevinroy/aligner_bench_live; \
    CB=/home/groups/larsms/users/kevinroy/anaconda3; P3=$CB/envs/pbsim3; cd $D; \
    PYTHONPATH=$D $CB/envs/rectify/bin/python scripts/benchmark/sim/live_roundtrip.py \
    --errhmm-model $P3/data/ERRHMM-ONT.model --pbsim-bin $P3/bin/pbsim \
    --minimap2 $P3/bin/minimap2 --samtools $P3/bin/samtools --out /tmp/pbsim_live_chk \
    --seed 7 --copies 60'` → exit 0 = green (junction recall reported; indel concordance
   is Tier-1-only, see addendum). NOTE: run with the **`rectify`** env python (full deps
   incl. numpy/pandas for `import rectify`), pbsim/minimap2/samtools from the **`pbsim3`**
   env. If the worktree code changed, re-`rsync` it first:
   `rsync -az --exclude data/ --exclude __pycache__/ --exclude '*.pyc' --exclude bin/ \
    rectify/ sherlock:$D/rectify/` and `rsync -az scripts/benchmark/ sherlock:$D/scripts/benchmark/`
   (mkdir the remote `$D/scripts` parent first — rsync won't create missing parents).
2. **NEXT INCREMENT — remaining SPEC strata** (M1-light Tier-1 controlled generators;
   schema already supports all): paralog/SMN, panel-failure/C5, coverage×Q, and the
   standing-variant **C6** generator (the one with first-class variant-adjacent junction
   FDR). Build in `scripts/benchmark/sim/controlled.py` (mirror the HP/STR/JUNCTION
   generators), add a smoke assertion per stratum, size cells ≥100 at `--reps 400`.
3. **NEXT INCREMENT — real Tier-2 transcriptome run** (Sherlock, `sherlock-sbatch`
   skill, owners partition, AVX-512 exclude; do NOT relay BAMs through the M1): drive
   `live_roundtrip.py`-style propagation over a real transcript panel (yeast saturation
   control + human SMN1/SMN2 + a NIC/NNC-rich set mirroring A549 chr5) for global
   junction recall/FDR + panel-failure tail sizing. Build the genome+transcript panel
   from `TranscriptModel`s off the real GFF; the projection path is verified.
4. **C1/member code is the NEXT GATED CYCLE — not here.** The arm-LAW vs arm-flat
   plug-in point is `run_flat_affine_arm` + `scorer.cigar_records_to_bam` (smoke);
   add `run_lengthlaw_arm` calling `align_exon_block_global(..., penalty_table=...)`
   once C1 Phase-1 exists.

## FILES (all NEW, benchmark-only paths — none touch shared hot files)
- `rectify/core/benchmark/__init__.py`
- `rectify/core/benchmark/truth_schema.py`   (component 2 — schema)
- `rectify/core/benchmark/scorer.py`          (component 3 — ambiguity-aware scorer)
- `scripts/benchmark/sim/transcript_model.py` (component 1 — truth backbone)
- `scripts/benchmark/sim/pbsim3_wrapper.py`   (component 1 — Tier-2 pbsim3 wrapper; gz/fq-robust)
- `scripts/benchmark/sim/live_roundtrip.py`   (Tier-2 LIVE mechanical-integration check — session 2)
- `scripts/benchmark/sim/controlled.py`       (component 1 — Tier-1 generators)
- `scripts/benchmark/smoke_roundtrip.py`      (end-to-end GATE smoke)
- `dev/SIMULATION_BENCHMARK_SPEC.md`          (UPDATED: simulator decision recorded)
