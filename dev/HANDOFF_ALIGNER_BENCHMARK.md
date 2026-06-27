# HANDOFF — Deliverable A: the simulation ground-truth benchmark (the GATE)

**Agent:** dedicated benchmark-builder (isolated worktree, branch
`worktree-agent-a25a2c1e784ad37dc`, based on `drs-validation-rebuild` so the reuse
primitives + design docs are present). **NEVER commit to `drs-validation-rebuild`.**
**Updated:** 2026-06-27 (session 4).

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
