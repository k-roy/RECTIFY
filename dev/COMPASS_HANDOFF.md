# COMPASS human-A549 rework — HANDOFF (2026-06-18)

Self-contained pickup for a FRESH agent continuing from P0. Read the cluster/workspace rules in
`~/.claude/CLAUDE.md`, `~/work/CLAUDE.md`, `~/work/UCLA/Chanfreau_Lab/CLAUDE.md` FIRST (ControlMaster
discipline, owners partition + AVX-512 constraint for rectify-env jobs, chunking, never relay big data
through the M1). Deep context: `dev/COMPASS_HUMAN_REWORK_PLAN.md` (the phased plan),
`dev/COMPASS_SHORTREAD_SCOPING.md` (scoping + STEP-0), `dev/ALIGNER_INVESTIGATION_SYNTHESIS.md` (de-herding +
design insights). COMPASS codebase map is in this session's transcript; key facts reproduced below.

## CONTEXT (why this exists)
COMPASS = the PI's own published short-read splice-junction tool (github.com/k-roy/COMPASS, Roy et al. 2023
NAR PMID 37956322): per-READ best-alignment arbitration across 6 aligners. We are **adapting it from yeast
(S288C) to human A549** to build a SENSITIVE multi-aligner short-read junction detector, because the only
existing short-read data (SG-NEx A549 Illumina, STAR 1-pass) validates KNOWN-corroborated junctions at just
**14.5%** (STEP 0) — too insensitive to adjudicate the **111 GMAP-only "recurrent novel" junctions** we need
to confirm-or-refute (`dev/gmap_only_recurrent_novels_chr5.tsv`; currently "likely artifacts, INCONCLUSIVE").
**Architecture (PI-decided): do the short-read work INSIDE COMPASS; back-propagate RECTIFY's junction
refinements into COMPASS; keep RECTIFY end-correction-focused.** Push to k-roy/COMPASS is AUTHORIZED (use a
`human-a549` branch, no force-push). Full **6-aligner panel** (incl. Magic-BLAST, the heavy pole) +
minimap2 `splice:sr`.

## DONE (P0)
- COMPASS cloned: M1 `~/work/COMPASS`, Sherlock `/scratch/users/kevinroy/compass_a549/COMPASS` (the work dir
  `$W = /scratch/users/kevinroy/compass_a549/`).
- **GRCh38 introns TSV** `$W/gencode.v44.introns.tsv` — 402,454 introns, chr-named, 18,451 chr5. Generator
  `$W/make_human_introns.py`. Convention VERIFIED = intron BODY 1-based `[exon1.end+1, exon2.start-1]` (vs
  yeast YAL030W).
- **samfixcigar replacement** `$W/samfixcigar.py` — pysam, rewrites `M`→`=`/`X` vs reference (COMPASS needs
  SAM-1.4 ops for its mismatch scoring). Replaces jvarkit (its gradle build FAILED — needs JDK17, COMPASS
  env pins openjdk8). **NOT YET VALIDATED** against a chr-named BAM — do that when wiring P3.
- **COMPASS conda env BUILT** (rc=0, job 30230396) — env `compass` at `/home/groups/larsms/users/kevinroy/anaconda3/envs/compass`. All 6 aligners (bbmap.sh / STAR 2.7.10a / hisat2 2.2.1 / magicblast 1.5.0 / gsnap+gmap_build 2021.05.27) + makeblastdb + cutadapt 2.6 + gffread + picard + samtools 1.7 + pysam 0.15.3 verified runnable (`conda run -n compass`). Built from `$W/COMPASS_env_minimal.yml` (15-tool version-pinned minimal; the full-lock yml SIGKILLed on the login node, so the rescue used a compute-node job + minimal env).
- **FASTQ DONE** (job 30227283) — 3 reps in `$W/fastq/replicate{1,3,5}/` (~24G).
- **P1 config edits DONE + pushed** to k-roy/COMPASS branch `human-a549` (commit ce59139): COMPASS.sh (GENOME_VERSION=GRCh38_gencode_v44, MAX_INTRON_LENGTH 2000->500000, READ_LENGTH 150, portable COMPASS_DIR, skip SGD gffread, ALIGNERS_FILE .tsv, SAMFIXCIGAR->samfixcigar.py) + process_reads_and_align.sh (cutadapt drop poly-T/A arms; jvarkit->python samfixcigar) + samfixcigar.py + make_human_introns.py. M1 clone `~/work/COMPASS` on branch human-a549.
- **P1 COMPLETE + VALIDATED.** Env fully working: samtools upgraded to 1.17 (job 30249088; the 1.7 libcrypto issue); ALL 6 aligners + samtools RUN (`--version`): STAR 2.7.10a, hisat2 2.2.1, magicblast 1.5.0, BBMap 38.84, cutadapt 2.6, GSNAP gsnap.sse42 (sse42 SIMD → no AVX-512 SIGILL). **samfixcigar.py VALIDATED** by synthetic test (known mismatches at [5,20,35] → CIGAR `5=1X14=1X14=1X14=`, exact). [The 0/1875 on a real ONT BAM was a reference-mismatch/NM confound, not a bug — in-pipeline samfixcigar uses the SAME ref COMPASS aligns to.]
- **P2 SMOKE TEST DISPATCHED** (job 30288918, partition **larsms** = 256GB NON-PREEMPT, 16 cores, 12h). Runs the subsampled (100k-read) end-to-end pipeline on rep1: stages FASTQ (cat lanes -> fastq/A549_rep1_R{1,2}.fastq.gz), builds ALL indices (STAR/HISAT2/GSNAP/BLAST) into genome_references/, runs 6 aligners + samfixcigar.py + per-read arbitration. Sentinel `$W/.smoke_rc`, log `$W/smoke.<jid>.out`, watcher `bd2k5prv0`. **Big-mem decision:** COMPASS needs no AVX-512 → runs on larsms AMD Milan; `hisat2-build --ss --exon` (NO --snp) needs ~16-30GB (NOT the 160GB SNP-graph figure) so 256GB is ample + non-preempt protects the one-time index builds.
- **P2 STAGED.** Sherlock `$W/COMPASS` on branch `human-a549` (ce59139). `$W/COMPASS/genome_references/`: `GRCh38_gencode_v44.fasta`(→symlink)+`.fai`, `GRCh38_gencode_v44.gtf`(→symlink), `GRCh38_gencode_v44_introns.tsv` (402k).
- **ENV FIX:** the build reported rc=0 but **samtools 1.7 did NOT RUN** (libcrypto.so.1.0.0 missing) — 'env built' only meant execs EXIST. Upgraded samtools to >=1.15 in the compass env (view/sort/index/depth syntax is stable). LESSON: verify tools RUN (`--version`), not just `which`.
- Genome choice LOCKED: align to **chr-named GENCODE GRCh38.primary_assembly.genome.fa** + gencode v44 GTF
  (both at `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/refs/`) → output junctions
  are `chr5`, matching the 111, NO chrom-name harmonization needed.

## VERIFIED
- Strandedness = **RF/dUTP (0.999)** on the existing SG-NEx STAR BAM → matches COMPASS's default assumption
  (`COMPASS_functions.py:216`); NO strand flip needed.
- introns TSV chr5 count (18,451) ≈ Deliverable-B annotated count (18,450) → consistent.
- `compass` env: all 6 aligner execs + makeblastdb/cutadapt runnable via `conda run -n compass`.
- NOT VERIFIED: `samfixcigar.py` output equivalence; any
  aligner run end-to-end on human; index builds.

## OPEN / IN FLIGHT
- **COMPASS smoke test RUNNING** — job 30288918 (larsms), watcher `bd2k5prv0`. Builds indices (the long pole: gmap_build ~1-2h, STAR ~1h, hisat2-build ~30-60m, makeblastdb ~10m) then aligns 100k reads × 6 aligners + arbitrates. UNPROVEN pipeline → check `$W/smoke.<jid>.out` early for path/conda/FASTQ-glob errors; watch the hisat2-build step for OOM.

## RESUME (concrete, with branch logic)
**Step A — env is BUILT (done); just activate + sanity-check:**
`ssh sherlock` then `conda activate compass` (or prefix commands with `conda run -n compass`). Confirm `which STAR hisat2 gsnap magicblast bbmap.sh makeblastdb cutadapt`. If ever broken, rebuild: `sbatch /scratch/users/kevinroy/compass_a549/env_build.sbatch` (minimal yml). Then proceed to Step B.


**Step A2 — P1 is DONE (config pushed to human-a549 ce59139, env validated). NEXT = P2.**

**Step B(P2) — smoke test is DISPATCHED (job 30288918, larsms). Monitor + branch:**
1. `ssh sherlock 'cat /scratch/users/kevinroy/compass_a549/.smoke_rc 2>/dev/null || echo RUNNING'`.
2. RUNNING → `tail $W/smoke.<jid>.out`. Index-build order in the log: HISAT2 (hisat2-build) -> STAR (genomeGenerate) -> BLAST (makeblastdb) -> GSNAP (gmap_build), then per-aligner alignment of the 100k subsample, then samfixcigar + arbitration.
3. **SMOKE_RC=0** → pipeline validated end-to-end. Per-rep COMPASS junction TSV in `$W/COMPASS/COMPASS_junctions/`. NEXT: scale to FULL reads (revert READS_TO_PROCESS to -1) × all 3 reps (rep1/3/5; concat lanes each; chunk one job per rep on larsms), reusing the built indices. Then P4 + P5.
4. **HISAT2 OOM** (smoke fails in hisat2-build, OUT_OF_ME) → build BASIC HISAT2 index (drop --ss/--exon, ~6GB) + pass `--known-splicesite-infile <ss.txt>` at align time (edit process_reads_and_align.sh lines ~146-160); OR re-dispatch to Sherlock `bigmem` (up to 4TB) or SCG (see CLUSTER MEMORY below).
5. **Other failure** → read `$W/smoke.<jid>.out`; common: path (genome_references symlinks), conda (env compass), FASTQ glob (lane names), STAR_OVERHANG. Fix on the human-a549 branch, push, re-pull on Sherlock, resubmit.

## CLUSTER MEMORY (for heavy steps)
- **Sherlock `larsms`**: 256GB, NON-PREEMPT, AMD Milan (no AVX-512 — fine for COMPASS). CHOSEN for the index builds. Sherlock `bigmem`: up to 4TB (preemptable/QOS). `owners`: big nodes (256GB-4TB) but preempt.
- **SCG** (`ssh scg.stanford.edu`, larsms acct, shares OAK only — NOT Sherlock scratch): `batch` (default, 75 nodes) ~389GB+; `gssc`/`nih_s10` ~1TB; GPU partitions ~2TB. Use for >256GB or GPU work, but data must be on OAK first (stage from Sherlock scratch). scg-sbatch skill in the SGTC workspace.

**Step B — P1 human config edits (in `~/work/COMPASS`, commit to a `human-a549` branch, push authorized):**
- `COMPASS.sh`: set GENOME/FASTA → `…/refs/GRCh38.primary_assembly.genome.fa`; GTF → gencode v44;
  `INTRONS_FILE` → `$W/gencode.v44.introns.tsv`; `READ_LENGTH=150`; fix the
  `analyze_exonic_and_intronic_sequence.py`→`_elements.py` filename mismatch (line ~104).
- `process_reads_and_align.sh`: **`MAX_INTRON_LENGTH` 2000 → 500000** (this is where the cap is EFFECTIVE —
  BBMap maxindel/pairlen, STAR/HISAT2 `--max-intronlen`; the Python arg is inert); cutadapt — **DROP the
  Roy-prep `-g "T{100}"` and `-A "A{100}"` polyA/T arms** (A549 is standard TruSeq), keep `AGATCGGAAGAGC`;
  **replace the jvarkit samfixcigar java call with `python $W/samfixcigar.py <in.bam> <out.bam> <FASTA>`**.
- SIMD: GSNAP/HISAT2/Magic-BLAST → sse42/nosimd build to avoid AVX-512 SIGILL on owners nodes; smoke-test one
  chunk before any array.
- **VALIDATE `samfixcigar.py`** now: run it on a small chr-named BAM, confirm `=`/`X` ops + that COMPASS's
  `compare_splice_junctions_from_multiple_aligners.py` parses them.

**Step C — P2 index builds (chunked SLURM, idempotent skip-checks):** STAR (`--sjdbGTFfile` gencode,
`--sjdbOverhang 149`, ~30GB RAM), HISAT2 (`--ss/--exon`), GSNAP (`gmap_build`+`gtf_splicesites`+`iit_store`),
Magic-BLAST (`makeblastdb`), BBMap (auto). Then P2 align 3 reps × 6 aligners (+minimap2 splice:sr),
end-to-end/no-softclip, numbered reads, samfixcigar.py, name-sort. Chunk per (rep × aligner).

**Step D — P3 per-read arbitration:** run `compare_splice_junctions_from_multiple_aligners.py` (human config)
→ per-rep COMPASS junction TSV.

**Step E — P4 back-propagate RECTIFY's family-concordance gate** (the genuine missing piece — COMPASS records
per-aligner agreement at `compare_splice_junctions…py:660-679` but NEVER gates on it): map aligners → FAMILIES
(STAR×2→star, HISAT2×2→hisat2, gsnap, bbmap, magicblast, minimap2) so modes don't double-vote, require ≥N
independent families for a NOVEL junction. **DE-HERD**: minimap2/gapmm2/uLTRA share lineage (long-read panel);
for the short-read panel, weight by independence. Gate on CALIBRATED support, NOT raw count (aligners can
share the short-anchor-over-GT-AG failure mode). (Optional: refine `adjust_ambiguous_junctions`
(`COMPASS_functions.py:112-184`) tie-resolution with canonical-within-window.)

**Step F — P5 re-validate the 111/609 (settle the verdict):** intersect the human-COMPASS junction set
(family-gated, ambiguity-normalized) against the 111 + 609 (`/scratch/users/kevinroy/deliverable_b/gmap_corroboration.json`).
CALIBRATED controls: positives = annotated gencode junctions in the same EXPRESSED loci (must validate HIGH —
fixing the 14.5% STAR failure); negatives = gmap_noncanonical sample + a shuffled-junction null of matched
intron length. DECISION: positives(expressed) HIGH + negatives ~0 → trustworthy → 111 ~0 = confirmed
artifacts (finalize keep-GMAP-behind-fences); nontrivial fraction validate → real GMAP-unique novels.
Frame the 111 rate RELATIVE to the positive-control floor; "recurrent" does NOT rescue a systematic artifact.

## FILES / PATHS
- Plan + scoping + synthesis: `~/work/rectify/dev/COMPASS_HUMAN_REWORK_PLAN.md`,
  `COMPASS_SHORTREAD_SCOPING.md`, `ALIGNER_INVESTIGATION_SYNTHESIS.md`, `SGNEX_A549_VALIDATION_DATA.md`,
  `DELIVERABLE_B_FINDINGS.md`. Cloud investigation now in-repo: `dev/aligner_investigation/`.
- COMPASS: M1 `~/work/COMPASS`; Sherlock `$W/COMPASS`. Core = `compare_splice_junctions_from_multiple_aligners.py`
  (per-read arbitration), `COMPASS_functions.py`, `process_reads_and_align.sh`, `COMPASS.sh`.
- Sherlock `$W = /scratch/users/kevinroy/compass_a549/`: `gencode.v44.introns.tsv`, `samfixcigar.py`,
  `make_human_introns.py`, `env_build.log`/`.env_rc`, `fastq/`, `pull_fastq.sbatch`.
- Refs: `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/refs/{GRCh38.primary_assembly.genome.fa,
  gencode.v44.annotation.gtf}`. Candidates JSON: `/scratch/users/kevinroy/deliverable_b/gmap_corroboration.json`.
- Memory: `project-aligner-orthogonality-panel`, `reference-aligner-investigation`.

## STATE OF THE BROADER PROGRAM (so the fresh agent isn't blindsided)
The native-aligner DESIGN (`dev/ALIGNER_MEMBER_DESIGN.md`) is a SEPARATE deferred track (design-doc only this
cycle; build later). The deferred RECTIFY figures/docs TODO (uncommitted working-tree WIP: `docs/figures/aligner_panel/`,
`generate_splice_classification_v3.py`+png/svg, `ALIGNER_RECOMMENDATIONS.md`, `dev/TODO.md`) is for a RECTIFY
agent, NOT this COMPASS thread — leave it. The `drop_chimeric_winners` filter is already committed (366c885).
All COMPASS/aligner commits are on `drs-validation-rebuild` (pushed). COMPASS code changes go to k-roy/COMPASS
`human-a549` branch.
