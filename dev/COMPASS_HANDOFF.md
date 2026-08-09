# COMPASS human-A549 rework — HANDOFF (2026-06-19, jobs submitted)

Self-contained pickup for a FRESH agent continuing from P0. Read the cluster/workspace rules in
`~/.claude/CLAUDE.md`, `~/work/CLAUDE.md`, `~/work/UCLA/Chanfreau_Lab/CLAUDE.md` FIRST (ControlMaster
discipline, owners partition + AVX-512 constraint for rectify-env jobs, chunking, never relay big data
through the M1). Deep context: `dev/COMPASS_HUMAN_REWORK_PLAN.md` (the phased plan),
`dev/COMPASS_SHORTREAD_SCOPING.md` (scoping + STEP-0), `dev/ALIGNER_INVESTIGATION_SYNTHESIS.md` (de-herding +
design insights). COMPASS codebase map is in this session's transcript; key facts reproduced below.

## ✅ PIVOT RESOLVED (2026-06-19) — Option B chosen; ACTIVE PLAN = `dev/RECTIFY_SHORTREAD_MODE_PLAN.md`
The user chose **Option B**: build the short-read multi-aligner junction caller as a **mode inside RECTIFY**
(skip 5'/3' end-correction); COMPASS algorithm survives as RECTIFY's `core/consensus/select.py` module.
Tiebreak = **configurable, default COMPASS published order**. The entire COMPASS-repo build below is now
**SUPERSEDED** (kept only for the human params/indices it produced, which the RECTIFY build REUSES). **Do all
new work per `dev/RECTIFY_SHORTREAD_MODE_PLAN.md`** (phased P0-P5, file:line targets, the paired-end + panel
work items, and the preserved Sherlock data/indices/env). Two RECTIFY infra investigations are DONE (their
maps are folded into that plan). Sherlock jobs cancelled; nothing running. The sections below are historical
COMPASS-build context only.

## ⛔ (HISTORICAL) STRATEGIC PIVOT — build path below was the COMPASS-repo approach, now superseded by Option B
The COMPASS-human build (everything below) is **paused pending an architecture decision**. Trigger: the user asked whether we run properly chunked job arrays (we do NOT — old `process_reads_and_align.sh` runs 7 aligners SEQUENTIALLY in one job, wall=sum, no failure isolation, plus FASTQ renumbering at lines 112-113) and directed that **RECTIFY's proven chunking/parallelization defaults + its "RN sidecar" read-name method (non-negotiable, optimized over months) be PORTED rather than re-solved**. The user is now seriously reconsidering the whole "build inside COMPASS" decision in favor of **adding a short-read MODE to RECTIFY (skipping 5'/3' end-refinement), reusing its infrastructure**, accepting loss of the standalone COMPASS repo (the COMPASS *algorithm* — per-read multi-aligner arbitration — survives as a named module/subcommand inside RECTIFY, staying citable).
**THE FORK:** (A) keep building in COMPASS, import RECTIFY infra; (B) short-read mode in RECTIFY, skip end-refinement, COMPASS-as-module *(user's lean)*; (C) extract shared infra (chunking/RN/parallel) into a common lib both import.
**Decision inputs being gathered (2 background investigations of `/Users/kevinroy/work/rectify`):** (1) what RECTIFY's RN sidecar is + chunking/parallelization defaults (file:line, concrete values); (2) how much of the COMPASS arbitration core is ALREADY in RECTIFY's junction-refiner/consensus-scoring, and how COUPLED RECTIFY's infra is to end-correction (can a short-read mode cleanly skip it / is infra cleanly extractable). **Sherlock jobs CANCELLED** (monolithic Stage A 30322715 + Stage B 30323309); numbered FASTQ + all genome indices PRESERVED, nothing burning.
**RESUME for the pivot:** if the 2 investigations' findings are not in hand, re-run them (prompts above). Then: run an Advisor pass on A-vs-B-vs-C, present the audit with a recommendation, and use AskUserQuestion to let the user pick the architecture BEFORE any further build. Do NOT resume the COMPASS monolithic build below until the user picks A/B/C. If they pick A, the build below + the chunked-array refactor (prep done→7-way parallel aligner array reusing the preserved numbered FASTQ→compare) applies. If B/C, the COMPASS short-read alignment + arbitration gets wired into RECTIFY's chunking/RN framework instead.

## ⭐ MAJOR FIND (2026-06-19) — PRIOR HUMAN COMPASS PIPELINE EXISTS; ADOPT IT (do not keep patching yeast scripts)
An Oak/scratch search found a COMPLETE, SUCCESSFUL prior human COMPASS run (PRPF18 RNA-seq, GRCh38, Jan-Feb 2022, all 7 aligners end-to-end through junctions + DESeq2). **We had been re-deriving it.**
- **Human scripts (Oak, persistent):** `/oak/stanford/groups/larsms/Users/kevinroy/projects/collaborations/COMPASS/main/scripts/` (+ `archived/`): `COMPASS_human.sh` (driver), `process_reads_and_align_human.sh`, `compare_splice_junctions_from_multiple_aligners_human.py`, `COMPASS_functions_human.py`, `add_exonic_intronic_sequence_human.py`, `COMPASS_filtering_human.R`, `compare_individual_aligners_human.R`, `COMPASS_combine_junction_tables_from_multiple_samples.py`. THESE handle human-scale coords/naming — START HERE.
- **Prebuilt GRCh38 indices (scratch, EXPIRE ~2026-07-02 purge):** `/scratch/users/kevinroy/COMPASS_alignments_archive/genome_references/GRCh38/` — bbmap/STAR(75bp)/HISAT2/GSNAP/BLAST + `GRCh38_latest_genomic.fasta/.gtf/.fai/.dict` (**NCBI RefSeq, NC_0000xx names**) + the HISAT2 `_splice_sites.txt` = the human INTRONS_FILE. `log/human/` has the exact resolved aligner commands.
- **CANONICAL human BBMap params (the fix):** `maxindel=200000 pairlen=200000 intronlen=20` (I used 500000 -> 37% spurious deletions). 2022 used BBMap DEFAULTS (ambig=best; NO ambig=random/xstag=us — my instinct there was wrong for matching the tested setup).
- **Human INTRONS_FILE = HISAT2 `hisat2_extract_splice_sites.py gencode.gtf > splice_sites.txt`** (NOT my make_human_introns.py GTF-derived TSV; different coord convention).
- KEY FORK (needs PI): the prebuilt indices are NCBI-RefSeq (NC_ names) + STAR@75bp. Our 111 candidates are GENCODE chr5. EITHER reuse NCBI indices (save hours; convert NC_000005->chr5 for the 111 intersection; STAR still rebuild @150bp) OR rebuild on GENCODE chr-named. Indices expire 2026-07-02 -> decide soon.

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

## ⭐ DECISION (2026-06-19) — chr5-SCOPE THE COMPARE STEP (advisor-vetted)
The full-genome compare_splice OOM (64G, rc137) was **NOT a read-volume problem** — it OOM'd even on a 100k-read subsample. Root cause (confirmed by code read + isolated startup test): `get_ambiguous_junctions_in_annotated_introns` (COMPASS_functions.py:35) runs an UNBOUNDED ambiguity walk (`get_ambiguous_junctions`, lines 15-28) over all **402k** introns; near N-blocks/homopolymers `'N'=='N'` loops and accumulates millions of tuples → 64G. The deliverable (the 111 GMAP novels) is **ALL chr5**, so we scope the COMPARE to chr5: chr5 introns TSV (18,451), chr5 FASTA, chr5-filtered BAMs. **Alignment stays GENOME-WIDE** (reads must map to true origin; we filter BAMs to chr5 only for the comparison). **VERIFIED** by isolated startup test (job 30321680, COMPLETED 12s): chr5 startup = **0.06 GB peak, 1.6s, rc=0**, `start` dtype int64 (no TypeError), ambiguous map bounded at 8,882 entries. The earlier str+int TypeError was a STALE full-TSV log; the current chr5 TSV parses clean.
*To redirect to genome-wide later:* either bigmem (4TB) OR add an N-guard / idx-cap to `get_ambiguous_junctions` (Advisor first — it's core logic). For the 111 question, chr5 is sufficient and correct.

## BBMap FIX — DONE
Non-negotiable workhorse restored. (1) Canonical 2022 params applied: `MAX_INTRON_LENGTH=200000` in COMPASS.sh (commit **6e346d5**, pulled on Sherlock) → bbmap `maxindel=200000 pairlen=200000 intronlen=20`, fixing the 37% spurious 500kb deletions. (2) bbmap re-added to `sample_aligner_info.tsv` (7-config panel). (3) The `chr1 1`→`chr1` @SQ inconsistency: genome FASTA cleaned (`sed '/^>/ s/ .*//'`) into a REAL file `$GR/GRCh38_gencode_v44.fasta` (was a symlink into the rectify workspace — do NOT edit that target), re-faidx'd; stale `$GR/bbmap/ref` deleted so bbmap rebuilds the index from the clean FASTA inside Stage A. samfixcigar.py is also hardened (split()[0] + per-read try/except, 0668db5) as belt-and-suspenders.

## ⭐ STAGE B CORRECTED (advisor, 2026-06-19) — Option A "normal-mode merge", NOT independent BAM filter
The FIRST Stage B (independent per-aligner chr5 BAM filter) was **WRONG**: filtering each aligner's BAM to chr5 separately gives the 7 BAMs DIFFERENT read sets, which (1) can desync compare_splice's per-read merge, and worse (2) HIDES a read's better off-chr5 placement → manufactures FALSE chr5 junctions (the exact artifact class we hunt). **Corrected approach (Option A):** run compare in COMPASS's NORMAL consistent-read-set mode — full genome-wide name_sorted BAMs (identical numbered-read sets across all 7 aligners), **chr5 introns** (this alone fixes the 64G OOM — it was the 402k-intron startup walk), **FULL clean FASTA** (on-demand fetch). Then filter the OUTPUT junction TSV to chr5. The observed-junction walk in the read loop is SAFE genome-wide because real reads never map into N-gaps (only the annotation hit that). Cost: the python read loop processes ALL genome-wide reads (slow — see timing test 30322625). *Future optimization (handoff TODO):* subset all BAMs to the UNION of read-names any aligner maps to chr5 (keeps read sets consistent, cuts reads ~20×) for a much faster compare — defer until the simple Option A is proven.

## LESSON (2026-06-19): `set -u` breaks conda activate
Do NOT put `set -u`/`set -uo pipefail` in a COMPASS sbatch: conda's `openjdk_activate.sh` references an unbound `target_platform` → job FAILS in 13s. Use `set -o pipefail` only, AFTER `conda activate`. (Stage A v1 job 30321881 died this way; fixed in v2.)

## OPEN / IN FLIGHT (2026-06-19) on Sherlock larsms (non-preempt)
- **Stage A = job 30322715** `cmp_rep1_align` (32c/128G/2day): full genome-wide alignment of rep1 (all 7 aligner configs) via `process_reads_and_align.sh` (READS_TO_PROCESS=-1) → `separate_alignments/<aligner>/A549_rep1_name_sorted.bam`. bbmap index rebuilds from clean FASTA here. Sentinel `$W/.rep1_align_rc` (`A:<rc>` or `A_MISSING:<aligners>`). Long pole = GSNAP/MagicBLAST. (v1 30321881 FAILED on `set -u`; this is v2.)
- **Stage B = job 30323309** `cmp_rep1_chr5` (16c/96G/**2day**, `--dependency=afterok:30322715`): the CORRECTED Option-A script `$W/run_rep1_chr5_compare.sbatch` — compare on full genome-wide BAMs (sample `A549_rep1`, full panel incl. fresh clean bbmap), chr5 introns + full FASTA, then filter output to chr5. Sentinel `$W/.rep1_chr5_rc` (`B:<rc>`). (First version 30321976, independent-filter, was CANCELLED as scientifically wrong.)
- **Timing/correctness test = job 30322625 PASSED** (`cmp_subB`, Option-A merge on 100k-subsampled BAMs, 6 aligners): **COMPARE_RC=0**, 4:24 wall, **357 MB** RSS, ~**19k reads/min** steady-state → full rep1 (~42M pairs) ≈ **37h** (under the 2-day limit, margin ~11h). Produced **30,093 junctions (1,476 chr5)**; positive control FIRES (3 chr5 junctions tagged `annotated_junction=True` — low only because 100k genome-wide → sparse chr5). Merge runs clean on consistent read sets = Option A validated.
- **OUTPUT FILE + COLUMNS (verified):** compare writes `$C/COMPASS_junctions/<sample>_COMPASS_splice_junctions.tsv` (NOT `_with_sequence_info.tsv` — that's a later analyze stage). Columns: 1=sample_name, **3=chrom**, 4=adj_start, 5=adj_stop, 6=five_SS, 7=three_SS, 8=RNA_strand, **10=annotated_junction(True/False)**, 16=COMPASS_counts, then per-aligner counts. Stage B filters chr5 on `$3` and counts annotated on `$10`.
- Scripts (on scratch, NOT in repo): `$W/run_rep1_align.sbatch`, `$W/run_rep1_chr5_compare.sbatch`, `$W/cmp_subB_test.sbatch`, `$W/test_compass_startup.py`. reps 3 & 5 fastqs NOT staged (only `$W/COMPASS/fastq/A549_rep1_R{1,2}.fastq.gz`).

## SUCCESS CRITERION (self-evaluable)
Stage B emits `$W/COMPASS/COMPASS_junctions/A549_rep1_chr5_COMPASS_splice_junctions.tsv`, **non-empty**, and the Stage B log prints `annotated_chr5_count` > 0 (chr5 junctions correctly tagged `annotated_junction=True` = positive control PASSES). At full rep1 depth expect FAR more than the test's 3 (sparse). "Job exited 0" alone is NOT sufficient.

## RESUME (concrete, with branch logic)
Pre-reqs ALL DONE (env, indices, config, BBMap fix, chr5 inputs, Option-A validated). BOTH jobs queued: Stage A 30322715 (RUNNING), Stage B 30323309 (PENDING afterok). Just monitor:
1. `ssh sherlock 'cat /scratch/users/kevinroy/compass_a549/.rep1_align_rc 2>/dev/null; echo ---; cat /scratch/users/kevinroy/compass_a549/.rep1_chr5_rc 2>/dev/null; echo ---; squeue -u kevinroy -o "%.10i %.14j %.9T %.10M %R"'`
2. **`.rep1_align_rc` absent** → Stage A still RUNNING/PENDING. `tail $W/logs/cmp_rep1_align_<jid>.out` and `tail $W/COMPASS/log/A549_rep1_process_reads_and_align.log` (per-aligner progress). If Stage A hits the 2-day TIMEOUT: completed-aligner BAMs persist; the bottleneck is GSNAP/MagicBLAST — the per-aligner resolved commands are in `$W/COMPASS/log/A549_rep1_process_reads_and_align.log` (and the 2022 `…/COMPASS/main/log/human/`); rerun the missing aligner(s) by hand from those, then submit Stage B manually (afterok not needed).
3. **`.rep1_align_rc` = `A:0`** and **`.rep1_chr5_rc` = `B:0`** → check the Stage B log `$W/logs/cmp_rep1_chr5_<jid>.out` for `JUNCTION_TSV_ROWS` and `annotated_tag_count`. If TSV non-empty AND annotated>0 → SUCCESS, go to P5 (the 111 intersection below).
4. **`.rep1_align_rc` = `A_MISSING:...`** → that aligner's BAM is empty; read its section in the process log. bbmap-specific empty BAM = re-check the clean-FASTA/index rebuild.
5. **`.rep1_chr5_rc` = `B:<nonzero>`** → Stage B failed AFTER alignment (BAMs are safe). Read `$W/COMPASS/log/A549_rep1_chr5_compare.log`. Re-run ONLY Stage B: `ssh sherlock 'cd $W && sbatch run_rep1_chr5_compare.sbatch'` (no dependency needed; reads the persisted name_sorted BAMs). **If the crash is bbmap-header-related** (`KeyError`/`fetch` on a `chr5 5`-style name — bbmap is the ONLY aligner not exercised by the test 30322625, which excluded it): first confirm the FRESH bbmap BAM `@SQ` is clean (`samtools view -H $C/separate_alignments/bbmap/A549_rep1_name_sorted.bam | grep '^@SQ' | head` should show `SN:chr5`, NOT `SN:chr5 5`); if dirty, the clean-FASTA/index rebuild didn't take — fix that. To UNBLOCK the run meanwhile, edit `run_rep1_chr5_compare.sbatch` to use `$W/aligners_no_bbmap.tsv` (6-aligner) and resubmit, then debug bbmap separately (non-negotiable, so don't drop it permanently).

## P5 — re-validate the 111 (the actual deliverable), once Stage B succeeds
Intersect the chr5 COMPASS junction TSV against the 111 GMAP-only recurrent novels.
- Candidates: `~/work/rectify/dev/gmap_only_recurrent_novels_chr5.tsv` (also `/scratch/users/kevinroy/deliverable_b/gmap_corroboration.json` for the 111+609).
- Junctions: `$W/COMPASS/COMPASS_junctions/A549_rep1_chr5_COMPASS_splice_junctions.tsv` (chr5, ambiguity-normalized by COMPASS). Coord columns (VERIFIED): chrom=`chrom`(col3), donor/acceptor=`adj_start`(col4)/`adj_stop`(col5), strand=`RNA_strand`(col8), `annotated_junction`(col10), per-aligner support in `*_counts` cols (16+), `COMPASS_counts`(col16). The 111's recurrence threshold (≥5 reads) maps to filtering on the per-aligner/COMPASS count columns.
- Intersection: inner-merge the 111 (chr5 start/stop) against the COMPASS chr5 junctions on (`adj_start`,`adj_stop`) with ±a few bp for residual ambiguity (COMPASS already collapses ambiguous junctions to a canonical rep). Report: how many of the 111 appear, with what multi-aligner support and annotated/canonical flags.
- CALIBRATED controls (per the original P5 plan, Step F below): positives = annotated gencode chr5 junctions in the SAME expressed loci (must validate HIGH — this is the cure for the 14.5% STAR floor); negatives ~0. Then: 111 ∩ COMPASS ≈ 0 with positives HIGH → confirmed artifacts; nontrivial fraction → real GMAP-unique novels.

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
