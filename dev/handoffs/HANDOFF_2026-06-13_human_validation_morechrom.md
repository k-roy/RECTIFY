# HANDOFF — Human DRS validation read set (A549), multi-chrom expansion in flight

**Date:** 2026-06-13 · **For:** a fresh agent resuming after Sherlock job **29466209** finishes.
**Read first:** `dev/handoffs/STATUS_human_validation_readset.md` (living doc — full background,
decisions, the 4 selector fixes, why multi-chrom). This file is the **resume runbook**.

---

## 1. What was done

- **Human penalty-table arm review** (separate sub-task, DONE, **[uncommitted]**): fixed doc drift
  in `rectify/data/genomes/homo_sapiens/penalty_tables/PENALTY_TABLE.md` + `PROVENANCE.json`
  (junction_overhang/cDNA/qsrev are bundled, were claimed absent); fixed a provenance bug in
  `rectify/core/splice/calibrate_junction_overhang.py` (header printed unclamped `n_concordant`);
  corrected the shipped human `junction_overhang_table.tsv` header `4`→`3`; added
  `tests/test_data_bundling_human.py` (19 tests incl. consumer-load smokes) +
  `tests/test_calibrate_overhang_header.py`. Added 2 TODOs to `dev/TODO.md` (yeast cat7 selector
  latent length-bug; Sherlock `correct_command.py` `--netseq` reconciliation).
- **Validation read set, chr5 (Sherlock)**: synced Sherlock `rectify/core/correct/walkback.py` to
  M1 HEAD `9f613a6` (was ~3 wks stale, missing the cat1 +strand fix; backup
  `walkback.py.bak_pre_head_sync_20260612`); re-aligned A549 chr5 @ `--max-intron 500000`
  (minimap2 only — uLTRA/deSALT are intron-model-independent, reused); full-depth `correct`
  (no `--annotation`); classify→select→IGV-vet with Kevin.
- **4 selector fixes** (in the SCRATCH `classify_candidates.py` / `select_validation.py`, NOT repo):
  cat7 length→support, cat5 same_intron→overlap, cat8 poly-A-run preference (not `polya_len`),
  global ≥90%-identity quality floor. **classify/select refactored chromosome-aware.**
- **Multi-chrom expansion LAUNCHED**: streamed chr1/11/17/19 from the same 4 A549 runs (2.42M
  reads), subsampled ~1M. First align job (29463887) FAILED at deSALT: **deSALT needs a pre-built
  index** (`deSALT index <ref.fa> $W/morechrom/desalt_index`) — minimap2/uLTRA auto-build, deSALT
  does NOT. Added the index build to `morechrom_align.sbatch` (skip-checked). **Resubmitted as
  29466209**; minimap2 (1.33 GB) + uLTRA (1.25 GB) BAMs survived, so this run only builds the
  deSALT index (~10–30 min) + aligns deSALT (~1–2 h).

## 2. What's verified

- chrom-aware classify/select **reproduce chr5 byte-for-byte** (cat5=24, all categories 8/8) — the
  refactor is safe.
- IGV (live, Kevin): cat5 genuine chimeric (overlap-fix works), cat8 real-CPA + poly-A run visible,
  cat4 legit (false N-junction over poly-A). After quality floor, **8/9 categories pick ≥90%
  identity**; cat5 the exception (genuine chimeric inherently low-quality → drove multi-chrom).
- M1 tests: `pytest tests/test_data_bundling_human.py tests/test_data_bundling.py tests/test_calibrate_overhang_header.py` → **40 passed**.
- morechrom job early stages clean: subsample = **999,228 reads**, trim-polya running, walkback md5
  `50ce85c` confirmed.
- **NOT VERIFIED:** morechrom align completion (running); whether multi-chrom fills cat5 +strand
  with ≥4 clean reads (the whole point); the committed artifact (not built); native-RNA004 table.

## 3. Open items (with why)

- **morechrom align (job 29466209) IN FLIGHT** — ETA ~4–8 h from 16:48 PDT 2026-06-13. owners is
  preemptable; per-aligner skip-checks make requeue cheap (resumes at the unfinished aligner).
- **walkback concurrent-session git race**: Sherlock walkback = `9f613a6`'s version (`50ce85c`); M1
  HEAD is now **`5f97af9`** ("docs(walkback): addendum — DRS policy is 100% non-A"). CHECK whether
  `5f97af9` touched walkback LOGIC or only docs: `git show 5f97af9 -- rectify/core/correct/walkback.py`.
  The whole fixture (chr5 + morechrom) is internally consistent on Sherlock's `50ce85c`. If `5f97af9`
  is logic, decide whether to re-sync + rebuild; if docs-only, ignore.
- **The 4 selector fixes live in SCRATCH scripts only** (`~/igv_data/a549_validation/scripts/`
  on M1 + `sgnex_a549/*.py` on Sherlock) — they are analysis scripts, NOT in the rectify repo. Do
  not "re-derive" them. The yeast-cat7-selector TODO in `dev/TODO.md` notes the same length-bug.
- **Committed artifact not built**: must mirror `rectify/data/validation/` (per-aligner BAMs in
  `aligners/`, `corrected_reads.tsv`, `validation_reads_dorado_source.bam` = untrimmed, trim
  metadata, PROVENANCE.json). Decisions locked: **multi-chromosome** (chr1/5/11/17/19); **preserve
  untrimmed dorado_source + trim metadata** (poly-A recoverable, mirrors yeast); **NO patient
  reads** (A549 is public; Sumner patient data is prototyping-only); R9.4 + M1-HEAD-sha in provenance.
- **Heavy uncommitted M1 working tree** — some is THIS session (penalty docs, tests, walkback sync,
  calibrate, TODO), some is OTHER sessions (`correct_command.py` netseq arm, `BUGS_TO_FIX.md`,
  `PLOT_SKILLS.md`, `dev/_*.py` probes). **Stage surgically by explicit path** (CLAUDE.md). Do NOT
  `git add -A`.

## 4. Resume command

**STATE 2026-06-13 18:48 PDT — READY FOR VETTING.** All mechanical pipeline DONE: align (29466209),
correct (29468007), classify+merge+select (29469951), extraction (29470587). **Combined IGV bundle
built and PULLED TO M1: `~/igv_data/a549_validation/combined/`** — `combined_ref.fa` (chr1/5/11/17/19),
`gencode.5chrom.sorted.gtf.gz`, `aligners/valc.{minimap2,uLTRA,deSALT}.bam`, `VETTING_all.tsv`
(72 reads, 8/cat, padded `igv_locus`), `README_VETTING.md` (per-category PASS/REJECT + cat5 decision),
`load_igv.sh`. **The next session just vets — no pipeline work remains.**

**⚠ cat5 +strand finding (resolve during vetting):** even across 5 chroms / ~1.3M reads, +strand
cat5 has only **2 reads ≥90% identity** (RPS15/chr19 91.5%, ERBB2/chr17 90.4%); the other 2 are
87.7% (UQCRQ) and 85.8% (SEC14L1) — near the 87.8% median. Genuine chimeric is fundamentally
rare-AND-noisy (aligner disagreement ⇔ hard-to-align read). −strand cat5 is all clean (90–91%).
Kevin to decide during vetting: accept the 2 median-quality +strand reads, or go asymmetric
(2+/4− = keep only ≥90%), or drop cat5. More chromosomes will NOT fix this (rate is fixed).

**Resume:** check extraction job 29470587 →
`ssh sherlock 'sacct -j 29470587 --format=State,Elapsed,ExitCode'`. If COMPLETED (or
`$W/combined/VETTING_all.tsv` exists): pull the bundle to M1 and start vetting (step 6). If
FAILED/RUNNING: tail `logs/extract_combined.29470587.*`; resubmit `sbatch $W/extract_combined.sbatch`
(idempotent). Earlier stages (align/correct/classify/select) are all DONE — do NOT redo them.

All paths under `W=/scratch/users/kevinroy/rectify_human_validation/sgnex_a549`. Env: `source
/home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh && conda activate rectify`.

1-4. **DONE.** Outputs on disk: `morechrom/correct_fulldepth/morechrom_corrected.tsv`,
   `morechrom_candidates.tsv`, `combined_candidates.tsv` (1.28M), `combined_selection.tsv` (72 reads).
2-4. **Classify morechrom + merge + select** — **job 29469951** (`morechrom_classify_select.sbatch`,
   compute node 64G). ⚠ classify at ~1M reads × 3 BAMs is **memory-heavy — MUST run on a compute
   node, NOT login** (login silently OOM-killed it the first time; chr5's 319k ran on login fine).
   ⚠ **chrom-label bug fixed 2026-06-13**: classify's output write line hardcoded the chrom column to
   `chr5` (the chrom-aware refactor missed it) → morechrom candidates were mislabeled `chr5`. Fixed to
   `reads[rid]['chrom']`. The chr5-only runs were unaffected (chr5 data → `chr5` is correct), so prior
   chr5 vetting stands. **Test-gap:** the chr5 regression test can't catch this — a multi-chrom fixture
   would. (Scratch scripts on M1 `~/igv_data/a549_validation/scripts/` + Sherlock are current.)
   The job classifies morechrom (chrom-aware) → `morechrom_candidates.tsv`, merges with chr5
   (`a549_candidates_fulldepth.tsv`) → `combined_candidates.tsv`, selects (quality floor, FULL gencode
   `sumner_lab/references/gencode.v44.basic.annotation.gtf.gz`) → `combined_selection.tsv`, and prints
   the cat5 picks with identity. Check `morechrom_candidates.tsv` exists; if classify failed, resubmit.
5. **Extract reads × 3 aligners from BOTH sources** (chr5 reads from `$W/alignments_500k/a549_chr5_trimmed.*.bam`;
   new-chrom reads from `$W/morechrom/alignments/morechrom_trimmed.*.bam`) → one combined per-aligner
   val BAM each (combined `@SQ` header = chr1/5/11/17/19). Pull to M1 `~/igv_data/a549_validation/combined/`.
6. **Vet ALL 9 categories with Kevin in IGV — every read, not just cat5.** (Kevin, 2026-06-13:
   "I would like to vet all categories because those weren't even vetted by me yet. These are public
   facing and must be high quality.") cat5 was the only category Kevin eyeballed so far; cat1/2/3/4/
   6/7/8/9 passed only the quality floor + agent checks, NOT his eyes. **Do not commit any category
   he hasn't personally signed off.** Setup (he runs IGV on M1; agent drives port 60151):
   - Build combined ref (chr1,5,11,17,19) from `error_model_gm12878/refs/GRCh38.primary_assembly.genome.fa`
     via `samtools faidx ... chr1 chr5 chr11 chr17 chr19`; pull to M1. Kevin loads it: **Genomes → Load
     Genome from File** (port has NO genome command). Sequence bases require the LOCAL fasta (hosted
     hg38 loads genes but NOT sequence → reads render as `=`).
   - Agent: `curl -s "http://localhost:60151/load?file=<gencode.sorted.gtf.gz>"`, then each BAM, then
     `curl -s "http://localhost:60151/goto?locus=<chrN:start-end>"`. Build a `VETTING_all.tsv` with an
     `igv_locus` (±500 bp) column for all 72 reads; step through **systematically, category by
     category**, present each read's signal + the per-category PASS/REJECT criteria (see the
     `README_VETTING.md` template), record Kevin's verdict per read.
   - For each REJECT: pull the next-best same-category/strand candidate from the merged inventory
     (respecting the quality floor), re-extract, re-present — iterate until Kevin signs off all 8/cat.
   - Name BAMs distinctly per re-pull (e.g. `valc.*`, `valc2.*`) to avoid IGV stale-track caching; or
     have Kevin refresh the session before reloading.
7. **Build the committed artifact ONLY after Kevin signs off all 9 categories.** Mirror §3 layout,
   apply **XV** (category) / **XG** (gene) tags, commit publicly. PROVENANCE records: A549 SG-NEx
   public, R9.4 caveat, chroms 1/5/11/17/19, minimap2 `--max-intron 500000` + uLTRA + deSALT,
   ~1M+319k reads, rectify M1 HEAD sha, the 4 selector-fix lineage, per-read Kevin-vetted flag.

## 5. Files touched

**rectify repo [uncommitted]:**
- `dev/handoffs/STATUS_human_validation_readset.md` [new, uncommitted] — living doc (read first)
- `dev/handoffs/HANDOFF_2026-06-13_human_validation_morechrom.md` [this file]
- `dev/TODO.md` — +2 TODOs (yeast cat7 selector; netseq reconcile)
- `rectify/data/genomes/homo_sapiens/penalty_tables/{PENALTY_TABLE.md,PROVENANCE.json,junction_overhang_table.tsv}` — doc/header fixes
- `rectify/core/splice/calibrate_junction_overhang.py` — header clamp fix
- `tests/test_data_bundling_human.py`, `tests/test_calibrate_overhang_header.py` [new]
- `rectify/core/correct/walkback.py` [ M — this is concurrent-session work, NOT mine; see §3 race]
- Other ` M`/`??` (correct_command.py, BUGS_TO_FIX.md, PLOT_SKILLS.md, dev/_*.py …) = OTHER sessions.

**Scratch (NOT in repo) — the 4 selector fixes + chrom-aware refactor live here:**
- M1: `~/igv_data/a549_validation/scripts/{classify_candidates,select_validation}.py`,
  `~/igv_data/a549_validation/` (IGV bundles, READMEs, load_igv.sh)
- Sherlock: `$W/{classify_candidates,select_validation,post_correct_pipeline}.py/.sh`,
  `$W/morechrom/{stream_chroms.sh,morechrom_align.sbatch,morechrom_correct.sbatch,morechrom_ref.fa,
  gencode.morechrom.gtf,polyasite_morechrom.bed}`

**Sherlock data:** `$W/correct_500k_fulldepth/` (chr5 corrected), `$W/alignments_500k/` (chr5 BAMs),
`$W/morechrom/` (new-chrom: merged.bam 2.42M, alignments/ being written, correct_fulldepth/ pending).
