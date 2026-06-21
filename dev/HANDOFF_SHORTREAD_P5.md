# HANDOFF — RECTIFY short-read P5 RUN on Sherlock (2026-06-19)

Picks up `dev/HANDOFF_SHORTREAD.md` (P1–P4 code done, committed, locally tested; pre-flight green).
This doc tracks the **P5 cluster run** of the paired short-read COMPASS pipeline on A549 to adjudicate
the **111 GMAP-only recurrent novels** (`dev/gmap_only_recurrent_novels_chr5.tsv`).

## Key correction to the original handoff's assumption
The original handoff said run rectify **in the compass env**. That is WRONG: the `compass` conda env is
**Python 3.7.12 / pysam 0.15.3**, and rectify requires **py≥3.8** (it dies at import — `config.py:33` uses
py3.9 generic subscription). Resolution:
- **Run rectify in the `rectify` conda env** (`/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify`,
  **py3.9.23 / pysam 0.23.3**) via `python -m rectify` (NOT a pip console-script — we do NOT `pip install`,
  to avoid mutating the shared env; `cd $RECTIFY_SRC` puts our branch on sys.path).
- The 4 specialized COMPASS aligners (**STAR, hisat2, magicblast, gsnap**) live ONLY in the `compass` env.
  They are exposed to the rectify env via **symlinks in `$HOME/.rectify/bin/`** (already on the generated
  script's PATH). Verified all 4 load via symlink with compass NOT on PATH (RPATH/$ORIGIN resolves). bbmap.sh,
  samtools, java come from the rectify env (fine).

## Done
- Deployed branch `drs-validation-rebuild` (local HEAD `b4d4b48`, incl. all P1–P5 commits through `2e0388c`)
  to Sherlock: `rsync --max-size=2M` (drops 2GB bundled yeast data) → `/scratch/users/kevinroy/compass_a549/rectify_src`.
- Bumped two generator constants in `rectify/core/commands/split_command.py` (local + redeployed):
  `SR_ALIGN_CORRECT_MEM_GB 32→64` (STAR loads a **29G** human index; sequential panel → peak ≈ one STAR;
  32G OOMs) and `SR_ALIGN_CORRECT_TIME '1:30:00'→'6:00:00'` (per-chunk 7-aligner index loads dominate;
  larsms is non-preempt). **NOT yet committed in git.**
- Generated + READ the smoke scripts (caught nothing wrong — the ANNOT path that looked off was just the
  COMPASS gtf symlink `.resolve()`-d to its real target `rectify_human_validation/.../gencode.v44.annotation.gtf`;
  identical file, 3,424,194 lines).

## Verified
- `rectify split --help` (deployed code) shows `-2/--read2 --read-length -n --short-read --generate-slurm
  --python-path --rectify-src --slurm-partition --slurm-account`.
- All 4 aligner symlinks load: STAR 2.7.10a, hisat2 2.2.1, magicblast 1.5.0, gsnap (runs; SSE4.2 note only).
- `for aligner in aligners:` in `multi_aligner.py:2802` is **sequential** → mem sizing = one STAR.
- Subsampled smoke input = **100k pairs** each (`A549_rep1_subsampled_{1,2}.fastq.gz`).

## CONSENSUS HANG — root-caused and FIXED (the big one; uncommitted)
Smoke attempt 2 ran the panel but then **hung ~2h in consensus** on 25k pairs. faulthandler traceback:
`select_best_alignment → score_alignment → _rescue_5prime_softclip → edit_distance`. Root cause
(`select.py:71-75`): the per-read 5' soft-clip rescue pool was built by scanning **every annotated
junction on the read's whole chromosome** (~tens of thousands) — O(reads × junctions_on_chrom), which is
fine for chr5-only long-read (deliverable_b) but hangs genome-wide at short-read scale (~50k reads). The
rescue only ever inspects junctions within `search_window_bp=300` of the read 5' boundary, so the chrom
scan was pure waste.
**Fix** (`select.py`): a memoized per-chrom **position index** (sorted start/end arrays + bisect) over
`annotated_junctions`; the candidate pool is now a locus window `[read_lo-350, read_hi+350]` query — a
**window-superset of the rescue's ≤300bp reach ⇒ byte-identical scoring**, O(log n) per read. The
`annotated_junctions` SET is kept intact for the `_n_annotated` tiebreak membership test. Helps the
long-read path identically.
**Verified in tight loop** (no SLURM burn): `repro2_consensus.py` on the 5 smoke BAMs went from >2h-hang →
**DONE 50.6s**, output BAM 64,586 records (= merge groups), **28,009 spliced N-op reads present**.

## gsnap + magicblast FAIL (panel currently 5/7) — env/tool issues, NOT code
- **gsnap**: compass env has GMAP workers (`gmap.sse42`,`gmap.nosimd`) but **NO gsnap worker binaries**
  (`gsnap.sse42`/`.nosimd` absent) — the `gsnap` dispatcher prints "does not exist" for every ISA and exits
  in 0.1s. The GSNAP half of the gmap/gsnap package was never installed. Fix = install complete gmap/gsnap
  (e.g. `conda install -n <env> -c bioconda gmap`, then re-point the `~/.rectify/bin/gsnap` symlink).
- **magicblast**: `.bam` has a valid `@HD` header but **0 records**; `samtools view failed ...magicblast.sam
  (exit 1)`. magicblast 1.5.0 (2019) emitted a SAM the rectify-env samtools rejected. Fix = inspect raw SAM /
  samtools-version compat / params in `run_magicblast`.
- 5 working = bbmap + STAR_default + STAR_noncanonical + HISAT2_default + HISAT2_noncanonical (the sensitive
  splice-aware workhorses). **USER CHOSE: fix env for full 7/7 before the run.**

### gsnap FIX (done, validated)
gmap/gsnap **2021-05-27** workers (`gsnap.sse42`,`gsnap.nosimd`,`gsnapl.*`) DO exist in compass/bin — the
dispatcher just looks for `gsnap.<isa>` **next to argv[0]** (= the `~/.rectify/bin` symlink dir), where only
`gsnap` itself was symlinked. Fix = symlinked the 4 workers into `~/.rectify/bin` too. `gsnap --version` via
that PATH now resolves `gsnap.sse42` cleanly (only a harmless `avx2 does not exist` note). **No env mutation.**
The GSNAP index was built with this same gmap 2021-05-27 → compatible.

### magicblast FIX (done, validated on subset; code change)
magicblast 1.5.0 under `-no_query_id_trim` keeps the whole multi-token FASTQ header (rectify-split injects
`RN:i:N` + the Casava comment) and spills those tokens into SAM cols 2-3, shifting the mandatory fields →
FLAG="RN:i:0" → samtools view exit 1 → 0 records. Dropping `-no_query_id_trim` makes magicblast fail
differently (rc=1). **Fix** (`multi_aligner.py::run_magicblast._ensure_plain`): rewrite every input-FASTQ
header to its bare QNAME (first whitespace token) before magicblast; RN is re-applied afterward from the
qname→RN map in `_finalize_short_read_bam` (keyed on bare qname). Verified on a 1k-pair subset: SAM columns
correct, `samtools view` rc=0, 2273 records. **Uncommitted.**

## ALIGNER-VERSION SENSITIVITY EXPERIMENT (user-requested 2026-06-20)
Goal: confidence aligner version plays a minimal role. Head-to-head on the SAME reads (100k-pair
A549_rep1 subsample), identical rectify pipeline, differing ONLY in aligner binaries/indices.
- Installed (COMPASS-pinned): STAR 2.7.10a, HISAT2 2.2.1, magicblast 1.5.0, gsnap 2021-05-27, bbmap 38.18.
  Latest: STAR 2.7.11b, HISAT2 2.2.1 (==, project dormant), magicblast 1.7.2, gmap/gsnap 2025.07.31, bbmap 39.79.
- Env `compass_latest` BUILT (via **micromamba** — classic conda solver OOMs on the login node; conda has
  no libmamba solver plugin and Sherlock curl is flaky, so micromamba binary was fetched on the M1 and
  rsync'd to `~/bin/micromamba`). Installed: STAR 2.7.11b, HISAT2 2.2.2, magicblast 1.6.0 (bioconda's
  newest is 1.7.0; solver picked 1.6.0), gsnap 2024-11-20, bbmap 39.26, samtools 1.21. Path:
  `/home/groups/larsms/users/kevinroy/anaconda3/envs/compass_latest`.
- **Invocation**: env has its own python 3.12, so run rectify via the ABSOLUTE rectify-env python with
  `PATH=compass_latest/bin:$PATH` (aligners + their wrapper-siblings — gsnap workers, hisat2-align-s,
  bbmap jar — all resolve from that one dir; no symlinks needed). Testing whether STAR 2.7.11b / gsnap
  2024 read the existing indices; rebuild only incompatible ones into `genome_references_latest/`.
- Baseline (COMPASS-pinned panel on the 100k subsample) → `$W/ver_cmp/compass/cmpver.rectified.bam`.
- Latest panel → `$W/ver_cmp/latest/latestver.rectified.bam`.
- Compare: `dev/compass_version_compare.py` (chr5 Jaccard, depth Spearman, novel-junction agreement,
  111∩each). High concordance ⇒ version minimal role. (Separate from the main full run job 30432422.)
- **Index compatibility tested**: STAR 2.7.11b READS the existing 2.7.4a index (reuse); HISAT2 2.2.2 reads
  2.2.1 index (reuse); BLAST db reused; bbmap builds on-the-fly. **gsnap 2024 CANNOT read the 2021 gmap
  index (rc=9)** → rebuilt with latest `gmap_build` (COMPASS recipe: gmap_build + `gtf_splicesites|awk
  '$4>9'|iit_store`). `genome_references_latest/` = symlinks to reusable + fresh GSNAP/.
- **Chained jobs**: `gsnap_latest_build` **30443075** (~1-2h, sentinel `$W/.gsnap_latest_build_rc`) →
  `latest_panel_cmp` **30443095** (`afterok` dep; runs latest panel on the 100k subsample then the
  comparison, writes `$W/ver_cmp/version_compare.json` + sentinel `$W/.latest_cmp_rc`). COMPASS-pinned
  baseline already done: `$W/ver_cmp/compass/cmpver.rectified.bam` (all 7 aligners, 311,614 reads).
- **RESUME (experiment)**: `cat $W/.latest_cmp_rc` → if `0`, read `$W/ver_cmp/version_compare.json`
  (jaccard_all, novel_jaccard, depth_spearman, 111_in_compass_pinned vs 111_in_latest, verdict). If a job
  failed: `sacct --name=gsnap_latest_build,latest_panel_cmp`; gsnap build logs `$W/logs/gsnap_latest_build_*`.
  Versions: STAR 2.7.10a→2.7.11b, HISAT2 2.2.1→2.2.2, magicblast 1.5.0→1.6.0, gsnap 2021-05-27→2024-11-20,
  bbmap 38.18→39.26.

## Test suite (user /goal: all tests pass) — GREEN after fixes
Full `pytest tests/` in the `rectify` env (py3.9), bundled data deployed (54M, indexes excluded):
- My 4 changes cause **zero regressions** — the 9 initial "failures" were all missing-bundled-data
  artifacts of the `--max-size=2M` deploy (yeast genome/gff); they pass once `rectify/data` is present.
- **Added regression tests**: `tests/test_consensus_locus_index.py` (window index correctness +
  far-junction scoring invariance), `tests/test_magicblast_header_strip.py` (bare-qname strip + gz +
  `check_aligner_available(None)` no-raise). Also hardened `check_aligner_available(None)→False` and
  refactored the magicblast strip into module-level `_write_bare_qname_fastq` (testable).
- **One pre-existing failure resolved, not mine**: `test_bam_parallel_state::…_deterministic` is a golden-
  hash test whose golden (recorded 2026-05-22 @55089f7) went stale when later 3'SS-rescue commits
  (bd20f9e/961c844/cf5ebb9/0c1773b) legitimately changed `process_bam_file_parallel` output. Proven not
  mine: baseline `select.py` yields the IDENTICAL observed hash; output is deterministic (stable across
  runs). Re-recorded the golden to `4f326833…` with documented rationale. (It was masked earlier because
  it's skipped when bundled data is absent.)

## ALIGNER-VERSION EXPERIMENT — RESULT (2026-06-20): version plays a NEGLIGIBLE role
Head-to-head on the same 100k-pair A549 subsample, identical rectify pipeline.
- **gsnap caveat**: gsnap 2024-11-20 removed `--ambig-splice-noclip` (rectify's shared wrapper passes it for
  gsnap 2021), so latest gsnap won't drop into the COMPASS flags unchanged → gsnap excluded from BOTH sides
  for a clean matched comparison (its version effect is UNtested, not shown-minimal). The other 4 distinct
  tools all upgraded cleanly: STAR 2.7.10a→2.7.11b, HISAT2 2.2.1→2.2.2, magicblast 1.5.0→1.6.0, bbmap 38.18→39.26.
- **CLEAN 6v6** (`$W/ver_cmp/version_compare_6v6.json`, bbmap+STAR×2+HISAT2×2+magicblast, pinned vs latest):
  - chr5 junctions 1661 vs 1662; **jaccard_all 0.9982**, **depth Spearman 0.9998**.
  - novel-junction Jaccard **0.9954 at ≥1 read**, **1.0 at ≥2/≥3/≥5/≥10 reads** (the single ≥1-read
    discordant junction is one 1-read call — detection-boundary noise, not algorithmic).
  - **111∩panel = 0 for BOTH** versions → the artifact conclusion is version-invariant.
- (As-is 7-pinned vs 6-latest, gsnap-confounded: jaccard 0.969, depth 0.971, novel 0.924 — conservative floor.)
- **Conclusion**: for STAR/HISAT2/magicblast/bbmap, latest vs COMPASS-pinned produce IDENTICAL well-supported
  junctions; aligner version does not change the 111 result.
- **IMPORTANT setup nuance**: `genome_references_latest/` fasta was a SYMLINK; rectify's
  `_compass_index_paths` does `Path(genome).resolve().parent`, and `.resolve()` followed the symlink back to
  the OLD dir — so the 6v6 "latest" panel actually used the OLD indices with the LATEST binaries (a clean
  binary-version isolation, still valid). Fixed by HARDLINKING the fasta+fai into `genome_references_latest/`
  so ref_dir resolves there (needed for gsnap's rebuilt index to be used).
- **gsnap 2024 — NOT a drop-in** (committed flag fix `_gsnap_supports_ambig_noclip`): (1) removed
  `--ambig-splice-noclip` (worked around); (2) requires a freshly-built index (rebuilt); (3) **~100x slower**
  — its new `localdb` mode allocs 12.4G and runs **13.6 queries/sec** (vs gsnap 2021 fast), so the full-100k
  rectify run hit a wall and exited nonzero at 512s. Standalone on 1k pairs it returns rc=0 / correct output,
  so the version is fine, just impractically slow as configured. Direct gsnap-2021-vs-2024 junction
  comparison on a 12k-pair subset: ATTEMPTED but **gsnap 2024 SIGSEGV**s on read
  `K00151:...:1112:3752:43128` (`Access_emergency_cleanup`), truncating output (12594 vs gsnap2021's 25828
  records; 55 vs 238 chr5 junctions) → comparison guardrail correctly refused. So gsnap 2024-11-20 is
  **broken on this A549 data** (crash bug), making its version effect untestable — which is itself the
  strongest reproducibility argument for KEEPING the COMPASS-pinned gsnap 2021-05-27.
- **FINAL ANSWER (user's question)**: aligner version plays a NEGLIGIBLE role for STAR/HISAT2/magicblast/bbmap
  (6v6: novel-Jaccard 1.0 at ≥2 reads, 111∩=0 on both). gsnap-latest is not a drop-in and crashes here, so
  pinning gsnap (as COMPASS does) is the right call; its junctions are corroborated by the other 6 aligners
  in consensus anyway. Net: the 111-artifact conclusion is robust to aligner version.

## Tool maintenance check (user-requested 2026-06-19) — both ACTIVELY maintained; ours are OLD
- **GMAP/GSNAP**: latest bioconda **2025.07.31** (Sep 2025); releases throughout 2024–2025; maintained by
  Thomas Wu / Genentech. **We run 2021-05-27** (~4 yr old). Functional (index built with same version).
- **Magic-BLAST**: GitHub commits into **Apr 2025**; latest tagged release **1.7.2 (Apr 2023)**. **We run
  1.5.0 (Aug 2019)** — ancient; its SAM-column mangling (now worked around) is plausibly fixed in 1.7.x.
- **Recommendation**: for the production/publishable run, consider upgrading both (esp. magicblast → 1.7.2)
  in a fresh env, then drop the magicblast header-strip workaround if 1.7.x handles multi-token headers.
  Not blocking — current fixes make the 2019/2021 versions produce correct 7/7 output.

## 7/7 END-TO-END VALIDATED (chunk-0, 2026-06-19 ~22:25)
Full `rectify align --short-read --read2 --aligners all` on smoke chunk 0 (25k pairs), 7m31s total:
- **All 7 aligners contribute**: By-aligner winners `{STAR_noncanonical 16822, STAR_default 1241,
  HISAT2_noncanonical 6067, bbmap 5206, gsnap 6583, magicblast 33459, HISAT2_default 5}` — gsnap+magicblast
  fixed and winning.
- **Consensus 187s for 69,383 (RN,mate) groups, no hang** (locus-index fix holds at full panel). 5' rescue
  intact (3,456 rescued, edit=0 matches in debug log).
- **Both mates survive**: records-per-RN dominated by 2 (10,223 RN); only 110 RN (0.5%) single-record →
  adversarial mate-drop bug absent. Rectified BAM 69,383 records (multimapper/secondary records inflate
  >2/RN; fine for junction extraction).
Output: `/scratch/users/kevinroy/compass_a549/chunk0_7aligner_test/chunk0test.rectified.bam`.

## Bug fixed before smoke could pass (commit-worthy, uncommitted)
First smoke (`30367486`) fast-failed: `align_command._run_one_aligner` sets `exec_path=None` for the
COMPASS panel (binary resolved inside the wrapper) but still called `check_aligner_available(None)` →
`shutil.which(None)` → TypeError. **Fix** (`align_command.py:475`): guard the check with
`exec_path is not None and not check_aligner_available(exec_path)`. This dispatch path had never executed
(handoff said so). Redeployed; COMPASS dispatch branches (583–629) otherwise present and correct.

## STATUS 2026-06-19 ~22:40 — pipeline 7/7 validated, tests green, NOT yet launched full run
Both smoke arrays (30367486, 30367768) were cancelled — the first hit the consensus hang; the second was
superseded by the direct chunk-0 7/7 validation (see below) once all bugs were fixed. **No cluster job is
currently running.** The pipeline is validated end-to-end at 7/7 on a real chunk; the FULL ~500-chunk A549
run has NOT been launched yet (awaiting go-ahead). The deployed code (`…/compass_a549/rectify_src`) carries
all fixes; `~/.rectify/bin` has all aligner symlinks incl. gsnap workers.

## FULL RUN LAUNCHED 2026-06-20 — split+chain job `30432422`
(First attempt 30431657 fast-failed on the `set -u` + conda-activate trap — conda's java_home.sh
references unbound JAVA_HOME; fixed sbatch to `set -o pipefail` only. Resubmitted as 30432422.)
Submitted `/scratch/users/kevinroy/compass_a549/cmp_sr_full_split.sbatch`. It:
1. `rectify split` the full 42M-pair A549_rep1 (R1/R2) into **500 chunks** → `$W/rectify_sr_full/`
   (~4h at observed throughput; writes sentinel `$W/.sr_full_split_rc` with the rc).
2. On rc==0, chains `bash $W/rectify_sr_full/submit_pipeline.sh` → submits the **500-task array**
   (`A549_rep1_sr`, 64G/6h/task, idempotent `.consensus.bam` skip) + `afterok` final merge.

### RESUME — concrete branch logic for the full run
SSH `sherlock` open; never tear down ControlMaster; retry transient sshd serially. Env: `rectify` conda
env + `export PATH=$PATH:$HOME/.rectify/bin`. `$W=/scratch/users/kevinroy/compass_a549`.
```
ssh sherlock "squeue -u kevinroy -o '%.14i %.16j %.8T %.10M %R'; echo ---; \
  cat $W/.sr_full_split_rc 2>/dev/null; echo '--- chunks done ---'; \
  ls $W/rectify_sr_full/chunk_outputs/*.consensus.bam 2>/dev/null | wc -l; \
  ls $W/rectify_sr_full/final/ 2>/dev/null"
```
- **split job 30432422 RUNNING** → wait (~4h).
- **`.sr_full_split_rc` absent & job gone** → split died; check `$W/logs/cmp_sr_split_*.{out,err}`.
- **`.sr_full_split_rc` == 0** → array was chain-submitted; find it: `squeue`/`sacct --name=A549_rep1_sr`.
  - array tasks idempotent (`.consensus.bam` skip + atomic copy) → safe to requeue / re-`sbatch
    rectify_sr_full/run_array_short_read.sh` if some failed.
  - **all 500 `.consensus.bam` present + merge COMPLETED** → final merged BAM in `rectify_sr_full/final/`.
    Proceed to adjudication (below).
- **`.sr_full_split_rc` != 0** → split failed; do NOT expect an array. Read split err.

### ADJUDICATION (P5 step 4) — tool now BUILT: `dev/compass_shortread_adjudicate_111.py`
Run AFTER the merged BAM exists (see that script's header for the exact command). Reports the 3 numbers:
positive control (annotated chr5 junctions HIGH), negative (~0), and `111 ∩ COMPASS`. Near-zero
intersection ⇒ the 111 are artifacts — ONLY valid if the positive control passed.

## (superseded) Resume — concrete branch logic
SSH ControlMaster `sherlock` is open; never tear it down; retry transient sshd errors serially.
Check the smoke:
```
ssh sherlock "sacct -j 30367486 -X -o JobID,State,Elapsed,MaxRSS,Start,End; \
  ls -la /scratch/users/kevinroy/compass_a549/rectify_sr_smoke/chunk_outputs/ 2>&1"
```
- **If all 4 tasks COMPLETED** → run the acceptance checks (per-chunk log shows all 7 aligners succeeding;
  merged BAM has **two records per RN** — the adversarial-review mate-drop fix; N-op junctions present):
  ```
  L=/scratch/users/kevinroy/compass_a549/rectify_sr_smoke/logs
  grep -iE "aligner|STAR|hisat|magicblast|gsnap|bbmap|reads|error|0 reads" $L/30367486_0.out | head -60
  # both-mates check on a chunk consensus BAM:
  source .../conda.sh && conda activate rectify
  samtools view .../chunk_outputs/A549_rep1_subsampled_chunk_000_of_004.consensus.bam | \
    grep -oE "RN:i:[0-9]+" | sort | uniq -c | awk '{print $1}' | sort | uniq -c
  ```
  If both-mates and 7-aligners pass → proceed to FULL run.
- **If any task FAILED** → read `$L/30367486_<task>.err`; most likely causes: an aligner not found on PATH
  (symlink/PATH), STAR `versionGenome … INCOMPATIBLE` (regenerate index, ~30min), or OOM (bump mem again).

### FULL run (after smoke passes)
The split of the FULL 42M-pair fastq (7.5GB gz) is heavy — run it AS A JOB, not on the login node.
Generate with the same flags but input `$W/COMPASS/fastq/A549_rep1_R{1,2}.fastq.gz` and **no `-n`**
(auto-size ~500 chunks via `--target-reads-per-chunk`). Then `bash submit_pipeline.sh`. Drop a sentinel,
refresh THIS handoff.

### Deliverable (P5 step 4) — TOOLING NOT YET WRITTEN
Junction extraction from the merged consensus BAM + the 111-adjudication (positive control: annotated chr5
junctions HIGH; negative ≈0; `111 ∩ COMPASS`) has **no code yet**. WRITE IT DURING the full run.
Report the three numbers together — near-zero `111 ∩ COMPASS` only means "artifact" if the positive control passed.

## Files
- Deployed code: `/scratch/users/kevinroy/compass_a549/rectify_src/` (rsync target; re-sync after local edits).
- rectify env python: `/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify/bin/python`.
- Aligner symlinks: `$HOME/.rectify/bin/{STAR,hisat2,magicblast,gsnap}` → compass env.
- Smoke: `/scratch/users/kevinroy/compass_a549/rectify_sr_smoke/`.
- Local generator edit (uncommitted): `rectify/core/commands/split_command.py` lines 108–110.
- Source of truth for aligner params: Sherlock `…/compass_a549/COMPASS/process_reads_and_align.sh`.
