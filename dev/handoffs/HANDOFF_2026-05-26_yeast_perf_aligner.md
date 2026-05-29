# Handoff — yeast perf + aligner-swap-in track (2026-05-26)

**Branch:** `drs-validation-rebuild` · **HEAD:** `1190129` (the *human-data track* keeps
advancing HEAD with branding/calibration/aligner-doc commits — that's fine).
**Sherlock** `$SR/rectify` (`/oak/.../software/rectify`) is synced to committed HEAD.

> ⚠️ **CONCURRENT HUMAN-DATA TRACK is active** (winnowmap2/minisplice on a549,
> `scripts/calibration/`, `docs/aligners/*`, and it edits `dev/PERF_AUDIT.md` +
> `dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`). **Do NOT edit/commit
> its files. Surgical staging only — `git add <explicit paths>`, never `git add -A`.**
> Sherlock ControlMaster: Kevin keeps it live; it goes stale intermittently
> (symptom: `Session open refused by peer` while `ssh -O check` says "running") —
> re-establish with `ssh sherlock -t 'bash --norc --noprofile'` (Duo).

## 1. What was done (this track)

- **rRNA/Pol III exclusion in `correct`** — committed **`cf5ebb9`**. Reads in
  rDNA/tRNA/SNR6/RDN5 skip the Module-2F 3'SS rescue (the `_hp_edit_distance`
  hotspot); default-on, `--include-rdna` to disable; threaded through the parallel
  worker, gated per-read in `correct_read_3prime`. Tests: `tests/test_correct_rdna_exclusion.py`
  (3) + 135 focused green.
- **Verified at correct-stage scale (DIAGG, job 26185659):** exclusion fires
  (610 regions skipped), **no OOM, advances** — but throughput is STILL limited by
  the large-5'-soft-clip archetype (NOT rDNA; chrI high-cov mRNA). The full
  manifest run-all (26167419) was **cancelled** — its blocker was mapPacBio-on-DRS
  slowness (separate issue), and DIAGG already proved the correct-stage fix.
- **Perf profiling #12/#13** — committed `befbc3a` (`dev/profile_correct_reads.py`,
  `dev/perf_panel_t13_results/`) + `c0702c4`/later (`PERF_AUDIT.md` Parts IV–V).
  Top-1% expensive reads are ~100% in the `_hp_edit_distance` DP. Archetypes
  (panel-time share): **large_5prime_clip 35.9% (#1)**, high_cigar_op 28.5%,
  productive 28.4%, over_cap_candidates 5.4% (17.4% gapmm2 → K=25 still binds),
  dead-end 0.8%, `excluded_BUG 0.0%` (positively confirms the rRNA gate).
- **Aligner swap-in eval (task #14)** — IN FLIGHT, see §3.

## 2. Verified

- Exclusion: 135 focused + 3 behavior tests; DIAGG no-OOM + advancing.
- **minisplice has only ONE model: `vi2-7k` (vertebrate+insect, NO fungal)** —
  confirmed via web (lh3/minisplice) + Sherlock (only `vi2-7k.kan` present). Yeast
  is out-of-distribution → expect minisplice to underperform on yeast; **winnowmap2
  is the cleaner yeast swap-in candidate.** minisplice/`--spsc` helps where splice
  signals are weak (metazoans); yeast's signals are strong → low marginal value.
- **NOT verified:** full-scale `correct` throughput (the large-clip cap isn't
  implemented); winnowmap2 yeast alignment (its meryl build was blocked — workaround
  in flight, see §3).

## 3. Open items

- **Aligner swap-in eval (task #14) — RESULTS IN for 4 aligners; winnowmap2 still
  blocked.** Eval dir `$SCRATCH/aligner_swapin_eval_yeast_20260526/`; comparison in
  `analyze_FULL.txt` (job `26231454`, COMPLETED). Per-aligner on upf1Δ common-read set:
    - **minisplice_mm2 ≈ minimap2 almost EXACTLY** (introns/read 0.167 vs 0.167;
      %annot 94.1% vs 93.8%; pairwise Jaccard 0.62 — the highest). With the vertebrate
      `vi2-7k` model on yeast, minisplice's `--spsc` barely shifts minimap2 →
      **minisplice is effectively redundant with minimap2 on yeast** (confirms the
      caveat; not a useful distinct panel member).
    - minimap2 cleanest (93.9% annotated), deSALT 91.3%, **gapmm2 noisy** (0.232
      introns/read but only 77.7% annot — 28% unique/novel).
    - **winnowmap2: STILL NO DATA — a SECOND failure.** The meryl workaround SUCCEEDED
      (`meryl_count_exit=0`, rep-kmers built, 1379 lines), but the winnowmap
      **alignment** step then failed (`winnowmap2_align_exit=1`, no BAM). Needs
      diagnosis: read `$SCRATCH/aligner_swapin_eval_yeast_20260526/winnow_26231454.err`
      (the winnowmap invocation itself, post-meryl). This is the key candidate — its
      result is the main open question for the eval.
    - uLTRA skipped (needs `--annotation`; reference aligner only).
  - NOTE: `analyze_swapin.py` needs the **conda env** (`conda activate rectify`) — a
    bare login-node `python` gives a SyntaxError; run it via sbatch or after activating.
- **[RELAY TO HUMAN TRACK] meryl bug** — `rectify/core/align/multi_aligner.py:1469`
  builds `meryl count k=15 output=<db> <genome>`; the installed meryl rejects
  `output=<db>` (wants `output <db>` separate). Breaks the FIRST from-scratch
  winnowmap2 build on any genome (their a549 only "works" via a cached DB). Fix:
  `f'output={meryl_db}'` → `'output', str(meryl_db)`. NOT fixed in code (human track
  owns it); workaround pre-builds the rep-kmers.
- **[DEFERRED — #13 top perf lever] Cap `rescue_seq` length fed to `_hp_edit_distance`
  (~60 bp)** — targets the large_5prime_clip class (35.9%). UNVERIFIED; in
  PERF_AUDIT Part IV/V. Correctness-sensitive (capping changes which rescues fire) →
  needs the verification discipline. **Awaiting Kevin's go to implement.**
- **[ORIGINAL GOAL — still pending] upf1Δ validation-set doubling.** 2 +/- reads per
  category into `rectify/data/validation/`, from upf1Δ. The `align_upf1d_nompb/` BAMs
  could feed a `correct` run (exclusion on) → harvest. cat5 hint: 2-intron genes
  YGL076C / YBR111W-A / YGL033W.
- #7 2F mapPacBio case-1 screening (deferred); #6 Commit-Zero-B profiling (pending).

## 4. Resume command

**First:** `ssh sherlock 'echo ok'` — if it fails with `Session open refused`, the CM
is stale; ask Kevin to re-establish before any cluster work.

1. **Harvest the aligner eval:**
   `ssh sherlock 'bash --norc --noprofile -c "sacct -j 26231454 --format=State,ExitCode,Elapsed; echo ---; grep -E \"meryl_count_exit|winnowmap2_align_exit|analyze_full_exit\" $SCRATCH/aligner_swapin_eval_yeast_20260526/winnow_*.out; echo ---; tail -60 $SCRATCH/aligner_swapin_eval_yeast_20260526/analyze_FULL.txt; echo ===PARTIAL===; tail -50 $SCRATCH/aligner_swapin_eval_yeast_20260526/analyze_partial.txt"'`
   - If `26231454` COMPLETED + `upf1d_rep1.winnowmap2.bam` present + `analyze_FULL.txt`
     populated → read it, conclude **winnowmap2 vs minisplice as mapPacBio swap-ins**
     (minisplice carries the vertebrate-model caveat → likely poor on yeast).
   - If winnowmap2 still failed → check `winnow_26231454.err` (meryl syntax); the
     workaround pre-builds `align_upf1d_nompb/S288C..._repetitive_k15.txt`.
   - Report the eval to Kevin.
2. **Then await Kevin's direction** on the deferred decisions: implement the
   `rescue_seq`-length cap (#13)? proceed to the upf1Δ validation doubling?

## 5. Files / paths

- **Committed (this track):** `rectify/core/bam/{bam_processor,parallel}.py`,
  `rectify/core/commands/correct_command.py` (exclusion `cf5ebb9`);
  `tests/test_correct_rdna_exclusion.py`; `dev/profile_correct_reads.py`,
  `dev/categorize_perf_panel.py`, `dev/perf_panel_t13_results/` (`befbc3a`);
  `PERF_AUDIT.md` Parts IV–V.
- **Sherlock eval dir:** `$SCRATCH/aligner_swapin_eval_yeast_20260526/` (sbatches,
  `analyze_swapin.py`, `align_upf1d_nompb/` BAMs, `analyze_partial.txt`,
  `analyze_FULL.txt` [pending], `drs_trim_upf1d/upf1d_rep1_trimmed.fastq.gz`).
- **Diag artifacts:** `$SCRATCH/rectify_runall_diag_20260524/`,
  `$SCRATCH/rectify_wt_by4742_rep1_25846844_0/`,
  `$SCRATCH/wt_upf1_manifest_test_20260523/`.
- **DO NOT TOUCH (human track):** `dev/PERF_AUDIT.md` (shared — append carefully if at
  all), `dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`,
  `scripts/calibration/`, `docs/aligners/`, the winnowmap2/minisplice code in
  `multi_aligner.py`, `dev/analyze_ab.py`.
- **Gotchas:** robust polls MUST distinguish ssh-failure (retry) from job-terminal
  (the false "GONE" bug burned a wrong "verified" claim). Unquoted `(` in `echo`
  inside nested SSH breaks the command. M1 has 8 GB RAM — heavy work on the cluster.
