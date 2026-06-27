# HANDOFF — DRS validation bundle: winner-selection bug + coherent rebuild

**Date:** 2026-06-26 · **Branch:** `drs-validation-rebuild` (SHARED — a concurrent COMPASS session also
commits here; user confirmed you own the validation work). **You are a fresh session: read this top-to-
bottom, then continue.** Read `~/.claude/CLAUDE.md` + `~/work/CLAUDE.md` + `CLAUDE.md` first (rm-safety,
M1 8GB discipline, ControlMaster, surgical `git add` — never `-A`).

## ⚑ UPDATE 2026-06-26 (#3) — WALKBACK extended-tail guard (force-aligned poly-A)

Figure review (cat1_minus_1) surfaced a real walkback gap. mapPacBio force-aligns the 3' poly-A tail
(no soft-clip, reference_start 9809 vs consensus 9831), and its corrected 3' came out **9811 vs the
consensus 9834** — it anchored on a coincidental read-A/genomic-A match at the minus-strand terminal.
Root cause: the tail-context false-stop guard's near window (`tail_context_k=4`) was CLEAN poly-T
(no mismatch) because the force-aligned-tail mismatch run only starts ~8 bp inward — so the guard
never fired, and the terminal non-stop match anchored the walk immediately. (The authors had already
added the leading-stop-match relaxation FOR this cat1_minus_1 pattern; the in-loop guard just couldn't
see far enough.)

**Fix** (`rectify/core/correct/walkback.py`, BOTH right- and left-side scans): added an extended
lookahead `tail_context_far=16`. When the near k-window is stop-base-only and mismatch-free, scan up to
16 aligned positions inward; if that wider run stays **stop-base-dominated (≥2/3)** AND is
**mismatch-rich (≥1/3)**, treat the candidate as still inside the force-aligned poly-A tail and keep
scanning. The ≥1/3 mismatch gate (the demonstrated case is ~56%) separates force-aligned tails from
genuine A/T-rich 3'UTRs with a lone basecall error (which would be ~6% and must NOT be over-walked).

**Verified (full pipeline, not the walkback-only proxy):** regen diff vs HEAD = exactly **1** change —
cat1_minus_1 mapPacBio 9811→**9834**, CONVERGING to the winner; **0 winner changes** in
corrected_reads.tsv. Tests: walkback units 54 passed; test_validation_reads 108 passed; full non-slow
suite (running). Figure: mapPacBio 3' now soft-clips the poly-A and lines up with the consensus.

## ⚑ UPDATE 2026-06-26 (#2) — POLY-A ROUND-TRIP: CORRECTED DIAGNOSIS (handoff #1's was WRONG)

The handoff's poly-A diagnosis ("parquet from a different build; aligned reads UNTRIMMED ==dorado len")
is **incorrect**. Verified the real state (the prior "412 vs 399" was read as a global parquet-vs-reads
skew; it is actually a localized dorado_source defect):

- **`val` reads ARE trimmed** (val ≠ dorado for all 36). **Parquet is coherent with `val`**, not stale:
  pq_trim == committed-val for **35/36** reads; pq_orig == build-X source for **36/36**.
- The earlier "val not a substring of dorado / independent basecalls" finding was an **artifact of the
  `=`-encoded SEQ** in `validation_reads.bam` (ref-relative; `query_sequence` returns `==A=C==…`). Decode
  vs reference before any seq comparison. (Minus-strand reads store the poly-A as leading poly-T — tail at
  the START — so `dor.startswith(val)` fails even for identical molecules.)
- **The authoritative build-X untrimmed source EXISTS**: `scripts/validation_data/rebuild_2026_05/
  combined_dorado_source.bam` (36 tagged reads, literal seqs) — its lengths == pq_orig for **36/36**, and
  committed-val == build-X-trimmed for 35/36. So parquet+val came from `combined_dorado_source.bam`.
  ⚠️ It is **NOT git-tracked** (scratch artifact) — copy its content into a tracked path to make durable.

### The bundle is build-X coherent EXCEPT two localized defects
1. **`validation_reads_dorado_source.bam`**: 6 reads were re-pulled from the WRONG basecall
   (`build_dorado_source.py` uses `wt_by4742_rep1.bam`), so their length ≠ pq_orig → **breaks test1**.
   Reads: cat2_plus_1(412/399), cat1_plus_2(449/451), cat1_minus_1(741/719), cat2_minus_1(656/655),
   cat1_minus_2(619/622), cat1_plus_1(187/167). FIX = source dorado_source from `combined_dorado_source.bam`.
2. **`validation_reads.bam` cat9_minus_2 (d3357db5)**: committed val=632 vs build-X-trimmed=640 (−8) →
   **breaks test2**. SEPARATE, val-side, 1 read; the dorado fix does NOT address it. Investigate before
   touching (changing val re-aligns → could shift a position).

### FIX PLAN (advisor-endorsed; does NOT touch val positions → cannot regress cat1_minus_1)
- Replace committed `validation_reads_dorado_source.bam` with `combined_dorado_source.bam` content
  (re-sorted/indexed, XV/XG preserved) → fixes test1 for all 6. Update `build_dorado_source.py` to use
  `combined_dorado_source.bam` as the authoritative source (it currently uses wt_by4742 → re-introduces skew).
- cat9_minus_2 test2: decide separately (re-derive that one val from build-X trimmed, or accept/xfail).
- The full re-trim+re-align ("trim-rebuild" that regressed cat1_minus_1) is NOT needed and is rejected.

### STATUS 2026-06-26 — POLY-A ROUND-TRIP: **FIXED** (all 3 tests pass; markers removed → live guards)
- **dorado_source swap DONE**: committed `validation_reads_dorado_source.bam` re-sourced from
  `combined_dorado_source.bam` (build X). 6 reads corrected, 30 byte-identical, 0 position changes. test1 passes.
- **cat9_minus_2 root cause = uLTRA read-trimming** (not build skew): uLTRA's `--reduce_read_ployA`
  (default 8) DROPS query bases from its emitted SEQ (640→632) despite help claiming otherwise. Per user's
  **no-trim policy** (aligners may soft/hard-clip — penalized, then rectify rescues — but may NOT delete
  bases), fixed in `rectify/core/align/multi_aligner.py::run_ultra`: inject `--reduce_read_ployA 10_000_000`
  (constant `_ULTRA_DISABLE_POLYA_REDUCE`) unless caller overrides. Verified: only extends terminal soft-clip;
  aligned core + reference_start + all other 35 reads byte-identical (blast radius = 1 read).
  - Regenerated committed `aligners/validation_reads.uLTRA.bam` cat9_minus_2 → 640 (via
    `update_validation_aligner_bams.py`, XV/XG preserved). Updated `validation_reads.bam` cat9_minus_2 → 640
    (copied fixed SEQ/CIGAR/quals onto bundle record, all consensus tags preserved). test2 passes.
- **Removed xfail markers** on both round-trip tests (assertions UNCHANGED; now live guards).
- **Reproducibility DONE** (user approved): `scripts/validation_data/combined_dorado_source.bam` now TRACKED
  (copied out of the gitignored rebuild_2026_05/, coord-sorted+indexed); `build_dorado_source.py` repointed to
  it as the authoritative source (old wt_by4742 path demoted to lazy fallback, only checked if a read is
  unresolved). Re-running the script is now idempotent + coherent (36 mapped, 0 missing; round-trip 3/3 pass).
- **VERIFIED**: full validation suite 109 passed / 0 failed; winner+position diff pre-vs-post uLTRA fix = 0 reads
  changed. cat9_minus_2 corrected_3prime + winner unchanged (restored bases are leading soft-clip, far from CPA).
- **Uncommitted change set** (mine — surgical add only; do NOT add the peer's dev/HANDOFF_SHORTREAD_P5.md or
  dev/HANDOFF_ALIGNER_BENCHMARK.md): multi_aligner.py, build_dorado_source.py, test_validation_polya_roundtrip.py,
  validation_reads.bam(.bai), validation_reads_dorado_source.bam(.bai), aligners/validation_reads.uLTRA.bam(.bai),
  scripts/validation_data/combined_dorado_source.bam(.bai), HANDOFF.md.
- **Full regen DONE** (regen_pa_rest_bundle.py): net committed-fixture change = 2 reads only —
  cat2_minus_2 128096→**128102** (now matches the test; fixes the real staleness) and cat7_plus_2 598986→598985
  (−1, no assertion). pA-rest restore acceptance gate PASSED: restored len == dorado len for all 36 (no
  double-count; the coherent parquet resolved the old "24 vs 12"). (Earlier 8-read estimate was inflated — that
  scratch merge omitted the penalty/overhang tables; the regen uses them and is authoritative.)
- **merge_corrected_tsvs idempotency BUG fixed** (`corrected_consensus.py`): the full non-slow suite surfaced
  5 failures (test_bam_writer_parallel_smoke ×4, test_run_single_sample_consensus ×1) — `KeyError
  'min_junction_anchor'`. Cause: the merge WRITES its internal `min_junction_anchor` column into
  corrected_reads.tsv (an omission — `hp_edit_distance`/`aligned_bases` were already dropped), so re-feeding
  that output collided (`_x`/`_y`) on the ed_df join. The stale HEAD fixture (pre-feature, 39 cols) masked it.
  Fix: (1) added `min_junction_anchor` to `drop_cols` (output back to 39-col schema); (2) defensive
  `all_df.drop([...3 computed cols...], errors='ignore')` before the ed_df merge → merge now idempotent.
  Both verified (10/10 in those two files pass). Re-running regen to rewrite corrected_reads.tsv as 39-col.
- **STILL OPEN**: (a) confirm full non-slow suite green after re-regen; (b) commit (bundled); (c) figure walkthrough.

## ⚑ UPDATE 2026-06-26 (#1) — NAMED WINNER BUG IS MOOT

Traced the winner-selection bug. **It does NOT reproduce on current HEAD.** Verified four ways:
1. **Fresh merge on committed BAMs:** cat2_minus_2→**128102**, cat4_minus_1→**128117** (both correct/test-expected).
2. **Committed deSALT BAM has NO 130168 mismap** — `9dbd37bf` aligns chrI:128112 (correct). The 130168
   only exists in the peer regen dir `scripts/validation_data/regen_2026_06_25/` (a different deSALT run).
3. **Even in that peer regen** the merge's own `comparison_summary.tsv` shows deSALT@130168 `_is_winner=False`;
   gapmm2@128102 (HP-ED 87.77) `_is_winner=True`. The handoff author read 130168 off deSALT's *per-aligner*
   `corrected_reads.region_000.tsv` and assumed it propagated — the merge had already rejected it.
4. **36-read invariant check** (fresh merge on committed bundle): winner == argmin(HP-ED) for **all 36 reads**;
   the chimera gate never overrode HP-ED. `pytest -k "cat2 or cat4"` = 19 passed.

→ **No merge code change needed.** `merge_corrected_tsvs` is sound. Do NOT touch the chimera sort.

### GENUINE remaining work (all user-gated — see OPEN DECISIONS below)
The only real discrepancy is the **stale committed fixture**. A fresh merge on the committed BAMs differs from
the committed `corrected_reads.tsv` for **8/36 reads** (all small shifts −83..+6 bp; no cross-chrom, no winner
swap). 3 move TO the test-expected value (cat2_minus_2 128096→128102, cat2_plus_2 8614→8605, cat4_minus_2
76251→76254); the other 5 (cat9_plus_2 −83, cat7_plus_2 −1, cat1_minus_2 +6, cat1_plus_2 −5, cat7_minus_1 +3)
have no test assertion. Tests pass because they run fresh; only the committed fixture + the figures built from
it are stale. **Decision pending: regenerate the committed fixture (+ figures) from a fresh merge?** This is the
"surgical" path — it does NOT touch BAMs or test assertions. cat1_minus_1 is NOT in this diff (stable on the
committed bundle; it only regressed under the *trim-rebuild*, which is a separate, rejected path).

---

## THE TASK (immediate — SUPERSEDED by the UPDATE above; kept for context)
Fix a winner-selection bug in `rectify/core/consensus/corrected_consensus.py::merge_corrected_tsvs`:
for some reads it selects a **spurious deSALT mismap with a far-WORSE HP-aware edit distance** over the
correct multi-aligner consensus. **User's governing principle: HP-aware ED is authoritative — lowest
HP-ED wins.** The bug is that the merge is NOT honoring that for these reads.

### Two exemplar reads (verified)
- **cat2_minus_2** (`9dbd37bf`): minimap2/gapmm2/uLTRA → chrI:128112, corr 3′=**128102**, HP-ED=**84.77**;
  deSALT → chrI:130168, HP-ED=**711.27**. Merge WRONGLY outputs 130168. Correct = **128102** (= test-expected).
- **cat4_minus_1** (`b956f764`): minimap2/gapmm2/mapPacBio/uLTRA → chrI:128117, corr 3′=**128117**,
  HP-ED=**11.42**; deSALT → chrI:130165, HP-ED=**334.26**. Merge WRONGLY outputs 130165. Correct = **128117**.

### Why deSALT's alignment is genuinely wrong (not a duplication)
deSALT's CIGAR at 130165 is BYTE-IDENTICAL to the consensus at 128117 (`66M8D2M1D47M1D52M2D45M1I87M2I11M
2I5M2D75M500N42M`, +2048bp), BUT it has **321/432 M-block mismatches** vs the consensus's **0/432**; the
two genome windows are only **26.1% identical**. A CIGAR encodes indel STRUCTURE, not base matches — deSALT
forced the same indel skeleton onto non-homologous sequence. So HP-ED 334/711 is DESERVED; lowest-ED →
consensus is correct.

### NEXT STEP
Trace `merge_corrected_tsvs` winner logic for `9dbd37bf` and find WHY it elevates the 711-ED deSALT over
the 85-ED consensus. Candidates: deSALT secondary record @128k scored while primary's 130k POSITION
reported; the chimera/`_n_agree`/drop-chimeric path overriding the ED sort; ED computed on wrong record.
Verify with: `_cigar_hp_edit_distance` (in that file) + `HpPenaltyTable.from_tsv` (penalty tables at
`rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/`). Fix so lowest-HP-ED wins → cat2_minus_2
→128102, cat4_minus_1→128117. Then `pytest tests/test_validation_reads.py -k cat2 -q` (slow ~2min).

## VERIFIED context
- **Root problem:** the validation FIGURES were built from STALE committed fixtures
  (`rectify/data/validation/rectified/per_aligner_summary.tsv` + `corrected_reads.tsv`, frozen 2026-05-18).
  Renderer made fresh-correct (commit `386814e`, `--fresh` opt-in + staleness banner). TESTS were always
  fresh (they run `rectify correct`+merge in tmp), so the regression suite is sound.
- **Poly-A round-trip is broken:** aligned reads (`validation_reads.bam`) are UNTRIMMED (==dorado_source
  len); trim parquet is from a DIFFERENT build (orig 412 vs actual 399); restore stacks parquet tail on
  the walkback clip (24 vs 12). Guard test committed (`57f3586`,
  `tests/test_validation_polya_roundtrip.py`, 2 xfail + 1 pass) + `regen_pa_rest_bundle.py --emit-merged-tsv` fix.
- **Coherent rebuild PROTOTYPED** (trim dorado_source → matching parquet → re-align trimmed reads):
  artifacts in `/tmp/cat2_rebuild/` (trimmed.bam/fastq.gz, trim_metadata.tsv [matching!], tagged BAMs,
  corrected/, merged_corrected.tsv) + Sherlock `/scratch/users/kevinroy/cat2_rebuild_desalt/`. Diff vs
  committed: 11/36 corrected_3prime changed.

## OPEN DECISIONS (user-facing; do NOT blind-decide)
1. **Surgical fix vs full trim-rebuild.** User leaning unclear. The trim-rebuild REGRESSES cat1_minus_1
   (9834→9852; test expects 9834) and leaves cat1_plus_2 wrong — so it's NOT a clean win. The deSALT
   winner bug is INDEPENDENT of trim (the mismap wins on the committed bundle too) → can be fixed
   surgically. **Recommended: (A) surgical = winner-selection fix + parquet round-trip fix on the
   existing bundle, NOT the full trim-rebuild.** Confirm with user.
2. **7 small shifts** (±1–18bp, legit trimmed-read diffs): cat2_plus_2(→8605✓test), cat4_minus_2(→76254
   ✓test) MATCH test (committed was wrong); cat1_minus_1 REGRESSES; others have no test assertion. Do NOT
   update `test_validation_reads.py` assertions without per-read user sign-off.
3. **cat2_minus_2 committed corrected_reads.tsv=128096 but test expects 128102** — even the old bundle was
   self-inconsistent (build skew). Target = 128102.

## RESUME
1. `cat HANDOFF.md` (this file). `ls /tmp/cat2_rebuild/` — if cleared (reboot), re-derive per "Coherent
   rebuild" (trim_drs_bam_polya on dorado_source → fastq → `rectify align --Scer --aligners minimap2
   mapPacBio gapmm2 --junction-aligners uLTRA --no-consensus`; deSALT via Sherlock sbatch pattern
   `/tmp/desalt_rebuild.sbatch` job 31409555; not strictly needed for the winner-bug trace).
2. Trace+fix the winner-selection bug (see NEXT STEP). This is the priority and is independent of /tmp.
3. Re-run `pytest tests/test_validation_reads.py -k "cat2 or cat4" -q`; expect 128102 / 128117.
4. Present any test-assertion changes to the user before committing.

## FILES / STATE
- Committed this thread: `57f3586` (round-trip guard + regen fix), `386814e` (--fresh opt-in). HEAD may be
  further (concurrent peer commits). Working tree clean except this HANDOFF.md.
- Env: M1 base conda `gapmm2==25.4.5` (do not upgrade — NEW-082). deSALT cannot run on M1 (x86); use Sherlock.
- `git add` EXPLICIT paths only. Don't touch `dev/COMPASS_*`, `scripts/validation_data/regen_2026_06_25/`
  (peer/other-thread).
