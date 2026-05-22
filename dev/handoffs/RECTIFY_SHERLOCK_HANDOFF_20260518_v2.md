# RECTIFY Sherlock Pipeline Handoff — 2026-05-18 (v2, post-afternoon session)

**Purpose:** State snapshot after agent session that completed the full pipeline submission for all set2 and set3 samples.

**Bypass mode:** `~/.claude/settings.json` has `"defaultMode": "bypassPermissions"`.

---

## TL;DR — What Happened This Session

All blockers from the previous handoff are resolved. **605 jobs are in the Sherlock queue.** The full RECTIFY pipeline (alignment → merge → prescan → correct → chunk_merge → final_merge → consensus) is now chained for all 23 samples in set2 and set3.

---

## Root Causes Fixed This Session

### 1. uLTRA not found on sh04 (THE real blocker)

`~/.rectify/bin/` existed but was empty — `rectify install-aligners` was never run on Sherlock. All run_array scripts have `$HOME/.rectify/bin` in PATH but uLTRA only exists in the `rectify` conda environment (`anaconda3/envs/rectify/bin/`), not in `anaconda3/bin` (base env).

**Fix:** `ln -sf /home/groups/larsms/users/kevinroy/anaconda3/envs/rectify/bin/uLTRA ~/.rectify/bin/uLTRA`

Confirmed working: `uLTRA --version` → `uLTRA 0.1`

Other aligners (mapPacBio.sh, gapmm2, minimap2) ARE in `anaconda3/bin` (base env). `deSALT` is in `$HOME/bin`. So only uLTRA needed the symlink.

### 2. PATH missing larsms home in set3 run_array scripts

Set2 scripts had the correct PATH (`/home/groups/larsms/...anaconda3/bin:...`), but set3 scripts only had `$HOME/bin:$HOME/.rectify/bin:$PATH` — missing the larsms anaconda bin. Also `wt_tfiiib_rep1` in set2 was missing it.

**Fix:** Patched 26 run_array scripts to include the larsms path as first PATH entry.

### 3. AVX-512 constraint not on downstream pipeline scripts

The previous session patched the 72 alignment + split + merge scripts, but the 240 downstream scripts (prescan, correct, chunk_merge, final_merge, consensus) had `--partition=larsms,owners` but no `--constraint` directive.

**Fix:** Added `#SBATCH --constraint="CPU_GEN:SIE|CPU_GEN:GEN|CPU_GEN:BGM|CPU_GEN:SKX"` to 240 downstream scripts.

### 4. rna15_rep1 merge: was in wrong directory

The handoff reported "no BAM found" for the rna15_rep1 merge. The audit script checked `merged_aligners/` but the `run_merge_aligners.sh` writes to `merged_bams/`. The merge completed successfully on 2026-05-18 12:50 (job 25330730, sh04-04n31). All 5 aligner BAMs are present.

---

## Current State: Pipeline Submission Complete

### Set2 (11 samples) — all have full pipeline chains

| Sample | Alignment | Full chain submitted |
|--------|-----------|---------------------|
| wt_rep1 | 9/9 mpb/mm2/gmm/dsl + uLTRA running | ✅ merge→consensus |
| wt_rep2 | 8/9 mpb (task 5 running) + uLTRA | ✅ |
| wt_rep3 | 6/6 all + uLTRA | ✅ |
| rna15_rep1 | ✅ ALL DONE including merge | ✅ prescan→consensus queued |
| rna15_rep2 | 4/4 + uLTRA | ✅ |
| rna15_rep3 | 3/4 deSALT (task 15 running) + uLTRA | ✅ |
| ysh1_rep1 | 2/4 mpb (0-1 running) + uLTRA | ✅ |
| ysh1_rep2 | 3/5 mpb (2-3 running) + uLTRA | ✅ |
| ysh1_rep3 | 4/4 + uLTRA | ✅ |
| wt_tfiiib_rep1 | FULL alignment running (24 chunks) | ✅ |
| wt_tfiiib_rep2 | 17/17 + uLTRA | ✅ |
| wt_tfiiib_rep3 | 15/16 deSALT (task 53 running) + uLTRA | ✅ |

### Set3 (12 samples) — all have full pipeline chains

All 12 samples had existing chunk FASTQs (splits were done in an earlier session, despite the failed re-submit jobs). Alignment now running for all 12.

| Sample | Chunks | Alignment | Full chain |
|--------|--------|-----------|------------|
| rrp6_rep1 | 20 ✓ | running | ✅ through run_merge_consensus.sh |
| rrp6_rep2 | 21 ✓ | running | ✅ |
| sen1_rep2 | 22 ✓ | running | ✅ |
| wt_tfiiib_rep1 | 24 ✓ | running | ✅ |
| wt_tfiiib_rep2 | 17 ✓ | running | ✅ |
| wt_tfiiib_rep3 | 16 ✓ | running | ✅ |
| prp28_rep1 | 24 ✓ | running | ✅ |
| prp28_rep2 | 21 ✓ | running | ✅ |
| prp28_rep3 | 20 ✓ | running | ✅ |
| sen1_rep3 | 15 ✓ | running | ✅ |
| sen1rrp6_rep2 | 21 ✓ | running | ✅ |
| sen1rrp6_rep3 | 15 ✓ | running | ✅ |

### rna15_rep1 pipeline chain (leading sample)

```
Prescan:         25339206
Correct (5×):    25339207–25339211
Chunk merge:     25339212
Final merge:     25339213
Consensus:       25339214
Merge consensus: 25339215
```

---

## Pipeline Stage Reference (set2 vs set3)

**set2** pipeline:
`align → merge_aligners → prescan → correct (×5) → chunk_merge → final_merge → consensus_per_chunk → merge_consensus_chunks`

**set3** pipeline (slightly different):
`align → merge_aligners → prescan → correct (×5) → chunk_merge → final_merge → merge_consensus`

Note: `run_array_consensus_per_chunk.sh` and `run_merge_consensus_chunks.sh` don't exist in set3. Set3 uses `run_merge_consensus.sh` directly.

---

## Key Infrastructure Fixes (Permanent)

- `~/.rectify/bin/uLTRA` symlink → `anaconda3/envs/rectify/bin/uLTRA`
- All 288+ Sherlock SBATCH scripts now have `--constraint="CPU_GEN:SIE|CPU_GEN:GEN|CPU_GEN:BGM|CPU_GEN:SKX"` (AVX-512 safe)
- All run_array scripts now have `export PATH="/home/groups/larsms/users/kevinroy/anaconda3/bin:..."` as first PATH entry
- `merged_bams/` (not `merged_aligners/`) is the correct directory for per-sample merged aligner BAMs

---

## What's Next

1. **Monitor jobs** — the alignment arrays will complete in 1-6 hours. Then merge_aligners jobs will auto-trigger, followed by prescan, correct, etc.
2. **Check for failures** — use `sacct` or check log files if any jobs fail. Key things to watch:
   - uLTRA correct jobs (new, first time running with symlink fix)
   - wt_tfiiib_rep1 full alignment (large: 24 chunks × 5 aligners = 120 tasks)
3. **rna15_rep1 will finish first** — prescan is starting now; pipeline should complete in ~2-4 hours
4. **Final analysis** — once all samples have `merge_consensus` output, run `rectify analyze` at the set level

---

## Key File Paths

```bash
BASE="/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429"
SET2="$BASE/set2_cpa_machinery"
SET3="$BASE/set3_anchor_away"
UTRA_SYMLINK="$HOME/.rectify/bin/uLTRA"  # → anaconda3/envs/rectify/bin/uLTRA

# Per-sample pipeline output
$BASE/<set>/<sample>/chunks/merged_bams/<sample>_trimmed.{mapPacBio,minimap2,gapmm2,uLTRA,deSALT}.bam
```

---

## Memory Notes

- Update `project_rectify.md` with set2/set3 pipeline status
- Update `feedback_sherlock_avx512_trap.md` — the canonical fix is the `--constraint=` line; the secondary issue (uLTRA not in PATH) is now fixed via symlink
- Note: `rectify install-aligners` ran but created an empty `~/.rectify/bin/` — it doesn't auto-install from conda envs; the symlink is the correct fix
