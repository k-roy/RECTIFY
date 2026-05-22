# RECTIFY Sherlock Pipeline Handoff — 2026-05-18

**Purpose:** Full state snapshot for agent handoff. This doc tells the next agent exactly
what happened, what was done, what failed, and what needs to happen next.

**Bypass mode:** `~/.claude/settings.json` has `"defaultMode": "bypassPermissions"` — no
permission prompts. New sessions pick it up automatically.

---

## TL;DR — Critical Blocker Before Anything Else

**The constraint fix is applied, but there is a new blocker:**
Aligner binaries (`gapmm2`, `minimap2`, `mapPacBio.sh`) live in
`/home/groups/larsms/users/kevinroy/anaconda3/bin/` — the larsms group home. This path
is **not accessible on sh04 (Genoa/Bergamo) nodes**, which are the only nodes where
anaconda3 Python runs without SIGILL. Until aligners are installed in a location accessible
on sh04, all alignment jobs will fail.

**Immediate action needed (requires Kevin interactively):**
```bash
srun --partition=owners --constraint="CPU_GEN:SIE" --time=00:05:00 --mem=1G --account=larsms \
  bash -c "echo HOST=\$(hostname); ls /home/groups/larsms/users/kevinroy/anaconda3/bin/gapmm2 2>/dev/null && echo LARSMS_OK || echo LARSMS_NOT_ACCESSIBLE"
```

If `LARSMS_NOT_ACCESSIBLE`: install aligners into `$HOME/bin` or `/oak/` and update scripts.
If `LARSMS_OK`: something else caused the "not found" — recheck PATH in run_array scripts.

---

## Background: What Happened in This Session

### Problem 1: AVX-512 SIGILL (exit code 132)

anaconda3 Python in `/home/groups/larsms/users/kevinroy/anaconda3/bin/python` uses
AVX-512 instructions. The following node families do NOT have AVX-512:
- **sh03-01 through sh03-06**: AMD Rome (EPYC 7502) — `CPU_GEN:RME`
- **sh03-07 through sh03-10**: AMD Milan (EPYC 7543) — `CPU_GEN:MLN`
- **sh02-01/02/04/05**: Intel Broadwell (E5-2640v4) — `CPU_GEN:BDW`

Previously, the only documented AMD Milan nodes were sh03-07n10 and sh03-08n22 (the larsms
partition). But when jobs ran on the owners partition (which includes all of sh03), hundreds
of tasks landed on sh03-01 through sh03-10 and sh02 Broadwell nodes — all failed with SIGILL.

**Fix applied:** Added `#SBATCH --constraint="CPU_GEN:SIE|CPU_GEN:GEN|CPU_GEN:BGM|CPU_GEN:SKX"`
to ALL SBATCH scripts (see section below for file count). This routes jobs to:
- sh04 AMD Genoa (8224P) — `CPU_GEN:SIE` ← primary target, has AVX-512
- sh04 AMD Genoa (9384X) — `CPU_GEN:GEN`
- sh04 AMD Bergamo (9754) — `CPU_GEN:BGM`
- sh02 Intel Skylake-SP — `CPU_GEN:SKX`

**Note:** sh03-11 to sh03-18 (Milan GPU nodes, CPU_GEN:MLN) also excluded by this constraint.
An older successful run on sh03-11n22 (May 15) ran Python fine — it's unknown if this is a
hardware anomaly or a different code path. The constraint is the safe choice.

### Problem 2: Aligner binaries not found on sh04

After constraining to sh04 nodes, jobs failed with exit 1 because `mapPacBio.sh`, `gapmm2`,
and `minimap2` were not found. Diagnosis from stderr logs:
```
WARNING - mapPacBio not found at mapPacBio.sh, skipping
WARNING - gapmm2 not found at gapmm2, skipping
ERROR - No aligners succeeded
```

These binaries ARE in `/home/groups/larsms/users/kevinroy/anaconda3/bin/` (confirmed with `ls`)
but apparently not accessible from sh04 nodes. `deSALT` IS in `$HOME/bin/deSALT` (user home,
accessible everywhere) — deSALT tasks on sh04 may have produced valid BAMs.

### What Was Patched

All 72 SBATCH scripts now have the constraint:
```
#SBATCH --constraint="CPU_GEN:SIE|CPU_GEN:GEN|CPU_GEN:BGM|CPU_GEN:SKX"
```
Patched scripts:
- 24 × `run_array_mapPacBio.sh`
- 24 × `run_array_others.sh`
- 24 × `run_merge_aligners.sh`
- 24 × `00_split.sh`

---

## Current State of Data

### Set3 Trimmed FASTQs — ALL COMPLETE

All 12 set3 samples now have trimmed FASTQs at
`/oak/.../v3_20260429/<sample>/drs_trim/<sample>_trimmed.fastq.gz`:

| Sample | FASTQ size |
|--------|-----------|
| rrp6_rep1 | 8.7G |
| rrp6_rep2 | 8.7G |
| sen1_rep2 | 5.9G |
| wt_tfiiib_rep1 | 6.2G |
| wt_tfiiib_rep2 | 3.8G |
| wt_tfiiib_rep3 | 3.6G |
| prp28_rep1 | 9.7G |
| prp28_rep2 | 8.5G |
| prp28_rep3 | 8.4G |
| sen1_rep3 | 4.0G |
| sen1rrp6_rep2 | 5.4G |
| sen1rrp6_rep3 | 3.6G |

### Set3 Split Status

| Sample | Chunk FASTQs |
|--------|-------------|
| rrp6_rep1 | 20 ✓ |
| rrp6_rep2 | 21 ✓ |
| sen1_rep2 | 22 ✓ |
| wt_tfiiib_rep1 | 0 — split FAILED (job 25330516, exit 132 on sh03-08n22) |
| wt_tfiiib_rep2 | 0 — split FAILED (job 25330521, exit 132 on sh03-07n10) |
| wt_tfiiib_rep3 | 0 — split FAILED (job 25330526, exit 132 on sh03-07n10) |
| prp28_rep1 | 0 — split FAILED (job 25330531, exit 132) |
| prp28_rep2 | 0 — split FAILED (job 25330537, exit 132) |
| prp28_rep3 | 0 — split FAILED (job 25330542, exit 132) |
| sen1_rep3 | 0 — split FAILED (job 25330548, exit 132 on sh03-08n22) |
| sen1rrp6_rep2 | 0 — split FAILED (job 25330553, exit 132 on sh03-07n10) |
| sen1rrp6_rep3 | 0 — split FAILED (job 25330558, exit 132 on sh03-07n10) |

All 9 split scripts have now been patched with the constraint.

### Set2 Aligner Chunks — Pre-session BAMs (valid, need audit)

These existed before the current session. BAMs were verified via @PG headers
in a previous session as matching current-HEAD aligner invocations.

| Sample | mapPacBio | minimap2 | gapmm2 | deSALT | uLTRA | Expected |
|--------|-----------|----------|--------|--------|-------|---------|
| wt_rep1 | 9 ✓ | 9 ✓ | 9 ✓ | 9 ✓ | 0 ✗ | 9 |
| wt_rep2 | 8 (1 missing) | 9 ✓ | 9 ✓ | 9 ✓ | 0 ✗ | 9 |
| wt_rep3 | 6 ✓ | 6 ✓ | 6 ✓ | 6 ✓ | 0 ✗ | 6 (verify) |
| rna15_rep1 | 4 ✓ | 4 ✓ | 4 ✓ | 4 ✓ | 4 ✓ | 4 |
| rna15_rep2 | 4 ✓ | 4 ✓ | 4 ✓ | 4 ✓ | 0 ✗ | 4 |
| rna15_rep3 | 4 ✓ | 4 ✓ | 4 ✓ | 3 (1 missing) | 0 ✗ | 4 |
| ysh1_rep1 | 2 (2 missing) | 4 ✓ | 4 ✓ | 4 ✓ | 0 ✗ | 4 |
| ysh1_rep2 | 3 (2 missing) | 5 ✓ | 5 ✓ | 5 ✓ | 0 ✗ | 5 |
| ysh1_rep3 | 4 ✓ | 4 ✓ | 4 ✓ | 4 ✓ | 0 ✗ | 4 |
| wt_tfiiib_rep1 | 0 | 0 | 0 | 0 | 0 | 24 |
| wt_tfiiib_rep2 | 17 ✓ | 17 ✓ | 17 ✓ | 17 ✓ | 0 ✗ | 17 (verify) |
| wt_tfiiib_rep3 | 16 ✓ | 16 ✓ | 16 ✓ | 15 (1 missing) | 0 ✗ | 16 (verify) |

**Note on uLTRA:** uLTRA is the 5th aligner (task indices 2×N_chunks to 3×N_chunks–1 in
run_array_others.sh). None of the set2 samples have uLTRA BAMs yet — the uLTRA jobs all
failed (they were split across sh03 and sh04 nodes, all failed).

### Set2 Merged Aligners: NONE

No sample has a merged BAM in `chunks/merged_aligners/`. The rna15_rep1 merge job
(25330730) appeared to complete in sacct but the output BAM was not found — needs
investigation (possibly wrote to SCRATCH and was deleted, or wrote to wrong location).

---

## Action Plan for Next Agent

### Step 0 (BLOCKER): Verify aligner accessibility on sh04

Kevin must run this check interactively from Sherlock:
```bash
# Start an interactive shell on a Genoa node
srun --partition=owners --constraint="CPU_GEN:SIE" --time=00:10:00 --mem=4G \
     --account=larsms --pty bash

# Inside the node:
hostname  # should be sh04-*
ls /home/groups/larsms/users/kevinroy/anaconda3/bin/gapmm2   # test larsms home
ls $HOME/bin/deSALT                                           # test user home
which mapPacBio.sh                                            # test PATH
```

**If `/home/groups/larsms/` is NOT mounted on sh04:**
Need to install gapmm2, minimap2, mapPacBio.sh into `$HOME/bin/` (user home) or `/oak/`.
Then update the PATH line in all run_array scripts:
```bash
# Current PATH in scripts (larsms home):
export PATH="/home/groups/larsms/users/kevinroy/anaconda3/bin:$HOME/bin:$HOME/.rectify/bin:$PATH"
# Fix: add /oak/ path or $HOME/bin where aligners were installed
```

**If `/home/groups/larsms/` IS accessible on sh04:**
The "not found" may have been from a different issue. Add `echo "which gapmm2=$(which gapmm2)"`
to the run_array script header as a diagnostic and resubmit one test task.

### Step 1: Fix split jobs for 9 missing set3 samples

After confirming aligners work on sh04, resubmit splits for the 9 failed samples.
The 00_split.sh scripts are now constraint-patched. Submit them:

```bash
BASE="/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set3_anchor_away"
for s in wt_tfiiib_rep1 wt_tfiiib_rep2 wt_tfiiib_rep3 prp28_rep1 prp28_rep2 prp28_rep3 sen1_rep3 sen1rrp6_rep2 sen1rrp6_rep3; do
  # Cancel old dependent pending split jobs first if any exist
  sbatch "$BASE/$s/00_split.sh"
done
```

The alignment arrays will be re-submitted manually after split completes (or chain with
`--dependency=afterok:$SPLIT_ID` by re-submitting from the submit_all_alignment_only.sh).

### Step 2: Fix set3 rrp6/sen1_rep2 alignment (already have chunks, just need arrays)

rrp6_rep1 (20 chunks), rrp6_rep2 (21 chunks), sen1_rep2 (22 chunks) have chunk FASTQs but
no alignment BAMs yet. Their run_array scripts are in:
`/oak/.../set3_anchor_away/<sample>/chunks/run_array_mapPacBio.sh`
`/oak/.../set3_anchor_away/<sample>/chunks/run_array_others.sh`

These scripts ARE constraint-patched. Once aligner availability on sh04 is confirmed, submit:
```bash
sbatch set3_anchor_away/rrp6_rep1/chunks/run_array_mapPacBio.sh
sbatch set3_anchor_away/rrp6_rep1/chunks/run_array_others.sh
# ... etc for rrp6_rep2, sen1_rep2
```

### Step 3: Fix set2 missing/incomplete alignment

For samples with incomplete mapPacBio/uLTRA chunks, resubmit specific array indices
(now safe with constraint patches). Samples needing attention:
- **wt_tfiiib_rep1**: 0 BAMs — needs full run_array_mapPacBio (0-23) and run_array_others (0-95)
- **wt_rep2**: missing mapPacBio chunk index 5 — `sbatch --array=5 run_array_mapPacBio.sh`
- **ysh1_rep1**: missing mapPacBio chunks 0-1 — `sbatch --array=0-1 run_array_mapPacBio.sh`
- **ysh1_rep2**: missing mapPacBio chunks 0-1 — `sbatch --array=0-1 run_array_mapPacBio.sh`
- **rna15_rep3**: missing deSALT chunk 15 — `sbatch --array=15 run_array_others.sh` (but array index ≠ chunk index; verify)
- **wt_tfiiib_rep3**: missing deSALT chunk — investigate which chunk is missing
- **uLTRA for all samples**: run_array_others.sh uLTRA tasks (indices 2×N to 3×N−1) for all set2 samples

### Step 4: Submit merge_aligners after alignment completes

Once all 5 aligner BAMs are present for a sample, submit:
```bash
sbatch /oak/.../set2_cpa_machinery/<sample>/chunks/run_merge_aligners.sh
```
Use `--dependency=afterok:<last_align_job>` if chaining from the resubmitted arrays.

### Step 5: Investigate rna15_rep1 merge output location

Job 25330730 (rna15_rep1 merge) completed (sacct exit 0) but no BAM found in
`chunks/merged_aligners/`. Check the run_merge_aligners.sh output path and SCRATCH
to see if the BAM was written there.

---

## Key File Paths

```
# Production root
BASE="/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429"

# Set2 samples
SET2="$BASE/set2_cpa_machinery"
# Samples: wt_rep1 wt_rep2 wt_rep3 rna15_rep1 rna15_rep2 rna15_rep3
#          ysh1_rep1 ysh1_rep2 ysh1_rep3 wt_tfiiib_rep1 wt_tfiiib_rep2 wt_tfiiib_rep3

# Set3 samples
SET3="$BASE/set3_anchor_away"
# Samples: rrp6_rep1 rrp6_rep2 sen1_rep2 wt_tfiiib_rep1 wt_tfiiib_rep2 wt_tfiiib_rep3
#          prp28_rep1 prp28_rep2 prp28_rep3 sen1_rep3 sen1rrp6_rep2 sen1rrp6_rep3

# Trimmed FASTQs (one level above set subdirectories)
# Format: $BASE/<sample>/drs_trim/<sample>_trimmed.fastq.gz

# Per-sample scripts
# $BASE/<set>/<sample>/00_split.sh
# $BASE/<set>/<sample>/chunks/run_array_mapPacBio.sh
# $BASE/<set>/<sample>/chunks/run_array_others.sh
# $BASE/<set>/<sample>/chunks/run_merge_aligners.sh
# $BASE/<set>/<sample>/chunks/logs/<jobid>_<taskid>.{out,err}

# Aligner binaries (confirm accessibility on sh04!)
# /home/groups/larsms/users/kevinroy/anaconda3/bin/ → gapmm2, minimap2, mapPacBio.sh
# $HOME/bin/ → deSALT (confirmed accessible all nodes)
```

---

## uLTRA Array Index Calculation

Run_array_others.sh covers 4 aligners: `["minimap2","gapmm2","uLTRA","deSALT"]`
For N_chunks chunks, total tasks = 4 × N_chunks:
- minimap2: tasks 0 to N−1
- gapmm2: tasks N to 2N−1
- **uLTRA: tasks 2N to 3N−1** ← these are the ones missing for all set2 samples
- deSALT: tasks 3N to 4N−1

---

## Node Architecture Reference (Sherlock)

| Node family | CPU | AVX-512 | larsms home accessible |
|-------------|-----|---------|------------------------|
| sh03-01 to sh03-10 | AMD Rome/Milan | NO | YES |
| sh04 (CBASE,CPERF,CSCALE) | AMD Genoa/Bergamo | YES | UNKNOWN — test! |
| sh02 Skylake-SP (SKX) | Intel Xeon 5118 | YES | YES |
| sh02 Broadwell (BDW) | Intel Xeon E5 | NO | YES |
| sh03-07n10, sh03-08n22 (larsms) | AMD Milan | NO | YES |

---

## Memory Notes (update if state changes)

- `project_rectify.md` memory has set2/set3 status — will be stale
- `feedback_sherlock_avx512_trap.md` memory only listed 2 nodes — needs update to say ALL sh03 CPUs (and sh02 BDW) lack AVX-512
- The constraint approach is now the canonical fix (not `--exclude`)
