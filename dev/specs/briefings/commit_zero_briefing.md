# Sonnet Briefing — Commit Zero (Profile + Audit)

**You are a Sonnet 4.6 subagent.** Your goal: produce `/Users/kevinroy/work/rectify/dev/specs/profile_results.md` — a profile attribution table for `write_corrected_bam` on Han wt_R1, plus per-stage wall-time timings for the full pipeline. This document gates ALL subsequent work in the parallelism refactor.

Read the parent spec at `/Users/kevinroy/work/rectify/dev/specs/parallel_refactor_plan.md` §3 ("Commit Zero — Profile + audit") before starting. This briefing fills in the operational details.

---

## 1. Critical rules (do not violate)

1. **You may run commands on H2 via `ssh hoffman2 'bash -lc "..."'`.** No 2FA required. Do NOT attempt to run on Sherlock (interactive auth blocks agents).
2. **Verify, don't recite.** Per `/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/CLAUDE.md`: prove every claim with a tool call (grep / import / --help). If you state "X is at line Y", you must have `Read` it. If you state "command Z takes N seconds", you must have timed it. No phantom-flag punch lists.
3. **No git commits.** This phase produces a results document only. Commit Zero changes no rectify source code.
4. **No `git add -A` or `git add .`** ever. Surgical staging only — but again, no commits this phase.
5. **Tight cluster polling.** When you submit any qsub job, poll every 270 seconds (`ScheduleWakeup delaySeconds=270` keeps the prompt cache warm). Do not walk away on a 3-hr run without polling — early-stage bugs surface within 5 min. See CLAUDE.md §"Cluster job monitoring — poll tight when the pipeline is unproven".
6. **Defer `set -u`** in any sbatch/qsub script: `source ~/.bashrc; set -u; ...` order is fatal on H2 too. Source FIRST, then set -u if you need it. See `feedback_sherlock_sbatch_set_u.md` (applies similarly on H2).
7. **PYTHONPATH must point to dev tree.** The installed `rectify` on H2 is `0.9.0` (public release); we need the dev branch `drs-validation-rebuild`. Always prepend `PYTHONPATH=/u/home/k/kevinroy/software/rectify` when invoking rectify on H2.

---

## 2. Confirmed environment state (do not re-derive)

These are verified facts. Do not waste cluster wall time re-checking them:

| Item | Value |
|---|---|
| M1 / GitHub / H2 HEAD | `7dbb1bd fix(ultra): ship SGD native GTF + fix normalizer for gffread-produced GTFs` |
| Branch | `drs-validation-rebuild` |
| Dev tree on H2 | `/u/home/k/kevinroy/software/rectify` |
| Conda env activation on H2 | `module load conda/23.11.0 && conda activate rectify` |
| Rectify binary in env | `/u/home/k/kevinroy/.conda/envs/rectify/bin/rectify` (version 0.9.0, public release — DO NOT USE; use PYTHONPATH override) |
| Canonical Han wt_R1 invocation | `/u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.1_multialigner/qsub_run_wt_R1.sh` |
| Han wt_R1 data dir | `/u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.1_multialigner/` |
| Manifest file | `manifest_wt_R1.tsv` in the same dir |
| Genome | `/u/project/guillom/shared/references/saccharomyces_cerevisiae/R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529.fasta` |
| Annotation | `/u/project/guillom/shared/references/saccharomyces_cerevisiae/R64-5-1_20240529/saccharomyces_cerevisiae_R64-5-1_20240529.gff` |
| Workload type | **QuantSeq REV / short-read** (`--short-read --dT-primed-cDNA`); NOT ONT DRS |
| Input read count | 12.49M short reads |
| Pre-correction BAM (input to `rectify correct`) | `results/wt_R1/wt_R1.rectified.bam` (453 MB) |
| Post-correction BAM | `results/wt_R1/wt_R1.rectified_corrected_3end.bam` (465 MB) |
| Bottleneck observation | `write_corrected_bam` runs single-threaded at 93.8% CPU on one core for ~2.5 hr after correction completes in ~5 min |

If any of these turns out to be wrong, **stop and report** — do not paper over a discrepancy by re-deriving in a different direction.

---

## 3. Phase Zero-A — 200k-read slice cProfile

**Goal:** decompose the per-read CPU cost of `write_corrected_bam` into four buckets so the architecture decision (§3.3 of parent spec) can be made. Wall: ~5 min.

### 3.1 Stage the slice on H2

```bash
ssh hoffman2 'bash -lc "
  set -e
  module load conda/23.11.0
  conda activate rectify
  source ~/.bashrc

  WORK=/u/scratch/k/kevinroy/profile_zero_a
  mkdir -p \$WORK
  cd \$WORK

  IN=/u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.1_multialigner/results/wt_R1/wt_R1.rectified.bam

  # 1.6% subsample -> ~200k reads. Use seed=42 for reproducibility.
  samtools view -b -s 42.016 \$IN > wt_R1_slice_200k.bam
  samtools index wt_R1_slice_200k.bam

  # Verify slice size
  N=\$(samtools view -c wt_R1_slice_200k.bam)
  echo \"Slice read count: \$N\"
  # Expect 180k-220k. If outside that range, report and stop.
"'
```

If `samtools view -c` returns outside [180000, 220000], **stop and report** — the subsample fraction may need adjustment.

### 3.2 Install py-spy in the rectify env (fallback to cProfile if it fails)

```bash
ssh hoffman2 'bash -lc "
  module load conda/23.11.0
  conda activate rectify
  pip install --quiet py-spy 2>&1 | tail -5
  which py-spy
"'
```

py-spy produces a flame-graph SVG that is much easier to read than cProfile's call tree. If install fails (network blocked, etc.), fall back to `python -m cProfile -o profile.prof -m rectify.cli correct ...` and analyze with `pstats`.

### 3.3 Run rectify correct on the slice with profiling

```bash
ssh hoffman2 'bash -lc "
  set -e
  source ~/.bashrc
  module load conda/23.11.0
  conda activate rectify
  export PYTHONPATH=/u/home/k/kevinroy/software/rectify

  WORK=/u/scratch/k/kevinroy/profile_zero_a
  cd \$WORK

  # Force single-threaded so attribution is clean. Match the canonical flags from
  # qsub_run_wt_R1.sh: --short-read --dT-primed-cDNA.
  GENOME=/u/project/guillom/shared/references/saccharomyces_cerevisiae/R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529.fasta
  GFF=/u/project/guillom/shared/references/saccharomyces_cerevisiae/R64-5-1_20240529/saccharomyces_cerevisiae_R64-5-1_20240529.gff

  # First: a single rectify correct invocation to get end-to-end timing AND the
  # output BAM that write_corrected_bam produces. Use the slice BAM as input.
  py-spy record --format flamegraph --output flame_correct.svg --rate 200 -- \
    python -m rectify.cli correct wt_R1_slice_200k.bam \
      --genome \$GENOME \
      --annotation \$GFF \
      --short-read \
      --dT-primed-cDNA \
      --threads 1 \
      --output-dir slice_correct_out/ 2>&1 | tee correct_run.log
"' 2>&1 | tail -50
```

If `rectify correct` does not have a `--threads 1` flag (verify with `--help`), use whatever the equivalent single-threaded flag is. Do NOT guess — read `rectify/core/commands/correct_command.py:create_correct_parser` if uncertain.

### 3.4 Decompose the flame graph into the four buckets

The parent spec §3.3 defines four buckets:
1. `realign_exon_blocks` (Numba global alignment) — under `rectify/core/correct/` or splice modules.
2. `set_tag` record rebuild — pysam's `pysam.libcalignedsegment.AlignedSegment.set_tag` + downstream `to_string` rewrites.
3. BGZF compression — pysam's `pysam.libchtslib.write_pysam` and below, OR `samtools` library calls.
4. Other (pysam fetch iteration, dict lookups, etc.).

Procedure:
- Download `flame_correct.svg` from H2: `scp hoffman2:/u/scratch/k/kevinroy/profile_zero_a/flame_correct.svg /tmp/`.
- For each bucket, sum the % time across matching stack frames. Use the `Read` tool on the SVG (it's text/XML) and grep for the function names.
- If py-spy unavailable, parse the cProfile output via `pstats` printed sort-by-cumulative.

Produce a table in `profile_results.md`:

```markdown
## Zero-A: cProfile attribution

Input: wt_R1_slice_200k.bam (N reads: <observed>)
Wall time: <observed seconds>
Profiling tool: py-spy <version> | cProfile

| Bucket | % CPU | Notes |
|---|---|---|
| realign_exon_blocks | XX.X | <top frames observed> |
| set_tag record rebuild | XX.X | <top frames observed> |
| BGZF compression | XX.X | <top frames observed> |
| Other (fetch, dict, ...) | XX.X | |
| **Total** | 100.0 | |

**Dominant bucket:** <name>, <%>.
**Decision per parent spec §3.3:** <Proceed with plan / Add Commit Zero-bis BGZF fix / Add batched set_tag to Commit A>.
```

### 3.5 Acceptance for Zero-A

- [ ] Slice exists, read count in expected range.
- [ ] py-spy or cProfile output exists at known path.
- [ ] Four-bucket attribution table written.
- [ ] Decision recorded.
- [ ] Output BAM `slice_correct_out/wt_R1_slice_200k.corrected.bam` exists (proves the pipeline ran to completion on the slice).

**If Zero-A FAILS** (rectify correct errors, py-spy crashes, etc.) — stop, report the error verbatim, do not proceed to Zero-B.

---

## 4. Phase Zero-B — Full Han wt_R1 stage timings

**Only proceed after Zero-A passes and the architecture decision is recorded.** Wall: ~3-4 hr cluster.

**Goal:** measure per-stage wall time and CPU% on the unprofiled, full-size Han wt_R1 run. Confirms the assumed bottleneck attribution holds at scale, and exposes any other stage (drs_trim, prescan, align, merge_corrected_tsvs, analyze) that might be the next ceiling.

### 4.1 Submit a timed qsub job

Create `qsub_zero_b_timing.sh` under `/u/scratch/k/kevinroy/profile_zero_b/`:

```bash
#!/bin/bash
#$ -N profile_zero_b
#$ -cwd
#$ -V
#$ -j y
#$ -o qsub_zero_b.log
#$ -l h_rt=8:00:00,h_data=8G
#$ -pe shared 16

source ~/.bashrc          # MUST come first
module load conda/23.11.0
conda activate rectify
export PYTHONPATH=/u/home/k/kevinroy/software/rectify

cd /u/scratch/k/kevinroy/profile_zero_b/

WT_DIR=/u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.1_multialigner
GENOME=/u/project/guillom/shared/references/saccharomyces_cerevisiae/R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529.fasta
GFF=/u/project/guillom/shared/references/saccharomyces_cerevisiae/R64-5-1_20240529/saccharomyces_cerevisiae_R64-5-1_20240529.gff

echo "=== Environment ==="
echo "Host:    $(hostname)"
echo "PID:     $$"
echo "Git SHA: $(cd /u/home/k/kevinroy/software/rectify && git log -1 --format='%h %s')"
echo "Started: $(date -Iseconds)"

# Background: per-second top sampler to capture CPU% on our process tree.
( while true; do
    date -Iseconds
    top -bn1 -u kevinroy -d 1 -p $$ -p $(pgrep -P $$ 2>/dev/null | tr '\n' ',' | sed 's/,$//')  2>/dev/null | head -30
    echo "---"
    sleep 5
  done ) > top_samples.log 2>&1 &
TOPPID=$!

# Stage timing wrapper
declare -A STAGE_START
declare -A STAGE_END

ts() { date -Iseconds; }

# Run run-all and rely on rectify's own [TIMING] log lines for stage attribution.
# Also wrap the whole thing in `time` for total wall.
ln -sf $WT_DIR/manifest_wt_R1.tsv .

time python -m rectify.cli run-all \
  --manifest manifest_wt_R1.tsv \
  --short-read \
  --dT-primed-cDNA \
  --genome $GENOME \
  --annotation $GFF \
  -o results_zero_b/ \
  --reference WT \
  --threads 16 \
  2>&1 | tee run_all.log

EXIT=$?

kill $TOPPID 2>/dev/null

echo "Exit:    $EXIT"
echo "Finished: $(date -Iseconds)"
```

Submit via `qsub qsub_zero_b_timing.sh`. Capture the job id.

### 4.2 Poll the job

Use `ScheduleWakeup delaySeconds=270` after submitting. On each wake:

1. `ssh hoffman2 'qstat -j <jobid> 2>&1 | head -20'` — running state.
2. `ssh hoffman2 'tail -30 /u/scratch/k/kevinroy/profile_zero_b/qsub_zero_b.log'` — most recent stdout.
3. `ssh hoffman2 'tail -10 /u/scratch/k/kevinroy/profile_zero_b/run_all.log'` — pipeline progress.
4. If the job has crashed (exit non-zero, qstat shows it gone but log shows error), STOP and report.
5. Otherwise, schedule the next wakeup.

Once the job completes successfully:

### 4.3 Extract per-stage timing

Grep the run-all log for `[TIMING]` lines (rectify's own per-stage instrumentation):
```bash
ssh hoffman2 'grep -E "\[TIMING\]|=== Step|Stage [0-9]" /u/scratch/k/kevinroy/profile_zero_b/run_all.log'
```

If the log doesn't have enough granularity (e.g., it lumps "correct" together rather than breaking out correction-phase vs write-bam vs sort vs index), report this gap — Commit A's tracker work may need to add more `[TIMING]` lines.

Produce the timing table in `profile_results.md`:

```markdown
## Zero-B: Per-stage wall-time on Han wt_R1

Job ID: <id>
Host: <n NNNN>
Total wall: <HH:MM:SS>
Exit status: 0
Git SHA: 7dbb1bd

| # | Stage | Wall (s) | % of total | CPU% (peak / mean) | Notes |
|---|---|---|---|---|---|
| 1 | drs_trim (if --drs) | — | — | — | not applicable for QuantSeq REV |
| 2 | prescan: variant scan | XX | X.X% | XX/XX | single-threaded today |
| 3 | prescan: junction pool | XX | X.X% | XX/XX | sequential across aligners today |
| 4 | align: per-aligner (bbmap, bwa, ...) | XX | X.X% | XX/XX | parallel across aligners today? |
| 5 | correct: region phase | XX | X.X% | XX/XX | already parallel; confirm wall ~5 min |
| 6 | correct: write_corrected_bam | XX | X.X% | XX/XX | **THE bottleneck**; confirm wall ~2.5 hr |
| 7 | restore_polya | XX | X.X% | XX/XX | applicable? |
| 8 | merge_corrected_tsvs | XX | X.X% | XX/XX | candidate next bottleneck |
| 9 | analyze: substep timings | XX | X.X% | XX/XX | bedgraph + cluster calling + DESeq2 if applicable |
| 10 | sort + index | XX | X.X% | XX/XX | |
| **Total** | | **<total_wall>** | 100% | | |
```

### 4.4 Acceptance for Zero-B

- [ ] qsub job completed with exit 0.
- [ ] Per-stage timing table populated.
- [ ] Any stage that takes >10% of total wall AND is currently single-threaded is flagged in a "next bottleneck candidates" subsection.
- [ ] Confirms or refutes the bottleneck attribution from the issue doc (correct: write_corrected_bam ~2.5 hr).

---

## 5. Final deliverable — `profile_results.md`

Write to `/Users/kevinroy/work/rectify/dev/specs/profile_results.md`. Structure:

```markdown
# Commit Zero — Profile Results

**Date:** 2026-05-19
**Git SHA profiled:** 7dbb1bd
**Sample:** Han wt_R1 (12.49M short reads, QuantSeq REV)
**Cluster:** Hoffman2, 16-core pe shared

---

## DECISION (top of doc, single line)

<one of:>
- **PROCEED with Commits A–F as written.** Bottleneck attribution: realign_exon_blocks XX.X%.
- **RESCOPE: add Commit Zero-bis (BGZF compresslevel=1) before Commit A.** Bottleneck attribution: BGZF XX.X%.
- **RESCOPE: add batched set_tag to Commit A scope.** Bottleneck attribution: set_tag XX.X%.
- **<other observed scenario>** — describe.

## Zero-A: cProfile attribution

<table from §3.4>

## Zero-B: Per-stage wall-time

<table from §4.3>

## Next-bottleneck candidates (>10% of total, currently single-threaded)

<list>

## Open observations / surprises

<anything that doesn't fit the planned architecture or that the parent spec missed>

## Artifacts produced (still on H2)

- /u/scratch/k/kevinroy/profile_zero_a/flame_correct.svg
- /u/scratch/k/kevinroy/profile_zero_a/correct_run.log
- /u/scratch/k/kevinroy/profile_zero_b/run_all.log
- /u/scratch/k/kevinroy/profile_zero_b/qsub_zero_b.log
- /u/scratch/k/kevinroy/profile_zero_b/top_samples.log
- /u/scratch/k/kevinroy/profile_zero_b/results_zero_b/wt_R1/

These are kept for re-analysis. Do not delete.
```

---

## 6. When to escalate back to the orchestrator (Opus)

Stop and hand back immediately if any of:

- `rectify correct` on the slice fails with any error not obviously a config typo.
- The slice read count is outside [180k, 220k] (sampling assumption wrong).
- The Han wt_R1 run takes >5 hr or hits the 8 hr h_rt cap.
- The bucket attribution doesn't cleanly fit the four §3.3 categories (e.g., a fifth dominant cost emerges — pysam fetch overhead, GFF reload, etc.).
- Any cross-cluster commit-status drift surfaces (M1, GitHub, H2 disagree on `drs-validation-rebuild` head).
- A stage you didn't expect dominates wall time.

Hand back via TaskUpdate marking your task `in_progress` with a comment about what blocked you. Do NOT continue past an ambiguous result.

---

## 7. After you complete

1. `TaskUpdate` task #3 → completed (Zero-A) and task #5 → completed (Zero-B).
2. Confirm `profile_results.md` exists and the DECISION line is filled in.
3. Report back to Opus with: (a) the decision, (b) any surprises, (c) the artifacts list.

The orchestrator will then update `parallel_refactor_plan.md` if the decision changes scope, and write the next commit's briefing.

**End of briefing.**
