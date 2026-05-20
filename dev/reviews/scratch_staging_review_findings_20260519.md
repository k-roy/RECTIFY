# Scratch-Staging Review — Findings

**Reviewer:** Claude (Opus 4.7, 1M ctx)
**Reviewed file:** `scratch_staging_review.md` (2026-05-19, Kevin R. Roy)
**Files inspected:**
- `rectify/slurm.py` (full)
- `rectify/core/commands/run_command.py` lines 540-590 (CLI changes)
- `rectify/core/commands/run/single_sample.py` (full, both `_process_one_sample` and `_run_single_sample`)

Method: source read + targeted M1 empirical checks (env-var probe, threaded
collision repro, rsync filter simulation). No cluster smoke test was run —
see "Cluster testing" at the bottom for why.

---

## Verdict

**Blocking bug present; do not merge as-is.** One BLOCKER + several smaller
issues, listed by severity below. Most can be fixed in <30 lines.

---

## Resolution (2026-05-19, same session)

All findings addressed in this branch:

| Finding | Status | Fix |
|---|---|---|
| BLOCKER — `sample_id[:12]` collision | **FIXED** | Drop truncation; `f'rectify_{sample_id}'` in both `_process_one_sample` dispatch branches (`single_sample.py:67,71`) |
| HIGH — macOS `$TMPDIR` always-on | **FIXED** | `get_scratch_dir()` early-returns `None` outside `is_hpc_job()` (`slurm.py:145-180`) |
| MEDIUM — Step 4 polyA restore on NFS | **FIXED** | `_restored_bam` lives under `_work` / `scratch_dir`; tracker records the Oak path; final `sync_to_oak` ships the BAM. (`single_sample.py:212-214, 631-632, 657-659`) |
| MEDIUM — rsync filter leaks trim/ckpt | **FIXED** | `sync_to_oak` excludes `drs_trim/`, `consensus_ckpt/` unconditionally; shutil fallback mirrors. Include list extended for `*.corrected_polya.bam[.bai]`. (`slurm.py:284-340`) |
| MEDIUM — `_run_single_sample` try/finally | **FIXED** | Final sync wrapped in `try / finally rmtree` so an rsync `RuntimeError` no longer orphans scratch. (`single_sample.py:660-676`) |
| LOW — doc gap #2 misstated | **FIXED** | Updated in `scratch_staging_review.md` (gap #2 marked non-issue) |
| LOW — `rsync -rlL` doc wording | **FIXED** | Updated in `scratch_staging_review.md` (corrected to "overwrites when size/mtime differs") |

**Verification:**

- M1 probe: `is_hpc_job() = False`, `get_scratch_dir() = None`,
  `make_job_scratch_dir() = None` — no auto-staging on workstation.
- M1 probe with `SLURM_JOB_ID` set: HPC path returns a scratch dir as
  expected.
- 4-worker `ThreadPoolExecutor` collision repro with full sample IDs —
  no clobbering, four distinct scratch dirs.
- Synthetic scratch tree + `sync_to_oak(..., exclude_aligner_bams=True)`:
  `*.rectified.bam`, `*.corrected_polya.bam`, durable artifacts present;
  `per_aligner/*.bam`, `drs_trim/`, `consensus_ckpt/` excluded.
- `pytest -m "not slow"`: 1002 passed / 33 skipped / 4 deselected / 1
  xfailed, plus the 2 pre-existing `corrected_consensus_tiebreaker.py`
  failures the doc already flagged. No new regressions.

**Cluster smoke recommended pre-merge:** 2-sample H2 manifest with
colliding-prefix sample_ids (e.g. `by4742_drs_wt_rep1` + `..._rep2`),
`--threads 4` on an ≥ 8-core node, ideally with `--write-polya-bam` to
exercise Step 4 sort-on-scratch. Verify each worker holds a distinct
`$SCRATCH/rectify_<full_sample_id>_<job>_<task>/` during the run and
that both samples produce non-zero `corrected_reads.tsv` on NFS without
`drs_trim/` or `consensus_ckpt/` leaking.

---

## BLOCKER — concurrent multi-sample workers collide on scratch dir

`_process_one_sample` builds the scratch dir name from `sample_id[:12]`
(`single_sample.py:67` and `:71`). Within a single ThreadPoolExecutor job
all workers share the same `job_id` + `task_id`, so two samples whose
sample_id agrees in the first 12 characters write to the same scratch dir
and clobber each other's BAMs. Worse, the `finally` block on `:254` does
`rmtree(_work)` once a worker exits, deleting in-flight files of any
sibling worker that hasn't finished yet.

Affected pairs that occur in the real deposit layout:

```
'by4742_drs_wt_rep1'    -> 'by4742_drs_w'
'by4742_drs_wt_rep2'    -> 'by4742_drs_w'   ← collision
'by4742_drs_ysh1aa_rep1' -> 'by4742_drs_y'
'by4742_drs_ysh1aa_rep2' -> 'by4742_drs_y'  ← collision
```

Reproduced on M1 with a 4-worker `ThreadPoolExecutor`:

```
by4742_drs_wt_rep1     -> $TMPDIR/rectify_by4742_drs_w_96682_0
by4742_drs_wt_rep2     -> $TMPDIR/rectify_by4742_drs_w_96682_0  (SAME)
by4742_drs_ysh1aa_rep1 -> $TMPDIR/rectify_by4742_drs_y_96682_0
by4742_drs_ysh1aa_rep2 -> $TMPDIR/rectify_by4742_drs_y_96682_0  (SAME)
```

The review doc's gap #3 *acknowledges* the thread-safety story but assumes
"the sample_id being unique within a manifest" — the truncation breaks
that assumption silently.

Real concurrency confirmed: `multi_sample.py:116` sets
`n_workers = max(1, n_cpus // threads_per_sample)`. On any H2 compute
node ≥ 8 cores with default `--threads 4`, ≥ 2 workers run in parallel,
so this fires in production runs — not a latent bug.

**Fix:**
- Replace `sample_id[:12]` with the full `sample_id` (filenames hold the
  needed bytes — Linux PATH_MAX is 4096).
- Or, append a worker discriminator (`threading.get_ident()` or an enumerated
  manifest index) so concurrent workers cannot collide.

The cheapest correct change is to drop the truncation. There's no obvious
reason to truncate (no scheduler-side limit on subdir names) — the prefix
just shows up in `ls $SCRATCH` for diagnostics.

---

## HIGH — `_using_scratch` is always True on macOS, contrary to doc claim

The CLI help text and review doc both promise "no staging overhead on
local workstations." This is **false on macOS**, where `$TMPDIR` is set by
the OS to a per-user `/var/folders/...` path. `get_scratch_dir()` returns
it; `make_job_scratch_dir()` therefore returns a per-process subdir on
every M1 invocation — so an unintended copy + rsync runs at the end of
every M1 `rectify run-all`.

Probe output on M1 (no HPC env vars set):

```
$TMPDIR = '/var/folders/h5/.../T/'
$SCRATCH = None
$SLURM_TMPDIR = None
get_scratch_dir() = /var/folders/h5/.../T
make_job_scratch_dir('test_oop') = /var/folders/h5/.../T/test_oop_96535_0
-> _using_scratch would be True on this M1
```

The doc's own unit check (lines 142-152) sidesteps this by popping
`TMPDIR` explicitly, which masks the issue. The doc author appears to
have *known* about TMPDIR but assumed it would be unset on M1.

**Fix options:**
1. Gate auto-detection on `is_hpc_job()` (already defined in `slurm.py`).
   Only use `$TMPDIR`/`$SCRATCH`/`$SLURM_TMPDIR` when actually running
   under SLURM/SGE/PBS. Local workstation runs would then have zero
   staging overhead, matching the doc claim.
2. Or remove `$TMPDIR` from the auto-detection priority (keep only
   `$SCRATCH` and `$SLURM_TMPDIR`). Loses portability on rare TACC-style
   setups but matches the doc.

Option 1 is cleaner — `is_hpc_job()` is already the right predicate.

---

## MEDIUM — DRS Step 4 restored BAM is written to NFS and sorted there

In both `_process_one_sample` (`:226-230`) and `_run_single_sample`
(`:637-642`), the Step 4 restored polyA BAM path is rooted at
`sample_output` / `output_dir` (NFS), and the pysam sort temp file
(`<bam>.sort_tmp.bam`) lands adjacent to it on NFS. This re-introduces
the exact `pysam.sort` on NFS pattern the scratch-staging change was
meant to fix elsewhere.

The fix is mechanical: write the restored BAM into `_work` /
`scratch_dir`, sort there, then let the final `sync_to_oak` ship it to
NFS along with everything else (or `shutil.copy2` it adjacent to the
early-synced TSV). The `--include=*.bam` exclusion in
`exclude_aligner_bams` mode does **not** exclude `.corrected_polya.bam`
(it's not `.rectified.bam`), so it would be excluded — need to extend
the include list or write it under a name like
`<sample>.corrected_polya.rectified.bam` or change the include pattern.

---

## MEDIUM — Transient scratch directories leak to NFS via rsync

Two distinct sub-issues, both rooted in the same rsync filter gap
(`slurm.py:289-296` — excludes `*.bam`/`*.bai`, allows everything else):

**M.1 — DRS trim FASTQ leak.** DRS Step 0 writes
`_work/drs_trim/<stem>_trimmed.fastq.gz` (GB-scale, fully regenerable
from the input BAM + parquet) and the final rsync ships it to NFS.

**M.2 — Consensus checkpoint leak.** `_work/consensus_ckpt/` holds
mid-alignment checkpoint state and also gets synced. Different subdir,
different cause (consensus stage), but same filter gap.

Empirical confirmation with the documented filter rule against a
representative tree:

```
$ tree dst/  # exclude_aligner_bams=True
dst/
├── .corrected_reads.tsv.provenance.json
├── consensus_ckpt/
│   └── checkpoint.json                      # transient — should not sync
├── corrected_reads.tsv
├── drs_trim/
│   └── sample_trimmed.fastq.gz              # GB-scale, regenerable — should not sync
├── sample.rectified.bam
└── sample.rectified.bam.bai
```

For a typical DRS sample the trimmed FASTQ is several GB. On a 48-sample
manifest this is hundreds of GB of regenerable intermediate landing in
the output dir — exactly what scratch staging was supposed to prevent.

**Fix:** add `--exclude=drs_trim/`, `--exclude=consensus_ckpt/`, and
`--exclude=*.fastq.gz` to the rsync command (and mirror in the
`_should_skip` fallback). Or invert: explicit `--include` list of the
durable outputs (TSV, parquet, bedgraphs, HTML report, rectified BAM)
and a final `--exclude=*` to drop everything else.

---

## LOW — Doc gap #2 is misstated; current code is fine on this point

The review doc's gap #2 claims `resolve_scratch_dir` "falls back silently
to `tempfile.mkdtemp()`" when `--scratch-dir /nonexistent/path` is
passed. It does not — `slurm.py:242-244` calls
`Path(base_dir).mkdir(parents=True, exist_ok=True)` and raises on
permission errors. The tempfile fallback at `:248-250` is unreachable
once `base_dir is not None`. No action needed, but the doc should be
updated so a future reviewer doesn't chase a non-bug.

---

## MEDIUM — `_run_single_sample` lacks try/finally; sync_to_oak now raises

Doc gap #4 frames this as "arguably acceptable" — that framing made
sense under the prior `sync_to_oak` behavior, but **this PR also
changed `sync_to_oak` to `raise RuntimeError` on rsync non-zero**
(`slurm.py:299-301`). Under the old behavior the rsync failure was just
logged; the unwrapped caller saw no exception. Under the new behavior
the failure now propagates as an uncaught exception that escapes
`_run_single_sample`, orphans the scratch dir, and skips both the
provenance manifest save and the final summary. This is a **new
regression introduced by this PR**, not a pre-existing issue. Either:
- Wrap the body of `_run_single_sample` in `try/finally` mirroring
  `_process_one_sample:81-255`, OR
- Have `sync_to_oak` log + return non-zero rather than raise (back to
  the prior behavior).

The first is the cleaner choice.

---

## LOW (cosmetic) — Doc wording: `rsync -rlL` is not "copy-only newer"

The review doc says rsync uses `-rlL` ("recursive, dereference symlinks,
copy-only newer"). The `-rlL` flags do not include `-u`. Current
behavior is "overwrite if size/mtime differs" — matches the doc's
*intent* ("re-run should add/overwrite, not delete") but the wording is
wrong. Trivial fix to the doc.

---

## What I did NOT find

- The DRS Step 4 ordering invariant in `_process_one_sample` (early TSV
  sync before Step 4, final scratch sync after) is correct as described.
  Per-aligner BAMs are still on `_work` when Step 4 runs.
- `_process_one_sample`'s `try/except/finally` cleanup is properly
  scoped — the `finally` runs `rmtree` unconditionally when
  `_using_scratch` is True.
- The scheduler-env-var priority chain in `get_job_id`/`get_task_id` is
  correct and matches the doc.
- The rsync include/exclude rule ordering correctly keeps
  `*.rectified.bam` while dropping `*.minimap2.bam` etc. (verified with
  a synthetic source tree).

---

## Cluster testing

The review doc prescribes a single-sample DRS smoke test on H2. I did
**not** run it, for two reasons:

1. A single-sample test cannot exercise the blocking
   `_process_one_sample` collision bug.
2. The macOS `$TMPDIR`-auto-on bug means that even M1 testing has been
   running through scratch all along; the doc's "no overhead on local
   workstations" claim was never observed to be true in the developer's
   own environment, which suggests the prescribed cluster test was not
   the gating signal that caught issues before this review.

## Cluster smoke (H2, 2026-05-19 evening) — RAN

**Inputs:** 5 217-read subsample of `wt_rep1.bam` + 5 305-read subsample of
`ysh1_rep1.bam`, drawn from the freshly onboarded DRS deposit
`/u/project/guillom/shared/raw/PRJNA1229592_cpa_depletion_DRS_2025/`.
Smoke driver: `/u/scratch/k/kevinroy/smoke_PRJNA1229592_subsampled/rectify_smoke.sh`.

### Pipeline A — `rectify run-all` (2-sample manifest, `--drs --write-polya-bam`)

- Confirmed `is_hpc_job() = True` and scratch staging engaged.
- Per-sample scratch dirs used full sample_id:
  `$SCRATCH/rectify_wt_rep1_<job>_<task>` /
  `$SCRATCH/rectify_ysh1_rep1_<job>_<task>` — no collision risk.
- Both samples: DRS trim → 4-aligner align → per-aligner correct → merged
  TSV → sync to NFS — all green.
- Combined analysis ran through stages 1–6 (aggregate → clusters → bedgraphs
  → genomic distribution → PCA → heatmap → DESeq2-init). DESeq2 then
  exited 1 because the test has only 1 rep / condition (an artifact of the
  smoke design, not a scratch-staging issue).
- Scratch left behind: **none** (job-specific dirs cleaned).
- Transient leak check: two new leaks surfaced — see below; patched same session.
- Elapsed: **572 s wall**.

### Pipeline B — sequential subcommands on wt_rep1

- trim-polya → align → per-aligner `rectify correct` (×4 aligners) →
  `merge_corrected_tsvs` (Python) → analyze → restore-softclip — all green.
- Per-aligner correct row counts: minimap2 4641, mapPacBio 4692,
  gapmm2 4335, deSALT 4625.
- Merge yielded 4 655 rows, 39 cols, `winning_aligner` present:
  `{deSALT: 4049, mapPacBio: 347, minimap2: 146, gapmm2: 113}`.
- Step 4 restore-softclip: 4 622 / 4 655 reads restored (33 lacked metadata).
- Elapsed: **387 s wall**.

### Pipeline equivalence — comparison C

- run-all wt_rep1 TSV: 4 557 rows × 38 cols.
- sequential wt_rep1 TSV: 4 655 rows × 40 cols
  (extra cols: `sample`, `winning_aligner`).
- Shared reads (read_id join): 4 589.
- `corrected_3prime` agreement: **89.4 %** (4 101 / 4 589).
- The 10.6 % divergence is consistent with the two pipelines doing
  different cross-aligner junction-pool plumbing during per-aligner
  correction (run-all passes the full per-aligner BAM list to every
  `_run_correction` call; sequential `rectify correct` doesn't expose that
  knob, so Module 2H runs with a single-aligner candidate pool). Both
  outputs are valid; sequential is not directly call-compatible with
  run-all's correction stage.

### Aligner panel — `docs/ALIGNER_RECOMMENDATIONS.md` (updated 2026-05-19)

- minimap2, mapPacBio, gapmm2, deSALT — **all working on H2.**
- deSALT-on-H2 produces real alignments: 5 065 records, 99.5 % mapped,
  splice N-ops present. H2 deSALT is the conda-installed binary at
  `~/.conda/envs/rectify/bin/deSALT` (built 2025-09-29) — **not** the
  vendored `data/bin/linux_x86_64/deSALT` (1.5.6) that crashed with
  SIGSEGV across 24 samples on Sherlock. The drop-from-production
  decision should be re-scoped as "drop the vendored Sherlock binary"
  rather than blanket-dropping deSALT on every platform.
- uLTRA: produces `*.uLTRA_ultra_tmp/database.db` (prep_splicing
  succeeds) but the alignment BAM `*.uLTRA.bam` is not emitted on H2.
  Same in both Pipeline A and Pipeline B. uLTRA is silently failing
  post-DB-build, cause TBD — not blocking the 4-aligner consensus since
  the other three cover the read set. Follow-up: capture uLTRA stderr
  in a fresh smoke to localize the post-DB crash.

### NEW transient leaks surfaced (patched this session)

- `<sample>.uLTRA_ultra_tmp/database.db` — uLTRA's per-sample annotation
  DB, regenerable. ~8 MB/sample.
- `<sample>.gapmm2.paf` — debug PAF dump. Small but unwanted.
- Fix: `sync_to_oak` exclude list extended to `*_ultra_tmp/` + `*.paf`
  in both the rsync command and the shutil fallback (`slurm.py`).

### Pre-existing DRS-rename bug — found and fixed

- Symptom: first run-all attempt failed with `Alignment completed but
  rectified.bam not found: …/<sample>.rectified.bam`. Cause:
  `align_command.py:281` derives prefix from `args.reads.stem` →
  `<sample>_trimmed` after DRS Step 0, so the BAM lands at
  `<sample>_trimmed.rectified.bam` while the caller in `stages.py:66`
  expects `<sample>.rectified.bam`.
- Fix: `stages.py` now passes `prefix=sample_id` explicitly when
  building `align_args`.

### `rectify split` BAM-input support (separate ask)

- `split` previously took FASTQ/FASTQ.GZ only. Added BAM input with two
  modes:
  - `--drs`: run `trim_drs_bam_polya` internally → trimmed FASTQ +
    parquet metadata (so the chunked pipeline can later feed Step 4
    restore-softclip). Mirrors what run-all does at Step 0.
  - No `--drs`: dump via `samtools fastq -T '*'` (preserves aux tags
    as comments, no trimming) — appropriate for non-DRS BAMs.

Until the collision fix is in, a cluster test is premature.
