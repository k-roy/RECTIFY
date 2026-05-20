# Session Handoff — BAM-provenance sidecar complete; H2 Han2023 ready to submit

**Date:** 2026-05-18
**Branch:** `drs-validation-rebuild` (HEAD `3f46a7f` — provenance commit)
**Prior handoff:** `handoffs/HANDOFF_2026-05-18_run_all_prelim_and_provenance.md`

---

## 1. What was done

- **BAM-provenance sidecar — complete implementation** (`3f46a7f`)
  - New module `rectify/utils/bam_provenance.py`: `compute_run_provenance`,
    `write_sidecar`, `read_sidecar`, `matches_strict`, `get_aligner_version`,
    `expected_provenance_for_aligner`, `get_rectify_git_sha`
  - `'consensus'` special-case in `get_aligner_version` returns
    `rectify.__version__` (no subprocess) — ensures write (consensus.py) and
    reuse gate (stages.py) see identical version strings
  - `--trust-existing-bams` flag added to `rectify run-all` (run_command.py)
  - Sidecar write wired into `rectify/core/consensus/consensus.py` (after
    coordinate-sort + index, wrapped in try/except)
  - Reuse gate wired into `rectify/core/commands/run/stages.py`
    (`_run_alignment`): strict match → reuse; mismatch → re-align
  - Provenance params plumbed through `rectify/core/commands/run/helpers.py`
    (`_collect_per_aligner_bams`) and `single_sample.py`
  - 21 unit tests in `tests/test_bam_provenance.py` (all pass)

- **H2 Han 2023 deposit created** (on Hoffman2, NOT committed — remote files)
  - `/u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.0/`
    — `manifest.tsv` (4 samples), `qsub_run.sh`, `README.md`, `results/` (empty)
  - Mode: `--short-read --dT-primed-cDNA`, R64-5-1 genome, HPC yeast annotation
  - **Status: ready to submit** (see Resume command below)

- **Sherlock DRS Set2 dry-run artifacts created** (on Sherlock, NOT committed)
  - `set2_cpa_machinery/dry_run_20260518/`: `manifest.tsv` (6 of 12 ready
    samples), `symlink_layout.sh`, `run_all.cmd`, `README.md`, `results/`
  - `run_all.cmd` includes `--trust-existing-bams` (BAMs predate sidecar)
  - **Status: ON HOLD** — do not submit until prescan ambiguity-window fix
    is merged (DRS_CPA_PROJECT_STATUS.md §7)

---

## 2. What's verified

- `pytest tests/test_bam_provenance.py` → **21 passed** (M1, 2026-05-18)
- `pytest -m "not slow"` → NOT VERIFIED this session (no full run after
  commit `3f46a7f`; last full run was 934 passed + 28 skipped before this
  branch; tests added since then: test_bam_provenance.py + xfail expansion)
- H2 manifest format verified via SSH:
  `ssh hoffman2 'cat .../han_2023_.../manifest.tsv'` → 4 correct
  tab-separated rows with correct FASTQ paths
- Sherlock sync verified: `git log --oneline -3` on Sherlock after
  fast-forward → confirmed at `3f46a7f`

NOT VERIFIED: end-to-end `rectify run-all` run with new provenance code
(requires cluster submit; deferred until Han2023 job runs on H2).

---

## 3. Open items

- **H2 Han2023 run not yet submitted.** Why deferred: provenance feature
  needed to land first so the run produces sidecars from the start. Feature
  is now at HEAD; submit is the immediate next step.

- **Sherlock Set2 on hold** (do not submit). Why: prescan ambiguity-window
  fix described in DRS_CPA_PROJECT_STATUS.md §7 must merge before `rectify
  correct` is authorized for DRS CPA samples. The hold is explicit.

- **`--skip-alignment` layout mismatch on Sherlock.** Why deferred: the
  dry-run is already on hold for an unrelated reason, so no urgency. But
  before lifting the hold, the manifest `path` column (currently sample
  directories) must be reconciled with how `detect_input_type` and
  `_collect_per_aligner_bams` work when `--skip-alignment` is set.
  Documented in `dry_run_20260518/README.md` under "To investigate before
  lifting the hold."

- **Retroactive sidecar audit tool** (deferred workstream from handoff
  `HANDOFF_2026-05-18_run_all_prelim_and_provenance.md` §3.3). Why
  deferred: not needed until there are existing production rectified BAMs
  that need sidecar back-fill. Set2 run hasn't happened yet.

- **`pytest -m "not slow"` full suite** — should be run on M1 or H2 before
  submitting the Han2023 job, to confirm no regressions from the provenance
  plumbing in consensus.py, stages.py, helpers.py, single_sample.py.

- **WIP working tree on `drs-validation-rebuild`** — many files staged but
  not committed (calibration scripts, validation BAMs, figure generators,
  etc.). These are pre-existing WIP from other workstreams; do NOT commit
  them as part of the provenance feature.

- **Sherlock stash** — `git stash list` on Sherlock has stash@{0} (pre-sync
  WIP from this session) plus older stashes @{1}-@{5}. Do NOT drop stash@{0}
  without reviewing it; it may contain unrelated work in progress.

---

## 4. Resume command

**H2 Han2023 (authorized — submit now):**
```
ssh hoffman2 'cd /u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.0 && qsub qsub_run.sh'
```
Then poll with:
```
ssh hoffman2 'qstat -u kevinroy'
```
When complete, check `results/` for per-sample output dirs and
`corrected_reads.tsv` files. If the job fails, pull the `.e<jobid>` error
log from the deposit dir.

**Before submitting**, optionally run `pytest -m "not slow"` on M1 or H2 to
confirm no regressions from the provenance plumbing:
```
cd /Users/kevinroy/work/rectify && pytest -m "not slow" 2>&1 | tail -5
```

**Sherlock Set2 (on hold — check gate first):**
Read `DRS_CPA_PROJECT_STATUS.md §7`. If the ambiguity-window fix is merged
and authorized, resolve the `--skip-alignment` layout question documented
in `dry_run_20260518/README.md`, update `run_all.cmd` if needed, then
`sbatch run_all.cmd` from the `dry_run_20260518/` directory on Sherlock.
If gate still closed, leave it.

---

## 5. Files touched

**Committed in `3f46a7f` (branch `drs-validation-rebuild`):**
- `rectify/utils/bam_provenance.py` — new module (complete provenance API)
- `tests/test_bam_provenance.py` — new file, 21 tests
- `rectify/core/commands/run_command.py` — `--trust-existing-bams` flag
- `rectify/core/commands/run/helpers.py` — provenance params + filtering
- `rectify/core/commands/run/stages.py` — reuse gate in `_run_alignment`
- `rectify/core/consensus/consensus.py` — sidecar write after rectified BAM
- `rectify/core/commands/run/single_sample.py` — plumb `trust_existing_bams`

**Remote only (Hoffman2 — NOT in git):**
- `/u/project/guillom/shared/processed/han_2023_ysh1aa_vs_wt_rectify_v0.9.0/manifest.tsv`
- `.../han_2023_.../qsub_run.sh`
- `.../han_2023_.../README.md`
- `.../han_2023_.../results/` (empty, output will land here)

**Remote only (Sherlock — NOT in git):**
- `set2_cpa_machinery/dry_run_20260518/manifest.tsv`
- `.../dry_run_20260518/symlink_layout.sh`
- `.../dry_run_20260518/run_all.cmd`
- `.../dry_run_20260518/README.md` [updated this session — added `--skip-alignment` layout caveat]
- `.../dry_run_20260518/results/` (empty)

**[uncommitted] — pre-existing WIP, not part of provenance feature:**
- `HANDOFF.md` (this file — modified)
- `CLAUDE.md`, `README.md` — docs WIP
- `rectify/core/commands/correct_command.py` — other workstream
- `rectify/core/splice/junction_refiner.py` — other workstream
- `rectify/data/*`, `scripts/*`, `tests/test_junction_refiner.py`,
  `docs/figures/*`, `scripts/calibration/*`, `scripts/validation_data/*`
  — see `git status -s` for full list
