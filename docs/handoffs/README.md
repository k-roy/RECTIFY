# Handoff Briefs

Scoped specs for the remaining single-file refactors after the
2026-05-16 package reorganization. Each brief describes a "split this
giant file into a subpackage" task that is too large and judgment-heavy
to be batched with the rest of the reorg.

## How to use

Each brief is self-contained — it describes the target file, the
proposed split, the migration plan, and the test gate. Pick up one
brief per CLI session. Do not run two splits in parallel: most of
these files import each other, and concurrent edits to call paths
will conflict.

Recommended order (most-actively-edited first):

1. ~~`bam_processor_split.md`~~ — **DONE 2026-05-16** (commit `0024fa3`)
2. ~~`cdna_correct_command_split.md`~~ — **DONE 2026-05-16** (commit `c411c80`)
3. `analyze_command_split.md` — `analyze/` subpackage already exists,
   so this is the cleanest carve
4. `bam_writer_split.md`
5. `consensus_split.md`
6. `junction_refiner_split.md`
7. `run_command_split.md` — the largest; do last, after the others
   are stable so the new structure is settled

## Conventions across briefs

- All file paths assume the post-reorg layout (commands/, bam/, etc.).
- Tests are at `tests/`; do not move them, just update imports.
- After each carve, run the broad sweep, fast variant:
  ```
  python -m pytest tests/ -m "not slow" -n auto
  ```
  This deselects the 4 chrI-pipeline smoke tests (5–10 min each) and
  parallelizes the rest across CPU cores. The full sweep used to take
  ~8 min sequentially; this drops it to ~1–2 min.
- For a final pre-merge gate that also runs the slow pipeline tests:
  ```
  python -m pytest tests/ -n auto
  ```
- Coverage is opt-in. For a coverage report:
  ```
  python -m pytest tests/ --cov=rectify --cov-report=html --cov-report=term
  ```
- For commits, follow the conventional-commits style already used
  in this repo: `refactor(<scope>): <one-line summary>`.
- Commit messages should include the
  `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`
  trailer.
