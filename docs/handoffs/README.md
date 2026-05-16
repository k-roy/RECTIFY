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

1. `bam_processor_split.md` — adjacent to the walkback follow-ups
2. `cdna_correct_command_split.md` — adjacent to cDNA pipeline work
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
- After each carve, run the full broad sweep:
  `python -m pytest tests/ --no-header --no-cov -q --ignore=tests/test_cdna_correct.py --ignore=tests/test_cdna_chain_canary.py`
- For commits, follow the conventional-commits style already used
  in this repo: `refactor(<scope>): <one-line summary>`.
- Commit messages should include the
  `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`
  trailer.
