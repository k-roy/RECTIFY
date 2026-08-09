# HANDOFF — Phase 0 (WIP landing + master sync) → Phase 1 (Rectify Re-aligner) — CURRENT 2026-08-09

**Branch:** `chore/vendor-desalt-chanfreau1` (main checkout). Kevin approved: land all dirty WIP,
fast-forward master, then build the Re-aligner consensus-triage layer from the validated
native-aligner work. Decision record: `dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md`.
Program state of record: `dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md` (restored copy — see Files).

## Done (this session, all commits on `chore/vendor-desalt-chanfreau1`)

- **`70e9664` fix(run-all)** — the three Chanfreau-630 manifest-era bugs (sample column into
  region TSVs; `_readable_corrected_tsv()`; empty-cluster schema; `_is_manifest` superset check).
- **`5792b22` feat(qc,browser-pack)** — 630's Rbrowse unit as run-all's fail-soft final stage.
- **`3d893ef` feat(ont-cdna)** — Path A (UMI-collapse to molecules) default, both entry points;
  chunked-alignment protocol-flag refusal; UMI neighbour-sets/star-split.
- **`b9ae0e1` perf(splice)** — `_hp_edit_distance` numba fast path (memory-gated OFF) + cutoff.
- **`7648725` feat(analyze)** — 377's "group 3": containment-first gene attribution (NEW DEFAULT),
  region_class, transcript model, ncRNA atlases. Pending since 2026-07-21.
- **`999ceb5` docs** — accumulated dev record (native-aligner docs, COMPASS/GMAP, reassessment).
- Inbox: all 6 messages actioned + archived to `.claude/inbox/.read/`. Coordination note sent to
  the Chanfreau workspace inbox (630/653/641/377) with commit hashes.

## Verified

- Every commit's staged .py content parse-checked; the three mixed files (single_sample /
  multi_sample / stages) were split hunk-level across units A/B/C and reassemble **byte-identical**
  to the pre-commit working tree (zero unstaged residue after the splits).
- 630 had verified the identical tree content end-to-end (full suite `1986 passed`, run-all
  single+multi green) before handing off; my own clean-tree suite run is the in-flight gate below.

## ⚠ Incidents / concurrent-agent facts (2026-08-09 afternoon)

1. **`.claude/worktrees/` (all 11 agent worktrees) was DELETED mid-session by an external cleanup**
   — not this session. All branches/commits survive (shared object store). Known uncommitted
   casualty: the benchmark worktree's state audit, restored verbatim as
   `dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md` (provenance note at top). If other worktrees held
   uncommitted work, it is gone — flag to Kevin.
2. **A live Chanfreau session (653) is editing this checkout** (`cdna/cluster.py`,
   `cdna/consensus.py`, `tests/test_cdna_consensus_653.py`, `tests/test_umi_clustering_cdist.py`
   — left uncommitted, theirs to land). Do not touch; do not `git checkout` branches in this tree.
3. **Shell gotcha (cost an hour):** this environment's `grep`→ugrep aliasing + git pipes returned
   FALSE ZEROS (`git grep`/`git diff | wc -l` claimed the worktree branch lacked motif_blind — it
   does not). For load-bearing checks, dump blobs to files (`git show ref:path > f`) and use
   `/usr/bin/grep` / `cmp` / `diff` on the files.

## Open / in flight

- **[IN FLIGHT] Clean-tree full suite** (background task `bzgwtfyf1`) on detached scratch worktree
  at `999ceb5` (excludes 653's live WIP):
  `<scratchpad>/suite_check`, log `<scratchpad>/../tasks/bzgwtfyf1.output`.
- **[BLOCKED on suite] master fast-forward** — ref-update only (master is checked out nowhere):
  `git fetch . chore/vendor-desalt-chanfreau1:master` (refuses non-ff by construction).
- Still-open bugs (filed, NOT this session's scope): AGENT_FIXES.md 2026-07-21 CRITICAL
  single-aligner path race (`align_command.py`); 2026-08-07 junction-pool-density cost scaling;
  Oak editable-install drift (Sherlock env runs uncommitted vintage).
- `feat/overhang-resolver-641`: 641/643's live branch (tip `b028e35`, 2 commits behind our line) —
  they land it themselves after cluster acceptance T3–T8; the triage layer adopts it as its
  overhang leg.

## Resume (concrete)

1. `pgrep -f suite_check || cat <scratchpad>/../tasks/bzgwtfyf1.output`
   - **PASS** (`… passed`, 0 failed) → `git fetch . chore/vendor-desalt-chanfreau1:master` →
     `git log -1 --oneline master` confirms → commit this HANDOFF refresh + audit-doc restore →
     start Phase 1 (next item).
   - **FAIL** → read the failing test; if it is in qc/browser/analyze/ont-cdna units, fix forward
     on this branch (the units are separable commits — worst case `git revert <unit>`); re-run.
     Do NOT ff master on a red suite.
2. **Phase 1.1 (port = MERGE):** `git worktree add <path> -b feat/realigner-triage master` (fresh
   worktree — NEVER switch branches in the main checkout, 653 is live) → `git merge 2b8d2ed`
   (= `worktree-agent-a25a2c1e784ad37dc`, the validated native re-aligner tip; Kevin already
   approved this merge once in July — it never executed). Expected conflicts are small and mapped:
   - `hp_penalty.py`: keep HEAD's `_hp_run_length` bounds-guard (1821d4d) + take wt's rate-table
     additions; the `from_tsv` return line combines both kwarg sets.
   - `junction_scoring.py`: wt side is additive (concat-DP `_USE_CONCAT_DP` default ON, flat-cost
     constants, full-run/refcol ins flags default OFF).
   - `local_aligner.py`: wt adds C1 knob (`penalty_table/lam/ins_lengthlaw`, None ⇒ byte-identical).
   - `junction_refiner.py` + `test_junction_refiner.py`: purely additive (motif_blind, guards
     dormant at 0.0, compensating-indel invariant e40ca00 always-on).
   Then full suite + the branch's fence suites (`test_junction_refiner`, `test_hp_drift_guard`,
   `test_microhom_drift_guard`, `test_c1_lengthlaw`) must be green; byte-identical-off fences are
   the acceptance bar. Port prep diffs: `<scratchpad>/port/port__*.diff`.
3. **Phase 1.2:** the `scoring.py::_count_junction_proximity_errors` post-N `prev_rp` surgical fix
   (audit §5.1) behind `discovery_tiebreak_probe.py` + `smoke_roundtrip.py`.
4. **Phase 1.3:** triage layer per `dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md` §3 (bypass /
   re-align / re-entry; pool-level discovery gate separate from read-level accuracy gate).

`<scratchpad>` = `/private/tmp/claude-501/-Users-kevinroy-work-rectify/12349562-4dab-426f-bb0c-0f5ca6d90c91/scratchpad`

## Files

- Decision + design: `dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md`
- Program state of record (restored): `dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md`
- Native-aligner sources: branch `worktree-agent-a25a2c1e784ad37dc` @ `2b8d2ed`; key commits
  `69a230f` (motif_blind), `e40ca00` (compensating-indel invariant), `dd257b8`/`e1ed90c`
  (concat-DP), guard-shelve verdict `d5b25d3`, Phase-6 overview `fc58950`.
- Prior CMA-session handoff: superseded; see `dev/CMA_PROGRESS.md` + git history of this file.
