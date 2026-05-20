# HANDOFF — Re-sync Sherlock PenaltyTableSet WIP from M1

**Date:** 2026-05-20
**Branch:** `drs-validation-rebuild`
**Discovered by:** another agent's session note (quoted in user message 2026-05-20)
**Estimated effort:** 30-60 min (depends on what's still uncommitted on M1 + on lab data-deposit conventions)
**Owner:** unassigned — Kevin should ideally drive this since it touches cluster sync + data deposit decisions

---

## Background

Another agent observed on Sherlock that a 177-line uncommitted diff existed which added:

1. A new `PenaltyTableSet` class — maps each read to a per-UMI-cluster-depth `HpPenaltyTable` based on the read's `XC` tag (UMI count). Bins: umi1, umi2, umi3plus, plus a pooled fallback for reads with no `XC` tag.
2. In `correct_command.py`: when `--dT-primed-cDNA` (cDNA mode) is set, looks up per-UMI-bin penalty TSVs from `rectify/data/genomes/<organism>/penalty_tables/` and instantiates a `PenaltyTableSet`.

That agent reported the diff was lost on Sherlock between sessions (presumably stashed/checkout-away during 7dbb1bd → cb2fe6c commits).

## Current state on M1 (verified 2026-05-20)

**`PenaltyTableSet` IS implemented and tested on M1 at HEAD `cb2fe6c`.** Specifically:

| Component | Location on M1 | git state |
|---|---|---|
| `PenaltyTableSet` class | `rectify/core/splice/junction_refiner.py:192` | committed (in HEAD) |
| `correct_command.py` cDNA wiring | `rectify/core/commands/correct_command.py:534-576` | committed (in HEAD) |
| `TestPenaltyTableSet` (7 tests) | `tests/test_junction_refiner.py:833-1031` | **uncommitted (modified)** |
| Per-UMI penalty TSVs (`penalty_scores_cdna_umi{1,2,3plus}.tsv`) | `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/` | **uncommitted (untracked)** |
| `penalty_scores_cdna.tsv` (pooled), `penalty_scores_qsrev.tsv` | same dir | **uncommitted (untracked)** |

All 7 `TestPenaltyTableSet` tests **pass** on M1: `pytest tests/test_junction_refiner.py::TestPenaltyTableSet` → 7 passed in 3.52s.

So the canonical PenaltyTableSet code is already committed at `cb2fe6c`. The tests and data files just need to be committed and pulled to Sherlock.

## Likely explanation of the Sherlock-side loss

The other agent's "177-line diff" was probably an **older draft** of `PenaltyTableSet` from a different session, which Kevin restructured and re-committed as `cb2fe6c`. When `cb2fe6c` was checked out on Sherlock (or when Kevin stashed Sherlock-side WIP to pull `cb2fe6c`), the older draft was discarded — superseded, not lost.

Verify by diffing what Sherlock had against what M1's `junction_refiner.py:192` has now (see "Step 2" below).

## What needs to happen

### Step 1 — Commit the M1 untracked / modified state

Commit just the PenaltyTableSet-related additions, surgically (per `CLAUDE.md` "Surgical staging" rule — no `git add -A`):

```bash
cd /Users/kevinroy/work/rectify

# 1a. The test class
git add tests/test_junction_refiner.py

# 1b. The per-UMI penalty TSVs (4 files for cerevisiae)
git add rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv
git add rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna_umi1.tsv
git add rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna_umi2.tsv
git add rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna_umi3plus.tsv
git add rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_qsrev.tsv

# 1c. Verify nothing else got staged
git status -s

# 1d. Verify tests still pass before committing
pytest tests/test_junction_refiner.py::TestPenaltyTableSet -v

# 1e. Commit
git commit -m "feat(cdna): commit PenaltyTableSet tests + per-UMI penalty TSVs

The PenaltyTableSet class was committed in cb2fe6c but its 7 unit tests
and the per-UMI penalty data files lived locally untracked. This commit
moves them into the repo so Sherlock + H2 can pull and run the cDNA
pipeline path that consumes the per-UMI tables."
```

### Step 2 — Push and pull on the clusters

```bash
# On M1
git push origin drs-validation-rebuild

# On Sherlock (Kevin must run; ControlMaster auth required)
ssh sherlock 'bash --norc --noprofile -c "
  cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify
  git status -s
  git stash save 'pre-pull-pentablset-resync'  # safety
  git fetch origin
  git merge origin/drs-validation-rebuild --ff-only
  git log --oneline -3
"'

# On H2 (no 2FA — agent can do this)
ssh hoffman2 'bash -lc "
  cd /u/home/k/kevinroy/software/rectify
  git fetch origin
  git merge origin/drs-validation-rebuild --ff-only
  git log --oneline -3
"'
```

After the pull, the agent on Sherlock that originally complained about a "discarded 177-line diff" should re-inspect HEAD's `junction_refiner.py:192` and confirm the M1 version is functionally equivalent to what they had. If equivalent: no action needed (the work wasn't lost, just re-committed). If non-equivalent: see Step 3.

### Step 3 (only if needed) — Recover the Sherlock-side draft

If Step 2 reveals real divergence (e.g., Sherlock's draft had a feature the M1 commit doesn't), recover from Sherlock-side stash or reflog:

```bash
ssh sherlock 'bash --norc --noprofile -c "
  cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify
  git stash list
  git reflog | head -20
"'
```

Look for an entry near the 7dbb1bd/cb2fe6c transition. If found, cherry-pick to a branch and diff against M1's HEAD. Surface any real differences to Kevin for review before merging.

## Acceptance

- [ ] M1 HEAD includes the test file + the 5 untracked penalty TSVs.
- [ ] `pytest tests/test_junction_refiner.py::TestPenaltyTableSet -v` passes (7 tests) on M1 after commit.
- [ ] Sherlock HEAD and H2 HEAD match M1 HEAD (`git log --oneline -1` returns same sha on all three).
- [ ] Sherlock-side agent confirms no real divergence between M1 PenaltyTableSet and their lost draft, OR documents any features that need porting forward.
- [ ] If divergence found, a Step-3 cherry-pick + Kevin review is recorded in a follow-up handoff.

## Dispatch options

- **Kevin runs Steps 1-2** directly — fastest, ~10 min, since Step 2's Sherlock side requires Kevin's 2FA Kerberos auth that an agent can't drive.
- **Agent drives Step 1** (M1 commit) + Step 2's H2 sync; Kevin manually does Step 2's Sherlock sync. Splits the work.
- **Bundle with the Commit A landing** so the cluster sync is a single pull rather than two — recommended sequencing.

## Why this is in its own commit (or paired with Commit A)

These files are not part of Commit A's parallelism refactor. They're an independent feature (cDNA UMI-stratified penalty tables) whose code shipped at `cb2fe6c` but whose tests + data didn't. Keep this commit (or set of commits) separate from the parallelism work for clean review and clean blame.
