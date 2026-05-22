# Documentation Consolidation Proposal — 2026-05-22

**Purpose:** Pre-publication cleanup of scattered MD docs in the RECTIFY repo.  
**Branch:** `drs-validation-rebuild`  
**Author:** Claude Code (review session 2026-05-22)  
**Status:** PROPOSAL — no files moved yet. Awaiting Kevin's approval.

---

## Scope and constraints

**Package boundary (pyproject.toml):** Only the `rectify/` Python subpackages + specific data files ship to PyPI. Root MDs, `docs/`, `dev/`, `scripts/`, `handoffs/` do NOT ship to package users. The clutter is entirely a GitHub / repo-visitor problem.

**Docs site boundary (mkdocs.yml nav):** The published ReadTheDocs nav is well-curated. Files in `docs/` but NOT in the nav are invisible to doc readers but still clutter the source tree.

**Active WIP — do not touch without Kevin's explicit sign-off:**
- `AGENT_FIXES.md` — uncommitted changes (modified)
- `HANDOFF.md` — uncommitted changes (modified)
- `README_KR_edits.md` — uncommitted changes (modified)

---

## Summary by bucket

| Bucket | Count | Action |
|--------|-------|--------|
| Public / external — keep as-is | 6 | Nothing to do |
| Root clutter — move to `dev/` | 6 | `git mv` |
| Root clutter — delete | 1 | `git rm` |
| Needs Kevin decision before action | 3 | Flag only |
| `docs/` orphans — move to `dev/handoffs/` | 15 | `git mv` |
| Figure scripts — move to `scripts/figures/` | 3 | `git mv` |
| `dev/` — already correct | all | Nothing to do |
| `handoffs/` root dir — consolidate | 25+ files | `git mv` into `dev/handoffs/` |

---

## Bucket 1: Public / external — keep at root, no changes

| File | Rationale |
|------|-----------|
| `README.md` | Public landing page; already polished |
| `CHANGELOG.md` | Keep a changelog format; referenced by `docs/changelog.md` snippet |
| `CONTRIBUTING.md` | Standard OSS file |
| `CODE_OF_CONDUCT.md` | Standard OSS file |
| `SECURITY.md` | Standard OSS file |
| `CITATION.cff` | Machine-readable citation; required for software papers |
| `CLAUDE.md` | Claude Code convention; increasingly common in OSS; content is dev guidance, not secrets. Keep. |
| `AGENTS.md` | Multi-agent startup doc; analogous to CLAUDE.md. Keep. |

---

## Bucket 2: Root clutter — move to `dev/`

These files are tracked in git, are internal/operational, and have no reason to be at the repo root. Moving them to `dev/` cleans the root without losing history.

| File | Move to | Rationale |
|------|---------|-----------|
| `TODO.md` | `dev/TODO.md` | Internal development backlog; not user-facing |
| `CODEX_AUDIT.md` | `dev/audits/CODEX_AUDIT_20260520.md` | Internal multi-agent audit; 422 lines of implementation notes |
| `PLOT_SKILLS.md` | `dev/PLOT_SKILLS.md` | Claude agent context for `rectify.visualize`; not user documentation |
| `pA_tail_DRS_citation_validation.md` | `dev/reviews/pA_tail_DRS_citation_validation.md` | Working analysis doc for the executive summary HTML; belongs in dev reviews |
| `RECTIFY_SHERLOCK_HANDOFF_20260518.md` | `dev/handoffs/RECTIFY_SHERLOCK_HANDOFF_20260518.md` | Old cluster handoff |
| `RECTIFY_SHERLOCK_HANDOFF_20260518_v2.md` | `dev/handoffs/RECTIFY_SHERLOCK_HANDOFF_20260518_v2.md` | Old cluster handoff |

**Commands (safe to run when ready):**
```bash
git mv TODO.md dev/TODO.md
git mv CODEX_AUDIT.md dev/audits/CODEX_AUDIT_20260520.md
git mv PLOT_SKILLS.md dev/PLOT_SKILLS.md
git mv pA_tail_DRS_citation_validation.md dev/reviews/pA_tail_DRS_citation_validation.md
git mv RECTIFY_SHERLOCK_HANDOFF_20260518.md dev/handoffs/RECTIFY_SHERLOCK_HANDOFF_20260518.md
git mv RECTIFY_SHERLOCK_HANDOFF_20260518_v2.md dev/handoffs/RECTIFY_SHERLOCK_HANDOFF_20260518_v2.md
```

Also update `CLAUDE.md`'s reference table to point `PLOT_SKILLS.md` → `dev/PLOT_SKILLS.md`.

---

## Bucket 3: Delete

| File | Rationale |
|------|-----------|
| `CODEX.md` | 8 lines that just say "read AGENTS.md first" and redirect to CLAUDE.md. Pure redundancy. |

```bash
git rm CODEX.md
```

---

## Bucket 4: Needs Kevin decision — flag only, do not touch

| File | Status | Question |
|------|---------|----------|
| `AGENT_FIXES.md` | **Uncommitted changes (modified)** | After committing/stashing WIP: move to `dev/AGENT_FIXES.md`. This is a 1,195-line internal coordination log — the largest non-public file at root. Key entries should be migrated to CHANGELOG.md as history. |
| `HANDOFF.md` | **Uncommitted changes (modified)** | After committing/stashing WIP: move to `dev/handoffs/HANDOFF_CURRENT.md`. Active session handoffs don't belong at root for a published repo. |
| `README_KR_edits.md` | **Uncommitted changes (modified)** | Which version is canonical — `README.md` or `README_KR_edits.md`? Once decided: reconcile the two (the KR version has better wording in several places), then `git rm` the draft. Don't touch until Kevin resolves. |

---

## Bucket 5: `docs/handoffs/` — move all to `dev/handoffs/`

14 handoff files live in `docs/handoffs/`. None appear in the mkdocs.yml nav. They are internal AI-agent session notes that pollute the `docs/` namespace (which GitHub visitors browse). Moving them to `dev/handoffs/` consolidates all handoffs in one place.

**Files:**
```
docs/handoffs/HANDOFF_2026-05-19_drs_validation_queue_audit.md
docs/handoffs/HANDOFF_2026-05-19_validation_qsrev_cdna_pass.md
docs/handoffs/HANDOFF_2026-05-20_corrected_consensus_nan_junctions_bug.md
docs/handoffs/HANDOFF_2026-05-20_gtf_feature_expansion.md
docs/handoffs/HANDOFF_2026-05-20_parallel_bam_writer_abort.md
docs/handoffs/HANDOFF_2026-05-20_parallel_checkpoint_resume_hardening.md
docs/handoffs/HANDOFF_2026-05-20_set2_minimap2_only_salvage.md
docs/handoffs/HANDOFF_2026-05-20_sherlock_penaltytableset_resync.md
docs/handoffs/anchored_prefix_calibration.md
docs/handoffs/cdna_umi_phase_d_resume.md
docs/handoffs/cdna_umi_stratified_calibration.md
docs/handoffs/debugger_queue.md
docs/handoffs/empirical_error_tables_cdna_qsrev.md
docs/handoffs/regression_resolution.md
```

Note: `docs/handoffs/HANDOFF_2026-05-20_sidecar_architecture.md` is **untracked** — stage it to `dev/handoffs/` before running the batch mv.

**Command:**
```bash
mkdir -p dev/handoffs
git mv docs/handoffs/*.md dev/handoffs/
# Then manually: git add dev/handoffs/HANDOFF_2026-05-20_sidecar_architecture.md
```

---

## Bucket 6: Consolidate `handoffs/` root directory into `dev/handoffs/`

The root `handoffs/` directory (25+ files across `handoffs/` and `handoffs/_archive/`) duplicates the function of `dev/handoffs/`. Consolidate:

```bash
git mv handoffs/HANDOFF_2026-05-18_*.md dev/handoffs/
git mv handoffs/HANDOFF_2026-05-19_*.md dev/handoffs/
git mv handoffs/HANDOFF_2026-05-20_*.md dev/handoffs/
mkdir -p dev/handoffs/archive
git mv handoffs/_archive/* dev/handoffs/archive/
# Then remove the now-empty handoffs/ dir (git does this automatically)
```

After this, all session handoffs live in exactly one place: `dev/handoffs/`.

---

## Bucket 7: Figure generation scripts — move to `scripts/figures/`

Three figure-generation scripts are tracked at the repo root:

| File | Move to |
|------|---------|
| `generate_figures.py` | `scripts/figures/generate_figures.py` |
| `generate_figures_v2.py` | `scripts/figures/generate_figures_v2.py` |
| `generate_cdna_figures.py` | `scripts/figures/generate_cdna_figures.py` |

```bash
mkdir -p scripts/figures
git mv generate_figures.py scripts/figures/
git mv generate_figures_v2.py scripts/figures/
git mv generate_cdna_figures.py scripts/figures/
```

Update `dev/figures/PROVENANCE.md` and `STYLE_GUIDE.md` references if they point to root-level scripts.

---

## Bucket 8: `docs/audit_history.md` — decide on fate

This file (951 lines) is not in the mkdocs.yml nav. It contains:
1. A 2026-05-21 documentation audit of 9 issues (H1–H3, M1–M4, L1–L2), most of which were fixed in commit `2b0b304`.
2. A pre-public version history section covering internal v2.x/v3.x development.

**Options:**
- A) **Add to mkdocs.yml nav** as a "Release History" page — gives it a home and makes it discoverable.
- B) **Move to `dev/audits/`** — treat it as internal history, not user-facing.
- C) **Split**: extract the coordinate convention statement (which is referenced by `dev/audits/bedgraph_coordinate_audit_20260520.md`) into `docs/coordinate_system.md` (already exists in nav), and move the rest to `dev/`.

Recommendation: **Option B** — move to `dev/audits/audit_history.md`. The pre-public version history mentioning "v2.x / v3.x — internal" is potentially confusing to external readers; CLAUDE.md already flags these as "internal pre-public versions" that should not be surfaced as canonical. The coordinate audit material is better placed in `docs/coordinate_system.md`.

---

## Bucket 9: `.github/ISSUE_walkback_integration.md`

This is a design-doc-style GitHub issue write-up, not a typical `.github/` file. Its acceptance criteria are all checked off (implemented and merged). 

**Options:**
- A) Convert to a proper GitHub issue (or verify it was already filed as one) and delete the file.
- B) Move to `dev/specs/walkback_integration_retrospective.md` as historical context.

Recommendation: **Option B** (safe, no GitHub action required). The document contains valuable algorithm description useful for future contributors.

---

## After all moves: update `CLAUDE.md` reference table

The `CLAUDE.md` READ FIRST table references `PLOT_SKILLS.md` and `HANDOFF.md` at their current root locations. Once the moves above are approved and executed, update the table entries to point to the new paths.

---

## After all moves: expected root-level MD files

```
AGENTS.md          ← AI agent startup (keep)
AGENT_FIXES.md     ← pending move to dev/ after WIP committed
CHANGELOG.md       ← public changelog
CITATION.cff       ← machine-readable citation
CLAUDE.md          ← Claude Code context (keep)
CODE_OF_CONDUCT.md ← standard OSS
CONTRIBUTING.md    ← standard OSS
HANDOFF.md         ← pending move to dev/handoffs/ after WIP committed
LICENSE            ← required
README.md          ← public landing page (resolve vs README_KR_edits.md first)
SECURITY.md        ← standard OSS
```

That is 11 files — down from 18 — with everything extraneous moved to `dev/` where it belongs.

---

## Prioritization

Execute in this order to minimize risk and get the biggest visual win first:

1. **Batch 1 (safe, no WIP conflict):** Bucket 2 moves + Bucket 3 delete + Bucket 5 + Bucket 6 + Bucket 7  
2. **Batch 2 (Kevin decides):** `README_KR_edits.md` reconciliation  
3. **Batch 3 (after committing WIP):** `AGENT_FIXES.md` and `HANDOFF.md` moves  
4. **Batch 4 (optional):** `docs/audit_history.md` fate, `.github/` issue  
