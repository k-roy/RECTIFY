# validation_read_review — status

This folder holds the rendered per-read review artifacts for the RECTIFY
validation bundle. The canonical operating reference lives **outside this
folder**:

> **`rectify/scripts/validation_data/PLOTTING.md`** — current renderer
> behavior, color conventions, CLI, bundle paths, toggles. Read first
> before touching anything visualization-related.

## What's in this folder

| Path | Role |
| --- | --- |
| `cat{1..9}_*_review.md` | Per-category review markdown — one per category, one read section per file. Regenerated from `generate_review_report.py --per-category`. |
| `cat{1..9}_*_review_pngs/` | The per-read PNGs (one per read) that the `.md` files embed. |
| `cat{1..5}_*_findings.md` | Plotter→debugger findings docs (open + resolved items, organized by category). These are the canonical channel for surfacing pipeline-side issues spotted during review. |
| `bundle_dorado_source_stale_findings.md` | RESOLVED. Kept as historical record (debugger fixed via `c3202f6` on Sherlock; M1 dorado-BAM scp'd in sync). |
| `polya_trim_audit.tsv` | 36-row audit output from `scripts/validation_data/audit_polya_trim.py`. Confirms the "trim to first non-A" invariant holds across the validation set. |
| `_archive/` | Previous session handoff docs (`PLOTTER_HANDOFF_2026-05-18*.md`). Useful for layout-decision history but no longer load-bearing — PLOTTING.md is the current source of truth. |

## Regenerating

```bash
cd /Users/kevinroy/work/rectify
# Full re-render of all 9 categories (~60–70 s for 36 reads on M1):
.venv/bin/python scripts/validation_data/generate_review_report.py \
    --arm drs --per-category --format md
# Single category for fast iteration:
.venv/bin/python scripts/validation_data/generate_review_report.py \
    --arm drs --per-category --category cat3_junction --format md
```

Output lands here under `cat{N}_*_review.md` + `cat{N}_*_review_pngs/`.

## When to update the findings docs

When reviewing a category and you notice a pipeline-side anomaly (a CIGAR
that doesn't match what the TSV claims, a missing N-op, an unexpected
soft-clip pattern, etc.) — log it in the matching `cat{N}_*_findings.md`
so the debugger session sees it. Don't edit `rectify/core/**` from the
plotter side; the role split is documented in PLOTTING.md.

## Last refresh

- PNGs regenerated 2026-05-19 (renderer commit pending). All 36 reads
  rendered cleanly; no glyph warnings.
- `polya_trim_audit.tsv` generated 2026-05-19 against the
  Sherlock-rebuilt dorado-source BAM (`c3202f6`). 36/36 trim_applied,
  invariant flag B = 0.
