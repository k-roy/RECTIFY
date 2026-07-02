# cDNA agent — scope & handoff (CI health + cDNA validation set)

**Created:** 2026-07-01 · by the DRS-validation-review session (branch `drs-validation-rebuild`, PR
k-roy/RECTIFY#3). **Owner:** the future cDNA agent. **Urgency:** LOW — RECTIFY is pre-release 0.9.0; there is
**no rush to green PR #3**. The governing principle: **fix more bugs than we introduce.** Carry the regression
discipline the DRS walkback work established (golden-hash test, clip-to-win gate, two-strand + before/after
audits, advisor cross-checks) into every cDNA change.

## Why this exists

The DRS-validation session touched only DRS code (`rectify/core/correct/walkback.py` + the DRS validation
bundle). When PR #3's CI ran for the FIRST time on this branch, it surfaced **three pre-existing** failures —
none caused by the DRS work or the doc merge. Two are trivial infra fixes (done/ready); the third is a
**cDNA-arm** issue that belongs with you, especially since you'll be building the cDNA validation read set
next. Bundle all of this into that effort rather than hot-patching CI now.

## CI failures on PR #3 (run 28497523365, HeadSHA 99cc38f) — status

1. **samtools missing in CI** — ✅ FIXED + pushed (`99cc38f`): added an apt `samtools` install step to
   `.github/workflows/tests.yml`. The bundled validation tests run live `rectify correct --aligner-bams`
   (Module 2H sorts/indexes BAMs) which needs the samtools binary.

2. **Python 3.8 collection errors** — ✅ FIXED locally in commit **`959eff8` (UNPUSHED, see risk below)**.
   `rectify/config.py` + `rectify/core/analyze/manifest.py` used PEP 585 builtin-generic annotations
   (`dict[...]`, `list[...]`, `tuple[...]`) without `from __future__ import annotations`; on 3.8 these
   evaluate at import and raise `TypeError: 'type' object is not subscriptable`. config.py is imported nearly
   everywhere → 29 collection errors. Fix = add the future-import to those two files (the other 4 branch
   files with that syntax already had it). Verified: no other 3.9+ runtime features (removeprefix,
   dict-union, match/case, graphlib) anywhere in `rectify/`. `pyproject` keeps `requires-python = ">=3.8"`
   (this was option A, restoring parity with the 3.8 classifier + CI matrix + master).

3. **`tests/test_validation_reads_cdna.py::TestCdnaWalkbackBaseline::test_anchor_and_tail` — ❌ YOUR TASK.**
   Env-sensitive baseline. 3 reads (`cdna_cat1_plus_2`, `cdna_cat1_minus_2`, `cdna_cat4_plus_2`) fail:
   **anchor matches, tail_len differs** (CI 41/35/56 vs baseline 23/18/36). It **PASSES locally**
   (numpy 1.23.4 / pysam 0.23.3) and **FAILS in CI** (numpy 2.0.2 / pysam 0.24.0). Pre-existing test
   (baselines locked 2026-05-16, commit c411c80). The tested code is `rectify/core/cdna/walkback.py::
   walk_back_anchor_and_tail` (the cDNA arm — untouched by the DRS session). Ruled out: edlib (test passes
   with/without it); the committed BAM + genome are byte-identical local vs CI. So the differ is **pysam 0.24
   or numpy 2** behavior vs the 2026-05-16 baseline-capture env.

### What to investigate for #3
- Pin down whether it's **pysam 0.24** (BAM/CIGAR/`query_sequence` parsing) or **numpy 2.0** — reproduce by
  installing those exact versions in a throwaway env and running the 3 reads through `walk_back_anchor_and_tail`;
  diff the intermediate `anchor_qp`, `seq` slice, and `_find_adapter_anchor_pos` result vs the old libs.
- Decide which tail_len is **biologically correct** (this is the real question — same discipline as the DRS
  cat2_plus_1 investigation: is the CI value right and the baseline stale, or vice versa?). `tail_len =
  tail_seg.count(tail_base)` where `tail_seg` is bounded by `anchor_qp` and the SSP/adapter anchor.
- Fix options (pick per your finding, not blindly): (a) if CI libs give the correct answer → **re-capture the
  baselines** via `scripts/validation_data/characterize_baseline.py --arm cdna` under the pinned CI libs, and
  make the capture env explicit; (b) if the old libs are correct → **pin `pysam<0.24` / `numpy<2`** in
  `pyproject`/CI; (c) make the test **tolerant** (assert anchor exactly, tail_len within a documented band) if
  the tail is legitimately lib-sensitive at the margin; (d) last resort: mark `slow` to keep it out of fast CI
  (loses coverage — avoid unless justified).
- **Regression gate:** whatever you choose, the 3 reads must become a *stable* baseline across the supported
  lib matrix, and you must confirm the fix doesn't silently change other cDNA reads. Add the reasoning to a
  comment on `EXPECTED_ANCHOR_AND_TAIL` (like the DRS upf1d cat1_minus_1 precedent).

## The cDNA validation read set (the larger deliverable)

Build the cDNA analogue of the DRS validation bundle once #3 is settled:
- **Model on the DRS bundle:** `rectify/data/validation/` (cat1–cat9 reads, per-aligner BAMs, `corrected_reads.tsv`)
  + `rectify/data/validation/CASE_STUDIES.md` (per-read mechanistic write-ups) +
  `scripts/validation_data/` (renderer + report generator).
- cDNA data already present: `rectify/data/validation/validation_reads_cdna.bam` (+ the
  `test_validation_reads_cdna.py` baselines). The cDNA pipeline is the 3-command sequence
  `correct-cdna → align → cdna-analyze` (NOT `run-all` — verified there is no `--protocol` flag; see
  `docs/quickstart_cdna.md`). Chain canary: `tests/test_cdna_chain_canary.py`.
- Produce a `CASE_STUDIES`-style mechanistic doc for the cDNA reads (UMI dedup-before-align, cross-orient
  merge with `--max-cross-orient-span`, `--per-cluster-cap`), and wire the review renderer for the cDNA arm.

## ⚠️ Risk on the unpushed py38 fix (`959eff8`)
It is committed **locally only** on `drs-validation-rebuild` (origin is at `99cc38f`). A local-only commit on a
SHARED branch (a COMPASS session also commits here) can be lost to a force-push/reset. When you take this up,
either **push it** (`git push origin HEAD:drs-validation-rebuild` after an FF-safety check that
`origin == 99cc38f`) or **cherry-pick it** onto your working branch so it isn't lost. It won't green CI alone
(#3 still fails), so it's fine to bundle with the #3 fix.

## Pointers
- CI config: `.github/workflows/tests.yml` (samtools step added; matrix `3.8`–`3.12`).
- DRS regression precedents to emulate: golden-hash `tests/test_bam_parallel_state.py`; the DRS walkback fix
  `fc44ee2`; the graded-clip audit `scripts/benchmark/graded_clip_audit.py` + `run_graded_clip_audit*.sh`
  (the Option-B "shelved: re-breaks clip-to-win" negative result is in `CASE_STUDIES.md` cat2_plus_1 — a
  cautionary tale on strand handling + calibrating on ambiguous reads).
- Full session state: root `HANDOFF.md`.
