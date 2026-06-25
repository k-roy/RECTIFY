# Handoff — TODO/BUGS cleanup sweep (session 2026-06-24/25)

**Date:** 2026-06-25
**Branch:** `drs-validation-rebuild` (26 commits ahead of origin at session start; this
session's work is UNCOMMITTED in the working tree — user has not asked to commit).
**Task:** "tackle all remaining TODOs in the rectify package." User chose all buckets;
cluster reachable; COMPASS short-read mode OUT of scope. **9 of 10 tracked items resolved;
GLASS aligner eval is the only one deferred.**

Recurring theme: most TODO/BUGS entries were **stale** — already fixed by the 26 unpushed
commits. Each was VERIFIED against current code/data, not taken at face value.

---

## 1. Done & verified this session

**Bucket A — M1-local code (all DONE, tests green):**
- **A1 cluster COM:** `cluster_com` (read-weighted centroid, floored, 0-based) added to both
  clusterers in `analyze/clustering.py` → `cpa_clusters.tsv` (std + manifest). Tests:
  `test_analyze.py::TestClusterCenterOfMass`.
- **A2 aligner stats:** `by_aligner_combo` surfaced in `*.consensus_aligner_stats.tsv`
  (`processing_stats.py`) + consensus HTML report (`summary.py`). New
  `tests/test_consensus_stats_surfacing.py`.
- **A3 manifest bedgraph/genomic-dist:** already implemented (Pass 3, `manifest.py:744`) — TODO stale.
- **A4 splice terminology:** dark figure regenerated; README + ARCHITECTURE.md relabeled with
  output→literature mapping + "novel = not annotated" caveat. Internal labels kept (non-breaking).
- **A5 indel-test verification:** the T=T(+)/A=A(-) inversion VERIFIED correct vs production
  (`walkback` stop_base); stale class docstring + NOTE block fixed.

**Bucket B:**
- **B1 deSALT Sherlock:** binary already = working bioconda build everywhere (md5 `e923d866…`);
  3 runtime mitigations in `run_desalt`. Empirically PASSED on Sherlock compute node (job
  31178641: 33,567-byte BAM, 34 mapped, no SIGSEGV).
- **B2 (NEW-082) gapmm2:** ROOT-CAUSED + FIXED + VERIFIED. gapmm2 25.4.13+'s
  `splice_aligner_minimap2` drops unspliced reads (requires minimap2 `ts:` tag, emitted only for
  spliced alignments). 25.8.12→10/35, 25.4.5→35/35. Fix: installer pinned `gapmm2==25.4.5`
  (`install_aligners_command.py`, matching pyproject); M1 env downgraded.
- **B3 validation bundle:** VERIFIED ALREADY CURRENT (not stale). cat2_plus_1 minimap2 BAM is
  through-aligned (`…49M9D39M8S`); gapmm2 36/36; `test_validation_reads.py` green.
  **Recommended NOT regenerating** (would inject align-env noise into a passing regression
  backbone) — user approved regen under the stale premise; surfaced for override.
- **B4 netseq reconcile:** DONE — Sherlock `correct_command.py` is byte-identical to M1 HEAD.
- **B4 aligner eval:** Winnowmap2 + Minisplice already wrapped; GMAP added (`7b32fa0`). GLASS = only
  unevaluated candidate (DEFERRED, see Open).

**Bucket C:**
- **C1 RN-sidecar:** leave archived (user decision). No action.
- **C2 Dorado pt:i:** IMPLEMENTED + tested. `--use-dorado-polya` (default off) threaded
  helper→`polya_trimmer`→`bam_processor`→`parallel`→`correct_command`. Always records
  `dorado_polya_length` in the result dict (TSV schema unchanged → golden hashes safe);
  authoritative `polya_length` when enabled. Tests:
  `test_polya_trimming.py::TestDoradoPolyaIntegration` (7).
- **C3 class rename:** keep internal labels + document mapping (done in A4).

Tests: full `not slow` suite green at session start (1577) and after Bucket A+B2 (1584). A
final post-C2 run is in flight (bg `bvfsmyi4j`) — CHECK IT (see Resume).

## 2. Open
- **GLASS aligner eval (only remaining item).** Not installed anywhere on Sherlock. Needs:
  install the 2025 graph-learning splice aligner (repo + weights + env), smoke on the validation
  reads, compare orthogonality vs the panel; add `run_glass()` + `SUPPORTED_ALIGNERS` ONLY if it
  earns its place. Open-ended research task — deliberately deferred.
- **B3 regen** — user approved but it's UNNEEDED (bundle already current). Only run if the user
  still wants a fresh regen for another reason; do it on H2/Sherlock (all 5 aligners), then
  `update_validation_aligner_bams.py` → `rectify correct` → `generate_review_report.py --arm drs`.

## 3. Resume
- **First:** `cat <scratchpad>/fulltest2.out` (bg `bvfsmyi4j`). If any FAILED, triage — most
  likely a C2-plumbing edge in `bam/parallel.py` (the use_dorado_polya threading) or a
  gapmm2-version-sensitive test. Targeted modules already passed (73 + golden hash).
- **Commit:** user has NOT asked. Surgical-staged in working tree. If asked: explicit `git add`
  only (§4), target `drs-validation-rebuild`. Do NOT `git add -A`.
- **GLASS:** the only forward work; treat as a fresh scoped task.

## 4. Files touched (uncommitted)
- Code: `analyze/clustering.py`, `analyze/summary.py`, `bam/processing_stats.py`,
  `bam/bam_processor.py`, `bam/parallel.py`, `polya/polya_trimmer.py`,
  `commands/install_aligners_command.py`, `commands/correct_command.py`.
- Tests: `test_analyze.py`, `test_indel_correction.py`, `test_consensus_stats_surfacing.py` (new),
  `test_polya_trimming.py`.
- Docs: `README.md`, `docs/ARCHITECTURE.md`, `dev/TODO.md`, `dev/BUGS_TO_FIX.md`,
  `docs/figures/splice_classification_dark.{png,svg}` (+ light re-touched).
- Env: M1 base conda `gapmm2` 25.8.12 → 25.4.5.
- Pre-existing WIP NOT mine: `dev/COMPASS_HANDOFF.md`, `docs/ALIGNER_RECOMMENDATIONS.md`,
  `docs/figures/generate_splice_classification_v3.py`, untracked `dev/_*.py`.
- Cluster: Sherlock `/scratch/users/kevinroy/desalt_b1_smoke/` (B1 confirmation job artifacts).
