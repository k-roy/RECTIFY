# HANDOFF 2026-05-19 — §3.2 cDNA pass + §3.3 QSrev pass

**Branch:** `drs-validation-rebuild`
**Top commit:** `245f486` — docs(handoffs): finish archiving stale session handoffs + add status preamble
**Date:** 2026-05-19
**Follows:** `HANDOFF_2026-05-19_drs_validation_queue_audit.md`

---

## Suite baseline (end of session)

| Suite | Command | Result |
| --- | --- | --- |
| DRS validation reads | `pytest tests/test_validation_reads.py` | **107 passed, 8 skipped** |
| cDNA walkback | `pytest tests/test_validation_reads_cdna.py` | **16 passed** |
| cDNA calibration + bundling | `pytest tests/test_data_bundling.py tests/test_profiler_umi_binning.py` | **31 passed** |
| QSrev walkback | `pytest tests/test_quantseq_rev_walkback.py` | **14 passed** |
| Splice junction | `pytest tests/test_splice_junction.py` | **89 passed, 1 xfailed** |
| Full non-slow | `pytest -m "not slow"` | **1002 passed, 33 skipped, 4 deselected, 1 xfailed, 2 pre-existing failures** |

Pre-existing failures remain in `test_corrected_consensus_tiebreaker.py`:
- `test_paralog_tiebreaker_picks_multi_aligner_consensus`
- `test_majority_consensus_picks_chrXIV_even_when_outlier_has_wider_span`

Both fail with `'float' object has no attribute 'split'` in `_eff_key` NaN path.

---

## §3.2 — cDNA validation pass

### 3.2.1 Walkback baseline — clean

`tests/test_validation_reads_cdna.py` — **16/16 pass** (2.15 s).
Tests exercise `walk_back_anchor_and_tail` on 12 pre-built validation reads
(cat1_cdna_polya, cat4_cdna_false_junc, cat5_cdna_chimeric; 4 reads each).

### 3.2.2 cDNA calibration — Phase A-C substantially complete (uncommitted)

Prior sessions completed the full UMI-bin-stratified penalty calibration.
State discovered at session start was ahead of the `cdna_umi_phase_d_resume.md`
handoff (which described the work as pending).

| Component | Status | Notes |
| --- | --- | --- |
| Phase A: profiler UMI-bin stratification | Done, uncommitted | `scripts/calibration/empirical_cigar_error_profiler.py` |
| Phase B: bundler + resolver bin-aware | Done, **committed** | `rectify/data/__init__.py` |
| Option 1: sidecar BAM lookup (fixes XC-tag drop) | Done, uncommitted | profiler + `run_profiler_qsrev_h2.sh` |
| Per-bin TSVs (umi1/umi2/umi3plus) | Done, **untracked** on disk | `penalty_tables/penalty_scores_cdna_umi*.tsv` |
| Phase C: `PenaltyTableSet` class | Done, **committed** | `junction_refiner.py:192` |
| Phase C: wired into `correct_command.py` | Done, **committed** | `correct_command.py:500-566` |
| Phase C smoke test | Done, uncommitted | `tests/test_data_bundling.py:191` |
| Convergence check (§4.9) | Done | `runs/error_profile_cdna_20260518/CONVERGENCE.md` |

**Convergence:** umi1 median |Δ%|=3.1%, umi2=3.0%, umi3plus=2.0% — all within
threshold. 8-chunk calibration is sufficient.

**What's still uncommitted (requires explicit commit request):**
- `scripts/calibration/empirical_cigar_error_profiler.py`
- `scripts/calibration/run_profiler_qsrev_h2.sh`
- `tests/test_data_bundling.py` (5 new test functions: Phase B/C smoke tests +
  per-bin flip test + sidecar test)

**What's untracked (also needs to be committed):**
- `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv`
- `penalty_scores_cdna_umi1.tsv`, `penalty_scores_cdna_umi2.tsv`,
  `penalty_scores_cdna_umi3plus.tsv`
- `docs/handoffs/cdna_umi_phase_d_resume.md`
- `docs/handoffs/cdna_umi_stratified_calibration.md`

### 3.2.3 DRS learnings applied to cDNA

- **HP-ED penalty table**: Phase C loads per-bin tables via `PenaltyTableSet`,
  selected by `XC` tag on each read. Committed and working.
- **Trimmed BAMs**: calibration used `stage1_consensus.bam` as sidecar (not
  pA_rest BAMs).
- **splice_aware_5prime for cDNA**: `rescue_3ss_truncation` has no
  `--ONT-cDNA` gate — it fires for all reads with annotated junctions.
  For orient=fwd cDNA reads, BAM 5' = RNA 5', so the rescue is appropriate.
  Full-pipeline verification (end-to-end `rectify correct --ONT-cDNA` on the
  12 validation reads) is not yet in the test suite; the current suite only
  tests the walkback baseline.

---

## §3.3 — QuantSeq REV pass

### 3.3.1 Walkback integration — VERIFIED ✅

`rectify/core/correct/protocols/quantseq_rev.py` exists and correctly implements:
- `is_reverse=True  → gene_strand='+'` (antisense chemistry)
- `is_reverse=False → gene_strand='-'`
- Delegates to `walkback_3prime_with_qpos` (protocol-agnostic core in `walkback.py`)

This matches the fix required in `ISSUE_walkback_integration.md` (line 208 strand
bug from `11_polya_walkback_recompute.py`). All acceptance criteria in the issue
are checked. Confirmed by `tests/test_quantseq_rev_walkback.py`: **14/14 pass**.

### 3.3.2 AG-priming threshold (Han QSrev — standard conditions) — OK

CLI default: `--ag-threshold 17.0` (argparse default in `correct_command.py:1216`).
Validated on Han 2023 wt_R1 50k subsample (per ISSUE doc): High 92.2% / Low 7.8%.
This is the correct rate for a yeast QuantSeq REV dataset with few internal-priming
artifacts. No recalibration needed.

Note: `parallel.py` still has function-level default `ag_threshold=0.65` — this is
shadowed by the value passed from `correct_command.py` (config `'ag_threshold': 17.0`).
The function-level default is never the operative value when called through the
normal CLI path. If future callers invoke the parallel functions directly, they
should pass the explicit 17.0.

### 3.3.3 Cold-shock DRS spot-check — deferred (future project)

The cold-shock DRS dataset (`by4742-wt-upf1D_cold_shock_DRS_16c_2026` on H2)
has never been analyzed with RECTIFY. TTTKGCAA validation and any cold-shock-
specific threshold calibration are future-project scope. The cold-shock DRS
analysis is separate from the Han QSrev analysis.

Deposit on H2:
`/u/project/guillom/shared/raw/by4742-wt-upf1D_cold_shock_DRS_16c_2026/`
(6 BAMs: wt_rep1-3, upf1d_rep1-3; chrI has ~135k reads in wt_rep1).

---

## Open items (carry forward)

### Genuinely open — code change needed

**xfail `test_plus_offset_junction_rescued` (Option C fix)** — `_off` search
window in `_rescue_3ss_truncation_body` lets shift=0 (in_amb=True, non-canonical
NN) score the same as shift=-3 (canonical GT, in_amb=False) when the window is
shifted by `_off=3`. Fix: track `_eff_intron_start -= _off` in the scoring loop.
Test is currently xfail with full diagnosis in the reason string.

**Two pre-existing `test_corrected_consensus_tiebreaker.py` failures** —
`_eff_key` NaN path; requires tracing NaN source in `merge_corrected_tsvs`.

### Uncommitted working-tree changes (need explicit commit request)

See §3.2.2 table above. Phase A/B/C cDNA calibration is done but not committed.
Commit order recommendation:
1. cDNA calibration scripts + test (Phase A + Option 1)
2. cDNA per-bin TSVs (data files)
3. cDNA handoff docs

### Design / deferred

**Ambiguity-window + motif-strength tiebreaker in `merge_corrected_tsvs`** —
see `debugger_queue.md` design note.

**Cold-shock DRS analysis** — first `rectify correct --direct-rna` run on
`by4742-wt-upf1D_cold_shock_DRS_16c_2026`; TTTKGCAA motif check deferred.

---

## Resume command

```bash
cd /Users/kevinroy/work/rectify
git log --oneline -5
pytest tests/test_validation_reads.py -v --tb=short 2>&1 | tail -5
pytest tests/test_validation_reads_cdna.py -q
pytest tests/test_quantseq_rev_walkback.py -q
```

Read this file and `docs/handoffs/debugger_queue.md` before making changes.
