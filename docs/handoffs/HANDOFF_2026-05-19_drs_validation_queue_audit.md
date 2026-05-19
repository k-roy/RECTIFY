# HANDOFF 2026-05-19 — DRS validation queue audit + in_amb priority tests

**Branch:** `drs-validation-rebuild`
**Top commit:** `cc4acbc` — docs(handoffs): mark cat1_plus_1, cat1_plus_2, cat2_minus_2 RESOLVED
**Date:** 2026-05-19
**Follows:** `HANDOFF_2026-05-18_validation_reads_drs_to_cdna_to_qsrev.md`

---

## Suite baseline (end of session)

| Suite | Command | Result |
| --- | --- | --- |
| Validation reads (DRS) | `pytest tests/test_validation_reads.py` | **107 passed, 8 skipped** |
| Splice junction | `pytest tests/test_splice_junction.py` | **89 passed, 1 xfailed** |
| Full non-slow | `pytest -m "not slow"` | **1000 passed, 35 skipped, 4 deselected, 1 xfailed, 2 pre-existing failures** |

The 2 pre-existing failures are in `test_corrected_consensus_tiebreaker.py`:
- `test_paralog_tiebreaker_picks_multi_aligner_consensus`
- `test_majority_consensus_picks_chrXIV_even_when_outlier_has_wider_span`

Both fail with `'float' object has no attribute 'split'` in the `_eff_key` NaN path inside
`merge_corrected_tsvs`. Pre-existing before this session; not caused by anything in this branch.

---

## This session's work

### 1. cat3_plus_2 HP-ED trace (closes §2.2 queue item)

Completed the per-op HP-ED trace for read `79f61403` (cat3_plus_2) using:
- Trimmed BAMs: `rectify/data/validation/rectified/per_aligner/{aligner}.trimmed.bam`
- `HpPenaltyTable` loaded from `penalty_tables/penalty_scores.tsv` + `str_penalty_scores.tsv`
- Trace script: `/tmp/trace_cat3plus2_hped_v2.py`

**Root cause of 1.39-pt gap** (deSALT=26.21, mpb=27.60): single D1 op at ref 142730 in a
1-bp HP context. `penalty_table.del_cost(1, base)` = 0.36; flat fallback = 1.0. deSALT calls
this position as `1D 1I` (1.25+0.36=1.61); mpb calls it as `3X` (3.0). HP-ED correctly
recognises the deSALT alignment as cheaper. Not a bug.

Both aligners produce identical TSV output: `corrected_3prime=143380`, `junctions=142253-142619`,
`effective_group=A`. The gap is body-only; winner selection is correct.

**Memory saved:** `feedback_hped_penalty_table_required.md` — always load `HpPenaltyTable`
and use trimmed BAMs when reproducing HP-ED from `per_aligner_summary.tsv`.

Queue entry updated in commit `6f39855`.

### 2. in_amb vs donor_ok priority — synthetic unit tests (§3.1 item 1)

New class `TestInAmbVsDonorOkPriority` appended to `tests/test_splice_junction.py`
in commit `bf1749e`. Two tests:

- **`test_plus_in_amb_beats_canonical_out_of_window`**: shift=0 (in_amb, non-canonical NN)
  wins over shift=-2 (out-of-amb, canonical GT) because Step 2 (in_amb) outranks Step 3
  (donor_ok). Tuple ordering: `(False,True,0) < (True,False,2)`.
- **`test_plus_canonical_wins_when_both_in_amb`**: shift=+2 (in_amb, canonical GT) wins
  over shift=0 (in_amb, non-canonical AA) because both in_amb → Steps 2 tie → Step 3
  (donor_ok) breaks tie. Tuple ordering: `(False,False,2) < (False,True,0)`.

These are the two cases that 99558c1 was designed to handle; they're now load-bearing tests.

### 3. + strand equivalence-extension (§3.1 item 2) — already resolved

Commit `acb508e` (prior session) had already implemented the fix. The handoff was written
in the window between the disable commit (`01c2a18`) and the re-enable. Updated queue in
commit `940d7cc`.

### 4. 15 pre-existing test_splice_junction failures (§3.1 item 3) — already resolved

Commit `d562395` (prior session) classified all 15 as stale-fixture failures (`MockRead`
missing `is_unmapped`/`query_name`). All now pass or carry the expected xfail marker.

### 5. debugger_queue.md staleness audit (two rounds)

**Round 1 (prior to user screenshot):** Found 2 stale "Deferred" entries that `acb508e`
had already resolved but the queue didn't reflect:
- `cat3_minus_2` "Still deferred" bullet → updated to "RESOLVED by acb508e" (commit `8327cd2`)
- `cat3_plus_2` subsection → updated to "RESOLVED (acb508e, 2026-05-18)" (commit `8327cd2`)

**Round 2 (triggered by cat1_plus_1 screenshot):** User showed cat1_plus_1 IGV figure
confirming mpb wins with corrected_3prime=10611 — looking correct. Queue said "currently
failing" for cat1_plus_1, cat1_plus_2, and cat2_minus_2 had an open design question.
Verified all tests pass; found three more stale entries:
- `cat1_plus_1`: RESOLVED by commit `09e4627` (early-exit window widened to ±20 bp)
- `cat1_plus_2`: RESOLVED (Phase A fix `a1728eb` was sufficient)
- `cat2_minus_2`: RESOLVED by commit `6943450` (2-bp del extension; endpoint=128096 per
  user directive 2026-05-18)

All three updated in commit `cc4acbc`.

---

## Validation read status (cat1–cat4 exact; cat5–cat9 by test class)

All tests in `test_validation_reads.py`: **107 passed, 8 skipped**.

### Cat1 — indel correction (all 4 PASS)

| Read | Locus | corrected_3prime | Notes |
| --- | --- | --- | --- |
| cat1_plus_1 | chrXIV:10435-10611 + | 10611 | mpb winner; early-exit widened |
| cat1_plus_2 | chrI:31118-31546 + | 31546 | |
| cat1_minus_1 | chrII:9826-10558 − | 9834 | |
| cat1_minus_2 | chrXII:15345-15964 − | {15345, 15351} | either accepted |

### Cat2 — soft-clip rescue (all 4 PASS)

| Read | Locus | corrected_3prime | Notes |
| --- | --- | --- | --- |
| cat2_plus_1 | chrI+ 23754 | 23754 | deSALT wins; calibration issue open (see below) |
| cat2_plus_2 | chrVI+ | 8605 | |
| cat2_minus_1 | chrV− | 186 | |
| cat2_minus_2 | chrI− | 128096 | 2-bp del extension; 6943450 |

### Cat3 — 5'-rescued splice junction (all 4 PASS)

| Read | Locus | 5'_position | Notes |
| --- | --- | --- | --- |
| cat3_plus_1 (0a28167d) | chrII:168808-169462 + | 168423 | |
| cat3_plus_2 (79f61403) | chrI:142618-143383 + | 142252 | HP-ED gap confirmed body-only |
| cat3_minus_1 (ac4db6da) | chrXV:900071-900767 − | 901193 | |
| cat3_minus_2 (28ea9379) | chrII:365845-366503 − | 366584 | |

### Cat4 — intron near 3' end (all 4 PASS)

| Read | Locus | corrected_3prime | Intron coords |
| --- | --- | --- | --- |
| cat4_plus_1 | chrXI:19592-22073 + | 20503 | 20527-22047 (1520 bp) |
| cat4_plus_2 | chrX:392246-393837 + | 393721 | 393725-393825 (100 bp) |
| cat4_minus_1 | chrI:128094-129063 − | 128117 | 128521-129021 (500 bp) |
| cat4_minus_2 | chrIX:76016-77313 − | 76254 | 76027-76250 (223 bp) |

### Cat5–Cat9 — summary by class (all passing or skipped-by-design)

- **Cat5** (chimeric reads): 4 reads on chrV+, chrII+, chrVII−, chrIII−. `test_category_tags` checks `Xz=1` for 3/4; cat5_plus_2 intron is confidently single-aligner. All pass.
- **Cat6** (intron near 5' end): 4 reads on chrII+(×2), chrII−, chrIV−. All pass; intron at annotated coords.
- **Cat7** (various): All pass (counts and shifts).
- **Cat8** (multi-peak / single-peak): Single-peak fraction tests pass.
- **Cat9** (Module-2H coverage gap): 4 reads skipped by design (`pytest.skip` when Module 2H refinement doesn't reach canonical output). These are the 4 skips in the 107+8 tally.

Winner counts from `regen_pa_rest_bundle.py`: `mm2=0, gapmm2=3, mpb=8, deSALT=19, uLTRA=6`.

---

## Open items (carry forward)

### Genuinely open — code change needed

**cat2_plus_1 calibration** — `del_cost(hp_len≈19, base='A')` in
`penalty_tables/penalty_scores.tsv` may be too expensive. deSALT (23754, 3× single-base
insertions in AAAT tetramers) wins over minimap2 (23759, 1× 8-bp HP deletion in a long
A-run). User diagnosis: minimap2's representation is more parsimonious; HP-ED should favor
it. Test currently hardcodes 23754 to accept the current (wrong) behavior — fix requires
recalibration of long-A-run del_cost and updating the test assertion to 23759. Cluster job
needed to extend empirical INS AT table for hp_len 4-20.

**xfail `test_plus_offset_junction_rescued` fix (option C)** — `_off` search window in
`_rescue_3ss_truncation_body` allows shift=0 (annotated NN, in_amb=True) to score the
same as shift=-3 (canonical GT, in_amb=False) because the window is shifted by `_off=3`.
Fix: track `_eff_intron_start -= _off` in scoring loop so recorded donor position matches
where the window actually landed.

### Design / deferred

**Ambiguity-window + motif-strength tiebreaker in `merge_corrected_tsvs`** — when two
aligners differ only by ±1 within the ambiguity window and motif strengths differ, the
current HP-ED sum tiebreaker may not surface the biologically correct winner. Deferred
until a real example surfaces.

**Two pre-existing `test_corrected_consensus_tiebreaker.py` failures** — `_eff_key` NaN
path in `merge_corrected_tsvs`; fix requires tracing the NaN source. Not caused by this
branch.

---

## Next protocol

### §3.2 — cDNA validation pass

```bash
cd /Users/kevinroy/work/rectify
pytest tests/test_validation_reads_cdna.py -v 2>&1 | tail -20
```

Read context files before starting:
- `memory/project_cdna_refactor.md`
- `memory/project_cdna_pipeline.md`
- Prior handoff on cDNA phase D (search `docs/handoffs/` for `cdna_umi_phase_d`)

DRS learnings to apply:
- HP-ED penalty table is required; load from `penalty_tables/`
- Use trimmed BAMs, not pA_rest BAMs
- Verify `splice_aware_5prime.py` + strand equivalence-extension fires for cDNA reads

### §3.3 — QuantSeq REV pass

- Walkback integration: verify `11_polya_walkback_recompute.py` strand-flip patch is
  captured in rectify integration (see `rectify_design_docs/ISSUE_walkback_integration.md`)
- TTTKGCAA validation: spot-check 3 internal-priming reads from `cold_shock_DRS_16c_2026`
- AG-priming threshold check: ensure `--ag-priming-threshold` default is calibrated for
  cold-shock vs standard conditions

---

## Resume command

```bash
cd /Users/kevinroy/work/rectify
git log --oneline -5
pytest tests/test_validation_reads.py -v --tb=short 2>&1 | tail -15
```

Read this file and `docs/handoffs/debugger_queue.md` before making any changes.
