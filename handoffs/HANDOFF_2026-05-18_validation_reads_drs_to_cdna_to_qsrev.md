# Handoff — validation read work (DRS → cDNA → QuantSeq REV)

**Date:** 2026-05-18 (late evening)
**Branch:** `drs-validation-rebuild`
**Top of branch at handoff:** `a133b54` (WIP — see `HANDOFF_2026-05-18_run_all_prelim_and_provenance.md` for the run-all+provenance work that's also in flight).
**Anchor commit for VALIDATION work specifically:** `e39089e` —
fix(rescue): propagate reanchor_clip_len → five_prime_soft_clip_length.
That's the last fully-tested validation-suite change. `a133b54` is
provenance scaffold (unrelated to validation reads).

---

## 1. Snapshot — where the suite is today

- `pytest tests/test_validation_reads.py` → **107 passed, 8 skipped.**
- `pytest tests/test_validation_reads.py -k cat4` → 11 passed, 1 skipped
  (false-N at 3' end).
- `false_junction_filter.validate_junction_filter()` → 3/3 synthetic cases
  pass.
- `pytest tests/test_corrected_consensus_tiebreaker.py` → 2 fail, 1 pass —
  **pre-existing**, in `_eff_key` NaN path; not regressions from anything
  this session shipped.
- `pytest tests/test_splice_junction.py` → 15 failures — **pre-existing**,
  all in `TestRescue3SSTruncation*` 5'-rescue unit tests. Confirmed
  pre-existing by stash-and-rerun on `e39089e`.

The validation bundle BAMs + TSVs at
`rectify/data/validation/rectified/` were regenerated against `e39089e`
and are committed.

---

## 2. This session's validation wins (cat3 cluster)

Three commits, all on this branch:

1. **`99558c1` fix(rescue): two-step scoring in `rescue_3ss_truncation`** —
   reordered the within-junction scoring tuple from
   `(not _donor_ok, not _in_amb, _shift_abs)` →
   `(not _in_amb, not _donor_ok, _shift_abs)` so that match-quality
   (in_amb) dominates signal-quality (canonical donor) when ED ties.
   Closed the queue top entry "Design: separate match-quality placement
   from canonical-signal slide in 5'-rescue."
   Pure structural cleanup on the current validation set — empty
   `per_aligner_summary.tsv` diff.
2. **`e39089e` fix(rescue): propagate reanchor_clip_len → five_prime_soft_clip_length** —
   the load-bearing fix this session. mpb runs in global mode so it
   never natively soft-clips; the reanchor pre-pass inside
   `rescue_3ss_truncation` collapses its `1X 2= 7I` 5' edge into a
   `10S` to enable downstream rescue. But the TSV's
   `five_prime_soft_clip_length` was still being computed from the
   *raw* (pre-reanchor) CIGAR via `extract_soft_clips`, so it stayed
   `0` for mpb. bam_writer's `extend_read_5prime_for_junction_rescue`
   gates on `soft_clip_len > 0` and returned False, leaving the
   reanchored 10S un-replaced in the final BAM. HP-ED 1.0/bp soft-clip
   penalty inflated mpb's score by 10–16 points on cat3_plus_1 /
   cat3_minus_1 despite the body alignment being structurally
   identical to deSALT's.
   Fix in `bam_processor.correct_read_3prime`: override
   `five_prime_soft_clip_len = _reanchor_clip_len` when reanchor fires
   materially. Concrete result: cat3_plus_1 mpb HP-ED **18.91 → 10.16**
   (mpb is now winner), BAM CIGAR head `10S 86=…` → `4M1I5M 384N 86=…`
   matching deSALT exactly.
3. **`2ee5d6c` docs(handoff)** — `e39089e` session wrap.

All four cat3 reads × 5 aligners preserve the documented post-regen
geometry. cat3_plus_2 verification target `14= 1D 9= 366N 50= …`
preserved on all 5 aligners.

---

## 3. Open items — finish DRS validation first

The user's framing: finish DRS validation, then move to cDNA, then
QuantSeq REV. Items below are tagged by protocol so the next session
can pick them in order.

### 3.1 DRS — open items (priority)

| Item | Where | Why it matters |
|---|---|---|
| 15 pre-existing failures in `test_splice_junction.py::TestRescue3SSTruncation*` | `tests/test_splice_junction.py` | 5'-rescue unit-level tests that pre-date the reanchor wiring. Need to be classified: stale-fixture, real regression, or known-skip with a reason. |
| **cat3_plus_2 HP-ED winner-selection** | `corrected_consensus.py` HP-ED scoring | mpb's raw alignment is already canonical; winner cluster has `1D 1=` exon-1 tail. HP-ED scores them within 1.4 points and picks the cluster (HP-ED 26.21 vs mpb 27.60). Earlier framing "HP-ED is undercounting the winner's split-D tail" hasn't been verified — trace per-op HP-ED contributions for cat3_plus_2 deSALT vs mpb and find which weight is responsible. Shares root cause with the open "Cat1 cluster (HP-ED metric)" entry. |
| **Cat1 HP-mode metric design** (architectural) | `docs/handoffs/debugger_queue.md` | mpb force-aligns 4+ mismatches into the body to anchor at a non-A past a poly-A tail. HP-ED currently undercounts this. Cross-cuts with cat3_plus_2. |
| **Defensive belt-and-suspenders** for bam_writer reanchor gate | `bam_writer.py` (3 write paths) | Optional hardening: gate `correction['reanchor_clip_len'] > 0` with `and correction.get('five_prime_rescued', False)`. Not done because the invariant is clear in `bam_processor.py:411-413`; the wasted call on a hypothetical regression is a no-op. |
| **+ strand 5'-rescue equivalence-extension proper-mirror fix** | `docs/handoffs/debugger_queue.md` "Bug note: + strand 5'-rescue equivalence-extension geometry inverted (DISABLED 2026-05-18)" | The + strand equivalence-extension was disabled because its geometry was inverted from the − strand mirror. Implementing the correct mirror is the structurally-symmetric undershoot case for cat3_plus_2 et al. |
| **Synthetic unit test for in_amb-vs-donor_ok priority** | new `tests/test_splice_5prime_priority.py` | Advisor-flagged after `99558c1`: the empty diff is two-way ambiguous (either the refactor is correct AND no real read exercises the new behavior, OR the refactor is a no-op AND no read exercises it). A synthetic test that asserts the in_amb-non-canonical vs out-of-amb-canonical ED-tie resolves to in_amb would make the priority-inversion claim load-bearing. |
| **Ambiguity-window + motif-strength tiebreaker for consensus** | `corrected_consensus.py` | Design note in `debugger_queue.md`. Symmetric slide already implemented in `junction_refiner.py:_apply_junction_replacement` lines 513-590; reuse that. |
| **Phase C 5' rescue calibration** | `scripts/calibration/` | Apply empirical HP penalty tables to the 5'-rescue scoring instead of the current 0.5/1.0 step function. The empirical profiler in `scripts/calibration/empirical_cigar_error_profiler.py` already exists. |

### 3.2 cDNA — when ready

The cDNA pipeline (ONT PCR-cDNA, PCB114 UMI architecture) has a
separate validation surface. Reference: memory `project_cdna_pipeline.md`
+ `project_cdna_refactor.md` + `tests/test_validation_reads_cdna.py`.

When DRS validation is at "no open regressions, all cat-categories
green," start the cDNA pass:

| First action | Goal |
|---|---|
| Run `pytest tests/test_validation_reads_cdna.py` | Establish cDNA-side baseline (parallel to the DRS baseline of 107/8). |
| Read `project_cdna_refactor.md` memory + `docs/handoffs/cdna_umi_phase_d_resume.md` (uncommitted in WIP, present locally) | Get full cDNA refactor status. Tasks 1-5 status, batch write_stage1_bam architecture, read_subtype naming. |
| Audit cDNA validation reads for the same fix categories | Apply learnings from DRS (reanchor, two-step scoring, equivalence extension). Some bugs may map; others won't because cDNA UMI architecture differs from DRS. |

### 3.3 QuantSeq REV — when ready

QuantSeq REV is dT-primed-cDNA, antisense, single-aligner BBMap path.
Smaller bug surface than DRS/cDNA but still has open items:

| Item | Where |
|---|---|
| Walkback integration into `rectify correct` | `rectify_design_docs/ISSUE_walkback_integration.md` — the `11_polya_walkback_recompute.py` post-hoc walkback should be moved into `correct`'s pipeline. |
| TTTKGCAA validation against an orthogonal yeast CPA dataset | Han 2023 prelim re-run (handoff #2) plus a TBD orthogonal dataset. |
| AG-priming threshold | `7bf0a1d/903cfc2` fixed default 0.65→17.0; verify no protocol still uses the wrong default. |

---

## 4. Resume command

```bash
cd /Users/kevinroy/work/rectify
git log --oneline HEAD~3..HEAD                       # sanity-check the chain
.venv/bin/pytest tests/test_validation_reads.py -q   # expect 107 passed, 8 skipped
.venv/bin/pytest tests/test_splice_junction.py -q    # expect 15 failed, 73 passed (pre-existing)
# Then pick an open item from §3.1 and dig in.
```

If the next session's first action is regen-the-bundle:

```bash
.venv/bin/python scripts/validation_data/regen_pa_rest_bundle.py
# Expected post-e39089e winner counts:
#   mm2=0  gapmm2=3  mpb=8  deSALT=19  uLTRA=6  (some samples flipped from pre-e39089e)
```

---

## 5. Cross-references

- Predecessor handoff: `HANDOFF_2026-05-18_reanchor_softclip_fix.md` —
  full verification trail for `e39089e`, including the open items in §4
  there (carried forward to §3.1 above).
- Two-step scoring handoff: `HANDOFF_2026-05-18_two_step_scoring.md` —
  for `99558c1`, including the advisor's empty-diff caveat.
- Run-all + provenance handoff (parallel workstream):
  `HANDOFF_2026-05-18_run_all_prelim_and_provenance.md`.

---

## 6. Memory notes worth checking before resuming

- `feedback_rectify_junction_slide.md` — junction_refiner's buggy
  realignment for the simple-slide case. Relevant if touching the
  `_apply_junction_replacement` path during 5'-rescue calibration work.
- `feedback_hp_edit_distance_semantics.md` — HP-ED is the metric the
  cat3_plus_2 winner-selection item is fighting. Soft-clip penalty
  is intentional (1.0/bp); the open question is whether the
  per-op HP-ED weights are correctly calibrated for cat3-class reads.
- `feedback_per_aligner_rescue_runs_first.md` — pipeline is
  correct-first: per-aligner `rectify correct` (incl. 5' rescue + 3'
  walkback) runs BEFORE `merge_corrected_tsvs`. Don't try to fix
  consensus-side issues from the per-aligner-correction side.
- `feedback_polya_tail_source.md` — when restoring poly-A tail as
  soft-clip, use parquet metadata, not the dorado_source softclip
  (which may contain bases the trimmer separated as non-poly-A).
