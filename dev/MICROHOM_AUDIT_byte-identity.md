# MICROHOM-DRIFT GUARD — BYTE-IDENTITY AUDIT

Auditor lens: **BYTE-IDENTITY**. Claim under audit: with `microhom_drift_margin=0.0` (default),
the microhomology-drift guard is truly inert — output is byte-identical to the pre-guard code.

## Commits
- Guard feat commit: `3a716aa` feat(refiner): microhomology-drift guard (general non-HP drift), flag-gated + byte-identical off
- Follow-up test/fix: `17e57ba` test(guard): microhomology-drift guard regression tests + detector geometry fix (both directions)
- Parent (pre-guard) for `junction_refiner.py`: `69a230f`

## Audit plan
1. (a) Read the guard code path; confirm the guard block is unreachable / a no-op when `microhom_drift_margin == 0.0`.
2. Confirm threading of `microhom_drift_margin` default through refine_read_junctions -> _run_sequential/_run_parallel -> refine_bam_junctions, and CLI default.
3. (b) Independent harness: refine a sim panel (mix_r3b_out / mix_fair_out aligned BAM) with guard OFF vs pre-guard commit; diff emitted BAMs byte-for-byte (js/je + CIGAR + read order).
4. (c) Full suite -m "not slow" at default must be green.
5. (d) Hunt for any FP-order / detector side-effect that changes behavior even at margin=0 (e.g., detector always runs, float compares, ordering of candidate evaluation, short-circuit changes).

## Checkpoints

### CK1 — (a) Code-path read: guard block gating at default. DONE.
Diff `69a230f..3a716aa` on `junction_refiner.py` is +156/-3 lines. The guard adds:
- New pure helpers `_hp_run_across`, `_frac_match`, `_move_microhomology` (only CALLED inside `margin>0.0` blocks).
- 5 new kwargs on `refine_read_junctions`, `_run_sequential`, `_run_parallel`, `refine_bam_junctions`, all
  defaulting `hold_margin=0.0, hp_drift_margin=0.0, hp_drift_min_run=4, microhom_drift_margin=0.0,
  microhom_threshold=0.5`, threaded through worker-pool kwargs dict and `_refine_read_batch`.

Hot-path logic (lines 684-804). UNCONDITIONAL new code at default:
- L685 `incumbent_score = None` (new local)
- L736-737 `if is_alt == 0: incumbent_score = score_cmp` (sets local only; does NOT touch best_tuple)
- L762 `best_score_cmp = best_tuple[0]` (new local)
- L765 `moves = (new_js != ns or new_je != ne)`
- L775 `eff_margin = hold_margin` (= 0.0 at default)
- L776 `if moves and hp_drift_margin > 0.0:` -> False at default -> `_hp_run_across` NOT called
- L789 `if moves and microhom_drift_margin > 0.0:` -> False at default -> `_move_microhomology` NOT called
- L794 `if moves and eff_margin > 0.0 and incumbent_score is not None:` -> `eff_margin=0.0` so `eff_margin>0.0`
  is False -> veto block SKIPPED -> `moves` never mutated
- L801 `if moves:` == old `if new_js != ns or new_je != ne:`

=> At default, the emitted `replacements` list is logically identical to pre-guard. Guard/detector functions
are never invoked. VERDICT (code read): inert at default. Independent BAM-diff harness pending (CK3).

### CK2 — Default propagation through all entry points. DONE.
Production caller `rectify/core/commands/correct_command.py:746` calls `refine_bam_junctions` WITHOUT passing
any of hold_margin/hp_drift_margin/microhom_drift_margin/microhom_threshold -> all take function defaults
(0.0/0.0/0.0/0.5). No CLI flag surfaces the guard yet (default-off, not user-reachable in production).
Worker path: `_refine_read_batch` reads `state['kwargs']` (built at L1372 incl. microhom_* = passed value)
and forwards `**eff_kw` to `refine_read_junctions`. At default these are 0.0/0.5. Sequential path forwards
identically. Both paths symmetric.
Benchmark scripts (scripts/benchmark/noncanon_sim/*) DO set hp_drift_margin / hold_margin / motif_blind, but
those are dev harnesses, not the production default path.

### CK3 — Confounder isolation. DONE (IMPORTANT).
The guard feat commit `3a716aa`'s DIRECT PARENT is `b6f07f7` (NOT `69a230f`). Between `b6f07f7` and current
HEAD, git confirms **ONLY `junction_refiner.py` changed** in the entire `rectify/` tree:
  - `git diff --stat b6f07f7 3a716aa` = junction_refiner.py +74 lines only.
  - `git diff --stat b6f07f7 HEAD -- rectify/core/splice/junction_scoring.py` = EMPTY (identical).
  - `git diff --name-only b6f07f7 HEAD -- rectify/` = junction_refiner.py ONLY.
NOTE a trap avoided: the range `69a230f..3a716aa` ALSO contains the concat-DP change (`_USE_CONCAT_DP`,
DEFAULT ON) in junction_scoring.py — a *separate* perf change. Diffing HEAD-default vs `69a230f` would
CONFLATE guard + concat-DP. Comparing against the DIRECT PARENT `b6f07f7` isolates the guard alone (scoring
modules byte-identical). Harness therefore compares HEAD vs `b6f07f7`.

Also: the follow-up commit `17e57ba` (detector geometry fix) touches ONLY the body of `_move_microhomology`,
which is called EXCLUSIVELY inside `if moves and microhom_drift_margin > 0.0:` — never reached at default. So
the whole guard (feat + geometry fix) is behind the margin>0 gate.

### CK4 — Independent BAM-diff harness (SEQUENTIAL, n_workers=1). DONE — BYTE-IDENTICAL.
Harness: /private/tmp/.../scratchpad/harness_byteident.py. Loads current
`rectify.core.splice.junction_refiner` (HEAD, guard present, microhom_drift_margin defaults 0.0) AND the
`b6f07f7` source as a sibling module (no guard code at all). Runs `refine_bam_junctions` on each panel's
`aligned.sorted.bam` (aligner_bams=[in_bam] -> real candidate pool -> real moves), sort_and_index=True,
n_workers=1, for two configs: `default` and `motif_blind=True` (the shipped re-placer config). Byte-compares
SHA256 over full SAM record text (read order + CIGAR + all fields).

RESULTS (all MATCH; `refined` counts identical HEAD vs PRE):
  mix_r3b_out/default    : refined 815/815, 5100 recs, sha256 MATCH
  mix_r3b_out/motifblind : refined 334/334, 5100 recs, sha256 MATCH
  mix_fair_out/default   : refined 585/585, 5800 recs, sha256 MATCH
  mix_fair_out/motifblind: refined 580/580, 5800 recs, sha256 MATCH
  fab_sweep/default      : refined 152/152, 2100 recs, sha256 MATCH
  fab_sweep/motifblind   : refined 136/136, 2100 recs, sha256 MATCH
OVERALL: BYTE-IDENTICAL across all panels/configs. Rich diff surface (up to 815 reads moved) exercised the
guard-off path against genuine junction moves. No js/je, CIGAR, or read-order divergence at default.

### CK5 — (d) FP-order / detector side-effect at default. NONE FOUND.
Key structural finding: the microhom feat commit `b6f07f7->3a716aa` has **ZERO removed/changed lines**
(git: `grep -cE '^-[^-]'` = 0; pure +74). The unconditional machinery people might worry about
(`moves = ...`, `incumbent_score`, `best_score_cmp`, `if moves:`) ALREADY EXISTED in the direct parent
`b6f07f7` (introduced by the earlier, already-shipped HP-drift guard). The microhom commit adds ONLY:
  1. pure helpers `_frac_match`, `_move_microhomology` — called nowhere at default;
  2. one gated block `if moves and microhom_drift_margin > 0.0:` — condition False at default (0.0), so
     `_move_microhomology` never runs and `eff_margin` is never incremented;
  3. kwargs with `=0.0`/`=0.5` defaults, threaded through worker-pool dict + `_refine_read_batch`.
=> No unconditional hot-path statement is added/altered by the microhom commit. `incumbent_score` is only
READ inside `if moves and eff_margin > 0.0 and incumbent_score is not None:` where `eff_margin>0.0` short-
circuits to False at default. No float compare, no candidate-eval reordering, no detector-side-effect, no
FP-order change. The `profile.inc('microhom_drift_flagged')` counters only fire inside the gated block, so no
side-effect on profile output either. Detector `_move_microhomology` never executes at default.

### CK6 — Parallel-path (n_workers=3) plumbing check. DONE — BYTE-IDENTICAL (order-normalized).
Exercises the worker-pool kwargs dict (`_WORKER_POOL_STATE['kwargs']` incl. microhom_*) + `_refine_read_batch`
`**eff_kw` forwarding. macOS notes: the parallel path is written for Linux `fork` (docstring L1334; state
inherited via fork, not pickled). Under the macOS default `spawn`, workers don't inherit `_WORKER_POOL_STATE`
(`KeyError: 'header'`) — a PLATFORM issue affecting HEAD and PRE identically, unrelated to the guard. Forcing
`mp.set_start_method('fork', force=True)` reproduces the production fork path.

SELF-CONSISTENCY diagnostic (harness_par_selfconsistency.py) — CRUCIAL:
  - HEAD-vs-HEAD, nw=3, RAW record order: **DIFFER**  (parallel batch-completion scheduling is non-deterministic
    in output order — this is inherent to the parallel path, NOT the guard; same code differs run-to-run).
  - HEAD-vs-HEAD, nw=3, POSITION-SORTED: **MATCH**.
  - nw=3 vs nw=1(sequential), position-sorted: **MATCH** (parallel == sequential content).
So the correct parallel comparison strips the guard-INDEPENDENT order noise via position-sort.

PARALLEL HEAD-vs-PRE (harness_par_final.py, nw=3, position-sorted digest):
  mix_r3b_out/default    : refined 815/815, MATCH
  mix_r3b_out/motifblind : refined 334/334, MATCH
  mix_fair_out/default   : refined 585/585, MATCH
  mix_fair_out/motifblind: refined 580/580, MATCH
  fab_sweep/default      : refined 152/152, MATCH
  fab_sweep/motifblind   : refined 136/136, MATCH
OVERALL: BYTE-IDENTICAL. Worker-pool kwargs plumbing carries the (default 0.0/0.5) microhom kwargs with no
effect. NOTE the raw-order non-determinism is a pre-existing parallel-path property (not introduced by the
guard); the sequential path (CK4) is fully deterministic AND byte-identical.

### CK7 — (c) Full suite -m "not slow". DONE — GREEN modulo one unrelated missing-fixture.
FIRST RUN `pytest -m "not slow" -q` -> **1653 passed, 39 skipped, 4 deselected, 1 xfailed, 1 error**
(1304s wall — slow only because the 4-worker parallel harness ran concurrently, CPU-contended).
EVIDENCE-INTEGRITY CORRECTION (auditor self-check caught): the first run captured RC via
`pytest ... | tail -35 > out; echo PYTEST_RC=$?` — in zsh with no `pipefail`, `$?` is the exit code of `tail`,
NOT pytest. pytest returns non-zero on a setup error, so the true RC of that run was ~1, not 0. The doc's
earlier "RC=0" was a broken self-check; corrected here.

The single ERROR is `tests/test_restore_polya_from_parquet.py::test_restore_cat3_plus_2` — a FIXTURE-SETUP
error: the validation-bundle TSV
`scripts/validation_data/rebuild_2026_05/trimmed/validation_reads_polya_trim_metadata.tsv` is not materialized
in THIS worktree (the `trimmed/` dir is absent). PROVENANCE: orthogonal to the guard —
  - the test file references junction_refiner/microhom/refine_bam **0 times**;
  - it was NOT touched in the `b6f07f7..HEAD` (guard) range;
  - it's a missing-data/environment condition (the TSV is git-tracked from CI commit 99cc38f but not present
    in this worktree's working tree), not a code regression.

CLEAN RE-RUN (no pipe; true RC): `pytest -m "not slow"
--deselect tests/test_restore_polya_from_parquet.py::test_restore_cat3_plus_2 -q; echo TRUE_RC=$?`
-> **1653 passed, 39 skipped, 5 deselected, 1 xfailed, TRUE_RC=0** (1133s). This RC IS pytest's own exit
code (no pipe swallowing it). GREEN. (CK7b)
Guard-specific confirmation: `tests/test_microhom_drift_guard.py` = **14 passed** (direct run). The refiner
test suite (test_junction_refiner*, test_hp_drift_guard, test_microhom_drift_guard, ...) is within the passed
set. => (c) satisfied: suite green apart from the one unrelated, pre-existing missing-data fixture.

### SCOPE / LENS notes
- Byte-identity verified at the RECORD level: SHA256 over pysam `to_string()` per read (js/je + CIGAR +
  FLAG/POS/MAPQ/tags + read order) — exactly the invariant the task named. NOT claiming BGZF-container byte
  equality (compression framing can differ for identical content; record-level is the meaningful level).
- The "threshold 0.5 set on 5 fab / 2 real (under-powered)" hint is an EFFICACY/POWER concern about the
  operating point when the guard is ON — it cannot affect inertness at margin=0 and is OUTSIDE the
  byte-identity lens. Noted so the verdict is not read as missing it.

---

## VERDICT (byte-identity lens)

**The microhomology-drift guard is TRULY INERT at the default `microhom_drift_margin=0.0`. NO fault found.**

Evidence:
- (a) Code path: the guard block is gated on `microhom_drift_margin > 0.0`; at default the block is skipped and
  the detector `_move_microhomology` is never called. The unconditional machinery (`moves`, `incumbent_score`,
  `best_score_cmp`, `if moves:`) is NOT introduced by the microhom commit — it pre-existed in the direct parent
  `b6f07f7` (the microhom feat commit has ZERO deleted lines; pure +74). `incumbent_score` is only read inside a
  condition that short-circuits on `eff_margin > 0.0` (False at default). No FP-order/candidate-eval/detector
  side-effect at default.
- Confounder correctly isolated by comparing against the DIRECT PARENT `b6f07f7` (not `69a230f`), so the
  concat-DP change (`_USE_CONCAT_DP`, default ON) in junction_scoring.py — present in the wider `69a230f..3a716aa`
  range — does NOT contaminate the comparison (junction_scoring.py identical b6f07f7..HEAD).
- (b) Independent BAM diff, SEQUENTIAL (n_workers=1): 3 sim panels x 2 configs, up to 815 reads moved,
  SHA256 over record text (js/je + CIGAR + flags/pos/tags + order) = MATCH in all 6. `refined` counts identical.
- (b) PARALLEL (n_workers=3, fork): byte-identical after normalizing the guard-INDEPENDENT parallel-scheduling
  order noise (proven guard-independent via HEAD-vs-HEAD self-consistency: raw-order DIFFER, pos-sorted MATCH).
- (c) `pytest -m "not slow"` = 1653 passed (clean re-run, TRUE pytest RC=0 with the unrelated broken fixture
  deselected). The one error in the first run is an unrelated missing validation-fixture TSV (test references
  junction_refiner 0 times, not in guard commit range). `test_microhom_drift_guard.py` 14/14.
  NOTE: the first run's "RC=0" was mis-captured (`$?` of `tail`, not pytest, in a zsh pipe w/o pipefail);
  corrected via a pipe-free re-run. Verdict unchanged.
- (d) No non-inertness at default. Detector never runs; profile counters never fire; no ordering/FP change.

NON-BLOCKING ASIDE (task named "read order"): the PARALLEL path (n_workers>1) emits records in a
non-deterministic RAW order due to batch-completion scheduling — this is PRE-EXISTING and GUARD-INDEPENDENT
(HEAD-vs-HEAD self-consistency: raw-order DIFFER, pos-sorted MATCH). So strict ORDER-LEVEL byte-reproducibility
holds for the SEQUENTIAL (n_workers=1) path; the parallel path is byte-identical only after coordinate/position
normalization. Not a fault of the guard; noted so a reviewer isn't surprised by a raw-order parallel diff.

SCOPE OF THIS BLESSING: this audit clears ONLY the byte-identity/inertness-at-default leg (microhom_drift_margin
=0.0 == pre-guard). It does NOT endorse flipping the default ON — the operating-point/threshold-0.5 power
question (set on 5 fab / 2 real, under-powered) is the efficacy leg's concern, outside this lens.

Byte-identity verified at the RECORD level (the invariant the task named); container/BGZF-framing equality not
claimed (compression framing may differ for identical content — record level is the meaningful level).

## Files
- Code audited: rectify/core/splice/junction_refiner.py (guard commit 3a716aa + geometry fix 17e57ba)
- Pre-guard reference: git b6f07f7:rectify/core/splice/junction_refiner.py
- Harnesses: /private/tmp/.../scratchpad/harness_byteident.py (sequential),
  harness_par_final.py (parallel, pos-sorted), harness_par_selfconsistency.py (parallel self-determinism)
- Test log: /private/tmp/.../scratchpad/pytest_out.txt

STATUS: COMPLETE.

