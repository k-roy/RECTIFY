# Audit — Lazy Corrected Consensus + Homopolymer Realign Speedup
Date: 2026-05-21
Auditor: independent review
Audits: `dev/specs/lazy_corrected_consensus_plan.md`
Scope: 85-file uncommitted WT on `drs-validation-rebuild`, M1 only.

## Status — 2026-05-21 (post-dispatch)

Three waves of Sonnet agents landed all P0/P1/P3 fixes plus one bug discovered during test authoring:

| Item | Status | Fix landed by |
| - | - | - |
| **P0** lazy fast-path unsoundness | ✅ FIXED | Wave 1 Agent A — `corrected_consensus.py:310-314` deleted; `_correction_requires_transient_edit` removed (no remaining callers); test suite 6/6 pass |
| **P1** silent fallback in `_chunk_index_from_path` | ✅ FIXED | Wave 1 Agent B — `sidecar.py:184-195` now raises `ValueError` on unparseable filenames |
| **P1** missing tests for targeted `realign_exon_blocks` | ✅ FIXED | Wave 2 Agent C — new `tests/test_realign_exon_blocks.py` with 13 tests covering all Stage 1 acceptance cases |
| **P1** degenerate equivalence test | ✅ FIXED | Wave 2 Agent C — `test_lazy_merge_distinct_aligners_select_same_winners` uses real perturbation (uLTRA `corrected_3prime` shifted ±50 bp on 18 of 36 reads); multi-aligner consensus test added |
| **P2** spec/code resolver mismatch | ✅ FIXED | Wave 3 Agent F — spec text in `lazy_corrected_consensus_plan.md` now matches actual 2-step resolver in code |
| **P3** dead `mv corrected.bam` line | ✅ FIXED | Wave 1 Agent B — `split_command.py:942-946` cleaned |
| **P3** no IGV opt-in flag | ✅ FIXED | Wave 1 Agent D — `--write-per-aligner-corrected-bams` CLI flag wired through `_correct_array_body` with samtools-sort+mv block on flag-on path |
| **NEW** `winning_aligner` column shadow in `merge_corrected_tsvs` | ✅ FIXED | Wave 3 Agent E — `_load_tsv` drops any pre-existing `winning_aligner` column on input; regression test `test_merge_corrected_tsvs_ignores_stale_winning_aligner_column` added |
| **P2** non-chunked `run` never calls `write_corrected_consensus_bam` | ✅ NOT-A-BUG (audit stale) | 2026-05-21 follow-up — `single_sample.py:253` (and `:641` for the post-hoc branch) already invokes `write_corrected_consensus_bam`, added in commit `22136af` (2026-05-18, "feat(run): thread per_aligner_corrected_bams through production HP-mode"). Comment at `single_sample.py:248` reads "Emit corrected_consensus.bam — symmetric with the chunked path." The audit predates verification; no code change needed |
| **P2** H2 rectify env lacks pyarrow | ✅ FIXED | 2026-05-21 — `conda install -c conda-forge pyarrow` into `/u/home/k/kevinroy/.conda/envs/rectify`; verified `import pyarrow` returns `24.0.0` |
| Phase C paired timing on real chunk data | ✅ DONE — see `dev/specs/lazy_consensus_timing_20260521.md` | 2026-05-22 Sherlock job 25662145 on `sh04-06n26` (AMD Milan, 4 CPUs). Dataset: `wt_tfiiib_rep3/chunk_008`, 3 aligners, ~15k corrected reads. **Lazy: 495.1 s wall, 286 MB RSS.** Full path: TIMEOUT (>81 min for minimap2 alone, 90-min limit exceeded). Lower-bound speedup: **> 9.9×**; true speedup estimated 20–30×. |
| Phase D reanchor optimization timing | ✅ DONE — see `dev/specs/lazy_consensus_timing_20260521.md#phase-d` | 2026-05-22 Sherlock job 25706050 on `sh04-11n19`. `_apply_reanchor_from_clip_len` (O(CIGAR-ops)) replaces per-base genome scan. **Total: 492.2 s (−2.9 s, 0.6% improvement).** Reanchor calls are fast in isolation (24.7 s for all 3 aligners) but the full merge is dominated by TSV/merge/UpSet overhead. Next bottleneck: make UpSet opt-in + profile merge internals + parallelize BAM write. |
| Phase E production timing + Lustre I/O diagnosis | ✅ DONE — see `dev/specs/lazy_consensus_timing_20260521.md#phase-e` | 2026-05-22 Sherlock job 25709013 on `sh04-11n19` (8 CPUs). `workers=3` (production config) reduces total from 492 s → **369.8 s (−25%)**. Root cause of remaining bottleneck identified: **OAK Lustre per-process client cache** — HP scoring subprocesses don't warm main process cache; cold BAM write = 146 s vs warm = 14 s. Fix: stage BAMs to `$SCRATCH` NVMe before merge → estimated total **~83 s**. ProcessPool shutdown also costs 22 s per chunk (fix: reuse pool). |

Test suite across affected files: **62 passed, 5 skipped** (skips are unrelated pyarrow streaming + pre-existing). Changes are in WT for user review; nothing committed.

The remainder of this audit document records the findings as originally written — kept for traceability. The status table above is authoritative for *current* state.

---

## P0 — Lazy path picks wrong winners (critical) — FIXED 2026-05-21

`tests/test_bam_writer_parallel_smoke.py::test_lazy_hp_edit_distance_matches_materialized_corrected_bam`
**fails on the WT right now**, on the same 36-read bundled fixture the spec's
checkpoint validation cites. 22/36 reads (61%) have HP-edit-distance values
that diverge from the materialized path by 20–160×. Spec's headline claim —
that lazy and materialized scoring choose the same winners — is empirically
false on the only fixture that exercises the new path.

**Root cause:** `corrected_consensus.py:310-314` is an unsound fast path:

```python
if _correction_requires_transient_edit(correction):
    apply_corrected_edits_to_read(read, correction, genome)
    edited += 1
else:
    fast_noedit += 1
```

`_correction_requires_transient_edit` (lines 222-259) inspects only the
correction-row dict. But `apply_corrected_edits_to_read` unconditionally calls
`realign_exon_blocks` (helper line 276) and `_hardclip_trailing_a_run` (line 335)
whenever `genome is not None`, regardless of correction-row state. Those edits
depend on **CIGAR + genome**, not the dict.

All 6 failing reads have `original_3prime=-1` (sentinel for "not recorded"),
which makes `_correction_requires_transient_edit` return False (line 257-258),
so the read is scored raw — even though the materialized path's writer
mutated its CIGAR via `realign_exon_blocks`. The reads diverge by exactly
the soft-clip-bases penalty (1.0/base) of the trailing A-run that
materialized hard-clipped but lazy did not.

**Blast radius:** path A only (lazy raw-BAM scoring,
`_read_hp_edit_distances_from_raw_bam`). Path B
(`_read_hp_edit_distances` reading a pre-materialized BAM)
reads CIGAR as-is and is unaffected. The lazy path is the spec's stated
default and the production policy explicitly disables corrected BAMs, so the
bug fires in the production-default code path.

**Why the spec's validation didn't catch it:** "145 passed locally" must
have excluded `test_bam_writer_parallel_smoke.py`, or the WT regressed
after the checkpoint write-up. The validation cited in the spec is not
load-bearing for the equivalence claim.

### Fix

Delete the outer fast path. `realign_exon_blocks` already has its own
internal X-at-HP pre-check (`read_edits.py:730-745`) that exits cheaply on
clean reads; the outer check was double-checking with the wrong predicate.

```diff
--- a/rectify/core/consensus/corrected_consensus.py
+++ b/rectify/core/consensus/corrected_consensus.py
@@ -307,11 +307,7 @@ def _read_hp_edit_distances_from_raw_bam(
                 correction = _lookup_read_correction(read, corrections, normalized_corrections)
                 if correction is None:
                     continue
                 seen_corrected_ids.add(_normalize_bam_read_name(str(correction.get('read_id', read_id))))
-                if _correction_requires_transient_edit(correction):
-                    apply_corrected_edits_to_read(read, correction, genome)
-                    edited += 1
-                else:
-                    fast_noedit += 1
+                apply_corrected_edits_to_read(read, correction, genome)
                 results[read_id] = (
                     _cigar_hp_edit_distance(read, genome, penalty_table),
                     _cigar_aligned_bases(read),
                 )
```

### Cost of the fix (fixture; lower bound on production scale)

| Path | Wall, 36 reads | Divergences |
| - | - | - |
| Materialized: full corrected-BAM write + scoring | 27.5 ms | 0 (baseline) |
| Lazy with buggy fast path (current WT) | 6.1 ms | **22/36** |
| Lazy after fix (always-apply) | **15.1 ms** | **0/36** |

The fix regresses fixture lazy timing by ~9 ms (2.5×) compared to the
buggy version. **It is still 45% faster than the materialized path** on the
fixture, and the materialized path's relative cost grows superlinearly
with read count (per-aligner `samtools sort` + index), so production-scale
savings are larger. Buggy-vs-fixed is not a useful comparison: the buggy
version is wrong.

Spec's bottleneck numbers (minimap2 +111.6 s, others >12.6 min on ysh1
chunk 000) for the corrected-BAM **write** step are not invalidated by the
bug — those are write/sort/index costs, not scoring. The lazy path still
avoids writing N full corrected BAMs. Net savings claim survives; the
absolute lazy timing in the spec is optimistic by an amount proportional
to how often the fast path fired.

`_correction_requires_transient_edit` and `fast_noedit`/`edited` counters
become unused — remove the dead helper as part of the same change.

---

## Audit findings by spec stage

| Stage | Spec claim | Audit result |
| - | - | - |
| Stage 0 — round-robin sidecar | "round-robin chunking maps local read index `j` in chunk `k` to `read_num = j * n_chunks + k`" | ✅ Formula is correct, verified by round-trip against mex67aa_rep1 chunk samples |
| Stage 1 — targeted `realign_exon_blocks` | "realigns only blocks containing homopolymer X" | ✅ Implementation correct (`read_edits.py:711-745`); **no unit tests** despite spec's "Acceptance: focused tests" |
| Stage 2 — shared edit helper | "single-threaded and parallel corrected-BAM writers call the same helper" | ✅ `apply_corrected_edits_to_read` defined once at `bam_writer.py:238`; 4 call sites funnel through it |
| Stage 3 — lazy HP scoring | "fixture test: full and raw-BAM lazy scoring choose same winners for every read" | ❌ See P0. The fixture test fails. |
| Stage 4 — lazy consensus BAM | "lazy final BAM is CIGAR/tag-equivalent to old pipeline" | ⚠️ Existing test passes for **single-aligner** (`minimap2` only); never exercised over a multi-aligner winner-flip scenario |

## Other findings (lower severity)

### P1 — Equivalence test surface too narrow

`test_lazy_merge_matches_materialized_corrected_bam_selection`
(`test_bam_writer_parallel_smoke.py:129-167`) passes the **same** BAM/TSV
for both aligners (`minimap2` and `uLTRA` both → `legacy_bam`/`fixture_bam`).
Winner selection is degenerate — every read ties — so the test cannot
catch a real divergence where mutation-in-memory and materialized produce
different CIGAR for the same read in one aligner only. Combined with the
36-read fixture size, the existing tests are insufficient surface area
for the spec's "every read" claim.

### P1 — Silent fallback in legacy sidecar reconstruction

`_chunk_index_from_path` (`sidecar.py:144-148`) returns `None` when the
filename doesn't match `_chunk_NNN_of_MMM`. `reconstruct_posthoc_sidecar_from_chunks`
then silently falls back to **iteration order** for the chunk index. For
the legacy CPA DRS / H2 mex67aa reconstruction this matters: a user passing
renamed or reordered chunk files (typo, glob expansion order, scp into a
flat dir) would produce silently-wrong `read_num` values. Strict mode does
not catch this.

```diff
--- a/rectify/core/chunking/sidecar.py
+++ b/rectify/core/chunking/sidecar.py
@@ -185,11 +185,11 @@ def reconstruct_posthoc_sidecar_from_chunks(
     indexed = []
     used_indices = set()
-    for fallback_idx, path in enumerate(chunk_paths):
+    for path in chunk_paths:
         idx = _chunk_index_from_path(path)
         if idx is None:
-            idx = fallback_idx
+            raise ValueError(
+                f"Cannot infer chunk index from filename: {path}. "
+                "Pass chunk_idx explicitly or rename to _chunk_NNN_of_MMM."
+            )
         if idx in used_indices:
             raise ValueError(f"Duplicate chunk index {idx} for {path}")
```

### P1 — Targeted-realign behavior is untested

Spec Stage 1 acceptance demands:
- "two short exon blocks where only one contains a homopolymer X; assert only that block changes"
- "neighboring clean short exon block remains unchanged"

Neither test exists anywhere in `tests/`. The code (`read_edits.py:711-745`)
looks correct on read, but a future edit that breaks `block_start` tracking
across N/S/H boundaries would not be caught.

### P2 — Identity resolver precedence: spec ≠ code

Spec says the lazy path's resolver should be RN-tag → sidecar QNAME →
legacy normalized QNAME → fail-loud. Code does legacy normalized QNAME + fail
only (`_lookup_read_correction`, `corrected_consensus.py:201-210`). No
correctness impact (TSV `read_id` IS the QNAME, suffix-strip aside) but
the spec is misleading. Either align the spec text to the code, or
implement RN-first if there's a perf reason.

### P2 — Spec acceptance for non-chunked `run` path not met for consensus BAM

`stages.py:240-243` correctly makes per-aligner corrected BAM opt-in. But
`single_sample.py` never invokes `write_corrected_consensus_bam` — the
non-chunked path emits a merged TSV but no final consensus BAM. Either
that's intentional (single-input runs don't need consensus) and the spec
should say so, or it's an omitted port and the chunked path's pattern
should be mirrored.

### P2 — H2 rectify env lacks pyarrow

`/u/home/k/kevinroy/.conda/envs/rectify` does not have pyarrow, so
`ReadNumSidecarWriter` falls back to a pickle file with a `.parquet`
extension. Downstream code that calls `pd.read_parquet(...)` on the
"sidecar parquet" crashes. Operational gap — not a code bug, but H2 is the
canonical location for legacy mex67aa reconstruction. Install pyarrow into
that env before running Stage 0 there.

### P3 — Dead `mv` line in chunk correction script

`split_command.py:943` does `[ -f "$SCRATCH_WORK/corrected.bam" ] && mv ...`,
but the chunk correction command no longer passes `--write-corrected-bam`,
so the file is never created. Cosmetic; remove for clarity.

### P3 — No user-visible flag for opting back into per-aligner corrected BAMs

The chunk-merge logic already falls back to materialized scoring when
`corrected_bams` are found on disk. But there's no CLI flag in
`rectify split --generate-slurm` to **request** them — a user wanting
IGV-loadable per-aligner BAMs must hand-edit the generated chunk script
to add `--write-corrected-bam`. Add a flag.

```diff
--- a/rectify/core/commands/split_command.py  (sketch)
+++ b/rectify/core/commands/split_command.py
@@  # CLI parser
+    parser.add_argument(
+        '--write-per-aligner-corrected-bams', action='store_true',
+        help='Also emit per-aligner corrected BAMs for IGV inspection '
+             '(default off; lazy consensus path does not need them).',
+    )
@@  # In _correct_array_body, when flag is set:
+    write_cb_line = (
+        '    --write-corrected-bam "$SCRATCH_WORK/corrected.bam" \\\n'
+        if write_per_aligner_corrected_bams else ''
+    )
@@ -928,7 +928,7 @@
 $PYTHON -m rectify correct \\
     "$LOCAL_BAM" \\
     --genome "{genome}" \\
     --annotation "{annot}" \\
 {jpt_line}
 {spt_line}
     $SCAN_CACHE_ARG \\
     --junction-pool-cache "$POOL_PKL" \\
 {aligner_bam_args}
+{write_cb_line}
     --threads "$CORRECT_CPUS" \\
     --streaming \\
     --checkpoint-dir "$CHECKPOINT_DIR" \\
     -o "$CHUNK_OUTDIR/corrected_reads.tsv"
```

---

## Phase C — timing notes

### What was measured (M1, 36-read fixture)

| Path | Wall | Divergences from materialized |
| - | - | - |
| Materialized (write + score) | 27.5 ms | 0 |
| Lazy current WT (buggy fast path) | 6.1 ms | 22/36 — invalid |
| Lazy fixed (always-apply) | 15.1 ms | 0 ✅ |

cProfile on lazy-fixed (50 iters × 36 reads = 1800 read-scans):

| Function | cumtime | % | per-read |
| - | - | - | - |
| `apply_corrected_edits_to_read` (helper) | 0.917 s | 71% | 510 µs |
| `apply_corrected_edits_to_read` self-time | 0.525 s | 41% | 290 µs |
| `realign_exon_blocks` | 0.197 s | 15% | 109 µs |
| `_cigar_hp_edit_distance` | 0.162 s | 13% | 90 µs |
| `align_exon_block_global` (inner DP) | 0.125 s | 10% | 50 calls / 1800 reads = **2.8% DP hit rate** |

**The realign DP fires on only 2.8% of fixture reads** because the
targeted pre-check rejects clean ones cheaply. Within
`apply_corrected_edits_to_read`, the dominant per-read cost on this
fixture is not the realign DP — it's `_decode_eq_seq_inplace` +
`clip_read_to_corrected_3prime` + `_hardclip_trailing_a_run`. Spec
implicitly assumed realign was the bottleneck; on this fixture it isn't.

### What was NOT measured

- **Production-scale paired timing on real chunk data.** I scoped this
  down once B1 revealed the bug — the absolute lazy number from the spec
  (and any version of it) is meaningless until the always-apply fix lands.
  Re-run with the fix on a 30k-read chr-contiguous subsample, both paths,
  same input, same commit. The blocking dependency is having corrected
  TSVs on disk for the subsample. wt_tfiiib_rep1 chunk_000 corrected TSVs
  on H2 are empty (failed run, per
  `HANDOFF_2026-05-20_parallel_bam_writer_abort.md`); ysh1 chunk_000 was
  not surveyed in this session. Sherlock has non-empty TSVs under
  `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/constructs/mdj1/processed_data/rectify/{r15_mdj1_rep1,r15_mdj1_rep2,r15_mdj1_rep3,wtaa_mdj1_rep1,wtaa_mdj1_rep3}/`
  but most appear to be `per_aligner_corrected/minimap2/` only — a true
  multi-aligner winner-flip comparison needs ≥2 aligners with completed
  TSVs at the same commit. Wider survey of the larsms tree warranted.
- **Targeted-realign 35–48% reduction claim** in the spec. Fixture has
  13 eligible short exon blocks; ysh1 has 184k per aligner. Fixture is
  too small to validate the magnitude. Targeted pre-check fires correctly
  (1/13 blocks realigned vs 13 under untargeted = 92% reduction at this
  scale).
- **Full mex67aa_rep1 sidecar reconstruction round-trip.** A 400-record
  smoke confirms the formula on real chunk samples; the 7.9M-record full
  run was not executed.

---

## Recommended order of operations

1. **Land the P0 fix** (delete fast-path lines 310-314, remove dead
   `_correction_requires_transient_edit` once unused). Re-run
   `tests/test_bam_writer_parallel_smoke.py` to confirm it passes (it
   did pass after the patch in this audit).
2. **Add the missing unit tests** Stage 1 acceptance demands: two short
   exon blocks, only one with X-at-HP, assert only that block changes.
3. **Strengthen the equivalence test** so the two-aligner path is
   genuinely two distinct BAMs/TSVs, not the same fixture aliased twice.
4. **Land the P1 sidecar fix** so the silent fallback becomes a loud
   error. Install pyarrow in the H2 rectify env.
5. **Re-run Phase C paired timing** on a 30k-read chr-contiguous slice
   of a real chunk (wt_tfiiib_rep1 chunk_000 chrIV has 54k reads in
   minimap2; chrII has 26k — both viable once a corrected TSV exists).
   Report the real lazy-vs-materialized ratio.
6. **Add the IGV opt-in flag** (P3 diff sketch above).

## Recommended intermediates to make opt-in / eliminate

Ranked by likely saved I/O × certainty of safe removal:

| Artifact | Current | Recommendation |
| - | - | - |
| Per-aligner `corrected.bam` (hard-clip) | opt-in via `args.write_corrected_bam`; chunk script never sets it → already off | Add CLI flag (above) for explicit IGV opt-in |
| Per-aligner `softclipped.bam` | opt-in via separate config; chunk script never sets → off | Same flag pattern |
| Per-aligner `rectified_pA_tail_trimmed.bam` | opt-in via `output_bam` config | Confirm not set by chunk script; already correct |
| `corrected_consensus.unsorted.bam` (disk spill) | always written then sorted in-place | Pipe to `samtools sort -` instead of spilling — saves a write |
| `corrected_polya.unsorted.bam` | same as above when PARQUET set | Same |
| WT untracked validation BAMs (`*.pA_rest.bam`, `*.softclipped.bam`, etc.) | written by `scripts/validation_data/*.py` only | Not pipeline outputs; out of scope |
| `corrected_reads_index.bed.gz`, `minimap2_unrectified.bam`, `corrected_reads_stats.tsv` (WT untracked) | no in-tree writer | Manual debug artifacts; clean up |
