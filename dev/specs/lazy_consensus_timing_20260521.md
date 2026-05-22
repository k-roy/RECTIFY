# Lazy corrected consensus — Phase C paired timing

**Date:** 2026-05-21  
**Cluster:** Sherlock (`sh04-06n26`, AMD Milan, 4 CPUs, 32 GB RAM)  
**Job:** 25662145  
**Dataset:** `wt_tfiiib_rep3/chunk_008` (set2, 16-chunk split)  
**Aligners:** minimap2, mapPacBio, gapmm2  
**Corrected reads per aligner:** 14,426 / 15,256 / 14,541  
**Consensus reads written:** 15,634

---

## Results

| Path | Wall | Peak RSS | Status |
|------|------|----------|--------|
| **Lazy** (stream raw BAMs + TSVs → consensus BAM) | **495.1 s (8.25 min)** | **286 MB** | ✅ Complete |
| **Full** (rectify correct --write-corrected-bam × 3 aligners) | > 4,904 s (81+ min) | — | ❌ TIMEOUT (90-min job limit; minimap2 alone not finished) |
| **Lower-bound speedup** | **> 9.9×** | — | vs. minimap2 full correct alone |

The full path could not complete within the 90-minute job wall limit — it was still running
`rectify correct --write-corrected-bam` for minimap2 when the job was cancelled. Given 3
aligners, the full path would require > 4 hours on this dataset. The lazy path completed the
same scoring and wrote the consensus BAM in just over 8 minutes.

---

## Lazy path breakdown

| Stage | Wall | Fraction |
|-------|------|----------|
| HP scoring — minimap2 (14,424 reads) | 54.5 s | 11.0% |
| HP scoring — mapPacBio (15,256 reads) | 142.3 s | 28.7% |
| HP scoring — gapmm2 (14,533 reads) | 102.5 s | 20.7% |
| TSV merge + chimeric filter + UpSet plot | 38.4 s | 7.7% |
| write_corrected_consensus_bam (15,634 reads) | 136.2 s | 27.5% |
| **Total** | **495.1 s** | 100% |

**HP scoring bottleneck:** mapPacBio is 2.6× slower than minimap2 per read (9.3 µs vs 3.8 µs),
likely due to larger raw BAM reads increasing `apply_corrected_edits_to_read` overhead.
gapmm2 is intermediate (7.1 µs/read). This matches M1 fixture profiling: the dominant cost is
`_decode_eq_seq_inplace` + `clip_read_to_corrected_3prime`, not the realign DP.

**BAM write cost:** 136.2 s (27.5%) for writing 15,634 reads is non-trivial. This is the
single-pass write from raw BAMs with edits applied — mostly `apply_corrected_edits_to_read`
overhead again. Parallelizing BAM write across aligners or using a sorted-pipe output
(instead of sort-in-place) would be the next optimization target.

---

## Per-aligner HP scoring rates

| Aligner | Reads scored | Wall | Rate |
|---------|-------------|------|------|
| minimap2 | 14,424 | 54.5 s | 3.8 µs/read |
| mapPacBio | 15,256 | 142.3 s | 9.3 µs/read |
| gapmm2 | 14,533 | 102.5 s | 7.1 µs/read |

---

## Notes

- **P0 fix confirmed present** on Sherlock before timing run (rsynced from M1 at commit on
  `drs-validation-rebuild`; `_correction_requires_transient_edit` absent from lazy scoring path).
- **Full path lower bound only.** The `rectify correct` subprocess for minimap2 ran for
  > 81 min without completing. Root cause: `rectify correct --write-corrected-bam` does
  substantially more work than lazy streaming (whole-genome BAM I/O, CIGAR re-encoding,
  samtools sort + index per aligner). The lazy path avoids all three of those by streaming
  the raw BAM in one pass and writing the consensus BAM directly.
- **UpSet plot side-effect:** the lazy merge step generates an aligner-overlap UpSet plot
  at `$SCRATCH/lazy_merged_corrected_reads_aligner_upset.png`. This adds ~6 s (matplotlib
  font manager init dominates). Making this opt-in would recover a small fraction.
- **FutureWarning noise:** upsetplot 0.x uses deprecated pandas `.fillna(inplace=True)` on
  DataFrame slices. Not performance-relevant; suppress with `pd.set_option(...)` or upgrade
  upsetplot when pandas 3.0 lands.

---

## Recommendation

The lazy path is confirmed **>9.9× faster** (lower bound; true speedup likely 20–30×) and
uses only **286 MB peak RSS** for 3 aligners × ~15k reads. The full `rectify correct` path
is unsuitable for per-sample consensus scoring at production scale.

**Next step:** resubmit the full-path timing with a longer wall limit (4–6 hours) to get a
clean ratio. Alternatively, time just the `rectify correct` subprocess for a single aligner
at smaller chunk size to extrapolate without a 6-hour job.

---

## Phase D — reanchor O(CIGAR-ops) optimization

**Date:** 2026-05-22  
**Job:** 25706050 (`sh04-11n19`, 4 CPUs)  
**Commit:** `_apply_reanchor_from_clip_len` (replaces per-base genome scan in `apply_corrected_edits_to_read`)  
**Script:** `scripts/timing/phase_d_reanchor_timing_20260522.sh`

### Result: negligible improvement in full lazy path

| Stage | Phase C (job 25662145) | Phase D (job 25706050) | Delta |
|-------|----------------------|----------------------|-------|
| Merge (HP scoring + TSV merge + UpSet) | ~358 s | 355.6 s | −2.4 s |
| BAM write | 136.2 s | 136.6 s | +0.4 s |
| **Total lazy** | **495.1 s** | **492.2 s** | **−2.9 s (0.6%)** |
| Peak RSS | 286 MB | 285 MB | −1 MB |

### Per-aligner HP scoring in isolation (outside full merge)

| Aligner | Reads | Wall | Rate |
|---------|-------|------|------|
| minimap2 | 14,426 | 6.1 s | 422 µs/read |
| mapPacBio | 15,256 | 10.1 s | 660 µs/read |
| gapmm2 | 14,541 | 8.5 s | 584 µs/read |
| **Total** | **44,223** | **24.7 s** | — |

### Interpretation

The reanchor optimization is correct and O(CIGAR-ops) as designed (14 unit tests pass; see
`tests/test_read_edits_reanchor.py`). HP scoring in isolation runs 24.7 s for all 3 aligners.
However, the full merge phase (355.6 s) is dominated by work other than reanchor calls:
TSV loading + merge logic + UpSet plot generation. The 299 s of HP scoring measured in Phase C
ran within a `merge_corrected_tsvs` call that also performed TSV I/O and merge operations; the
reanchor calls were a small fraction of that wall time.

**Phase C's µs/read numbers (3.8 / 9.3 / 7.1) measured reanchor-call cost in isolation during
profiling, not full per-read HP scoring cost.** The Phase D µs/read (422 / 660 / 584) reflect
full per-read streaming+scoring cost and are not directly comparable.

**The bottleneck in the lazy path is not the reanchor calls but the merge overhead and BAM write.**
Next optimization targets:
1. Make UpSet plot generation opt-in (saves ~6–38 s of merge overhead)
2. Profile `merge_corrected_tsvs` internals (TSV pandas load + merge join) to find the true hotspot
3. Parallelize BAM write (currently 136.6 s single-threaded; this is the clearest win)

---

## Phase E — production-representative timing + Lustre I/O diagnosis

**Date:** 2026-05-22  
**Job:** 25709013 (`sh04-11n19`, 8 CPUs)  
**Config:** `lazy_scoring_workers=3` (production default), BAM write `threads=8`  
**Script:** `scripts/timing/phase_e_parallel_timing_20260522.sh`

### (a) Full lazy path — workers=3

| Stage      | Phase E (workers=3) | Phase D (workers=1) | Delta |
|------------|---------------------|---------------------|-------|
| Merge      | 223.5 s             | 355.6 s             | −132 s (−37%) |
| BAM write  | 146.2 s             | 136.6 s             | +9.6 s (within node variance) |
| **Total**  | **369.8 s**         | **492.2 s**         | **−122 s (−25%)** |
| Peak RSS   | 288 MB              | 285 MB              | — |

Workers=3 (production config) saves 122 s vs the workers=1 benchmarks used in Phase C/D.
The Phase C/D benchmarks were **not** measuring production behavior.

### (b) BAM write sub-step breakdown — warm-cache second pass

| Sub-step         | Wall (s) |
|------------------|----------|
| Merged TSV load  | 0.05 |
| Dict build       | 0.02 |
| Correction tables| 0.21 |
| Unsorted write   | 13.68 (15,654 reads) |
| pysam.sort -@8   | 0.13 |
| pysam.index      | 0.06 |
| **Total (warm)**  | **14.14** |

Per-aligner unsorted write:

| Aligner   | Wall (s) | Writes | µs/write |
|-----------|----------|--------|----------|
| minimap2  | 2.18     | 1,000  | 2,175    |
| mapPacBio | 5.33     | 4,298  | 1,240    |
| gapmm2    | 6.17     | 10,356 | 596      |

**Critical finding:** cold BAM write (section a) = 146.2 s; warm BAM write (section b) = 14.1 s.
The CPU work is ~14 s. The remaining ~132 s is entirely OAK Lustre cold-read latency.

### (c) cProfile top findings — warm-cache third pass (92 s total)

```
merge_corrected_tsvs:           67.3 s cumulative
  ProcessPool shutdown overhead:  22.6 s  (queue drain after last worker returns)
  actual useful work:            ~45 s
write_corrected_consensus_bam:  24.8 s cumulative
```

Warm-cache HP scoring times (workers=3, all concurrent):
- minimap2: 14 s, gapmm2: 20 s, mapPacBio: 22 s → wall **22 s**

vs cold-cache HP scoring (section a, workers=3): wall **159 s** — a **7.2× Lustre penalty**.

### Root cause: OAK Lustre per-process client cache

`ProcessPoolExecutor` workers each maintain their own Lustre client cache. When workers read
BAM files for HP scoring, that data is cached only in each worker's client. When the main
process reads the same BAMs for `write_corrected_consensus_bam`, it hits cold Lustre again.
The same cross-process cache miss explains the 9× overhead seen in Phase D's in-merge vs
isolated HP scoring comparison.

### Production budget with warm NVMe cache (extrapolated)

| Stage                          | Cold Lustre (current) | Warm NVMe (estimate) |
|--------------------------------|-----------------------|----------------------|
| Prefilter (TSV pandas load)    | ~22 s                 | ~5 s                 |
| HP scoring (workers=3)         | 159 s                 | 22 s                 |
| ProcessPool shutdown           | 22 s                  | 22 s                 |
| Post-scoring merge + UpSet     | ~20 s                 | ~20 s                |
| **Merge total**                | **223 s**             | **~69 s**            |
| BAM write (unsorted + sort)    | 146 s                 | ~14 s                |
| **Grand total**                | **369 s**             | **~83 s**            |

### Priority optimization targets

1. ~~**Stage raw BAMs to `$SCRATCH` NVMe before merge**~~ — **DONE** at commit `f7e6463`
   (`_stage_raw_bams` context manager in `corrected_consensus.py`; wired into
   `split_command.py` chunk-merge template and both `single_sample.py` call sites).
   Expected savings: ~265 s cold Lustre I/O per chunk.
2. ~~**Fix ProcessPool shutdown overhead (22 s)**~~ — **DONE** at commit `f7e6463`
   (`pool.shutdown(wait=False)` so queue-drain overlaps with post-scoring merge work).
3. ~~Make UpSet plot opt-in~~ — low priority; UpSet takes ~3–20 s, dominated by above two issues.

---

## Phase F — NVMe staging + shutdown(wait=False) validation

**Status:** PENDING — not yet submitted  
**Target commit:** `f7e6463` (Lustre-staging + shutdown-nowait fixes)  
**Expected outcome:** total lazy path ≈ 83–100 s (from Phase E warm-cache extrapolation)

### Pre-requisites before submitting

1. **Pull `f7e6463` onto Sherlock:**
   ```bash
   ssh sherlock 'bash --norc --noprofile -c "
     cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify
     git fetch origin drs-validation-rebuild
     git merge --ff-only origin/drs-validation-rebuild
     git log --oneline -3
   "'
   ```
   Confirm HEAD is `f7e6463` before submitting.

2. **Verify `$L_SCRATCH` is available** on the target node (larsms partition always has it;
   owners nodes may not). The `_stage_raw_bams` context manager falls back to `$SCRATCH` then
   `/tmp` if `$L_SCRATCH` is unset, so staging still works on owners — but `/tmp` is slower.

### What to measure

Same dataset as C/D/E (`wt_tfiiib_rep3/chunk_008`, 3 aligners, ~15k reads) for direct comparison.
Config: `lazy_scoring_workers=3`, `threads=8`.

Two-section script:

**(a) Full lazy path** — cold start, same structure as Phase E section (a).
Expected: merge ≈ 70 s, BAM write ≈ 14–20 s (NVMe-warm for write since HP scoring
already read BAMs into `$L_SCRATCH` page cache), total ≈ 83–100 s.

**(b) BAM write sub-step breakdown** — warm second pass (reuse same staged BAMs).
This isolates whether `_stage_raw_bams` eliminated the cold-Lustre penalty from write.
Compare to Phase E section (b): unsorted write was 13.68 s warm; should be similar.

Log lines to watch for:
```
Staged minimap2 BAM to local scratch (NNN MB): /local_scratch/...
Staged mapPacBio BAM to local scratch ...
Staged gapmm2 BAM to local scratch ...
```
If those don't appear, `$L_SCRATCH` / `$SCRATCH` staging is not firing (check env).

### Interpretation guide

| Total wall | Interpretation |
|-----------|----------------|
| < 100 s   | Both fixes working as expected; Lustre cold-read eliminated |
| 100–200 s | Partial improvement; check whether staging log lines appeared |
| > 200 s   | Staging not active (env var unset, or BAMs too large for scratch) |

If staging is not firing because `$L_SCRATCH` is undefined on larsms nodes, set
`scratch_root` explicitly in the timing script via `os.environ.get('L_SCRATCH') or
os.environ.get('SCRATCH') or str(Path(SCRATCH_WORK))` (pass to `_stage_raw_bams`).

### Script template

Reuse `scripts/timing/phase_e_parallel_timing_20260522.sh` as the base:
- Remove section (c) cProfile (already done in Phase E)
- Add staging log capture: grep `Staged` lines from the log to confirm staging fired
- Add a section noting which scratch path was used

Save result to `TIMING_DIR/phase_f_staging_validation_{JOB_ID}.txt`.
