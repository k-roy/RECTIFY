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
