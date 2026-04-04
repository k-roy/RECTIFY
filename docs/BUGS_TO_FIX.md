# RECTIFY Bugs to Fix

## Last Updated: 2026-04-01

---

## Open

*(none)*

---

## Fixed

### ~~Bug 27 — `detect_partial_junction_crossings()` TypeError on SEQ=* reads~~ — Fixed 2026-04-03 (v2.7.5)
`terminal_exon_refiner.py`: `len(clip_seq)` crashed with `TypeError` when `five_prime_clip['sequence']` is `None` (unmapped reads with `SEQ=*` in BAM). Added `if clip_seq is None: continue` guard.
**File:** `rectify/core/terminal_exon_refiner.py`

---

### ~~Bug 26 — `generate_bedgraphs()` KeyError~~ — Fixed 2026-04-03 (v2.7.5)
`analyze_command.py`: Fallback column name `'position'` (non-existent) replaced with `'corrected_position'` (the actual column name in output TSVs).
**File:** `rectify/core/analyze_command.py`

---

### ~~Bug 24 — `_run_junction_aggregation()` KeyError + `.gff.gz` detection~~ — Fixed 2026-04-03 (v2.7.5)
`run_command.py`: Dict key `partial_results['summary']` replaced with `partial_results['stats']` (correct key from `detect_partial_junction_crossings()`). Also fixed `.suffix` → `.suffixes` so `.gff.gz` files are correctly detected as GFF.
**File:** `rectify/core/run_command.py`

---

### ~~Bug 23 — `deconvolve_region()` silent bad output~~ — Fixed 2026-04-03 (v2.7.5)
`netseq_deconvolution.py`: Inverted coordinates (`region_start ≥ region_end`) now raise `ValueError` instead of silently returning a zero-length result.
**File:** `rectify/core/netseq_deconvolution.py`

---

### ~~Bug 22 — `get_netseq_3prime_position()` crash on unmapped reads~~ — Fixed 2026-04-03 (v2.7.5)
`netseq_bam_processor.py`: Added `if read.is_unmapped or read.reference_end is None: raise ValueError(...)` guard to prevent `TypeError` on unmapped reads.
**File:** `rectify/core/netseq_bam_processor.py`

---

### ~~Bugs 20–21 — Silent fallbacks in `_stitch_group()`~~ — Fixed 2026-04-03 (v2.7.5)
`mpb_split_reads.py`: Added `logger.debug()` calls for (20) all-unmapped read groups and (21) cross-chromosome/cross-strand chunk fallbacks, so these events are visible at DEBUG level instead of silently dropped.
**File:** `rectify/core/mpb_split_reads.py`

---

### ~~Bug 19 — BAM file descriptor leak in `refine_terminal_exons()`~~ — Fixed 2026-04-03 (v2.7.5)
`terminal_exon_refiner.py`: Converted from `bam.open()/bam.close()` to `with pysam.AlignmentFile(...) as bam:` context manager, ensuring file handle is closed even on early exit.
**File:** `rectify/core/terminal_exon_refiner.py`

---

### ~~Bugs 17+11 — OOM and chrom normalization in `junction_validator.py`~~ — Fixed 2026-04-03 (v2.7.5)
`junction_validator.py`: `filter_cross_sample_junctions()` replaced `pd.concat` + `itertuples()` (OOM on large datasets) with streaming dict accumulation. `read.reference_name` now passes through `standardize_chrom_name()` before junction key construction.
**File:** `rectify/core/junction_validator.py`

---

### ~~Bug 15 — Negative `ambiguity_min` for reads near chromosome start~~ — Fixed 2026-04-03 (v2.7.5)
`bam_processor.py`: `ambiguity_min` is now clipped to `max(0, ...)` to prevent negative genomic coordinates for reads near the chromosome start.
**File:** `rectify/core/bam_processor.py`

---

### ~~Bug 7 — Missing timeouts on minimap2/mapPacBio/gapmm2~~ — Fixed 2026-04-03 (v2.7.5)
`multi_aligner.py`: `ALIGNER_TIMEOUT` (7200s) now applied to all five aligners (was only applied to uLTRA and deSALT). Added `try/except TimeoutExpired` with process kill for minimap2 pipe.
**File:** `rectify/core/multi_aligner.py`

---

### ~~Missing `logger` in `false_junction_filter.py`~~ — Fixed 2026-04-03 (v2.7.6)
`false_junction_filter.py`: `logger.debug()` was called on line 286 but `import logging` and `logger = logging.getLogger(__name__)` were never added. Would have raised `NameError` at runtime whenever a junction motif analysis exception occurred.
**File:** `rectify/core/false_junction_filter.py`

---

### ~~`NetSeqSignal` not picklable~~ — Fixed 2026-04-03 (v2.7.6)
`netseq_refiner.py`: Added `__getstate__()` and `__setstate__()` methods so `NetSeqSignal` can be pickled for use with `multiprocessing.Pool`. The methods exclude the `threading.Lock` and open `pyBigWig` file handles (which are not picklable) and recreate them on restore.
**File:** `rectify/core/netseq_refiner.py`

---

### ~~BAM file handle leak in `iter_netseq_reads()`~~ — Fixed 2026-04-03 (v2.7.6)
`netseq_bam_processor.py`: The read iteration loop in `iter_netseq_reads()` called `bam.close()` after the loop, but an early `break` (e.g. `max_reads` limit reached) or exception would skip the close. Wrapped the loop in `try: ... finally: bam.close()`.
**File:** `rectify/core/netseq_bam_processor.py`

---

### ~~Bug 1: SettingWithCopyWarning in clustering.py~~ — Fixed 2026-04-01

`clustering.py:473` warned about setting on a copy of a DataFrame slice.
Added explicit `.copy()` before the `_effective_count` assignment block.
**File:** `rectify/core/analyze/clustering.py`

---

### ~~Bug 2: Sample column naming inconsistency~~ — Fixed 2026-04-01

`_load_large_file_chunked` in `analyze_command.py` lacked the auto-detection
fallback that the small-file path already had. Now both paths try
`["sample", "replicate", "sample_id", "sample_name"]` and use the first match.
**Files:** `rectify/core/analyze_command.py`
**Tests:** `tests/test_sample_column_autodetect.py` (7 tests)

---

### ~~Bug 3: gapmm2 PAF→BAM drops read sequences~~ — Fixed 2026-04-01

`_paf_to_bam()` in `multi_aligner.py` never set `seg.query_sequence` because
PAF format omits the read sequence. Fixed by passing the source FASTQ dict
to the converter; sequence + quality are looked up by read_id and copied onto
every gapmm2 BAM record. `_process_and_write_batch` in `consensus.py` provides
the FASTQ dict; logs a warning rather than crashing when no sequence is found.
**Files:** `rectify/core/multi_aligner.py`, `rectify/core/consensus.py`
**Tests:** `tests/test_gapmm2_seq_restore.py`

---

### ~~Bug 4: Duplicate primary records in consensus BAM~~ — Fixed 2026-04-01

Two root causes:
1. `_paf_to_bam()` never checked the `tp:A:` PAF tag — `tp:A:S` secondary
   records got written as primaries. Now `tp:A:S` → `FLAG |= 0x100`.
2. The consensus BAM writer did not track which read_ids already had a primary
   written. Fixed to enforce one-primary-per-read-ID.
**Files:** `rectify/core/multi_aligner.py`, `rectify/core/consensus.py`
**Tests:** `tests/test_no_duplicate_primaries.py` (3 tests)

---

### ~~Bug 5: High-confidence huge introns~~ — Fixed 2026-04-01

Implemented COMPASS-style three-pass cross-sample junction validation:
- **Pass 1** `extract_sample_junctions()`: per-sample intron extraction with
  GT-AG / AT-AC splice motif detection
- **Pass 2** `filter_cross_sample_junctions()`: aggregate across samples, filter
  by min_samples ≥ 2, min_reads ≥ 3, max_intron ≤ 10,000 bp, canonical motifs
- **Pass 3** `apply_junction_filter()`: stream per-sample BAM, downgrade XC
  confidence for reads with unvalidated junctions (high→medium, medium→low)
**New file:** `rectify/core/junction_validator.py`
**Tests:** `tests/test_junction_validator.py` (44 tests)

---

### ~~Bug 6: XR flag inversion edge case~~ — Fixed 2026-04-01

Investigation found `get_softclip_lengths()` in `consensus.py` already handles
strand correctly: trailing clip used for minus-strand 5' end. The apparent
XR inversion in rna15_rep3 was a false positive — the function was correct.
Added tests to document and lock in the expected behavior for both strands.
**Tests:** `tests/test_xr_flag.py`

---

### ~~Chimeric assembly invalid CIGARs~~ — Fixed 2026-04-01

`build_chimeric_cigar()` produced phantom large insertions (e.g., `947I`, `557I`)
and implausible intron skips (e.g., `49426N`) when winning segments from different
aligners had a reference coordinate gap at the boundary. These malformed records
caused `malloc(): memory corruption` / SIGABRT in pysam workers during
region-boundary fetch calls (reproducibly crashed all chrI sub-regions in ysh1_rep2).

**Root cause:** No reference continuity tracking — gaps were silently absorbed as
phantom I operations instead of N bridges.

**Fix:** `rectify/rectify_beta_deploy/chimeric_consensus.py`
- `build_chimeric_cigar()` tracks `cur_ref`; inserts N bridges for gaps, returns
  `(None, [])` sentinel for reference regressions
- New `_validate_chimeric_cigar()`: rejects I > `read_len // 4` or N > `max_intron`
  (default 10,000 bp)
- `select_best_chimeric()` falls back to `_fallback_simple_selection` on either failure
- 15 unit tests in `rectify/tests/test_chimeric_consensus.py` — all passing

---

### ~~Parallel BAM processing crashes~~ — Fixed 2026-04-01 (v2.7.1)

Five root causes of `malloc(): memory corruption` / `SIGABRT` in
`process_bam_file_parallel()` diagnosed and fixed. See CHANGELOG [2.7.1] for full details.

- Per-task pickle flood (~17 GB for 235 regions) → worker initializer pattern
- Repeated BAM open/close cycles in workers → persistent `_WORKER_BAM` handle
- pysam called in main process before pool launch → `_count_reads_and_get_regions()` spawn subprocess
- Monolithic chromosome regions → `max_region_size` 500k → 100k + force-split fallback
- Unguarded `BrokenProcessPool` propagation → per-region `spawn`-context retry with 3 attempts

---

### ~~GFF3 parser improvements~~ — Fixed 2026-03-18

`_parse_gtf()` now handles both GTF and GFF3 formats:
- GFF3: Extracts `ID` → gene_id, `gene` → common name, `Name` → systematic name
- GTF: Extracts `gene_id` and `gene_name`

**Location:** `rectify/core/analyze_command.py:707-770`
