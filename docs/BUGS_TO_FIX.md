# RECTIFY Bugs to Fix

## Last Updated: 2026-04-01

---

## Open

*(none)*

---

## Fixed

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
