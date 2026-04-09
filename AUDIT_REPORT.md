# RECTIFY Codebase Audit Report

**Date:** 2026-04-08
**Version audited:** 2.7.7
**Scope:** Full codebase — 16 parallel agents covering all major modules
**Independent verification (round 1):** 2026-04-08 — 3 parallel checker agents; 2 findings REFUTED (C1, H1), 1 strengthened (H5), 1 nuanced (C3).
**Independent verification (round 2):** 2026-04-08 — 8 parallel checker agents re-read all remaining flagged lines; 1 additional finding REFUTED (C4), H5 narrowed to single call site, M5 partially already fixed. Net confirmed: **2 CRITICAL, 5 HIGH, 18 MEDIUM**.

---

## Executive Summary

16 agents audited all major modules in parallel. After two rounds of independent code verification, **2 CRITICAL bugs** and **5 HIGH-risk issues** were confirmed (3 initial findings were false positives). ~18 MEDIUM findings also confirmed. The most impactful bugs are a 2 bp tolerance error in 3'SS detection and the SLURM profile system silently ignoring `use_scratch`/`streaming`.

---

## CRITICAL

### ~~C1 — Minus-strand poly-A walk-back direction inverted~~ ✗ REFUTED
**File:** `rectify/core/indel_corrector.py:192-252`
**Agent:** 3 (3' Indel Correction) | **Checker verdict: REFUTED**

~~The `find_polya_boundary()` function correctly walks **backward** for plus-strand reads (finding the first non-A position moving left). For minus-strand reads it walks **forward**, returning a higher genomic position than `reference_start`. But minus-strand 3' ends are at `reference_start` (lowest coord); corrections must move **leftward** (lower coords). The current code shifts rightward — wrong direction.~~

**Checker analysis:** Walking FORWARD for minus-strand is correct. The minus-strand 3' end is at `reference_start` (leftmost coord). Genomic T-runs to the LEFT of `reference_start` can pull the apparent 3' end leftward; the true CPA is therefore at HIGHER genomic coordinates. Walking forward (rightward) through `aligned_positions` to find the first non-T position is the correct direction. The code comment ("Walk FORWARD from 3' end") is accurate. No bug.

---

### C2 — Double-subtraction of 3'SS position in terminal exon refiner
**File:** `rectify/core/terminal_exon_refiner.py:1072`
**Agent:** 8 (Terminal Exon Refiner)

`splice_index.three_ss` stores positions as `intron_end - 2` (the 'A' of the AG dinucleotide). Line 1072 then subtracts another `- 2`:

```python
if abs(read_end - pos - 2) <= tolerance and not has_junction:
```

Since `pos` is already `intron_end - 2`, this is equivalent to `abs(read_end - intron_end + 4)` instead of the intended `abs(read_end - intron_end + 2)`. The two-base double-subtraction may be numerically coincidental for some cases but is logically wrong and will misidentify truncated reads.

**Fix:** Remove the `- 2` from line 1072: `abs(read_end - pos) <= tolerance`.

---

### C3 — `use_scratch` and `streaming` silently dropped by `_apply_profile()`
**File:** `rectify/core/batch_command.py:452-477`
**Agent:** 12 (SLURM & HPC Integration) | **Checker verdict: CONFIRMED (nuanced)**

`_apply_profile()` merges YAML profile values into argparse args, but `use_scratch` and `streaming` are absent from the `profile_to_arg` mapping. Both flags bypass `_apply_profile()` entirely — they are handled in other functions via `getattr(args, 'use_scratch', False)` and `getattr(args, 'streaming', False)`. This means profile YAML values for these keys are silently ignored; only the hardcoded `False` defaults apply. Users who set `use_scratch: true` / `streaming: true` in their profile (as `sherlock_larsms.yaml` does) get neither scratch staging nor streaming. The result:
- Direct NFS I/O instead of high-bandwidth scratch (2–3× slower)
- Full BAM loaded into RAM (OOM risk for >3 GB BAMs)

The same function also silently drops `analyze_cpus`, `analyze_mem`, and `analyze_time` (see H3 below).

**Fix:** Add missing keys to `_apply_profile()`'s `profile_to_arg` dict and `defaults` dict:
```python
'use_scratch': 'use_scratch',
'streaming':   'streaming',
'analyze_cpus': 'analyze_cpus',
'analyze_mem':  'analyze_mem',
'analyze_time': 'analyze_time',
```

---

### ~~C4 — Unhandled `None` return from position functions in BAM processor~~ ✗ REFUTED
**File:** `rectify/core/bam_processor.py:191,198`
**Agent:** 2 (BAM Processing) | **Checker verdict: REFUTED**

~~`get_read_3prime_position()` can return `(None, strand)` for unmapped reads; `get_read_5prime_position()` can return `None`. At lines 191 and 198, these are unpacked and used without null checks. Lines 253–255, 285, 289–290 subsequently use these values as coordinates, which crashes atract_detector and ag_mispriming modules.~~

**Checker analysis:** All three code paths that call into `correct_read_3prime()` are guarded by `if read.is_unmapped: continue` at the loop level — lines 511–512, 652–654, and 979–980. Unmapped reads never reach the position functions. No fix needed.

---

## HIGH

### ~~H1 — Aligner-control flags undeclared in `run-all` parser~~ ✗ REFUTED
**File:** `rectify/core/run_command.py:591-603`
**Agent:** 1 (CLI & Command Wiring) | **Checker verdict: REFUTED**

~~None of these are declared in the `run-all` subparser in `cli.py`. They silently default to `None`/`False`, so passing `rectify run-all --parallel-aligners` has no effect.~~

**Checker analysis:** All 6 flags are properly declared in the `run-all` subparser in `cli.py`: `--parallel-aligners` (add_argument with `action='store_true'`), `--junction-aligners` (`nargs='+'`, `choices=['uLTRA', 'deSALT']`), `--chimeric-consensus`, `--ultra-path`, `--desalt-path`, `--mapPacBio-chunks` (`type=int`). No fix needed.

---

### H2 — Minus-strand 5'SS truncation check uses wrong coordinate
**File:** `rectify/core/terminal_exon_refiner.py:1085-1121`
**Agent:** 8 (Terminal Exon Refiner)

`detect_junction_truncated_reads()` checks `read_start` against the 5'SS position for minus-strand reads. But on minus strand, the 5' end is at `reference_end - 1` (highest coord), not `reference_start`. The check:

```python
if abs(read_start - pos) <= tolerance and not has_junction:
```

compares the 3' end (low coord) against the 5'SS (high coord). These are on opposite ends of the read; the detection will almost never trigger.

**Fix:** Use `read_end - 1` instead of `read_start` for the minus-strand 5'SS check.

---

### H3 — `analyze_cpus/mem/time` not applied from profiles
**File:** `rectify/core/batch_command.py:452-477`
**Agent:** 12 (SLURM & HPC Integration)

Same root cause as C3. The three analyze-stage resource parameters defined in all profile YAMLs are not transferred to args; the generated script uses hardcoded defaults (`8` CPUs, `64G` mem, `8:00:00` time) regardless of what the profile specifies.

**Fix:** Covered by the C3 fix (add `analyze_cpus`, `analyze_mem`, `analyze_time` to `_apply_profile()`).

---

### H4 — GFF-derived 3'SS positions not validated against genome
**File:** `rectify/core/terminal_exon_refiner.py:282-290`
**Agent:** 8 (Terminal Exon Refiner)

Splice site positions computed from GFF (`intron_end - 2`) are stored and used without verifying that the genomic sequence at those positions actually contains the expected AG dinucleotide. Off-by-one errors in the source annotation silently propagate.

**Fix:** After loading, fetch and assert `genome[chrom][pos:pos+2].upper() == 'AG'` (minus strand: reverse-complement). Log a warning and skip sites that fail.

---

### H5 — `run_command.py` doesn't pass config to provenance
**File:** `rectify/core/run_command.py:531`
**Agent:** 13 (Configuration & Provenance) | **Checker verdict: NARROWED**

Round-2 cross-check found that `correct_command.py:279` already calls `init_provenance(..., config=config)` and `batch_command.py:774` calls `init_provenance(..., config=vars(args))`. Only `run_command.py:531` omits the config argument, recording no parameter values for `run-all` runs.

**Fix:** Change `run_command.py:531` to `init_provenance(output_dir, description="...", config=vars(args))`.

---

### H6 — Zero unit tests for `terminal_exon_refiner.py`
**File:** `rectify/core/terminal_exon_refiner.py`
**Agent:** 14 (Test Coverage)

The module (1690 lines, multiple coordinate-sensitive code paths) has no dedicated test file. Given the CRITICAL and HIGH findings above (C2, H2, H4), this is the highest-priority test gap.

---

### H7 — `consensus.py` alignment selection only indirectly tested
**File:** `rectify/core/consensus.py`
**Agent:** 14 (Test Coverage)

Core per-read consensus selection logic is only exercised through two peripheral tests (`test_xr_flag.py`, `test_gapmm2_seq_restore.py`). Tie-breaking, penalty scoring, and multi-aligner merging have no direct unit tests.

---

## MEDIUM

### M1 — GC-AG missing from canonical motif set
**File:** `rectify/core/junction_validator.py:65`
**Agent:** 7 (Junction Validation)

```python
_CANONICAL_MOTIFS = frozenset({'GT-AG', 'AT-AC'})
```

`splice_motif.py` correctly lists GC-AG as canonical (~0.5–1% of all junctions), but the validator omits it. With `require_canonical=True`, real GC-AG junctions are silently filtered.

**Fix:** `_CANONICAL_MOTIFS = frozenset({'GT-AG', 'AT-AC', 'GC-AG'})`

---

### M2 — Position shift double-counting for A-tract reads
**File:** `rectify/core/processing_stats.py:152-168`
**Agent:** 13 (Configuration & Provenance)

`total_position_shifts` is incremented twice for reads with A-tract ambiguity resolution: once at line 154 (atract_ambiguity branch) and again at line 168 (any position change). Stats output overstates shifts by the count of A-tract-corrected reads.

**Fix:** Gate the line-168 increment on `'atract_ambiguity' not in corrections`, or use a single increment point.

---

### M3 — `--polya-model` trained model never used
**File:** `rectify/core/polya_model.py` + `correct_command.py:173-175`
**Agent:** 4 (Poly(A) Trimming & Model)

The `--polya-model` CLI flag, model training pipeline, and `load_model()` infrastructure exist, but `polya_model_path` is captured in config and then discarded. The hardcoded 80% A-richness threshold from `config.py` is always applied. Users who train a model and pass it via `--polya-model` will silently see no effect.

**Fix:** Either wire `polya_model_path` through `bam_processor.correct_read_3prime()` and use the trained threshold, or remove the flag and model loading code with a deprecation notice.

---

### M4 — Jensen-Shannon divergence epsilon regularization biased
**File:** `rectify/core/analyze/shift_analysis.py:158-179`
**Agent:** 9 (APA Detection & Clustering)

The JSD implementation adds `1e-10` uniformly before normalization, which biases estimates for genes with sparse counts. The correct approach normalizes first.

**Fix:** Replace the custom implementation with `scipy.spatial.distance.jensenshannon(p, q)`.

---

### M5 — Some APA clustering parameters not CLI-configurable
**File:** `rectify/core/analyze/clustering.py:26-31`
**Agent:** 9 (APA Detection & Clustering) | **Checker verdict: PARTIALLY CONFIRMED**

Round-2 cross-check found that `DEFAULT_CLUSTER_DISTANCE` and `DEFAULT_MIN_READS` are already exposed as `--cluster-distance` and `--min-reads` in `analyze_command.py`. The remaining three are not: `DEFAULT_MIN_PEAK_SEPARATION = 5`, `DEFAULT_MAX_CLUSTER_RADIUS = 10`, and `DEFAULT_MIN_SAMPLES = 2`. No signal smoothing before peak calling; single-read outliers in sparse regions can still be called as peaks.

**Fix:** Expose via CLI (`--min-peak-sep`, `--max-cluster-radius`, `--min-cluster-samples`); apply Gaussian smoothing (σ = 2–3 bp) before peak calling.

---

### M6 — A-richness threshold inconsistency
**File:** `rectify/core/false_junction_filter.py:67` vs `rectify/config.py:117`
**Agent:** 7 (Junction Validation)

`false_junction_filter.py` uses `A_RICHNESS_THRESHOLD = 0.7` while `config.py` uses `POLYA_RICHNESS_THRESHOLD = 0.8` for the same concept. Junctions downstream of 70–80% A-rich regions may be incorrectly flagged as poly(A) artifacts.

**Fix:** Import `POLYA_RICHNESS_THRESHOLD` from `config.py` in `false_junction_filter.py`; remove the local constant.

---

### M7 — `gene_attribution.py` silently uses un-converted GFF coords
**File:** `rectify/core/analyze/gene_attribution.py:189-196`
**Agent:** 15 (Coordinate System)

`build_cds_interval_tree()` assumes the input DataFrame is already 0-based, but has no explicit conversion and a comment noting "sometimes they aren't." A raw GFF input produces gene attribution off by 1 at every boundary with no error.

**Fix:** Add an explicit `start = int(start) - 1` step (or assert the DataFrame has already been converted) with a clear error message.

---

### M8 — Minus-strand ambiguity window direction incorrect
**File:** `rectify/core/bam_processor.py:401-407`
**Agent:** 2 (BAM Processing)

When computing the ambiguity window for minus-strand reads, the code adds `ambiguity_range` to `current_position` (moving rightward), but minus-strand corrections move leftward (lower coords). The window should extend in the correction direction.

**Fix:** For minus strand: `ambiguity_min = current_position - result['ambiguity_range']`, `ambiguity_max = current_position`.

---

### M9 — BAM handle not in context manager in `find_coverage_gaps()`
**File:** `rectify/core/bam_processor.py:826-852`
**Agent:** 2 (BAM Processing)

The BAM file is opened but only explicitly closed on the normal path. An exception during the read loop leaks the file descriptor. Under parallel workload this can exhaust system FD limits.

**Fix:** Wrap in `with pysam.AlignmentFile(bam_path, 'rb') as bam:`.

---

### M10 — Bare `except Exception: pass` swallows gene attribution errors
**File:** `rectify/core/bam_processor.py:276-282`
**Agent:** 2 (BAM Processing)

Gene attribution failures are silently caught and discarded. If `gene_attribution` is broken, the user has no indication that all `gene_id` fields will be empty.

**Fix:** Log the exception at WARNING level: `logger.warning("Gene attribution failed: %s", e)`.

---

### M11 — Profile fields not validated at load time
**File:** `rectify/core/batch_command.py:419-449`
**Agent:** 13 (Configuration & Provenance)

YAML profile files are loaded as flat dicts without field validation. A typo (e.g., `cpu:` instead of `cpus:`) silently uses the hardcoded default, generating SLURM scripts with wrong resource requests.

**Fix:** After loading, check for required keys (`partition`, `cpus`, `mem`, `time`) and log a warning for any unrecognized keys.

---

### M12 — `--continue-on-error` not wired to SLURM mode
**File:** `rectify/core/batch_command.py:1095-1098`
**Agent:** 1 (CLI & Command Wiring)

`--continue-on-error` is respected in interactive mode but ignored for SLURM-generated scripts, which always use `set -euo pipefail`. Users passing this flag with `--submit` get a false sense of resilience.

**Fix:** Either remove the flag from the batch parser (mark it interactive-only in the help text) or add error-handling wrappers to generated scripts when the flag is set.

---

### M13 — `--partition` not enforced when `--submit` is used
**File:** `rectify/core/batch_command.py:993-997`
**Agent:** 1 (CLI & Command Wiring)

The `--partition` help text says "required for --submit" but argparse does not enforce this; the script is generated and submitted with a comment instead of a real partition line.

**Fix:** Add a validation check in `_run_slurm_mode()` that raises an error if `--submit` is set and `partition` is `None`.

---

### M14 — Hardcoded `sherlock_larsms.yaml` references in docstrings/examples
**File:** `rectify/cli.py:503`, `rectify/core/batch_command.py:14,905`
**Agent:** 1 (CLI & Command Wiring)

Example invocations in docstrings still reference the lab-specific profile path, which is confusing for users on other clusters.

**Fix:** Replace with `--profile ./my_cluster.yaml` or reference `slurm_generic.yaml`.

---

### M15 — Silent data loss: reads without gene attribution not logged
**File:** `rectify/core/analyze/apa_detection.py:310`
**Agent:** 9 (APA Detection & Clustering)

Reads without a gene assignment are silently dropped. If gene attribution fails upstream, the user may be unaware that most reads are being excluded.

**Fix:** Log a warning if > 10% of reads lack gene attribution.

---

### M16 — Subprocess return code not checked for uLTRA/deSALT SAM→BAM
**File:** `rectify/core/multi_aligner.py:758-778, 914-928`
**Agent:** 5 (Consensus & Multi-Aligner)

`view_proc.returncode` is never checked after the `samtools view | samtools sort` pipeline for uLTRA and deSALT. A failed `samtools view` (e.g., invalid SAM) produces a corrupt or empty BAM that is silently passed to consensus selection.

**Fix:** Check `view_proc.returncode != 0` and `sort_proc.returncode != 0` after `communicate()`; raise `RuntimeError` on failure.

---

### M17 — SGE array task IDs are 1-based; script assumes 0-based
**File:** `rectify/core/batch_command.py:62-70`
**Agent:** 12 (SLURM & HPC Integration)

SLURM and PBS use 0-based array indices; SGE/UGE uses 1-based. The generated script uses `$RECTIFY_TASK_ID` directly to index into sample arrays without normalization, so SGE arrays would be off-by-one for all samples.

**Fix:** Add `if [ -n "$SGE_TASK_ID" ]; then RECTIFY_TASK_ID=$((SGE_TASK_ID - 1)); fi` normalization.

---

### M18 — Exclusion regions: 1-based coordinate detection is a warning, not an error
**File:** `rectify/core/exclusion_regions.py:347-355`
**Agent:** 16 (Chimeric & Special Case Handling)

When position values of `1` are detected (possible 1-based input), only a WARNING is logged; processing continues with the 0-based assumption, silently shifting all exclusion boundaries by 1.

**Fix:** Elevate to an error or add explicit GFF→0-based conversion with `start - 1`.

---

## LOW (selected)

| ID | File | Finding |
|----|------|---------|
| L1 | `bam_processor.py:693-695` | Dead `if True else` branch in `write_position_index()` |
| L2 | `bam_processor.py:321,341,348` | Unused `polya_shift` and `indel_shift` variables |
| L3 | `clustering.py:31` | `DEFAULT_MIN_SAMPLES = 2` defined but never used |
| L4 | `shift_analysis.py:147-148` | `'same_major': a == b` can be True when both are None |
| L5 | `provenance.py` | Timestamp format inconsistent: microseconds in sidecar, seconds-only in header |
| L6 | `multi_aligner.py:940` | Tiebreaker non-deterministic when all criteria tie (insertion-order dependent) |
| L7 | `splice_aware_5prime.py` | `first_exon_start` field has different semantics in junction vs annotation code paths |
| L8 | `cli.py:607-612` | `--no-polya-sequenced` in `run-all` is inverse of `--polya-sequenced` in `correct`; asymmetric naming |

---

## Modules That Passed

| Module | Agent | Result |
|--------|-------|--------|
| Visualization & strand handling | 11 | All CRITICAL strand requirements correct; comprehensive tests |
| DESeq2 & statistics | 10 | Correct matrix orientation, thread limits, PCA normalization; no issues |
| Chimeric & special case handling | 16 | No CRITICAL/HIGH findings; architecture sound |
| Splice motif coordinate fetching | 7 | Donor/acceptor dinucleotide fetching correct for both strands |
| Consensus subprocess safety | 5 | No shell injection risk; all paths use list-form subprocess calls |
| 5' splice rescue core logic | 6 | Soft-clip extraction and coordinate transformations correct |

---

## Disputed Finding (Agent 15 vs Agent 6)

Agent 15 flagged `splice_aware_5prime.py:270` as CRITICAL:
> "Minus strand sets `five_prime_corrected = last_junction_end` — should be `last_junction_start - 1`"

Agent 6 analyzed the same line as correct:
> "Sets corrected position to `last_junction_end` — the first base of the upstream exon (0-based). For minus strand, the intron `[start, end)` has the upstream exon at positions ≥ end; `end` is the correct boundary."

Agent 6's analysis is more thorough and geometrically consistent with the 0-based half-open convention. This is **not** a bug. The upstream exon on minus strand sits at genomic positions ≥ `intron_end`; `intron_end` is the correct 5' end correction target.

---

## Priority Fix Order

*(Refuted findings C1, C4, H1 excluded)*

| Priority | ID | Effort | Impact |
|----------|----|--------|--------|
| P1 | C3 | Small | Profile `use_scratch`/`streaming` silently ignored by `_apply_profile()` |
| P1 | C2 | Small | 3'SS double-subtraction gives 2 bp tolerance error in terminal exon refiner |
| P2 | H2 | Small | Minus-strand 5'SS truncation check always compares wrong end of read |
| P2 | H4 | Medium | GFF-derived 3'SS positions not validated against genome sequence |
| P2 | H5 | Trivial | `run_command.py` missing config in `init_provenance()` call |
| P3 | H6/H7 | Large | Add unit tests for `terminal_exon_refiner` and `consensus` selection |
| P3 | M1 | Trivial | Add GC-AG to `_CANONICAL_MOTIFS` in `junction_validator.py` |
| P3 | M2 | Small | Fix position shift double-counting in `processing_stats.py` |
| P3 | M4 | Small | Replace manual JSD with `scipy.spatial.distance.jensenshannon` |
| P3 | M6 | Trivial | Unify A-richness threshold: import from `config.py` |
| P4 | M3 | Medium | Wire `--polya-model` to correction pipeline or remove with deprecation notice |
| P4 | M5 | Small | Expose remaining 3 APA clustering params as CLI args |
| P4 | M8 | Small | Fix minus-strand ambiguity window direction in `bam_processor.py` |

---

---

## Legion Audit — Round 2 (2026-04-08)

**Scope:** 32 parallel discovery agents + 3-way independent verification per reported bug
**Coverage:** All 93 source files (~62k lines), with A01–A32 agents each covering a specific file slice or functional area
**Method:** Each unique bug received three simultaneous verifier verdicts — V1 (Code Confirmation), V2 (Test Coverage), V3 (Real-World Impact)
**Confirmation rule:** A bug is CONFIRMED if V1 = CONFIRMED and (V2 = NOT_TESTED or V3 ≠ UNREACHABLE)
**Exclusions:** Previously documented findings C1–C4, H1–H7, M1–M17 were explicitly excluded from agent scope
**Net findings:** 10 CRITICAL confirmed, 2 CRITICAL disputed, 18 HIGH confirmed, 1 HIGH disputed, 25 MEDIUM/LOW confirmed, 7 false positives

---

### Confirmed New Findings — CRITICAL

#### NEW-001 — `rescue_softclip_at_homopolymer` / `rescue_mismatch_inside_homopolymer` ignore corrected position; always re-anchor to `reference_end - 1`
**File:** `rectify/core/indel_corrector.py:1402, 1415`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

Both `rescue_softclip_at_homopolymer()` and `rescue_mismatch_inside_homopolymer()` internally compute their walk baseline as `raw_pos = read.reference_end - 1`, discarding any correction already applied by prior correction steps. In contrast, `rescue_polya_prefix_in_softclip()` at line 1388 correctly accepts `current_pos` as a parameter, allowing incremental correction.

When indel correction pipelines call these two functions after an initial 3' boundary correction, the re-anchor to `reference_end - 1` silently undoes the prior correction and applies its own walk starting from the original uncorrected position. The final reported position is wrong for any read where the pre-existing correction moved the position by more than 1 bp.

**Fix:** Add `current_pos` parameter to both functions; replace `raw_pos = read.reference_end - 1` with `raw_pos = current_pos`.

---

#### NEW-002 — Minus-strand exon sequence fetched from wrong genomic window in 5' splice rescue
**File:** `rectify/core/splice_aware_5prime.py:606`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

For minus-strand reads, the exonic sequence upstream of the intron is at genomic positions `[intron_end - rescue_len, intron_end)`. Line 606 extracts:

```python
exon_seq = genome_seq[intron_end:intron_end + rescue_len]
```

This fetches the sequence **downstream** of the intron (i.e., intronic sequence on the wrong side), instead of the upstream exonic sequence. The correct expression is:

```python
exon_seq = genome_seq[intron_end - rescue_len:intron_end]
```

All minus-strand 5' splice rescue comparisons are made against the wrong genomic sequence, causing incorrect rescue decisions.

**Fix:** Change line 606 to `genome_seq[intron_end - rescue_len:intron_end]`.

---

#### NEW-004 — `mapPacBio` timeout not caught; hangs indefinitely on stuck aligner
**File:** `rectify/core/multi_aligner.py:385`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`sort_proc.communicate(timeout=ALIGNER_TIMEOUT)` for mapPacBio is called without a surrounding `try-except subprocess.TimeoutExpired` block. The minimap2 path at line 181 correctly wraps this call with:

```python
try:
    stdout, stderr = proc.communicate(timeout=ALIGNER_TIMEOUT)
except subprocess.TimeoutExpired:
    proc.kill()
    raise
```

If `mapPacBio` hangs (e.g., on a malformed BAM or when alignment stalls), the entire rectify process hangs indefinitely with no kill or error. This is particularly severe in SLURM jobs where wall-time overruns cause downstream array tasks to be silently killed.

**Fix:** Wrap `sort_proc.communicate(timeout=ALIGNER_TIMEOUT)` in the same try-except pattern used by minimap2.

---

#### NEW-005 — `mapPacBio` `view_proc` created without `stderr=subprocess.PIPE`; stdout buffer deadlock possible
**File:** `rectify/core/multi_aligner.py` (mapPacBio subprocess setup)
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The `view_proc` subprocess for the `samtools view` step of the mapPacBio pipeline is spawned without `stderr=subprocess.PIPE`. If `samtools view` writes a large error message to stderr (e.g., on a corrupt BAM), the stderr pipe buffer fills. Because the parent process is blocked waiting on `sort_proc.communicate()`, neither pipe is drained, producing a classic deadlock. uLTRA and deSALT pipelines in the same file both set `stderr=subprocess.PIPE`.

**Fix:** Add `stderr=subprocess.PIPE` to the `view_proc` `Popen()` call.

---

#### NEW-006 — `build_query_ref_map` / `cigar_to_events` always anchor to `reference_start` regardless of strand
**File:** `rectify/core/chimeric_consensus.py:104-178`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

Both `build_query_ref_map()` (line 104) and `cigar_to_events()` (line 145) initialize the reference position walk as:

```python
rpos = read.reference_start
```

For minus-strand reads, CIGAR operations traverse from `reference_start` to `reference_end` in genomic space, but the query sequence is reverse-complemented. Failing to adjust the anchor for strand means query→reference position mappings are mirror-reversed for minus-strand reads. Chimeric read detection and split-read coordinate merging are consequently inverted for all minus-strand alignments.

**Fix:** For minus-strand reads, initialize `rpos = read.reference_end - 1` and walk CIGAR operations in reverse order, or account for the query reversal in the position map construction.

---

#### NEW-007 — `_SCRATCH_TEARDOWN_BLOCK` runs rsync unconditionally even when `SCRATCH_DIR` is empty
**File:** `rectify/core/batch_command.py:162`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The generated SLURM teardown block runs:

```bash
rsync -a "$SCRATCH_DIR/" "$OUTPUT_DIR/"
```

unconditionally at the end of each array task. If `SCRATCH_DIR` is empty (i.e., scratch staging was not actually used for a given task), rsync exits with a non-zero code or syncs an empty directory, potentially clobbering existing Oak output with an empty tree. The `set -euo pipefail` header in generated scripts means this can abort jobs that would otherwise succeed.

**Fix:** Guard the teardown block: `if [ -n "$SCRATCH_DIR" ] && [ -d "$SCRATCH_DIR" ]; then rsync ...; fi`.

---

#### NEW-008 — Bash array injection via unescaped sample IDs, BAM paths, and annotation arguments
**File:** `rectify/core/batch_command.py:506, 509, 524, 530, 542`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

Sample IDs and BAM paths are embedded in bash here-doc arrays with only double-quote wrapping:

```python
f'  "{sample_id}"'    # line 506
f'  "{bam_path}"'     # line 509
```

Organism and annotation arguments use f-strings directly:

```python
f'--organism "{args.organism}"'    # line 524
f'--annotation "{args.annotation}"'  # line 530
```

A sample ID containing `$(...)` or a backtick, or an organism name containing special characters, produces shell injection in the generated SLURM script. Paths or IDs sourced from manifest TSVs (user-controlled) are particularly vulnerable.

**Fix:** Use `shlex.quote()` for all values embedded in generated shell scripts.

---

#### NEW-011 — Donor and acceptor sequences swapped for minus-strand junction validation
**File:** `rectify/core/analyze/junction_validation.py:151-156`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The minus-strand validation block:

```python
donor_seq = _reverse_complement(acceptor_seq)
acceptor_seq = _reverse_complement(genome.fetch(chrom, donor, donor+2))
```

fetches the donor dinucleotide into `acceptor_seq` and derives `donor_seq` from the original `acceptor_seq`. The labels are inverted: what the code calls `donor_seq` is actually derived from the acceptor window, and vice versa. Legitimate minus-strand GT-AG junctions fail validation, while some invalid junctions pass, silently corrupting junction reliability scores in downstream APA and shift analyses.

**Fix:** Swap the assignment order and/or variable names so `donor_seq = _reverse_complement(genome.fetch(chrom, donor, donor+2))` and `acceptor_seq = _reverse_complement(genome.fetch(chrom, acceptor-2, acceptor))`.

---

#### NEW-012 — Junction aggregate never reverse-complements donor/acceptor dinucleotides for minus strand
**File:** `rectify/core/aggregate/junctions.py:156-164`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`five_ss_dinuc = genome_seq[intron_start:intron_start+2]` is computed with no strand check. For minus-strand junctions, the donor dinucleotide is on the reverse complement strand and must be fetched from the opposite end of the intron and reverse-complemented. Compare to `utils/splice_motif.py` which correctly handles strand. All minus-strand junction aggregates report `GT` where the actual motif is `CT` (complement), causing every minus-strand junction to be misclassified as non-canonical.

**Fix:** Add strand-aware dinucleotide extraction, matching the implementation in `utils/splice_motif.py`.

---

#### NEW-013 — NET-seq minus-strand 3' position double-subtracted by `n_trimmed`
**File:** `rectify/core/netseq_bam_processor.py:275, 416`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

Line 275 already subtracts `n_trimmed` from the minus-strand 3' position during initial extraction. Line 416 then subtracts `n_trimmed` a second time:

```python
three_prime_raw = three_prime_corrected - n_trimmed   # line 416, minus strand
```

The plus-strand path at line 418 correctly adds `n_trimmed` only once. The double-subtraction shifts all minus-strand NET-seq 3' positions by `2 × n_trimmed` bases upstream, corrupting Pol II occupancy profiles for minus-strand genes.

**Fix:** Remove the subtraction at line 416 for minus strand, or consolidate the two subtraction sites into a single canonical position.

---

### Confirmed New Findings — HIGH

#### NEW-018 — `None` values from failed position extraction propagate to TSV and arithmetic
**File:** `rectify/core/bam_processor.py` (position extraction paths)
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

When position extraction returns `None` for edge-case reads (e.g., zero-length alignments, reads with only soft-clip in CIGAR), downstream code writes `str(None)` to TSV columns and passes `None` into arithmetic operations that raise `TypeError`. The BAM processor does not universally guard against `None` position values after calling extraction helpers.

**Fix:** Add explicit `if pos is None: continue` guards after all position extraction calls, with a stat counter for skipped reads.

---

#### NEW-019 — `dist <= 0` vs `dist < 0` asymmetry excludes valid boundary junction on minus strand
**File:** `rectify/core/consensus.py:404`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The junction proximity check uses `dist <= 0` for minus strand but `dist < 0` for plus strand. A read whose corrected position falls exactly on a junction boundary (`dist == 0`) is excluded only for minus strand, silently discarding valid junction-adjacent reads from the consensus selection logic.

**Fix:** Use a consistent comparison: `dist < 0` for both strands (or `dist <= 0` for both if the semantics require it — but apply symmetrically).

---

#### NEW-020 — `candidate_junctions` in consensus has no strand; opposite-strand junctions corrupt 5' rescue
**File:** `rectify/core/consensus.py:370-379`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The candidate junction list is built from genomic coordinates without a strand field. On a locus where plus- and minus-strand genes overlap, junctions from the opposite-strand gene are included in the candidate set. The 5' rescue logic then snaps plus-strand reads to minus-strand junction boundaries (and vice versa), producing incorrect 5' end corrections for reads in overlapping gene regions.

**Fix:** Include strand in the junction key and filter `candidate_junctions` to matching strand before using in 5' rescue.

---

#### NEW-021 — `mapPacBio` `view_proc.returncode` never checked after `communicate()`
**File:** `rectify/core/multi_aligner.py` (mapPacBio SAM→BAM pipeline)
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

After `sort_proc.communicate()`, neither `view_proc.returncode` nor `sort_proc.returncode` is inspected for the mapPacBio pipeline. A failed `samtools view` (e.g., malformed SAM output from mapPacBio) produces an empty or corrupt BAM that is silently passed to consensus selection. uLTRA and deSALT pipelines in the same file both check return codes.

**Fix:** After `communicate()`, check both `view_proc.returncode != 0` and `sort_proc.returncode != 0`; raise `RuntimeError` with captured stderr on failure.

---

#### NEW-022 — `continue_on_error` never declared in `run-all` parser; always silently False
**File:** `rectify/core/run_command.py`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`getattr(args, 'continue_on_error', False)` is used to control error propagation in `_run_multi_sample()`, but `--continue-on-error` is not registered in the `run-all` subparser in `cli.py`. The attribute is always absent, so `getattr` always returns `False`. Users who pass `--continue-on-error` to `rectify run-all` receive no error message (argparse in partial-parse mode silently drops unknown flags), and all sample failures abort the entire run.

**Fix:** Declare `--continue-on-error` in the `run-all` subparser.

---

#### NEW-023 — `use_scratch` never in `run-all` parser; always evaluates to `True` via `getattr` default
**File:** `rectify/core/run_command.py`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`getattr(args, 'use_scratch', True)` is used with a default of `True`, but `--no-use-scratch` / `--use-scratch` are not declared in the `run-all` subparser. The attribute is always absent, so scratch staging is always attempted even when `$SCRATCH` is unavailable (e.g., non-Sherlock environments), causing `make_job_scratch_dir()` to fall back to `$TMPDIR` silently with no user control.

**Fix:** Declare `--use-scratch / --no-use-scratch` in the `run-all` subparser; default to `True` only on Sherlock.

---

#### NEW-024 — No duplicate `sample_id` detection in manifest; later samples silently overwrite earlier ones
**File:** `rectify/core/run_command.py` (manifest loading)
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

Manifest TSVs with duplicate `sample_id` rows are accepted without warning. In multi-sample mode, per-sample output directories are keyed by `sample_id`; a duplicate causes the second sample's outputs to overwrite the first. In two-pass streaming manifest mode, the duplicate sample is counted twice in the count matrix, inflating its apparent read depth.

**Fix:** After loading the manifest, check for duplicate `sample_id` values and raise an error (or at minimum a prominent warning) if any are found.

---

#### NEW-025 — `streaming: true` in SLURM profile YAML silently ignored for `batch` command
**File:** `rectify/core/batch_command.py` (profile application)
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`_apply_profile()` does not include `streaming` in its `profile_to_arg` mapping (related to, but distinct from, C3's `use_scratch`). When `sherlock_larsms.yaml` sets `streaming: true`, the profile value is loaded but never applied to the generated script's `--streaming` flag. The generated `rectify correct` commands run without `--streaming`, loading entire BAMs into RAM (~30–40 GB per BAM) instead of using the 4–5 GB streaming path.

**Fix:** Add `'streaming': 'streaming'` to `_apply_profile()`'s mapping (this is partially covered by the C3 fix but `batch_command`'s own invocation of profile application needs the same fix).

---

#### NEW-026 — No minimum-condition guard before DESeq2 in `analyze_command.py`
**File:** `rectify/core/analyze_command.py`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The shift analysis path in `analyze_command.py` correctly guards DESeq2-style comparisons with a minimum condition count check. The primary DESeq2 invocation does not perform this check. If the manifest contains only one condition (or only one sample per condition), `pydeseq2` raises an uninformative exception deep in the statistical fitting code rather than a clear user-facing error.

**Fix:** Add a pre-flight check: if fewer than 2 conditions or fewer than 2 samples per condition, log an error and skip DESeq2 with a descriptive message.

---

#### NEW-028 — `--polya-model` path captured but never forwarded to BAM processor
**File:** `rectify/core/correct_command.py:206`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`polya_model_path = args.polya_model` is assigned at line 206 but is then absent from both `process_bam_streaming()` and `process_bam_file_parallel()` call sites. The trained model is never used; the hardcoded 80% A-richness threshold from `config.py` always applies. This is a superset of M3 (from Round 1), which identified the same issue at a higher level; Round 2 pinpoints the exact missing argument in `correct_command.py`.

**Fix:** Pass `polya_model_path` through to both processor call sites and consume it inside the BAM processor to override the default threshold.

---

#### NEW-029 — Empty output TSV from `extract_command.py` is indistinguishable from success
**File:** `rectify/core/extract_command.py:454-455`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

When all reads are filtered out (e.g., wrong strand, empty BAM, overly strict quality thresholds), `extract_command.py` writes a header-only TSV and exits with code 0. Downstream pipeline steps that consume this TSV either crash with an unhelpful empty-dataframe error or silently produce empty results. No warning is emitted to the user indicating that zero reads passed filtering.

**Fix:** After filtering, if the output DataFrame is empty, log a WARNING with the filter breakdown and optionally exit with code 1 (or at minimum warn loudly before writing the empty file).

---

#### NEW-031 — Junction index built without strand; opposite-strand junctions merged
**File:** `rectify/core/analyze/junction_analysis.py` (`build_junction_index_from_records`)
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`build_junction_index_from_records()` aggregates junctions by `(chrom, junction_start, junction_end)` tuple with no strand component. On loci with overlapping plus- and minus-strand genes, same-coordinate junctions from opposite strands are merged into a single entry, conflating their read counts. Junction usage statistics (PSI, reliability scores) are incorrect for all overlapping gene pairs.

**Fix:** Add strand to the junction aggregation key: `(chrom, strand, junction_start, junction_end)`.

---

#### NEW-032 — GFF junction strand parsed but never stored in junction dict
**File:** `rectify/core/analyze/junction_analysis.py:473`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The GFF loading loop reads the strand field at line 473 but does not include `'strand'` in the junction dict that is created and appended to the junction list. All downstream junction operations that attempt to access `junction['strand']` raise `KeyError`, and any logic that defaults to `'+'` silently produces wrong results for minus-strand junctions.

**Fix:** Add `'strand': strand` to the junction dict construction at the GFF loading site.

---

#### NEW-034 — Exception in motif upstream-region extraction silently skips entire sequence
**File:** `rectify/core/analyze/motif_discovery.py:122`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The upstream sequence extraction loop contains:

```python
except Exception:
    continue
```

This swallows any exception (including `IndexError` from sequences shorter than the extraction window, or `KeyError` from missing chromosomes) and skips the entire downstream k-mer extraction for that sequence. Because no warning is logged, enrichment analysis silently runs on a smaller sequence set than expected, potentially with systematic bias toward longer chromosomes where extraction rarely fails.

**Fix:** Replace the bare `except` with specific exception types; log a WARNING per skipped sequence with the failure reason.

---

#### NEW-035 — `NCBI_TO_CHR` reversal overwrites `chrI` key with `I`, breaking chromosome normalization
**File:** `rectify/utils/chromosome.py:64`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`YEAST_CHR_TO_NCBI` maps both `'chrI'` and `'I'` to the same NCBI accession. The reversal at line 64:

```python
NCBI_TO_CHR = {v: k for k, v in YEAST_CHR_TO_NCBI.items()}
```

iterates in dict-insertion order; `'I'` is inserted after `'chrI'`, so it overwrites the mapping. `normalize_chromosome(ncbi_acc)` returns `'I'` instead of the expected `'chrI'` for all yeast chromosomes, causing chromosome name mismatches when comparing against GFF-sourced features that use the `'chrI'` convention.

**Fix:** Explicitly prefer `'chrI'`-form keys in the reversal: only insert `k` if it starts with `'chr'`, or build `NCBI_TO_CHR` manually.

---

#### NEW-036 — `sync_to_oak()` silently continues after partial rsync failure
**File:** `rectify/slurm.py` (`sync_to_oak`)
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`sync_to_oak()` catches `subprocess.CalledProcessError` from rsync, logs a warning, then falls back to `shutil.copytree()` without re-raising. If rsync partially transfers files before failing (e.g., disk quota exceeded on Oak), the fallback `shutil.copytree()` may succeed on the remaining files, leaving Oak in an inconsistent state with no indication that data is missing. The calling code receives no exception and treats the sync as successful.

**Fix:** If rsync fails, log the error and re-raise, allowing the caller to decide whether to abort or retry. Do not silently continue with `shutil`.

---

#### NEW-037 — `--aligner` choices inconsistent between `correct` and `run-all` subparsers
**File:** `rectify/cli.py`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

The `correct` subparser declares `choices=['minimap2', 'bwa', 'star', 'auto']`; the `run-all` subparser declares `choices=['minimap2', 'star', 'bowtie2', 'bwa']`. The `'auto'` option is only available in `correct` (not `run-all`), and `'bowtie2'` is only available in `run-all` (not `correct`). A user passing `--aligner auto` to `run-all` gets an argparse error; a user passing `--aligner bowtie2` to `correct` gets an argparse error. The `run-all` command calls `correct` internally, so the mismatched choices propagate to confusing failures.

**Fix:** Unify the `choices` list across both subparsers; use a shared constant.

---

#### NEW-038 — `n_reads_total` inflated by counting reads before `min_reads` filter
**File:** `rectify/core/train_polya_command.py:290, 292`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** HIGH_IMPACT

`n_reads_total += len(clips)` at line 290 increments before the `if len(clips) < min_reads: continue` filter at line 292. Training summary statistics therefore report the total number of reads seen before filtering, not the number actually used for model training. A model trained on 1,000 reads after a strict filter may report `n_reads_total = 50,000`, misleading users about actual training data volume.

**Fix:** Move the increment to after the min_reads guard, or use a separate counter for pre-filter reads.

---

### Confirmed New Findings — MEDIUM / LOW

| ID | File | Severity | Title | V2 | V3 |
|----|------|----------|-------|----|----|
| NEW-039 | `core/bam_processor.py` | MEDIUM | Dead code: `_df = pd.read_csv(...)` assigned but unused; file read twice | NOT_TESTED | HIGH_IMPACT |
| NEW-040 | `core/consensus.py:544` | MEDIUM | `max(five_clip, five_clip + terminal_end)` is always `five_clip + terminal_end`; adds index to clip length | NOT_TESTED | HIGH_IMPACT |
| NEW-041 | `core/consensus.py` | MEDIUM | Broad `except Exception` in `_get_effective_5prime_clip` silently returns default on malformed reads | NOT_TESTED | HIGH_IMPACT |
| NEW-042 | `core/polya_trimmer.py:532` | MEDIUM | `clip_len` can be 0; `a_frac = clip_seq.count('A') / clip_len` raises `ZeroDivisionError` | NOT_TESTED | EDGE_CASE |
| NEW-043 | `core/polya_trimmer.py:545` | MEDIUM | Second `a_frac` calculation has same missing guard for `clip_len == 0` | NOT_TESTED | EDGE_CASE |
| NEW-044 | `core/polya_model.py:130` | MEDIUM | `int(k)` on JSON position-profile key with no `try-except`; `ValueError` on corrupted model file | NOT_TESTED | HIGH_IMPACT |
| NEW-045 | `core/run_command.py` | MEDIUM | `except Exception: pass` swallows `FileNotFoundError` in BAM validation | NOT_TESTED | EDGE_CASE |
| NEW-046 | `core/run_command.py` | MEDIUM | `output_dir / sample['sample_id']` path construction without sanitization; path traversal via `../` in manifest | NOT_TESTED | HIGH_IMPACT |
| NEW-047 | `core/run_command.py` | MEDIUM | No file existence check for manifest BAM/FASTQ paths before spawning workers | NOT_TESTED | HIGH_IMPACT |
| NEW-048 | `core/batch_command.py` | MEDIUM | `streaming: true` in profile silently ignored by `generate_slurm_scripts()`; 30–40 GB RAM per BAM used | NOT_TESTED | HIGH_IMPACT |
| NEW-049 | `core/analyze_command.py` | MEDIUM | Sample with zero positions skipped in pass 1 but included as all-zeros column in pass 2 count matrix | NOT_TESTED | EDGE_CASE |
| NEW-050 | `core/analyze_command.py` | MEDIUM | Bedgraph written with direct `open()`/`write()` without temp-file + `os.replace` pattern; partial write on crash | NOT_TESTED | HIGH_IMPACT |
| NEW-051 | `core/analyze/apa_detection.py` | LOW | `isoform_counter` global across all genes; non-deterministic IDs if function called multiple times | NOT_TESTED | UNREACHABLE |
| NEW-052 | `core/analyze/apa_detection.py` | MEDIUM | Comment says "modal isoform" but code calls `np.median()` | NOT_TESTED | HIGH_IMPACT |
| NEW-053 | `core/analyze/atract_refiner.py:272` | MEDIUM | Tie-break on signal height selects upstream peak in coordinate order; should prefer downstream (closer to CPA) | NOT_TESTED | HIGH_IMPACT |
| NEW-054 | `core/analyze/gene_attribution.py` | LOW | `primary_fraction` inflated by double-counting spanning reads (display metadata only) | NOT_TESTED | UNREACHABLE |
| NEW-055 | `core/analyze/gene_attribution.py` | LOW | `round()` percentages in attribution summary do not sum to 100 for edge-case distributions | NOT_TESTED | EDGE_CASE |
| NEW-056 | `core/analyze/junction_analysis.py` | MEDIUM | Greedy junction clustering compares to running mean; enables transitivity merging of distant junctions | NOT_TESTED | HIGH_IMPACT |
| NEW-057 | `core/provenance.py` | MEDIUM | Sidecar `.prov.json` written with direct `open()` + `json.dump()`; crash mid-write produces corrupt file | NOT_TESTED | HIGH_IMPACT |
| NEW-058 | `rectify/utils/provenance.py` | LOW | `datetime.now()` (naive) in `utils/provenance.py` vs `datetime.now(timezone.utc)` in `provenance.py`; inconsistent timestamps | NOT_TESTED | EDGE_CASE |
| NEW-059 | `rectify/utils/provenance.py` | MEDIUM | `Path.is_relative_to()` requires Python ≥ 3.9; pyproject.toml declares `python_requires >= 3.8`; `AttributeError` on 3.8 | CAUGHT_BY_TESTS | HIGH_IMPACT |
| NEW-060 | `core/classify/full_length_classifier.py:150-153` | MEDIUM | Both if/else branches assign `five_prime_pos = record.five_prime_corrected`; minus-strand should use `record.three_prime_corrected` | NOT_TESTED | HIGH_IMPACT |
| NEW-061 | `core/validate_command.py:520-531` | MEDIUM | Ground-truth search anchored only to `corrected_position`; `improvement_bp` corrupted when correction overshoots | NOT_TESTED | HIGH_IMPACT |
| NEW-062 | `core/analyze/gene_attribution.py` | LOW | Non-deterministic gene ordering when `overlap_bp` tied; output ordering varies across Python versions | NOT_TESTED | EDGE_CASE |
| NEW-063 | `core/analyze/gene_attribution.py` | LOW | `row.get('strand', '+')` silent default; minus-strand reads silently attributed to plus-strand features on unstranded data | NOT_TESTED | EDGE_CASE |

---

### Disputed Findings

#### NEW-009 — `tes_tolerance == 0` causes `ZeroDivisionError` in APA detection
**File:** `rectify/core/analyze/apa_detection.py:293`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** UNREACHABLE (disputed)

`round(tes_modal / tes_tolerance) * tes_tolerance` raises `ZeroDivisionError` if `tes_tolerance == 0`. V1 confirms the code has no guard. V3 argues this is unreachable in practice: `tes_tolerance` is set from a CLI argument that defaults to a positive integer and argparse enforces `type=int` with no explicit `min=0` check, meaning zero is technically passable by a user. V3 verdict: UNREACHABLE only if users don't pass `--tes-tolerance 0`; a defensive guard is nonetheless warranted.

**Verifier verdicts:** V1=CONFIRMED | V2=NOT_TESTED | V3=UNREACHABLE
**Recommendation:** Add `if tes_tolerance == 0: raise ValueError("--tes-tolerance must be > 0")` at argument parsing time.

---

#### NEW-010 — Isoform grouping key missing strand; cross-strand isoforms merged
**File:** `rectify/core/analyze/apa_detection.py:297`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** UNREACHABLE (disputed)

The isoform grouping dict key omits strand. V1 confirms the omission. V3 argues this is unreachable because gene IDs in rectify's annotation system are already strand-unique (a gene on + strand has a different ID than any gene on − strand), so two isoforms from different strands would never share the same gene ID. However, if a user provides a custom GFF where the same gene name appears on both strands, incorrect merging would occur silently.

**Verifier verdicts:** V1=CONFIRMED | V2=NOT_TESTED | V3=UNREACHABLE
**Recommendation:** Add strand to the grouping key as a defensive measure; cost is trivial.

---

#### NEW-033 — BH FDR correction applied over all GO terms including zero-overlap terms (over-conservative)
**File:** `rectify/core/analyze/go_enrichment.py`
**V1:** CONFIRMED | **V2:** NOT_TESTED | **V3:** DISPUTED

`n_tests = len(go_terms)` includes GO terms with zero gene overlap, increasing the multiple-testing correction denominator. V1 confirms this makes the correction over-conservative (harder to reach significance). V3 disputes whether this is a bug: the Benjamini-Hochberg procedure is commonly applied over all tests including those with p=1 (from zero overlap), and filtering zero-overlap terms before correction would inflate false discovery rates under some interpretations of the method. This is a statistical philosophy disagreement.

**Verifier verdicts:** V1=CONFIRMED | V2=NOT_TESTED | V3=DISPUTED
**Recommendation:** Document the behavior explicitly; optionally offer `--go-filter-zero-overlap` flag to let users choose.

---

### False Positives

| ID | File | Agent Claim | Why Refuted |
|----|------|-------------|-------------|
| NEW-003 | `core/splice_aware_5prime.py:625` | Minus-strand 3'SS coordinate wrong | Line 625 is geometrically correct for 0-based half-open minus-strand coordinates; verified against convention in CLAUDE.md |
| NEW-014 | `core/netseq_refiner.py` | `aggregate_proportional_results` double-counts results | Function is dead code — never called in the active pipeline; `refine_batch` call pattern does not match the claimed code path |
| NEW-015 | `core/validate_command.py` | Header/data column count mismatch (`is_correct_5bp` missing) | Both header and data rows have 16 columns; `is_correct_5bp` is absent from both consistently — no mismatch |
| NEW-016 | `core/terminal_exon_refiner.py:869-870` | Raw alignment coordinates used for splice site arithmetic | Lines 869-870 use pre-loaded annotated splice site positions from `splice_index`, not raw alignment coordinates; no off-by-one |
| NEW-017 | `core/terminal_exon_refiner.py:1149` | Minus-strand position arithmetic wrong | Line 1149 correctly uses `read_end - 1` consistent with 0-based half-open minus-strand 5' end convention |
| NEW-027 | `core/analyze/motif_discovery.py` | Background/foreground sets swapped for depletion analysis | Swap is intentional: depletion analysis compares observed (foreground) to expected (background); the inversion is semantically correct |
| NEW-030 | `core/analyze/shift_analysis.py` | `jensenshannon()` epsilon padding biases sparse estimates | `scipy.spatial.distance.jensenshannon()` handles zero-count vectors correctly by design; no epsilon padding is applied |

---

### Round 2 Priority Fix Order

| Priority | ID | Effort | Impact |
|----------|----|--------|--------|
| P1 | NEW-011 | Small | Donor/acceptor swap corrupts all minus-strand junction validation |
| P1 | NEW-013 | Small | Double-subtraction corrupts all minus-strand NET-seq 3' positions |
| P1 | NEW-006 | Medium | Chimeric consensus coordinate maps inverted for all minus-strand reads |
| P1 | NEW-002 | Small | Minus-strand 5' splice rescue fetches wrong genomic window |
| P1 | NEW-012 | Small | Junction aggregate never RC's minus-strand dinucleotides |
| P1 | NEW-060 | Small | Full-length classifier both branches identical; minus-strand positions wrong |
| P2 | NEW-035 | Small | `NCBI_TO_CHR` reversal breaks yeast chromosome normalization |
| P2 | NEW-001 | Small | Softclip/mismatch rescue functions ignore corrected position baseline |
| P2 | NEW-019 | Trivial | `dist <= 0` vs `dist < 0` excludes boundary junctions on minus strand |
| P2 | NEW-020 | Small | Candidate junctions missing strand; opposite-strand junctions corrupt 5' rescue |
| P2 | NEW-008 | Small | Shell injection via unescaped sample IDs and paths in generated SLURM scripts |
| P2 | NEW-052 | Trivial | `np.median()` used where comment/semantics says "modal" |
| P2 | NEW-031 | Small | Junction index key missing strand; opposite-strand junctions merged |
| P2 | NEW-032 | Trivial | GFF junction strand parsed but never stored in junction dict |
| P3 | NEW-004 | Small | `mapPacBio` timeout not caught; hangs indefinitely |
| P3 | NEW-005 | Trivial | `mapPacBio` `view_proc` missing `stderr=PIPE`; deadlock risk |
| P3 | NEW-021 | Small | `mapPacBio` return code never checked; corrupt BAM silently passed downstream |
| P3 | NEW-007 | Small | Scratch teardown rsync runs unconditionally even when `SCRATCH_DIR` is empty |
| P3 | NEW-036 | Small | `sync_to_oak()` silently continues after rsync failure |
| P3 | NEW-025 | Small | `streaming: true` in SLURM profile silently ignored; 30–40 GB RAM per BAM |
| P3 | NEW-048 | Small | Same issue as NEW-025 in `generate_slurm_scripts()` path |
| P3 | NEW-056 | Medium | Junction clustering transitivity merges distant junctions |
| P3 | NEW-053 | Small | A-tract refiner tie-break selects wrong (upstream) peak |
| P3 | NEW-037 | Trivial | `--aligner` choices inconsistent between `correct` and `run-all` subparsers |
| P3 | NEW-059 | Trivial | `Path.is_relative_to()` Python 3.9+ only; breaks on Python 3.8 |
| P4 | NEW-018 | Small | `None` positions propagate to TSV writes and arithmetic |
| P4 | NEW-022 | Trivial | `--continue-on-error` undeclared in `run-all` parser |
| P4 | NEW-023 | Trivial | `--use-scratch` undeclared in `run-all` parser |
| P4 | NEW-024 | Small | No duplicate sample_id detection in manifest |
| P4 | NEW-026 | Small | No condition-count guard before DESeq2 |
| P4 | NEW-028 | Small | `--polya-model` path not forwarded to BAM processor (see also M3) |
| P4 | NEW-029 | Small | Empty extract output indistinguishable from success |
| P4 | NEW-034 | Small | Bare `except` in motif extraction silently skips sequences |
| P4 | NEW-038 | Trivial | Training read count inflated by counting before min_reads filter |
| P4 | NEW-040 | Trivial | `max()` tautology in consensus; always returns second operand |
| P4 | NEW-046 | Small | Path traversal via `../` in manifest sample_id |
| P4 | NEW-061 | Small | Ground truth search too narrow in validate_command; `improvement_bp` corrupted |
