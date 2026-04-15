# RECTIFY Bugs to Fix

## Last Updated: 2026-04-11 (Bugs 37, 38, 41, 55 fixed in v0.9.0)

---

## Open

---

### NEW-066 (MEDIUM) — uLTRA reuses stale cached database with empty genome, silently produces zero alignments

**File:** `rectify/core/multi_aligner.py`, `run_ultra()` (~line 848)

**Symptom:** uLTRA exits 0 but produces a SAM with only `@SQ` headers and no alignment records. Observed for chunks 0–6 of `wt_by4742_rep1_chunked_20260412` across three separate SLURM runs.

**Root cause:** uLTRA caches its genome index in `{output_bam.stem}_ultra_tmp/database.db`. When the first run fails mid-way (after the `ultra_tmp/` dir is created but before the genome is indexed), the `_tmp_dir` decompressed genome is cleaned up. On retry, uLTRA prints "Database found in directory — using this one" and skips genome indexing. `refs_sequences.fa` is 0 bytes → "Number of ref seqs in fasta: 0" → all chromosomes excluded → zero alignments. `run_ultra()` only checks `sam_path.exists() and sam_path.stat().st_size > 0`, but uLTRA writes a valid-looking SAM with just headers even when it aligned nothing.

**Fix:** After `ultra_out_dir.mkdir()`, check if `database.db` exists but `refs_sequences.fa` is empty or missing — this is the fingerprint of a failed prior indexing run. Remove the stale directory and recreate it so uLTRA re-indexes from scratch. A valid cache (non-empty `refs_sequences.fa`) should be preserved — the index is genome/GTF-derived and is intentionally shared across chunks. Fixed in `multi_aligner.py` `run_ultra()`.

---

### ~~NEW-061 (HIGH) — `clip_read_to_corrected_3prime` shrinks N ops instead of eliminating them — Fixed 2026-04-09 (v2.7.9)~~

When `polya_walkback` produces a corrected_3prime that falls *inside* a near-3' N op (intron artifact), `clip_read_to_corrected_3prime` partially trims the N rather than eliminating it. The walkback position is numerically correct (validated shifts of −45, −72 for cat4_plus_1/2) but a stub N op remains in the rectified CIGAR (1520→1501 bp, 100→40 bp). Any downstream junction caller will report this residual N as a spurious splice junction.

**Root cause:** The clipping loop walks backward consuming `n_ref_clip` reference bases; when it runs out of bases mid-N, it shrinks the N and stops. There is no special handling to detect "target is inside N → snap to junction_start−1."

**Fix (two options):**
- Option A: In `clip_read_to_corrected_3prime`, detect when walking lands inside an N op and snap the clip to the N's start (fully clip the N and all trailing ops). Update corrected_3prime in the TSV to `junction_start − 1`.
- Option B: In `bam_processor._run_correction`, when FJF detects a near-3' artifact junction, propagate FJF's `corrected_3prime` (junction_start−1) and do not override it with the walkback position when the walkback lands inside the artifact N span.

**Files:** `rectify/core/bam_processor.py` (`clip_read_to_corrected_3prime`), `rectify/core/false_junction_filter.py`
**Discovered:** 2026-04-09 validation audit (cat4_plus_1, cat4_plus_2)

---

### ~~NEW-062 (MEDIUM) — `five_prime_rescued` not written to TSV; `correction_applied` omits 5' rescue — Fixed 2026-04-09 (v2.7.9)~~

When `rescue_3ss_truncation` corrects a read's 5' end, `five_prime_rescued=True` is stored in the result dict (bam_processor.py:277) but never written to any TSV column. The `correction_applied` field for cat3 reads shows only `atract_ambiguity,polya_walkback` — the 5' junction rescue is invisible to downstream consumers. Users cannot identify rescued reads from the TSV without comparing `five_prime_position` against raw alignment coordinates.

**Fix:** Add `five_prime_rescued` (bool) to the TSV header and row in `write_results_tsv`, and append `five_prime_rescued` to `correction_applied` when True.

**Files:** `rectify/core/bam_processor.py:567–613` (`write_results_tsv`)
**Discovered:** 2026-04-09 validation audit (cat3 reads)

---

### ~~NEW-063 (MEDIUM) — Chimeric reads with 3'-hard-clipped alignment not exempt from poly-A walkback — Fixed 2026-04-09 (v2.7.9)~~

Reads with XK=1 (chimeric reconstruction) that have a hard-clip at the 3' end are subject to the full polya_walkback pipeline. For cat5_plus_3aligner (100H hard-clip), the unresolved 3' end causes an ambiguity window of 211 bp (positions 9753–9964), producing 5 TSV rows with walkback shifts of 6–185 bp. The correction pipeline has no short-circuit for chimeric reads whose 3' sequence is not present in the query.

**Fix:** In `bam_processor.correct_read_3prime()`, skip polya_walkback (or cap ambiguity window to 0) when the read has a hard-clip at its 3' end (indicates missing sequence; walkback is unanchored).

**Files:** `rectify/core/bam_processor.py` (`correct_read_3prime`), `rectify/core/atract_detector.py`
**Discovered:** 2026-04-09 validation audit (cat5_plus_3aligner)

---

### ~~NEW-064 (LOW) — `netseq_refinement` listed in `correction_applied` when ambiguity_range=1 — Fixed 2026-04-09 (v2.7.9)~~

When `ambiguity_range == 1`, there is only one candidate position; NET-seq refinement runs but is a no-op. Listing `netseq_refinement` in `correction_applied` is misleading — it implies the NET-seq signal was used to resolve ambiguity when no ambiguity existed. Seen in cat6_minus_single (ambiguity_range=1, fraction=1.0, correction_applied includes `netseq_refinement`).

**Fix:** In `bam_processor.correct_read_3prime()`, only add `netseq_refinement` to `correction_applied` when `ambiguity_range > 1` (i.e., when NET-seq signal was actually consulted to break a tie).

**Files:** `rectify/core/bam_processor.py` (correction_applied assembly)
**Discovered:** 2026-04-09 validation audit (cat6_minus_single)

---

### ~~Bug 37 (HIGH) — Zero unit tests for `terminal_exon_refiner.py` — Fixed 2026-04-11 (v0.9.0)~~

The module (1690 lines, multiple coordinate-sensitive code paths including Bugs 29, 33, 35) had no dedicated test file. Highest-priority test gap.

**Fix:** Added `tests/test_terminal_exon_refiner.py` with 51 tests covering: `SpliceSiteIndex` add/lookup/nearby, `load_splice_sites_from_gff` (plus and minus strand coordinate conversion, GFF position conflicts), `merge_splice_indices`, `detect_junction_truncated_reads` (both strands), `detect_partial_junction_crossings`, `get_soft_clip_info`, `simple_align`, `detect_mismatch_clusters`. Real-data class (`TestLoadSpliceSitesFromBundledGff`, `TestWithRealBam`) validates known S. cerevisiae intron positions (YAL030W, YAL001C) against the bundled R64-5-1 GFF and wt_by4742_rep1 BAM.
**File:** `tests/test_terminal_exon_refiner.py`

---

### ~~Bug 38 (HIGH) — `consensus.py` alignment selection only indirectly tested — Fixed 2026-04-11 (v0.9.0)~~

Core per-read consensus selection logic was exercised only through two peripheral tests (`test_xr_flag.py`, `test_gapmm2_seq_restore.py`). Tie-breaking, penalty scoring, and multi-aligner merging had no direct unit tests.

**Fix:** Added `tests/test_consensus_selection.py` with 40 tests covering: `extract_junctions_from_cigar` (multi-intron, soft-clip vs N op), `check_canonical_splice_sites` (GT/AG, GC/AG, non-canonical, real YAL030W intron), `score_alignment` (5' clip penalty −2/base, A-tract depth cap at 10, 3' clip penalty), `select_best_alignment` (winner selection, `was_5prime_rescued` flag, tiebreaker by annotated junction count / 3' agreement, `confidence` levels). Real-data class validates against wt_by4742_rep1 BAM at the YAL030W locus.
**File:** `tests/test_consensus_selection.py`

---

### ~~Bug 41 (MEDIUM) — Trained `--polya-model` is never used — Fixed 2026-04-11 (v0.9.0)~~

`--polya-model` flag, model training pipeline, and `load_model()` infrastructure existed, but `polya_model_path` was captured in config and then discarded. The hardcoded 80% A-richness threshold was always applied regardless.

**Fix:** Wired `polya_model_path` through `bam_processor.correct_read_3prime()` (both streaming and parallel modes). Added `pt_tag`, `polya_score`, and `polya_source` columns to `corrected_3ends.tsv`. Added `rectify tag-polya` subcommand for retroactive annotation of aligned BAMs. Added unaligned dorado BAM auto-detection and tag-preserving alignment pipeline in `preprocess.py`.
**Files:** `rectify/core/bam_processor.py`, `rectify/core/correct_command.py`, `rectify/core/tag_polya_command.py`, `rectify/core/preprocess.py`, `rectify/cli.py`

---

### ~~Bug 55 (MEDIUM) — Several APA clustering parameters not CLI-configurable — Fixed 2026-04-11 (v0.9.0)~~

`DEFAULT_MIN_PEAK_SEPARATION = 5`, `DEFAULT_MAX_CLUSTER_RADIUS = 10`, and `DEFAULT_MIN_SAMPLES = 2` had no corresponding CLI arguments.

**Fix:** Added `--min-peak-sep`, `--max-cluster-radius`, `--min-cluster-samples` to `create_analyze_parser()`. Both clustering call sites (single-sample and manifest mode) dispatch to `cluster_cpa_sites_adaptive()` when any non-default value is provided; otherwise fall through to the existing fixed-distance `cluster_cpa_sites()`.
**Files:** `rectify/core/analyze_command.py`

---

## Fixed (v2.7.8 — 2026-04-09)

---
### CRITICAL (Round 2 — Fixed v2.7.8)
---

### ~~NEW-011 (CRITICAL) — Donor and acceptor sequences swapped for minus-strand junction validation — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze/junction_validation.py:153-155`: the minus-strand block assigns `donor_seq = _reverse_complement(acceptor_seq)` and then `acceptor_seq = _reverse_complement(genome.fetch(chrom, donor, donor+2))`. The labels are inverted: what the code calls `donor_seq` is derived from the acceptor genomic window, and vice versa. Every minus-strand GT-AG junction fails validation; some invalid junctions pass.

**Fix:** Swap the assignments: `donor_seq = _reverse_complement(genome.fetch(chrom, donor, donor+2))` and `acceptor_seq = _reverse_complement(genome.fetch(chrom, acceptor-2, acceptor))`.
**File:** `rectify/core/analyze/junction_validation.py:151-156`

---

### ~~NEW-012 (CRITICAL) — Junction aggregate never reverse-complements donor/acceptor dinucleotides for minus strand — Fixed 2026-04-09 (v2.7.8)~~

`core/aggregate/junctions.py:158-159`: `five_ss_dinuc = genome_seq[intron_start:intron_start+2]` and `three_ss_dinuc = genome_seq[intron_end-2:intron_end]` with no strand check. For minus-strand junctions the donor is on the opposite strand and must be fetched from the reverse-complement end of the intron. Every minus-strand junction is misclassified as non-canonical.

**Fix:** Add strand-aware dinucleotide extraction matching `utils/splice_motif.py` which handles this correctly.
**File:** `rectify/core/aggregate/junctions.py:156-164`

---

### ~~NEW-013 (CRITICAL) — NET-seq minus-strand 3' position double-subtracted by `n_trimmed` — Fixed 2026-04-09 (v2.7.8)~~

`core/netseq_bam_processor.py:275` subtracts `n_trimmed` from the minus-strand 3' position during initial extraction. Line 416 subtracts it a second time: `three_prime_raw = three_prime_corrected - n_trimmed`. The plus-strand path correctly adds `n_trimmed` once to undo the subtraction. All minus-strand NET-seq 3' positions are shifted upstream by `2 × n_trimmed`.

**Fix:** Remove the duplicate subtraction at line 416 for minus strand, or consolidate both sites into one.
**File:** `rectify/core/netseq_bam_processor.py:275, 416`

---

### ~~NEW-001 (CRITICAL) — `rescue_softclip_at_homopolymer` / `rescue_mismatch_inside_homopolymer` ignore corrected position — Fixed 2026-04-09 (v2.7.8)~~

Both functions initialize `raw_pos = read.reference_end - 1`, discarding any correction already applied upstream. In contrast, `rescue_polya_prefix_in_softclip()` correctly accepts a `current_pos` parameter. When indel correction pipelines call these two functions after an initial correction, the re-anchor silently undoes the prior correction and applies its own walk from the original BAM position.

**Fix:** Add `current_pos` parameter to both functions; replace `raw_pos = read.reference_end - 1` with `raw_pos = current_pos`.
**File:** `rectify/core/indel_corrector.py:413-419, 654-659` (function definitions), called at lines ~1402, ~1415

---

### ~~NEW-002 (CRITICAL) — Minus-strand exon sequence fetched from wrong genomic window in 5' splice rescue — Fixed 2026-04-09 (v2.7.8)~~

`core/splice_aware_5prime.py:609`: for minus-strand reads, `exon_seq = genome_seq[intron_end - rescue_len:intron_end]` fetched sequence from the **end of the intron** instead of exon1. For minus-strand genes, exon1 (upstream in transcript order) is at **higher** genomic coordinates than the intron — `intron_end` (exclusive) marks where exon1 begins in genomic space. The soft-clip sequence stored in the BAM is in forward-strand orientation (RC of RNA), so no reverse-complement is needed; the correct window is `genome_seq[intron_end:intron_end + rescue_len]`. All minus-strand cat3 5' splice rescue calls compared the read tail against intronic sequence, causing the position to remain unchanged.

**Fix:** Changed line 609 from `genome_seq[intron_end - rescue_len:intron_end]` to `genome_seq[intron_end:intron_end + rescue_len]`.
**File:** `rectify/core/splice_aware_5prime.py:609`

---

### ~~NEW-008 (CRITICAL) — Shell injection via unescaped sample IDs and paths in generated SLURM scripts — Fixed 2026-04-09 (v2.7.8)~~

`core/batch_command.py:520-547`: sample IDs, BAM paths, organism, netseq-dir, and annotation arguments are embedded in generated bash scripts with only double-quote wrapping (`f'  "{sample_id}"'`), not `shlex.quote()`. A sample ID containing `$(...)` or backtick from a manifest TSV produces shell injection in the SLURM script.

**Fix:** Use `shlex.quote()` for all values embedded in generated shell scripts.
**File:** `rectify/core/batch_command.py:520-547`

---

### ~~NEW-035 (CRITICAL) — `NCBI_TO_CHR` reversal overwrites `chrI` form with bare `I` — Fixed 2026-04-09 (v2.7.8)~~

`utils/chromosome.py:64`: `YEAST_CHR_TO_NCBI` maps both `'chrI'` and `'I'` to the same NCBI accession. The dict-comprehension reversal `{v: k for k, v in ...}` iterates in insertion order; `'I'` is inserted after `'chrI'` and overwrites it. `normalize_chromosome(ncbi_acc)` returns `'I'` instead of `'chrI'`, causing chromosome name mismatches against GFF-sourced features.

**Fix:** Prefer `'chrI'`-form keys in the reversal: filter to only `k.startswith('chr')`, or build `NCBI_TO_CHR` manually.
**File:** `rectify/utils/chromosome.py:64`

---

### ~~NEW-007 (CRITICAL) — Scratch teardown rsync runs unconditionally even when `SCRATCH_DIR` is empty — Fixed 2026-04-09 (v2.7.8)~~

`core/batch_command.py:162`: `_SCRATCH_TEARDOWN_BLOCK` runs `rsync -a "${{SCRATCH_DIR}}/" "$OAK_OUTPUT_DIR/"` with no guard. If `SCRATCH_DIR` is unset or empty, rsync exits non-zero or clobbers Oak output with an empty tree. With `set -euo pipefail` in generated scripts, this aborts jobs that would otherwise succeed.

**Fix:** Guard the teardown: `if [ -n "$SCRATCH_DIR" ] && [ -d "$SCRATCH_DIR" ]; then rsync ...; fi`.
**File:** `rectify/core/batch_command.py:162`

---

### ~~NEW-004 (CRITICAL) — `mapPacBio` `sort_proc.communicate()` timeout not caught; hangs indefinitely — Fixed 2026-04-09 (v2.7.8)~~

`core/multi_aligner.py:385`: `sort_proc.communicate(timeout=ALIGNER_TIMEOUT)` for mapPacBio lacks a surrounding `try-except subprocess.TimeoutExpired`. The minimap2 path at line ~181 wraps this in try-except with `proc.kill()`. If mapPacBio hangs, the rectify process hangs indefinitely — severe in SLURM where wall-time overruns silently kill downstream tasks.

**Fix:** Wrap `sort_proc.communicate(timeout=ALIGNER_TIMEOUT)` in the same try-except pattern as minimap2.
**File:** `rectify/core/multi_aligner.py:385`

---

### ~~NEW-005 (CRITICAL) — `mapPacBio` `view_proc` deadlock — Fixed 2026-04-09 (v2.7.8)~~

`core/multi_aligner.py` (mapPacBio SAM→BAM setup): `view_proc` was spawned with `stderr=subprocess.PIPE` but no thread drained it. When samtools writes enough to fill the OS pipe buffer, `view_proc` blocks; since the parent is blocked on `sort_proc.communicate()`, neither pipe drains → classic deadlock producing a truncated BAM.

**Fix:** Changed `view_proc` to `stderr=subprocess.DEVNULL` (discards aligner stderr; error details come from `sort_proc`). Added `view_proc.kill()` in the timeout handler. `view_proc.wait()` + returncode check added after `sort_proc.communicate()`.
**File:** `rectify/core/multi_aligner.py`

---

### ~~NEW-003 (CRITICAL) — Bundled genome was regular gzip with stale bgzip index; pysam crash on fetch — Fixed 2026-04-09 (v2.7.8)~~

`data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz`: the file was replaced with a regular gzip archive on 2026-04-05, but the `.fai` and `.gzi` index files from the prior bgzip-compressed version (dated 2026-03-29) were left in place. pysam's `FastaFile()` opened without error (it uses the index to find chromosome positions) but `fasta.fetch(chrom)` raised `ValueError: failure when retrieving block` because random access requires bgzip BGZF blocks. This caused `rectify align` consensus selection to crash for every run that used the bundled genome.

**Fix:** Re-compressed with `bgzip -@ 4` and re-indexed with `samtools faidx`. Verified with `pysam.FastaFile().fetch('chrI', 0, 20)`.
**File:** `rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz`

---

### ~~NEW-006 (FEATURE) — `--write-corrected-bam`: hard-clip reads at corrected 3' position via CIGAR surgery — Added 2026-04-09 (v2.7.8)~~

`core/bam_processor.py` + `core/correct_command.py` + `cli.py`: the existing `--output-bam` flag only strips soft-clipped poly-A tails; reads with **aligned** A-mismatches (Cat1 poly-A walkback) still show the full alignment in IGV with A-mismatch tails. The new `--write-corrected-bam` flag applies CIGAR surgery to hard-clip every read at its corrected 3' position. Works for both strands: plus strand clips from the right (removes trailing ops, appends hard-clip, shrinks or removes the last CIGAR operations), minus strand clips from the left (removes leading ops, prepends hard-clip, updates `reference_start`). Output is coordinate-sorted and indexed.

**Implementation:** `clip_read_to_corrected_3prime()` + `write_corrected_bam()` in `core/bam_processor.py`; wired in `core/correct_command.py` after the correction TSV is written; `--write-corrected-bam PATH` argument added to `cli.py` correct subparser.

---
### HIGH (Round 2 — Fixed v2.7.8)
---

### ~~NEW-031 (HIGH) — Junction index built without strand; opposite-strand junctions merged — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze/junction_analysis.py:496`: junction index key is `(chrom, donor, acceptor)` with no strand. On loci where plus- and minus-strand genes overlap, same-coordinate junctions from opposite strands are merged into a single entry, conflating read counts. Junction PSI and reliability scores are incorrect for all overlapping gene pairs.

**Fix:** Add strand to the aggregation key: `(chrom, strand, donor, acceptor)`.
**File:** `rectify/core/analyze/junction_analysis.py`

---

### ~~NEW-032 (HIGH) — GFF junction strand parsed but never stored in junction dict — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze/junction_analysis.py:473`: `strand = parts[6]` is read but not included in the exon tuple (line 501: `exons.append((start, end, chrom))`). All downstream junction operations that access `junction['strand']` raise `KeyError`; logic that defaults to `'+'` silently produces wrong results for minus-strand junctions.

**Fix:** Add `'strand': strand` to the junction data structure at the GFF loading site.
**File:** `rectify/core/analyze/junction_analysis.py:473, 501`

---

### ~~NEW-018 (HIGH) — `None` positions from failed extraction propagate to TSV and arithmetic — Fixed 2026-04-09 (v2.7.8)~~

`core/bam_processor.py`: `get_read_3prime_position()` returns `(None, strand)` for zero-length/unmapped reads. The `None` position is stored directly to the result dict (`'original_3prime': original_position`) and written to TSV as the string `"None"`. Downstream arithmetic raises `TypeError`.

**Fix:** Add `if original_position is None: continue` guard after each position extraction call; increment a stat counter for skipped reads.
**File:** `rectify/core/bam_processor.py:191, 253-254`

---

### ~~NEW-019 (HIGH) — `dist <= 0` vs `dist < 0` asymmetry excludes valid boundary junction on minus strand — Fixed 2026-04-09 (v2.7.8)~~

`core/consensus.py`: plus-strand proximity check uses `if dist < 0` (line 378); minus-strand uses `if dist <= 0` (line 404). A read whose 5' end falls exactly on a junction boundary (`dist == 0`) is rescued on plus strand but excluded on minus strand.

**Fix:** Use consistent comparison: `dist < 0` for both strands (or `<= 0` for both — apply symmetrically).
**File:** `rectify/core/consensus.py:378, 404`

---

### ~~NEW-020 (HIGH) — `candidate_junctions` has no strand; opposite-strand junctions corrupt 5' rescue — Fixed 2026-04-09 (v2.7.8)~~

`core/consensus.py:899-907`: the candidate junction set is 3-tuples `(chrom, intron_start, intron_end)` with no strand. The rescue loop iterates these without strand filtering, potentially snapping plus-strand reads to minus-strand junction boundaries in overlapping gene regions.

**Fix:** Include strand in the junction key and filter `candidate_junctions` to matching strand before rescue.
**File:** `rectify/core/consensus.py:899-907`

---

### ~~NEW-021 (HIGH) — `mapPacBio` `view_proc.returncode` still not checked — Fixed 2026-04-09 (v2.7.8)~~

`core/multi_aligner.py` (mapPacBio SAM→BAM): `sort_proc.returncode` is checked after `communicate()` but `view_proc.returncode` is not. A failed `samtools view` (malformed SAM) produces a corrupt BAM silently passed to consensus. uLTRA and deSALT both call `view_proc.wait()` and check both return codes.

**Fix:** After `communicate()`, call `view_proc.wait()` and raise `RuntimeError` if `view_proc.returncode != 0`.
**File:** `rectify/core/multi_aligner.py` (mapPacBio section, ~line 387)

---

### ~~NEW-022 (HIGH) — `--continue-on-error` undeclared in `run-all` parser — Fixed 2026-04-09 (v2.7.8)~~

`run_command.py` uses `getattr(args, 'continue_on_error', False)` to control error propagation in `_run_multi_sample()`, but `--continue-on-error` is not registered in the `run-all` subparser in `cli.py`. The attribute is always absent; all sample failures abort the entire run even if the user passes the flag.

**Fix:** Declare `--continue-on-error` in the `run-all` subparser in `cli.py`.
**File:** `rectify/core/run_command.py`, `rectify/cli.py`

---

### ~~NEW-023 (HIGH) — `--use-scratch` undeclared in `run-all` parser; scratch always attempted — Fixed 2026-04-09 (v2.7.8)~~

`run_command.py:831`: `getattr(args, 'use_scratch', True)` defaults to `True`, but `--use-scratch` is not in the `run-all` subparser. Users cannot opt out; scratch staging is always attempted even in environments without `$SCRATCH`.

**Fix:** Declare `--use-scratch / --no-use-scratch` in the `run-all` subparser; default to `True` only when scratch is detected.
**File:** `rectify/core/run_command.py`, `rectify/cli.py`

---

### ~~NEW-024 (HIGH) — No duplicate `sample_id` detection in manifest — Fixed 2026-04-09 (v2.7.8)~~

`core/batch_command.py` `parse_manifest()`: duplicate `sample_id` rows are accepted without warning. Second sample's outputs overwrite the first, and in two-pass streaming mode the duplicate sample is counted twice in the count matrix, inflating apparent read depth.

**Fix:** After loading the manifest, check for duplicate `sample_id` values and raise an error (or prominent warning).
**File:** `rectify/core/batch_command.py` (`parse_manifest()`)

---

### ~~NEW-026 (HIGH) — No minimum-condition guard before primary DESeq2 invocation — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze_command.py:351-389`: the primary DESeq2 block relies only on `if args.run_deseq2:` without checking that there are ≥ 2 conditions or ≥ 2 samples per condition. The shift-analysis path at line 477 correctly guards with `len(sample_metadata['condition'].unique()) >= 2`. With a single-condition manifest, `pydeseq2` raises an uninformative deep exception.

**Fix:** Add a pre-flight check before the DESeq2 block: if fewer than 2 conditions or 2 samples per condition, log an error and skip DESeq2.
**File:** `rectify/core/analyze_command.py:351`

---

### ~~NEW-028 (HIGH) — `--polya-model` path captured but never forwarded to BAM processor — Fixed 2026-04-09 (v2.7.8)~~

`core/correct_command.py:206`: `polya_model_path = getattr(args, 'polya_model', None)` is assigned but absent from both `process_bam_streaming()` and `process_bam_file_parallel()` call sites. Both processor functions lack a `polya_model_path` parameter. The trained model is never used; hardcoded 80% A-richness threshold always applies.

**Fix:** Pass `polya_model_path` through to both processor call sites and consume it inside the BAM processor to override the default threshold. (Related to Bug 41.)
**File:** `rectify/core/correct_command.py:206`, `rectify/core/bam_processor.py`

---

### ~~NEW-029 (HIGH) — Empty output TSV from `extract_command.py` indistinguishable from success — Fixed 2026-04-09 (v2.7.8)~~

`core/extract_command.py:454-457`: when all reads are filtered, an empty DataFrame is written as a header-only TSV and the process exits with code 0. No warning is emitted. Downstream pipeline steps either crash on the empty file or silently produce empty results.

**Fix:** After filtering, if output is empty, log a WARNING with the filter breakdown before writing. Optionally exit with code 1.
**File:** `rectify/core/extract_command.py:454-457`

---

### ~~NEW-034 (HIGH) — Bare `except Exception` in motif extraction silently skips sequences — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze/motif_discovery.py:122`: `except Exception: extraction_errors += 1; continue` — no logging. Any exception (IndexError from short sequences, KeyError from missing chromosomes) causes silent skip. Enrichment analysis runs on a smaller, potentially biased sequence set with no diagnostic output.

**Fix:** Replace with specific exception types; log a WARNING per skipped sequence with the failure reason.
**File:** `rectify/core/analyze/motif_discovery.py:120-122`

---

### ~~NEW-036 (HIGH) — `sync_to_oak()` silently continues after partial rsync failure — Fixed 2026-04-09 (v2.7.8)~~

`slurm.py` `sync_to_oak()`: catches `subprocess.CalledProcessError` from rsync and falls back to `shutil.copytree()` without re-raising. If rsync partially transferred files before failing (e.g., disk quota exceeded), Oak is left in an inconsistent state. The caller receives no exception and treats the sync as successful.

**Fix:** Log the rsync failure and re-raise, allowing the caller to decide whether to abort or retry. Do not silently fall back to shutil.
**File:** `rectify/slurm.py` (`sync_to_oak`)

---

### ~~NEW-037 (HIGH) — `--aligner` choices inconsistent between `correct` and `run-all` subparsers — Fixed 2026-04-09 (v2.7.8)~~

`cli.py:106` (correct): `choices=['minimap2', 'bwa', 'star', 'auto']`. `cli.py:607` (run-all): `choices=['minimap2', 'star', 'bowtie2', 'bwa']`. `'auto'` is absent from `run-all`; `'bowtie2'` is absent from `correct`. Since `run-all` calls `correct` internally, the mismatch produces confusing argparse errors.

**Fix:** Unify the `choices` list via a shared constant.
**File:** `rectify/cli.py:106, 607`

---

### ~~NEW-038 (HIGH) — `n_reads_total` inflated by counting reads before `min_reads` filter — Fixed 2026-04-09 (v2.7.8)~~

`core/train_polya_command.py:290-292`: `training_data.n_reads_total += len(clips)` increments before `if len(clips) < min_reads: continue`. Sites that fail the threshold still inflate the total count. Training reports may show `n_reads_total = 50,000` when only 1,000 reads were actually used.

**Fix:** Move the increment to after the min_reads guard, or use a separate counter for pre-filter reads.
**File:** `rectify/core/train_polya_command.py:290, 292`

---
### MEDIUM (Round 2 — Fixed v2.7.8)
---

### ~~NEW-009 (MEDIUM) — `tes_tolerance == 0` causes `ZeroDivisionError` in APA detection — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze/apa_detection.py:304`: `_bin = round(tes_modal / tes_tolerance) * tes_tolerance` raises `ZeroDivisionError` if `--tes-tolerance 0` is passed. No guard at parsing or at the call site.

**Fix:** Add `if tes_tolerance <= 0: raise ValueError("--tes-tolerance must be > 0")` at argument-parsing time.
**File:** `rectify/core/analyze/apa_detection.py:304`

---

### ~~NEW-010 (MEDIUM) — Isoform grouping key missing strand; cross-strand isoforms merged — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze/apa_detection.py:307-308`: `key = (gene_id, gene_name, junction_sig, tes_key)` omits strand. For genomes where gene IDs are not strand-unique, isoforms from opposite strands with the same gene ID / junction signature / TES position are merged into one group, and the stored strand value is overwritten on each read.

**Fix:** Add strand to the key: `key = (gene_id, gene_name, junction_sig, tes_key, strand)`.
**File:** `rectify/core/analyze/apa_detection.py:307`

---

### ~~NEW-039 (MEDIUM) — Dead code: `_df` assigned but unused; input file read twice — Fixed 2026-04-09 (v2.7.8)~~

`core/bam_processor.py:694-707`: `_df = pd.read_csv(str(results_or_path), ...)` is assigned at line 694 with an always-true ternary condition (`if True else ...`) and never used. The file is then read again into `_df2` at lines 699 and 704.

**Fix:** Remove the `_df` assignment and the dead `if True else` ternary.
**File:** `rectify/core/bam_processor.py:694-695`

---

### ~~NEW-040 (MEDIUM) — `max(five_clip, five_clip + terminal_end)` is a tautology — Fixed 2026-04-09 (v2.7.8)~~

`core/consensus.py:544`: at this point `terminal_end > 0` is guaranteed (early return at line 537 handles zero, early return at line 542 exits if `total_errors < min_errors`). `max(five_clip, five_clip + terminal_end)` always returns `five_clip + terminal_end`. The `max()` adds confusion with no defensive value.

**Fix:** Replace with `return five_clip + terminal_end` directly.
**File:** `rectify/core/consensus.py:544`

---

### ~~NEW-041 (MEDIUM) — Broad `except Exception` in `_get_effective_5prime_clip` silently returns default — Fixed 2026-04-09 (v2.7.8)~~

`core/consensus.py`: `_get_effective_5prime_clip` wraps `read.get_aligned_pairs()` in a bare `except Exception: return five_clip`. Any unexpected failure (memory error, corrupted pysam object) is silently swallowed and the default clip length returned.

**Fix:** Narrow to specific exceptions (e.g., `ValueError`, `RuntimeError`); log at WARNING level before returning default.
**File:** `rectify/core/consensus.py` (`_get_effective_5prime_clip`)

---

### ~~NEW-042 (MEDIUM) — `clip_len == 0` causes `ZeroDivisionError` in poly(A) trimmer (plus strand) — Fixed 2026-04-09 (v2.7.8)~~

`core/polya_trimmer.py:532`: `a_frac = clip_seq.count('A') / clip_len` with no guard against `clip_len == 0`. A CIGAR soft-clip with length 0 (technically valid) causes `ZeroDivisionError`.

**Fix:** Add `if clip_len == 0: return read, 0` before the division.
**File:** `rectify/core/polya_trimmer.py:530-532`

---

### ~~NEW-043 (MEDIUM) — Same `clip_len == 0` ZeroDivisionError on minus strand — Fixed 2026-04-09 (v2.7.8)~~

`core/polya_trimmer.py:545`: `t_frac = clip_seq.count('T') / clip_len` — same unguarded division in the minus-strand branch.

**Fix:** Add `if clip_len == 0: return read, 0` before the division.
**File:** `rectify/core/polya_trimmer.py:543-545`

---

### ~~NEW-044 (MEDIUM) — `int(k)` on JSON model key with no try-except; crashes on corrupted model file — Fixed 2026-04-09 (v2.7.8)~~

`core/polya_model.py:128-131`: `{int(k): v for k, v in data['position_profile'].items()}` — any non-integer key in a corrupted JSON model file raises `ValueError` with no user-friendly error message.

**Fix:** Wrap the dict comprehension in a try-except; raise a descriptive `ValueError` naming the bad key.
**File:** `rectify/core/polya_model.py:128-131`

---

### ~~NEW-045 (MEDIUM) — `except Exception: pass` swallows `FileNotFoundError` in BAM MD-tag check — Fixed 2026-04-09 (v2.7.8)~~

`core/run_command.py` `_bam_has_md_tags()`: bare `except Exception: pass` returns `False` (no MD tags) when the BAM file does not exist. Callers can't distinguish "file missing" from "file has no MD tags."

**Fix:** Narrow the exception handler; let `FileNotFoundError` propagate or raise a descriptive error.
**File:** `rectify/core/run_command.py` (`_bam_has_md_tags`)

---

### ~~NEW-046 (MEDIUM) — Path traversal via `../` in manifest `sample_id` — Fixed 2026-04-09 (v2.7.8)~~

`core/run_command.py:569`: `sample_output = output_dir / sample_id` — a manifest with `sample_id = "../../etc/passwd"` resolves outside `output_dir`. No sanitization is performed.

**Fix:** Validate `sample_id` contains no path separator before constructing the output path (e.g., `if '/' in sample_id or sample_id.startswith('.'): raise ValueError`).
**File:** `rectify/core/run_command.py:569`

---

### ~~NEW-047 (MEDIUM) — No file existence check for manifest paths before spawning workers — Fixed 2026-04-09 (v2.7.8)~~

`core/run_command.py` (manifest dispatch): the code validates the `path` column exists but does not check that files are present before submitting all tasks to the thread pool. Missing files are discovered one-by-one as workers start, producing scattered error messages instead of a clear pre-flight failure.

**Fix:** Before spawning workers, check `Path(s['path']).exists()` for each manifest entry; collect all missing files and raise a single descriptive error.
**File:** `rectify/core/run_command.py` (manifest loading/dispatch section)

---

### ~~NEW-049 (MEDIUM) — Sample with zero positions skipped in pass 1 but added as all-zeros in pass 2 — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze_command.py` (two-pass streaming): pass 1 skips samples with no valid positions (`if not _agg: continue`). Pass 2 iterates the original manifest and fills missing sample columns with zeros, adding the zero-position sample back as an all-zeros column — inflating the sample count for DESeq2 and shift analysis.

**Fix:** Track which samples were skipped in pass 1; exclude them from both pass 2 and downstream tools.
**File:** `rectify/core/analyze_command.py` (two-pass streaming logic)

---

### ~~NEW-050 (MEDIUM) — Bedgraph written non-atomically; partial write on crash leaves corrupt file — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze_command.py:1673`: `with open(output_path, 'w') as f: f.write(...)` — no temp-file + `os.replace()` pattern. A crash mid-write leaves a partial bedgraph that is indistinguishable from a valid file on restart.

**Fix:** Write to a `.tmp` file alongside the target; call `os.replace(tmp, output_path)` on success.
**File:** `rectify/core/analyze_command.py:1673`

---

### ~~NEW-052 (MEDIUM) — Comment says "modal isoform" but code uses `np.median()` — Fixed 2026-04-09 (v2.7.8)~~

`core/analyze/apa_detection.py:324-325`: comment reads `# Calculate modal TES position` but the implementation calls `int(np.median(data['tes_positions']))`. Modal (most frequent) and median (middle value) can differ substantially for multi-modal distributions.

**Fix:** Either change the code to use the actual mode (`scipy.stats.mode` or `Counter.most_common`) or change the comment to say "median."
**File:** `rectify/core/analyze/apa_detection.py:324-325`

---

### ~~NEW-059 (MEDIUM) — `Path.is_relative_to()` requires Python 3.9; package declares `python_requires >= 3.8` — Fixed 2026-04-09 (v2.7.8)~~

`utils/provenance.py:138`: `filepath.is_relative_to(self.output_dir)` — `Path.is_relative_to()` was added in Python 3.9. `pyproject.toml` declares `requires-python = ">=3.8"`. Users on Python 3.8 receive `AttributeError`.

**Fix:** Replace with `try: filepath.relative_to(self.output_dir); is_relative = True except ValueError: is_relative = False`, which works on Python 3.8+.
**File:** `rectify/utils/provenance.py:138`

---

### ~~NEW-060 (MEDIUM) — Full-length classifier minus-strand branch uses same field as plus-strand — Fixed 2026-04-09 (v2.7.8)~~

`core/classify/full_length_classifier.py:150-153`: both `if record.strand == '+'` and `else` branches assign `five_prime_pos = record.five_prime_corrected`. The minus-strand branch should use a different coordinate to represent the actual 5' end.

**Fix:** For minus strand, use the coordinate that represents the true 5' (TSS) end per the coordinate convention (`reference_end - 1`).
**File:** `rectify/core/classify/full_length_classifier.py:150-153`

---

## Fixed (v2.7.7 — 2026-04-08)

### ~~Bug 29 (CRITICAL) — Double-subtraction of 3'SS position in terminal exon refiner~~ — Fixed 2026-04-08 (v2.7.7)
Removed the erroneous `- 2` from `abs(read_end - pos - 2)` → `abs(read_end - pos)`.
**File:** `rectify/core/terminal_exon_refiner.py:1072`

---

### ~~Bug 30 (CRITICAL) — `use_scratch` and `streaming` silently dropped by `_apply_profile()`~~ — Fixed 2026-04-08 (v2.7.7)
Added `use_scratch`, `streaming`, `analyze_cpus`, `analyze_mem`, `analyze_time` to `_apply_profile()` dispatch table and defaults dict.
**File:** `rectify/core/batch_command.py`

---

### ~~Bug 31 (HIGH) — `run_command.py` doesn't record resolved config in provenance~~ — Fixed 2026-04-08 (v2.7.7)
`_run_junction_aggregation()` now passes `config=vars(args)` to `init_provenance()`.
**File:** `rectify/core/run_command.py`

---

### ~~Bug 33 (HIGH) — Minus-strand 5'SS truncation check uses wrong coordinate~~ — Fixed 2026-04-08 (v2.7.7)
Changed to `read_5prime = read_end - 1` for minus-strand reads in `detect_junction_truncated_reads()`.
**File:** `rectify/core/terminal_exon_refiner.py`

---

### ~~Bug 34 (HIGH) — `analyze_cpus/mem/time` not applied from SLURM profiles~~ — Fixed 2026-04-08 (v2.7.7)
Resolved by Bug 30 fix.
**File:** `rectify/core/batch_command.py`

---

### ~~Bug 35 (HIGH) — GFF-derived 3'SS positions not validated against genome~~ — Fixed 2026-04-08 (v2.7.7)
Added optional `genome` parameter to `load_splice_sites_from_gff()`; validates AG/CT dinucleotide at each 3'SS and skips sites that fail.
**File:** `rectify/core/terminal_exon_refiner.py`

---

### ~~Bug 39 (MEDIUM) — GC-AG missing from canonical splice motif set~~ — Fixed 2026-04-08 (v2.7.7)
Added `'GC-AG'` to `_CANONICAL_MOTIFS`.
**File:** `rectify/core/junction_validator.py`

---

### ~~Bug 40 (MEDIUM) — Position shift double-counting for A-tract reads~~ — Fixed 2026-04-08 (v2.7.7)
`total_position_shifts` no longer increments for reads already counted under `atract_ambiguity`.
**File:** `rectify/core/processing_stats.py`

---

### ~~Bug 42 (MEDIUM) — Jensen-Shannon divergence epsilon regularization biased~~ — Fixed 2026-04-08 (v2.7.7)
Replaced manual JSD with `scipy.spatial.distance.jensenshannon() ** 2`.
**File:** `rectify/core/analyze/shift_analysis.py`

---

### ~~Bug 43 (MEDIUM) — A-richness threshold inconsistency~~ — Fixed 2026-04-08 (v2.7.7)
Local `A_RICHNESS_THRESHOLD = 0.7` removed; now uses `POLYA_RICHNESS_THRESHOLD` from `config.py`.
**File:** `rectify/core/false_junction_filter.py`

---

### ~~Bug 44 (MEDIUM) — `gene_attribution.py` silently uses un-converted GFF coords~~ — Fixed 2026-04-08 (v2.7.7)
Added heuristic warning in `build_cds_interval_tree()` when `start` column minimum is 1.
**File:** `rectify/core/analyze/gene_attribution.py`

---

### ~~Bug 45 (MEDIUM) — Minus-strand ambiguity window extends in wrong direction~~ — Fixed 2026-04-08 (v2.7.7)
Corrected to `ambiguity_min = current_position - range`, `ambiguity_max = current_position`.
**File:** `rectify/core/bam_processor.py`

---

### ~~Bug 46 (MEDIUM) — BAM handle not in context manager in `find_coverage_gaps()`~~ — Fixed 2026-04-08 (v2.7.7)
Wrapped in `with pysam.AlignmentFile(bam_path, 'rb') as bam:`.
**File:** `rectify/core/bam_processor.py`

---

### ~~Bug 47 (MEDIUM) — Bare `except Exception: pass` swallows gene attribution errors~~ — Fixed 2026-04-08 (v2.7.7)
Changed to `except Exception as _e: logger.warning(...)`.
**File:** `rectify/core/bam_processor.py`

---

### ~~Bug 48 (MEDIUM) — YAML profile fields not validated at load time~~ — Fixed 2026-04-08 (v2.7.7)
`load_slurm_profile()` now warns on unrecognized keys.
**File:** `rectify/core/batch_command.py`

---

### ~~Bug 49 (MEDIUM) — `--continue-on-error` not wired to SLURM mode~~ — Fixed 2026-04-08 (v2.7.7)
Generated scripts now use `set -uo pipefail` (no `-e`) when `--continue-on-error` is set.
**File:** `rectify/core/batch_command.py`

---

### ~~Bug 50 (MEDIUM) — `--partition` not enforced when `--submit` is used~~ — Fixed 2026-04-08 (v2.7.7)
Returns an error immediately if `--submit` is requested without `--partition`.
**File:** `rectify/core/batch_command.py`

---

### ~~Bug 51 (MEDIUM) — Subprocess return code not checked for uLTRA/deSALT SAM→BAM~~ — Fixed 2026-04-08 (v2.7.7)
Both pipelines now raise `RuntimeError` on non-zero exit from `samtools view` or `samtools sort`.
**File:** `rectify/core/multi_aligner.py`

---

### ~~Bug 52 (MEDIUM) — SGE array task IDs are 1-based; generated script assumes 0-based~~ — Fixed 2026-04-08 (v2.7.7)
Added `$((SGE_TASK_ID - 1))` normalization in the generated scheduler abstraction block.
**File:** `rectify/core/batch_command.py`

---

### ~~Bug 53 (MEDIUM) — Exclusion regions: 1-based input detected but not rejected~~ — Fixed 2026-04-08 (v2.7.7)
Elevated to `ERROR`-level log with fraction of affected rows; stronger diagnostic message.
**File:** `rectify/core/exclusion_regions.py`

---

### ~~Bug 54 (MEDIUM) — Silent data loss: reads without gene attribution not logged in APA detection~~ — Fixed 2026-04-08 (v2.7.7)
Logs a WARNING when >10% of records lack gene attribution.
**File:** `rectify/core/analyze/apa_detection.py`

---

## Fixed (prior releases)

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

### ~~Bug 4: Duplicate primary records in rectified BAM~~ — Fixed 2026-04-01

Two root causes:
1. `_paf_to_bam()` never checked the `tp:A:` PAF tag — `tp:A:S` secondary
   records got written as primaries. Now `tp:A:S` → `FLAG |= 0x100`.
2. The rectified BAM writer did not track which read_ids already had a primary
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
