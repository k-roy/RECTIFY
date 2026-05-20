# RECTIFY Bugs to Fix

## Last Updated: 2026-05-20 (NEW-076 through NEW-081 added — bedgraph coordinate audit + adjacent coordinate-convention findings from multi-agent audit run)

---

## Open

---

### NEW-076 (HIGH) — `rectify analyze` per-condition bedgraph: 1-bp left shift in `start`/`end` BED columns; wider audit of other bedgraph emitters pending

**Files:** confirmed in `rectify/core/analyze/bedgraph.py:99-100` (FIXED 2026-05-20, uncommitted). Audit pending in `rectify/core/bam/bedgraph_writers.py`, `rectify/core/netseq/netseq_output.py`, `rectify/core/commands/{consensus,run,export,analyze,correct,split}_command.py`, `rectify/core/analyze/manifest.py`, `rectify/data/validation/generate_igv_html.py`.

**Symptom:** Per-condition strand-specific bedgraphs from `rectify analyze` show their 3'-end pileup peak 1 bp to the LEFT of the true position. Caught 2026-05-20 in the He et al cross-modality TRT analysis: CST6 ysh1 peak in `han2023_analyze_wbfix_20260512/analyze/bedgraph/Ysh1AA_plus.bedgraph` at `chrIX 287748 287749` (IGV 287,749) vs the matching DRS bedtools-derived 3'-end peak at `chrIX 287749 287750` (IGV 287,750). Affects DRS, QuantSeq REV, and PCR-cDNA equally — the emitter is protocol-agnostic.

**Root cause:** the writer at `analyze/bedgraph.py:99-100` did `start = pos - 1; end = pos`, treating `corrected_3prime` as 1-based. It is **0-based-inclusive** (from `reference_end - 1` / `reference_start`). Correct: `start = pos; end = pos + 1`. See `dev/audits/bedgraph_coordinate_audit_20260520.md` for the full convention statement and the per-file audit table.

**Fix status:**
- `analyze/bedgraph.py:99-100` corrected on M1 working tree (uncommitted); 3 regression tests added in `tests/test_analyze.py::TestBedgraphCoordinates`. Verified by toggling the fix off/on (FAIL/PASS).
- Other emitters listed in the audit doc are NOT YET INSPECTED. Same bug pattern could exist in any of them. The audit doc has a verification recipe (trace each `start`/`end` assignment to its source column; if 0-based, must be `[pos, pos+1)`).

**Blast radius:** every prior `rectify analyze` output's `bedgraph/` subdir has 1-bp-left-shifted per-condition tracks. Clustering, shift analysis, and per-position attribution work from the `corrected_reads.tsv` / `corrected_3ends.tsv` position columns directly and are unaffected; only the bedgraph tracks (used for IGV inspection and cross-track comparison) need regeneration after fix lands.

**See:** `AGENT_FIXES.md` `[2026-05-20]` entry for the first emitter's fix detail; `dev/audits/bedgraph_coordinate_audit_20260520.md` for the audit plan.

---

### NEW-077 (HIGH) — `bam_processor.py` minus-strand artifact N snap can emit an intronic / skipped coordinate as the 3' end

**File:** `rectify/core/bam/bam_processor.py:826-836`.

**Symptom:** When `current_position` falls inside an artifact N span `[junction_start, junction_end)` on the minus strand, the snap routine sets `current_position = junction_start`. Because `junction_start` is the first base of the SKIPPED (N-op) span in half-open coordinates, the emitted 3' end coordinate lands ON a base that the alignment doesn't cover. Pathogenic for reads whose walkback lands inside an artifact N.

**Likely fix direction:** the target on the minus strand should be the first ALIGNED base on the 3' side (i.e. `junction_end` if interpreting half-open, or the aligned base immediately upstream of the N span — depends on intended geometry). Needs a careful review of what "3' end of the read" means semantically when the polished position falls inside an intron/N span.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### NEW-078 (HIGH/MEDIUM) — `bam_processor.py` minus-strand artifact ambiguity clipping leaves NET-seq refinement window touching the N span

**File:** `rectify/core/bam/bam_processor.py:847`.

**Symptom:** Minus-strand ambiguity-range clipping does `ambiguity_min <= junction_start` then sets `ambiguity_min = junction_start`. Because `[start, end)` is half-open, this leaves the refinement window starting at the first SKIPPED base of the N span. Downstream NET-seq refinement can still assign reads into the intronic / N-spanned region.

**Likely fix direction:** clip `ambiguity_min` to the first ALIGNED base outside the N span, not to `junction_start` (which is the boundary of the skipped span). Related to NEW-077 — same root confusion between half-open boundary semantics and aligned-base semantics on the minus strand.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### NEW-079 (MEDIUM) — `clustering.py` plus-strand `distance_to_gene_3prime` is off by 1

**File:** `rectify/core/analyze/clustering.py:532`.

**Symptom:** Uses `gene['end']` as the plus-strand gene 3' coordinate. Loaded annotations are 0-based half-open (per `rectify/core/analyze/loaders.py:411` which does `start - 1, end`), so the last aligned base of the gene is `gene['end'] - 1`, not `gene['end']`. Result: `distance_to_gene_3prime` is 1 bp too large for plus-strand clusters.

**Blast radius:** affects the `distance_to_gene_3prime` column in cluster output only; does NOT affect clustering itself or any downstream call that uses cluster identity rather than gene-3' distance. Low biological impact but easy to clean up — `gene_3prime_pos = gene['end'] - 1` for plus, `gene['start']` for minus.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### NEW-080 (MEDIUM) — `false_junction_filter.py` minus-strand wrong-flank lookup against artifact N target

**File:** `rectify/core/splice/false_junction_filter.py:282` (comment) vs `:300` (code).

**Symptom:** The comment at line 282 says the minus-strand poly-A side is BEFORE the N (i.e. on the `[start-20, start)` flank), but the code at line 300 always fetches `[end, end+20)` regardless of strand. Plausibly wrong-flank logic on the minus strand. Feeds into the artifact N snap issue tracked as NEW-077.

**Likely fix direction:** branch on strand in the flank lookup; minus strand should fetch the upstream-of-N flank, not the downstream-of-N flank.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### NEW-081 (LOW) — `analyze_command.py` / `analyze/manifest.py` rDNA exclusion uses inclusive `<= end` against half-open annotations

**Files:** `rectify/core/commands/analyze_command.py:130`; `rectify/core/analyze/manifest.py:217`.

**Symptom:** Both compare `corrected_position <= end` to test rDNA-region membership. If `rdna_regions` are half-open `[start, end)` (consistent with `loaders.py:411`), the `<= end` predicate excludes one EXTRA base at the right edge of each region. Low biological impact (off-by-one at region boundaries; rDNA is large enough that one extra base on each side is noise), but trivial to fix to `< end`.

**Source:** caught by the user's multi-agent coordinate-convention audit, 2026-05-20.

---

### NEW-075 (HIGH) — `correct`: dT-primed cDNA walkback is not integrated; post-hoc script emits BAM strand instead of gene/RNA strand

**File:** missing in `rectify/core/commands/correct_command.py`; current workaround is the standalone script at `projects/TRT/scripts/rectify/han2023/11_polya_walkback_recompute.py`.

**Symptom:** For dT-primed cDNA libraries (e.g. QuantSeq REV, the Han 2023 reanalysis), the corrected 3' end needs a polyA walkback step that scans backward from the alignment 3' end past stretches of genome-encoded A's to find the templated CPA. `rectify correct` does not perform this step. Users currently run `11_polya_walkback_recompute.py` post-hoc against the `corrected_3ends.tsv`. That standalone script has a strand bug: it emits the **BAM read strand** in the strand column, but for dT-primed antisense protocols the BAM strand is the **opposite** of the RNA/gene strand. Downstream clustering then sees flipped strands at convergent loci and the entire cluster taxonomy gets confused (caught 2026-05-12 during the Han 2023 cluster-revisit, fixed locally via a TSV strand-flip).

**Fix:** Integrate polyA walkback into `rectify correct` behind a `--dT-primed-cDNA` flag (or auto-detect from sample-sheet protocol field). Two requirements for the integration:

1. **Walkback algorithm:** scan upstream from the 3' end of each alignment past genome-encoded A-runs (parameter `--walkback-min-A-run`, default 4) and emit the templated CPA position. This is the existing logic in `11_polya_walkback_recompute.py` — port and harden.
2. **Strand convention:** the strand column in `corrected_3ends.tsv` must always be the **gene / RNA strand**, never the BAM-read strand. For dT-primed antisense protocols, flip the strand when emitting (BAM `is_reverse` → `+` RNA strand, and vice versa). Add a test fixture covering this case to prevent regression.

**Why HIGH:** the post-hoc workaround is fragile (lives in TRT-specific scripts, not rectify), the strand bug silently corrupts cluster output at convergent loci, and any future user reanalysing dT-primed cDNA libraries will hit the same trap.

**Discovered:** 2026-05-12 during Han 2023 cluster revisit; documented in `projects/TRT/analyses/han2023_cluster_revisit_20260512/`.

---

### NEW-074 (MEDIUM) — `consensus.py`: `shutil.copy2()` to NFS without `os.fsync()` risks partial writes

**File:** `rectify/core/consensus/consensus.py` (~line 1761)

**Symptom:** Sorted BAM copied from local scratch to NFS Oak with `shutil.copy2()`. Under NFS load, OS buffers the write; if the job crashes or the NFS connection drops before the buffer flushes, the destination BAM is silently truncated. Downstream readers crash with "truncated file".

**Fix:**
```python
shutil.copy2(src, dst)
with open(dst, 'r+b') as f:
    os.fsync(f.fileno())
```
Apply to both the sorted BAM copy and the `.bai` copy at the same location. Alternatively, switch to `rsync` (already used in `sync_to_oak()`) which performs checksum-verified transfer.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### NEW-073 (MEDIUM) — `bam_processor.py`: 17× redundant full BAM scans before correction begins

**File:** `rectify/core/bam/bam_processor.py` (lines 1335, 1352, 1768)

**Symptom:** Three independent full BAM scans run before any correction worker starts:
1. Line 1335: full linear scan for read-flag stats (`ProcessingStats`)
2. Line 1352: calls `get_processing_regions()` → `find_coverage_gaps()` → one `bam.fetch(chrom)` scan per chromosome (16 scans for yeast)
3. Line 1768: another full linear scan for the variant-aware prescan

For a 200k-read yeast BAM this consumes ~47 min of wall time before correction starts, contributing to the observed 6h timeout.

**Fix:** Add `prescan_bam_once(bam_path)` that performs a single `bam.fetch(until_eof=True)` pass, simultaneously accumulating flag stats, building per-chromosome coverage intervals for gap-finding, and optionally running `VariantAwareHomopolymerRescue.scan_read()` when variant-aware mode is enabled. Call this once before spawning workers; remove the redundant scans.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279); estimated ~47 min wall time savings per sample

---

### NEW-072 (MEDIUM) — `multi_aligner.py`: deSALT sort subprocess `stderr=PIPE` with no reader thread risks deadlock

**File:** `rectify/core/align/multi_aligner.py` (lines 1154–1161, `run_desalt()`)

**Symptom:** `sort_proc` spawned with `stderr=PIPE` (or inheriting parent stderr). If samtools sort writes ≥64 KB of warnings, the OS pipe buffer fills; since no thread drains it, `sort_proc.communicate()` deadlocks. (The minimap2 path already has a drain thread at lines 209–216; deSALT does not.)

**Note:** mapPacBio was already fixed (`stderr=subprocess.DEVNULL` with explanatory comment). This affects deSALT and uLTRA only.

**Fix:** Add a daemon stderr-drain thread matching the minimap2 pattern (lines 209–216):
```python
desalt_stderr_chunks = []
def _drain_stderr():
    desalt_stderr_chunks.append(sort_proc.stderr.read())
drain_thread = threading.Thread(target=_drain_stderr, daemon=True)
drain_thread.start()
sort_proc.wait()
drain_thread.join()
```

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### NEW-071 (MEDIUM) — `run_command.py`: `stderr=DEVNULL` on samtools fastq hides error context

**File:** `rectify/core/commands/run_command.py` (line 751, DRS trim → FASTQ conversion)

**Symptom:** `subprocess.run(['samtools', 'fastq', ...], check=True, stderr=subprocess.DEVNULL)`. If samtools fails, `check=True` raises `CalledProcessError` (error is caught), but all samtools stderr is discarded, making the failure impossible to diagnose. Line 1157 (same operation in single-sample path) correctly omits `DEVNULL` — only line 751 is affected.

**Fix:** Remove `stderr=subprocess.DEVNULL`. Capture or let stderr propagate; `check=True` already raises on non-zero exit.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### NEW-070 (MEDIUM) — `run_command.py`: non-atomic `unlink()` + `rename()` instead of `Path.replace()`

**File:** `rectify/core/commands/run_command.py` (lines 878–879, 1403–1404)

**Symptom:** Two locations use:
```python
_restored_bam.unlink(missing_ok=True)
Path(_restored_tmp).rename(_restored_bam)
```
If a crash or SIGKILL occurs between the two operations, `_restored_bam` is permanently absent (data loss). On NFS, `rename()` across directories can also fail.

**Fix:** Replace with the atomic POSIX `replace()` call:
```python
Path(_restored_tmp).replace(_restored_bam)
```
`Path.replace()` is atomic on POSIX filesystems and handles the unlink+rename as a single syscall.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### NEW-069 (MEDIUM) — `consensus.py`: `.done` checkpoint sentinel written without BAM validity check

**File:** `rectify/core/consensus/consensus.py` (lines 1667, 1698)

**Symptom:** After each batch BAM is written and `out_bam.close()` is called, the `.done` sentinel file is written immediately. If the batch BAM is corrupt (CIGAR/seq mismatch, partial write), the sentinel marks it complete. On resume, `pysam.cat()` or `pysam.sort()` reads the corrupt batch BAM and crashes.

**Fix:** Before writing the sentinel, call `_validate_bam_cigar(_cur_batch_path, sample_size=1000)` (the same function added for pre-sort validation). If validation fails, log the error and do NOT write the sentinel — the batch will be re-written on resume. A simple file-size check (`os.path.getsize > 0`) is insufficient: CIGAR-corrupted BAMs have non-zero size.

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### NEW-068 (MEDIUM) — `multi_aligner.py`: no assertion that minimap2 succeeded after aligner loop

**File:** `rectify/core/align/multi_aligner.py` (after aligner loop, ~line 1267)

**Symptom:** When all aligners are tried in a loop, exceptions are caught and logged per-aligner (`logger.error(f"Aligner {aligner} failed: {e}")`), but the function returns `results` with however many aligners succeeded — including zero. The caller receives an empty dict with no indication that minimap2 (the required baseline aligner) failed. Downstream consensus selection runs with no aligners and silently produces no output.

**Fix:** After the aligner loop, assert that minimap2 succeeded:
```python
if 'minimap2' not in results:
    raise RuntimeError(
        f"minimap2 (required baseline aligner) failed — "
        f"cannot proceed with consensus selection. "
        f"Succeeded: {list(results.keys())}"
    )
```

**Discovered:** 2026-04-22 post-mortem audit (job 22314279)

---

### NEW-067 (MEDIUM) — `align_command.py`: `calmd_bam.replace(rectified_bam)` before index creation leaves BAM without `.bai` on index failure

**File:** `rectify/core/commands/align_command.py` (lines 576–579)

**Symptom:** Current code:
```python
calmd_bam.replace(rectified_bam)          # destructive rename
_sp.run(['samtools', 'index', str(rectified_bam)], check=True)   # might fail
```
If `samtools index` fails (disk full, permissions, samtools crash), `rectified_bam` exists without a `.bai`. Downstream tools requiring an index fail with "index not found", and the skip-logic in `run_command.py` (now using `_validate_bam_integrity()`) will correctly detect the missing `.bai` and re-run alignment — but the original `calmd_bam` has been destroyed, so alignment must restart from scratch.

**Fix:** Create the index first on the temp BAM, then do two atomic renames:
```python
if result.returncode == 0 and calmd_bam.stat().st_size > 0:
    _sp.run(['samtools', 'index', str(calmd_bam)], check=True)
    calmd_bam.replace(rectified_bam)
    Path(str(calmd_bam) + '.bai').replace(Path(str(rectified_bam) + '.bai'))
    logger.info("  MD tags added successfully")
```

**Discovered:** 2026-04-22 post-mortem audit (job 22314279); not present in squishy-waddling-stream.md

---

### NEW-066 (MEDIUM) — uLTRA reuses stale cached database with empty genome, silently produces zero alignments

**File:** `rectify/core/align/multi_aligner.py`, `run_ultra()` (~line 848)

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

**Files:** `rectify/core/bam/bam_processor.py` (`clip_read_to_corrected_3prime`), `rectify/core/splice/false_junction_filter.py`
**Discovered:** 2026-04-09 validation audit (cat4_plus_1, cat4_plus_2)

---

### ~~NEW-062 (MEDIUM) — `five_prime_rescued` not written to TSV; `correction_applied` omits 5' rescue — Fixed 2026-04-09 (v2.7.9)~~

When `rescue_3ss_truncation` corrects a read's 5' end, `five_prime_rescued=True` is stored in the result dict (bam_processor.py:277) but never written to any TSV column. The `correction_applied` field for cat3 reads shows only `atract_ambiguity,polya_walkback` — the 5' junction rescue is invisible to downstream consumers. Users cannot identify rescued reads from the TSV without comparing `five_prime_position` against raw alignment coordinates.

**Fix:** Add `five_prime_rescued` (bool) to the TSV header and row in `write_results_tsv`, and append `five_prime_rescued` to `correction_applied` when True.

**Files:** `rectify/core/bam/bam_processor.py:567–613` (`write_results_tsv`)
**Discovered:** 2026-04-09 validation audit (cat3 reads)

---

### ~~NEW-063 (MEDIUM) — Chimeric reads with 3'-hard-clipped alignment not exempt from poly-A walkback — Fixed 2026-04-09 (v2.7.9)~~

Reads with XK=1 (chimeric reconstruction) that have a hard-clip at the 3' end are subject to the full polya_walkback pipeline. For cat5_plus_3aligner (100H hard-clip), the unresolved 3' end causes an ambiguity window of 211 bp (positions 9753–9964), producing 5 TSV rows with walkback shifts of 6–185 bp. The correction pipeline has no short-circuit for chimeric reads whose 3' sequence is not present in the query.

**Fix:** In `bam_processor.correct_read_3prime()`, skip polya_walkback (or cap ambiguity window to 0) when the read has a hard-clip at its 3' end (indicates missing sequence; walkback is unanchored).

**Files:** `rectify/core/bam/bam_processor.py` (`correct_read_3prime`), `rectify/core/polya/atract_detector.py`
**Discovered:** 2026-04-09 validation audit (cat5_plus_3aligner)

---

### ~~NEW-064 (LOW) — `netseq_refinement` listed in `correction_applied` when ambiguity_range=1 — Fixed 2026-04-09 (v2.7.9)~~

When `ambiguity_range == 1`, there is only one candidate position; NET-seq refinement runs but is a no-op. Listing `netseq_refinement` in `correction_applied` is misleading — it implies the NET-seq signal was used to resolve ambiguity when no ambiguity existed. Seen in cat6_minus_single (ambiguity_range=1, fraction=1.0, correction_applied includes `netseq_refinement`).

**Fix:** In `bam_processor.correct_read_3prime()`, only add `netseq_refinement` to `correction_applied` when `ambiguity_range > 1` (i.e., when NET-seq signal was actually consulted to break a tie).

**Files:** `rectify/core/bam/bam_processor.py` (correction_applied assembly)
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
**Files:** `rectify/core/bam/bam_processor.py`, `rectify/core/commands/correct_command.py`, `rectify/core/commands/tag_polya_command.py`, `rectify/core/align/preprocess.py`, `rectify/cli.py`

---

### ~~Bug 55 (MEDIUM) — Several APA clustering parameters not CLI-configurable — Fixed 2026-04-11 (v0.9.0)~~

`DEFAULT_MIN_PEAK_SEPARATION = 5`, `DEFAULT_MAX_CLUSTER_RADIUS = 10`, and `DEFAULT_MIN_SAMPLES = 2` had no corresponding CLI arguments.

**Fix:** Added `--min-peak-sep`, `--max-cluster-radius`, `--min-cluster-samples` to `create_analyze_parser()`. Both clustering call sites (single-sample and manifest mode) dispatch to `cluster_cpa_sites_adaptive()` when any non-default value is provided; otherwise fall through to the existing fixed-distance `cluster_cpa_sites()`.
**Files:** `rectify/core/commands/analyze_command.py`

---

## Historical bug fixes

Resolved bugs are no longer tracked here. See `../CHANGELOG.md` for the
release-by-release history.
