# AGENT_FIXES.md

Fast coordination log for active debugging sessions across M1 / H2 / Sherlock.
**Read this before touching pipeline code. Update it when you find a bug.**
Archive entries into CHANGELOG.md when the session wave is done.

---

## [2026-05-21] BAM-first 10k correct smoke + RN sidecar check

**Status:** FIRST PASS SUCCEEDED on H2 for TSV correction; follow-up 2H-enabled
pass needed because the first command omitted `--aligner-bams`.

**Why this was run:** the old chunk directories often have BAMs but zero-byte
`corrected_reads.tsv` files, so smoke-test candidates must be identified from
the BAMs themselves rather than from successful `rectify correct` outputs.

**Input panel:** H2
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/aligner_chunks`,
chunk `chunk_000`, all five aligners.

**BAM-derived candidate scan highlights:**
- Existing per-aligner chunk BAMs were present for `mapPacBio`, `minimap2`,
  `gapmm2`, `uLTRA`, and `deSALT`; their `corrected_reads.tsv` files were all
  0 bytes.
- The scan selected 10,000 normalized read IDs from BAM CIGAR/read-pattern
  triggers, including junction/N-op reads, terminal clips, terminal
  clip+junction combinations, short terminal exons, and terminal indel-rich
  reads.
- Reads shared by all five aligners in the source BAM panel: 346,618.

**Subset BAMs created:**
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/aligner_chunks/{aligner}/chunk_000/*.subset_10k.bam`

Primary read IDs retained:
- `mapPacBio`: 9,491 primary read IDs / 10,143 records
- `minimap2`: 8,704 primary read IDs / 10,168 records
- `gapmm2`: 8,366 primary read IDs / 8,372 records
- `uLTRA`: 8,697 primary read IDs / 8,706 records
- `deSALT`: 9,788 primary read IDs / 12,397 records

**First `rectify correct` smoke run:** H2 SGE array `13465761.1-5` succeeded
for all five aligners with `--streaming --emit-merged-tsv`, no corrected BAM
write, and no heap corruption. Outputs are under:
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/correct_outputs/{aligner}/chunk_000/`

Timings:
- `mapPacBio`: 9,504 reads processed; BAM processing 25.2s; total 27.4s
- `minimap2`: 8,707 reads processed; BAM processing 23.7s
- `gapmm2`: 8,372 reads processed; total completed successfully
- `uLTRA`: 8,700 reads processed; total completed successfully
- `deSALT`: 9,788 reads processed; total 33.3s

**Important caveat:** this first smoke run exercised 2F 3'SS rescue and core
correction, but Module 2H junction refinement was logged as disabled:
`Junction refinement: DISABLED (pass --aligner-bams to enable)`. A follow-up
10k run with the same subset BAM panel passed as repeated `--aligner-bams`
is needed before treating this as a full 2F+2H correction smoke.

**RN sidecar test:** wrote a post-hoc BAM-panel sidecar:
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/wt_tfiiib_rep2_chunk000_bam_panel_10k.read_num_sidecar.parquet`

Provenance:
`/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep2/chunks/correct_smoke_10k_from_bam/wt_tfiiib_rep2_chunk000_bam_panel_10k.read_num_sidecar.POSTHOC_PROVENANCE.json`

Validation:
- Sidecar rows: 10,000; missing selected reads: 0.
- Normalized QNAME lookup covered every primary record in all five subset
  BAMs: `mapPacBio` 9,504/9,504, `minimap2` 8,707/8,707, `gapmm2` 8,372/8,372,
  `uLTRA` 8,700/8,700, `deSALT` 9,788/9,788.
- Fingerprint verification is intentionally not a pass/fail criterion for this
  BAM-derived sidecar because aligners can hard-clip/trim/reorient query
  sequence and qualities. For production provenance, sidecars reconstructed
  from chunk FASTQs remain the correct fingerprint source.
- Local focused module tests passed:
  `tests/test_read_num_sidecar.py tests/test_consensus_tag_restoration.py`
  → `7 passed, 1 skipped`. H2 production env lacks `pytest`, so the cluster
  validation was the live sidecar open/lookup test above.

## [2026-05-21] Module 2H junction pool explosion — CANDIDATE LOOKUP + YEAST SIZE CAP IN WORKING TREE

**Status:** FIXED in working tree, uncommitted. Needs deployment to Sherlock/H2
RECTIFY checkouts before rerunning production chunk-correct jobs.

**Scope:** DRS/CPA `rectify correct` Module 2H junction refinement, especially
all-5-aligner per-chunk correction on Sherlock where `refine_bam_junctions`
was spending tens of minutes after loading a 53k min-support-filtered pool and
where the unfiltered all-5 pool reached ~5.0M junctions for `ysh1_rep2`.

**Root cause clarified:**
- Candidate lookup was already restricted to the same chromosome.
- There was a `--junction-search-radius` default of 5000 bp, but the active
  lookup call also used `start_radius=end_radius=--junction-max-boundary-shift`
  (default 50 bp). That means candidates were endpoint-bounded to ±50 bp, not
  a free 5 kb window.
- `_candidates_near()` still scanned the per-chromosome sorted list from the
  chromosome start for every N-op until `js > ns + radius`. With large noisy
  pools, this made the cost grow with chromosome-position × N-op count.
- Module 2H had no explicit max junction/intron-size cap in pool construction
  or candidate scoring. Yeast should use an organism-tuned cap; 10 kb is a
  conservative S. cerevisiae value.

**Fixes landed:**
- `junction_scoring.py`: `_candidates_near()` now uses `bisect_left` /
  `bisect_right` to scan only the intron-start window around the current N-op
  instead of repeatedly walking from chromosome start.
- `junction_scoring.py`: added optional `max_junction_size` filtering to
  `collect_junctions_from_bam()`, `collect_junction_counts_from_bam()`,
  `build_junction_pool()`, and `_candidates_near()`.
- `junction_scoring.py`: `build_junction_pool()` now supports
  `min_observed_support` and counts N-op support across aligner BAMs; annotated
  junctions are retained regardless of support or max-size cap.
- `junction_scoring.py`: pool scans skip secondary and supplementary records.
- `junction_refiner.py`: `max_junction_size` is propagated through
  `refine_bam_junctions()`, sequential/parallel workers, and
  `refine_read_junctions()`, so it applies both to newly built pools and to
  pre-built pools loaded from `rectify prescan`.
- `correct_command.py`: added `--junction-max-size BP`.
- `prescan_command.py`: added `--junction-min-support N` and
  `--junction-max-size BP`; both are recorded in `junction_pool.pkl` metadata.

**Recommended yeast command policy:**
- Prefer preserving all aligners, but pass `--junction-max-size 10000` to both
  `rectify prescan` and per-chunk `rectify correct`.
- Keep `--junction-min-support` available as an optional guardrail, but the
  10 kb cap plus bounded lookup may make min-support unnecessary for yeast.
- If using a pre-built pool, still pass `--junction-max-size 10000` to
  `rectify correct`; it filters oversized candidates at scoring time even if
  the pickle was built before the cap existed.

**Tests run:** focused splice/refiner suite passed with single-threaded BLAS
environment variables:
```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
pytest tests/test_junction_scoring_parallel.py tests/test_junction_refiner.py
```
Result: `41 passed, 17 skipped, 1 warning in 4.41s`.

**Caveat:** running the same pytest command without BLAS thread caps aborted
during NumPy's macOS import/runtime check before test collection. This appears
environmental and not specific to the junction-refiner patch.

**Operational note:** Sherlock test array `25600177_[0-4]` was still RUNNING at
25 min on the old code path after loading `junction_pool.min10.all5.pkl`
(53,392 total junctions). Replace/rerun the single-chunk all-5-aligner test
after syncing this patch, using `--junction-max-size 10000`.

**[2026-05-21 09:38 PDT] Cluster sync + first speed check:**
- M1 working-tree code/docs/tests were synced to Sherlock and H2 so the three
  checkouts are consistent. Large `rectify/data/` payloads were excluded.
- Backups were written before sync:
  - Sherlock: `.codex_backups/pre_m1_code_sync_20260521_092607`
  - H2: `.codex_backups/pre_m1_code_sync_20260521_092608`
- A first broad rsync briefly copied package/test/docs contents into the
  checkout root. Those accidental top-level copies were moved, not deleted:
  - Sherlock: `.codex_backups/accidental_top_level_20260521_093324`
  - H2: `.codex_backups/accidental_top_level_20260521_093330`
- Verification:
  - Sherlock production Python: `tests/test_junction_scoring_parallel.py`
    + `tests/test_junction_refiner.py` → `41 passed, 17 skipped`.
  - H2 production Python (`/u/project/guillom/shared/envs/rectify/bin/python`):
    import smoke passed and `rectify correct --help` shows
    `--junction-max-size`.
- Patched Sherlock benchmark `25601448_[0-4]` is running with the previous
  min10 all-aligner pool and `--junction-max-size 10000`.
  At 09:38 PDT, tasks were still in Module 2H:
  `mapPacBio` 12:15 elapsed, other aligners 7:35 elapsed, no outputs yet.
  That is not convincing speedup evidence. Since the 10 kb cap removes only
  1,342 / 53,392 junctions from the already-min10 pool (~2.5%), this benchmark
  mostly tests the bisect lookup. If it remains slow, the next bottleneck is
  per-candidate scoring and we should add scoring/memoization instrumentation
  rather than relying on lookup changes alone.

**[2026-05-21 11:27 PDT] Speed check result — not enough:**
- Patched benchmark `25601448_[0-4]` was stopped after proving the current
  optimization is insufficient. No final corrected TSV/BAM outputs were
  produced.
- Module 2H did complete for several aligners, but still took about an hour:
  - mapPacBio: 402,070 reads, 205,619 N-op reads, 101,183 refined;
    Module 2H timing `3915.4s` (~65.3 min), then timed out in downstream BAM
    parallel streaming before final output.
  - gapmm2: 330,634 reads, 211,002 N-op reads, 81,214 refined;
    Module 2H timing `3607.7s` (~60.1 min), then remained in downstream BAM
    streaming.
  - uLTRA: 340,584 reads, 209,726 N-op reads, 85,200 refined;
    Module 2H timing `3547.4s` (~59.1 min), then remained in downstream BAM
    streaming.
  - minimap2 was still in Module 2H at ~54 min when the benchmark was stopped.
  - deSALT log ended with a Python memory-map dump while still running; treat
    this test as failed/indeterminate.
- Interpretation: bisect lookup and a 10 kb max-size cap are necessary cleanup,
  but do not make all-5-aligner Module 2H tractable on this chunk. The dominant
  cost is now per-read/per-candidate scoring and/or too many N-op reads being
  refined. Next plan should instrument candidate counts and add a hard
  candidate cap/support-ranked candidate selection or a much narrower yeast
  junction pool before running production correction.

**[2026-05-21] Plan C implementation — profiling + candidate cap:**
- Added opt-in Module 2H profiling:
  - `rectify correct --junction-profile profile.json`
  - `--junction-profile-sample-rate N` (default 1)
- Profile JSON records aggregate timings, counts, raw/after-cap candidate
  histograms, p50/p90/p99/max candidate counts, `_score_junction` calls,
  `_score_hp_anchored` tier-1/tier-2 calls, boundary-filter counts, and final
  Module 2H stats.
- Added deterministic candidate cap:
  - `--junction-max-candidates-per-nop N`
  - ranking order: current junction first, annotated junctions next, then
    smallest boundary delta, smallest intron-length delta, coordinates.
  - Recommended first speed test value: `32`.
- Normal runs are unchanged unless either new flag is supplied.
- Focused local validation:
  ```bash
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  pytest tests/test_junction_refiner.py tests/test_junction_scoring_parallel.py -q
  ```
  Result: `43 passed, 17 skipped`.
- Synced to Sherlock and H2. Sherlock focused tests passed
  (`43 passed, 17 skipped`); H2 production Python import smoke passed.
- Sherlock Plan C benchmark submitted:
  - job array: `25624678_[0-4]`
  - sample/chunk: `ysh1_rep2`, chunk `000`
  - pool: previous `junction_pool.min10.all5.pkl`
  - flags: `--junction-max-size 10000`,
    `--junction-max-candidates-per-nop 32`,
    `--junction-profile <aligner>/chunk_000/junction_profile.json`
  - output root:
    `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/ysh1_rep2/chunks/planc_test`

**[2026-05-21] Bypass test hooks:**
- Added `rectify correct --skip-junction-refinement` to bypass Module 2H
  N-op boundary realignment/rewrite while still allowing an optional junction
  pool to be indexed for downstream Module 2F.
- Added `rectify correct --skip-3ss-rescue` to bypass Module 2F 3'SS
  truncation rescue entirely, including annotation, pool, and read-own N-op
  rescue triggers. Annotation-dependent gene attribution remains available.
- Planned benchmark: ysh1_rep2 chunk 000 all five aligners with both flags,
  to measure baseline correction runtime when both junction-specific rescue
  layers are disabled.

**[2026-05-21] ysh1 chunk 000 TSV-only vs corrected-BAM timing estimate:**
- In the no-2F/no-2H bypass benchmark, TSV correction finished quickly for
  the four usable aligners:
  - mapPacBio: BAM processing `93.0s`
  - minimap2: BAM processing `86.3s`
  - gapmm2: BAM processing `92.7s`
  - uLTRA: BAM processing `90.6s`
- The only fully completed corrected-BAM write was minimap2:
  `Corrected BAM write: 111.6s`, `Correction total: 202.6s`.
- mapPacBio/gapmm2/uLTRA had already emitted their corrected TSVs but were
  still in corrected-BAM writing when canceled at `14:35` elapsed, implying a
  lower bound of roughly `>12.6 min` extra BAM-writing time for each of those
  aligners on this chunk.
- Operational conclusion: corrected TSV manifests are sufficient for
  `rectify analyze`, so production correction should omit corrected BAMs by
  default. Generate corrected BAMs only for targeted IGV/debug chunks until
  the corrected-BAM writer is profiled/optimized.

**[2026-05-21] Lazy corrected consensus design note:**
- Dedicated plan written at `dev/specs/lazy_corrected_consensus_plan.md`.
- Core idea: keep homopolymer-aware corrected-CIGAR scoring and still emit a
  final corrected consensus BAM, but stop materializing full per-aligner
  corrected BAMs as an intermediate. Instead, score raw BAM records in memory
  using the per-aligner corrected TSV rows, select `winning_aligner`, then
  write only the winning corrected records into `corrected_consensus.bam`.
- Advisor review confirmed the main coupling points:
  `merge_corrected_tsvs()` uses per-aligner corrected BAMs only for
  `hp_edit_distance`/`aligned_bases`; generated split scripts use corrected
  BAMs both for that scoring and as the source for `corrected_consensus.bam`;
  non-chunked `run/stages.py` also unconditionally wires corrected-BAM output.
- Implementation should first target `realign_exon_blocks()` so it realigns
  only exon blocks containing homopolymer-position `X` ops, then extract the
  corrected-read edit sequence into a shared helper used by sequential writer,
  parallel writer, lazy HP scoring, and final consensus BAM generation.
- Important legacy-run caveat added to the plan: CPA DRS and H2 mex67aa
  predate RN-tagged chunk FASTQs / `*.read_num_sidecar.parquet`. Before lazy
  raw-BAM scoring is trusted on those outputs, reconstruct post-hoc read_num
  sidecars from existing chunk FASTQs using old round-robin chunk order
  (`read_num = local_index * n_chunks + chunk_index`) and write explicit
  `POSTHOC_PROVENANCE.json`. Use sidecar-backed normalized QNAME lookup first;
  inject `RN:i` into existing BAMs only if QNAME fallback coverage is
  insufficient.

**[2026-05-21] Lazy corrected consensus implementation checkpoint:**
- `realign_exon_blocks()` now targets only exon blocks that contain an `X` op
  at a homopolymer reference position, instead of realigning every eligible
  short exon block after a read-level trigger.
- Added shared `apply_corrected_edits_to_read()` in `bam_writer.py`; sequential
  and parallel corrected-BAM writers now use the same hard-clipped correction
  edit sequence that lazy scoring/final consensus BAM writing uses.
- `merge_corrected_tsvs()` now supports `per_aligner_raw_bams=...` for lazy
  HP-aware edit-distance scoring from raw BAM + corrected TSV, with strict
  identity checks by default. The old `per_aligner_corrected_bams` path remains
  supported and takes precedence when supplied.
- Added `write_corrected_consensus_bam()` to write the final corrected consensus
  BAM directly from raw BAMs + per-aligner corrected TSVs + merged
  `winning_aligner`, without materializing full per-aligner corrected BAMs.
- Generated chunk correction scripts now omit per-aligner
  `--write-corrected-bam` by default; chunk merge prefers lazy HP scoring and
  writes one `corrected_consensus.bam`.
- Added `reconstruct_posthoc_sidecar_from_chunks()` in
  `core/chunking/sidecar.py` for legacy round-robin chunk FASTQs, emitting
  `<sample>.POSTHOC_PROVENANCE.json`.
- Focused local validation:
  - `tests/test_bam_writer_parallel_smoke.py`,
    `tests/test_corrected_consensus_tiebreaker.py`,
    `tests/test_qname_sanitizer_and_validator.py`,
    `tests/test_read_num_sidecar.py`, `tests/test_splice_junction.py`,
    `tests/test_parallel_processing.py` →
    `145 passed, 1 skipped, 1 xfailed`.
  - Full `tests/test_validation_reads.py` did not run in the local sandbox:
    subprocess `rectify correct` calls aborted at OpenMP startup with
    `OMP: Error #179: Function Can't open SHM failed`, before exercising the
    assertions.
- Synced this implementation checkpoint to Sherlock and H2. Compile smoke
  passed on both. Sherlock focused test pass
  (`tests/test_bam_writer_parallel_smoke.py tests/test_read_num_sidecar.py`)
  passed: `11 passed`.

---

## [2026-05-20] CODEX audit bug sweep — CHECKPOINT/IO/PRESCAN FIXES IN WORKING TREE

**Status:** FIXED in working tree, uncommitted.

**Scope:** Follow-up fixes from `dev/BUGS_TO_FIX.md` NEW-067 through NEW-074 and
related `CODEX_AUDIT.md` durability/performance findings.

**Fixes landed:**
- `consensus.py`: checkpoint batches now fsync/validate BAMs before atomic
  `.done` sentinels; resume trusts only a contiguous valid batch prefix.
- `consensus.py`: scratch-to-output BAM and `.bai` copies now fsync destination
  files and parent dirs.
- `single_sample.py`: DRS `samtools fastq` failures preserve stderr/stdout in the
  raised error instead of discarding context.
- `single_sample.py`: DRS restored-polya BAM sort swaps now use `Path.replace()`
  instead of `unlink()` + `rename()`.
- `multi_aligner.py`: if minimap2 was requested but failed, the pipeline raises
  instead of returning a partial/empty aligner map into consensus.
- `align_command.py`: MD-tagged `.md.bam` is indexed before replacing
  `rectified.bam`; the `.bai` is atomically swapped only after index success.
- `bam/parallel.py`: parallel correction prescans the BAM once for read stats,
  coverage-region planning, and optional variant-aware rescue, removing the
  previous stats scan + per-chromosome coverage scans + variant scan fan-out.

**Tests run:** focused suite passed:
`tests/test_consensus_checkpoint_safety.py`, `tests/test_align_command_index_commit.py`,
`tests/test_multi_aligner_failures.py`, `tests/test_run_single_sample_safety.py`,
`tests/test_parallel_processing.py`, `tests/test_corrected_consensus_tiebreaker.py`,
`tests/test_analyze.py`, `tests/test_bam_writer.py`,
`tests/test_bam_writer_parallel_smoke.py` (`93 passed`). Earlier related slices
also passed for GTF feature expansion, checkpoint resume, and BAM writer smoke.

**Caveat:** full repository pytest has not been run in this session.

---

## [2026-05-20] mapPacBio QNAME sanitizer — ONT Dorado FASTQs

**Status:** FIXED (`e8c8070`)

**Affects:** any dataset sequenced with ONT Dorado (mex67aa, wtaa, and any
future DRS/cDNA deposits from the Nanopore). Does NOT affect cDNA pipeline
FASTQs produced by `rectify trim-polya` (those already have bare UUID QNAMEs).

**Symptom 1 — before fix 1:**
```
samtools view exit code 1: [E::sam_parse1] query name too long
mapPacBio failed after ~14400s
```
Dorado embeds full run metadata in the FASTQ header
(`@<uuid> runid=... ch=... flow_cell_id=... basecall_gpu=...`).
mapPacBio.sh copies the full header verbatim into SAM QNAME (346 chars);
SAM spec limit is 254. `samtools view -bS` exits 1 on every read.

**Symptom 2 — introduced by fix 1 (`838293c`), fixed in `382fcc7`:**
```
BBMap AssertionError: Error in mpb_san.fastq, line N, sequence line is blank
BBMap terminated in an error state
```
`split(' ', 1)[0]` on a no-space header (e.g. `__mpbsplit_*` sub-reads)
retains the trailing `\n`, so appending `\n` produced a double-newline
(blank sequence line).

**Symptom 3 — introduced by fix 2 (`382fcc7`), fixed in `e8c8070`:**
```
BBMap validateQualityLength: quality string length != sequence length
This can be bypassed with the flag 'tossbrokenreads' or 'nullifybrokenquality'
```
Line-by-line `startswith('@')` mis-identified quality score lines beginning
with `@` (valid FASTQ — Phred Q31) as headers and truncated them.

**Root fix (`e8c8070`):** parse sanitized FASTQ in 4-line blocks; only the
header line (block line 1) is ever touched. Quality lines are written as-is
regardless of content.

**Safe to pull:** strict correctness improvement; any FASTQ that worked
before still works.

**Downstream symptom in `rectify correct` (pre-fix):** if you didn't run
alignment yourself but inherited mapPacBio BAMs from a pre-`e8c8070` run,
`rectify correct` workers will hang/crash in the parallel-correction phase
even with the fix in place — because the bad BAMs already exist on disk.
Telltale:
```
INFO - rectify.core.bam.parallel - X regions across 8 workers
*** Error in `python': free(): invalid next size (fast): ... ***
*** Error in `python': free(): corrupted unsorted chunks: ... ***
```
followed by silent hang (Slurm reports RUNNING until walltime kill).
pysam/htslib's C code corrupts its buffer state when fetch() hits a
malformed QNAME or quality-length mismatch in the BAM. The bug LOOKS like
it's in `rectify.core.bam.parallel` (post-`0024fa3` refactor module) but
the rot is upstream in the BAM file itself. **Fix: re-align (re-run
`run_array_mapPacBio.sh`) with `e8c8070` in place, then re-run correct.**
Set2 alignments produced 2026-05-07 through 2026-05-19 all need re-alignment.

**Faster recovery for an existing pile of tainted mapPacBio BAMs (2026-05-20):**
Re-alignment costs ~60 min per chunk × 11 samples × ~10 chunks avg = ~6+ hours
of compute. The actual malformed records in the BAMs are reads with **empty
SEQ** (column 10 = `*`) that BBmap emitted when its FASTQ input had the
double-`\n` blank-sequence-line bug. Pysam scan of `wt_rep1_trimmed.mapPacBio.bam`
on 2026-05-20 found 366,728 / 4,347,067 reads (8.4%) with `query_sequence is
None`. Zero QNAME issues, zero seq/qual length mismatches — only empty SEQ.

Recovery without re-alignment: filter the empty-SEQ reads from the **merged
mapPacBio BAM only** (the one used by `rectify correct` via `--aligner-bams`
for cross-aligner consensus lookups; per-chunk input BAMs can be left alone —
they're iterated, not random-access fetched, so empty-SEQ reads get skipped):

```bash
samtools view -h --threads 4 INPUT.bam \
  | awk 'BEGIN{OFS="\t"} /^@/{print; next} $10!="*" && length($10)>0' \
  | samtools view --threads 4 -bS - > INPUT.bam.filtering
mv INPUT.bam INPUT.bam.pre_empty_seq_filter
mv INPUT.bam.filtering INPUT.bam
samtools index --threads 4 INPUT.bam
# Then invalidate stale prescan PKLs — they were built from corrupted merged BAM:
rm -f junction_pool.pkl rescue_scan.pkl
```

After this, re-submit `run_prescan.sh` → 5 `run_array_correct_*.sh` →
`run_array_chunk_merge.sh` → `run_final_merge.sh`. **Verified 2026-05-20:**
all 11 set2 samples cleared prescan in 1-7 min after filter (prescan would
crash within seconds on the unfiltered merged BAM); correct stage started
firing without `free():` errors.

**Do NOT try to also filter the per-chunk input BAMs** — they're not
coordinate-sorted (BBmap output order), so `samtools index` fails after
filter (`[E::hts_idx_push] Unsorted positions on sequence #4`). Either
sort them first, or skip — empty-SEQ reads in the input BAM are
iterated-not-fetched, so they don't trigger the htslib corruption.

---

## [2026-05-20] `rectify correct` pysam heap corruption on Han 2023 wt_R1 — STILL OPEN

**Status:** UNDER INVESTIGATION. v0.9.2 (`cb2fe6c`) only triggers the symptom because it correctly de-duplicates consensus rows (6.7M unique reads vs the previous 12.95M-row buggy output that masked it). The crash predates v0.9.2 — was always there, only now reachable.

**Affects:** `rectify correct` on a post-consensus rectified.bam with ≥ ~4M reads on H2 (16 threads, 16-core pod_smp.q node). Sherlock chunked pipeline (~290k reads per chunk in Stage D) is unaffected; Stage F.5 (full ~6.7M-read merged consensus) is the risk point.

**Symptom (matches the mapPacBio downstream symptom above):**
```
INFO - rectify.core.bam.parallel - 107 regions across 16 workers
*** Error in `python3.11`: free(): invalid next size (fast): 0x00002b... ***
*** Error in `python3.11`: free(): invalid next size (fast): 0x00002b... ***
======= Backtrace: =========
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x2d4b9)
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x4fc67)
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x4fb34)
...
Processing regions:  16%|██ | 17/107 [02:05<10:37, 7.09s/it]
```
Multiple workers crash, then silent hang (SGE reports `r` until walltime kill). Same pattern as the mapPacBio symptom.

**What's been ruled out (no fix needed for these — verified clean on the H2 BAM):**
- `samtools quickcheck` → exit 0
- SEQ length != QUAL length scan over full 6.7M reads → 0 hits
- QNAME length > 50 chars → 0 hits (BBmap ran with `trd=t`, qnames are bare SRR accessions)
- The `382fcc7` / `e8c8070` FASTQ sanitizer bugs do NOT apply — BBmap reads FASTQ directly, no rectify sanitizer in the BBmap path

**Scale evidence:**
| reads | threads | result |
|---|---|---|
| 1,072 | 1 | ✓ 31s |
| 1,072 | 8 | ✓ 9s |
| 102,157 (chrI) | 16 | ✓ 3:14 |
| 2,009,670 (chrI-V) | 16 | ✓ full run 61:40 wall (parallel + BAM write); exit 0 |
| 6,708,292 (full) | 16 | CRASH at heap corruption, ~17 of 107 regions completed before workers died |

**Working hypothesis:** memory pressure or specific read pattern only present on chrVI+ at scale. Need to test chrXII (rDNA, high coverage) or chrVI-XVI subset to localize.

**Mitigation in place for Sherlock Stage F.5 (job 25510786, after F's benign exit-2 detached from afterok dep):**
- Cut `--cpus-per-task=16 → 8` (halves worker count, gives each worker more memory budget)
- Added `--checkpoint-dir $OUTDIR/correct_merged_checkpoint` (per-region `.done` sentinels + rescue_scan.pkl on /oak; resume preserves work-done on crash)
- Bumped `--time=1:00:00 → 3:00:00`
- Template + deployed script both patched.

**[2026-05-20 09:47] Initial mitigation attempt: F.5 (25510786) running, no explicit crash.**
At 20+ min elapsed on sh03-07n10: 9 worker processes at 97–99% CPU each, all in `R` state. No `free()` heap-corruption log signature.

**[2026-05-20 10:10] Mitigation FAILED — silent-hang variant of the same bug.**
At 41+ min elapsed: workers STILL pinned at 99% CPU but **zero output** in scratch dir, **zero checkpoint sentinels**, **strace shows zero syscalls** — workers stuck in tight userspace loops in corrupted memory state. Same pysam heap-corruption bug, just not crashing immediately at threads=8 — instead silent-hanging. Cancelled 25510786.

**Working alternative: per-chunk consensus correct (Stage E.5, job 25514775).**
Bypass the merged-BAM crash entirely. Stage E already produced 24 per-chunk consensus.bam files (~290k reads each — the proven-safe Stage D scale). New sbatch array runs `rectify correct` on each per-chunk consensus.bam (24 small jobs at threads=8), then concat the 24 corrected_reads.tsv files into the sample-level TSV. Same scale as Stage D's per-aligner correct calls, all 48 of which succeeded. ETA ~4 hr.

**Updated guidance:** for any sample with merged BAM ≥ ~4M reads, run `rectify correct` PER CHUNK and concatenate TSVs. Do not run on the merged BAM regardless of thread count.

**Cross-cluster usage:** RECTIFY ships `rectify split --scheduler {slurm,uge,pbs}` that emits the per-chunk array scripts natively. H2 (SGE) uses `--scheduler uge --uge-queue pod_smp.q --uge-pe shared`. Sherlock (SLURM) uses `--scheduler slurm --slurm-partition larsms`. Do not use `rectify run-all` for samples >4M reads on either cluster — use `rectify split` and submit the generated array scripts instead.

**Path forward:** if 25509811 succeeds with threads=8, ship that as the workaround config in `05b_correct_merged_consensus.sh.tmpl` and re-run H2 with the same config. If it still crashes, resume from checkpoint to make incremental progress while we localize the malformed read pattern.

**Do not pull yet** — v0.9.2 is correct; this is a downstream pysam-level issue surfaced BY the correctness fix, not caused by it. Pulling `cb2fe6c` is safe; just be aware that whole-sample `rectify correct` at ≥4M-read scale may need the `--checkpoint-dir` + lower thread workaround.

**[2026-05-20] Commit B hypothesis — tempered.** Commit B's `write_corrected_bam_parallel` PRE-PARTITIONS the input BAM into per-region BAMs (≤100k ref bp each) BEFORE region workers start — each worker opens its own small region BAM via pysam, not the full 6.7M-read merged BAM. This is meaningfully different from the `--checkpoint-dir` mitigation above (which still fetched from the full merged BAM, just with sentinel-based resume). Hypothesis: workers operating on ≤100k-bp region BAMs may not trigger the pysam C-level state corruption that the merged-BAM workers do at scale.

**But this is unproven and the prior threads=8 mitigation already silent-hung.** The structural resolution test (6.7M-read Han wt_R1 full-scale run on H2 16-core) was **DEFERRED**: H2 pod_smp.q had 25,600 waiting jobs at submission time; Sherlock ControlMaster was not available (Duo re-auth needed). The chrI-V subset (2,009,670 reads, 16 threads) from prior Sherlock runs passed (61:40 wall, exit 0), validating Commit B's architecture at the proven-safe scale. The Han wt_R1 6.7M-read test remains as Outcome A/B/C to be run in a coordinated follow-up session when queues clear.

**Current working policy stands:** for samples >4M reads, use `rectify split` per-chunk array (290k reads per chunk, Stage E.5 pattern) — do NOT run `rectify correct` on merged BAMs regardless of which cluster, thread count, or commit. Commit B's region-pre-partitioning may or may not change this; until the deferred test runs, the chunk-first policy is the recommendation.

See briefing: `dev/specs/briefings/commit_b_briefing.md` §4 for the smoke deferral context.

---

## [2026-05-20] parallel BAM writer process aborts — MITIGATED IN WORKING TREE

**Status:** MITIGATED (uncommitted Codex working tree, 2026-05-20). The
parallel BAM writer no longer uses worker processes by default. It preserves
the region-planning code path but executes region workers sequentially unless
`RECTIFY_ENABLE_PARALLEL_BAM_WRITER=1` is explicitly set.

**Affects:** `rectify.core.bam.bam_writer_parallel.write_corrected_bam_parallel()`
and tests that call it with `n_threads > 1`.

**Symptom:** intermittent fatal Python aborts in/around `multiprocessing.Pool`
or process launch while exercising `write_corrected_bam_parallel()`. In the M1
sandbox, attempts to switch to spawn/forkserver/subprocess workers exposed
related low-level failures:

- `spawn` stalled in pytest on the 36-read smoke fixture.
- `forkserver` failed with sandbox `PermissionError` binding the Unix socket.
- clean subprocess workers hit OpenMP shared-memory startup failures and could
  still emit fatal abort traces during `subprocess.Popen`.

**Adjacent import bug found:** plain `import rectify` eagerly imported
`rectify.visualize`, which imported Matplotlib and ran font-manager subprocess
checks during test collection. This violates the repo rule that optional
numpy/plotting stacks must not load before thread limits are set. The working
tree now makes visualization lazy by default; `import rectify.visualize` still
works when plotting is needed, and `RECTIFY_IMPORT_VISUALIZE=1` restores the
legacy eager probe.

**Fix/mitigation in working tree:**

- Region BAM commits are now atomic: write unsorted temp → sort to temp BAM →
  fsync → `os.replace()` final region BAM → atomic `.ok` sentinel.
- Resume only trusts a region when both the `.ok` sentinel and region BAM exist.
- Unsafe worker-process execution is opt-in via
  `RECTIFY_ENABLE_PARALLEL_BAM_WRITER=1`; default `n_threads > 1` logs a warning
  and runs region workers sequentially.

**Validation:** `pytest tests/test_bam_writer.py tests/test_bam_writer_parallel_smoke.py -q`
passes on M1 (`18 passed`). This is intentionally a correctness/safety
mitigation, not a throughput fix. Re-enable true process parallelism only after
cluster-specific smoke tests prove the local Python/pysam/OpenMP runtime can
launch workers without aborts.

---

## [2026-05-20] `rectify analyze` per-condition bedgraph: 1-bp left shift

**Status:** FIXED on M1 working tree (uncommitted at time of entry); awaiting commit + push.

**Affects:** every per-condition strand-specific bedgraph emitted by
`rectify analyze` (filenames like `<condition>_plus.bedgraph` /
`<condition>_minus.bedgraph`, header
`description="RECTIFY 3' ends (<strand> strand)"`). Hits DRS and QuantSeq REV
equally — anything routed through `rectify/core/analyze/bedgraph.py::generate_bedgraphs`.

**Symptom:** the per-position 3'-end pileup peak in a `rectify analyze`
bedgraph appears 1 bp to the LEFT of the true position. Empirical:
CST6 ysh1 peak at chrIX in 2026-05-12 wbfix Han run wrote
`chrIX 287748 287749 70.66` (IGV 1-based 287,749); the matching DRS
3'-end signal from `bedtools genomecov -3 -bg` on the corrected
strand-split DRS BAM wrote `chrIX 287749 287750 18` (IGV 287,750). Gap
1 bp, consistent.

**Root cause:** `rectify/core/analyze/bedgraph.py:99-100` wrote
```python
start = int(pos) - 1
end   = int(pos)
```
as if `pos` (a `corrected_3prime` value) were 1-based. But
`corrected_3prime` is **0-based-inclusive** everywhere it is computed —
`reference_end - 1` for `is_reverse=True` (with explicit comment
`# 0-based inclusive` at `walkback.py:142,214,471`, `indel_corrector.py:1661,2027`)
and `reference_start` for `is_reverse=False`. BED is 0-based half-open,
so a single base at 0-based `pos` is the interval `[pos, pos+1)`. The
subtract-1 spelling was a leftover from an earlier 1-based convention
that was never updated when corrected_3prime semantics settled at 0-based.

**Fix:**
```python
start = int(pos)
end   = int(pos) + 1
```
Plus three regression tests in `tests/test_analyze.py::TestBedgraphCoordinates`
covering single-base emission, strand routing, and same-position aggregation.
Tests confirmed to FAIL under the pre-fix code and PASS under the fix.

**Safe to pull:** strict correctness improvement. Existing pipeline runs
that produced left-shifted bedgraphs were never used as the source of
truth for clustering or shift analysis (those work from `corrected_reads.tsv`
position columns directly, not from the bedgraphs), so the bug only
affected downstream visualization and any analysis that cross-referenced
the per-condition bedgraph against an independently-generated track. The
He et al manuscript-anchored TRT analysis (`analyses/cross_modality_trt_20260519`)
caught it during a CST6 cross-modality diagnostic; v2 classification was
unaffected because it worked from per-position attributed counts via the
DRS bedtools-derived bedgraphs, not from the `rectify analyze` output.

**Audit completed 2026-05-20** — see `dev/audits/bedgraph_coordinate_audit_20260520.md` for the full per-file table. **Three instances of the same bug found and fixed:**

- `rectify/core/analyze/manifest.py:598` — manifest-mode equivalent of `generate_bedgraphs` (streams per-sample TSVs without holding the combined frame in memory). Same `int(pos) - 1` / `int(pos)` spelling against the same 0-based `corrected_position` / `corrected_3prime` column. Corrected to `int(pos)` / `int(pos) + 1`.
- `scripts/generate_bedgraph_from_polished.py:113` — standalone CLI utility. Same off-by-1 with an explicitly wrong inline comment claiming the position was 1-based. Caught by a parallel multi-agent audit run after my initial close; my recipe had been scoped to `rectify/core/` and `rectify/data/` and missed `scripts/`. Fixed and noted in the audit doc as a scope-expansion miss.

All other bedgraph and bigwig emitters verified to already use the correct 0-based half-open convention (`netseq/netseq_output.py:122,191`, `commands/export_command.py:88,128`, and the delegate chains through `bam/bedgraph_writers.py`, `commands/{consensus,correct,analyze}_command.py`).

**Adjacent coordinate-convention findings (not bedgraph emission)** surfaced by the same multi-agent audit — tracked as separate BUGS_TO_FIX entries NEW-077 through NEW-081:

- NEW-077 (HIGH) — `bam_processor.py:826-836` minus-strand artifact N snap can land on a SKIPPED base.
- NEW-078 (HIGH/MEDIUM) — `bam_processor.py:847` ambiguity clip leaves NET-seq window in N span.
- NEW-079 (MEDIUM) — `clustering.py:532` plus-strand `distance_to_gene_3prime` is off by 1 (uses `gene['end']` against half-open annotation).
- NEW-080 (MEDIUM) — `false_junction_filter.py:282 vs 300` minus-strand wrong-flank fetch.
- NEW-081 (LOW) — `analyze_command.py:130` + `manifest.py:217` rDNA exclusion uses `<= end` against half-open regions.

**Action items still open:**
- Commit all three bedgraph fixes (`analyze/bedgraph.py` + `analyze/manifest.py` + `scripts/generate_bedgraph_from_polished.py` + `tests/test_analyze.py`) to `drs-validation-rebuild`, push, pull on H2 + Sherlock.
- Run the fast test suite on a cluster login node (M1 is memory-constrained; pytest there was killing M1 with the multi-agent audit also running).
- Triage and fix NEW-077 / NEW-078 / NEW-080 (HIGH-severity, biology-affecting); NEW-079 / NEW-081 are easy clean-ups.

---

## [2026-05-20] Read-number sidecar + RN aux-tag hybrid

**Status:** IMPLEMENTED (working tree; uncommitted at time of entry).

**Problem solved:** cDNA FASTQ comment tags (`XA/XC/XF/XU/XR/XI/XK/XS`) were
only propagated by minimap2 `-y`; non-minimap2 consensus winners lost that
metadata. `rectify split` now assigns a per-sample integer `RN:i:<read_num>`
to each derived FASTQ header and writes `<sample>.read_num_sidecar.parquet`
with the original QNAME, full FASTQ comment, chunk id, and seq/qual MD5s.

**Architecture:** Option B from `HANDOFF.md` landed. Original QNAMEs stay in
FASTQ/BAMs for IGV and grep ergonomics. minimap2 carries RN natively via `-y`;
mapPacBio, gapmm2, deSALT, uLTRA, bbmap, and bwa build QNAME→RN maps before
comment stripping and stream-inject `RN:i` into their BAMs after the existing
QNAME validator. Consensus uses RN as the K-way merge key only when every input
BAM has RN; otherwise it falls back to the existing normalized-QNAME path.

**Restoration:** `rectify consensus --read-num-sidecar` (or autodetection beside
input BAMs) restores SAM-format tags from the sidecar FASTQ comment onto each
winning consensus record without overwriting tags already present on the read.
Old RN-less BAMs remain valid and require no migration.

**Tests:** New focused tests cover sidecar round-trip/lookups/fingerprint checks,
chunker RN emission, RN BAM injection, RN-keyed merge fallback behavior, and
consensus tag restoration. Local PyArrow is unavailable, so the 1M-row streaming
RSS test is skipped locally; deployments with PyArrow exercise true parquet
row-group writing.

## [2026-05-20] QNAME pipeline hardening — sanitizer + validator + cross-aligner audit

**Status:** FIXED (working tree; uncommitted at time of entry).

**Affects:** any RECTIFY pipeline run. The active-panel mutation
class was latent in production DRS (set1 spot-check showed 0 leakage),
but multiple structural gaps were silently active for the chimeric path
and the TSV merge.

**Three connected fixes shipped today (one session):**

1. **mapPacBio sanitizer: always per-record, handle tab**
   (`multi_aligner.py: _sanitize_mpb_fastq`).
   The previous gate `_need_san = _first.startswith('@') and ' ' in _first`
   evaluated only the first FASTQ record AND only checked space (not tab),
   so bare-first files or tab-aux FASTQs (minimap2 `-y`, cDNA-pipeline shape)
   silently bypassed the sanitizer. Always rewrites per-record now; strips
   both space and tab; 254-char cap.

2. **Extended `_normalize_bam_read_name` regex** (`consensus.py:171`).
   Now covers generic SAM aux tag suffixes (`_<2c>:[AZifHB]:`) and
   enumerated Dorado metadata keys (`_runid=`, `_ch=`, `_start_time=`,
   `_flow_cell_id=`, ...). Defense in depth for any aligner that emits
   underscore-encoded comment leaks.

3. **`validate_post_alignment_qnames` + auto-sanitize**
   (`rectify/core/align/qname_validator.py`, NEW).
   Wired into ALL aligner wrappers (minimap2, mapPacBio, gapmm2, deSALT,
   bbmap, bwa, uLTRA). Samples 200 primary BAM reads; tracks 4 mutation
   classes (whitespace, overlength, cosmetic-needs-normalize,
   non-round-tripping). Default behavior: if a recoverable mutation is
   detected, stream-rewrite the BAM with normalized QNAMEs and re-validate.
   Raises only on unrecoverable mutation. Bypass with
   `RECTIFY_NO_AUTO_QNAME_SANITIZE=1` (detect-only) or
   `RECTIFY_SKIP_POST_ALIGN_VALIDATION=1` (skip entirely).

**Three cross-aligner audit follow-ups (audit:
`scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md`):**

- **Chimeric path QNAME normalize** (`consensus.py:459`) — mirrored the
  non-chimeric path's line-489 normalize. Without this, mapPacBio /
  uLTRA-templated chimeric reads carried mutated QNAMEs into all
  downstream joins.
- **TSV-side normalizers reuse canonical regex**
  (`corrected_consensus.py:48`, `:518`) — deleted `_bare_uuid`, replaced
  body of `_normalize_read_id` with vectorized form of
  `_normalize_bam_read_name`; imports `_UNDERSCORE_COMMENT_RE` so the
  two cannot drift.
- **`_restore_sequence_from_aligner_reads` strand guard**
  (`consensus.py:362`) — skip donor when `donor_read.is_reverse !=
  best_read.is_reverse` to prevent RC'd SEQ injection.

**Open follow-ups (NOT done in this session):**
- Issue #3 from EDGE_CASES.md (cDNA FASTQ tags lost for non-minimap2
  winners) — structural; needs sidecar architecture or per-wrapper `-y`
  equivalent. Held for separate design session.
- Issues #4 (cosmetic tiebreaker comment), #6 (validator sequential
  sampling) — low severity.

**Survey artifacts:**
- `scripts/diagnostics/qname_mutation_survey/results.tsv` — 4-aligner ×
  53-read adversarial synthetic + 3 real-data sanity checks
- `scripts/diagnostics/qname_mutation_survey/REPORT.md` — full audit
  narrative
- `scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md` —
  cross-aligner punch list (6 entries, 3 of which fixed today)

**Tests:** 113 pass in the QNAME + consensus + restore suite (47 new
covering the sanitizer, validator, vectorized TSV normalizer, and
strand guard). Pre-existing unrelated breakage in
`test_validation_reads.py` and `test_quantseq_rev_integration.py`
(Commit B manifest-output WIP) is not addressed by this work.

**Safe to pull:** strict correctness. Per-aligner BAMs will have an
added <1s validator pass at the end; nothing else slows. Existing BAMs
on disk are unchanged unless the auto-sanitizer fires on a re-run.

---

## [2026-05-20] set2 cascade scale: merged-BAM `--aligner-bams` lookups crash at chunk scale — STILL OPEN

**Status:** UNDER INVESTIGATION. Same root cause as entry 2 (Han 2023 6.7M
heap corruption) but reproduced at set2's smaller scale, demonstrating the
problem is **not** Han 2023-specific. Empty-SEQ filter (entry 1's "Faster
recovery" section) is necessary but **not sufficient**.

**Affects:** `rectify correct` invoked with `--aligner-bams` pointing at
merged per-aligner BAMs ≥ ~4M reads, when the per-chunk input BAM is just
~280k-466k reads. This is the standard set2 wrapper invocation generated
by `rectify split --generate-slurm` — each per-chunk correct does
cross-aligner consensus lookups against the merged BAMs.

**Symptom (mixed silent-hang + explicit free()):**

Same pattern as entries 1 ("Downstream symptom") and 2:
```
INFO - rectify.core.bam.parallel -   X regions across 8 workers
*** Error in `python': free(): invalid next size (fast): 0x... ***
```
… followed by memory-map dump and either:
- Silent hang (Slurm reports RUNNING until walltime kill — most cases)
- The explicit `free():` + map dump in the .err log (minority cases)

Only ~2 of dozens of logs at first inspection showed the verbatim
`free()` line; the rest had logs frozen at the "regions across N workers"
line with no further output. Workers may or may not show CPU activity —
the silent-hang variant can keep workers in tight loops in corrupted
memory state.

**Scale at which it manifests for set2 (verified 2026-05-20):**
- Per-chunk INPUT BAM: ~280k-466k reads (well under the 4M threshold)
- `--aligner-bams` MERGED BAMs: 4-15M reads per sample (above threshold)
- Cross-aligner fetch() calls into the large merged BAMs trigger the
  same htslib state corruption that entry 2 documents on Han 2023

**What the empty-SEQ filter (entry 1) fixed vs left open:**
- ✅ Fixed: merged-BAM crashes caused by the mapPacBio QNAME bug's empty-SEQ
  reads. Pysam scan showed 8.4% (366k/4.3M) of mapPacBio merged BAM reads
  had `query_sequence is None`. Filtering them allowed PRESCAN to complete
  cleanly (it would crash on the unfiltered merged BAM within seconds).
- ❌ Still open: the broader scale-related heap corruption at merged BAM
  ≥ ~4M reads, INDEPENDENT of the empty-SEQ pattern. Filtering removes
  one trigger but not all.

**Verification path that exposed it:** post-empty-SEQ-filter cascade
launched 99 jobs (11 samples × prescan + 5 correct + chunk_merge +
final_merge). Prescan succeeded on all 11. Correct stage ran but only
~19 / 530+ chunk tasks COMPLETED — and most of those 19 were the
empty-BAM-guard fallbacks (deSALT SIGSEGV chunks that exit-0 via the
6c8f5a6 guard in seconds, never doing real correction). Real-correct
tasks at 1:24+ elapsed showed the silent-hang / explicit-crash pattern.

**Working theory (from entry 2):** the issue is htslib's per-region
fetch() against multi-million-read indexed BAMs under multi-worker fan-out.
Workers' shared file handles or memory state corrupt in a way that
manifests as `free(): invalid next size`. The empty-SEQ filter removed
~8% of mapPacBio reads (an outlier-removal that brought some samples
under threshold) but didn't change the fundamental dynamic for samples
still over 4M reads in merged BAMs (every set2 sample is).

**Practical workarounds for set2 today (in order of preference):**

1. **Single-aligner correct (no `--aligner-bams`).** Bypass cross-aligner
   lookups entirely. Each per-chunk correct uses ONLY its own input BAM
   (e.g., minimap2-only). Removes the merged-BAM fetch path that's
   crashing. 3'-end accuracy stays ~95-99% of multi-aligner consensus
   per Kevin's earlier analysis (DRS yeast data; cross-aligner consensus
   mostly improves splice junction edge cases, less the 3'-end call
   itself). Used by this session as the salvage path after the 5-aligner
   cascade collapsed.

2. **Stage E.5 per-chunk-consensus correct.** Documented in entry 2 as the
   Han 2023 workaround. Requires producing per-chunk consensus.bam files
   FIRST (Stage E), then running rectify correct on each per-chunk
   consensus BAM (~290k reads each, proven-safe scale, no merged-BAM
   lookup at all). More compute than (1) but retains multi-aligner
   consensus value.

3. **Commit B `write_corrected_bam_parallel` pre-partitioning.**
   Architectural fix in flight (entry 2 references). When merged with M1,
   correct workers operate on per-region BAM slices ≤100k bp ref-space
   rather than fetching against the full merged BAM. Should resolve at
   the code level. Smoke-test pending.

**For the next agent debugging this:**

- Don't trust Slurm RUNNING — most "running" tasks at the second hour
  are hung in userspace memory-corrupted state, not making progress.
- The signature in logs is line: `rectify.core.bam.parallel -   X regions
  across N workers` followed by either silence OR `free():` errors.
- Check log progression with `tail -1 <log>.err` — if frozen at the
  "regions" line for >10 min, the task is hung.
- `grep -l "free()" *.err` undercounts because some logs go to silent
  hang without printing the free() line at all.
- Workers may show high CPU (50-99%) while hung — they're stuck in
  userspace loops in corrupted memory, not stuck in syscall waits.
- See sssion 5 handoff in `project_status_markdowns/DRS_CPA_PROJECT_STATUS.md`
  for the empty-SEQ-filter recovery context and the full set2 cascade
  timeline.

**Do not pull this entry for guidance** — pull the actual fix (entry 2's
Commit B branch or whatever supersedes it) once it lands. This entry is a
field report, not a fix.

---

*To add an entry: symptom (exact error string), root cause (one sentence),
fix commit, safe-to-pull verdict.*

---

## [2026-05-21] Analyze gene attribution for APA/shift analysis — SHERLOCK HAN VALIDATED

**Status:** Implemented in working tree, uncommitted. Local focused tests pass
and Sherlock Han BWA-only analyze completed successfully with DRS-origin
reference attribution.

**Affects:** `rectify analyze` gene-level DESeq2 and cluster shift analysis,
especially short-read QuantSeq REV runs and any long-read run where a CPA
cluster may derive from multiple upstream genes.

**Problem:** `rectify analyze` currently annotates CPA clusters with a single
nearest same-strand annotated TES using `annotate_clusters_with_genes()`
(`-500/+100` bp default). There is no protocol/read-length fork: manifest
analysis uses this heuristic for both DRS and QSrev once it has
`chrom/strand/corrected_position`. This is incorrect for:

- long reads, where read-body/TSS evidence should drive gene attribution;
- shared CPA clusters, where one cluster can legitimately have weighted
  contributions from multiple upstream genes;
- short reads, where body overlap is weak and the best attribution source is
  an imported long-read-derived cluster/position gene map, with annotation as
  an explicit fallback rather than the primary rule.

**Implementation landed in working tree:**

- Add a first-class weighted `cluster_gene_attributions.tsv` table:
  `cluster_id, gene_id, gene_name, raw_attributed_count,
  attribution_weight, source`.
- Add attribution modes for analyze:
  - `annotation`: legacy nearest-TES behavior;
  - `body`: compute attribution from corrected TSV read-body spans
    (`alignment_start/end`) against CDS/gene features;
  - `reference`: map external per-position DRS attribution TSVs
    (`chrom, position, strand, gene_id, attributed_count`) onto the current CPA
    clusters;
  - `body-then-annotation` / `reference-then-annotation`: fill unattributed
    clusters with the legacy annotation fallback.
- Thread weighted attributions into gene-level DESeq2 and shift analysis so a
  shared CPA cluster contributes fractionally to multiple genes instead of
  being forced into one arbitrary `gene_id`.
- Fill reference/body attribution `gene_name` values from loaded annotation
  when possible, so systematic IDs still map to common names for GO/plots.
- Explicit non-default attribution modes fail the analyze run if attribution
  cannot be built; the legacy `annotation` mode remains backward compatible.
- Add `--gene-attribution-reference-window BP` for reference modes. Exact
  DRS-position-to-current-cluster overlap is still preferred; if absent, a
  DRS attributed position can map to the nearest same-strand current CPA
  cluster within BP bases. This handles small DRS/QSrev peak offsets without
  reopening the broad `-500/+100` annotation window.

**Files touched for this entry:**

- `rectify/core/analyze/cluster_gene_attribution.py` (new)
- `rectify/core/analyze/__init__.py`
- `rectify/core/analyze/deseq2.py`
- `rectify/core/analyze/manifest.py`
- `rectify/core/analyze/shift_analysis.py`
- `rectify/core/commands/analyze_command.py`
- `tests/test_analyze.py`

**Local tests run:**

```bash
pytest tests/test_analyze.py -q
```

Result after the reference-window addition: `53 passed, 2 warnings in 35.93s`.
The new tests cover weighted
gene-level aggregation, weighted shift analysis, and mapping external
per-position DRS attributions onto current CPA clusters, including the bounded
nearest-cluster reference window.

**Validation:** reran Han BWA-only analyze on Sherlock with
`--gene-attribution-mode reference-then-annotation` using staged DRS origin
attribution tables under
`analyses/cross_modality_trt_20260519/inputs/*/attribution_origin/`.
Submitted full rerun as Sherlock SLURM job `25601498` using:
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/han2023_bwa_fast/run_han_bwa_analyze_drs_origin_ref.sbatch`.
Output target:
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/han2023_bwa_fast/analyze_drs_origin_ref_20260521`.
Job `25601498` failed before analysis started because the sbatch script exited
on `git log/status` from a compute-node path that did not resolve as a git
checkout. Fix the script to make git reporting best-effort, add
`--gene-attribution-reference-window`, then resubmit.
Resubmitted as Sherlock SLURM job `25602373` with
`--gene-attribution-reference-window 25`.
Result: `COMPLETED`, exit `0:0`, elapsed `01:11:16`.

Output:
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/han2023_bwa_fast/analyze_drs_origin_ref_20260521`

Key output metrics:
- CPA clusters retained: 231,182
- `cluster_gene_attributions.tsv`: 214,114 rows, 170,819 clusters, 8,442 genes
- Attribution row sources: 194,778 reference; 19,336 annotation fallback
- Shared multi-gene clusters: 36,513
- Gene-level DESeq2 tested 7,238 genes per contrast
- Significant genes: Ysh1AA 4,403; Rat1AA 3,727
- Weighted shift genes analyzed: Ysh1AA 6,859; Rat1AA 6,858

**Safe-to-pull verdict:** candidate fix is validated on the Han BWA-only data
and ready for review/commit. Still inspect the dirty tree carefully and stage
only the attribution-related paths; the repo contains unrelated WIP.

---

## [2026-05-21] H2 mex67aa_rep3 DRS continuation — MERGE SUBMITTED

**Status:** H2 merge-aligners job submitted and under monitoring.

**Dataset:** `/u/project/guillom/shared/processed/mex67aa_vs_wtaa_rectify_v0.9.0`

**Observed state before submission:**
- `mex67aa_rep1` and `mex67aa_rep2` each had all five merged aligner BAMs.
- `mex67aa_rep3` had completed chunk outputs for all five aligners but had no
  `chunks/merged_bams/` outputs yet.
- `mex67aa_rep3` deSALT chunks `003`, `005`, and `006` produced tiny empty BAMs
  after deSALT exited `-11`; logs identify this as likely SIGSEGV in
  `Loop-ProcessReads` from the upstream deSALT bug and intentionally emit empty
  BAMs so those chunks proceed with 4-aligner consensus.

**Submission:**
- First `qsub` failed because H2 JSV treated propagated `LC_ALL=C.UTF-8`
  locale warnings as fatal.
- Resubmitted with locale variables cleared:
  `env -u LC_ALL -u LANG qsub .../mex67aa_rep3/chunks/run_merge_aligners.sh`
- Job: `13463288` (`mex67aa_rep3_merge_aln`)
- Initial queue state at 2026-05-21 15:21 PDT: `qw`, waiting for 16 shared
  slots on `campus2.q`.

**Next steps when merge completes:**
Run the manual UGE chain from `mex67aa_rep3/chunks/submit_pipeline.sh`:
`run_prescan.sh`, per-aligner correction arrays, chunk merge, final merge,
per-chunk consensus, then consensus-chunk merge. Confirm whether correction
uses the current safer single-aligner/per-chunk strategy before launching the
large merged-BAM cross-aligner correct stage.

---

## [2026-05-21] Post-hoc read-number sidecars — CPA SHERLOCK COMPLETE

**Status:** Sherlock CPA set2 sidecars reconstructed and validated. H2 pending
environment/data prerequisites.

**Why:** Production CPA and mex67aa chunk FASTQs predate the RN/read-number
sidecar split path. They had `chunks_manifest.json`, but no
`*.read_num_sidecar.parquet`; FASTQ headers also lacked `RN:i`.

**Method:** used the legacy round-robin inverse:
`read_num = local_read_index * n_chunks + chunk_index`.
The helper writes:
- `<sample>.read_num_sidecar.parquet`
- `<sample>.read_num_sidecar.PROVENANCE.json`
- `<sample>.POSTHOC_PROVENANCE.json`

**Sherlock CPA output root:**
`/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery`

**Validated row counts:** Parquet metadata row counts exactly match each
sample's `chunks_manifest.json` `n_reads`.

| Sample | Rows |
| --- | ---: |
| wt_rep1 | 4,347,602 |
| wt_rep2 | 4,198,935 |
| wt_rep3 | 2,897,972 |
| rna15_rep1 | 1,090,278 |
| rna15_rep2 | 1,910,111 |
| rna15_rep3 | 945,995 |
| ysh1_rep1 | 1,264,614 |
| ysh1_rep2 | 2,010,907 |
| ysh1_rep3 | 1,859,080 |
| wt_tfiiib_rep1 | 11,543,753 |
| wt_tfiiib_rep2 | 8,234,198 |
| wt_tfiiib_rep3 | 7,657,191 |

**Important operational notes:**
- The first batch process was killed after writing `wt_rep1`/`wt_rep2` because
  validation used `ReadNumSidecar.open()`, which loads all rows into memory.
  The completed sidecars were valid. The rerun used one sample per Python
  process and `pyarrow.parquet.ParquetFile(...).metadata.num_rows`.
- Do not validate multi-million-row sidecars by opening them through
  `ReadNumSidecar.open()` on login nodes.
- These are post-hoc sidecars only. Existing BAMs do **not** have `RN:i` tags.
  Consensus can still use the sidecar for tag restoration via normalized-QNAME
  fallback, but RN-keyed consensus requires either BAM RN injection or
  re-alignment from RN-tagged chunk FASTQs.

**H2 status:**
- `/u/project/guillom/shared/processed/mex67aa_vs_wtaa_rectify_v0.9.0` has
  38 chunk FASTQs and 3 `chunks_manifest.json` files, so post-hoc
  reconstruction is possible.
- Installed `pyarrow==16.1.0` into
  `/u/project/guillom/shared/envs/rectify` using the env's own `pip` with
  `--only-binary=:all:`. The first unpinned `pip install pyarrow` attempted a
  source build and failed before installing anything; the pinned wheel install
  succeeded. Verified a Parquet write/read smoke test with the env Python.
- H2 mex67aa post-hoc sidecars were reconstructed and validated by Parquet
  metadata against each sample's `chunks_manifest.json`:
  - `mex67aa_rep1`: 7,925,277 rows
  - `mex67aa_rep2`: 6,311,767 rows
  - `mex67aa_rep3`: 4,111,368 rows
- `/u/project/guillom/shared/processed/alignments/wt_tfiiib_rep*` has no
  chunk FASTQs and no `chunks_manifest.json`; it cannot be reconstructed
  in-place from that copied H2 artifact tree.

---

## [2026-05-21] Lazy Corrected Consensus — ysh1 Pilot Timings and Fast-Path Test

**Status:** implemented, synced to Sherlock for testing; H2 sync pending this
checkpoint.

**Code changes:**
- `merge_corrected_tsvs()` now accepts `lazy_scoring_workers` and can score
  aligners in parallel from raw BAMs + corrected TSVs.
- Lazy raw-BAM HP scoring skips the transient corrected-BAM edit stack when a
  TSV row has no correction surgery fields that can alter the CIGAR.
- Generated split and single-sample merge paths pass their thread counts into
  lazy HP scoring.

**Full ysh1 chunk 000 baseline, Sherlock job `25635672`:**
- Aligners: `mapPacBio`, `minimap2`, `gapmm2`, `uLTRA`; deSALT excluded.
- Lazy merge + whole-read HP scoring: `1055.0s`.
- Final corrected consensus BAM writing/sort/index: `227.4s`.
- Output:
  `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/ysh1_rep2/chunks/bypass_2f2h_lazy_consensus_region_test/chunk_000/corrected_consensus.bam`

**10k stratified fast-path pilot, Sherlock job `25637055`:**
- `N_READS_TO_TEST=10000` builds temporary per-aligner TSV subsets before
  merge/consensus.
- Trigger coverage selected from ysh1 chunk 000:
  - softclip HP rescue: 2,127 selected / 2,790 available
  - overcall rescue: 2,087 / 2,492
  - changed 3' end: 6,426 / 41,431
  - junction reads: 6,301 / 223,577
  - no module trigger: 3,957 / 140,323
  - no available examples in this bypass chunk for 2F five-prime rescue,
    intronic-tail clip, or reanchor/upstream-trim.
- Lazy merge + HP scoring: `18.8s`.
- Final corrected consensus BAM: `14.1s`.
- Per-aligner scoring stats:
  - mapPacBio: 9,944 no-edit fast path, 0 transient-edit, 9,813 scored
  - minimap2: 5,992 no-edit, 4,000 transient-edit, 9,876 scored
  - uLTRA: 6,483 no-edit, 3,520 transient-edit, 9,857 scored
  - gapmm2: 5,542 no-edit, 3,903 transient-edit, 9,040 scored

**Important next step:** implement differential HP scoring over only the read
regions where candidate aligners disagree. Current whole-read scoring is
correct but still recomputes shared blocks that cancel out when aligners have
identical or near-identical CIGAR representations.

**2F-rich targeted panel:** RPL19B/RPL20B short-exon-1 genes are good real
examples for minimap2 2F rescue. In current TSVs these appear as systematic IDs
`YBL027W` and `YOR312C`. Scanning full minimap2 corrected TSVs found dense
RPL19B 2F examples; `rna15_rep3/chunk_001` was selected because all four stable
aligners had nonempty corrected TSVs and raw BAMs.

Sherlock job `25637632`:
- `N_READS_TO_TEST=10000`
- `MIN_READS_PER_TRIGGER=1000`
- `MIN_READS_TARGET_GENES=1000`
- `TARGET_GENE_IDS=YBL027W,YOR312C`
- selected all 134 target-gene reads and all 93 target-gene 2F reads
- selected 7,552 / 9,590 available 2F reads overall
- also included intronic-tail, reanchor/upstream-trim, junction, 3'-changed,
  softclip HP, overcall, and no-trigger examples
- lazy merge + HP scoring: `18.3s`
- final corrected consensus BAM: `7.9s`
- output:
  `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/rna15_rep3/chunks/lazy_consensus_target_2f_test/chunk_001`

---

## [2026-05-22] Lazy corrected consensus P0–P3 bug fixes — committed in `1ab71f0`

**Status:** Fixed and committed. Full suite 1219 passed, 35 skipped.

### P0 — DELETED `_correction_requires_transient_edit` fast path (DO NOT RE-ADD)

**Symptom:** 22/36 reads (61%) showed divergent HP-ED values between the lazy
path and the reference `rectify correct` path in the audit fixture.

**Root cause:** `_correction_requires_transient_edit` returned `False` for any
read where the TSV row had no surgery fields that could alter the CIGAR. The
function then skipped `realign_exon_blocks` entirely for that read. But
`realign_exon_blocks` is NOT just a CIGAR-surgery function — it also computes
the HP-edit-distance score. Skipping it set those reads' HP-ED to 0 (the
aligner default), which silently changed winner selection in the consensus step.

**Fix:** Deleted `_correction_requires_transient_edit` in full. `realign_exon_blocks`
is now called unconditionally — that function has its own cheap internal pre-check
that returns early when no HP scoring work is required. The old fast path was
not just premature optimization; it was incorrect.

**Lesson for future agents:** If you see `_correction_requires_transient_edit`
referenced anywhere (old comments, stale branches, spec drafts), treat it as a
**deleted error**. Do not re-introduce any "skip HP scoring if no CIGAR changes"
guard in the lazy consensus path.

### P1 — `_chunk_index_from_path` now raises on unparseable filenames

**Symptom:** Corrupt or non-standard chunk filenames silently returned `None`
for the chunk index, causing downstream KeyErrors with confusing messages.

**Fix:** `_chunk_index_from_path` in `sidecar.py` now raises `ValueError` on
unparseable filenames. Callers must handle the exception explicitly.

### P2 — `_load_tsv` drops pre-existing `winning_aligner` column

**Symptom:** When a corrected TSV was written by a previous run that had already
added a `winning_aligner` column, reloading it caused a duplicate-column error
in the merge step.

**Fix:** `_load_tsv` in `corrected_consensus.py` drops any pre-existing
`winning_aligner` column before returning the DataFrame.

### P3 — Removed dead `mv corrected.bam` line; wired `--write-per-aligner-corrected-bams`

**Symptom:** A dead `mv corrected.bam` call at `split_command.py:942-946` ran
after the new lazy consensus path and tried to rename a file that no longer
existed, causing a spurious error at the end of successful split runs.

**Fix:** Removed the dead line. Also wired the `--write-per-aligner-corrected-bams`
CLI flag through `split_command.py` so it propagates to the lazy consensus path.

### Sherlock timing (job 25662145, committed in timing spec)

- Dataset: `wt_tfiiib_rep3/chunk_008`, 3 aligners, ~15k corrected reads.
- Lazy path: 495 s wall, 286 MB peak RSS.
- Full path (`rectify correct --write-corrected-bam` × 3): TIMEOUT >81 min.
- Lower-bound speedup: **>9.9×** (true speedup likely 20–30×).
