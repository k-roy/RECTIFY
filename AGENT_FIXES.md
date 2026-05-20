# AGENT_FIXES.md

Fast coordination log for active debugging sessions across M1 / H2 / Sherlock.
**Read this before touching pipeline code. Update it when you find a bug.**
Archive entries into CHANGELOG.md when the session wave is done.

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

*To add an entry: symptom (exact error string), root cause (one sentence),
fix commit, safe-to-pull verdict.*
