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

**[2026-05-20] Commit B hypothesis:** Commit B's `write_corrected_bam_parallel` partitions the BAM into per-region BAMs (≤100k ref bp each) with idempotent `.ok` sentinels before region workers start. This is architecturally the same as the `--checkpoint-dir` mitigation above. The structural resolution test (6.7M-read Han wt_R1 full-scale run on H2 16-core) was **DEFERRED** — H2 pod_smp.q had 25,600 waiting jobs at submission time, and Sherlock ControlMaster was not available (required Duo re-auth). The chrI-V subset (2,009,670 reads, 16 threads) from prior Sherlock runs already passed (61:40 wall, exit 0), validating the Commit B architecture at this scale. The full 6.7M-read test remains as Outcome A/B/C to be run in a coordinated follow-up session when queues clear. See briefing: `dev/specs/briefings/commit_b_briefing.md` §4.

---

*To add an entry: symptom (exact error string), root cause (one sentence),
fix commit, safe-to-pull verdict.*
