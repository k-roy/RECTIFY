# Bug debrief: Axis 2 — Scale-induced htslib state corruption in `rectify correct` parallel workers

**Status (2026-05-20):** OPEN. Reproduced independently on two workstreams (Han 2023 wt_R1 6.7M-read run; set2 11-sample cascade with merged BAMs 3.5-15M reads). Field reports in `AGENT_FIXES.md` entries 2 and 6.

**Asking for:** root-cause diagnosis and a code-level fix. Workarounds exist (single-aligner correct without `--aligner-bams`; Stage E.5 per-chunk-consensus; Commit B's region-pre-partitioning hypothesis) but each compromises the multi-aligner consensus design or requires architectural changes. We want to keep multi-aligner consensus AND scale to multi-million-read merged BAMs without crashes.

---

## 1. Symptom signature

Every per-chunk `rectify correct` invocation that passes `--aligner-bams` flags pointing at multi-million-read merged BAMs eventually hits one of two failure modes:

### 1a. Explicit free() crash (minority)

```
INFO - rectify.core.bam.parallel - X regions across N workers
*** Error in `python': free(): invalid next size (fast): 0x... ***
*** Error in `python3.11': free(): invalid next size (fast): 0x... ***
======= Backtrace: =========
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x2d4b9)
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x4fc67)
libcalignedsegment.cpython-311-x86_64-linux-gnu.so(+0x4fb34)
...
Processing regions:  16%|██ | 17/107 [02:05<10:37, 7.09s/it]
======= Memory map: =========
[long mmap dump]
```

Workers print the glibc abort + map dump, then either die or zombie out. Parent process keeps Slurm thinking the job is still RUNNING.

### 1b. Silent hang (majority)

```
INFO - rectify.core.bam.parallel - X regions across N workers
[nothing more — log frozen for hours]
```

`ps` on the compute node shows worker processes still at 50–99% CPU but `strace` shows no syscalls; workers are spinning in userspace inside corrupted-memory tight loops. No `free()` line in the log. Slurm reports RUNNING until walltime kill. **This is the dominant failure mode at scale.**

Both modes are the same underlying bug per the backtraces and the fact that they coexist within the same job batches.

### 1c. Important non-symptom

`samtools quickcheck` passes on every input BAM. `samtools view -c` returns correct counts. `pysam.AlignmentFile(...).fetch(region)` works **single-threaded** on the same BAMs without error. The corruption emerges only when N pysam workers fan out and fetch concurrently against the same large merged BAMs.

---

## 2. Scope and reproducibility

### 2a. What triggers it (verified)

The triggering pattern is `rectify correct` invoked with:
- Per-chunk input BAM (any size; smallest seen was 887 reads in a smoke test).
- `--aligner-bams` flag(s) pointing at one or more merged BAMs **≥ ~4M reads each**.
- Default 8-worker parallel correction (some evidence the threshold lowers with more workers).

Each worker opens its own `pysam.AlignmentFile` handle to each of the N `--aligner-bams` and does indexed `fetch(chrom, start, end)` calls per read for cross-aligner consensus comparison. N concurrent workers × M merged BAMs × random-access fetch is the smoking gun.

### 2b. Scale evidence (compiled from both workstreams)

| Workload | Reads (input) | Workers | `--aligner-bams` merged-BAM size | Result |
|---|---|---|---|---|
| Han wt_R1 1,072 reads | 1,072 | 1 | (none — single BAM input) | ✓ 31 s |
| Han wt_R1 1,072 reads | 1,072 | 8 | (none) | ✓ 9 s |
| Han wt_R1 chrI | 102,157 | 16 | (none) | ✓ 3:14 |
| Han wt_R1 chrI–V | 2,009,670 | 16 | (none) | ✓ 61:40 |
| Han wt_R1 full | 6,708,292 | 16 | (none) | ✗ free() at ~17 / 107 regions |
| set2 wt_rep2 chunk_000 1k subsample | 887 | 8 | 5 × ~3.9M merged BAMs | ✓ 30 min, workers alive, no free() (but timed out at wall, possibly hung) |
| set2 wt_rep1 chunk_000 mapPacBio | 483,001 | 8 | 5 × ~4.3M merged BAMs (post-empty-SEQ-filter) | ✗ hung-or-crash at "X regions across 8 workers" line |
| set2 cascade ~530 chunks | 280k–466k each | 8 | 5 × 3.5–15M merged BAMs | ✗ ~0 of 530 chunks produced real corrected_reads.tsv; mix of silent-hang and explicit free() |

Inference: threshold appears to be ~4M reads of indexed-BAM cross-aligner data fetched concurrently. Below that, the worker fan-out is stable. Above that, corruption is reliable.

### 2c. Empty-SEQ records as a (separately documented) trigger amplifier

The set2 cascade also hit a SECOND independent corruption path from the pre-fix mapPacBio QNAME sanitizer producing 8.4% empty-SEQ unmapped records in the merged mapPacBio BAM (AGENT_FIXES entry 1). Filtering those (entry 1's "Faster recovery" section) was **necessary** to clear prescan but **not sufficient** for correct — the scale-related corruption (this debrief's subject) persisted after filtering. They are TWO bugs that compose in nasty ways:
- Empty-SEQ records make ANY pysam linear scan of the affected merged BAM crash near-instantly (prescan-class operations).
- The scale-related corruption hits cross-aligner `fetch()` against ANY large merged BAM regardless of content cleanliness.

The next agent debugging this should treat them separately. Filtering empty-SEQ first removes the noise; what remains is the scale issue.

---

## 3. What's been ruled out (verified clean by prior investigation)

For Han wt_R1's BAM (the cleanest available test case):
- `samtools quickcheck` → exit 0
- SEQ length != QUAL length scan over full 6.7M reads → 0 hits
- QNAME length > 50 chars → 0 hits (BBmap ran with `trd=t`, qnames are bare SRR accessions)
- The `382fcc7` / `e8c8070` FASTQ sanitizer bugs (AGENT_FIXES entry 1) do NOT apply — BBmap reads FASTQ directly via `trd=t`, no rectify sanitizer in the BBmap path
- `samtools view -c` returns expected counts
- Single-thread single-process `pysam.AlignmentFile().fetch()` succeeds on every region

For set2 merged BAMs (post-empty-SEQ-filter):
- Pysam scan for QNAME issues → 0 hits
- Pysam scan for SEQ/QUAL length mismatches → 0 hits
- Single-thread `prescan` (which does sequential pysam scan, not parallel fetch) succeeded on all 11 samples in 1–7 min each

So the records themselves are valid SAM/BAM (under static checks); the bug is in dynamic state across multi-worker concurrent fetch.

---

## 4. Working hypothesis (best theory)

**The corruption is in htslib's internal state when multiple processes share a parent's opened BAM file handle via Python's multiprocessing fork.**

Concrete chain of suspicion:

1. `rectify correct` opens `pysam.AlignmentFile(merged_bam, 'rb')` in the parent process (during config / region identification phase).
2. Workers spawn via `multiprocessing` — by default `fork()` on Linux. Child processes inherit the parent's open file descriptors and any mmap'd regions.
3. The htslib C-level `htsFile*` / `hts_idx_t*` structures inside libcalignedsegment.so contain mutable state (current file position, decompression buffers, index cache pages). After fork, both parent and children share this state via copy-on-write OR (worse) directly if mmap-shared.
4. Concurrent `fetch()` calls from multiple workers mutate this shared state simultaneously. For small BAMs, the operations may complete fast enough not to collide; for large BAMs, the seek+decompress+yield cycle is long enough that two workers can interleave inside the same C-level buffer.
5. The corruption manifests as glibc `free()` errors when one worker frees a buffer another worker is still reading, or as silent infinite loops when an index node pointer chain becomes circular.

This is **strongly consistent** with the scale dependency: large BAMs have larger compressed blocks and more index nodes, so concurrent decompress/seek operations span more time and more memory, increasing the collision window.

**Why it's HTSLIB-specific (not Python multiprocessing in general):** the backtraces show all faulting addresses inside `libcalignedsegment.cpython-311-x86_64-linux-gnu.so` — that's pysam's compiled C extension wrapping htslib. The corruption is in C, not in Python's pickling or queue.

### 4a. Why Commit B's per-region BAM pre-partitioning may fix it

Commit B (in flight on a separate branch) changes the architecture so that each worker reads from its own per-region BAM slice (≤100 kb ref-space), produced by the PARENT process before workers start, instead of all workers fetching from the same shared merged BAM handle. Each worker's pysam handle is independent and unique to it; no inherited state from parent; no concurrent access to the same C structures.

That should resolve the bug at the architectural level — but it's unproven for the full 6.7M Han wt_R1 run (the structural resolution test was deferred at the last attempt; see `dev/specs/briefings/commit_b_briefing.md` §4).

---

## 5. Investigation paths (concrete experiments)

These are ordered by expected information yield per hour of investigator time.

### 5a. Single-process baseline (1 hour)

Reproduce the Han wt_R1 6.7M-read crash with `--threads 1`. If it crashes single-threaded, the bug is in single-process pysam against a large BAM — much narrower search. If it succeeds, the bug is provably in worker fan-out and we narrow to that.

```bash
# Run the failing rectify correct invocation with --threads 1 instead of --threads 16:
python -m rectify correct <han_wt_r1_merged.bam> \
  --genome ... --annotation ... \
  --junction-pool-cache ... --variant-scan-cache ... \
  --aligner-bams "bbmap:<bbmap_merged>" --aligner-bams "bwa:<bwa_merged>" \
  --threads 1 \
  -o out.tsv --write-corrected-bam out.bam
```

Memory says this passed in 31 sec on 1,072 reads. Untested at 6.7M with `--threads 1`.

### 5b. Spawn vs fork (1 hour)

`multiprocessing` defaults to `fork` on Linux. Change to `spawn` — children get fresh interpreters with no inherited state.

Where to patch: `rectify/core/bam/parallel.py` worker pool initialization. Add `multiprocessing.set_start_method('spawn', force=True)` before pool creation. Re-run the failing test.

If `spawn` succeeds, the bug is fork-related state inheritance — supports the hypothesis.

### 5c. Per-worker pysam handle, not parent-inherited (2 hours)

Even with `fork`, parent shouldn't `open` the merged BAMs before workers start. Verify in `rectify/core/bam/parallel.py` that:
- Parent process does NOT instantiate `pysam.AlignmentFile(merged_bam)` before forking.
- Each worker, in its init function, opens its own `pysam.AlignmentFile(merged_bam)` (fresh fd, fresh htslib state).

If parent opens BAMs first and workers fetch via inherited handles, fix this. If workers already open their own, this hypothesis is weaker and we move on.

### 5d. mmap interaction (2 hours)

htslib uses mmap for BAM/BAI index access on some platforms. Check via `cat /proc/<worker_pid>/maps | grep .bam` while the job is running (or just-before-crash) — verify each worker has its own private mmap region per BAM or shared.

If shared: confirmed concurrent-mmap-mutation risk. Open `pysam.AlignmentFile` with `require_index=True` and check if the index loading goes through mmap.

### 5e. Glibc allocator vs jemalloc (2 hours)

Glibc `free()` is what's aborting. Test with `LD_PRELOAD=libjemalloc.so` or `LD_PRELOAD=libtcmalloc.so` — if these accept the same access pattern that glibc rejects, the corruption is a borderline-correct heap operation that glibc's hardened allocator catches and other allocators tolerate. Doesn't fix the bug but localizes it.

### 5f. ASAN / Valgrind on a smaller reproducer (4 hours)

Build a minimum repro: same BAMs, smallest input chunk that still triggers. Run under `LD_PRELOAD=libasan.so.6 python ...` or under `valgrind --tool=memcheck`. Should pinpoint the exact line in libcalignedsegment.

This is the highest-yield experiment if (5a)-(5e) don't already localize the bug.

### 5g. Upstream pysam / htslib bug report search (1 hour)

Search pysam issues + htslib issues for `concurrent fetch`, `fork multiprocessing`, `free invalid next size`, `large indexed BAM`. There may be a known bug + fix in a later pysam/htslib version. Verify the conda env's pysam version against latest stable; if behind, upgrade and re-test.

### 5h. Reduce worker count to characterize the threshold (1 hour)

Run the failing test with `--threads 16, 8, 4, 2, 1`. Map out the curve. Memory says 16 explicit-crashes, 8 silent-hangs. Does 4 silent-hang, complete, or also crash? Does 2? Threshold characterization narrows root cause.

---

## 6. Code locations to inspect

### 6a. RECTIFY code

- `rectify/core/bam/parallel.py` — worker pool, fan-out logic, the module whose log line "X regions across N workers" precedes every crash. This is THE module from `0024fa3 refactor(bam): split bam_processor.py into regions/output/variant_scan/parallel` (the refactor that surfaced the symptom; the bug pre-dates the refactor but moved here).
- `rectify/core/commands/correct_command.py` — caller of the parallel module, where `--aligner-bams` are parsed and passed in. Check whether merged BAMs are opened in parent before fork.
- `rectify/core/consensus/*.py` — cross-aligner consensus logic that the workers call into and that's why they need the merged BAMs at all.

### 6b. Compiled extensions in the backtrace

- `libcalignedsegment.cpython-311-x86_64-linux-gnu.so` — pysam's wrapper around htslib's `bam_*` API. Offsets shown in entry 2 backtrace: `+0x2d4b9`, `+0x4fc67`, `+0x4fb34`. Disassembly of these offsets against the conda env's pysam build would tell you the exact C function and line.

### 6c. Environment

- conda env at `/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify` on Sherlock and `/u/project/guillom/shared/envs/rectify` on H2.
- pysam version: check via `python -c "import pysam; print(pysam.__version__)"`.
- htslib version: bundled in pysam wheel; check via `python -c "import pysam; print(pysam.libchtslib.__file__)"` + ldd-style introspection.
- glibc version: `ldd --version` on the compute node.

---

## 7. Reproducer setup

### 7a. Sherlock — easiest

```
/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_tfiiib_rep1/chunks/
```

Has:
- 24 per-chunk per-aligner BAMs (input)
- 5 merged BAMs at 12-15M reads each (the failing input to `--aligner-bams`)
- prescan PKLs ready
- Slurm wrapper scripts in place

Empty-SEQ filter has been applied to merged mapPacBio BAM only (so it's the cleanest set2 sample for studying axis 2 specifically). Sample `run_array_correct_minimap2.sh` will trigger the bug within ~20-30 min of starting.

### 7b. Han 2023 — already has a chrVI–XVI subset

See `AGENT_FIXES.md` entry 2 for paths. Per the scale table, the chrI–V subset (~2M reads) passes, full 6.7M fails. Run the chrVI–XVI subset (~4.7M) to confirm threshold; bisect chromosomes to find the smallest crashing subset.

### 7c. Minimum data smoke

The 887-read subsample I made for the post-QNAME-fix smoke test is at:
- Sherlock: `/oak/stanford/groups/larsms/Users/kevinroy/correct_test_1k.sh` and the subsampled BAM in `/tmp/correct_test_1k_<jobid>/input_1k.bam` (regenerable from the script).

But: this DOESN'T reproduce axis 2 — workers stayed alive at high CPU for 30 min with no free() errors (timed out at wall). To trigger axis 2 reliably the merged BAMs must be ≥4M reads, regardless of input chunk size.

---

## 8. Definition of done

The fix is complete when:

1. `rectify correct` with `--aligner-bams` pointing at 4M+ read merged BAMs and `--threads ≥ 8` completes on:
   - The Han wt_R1 6.7M-read merged BAM (the original reproducer)
   - All 11 set2 cascade samples (3.5–15M-read merged BAMs)
   - The cdna/UMI workflow at its current scale (whatever that is — verify)
2. No `free()` errors in any log.
3. No silent hangs (workers either make progress or exit cleanly).
4. Output `corrected_reads.tsv` files have plausible non-empty content (validate against single-aligner output by sampling 1k reads and confirming cross-aligner consensus matches expected ratios).
5. Performance regression budget: < 2× the current per-chunk wall (single-threaded baseline is the lower bound).

If Commit B's region pre-partitioning is the resolution path, the same DoD applies — Commit B is just one possible implementation.

---

## 9. Related material

- `AGENT_FIXES.md` entry 1 — mapPacBio QNAME sanitizer + empty-SEQ filter recovery (axis 1, the noise to clear before debugging axis 2)
- `AGENT_FIXES.md` entry 2 — original Han 2023 wt_R1 reproducer field report
- `AGENT_FIXES.md` entry 6 — set2 cascade field report (same bug, smaller scale, more samples)
- `dev/specs/briefings/commit_b_briefing.md` — the in-flight architectural fix candidate
- `project_status_markdowns/DRS_CPA_PROJECT_STATUS.md` session 5 handoff — operational context

---

## 10. Out of scope for this briefing

- The empty-SEQ filter (axis 1) is solved in field; don't reinvestigate.
- The bedgraph 1-bp left shift in `rectify/core/analyze/bedgraph.py` is a separate axis with its own M1 WT fix waiting to be committed.
- The mapPacBio QNAME sanitizer fixes (`382fcc7`, `e8c8070`) are in place; further FASTQ-level sanitizer work is not in scope.
- The Sherlock owners-partition preemption thrashing is operational, not code.

Each of these may amplify or complicate the symptoms you'll see in the field, but none of them are the root cause of axis 2.
