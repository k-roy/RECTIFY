# deSALT: index requirement, duplicate alignments, deterministic crashes, threading

deSALT (De Bruijn graph / deBGA-based) is a high-sensitivity long-read splice
aligner. RECTIFY runs it via `run_desalt()` in `multi_aligner.py` as a Tier-2 /
junction aligner (`--junction-aligners deSALT`). It is sensitive but **fragile**:
it crashes deterministically on some inputs, duplicates every alignment, and
cannot be forked inside a multithreaded process. RECTIFY works around all of these;
this doc explains the quirks so the workarounds aren't mistaken for bugs.

---

## Vendored binary — Linux x86_64 only

deSALT is not pip-installable. RECTIFY ships a vendored binary at
`rectify/data/bin/linux_x86_64/deSALT` and resolves it via
`_get_vendored_desalt()` when `deSALT` is not on `PATH`. There is **no macOS
binary** — deSALT does not run on M1; it is a cluster-only aligner (Hoffman2 /
Sherlock). Install/expose with:

```bash
rectify install-aligners --desalt
```

On Sherlock it is exposed via a `~/bin/deSALT` symlink into the conda env.

---

## Pre-built RdBG index required

deSALT needs a De Bruijn graph index built ahead of time:

```bash
deSALT index <ref.fa> <index_dir>
```

`run_desalt()` looks for a `desalt_index/` directory adjacent to the genome if
`index_path` is not given. **An empty `desalt_index/` placeholder is a trap:** it
is picked up silently and deSALT then fails late with a confusing error because
the deBGA main reference output (the actual index files) is missing. Ensure the
index directory is fully built, not just present.

---

## Duplicate-alignment bug — RECTIFY dedups

deSALT outputs **each alignment N times**, where N is the number of secondary
alignment slots (`-N`, default 4). `_dedup_desalt_bam()` removes them by keeping
the first occurrence of each `(query_name, flag, reference_name, reference_start,
cigarstring)` tuple:

```
deSALT dedup: <n_total> → <n_kept> alignments (<n_removed> duplicates removed)
```

If you read deSALT BAMs that have **not** passed through `run_desalt()`, expect
~N× inflated alignment counts — dedup before any per-alignment counting.

---

## Deterministic SIGSEGV / OOM crashes — empty-BAM fallback

deSALT v1.5.6 segfaults (**SIGSEGV, exit 139**) or is **OOM-killed (exit 137)**
during its second-pass `Loop-ProcessReads` when certain pseudo-exon structures are
inferred from the input. The crash is **deterministic for a given input batch —
retries never succeed.** RECTIFY treats these exits as a recognized crash and
emits a valid empty name-sorted BAM so the merge proceeds with the other aligners:

```python
# multi_aligner.py
_DESALT_CRASH_EXITS = frozenset({139, 137, -11, -9})
# 139/137 = shell (128+signal); -11/-9 = Python subprocess negative signal
```

Consequently `align_command` **honors empty deSALT BAMs** (the usual
non-empty/size>2000 check is bypassed for deSALT) so a crash fallback is not
re-flagged as failure. Upstream bug: github.com/ydLiu-HIT/deSALT/issues/49.

A deSALT crash is therefore **not fatal** to a run — you simply lose deSALT's
evidence for that sample and the consensus proceeds on the remaining aligners.

---

## `-G` (GTF) flag triggers SIGSEGV on yeast GTF

Passing the annotation to deSALT via its `-G` flag causes a SIGSEGV on the SGD
yeast GTF, so `run_desalt()` does not feed the GTF through `-G`. deSALT still
performs splice-aware alignment from the RdBG index; annotation-guided behavior is
intentionally not used.

---

## Threading: "double free or corruption" when forked in a multithreaded process

deSALT crashes with `double free or corruption` when launched (forked) from inside
a multithreaded Python process. `align_command` therefore runs deSALT
**sequentially, after** the `ThreadPoolExecutor` parallel-aligner pool has closed —
never as one of the pooled futures:

```python
parallel_batch   = [a for a in remaining if a != 'deSALT']
sequential_batch = [a for a in remaining if a == 'deSALT']
# ... run parallel_batch in the pool, THEN deSALT alone, single-threaded launch
```

If you add deSALT to a custom parallel-dispatch path, keep it out of the thread
pool or it will double-free.

---

## Verifying deSALT works

```bash
# Cluster only (no macOS binary). Requires a pre-built desalt_index/.
rectify align --reads <reads.fastq.gz> --genome <genome.fa> \
    --annotation <annotation.gff.gz> --aligners deSALT \
    --output /tmp/desalt_smoke --threads 8
samtools flagstat /tmp/desalt_smoke/*.deSALT.bam
```

Expected: nonzero mapped count. An **empty** BAM means deSALT hit the deterministic
SIGSEGV/OOM crash (check logs for exit 139/137) — that is a tolerated fallback, not
a pipeline failure.

---

## Failure modes quick-reference

| Symptom | Cause | Fix |
|---------|-------|-----|
| `deSALT not found` and no alignment | not on PATH, no vendored binary for platform | `rectify install-aligners --desalt` (cluster only; no macOS binary) |
| Confusing late index error | empty/placeholder `desalt_index/` | rebuild fully: `deSALT index <ref.fa> <index_dir>` |
| ~N× inflated alignment counts | deSALT duplicate-output bug (`-N` slots) | use `run_desalt()` (auto-dedups) or dedup on (name,flag,chrom,pos,cigar) |
| Empty `.deSALT.bam`; exit 139/137 in logs | deterministic SIGSEGV/OOM in `Loop-ProcessReads` | tolerated fallback (other aligners proceed); upstream issue #49 |
| SIGSEGV when annotation passed | deSALT `-G` flag crashes on yeast GTF | RECTIFY omits `-G`; don't add it |
| `double free or corruption` | deSALT forked inside a multithreaded process | run deSALT sequentially, outside the thread pool |

---

## Forensic crash investigation (2026-05-18/19)

*Original doc: Drive `docs/desalt_crash_investigation_handoff.md`, now archived here.*


**Created:** 2026-05-19  
**Updated:** 2026-05-19 (post-minsearch)  
**Status:** Investigation complete — ready to file upstream  
**Upstream issue:** https://github.com/ydLiu-HIT/deSALT/issues/49  
**Binary:** `~/bin/deSALT` → vendored `he4a0461_5` bioconda build (v1.5.6); conda env copy at `anaconda3/envs/rectify/bin/deSALT`

---

## What we know

The conda binary (`he4a0461_5`) eliminates the SIGSEGV on the vast majority of chunks
but not all. Two chunks in set2 still crash:

| Sample | Chunk | FASTQ reads (total) | Clean reads (after `_clean_fastq`) | Exit code | Adjacent chunks |
|--------|-------|---------------------|------------------------------------|-----------|-----------------|
| wt_tfiiib_rep3 | chunk_005_of_016 | 478,575 | 420,447 | -11 (SIGSEGV) | chunk_004: 478,575 → PASS, chunk_006: 478,575 → PASS |
| rna15_rep3 | chunk_003_of_004 | 236,498 | 203,832 | -11 (SIGSEGV) | chunk_002: 236,499 → PASS |

Adjacent chunks have **identical read counts** — the crash is purely composition-dependent,
not a read-count threshold.

Both failures happened on May 18 during re-alignment runs (jobs 25388629_53 and
25388631_15). The `Loop-ProcessReads` SIGSEGV is documented in upstream deSALT#49 as
triggered by specific pseudo-exon graph construction.

For rna15_rep3/chunk_003, the original May 15 run produced a 444 MB `.deSALT.sam` that
is also truncated (CIGAR/query length mismatch), suggesting that run was killed by SIGTERM
(time limit) mid-write — this chunk has never produced a valid deSALT alignment.

---

## Current production state

Both chunks use the 4-aligner consensus fallback (mapPacBio + minimap2 + gapmm2 + uLTRA).
This is handled automatically by the correction pipeline's empty-BAM path. No manual
intervention is needed to proceed with the set2 correction runs.

---

## Investigation results (2026-05-18/19)

All experiments used `rna15_rep3/chunk_003` (203,832 clean reads) and
`~/bin/deSALT` (vendored, he4a0461_5). Work directory: `/scratch/users/kevinroy/desalt_bisect_33632/`.

### Finding 1: crash is data-dependent, NOT a thread race

Running `repro_14900.fastq` at **-t 1** (single thread) still crashes with SIGSEGV in
`Loop-ProcessReads` 1st loop. This rules out a thread-scheduling race condition —
the bug is deterministically triggered by specific read content.

### Finding 2: crash is non-monotonic with read count

The full clean FASTQ was split into prefix cuts (first N reads) and tested at -t 1.
Pass/fail by prefix length:

| Prefix (reads) | Result |
|----------------|--------|
| 0 – 14,850 | PASS |
| **14,860** | **CRASH** |
| 14,870 – 14,880 | PASS |
| **14,890 – 14,950** | **CRASH** |
| **15,000 – 19,000** | **CRASH** |

At coarser scale (-t 2):

| Prefix (reads) | Result |
|----------------|--------|
| 10,000 | PASS (empty input — file was 0 bytes; trivially passes) |
| **20,000** | **CRASH** |
| 30,000 | PASS |
| **40,000** | **CRASH** |
| **50,000–203,832** | **CRASH** (all tested lengths) |

The non-monotonic pattern — adding reads can BOTH trigger and cure the crash — means
the bug depends on **which reads are present together**, not just how many. Specific
reads force the pseudo-exon graph into a state that corrupts memory during second-pass
alignment. Other reads change the graph topology in a way that avoids the corrupt state.

### Finding 3: exact trigger reads identified — both map to convergent/divergent gene boundaries

Per-position testing of window 1 (14850–14860) and window 2 (14880–14890) shows a clean
boundary at each edge: c14851–c14859 all PASS; c14860 CRASH. c14881–c14889 all PASS;
c14890 CRASH. The crash is introduced by exactly **one specific read** at each boundary:

| Trigger | Read index (0-based) | UUID | Length | Locus | Strand |
|---------|---------------------|------|--------|-------|--------|
| #1 | 14859 | `d803d4d8-2a52-4aaa-9c14-adac9e05c377` | 1,080 bp | chrII:690,260–691,347 | minus |
| #2 | 14889 | `e46a457a-cfa7-46f6-8491-1c022966010e` | 1,326 bp | chrII:711,562–712,892 | plus |

**Both trigger reads map to compact intergenic boundaries between opposite-strand gene pairs:**

- **Trigger #1** spans the shared terminal region of **YBR235W** (+, 686,820–690,354) and
  **YBR236C** (−, 689,753–691,740). Their 3' ends overlap around 690,000–690,354 — a
  convergent arrangement where the read simultaneously provides exon-boundary evidence on
  both strands.

- **Trigger #2** starts in the 7-bp shared intergenic region of **YBR245C** (−, 5' end at
  711,539) and **YBR246W** (+, 5' end at 711,533). A divergent gene pair where both
  promoters are within 7 bp of each other.

Neither gene has annotated introns. The pattern: reads that span a strand-switching
boundary between two closely packed genes simultaneously suggest pseudo-exon split points
on both strands, likely forcing deSALT's pseudo-exon graph into a degenerate topology.

**Trigger reads do not crash in isolation.** Running only the 2 trigger reads through
deSALT exits 0. The crash requires the ~14,859-read background to accumulate enough
pseudo-exon graph state first; the trigger read then pushes it into the corrupt state.

**Minimum known crashing prefix:** `cleaned.fastq` reads 0–14859 (14,860 reads). The
`repro_14900.fastq` is slightly larger but also reliable.

Crash location in every repro run:
```
[Loop-ProcessReads] The 0th loop: 7,452 reads aligned ← always completes
[Loop-ProcessReads] The 1st loop: 7,448 reads ← SIGSEGV here
```

### Finding 4: conda binary also crashes — bug is not fixed in he4a0461_5

Running `repro_14900.fastq` through the conda binary
(`anaconda3/envs/rectify/bin/deSALT`) at -t 1 also exits 139 (SIGSEGV).
The he4a0461_5 build suppresses crashes on most chunks but does **not** fix
the underlying memory corruption for this specific read set.

### Finding 5: crash is content-dependent, not order-dependent

Shuffling `repro_14900.fastq` (seed=42, random.shuffle on 4-line read groups)
and running the vendored binary at -t 1 also exits 139. Any ordering of these
14,900 reads crashes. The bug is triggered by which reads are present together
— the specific pseudo-exon graph topology they force — not by the order they
arrive at `Loop-ProcessReads`.

Job 25452208 ran both tests on sh03-08n22 on 2026-05-19; logs at
`/scratch/users/kevinroy/desalt_bisect_33632/final_tests_25452208.{out,err}`.

### Finding 6: crash requires the specific original read set — random subsampling does not reproduce it

To find a smaller shareable reproducer, 30 tests were run (6 background sizes × 5 random
seeds): N ∈ {7000, 3000, 1000, 300, 100, 50} reads randomly sampled from the 14,859-read
background pool, each with trigger read `d803d4d8` appended. All 30 tests exited 0.

| Background N | Seeds crashed / 5 |
|-------------|------------------|
| 7,000 | 0/5 |
| 3,000 | 0/5 |
| 1,000 | 0/5 |
| 300 | 0/5 |
| 100 | 0/5 |
| 50 | 0/5 |

**Conclusion:** the crash is locked to the specific ordered composition of reads 0–14858.
A random draw of even 7,000 reads (≈half the background) never reproduced the crash across
5 seeds. The critical graph state depends on a precise combination of reads, not just any
large-enough random subset. The smallest confirmed reproducer remains `repro_14900.fastq`
(14,900 reads, gzipped to ~10 MB at
`/scratch/users/kevinroy/desalt_bisect_33632/desalt_minimal_reproducer.fastq.gz`).

Job 25457452 ran the minsearch on sh02-04n20 on 2026-05-19; logs at
`/scratch/users/kevinroy/desalt_bisect_33632/desalt_minsearch_25457452.{out,err}`.

### Finding 7: wt_tfiiib_rep3/chunk_005 confirmed crash at -t 1, with matching random and prefix search results

The second crashing chunk (wt_tfiiib_rep3/chunk_005, 420,447 clean reads) was investigated
in full on 2026-05-19 (job 25461125, sh03-08n22, AMD Milan).

**Crash confirmation at -t 1:** full cleaned FASTQ exits -11 (SIGSEGV). Original production
crash was at -t 2 (job 25388629_53). Confirmed single-threaded → not a race condition on
either chunk.

**Random subsampling (6 sizes × 5 seeds = 30 tests):** 0/30 crashed, identical to rna15.
The crash requires the specific composition of the full chunk.

**Prefix search (13 sizes, 1k–420k):**

| Prefix (reads) | Result |
|---------------|--------|
| 1,000 – 30,000 | PASS |
| **50,000** | **CRASH** |
| **75,000 – 420,000** | **CRASH** |

Non-monotonic: 30k PASS → 50k CRASH. Minimum crashing prefix: **50,000 reads**
(on AMD Milan sh03-08n22).

### Finding 8: wt_tfiiib crashes on both Intel and AMD; rna15 confirmed AMD-only

A targeted prefix test on an Intel Broadwell node (sh02-09n06, job 25470388) and a
refined AMD Milan test (sh03-08n22, job 25470332) compared crash thresholds across
architectures.

**wt_tfiiib_rep3/chunk_005 — crashes on BOTH architectures:**

| Node | Architecture | Sizes tested | Result |
|------|-------------|-------------|--------|
| sh02-09n06 | Intel Broadwell | 20k, 30k, 40k, 50k, 75k, 100k, 200k, 420k | ALL CRASH |
| sh03-08n22 | AMD Milan | 20k, 30k, 40k, 50k, 75k, 100k, 200k, 420k | 20k–30k PASS; 40k–420k CRASH |

Intel is **more sensitive**: crashes from ≤20,000 reads. AMD requires ≥40,000 reads (non-
monotonic: 20k–30k pass, 40k+ crash). Smallest confirmed cross-architecture reproducer:
**20,000 reads** from the cleaned wt_tfiiib FASTQ (crashes on both Intel and AMD).

**rna15_rep3/chunk_003 — confirmed AMD-only in all tested conditions:**
- CRASH on sh03-07n10 (AMD Milan, job 25452947) at 14,860 reads
- CRASH on sh03-08n22 (AMD Milan, job 25452208) at 14,900 reads
- PASS on sh02-02n12 (Intel, job 25461005) across all prefix sizes including 14,860
- PASS on sh02-04n20 (Intel, job 25457452) random subsampling
- PASS on sh03-08n22 (AMD Milan, job 25469142) at prefix sizes through 14,860 — sh03-08n22
  needs repro_14900 (14,900 reads) to crash, not just 14,860

**Interpretation:** the wt_tfiiib crash is a more severe memory corruption that causes SIGSEGV
on any tested architecture. The rna15 crash is marginal: the corrupted address only reaches
an unmapped page on specific AMD Milan nodes with specific read counts. Developers on Intel
hardware can use the wt_tfiiib 20k-read reproducer directly.

---

## What remains

**File upstream comment on deSALT#49** (see draft below — all experiments complete).

---

## Draft upstream comment for deSALT#49

> **Summary:** deSALT v1.5.6 (bioconda `he4a0461_5`) still crashes with SIGSEGV in
> `Loop-ProcessReads` on specific ONT DRS data, even at `-t 1`. We confirmed the crash
> on **two independent chunks from different experimental conditions**, on both AMD EPYC
> Milan and Intel Broadwell hardware (Sherlock HPC, Stanford). Both the manually compiled
> binary and a separately installed conda copy of the same version are affected.
>
> **Reproducers (happy to share either FASTQ directly):**
> - **rna15_rep3/chunk_003** — 203,832 clean reads total; bisected to a **14,860-read
>   minimum prefix** via per-read binary search. A 14,900-read prefix reproduces reliably.
>   Index: *S. cerevisiae* R64-5-1.
> - **wt_tfiiib_rep3/chunk_005** — 420,447 clean reads total; prefix search found crashes
>   from **50,000 reads up to the full chunk** (non-monotonic: 30k PASS → 50k CRASH).
>   Same index.
>
> These chunks come from different yeast experiments processed in the same pipeline run;
> 2 of ~200 chunks crash; the remaining 198 align cleanly.
>
> **Key findings (both chunks unless noted):**
>
> 1. **Two independent chunks, different conditions** — rna15_rep3 (RNA Pol III depletion)
>    and wt_tfiiib_rep3 (TFIIIb depletion) both crash on the same binary, different yeast
>    cultures, different read compositions. This is not an isolated edge case.
>
> 2. **Single-threaded** (`-t 1`) reproduces on both — not a thread-scheduling race. Bug
>    is data-dependent.
>
> 3. **Non-monotonic crash** — not a read-count threshold. rna15: first 14,850 reads exit 0,
>    first 14,860 crash, first 14,870 exit 0 again. wt_tfiiib: first 30,000 reads exit 0,
>    first 50,000 crash. Adding reads can both expose and hide the corruption.
>
> 4. **Content-dependent, not order-dependent** — shuffling the rna15 14,900-read FASTQ
>    (Python `random.shuffle`, seed=42) still crashes. However, random subsets do **not**
>    crash: 60 tests across both chunks (6 subset sizes × 5 seeds × 2 chunks) all exited 0.
>    The corruption requires the specific ordered composition of reads, not just a
>    large-enough random draw from the same pool.
>
> 5. **Crash location** — always in `Loop-ProcessReads` **1st loop** (rna15). The 0th
>    loop always completes cleanly.
>
> 6. **Exact trigger reads identified (rna15).** Per-read binary search pinpointed two
>    reads that each introduce the crash when appended to a passing prefix:
>    - `d803d4d8` (1,080 bp, chrII:690,260, minus strand) — maps to a region where the ends
>      of two opposite-strand genes overlap: **YBR235W** (+ strand, ending at 690,354) and
>      **YBR236C** (− strand, ending at 689,753) share ~600 bp of overlap at their 3' ends.
>      The read spans this overlap zone and runs along active gene sequence on both strands.
>    - `e46a457a` (1,326 bp, chrII:711,562, plus strand) — maps to the 7-bp gap between
>      **YBR245C** (− strand) and **YBR246W** (+ strand), a gene pair whose start sites face
>      each other only 7 bp apart. The read runs through that gap and into both flanking genes.
>
>    Neither locus has annotated introns. In both cases, a single read spans a narrow window
>    where gene sequence exists on **both DNA strands** — so deSALT must process it as
>    potentially informative for pseudo-exon construction on either strand. Our hypothesis:
>    these reads cause the 0th loop to write more pseudo-exon candidates than the graph
>    structure was allocated to hold (because both strands of the same window each contribute
>    candidates), producing a silent out-of-bounds write. The 1st loop then dereferences that
>    corrupted graph data → SIGSEGV. The non-monotonic pattern (adding reads can cure the
>    crash) is consistent with heap corruption: the overwritten bytes land at different
>    addresses depending on allocation order, sometimes in mapped memory (no crash) and
>    sometimes in an unmapped page (SIGSEGV).
>
>    Trigger reads do **not** crash in isolation — the ~14,859-read background must accumulate
>    graph state first.
>
>    **Suggested debugging path:** running the 14,900-read rna15 FASTQ under **AddressSanitizer**
>    (`-fsanitize=address`) should immediately identify the bad write. We'd expect it to be in
>    the 0th loop, in the code that inserts new pseudo-exon candidates — likely an array indexed
>    by node count or genomic position that lacks a bounds check when both strands of the same
>    narrow region each contribute candidates. A targeted fix: before inserting a new node into
>    the candidate array, verify the index is within the allocated size; if not, grow the array
>    dynamically or skip the insertion. The `.bin1pass_anchor` files captured at the crash point
>    are available for graph-state inspection if useful.
>
> 7. **Both binaries affected** — the manually compiled binary (`~/bin/deSALT`, he4a0461_5)
>    and a separately installed conda copy of the same version both exit 139. The he4a0461_5
>    build resolves 198/200 chunks but not these two.
>
> 8. **Architecture note** — crashes confirmed on both Intel Broadwell and AMD EPYC Milan:
>    - **wt_tfiiib_rep3/chunk_005** crashes on both architectures. Intel is *more* sensitive:
>      all tested prefix sizes from 20,000 reads crash on Intel Broadwell (sh02-09n06), while
>      AMD Milan (sh03-08n22) passes at 20k–30k but crashes from 40k up (non-monotonic). A
>      **20,000-read prefix** from the wt_tfiiib cleaned FASTQ crashes reliably on both.
>    - **rna15_rep3/chunk_003** is AMD-only in all tested conditions: SIGSEGV on two AMD Milan
>      nodes (sh03-07n10, sh03-08n22) but all-PASS across all Intel nodes tested. Consistent
>      with a marginal memory corruption where the SIGSEGV fires only when the corrupted
>      address lands in an unmapped page — AMD Milan's memory layout makes this more likely
>      for this read set.
>    - **Developers on Intel hardware** can reproduce directly using the wt_tfiiib 20,000-read
>      prefix; AMD hardware needed for rna15.
>
> First-pass `.bin1pass_anchor` files are preserved at the rna15 crash point if useful
> for graph-state inspection. Happy to share FASTQs or provide additional details.

---

## Files on Sherlock

| File | Path |
|------|------|
| Reproducer FASTQ (14,900 reads, raw) | `/scratch/users/kevinroy/desalt_bisect_33632/repro_14900.fastq` |
| Reproducer FASTQ (14,900 reads, gzipped) | `/scratch/users/kevinroy/desalt_bisect_33632/desalt_minimal_reproducer.fastq.gz` |
| First-pass bin (crash state) | `/scratch/users/kevinroy/desalt_bisect_33632/repro.bin1pass_anchor.{lines,pos}` |
| Full cleaned FASTQ (rna15_rep3/chunk_003) | `/scratch/users/kevinroy/desalt_bisect_33632/cleaned.fastq` |
| Bisect logs directory | `/scratch/users/kevinroy/desalt_bisect_33632/` |
| Failing FASTQ (wt_tfiiib_rep3) | `/oak/.../v3_20260429/set2_cpa_machinery/wt_tfiiib_rep3/chunks/wt_tfiiib_rep3_trimmed_chunk_005_of_016.fastq.gz` |
| Failing FASTQ (rna15_rep3) | `/oak/.../v3_20260429/set2_cpa_machinery/rna15_rep3/chunks/rna15_rep3_trimmed_chunk_003_of_004.fastq.gz` |
| wt_tfiiib cleaned FASTQ | `/scratch/users/kevinroy/desalt_tfiiib_bisect/cleaned.fastq` |
| wt_tfiiib investigation logs | `/scratch/users/kevinroy/desalt_tfiiib_bisect/wt_tfiiib_search_25461125.{out,err}` |
| Crash logs (production) | `/oak/.../wt_tfiiib_rep3/chunks/logs/25388629_53.err`, `.../rna15_rep3/chunks/logs/25388631_15.err` |
| deSALT index (yeast) | `/oak/.../software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/desalt_index` |

---

## Priority

Low — the empty-BAM fallback handles these chunks gracefully; 4-aligner consensus for
2 of ~200 chunks has negligible impact on the final output. All experiments are complete.
The only remaining step is posting the upstream comment (draft above) to deSALT#49.
