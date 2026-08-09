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

> **Current vendored build (2026-06).** The binary shipped in the repo is now the
> bioconda build `desalt-1.5.6-h577a1d6_7` (md5 `e923d86650fb6b677bb48f3fb53b0d5f`),
> reported md5-identical across the repo, Hoffman2, and Sherlock. It superseded the
> earlier `he4a0461_5` build referenced in the forensic crash investigation below;
> the SIGSEGV/OOM empty-BAM fallback (next section) still applies, so the crash logs
> recorded under the older build remain representative of current behaviour.

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

**Relation to the panel-wide duplicate hazard.** deSALT's `-N`-slot inflation is
the one case `rectify align` physically dedups (on `(name,flag,chrom,pos,cigar)`).
This is deSALT-specific and runs *inside* `rectify align` only — it does not
protect an external BAM handed straight to `rectify correct`. `rectify correct`
skips `is_secondary`/`is_supplementary` but does **not** dedup identical
**primary** records or honor `0x400`. See the canonical cross-aligner table and
the by4742 incident in [minimap2.md](minimap2.md#-duplicate-primary-alignments--2-double-counted-3-ends-external-bam-hazard).

---

## ✅ ROOT-CAUSED AND FIXED (2026-08-09) — `1.5.6-chanfreau1`, see `planning/639`

The deterministic SIGSEGV/SIGABRT documented below (upstream deSALT#49) has been
**root-caused and patched**. It is a **heap buffer overflow of the reference-window
buffer `tseq`** in `align_core_primary` (`src/aln_2pass.c`): `tseq` is sized once
per read from the pseudo-exon key span, but the per-anchor gap-fill loop recomputes
fill windows from the anchors themselves and overruns it; because `tseq` lives in a
kalloc pool the write lands on neighbouring pool allocations and surfaces later as
the composition/architecture-dependent crash. Two amplifiers: reference reads past
the packed-genome end, and window-length underflow on overlapping consecutive
anchors (the two-strand compact-locus geometry in the forensics below).

**The fix** is a 6-commit patch series on `chanfreau/fix-v1.5.6` (M1 `~/work/deSALT`,
off tag `v1.5.6`): grow `tseq`/`junc` before every fill (`TSEQ_FIT`), clamp all
reference reads at `reference_len`, guard degenerate windows, atomic-max the
`max_exon_num_per_read`/`readlen_max` seeding race, bound `fix_cigar`, and a final
CIGAR-vs-SEQ consistency gate that marks any inconsistent read unmapped rather than
emit malformed SAM. Version string bumped to `1.5.6-chanfreau1`.

**Verified** (`planning/639`): completes the full `ysh1_rep1` DRS sample
(1,224,892 records; vanilla crashed), byte-identical to bioconda on the
9,401/9,402 shared records of the rna15 reproducer's pre-crash output (extra
records are pure recovery), 0 malformed CIGARs, `samtools quickcheck` PASS, and a
non-empty `.deSALT.bam` through the real `run_desalt` chain. Builds and runs on H2,
Sherlock, and SCG. Portable build: `make CFLAGS="… -fcommon -gdwarf-4 …"
LIBS="-lm -lz -lpthread"` (the `-fcommon` works around header-defined globals; the
vestigial `-lgomp` is dropped).

⇒ With `1.5.6-chanfreau1` the empty-BAM crash fallback below **should no longer
fire** in practice; it remains as defence-in-depth. The `_DESALT_CRASH_EXITS` set
should also include `134`/`-6` (SIGABRT) — the whole-sample crash observed in
`[[632]]` exited 134, which the legacy set does not cover.

## Deterministic SIGSEGV / OOM crashes — empty-BAM fallback (legacy behaviour, pre-fix)

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
| `double free or corruption` at deSALT launch | deSALT forked inside a multithreaded process | run deSALT sequentially, outside the thread pool |
| `ValueError: file has no sequences defined (mode='rb')` in consensus after deSALT exit 139/137/-11/-9 | empty-BAM placeholder is `@HD`-only (no `@SQ`); consensus name-sort fails | **`_create_empty_name_sorted_bam` must include `@SQ` lines** — read `<genome>.fai` and synthesize. See "cDNA align: empty-BAM placeholder missing @SQ" below. |
| `double free or corruption (!prev)` in `rectify correct` on a deSALT BAM (human DRS) | malformed MD tag `…(N0)+^0` from chr4qter / chr7q telomeric N-padded regions; pysam aborts parsing it | Drop or sanitize MD tag from deSALT BAMs before `rectify correct` (or skip reads whose MD ends with `^0`). See "Malformed MD tags from subtelomeric N-runs" below. |

---

## cDNA align: empty-BAM placeholder missing @SQ — bug, not workaround

**Found 2026-05-28** while debugging human GM12878 PCR-cDNA alignment.

When deSALT crashes with one of the recognized exits (`139 / 137 / -11 / -9`),
`run_desalt()` correctly catches it and calls `_create_empty_name_sorted_bam`
(`multi_aligner.py`) to emit a placeholder so the consensus step can proceed
with the remaining aligners. The placeholder, however, contains **only one
header line:**

```
@HD	VN:1.6	SO:queryname
```

No `@SQ` lines for the reference contigs. On yeast that survived by accident
(`samtools sort -n` on the placeholder happened to succeed in the historical
path), but in `rectify align`'s current consensus pipeline the placeholder
flows into `consensus._ensure_name_sorted` → `samtools sort -n` → reopen with
pysam → **`ValueError: file has no sequences defined (mode='rb') — is it
SAM/BAM format? Consider opening with check_sq=False`**.

Reproduced 2026-05-28 on cDNA chunk 2 (one of 8 chunks): deSALT exited −11 in
the second-pass loop, the placeholder was emitted, consensus crashed.

**Fix (applied 2026-05-29):** `_create_empty_name_sorted_bam` now accepts an
optional `genome_path` argument; both `run_desalt()` call sites pass it. When
provided, the function reads `<genome_path>.fai` and synthesises
`@SQ\tSN:<seq>\tLN:<len>` lines into the placeholder header. The original
`@HD`-only behaviour is preserved when `genome_path` is omitted (backwards
compatible — no current caller omits it, but the parameter is `Optional` to
keep test/utility paths working). Implementation:

```python
def _create_empty_name_sorted_bam(output_bam: Path,
                                   genome_path: Optional[str] = None) -> None:
    header = [b'@HD\tVN:1.6\tSO:queryname']
    if genome_path:
        fai = Path(str(genome_path) + '.fai')
        if fai.exists():
            with open(fai) as fh:
                for line in fh:
                    parts = line.rstrip('\n').split('\t', 2)
                    if len(parts) >= 2:
                        header.append(f'@SQ\tSN:{parts[0]}\tLN:{parts[1]}'.encode())
    payload = b'\n'.join(header) + b'\n'
    result = subprocess.run(['samtools', 'view', '-bS', '-o', str(output_bam)],
                            input=payload, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(...)
```

This makes the placeholder a valid header-only BAM that pysam will accept with
default `check_sq=True`. The change is backwards-compatible: callers that don't
pass `genome_path` (none currently exist in the codebase, but future ones)
still get the historical `@HD`-only placeholder.

---

## Malformed MD tags from subtelomeric N-runs — `rectify correct` aborts

**Found 2026-05-28** while debugging human GM12878 DRS junction_overhang
calibration.

### Symptom

`rectify correct -j N` on a deSALT human DRS BAM aborts with
`*** Error in 'python': double free or corruption (!prev): 0x... ***` (SIGABRT,
exit 134). Symptom present at `-j 4` and `-j 1` — **not a multiprocessing or
spawn-pool issue.** Pure data-dependence in `rectify correct`'s pysam-using
hot path.

In the production junction_overhang run on whole-genome GM12878 IVT, 2 of 22
per-chrom tasks failed this way (chr4, chr7). chr8/15/16 and the other 17
chroms completed cleanly with the same code path and the same BAM source.

### Forensic bisection (2026-05-28, this session)

The full chr4 deSALT BAM (638,482 reads) was bisected by coordinate down to a
single triggering read. Every level halved/quartered the surviving range until
one read remained:

| Step | Range | Reads | Result |
|------|-------|-------|--------|
| Halves | chr4:0–95M | 317,444 | OK (11 min) |
| | chr4:95M–190M | 320k | CRASH |
| Quarters | chr4:95M–119M (q1) | 69k | OK |
| | chr4:119M–143M (q2) | clean | OK |
| | chr4:143M–167M (q3) | 69k | OK |
| | **chr4:167M–190M (q4)** | 70,088 | **CRASH @ 73 s** |
| q4 quarters | chr4:167M–172.8M (q4a) | 17k | OK |
| | chr4:172.8M–178.6M (q4b) | 16k | OK |
| | chr4:178.6M–184.4M (q4c) | 17k | OK |
| | **chr4:184.4M–190.214M (q4d)** | 19,320 | **CRASH @ 74 s** |
| q4d quarters | q4da/b/c | OK | |
| | **chr4:188.75M–190.214M (q4dd)** | 7,377 | **CRASH @ 18 s** |
| q4dd quarters | dd1/dd2/dd3 | OK | |
| | **chr4:189.848M–190.214M (q4dd4)** | 3,124 | **CRASH @ 26 s** |
| 4 × 91 kb | dd4a/b/c | OK | |
| | **chr4:190.1225M–190.2146M (dd4d)** | 909 | **CRASH** |
| 4 × 23 kb | dd4d1/dd4d2/dd4d3 | OK | |
| | **chr4:190.1915M–190.2146M (dd4d4)** | 145 | **CRASH** |
| 4 × 6 kb | dd4d4a/dd4d4b | OK | |
| | **chr4:190.203M–190.2088M (dd4d4c)** | 61 | CRASH |
| | **chr4:190.2088M–190.2146M (dd4d4d)** | **1** | **CRASH (single read!)** |

The 1-read BAM was the smoking gun.

### The trigger read

```
QNAME:  6bf557b1-019f-4432-be1f-a7498a493bab
RNAME:  chr4    POS: 190,203,757    MAPQ: 60    FLAG: 16 (reverse, primary)
SEQLEN: 2,822 bp
CIGAR:  1S 145M 31D 26M 5D 289M 1D 4M 1D 278M 1I 6M 1I 74M
        59453N 3D 21M 12D 21M 45N 1D 40M 1D 7M 5I 28M 91I 1M 50N
        10M 12D 22M 1751S
MD:     ...A103(N0){65}^0
```

The MD ends with **65 consecutive `N0` markers** (each = "an `N` ref base, then
0 matching bases") followed by **`^0`** — a deletion-marker with no base
letters (the SAM spec requires `^` to be followed by ACGT/N, not a digit). The
deSALT MD-emitter walks off the end of the genomic reference into the
chr4qter N-padded subtelomeric block and serializes the run incorrectly.

### Confirmation

Stripping the MD tag from this one read produces a BAM that `rectify correct
-j 1` processes cleanly in 1 second. Re-introducing the MD tag re-triggers the
crash deterministically. The CIGAR (with its 59-kb N op and 1,751-bp soft
clip) is **not** the trigger — the no-MD smoke succeeded with the same CIGAR.

### Hypothesis (not yet AddressSanitized)

Pysam's MD-tag parsing path (likely `get_aligned_pairs(with_seq=True)` or
`get_reference_sequence()`) allocates a buffer sized from the cumulative `N`
runs, then either over-writes when consuming the malformed `^0` (treating the
digit `0` as both the deletion-length and the start of the next op), or
double-frees the buffer when retrying. The crash is in glibc malloc, so
`fsanitize=address` on a pysam build should pinpoint the line within a single
run.

### Mitigations (pick one)

1. **Strip MD tags from deSALT BAMs before `rectify correct`.** Cheapest. The
   `rectify correct` walkback algorithm does not require MD — it reads
   sequence and reference independently. Drop with:
   ```bash
   samtools view -h foo.deSALT.bam | awk 'BEGIN{OFS="\t"} \
     /^@/{print; next} \
     {out=$1; for(i=2;i<=11;i++) out=out OFS $i; \
      for(i=12;i<=NF;i++) if (substr($i,1,5)!="MD:Z:") out=out OFS $i; \
      print out}' | samtools view -bS -o foo.deSALT.noMD.bam -
   ```
2. **Skip the bad reads at ingest time.** Filter records whose MD tag matches
   `/(N0){5,}\^0/` (or any `\^\d`) before passing to pysam.
3. **Recompute the MD with `samtools calmd`.** Produces a correct MD from the
   reference — but doesn't help if the reference itself contains the N-padded
   block (calmd will re-emit `N0…` markers, just without the trailing `^0`).
   Tested: untested in this session; expected to help only for the `^0` part.
4. **Patch pysam upstream.** Strongest, slowest. Requires building pysam with
   AddressSanitizer to confirm the over-write location.

In production we worked around the bug by dropping chr4 + chr7 deSALT
contributions from the WG junction_overhang calibration (calibrate script
already skips missing per-aligner TSVs; 20 of 22 chroms contributed 3-aligner
unanimous concordance, the other 2 contributed minimap2 + uLTRA).

### Files on Sherlock (preserved 2026-05-28)

| File | Path |
|------|------|
| 1-read trigger BAM | `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/desalt_smoke/chr4_dd4d4d.bam` |
| Full chr4 deSALT BAM | `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/junction_overhang/perchrom/chr4/chr4.deSALT.bam` |
| chr7 deSALT BAM (same failure mode) | `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/junction_overhang/perchrom/chr7/chr7.deSALT.bam` |
| Smoke crash log (-j 1) | `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/desalt_smoke/chr4_j1_26621424.out` |
| MD-stripped smoke (passes) | `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/desalt_smoke/chr4_dd4d4d_noMD_corrected.tsv` |
| Bisect sbatch scripts | `/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/desalt_smoke/chr4_*.sbatch` |

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
