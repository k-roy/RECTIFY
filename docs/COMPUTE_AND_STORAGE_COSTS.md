# Compute and Storage Costs

**What a RECTIFY run actually costs, measured on real production libraries — so you can choose an
aligner panel on evidence instead of defaults.**

> ## 🔴 The headline, because it is counter-intuitive
>
> **`rectify correct`, not alignment, is the pipeline. And `correct` runs once per aligner.**
>
> On a real production DRS library (1,264,614 reads, 2-aligner panel), the **entire alignment stage
> cost 2.50 CPU-hours** while **`correct`'s minimap2 arm alone cost 217.9 CPU-hours** — `correct` is
> **87×** alignment (195× counting CPU actually consumed including retries). **Alignment is 0.25–0.54%
> of the pipeline.**
>
> ⇒ **Dropping an aligner saves you its `correct` arm, not its alignment.** Choosing a panel by
> "how fast does this aligner run" optimises less than 1% of your bill. The relevant question is
> **how expensive is that aligner's BAM to correct** — and those costs differ by **7×** (§2b).

---

## 0. Scope — three different measurements, do not cross-quote them

| tag | what | scale |
|---|---|---|
| **PROD** | SCG production DRS `ysh1_rep1`, **2-aligner** (minimap2 + uLTRA), 32 chunks, `sacct` ledger | **1,264,614 reads** |
| **BENCH-DRS** | H2 DRS `mex67aa_rep1`, per-aligner instrumented, 4 threads | 10,000 reads |
| **BENCH-cDNA** | H2 UMI-collapsed ONT PCR-cDNA, 8 threads, **deSALT absent** | 33,923 reads |

Absolute throughput is protocol- and cluster-specific; the **ratios** are what transfer. Every
number below is read back from a log, `sacct`, `du`, or `samtools`, not estimated. Measured
2026-05→07, tabulated 2026-08-06.

⚠️ **`run-all` and `align` have different defaults.** `run-all` long-read = **5 aligners**
(minimap2 + mapPacBio + gapmm2 base, uLTRA + deSALT junction) **+ chimeric consensus**
(`run_command.py` 685–737). `rectify align` = **3** base aligners, **no** junction aligners,
chimeric **off** (`align_command.py` 126–171). A benchmark built on `align` does not exercise the
`run-all` panel.

⚠️ **And production used neither default** — all 35 samples of the SCG DRS fan-out ran
`ALIGNERS=("minimap2" "uLTRA")`. Three configurations are in play and the honest answer differs for
each.

---

## 1. Where the CPU actually goes — PROD, `sacct` TotalCPU

| stage | tasks | CPU-h (COMPLETED only) | CPU-h (consumed, incl. retries) |
|---|---:|---:|---:|
| `align` (minimap2 + uLTRA) | 64 | **2.50** | 2.50 |
| `merge_aln` | 1 | 0.051 | 0.051 |
| `prescan` | 1 | 0.030 | 0.030 |
| `chunk_merge` | 28/32 | 2.383 | 2.383 |
| **`correct` — minimap2 arm** | 32 | **217.9** | **486.7** |
| **`correct` — uLTRA arm** | 32 | **240.5** | **517.5** |
| **TOTAL** | | **463.4** | **1,009.2** |

**Alignment is 0.54% of the COMPLETED-only pipeline, 0.25% of consumed.**

⚠️ The correct arrays show 32 COMPLETED + 32 FAILED + 32 TIMEOUT per arm; the TIMEOUT rows came from
an earlier 90-minute submit config. Both ledgers are reported because they bracket the truth.

## 2. Per-aligner cost

### 2a. Alignment — mapPacBio dominates, *not* the junction aligners

BENCH-DRS, 10,000 reads, 4 threads:

| aligner | seconds | share of the projected 5-way align stage |
|---|---:|---:|
| minimap2 | 4.2 | 0.6% |
| **mapPacBio** | **533.2** | **75.2%** |
| gapmm2 | 26.7 | 3.8% |
| uLTRA | 56.3 | 7.9% |
| deSALT | 38.4 | 5.4% |
| chimeric consensus | 47.2 | 6.7% |

**uLTRA + deSALT together are 13.3%.** At PROD scale uLTRA costs **8.22×** minimap2 (4,379.3 s vs
532.5 s over 32 chunks).

**Scaling behaviour — of the ALIGNMENT stage only.** minimap2 runs **2,375 reads/s** (PROD) vs
**2,381 reads/s** (BENCH-DRS) — 0.3% agreement across two clusters and a 126× difference in library
size, so **alignment** extrapolates linearly in read count with confidence. mapPacBio fits
**0.0562 s/read** with intercept ≈ 0, i.e. essentially all marginal cost.

⚠️ **Correction (2026-08-07): mapPacBio DOES thread-scale, ~1.6× from 4→8 threads.** An earlier
version of this page said it did not, on the basis of 18.8 reads/s (BENCH-DRS, 4 thr) vs 18.1
reads/s (BENCH-cDNA, 8 thr). **That comparison was confounded** — different library *and* different
protocol. A same-library measurement shows ~1.6×.

🔴 **Do NOT extend this linearity to `rectify correct`.** Correct is driven by junction-pool density,
not read count — see **§8**. Read-count extrapolation of the *pipeline* is invalid.

### 2b. 🔴 Correction — this is the number that matters, and it varies 7×

Instrumented single-BAM `rectify correct`, BENCH-DRS, 1 thread:

| aligner's BAM | correct cost | relative |
|---|---:|---:|
| **minimap2** | **112.17 s** | **1.00× (cheapest)** |
| mapPacBio | 139.37 s | 1.24× |
| uLTRA | 291.54 s | 2.60× |
| deSALT | 611.73 s | 5.45× |
| **gapmm2** | **784.56 s** | **6.99× (most expensive)** |
| sum of all five | 1,939.37 s | 17.28× |

**minimap2's BAM is the cheapest of all five to correct.** A single-aligner run therefore keeps the
cheapest arm and drops the expensive ones — the saving compounds with §1.

## 3. Coverage — the panel adds NO reads

Read-level (distinct primary QNAMEs, `-F 0x900`), BENCH-cDNA, 33,923-read input:

| aligner | distinct reads | share | **reads unique to it** |
|---|---:|---:|---:|
| minimap2 | 33,923 | **100.00%** | **0** |
| mapPacBio | 33,923 | 100.00% | **0** |
| uLTRA | 33,904 | 99.94% | **0** |
| gapmm2 | 33,526 | 98.83% | **0** |
| **union of all four** | **33,923** | **100.00%** | — |

**The union of all four equals what minimap2 covers alone.** Not one aligner contributes a read the
others miss. **gapmm2's read set is a strict subset of minimap2's** (0 reads in gapmm2 not in
minimap2; 397 the other way) — as expected, since gapmm2 *is* a minimap2 wrapper adding edlib
terminal-exon refinement.

⇒ **A panel buys arbitration between competing placements of the same reads. It does not buy
coverage.**

> **⚠️ Coverage is a FLOOR, not the question.** That every aligner reports an alignment for nearly
> every read is close to trivial — the aligners differ in the **CIGAR** they produce, and that is
> what determines where 3′ ends, splice junctions and 5′ ends actually land. Nothing in this table
> tells you whether a reduced panel places reads *correctly*; it only tells you it does not lose
> them. **Do not read §3 as licence to reduce the panel.** The placement comparison — per-read CIGAR
> agreement, junction recovery and end-position deltas on corrected output at
> alternative-splicing-rich loci — is the test that settles it, and it is not yet done.

> ### 🔴 Never read the consensus log's `By aligner combo:` line as coverage
>
> It is a **per-record positional** census (40,448 records for 33,923 reads). A key like
> `'gapmm2': 6540` means *"gapmm2 placed those reads differently"*, **not** *"gapmm2 aligned reads
> no one else could"*. Misread as coverage it implies minimap2 misses ~16% of the data when minimap2
> in fact aligns **100%**. Use the read-level table above.
>
> Note also that `by_aligner_combo` is written by two code paths with **different semantics** —
> `consensus.py:805` = *aligners compared*; `consensus.py:778` = *segment winners* — into the same
> key. The table above is from a run with chimeric consensus **off**.

## 4. Storage

**PROD, measured `du -sh`** — whole `chunks/` tree = **9.5 GB for 1,264,614 reads = 7.51 KB/read**
(minimap2 2.1 GB · uLTRA 2.0 GB · merged_bams 1.5 GB · consensus 1.9 GB · input chunks 1.5 GB ·
`junction_pool.pkl` 12.4 MB · `rescue_scan.pkl` 0.51 MB).

**Against the raw BAM: 11.3×** for the 2-aligner config. ⚠️ **An in-house figure of "4–6× the input
BAM" circulates — it understates the real footprint by about 2×.** The `run-all` 5-aligner default
projects to order **~14.5 KB/read ≈ 22× raw BAM** *(projected, not measured)*.

**Two pure-overhead artifacts the panel creates:**

1. **A `namesorted` duplicate roughly doubles every aligner's BAM.** BENCH-DRS: minimap2 bam
   8,600,322 B + namesorted 8,692,894 B.
2. **`*.mpb_san.fastq` — uncompressed, and never deleted.** mapPacBio requires a QNAME-sanitised
   FASTQ (`mapPacBio.sh` copies the whole FASTQ header into the SAM QNAME, breaking the 254-char
   limit). It is created at `rectify/core/align/multi_aligner.py:1027` and — unlike every other
   intermediate in that module (`chunk_tmp_fq`, `mpb_split_fq`, `sam_path`, `mpb_pre_stitch`, all
   `unlink`ed) — **never removed**. In BENCH-cDNA it was **79,841,821 of 144,906,349 B = 55.1% of the
   entire panel output directory**.
   **If you run mapPacBio, delete `*.mpb_san.fastq` after alignment succeeds.** It is regenerable and
   nothing downstream reads it.

## 5. What a single-aligner run saves

Measured on PROD (`ysh1_rep1`, 1,264,614 reads):

| config | align | correct | other | **total CPU-h** |
|---|---:|---:|---:|---:|
| production as run (minimap2 + uLTRA) | 2.50 | 458.4 | 2.46 | **463.4** |
| **minimap2 only** | **0.27** | **217.9** | 2.41 | **220.6** |
| **saving** | 2.23 | **240.5** | 0.05 | **242.8 = 52.4%** |

**~99.5% of that saving is the second `correct` arm, not the alignment.**

| baseline | CPU saving | disk saving |
|---|---|---|
| production 2-aligner (minimap2 + uLTRA) | **~52%** *(measured)* | **~3.3×** *(measured)* |
| `run-all` 5-aligner default | **≥80%** *(floor; point estimate ~94% is **projected**)* | **~6.4×** *(projected)* |

The ≥80% floor is the defensible number: the 5-arm projection assumes correct cost is additive
across arms and uses relative BENCH-DRS ratios. The measured saving is also **biased downward** —
the run's `consensus` and `final_merge` stages never released (upstream failures), and a chimeric
consensus over one aligner is a no-op, so those would be pure additional single-aligner savings that
this ledger does not credit.

## 6. Practical guidance

- **Optimising alignment is optimising <1% of your bill.** Ask what an aligner's BAM costs to
  *correct* (§2b), not how fast it aligns.
- **Storage- or throughput-bound? Drop mapPacBio first** — 75.2% of the align stage, zero unique
  reads, and the source of the largest artifact on disk (§4).
- **Terminal-end work (poly(A) sites, APA, internal polyadenylation, TSS): minimap2 alone is
  defensible.** It aligns 100% of reads, its BAM is the cheapest of all five to correct, and the
  terminal-exon refinement gapmm2 adds at the aligner stage is the same class of correction
  `rectify correct` applies downstream via `splice/terminal_exon_refiner.py` and
  `splice_aware_5prime.py`. **gapmm2 is the single most expensive BAM to correct (6.99×)** — paying
  that to redo work the correction stage does anyway is hard to justify when ends are the
  deliverable. *(Not yet quantified: whether the two refinements agree read-for-read.)*
- **Novel-junction discovery** (cryptic splicing, unannotated introns)? **Keep the full panel.** The
  junction pool admits an observed junction by anchored support **or** cross-family concordance, and
  the concordance channel needs ≥2 independent algorithm families — **with one aligner it is
  structurally inert.** See [Multi-Aligner Consensus](algorithms/multi_aligner_consensus.md).
- **Always `du -sh` the output after your first sample**, before fanning out.
- 🔴 **Do not size a big run from a small-region `correct` benchmark** — cost appears superlinear in
  local junction-pool density (§8), so a small window under-states it badly.

## 7. What is NOT measured — read before extrapolating

1. **deSALT, mapPacBio and gapmm2 at production scale.** The SCG panel ran minimap2 + uLTRA only.
   Only 10 k-read points exist for the other three.
2. **A clean whole-pipeline `run-all` wall time.** The only in-corpus DRS `run-all` with `[TIMING]`
   prints **failed at the correction step**; its alignment numbers are sound, its total is not
   quotable.
3. ~~The `analyze` stage — never measured.~~ ✅ **MEASURED 2026-08-10 (planning/657, tree
   `999ceb5`, 9-library AA DRS panel, manifest mode, `--run-deseq2`): full-genome analyze =
   2,862 s wall ≈ 0.8 CPU-h ONCE per panel (~0.03% of the per-sample bill).** Per-stage: CPA
   clustering `[2/9]` **47.3%** (the only SUPERLINEAR stage: 52× wall at 17.5× clusters — the
   watch-item if panels get much deeper), genomic distribution 22.5%, count matrix 13.2%,
   DESeq2 3.9%. True max RSS **5.86 GB** — the 34.5 GB "maxvmem" seen in qacct is address
   space, not working set. cpu/wall ≈ 1.18 despite `--threads 8`: effectively serial, do not
   request 8 slots for it.
   **Consensus read selection is also now measured** (same unit, 31.9 k-read window, threads 1):
   raw `rectify consensus` 21–34 s at 2→5 aligners (sublinear in arm count, 1.5–1.6×),
   `--chimeric` **2.7×** dearer; corrected consensus (`merge_corrected_tsvs` +
   `write_corrected_consensus_bam`) 37→59.5 s at 2→5 arms (merge:write ≈ 3:1). Linear per-read
   projection to a 1.26 M-read PROD sample: raw ~22 min, corrected ~39 min, chimeric ~59 min,
   single-threaded ⇒ **~0.1% of the post-gate per-sample bill. Neither consensus nor analyze is
   a lever; `correct` remains ~99% even after the planning/649 pool gate.**
4. **Whether the panel's winning placements beat minimap2's** on reads all aligners placed. That
   needs a concordance comparison of *corrected* output.
5. 🔴 **A 12–90× discrepancy in per-read correct cost**: 0.95–2.86 CPU-s/read at PROD versus
   0.011–0.078 CPU-s/read in BENCH-DRS. **Partially accounted 2026-08-11 (planning/666,
   3-scale × 3-pool grid on the dense windows): startup amortisation and local pool density
   are now excluded — MARGINAL bench cost is 0.018 (gated pool) / 0.042 (ungated) CPU-s/read,
   leaving a REAL 23–68× residual vs PROD.** Leading suspects, testable in order: the
   production pool is panel-wide and genome-wide (~2 orders more junctions in the index),
   per-aligner arms with higher rescue rates (uLTRA), hardware. The closing experiment:
   `correct` the same window reads against the actual production `junction_pool.pkl`.
   **Until then, still plan with the PROD figure.**
   Three quotation rules from the same unit: (i) **quote MARGINAL per-read cost, never
   total/reads** — fixed startup is 10–28 s, so any per-read figure from a ≤500-read run is
   startup-dominated; (ii) `--threads 4` buys only **1.27–1.51×** (32–38% efficiency) —
   CPU-h is nearly threading-invariant, don't request wide allocations for `correct`;
   (iii) threaded output differs from serial on **~0.9% of reads** on pre-`305daff` trees
   (tie-break nondeterminism) — byte-identity gates must run single-threaded or on a tree
   with the coordinate-ordering fix.

## 8. ⚠️ `correct` cost scales with local junction-pool DENSITY, not with read count

**Status: the MECHANISM is verified (profile + code). The SCALING EXPONENT is NOT established.**
Stated here because sizing a large run on read count alone is the expensive mistake.

Same reads, same BAM, **same number of aligner BAMs** — only the junction pool differs. Two
measurements, and they do **not** agree on the exponent:

| observation | pool-density ratio | cost ratio |
|---|---:|---:|
| 3-region pass, run-all shape (partial — arm B never completed) | ~7.2× | **≥26×** |
| PTC7, full A/B under the single-`correct` design | **4.9×** | **3.7×** |

⚠️ **An earlier version of this page called this "roughly quadratic". That was over-claimed** — it
rested on the first row, where the slower arm had not finished (so its cost was a lower bound) and
was additionally inflated by a pipeline shape that ran `correct` once per aligner. The completed
head-to-head gives **roughly linear** in pool density. **Treat the exponent as unknown between
~linear and worse**, and do not extrapolate either way.

What *is* solid: **cost tracks pool density rather than reads**, for the reason profiled below.

**Root cause — profiled, not inferred.** `py-spy` on **4 hot workers across 2 jobs and 2 compute
nodes** landed in the same frame every time:

```
_hp_edit_distance             splice_aware_5prime.py:742/744/746/749
_rescue_3ss_truncation_body   splice_aware_5prime.py:1542, 2033
rescue_3ss_truncation         splice_aware_5prime.py:1147, 1155
correct_read_3prime           bam_processor.py:491
```

Module 2F's 3′SS rescue runs `for _shift … × for _off … × _hp_edit_distance` — an O(n·m) DP —
**once per candidate junction**. `_RESCUE_DP_CAP` (=100) bounds the DP *sequence length*;
**nothing bounds candidate multiplicity.** So cost grows with the number of candidates in the local
window, which is what the pool density measures.

An I/O confound was excluded by measurement, not argument: Δcpu/Δwall over 200 s = **2.00 cores** for
both jobs with `ioops` essentially frozen (+29, +18); `ps` showed 2 workers at 99% and 6 idle.

> ### 🐛 Related defect
> **`--junction-max-candidates-per-nop` does not gate the 3′SS rescue path.** It is plumbed from
> `correct_command.py:357` into `refine_bam_junctions` (`:764`) only; `rescue_3ss_truncation` in
> `splice_aware_5prime.py` takes no candidate cap (verified: no `max_candidates` symbol in that
> module). So the flag a user would reach for to bound this cost **does not bound it.**

**Consequences:**
- **A benchmark run over a small region under-states `correct` cost dramatically**, because its pool
  is small. This is the leading explanation for the 12–90× gap in §7 item 5.
- **Cost does not scale with read count alone.** A denser annotation, a repeat-rich locus, or a
  deeper library that surfaces more novel junctions all raise cost superlinearly.
- **A reduced panel may save more than linear accounting predicts**, since a single-aligner pool is
  smaller — but that is an *inference from the same unproven scaling*, so do not bank it.

**If you are sizing a large run**: measure `correct` on a region whose pool density matches
production, not on a convenient small window.

## 9. See also

- [Aligner Recommendations](ALIGNER_RECOMMENDATIONS.md) — which aligners to choose, and why
- [Multi-Aligner Consensus](algorithms/multi_aligner_consensus.md) — consensus and the junction pool
- [HPC / SLURM](user_guide/hpc_slurm.md) — chunking and resource requests
