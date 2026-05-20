# Commit Zero — Profile Results

**Date:** 2026-05-19
**Git SHA profiled:** `7dbb1bd fix(ultra): ship SGD native GTF + fix normalizer for gffread-produced GTFs`
**Branch:** `drs-validation-rebuild`
**Cluster:** Hoffman2 login node (compute-quality, profiles ran clean)
**Profiling tool:** py-spy 0.4.2, `--rate 250 --native --idle`, `--format raw` folded stacks
**Phase Zero-A:** ✅ DONE. Phase Zero-B (full Han wt_R1 stage timing on cluster): pending.

---

## DECISION (top of doc, governs all downstream commits)

**PROCEED with Commits A–F as written, with one protocol-conditional refinement: the success metrics + the relative weight of each commit's payoff need to be split between QSrev/short-read and DRS/long-read.** The architecture (region-parallel write + lazy TSV + concurrent aligners + analyze partial-streaming + provenance/resume) is correct for both arms, but the magnitude of each commit's win differs sharply by protocol.

**Refuted hypotheses:**
- The advisor's BGZF-compression hypothesis (`compresslevel=1` one-liner): **NO.** BGZF is 1.9% of WRITE phase for QSrev and 23.7% of (an already-tiny) WRITE phase for DRS. The cheap-win commit (Zero-bis) is **not needed**; do not add it.
- The "write_corrected_bam is THE bottleneck universally" framing from the original issue doc: **partially refuted.** True only for short-read protocols. For DRS/long-read, write_corrected_bam is a rounding error (0.5% of wall).
- The "batched set_tag" hypothesis: **NO.** SET_TAG is 0.2% (QSrev) and 2.0% (DRS) of WRITE phase. Not worth optimizing.

**Confirmed hypothesis:** `realign_exon_blocks` is the single dominant cost across BOTH protocols, but it runs in DIFFERENT phases — within `write_corrected_bam` for QSrev (88% of WRITE), within the correct-region loop for DRS (74% of CORRECT). Both phases get parallelized by the planned architecture (Commits A–B parallelize WRITE; the existing `process_bam_file_parallel` already parallelizes CORRECT).

**Per-protocol expected speedup (at 16 cores, scaled from this 1-core profile):**

| Sample | Today @ -j 16 | After Commits A-B | Speedup |
|---|---|---|---|
| QSrev (200k slice) | ~6.3 min (write dominates at 90%) | ~1 min (region-parallel write) | **6×** total, **~16× on write phase** |
| DRS (20k slice) | ~2.9 min (correct still 81%, write 6%) | ~2.6 min (only ~10% saved on already-small write) | ~1.1× total — marginal |

**Implication for the master spec:** §0 success metric (Han wt_R1 from 5+ hr → <75 min) is a QSrev story. For DRS, the bigger lever is reducing per-read realign cost itself (caching results across phases, or algorithmic optimization) — **outside the scope of the parallelism refactor.** I'll note this in §9 of the master plan as a follow-up.

---

## 1. Inputs profiled

| Sample | Protocol | Reads | Read length (p10/median/p90/max) | Wall @ -j 1 (with py-spy native overhead) |
|---|---|---|---|---|
| `wt_R1_slice_200k.bam` from Han 2023 multialigner v0.9.1 | QuantSeq REV (`--short-read --dT-primed-cDNA`) | 207,635 | ~50 bp short reads | 8m 42s |
| `drs_wt_rep1_20k.bam` from `nmd_as_audit/subsets/` (24% subsample) | ONT DRS (default flags) | 20,372 | 364 / 872 / 3195 / 4048 bp | 37m 54s |

Note: DRS at -j 1 takes far longer than QSrev despite being 10× fewer reads — `realign_exon_blocks` cost scales roughly quadratically in read length, and DRS reads are 17× longer at median. py-spy `--native` adds ~30-50% sampling overhead but is required to attribute pysam-C / BGZF / Numba-runtime frames.

## 2. Sample counts captured

| Sample | Total samples | Unique stacks | Wall (s) | Samples/sec (effective) |
|---|---|---|---|---|
| QSrev 200k | 130,597 | 8,849 | 522 | 250 |
| DRS 20k | 569,370 | 7,541 | 2,274 | 250 |

Sampling rate held at 250 Hz throughout; sample count ratio ≈ wall time ratio. py-spy reported 1 sampling error per run (probably a transient ptrace race) — negligible.

## 3. Phase breakdown

Phases are determined by the call-stack ancestors. Definitions:

- `WRITE` = stack contains `write_corrected_bam` or `bam_writer.py`.
- `CORRECT_REGION` = stack contains `_process_region_worker` or `process_bam_file_parallel` or `parallel.py`, AND NOT `bam_writer.py`.
- `VARIANT_SCAN` = stack contains `run_variant_aware_scan` or `variant_scan.py`.
- `STARTUP` = importlib-dominated stacks (genome load + module imports).
- `OTHER` = everything else.

| Phase | QSrev samples | QSrev % | DRS samples | DRS % |
|---|---|---|---|---|
| **WRITE** | **85,273** | **65.3%** | 3,059 | 0.5% |
| **CORRECT_REGION** | 37,737 | 28.9% | **560,125** | **98.4%** |
| STARTUP | 3,566 | 2.7% | 3,969 | 0.7% |
| OTHER | 2,486 | 1.9% | 1,698 | 0.3% |
| VARIANT_SCAN | 1,535 | 1.2% | 519 | 0.1% |

## 4. Bucket breakdown within WRITE phase

Buckets defined by regex match on ancestor frames (leaf-to-root, first match wins). Patterns documented at the top of `/tmp/parse_folded.py`.

| Bucket | QSrev % of WRITE | QSrev % of total | DRS % of WRITE | DRS % of total |
|---|---|---|---|---|
| REALIGN (realign_exon_blocks + align_exon_block_global) | **88.0%** | **57.4%** | 5.6% | 0.0% |
| OTHER (mostly `_decode_eq_seq_inplace` + `reroute_intronic_tail_5prime_via_junction`) | 5.5% | 3.6% | 15.7% | 0.1% |
| PYSAM_AC (pysam C-level accessors: cigartuples, query_sequence, etc.) | 3.2% | 2.1% | 40.2% | 0.2% |
| BGZF (`bgzf_*`, `libdeflate`, `htslib`) | 1.9% | 1.2% | 23.7% | 0.1% |
| PYSAM_IO | 0.6% | 0.4% | 2.5% | 0.0% |
| CIGAR_OPS (clip_read_to_corrected_3prime, _hardclip_trailing_a_run) | 0.3% | 0.2% | 4.0% | 0.0% |
| SET_TAG | 0.2% | 0.1% | 2.0% | 0.0% |
| NUMPY (array_subscr/array_new — Numba runtime) | 0.2% | 0.1% | 0.7% | 0.0% |
| WALKBACK | 0.0% | 0.0% | 5.6% | 0.0% |

**Key takeaway:** within the WRITE phase, QSrev work is overwhelmingly realign-driven (per-read CIGAR rewrite re-runs the HP-penalty global alignment). DRS WRITE is dominated by pysam accessors + BGZF, because by the time we're in write, the per-read mutations for DRS are mostly already done in correct (or are no-ops, given DRS reads have less to fix per read post-correction).

## 5. Bucket breakdown within CORRECT_REGION phase

| Bucket | QSrev % of CORRECT | QSrev % of total | DRS % of CORRECT | DRS % of total |
|---|---|---|---|---|
| **REALIGN** | 12.9% | 3.7% | **73.8%** | **72.6%** |
| **NUMPY** (Numba runtime: array_ass_subscr, array_subscr) | 6.4% | 1.9% | **17.9%** | **17.6%** |
| WALKBACK (walkback_3prime_guarded, rescue_3ss_truncation) | 37.1% | 10.7% | 7.4% | 7.3% |
| OTHER | 24.5% | 7.1% | 0.3% | 0.3% |
| PYSAM_AC | 14.1% | 4.1% | 0.3% | 0.3% |
| BGZF | 2.7% | 0.8% | 0.2% | 0.2% |
| PYSAM_IO | 2.3% | 0.7% | 0.0% | 0.0% |

For DRS, **REALIGN + NUMPY together = 91.7% of CORRECT phase = 90.2% of total wall.** Both are part of the Numba-compiled global alignment kernel. The single biggest lever for DRS speedup is reducing the per-read realign cost — either by caching results from prior phases or by algorithmic improvements (k-mer prefilter, banded DP at smaller bandwidth, etc.). **None of these are in scope for this refactor.**

For QSrev, the correct phase is dominated by WALKBACK (37%) — that's `walkback_3prime_guarded` and `rescue_3ss_truncation`. But correct is only 28.9% of QSrev wall, so this is 10.7% of total. Already region-parallelized via `process_bam_file_parallel`.

## 6. Top frames (leaf-to-root, ancestor-attribution; first match wins)

### QSrev top 10 within WRITE phase
1. `align_exon_block_global` (rectify/core/align/local_aligner.py, multiple lines around :600-:650) — ~700+ samples each line
2. `realign_exon_blocks` (rectify/core/bam/read_edits.py:700) — kernel entry
3. `_decode_eq_seq_inplace` (bam_writer.py:86) — `=` SEQ expansion, ~12 samples
4. Numba-internal allocator + libm `finite` — Numba runtime in realign kernel

### DRS top 10 within CORRECT_REGION phase
1. `realign_exon_blocks` (read_edits.py:700) — kernel entry, ~80%+ of CORRECT
2. `array_ass_subscr` (arraymodule.c) — Numba array writes within realign
3. `_int_malloc` (libc) — heap allocs during alignment DP
4. `_rescue_3ss_truncation_body` (splice_aware_5prime.py) — small contribution
5. `walkback_3prime_guarded` (walkback.py:498)

## 7. Implications for the master spec (`parallel_refactor_plan.md`)

Changes to roll into the spec immediately:

1. **§0 Goal — broaden the success metric.** Add a QSrev-vs-DRS split. QSrev: end-to-end < 75 min at 16 cores. DRS: existing correct-phase parallelism already does the bulk; the parallelism refactor adds Axis B (concurrent aligners) and resume, not a dramatic wall-time win.
2. **§2 Pipeline-stage audit — annotate rows with protocol applicability.** Mark rows 5–6 (`process_bam_file_parallel` + `write_corrected_bam`) as **dominant cost for QSrev/cDNA-short**. Mark CORRECT row as **dominant cost for DRS/cDNA-long**.
3. **§3 Commit Zero — mark DONE, link to this doc, drop the `compresslevel=1` cheap-win and the batched-set_tag escape hatches** (both refuted).
4. **§9 Risks — add a new risk:** "DRS expected speedup is small because the bottleneck is INSIDE the already-parallel correct loop. The parallelism refactor delivers the architecture (resume, sidecars, manifest TSVs, analyze partial-streaming) but does not materially change DRS wall time. Communicate this in CHANGELOG so users don't expect a 6× win for DRS like they will see for QSrev."
5. **§NEW Future work section — algorithmic redundancy elimination.** Caching realign output across the correct-phase TSV → write-phase BAM boundary would collapse QSrev's WRITE-phase realign (88% → ~0%). That alone gets QSrev another ~6× on top of the parallelism refactor. **OUT OF SCOPE for this PR.**

## 8. Open observations / surprises

- **The original issue doc's framing was QSrev-specific.** The "2.5 hr in write_corrected_bam" observation on Han wt_R1 (job 13446575) was a QuantSeq REV / short-read run. For DRS we won't see the same wall-time savings from write parallelization. This is fine — the architecture serves both arms — but the success-criteria messaging needs to be split.
- **Realign is being done twice** for QSrev: once in CORRECT_REGION (12.9% of correct = 3.7% of total wall) and once in WRITE (88% of write = 57% of total wall). The CIGAR computed in CORRECT is then re-recomputed in WRITE. Skipping the second pass would save QSrev ~57% of wall time on top of any parallelism. Worth investigating as a follow-up.
- **py-spy NUMPY bucket for DRS at 17.6% of total** is mostly Numba runtime (array_ass_subscr / array_new / array_dealloc inside compiled realign code). Numba's machine code is not symbolizable by py-spy without explicit symbol tables; this fraction is the "visible" portion of unsymbolized Numba work. The true REALIGN attribution is therefore a lower bound; including NUMPY, realign-rooted work is ≥91.7% of CORRECT phase for DRS.
- **DRS WRITE-phase BGZF at 23.7%** would be a meaningful target — but WRITE is only 0.5% of DRS total wall, so 23.7% × 0.5% = 0.1% absolute. Not a useful lever.

## 9. Methodology notes

- py-spy `--idle` captures stack samples even when the GIL is released or threads are sleeping. For these single-threaded -j 1 runs, idle samples should be near zero (the bottleneck is CPU-bound). They are.
- Bucket attribution walks the stack leaf-to-root and returns on the first regex match. Reason: leaf is what CPU is actually doing; ancestors give it semantic context. A `libc/_int_malloc` leaf with `align_exon_block_global` ancestor correctly attributes to REALIGN bucket.
- All py-spy data, the input BAMs, the output corrected BAMs, and the parse script live on H2 under `/u/scratch/k/kevinroy/profile_zero_a/` and `/u/scratch/k/kevinroy/profile_zero_a_drs/`. Parse script also copied to `/tmp/parse_folded.py` on M1.

## 10. Artifacts produced (kept on H2 for re-analysis)

| Path | Description | Size |
|---|---|---|
| `/u/scratch/k/kevinroy/profile_zero_a/wt_R1_slice_200k.bam` | QSrev input slice | 10.6 MB |
| `/u/scratch/k/kevinroy/profile_zero_a/wt_R1_slice_200k.corrected.bam` | QSrev write_corrected_bam output | 11.8 MB |
| `/u/scratch/k/kevinroy/profile_zero_a/folded_qsrev.txt` | QSrev folded stacks (130k samples) | 16 MB |
| `/u/scratch/k/kevinroy/profile_zero_a/flame_qsrev.svg` | QSrev flame graph (older, no --native) | 150 KB |
| `/u/scratch/k/kevinroy/profile_zero_a_drs/drs_wt_rep1_20k.bam` | DRS input slice | 11.3 MB |
| `/u/scratch/k/kevinroy/profile_zero_a_drs/drs_wt_rep1_20k.corrected.bam` | DRS write output | 11.1 MB |
| `/u/scratch/k/kevinroy/profile_zero_a_drs/folded_drs.txt` | DRS folded stacks (569k samples) | 16 MB |
| `/tmp/folded_qsrev.txt`, `/tmp/folded_drs.txt`, `/tmp/parse_folded.py` (M1) | Local copies for re-analysis | — |

Do not delete.

---

## Next step

Phase Zero-B (full Han wt_R1 stage timing on a 16-core qsub job, ~8 hr wall budget). The slice profile + projections above already make a confident case for the architecture decision, but Zero-B confirms:

1. The full-sample wall-time split between correct/write at -j 16 (validating the projection above).
2. Per-stage `[TIMING]` log granularity — if the existing log doesn't break out write-bam vs sort vs index, Commit A's tracker should add finer instrumentation.
3. Any unexpected stage (e.g., merge_corrected_tsvs at scale) that emerges as a next-bottleneck candidate after write parallelization.

The decision in §0 above does NOT depend on Zero-B; Zero-B is a confirmation run, not a gate.
