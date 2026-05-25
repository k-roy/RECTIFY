# PERF_AUDIT.md — per-read over-computation in RECTIFY (audit + resolution)

**Created:** 2026-05-24 (Kevin's request). **Status:** living doc.
**Merged 2026-05-25:** this doc now incorporates the former
`PERF_AUDIT_FINDINGS_2026-05-24.md` (10k timing benchmarks, expensive-read outcome
analysis, priority findings, output-necessity ledger, static scan trail) — see Part II
and Part III. That separate file has been deleted; this is the single canonical record.

**Why this exists:** the 3'SS-rescue bottleneck below (run-all manifest inline-correct
stall → OOM) was the worked example that triggered this audit, but it is not the only
place RECTIFY does far more work than necessary. This doc (a) records the 3'SS case and
its **resolution**, (b) records the empirical 10k timing audit that backs it, (c) gives
the next agent a methodology + grep targets + a prioritized open-findings list, and
(d) is the chain-of-custody for the perf work. **Read before any perf work.**

RECTIFY is v0.9.x; manifest mode is experimental. Most items here are perf bugs, not
correctness bugs — but at production scale (full junction pools, high-coverage loci,
human introns >100 kb) they make pipelines hang. Fix them, but verify they don't change
output (see Verification discipline below).

> **Headline (2026-05-25): the primary 3'SS-rescue / HP-edit-distance hotspot is
> RESOLVED.** dual-site interval-tree fetch (`961c844`) + pool anchor floor with
> cross-family concordance relaxation (`ed3df74`) collapse the rescue candidate set
> from ~390/read to **0.4–2.5/read across all 5 aligners (155×–911×)**, with **66–75%
> of reads doing zero rescue work**. No real-junction loss (verified). See Part I.

---

# Part I — Case study: 3'SS rescue (the worked example, now resolved)

## The bottleneck: O(reads × full-pool × ambiguity-window) HP-edit-distance

**Where:** `rectify/core/splice/splice_aware_5prime.py` (`_rescue_3ss_truncation_body`,
`_hp_edit_distance`) fed by `rectify/core/bam/bam_processor.py:correct_read_3prime`.

**Symptom:** `rectify run-all --manifest` (DRS, 5% wt+upf1 subset) hung at
`Processing 126 regions across 16 workers / 0%` for 6.5 h then OOM-killed at 64 GB
(Sherlock job 25846844). Align + trim + Module-2H refine all succeeded; the stall is the
inline `correct` parallel region step.

**Cost decomposition (py-spy, 29,999 samples, single-thread full-pool subset = job DIAGE):**
`_hp_edit_distance` (an O(n×m) DP, splice_aware_5prime.py:693-723) = **~86% of CPU**.
- Baseline rescue `_rescue_3ss_truncation_body` = **87%** of that. Two loops scan the
  candidate set and call the DP per candidate: the sequence-rescue loop (line ~1409) and
  the Case-4 intronic-snap loop (~1877/1878).
- The sequence-rescue loop is worse than "per candidate": for EACH candidate it runs
  `for _shift in [±~15]` × `for _off in [0..10]` ≈ **~341 DP calls per candidate**.
- Terminal peel `_terminal_peel_rescue` = 12.9% (its own per-depth scan, fixed first).

**Three compounding root causes — the generalizable lessons:**
1. **Per-read scan of a global collection.** `candidate_junctions` (the junction pool)
   can hold ~17k entries in run-all mode; the body looped over all of them per read.
2. **Vestigial oversized bound.** `_POOL_SEARCH_RADIUS` was **10000** bp, with a comment
   claiming it must cover `junction_proximity_bp` "default 5000" — but the real default
   is **10**. The 10k radius was sized for a parameter that no longer exists, so it pulled
   every intron in a ±10 kb window per read.
3. **Bound keyed on the wrong axis.** The pool index is sorted by `intron_start`; a radius
   on `intron_start` bounds the *intron length*. That happens to work for yeast (introns
   ≤~1.2 kb) but **breaks for human** (introns >100 kb): the donor sits 100 kb from the
   acceptor near the read, so it falls outside any sane radius. The correct axis is
   **splice-site proximity to the read's 5' boundary**, leaving the *other* site free.

## Fixes (status as of 2026-05-25)

- DONE (commit `8e8dc8c`): peel candidate-narrowing + depth cap. **Eliminated the OOM**
  (full-pool RSS 64 GB → ~14 GB, verified by live /proc probe).
- DONE (commit `bd20f9e`): body-level `_nearby_junctions` narrowing (union of the loops'
  own cheap gates, computed once) + `_POOL_SEARCH_RADIUS` 10k→2k. Correctness-clean.
  **Did NOT fix wall-time on its own**: DIAGF (full-pool, fixed code) still stalled at
  0/126 regions after 9.5 min — each surviving nearby candidate still triggered the
  ~341-call ambiguity block, and the full pool left many nearby candidates.
- DONE (commit `961c844`): **dual-site index + either-site fetch.** Each junction is
  indexed by BOTH splice sites (per-chrom `IntervalTree`); the fetch surfaces any junction
  with *either* site within the peeled-back genomic span + buffer, NOT a radius on one
  site. Intron-length-independent (human-ready) AND tighter. Localized to
  `_build_pool_chrom_index` + the `correct_read_3prime` fetch.
- DONE (commit `ed3df74`): **pool-build anchor floor + cross-family concordance
  relaxation** — cuts the candidate *count* at its source (complement to the fetch
  narrowing). A novel observed junction enters the shared rescue/refinement pool only if
  (a) some read crosses it with a clean ≥10 bp exon anchor on both flanks (after adjacent
  D/N-op collapse), OR (b) it is reported by ≥2 *independent algorithm families*
  (`ALIGNER_FAMILY`: minimap2/gapmm2 = one family since gapmm2 wraps minimap2, vs uLTRA /
  BBMap / deSALT). Drops single-family tiny-anchor splits (gapmm2 `4M4250N223M`-style)
  without losing real short-exon-1 junctions.
- DEFERRED to human DRS (handoff `dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`
  §7): a graded **periodicity/complexity** dimension for the consensus soft gate
  (`_add_chimera_flag`). Helpers built (`_seq_periodicity`,
  `_junction_worst_flank_periodicity`) but **not wired** — a no-op on yeast (≥0.9
  periodicity fired on 0/210 gapmm2-only and 0/139 real junctions; the gene-dense genome
  has no long low-complexity anchors). Must prove its worth on human before wiring.
  **No splice-motif gating** anywhere — pruning is anchor-quality / family-support only,
  to preserve unbiased non-canonical junction discovery (memory:
  `feedback_no_motifs_unbiased_discovery`).

## Verification of the resolution (committed code)

- **Candidate collapse, 10k yeast 5-aligner data:** rescue candidates/read dropped from
  ~390 (un-narrowed whole-chromosome pool) to **mean 0.4–2.5 across all aligners** —
  gapmm2 0.43 (911×), minimap2 0.84 (464×), uLTRA 0.82 (475×), deSALT 2.52 (155×),
  mapPacBio 2.38 (164×). **66–75% of reads get 0 candidates → instant rescue skip.**
- **No real-junction loss (pool-composition equivalence):** pool 648 (no floor) → 245
  (floor only) → 287 (floor+relaxation). The 361 net-dropped junctions are **all
  single-family** (zero ≥2-family losses); the 42 relaxation-recovered are **all
  ≥2-family**. The floor removes only junctions lacking independent corroboration.
- **800-read Sherlock end-to-end** correct verify (earlier): all dropped rescues spurious.
- **1315 fast tests green** (`pytest -m "not slow"`) at HEAD.
- STILL OPTIONAL: the named `dev/verify_dualsite_correct.sbatch` full-pool completion
  (<32 GB) run and the exact end-to-end "793 s → X" wall-time. The candidate collapse
  near-eliminates the candidate-driven HP-ED work (≈86% of correct CPU), so a large
  end-to-end speedup follows, but the precise full-pool wall/RSS was not measured.

**Note:** the advisor's standing guidance is *don't* touch the shift×offset ambiguity
windows themselves (intentional splice-site-slide search). Reduce the candidate *count*
that enters them — which is exactly what the dual-site fetch + anchor floor do.

---

# Part II — Empirical audit (2026-05-24 timing + 2026-05-25 follow-up)

Static scan of the run-all pipeline, aligner wrappers, correct command,
corrected-consensus merge/BAM writers, split/chunk helpers, and manifest analysis,
combined with bounded 10k-read timing runs on yeast and human DRS subsets. Treat the
findings below as prioritized output-necessity and over-computation items; un-resolved
ones still need focused tests, scale reruns, and output-equivalence checks before a fix
is declared done.

## 10k timing benchmarks

Jobs ran the evening of 2026-05-24 PDT with random seed `20260525`. A follow-up
junction-aligner pass for uLTRA and deSALT ran 2026-05-25 against the same trimmed 10k
read sets.

| Organism | Cluster | Random 10k input | Run details |
| --- | --- | --- | --- |
| Yeast | H2 | `mex67aa_rep1_chunk_000_of_016.fastq.gz` → `mex67aa_rep1_chunk000.random10k.seed20260525.fastq.gz` | Base-aligner stage job `13475812`; instrumented base correction job `13475879`; uLTRA/deSALT follow-up job `13477937`; 4 threads for alignment, 1 thread for instrumented `correct`. |
| Human | Sherlock | SG-NEx A549 DRS `a549_chr5.primary.bam` → `a549_chr5.primary.random10k.seed20260525.bam` | Base-aligner stage job `25988503` bounded/cancelled after 31:59 (gapmm2 zero-byte PAF); instrumented base correction job `25989452`; uLTRA/deSALT follow-up job `26037612`. |

Run directories:

- Yeast: `/u/project/guillom/shared/processed/rectify_perf_timing_10k_20260525/yeast_mex67aa_rep1_chunk000/runs/13475812_20260524_183850`
- Human: `/scratch/users/kevinroy/rectify_perf_timing_10k_20260525/human_sgnex_a549_chr5/runs/25988503_20260524_183435`
- Yeast junction-aligner follow-up: `.../runs/junction_aligners_13477937_20260525_081115`
- Human junction-aligner follow-up: `.../runs/junction_aligners_26037612_20260525_081135`

Copied timing artifacts live under `/private/tmp/rectify_perf_results/{yeast,human,yeast_junction,human_junction}/`
(`timing_summary.tsv`, `instrumented_correct_summary.tsv`, `instrumented_module_timing.tsv`,
`instrumented_top_reads.tsv`, and the `junction_*` equivalents).

### Stage timing

The instrumented correction runs are intentionally single-core so they show
per-read/subsection skew without hiding it behind worker parallelism.

| Organism | Stage / subsection | Wall / observed time | Notes |
| --- | --- | ---: | --- |
| Yeast | Poly(A) trimming | 9.9s user + 1.2s sys | 248 MB max RSS. |
| Yeast | Full `run-all` stage | ~17m58s, rc=1 | 4.28 GB max RSS; failure was a wrapper/output-contract issue after substantial work. |
| Yeast | Alignment total | 614.2s | minimap2 + mapPacBio + gapmm2. |
| Yeast | minimap2 alignment | 4.2s | Fast baseline. |
| Yeast | mapPacBio alignment | 533.2s | Dominant yeast alignment cost. |
| Yeast | gapmm2 alignment | 26.7s | Left a `.paf` intermediate; not needed after BAM conversion. |
| Yeast | uLTRA alignment | 66s step / 56.3s aligner | Follow-up, rc=0; max RSS 1.29 GB. |
| Yeast | deSALT alignment | 48s step / 38.4s aligner | Follow-up, rc=0; max RSS 2.26 GB. |
| Yeast | Raw consensus selection | 47.2s | Optional for normal correct-first output. |
| Yeast | Per-aligner correction inside `run-all` | 67.7s mm2; 77.7s mapPacBio; 279.8s gapmm2 | 4-thread stage timings. |
| Yeast | Fallback consensus-BAM correction | 32.0s | Unnecessary; triggered by the manifest output-contract bug. |
| Yeast | Instrumented `correct` | 134.1s mm2; 146.3s mapPacBio; 793.3s gapmm2 | Single-core per-read timing. |
| Yeast | Instrumented junction-aligner `correct` | 298.3s uLTRA; 630.2s deSALT | Single-core; both dominated by Module 2F. |
| Human | minimap2 alignment | 13.3s | Completed. |
| Human | mapPacBio alignment | 315.7s | Slow relative to minimap2. |
| Human | gapmm2 alignment | >31m59s, cancelled | Zero-byte PAF; `run-all` lacks a `--max-intron` passthrough. |
| Human | uLTRA alignment | 57s step / 50.8s aligner | Follow-up, rc=0; max RSS 2.93 GB. |
| Human | deSALT alignment | 18s step / 15.3s aligner | Follow-up, rc=0; max RSS 4.27 GB. |
| Human | Instrumented `correct` | 36.8s mm2; 19.2s mapPacBio | gapmm2 unprofiled (no BAM emitted). |
| Human | Instrumented junction-aligner `correct` | 86.6s uLTRA; 83.2s deSALT | Single-core; completed cleanly. |

### Read-time concentration

A small subset of reads accounts for most correction time, especially in yeast and
especially in Module 2F 3'SS rescue / HP-edit-distance. Percent columns are fractions of
instrumented per-read time.

| Organism | Aligner | Reads | Read s | p50 | p95 | p99 | Max | Top 1 | Top 10 | Top 100 | HP-ED | 2F |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Yeast | minimap2 | 8,932 | 99.2s | 1.23ms | 9.06ms | 0.107s | 9.23s | 9.3% | 49.4% | 74.1% | 65.1s | 88.5s |
| Yeast | mapPacBio | 9,775 | 122.1s | 1.14ms | 8.22ms | 0.108s | 11.23s | 9.2% | 51.4% | 77.4% | 85.5s | 110.3s |
| Yeast | gapmm2 | 7,925 | 771.7s | 1.59ms | 0.359s | 2.334s | 28.07s | 3.6% | 18.0% | 55.7% | 698.9s | 760.1s |
| Yeast | uLTRA | 8,932 | 281.8s | 1.44ms | 37.76ms | 0.405s | 22.13s | 7.9% | 38.5% | 79.0% | 105.4s | 269.9s |
| Yeast | deSALT | 9,641 | 610.5s | 1.60ms | 53.43ms | 1.334s | 20.73s | 3.4% | 16.7% | 78.7% | 177.1s | 599.4s |
| Human | minimap2 | 9,985 | 23.5s | 0.98ms | 4.15ms | 7.90ms | 1.85s | 7.9% | 38.1% | 42.6% | 8.85s | 10.3s |
| Human | mapPacBio | 5,156 | 9.1s | 1.08ms | 5.27ms | 10.97ms | 0.077s | 0.8% | 3.0% | 14.6% | 0s | 0s |
| Human | uLTRA | 9,998 | 71.9s | 1.39ms | 8.06ms | 19.03ms | 13.95s | 19.4% | 56.8% | 68.6% | 46.8s | 49.0s |
| Human | deSALT | 9,991 | 69.9s | 1.28ms | 7.14ms | 29.34ms | 3.56s | 5.1% | 36.0% | 70.4% | 48.2s | 52.5s |

- Yeast minimap2/mapPacBio correction is highly outlier-driven (top 10 ≈ half, top 100 ≈
  three quarters of read-time). Yeast gapmm2 has a broader tail (55.7% in top 100, 98.4%
  in top 1000) and is much slower overall. uLTRA/deSALT follow the same 2F-dominated shape.
- Human is far lighter; mapPacBio did **zero** 2F/HP-ED (rescue never engaged), and gapmm2
  is absent (alignment stalled). So the human profile is *blind* to the worst-case 2F cost.

### Expensive-read outcome analysis

Full per-read `read_times.tsv` + `corrected_reads.tsv` were copied under
`/private/tmp/rectify_perf_full/{yeast_base,yeast_junction,human_base,human_junction}/`.
Slowest-1%-per-aligner reads classified by outcome (No-op/dead-end; Direct 3' shift;
5' only; Other flagged, no 3' shift):

| Organism | Aligner | Top 1% reads | Top 1% time | No-op in top 1% | Direct 3' | 5' only | Other, no 3' | HP time | Median cand |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Yeast | minimap2 | 90/8932 | 72.5s (73.1%) | 47 / 43.8s (60.5%) | 5 | 27 | 11 | 52.4s | 389 |
| Yeast | mapPacBio | 98/9775 | 94.3s (77.2%) | 48 / 47.1s (50.0%) | 7 | 37 | 6 | 72.1s | 389 |
| Yeast | gapmm2 | 80/7925 | 386.8s (50.1%) | 45 / 173.3s (44.8%) | 19 | 7 | 9 | 338.6s | 388 |
| Yeast | uLTRA | 90/8932 | 218.7s (77.6%) | 17 / 26.8s (12.2%) | 25 | 39 | 9 | 61.3s | 398 |
| Yeast | deSALT | 97/9641 | 476.4s (78.0%) | 10 / 34.4s (7.2%) | 44 | 36 | 7 | 75.7s | 396 |
| Human | minimap2 | 100/9985 | 10.0s (42.6%) | 49 / 0.6s (5.9%) | 38 | 6 | 7 | 8.8s | 3 |
| Human | mapPacBio | 52/5156 | 0.9s (9.4%) | 23 / 0.4s (49.2%) | 14 | 0 | 15 | 0.0s | 0 |
| Human | uLTRA | 100/9998 | 49.4s (68.6%) | 32 / 0.9s (1.8%) | 31 | 9 | 28 | 46.8s | 7 |
| Human | deSALT | 100/9991 | 49.2s (70.4%) | 4 / 0.3s (0.7%) | 37 | 40 | 19 | 47.3s | 5 |

No-op burden across the full profiled read set:

| Organism | Aligner | All no-op reads | All no-op read-time | Top-1% no-op share of all read-time |
| --- | --- | ---: | ---: | ---: |
| Yeast | minimap2 | 7072/8932 | 63.8s (64.3%) | 44.2% |
| Yeast | mapPacBio | 7669/9775 | 67.3s (55.2%) | 38.6% |
| Yeast | gapmm2 | 5001/7925 | 406.7s (52.7%) | 22.5% |
| Yeast | uLTRA | 6584/8932 | 61.4s (21.8%) | 9.5% |
| Yeast | deSALT | 6618/9641 | 89.6s (14.7%) | 5.6% |
| Human | minimap2 | 4253/9985 | 6.5s (27.5%) | 2.5% |
| Human | mapPacBio | 2082/5156 | 3.7s (41.0%) | 4.6% |
| Human | uLTRA | 4060/9998 | 10.7s (14.9%) | 1.3% |
| Human | deSALT | 4091/9991 | 8.8s (12.6%) | 0.5% |

Slowest no-op examples (read_id, no-op seconds, 3'SS candidates, HP calls):

| Organism | Aligner | Slowest no-op read | s | cand | HP |
| --- | --- | --- | ---: | ---: | ---: |
| Yeast | minimap2 | `24dd3f3c-2ab0-46ce-bc9f-060a2ff1740b` | 9.23 | 389 | 492 |
| Yeast | mapPacBio | `66f641d8-fc1d-404a-ab5a-b8abe1e491c9` | 11.23 | 390 | 615 |
| Yeast | gapmm2 | `0811b519-671f-4041-ba2c-fab96b9e8355` | 12.85 | 395 | 865 |
| Yeast | uLTRA | `66f641d8-fc1d-404a-ab5a-b8abe1e491c9` | 9.05 | 389 | 492 |
| Yeast | deSALT | `132b708e-e02b-4b35-9a1b-464fbf92e1d8` | 11.66 | 392 | 625 |
| Human | minimap2 | `55d9c950-ea99-4f87-ae6c-c01acff8abfb` | 0.07 | 4 | 0 |
| Human | mapPacBio | `67ca3eb2-02ad-4da2-89fb-455ded639fd9` | 0.08 | 0 | 0 |
| Human | uLTRA | `174b622c-e217-4ba8-8a63-ab7eb88f616e` | 0.09 | 14 | 0 |
| Human | deSALT | `39428054-d4f7-4680-b3d6-831bf5db0fac` | 0.11 | 4 | 0 |

**Outcome interpretation (and how the resolution addressed it):**
- The slow tail is **not uniformly waste** — yeast/human uLTRA/deSALT slow reads usually
  emit a real 5' rescue or 3' shift, so a blanket skip of the top-cost pattern would
  discard useful output. But yeast base-aligner profiles do carry clear dead-end compute
  (top-1% no-ops alone = 22–60% of read-time).
- The audit's own caution: **candidate count is not a safe early-exit predicate** — no-op
  and useful reads both commonly had ~388–398 candidates, so a hard HP-call/candidate
  budget risks dropping real corrections.
- **The resolution sidesteps that risk.** Rather than a heuristic budget, the dual-site
  fetch narrows candidates by *geometric* splice-site proximity and the anchor floor
  prunes the *pool* — both content-true, not a guess. The resulting 66–75% zero-candidate
  reads are reads with genuinely no junction near their 5' boundary (nothing to rescue
  against), not heuristically-skipped reads. That is the safe form of the "flag dead-ends
  early" idea the audit was circling.

### Instrumented subsection timing

`_hp_edit_distance` is nested inside `2F_3ss_rescue` (not an independent stage total).

| Organism | Aligner | Dominant instrumented subsections |
| --- | --- | --- |
| Yeast | minimap2 | `2F_3ss_rescue` 88.5s; `_hp_edit_distance` 65.1s; `2H_refine_bam_junctions` 6.4s; `2H_score_junction` 3.1s; `2E_polya_walkback_drs` 2.8s; gene-tree build 2.7s. |
| Yeast | mapPacBio | `2F_3ss_rescue` 110.3s; `_hp_edit_distance` 85.5s; `2H_refine_bam_junctions` 10.0s; `2H_score_junction` 6.8s; `2E_polya_walkback_drs` 3.2s. |
| Yeast | gapmm2 | `2F_3ss_rescue` 760.1s; `_hp_edit_distance` 698.9s; `2H_refine_bam_junctions` 8.5s; `2C_indel_correction` 3.4s; `2E_polya_walkback_drs` 2.9s. |
| Yeast | uLTRA | `2F_3ss_rescue` 269.9s; `_hp_edit_distance` 105.4s; `2H_refine_bam_junctions` 5.3s; `2C_indel_correction` 3.3s. |
| Yeast | deSALT | `2F_3ss_rescue` 599.4s; `_hp_edit_distance` 177.1s; `2H_refine_bam_junctions` 6.6s; `2H_score_junction` 5.2s. |
| Human | minimap2 | `2F_3ss_rescue` 10.3s; `_hp_edit_distance` 8.85s; `2C_indel_correction` 6.2s; gene-tree build 3.9s; `2H_refine_bam_junctions` 3.4s. |
| Human | mapPacBio | `2C_indel_correction` 5.2s; gene-tree build 4.0s; `2H_refine_bam_junctions` 1.9s; `2E_polya_walkback_drs` 1.5s. |
| Human | uLTRA | `2F_3ss_rescue` 49.0s; `_hp_edit_distance` 46.8s; `2C_indel_correction` 13.6s; `2H_refine_bam_junctions` 5.7s. |
| Human | deSALT | `2F_3ss_rescue` 52.5s; `_hp_edit_distance` 48.2s; `2C_indel_correction` 9.1s; `2H_refine_bam_junctions` 4.6s. |

Conclusions: Module **2F** (not 2H) is the correction hotspot; it is strongly tied to HP
edit-distance (yeast gapmm2 = 227,600 HP calls, 698.9s; deSALT = 1.47M HP calls on yeast).
**This is the hotspot resolved in Part I.** Secondary, still-open observations: constant
per-aligner setup (annotation/gene-tree, pool/refinement) is rebuilt per aligner (Finding
#3); `2C_indel_correction` is a visible per-read cost on human; the `run-all`
manifest output-contract bug forced a fallback consensus correction; human gapmm2 stalls
without a `--max-intron` passthrough.

Timing limitations: human used a chr5 subset (not full-genome); human gapmm2 correction
was never profiled (no BAM); Sherlock lacked `pyarrow` (DRS poly(A) metadata fell back to
TSV); 10k samples locate bottleneck *shape*, not production-scale absolute wall-time.

---

# Part III — Priority findings (open unless marked resolved) + output-necessity ledger

1. **P0: Fresh correct-first `run-all` pays for a raw alignment consensus it normally ignores.**
   `_run_alignment` calls `run_align(..., no_consensus=False)` → pre-correction consensus
   selection, `<sample>.rectified.bam`, raw stats/report, `samtools calmd`. Correct-first
   then uses per-aligner BAMs anyway. **Necessity:** optional/debug/fallback.
   **Fix:** align-only multi-aligner output, or gate raw consensus behind
   `--emit-raw-consensus` / diagnostics, preserving fallback when per-aligner correction
   is absent/fails.

2. **P0: Manifest analysis defaults emit convenience outputs that violate the documented memory contract.**
   `docs/architecture/multi_sample_pipeline.md` says manifest mode skips bedgraph +
   genomic-distribution, but `_run_analysis_manifest` passes `no_bedgraph=False` and
   `no_genomic_distribution=False`. **Necessity:** QC/convenience only.
   **Fix:** make opt-in for production manifest runs, or stream-aggregate genomic
   distribution before enabling by default.

3. **P1: Per-aligner correction rebuilds constant-per-sample junction context.**
   `_run_correction_per_aligner` runs `correct` once per aligner; each process reloads
   annotation, builds gene interval trees, loads annotated junctions, and builds/loads the
   same cross-aligner pool. `prescan` materializes `rescue_scan.pkl` / `junction_pool.pkl`,
   but fresh run-all paths don't consistently reuse them. **Fix:** auto-prescan or share a
   per-sample cache before the per-aligner loop. (Corroborated by the per-aligner setup
   cost in Part II's instrumented timing.)

4. **P1: Position-index writers are duplicated and schema-drifted.**
   Canonical `write_position_index` emits `chrom, corrected_3prime, strand, count,
   count_ag_rich`; streaming paths hand-write `*_index.bed.gz` with only the first four,
   so manifest loaders silently default `count_ag_rich=0` (lost AG-rich annotation).
   **Fix:** route streaming writers through the canonical writer, or track `(total,
   ag_rich)` in streaming accumulators — one schema.

5. **P1 — ✅ RESOLVED (`961c844` + `ed3df74`): the 3'SS pool lookup needed the dual-site index.**
   The single-coordinate index was yeast-safe only by accident (human introns put the
   opposite splice site 100s of kb away). **Resolved** by the dual-site interval-tree
   either-site fetch (`961c844`) + the pool anchor floor / cross-family relaxation
   (`ed3df74`). See Part I (candidate collapse 155×–911×, no real-junction loss).

6. **P1: Raw-consensus 5' softclip rescue scans chromosome-wide annotated junctions.**
   `select_best_alignment` collects same-chromosome annotated junctions, then
   `_rescue_5prime_softclip` does per-read distance + edit-distance work — same structure
   as the (now-resolved) 3'SS bug, for human-scale annotations. **Fix:** pass a
   per-chrom/per-strand splice-site index instead of every same-chromosome junction.
   (Still open; same dual-site treatment applies.)

7. **P2: Corrected-consensus diagnostics are emitted as if core.** `merge_corrected_tsvs`
   writes `comparison_summary.tsv` + the aligner UpSet plot; run-all writes
   `cat5_candidates.tsv`. No core step consumes them. **Fix:** gate behind a diagnostics
   flag / profile level.

8. **P2: `align` generates a dead duplicate junction BED.** `align_command.py` writes
   `<prefix>_junctions.bed` but minimap2 uses `annotation.junc.bed` via
   `get_minimap2_junc_args`; `junc_bed_path` is unconsumed. **Fix:** remove or gate.

9. **P2: gapmm2 leaves a PAF conversion sidecar.** `run_gapmm2` writes `<output>.paf` with
   no downstream consumer. **Fix:** delete by default; keep under `--keep-intermediates`.

10. **P2: `correct --emit-manifest` routes through a monolithic TSV.** Writes
    `corrected_reads.tsv`, hashes/counts the whole file, then renames to
    `..region_000.tsv` unless `--emit-merged-tsv`. **Fix:** promote true per-region
    manifest emission from streaming/checkpoint workers; reserve the monolithic TSV for
    small/legacy runs.

11. **P3: Non-streaming correction buffers all per-read results in memory.**
    `process_bam_file_parallel` accumulates all region results before writing.
    **Fix:** prefer streaming correction by default, or auto-switch above a size threshold.

## Output necessity ledger

| Artifact / output | Needed for core pipeline? | Notes / action |
| --- | --- | --- |
| Per-aligner raw BAMs (`<sample>.<aligner>.bam` + `.bai`) | Yes, until finalization | Needed for Module 2H pool/refinement, lazy HP scoring, `corrected_consensus.bam` surgery. Archive/remove only after final TSV/BAM built. |
| Raw alignment consensus (`<sample>.rectified.bam`, raw stats/report, raw `calmd`) | No for normal correct-first | Keep as fallback/debug; skip by default when per-aligner correction + corrected consensus will run. |
| Per-aligner corrected TSVs / manifests | Yes, until merge/final BAM | Inputs to merge + BAM surgery. Per-aligner corrected BAMs not required for the lazy final-BAM path. |
| Per-aligner corrected BAMs | Optional | IGV/debug; corrected-consensus writer can use raw BAMs + correction TSVs. Gate for production. |
| Final `corrected_reads.tsv` / `.manifest.tsv` + shards | Yes | Primary corrected-read table. Prefer manifest/shards at scale. |
| `corrected_reads_index.bed.gz` | Yes for manifest mode | Cluster/count loading without per-read scans. Must use the AG-rich-aware schema everywhere. |
| `corrected_consensus.bam` + `.bai` | Yes for correct-first run-all | Final inspectable consensus BAM. |
| `junction_pool.pkl` / `rescue_scan.pkl` | Internal yes | Required for chunked/prescan; should be the reusable per-sample cache in run-all (Finding #3). |
| `<prefix>_junctions.bed` from `align_command.py` | No | Redundant next to minimap2's `annotation.junc.bed`; remove or gate (Finding #8). |
| gapmm2 `.paf` intermediate | No after BAM conversion | Delete by default (Finding #9). |
| `comparison_summary.tsv`, aligner UpSet PNG, `cat5_candidates.tsv` | Optional | Diagnostics; gate (Finding #7). |
| Manifest bedgraphs + genomic-distribution plots | Optional | QC/convenience; opt-in or stream-aggregate (Finding #2). |
| Poly(A) trim metadata parquet | Yes for DRS Step 4 | Needed to restore trimmed sequence/quality context. |
| Stage provenance/checkpoint sidecars | Yes for resume/audit | Keep; avoid sidecars that duplicate recorded info or force full-file reads. |

---

# Anti-patterns to hunt for elsewhere

1. **Per-item scan of a global collection** — `for x in <big_set>:` inside a per-read
   function where a cheap gate rejects most iterations. Fix: pre-narrow ONCE via a
   sorted/bisect index, interval tree, or single-pass filter, then loop the small subset.
2. **Nested windows multiplying an expensive op.** An O(n×m) DP / edit distance inside ≥2
   nested loops (shift × offset × candidate). Cut the outer count, not the inner search.
3. **Vestigial / oversized constants.** A magic number whose comment cites a default that
   no longer exists (the 10k-for-5000 case). Re-derive against the *current* defaults.
4. **Bound keyed on the wrong axis.** A radius/window scaling with the wrong quantity
   (intron length vs splice-site proximity). Ask: does it break for human / long features?
5. **Recompute-per-item what's constant-per-read/region.** Hoist read-level work out of
   per-candidate loops (and per-aligner setup out of the per-aligner loop — Finding #3).

---

# How to hunt (methodology — profile, don't guess)

We theorized the cause twice and were wrong both times; only the sampling profile settled
it. Always measure.

- **py-spy (Sherlock rectify env).** `gdb py-bt` is UNUSABLE here (DWARF v5 vs gdb v2-4).
  Launch single-threaded (`--legacy-single-threaded -j 1`), `sleep 25` past setup, then
  `py-spy record -f raw -o stacks.folded --pid $PID --duration 300 --rate 100 --nonblocking`.
  Aggregate by leaf frame:
  `awk '{c=$NF; s=$0; sub(/ [0-9]+$/,"",s); n=split(s,a,";"); sum[a[n]]+=c} END{for(k in sum) print sum[k]"\t"k}' stacks.folded | sort -rn | head`.
- **faulthandler for hangs:** `timeout -s ABRT <sec> python -X faulthandler -m rectify.cli ...`.
- **Reproduce at scale on a dense input.** Synthetic unit tests miss this (tiny pools).
  Diag artifacts: `$SCRATCH/rectify_runall_diag_20260524/` + trimmed/aligned BAMs in
  `$SCRATCH/rectify_wt_by4742_rep1_25846844_0/`. 10k timing run dirs: see Part II.
- **For pool/candidate-count questions, you can measure on M1** with committed code +
  the 5-aligner BAMs (build the pool both ways, replicate the `bam_processor` fetch, count
  candidates/read) — no cluster needed. This is how the 155×–911× collapse was confirmed.
- **grep targets:** `for .* in candidate_` / `.*junctions` / `.*pool`; nested `for _shift`
  / `for _off`; `_hp_edit_distance` / `*edit_distance*` / `*align*` call sites;
  `_SEARCH_RADIUS` + large integer window/radius literals; `.get_aligned_pairs()` /
  genome slicing inside per-candidate loops; intervaltree `.overlap` per read.

---

# Suspect hotspots to audit next (still open)

- **`walkback_3prime_guarded` / `walkback_drs_full`** (`core/correct/walkback.py`) —
  per-read poly-A walkback; showed in the DIAGE + instrumented profiles (~2.8–3.2s yeast).
- **gene-attribution interval-tree** — `search_point` per read; per-aligner rebuild is
  Finding #3.
- **`indel_corrector`** (`2C_indel_correction`) — a notable per-read cost on human
  (5–14s in the instrumented timing); audit at human scale.
- **`junction_refiner`** (Module 2H) — shares the pool (now floored/relaxed); re-audit the
  per-N-op candidate search for the same radius-vs-axis issue at human scale.
- **Raw-consensus 5' softclip rescue** (Finding #6) — same dual-site treatment applies.
- **The 28 s non-HP-ED yeast outlier** (gapmm2 rank-1 read) — its cost is NOT HP-ED; likely
  per-candidate softclip realignment over a long read. Candidate narrowing helps but won't
  fully address it.
- **`_real_junctions` pool-floor bypass** — a read's own N-ops enter rescue candidates
  without the anchor floor; consensus catches them softly. Documented open item.

---

# Verification discipline (REQUIRE before declaring any perf fix "done")

A passing unit test is NOT evidence a perf fix works — tests use tiny inputs. Require:
1. Focused blast-radius tests green, then `pytest -m "not slow"`.
2. **Scale re-run COMPLETES** on a dense representative input AND peak memory bounded
   (target <32 GB for the 5% subset).
3. **Measured** wall-time / call-count (or candidate-count) drop — not "it feels faster."
4. **Output unchanged** — rescue counts / corrected-TSV within tolerance of the pre-fix
   run (no over-narrowing that silently drops real rescues). For the pool floor this was
   the pool-composition equivalence + 800-read end-to-end check in Part I.

---

# Static scan trail (2026-05-24 audit provenance)

Files/functions reviewed in the static pass that produced the Part III findings:
`core/commands/run/stages.py`, `run/single_sample.py`, `run/multi_sample.py`,
`core/commands/align_command.py`, `core/align/multi_aligner.py`,
`core/consensus/select.py`, `consensus/scoring.py`, `consensus/corrected_consensus.py`,
`core/commands/correct_command.py`, `core/bam/parallel.py`, `core/position_index.py`,
`core/analysis/manifest.py`, `analysis/loaders.py`, `core/commands/prescan_command.py`,
`core/commands/split_command.py`, and the chunked/split generated script paths.

Greps used: `no_consensus`, `junc_bed_path`, `comparison_summary`, `cat5_candidates`,
`write_position_index`, `count_ag_rich`, `no_bedgraph`, `no_genomic_distribution`,
`.paf`, `candidate_junctions`, `annotated_junctions`.

---

# Pointers

- `AGENT_FIXES.md` → `[2026-05-24] ... 3'SS rescue` entry (full diagnosis + commands).
- Resolution commits: `8e8dc8c` (peel/OOM), `bd20f9e` (narrowing), `961c844` (dual-site),
  `ed3df74` (anchor floor + cross-family relaxation), `93d99b8` (this doc + human handoff).
- Human-readiness follow-ups: `dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`
  (periodicity validation, `--annotation` tractability, A549 equivalence).
- Diagnostic jobs A–F + py-spy folded stacks: `$SCRATCH/rectify_runall_diag_20260524/`.
- 10k timing run dirs + copied artifacts: see Part II.
