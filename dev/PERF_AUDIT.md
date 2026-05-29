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

---

# Part IV — Per-read profile POST all 2026-05-26 fixes (task #12, 2026-05-26)

**Baseline: HEAD `c1d0e83`** (NOT `cf5ebb9` as the task brief stated). Between
`cf5ebb9` and `c1d0e83` the loop-1 shift×offset cost the brief said was "deliberately
left untouched" was in fact ALREADY narrowed by four commits: `da18965`
(`align_clip_to_exon` 500 bp guard), `fe1bb69` (terminal-peel K=20 cap), `25d7a30`
(baseline `_rescue_3ss_truncation_body` K=25 candidate cap), `5ceb243` (offset-loop ED=0
early-exit), plus `cf5ebb9` (rRNA/Pol III exclusion). The "soft-hung at region 29/126 /
216 s per region" anecdote in the brief predates the K-caps — do NOT expect to reproduce
it on current code. The right question, answered here: **given all those fixes, what STILL
dominates `correct_read_3prime` and which reads are dead ends?**

## Method

Standalone profiler `dev/profile_correct_reads.py` (NOT a hot-path change) replicates the
`correct` setup exactly: load genome + GFF, build the cross-aligner junction pool from the
4 staged per-aligner BAMs, build `pool_chrom_index` via `_build_pool_chrom_index`, build
the `ExclusionRegionDetector` (rDNA/Pol III, default-on), build gene interval trees; then
time `correct_read_3prime` per read recording chrom/pos/strand/mapq, CIGAR complexity
(n_ops, n_S/I/D/N/X, ref_len), 5'-softclip len, `n_pool_nearby` (pool candidates in the
real fetch window), `n_op_intervals`, `is_excluded`, and the emitted outcome
(five_rescued / three_shifted / no-op). Sorts by elapsed, writes top-1% + all rows to TSV.

- **Input:** Sherlock staged `wt_by4742_rep1.{minimap2,mapPacBio,gapmm2}.bam` (yeast DRS),
  6000 reads each, single-thread. Job 26151151 (profile) + 26152138 (py-spy + cross-chrom).
- **py-spy** (job 26152138) on the single worst read of each aligner, 200 Hz × 25 s.
- Output TSVs copied to M1: `dev/perf_task12_results/profile_{aligner}.{top1pct,all}.tsv`.
- **Cross-chrom corroboration (job 26152138):** a stride-40 (1-of-40) pass over the WHOLE
  coord-sorted BAM (8000 sampled reads spanning all 16 chroms) was run for gapmm2 + minimap2
  to test the chrI-only limitation. It is decisive — see Headline (B). TSVs:
  `dev/perf_task12_results/profile_{gapmm2,minimap2}_xchrom.{top1pct,all}.tsv`.
- **Limitations:** (1) yeast only — human chr5 dense-splice loci (per the
  `[2026-05-26] A549 _hp_edit_distance slow` AGENT_FIXES entry) likely look like the
  cross-chrom (B) profile or worse, not the chrI (A) profile; (2) the **chrI 6k pass is a
  splice-sparse LOWER bound** (`fetch(until_eof=True)` starts at chrI) — always read it
  alongside the cross-chrom (B) numbers, which are the representative picture; (3)
  `apply_indel_correction` was left
  **OFF** so the profile isolates the 3'SS-rescue cost — production `correct` enables 2C by
  default (`not skip_indel_correction and not is_short_read`). The production 2C cost adds
  on top of the 3'SS-rescue cost shown here (per Part II, `2C_indel_correction` was 5–14 s
  on human, smaller on yeast); it is a separate, already-logged open finding and does not
  affect the 3'SS-rescue conclusions below.

## Headline

**The pre-fix per-read pool *scan* is gone, but on splice-DENSE chromosomes the
candidate-driven HP-ED cost RESURGES — and a few catastrophic reads (one at 292 SECONDS)
dominate.** The chrI-only pass looked clean (top-1% median candidates 0–4); the cross-chrom
stride pass (below) shows the chrI numbers were a splice-sparse lower bound. Two views:

**(A) chrI 6k subset (splice-sparse — LOWER BOUND):** candidate counts look collapsed
(top-1% median `n_pool_nearby` 0–4), cost is a handful of outliers at 1.8 / 6.5 / 12.3 s.

| Aligner (chrI 6k) | total | p50 | p99 | max | top-1% share | top-1% no-ops | top-1% rescued |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| minimap2 | 6.15 s | 0.58 ms | 1.56 ms | **1761 ms** | 38.5% | 52/60 (12.7% of top time) | 4 |
| mapPacBio | 21.0 s | 0.49 ms | 1.55 ms | **6476 ms** | 83.8% | 29/60 (1.6%) | 25 |
| gapmm2 | 29.85 s | 0.60 ms | 2.14 ms | **12349 ms** | 86.7% | 33/60 (20.5%) | 23 |

**(B) cross-chrom stride-40 8k subset (spans all chroms — the REAL picture):** total time is
**20–100× higher** than chrI; the worst single read is **292 s** (gapmm2 chrIV:1360492);
top-1% candidate medians rise to **6** and the worst reads carry **3–28 nearby candidates**
(some ABOVE the K=25 cap); **most expensive reads are PRODUCTIVE rescues, not dead ends**
(72–74 of top-80 five_rescued).

| Aligner (xchrom 8k) | total | p50 | p95 | p99 | max | top-1% share | top-1% no-ops | top-1% median cand |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| minimap2_xchrom | **618 s** | 0.81 ms | 84.6 ms | 331 ms | **87194 ms** | 83.9% | 4/80 (0.9%) | 6 |
| gapmm2_xchrom | **595 s** | 0.76 ms | 90.0 ms | 296 ms | **292124 ms** | 84.0% | 7/80 (1.4%) | 6 |

The cross-chrom disjoint classes confirm the candidate-driven cost is back at scale:
- **minimap2_xchrom:** `nearby_cand_productive` = **61.5%** (1533 reads), high-op 17.6%,
  large-clip 14.0%, dead-end-noops 6.5%. **235 reads exceed the K=25 cap (45 s).**
- **gapmm2_xchrom:** `large_subguard_5clip` = **62.6%** (just 7 reads — the 292/45/31 s
  monsters), `nearby_cand_productive` 22.0%, dead-end-noops 7.1%, bigN 2.9%. **167 reads
  exceed the K=25 cap = 300 s = HALF the total.**

So the K=25 cap (`25d7a30`) helps but is *insufficient* on dense loci: reads with 27+
candidates still hit the ceiling, and the worst reads combine **large 5' clip + many
candidates + high op count** into a single catastrophic DP volume (the 292 s read has
five_clip=131, n_cigar_ops=100, n_pool_nearby=27).

## py-spy: the residual cost is uniformly `_hp_edit_distance`

py-spy on the single worst read of **all three** aligners (the 1761 ms minimap2 / 6476 ms
mapPacBio / 12349 ms gapmm2 reads) showed **~100% of samples in `_hp_edit_distance`,
lines 727–734 = the inner O(n×m) DP cell loop.** NOT `align_clip_to_exon` (Gotoh), NOT
`get_aligned_pairs`/genome-slicing, NOT 2C. So even after the candidate cap, the leaf cost
is the same DP as before — just reached via a different driver. The DP is reached
`(candidates) × (shift range ~21) × (offset range ≤ junction_proximity_bp)` times per read,
each call up to the 200×200 `_HP_ED_MAX_LEN` cap. The K=25 cap bounds *candidates*; it does
**not** bound the shift×offset product per candidate, and the ED=0 early-exit (`5ceb243`)
only fires when a perfect match exists — **dead-end reads (no good match) never hit ED=0,
so they pay the full shift×offset×DP cost.**

Why these specific reads stay expensive even with `n_pool_nearby=0`:
- **Annotated junctions + the read's own `_real_junctions` bypass the K=25 pool cap**
  (`bam_processor.py` ~L429–432 add them to `_ss_junctions` unconditionally; the cap is on
  the pool-fetched subset). A read near several annotated 3'SS still feeds many candidates.
- The **mapPacBio leading-D N-op-matched** junctions are *deliberately preserved* past the
  cap (correctness — they are the evidence basis), so garbled high-op mapPacBio reads run
  the full DP on every N-op-matched junction.
- A **large 5' soft-clip** (e.g. the 456 bp minimap2 outlier) yields a long `rescue_seq`
  (capped at 200), making each DP call hit the 200×200 ceiling — maximal per-call cost.

## Disjoint cost-class attribution (chrI 6k, by precedence)

Classes are assigned by precedence (bigN > large-clip > high-op > dead-end-noop > productive
> other) so percentages are disjoint and sum to ~100%. ("other" = the long tail of normal
reads at baseline ~0.5 ms each — not a hotspot.)

| Class (precedence) | minimap2 | mapPacBio | gapmm2 |
| --- | ---: | ---: | ---: |
| `bigN_artifact_intron` (n_N>2 kb) | 0.0% (1 rd) | 0% (0) | **49.9%** (28 rd) |
| `large_subguard_5clip` (100–499 bp clip) | **32.8%** (23 rd) | 0% (0) | 15.1% (12 rd) |
| `high_cigar_op` (≥50 CIGAR ops) | 9.9% (571 rd) | **83.1%** (1500 rd) | 8.3% (1796 rd) |
| `deadend_noop_with_cands` (no-op, n_pool_nearby>0) | 7.7% (813 rd) | 1.3% (310 rd) | **16.8%** (754 rd) |
| `nearby_cand_productive` | 0.9% | 5.7% | 3.0% |
| `other` (baseline tail) | 48.6% | 9.9% | 6.9% |

Per-aligner dominant driver:
- **minimap2** — large (but sub-500-guard) 5' soft-clips: 23 reads = 33% of total. The
  `da18965` 500 bp `align_clip_to_exon` guard does NOT catch these (the worst is 456 bp).
- **mapPacBio** — garbled high-CIGAR-op alignments (≥50 ops): 1500 reads = 83%. These are
  the BBMap PacBio-error-model artifacts on ONT (see the mapPacBio `intronlen`/`maxindel`
  AGENT_FIXES entries) — many N-op-matched candidates each running the full DP.
- **gapmm2** — `bigN_artifact_intron` (a single read with an **88 kb spurious N-op**,
  n_N=88240, alone = 12.3 s ≈ 41% of total), plus dead-end no-ops (16.8%) and large clips
  (15.1%). The 88 kb N is a gapmm2 mis-relabel; the cost is still the DP volume, not the N.

**Cross-chrom (xchrom 8k) disjoint classes — the representative picture (see Headline B):**
on dense chroms the dominant class shifts to **productive candidate-driven rescue**, NOT
dead ends: minimap2_xchrom `nearby_cand_productive` 61.5% / high-op 17.6% / large-clip 14.0%
/ dead-end-noop 6.5%; gapmm2_xchrom `large_subguard_5clip` 62.6% (7 monster reads) /
`nearby_cand_productive` 22.0% / dead-end-noop 7.1%. **167 gapmm2 + 235 minimap2 reads
exceed the K=25 cap (300 s + 45 s respectively)** — i.e. half of gapmm2's total cross-chrom
cost is in reads that hit the candidate ceiling. The cost is genuine rescue work on
candidate-rich reads, so the lever is per-candidate DP-skipping (#1) and a tighter/unified
candidate cap (#3), NOT blanket read skipping.

## Dead-end classes (reads where rescue can never succeed → should early-skip)

The cleanest flaggable dead-end class: **`is_noop=1 AND (n_pool_nearby>0 OR n_op_intervals>0)`**
— reads that RUN the rescue DP (had candidates) but emit zero correction (no five_prime
rescue, no 3' shift, no category). Whole-set burden:

| Aligner | dead-end-with-work reads | their time | % of total |
| --- | ---: | ---: | ---: |
| gapmm2 (chrI) | 862 | 5.64 s | **18.9%** |
| minimap2 (chrI) | 838 | 0.49 s | 7.9% |
| mapPacBio (chrI) | 334 | 0.29 s | 1.4% |
| gapmm2_xchrom | 2274 | 42.2 s | 7.1% |
| minimap2_xchrom | 2248 | 40.0 s | 6.5% |

These are reads near a junction (or with an N-op) whose `rescue_seq` matches NO candidate
exon-1 window well enough — they exhaust the full shift×offset×DP search and return "none".
They are genuine dead ends: the DP can be skipped entirely if a *cheap* pre-gate proves no
candidate can match. NOTE: on the representative cross-chrom set the dead-end fraction is
SMALLER (6.5–7.1%) than on chrI (because dense chroms have more *productive* candidate-rich
reads inflating the denominator) — but in absolute terms it is still 40+ s/8k reads, and the
same cheap pre-gate that skips dead-ends (#1) also skips the hopeless *candidates within*
productive reads, so it attacks both the dead-end and the productive-but-bloated classes.

## Concrete speedup + dead-end-flagging recommendations (UNVERIFIED — propose, don't ship)

Per the standing advisor guidance, do NOT touch the shift×offset windows themselves; cut the
*number of DP calls* or *skip provably-hopeless reads*. Candidates, in priority order:

1. **Cheap pre-gate before the shift×offset DP (kills the dead-end class).** For each
   candidate junction, before the `for _shift … for _off … _hp_edit_distance` block, compute
   a *cheap* lower-bound on achievable ED — e.g. a k-mer (k=6–8) presence/Jaccard between
   `rescue_seq` and the exon-1 window at the unshifted position, or a single Hamming at
   shift=0. If the best cheap score across all shifts of that candidate already exceeds
   `max_edit_frac * rescue_len`, skip the full DP for that candidate. Dead-end reads (no
   candidate can match) then skip ALL DP calls. **Risk:** a too-aggressive bound drops real
   rescues — must verify rescue-count equivalence (Part I discipline). This is the
   "flag dead-ends early" idea the audit was circling, now grounded: the 862 gapmm2 dead-end
   reads (18.9%) are the target.

2. **Cap `rescue_seq` length used in the DP (helps the large-clip class).** A 456 bp clip is
   truncated to `five_clip` then the DP caps at 200 — but a real exon-1 boundary match needs
   only the ~20–40 bp adjacent to the splice site. Truncating `rescue_seq` to e.g. the 60 bp
   nearest the 5' boundary before the DP would shrink each call from 200×200 to ~60×~80.
   **Risk:** changes which window scores best on legitimately long short-exon-1 reads —
   verify on the validation set.

3. **Unify + tighten the K cap (the cross-chrom data shows K=25 is insufficient on dense
   loci).** Today the K=25 cap (`_rescue_3ss_truncation_body`) bounds only the pool-fetched
   subset; annotated + read-own N-op junctions are added unconditionally. The cross-chrom run
   is hard evidence this matters: **167 gapmm2 reads exceed 25 nearby candidates and account
   for 300 s = HALF of gapmm2's 595 s total**; 235 minimap2 reads (45 s). The worst read (292 s)
   has n_pool_nearby=27. Apply the edge-distance K cap to the *union* of annotated + pool +
   N-op-matched candidates (preserving N-op-matched ones, as the existing partition does), and
   consider lowering K (e.g. 10–15) on dense loci. **Risk:** could drop a real rescue if the
   correct junction is the 26th-closest — bound the risk by ranking annotated junctions ahead
   of novel ones and verifying rescue-count equivalence. This is independent of #1 and stacks
   with it.

4. **Leave the `align_clip_to_exon` 500 bp guard (`da18965`) alone — it does NOT apply to
   this hotspot.** The 456 bp minimap2 outlier sits just under the guard, so it looked like a
   candidate for tightening 500→150 bp — but py-spy showed its cost is in `_hp_edit_distance`,
   NOT the Gotoh `align_clip_to_exon`. Tightening the guard would not touch the large-clip
   hotspot; recommendation #2 (cap `rescue_seq` length fed to the DP) is the actual lever for
   that class. Recorded here only to prevent a future agent from chasing the guard.

5. **(mapPacBio-specific, low ROI) reconsider mapPacBio in the ONT panel.** The high-op class
   is 83% of mapPacBio's correct-stage cost and stems from BBMap's PacBio error model
   producing garbled ONT alignments (Křižanović 2018; see mapPacBio AGENT_FIXES). The cost
   is real rescue work, not waste — but it is the worst aligner-per-read profile. Per the
   `[2026-05-25]` entry, whether mapPacBio belongs in an ONT splice panel at all is the
   higher-leverage question than micro-optimizing its rescue.

**None of these are shipped.** Each needs the Part I verification discipline (focused tests
+ `pytest -m "not slow"` + scale rerun + rescue-count/output equivalence) before landing.

## Artifacts

- Profiler: `dev/profile_correct_reads.py` (M1 + Sherlock
  `…/software/rectify/dev/`). Sbatch: `dev/profile_correct_task12.sbatch`,
  `dev/profile_pyspy_task12.sbatch`.
- Top-1% + all-row TSVs (M1): `dev/perf_task12_results/profile_{minimap2,mapPacBio,gapmm2}.{top1pct,all}.tsv`
  (chrI 6k) and `profile_{gapmm2,minimap2}_xchrom.{top1pct,all}.tsv` (cross-chrom 8k).
- py-spy folded stacks (Sherlock): `…_pyspy_…/stacks_{minimap2,mapPacBio,gapmm2}.folded`
  (all three → ~100% `_hp_edit_distance` inner DP).
- Sherlock run dirs: `$SCRATCH/rectify_perf_task12_20260526_114334/` (chrI profile),
  `$SCRATCH/rectify_perf_task12_pyspy_20260526_115051/` (py-spy `.folded` + cross-chrom).
- Jobs: 26151151 (profile), 26152138 (py-spy + cross-chrom stride). Sherlock rectify synced
  to `c1d0e83` (`bam_processor.py` was the only file behind; md5-verified) before running.

---

# Part V — Genome-wide expensive-read panel + cost archetypes (task #13, 2026-05-26)

**Baseline: HEAD `c52966c`** (3 commits past Part IV's `c1d0e83`; only branding/docs commits
between them — no rescue-path code change, so the cost *shape* matches Part IV (B), the
cross-chrom picture). This part REPLACES Part IV's top-1%/stride-40 sampling with a
**genome-wide stride-5 sweep across all 5 aligners, collecting EVERY read above a fixed
500 ms cost threshold** into a panel, then categorizing the panel into disjoint cost
archetypes by precedence.

## Method

`dev/profile_correct_reads.py` gained a **panel mode** (`--cost-threshold-ms`): it streams
every read whose `correct_read_3prime` elapsed_ms ≥ threshold to `<prefix>.panel.tsv` while
accumulating whole-set aggregates (reads scanned, total time, a log-spaced time-share
histogram) in `--summary-out` — so the whole BAM is swept without buffering ~40k rows. The
`correct` setup is replicated exactly (genome + GFF, 5-BAM cross-aligner pool = 3,459
junctions / 385 annotated, `pool_chrom_index`, rDNA/Pol III `ExclusionRegionDetector` = 610
regions, gene interval trees); `apply_indel_correction` is OFF to isolate the 3'SS-rescue
cost (matches Part IV; 2C is a separate logged finding).

- **Input:** Sherlock `wt_by4742_rep1.{minimap2,mapPacBio,gapmm2,deSALT,uLTRA}.bam` (yeast DRS).
- **Sweep:** stride-5 (1-of-5 primary reads, all 16 chroms), cap 120k reads/aligner, threshold
  **500 ms**. SLURM array job **26175498** (one task per aligner, single-thread).
- **Threshold calibration:** Part IV `c1d0e83` xchrom stride-40 8k → 42 (gapmm2) / 59 (minimap2)
  reads >500 ms. Scaling to stride-5 projected ~1000–1500 pooled — confirmed in range.
- **Categorizer:** `dev/categorize_perf_panel.py` (disjoint precedence; count + panel-time share).
- **py-spy:** `dev/profile_pyspy_task13.sbatch` on the worst read/aligner (job 26181308).

## ⚠️ Caveat — PARTIAL SWEEP (chrI–chrVII)

The genome-wide array is **single-thread on a shared larsms node and crawls** (~80 min wall to
reach chrVII of XVI; the heavy per-read tail dominates). The numbers below are the **first 200
>500 ms reads flushed per aligner = 1000 pooled reads spanning chrI–chrVII**, which **includes
the dense, expensive chrIV** (where Part IV found its 292 s monster). Reads 201+ are buffered in
Python and flush only at `close()`; the **full-genome panel + `summary_*.tsv` (histogram +
panel-share-of-total denominator) are pending job completion** (background array 26175498). The
remaining chroms (VIII–XVI) add chrXII — largely rDNA-**excluded** — and the smaller chromosomes;
the disjoint archetype *ratios* below are decisive and are not expected to invert with their
addition. Treat the table as a representative chrI–VII picture, not the final full-genome panel.

## Panel sizing + archetype breakdown (chrI–VII, 1000 pooled reads, 6,688 s panel time)

Per-aligner each contributed its first 200 >500 ms reads. The worst single reads are extreme:
**gapmm2 320 s** (chrIV, n_pool=31>K25, five_clip=93, n_ops=109), **minimap2 300 s** (chrIV,
five_clip=131, n_pool=31), **uLTRA 181 s** (five_clip=**497**), mapPacBio 154 s (n_ops=129),
deSALT 114 s (n_ops=86).

| Archetype (precedence) | predicate | POOLED reads | POOLED % of panel time |
| --- | --- | ---: | ---: |
| `excluded_BUG` | `is_excluded==1` | **0** | **0.0%** ✅ (rRNA/PolIII gate clean) |
| `bigN_artifact_intron` | `n_N > 2000` | 25 | 1.1% |
| `large_5prime_clip` | `five_clip >= 100` | 106 | **35.9%** |
| `over_cap_candidates` | `n_pool_nearby > 25` | 28 | 5.4% |
| `high_cigar_op` | `n_cigar_ops >= 50` | 197 | **28.5%** |
| `dead_end_with_work` | `is_noop==1 AND (n_pool_nearby>0 OR n_op_intervals>0)` | 66 | 0.8% |
| `productive_candidate` | `(5'rescued OR 3'shifted) AND n_pool_nearby>0` | 576 | 28.4% |
| `other` | baseline tail | 2 | 0.1% |

**Per-aligner dominant class differs sharply** (one optimization won't help all equally):

| Aligner | panel time (s) | dominant class (% of its panel time) |
| --- | ---: | --- |
| minimap2 | 1706 | `large_5prime_clip` 58.8% |
| gapmm2 | 1884 | spread: `high_cigar_op` 22.4%, `productive` 22.7%, `over_cap` 17.4%, `large_clip` 35.2%* |
| mapPacBio | 1039 | `high_cigar_op` **70.6%** |
| deSALT | 884 | `productive_candidate` 49.0% |
| uLTRA | 1175 | `large_5prime_clip` 43.7% |

*gapmm2 large_clip+high_op+over_cap together ≈ 75% — the catastrophic-combo aligner.

Feature medians confirm the classes: `large_5prime_clip` median elapsed **4636 ms**, median
five_clip 244 bp; `over_cap_candidates` median n_pool_nearby **33** (binding above K=25);
`high_cigar_op` median n_ops 86; `bigN_artifact_intron` median n_N ~100 kb.

## What CHANGED vs Part IV (the panel overrules Part IV's priorities)

1. **`dead_end_with_work` collapsed from Part IV's chrI 18.9% to 0.8% here.** The genome-wide
   sample has far more *productive* candidate-rich reads inflating the denominator, and the
   post-`c1d0e83` code reduced absolute dead-end cost. **Recommendation #1 (cheap pre-gate to
   skip dead-ends) now targets a much smaller class** and is re-ranked BELOW the large-clip and
   K-cap levers.
2. **`large_5prime_clip` is now the single biggest class (35.9%)** — directly actionable by
   capping `rescue_seq` length fed to the DP (#2). The uLTRA 497 bp / minimap2 131 bp / gapmm2
   93 bp worst-read clips are the exemplars.
3. **`over_cap_candidates` = 5.4% pooled but 17.4% on gapmm2, median n_pool_nearby=33** — hard
   empirical evidence K=25 is still binding on a real subset. Recommendation #3 has support.
4. **`excluded_BUG = 0` is a positive VERIFICATION** — the rRNA/Pol III exclusion gate
   (`bam_processor.py:423-426`, `cf5ebb9`) correctly keeps excluded reads out of the >500 ms
   panel. No exclusion-gate leak.

## py-spy corroboration — CONFIRMED ~100% `_hp_edit_distance`

py-spy on the worst read of each aligner (`dev/profile_pyspy_task13.sbatch`, 200 Hz × 25 s, job
26181577; `26181308` first failed in 4 s with exit 141/SIGPIPE from a `... | head -1` under
`pipefail` — fixed by dropping `set -e`/`pipefail`). Profiled the panel's worst reads (gapmm2
320 s, minimap2 300 s, mapPacBio 154 s). **Result is decisive — the residual cost is the inner
O(n×m) DP cell loop, nothing else:**

| Aligner (worst read) | `_hp_edit_distance` share | top leaf frames |
| --- | ---: | --- |
| gapmm2 (320 s, chrIV) | **99.9%** (4992/4999) | `splice_aware_5prime.py:731,734,732` |
| minimap2 (300 s, chrIV) | **99.7%** (4984/4999) | `splice_aware_5prime.py:731,734,732` |
| mapPacBio (154 s, chrV) | **99.7%** (4982/4999) | `splice_aware_5prime.py:731,734,732` |

Lines 727–734 are the inner DP cell loop of `_hp_edit_distance`. NOT `align_clip_to_exon`
(Gotoh), NOT `get_aligned_pairs`/genome-slicing, NOT 2C — matching Part IV exactly. So even
after the candidate cap + early-exits, the catastrophic-read leaf cost is the same DP, reached
via (candidates) × (shift range ~21) × (offset range ≤ `junction_proximity_bp`) calls per read,
each up to the 200×200 `_HP_ED_MAX_LEN` cap. This is what recommendations #1 (truncate
`rescue_seq` → smaller DP) and #2/#3 (fewer DP calls) attack. *(deSALT + uLTRA were still
profiling when harvested; the three above are unanimous and representative.)* Folded stacks:
`dev/perf_panel_t13_results/stacks_{gapmm2,minimap2,mapPacBio}.folded`.

## Optimization lessons for the expensive top-1% (UNVERIFIED — propose, don't ship)

Standing advisor guidance: do NOT touch the shift×offset windows (intentional splice-slide
search); cut the *number of DP calls* or *skip provably-hopeless reads*. **Code-accurate framing
(verified against `c52966c`):** the K=25 cap in `_rescue_3ss_truncation_body`
(`splice_aware_5prime.py:1298`) IS already unified across the annotated + read-own
(`_real_junctions`) + pool candidate union for the *proximity-gated* portion — all flow through
`_ss_junctions` (`bam_processor.py:428-449`) into `_nearby_junctions` then the sort-and-slice.
Only **N-op-matched** junctions are intentionally preserved past the cap. So the lever is NOT
"unify the cap" (Part IV's framing — inaccurate) but the items below. Priority is set by the
panel time-share (large-clip > high-op ≈ productive > over-cap > dead-end):

1. **Cap `rescue_seq` length fed to the DP — TOP lever (targets `large_5prime_clip`, 35.9%).**
   A real exon-1 boundary match needs only ~20–60 bp adjacent to the splice site, but a long 5'
   clip yields a `rescue_seq` up to `_HP_ED_MAX_LEN`=200, so every DP call at
   `splice_aware_5prime.py:1496` hits the 200×200 ceiling. Truncate `rescue_seq` to e.g. the
   60 bp nearest the 5' boundary before the shift/offset loop → each call drops to ~60×~80
   (~6× cheaper). The worst panel reads (497/131/93 bp clips) are the targets. **Risk:** changes
   which window scores best on legitimately long short-exon-1 reads — verify rescue-count /
   corrected-3end equivalence on the validation set.

2. **Lower K and/or narrow the proximity window feeding the K=25 cap (targets
   `over_cap_candidates`, 5.4% pooled / 17.4% gapmm2; two worst reads at n_pool=31).** Options:
   lower K to 10–15 ranking annotated junctions ahead of novel ones; OR shrink
   `_POOL_FETCH_HALF_WINDOW` / `junction_proximity_bp`. **Risk:** drops a real rescue if the
   correct junction is the 26th-closest — bound by annotated-first ranking + rescue-count
   equivalence. Independent of #1, stacks with it.

3. **Cheap per-candidate pre-gate before the shift×offset DP (targets the now-small
   `dead_end_with_work` 0.8% AND the hopeless candidates inside `productive_candidate`).** Before
   the `for _shift … for _off … _hp_edit_distance` block (L1475-1525), compute a cheap ED lower
   bound — k-mer (k=6–8) Jaccard or single Hamming at shift=0 between `rescue_seq` and the
   unshifted exon-1 window. If the best cheap score already exceeds `max_edit_frac * rescue_len`,
   skip the full DP for that candidate. The ED=0 early-exit (`5ceb243`) only fires on a *perfect*
   match, so dead ends pay the full cost today. **Re-ranked below #1/#2** — the dead-end class
   shrank from Part IV's 18.9% to 0.8%, though the pre-gate also trims hopeless candidates inside
   productive reads. **Risk:** a too-loose bound drops real rescues — verify equivalence.

4. **`high_cigar_op` (28.5% pooled, 70.6% of mapPacBio) is mostly real candidate-driven work on
   garbled alignments, not a clean skip.** It stems from BBMap's PacBio error model on ONT
   (Křižanović 2018; mapPacBio `intronlen`/`maxindel` AGENT_FIXES entries). #1+#3 help indirectly
   (these reads have many candidates each running the DP); the higher-leverage question is
   whether mapPacBio belongs in an ONT splice panel at all (`[2026-05-25]` entry).

5. **Leave the `align_clip_to_exon` 500 bp guard (`da18965`) alone** — py-spy (Part IV) showed
   the large-clip cost is in `_hp_edit_distance`, NOT the Gotoh `align_clip_to_exon`. #1 is the
   lever; recorded to prevent a future agent chasing the guard.

**None shipped.** Each needs the Part I verification discipline (focused tests + `pytest -m "not
slow"` + scale rerun + rescue-count/output equivalence) before landing.

## Artifacts

- Profiler: `dev/profile_correct_reads.py` (panel mode, `--cost-threshold-ms` + `--summary-out`).
  Sbatch: `dev/profile_panel_task13.sbatch`, `dev/profile_pyspy_task13.sbatch`. Categorizer:
  `dev/categorize_perf_panel.py`.
- Panel TSVs (M1, chrI-VII): `dev/perf_panel_t13_results/profile_<aligner>.panel.tsv`;
  provenance `dev/perf_panel_t13_results/PROVENANCE.json`.
- Sherlock run dir: `$SCRATCH/rectify_perf_panel_t13_26175498/` (full-genome panel +
  `summary_*.tsv` land here on array completion). Jobs: 26175498 (panel array), 26181308 (py-spy).
- Sherlock rectify md5-verified at `c52966c` (splice_aware_5prime / bam_processor /
  junction_scoring all match M1) before profiling.
