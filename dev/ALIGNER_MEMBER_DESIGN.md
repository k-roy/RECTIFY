# RECTIFY native-aligner member — design doc (2026-06-18)

> Source: crafter fan-out `wfpiu0qap` (5 concept-architects → architect → adversarial redteam).
> Redteam verdict: **sound_with_fixes**. Cycle scope: validated simulation benchmark + THIS design doc only;
> member build deferred and gated on the benchmark. Every claim is a hypothesis-pending-ablation;
> fitness = the simulation truth set, NEVER the internal score.

## ⚠ REDTEAM MUST-FIX (apply before any member code — corrections to the body below)
1. **DP CONFLATION (blocker, was seeded in the brief):** `align_exon_block_global` (local_aligner.py:522)
   is called at **read_edits.py:811 ONLY** (exon-block realignment). The 5'-rescue DP is a DIFFERENT
   function, `align_clip_to_exon` (local_aligner.py:716), called at junction_refiner.py:717 +
   splice_aware_5prime.py:1930/2052. They are two distinct live affine DPs — the "ONE DP substrate" claim
   must be re-scoped per-DP. C1 extends `align_exon_block_global`; do not present them as one recurrence.
2. **Prior-art attribution (reframe novelty as integration/orthogonality, not method invention):**
   C3 = pair-HMM (Durbin) / Dindel (Albers 2011) / GATK PairHMM; C5 = FracMinHash (sourmash, Irber 2022) /
   Mash containment (Ondov 2016, Koslicki-Brown); C4 = RATTLE / isONcorrect / FLAIR pool→POA→project recipe.
   The novelty is applying these as an ORTHOGONAL member to a seed-and-chain panel, not the methods.
3. **"integer" → "flat float":** the current selection path is float-valued (score_segment uses 0.5/3.0
   multipliers; junction_score is float, extract.py:77). Replace every "integer-argmax/integer path" with
   "flat float / flat-affine discrete-tie-prone". The tie-broken-by-LLR argument survives; the word does not.
4. **C4 ablation must be position-exact:** net change in EXACT-INDEL-POSITION concordance (ambiguity-aware)
   strictly positive AND zero minority-isoform reassignments — not a raw "errors introduced < removed" count.
5. **Discovery-FDR ↔ Q-provenance gate:** a stale/chemistry-mismatched Q-recalibration stamp must force the
   non-canonical discovery LLR to ABSTAIN (the de-novo non-canonical LLR depends on recalibrated ONT Q).

Other minor corrections in the redteam appendix at the end of this file.

---

# RECTIFY Native Aligner Member — Design Doc (DESIGN ONLY; member built later, gated on the benchmark)

> Scope of this cycle: a VALIDATED SIMULATION BENCHMARK + this design doc. No member code ships this cycle. Every primitive below is a **hypothesis pending ablation on the benchmark**, never a finding. The framing metric for every ablation is uniform: **exact indel-position concordance with truth, ambiguity-aware, never edit distance** (edit distance ties by construction at the contested positions).

---

## 1. Architecture: ONE native member, three layers, one substrate

The five grounded concepts are **not five members**. They are facets of a **single native member** that plugs into the retained 5-aligner consensus fan-out as an orthogonal arbitration/refinement layer. The member has three layers over one DP substrate:

```
                 ┌─────────────────────────────────────────────────────────────┐
  LOCALIZATION   │  inherited panel-union window (SETTLED)                       │
                 │     └─ gated FracMinHash containment fallback  ── C5          │
                 │        (fires ONLY in the all-5-unplaced complement set)      │
                 └──────────────────────────┬──────────────────────────────────┘
                                            │  proposes a reference window
                 ┌──────────────────────────▼──────────────────────────────────┐
  REFINEMENT     │  POA-pooled junction-homogeneous subpool consensus ── C4     │
  CORE           │     (UPSTREAM; emits a denoised exon-block per subpool)       │
                 │                          │                                    │
                 │  ┌───────────────────────▼─────────────────────────────────┐ │
                 │  │  ONE calibrated-NLL exon-block DP                        │ │
                 │  │    align_exon_block_global (local_aligner.py:522)        │ │
                 │  │    · per-position HP-length-law gap cost   ── C1         │ │
                 │  │    · forward-backward posterior on SAME lattice ── C3    │ │
                 │  │    · CPA change-point decoder shares length table ── C2  │ │
                 │  └─────────────────────────────────────────────────────────┘ │
                 └──────────────────────────┬──────────────────────────────────┘
                                            │  posteriors + runner-up per decision
                 ┌──────────────────────────▼──────────────────────────────────┐
  ARBITRATION    │  calibrated log-likelihood-RATIO consensus  ── C3            │
                 │    select.py:65 / chimeric_consensus.py:1135 integer-argmax  │
                 │    → LLR + abstain band; integer path retained as fenced     │
                 │      fallback; panel diversity preserved (no collapse)       │
                 └─────────────────────────────────────────────────────────────┘
```

**The unifying claim.** The 5-aligner panel (minimap2, gapmm2, uLTRA, GMAP, mapPacBio) shares ONE error family: a **flat-affine, quality-blind, run-length-agnostic** gap model that consumes neither per-base phred (`query_qualities` propagated but unconsumed: `chimeric_consensus.py:952`, `consensus.py:532`) nor the empirical per-length indel costs that already exist (`HpPenaltyTable.del_cost/ins_cost`, `hp_penalty.py:261/270`) but which the live recurrence never consults (`local_aligner.py:624-646`). Because the panel shares that bias, on HP/poly-A undercalls it **converges on the same wrong run length** — consensus is unanimous-but-wrong and cannot vote its way out. The native member's value is that it scores on a **different axis** (calibrated −log P from a held-out-validated length law + recalibrated Q), so its error family (mis-calibration) is statistically independent of the panel's shared flat-cost-tie family. This is simultaneously the **orthogonality source** and the **defense against the scoring-artifact family** that already bit the project (a pure-scoring re-weighting flipped GMAP's annotated win/loss 0.09→1.07 with no aligner change). The goal is to **add an orthogonal member and cut external deps 5→2-3 (keep minimap2)** — never replace the panel, never 5→0.

**The 87-90% anchor.** ~87-90% of contested per-segment decisions are INDEL/residue micro-arbitration at HP/repeat/poly-A boundaries — NOT splice biology. Each layer names its slice: C1 = the dominant HP/poly-A/STR micro-indel slice (the highest-leverage facet); C2 = the 3'/CPA terminus slice; C3 = the cross-decision arbitration scale (the meta-lever); C4 = the pooled-denoising slice; C5 = the residual ~12% panel-unplaced tail (localization only, then hands the DP the window so the standard concordance metric applies).

---

## 2. The shared substrate: the calibrated-NLL exon-block DP

Today `align_exon_block_global` (`local_aligner.py:522-688`) is a Gotoh three-matrix (H/D/I) affine-gap NW that **maximizes** a flat score (`_MATCH=2,_MISMATCH=-4,_GAP_OPEN=-4,_GAP_EXTEND=-1` at `:73-76`). Its only HP-awareness is the substitution branch: a binary `homo_mask` swaps `_MISMATCH=-4 → homo_mismatch=-2` at ref positions inside a run ≥ `min_run` (`:601-606`). The **gap** recurrence (D at `:624-634`, I at `:636-646`) is run-length-AGNOSTIC: a k-base deletion costs `gap_open + k*gap_extend` whether it sits in a 12-base poly-A or in unique sequence.

The member re-expresses this single recurrence in **NLL / minimization space**: every op contributes `−log P` from a calibrated emission, so match, sub, and indel|run sit on ONE validated scale (discharging the meta-lever). Scale/sign is load-bearing and NOT interchangeable: `align_exon_block_global` maximizes affine scores, while `del_cost/ins_cost` are POSITIVE costs on the minimization scale used by the numba HP DP (`hp_penalty.py:69`, substitution=1.0). Concretely:

- match → `−log P(match | run, base)`
- sub → `−log P(sub | run, base)` (subsuming the binary `homo_mismatch` knob — it becomes an ablation arm, not a co-resident heuristic)
- deletion ending a true-length-L run → `−log P(obs=L−k | true=L, base)`
- insertion → `−log P(obs=L+k | true=L, base)`

read off the empirical `HpPenaltyTable` (already a −log-likelihood-style penalty, populated by `scripts/calibration/empirical_cigar_error_profiler.py`).

**The substrate has two tiers** that define the member's MVP→full fork (C1) and host every other facet:

- **Tier A (per-position gap cost, MVP):** keep the Gotoh H/D/I structure; make the gap term position-dependent. Reuse `_precompute_del_costs(...)` (`hp_penalty.py:315-371`), which already emits a per-ref-position del-cost vector keyed on run length (and STR via `_str_repeat_info`); add a symmetric `_precompute_ins_costs`. Inject `−del_cost_vec[j-1]` at `:628/:633` and `−ins_cost_vec[i-1]` at `:640/:645`. Gap-OPEN stays a fixed prior (open/extend separation preserved). Minimal, contained edit.
- **Tier B (run-unit length-transition emission, FULL):** run-segment the ref once, then a within-run length-DP scores an entire HP/STR run as ONE `−log P(obs_run | true_run, base)` (correctly modelling A12→A10 as a single length error, not two iid deletions) and emits a posterior over candidate run lengths.

**Boundary invariance (load-bearing for FDR scoping):** `align_exon_block_global` is global (both spans fully consumed); the caller hard-asserts span preservation (`read_edits.py:818-828`) and exon-block boundaries (the adjacent N-op/junction edges) are FIXED by block segmentation (`read_edits.py:733-784`, op in {3,4,5} ends a block). The length-law therefore **re-attributes indel bases WITHIN a block** and can never create, delete, or shift a junction — so the C1 facet is structurally discovery-neutral (its FDR axis is "false indels", not "false junctions").

### Substrate integration points (C1 facet)
| file:line | change |
|---|---|
| `local_aligner.py:522-688` (`align_exon_block_global`) | add `penalty_table: Optional[HpPenaltyTable]=None`. None ⇒ byte-identical flat-affine (backward-compat fence). Provided ⇒ D/I gap recurrence (`:624-646`) consumes per-position cost vectors; H/D/I re-expressed on the NLL/min scale. |
| `local_aligner.py:73-76` + `:601-606` | flat `_GAP_OPEN/_GAP_EXTEND` → fixed open-prior + table-driven extend; the binary `homo_mismatch=-2` branch subsumed by `−log P(sub\|run,base)` — the law must **beat, not tie** the existing knob. |
| `hp_penalty.py:315-371` (`_precompute_del_costs`) | reused as-is for the MVP per-position del-cost vector; add symmetric `_precompute_ins_costs` for the I state. No new estimation code. |
| `read_edits.py:811-816` (`realign_exon_blocks`) | thread `penalty_table` through, sourced from `correct_command.py:707` (already imports `HpPenaltyTable`) or per-read via `PenaltyTableManager.for_read` (`junction_refiner.py:349`). |
| `hp_penalty.py:42-80` (`_score_hp_dp_numba`) | pattern to mirror for the numba-fast path of the new recurrence (`njit, cache=True`). **MUST pin AVX2/scalar (`NUMBA_CPU_FEATURES`) to avoid the AVX-512 SIGILL build trap.** |

---

## 3. Facet C2 — channel-aware 3'/poly-A CPA change-point decoder

**Slice:** the 3'/CPA terminus (a joint **localization+refinement** facet — you cannot inherit-then-refine here because the inherited END is itself the drifted quantity). Today the CPA is decided by a deterministic walkback that BREAKS at the first read==ref non-stop base (`walkback.py:178-186`), delegated via `find_polya_boundary` (`indel_corrector.py:98-164`); the A-tract detector emits only a fixed-shift window from `SHIFT_CORRECTIONS_BY_ACOUNT` (`config.py:92-104`); and `scoring.py:684-732` deliberately does NOT score the 3' terminus at all.

**Mechanism:** reframe the cut as a **change-point** over the inherited candidate window (deepest panel raw end inward to `find_atract_boundaries`' `first_non_a`, `atract_detector.py:101-206`). Each candidate cut k partitions the read into a **BODY** (templated genomic) and **TAIL** (untemplated poly-A) segment. For each k compute a log-likelihood: BODY = templated emission log-probs (match/mismatch as calibrated log-ratios + per-length indel costs from the SHARED `HpPenaltyTable`) **plus a per-base phred term** (the Q finally consumed); TAIL = `log P(base | poly-A process)` from `PolyAModel.non_a_frequencies` (`polya_model.py:88-90`). The two are tied by an **A-run length law** P(templated genomic-A-run length | downstream genomic-A count g), with g from `_count_downstream_as` (`atract_detector.py:49-98`). Return `argmax_k` posterior as `corrected_pos` AND the full posterior (credible interval over k) for consensus arbitration.

**Cross-facet wiring:** C2 **shares the length table** with C1/C3 (BODY indel log-likelihoods) and the Q-recalibration table with C3/C4. It consumes the same calibration split (§7).

### Integration points (C2)
| file:line | change |
|---|---|
| `walkback.py:178-186` (hard break-loop) | the primitive being replaced; retained as MVP/fallback comparator. Decoder scores ALL k, returns argmax + posterior. |
| `indel_corrector.py:1622` + `:98-164` (`find_polya_boundary`) | opt-in decoder mode returning `corrected_pos + cpa_posterior`; inherits the guarded walkback's artifact clipping (large-del/intron/tail-context guards, `walkback.py:274-300`) to build the candidate window. |
| `atract_detector.py:209-331` + `config.py:92-104` | fixed-shift point estimate → run-length-law MAP/posterior-mean; `has_ambiguity/ambiguity_range` → posterior credible interval. |
| `scoring.py:684-732` + `consensus.py:532` | decoded CPA posterior becomes the calibrated, aligner-agnostic terminus signal (meta-lever); `query_qualities` consumed via the recalibration map. |
| `hp_penalty.py:151-154,261-270` | `del_cost/ins_cost` supply BODY-state indel log-likelihoods at the HP/poly-A boundary. |

**CPA-specific FDR (first-class, separate from junction FDR):** before accepting any cut that MOVES the terminus, screen the downstream window with the existing weighted AG-mispriming score (`ag_mispriming.py:43-70`, threshold>15) — a decoded site over an A/G-rich stretch is flagged likely internal-priming, posterior penalized / deferred to fallback. The g=0 stratum returns EXACTLY the walkback anchor (anti-drift fence, §6).

---

## 4. Facet C3 — posterior-emitting refinement + LLR consensus arbitration (the meta-lever)

**Slice:** the cross-decision arbitration scale. This facet is what makes the member "calibrated-likelihoods-and-posteriors" rather than flat costs, and it is the **structural defense against the scoring-artifact family**.

**(A) Emit (the refiner).** Add a sibling **forward-backward (sum-product)** pass over the SAME bounded H/D/I lattice already allocated by `align_exon_block_global` / `score_left_anchored` (`:571-577`). Because the window is the inherited cluster (tens-to-low-hundreds of ref bp), this is O(Q·R) on cells already allocated — cheap CPU, AVX2/scalar. Emit, per contested decision, a **posterior** (normalized marginal of that op-at-that-position over all tracebacks) + the **runner-up** op/position and its posterior. Scores are interpreted as log-likelihoods on the SAME calibrated scale as C1 (recalibrated Q + `HpPenaltyTable`). The posterior is a temperature-scaled log-sum-exp; temperature fit on held-out truth.

**(B) Arbitrate (whole-read).** `select_best_alignment` (`select.py:65-91`) computes integer `max(junction_score)` + equality tiebreak. Replace with a calibrated **LLR** between the top-2 candidates: `LLR = Σ_contested [logL(A) − logL(B)]`. Winner decided by whether `|LLR|` exceeds a calibrated **abstain band**; inside the band → reported tie/low-confidence, not force-picked. **INERT when integer and LLR agree** (the common case), so it cannot regress clean reads.

**(C) Per-segment (chimeric).** `select_best_chimeric` picks each segment by `score > best_score` (`chimeric_consensus.py:1135`) over `SegmentScore` (`:173-194`). Attach `.posterior/.runner_up_aligner/.runner_up_posterior`; arbitrate by LLR at the UNCHANGED sync points (`find_sync_points :317`).

**(D) Propagate.** Surface `best_posterior/runner_up/abstained` on `ConsensusResult` (`extract.py:77,107-118`) so downstream FDR + NIC/NNC discovery consume a calibrated number.

### Integration points (C3)
| file:line | change |
|---|---|
| `local_aligner.py:522-620` + `:410-470` | sibling forward-backward over the SAME H/D/I lattice → per-position marginal + runner-up; gated behind a flag so Viterbi CIGAR is unchanged when off. AVX2/scalar. |
| `local_aligner.py:73-76` | flat constants become the FALLBACK when no calibrated table is loaded; emission path reads calibrated per-base log-likelihoods (Q + `HpPenaltyTable`), integer path bit-identical for regression. |
| `select.py:65-91` | integer-argmax + tiebreak → calibrated-LLR top-2 + abstain band; integer path retained as fallback and inert baseline. |
| `chimeric_consensus.py:1119-1136` + `:173-194` | add posterior/runner-up to `SegmentScore`; LLR arbitration at unchanged sync points. |
| `extract.py:77,107-118` | add `best_posterior/runner_up/abstained` to `ConsensusResult`. |
| `chimeric_consensus.py:952` + `consensus.py:532` | wire (recalibrated) Q into the emission calibration. |

**The single highest cross-concept regression risk — GMAP fence vs LLR.** `test_gmap_fence_regression.py` documents that GMAP earns its panel seat ONLY because the integer `score_segment` scoring suppresses its 878k spurious non-canonical junctions (55.6% of its calls) while keeping its ~100 genuine novel junctions + annotated parity. C3 replaces integer-argmax at exactly that selection path. **The doc's hard requirement: the new non-canonical FDR guard (§8) must demonstrably subsume Fence-1's suppression role BEFORE the integer path is retired**, or GMAP's noise leaks back through the new arbitrator. Resolution adopted here: **`score_segment` and its three fences are RETAINED as the fenced fallback** (the integer arbitration path stays behind a flag, §6 inertness fence pins it byte-identical). The LLR arbitrator is gated ON only after Phase-2 ablation shows LLR's non-canonical FDR (measured against NNC truth) is ≤ the integer path's at matched panel placements. Until then the fences are re-expressed as LLR-level invariants (a non-canonical de-novo placement must score strictly worse than a clean canonical one on the calibrated scale, and a short-anchor forced-canonical placement earns no posterior bonus) and BOTH the integer fences and the LLR fences run in CI.

---

## 5. Facet C4 — POA-pooled junction-homogeneous subpool refinement (UPSTREAM of the DP)

**Slice:** the same 87-90% HP/indel residue, but the refinement UNIT becomes a per-locus read POOL — pooled ONLY within a junction-homogeneous subpool, so majority denoising cancels stochastic per-read HP miscalls and **never votes a junction out** (discovery-suppression designed out, not mitigated later).

**Pipeline (AVX2/scalar, KB-scale C lib only):**
- **(a) BIN** = the inherited panel-union cluster. Reuse `build_junction_pool` (`junction_scoring.py:489`) + `collect_junction_counts_from_bam` (`:473`) for per-junction `{raw_count, anchored_count, read_ids}` membership — no new clusterer.
- **(b) SUBPOOL** by normalized junction-set signature (`normalize_junction :59`, `junction_ambiguity_window :78`, `_normalized_annotation_set :136`). A junction 10% of the bin uses forms its OWN subpool and survives.
- **(c) POA** per subpool, per exon block (abPOA/spoa) → denoised block consensus.
- **(d) ALIGN** the consensus ONCE to genome via the SHARED `align_exon_block_global` (`:522`).
- **(e) PROJECT-BACK** = explicit composition `(read→consensus) ∘ (consensus→genome)`: the POA per-read traceback composed with the single genome CIGAR. Indel POSITIONS inherited from the majority-denoised consensus; read-specific true differences (real SNV/minority indel) preserved as divergent columns. This composition is the buildable core (the part easiest to hand-wave).
- **(f) EMIT a POSTERIOR** (subpool depth, column entropy at the contested HP boundary, fraction agreeing) — ONE candidate per read arbitrated ALONGSIDE the panel + the native refiner via C3, not a replacement.

**DISPERSION GATE (first-class):** intra-subpool entropy/paralog-divergence. If a subpool exceeds a calibrated dispersion threshold (merged paralogs / mixed isoforms), SPLIT and recurse; if it cannot be cleanly split, **abstain to per-read placement** for those reads.

**Cross-facet placement:** C4 sits UPSTREAM of the §2 DP (feeds it a pooled exon-block) and feeds its posterior INTO the §4 arbitration. It is the only facet with a new hard external dep.

### Integration points (C4)
| file:line | change |
|---|---|
| `junction_scoring.py:473` + `:489` | reused read-only for the bin + per-junction `read_ids`. No edit to the pool builder. |
| `chimeric_consensus.py:59,78,136` | reused to compute the junction-set signature (subpool key). No edit. |
| `local_aligner.py:522` | called ONCE per subpool exon block (consensus→genome). No recurrence/constant change. |
| **NEW** `rectify/core/consensus/poa_pool.py` | bin→subpool partition; POA invocation (abPOA/spoa binding); per-read traceback capture; project-back composition → per-read CIGAR + confidence; dispersion gate. Only net-new file in the MVP. |
| `consensus.py:685` (`run_consensus_selection`) | FULL only: a PRECEDING locus-grouped pass produces a `read_id→(projected_cigar, confidence)` sidecar; `select_best_alignment` weighs it (not replaced). |
| `chimeric_consensus.py:952` + `consensus.py:532` | FULL: per-base Q → POA column confidence (ablated separately). |

---

## 6. Facet C5 — FracMinHash containment fallback localizer (gated, panel-failure tail)

**Slice:** the residual ~12% with NO acceptable panel window — all-5-unmapped/GMAP-absent or all-5-misplaced (empty/sub-threshold union). This is the ONLY place the independent fallback localizer earns its dep, **gated behind a measured depletion trigger**. It is NOT a 6th placer: it proposes ONLY a reference window and hands off to the existing DP.

**Mechanism:** FracMinHash sketch (hash every k-mer, keep hash < H_max/s); `Containment(read in W) = |sketch(read) ∩ sketch(W)| / |sketch(read)|` — an unbiased containment estimate from surviving k-mers, **no contiguous chain required**, which is why it survives where minimizer chains vanish at 5-10% ONT error (the panel's CORRELATED failure: all 5 share the seed-and-chain assumption). Candidate refs pre-sketched per fixed-stride tile; an integer inverted index (`sketch-hash → tile-ids`) yields top-containment windows in one pass. The window (with a small splice/indel pad) feeds the EXISTING DP (`align_exon_block_global :522` / `align_clip_to_exon :716`). k, s, and the accept threshold are EXPLICIT ablatable hypotheses (error-free k-mer fraction ≈ (1−e)^k: at e=0.07, k=15→~0.34, k=21→~0.22), co-tuned on the held-out split, never on internal score.

**Two sub-cases:**
- **(a) all-5-unmapped:** these reads are ABSENT from every iterator (`_filtered_read_iterator` drops unmapped/secondary/supplementary, `consensus.py:210-214`), so the `if not aligner_reads` guard (`chimeric_consensus.py:1033-1034`) is DEAD for this path. Hook = a **NEW separate fallback pass after the main batch loop** (`consensus.py` post-`:1005-1012` flush) whose universe = (FASTQ/sidecar read-ID set) MINUS (placed IDs). Re-source the raw sequence from the unfiltered BAM / sidecar / FASTQ. **Ships standalone.**
- **(b) all-5-misplaced:** read DOES enter the loop with a non-empty group; trigger = uniformly-low POST-REFINEMENT likelihood. **CROSS-CONCEPT DEPENDENCY: gated on C3's calibrated posterior existing** (it supplies the trigger), so it is a later phase.

### Integration points (C5)
| file:line | change |
|---|---|
| `consensus.py:210-214` + `:446-475` | DIAGNOSIS only (the all-unmapped tail manifests as ABSENT reads). No edit. |
| `consensus.py:685` post-loop | NEW fallback pass for sub-case (a): complement-set → sketch → containment → top window → existing DP → emit primary (or leave unplaced). |
| `consensus.py:549` + `chimeric_consensus.py:999-1034` | sub-case (b): containment re-localization adjacent to `select_best_chimeric`, gated on C3. |
| `local_aligner.py:522` + `:716` | REUSED unchanged as the hand-off target. |
| `junction_scoring.py:473-485` + `chimeric_consensus.py:59-155` | reused by the FDR guard to normalize/ambiguity-window fallback-localized junctions, counted separately from panel-corroborated ones. |
| **NEW** `rectify/core/align/containment_localizer.py` | pure-Python (MVP) then integer-inverted-index FracMinHash. AVX2/scalar, NO AVX-512. |

**Dep commitment is GATED:** the index dep is adopted ONLY if Phase-1 measures the tail large enough to matter (recommend skip if tail < ~1-2%).

---

## 7. The shared calibration regime (ONE section — all facets share it)

The `HpPenaltyTable` length law is consumed by C1, C2, C3, C4; the Q-recalibration table by C2, C3, C4. The held-out / maintenance discipline is stated ONCE here, not repeated per facet:

1. **ONE genomic-region-disjoint train/test split, shared across ALL facets.** A locus or region used to estimate ANY calibration table for ANY facet is in the TRAIN partition; ALL ablation scoring is on the disjoint TEST partition. The split must be the SAME across concepts — otherwise concept A's test region is concept B's train region and the win leaks across concepts even when each is internally clean.
2. **Per-chemistry/basecaller re-estimation + provenance stamp.** Every table (length law, A-run law, Q-recalibration) is a SILENT MAINTENANCE DEP: it must be re-estimated per chemistry/basecaller version, each carrying a version stamp; a stale table degrades the law silently toward flat. Cost this as a recurring recalibration, not a one-time fit.
3. **`min_count=100` sparsity tail (`hp_penalty.py:184,213`).** Below the gate the table row is dropped and the law falls back to flat. The benchmark MUST size every `(run_length, base_class[AT/CG])` cell and every STR `(unit, n_copies)` cell to clear the gate, OR the law silently nullifies and produces a FALSE REFUTE. The long-HP tail (≥8-10) and rare STR units may never clear even a large benchmark → a smoothed/parametric length-law prior for the tail is a named (tiny, KB-scale) extension.
4. **Free parameters tuned on TRAIN, frozen, never on TEST/truth:** the NLL gap-open prior (C1), the temperature + Q-recalibration (C3), the A-run length law + Q-recal (C2), the dispersion/confidence thresholds (C4), and k/s/threshold (C5).

**Dependency ledger (5→2-3).** External ALIGNERS drop (the member covers retired members' contributions — see open decisions). The member internally adds only: **abPOA/spoa** (C4, the one hard new external C lib, AVX2/scalar) and a **conditional inverted index** (C5, adopted only if the tail justifies it). Everything else — `HpPenaltyTable`, `_precompute_del_costs`, the numba pattern, the forward-backward over the existing lattice, the FracMinHash sketch (reimplementable dependency-free) — is in-tree or pure-Python. **No torch, no basecaller, no large weights. AVX2/scalar + `NUMBA_CPU_FEATURES` pin everywhere.**

---

## 8. Discovery / FDR strategy (CO-PRIMARY, first-class — one unified strategy)

De-novo NON-canonical junction discovery is a CO-PRIMARY deliverable alongside indel/CPA refinement. It RAISES FDR risk and is ANTI-CORRELATED with annotation reward, so the guard is a first-class design element. C3, C4, and C5 each touch discovery; their three guards are UNIFIED into one strategy:

1. **Separate canonical/non-canonical FDR tracks — never one pooled number.** The calibrated LLR (and every discovery threshold) is estimated SEPARATELY for canonical (GT-AG within the ambiguity window, via `_canonical_within_window` / `_normalized_annotation_set`, `chimeric_consensus.py:113-155`; annotated set via `load_annotated_junctions`, `consensus.py:1222`) vs non-canonical placements. Two FDR curves are always reported. A non-canonical de-novo win must clear a **non-canonical-specific** threshold set to hold non-canonical FDR at a pre-registered level on held-out TEST against NIC/NNC truth — NOT the canonical threshold (canonical placements get free annotation/motif support the de-novo ones must not borrow).
2. **Anti-minisplice resolution (the head-on tension): GT-AG enters as a SOFT posterior/likelihood term, NEVER a hard filter — the SAME rule in C3, C4, and C5.** Baking a GT-AG/annotated prior into the per-base emission would let a learned canonical prior fight de-novo discovery. So motif/annotation stays a downstream TIEBREAKER (as today: `scoring.py:680-682` / `select.py:79-89` / `SegmentScore.n_novel_canonical_junctions :190`), recoverable by strong read/subpool support. A genuine novel non-canonical junction with strong calibrated likelihood is NOT suppressed for being absent from the catalog; a short-anchor forced-canonical placement (`n_annotated_unsupported :193`) earns no posterior bonus. C4 reinforces this structurally: discovery support is required from a junction-homogeneous SUBPOOL (reads independently sharing the novel N-op), never from a consensus majority — a minority-subpool junction is never voted out because subpools never share a POA graph.
3. **The abstain band is the shared safety valve.** Low-LLR / low-support de-novo calls are reported as ABSTAINED, not asserted — capping false discovery across all three facets.
4. **Fallback-only junctions (C5, zero panel corroboration) clear a STRICTLY higher bar:** multiple agreeing fallback reads at the same ambiguity-normalized junction AND a minimum containment quality gate; admitted only through the COMPASS 3-pass `junction_validator.py` path with canonical/non-canonical tracked separately. The canonical/non-canonical split may force suppression of non-canonical FALLBACK calls entirely (capping that discovery upside) — an accepted, measured outcome, not a silent one.
5. **C1 is structurally excluded from junction FDR** (§2 boundary invariance): it cannot create/move a junction, so its FDR axis is false-indels only.
6. **The benchmark MUST carry NIC/NNC truth labels from Phase 0** or none of this is measurable.

**Subsumption requirement (re-stated):** before the integer `score_segment` path is retired, the non-canonical FDR guard above must demonstrably hold GMAP's non-canonical noise at or below the level the integer fences achieve, on held-out TEST. Until then, integer fences + LLR fences both run (§4).

---

## 9. Benchmark coupling (the integrating artifact — gates everything)

The simulation benchmark is the union of all five `benchmark_requirement_consumed` blocks; it is the system of record for fitness (NEVER the internal score). It MUST CONTAIN / MEASURE:

**Per-read truth & metric machinery (all facets):**
- exact-truth indel/junction POSITIONS at HP/poly-A/repeat boundaries, with the **precomputed ambiguity-equivalence set** of position-tied CIGARs per locus, so the harness measures indel-COUNT concordance (and indel-in-run vs mismatch-after-run attribution) without charging uncharged within-run shifts. **Framing metric = exact indel-position concordance, ambiguity-aware, never edit distance** (a call one bp into a donor/acceptor/poly-A repeat is NOT an FP).
- per-read which-aligners-placed labels (so "panel-unplaced" and "panel placement at each contested decision" are decidable).
- the **single shared genomic-region-disjoint TRAIN/TEST split tag** (§7).
- NIC/NNC junction truth labels (§8).

**C1 (HP-length law):** per HP/poly-A/STR locus — true run length, an error-injected observed run from a realistic DRS length-error model, a canonical left-aligned truth CIGAR, AND **clean (no-error) runs** to measure the false-indel rate. Every `(run_length, base_class)` and STR `(unit, n_copies)` cell sized to clear `min_count=100`.

**C2 (CPA):** per-read ground-truth CPA position + the stratification context — downstream genomic-A count g, abutting-A-tract length/boundary — with the **tail-over-genomic-A-tract class OVERSAMPLED** (the drift-prone regime). Emit per read the 5 panel raw ends + walkback + atract estimates so all comparators score identical reads.

**C3 (LLR):** ≥1 constructed family where two aligners TIE under flat-affine integer scoring but DIFFER under truth (to exercise LLR vs integer); metadata to REPLAY a cost re-weighting (the 0.09→1.07 scenario) at fixed placements.

**C4 (POA):** PARALOG loci (≥2 near-identical copies), ISOFORM-MIXTURE loci (major + real minor junction sharing a window), **per-read READ-OF-ORIGIN labels** (so mis-clustering-induced errors are MEASURED not assumed), and per-locus COVERAGE (for the ≥5-read gate + stratification).

**C5 (fallback):** a measurable, labeled population of PANEL-UNPLACED reads (empty union OR all-misplaced), each RETAINING ground-truth origin window + indel positions. This population will NOT appear at low simulated error; the benchmark must GUARANTEE it (inject a tunable sub-population at elevated error / in repeat / low-complexity / novel-junction contexts) and MEASURE its size — the tail size is the C5 dep-commit decision gate.

**The benchmark is the gating first step of the roadmap.** No member facet is built before its `benchmark_requirement_consumed` is met and its ablation is runnable.

---

## 10. Phased roadmap (sequenced by the dependency DAG, then leverage × risk)

The hard ordering constraint: **C2 sub-case (b) and C5 sub-case (b) are gated on C3's posterior existing.** That drives the sequence. (Exit criteria are pulled from each concept's `build_phases`, not re-derived.)

1. **Benchmark** — the gate; co-deliverable of the cycle. Exit: union of all five requirement blocks met, NIC/NNC + split tags present, every C1 cell clears `min_count`, panel-unplaced population guaranteed + measured.
2. **C1 Phase 0 (equivalence harness)** — `penalty_table=None` byte-identical fence + CIGAR left-aligner + ambiguity-set generator. Exit: 0 CIGAR diffs over the corpus.
3. **C1 Phase 1 (per-position gap cost MVP)** — highest leverage (the 5-panel-CONVERGENT bias, the dominant slice), lowest risk (contained D/I edit). Exit: arm-law STRICTLY beats BOTH arm-flat AND arm-knob on held-out position-exact concordance; REFUTED if it only ties arm-knob or the win vanishes on train/test swap.
4. **C3 Phase 0-1 (emit + calibrate)** — forward-backward on C1's lattice (cheap, reuses cells); the meta-lever; unblocks 2b/5b. Exit: marginal matches brute-force enumeration to 1e-6 with flag off byte-identical; on TEST posterior calibrated (ECE under threshold) and beats uncalibrated softmax.
5. **C2 (CPA decoder)** — shares the length table; narrower slice. Exit: MVP reduces median |est−true CPA| AND into-tract drift vs walkback with NO g=0 regression; FULL strictly dominates MVP on held-out with calibrated posterior coverage.
6. **C3 Phase 2-3 (LLR arbitrate + artifact-replay)** — Exit: LLR concordance ≥ integer at matched placements with verified inertness on the clean subset; the 0.09→1.07 re-weighting flips the integer winner but NOT the LLR winner; non-canonical FDR ≤ integer path (subsumption of Fence-1) before integer retirement.
7. **C5 MVP + Phase 1 (kill-gate, then size the tail)** — Exit: containment beats the random-window null on held-out panel-unplaced reads (else KILL before any index); then commit the index dep ONLY if measured tail ≥ pre-set fraction AND recovered-read concordance ≥ panel baseline.
8. **C4 (POA)** — last of the core; only hard external dep + new failure family (mis-clustering). Exit: pooled-then-projected beats per-read on HP-boundary position concordance at ≥5-read loci; mis-clustering errors INTRODUCED stay strictly below indel errors REMOVED; minority-isoform reads never collapsed.
9. **C2b / C5b (posterior-gated)** — built last, gated on C3's posterior. Exit: each recovers additional correctly-localized/refined reads at controlled (separate canonical/non-canonical) FDR.

---

## 11. Regression-guard plan (extends the `test_gmap_fence_regression.py` discipline)

The existing fence test pins **behavioral invariants (relative orderings), NOT exact constants** — its header literally says "if any of these orderings inverts… Break one ⇒ revisit the decision, do not just update the assertion." Every new guard inherits that discipline. **Critically, to prevent a calibration table from being overfit INTO the guards, every golden-set guard pins a RELATIVE invariant on held-out loci (e.g. arm-law concordance > arm-knob > arm-flat), NEVER an absolute concordance number** — an absolute-number golden set IS the overfit trap (the table could be fit to reproduce it). Guards by facet:

- **Backward-compat / inertness fences (all facets):** `penalty_table=None` ⇒ byte-identical CIGARs (C1); calibration-flag-OFF ⇒ byte-identical CIGAR + selection output (C3); fallback-flag-OFF ⇒ no panel-placed read is ever re-localized/altered (C5). Fail loudly if a future change re-couples arbitration to the flat integer scale.
- **Span invariants:** `sum(M,I)==Q_block` and `sum(M,D)==R_block` for every emitted block (C1, mirrors `read_edits.py:819-821`); the same on every projected per-read CIGAR (C4).
- **Frozen HP-boundary golden set (C1):** pinned poly-A/HP loci whose calibrated concordance must never regress BELOW the established arm-law baseline — but expressed as the relative ordering arm-law ≥ arm-knob ≥ arm-flat on HELD-OUT loci, guarding against a future table re-estimation silently degrading the law without becoming an absolute number to fit.
- **Anti-drift fence (C2):** in the g=0 stratum the decoder returns EXACTLY the walkback anchor, and decoded cuts NEVER leave the inherited window into gene body.
- **Calibration/overfit fence (C2/C3):** the held-out offset/ECE is recomputed in CI and FAILS if it degrades vs the MVP baseline — pinning that the win is not train-split hill-climbing into the simulation truth set.
- **Membership-preserving fence (C4):** for a synthetic 90%-major / 10%-minor bin, project-back NEVER reassigns any minority read's junction to the majority coordinates (fences the discovery-suppression mode permanently).
- **GMAP-fence subsumption / artifact-replay fence (C3 — the highest cross-concept risk):** the three existing `score_segment` fences are RETAINED and keep running; ADD an LLR-level re-expression (non-canonical < canonical on the calibrated scale; short-anchor forced-canonical earns no bonus) and the Phase-3 replay fence (LLR winner INVARIANT under the 0.09→1.07 re-weighting at fixed placements; integer path reproduces the flip). Both suites in CI; the integer path is not retired until the LLR non-canonical FDR demonstrably subsumes Fence-1's suppression on held-out TEST.
- **No-fabrication fence (C5):** the fallback emits a record only when containment AND post-DP score clear threshold; else the read stays unplaced.

---

## 12. Risk register

(See structured `risk_register` field.)

---

## 13. Open decisions for the user

(See structured `open_decisions_for_user` field — chief among them: WHICH 2-3 aligners survive is benchmark-earned, not decidable now.)


---

# APPENDIX — Redteam punch list (full)

**Verdict:** sound_with_fixes

**Unsupported / mischaracterized file:lines found:**
- VERIFIED GROUNDING + doc §2/§C5: 'align_exon_block_global ... called in the live 5' rescue path (read_edits.py:811, junction_refiner.py:717, splice_aware_5prime.py:1930/2052)' — WRONG: align_exon_block_global (local_aligner.py:522) is called ONLY at read_edits.py:811; junction_refiner.py:717 and splice_aware_5prime.py:1930/2052 call align_clip_to_exon (local_aligner.py:716), a different DP with a different signature.
- §4 + §C3(D): 'extract.py:77,107-118' for ConsensusResult — extract.py:77 is AlignmentInfo.junction_score (class at :44), not ConsensusResult (class at :107). The :77 cite is mislabeled; :107 is correct.
- §2 + §11: 'caller hard-asserts span preservation (read_edits.py:818-828)' / 'mirrors read_edits.py:819-821' — :818-828 is a logger.warning + keep-original soft fallback, NOT a hard assert.
- §4/§3/§10/§11: 'integer max(junction_score) / integer-argmax / integer path' — junction_score is float (extract.py:77) and score_segment (chimeric_consensus.py:553-697) uses non-integer multipliers (ev.length * 0.5, * 3.0). Scores are flat float-valued, not integer-valued.

**Punch list:**

- **[BLOCKER / unsupported_file_line]** VERIFIED GROUNDING block + doc §2/§C5 ("the live 5' rescue path (read_edits.py:811, junction_refiner.py:717, splice_aware_5prime.py:1930/2052)")
  - issue: The grounding asserts align_exon_block_global is 'imported + called in the live 5' rescue path' at those four sites. Verified against source: align_exon_block_global (local_aligner.py:522) is called at exactly ONE site, read_edits.py:811 (the exon-block REALIGNER, not 5' rescue). junction_refiner.py:717 and splice_aware_5prime.py:1930/2052 import and call a DIFFERENT function, align_clip_to_exon (local_aligner.py:716), which has a different signature (takes intron_start/intron_end/strand) and different semantics (clip-to-exon, not global-block). The doc's central 'ONE DP substrate' unification (§2) and C5's hand-off claim ('align_exon_block_global :522 / align_clip_to_exon :716') treat these as interchangeable substrates, but they are two distinct live affine DPs. The 'single recurrence' framing must be re-scoped to whichever DP each facet actually extends.
  - fix: Correct the grounding: align_exon_block_global is called only at read_edits.py:811 (exon-block realignment). The 5' rescue DP is align_clip_to_exon (local_aligner.py:716), called at junction_refiner.py:717 and splice_aware_5prime.py:1930/2052. Decide explicitly whether the calibrated-NLL substrate extends ONE or BOTH DPs, and write per-DP integration points; do not present them as one recurrence.

- **[MAJOR / overclaim_vs_prior_art]** §4 Facet C3 (the meta-lever) and §1 'orthogonality source'
  - issue: C3 = 'forward-backward (sum-product) over the H/D/I lattice + calibrated -log P emissions + posterior/runner-up decoding + LLR arbitration' is precisely a pair-HMM (Durbin et al.), as instantiated in Dindel (Albers 2011) and GATK HaplotypeCaller's PairHMM. The doc frames calibrated-likelihood-DP-with-posteriors as THE novel orthogonality source and meta-lever with zero citation of pair-HMM/Dindel/GATK. NOTE: the orthogonality claim itself (a pair-HMM/mis-calibration error family is statistically independent of the panel's shared flat-affine seed-and-chain family) is SOUND - the issue is presenting a 15-year-old standard method as novel insight, which is exactly the overclaim a redteam exists to catch.
  - fix: Cite Dindel (Albers 2011), GATK HaplotypeCaller PairHMM, and Durbin pair-HMM. Reframe C3's contribution as 'apply the established pair-HMM forward-backward/posterior machinery, with an ONT-recalibrated emission, as an orthogonal member to a seed-and-chain panel' - the novelty is the integration/orthogonality, not the DP.

- **[MAJOR / overclaim_vs_prior_art]** §6 Facet C5 (FracMinHash containment fallback)
  - issue: The doc uses the literal coinages 'FracMinHash' (sourmash, Irber et al. 2022) and 'Containment(read in W) = |sketch(read) ∩ sketch(W)| / |sketch(read)|' (Mash containment, Ondov 2016; Koslicki & Brown 2019) and presents 'survives where minimizer chains vanish at 5-10% ONT error' as its own derivation. That containment-beats-chaining-under-high-error property is the well-known sourmash/Mash selling point, not a novel result, and the error-free k-mer fraction (1-e)^k is the standard Mash derivation. No citation of sourmash/Mash anywhere. Concept is sound and the gating discipline is good; the framing overclaims.
  - fix: Attribute FracMinHash to sourmash (Irber 2022) and containment estimation to Mash (Ondov 2016) / Koslicki-Brown. State the contribution as 'apply a known sketch-containment localizer to the panel-unplaced tail', and note the option to simply reuse sourmash rather than reimplement (costs the dep honestly).

- **[MAJOR / overclaim_vs_prior_art]** §5 Facet C4 (POA-pooled junction-homogeneous subpool)
  - issue: pool-by-cluster -> POA (abPOA/spoa) -> align-the-consensus-once-to-genome -> project back per-read is the established long-read isoform/consensus-correction recipe (RATTLE, isONcorrect, FLAIR's correction step). The doc presents the pipeline (a)-(f) as its own design with no acknowledgment of this prior art. The genuinely defensible novel elements are the junction-homogeneous SUBPOOL keying (discovery-non-suppression) and the project-back composition into RECTIFY's arbitration - those should be foregrounded.
  - fix: Cite RATTLE / isONcorrect / FLAIR-correct as the pooling-POA precedent and scope C4's novelty to the junction-homogeneous subpool key + the dispersion gate + feeding a posterior into the C3 arbitrator, not the pool->POA->align pattern itself.

- **[MAJOR / confabulation]** §4 (B)/(C), §3, §10, §11 - repeated 'integer max(junction_score) / integer-argmax / integer arbitration path / integer path retained as fenced fallback'
  - issue: The doc characterizes the current selection path as 'integer' ~10 times and builds the orthogonality framing on 'integer path vs calibrated LLR'. Verified false: select.py:65 maxes junction_score which is typed float (extract.py:77) and score_segment (chimeric_consensus.py:553-697) accumulates non-integer terms (ev.length * 0.5, * 3.0, etc. at :606/:614/...). The scores are flat float-valued, not integer-valued. 'Integer-argmax' is a mischaracterization; the substantively correct critique (which the doc itself half-concedes as 'flat-affine') is 'flat / quality-blind / discrete-tie-prone float scores'.
  - fix: Replace every 'integer' with 'flat float' (or 'flat-affine discrete-tie-prone'). The argument that the panel ties under the flat scale and LLR breaks the tie survives the correction; the 'integer' word does not and undercuts credibility.

- **[MINOR / overclaim_vs_prior_art]** §2 / §1 presentation of C1 (per-position HP-length gap cost)
  - issue: C1 (restructure the gap recurrence to look up per-ref-position del cost as a function of HP run length rather than per-gap-event) is presented as the architect's highest-leverage insight, but local_aligner.py:30-60 docstring ALREADY documents this exact restructuring as a known TODO: '...not directly interchangeable. If HP-specific gap costs are needed here, the DP would need to be restructured to look up del cost per ref position rather than per gap-event. See junction_refiner._score_hp_anchored for the per-position HP del cost model.' This is good (C1 is well-grounded, not confabulated) but the doc should credit the source's own note rather than imply de-novo ideation. Substance is sound.
  - fix: Add one line acknowledging local_aligner.py:30-60 already names the per-position-vs-per-gap restructuring and that _score_hp_anchored / _precompute_del_costs already implement the per-position model on the minimization scale; frame C1 as 'wire the documented restructuring into align_exon_block_global', which it is.

- **[MINOR / unsupported_file_line]** §4 integration table and §C3(D): 'extract.py:77,107-118' for adding best_posterior/runner_up/abstained to ConsensusResult
  - issue: extract.py:77 (junction_score: float) is a field of dataclass AlignmentInfo (class at :44), NOT ConsensusResult (class at :107). The doc conflates the two classes in a single citation. The actual target for the new posterior fields (ConsensusResult, :107) is correct; the :77 line is mislabeled.
  - fix: Cite ConsensusResult at extract.py:107 for the new fields; drop or relabel :77 (it is AlignmentInfo.junction_score).

- **[MINOR / unsupported_file_line]** §2 'caller hard-asserts span preservation (read_edits.py:818-828)' and §11 'mirrors read_edits.py:819-821'
  - issue: read_edits.py:818-828 is NOT a hard assert: it computes new_q/new_r and on mismatch emits logger.warning(...) then 'keeps original' (continue) - a soft check with graceful fallback, not an assertion that aborts. The boundary-invariance / FDR-scoping argument (C1 cannot create/move a junction) still holds because of block segmentation at op in {3,4,5} (read_edits.py:733-784, verified), but 'hard-asserts' overstates the runtime guarantee.
  - fix: Change 'hard-asserts span preservation' to 'checks span preservation and reverts to the original block on mismatch (warn + keep-original, read_edits.py:818-828)'. The discovery-neutrality argument is unaffected.

- **[MINOR / untestable_ablation]** §8.1 / §10 step 8 / §11 C4 membership fence: 'mis-clustering errors INTRODUCED stay strictly below indel errors REMOVED'
  - issue: This C4 exit criterion is stated as a net aggregate over 'errors' without pinning that 'errors' = exact-indel-position concordance. A net 'introduced < removed' count can hide position SHIFTS that cancel (a removed false-indel and an introduced position-shift net to zero) and risks bottoming out in indel-count/edit-distance rather than the mandated position-exact, ambiguity-aware metric. §9 C4 does specify 'HP-boundary position concordance', so the intent is right; the §8/§10 aggregate must inherit it explicitly.
  - fix: Restate the C4 exit criterion as: net change in EXACT-INDEL-POSITION concordance (ambiguity-aware) is strictly positive, AND no minority-isoform read's junction is reassigned - not a raw introduced-vs-removed error count.

- **[MINOR / fdr_risk]** §4 C3(A) emission uses recalibrated Q + §7 (ONT Q miscalibration) feeding the non-canonical discovery LLR in §8
  - issue: C3's posterior/LLR is the discovery signal for de-novo non-canonical junctions (§8), and its emission depends on recalibrated ONT Q. The doc flags ONT Q miscalibration and asks to ablate recalibration separately (good), but does not close the loop on the FDR consequence: a stale or regionally-biased Q-recalibration table (a silent maintenance dep, §7.2) directly inflates or deflates the non-canonical LLR and therefore the co-primary discovery FDR, with no stated guard tying Q-table staleness to the non-canonical FDR track specifically.
  - fix: Add to §8 an explicit dependency: the non-canonical FDR track must be re-validated whenever the Q-recalibration table version changes, and the per-chemistry provenance stamp (§7.2) must gate discovery output - i.e. a stale Q stamp forces the discovery LLR to abstain rather than assert.


# APPENDIX — Open decisions for user

1. WHICH 2-3 external aligners survive the 5->2-3 cut is NOT decidable now and is benchmark-EARNED, not assumed: the member must demonstrate on the benchmark that it covers each retired aligner's unique contribution (e.g. GMAP's ~100 genuine novel junctions + annotated parity per the fence test) before that aligner is dropped. minimap2 is kept by decree; the other 1-2 keepers are an output of the validation, not an input.

2. Abstain-band width and the non-canonical FDR target level are pre-registered knobs that trade discovery sensitivity against false-discovery: set them on the held-out split to a chosen FDR ceiling. The user must pick the operating point (e.g. control non-canonical FDR at 5% vs 1%) since it directly caps the co-primary discovery upside.

3. Should the integer score_segment path be RETIRED once LLR subsumes Fence-1, or kept permanently as a redundant fenced fallback? Retiring simplifies the code but removes a proven suppressor; keeping it costs a flag + dual-CI but is strictly safer. Recommend keep-as-fallback until multiple chemistries validate.

4. C5 index dependency commitment: confirm the tail-size threshold below which the FracMinHash index is NOT adopted (proposed ~1-2% of reads). This is a one-way-ish infra decision tied to the measured panel-unplaced fraction on the benchmark.

5. Per-read UMI-binned penalty tables (PenaltyTableManager.for_read) vs a single pooled table for the exon-block path: the per-read branch adds hot-loop cost; decide whether the binning win justifies it or default to pooled until measured.

6. Q-recalibration scope: ship Q-blind first and add the recalibration table only where its ablation shows non-negative lift (C2/C3/C4 each ablate it separately). Confirm this conservative default vs wiring Q everywhere up front given the per-basecaller maintenance burden.



---

# ADDENDUM — two user-proposed novelty adds (2026-06-18)

Two PI-proposed adds, assessed + grounded against the live tree. Both fit the architecture; both
remain hypotheses-pending-ablation on the benchmark.

## Add (a) — generalize the empirical error table into the emission model (STRENGTHENS C1/C3)
**Status: already the spine; the add is a generalization, not a new concept.** C1/C3's meta-lever already
wires `HpPenaltyTable` (`−log P` from `scripts/calibration/empirical_cigar_error_profiler.py`) into the
DP. The profiler currently keys emissions on `(op_type, ref_base, hp_run_length, strand)` — i.e. the
empirical model is **homopolymer/STR-conditioned only**. The add: broaden the conditioning context so the
WHOLE emission (match/mismatch/indel) is empirical, not just the in-run gap term — candidate contexts:
short k-mer / dinucleotide neighborhood, position-in-read (basecaller drift), and base-quality decile.
- **Honest novelty:** LOW as method (this is a context-richer pair-HMM emission — see C3 prior-art:
  Dindel/GATK). The orthogonality is the **ONT-chemistry-specific empirical table**, not the DP.
- **Integration:** extend `empirical_cigar_error_profiler` key tuple + `HpPenaltyTable` lookup signature;
  the C1 recurrence already consumes a per-position cost vector, so richer context is a table swap, not a
  recurrence change. Strict backward-compat: HP-only context = today's behavior (an ablation arm).
- **Guardrail:** more context cells ⇒ the `min_count=100/cell` benchmark floor bites harder (sparse
  high-order contexts silently nullify to flat → false refute). Add context cells ONLY where the cell
  count clears the floor on TRAIN; fall back to the marginal (HP-only) cost otherwise (hierarchical
  back-off). This is the FDR-safe way to add context.
- **Ablation:** richer-context emission must beat HP-only-context on position-exact (ambiguity-aware)
  indel concordance on held-out TEST; refuted if back-off-to-marginal is never beaten.

## Add (b) — variant / haplotype-aware alignment (NEW CONCEPT C6) — defensively load-bearing for discovery
**Status: partially already built.** RECTIFY already ships `run_variant_aware_scan` (`variant_scan.py`) →
`VariantAwareHomopolymerRescue` (`indel_corrector.py`): a single-pass position-level **variant-frequency
map** that distinguishes a basecaller HP error from a real variant (`min_variant_fraction=0.8`,
`min_reads_for_variant_call=5`). Today that map only gates **HP rescue**. C6 extends it.
- **Why it matters for CO-PRIMARY discovery (the load-bearing reason):** a standing germline/somatic
  variant (SNP or indel) near a splice site can shift the local alignment and **fabricate a "novel"
  junction**, OR a variant in the exon can be mis-charged as a mismatch that the DP "repairs" with a
  spurious micro-indel. Both directly **inflate the de-novo (esp. non-canonical) discovery FDR** that §8
  is built to control. A variant-aware reference reclassifies variant-induced pseudo-junctions BEFORE the
  discovery LLR sees them. C6 is therefore not just an accuracy add — it is a first-class **discovery-FDR
  guard**, complementary to §8's abstain band.
- **Two tiers (MVP → frontier):**
  - **C6-MVP (variant-aware emission, dependency-light):** reuse the existing variant-frequency scan to
    build a per-position alt-allele set; in the C1/C3 emission, a read base matching a **supported alt
    allele** is scored as a MATCH-to-variant (`−log P(match)`), not a mismatch — so the DP stops
    "repairing" real variants with micro-indels/junctions. The discovery track must mark any junction
    whose donor/acceptor falls within `±k` of a supported variant as **variant-adjacent → abstain unless
    independently corroborated** (ties into §8 + the Deliverable-B corroboration logic). No new external
    dep (in-house pileup-style scan already exists).
  - **C6-FULL (haplotype phasing + ploidy inference — the frontier, genuinely novel):** phase the
    supported variants into haplotypes and realign against the **allele-aware (diploid/polyploid)**
    reference. A549 is an **aneuploid** cancer line, so this needs **ploidy / copy-number-aware allele
    inference** (allele fractions are not 0.5/1.0) — this is the genuinely novel, research-grade piece and
    where the orthogonality is strongest (NO panel aligner is haplotype-aware). Gate FULL behind MVP
    proving variant-induced FDR is non-trivial on the benchmark.
- **Orthogonality:** all 5 panel aligners align to a flat haploid reference → variant-induced errors are a
  CORRELATED family they all share; a variant-aware member's error mode (variant-call FDR / mis-phasing)
  is independent.
- **Dependency cost:** MVP light (reuses in-house scan; no torch/no external caller — explicitly NOT
  clair3/deepvariant). FULL adds a phasing step (in-house or a light tool); ploidy inference is compute,
  not a heavy dep.
- **Benchmark requirement (MUST land in P0 — see SIMULATION_BENCHMARK_SPEC.md):** a **standing-variant
  stratum** — simulated reads carrying known SNPs/indels at known positions, het vs hom, near vs far from
  junctions/CPA, and at non-Mendelian VAFs (to mimic aneuploid A549) — so we can MEASURE (1) how much
  variant-induced junction/indel FDR the current flat-reference panel suffers, and (2) C6's reduction of
  it. Without this stratum C6 is unfalsifiable.
- **Ablation:** with the variant stratum, exact junction/indel FDR with vs without C6-MVP variant-aware
  emission; FULL must additionally reduce error at heterozygous/allele-specific loci over MVP. Refuted if
  variant-induced FDR is negligible on the benchmark (then C6 is descoped — a measured decision).

## Roadmap placement
- (a) folds into **C1 Phase 2 / C3** (table generalization with hierarchical back-off) — no new phase.
- (b) C6-MVP is a **new P-phase after C3 (LLR) and before/with C4**, because it depends on the calibrated
  emission (C3) and feeds the discovery-FDR track (§8). C6-FULL (phasing+ploidy) is **deferred to a later
  cycle**, gated on C6-MVP showing variant-induced FDR is material on the benchmark.


---

# INVESTIGATION INTEGRATION (2026-06-18) — see dev/ALIGNER_INVESTIGATION_SYNTHESIS.md

The cloud-agent aligner investigation (re-verified vs the build via CORRECTIONS) converges on this design
and adds critical refinements. Headline flip: **junction quality ALREADY drives production selection**
(HP-edit-distance Path A is live) → we are sharpening a live metric. Critical ADDs to fold in:
- **#1 HP-length-law DP MUST charge each N-op a CALIBRATED −logP(novel intron), not 0** (the live
  free-false-intron exploit, J1a). NOT a fixed penalty (that re-introduces a gate, violating the PERMANENT
  no-candidate-gate policy). + add a **positive splice-strength term** (yeast PWM/MaxEnt, swappable iface) —
  RECTIFY has none today (only canonical/non-canonical tiers).
- **#3 LLR/posterior arbitration MUST de-herd correlated aligners** — minimap2/gapmm2/uLTRA share lineage
  (~3× over-count); only deSALT + mapPacBio are independent → lineage/empirical-ρ weights. Add a calibrated
  per-read `selection_confidence` + `_abstain` flag.
- **#4 POA over-homogenization guards are MANDATORY** (per-read veto + tight max_boundary_shift +
  shared-boundary-anchoring) or it erases alt-SS isoforms <50bp apart; validate vs orthogonal truth only.
- **#2 CPA decoder scope DOWN to the realizable** (`pt:i` tail-length Bayesian prior + 3′-tight/5′-loose DP
  asymmetry + cDNA mispriming veto); full dwell/Remora is blocked on POD5 squiggle RECTIFY doesn't retain.
- **#5 FracMinHash** has no investigation analogue — genuinely novel, keep.
