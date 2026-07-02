# INTEGRATION REVIEW — DIRECTOR_ALGO_EVAL_SYNTHESIS.md
*Sonnet reviewer (agent a23805945182aa994). Task: what in the synthesis is already-in-plan vs
genuinely new, and how should it re-sequence the plan. Report-back only — no commits.*
*Date: 2026-07-01.*

---

## Q1. What's already in plan (confirmation, not new)

| Synthesis claim | Already-in-plan as |
|---|---|
| "Union-of-aligners floor control" (PANEL column) | `panel_blindspot_score.py` PANEL column — BUILT. In-flight as job 32422876. |
| "Motif-blind -logP realignment member" as leading concept | `dev/NOVEL_JUNCTION_BLINDSPOT.md` §Addressability + `dev/RECTIFY_STRATEGIC_FRAME.md` — already the named leading member. |
| "Run the real panel over the ladder" | JOB 32422876 IN FLIGHT (minimap2 + mapPacBio + deSALT + uLTRA on the NIC/NNC ladder, `run_panel_blindspot.sbatch`). |
| Decision threshold ~60% → build / ~80% → pivot | `dev/NOVEL_JUNCTION_BLINDSPOT.md` "Cross-model correction" section and `HANDOFF_ALIGNER_BENCHMARK.md` line 1406. |
| FDR guard for novel calls | `scripts/benchmark/novel_support_probe.py` (WS-3) — built, committed. |
| deSALT best novelty coverage; uLTRA annotation-snapping; gapmm2 near-redundant | `dev/NATIVE_ALIGNER_OVERVIEW.md` + handoff + committed drop of gapmm2 from human-DRS panel. |
| Isoform assemblers consume minimap2 BAMs → correlated → not discovery engines | `dev/NATIVE_ALIGNER_OVERVIEW.md` + SPEC. |
| Consensus items 1-11 broadly | These faithfully restate priors in NATIVE_ALIGNER_OVERVIEW + SPEC — the synthesis itself flags many as "low-surprise restatements." |

**The "union-of-aligners floor control" = PANEL column = already the exact column in `panel_blindspot_score.py`.** The synthesis is confirming the design, not adding a new column.

---

## Q2. Genuinely new and actionable (ranked)

### #1 — `_CANONICAL_HP_PRIOR` / `is_novel` CORRECTOR ablation [NEW, high-leverage]

**The finding:** `junction_refiner.py:647` applies a `_CANONICAL_HP_PRIOR=0.5` discount during
per-read CORRECTION (tier < 4: a non-canonical junction must outscore the best canonical
alternative by ~one HP=1 deletion unit), plus an `is_novel` tiebreaker
(`junction_refiner.py:638,656–663`). **This corrector-step canonical preference has never been
measured for its cost to non-canonical discovery.** Neither the HANDOFF's inline C-gate list, nor
DISCOVERY_TIEBREAK.md (which covers `select.py` arbiter tiebreaker — a different code path), nor
any planning doc proposes this ablation. The HANDOFF at line 1416-1418 explicitly flags it as "a
NEW verified finding" surfaced BY the synthesis.

**Why distinct from the panel blindspot job (job 32422876):** Job 32422876 measures whether any
aligner **generates** the true novel site (aligner-generation failure). The corrector ablation
measures whether the corrector **discards correctly-placed** non-canonical junctions (corrector-
flattening failure). A read where deSALT placed the AT-AC junction correctly can still be
re-snapped to canonical by the corrector. Same corpus, different failure mode, different code path.

**Reuses existing harness:** same NIC/NNC ladder corpus. Needs a new probe script (~50 lines) that
runs the shipped `_junction_refiner` on per-rung aligner output with `canonical_discount=0.5`
(ON) vs `0.0` (OFF) and reports corrector-step recall per rung. M1-local, fast.

**Gate outcome if cost > ~5% recall on AT-AC rungs:** ship the evidence-scaled variant (reads the
discount off -logP tables) immediately — a zero-new-member win. If < 2%: confirmed precision knob,
no change.

### #2 — cDNA-UMI-stratified penalty tables as per-read cost signal IN the realigner [NEW design spec]

**The finding:** `penalty_scores_cdna_umi{1,2,3plus}.tsv` exist (verified). These are per-UMI-
family-size cost tables — exactly the per-read reliability signal that could weight the motif-blind
realigner's DP by read redundancy (singleton reads down-weighted vs umi3plus). Neither the
NOVEL_JUNCTION_BLINDSPOT.md nor the STRATEGIC_FRAME mentions per-UMI-bin costs as a design input
to the realigner; WS-3's novel_support_probe.py uses a Binomial p-value, not UMI-bin tables.

**Not a standalone experiment:** this is a design requirement to fold into the motif-blind
realigner design (step 2 of the build sequence). It doesn't need its own harness; it reuses
the same corpus. Rank it below the double-prior ablation because it depends on the realigner being
built.

### #3 — Chimera/trans-splicing stratum as a named PLANNING BLOCKER [NEW planning item]

**The finding:** `chimeric_consensus.py` exists in the aligner but there is no chimera/trans-
splicing gate stratum in the benchmark. Per the standing discipline (no stratum = not yet a
proposal), any chimera-member proposal is blocked. This is not in the HANDOFF or planning docs as
an explicit item.

**Not an immediate priority:** this is a future/blocking dependency, not a top-3 build item. Note
it on the roadmap; do not build it until the top-3 are resolved.

---

## Q3. Updated build-and-gate-first sequence

**Sequence is the synthesis's top-3 with one insertion and two reuse/new-code clarifications:**

### 1. Panel-native measurement + double-prior ablation (co-first, shared corpus)

- **Panel job (32422876):** already in flight. On completion, read PANEL column per rung (AT-AC
  recovery). This is the pivot: < ~60% → member justified; > ~80% → pivot to correct-step + C6.
- **Double-prior ablation (NEW, add to this step):** write a new probe (M1, ~50 lines) running
  `junction_refiner` with prior ON vs OFF on the same ladder corpus. Separate from the aligner-
  generation measurement. M1-inline; does NOT need to wait for cluster job results.
- **Reuses:** existing NIC/NNC ladder corpus (`novel_junction_blindspot.py --emit-corpus` output).
- **New code:** one probe script (~50 lines, mirrors c3_junction_headroom.py pattern).

### 2. Motif-blind empirical-logP local realignment

- **Gated on step 1:** panel AT-AC recovery < ~60% AND true site addressable on hp_penalty -logP.
- **Reuses:** C1 engine (`local_aligner.py`), shipped penalty tables, NIC/NNC ladder corpus. Add
  cDNA-UMI-bin tables as per-read cost signal (the #2 new design spec above).
- **New code:** member wrapper + integration with junction_refiner (the motion the step-1 ablation
  will inform).
- Sim win is necessary-not-sufficient → transfer-check on real SIRV / SG-NEx after sim gate.

### 3. Variant-aware track (C6 + SMN1/SMN2 paralog)

- **C6 inline gate:** already in the HANDOFF inline list. Test winnowmap2 (already wrapped) +
  abPOA-pooling first as cheap baseline before proposing a local 2-copy graph.
- **SMN2 framing:** VCF-conditioned emission for c.840C>T + joint which-copy × which-junction
  readout. PI-redirect, Sumner-facing.
- **Reuses:** existing C6 stratum, `fp_variant_adjacent` gate, `junction_refiner` with VCF hook.
- **New code:** optional VCF emission context in walkback/refiner (the C6 extension).

**Chimera stratum:** goes on roadmap as a named blocker before any chimera member, not in top-3.

---

## Summary: what's genuinely new vs already-had

**Already had (synthesis confirms, doesn't add):** PANEL union column, motif-blind member concept,
~60% decision threshold, FDR guard.

**Top-3 genuinely new:**
1. `_CANONICAL_HP_PRIOR`/`is_novel` corrector ablation — distinct from aligner-generation
   measurement; same corpus, new code, M1-inline, potentially zero-member-win.
2. cDNA-UMI-stratified tables as per-read cost signal in the realigner design — new design spec,
   not standalone experiment.
3. Chimera stratum as a named planning blocker — new planning item, not immediate.

**Build sequence change:** the double-prior ablation becomes co-first with the in-flight panel
job (not after it), because it's M1-inline and tests a distinct failure mode the panel job cannot.
The rest of the sequence is unchanged.

*Reviewer: Sonnet 4.6 (agent a23805945182aa994). Code claims verified against repo (junction_refiner.py,
junction_scoring.py, panel_blindspot_score.py, penalty_tables/). Negative verified: planning docs
(dev/) contained NO prior proposal for the corrector ablation before the synthesis introduced it.*
