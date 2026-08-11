# CODE-VERIFICATION REVIEW — DIRECTOR_ALGO_EVAL_SYNTHESIS.md
*Reviewer: independent Sonnet (claude-sonnet-4-6). Date: 2026-07-01.*
*Branch: worktrees/agent-a25a2c1e784ad37dc. Read-only; not committed.*

---

## Verdict table

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1a | `_CANONICAL_HP_PRIOR = 0.5` at junction_scoring.py:293 | **CONFIRMED** | junction_scoring.py:293 exact match |
| 1b | Calibrated so "one HP deletion (del_cost ≈ 0.44–0.85 at HP=1) gives canonical junctions tie-breaking advantage" | **CONFIRMED** | docs/penalty_tables_quickref.md:61-65 |
| 1c | `canonical_discount = _CANONICAL_HP_PRIOR if tier < 4 else 0.0` at junction_refiner.py:647 | **CONFIRMED** | junction_refiner.py:647 exact match |
| 1d | Tier-gating: "disabled at deep deviation" (non-canonical candidate, tier≥4) | **CONFIRMED — one nuance** | junction_refiner.py:616,646-650; `tier` = CANDIDATE's tier; WHOLE block gated on `tier_beats_alt = current_tier >= 4` (current junction must itself be non-canonical to enter this path) |
| 1e | Tier semantics: what tier<4 covers | **CONFIRMED — nuance** | junction_scoring.py:349-374: tier 0 = GT/GC+YAG, tier 1 = GT/GC+RAG, tier 2 = GT/GC+NBG or non-canon-5'SS+YAG, tier 3 = semi-canon, tier 4 = non-canonical at both sites. Tier<4 is wider than just GT-AG (it includes GC-AG and similar) |
| 2 | `is_novel` tiebreaker at junction_refiner.py:638,656-663; "annotated preferred over novel"; priority 4 | **CONFIRMED** | 638: `is_novel = 0 if (chrom, js, je) in annotated_set else 1`; 660/663: position 4 (0-indexed: score_cmp, pri2, pri3, is_novel, abs_delta, …); is_novel=0 for annotated, lower=better |
| 2b | is_novel can "actually flip a novel-vs-annotated call" | **CONFIRMED** | Yes — when score_cmp and tier/is_alt and is_alt/tier are all equal between candidates, is_novel is the deciding field |
| 3 | cDNA-UMI penalty tables `penalty_scores_cdna_umi{1,2,3plus}.tsv` exist | **CONFIRMED** | rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/ contains all three; columns: op_type, base_class, hp_length, rate_mean, count_total, penalty_score, low_count |
| 4 | "code wraps 8 long-read aligners (minimap2, mapPacBio, gapmm2, deSALT, uLTRA, gmap, winnowmap2, minisplice_mm2)" | **CONFIRMED** | multi_aligner.py has run_* functions for all 8; align_command.py:419-420 dispatches winnowmap2 and minisplice_mm2 as opt-in extras; align_command.py:461ff dispatches uLTRA, deSALT, gmap |
| 4b | "winnowmap2 and minisplice are ALREADY present" | **CONFIRMED** | align_command.py:111 lists both as valid --aligners choices; 419-420 dispatches them; 241ff provides --winnowmap-repetitive-kmers and --minisplice-model flags |
| 4c | "the active set is config-dependent (gapmm2 dropped for human-DRS)" | **CONFIRMED by docs** | docs/ALIGNER_RECOMMENDATIONS.md; MultiAlignerConfig enables minimap2/mapPacBio/gapmm2 by default (multi_aligner.py:291-305) |
| 4d | ORIGINAL DISPUTE: "quick grep found only minimap2/mapPacBio/gapmm2 via AlignerConfig name= pattern" | **DISPUTE WAS TOO NARROW** | MultiAlignerConfig only lists 3 Tier-1 defaults; the full panel (all 8) is dispatched via align_command.py, not multi_aligner.run_multi_aligner |
| 5a | gapmm2 "IS minimap2 + edlib terminal rescue that forces GT-AG" | **CONFIRMED** | docs/aligners/gapmm2.md:137 "gapmm2 = minimap2 mapping + a single-threaded Python edlib loop"; ALIGNER_RECOMMENDATIONS.md:85 "forces GT-AG via edlib" |
| 5b | "≥85% of its unique introns are forced-canonical artifacts" | **CONFIRMED — sourced** | docs/aligners/gapmm2.md:194,207: "≥85 % of gapmm2's unique introns" and "~85 % of those uniques are novel/unannotated with an artifact signature (clustered forced GT-AG terminal rescues)". Applied to UNIQUE introns only; 1.3-1.8% of all gapmm2 introns are unique to gapmm2. |
| 5c | minisplice "scores GT/AG only" | **CONFIRMED** | docs/aligners/EVALUATED_AND_DISQUALIFIED.md:69 "A 1D-CNN that scores canonical GT/AG sites genome-wide"; ALIGNER_RECOMMENDATIONS.md "the prior is GT-AG-specific" |
| 5d | fp_variant_adjacent FDR=1.0 on ≥40bp GT..AG deletions | **CONFIRMED — measured result** | dev/C6_DESIGN.md:78 "HEADROOM=1.0000 (180/180)"; C6_DESIGN.md:85 "minimap2 fabricates a variant-adjacent FP intron on 100% of the ≥40 bp drivers"; also SIMULATION_BENCHMARK_SPEC.md:162; DRIVER_LENS=[40,60,100] in c6_headroom.py:87 |

---

## Three claims warranting precision fixes (non-blocking but affect planning)

### A. The tier-gating and candidate window scope — "discovery tax" framing is OVERSTATED

The synthesis says `_CANONICAL_HP_PRIOR` applies as a "double-prior penalty" that creates a "discovery tax" for non-canonical junctions. This framing overstates the scope. Two corrections:

**TWO gates, not one.** The mechanism has an outer gate AND an inner gate:

1. **Outer gate** (`tier_beats_alt = current_tier >= 4`, junction_refiner.py:616): the ENTIRE canonical-discount block only fires when the **current N-op** is itself non-canonical. Reads already placed at canonical junctions never enter this path.
2. **Inner gate** (`tier < 4`, line 647): within the outer gate, the discount applies only to **canonical CANDIDATE** alternatives. Non-canonical candidates get no extra penalty.

**±50bp repositioning window, not global discovery.** The candidates for scoring are drawn from `_candidates_near()` (junction_scoring.py:630), which applies `start_radius=max_boundary_shift` (default 50 bp) and `end_radius=max_boundary_shift` (default 50 bp) per endpoint. The wider `search_radius=5000 bp` (discovery radius) is irrelevant for the endpoint-shift constraint because `start_window = min(radius, start_radius) = min(5000, 50) = 50` (junction_scoring.py:676). The pool includes both annotated junctions and observed (novel) junctions from all aligner BAMs, but any candidate must have both endpoints within ±50 bp of the current N-op endpoints. Note: `max_slide=10` in the function signature is "retained for API compatibility; unused by the current HP-anchored semi-global implementation" (junction_refiner.py:490-491) — the actual constraint is ±50bp.

**Impact on the ablation experiment:** Toggling `_CANONICAL_HP_PRIOR` ON vs OFF measures how often the 0.5 discount tips junction SELECTION between a canonical and a non-canonical alternative within ±50bp of the current endpoints — for reads ALREADY placed at non-canonical junctions. This is a real and measurable effect, but it is NOT "non-canonical discovery across the transcriptome." The synthesis's #1-ranked ablation experiment is valid in design, but the scope should be stated as "repositioning of non-canonically-placed reads toward nearer canonical alternatives," not a global discovery penalty.

### B. 8-aligner claim requires context: "wired" means the CLI, not run_multi_aligner

The synthesis says "the code wraps 8 long-read aligners." All 8 have `run_*` functions and CLI access via `align_command.py`, but `run_multi_aligner()` in `multi_aligner.py` dispatches only 3 defaults (minimap2/mapPacBio/gapmm2) plus short-read alternatives. deSALT/uLTRA/gmap are labeled "Tier 2 opt-in" in the module docstring and dispatched via `align_command.py` separately from `run_multi_aligner`. Winnowmap2 and minisplice are further opt-in extras ("not in 'all'", per align_command.py:119). The synthesis correctly says "active set is config-dependent" and the claim is accurate, but the term "wraps" should be understood as "has a CLI-accessible wrapper function," not "is in the default multi-aligner pipeline."

### C. fp_variant_adjacent "FDR=1.0" is a measured simulation result, not a static code assertion

The "FDR=1.0" value comes from running `c6_headroom.py` with reps=60, seed=7 (documented in dev/C6_DESIGN.md). The static code only defines the drivers (DRIVER_LENS=[40,60,100]) and the measurement framework; the 1.0 value is a stored result. The synthesis correctly attributes it ("scored on the C6/VARIANT stratum"), and the result IS documented in the repo (dev/C6_DESIGN.md:78). No correction needed, but reviewers should understand this is "measured and stored" not "hardcoded."

---

## Bottom line

All five claim categories confirmed. One claim needed precision correction (A above): the synthesis's "double-prior discovery tax" framing for `_CANONICAL_HP_PRIOR` is **OVERSTATED** — the mechanism is a ±50bp repositioning knob for reads already placed at non-canonical junctions, not a genome-wide discovery penalty. The ablation experiment is valid in design but its scope should be stated precisely.

The "DISPUTED" aligner count was an incorrect dispute: winnowmap2 and minisplice are genuinely wired as CLI opt-in extras via align_command.py, not `MultiAlignerConfig` defaults. The synthesis is right to treat them as "measure the gap, not propose to add."

The FDR=1.0 (HEADROOM=1.0000) figure is a stored simulation result at dev/C6_DESIGN.md:78, not a static code assertion. Reviewers should know it requires re-running c6_headroom.py to reproduce.
