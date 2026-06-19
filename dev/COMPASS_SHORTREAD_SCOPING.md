# COMPASS short-read junction discovery — scoping synthesis + architecture call (2026-06-18)

From the scoping team (`wbb7hmt0h`, 4 grounded agents + synthesis). RE-LENSED per the PI's architecture
call: **back-propagate RECTIFY's junction refinements INTO COMPASS; keep RECTIFY end-correction-focused;
short-read junction discovery lives in COMPASS.** (The team's synthesis packaged it as a new RECTIFY module
— overridden here; its technical findings are direction-agnostic and stand.)

## What the grounding established (corrects the prompt's framing)
- **COMPASS = the PI's own repo** (github.com/k-roy/COMPASS, "Comparison Of Multiple alignment Programs for
  Alternative Splice Site discovery"). Algorithm = per-READ best-alignment ARBITRATION across SIX aligners
  (BBMap, STAR×2, HISAT2×2, Magic-BLAST, GSNAP): per read pick fewest-mismatch, tie-break
  ungapped>gapped>annotated>shorter-intron. (PI cites a 2023 NAR paper — the agent saw only the README;
  trust the PI for the methods/thresholds.) The "Pleiss/Chanfreau lineage" framing in my prompt was wrong.
- **RECTIFY's in-codebase "COMPASS 3-pass"** (`junction_validator.py`) is **cross-SAMPLE, and ORPHANED**
  (no CLI import — verified dead code).
- **The real cross-ALIGNER engine already in RECTIFY = `build_junction_pool`** (`junction_scoring.py:489`):
  admits a junction on anchored support OR ≥2 independent FAMILIES, with `ALIGNER_FAMILY` (:84-90)
  deduplicating wrappers. **One-line gap:** `ALIGNER_FAMILY` has only the 5 long-read tokens → short-read
  families (star/hisat2/gsnap) map to "?" and concordance silently no-ops.
- **Ambiguity-aware matching** (`normalize_junction`/`_canonical_within_window`, `chimeric_consensus.py:59-136`)
  is self-contained (coords+seq) → clean to back-propagate.

## Back-propagation payload (RECTIFY → COMPASS)
1. **Ambiguity-aware junction normalization** (left-slide + canonical-within-window) — clean, self-contained.
2. **Cross-FAMILY concordance gate with wrapper dedup** (`build_junction_pool` logic + `ALIGNER_FAMILY`) —
   a SECOND FP-control layer ON TOP of COMPASS's per-read arbitration (they are distinct: arbitration picks
   the best alignment per read; concordance requires independent families to agree on a junction).
3. **FDR calibration discipline** (below) — controls + a shuffled-junction null.
RECTIFY stays end-correction. Optional long-term: a shared `junction_core` mini-lib both import, to prevent
drift (decide consciously vs a one-time back-port).

## Recommended short-read panel (all availability VERIFIED on Sherlock, no installs except HISAT2 via Lmod)
**STAR (single-pass) + GSNAP + HISAT2 + BBMap + minimap2 `splice:sr`** = 5 distinct families.
- STAR single-pass ONLY (PI: 2-pass raises novel FDR); one vote, never sole arbiter; do NOT drop it.
- **minimap2 `splice:sr`** (preset added v2.29, verified in installed 2.30) DOES emit N-ops — overturns the
  earlier "-ax sr is unspliced" exclusion (that was plain `-ax sr`).
- GSNAP/HISAT2 indices do NOT exist on Sherlock (~1h `gmap_build`/`hisat2-build` each). GSNAP/HISAT2 SIMD
  can SIGILL on AVX-512-excluded nodes → pick sse42/nosimd variant + smoke-test.
- Skip STAR-noncanonical/HISAT2-noncanonical for THIS validation (the 111 are all GT-AG → noncanonical
  modes add FDR for no gain).

## ⚠ The decision-critical correction (owning my earlier overreach)
**The prior STAR-1-pass run's POSITIVE control failed: corrob_609 = 87/609 = 14.3%.** A detector that
validates KNOWN-corroborated junctions at 14.3% is **too insensitive to declare anything an artifact.**
Declaring the 111 "refuted" from it repeats the exact error that triggered this exercise. **The verdict is
THREE-WAY: corroborated / refuted / INCONCLUSIVE — 0/111 must NOT collapse to "artifact."** My earlier
"likely artifacts" was premature. (The STAR-independent coverage test still LEANS artifact — all 111
introns not spliced out — but has an intron-retention confound and cannot be conclusive alone.)

**Also (risk #7): family-concordance is NOT self-evidently FDR-controlling** — short-read aligners can SHARE
the failure mode (short anchor across a genomic GT-AG), the very thing that could have produced the 111. The
gate MUST be calibrated against the gmap_noncanonical negative control + a shuffled-junction null, not asserted.

## STEP 0 — the no-regret immediate move (NO compute, EXISTING BAMs)
Fix the denominator: per junction, compute Illumina flank coverage from the 3 existing STAR BAMs; define
"expressed" = both flanks mean cov ≥10 in ≥2 reps. Report **validated / EXPRESSED**, not validated/609. This
decides whether STAR *sensitivity* (not the 111) is the bottleneck — and whether the full multi-aligner
re-alignment is even needed. Cheap, decisive, uses data already on disk.

## Decision tree (after the full multi-aligner panel + calibration)
- Calibrate on controls FIRST. If positive(expressed) HIGH (>~70%) AND negative ~0 → panel trustworthy:
  then 111 still ~0 → CONFIRMED artifacts (finalize keep-GMAP-behind-fences); nontrivial fraction validate
  → REAL GMAP-unique novels (revise Deliverable B).
- If positive(expressed) stays LOW even with the full panel → INCONCLUSIVE; escalate to the P0 sim
  benchmark + the C6 variant-aware check on the clustered loci (chr5:171388xxx / 181237xxx). KEEP-GMAP
  stays untouched pending the benchmark.

## Open questions for the PI (from the synthesis)
1. Reuse the existing SG-NEx STAR 2.6.0c BAMs, or re-run STAR 2.7.11b single-pass for consistency?
2. Panel = the 5-family no-install set, or the full 6-aligner COMPASS (adds subread/Magic-BLAST, need installs)?
3. PacBio A549 (1 SMRTcell) as a THIRD orthogonal modality in the per-candidate verdict, or separate axis?
4. Per-read COMPASS arbitration (Layer B) in scope, or is cross-family concordance (Layer A) the deliverable?
   (The immediate experiment needs only Layer A.)
5. Sign-off on the "inconclusive" escape hatch so a 0/111 is not over-read as "artifact."
6. Architecture: confirm back-propagate-to-COMPASS (vs a RECTIFY module) and whether a shared junction lib is wanted.

---

# PI DECISIONS (2026-06-18) + STEP 0 RESULT

**PI answers to the open questions:** (1) run STEP 0 — done; (2) the COMPASS paper = **Roy et al. 2023
NAR, PMID 37956322** (note: the *paper* describes COMPASS as STAR+BBMap; the *repo* has since expanded to
6 aligners — BBMap, STAR×2, HISAT2×2, Magic-BLAST, GSNAP); (3) **wire up FULL per-read COMPASS arbitration
(Layer B)**, not just cross-family concordance; (4) **do it INSIDE the COMPASS repo** (github.com/k-roy/COMPASS,
cloned to ~/work/COMPASS) — confirms the back-propagate-to-COMPASS architecture.

**STEP 0 result (existing STAR BAMs, no compute) — STAR sensitivity IS the bottleneck, not expression:**
| set | total | validated/total | EXPRESSED | validated/EXPRESSED |
| --- | --- | --- | --- | --- |
| the_111 | 111 | 0 (0.0%) | 111/111 (ALL expressed) | 0 (0.0%) |
| corrob_609 | 609 | 87 (14.3%) | 567 | 82 (14.5%) |
The expressed-denominator fix barely moved the positive control (14.3%→14.5%) → STAR is genuinely
insensitive to these (mostly novel) junctions; a 14.5%-on-known-positives detector CANNOT adjudicate the
111. **Verdict stays INCONCLUSIVE; the multi-aligner COMPASS rework is REQUIRED** to build a sensitive,
calibrated detector before any keep-GMAP / artifact call. The 111 are all expressed (not a coverage issue),
which keeps them genuinely open rather than dismissible.

**COMPASS core algorithm (README + paper, definitive):** per-read best alignment = FEWEST MISMATCHES vs
reference; ties (same score, different junction) broken by ungapped>gapped → annotated>unannotated →
shorter>longer intron; **soft-clipping DISABLED** (force full-read mapping so mismatches can't hide);
cutadapt 3' trim; samfixcigar CIGAR reformatting; then ambiguous-junction adjustment to species splice
signals + branchpoint selection + unspliced-read splicing-efficiency.

**Deferred RECTIFY TODO (for the next RECTIFY agent — PI said leave it):** uncommitted RECTIFY working-tree
changes NOT pushed — the `docs/figures/aligner_panel/` schematics (aligner_schematics.py + fig1/2/3
png/svg, exclude __pycache__), `docs/figures/generate_splice_classification_v3.py` + regenerated
splice_classification.png/svg, `docs/ALIGNER_RECOMMENDATIONS.md` (+4/-2), `dev/TODO.md` (+79). The
`drop_chimeric_winners` filter + test were already committed (366c885). Scratch to leave: `dev/_bb_*`,
`_oc_*`, `_smoke_*`, `_spotcheck_*`, `README_preview.html`, `scripts/validation_data/upf1d_2026_05/stage/`.
