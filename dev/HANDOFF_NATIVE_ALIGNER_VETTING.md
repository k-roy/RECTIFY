# HANDOFF — native aligner vetting via a ground-truth prp18Δ sim (PI-directed pivot)

**Branch:** `worktree-agent-a25a2c1e784ad37dc` (the agent worktree). **Written:** 2026-07-03.
Supersedes the A549 real-data snap-or-hold track. Report to the PI (Kevin); commit surgically; NEVER
commit to `drs-validation-rebuild`.

## What changed this session (the load-bearing corrections)

1. **C1 (inner −logP scorer) is ALREADY BUILT + validated** on this branch (`dev/C1_DESIGN.md`,
   `align_exon_block_global(penalty_table, lam, ins_lengthlaw)`, `test_c1_lengthlaw.py`). The handoff that
   said "build C1" was stale; a wrong-design reimplementation was built on the WRONG branch
   (drs-validation-rebuild, which has C1 removed), caught before commit, discarded. See
   `dev/C1_BASELINE_CORRECTION_FINDINGS.md`.
2. **The 3 A549 "non-canonical anchors" don't reproduce as a do-no-harm substrate** — both short+long reads
   place them at canonical, ANNOTATED GT-AG introns; the harness coords are 175–1017 bp off (a cross-pipeline
   placement disagreement, unresolvable from surviving data). See `dev/SNAP_OR_HOLD_ANCHOR_COORD_FINDING.md`.
3. **PI decision (confirmed):** BUILD the native junction re-aligner. Vetting moves from A549 real data to a
   **ground-truth simulation**, because non-canonical splicing is a MUTANT phenomenon, not WT — A549 (WT-like)
   may lack it. Biological grounding = **Roy et al. 2023 NAR** (Prp18p promotes 3′SS fidelity; prp18Δ
   activates **non-YAG 3′SS**, enhanced by **upstream poly(U) tracts** + first-intronic-nt complementarity).

## The plan (advisor-vetted; two-advisor consult done)

Build a **yeast prp18Δ-grounded ground-truth panel** + vet arm-A/B/C:
- **Panel:** non-YAG 3′SS motif gradient (YAG canonical control → RAG → non-YAG "unusual") × upstream
  poly(U)-tract length × flanking exon sizes (yeast-scale) × first-intronic-nt complementarity. Each non-YAG
  acceptor has a canonical YAG a few bp away = the SNAP target. By-construction junction truth.
- **Arms:** A = motif-biased incumbent; **B = motif-blind** (canonical/annotation priors neutralized);
  C = B + empirical penalty table (the C1 −logP law).
- **PRIMARY METRIC = the recovery-vs-false-junction-FDR TRADE curve**, stratified by motif-deviation rung ×
  poly(U) length. Do NOT lead with the recovery number — arm-B recovers the injected non-YAG by construction
  (tautological). The aligner earns its keep only if recovery doesn't buy proportional fabrication at noisy
  canonical sites (INTRONFREE + YAG-canonical controls).

### Two hard constraints (advisor)
- **Error-model independence is STRUCTURAL (verified):** `error_injector.py` uses placeholder marginal rates
  (NOT `penalty_scores.tsv`); `scorer.py` never reads the penalty table; pbsim3/controlled carry external
  truth. So arm-C's win cannot be circular by construction — but do NOT recalibrate `penalty_scores` from sim
  output.
- **HP-context is load-bearing for arm-A/B, not just arm-C:** the prp18Δ sites are poly(U)-activated, and the
  elevated deletion IN that poly(U) tract is what makes the non-YAG-vs-YAG placement ambiguous (the snap-vs-
  hold decision). `error_injector.py` is HP-FLAT → arm-B would hold trivially → won't transfer to real
  prp18Δ. **An independent HP-context-aware injector is a NEAR-TERM dependency.** Build plumbing on the flat
  injector first; any arm-A/B conclusion is PROVISIONAL until the poly(U) carries realistic, independent-
  source deletion structure.

## DONE this session
- **Arm-B motif-blind toggle wired into `junction_refiner.py`** (`motif_blind: bool = False`), threaded through
  `refine_read_junctions` (neutralizes `tier_beats_alt`, the `_CANONICAL_HP_PRIOR` discount, and drops
  `tier`/`is_novel` from the tie-break sort), `_run_sequential`, `_run_parallel` (via the worker kwargs dict),
  and `refine_bam_junctions`. **Byte-identical when off** — all 59 `test_junction_refiner`/`false_junction`/
  `anchor_gate`/`simple_slide` tests pass. **UNCOMMITTED WIP** (hold commit until the panel exercises the ON
  path). Files changed: `rectify/core/splice/junction_refiner.py` only.

## Code surface (verified this session)
- Arm-B neutralization site: `junction_refiner.py` ~614 (`tier_beats_alt`), ~646 (`canonical_discount`),
  ~658–665 (sort tuples; now use `tier_key`/`novel_key` = 0 when motif_blind).
- `_CANONICAL_HP_PRIOR=0.5` def `junction_scoring.py:293`, consumed `junction_refiner.py:647`.
- Existing sim infra (reuse): `scripts/benchmark/human_noncanon_sim.py` (by-construction junction truth, OBS/
  LADDER/INTRONFREE arms — human-scale, hardcoded panel → generalize to yeast + `--motif-list` + poly(U)),
  `scripts/benchmark/sim/{error_injector,controlled,transcript_model,pbsim3_wrapper}.py`,
  `rectify/core/benchmark/scorer.py` (ambiguity-aware position-exact concordance, no penalty-table use),
  `scripts/benchmark/error_realism_validate.py`, `dev/SIMULATION_BENCHMARK_SPEC.md`.

## RESUME — next steps (in order)
1. **Yeast prp18Δ panel generator** (task 10): generalize `human_noncanon_sim.py` to yeast + parameterize the
   motif panel (`--motif-list`: non-YAG gradient) + add the poly(U)-tract context axis + yeast-scale
   exon/intron sizes. By-construction truth (donor GT fixed; acceptor = non-YAG with a canonical YAG neighbor).
2. **End-to-end ~5 constructed cases**: generate → align (minimap2) → `refine_bam_junctions(motif_blind=False)`
   (arm-A) vs `(motif_blind=True)` (arm-B) → `scorer.py`. Confirm plumbing + that arm-A snaps non-YAG→YAG while
   arm-B holds. This is also the arm-B ON-path behavior test.
3. **Independent HP-context-aware injector** (task 11): per-(HP-length,base) deletion structure from a source
   independent of `penalty_scores.tsv` (SIRV-fit / pbsim3 / held-out raw-read HP profile). Re-run the panel;
   only THEN report arm-A/B conclusions as non-provisional.
4. **Arm-C** = arm-B + penalty table; then the recovery-vs-FDR trade curve across all arms.
5. Commit arm-B once (2) exercises it; keep `test_junction_refiner` green.

## Repro / env
- Tests (agent worktree): `cd <agent-worktree> && /Users/kevinroy/miniconda3/bin/python -m pytest tests/test_junction_refiner.py -q`.
- Yeast genome+annotation bundled: `rectify/data/genomes/saccharomyces_cerevisiae/`.
- Sherlock (for cluster runs): `source …/anaconda3/etc/profile.d/conda.sh; conda activate rectify; PYTHONPATH=…/nj_panel_code`.
- Prior anchor-investigation scratch (for-record): `/scratch/users/kevinroy/c1_realdata_dB/`.

## IN FLIGHT (2026-07-03)
Workflow `build-yeast-noncanon-sim` (5 Opus/xhigh agents) building the panel+reads+arms+score pipeline under
`scripts/benchmark/noncanon_sim/` against `scripts/benchmark/noncanon_sim/SPEC.md`, then an integration smoke
run → `scripts/benchmark/noncanon_sim/smoke_out/trade_curve.json`. PI-corrected biology: poly(U) is intronic
(spliced out, NOT in the read); model non-YAG ACCEPTOR + canonical-YAG decoy + exonic-context ambiguity.
Goal = simulated non-canonical LONG reads + the recovery-vs-false-junction-FDR trade curve for arm-A/B/C.
On resume: check the workflow result / `smoke_out/`; if trade_curve.json exists, review recovery-vs-FDR;
else re-run the pipeline scripts (build_panel → gen_reads → run_arms → score_trade) per SPEC.

## v1 RESULT + v2 IN FLIGHT (2026-07-03, later)
v1 pipeline PASS end-to-end (smoke_out/trade_curve.json). Independent error model = fallback (pbsim3 absent),
deletion-dominant/HP-elevated, verified NOT from penalty_scores.tsv. KEY FINDING: arm-A's snap is real but
GATED by junction-pool admission — build_junction_pool needs clean-anchor support; smoke panel simulated ONLY
non-canonical reads so the canonical decoy was never pooled → arm-A==arm-B. Forced-in probe: arm-A snaps 53/53,
arm-B 0/53. FIX (v2, agent ab931dea in flight): model prp18Δ WT+cryptic MIXTURE — abundant canonical(WT) reads
+ minor non-canonical(cryptic) reads per locus → canonical decoy becomes an observed pool member → arm-A
flattens cryptic onto WT, arm-B holds. Output → scripts/benchmark/noncanon_sim/mix_out/trade_curve.json.
Success = recovery(cryptic) arm-A < arm-B, controls FDR ~0. Also fixing: R2 tier unrealizable, control-row TSV
alignment, R3-HP fjFDR artifact; trying to install pbsim3.
On resume: check mix_out/trade_curve.json for the A<B recovery gap on cryptic cells.

## ★ v2 RESULT — the demonstration LANDED (2026-07-03)
mix_out/trade_curve.json (WT+cryptic mixture; the agent stalled → I implemented the co-located WT isoform in
build_panel.py myself). Ground-truth, BAM-verified:
- arm-A (incumbent, motif-biased): FLATTENS 100% of cryptic non-canonical reads onto the WT canonical decoy
  (chrSIM_2: 98/98 cryptic reads moved 180→183). recovery(R3)=0.000.
- arm-B (motif-blind): HOLDS the true non-canonical — 75/98 at 180 (per-cell recovery 0.61–0.63; hold 76%).
- arm-C (−logP): slightly better — 81/98 at 180 (recovery 0.63–0.68; hold 83%).
- false-junction FDR = 0.000 for ALL arms on ALL controls (INTRONFREE, R0+decoy, WT, YAG-no-decoy); WT
  canonical reads untouched by every arm (all stay at 183).
So: motif-blind/−logP re-aligner recovers ~62–83% of non-canonical junctions the incumbent flattens 100%,
at ZERO false-junction cost on clean sites. First honest ground-truth evidence the native aligner works.
CAVEATS: smoke panel (small, N~100/cell); independent FALLBACK error model (pbsim3 NOT installed, HP-elevated
but placeholder magnitude); arm-C increment over arm-B modest (provisional until independent HP-aware injector,
task 11); ~15% of arm-B/C R3 reads land 1bp off (within-cell cost to characterize).
FILES CHANGED (uncommitted): rectify/core/splice/junction_refiner.py (motif_blind toggle),
scripts/benchmark/noncanon_sim/{SPEC.md,build_panel.py,gen_reads.py,run_arms.py,score_trade.py}, mix_out/.
NEXT: install pbsim3 for realistic ONT errors; fuller panel + realistic mixture fractions; independent
HP-aware injector (task 11); then commit the milestone.
On resume: mix_out/trade_curve.json holds the result; re-run = build_panel --out-dir mix_out; gen_reads
--out-dir mix_out; run_arms --reads mix_out/reads.fastq --sim-ref mix_out/sim_ref.fa --outdir mix_out;
score_trade --work-dir mix_out --arm A=... B=... C=...
