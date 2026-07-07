═══════════════════════════════════════════════════════════════════════════════
# ★ DURABLE HANDOFF (current state) — refreshed 2026-07-06
**Branch:** `worktree-agent-a25a2c1e784ad37dc` (agent worktree). Commit surgically (`git add <paths>`,
never `-A`); NEVER commit to `drs-validation-rebuild`; don't push without asking. The detailed
blow-by-blow log is BELOW this section; this is the clean summary.

## THE OUTCOME (what to ship)
The native junction re-aligner = **motif-blind refinement (`motif_blind=True`) + the targeted HP-drift
guard (`hp_drift_margin ≈ 3.0`)**. Motif-blind gives non-canonical discovery sensitivity; the guard adds
homopolymer-specific specificity with **zero discovery cost**. The empirical −logP `penalty_table` does NOT
belong in the re-placer (net-harmful in the boundary search); it keeps its validated home in **consensus
ranking + exon-interior indel scoring**. Reconciliation: consensus RANKS fixed motif-anchored alignments
(safe); re-placement SEARCHES boundaries (the intron-absorption degeneracy needs a prior, not the likelihood).

## DONE (committed, this arc)
- `motif_blind` toggle in junction_refiner (byte-identical off).
- Ground-truth yeast prp18Δ long-read sim + vetting harness (`scripts/benchmark/noncanon_sim/`), pbsim3 +
  independent-fallback error models.
- arm-A/B/C vetting → the −logP penalty table HURTS junction re-placement (over-shifts into HP; on realistic
  majority-undercall reads it's net-negative). Blunt hold-margin can't salvage it (no sweet spot).
- **HP-drift guard** (`_hp_run_across` + `hp_drift_margin`/`hp_drift_min_run`), threaded through all 4 refine
  fns, byte-identical at 0. `run_arms` arm-E = motif_blind + guard, no table.
- Tuned `hp_drift_margin ≈ 3.0` (drift-fix plateaus 3–4; discovery flat at all margins — decoupled).
- 11-test regression suite `tests/test_hp_drift_guard.py` (full refiner suite: 70 passed).
- Commits: 69a230f 07a712c 7ceed77 0984693 40b865a 796226b 92937d5 d70b036 (see `git log`).

## VERIFIED (numbers)
- Guard FIXES canonical HP-drift: D0 recovery 0.31→0.91 (margin 2) / 0.99 (margin 3); refines 215 vs arm-C 1085.
- Guard PRESERVES non-canonical discovery: realistic R3-HP arm-B 0.284 == arm-E 0.284 (touches 0 reads;
  acceptor is a transition, spared). R3-plain preserved.
- Byte-identical when guard off (59 refiner tests + 11 new guard tests green).

## OPEN / IN-FLIGHT
- **Real-DRS transfer test — sbatch `32967406` on Sherlock (PENDING in queue as of 2026-07-06).** Refines the
  real wt_by4742_rep1 DRS minimap2 BAM (309MB) twice (arm-B margin0 / arm-Bguard margin3, same pool) + compares
  at 119 HP-abutting annotated junctions (truth = SGD annotations). PI decisions: dataset = BY4742 DRS;
  discovery half = sim-proven + real do-no-harm.
- HELD (task 9): outer junction-enumeration loop + arm-C (−logP) — NOT pursued; the table doesn't earn its keep
  in the re-placer, so this is effectively closed unless the coherent `del_open_delta` law (untested) is tried.

## RESUME (concrete)
1. **Real-DRS job:** `ssh sherlock 'sacct -j 32967406 -X -o State%14; cat /scratch/users/kevinroy/real_drs_out/.real_drs_rc 2>/dev/null'`
   - sentinel `.real_drs_rc == 0` → `ssh sherlock 'python -c "import json;print(open(\"/scratch/users/kevinroy/real_drs_out/real_drs_hp_drift.json\").read())"'`; interpret vs EXPECT:
     Bguard `annotated_match_at_hp_abutting` > B (drift fixed on REAL undercalls); B/Bguard overall match not
     lower (do-no-harm); `reads_differing_between_arms ≈ reads_differing_at_hp_abutting` (guard is HP-specific).
   - rc != 0 OR sacct FAILED/TIMEOUT/OOM → `ssh sherlock 'tail -40 /scratch/users/kevinroy/real_drs_out/slurm-32967406.log'`, fix, `sbatch /scratch/users/kevinroy/real_drs_run.sbatch`.
   - still PENDING/RUNNING → poll again later (Sherlock queue wait).
   - Deployed guard code lives at `/scratch/users/kevinroy/rectify_guard` (my package overlaid; production
     `/oak/.../software/rectify` UNTOUCHED). To redeploy after a code change: rsync the changed .py there.
2. **Re-run any sim panel locally (M1):** `cd scripts/benchmark/noncanon_sim`; build_panel → gen_reads
   (`--force-fallback --hp-del-mult 8` for realistic undercall, or pbsim3 default) → run_arms
   (`--refine-workers 4 --hp-drift-margin 3.0`) → paired_arm_test (`--ref B --test E`). PY=/Users/kevinroy/miniconda3/bin/python.
   NOTE: M1 kills long/heavy refine jobs — keep reads ≤ ~5k or use `_make_arm_e.py` (single-arm on an existing bam).

## FILES
- Core: `rectify/core/splice/junction_refiner.py` (motif_blind, hold_margin, hp_drift_margin, `_hp_run_across`).
- Harness: `scripts/benchmark/noncanon_sim/{SPEC.md,build_panel.py,gen_reads.py,run_arms.py,paired_arm_test.py,
  real_drs_hp_drift.py,_make_arm_e.py,_sweep_refine.py}` (output dirs gitignored).
- Tests: `tests/test_hp_drift_guard.py`.
- This handoff (detailed log below).
═══════════════════════════════════════════════════════════════════════════════

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

## COMMITTED + v3 IN FLIGHT (2026-07-03)
Milestone committed: 69a230f "feat(benchmark): motif-blind arm-B re-aligner toggle + yeast non-canonical
long-read vetting sim" (junction_refiner motif_blind toggle + scripts/benchmark/noncanon_sim/ + handoff).
v3 Workflow `noncanon-sim-v3-pbsim3-fullpanel` (Opus/xhigh: 2 builders + integrator): install pbsim3 (realistic
ONT errors) + fuller panel (R0/R1/R3 × plain/HP × decoy 2-4 × exon/intron sizes) + realistic WT/cryptic mixture
(n_reads column, --cryptic-frac) + fix R2 tier + control TSV; then run full pipeline → mix_full_out/trade_curve.json
+ deliver the RIGOROUS arm-C-vs-arm-B verdict under realistic errors (does −logP beat motif-blind, esp. in HP cells).
On resume: check mix_full_out/trade_curve.json + the workflow verdict.

## ★★ v3 RIGOROUS RESULT (2026-07-03) — pbsim3 realistic errors, fuller panel; committed 07a712c
pbsim3 (ERRHMM-ONT, source-built, independent of penalty_scores) is the default read model. Fuller panel
R0/R1/R3 × plain/HP × decoy, WT+cryptic mixture (cryptic-frac 0.2), 4650 reads, N>=200/cell → mix_full_out/.
HEADLINE (robust): motif-blind arm-B recovers 0.788 of non-canonical junctions vs incumbent arm-A 0.568
(+22pp) at EQUAL false-junction FDR (~0.006) and 0.000 on INTRONFREE. arm-A flattens ONLY when a canonical
decoy is present (R3-plain-no-decoy: arm-A holds 0.835) — decoy-specific, correct.
arm-C vs arm-B VERDICT (make-or-break): TIE overall (C 0.781 vs B 0.788). arm-C's −logP shows a directional
edge ONLY in the HP-context non-canonical cell (R3-HP-decoy: rec 0.745 vs 0.710, FDR 0.015 vs 0.045) — where
its HP-length law applies — but WITHIN NOISE at n=200; slightly WORSE in R1 semi-canonical cells. So the
−logP increment is NOT demonstrated significant at this scale; motif-blindness (arm-B) is the robust win.
NEXT (to resolve arm-C): larger HP-focused panel / more reps to test if the HP-cell edge is real; realistic
cryptic-fraction sweep; decoy-offset/exon-size sweeps; then decide whether arm-C ships or arm-B alone.
On resume: mix_full_out/trade_curve.json holds it; re-run per the v3 pipeline (pbsim3 default).

## ★★★ arm-C VERDICT (2026-07-05) — the −logP table HURTS junction placement near HP runs
Advisor-designed hp_dist panel (canonical junction × flanking-A-run DISTANCE sweep, one-sided acceptor
primary + donor mirror + base-class + degenerate controls), pbsim3 reads, 700 reads/cell, paired McNemar +
bootstrap (paired_arm_test.py → mix_hpd_out/). Metric verified sound (normalize-equality, GT/AG-anchored;
mis-placements incl. 8bp NOT counted recovered). Fork-parallel refine (run_arms --refine-workers; forced
mp fork for macOS).
RESULT — OPPOSITE of the pre-registered hump: arm-C is ≤ arm-B EVERYWHERE (no niche where −logP wins).
 - ACC_A distance sweep: dC-B = -0.067(D0), -0.104(D1), -0.020(D2), then ≈-0.01/0 at D3-D10. Worst at
   SHORT distance (HP near the junction). BOT_A_D1 -0.050. poly-T ACC_T_D3 -0.081 (9× worse than poly-A
   at same distance — base-class calibration doesn't transfer A→T).
 - arm-C refined 442 reads vs arm-A/arm-B 164 — arm-C OVER-MOVES.
 - MECHANISM (BAM-verified, ACC_A_D1): arm-C shifts 73 reads 2bp INTO the A-run (180→182), absorbing A's
   into the intron; arm-B holds 606/612 at the true acceptor. The HP-cheap-del table makes the wrong
   HP-absorbing shift cheap. Exactly the PI's "moves an A into the intron", and arm-C is MORE prone to it.
VERDICT: arm-B (motif-blind, flat cost) is the win. arm-C's −logP penalty table does NOT earn its keep for
junction RE-PLACEMENT — it ties arm-B on non-canonical recovery (v3) and LOSES on canonical placement
precision near HP runs (this test). The −logP law's validated value is its ORIGINAL exon-block indel
attribution (the C1 false-indel test), NOT junction re-placement. SHIP arm-B; do NOT ship arm-C's table in
the junction scorer. Sub-finding: the AT/CG base-class table mis-transfers poly-A→poly-T (follow-up).
NEXT (per PI's approved sequence): human-transcript transfer check (human geometry + human penalty table;
then real gencode contexts).

## ★★★★ CORRECTED arm-C VERDICT (2026-07-05) — red team + PI overturned "arm-C hurts / don't ship"
The prior "arm-C hurts, never wins, don't ship" is WRONG on the "never wins" leg and OVERSTATED on harm.
Independent 3-lens Opus red team (design/impl/stats) + PI's logical-error catch:
1. STATS lens REFUTES "never wins": disaggregating mix_full non-canonical by rung — on R3 (deepest non-YAG,
   the prp18Δ biology) arm-C UNIQUELY rescues c=13 reads arm-B misses vs b=4 lost, ALL in the HP context
   (R3×HP c=12 b=4 +2pp; R3×plain null). arm-C is DOUBLE-EDGED, not purely harmful.
2. IMPL lens — real mechanism (corrects earlier framing): on the harmed canonical reads, arm-B scores
   truth==shifted EXACTLY 2.5000 (a TIE) and is_alt holds truth (right); arm-C's empirical fractional costs
   BREAK the tie toward the biologically-wrong intron-growth (truth 1.879 vs shifted 1.237). So harm =
   turning an ambiguous tie into a confident wrong answer. SAME machinery that helps R3. Metric verified
   correct (length-changing shift, not a normalize-slide). PI's "undecidable HP-abutting" intuition VALIDATED
   in spirit (arm-B scores exact ties there) — but arm-B handles the ambiguity right, arm-C doesn't.
3. DESIGN lens qualifications: (a) we tested the `penalty_score` reciprocal-rate COLUMN, NOT the coherent
   -logP law (`del_open_delta`=λ·ln(rate(hp)/rate(1)), no re-placer hook, UNTESTED); (b) arm-C's UPSIDE is
   under-tested — pbsim@0.95 reads have NO HP-length-dependent deletion structure (flat ~1.8%), so arm-C's
   raison d'être (absorbing a genuinely HP-collapsed junction) was never presented. Real DRS has ~4x more HP
   miscalls → harm scales UP, but upside also untested.
4. CHEAP FIX (impl lens): motif_blind has NO hold-margin/canonical guard — arm-C displaces a correct junction
   on a MARGINAL score win. A hold-margin (alt must beat incumbent by a threshold) would veto the harmful
   shifts while keeping the R3 rescues → potentially arm-C >= arm-B.
CORRECTED VERDICT: arm-C is DOUBLE-EDGED. Harm (canonical over-shift near HP) is real but a fixable
tie-breaking issue; benefit (rescue deep non-YAG near HP) is real but under-tested. NOT "don't ship" — it's
"needs the coherent del_open_delta law + a hold-margin + reads with real HP-length-collapse + ultimately real
DRS truth." arm-B remains the robust core (+22pp over incumbent); arm-C is a refinement to earn or reject on a
FAIR test. Redesign per PI: drop D0 abutting (degenerate), HP at distance + noisy flank, + arm-C's upside case.

## ★★★★★ PLAN 1+3 FINAL VERDICT (2026-07-05) — hold-margin can't salvage arm-C; ranking≠search reconciliation
Plan 1 (hold_margin, threaded through junction_refiner, byte-identical@0, 59 tests green; arm-D = arm-C+margin
in run_arms) + Plan 3 (realistic majority-undercall reads: fallback hp_del_mult=8 → 93% of reads undercall).
TRADE-OFF (paired McNemar, realistic reads):
- arm-C (no margin) under realistic reads = NET-NEGATIVE: large canonical harm (over-shift into undercalled
  runs, refines 1085 vs B 580), tiny R3 benefit (+0.012 ns; earlier clean-read +2pp evaporated).
- hold_margin=1.0: FIXES canonical harm (arm-D≈/>arm-B, D2 -0.023) but KILLS R3 benefit (R3-HP recD 0.329 vs
  0.665, -0.336).
- hold_margin=0.5: PRESERVES R3 (R3-HP -0.015) but does NOT fix canonical harm (D2 -0.115).
- NO SWEET SPOT: canonical harm needs margin≥~1.0, R3 needs margin≤~0.5 — requirements don't overlap. A global
  "don't move" prior can't separate harmful (canonical over-shift) from beneficial (R3 rescue) moves — same
  score magnitude.
SIDE-FINDING (positive, separate from discovery goal): a hold-margin DRAMATICALLY improves CANONICAL placement
near HP (D0: arm-D 0.92 vs arm-B 0.31, +0.6) — because it holds minimap2's motif-anchored placement instead of
drifting. Motif-blindness is a liability for canonical-near-HP; a prior fixes it.
RECONCILIATION (PI's Q — why is the SAME table fine for consensus?): consensus RANKS fixed, motif-anchored
per-aligner alignments (the degenerate shift-into-run isn't a candidate) — max-likelihood ranking, correct +
validated (C1 false-indel test). Re-placement SEARCHES over junction boundaries → the intron-absorption
degeneracy (HP-del-at-truth vs shift-into-run, equal likelihood) → pure likelihood picks wrong without a prior.
The table is a good LIKELIHOOD model; ranking it is right, searching boundaries with it (no prior) is not.
FINAL: arm-B (plain motif-blind) is the non-canonical-discovery tool. The penalty table stays in its validated
role (consensus ranking + exon-interior indel scoring), NOT junction re-placement search. Do NOT ship
penalty_table_path in the re-placer. Genuinely-open next steps: (2) coherent del_open_delta law (untested,
different quantity), (4) real DRS truth. Low-risk audit worth doing: does consensus ever pick a shift-into-HP
over a motif-anchored alt (mapPacBio the suspect)?

## ★★★★★★ HP-DRIFT GUARD (arm-E) WORKS (2026-07-06) — motif-blind discovery + HP-specific specificity, no trade-off
The PI-designed TARGETED guard succeeds where the penalty table (arm-C, hurt) and blunt hold-margin (arm-D,
killed discovery) failed. junction_refiner: `_hp_run_across` flags a move that lands a boundary INSIDE a
homopolymer run (≥ hp_drift_min_run); `hp_drift_margin` applies extra evidence margin ONLY to such into-HP
moves, sparing moves to genuine sequence transitions (real non-canonical acceptors). Byte-identical at 0
(59 tests green). arm-E = motif_blind + hp_drift_margin, NO penalty table.
RESULT (realistic majority-undercall reads, hp_drift_margin=2.0):
- Canonical HP-drift (mix_fair_hpg, arm-B vs arm-E): arm-E >= arm-B EVERYWHERE — D0 +0.603 (0.307->0.910),
  D2 +0.068, tied D1/D3/D5/D8. Fixes the false-non-canonical fabrication; refines 215 vs arm-C's 1085.
- Non-canonical discovery PRESERVED: R3-plain (non-canon AC at transition) +0.000; REALISTIC R3-HP (non-canon
  AC at transition + 8-A run DOWNSTREAM) arm-B 0.284 == arm-E 0.284 (+0.000, b=0 c=0) — guard touches 0 reads.
- Mechanism verified: _hp_run_across(true non-canon acceptor)=0 (transition, spared); =8 into the downstream
  run (drift vetoed). Earlier R3-HP -0.336 was a DEGENERATE test (acceptor placed inside an A-run, which no
  defined non-canonical dinucleotide is) — fixed build_panel R3-HP to a transition acceptor + downstream run.
VERDICT: SHIP arm-B + the HP-drift guard (hp_drift_margin) as the native re-placer: motif-blind sensitivity +
homopolymer-specific specificity, no discovery cost. The −logP penalty table stays in consensus/exon-indel
scoring (ranking), NOT the re-placer search. Params: hp_drift_margin (~2.0), hp_drift_min_run (4). Next: tune
hp_drift_margin, add a regression test for the guard, and the real-DRS transfer test (plan 4).

## hp_drift_margin TUNED (2026-07-06) — decoupled: drift-fix scales with margin, discovery flat at all margins
Sweep (_sweep_refine.py reuses aligned bam + pool; refine arm-E per margin, no re-align):
- CANONICAL drift-fix (D0 recovery): 0.307(m0.5) → 0.515(1.0) → 0.910(2.0) → 0.993(3.0) → 1.000(4.0). Plateaus ~3-4.
- R3-HP DISCOVERY (arm-E vs arm-B): 0.284 == 0.284 at EVERY margin (0.5/1.0/2.0), Δ+0.000 — guard never fires
  on the transition-site non-canonical acceptor, so discovery cost is ZERO regardless of margin.
=> RECOMMENDED hp_drift_margin ≈ 3.0 (D0 0.993, near-complete drift-fix, refines 171; conservative). Can go to
4.0 for a complete fix at no discovery cost. Units = edit-distance; an into-HP move must beat the incumbent by
>= this to be accepted (HP undercalls give only ~1-2-edit spurious advantages, so 3 cleanly vetoes them).
Loose ends: add a unit/regression test for _hp_run_across + the guard; real-DRS-with-truth transfer test.

## REAL-DRS TRANSFER TEST — SCOPED + READY (2026-07-06); execution = next milestone
PI decisions: dataset = BY4742 R10.4 DRS (penalty-table calibration data); discovery half = sim-proven + real
do-no-harm (no real non-canonical DRS truth exists; prp18Δ is short-read). Guard regression test DONE (11 tests,
committed 92937d5). Comparison tool written: scripts/benchmark/noncanon_sim/real_drs_hp_drift.py.
DATA (Sherlock, verified reachable): real minimap2 DRS BAM
  /scratch/users/kevinroy/rectify_wt_by4742_rep1_26167419_0/wt_by4742_rep1.minimap2.namesorted.bam
  yeast genome+GFF on Oak /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/
TEST (truth=annotations): flag annotated junctions that ABUT a homopolymer (a boundary can drift <=3bp into a
run>=4); refine the real BAM twice (arm-B motif_blind hp_drift_margin=0, arm-Bguard hp_drift_margin=3, same
pool); measure annotated-match rate at HP-abutting junctions per arm. EXPECT: Bguard match@HP > B (drift fixed
on real undercalls); overall match not lower (do-no-harm); arms differ ~only at HP-abutting junctions (guard
HP-specific). real_drs_hp_drift.py computes all of this -> real_drs_hp_drift.json.
RESUME (deploy guard code to Sherlock + run — Sherlock rectify install at /oak/.../software/rectify may lack the
guard; do NOT overwrite the production install or push without asking):
  1. Make a scratch rectify copy: cp -r /oak/.../software/rectify /scratch/users/kevinroy/rectify_guard/ (or reuse).
  2. rsync (M1->Sherlock, tiny) the guard file + tool:
       rectify/core/splice/junction_refiner.py -> /scratch/users/kevinroy/rectify_guard/rectify/core/splice/
       scripts/benchmark/noncanon_sim/real_drs_hp_drift.py -> /scratch/users/kevinroy/rectify_guard/
  3. Run (sbatch or interactive; PYTHONPATH=/scratch/users/kevinroy/rectify_guard):
       python real_drs_hp_drift.py --bam <the BAM> --genome <yeast.fsa(.gz)> --gff <yeast.gff.gz>
         --hp-drift-margin 3.0 --workers 8 --outdir /scratch/users/kevinroy/real_drs_out
  4. Read real_drs_out/real_drs_hp_drift.json; interpret vs the EXPECT above.
NOTE: real_drs_hp_drift.py syntax-checked on M1; not yet run on real data. Refine on a full DRS BAM is many
reads -> cluster only (n_workers>1; fork ok on Linux). If the BAM is huge, subset to a few chroms first.

## ★ REAL-DRS RUN SUBMITTED (2026-07-06) — sbatch 32967406 on Sherlock (job in flight)
Deployed my rectify package (with the guard) to /scratch/users/kevinroy/rectify_guard (production install
/oak/.../software/rectify UNTOUCHED; nothing pushed). Setup verified on real inputs: 17 chroms, 385 annotated
introns, 119 HP-abutting. Job refines wt_by4742_rep1 DRS minimap2 BAM (309MB) twice (arm-B margin0 / arm-Bguard
margin3, same pool) + compares -> /scratch/users/kevinroy/real_drs_out/real_drs_hp_drift.json (+ sentinel
.real_drs_rc). A local watcher polls it.
RESUME: ssh sherlock 'sacct -j 32967406 -X -o State%14; cat /scratch/users/kevinroy/real_drs_out/.real_drs_rc 2>/dev/null'
  - .real_drs_rc==0 -> read real_drs_out/real_drs_hp_drift.json; interpret vs EXPECT (Bguard match@HP-abutting > B
    = drift fixed on REAL undercalls; overall annotated match NOT lower = do-no-harm; reads_differing_between_arms
    ~= reads_differing_at_hp_abutting = guard is HP-specific).
  - rc!=0 or sacct FAILED/TIMEOUT/OOM -> read real_drs_out/slurm-32967406.log, fix, resubmit real_drs_run.sbatch.
  - PENDING/RUNNING -> wait.
