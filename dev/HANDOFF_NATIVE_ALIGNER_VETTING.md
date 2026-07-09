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

## ★★★★★★★ REAL-DRS TRANSFER TEST — PASS (2026-07-06, sbatch 32967406, rc=0)
Real BY4742 wt DRS minimap2 BAM (309MB), refined twice (arm-B margin0 / arm-Bguard margin3), truth=SGD
annotations. Result (real_drs_out/real_drs_hp_drift.json + a corrected fix/harm pass):
- DO-NO-HARM: overall annotated-match 0.9882 (B) -> 0.9884 (Bguard) over 94,907 junctions — NOT lowered
  (slightly higher). Guard changed only 20/94907 placements (0.02%).
- HP-SPECIFIC: 16/20 differing reads are at HP-abutting junctions.
- DRIFT-FIX TRANSFERS: of the changed reads, 17 FIX (junction moved back TO the annotation, e.g. 223786->223781),
  0 HARM, 3 neutral. Every guard change on real data fixed a genuine drift; zero harm.
- Magnitude is SMALL on WT (only ~20 reads drift — WT reads at annotated junctions are ~98.8% correctly placed;
  the into-HP drift is rare in WT). But when the motif-blind refiner DOES drift on real undercalls, the guard
  catches it 100% correctly (17/17). Larger benefit expected in heavy-undercall / prp18d-type contexts.
METRIC BUG (cosmetic, in the committed real_drs_hp_drift.py): `annotated_match_at_hp_abutting` divides by ALL
junctions (n_at_hp_junctions=94907), so it reports the FRACTION of reads at HP sites (~0.30), not the match
RATE there. The decisive metric is the fix/harm classification above (17/0). Fix the tool before reuse.
=> REAL-DRS GATE PASSED. Whole arc lands: SHIP motif-blind refinement + HP-drift guard (hp_drift_margin ~3.0).

## IN-FLIGHT (2026-07-06 session cont.) — human transfer + del_open_delta triple-check
STATE: del_open_delta REJECTED (committed b4cab89; finding dev/DEL_OPEN_DELTA_FINDING.md + wiring preserved
dev/arm_f_del_open_delta.patch, reverted from production). Scope of the guard distilled into the layman doc
(cfa08b1 + the "narrow footprint" note). Metric bug in real_drs_hp_drift.py FIXED (fix/harm classifier now
the headline).

IN FLIGHT #1 — del_open_delta TRIPLE ADVERSARIAL PASS (Workflow whf5mx4qj, 3 Opus xHigh arbiters + synth).
  RESUME: watch /workflows or await the completion notification. If synth.all_clear -> del_open_delta is
  DROPPED DEFINITIVELY (close task #9). If a needs_rerun is flagged (most likely the untested arm-F+guard
  combo) -> apply dev/arm_f_del_open_delta.patch, run scripts/benchmark/noncanon_sim/_make_arm_f.py WITH
  hp_drift_margin=3.0, compare to arm-E, then revert the patch.

IN FLIGHT #2 — HUMAN chr5 TRANSFER TEST (sbatch 32972957 on Sherlock larsms).
  RESUME: ssh sherlock 'sacct -j 32972957 -X -o State%14; cat /scratch/users/kevinroy/human_drs_out/.human_rc 2>/dev/null'
   - .human_rc==0 -> ssh sherlock 'cat /scratch/users/kevinroy/human_drs_out/human_rep1_chr5.json'; interpret:
     guard_changes fix >> harm (drift-fix transfers to human); overall annotated_match_rate not lower
     (do-no-harm); match_rate_at_hp_abutting Bguard >= B; reads_differing ~= at_hp_abutting (HP-specific).
   - rc!=0 or FAILED/TIMEOUT -> ssh sherlock 'tail -40 /scratch/users/kevinroy/human_drs_out/slurm-32972957.log'; fix; sbatch /scratch/users/kevinroy/human_drs_run.sbatch
   - PENDING/RUNNING -> poll again (larsms queue).
  SETUP (done): BAM reheadered 5->chr5 at /scratch/users/kevinroy/human_drs_out/rep1.chr5.rehdr.bam; genome
   /scratch/users/kevinroy/compass_a549/COMPASS/genome_references/GRCh38_gencode_v44_chr5.fasta; GTF
   /scratch/users/kevinroy/sumner_lab/references/gencode.v44.basic.chr5.gtf; 13548 chr5 introns, 3007 HP-abutting.
  GOTCHA (roadmap): standardize_chrom_name is YEAST-CENTRIC (chr5->chrV; chr10->chrX collision). register_
   genome_contigs_from_fasta fixes it for a chr5-ONLY run (all reconcile to chr5); a FULL-GENOME human run
   needs standardize_chrom_name generalized first.

## UPDATE (2026-07-06 cont.2) — del_open_delta triple pass = HOLD; faithful re-run + human both in flight
del_open_delta TRIPLE ADVERSARIAL PASS done (Workflow whf5mx4qj): NOT all-clear.
 - Arbiter2 (mechanism): NO fault; CLOSED the over-shift axis analytically (P_shift into-run = emission-identical
   to truth except free-intron-N vs penalized-deletion-D label → P_shift=0 global min, del-cost-INVARIANT; λ=0≡arm-B
   caps axis [580 over-shifts,654] strictly above guard's 171). Over-shift axis SETTLED.
 - Arbiter1 (wiring): REAL staircase bug in _score_hp_affine_del — del→match(0)→del re-collects the one-time OPEN
   discount per in-run base, so arm-F degenerated to per-base (arm-C-like) IN RUNS. RECOVERY axis (D0/D2) never
   faithfully tested. Verdict survives (580 core is cost-free) but NOT droppable-as-definitive yet.
 - Arbiter3 (results): STALLED (API), no return.

IN FLIGHT #A — FAITHFUL arm-F RE-RUN (background Agent ada70ae3ca5401b1e). Fix staircase→true one-time-per-run
 discount, re-derive λ≈0.05-0.1 (graded, no floor), test arm-Ff+guard vs arm-E on mix_fair_out (D0/D2 recovery +
 over-shift) & mix_r3b_out (R3). DECISION: arm-Ff+guard fails to Pareto-beat arm-E → DROP del_open_delta
 DEFINITIVELY (close #9); beats it → coherent law belongs. RESUME: await agent completion notification; read its
 report (staircase-fix proof k=2,3, λ curve, arm-B/E/Ff/Ff+guard numbers, Pareto verdict). Both arbiters expect DROP.
IN FLIGHT #B — HUMAN chr5 (sbatch 32972957, RUNNING; bounded poll bu0u3tfob). RESUME as in the prior block
 (.human_rc sentinel → human_rep1_chr5.json; interpret guard_changes fix>>harm + do-no-harm).

## ★ HUMAN chr5 TRANSFER TEST — PASS (2026-07-07, sbatch 32972957, rc=0)
A549 chr5 direct-RNA, GENCODE-basic truth (13547 introns, 3007 HP-abutting), 69610 junctions scored/arm.
- DO-NO-HARM: overall annotated match 0.7777 (B) -> 0.7914 (Bguard) — NOT lowered (RAISED +0.0137 ≈ +956
  junctions moved back onto annotations). At HP-abutting: 0.810 -> 0.818. 0 harm.
- DRIFT-FIX TRANSFERS TO HUMAN: the guard held ~956 junctions that arm-B drifted OFF annotations into runs ->
  net +956 annotation matches, 0 harm. 950 reads changed (25x the yeast substrate); 777/950 (82%) at HP-abutting.
- WRINKLE (honest): the per-READ fix/harm classifier shows only 6 fix / 0 harm / 944 NEUTRAL — because (a) human
  reads are MULTI-JUNCTION (a read matches via any of its many junctions, so a single fixed junction reads as
  "neutral"), and (b) ~944 of the guard's changes are at NON-ANNOTATED (novel / cancer-specific) junctions that
  GENCODE-basic can't judge. The junction-level match rate (+956, all positive, 0 harm) is the reliable metric
  on human; the per-read classifier should be made per-JUNCTION for multi-junction organisms (follow-up).
  The guard's activity on NOVEL junctions (the discovery territory) needs the short-read (COMPASS) cross-check.
=> HUMAN do-no-harm + drift-fix CONFIRMED. Full novel-junction picture -> COMPASS short-read validation (next).
PERF (human-readiness): refiner ~0.24s/read on human (build_junction_pool 1.7s is fine); chr5 ~10min/arm at
8 workers. Full-genome human needs refiner perf work + standardize_chrom_name generalized (chr5->chrV; chr10->chrX collision).

## RESOLVED (2026-07-07) — del_open_delta DROPPED DEFINITIVELY (PI accepted)
Triple-pass caught a real staircase bug; faithful re-run fixed+re-tested; REJECT confirmed 3 ways (faithful
re-run + independent re-scoring arm_Ffg≡arm_E + Arbiter2 del-cost-invariance theorem). PI accepted the drop
without a further triple pass (corroboration = impossibility theorem + faithful empirical + independent verify).
The guard ALONE is the native re-aligner's answer; the empirical table stays OUT of the boundary search in
every form. Finding: dev/DEL_OPEN_DELTA_FINDING.md. Faithful wiring: dev/arm_f_del_open_delta_faithful.patch.
Task #9 CLOSED. Human chr5 transfer (#18) PASSED. Next milestone: COMPASS short-read cross-check for the
novel-junction activity annotation-truth can't score; then per-junction fix/harm metric; refiner perf +
standardize_chrom_name generalization for full-genome human.

## SHORT-READ VALIDATION SCOPED + Q1 STARTED (2026-07-07 cont.)
Scope VETTED (dev/SHORTREAD_DISCOVERY_VALIDATION_SCOPE.md, 3-critic pass): Q1 (guard/HP-drift on novels) VALID
after 5 must-fixes; Q2 (non-canonical discovery) RIGGED as drafted (short-read aligners share minimap2 motif-
snapping) — GATED on reviving COMPASS-A549 (NOT STAR-2pass). Q1-first.
IN FLIGHT:
 #1 Q1 POSITIVE CONTROL (sbatch 33013531, larsms) — base-exact Illumina split-read concordance at the 3007
    annotated HP-abutting chr5 junctions (+ all-annotated canonical floor). Script rectify_guard/
    q1_illumina_hp_concordance.py. RESUME: ssh sherlock 'sacct -j 33013531 -X -o State; cat /scratch/users/
    kevinroy/human_drs_out/.q1pc_rc 2>/dev/null; tail -30 /scratch/users/kevinroy/human_drs_out/q1_pc-33013531.log'.
    GO/NO-GO: HP-abutting base-exact concordance >0.3 (K>=1, >=2/3 reps) -> Q1 answerable, report vs this floor;
    ~0 while canonical floor high -> Illumina can't resolve HP -> Q1 needs COMPASS. Illumina BAMs verified:
    SN:5 (Ensembl), coord-sorted, indexed, ample split reads. rep1/3/5.
 #2 WEB SURVEY (background Agent abafcd2c8d191150e) — exhaustive prior non-canonical junction benchmarks +
    methods papers + spike-in prior art. RESUME: await completion notification; fold findings into the
    spike-in design + the benchmark plan.
 #3 SPIKE-IN TRACK (task #20, scoped not built) — spike noncanon_sim's non-canonical reads (as SYNTHETIC
    CONTIGS added to GRCh38, SIRV-model) into real A549 DRS -> ground truth by construction, sidesteps the Q2
    circularity. Measure recall/precision on the spiked junctions. Design after survey + Q1 PC land.

## Q1 POSITIVE CONTROL = GO; ADJUDICATION IN FLIGHT (2026-07-07)
Q1 PC (sbatch 33013531, DONE): base-exact Illumina split-read concordance HP-abutting=0.517 ≈ canonical=0.512
→ Illumina IS competent at HP base-exact placement (vet's sharpest risk REFUTED); ~52% decidable floor. Q1
answerable. (My proper split-counter beats the old 14.5% STAR-1pass figure.)
IN FLIGHT:
 #1 Q1 ADJUDICATION (sbatch 33013750). RESUME: ssh sherlock 'cat /scratch/users/kevinroy/human_drs_out/.q1adj_rc;
    tail -40 /scratch/users/kevinroy/human_drs_out/q1_adj-33013750.log'. Reports (per distinct guard MOVE at
    non-annotated junctions): annotated_in_comprehensive (adjudicated-free), truly_novel_residual, resolvable
    (crosses ambiguity window), FIX_guard_supported vs HARM_armB_supported (base-exact Illumina support-ratio,
    K=2 ratio=3, over decidable), inconclusive. EXPECT FIX>>HARM (do-no-harm into novels); HARM = investigation
    set (one-directional). Script rectify_guard/q1_adjudicate.py. If OOM/slow on the 1.5GB comprehensive GTF,
    subset it to chr5 first.
 #2 WEB SURVEY (Agent abafcd2c8d191150e) — non-canonical junction benchmarks + spike-in prior art; await notif.
 #3 SPIKE-IN track (#20) — scoped; design after survey.

## Q1 RESULT (re-partition) + SURVEY DONE (2026-07-07)
Q1 ADJUDICATION (sbatch 33013750): the vet's comprehensive re-partition transformed it — 944 per-read changes
= 173 DISTINCT junction moves; 168/173 (97%) annotated in COMPREHENSIVE gencode v44 (only "novel" vs basic);
153 clean FIXES (guard->comp-annotated arm-B wasn't at); 5 truly-novel (all zero Illumina coverage). Guard's
changes are overwhelmingly fixes onto REAL junctions. Harm-breakdown re-run (sbatch 33013898) IN FLIGHT to
count comp_HARM (guard-off-annotated) — expect ~0. RESUME: ssh sherlock 'tail -22 /scratch/users/kevinroy/
human_drs_out/q1_adj-33013898.log'. Scripts committed (q1_illumina_hp_concordance.py, q1_adjudicate.py).
SURVEY DONE (dev/NONCANON_BENCHMARK_SURVEY.md): PI's recalled study = LRGASP. Spike-in (#20) has a published
precedent (SQANTI-SIM: novel junctions by transcript-removal + NanoSim/pbsim3); the SPECIFIC combo (unannotated
non-canonical in real ONT-dRNA background) is GENUINELY NOVEL design space = the defensible contribution. Borrow
SIRV-in-real-dRNA (SG-NEx RNA002, LongBench RNA004=target chem). Adopt SQANTI3 NIC/NNC motif-stratified metrics.
NEXT: finalize Q1 harm count; then lock spike-in design (LongBench RNA004 bg + SQANTI-SIM-style injection).
Q2 still gated on COMPASS revival.

## SPIKE-IN DESIGN + SUMNER TEST (2026-07-07)
Q2 path = "both, spike-in first" (PI). SPIKE-IN design LOCKED (dev/SPIKEIN_DESIGN.md): SQANTI-SIM-style
unannotated non-canonical synthetic contigs (noncanon_sim) spiked into a real ONT-dRNA background; panel align
-> refine arm-A/B/Bguard; SQANTI3 NIC/NNC motif-stratified recall(flattening)+precision(fabrication). PENDING
DECISION: background dataset (LongBench RNA004 target-chem download [my lean] vs reuse SG-NEx/BY4742).
SUMNER SMN test (#21, SEPARATE JHU collab, dual-purpose per PI): data /scratch/users/kevinroy/sumner_lab/
chr5_bams/ (4 CNTL + 7 SMA, chr5-named, no reheader). BIOLOGY FIRST LOOK (raw junctions, chr5:69.9-71Mb): SMN1
(~70.94Mb) junctions strongly DOWN in SMA, SMN2 (~70.05-70.08Mb) UP — the SMN1-loss/SMN2-compensation mechanism,
recovered. IN FLIGHT: RECTIFY refine (sbatch 33015742) on SMA_GSB2945 + CNTL_HB53 (arm-B/Bguard do-no-harm +
arm bams for non-canonical). RESUME: ssh sherlock 'cat /scratch/users/kevinroy/sumner_out/.sumner_rc; for S in
SMA_GSB2945 CNTL_HB53; do cat /scratch/users/kevinroy/sumner_out/$S/${S}_drift.json; done'. NEXT: SMN2 exon 7
skip ratio SMA-vs-CNTL + non-canonical junctions the refine reveals. CAVEATS: SMN1/SMN2 paralog confound
(within-locus assignment provisional; locus-level robust); real data = discovery not recall/precision.
Q1 DONE (153 FIX/0 HARM). Q2 substrates: spike-in (ground truth) + Sumner (real biology) + COMPASS-A549 (later).

## GENOME-WIDE SUMNER — blocker resolved; awaiting ssh (2026-07-07)
PI: fan out Sumner to whole genome (chr5 had few novel events). Genome-wide BAMs ARE on cluster:
/scratch/users/kevinroy/sumner_lab/full_genome_bams/ (full SMA + WT panel, some reps; WT=control naming).
★ standardize_chrom_name "blocker" is a NON-ISSUE (verified locally): unregistered it collides chr10->chrX==chrX,
but register_genome_contigs_from_fasta(FULL GRCh38) returns ALL human chroms VERBATIM (chr10->chr10, no
collision). NO CODE FIX — just pass the FULL genome (GRCh38.primary_assembly.genome.fa, staged at
error_model_gm12878/refs/) + comprehensive GTF (gencode.v44.annotation.gtf). Corrects the earlier roadmap note.
BLOCKED ON: ssh ControlMaster DIED (falling to password/2FA) — asked Kevin to `! ssh sherlock` to reopen (never
reopen it myself). RESUME once ssh back: (1) verify full_genome_bams @SQ naming (chr-named? no reheader);
(2) fan out refine (arm-B/Bguard) per sample with FULL GRCh38 + comprehensive GTF -> genome-wide junctions +
do-no-harm; (3) novel/non-canonical junctions genome-wide, SMA vs WT. Chunk across ~12 samples (perf ~0.24s/read).
chr5 refine (33015742) do-no-harm result not yet retrieved (ssh blip) — grab when back.

## GENOME-WIDE SUMNER FAN-OUT — downsampled first pass (2026-07-07)
Cluster = SHERLOCK (data+tools there; SCG sees only Oak). 15 WGS BAMs (8 SMA + 7 WT), chr-named, indexed,
~13.7M reads each, ~51% junction reads -> ~6.9M junction reads/sample. PERF WALL: refiner ~0.34s/junction-read
=> full-depth = ~83hr/arm/sample (~10k CPU-hr) INFEASIBLE. PI chose: DOWNSAMPLED (~5%) genome-wide first pass,
motif-blind arm, 15 parallel. Tool: rectify_guard/sumner_gw_discover.py (downsample -> refine motif-blind+guard
w/ FULL GRCh38 [chrom names register verbatim] + comprehensive gencode.v44.annotation.gtf -> classify raw vs
refined junctions annot/novel x canon/noncanon -> per-sample REVEALED novel-non-canonical). Genome
/scratch/users/kevinroy/rectify_human_validation/error_model_gm12878/refs/GRCh38.primary_assembly.genome.fa (+GTF).
IN FLIGHT: SMOKE test (sbatch 33047424, SMA_7.12 @0.5%). RESUME: ssh sherlock 'cat /scratch/users/kevinroy/
sumner_gw/.smoke_rc; tail -20 /scratch/users/kevinroy/sumner_gw/smoke-33047424.log; cat /scratch/users/kevinroy/
sumner_gw/smoke/SMA_7.12_smoke.summary.tsv'. IF SMOKE_RC=0 + sensible summary -> launch 15-sample array at
frac=0.05 (one task/sample, --array=1-15%15, 8c/20G/8h each, sample list from full_genome_bams/*.bam). Then
aggregate SMA vs WT revealed-novel-noncanon genome-wide (positive control = the SMN1-down/SMN2-up already seen
raw). IF smoke fails -> read the log, fix sumner_gw_discover.py, redeploy, re-smoke.
NOTE: standardize_chrom_name collision is a NON-ISSUE with the full genome (verified). Q1 done (153 FIX/0 HARM).
Spike-in design locked (bg decision pending). Survey done.

## GENOME-WIDE SUMNER ARRAY LAUNCHED (2026-07-07) — smoke PASSED, rich yield
Smoke (33047424, SMA_7.12 @0.5%) PASSED: pipeline works end-to-end (full GRCh38 194 contigs, 402k comprehensive
annotated junctions, chrom names verbatim). At 0.5% of ONE sample the refine REVEALED 345 novel non-canonical
junctions minimap2 flattened -> genome-wide >> chr5's handful. 15-sample ARRAY LAUNCHED: sbatch 33047654
(--array=1-15%15, frac 0.05, 8c/20G/8h). RESUME: ssh sherlock 'wc -l /scratch/users/kevinroy/sumner_gw/.panel_done'
(15 lines = done). WHEN DONE: ssh sherlock 'source ...conda; conda activate rectify; export PYTHONPATH=
/scratch/users/kevinroy/rectify_guard; python /scratch/users/kevinroy/rectify_guard/sumner_gw_aggregate.py
/scratch/users/kevinroy/sumner_gw/panel' -> per-sample yield + recurrence (>=N samples) + SMA-vs-WT specificity.
IF a task fails: read gw-33047654_<task>.log, fix sumner_gw_discover.py, redeploy, resubmit that task.
★ HONEST CAVEAT: 'revealed novel non-canonical' = DISCOVERY YIELD (real flattening-recovery + FABRICATION mixed);
no ground truth on real data. Precision is the SPIKE-IN track's job; recurrence + SMA-vs-WT specificity are the
real-data confidence signals; SMN region (chr5:70.0-70.95Mb) = built-in positive control (SMN1 down/SMN2 up seen raw).
Tools: sumner_gw_discover.py, sumner_gw_aggregate.py (committed).

## REFINER PERF INVESTIGATION (2026-07-07, task #22) — parallel with the discovery array
PI circled back to option 3: attack ~0.34s/junction-read (enables full-depth genome-wide human). Candidate
lookup RULED OUT (_candidates_near is bisect-indexed O(log pool), junction_scoring.py:630). py-spy blocked
(ptrace). IN FLIGHT: cProfile (sbatch 33048402) on 300 human junction reads -> top-by-tottime hotspot. RESUME:
ssh sherlock 'cat /scratch/users/kevinroy/sumner_gw/.prof_rc; grep -E "ms/read|tottime|_score|_apply|numba|refine_read" /scratch/users/kevinroy/sumner_gw/prof-33048402.log'. HYPOTHESIS: scoring DP (_score_hp_anchored) —
is the production numba kernel hit at human scale, or a pure-Python fallback (the del_open work noted a pure-Py
affine path)? THEN: implement fix, verify BYTE-IDENTICAL (pytest tests/test_junction_refiner.py
tests/test_hp_drift_guard.py = 46), measure speedup. TWO in-flight: discovery array 33047654 (RESUME: wc -l
/scratch/users/kevinroy/sumner_gw/.panel_done ==15 -> run sumner_gw_aggregate.py) + this profile 33048402.

## PERF INVESTIGATION (2026-07-08) — concat-DP REFUTED; full-run ins_cost under investigation
concat/2-pass DP speedup = PROVEN IMPOSSIBLE byte-identically (triple audit, committed): OLD _score_junction
computes ins_cost on the PER-K TRUNCATED substring -> splits homopolymers to exploit the NON-MONOTONIC DRS ins
table (witness A*12: old=1.7604 vs concat=8.2584). Durable audit records dev/CONCAT_DP_*.md.
PI-approved: investigate FULL-RUN ins_cost (compute on the read's full HP run, cut-independent) = sounder model
+ removes fabrication + UNLOCKS the single-pass DP. IN FLIGHT: Workflow wc2hfx0dx (durable-output baked in) —
investigate (flag _USE_FULL_RUN_INS, measure) -> 3 audits (model-correctness/impact-validity/revalidation-
completeness) -> synthesize. RESUME: await completion notif; read dev/INSCOST_INVESTIGATION.md + INSCOST_AUDIT_*.md
+ INSCOST_SYNTHESIS.md (durable). PRELIM (durable record): full-run is ANTI-FABRICATION + net-positive (fair panel
0.75% change, 35 noncanon->canon vs 6; R3 discovery preserved; guard-compatible). Auditors must check: does it drop
REAL non-canonical recoveries? does the guard tuning (m=3.0) still hold? could it flip the del_open REJECT verdict?
I review synthesis before any switch. Durable-output policy VERIFIED working (plan-first + 10 checkpoints).
IN FLIGHT #2: SUMNER genome-wide discovery array 33047654 (15/15 done). RESUME: when 15/15 ->
ssh sherlock 'conda activate rectify; export PYTHONPATH=/scratch/users/kevinroy/rectify_guard; python
/scratch/users/kevinroy/rectify_guard/sumner_gw_aggregate.py /scratch/users/kevinroy/sumner_gw/panel'.

## ins_cost A/B: REFERENCE-COLUMN vs full-run (2026-07-08) — PI chose the audit's path
Prior verdict: per-cut run-splitting = real fabrication artifact; full-run fixes it + unlocks the single-pass DP,
but full-run has a 4-6/panel canonical-DEMOTION counter-signal, and the model-correctness auditor recommended
REFERENCE-COLUMN ins indexing (charge over-call against the GENOME HP context, mirror _precompute_del_costs) as
likely-better (keeps the ~35 anti-fab wins, drops the demotions). PI: A/B ref-column vs full-run first.
IN FLIGHT: Workflow wb7dhcdhn (durable-output baked in) — build ref-column (flag _USE_REFCOL_INS) + 3-way A/B
(per-cut vs full-run vs ref-column on mix_fair_out + mix_r3b_out) -> 3 audits (refcol-model-correctness / ab-
measurement-validity / revalidation-and-coupling) -> synthesize (winner + re-validation checklist). RESUME: await
notif; read dev/INSCOST_REFCOL_BUILD.md + INSCOST_REFCOL_AUDIT_*.md + INSCOST_REFCOL_SYNTHESIS.md (durable). Key
question: does ref-column keep full-run's anti-fab wins AND eliminate the canonical-demotion losses? Then re-validate
the winner (8-item checklist in dev/INSCOST_SYNTHESIS.md: coupled-constant re-tune _CANONICAL_HP_PRIOR + hp_drift_margin,
del_open re-run, human transfer). I review synthesis before any default flip. NOTE flags uncommitted on THIS benchmark
worktree (has motif_blind); master has neither — merge-target gap in the checklist.
Also DONE this session: Sumner genome-wide discovery (219k revealed, dev/SUMNER_GW_DISCOVERY_RESULT.txt); the
ins_cost fix will clean that fabrication-dominated yield.

## FULL-RUN RE-VALIDATION IN PROGRESS (2026-07-08) — PI: run the 8-item checklist
Reference-column REFUTED; full-run is the winner (baseline committed 56addde: flags _USE_FULL_RUN_INS/
_USE_REFCOL_INS default OFF = byte-identical per-cut, 70 passed both ways). TOPOLOGY RESOLVED: whole arc
(motif_blind+guard+ins flags) is on THIS branch ONLY; master lacks it -> re-validation is apples-to-apples
HERE (all prior results computed here); master-merge = separate future step. env: RECTIFY_FULL_RUN_INS=1.
CHECKLIST STATUS:
 [x] 1 merge-target/baseline — DONE (topology resolved, flag committed 56addde).
 [x] 2a/2b narrow suite flag off+on — 70 passed both.
 [~] 2c FULL suite (not slow) flag-ON — RUNNING (bg b7zchzgow; log scratchpad/suite_flagon.log). RESUME: read that log; expect green, any fail = a test pinning old-scale scores/CIGARs -> triage.
 [ ] 4 _CANONICAL_HP_PRIOR re-confirm on full-run scale (0.5 anchor = 1 HP-DEL; DEL unchanged so anchor holds; confirm noise floor + no discovery regression fair+r3b no-guard).
 [ ] 5 hp_drift_margin re-sweep (task #16 redo under full-run; smallest zero-discovery-cost margin may DROP below 3.0). _make_arm_e.py / _sweep_refine.py pattern on mix_fair_out+mix_r3b_out.
 [ ] 6 del_open/arm-F verdict RE-RUN under full-run scale (could flip; scale x4.69). dev/arm_f_del_open_delta_faithful.patch + _make_arm_ff.py.
 [ ] 3 numba-ON DP path flag-ON on CLUSTER build (reversed-slice list->float64 -> _score_hp_dp_numba).
 [ ] 7 yeast DRS transfer (Sherlock, task #17 redo flag-on).
 [ ] 8 human ONT DRS transfer (task #18 redo flag-on) — DECISIVE / most-likely-to-move.
Flip default-ON ONLY after 2c/3/4 green + 7/8 no-regression. Each sub-agent brief MUST carry the durable-output policy.

## TABLE-FREE CONCAT-DP PORTED (2026-07-08) — flag _USE_CONCAT_DP, verifying
Ported the vectorized DP into _score_junction (junction_scoring.py), flag _USE_CONCAT_DP (env RECTIFY_CONCAT_DP,
default OFF), gated on penalty_table is None. Helper _all_suffix_scores (query-suffix reversal DP). Replaces ~60
per-k calls with 2 passes, 14.3x. Standalone proto 8000/8000 exact; suites 70 passed flag-ON / 46 flag-OFF.
IN FLIGHT: (1) in-repo 20k byte-identity harness (bg b45skijqt, log scratchpad/port_verify.log) — expect 20000/20000.
(2) end-to-end BAM diff (bg bn5ksjo6y, log scratchpad/e2e_diff.log) — refine mix_fair_out flag off vs on, expect 0
reads differ (same js/je+CIGAR). RESUME: read those logs. If both clean -> 1 adversarial audit -> flip default +
cluster deploy + re-run Sumner genome-wide/human at 14x. If any mismatch -> the fast path has a bug, keep default off.
ALSO: full-run transfers (33182801) confirm the table-free no-op premise (ssh was down; re-check flag-on==flag-off:
yeast 0.9884 / human 0.7914). Full-run re-val items 4-8 MOOT for table-free (this DP is the correct perf answer).

## TABLE-VS-FLAT ADVERSARIAL CHALLENGE (2026-07-08) — PI not convinced dropping the table is right
Concat-DP perf DONE (committed e1ed90c, 4 byte-identity proofs + audit HELD, 14.3x, flag _USE_CONCAT_DP default OFF).
PI CHALLENGE: is dropping the empirical DRS table from the re-placer (for FLAT hand-set costs sub=1.0/del_normal=1.0/
del_hp=0.5/ins=1.25, NOT calibrated) justified? DISCONNECT FOUND: production (correct_command.py:746) uses the table
+ motif_blind=False (incumbent); table-free is the NATIVE re-aligner (motif-blind+guard) research config only. What
was rejected is NARROW (arm-C penalty_score column; del_open log-odds delta — both over-shift by cheapening in-run
indels). NOT tested: table for non-drift-axis costs (sub/mismatch); table+guard vs flat+guard; flat-cost placement
ACCURACY vs calibrated; where the guard doesn't reach.
IN FLIGHT: Workflow wny725kc1 (durable-output baked in) — 3 orthogonal adversarial angles (scope-of-rejection /
flatcost-calibration / guard-blind-spots) + synthesis. RESUME: await notif; read dev/TABLE_VS_FLAT_*.md (durable).
Key Q: does ANY table use / cost calibration beat flat+guard on ground-truth recovery/accuracy? If yes, "drop the
table" is wrong. NOTE: ssh master DOWN (cluster confirmation of concat-DP transfers + deploy pending reopen).

## TABLE-VS-FLAT RESOLVED + flat-cost sweep in flight (2026-07-08)
Adversarial 3-angle attack (wzhvi1its, survived 2 API outages via durability): DROP the table = JUSTIFIED
(measured). (1) no sub row -> table can't touch mismatch; (2) HP axis = arm-C drift-grease, guard owns it;
(3) STR-del (only guard-blind table axis) RANK-EQUIVALENT to flat (compresses margins, never flips argmin winner),
NET-NEGATIVE (=rejected arm-C heuristic, greases 2.3x>forbids), NEAR-EMPTY (0.29% real Scer junctions). Verdict +
records: dev/TABLE_VS_FLAT_*.md. Flat+guard already validated real DRS (#17) + human (#18).
IN FLIGHT: flat-cost optimality sweep (bg bb9je2wn0) — the PI's 'are 0.5/1.25 right'. 7 refines (del_hp x ins around
0.5/1.25) on mix_fair_out guard-ON, D0/D2 recovery vs arm-B. Made flat costs env-tunable (RECTIFY_FLAT_DEL_HP/INS,
committed; defaults preserved). RESUME: read scratchpad task bb9je2wn0 / the printed table. Plateau near (0.5,1.25)
-> costs validated as MEASURED (not just ordinal); peak elsewhere -> flat under-tuned (still not a table argument).
PERF: concat-DP DONE (14.3x byte-identical, flag _USE_CONCAT_DP default OFF, committed e1ed90c); full-run transfers
confirmed the table-free no-op end-to-end. ssh back. NEXT after sweep: PI decision on flipping _USE_CONCAT_DP default
+ cluster-deploy + re-run Sumner genome-wide at 14x.
