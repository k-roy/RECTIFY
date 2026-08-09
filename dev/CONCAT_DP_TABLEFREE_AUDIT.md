# CONCAT-DP TABLE-FREE FAST-PATH — ADVERSARIAL BYTE-IDENTITY AUDIT

**Auditor:** single adversarial auditor (worktree agent-a25a2c1e784ad37dc)
**Target:** `_all_suffix_scores` + fast-path block in `_score_junction`
(`rectify/core/splice/junction_scoring.py`), gated on
`_USE_CONCAT_DP (env RECTIFY_CONCAT_DP=1) AND penalty_table is None`.
**Goal:** find ONE input where fast path (flag ON, table=None) != legacy
path (flag OFF, table=None), EXACT float. Attack weak spots of random testing.
**Constraint:** numba NOT installed (confirmed — import errors), so legacy uses
the pure-Python one-row DP in `_score_hp_anchored`. Do NOT edit junction_scoring.py.

## Method
- In-process toggle: set `os.environ["RECTIFY_CONCAT_DP"]`, `importlib.reload`
  the module, call `J._score_junction(...)` on IDENTICAL inputs, assert `==`
  (exact float, `repr()` compare + `==`).
- Harness reloads the module fresh for OFF and ON on each identical input to
  avoid stale module-level flag capture.

## Attack plan / case classes
1. **MECH1 boundary/reversal indexing** — k=L exclusion; reverse(rescue)/
   reverse(ref)/reverse(del_costs) alignment; t2_vec[L-k] index; t2(0)=0.
   Construct L=1, L=2 rescues where the k index arithmetic is fragile.
2. **EDGE cases random rarely hits:**
   - L=1 rescue.
   - rescue longer than ref_exon2 / ref_intron_end (R < L).
   - intron_end within _BUF of a genome END → ref windows TRUNCATE to <_BUF,
     possibly R=0. Does `_all_suffix_scores` handle R=0 (all-insertion branch)?
   - intron_end == 0 (ref_intron_end empty → R=0 on the rev side).
   - intron_end == gs (ref_exon2 empty → R=0 on fwd side).
   - rescue with N / non-ACGT bases (uppercase already applied).
   - all-homopolymer rescue AND all-homopolymer ref.
   - genome all one base.
3. **del_costs edges** — STR/HP context at reversed-window boundaries
   (ref_genome_rev=True path); del_hp vs del_normal transitions at ref ends.
   Since fast path and legacy SHARE the same precomputed del_costs arrays,
   this tests whether `_all_suffix_scores`'s `dc = del_costs[::-1]` reversal
   matches the legacy per-k `precomputed_del_costs` usage.
4. **GATING** — fast path fires ONLY when (flag ON and table is None). With a
   penalty_table it MUST fall through to legacy: verify flag ON + table ==
   flag OFF + table. And flag OFF is completely inert.
5. **FP order** — 1-ULP divergence between 2-pass and per-k accumulation. ins,
   del, sub all dyadic (1.25, 1.0, 0.5). Argue whether summation ORDER diverges.

## Checkpoints (appended as each class completes)

- **CLASS 3 (asymmetric non-uniform del_costs, HP dropout) -> IDENTICAL.**
  3 targeted cases (fwd A-run, rev T-run, both-sides) with rescue HAVING FEWER
  run bases than genome so del_hp=0.5 actually enters. Non-trivial scores
  6.5 / 3.0 / 3.5 (proves del_hp exercised + del_costs[::-1] reversal exercised).
  off==on exact for all 3. This is the case random ACGT genomes never hit
  (they give uniform del_costs=1.0 -> reversal is a no-op).

- **del_costs non-uniformity PROVEN:** class-3 fwd window del_costs_fwd =
  [0.5]*8 + [1.0]*... (unique {0.5,1.0}); reversed != forward = True. The
  del_costs[::-1] reversal in _all_suffix_scores was genuinely exercised
  (NOT the random-ACGT no-op).

- **CLASS 1 (MECH1 boundary/reversal) -> IDENTICAL.** L=1 (k=0 only, t2 unused),
  L=2, L=2 HP tie, L=3 with-N. off==on exact (1.0/1.0/2.0/1.0).

- **CLASS 2 (truncation/edge) -> IDENTICAL (all 8).** fwd R<L; fwd R=0
  (intron_end==gs -> ref_exon2 empty, score 1.25 = all-insertion); rev R=0
  (intron_end==0); rev R<L; N + lowercase query; all-A genome + all-A rescue
  (0.0); all-A genome + diverse rescue (7.0); junk/non-ACGT genome. off==on exact.

- **CLASS 4 GATING -> CORRECT.** Real HpPenaltyTable (from penalty_scores.tsv,
  NON-dyadic values e.g. 1.1615, 0.1938): flag ON + table == flag OFF + table
  EXACT (3 cases). Fast path correctly falls through to legacy with a table.
  Flag-OFF inertness: env unset / "0" / "2" all -> _USE_CONCAT_DP False (only
  "==1" enables); results identical. Fast path cannot touch the production
  table path.

- **BROAD SWEEP (asymmetric-HP-seeded random) -> IDENTICAL.** 6000 genomes each
  SEEDED with HP runs (len 4-10) on BOTH flanks of intron_end + random rescue
  dropout — the exact regime plain-random ACGT misses. 0/6000 DIVERGE.

- **CLASS 5 FP-ORDER -> NO DIVERGENCE POSSIBLE (by argument).** All live DP costs
  dyadic: ins=1.25, del in {1.0,0.5}, sub=1.0; scores bounded ~30*1.25=37.5.
  Every partial sum is exactly representable in float64 -> zero rounding occurs
  anywhere, so accumulation ORDER (2-pass vs per-k) is irrelevant (exact
  arithmetic is associative): 0 ULP GUARANTEED, not merely observed. The
  hand-rolled `min` (diag<above else above; if left<curr) vs min(diag,above,left)
  differ only in WHICH tied arg is selected, never the numeric value.

## RESIDUAL RISK (one)
`ins=1.25` and `sub=1.0` are HARDCODED INDEPENDENTLY as defaults in BOTH
`_all_suffix_scores` (params `ins=1.25, sub=1.0`) AND `_score_hp_anchored`
(params `ins=1.25, sub=1.0`), with NO shared named constant. Byte-identity holds
today only because they COINCIDE. A future edit to one path's default (e.g.
retuning the flat insertion cost) silently breaks identity, and since the flag
defaults OFF the divergence would ship unnoticed. Recommend a shared module-level
constant (e.g. `_FLAT_INS = 1.25`, `_FLAT_SUB = 1.0`) referenced by both. This is
a coupling/maintenance hazard, NOT a current-code divergence.

- **L=0 EARLY-RETURN + MID-K t2-participation -> IDENTICAL.** q_split past query
  end and empty query both hit the L==0 early return (0.0,0) before the fast path
  — match. Plus a 4000-case dense-HP sweep engineered so the winning k is
  frequently MID (t2 genuinely enters the argmin, not just k=0) -> 0 DIVERGE.
  This confirms t2_vec[L-k] indexing is correct when it actually decides the score.

## VERDICT
Byte-identity HELD across every attacked case class (targeted MECH1/edge/del-costs
+ 10000 asymmetric/dense-HP-seeded random + gating). No witness found. The identity is
exact by construction (alignment-reversal duality; all-dyadic costs). Gating is
sound: the fast path is provably inert unless (RECTIFY_CONCAT_DP==1 AND
penalty_table is None). Recommend keeping default OFF and, before any promotion,
factoring the flat ins/sub constants into shared names to prevent silent drift.

STATUS: COMPLETE.
