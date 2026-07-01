# C1 Claim B — held-out injection simulator (non-circular length-SHAPE test)

Status: **SCAFFOLD + MECHANISM RESULT (Claim B UNVALIDATED — pre-committed null
confirmed on synthetic substrate; genuine Claim B blocked on the SIRV curve).**
Agent C1-ClaimB, worktree agent-a00cde628f2d37a20, off a25 tip 7a05a2c.
This note + `scripts/benchmark/c1_claimb_injection.py` are the durable deliverable.

## Goal
Validate that the C1 **per-length SHAPE** (learned from the Scer `rate_mean` penalty
table) TRANSFERS to an INDEPENDENT error source. This is **Claim B** — distinct from
Claim A ("an in-run discount improves placement", already PASS on the flat-in-L sim,
`dev/C1_DESIGN.md` §RESULT). The flat-in-L `c1_ablation.py` sim CANNOT test Claim B
(the cost and the test data would share a source → law wins by construction). Real
SIRV placement CANNOT either (underpowered long runs + truth-confounded "iron
triangle"; `dev/C1_DESIGN.md` §"Claim B — VETTED CONCLUSION").

## The non-circular design (dev/C1_DESIGN.md lines 146-152)
A held-out INJECTION simulator:
1. **Templates**: natural-sequence windows with ABUNDANT LONG homopolymers.
2. **Inject** length-correlated deletions **+ boundary substitutions** at a rate
   **MEASURED FROM REAL SIRV** — table-INDEPENDENT of the Scer `rate_mean` cost.
3. Injected edits = KNOWN per-read truth.
4. Ablate **flat vs B0 vs law** (mirror `c1_ablation.py`:
   `align_exon_block_global(..., penalty_table, lam)`, same 3 arms).
5. **law > B0 on the KNOWN edits ⇒ the Scer-learned shape transfers ⇒ Claim B.**

Non-circularity hinges on TWO independent sources:
- the DP cost = Scer `penalty_scores.tsv` (`rate_mean`);
- the injection rate = SIRV-observed length→deletion curve.
If law>B0 here, the shape the DP learned from Scer matched the shape SIRV actually
has — that is genuine transfer, not fitting-to-self.
**MUST inject boundary-subs at the SIRV-measured rate, NEVER from the Scer table**
(else it reproduces the original circular trap).

## Template source (DECIDED)
Synthesized long-HP loci via the VETTED `controlled.gen_hp_hard_stratum` (planted
runs len 4/6/8/10/12, 4 bases, unique non-run flanks). Advisor-endorsed synthesis
over yeast-genome extraction (yeast is HP-poor at long lengths — the same underpower
that killed the real-SIRV substrate, C1_DESIGN line 124). Reusing the vetted
generator also minimizes bug surface and keeps the construction identical to the
one behind Claim A's real headline.

## SIRV-rate source — NOT AVAILABLE (the honest blocker)
Genuine Claim B needs the injection rate to be a SIRV-MEASURED per-HP-run-length
del-rate curve. **That curve was never recorded as numbers.** The SPEC §SIRV
ABSOLUTE-TRUTH RESULT records error-STRUCTURE stats only (`indel_run>=2`=0.36,
overdispersion, autocorr r~0.20, gap-excess) — NONE is a del-rate-by-runlen shape.
Using `indel_run>=2` (an indel-LENGTH stat) to fake a per-run-LENGTH rate would
reintroduce circularity by misrepresentation (advisor). So this run uses the
FLAT-in-L `K_DIST` from `controlled.py` and is labeled SYNTHETIC, never
SIRV-measured. Getting the real curve = the deferred multi-night SIRV-BAM track
(see RESUME).

## Method (arms + metrics — mirror c1_ablation.py)
- Arms: `flat` (penalty_table=None), `B0` (ConstantDiscount = mean law delta),
  `law` (HpPenaltyTable from Scer penalty_scores.tsv). lam=1.0.
- Align each injected read to its template with `align_exon_block_global(seq, ref,
  chrom_ref=ref, penalty_table=pt, lam=lam)`.
- Score vs KNOWN injected truth: position-exact indel concordance, stratified by
  the reference run length of the edited run (the SHAPE axis).
- **Verdict**: law > B0 on the LONG-run bins (with the injection rate length-
  correlated per SIRV) ⇒ Claim B CONFIRMED. law ~ B0 ⇒ still unvalidated. law < B0
  ⇒ anti-helpful (shape mis-transfers).

## Result (2026-06-30, reps=60/locus, seed=7 — MECHANISM, not Claim B)
Substrate = the VETTED `gen_hp_hard_stratum` boundary_sub construction (the exact
generator behind Claim A's 0.00->0.55/0.78 headline), scored flat/B0/law by run
length L. DP cost = Scer `penalty_scores.tsv`; nothing calibrated on the sim.

| run_len L | n   | flat   | B0     | law    | law-B0 |
| --- | --- | --- | --- | --- | --- |
| 4   | 120 | 0.0000 | 0.7917 | 0.0000 | -0.7917 |
| 6   | 120 | 0.0000 | 0.8000 | 0.3917 | -0.4083 |
| 8   | 120 | 0.0000 | 0.8000 | 0.8000 |  0.0000 |
| 10  | 120 | 0.0000 | 0.8167 | 0.8167 |  0.0000 |
| 12  | 120 | 0.0000 | 0.7333 | 0.7333 |  0.0000 |
| **ALL (flat)** | 600 | 0.0000 | **0.7883** | **0.5483** | -0.2400 |
| reweight-L (up-weight long) | — | 0.0000 | 0.7833 | 0.6429 | -0.1404 |

### Second stratum — `noisy` (the GRADED case, `--stratum noisy`, reps=60)
The `boundary_sub` case above is a pure THRESHOLD (flat=0 everywhere). The `noisy`
sub-case is GRADED (flat near ceiling in Claim A), so it is the one place a
per-length law-vs-B0 signal could actually surface. It does NOT:

| run_len L | n | flat | B0 | law | law-B0 |
| --- | --- | --- | --- | --- | --- |
| 4  | 120 | 0.9667 | 0.9833 | 0.9667 | -0.0167 |
| 6  | 120 | 0.9583 | 0.9833 | 0.9667 | -0.0167 |
| 8  | 120 | 0.9750 | 0.9917 | 0.9917 |  0.0000 |
| 10 | 120 | 0.9500 | 0.9833 | 0.9833 |  0.0000 |
| 12 | 120 | 0.9750 | 1.0000 | 1.0000 |  0.0000 |
| **ALL** | 600 | 0.9650 | **0.9883** | **0.9817** | -0.0067 |

flat=0.965 reproduces Claim A's noisy flat (0.962). law again TIES B0 at L>=8 and
slightly loses at short L — law NEVER exceeds B0 on the graded stratum either.
**So BOTH the threshold (boundary_sub) AND the graded (noisy) sub-cases fail to
reveal any beneficial length-shape** — the synthetic-sim-cannot-validate-Claim-B
conclusion is airtight, resting on both strata, not one.

**Findings:**
1. **Harness-sanity PASS.** Aggregate B0=0.7883 >= law=0.5483 reproduces Claim A
   (0.78 / 0.55) to 2 decimals. The harness is NOT rigged to favor law.
2. **law NEVER exceeds B0 at ANY length** — it TIES at long runs (L>=8, where its
   log-odds delta saturates high: delta@L8(A)=2.60 ~ B0 const 2.42) and LOSES at
   short runs (L=4: law=0.00, collapsing to flat, because its delta there is below
   the placement threshold). Even up-weighting long runs keeps law<B0 (0.64<0.78).
3. **WHY (a mechanism the flat-in-L sim exposes):** this placement task is
   THRESHOLD-like, not graded — flat (no discount) fails everywhere (0.00); ANY
   discount clearing the threshold succeeds (~0.80). B0's constant 2.42 clears it
   at every L; law only clears it at L>=8. So a bigger-when-longer discount buys
   nothing once the threshold is met, and a graded SHAPE cannot be rewarded here.
   This is a SECOND, independent reason (on top of "K_DIST is flat-in-L") that the
   synthetic sim cannot validate Claim B.

**Verdict: pre-committed NULL confirmed (B0 >= law on flat-in-L), with a per-length
mechanism.** NOT a failure of the law and NOT Claim B — the length-SHAPE is
unvalidatable on this substrate. lambda was NOT cranked (C1_DESIGN 75-77).

Reproduce: `PYTHONPATH=. python scripts/benchmark/c1_claimb_injection.py --out
/tmp/cb_run --reps 60 --reweight-long --penalty-table
rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv`

## RESUME (next session — the genuine Claim B, multi-night)
The synthetic-injection track is EXHAUSTED for Claim B: any flat-in-L placement
task is threshold-like and cannot reward the shape. Genuine Claim B needs a
graded, length-CORRELATED error source, i.e. the SIRV-measured per-run-length
del-rate curve, which is NOT yet recorded. Concrete next step:
1. On Sherlock, from the SIRV BAMs `c1_real_sirv_ablation.py` already opens
   (LRGASP RNA002 + SG-NEx), enumerate reference HP runs and measure
   P(net deletion) and mean deletion depth as a function of reference run length
   -> the SIRV del-rate-by-runlen curve (table-INDEPENDENT of Scer).
2. Feed that curve as the injection rate here (replace the flat K_DIST with a
   length-correlated k-draw), keeping boundary-subs at the SIRV-measured rate.
3. Re-run flat vs rising; law>B0 in the LONG bin under the SIRV curve = Claim B.
Also test a GRADED placement metric (partial credit for how close the deletion
lands to the run) so the shape has a continuous axis to be rewarded on.

CAVEATS the RESUME must inherit (do NOT execute without accounting for these):
- **SIRV substrate is underpowered exactly where shape lives:** the SIRV refs have
  only ~9 distinct runs >=9 (C1_DESIGN line 124), so a measured del-rate-by-runlen
  curve is noisy/absent at long L — the very regime the law's shape distinguishes.
- **The law's delta SATURATES at L8** in the table (CG L8=L9=L10=+3.25), so there
  is little shape to reward above L8 regardless of the injection — consistent with
  the L>=8 ties observed here. Any Claim B win must come from the L5-8 gradient.
- **B0 is a single scalar across AT and CG**, so the measured law-B0 gap captures
  the full (run_length x base_class) shape, not length alone (code-review nit); a
  base-class-matched B0 would isolate length more cleanly if that distinction matters.

## Files
- `scripts/benchmark/c1_claimb_injection.py` — the simulator + ablation.
- `dev/C1_DESIGN.md` — the locked design + Claim A result + the Claim B path spec.
- `scripts/benchmark/c1_ablation.py` — the Claim A harness this mirrors.
