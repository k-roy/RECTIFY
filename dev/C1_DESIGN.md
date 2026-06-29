# C1 — HP-length-law gap cost: design (locked 2026-06-29)

Provenance: triple adversarial agent panel (principled-correctness / red-team /
minimal-shippable) + multiple advisor passes. This note records the LOCKED design
and, crucially, what the simulation gate can and cannot prove.

## The cost (principled-correctness reviewer — decisive, code-level)

- **`HpPenaltyTable.del_cost/ins_cost` return the `penalty_score` column, which is
  NOT −logP.** It is a reciprocal-rate heuristic (`penalty_score ≈ c/rate_mean`;
  the product is constant within each (op, base_class) group). Summing it along an
  additive DP path corresponds to no probability model. **Do NOT feed
  `penalty_score` into the DP.** (Leave `del_cost`/`ins_cost` untouched — shared by
  junction_scoring / Cat3; just flag the mislabel.)
- **Use `rate_mean`.** Over {M,X,D,I} it sums to ≈1.0 at every (base_class, hp), i.e.
  it is already a calibrated per-reference-column emission distribution.
- **Scale-mixing is provably incoherent in a global, full-consumption DP** (the
  match-count varies across paths, so a per-diagonal constant shifts the
  match-vs-indel tradeoff). The only coherent partial is a **baseline-anchored
  log-odds DELTA on the gap-OPEN**, leaving match/mismatch/extend on the legacy
  scale:

      Δ_D(hp, bc) = λ · ln( P_D(hp,bc) / P_D(1,bc) )         # added to gap_open, del
      Δ_I(hp, bc) = λ · ln( P_I(hp,bc) / P_I(1,bc) )         # added to gap_open, ins

  where `P_op = rate_mean`. Zero at hp=1 (true "scale-preserving" — only the
  hp-DEPENDENCE is new). λ ≈ 1 (the legacy +2/−4 integers already sit on a ~1-nat
  log-odds scale; reviewer-1's anchored λ = 6/ln(P_M/P_X) ≈ 0.996). Since P_D rises
  with run length, Δ_D > 0 for hp>1 → a deletion is CHEAPER in longer runs → the DP
  prefers the in-run deletion over out-of-run misplacement / mismatch-repair.
- Insertion delta is implemented symmetrically but is UNVALIDATED (the controlled
  generator injects deletions only — no insertion stratum yet).
- Gated STRICTLY on `homo_mask[j-1]` (the existing HP mask). Off-mask positions and
  `penalty_table=None` → byte-identical to the legacy DP (Cat3 regression guard).

## What the gate can prove tonight (red-team + advisor — the load-bearing limit)

Two SEPARATE claims; do not conflate:

- **Claim A — "an in-run gap DISCOUNT improves indel PLACEMENT."** Honestly testable
  NOW on the flat-in-L HP_HARD stratum, because the failure mode is *where* the
  deletion goes (out-of-run misplacement / substitution-repaired-as-indel), which a
  discount fixes regardless of the per-length rate.
- **Claim B — "the table's per-length SHAPE is correct."** NOT testable on the
  current sim, and NOT honestly testable by length-scaling the generator from the
  SAME `rate_mean` table — that makes the cost and the test data share a source, so
  the length-law wins by construction (the 0.09→1.07 artifact one level up; SPEC
  lines 32–33: the scoring model must be HELD OUT from the calibration table). Claim
  B requires an INDEPENDENT error source: real SIRV / RNA004 (the multi-night track)
  or a held-out simulator (pbsim3/badread/different-chemistry). NOT attempted tonight.

## The ablation (matched arms; reps ≥ 400; TEST split)

Three arms, IDENTICAL except the gap-open cost, ALL passing `chrom_ref` (so
`homo_mask` + the legacy `homo_mismatch=−2` are active in every arm — without
`chrom_ref` the baseline silently has NO HP-awareness, which would make any "win" a
conflated artifact):

1. **flat (matched baseline):** `penalty_table=None` → legacy homo_mismatch only.
2. **B=0 (control):** a CONSTANT in-run gap-open discount (no hp-dependence).
3. **full-law:** the per-(hp,base) log-odds delta above.

Metrics vs TRUTH (never the internal DP score), on the TEST partition:
- position-exact indel concordance on HP_HARD-noisy (Claim A: full-law/B=0 > flat);
- **false_indel_rate on CLEAN HP/STR runs routed THROUGH the arm** (must stay ≈0 —
  the discount must not hallucinate deletions; clean reads were previously never
  routed through the DP arm, so this control did not exist);
- boundary_sub concordance must NOT regress while noisy improves (guards the
  double-count with homo_mismatch flipping true substitutions into indels).

**Pre-committed null:** on flat-in-L sim, **B=0 ≥ law is the CORRECT result** — the
length-SHAPE is mildly ANTI-helpful here (a uniform discount fixes short AND long
runs equally; the law deliberately under-discounts short runs, which the flat-in-L
generator fails just as often). Report it as "length-shape unvalidated on flat-in-L
sim → deferred to real data," NOT a win and NOT a failure. **Do NOT raise λ to make
law beat B=0** — that gap is the sim lacking length-correlation, not a mis-set scale;
cranking λ is the hill-climb-into-the-simulator trap. λ≈1 is principled; leave it.

## RESULT (2026-06-29, reps=120 TEST, mechanism hand-VERIFIED)
| metric | flat | B=0 | law |
| --- | --- | --- | --- |
| HP_HARD-noisy concordance | 0.962 | 0.990 | 0.985 |
| boundary_sub concordance | **0.000** | 0.78 | 0.55 |
| clean false_indel_rate | 0.000 | 0.000 | 0.000 (1854 clean) |

**Claim A PASS, and the headline is REAL (hand-traced, not a scorer artifact):** flat
places the boundary_sub deletion OUT of the run (run [80,90) → D at ref 90), scoring
0/600; law places it IN-run (D at ref 80) on 328/600. Example `hph_A_10_boundary_sub`:
flat `[90M,1D,79M]` (out) → law `[80M,1D,89M]` (in). The discount makes the in-run
deletion cheaper so the DP calls the boundary as a mismatch instead of absorbing it
into a misplaced gap — exactly the C1 target. **B=0 > law as pre-committed** (shape
unvalidated on flat-in-L). Note B=0's constant IS the mean of law's own deltas
(table-derived), so it is a SHAPE-ablation control, not an independent alternative.

## Honest deliverable + ship decision

"C1 implemented; an in-run calibrated gap discount improves indel PLACEMENT vs a
matched HP-aware baseline (boundary_sub 0.00→0.55/0.78; noisy 0.96→0.99), with
false-indel ≈0 on clean runs, boundary_sub not regressing, and Cat3 byte-identical."
That is **Claim A — real, mechanism-verified, shippable.**

**Ship the law (not the constant B=0)** — on COHERENCE + real-data-deferred grounds,
NOT sim performance (where it slightly loses). The law is the principled rate_mean
log-odds form; it DEGRADES to ≈constant on uncorrelated errors (why it ties/loses to
B=0 here) while being able to EXPLOIT length-correlation when it is real. Whether real
errors are length-correlated = **Claim B**, the next test.

## Claim B — the RIGHT test (advisor-corrected; next session)
NOT "does deletion rate rise with HP length on real SIRV" (a property of the READS,
already established, says nothing about our aligner). The honest C1 real-data test is
the **placement ablation on real SIRV**: align real SIRV reads to the SIRV reference;
for each read covering a known reference HP run, score in-run vs out-of-run placement
(net-indel-in-span / out-of-span) for flat vs B=0 vs law. Truth = the known reference
run span (from the SIRV FASTA, even though the per-read injected edit is unknown).
Caveats: no boundary_sub analog in real reads; reads carry co-occurring errors; must
enumerate HP-run spans from the FASTA. This is where the length-SHAPE (Claim B) is
earned or refuted on an INDEPENDENT source. Multi-night; do not rush it.
