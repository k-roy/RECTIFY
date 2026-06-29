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

**Pre-committed null:** on flat-in-L sim, **full-law ≈ B=0 is the CORRECT result**,
reported as "length-shape inconclusive on flat-in-L sim → deferred to SIRV/RNA004,"
NOT a win and NOT a failure. Reaching for a rate_mean-scaled generator to make
full-law beat B=0 is the tell that the win is being manufactured.

## Honest deliverable

"C1 implemented; an in-run calibrated gap discount improves indel placement vs a
matched HP-aware baseline, with false-indel ≈0 on clean runs, boundary_sub not
regressing, and Cat3 byte-identical; the per-length SHAPE validation is gated to
real SIRV/RNA004 data." That is Claim A — real and shippable. Claim B is multi-night.
