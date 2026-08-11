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

### Claim B — VETTED CONCLUSION (2026-06-29, coordinator agent + 3-way fan-out + 2 advisor passes)
**The real-SIRV PLACEMENT ablation CANNOT honestly validate the per-length SHAPE.** Two
independent, empirically-checked reasons:
1. **Underpowered substrate.** Enumerated HP runs in the actual refs (SIRVomeERCCome ==
   sirv4.fasta sequences): len4=4304, 5=1447, 6=367, 7=93, 8=35, 9=7, 10=2 → only ~137
   distinct runs ≥7, **9 distinct runs ≥9**. The shape signal lives at long runs; SIRV is
   engineered to avoid them. The in-flight LRGASP `sirv.sorted.bam` won't help (same
   sequences, same ceiling — it adds reads-per-run, never distinct runs). Bonus cap: the
   law's delta SATURATES at L8 in the table (CG L8=L9=L10=+3.25), so long-run shape isn't
   even encoded.
2. **Truth-confounded ("iron triangle").** On real single reads:
   `discriminating (flat CIGAR≠law CIGAR) ⟺ a co-occurring boundary mismatch ⟺ per-read
   truth unknowable` (the data is equally consistent with in-run-del+boundary-sub vs
   last-run-base-sub+flank-del; truth_net identical). Crediting "in-run = correct" on
   ambiguous reads is a UNIDIRECTIONAL law-favoring artifact through the SAME length-scaled
   mechanism as the genuine effect → the LONG-bin `law>B0` headline is CONFOUNDED, not just
   noisy. A tie-detector that excludes ambiguous trials → expected NULL. Changing to a
   natural transcriptome fixes POWER but NOT truth (triangle persists on any real reads).

**What SURVIVES as a real-SIRV deliverable (winnable, non-circular, well-powered):** the
**over-call / false-indel-rate on CLEAN (truth_net==0) HP runs, stratified by run length**
— KNOWN truth (no indel), not subject to the in-run-credit circularity, ~6141 clean SHORT
trials in 4000 reads. A standalone real-data SAFETY result complementing the sim's clean
`false_indel_rate`. (The vetted metric backbone — per-(read,run) local window, base-count
`truth_net`, paired arms, run-clustered bootstrap, tie-detector, wrong-monotone competitor
arm — is sound and reusable.)

**The GENUINE non-circular Claim B path (deferred, multi-night):** a held-out INJECTION
simulator on natural-sequence templates with abundant long HPs, injecting length-correlated
deletions **+ boundary substitutions** at a rate **MEASURED FROM REAL SIRV** (the observed
forced-net curve — table-INDEPENDENT of the Scer `rate_mean` cost). Injected = known per-read
truth; SIRV-measured rate = independent shape. `law>B0` on the KNOWN edits ⇒ the Scer-learned
shape TRANSFERS to the SIRV-measured shape = genuine Claim B. MUST inject boundary-subs at
the SIRV-measured rate (never the Scer table) or it reproduces the original trap.

**Build prerequisite:** the C1 code (`penalty_table`/`lam`/`_homopolymer_run_len`/`del_open_delta`)
must be synced to Sherlock — `aligner_bench_live` was PRE-C1.

### OVER-CALL RESULT (2026-06-29) — real SIRV caught a bug: the INSERTION discount is harmful → GATED OFF
The surviving winnable deliverable (`scripts/benchmark/c1_real_sirv_ablation.py`) RAN on both real SIRV
BAMs (LRGASP RNA002 + SG-NEx) and produced a decisive, replicated finding. Honest control = the **sub-only
stratum** (same-length windows with 1-2 mismatches → correct alignment is provably indel-free, so any
indel is a TRUE hallucination; the `truth_net==0` ALL-clean table is contaminated by balanced real
indel pairs and is context-only). Sub-only over-call:
| BAM | bin | flat | B0 | full law | law (del-only) |
| --- | --- | --- | --- | --- | --- |
| LRGASP RNA002 | SHORT | 0.00 | 0.00 | **0.034** | 0.00 |
| | MID | 0.00 | 0.00 | **0.070** | 0.00 |
| SG-NEx | SHORT | 0.00 | 0.00 | **0.032** | 0.00 |
| | MID | 0.00 | 0.00 | **0.041** | 0.00 |
**flat and B0 hallucinate 0%; the FULL law hallucinates 3-7%, GROWING with run length; the `--zero-ins`
diagnostic isolates the cause ENTIRELY to `ins_open_delta`** (a cheap insertion lets the DP rewrite a
length-preserving substitution as a spurious D+I). Reproduced locally:
`flat [28M]` → `full-law [10M,1D,4M,1I,13M]` on a one-substitution A8 window.
**FIX:** `ins_open_delta` is now **GATED OFF by default** (`align_exon_block_global(..., ins_lengthlaw=False)`)
— the insertion discount was already flagged UNVALIDATED (the generator injects deletions only) and real
SIRV proved it actively harmful. The **deletion-only law is SAFE (0% over-call = flat) AND Claim-A-proven**;
that is what ships. Re-enable `ins_lengthlaw=True` only after the injection-simulator validates it.
Regression test: `test_ins_discount_gated_off_by_default_no_hallucination`. This is the spike-in grounding
working exactly as intended — real data caught an implementation bug the sim could not (the sim has no
length-preserving substitution-vs-indel ambiguity at scale).
