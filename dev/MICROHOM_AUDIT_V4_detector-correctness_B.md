# MICROHOM AUDIT V4 — detector-correctness — AUDITOR B (independent)

**Task:** detector-correctness. **Auditor:** B. **Commit:** 05664bc. **Mode:** READ-ONLY.
**Lens:** are `_move_microhomology` / `_frac_match` / `_effective_veto_margin` (junction_refiner.py)
INCORRECT for some input — scoring WRONG in a way that vetoes a REAL move (discovery loss) or
lets a fab through (specificity loss)?

Independent of Auditor A. Fresh harness. My own verdict.

---

## Guard mechanics (verified in rectify/core/splice/junction_refiner.py @ 05664bc)

- `_frac_match(a,b)` (L491): fraction of equal positions; **`x == y`** with no ACGT filter →
  **N==N counts as a match**, lowercase==lowercase counts as a match, but 'a'!='A'.
- `_move_microhomology(seq,ns,ne,js,je)` (L499): for the acceptor leg (je!=ne) compares
  `seq[lo:hi]` vs `seq[hi:hi+k]`; donor leg (js!=ns) compares `seq[lo-k:lo]` vs `seq[lo:hi]`.
  Returns **max over the two legs**. Guarded for off-end (`hi+k<=n`, `lo-k>=0`).
- `_effective_veto_margin(hold,eff,cap)` (L535): `max(hold, min(eff,cap))` when `cap>0 and eff>hold`,
  else `eff`.
- Call site (L802-841): `eff_margin = hold_margin`; if move lands in HP (`_hp_run_across>0`) add
  `hp_drift_margin`; if `_move_microhomology >= microhom_threshold` (default 0.5) add
  `microhom_drift_margin`; then `veto_margin = _effective_veto_margin(...)`; **veto fires when
  `best_score_cmp > incumbent_score - veto_margin`** (score LOWER=better → larger veto_margin ⇒
  move held at incumbent ⇒ discovery loss if the move was real).
  A false-POSITIVE microhomology (helper ≥0.5 when it should be <0.5) ⇒ margin added ⇒ real move vetoed.

## Named priority risks
- **A5**: `_frac_match` scores N==N (any non-ACGT==non-ACGT) as MATCH → phantom microhomology ≥0.5
  can veto a REAL move.
- **A8**: max-over-both-boundaries → benign transition on one boundary masked by unrelated repeat on
  the other trips ≥0.5, vetoes a move that should be spared.

---

## CHECKPOINT 1 — orientation done, harness reviewed, genome alphabet measured

- Prior attempt-1 left a scratch harness (harness.py) + unit_results.json but **NO durable record**.
  Reviewed for orientation; building my OWN harness fresh (not resuming) per mandate.
- **A5 unit-level CONFIRMED fires** (prior + will re-derive independently):
  `_frac_match("NNNNNN","NNNNNN") = 1.0` (correct N-excluding = 0.0);
  tip case `_frac_match("ACNTTT","ACNGGG") = 0.5` TRIPS threshold, real-base agreement = 2/6 = 0.333
  does NOT trip. So a genuine sequence TRANSITION (0.333, should be a discoverable move) is flagged as
  microhomology purely because of one N==N pair.
- **A8 unit-level fires** (donor repeat 1.0 masks benign acceptor 0.0 → max=1.0).
- **`_effective_veto_margin` truth table: ALL PASS** (10 regimes incl. hold>cap→hold, eff==hold no-op,
  cap≤0 disabled, cap≥eff inert). Invariant hold≤got≤eff holds everywhere. → EVM CLEAR.

### Decisive real-world facts measured
- **Bundled S. cerevisiae reference (S288C R64-5-1) alphabet = pure ACGT, ZERO N**
  (A:3766349 C:2320576 G:2317100 T:3753080; total 12,157,105; non-ACGT fraction = 0.0).
  ⇒ **On the PRIMARY (yeast DRS) substrate, A5 CANNOT fire — bounded to zero.**
- The transfer/validation tests (#18) use **human ONT DRS** (GRCh38), which contains ~150 Mb of N
  (assembly gaps/centromeres/telomeres). That is the substrate where A5 could become reachable —
  IF the genome loader preserves N (`upper("N")=="N"`). **Loader trace = next step.**

## Verdict framing (per advisor, split A5 vs A8)
- **A5 IS a genuine detector-correctness DEFECT** (wrong number: N-vs-N is not a tandem repeat).
  Whether it is a **HOLDING** fault requires BOTH: (1) loader preserves N into `genome_seq`, AND
  (2) an end-to-end decision-flip (real move vetoed by N-inflation that N-correct handling emits).
- **A8 is NOT a second fault**: `max()` faithfully reports the max genomic microhomology across the
  two shifted boundaries of an ATOMIC move — a correct number. "Masks a benign acceptor" = read-blind
  genomic context, which is the V2 STRUCTURAL claim already CLEARED + explicitly out of scope. A8 blocks
  nothing (acceptor-only move stays independently reachable when it scores best).
- Soft-mask/case (lowercase mismatch) is real but OPPOSITE sign (false-NEGATIVE → guard misses drift),
  gated on the same loader fact.

**STATUS: in progress — pending (a) loader trace, (b) e2e decision-flip. Verdict not yet final.**

---

## CHECKPOINT 2 — LOADER TRACE DONE (decisive: N IS preserved into genome_seq)

Traced the genome dict that reaches `refine_bam_junctions(genome=...)`:
`correct_command.py:697-698` → `rectify.utils.genome.load_genome` →
```python
genome[record.id] = str(record.seq).upper()   # genome.py:152 / :155
```
Empirically verified (tiny FASTA `ACGTNNNNNNNNNNacgtACGT...`):
- **N SURVIVES** the `.upper()` (14 N's in → 14 N's out; `"n".upper()=="N"`). ⇒ **A5 is REACHABLE**:
  on a human/GRCh38 substrate (the transfer test #18, ~150 Mb of literal N at gaps/centromeres/telomeres),
  `genome_seq` passed to `_move_microhomology` contains N. `_frac_match` then scores N==N as a match.
- **Soft-masking (lowercase) is DESTROYED** by `.upper()` → the case-sensitivity/false-NEGATIVE axis
  (lowercase mismatch making `_frac_match` MISS microhomology) is **NOT reachable** through the real
  pipeline. Closed CLEAR-by-loader. (Only a caller that bypassed `load_genome` and fed lowercase would
  hit it; the sole production loader uppercases.)
- All refiner call sites (correct_command single-shot + streaming `_refine_read_batch` at ~1479/1548/1861)
  pull `genome.get(chrom,'')` from this same dict — no alternate N-stripping path.

⇒ The advisor's "IF (loader preserves N)" gate is **RESOLVED TRUE** for the human-transfer substrate.
Remaining gate for HOLD: the **end-to-end decision-flip** (does N-inflation actually veto a real move).

---

## CHECKPOINT 3 — END-TO-END DECISION-FLIP LANDED (A5 = HOLDING FAULT)

Harness: `.../detector-correctness_B/b_e2e_independent.py` → `b_e2e_independent.json` (my own,
built fresh; NOT attempt-1's). Drives the real `refine_read_junctions` (motif_blind, pool =
{incumbent (D,A), drift (D,A+6)}), read constructed so exon2 honestly encodes X-ops at the
incumbent and the query supports the +6 acceptor move (delta_improve > 0 — a GENUINE move).
Only `microhom_drift_margin` (0.0 off → 8.0 on) changes between the two runs. Profile counters
prove causation.

ACCEPTOR-leg grid (accepted `je`; A=170, A+6=176):

| case | flank(seg_true·seg_drift) | actual_mh | REAL-base_mh | guard OFF | guard ON | flagged | vetoed | FAULT |
|---|---|---|---|---|---|---|---|---|
| acc_A pos-control (real repeat) | ACGTAC·ACGTAC | 1.0 | 1.0 | 176 moved | 170 held | 1 | 1 | correct (real mh≥0.5) |
| acc_B **pure-ACGT control** | ACGTGA·TCAGCT | 0.167 | 0.167 | 176 | **176 moved** | 0 | 0 | spared (isolates N) |
| acc_C N-phantom partial | ACNTGA·ACNCTT | **0.5** | 0.333 | 176 | **170 held** | 1 | 1 | **TRUE** |
| acc_D N-phantom **full** | NNNNNN·NNNNNN | **1.0** | **0.0** | 176 | **170 held** | 1 | 1 | **TRUE** |

**Decisive isolation:** acc_B and acc_D have the SAME move geometry and same read shape; the ONLY
difference is the flank sequence. acc_B (genuine ACGT transition, real mh 0.167) is EMITTED with the
guard on; acc_D (flank = N, real mh **0.0** — a maximal genuine transition the guard is DESIGNED to
spare) is VETOED, and `microhom_drift_flagged=1`, `move_margin_vetoes=1` confirm the veto is the
MICROHOMOLOGY path caused by N==N inflating `_move_microhomology` 0.0→1.0. A correct (N-excluding)
detector returns 0.0 for acc_D and 0.333 for acc_C → would NOT veto either. One-line fix: in
`_frac_match`, require `x == y and x in "ACGT"` (or `.upper() in ACGT`).

**Near-tie cap:** `drift_near_tie_cap=2.0` RECOVERS the move (accepted je → 176) in every N case —
confirming the cap structurally bounds this loss. BUT the cap DEFAULTS to 0.0 (disabled), so in the
shipped default the phantom veto STANDS. The cap being off-by-default is exactly the V2 established
scope; the A5 detector wrongness is orthogonal and independently real.

**VERDICT (A5): HOLD.** A5 is a genuine detector-correctness fault that flips a real discovery
decision on the human/GRCh38 substrate (N present, loader preserves it). Reproducing case = acc_D
(+ acc_C). A8/EVM below.

**STATUS: A5 HOLD proven. Remaining: donor-leg + minus-strand robustness (not verdict-changing);
A8 disposition; EVM close; realism bound.**

---

## CHECKPOINT 4 — UNIT MATRIX + A8/EVM disposition (independent re-derivation)

Harness `.../detector-correctness_B/b_units_and_donor.py` → `b_units_and_donor.json`.

**A5_unit** (re-derived from scratch):
- `_frac_match("NNNN","NNNN") = 1.0` (correct N-excluding = 0.0).
- `_frac_match("RYRY","RYRY") = 1.0` — the bug is NOT N-specific: ANY non-ACGT==non-ACGT
  (IUPAC ambiguity R/Y/K/M/S/W, which DO occur in some assemblies) phantom-matches. Broadens A5.
- tip `_frac_match("ACNTTT","ACNGGG") = 0.5` vs real 0.333.
- `_frac_match("acgt","ACGT") = 0.0` (case-sensitive miss) — but `load_genome` uppercases, so this
  false-NEGATIVE axis is UNREACHABLE via the pipeline (documented, closed).

**A8_unit — NOT a fault (confirmed):** with a donor real repeat (real mh 1.0) masking a benign
acceptor (real mh 0.0), `_move_microhomology` returns 1.0, which **EQUALS `max(donor_realbase,
acc_realbase) = 1.0`** (`actual_equals_correct_max = True`). max() faithfully reports the maximum
GENOMIC microhomology across the two shifted boundaries of an ATOMIC move — a correct number. The
"benign acceptor masked" objection is that genomic context is read-blind (can't tell a real compound
move from one-real-one-drift) = the **V2 STRUCTURAL claim, already CLEARED and explicitly out of
scope**. A8 blocks nothing (an acceptor-only move stays independently reachable when it scores best).
→ **A8 is not a second holding fault.**

**geometry/edges:** all off-end guards return 0.0 (`hi+k<=n`, `lo-k>=0`); k=0 no-move → 0.0; k=1
into HP → 1.0; upstream/downstream acceptor shifts handled. `bigk_Nflank` (k=20 N run) → actual 0.75
vs correct 0.0 (still trips ≥0.5) — **A5 scales with drift distance k**, not just k=6.

**EVM truth table: 12/12 PASS** (incl. cap<hold<eff→hold, cap<0→eff, eff==hold→eff no-op,
cap>=eff→eff inert, cap==hold<eff→hold, eff==hold>cap→eff). Invariant `hold<=got<=eff` holds
everywhere. → **`_effective_veto_margin` is CORRECT; CLEAR.**

**Donor leg (helper-level):** donor_N_full → actual 1.0 vs correct 0.0; donor_N_partial → 0.5 vs
0.333; donor_pureACGT → 0.167 (spared); donor_real_repeat → 1.0 (correct). Same phantom on the donor
leg.

**Donor FULL pipeline flip attempt** (`b_donor_pipeline.py` → `b_donor_pipeline.json`): my donor read
construction did NOT produce a scored donor move — the `(D-K,A)` candidate was never emitted even with
the guard OFF (`donor_OFF = D`, `move_emitted_off = false`, `microhom_flagged = 0` in all four cells).
So this harness does not exercise the guard on the donor leg (a read-construction limitation, NOT a
finding — donor-shift scoring anchors differently than my mismatch design assumed). **This does not
weaken A5:** (a) the acceptor-leg pipeline flip (acc_D) is already airtight + isolated; (b) the veto
code at junction_refiner.py:816-818 is **boundary-agnostic** — `if _move_microhomology(...) >=
microhom_threshold: eff_margin += microhom_drift_margin` adds the margin whenever the MAX over legs
trips, with no donor-vs-acceptor distinction — so a donor move that did score better would be vetoed
identically; (c) the donor helper value confirms the same phantom number arises on the donor leg.
A5 is not acceptor-specific.
