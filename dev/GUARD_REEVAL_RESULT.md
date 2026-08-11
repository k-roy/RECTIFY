# Guard / positional-close re-eval on the corrected residual — VERDICT: SHELVE (2026-07-17)

## Question
After the compensating-indel fix (commit e40ca00) removed ~91% of arm-B fabrication, is the microhomology
guard / positional close (`_positional_signal`, `drift_positional_gate`, `microhom_drift_margin`) still
justified on the corrected residual, or should the dormant machinery be shelved to simplify?

Decision rule (advisor-set): the close earns its keep ONLY if `_positional_signal` SEPARATES the residual
recovery (genuine move → read matches moved acceptor → signal >0) from drift (fabrication → read matches
incumbent → signal ≤0). If the distributions overlap, it does not discriminate on real data → shelve.

## Evidence 1 — residual drift is DIFFUSE
`fixed_4b_dump.tsv` (§4b FIXED, ≥2 samples): 3488 drift, 1462 recovery, 18327 inconclusive. The 3488 drift
junctions are spread across **2368 distinct donor-kb loci on every chromosome** (~1.5/locus). Diffuse,
genome-wide low-level artifact — NOT concentrated real biology. Weak case for a guard on its own.

## Evidence 2 — `_positional_signal` does NOT separate recovery from drift (`pos_signal_eval.py`)
Per-read `_positional_signal(genome, q, q_split, ne=raw_acc, new_je=refined_acc)` over ALL arm-B
acceptor-only moves in the residual set, 13 BAM pairs:

| category | level | n | median | frac>0 | frac≤0 |
|---|---|---|---|---|---|
| **recovery** (genuine) | per-read | 3342 | 0.00 | **48.3%** | 51.7% |
| recovery | per-junction | 188 | 0.00 | 40.4% | 59.6% |
| **drift** (fabrication) | per-read | 72949 | −3.00 | 15.0% | 85.0% |
| drift | per-junction | 1543 | 0.00 | 33.3% | 66.7% |

**Interpretation:**
- DRIFT reads cluster negative (85% ≤0, median −3): the close CAN detect "this read matches its incumbent."
- RECOVERY reads are essentially **NOISE** (48.3% >0 ≈ coin flip, median 0): the close CANNOT detect
  "this read genuinely supports the move." A true recovery should be strongly >0; it is not — ONT error +
  the W=28 window wash out the discriminating signal.
- The close's PURPOSE is to SPARE genuine moves from the drift veto. A `signal>0` spare gate would spare
  only **48% of real recoveries** (losing 52%) while sparing 15% of drift. At the junction level recovery
  (40.4% >0) and drift (33.3% >0) are barely 7 points apart — no usable separation.

**The eval is self-validating (not a broken metric).** The drift channel shows clean, directionally-correct
signal (85% ≤0, median −3) — proving `_positional_signal` is computed correctly. So recovery-at-noise is a
REAL finding, not a computation bug: on real ONT data a genuine acceptor move does not produce a positive
positional signal.

**The close's addressable set is even smaller than the drift total.** `_positional_signal` returns None for
non-acceptor moves, so the gate can only ever fire on the **acceptor-only** subset: ~**1543** drift junctions
(per-junction row above), not 3488. Its ceiling benefit is ~1543 diffuse junctions it cannot cleanly
separate — weaker still.

**Framed against the project's purpose (discovery of real cryptic junctions):** the fix already bought the
specificity win (−91% fabrication) at ZERO recall cost. Enabling the guard would give back ~52% of genuine
recoveries to chase a diffuse residual — it is anti-discovery. That is the decision in one sentence.

**The "recovery" label is imperfect but it doesn't change the verdict:** recovery = SR-exact-support
junction, but a given arm-B read assigned there may not actually splice there (false recovery). Both readings
land on shelve — if recoveries are real, the close fails to protect them; if false, there is nothing to protect.

## Verdict: SHELVE the guard machinery; ship the fix alone
- The compensating-indel fix (always-on invariant, no tuning) does the real work: −91% fabrication,
  precision up ~2× at high confidence, SR-canonical control untouched.
- The residual drift (3488 junctions) is diffuse and small, and the positional close cannot discriminate it
  from genuine recoveries without killing ~half the recoveries. The guard was originally tuned against a
  signal that was 80–90% CIGAR artifact — with that removed, it has no real-data footing.
- The guard params (`microhom_drift_margin`, `drift_near_tie_cap`, `drift_positional_gate`) all DEFAULT to
  0.0 (dormant / byte-identical). Recommendation: **leave them OFF in production** (do not enable/tune),
  and consider removing the dormant machinery in a future cleanup pass for simplicity. No code change
  needed to ship — the fix is the production behavior; the guard is already off by default.

## Files
- Sherlock `sma_recall/`: `fixed_4b_dump.tsv`, `pos_signal_eval.py`, `pos_signal_eval.{tsv,summary.txt}`,
  `recall_result_{OLD,FIXED}_matched.txt`.
- Worktree: `dev/REPLACER_COMPENSATING_INDEL_BUG.md` (fix + §4b), this verdict.
