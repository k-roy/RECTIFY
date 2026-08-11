# MICROHOM AUDIT V4 — Detector Correctness (Auditor A)

**Task:** detector-correctness | **Auditor:** A | **Commit:** 05664bc
**Target helpers:** `_move_microhomology(seq, ns, ne, js, je)` + `_effective_veto_margin(hold, eff, cap)`
in `rectify/core/splice/junction_refiner.py`.
**Lens:** detector/helper CORRECTNESS — inputs that score WRONG. Prioritize A5 (`N==N` false-positive
microhomology) + A8 (max-over-both-boundaries masking).

Status: LIVE. Independent of Auditor B. READ-ONLY on repo.

---

## CHECKPOINT 0 — mechanics verified (from source)

`_frac_match(a, b)`: returns 0.0 if empty or unequal length; else `sum(x==y) / len(a)`. **Pure byte equality
`x == y`** — so `N==N`, `n==n`, `X==X`, any char equal to itself, counts as a MATCH.

`_move_microhomology(seq, ns, ne, js, je)`:
- `best = 0.0`
- Acceptor shift (`je != ne`): `lo,hi = sorted(ne, je)`; `k = hi-lo`; if `hi+k <= n`: `best = max(best, _frac_match(seq[lo:hi], seq[hi:hi+k]))`.
- Donor shift (`js != ns`): `lo,hi = sorted(ns, js)`; `k = hi-lo`; if `lo-k >= 0`: `best = max(best, _frac_match(seq[lo-k:lo], seq[lo:hi]))`.
- returns `best` (max over the two boundaries).

Caller (line ~816): if `microhom_drift_margin>0` AND `_move_microhomology(...) >= microhom_threshold(0.5)`:
`eff_margin += microhom_drift_margin`; then veto (line ~834): `if best_score_cmp > incumbent_score - veto_margin: moves=False`.
**Guard firing RAISES the veto bar → HOLDS the move at incumbent.**
- False-POSITIVE frac (spuriously >= 0.5): vetoes a move that should be allowed → **discovery loss on a REAL move**.
- False-NEGATIVE (spuriously < 0.5): fails to guard a fabricated drift → guard silently ineffective.

`_effective_veto_margin(hold, eff, cap)`: `if cap>0 and eff>hold: return max(hold, min(eff, cap)); else return eff`.

**A5 hypothesis:** genomic `N`-run (ambiguity / masked repeat) makes `_frac_match` of two all-N k-mers = 1.0,
tripping guard on ANY move landing boundaries in the N region — INCLUDING a real transition move. FALSE POSITIVE.

**A8 hypothesis:** `max` over both boundaries. A move that shifts BOTH donor and acceptor: acceptor side is a
genuine transition (frac ~0.25, should be spared) but donor side has an unrelated repeat (frac >=0.5) → max
trips guard → real acceptor cryptic vetoed. Masking = the benign boundary's low frac hidden by the other's repeat.

---

## CHECKPOINT 1 — harness written
(pending)

## CHECKPOINT 1 — harness written + FIRST NUMBERS (all probes ran, 0 unexpected)

Harness: `/private/tmp/.../audit_v4/detector-correctness_A/harness.py`. Direct function calls.

**A5 CONFIRMED (unit level):** `FM('NNNN','NNNN')=1.0`, `FM('nnnn','nnnn')=1.0`, `FM('XXXX','XXXX')=1.0`.
Constructed masked-acceptor move: `MH=1.0` (>=0.5 -> VETO fires). Identical geometry with the REAL distinct
underlying sequence (unmasked): `MH=0.0` (spared). So an N-masked region FLIPS a real transition move from
spared -> vetoed. `_hp_run_across` ALSO counts an N-run (HP('NNNNN',2,4)=5) — so hp_drift guard shares the flaw.

**A8 CONFIRMED (unit level):** move shifting BOTH donor(+4) and acceptor(+4): donor-side frac=1.0 (upstream
repeat), acceptor-side frac=0.0 (genuine transition). Combined `max` MH=1.0 -> VETO fires. Acceptor-only
(donor held) MH=0.0 -> would be spared. So an UNRELATED donor repeat masks a genuine acceptor cryptic.

**EVM: all 12 regime probes PASS** — cap<=0 / eff==hold no-ops byte-identical; active returns
max(hold,min(eff,cap)); veto never below hold; monotone in cap; range invariant [hold,eff]. **EVM is CORRECT.**

**Geometry/edges/k:** direction-symmetric; edges (lo-k<0, hi+k>n) safely skip (MH=0, no IndexError);
k=0 no-move=0; k=1 into HP=1.0; huge-k-fits=1.0. All as expected — no crash, no wrong slice.

### OPEN: realism gate (decides HOLDING vs non-holding)
A5/A8 are unit-level false positives. Must verify they REACH the guard in the real pipeline:
- A5: does yeast/human genome_seq contain N at positions where the refiner proposes a move? Does refiner
  even enumerate candidate junctions whose boundary lands in an N-run?
- A8: does the refiner ever propose a SIMULTANEOUS donor+acceptor shift (both js!=ns AND je!=ne)?
- Both: guard is OFF by default (microhom_drift_margin=0.0). Is it ever enabled in shipped configs/presets?

## CHECKPOINT 2 — REALISM GATE (reachability in the real pipeline)

Re-ran prior unit harness: **0 FAILS** (A5, A8, EVM, geometry, edges, k, repeats all as prior). Confirmed stable.

**R1 — Guard is UNWIRED in every shipped command path.** `microhom_drift_margin`, `hp_drift_margin`,
`drift_near_tie_cap` appear ONLY in `junction_refiner.py` + `tests/test_microhom_drift_guard.py`.
- `grep` across `rectify/core/commands/`, `rectify/core/config.py`, `rectify/slurm_profiles/`: **0 hits.**
- `correct_command.py:746 refine_bam_junctions(...)` passes NO drift kwargs → all default 0.0 → guard NEVER fires
  in `rectify correct`. No CLI flag, no config key, no preset enables it. **Guard reachable only via direct
  Python API / the dedicated test.** So A5/A8 CANNOT manifest through any user-facing command today.

**R2 — genome_seq PRESERVES N.** `load_genome` (rectify/utils/genome.py:152,155) = `str(record.seq).upper()`
— uppercases but does NOT strip/replace N. So genomic N (assembly gaps / masking) survives into `genome_seq`
that the guard reads. => A5's N-run precondition is satisfiable IF a candidate boundary lands in an N-run.

**R3 — simultaneous donor+acceptor shift is enumerable.** `_candidates_near` (junction_scoring.py:754-759)
emits arbitrary `(js,je)` pairs from the junction index; the move gate sets `new_js,new_je` from the best
tuple with NO constraint that only one boundary differ. So `new_js!=ns AND new_je!=ne` is reachable => A8's
both-boundary geometry is reachable IF guard enabled.

### Verdict axis so far
- Helper CORRECTNESS: A5 + A8 are GENUINE detector defects (unit-confirmed false positives). EVM CORRECT.
- HOLDING-ness hinges on the guard being OFF-by-default AND UNWIRED → today zero user impact.
- REMAINING to nail: (a) does bundled yeast genome actually contain N at plausible boundaries? (b) exact args
  to _move_microhomology (ns,ne vs new_js,new_je) — any incumbent/target arg-swap bug?

## CHECKPOINT 3 — N-in-genome measured; A5 realism narrowed

**Bundled YEAST S288C (R64-5-1): total=12,157,105 bp, N=0, n_runs=0.** The default `--Scer` genome has NO N
at all → A5 can NEVER fire on yeast. (Finished/telomere-to-telomere-quality assembly.)
**No bundled HUMAN FASTA** (only penalty tables under data/genomes/homo_sapiens/) → human genome is
user-supplied at runtime. Real GRCh38 primary assembly does contain ~150 Mb N (gaps/masking).

**A5 realism chain (to be a HOLDING fault) requires ALL of:**
1. guard enabled (microhom_drift_margin>0) — NOT wired in any command (R1) → fails immediately today.
2. genome_seq contains an N-run — false on bundled yeast; true only on user human/other with gaps.
3. a *candidate junction boundary* (js,je from the aligner-observed pool) lands INSIDE that N-run — but
   reads don't align across assembly-gap N (no bases to match), so the observed junction pool essentially
   never places a boundary in a gap-N-run. Soft-masked lowercase is UPPERCASED to ACGT by load_genome
   (`.upper()`), NOT to N — so masking does NOT create N. The only N reaching the guard = hard-masked/gap N,
   which is exactly where reads don't map.
=> A5 is a REAL detector defect (N==N mis-scored as match) but its pipeline-realization path is effectively
   empty: default genome has no N, and the N that exists elsewhere sits where no candidate boundary lands.

**A8 realism:** both-boundary shift IS enumerable (R3). The masking scenario needs one boundary to be a
genuine transition AND the other to sit in a >=0.5 microhomology repeat AND the best-scoring candidate to be
that simultaneous shift AND guard enabled. Reachable in principle on real data IF guard on; unit-confirmed.

## CHECKPOINT 4 — VERDICT CORRECTION: HOLD (A5 + A8 are real detector faults)

**Task decision rule (verbatim): "is either helper INCORRECT? If A5/A8 are real, they are HOLDING faults."**
Both are real. => VERDICT = **HOLD**. The off-by-default/unwired facts are the DISCOVERY-LOSS auditor's
scope (structural leg, explicitly CLEARED / "do not re-litigate" in my brief); they are severity/scope
context here, NOT grounds to clear a confirmed detector-correctness fault. My prior "HOLDING-ness hinges on
unwired -> zero impact" framing was WRONG and is retracted.

### A8 — CONFIRMED HOLDING, needs NO N (the load-bearing correction)
Pure-ACGT, k=3 (well within default `junction_max_boundary_shift=50`), both-boundary shift:
- donor upstream `(CAG)` tandem repeat -> donor frac = **1.0**
- acceptor = genuine sequence transition (a REAL non-canonical acceptor the read distinguishes) -> frac = **0.0**
- combined `max` MH = **1.0** >= 0.5 -> **VETO fires** (move held at incumbent = discovery loss)
- acceptor-only (donor held) MH = 0.0 -> would be SPARED.
=> an UNRELATED donor `(CAG)n` microsatellite (ubiquitous at real splice boundaries) masks a genuine acceptor
   cryptic via the max-over-boundaries. **Violates the detector's own docstring contract** ("A move to a
   genuine sequence transition ... is untouched, preserving motif-blind discovery"). No N required, no
   yeast/human distinction — fires on the DEFAULT --Scer genome. The yeast-N=0 result clears at most A5-via-N;
   it does NOTHING for A8.
Also: single-boundary acceptor into a 2bp coincidental repeat (`AT|AT`, random freq ~1/16) -> MH=1.0 -> VETO.
So even without both-boundary geometry, short coincidental direct repeats over-trip the 0.5 threshold.

### A5 — real detector fault; pipeline-realization narrow but non-empty
`_frac_match` scores N==N (any x==x incl. lowercase/X) as a match -> all-N k-mers = 1.0 microhomology.
`_hp_run_across` shares it (`HP('NNNNN',2,4)=5`). Realization requires genome_seq N at a candidate boundary:
default yeast has 0 N; user human GRCh38 has ~150 Mb N but at gaps where reads don't map. So A5 is a genuine
SCORING defect with a near-empty realistic trigger — reported as a fault with narrow scope, secondary to A8.

### near-tie cap does NOT rescue either
The cap is on the SAME score axis as delta_improve; a FALSE-POSITIVE MH still adds the drift margin, so
`_effective_veto_margin(0.5, 1.5, cap)` returns a nonzero veto_margin for every cap (1.5/0.5/1.0/1.5 at
cap=0/0.25/1.0/5.0). The cap can only let a move through if the read's delta_improve>=cap; it cannot correct
the wrong frac SCORE. So the 05664bc cap does not fix detector correctness.

### EVM — CLEAR (within the HOLD)
All 12 regime probes pass (cap<=0 / eff==hold no-ops byte-identical; active = max(hold,min(eff,cap)); veto
never below hold; monotone in cap; range [hold,eff]). `_effective_veto_margin` is CORRECT. One broken helper
(_move_microhomology / _frac_match) does not contaminate the other.

### FIX (recommendation)
1. `_frac_match`: treat any non-ACGT (N/n/X/etc.) as a NON-match (`x==y and x in 'ACGT'`), killing the N==N
   false positive in BOTH `_move_microhomology` and (analogously) guarding `_hp_run_across` against N-runs.
2. Evaluate the two boundaries INDEPENDENTLY (guard each moved boundary on its own frac) instead of `max`
   over both — so an unrelated repeat on the held/other boundary can't mask a genuine transition (fixes A8).
3. Consider a minimum k (>=3-4) and/or requiring the repeat to be the boundary that ACTUALLY moved, to stop
   2bp coincidental direct repeats from tripping 0.5.
4. Scope note (for the discovery-loss auditor, not a clearance here): guard is currently unwired in every
   shipped command (correct_command passes no drift kwargs; 0 CLI/config/profile hits), so NO current
   user-facing impact — but the feature is flag-gated to be enabled, and detector correctness is intrinsic.

**VERDICT: HOLD** — `_move_microhomology`/`_frac_match` are INCORRECT (A8 unconditionally on default genome;
A5 on N-containing genomes). `_effective_veto_margin` is CORRECT.
