# MICROHOM DRIFT GUARD — Adversarial Audit V2: DETECTOR / HELPER CORRECTNESS

**Auditor lens:** correctness of the two pure decision functions —
`_move_microhomology(seq, ns, ne, js, je)` and
`_effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap)` —
in `rectify/core/splice/junction_refiner.py`.
**Mode:** READ-ONLY. I may write ONLY this file. Python: `/Users/kevinroy/miniconda3/bin/python`.
**Predecessor:** `dev/MICROHOM_AUDIT_detector-correctness.md` STALLED at PLAN (zero checkpoints). I resume its plan but build every case independently.
**Job:** FIND THE FAULT, not confirm.

## Exact code under audit (transcribed from lines 491-558, verified against source)
```python
def _frac_match(a, b):
    if not a or len(a) != len(b):
        return 0.0
    return sum(1 for x, y in zip(a, b) if x == y) / len(a)

def _move_microhomology(seq, ns, ne, js, je):
    n = len(seq); best = 0.0
    if je != ne:                                   # ACCEPTOR shift (intron END ne->je)
        lo, hi = (ne, je) if ne < je else (je, ne)
        k = hi - lo
        if hi + k <= n:
            best = max(best, _frac_match(seq[lo:hi], seq[hi:hi + k]))
    if js != ns:                                   # DONOR shift (intron START ns->js)
        lo, hi = (ns, js) if ns < js else (js, ns)
        k = hi - lo
        if lo - k >= 0:
            best = max(best, _frac_match(seq[lo - k:lo], seq[lo:hi]))
    return best

def _effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap):
    if drift_near_tie_cap > 0.0 and eff_margin > hold_margin:
        return max(hold_margin, min(eff_margin, drift_near_tie_cap))
    return eff_margin
```
Guard wiring (line 816-818): `if moves and microhom_drift_margin>0: if _move_microhomology(genome, ns, ne, new_js, new_je) >= microhom_threshold: eff_margin += microhom_drift_margin`.
Veto (line 834-839): `if moves and eff_margin>0 and incumbent_score is not None: veto_margin=_effective_veto_margin(...); if best_score_cmp > incumbent_score - veto_margin: moves=False`.
Semantics: `(ns,ne)`=incumbent intron [start,end); `(new_js,new_je)`=candidate. Returns MAX over the two moved boundaries. Guard is strand-SYMMETRIC (no strand arg passed).

## PLAN (case classes — adversarial, independent)
A. `_move_microhomology`:
   A1. Baseline: acceptor downstream drift into tandem repeat -> 1.0; transition -> low. Reproduce.
   A2. DONOR-side drift geometry: upstream repeat -> fires; verify lo-k slice is the correct upstream window.
   A3. Genome-edge OOB: acceptor hi+k>n, donor lo-k<0 -> must return 0 (no crash, no false score).
   A4. k=0 (no move on that boundary) ; k=1 (single-base slide) ; huge k.
   A5. non-ACGT / N / lowercase bases -> _frac_match treats by char equality; N==N scores as match (FALSE-POS risk).
   A6. Overlapping / nested repeats; period-mismatch (period p != k).
   A7. BOTH boundaries shift simultaneously (donor AND acceptor move).
   A8. **max-over-both masking**: a benign transition shift on boundary X, masked by an unrelated
       repeat on boundary Y that ALSO moved -> spurious veto? (the headline concern).
   A9. Strand: is acceptor=downstream / donor=upstream correct for BOTH + and - reads? (guard passes no strand.)
   A10. Self-overlap when k > drift gap: seq[lo:hi] and seq[hi:hi+k] disjoint by construction, but
        seq[lo-k:lo] vs seq[lo:hi] — check the k>span cases and near-edge partial slices.
   A11. FALSE-NEG hunt: a real fab drift the detector MISSES (returns <0.5) -> guard fails to fire.
   A12. FALSE-POS hunt: a real cryptic transition the detector FLAGS (>=0.5) -> guard spuriously vetoes.
B. `_effective_veto_margin`: exhaustive regime grid (hold vs eff vs cap orderings), esp.
   hold>cap (drift move with delta in (cap,hold): still vetoed by hold? must be), cap<=0 no-op,
   eff==hold no-op, negative/nonsensical output hunt, and behavior-at-default (cap=0) byte-identity.
C. Reproduce tests/test_microhom_drift_guard.py; try to BREAK each with a case it misses.

## CHECKPOINTS
(appended below as each sub-step completes — persist every number the moment it lands)
