# MICROHOM DRIFT GUARD — Adversarial Audit: DETECTOR CORRECTNESS

Auditor lens: **detector correctness** of `_move_microhomology` (rectify/core/splice/junction_refiner.py).
Read-only audit. Python: /Users/kevinroy/miniconda3/bin/python.

## Target under audit (exact code)
```python
def _frac_match(a, b):
    if not a or len(a) != len(b):
        return 0.0
    return sum(1 for x, y in zip(a, b) if x == y) / len(a)

def _move_microhomology(seq, ns, ne, js, je):
    n = len(seq); best = 0.0
    # ACCEPTOR shift: downstream tandem repeat seq[lo:hi] ~ seq[hi:hi+k]
    if je != ne:
        lo, hi = (ne, je) if ne < je else (je, ne)
        k = hi - lo
        if hi + k <= n:
            best = max(best, _frac_match(seq[lo:hi], seq[hi:hi + k]))
    # DONOR shift: upstream tandem repeat seq[lo-k:lo] ~ seq[lo:hi]
    if js != ns:
        lo, hi = (ns, js) if ns < js else (js, ns)
        k = hi - lo
        if lo - k >= 0:
            best = max(best, _frac_match(seq[lo - k:lo], seq[lo:hi]))
    return best
```
Guard: `if moves and microhom_drift_margin>0: if _move_microhomology(genome,ns,ne,new_js,new_je) >= 0.5: eff_margin += margin`.
Semantics: `(ns,ne)` = incumbent intron [start,end); `(new_js,new_je)` = candidate. Acceptor = intron END (`ne->je`);
donor = intron START (`ns->js`). Detector returns MAX over the two moved boundaries; guard is symmetric (no strand).

## PLAN (case classes)
1. Baseline sanity: fab-style downstream-acceptor drift -> detector fires (>=0.5). Reproduce Phase-1 numbers.
2. (a) DONOR-side drift fab: build donor-drift near-tandem repeat, confirm detector fires + guard vetoes.
3. (b) Edge cases: genome ends (lo-k<0 / hi+k>n), k=0, k=1, huge k, N/non-ACGT, overlapping/nested repeats, BOTH shift.
4. (c) max-over-both masking: benign shift on one boundary masked by unrelated repeat on the other.
5. (d) acceptor-vs-donor asymmetry: is the geometry (downstream vs upstream repeat) correct for each direction?
6. Hunt for a case the detector scores WRONG (false-neg = misses real drift; false-pos = flags real discovery).

## CHECKPOINTS
(appended below as each sub-step completes)
