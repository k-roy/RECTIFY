# MICROHOM AUDIT V5 — signal-correctness (Auditor A)

**Task:** Adversarially test `_semiglobal_ed` + `_positional_signal` correctness in
`rectify/core/splice/junction_refiner.py` (worktree branch `worktree-agent-a25a2c1e784ad37dc`,
HEAD a97ff6d, includes CLOSE commit cd3de46).

**Mode:** READ-ONLY. Independent of Auditor B. Verdict CLEAR/HOLD with reproducing cases.

**Python:** `/Users/kevinroy/miniconda3/bin/python`
**Import verified from worktree:** `jr.__file__ = /Users/.../worktrees/agent-a25a2c1e784ad37dc/rectify/core/splice/junction_refiner.py`

---

## CHECKPOINT 0 — orientation + geometry (DONE)

### The two functions under audit (verbatim)

```python
def _semiglobal_ed(query, ref) -> int:   # lines 580-597
    m, n = len(query), len(ref)
    if m == 0: return 0
    if n == 0: return m
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        cur = [i] + [0] * n
        qi = query[i - 1]
        for j in range(1, n + 1):
            cost = 0 if qi == ref[j - 1] else 1
            cur[j] = min(prev[j-1]+cost, prev[j]+1, cur[j-1]+1)
        prev = cur
    return min(prev)   # min over LAST ROW => free ref SUFFIX; row0 = range() => free ref PREFIX too?!

def _positional_signal(genome_seq, q, q_split, ne, new_je, W=28) -> Optional[int]:   # 600-629
    if new_je == ne: return None
    rescue = q[q_split:q_split + W]
    if not rescue: return None
    n = len(genome_seq)
    ref_inc = genome_seq[ne:min(n, ne + W + 6)]
    ref_mov = genome_seq[new_je:min(n, new_je + W + 6)]
    return _semiglobal_ed(rescue, ref_inc) - _semiglobal_ed(rescue, ref_mov)
```

### Geometry established (verify not trust):
- `q = read.query_sequence` (line 703) — pysam returns the read **reference-oriented**
  (revcomp of original for is_reverse). SAME coordinate frame as `genome_seq` for BOTH strands.
- `q_split` from `_iter_n_ops` (366-401): query bases consumed before the N-op, walking CIGAR
  **left-to-right in reference-coordinate order**. So `q[q_split:]` = read's exon2 in GENOME coords.
- `ne` = incumbent intron_end (genome); `new_je` = candidate acceptor intron_end (genome).
  `genome_seq[ne:]` / `genome_seq[new_je:]` = incumbent / moved reference exon2. Coordinate-consistent
  with `q[q_split:]` on both strands (all three in reference frame). **=> minus-strand geometry
  appears CORRECT a priori — must confirm empirically.**
- Wiring (915-930): drift-flagged would-be-veto SPARED iff `psig is not None and psig >= gate`
  (gate default 0.0 = OFF). Acceptor moves only (`new_je==ne` → None → not spared → veto stands).
- The "free-k soft-clip escape": `_score_junction` (junction_scoring.py:983) loops `k in [0,L)`,
  scores `rescue[k:]` left-anchored + `rescue[:k]` right-anchored; picking a clipping k lets the
  INCUMBENT also score ~0 → delta_improve ~0. `_semiglobal_ed` hard-anchors query[0]@ref[0] → no
  clip escape → exposes discriminating exon2 bases. THIS is the design claim to verify.

### FIRST SUSPICION (to test): docstring says `_semiglobal_ed` is "HARD-anchored at ref[0]
— no free query/ref prefix". But `prev = list(range(n+1))` initializes row 0 to [0,1,2,...,n],
which is the FREE-PREFIX initialization (aligning empty query to any ref prefix costs 0... no wait,
range gives cost j for consuming j ref chars). Need to trace: is the ref PREFIX actually free or
not? If `cur[j-1]+1` lets the query slide right (skip ref prefix) cheaply via the first row, the
"hard-anchored" claim could be FALSE. Build a case: query="XYZ", ref="AAAXYZ" — if hard-anchored,
ed should be high (3 mism at ref[0:3]); if free-prefix, ed=0. TEST THIS FIRST.

---
(checkpoints appended below as work proceeds)
