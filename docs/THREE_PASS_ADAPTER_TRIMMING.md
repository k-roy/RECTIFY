# DRS Poly(A) + Adapter Trimming: Three-Pass Algorithm

**Relevant source:** `rectify/core/drs_trim_command.py` — `find_polya_and_adapter()`

---

## Purpose

Before re-alignment, `rectify trim-polya` strips the poly(A) tail and any residual
adapter stub from each DRS read. The algorithm operates on the **last 150 bases of
the read in RNA 5'→3' orientation** and must locate the exact position where the
transcript body ends and the poly(A) tail begins.

The boundary is non-trivial because:
- Nanopore basecallers miscall A's within the tail (A→T or A→C errors)
- The ONT adapter isn't always fully stripped by Dorado
- Basecalling errors can make the adapter stub unrecognisable to a simple regex

The three passes handle these cases in order of increasing difficulty. Passes 0 and 1
cover >95% of reads cheaply; Pass 2 activates only for the hard cases.

---

## The Adapter Regex

```python
_ADAPTER_RE = re.compile(r'T[CT]{0,10}$')
```

The ONT cDNA adapter sequence ends in a T-rich, C-containing stub.
This regex matches runs of 1–11 T/C characters anchored at the very end of the
window. Gene-body sequence (containing G and A) will not match.

---

## Pass 0 — Clean read, no adapter stub

**Trigger:** The regex finds nothing at the end of the window.

**Situation:** Dorado already stripped the adapter completely. The read ends cleanly
in A's — the poly(A) tail runs all the way to the last base.

**What happens:** `_scan_polya()` runs directly on the raw window, walking rightward
from the transcript body into the tail until it hits a non-A base.

```
Read (last 30 bases, RNA orientation):

  ...G C U A A G C A A A A A A A A A A A A A A A A A A
                    |<-------- poly(A) tail, 16 bp ------->|

Regex:   no match
Scan:    walks left from final A, counts 16 consecutive A's
Result:  polya_len=16, adapter_seq='', pass=0
```

---

## Pass 1 — Adapter stub cleanly detectable by regex

**Trigger:** The regex matches a T[CT]+ run at the end of the window.

**Situation:** Dorado left a residual adapter stub after the poly(A). The stub
consists only of T and C characters, so the regex identifies and removes it cleanly.
After stripping the stub, `_scan_polya()` finds the tail in what remains.

```
Read (last 35 bases, RNA orientation):

  ...G C U A A G C A A A A A A A A A A A A A A T C T T C T
                    |<------ poly(A), 14 bp ------>| |stub|

Regex T[CT]{0,10}$:  matches "TCTTCT" (6 chars) at the end
Strip stub:          window becomes ...GCUAAGCAAAAAAAAAAAAAAAAA
Scan:                finds 14 A's
Result:              polya_len=14, adapter_seq='TCTTCT', pass=1
```

### Why T[CT]+ and not just T+?

ONT adapter sequences contain both T and C. Allowing C in the pattern captures real
adapter stubs that would be missed by a pure `T+` regex. G and A are excluded because
they signal gene-body sequence, not adapter.

---

## Pass 2 — Adapter stub with internal basecalling errors

**Trigger:** After Pass 1 (regex hit or miss) and `_scan_polya()`, the detected
`polya_len` is still below `min_polya` (default 1), **and** the last base of the
remaining window is not `'A'`.

**Situation:** A basecalling error at the poly(A)/adapter boundary scrambles the
stub so the regex either doesn't match or strips the wrong amount. For example, an
A→G miscall inside the stub introduces a G, which breaks `T[CT]{0,10}$`:

```
True sequence:   ...AAAAAAAAAAAAAAAAA | T C T T C T
Called sequence: ...AAAAAAAAAAAAAAAAA | T C G T C T
                                            ^ A→G miscall breaks regex
```

**What happens:** The algorithm peels 1, 2, 3, … up to 15 bases from the right of
the current window. For each peel depth `k`:

1. `candidate = window[:-k]`  (everything except the last k bases)
2. `stub = window[-k:]`       (the k peeled bases)
3. **Skip** if `candidate[-1] != 'A'` — no poly(A) tail could end here
4. **Skip** if `'T' not in stub` — no adapter signal; likely gene-body sequence
5. Run `_scan_polya(candidate)` — if it finds ≥ 5 A's, accept this boundary

```
Read (last 32 bases, RNA orientation):

  ...G C U A A G C A A A A A A A A A A A A A A A T C G T C T
                    |<-------- poly(A), 15 bp ------->|  |stub|

Regex T[CT]{0,10}$:  NO MATCH  (G breaks the pattern)
Scan on raw window:  finds 0 A's (last base is T; scan stops immediately)
polya_len=0, last_base='T', last_base != 'A'  →  Pass 2 activates

k=1: strip 'T'  → candidate ends in 'C'  → skip (not A)
k=2: strip 'CT' → candidate ends in 'C'  → skip (not A)
k=3: strip 'TCT'→ candidate ends in 'G'  → skip (not A)
k=4: strip 'GTCT'→ candidate ends in 'C' → skip (not A)
k=5: strip 'CGTCT'→ candidate ends in 'C'→ skip (not A)
k=6: strip 'TCGTCT'→ candidate ends in 'A' ✓, stub contains 'T' ✓
     scan candidate: finds 15 A's  ≥ 5  ✓
Result:  polya_len=15, adapter_seq='TCGTCT', pass=2
```

### The 'T' not in stub guard

If you peel `k` bases and none of them are T, you probably haven't found an adapter
fragment — you've peeled into gene-body sequence. The guard prevents the algorithm
from misidentifying a read that simply ends in a non-A gene-body base as having a
poly(A) tail.

```
Read ending in gene body (no tail):

  ...G C U A A G C U A G C U A G C U A G C U  (ends in U/T in DNA space)

k=1: strip 'T' → candidate ends in 'C'  → skip (not A)
k=2: strip 'CT'→ candidate ends in 'G'  → skip (not A)
...
No k produces: candidate[-1]=='A' AND 'T' in stub AND polya_len >= 5
Result: polya_len=0  →  read passes through unchanged
```

---

## Summary Table

| Pass | Trigger | Adapter stub state | Cost |
|------|---------|-------------------|------|
| 0 | Regex finds nothing; scan succeeds | None (Dorado already stripped) | 1 scan |
| 1 | Regex matches T[CT]+ at end | Present, recognisable | 1 regex + 1 scan |
| 2 | After regex+scan, polya_len < min_polya and last_base != 'A' | Present but scrambled by basecalling errors | Up to 15 scans |

---

## Return Value

```python
(polya_len, adapter_seq, last_base, adapter_pass)
```

| Field | Meaning |
|-------|---------|
| `polya_len` | Number of bases belonging to the poly(A) tail (0 = no tail detected) |
| `adapter_seq` | The stripped adapter stub sequence (empty string for Pass 0) |
| `last_base` | The last base of the poly(A)-trimmed sequence (sanity check: should be non-A) |
| `adapter_pass` | 0, 1, or 2 — which pass found the boundary |

`polya_len=0` means no tail was detected; the read is emitted to the output BAM
unchanged.

---

## Reads with no poly(A) tail

These pass through all three passes without a match:

```
Read (last 20 bases):  ...G C U A G C U A G C U A G C U A G C U A

Regex:   no match
Scan:    0 A's (scan stops at first non-A from the right)
Pass 2:  no k satisfies both candidate[-1]=='A' and polya_len >= 5
Result:  polya_len=0  →  read emitted unchanged
```

Typical cases: pre-mRNA intermediates, degraded RNA, reads from non-polyadenylated
species.

---

## Related documentation

- `docs/DRS_POLYA_ADAPTER_ANALYSIS.md` — empirical characterisation of adapter stubs
  across multiple DRS datasets; motivation for the regex pattern and Pass 2 threshold
- `rectify/core/drs_trim_command.py` — full implementation including `_scan_polya()`,
  `find_polya_and_adapter()`, and the BAM I/O layer
- `rectify/core/restore_polya_command.py` — Step 4: re-attaches trimmed sequence as
  soft-clips for IGV tail-length visualization
