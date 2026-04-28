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

The ONT cDNA adapter sequence ends in a T-rich, C-containing stub. This regex
matches runs of 1–11 T/C characters anchored at the very end of the window.
Gene-body sequence (which contains G and A) will not match.

---

## Pass 0 — Clean read, no adapter stub

**Trigger:** The regex finds nothing at the end of the window.

**Situation:** Dorado already stripped the adapter completely. The read ends cleanly
in A's — the poly(A) tail runs all the way to the last base.

**What happens:** `_scan_polya()` walks right-to-left from the last base, counting
consecutive A's until it exits the tail into the transcript body.

```
  ...GCUAAGCAAAAAAAAAAAAAAAA
            ^^^^^^^^^^^^^^^^
            poly(A), 16 bp

  Regex:  no match (window does not end in T or C)
  Scan:   16 consecutive A's from the right
  Result: polya_len=16, adapter_seq='', pass=0
```

---

## Pass 1 — Adapter stub cleanly detectable by regex

**Trigger:** The regex matches a `T[CT]+` run at the end of the window.

**Situation:** Dorado left a residual adapter stub after the poly(A). The stub
consists only of T and C characters, so the regex identifies and removes it cleanly.
After stripping the stub, `_scan_polya()` finds the tail boundary in what remains.

```
  ...GCUAAGCAAAAAAAAAAAAAATCTTCT
            ^^^^^^^^^^^^^^~~~~~~
            poly(A), 14 bp  stub

  Regex T[CT]{0,10}$:  matches "TCTTCT" (6 chars) at end
  Strip stub:          window → ...GCUAAGCAAAAAAAAAAAAAA
  Scan:                14 consecutive A's
  Result:              polya_len=14, adapter_seq='TCTTCT', pass=1
```

### Why `T[CT]+` and not just `T+`?

ONT adapter sequences contain both T and C. Allowing C captures real stubs that
a pure `T+` regex would miss. G and A are excluded because they signal gene-body
sequence, not adapter.

---

## Pass 2 — Adapter stub with internal basecalling errors

**Trigger:** After the regex attempt and `_scan_polya()`, `polya_len < min_polya`
**and** the last base of the window is not `'A'`.

**Situation:** A basecalling error inside the adapter stub introduces a non-T/C base,
breaking the regex pattern entirely. For example, a T→G miscall in the stub:

```
  True:    ...GCUAAGCAAAAAAAAAAAAAAATCTTCT
  Called:  ...GCUAAGCAAAAAAAAAAAAAAATCGTCT
                                      ^
                                    T→G miscall in adapter stub
```

The G makes `T[CT]{0,10}$` fail — `TCGTCT` contains G so the regex cannot match.
The raw window then ends in `T`, so `_scan_polya()` finds 0 A's. With
`polya_len=0` and `last_base='T'` (not A), Pass 2 activates.

**What happens:** The algorithm peels 1, 2, 3, … up to 15 bases from the right of
the window, testing each peel depth `k`:

1. `candidate = window[:-k]`  (everything left after peeling)
2. `stub      = window[-k:]`  (the k peeled bases)
3. **Skip** if `candidate[-1] != 'A'` — no poly(A) tail could end here
4. **Skip** if `'T' not in stub` — no adapter T signal; probably gene-body sequence
5. Run `_scan_polya(candidate)` — if ≥ 5 A's found, accept this boundary

Tracing the example above (`polya_len=0`, window ends in `...AAAAAAAAAAAAAAATCGTCT`):

```
  k=1  strip T      last of candidate = C  → skip (not A)
  k=2  strip CT     last of candidate = C  → skip (not A)
  k=3  strip TCT    last of candidate = G  → skip (not A)
  k=4  strip GTCT   last of candidate = C  → skip (not A)
  k=5  strip CGTCT  last of candidate = T  → skip (not A)
  k=6  strip TCGTCT last of candidate = A  ✓
                    stub = 'TCGTCT', 'T' in stub ✓
                    scan candidate: 15 A's ≥ 5 ✓ → ACCEPT

  Result: polya_len=15, adapter_seq='TCGTCT', pass=2
```

### The `'T' not in stub` guard

If the peeled chunk contains no T at all, it is almost certainly gene-body sequence,
not adapter. The guard prevents mis-trimming reads that simply end in a non-A
transcript body base.

Counter-example — a read with no tail that ends in `...CUAGCU`:

```
  k=1  strip U/T   last of candidate = C → skip (not A)
  k=2  strip CU    last of candidate = G → skip (not A)
  k=3  strip GCU   last of candidate = A ✓, but 'T' not in 'GCU' → skip
  k=4  strip AGCU  last of candidate = U → skip (not A)
  ...
  No k produces: last=='A' AND 'T' in stub AND scan ≥ 5 A's
  Result: polya_len=0  →  read passes through unchanged
```

---

## Reads with no poly(A) tail

Any read that exits all three passes with `polya_len=0` is emitted to the output
BAM unchanged.

```
  ...GCUAGCUAGCUAGCUAGCUA   (transcript body, ends in A-containing gene sequence)

  Regex:  no match
  Scan:   0 A's (scan stops at first non-A from the right)
  Pass 2: no k satisfies candidate[-1]=='A' AND 'T' in stub AND polya_len ≥ 5
  Result: polya_len=0  →  read emitted unchanged
```

Typical cases: pre-mRNA intermediates, degraded RNA, non-polyadenylated transcripts.

---

## Summary

| Pass | Trigger condition | Adapter state | Per-read cost |
|------|-------------------|---------------|---------------|
| 0 | Regex: no match; scan finds ≥1 A | None — Dorado already stripped | 1 regex + 1 scan |
| 1 | Regex: matches `T[CT]+$`; scan finds ≥1 A after stripping | Present, cleanly recognisable | 1 regex + 1 scan |
| 2 | polya_len=0 after Pass 1; last base ≠ A | Present but scrambled by basecalling errors | 1 regex + up to 15 scans |

---

## Return value

```python
(polya_len, adapter_seq, last_base, adapter_pass)
```

| Field | Meaning |
|-------|---------|
| `polya_len` | Bases belonging to the poly(A) tail (0 = no tail detected) |
| `adapter_seq` | Stripped adapter stub sequence (empty string for Pass 0) |
| `last_base` | Last base of the trimmed body sequence (sanity: should be non-A) |
| `adapter_pass` | 0, 1, or 2 — which pass found the boundary |

---

## Related documentation

- `docs/DRS_POLYA_ADAPTER_ANALYSIS.md` — empirical characterisation of adapter stubs
  across multiple DRS datasets; motivation for the regex pattern and Pass 2 threshold
- `rectify/core/drs_trim_command.py` — full implementation: `_scan_polya()`,
  `find_polya_and_adapter()`, and the BAM I/O layer
- `rectify/core/restore_polya_command.py` — Step 4: re-attaches trimmed sequence as
  soft-clips for IGV tail-length visualization
