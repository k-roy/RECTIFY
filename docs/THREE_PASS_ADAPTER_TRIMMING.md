# DRS Poly(A) + Adapter Trimming: Three-Pass Algorithm

**Relevant source:** `rectify/core/commands/drs_trim_command.py` — `find_polya_and_adapter()`

<figure markdown>
  ![Three-pass pre-trim schematic: Pass 1 scans the last 150 bases for a high-confidence adapter motif and trims poly(A) + adapter; Pass 2 handles reads where the adapter is degraded or absent by anchoring on the poly(A) run; Pass 3 catches the remaining edge cases by looking for an A-rich window without requiring an adapter.](figures/polya_pretrim.png){ width="680" }
  <figcaption>The three-pass pre-trim ensures that the poly(A) tail and any residual adapter are removed before alignment, so downstream walk-back operates on a read whose 3' boundary is genomic sequence rather than artefactual tail.</figcaption>
</figure>

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
cover >95% of reads cheaply; Pass 2 activates only for the hard cases. Between the
regex step and the poly(A) scan, the algorithm also peels a terminal triplet-repeat
basecaller artifact when present (see "Terminal repeat strip" below).

---

## The Adapter Regex

```python
_ADAPTER_RE = re.compile(r'[TU][CTU]{0,10}$')
```

The ONT SQK-RNA004 DNA adapter (`GGCTTCTTCTT`) leaves a T/C-rich residual stub
after Dorado pre-trimming. This regex matches a run of 1–11 T/C characters
anchored at the very end of the window, accepting `U` as well as `T` so it works
on RNA-alphabet reads. Gene-body sequence (which contains G and A) will not
match. See `docs/DRS_POLYA_ADAPTER_ANALYSIS.md` for the empirical characterisation
of this stub.

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

**Trigger:** The regex matches a `[TU][CTU]+` run at the end of the window.

**Situation:** Dorado left a residual adapter stub after the poly(A). The stub
consists only of T and C characters, so the regex identifies and removes it cleanly.
After stripping the stub, `_scan_polya()` finds the tail boundary in what remains.

```
  ...GCUAAGCAAAAAAAAAAAAAATCTTCT
            ^^^^^^^^^^^^^^~~~~~~
            poly(A), 14 bp  stub

  Regex [TU][CTU]{0,10}$:  matches "TCTTCT" (6 chars) at end
  Strip stub:             window → ...GCUAAGCAAAAAAAAAAAAAA
  Scan:                   14 consecutive A's
  Result:                 polya_len=14, adapter_seq='TCTTCT', pass=1
```

### Why `[TU][CTU]+` and not just `T+`?

ONT adapter sequences contain both T and C. Allowing C captures real stubs that
a pure `T+` regex would miss. G and A are excluded because they signal gene-body
sequence, not adapter.

---

## Terminal repeat strip — `(AAG/GAA)n` basecaller artifact

After the regex step and **before** the poly(A) scan, `find_polya_and_adapter()`
peels a terminal triplet-repeat block (`_find_terminal_repeat_block`) when one is
present. Some reads carry a periodic `(AAG)n` / `(GAA)n` artifact 3' of the genuine
poly(A) tail. Because this block has no terminal A-run, `_scan_polya()` cannot see
past it; if left in place, the artifact would be force-aligned and the true tail
boundary lost. Peeling it first exposes any genuine poly(A) run that lies 5' of the
artifact for the scan.

The block is peeled only when it is at least `_REPEAT_MIN_LEN` = 15 bp (≥ 5 tandem
copies of a triplet) with canonical-motif purity ≥ `_REPEAT_MIN_FRAC` = 0.8. The
peeled length and motif are returned as `repeat_len` and `repeat_motif`; callers
fold the length into the total trim: `total_trim = repeat_len + polya_len + len(adapter_seq)`.
Stripping is controlled by the `strip_repeat` argument (on by default).

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

The G makes `[TU][CTU]{0,10}$` fail — `TCGTCT` contains G so the regex cannot match.
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
| 1 | Regex: matches `[TU][CTU]+$`; scan finds ≥1 A after stripping | Present, cleanly recognisable | 1 regex + 1 scan |
| 2 | polya_len=0 after Pass 1; last base ≠ A | Present but scrambled by basecalling errors | 1 regex + up to 15 scans |

(A terminal `(AAG/GAA)n` repeat block, when present, is peeled between the regex
and the scan in every pass — see "Terminal repeat strip" above.)

---

## Return value

```python
(polya_len, adapter_seq, last_base, adapter_pass, repeat_len, repeat_motif)
```

| Field | Meaning |
|-------|---------|
| `polya_len` | Bases belonging to the poly(A) tail (0 = no tail detected), measured on the window *after* the repeat block is peeled |
| `adapter_seq` | Stripped adapter stub sequence (empty string for Pass 0) |
| `last_base` | Last base of the trimmed body sequence (sanity: should be non-A) |
| `adapter_pass` | 0, 1, or 2 — which pass found the boundary |
| `repeat_len` | Length of the peeled terminal `(AAG/GAA)n` artifact block (0 if none / strip disabled); sits 3' of `polya_len` |
| `repeat_motif` | Canonical motif of the peeled block (`''` if none) |

---

## Related documentation

- `docs/DRS_POLYA_ADAPTER_ANALYSIS.md` — empirical characterisation of adapter stubs
  across multiple DRS datasets; motivation for the regex pattern and Pass 2 threshold
- `rectify/core/commands/drs_trim_command.py` — full implementation: `_scan_polya()`,
  `find_polya_and_adapter()`, and the BAM I/O layer
- `rectify/core/commands/restore_polya_command.py` — Step 4: re-attaches trimmed sequence as
  soft-clips for IGV tail-length visualization
