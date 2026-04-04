# 3' End Indel Correction

The core RECTIFY algorithm corrects systematic alignment artifacts that arise when poly(A) tails land on genomic A-tracts (adenosine-rich sequences downstream of CPA sites).

**Implementation:** `rectify/core/indel_corrector.py` — `find_polya_boundary()`

---

## The problem

### What nanopore aligners do wrong

When a long-read poly(A) tail meets a genomic A-tract, the aligner cannot distinguish tail A's from genomic A's. To maximize the alignment score, it:

1. Extends the aligned region into the A-tract
2. Introduces deletions in the genomic A-tract to "use up" the poly(A) bases
3. Places the 3' end **downstream** of the true CPA site

```
True CPA site
     ↓
─────|AAAAAAA|─────────────────────────────────────────────── Genome
     [TTTTTTTTTTTTTTTTTT]AAAAGGCC                             Read

              ↓ Aligner (wrong):
─────|AAAAAAA|─────────────────────────────────────────────── Genome
     [TTTTTTTTTTTTTTTTTT]AAAAGGCC
 ════════════════════════
         Alignment ends here  ← apparent 3' end (wrong)

              ↓ RECTIFY (correct):
─────|AAAAAAA|─────────────────────────────────────────────── Genome
     [TTTTTTTTTTTTTTTTTT]
 ════════════════
      Alignment ends here ← corrected 3' end (correct)
```

---

## The walk-back algorithm

`find_polya_boundary()` starts at the soft-clip boundary (the raw mapped 3' end) and steps upstream, stopping at the first position where:

1. The reference base is **not A** (plus strand) or **not T** (minus strand), AND
2. The read base at that position **agrees with the reference**

At each step, it skips:
- **A's and T's** — poly(A) territory
- **Deletions** in the CIGAR — inserted by the aligner to span the A-tract
- **T errors** — common sequencing errors in poly(A) runs (T→ miscall)

```python
def find_polya_boundary(read, strand, genome, min_polya_len=5):
    """
    Walk upstream from 3' end to find true CPA site.

    Returns dict with:
      corrected_pos: int     — new 3' end position (0-based)
      original_pos:  int     — raw aligned 3' end
      polya_aligned: int     — A's within alignment
      correction_bp: int     — shift applied (negative = upstream)
    """
```

### Plus strand example

```
Genome: ...GCTT|AAAAAAA|GCATGG...     (| marks A-tract)
Read:   ...GCTTAAAAAAAAGC              (poly(A) tail extends into A-tract)

Walk from right:
  pos 20: G — reference G, read A → mismatch, skip
  pos 19: A — reference A, skip
  pos 18: A — reference A, skip
  ...
  pos 14: T — reference T, read T → AGREE, non-A → STOP

Corrected 3' end: position 14
```

### Minus strand

For minus strand genes, the poly(A) tail appears as a **poly(T) prefix** on the left end of the alignment. The walk-back runs right-to-left from `reference_start`, looking for the first non-T position where read and genome agree.

---

## Soft-clip rescue at homopolymer boundaries

A related correction handles basecaller under-calling at homopolymer runs. When the basecaller miscounts T's in an upstream T-tract, the aligner soft-clips the uncounted bases instead of placing them. RECTIFY rescues these:

```
Genome: ...TTTTT|GCATGG...
Read:   ...TTT[GCA...]         ← 2 T's under-called; GCA soft-clipped

RECTIFY: skip remaining reference T's, match soft-clip GCA → extend alignment
```

**Implementation:** `rescue_softclip_at_homopolymer()`

---

## A-tract ambiguity detection

After walk-back, RECTIFY classifies each read:

| A-count downstream | Classification | Action |
|--------------------|---------------|--------|
| 0 | **Clear** — no A-tract | No correction needed |
| 1–5 | **Minor ambiguity** | Correction applied, `MEDIUM` confidence |
| 6–15 | **Moderate ambiguity** | Correction applied, `LOW` confidence |
| >15 | **Severe ambiguity** | NET-seq refinement attempted; `ATRACT_AMBIGUOUS` flag |

**Implementation:** `rectify/core/atract_detector.py`

---

## Statistics

In typical *S. cerevisiae* direct RNA experiments:

- ~52% of reads require no correction (position already correct)
- ~32% shift 1–5 bp upstream
- ~13% shift 6–15 bp upstream
- ~3% shift > 15 bp (high A-tract depth)
- Mean shift: −3 to −7 bp

---

## See also

- [NET-seq Refinement](netseq_refinement.md) — resolves remaining A-tract ambiguity
- [False Junction Handling](false_junction_handling.md) — removes spurious junctions created by poly(A) aligning to downstream A-tracts
- [Coordinate System](../coordinate_system.md) — 0-based coordinate convention used throughout
