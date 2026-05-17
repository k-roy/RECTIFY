# 3' End Indel Correction

The core RECTIFY algorithm corrects systematic alignment artifacts that arise when poly(A) tails land on genomic A-tracts (adenosine-rich sequences downstream of CPA sites).

**Implementation:** `rectify/core/correct/indel_corrector.py` — `find_polya_boundary()`

<figure markdown>
  ![Three-row schematic showing genome, aligned read, and corrected 3' end. The aligned read extends past the true CPA into a downstream A-tract with deletions; the walk-back arrow steps back from the soft-clip boundary, skipping A's and a sequencing-error T, and lands on the CPA position one base 5' of the A-tract.](../figures/indel_correction.png){ width="680" }
  <figcaption>The walk-back algorithm steps inward from the alignment's 3' boundary, skipping A's, deletions, and T sequencing errors, and stops at the first position where the read base matches the (non-A) reference base — the true CPA site.</figcaption>
</figure>

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

<figure markdown>
  ![Two-panel comparison contrasting reference-only A-tract detection (which stops at the first genomic A and over-corrects) with the read-vs-reference walkback (which walks past every A regardless of the genome and stops at the first non-A position where read and reference agree).](../figures/walkback_readvsref.png){ width="680" }
  <figcaption>The walk-back compares <em>read</em> and <em>reference</em> base by base, not reference alone — so genomic A's on the way to the CPA are walked past correctly. The terminal-gate is a single rule: "terminal base is not A and matches the genome → already correctly anchored," with no exception for genomic-A-matching terminal A's.</figcaption>
</figure>

At each step, it skips:
- **A's and T's** — poly(A) territory
- **Deletions** in the CIGAR — inserted by the aligner to span the A-tract
- **N operations (skips)** — treated as A-territory when flanked by A's (handles false junctions in poly-A regions)
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

### Guard: large-deletion pre-scan (v2.9.3)

Before the main walk-back, `find_polya_boundary` scans for large deletions (≥ 5 bp) within 50 bp of the RNA 3' end when the alignment starts in a poly-T/poly-A context. Minimap2 sometimes bridges a mismatching region with a large deletion when extending a poly-A tail alignment, causing false-positive non-A/non-T matches before the deletion that stop the walk-back too early (under-correction by tens of bp). The pre-scan detects and skips these artifactual large-deletion spans.

### Guard: N-op boundary for spliced minus-strand reads (v2.9.4)

For minus-strand reads, `find_polya_boundary` records the first N-op (intron skip) boundary encountered while parsing the CIGAR and limits the forward scan to positions before that boundary. Without this guard, the scan could silently cross a spliced intron and find a spurious non-T match in a downstream exon, producing an erroneous correction spanning the entire intron. If no non-T match is found before the N boundary, `first_n_start` is returned as the CPA (the exon-intron boundary is the natural terminal-exon end).

### Guard: trailing-base false-stop (v3.0.3)

Before accepting a candidate stop position (read base matches reference base at a non-A/non-T position), `find_polya_boundary` inspects the K = 4 positions immediately toward the poly-A tail. If all K positions have `rb = 'A'` (plus strand) or `rb = 'T'` (minus strand) AND at least one has a mismatching genomic base, the candidate is rejected and scanning continues. This prevents the trailing base of a poly-A tail (e.g., a T that coincidentally matches a genomic T at the alignment boundary) from causing premature termination inside the poly-A zone.

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

## Soft-clip rescue at homopolymer boundaries (Module 2G)

A separate module (`rescue_softclip_at_homopolymer`, Module 2G) handles basecaller under-calling at homopolymer runs. When the basecaller miscounts T's in an upstream T-tract, the aligner soft-clips the uncounted bases instead of placing them. RECTIFY rescues these by extending the 3' end outward past the remaining reference homopolymer bases to include the soft-clipped match.

```
Genome: ...TTTTT|GCATGG...
Read:   ...TTT[GCA...]         ← 2 T's under-called; GCA soft-clipped

RECTIFY: skip remaining reference T's, match soft-clip GCA → extend alignment
```

Module 2G runs for **all protocols** (DRS and dT-primed cDNA). It runs after the poly-A walk-back in the correction order and takes priority over opposite-direction corrections to prevent cancellation.

**Implementation:** `rectify/core/correct/indel_corrector.py` — `rescue_softclip_at_homopolymer()`

---

## A-tract ambiguity detection

After walk-back, RECTIFY classifies each read:

| A-count downstream | Classification | Action |
|--------------------|---------------|--------|
| 0 | **Clear** — no A-tract | No correction needed |
| 1–5 | **Minor ambiguity** | Correction applied, `MEDIUM` confidence |
| 6–15 | **Moderate ambiguity** | Correction applied, `LOW` confidence |
| >15 | **Severe ambiguity** | NET-seq refinement attempted; `ATRACT_AMBIGUOUS` flag |

**Implementation:** `rectify/core/polya/atract_detector.py`

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
