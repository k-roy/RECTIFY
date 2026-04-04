# False Junction Handling

Poly(A) tails can create spurious splice junctions (N operations in the CIGAR) when the aligner introduces a skip to align tail bases to a downstream A-tract.

**Implementation:** `rectify/core/false_junction_filter.py`, integrated into `rectify/core/bam_processor.py`

---

## The problem

When a long poly(A) tail aligns beyond the true CPA site and into a downstream A-tract, some aligners insert an `N` (skip) operation to bridge the non-A sequence between two A-runs. This appears as a spurious intron in the alignment.

```
Genome:  ...GCTT|AAAAAAA|CCGG|AAAAAAA|GCATGG...
Read:    ...GCTTAAAAAAAAAAAAAAA
                ↑                ↑
           True CPA          Spurious N operation here

CIGAR:  50M 4N 12M       ← the 4N is not a real intron!
```

---

## The walk-back solution

RECTIFY's walk-back algorithm handles false junctions **automatically** without special detection:

When stepping upstream through the CIGAR, the walk-back eats through:

- All aligned A's
- All deletions (D operations) in the A-region
- **All N operations (skips) in the A-region**

Because false junctions appear within the poly(A) territory (all surrounding bases are A's), the walk-back traverses through them and finds the true CPA site at the first non-A agreement between read and genome.

```python
# In find_polya_boundary():
while pos > start:
    cigar_op = cigar_at(pos)
    ref_base = genome[chrom][pos]

    if cigar_op == 'N':
        # Skip operation — treat as A-territory if flanked by A's
        pos -= 1
        continue

    if ref_base in 'aAtT' or cigar_op == 'D':
        pos -= 1
        continue

    # Non-A reference, non-deletion → check agreement
    if read_base_at(pos) == ref_base:
        return pos   # True CPA site
```

---

## False junction filter module

For cases where false junctions appear further from the poly(A) tail and are not caught by the walk-back, `false_junction_filter.py` provides an explicit filter:

```python
def filter_false_polya_junctions(read, genome, threshold_a_fraction=0.8):
    """
    Remove N operations flanked by high-A-fraction sequence on both sides.

    An N is a false junction if:
    - The 20 bp upstream are > threshold_a_fraction A's
    - The 20 bp downstream are > threshold_a_fraction A's
    """
```

---

## AG mispriming

A related artifact is **AG mispriming**: oligo(dT) primers can prime internally at AG-rich sequences, producing reads that appear to end at non-CPA sites. RECTIFY flags these reads:

```python
def screen_ag_mispriming(read, genome, threshold=0.65):
    """
    Flag reads where the 3' soft-clip region is AG-rich.

    threshold: fraction of A+G bases in the 20 bp upstream of 3' end
    """
```

Flagged reads receive the `AG_RICH` QC flag. They are still reported but should be interpreted with caution for CPA site mapping.

**Implementation:** `rectify/core/ag_mispriming.py`

---

## See also

- [3' End Indel Correction](3prime_indel_correction.md) — the walk-back algorithm
- [Output Formats](../user_guide/output_formats.md) — `qc_flags` column descriptions (`AG_RICH`, `ATRACT_AMBIGUOUS`)
