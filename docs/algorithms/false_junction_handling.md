# False Junction Handling

Poly(A) tails can create spurious splice junctions (N operations in the CIGAR) when the aligner introduces a skip to align tail bases to a downstream A-tract.

**Implementation:** `rectify/core/splice/false_junction_filter.py`, integrated into `rectify/core/bam/bam_processor.py`

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

**Implementation:** `rectify/core/polya/ag_mispriming.py`

---

## Post-consensus N-op junction refinement (Module 2H)

**Implementation:** `rectify/core/splice/junction_refiner.py` — `refine_read_junctions()`

Module 2H runs **post-consensus** on the selected consensus read set. For every N-op in every read, it tests all candidate junctions within a search radius and replaces imprecise N-op boundaries with the best sequence-supported alternative.

### Scoring: sequence-first, canonical/annotated as tie-breakers only

Scoring is **HP-aware edit distance primary**. Canonical GT-AG status and annotated-junction status are tie-breakers only — annotation NEVER overrides a better-scoring junction. The scoring tuple (lower = better) is:

```
(hp_edit_score, canonical_tier, is_alternative, ...)
```

When the current N-op is non-canonical (tier ≥ 4), a canonical-tier alternative receives a 0.5 edit-distance discount (`_CANONICAL_HP_PRIOR = 0.5`) — equivalent to one expected Nanopore HP deletion — so canonical junctions win within the sequencing noise floor.

### `--aligner-bams`: candidate junction pool, not per-aligner correction

`--aligner-bams` provides a **pool of candidate junction coordinates** collected from all aligner BAMs. Module 2H scores every junction in this pool against each read using HP-aware edit distance. The aligner BAMs do not drive per-aligner correction; they expand the candidate set beyond annotated junctions.

### No-candidate-guard policy (PERMANENT)

All junctions in the candidate pool are scored. Non-canonical, non-annotated (novel) alternatives are **not filtered out** before scoring. A novel non-canonical junction can win only if it has a strictly better HP-aware edit-distance score than all canonical and annotated alternatives. Annotation and canonical tier operate as tie-breakers, never as gates. This policy must not be reversed.

### Fast path

Reads already at an annotated canonical-tier-0 junction skip scoring entirely (identical junction coordinates at both boundaries), providing a substantial speedup on typical datasets where most reads are at correct annotated GT-AG junctions.

---

## Cat3 5' junction rescue: full CIGAR surgery is post-consensus

Module 2F (`rescue_3ss_truncation` in `splice_aware_5prime.py`) performs full exon CIGAR surgery **post-consensus**, not during aligner consensus selection. During consensus, only a lightweight 5' soft-clip proximity check is used to score aligners. The full semi-global NW alignment (`local_aligner.py`, affine gap Gotoh 1982) that generates the `five_prime_exon_cigar` column runs after the consensus read is selected.

---

## See also

- [3' End Indel Correction](3prime_indel_correction.md) — the walk-back algorithm
- [Output Formats](../user_guide/output_formats.md) — `qc_flags` column descriptions (`AG_RICH`, `ATRACT_AMBIGUOUS`)
