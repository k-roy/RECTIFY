# False Junction Handling

Poly(A) tails can create spurious splice junctions (N operations in the CIGAR) when the aligner introduces a skip to align tail bases to a downstream A-tract.

**Implementation:** `rectify/core/splice/false_junction_filter.py`, integrated into `rectify/core/bam/bam_processor.py`

<figure markdown>
  ![False-junction walk-back schematic: an N-operation (genomic skip) bridges the poly(A) tail to a downstream genomic A-tract; the false-junction filter recognises this pattern (N-op flanked by A's on both sides) and walks back through the N-op as A-territory, recovering the true CPA.](../figures/false_junction_walkback.png){ width="660" }
  <figcaption>When an aligner splices the poly(A) tail to a downstream genomic A-tract, the walk-back treats the N-operation as A-territory (a "false junction") and recovers the true CPA. The companion splice-classification figure (below) shows the full decision tree.</figcaption>
</figure>

<figure markdown>
  ![Splice classification decision tree: every N-operation in the CIGAR is classified as canonical (GT/AG-bounded, both flanks non-A), false-junction (A-flanked on both sides → poly(A) skip), partial-evidence (one A flank), or unsplice (no canonical motif). Each class is handled by a dedicated module.](../figures/splice_classification.png){ width="620" }
  <figcaption>Every N-operation is classified by the dinucleotides at its flanks. False junctions (A-flanked on both sides) are walked through during 3'-end correction; canonical junctions are preserved; partial-evidence junctions are handled by the partial-rescue module.</figcaption>
</figure>

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

For cases where false junctions appear further from the poly(A) tail and are not caught by the walk-back, `false_junction_filter.py` provides an explicit filter. `analyze_junction_for_artifact()` flags an N-op as a poly(A) artifact when the flanking sequence is A-rich (plus strand) / T-rich (minus strand), and `filter_polya_artifact_junctions()` removes the flagged junctions; `correct_3prime_for_artifact_junctions()` then re-derives the 3' end. A-richness is measured by `calculate_a_richness()` / `calculate_t_richness()` over the windows immediately up- and downstream of the junction.

---

## AG mispriming

A related artifact is **AG mispriming**: oligo(dT) primers can prime internally at A/G-rich sequences, producing reads that appear to end at non-CPA sites. RECTIFY screens for this with the weighted A/G score of Roy & Chanfreau 2019 (PMID 31128237):

```python
def screen_ag_mispriming(genome, chrom, position, strand,
                         window=19, threshold=15):
    """
    Score the genomic window immediately downstream (in mRNA orientation) of
    the called 3' end for A/G enrichment, matching the oligo-dT19 primer.
    A weighted A/G score above `threshold` (Youden-optimal cutoff) flags
    likely internal priming.
    """
```

`get_ag_qc_flag()` then labels the read `AG_RICH` (score above threshold),
`PASS` (below), or `INSUFFICIENT_DATA` (window shorter than `min_window`, e.g.
near a chromosome end). `AG_RICH` reads are still reported but should be
interpreted with caution for CPA-site mapping; the read's `confidence` is
downgraded to `low`.

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
- [Output Formats](../user_guide/output_formats.md) — `qc_flags` column descriptions (`AG_RICH`, `PASS`, `INSUFFICIENT_DATA`)
