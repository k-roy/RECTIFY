# Soft-Clip Rescue

RECTIFY rescues two types of soft-clipped bases that represent real sequence but are left unaligned by the aligner.

---

## 1. 5' Junction Rescue (Module 2F)

**Implementation:** `rectify/core/splice/splice_aware_5prime.py` — `rescue_3ss_truncation()`

### The problem

Long reads spanning splice junctions often have soft-clipped bases at the 5' end. These bases match the upstream exon, but the aligner couldn't extend the alignment through the junction because it lacked global junction awareness.

```
Exon 1        Intron         Exon 2
──────────────|>>>>>>>>>>>>>|──────────────
         CATGGG|GTAAGC...|AGACGT

Read:  [ATGGG]ACGT...
        ^^^^^
        Soft-clipped (the aligner started at Exon 2)
        These bases match the end of Exon 1!
```

### RECTIFY's solution

Full CIGAR surgery for 5' junction rescue runs **post-consensus** (not during consensus). During consensus, only a lightweight proximity check is performed. The full rescue (`rescue_3ss_truncation`) runs after aligner selection on the consensus read set.

For each read with a 5' soft-clip near a splice junction:

1. Look up the nearest annotated splice junction upstream of the alignment start
2. Fetch the upstream exon sequence (the last N bp, where N = soft-clip length)
3. Align the soft-clipped sequence against the exon sequence using semi-global NW (affine gap, Gotoh 1982) via `local_aligner.py`
4. If the alignment is within edit-distance threshold: extend the alignment to the exon boundary and record the `five_prime_exon_cigar` (v2.8.0)

```python
def rescue_3ss_truncation(read, genome, annotated_junctions):
    """
    Post-consensus 5' soft-clip rescue at 3' splice sites.

    Returns rescue dict including five_prime_exon_cigar if rescued, else rescued=False.
    """
```

### Scoring in consensus selection

During multi-aligner consensus, soft-clipped 5' bases are penalized at **−2 bp per effective clip**. This lightweight score drives aligner selection only; full CIGAR surgery for the rescued exon segment runs post-consensus.

---

## 2. 3' Soft-Clip Rescue at Homopolymer Boundaries

**Implementation:** `rectify/core/correct/indel_corrector.py` — `rescue_softclip_at_homopolymer()`

### The problem

Nanopore basecallers systematically **under-call homopolymer runs**. When the read's T-tract upstream of the CPA site is shorter than the reference T-tract, the aligner places the alignment start at the wrong position and soft-clips the "extra" non-T bases.

```
Genome:  ...TTTTTTTTT|GCATGG...    (9 T's)
Read:    ...TTTTTT   [GCATGG]      (6 T's — under-called; GCA soft-clipped)
                      ^^^^^^
                      These match the reference! The alignment end is wrong.
```

### RECTIFY's solution

1. Identify reads where the aligned region ends inside a T-tract (plus strand: A-tract downstream; minus strand: T-tract upstream)
2. Check if soft-clipped sequence matches the reference at the expected position after the homopolymer
3. If yes: extend the 3' end past the remaining reference homopolymer bases to include the soft-clipped match

```python
def rescue_softclip_at_homopolymer(read, strand, genome, max_rescue_len=10):
    """
    Extend 3' end past reference homopolymer to include matching soft-clipped bases.
    """
```

---

## Effect on 5' end accuracy

!!! note "5' ends and TSS accuracy"
    Due to 5'-to-3' degradation in direct RNA sequencing, the read's 5' end often does **not** represent the true transcription start site (TSS). RECTIFY reports the corrected 5' end but this should be interpreted as the 5'-most position of the captured RNA fragment, not necessarily the TSS.

---

## See also

- [Multi-Aligner Consensus](multi_aligner_consensus.md) — how soft-clip scores feed into aligner selection
- [3' End Indel Correction](3prime_indel_correction.md) — the walk-back algorithm for poly(A) artifacts
