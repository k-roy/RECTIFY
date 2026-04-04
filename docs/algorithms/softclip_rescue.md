# Soft-Clip Rescue

RECTIFY rescues two types of soft-clipped bases that represent real sequence but are left unaligned by the aligner.

---

## 1. 5' Junction Rescue

**Implementation:** `rectify/core/splice_aware_5prime.py`
**Also in:** `rectify/core/consensus.py` — `_rescue_5prime_softclip()`

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

For each read with a 5' soft-clip:

1. Look up the nearest annotated splice junction upstream of the alignment start
2. Fetch the upstream exon sequence (the last N bp, where N = soft-clip length)
3. Align the soft-clipped sequence against the exon sequence using edit distance
4. If edit distance ≤ 20% mismatches: extend the alignment to the exon boundary

```python
def _rescue_5prime_softclip(read, genome, annotated_junctions):
    """
    Sequence-based 5' soft-clip rescue.

    Returns extended five_prime_position if rescued, else original.
    """
```

### Scoring in consensus selection

During multi-aligner consensus, soft-clipped 5' bases are penalized at **−2 bp per effective clip**. After rescue, the effective clip is the remaining unrescued length. An aligner with fewer unrescued soft-clips scores higher.

---

## 2. 3' Soft-Clip Rescue at Homopolymer Boundaries

**Implementation:** `rectify/core/indel_corrector.py` — `rescue_softclip_at_homopolymer()`

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

## 3. 3'SS Truncation Rescue

**Implementation:** `rectify/core/bam_processor.py` — Module 2F

Some reads are truncated or soft-clipped at the 3' splice site (the exon 2 / 3'SS boundary). RECTIFY identifies these by checking for a canonical AG dinucleotide at the alignment start and attempts to extend the 5' alignment upstream to the annotated 3'SS.

This runs post-consensus on reads with real (non-poly-A-artifact) junctions.

---

## Effect on 5' end accuracy

!!! note "5' ends and TSS accuracy"
    Due to 5'-to-3' degradation in direct RNA sequencing, the read's 5' end often does **not** represent the true transcription start site (TSS). RECTIFY reports the corrected 5' end but this should be interpreted as the 5'-most position of the captured RNA fragment, not necessarily the TSS.

---

## See also

- [Multi-Aligner Consensus](multi_aligner_consensus.md) — how soft-clip scores feed into aligner selection
- [3' End Indel Correction](3prime_indel_correction.md) — the walk-back algorithm for poly(A) artifacts
