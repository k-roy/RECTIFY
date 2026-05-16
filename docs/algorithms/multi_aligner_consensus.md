# Multi-Aligner Consensus

RECTIFY runs up to five aligners in parallel and selects the best alignment per read using a composite score.

**Implementation:** `rectify/core/align/multi_aligner.py`, `rectify/core/consensus/consensus.py`

---

## Why multiple aligners?

Different aligners make different tradeoffs:

| Aligner | Strength | Weakness |
|---------|---------|---------|
| minimap2 | Fast; good general splice junction detection | Can be aggressive about junction placement |
| mapPacBio | Conservative at junctions (forces mismatches) | Slow (~10× minimap2) |
| gapmm2 | Handles large indels well | Less accurate at short exons |
| uLTRA | Annotation-guided; best for small exons (11–20 nt) | Requires annotation; slow |
| deSALT | De Bruijn graph; can resolve complex junctions | Memory intensive; known duplication bug |

By running all available aligners and selecting per-read, RECTIFY avoids systematic biases that would arise from committing to one aligner.

---

## Two-phase scheduling

mapPacBio is ~10× slower than the other aligners. RECTIFY uses a two-phase schedule to avoid resource contention:

- **Phase 1:** mapPacBio runs alone with all allocated threads
- **Phase 2:** remaining aligners run in parallel, splitting threads equally

This is faster than running all aligners simultaneously.

!!! warning "deSALT and fork-safety"
    deSALT crashes with "double free or corruption" when forked from a multithreaded Python process. It runs sequentially after the parallel pool in phase 2.

---

## Scoring

Each aligner's output is represented as an `AlignmentInfo` dataclass. The composite score is:

```python
score = (
    - 2 * effective_5prime_clip    # Unrescued 5' soft-clips (-2/bp)
    - 1 * downstream_a_count       # A-tract depth penalty (-1/A, capped at 10)
    - 2 * terminal_errors          # Non-poly(A) terminal errors near 3' end (-2/bp, capped at 10)
    # NOTE: canonical GT-AG and annotated junction counts are tiebreakers only,
    # not additive components of the primary score.
)
```

### 5' soft-clip penalty

The soft-clip penalty uses **sequence-based rescue** (v2.7.5+): soft-clipped bases are aligned against the upstream exon sequence (edit distance ≤ 20% mismatches). Rescued bases carry zero penalty; unrescued bases cost −2/bp.

```python
@dataclass
class AlignmentInfo:
    aligner: str
    five_prime_softclip_seq: str        # Full soft-clipped sequence
    effective_five_prime_clip: int      # After rescue: remaining unrescued length
    effective_three_prime_clip: int     # mapPacBio correction for forced mismatches
    n_junctions: int
    canonical_junctions: int
    annotated_junctions: int
    three_prime_a_count: int            # Downstream A's (from A-tract detector)
```

### mapPacBio effective clip correction

mapPacBio forces mismatches at splice junction boundaries instead of soft-clipping (a deliberate design choice). This makes its 5' junction representation functionally equivalent to a soft-clip but encoded differently in CIGAR. `effective_five_prime_clip` normalises this so mapPacBio is scored fairly against aligners that do soft-clip.

### 3' non-poly(A) terminal error penalty

`effective_three_prime_clip` captures non-A (plus strand) or non-T (minus strand) terminal errors near the 3' end via a greedy sliding-window scan. This penalises aligners that stop before the true 3' end with non-poly(A) content at the tail. Penalty: −2/bp, capped at 10. This is distinct from the A-tract depth penalty, which penalises going too far *into* a downstream genomic A-tract.

---

## Selection logic

```python
def select_best_alignment(read_group, annotated_junctions=None):
    """
    Per-read selection. Consensus operates on RAW (uncorrected) BAMs.
    rectify correct runs ONCE after consensus, on the winning aligner's BAM.

    1. Score each aligner using composite penalty (5' clip, A-tract depth,
       3' non-poly(A) terminal errors) — all on the raw BAM records
    2. Select highest scorer
    3. Tiebreaker 1: canonical GT-AG junction count
    4. Tiebreaker 2: annotated junction count
    5. Tiebreaker 3: majority 3' position vote across aligners
    6. Tiebreaker 4: wider alignment span
    """
```

The winning aligner's BAM record is written to the rectified BAM.

---

## Chimeric consensus (experimental)

With `--chimeric-consensus`, RECTIFY uses a more sophisticated assembly:

1. Find "sync points" where all aligners agree on the alignment
2. Divide the read into segments at sync points
3. Independently select the best aligner per segment
4. Assemble a chimeric CIGAR from the best segments

This can improve accuracy for reads where different aligners are better at different parts, but it is experimental and requires further validation. Off by default.

**Implementation:** `rectify/core/consensus/chimeric_consensus.py`

---

## minimap2 parameters

```bash
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD \
    -t <threads> \
    --junc-bed annotation.junc.bed \
    --junc-bonus 9 \
    genome.fa.gz reads.fastq.gz
```

Key flags:
- `-uf`: forward-strand only (direct RNA / sense-strand cDNA)
- `-k14`: smaller k-mer for noisy nanopore reads
- `-G 5000`: max intron (yeast; increase for other organisms)
- `--splice-flank=no`: disables GT-AG bonus within minimap2 (RECTIFY does its own scoring)

---

## What consensus does NOT do

**Module 2F (Cat3 / `splice_aware_5prime.py`):** During consensus, only a lightweight 5' junction rescue *score* is computed to penalise (or not) each aligner's soft-clip. The full CIGAR surgery — semi-global NW alignment via `local_aligner.py`, exon CIGAR reconstruction — runs entirely **post-consensus** in `bam_processor.py` as Module 2F. Consensus never rewrites CIGARs.

**Module 2H (`junction_refiner.py`):** Post-consensus N-op junction refinement. The `--aligner-bams` option supplies all per-aligner raw BAMs as a **candidate junction pool** only — Module 2H does not re-score or re-correct individual aligners. Scoring is sequence-first (HP-aware edit distance); canonical GT-AG and annotation are tie-breakers, never gates.

**`rectify correct`** runs **once** after consensus on the winning aligner's BAM. It is not part of consensus selection.

---

## See also

- [Aligner Recommendations](https://github.com/k-roy/RECTIFY/blob/main/docs/ALIGNER_RECOMMENDATIONS.md) — guidance on when to add uLTRA or deSALT
- [Soft-Clip Rescue](softclip_rescue.md) — how 5' soft-clips are rescued before scoring
