# Multi-Aligner Consensus

RECTIFY runs up to five aligners in parallel and selects the best alignment per read using a composite score.

**Implementation:** `rectify/core/multi_aligner.py`, `rectify/core/consensus.py`

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
    - 2 * effective_5prime_clip    # Unrescued 5' soft-clips
    - 1 * downstream_a_count       # A-tract depth penalty (capped at 10)
    - 2 * terminal_errors          # Non-A errors near 3' end (capped at 10)
    + canonical_junction_bonus     # +1 per canonical GT-AG junction
    + annotated_junction_bonus     # +2 per junction matching annotation
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

mapPacBio forces mismatches at splice junction boundaries (a deliberate design choice for accuracy). This inflates the apparent 3' soft-clip length. RECTIFY uses `effective_three_prime_clip` to subtract these forced-mismatch bases from the penalty.

---

## Selection logic

```python
def select_best_alignment(read_group, annotated_junctions=None):
    """
    Per-read selection.

    1. Score each aligner
    2. Select highest scorer
    3. Tiebreaker: canonical GT-AG count
    4. Tiebreaker: annotated junction count
    5. Tiebreaker: majority 3' position vote
    """
```

The winning aligner's BAM record is written to the consensus BAM.

---

## Chimeric consensus (experimental)

With `--chimeric-consensus`, RECTIFY uses a more sophisticated assembly:

1. Find "sync points" where all aligners agree on the alignment
2. Divide the read into segments at sync points
3. Independently select the best aligner per segment
4. Assemble a chimeric CIGAR from the best segments

This can improve accuracy for reads where different aligners are better at different parts, but it is experimental and requires further validation. Off by default.

**Implementation:** `rectify/core/chimeric_consensus.py`

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

## See also

- [Aligner Recommendations](https://github.com/k-roy/RECTIFY/blob/main/docs/ALIGNER_RECOMMENDATIONS.md) — guidance on when to add uLTRA or deSALT
- [Soft-Clip Rescue](softclip_rescue.md) — how 5' soft-clips are rescued before scoring
