# Multi-Aligner Consensus

RECTIFY runs the default Tier-1 trio of aligners (and, when opted in, two
Tier-2 aligners) and selects the best alignment per read using a composite
score.

**Implementation:** `rectify/core/align/multi_aligner.py`, `rectify/core/consensus/` (`extract.py`, `scoring.py`, `select.py`)

<figure markdown>
  ![Multi-aligner consensus flowchart: a read is aligned by minimap2, mapPacBio, gapmm2 (and optionally uLTRA + deSALT); each alignment is scored on clip and junction-proximity metrics; the highest-scoring alignment is selected and a chimeric CIGAR can optionally be assembled from per-segment winners.](../figures/multi_aligner_consensus.png){ width="660" }
  <figcaption>The Tier-1 aligners run per read (with optional Tier-2 aligners); each alignment is scored on its 5' clip penalty, 3' non-poly(A) terminal errors, and junction-proximity errors. Canonical splice sites and annotation matches act only as tiebreakers. The highest-scoring alignment wins per read; with <code>--chimeric-consensus</code> the per-segment winners are stitched together into a chimeric CIGAR.</figcaption>
</figure>

---

## Why multiple aligners?

Different aligners make different tradeoffs:

| Tier | Aligner | Strength | Weakness |
|------|---------|---------|---------|
| 1 (default) | minimap2 | Fast; good general splice junction detection | Can be aggressive about junction placement |
| 1 (default) | mapPacBio | Conservative at junctions (forces mismatches) | Slower than minimap2 |
| 1 (default) | gapmm2 | Terminal-exon homopolymer refinement | Less accurate at very short exons |
| 2 (opt-in) | deSALT | High-sensitivity splice aligner | SIGSEGV-on-certain-inputs fallback; NFS/allocator gotchas |
| 2 (opt-in) | uLTRA | Annotation-guided; good for small exons | Requires annotation; slow |

Tier-1 always runs. The splice-aware aligners uLTRA + deSALT default differently
by command: for **long-read `rectify run-all` they are ON by default** (the
default panel is all 5 aligners; disable with `--no-junction-aligners`), whereas
for the bare `rectify align` command `--junction-aligners` defaults to `[]` and
they must be opted in. `gmap` (and winnowmap2 / minisplice) remain opt-in via
`--junction-aligners` everywhere. By running multiple aligners
and selecting per-read, RECTIFY avoids systematic biases that would arise from
committing to one aligner.

---

## Aligner execution

`run_multi_aligner()` launches each enabled aligner in turn, giving each the
full per-job thread budget; there is no in-process parallel pool across
aligners. (Coarse-grained parallelism comes from the chunked SLURM array
pipeline, which fans out chunks and aligners across array tasks — see
[Pipeline order and I/O](../architecture/pipeline_and_io.md).)

!!! warning "deSALT fork-safety"
    deSALT can emit "double free or corruption" under a conda allocator or
    when its scratch file is memory-mapped over NFS. `run_desalt()` mitigates
    this by stripping `LD_LIBRARY_PATH` (forcing system glibc) and placing the
    temp file on a local filesystem; it also falls back to an empty BAM on a
    SIGSEGV so the rest of the panel still completes.

---

## Scoring

Each aligner's output is represented as an `AlignmentInfo` dataclass
(`consensus/extract.py`), scored by `score_alignment()` in
`consensus/scoring.py`. The composite score (higher is better, starting at 0)
is:

```python
score = (
    - 2 * effective_five_prime_clip   # Unrescued 5' clip/terminal-error region (-2/bp; 0 if rescued)
    - min(2 * effective_three_prime_clip, 10)   # Non-poly(A) 3' terminal errors (-2/bp, capped at 10)
    - min(junction_proximity_errors, 10)        # Mismatch/indel within a few bp of a junction (-1/error, capped at 10)
)
# NOTE: three_prime_atract_depth is intentionally NOT scored — an aligner that
# extends further into a genomic A-tract is more informative, not worse; the
# walk-back assigns the CPA regardless. Canonical GT-AG and annotated-junction
# counts are tiebreakers only, never additive score components.
```

### 5' soft-clip penalty

The soft-clip penalty uses **sequence-based rescue**: the clipped bases are
aligned against the upstream exon sequence, and the clip penalty is waived
when `edit_distance / clip_length ≤ max_edit_frac` (default 0.2, i.e. ≤20%).
Rescued clips carry zero penalty; unrescued bases cost −2/bp.

```python
@dataclass
class AlignmentInfo:
    aligner: str
    junctions: List[Tuple[int, int]]        # (start, end) intron pairs
    five_prime_softclip_seq: str            # Explicit 5' soft-clip sequence (for rescue)
    effective_five_prime_clip_seq: str      # mapPacBio forced-mismatch terminal region (rescue fallback)
    effective_five_prime_clip: int          # max(explicit clip, terminal-error region)
    effective_three_prime_clip: int         # Non-poly(A) 3' terminal errors
    junction_proximity_errors: float        # Mismatch/indel within junction_window_bp of a junction
    three_prime_atract_depth: int           # Downstream A's (informational; NOT scored)
    canonical_count: int                    # Canonical GT-AG junctions (tiebreaker only)
    non_canonical_count: int
    corrected_3prime: Optional[int]         # Genome-only 3' estimate (tiebreaker only)
```

### mapPacBio effective clip correction

mapPacBio forces mismatches at splice junction boundaries instead of soft-clipping (a deliberate design choice). This makes its 5' junction representation functionally equivalent to a soft-clip but encoded differently in CIGAR. `effective_five_prime_clip` normalises this so mapPacBio is scored fairly against aligners that do soft-clip.

### 3' non-poly(A) terminal error penalty

`effective_three_prime_clip` captures non-A (plus strand) or non-T (minus strand) terminal errors near the 3' end via a greedy sliding-window scan. This penalises aligners that stop before the true 3' end with non-poly(A) content at the tail. Penalty: −2/bp, capped at 10. This is distinct from the A-tract depth penalty, which penalises going too far *into* a downstream genomic A-tract.

---

## Selection logic

```python
def select_best_alignment(alignments, genome, annotated_junctions=None, tiebreak='rectify'):
    """
    Per-read selection. Consensus operates on RAW (uncorrected) BAMs.
    rectify correct runs ONCE after consensus, on the winning aligner's BAM.

    1. Score each aligner using the composite penalty (5' clip, 3' non-poly(A)
       terminal errors, junction-proximity errors) — all on the raw BAM records
    2. Select the highest scorer
    3. On a score tie, the default 'rectify' (long-read) order applies, in
       priority:
         Tiebreaker 1: most 3' agreement — corrected_3prime matching the
                       majority across aligners
         Tiebreaker 2: most annotated junction matches
         Tiebreaker 3: most canonical GT-AG splice sites

    (The short-read 'compass' tiebreak order — ungapped > more-annotated >
    shorter-intron > canonical — is used only by the short-read entrypoint.)
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

- [Aligner Recommendations](../ALIGNER_RECOMMENDATIONS.md) — guidance on when to add uLTRA or deSALT
- [Soft-Clip Rescue](softclip_rescue.md) — how 5' soft-clips are rescued before scoring
