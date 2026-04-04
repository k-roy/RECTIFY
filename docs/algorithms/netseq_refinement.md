# NET-seq Refinement

For reads in A-tract regions where the walk-back algorithm cannot uniquely resolve the true CPA position, RECTIFY uses NET-seq (Native Elongating Transcript sequencing) data to guide the final position assignment.

**Implementation:** `rectify/core/netseq_refiner.py`, `rectify/core/netseq_deconvolution.py`

---

## Why NET-seq?

NET-seq captures nascent RNA 3' ends with single-nucleotide resolution. Because the reads come from actively transcribing polymerases (not polyadenylated mRNAs), they are not subject to the A-tract ambiguity that plagues poly(A)-enriched libraries.

However, NET-seq has its own artifact: **oligo(A)-spreading**. Nuclear RNA undergoes oligoadenylation (addition of short A tails), which creates downstream signal spreading around CPA sites.

---

## Oligo(A) spreading artifact

Raw NET-seq signal around a CPA site looks like:

```
True CPA
    ↓
────|────────────────────────────────────────── Genome
    ↑↑↑↑↑↑↑↑↑↑
NET-seq signal spreading into A-tract
(each A is a possible oligo-A extension endpoint)
```

The spreading width depends on the downstream A-tract length. The point-spread function (PSF) describes this spreading distribution.

---

## NNLS deconvolution

RECTIFY derives the PSF from 5,000+ **zero-A control sites** — positions where the genome immediately downstream has no A's. At these sites, all NET-seq signal must originate from the true CPA (no spreading possible), so the signal width measures the PSF directly.

The observed NET-seq signal `y` is modeled as:

```
y = A · x + ε
```

where:
- `y` = observed signal vector (counts at each genomic position)
- `A` = convolution matrix derived from the PSF
- `x` = true CPA signal (what we want to recover)
- `ε` = noise

RECTIFY solves for `x` using **Non-negative Least Squares (NNLS)**, which enforces `x ≥ 0` (counts can't be negative):

```python
from scipy.optimize import nnls

def deconvolve_signal(observed, psf, region_start, region_end):
    """
    NNLS deconvolution of NET-seq signal.

    Returns: deconvolved_signal, residual
    """
    A = build_convolution_matrix(psf, len(observed))
    x, residual = nnls(A, observed)
    return x, residual
```

---

## Refinement workflow

For reads flagged as `ATRACT_AMBIGUOUS` (A-count downstream > 15):

1. Look up NET-seq signal in the ±50 bp window around the candidate CPA positions
2. Apply NNLS deconvolution
3. Find the peak of the deconvolved signal
4. If the deconvolved peak is within the ambiguous window: set corrected position to peak
5. If no clear peak: retain the walk-back position with `LOW` confidence

```python
class NetseqRefiner:
    """
    Resolves A-tract ambiguity using pre-loaded NET-seq data.

    The refiner is initialized once with the NET-seq BigWig files
    and called per-read. Uses an LRU cache (default 10,000 entries)
    to avoid redundant BigWig lookups.
    """

    def refine_position(self, chrom, strand, candidate_positions, genome):
        """
        Select best CPA position from candidates using NET-seq.

        Returns: (position, confidence_boost)
        """
```

---

## Bundled data (yeast)

For *S. cerevisiae*, RECTIFY bundles pre-processed WT NET-seq BigWigs:

- `netseq/wt_plus.bw` — plus strand
- `netseq/wt_minus.bw` — minus strand

These are deconvolved from oligo(A)-spreading and represent the true CPA signal. Auto-detected when `--Scer` is set.

For custom NET-seq data (other organisms or mutant conditions):

```bash
rectify netseq my_netseq.bam \
    --genome genome.fa.gz \
    --gff genes.gff.gz \
    -o netseq_processed/

# Then use the processed BigWigs:
rectify correct reads.bam \
    --genome genome.fa.gz \
    --netseq-dir netseq_processed/ \
    -o corrected.tsv
```

---

## PSF parameters

The PSF is characterized by:

- **σ (sigma)**: standard deviation of the Gaussian spreading (default: 25 bp)
- **Window**: ±3σ = ±75 bp around each CPA site
- **Min signal**: minimum NET-seq counts to accept a deconvolution result

The PSF can be updated for a new organism by training on zero-A control sites:

```python
from rectify.core.analyze.deconvolution import fit_psf_from_controls

psf = fit_psf_from_controls(netseq_bigwigs, zero_a_sites, genome)
```

---

## See also

- [rectify netseq command](../user_guide/commands/netseq.md) — process raw NET-seq BAMs
- [3' End Indel Correction](3prime_indel_correction.md) — the primary correction that NET-seq refines
