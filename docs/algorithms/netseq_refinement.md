# NET-seq Refinement

NET-seq refinement is an **optional** post-correction step that snaps 3' ends to NET-seq signal peaks within the read's ambiguity window. It is only active when `--netseq-dir` is provided to `rectify correct`; without it the step is skipped entirely.

NET-seq runs as the final refinement step, after all main 3' end correction modules (2E pre → 2F → 2A → 2B → 2C → 2D → 2G → 2E main → 2H) have completed. For reads in A-tract regions where the walk-back algorithm cannot uniquely resolve the true CPA position, NET-seq signal is used to guide the final position assignment.

**Implementation:** `rectify/core/netseq/netseq_refiner.py` (`NetseqLoader`,
`refine_with_netseq`), `rectify/core/analyze/deconvolution.py`
(`DEFAULT_PSF`, `load_psf`, `build_convolution_matrix`, `deconvolve_signal`).

---

## Why NET-seq?

NET-seq captures nascent RNA 3' ends with single-nucleotide resolution. Because the reads come from actively transcribing polymerases (not polyadenylated mRNAs), they are not subject to the A-tract ambiguity that plagues poly(A)-enriched libraries.

However, NET-seq has its own artifact: **oligo(A)-spreading**. Nuclear RNA undergoes oligoadenylation (addition of short A tails), which creates downstream signal spreading around CPA sites.

<figure markdown>
  ![Oligo(A) spreading artifact schematic: a true CPA produces a sharp NET-seq peak; oligoadenylation appends short A tails to the nascent RNA, which then get sequenced at offsets +1, +2, +3, … downstream of the true CPA, producing a characteristic spreading footprint.](../figures/oligo_a_spreading.png){ width="500" }
  <figcaption>NET-seq signal at a CPA is not a single sharp peak — oligoadenylation of nascent RNA produces a downstream "spread" of secondary peaks at +1, +2, +3 nt offsets that must be deconvolved before the signal can be used as ground truth.</figcaption>
</figure>

<figure markdown>
  ![NNLS deconvolution schematic: the spread NET-seq signal at a locus is decomposed into a sum of point-spread-function-convolved peak intensities; the deconvolved peaks recover the true single-nucleotide CPA positions used for walk-back tie-breaking.](../figures/oligo_a_deconvolution.png){ width="680" }
  <figcaption>RECTIFY's NNLS deconvolution decomposes the spread NET-seq signal into a sum of point-spread-function-convolved peaks, recovering sharp single-nucleotide CPA positions. These deconvolved peaks then guide walk-back tie-breaking within A-tract regions.</figcaption>
</figure>

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

The point-spread function is an **empirical 0A-site PSF**: it is the observed
offset distribution `P(observed_position = j | true_CPA = i)` measured at
**zero-A control sites** — positions where the genome immediately downstream
has no A's, so all NET-seq signal must originate from the true CPA (no spreading
possible). RECTIFY ships a `DEFAULT_PSF` derived from WT NET-seq (GSE25107 +
GSE159603), peaking at ~54% on the true position; `load_psf()` can load a
sample- or organism-specific PSF table instead.

The observed NET-seq signal `y` is modeled as:

```
y = A · x + ε
```

where:
- `y` = observed signal vector (counts at each genomic position)
- `A` = convolution matrix built from the PSF (`build_convolution_matrix`)
- `x` = true CPA signal (what we want to recover)
- `ε` = noise

RECTIFY solves for `x` using **Non-negative Least Squares (NNLS)** with L2
regularization, which enforces `x ≥ 0` (counts can't be negative):

```python
def deconvolve_signal(observed, convolution_matrix=None, tract_length=None,
                      psf=None, regularization=None):
    """
    NNLS deconvolution of NET-seq signal (L2-regularized).

    Returns the deconvolved signal array (recovered per-position peak
    intensities). Implemented in rectify/core/analyze/deconvolution.py.
    """
```

---

## Refinement workflow

For reads whose corrected 3' end carries a non-zero A-tract `ambiguity_range`
(i.e. it sits inside a genomic A/T-tract), `refine_with_netseq()` operates over
the `[ambiguity_min, ambiguity_max]` window:

1. Extract the NET-seq signal across the ambiguity window
2. Optionally deconvolve it with NNLS to remove the oligo(A)-spreading artifact
3. Find peaks in the (deconvolved) signal
4. Return either a winner-take-all position or a proportional split of the read
   across the most likely true CPA sites (`proportional_split=True`)
5. The resulting `netseq_confidence` can downgrade the read's `confidence`
   (e.g. to `low`) when no clear peak is found

```python
def refine_with_netseq(netseq_loader, chrom, ambiguity_min, ambiguity_max,
                       strand, original_position, use_deconvolution=True,
                       psf=None, proportional_split=True):
    """
    Refine a CPA position using NET-seq signal across the ambiguity window.
    Returns a single assignment dict, or a list of proportional assignments.
    """
```

Signal access goes through `NetseqLoader`, initialized once with NET-seq signal
(BigWig files via `load_bigwig`/`load_directory`, or the bundled per-organism TSV
arrays via `load_bundled`) and reused per read; it keeps a thread-safe LRU cache
(`MAX_CACHE_SIZE = 10000`) to avoid redundant signal lookups.

---

## Bundled data (yeast)

For *S. cerevisiae*, RECTIFY bundles pre-processed NET-seq 3'-end references as
gzipped TSVs (auto-detected when `--Scer` is set):

- `saccharomyces_cerevisiae_netseq_wt.tsv.gz` — wild-type NET-seq 3' ends
- `saccharomyces_cerevisiae_netseq_pan.tsv.gz` — pan-mutant reference (broader coverage)
- `saccharomyces_cerevisiae_atract_netseq.tsv.gz` — A-tract-focused reference

To process custom NET-seq data (other organisms or mutant conditions), run
`rectify netseq` on the NET-seq BAM(s) — it extracts 3' ends, applies NNLS
deconvolution, and writes per-read parquet plus RPM-normalized bedgraphs:

```bash
rectify netseq my_netseq.bam \
    --genome genome.fa.gz \
    --gff genes.gff.gz \
    -o netseq_processed/

# Then point correction at the processed NET-seq directory:
rectify correct reads.bam \
    --genome genome.fa.gz \
    --netseq-dir netseq_processed/ \
    -o corrected_reads.tsv
```

---

## PSF representation

The PSF is **not** a parametric Gaussian — it is an empirical offset
distribution indexed by `offset = observed − true` and normalized to sum to 1.
The bundled `DEFAULT_PSF` (in `rectify/core/analyze/deconvolution.py`) places
the majority of mass at offset 0 with a downstream tail reflecting
oligo(A)-spreading. `load_psf(filepath)` reads an alternative PSF table (the 0A
slice, `acount == 0`) when adapting to a new chemistry or organism, and
`deconvolve_signal(..., regularization=...)` controls the L2 strength applied
during NNLS.

---

## See also

- [rectify netseq command](../user_guide/commands/netseq.md) — process raw NET-seq BAMs
- [3' End Indel Correction](3prime_indel_correction.md) — the primary correction that NET-seq refines
