# RECTIFY: Unified RNA 3' End Correction Framework

[![PyPI version](https://badge.fury.io/py/rectify-rna.svg)](https://pypi.org/project/rectify-rna/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-192%20passing-brightgreen.svg)](tests/)

**RECTIFY** (**R**NA 3' **E**nd **C**orrection **T**ool **I**ntegrating **F**alse-priming and pol**y**(A) ambiguity) is a unified framework for correcting 3' end mapping artifacts in poly(A)-tailed RNA sequencing data.

## Overview

RECTIFY addresses two fundamental problems affecting RNA 3' end mapping:

1. **A-tract Ambiguity (Universal)**: Genomic A-tracts near true 3' ends create positional uncertainty affecting ALL poly(A)-tailed RNA-seq technologies
2. **Technology-Specific Artifacts**:
   - **AG mispriming**: Internal priming on A/G-rich regions (oligo-dT methods)
   - **Poly(A) tail alignment**: Tail bases align to genomic A-tracts creating systematic shifts (when poly(A) is sequenced)

### Key Features

- **Modular correction strategies** that apply based on sequencing technology
- **Universal A-tract ambiguity detection** for all poly(A)-tailed RNA-seq
- **AG mispriming screening** (from original RECTIFY, Roy & Chanfreau 2019)
- **Poly(A) tail trimming and indel artifact correction** (for direct RNA-seq: nanopore, Helicos, QuantSeq)
- **NET-seq refinement** (optional, technology-independent)
- **Unified output format** with confidence scores and QC flags

## How It Works

RECTIFY corrects 3' end mapping artifacts by walking a read through multiple correction steps.
Here's a single Nanopore read traversing the full pipeline:

```
===============================================================================
STEP 1: RAW ALIGNMENT
===============================================================================

The aligner maps a Nanopore direct RNA read to the genome. The poly(A) tail
is soft-clipped where the genomic A-tract ends, but the aligned region of the tail
contains additional errors due to alignment and sequencing errors.

genome: 5'..CTAGTGACAGTCAAAAAAAA-AAACAAAAGTAAAAAAAAAAAA|CTAGCGATC..3'
                                                       |
                                              CPA site (true 3' end)

read:   5'..CTAGTGACAGTCAAAAAAAATAAA-AAAAA--AAAAAAAAAAAAAAAAAAAAA..
                                 ^      ^^  |___________________|
                              T error  dels   soft-clipped tail

Key observations:
  - Soft-clip boundary placed at mapped 3' end
  - Deletions (-) = aligner attempts to maximize alignment of the poly(A) tail
    at expense of indels
  - 'T' in the poly(A) tail = Nanopore sequencing error in the homopolymer

===============================================================================
STEP 2: INDEL CORRECTION - Walk Backwards to Find True 3' End
===============================================================================

Starting from soft-clip boundary, walk upstream comparing genome vs read
to find the first position where both agree on a non-A base:

genome: ..CTAGTGACAGTCAAAAAAAA-AAACAAAAGTAAAAAAAAAAAA|..
read:   ..CTAGTGACAGTCAAAAAAAATAAA-AAAAA--AAAAAAAAAA.|..
                     ^         ^      ^^
                  AGREE     T error  dels
                  (C=C)    (ignore) (count)

RECTIFY walks back from soft-clip boundary:
  - Skip all A's (ambiguous with poly(A) tail)
  - Ignore deletions where mRNA has the poly(A) tail: 3 bp total
  - Ignore T (likely sequencing error in homopolymer)
  - Find first non-A agreement: 'C' at upstream position
  - Result: ambiguity window spanning the A-tract from upstream C to downstream C

===============================================================================
STEP 3: NET-seq REFINEMENT (Optional)
===============================================================================

For species with NET-seq data, we can resolve the ambiguity window.

NET-seq captures nascent RNA bound to RNA Pol II. Cleavage intermediates
are oligo-adenylated before release, with a distribution of tail lengths:

                           CPA site (cleavage here)
                                  |
                                  v
    ___________                  ✂️
   |  Pol II   |====~~~CGUACGUAG*AAA...     <- oligo-A tail added
   |___________|                 3'            after cleavage

Oligo(A) tail length distribution (typical):
   1A:  ████████████  ~12%
   2A:  ██████████    ~10%
   3A:  █████████     ~9%
   4A:  ████████      ~8%
   5A:  ███████       ~7%
   6A:  ██████        ~6%    mean ~6.6 A's
   7A:  █████         ~5%
   8A:  ████          ~4%
   9A:  ███           ~3%
  10A:  ██            ~2%
  ...decreasing...

When oligo-adenylated intermediates align to a genome with downstream A's,
the oligo-A tail extends the alignment into the genomic A-tract:

                        true CPA
                           |
genome:       ...CGTACGTAG|AAAAAAAA|GTCACC...
                          |^^^^^^^^|
                          genomic A-tract
                                   |
aligned:      ...CGTACGTAG*AAAAAAA      <- oligo-A aligns to genomic A's
                          |<------>|
                           SHIFT (apparent 3' end moves downstream)

This creates a spreading artifact: signal is shifted downstream.

![Oligo(A) Spreading Artifact](docs/figures/oligo_a_spreading.png)
*Figure: Oligo(A) spreading artifact showing three CPA sites within a 12bp A-tract. Each peak has a right-tailed distribution where signal "spills" downstream. Signal beyond the A-tract boundary is soft-clipped.*

RECTIFY uses NNLS (Non-Negative Least Squares) deconvolution to remove
the spreading artifact and recover true peak positions:

1. Build convolution matrix from 0A PSF (Point-Spread Function)
   - ~54% of signal stays at true position
   - ~46% spreads downstream (mean ~3bp)

2. Solve: observed = A @ true_peaks (regularized NNLS)

3. Assign reads PROPORTIONALLY to deconvolved peaks:
   - If 2 peaks with 70%/30% signal → assign 0.7 and 0.3 reads

Confidence assignment:
  - HIGH:   Single dominant peak (>90% of signal)
  - MEDIUM: Dominant peak (>70% of signal)
  - SPLIT:  Multiple significant peaks (no dominant)
  - LOW:    Weak signal or no NET-seq data

![Deconvolution](docs/figures/oligo_a_deconvolution.png)
*Figure: NNLS deconvolution "adds back" the spread oligo-A tails to their true CPA positions. Left: observed signal with spreading. Right: deconvolved signal with all reads assigned to true peaks.*

===============================================================================
STEP 4: FINAL OUTPUT (with proportional apportionment)
===============================================================================

When NET-seq reveals multiple peaks within an ambiguity window, reads are
APPORTIONED PROPORTIONALLY to each peak position. This preserves quantitative
accuracy for downstream analysis.

Example 1: Single dominant peak (HIGH confidence)
-------------------------------------------------
Ambiguity window [31, 42], NET-seq shows single peak at position 35:

                 31              35                          42
                 |               |                           |
Ambiguity:       |===============================================|
NET-seq peak:                    #

Output: read001 → position 35, weight 1.0, confidence HIGH

Example 2: Multiple peaks (SPLIT confidence)
--------------------------------------------
A SINGLE Nanopore read with ambiguity window [31, 42] can be refined using
NET-seq to reveal THREE distinct CPA sites at positions 33, 36, and 40:

                 31   33   36        40                       42
                 |    |    |         |                        |
Ambiguity:       |==============================================|
NET-seq peaks:        ###  ##        #
                      50%  30%       20%

Output: read001 → split into THREE output rows (one read becomes three):
  - read001 at position 33, weight 0.5
  - read001 at position 36, weight 0.3
  - read001 at position 40, weight 0.2

This preserves quantitative accuracy: if 100 reads map to this ambiguous
region, the output will have ~50 reads at pos 33, ~30 at pos 36, ~20 at pos 40.

Final output format:
+----------+----------+--------+------------+-----------------+
| read_id  | position | weight | confidence | ambiguity_range |
+----------+----------+--------+------------+-----------------+
| read001  |    33    |  0.5   |   SPLIT    |       11        |
| read001  |    36    |  0.3   |   SPLIT    |       11        |
| read001  |    40    |  0.2   |   SPLIT    |       11        |
+----------+----------+--------+------------+-----------------+

Without NET-seq: Reports ambiguity window [31, 42] and uses leftmost
position (most conservative estimate), weight 1.0.
```

## Installation

### From source (development)

```bash
git clone https://github.com/k-roy/RECTIFY.git
cd RECTIFY
pip install -e .
```

### From PyPI

```bash
pip install rectify-rna
```

## Quick Start

### QuantSeq (oligo-dT short-read)

```bash
rectify correct quantseq.bam \
  --genome sacCer3.fa \
  --annotation genes.gtf \
  --polya-sequenced \
  --output corrected_3ends.tsv
```

### Nanopore direct RNA-seq with NET-seq refinement

```bash
rectify correct nanopore.bam \
  --genome sacCer3.fa \
  --annotation genes.gtf \
  --polya-sequenced \
  --aligner minimap2 \
  --netseq-dir churchman_bigwigs/ \
  --output corrected_3ends.tsv
```

## Output Format

RECTIFY produces two output files:

### 1. Corrected Positions TSV (`output.tsv`)

Per-read corrected 3' end positions with QC metrics:

```
read_id  chrom  strand  original_3prime  corrected_3prime  ambiguity_min  ambiguity_max  ambiguity_range  polya_length  correction_applied  confidence  qc_flags
read001  chrI   +       147588           147585            147583         147588         5                42            atract_ambiguity    high        PASS
read002  chrI   +       147593           147591            147591         147593         2                38            atract_ambiguity    medium      AG_RICH
```

**Key columns:**
- `polya_length`: Total observed poly(A) tail length = aligned A's (from A-tract) + soft-clipped A's
- `ambiguity_range`: Positional uncertainty due to poly(A) tail aligning to genomic A-tract
- `correction_applied`: `atract_ambiguity` (position correction), `indel_correction`, `netseq_refinement`

### 2. Processing Statistics TSV (`output_stats.tsv`)

Comprehensive read flow statistics through each filtering and correction stage:

```
metric                       count      percent   description
total_reads_in_bam           1000000    100.00    Total reads in BAM file
reads_unmapped               5000       0.50      Unmapped reads (skipped)
reads_secondary              2000       0.20      Secondary alignments (skipped)
reads_processed              993000     99.30     Reads with 3' ends corrected
ends_with_downstream_A       450000     45.32     3' ends with >=1 A immediately downstream
ends_ambiguous_atract        180000     18.13     3' ends in A-tract (ambiguity range > 0)
ends_shifted_atract_walking  85000      8.56      3' ends shifted by A-tract walking
reads_with_polya             890000     89.63     Reads with detected poly(A) tail
polya_length_mean            42.3       -         Mean poly(A) tail length (bp)
polya_length_max             185        -         Maximum poly(A) tail length (bp)
confidence_high              750000     75.53     High confidence assignments
confidence_medium            200000     20.14     Medium confidence assignments
confidence_low               43000      4.33      Low confidence assignments
```

## Module Architecture

RECTIFY applies corrections modularly based on your data:

1. **Module 1: A-tract Ambiguity** (always applied)
   - Identifies genomic A-tracts near 3' ends
   - Calculates ambiguity windows

2. **Module 2A: AG Mispriming** (when oligo-dT priming used)
   - Screens for downstream AG-richness
   - Flags likely misprimed reads

3. **Module 2B: Poly(A) Detection** (when poly(A) IS sequenced)
   - Measures observed poly(A) tail length (aligned A's + soft-clipped A's)
   - Reports tail length for QC (does NOT correct position - that's handled by A-tract detection)

4. **Module 2C: Indel Correction** (when poly(A) IS sequenced)
   - Detects and corrects indel artifacts in aligned region

5. **Module 3: NET-seq Refinement** (optional)
   - Resolves ambiguity using NET-seq data via NNLS deconvolution
   - Apportions reads PROPORTIONALLY to multiple peaks when ambiguity exists
   - Assigns confidence scores (HIGH/MEDIUM/SPLIT/LOW)

## Downstream Analysis

RECTIFY includes a comprehensive `analyze` command for downstream analysis of corrected 3' ends:

```bash
rectify analyze corrected_3ends.tsv \
  --annotation genes.gtf \
  --output-dir analysis_results/
```

### Analysis Modules

| Module | Description | Output |
|--------|-------------|--------|
| **Clustering** | Group nearby CPA sites into clusters | `clusters.tsv`, count matrices |
| **Differential Expression** | DESeq2-based gene and cluster-level analysis | `deseq2_results.tsv` |
| **PCA** | Sample quality control and batch effect detection | `pca_plot.png` |
| **Shift Analysis** | Detect condition-specific CPA site usage changes | `shift_results.tsv`, browser plots |
| **GO Enrichment** | Functional enrichment of genes with shifted 3' ends | `go_enrichment.tsv` |
| **Motif Discovery** | Identify sequence motifs near CPA sites | `motif_results/` |

### Browser-Style Visualization

The shift analysis module generates publication-ready genome browser plots showing:
- Per-condition CPA site usage distributions
- Gene track with CDS boxes as directional pentagon arrows
- CPA cluster markers

## Citation

If you use RECTIFY, please cite:

**Original RECTIFY (AG mispriming correction):**
> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020 Apr 1;176:4-13. doi: 10.1016/j.ymeth.2019.05.016. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.0 (unified framework):**
> Manuscript in preparation

## License

MIT License - See [LICENSE](LICENSE) for details

## Contact

- Kevin R. Roy - [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com)
- GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)

## Acknowledgments

- Original RECTIFY development supported by Chanfreau Lab, UCLA
- NET-seq data from Churchman Lab, Harvard Medical School
