# RECTIFY: Unified RNA 3' End Correction Framework

[![PyPI version](https://badge.fury.io/py/rectify-rna.svg)](https://pypi.org/project/rectify-rna/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-192%20passing-brightgreen.svg)](tests/)

**Correct poly(A) tail artifacts and perform comprehensive differential expression analysis with a single command.**

---

## Quick Start

### Install

```bash
pip install rectify-rna
```

### Run

```bash
# Simplest: FASTQ input with bundled yeast genome (no external files needed!)
rectify correct reads.fastq.gz --organism yeast -o corrected.tsv

# All-in-one: correct + analyze (uses bundled WT NET-seq for yeast by default)
rectify run reads.bam --genome genome.fa --annotation genes.gtf --output-dir results/
```

**New in v2.2.0:** RECTIFY now includes:
- **Bundled yeast genome** (S288C R64-5-1) and GFF annotations from SGD
- **Pre-processed WT NET-seq data** (Churchman lab)
- **FASTQ support** with automatic minimap2 alignment
- **CPA motif database** for transcription factor binding site matching

<details>
<summary>Or run steps separately / use custom NET-seq data</summary>

```bash
# 1. Correct 3' end positions (bundled WT NET-seq for yeast)
rectify correct reads.bam --genome genome.fa --organism yeast --output corrected.tsv

# With custom/mutant NET-seq data (overrides bundled WT data)
rectify correct reads.bam --genome genome.fa --netseq-dir my_mutant_netseq/ --output corrected.tsv

# 2. Analyze results (clustering, differential expression, GO enrichment, enriched motifs around cluster peaks)
rectify analyze corrected.tsv --annotation genes.gtf --output-dir results/

# With explicit sample metadata (optional - conditions auto-inferred from sample names like WT_rep1, KO_rep2)
rectify analyze corrected.tsv --annotation genes.gtf --manifest samples.tsv --reference WT --output-dir results/
```

</details>

---

## What You Get

### Corrected 3' End Positions

Each read gets a corrected position with confidence scores and QC metrics:

```
┌───────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│  read_id   │ chrom │ strand │ original_3prime │ corrected_3prime │ shift │ confidence │ polya_length │ qc_flags │
├───────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│  read001   │ chrI  │   +    │     147592      │     147585       │  -7   │    HIGH    │      42      │   PASS   │
│  read002   │ chrI  │   +    │     147594      │     147591       │  -3   │   MEDIUM   │      38      │   PASS   │
│  read003   │ chrI  │   -    │     147602      │     147602       │   0   │    HIGH    │      55      │   PASS   │
│  read004   │ chrII │   +    │     283109      │     283104       │  -5   │    LOW     │      31      │ AG_RICH  │
└───────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

**Key fields:**
- `original_3prime`: Mapped 3' end position from aligner
- `corrected_3prime`: True CPA site (shifted upstream to first non-A)
- `shift`: Correction applied (negative = shifted upstream through A-tract)
- `confidence`: HIGH (single peak), MEDIUM (dominant peak), SPLIT (multiple peaks), LOW (uncertain)
- `polya_length`: Measured poly(A) tail length (aligned + soft-clipped A's)
- `qc_flags`: PASS, AG_RICH (possible mispriming), or other warnings

### Processing Statistics

Comprehensive QC report showing read flow through each correction stage:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  RECTIFY Processing Summary                                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│  Total reads processed         │  993,000   │  99.3%                        │
│  Reads with poly(A) detected   │  890,000   │  89.6%                        │
│  Positions corrected           │  180,000   │  18.1%                        │
│  Mean poly(A) length           │     42.3 bp                                │
├─────────────────────────────────────────────────────────────────────────────┤
│  Confidence Distribution                                                    │
│    HIGH                        │  750,000   │  ████████████████████  75.5%  │
│    MEDIUM                      │  200,000   │  █████                 20.1%  │
│    LOW                         │   43,000   │  █                      4.3%  │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Analysis Outputs (`rectify analyze`)

**Cluster-level counts** for differential expression:

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  cluster_id  │  gene   │ chrom │ strand │  start  │   end   │ sample_A │ sample_B │ ... │
├────────────────────────────────────────────────────────────────────────────────────────┤
│  cluster_1   │  ACT1   │ chrVI │   +    │ 54696   │ 54702   │   1523   │   1891   │     │
│  cluster_2   │  ACT1   │ chrVI │   +    │ 54715   │ 54718   │    342   │    128   │     │
│  cluster_3   │  PGK1   │ chrIII│   -    │ 137218  │ 137225  │   2104   │   2087   │     │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

**DESeq2 differential expression** (gene and cluster level):

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│   gene   │ baseMean │ log2FoldChange │  lfcSE  │  stat   │  pvalue  │    padj     │
├──────────────────────────────────────────────────────────────────────────────────────┤
│  RPL3    │  8542.3  │     1.82       │  0.12   │  15.2   │ 2.1e-52  │  4.8e-49    │
│  HSP82   │  3201.7  │    -2.14       │  0.18   │ -11.9   │ 1.3e-32  │  1.5e-29    │
│  ENO2    │  5876.2  │     0.92       │  0.09   │  10.2   │ 1.8e-24  │  1.4e-21    │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

**3' UTR shift analysis** (alternative polyadenylation):

```
┌──────────────────────────────────────────────────────────────────────────────────────────┐
│   gene   │ proximal_cluster │ distal_cluster │ shift_score │ direction │    padj     │
├──────────────────────────────────────────────────────────────────────────────────────────┤
│  FAS1    │    cluster_12    │   cluster_14   │    0.42     │  SHORTEN  │  3.2e-08    │
│  CDC19   │    cluster_28    │   cluster_31   │   -0.38     │  LENGTHEN │  1.7e-06    │
│  TDH3    │    cluster_45    │   cluster_47   │    0.25     │  SHORTEN  │  4.1e-04    │
└──────────────────────────────────────────────────────────────────────────────────────────┘
```

**Visualization outputs:**
- `pca_plot.png` - Sample clustering and batch effect detection
- `heatmap.png` - Gene expression heatmap with hierarchical clustering
- `shift_browser_*.png` - Genome browser plots showing CPA site usage per condition
- `motif_results/` - Sequence motifs enriched near CPA sites (MEME format)

---

## Overview

RECTIFY addresses two fundamental problems affecting RNA 3' end mapping:

1. **A-tract Ambiguity (Universal)**: Genomic A-tracts near true 3' ends create positional uncertainty affecting ALL poly(A)-tailed RNA-seq technologies
2. **Technology-Specific Artifacts**:
   - **AG mispriming**: Internal priming on A/G-rich regions (oligo-dT methods)
   - **Poly(A) tail alignment**: Tail bases align to genomic A-tracts creating systematic shifts (when poly(A) is sequenced)

### Indel Correction in A-rich Regions

When poly(A) tails align to genomic A-tracts, aligners like minimap2 introduce indel artifacts to maximize alignment score. RECTIFY detects and corrects these artifacts.

![Indel Correction](docs/figures/indel_correction.png)

**Key Insight for IGV Users:** When viewing minus strand reads in IGV:
- The poly(A) tail appears as **poly(T)** (because IGV shows the read sequence, not the RNA)
- The tail extends **leftward** (toward lower genomic coordinates)
- RECTIFY corrects by shifting **rightward** (toward the true CPA site)

**The Problem:** Minimap2 introduces deletions to extend alignment of poly(A) tail bases into the genomic A-tract. This shifts the apparent 3' end **downstream** of the true CPA site.

**RECTIFY's Solution:** Walk backwards from the soft-clip boundary to find the first non-A position where genome and read agree:

1. Start at soft-clip boundary (poly(A) tail begins)
2. Walk upstream through aligned region:
   - Skip A/T positions (ambiguous with poly(A) tail)
   - Absorb deletions (D) - these are alignment artifacts
   - Absorb T errors - likely sequencing errors in homopolymer
   - STOP at first non-A/T agreement
3. Calculate corrected position and total poly(A) length

**Result:** The true CPA site is shifted **upstream** by the number of absorbed A's and deletions. Reads that appeared to end at different positions within the A-tract are unified to the true cleavage site.

### NET-seq Refinement

For organisms with available NET-seq data, RECTIFY can resolve A-tract ambiguity by using nascent RNA 3' end positions to infer true CPA sites. The bundled WT NET-seq data is auto-detected based on your genome:

| Organism | NET-seq Data | Auto-detected from |
|----------|--------------|-------------------|
| *S. cerevisiae* | Churchman lab (bundled) | Genome size ~12 Mb, chrI-XVI |
| Other eukaryotes | User-provided | Use `--netseq-dir` |

**Without NET-seq**: RECTIFY still corrects A-tract artifacts and reports ambiguity windows - downstream analysis remains valid, just with wider confidence intervals.

**Custom NET-seq**: Provide your own NET-seq data (e.g., from mutant strains) with `--netseq-dir` to override the bundled WT data for condition-specific refinement.

![Oligo(A) Spreading Artifact](docs/figures/oligo_a_spreading.png)

**Poly(A) tails cause systematic 3' end mapping errors.** When poly(A) tails align to genomic A-tracts, the apparent 3' end shifts downstream. RECTIFY corrects this artifact using NNLS deconvolution:

![Deconvolution](docs/figures/oligo_a_deconvolution.png)

---

## Features

| Feature | Description |
|---------|-------------|
| **A-tract Correction** | Detects genomic A-tracts and calculates ambiguity windows |
| **Poly(A) Measurement** | Reports observed tail length (aligned + soft-clipped) |
| **Indel Correction** | Fixes alignment artifacts in poly(A) regions |
| **AG Mispriming Detection** | Flags likely internal priming events (oligo-dT methods) |
| **NET-seq Refinement** | Resolves ambiguity using nascent RNA data (optional) |
| **Proportional Assignment** | Splits reads across multiple CPA sites when appropriate |
| **FASTQ Support** | Direct FASTQ input with automatic minimap2 alignment |
| **Bundled Genomes** | Yeast S288C genome and annotations included |
| **Motif Database** | CPA-related TF binding sites for motif matching |

### Bundled Data (New in v2.2.0)

For yeast (*S. cerevisiae*), RECTIFY includes:

| Data | Description | Size |
|------|-------------|------|
| **Genome** | S288C R64-5-1 reference | 3.7 MB |
| **Annotation** | SGD GFF3 with gene names | 5.0 MB |
| **GO Annotations** | Gene Ontology from SGD | 400 KB |
| **NET-seq** | WT pan-mutant consensus | ~1 MB |
| **Motif Database** | CPA factors, NNS pathway, general TFs | 10 KB |

### Supported Technologies

- **Nanopore direct RNA-seq** (minimap2)
- **QuantSeq** (oligo-dT short-read)
- **Helicos** (single-molecule)
- **PacBio Iso-Seq**
- **Any poly(A)-tailed RNA-seq**

---

## Export to bedGraph/bigWig

Generate genome browser tracks from corrected 3' ends:

```bash
# Export per-replicate and per-condition bigWig files
rectify export corrected.tsv -o tracks/ --genome genome.fa

# Or use bedGraph format
rectify export corrected.tsv -o tracks/ --format bedgraph

# Per-replicate only
rectify export corrected.tsv -o tracks/ --per-replicate

# Per-condition summed only
rectify export corrected.tsv -o tracks/ --per-condition
```

**Output structure:**
```
tracks/
├── per_replicate/
│   ├── wt_rep1.plus.bw
│   ├── wt_rep1.minus.bw
│   ├── wt_rep2.plus.bw
│   └── ...
└── per_condition/
    ├── wt.plus.bw
    ├── wt.minus.bw
    ├── mutant.plus.bw
    └── ...
```

---

## Downstream Analysis

RECTIFY includes a comprehensive `analyze` command:

```bash
rectify analyze corrected.tsv --annotation genes.gtf --output-dir results/
```

| Module | Description | Output |
|--------|-------------|--------|
| **Clustering** | Group nearby CPA sites | `clusters.tsv`, count matrices |
| **Differential Expression** | DESeq2-based analysis | `deseq2_results.tsv` |
| **PCA** | Sample QC and batch effects | `pca_plot.png` |
| **Shift Analysis** | Condition-specific CPA usage | `shift_results.tsv` |
| **GO Enrichment** | Functional enrichment | `go_enrichment.tsv` |
| **Motif Discovery** | Sequence motifs near CPA sites | `motif_results/` |

---

## Installation Options

### From PyPI (recommended)

```bash
pip install rectify-rna
```

### From Conda (recommended)

```bash
# Bioconda (recommended - includes MEME Suite for motif analysis)
conda install -c conda-forge -c bioconda rectify-rna

# Or with mamba (faster)
mamba install -c conda-forge -c bioconda rectify-rna
```

> **Note**: The conda installation includes MEME Suite (MEME, STREME, FIMO) for motif discovery.
> If using pip, install MEME separately: `conda install -c bioconda meme>=5.4.0`

### From source (development)

```bash
git clone https://github.com/k-roy/RECTIFY.git
cd RECTIFY
pip install -e .
```

> **Note**: This is a beta release. Please report issues at [GitHub Issues](https://github.com/k-roy/RECTIFY/issues).

---

## Command Reference

### Basic Usage

```bash
# Simplest: FASTQ input with bundled yeast genome
rectify correct reads.fastq.gz --organism yeast -o corrected.tsv

# BAM input with custom genome
rectify correct reads.bam --genome genome.fa --output corrected.tsv

# With annotation (recommended)
rectify correct reads.bam --genome genome.fa --annotation genes.gtf --output corrected.tsv

# Nanopore with NET-seq refinement
rectify correct reads.bam \
  --genome genome.fa \
  --annotation genes.gtf \
  --aligner minimap2 \
  --netseq-dir netseq_bigwigs/ \
  --output corrected.tsv

# QuantSeq (oligo-dT)
rectify correct reads.bam \
  --genome genome.fa \
  --annotation genes.gtf \
  --polya-sequenced \
  --output corrected.tsv
```

### Key Options

| Option | Description |
|--------|-------------|
| `input` | Input BAM or FASTQ file (FASTQ auto-aligned with minimap2) |
| `--genome` | Reference genome FASTA (optional with `--organism yeast`) |
| `--organism` | Use bundled genome/annotation for organism (e.g., `yeast`) |
| `--annotation` | Gene annotation GTF/GFF |
| `--netseq-dir` | Directory with NET-seq bigWig files |
| `--aligner` | Aligner used: `minimap2`, `star`, `bowtie2` |
| `--polya-sequenced` | Poly(A) tail was sequenced (not just primed) |
| `--threads` | Number of threads (default: 4) |

---

## How It Works (Technical Details)

<details>
<summary>Click to expand detailed algorithm description</summary>

### The Correction Pipeline

RECTIFY corrects 3' end mapping artifacts by walking each read through multiple correction steps:

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

When poly(A) tails align to genomic A-tracts, aligners like minimap2 introduce
indels (insertions/deletions) to maximize the alignment score. These artifacts
shift the apparent 3' end downstream.

Real Example: Read SRR32518284.2355129 at chrI:31,393 (IGV screenshot)
------------------------------------------------------------------------
This minus-strand read shows typical minimap2 indel artifacts in an A-rich region:

Position:      31,395  31,400  31,405  31,410  31,415  31,420
               |       |       |       |       |       |
Genome (ref):  GAAT TTTTT T CTCAT AAA GAAA AA T AAA G...
                        |       |||
Read aligned:  GAAT-TTTTT T CTCAT AAA GAAA--AA T AAA G...
                   ^           ^^^       ^^
                   1D          3bp       2D
               (deletion)   (aligned) (deletions)

The aligner introduces deletions (shown as '-' or numbers in IGV) to
maximize alignment of the poly(A) tail bases to genomic A's.

RECTIFY's correction algorithm:
------------------------------
Starting from the soft-clip boundary, walk UPSTREAM comparing genome vs read
to find the first position where both agree on a non-A/non-T base:

STEP A: Start at soft-clip boundary (rightmost aligned position)
        Count soft-clipped A's as part of poly(A) tail

STEP B: Walk upstream through aligned region:
        - Skip positions where genome = A (or T for minus strand)
        - Skip deletions (D in CIGAR) - these are tail alignment artifacts
        - Skip insertions (I in CIGAR) of A/T - these are sequencing errors
        - Stop when genome and read AGREE on a non-A/T base

STEP C: Calculate total poly(A) length:
        = soft-clipped A's + aligned A's + deletion-absorbed A's

Example walkback:
                                          soft-clip boundary
                                                   |
Genome: ...CTAGTGACAGTCAAAAAAAA-AAACAAAAGTAAAAAAAAAAAA|CTAGCGATC...
Read:   ...CTAGTGACAGTCAAAAAAAATAAA-AAAAA--AAAAAAAAAA.|AAAAAAAAAAA
                     ^         ^      ^^             |<-------->
                  AGREE     T error  dels          soft-clipped
                  (C=C)    (ignore) (count)           tail

Walk back from '|':
  pos-1: A (genome) = A (read) → continue (ambiguous)
  pos-2: A (genome), deletion → count as absorbed
  ...
  pos-N: C (genome) = C (read) → STOP (first non-A agreement)

Result: True 3' end at position of C
        Poly(A) length = 11 soft-clipped + 3 dels + 12 aligned A's = 26 bp

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

RECTIFY uses NNLS (Non-Negative Least Squares) deconvolution to remove
the spreading artifact and recover true peak positions:

1. Build convolution matrix from oligo(A) PSF (Point-Spread Function)
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

===============================================================================
STEP 4: FINAL OUTPUT (with proportional apportionment)
===============================================================================

When NET-seq reveals multiple peaks within an ambiguity window, reads are
APPORTIONED PROPORTIONALLY to each peak position. This preserves quantitative
accuracy for downstream analysis.

Example: Single dominant peak (HIGH confidence)
-------------------------------------------------
Ambiguity window [31, 42], NET-seq shows single peak at position 35:

                 31              35                          42
                 |               |                           |
Ambiguity:       |===============================================|
NET-seq peak:                    #

Output: read001 → position 35, weight 1.0, confidence HIGH

Example: Multiple peaks (SPLIT confidence)
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

Without NET-seq: Reports ambiguity window [31, 42] and uses leftmost
position (most conservative estimate), weight 1.0.
```

### Module Architecture

RECTIFY applies corrections modularly based on your data:

1. **Module 1: A-tract Ambiguity** (always applied)
   - Identifies genomic A-tracts near 3' ends
   - Calculates ambiguity windows

2. **Module 2A: AG Mispriming** (when oligo-dT priming used)
   - Screens for downstream AG-richness
   - Flags likely misprimed reads

3. **Module 2B: Poly(A) Detection** (when poly(A) IS sequenced)
   - Measures observed poly(A) tail length (aligned A's + soft-clipped A's)
   - Reports tail length for QC

4. **Module 2C: Indel Correction** (when poly(A) IS sequenced)
   - Detects and corrects indel artifacts in aligned region

5. **Module 3: NET-seq Refinement** (optional)
   - Resolves ambiguity using NET-seq data via NNLS deconvolution
   - Apportions reads PROPORTIONALLY to multiple peaks
   - Assigns confidence scores (HIGH/MEDIUM/SPLIT/LOW)

</details>

---

## Citation

If you use RECTIFY, please cite:

**Original RECTIFY (AG mispriming correction):**
> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020 Apr 1;176:4-13. doi: 10.1016/j.ymeth.2019.05.016. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.0 (unified framework):**
> Manuscript in preparation

---

## License

MIT License - See [LICENSE](LICENSE) for details

## Contact

- Kevin R. Roy - [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com)
- GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)

## Acknowledgments

- Original RECTIFY development supported by Chanfreau Lab, UCLA
- NET-seq data from Churchman Lab, Harvard Medical School
