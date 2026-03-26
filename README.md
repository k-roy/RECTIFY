# RECTIFY: Precision Transcript Structure Mapping for Direct RNA Nanopore Sequencing

[![PyPI version](https://badge.fury.io/py/rectify-rna.svg)](https://pypi.org/project/rectify-rna/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**Accurate 5' ends, 3' ends, and splice junctions through artifact correction and optional NET-seq refinement.**

---

## Quick Start

```bash
pip install rectify-rna

# FASTQ input with bundled yeast genome (no external files needed)
rectify correct reads.fastq.gz --organism yeast -o corrected.tsv

# Full pipeline: correct + analyze
rectify run reads.bam --genome genome.fa --annotation genes.gtf --output-dir results/
```

---

## Key Features

| Feature | Description |
|---------|-------------|
| **A-tract Correction** | Detects genomic A-tracts and resolves positional ambiguity |
| **Indel Correction** | Fixes minimap2 alignment artifacts in poly(A) regions |
| **Poly(A) Measurement** | Reports tail length (aligned + soft-clipped) |
| **AG Mispriming Detection** | Flags internal priming on A/G-rich regions |
| **NET-seq Refinement** | Resolves ambiguity using nascent RNA data (optional) |
| **Adaptive Clustering** | Groups CPA sites with valley-based algorithm |
| **Dual-Resolution DESeq2** | Gene-level and cluster-level differential expression |
| **APA Shift Analysis** | Detects proximal/distal CPA site usage changes |
| **Visualization** | Metagene plots, genome browser figures (`pip install rectify-rna[visualize]`) |

### Bundled Data (Yeast)

For *S. cerevisiae*, RECTIFY includes the S288C genome, SGD annotations, GO terms, WT NET-seq data, and 64K pre-computed A-tract CPA sites—no external files needed.

---

## What RECTIFY Corrects

### The Problem

When poly(A) tails align to genomic A-tracts, aligners introduce indels that shift the apparent 3' end downstream:

```
Genome: ...CTAGTGACAGTC|AAAAAAAAAAA|CTAGCGATC...
                       |           |
                  true CPA    genomic A-tract

Read:   ...CTAGTGACAGTCAAAAAAAATAAA-AAAAA--AAAAAAAAAA...
                                ↑      ↑↑
                             T error  deletions (artifacts)
```

### The Solution

RECTIFY walks upstream from the soft-clip boundary, absorbing A's and deletions until finding the first non-A agreement between genome and read. This recovers the true CPA position.

For ambiguous cases, optional NET-seq refinement uses NNLS deconvolution to resolve multiple peaks and assign reads proportionally.

---

## Output

Each read gets a corrected position with confidence scores:

```
read_id   │ chrom │ strand │ original │ corrected │ shift │ confidence │ polya_len │ qc_flags
read001   │ chrI  │   +    │  147592  │   147585  │  -7   │    HIGH    │    42     │   PASS
read002   │ chrI  │   +    │  147594  │   147591  │  -3   │   MEDIUM   │    38     │   PASS
read003   │ chrII │   +    │  283109  │   283104  │  -5   │    LOW     │    31     │ AG_RICH
```

The `rectify analyze` command produces:
- `clusters.tsv` - CPA site clusters with read counts
- `deseq2_gene_results.tsv` - Gene-level differential expression
- `deseq2_cluster_results.tsv` - Cluster-level differential expression
- `shift_results.tsv` - Genes with APA shifts
- `go_enrichment.tsv` - GO enrichment for DE genes
- `motif_results/` - Enriched sequence motifs near CPA sites

---

## Installation

```bash
# PyPI
pip install rectify-rna

# With visualization support
pip install rectify-rna[visualize]

# Conda (includes MEME Suite for motif discovery)
conda install -c conda-forge -c bioconda rectify-rna
```

---

## Usage

```bash
# Basic correction
rectify correct reads.bam --genome genome.fa --output corrected.tsv

# With custom NET-seq data
rectify correct reads.bam --genome genome.fa --netseq-dir my_netseq/ --output corrected.tsv

# Analysis pipeline
rectify analyze corrected.tsv --annotation genes.gtf --output-dir results/

# Export bigWig tracks
rectify export corrected.tsv -o tracks/ --genome genome.fa
```

---

## Supported Technologies

- Nanopore direct RNA-seq (minimap2)
- QuantSeq (oligo-dT short-read)
- PacBio Iso-Seq
- Any poly(A)-tailed RNA-seq

---

## Citation

**Original RECTIFY:**
> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020;176:4-13. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.0:** Manuscript in preparation

---

## License

MIT License - See [LICENSE](LICENSE) for details

## Contact

- Kevin R. Roy - [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com)
- GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)
