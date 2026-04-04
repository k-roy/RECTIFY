# RECTIFY

<p align="center">
  <strong>R</strong>NA 5' and 3' <strong>E</strong>nd <strong>C</strong>orrection <strong>T</strong>ool with <strong>I</strong>ntron re<strong>F</strong>inement and ambiguit<strong>Y</strong> resolution
</p>

<p align="center">
  <a href="https://pypi.org/project/rectify-rna/"><img src="https://img.shields.io/pypi/v/rectify-rna?color=blue&label=PyPI" alt="PyPI"></a>
  <a href="https://rectify-rna.readthedocs.io"><img src="https://img.shields.io/readthedocs/rectify-rna?label=docs" alt="Docs"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.8%2B-blue.svg" alt="Python 3.8+"></a>
</p>

---

Precision transcript structure mapping for direct RNA nanopore sequencing. RECTIFY provides accurate 5' ends, 3' ends, and splice junctions through multi-aligner consensus, artifact correction, and optional NET-seq refinement.

## Key Features

| Feature | Description |
|:--------|:------------|
| **Multi-Aligner Consensus** | Runs minimap2, pbmm2, gapmm2 (optionally uLTRA + deSALT) and selects best junction set per read |
| **5' End Junction Recovery** | Rescues soft-clipped bases by extending alignments through annotated splice junctions |
| **3' End Indel Correction** | Walk-back algorithm fixes alignment artifacts where poly(A) tails land on genomic A-tracts |
| **False Junction Handling** | Removes spurious N-operations created by poly(A) tail alignment into downstream A-tracts |
| **Junction Ambiguity Resolution** | Resolves reads matching multiple junctions using proportional assignment |
| **Poly(A) Measurement** | Reports tail length (aligned + soft-clipped A's) per read |
| **NET-seq Refinement** | Optional NNLS deconvolution resolves remaining A-tract ambiguity using nascent RNA data |
| **Adaptive Clustering** | Groups CPA sites using valley-based algorithm; handles multi-site genes cleanly |
| **Dual-Resolution DESeq2** | Gene-level and cluster-level differential expression in one command |
| **APA Shift Analysis** | Detects proximal/distal CPA site usage changes across conditions |
| **Visualization** | Metagene plots, genome browser figures, heatmaps (`pip install rectify-rna[visualize]`) |

**Bundled data for yeast:** For *S. cerevisiae*, RECTIFY includes the S288C genome, SGD annotations, GO terms, WT NET-seq data, and 64K pre-computed A-tract CPA sites — no external files needed.

---

## Quick Start

```bash
pip install rectify-rna

# Single sample — bundled yeast genome (no external files needed)
rectify correct reads.fastq.gz --organism yeast -o corrected.tsv

# Full pipeline: correct + multi-sample analysis
rectify run reads.bam --genome genome.fa --annotation genes.gtf --output-dir results/
```

See the [Quick Start guide](quickstart.md) for a step-by-step walkthrough.

---

## How It Works

### 1. Multi-Aligner Consensus

Different aligners make different tradeoffs at splice junctions. RECTIFY runs up to five aligners in parallel, attempts to rescue soft-clips through known junctions, scores alignments by canonical GT-AG sites and annotation matches, and selects the best alignment per read.

### 2. 3' End Indel Correction

When poly(A) tails align to genomic A-tracts, aligners introduce indel artifacts that shift the apparent 3' end **downstream** of the true cleavage site. RECTIFY's walk-back algorithm starts at the soft-clip boundary and steps upstream, skipping A's, deletions, and T sequencing errors, stopping at the first non-A/T position where the read and genome agree.

### 3. Optional NET-seq Refinement

For organisms with NET-seq data (bundled for yeast), RECTIFY resolves remaining A-tract ambiguity using NNLS deconvolution with a point-spread function derived from 5,000+ zero-A control sites.

### 4. Adaptive Clustering and Differential Expression

After correction, CPA sites are clustered using a valley-based algorithm, and DESeq2 is run at both gene and cluster resolution. Cluster-level analysis detects isoform-specific changes that gene-level counts would miss.

---

## Supported Technologies

Nanopore direct RNA-seq · QuantSeq (oligo-dT short-read) · PacBio Iso-Seq · NET-seq · Any poly(A)-tailed RNA-seq

---

## Citation

> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020;176:4-13. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.0:** Manuscript in preparation.

---

## License

MIT — see [LICENSE](https://github.com/k-roy/RECTIFY/blob/main/LICENSE) for details.

## Contact

Kevin R. Roy · [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com) · GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)
