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
| **3' End Indel Correction** | Fixes alignment artifacts where poly(A) tails align to genomic A-tracts |
| **5' End Junction Correction** | Recovers true splice sites from soft-clipped junction reads |
| **Poly(A) Measurement** | Reports tail length (aligned + soft-clipped) |
| **AG Mispriming Detection** | Flags internal priming on A/G-rich regions |
| **NET-seq Refinement** | Resolves A-tract ambiguity using nascent RNA data (optional) |
| **Adaptive Clustering** | Groups CPA sites with valley-based algorithm |
| **Dual-Resolution DESeq2** | Gene-level and cluster-level differential expression |
| **APA Shift Analysis** | Detects proximal/distal CPA site usage changes |
| **Visualization** | Metagene plots, genome browser figures (`pip install rectify-rna[visualize]`) |

### Bundled Data (Yeast)

For *S. cerevisiae*, RECTIFY includes the S288C genome, SGD annotations, GO terms, WT NET-seq data, and 64K pre-computed A-tract CPA sites—no external files needed.

---

## 3' End Correction: Indel Artifacts in Poly(A) Regions

When poly(A) tails align to genomic A-tracts, aligners like minimap2 introduce indel artifacts to maximize alignment score. This shifts the apparent 3' end **downstream** of the true cleavage site.

```
===============================================================================
THE PROBLEM: Poly(A) tail alignment artifacts
===============================================================================

                                              true CPA site
                                                   ↓
Genome: 5'..CTAGTGACAGTCAAAAAAAA-AAACAAAAGTAAAAAAAAAAAA|CTAGCGATC..3'
                                                       |
                                              genomic A-tract ends here

Read:   5'..CTAGTGACAGTCAAAAAAAATAAA-AAAAA--AAAAAAAAAA.|AAAAAAAAAAA
                                 ↑      ↑↑            |<--------->
                              T error  deletions    soft-clipped tail
                              (seq err) (artifacts)

The aligner introduces deletions to extend alignment of poly(A) tail bases
into the genomic A-tract. The apparent 3' end shifts downstream.

===============================================================================
RECTIFY'S SOLUTION: Walk upstream to find true CPA
===============================================================================

Starting from the soft-clip boundary, walk UPSTREAM through the aligned region:

  1. Skip positions where genome = A (ambiguous with poly(A) tail)
  2. Skip deletions (D in CIGAR) - these are alignment artifacts
  3. Skip T sequencing errors in the tail
  4. STOP at first non-A/T agreement between genome and read

                                          soft-clip boundary
                                                   ↓
Genome: ...CTAGTGACAGTCAAAAAAAA-AAACAAAAGTAAAAAAAAAAAA|CTAGCGATC...
Read:   ...CTAGTGACAGTCAAAAAAAATAAA-AAAAA--AAAAAAAAAA.|AAAAAAAAAAA
                     ↑         ↑      ↑↑             |
                  STOP here  T err   dels         soft-clip
                  (C = C)   (skip)  (absorb)

Result: True 3' end at the C position
        Poly(A) length = soft-clipped + aligned A's + absorbed deletions
```

**Key insight for IGV users:** For minus strand reads, the poly(A) tail appears as poly(T) extending leftward. RECTIFY corrects by shifting rightward toward the true CPA.

---

## 5' End Correction: Splice Junction Soft-Clips

Long reads spanning splice junctions often have soft-clipped bases at the 5' end where the aligner fails to find the exact junction boundary. RECTIFY recovers the true splice site.

```
===============================================================================
THE PROBLEM: Soft-clipped bases hide the true 5' junction
===============================================================================

                    true 5' splice site
                           ↓
Genome:     ...EXON1|gt----intron----ag|EXON2...
                    ↑                  ↑
                  donor              acceptor

Read:       NNNNNN|====================|===================>
            ↑     ↑
         soft-    aligned portion starts here
         clipped  (aligner missed the exact junction)

The soft-clipped bases (N's) actually match EXON1, but the aligner
couldn't extend through the splice junction.

===============================================================================
RECTIFY'S SOLUTION: Extend alignment to known junction
===============================================================================

  1. Identify reads with 5' soft-clips near annotated splice sites
  2. Check if soft-clipped sequence matches upstream exon
  3. Extend the alignment to the canonical splice donor (GT)

Before:  NNNNNN|====================...
              ↑
           read starts here (wrong)

After:   |=========================...
         ↑
      read starts at true exon boundary

Result: Accurate 5' end for TSS analysis and full-length read classification
```

---

## Adaptive Clustering and Differential Expression

After correction, RECTIFY groups nearby CPA sites into clusters using a valley-based algorithm, then runs DESeq2 at both gene and cluster resolution.

```
===============================================================================
ADAPTIVE CLUSTERING: Valley-based CPA site grouping
===============================================================================

Signal:
            ▲                    ▲▲
           ███                  ████
          █████      ▲         ██████
         ███████    ███       ████████
        █████████  █████     ██████████
       ─────┴─────┴─────────┴──────────────> position
            │  ↑  │    ↑    │     ↑
         cluster1  valley  cluster2  cluster3

Algorithm:
  1. Find peaks (local maxima in 3' end signal)
  2. Find valleys (local minima between peaks)
  3. Set boundaries at midpoint between peak and valley (capped at ±10bp)

Why cluster-level analysis matters:
  - Genes often have MULTIPLE CPA sites (alternative polyadenylation)
  - Conditions may shift usage between proximal/distal sites
  - Cluster-level DESeq2 detects isoform-specific changes that gene-level misses
```

**Dual-resolution output:**

| Level | Detects | Example |
|-------|---------|---------|
| **Gene** | Total expression changes | HSP82 is 2-fold down in heat shock |
| **Cluster** | CPA site usage changes | FAS1 shifts from distal to proximal site |

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
- `clusters.tsv` - CPA site clusters with read counts per sample
- `deseq2_gene_results.tsv` - Gene-level differential expression
- `deseq2_cluster_results.tsv` - Cluster-level differential expression
- `shift_results.tsv` - Genes with significant APA shifts
- `go_enrichment.tsv` - GO enrichment for DE genes
- `motif_results/` - Enriched sequence motifs near CPA sites

---

## NET-seq Refinement (Optional)

For organisms with NET-seq data, RECTIFY can resolve remaining ambiguity within A-tracts. Nascent RNA 3' ends from NET-seq are oligo-adenylated, creating a characteristic spreading pattern. RECTIFY uses NNLS deconvolution to recover true CPA positions and assign reads proportionally to multiple peaks.

Bundled WT NET-seq data for yeast is auto-detected. For other organisms or mutant conditions, provide NET-seq bigWigs with `--netseq-dir`.

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
