<p align="center">
  <h1 align="center">RECTIFY</h1>
  <p align="center">
    <b>R</b>NA 5' and 3' <b>E</b>nd <b>C</b>orrection <b>T</b>ool with <b>I</b>ntron re<b>F</b>inement and ambiguit<b>Y</b> resolution
  </p>
  <p align="center">
    <a href="https://pypi.org/project/rectify-rna/"><img src="https://img.shields.io/pypi/v/rectify-rna?color=blue&label=PyPI" alt="PyPI"></a>
    <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT"></a>
    <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.8%2B-blue.svg" alt="Python 3.8+"></a>
  </p>
</p>

Precision transcript structure mapping for direct RNA nanopore sequencing. RECTIFY provides accurate 5' ends, 3' ends, and splice junctions through multi-aligner consensus, artifact correction, and optional NET-seq refinement.

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
|:--------|:------------|
| **Multi-Aligner Consensus** | Runs minimap2, pbmm2, gapmm2 and selects best junction set per read |
| **5' End Junction Recovery** | Rescues soft-clipped bases by extending alignments through splice junctions |
| **3' End Indel Correction** | Fixes alignment artifacts where poly(A) tails align to genomic A-tracts |
| **3' False Junction Handling** | Walk back correction removes spurious junctions from poly(A) artifacts |
| **Junction Ambiguity Resolution** | Resolves reads matching multiple junctions using proportional assignment |
| **Poly(A) Measurement** | Reports tail length (aligned + soft-clipped) |
| **NET-seq Refinement** | Resolves A-tract ambiguity using nascent RNA data (optional) |
| **Adaptive Clustering** | Groups CPA sites with valley-based algorithm |
| **Dual-Resolution DESeq2** | Gene-level and cluster-level differential expression |
| **APA Shift Analysis** | Detects proximal/distal CPA site usage changes |
| **Visualization** | Metagene plots, genome browser figures (`pip install rectify-rna[visualize]`) |

**Bundled data for yeast:** For *S. cerevisiae*, RECTIFY includes the S288C genome, SGD annotations, GO terms, WT NET-seq data, and 64K pre-computed A-tract CPA sites — no external files needed.

---

## How It Works

### 3' End Indel Correction

When poly(A) tails align to genomic A-tracts, aligners introduce indel artifacts to maximize alignment score, shifting the apparent 3' end **downstream** of the true cleavage site. RECTIFY walks upstream from the soft-clip boundary, skipping A's, deletions, and T sequencing errors, stopping at the first non-A/T agreement between genome and read.

<p align="center">
  <img src="docs/figures/indel_correction.png" alt="3' End Indel Correction" width="680">
</p>

**Why not just trim poly(A)?** Simple poly(A) trimming fails at A-tracts because the boundary between genomic A's and tail A's is ambiguous. RECTIFY's walk-back algorithm handles deletions and T sequencing errors within the A-tract, finding the true CPA site even when the aligner has spread the tail across multiple genomic A-runs. For minus-strand genes, the poly(A) tail appears as a poly(T) prefix extending leftward — RECTIFY applies the same walk-back logic in the opposite direction.

### Soft-Clip Rescue at Homopolymer Boundaries

Nanopore basecallers systematically under-call homopolymer runs. When this happens at CPA sites with upstream T-tracts, the aligner soft-clips non-T bases instead of placing them correctly. RECTIFY skips over remaining reference homopolymer bases and matches the soft-clipped sequence to the downstream reference.

<p align="center">
  <img src="docs/figures/softclip_rescue.png" alt="Soft-Clip Rescue" width="620">
</p>

### 3' False Junction Handling

Poly(A) tails can create spurious "junctions" when the aligner introduces a skip (N) operation to align tail bases to a downstream A-tract. The walk back algorithm handles this automatically — it eats through all aligned A's and discards any N operations it encounters, finding the true CPA site without needing special false junction detection.

<p align="center">
  <img src="docs/figures/false_junction_walkback.png" alt="False Junction Walk Back" width="660">
</p>

### 5' End Correction

Long reads spanning splice junctions often have soft-clipped bases at the 5' end. These bases actually match the upstream exon, but the aligner couldn't extend through the junction. RECTIFY identifies these reads, checks for upstream exon matches, and extends the alignment to the canonical splice donor.

<p align="center">
  <img src="docs/figures/5prime_junction_rescue.png" alt="5' Junction Rescue" width="660">
</p>

> **Note:** Due to 5'-to-3' degradation in direct RNA sequencing, the read's 5' end often does not represent the true TSS.

### Multi-Aligner Consensus

Different aligners make different tradeoffs at splice junctions. RECTIFY runs three aligners in parallel, attempts to rescue soft-clips through known junctions, scores the resulting alignments by canonical splice sites and annotation matches, and selects the best alignment per read.

<p align="center">
  <img src="docs/figures/multi_aligner_consensus.png" alt="Multi-Aligner Consensus" width="660">
</p>

**Scoring criteria:** Rescued alignments are scored by the number of canonical GT-AG splice junctions, matches to annotated junctions in the provided GFF/GTF, and remaining soft-clip length. The alignment with the highest composite score is selected per read and written to the rectified BAM.

```bash
# Multi-aligner consensus alignment (default)
rectify align reads.fastq.gz --genome genome.fa --annotation genes.gff -o aligned.bam

# Single aligner mode (faster, less accurate)
rectify align reads.fastq.gz --genome genome.fa --aligner minimap2 -o aligned.bam
```

### Adaptive Clustering and Differential Expression

After correction, RECTIFY groups nearby CPA sites into clusters using a valley-based algorithm (find peaks → find valleys → set boundaries), then runs DESeq2 at both gene and cluster resolution. Cluster-level analysis detects isoform-specific changes that gene-level analysis misses.

<p align="center">
  <img src="docs/figures/adaptive_clustering.png" alt="Adaptive Clustering" width="600">
</p>

| Level | Detects | Example |
|:------|:--------|:--------|
| **Gene** | Total expression changes | HSP82 is 2-fold down in heat shock |
| **Cluster** | CPA site usage changes | FAS1 shifts from distal to proximal site |

---

## NET-seq Refinement

For organisms with NET-seq data, RECTIFY resolves remaining ambiguity within A-tracts. Nascent RNA 3' ends from NET-seq are oligo-adenylated, creating downstream signal spreading. RECTIFY uses NNLS deconvolution with a point-spread function derived from 5000+ zero-A sites to recover true CPA positions.

<p align="center">
  <img src="docs/figures/oligo_a_spreading.png" alt="Oligo(A) Spreading Artifact" width="500">
</p>

<p align="center">
  <img src="docs/figures/oligo_a_deconvolution.png" alt="Oligo(A) Deconvolution" width="680">
</p>

Bundled WT NET-seq data for yeast is auto-detected. For other organisms or mutant conditions, provide NET-seq bigWigs with `--netseq-dir`.

```bash
# Basic NET-seq processing
rectify netseq input.bam --genome genome.fa --gff genes.gff -o output/

# With exclusion region control
rectify netseq input.bam --genome genome.fa --gff genes.gff \
    --include-rdna --include-pol3 --exclude-mito -o output/
```

---

## Output

Each read gets a corrected position with confidence scores:

```
read_id   │ chrom │ strand │ original │ corrected │ shift │ confidence │ polya_len │ qc_flags
read001   │ chrI  │   +    │  147592  │   147585  │  -7   │    HIGH    │    42     │   PASS
read002   │ chrI  │   +    │  147594  │   147591  │  -3   │   MEDIUM   │    38     │   PASS
read003   │ chrII │   +    │  283109  │   283104  │  -5   │    LOW     │    31     │ AG_RICH
```

The `rectify analyze` command produces: `clusters.tsv` (CPA site clusters with read counts), `deseq2_gene_results.tsv` and `deseq2_cluster_results.tsv` (differential expression at both resolutions), `shift_results.tsv` (genes with significant APA shifts), `go_enrichment.tsv` (GO enrichment), and `motif_results/` (enriched sequence motifs near CPA sites).

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

## Commands

| Command | Description |
|:--------|:------------|
| `rectify correct` | Correct 3' end positions (indel correction, A-tract resolution) |
| `rectify analyze` | Downstream analysis (clustering, DESeq2, GO, motifs) |
| `rectify export` | Export corrected positions to bigWig/bedGraph tracks |
| `rectify extract` | Extract per-read info from BAM to TSV (5'/3' ends, junctions) |
| `rectify aggregate` | Aggregate reads into 3' end, 5' end, and junction datasets |
| `rectify align` | Align FASTQ with multi-aligner consensus (DRS-optimized) |
| `rectify netseq` | Process NET-seq BAM files (3' end extraction, deconvolution) |
| `rectify run` | Full pipeline: align (if FASTQ) → correct → analyze |

<details>
<summary><b>Usage examples</b></summary>

```bash
# Correct 3' ends (bundled yeast genome)
rectify correct reads.fastq.gz --organism yeast -o corrected.tsv

# Correct with custom genome and NET-seq
rectify correct reads.bam --genome genome.fa --netseq-dir my_netseq/ -o corrected.tsv

# Extract per-read features from BAM
rectify extract reads.bam -o reads.tsv --genome genome.fa --annotation genes.gff

# Aggregate into 3'/5'/junction datasets
rectify aggregate reads.bam -o aggregated/ --annotation genes.gff --mode all

# Differential expression analysis
rectify analyze corrected.tsv --annotation genes.gtf --output-dir results/

# Export bigWig tracks
rectify export corrected.tsv -o tracks/ --genome genome.fa

# Full pipeline
rectify run reads.bam --genome genome.fa --annotation genes.gtf --output-dir results/

# Process NET-seq data
rectify netseq netseq.bam --genome genome.fa --gff genes.gff -o netseq_output/
```

</details>

---

## Supported Technologies

Nanopore direct RNA-seq · QuantSeq (oligo-dT short-read) · PacBio Iso-Seq · NET-seq · Any poly(A)-tailed RNA-seq

---

## Citation

> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020;176:4-13. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.0:** Manuscript in preparation.

---

## License

MIT — see [LICENSE](LICENSE) for details.

## Contact

Kevin R. Roy — [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com) · GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)
