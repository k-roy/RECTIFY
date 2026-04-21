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

---

## Overview

Nanopore direct RNA sequencing offers unprecedented read lengths, but accurate transcript structure mapping requires solving five intertwined problems: poly(A) tails that mislead aligners into spurious splice junctions and indels (requiring pre-trimming for DRS data), residual 3' end artifacts in A-tract regions, imprecise 5' ends and junction boundaries from soft-clipping, homopolymer-driven 3' soft-clipping at CPA sites, and conflicting junction calls between aligners. **RECTIFY** solves all five through pre-alignment poly(A) removal, multi-aligner rectification, artifact-aware corrections, sequence-evidence-first junction refinement, and optional NET-seq refinement — delivering nucleotide-precision 5' and 3' end coordinates and splice junction sets.

**Use RECTIFY when you need:**
- Accurate cleavage and polyadenylation (CPA) site mapping from DRS data
- Correction of poly(A) misalignment artifacts in A-tract regions
- Robust splice junction calls from reads spanning multiple exons
- Detection of alternative polyadenylation (APA) with cluster-level resolution
- Differential expression analysis at gene and isoform levels
- Optional NET-seq-informed refinement for A-tract ambiguity

---

## Quick Start

### Installation

```bash
# Via PyPI
pip install rectify-rna

# With visualization support (metagene plots, genome figures)
pip install rectify-rna[visualize]

# Via Conda (includes MEME Suite for motif discovery)
conda install -c conda-forge -c bioconda rectify-rna
```

### Basic Usage

```bash
# Correct 3' ends from FASTQ (bundled yeast genome — no external files needed)
rectify correct reads.fastq.gz --organism yeast -o corrected.tsv

# Full pipeline: alignment → correction → analysis
rectify run-all reads.bam --genome genome.fa --annotation genes.gtf --output-dir results/

# Process NET-seq data (nascent RNA 3' ends)
rectify netseq netseq.bam --genome genome.fa --gff genes.gff -o netseq_output/
```

---

## How It Works

RECTIFY reconstructs true RNA 3' and 5' ends through a pre-alignment trimming step and four sequential corrections, each addressing a specific alignment artifact.

### 0. DRS Pre-Trimming: Removing Poly(A) Tails Before Alignment

For direct RNA sequencing (DRS), poly(A) tails are present in the read sequence. When left untrimmed, they mislead aligners into creating spurious splice junctions and indels as aligners try to force-align A-rich tail bases against genomic A-tracts — the very artifacts that downstream correction modules must then repair. **RECTIFY eliminates the source** with `rectify trim-polya`, which strips the poly(A) tail and adapter stub before re-alignment, producing clean reads that align correctly.

<p align="center">
  <img src="docs/figures/polya_pretrim.png" alt="DRS Poly(A) Pre-Trimming" width="680">
</p>

The trimmer uses a three-pass algorithm in RNA 5'→3' orientation:

- **Pass 0** — Adapter stub removal: regex `T[CT]{0,10}$` strips the 3'-terminal adapter stub in a single scan.
- **Pass 1** — Pure-A tail scan: slides leftward from the stripped boundary counting consecutive A-rich bases (strict mode: 0% non-A allowed) to locate the poly(A) / transcript body junction.
- **Pass 2** — Iterative peel: handles ambiguous boundaries where Nanopore basecalling errors (T calls within the tail) confuse the pure scan; peels the tail one base at a time until A-richness drops below threshold.

Reads with no detectable tail pass through unchanged. Output: an **unaligned BAM** with poly(A)-free reads for re-alignment, plus a **per-read metadata parquet** recording tail length, adapter sequence, and pass number — used by `rectify restore-softclip` (Step 4) to re-attach the original tail to the IGV softclip BAM for tail-length visualization.

### 1. 3' End Walk-Back: Recovering the True CPA Site

When poly(A) tails align to genomic A-tracts, aligners introduce indels and spurious splice junctions (N operations) to maximize alignment score, **shifting the apparent 3' end far downstream** of the true cleavage site. RECTIFY walks backward from the soft-clip boundary, skipping A's, deletions, T sequencing errors, and any intron-skip (N) operations it encounters, until it finds the first non-A/T agreement between genome and read — the true CPA site.

<p align="center">
  <img src="docs/figures/indel_correction.png" alt="3' End Walk-Back Correction" width="680">
</p>

**Why simple poly(A) trimming fails:** The boundary between genomic A's and tail A's is ambiguous in A-tract regions. RECTIFY's walk-back algorithm handles deletions, T sequencing errors, and false splice junctions within the A-tract, recovering the true CPA position even when the aligner has spread the poly(A) signal across multiple genomic A-runs or introduced spurious N operations to reach downstream A-tracts. For minus-strand genes, the poly(A) tail appears as a poly(T) prefix extending leftward — RECTIFY applies identical logic in reverse orientation.

**False junction cleanup is built-in:** Poly(A) tails can cause aligners to introduce skip (N) operations to reach downstream A-tracts, creating spurious splice junctions. The same walk-back that corrects indel artifacts transparently absorbs these N operations — they require no separate detection step.

<p align="center">
  <img src="docs/figures/false_junction_walkback.png" alt="False Junction Walk-Back" width="680">
</p>

### 2. Unified 5' End and Junction Correction: Rescue and N-Op Refinement

RECTIFY applies two complementary modules that together deliver nucleotide-precision 5' ends and splice junction boundaries.

**Part 1 — 5' Soft-Clip Junction Rescue (Cat3):** Nanopore reads that begin near a splice junction frequently have their 5'-most bases soft-clipped rather than placed in the upstream exon, because the short exon fragment is too ambiguous for the aligner to place confidently. RECTIFY identifies these soft-clipped sequences, locates the nearest annotated donor site, performs a semi-global Needleman-Wunsch alignment (affine gap, Gotoh 1982) between the soft-clip bases and the upstream exon reference, and extends the alignment through the intron with a proper M/I/D CIGAR — recovering the true transcription start position with per-base resolution.

**Part 2 — Post-Consensus N-Op Junction Refinement (Module 2H):** After consensus aligner selection, every N-op (splice junction) in every read is re-scored. For each N-op, RECTIFY collects all candidate junctions within a search radius and scores each with a homopolymer-aware semi-global alignment: the rescue sequence (bases downstream of the current N-op split) is aligned left-anchored to the candidate intron end, with HP-aware linear gap costs (HP deletions cost 0.5, non-HP deletions 1.0, insertions 1.25). The winner is selected by strict priority: sequence match score first, then current-junction stability (equal-scoring candidates that match the existing N-op are never displaced), then canonical GT-AG motif, then annotation status, then boundary shift distance — **sequence evidence always overrides annotation**.

<p align="center">
  <img src="docs/figures/5prime_junction_rescue.png" alt="Unified 5' End and Junction Correction" width="680">
</p>

A fast path skips scoring for reads already at an annotated canonical-tier-0 junction, providing a 255× speedup. The `max_boundary_shift` parameter (default 50 bp) prevents false matches from junctions in neighbouring genes; `search_radius` (default 5000 bp) controls candidate discovery.

### 3. 3' End Soft-Clip Rescue at Homopolymer Boundaries

Nanopore basecallers systematically under-call homopolymer runs. At cleavage and polyadenylation (CPA) sites with upstream T-rich genomic regions, this causes the aligner to terminate the alignment prematurely, soft-clipping non-T bases that actually belong to the transcript body. RECTIFY identifies soft-clipped sequences, skips remaining reference homopolymer bases, and matches them to downstream reference positions.

<p align="center">
  <img src="docs/figures/softclip_rescue.png" alt="Soft-Clip Rescue at Homopolymer Boundaries" width="680">
</p>

This correction is especially critical for detecting true 3' ends in regions where weak basecalling and homopolymer under-calling create false soft-clip boundaries.

### 4. Multi-Aligner Rectification: Selecting the Optimal Junction Set

Different aligners make different tradeoffs at splice junctions. RECTIFY solves this in three stages:

**Stage 1 — Per-aligner rectification:** `rectify correct` is applied independently to each aligner's BAM (minimap2, mapPacBio, gapmm2). Every correction module (3' walk-back, 5' junction rescue, soft-clip rescue, false-junction filter) runs on each aligner's output, producing a separate corrected TSV per aligner.

**Stage 2 — Consensus selection:** The per-aligner corrected TSVs are merged by `rectify consensus`, which selects the winning aligner per read using post-rectification features (in priority order): (1) `five_prime_rescued` — prefer the aligner where Cat3 5' rescue fired; (2) `confidence` — high > medium > low; (3) corrected_3′ agreement — prefer positions agreed on by most aligners; (4) alignment span — prefer wider reference span; (5) `n_junctions` — prefer more completely spliced.

**Stage 3 — Chimeric reconstruction:** For reads where two or more aligners each uniquely contribute a junction not present in the other's corrected output, `rectify consensus` can optionally stitch the complementary junctions into a single chimeric alignment — recovering reads that no single aligner handles completely.

<p align="center">
  <img src="docs/figures/multi_aligner_consensus.png" alt="Multi-Aligner Rectification Pipeline" width="680">
</p>

```bash
# Align and rectify with all three aligners (default, DRS-optimized)
rectify align reads.fastq.gz --genome genome.fa --annotation genes.gff -o aligned_dir/

# Merge per-aligner corrected TSVs into a single consensus result
rectify consensus minimap2:aligned_dir/minimap2/corrected_3ends.tsv \
                  mapPacBio:aligned_dir/mapPacBio/corrected_3ends.tsv \
                  gapmm2:aligned_dir/gapmm2/corrected_3ends.tsv \
                  -o corrected_3ends.tsv

# Single-aligner mode (faster, less accurate)
rectify align reads.fastq.gz --genome genome.fa --aligner minimap2 -o aligned.bam
```

---

## Key Features

| Feature | Benefit |
|:--------|:--------|
| **DRS Poly(A) Pre-Trimming** | Three-pass algorithm (adapter stub strip → pure-A scan → iterative peel) removes poly(A) + adapter before re-alignment; tail lengths and adapter sequences stored in metadata parquet for downstream restoration |
| **Multi-Aligner Rectification** | Rectifies each aligner independently, then selects the winning aligner per read using post-rectification features (5' rescue > confidence > agreement > span > junctions); optionally stitches complementary junctions from two aligners (chimeric reconstruction) |
| **Unified 5' Rescue and Junction Refinement** | Cat3 rescues 5'-soft-clipped reads via semi-global NW alignment to the upstream exon; Module 2H refines every N-op boundary post-consensus using homopolymer-aware split-alignment where sequence evidence always overrides annotation |
| **3' End Walk-Back** | Walks backward from soft-clip boundary to recover true CPA site, transparently absorbing indels, T sequencing errors, and spurious splice junctions (N ops) in a single pass |
| **Junction Ambiguity Resolution** | Resolves reads matching multiple junctions using proportional assignment |
| **Poly(A) Measurement** | Reports tail length including both aligned and soft-clipped bases |
| **NET-seq Refinement** | Uses nascent RNA 3' ends to deconvolve A-tract ambiguity (optional) |
| **Adaptive Clustering** | Groups nearby CPA sites using valley-based peak detection |
| **Dual-Resolution Differential Expression** | DESeq2 at both gene level and cluster (isoform) level |
| **APA Shift Analysis** | Detects significant proximal/distal CPA site usage changes |
| **Visualization** | Metagene plots and genome browser figures (`pip install rectify-rna[visualize]`) |
| **Bundled Yeast Data** | S288C genome, SGD annotations, GO terms, WT NET-seq, 64K pre-computed A-tract CPA sites |

---

## Output and Results

Each read receives a corrected position with confidence scoring:

```
read_id   │ chrom │ strand │ original │ corrected │ shift │ confidence │ polya_len │ qc_flags
read001   │ chrI  │   +    │  147592  │   147585  │  -7   │    HIGH    │    42     │   PASS
read002   │ chrI  │   +    │  147594  │   147591  │  -3   │   MEDIUM   │    38     │   PASS
read003   │ chrII │   +    │  283109  │   283104  │  -5   │    LOW     │    31     │ AG_RICH
```

The `rectify analyze` command produces:
- **clusters.tsv** — CPA site clusters with read counts per condition
- **deseq2_gene_results.tsv** — Differential expression at gene level
- **deseq2_cluster_results.tsv** — Differential expression at cluster (isoform) level
- **shift_results.tsv** — Genes with statistically significant APA shifts
- **go_enrichment.tsv** — GO term enrichment on shifted genes
- **motif_results/** — Enriched sequence motifs near CPA sites

---

## NET-seq Refinement (Optional)

For organisms with nascent RNA (NET-seq) data, RECTIFY resolves remaining ambiguity within A-tracts. NET-seq samples RNA still attached to polymerase, providing a reference for true CPA positions. Since nascent RNA is oligo-adenylated post-capture, RECTIFY uses NNLS deconvolution with a point-spread function derived from 5000+ zero-A calibration sites to recover true CPA positions.

<p align="center">
  <img src="docs/figures/oligo_a_spreading.png" alt="Oligo(A) Spreading Artifact" width="500">
</p>

<p align="center">
  <img src="docs/figures/oligo_a_deconvolution.png" alt="Oligo(A) Deconvolution" width="680">
</p>

**For *S. cerevisiae***, bundled WT NET-seq data is auto-detected. For other organisms or mutant conditions, provide NET-seq bigWigs with the `--netseq-dir` flag.

---

## Commands Reference

**Pipeline steps (DRS direct RNA-seq):**

| Step | Command | Purpose |
|:-----|:--------|:--------|
| 0 *(DRS only)* | `rectify trim-polya` | Trim poly(A) tail + adapter stub from Dorado-aligned BAM; writes unaligned BAM + per-read metadata parquet |
| 1 | `rectify align` | Align FASTQ with multiple aligners (minimap2, mapPacBio, gapmm2, deSALT, uLTRA) |
| 2 | `rectify consensus` | Select best aligner per read; tags: XA (winner), XC (confidence), XN (agreement count) |
| 3 | `rectify correct` | Correct 3' ends (indel correction + A-tract resolution + junction rescue); writes cp:i: tag |
| 4 *(DRS only)* | `rectify restore-softclip` | Re-attach trimmed poly(A)+adapter to softclip BAM for IGV visualization of tail lengths |
| 5 | `rectify analyze` | Downstream analysis (CPA clustering, DESeq2, GO enrichment, motif discovery) |

For oligo-dT cDNA sequencing, skip Steps 0 and 4 — start at Step 1 directly from your FASTQ.

**Other commands:**

| Command | Purpose |
|:--------|:--------|
| `rectify export` | Export corrected positions to bigWig/bedGraph tracks |
| `rectify extract` | Extract per-read 5'/3' ends and junctions to TSV |
| `rectify aggregate` | Group reads into CPA / TSS / junction datasets |
| `rectify netseq` | Process NET-seq BAM files (3' extraction + deconvolution) |
| `rectify run-all` | Full pipeline (Steps 1–5) with provenance tracking and step-skip |

<details>
<summary><b>Usage examples</b></summary>

```bash
# Correct 3' ends (bundled yeast genome, no external files needed)
rectify correct reads.fastq.gz --organism yeast -o corrected.tsv

# Correct with custom genome and optional NET-seq deconvolution
rectify correct reads.bam --genome genome.fa --netseq-dir my_netseq/ -o corrected.tsv

# Extract per-read features (5'/3' ends, junctions) to TSV
rectify extract reads.bam -o reads.tsv --genome genome.fa --annotation genes.gff

# Aggregate into separate 3'/5'/junction datasets by condition
rectify aggregate reads.bam -o aggregated/ --annotation genes.gff --mode all

# Differential expression analysis (gene and cluster level)
rectify analyze corrected.tsv --annotation genes.gtf --output-dir results/

# Export corrected positions as genome browser tracks
rectify export corrected.tsv -o tracks/ --genome genome.fa

# Complete pipeline from reads to differential expression
rectify run-all reads.bam --genome genome.fa --annotation genes.gtf --output-dir results/

# Process NET-seq data (nascent RNA 3' ends for A-tract refinement)
rectify netseq netseq.bam --genome genome.fa --gff genes.gff -o netseq_output/
```

</details>

---

## Supported Technologies

**Direct RNA sequencing:** Nanopore direct RNA-seq (DRS)
**Short-read quantification:** QuantSeq (oligo-dT), PacBio Iso-Seq, NET-seq
**General:** Any poly(A)-tailed RNA-seq platform

---

## Citation

Please cite RECTIFY if you use it in your research:

> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020;176:4-13. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.0:** Manuscript in preparation.

---

## License

MIT — see [LICENSE](LICENSE) for details.

## Contact

**Kevin R. Roy**
Email: [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com)
GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)
