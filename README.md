# RECTIFY

**R**NA **E**nd **C**orrection **T**ool for **I**ntron re**F**inement with ambiguit**Y** resolution

[![PyPI version](https://badge.fury.io/py/rectify-rna.svg)](https://pypi.org/project/rectify-rna/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**Precision transcript structure mapping for direct RNA nanopore sequencing.** Accurate 5' ends, 3' ends, and splice junctions through multi-aligner consensus, artifact correction, and optional NET-seq refinement.

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
| **Multi-Aligner Consensus** | Runs minimap2, mapPacBio, gapmm2 and selects best junction set per read |
| **5' End Junction Recovery** | Rescues soft-clipped bases by extending alignments through splice junctions |
| **3' End Indel Correction** | Fixes alignment artifacts where poly(A) tails align to genomic A-tracts |
| **3' False Junction Handling** | Walk back correction eats through spurious junctions from poly(A) artifacts |
| **Junction Ambiguity Resolution** | Resolves reads matching multiple junctions using proportional assignment |
| **Poly(A) Measurement** | Reports tail length (aligned + soft-clipped) |
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

## Multi-Aligner Consensus Pipeline

Different aligners make different tradeoffs at splice junctions. RECTIFY runs three aligners in parallel and selects the best alignment per read.

```
===============================================================================
THE PROBLEM: Aligner disagreement at splice junctions
===============================================================================

minimap2:  ...EXON1|=====N=====|EXON2...  (spliced, but 5' soft-clipped)
              SSSSS↑

mapPacBio: ...EXON1|=====N=====|EXON2...  (spliced, no soft-clip)
                   ↑

gapmm2:    ...EXON1|=====N=====|EXON2...  (spliced, terminal exon refined)
                   ↑

Each aligner may find the correct junction, miss it, or soft-clip at the boundary.

===============================================================================
RECTIFY'S SOLUTION: Multi-aligner consensus with 5' rescue
===============================================================================

  1. Run all 3 aligners on the same reads
  2. For each read, compare junction sets across aligners
  3. Prefer alignments that:
     a. Splice through known junctions rather than soft-clipping (5' rescue)
     b. Use canonical splice sites (GT-AG)
     c. Are supported by multiple aligners (high confidence)
  4. Output: Single consensus BAM with best alignments

Note: 3' false junctions from poly(A) artifacts are handled separately by
walk back correction (see "3' False Junction Handling" below).

Confidence scoring:
  - HIGH:   3/3 aligners agree on junction
  - MEDIUM: 2/3 aligners agree
  - LOW:    1/3 aligners (only one found it)
```

**Usage:**

```bash
# Multi-aligner consensus alignment (default)
rectify align reads.fastq.gz --genome genome.fa --annotation genes.gff -o aligned.bam

# Single aligner mode (faster, less accurate)
rectify align reads.fastq.gz --genome genome.fa --aligner minimap2 -o aligned.bam
```

---

## 3' False Junction Handling

Poly(A) tails can create spurious "junctions" when the aligner introduces a skip (N) operation to align tail bases to a downstream genomic A-tract. **RECTIFY's walk back correction completely handles this artifact.**

```
===============================================================================
THE PROBLEM: Poly(A) tails create false junctions
===============================================================================

                              true CPA    false "junction"
                                  ↓              ↓
Genome:  ...EXON|AAAAA|----N----|AAAAA|...
Read:    ...EXON|AAAAAAAAAAAAAAAAAAAAA.|AAAAA (poly(A) tail)
                    ↑                        ↑
              aligned A's              soft-clipped tail

The aligner introduces an N (skip) to extend the alignment of poly(A) tail
bases into a downstream genomic A-tract, creating a spurious junction.

===============================================================================
RECTIFY'S SOLUTION: Walk back eats through false junctions
===============================================================================

The walk back algorithm finds the true 3' end by walking upstream through
ALL aligned A's until it finds the first non-A agreement between genome
and read. Crucially, it DISCARDS any N (skip) operations it encounters:

         walk back eats through everything
         <─────────────────────────────────
                                          ↓
Genome:  ...EXON|AAAAA|----N----|AAAAA|...
Read:    ...EXON|AAAAAAAAAAAAAAAAAAAAA.|AAAAA
              ↑                             ↑
           STOP here                  walk back starts
           (non-A match)              (eats A's, discards N)

Result:
  - Walk back finds the true CPA at the EXON/A boundary
  - The false junction (N operation) is completely ignored
  - No special false junction detection needed

This means 3' false junctions are NOT a problem for downstream analysis—
walk back correction handles them automatically as part of finding the
true cleavage and polyadenylation site.
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

## Commands

| Command | Description |
|---------|-------------|
| `rectify correct` | Correct 3' end positions (indel correction, A-tract resolution) |
| `rectify analyze` | Downstream analysis (clustering, DESeq2, GO, motifs) |
| `rectify export` | Export corrected positions to bigWig/bedGraph tracks |
| `rectify extract` | Extract per-read info from BAM to TSV (5'/3' ends, junctions) |
| `rectify aggregate` | Aggregate reads into 3' end, 5' end, and junction datasets |
| `rectify align` | Align FASTQ with minimap2 (DRS-optimized settings) |
| `rectify netseq` | Process NET-seq BAM files (3' end extraction, deconvolution) |
| `rectify run` | Full pipeline: align (if FASTQ) → correct → analyze |

### Examples

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

---

## NET-seq Processing

RECTIFY includes a dedicated pipeline for processing NET-seq (Native Elongating Transcript sequencing) data. NET-seq captures nascent RNA 3' ends, which undergo oligo-adenylation during library preparation, creating characteristic signal spreading at A-tract regions.

```
===============================================================================
NET-SEQ SIGNAL SPREADING: Oligo(A) tails cause position ambiguity
===============================================================================

True CPA site at position P:
                     ↓
Genome:    ...GCTA|AAAAAAAA|TGCG...
                  └────────┘
                   A-tract region

With oligo(A) tails (~5.5 bp mean), reads can prime at ANY downstream A:

             Position:  P   P+1  P+2  P+3  P+4  P+5  P+6  P+7
             Signal:   54%  0.5% 1.5% 3.2% 5.2% 6.7% 7.3% 6.7% ...

The signal "spreads" downstream into the A-tract, obscuring the true CPA site.

===============================================================================
RECTIFY'S SOLUTION: NNLS deconvolution with empirical PSF
===============================================================================

Using the Point-Spread-Function (PSF) derived from 5000+ 0A sites (sites with
no downstream genomic A's, where the true position is known):

  1. Build convolution matrix: A[i,j] = P(observe at j | true peak at i)
  2. Solve NNLS with L2 regularization: min ||Ax - observed||² + λ||x||²
  3. Recover true peak positions from deconvolved signal

Result: Sharper, more accurate 3' end signal with A-tract ambiguity resolved.
```

### NET-seq Command

```bash
# Basic usage
rectify netseq input.bam --genome genome.fa --gff genes.gff -o output/

# With exclusion region control
rectify netseq input.bam --genome genome.fa --gff genes.gff \
    --include-rdna \        # Don't exclude rDNA locus
    --include-pol3 \        # Don't exclude Pol III genes (tRNAs)
    --exclude-mito \        # Exclude mitochondrial genome
    -o output/

# Disable deconvolution (raw 3' ends only)
rectify netseq input.bam --genome genome.fa --no-deconvolution -o output/

# Process multiple samples
rectify netseq sample1.bam sample2.bam sample3.bam \
    --genome genome.fa --gff genes.gff -o output/
```

### Output Files

| File | Description |
|------|-------------|
| `{sample}.unified_reads.parquet` | Per-read records (25 columns, same schema as nanopore) |
| `{sample}.raw.plus.bedgraph` | Raw 3' ends, plus strand, RPM-normalized |
| `{sample}.raw.minus.bedgraph` | Raw 3' ends, minus strand, RPM-normalized |
| `{sample}.deconv.plus.bedgraph` | Deconvolved signal, plus strand, RPM-normalized |
| `{sample}.deconv.minus.bedgraph` | Deconvolved signal, minus strand, RPM-normalized |

### Exclusion Regions

By default, RECTIFY excludes regions with non-standard transcription:

- **rDNA locus** (chrXII ~450,000-490,000 in yeast): Highly repetitive, Pol I transcribed
- **Pol III genes** (tRNAs, SNR6, RDN5, RPR1, SCR1): Different transcription termination mechanism
- **Flanking regions** (100 bp by default): Buffer around excluded genes

These regions are auto-detected from GFF annotation. Use `--include-rdna` and `--include-pol3` flags to include them if needed.

### Strand-Aware Coordinate Mapping

NET-seq reads are short (~40-76 bp) and represent the 3' end of nascent RNA:

| Strand | 3' end position | 3' soft-clip | Oligo(A) detection |
|--------|-----------------|--------------|-------------------|
| **+** | `reference_end - 1` (rightmost) | Right clip | Count A's |
| **-** | `reference_start` (leftmost) | Left clip | Count T's (= RNA A's) |

---

## Supported Technologies

- Nanopore direct RNA-seq (minimap2)
- QuantSeq (oligo-dT short-read)
- PacBio Iso-Seq
- NET-seq (nascent RNA 3' end sequencing)
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
