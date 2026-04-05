# RECTIFY

**R**NA 5' and 3' **E**nd **C**orrection **T**ool with **I**ntron re**F**inement and ambiguit**Y** resolution

[![PyPI version](https://badge.fury.io/py/rectify-rna.svg)](https://pypi.org/project/rectify-rna/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**Precision transcript structure mapping for direct RNA nanopore sequencing data.** Accurate 5' ends, 3' ends, and splice junctions through refinement of alignments at transcript termini, correction of artifacts due to sequencing errors, and selection of the optimal alignments from a panel of aligners (minimap2, gapmm2, and mapPacBio).

---

## Quick Start

```bash
pip install rectify-rna

# Single sample — bundled yeast genome, no external files needed
rectify run-all reads.fastq.gz --Scer -o results/

# Multiple samples via manifest (typical usage)
rectify run-all --manifest samples.tsv --Scer -o results/
```

The manifest `samples.tsv` is a tab-separated file (no header):

```
wt_rep1.fastq.gz      WT    rep1
wt_rep2.fastq.gz      WT    rep2
wt_rep3.fastq.gz      WT    rep3
rna15_rep1.fastq.gz   rna15 rep1
rna15_rep2.fastq.gz   rna15 rep2
rna15_rep3.fastq.gz   rna15 rep3
```

Columns: `filename` (resolved relative to the manifest), `group` (condition label for DESeq2 contrasts), `bio_rep`.

**Outputs** (`results/`):

```
results/
├── <sample_id>/                        # Per-sample
│   ├── <sample_id>.rectified.bam       # Optimal aligner selection + 3' end rectification
│   ├── corrected_3ends.tsv             # Per-read corrected 3' ends, confidence, poly(A) length, fraction
│   ├── corrected_3ends_stats.tsv       # Correction statistics
│   ├── corrected_3ends_report.html     # Per-sample QC report
│   ├── junctions/junctions.tsv         # Splice junctions with partial-rescue evidence
│   └── PROVENANCE.json
│
└── combined/                           # Cross-sample analysis (≥2 samples in manifest)
    ├── cpa_clusters.tsv                # CPA site clusters with per-sample read counts
    ├── tables/
    │   ├── deseq2_genes_*.tsv          # Gene-level differential expression (≥2 conditions)
    │   ├── deseq2_clusters_*.tsv       # Cluster-level differential expression (≥2 conditions)
    │   └── shift_results.tsv           # APA shift analysis (≥2 conditions)
    ├── go_enrichment.tsv               # GO enrichment for DE genes (≥2 conditions)
    ├── motif_results/                  # Enriched sequence motifs near CPA sites
    ├── report.html                     # Combined QC and results report
    └── PROVENANCE.json
```

---

## Key Features

| Feature | Description |
|---------|-------------|
| **Spike-in Filtering** | Removes synthetic spike-in reads (e.g. ENO2) by k-mer sequence matching, restricted to the spike-in gene locus to prevent false positives |
| **5' End Junction Recovery** | Attempts to rescue 5' soft-clipped bases in each aligner's output by extending through splice junctions |
| **3' End A-tract Estimation** | Estimates true CPA position for each aligner using downstream A-tract depth; used to score and break ties in consensus selection |
| **Multi-Aligner Consensus** | Scores rescued alignments from minimap2, mapPacBio, and gapmm2 — penalizing 5' soft-clips and 3' A-tract depth — and selects the best per read |
| **3' End Correction (Walk-back)** | Refines 3' end on the rectified BAM: fixes CIGAR deletion artifacts and walks upstream to the true CPA site; false splice junctions (N operations) are discarded automatically |
| **Poly(A) Measurement** | Reports tail length (aligned + soft-clipped) |
| **Junction Ambiguity Resolution** | Resolves reads matching multiple junctions using proportional assignment |
| **NET-seq Refinement** | Resolves A-tract ambiguity using nascent RNA data; reads assigned proportionally across peaks |
| **CPA Clustering** | Groups CPA sites with fixed-distance clustering (default) or adaptive valley-based algorithm (optional) |
| **Dual-Resolution DESeq2** | Gene-level and cluster-level differential expression |
| **APA Shift Analysis** | Detects proximal/distal CPA site usage changes |
| **Visualization** | Metagene plots, genome browser figures (`pip install rectify-rna[visualize]`) |

### Bundled Data (Yeast)

For *S. cerevisiae*, RECTIFY includes the S288C genome, a merged annotation (see below), GO terms, and pre-deconvolved WT NET-seq data (pan-mutant consensus, NNLS-deconvolved offline)—no external files needed.

**Bundled annotation** (`--Scer`) fuses four sources into a single GFF3:

| Source | Features | Reference |
|--------|----------|-----------|
| SGD R64-5-1 | All protein-coding genes, tRNAs, snoRNAs, rDNA, etc. | SGD (yeastgenome.org) |
| CUTs (925) | Cryptic Unstable Transcripts | Xu et al. 2009, *Nature* 457:1033 ([PMID:19169243](https://pubmed.ncbi.nlm.nih.gov/19169243/)) |
| SUTs (847) | Stable Unannotated Transcripts | Xu et al. 2009, *Nature* 457:1033 ([PMID:19169243](https://pubmed.ncbi.nlm.nih.gov/19169243/)) |
| XUTs (1,658) | Xrn1-sensitive Unstable Transcripts | van Dijk et al. 2011, *Nature* 475:114 ([PMID:21697827](https://pubmed.ncbi.nlm.nih.gov/21697827/)) |

CUT/SUT coordinates are SGD-curated and lifted to R64-1-1. XUT coordinates are from the original van Dijk 2011 supplementary data (R64/sacCer3). All features use Roman numeral chromosome names (`chrI`–`chrXVI`, `chrMito`).

---

## 3' End Correction: Indel Artifacts in Poly(A) Regions

When poly(A) tails align to genomic A-tracts, aligners like minimap2 introduce indel artifacts to maximize alignment score. This shifts the apparent 3' end **downstream** of the true cleavage site.

![3' End Indel Correction](docs/figures/indel_correction.png)

**The Problem:** The aligner introduces deletions to extend alignment of poly(A) tail bases into downstream genomic A-tracts, shifting the apparent 3' end.

**RECTIFY's Solution:** Walk upstream from the soft-clip boundary through the aligned region:
1. Skip positions where genome = A (ambiguous with poly(A) tail)
2. Skip deletions (D in CIGAR) - these are alignment artifacts
3. Skip T sequencing errors in the tail
4. **STOP** at first non-A/T agreement between genome and read

**Result:** True 3' end at the correct position. Poly(A) length = soft-clipped + aligned A's + absorbed deletions.

**Key insight for IGV users:** For minus strand reads, the poly(A) tail appears as poly(T) extending leftward. RECTIFY corrects by shifting rightward toward the true CPA.

---

## Soft-Clip Rescue at Homopolymer Boundaries

Nanopore basecallers systematically under-call homopolymer runs (e.g., calling 8 U's instead of 10). When this happens at CPA sites with upstream T-tracts, the aligner soft-clips the non-T bases instead of placing them correctly.

![Soft-Clip Rescue](docs/images/softclip_rescue.png)

**The Problem:** The basecaller under-calls the T-tract, so when the aligner reaches the first non-T base, it can't fit it into the shortened homopolymer and soft-clips it instead.

**RECTIFY's Solution:**
1. Detect soft-clips adjacent to homopolymer boundaries
2. Skip over remaining reference homopolymer bases (the under-called T's)
3. Match soft-clipped bases to reference sequence beyond the homopolymer
4. Extend the 3' end to include the rescued bases

---

## 5' End Correction: Splice Junction Soft-Clips

Long reads spanning splice junctions often have soft-clipped bases at the 5' end where the aligner fails to find the exact junction boundary. RECTIFY recovers the true splice site.

![5' Junction Rescue](docs/figures/5prime_junction_rescue.png)

**The Problem:** Soft-clipped bases at the 5' end actually match the upstream exon, but the aligner couldn't extend through the splice junction.

**RECTIFY's Solution:**
1. Identify reads with 5' soft-clips near annotated splice sites
2. Check if soft-clipped sequence matches upstream exon
3. Extend the alignment to the canonical splice donor (GT)

**Result:** Accurate read 5' end and recovered splice junction. Note: Due to 5'-to-3' degradation in direct RNA sequencing, the read's 5' end often does not represent the true TSS.

---

## Multi-Aligner Rectification Pipeline

Different aligners make different tradeoffs at splice junctions. RECTIFY runs minimap2, mapPacBio, and gapmm2, selects the optimal alignment per read, then applies 3' end rectification.

![Multi-Aligner Rectification](docs/figures/multi_aligner_consensus.png)

**The Problem:** Different aligners handle the same read differently. Some soft-clip at splice boundaries while others find the junction.

**RECTIFY's Solution:**
1. Run all 3 aligners (minimap2, mapPacBio, gapmm2) on the same reads
2. **Attempt to rescue** each alignment's 5' soft-clips by extending through known junctions
3. Score rescued alignments: penalize remaining 5' soft-clips and 3' A-tract depth (GT-AG motifs serve as tiebreakers only — deliberately excluded from primary score to avoid penalizing novel non-canonical junctions)
4. Select the **optimal** alignment per read
5. Apply 3' end rectification; output **rectified BAM**

**Note:** 3' false junctions from poly(A) artifacts are handled separately by walk back correction (see "3' False Junction Handling" below).

**Usage:**

```bash
# Multi-aligner rectified alignment (default, aligners run sequentially)
rectify align reads.fastq.gz --genome genome.fa --annotation genes.gff -o aligned/

# Parallel aligners with proportional thread allocation (minimap2 gets fewer threads
# since it's faster; mapPacBio and gapmm2 get more to finish at the same time)
rectify align reads.fastq.gz --genome genome.fa --annotation genes.gff \
    --parallel-aligners --threads 16 -o aligned.bam

# Single aligner mode — outputs optimal alignment from one aligner only
rectify align reads.fastq.gz --genome genome.fa --aligners minimap2 -o aligned/
```

---

## 3' False Junction Handling

Poly(A) tails can create spurious "junctions" when the aligner introduces a skip (N) operation to align tail bases to a downstream genomic A-tract. **RECTIFY's walk back correction completely handles this artifact.**

![False Junction Walk Back](docs/figures/false_junction_walkback.png)

**The Problem:** The aligner introduces an N (skip) to extend the poly(A) tail alignment into a downstream A-tract, creating a spurious junction that doesn't exist in the transcript.

**RECTIFY's Solution:** The walk back algorithm finds the true 3' end by walking upstream through ALL aligned A's until it finds the first non-A agreement between genome and read. N (skip/junction) operations are excluded from the walk-back position list during CIGAR parsing — so they are never traversed, and any spurious junction is simply absent from the corrected output.

**Result:**
- Walk back finds the true CPA at the EXON/A boundary
- The false junction (N operation) is completely ignored
- No special false junction detection needed

**This means 3' false junctions are NOT a problem for downstream analysis**—walk back correction handles them automatically as part of finding the true cleavage and polyadenylation site.

---

## Adaptive Clustering and Differential Expression

After correction, RECTIFY groups nearby CPA sites into clusters, then runs DESeq2 at both gene and cluster resolution. Two modes: **fixed-distance** (default — merge sites within a set bp window) and **adaptive valley-based** (optional, `--adaptive-clustering` — find signal peaks then split at valleys).

![Adaptive Clustering](docs/figures/adaptive_clustering.png)

**Adaptive valley-based algorithm** (`--adaptive-clustering`):
1. Find peaks (local maxima in 3' end signal)
2. Find valleys (local minima between peaks)
3. Set boundaries at midpoint between peak and valley (capped at ±10bp)

**Why cluster-level analysis matters:**
- Genes often have MULTIPLE CPA sites (alternative polyadenylation)
- Conditions may shift usage between proximal/distal sites
- Cluster-level DESeq2 detects isoform-specific changes that gene-level misses

**Dual-resolution output:**

| Level | Detects | Example |
|-------|---------|---------|
| **Gene** | Total expression changes | HSP82 is 2-fold down in heat shock |
| **Cluster** | CPA site usage changes | FAS1 shifts from distal to proximal site |

---

## Output

Each read gets a corrected position with confidence scores:

```
read_id   │ chrom │ strand │ original_3prime │ corrected_3prime │ confidence │ polya_length │ fraction │ qc_flags
read001   │ chrI  │   +    │    147592       │     147585       │    high    │      42      │  1.0000  │   PASS
read002   │ chrI  │   +    │    147594       │     147591       │   medium   │      38      │  1.0000  │   PASS
read003   │ chrII │   +    │    283109       │     283104       │    low     │      31      │  0.6500  │ AG_RICH
```

The `fraction` column reflects proportional NET-seq assignment: when a read falls in an A-tract with multiple NET-seq peaks, it is split across peaks rather than snapped to a single position. Fractions sum to 1.0 per input read.

The `rectify analyze` command produces:
- `cpa_clusters.tsv` - CPA site clusters with read counts per sample
- `deseq2_genes_{condition}.tsv` - Gene-level differential expression
- `deseq2_clusters_{condition}.tsv` - Cluster-level differential expression
- `shift_analysis_{condition}.tsv` - Genes with significant APA shifts
- `go_enrichment_up_{condition}.tsv` / `go_enrichment_down_{condition}.tsv` - GO enrichment for up/down DE genes
- `motif_results/` - Enriched sequence motifs near CPA sites

#### Genomic category distribution

The primary per-sample output includes a horizontal bar chart (`genomic_distribution.png`) showing the fraction of reads and number of unique clusters falling into each genomic category — directly comparable to the style of Xu et al. 2009 Figure 1B/C:

```
3' UTR            ████████████████████████████████████████  ~90%
snoRNA +/- 300 bp ██  ~2%
intergenic        ███  ~3%
CUTs              █  ~1%
SUTs / XUTs       █  ~1%
5' UTR / CDS      ██  ~2%
antisense CDS     █  <1%
```

The x-axis uses a broken scale so the dominant 3' UTR bar and the minority categories are both visible. Multiple conditions (e.g. WT vs *rrp6*Δ) are overlaid as separate bars.

**Category definitions and classification priority:**

| Category | Feature types | Source |
|----------|--------------|--------|
| 3' UTR | `UTR3` annotation (or 3' 10% heuristic) | SGD |
| snoRNA +/- 300 bp | `snoRNA_gene` ± 300 bp on either strand | SGD |
| CUTs | `CUT` features | Xu et al. 2009 |
| SUTs / XUTs | `SUT` or `XUT` features | Xu 2009 / van Dijk 2011 |
| 5' UTR / CDS | `UTR5` or `CDS` annotation | SGD |
| antisense CDS | Position overlaps CDS on opposite strand | SGD |
| intergenic / intronic | None of the above | — |

When a position matches multiple categories, the highest-priority rule wins (order as listed above).

The companion `genomic_distribution_3prime_summary.tsv` table contains:

```
condition  category  display_label        reads   reads_pct  clusters
WT         UTR3      3' UTR               182034  91.2       4812
WT         snoRNA    snoRNA +/- 300 bp    3940    1.97       208
...
```

#### Transcript body distribution

A second figure (`genomic_distribution_body.png`) classifies each read by RNA biotype based on the majority-overlap of its full alignment span with annotated features. Replicates are merged by condition; one pie chart is shown per condition.

**Biotype assignments and priority:**

| Biotype | Feature types matched | Source |
|---------|-----------------------|--------|
| protein-coding | `mRNA`, `CDS` | SGD |
| CUTs | `CUT` | Xu et al. 2009 |
| SUTs / XUTs | `SUT`, `XUT` | Xu 2009 / van Dijk 2011 |
| snoRNA | `snoRNA_gene`, `snoRNA` | SGD |
| tRNA | `tRNA_gene`, `tRNA` | SGD |
| rRNA | `rRNA_gene`, `rRNA` | SGD |
| LTR / retrotransposon | `LTR_retrotransposon`, `transposable_element_gene` | SGD |
| pseudogene | `pseudogene` | SGD |
| antisense | Read span overlaps an annotated feature on opposite strand | — |
| intergenic | No annotated feature overlaps the alignment span | — |

The read is assigned to whichever overlapping feature has the most base-pair overlap. Ties are broken by the priority order above (protein-coding > CUT > SUT_XUT > ...).

The companion `genomic_distribution_body_summary.tsv` table contains:

```
condition  category      display_label    reads   reads_pct
WT         protein_coding  protein-coding  175000  87.5
WT         CUT             CUTs             4500   2.25
...
```

---

## NET-seq Refinement (Optional)

For organisms with NET-seq data, RECTIFY resolves remaining A-tract ambiguity by assigning reads proportionally across NET-seq peaks rather than snapping to a single position.

**Bundled yeast data** (`--Scer`): pre-deconvolved pan-mutant NET-seq (WT + DST1Δ consensus, NNLS deconvolution applied once offline). No runtime deconvolution — reads are assigned directly using the pre-computed signal.

**Custom data**: provide raw NET-seq bigWigs with `--netseq-dir`. NNLS deconvolution is applied at runtime to recover true CPA positions from the oligo-adenylation spreading pattern.

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

### All-in-one

| Command | Description |
|---------|-------------|
| `rectify run-all` | Full pipeline: align (if FASTQ) → correct → analyze → aggregate junctions. Skips completed steps automatically on re-run. |

### Individual steps

Run steps independently to re-process from any point in the pipeline.

| Command | Description |
|---------|-------------|
| `rectify align` | Align FASTQ with the aligner panel (minimap2, gapmm2, mapPacBio; all with DRS-optimized settings) select the optimal alignment per read, and output a rectified BAM |
| `rectify correct` | Correct 5' ends, 3' ends, and junctions — indel correction at poly(A) boundaries, A-tract ambiguity resolution, NET-seq refinement |
| `rectify analyze` | Downstream analysis: CPA clustering, DESeq2, GO enrichment, motif discovery |
| `rectify export` | Export corrected positions to bigWig/bedGraph tracks |
| `rectify extract` | Extract per-read features from BAM to TSV (5'/3' ends, junctions, poly(A) length) |
| `rectify aggregate` | Aggregate reads into 3' end, 5' end, and junction-centered datasets |
| `rectify netseq` | Process NET-seq BAM files (3' end extraction, NNLS deconvolution; for standalone analysis or for assigning DRS 3' ends in A-tracts to their likely CPA sites) |

### Examples

```bash
# Multiple samples via manifest — typical multi-condition experiment
rectify run-all --manifest samples.tsv --Scer --filter-spikein ENO2 -o results/

# Single sample, bundled yeast genome
rectify run-all reads.fastq.gz --Scer -o results/

# Non-DRS protocol where poly(A) tail is not sequenced
rectify run-all reads.fastq.gz --genome genome.fa --annotation genes.gff \
    --no-polya-sequenced -o results/

# Re-run correction only (alignment already done)
rectify correct reads.bam --genome genome.fa --netseq-dir my_netseq/ -o corrected.tsv

# Re-run analysis only (correction already done)
rectify analyze corrected.tsv --annotation genes.gff --output-dir results/

# Process NET-seq data
rectify netseq netseq.bam --genome genome.fa --gff genes.gff -o netseq_output/
```

---

## NET-seq Processing

RECTIFY includes a dedicated pipeline for processing NET-seq (Native Elongating Transcript sequencing) data. NET-seq captures nascent RNA 3' ends, which undergo oligo-adenylation during library preparation, creating characteristic signal spreading at A-tract regions.

### The Problem: Oligo(A) Signal Spreading

![Oligo(A) Spreading](docs/figures/oligo_a_spreading.png)

With oligo(A) tails (~5.5 bp mean), reads can prime at ANY downstream A in genomic A-tracts. This causes signal to "spread" downstream, obscuring the true CPA site.

### RECTIFY's Solution: NNLS Deconvolution

![Oligo(A) Deconvolution](docs/figures/oligo_a_deconvolution.png)

Using the Point-Spread-Function (PSF) derived from 5000+ 0A sites (sites with no downstream genomic A's, where the true position is known):

1. Build convolution matrix: A[i,j] = P(observe at j | true peak at i)
2. Solve NNLS with L2 regularization: min ||Ax - observed||² + λ||x||²
3. Recover true peak positions from deconvolved signal

**Result:** Sharper, more accurate 3' end signal with A-tract ambiguity resolved.

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

## Visualization

Install with: `pip install rectify-rna[visualize]`

### Metagene: single condition

Aggregate 3' end signal around a set of loci (TRT sites, CPA clusters, gene TES/TSS, motif matches).

```python
import numpy as np
import matplotlib.pyplot as plt
from rectify.visualize import (
    MetagenePipeline,
    position_index_from_tsv,
    loci_from_pickle,
    loci_from_gff,
    plot_metagene_line,
    add_metagene_annotations,
    set_publication_style,
    despine,
    WONG_COLORS,
)

# Build a position index from RECTIFY corrected 3' end TSVs
# Multiple replicates merge automatically
index, total_reads = position_index_from_tsv(
    ["wt_rep1/corrected_3ends.tsv",
     "wt_rep2/corrected_3ends.tsv",
     "wt_rep3/corrected_3ends.tsv"],
    position_col='corrected_3prime',
)

# Load loci — choose the loader that matches your input
trt_loci = loci_from_pickle("cache/trt_signals_v3.pkl")            # TRT/CPA pkl cache
tes_loci = loci_from_gff("genes.gff", feature_type='gene', center='end')   # gene 3' ends
tss_loci = loci_from_gff("genes.gff", feature_type='gene', center='start') # gene 5' ends

# Compute metagene (strand coordinate transformation + reversal handled internally)
pipeline = MetagenePipeline()
result = pipeline.compute_center_profile(
    trt_loci, index,
    window=(-50, 50),
    total_reads=total_reads,
    verify_strands=True,   # raises StrandOrientationError on strand bug — never skip this
    cap_percentile=50,     # suppress outlier loci (e.g. TDH3) before peak detection
)

set_publication_style()
fig, ax = plt.subplots(figsize=(5, 3), constrained_layout=True)

x    = result['x']
mean = np.mean(result['profile_matrix'], axis=0)
sem  = np.std(result['profile_matrix'],  axis=0) / np.sqrt(result['n_loci'])

plot_metagene_line(ax, x, mean, sem,
                   color=WONG_COLORS['blue'],
                   label=f"WT (n={result['n_loci']})")
add_metagene_annotations(ax, positions=[0], labels=["TRT start"])
ax.set_xlabel("Position (bp from TRT start)")
ax.set_ylabel("RPM")
despine(ax)
ax.legend(fontsize=8)
fig.savefig("trt_metagene_wt.png", dpi=150, bbox_inches='tight')
```

**Always set `cap_percentile=50`** with window ≥ 100 bp — a single highly expressed gene
can dominate the aggregate and mask the biological peak.

---

### Metagene: multi-condition ridge plot

Compare WT vs mutants at TRT sites, stacked as a ridge plot (conditions offset vertically
so shapes can be compared without overlap).

```python
import numpy as np
import matplotlib.pyplot as plt
from rectify.visualize import (
    MetagenePipeline,
    position_index_from_tsv,
    loci_from_pickle,
    plot_ridge_profiles,
    apply_window_sum_capping,
    set_publication_style,
    WONG_COLORS,
)

trt_loci = loci_from_pickle("cache/trt_signals_v3.pkl")
pipeline  = MetagenePipeline()

# Build per-condition profile matrices
CONDITIONS = {
    'wt':     ["wt_rep1/corrected_3ends.tsv",     "wt_rep2/corrected_3ends.tsv"],
    'rna15':  ["rna15_rep1/corrected_3ends.tsv",  "rna15_rep2/corrected_3ends.tsv"],
    'ysh1':   ["ysh1_rep1/corrected_3ends.tsv",   "ysh1_rep2/corrected_3ends.tsv"],
}

profiles = {}
for cond, paths in CONDITIONS.items():
    idx, n = position_index_from_tsv(paths, position_col='corrected_3prime')
    r = pipeline.compute_center_profile(
        trt_loci, idx, window=(-50, 50), total_reads=n,
        verify_strands=True, cap_percentile=50,
    )
    profiles[cond] = r['profile_matrix']

# Optional: cap outliers consistently across conditions before ridge plot
profiles_capped = {k: apply_window_sum_capping(v, percentile=90)[0]
                   for k, v in profiles.items()}

x = np.arange(-50, 51)
COLORS = {'wt': WONG_COLORS['blue'], 'rna15': WONG_COLORS['vermillion'],
          'ysh1': WONG_COLORS['green']}
LABELS = {'wt': 'WT', 'rna15': 'RNA15-AA', 'ysh1': 'YSH1-AA'}

set_publication_style()
fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
plot_ridge_profiles(ax, x, profiles_capped, colors=COLORS, labels=LABELS,
                    order=['wt', 'rna15', 'ysh1'])
ax.set_xlabel("Position (bp from TRT start)")
ax.axvline(0, color='#CC0000', linestyle='--', linewidth=0.8, alpha=0.7)
fig.savefig("trt_metagene_ridge.png", dpi=150, bbox_inches='tight')
```

---

### Stacked read browser

Genome-browser-style view of individual nanopore reads, colored by transcript category.
Uses `LineCollection` for fast batch rendering — 400 reads × 5 conditions renders as
~30 draw calls instead of 2000.

```python
import pandas as pd
import matplotlib.pyplot as plt
from rectify.visualize import plot_stacked_read_panel

# ── Read classification is your responsibility (domain-specific) ───────────────
CAT_COLORS = {
    'CDS-internal':   '#d62728',  # red
    '3UTR-premature': '#ff9896',  # coral
    'intronic':       '#9467bd',  # purple
    'canonical':      '#2ca02c',  # green
    'readthrough':    '#e6ac00',  # gold
    'other':          '#aaaaaa',  # gray
}

# Load reads from RECTIFY TSV, filter to gene region
reads_df = pd.read_csv("wt_rep1/corrected_3ends.tsv", sep='\t',
                       usecols=['read_id', 'chrom', 'strand',
                                'corrected_3prime', 'five_prime_position',
                                'alignment_start', 'alignment_end',
                                'qc_flags', 'junctions'])
reads_df = reads_df[
    (reads_df['qc_flags'] == 'PASS') &
    (reads_df['chrom'] == 'chrXIV') &
    (reads_df['strand'] == '+') &
    (reads_df['alignment_start'] < 417000) &
    (reads_df['alignment_end']   > 413000)
]

# Sort by 3' end so reads with the same termination site appear contiguous
reads_df = reads_df.sort_values('corrected_3prime')

# Classify and map to color (your logic here)
reads_df['color'] = reads_df.apply(classify_read, axis=1).map(CAT_COLORS)

# ── Draw ───────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    3, 1,                                         # gene track + 2 conditions
    figsize=(14, 8),
    gridspec_kw={'height_ratios': [1, 3, 3]},
    constrained_layout=True,                      # NOT tight_layout()
)

# condition panel
n_rows = plot_stacked_read_panel(
    axes[1], reads_df,
    color_col='color',
    start_col='alignment_start',
    end_col='alignment_end',
    junction_col='junctions',  # "start-end,start-end" RECTIFY format
    gap=50,
)

# mark internal pA sites
for pos in [414500, 415200]:
    axes[1].axvline(pos, color='#d62728', linestyle='--', linewidth=1.0, alpha=0.8)

axes[1].set_xlim(413000, 417000)
axes[1].set_ylabel("WT", rotation=0, ha='right', va='center')
```

---

### Locus loaders

| Function | Input | Use for |
|----------|-------|---------|
| `loci_from_pickle(path)` | `.pkl` cache | TRT/CPA caches |
| `loci_from_gff(path, feature_type, center)` | GFF3 | TSS (`center='start'`), TES (`center='end'`) |
| `loci_from_motif_scan(sequences, motif)` | dict of sequences | Regex motif, both strands (e.g. `T{6,}`) |
| `loci_from_bed(path, center)` | BED file | Custom interval lists |
| `loci_from_tsv(path)` | TSV with chrom/strand/center | Pre-computed loci |
| `position_index_from_tsv(paths)` | RECTIFY `corrected_3ends.tsv` | 3' end signal |
| `position_index_from_bigwig(path)` | bigWig | NET-seq, PAR-CLIP, any coverage |

---

### Figure utilities

```python
from rectify.visualize import set_publication_style, save_multi_format, WONG_COLORS

set_publication_style()   # consistent fonts, tick sizes, rc params for all figures

# Color-blind-safe palette (use for all new figures)
blue       = WONG_COLORS['blue']        # #0072B2
orange     = WONG_COLORS['orange']      # #E69F00
green      = WONG_COLORS['green']       # #009E73
vermillion = WONG_COLORS['vermillion']  # #D55E00

# Save PNG + PDF + SVG in one call
save_multi_format(fig, "figures/fig1_trt_metagene")
# → figures/fig1_trt_metagene.png, .pdf, .svg
```

See [PLOT_SKILLS.md](PLOT_SKILLS.md) for the full API reference and list of pitfalls.

---

## Supported Technologies

- Nanopore direct RNA-seq (minimap2)
- QuantSeq (oligo-dT short-read)
- PacBio Iso-Seq (experimental — mapPacBio/gapmm2 aligners available)
- NET-seq (nascent RNA 3' end sequencing)
- Any poly(A)-tailed RNA-seq

---

## Citation

**Original RECTIFY:**
> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020;176:4-13. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.x (current):** Manuscript in preparation

---

## License

MIT License - See [LICENSE](LICENSE) for details

## Contact

- Kevin R. Roy - [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com)
- GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)
