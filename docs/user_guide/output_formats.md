# Output Formats

## Correction outputs

### `corrected_3ends.tsv`

The primary per-read output. One row per read.

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | str | Read name from BAM |
| `chrom` | str | Chromosome (UCSC format, e.g. `chrI`) |
| `strand` | str | `+` or `-` |
| `original_3prime` | int | Raw 3' end before correction (0-based) |
| `corrected_3prime` | int | Corrected 3' end after walk-back (0-based) |
| `correction_modules` | str | Comma-separated list of modules that modified this read |
| `confidence` | str | `HIGH`, `MEDIUM`, or `LOW` |
| `polya_length` | int | Total poly(A) tail length (aligned + soft-clipped A's) |
| `polya_source` | str | Source of poly(A) length estimate (`aligned`, `softclip`, `pt_tag`, etc.) |
| `five_prime_rescued` | int | `1` if 5' end was rescued by Cat3 splice-site truncation rescue |
| `five_prime_exon_cigar` | str | SAM CIGAR string for the rescued exon segment (e.g. `8M1D3M`); empty if no rescue |
| `sc_homopolymer_extension` | int | Number of bases extended into genomic homopolymer (Module 2G) |
| `sc_rescued_seq` | str | Soft-clipped sequence rescued by Module 2G |
| `sc_original_softclip_len` | int | Original soft-clip length before Module 2G rescue |
| `best_aligner` | str | Winning aligner from consensus selection |
| `qc_flags` | str | `PASS`, `AG_RICH`, `ATRACT_AMBIGUOUS`, `LOW_MAPQ`, etc. |

!!! note "Coordinate convention"
    All positions are **0-based, half-open** (BED/pysam convention). See the [Coordinate System](../coordinate_system.md) page for strand-specific definitions.

---

### `corrected_3ends_index.bed.gz`

A pre-aggregated position count file — approximately 300× smaller than the full TSV. Used by manifest-mode analysis to speed up multi-sample clustering.

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome |
| `corrected_3prime` | Corrected 3' end position (0-based) |
| `strand` | `+` or `-` |
| `count` | Number of reads at this position |

---

### BAM outputs from `rectify correct`

| File | Description |
|------|-------------|
| `rectified_corrected_3end.bam` | Corrected BAM with reads hard-clipped to the corrected 3' end |
| `rectified_pA_tail_trimmed.bam` | Soft-clipped BAM showing the poly(A) tail as soft-clips at the 3' end |
| `rectified_pA_tail_soft_clipped.bam` | DRS Step 4 output: poly(A) tail restored as soft-clips after `rectify restore-softclip` |

!!! note "DRS Step 4"
    `rectified_pA_tail_soft_clipped.bam` is produced only when `--drs` is passed to `rectify run-all`. It is written by Step 4 (`rectify restore-softclip`) and is sorted and indexed.

---

### `alignment_features.tsv`

Per-read alignment metadata.

| Column | Description |
|--------|-------------|
| `read_id` | Read name |
| `mapq` | Mapping quality |
| `cigar_summary` | Compact CIGAR description |
| `alignment_identity` | Fraction of aligned bases matching reference |
| `five_prime_soft_clip` | Soft-clipped bases at 5' end |
| `three_prime_soft_clip` | Soft-clipped bases at 3' end |

---

### `processing_stats.tsv`

Per-sample QC summary.

| Column | Description |
|--------|-------------|
| `total_reads` | Total reads processed |
| `corrected` | Reads with position correction applied |
| `ag_rich_flagged` | Reads flagged for AG-rich context |
| `atract_ambiguous` | Reads in A-tract ambiguous region |
| `low_mapq_filtered` | Reads filtered for low mapping quality |
| `spikein_filtered` | Spike-in reads removed |
| `mean_shift_bp` | Mean shift (usually −3 to −7 bp) |
| `median_polya_length` | Median poly(A) tail length |

---

## Analysis outputs

### `cpa_clusters.tsv`

CPA site clusters from adaptive clustering.

| Column | Description |
|--------|-------------|
| `cluster_id` | Unique cluster identifier |
| `chrom` | Chromosome |
| `strand` | `+` or `-` |
| `start` | Cluster start (0-based) |
| `end` | Cluster end (exclusive) |
| `peak_pos` | Most-used CPA position within cluster |
| `total_count` | Total reads across all samples |

---

### `tables/deseq2_genes_{condition}.tsv`

Gene-level differential expression results.

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier |
| `baseMean` | Average normalized expression |
| `log2FoldChange` | Log2 fold change (condition vs reference) |
| `lfcSE` | Standard error of fold change |
| `stat` | Wald statistic |
| `pvalue` | Nominal p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |

---

### `tables/deseq2_clusters_{condition}.tsv`

Cluster-level (isoform-resolution) differential expression. Same columns as gene-level, but `gene_id` is replaced by `cluster_id`.

---

### `tables/shift_results.tsv`

Genes with significant APA site usage shifts.

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier |
| `js_divergence` | Jensen-Shannon divergence between conditions |
| `direction` | `proximal_shift` or `distal_shift` |
| `proximal_delta` | Change in proximal site usage fraction |
| `distal_delta` | Change in distal site usage fraction |
| `padj` | Adjusted p-value |

---

### `tables/go_enrichment_{condition}.tsv`

GO term enrichment.

| Column | Description |
|--------|-------------|
| `go_term` | GO term ID |
| `go_name` | GO term description |
| `namespace` | `biological_process`, `molecular_function`, `cellular_component` |
| `n_genes` | Number of genes in term |
| `n_sig` | Significant genes in term |
| `fold_enrichment` | Observed / expected |
| `pvalue` | Hypergeometric p-value |
| `padj` | Benjamini-Hochberg adjusted |

---

### `tables/motif_summary_{condition}.tsv`

Enriched sequence motifs near CPA sites (from STREME/MEME).

| Column | Description |
|--------|-------------|
| `motif_id` | STREME/MEME motif ID |
| `consensus` | Consensus sequence |
| `evalue` | E-value |
| `n_sites` | Number of matched sites |
| `offset` | Distance from CPA site |

Full motif files (logo images, MEME format) are in `motifs/{condition}/`.

---

## Export outputs

### `{prefix}.bw`

BigWig format — per-base 3' end coverage, strand-separated. Load directly in IGV, UCSC Genome Browser, or pyGenomeTracks.

### `{prefix}.bedgraph`

BedGraph format — same data as bigWig, plain text.

---

## Visualization outputs

### `plots/pca_samples.png`

PCA of samples by CPA cluster usage (first 2 PCs).

### `plots/heatmap_top_clusters.png`

Heatmap of top differentially used CPA clusters across samples.

### `analysis_report.html`

Comprehensive HTML report including QC stats, PCA, heatmaps, top DE genes, and motif logos.
