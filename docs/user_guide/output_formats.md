# Output Formats

## Correction outputs

### `corrected_reads.tsv`

The primary per-read output of `rectify correct` (and the per-sample output of
`rectify run-all`). One row per read; columns follow the canonical
`CORRECTION_TSV_HEADER`. The most commonly used columns:

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | str | Read name from BAM |
| `chrom` | str | Chromosome (UCSC format, e.g. `chrI`) |
| `strand` | str | `+` or `-` |
| `original_3prime` | int | Raw 3' end before correction (0-based) |
| `corrected_3prime` | int | Corrected 3' end after walk-back (0-based) |
| `five_prime_position` | int | 5' end position (0-based) |
| `five_prime_rescued` | int | `1` if the 5' end was rescued by Cat3 splice-site truncation rescue |
| `five_prime_exon_cigar` | str | SAM CIGAR for the rescued exon segment (e.g. `8M1D3M`); empty if no rescue |
| `alignment_start` / `alignment_end` | int | Reference alignment span (0-based half-open) |
| `ambiguity_min` / `ambiguity_max` / `ambiguity_range` | int | A-tract ambiguity window bounds and width |
| `polya_length` | int | Poly(A) tail length (aligned + soft-clipped A's) |
| `aligned_a_length` / `soft_clip_a_length` | int | Poly(A) split into aligned vs soft-clipped components |
| `junctions` / `n_junctions` | str / int | Splice junctions in the read and their count |
| `five_prime_soft_clip_length` / `three_prime_soft_clip_length` | int | Terminal soft-clip lengths |
| `mapq` | int | Mapping quality |
| `correction_applied` | int/bool | Whether a position correction was applied |
| `confidence` | str | `HIGH`, `MEDIUM`, or `LOW` |
| `qc_flags` | str | e.g. `PASS`, `AG_RICH`, `ATRACT_AMBIGUOUS`, `LOW_MAPQ` |
| `gene_id` | str | Assigned gene (annotation overlap) |
| `pt_tag` | int | Dorado `pt:i` poly(A) estimate, if present |
| `polya_score` | float | Poly(A) model score |
| `polya_source` | str | Source of the poly(A) estimate (`aligned`, `softclip`, `pt_tag`, etc.) |
| `sc_homopolymer_extension` | int | Bases extended into a genomic homopolymer (Module 2G) |
| `sc_rescued_seq` | str | Soft-clipped sequence rescued by Module 2G |
| `sc_original_softclip_len` | int | Original soft-clip length before Module 2G rescue |
| `oc_homopolymer_extension` / `oc_overcall_count` / `oc_terminal_base` | int/str | Over-call (over-extension) diagnostics |
| `five_prime_intron_clip_pos` / `five_prime_upstream_trim` / `reanchor_clip_len` | int | 5' rescue / re-anchoring diagnostics |
| `strand_evidence` | str | How `strand` was decided; ONT PCR-cDNA only (`polyA_3p`, `polyT_5p`, `gene_overlap`, `unassigned`) |
| `consensus_aligner` | str | `Xa` — aligner the multi-aligner consensus chose for this read (comma list for a chimeric stitch) |
| `consensus_confidence` | str | `Xc` — consensus confidence for the choice |
| `consensus_n_agree` | int | `Xn` — how many aligners agreed |
| `consensus_tied` | str | `Xt` — comma-separated tied aligners; empty unless the vote tied |

The full 43-column header is defined by `CORRECTION_TSV_HEADER`
(`rectify/core/bam/output.py`). After consensus merge (`rectify run-all`,
`rectify consensus`), a `winning_aligner` column is appended recording the
aligner selected per read.

!!! note "`consensus_*` vs `winning_aligner`"
    The four `consensus_*` columns are the input BAM's `Xa`/`Xc`/`Xn`/`Xt`
    tags, which the **alignment** stage writes onto
    `<sample>.multialigned.bam`. `rectify correct` copies them straight
    through when it reads a BAM that carries them. In the default order it
    does not: it corrects each per-aligner BAM, and those carry no tags — so
    the **merge** joins the tags back in by read name from
    `<sample>.multialigned.bam`, which it finds next to the per-aligner BAMs
    (override with `--consensus-bam`; the path used is logged). Only empty
    cells are filled, and the columns stay empty when no such BAM exists.
    `winning_aligner` is the separate, later verdict of the merge itself —
    these columns are never derived from it.

!!! note "Coordinate convention"
    All positions are **0-based, half-open** (BED/pysam convention). See the [Coordinate System](../coordinate_system.md) page for strand-specific definitions.

!!! note "`--use-dorado-polya`"
    By default, RECTIFY measures poly(A) length from the alignment/soft-clip
    and records Dorado's `pt:i` estimate in `pt_tag` for reference. Passing
    `--use-dorado-polya` to `rectify correct` makes Dorado's `pt:i` value the
    authoritative `polya_length` instead.

---

### BAM outputs from `rectify correct`

These BAMs are optional and written to the paths you pass on the command line:

| Flag | Description |
|------|-------------|
| `--write-corrected-bam BAM` | Corrected BAM with each read hard-clipped at its corrected 3' end |
| `--write-softclipped-bam BAM` | Corrected BAM with poly(A) tail bases soft-clipped (rather than hard-clipped) at the corrected 3' end |
| `--output-bam BAM` | The corrected alignments written back as a BAM (see `rectify correct --help`) |

In `rectify run-all`, the per-sample corrected consensus alignments are written
as `corrected_consensus.bam` (sorted + indexed). DRS Step 4
(`rectify restore-softclip`) re-attaches the trimmed poly(A) tail as a 3'
soft-clip for IGV inspection — produced only when `--drs` is passed to
`rectify run-all`, or by invoking `rectify restore-softclip` directly.

---

### `<output>_stats.tsv`

Per-sample QC summary, written alongside the corrected TSV (the `.tsv`
extension of `-o`/`--output` is replaced with `_stats.tsv`). Each row is a
`metric  count  percent  description` record; common metrics include:

| Metric | Description |
|--------|-------------|
| `total_reads_in_bam` | Total reads in the input BAM |
| `reads_processed` | Reads with 3' ends corrected |
| `spikein_reads_filtered` | Spike-in reads filtered |
| `ends_corrected_indel` | 3' ends corrected for indel artifacts |
| `ends_walkback_readgenome` | 3' ends corrected by read-vs-reference walkback (`--dT-primed-cDNA`) |
| `total_position_shifts` | Total reads whose 3' position changed |

---

## Analysis outputs

### `cpa_clusters.tsv`

CPA site clusters from adaptive clustering.

| Column | Description |
|--------|-------------|
| `cluster_id` | Unique cluster identifier (e.g. `cluster_000123`) |
| `chrom` | Chromosome |
| `strand` | `+` or `-` |
| `start` | Cluster start (0-based) |
| `end` | Cluster end (exclusive) |
| `modal_position` | Most-used CPA position within the cluster |
| `cluster_com` | Read-weighted center-of-mass (signal centroid) of the cluster |
| `n_positions` | Number of distinct CPA positions in the cluster |
| `n_reads` | Total reads in the cluster (across all samples) |
| `cluster_width` | Cluster width in bp (`end - start + 1`) |

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

### `tables/shift_analysis_{condition}.tsv`

Genes with APA site-usage shifts between conditions. Columns include a
per-gene Jensen-Shannon divergence and proximal/distal usage deltas with an
adjusted p-value.

---

### `tables/go_enrichment_up_{condition}.tsv` / `tables/go_enrichment_down_{condition}.tsv`

GO term enrichment, split into up- and down-regulated gene sets.

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

`rectify export` writes strand-separated tracks, named per replicate and per
summed condition (e.g. `wt_by4742_rep1.plus.bw`, `wt_by4742.minus.bw`).

### `{name}.plus.bw` / `{name}.minus.bw`

BigWig format — per-base 3' end coverage, one file per strand. Load directly in IGV, the UCSC Genome Browser, or pyGenomeTracks.

### `{name}.plus.bedgraph` / `{name}.minus.bedgraph`

BedGraph format — same data as bigWig, plain text (written with `--format bedgraph`).

---

## Visualization outputs

### `plots/pca_samples.png`

PCA of samples by CPA cluster usage (first 2 PCs).

### `plots/sample_heatmap.png`

Sample-clustering heatmap of CPA cluster usage across samples (≥ 2 samples required).

### `report.html`

Comprehensive HTML report including QC stats, the per-aligner consensus
breakdown, PCA, heatmaps, top DE genes, and motif logos.

---

## Consensus selection stats

### `<prefix>.consensus_aligner_stats.tsv`

Written by `rectify align` and `rectify consensus`. Each row is a
`metric  count  percent  description` record summarizing aligner selection:

| Metric (row) | Description |
|--------------|-------------|
| `total_reads` | Total reads written to the consensus BAM |
| `aligner_wins_<aligner>` | Reads where `<aligner>` was selected |
| `confidence_high` / `confidence_medium` / `confidence_low` | Selection-confidence breakdown |
| `5prime_rescued` | Reads rescued via 5' soft-clip splice detection |
| `tied_score` | Reads requiring a tiebreaker (equal top junction score) |
| `aligner_combo_<panel>` | Reads whose available/compared aligner panel was `<panel>` (e.g. `gapmm2+minimap2`) |

The per-aligner-combination and tied-score breakdowns also surface in the
consensus section of the HTML report.
