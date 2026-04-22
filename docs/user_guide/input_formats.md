# Input Formats

## Reads

RECTIFY accepts reads in two formats:

### BAM

Pre-aligned BAM files. Used directly by `rectify correct` (no re-alignment). Must be sorted and indexed.

```bash
rectify correct reads.bam --genome genome.fa -o corrected.tsv
```

### FASTQ / FASTQ.GZ

Unaligned reads. `rectify run` and `rectify align` will align these using the multi-aligner consensus pipeline before correction.

```bash
rectify run reads.fastq.gz --genome genome.fa --annotation genes.gff -o results/
```

!!! tip
    Pass `--skip-alignment` to `rectify run` if you have a FASTQ but want to skip alignment (unusual).

---

## Reference genome (FASTA)

Standard FASTA format, optionally gzip-compressed:

```
--genome genome.fa
--genome genome.fa.gz
```

For *S. cerevisiae*, skip this entirely with `--Scer` — the S288C genome is bundled.

---

## Gene annotation (GFF/GTF)

Used for:
- Junction BED generation for minimap2 (`--junc-bed`)
- AG-mispriming context
- Gene-level attribution in `rectify analyze`
- DESeq2 gene-level count matrix

Supported formats:

```
--annotation genes.gff
--annotation genes.gff.gz
--annotation genes.gtf
--annotation genes.gtf.gz
```

!!! note
    RECTIFY auto-detects format from the file extension. For uLTRA (optional 5th aligner), `.gff` files are auto-converted to `.gtf` via `gffread`.

---

## NET-seq data (BigWig)

Optional. Provides nascent RNA 3' end positions for NNLS deconvolution within A-tracts.

```
--netseq-dir /path/to/netseq/bigwigs/
```

The directory should contain one or more `.bw` files:

```
netseq_dir/
├── wt_plus.bw
├── wt_minus.bw
└── ...
```

For *S. cerevisiae* with `--Scer`, WT NET-seq is bundled and auto-detected.

---

## Manifest TSV (multi-sample)

A tab-separated file with three required columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier (used for output directory names) |
| `path` | Absolute path to input BAM or FASTQ.GZ file |
| `condition` | Condition label (e.g. `wt`, `ko`, `heat_shock`) |

Example:

```tsv
sample_id	path	condition
wt_rep1	/data/wt_rep1.fastq.gz	wt
wt_rep2	/data/wt_rep2.fastq.gz	wt
ko_rep1	/data/ko_rep1.fastq.gz	ko
ko_rep2	/data/ko_rep2.fastq.gz	ko
```

!!! tip
    The `--reference` flag to `rectify run` / `rectify analyze` is matched **case-insensitively** against the `condition` column, so `--reference wt` will match `WT`, `Wt`, `wt`.

---

## Poly(A) model (JSON)

A trained poly(A) tail model produced by `rectify train-polya`. If omitted, the built-in model is used.

```
--polya-model model.json
```

---

## GO annotations (TSV)

For GO enrichment analysis in `rectify analyze`. Tab-separated with `gene_id` and `go_term` columns. For yeast, bundled with `--Scer`.

```
--go-annotations go_annotations.tsv.gz
```

---

## SLURM profile (YAML)

Configuration for SLURM job submission via `rectify run --profile` or `rectify batch --profile`. See the [HPC / SLURM guide](hpc_slurm.md) for the full profile format.

```
--profile slurm_profiles/hpc_cpu.yaml
```
