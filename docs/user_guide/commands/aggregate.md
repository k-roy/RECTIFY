# rectify aggregate

Aggregate reads from a BAM file into position-level datasets for 3' ends, 5' ends, and splice junctions.

Produces count files suitable for genome browser visualization, metagene analysis, or downstream statistical analysis.

---

## Usage

```bash
rectify aggregate <input.bam> [options] -o <output_dir>
```

## Examples

```bash
# All three datasets
rectify aggregate reads.bam \
    --annotation genes.gff.gz \
    --mode all \
    -o aggregated/

# 3' ends only
rectify aggregate reads.bam \
    --annotation genes.gff.gz \
    --mode 3prime \
    -o aggregated/

# Bundled yeast data
rectify aggregate reads.bam --Scer --mode all -o aggregated/
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `bam` | — | Input BAM file (sorted and indexed) |
| `-o, --output-dir` | — | Output directory |
| `--genome` / `--annotation` | — | Reference genome FASTA / gene annotation GFF/GTF |
| `--Scer` / `--organism` | — | Use bundled *S. cerevisiae* data |
| `--mode` | `all` | What to aggregate: `3prime`, `5prime`, `junctions`, or `all` |
| `--format` | `tsv` | Output format: `tsv` or `parquet` |
| `--include-read-ids` | off | Include contributing read IDs (larger files) |
| `--prefix` | input stem | Output file-name prefix |

---

## Output files

Files are named `<prefix>.<format>` (the `<prefix>` defaults to the input BAM stem):

| File | Description |
|------|-------------|
| `<prefix>_3prime_clusters.{tsv,parquet}` | Per-position 3' end counts |
| `<prefix>_5prime_clusters.{tsv,parquet}` | Per-position 5' end counts |
| `<prefix>_junctions.{tsv,parquet}` | Per-junction read counts |

---

## Notes

- Use `rectify correct` + `rectify export` for corrected 3' end bigWig tracks
- This command produces raw (uncorrected) position counts from the BAM
- Junction counts from this command can feed into `rectify validate` for junction-level QC
