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
| `input` | — | Input BAM file |
| `-o, --output-dir` | — | Output directory |
| `--annotation` | — | Gene annotation GFF/GTF |
| `--Scer` | — | Bundled *S. cerevisiae* data |
| `--mode` | `3prime` | What to aggregate: `3prime`, `5prime`, `junctions`, or `all` |
| `-j, --threads` | auto | Number of threads |

---

## Output files

| File | Description |
|------|-------------|
| `3prime_ends.tsv` | Per-position 3' end counts (chrom, pos, strand, count) |
| `5prime_ends.tsv` | Per-position 5' end counts |
| `junctions.tsv` | Per-junction read counts (donor_pos, acceptor_pos, strand, count) |

---

## Notes

- Use `rectify correct` + `rectify export` for corrected 3' end bigWig tracks
- This command produces raw (uncorrected) position counts from the BAM
- Junction counts from this command can feed into `rectify validate` for junction-level QC
