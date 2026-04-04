# rectify correct

Correct 3' (and 5') end positions in a BAM or FASTQ file.

This is the core RECTIFY command. It applies the walk-back indel correction algorithm, A-tract ambiguity detection, AG-mispriming screening, and optionally NET-seq refinement.

---

## Usage

```bash
rectify correct <input> [options] -o <output>
```

## Examples

```bash
# Bundled yeast data — no external files needed
rectify correct reads.bam --Scer -o corrected.tsv

# Custom genome + annotation
rectify correct reads.bam \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o corrected.tsv

# With NET-seq refinement
rectify correct reads.bam --Scer --netseq-dir my_netseq/ -o corrected.tsv

# Streaming mode for large BAMs (>2 GB)
rectify correct large.bam --Scer --streaming -o corrected.tsv

# Remove spike-in reads
rectify correct reads.bam --Scer --filter-spikein ENO2 -o corrected.tsv
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `input` | Input BAM or FASTQ/FASTQ.GZ file |
| `-o, --output` | Output TSV file path |

### Reference data

| Argument | Description |
|----------|-------------|
| `--genome` | Reference genome FASTA (optionally gzipped) |
| `--annotation` | Gene annotation file (GTF or GFF, optionally gzipped) |
| `--Scer` | Use bundled *S. cerevisiae* S288C data (genome + annotation + NET-seq) |
| `--organism` | Organism name for bundled data (e.g. `yeast`) |

### NET-seq refinement

| Argument | Default | Description |
|----------|---------|-------------|
| `--netseq-dir` | — | Directory of NET-seq BigWig files for A-tract refinement |
| `--netseq-samples` | all | Specific NET-seq samples to use |

### Poly(A) handling

| Argument | Default | Description |
|----------|---------|-------------|
| `--polya-sequenced` | off | Enable poly(A) trimming and indel correction (use for direct RNA / QuantSeq) |
| `--polya-model` | built-in | Pre-trained poly(A) tail model (JSON from `rectify train-polya`) |

### Module selection

| Argument | Description |
|----------|-------------|
| `--skip-atract-check` | Skip A-tract ambiguity detection |
| `--skip-ag-check` | Skip AG-mispriming screening |
| `--skip-polya-trim` | Skip poly(A) tail trimming |
| `--skip-indel-correction` | Skip indel artifact correction |
| `--skip-variant-aware` | Skip variant-aware position rescue |

### Filtering

| Argument | Description |
|----------|-------------|
| `--filter-spikein GENE [GENE ...]` | Remove reads from spike-in genes by name |

### Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--ag-threshold` | 0.65 | AG-richness threshold (0–1) for mispriming flag |

### Output options

| Argument | Description |
|----------|-------------|
| `--report` | Write QC report to this path (`.html` or `.pdf`) |

### Performance

| Argument | Default | Description |
|----------|---------|-------------|
| `-j, --threads` | auto | Number of processing threads (0 = auto-detect) |
| `--streaming` | off | Streaming output mode — keeps peak RAM ~4–5 GB for any BAM size |
| `--chunk-size` | 10000 | Reads per output chunk (streaming mode only) |

---

## Output files

| File | Description |
|------|-------------|
| `<output>.tsv` | Per-read corrected positions (see [Output Formats](../output_formats.md)) |
| `<output>_index.bed.gz` | Pre-aggregated position counts (~300× smaller) |
| `<output>_alignment_features.tsv` | Per-read alignment metadata |
| `<output>_stats.tsv` | Processing QC summary |

---

## Notes

- For FASTQ input, pre-align first with `rectify align`, then pass the BAM to `rectify correct`
- To run alignment + correction in one step, use `rectify run`
- `--streaming` is recommended for BAMs larger than 2 GB; it is the default in the bundled SLURM profiles
- The output index file (`_index.bed.gz`) is used by manifest-mode analysis; generate it for all samples before running `rectify analyze --manifest`
