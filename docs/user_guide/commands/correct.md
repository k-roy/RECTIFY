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

### Technology / protocol

| Argument | Default | Description |
|----------|---------|-------------|
| `--dT-primed-cDNA` | off | Input was generated with oligo-dT priming (QuantSeq, dT-primed cDNA). Poly(A) tail is NOT in the read. Enables indel artifact correction and poly(A) trimming modules. Do NOT use for ONT direct RNA. |
| `--ONT-cDNA` | off | Input is Oxford Nanopore PCR-cDNA (e.g. SQK-PCB114). Poly(A) tail IS present as a 3' soft-clip. Enables poly(A) trimming + indel correction; disables AG-mispriming. Do NOT combine with `--dT-primed-cDNA`. |
| `--short-read` | off | Input is short-read data (Illumina/Aviti ≤150 bp). Disables poly(A)-tail trimming, A-tract correction, and indel modules. Pair with `rectify align --short-read` (bbmap + bwa). |

*(deprecated aliases retained for backwards compatibility: `--polya-sequenced` / `--no-polya-sequenced`. Use the current flags above for new pipelines.)*

### Poly(A) handling

| Argument | Default | Description |
|----------|---------|-------------|
| `--polya-model` | built-in | Pre-trained poly(A) tail model (JSON from `rectify train-polya`) |
| `--min-polya-score` | — | Minimum poly(A) model confidence score (0–1) to mark `polya_pass=1`; reads below are flagged but still written. Requires `--polya-model` or `--dT-primed-cDNA`. |

### Module selection

| Argument | Description |
|----------|-------------|
| `--skip-atract-check` | Skip A-tract ambiguity detection |
| `--skip-ag-check` | Skip AG-mispriming screening |
| `--skip-polya-trim` | Skip poly(A) tail trimming |
| `--skip-indel-correction` | Skip indel artifact correction |
| `--skip-variant-aware` | Skip variant-aware homopolymer rescue (two-pass) |
| `--min-mapq` | Minimum MAPQ to include a read (default 0; use 1 to drop multi-mappers) |
| `--min-aligned-length` | Minimum reference bases consumed by the alignment (default 0; 30 recommended for dT-primed short-read) |

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
| `--output-bam` | Write a poly(A)-trimmed BAM alongside the corrected TSV (requires `--dT-primed-cDNA`) |
| `--write-corrected-bam` | Write a corrected BAM hard-clipped at every read's corrected 3' end (sorted + indexed) |
| `--write-softclipped-bam` | Like `--write-corrected-bam` but poly(A) bases are soft-clipped (visible in IGV) |
| `--write-bedgraph PREFIX` | Write strand-specific bedGraph files (`PREFIX.plus.bedgraph` / `PREFIX.minus.bedgraph`) for NET-seq fractional pileups |
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
- To run alignment + correction in one step, use `rectify run-all`
- `rectify correct` runs post-consensus on the winning aligner's BAM; use `--aligner-bams` to supply per-aligner BAMs as a junction candidate pool for Module 2H (post-consensus N-op refinement)
- `--streaming` is recommended for BAMs larger than 2 GB; it is the default in the bundled SLURM profiles
- The output index file (`_index.bed.gz`) is used by manifest-mode analysis; generate it for all samples before running `rectify analyze --manifest`
