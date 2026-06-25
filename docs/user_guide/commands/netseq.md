# rectify netseq

Optional refinement step. Process NET-seq BAM files: extract 3' end positions, apply oligo(A)-spreading deconvolution, and write strand-specific BigWig files for use with `--netseq-dir` in `rectify correct`.

---

## Usage

```bash
rectify netseq <input.bam> [options] -o <output_dir>
```

## Examples

```bash
# Basic NET-seq processing
rectify netseq netseq.bam \
    --genome genome.fa.gz \
    --gff genes.gff.gz \
    -o netseq_output/

# With exclusion regions
rectify netseq netseq.bam \
    --genome genome.fa.gz \
    --gff genes.gff.gz \
    --exclude-mito \
    -o netseq_output/

# Include rDNA and Pol III genes
rectify netseq netseq.bam \
    --Scer \
    --include-rdna \
    --include-pol3 \
    -o netseq_output/
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `input` | Input NET-seq BAM file |
| `-o, --output-dir` | Output directory |

### Reference data

| Argument | Description |
|----------|-------------|
| `--genome` | Reference genome FASTA |
| `--gff` | Gene annotation GFF file |
| `--Scer` | Bundled *S. cerevisiae* data |

### Exclusion regions

| Argument | Default | Description |
|----------|---------|-------------|
| `--exclude-mito` | off | Exclude mitochondrial reads |
| `--include-rdna` | off | Include rDNA loci (excluded by default) |
| `--include-pol3` | off | Include Pol III transcripts (excluded by default) |
| `--pol3-flanking` | 100 | Flanking bp around Pol III genes to exclude |

### Deconvolution

| Argument | Default | Description |
|----------|---------|-------------|
| `--no-deconvolution` | off | Disable NNLS deconvolution (output raw 3' ends only) |
| `--min-atract-length` | 3 | Minimum downstream A's to trigger deconvolution |

### Output options

| Argument | Default | Description |
|----------|---------|-------------|
| `--output-format` | parquet bedgraph | One or more of `parquet`, `bedgraph`, `bigwig`, `tsv` |
| `--no-rpm-normalize` | off | Disable RPM normalization for bedgraph/bigwig output |

### Processing

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-mapq` | 0 | Minimum mapping quality |
| `--max-reads` | — | Cap total reads processed (for testing) |
| `-v, --verbose` | off | Verbose output |

---

## Output files

Defaults: per-read records as parquet plus per-position signal as bedgraph (one strand-specific file each, raw and deconvolved). Use `--output-format bigwig` to produce BigWig files in addition (requires `pyBigWig`); `--output-format tsv` writes a plain-text equivalent of the parquet.

| File | Description |
|------|-------------|
| `{sample}.unified_reads.parquet` | Per-read 3' end records (default) |
| `{sample}.raw.plus.bedgraph` | Plus strand 3' end signal (default; RPM-normalized unless `--no-rpm-normalize`) |
| `{sample}.raw.minus.bedgraph` | Minus strand 3' end signal (default) |
| `{sample}.deconv.plus.bedgraph` | NNLS-deconvolved plus strand signal (default) |
| `{sample}.deconv.minus.bedgraph` | NNLS-deconvolved minus strand signal (default) |
| `{sample}.{raw,deconv}.{plus,minus}.bw` | BigWig variants of the above (only with `--output-format bigwig`) |

---

## Oligo(A) deconvolution

NET-seq 3' ends are oligo-adenylated, creating downstream signal spreading. RECTIFY applies NNLS deconvolution using a point-spread function (PSF) derived from 5,000+ zero-A control sites to recover true CPA positions.

The deconvolved BigWigs are used by `rectify correct --netseq-dir` for A-tract ambiguity resolution in direct RNA experiments.

See [NET-seq Refinement](../../algorithms/netseq_refinement.md) for details on the deconvolution algorithm.

---

## Notes

- rDNA and Pol III genes are excluded by default because their extremely high signal density distorts the deconvolution PSF
- The PSF is derived from positions with zero downstream A's — sites where the true CPA is unambiguous
- Use `--include-rdna` and `--include-pol3` only if you specifically want to analyze those loci
