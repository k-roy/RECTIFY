# rectify netseq

Process NET-seq BAM files: extract 3' end positions, apply oligo(A)-spreading deconvolution, and write strand-specific BigWig files.

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

---

## Output files

| File | Description |
|------|-------------|
| `netseq_plus.bw` | Plus strand 3' end signal (BigWig) |
| `netseq_minus.bw` | Minus strand 3' end signal (BigWig) |
| `netseq_plus_deconvolved.bw` | Deconvolved plus strand signal |
| `netseq_minus_deconvolved.bw` | Deconvolved minus strand signal |
| `netseq_stats.tsv` | Processing statistics |

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
