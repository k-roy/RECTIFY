# rectify align

Multi-aligner consensus alignment for direct RNA and long-read sequencing.

Runs up to five aligners in parallel, attempts to rescue soft-clips through annotated splice junctions, scores alignments by canonical GT-AG sites and annotation matches, and writes the best alignment per read to a consensus BAM.

---

## Usage

```bash
rectify align <reads> [options] -o <output.bam>
```

## Examples

```bash
# Default 3-aligner consensus
rectify align reads.fastq.gz \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o aligned.bam

# Bundled yeast data
rectify align reads.fastq.gz --Scer -o aligned.bam

# 5-aligner consensus (add uLTRA + deSALT)
rectify align reads.fastq.gz \
    --Scer \
    --junction-aligners uLTRA deSALT \
    -o aligned.bam

# Single aligner (faster, less accurate)
rectify align reads.fastq.gz --Scer --aligner minimap2 -o aligned.bam
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `reads` | — | Input FASTQ or FASTQ.GZ file |
| `-o, --output` | — | Output BAM file |
| `--genome` | — | Reference genome FASTA |
| `--annotation` | — | Gene annotation GFF/GTF (for junction BED) |
| `--Scer` | — | Bundled *S. cerevisiae* data |
| `--aligner` | minimap2 | Primary aligner choice |
| `--junction-aligners` | — | Add optional junction-mode aligners: `uLTRA`, `deSALT` |
| `--chimeric-consensus` | off | Use chimeric CIGAR assembly (experimental) |
| `--ultra-path` | auto | Path to uLTRA executable |
| `--desalt-path` | auto | Path to deSALT executable |
| `-j, --threads` | auto | Threads for alignment |

---

## Aligners

### Default (always run)

| Aligner | Strengths |
|---------|-----------|
| **minimap2** | Fast, splice-aware; uses junction BED annotation for improved accuracy |
| **mapPacBio** (pbmm2) | PacBio RNA mode; forces mismatches at splice junctions for fair scoring |
| **gapmm2** | Gap-aware minimap2 variant; handles reads with large indels |

### Optional (add with `--junction-aligners`)

| Aligner | When to use |
|---------|-------------|
| **uLTRA** | Small exons (11–20 nt); annotation-guided collinear chaining |
| **deSALT** | Additional De Bruijn graph mapper; can resolve some complex junctions |

!!! note "Two-phase scheduler"
    mapPacBio is ~10× slower than other aligners. RECTIFY runs it alone with all threads first (phase 1), then runs the remaining aligners in parallel (phase 2). This prevents resource contention and is faster than running all aligners together.

---

## Consensus scoring

Per read, RECTIFY scores each aligner's output on:

1. **5' soft-clip penalty** — sequence-based: aligns clipped bases against upstream exon (edit distance ≤ 20% mismatches); rescued clips carry no penalty
2. **3' A-tract depth penalty** — penalizes aligners that land further into downstream A-tracts
3. **Canonical GT-AG splice sites** — more canonical junctions = higher score
4. **Annotated junction support** — bonus for junctions matching the provided annotation

The aligner with the highest composite score is written to the output BAM.

---

## minimap2 parameters

```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD
    -t <threads>
    --junc-bed annotation.junc.bed
    --junc-bonus 9
    genome.fa.gz reads.fastq.gz
```

Key flags:
- `-uf`: forward-strand only (correct for direct RNA)
- `-k14`: smaller k-mer for sensitivity on noisy nanopore reads
- `-G 5000`: max intron size (for yeast)
- `--splice-flank=no`: disables GT-AG bonus (important for 3' end accuracy)
- `--MD`: required for indel artifact correction downstream

---

## Notes

- The junction BED file is cached as `annotation.junc.bed` in the output directory on first run
- deSALT has a known output-duplication bug; RECTIFY deduplicates on (read_name, flag, chrom, pos, cigar) automatically
- uLTRA requires an uncompressed genome FASTA; RECTIFY decompresses gzipped genomes to a temp file automatically
- Dorado RNA004 basecalled reads use `U` instead of `T`; RECTIFY normalizes `U→T` before passing to mapPacBio
