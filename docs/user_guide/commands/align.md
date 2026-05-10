# rectify align

Multi-aligner consensus alignment for direct RNA and long-read sequencing.

Runs up to five aligners in parallel, attempts to rescue soft-clips through annotated splice junctions, scores alignments by canonical GT-AG sites and annotation matches, and writes the best alignment per read to a rectified BAM.

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

Per read, RECTIFY scores each aligner's output on five signals:

1. **Effective 5' clip penalty** — −2 per unrescued bp of 5' soft-clip
2. **5' rescue attempt** — clipped bases are aligned (NW) against the upstream exon; a successful rescue cancels the 5' clip penalty
3. **A-tract 3' depth** — −1 per A extending into a downstream A-tract, capped at 10
4. **3' non-poly(A) terminal errors** — −2 per non-poly(A) base error at the 3' end, capped at 10
5. **Junction-proximity errors** — −1 per error within a window around splice junctions, capped at 10

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

## Parallel alignment (large datasets)

For datasets where a single 5-aligner run exceeds available wall time, use
`rectify split` to divide the FASTQ into N chunks and run each (chunk, aligner)
pair as a SLURM array task. See [`rectify split`](split.md) for the complete workflow.

```bash
# Quick summary: 16 chunks × 5 aligners = 80 array tasks
rectify split reads.fastq.gz -n 16 -o chunks/ \
    --generate-slurm \
    --aligners minimap2 mapPacBio gapmm2 uLTRA deSALT \
    --genome genome.fa --annotation genes.gff
sbatch chunks/run_array_align.sh          # submit 80-task array
bash  chunks/run_merge_and_consensus.sh   # run after array completes
```

---

## Notes

- The junction BED file is cached as `annotation.junc.bed` in the output directory on first run
- deSALT has a known output-duplication bug; RECTIFY deduplicates on (read_name, flag, chrom, pos, cigar) automatically
- deSALT: a pre-built Linux/x86_64 binary (v1.5.6) is bundled with RECTIFY and used automatically when `deSALT` is not on `PATH`. Use `rectify install-aligners --desalt` to install to `~/.rectify/bin/` for other platforms.
- uLTRA requires an uncompressed genome FASTA; RECTIFY decompresses gzipped genomes to a temp file automatically
- Dorado RNA004 basecalled reads use `U` instead of `T`; RECTIFY normalizes `U→T` before passing to mapPacBio
