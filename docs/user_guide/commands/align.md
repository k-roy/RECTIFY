# rectify align

Multi-aligner consensus alignment for direct RNA and long-read sequencing.

By default it runs the three-aligner long-read panel (minimap2, mapPacBio, gapmm2); `--junction-aligners` adds the opt-in splice-aware aligners uLTRA and deSALT for the full five-aligner panel. It attempts to rescue soft-clips through annotated splice junctions, scores alignments by canonical GT-AG sites and annotation matches, and writes the best alignment per read to a rectified BAM.

---

## Usage

```bash
rectify align <reads> [options] -o <output_dir>
```

The output is a **directory**: per-aligner BAMs are written as
`<prefix>.<aligner>.bam`, and the consensus result as
`<prefix>.rectified.bam` (the `<prefix>` defaults to the input file stem; set
it with `--prefix`).

## Examples

```bash
# Default 3-aligner consensus (minimap2 + mapPacBio + gapmm2)
rectify align reads.fastq.gz \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o aligned/

# Bundled yeast data
rectify align reads.fastq.gz --Scer -o aligned/

# 5-aligner consensus (add uLTRA + deSALT)
rectify align reads.fastq.gz \
    --Scer \
    --junction-aligners uLTRA deSALT \
    -o aligned/

# Single aligner (faster, less accurate)
rectify align reads.fastq.gz --Scer --aligners minimap2 -o aligned/
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `reads` | — | Input FASTQ or FASTQ.GZ file |
| `-o, --output-dir` | — | Output directory for BAM files |
| `--prefix` | input stem | Output file-name prefix |
| `--genome` | — | Reference genome FASTA |
| `--annotation` | — | Gene annotation GFF/GTF (for junction BED) |
| `--Scer` | — | Bundled *S. cerevisiae* data |
| `--aligners` | all | Base aligners to run (`all` = minimap2 + mapPacBio + gapmm2 for long-read); pass a subset (e.g. `minimap2`) or `none` |
| `--junction-aligners` | (none) | Add opt-in splice-aware aligners: `uLTRA`, `deSALT`, `gmap` |
| `--chimeric-consensus` | off | Use chimeric CIGAR assembly (experimental) |
| `--ultra-path` | `uLTRA` | Path to uLTRA executable |
| `--desalt-path` | `deSALT` | Path to deSALT executable |
| `--parallel-aligners` | off | Run base aligners in parallel (phase 2) |
| `-t, --threads` | 8 | Threads per aligner |

---

## Aligners

### Default base panel (long-read)

| Aligner | Strengths |
|---------|-----------|
| **minimap2** | Fast, splice-aware; uses junction BED annotation for improved accuracy |
| **mapPacBio** (BBMap's `mapPacBio.sh`) | PacBio long-read RNA mode |
| **gapmm2** | Gap-aware minimap2 wrapper (pinned to `==25.4.5`); handles reads with large indels |

### Optional (add with `--junction-aligners`)

| Aligner | When to use |
|---------|-------------|
| **uLTRA** | Small exons; annotation-guided collinear chaining (requires `--annotation`) |
| **deSALT** | Additional reference-graph mapper; can resolve some complex junctions |
| **gmap** | Splice-aware mapper; requires a pre-built database (see `--gmap-db`) |

Benchmark any opt-in junction aligner against your data before relying on it in production.

!!! note "Two-phase scheduler"
    mapPacBio is ~10× slower than the other aligners. With `--parallel-aligners`, RECTIFY runs mapPacBio alone with all threads first (phase 1), then runs the remaining aligners in parallel (phase 2). This avoids resource contention and is faster than running every aligner at once.

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
- `--splice-flank=no`: disables the GT-AG splice-flank model (set for compatibility / 3' end accuracy)
- `--MD`: required for indel artifact correction downstream

---

## Parallel alignment (large datasets)

For datasets where a single 5-aligner run exceeds available wall time, use
`rectify split` to divide the FASTQ into N chunks and run each (chunk, aligner)
pair as a SLURM array task. See [`rectify split`](split.md) for the complete workflow.

```bash
# Quick summary. mapPacBio gets its own array script; the remaining
# aligners (default --other-aligners = minimap2 gapmm2 uLTRA deSALT) share
# the "others" array.
rectify split reads.fastq.gz -n 16 -o chunks/ \
    --generate-slurm \
    --genome genome.fa --annotation genes.gff
bash chunks/submit_pipeline.sh   # submits the full DAG with dependencies
```

---

## Notes

- The junction BED file is cached as `annotation.junc.bed` in the output directory on first run
- deSALT has a known output-duplication bug; RECTIFY deduplicates on (read_name, flag, chrom, pos, cigar) automatically
- deSALT: a pre-built Linux/x86_64 binary (v1.5.6) is bundled with RECTIFY and used automatically when `deSALT` is not on `PATH`. Use `rectify install-aligners --desalt` to install to `~/.rectify/bin/` for other platforms.
- uLTRA requires an uncompressed genome FASTA; RECTIFY decompresses gzipped genomes to a temp file automatically
- Dorado RNA004 basecalled reads use `U` instead of `T`; RECTIFY normalizes `U→T` before passing to mapPacBio
