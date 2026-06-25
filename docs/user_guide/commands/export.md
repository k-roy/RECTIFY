# rectify export

Export corrected 3' end positions to bigWig or bedGraph format for genome browser visualization.

---

## Usage

```bash
rectify export <input.tsv> [options] -o <output_dir>
```

## Examples

```bash
# BigWig output (default)
rectify export corrected_reads.tsv \
    --genome genome.fa.gz \
    -o tracks/

# BedGraph output
rectify export corrected_reads.tsv \
    --genome genome.fa.gz \
    --format bedgraph \
    -o tracks/

# Using a chromosome sizes file
rectify export corrected_reads.tsv \
    --chrom-sizes chrom.sizes \
    -o tracks/

# Per-replicate and per-condition tracks
rectify export corrected_reads.tsv \
    --genome genome.fa.gz \
    --per-replicate \
    --per-condition \
    -o tracks/

# Bundled yeast genome (auto-fills --genome / chromosome sizes)
rectify export corrected_reads.tsv --Scer -o tracks/
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `input` | — | Corrected 3' end TSV file |
| `-o, --output-dir` | — | Output directory |
| `--format` | `bigwig` | Output format: `bigwig` or `bedgraph` |
| `--genome` | — | Reference genome FASTA (for chromosome sizes) |
| `--chrom-sizes` | — | Chromosome sizes file (alternative to `--genome`) |
| `--position-col` | `position` | Column name containing the corrected position |
| `--per-replicate` | off | Write one file per replicate |
| `--per-condition` | off | Write one summed file per condition |

---

## Output files

Tracks are strand-separated, named per replicate and/or per summed condition
(e.g. `wt_by4742_rep1.plus.bw`, `wt_by4742.minus.bw`):

| File | Description |
|------|-------------|
| `{name}.plus.bw` / `{name}.minus.bw` | BigWig — per-base 3' end coverage, one file per strand (`--format bigwig`, default) |
| `{name}.plus.bedgraph` / `{name}.minus.bedgraph` | BedGraph — plain text equivalent (`--format bedgraph`) |

---

## Loading in genome browsers

**IGV:**

1. File → Load from File → select `*.bw`
2. Right-click track → Set Data Range

**UCSC Genome Browser:**

```
track type=bigWig name="RECTIFY 3' ends" bigDataUrl=https://your-server/sample.bw
```

**pyGenomeTracks:**

```ini
[rectify_3prime]
file = sample.bw
type = bigwig
color = steelblue
```

---

## Notes

- The input TSV must contain `chrom`, `strand`, and the position column named by
  `--position-col` (default `position`). To export from `corrected_reads.tsv`,
  point at its corrected-3'-end column: `--position-col corrected_3prime`.
- Positions in the output represent 3' end (CPA) positions, not read coverage
- For reads on the minus strand, the position is at `reference_start` (leftmost coordinate, which is the 3' end)
- Generate tracks per sample first, then optionally use `bedGraphToBigWig` to merge conditions
