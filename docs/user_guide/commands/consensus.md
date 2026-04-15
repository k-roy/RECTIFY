# rectify consensus

Run consensus aligner selection on pre-built per-aligner BAMs.

This command is the final step in the chunked parallel alignment workflow: after
`samtools merge` produces one sorted BAM per aligner, `rectify consensus` reads all
BAMs simultaneously, scores each aligner's alignment per read, and writes the
best-aligner-per-read result to a single rectified BAM.

It uses the same scoring logic as `rectify align`: 5' soft-clip quality, 3' A-tract
depth, canonical GT-AG splice sites, and annotated junction bonuses.

---

## Usage

```bash
rectify consensus ALIGNER:BAM [ALIGNER:BAM ...] \
    --genome <genome.fa> \
    -o <output_dir> \
    [--annotation <genes.gff>]
```

Each positional argument is an `aligner:path` pair. Accepted aligner names:
`minimap2`, `mapPacBio`, `gapmm2`, `uLTRA`, `deSALT`.

## Examples

```bash
# 5-aligner consensus after chunked alignment
rectify consensus \
    minimap2:merged/sample.minimap2.sorted.bam \
    mapPacBio:merged/sample.mapPacBio.sorted.bam \
    gapmm2:merged/sample.gapmm2.sorted.bam \
    uLTRA:merged/sample.uLTRA.sorted.bam \
    deSALT:merged/sample.deSALT.sorted.bam \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o results/sample/ \
    --prefix sample

# 3-aligner consensus (no junction-mode aligners)
rectify consensus \
    minimap2:sample.minimap2.sorted.bam \
    mapPacBio:sample.mapPacBio.sorted.bam \
    gapmm2:sample.gapmm2.sorted.bam \
    --genome genome.fa.gz \
    -o results/sample/
```

---

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `ALIGNER:BAM` | Yes (≥2) | Per-aligner BAMs as `name:path` pairs |
| `--genome` | Yes | Reference genome FASTA (for MD-tag recalculation) |
| `-o, --output-dir` | Yes | Output directory |
| `--prefix` | No | Output file prefix (default: derived from first BAM name) |
| `--annotation` | No | Gene annotation GFF/GTF for junction scoring |
| `--chimeric` | No | Use chimeric CIGAR assembly (experimental) |
| `--verbose` | No | Verbose logging |

---

## Output

```
output_dir/
└── <prefix>.consensus.bam      # coordinate-sorted, indexed, MD-tagged
└── <prefix>.consensus.bam.bai
```

The output BAM carries per-read tags:

| Tag | Type | Description |
|-----|------|-------------|
| `XA` | String | Winning aligner name (`minimap2`, `mapPacBio`, etc.) |
| `XC` | Integer | Confidence level (0=low, 1=medium, 2=high) |
| `XN` | Integer | Number of aligners agreeing with winner |
| `XR` | Integer | 5' rescue flag (1 = soft-clip rescued through splice junction) |

---

## Relationship to `rectify align`

`rectify align` performs alignment **and** consensus selection in one step.
`rectify consensus` performs consensus selection only, on BAMs you provide.

Use `rectify consensus` when:
- You ran alignment via `rectify split` + SLURM array and need to combine the results
- You have pre-existing per-aligner BAMs from an earlier run
- You want to re-run consensus scoring with different parameters without re-aligning

---

## Notes

- At least 2 aligner BAMs are required. For a single-aligner output, use `rectify align --aligners minimap2 --no-consensus` instead.
- All input BAMs must be coordinate-sorted and indexed (`.bai` file present).
- The genome path is used only for `samtools calmd` (MD-tag recalculation). It must be bgzip-compressed or uncompressed — a gzip-only file will be auto-converted.
