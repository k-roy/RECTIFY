# rectify extract

Extract per-read features from a BAM file to a TSV — 5' ends, 3' ends, splice junctions, alignment quality, poly(A) length, and soft-clip sequences.

Useful when you want per-read metadata without running the full correction pipeline.

---

## Usage

```bash
rectify extract <input.bam> [options] -o <output.tsv>
```

## Examples

```bash
# Basic feature extraction
rectify extract reads.bam \
    --genome genome.fa.gz \
    -o features.tsv

# With gene annotation (adds junction validation)
rectify extract reads.bam \
    --genome genome.fa.gz \
    --annotation genes.gff.gz \
    -o features.tsv

# Bundled yeast data
rectify extract reads.bam --Scer -o features.tsv
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `input` | — | Input BAM file |
| `-o, --output` | — | Output TSV file |
| `--genome` | — | Reference genome FASTA |
| `--annotation` | — | Gene annotation GFF/GTF |
| `--Scer` | — | Bundled *S. cerevisiae* data |
| `-j, --threads` | auto | Number of threads |

---

## Output columns

| Column | Description |
|--------|-------------|
| `read_id` | Read name |
| `chrom` | Chromosome |
| `strand` | `+` or `-` |
| `five_prime_raw` | Raw 5' end position (0-based) |
| `three_prime_raw` | Raw 3' end position (0-based) |
| `alignment_start` | Leftmost aligned position (0-based) |
| `alignment_end` | Rightmost aligned position + 1 (exclusive) |
| `n_junctions` | Number of splice junctions |
| `junctions` | Comma-separated list of `start:end` junction coordinates |
| `five_prime_soft_clip` | Soft-clipped bases at 5' end |
| `three_prime_soft_clip` | Soft-clipped bases at 3' end |
| `five_prime_soft_clip_seq` | Sequence of 5' soft-clipped bases |
| `three_prime_soft_clip_seq` | Sequence of 3' soft-clipped bases |
| `polya_length` | Total poly(A) length (aligned + soft-clipped) |
| `mapq` | Mapping quality |
| `alignment_identity` | Fraction of aligned positions matching reference |

---

## Notes

- Coordinate convention: all positions are 0-based, half-open (BED/pysam convention)
- For the strand-specific 5' and 3' end definitions, see [Coordinate System](../../coordinate_system.md)
- For corrected positions (after indel correction), use `rectify correct` instead
