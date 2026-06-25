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
| `bam` | — | Input BAM file |
| `-o, --output` | — | Output TSV file |
| `--genome` | — | Reference genome FASTA |
| `--annotation` | — | Gene annotation GFF/GTF |
| `--Scer` / `--organism` | — | Use bundled *S. cerevisiae* data |
| `--include-sequence` | off | Add the read `sequence` column |
| `--include-quality` | off | Add the per-base `quality` column |
| `--include-junctions` | on | Add `n_junctions` / `junctions` / `is_spliced` |
| `--include-softclips` | on | Add soft-clip length columns |
| `--include-context` | off | Add downstream genomic context (requires `--genome`; with `--context-length`) |
| `--spliced-only` | off | Emit only spliced reads |
| `--streaming` | off | Stream output to file (lower memory) |
| `--chunk-size` | — | Reads per chunk in streaming mode |

---

## Output columns

Always present:

| Column | Description |
|--------|-------------|
| `read_id` | Read name |
| `chrom` | Chromosome |
| `strand` | `+` or `-` |
| `reference_start` | Leftmost aligned position (0-based) |
| `reference_end` | Rightmost aligned position + 1 (exclusive) |
| `five_prime` | 5' end position (0-based, strand-aware) |
| `three_prime` | 3' end position (0-based, strand-aware) |
| `aligned_length` | Reference bases consumed by the alignment |
| `mapq` | Mapping quality |

Optional (added by the matching `--include-*` flag):

| Column | Added by |
|--------|----------|
| `n_junctions`, `junctions`, `is_spliced` | `--include-junctions` (on by default) |
| `left_softclip`, `right_softclip`, `five_prime_softclip`, `three_prime_softclip` | `--include-softclips` (on by default) |
| `downstream_context` | `--include-context` |
| `sequence` | `--include-sequence` |
| `quality` | `--include-quality` |

---

## Notes

- Coordinate convention: all positions are 0-based, half-open (BED/pysam convention)
- For the strand-specific 5' and 3' end definitions, see [Coordinate System](../../coordinate_system.md)
- For corrected positions (after indel correction), use `rectify correct` instead
