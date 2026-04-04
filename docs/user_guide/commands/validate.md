# rectify validate

Validate corrected 3' end positions against ground truth data — NET-seq signal, gene annotations, or known positions.

---

## Usage

```bash
rectify validate <corrected.tsv> [options] -o <output.tsv>
```

## Examples

```bash
# Validate against bundled WT NET-seq
rectify validate corrected.tsv --Scer -o validation.tsv

# Validate against custom NET-seq
rectify validate corrected.tsv \
    --netseq-dir netseq_bigwigs/ \
    -o validation.tsv

# Validate against known positions file
rectify validate corrected.tsv \
    --ground-truth known_cpa_sites.tsv \
    -o validation.tsv

# Multiple ground truth sources
rectify validate corrected.tsv \
    --netseq-dir netseq/ \
    --annotation genes.gff.gz \
    --ground-truth known_sites.tsv \
    -o validation.tsv
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `corrected` | Corrected 3' end TSV (from `rectify correct`) |
| `-o, --output` | Output validation results TSV |

### Ground truth (at least one required)

| Argument | Description |
|----------|-------------|
| `--netseq-dir` | NET-seq BigWig directory |
| `--netseq-samples` | Specific NET-seq samples to use |
| `--annotation` | Gene annotation GTF/GFF (validates against annotated 3' UTR ends) |
| `--ground-truth` | TSV file with known CPA positions (columns: chrom, strand, pos) |
| `--Scer` | Bundled *S. cerevisiae* data |

### Validation parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--tolerance` | 1 | Position tolerance in bp for "exact match" scoring |
| `--min-signal` | 0.5 | Minimum NET-seq signal to consider a position supported |
| `--search-window` | 10 | Window (±bp) for finding nearest ground truth position |

---

## Output

### `validation.tsv`

Per-read validation results:

| Column | Description |
|--------|-------------|
| `read_id` | Read name |
| `corrected_position` | Position being validated |
| `nearest_ground_truth` | Closest ground truth position |
| `distance_bp` | Distance from ground truth (signed) |
| `netseq_support` | NET-seq signal at corrected position |
| `validated` | Boolean — within tolerance |
| `confidence` | `HIGH`, `MEDIUM`, `LOW` |

### Summary statistics (printed to stdout)

```
Validation summary:
  Total reads: 45,231
  Exact match (±1 bp):   88.3%
  Match at ±2 bp:        91.7%
  Match at ±5 bp:        94.2%
  Mean distance: -0.3 bp (negative = upstream of ground truth)

By correction type:
  No correction needed:  52.1%  — 99.1% at ±1 bp
  Shifted 1–5 bp:        31.4%  — 95.3% at ±1 bp
  Shifted 6–15 bp:       12.8%  — 78.2% at ±1 bp
  Shifted >15 bp:         3.7%  — 61.0% at ±1 bp
```

---

## Notes

- Use `--Scer` to validate against bundled yeast WT NET-seq without providing any files
- Reads with `qc_flags = AG_RICH` are reported separately (lower expected accuracy at AG-rich sites)
- The `--tolerance` parameter sets what counts as "correct" for accuracy metrics; 1 bp is standard for single-nucleotide resolution claims
