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
| `chrom` | Chromosome |
| `strand` | `+` or `-` |
| `original_3prime` | 3' end before correction |
| `corrected_3prime` | 3' end after correction (the position validated) |
| `ground_truth_position` | Nearest ground truth position |
| `ground_truth_source` | `netseq`, `annotation`, etc. |
| `distance_from_truth` | Corrected-position distance from ground truth (signed) |
| `original_distance` | Original-position distance from ground truth |
| `is_correct_exact` / `is_correct_1bp` / `is_correct_2bp` | Match flags at increasing tolerance |
| `improvement_bp` | Distance improvement from correction |
| `validation_status` | Per-read validation outcome |
| `correction_applied` | Whether a correction was applied |
| `confidence` | `HIGH`, `MEDIUM`, `LOW` |

### Summary statistics (printed to stdout)

The summary reports corrected-position accuracy, the original (pre-correction)
accuracy for comparison, and the net improvement — for example:

```text
Corrected Position Accuracy:
  Exact match:            ...
  Within 1 bp:            ...
  Within 2 bp:            ...
  Within 5 bp:            ...

Original Position Accuracy (Before Correction):
  Exact match:            ...
  Within 1 bp:            ...
  Within 2 bp:            ...

Improvement Over Original:
  Improved:               ...
  Unchanged:              ...
  Worsened:               ...
  Net ±1bp improvement:   ...
```

---

## Notes

- Use `--Scer` to validate against bundled yeast WT NET-seq without providing any files
- Reads with `qc_flags = AG_RICH` are reported separately (lower expected accuracy at AG-rich sites)
- The `--tolerance` parameter sets what counts as "correct" for accuracy metrics; 1 bp is standard for single-nucleotide resolution claims
