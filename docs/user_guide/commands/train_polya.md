# rectify train-polya

Train a custom poly(A) tail model from control sites with zero downstream A's.

The built-in model was trained on *S. cerevisiae* WT data. Use this command to train a model for a different organism or a different sequencing protocol.

---

## Usage

```bash
rectify train-polya <input.bam> [options] -o <model.json>
```

## Examples

```bash
# Train from control CPA sites
rectify train-polya reads.bam \
    --genome genome.fa.gz \
    --control-sites control_cpa_sites.tsv \
    -o my_polya_model.json

# With custom minimum reads per site
rectify train-polya reads.bam \
    --genome genome.fa.gz \
    --control-sites control_cpa_sites.tsv \
    --min-reads 20 \
    -o my_polya_model.json
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `bam` | — | Input BAM file (aligned reads) |
| `--genome` | — | Reference genome FASTA |
| `--control-sites` | — | TSV file with control CPA sites (zero downstream A's) |
| `-o, --output` | — | Output model file (JSON) |
| `--min-reads` | 10 | Minimum reads per control site |

---

## Control sites file format

A tab-separated file with at minimum `chrom`, `strand`, and `pos` columns:

```tsv
chrom	strand	pos
chrI	+	34521
chrI	+	89201
chrII	-	12043
```

Control sites are positions where the genome sequence immediately downstream has zero A's (for plus strand) or zero T's (for minus strand). These are used to learn the poly(A) tail detection model without A-tract confounders.

---

## Using the trained model

Pass the model JSON to `rectify correct`:

```bash
rectify correct reads.bam \
    --genome genome.fa.gz \
    --polya-model my_polya_model.json \
    -o corrected.tsv
```

---

## What the model learns

The `PolyAModel` captures:

- **A-richness distribution** at poly(A) tail positions
- **Tail length distribution** (mean, variance by site type)
- **Position-dependent A-frequency profiles** (first 50 bp upstream/downstream of CPA)

These parameters inform the A-tract ambiguity detection and walk-back algorithm thresholds.
