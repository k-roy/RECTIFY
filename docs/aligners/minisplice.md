# Minisplice (minisplice_mm2): deep-learning splice-site scoring for minimap2

Minisplice trains a neural network model on splice-site sequence context and uses
its predictions to score candidate splice junctions, which are then passed to
minimap2 via `--spsc` to bias alignment toward biologically plausible junctions.

**Benchmark result (human chr5, 10k reads, 2026-05-25):** 96.0% GT-AG vs
minimap2's 91.6% baseline, at 1.8× minimap2's wall time (16 s vs 9 s for 10k
reads). The improvement is consistent and meaningful for splice-junction analysis.

See [ALIGNER_RECOMMENDATIONS.md](../ALIGNER_RECOMMENDATIONS.md) for the full
benchmark table.

---

## Installation

```bash
conda install -c bioconda minisplice
```

Minisplice installs:
- `minisplice` — the predict binary
- A patched `minimap2` (or the `--spsc` flag may be supported in minimap2 ≥2.29)

Download the pre-trained model (e.g., `vi2-7k.kan` for ONT direct RNA):
see the [minisplice releases page](https://github.com/comprna/minisplice/releases).
A calibration file (`vi2-7k.kan.cali`) is used with `-c` for better-calibrated
probability scores; it is optional.

---

## How RECTIFY runs it

`run_minisplice_mm2()` in `multi_aligner.py`:

1. **Predicts splice site scores** (one-time per genome; cached):
   ```bash
   minisplice predict -t <threads> [-c <model.cali>] <model> <genome> \
       > <genome_stem>_splice_scores.tsv
   ```
   The scores TSV is cached in `.minisplice_cache/` adjacent to the genome (or
   in `cache_dir`). Subsequent runs for the same genome skip this step.

2. **Aligns** using minimap2 with `--spsc`:
   ```bash
   minimap2 --spsc=<splice_scores.tsv> -ax splice -uf -k14 \
       -G <max_intron> --splice-flank=no --secondary=no --MD -y \
       -t <threads> <genome> <reads> | samtools sort -n ...
   ```

3. Applies `samtools calmd -e`, validates QNAMEs, and injects `RN` tags.

---

## CLI usage

```bash
rectify align <reads.fastq.gz> --genome <genome.fa> \
    --aligners minisplice_mm2 \
    --minisplice-model /path/to/vi2-7k.kan \
    --minisplice-model-cali /path/to/vi2-7k.kan.cali \
    --max-intron 500000 \
    -o /tmp/ms_out -t 8
```

If splice scores are pre-computed (e.g., shared across samples for the same genome):

```bash
rectify align <reads.fastq.gz> --genome <genome.fa> \
    --aligners minisplice_mm2 \
    --minisplice-scores /path/to/splice_scores.tsv \
    -o /tmp/ms_out -t 8
```

`minisplice_mm2` is not in the default `"all"` aligner set — it must be requested
explicitly. To combine with the default panel:

```bash
rectify align <reads.fastq.gz> --genome <genome.fa> \
    --aligners minimap2 mapPacBio gapmm2 minisplice_mm2 \
    --minisplice-scores /path/to/genome_splice_scores.tsv \
    --junction-aligners uLTRA deSALT \
    -o /tmp/out -t 8
```

---

## Performance notes

- Predict step: ~10–30 min for yeast; ~60–120 min for GRCh38 (cached after first run;
  amortized across all samples aligned to the same genome).
- Alignment: ~1.8× minimap2 wall time (16 s vs 9 s for 10k chr5 reads).
- Memory: similar to minimap2.

Pre-computing splice scores once per genome and sharing them across samples
(via `--minisplice-scores`) eliminates the predict overhead from per-sample runs.

---

## Failure modes quick-reference

| Symptom | Cause | Fix |
|---------|-------|-----|
| `minisplice not found` | minisplice not installed | `conda install -c bioconda minisplice` |
| `minisplice_mm2 requires --minisplice-model` | neither `--minisplice-model` nor `--minisplice-scores` given | provide one of the two |
| `minimap2 failed` with `--spsc` not recognized | minimap2 version doesn't support `--spsc` | use the minisplice-patched minimap2 or upgrade minimap2 |
| Low GT-AG (same as plain minimap2) | model mismatch (e.g., PCR-cDNA model used for DRS) | use a DRS-appropriate model (e.g., `vi2-7k.kan`) |
