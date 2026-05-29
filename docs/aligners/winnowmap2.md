# Winnowmap2: weighted minimizers for repetitive regions

Winnowmap2 suppresses false alignments in repetitive regions by down-weighting
high-frequency minimizers (identified via a `meryl` k-mer frequency database).
On ONT DRS data it achieves the highest GT-AG rate of any aligner tested, at the
same speed as minimap2.

**Benchmark result (human chr5, 10k reads, 2026-05-25):** 99.0% GT-AG, 0.4%
non-canonical, 9 s — compared to minimap2's 91.6% GT-AG at 9 s.
Weighted minimizers are especially useful at the SMN1/SMN2 locus on chr5.

See [ALIGNER_RECOMMENDATIONS.md](../ALIGNER_RECOMMENDATIONS.md) for the full
benchmark table.

---

## Installation

```bash
conda install -c bioconda winnowmap meryl
```

Winnowmap installs as the binary `winnowmap` (not `winnowmap2`). RECTIFY's
`run_winnowmap2()` checks `winnowmap` first, then falls back to `winnowmap2`.
`meryl` is required to build the repetitive k-mer list.

---

## How RECTIFY runs it

`run_winnowmap2()` in `multi_aligner.py`:

1. **Builds (or reuses) a repetitive k-mers file** via meryl:
   ```bash
   meryl count k=15 output <merylDB> <genome>
   meryl print greater-than distinct=0.9998 <merylDB> > <repetitive_k15.txt>
   ```
   The meryl database and repetitive k-mers file are cached in
   `.winnowmap_cache/` adjacent to the genome (or in `cache_dir`). Subsequent
   runs for the same genome skip this step.

2. **Aligns** with the standard splice-aware minimap2 flags plus `-W`:
   ```bash
   winnowmap -W repetitive_k15.txt -ax splice -uf -k14 \
       -G <max_intron> --secondary=no --MD -y -t <threads> \
       <genome> <reads> | samtools sort -n ...
   ```

3. Applies `samtools calmd -e`, validates QNAMEs, and injects `RN` tags.

---

## CLI usage

```bash
rectify align <reads.fastq.gz> --genome <genome.fa> \
    --aligners winnowmap2 \
    --max-intron 500000 \
    -o /tmp/wm2_out -t 8
```

Pass a pre-computed repetitive k-mers file to skip the meryl build step:

```bash
rectify align <reads.fastq.gz> --genome <genome.fa> \
    --aligners winnowmap2 \
    --winnowmap-repetitive-kmers /path/to/repetitive_k15.txt \
    -o /tmp/wm2_out -t 8
```

`winnowmap2` is not in the default `"all"` aligner set — it must be requested
explicitly. To combine with the default panel:

```bash
rectify align <reads.fastq.gz> --genome <genome.fa> \
    --aligners minimap2 mapPacBio gapmm2 winnowmap2 \
    --junction-aligners uLTRA deSALT \
    -o /tmp/out -t 8
```

---

## Performance notes

- Meryl build time: ~1–5 min for yeast; ~15–30 min for GRCh38 (cached after first run).
- Alignment: same speed as minimap2 (~9 s for 10k reads on chr5).
- Memory: slightly higher than minimap2 due to k-mer weight table.

---

## Failure modes quick-reference

| Symptom | Cause | Fix |
|---------|-------|-----|
| `meryl not found` | meryl not installed | `conda install -c bioconda meryl` |
| `winnowmap not found` | winnowmap not installed | `conda install -c bioconda winnowmap` |
| `meryl count failed` | gzipped genome (some meryl versions) | Decompress genome first: `zcat genome.fa.gz > genome.fa` |
| Empty BAM | Genome path mismatch | Confirm genome path is uncompressed FASTA |

---

## Primary-alignment & duplicate handling

`rectify align` passes **`--secondary=no`** to winnowmap (`multi_aligner.py:1531`;
see the `-W ... --secondary=no --MD -y` line under "How RECTIFY runs it"), so
rectify-produced winnowmap2 BAMs carry no secondary records. Supplementary records
from chimeric/split reads can still occur (winnowmap, like minimap2, has no flag to
suppress them) but are dropped by `rectify correct`'s `is_supplementary` filter.

The unguarded hazard for *all* aligners is duplicate **primary** records in an
external BAM (e.g. a doubled input FASTQ), which `rectify correct` does **not**
dedup → 2× double-counted 3′ ends. See the canonical writeup and cross-aligner
table in [minimap2.md](minimap2.md#-duplicate-primary-alignments--2-double-counted-3-ends-external-bam-hazard).
