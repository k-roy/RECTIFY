# rectify correct-cdna

UMI-aware Stage 1 of the ONT PCR-cDNA pipeline (PCB114.24 chemistry).

`rectify correct-cdna` consumes a pre-aligned PCR-cDNA BAM and collapses PCR
siblings of the same starting RNA molecule into a single consensus record.
Per-cluster tags are written on TAB-separated FASTQ comments so that
`rectify align -y` propagates them to the post-alignment BAM for downstream
analysis by [`rectify cdna-analyze`](cdna_analyze.md).

For pipeline context, see the [ONT PCR-cDNA pipeline overview](correct_cdna_overview.md).

---

## Usage

```bash
rectify correct-cdna <input.bam> --reference <genome.fa> -o <outdir> [options]
```

## When to use

Use `rectify correct-cdna` when your input is an Oxford Nanopore PCR-cDNA BAM
built with the SQK-PCB114.24 kit (or a chemistry that places a 27-nt structured
UMI between the strand-switching primer (SSP) and a GGG template-switching
bridge). The command:

- Extracts the 27-nt UMI from each read by regex-matching the SSP motif
  (`TTTCTGTTGGTGCTGATATTGCT`) on either orientation.
- Clusters reads by `(chrom, anchor_window, orient, UMI)` using **directional
  clustering** with a 2x count rule (default; see `--umi-clustering`).
- Builds an abPOA consensus per multi-read cluster (or passes a singleton
  through unchanged).
- Strips SSP / UMI / GGG bridge at the 5' end and any basecalled poly-A tail
  at the 3' end so the downstream aligner receives a clean mRNA body.

Do **not** use this command for direct RNA-seq (DRS) or QuantSeq REV input —
those have no UMI and route through `rectify correct` instead.

## Examples

```bash
# Single sample, full genome
rectify correct-cdna pcb114.bam \
    --reference genome.fa \
    --gff genes.gff \
    -o out/

# Single-chromosome test run
rectify correct-cdna pcb114.bam \
    --reference genome.fa \
    --gff genes.gff \
    --region chrI \
    -o out/

# v1.18 strand-aware POA (cancels strand-specific basecaller bias)
rectify correct-cdna pcb114.bam \
    --reference genome.fa \
    --gff genes.gff \
    --strand-aware-consensus \
    -o out/
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `bam` | Input BAM — aligned PCR-cDNA reads, indexed (`.bai` alongside) |
| `-o, --out` | Output directory (will contain `stage1_consensus.fastq.gz`) |

### Reference data

| Argument | Description |
|----------|-------------|
| `--reference` | Genome FASTA (gzip OK). Required for walk-back anchor canonicalisation and pre-align poly-A tail-length measurement during UMI extraction. Without it, the legacy aln-end anchor is used. |
| `--gff` | Genome annotation GFF3 (gzip OK). Required for sense / antisense `XS` classification, `XG` gene-name tagging, isoform clustering, and rDNA masking. |

### UMI clustering

| Argument | Default | Description |
|----------|---------|-------------|
| `--umi-edit-distance` | 3 | Max Levenshtein distance between UMIs in the same cluster |
| `--anchor-window-bp` | 25 | Window around alignment-end anchor for same-locus clustering (bp) |
| `--per-cluster-cap` | 200 | Max reads per cluster — guards against polyA-pileup hot-spots (chrXII rDNA, etc.) |
| `--umi-clustering` | directional | UMI clustering method: `directional` (umi_tools 2x rule) or `components` (connected components) |

### Consensus

| Argument | Default | Description |
|----------|---------|-------------|
| `--no-poa` | off | Force pileup-only consensus even if abPOA is available |
| `--strand-aware-consensus` | off | v1.18: split reads by `is_reverse`, build per-strand sub-consensuses, then merge — cancels strand-specific systematic basecaller errors |

### Filtering

| Argument | Description |
|----------|-------------|
| `--region CHROM` | Restrict to one BAM region (e.g. `chrI`) — useful for testing |
| `--no-mask-rdna` | Disable rDNA masking. By default, reads overlapping `rRNA_gene` loci in `--gff` are excluded to prevent the O(n^2) UMI bottleneck on chrXII tandem repeats. |

### Logging

| Argument | Description |
|----------|-------------|
| `-v, --verbose` | Enable DEBUG-level logging |

---

## Output files

| File | Description |
|------|-------------|
| `<out>/stage1_consensus.fastq.gz` | One FASTQ record per UMI cluster. Per-cluster SAM-format tags are appended as TAB-separated comment fields so that `minimap2 -y` propagates each `XX:T:value` into its own BAM aux entry. |

The FASTQ is the single canonical output. Feed it directly to
[`rectify align`](align.md) (Stage 2), then to
[`rectify cdna-analyze`](cdna_analyze.md) (Stage 3) on the post-aligned BAM.

### Per-cluster tag glossary

The following tags are emitted on every FASTQ record (and survive through to
the post-alignment BAM via `minimap2 -y`):

| Tag | Type | Meaning |
|-----|------|---------|
| `XU:Z` | string | Canonical UMI (representative or pileup-consensus UMI) |
| `XO:Z` | string | Orientation: `fwd` (SSP found on forward strand of BAM SEQ) or `rev` (SSP found as reverse complement) |
| `XC:i` | int | Cluster size — number of raw reads merged into this consensus |
| `XR:Z` | string | Comma-separated list of input read IDs that contributed |
| `XM:Z` | string | Consensus method: `poa`, `strand_split_poa`, or `pileup` |
| `XF:i` | int | Full-length tier — `0` = not detected, `1` = unanchored polyA/T match, `2` = anchored polyA/T at adapter (HIGH confidence) |
| `XA:i` | int | Pre-align poly-A tail length (median over cluster); recomputed by `cdna-analyze` on post-align CIGAR |
| `XT:i` | int | Read type — `1` = SSP+UMI captured, `2` = SSP-less (5'-truncated; not UMI-deduplicated) |
| `XY:Z` | string | Read subtype — `1a` (Type-1, orient=fwd), `1b` (Type-1, orient=rev), `2` (Type-2, no UMI) |
| `XQ:i` | int | 5' pre-trim bases stripped (SSP+UMI+GGG for Type-1, polyT for orient=rev) |
| `XK:i` | int | 3' pre-trim bases stripped (polyA for orient=fwd, SSP_RC suffix for orient=rev) |
| `XB:Z` | string | Strand-split count `n_top/n_bottom` (only meaningful with `--strand-aware-consensus`) |

> **Tag namespace.** `X[upper]` tags are persistent user-visible metadata
> owned by the cDNA pipeline. `rectify align`'s internal aligner-selection
> bookkeeping uses `X[lower]` (`Xa`/`Xc`/`Xn`/`Xj`/...) — those are debug
> metadata and are not stable for downstream consumers.

### Read-type classification

| `XY:Z` | `XT:i` | When | Deduplication |
|--------|--------|------|---------------|
| `1a` (`umi_captured_fwd`) | 1 | SSP + UMI captured at basecalled 5' (full-length sense molecule) | UMI-anchored |
| `1b` (`umi_captured_rev`) | 1 | SSP_RC + UMI_RC captured at basecalled 3' (pA-first read traveled far enough) | UMI-anchored |
| `2` (`umi_not_captured`) | 2 | pA-first read truncated before reaching SSP/UMI | Not merged — each read is one observation |

---

## Notes

- The input BAM must be **pre-aligned** — typically with `minimap2 -ax splice`
  using a yeast or other reference. `rectify correct-cdna` consumes its
  positions as anchors for UMI bucketing; it does not align reads itself.
- The `XA` (pre-align tail length) tag emitted here is a coarse estimate used
  only for the Stage 1 full-length tier (`XF`); `rectify cdna-analyze`
  recomputes `XA` on the post-align CIGAR, which is the value that should be
  consumed downstream.
- The first run on a fresh BAM with many polyA-pileup hot-spots can be slow
  if rDNA masking is disabled — keep `--no-mask-rdna` off unless you have
  manually filtered chrXII.

---

## Next step

```bash
rectify align out/stage1_consensus.fastq.gz \
    --genome genome.fa \
    --prefix stage1 \
    -o out/
```

Then proceed to [`rectify cdna-analyze`](cdna_analyze.md).

## See also

- [ONT PCR-cDNA pipeline overview](correct_cdna_overview.md)
- [`rectify align`](align.md) — Stage 2 multi-aligner consensus
- [`rectify cdna-analyze`](cdna_analyze.md) — Stage 3 post-align isoform clustering
- [`rectify correct`](correct.md) — for DRS and QuantSeq REV inputs
- [cDNA correct design spec](../../algorithms/cdna_correct.md)
