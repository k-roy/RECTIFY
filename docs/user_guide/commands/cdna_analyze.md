# rectify cdna-analyze

Stage 3 of the ONT PCR-cDNA pipeline — post-alignment walkback, gene
assignment, isoform clustering, and Type-1 / Type-2 same-molecule pairing
on the multi-aligner consensus BAM produced by `rectify align`.

For pipeline context, see the [ONT PCR-cDNA pipeline overview](correct_cdna_overview.md).

---

## Usage

```bash
rectify cdna-analyze <consensus.bam> \
    --reference <genome.fa> \
    --gff <genes.gff> \
    -o <outdir> \
    [options]
```

`cdna-analyze` consumes the post-aligned consensus BAM from
[`rectify align`](align.md) on the Stage 1 FASTQ. Each BAM record is a
single UMI-deduplicated molecule whose per-cluster tags
(`XU` / `XC` / `XR` / `XO` / `XT` / `XY` / `XF` / `XB`, etc.) rode along
from the FASTQ comment via `minimap2 -y`. Stage 3 then:

1. Recomputes the corrected 3' end (poly-A tail length → `XA`) by walking
   back across the post-align CIGAR.
2. Recomputes the corrected 5' end (TSS via bridge-G walk-forward).
3. Assigns each cluster to a gene (`XG`) and classifies sense / antisense
   (`XS`).
4. Groups clusters into isoforms — Type-1 uses both 5' and 3' positions;
   Type-2 uses 3' only because the 5' end is random truncation noise.
5. Links Type-1 and Type-2 clusters from the same molecule at the same
   gene + orient with `|d5'| <= 5` and `|d3'| <= 5`, writing the partner
   cluster id to each record's `XL` tag.

Why post-align rather than pre-align coordinates: empirically, the
chimeric-consensus multi-aligner output gives 5' and 3' positions
accurate to within a few bp, which is required for the `tol-5 = 5`
and `tol-3 = 5` grouping defaults; pre-align positions are too coarse.

## Examples

```bash
# Standard run (defaults: tol-5 = tol-3 = 5)
rectify cdna-analyze out/stage1.rectified.bam \
    --reference genome.fa \
    --gff genes.gff \
    -o out/

# Loosen isoform-grouping tolerance (more lumping)
rectify cdna-analyze out/stage1.rectified.bam \
    --reference genome.fa \
    --gff genes.gff \
    --isoform-tol-5 10 --isoform-tol-3 10 \
    -o out/

# Tighter sense/antisense XS classifier (require larger reciprocal overlap)
rectify cdna-analyze out/stage1.rectified.bam \
    --reference genome.fa \
    --gff genes.gff \
    --min-gene-frac 0.5 --min-read-frac 0.9 \
    -o out/
```

---

## Arguments

### Required

| Argument | Description |
|----------|-------------|
| `bam` | Consensus BAM from [`rectify align`](align.md) on the `correct-cdna` FASTQ |
| `-o, --out` | Output directory (will contain `clusters.tsv`, `corrected_reads.tsv`, `isoforms.tsv`, `t1t2_pairs.tsv`, `consensus_tagged.bam`) |
| `--reference` | Genome FASTA — required for walkback (XA tail-length recomputation) |
| `--gff` | Genome annotation GFF3 — required for gene tree and `XS` classification |

### Sense / antisense classifier

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-gene-frac` | 0.3 | `gene_overlap / gene_length` threshold for the `XS` classifier |
| `--min-read-frac` | 0.8 | `gene_overlap / read_aln_length` threshold for the `XS` classifier |

### Isoform clustering

| Argument | Default | Description |
|----------|---------|-------------|
| `--isoform-tol-5` | 5 | bp tolerance on the 5' axis (Type-1 only; Type-2 ignores 5') |
| `--isoform-tol-3` | 5 | bp tolerance on the 3' axis |

### Type-1 / Type-2 reconciliation

| Argument | Default | Description |
|----------|---------|-------------|
| `--t1t2-tol-5` | 5 | bp tolerance for matching `pos5` between Type-1 and Type-2 clusters |
| `--t1t2-tol-3` | 5 | bp tolerance for matching `pos3` between Type-1 and Type-2 clusters |

---

## Output files

| File | Description |
|------|-------------|
| `<out>/clusters.tsv` | Per-cluster manifest (one row per UMI cluster) |
| `<out>/corrected_reads.tsv` | Per-molecule (UMI-consensus) corrected 3' ends, with DRS-identical columns so a single loader works across modalities |
| `<out>/isoforms.tsv` | Isoform-level aggregation (one row per isoform after `tol-5`/`tol-3` grouping) |
| `<out>/t1t2_pairs.tsv` | Type-1 / Type-2 reconciliation pairs (one row per same-molecule link) |
| `<out>/consensus_tagged.bam` | Input BAM rewritten with the new `XA` / `XG` / `XS` / `XI` / `XL` tags added per record. Indexed. |

### Schema — `clusters.tsv`

| Column | Type | Meaning |
|--------|------|---------|
| `cluster_id` | int | Internal cluster index (0-based) |
| `chrom` | string | Reference chromosome |
| `anchor` | int | Walkback-canonicalised cleavage anchor (post-align coordinate) |
| `orient` | string | `fwd` or `rev` (from Stage 1 SSP orientation) |
| `n_reads` | int | Number of input reads merged into this consensus (= `XC`) |
| `umi_canonical` | string | Canonical UMI for the cluster (= `XU`) |
| `xs` | string | Sense / antisense classification: `sense`, `antisense`, or `unannotated` |
| `xf` | int | Full-length tier from Stage 1: `0` / `1` / `2` (see [`correct-cdna`](correct_cdna.md)) |
| `xt` | int | Read type: `1` (SSP+UMI captured) or `2` (SSP-less) |
| `read_subtype` | string | `1a` / `1b` / `2` (= `XY`) |
| `tail_len` | int | Corrected poly-A tail length (post-align walkback; emitted as `XA` on the tagged BAM) |
| `gene` | string | Assigned gene name (= `XG`); empty for unannotated clusters |
| `isoform_id` | string | `<gene>_t<1\|2>_5g<i>_3g<j>` for Type-1, `<gene>_t2_3g<j>` for Type-2 (= `XI`); empty when no gene was assigned |
| `aln_start` | int | Alignment start (post-align reference_start) |
| `aln_end` | int | Alignment end (post-align reference_end) |
| `read_ids` | string | Comma-separated list of contributing input read IDs |

### Schema — `isoforms.tsv`

| Column | Type | Meaning |
|--------|------|---------|
| `isoform_id` | string | Same key as the `isoform_id` column in `clusters.tsv` |
| `gene` | string | Gene name (same for all clusters in the isoform) |
| `chrom` | string | Reference chromosome |
| `orient` | string | `fwd` or `rev` |
| `read_type` | int | `1` (Type-1 isoform, 5'+3' grouped) or `2` (Type-2 isoform, 3'-only) |
| `read_subtype` | string | `1a` / `1b` / `2` — from the first cluster in the isoform |
| `n_clusters` | int | Number of Stage-1 clusters folded into this isoform |
| `n_reads_total` | int | Sum of `n_reads` across all clusters in the isoform |
| `pos5_modal` | int | Modal 5' position across the isoform's clusters; `-1` for Type-2 isoforms (no 5' grouping) |
| `pos3_modal` | int | Modal 3' position across the isoform's clusters |
| `tail_len_median` | int | Median `tail_len` across the isoform's clusters |
| `cluster_ids` | string | Comma-separated `cluster_id` list of the constituent clusters |

### Schema — `t1t2_pairs.tsv`

One row per Type-1 / Type-2 pairing within the same `(gene, orient)` group
where both 5' and 3' termini fall within `--t1t2-tol-5` / `--t1t2-tol-3`.

| Column | Type | Meaning |
|--------|------|---------|
| `t1_cid` | int | `cluster_id` of the Type-1 partner |
| `t2_cid` | int | `cluster_id` of the Type-2 partner |
| `gene` | string | Shared gene assignment |
| `orient` | string | Shared orientation (`fwd` or `rev`) |
| `t1_pos5` | int | Type-1 corrected 5' position |
| `t1_pos3` | int | Type-1 corrected 3' position |
| `t2_pos5` | int | Type-2 5' position (raw aln-side; noisy) |
| `t2_pos3` | int | Type-2 corrected 3' position |
| `d5` | int | `\|t1_pos5 - t2_pos5\|` (always <= `--t1t2-tol-5`) |
| `d3` | int | `\|t1_pos3 - t2_pos3\|` (always <= `--t1t2-tol-3`) |
| `t1_n_reads` | int | Number of reads in the Type-1 cluster |
| `t2_n_reads` | int | Number of reads in the Type-2 cluster |
| `t1_umi` | string | Canonical UMI of the Type-1 cluster |

Both partner records also receive an `XL:i:<partner_cid>` tag in
`consensus_tagged.bam`.

### Schema — `consensus_tagged.bam`

The input BAM rewritten record-by-record. Records that survived the
walkback + classification gates receive these added tags (others pass
through unmodified):

| Tag | Type | Source |
|-----|------|--------|
| `XA` | int | Post-align corrected poly-A tail length |
| `XS` | string | `sense` / `antisense` / `unannotated` |
| `XG` | string | Gene name (only set if assigned) |
| `XI` | string | `isoform_id` (only set if assigned) |
| `XL` | int | Partner `cluster_id` from the T1 / T2 pairing (only set if matched) |

---

## Notes

- The Stage 1 `XA` (pre-align tail length) is replaced by the post-align
  walkback estimate. The Stage 1 value is more conservative because it
  cannot see post-correction CIGAR edits.
- Type-1 isoforms always have `pos5_modal != -1`; Type-2 isoforms always
  have `pos5_modal = -1` (no 5' grouping is performed).
- `XL` indicates same-molecule pairing — if a Type-2 cluster has `XL`
  set, its CPA call is corroborated by a Type-1 cluster on the opposite
  strand of the same dsDNA molecule.
- The classifier honours strict-overlap thresholds. To tighten what
  counts as `sense` / `antisense`, raise `--min-gene-frac` and
  `--min-read-frac`. Reads that fail both checks are tagged
  `unannotated`.

## See also

- [ONT PCR-cDNA pipeline overview](correct_cdna_overview.md)
- [`rectify correct-cdna`](correct_cdna.md) — Stage 1 (UMI extraction + consensus)
- [`rectify align`](align.md) — Stage 2 (multi-aligner consensus)
- [cDNA correct design spec](../../algorithms/cdna_correct.md)
