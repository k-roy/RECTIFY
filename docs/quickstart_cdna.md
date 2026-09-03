# Quick Start — ONT PCR-cDNA (SQK-PCB114.24)

End-to-end walkthrough for processing Oxford Nanopore PCR-cDNA data from the
SQK-PCB114.24 kit with RECTIFY. For DRS data, see the
[DRS quick start](quickstart.md); for QuantSeq REV, see the
[QuantSeq REV quick start](quickstart_quantseq_rev.md).

---

## What is SQK-PCB114.24

PCB114.24 is Oxford Nanopore's current PCR-based cDNA kit (native cDNA via
reverse transcription, amplified by PCR, long-read sequenced). Each starting
RNA molecule is reverse-transcribed and PCR-amplified, so a single nanopore
run contains **multiple PCR copies of every starting molecule**. Unless you
deduplicate by unique molecular identifier (UMI), every copy is counted
separately, inflating abundance estimates and CPA-site precision.

### The UMI architecture

PCB114.24 places a **27-nt structured UMI** immediately downstream of the
strand-switching primer (SSP) motif:

```
5'— [SSP motif: TTTCTGTTGGTGCTGATATTGCT] — [27-nt UMI: (TT-VVVV)×4 + TTT] — GGG — [mRNA body] — [polyT/VNP] —3'
```

The 16 informative V positions supply ~43 M unique barcodes. Key empirical
properties (validated on *S. cerevisiae* `wt_rep1` chrI, 40K reads):

- **PCR siblings are at Hamming=0.** Basecaller accuracy on the structured UMI
  is remarkably high; there is essentially no d=1–3 substitution smear. Indels
  in the T-runs do occur, so **Levenshtein distance, not Hamming, is mandatory**
  for correct clustering.
- **Sense and antisense reads of the same molecule share the same UMI.**
  Dorado does *not* reverse-complement antisense reads, so both orientations of
  a dsDNA PCR copy are present in the BAM. `rectify correct-cdna` handles both.
- **~81.7 % of reads carry an extractable UMI.** The remaining ~18 % had
  basecaller errors in the SSP region or were truncated before reaching it;
  they are classified as **Type-2 reads** and handled separately (see below).

---

## Three-stage pipeline overview

| Stage | Command | What it does |
|-------|---------|--------------|
| **1** | `rectify correct-cdna` | Extract UMIs, directional-cluster PCR siblings, build abPOA consensus per cluster, strip adapters. Emits one FASTQ record per UMI cluster with per-cluster tags. |
| **2** | `rectify align` | Multi-aligner consensus on the Stage 1 FASTQ; propagates Stage 1 tags into the BAM via `minimap2 -y`. |
| **3** | `rectify cdna-analyze` | Post-align walkback + walk-forward, gene assignment, isoform clustering, Type-1 ↔ Type-2 same-molecule pairing. |

The critical design choice: **UMI deduplication runs before alignment (Stage 1),
not after.** This means the aligner in Stage 2 sees one canonical read per
starting molecule with adapter sequences already stripped, giving it the best
possible input for splice-junction calling.

---

## Installation

UMI clustering (`rapidfuzz`) and per-cluster re-alignment (`mappy`) ship with
the core `rectify-rna` package. For the best consensus quality, install the
`[cdna-correct]` extra, which adds fuzzy adapter anchoring (`edlib`) and
abPOA consensus (`pyabpoa`):

```bash
pip install "rectify-rna[cdna-correct]"
```

| Package | Role | Tier |
|---------|------|------|
| `rapidfuzz` | Fast Levenshtein distance for UMI clustering | core |
| `mappy` | In-process minimap2 for per-cluster consensus re-alignment | core |
| `edlib` | Fuzzy (Lev≤2) adapter-anchor alignment | `[cdna-correct]` |
| `pyabpoa` | Partial-order-alignment (POA) consensus across PCR siblings | `[cdna-correct]` |

Without the `[cdna-correct]` extras, `correct-cdna` falls back to exact-match
adapter anchors and pileup-style consensus — still functional, just less robust.

---

## Input requirements

| Input | Notes |
|-------|-------|
| Pre-aligned BAM | Aligned with any long-read splice-aware aligner (minimap2 `-ax splice` recommended). Must be sorted and indexed (`.bai` alongside). |
| Reference genome | FASTA (gzip OK). Required for anchor canonicalisation and poly-A tail estimation during UMI extraction. |
| Annotation | GFF3 (gzip OK). Required for sense/antisense classification, gene assignment (`XG`), and rDNA masking. |

> **Pre-alignment note.** `rectify correct-cdna` uses the pre-aligned positions
> only as anchors for grouping PCR siblings into locus buckets. It does not
> trust the alignment for end-position calling — that happens post-consensus in
> Stage 3. Any splice-aware long-read aligner works; minimap2 is the default.

---

## Step by step

The cDNA pipeline runs as three explicit stages. Because Stage 1 deduplicates
by UMI *before* alignment, there is no single `run-all` shortcut for cDNA — run
the three commands in sequence.

### Stage 1 — UMI extraction + abPOA consensus

```bash
rectify correct-cdna pcb114.bam \
    --reference genome.fa \
    --gff genes.gff \
    -o out/
```

**What happens internally:**

1. **UMI extraction.** For each read, search for the SSP motif
   (`TTTCTGTTGGTGCTGATATTGCT`) on either strand. If found, extract the 27 bases
   immediately downstream as the UMI. Reads where the SSP is found at the
   basecalled 5' end are tagged `orient=fwd`; reads where the SSP_RC is found
   at the basecalled 3' end are tagged `orient=rev`. Reads with no extractable
   UMI are classified as **Type-2** and handled separately.
2. **Anchor bucketing.** Group reads by `(chrom, anchor_window ±25 bp, orient)`.
   The anchor is the basecalled-3' genomic position (the end opposite the
   SSP/UMI side).
3. **Directional UMI clustering.** Within each anchor bucket, cluster UMIs by
   Levenshtein ≤ 3 using the `umi_tools` directional 2x rule. Each connected
   component is one starting-molecule cluster.
4. **abPOA consensus.** For multi-read clusters, build a partial-order-alignment
   consensus across the PCR siblings' cDNA inserts (SSP/UMI/GGG already
   stripped). Singleton clusters pass through unchanged.
5. **Output.** Writes `out/stage1_consensus.fastq.gz` — one record per UMI
   cluster, with per-cluster metadata on TAB-separated FASTQ comment fields so
   `minimap2 -y` can propagate them as BAM aux fields in Stage 2.

**Key flags:**

| Flag | Default | Notes |
|------|---------|-------|
| `--umi-edit-distance` | 3 | Levenshtein threshold for UMI clustering. 3 is empirically validated; ≥ 5 starts admitting random collisions |
| `--anchor-window-bp` | 25 | ±bp window on basecalled-3' for same-locus grouping |
| `--per-cluster-cap` | 200 | Drops anchor buckets with >N reads. Required: polyA-pileup hot-spots (rDNA loci on chrXII) have thousands of distinct molecules at the same anchor and are computationally expensive |
| `--region CHROM` | — | Restrict to one chromosome for smoke-testing |
| `--no-poa` | off | Skip abPOA; use pileup-style consensus instead (faster, lower quality) |
| `--no-mask-rdna` | off | Disable rDNA masking (see Gotchas — not recommended) |

### Stage 2 — Multi-aligner consensus

```bash
rectify align out/stage1_consensus.fastq.gz \
    --genome genome.fa \
    --prefix stage1 \
    -o out/
```

Runs the standard RECTIFY multi-aligner panel (minimap2 + uLTRA + deSALT,
with the overhang resolver on the minimap2 arm) on the per-cluster FASTQ
from Stage 1. The
`--prefix stage1` flag ensures the output BAM is named `out/stage1.rectified.bam`.
`minimap2 -y` propagates each `XX:T:value` tag from the FASTQ comment into its
own BAM auxiliary field — the Stage 1 metadata (`XU`, `XC`, `XO`, `XT`, `XY`,
`XF`, etc.) rides through intact to Stage 3 without any custom BAM rewriting.

> **Do not feed the original raw BAM to `rectify align`.** The raw BAM contains
> unseparated PCR siblings with adapter sequences attached. Always align the
> `stage1_consensus.fastq.gz` output from Stage 1.

### Stage 3 — Post-align walkback + isoform clustering

```bash
rectify cdna-analyze out/stage1.rectified.bam \
    --reference genome.fa \
    --gff genes.gff \
    -o out/
```

**What happens internally:**

1. **Post-align walkback.** Recomputes the corrected 3' end (poly-A tail
   length → `XA`) by walking back across the post-alignment CIGAR. The Stage 1
   `XA` estimate (pre-alignment) is replaced here — post-align is the
   authoritative value because it sees the corrected splice-junction CIGAR.
2. **Walk-forward.** Recomputes the corrected 5' end (TSS via template-switching
   bridge-G walk-forward) for Type-1 reads.
3. **Gene assignment + sense/antisense.** Assigns each cluster to a gene (`XG`)
   via interval overlap with the GFF annotation. Classifies as `sense`,
   `antisense`, or `unannotated` (`XS`).
4. **Isoform clustering.** Type-1 reads group by `(gene, orient, pos5 ±tol-5,
   pos3 ±tol-3)`; Type-2 reads group by `(gene, orient, pos3 ±tol-3)` only
   (their 5' end is stochastic truncation noise, not a real TSS). Default
   tolerance: ±5 bp on each axis (configurable via `--isoform-tol-5` and
   `--isoform-tol-3`).
5. **Type-1 / Type-2 pairing.** Links each Type-2 cluster to the matching
   Type-1 cluster (same gene, orient, and both termini within tolerance) by
   writing the partner cluster id to `XL`. This recovers 3'-end CPA calls from
   the deeper Type-2 pool, corroborated by a full-length Type-1 read on the
   opposite dsDNA strand.

**Key flags:**

| Flag | Default | Notes |
|------|---------|-------|
| `--isoform-tol-5` | 5 | bp tolerance on 5' axis for isoform grouping (Type-1 only) |
| `--isoform-tol-3` | 5 | bp tolerance on 3' axis for isoform grouping |
| `--t1t2-tol-5` | 5 | bp tolerance for T1 ↔ T2 pairing on 5' axis |
| `--t1t2-tol-3` | 5 | bp tolerance for T1 ↔ T2 pairing on 3' axis |
| `--min-gene-frac` | 0.3 | Min `gene_overlap / gene_length` for sense/antisense classification |
| `--min-read-frac` | 0.8 | Min `gene_overlap / read_length` for sense/antisense classification |

---

## Output files

After all three stages:

```
out/
├── stage1_consensus.fastq.gz   # Stage 1: one record per UMI cluster (→ input to Stage 2)
├── stage1.rectified.bam        # Stage 2: multi-aligner consensus BAM (→ input to Stage 3)
├── stage1.rectified.bam.bai
├── clusters.tsv                 # Stage 3: per-cluster manifest
├── corrected_reads.tsv          # Stage 3: per-molecule 3' ends, DRS-equivalent columns
├── isoforms.tsv                 # Stage 3: isoform-level aggregation
├── t1t2_pairs.tsv               # Stage 3: Type-1 / Type-2 same-molecule links
└── consensus_tagged.bam         # Stage 3: BAM with XA / XG / XS / XI / XL tags added
```

The `corrected_reads.tsv` emitted here uses the same per-read column schema as
the DRS `rectify correct` output, so a single loader works across modalities.

### `clusters.tsv` — one row per UMI cluster

Key columns: `cluster_id`, `chrom`, `anchor` (walkback-corrected 3' position),
`orient` (`fwd`/`rev`), `n_reads` (PCR siblings merged), `umi_canonical`,
`xs` (sense/antisense/unannotated), `xf` (full-length tier 0/1/2), `xt`
(read type 1/2), `tail_len` (corrected poly-A tail length), `gene`, `isoform_id`.

### `isoforms.tsv` — one row per isoform

Key columns: `isoform_id`, `gene`, `n_clusters` (UMI clusters), `n_reads_total`,
`pos5_modal`, `pos3_modal`, `tail_len_median`. Type-2 isoforms have
`pos5_modal = -1` (no 5' grouping).

### `t1t2_pairs.tsv` — Type-1 / Type-2 pairs

One row per same-molecule link. Columns: `t1_cid`, `t2_cid`, `gene`, `orient`,
`t1_pos5`, `t1_pos3`, `t2_pos3`, `d5`, `d3`, `t1_n_reads`, `t2_n_reads`.

---

## Read types and orientations

Understanding Type-1 vs Type-2 and `fwd` vs `rev` is essential for
interpreting the output.

### Type-1 vs Type-2

| | Type-1 (`XT=1`) | Type-2 (`XT=2`) |
|---|---|---|
| SSP found? | Yes | No |
| UMI available? | Yes | No |
| Deduplication | UMI-anchored | **None — each read is one observation** (enforced; see note) |
| 5' end in output | Corrected (walk-forward) | Stochastic truncation — unreliable |
| 3' end in output | Corrected (walkback) | Corrected (walkback) |
| Isoform grouping | By 5' + 3' | By 3' only |

> 🔴 **Type-2 reads are never deduplicated, and this is now enforced in code.** Until the
> `--type2-collapse` flag was added, Stage 1 silently collapsed Type-2 reads that shared an
> exact `(aln_start, aln_end)` and counted them as PCR duplicates. With no UMI there is no
> evidence two Type-2 reads are the same molecule, so that measured **positional
> concentration, not amplification** — on one 18-library cohort it removed **51.5 % of all
> Type-2 reads**, scaling with depth (4-6 % at ~50 k reads, 44-57 % at multi-million) while
> true UMI-measured duplication on the same libraries was 24-41 %. The default is now
> `--type2-collapse none`; `position` reproduces the old behaviour and should be used only
> for that.

Type-2 reads are NOT discarded: they carry a valid, corrected 3' end and are
clustered separately by CPA site **in Stage 3** — that is isoform/CPA clustering, which is a
different operation from deduplication and is labelled as such. If a matching Type-1 cluster exists, the pair
is recorded in `t1t2_pairs.tsv` and the Type-2 cluster gains the
`XL` same-molecule cross-reference — its CPA call is then corroborated by a
full-length read.

### `fwd` vs `rev` orientation

These label which strand of the dsDNA cDNA molecule the read came from, **not**
the gene strand:

| `XO` tag | What it means |
|----------|---------------|
| `fwd` | SSP motif found at the basecalled 5' end → polyT/VNP is at the basecalled 3' end |
| `rev` | SSP_RC found at the basecalled 3' end → the read traversed the molecule in the other direction |

Both orientations of the same molecule share the same UMI. `rectify correct-cdna`
clusters `fwd` and `rev` reads separately (each within its own anchor bucket);
same-molecule reconciliation across read types happens later, in the Stage 3
Type-1 ↔ Type-2 pairing.

> `fwd`/`rev` does **not** correlate directly with BAM `FLAG is_reverse`. A
> `fwd`-orient read (SSP at basecalled 5') may map to either BAM strand depending
> on which RNA strand the gene is on.

---

## QC — Stage 1 summary

Stage 1 (`rectify correct-cdna`) prints a per-run QC summary to the console and
writes a provenance sidecar (JSON) alongside the FASTQ output. The console
summary reports:

- Total reads processed and the consensus-method breakdown (singletons passed
  through, multi-read POA/pileup, fallback).
- **Read-type breakdown** — Type-1 (SSP + UMI captured) vs Type-2 (SSP-less,
  5'-truncated), with per-type read and cluster counts. ~80 % Type-1 is typical
  for clean PCB114 yeast data; a much lower fraction suggests a chemistry
  mismatch or heavy truncation (see Gotchas).
- **XF tier breakdown** — full-length confidence tier counts.
- polyA-pileup buckets dropped by `--per-cluster-cap`.

---

## Gotchas

### rDNA pileup on chrXII (yeast)

The 100+ tandem rDNA repeats on yeast chrXII produce enormous polyA-pileup
clusters at rRNA gene loci. These have thousands of distinct molecules at the
same genomic anchor, and naive UMI clustering inside them is computationally
expensive. By default, `rectify correct-cdna` masks reads overlapping
`rRNA_gene` loci in the GFF annotation. **Do not disable this (`--no-mask-rdna`)
unless you have manually filtered chrXII reads beforehand.**

### Cross-orientation span filter

The Stage 2 cross-orient merge pairs `fwd` and `rev` consensuses by shared UMI
**within a genomic span filter** (`--max-cross-orient-span`, default 3000 bp for
yeast). This filter is mandatory: empirically, ~95% of cross-orient UMI matches
are random collisions; only pairs within the expected mRNA span are true
same-molecule pairs. Without the filter, thousands of unrelated molecules would
be incorrectly merged.

For organisms with larger average mRNA span (e.g., *H. sapiens*), increase this
value:

```bash
rectify correct-cdna pcb114.bam --reference genome.fa --gff genes.gff \
    --max-cross-orient-span 200000 \
    -o out/
```

### `--per-cluster-cap` protects against O(N²) hot-spots

Genomic loci with very high read pileup (rDNA after masking fails to catch a
locus, or exceptionally highly expressed genes) can have hundreds of reads per
anchor bucket. The pairwise Levenshtein clustering inside these buckets is O(K²).
The default cap of 200 reads per bucket prevents runaway runtimes on edge cases.
If you see unexpectedly slow Stage 1 on a specific chromosome, check which locus
is the bottleneck with `--verbose` and either add it to the GFF rRNA mask or
reduce `--per-cluster-cap`.

### ~18 % of reads have no extractable UMI

Reads where the SSP motif is not found (basecaller errors in the structured
adapter, or reads truncated before reaching the SSP) are classified as Type-2.
They are **not dropped** — they contribute valid corrected 3' ends but are
excluded from UMI-anchored deduplication. This is expected behaviour.
Check the Type-1 fraction in the Stage 1 console summary; a Type-1 fraction
well below ~70 % suggests a chemistry mismatch (wrong SSP motif or UMI length)
or very heavy truncation.

### Do not mix Type-1 and Type-2 for isoform-level analyses

For TSS analysis (`pos5_modal`) and full-length isoform characterization, use
only Type-1 clusters (`xt == 1`). Type-2 5' positions are random truncation
endpoints, not true TSSs. The `isoforms.tsv` file already marks Type-2 isoforms
with `read_type = 2` and `pos5_modal = -1`.

### Strand-split consensus (`--strand-aware-consensus`)

By default, abPOA builds consensus across all reads in a cluster regardless of
their BAM strand. For clusters that contain both FLAG=0 and FLAG=16 reads, the
consensus may be blurred by strand-specific systematic basecaller errors. The
`--strand-aware-consensus` flag builds separate per-strand sub-consensuses
first and then merges them, cancelling this bias. It is slower
(~2× Stage 1 runtime) and recommended for deeply-sequenced samples or when
per-strand error profiles are expected to differ substantially.

---

## Multi-sample analysis

The cDNA pipeline has no single-command multi-sample mode. Run the three
stages per sample (e.g. in a shell loop), build an analyze manifest whose
`path` column points at each sample's Stage 3 `corrected_reads.tsv`, then run
`rectify analyze --manifest` for the combined differential analysis.

```bash
# Input manifest (input.tsv): columns sample_id, bam, condition
# Process each sample through all three stages.
printf 'sample_id\tpath\tcondition\n' > corrected_manifest.tsv
while IFS=$'\t' read -r sample bam condition; do
    rectify correct-cdna "$bam"                                       --reference genome.fa --gff genes.gff -o "results/$sample/"
    rectify align         "results/$sample/stage1_consensus.fastq.gz" --genome genome.fa --prefix stage1 -o "results/$sample/"
    rectify cdna-analyze  "results/$sample/stage1.rectified.bam"      --reference genome.fa --gff genes.gff -o "results/$sample/"
    # Record the per-sample corrected TSV for the combined analyze step
    printf '%s\t%s\t%s\n' "$sample" "results/$sample/corrected_reads.tsv" "$condition" >> corrected_manifest.tsv
done < <(tail -n +2 input.tsv)

# Combined analysis across conditions
rectify analyze /dev/null \
    --manifest corrected_manifest.tsv \
    --annotation genes.gff \
    --reference wt \
    --run-deseq2 \
    -o results/combined/
```

The `combined/` directory aggregates CPA-site clusters and runs DESeq2 on
CPA-site counts across conditions.

---

## HPC / large datasets

Stage 1 (`rectify correct-cdna`) is the most memory-intensive step — plan for
~4–8 GB RAM per sample at yeast scale. Stage 2 (`rectify align`) fans out to
the multi-aligner panel and benefits most from parallelism. For SLURM, run one
array task per sample with the three stages chained inside each task (see the
[HPC SLURM guide](user_guide/hpc_slurm.md) for the task-array template and
per-sample memory sizing).

---

## Tag namespace summary

The cDNA pipeline writes `X[upper]` BAM tags — persistent, user-visible metadata
that survives through all three stages:

| Tag | Stage set | Meaning |
|-----|-----------|---------|
| `XU:Z` | 1 | Canonical UMI |
| `XO:Z` | 1 | Orientation: `fwd` or `rev` |
| `XC:i` | 1 | Cluster size (PCR siblings merged) |
| `XR:Z` | 1 | Comma-separated input read IDs |
| `XM:Z` | 1 | Consensus method: `poa`, `strand_split_poa`, or `pileup` |
| `XF:i` | 1 | Full-length tier: `0` none, `1` unanchored, `2` HIGH confidence |
| `XT:i` | 1 | Read type: `1` (Type-1, UMI captured) or `2` (Type-2, no UMI) |
| `XY:Z` | 1 | Read subtype: `1a` (fwd), `1b` (rev), `2` |
| `XQ:i` | 1 | 5' pre-trim bases stripped |
| `XK:i` | 1 | 3' pre-trim bases stripped |
| `XA:i` | 1 (pre), **3 (post)** | Poly-A tail length. Stage 1 = coarse pre-align estimate; Stage 3 overwrites with post-align walkback value (use this one) |
| `XG:Z` | 3 | Assigned gene name |
| `XS:Z` | 3 | `sense` / `antisense` / `unannotated` |
| `XI:Z` | 3 | Isoform ID |
| `XL:i` | 3 | Partner cluster ID (Type-1 ↔ Type-2 same-molecule pairing) |

`X[lower]` tags (`Xa`, `Xc`, `Xn`, etc.) are internal `rectify align` aligner-
selection bookkeeping — not stable for downstream consumers.

---

## Performance

Typical runtimes on M1 (8 cores), *S. cerevisiae* single sample (~1.5 M reads):

| Stage | Typical time |
|-------|-------------|
| Stage 1 (`correct-cdna`) | ~90 s |
| Stage 2 (`rectify align`) | ~3–5 min (multi-aligner panel) |
| Stage 3 (`cdna-analyze`) | ~30 s |

---

## See also

- [`rectify correct-cdna`](user_guide/commands/correct_cdna.md) — Stage 1 full argument reference
- [`rectify align`](user_guide/commands/align.md) — Stage 2 aligner panel
- [`rectify cdna-analyze`](user_guide/commands/cdna_analyze.md) — Stage 3 full argument reference
- [ONT PCR-cDNA pipeline overview](user_guide/commands/correct_cdna_overview.md) — architecture diagram + tag glossary
- [cDNA correct design spec](algorithms/cdna_correct.md) — empirical UMI validation and algorithm details
- [DRS quick start](quickstart.md) — for direct RNA-seq data
- [QuantSeq REV quick start](quickstart_quantseq_rev.md) — for short-read 3'-end data
