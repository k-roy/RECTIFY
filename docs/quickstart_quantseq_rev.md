# Quick Start — QuantSeq 3' REV

End-to-end walkthrough for processing QuantSeq 3' REV (dT-primed antisense
short-read) data with RECTIFY. For DRS data, see the [DRS quick start](quickstart.md);
for ONT PCR-cDNA, see the [ONT PCR-cDNA pipeline overview](user_guide/commands/correct_cdna_overview.md).

---

## What is QuantSeq REV

Lexogen QuantSeq 3' mRNA-Seq REV is a short-read 3'-end sequencing protocol
that uses **oligo-dT priming** at the poly-A tail of mature mRNAs. The library
captures the 3' end of every transcript in a strand-specific, 3'-end-focused
read population.

Key chemistry consequences (versus DRS / ONT cDNA):

- Each read maps to **one short region** near a cleavage-and-polyadenylation
  site (CPA); read body coverage is not meaningful.
- The poly-A tail is **NOT** in the read (oligo-dT primes against it and
  the tail is consumed during second-strand synthesis).
- Reads are **antisense** to the mRNA — see the strand-convention note
  below.

## Why the strand convention is inverted

QuantSeq REV is **antisense**: the sequencing read maps to the strand opposite
the mRNA. In BAM terms:

- `is_reverse=False` (BAM "+" strand) reads come from `-` strand genes
- `is_reverse=True`  (BAM "-" strand) reads come from `+` strand genes
- The 3' end of the transcript is the **LEFT** side of the BAM alignment
  (not the right, as in DRS / ONT cDNA).

RECTIFY's QuantSeq REV protocol wrapper handles this inversion internally —
the walkback routine uses `three_prime_side = "left"` and maps the BAM
`is_reverse` flag to the **opposite** gene strand. From the user's point of
view, the `corrected_3ends.tsv` `strand` column is the gene / RNA strand, not
the BAM strand.

This is the source of the most common QuantSeq REV pitfall: if you split a
QuantSeq BAM for IGV using `samtools view -f 16` / `-F 16` you will see the
two gene strands **inverted** relative to what you expect — see the
gotchas section below.

---

## Input requirements

| Input | Notes |
|-------|-------|
| Aligned BAM | Pre-aligned with BWA-MEM, BBMap, or another short-read aligner that handles spliced reads. RECTIFY does not currently align QuantSeq REV reads itself (see [target state](#target-state-multi-aligner-consensus) below). |
| Reference genome | FASTA (gzip OK). Required for read-vs-reference walkback. |
| Annotation | GFF or GTF (gzip OK). Required for gene assignment in `rectify analyze`. |

## Step-by-step

### 1. Pre-align with BWA-MEM or BBMap

QuantSeq REV reads are short (typically 75–150 bp) and 3'-end focused.
RECTIFY does not bundle a short-read aligner — pre-align with your tool of
choice and write the result to a sorted, indexed BAM.

```bash
# Example: BWA-MEM
bwa index genome.fa
bwa mem -t 8 genome.fa reads.fastq.gz \
    | samtools sort -@ 8 -o aligned.bam -
samtools index aligned.bam
```

```bash
# Example: BBMap (better splice sensitivity for yeast)
bbmap.sh ref=genome.fa in=reads.fastq.gz \
    out=aligned.sam intronlen=10 maxindel=200000 \
    threads=8
samtools sort -@ 8 -o aligned.bam aligned.sam
samtools index aligned.bam
```

### 2. Run `rectify correct` with the QuantSeq REV protocol wrapper

```bash
rectify correct aligned.bam \
    --genome genome.fa \
    --annotation genes.gff \
    --dT-primed-cDNA \
    --short-read \
    -o corrected.tsv
```

What each flag does:

- `--dT-primed-cDNA` activates the QuantSeq / dT-primed-cDNA mode. The
  poly-A tail is **not** in the read; RECTIFY engages indel-artifact
  correction at internal genomic A-tracts where the aligner over-walked
  into the priming site.
- `--short-read` disables long-read-only modules (poly-A tail trimming,
  full-read indel correction) and uses short-read tolerances. Pairs with
  `rectify align --short-read` (bbmap + bwa) once that path is wired
  into the CLI.

Internally this dispatches to the `walkback_quantseq_rev()` protocol
wrapper (`rectify/core/correct/protocols/quantseq_rev.py`):

| Parameter | Value |
|-----------|-------|
| `three_prime_side` | `"left"` (not `"right"` as in DRS) |
| `stop_base` | `"A"` (walk inward from the left until a non-A match) |
| Strand handling | BAM `is_reverse` mapped to opposite gene strand |

The output `corrected.tsv` `strand` column is the **gene / RNA strand**,
not the BAM `is_reverse` strand.

### 3. Downstream analysis with `rectify analyze`

```bash
rectify analyze corrected.tsv \
    --annotation genes.gff \
    --output-dir results/
```

`rectify analyze` is protocol-agnostic on the corrected TSV — it consumes
gene-strand positions and runs clustering, DESeq2, GO enrichment, and motif
discovery in the same way as DRS. See [`rectify analyze`](user_guide/commands/analyze.md).

---

## Strand convention gotcha — IGV BAM-split tip

Because QuantSeq REV reads are antisense, the canonical "filter by
`is_reverse`" trick splits the BAM by **DNA** strand — which is the
opposite of what you usually want when looking at gene-strand expression
in IGV. To get gene-strand views in IGV:

```bash
# Gene + strand (RNA +) = DNA - reads in the BAM
samtools view -h -f 16 aligned.bam | samtools sort -o sample.rna_plus.bam -
samtools index sample.rna_plus.bam

# Gene - strand (RNA -) = DNA + reads in the BAM
samtools view -h -F 16 aligned.bam | samtools sort -o sample.rna_minus.bam -
samtools index sample.rna_minus.bam
```

Loading these split BAMs into IGV gives proper strand-coloured tracks
that match the gene orientation. Without this split, convergent loci
(common in yeast Ysh1AA-depletion data, for example) will look
strand-confused because both signals overlay.

> **Pileup, not coverage.** For 3'-end sequencing, the meaningful
> per-position signal is the count of reads whose corrected 3' end
> maps to a base — **not** the read-body coverage. Reach for
> [`rectify export`](user_guide/commands/export.md) (which emits proper
> 3'-end pileup bedgraphs / bigWigs) rather than `samtools depth` or
> `bedtools genomecov`, which both produce smeared read-body coverage
> for QuantSeq / ONT 3'-end data.

---

## Target state — multi-aligner consensus

The QuantSeq REV path is currently single-aligner: the user pre-aligns
with one short-read aligner (BWA-MEM or BBMap) and feeds the BAM in.
**Target state** is a BWA + BBMap consensus through the existing
`chimeric_consensus.select_best()` infrastructure already used for the
DRS / ONT cDNA tracks — both aligners produce SAM/BAM that the existing
selector can consume, and short-read splice calls are notably more
robust when two aligners agree. Until the multi-aligner shim lands in
`rectify/core/correct_command.py` (the `--short-read --dT-primed-cDNA`
branch), pick whichever aligner gives the best splice sensitivity for
your organism.

---

## See also

- [`rectify correct`](user_guide/commands/correct.md) — full argument reference
- [`rectify analyze`](user_guide/commands/analyze.md) — clustering + DESeq2 + motifs
- [`rectify export`](user_guide/commands/export.md) — 3'-end pileup bedgraphs / bigWigs
- [Coordinate System](coordinate_system.md) — strand and 3'-end conventions across tracks
- DRS [Quick Start](quickstart.md)
- ONT PCR-cDNA [pipeline overview](user_guide/commands/correct_cdna_overview.md)
