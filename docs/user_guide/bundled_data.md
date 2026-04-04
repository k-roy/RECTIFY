# Bundled Data (Yeast)

For *Saccharomyces cerevisiae*, RECTIFY bundles all reference data required for a complete run. No external files are needed.

---

## Activating bundled data

Use either flag:

```bash
--Scer              # Shorthand
--organism yeast    # Long form (also accepts: saccharomyces_cerevisiae)
```

Both flags are recognized by `rectify correct`, `rectify run`, `rectify align`, `rectify analyze`, and `rectify netseq`.

---

## What is bundled

| Data | Source | Notes |
|------|--------|-------|
| **Genome** | *S. cerevisiae* S288C R64-5-1 | FASTA (gzipped); includes all 16 chromosomes + mitochondria |
| **Gene annotations** | SGD (Saccharomyces Genome Database) | GFF format; 6,000+ genes, UTRs, introns |
| **GO annotations** | SGD GO slim | 5,000+ gene–GO term associations |
| **WT NET-seq** | Published dataset | Strand-specific BigWig; pre-deconvolved |
| **A-tract CPA sites** | Pre-computed | 64K known CPA sites in A-tract regions (for model validation) |
| **Poly(A) model** | Trained on yeast WT | JSON; used by `rectify correct` for A-richness scoring |
| **Spike-in sequences** | ENO2 k-mers | Used by `--filter-spikein` |

### Bundled data location

Data lives in the `rectify/data/` directory inside the installed package:

```
rectify/data/
├── genomes/
│   └── saccharomyces_cerevisiae/
│       ├── sacCer3.fa.gz           # Genome FASTA
│       ├── SGD_genes.gff.gz        # Gene annotations
│       └── bbmap_index/            # Pre-built BBMap index for mapPacBio
├── models/
│   └── polya_tail_model.json       # Trained poly(A) model
├── go_annotations/
│   └── yeast_go_annotations.tsv.gz
├── netseq/
│   ├── wt_plus.bw                  # WT NET-seq (+ strand)
│   └── wt_minus.bw                 # WT NET-seq (- strand)
└── atract_sites/
    └── yeast_atract_cpa_sites.tsv.gz  # 64K pre-computed CPA positions
```

---

## Chromosome naming

Bundled data uses UCSC-style chromosome names: `chrI`, `chrII`, ..., `chrXVI`, `chrM`.

RECTIFY auto-converts between UCSC and NCBI formats:

| Format | Example |
|--------|---------|
| UCSC (preferred) | `chrI` |
| NCBI RefSeq | `ref\|NC_001133\|` |
| Short | `I` |

If your BAM uses a different naming scheme, RECTIFY normalizes chromosome names internally. Check `processing_stats.tsv` for any chromosome name warnings.

---

## NET-seq data

The bundled WT NET-seq is used by `rectify correct` for **A-tract ambiguity resolution** when `--Scer` is set. It represents WT (BY4742) nascent RNA 3' end positions, deconvolved from oligo(A)-spreading artifacts.

For mutant conditions or custom NET-seq, provide your own BigWigs:

```bash
rectify correct reads.bam \
    --Scer \
    --netseq-dir /path/to/mutant_netseq/ \
    -o results/
```

---

## Other organisms

For non-yeast organisms, provide all reference files explicitly:

```bash
rectify run reads.fastq.gz \
    --genome /path/to/genome.fa.gz \
    --annotation /path/to/genes.gff.gz \
    --go-annotations /path/to/go.tsv.gz \
    -o results/
```

NET-seq refinement is skipped if no NET-seq data is provided.
