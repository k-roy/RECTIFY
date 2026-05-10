# Bundled Data (Yeast)

For *Saccharomyces cerevisiae*, RECTIFY bundles all reference data required for a complete run. No external files are needed.

---

## Activating bundled data

Use either flag:

```bash
--Scer              # Shorthand
--organism yeast    # Long form (also accepts: saccharomyces_cerevisiae)
```

Both flags are recognized by `rectify correct`, `rectify run-all`, `rectify align`, `rectify analyze`, and `rectify netseq`.

---

## What is bundled

| Data | Source | Notes |
|------|--------|-------|
| **Genome** | *S. cerevisiae* S288C R64-5-1 | FASTA (`.fsa`); includes all 16 chromosomes + mitochondria |
| **Gene annotations** | SGD (Saccharomyces Genome Database) | GFF and GTF formats; 6,000+ genes, UTRs, introns |
| **Junction BED** | Derived from annotations | Pre-built for minimap2 `--junc-bed`; improves splice accuracy |
| **WT NET-seq** | Published dataset | Pre-aggregated TSV tables (A-tract, pan, and WT-only); used for A-tract ambiguity resolution |
| **Poly(A) model** | Trained on yeast WT | Model weights; used by `rectify correct` for poly(A) scoring |

### Bundled data location

Data lives in the `rectify/data/` directory inside the installed package:

```
rectify/data/
├── S288C_reference_sequence_R64-5-1_20240529.fsa       # Genome FASTA (S288C R64-5-1)
├── S288C_reference_sequence_R64-5-1_20240529.fsa.fai   # samtools FASTA index
├── S288C_reference_sequence_R64-5-1_20240529.pkl       # Pre-loaded genome dict (fast access)
├── saccharomyces_cerevisiae_R64-5-1_20240529.gff       # Gene annotations (GFF)
├── saccharomyces_cerevisiae_R64-5-1_20240529.gtf       # Gene annotations (GTF)
├── saccharomyces_cerevisiae_R64-5-1_20240529.junc.bed  # Pre-built junction BED for minimap2 --junc-bed
├── saccharomyces_cerevisiae_atract_netseq.tsv.gz       # A-tract NET-seq signal table
├── saccharomyces_cerevisiae_netseq_pan.tsv.gz          # Pan NET-seq table
├── saccharomyces_cerevisiae_netseq_wt.tsv.gz           # WT NET-seq table
├── models/                                              # Poly(A) model weights
├── motif_databases/                                     # JASPAR-format motif databases
├── genomes/saccharomyces_cerevisiae/
│   └── saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz  # GFF.GZ copy (for rectify prescan)
└── bin/
    └── linux_x86_64/
        └── deSALT                                       # Vendored deSALT v1.5.6 binary
```

To reference bundled data paths programmatically:

```python
import rectify
from pathlib import Path
DATA = Path(rectify.__file__).parent / 'data'

GENOME     = DATA / 'S288C_reference_sequence_R64-5-1_20240529.fsa'
ANNOTATION = DATA / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff'
JUNC_BED   = DATA / 'saccharomyces_cerevisiae_R64-5-1_20240529.junc.bed'
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

The bundled WT NET-seq is used by `rectify correct` for **A-tract ambiguity resolution** when `--Scer` is set. It represents WT (BY4742) nascent RNA 3' end positions in pre-aggregated TSV table format (`saccharomyces_cerevisiae_netseq_wt.tsv.gz`, `saccharomyces_cerevisiae_netseq_pan.tsv.gz`, `saccharomyces_cerevisiae_atract_netseq.tsv.gz`).

For mutant conditions or custom NET-seq, provide your own BigWig directory:

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
rectify run-all reads.fastq.gz \
    --genome /path/to/genome.fa.gz \
    --annotation /path/to/genes.gff.gz \
    --go-annotations /path/to/go.tsv.gz \
    -o results/
```

NET-seq refinement is skipped if no NET-seq data is provided.
