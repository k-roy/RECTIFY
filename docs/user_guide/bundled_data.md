# Bundled Data (Yeast)

For *Saccharomyces cerevisiae*, RECTIFY bundles all reference data required for a complete run. No external files are needed.

---

## Activating bundled data

Use either flag:

```bash
--Scer              # Shorthand
--organism yeast    # Long form (also accepts: saccharomyces_cerevisiae)
```

Both flags are recognized by every reference-aware subcommand: `align`, `analyze`, `aggregate`, `extract`, `validate`, `netseq`, `consensus`, `train-polya`, `prescan`, `split`, `export`, `correct`, `run-all`, and `batch`. When set, `--genome` and `--annotation` become optional — the bundled paths are auto-filled.

---

## What is bundled

| Data | Source | Notes |
|------|--------|-------|
| **Genome** | *S. cerevisiae* S288C R64-5-1 (SGD) | Bgzipped FASTA (`.fsa.gz` + `.fai` + `.gzi`); all 16 chromosomes + mitochondria. NCBI-style chromosome names available in a parallel `.fsa.ncbi.gz`. |
| **Gene annotation** | SGD + ncRNA additions (CUT/SUT/XUT) | Bgzipped GFF3 (`.gff.gz`); ~6 K protein-coding genes plus 925 CUTs, 847 SUTs, 1658 XUTs |
| **GO annotations** | SGD | Tab-separated, gzipped — used by `rectify analyze --go-annotations` for enrichment |
| **NET-seq tables** | Churchman & Weissman 2011 (Nature) and Harlen et al. 2016 (Cell) | Three pre-aggregated TSV.GZ tables (`netseq_pan` = 6 reps WT + DST1D; `netseq_wt` = 3 reps WT; `atract_netseq` = A-tract preview) used for A-tract ambiguity resolution |
| **Motif database** | Curated yeast CPA factors + NNS pathway + general TFs | MEME-format file used by `rectify analyze --run-motif` |
| **Validation reads** | 33 hand-curated reads (DRS + cDNA + QuantSeq REV) | BAM + FASTQ in `rectify/data/validation/` — used by `rectify test` smoke check |
| **deSALT binary** | v1.5.6 Linux/x86_64, vendored | Used automatically when `deSALT` is not on `PATH` (Linux/x86_64 only) |

### Bundled data layout

Data lives in the `rectify/data/` directory inside the installed package:

```
rectify/data/
├── __init__.py                                         # Loaders: ensure_reference_data,
│                                                       #          get_bundled_genome_path,
│                                                       #          resolve_reference_paths, ...
├── saccharomyces_cerevisiae_netseq_pan.tsv.gz         # Pan-mutant NET-seq table
├── saccharomyces_cerevisiae_netseq_wt.tsv.gz          # WT-only NET-seq table
├── saccharomyces_cerevisiae_atract_netseq.tsv.gz      # A-tract preview NET-seq table
├── motif_databases/
│   └── scerevisiae_tf_motifs.meme                     # CPA factors + NNS + TFs (MEME format)
├── genomes/saccharomyces_cerevisiae/
│   ├── S288C_reference_sequence_R64-5-1_20240529.fsa.gz      # Bgzipped genome (chrI..chrXVI + chrM)
│   ├── S288C_reference_sequence_R64-5-1_20240529.fsa.gz.fai  # samtools FASTA index
│   ├── S288C_reference_sequence_R64-5-1_20240529.fsa.gz.gzi  # bgzip block index
│   ├── S288C_reference_sequence_R64-5-1_20240529.fsa.ncbi.gz # NCBI-named copy
│   ├── S288C_reference_sequence_R64-5-1_20240529.fsa.pkl     # Pre-loaded genome dict (fast access)
│   ├── saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz      # SGD + ncRNA GFF3
│   ├── go_annotations.tsv.gz                          # GO annotations
│   ├── bbmap_index/                                   # mapPacBio (BBMap) index (generated lazily)
│   └── desalt_index/                                  # deSALT index (generated lazily)
├── validation/                                         # 33-read validation set (DRS + cDNA + QuantSeq REV)
│   ├── validation_reads.bam
│   ├── validation_reads_cdna.bam
│   ├── validation_reads_quantseq_rev.bam
│   ├── corrected_reads.tsv
│   └── aligners/                                       # Per-aligner BAMs for consensus regression tests
└── bin/
    └── linux_x86_64/
        ├── deSALT                                       # Vendored deSALT v1.5.6 binary
        └── namfinder                                    # Helper used by deSALT
```

To reference bundled data paths programmatically, use the helpers in
`rectify.data` rather than hardcoded paths — the layout may evolve:

```python
from rectify.data import (
    get_bundled_genome_path,
    get_bundled_annotation_path,
    get_bundled_go_annotations_path,
    get_bundled_netseq_path,
    get_motif_database_path,
)

GENOME     = get_bundled_genome_path('saccharomyces_cerevisiae')        # Path to .fsa.gz
ANNOTATION = get_bundled_annotation_path('saccharomyces_cerevisiae')    # Path to .gff.gz
GO         = get_bundled_go_annotations_path('saccharomyces_cerevisiae')
NETSEQ     = get_bundled_netseq_path('saccharomyces_cerevisiae')        # netseq_pan by default
MOTIFS     = get_motif_database_path('saccharomyces_cerevisiae', 'cpa_factors')
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
