# RECTIFY — Project Status

**Version:** 2.7.8 (released 2026-04-09)
**In-development:** 2.7.9 (post-tag commits pending release)
**Last updated:** 2026-04-11

---

## What rectify does

RECTIFY is a framework for correcting systematic 3' and 5' end positioning
errors in nanopore direct RNA-seq data. Errors arise from poly(A) homopolymer
ambiguity, AG mispriming, and sequencing indels near cleavage/polyadenylation
(CPA) sites. The corrected positions enable accurate CPA site mapping,
differential CPA usage analysis, and splice junction characterization across
conditions.

---

## Distribution

| Channel | Version | URL |
|---------|---------|-----|
| PyPI | 2.7.8 | https://pypi.org/project/rectify-rna/2.7.8/ |
| Anaconda.org | 2.7.8 (noarch) | https://anaconda.org/kevinrjroy/rectify-rna |
| GitHub | tag `v2.7.8` | https://github.com/k-roy/RECTIFY |

```bash
pip install rectify-rna==2.7.8
conda install -c kevinrjroy -c conda-forge -c bioconda rectify-rna
```

---

## Test suite

**472 tests passing** (0 failures, 2 warnings) — ~81 s
**Coverage:** 25% overall; core pipeline modules ~40–60%

Coverage gaps: `visualize/` (ridge, vep_panels, read_browser) at 0–16%.

---

## Core pipeline — module status

| Module | Status | Notes |
|--------|--------|-------|
| A-tract ambiguity correction | Stable | NET-seq-guided repositioning |
| Poly(A) soft-clip rescue | Stable | Sequence-based (v2.7.5+) |
| Indel artifact correction | Stable | Preserves prior corrections (NEW-001, v2.7.8) |
| AG mispriming detection | Stable | |
| Variant-aware rescue | Stable | |
| NET-seq refinement | Stable | Bundled S. cerevisiae WT data |
| Spike-in filter | Stable | |
| N-base snapping | Stable | NEW-061 fixed (v2.7.9-dev) |
| Corrected BAM output | Stable | `--write-corrected-bam` (v2.5.0+) |
| 5' end correction | Stable | Sequence-based rescue; minus-strand fix NEW-002 |
| Full-length classifier | Stable | Minus-strand coordinate fix NEW-060 (v2.7.8) |

### Multi-aligner consensus

Three aligners run in parallel (minimap2, mapPacBio/BBMap, gapmm2) and vote
on the best alignment per read. Junction tuples now carry strand (NEW-020);
minus-strand donor/acceptor assignments corrected throughout (NEW-011/012).
mapPacBio subprocess handling hardened: timeout, stderr=DEVNULL, returncode
checks (NEW-004/005/021).

### Streaming memory pipeline

Two-pass manifest mode (`rectify analyze --manifest`) keeps peak RAM at
O(clusters × samples) regardless of read depth. Validated at 21 samples /
150 M reads on a 16 GB node. Position index files
(`corrected_3ends_index.bed.gz`) generated alongside every
`corrected_3ends.tsv` for near-instant second-pass loading.

### Analysis outputs

| Output | Status |
|--------|--------|
| CPA clusters (`cpa_clusters.tsv`) | Stable |
| Count matrix (`cluster_counts.tsv`) | Stable; zero-sample pre-init (NEW-049) |
| Bedgraphs | Stable; atomic writes (NEW-050) |
| PCA / heatmap / shift analysis | Stable |
| DESeq2 gene + cluster results | Stable; condition guard added (NEW-026) |
| GO enrichment | Stable |
| De novo motif discovery | Stable; exceptions logged (NEW-034) |
| 5' / 3' / transcript-body genomic distribution | Stable (v2.6.0+) |
| HTML report | Stable |

---

## Bundled data (S. cerevisiae)

| Resource | Source | Notes |
|----------|--------|-------|
| Genome (R64-5-1) | SGD 2024-05-29 | bgzip + samtools faidx indexed |
| GFF annotation | SGD R64-5-1 | SGD-curated ncRNAs; functionally characterized CUT/SUT/XUT only |
| GO annotations | SGD R64-5-1 | |
| WT NET-seq | Churchman 2011 / Harlen 2016 | Pre-processed TSV |
| TF motif database | CPA factors + NNS pathway | MEME format |

**Full CUT/SUT/XUT sets** (not in bundled GFF — from original papers):
- `common/annotation_files/non_coding_rna/CUTs_Xu2009.gff` — 925 loci
- `common/annotation_files/non_coding_rna/SUTs_Xu2009.gff` — 847 loci
- `common/annotation_files/non_coding_rna/XUTs_VanDijk2011.gff` — 1658 loci

---

## Release history (2026)

| Version | Date | Highlights |
|---------|------|-----------|
| **2.7.9-dev** | unreleased | NEW-061–064: N-snap, 5' TSV column, hard-clip walkback, netseq flag |
| **2.7.8** | 2026-04-09 | 35 Round 2 audit bugs; strand fixes throughout; shell injection fix |
| **2.7.7** | 2026-04-08 | 17 Round 1 audit bugs; 3'SS off-by-one; JSD sqrt; SGE 1-based fix |
| **2.7.6** | 2026-04-03 | NetSeqSignal pickle; BAM handle leak; false junction filter logger |
| **2.7.5** | 2026-04-03 | Sequence-based 5' soft-clip rescue; junction filtering; genome cache |
| **2.6.0** | 2026-04-03 | 5' genomic distribution analysis; dinucleotide 3'SS validation |
| **2.5.0** | 2026-04-02 | `--write-corrected-bam`; NET-seq bedgraph export |
| **2.4.0** | 2026-04-01 | Full-length classifier; VEP panels; read browser |
| **2.3.0** | 2026-03-30 | Streaming two-pass manifest analysis; position index files |
| **2.2.0** | 2026-03-28 | Bundled yeast genome/annotation; FASTQ input support |
| **2.1.0** | 2026-03-16 | Multi-aligner consensus (mapPacBio + gapmm2); HPC batch mode |
| **2.0.0** | 2026-03-09 | Complete rewrite: streaming pipeline, DESeq2, motif discovery |

---

## In-development (post-v2.7.8, not yet tagged)

Four bugs fixed from a validation audit — staged in master, tag `v2.7.9` and
PyPI/conda publish not yet done:

| Bug | Module | Description |
|-----|--------|-------------|
| NEW-061 | `indel_corrector.py` | N-absorption: N CIGAR ops shifted corrected position into poly-A region |
| NEW-062 | `bam_processor.py` | `five_prime_rescued` column missing from streaming-mode TSV output |
| NEW-063 | `consensus.py` | Hard-clip walkback incorrectly applied to chimeric reads |
| NEW-064 | `netseq_bam_processor.py` | NET-seq flag guard missing for unstranded BAMs |

---

## Open bugs (tracked in docs/BUGS_TO_FIX.md)

| ID | Priority | Description |
|----|----------|-------------|
| Bug 37 | HIGH | Zero unit tests for `terminal_exon_refiner.py` (1690 lines) |
| Bug 38 | HIGH | `consensus.py` alignment selection only indirectly tested |
| Bug 41 | MEDIUM | `--polya-model` now forwarded but BAM processor doesn't consume it end-to-end |
| Bug 55 | MEDIUM | APA clustering params `--min-peak-sep`, `--max-cluster-radius`, `--min-cluster-samples` not CLI-configurable |

---

## Active experiment: roadblocks project

21-sample nanopore run analyzing transcription termination roadblocks in
*S. cerevisiae*. All per-sample corrections complete.

**Location:** `projects/roadblocks/rectify_alpha_output/minimap2_run/`

### Analysis split

| Task | Samples | Background | Output dir |
|------|---------|------------|------------|
| 1 | rna15_rep1–3, wt_rep1–3, ysh1_rep1–3 | WT S288C | `combined_cpa/` |
| 2 | wt_by4742_rep1–3, dst1d_rep1–3, ski7d_rep1–3, 4nqo_rep1–3 | BY4742 | `combined_by4742/` |

Script: `minimap2_run/run_combined_analysis.sh` (manifest mode, 16 GB RAM)

---

## Key file locations

| Resource | Path |
|----------|------|
| rectify source | *(your local clone)* |
| roadblocks run dir | `projects/roadblocks/rectify_alpha_output/minimap2_run/` |
| CUT/SUT/XUT GFFs | `common/annotation_files/non_coding_rna/` |
| Generic SLURM profile | `rectify/slurm_profiles/slurm_generic.yaml` |
| Bundled annotation | `rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz` |
