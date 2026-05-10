# RECTIFY — Project Status

**Version:** 3.2.5 (released 2026-04-24)
**Last updated:** 2026-04-30

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
| PyPI | 3.2.5 | https://pypi.org/project/rectify-rna/3.2.5/ |
| Anaconda.org | 3.2.5 (noarch) | https://anaconda.org/kevinrjroy/rectify-rna |
| GitHub | tag `v3.2.5` | https://github.com/k-roy/RECTIFY |

```bash
pip install rectify-rna==3.2.5
conda install -c kevinrjroy -c conda-forge -c bioconda rectify-rna
```

---

## Test suite

**708 tests passing** (0 failures, 4 skipped) — ~120 s
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

**Tier 1 (default):** minimap2, mapPacBio/BBMap, and gapmm2 run in parallel and vote
on the best alignment per read. **Tier 2 (opt-in):** deSALT and uLTRA add additional
splice-site sensitivity. Short-read mode (`--short-read`): bbmap and bwa replace
the Tier 1 long-read panel for Illumina/Aviti ≤150 bp reads. Junction tuples
now carry strand (NEW-020); minus-strand donor/acceptor assignments corrected
throughout (NEW-011/012). mapPacBio subprocess handling hardened: timeout,
stderr=DEVNULL, returncode checks (NEW-004/005/021).

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
| **3.2.5** | 2026-04-24 | Validation Cat1/Cat2 replaced with DRS-trimmed examples; 708 tests |
| **3.2.2** | 2026-04-22 | XV/XG tag fix; N-op fix for aligner BAMs; 36/36 validation reads certified |
| **3.2.1** | 2026-04-22 | `rescue_3ss_truncation` soft-clip exon CIGAR body-borrowing fix |
| **3.2.0** | 2026-04-22 | Validation aligner BAMs rebuilt after DRS rebuild |
| **3.1.9** | 2026-04-22 | Step 4 `rectified_pA_tail_soft_clipped.bam`; sort+index both paths |
| **3.1.8** | 2026-04-22 | Validation Cat6/Cat7/Cat9 updated to DRS pre-trim mapPacBio alignments |
| **3.1.7** | 2026-04-21 | Module 2H bilateral t2, no-candidate-guard policy, canonical HP prior |
| **3.1.6** | 2026-04-21 | `--drs` flag wired into `run-all`; ARCHITECTURE.md updated |
| **3.1.5** | 2026-04-21 | `--checkpoint-dir` two-level checkpoint/resume for `rectify correct` |
| **3.1.4** | 2026-04-21 | Module 2H candidate guard + adaptive tie-break |
| **3.1.3** | 2026-04-21 | `--aligner-bams` `aligner:path` prefix stripping fix |
| **3.1.2** | 2026-04-21 | Cat9 validation reads added (Module 2H junction refinement) |
| **3.1.1** | 2026-04-21 | `_score_junction` range(L) fix + `is_alt` tiebreaker |
| **3.1.0** | 2026-04-20 | Module 2H — post-consensus N-op junction refinement |
| **3.0.4** | 2026-04-20 | `rescue_3ss_truncation` minus-strand soft-clip truncation fix |
| **3.0.3** | 2026-04-16 | `find_polya_boundary` trailing-base false-stop guard |
| **3.0.2** | 2026-04-16 | `clip_read_to_corrected_3prime` terminal D/N stripping |
| **3.0.1** | 2026-04-15 | `clip_intronic_tail_5prime` off-by-one + trailing I/S stripping |
| **3.0.0** | 2026-04-15 | `clip_intronic_tail_5prime` BAM sequence trimming fix |
| **2.9.9** | 2026-04-15 | HP-aware edit distance + 3'SS acceptor tiebreaker |
| **2.9.8** | 2026-04-15 | Two-phase discovery + canonical refinement, ±5 bp baseline |
| **2.9.7** | 2026-04-14 | Natural sort K-way merge fix in `consensus.py` |
| **2.9.6** | 2026-04-14 | Splice-site ambiguity resolution via data-driven shift range |
| **2.9.5** | 2026-04-14 | Case 4 intronic-snap rescue in `rescue_3ss_truncation` |
| **2.9.4** | 2026-04-14 | N-op boundary guard for spliced minus-strand reads |
| **2.9.3** | 2026-04-14 | Large-deletion pre-scan for poly-A over-calling artifacts |
| **2.9.2** | 2026-04-14 | Cat3 5' rescue mapPacBio intrusion fix |
| **2.9.1** | 2026-04-12 | Cat2 soft-clip rescue CIGAR surgery |
| **v0.9.0-dev** | 2026-04-11 | Bugs 37/38/41/55 fixed; `rectify split/consensus/install-aligners` |
| **2.8.1** | 2026-04-11 | Module 2G standalone; `--dT-primed-cDNA` rename |
| **2.8.0** | 2026-04-11 | Cat3 semi-global NW local aligner; chimeric stitch D/N fix |
| **2.7.9** | 2026-04-09 | NEW-061–064: N-snap, 5' TSV column, hard-clip walkback, netseq flag |
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

## Open bugs (tracked in docs/BUGS_TO_FIX.md)

All previously tracked HIGH-priority bugs (37, 38, 41, 55) resolved in v0.9.0-dev / v3.x series.

Current open items are MEDIUM priority (8 items, NEW-067 through NEW-074) — see
`docs/BUGS_TO_FIX.md` for full details. No CRITICAL or HIGH open issues.

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
