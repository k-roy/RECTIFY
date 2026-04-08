# RECTIFY Architecture

## What is this document?

This is a narrative architecture guide for developers, collaborators, and
automated agents. It covers the pipeline design, module responsibilities,
and data flow at a level deeper than the README but focused on the "why" —
not a function-by-function reference.

For auto-generated API reference (every public function, parameter, and
return type), see the ReadTheDocs site built from docstrings:
[rectify-rna.readthedocs.io](https://rectify-rna.readthedocs.io).
RTD answers *what does this function do?* — this document answers
*how do the modules fit together?*

---

## Pipeline overview

RECTIFY corrects systematic errors in poly(A)-tailed RNA-seq data (nanopore
Direct RNA-seq, QuantSeq, Helicos, PacBio Iso-Seq) and then performs
differential 3' end usage analysis across conditions.

```
FASTQ / BAM
    │
    ▼
[Step 0: Alignment]        multi-aligner rectification pipeline
    │                      Tier 1 (default): minimap2 + mapPacBio + gapmm2
    │                      Tier 2 (opt-in):  + deSALT, + uLTRA
    │                      All enabled aligners run in parallel subprocesses.
    │                      Chimeric consensus selection per read:
    │                        find sync-points where aligners agree →
    │                        score each inter-sync segment independently →
    │                        stitch best segment from each aligner
    ▼
    <sample>.rectified.bam
    │
    ▼
[Step 1: Correction]       per-read 3' and 5' end correction
    │                      ① A-tract ambiguity detection (universal)
    │                      ② AG mispriming screen (oligo-dT only)
    │                      ③ Poly(A) tail trimming (if sequenced)
    │                      ④ Indel artifact correction
    │                      ⑤ False junction filtering (poly(A)-artifact N-ops)
    │                      ⑥ 5' soft-clip junction rescue
    │                      ⑦ NET-seq refinement (optional, resolves ambiguous windows)
    │                      ⑧ Spike-in read filtering
    ▼
    corrected_3ends.tsv           per-read corrected positions; each row is one read
                                  columns include: chrom, strand, original_3prime,
                                  corrected_3prime, five_prime_position, polya_length,
                                  junctions_str, n_junctions, confidence, …
    corrected_3ends_index.bed.gz  position-count summary: one row per unique
                                  (chrom, corrected_3prime, strand) with read count.
                                  ~300× smaller than the per-read TSV; used by
                                  manifest-mode analysis to skip re-reading every read.
    <sample>.stats.tsv            per-sample QC report (reads processed, corrections
                                  applied per module, confidence distribution)
    │
    ▼
[Step 2: Analysis]         multi-sample downstream analysis
    │                      ① CPA site clustering
    │                      ② Gene attribution (which gene owns each cluster?)
    │                      ③ DESeq2 differential usage (multi-sample only):
    │                           — cluster level: per-CPA-site read counts
    │                           — gene level: summed CPA counts across all clusters
    │                           each produces deseq2_{clusters,genes}_{condition}.tsv
    │                      ④ 3' UTR shift analysis: gene-level weighted-mean CPA
    │                           position shift in bp between conditions — detects
    │                           global 3' UTR lengthening/shortening per gene
    │                      ⑤ APA isoform detection: read-level grouping by
    │                           (gene, junction signature, 3' cluster) to identify
    │                           distinct isoforms and proximal:distal usage ratios
    │                      ⑥ De novo motif discovery: STREME/MEME on sequences
    │                           around DESeq2-enriched (padj<0.05, log2FC>1) and
    │                           -depleted (log2FC<-1) CPA clusters — identifies
    │                           differential poly(A) signal or downstream elements
    │                      ⑦ GO enrichment: Fisher's exact test on genes with
    │                           significantly increased or decreased CPA cluster
    │                           usage (padj<0.05, |log2FC|>1 at gene level)
    │                      ⑧ Genomic region distribution (two analyses):
    │                           — 3' end distribution: classify corrected 3' end
    │                             positions by feature (3'UTR > snoRNA > CUT >
    │                             SUT/XUT > 5'UTR/CDS > antisense > intergenic)
    │                           — Transcript body distribution: classify each read
    │                             by the biotype of the feature its full alignment
    │                             span (alignment_start → alignment_end) overlaps most
    ▼
    cpa_clusters.tsv
    deseq2_genes_{condition}.tsv
    deseq2_clusters_{condition}.tsv
    shift_analysis_{condition}.tsv
    go_enrichment_{up,down}_{condition}.tsv
    plots/
```

Steps 0–2 are orchestrated by `rectify run-all`. Each step can also be
invoked independently through its own subcommand.

---

## CLI subcommand dispatch

Entry point: `rectify.cli:main` → `create_parser()` → per-subcommand
`create_*_parser()` → dispatches to the corresponding `*_command.py`.

| Subcommand | Module | Purpose |
|---|---|---|
| `run-all` | `core/run_command.py` | Full end-to-end pipeline (Steps 0–2) |
| `align` | `core/align_command.py` | Multi-aligner alignment only (Step 0) |
| `correct` | `core/correct_command.py` | 3' end correction only (Step 1) |
| `analyze` | `core/analyze_command.py` | Downstream analysis only (Step 2) |
| `batch` | `core/batch_command.py` | Parallel correction across samples; interactive mode auto-sizes to available CPUs, HPC mode generates array job scripts for SLURM, PBS/Torque, or UGE/SGE via a portable scheduler abstraction |
| `train-polya` | `core/train_polya_command.py` | Train a poly(A) tail model from calibration data (see below) |
| `validate` | `core/validate_command.py` | Post-correction quality check against NET-seq or known CPA sites (see below) |
| `export` | `core/export_command.py` | Export corrected 3' ends to bedGraph/bigWig |
| `extract` | `core/extract_command.py` | Extract per-read info from BAM to TSV |
| `netseq` | `core/netseq_command.py` | Process NET-seq BAM files with deconvolution |
| `aggregate` | `core/aggregate_command.py` | Aggregate reads into CPA / TSS / junction datasets |

> **Note on `train-polya`:** Run once per sequencing technology/chemistry to
> fit A-richness thresholds. Operates on reads mapping to control sites that
> have zero downstream genomic A's — any A-rich soft-clipped sequence at such
> sites must be a genuine poly(A) tail, not a genomic artifact. Outputs a JSON
> model loaded by `rectify correct --polya-model`. Bundled models for common
> platforms are included; re-run only when adapting to a new chemistry.

> **Note on `validate`:** Optional post-correction QC. Compares corrected 3'
> end positions against NET-seq bigWig files or a known-CPA-positions TSV and
> reports accuracy improvement over raw positions. Used during development and
> benchmarking; not needed for routine production runs.

> **Note on `run_all_command.py`:** This file contains an alternative
> orchestrator (align → correct → aggregate → export → analyze) with
> provenance-based step-skipping built into the chain. It is **not wired to
> the CLI** — the active `run-all` dispatcher is `run_command.py`. The file
> should either be promoted to replace `run_command.py` (once the aggregate
> step is validated end-to-end) or deleted. Keeping both is a maintenance
> liability.

---

## Directory structure

```
rectify/                              ← git repo root
├── rectify/                          ← importable Python package
│   │
│   ├── __init__.py                   version string, public API re-exports
│   ├── __main__.py                   allows `python -m rectify`
│   ├── cli.py                        argparse entry point; dispatches subcommands
│   ├── config.py                     all constants (chroms, thresholds, shifts)
│   ├── slurm.py                      HPC utilities: thread limits, scratch staging
│   ├── slurm_profiles/               YAML configs for Sherlock partitions
│   │   ├── sherlock_larsms.yaml      CPU partition (streaming + scratch on by default)
│   │   └── sherlock_gpu.yaml         GPU partition
│   ├── provenance.py                 SHA-256 provenance tracking; skip-if-unchanged logic
│   │
│   ├── core/                         pipeline step implementations
│   │   │
│   │   ├── run_command.py            Step 0+1+2 orchestrator (the "run-all" dispatcher)
│   │   ├── run_all_command.py        ⚠️ experimental redesign — NOT active
│   │   ├── align_command.py          Step 0 CLI wrapper
│   │   ├── correct_command.py        Step 1 CLI wrapper
│   │   ├── analyze_command.py        Step 2 CLI wrapper + GFF/GTF parsing
│   │   ├── batch_command.py          parallel/SLURM batch correction
│   │   │
│   │   ├── multi_aligner.py          Tier 1: minimap2+mapPacBio+gapmm2; Tier 2: +deSALT+uLTRA
│   │   ├── chimeric_consensus.py     chimeric alignment stitching from sync-points
│   │   ├── consensus.py              per-read optimal aligner selection (non-chimeric fallback)
│   │   ├── bam_processor.py          correction pipeline orchestrator (Steps 1①–⑧)
│   │   │
│   │   ├── atract_detector.py        A-tract boundary detection + ambiguity window
│   │   ├── ag_mispriming.py          AG-richness mispriming screen (Roy & Chanfreau 2019)
│   │   ├── polya_trimmer.py          poly(A) tail detection and trimming
│   │   ├── polya_model.py            JSON poly(A) model (A-richness thresholds)
│   │   ├── indel_corrector.py        indel artifact correction; variant-aware rescue
│   │   ├── false_junction_filter.py  poly(A)-artifact N-op (splice junction) removal
│   │   ├── splice_aware_5prime.py    5' soft-clip junction rescue (intron-displaced TSS)
│   │   ├── netseq_refiner.py         NET-seq NNLS deconvolution for CPA localization
│   │   ├── netseq_deconvolution.py   low-level NNLS / PSF math
│   │   ├── netseq_command.py         `rectify netseq` CLI
│   │   ├── netseq_bam_processor.py   NET-seq BAM → 3' end TSV
│   │   ├── netseq_output.py          NET-seq output formatting
│   │   ├── terminal_exon_refiner.py  terminal exon boundary refinement
│   │   ├── spikein_filter.py         spike-in construct detection and filtering
│   │   ├── exclusion_regions.py      blacklist regions (repetitive elements, etc.)
│   │   ├── junction_validator.py     cross-sample COMPASS-style junction validation
│   │   ├── mpb_split_reads.py        mapPacBio long-read splitting and stitching
│   │   ├── preprocess.py             input detection (FASTQ vs BAM), bundled genome prep
│   │   ├── unified_record.py         unified read record dataclass
│   │   ├── processing_stats.py       per-sample QC stat accumulation and reporting
│   │   ├── export_command.py         bedGraph/bigWig export
│   │   ├── extract_command.py        per-read BAM → TSV extraction
│   │   ├── validate_command.py       correction validation vs NET-seq / known sites
│   │   ├── train_polya_command.py    poly(A) model training
│   │   ├── aggregate_command.py      3'/5'/junction aggregation dispatcher
│   │   │
│   │   ├── analyze/                  downstream statistical analysis modules
│   │   │   ├── clustering.py         CPA site clustering (fixed-distance + adaptive)
│   │   │   ├── gene_attribution.py   read-body attribution of 3' clusters to genes
│   │   │   ├── deseq2.py             DESeq2 differential cluster + gene usage
│   │   │   ├── shift_analysis.py     3' UTR shift scoring and visualization
│   │   │   ├── apa_detection.py      APA isoform quantification (Isosceles-style)
│   │   │   ├── motif_discovery.py    de novo motif discovery via STREME/MEME
│   │   │   ├── go_enrichment.py      GO enrichment via Fisher's exact test
│   │   │   ├── genomic_distribution.py   read classification by genomic feature
│   │   │   ├── junction_analysis.py  splice junction statistics
│   │   │   ├── junction_validation.py    annotation-based junction validation
│   │   │   ├── deconvolution.py      generalized deconvolution utilities
│   │   │   ├── pan_mutant_refiner.py pan-mutant NET-seq CPA refinement
│   │   │   ├── atract_refiner.py     A-tract-aware CPA position refinement
│   │   │   ├── heatmap.py            heatmap visualization of cluster usage
│   │   │   ├── pca.py                PCA of per-sample cluster counts
│   │   │   └── summary.py            per-run summary statistics
│   │   │
│   │   ├── aggregate/                per-end aggregation (output of `rectify aggregate`)
│   │   │   ├── three_prime.py        cluster by 3' CPA end; attribute to gene via read body
│   │   │   ├── five_prime.py         cluster by 5' TSS end; Bambu-style full-length filter
│   │   │   └── junctions.py          aggregate splice junctions with partial evidence
│   │   │
│   │   └── classify/
│   │       └── full_length_classifier.py   Bambu-style full-length vs truncated classifier
│   │
│   ├── utils/                        shared low-level utilities
│   │   ├── genome.py                 load_genome(), fetch_genomic_sequence(), A-tract utils
│   │   ├── alignment.py              CIGAR parsing, soft-clip extraction, coordinate ops
│   │   ├── chromosome.py             chromosome name normalization helpers
│   │   ├── stats.py                  confidence scoring, QC metrics
│   │   ├── junction_bed.py           minimap2 junction BED generation from annotation
│   │   ├── splice_motif.py           GT-AG and other splice site motif scoring
│   │   └── provenance.py             (re-exports from rectify.provenance)
│   │
│   ├── visualize/                    plotting layer
│   │   ├── config.py                 plot-wide style constants (colors, fonts, DPI)
│   │   ├── figure_utils.py           shared figure helpers (axes, colorbars, layout)
│   │   ├── gene_track.py             gene structure panel (box-arrow CDS glyphs)
│   │   ├── coverage.py               BAM coverage extraction and filled-area plots
│   │   ├── metagene.py               metagene signal aggregation (PositionIndex)
│   │   ├── metagene_loaders.py       convenience loaders for metagene workflows
│   │   ├── multi_track.py            multi-panel genome browser layout
│   │   ├── read_browser.py           stacked read browser (LineCollection-based)
│   │   ├── ridge.py                  ridge/joy-division plots of CPA distributions
│   │   └── vep_panels.py             variant effect prediction visualization panels
│   │
│   ├── data/                         bundled reference data
│   │   ├── __init__.py               data-loading API (ensure_reference_data(), etc.)
│   │   ├── genomes/saccharomyces_cerevisiae/
│   │   │   ├── S288C_reference_sequence_R64-5-1_20240529.fsa.gz   genome FASTA (bgzipped)
│   │   │   ├── *.fai, *.gzi, *.pkl                                 index and pickle cache
│   │   │   ├── S288C_reference_sequence_R64-5-1_20240529.ncbi.gz   NCBI-named version (legacy)
│   │   │   ├── saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz   gene annotation (GFF3)
│   │   │   ├── go_annotations.tsv.gz                               SGD GO term assignments
│   │   │   └── bbmap_index/                                        mapPacBio index
│   │   ├── saccharomyces_cerevisiae_netseq_wt.tsv.gz               WT NET-seq reference
│   │   ├── saccharomyces_cerevisiae_netseq_pan.tsv.gz              pan-mutant NET-seq reference
│   │   ├── saccharomyces_cerevisiae_atract_netseq.tsv.gz           A-tract NET-seq reference
│   │   ├── motif_databases/                                        MEME-format motif databases
│   │   ├── models/                                                  trained poly(A) models (JSON)
│   │   ├── validation_reads.bam / .bai                             test reads for CI/validation
│   │   ├── rectified.bam / .bai                                    expected output for CI
│   │   └── aligners/                                               per-aligner test BAMs for CI
│   │       ├── minimap2.bam, mapPacBio.bam, gapmm2.bam, deSALT.bam, uLTRA.bam
│   │
│   ├── calibration/
│   │   └── calibrate_shift_corrections.py   derive A-count→shift table from NET-seq data
│   │
│   └── scripts/                      one-off database construction scripts
│       ├── build_pan_mutant_netseq_database.py   build pan-mutant NET-seq TSV
│       └── create_atract_netseq_reference.py     build A-tract NET-seq TSV
│
├── tests/                            pytest suite (~25 test files, one per module)
├── docs/                             MkDocs source for RTD
├── conda-recipe/                     conda-forge recipe (meta.yaml)
├── examples/                         usage notebooks and scripts
├── scripts/                          helper shell/Python scripts (not part of package)
├── slurm_profiles/ → rectify/slurm_profiles/  (profiles live inside the package)
├── pyproject.toml                    pip packaging metadata + dependency list
├── mkdocs.yml                        ReadTheDocs / MkDocs configuration
├── ARCHITECTURE.md                   this file
├── CLAUDE.md                         AI agent / developer context
└── README.md                         user-facing overview
```

---

## Layer-by-layer module descriptions

### Layer 1: Entry point

**`cli.py`** — Builds the top-level `argparse` parser and 11 subparsers.
Each subparser delegates to a `create_*_parser()` function in the
corresponding `*_command.py`. No business logic lives here.

**`config.py`** — Single source of truth for all numeric constants:
chromosome name maps (chrI ↔ NCBI), A-count→expected-shift calibration
table, poly(A) detection thresholds, NET-seq signal thresholds, indel
detection windows. Import from here; never hard-code thresholds elsewhere.

**`slurm.py`** — SLURM-aware utilities called at process startup (before
numpy is imported): `set_thread_limits()` sets `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `LOKY_MAX_CPU_COUNT` to
prevent thread oversubscription. Also provides `make_job_scratch_dir()` and
`sync_to_oak()` for scratch-staging patterns.

**`provenance.py`** — Lightweight Snakemake-style provenance: SHA-256
hashes input/output files and writes sidecar `.provenance.json`. The
`step_needed()` method enables skip-if-unchanged logic without Snakemake.

---

### Layer 2: Pipeline orchestration

**`core/run_command.py`** — The full end-to-end orchestrator. Handles
single-sample (FASTQ or BAM input) and multi-sample (manifest TSV) modes.
In multi-sample mode, correction runs per-sample in parallel, then analysis
runs once across all corrected outputs. DESeq2, GO enrichment, and motif
discovery only fire in multi-sample mode (need multiple conditions for
statistics). Also resolves bundled genome/annotation paths transparently.

**`core/batch_command.py`** — Parallel batch correction across samples.
In interactive mode, auto-sizes a `ThreadPoolExecutor` to available CPUs
and runs `rectify correct` per sample in parallel. In HPC mode, reads a
profile YAML and generates array job scripts. Generated scripts include a
portable scheduler abstraction block that normalises CPU count, job ID, and
task ID from SLURM (`$SLURM_CPUS_PER_TASK`), PBS/Torque (`$PBS_NUM_PPN`),
and UGE/SGE (`$NSLOTS`) environment variables, so the same script body
works on all three schedulers without modification.

**`core/align_command.py`** — Thin CLI wrapper around `multi_aligner.py`.

**`core/correct_command.py`** — Step 1 orchestrator: validates inputs,
sets thread limits, calls `bam_processor.process_bam_file_parallel()` or
`process_bam_streaming()`, writes `corrected_3ends.tsv` and stats report.

**`core/analyze_command.py`** — Step 2 orchestrator: parses GFF/GTF,
calls analysis modules in sequence, writes all output TSVs and plots.
Handles both single-sample (no DESeq2) and manifest mode (full analysis).

---

### Layer 3: Alignment (Step 0)

**`core/multi_aligner.py`** — Runs all enabled aligners in parallel
subprocesses. Tier 1 (default): minimap2, mapPacBio, gapmm2. Tier 2
(opt-in via `--aligners`): deSALT (high-sensitivity splice aligner) and
uLTRA (annotation-guided graph aligner, requires GTF). Each aligner
produces a sorted, indexed BAM. Junction annotations from GFF are passed
to minimap2 via `--junc-bed` to improve splice site accuracy while keeping
scoring annotation-blind (novel junctions are still discoverable). Returns
per-aligner BAM paths.

**`core/chimeric_consensus.py`** — Chimeric rectification: finds "sync
points" where two or more aligners agree on query→reference mapping, then
scores segments between sync points independently. The best segment from
each aligner can be combined into one chimeric alignment. This is the key
innovation that handles DRS reads where no single aligner wins everywhere:
one may handle the 5' splice-through better, another the 3' poly(A)
boundary.

**`core/consensus.py`** — Simpler per-read optimal selection (non-chimeric):
scores each aligner's full alignment and picks the best one. Used as
fallback when chimeric stitching is disabled or sync points are not found.

**`core/mpb_split_reads.py`** — mapPacBio can fail on reads >100 kb.
Splits long reads, aligns each chunk, then stitches the resulting BAMs
back together with corrected CIGAR strings.

---

### Layer 4: Correction (Step 1)

**`core/bam_processor.py`** — The correction pipeline hub. Reads the BAM,
dispatches per-read correction through modules ①–⑧ in sequence, and writes
results. Supports two execution modes:
- `process_bam_file_parallel()`: per-chromosome parallelism, accumulates
  all results in memory before writing (use for small genomes or low RAM).
- `process_bam_streaming()`: 10,000-read chunks, writes incrementally.
  Peak RAM ~4–5 GB regardless of BAM size. **Use this for SLURM jobs.**

Correction modules are called in this order:
1. **`atract_detector.py`** — Counts A's (or T's for minus strand) in the
   downstream genomic window, looks up expected alignment shift from the
   calibration table in `config.py`, and computes an `ambiguity_min`/
   `ambiguity_max` window. This correction is UNIVERSAL — it applies to
   all poly(A)-tailed RNA-seq, not just direct RNA.

2. **`ag_mispriming.py`** — For oligo-dT libraries: screens downstream
   sequence for AG-richness. High AG-richness flags a read as likely
   internally primed (misprimed on a genomic A/G run). The original
   RECTIFY algorithm from Roy & Chanfreau 2019.

3. **`polya_trimmer.py`** — For reads that directly sequence the poly(A)
   tail: detects A-rich soft-clipped sequence, calculates `polya_length`,
   and back-calculates the pre-tail cleavage position. Uses a trained
   model (JSON) for A-richness thresholds.

4. **`indel_corrector.py`** — Detects deletion/insertion artifacts that
   arise when aligners force-align poly(A) tails to genomic A-tracts,
   creating spurious CIGAR gaps. The `VariantAwareHomopolymerRescue` class
   also handles cases where a SNP in an A-tract could be misinterpreted
   as an artifact.

5. **`false_junction_filter.py`** — Removes spurious N operations (introns)
   in CIGAR strings near the 3' end that arise from poly(A) tail
   misalignment to distant A-tracts. A junction is flagged if it is within
   a configurable window of the 3' end AND the downstream region is
   highly A-rich.

6. **`splice_aware_5prime.py`** — 5' soft-clip junction rescue. When the
   first exon is very short, aligners may place the read's 5' end (TSS)
   within an annotated intron rather than at the true TSS. This module
   checks whether the 5' end falls inside an intron and, if so, shifts it
   to the nearest upstream exon boundary. Uses both annotated junctions and
   real junctions observed in the read's own CIGAR string.

7. **`netseq_refiner.py`** — For reads whose corrected position falls in
   an ambiguous window: queries NET-seq signal, deconvolves the PSF spread
   caused by short poly(A) tails in NET-seq, and re-assigns reads
   proportionally to deconvolved peak intensities. Converts positional
   uncertainty into a probabilistic split across the most likely true CPA
   sites. Requires NET-seq bigWig or TSV data.

8. **`spikein_filter.py`** — Detects and removes reads from synthetic
   spike-in constructs by comparing 3' UTR sequence to the genomic
   reference; spike-ins have a synthetic 3' UTR not present in the genome.

Support modules called by `bam_processor`:
- **`polya_model.py`** — Loads/saves the JSON poly(A) model.
- **`terminal_exon_refiner.py`** — Refines terminal exon boundaries.
- **`junction_validator.py`** — Cross-sample COMPASS-style junction
  validation: aggregates junctions across samples, applies minimum-support
  and splice-motif filters, then downgrades low-confidence junctions.
- **`exclusion_regions.py`** — Reads a BED blacklist and marks reads
  overlapping excluded regions.

---

### Layer 5: Analysis (Step 2)

**`core/analyze/clustering.py`** — Groups corrected 3' positions into CPA
site clusters. Two algorithms:
- Fixed-distance (default, `--cluster-distance 25`): reads within 25 bp
  on the same strand are merged. O(n log n).
- Adaptive valley-based (`--adaptive-clustering`): finds valleys in the
  read-depth histogram to split clusters at natural breakpoints.
After clustering, `annotate_clusters_with_genes()` uses an IntervalTree to
attribute each cluster to the gene whose 3' region it falls in.

**`core/analyze/gene_attribution.py`** — Attributes 3' end clusters to
genes based on read bodies (not just proximity). For each cluster, tallies
what fraction of reads whose 3' end falls in the cluster also have their
read body overlapping a given CDS/ncRNA feature. This handles the common
case of two adjacent genes where the read bodies disambiguate which gene
the 3' end belongs to.

**`core/analyze/deseq2.py`** — Wraps pyDESeq2 for differential CPA usage.
Runs two independent analyses per condition pair:
- **Cluster level**: per-CPA-site read counts across samples → identifies
  which individual CPA sites are gained or lost between conditions.
  Output: `deseq2_clusters_{condition}.tsv`
- **Gene level**: CPA counts summed across all clusters within each gene →
  identifies genes whose overall 3' end usage changes.
  Output: `deseq2_genes_{condition}.tsv`
Both analyses output padj, log2FoldChange, and baseMean columns. Results
from the cluster-level analysis feed motif discovery; gene-level results
feed GO enrichment.

**`core/analyze/shift_analysis.py`** — Gene-level 3' UTR shift scoring.
For each gene, computes the weighted mean CPA position (weights = read
counts per cluster) in condition A and condition B, and reports the
difference in bp. A positive shift means the 3' ends moved downstream
(3' UTR lengthening); negative means upstream (3' UTR shortening). This
is a single scalar per gene — it summarizes the overall directionality of
CPA redistribution without specifying which isoforms are responsible.

**`core/analyze/apa_detection.py`** — Read-level APA isoform detection
(Isosceles-style). Where shift analysis gives one number per gene, APA
detection identifies *which isoforms* drive the change. Groups reads by
the triplet `(gene, junction_signature, 3' cluster)` — reads with the
same gene body splicing pattern and the same CPA cluster are counted as
one isoform. Computes proximal:distal CPA usage ratio and flags genes
where this ratio changes between conditions.

**`core/analyze/motif_discovery.py`** — De novo motif discovery via STREME
or MEME. Operates on DESeq2 cluster-level results: extracts sequences from
±100/50 bp windows around enriched CPA clusters (padj<0.05, log2FC>1) and
depleted clusters (log2FC<−1), then runs motif discovery on each set. This
identifies poly(A) signals, downstream elements, or other sequence features
that distinguish gained from lost CPA sites. Writes MEME-format output and
a summary TSV per condition.

**`core/analyze/go_enrichment.py`** — GO term enrichment via Fisher's exact
test. Operates on DESeq2 gene-level results: collects genes with
significantly increased CPA usage (padj<0.05, log2FC>1) and genes with
significantly decreased usage (log2FC<−1) separately, then tests for GO
term over-representation in each set. Outputs
`go_enrichment_up_{condition}.tsv` and `go_enrichment_down_{condition}.tsv`.

**`core/analyze/genomic_distribution.py`** — Two complementary distribution
analyses, each producing a horizontal bar chart across conditions:
- **3' end distribution**: classifies each corrected 3' end *position* by
  the genomic feature it falls in. Priority order: 3'UTR > snoRNA±300bp >
  CUTs > SUTs/XUTs > 5'UTR/CDS > antisense CDS > intergenic. Mirrors
  panels B/C of Xu et al. 2009 (Nature 457:1033).
- **Transcript body distribution**: classifies each *read* by the RNA
  biotype of the feature whose bp overlap its full alignment span
  (alignment_start → alignment_end) is greatest. Categories: protein_coding,
  CUT, SUT/XUT, snoRNA, tRNA, rRNA, LTR, pseudogene, antisense, intergenic.

Other analysis modules: `junction_analysis.py` (splice stats),
`junction_validation.py` (annotation-based validation), `heatmap.py`
(cluster usage heatmaps), `pca.py` (PCA of count matrix), `summary.py`
(run-level summary statistics), `pan_mutant_refiner.py` (pan-mutant
NET-seq CPA refinement), `atract_refiner.py` (A-tract-aware refinement),
`deconvolution.py` (generalized NNLS math).

---

### Layer 6: Utilities

**`utils/genome.py`** — `load_genome()`: reads a FASTA (gzipped or plain),
normalizes all chromosome names to chrI/chrII/… format at load time using
`standardize_chrom_name()` (which maps NCBI `ref|NC_001133|` → `chrI`), and
caches as a pickle for fast subsequent loads. All downstream code receives a
`{chrI: sequence, ...}` dict and never needs NCBI-format lookups.
Also provides `fetch_genomic_sequence()`, `reverse_complement()`,
`get_downstream_sequence()`, and A-tract detection helpers.

**`utils/alignment.py`** — CIGAR operation parsing: `extract_junctions_simple()`
(splice-aware N-op extraction), `extract_soft_clips()`, `get_cigar_stats()`,
`format_junctions_string()`. All use 0-based coordinates consistent with pysam.

**`utils/chromosome.py`** — `standardize_chrom_name()` and friends; called
by `load_genome()` at load time so translation is a one-time cost.

**`utils/stats.py`** — `calculate_confidence()` assigns high/medium/low
confidence labels based on the number of corrections applied, degree of
ambiguity, and NET-seq agreement.

**`utils/junction_bed.py`** — Generates the junction BED file from GFF
annotation for `minimap2 --junc-bed`. Cached per sample.

**`utils/splice_motif.py`** — GT-AG, AT-AC, and other splice site motif
classification. Used as a tiebreaker in consensus scoring (deliberately
excluded from primary score to avoid penalizing novel junctions).

---

### Layer 7: Visualization

All plots use matplotlib. Visualization modules are stateless — callers
pass DataFrames and axes, modules draw and return.

**`visualize/metagene.py`** — `PositionIndex`: O(1) read-depth lookups
using a `{chrom: {strand: {position: count}}}` dict. Supports per-locus
normalization, percentile capping, and trimmed-mean aggregation across
loci. `MetagenePipeline` provides a one-call interface for the full
metagene workflow including strand verification.

**`visualize/gene_track.py`** — Draws gene structure panels with pentagon/
arrow CDS glyphs indicating strand. Handles overlapping genes via vertical
staggering. Supports both common names and systematic names.

**`visualize/coverage.py`** — Extracts coverage from BAM files and draws
strand-specific filled-area plots.

**`visualize/read_browser.py`** — Stacked read browser using
`LineCollection` for batch rendering (~6 draw calls for 400 reads vs
2000 `Line2D` objects). Supports junction-spanning reads (exon blocks
connected by thin intron lines).

**`visualize/ridge.py`** — Ridge/joy-division plots of CPA position
distributions per gene or cluster, useful for showing condition-specific
shifts in 3' end usage.

**`visualize/multi_track.py`** — Multi-panel genome browser layout
combining coverage, gene track, and read browser panels with shared
x-axis.

**`visualize/metagene_loaders.py`** — Convenience loaders that build a
`PositionIndex` from corrected TSV files.

**`visualize/figure_utils.py`** — Shared helpers: axis formatting,
colorbar placement, figure sizing, color palette definitions.

**`visualize/vep_panels.py`** — Variant effect prediction visualization:
per-variant effect score panels aligned to genomic tracks.

---

## Bundled data

The package ships with reference data for *Saccharomyces cerevisiae* (S288C
R64-5-1). Other organisms require user-supplied files.

| File | Purpose |
|---|---|
| `genomes/saccharomyces_cerevisiae/*.fsa.gz` | Genome FASTA (bgzipped, indexed) |
| `genomes/saccharomyces_cerevisiae/*.pkl` | Pickle cache for fast genome loading |
| `genomes/saccharomyces_cerevisiae/*.gff.gz` | GFF3 gene annotation |
| `genomes/saccharomyces_cerevisiae/go_annotations.tsv.gz` | SGD GO term assignments |
| `genomes/saccharomyces_cerevisiae/bbmap_index/` | mapPacBio alignment index |
| `netseq_wt.tsv.gz` | Wild-type NET-seq 3' end reference (for refinement) |
| `netseq_pan.tsv.gz` | Pan-mutant NET-seq reference (broader coverage) |
| `atract_netseq.tsv.gz` | A-tract-focused NET-seq reference (highest resolution at A-tracts) |
| `motif_databases/` | MEME-format poly(A) signal motif databases |
| `models/` | Trained poly(A) tail JSON models |
| `validation_reads.bam` | Small test BAM for CI/regression testing |
| `rectified.bam` | Expected output BAM for CI comparison |
| `aligners/*.bam` | Per-aligner test BAMs (minimap2, mapPacBio, gapmm2, deSALT, uLTRA) |

---

## Coordinate conventions

All coordinates are **0-based, half-open** (pysam / BED convention).

| Strand | 5' end (TSS) | 3' end (CPA) |
|---|---|---|
| `+` | `read.reference_start` | `read.reference_end - 1` |
| `-` | `read.reference_end - 1` | `read.reference_start` |

GFF/GTF files are 1-based inclusive. When loading:
```python
start = int(fields[3]) - 1   # GFF 1-based → 0-based
end   = int(fields[4])        # 1-based inclusive → 0-based exclusive (half-open)
```

Chromosome names are normalized to `chrI`, `chrII`, … at genome load time
in `load_genome()`. NCBI-format names (`ref|NC_001133|`) in FASTA headers
are converted transparently; all downstream code uses the canonical format.

---

## Key design principles

**1. Universal A-tract correction, technology-specific stack.**
The A-tract ambiguity correction in `atract_detector.py` applies to every
poly(A)-tailed technology. Technology-specific corrections (AG mispriming
for oligo-dT, poly(A) trimming for direct RNA) stack on top and are
enabled/disabled by flags.

**2. Annotation-aware but annotation-blind scoring.**
Junction BED files are passed to aligners to improve splice site detection,
but the primary alignment scoring does not reward canonical GT-AG motifs.
GT-AG scoring is used only as a tiebreaker between otherwise-equal alignments.
This ensures novel splice sites are not penalized.

**3. Memory-efficient streaming for large datasets.**
`--streaming` mode processes reads in 10,000-read chunks with O(1) peak RAM
relative to BAM size. The `corrected_3ends_index.bed.gz` position count file
(~300× smaller than per-read TSV) further accelerates multi-sample analysis
by making Pass 1 and Pass 2 near-instant.

**4. HPC-aware from the start.**
Thread limits are set before numpy is imported. Scratch staging is built
into `run_command.py`. SLURM profile YAMLs live inside the package.
All SLURM scripts generated by `batch_command.py` include explicit conda
Python paths (no `conda activate`).

**5. Provenance tracking without Snakemake.**
Every output file gets a sidecar `.provenance.json` with SHA-256 hashes
of its inputs. `step_needed()` enables skip-if-unchanged re-runs without
requiring a workflow manager.

**6. Single genome normalization point.**
All FASTA chromosome name translation happens in `load_genome()` and never
again inside the pipeline. Any `CHROM_TO_GENOME` / `GENOME_TO_CHROM` lookup
outside of `load_genome()` is dead code (legacy from pre-2.x).
