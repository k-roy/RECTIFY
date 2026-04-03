# Changelog

All notable changes to RECTIFY will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.6.0] - 2026-04-03

### Fixed
- **mapPacBio AssertionError on RNA004/Dorado basecalled reads** (`mpb_split_reads.py`):
  Dorado basecaller using the RNA004 model outputs `U` (uracil) instead of `T` in FASTQ
  sequences. BBMap v39.26 does not handle `U` and throws a Java `AssertionError` in
  `genMatchStringForSite()` on any such read, crashing the entire alignment job.
  Fix: normalize `U→T` in `split_long_reads()` before writing the split FASTQ, so all
  reads reaching mapPacBio are in standard DNA notation regardless of basecaller.

- **Missing `_GenomeDictReference` adapter in `false_junction_filter`** (`false_junction_filter.py`):
  `bam_processor.py` called `_fjf._GenomeDictReference(genome)` which existed in beta
  but was never ported to production, causing a crash during the correction step for
  every sample: `AttributeError: module 'rectify.core.false_junction_filter' has no
  attribute '_GenomeDictReference'`. Ported the thin adapter class (wraps a genome dict
  to present a `pysam.FastaFile`-compatible interface) from beta.

## [2.5.0] - 2026-04-02

### Added
- **uLTRA aligner support** (`multi_aligner.py`, `align_command.py`, `run_command.py`, `cli.py`):
  Annotation-guided collinear chaining aligner optimized for small exons (11–20 nt).
  Handles GFF→GTF conversion automatically (looks for sibling `.gtf`, or converts via
  `gffread`). Decompresses gzip genomes to a temp file (uLTRA cannot read gzip).
  Enabled via `--junction-aligners uLTRA` (opt-in, requires `--annotation`).

- **deSALT aligner support** (`multi_aligner.py`, `align_command.py`, `run_command.py`, `cli.py`):
  De Bruijn graph mapper. Handles deSALT's known output-duplication bug by deduplicating
  on (read_name, flag, chrom, pos, cigar) via `_dedup_desalt_bam()`. Strips
  `LD_LIBRARY_PATH` from deSALT's environment (prevents library conflicts). Resolves
  the deSALT index directory (`desalt_index/`) adjacent to the genome FASTA.
  Enabled via `--junction-aligners deSALT` (opt-in).

- **5-aligner consensus pipeline**: `--junction-aligners uLTRA deSALT` adds both
  junction-mode aligners to the 3-aligner default (minimap2 + mapPacBio + gapmm2),
  producing a 5-aligner consensus BAM.

- **`run-all` CLI exposes junction aligner flags**: `--junction-aligners`, `--ultra-path`,
  `--desalt-path`, `--chimeric-consensus` added to the `run-all` subcommand.

### Changed
- **Two-phase scheduler updated for deSALT fork-safety**: deSALT crashes with
  "double free or corruption" when forked from a multithreaded Python process.
  The scheduler now runs deSALT sequentially after the parallel pool completes
  (`parallel_batch = [a for a in remaining if a != 'deSALT']`).

- **`_run_alignment()` in `run_command.py`** now accepts and threads through
  `junction_aligners`, `chimeric_consensus`, `ultra_path`, `desalt_path`.

## [2.4.0] - 2026-04-01

### Added
- **Chimeric consensus selection** (`chimeric_consensus.py`):
  Divides a read into segments at "sync points" where all aligners agree, then
  independently selects the best aligner per segment and assembles a chimeric CIGAR.
  Wired into the pipeline behind `--chimeric-consensus` flag (default off, experimental —
  requires further validation before enabling by default).

- **Test suite ported from beta**: Unit and integration tests for consensus selection,
  chimeric consensus, junction validation, false junction filtering, and aligner wrappers.

### Changed
- **Visualization module updated**: Additional plot types and improved styling.

## [2.3.0] - 2026-03-30

### Added
- **Module 2F: 3'SS truncation rescue** (`bam_processor.py`):
  Corrects 5' positions for reads truncated or soft-clipped at the exon 2 / 3' splice
  site boundary. Runs post-consensus on reads with real (non-poly-A-artifact) junctions.

- **Memory-efficient streaming pipeline for multi-sample analysis** (`analyze_command.py`):
  Two-pass streaming over per-sample corrected TSVs using a manifest; never loads more
  than one sample into RAM. Position index (`corrected_3ends_index.bed.gz`) support for
  ~300× faster Pass 1 and Pass 2. Peak RAM is O(clusters × samples) regardless of depth.

- **Pre-consensus 3' A-tract scoring**: A-tract ambiguity is now scored before consensus
  selection so it can inform aligner selection, not just final position correction.

- **Position index writing**: `rectify correct` now writes `corrected_3ends_index.bed.gz`
  alongside `corrected_3ends.tsv` for use by manifest-mode analysis.

### Fixed
- **Spike-in filter false positives** (`correct_command.py`): Fixed over-aggressive
  filtering that was removing valid reads near spike-in sequences.

- **Bug 14: Silent exception swallow in `false_junction_filter.py`**: Exceptions in
  junction motif analysis were silently swallowed (bare `except: pass`), hiding real
  errors. Now logs at DEBUG level with the exception message.

- **Bug 16: CIGAR/seq mismatch for hard-clipped donors** (`consensus.py`):
  `_restore_sequence_from_aligner_reads()` copied donor sequence without checking length
  against CIGAR-implied query length. Hard-clipped donors (deSALT) caused
  `samtools sort` crash: "CIGAR and query sequence lengths differ". Fixed with
  `_cigar_query_length()` helper; donors whose sequence length doesn't match are skipped.

- **`ref_end is None` guards in `bam_processor.py`**: Unmapped reads (no `reference_end`)
  caused `TypeError: unsupported operand type(s) for -: 'NoneType' and 'int'` in
  `get_read_3prime_position()` and `get_read_5prime_position()`. Now returns `None`
  for unmapped reads.

- **`itertuples()` OOM in `junction_validator.py`**: `filter_cross_sample_junctions()`
  used `itertuples()` on a large aggregated DataFrame, causing OOM on large datasets.
  Replaced with numpy array iteration.

### Changed
- **Two-phase aligner scheduler**: mapPacBio (~10× slower) runs alone with all threads
  (phase 1); remaining aligners split threads equally in parallel (phase 2). Replaces
  the old `_THREAD_WEIGHTS` dict approach.

- **Per-aligner BAMs name-sorted directly**: Skips intermediate coordinate sort and
  index, reducing disk I/O during the alignment phase.

## [2.2.0] - 2026-03-28

### Added
- **`run-all` command** (`run_command.py`): Complete end-to-end pipeline
  (align → correct → analyze) for single samples and multi-sample manifests.
  Single-sample: steps 0–3 (alignment, correction, analysis, junction aggregation).
  Multi-sample: parallel per-sample correction + combined DESeq2/GO/motif analysis.

- **Scratch staging** (`slurm.py`): SLURM array jobs automatically stage I/O through
  `$SCRATCH` (75 GB/s) to avoid Oak NFS contention. Final outputs sync'd back to Oak
  via `rsync`. `make_job_scratch_dir()` / `sync_to_oak()` utilities.

- **Infrastructure ported from beta**:
  - `junction_validator.py`: Cross-sample junction validation with configurable thresholds
  - `slurm.py`: CPU detection, thread limits, scratch staging utilities
  - Sherlock SLURM profile (`sherlock_larsms.yaml`) with `use_scratch: true` and
    `streaming: true` defaults

## [2.1.0] - 2026-03-16

### Added
- **train-polya command**: Train poly(A) tail models from control sites with 0 downstream A's
  - Learns A-richness distribution, length distribution, position profiles
  - Outputs JSON model file for use with `correct` command
  - Supports minimum read count thresholds and per-site statistics

- **validate command**: Comprehensive validation against ground truth
  - Multiple ground truth sources: NET-seq BigWig, gene annotations (GTF/GFF), known positions TSV
  - All sources optional but at least one required
  - Calculates accuracy at multiple tolerances (exact, 1bp, 2bp, 5bp)
  - Per-correction-type breakdown and confidence calibration metrics

- **PolyAModel class**: Dataclass-based poly(A) tail model with JSON serialization
  - Training statistics tracking
  - Position-dependent A-frequency profiles
  - Score sequences against learned parameters

- **LRU cache in NetseqLoader**: Bounded cache (default 10,000 entries) prevents memory issues

- **GitHub Actions CI/CD**: Automated testing across Python 3.8-3.12, linting, and package building

### Changed
- **A-tract ambiguity detection**: Now correctly counts downstream A's in window (not just tract length)
  - Added `_get_sequence()` helper for chromosome name conversion
  - Added `_count_downstream_as()` for proper window-based counting
  - Returns both `downstream_a_count` and `expected_shift` fields

- **Performance optimizations**:
  - `find_polya_boundary()`: O(n) sliding window instead of recalculating per position
  - `complement()`: Uses `str.maketrans`/`translate` instead of per-character dict lookup
  - `generate_summary_report()`: Single-pass accumulation instead of 5 iterations
  - `find_coverage_gaps()`: Avoid O(n) `.index()` call in chromosome lookup

- **Consolidated `find_atract_boundaries`**: Removed duplicate implementation from `genome.py`

### Fixed
- Version mismatch between config.py and pyproject.toml
- Test expectations aligned with actual algorithm behavior:
  - Insertions do NOT change reference coordinates (counted for QC only)
  - Soft-clipped poly(A) tails don't shift genomic position

### Dependencies
- Added `tqdm>=4.60.0` for progress bars

## [2.0.0] - 2026-03-09

### Added
- Initial unified RECTIFY framework
- A-tract ambiguity detection module
- Poly(A) tail trimming with adapter detection
- Indel artifact correction
- NET-seq refinement for position validation
- AG mispriming screening
- SLURM-aware thread management
- Parallel processing with coverage-gap-based region splitting

### Core Modules
- `atract_detector.py`: Universal A-tract ambiguity detection
- `polya_trimmer.py`: Poly(A) tail detection and scoring
- `indel_corrector.py`: Small indel artifact correction
- `netseq_refiner.py`: NET-seq-based position refinement
- `ag_mispriming.py`: AG-rich region flagging
- `bam_processor.py`: Parallel BAM processing pipeline

### Commands
- `rectify correct`: Main correction pipeline
- `rectify train-polya`: (stub) Train poly(A) model
- `rectify validate`: (stub) Validate corrections
