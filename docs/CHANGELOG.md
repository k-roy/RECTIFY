# Changelog

All notable changes to RECTIFY will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
