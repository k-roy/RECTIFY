# RECTIFY Nanopore Error Profiles + Penalty Tables

## Purpose
Empirical error profiles and penalty score tables derived from BY4742 DRS alignments,
used by RECTIFY's local re-alignment step to score indel corrections. Contains two
independently derived profiles (April 2026 and May 2026) plus strand-stratified rates.
These are internal lab resources, not public datasets.

## Contents
```
empirical_cigar_error_profiler.py     Profiler script (run on DRS BAM to generate profiles)
combined_100k_counts.tsv              Aggregate error counts (100k reads)
penalty_scores.tsv                    Canonical penalty table (symmetric, used in production)
hp_ambiguity_diagnostic.py            Diagnostic: homopolymer ambiguity analysis
plot_error_profiles.py                Visualization script
run_profiler_20260422.sh / _20260501.sh  SGE/SLURM submission scripts

error_profile_20260422/               April 2026 profile (BY4742 DRS wt_rep1)
  penalty_scores.tsv                  Per-base indel penalties
  str_penalty_scores.tsv              STR-aware penalties (homopolymer runs)
  error_rates.tsv / error_counts.tsv  Raw rates and counts
  union_*                             Union-read rates (inclusive counting)
  diagnostics/                        12 diagnostic PNGs
  plots/                              Summary figures

error_profile_strand_20260501/        May 2026 strand-stratified profile
  penalty_scores.tsv                  Strand-aware penalties
  strand_error_rates_by_base.tsv      Per-strand, per-base error rates
  plots/                              Summary figures

overhang_table_20260502/
  junction_overhang_table.tsv         5' overhang distribution table (splice junction QC)
```

## Size
~2.8 MB

## Source
Internal — derived from: `/u/project/guillom/shared/processed/by4742_polya_drs_2025/`
(aligned BAMs on H2). The profiler script is self-contained and can be re-run on any
DRS BAM to regenerate or update error profiles.

## Retrieval / Regeneration
```bash
# Re-run profiler on a DRS BAM (produces all TSVs):
python empirical_cigar_error_profiler.py \
  --bam /u/project/guillom/shared/processed/by4742_polya_drs_2025/<sample>.bam \
  --genome /u/project/guillom/kevinroy/common/references/S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529.fsa \
  --out error_profile_$(date +%Y%m%d)/

# Use qsub wrapper:
qsub run_profiler_20260422.sh  # (update paths as needed)
```

## Migration
2026-05-13 from Sherlock: `/oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/`
