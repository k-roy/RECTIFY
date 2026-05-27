# Homo sapiens empirical penalty tables

Human analog of the yeast (`saccharomyces_cerevisiae`) empirical HP/STR penalty
tables, consumed by `rectify correct` / `run` junction refinement (Module 2H)
and 3'SS rescue scoring. Auto-resolved for `--organism homo_sapiens` (aliases:
`human`, `hg38`, `grch38`) via `rectify/data/__init__.py` `BUNDLED_GENOMES`.

## Files (DRS)

| File | Contents |
| --- | --- |
| `penalty_scores.tsv` | HP-context deletion/substitution/insertion penalties, AT/CG base-class split, isotonic-smoothed |
| `str_penalty_scores.tsv` | STR (di/tri-nucleotide repeat) slippage penalties by (unit, n_copies) |
| `error_rates.tsv`, `error_counts.tsv` | Raw rates + counts the penalties derive from |
| `PROVENANCE.json` | Machine-readable provenance |

`junction_overhang_table.tsv` is intentionally **not** bundled for human — its
calibration needs `rectify correct` output (`corrected_reads.tsv`), which carries
the documented human-`correct` performance caveat; the resolver falls back to the
DRS default. cDNA (`_cdna` + per-UMI-depth `_umi1/2/3plus`) and QuantSeq
(`_qsrev`) tables are pending separate human data (see project notes).

## Derivation (2026-05-27)

- **Sample:** GM12878 (= NA12878 = GIAB HG001), whole-transcriptome **IVT**
  (in-vitro-transcribed, unmodified RNA) — ENA `ERR15839422` (`PRJEB102644`),
  31,086,477 reads.
- **Chemistry:** ONT direct RNA **SQK-RNA004**; dorado 0.8.1, model
  `rna004_130bps_sup@v5.1.0`.
- **Reference:** GRCh38 primary assembly (GENCODE v44).
- **Aligners:** minimap2 + uLTRA + deSALT, strict 3-aligner consensus,
  `--max-intron 100000`. gapmm2 dropped (too slow on human); mapPacBio dropped
  (games HP-ED consensus). ~8.4M reads common to all 3 aligners (chr1–22).
- **Masking:** GIAB HG001 NISTv4.2.1 small-variant VCF (3.89M records, indels
  expanded) excluded + restricted to the high-confidence-regions BED (chr1–22);
  `--exclude-3prime-bp 50` (poly-A boundary; IVT poly-A ~15 bp ≈ trim-equivalent
  for internal-HP rates); `--isolation-flank 10`.
- **Compute:** 45-chunk align array → per-aligner merge → 22-chromosome parallel
  profiler (`empirical_cigar_error_profiler.py` + GIAB masking) → summed
  `error_counts` → `compute_rates` / `derive_penalties` (HP) and
  `_write_str_outputs` (STR).

## Validation

HP penalties match the yeast R10.4.1 reference closely (D AT HP1 0.46 vs 0.44;
HP4 0.18 vs 0.17; HP8 0.032 vs 0.033; D CG HP1 0.76 vs 0.85; X AT/CG HP1 = 1.0
baseline), well-sampled (HP=1 denominator ~2e9 positions, all cells through
HP=10 above the low-count threshold).

## Caveats

- **IVT, not native.** Avoids A-to-I-editing inflation of the AT substitution
  anchor, but does **not** capture modification-induced miscalls present in
  native RNA004 production data. Validate against native RNA004 before driving
  native production. Carries a T7/RT error signature.
- Model-version note: reference is `@v5.1.0`; the Sumner production cohort is
  `@v5.0.0` / `@v5.2.0` (same RNA004 sup-model family).
