# Empirical Penalty Tables (`--junction-penalty-table`)

Quick reference. For algorithm theory + design rationale see
`docs/EMPIRICAL_HP_PENALTY_SCORING.md`.

Production tables: `common/scripts/nanopore/error_profile_20260422/`
(generated 2026-04-22).

## Format: AT/CG base-class split

`penalty_scores.tsv` groups reference bases into two classes: **AT** and
**CG**. AT runs have ~10–20% higher deletion rates than CG runs at the
same HP length (Nanopore pore-ratcheting asymmetry).

Columns: `op_type`, `base_class`, `hp_length`, `rate_mean`,
`count_total`, `penalty_score`, `low_count`.

Key numbers:

| op_type | base_class | HP=1 | HP=4 | HP=8 |
|---|---|:---:|:---:|:---:|
| D | AT | 0.37 | 0.33 | 0.28 |
| D | CG | 0.58 | 0.42 | 0.37 |
| X | AT | 1.00 (ref) | 1.55 | 10.0 (cap) |
| X | CG | 1.00 (ref) | 2.56 | 10.0 (cap) |

Insertions: all `low_count=True`, treat all as 10.0.

## Usage

```bash
rectify correct \
    --junction-penalty-table .../error_profile_20260422/penalty_scores.tsv \
    --str-penalty-table      .../error_profile_20260422/str_penalty_scores.tsv \
    ...
```

## Regeneration

See `common/scripts/nanopore/PENALTY_TABLE.md` for the full regeneration
command, known caveats (deSALT failures on chunks 2/11/15), and
diagnostic plot instructions.

## Design notes

- Isotonic smoothing ensures deletion penalties are monotone
  non-decreasing with HP length.
- `low_count=True` rows (count < 100) should not be used as reliable
  estimates.
- Tables are S. cerevisiae R10.4.1-specific — regenerate for other
  organisms/chemistries.
- The `_CANONICAL_HP_PRIOR = 0.5` in `junction_refiner.py` was
  calibrated so that one HP deletion (`del_cost ≈ 0.37–0.58` at HP=1)
  gives canonical junctions the tie-breaking advantage. Works correctly
  with both default heuristics and empirical tables.
