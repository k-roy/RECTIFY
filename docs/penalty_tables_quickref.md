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
| D | AT | 0.44 | 0.17 | 0.03 |
| D | CG | 0.85 | 0.35 | 0.03 |
| X | AT | 1.00 (ref) | 0.64 | 1.24 |
| X | CG | 1.00 (ref) | 1.77 | 1.26 |

Deletion penalty drops sharply with HP length (deletions become nearly
free in long runs). Substitution penalty is non-monotone: it dips below
the HP=1 reference mid-tract (subs more expected there) then rises again
at long runs where substitutions are rare.

Insertions: populated and loaded (`low_count=False`) at short HP, not a
uniform cap. I AT HP=1 = 1.25, falling to ~0.15 by HP=6; only the
sparse long-HP tail (HP≳15) is `low_count=True` / near the 10.0 cap.

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
- `low_count=True` rows (the profiler flags `count_total < 50`, its
  `--min-count-flag` default) should not be used as reliable estimates.
  Separately, the loader (`HpPenaltyTable.from_tsv`) skips rows below its
  own `min_count=100` runtime threshold.
- Tables are S. cerevisiae R10.4.1-specific — regenerate for other
  organisms/chemistries.
- `_CANONICAL_HP_PRIOR = 0.5` (defined in `junction_scoring.py`, applied
  in `junction_refiner.py`) was calibrated so that one HP deletion
  (`del_cost ≈ 0.44–0.85` at HP=1) gives canonical junctions the
  tie-breaking advantage. Works correctly with both default heuristics
  and empirical tables.
