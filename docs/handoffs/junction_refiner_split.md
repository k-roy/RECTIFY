# Split `rectify/core/splice/junction_refiner.py`

**Current:** 2,136 lines. Implements junction refinement: for every
N-op in every consensus read, test all candidate junctions within
a search radius and replace the imprecise N-op boundary with the
best sequence-supported junction.

## Goal

Separate the scoring/penalty machinery (the `HpPenaltyTable` and edit
distance algorithms) from the refinement orchestration (which iterates
N-ops and applies replacements). The scoring layer is computationally
intensive and well-defined enough to live in its own module.

## Proposed split

| New file | Contents | Approx lines |
|---|---|---|
| `splice/junction_refiner.py` (kept) | `refine_read_junctions`, `refine_bam_junctions`, `_iter_n_ops`, `_has_boundary_error`, `_apply_junction_replacement`, `_apply_replacements_to_read`, `_run_sequential`, `_run_parallel`, `_refine_read_batch`, `_warmup_numba_dp`, `_sort_and_index` (the refinement orchestrator) | ~1,000 |
| `splice/junction_scoring.py` (new) | `_score_hp_anchored`, `_score_junction`, `_score_junction_fallback`, `_3ss_tier_from_rna_trinucleotide`, `_canonical_tier`, `collect_junctions_from_bam`, `build_junction_pool`, `_build_junction_index`, `_candidates_near` (per-junction scoring + candidate pool) | ~600 |
| `splice/hp_penalty.py` (new) | `HpPenaltyTable`, `_hp_run_length`, `_str_repeat_info`, `_precompute_del_costs`, `_hp_edit_distance` (homopolymer-aware edit distance) | ~440 |

## Migration plan

1. `hp_penalty.py` first — pure algorithmic module, no external
   state. `HpPenaltyTable` is the central dataclass; everything else
   here is helper functions for it.
2. `junction_scoring.py` next. It imports `HpPenaltyTable` and uses
   it to score candidate junctions.
3. The refiner stays in `junction_refiner.py` and imports from the
   new modules.
4. Update tests in `tests/test_junction_refiner.py` and
   `tests/test_junction_validator.py` to point at the new locations
   where they patch internals.

## Test gate

- `tests/test_junction_refiner.py` — 41 tests, heavy
- `tests/test_junction_validator.py`
- `tests/test_terminal_exon_refiner.py`
- `tests/test_splice_junction.py`
- `tests/test_splice_summary.py`

## Critical pitfalls

- `_hp_edit_distance` is compiled by numba on first call; the
  `_warmup_numba_dp` function pre-triggers the compilation. If you
  move them apart, ensure `_warmup_numba_dp` still triggers the
  compile on the right function (it imports and calls it).
- `HpPenaltyTable` is loaded lazily from a JSON file (or built
  empirically); the load path is referenced from `commands/correct_command.py`.
  Check that path after the move.
- `build_junction_pool` is also called from
  `commands/run_command.py` and `commands/aggregate_command.py`
  for the cross-sample junction aggregation. Verify imports there.
- `_run_parallel` uses `concurrent.futures.ProcessPoolExecutor` and
  pickles `_refine_read_batch` — verify the module path stays
  importable from worker processes after the split.
