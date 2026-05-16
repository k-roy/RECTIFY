# Split `rectify/core/commands/cdna_correct_command.py`

**Current:** 2,039 lines. The `rectify correct-cdna` (Stage 1)
implementation: UMI extraction, directional clustering, abPOA
consensus, walkback, pretrim, FASTQ emission.

## Goal

The file is half "command entry point" and half "ONT PCR-cDNA
algorithms." Carve the algorithms into a new `rectify/core/cdna/`
subpackage; leave only the CLI wiring + orchestration in
`commands/cdna_correct_command.py`.

## Proposed split

| New file | Contents | Approx lines |
|---|---|---|
| `commands/cdna_correct_command.py` (kept) | `run()`, `create_correct_cdna_parser()`, thin orchestration calling into the cdna/ modules | ~250 |
| `cdna/read_info.py` (new) | `ReadInfo` dataclass, `detect_full_length_tier`, `detect_other_end_adapter`, `revcomp`, `_find_adapter_anchor_pos`, `_find_anchor_fuzzy`, `extract_read_info` | ~250 |
| `cdna/walkback.py` (new) | `walk_back_anchor_and_tail`, `walk_forward_tss`, `_score_boundary_window`, `_find_boundary_match` (currently thin wrappers around `core/correct/walkback`) | ~250 |
| `cdna/umi.py` (new) | `umi_components`, `umi_components_directional`, `position_components`, `position_components_directional`, `canonical_umi`, `umi_canon_placeholder` | ~280 |
| `cdna/cluster.py` (new) | `cluster_reads`, `pick_representative` | ~80 |
| `cdna/consensus.py` (new) | `pileup_consensus`, `poa_consensus`, `poa_consensus_from_strings`, `poa_consensus_strand_aware`, `pretrim_consensus` | ~280 |
| `cdna/isoforms.py` (new) | `assign_isoforms`, `reconcile_t1_t2_pairs`, `_emit_isoform` | ~190 |
| `cdna/gff.py` (new) | `_parse_gff_gene_name`, `load_rdna_intervals`, `load_gff_genes`, `classify_sense_antisense` | ~140 |
| `cdna/io.py` (new) | `stream_reads`, `write_stage1_fastq` | ~270 |

This is a LOT of files. The alternative is one big `cdna/internals.py`
(~1,800 lines) — choose based on your taste. The above is the more
ambitious version that mirrors how `analyze/` and `splice/` are
already organized.

## Migration plan

1. Create `rectify/core/cdna/` with `__init__.py`.
2. Move modules one at a time, smallest-first. Suggested order:
   `gff.py` → `io.py` → `read_info.py` → `walkback.py` → `umi.py`
   → `cluster.py` → `consensus.py` → `isoforms.py`.
3. For each move: copy the function(s), fix imports (this file's
   relatives become `from ..correct.walkback`, etc.), add re-export
   stubs to `cdna_correct_command.py` so `run()` continues to work.
4. After all moves, delete the re-exports from
   `cdna_correct_command.py` and update `run()` to import from
   `..cdna.X` directly.
5. Verify test_cdna_correct.py, test_cdna_analyze.py,
   test_cdna_chain_canary.py, test_validation_reads_cdna.py.

## Test gate

- `tests/test_cdna_correct.py` (slow: includes 5-min pipeline smoke)
- `tests/test_cdna_analyze.py`
- `tests/test_cdna_chain_canary.py` (slow: full chain canary)
- `tests/test_validation_reads_cdna.py`

Plan ~15 minutes per pipeline test run; budget accordingly.

## Status

**Completed 2026-05-16.** Final layout matches the table above with these
deviations:

- A small `cdna/_constants.py` was added for the chemistry constants
  (`SSP_FWD`, `SSP_RC`, `UMI_LEN`, `BRIDGE_LEN`, `ANCHOR_*`, polyA regex,
  `COMPLEMENT_TABLE`). Multiple target modules need them; putting them in
  any one of them creates a hub-and-spoke import shape that this small file
  avoids.
- `_find_adapter_anchor_pos` lives in `cdna/walkback.py`, not
  `cdna/read_info.py` as the table above suggests. `extract_read_info`
  (read_info.py) calls `walk_back_anchor_and_tail` (walkback.py), and
  `walk_back_anchor_and_tail` calls `_find_adapter_anchor_pos` — putting the
  anchor finder in read_info would create a read_info ↔ walkback cycle.
  `consensus.pretrim_consensus` imports it from walkback alongside the
  walkback wrappers.
- `cdna_correct_command.py` ended at 330 lines (vs ~250 target). The parser
  function `create_correct_cdna_parser` is ~150 lines and stays in the command
  file by design (matches how other commands are laid out).
- Test/script imports were updated to point at the new module locations
  (no re-export shims left in `cdna_correct_command.py`):
  `tests/test_cdna_correct.py`, `tests/test_validation_reads_cdna.py`,
  `scripts/validation_data/characterize_baseline.py`.

Test gate: all four target tests pass (test_cdna_correct 25/25 incl. the
5-min smoke; test_cdna_analyze 3/3; test_cdna_chain_canary 1/1 in 4:21;
test_validation_reads_cdna 16/16).

## Critical pitfalls

- `ReadInfo` is referenced by name in `cdna_analyze_command.py`
  via `from rectify.core.commands.cdna_correct_command import (...)`
  — when you move it to `cdna/read_info.py`, update that import.
- `walk_back_anchor_and_tail` delegates to
  `rectify.core.correct.walkback.walkback_3prime_with_qpos`; the
  delegation pattern stays the same after the move.
- The `_BOUNDARY_PATTERN_FWD` / `_BOUNDARY_PATTERN_REV` module-level
  constants belong with `walk_forward_tss` — they're used only by it.
- `pyabpoa` is an optional dep (loaded lazily). Preserve the
  `try: import pyabpoa except: HAS_PYABPOA = False` pattern when
  splitting `consensus.py` out.
- `intervaltree` and `rapidfuzz` are also imported at module top; they
  belong with the modules that use them (`gff.py` for intervaltree,
  `read_info.py` for rapidfuzz).
