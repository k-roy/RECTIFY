# Split `rectify/core/commands/analyze_command.py`

**Current:** 2,093 lines. Drives `rectify analyze`. Already has a
peer `rectify/core/analyze/` subpackage with the actual analysis
modules (clustering, motif_discovery, apa_detection, etc.).

## Goal

Make `analyze_command.py` a thin orchestrator. Move the data-loading
helpers and the manifest-mode runner into `analyze/`.

## Proposed split

| New file | Contents | Approx lines |
|---|---|---|
| `commands/analyze_command.py` (kept) | `run_analyze()` (top-level orchestrator), `create_analyze_parser()` | ~450 |
| `analyze/loaders.py` (new) | `load_corrected_positions`, `_load_large_file_chunked`, `load_annotation`, `load_position_index`, `_parse_gtf` | ~330 |
| `analyze/exclusions.py` (new) | `_detect_exclusion_regions` (this might overlap with `core/exclusion_regions.py` — investigate whether they should be merged) | ~90 |
| `analyze/manifest.py` (new) | `_run_analyze_manifest` (currently 780 lines — by far the biggest function in this file) | ~780 |
| `analyze/bedgraph.py` (new) | `generate_bedgraphs` | ~90 |

This is the cleanest of the giants because the destination subpackage
already exists.

## Migration plan

1. Move `_run_analyze_manifest` first (it's the biggest, and moving
   it gives the biggest size reduction). It's a self-contained
   function — most of its size comes from inline logic that could
   itself be broken up later, but don't try in this PR.
2. Move the loaders (`load_corrected_positions`, `load_annotation`,
   etc.) into `analyze/loaders.py`. Each is independent of the others.
3. Move `_detect_exclusion_regions` into `analyze/exclusions.py`.
   Note: there's a top-level `rectify/core/exclusion_regions.py` —
   check whether they should be consolidated.
4. Move `generate_bedgraphs` into `analyze/bedgraph.py`.
5. Update `run_analyze` to import from the new locations.

## Test gate

- `tests/test_analyze.py` — extensive analyze tests
- `tests/test_apa_detection.py`
- `tests/test_splice_summary.py`
- `tests/test_sample_column_autodetect.py`

## Critical pitfalls

- `_run_analyze_manifest` shares a lot of variables with `run_analyze`
  (config, output dirs, intermediate dataframes). When you move it,
  add explicit parameters for everything `run_analyze` was passing
  via closure.
- `generate_bedgraphs` is called from `run_analyze` near the end; it
  uses several module-level constants that need to migrate with it.
- Some `analyze/` submodules import back from `analyze_command.py`
  (check via `grep -rn analyze_command rectify/core/analyze/`). If
  so, these are circular-import risks — break by moving the shared
  helper into a new `analyze/_common.py`.
