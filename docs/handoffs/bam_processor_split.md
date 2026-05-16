# Split `rectify/core/bam/bam_processor.py`

**Current:** 2,308 lines. Imported by `commands/correct_command.py`,
`commands/cdna_analyze_command.py`, `commands/run_command.py`, tests.

## Goal

Reduce the file to the per-read correction core. Move the surrounding
scaffolding (file orchestration, parallel/streaming dispatch, output
writing, region planning) into peer modules under `bam/`.

## Proposed split

| New file | Contents | Approx lines |
|---|---|---|
| `bam/bam_processor.py` (kept) | `correct_read_3prime` (the 665-line workhorse), `get_read_3prime_position`, `get_read_5prime_position`, `_build_pool_chrom_index`, `_load_netseq`, `process_bam_file` | ~900 |
| `bam/regions.py` (new) | `find_coverage_gaps`, `get_processing_regions` | ~140 |
| `bam/parallel.py` (new) | `_process_region_worker`, `process_bam_file_parallel`, `process_bam_streaming`, `process_bam_streaming_parallel`, `_rebuild_pos_counts_from_partial`, `_write_results_chunk` | ~900 |
| `bam/output.py` (new) | `write_output_tsv`, `generate_summary_report`, `generate_summary_from_stats` | ~270 |
| `bam/variant_scan.py` (new) | `run_variant_aware_scan` | ~85 |

The 665-line `correct_read_3prime` is the load-bearing per-read
function — do NOT split it internally without a separate review.

## Migration plan

1. Create empty `bam/regions.py`, `bam/parallel.py`, `bam/output.py`,
   `bam/variant_scan.py` (with module docstrings).
2. For each new file: copy the relevant function(s) from
   `bam_processor.py` and resolve their imports (most use `..polya`,
   `..correct`, `..splice`, `..netseq`, `..utils` already from the prior
   reorg — verify these still work).
3. In `bam_processor.py`, replace each moved function with `from .X
   import Y` re-exports so existing callers keep working. Or — preferred
   — update all callers to import from the new location directly.
4. Delete the moved bodies from `bam_processor.py`.
5. Run test gate (below).
6. Drop the re-exports in `bam_processor.py` once you've confirmed
   nothing imports through it.

## Test gate

- `tests/test_bam_writer.py` — exercises BAM I/O round-trip
- `tests/test_correct_command_drs.py` — drives the correct command (subprocess)
- `tests/test_run_command_wiring.py` — pipeline wiring
- `tests/test_validation_reads.py` — full DRS validation suite
- `tests/test_parallel_processing.py` — parallel correction
- Full broad sweep at the end

## Critical pitfalls

- `process_bam_streaming_parallel` shares pre-allocated scratch arrays
  with `correct_read_3prime` via module-level state. If you move the
  parallel functions to `bam/parallel.py`, those scratch arrays stay
  imported from `bam_processor.py` — multiprocessing forks copy them
  per worker, so the import path matters.
- `_load_netseq` constructs a `NetseqLoader` from
  `..netseq.netseq_refiner`; some workers re-import it. Verify the
  refactored module path is referenced by string anywhere
  (e.g. dynamic imports, multiprocessing pickling).
- The DRS-specific guarded walkback call at `bam_processor.py:679`
  (`walkback_drs_full(read, _chrom_seq)`) is the production wiring
  for `correct_read_3prime`. Don't touch it; it's the entry point
  the other agent's Gap 2 work hooks into.
