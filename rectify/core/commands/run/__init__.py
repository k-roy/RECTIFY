"""
Sub-subpackage for ``rectify run-all`` pipeline runners.

Each module isolates one slice of the orchestrator:

- ``helpers``       — reference-path resolution, per-aligner BAM collection,
                      BAM integrity / MD-tag checks, canonical rectified-BAM path.
- ``stages``        — per-stage runners called by sample-level dispatch
                      (alignment, correction, per-aligner correction,
                      corrected-TSV combination, analysis, junction aggregation).
- ``single_sample`` — single-sample pipeline (``_run_single_sample``) and the
                      shared per-sample worker (``_process_one_sample``).
- ``multi_sample``  — manifest-driven multi-sample pipeline + manifest-mode
                      analysis runner.
- ``chunked_batch`` — generator for scheduler-array shell scripts
                      (``--chunked-alignment``); emits chunk-split/align/merge/
                      correct-analyze submission chains for single-sample or
                      manifest input.
"""
