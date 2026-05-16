# Split `rectify/core/bam/bam_writer.py`

**Current:** 2,187 lines. Houses every "write a corrected BAM" path
that the pipeline produces, plus several read-manipulation primitives
(clip, softclip, fix-homopolymer, re-align, extend, hardclip).

## Goal

Separate the small per-read primitives (functional, testable in
isolation) from the file-level writers (which orchestrate multiple
primitives over a whole BAM). The primitives form a "BAM read editing
toolkit" that other modules can call independently.

## Proposed split

| New file | Contents | Approx lines |
|---|---|---|
| `bam/bam_writer.py` (kept) | `write_corrected_bam`, `write_softclipped_bam`, `write_dual_bam`, `write_polya_trimmed_bam`, `_load_corrections_from_tsv` (the file-level writers) | ~700 |
| `bam/read_edits.py` (new) | `clip_read_to_corrected_3prime`, `softclip_read_to_corrected_3prime`, `fix_homopolymer_mismatches`, `realign_exon_blocks`, `_normalize_cigar_ops`, `softclip_intronic_tail_5prime`, `reroute_intronic_tail_5prime_via_junction`, `_cigar_ref_end`, `extend_read_5prime_for_junction_rescue`, `extend_read_3prime_for_softclip_rescue`, `_hardclip_trailing_a_run` (per-read CIGAR/seq surgery primitives) | ~1,200 |
| `bam/bedgraph_writers.py` (new) | `write_netseq_assigned_bedgraph`, `write_corrected_reads_bedgraph` | ~190 |

The per-read primitives in `read_edits.py` are individually small
(~100-200 lines each) but together form a coherent toolkit. Several
of them are unit-tested in `tests/test_bam_writer.py`.

## Migration plan

1. Create `bam/read_edits.py` and move all per-read primitives.
   Most are pure functions of `(pysam.AlignedSegment, args...) ->
   pysam.AlignedSegment` so they have no shared state.
2. Create `bam/bedgraph_writers.py` and move the two bedgraph writers.
3. Update the file-level writers (`write_corrected_bam`, etc.) to
   import from `..read_edits`.
4. Update external callers — chiefly
   `commands/correct_command.py`, `commands/restore_polya_command.py`,
   `commands/split_command.py`, `commands/run_command.py`, and
   the chain of tests.

## Test gate

- `tests/test_bam_writer.py` — direct tests of the primitives
- `tests/test_validation_reads.py` — exercises the full writer path
- `tests/test_xr_flag.py`, `tests/test_chimeric_consensus.py`,
  `tests/test_gapmm2_seq_restore.py` — call the writer's edit
  primitives indirectly via consensus selection

## Critical pitfalls

- `realign_exon_blocks` uses `parasail` (optional dependency, loaded
  lazily). Preserve the lazy-import pattern when moving it.
- `fix_homopolymer_mismatches` uses module-level scratch buffers for
  performance — these need to move with it. (Same kind of pattern as
  `bam_processor.py`.)
- `write_corrected_bam` writes BOTH the corrected BAM and the
  trimmed-polyA BAM in some configurations. Don't accidentally split
  the trimmed-polyA case into a separate writer; it's path-dependent.
- The two `_bedgraph` writers consume `corrected_reads.tsv` files and
  emit per-strand bedgraphs. They share file-naming conventions with
  `commands/correct_command.py`; double-check the suffix conventions
  if you change the import path.
