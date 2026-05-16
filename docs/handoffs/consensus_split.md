# Split `rectify/core/consensus/consensus.py`

**Current:** 2,033 lines. Houses the multi-aligner consensus selection
logic — `select_best_alignment`, `extract_alignment_info`,
`score_alignment`, plus per-read scoring helpers and the streaming
multi-BAM iterator.

## Goal

Group the file into three layers: scoring primitives, alignment-info
extraction, and the streaming selection orchestrator. The primitives
deserve their own module so they can be unit-tested without setting
up a streaming pipeline.

## Proposed split

| New file | Contents | Approx lines |
|---|---|---|
| `consensus/consensus.py` (kept) | `run_consensus_selection`, `_process_and_write_batch`, the streaming/name-grouping iterator helpers (`_ensure_name_sorted`, `_read_id_hash`, `_filtered_read_iterator`, `_natural_sort_key`, `_normalize_bam_read_name`, `_iter_name_grouped_bams`, `_cigar_query_length`, `_restore_sequence_from_aligner_reads`) | ~700 |
| `consensus/scoring.py` (new) | `score_alignment`, `_cigar_terminal_errors`, `_get_effective_5prime_clip`, `_get_effective_3prime_clip`, `_count_junction_proximity_errors`, `_is_homopolymer_position`, `check_canonical_splice_sites`, `detect_false_3prime_junction`, `_rescue_5prime_softclip` (per-alignment scoring primitives + softclip rescue) | ~900 |
| `consensus/extract.py` (new) | `AlignmentInfo` dataclass, `ConsensusResult` dataclass, `extract_junctions_from_cigar`, `get_softclip_lengths`, `extract_alignment_info` (CIGAR → AlignmentInfo conversion) | ~300 |
| `consensus/select.py` (new) | `select_best_alignment` (the actual selection logic that compares scored AlignmentInfos) | ~120 |

## Migration plan

1. Move the dataclasses + `extract_alignment_info` to `extract.py`
   first — they're the type system everything else uses.
2. Move the scoring primitives to `scoring.py`. They depend only on
   `AlignmentInfo` (now in `extract.py`) and stdlib.
3. Move `select_best_alignment` to `select.py`. It depends on both
   `extract.py` and `scoring.py`.
4. Leave streaming/orchestration in `consensus.py`. Update its
   imports from `from .scoring import ...`, etc.
5. Add re-export stubs in `consensus.py` for `AlignmentInfo` and
   `extract_alignment_info` since many external callers import those
   from `consensus.consensus`. Or — preferred — update all callers
   to import from `consensus.extract` directly.

## Test gate

- `tests/test_consensus_selection.py` — tests `select_best_alignment` directly
- `tests/test_chimeric_consensus.py` — tests `select_best_chimeric` (in
  `chimeric_consensus.py` — separate file, but consumes `AlignmentInfo`)
- `tests/test_no_duplicate_primaries.py`
- `tests/test_xr_flag.py`
- `tests/test_gapmm2_seq_restore.py`
- `tests/test_splice_junction.py`

## Critical pitfalls

- `tests/test_gapmm2_seq_restore.py` patches the logger
  `"rectify.core.consensus.consensus"` (the module path). If you
  move `_restore_sequence_from_aligner_reads` out, the logger name
  changes — update the test in lockstep, OR keep that function in
  `consensus.py`.
- `extract_alignment_info` is called once per read per aligner —
  performance-sensitive path. Don't add extra indirection in its
  internals.
- `_rescue_5prime_softclip` writes to the read it's given (mutates
  the `pysam.AlignedSegment` in place). When you move it, document
  this side effect at the top of the new module.
- `chimeric_consensus.py` imports from `consensus.py`; verify no
  circular import when you split. If `consensus.consensus` ends up
  importing from `consensus.scoring` which ends up importing from
  `consensus.chimeric_consensus`, you'll have a problem.
