# Lazy Corrected Consensus + Homopolymer Realign Speedup Plan

Date: 2026-05-21
Status: implementation checkpoint in working tree; synced to Sherlock/H2

## Goal

Preserve TRT-relevant homopolymer correction and still produce a final
corrected consensus BAM, without requiring full per-aligner corrected BAMs for
every aligner/chunk.

The current path materializes corrected BAMs for all aligners, then reads them
back for two downstream purposes:

1. `merge_corrected_tsvs()` computes HP-aware edit distance from each
   aligner's final corrected CIGAR and uses it for winner selection.
2. The chunk script builds `corrected_consensus.bam` by copying winning reads
   from the per-aligner corrected BAMs.

That is functionally clean but expensive. For ysh1 chunk 000, minimap2 finished
corrected-BAM writing in 111.6s, while mapPacBio/gapmm2/uLTRA were still in the
writer after >12.6 min extra each despite TSVs already being complete.

## Key Insight

The biology requires corrected-CIGAR scoring and a final corrected consensus
BAM. It does not require writing a full corrected BAM for every aligner first.

We can stream raw per-aligner BAMs, apply the same correction/CIGAR surgery in
memory, and either:

- compute HP-aware edit distance for winner selection, or
- write the corrected record only if that aligner won for that read.

This turns "write N full corrected BAMs, then select" into "score N aligners,
then write one final consensus BAM".

## Critical Caveat: Read Identity for Legacy Runs

Lazy consensus depends on joining three objects by read identity:

1. per-aligner corrected TSV rows,
2. raw per-aligner BAM records,
3. optional trim/read metadata sidecars.

Newer RECTIFY runs solve this with RN-tagged chunk FASTQs,
`*.read_num_sidecar.parquet`, and stage provenance sidecars. The active CPA DRS
and H2 mex67aa runs predate that machinery. Several aligners also mutate QNAMEs
or strip FASTQ comments, so a lazy raw-BAM scorer must not assume raw
`read.query_name` is always the canonical key.

For legacy runs, add a Stage 0 identity bridge before lazy consensus:

1. **Post-hoc sidecar reconstruction from existing chunk FASTQs**
   - no re-alignment
   - reconstruct `read_num`, `original_qname`, `fastq_comment`, `chunk_id`,
     `seq_md5`, and `qual_md5`
   - old round-robin chunking maps local read index `j` in chunk `k` to:
     `read_num = j * n_chunks + k`
   - write an explicit `POSTHOC_PROVENANCE.json` describing source FASTQs,
     manifest, read counts, hash policy, and reconstruction formula

2. **Use sidecar-assisted QNAME fallback first**
   - for CPA DRS UUID reads, normalized QNAME lookup should be sufficient in
     most cases because the original read IDs are unique
   - this improves provenance and metadata restoration without touching BAMs

3. **Inject `RN:i` into existing BAMs only if needed**
   - stream-rewrite BAMs using `qname -> read_num` from the reconstructed
     sidecar
   - validate with read counts, `samtools quickcheck`, indexability, and RN tag
     coverage
   - higher risk because it touches hundreds of BAMs, but still avoids
     re-alignment

4. **Clean re-run remains best provenance**
   - current `rectify split` can emit RN-tagged chunks and sidecars directly
   - cost is redoing alignment, so use only if the post-hoc bridge fails

Lazy consensus implementation uses a read-identity resolver with this
precedence (implemented in `_lookup_read_correction`,
`rectify/core/consensus/corrected_consensus.py:201–210`):

1. **Normalized QNAME match** — raw `read.query_name` is looked up first
   with an exact match against the corrections dict, then via aligner-suffix
   stripping (`_normalize_bam_read_name`, `rectify/core/consensus/consensus.py`)
   with collision exclusion (`_normalized_correction_lookup`,
   `corrected_consensus.py:188–198`). This handles mapPacBio and dorado_source
   suffix variants without a sidecar.
2. **Fail loudly** — when `strict_lazy_identity` is set and a read cannot be
   matched, the resolver raises rather than silently dropping the read.

Note: `RN:i`-tag lookup and sidecar-backed QNAME resolution exist in
`rectify/core/consensus/consensus.py` for the *raw* consensus path, but are
not needed in the *corrected* consensus path. Both the TSV `read_id` column
and the BAM `query_name` are derived from the same canonical QNAME; the only
divergence is aligner-appended suffixes, which `_normalize_bam_read_name`
already strips. Adding sidecar lookup here would add complexity without
fixing any real mismatch.

> **Future-work note:** if RN-tag-first lookup becomes desirable for
> performance reasons (e.g. sidecar-indexed O(1) lookup on very large chunks),
> the wiring point is `_lookup_read_correction` in `corrected_consensus.py`.

## Current Coupling Points

- `rectify/core/consensus/corrected_consensus.py`
  - `merge_corrected_tsvs(..., per_aligner_corrected_bams=...)`
  - if corrected BAMs are supplied, `_read_hp_edit_distances()` scans them and
    joins `hp_edit_distance` / `aligned_bases` into the per-aligner TSV table.
  - if not supplied, winner selection falls back to legacy TSV-sort keys.

- `rectify/core/commands/split_command.py`
  - generated chunk merge script discovers per-aligner `corrected.bam`.
  - passes those BAMs into `merge_corrected_tsvs()`.
  - builds `corrected_consensus.bam` by copying winning primary reads from the
    per-aligner corrected BAMs.
  - generated per-aligner chunk correction currently passes
    `--write-corrected-bam "$SCRATCH_WORK/corrected.bam"` by default.

- `rectify/core/commands/run/stages.py`
  - `_run_correction()` currently constructs a corrected-BAM destination and
    includes it in `correct_args`, even though the `rectify correct` CLI itself
    treats corrected BAM writing as optional. This is a non-chunked run-path
    coupling that must be updated separately from generated split scripts.

- `rectify/core/commands/restore_polya_command.py`
  - uses merged TSV `winning_aligner` plus raw aligner BAMs to restore trimmed
    poly(A) sequence. This already proves a raw-BAM + winning-aligner pattern
    is workable for post-merge BAM generation.

- `rectify/core/bam/bam_writer.py` and
  `rectify/core/bam/bam_writer_parallel.py`
  - the canonical corrected-read edit sequence is duplicated between the
    sequential writer and the parallel writer's private
    `_apply_corrections_to_read()` helper. Lazy scoring/final-BAM writing
    needs this to become a shared helper so all modes mutate reads identically.

## Proposed Architecture

### Stage 0: Post-Hoc Legacy Identity Bridge

Add a utility/command or documented script to reconstruct read sidecars from
existing chunk FASTQs and `chunks_manifest.json` for Sherlock CPA and H2 mex67aa.

Minimum output:

- `<sample>.read_num_sidecar.parquet`
- `<sample>.POSTHOC_PROVENANCE.json`
- per-chunk read-count/hash validation report

Acceptance:

- reconstructed global read count equals the sum of chunk FASTQ records
- `read_num` is unique and contiguous under the old round-robin formula
- `original_qname` is unique for CPA DRS UUID reads
- sampled `seq_md5`/`qual_md5` values match source chunk FASTQs
- sidecar can resolve all corrected TSV `read_id`s by normalized QNAME for the
  four stable aligners on a ysh1 chunk

### Stage 1: Target `realign_exon_blocks()`

Problem: current `realign_exon_blocks()` uses a read-level trigger: once a read
has a homopolymer `X`, it walks all eligible short exon blocks and realigns
them. For ysh1 chunk 000:

| aligner | eligible short blocks | blocks with X | target fraction |
| --- | ---: | ---: | ---: |
| mapPacBio | 184,757 | 106,991 | 57.9% |
| gapmm2 | 193,167 | 100,238 | 51.9% |
| uLTRA | 185,552 | 104,543 | 56.3% |

Change `realign_exon_blocks()` so it records which exon blocks contain
homopolymer-position `X` ops during the pre-check, then realigns only those
blocks. This should preserve intended behavior: the global DP is only needed
where a homopolymer mismatch exists inside that exon block.

Acceptance:

- Existing BAM writer tests still pass.
- Add a focused test with two short exon blocks where only one contains a
  homopolymer `X`; assert only that block changes.
- Add a test where a neighboring clean short exon block remains unchanged.
- Add debug/profile counters for reads inspected, blocks eligible,
  blocks realigned, and blocks skipped.

### Stage 2: Add an In-Memory Corrected-Record Function

Extract the shared per-read mutation sequence from
`bam_writer_parallel._apply_corrections_to_read()` / `bam_writer.write_corrected_bam()`
into a reusable helper, for example:

```python
apply_corrected_bam_edits(read, correction, genome, *, hard_clip=True) -> bool
```

The helper should:

- decode `=` SEQ
- apply 5' reanchor when needed
- apply targeted `realign_exon_blocks()`
- apply 5' rescue CIGAR surgery
- apply 3' softclip/overcall rescue
- apply final 3' hard-clip
- set `cp:i`

It should not write the read.

Acceptance:

- Single-threaded and parallel corrected-BAM writers call the same helper.
- Existing equivalence tests between sequential and parallel BAM writers still
  pass.
- New unit test: helper-mutated record is CIGAR/tag-equivalent to the record
  produced by `write_corrected_bam()` for the same input/TSV row.

### Stage 3: Score HP Edit Distance Without Corrected BAMs

Add a new scoring path to `merge_corrected_tsvs()`:

```python
per_aligner_raw_bams: Optional[Dict[str, str]]
score_from_raw_bams: bool = False
```

When raw BAMs are provided, compute HP-aware edit distances by:

1. Loading each aligner's corrected TSV/manifest into correction rows.
2. Streaming that aligner's raw BAM.
3. For each primary read with a correction row, mutate a copy in memory with
   `apply_corrected_bam_edits()`.
4. Compute `_cigar_hp_edit_distance()` and `_cigar_aligned_bases()` from the
   mutated read.
5. Join scores into `merge_corrected_tsvs()` exactly like the current
   corrected-BAM path.

Missing raw BAM policy: in lazy HP mode, missing raw BAMs should be loud and
deterministic. Prefer failing the merge for production unless an explicit
fallback flag is supplied. Silent `inf` scoring is acceptable for backward
compatibility with the existing corrected-BAM path, but it is not safe as the
default for a new production mode because it can change winners without making
the run obviously incomplete.

The output `aligner_summary.tsv` should be identical, or near-identical where
sort/index order previously created incidental differences, to the full
per-aligner corrected-BAM path.

Acceptance:

- Fixture test: full corrected-BAM scoring and raw-BAM lazy scoring choose the
  same `winning_aligner` for every read.
- The summary TSV includes `hp_edit_distance` and `aligned_bases` in both
  modes.
- Missing raw BAMs fall back with a loud warning rather than silently switching
  scoring modes.

### Stage 4: Build Only the Final Corrected Consensus BAM

After `merge_corrected_tsvs()` selects `winning_aligner`, build the final
`corrected_consensus.bam` directly from raw aligner BAMs:

1. Load merged TSV as `read_id -> winning_aligner`.
2. Load per-aligner correction rows for winning reads.
3. Stream each raw BAM once.
4. If read is primary and that aligner won, mutate read in memory using
   `apply_corrected_bam_edits()`.
5. Write to a single unsorted consensus BAM, then sort/index once.

This replaces the current "copy from full per-aligner corrected BAMs" step in
generated chunk scripts.

Acceptance:

- On a small fixture, lazy final consensus BAM is query-name-equivalent and
  CIGAR/tag-equivalent to the old full-corrected-BAM pipeline.
- Final BAM contains one primary record per winning read.
- Final BAM has `cp:i`, corrected 3' hard-clips, 5' rescue edits, and targeted
  homopolymer CIGAR surgery.
- Sort/index/calmd behavior remains unchanged.

## Expected Savings

For ysh1 chunk 000, corrected TSV generation was ~90s/aligner in the no-2F/no-2H
baseline. Corrected-BAM writing added:

- minimap2: +111.6s
- mapPacBio/gapmm2/uLTRA: >12.6 min each, censored by cancellation

The lazy path should avoid writing/sorting/indexing losing aligner corrected
BAMs entirely. It will still pay in-memory corrected-CIGAR scoring for all
aligners and corrected-record writing for winners, so it is not "free"; the
savings come from eliminating unnecessary BAM I/O and avoiding repeated
materialization of non-winning records.

Targeted `realign_exon_blocks()` should independently reduce realign work by
roughly 35-48% in mapPacBio/gapmm2/uLTRA on ysh1 chunk 000 based on block counts.

## Risks

- Corrected-BAM scoring and lazy scoring must use exactly the same correction
  sequence, or winner selection may drift.
- Pysam record mutation is in-place. The lazy scorer should mutate a copied
  record or ensure records are not reused after scoring.
- `=` SEQ decoding must happen before any sequence-vs-reference scoring, just
  as the writer does today.
- Secondary/supplementary/unmapped handling must mirror the corrected-BAM
  writer and current consensus-builder assumptions.
- `winning_aligner` matching must use the same QNAME normalization as current
  corrected-BAM scoring.
- Read ID normalization is a real hazard: TSV loading normalizes `read_id`,
  while BAM writers historically index by raw `read.query_name`. Lazy code must
  normalize BAM names consistently enough to survive Dorado comments,
  mapPacBio suffixes, BBMap suffixes, and the RN-sidecar transition.
- NET-seq/multi-row reads currently collapse to the first correction row in BAM
  writers. Lazy scoring should preserve that contract unless a separate NET-seq
  redesign is intended.
- deSALT remains unsafe until the separate htslib/heap-corruption issue is
  isolated. The lazy path should not silently re-enable deSALT production.

## Implementation Order

1. Implement targeted `realign_exon_blocks()` and tests.
2. Extract shared `apply_corrected_bam_edits()` helper; wire writers through it.
3. Add lazy HP-score-from-raw-BAM mode to `merge_corrected_tsvs()`.
4. Add lazy final consensus BAM builder.
5. Add post-hoc sidecar reconstruction for legacy Sherlock CPA / H2 mex67aa.
6. Update generated chunk scripts and non-chunked run stages to prefer:
   - per-aligner corrected TSVs
   - raw aligner BAMs
   - sidecar/RN-aware identity resolution
   - lazy HP scoring
   - one final corrected consensus BAM
7. Keep `--write-corrected-bam` available for debug/IGV per-aligner outputs.

## Production Policy If Implemented

For CPA DRS production:

- Run all stable aligners, likely excluding/quarantining deSALT until its BAM
  corruption is solved.
- Keep Plan C capped Module 2H.
- Emit per-aligner corrected TSVs and indexes.
- Do not emit per-aligner corrected BAMs by default.
- Use lazy HP scoring for multi-aligner winner selection.
- Emit one final corrected consensus BAM per chunk/sample.

## Implementation Checkpoint

Implemented in the working tree on 2026-05-21:

- targeted `realign_exon_blocks()`
- shared hard-clipped corrected-read edit helper
- lazy HP edit-distance scoring from raw BAMs + corrected TSVs
- final `write_corrected_consensus_bam()` writer
- generated chunk-script wiring to omit per-aligner corrected BAMs by default
- post-hoc round-robin chunk sidecar reconstruction API

Validation:

- Local focused suite:
  `145 passed, 1 skipped, 1 xfailed`
- Sherlock focused suite:
  `tests/test_bam_writer_parallel_smoke.py tests/test_read_num_sidecar.py`
  passed with `11 passed`
- Full local `tests/test_validation_reads.py` was blocked by an environmental
  OpenMP shared-memory startup error in subprocess `rectify correct` calls.

## Timing Checkpoint: ysh1 Chunk 000

Sherlock pilot using ysh1 chunk 000, four stable aligners
(`mapPacBio`, `minimap2`, `gapmm2`, `uLTRA`), deSALT excluded:

- Full-chunk lazy merge + whole-read HP scoring: `1055.0s`
- Final corrected consensus BAM writing/sort/index: `227.4s`
- Output: `corrected_reads.tsv`, `aligner_summary.tsv`,
  `corrected_consensus.bam`, and `.bai`

Follow-up fast-path implementation:

- `merge_corrected_tsvs(..., lazy_scoring_workers=N)` parallelizes independent
  aligner-level HP scoring scans.
- Lazy raw-BAM scoring now skips the transient corrected-BAM edit stack for rows
  where the TSV indicates no correction surgery can change the CIGAR.
- Generated split and single-sample merge paths pass available thread counts to
  `lazy_scoring_workers`.

Fast iteration pilot:

- Added a chunk-test `N_READS_TO_TEST` workflow in the ysh1 pilot script.
- `N_READS_TO_TEST=10000` creates temporary per-aligner TSV subsets across
  Module 2-related triggers, then runs the same lazy merge and final consensus
  BAM writer.
- The pilot sampler now supports `MIN_READS_PER_TRIGGER` so small tests first
  reserve a quota for each available trigger class before filling remaining
  reads round-robin.
- ysh1 chunk 000 stratified 10k result:
  - selected 10,000 read IDs
  - lazy merge + HP scoring: `18.8s`
  - final corrected consensus BAM: `14.1s`
  - available trigger coverage in this chunk: softclip HP rescue, overcall
    rescue, changed 3' end, junction reads, and no-trigger reads
  - unavailable in this bypass chunk: 2F five-prime rescue, intronic-tail clip,
    and reanchor/upstream-trim examples

2F-rich targeted pilot:

- RPL19B/RPL20B short-exon-1 genes are good real-data 2F stress cases.
  In current TSVs they appear by systematic IDs:
  - RPL19B: `YBL027W`
  - RPL20B: `YOR312C`
- A scan of full minimap2 corrected TSVs found dense RPL19B 2F examples.
  `rna15_rep3/chunk_001` was selected because all four stable aligners had
  nonempty corrected TSVs and raw BAMs.
- `N_READS_TO_TEST=10000`, `MIN_READS_PER_TRIGGER=1000`,
  `MIN_READS_TARGET_GENES=1000`, `TARGET_GENE_IDS=YBL027W,YOR312C`.
- Trigger coverage in the selected panel:
  - target short-exon genes: 134 / 134 available
  - target-gene 2F: 93 / 93 available
  - 2F five-prime rescue overall: 7,552 / 9,590 available
  - softclip HP rescue: 122 / 122 available
  - intronic-tail clip: 7,277 / 9,227 available
  - overcall rescue: 69 / 69 available
  - reanchor/upstream-trim: 3,183 / 3,279 available
  - changed 3' end: 1,335 / 1,335 available
  - junction reads: 7,928 / 10,647 available
  - no-trigger reads: 2,084 / 4,490 available
- Timing:
  - lazy merge + HP scoring: `18.3s`
  - final corrected consensus BAM: `7.9s`
- Output:
  `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/rna15_rep3/chunks/lazy_consensus_target_2f_test/chunk_001`

Next optimization target:

- Implement differential HP scoring over only the alignment regions where
  candidate aligners disagree. Whole-read HP scores remain correct but waste
  work on shared blocks that cancel out across aligners. A safe implementation
  should first fingerprint final corrected alignments per read, reuse scores for
  identical fingerprints, then restrict penalty calculation to disagreement
  intervals for non-identical candidates.
