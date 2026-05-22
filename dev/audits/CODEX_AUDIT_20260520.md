# Codex Audit: Rectify Coordinate and Pipeline Review

Date: 2026-05-20

Scope: static multi-agent review of Rectify coordinate handling and adjacent
pipeline logic, covering 0-/1-based conventions, clustering, CIGAR surgery,
intron and junction handling, 3' walkback, 5' rescue, consensus selection,
manifest/index analysis paths, bedgraph output, and parallel/streaming modes.

No source files were edited during this audit.

## Coordinate Conventions Observed

- `corrected_3prime` is treated throughout the correction path as a 0-based,
  inclusive aligned-base coordinate.
- Alignment spans from pysam/BAM are 0-based half-open:
  `reference_start` is inclusive and `reference_end` is exclusive.
- GFF/GTF loading converts 1-based inclusive annotation coordinates to
  0-based half-open intervals with `start - 1, end`.
- BED/bedGraph output should therefore emit a corrected 3' base as
  `start = corrected_3prime`, `end = corrected_3prime + 1`.
- CIGAR `N` operations are represented as half-open skipped spans
  `[junction_start, junction_end)`. Boundary fixes must be explicit about
  whether the desired value is the skipped-span boundary or the nearest
  aligned exon base.

## Highest-Priority Findings

### 1. Minus-strand artifact `N` correction can emit intronic coordinates

Refs:

- `rectify/core/bam/bam_processor.py:836`
- `rectify/core/bam/bam_processor.py:850`
- `rectify/core/splice/false_junction_filter.py:282`
- `rectify/core/splice/false_junction_filter.py:300`
- `rectify/core/splice/false_junction_filter.py:423`

For minus-strand false-junction/artifact `N` handling, production code snaps
and clips ambiguity bounds using `junction_start`. Under the half-open
`[junction_start, junction_end)` convention, `junction_start` is inside or at
the skipped span boundary, while the aligned gene-body side for a minus-strand
3' artifact is plausibly `junction_end`.

Impact: `corrected_3prime` or NET-seq refinement windows can land inside an
artifact gap, causing downstream BAM clipping or CPA assignment in an intronic
or skipped region.

Suggested tests:

- Minus-strand artifact `N` near the 3' end where walkback lands inside the
  skipped span.
- NET-seq refinement constrained by the artifact ambiguity window.
- Direct comparison of `false_junction_filter.correct_3prime_for_artifact_junctions`
  and `bam_processor` boundary choices.

### 2. Generic walkback real-`N` fallback can return the first intron base

Ref:

- `rectify/core/correct/walkback.py:756`

The left-side real-`N` fallback sets `true_cpa = first_real_n[0]`. That value is
the first skipped base in a half-open CIGAR `N` span. Since `corrected_3prime`
is an inclusive aligned-base coordinate, this can return a position inside a
real intron rather than the last pre-`N` exon base.

Suggested tests:

- Left-side/negative-strand real-`N` fallback with an all-stop-base terminal
  exon and no non-stop anchor before the `N`.
- Assert that the returned CPA is an aligned exon base, not a skipped base.

### 3. Streaming/parallel TSV schema drops fields required for BAM CIGAR surgery

Refs:

- `rectify/core/bam/parallel.py:42`
- `rectify/core/bam/parallel.py:535`
- `rectify/core/bam/parallel.py:841`
- `rectify/core/bam/bam_writer.py`

Streaming and parallel TSV writers omit fields that `bam_writer` uses to apply
CIGAR edits, including some of:

- `oc_homopolymer_extension`
- `oc_overcall_count`
- `oc_terminal_base`
- `five_prime_upstream_trim`
- `reanchor_clip_len`
- `five_prime_intron_clip_pos`

Because `bam_writer._load_corrections_from_tsv()` only enables these edits when
the columns exist, the same reads can produce different BAM CIGARs depending on
whether correction used regular, streaming, or parallel output.

Suggested tests:

- Schema parity test across `output.write_output_tsv()` and every
  streaming/parallel TSV writer.
- Integration fixture with overcall rescue, upstream trim, reanchor rescue, and
  intron-clip cases, followed by BAM-writing assertions.

### 4. `write_dual_bam()` is not equivalent to the standalone BAM writers

Refs:

- `rectify/core/bam/bam_writer.py:623`
- `rectify/core/bam/bam_writer.py:658`
- `rectify/core/bam/bam_writer.py:690`

`write_dual_bam()` omits behavior present in standalone hardclip/softclip
writers:

- It calls `extend_read_5prime_for_junction_rescue()` without `upstream_trim`.
- It does not call `extend_read_3prime_for_overcall_rescue()` on either output.

Impact: normal paired-output runs can produce CIGARs that differ from running
the hardclip and softclip writers separately.

Suggested tests:

- Equivalence test: `write_dual_bam()` vs `write_corrected_bam()` plus
  `write_softclipped_bam()` for synthetic reads with `five_prime_upstream_trim`
  and `oc_terminal_base`.

### 5. Opposite-strand junctions can drive 5' rescue or consensus snaps

Refs:

- `rectify/core/consensus/select.py:50`
- `rectify/core/consensus/scoring.py:194`
- `rectify/core/splice/splice_aware_5prime.py:1014`
- `rectify/core/splice/splice_aware_5prime.py:1555`
- `rectify/core/splice/junction_scoring.py:218`

Candidate junction pools preserve or receive strand in some places, but rescue
and refinement paths often ignore or drop strand. A plus-strand read can
therefore potentially use a nearby minus-strand junction from an overlapping
antisense gene to cancel a 5' clip penalty, snap a boundary, or drive a rescue.

Suggested tests:

- Same-chromosome, near-boundary opposite-strand annotated junction where the
  clip sequence matches.
- Assert no rescue unless the candidate junction strand matches the read/gene
  strand.
- Include both annotation-derived and pool-derived junctions.

### 6. Consensus canonical splice-site counting is wrong for minus-strand junctions

Refs:

- `rectify/core/consensus/extract.py:175`
- `rectify/core/consensus/extract.py:267`
- `rectify/core/consensus/select.py:79`
- `rectify/core/splice/junction_validator.py:189`

`check_canonical_splice_sites()` has no strand parameter and always checks
genomic `GT/GC...AG`. A canonical minus-strand intron appears as genomic
`CT...AC` under the same coordinate convention. The consensus code can
therefore count canonical minus-strand junctions as non-canonical, weakening
tie-breaks in `select_best_alignment()`.

Suggested tests:

- Strand-aware canonical motif checks for plus and minus introns.
- Consensus tie-break fixture where a minus-strand canonical alignment should
  beat a non-canonical alternative.

### 7. Plus-strand 5' `intronic_snap` is off by one

Refs:

- `rectify/core/splice/splice_aware_5prime.py:1573`
- `rectify/core/splice/splice_aware_5prime.py:253`
- `rectify/core/splice/splice_aware_5prime.py:1351`
- `rectify/core/splice/splice_aware_5prime.py:1532`
- `rectify/core/bam/read_edits.py:1215`

For plus strand, `five_prime_position` is consistently treated elsewhere as
the last upstream-exon base, `intron_start - 1`. The Case 4 `intronic_snap`
path returns `intron_start`, which is inside the intron under the half-open
convention.

Suggested tests:

- Plus-strand Case 4 intronic snap with `align_5prime` inside
  `[intron_start, intron_end)`.
- Assert `five_prime_corrected == intron_start - 1` and validate the final
  CIGAR/intron length.

### 8. Analyze mode can silently choose different coordinate columns

Refs:

- `rectify/core/analyze/loaders.py:82`
- `rectify/core/analyze/loaders.py:181`
- `rectify/core/analyze/manifest.py:137`
- `rectify/core/analyze/manifest.py:471`
- `rectify/core/position_index.py:70`

Standard and manifest full-TSV paths prefer `corrected_position` when multiple
position columns exist. Position-index mode requires/uses `corrected_3prime`.
If a TSV contains both columns and they differ, standard/manifest fallback and
index mode can cluster different coordinates.

Suggested tests:

- TSV fixture with both `corrected_position` and `corrected_3prime` populated
  with intentionally different values.
- Assert standard mode, manifest fallback, and manifest index mode either
  produce identical results or fail loudly.

### 9. Fraction/count semantics differ across analyze paths

Refs:

- `rectify/core/commands/analyze_command.py:227`
- `rectify/core/analyze/clustering.py:202`
- `rectify/core/analyze/manifest.py:146`
- `rectify/core/position_index.py:84`
- `rectify/core/analyze/clustering.py:469`

Small-file clustering and manifest fallback can ignore `fraction`, while index
mode sums it. Count matrices may then use fractional weights even when cluster
formation did not. This makes cluster formation and `min_reads` thresholds
dependent on file size/index availability.

Suggested tests:

- Split-read rows with fractions summing to 1.0.
- Run through standard small-file, large/chunked, manifest fallback, and
  manifest-index paths with identical expected clusters/counts.

### 10. Adaptive clustering can lose valley positions

Refs:

- `rectify/core/analyze/clustering.py:307`
- `rectify/core/analyze/clustering.py:340`
- `rectify/core/analyze/clustering.py:348`
- `rectify/core/analyze/clustering.py:404`

`_find_valleys_between_peaks()` identifies a valley coordinate, but
`_get_adaptive_cluster_boundaries()` places boundaries halfway between peak and
valley. For example, peaks at 100 and 106 with valley 103 can create boundaries
that leave position 103 assigned to no cluster.

Suggested tests:

- Two-peak/one-valley synthetic input.
- Assert every position between retained peaks is intentionally assigned or
  explicitly reported as discarded.

## Additional Findings

### Bedgraph and annotation boundary issues

Refs:

- `rectify/core/analyze/bedgraph.py:70`
- `rectify/core/analyze/bedgraph.py:76`
- `rectify/core/analyze/bedgraph.py:99`
- `rectify/core/analyze/manifest.py:506`
- `rectify/core/analyze/manifest.py:593`
- `scripts/generate_bedgraph_from_polished.py`
- `rectify/core/analyze/clustering.py:531`
- `rectify/core/commands/analyze_command.py:130`
- `rectify/core/analyze/manifest.py:217`

Observations:

- Some bedgraph off-by-one issues have been documented/fixed in adjacent notes,
  but all emitters should be kept under a shared regression test.
- Bedgraph output can ignore pre-aggregated `count` or `fraction` columns in
  some paths.
- Hard-coded yeast chromosome order can drop non-yeast or scaffold contigs.
- Plus-strand gene 3' distance uses half-open `gene['end']` instead of
  `gene['end'] - 1`.
- rDNA exclusion uses `<= end` against half-open annotation intervals and
  excludes one extra base at the right edge.

### Manifest/index analysis divergence

Refs:

- `rectify/core/analyze/manifest.py:321`
- `rectify/core/analyze/manifest.py:452`
- `rectify/core/analyze/manifest.py:466`
- `rectify/core/analyze/manifest.py:638`
- `rectify/core/analyze/manifest.py:646`
- `rectify/core/analyze/genomic_distribution.py:1079`
- `rectify/core/analyze/genomic_distribution.py:1107`

Observations:

- Manifest genomic-distribution output is run per condition but writes fixed
  filenames, so later conditions overwrite earlier ones.
- Manifest index mode skips TSS count assignment for indexed samples.
- `--min-cluster-samples` is parsed/read but not applied.

### Annotation feature loss

Refs:

- `rectify/core/analyze/loaders.py:407`
- `rectify/core/analyze/loaders.py:445`
- `rectify/core/analyze/genomic_distribution.py:146`

`_parse_gtf()` keeps only `gene` rows and does not emit the feature-level data
expected by genomic-distribution code (`UTR3`, `UTR5`, `CDS`, `CUT`, `SUT`,
`XUT`, `snoRNA`, etc.). Parsed GTF/GFF therefore falls back to broader gene
body heuristics.

### Chromosome aliasing can disable walkback or `=` decoding

Refs:

- `rectify/core/bam/bam_processor.py:279`
- `rectify/core/bam/bam_processor.py:767`
- `rectify/core/bam/bam_processor.py:790`
- `rectify/core/bam/bam_processor.py:947`

Walkback and `=` sequence decoding only try limited contig-name lookups in some
paths. A BAM using `chrI` with an NCBI-keyed FASTA can silently disable
DRS/QuantSeq REV walkback or leave `=` characters undecoded in query sequence.
The post-walkback NET-seq poly-A-base filter can also be skipped when walkback
uses a raw contig key rather than `chrom_std`.

### Non-walkback ambiguity invariants can become inconsistent

Ref:

- `rectify/core/bam/bam_processor.py:902`

Non-walkback position shifts rewrite `ambiguity_min`/`ambiguity_max` but do not
always recompute `ambiguity_range`. This can skip NET-seq refinement despite a
nonzero window, or produce invalid min/max ordering in some softclip/overcall
rescue cases.

### Parallel checkpoint/resume can duplicate or drop rows

Refs:

- `rectify/core/bam/parallel.py:880`
- `rectify/core/bam/parallel.py:889`

Region rows are written to the shared output before the `.done` sentinel is
touched, without atomic per-region output or fsync. A crash between those steps
can duplicate rows on resume; a crash after `.done` but before durable data can
drop rows.

### Parallel BAM writer instability

Ref:

- `rectify/core/bam/bam_writer_parallel.py:389`

During this audit, `pytest tests/test_bam_writer.py tests/test_bam_writer_parallel_smoke.py -q`
passed low-level writer tests and the `n_threads=1` smoke case, then repeatedly
hit fatal Python abort traces in `multiprocessing.Pool` while exercising
`write_corrected_bam_parallel()`. The run was interrupted after 16 passes.

### Corrected-consensus effective grouping ignores chrom/strand

Refs:

- `rectify/core/consensus/corrected_consensus.py:695`
- `rectify/core/consensus/corrected_consensus.py:801`

`_eff_key()` groups by `(corrected_3prime, junctions)` only. Two aligners on
different chromosomes or opposite strands can be reported as effectively
matching the winner when their numeric CPA and junction string collide. This
may not select the winner, but it corrupts summary interpretation.

### Consensus docs still mention A-tract depth as a scoring factor

Refs:

- `rectify/core/consensus/consensus.py:42`
- `docs/user_guide/commands/consensus.md:11`
- `rectify/core/consensus/scoring.py:724`

Docs/header text still describe A-tract depth as a consensus scoring factor,
while `score_alignment()` explicitly removed it from scoring. This is a user
contract mismatch.

## Suggested Fix Order

1. Normalize all `N`/junction boundary conventions in correction and walkback.
   Start with minus-strand artifact `N` handling and real-`N` fallback.
2. Make output TSV schemas identical across regular, streaming, and parallel
   correction paths.
3. Make `write_dual_bam()` call the same CIGAR edit logic as the standalone
   hardclip and softclip writers.
4. Enforce strand-aware junction matching for consensus, 5' rescue, 3'SS rescue,
   and annotation-derived junction pools.
5. Add strand-aware canonical splice-site scoring to consensus extraction.
6. Lock analyze standard/manifest/index equivalence with shared fixtures for
   `corrected_3prime`, `corrected_position`, `fraction`, and `count`.
7. Add regression tests for bedgraph coordinate output and annotation
   half-open boundary handling.
8. Add crash/resume and multiprocessing smoke tests for parallel correction and
   parallel BAM writing.

## Test Coverage To Add

- Boundary tests for `corrected_3prime` on both strands around CIGAR `N`, `D`,
  and terminal soft/hard clips.
- Minus-strand false-junction artifact classification and correction.
- Real-intron walkback fallback with no non-stop anchor before/after `N`.
- Strand-aware canonical splice motif tests for plus and minus introns.
- Opposite-strand junction-rescue rejection tests.
- `write_dual_bam()` equivalence tests against standalone BAM writers.
- TSV schema parity tests for all correction output modes.
- Analyze equivalence tests across standard, large/chunked, manifest fallback,
  and manifest-index modes.
- Bedgraph tests for 0-based corrected positions, `count`, `fraction`, and
  non-yeast contigs.
- Manifest genomic-distribution multi-condition tests to prevent overwrite.
- Parallel checkpoint crash-injection tests.
