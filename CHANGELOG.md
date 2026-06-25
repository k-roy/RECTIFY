<!-- markdownlint-disable MD022 MD024 MD025 MD032 -->
# Changelog

All notable changes to RECTIFY will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`rectify cdna-analyze` now emits a per-molecule `corrected_reads.tsv`**
  (`cdna_analyze_command.py`) — one row per UMI-consensus molecule, using the
  same per-read column schema as the DRS `rectify correct` output, so a single
  loader works across modalities.
- **`--use-dorado-polya` flag** (`rectify correct`, default off): when set, the
  Dorado `pt:i` tag is used as the authoritative poly(A) length. The Dorado
  estimate is always recorded as `dorado_polya_length` regardless of this flag
  (`polya_trimmer.py`, `correct_command.py`).
- **Consensus aligner-selection stats surface `by_aligner_combo` and
  `tied_score`** in `*.consensus_aligner_stats.tsv` and the HTML consensus
  report (`consensus/`, `analyze/summary.py`).
- **`rectify analyze` reports a read-weighted center-of-mass (`cluster_com`)
  per CPA cluster** in `cpa_clusters.tsv` (`analyze/clustering.py`).
- **DRS pre-trim now strips the terminal `(AAG/GAA)ₙ` triplet-repeat basecaller
  artifact** (`rectify trim-polya`). Dorado v5.2.0 / RNA004 mis-basecalls the poly-A
  homopolymer tail as a low-period multi-base repeat with no terminal A-run, so the
  poly-A scan misses it and the artifact is force-aligned downstream. The trimmer now
  peels a terminal low-period **multi-base** repeat block before the poly-A scan
  (`_find_terminal_repeat_block`), exposing any genuine poly-A behind it. The boundary
  is anchored on the repeat's non-A phase, so a poly-A run 5′ of the artifact is
  preserved (poly-A homopolymers are period-consistent at every period and would
  otherwise be merged into the block). Reuses the existing motif-agnostic
  `repeat_expansion` detector — strict no-op on clean poly-A tails (no motif gating).
  - New flags: `--strip-repeat-expansion` / `--no-strip-repeat-expansion` (default ON),
    `--repeat-min-len` (default 15), `--repeat-min-frac` (default 0.8).
  - New trim-metadata columns: `repeat_len`, `repeat_motif`.
  - Measured on Sumner SMA `CNTL_21.8` (RNA004) at the default `min_len=15`: terminal-GAA
    fraction on trimmed 3′ ends 42.6% → 3.5%; strict no-op on the yeast Cat1–9 validation
    reads (0/36 fire).

### Fixed

- **DRS 3' ends now land on a non-A 100% of the time** (walkback guard-refactor +
  soft-clip-rescue stop-base fix). RECTIFY's DRS policy left-shifts every corrected
  3' end to the first non-A (the pre-trim removes the poly-A tail, so a corrected end
  parked on a genomic A/T — an RNA A — is always a bug). Three changes:
  - **Removed the `early_exit_homopolymer_check` suppressor** from the poly-A walkback
    (`walkback.py`). Its "no genomic A-tract near the 3' end ⇒ skip the walk" heuristic
    was invalid (the poly-A can be basecall-truncated) and was the recurring source of
    the +strand lone-A regression. The walkback now always scans; the artifact guards
    only *restrict* where the anchor may land, never suppress the walk.
  - **Added an `_enforce_non_stop_anchor` post-condition** that re-runs the unguarded
    core (uncapped) if the guarded walk would leave the end on a stop base and no
    artifact guard clipped the scan — making the on-stop-base regression class
    structurally impossible while preserving every artifact rejection.
  - **Removed a soft-clip-rescue special case** (`indel_corrector.rescue_softclip_at_
    homopolymer`, −strand) that could extend a 3' end through the poly-A tail onto a
    genomic T (RNA A).
  - Added CI gate `TestCorrectedEndsAreNonA` asserting zero corrected DRS 3' ends on a
    gene-strand stop base (strand-skew enumerated as the regression canary).
  - **Universal non-A post-condition in `bam_processor`**: after every correction module
    settles the 3' end, if it lands on a gene-strand stop base it is re-anchored via the
    guard-respecting walkback (or reverted to the raw 3' end), with rescue CIGAR-surgery
    metadata reverted for coherence. This catches the rescue modules that previously
    extended onto a genomic A by design (notably variant-aware `homopolymer_rescue`).
    Verified on human DRS (Sumner chr5, 11 samples): corrected-on-A 7.09% → 0.0056%.

## [0.9.0] - 2026-05-16

**First public release of RECTIFY.** Prior 2.x and 3.x versions were internal
pre-public development; their history is preserved below for completeness but
should be read as "what landed before 0.9.0 shipped publicly."

### Changed

- **Universal `--Scer` / `--organism` shorthand** (commit `3454418`): every
  reference-aware subcommand now accepts `--Scer` (= `--organism
  saccharomyces_cerevisiae`) and resolves the bundled genome / GFF / GO
  annotations automatically. A single shared helper
  (`rectify.data.add_organism_args`) attaches the argparse pair, and a single
  CLI hook (`rectify.data.resolve_reference_paths`, called once in
  `rectify.cli.main()` between `parse_args()` and dispatch) fills `args.genome`
  / `args.annotation` / `args.go_annotations` from the bundled paths. Covers:
  `align`, `analyze`, `aggregate`, `extract`, `validate`, `netseq`, `consensus`,
  `train-polya`, `prescan`, `split`, `export`, `correct`, `run-all`, `batch`.
  `--genome` is relaxed from `required=True` to `default=None` on the commands
  where it was strict (`align`, `consensus`, `netseq`, `prescan`,
  `train-polya`). 46 new tests in `tests/test_organism_helper.py` cover the
  helper, each subcommand's `--Scer` wiring, and the CLI main() hook
  end-to-end.
- **Internal documentation restructure**: separated user-facing docs (`docs/`)
  from contributor-only material (new `dev/` tree). Promoted `CHANGELOG.md` to
  repo root per OSS convention; mkdocs site re-exposes it via a snippet stub.
  Pruned 10 stale files (7 completed refactor briefs + their index + the
  workshop walkthrough + a stale README clone).
- **Internal refactor — split seven 2000+-line files** into focused sub-modules
  under `bam/`, `cdna/`, `commands/run/`, `consensus/`, `splice/`. No public
  API change; old paths preserved as re-exports. Commits: `0024fa3`
  (bam_processor), `c411c80` (cdna), `49478a0` (analyze), `77d3828` (bam_writer),
  `33c683c` (consensus), `bf55e96` (junction_refiner), `b3123fe` (run_command).

### Added

- **Poly(A) evidence tracking** (`core/bam_processor.py`, `core/correct_command.py`): `corrected_3ends.tsv` now includes three new columns — `pt_tag` (dorado signal-level poly(A) estimate, integer or empty), `polya_score` (poly(A) model confidence 0–1, four decimal places, empty if model not loaded), and `polya_source` (`'pt_tag'` | `'model'` | `'none'`). Both streaming and parallel BAM processing modes accept `--polya-model`; model is loaded once and applied to all reads.

- **`rectify tag-polya` subcommand** (`cli.py`, `core/tag_polya_command.py`): retroactively annotates an aligned BAM with `ps:f` (poly(A) model score) and optionally `pt:i` (when absent). Preserves pre-existing dorado `pt:i` values.

- **Unaligned dorado BAM auto-detection** (`core/preprocess.py`): when the user supplies an unaligned BAM (all primary reads unmapped, as produced by `dorado basecall`), `prepare_input()` automatically routes to `align_ubam_to_genome()` which runs `samtools fastq -T '*' | minimap2 -y | samtools sort` to preserve aux tags (including `pt:i`) through alignment.

- **`--max-cluster-radius`, `--min-peak-sep`, `--min-cluster-samples` CLI flags** (`core/analyze_command.py`): the three remaining APA clustering constants are now user-configurable. Passing any non-default value activates `cluster_cpa_sites_adaptive()` (valley-based) in place of the fixed-distance `cluster_cpa_sites()`.

- **`tests/test_terminal_exon_refiner.py`** (Bug 37): 51 tests covering `SpliceSiteIndex`, `load_splice_sites_from_gff` (plus/minus strand coordinate conversion, GFF position conflicts), `merge_splice_indices`, `detect_junction_truncated_reads`, `detect_partial_junction_crossings`, `get_soft_clip_info`, `simple_align`, `detect_mismatch_clusters`. Real-data classes validate YAL030W and YAL001C intron positions from the bundled R64-5-1 GFF and wt_by4742_rep1 BAM.

- **`tests/test_consensus_selection.py`** (Bug 38): 40 tests covering `extract_junctions_from_cigar`, `check_canonical_splice_sites` (GT/AG, GC/AG, real YAL030W intron), `score_alignment` (5' clip −2/base, A-tract depth cap, 3' clip penalty), `select_best_alignment` (winner selection, `was_5prime_rescued`, tiebreakers, confidence levels). Real-data class validates junction extraction and aligner selection at the YAL030W locus.

### Fixed

- **`rescue_3ss_truncation` Bug 3: mapPacBio terminal-D overshoot, Pattern 2** (`core/splice_aware_5prime.py`): added `_leading_del` detection for CIGAR pattern `… N D =` (reversed 5' ops: `= D N`), where mapPacBio encodes a few missed exon bases as deletions between the terminal match block and the intron N-op. The new `_5p_ops` collector reads up to 3 non-H ops from the 5' end; `_leading_del=True` fires when `op[0] ∈ {M,=,X}`, `op[1]==D`, `op[2]==N`. Pattern 1 (plain terminal D at 5' end, `op[0]==D`) was already handled but is now detected via the same collector for consistency. When both `_n_match` and `_leading_del`, the junction is recorded as `_forced_snap_junction`; `n_op_snap` fires if all sequence-based rescue attempts fail. Verified on real read `de69692f-04e9` (chrII minus strand, N-op `[89132,89436)`): `five_prime_corrected=89440`.

- Three new regression tests in `TestRescue3SSTruncation`: `test_minus_n_op_snap_pattern1_terminal_d`, `test_minus_n_op_snap_pattern2_d_between_n_and_match`, `test_minus_clean_5prime_no_snap`.

- **`_score_junction` degenerate k=L case** (`core/splice/junction_refiner.py`): the k-loop previously ran `range(L+1)`, where k=L gives an empty `rescue[L:]`, making the score 0.0 for every candidate regardless of junction quality — a completely non-discriminating case that caused the wrong junction to win whenever candidates tied at 0.0. Fixed by changing to `range(L)` so `q1 = rescue[k:]` always contains at least one base.

- **`refine_read_junctions` stability tiebreaker** (`core/splice/junction_refiner.py`): added `is_alt` as the second element of the candidate tuple (after `score`, before `canonical_tier`). `is_alt=0` when the candidate matches the current N-op boundaries exactly; `is_alt=1` otherwise. This prevents equal-scoring candidates from spuriously displacing already-correct junctions — critical when multiple candidates share the same `intron_end` (e.g. TFC3 annotated junction and a nearby alternative both score 0.0 at the same `je`).

- **Impact**: 9 previously failing tests in `tests/test_junction_refiner.py` now pass. All 698 tests pass (4 skipped). Key reads fixed: RPL20B read `0b3b593b` now correctly refined from `[900758,901189)` to `[900767,901193)`; TFC3 annotated reads no longer displaced by same-score alternative at `[150989,151096)`.

---

## Pre-public internal history

Earlier 2.x and 3.x version entries describe internal pre-public development.
Module paths cited in those entries refer to the layout that existed at the time
of each fix; some paths have moved to subpackages under `core/` since (see the
0.9.0 refactor note above). The bug rationales and behaviour descriptions remain
accurate; only the file paths may be stale.

## [3.2.5] - 2026-04-24

### Changed
- Validation Cat1/Cat2 reads replaced with DRS-trimmed examples showing genuine correction artifacts
- Cat1 (4 reads): DRS minimap2 reads on chrXIV/chrI/chrII/chrXII with `atract_ambiguity,indel_correction,polya_walkback`; shifts: −16, −4, +3, +3 bp
- Cat2 (3 replaced, 9dbd37bf kept): DRS minimap2 reads on chrI/chrVI/chrV with `softclip_rescue`; shifts: +11, +4, −8 bp
- **708 tests passing** (Cat1/Cat2 expected positions updated to new chromosomes/positions)

## [3.2.3] - 2026-04-24

### Fixed
- Validation reads Cat5–9 sequences synced to DRS-trimmed run (v3.1.8 updated CIGAR/boundaries but left old chunked sequences)
- cat5_minus_1 `SEQ='*'` filled with chunked consensus sequence (chimeric CIGAR incompatible with DRS seqlen)
- Cat6/Cat7/Cat9 sequences + CIGAR from DRS mapPacBio merged BAM

## [3.2.2] - 2026-04-22

### Fixed
- XV/XG tag fix: both tags now present on all 36 reads in `validation_reads.bam`; XV required by test fixture, XG drives category assertions
- N-op fix for 3 replacement reads (f8050895, 7d5e8dc2, 72557a9a) re-sourced directly from dev run BAMs; correct intervals: f8050895 all-aligners (45644–45977), 7d5e8dc2 mapPacBio-only (60193–60697), 72557a9a mapPacBio-only (104435–104495)
- Full pipeline re-run: `corrected_3ends.tsv`, rectified BAMs, `rectified_pA_tail_soft_clipped.bam` all regenerated from fixed dataset
- Independent certification: 4 independent agent reviews confirmed 36/36 PASS

## [3.2.1] - 2026-04-22

### Fixed
- `rescue_3ss_truncation` soft-clip exon CIGAR body-borrowing: for `softclip` rescues, use `rescue_seq` (exact soft-clip) instead of `_intronic_seq` as alignment query; cat3_plus_2 (79f61403, YAL003W) now produces `14M1D7M1D1M` instead of flat `22M`

## [3.2.0] - 2026-04-22

### Fixed
- Validation aligner BAMs rebuilt (5 BAMs) after DRS rebuild removed 3 old Cat6/Cat7 reads but left aligner BAMs stale

## [3.1.9] - 2026-04-22

### Added
- `rectified_pA_tail_soft_clipped.bam` generated as default Step 4 output
- Sort + index added to both `_process_one_sample` and `_run_single_sample` Step 4 code paths (previously written unsorted/unindexed, unusable in IGV)
- `TestPolyASoftClippedBam` added (4 tests)

## [3.1.8] - 2026-04-22

### Changed
- Validation Cat6/Cat7/Cat9 (12 reads) updated to DRS pre-trim mapPacBio alignments
- 3 reads replaced entirely (pre-DRS alignments lost N-op after trimming): ba761413→f8050895, 64f4da08→7d5e8dc2, 5c59f0bc→72557a9a
- BAM header reduced to 16 nuclear chromosomes (chrI–chrXVI)

## [3.1.7] - 2026-04-21

### Changed
- Module 2H `refine_read_junctions`: bilateral t2 scoring — `score(k) = t1(k) + t2(k)`, both anchored to `intron_end`; eliminates degenerate k=L-1 coincidental single-base matches
- **No candidate guards (permanent policy)**: removed `if is_alt==1 and tier>=4 and is_novel==1: continue`; annotation/canonical tier are tiebreakers only, never gates
- Canonical HP prior (`_CANONICAL_HP_PRIOR = 0.5`): replaces fragile `int(score)` floor; canonical junctions win within one Nanopore HP deletion noise floor regardless of penalty table

## [3.1.6] - 2026-04-21

### Added
- `--drs` flag wired into `rectify run-all` (cli.py); enables Step 0 (poly-A pre-trim) and Step 4 (poly-A restore) for Dorado BAM inputs

### Fixed
- ARCHITECTURE.md mermaid diagram: Step 3=analysis, Step 4=restore-softclip (was Step 4=analysis, Step 5=restore)
- README.md: `rectify run` → `rectify run-all`; Module 2H scoring description updated to HP-anchored semi-global DP

## [3.1.5] - 2026-04-21

### Added
- `--checkpoint-dir DIR` flag for `rectify correct --streaming`; two-level checkpoint/resume: scan-phase pickle + per-region `.done` sentinels

## [3.1.4] - 2026-04-21

### Fixed
- Module 2H candidate guard: non-canonical + non-annotated alternatives excluded during candidate *selection* instead of post-selection (preventing spurious wins)
- Adaptive tiebreak ordering: `(score, tier, is_alt, ...)` when current junction is non-canonical; `(score, is_alt, tier, ...)` when canonical

## [3.1.3] - 2026-04-21

### Fixed
- `--aligner-bams` `aligner:path` prefix stripping: `_strip_aligner_prefix()` removes `minimap2:` etc. before pysam path lookup; previously caused "Protocol not supported" error

## [3.1.2] - 2026-04-21

### Added
- Cat9 validation reads (4 reads): Module 2H junction refinement test cases; `TestCategory9JunctionRefinement`
- `validation_reads.bam` expanded 32 → 36 reads; all 5 aligner BAMs updated

## [3.1.1] - 2026-04-21

### Fixed
- `_score_junction` k-loop `range(L+1)` → `range(L)`: eliminates degenerate k=L empty-query scores
- `is_alt` tiebreaker: existing junction preferred over equal-scoring alternatives

## [3.1.0] - 2026-04-20

### Added
- **Module 2H** (`junction_refiner.py`): post-consensus N-op junction refinement
  - HP-aware split-alignment scores all candidates within `search_radius`; canonical/annotated are tiebreakers only
  - Fast path: tier-0 annotated canonical reads skip scoring (255× speedup)
  - CIGAR surgery encodes boundary changes as I/D ops; MD/cs tags stripped
  - Wired into `rectify correct` via `--aligner-bams`
  - Flags: `--junction-hp-pen`, `--junction-search-radius`, `--junction-window`, `--junction-max-slide`, `--junction-max-boundary-shift`
  - 41 tests in `tests/test_junction_refiner.py`

## [3.0.4] - 2026-04-20

### Fixed
- `rescue_3ss_truncation` minus-strand soft-clip truncation: `rescue_seq[-five_clip:]` (was `[:five_clip]`); fixes Cases 1/2 where rightmost bases are the actual soft-clip

## [3.0.3] - 2026-04-16

### Fixed
- `find_polya_boundary` trailing-base false-stop guard: before accepting a stop where terminal base matches genome, check K=4 upstream positions; if all poly-A context, continue scanning

## [3.0.2] - 2026-04-16

### Fixed
- `clip_read_to_corrected_3prime` / `softclip_read_to_corrected_3prime`: strip trailing D/N ops (plus strand) or leading D/N ops (minus strand) before appending H/S clip; fixes `4D6H` invalid CIGAR when corrected position falls inside deletion span
- `tests/test_bam_writer.py` added (11 tests)

## [3.0.1] - 2026-04-15

### Fixed
- `clip_intronic_tail_5prime` off-by-one (minus strand): `<= clip_boundary + 1` → `<= clip_boundary`
- Trailing I/S stripping at 5' end before main trim loop
- Existing H preservation: pre-existing H extracted and merged back after surgery
- `_MIN_SC_FOR_JUNCTION_EXTENSION = 3`: 1–2 bp soft clips fall through to intronic-tail clip instead of spurious Cat3 extension

## [3.0.0] - 2026-04-15

### Fixed
- `clip_intronic_tail_5prime`: `query_sequence` and `query_qualities` now trimmed to match new CIGAR; previously only CIGAR updated → malformed BAM records → `samtools sort: truncated file`
- `five_prime_intron_clip_pos` assignment generalised to all rescues (Cases 1/2/4) where alignment falls inside rescued intron with no soft clip

## [2.9.1] — 2026-04-12

### Fixed

- **Cat2 soft-clip rescue stops at poly-A tail bases** (`core/indel_corrector.py`): `rescue_softclip_at_homopolymer()` now halts rescue when it encounters an `A` base in the soft-clip (plus strand) or a `T` base (minus strand). Previously it would match poly-A tail bases to genomic A-runs, causing the corrected 3' position to be shifted downstream into the poly-A tail. Fixes shifted corrected positions for affected Cat2 reads (e.g. +10→+9, -17→-10, -12→-11 in validation set).

- **`extend_read_3prime_for_softclip_rescue()`** (`core/bam_writer.py`): new function converts the 3' soft-clip of a Cat2-rescued read from `{M}S` to `{D}D{M}M{poly-A}H|S`, making the true RNA 3' end visible in IGV. The deletion operations represent the skipped reference homopolymer; the match operations represent the rescued genomic bases.

- **Cat2 rescue metadata in TSV** (`core/bam_processor.py`): three new columns written to `corrected_3ends.tsv` — `sc_homopolymer_extension` (bases skipped in reference homopolymer), `sc_rescued_seq` (rescued bases), `sc_original_softclip_len` (original soft-clip length before rescue).

- **Bundled BAMs renamed** (`rectify/data/`): `rectified.bam` → `rectified_corrected_3end.bam`; `rectified_softclip.bam` → `rectified_pA_tail_trimmed.bam`. Names now reflect the poly-A handling mode rather than the correction step.

## [2.9.0] — 2026-04-12

### Added

- **`rectify split`** (`core/split_command.py`): splits a FASTQ/FASTQ.GZ into N equal chunks at read boundaries using round-robin interleaving, giving each chunk an even read-length distribution. `--generate-slurm` writes `run_array_align.sh` (SLURM `--array=0-{N×M-1}`) and `run_merge_and_consensus.sh` alongside the chunk files. Task IDs decode to `(chunk_idx, aligner)` via modulo/division; thread limits (`OMP_NUM_THREADS`, `LOKY_MAX_CPU_COUNT`, etc.) are set automatically.

- **`rectify consensus`** (`core/consensus_command.py`): runs consensus aligner selection on pre-built per-aligner BAMs. Accepts `aligner:path` positional arguments (e.g. `minimap2:sample.minimap2.sorted.bam`). Used as the final step after `samtools merge` in the chunked-alignment workflow. Writes a coordinate-sorted, indexed, MD-tagged `<prefix>.consensus.bam`.

- **`rectify install-aligners`** (`core/install_aligners_command.py`): downloads and installs external aligners. `--check` reports PATH availability of all five aligners including the vendored deSALT binary. `--all` installs everything possible. `--desalt` copies the vendored Linux/x86_64 binary to `~/.rectify/bin/` or builds from GitHub source for other platforms. `--minimap2` downloads the pre-built Linux/x86_64 release binary. `--gapmm2`/`--ultra` install via pip.

- **Vendored deSALT binary** (`rectify/data/bin/linux_x86_64/deSALT`): deSALT v1.5.6 (bioconda build `h577a1d6_7`, 837 KB, glibc 2.6.32+) is bundled with the package. `_get_vendored_desalt()` in `multi_aligner.py` resolves the binary automatically when `deSALT` is not on `PATH`, matching platform and architecture. Error message now points to `rectify install-aligners --desalt` for unsupported platforms.

## [2.7.8] - 2026-04-09

### Fixed

- **NEW-011 — Donor/acceptor labels swapped for minus-strand junction validation** (`core/analyze/junction_validation.py`): the minus-strand block derived `donor_seq` from the acceptor genomic window and vice versa, causing every minus-strand GT-AG junction to fail validation. Fixed by re-fetching both sequences from their correct genomic coordinates before reverse-complementing.

- **NEW-012 — Junction aggregate never RC'd donor/acceptor dinucleotides for minus strand** (`core/aggregate/junctions.py`): `five_ss_dinuc` and `three_ss_dinuc` were extracted from raw genomic coordinates with no strand check, misclassifying every minus-strand junction as non-canonical. Added a `_reverse_complement` helper and strand-aware extraction (coordinate swap + RC for minus strand).

- **NEW-013 — NET-seq minus-strand 3' position double-subtracted by `n_trimmed`** (`core/netseq_bam_processor.py`): line 416 subtracted `n_trimmed` a second time from the already-corrected minus-strand position, shifting all NET-seq minus-strand 3' ends upstream by `2 × n_trimmed`. Fixed by changing the second operation to `+ n_trimmed`.

- **NEW-001 — Indel rescue functions ignore prior corrected position** (`core/indel_corrector.py`): `rescue_softclip_at_homopolymer` and `rescue_mismatch_inside_homopolymer` always re-anchored to `read.reference_end - 1`, discarding upstream corrections. Added `current_pos: Optional[int] = None` parameter to both functions; call sites now forward the previously corrected position.

- **NEW-002 — Minus-strand exon sequence fetched from wrong genomic window in 5' splice rescue** (`core/splice_aware_5prime.py:609`): `genome_seq[intron_end - rescue_len:intron_end]` fetched the end of the intron instead of exon1. For minus-strand genes, exon1 (upstream in transcript) is at higher genomic coordinates; `intron_end` marks where exon1 begins in genomic space. The BAM soft-clip for minus-strand reads is in forward-strand orientation (no RC needed), so the correct window is `genome_seq[intron_end:intron_end + rescue_len]`. Fixes cat3_minus_1/cat3_minus_2 test failures where 5' position was left unchanged after splice rescue.

- **NEW-008 — Shell injection via unescaped sample IDs/paths in generated SLURM scripts** (`core/batch_command.py`): sample IDs and BAM paths were embedded in generated bash scripts without `shlex.quote()`, allowing metacharacters from manifest TSVs to inject arbitrary shell commands. All manifest-derived values in generated scripts now use `shlex.quote()`.

- **NEW-035 — `NCBI_TO_CHR` reversal overwrites `chrI` form with bare `I`** (`utils/chromosome.py:64`): the dict-comprehension reversal of `YEAST_CHR_TO_NCBI` iterated over both `'chrI'` and `'I'` forms; `'I'` overwrote `'chrI'` in the result. Fixed by filtering to only `chr`-prefixed keys in the reversal.

- **NEW-007 — Scratch teardown rsync runs unconditionally when `SCRATCH_DIR` is empty** (`core/batch_command.py:162`): the teardown block rsync'd even when scratch was not created, potentially clobbering Oak output with an empty tree. Added `if [ -n "${SCRATCH_DIR}" ] && [ -d "${SCRATCH_DIR}" ]` guard.

- **NEW-004 — mapPacBio `sort_proc.communicate()` timeout uncaught; process hangs indefinitely** (`core/multi_aligner.py:385`): wrapped `communicate(timeout=ALIGNER_TIMEOUT)` in `try/except subprocess.TimeoutExpired` with kill + re-raise, matching the existing minimap2 pattern.

- **NEW-005 — mapPacBio `view_proc` missing `stderr=subprocess.PIPE`; deadlock risk** (`core/multi_aligner.py`): added `stderr=subprocess.PIPE` to the mapPacBio `view_proc` Popen call to prevent stderr buffer fill deadlock.

- **NEW-031 — Junction index missing strand; opposite-strand junctions merged** (`core/analyze/junction_analysis.py`): junction index key `(chrom, donor, acceptor)` promoted to `(chrom, strand, donor, acceptor)`. Updated return type annotation.

- **NEW-032 — GFF junction strand parsed but not stored in junction dict** (`core/analyze/junction_analysis.py:473,501`): `strand = parts[6]` was discarded; exon tuple now includes strand as 4th element, carried through to the index key.

- **NEW-018 — `None` positions from failed extraction propagate to TSV and arithmetic** (`core/bam_processor.py:191,253`): `get_read_3prime_position()` and `get_read_5prime_position()` can return `None` for unmapped reads. Added guards after each call; reads with `None` positions are now skipped with a warning.

- **NEW-019 — `dist <= 0` vs `dist < 0` asymmetry in consensus junction rescue** (`core/consensus.py:378,404`): minus-strand used `dist <= 0`, excluding valid boundary junctions where `dist == 0`. Changed to `dist == 0` (exact match) for both strands for clarity; `dist < 0` is used to reject junctions on the wrong side.

- **NEW-020 — Candidate junctions missing strand; opposite-strand junctions corrupt 5' rescue** (`core/consensus.py:899-907`): junction tuples promoted from `(chrom, start, end)` to `(chrom, start, end, strand)` throughout the pipeline. Rescue loop now filters to matching strand.

- **NEW-021 — mapPacBio `view_proc.returncode` unchecked after sort completes** (`core/multi_aligner.py`): added `view_proc.wait()` followed by returncode check after `sort_proc.communicate()`, raising `RuntimeError` on failure.

- **NEW-022 — `--continue-on-error` undeclared in `run-all` subparser** (`cli.py`): added `--continue-on-error` argument to the `run-all` subparser so the flag is properly parsed and forwarded to `_run_multi_sample()`.

- **NEW-023 — `--use-scratch` undeclared in `run-all` subparser** (`cli.py`): added `--use-scratch` argument to the `run-all` subparser; scratch staging can now be explicitly controlled.

- **NEW-024 — No duplicate `sample_id` detection in manifest** (`core/batch_command.py`): `parse_manifest()` now raises `ValueError` listing all duplicated `sample_id` values before any processing begins.

- **NEW-026 — No minimum-condition guard before primary DESeq2 invocation** (`core/analyze_command.py:351`): added `condition.nunique() >= 2` pre-flight check in both `run_analyze()` and `_run_analyze_manifest()`; single-condition runs log a warning and skip DESeq2.

- **NEW-028 — `--polya-model` path captured but never forwarded to BAM processor** (`core/correct_command.py:206`): `polya_model_path` is now passed to both `process_bam_streaming()` and `process_bam_file_parallel()` call sites.

- **NEW-029 — Empty extract output indistinguishable from success** (`core/extract_command.py:454-457`): `run_batch_extract()` now emits `logger.warning(...)` when the output DataFrame is empty, prompting the user to check filters and input data.

- **NEW-034 — Bare `except Exception` in motif extraction silently skips sequences** (`core/analyze/motif_discovery.py:122`): replaced with `except Exception as e: logger.warning(...)` so sequence-level failures are visible in logs.

- **NEW-036 — `sync_to_oak()` silently falls back after partial rsync failure** (`slurm.py`): `CalledProcessError` from rsync now logs an error and raises `RuntimeError`; `FileNotFoundError` (rsync absent) falls back to shutil but logs a `WARNING`.

- **NEW-037 — `--aligner` choices inconsistent between `correct` and `run-all` subparsers** (`cli.py:106,607`): unified to `['minimap2', 'bwa', 'star', 'auto']` with default `'auto'` in both subparsers.

- **NEW-038 — `n_reads_total` inflated by counting reads before `min_reads` filter** (`core/train_polya_command.py:290`): counter increment moved to after the `min_reads` guard so only reads from sites that pass the filter are counted.

- **NEW-009 — `tes_tolerance == 0` causes `ZeroDivisionError` in APA detection** (`core/analyze/apa_detection.py:304`): added guard; when `tes_tolerance == 0`, the exact position is used as the bin key and a warning is logged.

- **NEW-010 — Isoform grouping key missing strand; cross-strand isoforms merged** (`core/analyze/apa_detection.py:307`): added `chrom` and `strand` to the isoform group key.

- **NEW-039 — Dead `_df` assignment causes input file to be read twice** (`core/bam_processor.py:694`): removed the dead `_df = pd.read_csv(...)` assignment; file is now read only once.

- **NEW-040 — `max(five_clip, five_clip + terminal_end)` tautology in consensus** (`core/consensus.py:544`): replaced with `return five_clip + terminal_end` directly (the max is always the second argument when `terminal_end > 0`).

- **NEW-041 — Bare `except Exception` in `_get_effective_5prime_clip`** (`core/consensus.py`): changed to `except Exception as e: logger.warning(...)` so unexpected failures are logged before returning the default clip length.

- **NEW-042/NEW-043 — `clip_len == 0` causes `ZeroDivisionError` in poly(A) trimmer** (`core/polya_trimmer.py:532,545`): added `if clip_len == 0: return read, 0` guard in both plus-strand and minus-strand branches.

- **NEW-044 — `int(k)` without try-except on poly(A) model JSON keys** (`core/polya_model.py:130`): replaced dict comprehension with an explicit loop that skips non-integer keys with a DEBUG-level log.

- **NEW-045 — Bare `except` swallows `FileNotFoundError` in BAM MD-tag check** (`core/run_command.py`): `FileNotFoundError` now propagates; other exceptions log a warning and return `False`.

- **NEW-046 — Path traversal via `sample_id` from manifest** (`core/run_command.py:569`): `sample_id` values are now sanitized with `re.sub(r'[^a-zA-Z0-9_\-.]', '_', ...)` before constructing output paths.

- **NEW-047 — No pre-flight file existence check for manifest paths** (`core/run_command.py`): added pre-flight loop in `_run_single_sample()` that raises `FileNotFoundError` immediately for any missing input file.

- **NEW-049 — Zero-position sample missing from pass 2 count matrix** (`core/analyze_command.py`): all `(cluster_id, sample_id)` combinations are now pre-initialized to 0 in the count accumulator before streaming.

- **NEW-050 — Bedgraph written non-atomically** (`core/analyze_command.py:1673`): output is now written to a `.tmp` file then atomically renamed via `os.rename()`.

- **NEW-052 — Comment says "modal" but code uses `np.median()`** (`core/analyze/apa_detection.py:324`): corrected comment to "median" (the implementation is correct for noisy long-read data).

- **NEW-059 — `Path.is_relative_to()` requires Python 3.9 but package declares `python_requires >= 3.8`** (`utils/provenance.py:138`): added `_is_relative_to()` helper using `try/except ValueError` around `path.relative_to()`, compatible with Python 3.8+.

- **NEW-060 — Full-length classifier minus-strand branch uses same field as plus-strand** (`core/classify/full_length_classifier.py:150-153`): both branches were assigning `record.five_prime_corrected`; minus-strand branch now correctly uses `record.reference_end - 1`.

- **NEW-005b — mapPacBio `view_proc` stderr changed to DEVNULL to prevent deadlock** (`core/multi_aligner.py`): changed `stderr=subprocess.PIPE` to `stderr=subprocess.DEVNULL` for `view_proc` to prevent buffer-fill deadlock when no thread drains stderr; `sort_proc` retains `stderr=subprocess.PIPE` for error reporting. `view_proc.kill()` added to timeout handler.

- **`--write-corrected-bam` flag** (`core/bam_processor.py`, `core/correct_command.py`, `cli.py`): new output option that hard-clips every read at its corrected 3' position via CIGAR surgery. Unlike `--output-bam` (which only strips soft-clipped poly-A), this flag modifies aligned CIGAR operations so reads visually end at the corrected position in IGV. Works for both plus strand (clip from right) and minus strand (clip from left, updates `reference_start`). Output is coordinate-sorted and indexed.

- **Bundled genome re-bgzipped** (`data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz`): the bundled yeast reference genome was replaced with a regular gzip file on 2026-04-05, leaving stale `.fai`/`.gzi` index files from the prior bgzip version. pysam opened the file successfully but raised `ValueError` on the first `fasta.fetch()` call, causing `rectify align` consensus selection to crash for all users. Re-compressed with `bgzip` and re-indexed with `samtools faidx`.

- **4-tuple junction unpack** (`core/splice_aware_5prime.py`): two loops in `rescue_3ss_truncation()` unpacked annotated junctions as 3-tuples `(j_chrom, intron_start, intron_end)` but `load_annotated_junctions()` returns 4-tuples including strand. Fixed both sites (different indentation levels) to use `j_entry[0:3]` positional access.

## [2.7.7] - 2026-04-08

### Fixed

- **C2 — 3'SS proximity tolerance off-by-one** (`terminal_exon_refiner.py`): `abs(read_end - pos - 2)` corrected to `abs(read_end - pos)`, eliminating a systematic 2-bp error in soft-clip proximity detection and the matching `distance_from_ss` diagnostic field.

- **C3 — `_apply_profile()` missing keys** (`batch_command.py`): Added `use_scratch`, `streaming`, `analyze_cpus`, `analyze_mem`, and `analyze_time` to the `profile_to_arg` dispatch table and defaults dict. These keys were silently ignored when passed via a profile YAML, so `--streaming`, `--use-scratch`, and analyze-step resource settings had no effect in batch mode.

- **H2 — Minus-strand 5'SS truncation check used wrong coordinate** (`terminal_exon_refiner.py`): Proximity check compared `read_start` (3' end of minus-strand read) to the 5'SS position instead of `read_end - 1` (true 5' end). Corrected to `read_5prime = read_end - 1`.

- **H4 — No dinucleotide validation for GFF-derived 3'SS positions** (`terminal_exon_refiner.py`): Added optional `genome` parameter to `load_splice_sites_from_gff()`; when provided, validates that each 3'SS has the expected AG (plus-strand) or CT (reverse complement of AG on minus-strand) dinucleotide and skips introns with annotation coordinate errors.

- **H5 — `init_provenance()` called without `config`** (`run_command.py`): `_run_junction_aggregation()` now accepts and forwards `config=vars(args)` to `init_provenance()`, so provenance records include the full run configuration.

- **M1 — GC-AG canonical splice motif missing** (`junction_validator.py`): Added `'GC-AG'` to `_CANONICAL_MOTIFS`. This ~0.5–1% of junctions was previously classified as non-canonical, inflating false-positive junction filter rates.

- **M2 — Position shift double-counting** (`processing_stats.py`): `total_position_shifts` counter no longer increments for reads already counted under `atract_ambiguity`, preventing the same correction event from being reported in two categories.

- **M4 — Jensen-Shannon divergence computed as sqrt(JSD)** (`analyze/shift_analysis.py`): Replaced manual JSD implementation with `scipy.spatial.distance.jensenshannon() ** 2`. The scipy function returns the square root of JSD by design; squaring restores the [0, 1] scale expected by downstream shift detection.

- **M6 — Duplicate A-richness threshold** (`false_junction_filter.py`): Local constant `A_RICHNESS_THRESHOLD = 0.7` replaced with the canonical `POLYA_RICHNESS_THRESHOLD` (0.8) imported from config, ensuring consistent filtering behavior.

- **M7 — Missing 0-based coordinate warning** (`analyze/gene_attribution.py`): `build_cds_interval_tree()` now logs a warning when the `start` column minimum is 1, indicating likely un-converted GFF 1-based coordinates.

- **M8 — Minus-strand ambiguity window direction inverted** (`bam_processor.py`): Ambiguity window for minus-strand reads now correctly extends toward lower coordinates (`ambiguity_min = current_position - range`, `ambiguity_max = current_position`) instead of higher.

- **M9 — BAM handle leak in `find_coverage_gaps()`** (`bam_processor.py`): Wrapped `pysam.AlignmentFile` in a `with` context manager; removed manual `.close()` calls and fixed the early-return path for missing chromosomes.

- **M10 — Silent gene attribution exceptions** (`bam_processor.py`): `except Exception: pass` replaced with `except Exception as _e: logger.warning(...)` so attribution failures appear in logs.

- **M11 — Unrecognized SLURM profile keys silently ignored** (`batch_command.py`): `load_slurm_profile()` now checks keys against a known-good set and logs a `WARNING` for any unrecognized keys (typo detection).

- **M12 — `--continue-on-error` not wired to generated SLURM scripts** (`batch_command.py`): Both `SLURM_CORRECT_TEMPLATE` and `SLURM_ANALYZE_TEMPLATE` now use a `{error_handling}` placeholder. When `--continue-on-error` is set, generated scripts use `set -uo pipefail` (no `-e`) so individual task failures do not abort the array job.

- **M13 — No partition check before `--submit`** (`batch_command.py`): `run_batch_command()` now returns an error immediately if `--submit` is requested without `--partition`, preventing silent submission to the cluster default partition.

- **M14 — Cluster-specific example paths in public docstrings** (`batch_command.py`, `cli.py`, `run_command.py`): Replaced cluster-specific profile names in usage examples with the generic `my_cluster.yaml`.

- **M15 — Silent gene attribution drop in APA isoform building** (`analyze/apa_detection.py`): `build_apa_isoforms()` now logs a WARNING when more than 10% of input records lack gene attribution, with the exact count and fraction.

- **M16 — SAM→BAM pipeline return codes unchecked** (`multi_aligner.py`): Both the uLTRA and deSALT `samtools view | samtools sort` pipelines now call `view_proc.wait()` and raise `RuntimeError` on non-zero exit codes from either process.

- **M17 — SGE array indices not normalized to 0-based** (`batch_command.py`): In the portable scheduler abstraction block, `SGE_TASK_ID` is now decremented by 1 (`$((SGE_TASK_ID - 1))`) to match SLURM's 0-based `SLURM_ARRAY_TASK_ID` convention.

- **M18 — Exclusion regions coordinate check emitted only a debug warning** (`exclusion_regions.py`): `filter_dataframe()` now emits an `ERROR`-level log (with fraction of affected rows) when the position column minimum is 1 and >1% of rows are at position 1, strongly indicating un-converted 1-based GFF coordinates.

## [2.7.6] - 2026-04-03

### Fixed

- **Missing `logger` in `false_junction_filter.py`**: `logger.debug()` was called but `import logging` / `logger = logging.getLogger(__name__)` were absent — would raise `NameError` at runtime on any junction motif exception.

- **`NetSeqSignal` not picklable** (`netseq_refiner.py`): Added `__getstate__()` / `__setstate__()` to exclude the `threading.Lock` and open `pyBigWig` handles from pickle state and recreate them on restore. `NetSeqSignal` is now safe to use with `multiprocessing.Pool`.

- **BAM handle leak in `iter_netseq_reads()`** (`netseq_bam_processor.py`): Wrapped read iteration loop in `try/finally` so `bam.close()` is always called, even when `max_reads` limit causes an early `break`.

## [2.7.5] - 2026-04-03

### Added

- **Sequence-based 5' soft-clip rescue** (`consensus.py`): `_rescue_5prime_softclip()` aligns
  soft-clipped bases against annotated upstream exon sequences (edit distance ≤ 20% mismatches),
  replacing the prior length-only scoring. Added `_get_effective_5prime_clip()` and
  `_get_effective_3prime_clip()` for fair mapPacBio scoring (accounts for forced mismatches).
  Expanded `AlignmentInfo` dataclass with `five_prime_softclip_seq`, `effective_five_prime_clip`,
  `effective_three_prime_clip` fields.

- **Junction filtering utilities** (`aggregate/junctions.py`): `filter_junctions()` removes
  junctions below read-count and canonical-motif thresholds; `resolve_homopolymer_ambiguity()`
  assigns ambiguous junctions flanking homopolymers to the best-supported position.

- **Genome pickle cache** (`utils/genome.py`): `load_genome()` writes a `.pkl` sidecar on first
  load; subsequent loads use the cache, saving 10–120s per sample.

### Fixed

- **Bug 27 — `detect_partial_junction_crossings()` TypeError on SEQ=* reads**
  (`terminal_exon_refiner.py`): `len(clip_seq)` crashed with `TypeError` when
  `five_prime_clip['sequence']` is `None`. Added `if clip_seq is None: continue` guard.

- **Bug 26 — `generate_bedgraphs()` KeyError** (`analyze_command.py`): Fallback column
  name `'position'` (non-existent) replaced with `'corrected_position'`.

- **Bug 24 — `_run_junction_aggregation()` KeyError** (`run_command.py`): Dict key
  `partial_results['summary']` replaced with `partial_results['stats']` (correct key from
  `detect_partial_junction_crossings()`). Also fixed `.suffix` → `.suffixes` check for
  `.gff.gz` annotation detection.

- **Bug 23 — `deconvolve_region()` silent bad output** (`netseq_deconvolution.py`): Inverted
  `region_start ≥ region_end` coordinates now raise `ValueError` instead of returning a
  zero-length result silently.

- **Bug 22 — `get_netseq_3prime_position()` crash on unmapped reads** (`netseq_bam_processor.py`):
  Added `if read.is_unmapped or read.reference_end is None: raise ValueError(...)` guard.

- **Bugs 20–21 — Silent fallbacks in `_stitch_group()`** (`mpb_split_reads.py`): Added
  `logger.debug()` calls for all-unmapped reads and cross-chromosome/cross-strand fallbacks
  so these events appear in logs at DEBUG level.

- **Bug 19 — BAM file descriptor leak** (`terminal_exon_refiner.py`):
  `refine_terminal_exons()` converted from `bam.open()/bam.close()` to
  `with pysam.AlignmentFile(...) as bam:` context manager.

- **Bugs 17+11 — OOM and chrom normalization in `junction_validator.py`**:
  `filter_cross_sample_junctions()` replaced `pd.concat` + `itertuples()` (OOM on large
  datasets) with streaming dict accumulation. `read.reference_name` now passes through
  `standardize_chrom_name()` before junction key construction.

- **Bug 15 — Negative `ambiguity_min` for reads near chromosome start** (`bam_processor.py`):
  `ambiguity_min` clipped to `max(0, ...)` to prevent negative coordinates.

- **Bug 7 (partial) — Missing timeouts on minimap2/mapPacBio/gapmm2** (`multi_aligner.py`):
  `ALIGNER_TIMEOUT` (7200s) now applied to all five aligners (was only applied to uLTRA
  and deSALT). Added `try/except TimeoutExpired` with process kill for minimap2 pipe.

### Ported from rectify-beta

- `docs/IMPLEMENTATION_SUMMARY.md` — technical architecture and validation suite
- `docs/ALIGNER_RECOMMENDATIONS.md` — guidance on expanding aligner panel

## [2.6.0] - 2026-04-03

### Fixed
- **mapPacBio AssertionError on RNA004/Dorado basecalled reads** (`mpb_split_reads.py`):
  Dorado basecaller using the RNA004 model outputs `U` (uracil) instead of `T` in FASTQ
  sequences. BBMap v39.26 does not handle `U` and throws a Java `AssertionError` in
  `genMatchStringForSite()` on any such read, crashing the entire alignment job.
  Fix: normalize `U→T` in `split_long_reads()` before writing the split FASTQ, so all
  reads reaching mapPacBio are in standard DNA notation regardless of basecaller.

- **Missing `_GenomeDictReference` adapter in `false_junction_filter`** (`false_junction_filter.py`):
  `bam_processor.py` called `_fjf._GenomeDictReference(genome)` which existed in beta
  but was never ported to production, causing a crash during the correction step for
  every sample: `AttributeError: module 'rectify.core.splice.false_junction_filter' has no
  attribute '_GenomeDictReference'`. Ported the thin adapter class (wraps a genome dict
  to present a `pysam.FastaFile`-compatible interface) from beta.

## [2.5.0] - 2026-04-02

### Added
- **uLTRA aligner support** (`multi_aligner.py`, `align_command.py`, `run_command.py`, `cli.py`):
  Annotation-guided collinear chaining aligner optimized for small exons (11–20 nt).
  Handles GFF→GTF conversion automatically (looks for sibling `.gtf`, or converts via
  `gffread`). Decompresses gzip genomes to a temp file (uLTRA cannot read gzip).
  Enabled via `--junction-aligners uLTRA` (opt-in, requires `--annotation`).

- **deSALT aligner support** (`multi_aligner.py`, `align_command.py`, `run_command.py`, `cli.py`):
  De Bruijn graph mapper. Handles deSALT's known output-duplication bug by deduplicating
  on (read_name, flag, chrom, pos, cigar) via `_dedup_desalt_bam()`. Strips
  `LD_LIBRARY_PATH` from deSALT's environment (prevents library conflicts). Resolves
  the deSALT index directory (`desalt_index/`) adjacent to the genome FASTA.
  Enabled via `--junction-aligners deSALT` (opt-in).

- **5-aligner consensus pipeline**: `--junction-aligners uLTRA deSALT` adds both
  junction-mode aligners to the 3-aligner default (minimap2 + mapPacBio + gapmm2),
  producing a 5-aligner rectified BAM.

- **`run-all` CLI exposes junction aligner flags**: `--junction-aligners`, `--ultra-path`,
  `--desalt-path`, `--chimeric-consensus` added to the `run-all` subcommand.

### Changed
- **Two-phase scheduler updated for deSALT fork-safety**: deSALT crashes with
  "double free or corruption" when forked from a multithreaded Python process.
  The scheduler now runs deSALT sequentially after the parallel pool completes
  (`parallel_batch = [a for a in remaining if a != 'deSALT']`).

- **`_run_alignment()` in `run_command.py`** now accepts and threads through
  `junction_aligners`, `chimeric_consensus`, `ultra_path`, `desalt_path`.

## [2.4.0] - 2026-04-01

### Added
- **Chimeric consensus selection** (`chimeric_consensus.py`):
  Divides a read into segments at "sync points" where all aligners agree, then
  independently selects the best aligner per segment and assembles a chimeric CIGAR.
  Wired into the pipeline behind `--chimeric-consensus` flag (default off, experimental —
  requires further validation before enabling by default).

- **Test suite ported from beta**: Unit and integration tests for consensus selection,
  chimeric consensus, junction validation, false junction filtering, and aligner wrappers.

### Changed
- **Visualization module updated**: Additional plot types and improved styling.

## [2.3.0] - 2026-03-30

### Added
- **Module 2F: 3'SS truncation rescue** (`bam_processor.py`):
  Corrects 5' positions for reads truncated or soft-clipped at the exon 2 / 3' splice
  site boundary. Runs post-consensus on reads with real (non-poly-A-artifact) junctions.

- **Memory-efficient streaming pipeline for multi-sample analysis** (`analyze_command.py`):
  Two-pass streaming over per-sample corrected TSVs using a manifest; never loads more
  than one sample into RAM. Position index (`corrected_3ends_index.bed.gz`) support for
  ~300× faster Pass 1 and Pass 2. Peak RAM is O(clusters × samples) regardless of depth.

- **Pre-consensus 3' A-tract scoring**: A-tract ambiguity is now scored before consensus
  selection so it can inform aligner selection, not just final position correction.

- **Position index writing**: `rectify correct` now writes `corrected_3ends_index.bed.gz`
  alongside `corrected_3ends.tsv` for use by manifest-mode analysis.

### Fixed
- **Spike-in filter false positives** (`correct_command.py`): Fixed over-aggressive
  filtering that was removing valid reads near spike-in sequences.

- **Bug 14: Silent exception swallow in `false_junction_filter.py`**: Exceptions in
  junction motif analysis were silently swallowed (bare `except: pass`), hiding real
  errors. Now logs at DEBUG level with the exception message.

- **Bug 16: CIGAR/seq mismatch for hard-clipped donors** (`consensus.py`):
  `_restore_sequence_from_aligner_reads()` copied donor sequence without checking length
  against CIGAR-implied query length. Hard-clipped donors (deSALT) caused
  `samtools sort` crash: "CIGAR and query sequence lengths differ". Fixed with
  `_cigar_query_length()` helper; donors whose sequence length doesn't match are skipped.

- **`ref_end is None` guards in `bam_processor.py`**: Unmapped reads (no `reference_end`)
  caused `TypeError: unsupported operand type(s) for -: 'NoneType' and 'int'` in
  `get_read_3prime_position()` and `get_read_5prime_position()`. Now returns `None`
  for unmapped reads.

- **`itertuples()` OOM in `junction_validator.py`**: `filter_cross_sample_junctions()`
  used `itertuples()` on a large aggregated DataFrame, causing OOM on large datasets.
  Replaced with numpy array iteration.

### Changed
- **Two-phase aligner scheduler**: mapPacBio (~10× slower) runs alone with all threads
  (phase 1); remaining aligners split threads equally in parallel (phase 2). Replaces
  the old `_THREAD_WEIGHTS` dict approach.

- **Per-aligner BAMs name-sorted directly**: Skips intermediate coordinate sort and
  index, reducing disk I/O during the alignment phase.

## [2.2.0] - 2026-03-28

### Added
- **`run-all` command** (`run_command.py`): Complete end-to-end pipeline
  (align → correct → analyze) for single samples and multi-sample manifests.
  Single-sample: steps 0–3 (alignment, correction, analysis, junction aggregation).
  Multi-sample: parallel per-sample correction + combined DESeq2/GO/motif analysis.

- **Scratch staging** (`slurm.py`): SLURM array jobs automatically stage I/O through
  `$SCRATCH` (75 GB/s) to avoid Oak NFS contention. Final outputs sync'd back to Oak
  via `rsync`. `make_job_scratch_dir()` / `sync_to_oak()` utilities.

- **Infrastructure ported from beta**:
  - `junction_validator.py`: Cross-sample junction validation with configurable thresholds
  - `slurm.py`: CPU detection, thread limits, scratch staging utilities
  - HPC SLURM profile (`hpc_cpu.yaml`) with `use_scratch: true` and
    `streaming: true` defaults

## [2.1.0] - 2026-03-16

### Added
- **train-polya command**: Train poly(A) tail models from control sites with 0 downstream A's
  - Learns A-richness distribution, length distribution, position profiles
  - Outputs JSON model file for use with `correct` command
  - Supports minimum read count thresholds and per-site statistics

- **validate command**: Comprehensive validation against ground truth
  - Multiple ground truth sources: NET-seq BigWig, gene annotations (GTF/GFF), known positions TSV
  - All sources optional but at least one required
  - Calculates accuracy at multiple tolerances (exact, 1bp, 2bp, 5bp)
  - Per-correction-type breakdown and confidence calibration metrics

- **PolyAModel class**: Dataclass-based poly(A) tail model with JSON serialization
  - Training statistics tracking
  - Position-dependent A-frequency profiles
  - Score sequences against learned parameters

- **LRU cache in NetseqLoader**: Bounded cache (default 10,000 entries) prevents memory issues

- **GitHub Actions CI/CD**: Automated testing across Python 3.8-3.12, linting, and package building

### Changed
- **A-tract ambiguity detection**: Now correctly counts downstream A's in window (not just tract length)
  - Added `_get_sequence()` helper for chromosome name conversion
  - Added `_count_downstream_as()` for proper window-based counting
  - Returns both `downstream_a_count` and `expected_shift` fields

- **Performance optimizations**:
  - `find_polya_boundary()`: O(n) sliding window instead of recalculating per position
  - `complement()`: Uses `str.maketrans`/`translate` instead of per-character dict lookup
  - `generate_summary_report()`: Single-pass accumulation instead of 5 iterations
  - `find_coverage_gaps()`: Avoid O(n) `.index()` call in chromosome lookup

- **Consolidated `find_atract_boundaries`**: Removed duplicate implementation from `genome.py`

### Fixed
- Version mismatch between config.py and pyproject.toml
- Test expectations aligned with actual algorithm behavior:
  - Insertions do NOT change reference coordinates (counted for QC only)
  - Soft-clipped poly(A) tails don't shift genomic position

### Dependencies
- Added `tqdm>=4.60.0` for progress bars

## [2.0.0] - 2026-03-09

### Added
- Initial unified RECTIFY framework
- A-tract ambiguity detection module
- Poly(A) tail trimming with adapter detection
- Indel artifact correction
- NET-seq refinement for position validation
- AG mispriming screening
- SLURM-aware thread management
- Parallel processing with coverage-gap-based region splitting

### Core Modules
- `atract_detector.py`: Universal A-tract ambiguity detection
- `polya_trimmer.py`: Poly(A) tail detection and scoring
- `indel_corrector.py`: Small indel artifact correction
- `netseq_refiner.py`: NET-seq-based position refinement
- `ag_mispriming.py`: AG-rich region flagging
- `bam_processor.py`: Parallel BAM processing pipeline

### Commands
- `rectify correct`: Main correction pipeline
- `rectify train-polya`: (stub) Train poly(A) model
- `rectify validate`: (stub) Validate corrections
