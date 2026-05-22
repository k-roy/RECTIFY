# Audit history & pre-public version log

Two rounds of systematic codebase audits were completed 2026-04-08. All
findings are tracked in `dev/BUGS_TO_FIX.md`. **No open bugs remaining.**

---

## Documentation audit — 2026-05-21

**Baseline:** HEAD commit `6813bd1` (branch `drs-validation-rebuild`).
**Method:** Cross-checked all documentation against `python -m rectify <cmd> --help` output
(stashed WIP before each invocation), file existence checks, and grep verification.
63 files modified / 6078 lines inserted in WIP not audited (deferred to WIP resolver).

### HIGH severity

**H1 — `docs/user_guide/commands/correct.md`: `--ag-threshold` range and default wrong**

- Doc states: "AG-richness threshold (0–1); default 0.65"
- Actual CLI: "AG-richness weighted score threshold for mispriming flagging (0.0–34.5; default 17.0)"
- The threshold is a weighted composite score, not a simple [0,1] fraction.
  Users reading the docs who set `--ag-threshold 0.5` would be supplying a value
  far below the noise floor and flag almost nothing.

**H2 — `docs/user_guide/commands/correct.md`: Entire "Junction refinement" section absent**

The following eight flags exist in `rectify correct` (Module 2H, post-consensus N-op
refinement) and are not documented in `correct.md`:

| Flag | Default |
|------|---------|
| `--aligner-bams BAM` | (repeatable; enables Module 2H) |
| `--junction-hp-pen FLOAT` | 0.25 |
| `--junction-search-radius BP` | 5000 |
| `--junction-window BP` | 15 |
| `--junction-max-slide BP` | 10 |
| `--junction-max-boundary-shift BP` | 50 |
| `--junction-penalty-table PATH` | (empirical HP costs) |
| `--str-penalty-table PATH` | (empirical STR costs) |

Module 2H is the primary mechanism for correcting off-by-a-few-bp junction boundaries
and is the only user-configurable path for empirical penalty tables. Having no `correct.md`
entry makes these flags effectively invisible outside `--help`.

**H3 — `docs/user_guide/commands/correct.md`: "Resume / sidecar" and other operational flags absent**

The following flags exist in `rectify correct` but have no `correct.md` entry:

| Flag group | Flags |
|------------|-------|
| Checkpoint/resume | `--checkpoint-dir DIR`, `--tmp-dir DIR`, `--variant-scan-cache PKL`, `--junction-pool-cache PKL` |
| Force/sidecar | `--force-all`, `--force-stage NAME[,NAME...]`, `--accept-prior-provenance`, `--dry-run-resume` |
| Output | `--emit-merged-tsv`, `--legacy-single-threaded` (DEPRECATED) |

These are operational flags most users will need when running `rectify correct` in HPC batch
pipelines. `--checkpoint-dir` in particular is critical for long streaming runs (>10 GB BAMs).

### MEDIUM severity

**M1 — `docs/user_guide/commands/run.md`: ~11 new flags absent**

The following flags exist in `rectify run-all` but are not in `run.md`:

| Flag | Description |
|------|-------------|
| `--aligner-concurrency {auto,1,N}` | Parallel-aligner correction worker budget |
| `--junction-overhang-table PATH` | Empirical overhang filter (v3.3.0 feature) |
| `--write-softclip-bam` | Write soft-clipped corrected BAM alongside primary |
| `--write-polya-bam` | Write poly(A)-restored BAM for IGV validation |
| `--trust-existing-bams` | Reuse pre-existing BAMs without sidecar check |
| `--bam-dir DIR` | Redirect alignment BAM output directory |
| `--keep-aligner-bams` | Retain per-aligner BAMs after consensus selection |
| `--scratch-dir DIR` | HPC scratch for intermediate BAM I/O |
| `--continue-on-error` | Continue remaining samples if one fails |
| `--junction-penalty-table PATH` | Pass-through to `rectify correct` |
| `--str-penalty-table PATH` | Pass-through to `rectify correct` |

**M2 — `docs/user_guide/commands/correct.md` L72: deprecated alias `--no-polya-sequenced` incorrectly attributed**

`correct.md` states: *"deprecated aliases retained for backwards compatibility:
`--polya-sequenced` / `--no-polya-sequenced`"*

Verified in source: `--polya-sequenced` exists in `correct_command.py` (line 1370).
`--no-polya-sequenced` does **not** exist in `correct_command.py` — it is in
`run_command.py` only (line 413). The `correct.md` entry overstates the alias
coverage for this subcommand.

**M3 — `docs/ALIGNER_RECOMMENDATIONS.md`: Dead internal link (×2)**

Lines 72 and 172 both reference `docs/desalt_crash_investigation_handoff.md`.
That file does not exist. Verified with `ls docs/ | grep desalt` → no output.
The link should either be removed, replaced with a pointer to the relevant
`AGENT_FIXES.md` entry, or the handoff doc should be created.

**M4 — `docs/user_guide/commands/consensus.md`: `--no-bedgraph` flag absent**

`rectify consensus --help` includes `--no-bedgraph` (skips bedGraph/bigWig
generation post-consensus). Not documented in `consensus.md`.

### LOW severity

**L1 — `docs/user_guide/commands/run.md` L147: `--threads` default says "4"**

- Doc: "default 4"
- Actual CLI: "0 (auto-detect from SLURM_CPUS_PER_TASK or CPU count)"

**L2 — Memory index `reference_rectify_5aligner_status.md` is stale**

Entry says `docs/ALIGNER_RECOMMENDATIONS.md` "is stale (says SIGSEGV — needs update)".
The file was updated 2026-05-19 and now correctly documents the SIGSEGV as a resolved
historical issue with the vendored binary. Memory entry should be updated to reflect
that the doc is current.

### Confirmed current

The following were checked and found accurate against HEAD:

| Document | Status |
|----------|--------|
| `docs/ALIGNER_RECOMMENDATIONS.md` | Current (updated 2026-05-19; SIGSEGV documented as resolved) |
| `docs/algorithms/overview.md` | Current; Module 2H described; HP penalty calibration accurate |
| `docs/ARCHITECTURE.md` | Current; DRS Steps 0–5, cDNA pipeline, consensus scoring all match code |
| `docs/quickstart.md` | Current; all example commands use valid flags |
| `docs/index.md` | Current; `mapPacBio (BBMap)` label is correct (`mapPacBio.sh` is from BBMap/BBTools) |
| `docs/user_guide/commands/split.md` | Current; `--junction-overhang-table` documented |
| `docs/user_guide/commands/consensus.md` | Mostly current; `--read-num-sidecar` and `--chimeric` documented |
| `CHANGELOG.md` | Current for 0.9.0 public release |
| `docs/user_guide/output_formats.md` | Current; columns match `corrected_3ends.tsv` schema |
| `docs/user_guide/commands/run.md` | Mostly current; pipeline steps, DRS workflow, manifest format accurate |

### Recommended actions (priority order)

1. **Fix `--ag-threshold` range and default in `correct.md`** (H1) — one-line edit, high
   user-facing impact.
2. **Add "Junction refinement" section to `correct.md`** (H2) — Module 2H is production-grade
   for yeast DRS runs; currently invisible to docs readers.
3. **Add "Checkpoint / resume / sidecar" section to `correct.md`** (H3) — required knowledge
   for any HPC deployment.
4. **Update `run.md` with 11 missing flags** (M1) — `--aligner-concurrency` and
   `--junction-overhang-table` in particular affect production run quality.
5. **Fix deprecated-alias attribution in `correct.md`** (M2) — remove `--no-polya-sequenced`
   from the `correct` entry.
6. **Resolve dead link to `desalt_crash_investigation_handoff.md`** (M3) — either create
   the stub file or change the reference.
7. **Add `--no-bedgraph` to `consensus.md`** (M4) — minor but complete.
8. **Fix `--threads` default in `run.md`** (L1) — trivial.
9. **Update memory entry `reference_rectify_5aligner_status.md`** (L2) — stale agent context.

---

This file collects per-version notes for changes prior to the 0.9.0
public release. New entries should go in `CHANGELOG.md`; this file is
for historical reference only.

---

**v2.7.8 (2026-04-09):** Fixed 35 Round 2 findings including all
confirmed CRITICALs and HIGHs. Most impactful fixes:
- 6 minus-strand correctness bugs (junction validation donor/acceptor
  swap, junction aggregate strand handling, NET-seq double-subtraction,
  chimeric consensus CIGAR anchoring, 5' splice rescue window,
  full-length classifier).
- Shell injection in generated SLURM scripts (`shlex.quote()` added).
- `streaming: true` / `use_scratch: true` silently ignored by profile
  apply.
- `NCBI_TO_CHR` reversal overwrote `chrI` with `I`.
- `mapPacBio` timeout/deadlock/return-code issues.

**v2.7.9 (2026-04-09):** Fixed 4 additional findings from validation
audit (N-op absorption in CIGAR surgery, `five_prime_rescued` TSV gap,
chimeric walkback exemption, NET-seq flag guard).

**v2.8.0 (2026-04-11):** Cat3 5' junction rescue now uses semi-global NW
local alignment (`local_aligner.py`) for the exon CIGAR instead of a
flat `nM` block. New `five_prime_exon_cigar` column in
`corrected_reads.tsv`. Chimeric stitch gaps ≤ 10 bp use CIGAR op `D`
instead of `N`. Local aligner upgraded to affine gap scoring (Gotoh
1982, gap_open=-4, gap_extend=-1) to prevent staircase artifacts where
isolated 1-base deletions outscore a single consolidated deletion.

**v2.8.1 (2026-04-11):** New Module 2G: standalone soft-clip rescue at
homopolymer boundaries (`rescue_softclip_at_homopolymer`) runs for all
protocols (not just `--dT-primed-cDNA`). Renamed `--polya-sequenced` →
`--dT-primed-cDNA` (deprecated alias preserved). Cat2 validation reads
replaced with T-tract homopolymer examples. Direct RNA / dT-primed cDNA
protocol distinction clarified.

**v0.9.0-dev (2026-04-11):** All 4 previously open items resolved:
- Bug 37 (HIGH): `tests/test_terminal_exon_refiner.py` added (51 tests,
  real wt_by4742_rep1 data).
- Bug 38 (HIGH): `tests/test_consensus_selection.py` added (40 tests,
  real wt_by4742_rep1 data).
- Bug 41 (MEDIUM): `--polya-model` wired through correction pipeline;
  `pt_tag`/`polya_score`/`polya_source` columns in TSV; `rectify
  tag-polya` subcommand; unaligned dorado BAM auto-detection.
- Bug 55 (MEDIUM): `--max-cluster-radius`, `--min-peak-sep`,
  `--min-cluster-samples` added to `analyze` subparser.

**v3.1.7 (2026-04-21):** `refine_read_junctions` — bilateral t2, no-
candidate-guard policy, canonical HP prior:

- **No candidate guards**: all junctions in the candidate pool are
  scored; non-canonical, non-annotated (novel) alternatives are no
  longer filtered out before scoring. The previous guard (`if is_alt==1
  and tier>=4 and is_novel==1: continue`) was removed because it would
  silently discard reads that genuinely belong at non-canonical
  junctions — e.g. when many reads from the same splice isoform all
  score perfectly at a novel non-canonical site. Annotation and
  canonical tier remain as TIE-BREAKERS only, never as gates. This
  policy is PERMANENT and must not be re-introduced.

- **Bilateral t2 scoring**: `_score_junction` now uses `score(k) = t1(k)
  + t2(k)` where:
  - `t1(k) = _score_hp_anchored(rescue[k:],       g[je      : je+buf ])` — intron_end-anchored
  - `t2(k) = _score_hp_anchored(rescue[:k][::-1], g[je-buf  : je     ][::-1])` — intron_end-proximal
  Both anchored to `intron_end` (je). t2 prevents degenerate k=L-1
  coincidental single-base matches from scoring 0 at the wrong junction
  without any pre-filtering.

- **Canonical HP prior** (replaces `int(score)` floor binning): when the
  current N-op is non-canonical (`tier_beats_alt=True`), canonical-tier
  alternatives (tier < 4) receive a 0.5 edit-distance discount
  (`_CANONICAL_HP_PRIOR = 0.5`). This equals one Nanopore HP deletion
  equivalent — the expected noise floor — and ensures canonical
  junctions win within the noise floor regardless of which penalty
  table is in use. Non-canonical alternatives must exceed the canonical
  score by >0.5 to win. The old `int(score)` floor was fragile: it
  worked only when both scores happened to fall in the same [n,n+1)
  integer bucket, which fails when the empirical del_cost(1)≈0.43 and
  scores straddle an integer boundary.

- **Impact**: all 146 tests pass with both default and empirical
  penalty table (`penalty_scores.tsv`); all 4 cat9 reads correctly
  refined with `--junction-penalty-table`.

**v3.3.0 (2026-05-02):** Chimeric-junction overhang filter
(`--junction-overhang-table`) + two critical `corrected_consensus.py`
correctness fixes:

- **New module: `calibrate_junction_overhang.py`** — empirically
  calibrates the minimum-overhang threshold per intron-size bin from
  concordant multi-aligner reads (≥4 aligners agree on the junction).
  Produces `junction_overhang_table.tsv`. Uses isotonic regression
  (PAV) to enforce monotone non-decreasing thresholds. Example
  production table in
  `common/scripts/nanopore/overhang_table_20260502/`.
- **New chimera filter in `merge_corrected_tsvs`** —
  `_add_chimera_flag(rep_df, overhang_table, junction_stats)` sets
  `_chimera_ok=1` for aligner rows whose junction has insufficient
  flanking overhang for its intron size. Chimeric aligners are sorted
  AFTER non-chimeric ones; if all aligners are chimeric, HP edit
  distance resolves the tie (graceful degradation). Short-intron
  relaxation (< 500 bp, ≥5 cross-read support, max_overhang ≥ 10 bp)
  allows 1 bp minimum overhang for well-supported junctions. Wired
  into `rectify run-all` via `--junction-overhang-table`.
- **FIX: `index_col=False` in `_load_tsv`** — minimap2/gapmm2/uLTRA
  corrected TSVs are written with `index=True` (32 header columns, 33
  data columns). Without `index_col=False`, pandas auto-uses the first
  data column (UUID) as the row index and shifts ALL column names one
  right: `read_id` → chromosome, `chrom` → strand, etc. With
  `index_col=False`, column mapping is correct. The root cause: these
  three aligners write TSVs via `pd.DataFrame.to_csv(index=True)` with
  an unnamed index, which pandas reads back as a header/data length
  mismatch. deSALT and mapPacBio write with `index=False` (or with a
  named index column) and are unaffected. **Impact**: in any merge that
  ran before this fix, gapmm2/minimap2/uLTRA read_ids were chromosome
  names (~17 unique values), so those aligners NEVER properly competed
  with mapPacBio/deSALT in winner selection. The original
  chunk003_correct_test merged TSV (pre-fix) was generated with this
  bug active.
- **FIX: `_normalize_read_id` handles both `_pt:i:N` forms** — the
  mapPacBio corrected TSV stores read names as `UUID_pt:i:N`
  (underscore-separated, from `samtools sort`), but the previous
  `_normalize_read_id` only stripped the space-separated form (`UUID
  pt:i:N`). After the fix, both space and underscore forms are
  stripped to bare UUID before the merge join. Without this fix,
  mapPacBio rows would have `UUID_pt:i:N` read_ids that didn't match
  any other aligner's `UUID` read_ids, causing mapPacBio to win 100%
  of reads from its own isolated groups (no competition) rather than
  competing on HP edit distance.

**v3.2.5 (2026-04-24):** Validation Cat1/Cat2 reads replaced with
DRS-trimmed examples showing actual correction artifacts:
- **Cat1** (4 reads): replaced old chunked-pipeline reads
  (9747f421/0821dc9e/058438bb/4e5e6eee, chrVIII/chrIX) with DRS
  minimap2 reads (0cb5a111/a146838d/77b392d9/34ba198b,
  chrXIV/chrI/chrII/chrXII) that genuinely show
  `atract_ambiguity,indel_correction,polya_walkback`. Shifts: −16, −4,
  +3, +3 bp.
- **Cat2** (3 reads replaced, 9dbd37bf kept): replaced old chunked
  reads (6ad42e7a/5bd31a5e/8786c81f, chrI) with DRS minimap2 reads
  (61b0c014/88953e9c/b313b50d, chrI/chrVI/chrV) that genuinely show
  `softclip_rescue`. Shifts: +11, +4, −8 bp; sc_rescued_seq:
  TT/AGT/CCTAG.
- **Rebuild method**: `replace_cat1_reads.py` + `replace_cat2_reads.py`
  sourced from `wt_by4742_rep1_drs_trim_20260417` merged BAMs.
  In-place pysam BAM rebuild corrupted Cat4/Cat6/Cat7/Cat9 reads (N
  ops lost, CIGARs changed). Fixed by `rebuild_from_committed.py`
  which merged Cat1/Cat2 from working tree with Cat3–9 from `git HEAD`.
  Same fix applied to all 5 aligner BAMs via `rebuild_aligner_bams.py`.
  See CAUTION below.
- **CAUTION — in-place pysam rebuild**: the `read→.tmp→pysam.sort→overwrite`
  pattern silently converts `=`/`X` CIGAR ops to `M` and repositions
  reads when the source BAM uses mixed CIGAR op styles. Always verify
  N-ops and read IDs for ALL categories after any BAM rebuild, not
  just the target category.
- Test updates: Cat1 expected positions changed (chrVIII/chrIX→chrXIV/
  chrI/chrII/chrXII); Cat2 expected positions changed
  (chrI→chrI/chrVI/chrV). All 708 tests pass.

**v3.2.3 (2026-04-24):** Validation read sequences synced to DRS-trimmed
run for Cat5-9:
- **Root cause (v3.1.8 partial update)**: `update_validation_drs.py`'s
  "update in place" path updated CIGAR/N-op boundaries from the DRS-
  trimmed mapPacBio run but left read *sequences* from the old
  chunked-consensus BAM. This caused visual discrepancies in IGV
  between `rectified_corrected_3end.bam` (corrected from chunked
  sequences) and `validation_reads.mapPacBio.bam` (aligner BAM with
  DRS sequences). cat5_minus_1 retained `SEQ='*'` from the gapmm2 PAF
  issue (fix was never committed).
- **Fix**: targeted 9-read update preserving all XV/XG/XU tags and
  MD/cs/NM stripped:
  - cat5_minus_1: filled 561-base chimeric sequence from chunked
    consensus (chimeric CIGAR qlen=561; cannot use DRS mapPacBio
    seqlen=538 without breaking the chimeric structure).
  - cat6_plus_1, cat6_minus_2: seq+CIGAR+start from DRS mapPacBio
    merged BAM.
  - cat7_minus_1, cat7_plus_1, cat7_plus_2: same.
  - cat9_plus_1, cat9_plus_2, cat9_minus_1: from DRS per-chunk
    mapPacBio BAMs.
  - **Cat1–4, Cat8 intentionally kept on chunked sequences**: a full
    swap to DRS-trimmed caused Cat4 reads to map to wrong genomic
    regions (minimap2 DRS alignment chose a different locus for
    cat4_plus_1: 19589 vs expected 22072), and Cat1–2 reads no longer
    exhibited the indel/soft-clip artifacts required for their
    validation scenarios.
- **Test updates**: Cat7 EXPECTED_JUNCTIONS shifted 3–8 bp (cat7_plus_1:
  138856→138864, cat7_plus_2: 595736→595739); Cat9 RAW_JUNCTIONS shifted
  1–3 bp (cat9_plus_1 start: 555824→555825; cat9_plus_2 end:
  439324→439321) — same GT-AG annotated junctions, same correction
  outcomes.
- All 110 tests pass (3 skipped).

**v3.2.2 (2026-04-22):** Validation dataset fully certified for production
— XV/XG tag fix + aligner BAM N-op correction:
- **Root cause**: `rebuild_aligner_bams_cat679.py` (v3.2.0) sourced
  replacement read alignments from `validation_reads.bam` before that
  BAM's replacement reads had proper N-ops. Additionally, the 3
  replacement reads (f8050895, 7d5e8dc2, 72557a9a) were missing XV and
  XG auxiliary tags because they were inserted with full UUIDs while
  the source (stale BAM) had 8-char prefix names — the tag-copy logic
  couldn't match them.
- **XV/XG tag fix**: both tags are now present on all 36 reads in
  `validation_reads.bam`. XV (read label, e.g. `cat6_plus_2`) is
  required by `load_reads()` in the test fixture; XG (category, e.g.
  `cat6_chimeric`) drives category-level test assertions.
- **N-op fix for aligner BAMs**:
  `dev_runs/wt_by4742_rep1_drs_trim_20260417/fix_replacement_read_cigars.py`
  re-sourced correct alignments directly from dev run BAMs —
  `wt_by4742_rep1_drs_trim_20260417` (mapPacBio) and
  `wt_by4742_rep1_chunked_20260412` (deSALT/gapmm2/minimap2/uLTRA) —
  preserving full UUID QNAMEs from aligner BAMs and stripping
  `_pt:i:N` suffixes when matching. Correct N-op intervals: f8050895
  all-aligners (45644–45977); 7d5e8dc2 mapPacBio-only (60193–60697);
  72557a9a mapPacBio-only (104435–104495).
- **Full pipeline re-run**: `corrected_reads.tsv` (36 rows: 32 high + 4
  chimeric), `rectified_corrected_3end.bam`,
  `rectified_pA_tail_trimmed.bam`, `rectified_pA_tail_soft_clipped.bam`
  all regenerated from the fixed validation_reads.bam.
- **Independent certification**: 4 independent agent reviews
  (sequence-level genomic context, structural checks, category-specific
  assertions, cross-agent consensus) all confirmed 36/36 PASS. Cat1/Cat2
  shifts match VALIDATION_READS.md exactly; Cat3 five_prime_rescued=1;
  Cat4 n_junctions=1; Cat5 correction=none+chimeric; Cat6/7 XU=1; Cat9
  junction refinement corrects all 4 reads.
- All 708 tests pass.

**v3.2.1 (2026-04-22):** `rescue_3ss_truncation` — soft-clip exon CIGAR
body-borrowing fix:
- `_get_intronic_query_bases` includes one extra body base when the
  alignment ends exactly at `intron_start` (boundary straddle: last M
  op's final ref base == intron_start). For soft-clip rescues this
  produced a CIGAR with one more query-consuming op than
  `five_prime_soft_clip_length` (e.g. `14M1D9M` for 79f61403,
  query_span=23 vs soft_clip=22). The `bam_writer` guard then fell
  back to a flat `22M` block, applying wrong exon geometry in IGV.
- Fix: for `rescue_type_candidate == 'softclip'`, use `rescue_seq`
  (exactly the soft-clip, already truncated to `five_clip` bases) as
  `_align_seq` instead of `_intronic_seq`. For `mpb_mismatch` rescues
  `_intronic_seq` is still used (it correctly excludes exon-2 body
  bases beyond the intron boundary).
- Impact: 79f61403 (`cat3_plus_2`, YAL003W) now correctly produces
  `14M1D7M1D1M` (query_span=22) instead of flat `22M`; exon geometry
  in IGV matches mapPacBio alignment. 28ea9379 (`cat3_minus_2`,
  YBR062C) was unaffected (its `_intronic_seq` was already empty in the
  current validation data).
- All 708 tests pass.

**v3.2.0 (2026-04-22):** Validation aligner BAMs rebuilt to match
`validation_reads.bam` after DRS rebuild (v3.1.8):
- The DRS rebuild (v3.1.8) replaced 3 old Cat6/Cat7 reads in
  `validation_reads.bam` but did NOT update the individual aligner
  BAMs (`aligners/validation_reads.*.bam`). All 5 aligner BAMs still
  contained the old reads (ba761413, 64f4da08, 5c59f0bc) and lacked
  the replacements (f8050895, 7d5e8dc2, 72557a9a).
- `dev_runs/wt_by4742_rep1_drs_trim_20260417/rebuild_aligner_bams_cat679.py`
  removes the 3 old reads from each aligner BAM and adds the 3 new
  reads. NOTE: this script sourced alignments from
  `validation_reads.bam` which at the time also lacked N-ops for the
  replacement reads — the N-ops were fixed in a subsequent pass by
  `fix_replacement_read_cigars.py` (see v3.2.2).
- This inconsistency did not break any tests (Cat6 tests use only the
  basic `corrected` fixture; Cat9 tests use `--aligner-bams` but Cat9
  reads were not changed). It was confusing because ba761413 appeared
  in IGV from aligner BAMs but not from the main rectify list.
- All 708 tests pass.

**v3.1.9 (2026-04-22):** `rectified_pA_tail_soft_clipped.bam` added as
default Step 4 output; sort+index added to both Step 4 code paths:
- `run_command.py` Step 4 in both `_process_one_sample` (multi-sample)
  and `_run_single_sample` (single-sample): after
  `restore_polya_softclips()` completes, the output BAM is now sorted
  (via `pysam.sort`) and indexed (via `pysam.index`). Without
  sort+index, the BAM was written unsorted and had no `.bai` —
  unusable in IGV/samtools.
- Multi-sample path also gained the missing
  `threads=getattr(args,'threads',4)` argument to
  `restore_polya_softclips()`.
- Bundled validation data: `rectified_pA_tail_soft_clipped.bam`
  generated from `rectified_pA_tail_trimmed.bam` +
  `wt_by4742_rep1_polya_trim_metadata.parquet`. 32/36 reads have
  poly-A/adapter soft-clips restored; 3 replacement reads (f8050895,
  7d5e8dc2, 72557a9a) have no trim metadata (new UUIDs) and are left
  unchanged.
- `TestPolyASoftClippedBam` added to `tests/test_validation_reads.py`
  (4 tests): checks 36 primary reads, coordinate-sorted+indexed,
  reads-with-metadata have ≥1 soft-clip base, replacement reads
  present.
- All 708 tests pass.

**v3.1.8 (2026-04-22):** Validation reads updated to DRS pre-trim
mapPacBio alignments for Cat6/Cat7/Cat9:
- `validation_reads.bam`: all 12 Cat6/Cat7/Cat9 reads now use
  DRS-trimmed mapPacBio alignments. Three reads replaced entirely
  because their pre-DRS mapPacBio alignments lost the N op after DRS
  trimming: `ba761413→f8050895` (cat6_plus_2, chrII:45644–45977),
  `64f4da08→7d5e8dc2` (cat6_minus_1, same intron 60193–60697),
  `5c59f0bc→72557a9a` (cat7_minus_2, chrIII:104435–104495). Nine reads
  updated in place.
- BAM header reduced to 16 nuclear chromosomes only (chrI–chrXVI);
  `chrmt`/`chrMito` removed.
- `corrected_reads.tsv` and rectified BAMs regenerated from DRS-trimmed
  validation set.
- Test updates: Cat7 `EXPECTED_JUNCTIONS` updated (cat7_plus_1:
  138856→138864 start; cat7_plus_2: 595736→595739; cat7_minus_2:
  882352-882702→104435-104495 new read); Cat9 `RAW_JUNCTIONS` updated
  for plus-strand reads (cat9_plus_1: 555824→555825; cat9_plus_2:
  439324→439321). Minus-strand raw junctions unchanged.
- Root cause of mapPacBio BAM unmapped-duplicate issue: mapPacBio
  emits unmapped (flag=4) copies of reads alongside the primary
  alignment; these pass `not is_secondary and not is_supplementary`
  checks and must be explicitly skipped with `read.is_unmapped` guard.
- All 704 tests pass.

**v3.1.6 (2026-04-21):** DRS pre-trimming integrated as default Step 0+4
in `run-all`:
- `--drs` flag added to `rectify run-all` (cli.py). When set with a BAM
  input, Step 0 (`trim_drs_bam_polya`) runs before alignment and Step
  4 (`restore_polya_softclips`) runs after correction. FASTQ inputs
  are unaffected.
- `_run_single_sample()` docstring updated to list all 5 steps (0–4).
- `ARCHITECTURE.md` updated: mermaid diagram corrected (Step 3=analysis,
  Step 4=restore-softclip), `run-all` description expanded with DRS
  wiring details, `junction_refiner.py` added to directory listing
  and Layer 4 narrative (Module 2H, ⑦).
- `README.md` accuracy fixes: `rectify run` → `rectify run-all`
  (command didn't exist); Module 2H scoring description updated
  (W/max_slide unused, HP-anchored semi-global DP replacing
  "split-alignment"); tiebreaker priority corrected (score → canonical
  → annotation → shift, not score → shift → canonical → annotation).
- `junction_refiner.py` `_score_junction` docstring Returns line
  corrected to match tier-1-only implementation (`min_k t1(k)` not
  `t1+t2`).

**v3.1.5 (2026-04-21):** `process_bam_streaming_parallel` — two-level
checkpoint/resume for `rectify correct`:
- New `--checkpoint-dir DIR` flag (performance group in `cli.py`). Has
  no effect without `--streaming`.
- **Scan-phase checkpoint**: after `run_variant_aware_scan` (~30–60
  min), the `VariantAwareHomopolymerRescue` object is pickled to
  `<checkpoint_dir>/rescue_scan.pkl`. On resume the pickle is loaded
  and the scan is skipped entirely.
- **Region-phase checkpointing**: after each region batch completes, a
  `region_NNNN.done` sentinel file is written to the checkpoint dir.
  On resume, done regions are filtered from `_regions_to_run` and
  `_orig_idxs`; the partial output TSV is opened in append mode; the
  header line is not re-written;
  `_rebuild_pos_counts_from_partial` re-reads the partial TSV to
  restore `_pos_counts` for the final position index.
- **Map function**: `pool.imap` (ordered) is used when
  `checkpoint_dir` is set so that `_batch_num → _orig_idxs[_batch_num]`
  is stable; `pool.imap_unordered` is used otherwise (unchanged).
- **Failure behaviour**: when `--checkpoint-dir` is set, partial output
  and sentinels are preserved on failure (warning logged) instead of
  being deleted. Re-running the identical command resumes from the
  last completed region.
- Wired through `correct_command.py` via `getattr(args,
  'checkpoint_dir', None)`.

**v3.1.4 (2026-04-21):** `refine_read_junctions` — candidate guard +
adaptive tie-break scoring:
- **Candidate guard**: non-canonical (tier ≥ 4) AND non-annotated
  (novel) junctions are now excluded as replacement candidates during
  scoring. They were previously considered and could "win" by scoring
  0.0 (aligner placed a non-canonical junction with a perfect
  split-alignment score), then get blocked by a post-selection guard
  that fired `continue` on the N-op, skipping all other candidates
  including the correct annotated GT-AG junction. The guard is now
  applied during candidate selection instead: `if is_alt==1 and
  tier>=4 and is_novel==1: continue`.
- **Adaptive tie-break ordering**: when the current N-op junction is
  non-canonical (tier ≥ 4), the scoring tuple uses `(score, tier,
  is_alt, ...)` so a canonical alternative wins at equal
  edit-distance. When the current junction is acceptably canonical
  (tier < 4), the tuple uses `(score, is_alt, tier, ...)` so the
  current junction is preferred at equal score (prevents displacing
  correctly-placed reads at annotated non-GT-AG junctions like TFC3's
  RAG site). Controlled by `tier_beats_alt = current_tier >= 4`.
- **Impact**: `cat9_plus_2` (read 00a1e01e, chrVII RPL30, current
  non-canonical GA-GG tier=6) now correctly refined to annotated GT-AG
  (439093,439323). TFC3 reads at annotated RAG 3'SS (151006,151096,
  tier=1) remain stable against the nearby YAG alternative
  (150989,151096, tier=0). All 704 tests pass.

**v3.1.3 (2026-04-21):** `--aligner-bams` `aligner:path` prefix stripping
bug fix:
- `correct_command.py` `_strip_aligner_prefix()`: new private helper
  that strips the `aligner:` prefix (e.g.
  `"minimap2:/path/to/file.bam"` → `"/path/to/file.bam"`) before
  storing BAM paths in `config['aligner_bams']`. Without the fix, pysam
  received `"minimap2:/path/..."` as a URL scheme → "[Errno 93]
  Protocol not supported" → 0 novel junctions loaded from aligner BAMs.
- Root cause: `--aligner-bams` accepts the same `aligner:path` format
  as `rectify consensus`, but `collect_junctions_from_bam` needs a
  plain file path. The prefix was never stripped before pysam.
- Impact: With the fix, `rectify correct --aligner-bams minimap2:/path
  gapmm2:/path` loads all aligner-specific novel junctions into the
  pool (377 vs 362 for the validation dataset). Previously
  annotation-only junctions were used even when aligner BAMs were
  specified.
- `_strip_aligner_prefix` checks `'/' not in prefix` to distinguish
  aligner names from Windows-style drive letters (future-proofing);
  plain paths without `:` pass through unchanged.

**v3.1.2 (2026-04-21):** Cat9 validation reads added (Module 2H junction
refinement):
- `validation_reads.bam` expanded from 32 → 36 reads; aligner BAMs
  updated to 36 reads each (uLTRA/deSALT stay at 34/36 due to dev-run
  coverage).
- Four Cat9 reads: `cat9_plus_1` (00a1c9b3, chrVII:555824→555830, +),
  `cat9_plus_2` (00a1e01e, chrVII:439089→439093, +), `cat9_minus_1`
  (0b3b593b, chrXV:900760→900767, −), `cat9_minus_2` (d3357db5,
  chrXV:900760→900767, −). Each has `XG=cat9_junction_refine`.
- `TestCategory9JunctionRefinement` class in
  `tests/test_validation_reads.py`: `corrected_with_aligner_bams`
  fixture runs `rectify correct --aligner-bams --annotation`; verifies
  corrected junction ≠ wrong junction for all 4 reads.
- `TestBamIntegrity` updated: 32→36 reads, cat9 labels added to
  `test_all_labels_present`, `cat9→cat9_junction_refine` in
  `test_category_tags`, cat9 added to `test_strand_balance`.
- Bundled validation outputs regenerated: `corrected_reads.tsv`,
  rectified BAMs, PROVENANCE.json all updated to reflect 36 reads.

**v3.1.1 (2026-04-21):** `_score_junction` + `refine_read_junctions` —
two correctness fixes:
- **`range(L)` fix**: the k-loop in `_score_junction` previously ran
  `range(L+1)`. At k=L, `q1 = rescue[L:]` is empty →
  `_score_hp_anchored` returns 0.0 for every candidate regardless of
  junction quality. All candidates tied at 0.0, making selection
  arbitrary. Changed to `range(L)`: q1 always has ≥1 base; scoring is
  always discriminating.
- **`is_alt` tiebreaker in `refine_read_junctions`**: added `is_alt`
  (0 if candidate = current N-op, 1 otherwise) as the second tuple
  element after `score`. When multiple candidates tie (e.g. two
  junctions share the same `intron_end` and both score 0.0), the
  existing junction is preferred over alternatives, preventing
  spurious displacement of already-correct reads.
- **Tests**: 9 previously failing tests now pass; all 41 tests in
  `test_junction_refiner.py` pass. Total suite: 698 tests pass, 4
  skipped. Key reads: RPL20B `0b3b593b` corrected to `[900767,901193)`;
  TFC3 annotated junction stable against alternative
  `[150989,151096)`.

**v3.1.0 (2026-04-20):** Module 2H — post-consensus N-op junction
refinement (`junction_refiner.py`):
- New module `rectify/core/splice/junction_refiner.py`: for every N-op
  in every consensus read, tests all candidate junctions within a
  search radius and replaces imprecise N-op boundaries with the best
  sequence-supported junction.
- Scoring is **sequence-first**: hp_score (split-alignment with
  homopolymer-aware edit distance) is the primary criterion;
  canonical GT-AG and annotated status are tie-breakers only.
  Annotation NEVER overrides a better-scoring junction.
- Split-alignment: query window of ±W bp around the current N-op
  split point; query split can slide ±max_slide bp; scores both exon2
  and exon1 genomic context simultaneously.
- Fast path: reads already at an annotated canonical-tier-0 junction
  skip scoring entirely (255× speedup over naive approach).
- `max_boundary_shift=50`: prevents false-positive matches from
  junctions in neighbouring genes. The `search_radius=5000` discovers
  candidates; `max_boundary_shift` constrains individual endpoint
  shifts.
- CIGAR surgery encodes boundary changes as I/D ops (not M
  adjustments) to preserve both ref and query spans. MD/cs tags
  stripped from modified reads to avoid pysam arithmetic errors on
  CentOS 7.
- Wired into `rectify correct` via `--aligner-bams`; also exposes
  `--junction-hp-pen`, `--junction-search-radius`, `--junction-window`,
  `--junction-max-slide`, `--junction-max-boundary-shift`.
- Key fixed reads: b2c4d195 (27S → 426N+exon1), f6590560 (422N1M1D1M
  → 426N2M), 0b3b593b (431N wrong boundaries → 430N correct exon1
  end).
- 41 tests in `tests/test_junction_refiner.py` covering RPL20B, GCR1,
  TFC3, RPL22B, SRC1. (9 were failing at initial merge; fixed in
  v3.1.1.)

**v3.0.4 (2026-04-20):** `rescue_3ss_truncation` — minus-strand
soft-clip truncation bug fix:
- `_extract_5prime_rescue_seq` returns `query_seq[n - last_imp_q:]`
  for minus-strand reads (the rightmost portion). The truncation to
  `five_clip` bases was using `[:five_clip]` (leftmost), taking
  aligned exon2 body instead of the actual soft-clip bases.
- Fix: `rescue_seq[-five_clip:]` for minus strand,
  `rescue_seq[:five_clip]` for plus strand (unchanged).
- Impact: b2c4d195 (27S), 7bf94550 (1S), f6590560 (2S), 08f6ddf7 (5S)
  all now correctly rescued to `five_prime_corrected=901193`.
  Previously all fell through to `rescue_type='proximity'` with
  `rescued=False`.

**v3.0.3 (2026-04-16):** `find_polya_boundary` — poly-A tail
trailing-base false-stop guard:
- `find_polya_boundary` in `indel_corrector.py` now detects false stops
  where the trailing base of a poly-A tail (e.g., a T at the very end
  of `...AAAAAAAAAAAAAAAAT`) coincidentally matches a genomic base
  (T=T) at the alignment boundary, causing the backward scan to stop
  prematurely at the poly-A/exon junction rather than continuing to
  the true exon body.
- Fix (plus strand): Before accepting a candidate stop (`rb==gb,
  gb!='A'`), inspect the K=4 positions to the left. If all K have
  `rb='A'` AND at least one has `gb≠'A'` (unmistakably poly-A tail
  context), the candidate is skipped and scanning continues leftward.
- Fix (minus strand): Symmetric check — inspect K=4 positions to the
  right; if all have `rb='T'` AND at least one has `gb≠'T'`, continue
  scanning rightward.
- Root cause in cat6_plus_2 (RPL19B, ba761413): alignment ends
  `...AAAAAAAAAAAAAAAAT` at 169492 (T=T exact match). The 14 A
  positions before it (169478–169491) all have `rb='A'` with several
  mismatches (`gb≠'A'`) — pure poly-A tail. With the guard, the scan
  continues to 169476 (G=G, true exon end). `corrected_3prime` moves
  16 bp upstream from 169492 to 169476; NET-seq refinement stays at
  169476 (signal=75.0). Bedgraph now shows signal at 169476 rather
  than 169491.
- All 659 tests pass.

**v3.0.2 (2026-04-16):** `clip_read_to_corrected_3prime` /
`softclip_read_to_corrected_3prime` — terminal D/N stripping:
- Both functions now strip trailing D/N ops (plus strand) or leading
  D/N ops (minus strand) that are left dangling after the CIGAR walk
  loop. The bug arose when `corrected_3prime` fell within a deletion
  span: the loop removed all query-consuming ops to the right of the
  D (satisfying `n_ref_removed >= n_ref_clip`) and exited before
  reaching the D, then appended `H`/`S` directly after it, producing
  invalid CIGAR strings like `4D6H`.
- Example: read 299e1402 (chrII plus strand) CIGAR previously ended
  `5X3=4D6H` — a terminal deletion before a hard-clip that violated
  SAM spec. Now the 4D is stripped before the H is appended.
- Root cause: `corrected_3prime` from `find_polya_boundary` can land
  inside a deletion span (a reference position with no corresponding
  query base), which is a valid CPA but requires CIGAR surgery to
  snap to the last real query base.
- Fix applied in four locations: `clip_read_to_corrected_3prime` (plus
  + minus strand paths) and `softclip_read_to_corrected_3prime` (plus
  + minus strand paths).
- New test file: `tests/test_bam_writer.py` (11 tests covering normal
  clipping, terminal-D stripping, in-deletion-span corrected positions,
  and N-op stripping for spliced reads).

**v3.0.1 (2026-04-15):** `clip_intronic_tail_5prime` — off-by-one fix +
trailing-I/S stripping + existing-H handling +
`_MIN_SC_FOR_JUNCTION_EXTENSION` guard:
- **Off-by-one fix (minus strand)**: Exit condition changed from `<=
  clip_boundary + 1` → `<= clip_boundary`. Previously reads ending at
  `reference_end = intron_start + 1` (last mapped base = `intron_start`
  = first intron base) were left unclipped, appearing as red T→C
  mismatches in IGV. Now they are trimmed so the last mapped base is
  `intron_start - 1` (last exon base).
- **Trailing I/S stripping**: Before the main ref-consuming trim loop,
  `clip_intronic_tail_5prime` now explicitly strips trailing I
  (insertion) and S (soft-clip) ops from the 5' end. These bases lie
  at/past the intron boundary and were silently left in the CIGAR
  under the old code when `reference_end` was already at the boundary.
- **Existing H handling**: Any pre-existing trailing H op is extracted
  before processing and merged back at the end. This prevents
  double-counting query bases that were already removed from
  `query_sequence` in a prior call.
- **`_MIN_SC_FOR_JUNCTION_EXTENSION = 3`**: New module-level constant
  in `bam_writer.py`. All three write functions (`write_corrected_bam`,
  `write_softclipped_bam`, `write_dual_bam`) now only call
  `extend_read_5prime_for_junction_rescue` when the 5' soft clip is ≥
  3 bp. Reads with 1-2 bp soft clips fall through to
  `clip_intronic_tail_5prime` instead, preventing spurious Cat3
  exon-N-exon CIGAR surgery on single-base alignment artefacts. An
  `_extended` flag prevents `clip_intronic_tail_5prime` from undoing
  a successful Cat3 extension.
- **`bam_processor.py` `_sc_at_5p == 0` restriction removed**:
  `five_prime_intron_clip_pos` is now set for ALL `_in_intron` reads
  regardless of whether a soft clip is present at the 5' end. Reads
  with small soft clips now get the clip applied rather than being
  directed to Cat3 extension.
- Impact: zero reads with last base at intron_start (X mismatch),
  zero trailing I ops at boundary, 292 reads properly clipped to
  `ref_end = intron_start` (last base = last exon base); all 2,073
  reads pass CIGAR/sequence-length validation.

**v3.0.0 (2026-04-15):** `clip_intronic_tail_5prime` — BAM sequence
trimming fix + generalised intron-clip trigger:
- `clip_intronic_tail_5prime` in `bam_writer.py` now trims
  `query_sequence` and `query_qualities` to match the new CIGAR after
  adding H ops. Previously only the CIGAR was updated, leaving
  `query_sequence` too long → pysam wrote malformed BAM records →
  `samtools sort: truncated file` error and silent loss of all
  intronic-snap clip corrections.
  - Minus strand (clip from right): `read.query_sequence =
    seq[:-clipped_query_bases]`.
  - Plus strand (clip from left): `read.query_sequence =
    seq[clipped_query_bases:]`.
- `bam_processor.py` `five_prime_intron_clip_pos` assignment
  generalised: previously only set for `rescue_type == 'intronic_snap'`
  (Case 4). Now set for ALL rescues (Cases 1/2/4) where the
  alignment's 5' end falls inside the rescued intron AND there is no
  5' soft-clip. Case 1/2 reads rescued with no soft clip (e.g.
  d36f9748: `16M` exon match, alignment ending 16 bp inside intron)
  also need the intronic tail hard-clipped.
- Impact: clipped BAM count rose from 1,287 → 1,355; all 2,073 reads
  pass CIGAR/sequence-length validation.

**v2.9.1 (2026-04-12):** Cat2 soft-clip rescue CIGAR surgery:
- `rescue_softclip_at_homopolymer` now stops at `A` (plus strand) or
  `T` (minus strand) in the soft-clip, preventing poly-A tail bases
  from matching genomic A-runs. Fixes shifted corrected positions for
  cat2_plus_2 (+10→+9), cat2_minus_1 (-17→-10), cat2_minus_2
  (-12→-11).
- New `extend_read_3prime_for_softclip_rescue()` in `bam_writer.py`:
  converts the 3' soft-clip to `{D}D{M}M{poly-A}H|S`, making the true
  RNA 3' end visible in IGV.
- Cat2 rescue metadata (`sc_homopolymer_extension`, `sc_rescued_seq`,
  `sc_original_softclip_len`) now written to corrected_reads.tsv and
  read back by bam_writer.
- Bundled BAMs renamed: `rectified.bam` → `rectified_corrected_3end.bam`;
  `rectified_softclip.bam` → `rectified_pA_tail_trimmed.bam`.

**v2.9.7 (2026-04-14):** `_iter_name_grouped_bams` — natural sort K-way
merge fix:
- The K-way merge in `consensus.py` used Python `min()` (lexicographic)
  to pick the next read ID, but the name-sorted BAMs use
  `SS:queryname:natural` (samtools natural sort, digit runs compared
  as integers). For UUID-format read names, these orderings diverge:
  `98297e97` (key=98297) sorts BEFORE `0633141e` (key=633141) in
  natural sort, but AFTER it lexicographically ('9' > '0').
- When a read present only in aligner A sits between two reads in
  natural sort, the lexicographic merge incorrectly pulled the shared
  read from aligner B first (as a single-aligner group), then later
  processed aligner A's copy as a separate single-aligner group. Both
  single-aligner groups "won" by default → duplicate records with
  different aligners/strands in the output BAM.
- Fix: added `_natural_sort_key()` (splits on `(\d+)`, converts digit
  tokens to `int`) and passed it as `key=` to the `min()` call in
  `_iter_name_grouped_bams`. The merge now uses the same natural sort
  as the BAMs.
- Impact on RPL19B+RPL20B test: mixed-strand duplicate QNAMEs in
  consensus BAM dropped from 51 → 1 (the 1 remaining is an inherent
  FASTQ-duplicate edge case where copies map to opposite strands).
  Plus-strand soft-clipped reads on chrXV dropped from 53 → 9
  (residual are other-region reads + FASTQ-dup pairs).

**v2.9.9 (2026-04-15):** `rescue_3ss_truncation` — homopolymer-aware
edit distance + 3'SS acceptor tiebreaker:
- New `_hp_edit_distance(s1, s2) -> float`: indels within a
  homopolymer run (indel base equals its immediate neighbour in the
  same string) cost 0.5 instead of 1.0; substitutions always cost
  1.0. Used for all four edit-distance calls in
  `rescue_3ss_truncation` (inner-loop candidate scoring, outer
  `ed_exon`, and intron comparison `ed_intron`). Replaces
  `_edit_distance` (integer Levenshtein) for all splice-rescue
  scoring; `_edit_distance` is retained for other callers.
- New `_ACCEPTOR_PRIORITY_PLUS` and `_ACCEPTOR_PRIORITY_MINUS`
  module-level dicts encode 3'SS acceptor dinucleotide priority:
  AG=0, CG=1, TG=2, AT=3, other=4 (plus strand, last 2 bases of
  intron); minus strand uses RC equivalents (CT=0, CG=1, CA=2,
  AT=3).
- `_acceptor_priority` is computed once per junction (outside the
  shift loop) and added as the 5th element of the outer tiebreaking
  tuple: `(ed_exon, not_canonical_donor, not_in_amb, shift_abs,
  acceptor_priority)`.
- Outer tracking now also records `best_in_amb`, `best_shift_abs`,
  `best_acceptor_priority` alongside the existing `best_ed` and
  `best_is_canonical`.

**v2.9.8 (2026-04-15):** `rescue_3ss_truncation` — two-phase discovery
+ canonical refinement, ±5 bp baseline for both strands:
- Plus and minus strand branches now both use `max(5, amb)` as the
  shift range floor (previously `max(1, amb)`). The ±5 bp baseline
  catches imprecise aligner annotations (e.g. off by 3-4 bp) even
  when local sequence ambiguity is only 1 bp.
- Minus strand canonical check expanded: `in ('AC', 'GC')` (RC of
  GT/GC) instead of `== 'AC'` only. Mirrors plus strand which uses
  `in ('GT', 'GC')`.
- Minus strand now uses the same tuple scoring as plus strand:
  `(not_canonical, not_in_amb, shift_abs)` — canonical wins first,
  then within-ambiguity-window shift, then smallest |shift|.
  Previously used a simpler chain lacking `_in_amb` refinement.
- Three production-data failed rescues motivating this fix (not yet
  in validation set): chrI:142,611-142,653 (EFB1/YAL003W, plus
  strand), chrII:366,454-366,654 (YBR062C, minus strand),
  chrXI:20,510-20,584 (plus strand, poly-A run).

**v2.9.6 (2026-04-14):** `rescue_3ss_truncation` — splice-site
ambiguity resolution via data-driven shift range and canonical GT-AG
preference:
- For each candidate junction, the exon-end scan tries a range of
  donor/acceptor position shifts derived from the actual run-length
  of matching bases at the boundary: `r_amb` = consecutive intron
  bases (going rightward from the annotated donor) that equal the
  last exon base; `l_amb` = consecutive exon bases (going leftward)
  that equal the first intron base. The shift range is `[-max(5,
  l_amb), max(5, r_amb)]`, capped at ±15. The minimum ±5 ensures
  off-by-one errors (e.g. 1-based GFF vs 0-based) and imprecise
  aligner annotations are always caught.
- Scoring priority: lower edit distance wins first; among ties,
  canonical GT/GC donor (plus strand) or AC/GC (minus strand
  genomic) wins over non-canonical; among remaining ties,
  within-ambiguity-window shifts win over outside; then smallest
  |shift| wins.
- The data-driven range is critical for correctness: a fixed wide
  range (e.g. ±15) would spuriously try shifts deep into the intron,
  finding false perfect matches when exon1-end + intron-start looks
  identical to intron-end + exon2-start (e.g. `TTTCAG-GTA` at both
  the 5'SS and 3'SS). The run-length bound prevents this.
- For long homopolymer runs at the junction, the range automatically
  extends to cover the full ambiguous window (up to 15 bp). For
  typical 1-bp repeats (like G|G at YAL003W), the range is ±5
  (dominated by the baseline).
- Example: YAL003W (chrI) read 79f61403. ref[142252]=G (last exon
  base), ref[142253]=G (first intron base, canonical GT donor at
  142253-142254). With the wrong junction (intron_start=142254),
  shift=-1 gives intron_start=142253 (GT canonical) with the same
  edit distance → canonical wins → five_prime_corrected=142252,
  exon_cigar=14M1D7M1D1M. Deterministic across all set iteration
  orders.

**v2.9.5 (2026-04-14):** `rescue_3ss_truncation` — Case 4 intronic-snap
rescue:
- New `_get_n_op_intervals()` helper extracts (start, end) genomic
  intervals for every N-op in a CIGAR.
- New Case 4 in `rescue_3ss_truncation()`: fires when `align_5prime`
  is strictly inside `[intron_start, intron_end)` AND no existing
  N-op already covers that intron (checked with
  ±`junction_proximity_bp` tolerance). Cases 1–2 already handle this
  via terminal-error sequence matching; Case 4 handles the remainder
  (clean alignments with D-ops or plain truncation inside the intron,
  no detectable terminal errors).
- Action: snaps `five_prime_corrected` to `intron_end` (minus strand)
  or `intron_start` (plus strand) — the exon-1-side splice donor
  boundary. Records `rescue_type = 'intronic_snap'`, `rescued = True`.
- Works for all aligner representations: does not depend on CIGAR op
  type, only on final `align_5prime` position.
- Impact on RPL20B test set: reduced reads with `five_prime_position`
  inside annotated intron from 170 → 0; added 281 `intronic_snap`
  rescues at `intron_end = 901193`.

**v2.9.4 (2026-04-14):** `find_polya_boundary` — N-op boundary guard
for spliced minus-strand reads:
- `find_polya_boundary` in `indel_corrector.py` now records the first
  N-op (intron skip) boundary encountered while parsing the CIGAR
  for minus-strand reads and limits the forward scan to positions
  before that boundary (`scan_limit`). Previously the scan silently
  skipped N ops (they are absent from `aligned_positions`) and could
  find a spurious non-T match in a downstream exon, producing a huge
  erroneous correction that spanned the intron.
- If no non-T match is found before the N boundary and the last
  pre-N position is T (poly-A zone), `first_n_start` is returned as
  the CPA — the intron boundary is the natural exon-end for the
  terminal exon.
- Plus-strand backward scan is unchanged: for plus-strand reads the
  3' terminal exon is to the right and crossing N ops leftward into
  earlier exons is legitimate (the entire last exon may be in the
  poly-A zone, requiring the scan to continue into exon N-1 to find
  the true non-A boundary).
- Example: YIL145C (chrIX) minus-strand read a9706bbe with CIGAR
  `1S11M223N1D34M...`. All 11M positions (76016-76026) are T's; scan
  previously crossed the 223N and found seq[12]='G'==genome[76251]='G',
  giving a 235bp correction. With the fix, scan stops at scan_limit
  (exon 1 positions only), finds no non-T match, and returns
  first_n_start=76027 as the CPA (11bp correction). The 12-row
  ambiguity problem is resolved.

**v2.9.3 (2026-04-14):** `find_polya_boundary` — large-deletion
pre-scan for poly-A over-calling artifacts:
- `find_polya_boundary` in `indel_corrector.py` now detects large
  deletions (≥5bp) within 50bp of the RNA 3' end that are poly-A
  over-calling artifacts. The artifact: minimap2 over-extends the
  poly-A tail alignment into a mismatching region, then uses a large
  deletion to bridge back to the next exon match. The false-positive
  non-T (minus strand) or non-A (plus strand) matches before the
  deletion caused the walkback to stop too early (under-correcting
  by tens of bp).
- Pre-scan is gated on whether the alignment starts in a poly-T/
  poly-A context (first real aligned position has `gb='T'` for minus
  strand, `gb='A'` for plus strand), preventing false activation for
  reads whose 3' end is NOT in a poly-A zone (e.g., exon body
  deletions).
- Example: RPL20B (YOR312C) minus strand reads with CIGAR
  `3=1X1=2X3=18D...` at chrXV:900,120 — previous correction was 2bp
  (stopped at coincidental G match before the 18D), now correctly
  28bp (stops at first exon match after the 18D, within 2bp of the
  annotated CPA at 900,150).

**v2.9.2 (2026-04-14):** Cat3 5' rescue — mapPacBio intrusion fix:
- `rescue_3ss_truncation()` in `splice_aware_5prime.py` now handles
  the case where mapPacBio's alignment extends INTO the upstream
  intron (`align_5prime` lands between `intron_start` and
  `intron_end`). Previously the `dist < 0` guard skipped these
  junctions entirely before any sequence comparison; now they are
  treated as `dist = 0` (touching the boundary) and proceed to
  exon-sequence alignment. Applies symmetrically for both strands.

**v0.9.0-dev (2026-04-12):** Parallel alignment infrastructure:
- `rectify split` — splits FASTQ.GZ into N equal chunks (round-robin
  interleave) for SLURM array jobs; `--generate-slurm` emits
  `run_array_align.sh` (array 0-79) and
  `run_merge_and_consensus.sh`.
- `rectify consensus` — aligner selection on pre-built BAMs; accepts
  `aligner:path` pairs; used after merging per-chunk per-aligner
  BAMs from array jobs.
- `rectify install-aligners` — downloads/installs external aligners;
  `--check` shows status; `--all` installs everything; deSALT gets
  vendored binary (Linux/x86_64) or source build fallback.
- deSALT vendored binary bundled at
  `rectify/data/bin/linux_x86_64/deSALT` (v1.5.6, 773 KB);
  `_get_vendored_desalt()` in `multi_aligner.py` resolves
  automatically when deSALT not on PATH.
