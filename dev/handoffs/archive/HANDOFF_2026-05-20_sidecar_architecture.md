# Session Handoff Update — read_num sidecar implemented

**Date:** 2026-05-20
**Branch:** `drs-validation-rebuild`
**Status:** Option B (RN aux tag hybrid) is implemented in the working tree.

Implemented pieces:
- `rectify split` emits RN-tagged chunk FASTQs and `<sample>.read_num_sidecar.parquet` plus provenance JSON.
- Comment-stripping aligner wrappers inject `RN:i` into BAM records after QNAME validation.
- Consensus uses RN keys only when all BAM inputs have RN, with normalized-QNAME fallback for old BAMs.
- `rectify consensus --read-num-sidecar` restores original FASTQ-comment aux tags without overwriting existing tags.

Remaining validation gap: full production end-to-end smoke on the validation fixture with real aligners was not run in this local pass; focused unit/regression tests are green.

---

# Session Handoff — `read_num + sidecar` architecture design

**Date:** 2026-05-20
**Branch:** `drs-validation-rebuild`
**HEAD:** `45cbc13` — `fix(qname): cross-aligner QNAME hardening + auto-sanitize validator`
**Scope:** This handoff is *specifically* for the deferred sidecar
architecture design discussion held this session. Today's actual code
work (the defensive QNAME layer) is landed and tested; what's deferred
is the next-generation cross-aligner identity scheme.

---

## 1. What was done

- Cross-aligner QNAME hardening landed end-to-end (`45cbc13`) — see commit
  message for the full file list. 113 tests pass in the QNAME + consensus
  + restore suite. Three of six issues from
  `scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md` fixed in the
  same commit (#1 chimeric normalize, #2 TSV-side normalizer parity,
  #5 strand guard).
- Cross-aligner audit run: `scripts/diagnostics/qname_mutation_survey/`
  contains the synthetic-FASTQ survey (REPORT.md), the cross-aligner
  punch list (EDGE_CASES.md), and the 449-row characterization TSV.
- The sidecar architecture design was discussed but **not implemented**.
  This handoff captures the design conclusion so the next session can
  resume without re-deriving it.

---

## 2. What's verified

- `git log --oneline -1` → `45cbc13 fix(qname): cross-aligner QNAME hardening + auto-sanitize validator`
- `pytest tests/test_qname_sanitizer_and_validator.py tests/test_normalize_read_name.py tests/test_consensus_selection.py tests/test_corrected_consensus_tiebreaker.py tests/test_chimeric_consensus.py tests/test_parallel_aligner_schedule.py tests/test_gapmm2_seq_restore.py -q`
  → 113 passed, 4 skipped, ~14 s
- Production DRS spot-check confirmed bug-free at HEAD: set1 wt_rep1
  (Sherlock) had 0 occurrences of `_pt:i:` in first 1000 mapPacBio
  QNAMEs. Today's work is preventive for the in-place flows.

**NOT VERIFIED:**
- No sidecar prototype written — design only.
- The "cross-aligner tag drift" hypothesis (EDGE_CASES.md issue #3) is
  backed by code inspection, not a production read-count diff. If you
  want a quantitative impact number before committing to the sidecar
  work, run `samtools view <consensus.bam> | awk '$0 !~ /XA:Z/'` on a
  recent cDNA production output and count records — every record
  without an XA:Z: tag is a winner from a non-minimap2 aligner that lost
  its tag set.

---

## 3. Open items

### The sidecar architecture (THE deferred item this handoff is for)

**Motivating problem:** EDGE_CASES.md issue #3 — cDNA FASTQ tags
(XA/XC/XF/XU/XR) are propagated to BAM aux records *only* by minimap2
(via its `-y` flag). Every other aligner wrapper strips FASTQ comments
before invoking the aligner (`_sanitize_mpb_fastq` for mapPacBio,
`_clean_fastq` for deSALT/gapmm2; uLTRA and bwa have no equivalent
mechanism). When mapPacBio/deSALT/gapmm2/uLTRA wins consensus for a
read, the output BAM record has **no cDNA metadata tags**. Downstream
cDNA analyses (tail length, cluster size, full-length tier, UMI) read
missing values for any read where minimap2 was not the winner. This is
structural; it cannot be patched piecewise without changing how the
pipeline carries per-read metadata across the aligner stage.

**The proposal (user's framing, refined this session):**

1. At chunking time, assign each FASTQ record a sequential integer
   `read_num`. Globally unique per sample — either via per-chunk prefix
   (`c01_000001`) or via a globally synchronized counter.
2. Write a parquet sidecar mapping `read_num → original_qname,
   fastq_comment, chunk_id`. Original FASTQ stays untouched as the
   source of truth.
3. Use the `read_num` as the cross-aligner primary key.
4. `(read_num + seq + qual)` is the verification fingerprint —
   collision-free by construction because `read_num` is bijective.
   `(seq + qual)` alone collides on PCR duplicates; with `read_num`,
   PCR duplicates are still distinguishable.
5. At consensus emission, restore the original FASTQ tags from the
   sidecar onto the winning BAM record, regardless of which aligner won.

**Why this is the right architectural target:**
- Solves issue #3 cleanly. Tag restoration is centralized at consensus
  emission, not scattered per-wrapper.
- Retroactively obsoletes most of today's defensive machinery — the
  QNAME-mutation class can't exist if QNAMEs are integers (or carried
  in an aux tag that aligners can't mutate). The validator and
  sanitizer become belt-and-suspenders rather than primary.
- Centralizes the cross-aligner identity problem in one place
  (chunker + consensus restore) instead of N per-aligner sanitizers.

### Two implementation variants discussed (decide first thing next session)

**A — Full anonymization.** Chunker writes derived FASTQs with integer
QNAMEs only; aligners see integers; consensus restores original QNAMEs
from sidecar. Maximally clean cross-aligner story.

*Costs:*
- Per-aligner BAMs have anonymous QNAMEs. IGV/grep ergonomics suffer.
- Migration burden for every existing UUID-QNAME BAM on disk
  (set1/set2/set3 production runs).
- Chunker refactor touches more code paths (`rectify split`,
  `rectify run-all`, the cDNA pipeline's own chunker per
  `project_cdna_pipeline.md` memory).
- Test fixture migration across hundreds of tests that hard-code UUID
  QNAMEs.
- Restoration step at consensus emission must be reliable — if the
  sidecar is unreachable, output is anonymous forever.

**B — RN aux tag hybrid (RECOMMENDED).** Chunker keeps original QNAMEs,
attaches `RN:i:<num>` as a FASTQ-comment-side aux tag (alongside the
existing cDNA pipeline X-tags). Aligners that natively support tag
passthrough (minimap2 with `-y`) carry RN automatically. Others need a
small post-alignment tag-injection pass keyed on QNAME. K-way merge
prefers RN when available, falls back to normalized QNAME otherwise.

*Costs:*
- One small post-alignment tag-injection pass for mapPacBio / deSALT /
  gapmm2 / uLTRA wrappers (~50 lines per wrapper).
- The sidecar still has to be written at chunking time.

*Benefits:*
- Backwards compatible with existing UUID-QNAME BAMs (RN-less BAMs use
  the QNAME fallback path).
- IGV/samtools/grep all keep working unchanged.
- Today's validator + auto-sanitize stays as belt-and-suspenders
  defense in depth.
- Closes issue #3 (sidecar carries cDNA tags; consensus restoration
  applies them regardless of winner).

**Recommendation: B.** ~80% of the value with ~20% of the disruption.
A is a v3.0 rewrite worth considering only if cross-aligner sanitizer
maintenance becomes intolerable.

### Other open items (not blocking the sidecar work)

- **EDGE_CASES.md #4** (cosmetic tiebreaker comment claims MAPQ but
  code uses only N-ops). Trivial fix; deferred because it's a comment
  drift, not a bug.
- **EDGE_CASES.md #6** (validator samples first N reads sequentially,
  not randomly). Robust enough for current aligner panel; revisit if
  new aligners are added with position-dependent mutation patterns.
- **Pre-existing WIP in working tree** (not touched today): Commit B
  manifest split (`correct_command.py`, `cli.py`,
  `export_merged_tsv_command.py`) currently breaks
  `test_quantseq_rev_integration.py` (6 errors) and
  `test_validation_reads.py` (83 errors). All errors are the same
  pattern: integration tests call `rectify correct -o foo.tsv` then
  open `foo.tsv`, but Commit B emits `foo.manifest.tsv` + per-region
  TSVs. Tests need to call `rectify export-merged-tsv` after `correct`,
  or a back-compat `--legacy-monolithic-tsv` flag needs to land before
  these tests pass again. Out of scope for the QNAME work and the
  sidecar design.

---

## 4. Resume command

**Resume:** read this handoff plus EDGE_CASES.md issue #3 (in
`scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md`) and the
AGENT_FIXES.md entry titled "QNAME pipeline hardening" — together they
are the full context. Then:

1. **Decide A vs B** (recommendation in §3 is B). If picking B, proceed.
   If picking A, write a separate design doc covering the migration
   path for existing BAMs and the test-fixture rewrite before any code.
2. **Locate the chunker entry points** with
   `grep -rn "def.*chunk\|extract_fastq_chunk\|split_long_reads" rectify/ scripts/`
   then walk `rectify/core/commands/split_command.py` and
   `rectify/core/commands/run_all_command.py`. The cDNA pipeline has
   its own chunker — `project_cdna_pipeline.md` memory points at it.
3. **Sketch the sidecar parquet schema** as a code comment in a new
   file `rectify/core/chunking/sidecar.py` BEFORE writing code:
   - `read_num` (int64, primary key)
   - `original_qname` (str)
   - `fastq_comment` (str — full tab-separated comment string, or empty)
   - `chunk_id` (str)
   - `seq_md5` (16 bytes — for the fingerprint verification path)
   - `qual_md5` (16 bytes)
4. **Identify the RN injection points**:
   - minimap2 (`-y`): RN flows automatically if present in FASTQ comment.
     Today's chunker can be extended to write
     `@<qname>\tRN:i:<num>\t<existing-tags>` and minimap2 just propagates.
   - mapPacBio: `_sanitize_mpb_fastq` currently strips ALL comments.
     Need to capture RN before strip and re-inject into the BAM with
     `pysam.AlignmentFile` post-run. Pattern identical to today's
     auto-sanitize streaming rewrite — look at
     `validate_post_alignment_qnames` for the template.
   - deSALT/gapmm2: `_clean_fastq` strips comments. Same pattern as
     mapPacBio.
   - uLTRA/bwa: no FASTQ-comment passthrough mechanism. Same as above.
5. **Sketch the consensus restoration pass** in `consensus.py` near the
   line-489 normalize site: for each winning read, read the `RN:i` tag
   (or fall back to QNAME-based lookup if no RN), fetch the original
   tags from the sidecar, set them on the output record. Backwards-
   compatible: if no sidecar, behavior is identical to today.
6. **Write an integration test** before any production code: a tiny
   chunker → 5-aligner-mock → consensus pipeline that exercises the
   tag restoration. Should pass for both the new-style (with sidecar)
   and old-style (no sidecar, falls through to QNAME normalize) paths.

If the next session is short on time, an acceptable scope-reduction is
to implement the chunker + sidecar emission first (read-only — no
restoration pass yet) and run one production sample through to confirm
the sidecar parquet is correct before touching the consensus side.

---

## 5. Files touched

This session committed the entire QNAME defensive layer at `45cbc13`.
The sidecar work touched zero source files. The following references
will be useful next session:

- `scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md` — issue #3
  is the spec for what the sidecar must solve
- `scripts/diagnostics/qname_mutation_survey/REPORT.md` — full survey
  context (per-aligner tag-stripping behavior is documented per wrapper)
- `rectify/core/align/qname_validator.py` — the existing defensive
  layer that stays as belt-and-suspenders. The streaming-rewrite
  pattern at `_rewrite_bam_normalize_qnames` is the template for the
  RN injection pass.
- `rectify/core/consensus/consensus.py:171` — canonical
  `_UNDERSCORE_COMMENT_RE`. The QNAME fallback path in option B will
  reuse this exact regex.
- `rectify/core/align/multi_aligner.py` —
  `_sanitize_mpb_fastq` (lines ~494–525) and `_clean_fastq` (lines
  ~51 onward) are the current comment-strip points. The sidecar work
  needs to capture comments *before* these strip them.
- `AGENT_FIXES.md` — the "QNAME pipeline hardening" entry sets the
  context for why today's work happened and what was explicitly
  deferred.

**No `[uncommitted]` files from the sidecar work** — this handoff
exists precisely because there is no code to commit. The pre-existing
WIP in the working tree (Commit B manifest split, validation BAM
regen, figures, docs, handoffs) is out of scope for the sidecar work
and was not touched.
