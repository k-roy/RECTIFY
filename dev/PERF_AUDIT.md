# PERF_AUDIT.md — hunting per-read over-computation in RECTIFY

**Created:** 2026-05-24 (Kevin's request). **Status:** living doc.

**Why this exists:** the 3'SS-rescue bottleneck below (run-all manifest inline-correct
stall → OOM) is almost certainly *not* the only place RECTIFY does far more work than
necessary. This doc (a) records that case as a worked example, and (b) gives the next
agent a concrete methodology + grep targets + suspect list to find and fix similar
bottlenecks. **Read this before any perf work; add findings as you go.**

RECTIFY is v0.9.x; manifest mode is experimental. These are perf bugs, not correctness
bugs — but at production scale (full junction pools, high-coverage loci, eventually
human introns >100 kb) they make pipelines hang. Fix them, but verify they don't change
output (see Verification below).

---

## Case study: 3'SS rescue = O(reads × full-pool × ambiguity-window) HP-edit-distance

**Where:** `rectify/core/splice/splice_aware_5prime.py` (`_rescue_3ss_truncation_body`,
`_hp_edit_distance`) fed by `rectify/core/bam/bam_processor.py:correct_read_3prime`.

**Symptom:** `rectify run-all --manifest` (DRS, 5% wt+upf1 subset) hung at
`Processing 126 regions across 16 workers / 0%` for 6.5 h then OOM-killed at 64 GB
(Sherlock job 25846844). Align + trim + Module-2H refine all succeeded; the stall is the
inline `correct` parallel region step.

**Cost decomposition (py-spy, 29,999 samples, single-thread full-pool subset = job DIAGE):**
`_hp_edit_distance` (an O(n×m) DP, splice_aware_5prime.py:693-723) = **~86% of CPU**.
- Baseline rescue `_rescue_3ss_truncation_body` = **87%** of that. Two loops scan the
  candidate set and call the DP per candidate: the sequence-rescue loop (line ~1409) and
  the Case-4 intronic-snap loop (~1877/1878).
- The sequence-rescue loop is worse than "per candidate": for EACH candidate it runs
  `for _shift in [±~15]` × `for _off in [0..10]` ≈ **~341 DP calls per candidate**.
- Terminal peel `_terminal_peel_rescue` = 12.9% (its own per-depth scan, fixed first).

**Three compounding root causes — the generalizable lessons:**
1. **Per-read scan of a global collection.** `candidate_junctions` (the junction pool)
   can hold ~17k entries in run-all mode; the body looped over all of them per read.
2. **Vestigial oversized bound.** `_POOL_SEARCH_RADIUS` was **10000** bp, with a comment
   claiming it must cover `junction_proximity_bp` "default 5000" — but the real default
   is **10**. The 10k radius was sized for a parameter that no longer exists, so it pulled
   every intron in a ±10 kb window per read.
3. **Bound keyed on the wrong axis.** The pool index is sorted by `intron_start`; a radius
   on `intron_start` bounds the *intron length*. That happens to work for yeast (introns
   ≤~1.2 kb) but **breaks for human** (introns >100 kb): the donor sits 100 kb from the
   acceptor near the read, so it falls outside any sane radius. The correct axis is
   **splice-site proximity to the read's 5' boundary**, leaving the *other* site free.

**Fixes (status as of 2026-05-24):**
- DONE (commit 8e8dc8c): peel candidate-narrowing + depth cap. **Eliminated the OOM**
  (full-pool RSS 64 GB → ~14 GB, verified by live /proc probe).
- DONE (uncommitted): body-level `_nearby_junctions` narrowing (union of the loops' own
  cheap gates, computed once) + `_POOL_SEARCH_RADIUS` 10k→2k. Correctness-clean
  (246 focused tests). **Did NOT fix wall-time**: DIAGF (full-pool, fixed code) still
  stalled at 0/126 regions after 9.5 min — because each surviving nearby candidate still
  triggers the ~341-call ambiguity block, and the full pool leaves more nearby candidates
  than GFF-only.
- PENDING (Kevin's design, the real fix): **dual-site index + either-site fetch.** Index
  each junction by BOTH splice sites; fetch any junction with *either* site within the
  peeled-back genomic span + ~10 bp buffer (~60 bp window), NOT a radius on one site. This
  is intron-length-independent (human-ready) AND tighter (far fewer nearby candidates →
  far fewer 341-blocks). The `correct_read_3prime` fetch is the ONLY consumer that reads
  the index internals, so the restructure is localized to `_build_pool_chrom_index` +
  that fetch.

**Note:** the advisor's standing guidance is *don't* touch the shift×offset ambiguity
windows themselves (intentional splice-site-slide search). Reduce the candidate *count*
that enters them.

---

## Comprehensive pipeline audit — 2026-05-24 static pass

**Scope:** static scan of the main run-all pipeline, aligner wrappers, correct command,
corrected-consensus merge/BAM writers, split/chunk helpers, and manifest analysis. No new
cluster-scale timing run was performed for this section; treat these as prioritized
output-necessity and over-computation findings that still need the Verification section
before any code fix is declared done.

### Priority findings

1. **P0: Fresh correct-first `run-all` pays for a raw alignment consensus that it normally
   ignores.** `_run_alignment` calls `run_align(..., no_consensus=False)`, so `align` runs
   pre-correction consensus selection, writes `<sample>.rectified.bam`, writes raw
   consensus stats/report, and runs `samtools calmd`. Correct-first `run-all` then uses
   the per-aligner BAMs for per-aligner correction and builds `corrected_consensus.bam`
   from corrected TSVs plus raw aligner BAMs. The raw rectified BAM is only a fallback if
   per-aligner correction is unavailable, or a reuse target for `--skip-alignment`.
   **Necessity:** optional/debug/fallback, not required for the normal correct-first
   pipeline. **Fix direction:** let correct-first callers request align-only
   multi-aligner output (or gate raw consensus behind `--emit-raw-consensus` /
   diagnostics), while preserving fallback behavior when per-aligner correction is not
   requested or fails.

2. **P0: Manifest analysis defaults still emit convenience outputs that violate the
   documented memory contract.** `docs/architecture/multi_sample_pipeline.md` says
   manifest mode skips bedgraph and genomic-distribution outputs, but
   `_run_analysis_manifest` passes `no_bedgraph=False` and
   `no_genomic_distribution=False`. The manifest analyzer then accumulates bedgraph
   per-condition counts and, for genomic distribution, reads per-sample TSV columns into
   DataFrames and concatenates per condition. **Necessity:** bedgraphs and genomic
   distribution plots are QC/convenience outputs; DESeq2/PCA/GO/motif/cluster tables do
   not require them. **Fix direction:** make these opt-in for production manifest runs, or
   rework genomic distribution into true streaming aggregation before keeping it default.

3. **P1: Per-aligner correction rebuilds constant-per-sample junction context.**
   `_run_correction_per_aligner` invokes `correct` once per aligner; each process can
   reload annotation, build gene interval trees, load annotated junctions, and build/load
   the same cross-aligner junction pool. `prescan` already materializes
   `rescue_scan.pkl` and `junction_pool.pkl` for chunked pipelines, but fresh run-all
   paths do not consistently create/reuse those caches once per sample. **Necessity:**
   the junction-pool/rescue-scan caches are valid internal pipeline artifacts and are
   essential for chunked runs; they should also be the single source for per-aligner
   correction inside run-all. **Fix direction:** auto-prescan or share a per-sample cache
   before the per-aligner correction loop.

4. **P1: Position-index writers are duplicated and schema-drifted.** The canonical
   `write_position_index` emits `chrom, corrected_3prime, strand, count, count_ag_rich`,
   but streaming correction paths hand-write `*_index.bed.gz` with only
   `chrom, corrected_3prime, strand, count`. Manifest loaders tolerate missing
   `count_ag_rich` by defaulting to 0, silently losing AG-rich cluster annotation for the
   streaming path. **Necessity:** `*_index.bed.gz` is required for manifest-scale analysis;
   duplicate writers are not. **Fix direction:** route streaming writers through the
   canonical writer, or track `(total, ag_rich)` in streaming accumulators and keep one
   schema.

5. **P1: The 3'SS pool lookup still needs the pending dual-site index.** The current WIP
   narrows the full pool and reduces `_POOL_SEARCH_RADIUS`, but the index is still keyed
   on one intron coordinate. That is yeast-safe only by accident: human introns can put
   the opposite splice site tens or hundreds of kb away. **Necessity:** the pool itself is
   required for rescue, but the current index shape is a scaling/correctness hazard for
   human. **Fix direction:** index both splice sites and fetch junctions with either site
   within the peeled-back genomic span plus the biological buffer.

6. **P1: Raw-consensus 5' softclip rescue scans chromosome-wide annotated junctions.**
   `select_best_alignment` collects same-chromosome annotated junctions, then
   `_rescue_5prime_softclip` performs per-read distance checks and edit-distance work.
   This is small in yeast but has the same structure as the 3'SS bug for human-scale
   annotations. **Necessity:** no extra output; pure compute risk. **Fix direction:** pass
   a per-chrom/per-strand splice-site index or query a nearby-site index instead of
   handing every same-chromosome junction to every read.

7. **P2: Corrected-consensus diagnostics are emitted as if they were core outputs.**
   `merge_corrected_tsvs` loads all per-aligner TSVs into DataFrames, writes
   `comparison_summary.tsv` when requested by run-all, and always writes the aligner
   UpSet plot. Run-all also writes `cat5_candidates.tsv`. These are useful review
   artifacts, but no downstream core step consumes them. **Necessity:** optional
   diagnostics. **Fix direction:** gate behind a diagnostics flag or profile level; keep
   core merged TSV/manifest and final corrected-consensus BAM unconditional.

8. **P2: `align` generates a dead duplicate junction BED.** `align_command.py` writes
   `<prefix>_junctions.bed` when annotation is present, but minimap2 independently uses
   `annotation.junc.bed` through `get_minimap2_junc_args`, and the local
   `junc_bed_path` is not consumed. **Necessity:** dead output in the current align path.
   **Fix direction:** remove it or only emit it when an explicit export/debug flag asks
   for it.

9. **P2: gapmm2 leaves a PAF conversion sidecar.** `run_gapmm2` writes
   `<output>.paf` as an intermediate before BAM conversion, with no downstream consumer.
   **Necessity:** debug artifact only after BAM conversion succeeds. **Fix direction:**
   delete it by default, preserve it under `--keep-sam`/diagnostics or a dedicated
   `--keep-intermediates` flag.

10. **P2: `correct --emit-manifest` still routes through a monolithic TSV.** The command
    writes `corrected_reads.tsv`, hashes/counts the whole file, then renames it to
    `corrected_reads.region_000.tsv` unless `--emit-merged-tsv` is requested. The region
    checkpoint machinery already knows how to write region TSVs. **Necessity:** the
    manifest is necessary; a single-entry manifest around a full TSV provides little
    sharding value and adds a full-file read/hash pass. **Fix direction:** promote true
    per-region manifest emission from streaming/checkpoint workers and reserve the
    monolithic TSV for small or legacy runs.

11. **P3: Non-streaming correction buffers all per-read results in memory.**
    `process_bam_file_parallel` accumulates all region results before writing the TSV and
    index. This is acceptable for small tests and maybe yeast subsets, but production
    run-all should prefer streaming correction by default or switch automatically above a
    BAM/region threshold. **Necessity:** the final per-read TSV is core output; retaining
    every row in RAM before writing is not.

### Output necessity ledger

| Artifact / output | Needed for core pipeline? | Notes / action |
| --- | --- | --- |
| Per-aligner raw BAMs (`<sample>.<aligner>.bam` + `.bai`) | Yes, until finalization | Required for Module 2H pool/refinement, lazy HP scoring, and `corrected_consensus.bam` surgery. Can be archived or optionally removed only after final TSV/BAM are built and audit needs are satisfied. |
| Raw alignment consensus (`<sample>.rectified.bam`, raw stats TSV, raw HTML report, raw `calmd`) | No for normal correct-first run-all | Keep as fallback/debug/legacy output; skip by default when per-aligner correction and corrected consensus will run. |
| Per-aligner corrected TSVs / manifests | Yes, until merge/final BAM | Inputs to corrected-consensus merge and BAM surgery. Keep for reproducibility/debug, but do not require per-aligner corrected BAMs for the lazy final-BAM path. |
| Per-aligner corrected BAMs | Optional | Useful for IGV/debug; the corrected-consensus BAM writer can use raw BAMs plus correction TSVs. Gate for production if storage/runtime matters. |
| Final `corrected_reads.tsv` or `corrected_reads.manifest.tsv` + shards | Yes | Primary corrected-read table. Prefer manifest/shards for scale; keep monolithic TSV for small/legacy outputs. |
| `corrected_reads_index.bed.gz` | Yes for manifest mode | Needed for cluster/count loading without per-read TSV scans. Must use the canonical AG-rich-aware schema in every correction mode. |
| `corrected_consensus.bam` + `.bai` | Yes for correct-first run-all | Final inspectable consensus BAM after correction. |
| `junction_pool.pkl` / `rescue_scan.pkl` | Internal yes | Required for chunked/prescan and should become a reusable per-sample cache in run-all. Not a final user-facing result. |
| `<prefix>_junctions.bed` from `align_command.py` | No | Redundant/dead next to minimap2's `annotation.junc.bed`; remove or gate. |
| gapmm2 `.paf` intermediate | No after BAM conversion | Delete by default unless intermediate retention is requested. |
| `comparison_summary.tsv`, aligner UpSet PNG, `cat5_candidates.tsv` | Optional | Diagnostics only; gate behind an emit-diagnostics/profile flag. |
| Manifest bedgraphs and genomic-distribution plots | Optional | QC/convenience outputs; make opt-in or stream-aggregate before enabling at production scale. |
| Poly(A) trim metadata parquet | Yes for DRS Step 4 | Needed to restore trimmed sequence/quality context. |
| Stage provenance/checkpoint sidecars | Yes for resume/audit | Keep, but avoid sidecars that duplicate already-recorded information or force full-file reads. |

### Static scan trail

Files/functions reviewed in this pass: `core/commands/run/stages.py`,
`core/commands/run/single_sample.py`, `core/commands/run/multi_sample.py`,
`core/commands/align_command.py`, `core/align/multi_aligner.py`,
`core/consensus/select.py`, `core/consensus/scoring.py`,
`core/consensus/corrected_consensus.py`, `core/commands/correct_command.py`,
`core/bam/parallel.py`, `core/position_index.py`,
`core/analysis/manifest.py`, `core/analysis/loaders.py`, `core/commands/prescan_command.py`,
`core/commands/split_command.py`, and chunked/split generated script paths.

Useful greps from this pass: `no_consensus`, `junc_bed_path`, `comparison_summary`,
`cat5_candidates`, `write_position_index`, `count_ag_rich`, `no_bedgraph`,
`no_genomic_distribution`, `.paf`, `candidate_junctions`, and `annotated_junctions`.

---

## Anti-patterns to hunt for elsewhere

1. **Per-item scan of a global collection** — `for x in <big_set>:` inside a per-read /
   per-record function, where a cheap gate (`if dist > ...: continue`, containment) rejects
   most iterations. Fix: pre-narrow ONCE via a sorted/bisect index, interval tree, or a
   single-pass filter, then loop the small subset.
2. **Nested windows multiplying an expensive op.** An O(n×m) DP / local alignment / edit
   distance inside ≥2 nested loops (shift × offset × candidate). The multiplier hides the
   real cost; cut the outer count.
3. **Vestigial / oversized constants.** A magic number whose comment cites a default that
   no longer exists (the 10k-for-5000 case). Grep magic constants; re-derive against the
   *current* defaults.
4. **Bound keyed on the wrong axis.** A radius/window that scales with the wrong quantity
   (intron length instead of splice-site proximity). Ask: does this break for human / long
   features? If the bound is on coordinate A but the biology cares about coordinate B,
   it's wrong.
5. **Recompute-per-item what's constant-per-read/region.** Hoist read-level work out of
   per-candidate loops (RECTIFY already did this once — the `_leading_del` flag hoist,
   "Bug 3 patch perf fix" — look for more).

---

## How to hunt (methodology — profile, don't guess)

We theorized the cause twice and were wrong both times; only the sampling profile settled
it. Always measure.

- **py-spy (installed in the Sherlock rectify env).** `gdb py-bt` is UNUSABLE on this
  conda python3.9 (DWARF v5 vs gdb v2-4). Recipe:
  - Launch the target single-threaded (`--legacy-single-threaded -j 1`) in the background,
    `sleep 25` to pass setup, then
    `py-spy record -f raw -o stacks.folded --pid $PID --duration 300 --rate 100 --nonblocking`.
  - Aggregate by leaf frame:
    `awk '{c=$NF; s=$0; sub(/ [0-9]+$/,"",s); n=split(s,a,";"); sum[a[n]]+=c} END{for(k in sum) print sum[k]"\t"k}' stacks.folded | sort -rn | head`.
  - Split by caller: grep the folded lines for the suspect frame, sum counts with/without
    a calling-frame substring.
- **faulthandler for hangs:** `timeout -s ABRT <sec> python -X faulthandler -m rectify.cli ...`
  prints the current stack on the SIGABRT. Good for "where is it stuck right now."
- **Reproduce at scale on a dense input.** Synthetic unit tests miss this entirely (they
  use tiny pools). Use a real high-coverage region. Diag artifacts (subset BAM, 5 aligner
  BAMs, py-spy folded stacks, sbatch templates) live in
  `$SCRATCH/rectify_runall_diag_20260524/` and the trimmed/aligned BAMs in
  `$SCRATCH/rectify_wt_by4742_rep1_25846844_0/` on Sherlock.
- **grep targets to start from:**
  - `for .* in candidate_`, `for .* in .*junctions`, `for .* in .*pool`
  - nested `for _shift` / `for _off`; call sites of `_hp_edit_distance`,
    `align_clip_to_exon`, anything named `*edit_distance*` / `*align*`
  - `_SEARCH_RADIUS`, large integer literals as window/radius bounds
  - `.get_aligned_pairs()` or genome slicing inside per-candidate loops
  - intervaltree `.search` / `.overlap` / `search_point` called per read

---

## Suspect hotspots to audit next (from the DIAGE profile + code structure)

These appeared in the profile (minor here, but per-read and worth checking at other scales
/ organisms) or are structurally similar to the 3'SS case:

- **`walkback_3prime_guarded` / `walkback_drs_full`** (`core/correct/walkback.py`) — showed
  up in the DIAGE profile; per-read poly-A walkback.
- **gene-attribution interval-tree** — `search_point (intervaltree/node.py)` appeared in
  the profile; per-read interval-tree lookups for gene assignment (`gene_interval_trees`).
- **`indel_corrector`** (`detect_indel_artifacts`, `extract_insertions/deletions`,
  `utils/alignment.py`) — per-read, in the profile tail.
- **`junction_refiner`** (Module 2H) — builds + scans the junction pool; already partly
  bisect-optimized, but re-audit the per-N-op candidate search for the same radius-vs-axis
  issue, especially for human.
- **variant-aware rescue scan** and **per-aligner consensus scoring** — any path that
  receives `annotated_junctions` / the pool in full and loops it per read.

---

## Verification discipline (REQUIRE before declaring any perf fix "done")

A passing unit test is NOT evidence a perf fix works — the tests use tiny inputs. Require:
1. Focused blast-radius tests green, then `pytest -m "not slow"`.
2. **Scale re-run COMPLETES** on the dense representative input AND peak memory stays
   bounded (target <32 GB for the 5% subset).
3. **Measured** wall-time / call-count drop (e.g. count `_hp_edit_distance` calls before vs
   after; must drop ~10× in a dense region) — not "it feels faster."
4. **Output unchanged** — rescue counts / corrected-TSV within tolerance of the pre-fix
   run (no over-narrowing that silently drops real rescues).

---

## Pointers

- `AGENT_FIXES.md` → `[2026-05-24] ... 3'SS rescue` entry (full diagnosis + commands).
- Commit `8e8dc8c` (peel fix, OOM eliminated). Task #10 (baseline narrowing, in progress).
- Diagnostic jobs A–F + py-spy folded stacks: `$SCRATCH/rectify_runall_diag_20260524/`.
