# Handoff: `rectify run-all` manifest vs single-sample correction-path divergence

**Date:** 2026-05-19
**Author:** Claude (Opus 4.7), based on smoke testing run by Kevin Roy on H2
**Status:** FIXED — `_process_one_sample` now uses per-aligner correct + merge (2026-05-19)
**Fix commit:** pending — change is in `rectify/core/commands/run/single_sample.py`
**Remaining:** smoke re-run on H2 to verify manifest vs single-sample output now agrees
**Brief to the reader:** read this end-to-end before writing code. The diagnosis
is solid but the decision about which path is canonical is yours — there's
likely intent behind the two-path design that this doc doesn't fully recover.

---

## TL;DR

`rectify run-all` has two separate correction code paths that produce
**structurally different** corrected_reads.tsv output for the same input:

| Invocation | Correction path | Output schema | `winning_aligner` col |
|---|---|---|---|
| `rectify run-all <input> ...` (single positional arg) | per-aligner correct → `merge_corrected_tsvs` | 40 columns | **yes** |
| `rectify run-all --manifest manifest.tsv ...` | single `_run_correction` on the consensus rectified.bam | 38 columns | **no** |

On a 5 217-read DRS subsample of `wt_rep1.bam` (PRJNA1229592), the two paths
disagree on `corrected_3prime` for **488 / 4 589 shared reads (10.6 %)**,
with 31 % of disagreements >100 bp and 5 % mapping to different chromosomes
entirely. The divergence is not a bug in either path's internal logic — it
reflects a fundamental design difference that's not documented anywhere.

---

## The two code paths, with file:line citations

### Manifest path: `_process_one_sample` (used by `--manifest`)

`rectify/core/commands/run/single_sample.py:185-205`

```python
print(f"  [{sample_id}] Correcting 3' ends…", flush=True)
try:
    _run_correction(
        bam_path=bam_to_correct,        # = rectified.bam from consensus selection
        output_dir=_work,
        genome_path=genome_path,
        annotation_path=annotation_path,
        args=args,
    )
    ...
```

A **single** `_run_correction` call on the consensus rectified.bam. No
per-aligner correction. No merge step. Output has 38 columns; no
`winning_aligner` column because there's no winner selection happening
post-consensus — the consensus already chose at align-time.

### Single-sample path: `_run_single_sample` (used by single positional arg)

`rectify/core/commands/run/single_sample.py:486-528`

```python
if per_aligner_bams:
    from ...consensus.corrected_consensus import merge_corrected_tsvs, identify_cat5_candidates
    print(f"    Running per-aligner correction ({len(per_aligner_bams)} aligners)...")
    per_aligner_tsvs, per_aligner_corrected_bams = _run_correction_per_aligner(
        per_aligner_bams=per_aligner_bams,
        output_dir=work_dir,
        genome_path=genome_path,
        annotation_path=annotation_path,
        args=args,
    )
    if per_aligner_tsvs:
        ...
        corrected_tsv = merge_corrected_tsvs(
            per_aligner_tsvs=per_aligner_tsvs,
            output_tsv=work_dir / 'corrected_reads.tsv',
            summary_tsv=_summary_tsv,
            per_aligner_corrected_bams={
                a: str(p) for a, p in per_aligner_corrected_bams.items()
            } if per_aligner_corrected_bams else None,
            overhang_table=_overhang_table,
        )
```

Per-aligner correct on each of the N raw aligner BAMs, then HP-aware
edit-distance merging in `merge_corrected_tsvs`. Output has 40 columns
including `winning_aligner` and `sample`. Comment at L508-511 calls this
"the validated correct-first path per CLAUDE.md PIPELINE ORDER" — strongly
suggesting this path is the *intended* / canonical one.

### Supporting machinery

- `_run_correction_per_aligner`: defined at
  `rectify/core/commands/run/stages.py:294-416`.
- `_run_correction`: defined at `rectify/core/commands/run/stages.py:188-293`.
- `merge_corrected_tsvs`: defined at
  `rectify/core/consensus/corrected_consensus.py:577`, docstring documents
  HP-ED winner selection when `per_aligner_corrected_bams` is provided.

CLAUDE.md (project-level instructions) references "correct-first pipeline
order" — see `feedback_per_aligner_rescue_runs_first.md` memory and the
in-code comment at L508-511 — both reinforce that **per-aligner correction
is the documented intent**. The manifest path appears to predate this
shift and was never updated.

---

## Empirical evidence (smoke 13445930, H2, 2026-05-19)

### Schema diff

```
=== run-all (manifest mode) corrected_reads.tsv columns ===
[1..38]: read_id, chrom, strand, ..., reanchor_clip_len
(38 columns, no winning_aligner)

=== sequential (mimicking single-sample mode) corrected_reads.tsv columns ===
[1..38]: same as above, plus:
[39]: winning_aligner
[40]: sample
```

### Per-row corrected_3prime divergence

| Test config | --aligner-bams | corrected BAM artifact | Divergence |
|---|---|---|---|
| Sequential per-aligner correct, no Module 2H | ✗ | `--output-bam` (wrong artifact) | 488 / 4 589 = **10.6 %** |
| + `--aligner-bams` to each correct call | ✓ | `--output-bam` (wrong artifact) | 483 / 4 589 = **10.5 %** |
| + `--write-corrected-bam` (correct artifact for merge) | ✓ | `--write-corrected-bam` | 491 / 4 589 = **10.7 %** |

**The plumbing dimension (Module 2H, prescan, write-corrected-bam) explains
~0 % of the divergence.** All three configs give the same answer within
sampling noise. The divergence is an *invariant* property of comparing the
two code paths.

### Breakdown of disagreeing rows (483 of 4 589)

| Field | % differ in disagreeing rows |
|---|---:|
| `alignment_end` | 65.4 % |
| `correction_applied` | 64.3 % |
| `confidence` | 51.0 % |
| `five_prime_position` | 43.9 % |
| `alignment_start` | 43.4 % |
| `polya_length` | 43.0 % |
| `junctions` | 38.7 % |
| `n_junctions` | 32.4 % |
| **`chrom`** | **4.9 %** |
| **`strand`** | **2.5 %** |
| `five_prime_rescued` | 7.0 % |

### Shift-magnitude distribution

| Magnitude | Count |
|---|---:|
| `|delta| <= 2 bp` | 132 (27 %) |
| `|delta| <= 10 bp` | 240 (49 %) |
| `|delta| <= 100 bp` | 333 (69 %) |
| `|delta| > 100 bp` | 150 (31 %) |
| `|delta| > 10 kb` (different chromosome territory) | 24 (5 %) |

So divergence is bimodal: a tail of small slides (1-10 bp, junction-shift
territory) plus a tail of complete-aligner-swap cases (chromosome-flipping
reads where the two paths pick entirely different alignments).

### Sequential winning_aligner distribution

```
{'deSALT': 4070, 'mapPacBio': 343, 'minimap2': 199, 'gapmm2': 43}
```

deSALT wins ~87 % of reads under HP-ED on this dataset. The manifest path,
by working from the consensus rectified.bam (which selects per-alignment
from the same aligners but with a different scoring rubric), would silently
overwrite many of these deSALT choices.

---

## What I ruled out

Don't redo these — I already burned cluster time on them.

1. **Module 2H plumbing**: passing `--aligner-bams` to every per-aligner
   `rectify correct` did not move divergence (488 → 483).
2. **Variant prescan**: there's no separate prescan cache being passed in
   the run-all manifest mode either; both invocations use the in-process
   variant scan over the same input BAM.
3. **`--write-corrected-bam` vs `--output-bam` artifact mismatch**: I
   initially thought sequential was feeding the wrong BAM to
   `merge_corrected_tsvs`. Fixed it; divergence didn't move (483 → 491).
4. **Trimmed FASTQ ordering**: both invocations call `trim_drs_bam_polya`
   internally with identical args; no observed byte-level diff.
5. **uLTRA failure**: uLTRA produces no BAM in the run-all manifest mode
   but produces a real 2.9 MB BAM in fresh-scratch isolated invocations.
   This is a stale-cache symptom in the uLTRA tmp dir; uLTRA itself is
   working on H2. Document but not relevant to the manifest/single-sample
   divergence.
6. **deSALT identity**: H2 uses the conda binary
   (`~/.conda/envs/rectify/bin/deSALT`) which works; Sherlock previously
   used the vendored binary (SIGSEGV) but per the freshly-updated
   `reference_rectify_5aligner_status.md` memory entry, the conda symlink
   at `~/bin/deSALT` is now in place there too. deSALT identity is not
   the divergence source.

---

## Open architectural questions

These are the questions a fresh-eyes agent should answer before writing
code. Read git history (`git log --follow rectify/core/commands/run/`) to
recover intent.

1. **When did the two paths diverge?** Was the manifest path written first
   (when there was only single-BAM correct) and never updated when
   per-aligner correction landed? Or was per-aligner correction added
   *only* to the single-sample path on purpose (e.g., for benchmarking),
   and manifest mode is the production-intended path?

2. **Which path is canonical?** The in-code comment at
   `single_sample.py:508-511` calls per-aligner+merge "the validated
   correct-first path per CLAUDE.md PIPELINE ORDER" — but if that's
   canonical, why does the manifest path (the one used in production for
   multi-sample experiments) not do it?

3. **Is the schema gap intentional?** Should `winning_aligner` be present
   in manifest-mode output for downstream consumers, or is its absence
   load-bearing for something (e.g., DESeq2 step expecting a specific
   column set)?

4. **Test coverage gap.** Are there tests that pin either path's behavior?
   Grep `tests/` for `_process_one_sample`, `_run_correction_per_aligner`,
   `merge_corrected_tsvs`. If both paths have explicit tests, aligning
   them may break tests on purpose; if only one has tests, that hints at
   which is canonical.

5. **Production blast radius.** Look at recent production runs (set1/2/3
   on Sherlock per `project_rectify.md` memory). Did they use manifest
   or single-sample mode? If manifest, the winning_aligner-derived
   downstream analyses have been silently missing a column for the entire
   production history. If single-sample, the discrepancy is theoretical
   only.

---

## Recommended approach

1. **Investigate first, code second.**
   - `git log --follow rectify/core/commands/run/single_sample.py` —
     when did per-aligner+merge land? Was manifest mode considered?
   - `grep -rn "winning_aligner" rectify/ tests/` — who consumes this
     column? If downstream code depends on it, manifest-mode output is
     incomplete for those consumers.
   - Check `rectify/core/commands/run/multi_sample.py` to see how
     the manifest dispatch ties into combined analysis — does combined
     analysis read winning_aligner?

2. **Decide which path is canonical.** Best guess from the evidence:
   single-sample (per-aligner + merge) is canonical, manifest is a stale
   workflow. But the fresh-eyes agent should confirm this independently.

3. **Align the two paths.** If single-sample is canonical:
   - Refactor `_process_one_sample` to call the same
     `_run_correction_per_aligner` + `merge_corrected_tsvs` flow as
     `_run_single_sample` does.
   - This unblocks consistent `winning_aligner` for downstream analyses.
   - Watch out for the scratch-staging plumbing — the per-aligner BAMs
     produced by alignment must be available to the correction step,
     which requires careful ordering with the per-worker scratch cleanup.
     (My session already added that scratch staging; see
     `scratch_staging_review_findings_20260519.md`.)

4. **Add a regression test.** Two-sample manifest with a one-sample
   single-positional invocation on the same data should produce
   byte-identical per-sample `corrected_reads.tsv`.

5. **Doc updates.** Once the paths agree, document the consolidated
   pipeline in `docs/user_guide/commands/run.md` and remove the implicit
   manifest-vs-single-sample inconsistency from `docs/ARCHITECTURE.md`.

---

## Reproduction artifacts (on H2)

The smoke driver and outputs are preserved at:

- **Driver:**
  `/u/scratch/k/kevinroy/smoke_PRJNA1229592_subsampled/rectify_smoke.sh`
- **Subsampled inputs:**
  `/u/scratch/k/kevinroy/smoke_PRJNA1229592_subsampled/{wt,ysh1}_rep1.bam`
- **Last smoke log (job 13445930):**
  `/u/scratch/k/kevinroy/smoke_results/rectify_smoke.log`
- **run-all output:**
  `/u/scratch/k/kevinroy/smoke_results/runall/wt_rep1/corrected_reads.tsv`
  (38 cols, manifest mode, no winning_aligner)
- **Sequential output:**
  `/u/scratch/k/kevinroy/smoke_results/sequential/wt_rep1/correct/corrected_reads.tsv`
  (40 cols, per-aligner+merge, has winning_aligner)
- **Per-aligner corrected BAMs:**
  `/u/scratch/k/kevinroy/smoke_results/sequential/wt_rep1/correct/<aligner>/<aligner>.rectified_corrected_3end.bam`

The smoke takes ~10 min wall on `-pe shared 8 -l h_data=4G` and is
deterministic across reruns. Resubmit with `qsub rectify_smoke.sh`.

Canonical deposit used for the smoke:
`/u/project/guillom/shared/raw/PRJNA1229592_cpa_depletion_DRS_2025/`
(15 GB, README + sha256, read-only).

---

## Related findings & context from this session

- `rectify/dev/reviews/scratch_staging_review_findings_20260519.md` —
  the upstream scratch-staging review whose smoke test surfaced this
  divergence. The collision BLOCKER fix, scratch-on-HPC gating, rsync
  filter, Step 4 routing through scratch, and split BAM-input rework
  are all already landed.
- `docs/ALIGNER_RECOMMENDATIONS.md` — production panel (4 vs 5 aligners,
  deSALT status per platform) is being revised; the memory entry
  `reference_rectify_5aligner_status.md` was updated 2026-05-19 to
  reflect that deSALT is working everywhere via conda symlink.
- Memory entries to be aware of when reading the code:
  - `feedback_per_aligner_rescue_runs_first.md` — "correct-first" intent
  - `feedback_hp_edit_distance_semantics.md` — HP-ED scoring rules
  - `feedback_rectify_junction_slide.md` — Module 2H junction-slide
    handling
