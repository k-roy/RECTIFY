# WIP Disposition — 2026-05-20

Purpose: stop future agents from treating the dirty tree as one undifferentiated
blob. This file groups the current WIP into reviewable batches and separates
source/docs worth committing from generated artifacts that should be archived or
ignored.

Current branch at triage time: `drs-validation-rebuild`, head `a4e2933`.

Hard rule: do not use `git add -A` or `git add .`. Stage explicit paths only.

---

## Summary

The dirty tree is not one workstream. It is at least seven:

1. Documentation refresh and generated figures.
2. Architecture / aligner / protocol documentation.
3. Validation bundle rebuild and pA-rest renderer work.
4. pA-tail restore script rename/migration.
5. QNAME diagnostics scratch output.
6. Calibration run outputs.
7. Handoff / review / local notes.

Largest generated directories:

| Path | Size | Recommendation |
|---|---:|---|
| `scripts/validation_data/rebuild_2026_05/` | 306 MB | archive or ignore; commit only selected scripts/reports if needed |
| `scripts/diagnostics/` | 58 MB | commit reports/scripts only; ignore BAM/PAF/FASTQ scratch outputs |
| `scripts/calibration/runs/` | 3.6 MB | archive or ignore unless these are intended bundled calibration artifacts |
| `_h2_rebaseline/` | 844 KB | archive or ignore; not source |
| `rectify/data/validation/rectified/` untracked additions | 1.4 MB | decide whether these are bundled fixtures |

---

## Batch A — Operational Docs, Commit Separately

Paths:

- `CLAUDE.md`
- `HANDOFF.md`
- `TODO.md`
- `CODEX_AUDIT.md`
- `dev/specs/briefings/read_num_sidecar_briefing.md`

Recommendation:

- Review manually before commit because these steer agent behavior.
- Keep `CLAUDE.md` separate: it is a large rewrite (`~1077` changed lines).
- Keep `HANDOFF.md` separate or regenerate it last, after all cleanup decisions.

Suggested commits:

```bash
git add CLAUDE.md
git commit -m "docs(agent): refresh RECTIFY operating guidance"

git add TODO.md CODEX_AUDIT.md dev/specs/briefings/read_num_sidecar_briefing.md
git commit -m "docs(dev): record active audit and sidecar follow-ups"
```

Do not commit `HANDOFF.md` until the end unless Kevin wants the current handoff
as-is.

---

## Batch B — Architecture / Protocol Docs, Good Commit Candidate

Untracked source docs:

- `docs/aligners/gapmm2.md`
- `docs/aligners/minimap2.md`
- `docs/aligners/ultra.md`
- `docs/architecture/multi_sample_pipeline.md`
- `docs/architecture/pipeline_and_io.md`
- `docs/protocols/dt_primed_cdna.md`
- `docs/penalty_tables_quickref.md`
- `docs/audit_history.md`
- `dev/reviews/run_all_path_divergence_handoff_20260519.md`
- `dev/reviews/scratch_staging_review.md`
- `dev/reviews/scratch_staging_review_findings_20260519.md`

Recommendation:

- These are human-readable docs and should be easy to review.
- Commit before the README/figure rewrite so docs history stays readable.

Suggested commit:

```bash
git add \
  docs/aligners/gapmm2.md \
  docs/aligners/minimap2.md \
  docs/aligners/ultra.md \
  docs/architecture/multi_sample_pipeline.md \
  docs/architecture/pipeline_and_io.md \
  docs/protocols/dt_primed_cdna.md \
  docs/penalty_tables_quickref.md \
  docs/audit_history.md \
  dev/reviews/run_all_path_divergence_handoff_20260519.md \
  dev/reviews/scratch_staging_review.md \
  dev/reviews/scratch_staging_review_findings_20260519.md
git commit -m "docs: add architecture, protocol, and aligner notes"
```

Validation:

- `python -m pytest -q tests/test_run_command_wiring.py` is optional; docs only.
- Link/path sanity check with `rg '\]\([^)]*\)' docs/architecture docs/aligners docs/protocols`.

---

## Batch C — README + Figure Refresh, Review As One Unit

Tracked modified:

- `README.md`
- `docs/ALIGNER_RECOMMENDATIONS.md`
- `docs/troubleshooting.md`
- `docs/figures/*.png`
- `docs/figures/*.svg`
- deleted `docs/images/*.png`

Untracked figure generation scripts:

- `docs/figures/generate_5prime_junction_v3.py`
- `docs/figures/generate_adaptive_clustering_v3.py`
- `docs/figures/generate_cdna_isoform_v3.py`
- `docs/figures/generate_cdna_umi_v3.py`
- `docs/figures/generate_multi_aligner_v3.py`
- `docs/figures/generate_pipeline_overview_v3.py`
- `docs/figures/generate_polya_pretrim_v3.py`
- `docs/figures/generate_splice_classification_v3.py`
- `docs/figures/generate_walkback_v3.py`

Untracked backups:

- `docs/figures/*.v2-backup.png`
- `docs/figures/*.v2-backup.svg`

Recommendation:

- Commit the README/docs/figure refresh and generation scripts together.
- Do not commit `*.v2-backup.*`; archive outside the repo or add ignore rules.
- The deleted `docs/images/*.png` look intentional if README no longer links
  them. Verify with `rg 'docs/images|5prime_distribution|softclip_rescue'`.

Suggested commit:

```bash
git add README.md docs/ALIGNER_RECOMMENDATIONS.md docs/troubleshooting.md
git add docs/figures/*.png docs/figures/*.svg docs/figures/generate_*_v3.py
git add docs/images/5prime_distribution_pie_by4742.png \
        docs/images/genomic_distribution_pie_by4742.png \
        docs/images/softclip_rescue.png \
        docs/images/transcript_body_orf_distribution_by4742.png \
        docs/images/transcript_body_orf_distribution_set1_grid.png
git commit -m "docs: refresh README and protocol figures"
```

Suggested ignore addition:

```gitignore
docs/figures/*.v2-backup.png
docs/figures/*.v2-backup.svg
```

---

## Batch D — pA-tail Restore Script Migration, Source Commit Candidate

Tracked modified/deleted:

- `scripts/validation_data/generate_review_report.py`
- `scripts/validation_data/regen_pa_rest_bundle.py`
- `scripts/validation_data/render_read_alignment.py`
- deleted `scripts/validation_data/splice_polya_from_parquet.py`
- deleted `tests/test_splice_polya_from_parquet.py`

Untracked new:

- `scripts/validation_data/restore_polya_from_parquet.py`
- `tests/test_restore_polya_from_parquet.py`
- `scripts/validation_data/PLOTTING.md`
- `scripts/validation_data/audit_polya_trim.py`
- `scripts/validation_data/diagnose_cat3.py`
- `scripts/validation_data/diagnose_cat3_rescue_3ss.py`
- `scripts/validation_data/diff_fixture_vs_bundle.py`
- `scripts/validation_data/render_option_a_v2.py`

Recommendation:

- This is source work, but it is large (`render_read_alignment.py` changed by
  ~2293 lines). Review and test before committing.
- Keep it separate from regenerated BAM fixtures.

Suggested tests:

```bash
pytest -q tests/test_restore_polya_from_parquet.py
pytest -q tests/test_validation_reads.py -k 'polya or restore or render'
```

Suggested commit after tests:

```bash
git add scripts/validation_data/generate_review_report.py \
        scripts/validation_data/regen_pa_rest_bundle.py \
        scripts/validation_data/render_read_alignment.py \
        scripts/validation_data/restore_polya_from_parquet.py \
        scripts/validation_data/PLOTTING.md \
        scripts/validation_data/audit_polya_trim.py \
        scripts/validation_data/diagnose_cat3.py \
        scripts/validation_data/diagnose_cat3_rescue_3ss.py \
        scripts/validation_data/diff_fixture_vs_bundle.py \
        scripts/validation_data/render_option_a_v2.py \
        tests/test_restore_polya_from_parquet.py \
        scripts/validation_data/splice_polya_from_parquet.py \
        tests/test_splice_polya_from_parquet.py
git commit -m "refactor(validation): rename pA-tail restore tooling"
```

---

## Batch E — Validation Bundle Rebuild, Needs Decision

Tracked modified:

- `rectify/data/validation/PROVENANCE.json`
- `rectify/data/validation/README.md`
- `rectify/data/validation/VALIDATION_READS.md`
- `rectify/data/validation/aligners/validation_reads.*.bam`
- `rectify/data/validation/aligners/validation_reads.*.bam.bai`

Untracked fixture-like outputs:

- `rectify/data/validation/corrected_reads.tsv.pre_regen`
- `rectify/data/validation/corrected_reads_index.bed.gz`
- `rectify/data/validation/corrected_reads_stats.tsv`
- `rectify/data/validation/rectified/per_aligner/*`
- `rectify/data/validation/rectified/*pre_regen`

Recommendation:

- Commit only if these are intentionally replacing the packaged validation
  bundle.
- Do not commit `*.pre_regen` backups.
- Decide whether `rectify/data/validation/rectified/per_aligner/*.bam` should
  become bundled fixtures. They are small enough, but they add maintenance cost.

Suggested tests before commit:

```bash
pytest -q tests/test_validation_reads.py
pytest -q tests/test_gapmm2_seq_restore.py tests/test_normalize_read_name.py
```

Suggested ignore additions if backups remain local:

```gitignore
rectify/data/validation/**/*.pre_regen
```

---

## Batch F — Diagnostics: Commit Reports/Scripts, Ignore Outputs

Directory:

- `scripts/diagnostics/qname_mutation_survey/` (~58 MB)

Likely commit:

- `EDGE_CASES.md`
- `REPORT.md`
- `analyze_qnames.py`
- `make_synthetic_fastq.py`
- `rebuild_results_tsv.py`
- `results.tsv`
- `test_cases.tsv`

Likely ignore/archive:

- `*.bam`
- `*.paf`
- `*.fastq`
- `*.sam`
- `*.stderr`
- `qnames_*.txt`
- `ref/`

Suggested ignore additions:

```gitignore
scripts/diagnostics/**/*.bam
scripts/diagnostics/**/*.bai
scripts/diagnostics/**/*.paf
scripts/diagnostics/**/*.fastq
scripts/diagnostics/**/*.sam
scripts/diagnostics/**/*.stderr
scripts/diagnostics/**/qnames_*.txt
scripts/diagnostics/**/ref/
```

---

## Batch G — Generated Runs / Scratch, Archive Or Ignore

Do not commit by default:

- `_h2_rebaseline/`
- `scripts/validation_data/rebuild_2026_05/`
- `scripts/calibration/runs/`

Possible exception:

- Selected scripts from `scripts/validation_data/rebuild_2026_05/`:
  - `attach_xv_xg_tags.py`
  - `build_combined_dorado_source.py`
  - `option_b_recover_body_bases.py`
  - `fix_run_imports.py`
  - `run_canonical_correct.py`

Recommendation:

- Move these generated directories outside the repo or add ignore rules.
- If preserving exact run outputs matters, tar them under an external archive
  path, not in the git tree.

Suggested ignore additions:

```gitignore
_h2_rebaseline/
scripts/validation_data/rebuild_2026_05/
scripts/calibration/runs/
```

---

## Batch H — Loose Notes / Handoffs, Needs Kevin Decision

Untracked:

- `README_KR_edits.md`
- `RECTIFY_SHERLOCK_HANDOFF_20260518.md`
- `RECTIFY_SHERLOCK_HANDOFF_20260518_v2.md`
- `docs/handoffs/HANDOFF_2026-05-19_validation_qsrev_cdna_pass.md`
- `docs/handoffs/cdna_umi_phase_d_resume.md`
- `docs/handoffs/cdna_umi_stratified_calibration.md`
- `handoffs/HANDOFF_2026-05-18_bam_provenance_h2_han2023.md`
- `handoffs/HANDOFF_2026-05-19_245f486.md`
- `handoffs/HANDOFF_2026-05-19_manifest_correction_divergence.md`
- `handoffs/HANDOFF_2026-05-20_ba61d31_commits_zero_a_b.md`
- `pA_tail_DRS_adapter_executive_summary.html`
- `pA_tail_DRS_citation_validation.md`

Recommendation:

- If these are useful historical handoffs, commit under `docs/handoffs/` or
  `handoffs/` consistently.
- If they are personal scratch notes, archive outside the repo.

---

## Proposed Cleanup Order

1. Commit Batch B docs.
2. Commit Batch F reports/scripts and add ignore rules for diagnostics outputs.
3. Commit Batch D pA-tail restore migration after focused tests.
4. Commit Batch C README/figures after link sanity check.
5. Decide Batch E validation bundle.
6. Add ignore rules for Batch G generated runs and backup artifacts.
7. Decide Batch A/H operational notes and handoffs.
8. Regenerate final `HANDOFF.md`.

After each commit:

```bash
git status --short
git push origin drs-validation-rebuild
ssh hoffman2 'bash -lc "cd /u/home/k/kevinroy/software/rectify && git pull --ff-only"'
ssh sherlock 'bash --norc --noprofile -c "cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify && git pull --ff-only"'
```

