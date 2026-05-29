# Handoff — upf1Δ validation review PNGs (task #10) · 2026-05-29

Branch `drs-validation-rebuild`. Origin tip after this session: **`e034ca7`**
(pushed). The upf1Δ validation set is built, green, and committed; the only
remaining piece is the **per-category review PNGs**.

## 1. What was done

- **upf1Δ parallel validation set — 32 reads, 8 categories × 4 (2+/2−)** (`de916fc`).
  `validation_reads_upf1d.bam` + `aligners/validation_reads_upf1d.{minimap2,gapmm2,mapPacBio,deSALT,uLTRA}.bam`
  + dorado-source archive + `VALIDATION_READS_upf1d.md` + `tests/test_validation_reads_upf1d.py`
  + `scripts/validation_data/upf1d_2026_05/{build_upf1d_validation.py,manifest.tsv}`.
  Categories cat1,2,3,5,6,7,8,9. cat4 omitted by design (trimmed panel yields zero
  `false_junction` artifacts — needs the untrimmed tail; not a real-pipeline scenario).
- **Cross-session WIP preservation** (`e034ca7`): committed the concurrent human
  track's uncommitted aligner docs / perf audit / branding / `dev/validation_review/`
  + the assistant's mapPacBio multi-primary note and over-call/ED spec. Excluded
  `dev/_*.py`/`_bb_*.qsub` scratch and build intermediates.

## 2. What's verified

- `pytest tests/test_validation_reads_upf1d.py -q` → **77 passed, 14 warnings, 412 s**
  (runs both no-Module-2H and Module-2H correct-first pipelines on M1).
- Bundle integrity (pysam): 32 primary reads, 8 categories × 4, **16+/16−**, all XV
  labels unique, cat5/cat6 stored alignments carry an N-op (Xg≥2), zero strand
  mismatches vs label.
- Ground truth biologically verified: cat3/5/6/9 land on annotated GT-AG introns;
  cat7 motifs are non-canonical & >50 bp from any annotated boundary; cat8 corrected
  3′ ends sit on a non-A genomic base.
- Push: `git push origin drs-validation-rebuild` → `726ebe1..e034ca7`, fast-forward
  (origin had not diverged).
- **NOT VERIFIED:** the review PNGs (not yet generated — task #10, below). No
  `rectified_upf1d/` dir exists yet.

## 3. Open items

- **Review PNGs (task #10) — DEFERRED, the subject of this handoff.** Why not done
  in-session: the renderer needs a `rectified_upf1d/` artifact dir (corrected +
  pA-restore BAMs + bedgraphs + per-aligner summary) that doesn't exist yet, and the
  batch driver `scripts/validation_data/generate_review_report.py` **hardcodes the WT
  bundle paths** (`VAL_DIR/rectified/...`, lines 57-61) — it must be parameterized for
  the upf1d bundle. Resume block has the plan.
- **Over-call rescue + raw-ED metric** — spec at `dev/specs/SPEC_overcall_rescue_and_ed_metric_20260529.md`,
  flagged for a **dedicated session** (touches `rectify/core/correct/` production code).
  Note the appended INVESTIGATION FINDINGS: the flat S=X cost is display-only;
  production `score_alignment` already defeats the cat1_minus_1 gaming pattern, so
  most of Item 1's assumed work is unnecessary — Item 2 (over-call) needs error-table
  insertion profiling first. Do NOT fold into validation-data commits.
- **`docs/aligners/mapPacBio.md` entanglement** — committed in `e034ca7` with the
  assistant's multi-primary note + the human track's WIP together (couldn't split
  without interactive `add -p`). If the human session resurfaces, reconcile there.

## 4. Resume command

**Resume:** generate the upf1Δ review PNGs. Bundle is committed at `e034ca7`;
per-aligner BAMs live at `rectify/data/validation/aligners/validation_reads_upf1d.*.bam`.

1. **Produce `rectify/data/validation/rectified_upf1d/`** (the analog of the WT
   `rectified/` dir). The renderer consumes: `rectified_pA_tail_trimmed.bam`
   (corrected), `rectified_pA_tail_soft_clipped.bam` (pA-restored),
   `corrected_3ends.plus.bedgraph`, `corrected_3ends.minus.bedgraph`,
   `per_aligner_summary.tsv`, `corrected_reads.tsv`.
   - Cleanest: `rectify run-all --drs rectify/data/validation/validation_reads_upf1d.bam --Scer --aligner-bams <the 5 aligners/ BAMs> -o rectify/data/validation/rectified_upf1d/` — **BUT** watch the documented 0.9.x manifest-hang (memory `project_rectify_release_strategy`); if it hangs, fall back to the per-aligner correct-first pipeline (mirror `_run_correct_first_pipeline` in the test) to make the corrected BAM, then `rectify analyze` to emit bedgraphs + the pA-restore BAM.
   - Run on M1 if swap is healthy (`sysctl vm.swapusage`), else Sherlock (push the
     bundle aligner BAMs first; splice fix already deployed there).
2. **Parameterize the driver** `scripts/validation_data/generate_review_report.py`
   for the upf1d bundle (it hardcodes WT `VAL_DIR/rectified/...` at lines 57-61 and
   `validation_reads.bam`). Add a `--bundle upf1d` switch (or a sibling script) that
   points `BAM_XV_MAP`→`validation_reads_upf1d.bam`, the `rectified/`→`rectified_upf1d/`,
   and the aligner-BAM glob → `validation_reads_upf1d.*`.
3. **Render per category:** `python scripts/validation_data/generate_review_report.py
   --per-category` (with the upf1d switch) → emits one PNG per XV label into
   `dev/validation_review/validation_read_review/cat*_upf1d_review_pngs/` (mirror the
   WT `cat*_review_pngs/` layout). Per-read coords come from
   `rectified_upf1d/corrected_reads.tsv` (chrom/start/end/orig_3p/corr_3p/five_prime).
4. **Sanity-check 2-3 PNGs** against `VALIDATION_READS_upf1d.md` (e.g. cat3_minus_1
   chrIV rescued 5′ = 307765; cat7_plus_1 chrII junction 246695-246734) and surface
   them to the user.
5. Commit surgically (`git add` explicit PNG dirs + `rectified_upf1d/` + the driver
   change) and push. Mark task #10 complete.

If `run-all` hangs >10 min on the manifest stage, kill it and use the correct-first
fallback — do not wait it out.

## 5. Files touched (this session)

Committed (`de916fc`):
- `rectify/data/validation/validation_reads_upf1d.bam` (+`.bai`), `_dorado_source.bam`
- `rectify/data/validation/aligners/validation_reads_upf1d.{5 aligners}.bam` (+`.bai`)
- `rectify/data/validation/VALIDATION_READS_upf1d.md` — new reference doc
- `tests/test_validation_reads_upf1d.py` — new, 77 tests
- `scripts/validation_data/upf1d_2026_05/build_upf1d_validation.py`, `manifest.tsv`

Committed (`e034ca7`): cross-session WIP — see commit body. Includes
`dev/specs/SPEC_overcall_rescue_and_ed_metric_20260529.md` and the
`docs/aligners/*` edits.

`[uncommitted]` in working tree (intentionally not committed):
- `scripts/validation_data/upf1d_2026_05/stage/` — pulled build-input BAMs (reproducible)
- `dev/_*.py`, `dev/_bb_*.qsub`, `README_preview.html` — human-track scratch / artifacts

Sherlock scratch (not in repo): `$SCRATCH/upf1d_fullpanel_20260527/` — full 5-aligner
panel, harvest correct outputs, scan/curate/subset scripts, `build32/` (the staged
32-read subset). Reusable if PNGs need re-derivation.
