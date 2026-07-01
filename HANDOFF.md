# HANDOFF — DRS validation bundle: winner-selection bug + coherent rebuild

**Date:** 2026-06-26 · **Branch:** `drs-validation-rebuild` (SHARED — a concurrent COMPASS session also
commits here; user confirmed you own the validation work). **You are a fresh session: read this top-to-
bottom, then continue.** Read `~/.claude/CLAUDE.md` + `~/work/CLAUDE.md` + `CLAUDE.md` first (rm-safety,
M1 8GB discipline, ControlMaster, surgical `git add` — never `-A`).

## ⚑ UPDATE 2026-06-29 (#5) — walkback homopolymer-undercall fix (CORE change, IN FLIGHT)

Triggered by cat2_plus_1 review: a real DRS long-homopolymer undercall was being eaten by the walkback.

**Done (uncommitted):**
- `rectify/core/correct/walkback.py`: large-deletion pre-scan (BOTH right-side ~line 600 and left-side
  ~line 815) now gated by a **homopolymer-undercall guard** — it will NOT skip past a ≥`large_del_min_bp`
  deletion if ≥`_MIN_GENOMIC_ANCHOR_3P` (=5) non-stop-base read=ref MATCHES sit 3' of it (those are genuine
  templated residues = a homopolymer undercall flanked by real genomic seq, not a force-aligned poly-A tail).
  New module constant `_MIN_GENOMIC_ANCHOR_3P = 5` near top (after `APPLIED_NONE`).
- `tests/test_validation_reads_upf1d.py`: cat1_minus_1 expected corrected_3prime `1162861 → 1162817`
  (same bug class; comment added). `test_indel_correction_applied` still passes (label still present).

**Verified:**
- ROOT CAUSE proven via `/private/tmp/.../scratchpad/wb_probe5.py` + `wb_probe7.py`: must call
  `_decode_eq_seq_inplace(read, genome)` BEFORE the walkback (uLTRA/minimap2 emit `=`-SEQ; undecoded → walkback
  sees all-mismatch → returns None, masking the bug). After decode: cat2_plus_1 minimap2/uLTRA walked 23759→23711
  (correction_bp=48 = 9D+39M, eating genuine `39=` over AT-rich genomic incl. matched T's). deSALT (no 9D)→None.
- POST-FIX: cat2_plus_1 minimap2/uLTRA → corrected 23758 (bp=1, KEEPS the 9D+match); deSALT unchanged.
  upf1d cat1_minus_1 → 1162817 (37 bp perfect C/G-rich match 3' of a 6D; old 1162861 over-walked it).
- Targeted suite (walkback/atract/indel/correct/polya): was 1 fail (the upf1d test, now updated).
- Advisor consulted twice (approach + reconcile); 1162861 confirmed a regression snapshot (git blame
  de916fcf, no orthogonal-truth comment), not external ground truth.

**Open / NOT YET DONE:**
- FULL `pytest -m "not slow"` running in background (task `bovc9ov6z`; result tail in
  `scratchpad/pytest_full.txt`) — CHECK RESULT. (Earlier runs `b3upqsr0d` were killed by M1 pressure.)
- FULL suite result: **1605 passed, upf1d now green**; 1 EXPECTED fail =
  `test_bam_parallel_state.py::test_process_bam_file_parallel_deterministic` (a GOLDEN-HASH test that pins
  correction output — my fix changes it intentionally; re-record via `RECTIFY_RECORD_GOLDEN=1` and paste into
  `GOLDEN_HASH_VALIDATION_MINIMAP2_NT2`). This is the expected signal, not a bug.
- 36-read blast-radius diff: **DONE** (regen task `bwyhclvz0` exit 0). **RESULT: 0 changes** to any merged
  winner/corr3 across all 36 reads (diff of regenerated `corrected_reads.tsv` vs
  `scratchpad/BEFORE_corrected_reads.tsv`). The fix changes ONLY per-aligner diagnostic BAMs: cat2_plus_1
  minimap2/gapmm2/uLTRA rectified rows now keep `9D38(=|M)9S` (the 9D), instead of clipping to `47S`.
  cat2_plus_1 winner UNCHANGED = deSALT @ 23754 (uLTRA now @ 23758 with clean 49=9D38= but does NOT win —
  a SEPARATE EER-ED/winner-selection question the user was asked about).
- Fixtures regenerated in the working tree (UNCOMMITTED, `git status` shows `rectified/` BAMs + CASE_STUDIES.md
  modified). Golden hash NOT yet re-recorded.

**EER-ED follow-up (user chose to pursue it FIRST, before committing the walkback fix):**
- Diagnosed why deSALT still wins cat2_plus_1 over uLTRA's clean 9D. Per-op ED breakdown (probe
  `scratchpad/ed2.py`, uses bundled `penalty_tables/penalty_scores.tsv`+`str_penalty_scores.tsv`):
  uLTRA/minimap2/gapmm2 ED=12.99 vs deSALT ED=9.90. **The 9D is correctly CHEAP = 0.078**
  (`del_cost(hp=24,'A')=0.0087`/base) — the user's intuition was right. **The blocker is the poly-A-tail
  CLIP penalty**: uLTRA honestly clips the 9-base tail (`9H`→9.0 flat in `_cigar_hp_edit_distance`,
  corrected_consensus.py line ~139-143), deSALT force-aligns the tail (insertions=6.8, no clip) and wins.
  Merge scores the `.trimmed` (HARD-clipped) per-aligner BAMs (regen passes `per_aligner_corrected_bams=
  per_aligner_trimmed`). uLTRA also has the MORE-correct 3' (23758 vs deSALT 23754; genomic match ends 23759).
- Clip penalty exists on purpose (stop walkback+clip from beating an inline deletion — cat1_minus_2), so
  the fix must NOT just remove it. Recommended: **tail-aware clip penalty** — subtract the known poly-A tail
  length (TSV `three_prime_soft_clip_length`/`polya_length`) from the clip cost; keep penalizing genuine
  genomic clips. Wrinkle: `_cigar_hp_edit_distance` sees the BAM read only (tail H-clipped → bases gone),
  so the tail length must be threaded in from the TSV.

**Tail reconcile (probe `scratchpad/tail_reconcile.py`):** uLTRA clips final 9 read bases `AAATAAAAT`;
genome past its endpoint (23759+) is `ACAATGAT` → the clip is genuine non-templated tail, NOT discarded
genomic. uLTRA cleaner alignment (3' 23758 at genomic-match boundary) but CPA genuinely ambiguous (tail is
AT-rich/degenerate). User calls it a toss-up.

**Two independent advisor agents** (neutral identical briefing, A=leave winner / B=tail-aware clip penalty):
- Advisor 2 (DONE): **Option A now, Option B as a separate validated project.** It IS a systematic flaw
  (flat clip vs near-free HP-del rewards force-aligning the tail), but the clip penalty is load-bearing
  across all reads; a tail-aware spec needs a robust poly-A/AT-tail detector (this tail is degenerate AT-rich
  = where naive detectors misfire) + full before/after winner audit + regression gates. Don't retune a core
  term on one ambiguous read where the "wrong" winner still lands on a defensible CPA.
- Advisor 1 (DONE): **Option A.** "A tidier alignment is not evidence of a more correct CPA when the CPA is
  ambiguous" — preferring uLTRA = confirmation bias; calibrating a load-bearing cross-cutting cost on one
  read with no known-correct answer = overfitting. B = flagged hypothesis for a truth-anchored study.
- BOTH advisors + assistant agree: **ship A now, defer B to a validated effort.**

**Option B refined (user's idea, 2026-06-29) — this is the deferred-B spec:** instead of binary tail-length
subtraction, make the 3' clip penalty GRADED by a probability the clipped tail is poly-A-derived (assess
A-richness; e.g. 9S=AAATAAAAT = 7A/9 → high P → discount the 9.0 toward ~0; low-A clip keeps full penalty).
Elegant: A/T-richness may AUTO-preserve the clip penalty's load-bearing role (genomic discards = low-A = full
penalty (cat1_minus_2); real tails = high-A = forgiven). Constraints to handle in the B effort:
  - STRAND-AWARE: A-rich on plus, T-rich on minus (key off per-strand stop_base, not literal A).
  - PLUMBING: merge scores `.trimmed` (HARD-clipped) BAMs → clipped bases absent; must score `.softclipped`
    variant OR carry clipped-tail composition via TSV. (`_cigar_hp_edit_distance` in corrected_consensus.py.)
  - CALIBRATION: the P→penalty curve can't be fit on cat2_plus_1 (unknowable CPA). Validate on reads with
    UNAMBIGUOUS CPA (clean poly-A, ideally orthogonal 3'-truth) + before/after 36-read winner reshuffle +
    regression gate (cat1_minus_2-class must still lose). cat2_plus_1→uLTRA should be a CONSEQUENCE, not target.

**DONE: walkback fix committed `fc44ee2`** (Option A shipped). Surgical add: walkback.py + 2 tests
(upf1d expectation + golden hash a41ec734…) + CASE_STUDIES.md cat2_plus_1 entry + regenerated `rectified/`.
Reverted dorado_source.bam (re-serialization churn, identical CIGARs). `pytest -m "not slow"` = 1606 passed.

**NOW STARTING — Option B effort (user approved 2026-06-29): graded A/T-richness clip penalty.**
User guidance: "test the graded penalty on a subset of our DRS datasets to see the blast radius" (i.e. go
beyond the 36 validation reads — use real DRS data for the winner-reshuffle audit). Spec lives in
CASE_STUDIES cat2_plus_1 "Open EER-ED hypothesis". Resume:
1. Implement graded clip cost in `_cigar_hp_edit_distance` (corrected_consensus.py ~line 139-143):
   3'-clip penalty scaled by P(tail is poly-A-derived) = strand-aware stop-base-richness of the clipped run
   (plus→A, minus→T). PLUMBING: merge scores `.trimmed` (hard-clipped) BAMs → clipped bases absent; score
   the `.softclipped` variant instead, or carry clipped-tail composition via the per-aligner TSV.
2. Find a DRS subset to test on (local validation bundle is only 36 reads; real DRS lives on cluster H2 —
   ask user for path or use an h2-qsub chunk). Run before/after winner reshuffle; REGRESSION GATE: the
   cat1_minus_2-class clip-to-win must still LOSE; cat2_plus_1→uLTRA must emerge as a consequence.
3. Calibrate the P→discount curve on UNAMBIGUOUS-CPA reads (clean poly-A), not cat2_plus_1.

**B prototype DONE (local sanity gate, non-invasive)** — `scratchpad/graded_proto.py` recomputes per-aligner
ED over the 36-read `.softclipped` BAMs with a graded trailing-clip cost = `len × (1 − stop_base_richness)`
(plus→A, minus→T). Results: **2/36 winners flip** (min-ED approx, NO merge tiebreakers — approximate!):
  - cat2_plus_1: deSALT(9.90) → minimap2/gapmm2/uLTRA (tie 12.99, clean 9D38M) ✓ target flips.
  - cat1_minus_2: uLTRA → uLTRA ✓ regression case holds.
  - 00a1e01e (chrVII, huge-ED messy read): mapPacBio→uLTRA — UNVETTED second flip, investigate.
Caveats: min-ED only (real merge differs), linear (1−r) curve uncalibrated, 36 reads = curated edge cases
not representative DRS.

**00a1e01e INVESTIGATED (false positive of naive metric):** its 76-bp clip `TACAAAATGAAAT…` is NOT poly-A
(whole-clip A=0.36, terminal-12 A=0.25) — a discarded complex/2nd-exon segment. mapPacBio force-aligns it as
insertion garbage (no clip). The naive whole-clip-average linear `(1−r)` over-discounts it → false flip.
**REFINED METRIC (key result):** assess stop-base richness in the **terminal window** (3'-most ~12 bp where
the tail sits), NOT whole-clip average. cat2_plus_1 terminal-12 A=0.78 vs 00a1e01e 0.25 → cleanly separates.
A hard ≥0.8 threshold FAILS (cat2's degenerate tail is 0.78) → must be a GRADED score on the terminal window.

**DRS audit data LOCATED (user approved wt_by4742_rep1):** per-aligner RAW BAMs on Sherlock scratch
`/scratch/users/kevinroy/rectify_single_22321416_0/wt_by4742_rep1.{minimap2,gapmm2,mapPacBio,deSALT,uLTRA}.bam`
(+ `.rectified.bam` merged). ⚠️ Sherlock scratch is EPHEMERAL (purge policy) — re-confirm/re-stage if gone
(raw source: `/oak/…/larsms/Users/kevinroy/projects/TRT/raw_data/nanopore/inhouse_by4742_dst1_4nqo_ski7/`,
but that dir has dst1d/ski2d/4nqo — wt_by4742_rep1 itself via the build_dorado fallback path / the scratch run).
Oak rectify install: `/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/` (may be STALE — needs the
`fc44ee2` walkback fix deployed before the audit, else baseline is pre-fix).

**AUDIT: Sherlock job `32183378`** (larsms; RESUBMIT of 32181800 which FAILED rc=3 "fetch on bamfile without
index" — fixed by sort+indexing raw BAMs first, commit `<sbatch fix>`). Sentinel + report paths use the CURRENT
job id: `/scratch/users/kevinroy/graded_clip_audit/.audit_rc` and `graded_clip_reshuffle_32183378.tsv`. M1
background watchers keep getting killed — rely on the sentinel + this check logic, not a live watcher.
CHECK: `ssh sherlock 'cat /scratch/users/kevinroy/graded_clip_audit/.audit_rc; sacct -j 32183378 -X -o State'`.
If rc=0 → fetch+read `graded_clip_reshuffle_32183378.tsv`. If rc≠0 → `tail slurm_32183378.err`.

**FIRST-PASS AUDIT RESULT (job 32183378, rc=0, 2026-06-30):** 1929 real DRS reads → **3 winners flip (0.16%),
ALL tail-flips (tail_share=1.0, minimap2↔uLTRA), 0 SUSPICIOUS, 0 regressions.** The graded terminal-window
clip penalty has a tiny, clean blast radius on real data — de-risks the Option-B core change. CAVEAT: stale
Oak rectify (no fc44ee2 guard, confirmed in log) → cat2_plus_1-class 9D-preservation synergy UNDERCOUNTED;
expect a few more legit tail-flips on a fix-deployed pass. Reshuffle TSV saved:
`/scratch/users/kevinroy/graded_clip_audit/graded_clip_reshuffle_32183378.tsv` (3 rows).

**NEXT STEPS for Option B:**
1. ✅ DONE — fix-deployed FRESH tree: `fc44ee2` code rsync'd to `/scratch/users/kevinroy/rectify_fc44ee2/`
   (176 py files, imports as `rectify`, guard=True; uses Oak data via absolute paths). fc44ee2 is NOT on
   origin (origin head 5b44b31 = shared/COMPASS) so do NOT push/clone — rsync code from M1 is the deploy path.
2. **IN FLIGHT — 2nd pass + calibration: Sherlock job `32193027`** (`run_graded_clip_audit_p2_sherlock.sh`,
   PYTHONPATH=fresh tree). Corrects once, then sweeps graded_clip_audit.py over tail_frac{0.5,0.6,0.7,0.8}
   × term_window{8,12,16}. Sentinel `/scratch/users/kevinroy/graded_clip_audit/.audit_p2_rc`.
   CHECK: `ssh sherlock 'cat /scratch/users/kevinroy/graded_clip_audit/.audit_p2_rc; sacct -j 32193027 -X -o State'`.
   If rc=0 → read `calibration_summary_32193027.tsv` (flip/tailflip/SUSPICIOUS per param) +
   `graded_clip_reshuffle_p2_32193027.tsv` (default tf0.6/tw12); compare to 1st-pass (3 flips). If rc≠0 →
   `tail slurm_p2_32193027.err`. (M1 watchers keep dying — sentinel is source of truth.)
3. ✅ DONE — 2nd pass + calibration (job 32193027, rc=0, guard active). Sweep result
   (`calibration_summary_32193027.tsv`, 1929 reads): tf=0.5 → 11-13 flips but **1 SUSPICIOUS** (false-pos
   creeps in); tf=0.6 → 3-4 flips, **0 SUSPICIOUS**; tf=0.7/0.8 → 1-2 flips, 0 SUS (over-conservative).
   term_window 8≈12, 16 = +1 legit flip. **CALIBRATION DECISION: tail_frac=0.6, term_window=12** (knee: max
   legit tail-flips at 0 SUSPICIOUS; 0.5 is unsafe, ≥0.7 misses flips). tw=16 = higher-recall alt (4 vs 3).
   Fix-deployed reshuffle = SAME 3 flips as 1st pass (902999d7/b164de40/de067271, all minimap2↔uLTRA
   tail_share=1.0) → walkback fix doesn't destabilize winner selection; cat2_plus_1's 24-A pattern is rare
   & absent from this 1000-read chunk. Blast radius stable: 3/1929, 0 suspicious.
4. ❌ CORE WIRING ABANDONED (2026-07-01) — **Option B is unsafe, do NOT wire.** The clean first audit was a
   STRAND-BUG artifact (grader only did the trailing clip; minus-strand tail is the LEADING clip → no-op on
   minus + spurious T-rich-5'-fragment flips = 2 of 3 first-pass flips). Strand-fixed the audit tool
   (commit: graded_clip_audit.py from_start on minus) → it then flips **cat1_minus_2 (the regression case):
   uLTRA @ CPA 15345 correct → minimap2 @ 15351 WRONG (6 bp)**. minimap2's leading clip CCTTT = AAAGG on the
   RNA = the GG genomic landmark uLTRA threads, but T-rich on the forward strand → a pure A/T-richness
   estimator can't tell mis-clipped genomic evidence from a real tail → re-breaks clip-to-win. Verdict +
   full arc in CASE_STUDIES cat2_plus_1. A future Option-B needs a stronger tail-vs-genomic discriminator
   (cross-aligner 3' consensus, or forbid forgiving a clip abutting a genomic non-stop landmark), validated
   on BOTH strands with cat1_minus_2 as a hard gate. cat2_plus_1 stays deSALT-winner. **Walkback fix fc44ee2
   stands independently — the ONLY shipped production change from this whole thread.** Scripts deployed to
`/scratch/users/kevinroy/graded_clip_audit/`. Sentinel `= /scratch/users/kevinroy/graded_clip_audit/.audit_rc`
(rc inside; 0=ok). Report `= /scratch/users/kevinroy/graded_clip_audit/graded_clip_reshuffle_32181800.tsv`.
RESUME/CHECK: `ssh sherlock 'cat /scratch/users/kevinroy/graded_clip_audit/.audit_rc 2>/dev/null; sacct -j
32181800 -X -o State'` — if rc=0 → fetch+read the reshuffle TSV (flip count, tail-flip vs SUSPICIOUS,
regression gate); if rc≠0 or FAILED → read `slurm_32181800.err` in that dir. First-pass caveat: deployed Oak
rectify is STALE (job log prints whether walkback has the fc44ee2 guard — expect False).

**AUDIT TOOLING BUILT + COMMITTED `2d0acfa`:** `scripts/benchmark/graded_clip_audit.py` (graded vs flat
winner reshuffle; terminal-window tail estimator; tail-flip/SUSPICIOUS classifier; local 36-read: 1 flip =
cat2_plus_1 tail-flip, 0 SUSPICIOUS, 00a1e01e no longer flips) + `scripts/benchmark/run_graded_clip_audit_
sherlock.sh` (larsms sbatch). **BLOCKED: SSH ControlMaster to `sherlock` expired → needs Kevin's interactive
2FA re-open** (do NOT force; BatchMode can't re-auth). Remaining deploy+submit (run AFTER Kevin re-opens ssh,
e.g. he types `! ssh sherlock true`):
```
scp scripts/benchmark/graded_clip_audit.py scripts/benchmark/run_graded_clip_audit_sherlock.sh \
    sherlock:/scratch/users/kevinroy/graded_clip_audit/
ssh sherlock 'cd /scratch/users/kevinroy/graded_clip_audit && sbatch run_graded_clip_audit_sherlock.sh'
# watch sentinel: /scratch/users/kevinroy/graded_clip_audit/.audit_rc (rc inside); report
#   /scratch/users/kevinroy/graded_clip_audit/graded_clip_reshuffle_<jobid>.tsv
# verify raw BAMs still on scratch first (Apr-22, ephemeral): ls $SRC/wt_by4742_rep1.*.bam
```

**Audit plan (sherlock-sbatch, `larsms` acct; chunk it):** (1) deploy `fc44ee2` to Oak rectify; (2) `rectify
correct` per aligner on a wt_by4742 read subset → per-aligner corrected BAMs+TSVs; (3) `merge_corrected_tsvs`
with current ED = baseline winners; (4) recompute ED with the GRADED terminal-window clip penalty → graded
winners; (5) report winner reshuffle: legit tail-clip flips vs 00a1e01e-class false positives; REGRESSION
GATE: clip-to-win (cat1_minus_2-class) must still lose. Then calibrate the terminal-window P→discount curve.
**Do NOT touch core `_cigar_hp_edit_distance` until this audit + gate pass.** Walkback fix `fc44ee2` is
shipped & independent. Prototype: `scratchpad/graded_proto.py`, `terminus.py`, `probe_00a1.py`.

**Resume (after user answers):**
1. If user says commit: re-record golden hash —
   `RECTIFY_RECORD_GOLDEN=1 conda run -n base python -m pytest tests/test_bam_parallel_state.py::test_process_bam_file_parallel_deterministic -s 2>&1 | grep "hash ="`
   then paste into the `GOLDEN_HASH_VALIDATION_MINIMAP2_NT2` constant in `tests/test_bam_parallel_state.py`.
2. Re-run that one test to confirm green, then `pytest -m "not slow"` sanity (expect all pass).
3. Add the `cat2_plus_1` case study to `rectify/data/validation/CASE_STUDIES.md` (homopolymer-undercall +
   the walkback-guard fix; note winner still deSALT pending EER-ED follow-up).
4. Surgical commit: `git add rectify/core/correct/walkback.py tests/test_validation_reads_upf1d.py
   tests/test_bam_parallel_state.py rectify/data/validation/CASE_STUDIES.md` + the regenerated
   `rectify/data/validation/rectified/**` + `rectify/data/validation/corrected_reads.tsv`
   (NEVER `git add -A`; do NOT touch COMPASS peer files: regen_2026_06_25/, compass_*, upf1d_2026_05/stage/).
   If user says DON'T commit fixtures: `git checkout -- rectify/data/validation/rectified rectify/data/validation/corrected_reads.tsv` to revert the regen.

**Files:** `rectify/core/correct/walkback.py` (fix), `tests/test_validation_reads_upf1d.py` (expectation),
`rectify/data/validation/CASE_STUDIES.md` (add cat2_plus_1), scratchpad `wb_probe5/6/7.py` + `BEFORE_corrected_reads.tsv`.

## ⚑ UPDATE 2026-06-27 (#4) — figure-review session (per-read renderer overhaul)

Interactive figure-by-figure review of the DRS validation bundle with the user. All changes are to
`scripts/validation_data/render_read_alignment.py` (+ `generate_review_report.py`), **render-only — no
test imports the renderer; full validation/poly-A/walkback suites remain green**. Committed in order:

- `a389793` unrectified-row fix: the "minimap2 (unrectified)" overview/per-base row prefers the dorado
  source, but the 7 build-X re-sourced reads (cat1/cat2/cat9) have a placeholder all-M CIGAR → painting it
  per-base frameshifts into an all-pink row. Fall through to `minimap2_unrectified.bam` when the dorado read
  is single-all-M.
- `eb54a81` (+ several `feat(validation-render)` commits squashed/amended along the way) — the per-read
  figure now shows, per aligner row:
  • **green ▲** = that aligner's RECTIFY-corrected 3' end (on EVERY row incl. the unrectified baseline);
    the unrectified row's ▲ marks the NAIVE (samtools, no-RECTIFY) 3' end for contrast.
  • **tail cells** decomposed non-overlapping: **crimson** = force-aligned-then-walkback-removed,
    **grey** = native aligner soft-clip, **green** = parquet pA tail.
  • ref row CLEANED — old `orig=corr`/`samtools` markers removed (redundant with per-row ▲).
  • the old vacuous bundle-bedgraph "3' pileup" top track was repurposed into an **"aligner 3' agreement"**
    track (each aligner's corrected 3' stacked as one labeled cell; consensus = one tall stack) but is now
    **opt-in, OFF by default** (`--agreement-track` / `render(show_agreement_track=True)`) — redundant with
    the per-row ▲s on the validation bundle; kept for real multi-read pileup scenarios.
  • `--dpi` option added (default 150).
  • **zoom indicator** (`6115776`): faint yellow box marks the detailed sub-region in the whole-read
    overview + two dotted fan lines down to the ref row, so the reader sees the per-base view is the 3' end.

**Walkthrough position: read 3/36 reviewed** (cat1_minus_1, _minus_2, _plus_1). NEXT = cat1_plus_2, then
cat2…cat9. Latest full render: `scratchpad/figs_v9/` (all 36). The renderer feature-set is settled per the
user's iterative feedback; resume the read-by-read review at cat1_plus_2.

## ⚑ UPDATE 2026-06-26 (#3) — WALKBACK extended-tail guard (force-aligned poly-A)

Figure review (cat1_minus_1) surfaced a real walkback gap. mapPacBio force-aligns the 3' poly-A tail
(no soft-clip, reference_start 9809 vs consensus 9831), and its corrected 3' came out **9811 vs the
consensus 9834** — it anchored on a coincidental read-A/genomic-A match at the minus-strand terminal.
Root cause: the tail-context false-stop guard's near window (`tail_context_k=4`) was CLEAN poly-T
(no mismatch) because the force-aligned-tail mismatch run only starts ~8 bp inward — so the guard
never fired, and the terminal non-stop match anchored the walk immediately. (The authors had already
added the leading-stop-match relaxation FOR this cat1_minus_1 pattern; the in-loop guard just couldn't
see far enough.)

**Fix** (`rectify/core/correct/walkback.py`, BOTH right- and left-side scans): added an extended
lookahead `tail_context_far=16`. When the near k-window is stop-base-only and mismatch-free, scan up to
16 aligned positions inward; if that wider run stays **stop-base-dominated (≥2/3)** AND is
**mismatch-rich (≥1/3)**, treat the candidate as still inside the force-aligned poly-A tail and keep
scanning. The ≥1/3 mismatch gate (the demonstrated case is ~56%) separates force-aligned tails from
genuine A/T-rich 3'UTRs with a lone basecall error (which would be ~6% and must NOT be over-walked).

**Verified (full pipeline, not the walkback-only proxy):** regen diff vs HEAD = exactly **1** change —
cat1_minus_1 mapPacBio 9811→**9834**, CONVERGING to the winner; **0 winner changes** in
corrected_reads.tsv. Tests: walkback units 54 passed; test_validation_reads 108 passed; full non-slow
suite (running). Figure: mapPacBio 3' now soft-clips the poly-A and lines up with the consensus.

## ⚑ UPDATE 2026-06-26 (#2) — POLY-A ROUND-TRIP: CORRECTED DIAGNOSIS (handoff #1's was WRONG)

The handoff's poly-A diagnosis ("parquet from a different build; aligned reads UNTRIMMED ==dorado len")
is **incorrect**. Verified the real state (the prior "412 vs 399" was read as a global parquet-vs-reads
skew; it is actually a localized dorado_source defect):

- **`val` reads ARE trimmed** (val ≠ dorado for all 36). **Parquet is coherent with `val`**, not stale:
  pq_trim == committed-val for **35/36** reads; pq_orig == build-X source for **36/36**.
- The earlier "val not a substring of dorado / independent basecalls" finding was an **artifact of the
  `=`-encoded SEQ** in `validation_reads.bam` (ref-relative; `query_sequence` returns `==A=C==…`). Decode
  vs reference before any seq comparison. (Minus-strand reads store the poly-A as leading poly-T — tail at
  the START — so `dor.startswith(val)` fails even for identical molecules.)
- **The authoritative build-X untrimmed source EXISTS**: `scripts/validation_data/rebuild_2026_05/
  combined_dorado_source.bam` (36 tagged reads, literal seqs) — its lengths == pq_orig for **36/36**, and
  committed-val == build-X-trimmed for 35/36. So parquet+val came from `combined_dorado_source.bam`.
  ⚠️ It is **NOT git-tracked** (scratch artifact) — copy its content into a tracked path to make durable.

### The bundle is build-X coherent EXCEPT two localized defects
1. **`validation_reads_dorado_source.bam`**: 6 reads were re-pulled from the WRONG basecall
   (`build_dorado_source.py` uses `wt_by4742_rep1.bam`), so their length ≠ pq_orig → **breaks test1**.
   Reads: cat2_plus_1(412/399), cat1_plus_2(449/451), cat1_minus_1(741/719), cat2_minus_1(656/655),
   cat1_minus_2(619/622), cat1_plus_1(187/167). FIX = source dorado_source from `combined_dorado_source.bam`.
2. **`validation_reads.bam` cat9_minus_2 (d3357db5)**: committed val=632 vs build-X-trimmed=640 (−8) →
   **breaks test2**. SEPARATE, val-side, 1 read; the dorado fix does NOT address it. Investigate before
   touching (changing val re-aligns → could shift a position).

### FIX PLAN (advisor-endorsed; does NOT touch val positions → cannot regress cat1_minus_1)
- Replace committed `validation_reads_dorado_source.bam` with `combined_dorado_source.bam` content
  (re-sorted/indexed, XV/XG preserved) → fixes test1 for all 6. Update `build_dorado_source.py` to use
  `combined_dorado_source.bam` as the authoritative source (it currently uses wt_by4742 → re-introduces skew).
- cat9_minus_2 test2: decide separately (re-derive that one val from build-X trimmed, or accept/xfail).
- The full re-trim+re-align ("trim-rebuild" that regressed cat1_minus_1) is NOT needed and is rejected.

### STATUS 2026-06-26 — POLY-A ROUND-TRIP: **FIXED** (all 3 tests pass; markers removed → live guards)
- **dorado_source swap DONE**: committed `validation_reads_dorado_source.bam` re-sourced from
  `combined_dorado_source.bam` (build X). 6 reads corrected, 30 byte-identical, 0 position changes. test1 passes.
- **cat9_minus_2 root cause = uLTRA read-trimming** (not build skew): uLTRA's `--reduce_read_ployA`
  (default 8) DROPS query bases from its emitted SEQ (640→632) despite help claiming otherwise. Per user's
  **no-trim policy** (aligners may soft/hard-clip — penalized, then rectify rescues — but may NOT delete
  bases), fixed in `rectify/core/align/multi_aligner.py::run_ultra`: inject `--reduce_read_ployA 10_000_000`
  (constant `_ULTRA_DISABLE_POLYA_REDUCE`) unless caller overrides. Verified: only extends terminal soft-clip;
  aligned core + reference_start + all other 35 reads byte-identical (blast radius = 1 read).
  - Regenerated committed `aligners/validation_reads.uLTRA.bam` cat9_minus_2 → 640 (via
    `update_validation_aligner_bams.py`, XV/XG preserved). Updated `validation_reads.bam` cat9_minus_2 → 640
    (copied fixed SEQ/CIGAR/quals onto bundle record, all consensus tags preserved). test2 passes.
- **Removed xfail markers** on both round-trip tests (assertions UNCHANGED; now live guards).
- **Reproducibility DONE** (user approved): `scripts/validation_data/combined_dorado_source.bam` now TRACKED
  (copied out of the gitignored rebuild_2026_05/, coord-sorted+indexed); `build_dorado_source.py` repointed to
  it as the authoritative source (old wt_by4742 path demoted to lazy fallback, only checked if a read is
  unresolved). Re-running the script is now idempotent + coherent (36 mapped, 0 missing; round-trip 3/3 pass).
- **VERIFIED**: full validation suite 109 passed / 0 failed; winner+position diff pre-vs-post uLTRA fix = 0 reads
  changed. cat9_minus_2 corrected_3prime + winner unchanged (restored bases are leading soft-clip, far from CPA).
- **Uncommitted change set** (mine — surgical add only; do NOT add the peer's dev/HANDOFF_SHORTREAD_P5.md or
  dev/HANDOFF_ALIGNER_BENCHMARK.md): multi_aligner.py, build_dorado_source.py, test_validation_polya_roundtrip.py,
  validation_reads.bam(.bai), validation_reads_dorado_source.bam(.bai), aligners/validation_reads.uLTRA.bam(.bai),
  scripts/validation_data/combined_dorado_source.bam(.bai), HANDOFF.md.
- **Full regen DONE** (regen_pa_rest_bundle.py): net committed-fixture change = 2 reads only —
  cat2_minus_2 128096→**128102** (now matches the test; fixes the real staleness) and cat7_plus_2 598986→598985
  (−1, no assertion). pA-rest restore acceptance gate PASSED: restored len == dorado len for all 36 (no
  double-count; the coherent parquet resolved the old "24 vs 12"). (Earlier 8-read estimate was inflated — that
  scratch merge omitted the penalty/overhang tables; the regen uses them and is authoritative.)
- **merge_corrected_tsvs idempotency BUG fixed** (`corrected_consensus.py`): the full non-slow suite surfaced
  5 failures (test_bam_writer_parallel_smoke ×4, test_run_single_sample_consensus ×1) — `KeyError
  'min_junction_anchor'`. Cause: the merge WRITES its internal `min_junction_anchor` column into
  corrected_reads.tsv (an omission — `hp_edit_distance`/`aligned_bases` were already dropped), so re-feeding
  that output collided (`_x`/`_y`) on the ed_df join. The stale HEAD fixture (pre-feature, 39 cols) masked it.
  Fix: (1) added `min_junction_anchor` to `drop_cols` (output back to 39-col schema); (2) defensive
  `all_df.drop([...3 computed cols...], errors='ignore')` before the ed_df merge → merge now idempotent.
  Both verified (10/10 in those two files pass). Re-running regen to rewrite corrected_reads.tsv as 39-col.
- **STILL OPEN**: (a) confirm full non-slow suite green after re-regen; (b) commit (bundled); (c) figure walkthrough.

## ⚑ UPDATE 2026-06-26 (#1) — NAMED WINNER BUG IS MOOT

Traced the winner-selection bug. **It does NOT reproduce on current HEAD.** Verified four ways:
1. **Fresh merge on committed BAMs:** cat2_minus_2→**128102**, cat4_minus_1→**128117** (both correct/test-expected).
2. **Committed deSALT BAM has NO 130168 mismap** — `9dbd37bf` aligns chrI:128112 (correct). The 130168
   only exists in the peer regen dir `scripts/validation_data/regen_2026_06_25/` (a different deSALT run).
3. **Even in that peer regen** the merge's own `comparison_summary.tsv` shows deSALT@130168 `_is_winner=False`;
   gapmm2@128102 (HP-ED 87.77) `_is_winner=True`. The handoff author read 130168 off deSALT's *per-aligner*
   `corrected_reads.region_000.tsv` and assumed it propagated — the merge had already rejected it.
4. **36-read invariant check** (fresh merge on committed bundle): winner == argmin(HP-ED) for **all 36 reads**;
   the chimera gate never overrode HP-ED. `pytest -k "cat2 or cat4"` = 19 passed.

→ **No merge code change needed.** `merge_corrected_tsvs` is sound. Do NOT touch the chimera sort.

### GENUINE remaining work (all user-gated — see OPEN DECISIONS below)
The only real discrepancy is the **stale committed fixture**. A fresh merge on the committed BAMs differs from
the committed `corrected_reads.tsv` for **8/36 reads** (all small shifts −83..+6 bp; no cross-chrom, no winner
swap). 3 move TO the test-expected value (cat2_minus_2 128096→128102, cat2_plus_2 8614→8605, cat4_minus_2
76251→76254); the other 5 (cat9_plus_2 −83, cat7_plus_2 −1, cat1_minus_2 +6, cat1_plus_2 −5, cat7_minus_1 +3)
have no test assertion. Tests pass because they run fresh; only the committed fixture + the figures built from
it are stale. **Decision pending: regenerate the committed fixture (+ figures) from a fresh merge?** This is the
"surgical" path — it does NOT touch BAMs or test assertions. cat1_minus_1 is NOT in this diff (stable on the
committed bundle; it only regressed under the *trim-rebuild*, which is a separate, rejected path).

---

## THE TASK (immediate — SUPERSEDED by the UPDATE above; kept for context)
Fix a winner-selection bug in `rectify/core/consensus/corrected_consensus.py::merge_corrected_tsvs`:
for some reads it selects a **spurious deSALT mismap with a far-WORSE HP-aware edit distance** over the
correct multi-aligner consensus. **User's governing principle: HP-aware ED is authoritative — lowest
HP-ED wins.** The bug is that the merge is NOT honoring that for these reads.

### Two exemplar reads (verified)
- **cat2_minus_2** (`9dbd37bf`): minimap2/gapmm2/uLTRA → chrI:128112, corr 3′=**128102**, HP-ED=**84.77**;
  deSALT → chrI:130168, HP-ED=**711.27**. Merge WRONGLY outputs 130168. Correct = **128102** (= test-expected).
- **cat4_minus_1** (`b956f764`): minimap2/gapmm2/mapPacBio/uLTRA → chrI:128117, corr 3′=**128117**,
  HP-ED=**11.42**; deSALT → chrI:130165, HP-ED=**334.26**. Merge WRONGLY outputs 130165. Correct = **128117**.

### Why deSALT's alignment is genuinely wrong (not a duplication)
deSALT's CIGAR at 130165 is BYTE-IDENTICAL to the consensus at 128117 (`66M8D2M1D47M1D52M2D45M1I87M2I11M
2I5M2D75M500N42M`, +2048bp), BUT it has **321/432 M-block mismatches** vs the consensus's **0/432**; the
two genome windows are only **26.1% identical**. A CIGAR encodes indel STRUCTURE, not base matches — deSALT
forced the same indel skeleton onto non-homologous sequence. So HP-ED 334/711 is DESERVED; lowest-ED →
consensus is correct.

### NEXT STEP
Trace `merge_corrected_tsvs` winner logic for `9dbd37bf` and find WHY it elevates the 711-ED deSALT over
the 85-ED consensus. Candidates: deSALT secondary record @128k scored while primary's 130k POSITION
reported; the chimera/`_n_agree`/drop-chimeric path overriding the ED sort; ED computed on wrong record.
Verify with: `_cigar_hp_edit_distance` (in that file) + `HpPenaltyTable.from_tsv` (penalty tables at
`rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/`). Fix so lowest-HP-ED wins → cat2_minus_2
→128102, cat4_minus_1→128117. Then `pytest tests/test_validation_reads.py -k cat2 -q` (slow ~2min).

## VERIFIED context
- **Root problem:** the validation FIGURES were built from STALE committed fixtures
  (`rectify/data/validation/rectified/per_aligner_summary.tsv` + `corrected_reads.tsv`, frozen 2026-05-18).
  Renderer made fresh-correct (commit `386814e`, `--fresh` opt-in + staleness banner). TESTS were always
  fresh (they run `rectify correct`+merge in tmp), so the regression suite is sound.
- **Poly-A round-trip is broken:** aligned reads (`validation_reads.bam`) are UNTRIMMED (==dorado_source
  len); trim parquet is from a DIFFERENT build (orig 412 vs actual 399); restore stacks parquet tail on
  the walkback clip (24 vs 12). Guard test committed (`57f3586`,
  `tests/test_validation_polya_roundtrip.py`, 2 xfail + 1 pass) + `regen_pa_rest_bundle.py --emit-merged-tsv` fix.
- **Coherent rebuild PROTOTYPED** (trim dorado_source → matching parquet → re-align trimmed reads):
  artifacts in `/tmp/cat2_rebuild/` (trimmed.bam/fastq.gz, trim_metadata.tsv [matching!], tagged BAMs,
  corrected/, merged_corrected.tsv) + Sherlock `/scratch/users/kevinroy/cat2_rebuild_desalt/`. Diff vs
  committed: 11/36 corrected_3prime changed.

## OPEN DECISIONS (user-facing; do NOT blind-decide)
1. **Surgical fix vs full trim-rebuild.** User leaning unclear. The trim-rebuild REGRESSES cat1_minus_1
   (9834→9852; test expects 9834) and leaves cat1_plus_2 wrong — so it's NOT a clean win. The deSALT
   winner bug is INDEPENDENT of trim (the mismap wins on the committed bundle too) → can be fixed
   surgically. **Recommended: (A) surgical = winner-selection fix + parquet round-trip fix on the
   existing bundle, NOT the full trim-rebuild.** Confirm with user.
2. **7 small shifts** (±1–18bp, legit trimmed-read diffs): cat2_plus_2(→8605✓test), cat4_minus_2(→76254
   ✓test) MATCH test (committed was wrong); cat1_minus_1 REGRESSES; others have no test assertion. Do NOT
   update `test_validation_reads.py` assertions without per-read user sign-off.
3. **cat2_minus_2 committed corrected_reads.tsv=128096 but test expects 128102** — even the old bundle was
   self-inconsistent (build skew). Target = 128102.

## RESUME
1. `cat HANDOFF.md` (this file). `ls /tmp/cat2_rebuild/` — if cleared (reboot), re-derive per "Coherent
   rebuild" (trim_drs_bam_polya on dorado_source → fastq → `rectify align --Scer --aligners minimap2
   mapPacBio gapmm2 --junction-aligners uLTRA --no-consensus`; deSALT via Sherlock sbatch pattern
   `/tmp/desalt_rebuild.sbatch` job 31409555; not strictly needed for the winner-bug trace).
2. Trace+fix the winner-selection bug (see NEXT STEP). This is the priority and is independent of /tmp.
3. Re-run `pytest tests/test_validation_reads.py -k "cat2 or cat4" -q`; expect 128102 / 128117.
4. Present any test-assertion changes to the user before committing.

## FILES / STATE
- Committed this thread: `57f3586` (round-trip guard + regen fix), `386814e` (--fresh opt-in). HEAD may be
  further (concurrent peer commits). Working tree clean except this HANDOFF.md.
- Env: M1 base conda `gapmm2==25.4.5` (do not upgrade — NEW-082). deSALT cannot run on M1 (x86); use Sherlock.
- `git add` EXPLICIT paths only. Don't touch `dev/COMPASS_*`, `scripts/validation_data/regen_2026_06_25/`
  (peer/other-thread).
