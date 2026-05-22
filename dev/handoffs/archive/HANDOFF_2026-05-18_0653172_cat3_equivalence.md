# Session Handoff — Cat3 5'-rescue equivalence-extension (parallel to UMI-bin work)

**Date:** 2026-05-18 (afternoon → late evening)
**Branch:** `drs-validation-rebuild` (off `master`)
**Top of branch:** `0653172` (fix: equivalence-extension slice + + strand mirror + k-sweep)
**Parallel handoff:** the in-tree `HANDOFF.md` describes a *separate* UMI-bin
Phase D calibration work track. **This handoff is independent of that** —
different code paths, different files, different test surface. The next
agent for cat3 work should read THIS file; the next agent for the calibration
work should read the in-tree `HANDOFF.md`.

**Scope of this session:** shipped a 5'-rescue equivalence-extension that
eliminates the `2D-N` cosmetic CIGAR for cat3-style rescued reads on - strand,
plus a follow-up that fixed a slice-geometry bug, added + strand mirror
plumbing + trigger, and added a k-sweep. Five commits total.

---

## 1. What was done

- **Doc cleanup pass + cat2/cat3 re-classification** (`540f73c`, `30dbcf9`).
  Re-investigated `cat2_minus_2` "Module 2G lookahead" — found the queue
  description was geometrically imprecise (intervening bases A+T not pure
  HP-T; 4 not 3 matches; ambiguous endpoint). Re-investigated cat3_minus_2 —
  found the TSV-level output was already canonical (`junctions=366502-366584`
  for all 5 aligners); the visible issue was BAM-CIGAR cosmetic only. Marked
  bedgraph-regen and per-aligner-effective-utility items as RESOLVED by
  prior commit `75b0338`.

- **Reviewed `validation_read_review/cat3_junction_findings.md`** — found the
  plotter's primary recommendation was actually the *chimera-exemption fix*
  in `merge_corrected_tsvs`, not the 5'-rescue extension. Blast-radius
  diagnostic showed the chimera fix would flip only 1 of 36 validation
  reads (cat3_minus_2: gapmm2 vs deSALT, hp_ed 20.63 vs 21.49). User
  pivoted to option B (5'-rescue extension) for broader real-world benefit.

- **Implemented 5'-rescue equivalence-extension** (`6d2cf59`). When body
  alignment overshoots annotated intron boundary by *k* bases AND borrowed
  query bases match both the old and new positions, prepend (− strand) /
  append (+ strand) the borrowed bases to `_align_seq` so
  `align_clip_to_exon` produces a clean M op covering both regions. New
  TSV column `five_prime_upstream_trim`. Threaded through 5 files. Bundle
  regenerated.

- **Refreshed handoff + queue docs** (`accc7df`). Recorded the cat3 cleanup
  outcome + added a design note at the top of `debugger_queue.md` for the
  splice-junction ambiguity-window + motif-strength tiebreaker (queued as
  next-session work; structurally a consensus-selection concern, not a
  rescue concern).

- **Slice bug fix + + strand mirror + k-sweep** (`0653172`). Three things:
  - The initial commit's `_ref_left` slice was off by k positions
    (`genome[intron_start - k : intron_start]` instead of the geometrically
    correct `genome[ref_end - k : ref_end]`). For cat3_minus_2 these
    slices both happened to be `CT` because chrII at 366500-366503 has a
    CT-CT repeat — the test passed by coincidence, not by correctness.
    Fixed.
  - Added + strand overshoot trigger with the symmetric criterion
    `genome[ref_start : ref_start + k] == genome[intron_start - k : intron_start]`.
    Plumbing through `extend_read_5prime_for_junction_rescue` already
    handled + strand symmetrically; only the rescue-side detection was
    missing.
  - Added k-sweep over `k = min(overshoot, 10)` downward, accepting the
    largest k where equivalence holds. Enables partial absorption when
    full-k doesn't qualify.

---

## 2. What's verified

- `pytest tests/test_validation_reads.py --no-header -q` → **106 passed, 8
  skipped** (was 105/8 baseline; +1 from the new
  `test_cat3_minus_2_rescued_aligners_have_clean_intron_cigar`). Run twice:
  after `6d2cf59` and after `0653172`.
- cat3_minus_2 per-aligner trimmed BAM CIGAR (post-regen, verified via
  direct `pysam` inspection):
  - minimap2: `… D1 M256 N(366502-366584) M15` (canonical)
  - gapmm2:   `… D1 M256 N(366502-366584) M15` (canonical, unchanged)
  - mapPacBio: `… D1 =256 N(366502-366584) =15` (canonical, unchanged)
  - deSALT:   `… D1 M256 N(366502-366584) M15` (canonical, was `M258 N(80) D2 M13`)
  - uLTRA:    `… D1 =256 N(366502-366584) M15` (canonical, was `=258 N(80) D2 M13`)
- corrected_reads.tsv row for cat3_minus_2 records
  `five_prime_upstream_trim=2` + `five_prime_exon_cigar=15M`.
- All 5 aligners now share `hp_edit_distance = 20.6323` for cat3_minus_2
  (was 20.63 / 20.83 / 21.49 split before the fix).

**NOT VERIFIED:**

- **+ strand mirror code path has no test that exercises it.** The criterion
  `genome[ref_start : ref_start + k] == genome[intron_start - k : intron_start]`
  was derived from scratch by mirror reasoning from the - strand case. Type
  checks pass and pytest doesn't regress, but no synthetic or real + strand
  cat3 read in the bundle has the overshoot pattern needed to trigger it.
  cat3_plus_2 is an *undershoot* pattern (`ref_start < intron_end`), so it
  doesn't exercise this mirror either.
- **k-sweep at k > 2 has no test exercise.** The only test case
  (cat3_minus_2) has overshoot k=2 and the equivalence holds at k=2 on the
  first try. Edge cases not tested: overshoot=0 (sweep should no-op),
  overshoot > _MAX_K (should clamp), k=2 fails but k=1 succeeds (partial
  absorption).
- **Interaction with the existing shift loop in `rescue_3ss_truncation`
  lines 1191-1218.** That loop already searches over `_shift` values to
  find the best junction placement. My equivalence-extension runs AFTER
  it (after `best_junction` is committed). Not deeply traced; no
  regression observed but worth a careful read.

---

## 3. Open items

- **cat3_plus_2 "off-by-1 acceptor" undershoot pattern.** 4 of 5 aligners
  on cat3_plus_2 place the intron 1 bp short of canonical: `N(142253,142618)`
  instead of `N(142253,142619)`. This is an UNDERSHOOT on the acceptor side
  (`ref_start < intron_end`), not an overshoot. My fix only handles
  overshoot. **Why deferred:** the undershoot is a structurally different
  transformation (extend body M by k bases, taking from the END of the
  upstream softclip-aligned exon ops, rather than absorbing the trailing D)
  — needs its own equivalence criterion and code path. Detection logic for
  + strand undershoot would go in the same `rescue_3ss_truncation` block
  alongside the overshoot triggers.

- **Splice-junction ambiguity-window + motif-strength tiebreaker.** Full
  4-commit design sketch at the top of `docs/handoffs/debugger_queue.md`.
  This is the consensus-selection concern that the user originally
  motivated with "We check upstream of the called 5' SS and 3' SS for the
  same nt(s)..." — it goes in `merge_corrected_tsvs`, NOT in the rescue
  function. **Why deferred:** structurally separate scope (winner selection,
  not rescue); needs its own session with fresh context on the chimera
  exemption logic at `corrected_consensus.py:760-763`. The advisor flagged
  earlier this session that the blast radius of any tiebreaker change is
  narrow (only cat3_minus_2 in the validation set would flip; check
  production scope before assuming impact).

- **Cat1 HP-mode metric design** — architectural, still deferred from
  prior sessions. See `docs/handoffs/debugger_queue.md` Cat1 entries.

- **Phase C 5' rescue calibration application.** The report at
  `docs/handoffs/anchored_prefix_calibration.md` recommends a threshold of
  "≥4 consecutive matches in first 6 bp." But the advisor flagged early
  this session that the calibration profiled random soft-clip-vs-upstream-
  genome (the wrong null) — translating that recommendation to the existing
  `_rescue_5prime_softclip`'s junction-specific threshold needs a
  precision/recall measurement on actual rescue cases first. Don't ship
  the threshold change without that intermediate measurement.

- **`[uncommitted]` raw aligner BAM WIP** —
  `rectify/data/validation/aligners/validation_reads.{minimap2,uLTRA}.bam`
  have content changes from before this session (mtime 2026-05-18 10:07).
  Not introduced by this session, intentionally not committed. Next session
  should ignore unless the user surfaces them.

---

## 4. Resume command

**Resume:** double-check the cat3 equivalence-extension by running:

```bash
cd /Users/kevinroy/work/rectify
git log --oneline 0653172^..0653172  # the slice-fix + + strand + k-sweep
git show 0653172 -- rectify/core/splice/splice_aware_5prime.py | less
```

Then:

1. **Manually trace the + strand mirror geometry** against a hypothetical
   cat3_plus_X read where the body M starts at `intron_end + k` (overshoot
   on acceptor side). Verify the equivalence criterion
   `genome[ref_start : ref_start + k] == genome[intron_start - k : intron_start]`
   is right. The corresponding read-side check is that the FIRST k query
   bases of body M (just after the 5' softclip) equal both genome slices.

2. **Construct a synthetic test** in `tests/test_validation_reads.py`
   that exercises the + strand overshoot path. The validation bundle has
   no real + strand cat3 read with this pattern (cat3_plus_2 is undershoot;
   cat3_plus_1 has ref_start = intron_end exactly = no overshoot). A
   synthetic-BAM unit test in a fresh tests/test_splice_equivalence.py
   would be cleaner than adding to the validation suite.

3. **Trace the k-sweep edge cases.** With my code:
   - `overshoot=0` → `if 0 < _overshoot` fails → no-op ✓ (verified by
     reading code; not test-exercised)
   - `overshoot > 10` → `_max_k = min(overshoot, 10) = 10` → sweep
     1..10 ✓ (verified by reading; not test-exercised)
   - `overshoot=2, k=2 fails, k=1 succeeds` → sweep should land at k=1
     for partial absorption. NOT test-exercised; verify by reading the
     loop body.

4. **Verify the `_ref_old == _borrowed AND _ref_new == _borrowed` joint
   check.** Two separate conditions:
   - `_borrowed == _ref_old`: the body M was a REAL match at those
     positions (filters out cases where body had op-code M but with
     mismatched bases — i.e., minimap2's encoded mismatch).
   - `_borrowed == _ref_new`: the slide-to-new-position would also be a
     valid match.
   - Both required → the equivalence-extension is safe to apply.

5. **Check interaction with the existing shift loop** in
   `rescue_3ss_truncation` lines 1191-1218. That loop searches for the
   best junction placement within `(_shift_lo, _shift_hi)`. My
   equivalence-extension runs AFTER `best_junction` is committed
   (lines 1311 onward). Verify there's no case where the shift loop's
   chosen `_intron_end` is non-annotated (then my equivalence check
   against the chosen `_intron_end` might fire on a wrong placement).

**If all 5 checks pass:** the equivalence-extension is in good shape.
Optional next step: cat3_plus_2 undershoot pattern (item 1 in §3).

**If the user instead wants the ambiguity-window + motif tiebreaker:**
read the design note at the top of `docs/handoffs/debugger_queue.md` and
execute the 4-commit sketch there. That work is independent of cat3
equivalence-extension review.

---

## 5. Files touched

Across this session's commits (`540f73c`, `30dbcf9`, `accc7df`, `6d2cf59`,
`0653172`):

**CODE:**

- `rectify/core/splice/splice_aware_5prime.py` — equivalence-extension
  logic in `rescue_3ss_truncation` (commits `6d2cf59`, `0653172`). The
  block to read carefully is approximately lines 1298-1410 of that
  function (after `best_junction` is set, before the `return` dict is
  built).
- `rectify/core/bam/bam_processor.py` — thread `five_prime_upstream_trim`
  from `_3ss_result` into the per-read TSV dict.
- `rectify/core/bam/output.py` — new TSV column header + write.
- `rectify/core/bam/bam_writer.py` — read TSV column, pass to extend
  function (both streaming and parallel code paths).
- `rectify/core/bam/read_edits.py` — `extend_read_5prime_for_junction_rescue`
  accepts `upstream_trim`, applies the body-M trim symmetrically for both
  strands. Both branches updated.

**TEST:**

- `tests/test_validation_reads.py`:
  - `_run_correct_first_pipeline` now returns
    `(merged_tsv_path, per_aligner_corrected_bams)`.
  - `_correction_outputs` fixture returns BAM paths too.
  - New fixture `per_aligner_corrected_bams_no_2h`.
  - New test method
    `test_cat3_minus_2_rescued_aligners_have_clean_intron_cigar` in
    `TestCategory3JunctionRescue`.

**DOCS:**

- `docs/handoffs/debugger_queue.md` — cat2_minus_2 re-classification, cat3
  RESOLVED status with geometric details, ambiguity-window design note at
  the top, original-queue section preserved below.
- This file (`handoffs/HANDOFF_2026-05-18_0653172_cat3_equivalence.md`).

**DATA (bundle regen output from `6d2cf59`):**

- `rectify/data/validation/corrected_reads.tsv` — new
  `five_prime_upstream_trim` column.
- `rectify/data/validation/rectified/per_aligner_summary.tsv` — re-merged.
- `rectify/data/validation/rectified/corrected_3ends.{plus,minus}.bedgraph`
  — refreshed counts.
- `rectify/data/validation/rectified/corrected_reads.{plus,minus}.bedgraph`
  — new files (bedgraph writer output, parallel to `corrected_3ends.*`).
- `rectify/data/validation/rectified/per_aligner/{5 aligners}.{trimmed,
  softclipped}.bam(.bai)` — new files.
- `rectify/data/validation/rectified/rectified_corrected_3end.bam(.bai)`
  — refreshed (cat3 canonical CIGAR visible here).
- `rectify/data/validation/rectified/rectified_pA_tail_{trimmed,
  soft_clipped}.bam(.bai)` — refreshed.
- `rectify/data/validation/validation_reads.bam(.bai)` — consensus output
  refreshed.

**`[uncommitted]` — working tree items present BEFORE this session and not
changed here:**

- `HANDOFF.md` — replaced by user with UMI-bin Phase D content. A
  *separate* work track. Don't conflate.
- `rectify/data/validation/aligners/validation_reads.{minimap2,uLTRA}.bam`
  — pre-existing WIP modifications to raw aligner BAMs from before this
  session (mtime 2026-05-18 10:07). Not from my session; not committed.
- `CLAUDE.md`, `README.md`, `docs/figures/*` — pre-existing WIP, unrelated
  to this session.
