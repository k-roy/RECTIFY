# Session Handoff — cat3_plus_2 + cat2_minus_2 shipped; mpb reanchor deferred

**Date:** 2026-05-18 (evening, continuing from `01c2a18_plus_strand_disabled`)
**Branch:** `drs-validation-rebuild`
**Top of branch:** `6943450` (cat2_minus_2 2-bp del + reanchor sketch)
**Predecessor handoffs:**
- `handoffs/HANDOFF_2026-05-18_0653172_cat3_equivalence.md` — original cat3_minus_2 - strand work
- `handoffs/HANDOFF_2026-05-18_01c2a18_plus_strand_disabled.md` — verification that caught the + strand geometry inversion

---

## 1. What was done

Three substantive commits, two shipped fixes + one shipped-but-unhooked sketch.

### `acb508e` — fix(splice): cat3_plus_2 + strand mirror

Two coupled fixes that together produce canonical 366 bp intron for any
+ strand 5'-rescue read whose body M starts INSIDE the canonical intron
(the "off-by-1 acceptor" pattern affecting 4 of 5 aligners on cat3_plus_2):

- **`rescue_3ss_truncation` + strand equivalence-extension re-enabled
  with corrected trigger.** Replaces the conservative no-op from `01c2a18`.
  Correct trigger is `_intron_end - read.reference_start > 0` (undershoot:
  body M starts INSIDE the intron from the high-coord edge), NOT
  `read.reference_start - _intron_end > 0` (the inverted sign that
  produced `canonical + 2k` intron). Surgery routine in
  `extend_read_5prime_for_junction_rescue` was already correct for this
  case: `intron_len` collapses to `intron_end - intron_start`.

- **`reroute_intronic_tail_5prime_via_junction` + strand off-by-one.**
  The + strand branch was treating `five_prime_position` as `intron_start`
  when it is actually `intron_start - 1` (the LAST genomic position of
  the upstream exon). Fix: `intron_len = clip_boundary - five_prime_position - 1`
  and `new_ref_start = five_prime_position - exon_ref_span + 1`. The
  - strand branch was correct because - strand `five_prime_position` IS
  `intron_end` directly.

  This bug was previously masked because the equivalence-extension didn't
  fire on + strand (was disabled). Once `acb508e`'s first half re-enabled
  it, the buggy reroute partially fired AFTER my equivalence-extension
  and re-introduced a +1 to `intron_len` and a -1 to `new_ref_start`,
  producing `N(367)` instead of canonical `N(366)`. Diagnostic trail
  identified this via per-step file logging in `/tmp/cat3p2_dbg.log`.

  Added regression test
  `test_cat3_plus_2_rescued_aligners_have_clean_intron_cigar` in
  `tests/test_validation_reads.py:681`. Mirrors the existing cat3_minus_2
  cigar test.

### `6943450` — fix(rescue): cat2_minus_2 2-bp del extension + reanchor sketch

**cat2_minus_2:** Post-pass in `rescue_softclip_at_homopolymer`'s
left-side branch (minus-strand 3'-end rescue): when the main walk
terminates on a softclip-vs-ref mismatch (not poly-T entry /
chromosome boundary), try absorbing 2 more ref positions into the
deletion and walking outward until first true mismatch. Accepted only
when ≥3 additional bases match (the TGC threshold per user directive
2026-05-18: "Allowing for 3 straight nt to match should outweigh the
single T deletion"). The 2-bp non-HP skip is rolled into
`homopolymer_extension` so the downstream BAM surgery in
`extend_read_3prime_for_softclip_rescue` produces `D(HP+2) + M(rescued)`
without needing a new plumbing path.

For cat2_minus_2 (chrI, raw 128113, minimap2 winner): corrected_3p
moves from **128102** (1-base A@128102 match only) to **128096** (TTGC
motif at 128096-128099 + 2-bp skip past A@128100 / T@128101). Test
updated: `test_3prime_exact_position[cat2_minus_2-128096]` at
`tests/test_validation_reads.py:497`.

**Reanchor sketch (DEFERRED):** `reanchor_5prime_for_rescue` in
`rectify/core/bam/read_edits.py:38` is functionally complete but NOT
wired. Strand-aware logic walks the CIGAR from the 5' edge inward to
find the first sustained match run of length ≥10 (M ops checked
per-base against genome; =/X read directly; S/I/D/N break the run)
and collapses everything upstream into a leading soft-clip. See §3 for
why the hooks were reverted.

### `5c13637` (predecessor session) — handoff doc for the disable-and-document `01c2a18` fix

For context only; not introduced this session.

---

## 2. What's verified

- `pytest tests/test_validation_reads.py --no-header -q` →
  **107 passed, 8 skipped** (was 106/8 at start of session; +1 from the
  new cat3_plus_2 cigar test).
- cat3_plus_2 per-aligner CIGAR (post-`acb508e`): all 4 off-by-1 aligners
  (minimap2, deSALT, uLTRA, gapmm2) now share `N(142253, 142619)`
  intron of length 366 (canonical). mapPacBio unchanged.
- cat3_minus_2 unchanged from `0653172` baseline (- strand equivalence-
  extension still fires; canonical M256 N(366502-366584) M15).
- cat2_minus_2 corrected_3p = 128096 (was 128102).
- Bundle BAMs **NOT regenerated this session** — the user's image
  showing the old "winner cluster has 2D, mapPacBio has 1D" geometry
  reflects the pre-`acb508e` bundle. Bundle regen needed to surface
  the new HP-ED-driven consensus pick (expected: mapPacBio wins on
  body-quality once all 5 aligners have canonical intron).

**NOT VERIFIED:**

- **mpb 5'-rescue for cat3_minus_1 / cat3_plus_1.** Reanchor function
  works in isolation (verified after the `run_end_off = k` fix —
  initial code had `run_end_off = length - 1 - k` which dropped
  post-run ops entirely). But hooks were reverted (see §3).
- **Cat3_plus_2 consensus pick.** Per the user's image, the winner
  cluster is currently `deSALT, gapmm2, minimap2, uLTRA` (not
  mapPacBio). User predicts HP-ED will naturally favor mapPacBio after
  bundle regen because mapPacBio's exon-1 has a single 1D vs the
  winner cluster's 2D. NOT confirmed empirically this session.
- **Pre-existing tiebreaker regressions in
  `test_corrected_consensus_tiebreaker.py`** —
  `test_paralog_tiebreaker_picks_multi_aligner_consensus` and
  `test_majority_consensus_picks_chrXIV_even_when_outlier_has_wider_span`
  fail with `'float' object has no attribute 'split'` in `_eff_key`
  (`corrected_consensus.py:830`). Confirmed pre-existing on master
  (reproduced by stashing this session's changes). Not introduced by
  any of `acb508e` / `6943450`.

---

## 3. Open items

### Reanchor unhook (highest-priority deferred work)

**Hook attempt summary:** Added `reanchor_5prime_for_rescue` calls in:
- `bam_processor.py` before `_rescue_3ss` (line ~407)
- `bam_writer.py` inside the `if correction['five_prime_rescued']`
  conditional, before `extend_read_5prime_for_junction_rescue` (line ~270)

**Result:** 6 regressions: cat3_minus_2_rescued_aligners_have_clean_intron_cigar,
cat3_plus_2_rescued_aligners_have_clean_intron_cigar, cat4_plus_2-393721,
cat7_plus_1, cat7_plus_2. cat3_plus_2 actually passes in isolation —
cross-test interference suggests module-scope state pollution between
fixture runs.

**Diagnosed root cause:** Cross-stage CIGAR inconsistency. Workflow:
1. `bam_processor.run_correct` reads the raw input BAM and calls
   `rescue_3ss_truncation` to compute the TSV row (including
   `five_prime_exon_cigar` and `five_prime_upstream_trim`). With the
   bam_processor reanchor hook, the read's in-memory CIGAR is mutated
   before `_rescue_3ss` runs, so the TSV is sized to the reanchored
   geometry.
2. `bam_writer.write_corrected_bam` reads the raw input BAM AGAIN
   (fresh pysam iterator), receives the ORIGINAL un-reanchored CIGAR.
   With the bam_writer reanchor hook, the CIGAR is reanchored again
   before surgery. SHOULD produce the same geometry as bam_processor
   saw.
3. Surgery uses `actual_sc = cigar.pop(0)[1]` from the live CIGAR
   (post-reanchor) and `exon_cigar_str` from the TSV. If the two
   reanchor results diverge by even 1 base, `exon_query_span` vs
   `expected_query` mismatches and the surgery falls back to flat M
   (line `if exon_query_span != expected_query: exon_ops = [(0, expected_query)]`),
   corrupting the read.

The reanchor function is deterministic given (read, genome,
anchor_min_run), so calls (1) and (2) SHOULD produce the same result.
But there are intervening steps in `write_corrected_bam` (notably
`realign_exon_blocks` at line 266) that mutate the CIGAR between the
fresh-from-disk read and the reanchor call. realign_exon_blocks does
NOT run in bam_processor. So the bam_writer's reanchor sees a different
CIGAR than the bam_processor's reanchor saw → divergence → corruption.

**Unhook path (a):** Reanchor in bam_processor BEFORE
`rescue_3ss_truncation`, AND in bam_writer BEFORE
`realign_exon_blocks`. Risk: realign might still alter the post-reanchor
CIGAR in ways that don't match what bam_processor saw (bam_processor
doesn't run realign at all). Need to confirm realign is a no-op on
already-reanchored reads, or move realign to bam_processor too.

**Unhook path (b) — RECOMMENDED:** Persist reanchor as a TSV column.
- Add `reanchor_clip_len: int` to the result dict from
  `rescue_3ss_truncation` (or compute in bam_processor adjacent to the
  rescue call).
- Add column to `rectify/core/bam/output.py:50` schema.
- Read it in `bam_writer.py:133` and pass as a new parameter to
  `extend_read_5prime_for_junction_rescue`.
- Function signature gains `reanchor_clip_len: int = 0`. When > 0,
  the function knows to look for the leading S of THAT length (not
  just the leading S op as-is), accounting for the synthetic S that
  reanchor produces.

This avoids any in-memory CIGAR mutation in the correction pass. The
raw BAM stays raw; only the surgery output reflects the reanchor.

Pseudocode reference for reanchor at
`validation_read_review/cat3_junction_findings.md:252-282`. Function
implementation at `rectify/core/bam/read_edits.py:38`.

### Bundle regeneration

The validation bundle in `rectify/data/validation/rectified/` was not
regenerated this session. After `acb508e`, cat3_plus_2 per-aligner BAMs
will show canonical 366 bp intron for all 5 aligners; the consensus
summary should re-rank mapPacBio first on body-quality HP-ED.

Run: `scripts/validation_data/regen_pa_rest_bundle.py` (or equivalent
per `methods/rectify_quickstart_h2.md`). Verify
`rectify/data/validation/rectified/per_aligner_summary.tsv` row for
cat3_plus_2 (`79f61403-cf63-4522-b555-569590dc4304`) shows
`_winning_aligner=mapPacBio`.

### Pre-existing tiebreaker test failures (separate issue)

`test_paralog_tiebreaker_picks_multi_aligner_consensus` and
`test_majority_consensus_picks_chrXIV_even_when_outlier_has_wider_span`
in `tests/test_corrected_consensus_tiebreaker.py` fail with
`'float' object has no attribute 'split'` at
`rectify/core/consensus/corrected_consensus.py:830`. The `junctions`
column has a NaN (float) for at least one row, and `_eff_key` calls
`.split(';')` on it unconditionally.

Likely fix: `_juncs = _row.get('junctions', '') or ''` is already
defensive but pandas `apply` may pass NaN-as-float through `.get` if
the column literally holds `np.nan`. Try
`_juncs = _row.get('junctions') or ''` after coercing NaN→'' upstream,
or `_juncs = '' if pd.isna(_row.get('junctions')) else _row.get('junctions')`.

### Open items from prior handoffs (still open)

- **Cat3_plus_1 chimera exemption** (per `cat3_junction_findings.md:90`)
  — separate from the reanchor work; affects whether mapPacBio wins
  cat3_minus_2 consensus. Per user 2026-05-18: don't need a hack; once
  CIGARs converge after canonical-intron fixes, HP-ED naturally
  decides. Verify after bundle regen.
- **Cat2_minus_2 and other cat2 entries** in `debugger_queue.md` —
  cat2_minus_2 is now SHIPPED. Other cat2 design notes remain open
  (cat2_plus_1 HP-ED del_cost calibration at line 313 of queue;
  cat2_plus_2 effective-utility feature request RESOLVED at line 448).
- **Ambiguity-window + motif-strength tiebreaker** for consensus
  selection (queued at top of `debugger_queue.md`).
- **Cat1 HP-mode metric design** (architectural, deferred from prior
  sessions).
- **Phase C 5' rescue calibration application** — flagged in `0653172`
  handoff §3; requires precision/recall measurement before shipping
  the threshold change.
- **Cat9 Module 2H** — user-deferred per queue line 657; no action
  this session.

---

## 4. Resume command

For unhook path (b) — the recommended next-session work:

```bash
cd /Users/kevinroy/work/rectify
git log --oneline 6943450^..6943450  # see what was deferred
sed -n '38,200p' rectify/core/bam/read_edits.py  # reanchor function + DEFERRED docstring
head -75 docs/handoffs/debugger_queue.md  # full unhook sketches
```

Then:

1. **Add `reanchor_clip_len` to `rescue_3ss_truncation` result dict.**
   `rectify/core/splice/splice_aware_5prime.py:1416-1425` (the `return
   {'rescued': True, ...}` block). Compute by calling
   `reanchor_5prime_for_rescue` against a COPY of the read (or compute
   the would-be S length without mutating), and emit the new field.

2. **Plumb through bam_processor.py.** `rectify/core/bam/bam_processor.py:412`
   currently reads `five_prime_upstream_trim` from `_3ss_result`. Add a
   parallel line for `reanchor_clip_len`. Then thread into the result
   dict at lines 522-528 (which becomes a TSV row column).

3. **Add TSV column to output.py.**
   `rectify/core/bam/output.py:50` (header schema) and `:95` (write).
   Place it alongside `five_prime_upstream_trim`.

4. **Read in bam_writer.py.** `rectify/core/bam/bam_writer.py:133`
   (header parsing) and `:180` (correction dict population). Pass to
   `extend_read_5prime_for_junction_rescue` as a new kwarg
   `reanchor_clip_len=...`.

5. **Update `extend_read_5prime_for_junction_rescue`** to accept the
   new kwarg. When > 0 and the read's leading S is shorter than
   `reanchor_clip_len`, call `reanchor_5prime_for_rescue` on the read
   FIRST (in bam_writer's processing) so the live CIGAR has the
   correct S, then proceed with surgery. OR: trust that the leading S
   in the BAM is what reanchor produced and use it as-is (requires
   bam_writer to have already reanchored).

6. **Verify cat3_minus_1 and cat3_plus_1.** After bundle regen, both
   reads' mapPacBio should produce clean `S(reanchor_clip_len) {body_M
   ... clean N ...}` CIGAR matching the winner cluster's geometry.
   No new test needed if the existing `test_5prime_present` and
   `test_5prime_exact_position` cover them (they do —
   `tests/test_validation_reads.py:516-542`).

7. **Run validation suite.** `pytest tests/test_validation_reads.py
   --no-header -q`. Expect 107 passed (no regressions) and possibly +1
   if you add a new test asserting cat3_minus_1 / cat3_plus_1 winner
   identity.

**If the user instead wants bundle regen first to confirm cat3_plus_2
consensus picks mapPacBio:** that's a single command per
`methods/rectify_quickstart_h2.md` and doesn't gate the reanchor work.

---

## 5. Files touched

**CODE (committed):**

- `rectify/core/splice/splice_aware_5prime.py` — re-enabled + strand
  equivalence-extension with correct undershoot trigger
  (`_intron_end - read.reference_start > 0`). Lines 1374-1416.
- `rectify/core/bam/read_edits.py` — two changes:
  - Added `reanchor_5prime_for_rescue` at line 38 (DEFERRED docstring
    notes it is not wired).
  - Fixed `reroute_intronic_tail_5prime_via_junction` + strand
    off-by-one in `intron_len` and `new_ref_start` (lines ~1014-1020
    of the - strand vs + strand symmetric block).
- `rectify/core/correct/indel_corrector.py` — added 2-bp deletion
  extension post-pass to `rescue_softclip_at_homopolymer`'s left-side
  branch (after the existing strip-trailing-T cleanup, before the
  result dict).

**TEST (committed):**

- `tests/test_validation_reads.py` — two changes:
  - New test
    `test_cat3_plus_2_rescued_aligners_have_clean_intron_cigar` at
    line 681, mirroring the cat3_minus_2 cigar test.
  - Updated `test_3prime_exact_position[cat2_minus_2-...]` from 128102
    to 128096 at line 497.

**DOCS (committed):**

- `docs/handoffs/debugger_queue.md` — two top-level entries added:
  - "Bug note: + strand 5'-rescue equivalence-extension geometry
    inverted (DISABLED 2026-05-18)" (predecessor session).
  - "Deferred: mapPacBio 5'-edge reanchor for rescue (function
    shipped, hooks reverted)" — this session. Documents the cross-
    stage CIGAR inconsistency diagnosis and two unhook paths.
- This file (`handoffs/HANDOFF_2026-05-18_6943450_cat3p2_cat2m2_reanchor.md`).

**`[uncommitted]` — working tree items present BEFORE this session and
not changed here:**

- `CLAUDE.md`, `HANDOFF.md`, `README.md`, `docs/figures/*` —
  pre-existing WIP.
- `rectify/core/commands/correct_command.py`,
  `rectify/core/splice/junction_refiner.py` — pre-existing WIP from the
  parallel UMI-bin Phase D work track.
- `rectify/data/validation/aligners/validation_reads.{deSALT,minimap2,uLTRA}.bam`
  — pre-existing WIP from before the `0653172` session.
- Various data/validation/PROVENANCE.json, README.md, and bundle
  output files — touched by parallel work tracks, not this session.
