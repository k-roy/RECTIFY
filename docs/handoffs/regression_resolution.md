# Validation suite regression resolution log

Tests in `tests/test_validation_reads.py` were written against a bundle whose
consensus output was produced by an older pipeline. The new HP-mode merge
(`per_aligner_corrected_bams` threaded through production, landed in
`22136af`) picks different aligners for some reads — usually equivalent
biological outcomes via a different mechanism. This file documents which
tests were updated and why, per CONTINUE.md §6 ground rule 4 ("if a
regression is genuinely 'the new pipeline is correct, the old bundle was
wrong' — update the test, but also document why").

## Cat6 — `test_no_five_prime_rescue` replaced with `test_intron_present_at_expected_coords`

**Old assertion:** `row['five_prime_rescued'] == '0'`.

**Old rationale (per docstring):** "Cat6 reads come from a single winning
aligner (mapPacBio), so... no 5' rescue should be needed." The bundle's
consensus pre-HP-mode picked mapPacBio (native intron span). mapPacBio
spans the 5' intron with an explicit N-op; rescue isn't needed.

**New behavior:** HP-mode consensus picks deSALT (or another rescue-applied
aligner) for all 4 Cat6 reads. Those aligners had a 5' soft-clip in their
raw CIGAR; `rescue_3ss_truncation` correctly recognized the annotated 5'
intron and materialized the N-op via `extend_read_5prime_for_junction_rescue`.
The corrected BAM ends up with the **same** N-op at the **same** coordinates
— just produced by a different path (rescue surgery instead of native
intron span).

**Per-aligner inspection (all 4 Cat6 reads):**

```
cat6_plus_1 (875a773c) at chrII:125154-125270:
  minimap2/uLTRA/deSALT  → 5p_resc=1 n_junc=0+1(post-fix)  five_prime_rescued
  mapPacBio              → 5p_resc=0 n_junc=1              no rescue (native N)
  MERGED winner: deSALT (rescue-applied)
```

(Pattern is identical for cat6_plus_2, cat6_minus_1, cat6_minus_2.)

**Why the new behavior is correct:** the rescue-applied corrected BAM
contains the same N-op as mapPacBio's native span. The annotated 5' intron
is recovered. `five_prime_rescued` records *how* the N-op got there (path
metadata), not whether the read has the right biology. The test should
check the BIOLOGICAL OUTCOME, not the path.

**Fix:** `test_no_five_prime_rescue` replaced with
`test_intron_present_at_expected_coords` which asserts the corrected output
contains the annotated 5' intron at the expected coordinates. Both rescue-
applied and native-span winners satisfy this — path-agnostic.

A complementary fix in `bam_processor.py` (commit `dc5591e`) adds the
rescued junction to the TSV's `junctions` list so `n_junctions` reflects
post-rescue state, which the new test (and `test_has_one_junction`)
depend on.

## Cat6 + Cat7 — `test_xu_tag` relaxed from `Xm == 1` to `Xm >= 1`

**Old assertion:** `r.get_tag('Xm') == 1` (for source-BAM reads).

**Old rationale:** "Cat6/Cat7 reads come from a single winning aligner
(mapPacBio), so Xm=1." The legacy chimeric consensus always picked one
aligner per read in cases where mapPacBio's intron span was strong.

**New behavior:** The current chimeric consensus stitches multiple
aligners when their alignments are equivalent — Xm == 2 is common for
reads where minimap2 and mapPacBio agree on intron coordinates and both
contribute to the multi-aligner stitch. Both single-aligner and multi-
aligner stitches are biologically correct.

**Per-read observation** (cat6_plus_1 bundle source BAM): `Xm == 2` —
two aligners contributed to the stitch.

**Fix:** Relax to `xm >= 1` (presence + positive count). The semantic
meaning ("at least one aligner contributed") is preserved without
over-constraining on the exact contributor count.
