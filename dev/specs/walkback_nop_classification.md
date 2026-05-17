# Symmetric N-op handling in `walkback_3prime_guarded`

**Status:** open. Predicate work is committed (Gap 2 = QSrev routing +
opt-in right-side bridging guard). This brief replaces that hack with a
principled, side-symmetric design.

**Origin:** raised 2026-05-16 during Gap 2 review. The user
("Why are left and right sides treated differently? Sounds like a
fragile hack to skirt around the issue by just getting by these test
cases.") was right — the current asymmetry is data-driven heuristic,
not biology. This brief is the cleanup.

## Problem

After Gap 2, `walkback_3prime_guarded` has two N-op behaviors:

- **LEFT side (`intron_boundary_guard=True`, default on).** Scan is
  unconditionally clipped at the **first** N-op. Justification: small
  terminal exon fully consumed by poly-T tail would otherwise let the
  scan anchor spuriously in the next exon.
- **RIGHT side (`right_side_bridging_guard=False`, default off, opt-in
  via QSrev wrapper).** When enabled, scan is clipped at the
  **rightmost** N-op end if it lies within `poly_noise_window` of the
  3' end. DRS keeps this off so Cat4-style false-N reads
  (aligner-inserted intron in the polyA noise window) still get
  absorbed and walkback can anchor a few bp upstream in real gene-body
  sequence.

The two cases are mirror images alignment-wise but have **opposite**
desired walkback behavior because:

1. Real intron with poly-A-consumed terminal exon → don't cross; anchor
   at boundary.
2. Aligner-inserted false junction → cross; anchor in pre-intron
   gene-body sequence.

The walkback can't tell these apart from the read alone, so the current
code outsources the decision to the caller's *protocol identity*
(`right_side_bridging_guard=True` for QSrev, off for DRS). That's the
fragile hack.

## The principled fix

Move the artifact-vs-real classification **upstream of the walkback**
(it already lives at `rectify/core/splice/false_junction_filter.py::
filter_polya_artifact_junctions`) and pass the classification *into*
the walkback as data. The walkback then handles both sides
symmetrically:

- **Artifact N-ops** → cross silently (scratch array already skips them
  transparently; current right-side default).
- **Real N-ops** → clip at boundary (current left-side default).

Same code path, same biology, both sides.

## Proposed API

```python
def walkback_3prime_guarded(
    read: pysam.AlignedSegment,
    chrom_seq: str,
    three_prime_side: str,
    stop_base: str = "A",
    *,
    artifact_n_ref_starts: Optional[set[int]] = None,
    # ↑ ref_starts of N-ops the caller has classified as artifacts.
    #   None  → behave as if no classification was done; clip at every
    #           N-op on side='left' (preserves legacy intron_boundary_guard
    #           behavior) and at every N-op within poly_noise_window on
    #           side='right' (preserves Gap 2's right_side_bridging_guard
    #           behavior when it would have been enabled).
    #   set() → treat every N-op as REAL → clip at every boundary on
    #           both sides.
    #   {a,b,c} → these N-op ref_starts are artifacts (skip transparently);
    #           any other N-op is real (clip at boundary).
    large_del_min_bp: int = 5,
    poly_noise_window: int = 50,
    tail_context_k: int = 4,
    max_scan_depth: int = 1000,
    early_exit_homopolymer_check: bool = True,
    early_exit_min_homopolymer_len: int = 4,
) -> Optional[dict]:
```

Retire these parameters:

- `intron_boundary_guard` — semantics replaced by `artifact_n_ref_starts`.
- `right_side_bridging_guard` — same.

Keep a one-release deprecation alias if you want; the user has been
moving fast and may be fine with a hard cut since the function is only
called from two wrappers.

## What "clip at boundary" means, symmetrically

Both sides:

- For each **real** N-op encountered while building the scratch array,
  remember `(n_ref_start, n_ref_end, scratch_idx_just_after)`.
- `side='left'` scan: clip `scan_hi` to the smallest
  `scratch_idx_just_after` among real N-ops (= first real-intron exit
  on the leftmost-exon's right side). Fallback anchor at that N-op's
  `n_ref_start` if the pre-N exon was all stop-base — current behavior.
- `side='right'` scan: clip `scan_lo` to the largest
  `scratch_idx_just_after` among real N-ops within
  `poly_noise_window` of the 3' end (= last real-intron exit on the
  rightmost-exon's left side). Fallback anchor at that N-op's
  `n_ref_end` if the post-N exon was all stop-base.

Artifact N-ops contribute *nothing* to the scratch metadata — they're
skipped during scratch building (already the case) and ignored during
pre-scan.

## Wrapper changes

### `walkback_drs_full` (DRS production)

Currently in `rectify/core/correct/walkback.py:walkback_drs_full`:

```python
return walkback_3prime_guarded(read, chrom_seq, side, stop_base=stop_base, ...)
```

Becomes:

```python
artifact_starts = _classify_artifact_nops(read, chrom_seq, gene_strand)
return walkback_3prime_guarded(
    read, chrom_seq, side,
    stop_base=stop_base,
    artifact_n_ref_starts=artifact_starts,
    ...
)
```

where `_classify_artifact_nops` is a thin wrapper around
`filter_polya_artifact_junctions(read, _GenomeDictReference({chrom: chrom_seq}), gene_strand)`
that returns `{analysis.junction_start for analysis in artifact_analyses}`.

### `walkback_quantseq_rev`

Same change. QSrev's dT-primer bridging artifacts will be detected by
`filter_polya_artifact_junctions` as artifacts (their downstream is
A-rich, target is A-rich, no canonical motif). So they end up in
`artifact_n_ref_starts` and the walkback skips them — same outcome as
the current `right_side_bridging_guard=True` path, but reached by
classification, not by protocol identity.

### `bam_processor.py:correct_read_3prime`

Currently calls `filter_polya_artifact_junctions` at line 355 and
walkback (DRS or QSrev) ~300 lines later. The classification is
*already* available before the walkback call. Plumbing options:

- **Cleanest:** pass `_artifact_analyses` into the walkback wrappers
  (`walkback_drs_full(read, _chrom_seq, artifact_analyses=_artifact_analyses)`).
  Wrappers extract `junction_start` and pass to the guarded core.
- **Simplest:** keep wrappers self-contained — they call
  `filter_polya_artifact_junctions` themselves. Slightly redundant
  work in `bam_processor` (filter runs twice per read) but no API
  contortions. Probably fine — the filter is cheap relative to
  walkback scratch building.

Pick "cleanest" if it doesn't require changing 4+ call sites; "simplest"
otherwise.

## Test gate

The classification-based behavior must produce **the same** corrected
positions as today on:

- `tests/test_validation_reads.py` — all 113 tests, especially:
  - `TestCategory4FalseJunction` — Cat4 artifact-N must still be crossed.
    cat4_plus_2 → 393721 (cross 100 bp artifact-N at 393725–393825).
  - `TestCategory2InternalPriming` — normal walkback (no N-ops involved).
  - `TestCategory5LongIntron` — real introns, walkback must NOT cross.
- `tests/test_quantseq_rev_walkback.py` — QSrev wrapper.
  - `TestQuantSeqRevGuardedRouting::test_nop_bridging_returns_no_correction`
    — QSrev bridging N must still be blocked.
  - `test_strand_inverted_gapped_read_walkback` — synthetic intron with
    non-A post-intron content, anchor must land in post-intron at 1059.
    (This test uses synthetic genome `X*1000 + ACGTC + X*50 + ACGTC + A*8 + X*100`,
    so `filter_polya_artifact_junctions` won't classify it cleanly —
    consider mocking the classifier or adjusting the test to call the
    guarded core directly with a chosen artifact set.)
- `tests/test_walkback_readvsref.py::TestGuardedParityWithFindPolyaBoundary`
  — DRS validation BAM parity. Must remain byte-identical.

The Cat4 tests are the load-bearing gate. They were the original
motivator and the one that caught Gap 2's regression.

## Critical pitfalls

1. **Classifier reliability.** `filter_polya_artifact_junctions` uses
   downstream A-richness, target A-richness, and canonical motif
   (GT-AG). It can misclassify:
   - **False positive** (real intron tagged as artifact): walkback
     crosses where it shouldn't. Result: walkback anchors too far
     upstream. Mitigation: the classifier already returns conservative
     results — a real intron without A-rich downstream stays
     classified as real.
   - **False negative** (artifact tagged as real): walkback blocks
     where it should cross. Result: no correction applied (corrected_pos
     stays at original). This is the *conservative* failure mode.
   Both modes were already in play under the existing Gap 2 design —
   the new design isn't introducing new failure modes, just redirecting
   their source from "protocol flag" to "classifier accuracy."

2. **LEFT-side behavior change.** The current LEFT-side guard fires
   unconditionally at the first N-op. With classification, a LEFT-side
   N-op classified as artifact would be crossed. **Verify Cat4-minus**
   (cat4_minus_1, cat4_minus_2) tests still pass. If LEFT-side artifact
   crossing is undesired biologically (no known production scenario),
   it's defensible to keep LEFT clip-at-all-N-ops by passing an empty
   artifact set on LEFT — but that's the asymmetry-revival path.

3. **Synthetic test reads.** Tests in `test_walkback_readvsref.py` and
   `test_quantseq_rev_walkback.py` use synthetic genomes
   (`X * 1000 + ACGTC + X*50 + ...`). The classifier's A-richness/motif
   checks may not behave usefully on these. Options:
   - Call the guarded core directly with a chosen
     `artifact_n_ref_starts` (bypass the classifier).
   - Construct synthetic reads with realistic flanking A-tracts so the
     classifier gives the intended verdict.
   The first option is preferred for unit-level tests; the second is
   more end-to-end but harder to control.

4. **The `intron_boundary_guard` fallback (LEFT side anchor at
   first_n_start when terminal poly-T zone has no non-T match).** This
   fallback should still apply when the LEFT-side N-op is classified
   as real. Keep the fallback wired in the new design — it's
   independent of the artifact-vs-real distinction.

5. **`walkback_3prime_with_qpos`** (line 189) is the cDNA-pipeline
   variant. It doesn't currently route through the guarded core. Out
   of scope for this brief, but worth noting: the cDNA path will keep
   its current "no guard" behavior until someone unifies it with the
   guarded core. The user has flagged that as separate cdna walkback
   work.

6. **Migration order.** Do this *after* `bam_processor_split.md` if
   you're running both. The classification plumbing is much cleaner
   when `correct_read_3prime` is already extracted into its own file.
   But not strictly required — both refactors are mostly orthogonal.

## Files touched

- `rectify/core/correct/walkback.py` — `walkback_3prime_guarded`
  signature + scratch metadata + pre-scan logic. Retire
  `intron_boundary_guard` and `right_side_bridging_guard`. ~150 lines
  changed.
- `rectify/core/correct/walkback.py` — `walkback_drs_full` and
  `walkback_drs` wrappers: add classification call. ~40 lines.
- `rectify/core/correct/protocols/quantseq_rev.py` — same. Drop the
  `right_side_bridging_guard=True` argument. ~15 lines.
- `rectify/core/correct/indel_corrector.py` — `find_polya_boundary` shim:
  update to new signature. ~5 lines.
- `rectify/core/bam/bam_processor.py` — if you go "cleanest" plumbing,
  pass `_artifact_analyses` into the walkback wrappers. ~10 lines.
  Otherwise no change.
- `tests/test_walkback_readvsref.py` — update any test that called
  `walkback_3prime_guarded` with `intron_boundary_guard=` or
  `right_side_bridging_guard=` directly. ~30 lines.
- `tests/test_quantseq_rev_walkback.py` — same. ~15 lines.

## Commit shape

One commit, conventional-commits style:

```
refactor(walkback): replace side-asymmetric N-op flags with artifact-set classification

Retires `intron_boundary_guard` (always-on LEFT) and
`right_side_bridging_guard` (opt-in RIGHT) in favor of a single
`artifact_n_ref_starts` parameter that both sides consume
symmetrically. Caller classifies each N-op via
filter_polya_artifact_junctions; walkback transparently skips artifact
N-ops and clips at real ones — same code path on both sides.

Cat4 false-junction reads (artifact-N within poly-A noise window) and
QSrev dT-primer bridging are now distinguished by biology, not by
protocol identity.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```

## Estimated scope

~250 lines net change across 5–6 files. One-CLI-session task assuming
the test gate doesn't surface surprises in the classifier's
production-data behavior. If `TestGuardedParityWithFindPolyaBoundary`
diverges on validation reads, that's the signal that the classifier
is non-trivially affecting DRS behavior and the brief becomes a
2-session task (one for the refactor, one for re-tuning classifier
thresholds against validation_reads.bam).
