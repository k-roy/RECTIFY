# Validation suite debugger queue

User-supplied findings for follow-up sessions. These are real bugs/design
questions that surfaced during Phase D investigation but were deferred to
keep the 40 → 0 push moving.

## Foundational policy: reads never end in A

Per user (2026-05-18, plotter session via
`validation_read_review/cat1_walkback_findings.md`):

> Reads can NEVER end in A in RECTIFY, as per policy we always walkback to
> the first non-A due to the inherent ambiguity at the genomic-A /
> poly(A)-tail boundary.

This is stronger than "leftmost-possible-CPA" — it is an *absolute terminal
base constraint* on the corrected 3' end. The base at the corrected
position must be non-A on the **RNA strand**:
- For ``+`` strand reads: ``genome[corrected_3prime]`` must not equal ``'A'``.
- For ``−`` strand reads: ``genome[corrected_3prime]`` must not equal ``'T'``
  (because RNA = revcomp of genome, so a genomic T on the minus-strand
  alignment is an A on the RNA).

`walkback_3prime_guarded` should be the single enforcement point for this
policy. If the walkback exits early (e.g. the 4-A homopolymer check at
walkback.py:437-449 doesn't trigger), the corrected position may end up at
an A on RNA — the regression on cat1_minus_1 in the
`rectified/per_aligner_summary.tsv` bundle reflected exactly this:
walkback short-circuited and left corr_3p on a chrII T (RNA A).

## Cat1 cluster (HP-mode metric / walkback design)

### cat1_plus_1 — mapPacBio force-aligns 4+ mismatches into the body

**Status:** currently failing `test_3prime_shifted` + `test_3prime_exact_position[cat1_plus_1-10611]`.

**Panel finding (user, 2026-05-18):** mapPacBio winner extends to chrXIV:10617
(vs 10611 for minimap2/gapmm2/uLTRA via `overcall_rescue`). Genome around:
`chrXIV[10600..10620] = AAAAAAAAAAATAGCTCTAT`. The body alignment ends with
poly-A tail crossing the genomic boundary at ~10611; mapPacBio forces 4+
mismatches `GCTC → ATAA` to anchor at a "non-A" past the tail.

**User's diagnosis (plotter session, verbatim):**
> "We are matching terminal non-A residues at the expense of too many
> mismatches. Over/under-calls of A-tracts are allowed, with perhaps a
> mismatch or two, but in this case we need to compare what's more likely:
> (a) two A → T nanopore sequencing errors that would mean the pA tail is
> simply longer, or (b) four consecutive errors GCTC → ATAA."

(a) is the more parsimonious answer; mapPacBio's choice of (b) loses.

**Two related bugs:**
1. **HP-ED scoring shouldn't reward mapPacBio's over-extension.** The H/S
   penalty doesn't catch mismatches inside `=`/`X` ops as expensively as a
   soft-clip extension would, so mapPacBio's `…1X13=4X1=` tail wins HP-ED
   over the rescue-applied `…{hp_ext}D{overcall_count}I 1=` tail.
2. **mapPacBio's extension itself shouldn't be produced.** Even if winner
   selection were fixed, mapPacBio independently produces this over-extended
   alignment — a real-data correctness issue beyond the validation suite.

**Walkback investigation (2026-05-18 from a debugger session):**

Tracing `walkback_3prime_guarded` on mapPacBio cat1_plus_1 shows the
**early-exit homopolymer check at walkback.py:437-449 returns None before
the main scan runs.** The check looks for ``stop_base × early_exit_min_
homopolymer_len`` (4 A's by default) within a narrow window:
``[_raw_3p - 1, _raw_3p + early_exit_min_homopolymer_len + 1)`` = 6 bp.
For mapPacBio's 3' end at 10617, the window is chrXIV[10616..10622] =
`CTATTC` — no `AAAA` → walkback exits early.

The genomic A-tract is at chrXIV[10600..10610] (7 bp upstream of the 3'
end). The narrow early-exit window misses it.

**If the early-exit were removed/widened, the main scan would handle this
correctly:**
- At i=10617 (rb=T, gb=T): candidate match. Tail-context-false-stop guard
  (lines 578-597) looks leftward: positions 10613-10616 have rb=A (stop_base)
  and rb != gb → 4 stop-base mismatches → `ctx_n >= tail_context_k (4)` AND
  `ctx_all_stop` AND `ctx_has_mismatch` → reject as false stop, continue.
- Scan continues leftward through the 4X region (no candidates), through the
  13= region (rb=gb=A matches; not anchors), reaches i=10611 (rb=T, gb=T,
  gb != A) → MATCH. Tail-context guard at 10611 looks leftward at A-run
  positions: `ctx_has_mismatch=False` → guard accepts → anchor at 10611.

**Proposed fix:** widen the early-exit window to cover the A-tract region
broadly (e.g. ±20 bp from the 3' end), OR add a SECONDARY entry condition
that triggers walkback when the trailing ≥k CIGAR positions have
`rb=stop_base` and `rb != gb` (force-alignment-past-pA-tail signature).

**Location:** `rectify/core/correct/walkback.py:437-449`.

### cat1_plus_2 — HP-aware insertion cost may not fire

**Status:** currently failing `test_3prime_shifted` + `test_3prime_exact_position[cat1_plus_2-31546]` (wait — Phase A fix should have addressed this; verify the latest pytest).

**Panel finding (user, 2026-05-18):** `corr_3p = 31545` on chrI. Genome
around: `chrI[31535..31560] = GTCACCGAAAAGAAAAGGTAAAAAG`. Genuinely ambiguous:
(a) `GAAAAG` → `GAAAA(A)G` with 1 bp HP insertion + 1 G→A mismatch, alignment
stops here; (b) extend further into `AAAAA G` with multiple errors. Treating
the first surplus A as `1I` (insertion into the HP run) is what the HP-aware
penalty table should favor.

**Proposed fix:** verify `indel_corrector` produces `…GAAAA1IG…`-style CIGARs
in cases like this. The HP-aware `ins_cost` lookup may not be firing where
expected.

**Location:** `rectify/core/correct/indel_corrector.py:correct_indels_from_read`
(insertion-detection branch). Check the penalty table lookup path.

### cat1_minus_1 — resolved by Phase A fix; bundle regen pending

**Status:** RESOLVED in the current code (Phase A commit a1728eb). The
plotter's panel finding ("walkback under-shoots by 1 bp; corr_3p=9831 lands
on chrII T = RNA A") reflected pre-Phase-A bundle data in
`rectify/data/validation/rectified/per_aligner_summary.tsv`. Pytest fixture
runs against the current code show `corr_3p = 9834` (chrII[9834] = C =
RNA G) — policy satisfied.

The `_decode_eq_seq_inplace` call at the start of `correct_read_3prime`
(Phase A) was the load-bearing change: minimap2/gapmm2/deSALT propagate the
SAM-spec `=` shorthand in SEQ at match positions, and walkback was comparing
`'=' != 'T'` as a mismatch, falsely bailing. Decoding at intake makes
walkback see the actual T bases and walk back the 3 bp correctly.

**Follow-up:** regenerate `rectify/data/validation/rectified/` so the
plotter's panel data reflects the new positions. The current
`per_aligner_summary.tsv` is stale.

### cat1_minus_2 — positive control

Correctly handled (GG rescued via 1bp del in the A4 tract). Use as a positive
control while debugging the other Cat1 reads. The corrected 3' end lands at
chrXII[15345] = C on plus = G on RNA → policy satisfied.

---

## Cat6 / Cat7 (HP-mode metric — Xm semantics)

Deferred — see `regression_resolution.md` for the test-side relaxation. The
underlying question (when should chimeric_consensus stitch 1 vs 2 aligners?)
is architectural and warrants its own session.

---

## Cat4 (false-junction filter — no exemplar)

Per-aligner BAMs in the current bundle don't reproduce the false-N pattern
the test was designed for. Either re-align with a configuration that produces
the false N (so the filter has work to do), or retire Cat4 from the suite
and replace with a new exemplar. Defer until Phase B bundle regen and a
fresh per-aligner inventory.

---

## Cat9 (Module 2H contract)

The legacy bundle had aligners producing imprecise raw N-ops that Module 2H
corrected. The current aligner pool produces canonical N-ops natively, so
the "before 2H" / "after 2H" distinction the tests check is no longer
exercised. Defer until Module 2H scope is revisited (user-deferred per the
plan).
