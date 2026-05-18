# Validation suite debugger queue

User-supplied findings for follow-up sessions. These are real bugs/design
questions that surfaced during Phase D investigation but were deferred to
keep the 40 → 0 push moving.

## Cat1 cluster (HP-mode metric / walkback design)

### cat1_plus_1 — mapPacBio force-aligns 4+ mismatches into the body

**Status:** currently failing `test_3prime_shifted` + `test_3prime_exact_position[cat1_plus_1-10611]`.

**Panel finding (user, 2026-05-18):** mapPacBio winner extends to chrXIV:10617
(vs 10611 for minimap2/gapmm2/uLTRA via `overcall_rescue`). Genome around:
`chrXIV[10600..10620] = AAAAAAAAAAATAGCTCTAT`. The body alignment ends with
poly-A tail crossing the genomic boundary at ~10611; mapPacBio forces 4+
mismatches `GCTC → ATAA` to anchor at a "non-A" past the tail. Better answer:
accept 2 A→T sequencing errors (poly-A tail is simply longer) and stop the
alignment shorter. HP-ED currently rewards mapPacBio for getting more aligned
bases despite the mismatches because the trailing run isn't being recognized
as poly-A tail.

**Proposed fix:** walkback / 3'-end extension heuristic — bound the number of
consecutive mismatches accepted when "extending toward a non-A anchor". When
mismatches > 2–3 and the read bases are all A, treat as poly-A tail
continuation, not as true body extension.

**Location candidates:**
- `rectify/core/correct/walkback.py:walkback_3prime_guarded` (the right-side
  scan that finds the first read-vs-ref non-stop-base match).
- `rectify/core/consensus/corrected_consensus.py:_cigar_hp_edit_distance`
  (the metric that lets force-aligned X-runs win).

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

### cat1_minus_2 — positive control

Correctly handled (GG rescued via 1bp del in the A4 tract). Use as a positive
control while debugging the other Cat1 reads.

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
