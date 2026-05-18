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

## Cat2 cluster (HP-ED penalty calibration + soft-clip rescue lookahead)

Plotter findings, 2026-05-18, from
`validation_read_review/cat2_softclip_findings.md`.

### cat2_plus_1 — HP-ED penalty calibration: del_cost(long-A-run) may be too expensive

**Panel (verbatim):**
```
aligner     5'    Δ   3'    Δ    HP-ED  span  pick
minimap2    23362 0   23759 0    20.1   398
mapPacBio   23362 0   23754 0    16.9   393
deSALT      23362 0   23754 0    15.6   393   ← winner
uLTRA       23362 0   23759 0    20.1   398
gapmm2      23362 0   23759 0    20.1   398
```

**Two-cluster divergence:** minimap2/gapmm2/uLTRA at 3'=23759 (Cluster A);
mapPacBio/deSALT at 3'=23754 (Cluster B). HP-ED picks deSALT (Cluster B).

**User's diagnosis (verbatim):**
> "For cat2 plus1, I actually think minimap2 did the right thing. The
> undercalling of the A-tract is the most plausible alignment here, and it
> allows everything else to cleanly align. The extra AAAT repeats are
> classic nanopore over-calling of tetramer repeats. The hp-aware edit
> distance should be allowing minimap2 to win this one outright."

**Representations:**
- minimap2: 1× 8-bp HP deletion in a long A-tract (1 biological event).
- deSALT: 3× single-base insertions in AAAT tetramer repeats (3 events).

**Hypothesis:** `del_cost(hp_len ≈ 19, base='A')` in
`rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv`
may be too expensive — empirically nanopore DRS undercalls long A-runs at
high frequency, so cost should be very cheap. If del_cost were correct,
minimap2's HP-ED would dip below deSALT's.

Out-of-scope ideas:
- HP-ED could normalize per "biological event" rather than per CIGAR-op
  base (one 8-bp HP under-call is one event, not eight).
- Correlated insertion penalty for AAAT-tetramer over-calls (3 events at
  the same boundary are correlated, not independent).

### cat2_minus_2 — soft-clip rescue stops 3 nt short of a TGC body match

**Panel summary:** 4 aligners agree corr_3p = 128102 (T on RNA, non-A,
policy satisfied). mapPacBio uniquely stops at 128106.

**User's diagnosis (verbatim):**
> "cat2 minus 2 is tricky, but I suspect the right call would be to allow
> the TGC to bind and accept a T deletion as part of the massive A
> homopolymer undercall for this read. Allowing for 3 straight nt to match
> should outweigh the single T deletion here."

**Pattern:** trailing soft-clip `GCAAT` has a `TGC` that maps to upstream
genome at ~128097-128099. To incorporate those 3 matches, the rescue
would need to accept ONE additional T deletion in the large A-HP
undercall. Parsimony: 3 new matches − 1 new del = +2 net.

**Proposed fix:** extend Module 2G (soft-clip rescue,
`rectify/core/correct/indel_corrector.py:rescue_softclip_at_homopolymer`)
with a "cheap-HP-del lookahead":
1. After the standard rescue completes, take the soft-clip suffix.
2. Try to align the first ~10 bp of the soft-clip against the next
   ~10–20 bp of upstream genome (allowing a single HP del slip).
3. Accept the extension if it produces ≥3 consecutive matches with ≤1
   cheap HP del. Reject if the cost exceeds 1.5× the gain.

### cat2_plus_2 — clean; surfaces an effective-utility feature request

Correction looks correct. The interesting observation is a **two-cluster
pattern** within the per-aligner picks:
- Cluster A (3'=8614, Δ=0): mapPacBio + uLTRA.
- Cluster B (3'=8605, Δ=+4): minimap2/gapmm2/deSALT.

uLTRA wins HP-ED within Cluster A. Both clusters have valid biology.

### Feature request: per-aligner effective-utility tracking

**User's verbatim ask (2026-05-18):**
> "I'm wondering if we should have a separate column that acknowledges
> which aligners effectively got the same 5′, 3′ and junctions (if
> applicable). This is useful for knowing the 'effective' utility of each
> aligner for rectify. Presumably, uLTRA did a cleaner job of aligning the
> mRNA body that was more favorable with our empirical penalty table.
> Both analyses are useful. The 'effective' utility should be something
> that is tracked at read-level and sample-wide as well."

**Implementation sketch:**
- **Read-level**: a new `effective_group` column in
  `merge_corrected_tsvs`'s summary TSV. Cluster per-aligner rows by
  `(corrected_5p, corrected_3p, frozenset(corrected_junctions))`; assign
  each distinct tuple a letter (A, B, C, …).
- **Sample-wide**: aggregate per-aligner `effectively_matched_winner`
  boolean — "this aligner's 5'/3'/junctions matched the winning
  cluster's." Surface via a `--summary-stats` flag on
  `merge_corrected_tsvs` that prints the rollup.

This is genuinely interesting because the current `winning_aligner`
column hides the case where multiple aligners are in the right biological
cluster but lose HP-ED tiebreaks within it.

---

## Cat3_plus_1 reproduction — reroute_intronic_tail_5prime_via_junction silent False (HANDOFF Open §3)

**Reproduced 2026-05-18:**

```
TSV: 5p_rescued=1 fpp=168423 exon_cigar=4M1I5M icp=168808
     corrected_3prime=169459 n_junctions=1 junctions=168424-168808
BAM: cigar head=1X2=7I86=2I65=1X210=... (NO N-op)
     Xj tag: NOT set
```

mapPacBio raw CIGAR has no soft-clip → `extend_read_5prime_for_junction_
rescue` skips (requires `cigar[0][0]==4`). Should fall into
`reroute_intronic_tail_5prime_via_junction` (read_edits.py:816).
The function has multiple silent `return False` paths; need to add
logging to identify which precondition fails for cat3_plus_1.

**Likely culprit (HANDOFF):** query-length sanity (line 898:
`if n_intronic_q != exon_q_bases: return False`) or boundary check
(line 941: `if read.reference_start >= clip_boundary: return False`).

**RESOLVED 2026-05-18 (commit `62241ea`):** mapPacBio cat3_plus_1's
raw CIGAR is `1X 2= 7I 86=...`. The intronic head walk hits
`break` at `ref_pos >= clip_boundary` before consuming the 7I op
at the boundary, so `n_intronic_q=3` while exon_cigar `4M1I5M`
needs 10 query bases. Fix: after the loop hits the clip boundary,
also consume any boundary-adjacent query-only ops (I, S) so the
strict `n_intronic_q == exon_q_bases` check passes. Applied to
both plus and minus strand branches. Verified: cat3_plus_1
corrected CIGAR now `4M1I5M385N86=2I65=...` with the N-op
canonical at (168424, 168808).

## Cat3_plus_2 downstream CIGAR mangling — RESOLVED transitively

HANDOFF §3 NEW BUG flagged that refiner produces clean
`14=1D9=366N50=...` but downstream pipeline mangles upstream to
`22D1M21I1M`. Verified 2026-05-18 that with `--aligner-bams` the
corrected BAM is clean `14=1D9=366N50=1I2=...` — N-op canonical at
(142253, 142619). Likely fixed transitively by Phase A's
`_decode_eq_seq_inplace` (which keeps indel_corrector and friends
from re-rewriting CIGARs on `=`-encoded SEQ reads).

## Penalty table calibration — INS AT hp_len ≥ 4 missing (Cat2 cluster)

Investigation 2026-05-18 on `rectify/data/genomes/saccharomyces_
cerevisiae/penalty_tables/penalty_scores.tsv`:

```
I AT 1  0.0000  count=1.0   penalty=1.2500  low_count=True
I AT 2  0.0000  count=3.0   penalty=0.3805  low_count=True
I AT 3  0.0000  count=1.0   penalty=0.7156  low_count=True
(no entries for hp ≥ 4 — falls back to hp=3 value via
 PenaltyTable.ins_cost line 274)
D AT 15 0.9951  count=34922 penalty=0.2736  low_count=False
D AT 16 0.9951  count=25418 penalty=0.2736  low_count=False
...
D AT 18 0.9951  count=2686  penalty=0.2736  low_count=False
```

**Observations:**
1. INS AT hp ≥ 4 has NO entries. `ins_cost(hp=4, 'A')` falls back to
   `table[3] = 0.7156` (which itself is `low_count=True`, count=1).
   The cat2_plus_1 deSALT case has 3 single-base insertions in AAAT
   tetramer regions — each costs ~0.72 → total ~2.16.
2. DEL AT hp 15-20 has solid empirical entries (count > 2K) at
   0.27/base. minimap2's 8-bp deletion in a 19-A run costs ~2.19.
3. The two are roughly equal — but the cat2_plus_1 panel shows
   deSALT winning by 4.5 hp_ed (15.6 vs minimap2 20.1). The
   remaining gap (~4-5) is likely from elsewhere (hardclip / X-ops
   in the body alignment), not the indel cost itself.

**Follow-up:** extend the empirical INS AT table to cover hp_len
4-20 (recalibration job). Until then, INS-heavy alignments at
moderate-length HP runs may be either over- or under-penalized
relative to the empirically-calibrated DEL cost. This is a known
gap rather than a code bug.

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
