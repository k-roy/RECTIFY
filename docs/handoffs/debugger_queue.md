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

### cat2_minus_2 — soft-clip rescue: extension geometry is ambiguous; user input needed before code

**Panel summary:** 4 aligners agree corr_3p = 128102 (T on RNA, non-A,
policy satisfied). mapPacBio uniquely stops at 128106.

**User's diagnosis (verbatim, 2026-05-18):**
> "cat2 minus 2 is tricky, but I suspect the right call would be to allow
> the TGC to bind and accept a T deletion as part of the massive A
> homopolymer undercall for this read. Allowing for 3 straight nt to match
> should outweigh the single T deletion here."

**Genome + RNA geometry (investigated 2026-05-18 next session, pre-code):**

Genome + strand, chrI[128095:128117]:
```
pos    : 128095 96 97 98 99 100 101 102 103 ... 115 116
+strand:    T   T  T  G  C  A   T   A  T-tract  T  G
```

The long T-tract is 128103-128115 (13 bp), bracketed by A(128102) and
G(128116). RNA = revcomp, so the RNA-orientation 3'-end region reads
`...AAAAAAAAAAAAA-T-A-T-G-C-A-A-A-A...` going 5'→3' (from genome high
position downward).

Raw alignment (minimap2/gapmm2/deSALT/uLTRA): rs=128112, 5p softclip = 6
bases `ATTGCA` in BAM orientation (= `TGCAAT` in RNA 5'→3' orientation).

Current rescue trace (`rescue_softclip_at_homopolymer`, left side):
- `boundary_base = T` (locked to `genome[raw_pos]=128112` = T); the loop
  walks leftward as long as `genome == T`, so it cannot slip past the
  A at 128102. `homopolymer_extension = 9` (positions 128103-128111).
- `match_start = raw_pos - 9 - 1 = 128102`.
- Walk softclip innermost→outermost (BAM indices 5→0): softclip[5]=A vs
  genome[128102]=A → MATCH (rescue=1). softclip[4]=C vs genome[128101]=T
  → MISMATCH → terminate.
- Result: 1M at 128102, 9D over T-tract, 7M from 128112 onward.
  corrected_3p = 128102. Hard-clip retains softclip indices 0-4.

**What a lookahead-past-HP-del would actually capture:**

If the rescue could "skip" 2 ref bases (genome[128101]=T and genome[128100]=A)
after the matched A at 128102, the next softclip bases would line up:
- softclip[4]=C vs genome[128099]=C → match
- softclip[3]=G vs genome[128098]=G → match
- softclip[2]=T vs genome[128097]=T → match
- softclip[1]=T vs genome[128096]=T → match
- softclip[0]=A vs genome[128095]=T → mismatch (terminate)

That's **4 more matches, with a 2-bp ref skip** (one T at 128101 + one A
at 128100). corrected_3p would become **128096** (or 128099 if we anchor
on the TGC alone). NOT 128097-128099 (those are the matching positions,
not the corrected endpoint).

**Why the queue note's wording was imprecise:**

- "ONE additional T deletion": the 2-bp ref skip includes a non-HP A at
  128100, not just one extra HP-T. Framing it as "more HP-undercall" is
  only correct if we accept that the user's mental model lumps the A at
  128100 into the homopolymer for parsimony purposes.
- "3 straight nt to match": 4 softclip bases (TTGC) line up cleanly, not 3.
- "the TGC binds at 128097-128099": correct, but that's the inner part of
  the bound region. The full bound region is 128096-128099 (TTGC) with
  the outermost A at softclip[0] still unmatched.

**Open design question for the user:**

This shifts the consensus past the test's hardcoded 128102. Three
candidate endpoints, depending on how aggressive we let the rescue be:
1. **128099** (TGC only, 3 matches + 2-bp skip; the user's literal
   description). Net: +1.
2. **128096** (TTGC, 4 matches + 2-bp skip). Net: +2. Maximal HP-del
   parsimony.
3. **128095** (TTGC + try 1 mismatch extension). Probably not, since the
   T at 128095 mismatches the softclip A.

Before coding, need to know:
- Which endpoint does the user actually want?
- Is the rescue allowed to consume non-HP ref bases (the A at 128100) as
  part of a "cheap del" run, or strictly HP-T only? If strict, the 2-bp
  skip in this case isn't allowed and the answer is "current 128102 stands."
- Does this change the test (`test_3prime_exact_position[cat2_minus_2-128102]`)
  or do we want to gate the new behavior behind a flag?

**Files implicated when implementing:**
`rectify/core/correct/indel_corrector.py:rescue_softclip_at_homopolymer`
(both right-side and left-side branches need symmetric extension).
Test: `tests/test_validation_reads.py` line 487.

### cat2_plus_2 — clean; surfaces an effective-utility feature request

Correction looks correct. The interesting observation is a **two-cluster
pattern** within the per-aligner picks:
- Cluster A (3'=8614, Δ=0): mapPacBio + uLTRA.
- Cluster B (3'=8605, Δ=+4): minimap2/gapmm2/deSALT.

uLTRA wins HP-ED within Cluster A. Both clusters have valid biology.

### Feature request: per-aligner effective-utility tracking — RESOLVED (commit 75b0338)

**Status:** shipped 2026-05-17 in commit `75b0338` ("feat: bedgraph regen +
per-aligner effective-utility column"). The summary TSV now has
`effective_group` (read-level cluster letter) and `effectively_matched_winner`
(boolean per aligner-row); sample-wide rollup is emitted at `logger.info` from
`merge_corrected_tsvs`. The original ask, retained for context:

> "I'm wondering if we should have a separate column that acknowledges
> which aligners effectively got the same 5′, 3′ and junctions (if
> applicable). This is useful for knowing the 'effective' utility of each
> aligner for rectify. Presumably, uLTRA did a cleaner job of aligning the
> mRNA body that was more favorable with our empirical penalty table.
> Both analyses are useful. The 'effective' utility should be something
> that is tracked at read-level and sample-wide as well."

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

## Cat3 cluster (plotter session 2026-05-18)

Source: `validation_read_review/cat3_junction_findings.md`.

### cat3_minus_2 + cat3_plus_2 — asymmetric 2-bp slide producing `M 2D N` instead of clean `N` (PARTIAL: cat3_plus_2 already clean post-Phase-A)

**User's diagnosis (verbatim):**
> "We should be capturing the correct intron with a clean N-op, rather
> than 2D-Nop. All we need to do is move the AG aligned bases to exon 1.
> Also, it seems our 3' pileup signal is stale."

**Investigation 2026-05-18:** post-Phase-A/E.1 bundle regen via
`regen_pa_rest_bundle.py`. Per-aligner picture for cat3_minus_2:

```
minimap2/uLTRA/deSALT  N(366504,366584)  ... N80 2D 13M ...
mapPacBio              N(366502,366584)  ... N82 15= ...     (clean)
gapmm2                 N(366502,366584)  ... N82 15M ...     (clean)
```

3 of 5 aligners place the intron at (366504, 366584) with a trailing
`2D 13M`; 2 aligners place it at (366502, 366584) with clean `15=`/`15M`.
The legacy simple-slide fast path
(`junction_refiner.py:_apply_junction_replacement`, lines 513-590) only
handles the **symmetric** case `delta_start == delta_end != 0` (pure
slide, intron length preserved). For cat3_minus_2 the canonical placement
needs `delta_start = -2, delta_end = 0` — an **asymmetric** slide where
the intron *length changes* by 2.

cat3_plus_2 is already clean post-Phase-A:
```
mapPacBio    N(142253,142619)  14=1D9= N366 50=1I2=  (clean — slid by Bug 1)
others       N(142253,142618)  ...
```

So 4 of 5 aligners on cat3_plus_2 are off by 1 on the acceptor — but
mapPacBio (winner under HP-ED) has the canonical placement.

**Proposed fix (deferred):** add a new fast path for the asymmetric-slide-
with-length-change case. When the post-N CIGAR has trailing D ops with
no query consumption that could be absorbed into the N (extending the
intron) AND the k bases at the would-be new boundary equivalence-match
the bases at the current boundary, emit a length-extension CIGAR
transformation rather than letting the general path produce `M 2D N`.

The general criterion (user, verbatim):
> "`genome[old_donor:old_donor+k] == genome[old_acceptor:old_acceptor+k]`
> for k=1, 2, 3, …; when the equivalence holds, emit a length-adjustment
> CIGAR fix rather than a local realignment that produces `kD N`."

This is a meaningful architectural change to the refiner. Out of scope
for the current session.

### Stale 3' pileup bedgraphs (mechanical follow-up) — RESOLVED (commit 75b0338)

**Status:** shipped 2026-05-17 in commit `75b0338`. `regen_pa_rest_bundle.py`
now refreshes `rectified/corrected_3ends.{plus,minus}.bedgraph` after the
corrected BAMs + TSV are written, using `bedgraph_writers.py`.

### Summary TSV `junctions` column — RESOLVED transitively

Plotter found: rescue-applied rows showed `junctions=1` instead of
`donor-acceptor` coords. Resolved by my Phase D fix (commit `dc5591e`)
which appends the rescued junction to the per-aligner TSV's `junctions`
list. Current bundle's `rectified/per_aligner_summary.tsv` shows
`junctions=142253-142619` for all aligners on cat3_plus_2 (including
the rescued ones).

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
