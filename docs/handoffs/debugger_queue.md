# Validation suite debugger queue

User-supplied findings for follow-up sessions. These are real bugs/design
questions that surfaced during Phase D investigation but were deferred to
keep the 40 → 0 push moving.

---

## Bug note: + strand 5'-rescue equivalence-extension geometry inverted (DISABLED 2026-05-18)

**Scope:** `rectify/core/splice/splice_aware_5prime.py`, `rescue_3ss_truncation`,
the equivalence-extension block added in commit `0653172`.

**What was wrong:** The + strand branch triggered on
`_overshoot = read.reference_start - _intron_end > 0` (body M starts PAST
canonical intron_end, i.e., a gap between intron and body M). The
structural mirror of the - strand `ref_end > intron_start` overshoot is
the OPPOSITE sign: `_intron_end - read.reference_start > 0` (body M
extends INTO the intron from the high-coord edge).

**Why it produced wrong BAMs (if it had ever fired):** When the + strand
branch triggered, the downstream surgery in
`extend_read_5prime_for_junction_rescue` (read_edits.py) computed:

```
intron_len = (read.reference_start + effective_trim) - five_prime_position - 1
           = (intron_end + k + k) - (intron_start - 1) - 1
           = canonical_intron_len + 2k
```

— a non-canonical intron 2k bases longer than the annotated splice junction,
spanning `[intron_start, intron_end + 2k)`. The - strand path's arithmetic
collapses cleanly to canonical (`canonical_intron_len`); the + strand path
doesn't because the geometry is inverted.

**Independent giveaway:** The + strand case the code handles has no D op
to absorb. Body M's first k bases at `[intron_end + k, intron_end + 2k)`
are pure matches inside the downstream exon. The transformation just
moves k bases from body M to rescue M and lengthens the intron — no
D-elimination benefit, which is the whole point of equivalence-extension.

**Verification trail:** No real or synthetic + strand cat3 read in the
validation bundle has the overshoot pattern needed to trigger the buggy
path (cat3_plus_2 is undershoot; cat3_plus_1 has `ref_start == intron_end`
exactly). Pytest passed because the path was never exercised. The
handoff `HANDOFF_2026-05-18_0653172_cat3_equivalence.md` §2 already
flagged this as "+ strand mirror code path has no test that exercises
it" — the verification follow-up caught the geometric inversion.

**Current state:** The + strand `elif` branch is replaced with an inline
explanatory comment; behavior is now identical to pre-`6d2cf59` on
+ strand rescues. Validation suite still 106 passed / 8 skipped.

**To fix properly:** Implement the correct mirror, which is structurally
the same transformation as the deferred cat3_plus_2 "off-by-1 acceptor"
undershoot case (§3 of the cat3 handoff). Sketch:

- Trigger: `strand == '+' AND _intron_end - read.reference_start > 0`
  (body M starts k bases short of canonical intron_end, INSIDE intron).
- Borrowed query bases: first `_k_try` body-M query bases (same as
  before): `_q[_scl : _scl + _k_try]`.
- `_ref_old = genome[read.reference_start : read.reference_start + _k_try]`
  — currently aligned at `[intron_end - k, intron_end)` inside intron region.
- `_ref_new = genome[intron_end : intron_end + _k_try]` — downstream
  exon start (structural mirror of - strand's `_ref_new`).
- Surgery: the borrowed bases need to move from intron-region to
  downstream-exon-start. That is NOT the same operation as appending to
  rescue M (which is on the upstream side). On + strand, "move to
  downstream exon start" means shifting body M's ref_start LEFT by k
  (from `intron_end - k` to `intron_end - 2k`? — re-derive carefully).
  This may need a different BAM-writer code path than the - strand
  case, since both - and + strand currently route the "borrowed bases"
  into the rescue M, which is the upstream side for + strand but the
  downstream side for - strand. The asymmetry probably needs a separate
  surgery routine.

**Blast radius:** The handoff §3 cat3_plus_2 case (4/5 aligners) is
exactly this undershoot pattern. The proper + strand fix would clean
up its cosmetic CIGAR the same way - strand's was cleaned for
cat3_minus_2.

---

## Design note: splice-junction ambiguity window + motif-strength tiebreaker (next session)

**User ask (2026-05-18 evening):** "We check upstream of the called 5' SS
and 3' SS for the same nt(s). We also check downstream for the same nt(s).
We report back those numbers as the ambiguity window, and use the splice
site signal motif strengths (defined in one of our rectify docs) if the
window is >0 for tiebreaking purposes."

**Geometric criterion (symmetric slide — already what the existing fast
path in `junction_refiner.py:_apply_junction_replacement` lines 513-590
implements):**

For an intron at (intron_start, intron_end):
- `up_amb(k)` holds when `genome[intron_start - k : intron_start] ==
  genome[intron_end - k : intron_end]` (upstream of each SS matches).
- `down_amb(k)` holds when `genome[intron_start : intron_start + k] ==
  genome[intron_end : intron_end + k]` (downstream of each SS matches).
- The ambiguity window is the largest k where either holds → intron can
  be slid by up to k bp in that direction without changing length.

(This is distinct from the asymmetric-slide criterion
`genome[old_ns - k : old_ns] == genome[old_ne : old_ne + k]` shipped this
session in `rescue_3ss_truncation`. That one changes intron *length*.)

**Existing infrastructure (already in tree):**

- `rectify/utils/splice_motif.py`:
  - `SpliceMotifScorer.score_five_ss(seq) -> float` — yeast penalty matrix
    `[4, 3, 1, 1, 2, 1]` against consensus `GTATGT`. Lower = better.
  - `SpliceMotifScorer.score_three_ss(seq)` — `[1, 1, 1, 3, 4]` against
    `YYYAG`.
  - `get_splice_site_dinucleotides(genome, chrom, intron_start, intron_end,
    strand)` — returns (five_dinuc, three_dinuc) in transcript orientation.
  - `get_splice_site_sequences(genome, chrom, intron_start, intron_end,
    strand, context=6)` — returns the consensus-length sequences for
    scoring.
- `rectify/core/consensus/consensus.py:load_annotated_junctions(annotation_path)`
  — returns `Set[(chrom, intron_start, intron_end, strand)]`.

**CLAUDE.md line 202-204** already prefigures this work:
> "Future enhancement: pre-compute `up_amb` / `down_amb` fields on
> annotated junctions so soft-clip rescue can flex the match length
> within the ambiguity window."

**Implementation sketch (suggested 4 commits):**

1. **Per-junction ambiguity precompute** (no behavior change yet):
   - In `consensus.py:load_annotated_junctions`, augment the return type
     from `Set[tuple]` to `Dict[tuple, dict]` where each value carries
     `up_amb`, `down_amb`, `five_ss_motif_score`, `three_ss_motif_score`.
   - Compute by walking outward from each SS while genome bases match
     between the two SS (cap at, say, k=10).
   - Backward-compatible callers: those that just need the set of
     coordinates can still iterate keys.

2. **TSV reporting columns** (observability):
   - Add `junction_up_amb`, `junction_down_amb`, `junction_five_ss_score`,
     `junction_three_ss_score` to `corrected_reads.tsv` (one set of values
     per row; for multi-junction reads, use the worst-scoring junction or
     join all). Plumbing parallels what this session did for
     `five_prime_upstream_trim` — through `bam_processor.py`/`output.py`.

3. **Tiebreaker integration** (winner selection):
   - In `corrected_consensus.merge_corrected_tsvs`, when the HP-edit
     distances are tied AND the tied aligners disagree on junction
     placement within an ambiguity window of an annotated junction, prefer
     the placement with the lower motif penalty score. Concretely: extend
     the sort key tuple at line 765 with a motif-score component that
     activates only on `up_amb + down_amb > 0`.
   - The advisor noted (2026-05-18 session) that the blast radius of any
     tiebreaker change is narrow (1-3 reads in validation; check production
     before assuming impact).

4. **Tests** in `tests/test_validation_reads.py`:
   - `test_annotated_junction_carries_ambiguity_fields` — every row whose
     `junctions` column refers to an annotated junction should have non-empty
     ambiguity/motif columns.
   - `test_tied_hp_ed_resolves_by_motif_score` — for a curated read where
     two aligners tie on HP-ED but place the junction at different points
     in the ambiguity window, the higher-motif-score placement wins.

**Out of scope for this design note:**

- Whether to use motif-strength for ANY tiebreaking (not just ambiguity-
  window-driven). The user's verbatim was specific to "if the window is >0
  for tiebreaking purposes" — keeping that scope is safe.
- Whether to use motif-strength as a NOVEL junction filter. That's a
  separate concern that interacts with the chimera-exemption discussion
  in the cat3 entry below.

---

## Original queue (pre-2026-05-18 evening)

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

### cat3_minus_2 — RESOLVED (5'-rescue equivalence-extension; iterated 2026-05-18)

**Status:** shipped + iterated this session (commits 6d2cf59, then a follow-up
that fixes the slice geometry, adds the + strand trigger, and adds a k-sweep).
The rescue_3ss_truncation function detects when the body alignment overshoots
the annotated intron boundary by *k* bases AND the borrowed read bases at
those positions equal BOTH the original genome slice AND the would-be new
slice. The intron length changes accordingly and what would otherwise be a
trailing/leading kD op gets absorbed into the rescued M op.

**Correct geometric criteria (derived empirically — these are NOT the
ambiguity-window criteria; that's separate symmetric-slide work in the
design note above):**

  - strand (body overshoots intron_start):
    ```
    genome[ref_end - k : ref_end] == genome[intron_end : intron_end + k]
    AND read[last k of body M] equals both slices
    ```

  + strand (body overshoots intron_end):
    ```
    genome[ref_start : ref_start + k] == genome[intron_start - k : intron_start]
    AND read[first k of body M, after softclip] equals both slices
    ```

**k-sweep**: tries k = min(overshoot, 10) downward and accepts the largest
k where the equivalence holds. Allows partial absorption when full-overshoot k
doesn't qualify.

**Initial (commit 6d2cf59) bug**: the - strand `_ref_left` slice used
`genome[intron_start - k : intron_start]` instead of the correct
`genome[ref_end - k : ref_end]`. For cat3_minus_2 these slices happened to
both be `CT` (chrII at positions 366500-366503 has a CT-CT repeat pattern),
so the test passed by coincidence. The follow-up commit corrects the slice
to the geometrically correct position.

**Still deferred (next session work):**

- **cat3_plus_2 "off-by-1 acceptor" pattern**: this is an UNDERSHOOT, not
  overshoot — body M starts at `ref_start = intron_end - 1` (1 bp short of
  canonical). The current overshoot-only detection doesn't fire for this
  pattern. Fix would be the mirror: extend body M by k bases (taking from
  the END of the upstream softclip-aligned region) when ref_start <
  intron_end. Same equivalence concept but rotated. Untouched this session.
- **Symmetric-slide ambiguity window + motif-strength tiebreaker at the
  consensus-selection level**: see "Design note" at the top of this file.
  The existing rescue function already does ambiguity-window SEARCH for
  best junction placement (lines 1191-1218); what's missing is the
  motif-based tiebreaker in `merge_corrected_tsvs`.

### cat3_plus_2 — `M 2D N` vs clean `N` (+ strand follow-up)

Per-aligner picture (current bundle):
```
mapPacBio  N(142253,142619)  14=1D9= N366 50=1I2=  (canonical placement;
                                                    HP-ED 27.60 — WORST of 5)
others     N(142253,142618)  ...  (off-by-1 on acceptor, but cleaner HP-ED
                                   so they win consensus)
```

Different geometric pattern than cat3_minus_2 — here 4 of 5 aligners place
the intron's RIGHT boundary 1 bp earlier than canonical, rather than the
body alignment overshooting the donor on the - strand body side. The
- strand fix's body-overshoot detection won't fire because there isn't an
overshoot.

The + strand mirror of the - strand fix (trim N bases off the START of the
body M after the softclip) is plumbed through `extend_read_5prime_for_junction_rescue`
but `rescue_3ss_truncation` only emits `five_prime_upstream_trim > 0` for
- strand reads. For + strand cat3_plus_2's "off-by-1 acceptor" pattern, the
detection logic needs to be different — possibly checking whether the
annotated intron's acceptor is 1 bp later than the aligner's placed
acceptor, AND the equivalence criterion still holds. Deferred.

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
