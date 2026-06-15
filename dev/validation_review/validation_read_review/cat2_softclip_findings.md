# cat2_softclip — findings surfaced by the End-correction panel

Written by the **plotter session** for the **debugger session** to consume.
Plotter does not own walkback / indel_corrector / scoring code; these are
observations only, with the user's verbatim diagnosis preserved.

**Source bundle**: `rectify/data/validation/rectified/` regenerated 2026-05-18.

---

## cat2_plus_1 (61b0c014, chrI:23704–23804, + strand) — HP-ED ranks parsimony-wrong

**Panel summary (verbatim from rendered PNG):**

```
aligner       5' orig  5' corr  Δ    3' orig  3' corr  Δ    HP-ED  span  chim  pick
minimap2      23362    23362    0    23759    23759    0    20.1   398   ✓
gapmm2        23362    23362    0    23759    23759    0    20.1   398   ✓
mapPacBio     23362    23362    0    23754    23754    0    16.9   393   ✓
deSALT        23362    23362    0    23754    23754    0    15.6   393   ✓   ← winner
uLTRA         23362    23362    0    23759    23759    0    20.1   398   ✓
```

**Two-cluster divergence**: minimap2/gapmm2/uLTRA agree on `3'=23759` (5 bp
downstream of deSALT/mapPacBio's `3'=23754`). HP-ED ranks deSALT lowest →
deSALT wins the bundle. The user argues this ranking is biologically wrong.

**Genome at the divergence (from PNG ref row):**
```
…23722         23752 23759            
A-tract       AATTAAATAAATAAATAAATA CAAT GATAT…
                ↑ deSALT stops here    ↑ minimap2 stops here
                  (23754)                (23759)
```

The 5-bp gap is the `AAATA` between deSALT's 3'-end and minimap2's 3'-end.

**User's diagnosis (verbatim, 2026-05-18)**:

> "For cat2 plus1, I actually think minimap2 did the right thing. The
> undercalling of the A-tract is the most plausible alignment here, and it
> allows everything else to cleanly align. The extra AAAT repeats are
> classic nanopore over-calling of tetramer repeats. The hp-aware edit
> distance should be allowing minimap2 to win this one outright."

**Reading the rendered tracks (PNG):**

- **minimap2 row**: 8 deletions early in the A-tract (orange `------` block).
  Body alignment extends cleanly to 23759, then a short softclip tail
  (`AATAAAAT`) past the body.
- **corrected / pA-rest tracks** (built from deSALT-as-winner): show three
  separate 1I insertion pills (purple ↓ markers labeled `+1` with base
  `T`) sitting in the `AAATAAATAAATAAAT` AAAT-repeat region. Body alignment
  stops at 23754.

The pattern matches the user's parsimony argument:
- minimap2's representation: **1× 8-bp HP deletion** in a long A-tract
  (nanopore HP under-call) — biologically expected, should be HP-cheap.
- deSALT's representation: **3× single-base insertions** in AAAT tetramer
  repeats (nanopore tetramer over-call) — also expected, but 3 independent
  events. Plus a body alignment that stops 5 bp short, leaving more bases
  as soft/hard clip.

**Hypothesis for the debugger**:

If `del_cost(hp_len ≈ 19, base='A')` were sufficiently cheap (per the
empirical penalty table at
`rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv`),
minimap2's HP-ED would dip below deSALT's. Two possibilities:

1. **Penalty table miscalibration**: `del_cost` for long A-runs isn't as
   cheap as it should be empirically. Worth checking the row for `del`,
   base `A`, hp_len 15–20 in the TSV. Empirically nanopore DRS undercalls
   long A-runs at high frequency; this should be very cheap.
2. **Indep. vs. correlated penalty**: deSALT's 3 separate insertions
   probably get scored at `3 × ins_cost(hp_len=4, base='T')` even though
   they're correlated errors at AAAT tetramer boundaries. If the penalty
   model assumes independence, it under-counts the joint cost. Out of
   scope for the immediate fix.
3. **Soft/hard-clip penalty mismatch**: deSALT covers 5 fewer aligned bases
   (393 vs 398). Those 5 bp are soft/hard-clipped at 1.0/base = +5.0 ED.
   So deSALT's "win" margin of 4.5 ED (20.1 − 15.6) actually comes from
   ~3 cheap insertions vs minimap2's 8 expensive deletions — the gap is
   _entirely_ in the per-event ED ratio.

**Suggested debugger workflow**:

1. Dump the per-op ED contribution for cat2_plus_1's minimap2-corrected vs
   deSALT-corrected BAMs (use the diagnostic from the earlier `=`-SEQ
   investigation). Confirm where the 4.5-ED gap is concentrated.
2. If it's in the long-A-run `del_cost`, revisit the penalty table calibration.
3. Consider whether HP-ED should normalize per "biological event" rather
   than per "CIGAR-op base" — a single 8-bp HP under-call is one event,
   not eight.

---

## cat2_minus_1 — clean

Reviewed by user 2026-05-18: looks correct as rendered. No debugger action.

---

## cat2_minus_2 (b313b50d, chrI:128052–128152, − strand, YAL014C) — body extension blocked by single T-del; soft-clip retains 3 matchable nt

**Panel summary**:

```
aligner       5' orig   5' corr   Δ   3' orig   3' corr   Δ    HP-ED  span  pick
minimap2      129051    129051    0   128112    128102   +10   86.9   940
gapmm2        129051    129051    0   128112    128102   +10   86.9   940
mapPacBio     129051    129051    0   128106    128106    0    88.1   946
deSALT        129051    129051    0   128112    128102   +10   86.9   940   ← winner
uLTRA         129051    129051    0   128112    128102   +10   86.9   940
```

Four of five aligners agree the corrected 3' lands at 128102 (the T position
on the RNA strand). mapPacBio uniquely stops 4 bp short at 128106 with Δ=0.

**What the rendered tracks show** (corrected + pA-rest):

- Alignment body ends at the `T` (RNA-direction terminal non-A) at 128102.
- The corrected track then has a long `-----------` deletion (orange) where
  the read undercalls a long genomic A-tract.
- The pA-rest track shows the trailing `GCAATAAA` as soft-clipped (faded
  grey) — these are unmatched bases past the body end.

**User's diagnosis (verbatim, 2026-05-18)**:

> "cat2 minus 2 is tricky, but I suspect the right call would be to allow
> the TGC to bind and accept a T deletion as part of the massive A
> homopolymer undercall for this read. Allowing for 3 straight nt to match
> should outweigh the single T deletion here."

**Interpretation**: the trailing `GCAAT` soft-clip contains a `TGC` (RNA
direction, read-3'→5') that maps to an upstream `TGC` in the genome at
~128097–128099. To incorporate those 3 matches into the body alignment,
the algorithm would have to accept ONE additional `T` deletion as part of
the (already large) A-homopolymer undercall the read is exhibiting. The
parsimony argument: **3 new matches − 1 new del = +2 informational gain**,
should always be a net win for any sensible scoring.

**Hypothesis for the debugger**:

Cat2 soft-clip rescue (Module 2G) is the obvious candidate. The current
behavior: the rescue extended into the A-tract but stopped at the first T
(128102) and left the `GCAAT` past it as soft-clip. Specifically:
- The current rescue logic appears to terminate when it encounters a
  non-A base. It should instead consider: "is there a downstream non-A
  read=ref match that's reachable through a cheap HP del?"
- The lookahead would be: starting from the current rescue endpoint, scan
  the next ~10 bp of soft-clip for a 2-3-bp read=ref match, and if found,
  emit a HP del + the matched bases.

Implementation sketch (debugger's call to make):
1. In `extend_read_3prime_for_softclip_rescue` (or a successor pass), after
   the current rescue completes, take the soft-clip suffix.
2. Try to align the first ~10 bp of the soft-clip against the next ~10–20
   bp of upstream genome (allowing a single HP del slip).
3. Accept the extension if it produces ≥3 consecutive matches with ≤1
   cheap HP del. Reject if the cost exceeds 1.5× the gain (suppress
   speculative extensions).

mapPacBio's behavior here (Δ=0, body at 128106) is the conservative end:
it never extended the rescue at all. That's the wrong direction for the
test scenario this read encodes (Cat2 = soft-clip rescue should fire).

---

## cat2_plus_2 (88953e9c, chrVI:8564–8664, + strand) — clean correction; surfaces an "effective-utility" feature request

**Panel summary**:

```
aligner       5' orig  5' corr  Δ   3' orig  3' corr  Δ    HP-ED  span  pick
minimap2      7832     7832     0   8601     8605    +4    21.3   770
gapmm2        7832     7832     0   8601     8605    +4    21.3   770
mapPacBio     7832     7832     0   8614     8614     0    20.4   783
deSALT        7832     7832     0   8601     8605    +4    21.3   770
uLTRA         7832     7832     0   8614     8614     0    19.2   783   ← winner
```

Correction itself is fine (user: "looks great"). The interesting finding
is a **two-cluster pattern** in the per-aligner 3' calls:

- **Cluster A** (3'=8614, Δ=0): mapPacBio, uLTRA. The body alignment
  already ended at the correct non-A in the aligner output; no walkback
  needed.
- **Cluster B** (3'=8605, Δ=+4): minimap2, gapmm2, deSALT. These three
  needed a +4 bp walkback to arrive at 8605, a *different* corrected
  position than cluster A's 8614.

uLTRA wins HP-ED within cluster A (19.2 vs 20.4 — both more biologically
correct than cluster B, presumably). The user observes: "uLTRA did a
cleaner job of aligning the mRNA body that was more favorable with our
empirical penalty table."

---

## Update — 2026-05-18, formal cat1–3 review

The freshly regenerated bundle (commit `5df89c9` + walkback fixes) still
shows the cat2_minus_2 issue. The user has restated the request as a
**generalized "long-deletion-extension" policy**:

**User's verbatim ask (2026-05-18)**:

> "Let's also add a note to handle cases like cat2 minus2, where simply
> extending the deletion by 2 bases allows for a full TGC match as
> opposed to just the T. Makes much more sense to me to maximize the
> ending matches when an already long deletion is simply extended."

### Proposed policy

When a CIGAR contains an **already-long deletion** (≥ ~5 bp) and the
adjacent soft-clip on the same side has a short matching run with the
reference at positions reachable by extending the deletion:

1. Quantify the trade — for each candidate extension length `k`:
   - Cost: `k × del_cost(hp_context, base)` extra deletion bases.
   - Gain: `n_matches × 1.0` (soft-clip `S` reclassified as match `=`).
2. Accept the extension if `gain > cost` (i.e. each additional deletion
   base pays for itself by enabling at least ~2 matches in HP context).
3. Cap the extension at a reasonable max (e.g. `k ≤ 5`) to prevent
   pathological deletion runs.

### Worked example (cat2_minus_2)

```
current: …CCCAAAAA --------- T GCAATAAA(softclip)
                  9-bp del   T = 1 match, then 8-bp soft-clip
                              ED contrib: 9 × del_cost(HP-A,A) + 8 × 1.0
                              ≈ 9 × 0.43 + 8.0 = 11.87 ED

proposed: …CCCAAAAA ----------- TGC AATAAA(softclip)
                   11-bp del    3 matches,    6-bp soft-clip
                                ED contrib: 11 × del_cost + 0 + 6 × 1.0
                                ≈ 11 × 0.43 + 6.0 = 10.73 ED
```

Net improvement: **−1.14 ED** (assuming flat HP-A `del_cost` ≈ 0.43).
Also biologically more parsimonious — the read truly does match `TGC`
at those positions; the long A-tract deletion bridges to it.

### Relation to the cat3 mpb 5'-anchor finding (task #12)

Both findings share the underlying principle: **when the existing
alignment has a noisy edge (long deletion OR mismatch/indel cluster),
extending the noise slightly to reach a longer clean-match region is
net informational gain**. The two cases differ only in geometry:

- cat2_minus_2 (this finding): long deletion in BODY interior →
  extend to absorb adjacent soft-clip matches into body.
- cat3 mpb (task #12): mismatch/indel cluster at 5' EDGE → reanchor
  past the cluster to enable 5'-rescue.

Both could be unified into a single "extend-noisy-region-toward-clean-
match-anchor" pass in the correction pipeline, parameterized by side
(5' vs 3' / body-interior vs body-edge).

---

## Feature request: per-aligner "effective utility" tracking

**User's verbatim ask (2026-05-18)**:

> "I'm wondering if we should have a separate column that acknowledges
> which aligners effectively got the same 5′, 3′ and junctions (if
> applicable). This is useful for knowing the 'effective' utility of each
> aligner for rectify. Presumably, uLTRA did a cleaner job of aligning the
> mRNA body that was more favorable with our empirical penalty table. Both
> analyses are useful. The 'effective' utility should be something that is
> tracked at read-level and sample-wide as well."

**Read-level (plotter scope)**: add an "eff. group" column to the
End-correction panel between `chim` and `pick`. Cluster the per-aligner
rows by `(corrected_5p, corrected_3p, frozenset(corrected_junctions))`,
assign each distinct tuple a letter (A, B, C, …), and label each row.
Aligners in the same cluster ARE biologically interchangeable for RECTIFY,
even if their CIGARs differ in mismatch/indel placement (which is what
HP-ED then ranks within the cluster).

**Sample-wide (debugger scope)**: a per-aligner aggregate of "how often
this aligner's 5'/3'/junctions matched the winning cluster". Suggested
output:

```
sample-wide effective utility (out of N reads):
  minimap2:   matched winner in 7,234 / 9,800 reads (73.8%)
  gapmm2:     matched winner in 6,891 / 9,800 reads (70.3%)
  mapPacBio:  matched winner in 8,123 / 9,800 reads (82.9%)
  deSALT:     matched winner in 7,512 / 9,800 reads (76.7%)
  uLTRA:      matched winner in 8,401 / 9,800 reads (85.7%)
```

Could live as a new column in the per_aligner_summary.tsv (per-read
boolean `effectively_matched_winner`) plus a `--summary-stats` flag on
`merge_corrected_tsvs` that prints the rollup.

This is genuinely interesting information that the current
"winning_aligner" column hides: an aligner can lose HP-ED inside the
correct cluster (cluster A) and still be carrying the right biological
answer. That's a different signal than an aligner that ends up in the
WRONG cluster.

**Plotter could ship the read-level column now if you want it before the
sample-wide piece lands.** Just say the word.

---

## ✅ RESOLUTION — cat2_plus_1 (2026-06-01)

**The minimap2 "clip" (EER 56.6, corrected 3′=23711) is a STALE-BUNDLE ARTIFACT, not a real
behavior, and needs no walkback/trim/softclip-rescue/cigar-surgery change.** Full root-cause:

- `polya_walkback` is **not** involved — `walkback_drs_full` returns `None` for this read.
- `trim-polya` is correct — removed only the 12-base tail (`trimmed_3prime_seq`=12 A's).
- **Every current minimap2 aligns this read THROUGH with a clean single 9-bp deletion**
  (`…49M 9D 39M 8S`, 3′=23759): M1 binary 2.28, M1 mappy 2.30, **H2 & Sherlock production
  2.30-r1287**, and a fresh M1 conda 2.30 build. Confirmed re-align: minimap2 & uLTRA → 9D through.
- Controlled for version (changelog 2.28→2.30 shows no relevant change), arch, conda build, flags,
  junc-bed (the bundle's was EMPTY), and read seq — the 47H clip reproduced in **none**. It was a
  one-off artifact of the original bundle build (only unrecoverable variable: the since-deleted exact
  input fastq). Production aligns these A-tract reads correctly.

**Figure-refresh status (deferred):** regenerating the committed bundle to show the through-alignment
is a maintenance task blocked by (a) M1 gapmm2 v25.8.12 dropping 30/36 reads, (b) H2 `rectify align`
env/aligner-path setup, and (c) `tests/test_validation_reads.py` assertions on this read's 3′ end
(23711→~23758) needing a deliberate update. Do it as a focused pass once gapmm2 is fixed, not ad hoc.
Full detail: memory `project_cat2_plus1_minimap2_clip`.
