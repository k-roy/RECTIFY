# 644h — the phase-2 overhang-quality likelihood, measured: it upgrades the canonical track, rank-orders everything, and CANNOT alone purify the non-canonical track

**Date:** 2026-08-11. **Question (the PI's phase-2 metric, STATIONC_MAPPACBIO_HARVEST §PI-reframe):
does the longest high-quality SHORT-EXON-side overhang — complexity × length × low error —
separate the 37 residual gold from the junk that support-recurrence and motif could not (644f/644g)?**

**Method:** `scripts/benchmark/644h_phase2_overhang_likelihood.py` (H2 `644_accept/`; out
`t3/full/644h_phase2_overhang.json`). Per supporting read of each mapPacBio
beyond-(mm2 ∪ resolver_v5) junction: short-exon side = the side with the smaller contiguous
aligned anchor; overhang = ≤60 junction-adjacent query bases; I_eff = 641's
`effective_information_bits` verbatim (the resolver's currency); errors = mismatches + indel
bases in the window vs genome (the mapPacBio BAM's `=`-compressed SEQ decoded against the
reference — the smoke caught this: all-zero scores until decoded). **q_read = max(0, I_eff −
2·errors); junction q_max = max over reads.** Calibration positive class = annotated junctions
in the SAME arm (255 junctions, ≤10 reads each; q_max p25/p50 = 91/97 bits). Frames reproduce
644f/644g exactly (37 gold / 5,414 junk; phase-1 survivors 35/5,068; recurrent 17/1,737).

## Results

| frame / track | gold | junk | AUC (q_max, gold vs junk) |
|---|---:|---:|---:|
| all / canonical-in-class | 23 | 62 | 0.78 |
| all / non-canonical | 14 | 5,352 | 0.77 |
| phase1-surv≥2 / canonical | 12 | 12 | **0.91** |
| phase1-surv≥2 / non-canonical | 5 | 1,725 | 0.69 |

**Operating points that matter:**

- **Canonical track (surv≥2): q≥40 keeps 11/12 gold vs 3/12 junk** — the 644f support-only
  shortlist (~48% precision) upgrades to **~79% precision at 92% recall**. The overhang
  likelihood is ready to be the canonical track's adjudication axis right now.
- **Non-canonical track (surv≥1, the Gould class): q≥80 keeps 11/14 gold vs 1,357 junk**
  (rDNA-wide excluded; 3.7× junk cut at 79% Gould recall — a stronger cut than support≥2's 3×
  at 43% recall, and composable with it). But composing q≥80 AND surv≥2 yields 4/14 gold vs
  675 junk — **~0.6% precision. The metric alone cannot purify this track at any threshold.**
- The gold that scores low is honest: the bottom-5 gold (q 8–16) are single reads whose
  short-side anchor is 9–14 bp — genuinely indistinguishable evidence per-read (they are
  canonical-track junctions, admitted by track-1 support+class, not by phase 2).

## The mechanism finding: clean overhangs cannot see COPY misplacement

The top-scoring junk is the tell: rDNA-array junctions (chrXII 451–490 kb; 106–109 bits,
60 bp clean) and boundary-jitter twins of real junctions (chrVIII:498706-498790 /
498719-498790 = jittered variants of the #1 gold Gould junction chrVIII:498706-498786,
outside the ambiguity class). mapPacBio's systematic junk is misalignment between
NEAR-IDENTICAL sequence copies — so the overhang matches the (wrong) reference perfectly.
**Local overhang quality measures "is this placement locally well-supported"; it is
structurally blind to "is this placement unique".** rDNA explains only ~12% of the ≥80-bit
non-canonical junk (180/1,537 surv≥1; `RECTIFY_SKIP_REGIONS=yeast-rdna` removes it anyway) —
the rest is dispersed homology, Ty/paralog-class.

## Consequences for Station C (the spec, updated with numbers)

1. **Phase-2 scorer = overhang quality × placement UNIQUENESS × recurrence.** The missing leg
   is uniqueness/mappability of the overhang+junction context — the same T6 decoy check named
   for the resolver, now measured as the binding constraint on the non-canonical track. The
   644f "orthogonal evidence" requirement stands, sharpened: uniqueness is the orthogonal axis
   to build FIRST (cross-sample recurrence and short-read corroboration remain).
2. **Adopt the canonical-track upgrade now:** support≥2 + canonical-in-class + q≥40 → ~79%
   precision, 11/12 gold — a reviewable shortlist that no longer needs manual triage.
3. **The jitter-twin observation feeds the two-sided enumeration design**
   (`644h_realigner_integration.md` §3): near-boundary junk variants of a real junction are
   VOTES for the locus under one-sided/two-sided candidate enumeration, not independent
   junctions — locus-scoped aggregation will both purify junk and strengthen gold recurrence.
4. **The per-read machinery transfers as-is** to the triage clip legs: q_read is computed from
   the BAM + genome alone (no truth set), so the same scorer prices one-sided-anchor
   candidates and terminal-clip rescues at adjudication time.

**Files:** tool `scripts/benchmark/644h_phase2_overhang_likelihood.py` (this branch); data
`t3/full/644h_phase2_overhang.json` (H2); builds on 644f (`STATIONC_MAPPACBIO_HARVEST`), 644g
(`PHASE1_CONTEST_644G`), and the 641→Re-aligner integration handoff
(`~/work/UCLA/Chanfreau_Lab/planning/644h_realigner_integration.md`).
