# 644g — the phase-1 per-read contest: both comfortable hypotheses fail; the burden lands on junction-aware scoring and phase 2

**Date:** 2026-08-10 (late night). **Context:** the PI's two-phase reframe — phase 1 (per-read
consensus) should refute mapPacBio's junk read-by-read; a 5′-clip routing scheme could buy its
unique gold at a fraction of the compute. **Both quantified here; both fail in their naive
form — while gold recall through consensus is confirmed.**

**Method:** `644g_phase1_contest.py` (H2 `644_accept/`; out `t3/full/644g_phase1_contest.json`).
For every read supporting a mapPacBio beyond-(mm2 ∪ resolver_v5) junction, contest mapPacBio's
vs minimap2's alignment of the SAME read under a flat score `NM + soft/hard-clip bases` (N-ops
free; both BAMs carry NM; lower wins). ⚠ This is a RAW-alignment flat proxy — production phase 1
contests CORRECTED alignments under `score_segment`/hp_ed with junction-aware fences. That
difference turns out to be the finding.

## Results

**Junk (5,414 junctions · 50,190 supporting read slots):**

| verdict | slots | share |
|---|---:|---:|
| mapPacBio wins | 42,659 | **85.0%** |
| tie | 4,491 | 8.9% |
| minimap2 wins | 2,971 | **5.9%** |
| mm2 absent (mpb wins by default) | 69 | 0.1% |

Junk junctions retaining ≥1 win-or-tie read: **5,068/5,414 (94%)**; ≥2: 1,737.

**Residual gold (37 junctions · 217 read slots):** mapPacBio wins 199 (92%) + 9 ties → **96%
win-or-tie**; minimap2 wins only 9. **35/37 gold junctions survive with ≥1 winning read**
(top alt-introns survive with 18–34 winning reads each).

**Routing (the 5′-clip idea):** of the 217 residual-gold read slots, only **20 (9%) are
5′-clipped ≥8 bp in minimap2's alignment** (6 at ≥15 bp). Budget side: 5′-clip ≥8 routing would
send 650,897/4,571,871 reads (14%). ⇒ routing 14% of compute **misses 91% of the reads carrying
mapPacBio's unique gold**. Any-clip ≥8 sends 18% of reads; the gold reads are simply not
clipped in minimap2 — it aligns them *confidently and differently* (the SRC1 mechanism again).

## The three conclusions

1. **"Most junk is refuted by another aligner scoring better on that read" is REFUTED under a
   naive alignment score — and the failure mode is exactly the free-N-op gaming the PI named.**
   A fabricated splice converts clip/mismatch cost into a free N: mapPacBio's junk alignments
   BEAT minimap2's honest clip on the same read 85% of the time. A flat alignment-quality
   contest cannot tell mapPacBio's real splices from its fabrications (junk wins 94%
   win-or-tie; gold wins 96% — **no discrimination**). This is the C3/GMAP lesson re-measured:
   the production `score_segment` fences (junction-proximal error charges, unsupported-junction
   terms, the canonical tiebreak) are not decoration — they are the ONLY thing standing between
   phase 1 and 5,068 fabricated junctions. The definitive phase-1 number (does the PRODUCTION
   corrected+fenced contest refute the junk?) needs a `correct`+consensus run over both arms on
   a gold-window subset — the natural next measurement, and a direct test of whether the fences
   built for GMAP's 878k-junction noise transfer to mapPacBio's.
2. **Gold recall through consensus is CONFIRMED.** 96% of residual-gold read slots win the
   contest; 35/37 junctions survive. With mapPacBio in the panel, phase 1 does NOT discard its
   true positives — the "pull it back into the panel" alarm does not fire on recall grounds.
   The cost of that recall is the 5,068-junction junk caseload it admits alongside — which is
   phase 2's job, exactly as the PI's framing anticipated ("the all-important phase 2").
3. **5′-clip routing is REFUTED as the economy measure** (9% gold-read recall at 14% compute).
   The reads mapPacBio uniquely rescues are unclipped, confidently-misplaced reads in minimap2 —
   invisible to any clip trigger. Viable alternatives, in cost order: (a) route on
   junction-proximal distress (M-D-N signatures + proximity-error density — the 644d census
   already finds the alt-intron class this way, from minimap2's own records, at zero mapPacBio
   compute); (b) route on pool-level disagreement loci (Station C scout); (c) full-depth
   mapPacBio in the panel, accepted as ~80 CPU-h, with phase 2 carrying the FDR. Given (a)
   recovered ~half of mapPacBio's contribution with NO mapPacBio at all (644d), the economical
   scheme is distress-routing, not clip-routing.

## Consequence for the Re-aligner program

Phase 2 (per-junction likelihood; the PI's metric: longest high-quality short-exon-side
overhang — complexity × length × low error) is not an enhancement, it is THE decision layer:
phase 1 either admits fabrications (naive score) or leans on fences whose transfer to
mapPacBio is unmeasured. Station C = phase 2, and its caseload is now sized: ~1,700–5,100
surviving junk junctions vs 35 recoverable gold under the naive bound; the fenced production
number will be smaller but not zero. Next measurements: (i) production-fenced phase-1 contest
on a gold-window subset (correct + consensus over both arms); (ii) the overhang-quality
likelihood prototype scored against this same census (does the PI's metric separate the 35
gold from the 1,737 recurrent junk where support and motif could not?).
