# 644f — Station-C harvest census: how much of mapPacBio can pool-level arbitration keep?

**Date:** 2026-08-10. **Question (Kevin):** of the gold junctions the resolver misses, how many
would Station C flag from a mapPacBio alignment — and how many junk junctions come with them?
Is mapPacBio's discovery potential recoverable under judicious arbitration?

**Method:** `644f_stationc_gate_census.py` (H2 `644_accept/`; output
`t3/full/644f_stationc_gate.json`). Same census + ambiguity-canonicalised matching as the 644
T3 score (primary reads, min anchor 8); per-junction read support summed over exact variants of
each canonical junction. Frames: mapPacBio-beyond-mm2 (the +105/5,456 set) and
**mapPacBio-beyond-(mm2 ∪ resolver_v5)** — the true Station-C residual. Sweep the support gate;
split junk and gold by canonical-in-class (GT/GC..AG within the ambiguity window).

## Results

**True residual (mapPacBio-only, absent from mm2 AND resolver v5.1): 37 gold / 5,414 junk.**
(The planning/644b "22-junction residual" was category arithmetic; the set-difference is 37 —
the resolver also holds ~15 gold junctions mapPacBio lacks, mostly protointrons.)

| support ≥ | gold kept (of 37) | junk kept (of 5,414) | junk canonical-class |
|---:|---:|---:|---:|
| 1 | 37 | 5,414 | 62 |
| 2 | 18 | 1,792 | 13 |
| 3 | 14 | 1,248 | 6 |
| 5 | 11 | 824 | 3 |

**The class split is the finding.** The 37 residual gold decompose into:
- **23 canonical-in-class** (20 alt introns + 1 protointron + 2 standard) — 12 survive ≥2;
- **14 Gould-S6 novels — 0/14 canonical-in-class** (the entire class the splice-site index
  cannot see is non-canonical), 6 survive ≥2.

### Two-track gate performance (the only honest framing)

- **Canonical track** (support ≥2 AND canonical-in-class): admits **12 gold + 13 junk = 25
  candidates genome-wide** → ~48% precision, up from 0.68% raw — a ~70× purification that
  harvests half the canonical residual as a trivially reviewable shortlist. (And "junk" =
  not-in-catalogue; per the GCR1-alt lesson some of the 13 may be real uncatalogued junctions,
  so 48% is a lower bound.)
- **Non-canonical track** (support ≥2, not canonical-class): 6 Gould gold vs **1,779 junk —
  ~0.3% precision. Support recurrence CANNOT purify this track**, at any threshold
  (≥5: 2 gold / 821 junk).

### ⚠ The scout-mode statistics correction (honest, load-bearing)

**mapPacBio's junk is NOT per-read-random: 33% of its junk junctions carry ≥2 reads**
(1,792/5,414). The "fabrication doesn't recur, genuine signal does" argument (landscape doc
§3b scout mode) holds for random fabrication but NOT for mapPacBio's systematic
misalignment (repeat/homology-driven, reproduced read after read). Recurrence remains a
useful FIRST gate (3× reduction) — it is nowhere near a sufficient one for the
non-canonical track.

## Answer to the PI's question

1. **Flagged residual gold:** a bare support≥2 gate flags **18/37 (49%)** of the
   mapPacBio-only gold. With the class split: **12 canonical residuals are immediately
   harvestable at ~48% precision** (25-candidate shortlist → hp_ed/COMPASS adjudication);
   **6 recurrent Gould novels are flagged but arrive buried in 1,779 recurrent junk.**
2. **Junk flagged alongside:** 25 candidates on the canonical track (manageable); 1,779 on the
   non-canonical track (not manageable by support alone).
3. **The verdict on "enormous discovery potential, judiciously handled": half right, and the
   half matters.** The canonical residual — the alt-intron edge mapPacBio still holds over the
   resolver — is fully harvestable by a cheap pool gate + the existing arbiter. The genuinely
   unique class (Gould novels: non-canonical, index-invisible, reachable by NO other arm) is
   real discovery potential but needs **orthogonal evidence** to admit: short-read/COMPASS
   corroboration, cross-sample recurrence, co-localized mm2-side distress (M-D-N signatures,
   junction-proximal error density), mappability. This is precisely the pre-registered
   two-FDR-track discipline (canonical and non-canonical tracks never share a threshold) from
   the member-design era — now with measured numbers attached.

**Station-C spec consequence:** mapPacBio probe admission = two tracks. Track 1 (canonical):
support ≥2 + canonical-in-class → straight to hp_ed re-entry/adjudication. Track 2
(non-canonical): support ≥2 is a PRE-filter only; admission requires ≥1 orthogonal evidence
class; abstain (report, don't assert) otherwise. Next measurement worth its cost: how many of
the 1,779 track-2 junk (and the 6 Gould) show co-localized minimap2-side distress — that is the
scout-free trigger's true precision, and the 644d D-signature census machinery already exists.
