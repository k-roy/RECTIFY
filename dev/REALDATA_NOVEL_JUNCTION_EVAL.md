# Real-data per-aligner novel-junction precision — SG-NEx A549 chr5 ONT-DRS (2026-07-02)

The real-data adjudication (#2) of the mapPacBio question the SYNTHETIC gate was too benign to settle.
Agent a7427e3abf13d798b; results /scratch/users/kevinroy/nativealigner_realeval/jxn_eval_results.json.
PRELIMINARY (a corrected re-run 32528966 is pending — re-verify annotation-frame + the real-vs-artifact split).

## Support>=2 profile (real junctions; singleton noise filtered)
| aligner | distinct jxns | non-canon% | unique% | unique-non-canon n |
| --- | --- | --- | --- | --- |
| minimap2 | 8205 | 1.6 | 0.9 | 23 |
| uLTRA | 22614 | 2.2 | 0.9 | 45 |
| winnowmap2 | 19495 | 3.8 | 0.9 | 47 |
| deSALT | 25222 | 14.0 | 3.1 | 214 |
| minisplice_mm2 | 40257 | 13.8 | 14.6 | 1169 |
| mapPacBio | 28920 | 94.5 | 47.9 | 2519 |
| GMAP | 543378 | 88.8 | 66.3 | 59775 |
(mapPacBio only 165,846 primary reads vs ~318k others — it DROPS half the reads, real-data fragility.)

## VERDICT — the native member's TARGET is REAL on real data (reverses the inconclusive synthetic rgfdr)
1. WORKHORSES FLATTEN: minimap2/uLTRA/winnowmap2 are ~97-98% canonical even at support>=2 — they emit almost no
   novel non-canonical. deSALT is the least-biased precise member (14% non-canon, 3.1% unique).
2. mapPacBio + GMAP FABRICATE: 89-95% non-canonical, dominated by aligner-UNIQUE junctions (GMAP 59,775 unique
   non-canon), even at support>=2. The real-human pathology the synthetic gate MISSED — confirmed. (PI was right.)
3. THE GAP IS REAL: unannotated non-canonical junctions with real-looking support (e.g. 2590 reads) that ONLY the
   fabricators (GMAP+mapPacBio) find and the PRECISE members MISS. No member is BOTH sensitive to non-canonical
   AND precise => the calibrated-empirical native member's exact target.
CAVEATS (discipline): (a) corrected run pending; (b) recall-gap candidates found ONLY by GMAP+mapPacBio may be
SHARED ARTIFACTS of two permissive aligners — needs COMPASS short-read cross-platform validation to split
real-novel from artifact; (c) this confirms the TARGET is real; it does NOT prove the CURE (a calibrated member
recovers these AT PRECISION — the addressability ablation, still unproven, remains the make-or-break before/while building).
NEXT: recover the corrected run + the 3 COMPASS-anchor per-aligner placements (SLC35A4/TMED9/SQSTM1); cross-
validate the recall-gap candidates vs COMPASS short reads; THEN the addressability ablation defines the member spec.
