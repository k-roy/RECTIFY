# Deliverable B — GMAP novel-junction read-level corroboration (A549 chr5) — 2026-06-18

Aligner-stratified, ambiguity-normalized, anchored (≥10bp) corroboration of every junction GMAP places
on A549 chr5, against the other 4 panel aligners. Driver `dev/deliverable_b_gmap_corroboration.py`
(commit 0326f93); Sherlock job 30205129 (COMPLETED, 2:22). Genome contigs registered (no chr5→chrV).

## Per-aligner anchored junctions (chr5, distinct after ambiguity-normalization)
| aligner | anchored junctions |
| --- | --- |
| minimap2 | 28,031 |
| uLTRA | 54,156 |
| deSALT | 75,124 |
| mapPacBio | 3,558 |
| **GMAP** | **221,354** |
annotated catalog (gencode v44, chr5): 18,450 junctions.

## GMAP's junctions, classified
| bucket | count | interpretation |
| --- | --- | --- |
| gmap_annotated | 2,785 | GMAP anchors a catalogued junction |
| **gmap_noncanonical** | **198,144** | non-GT-AG within ambiguity window → **the noise the fences suppress** |
| independently_corroborated (novel-canonical) | 609 | novel GT-AG junction GMAP shares with ≥1 other aligner |
| **gmap_only_recurrent (novel-canonical, ≥5 reads)** | **111** | **GMAP's genuine UNIQUE contribution** (was loosely "~127") |
| gmap_only_singleton (novel-canonical, <5 reads) | 2,300 | likely noise |

## Headline — the "fenced seat" verdict, now quantified
- **89.5% of GMAP's anchored junctions (198,144 / 221,354) are non-canonical** — even the 10bp anchor gate
  does not filter GMAP's noise; only the scoring fences (non-canonical −3 + support gate + ambiguity
  match) do. This is the strongest quantitative case yet that GMAP is net-positive ONLY behind the fences.
- **GMAP's real unique value = 111 recurrent canonical novel junctions** no other aligner anchors — buried
  in ~200k non-canonical + 2,300 singleton noise. The fences extract the 111 while suppressing the rest.
- Independently-corroborated novel-canonical junctions (609) by # corroborating aligners: {1:173, 2:259,
  3:157, **4:20**}. The 20 corroborated by ALL 4 others are the highest-confidence novel junctions on chr5.

## The 111 candidates → `dev/gmap_only_recurrent_novels_chr5.tsv`
Top by GMAP read depth: chr5:171388191-171388472 (109 reads), chr5:181237624-181238158 (61),
chr5:179821705-179829510 (52, ~7.8kb intron). NOTE clustering at chr5:171388xxx and chr5:181237xxx
(multiple novel junctions each) — either alt-splicing hotspots OR a GMAP systematic artifact at those
loci. **Flag for C6 variant-aware check:** are these loci variant/repeat-rich (variant-induced pseudo-
junctions)? This is exactly the discovery-FDR failure mode C6 guards.

## Validity caveat (the two-claims rule)
The 111 are recurrent + GT-AG + anchored — STRONGER evidence, NOT proof. By construction they have ZERO
independent aligner corroboration, so their realness rests on recurrence + canonical motif only. Orthogonal
validation still owed: the P0 simulation benchmark, A549 short-read RNA-seq junction support, or
curated splice DBs. recurrence + GT-AG ≠ proof.

Raw 39 MB JSON stays on Sherlock scratch (`/scratch/users/kevinroy/deliverable_b/gmap_corroboration.json`,
purgeable; regenerable from the committed driver). Only this doc + the 111-row TSV are kept durably.
