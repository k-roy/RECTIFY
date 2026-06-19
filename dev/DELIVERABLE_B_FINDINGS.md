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

---

# SG-NEx Illumina short-read VALIDATION of the candidates (2026-06-18) — VERDICT-CHANGING

Job chain 30208125 (pull) → 30208350 (validate). Driver `dev/sgnex_illumina_validate.py`. A junction is
"validated" with ≥2 anchored split reads in ≥2 of the 3 A549 Illumina replicates (ambiguity-normalized).

| set | validated | rate | any Illumina read (≥1) |
| --- | --- | --- | --- |
| **the_111** (GMAP-only recurrent novels) | **0 / 111** | **0.0%** | **1 / 111** |
| corrob_609 (multi-aligner novels, POSITIVE control) | 87 / 609 | 14.3% | 163 / 609 (100 with ≥5 reads, 62 in 3/3 reps) |
| noncanon 3000 (NEGATIVE control) | 21 / 3000 | 0.7% | — |

## Verdict: the 111 GMAP-only recurrent novels are GMAP ARTIFACTS, NOT novel biology
The method works (positive control shows a clear graded real signal: 100 junctions with ≥5 reads, 62 in
3/3 reps; negative control at the 0.7% noise floor). The 111 are **NOT low-coverage/undetectable** — their
loci are heavily expressed and heavily spliced in Illumina:
- 5:171388191-171388472 (GMAP 109 reads): real Illumina junction **171388006-171390050 = 14,931 reads**,
  acceptor **+1,578 bp** away; GMAP's exact coord = 0.
- 5:181237624-181238158 (GMAP 61): real **181237719-181238098 = 8,594 reads**, donor +95 / acceptor −60 bp
  off (NOT an ambiguity-window shift).
- 5:179821705-179829510 (GMAP 52, ~7.8kb): real **179821141-179822957 = 1,721 reads**, acceptor **−6,553 bp**
  (GMAP fabricated a spurious long-range splice).
So GMAP places junctions tens-of-bp to kilobases from the real, heavily-supported splice sites. recurrence +
GT-AG-within-window was NOT evidence of truth — these were recurrent ONT-platform artifacts.

## Consequences
1. **GMAP's "fenced seat" is in serious doubt.** Its headline unique value — the ~111/127 recurrent novel
   junctions — does NOT survive orthogonal validation. Combined with 89.5% non-canonical noise, GMAP's
   UNIQUE contribution is artifact; any real value is only via the MULTI-aligner-corroborated set (where the
   87 short-read-validated junctions live), not its solo calls. → revisit the keep-GMAP decision.
2. **Vindicates the two-validity-claims rule** (advisor) and the PI's instinct to pull short reads: ONT
   recurrence (placement mechanics within one platform) ≠ error-model realism. The orthogonal gold standard
   caught artifacts the ONT-only analysis rated as promising.
3. **Design-doc impact:** the de-novo discovery track must require orthogonal/short-read corroboration, NOT
   recurrence+motif (which passed all 111 false positives). Strengthens the §8 abstain band + the C6
   variant-aware guard + LLR arbitration. "Recurrent + GT-AG" is necessary, not sufficient.
4. **The real novels are the 87** short-read-validated members of the 609 multi-aligner set (esp. the 62 in
   3/3 reps) — a high-confidence chr5 novel-junction set worth carrying forward (some may be genuinely
   novel vs annotation; intersect with gencode to confirm).

Caveat: spot-checked the top 3 of the 111 in coordinate detail; the class-level signal (110/111 zero
Illumina reads at exact coords) + these three make the artifact verdict robust. Single marginal case had
2-4 reads in 1 rep.

---

# CAVEATS + STAR-independent coverage check (2026-06-18, after PI challenge on the gold-standard assumption)

PI rightly challenged trusting the SG-NEx Illumina BAMs as a gold standard. Provenance + checks:
- **I did NOT align these reads — SG-NEx did.** `@PG`: **STAR 2.6.0c, SINGLE-PASS** (no `--twopassMode`),
  Ensembl-91 GRCh38 index, ~74M primary reads/rep. **STAR is well-documented as poor for NOVEL junctions
  (high FDR; 2-pass can make it worse), so the split-read 0/111 result is confounded by STAR's novel-junction
  insensitivity** and must NOT be treated as definitive on its own.
- Triage of the 111 (3 reps): **0 in low-coverage, 0 in long-intron (>590kb alignIntronMax)** — all 111 loci
  ARE expressed/covered by Illumina (often thousands of reads). So "undetectable due to low expression" is
  NOT the explanation. 20/111 have a competing junction within 50bp (clear artifacts); 91 "covered, no
  competing-within-50bp" (but the top-3 spot-checks had real junctions 185–6553bp away → the 50bp window
  undercounts artifacts).
- **STAR-INDEPENDENT coverage test (does GMAP's intron actually splice out?):** inside/flank coverage ratio
  for ALL 111 is **>0.5** (min 0.785; many >1.0, e.g. 181237516-181239425 ratio 1.46). A real intron is
  DEPLETED inside (spliced out); these are NOT depleted → the "intron" sequence is present in the
  transcriptome → consistent with GMAP fabricating splices in expressed/contiguous sequence. This test does
  NOT rely on STAR detecting junctions (intron-interior coverage = contiguous reads, not junction-spanning),
  so it survives the STAR critique. Residual alt-explanation: intron retention with a minor spliced form
  (the split-read test argues against a spliced subpopulation, but that test is STAR-confounded).

**STATUS: the 111 are LIKELY artifacts (coverage test is strong + STAR-independent), but the split-read
validation must be REDONE properly** with a multi-aligner short-read consensus (bbmap + minimap2 + a third),
COMPASS-integrated — NOT STAR alone. The positive-control 14.3% and the whole "validated 87" set ALSO need
re-deriving under the proper pipeline before being trusted. → see the COMPASS-short-read scoping effort.
