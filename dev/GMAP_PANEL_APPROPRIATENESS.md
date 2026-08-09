# Is GMAP appropriate for the RECTIFY long-read panel? — evidence synthesis (2026-06-29)

Track goal #2 ("validate GMAP's appropriateness for the panel"). Synthesizes three independent lines of
evidence: (A) the per-bucket GMAP junction census on A549 chr5 (deliverable_b), (B) the deep
multi-platform + WGS validation of GMAP's UNIQUE novel calls (the 111, this session), (C) the
aligner-version experiment + the live `test_gmap_fence_regression` rationale.

## A. GMAP junction census — A549 chr5 long-read (deliverable_b `gmap_corroboration.json`)

| bucket | count | % of GMAP chr5 calls |
| --- | ---: | ---: |
| **non-canonical** (not GT-AG/GC-AG within ambiguity window) | **198,144** | **~97.2 %** |
| annotated (canonical, in catalog) | 2,785 | 1.4 % |
| independently-corroborated novels (≥1 OTHER aligner agrees) | 609 | 0.3 % |
| gmap-only recurrent GT-AG novel (the "111") | 111 | 0.05 % |
| gmap-only singleton (canonical, <5 reads, no corroboration) | 2,300 | 1.1 % |
| **total GMAP calls** | **~203,949** | 100 % |

**~97 % of GMAP's raw junction output is non-canonical** — overwhelmingly spurious (the
`test_gmap_fence_regression` doc independently notes GMAP emits ~878k spurious non-canonical junctions =
55.6 % of its calls genome-wide; 198k on chr5 alone here). GMAP is a HIGH-SENSITIVITY, HIGH-NOISE member.

## B. GMAP's UNIQUE contribution is tiny and mostly artifact (this session's deep validation)

The 111 gmap-only recurrent GT-AG novels are GMAP's BEST unique-novel candidates (recurrent ≥5 reads,
canonical, no other aligner). Deep validation (ambiguity-tolerant short-read COMPASS panel + cross-aligner
long-read placement + A549 WGS deletion/rearrangement/CN tests + circRNA/fusion checks):

- **~107 / 111 = artifacts** (no independent short-read support even under ambiguity-tolerant matching).
- **3 / 111 = real intragenic novel RNA junctions** (SQSTM1/TMED9/SLC35A4; DNA-confirmed genomic, not
  deletion/rearrangement/circRNA/fusion) — but their dominant coord is NON-canonical with a canonical motif
  1-4 bp away (genuinely-non-canonical vs canonical-misplaced is open → C1 realign / RT-PCR).
- **1 / 111 = inconclusive** alt-acceptor (J2).
- 2,300 gmap-only singletons are below the recurrence gate → treated as noise.

So GMAP's GENUINELY-UNIQUE real-junction yield on chr5 is **~3-4 junctions**, against ~107 unique artifacts
(of the recurrent set) + 2,300 singletons + 198,144 non-canonical. **GMAP's panel value is NOT as a
unique-novel finder.** Its real value is as a CORROBORATOR (609 novels it shares with independent aligners)
and annotated-recoverer (2,785) — contributions that, by definition of "independently corroborated", the
OTHER panel members also make.

## C. Version + fence context

- **Version is negligible** (prior experiment): STAR/HISAT2/magicblast/bbmap latest-vs-pinned give
  novel-Jaccard 1.0 @≥2 reads; GMAP 2024-11-20 is NOT a drop-in (removed `--ambig-splice-noclip`, needs index
  rebuild, ~100x slower, and SIGSEGVs on this A549 data) → pinning gmap 2021-05-27 (as COMPASS does) is correct.
- **The fence is load-bearing:** `test_gmap_fence_regression` documents that GMAP earns its seat ONLY because
  the integer `score_segment` scoring SUPPRESSES its spurious non-canonical mass while keeping its genuine
  novels + annotated parity. Remove the fences → GMAP floods consensus with false non-canonical junctions.
  Fence-test figure (A549 **chr5**, `test_gmap_fence_regression.py:5-8`): GMAP = 55.6 % non-canonical
  (878k of 1.58M) vs 1.1 % for minimap2; "~100 genuine unique novel junctions."

- **RECONCILIATION (two count systems — do NOT subtract across them).** The fence test's chr5 totals
  (1.58M junctions, 878k non-canonical) are ~8x my census's (204k, 198k) although BOTH are "A549 chr5" —
  because they COUNT DIFFERENTLY (fence test ≈ read-supported junction OBSERVATIONS / all raw N-ops; the
  deliverable_b census ≈ UNIQUE normalized junctions after recurrence/anchor gates). The absolute numbers are
  therefore NOT subtractable. What IS consistent across both: GMAP is dominated by non-canonical calls
  (55.6 %–97 % depending on counting) — the qualitative noise picture is robust.

- **The "~100 genuine novels" claim needs a GENOME-WIDE validated count before it drives any decision.** The
  fence test's "~100 genuine unique novel junctions" is chr5 and most plausibly ≈ the 111 gmap-only-recurrent
  CANDIDATES — i.e. an UNvalidated candidate count, which this session reduces to **~3-4 actually real on
  chr5**. BUT (caution, per review): (a) it is NOT confirmed the fence-test "~100" is exactly the 111 set vs
  a broader definition; (b) chr5 is only ~6 % of the genome, so **~3-4 real on chr5 scales to ~50-65 real
  unique novels GENOME-WIDE** if chr5 is representative — a NON-negligible contribution. So GMAP's unique
  real-novel yield is SMALL PER-CHROMOSOME but plausibly tens genome-wide; "marginal" is not established
  without a genome-wide validation. The fence test's behavioral assertions stand regardless; only its
  motivating "~100 genuine" comment should be annotated "candidate count; chr5; ~3-4 validated-real (v3)".

## D. 5-CHROMOSOME extension (2026-06-29, `gmap_validate_harness.py` on chr5+chr1/11/17/19)

Ran the validation harness on 5 chromosomes (chr5 from the full panel; chr1/11/17/19 by adding GMAP-2021 to
the DRS agent's `morechrom` minimap2/deSALT/uLTRA panel, read-only). GMAP aligned 999,228 primary reads.

| chrom | gmap-only-recurrent-canonical-novel candidates | SUPPORTED (indep. short-read) | robustly (sr≥30) | non-canonical at supported coord |
| --- | ---: | ---: | ---: | --- |
| chr5  | 111 | 4 | 3 | 3 of 4 (the canonical one = J2, sr3, inconclusive) |
| chr1  | 54  | 2 | 1 | chr1:19219782 sr99 CA..AC (noncanon); chr1:223712876 sr4 GT-AG |
| chr11 | 43  | 2 | 0 | both sr≤2 (GT-AG, but detection-boundary noise) |
| chr17 | 57  | 2 | 2 | both non-canonical (GA..CC sr51, AG..CA sr32) |
| chr19 | 84  | 3 | 1 | chr19:797534 sr2643 CC..CC (noncanon); 2× sr2 GC-AG |
| **TOTAL** | **349** | **13** | **7** | **high-depth SUPPORTED are mostly NON-canonical** |

**Pattern generalizes from chr5.** Of 349 candidates, ~336 (96 %) are artifacts (no independent support).
Of the 13 SUPPORTED, only ~7 have real short-read depth (sr≥30); the other 6 sit at sr≤4 (detection-boundary
noise). And the high-depth SUPPORTED are dominated by NON-canonical / complex-locus signatures (CA..AC,
CC..CC, GA..CC — like the chr5 7805 bp "deletion-lookalike", which WGS then showed is real genomic sequence
but NOT clean canonical splicing), while the canonical-at-coord calls are all low-depth (sr≤4). So GMAP's
per-chromosome yield of well-supported, canonical, unique novel junctions is **~0-1 per chromosome** — the
chr5 "3-4 real" was, if anything, the high end. Extrapolating 349→13 across ~19 remaining autosomes: order
~1-2k candidates genome-wide, ~few-dozen SUPPORTED, of which a minority are well-supported canonical novels.
CAVEAT: these 13 have NOT had the full chr5-grade confirm (cross-aligner + WGS motif adjudication); "SUPPORTED"
= independent short-read corroboration only. A definitive genome-wide real-count needs that confirm step on
the 13 (cheap now — feed them to `dev/dna_split.py` + the WGS + `lr_probe`).

## E. CONFIRM STEP on the 13 SUPPORTED (2026-06-30) — WGS + reliable cross-aligner

Ran the chr5-grade confirm on all 13 SUPPORTED: WGS interior/flank DNA coverage (deletion test) + cross-aligner
long-read placement via `collect_junction_counts_from_bam` (the validated method; a first manual-cigar-walk
attempt was buggy — undercounted long introns). Data: `$W/confirm_13.json`, `dev/confirm13*.py`.

- **WGS: 13/13 present at NORMAL copy** (interior/flank ratio 0.72–1.33; a het del would be ~0.5, homo ~0).
  **NONE of the 13 is a deletion or SV** — every one is RNA-level (splicing), like the chr5 three. The
  "structural-deletion lookalike" worry (from the chr5 7805bp) does NOT generalize into real deletions.
- **⚠ Cross-aligner REFRAME — the harness's "gmap-only" is INFLATED by the SAME normalization gap.** Checking
  the 13 with ambiguity-tolerant (length-matched ±8) placement: **10 of 13 are actually placed by ≥2
  INDEPENDENT long-read aligners** (minimap2/deSALT/uLTRA) — a few bp off GMAP's exact coordinate, so the
  harness's exact-match STEP-1 mislabeled them "gmap-only." They are real multi-aligner junctions, NOT unique
  to GMAP → they do NOT count toward GMAP's unique value.
- **Genuinely gmap-only in long reads = only 3 of 13:** chr1:19219782 (sr99, non-canon CA..AC — STRONG
  short-read support, a real junction only GMAP found in the LR panel = GMAP earning its seat); chr11:308314
  (sr1, canon GT-AG — thin); chr19:48965691 (sr2, canon GC-AG, GMAP-LR 185 but others ~0 — GMAP-specific,
  thin short-read). So GMAP's genuinely-unique, well-corroborated real-junction yield across 5 chromosomes is
  **~1 (chr1:19219782), maybe 2** — even lower than §D's per-chromosome estimate.

**Two consequences:** (1) GMAP's true unique value is even smaller than the raw "13 SUPPORTED" implied — most
of those are shared junctions. (2) The harness itself needs the ambiguity-tolerant fix applied to its STEP-1
"others" set (currently exact-normalized), or it over-counts gmap-only candidates. Fix before the genome-wide run.

## F. CORRECTED 5-chrom re-validation (2026-07-01, ambiguity-tolerant gmap-only test — job 32220563)

Re-ran with the harness STEP-1 fix (gmap-only = exact OR same-length within ±win, not exact-normalized).
Candidate counts drop ~42 % vs the exact-match version (which inflated "gmap-only" via the normalization gap):

| chrom | candidates (corrected) | SUPPORTED | was exact-match |
| --- | ---: | ---: | --- |
| chr5  | 61  | 1 (J2, sr3, inconclusive) | 111 / 4 |
| chr1  | 23  | 2 (19219782 sr99 non-canon; 223712876 sr4 GT-AG) | 54 / 2 |
| chr11 | 33  | 1 (308314 sr1 GT-AG) | 43 / 2 |
| chr17 | 34  | 0 | 57 / 2 |
| chr19 | 50  | 1 (41860358 sr2 GC-AG) | 84 / 3 |
| **TOTAL** | **201** | **5** | **349 / 13** |

On the CORRECTED (non-inflated) candidate set, GMAP's genuinely-unique short-read-SUPPORTED yield across 5
chromosomes is **5**, of which only **chr1:19219782 (sr99)** is well-supported; the rest are sr≤4
(detection-boundary). So GMAP's genuinely-unique, well-corroborated real-junction yield ≈ **~1 across 5
chromosomes** — the fix STRENGTHENS "GMAP unique value is very low." (Note: §D/§E's 349/13 numbers were on the
inflated exact-match candidate set; §F supersedes them for the candidate/SUPPORTED counts. The §E confirm
conclusions — none are SV; most "gmap-only" are actually multi-aligner — are consistent and reinforced.)

**KEY CORRECTNESS POINT (verified 2026-07-02):** on the corrected chr5, the 3 REAL novel junctions
(SLC35A4:140564954, TMED9:177592500, SQSTM1:179824400) are **EXCLUDED from the gmap-only candidate set
entirely** — because deSALT/uLTRA place them too (ambiguity-tolerant). Only J2 (inconclusive alt-acceptor)
remains a genuinely gmap-only chr5 candidate. **So the 3 real junctions were NEVER GMAP's unique contribution
— they are real MULTI-ALIGNER junctions the exact-coordinate matching mislabeled "gmap-only."** This means
the original adjudication's very premise ("111 GMAP-only recurrent novels") was itself inflated by the
normalization gap: some of the 111 are real junctions other aligners also find. GMAP's TRUE unique yield is
the corrected §F set (201 cand → 5 SUPPORTED → ~1 well-supported), NOT the "3-4 real" I earlier attributed to
it. `adjudication_111_v3.json` should note: the 3 "real intragenic junctions" are real but multi-aligner, not
gmap-unique (does not change their reality, only the gmap-credit).

## VERDICT

**GMAP is appropriate for the panel in its current fenced, scoring-suppressed configuration — and its genuinely-
unique value is VERY LOW (~1 well-supported novel junction across 5 chromosomes).** Robustly established across
chr5 + chr1/11/17/19 with WGS + cross-aligner confirm (§A–F):
- ~96-97 % of GMAP's unique-junction output is non-canonical noise; the `score_segment` fences neutralize it
  at no accuracy cost (the fence test's primary role — VALID and load-bearing).
- On the AMBIGUITY-TOLERANT-corrected candidate set (§F, the numbers to trust): 201 gmap-only candidates across
  5 chroms → **5 short-read-SUPPORTED → ~1 well-supported** (chr1:19219782, sr99). The exact-match "349/13"
  (§D/§E) OVER-counted gmap-only via the normalization gap.
- WGS confirm: none of the SUPPORTED are deletions/SV — all real RNA splicing (§E).
- **The 3 "real intragenic novel junctions" from the original 111 adjudication (SLC35A4/TMED9/SQSTM1) are REAL
  but MULTI-ALIGNER, not gmap-unique** (deSALT/uLTRA place them; excluded from the corrected gmap-only set,
  verified §F). So GMAP earned essentially NO unique credit for them.
- Net: GMAP is a LOW-unique-value member, kept ONLY because its noise is fully fenced and it occasionally
  (~1/5-chrom) surfaces a real short-read-corroborated junction the other LR aligners miss (chr1:19219782).
  The earlier "~3-4 real on chr5 / needs genome-wide count" framing is SUPERSEDED — those 3 weren't gmap-unique.

## RECOMMENDATION (informing the native-aligner 5→2-3 cut; benchmark-earned, not a directive)

1. **Keep GMAP, fenced** — its noise is controlled; do NOT relax the `score_segment` non-canonical fences
   (the fence test's behavioral invariants must keep passing).
2. **GMAP is now a well-supported DROP candidate once C1 lands.** Corrected 5-chrom evidence shows its unique
   real-junction yield is ~1 per 5 chromosomes (~a handful genome-wide) — and even that (chr1:19219782) is a
   non-canonical junction a calibrated empirical-penalty realigner (C1) is designed to recover. Annotate
   `test_gmap_fence_regression.py:8`'s "~100 genuine" as an UNVALIDATED chr5 candidate count (corrected: ~1
   genuinely-unique well-supported). Optional genome-wide confirm: the harness is ready (`--chrom all`), but
   needs a genome-wide GMAP+independent LR panel first (current BAMs cover only chr5 + morechrom chr1/11/17/19);
   the 5-chrom result is already decision-grade, so a full genome-wide run is low marginal value.
3. **The 3 real (multi-aligner) non-canonical junctions are a concrete acceptance test** for whether the panel
   — especially the C1 realigner — can recover genuine non-canonical human splicing. Handed to the C1 agent
   (`.claude/inbox/`); in `rectify/data/validation/CASE_STUDIES.md` (de-novo aligner acceptance criteria).

## KEEP GMAP FOR TESTING, treat as a production DROP candidate (2026-07-02)

GMAP's residual value has shifted from "production consensus member" to "test/dev reference." Keep it in the
DEV/TEST panel even if it is dropped from production, because its distinctive properties make it uniquely
useful for validation:
1. **Sensitivity ceiling.** GMAP is BLAST-like (no seed-and-chain requirement) → it maps divergent/spliced
   reads the seed-and-chain members (minimap2 lineage incl. uLTRA/gapmm2; deSALT; mapPacBio) cannot. That is
   simultaneously its 97 %-noise source AND the only view of the discovery frontier — you cannot know what the
   panel MISSES without a maximally-sensitive reference.
2. **Fence adversary.** GMAP's ~878k-non-canonical mass is the real-world negative control that proves the
   `score_segment` fences work (`test_gmap_fence_regression.py`). Drop GMAP from the test panel and you lose
   the adversary that catches a future scoring change leaking noise into consensus.
3. **The C1 drop-test.** The clean way to retire GMAP from production is to show C1 + minimap2 +
   (deSALT|mapPacBio) recovers GMAP's genuine unique finds (chr1:19219782-class) WITHOUT its noise. That test
   REQUIRES GMAP in the test panel as the thing C1 must subsume.
4. **Benchmark FDR calibration.** GMAP's non-canonical output is a known FP generator → ideal for validating
   the simulation benchmark's non-canonical FDR control against a real adversary, not just synthetic noise.

**Cost side (why NOT production):** GMAP is the slowest, noisiest member, pinned to 2021-05-27 (2024 is broken
— SIGSEGV + ~100× slower). Paying that per-read on every production run to net ~1 unique real junction / 5
chromosomes is a poor trade; the fences spend effort suppressing a mostly-noise member.

**STATED LIMITATION of this whole validation (important):** "real" was gated on INDEPENDENT SHORT-READ
corroboration. A genuine LONG-READ-ONLY GMAP junction (a real isoform in a region short reads cover poorly, or
a true junction the short-read panel also misses) is scored as "artifact" here because we cannot independently
confirm it. So the "~1 unique real / 5 chrom" figure is a LOWER BOUND biased toward what short reads see — NOT
a hard ceiling. This is itself a reason to keep GMAP for testing (to probe exactly that long-read-only tail
once C1 / the benchmark can adjudicate without short-read dependence).

**Net recommendation:** keep GMAP fenced + pinned in the DEV/TEST panel; execute the production drop ONLY after
C1 demonstrably recovers its unique finds on the benchmark. Until C1 exists, status quo (GMAP fenced, in the
production panel) is fine — its noise is controlled and we cannot yet prove nothing real is lost.

## Sources / provenance
- A549 chr5 GMAP buckets: `/scratch/users/kevinroy/deliverable_b/gmap_corroboration.json` (deliverable_b).
- 111 deep validation: `dev/COMPASS_2corroborated_CROSSPLATFORM.md`, `$OAK/adjudication_111_v3.json`.
- Version experiment + fences: `dev/HANDOFF_SHORTREAD_P5.md`, `tests/test_gmap_fence_regression.py`.
