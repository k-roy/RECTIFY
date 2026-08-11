# SPIKE-IN PRECISION BENCHMARK — RESULT (task #20)

**Agent:** spike-in agent · branch `worktree-agent-a25a2c1e784ad37dc` · started 2026-07-09
**Mission:** the ONLY unbiased ground-truth measure of the native re-aligner's non-canonical
junction PRECISION + RECALL, motif- and decoy-distance-stratified across the four arms
(raw-mm2 / arm-A motif-biased / arm-B motif-blind / arm-Bguard). Snaptron/recount3 is CIRCULAR
(STAR-built, motif-filtered) so it cannot prove fabrication — the spike-in gives ABSOLUTE truth
because we build the junctions.

**Why now (the motivating finding, dev/SUMNER_SMA_LEADS_B_splicing.md):** on real SMA data the
re-placer DRIFTS/mis-places real junctions at NON-HOMOPOLYMER positions 6-28bp (SNRPN ~6bp,
UBA1 ~26bp, PCBP2 ~28bp). The HP-drift guard (arm-E) MISSES this general non-HP drift. The
spike-in must measure exactly this: RECALL (does motif-blind recover non-canonical junctions
mm2 flattens) and PRECISION/FABRICATION (junctions each arm calls that are NOT true, per
decoy-distance, especially NON-HP), and whether the guard helps or misses the non-HP drift.

---

## PLAN (written before compute)

Reuse the SIRV-model harness in `scripts/benchmark/noncanon_sim/`:
- `build_panel.py` — synthetic contigs carrying KNOWN non-canonical junctions + canonical decoys
  at graded distances. NON-HP DECOY-DISTANCE is the priority axis (the non-HP drift is the real
  failure). Need to add a decoy-distance sweep panel (non-HP R3 non-canonical acceptor with the
  canonical AG decoy at k = 2,4,6,8,12,20 bp — spanning the 6-28bp real-data drift range).
- `gen_reads.py` — pbsim3 ERRHMM-ONT error model (validated), full-length templ mode.
- `run_arms.py` — align (minimap2 -ax splice) then refine to arms A/B/C/E. arm-A (motif-biased,
  incumbent) is the MANDATORY control: motif-blind contribution = arm-B − arm-A.
- SCORE on the synthetic contigs (known truth): recall + precision/FDR + canonical-snap
  fabrication via `rectify/core/benchmark/scorer.py` (already ambiguity-aware, already stratifies
  by class + canonicity + fp_canonical_snap). Stratify additionally by decoy-distance rung.

### The four arms (SPEC-mandated ladder)
- **raw-mm2** = `aligned.sorted.bam` (minimap2 -ax splice, no refinement) — the flattening baseline.
- **arm-A** = refine motif_blind=False (incumbent motif-biased re-placer).
- **arm-B** = refine motif_blind=True (motif-blind — the anti-flattening lever).
- **arm-Bguard (arm-E)** = motif_blind=True + hp_drift_margin>0 (the shipped guard).

### Headline metrics per arm, per decoy-distance rung, motif-stratified
- **RECALL** = fraction of spiked non-canonical junctions placed at the TRUE site (normalize_junction
  ambiguity-aware) vs flattened to the canonical decoy. Curve: raw-mm2 (worst) → arm-A → arm-B (best).
- **PRECISION / FABRICATION** = junctions called that are NOT true, per decoy-distance. The
  canonical-snap FP (`fp_canonical_snap`) = the flattening fabrication signature. Does arm-Bguard
  cut fabrication vs arm-B WITHOUT cutting recall, and does it MISS the non-HP drift (Snaptron finding)?

### Decision (background choice)
SPEC offers LongBench RNA004 / SG-NEx / reuse-BY4742. Per SPEC note: "the background genome is
irrelevant to the spiked-contig truth" — the spiked reads carry SIMULATED (validated) error, the
background only supplies run-context realism. FIRST PASS = pure synthetic contigs (SIRV-model,
clean truth separation), no real background mixed in yet — this isolates the precision/recall
signal on the spiked junctions with zero confound. (Real-background mixing = external-validity
add-on, notably NOT needed for the truth measurement itself.)

---

## ★ ADVISOR CORRECTION (2026-07-09) — the recall panel cannot measure fabrication

The R3(non-canonical)+canonical-decoy panel measures RECALL (arm-B recovers what arm-A flattens)
— keep it, mission-required. But it STRUCTURALLY CANNOT exhibit arm-B fabrication: reads are
spliced at the true non-canonical site, so the true placement matches perfectly and arm-B holds.
On that panel arm-B fabrication ~= 0 and arm-E == arm-B (guard never fires). A clean precision
table there means "the panel contains no failure," NOT "arm-B is precise" — a FALSE all-clear.

The Snaptron finding is the OPPOSITE direction: a REAL (canonical/dominant) junction gets DRIFTED
to a false NON-CANONICAL site at a non-HP position (6-28bp). Truth=canonical; false call=non-canonical.
To measure fabrication I need the MIRROR panel:
  - canonical TRUE junction, reads generated from it;
  - a nearby non-canonical DRIFT TARGET at graded NON-HP distance (matching SNRPN/UBA1/PCBP2 6-28bp);
  - the drift target must NOT be registered as a true isoform (else a drifted call scores TP not FP);
  - FABRICATION metric = NON-CANONICAL FP (by_canon["noncanonical"]["fp"]), NOT fp_canonical_snap
    (fp_canonical_snap is arm-A's FLATTENING signature; arm-B's fabrication is a non-canonical FP
    near a canonical truth — a DIFFERENT number).

MECHANISM (how fabrication arises, from the refiner code): the candidate pool = junctions OTHER
reads placed (within max_boundary_shift=50bp) + annotated. minimap2 error-driven boundary variance
populates the pool with slightly-shifted junctions; arm-B re-places each read to the best-scoring
pool member by its OWN error-laden evidence -> can land on a shifted non-canonical position = the drift.

REVISED PLAN: run a cheap MECHANISM PROBE FIRST (canonical-truth contig, error-laden pbsim reads,
arms A/B/E) — does arm-B actually emit non-canonical FP off the true site? If arm-B holds even under
elevated error, THAT is a first-order finding (pure sim under-represents real non-HP drift; real-
error-on-non-canonical needs wet-lab SIRVs, out of scope). Only build the full decoy-distance sweep
once arm-B fabricates at least once. Score with scorer.py:score_bam (build ReadTruth via
JunctionTruth.from_intron) — gives ambiguity-aware precision/FDR by canonicity = the fabrication metric.

## CHECKPOINTS (append-only; every number persisted the moment it lands)

- [x] harness + refiner + scorer read; advisor flagged the panel-direction flaw (FIXED above)
- [ ] MECHANISM PROBE: canonical-truth contig + drift target, arms A/B/E — does arm-B fabricate?
- [ ] RECALL panel (R3+decoy): recall per arm, decoy-distance + motif stratified
- [ ] FABRICATION panel (canonical-truth + non-HP drift target sweep): non-canonical FP per arm per distance
- [ ] scored via scorer.py:score_bam; arm ladder raw-mm2/A/B/E
- [ ] verdict: does arm-B recover (recall win) & at what fabrication cost & does the guard miss non-HP drift

## CHECKPOINT 1 (2026-07-09) — mechanism probe on FEATURELESS panel: arm-B barely fabricates

Built fab_probe (canonical GT-AG truth at (120,240), non-canonical drift target TC at +10bp=250,
HP-scrubbed to non-HP everywhere), pbsim3 ERRHMM-ONT reads.
- acc=0.95: raw-mm2 280 TRUE / 0 FP; arm-A 6 FP@(120,251); arm-B 1 FP@drift + 5@(120,241); arm-E == arm-B.
- acc=0.88: raw-mm2 232 TRUE / 0 FP; arm-A == arm-B == arm-E, all 7 reads -> (120,241) [1bp shift].
- SEED reads (spliced at true drift 250) were placed by minimap2 at (120,248) — mm2 PRE-FLATTENED the
  non-canonical drift target 2bp. So the pool's drift member was at 248, not 250.

★ DIAGNOSIS (advisor-confirmed, I OWN the correction): on RANDOM / HP-scrubbed sequence a read spliced
at canonical A CANNOT re-place to A+k without ~k mismatches (align exon2 genome[A:A+e3] starting at A+k
costs the mismatches between genome[A:A+k] and genome[A+k:A+2k]). On featureless sequence that's ~k
mismatches -> true site ALWAYS wins, no error rate tips it. arm-A==arm-B==arm-E all agreeing (the
project's own ceiling/zero-evidence trap) = the panel has NO discriminating stratum. Error alone can't
produce non-HP drift; the ENABLING variable is LOCAL SEQUENCE SIMILARITY (microhomology/near-repeat) at
the drift distance — exactly what highly-expressed/ribosomal genes (the SMA fabrication-enriched set)
are rich in, and exactly what my `_no_run` scrub removed. NOT a scope finding — a featureless-panel artifact.

FIX: add microhomology to spikein_fab.py — make exon2[0:k] a NEAR-repeat (period-k, imperfect, non-HP)
of genome[A+k:A+2k], drift dinuc forced non-canonical. Verify normalize_junction(D,A) != (D,A+k) (genuine
FP, not ambiguity slide — needs IMPERFECT similarity). Let pbsim error decide which reads drift; count any
non-true non-ambiguity N-op as non-canonical FP (score_bam by_canon["noncanonical"]["fp"]). If arm-B then
fabricates & arm-E==arm-B (guard blind, non-HP) while arm-A differs -> real benchmark; sweep 6-28bp.

## CHECKPOINT 2 (2026-07-09) — microhomology ENABLES the drift; but arms don't differ (reads start canonical)

Added --microhom-mismatch to spikein_fab.py: genome[A:A+k] := imperfect near-copy of genome[A+k:A+2k]
(installs the local sequence similarity real drift-prone loci carry). Genuine-FP guard: verified
normalize_junction(D,A) != (D,A+k) (not an ambiguity slide). Panel: drift=10bp, microhom_mismatch=0.2
(7/10 match), seed 5%, pbsim acc 0.92.

RESULT (scored, true-junction reads only, seed reads excluded):
  raw-mm2 : TRUE 352 / DRIFT@+10 non-canon FP = 0 / other 0
  arm-A   : TRUE 335 / DRIFT@+10 = 4 / other 13 @(120,241)
  arm-B   : TRUE 335 / DRIFT@+10 = 4 / other 13 @(120,241)   [IDENTICAL to A]
  arm-E   : TRUE 335 / DRIFT@+10 = 4 / other 13 @(120,241)   [IDENTICAL to B — guard blind, non-HP]

★ FINDINGS:
1. Microhomology WORKS: arm-B now fabricates 4 non-canonical FP @+10bp that raw-mm2 (0 FP) does NOT.
   Confirmed the drift target (120,250) is in the pool (13 seed reads land there). The re-placer DOES
   manufacture non-canonical junctions off a canonical truth — the Snaptron mechanism reproduced in sim.
2. arm-A == arm-B == arm-E IDENTICALLY. The fabrication is NOT motif-blind-specific here: because the
   reads START at a canonical junction (tier<4), the refiner's `tier_beats_alt` canonical-preference
   NEVER engages (it only fires when the CURRENT junction is non-canonical, tier>=4). So arm-A's motif
   bias is dormant for canonical-origin reads -> A behaves like B. The drift is driven purely by
   read-evidence competitiveness (the microhomology), which is arm-invariant.
3. arm-E (guard) == arm-B: the guard is non-HP-blind BY CONSTRUCTION — it only vetoes moves INTO a
   homopolymer run; this drift is at a scrubbed NON-HP position, so the guard never fires. This
   DIRECTLY confirms the mission's key concern: the HP-drift guard MISSES the general non-HP drift.

IMPLICATION for the arm ladder: the mission's headline (motif-blind CONTRIBUTION = arm-B - arm-A) is
about RECALL (does motif-blindness recover non-canonical junctions arm-A flattens). Fabrication, by
contrast, is a property of the re-placement ENGINE (shared A/B/E) on canonical-origin reads, NOT of
motif-blindness. So the deliverable needs BOTH: the RECALL panel (reads at non-canonical truth ->
arm-A flattens, arm-B recovers = the motif-blind win) AND this FABRICATION panel (reads at canonical
truth -> the engine drifts to non-canonical = the precision cost, arm-invariant, guard-blind).

NEXT: (a) sweep drift distance 6-28bp x microhom strength on the fabrication panel -> fabrication rate
curve; (b) run the RECALL panel (existing hp_power/full harness) for the arm-B-recovers-flattened
number; (c) integrate via scorer.py:score_bam for the ambiguity-aware by-canonicity FP metric.

## CHECKPOINT 3 (2026-07-09) — FABRICATION SWEEP (drift 6-28bp): the precision cost number

Panel fab_sweep: 5 contigs (drift 6,10,14,20,28bp), microhom_mismatch=0.2, seed 5%, pbsim acc 0.92,
400 true-junction reads/contig (2000 scored). Scored via scorer.py:score_bam (ambiguity-aware, only
true canonical junctions registered -> a drifted call = non-canonical FP).

OVERALL (scorer.py, by-canonicity):
  arm       recall  prec    FDR     TP    FP  FP_noncanon FP_canon  FN
  raw-mm2   0.9325  0.9941  0.0059  1865  11  0           11        135
  arm-A     0.9050  0.9648  0.0352  1810  66  63          3         190
  arm-B     0.9010  0.9606  0.0394  1802  74  74          0         198
  arm-E     0.9010  0.9606  0.0394  1802  74  74          0         198

PER-DRIFT-DISTANCE (reads landing at the exact +k non-canonical target / total):
  dist=6:  raw 0 drift, A/B/E all 5 @+6   (fab_rate raw 1.6% -> A 2.3% / B 3.9%)
  dist=10: raw 0, A/B/E 7 @+10            (raw 0.3% -> A 4.0% / B 4.3%)
  dist=14: raw 0, A 0 @+14 / B/E 4 @+14   (raw 0.8% -> ~3.8% ; A deflects the exact target)
  dist=20: raw 0, A/B/E 3 @+20            (raw 0.3% -> A 4.2% / B 4.5%)
  dist=28: raw 0, A/B/E 5 @+28            (raw 0.0% -> A/B/E 3.2%)

★★ HEADLINE FABRICATION FINDINGS:
1. raw-mm2 (minimap2) is NEARLY PERFECT: FDR 0.6%, ZERO non-canonical FP. It places canonical-origin
   reads at the true site; it does NOT fabricate non-canonical junctions. (This is the flip side of
   the recall story: minimap2's precision comes from its canonical bias / motif snapping.)
2. THE RE-PLACER FABRICATES at ~3.5-4% junction-FDR across ALL drift distances 6-28bp. At every
   distance some canonical-origin reads drift to the EXACT non-canonical target. The Snaptron real-data
   drift (SNRPN/UBA1/PCBP2, 6-28bp) is REPRODUCED in ground-truth sim — it is a genuine, unbiased,
   non-circular measurement of the fabrication the motif-agnostic short-read oracle could NOT prove.
3. arm-E (HP-drift guard) == arm-B EXACTLY at every distance: the guard reduces the non-HP drift by
   ZERO. Ground-truth confirmation of the mission's key concern — the HP-drift guard MISSES the general
   non-HP drift (it only vetoes moves INTO a homopolymer run; this drift is non-HP by construction).
4. arm-A (motif-biased) fabricates SLIGHTLY LESS than arm-B (63 vs 74 non-canonical FP; at dist=14 arm-A
   places 0 at the exact +14 target vs arm-B's 4). arm-A's canonical tiebreaker occasionally deflects a
   read off the non-canonical target -> the motif bias buys a LITTLE precision. But it is a WEAK effect
   (both ~4% FDR) because canonical-origin reads keep tier<4 so the canonical preference is mostly dormant.
   => motif-blindness (arm-B) is NOT the primary fabrication source; the re-placement ENGINE is.

## CHECKPOINT 4 (2026-07-09) — RECALL PANEL: the motif-blind win (arm-B − arm-A)

Panel recall_panel: build_panel --panel full (R0/R1/R3 x context x decoy), e5/e3=80, intron 100,
cryptic 400 reads/R3 cell, cryptic-frac 0.3 (abundant WT), pbsim acc 0.92 -> 6166 reads. R3 = non-
canonical acceptor truth + co-located canonical decoy = the FLATTENING case. Arms A/B/C/E refined
(run_arms bg rc=0, all quickcheck OK).

RECALL of R3 NON-CANONICAL reads (placed at the TRUE non-canonical site, normalize_junction
ambiguity-aware, vs snapped to the canonical decoy):
  raw-mm2 : 1188/1433 = 0.8290
  arm-A   :  344/1433 = 0.2401   <-- motif-biased re-placer FLATTENS (worse than raw mm2!)
  arm-B   : 1183/1433 = 0.8255   <-- motif-blind RECOVERS (restores raw-mm2 recall)
  arm-E   : 1183/1433 = 0.8255   <-- guard doesn't touch discovery

★★ HEADLINE RECALL FINDING — the motif-blind contribution (arm-B − arm-A) = +58.5 percentage points
(0.826 − 0.240). The INCUMBENT motif-biased refiner (arm-A) actively SNAPS non-canonical junctions onto
the canonical decoy, collapsing recall from raw-mm2's 82.9% to 24.0% — it is WORSE than doing nothing
(raw minimap2) on non-canonical discovery. Motif-BLIND re-placement (arm-B) undoes that flattening and
restores recall to raw-mm2's level. This is the make-or-break motif-blind win, on ground truth immune to
the STAR/motif circularity that rigs the short-read Q2. arm-E==arm-B confirms the guard preserves
discovery at zero cost (mission's "guard should not cut recall").

## CHECKPOINT 5 (2026-07-09) — NO-SEED fabrication control [SUPERSEDED by CHECKPOINT 9's decomposition]

Panel fab_noseed: identical to fab_sweep but seed-frac=0.0 (NO contaminant reads priming the pool at the
drift target). Question: does the re-placer fabricate WITHOUT me seeding the pool?

  arm       recall  prec    FDR     FP_noncanon  (vs seeded fab_sweep FP_noncanon)
  raw-mm2   0.9365  0.9936  0.0064  0            (0)
  arm-A     0.9195  0.9756  0.0244  42           (63)
  arm-B     0.9155  0.9714  0.0286  54           (74)
  arm-E     0.9155  0.9714  0.0286  54           (74)

INITIAL (over-broad) reading — CORRECTED in CHECKPOINT 9: even with ZERO pool seed the total non-canonical
FP stays 42-54 (FDR 2.4-2.9%) vs raw-mm2's 0. BUT CHECKPOINT 9's decomposition shows this residual is
ENTIRELY the generic 1bp boundary noise; the SMA-relevant k-bp drift is 0.00% without the seed (POOL-
DEPENDENT). So "intrinsic" applies ONLY to the 1bp noise, NOT to the 6-28bp SMA drift. arm-E==arm-B (guard
blind, non-HP).

## ★ ADVISOR CORRECTION 2 (2026-07-09) — raw-mm2 >= arm-B everywhere; the "win" is arm-B undoing arm-A

CRITICAL framing fix. In EVERYTHING measured, raw minimap2 >= arm-B:
 - recall panel: raw-mm2 0.829 vs arm-B 0.826 (arm-B MATCHES, does not beat mm2).
 - fabrication panel: raw-mm2 prec 0.994/rec 0.933 vs arm-B 0.961/0.901 (raw-mm2 BEATS arm-B).
The +58pt "motif-blind win" is arm-B vs the INCUMBENT REFINER (arm-A), NOT vs minimap2. arm-B undoes
damage arm-A inflicts and lands back at minimap2's baseline. Honest statement: "the incumbent motif-
biased refiner is net-HARMFUL to non-canonical recall; motif-blind repair restores it to raw-mm2." It is
NOT "arm-B recovers junctions minimap2 flattens" — minimap2 barely flattens at these settings (~17-20%).

## CHECKPOINT 6 (2026-07-09) — Check 1: recall stratified by rung x decoy-presence (no recompute)

recall_panel R3/R0/R1 x {decoy, nodecoy}, recovery of the TRUE junction:
  cell            raw-mm2   arm-A    arm-B
  R0 decoy         0.981    0.976    0.972
  R0 nodecoy       0.901    0.960    0.885
  R1 decoy         0.955    0.976    0.950
  R1 nodecoy       0.967    0.974    0.952
  R3 decoy         0.796    0.000    0.796   <-- arm-A TOTAL flatten (0%!) w/ canonical decoy; arm-B == raw
  R3 nodecoy       0.862    0.483    0.855   <-- arm-A loses HALF even w/o decoy; arm-B == raw

MECHANISM (advisor's suspicion confirmed + refined): arm-A's flattening is decoy-DRIVEN (0% R3 recall
when a canonical decoy exists — it snaps EVERY non-canonical junction onto the decoy) but ALSO loses
no-decoy R3 (0.483) — so it is a broad "prefer canonical" bias, not purely decoy-snapping. arm-B recovers
to raw-mm2's level in BOTH R3 cells. Canonical R0/R1 unaffected across arms (~0.95-0.98). The prior R0-HP
cells in mix_full_out/trade_curve.json also show A==B==C ~0.87 (canonical, no snap) — consistent.

NEXT: Check 2 — sweep decoy-offset {3,6,10,14,20,28} to find whether minimap2 EVER flattens hard enough
that arm-B BEATS raw-mm2 (the regime the mission's headline needs). If none, the premise "motif-blind
recovers what minimap2 flattens" is challenged and must be reported as such.

## CHECKPOINT 7 (2026-07-09) — Check 2: recall vs decoy DISTANCE sweep (does arm-B EVER beat raw-mm2?)

6 panels, decoy-offset {3,6,10,14,20,28}, R3 non-canonical + canonical decoy, pbsim acc 0.92, ~250
cryptic reads/cell. R3 non-canonical recall:
  decoy_off  raw-mm2  arm-A   arm-B   B-raw    B-A
     3        0.818    0.233   0.813   -0.005   +0.580
     6        0.861    0.233   0.851   -0.010   +0.618
    10        0.837    0.232   0.830   -0.008   +0.598
    14        0.837    0.232   0.833   -0.004   +0.600
    20        0.837    0.235   0.834   -0.003   +0.598
    28        0.838    0.261   0.831   -0.007   +0.570

★★ DECISIVE: arm-B NEVER beats raw-mm2 at ANY decoy distance (B-raw = -0.003 to -0.010, always ~1pt
BELOW). raw-mm2 recall is a flat ~0.82-0.86 across ALL distances — minimap2 does NOT meaningfully flatten
non-canonical junctions in THIS simulator (it places ~83% at the true site regardless of decoy proximity).
arm-A flattens to ~0.23-0.26 everywhere; arm-B recovers to the raw-mm2 baseline. The motif-blind win
(B-A = +0.57 to +0.62) is large and consistent — but it is entirely arm-B repairing the incumbent
refiner's self-inflicted damage, NOT arm-B recovering junctions minimap2 flattened.

INTERPRETATION: the mission's literal Q ("does motif-blind recover non-canonical junctions minimap2
FLATTENS") is NOT demonstrable on this pbsim sim because minimap2 barely flattens here. This is the
known simulator-realism gap (NONCANON_BENCHMARK_SURVEY: "pbsim over-dispersion ~13x too low vs real
SIRV; the 0.09->1.07 trap; sim wins necessary not sufficient -> must transfer to SIRV/real"). The prior
project claim "minimap2 flattens 40-100%" came from ERROR-FREE reads / a different geometry (min-align
seeding), not this pbsim-error regime. HONEST verdict below states this plainly.

## CHECKPOINT 8 (2026-07-09) — microhomology-strength robustness (fabrication is not a single-parameter artifact)

drift=10, arm-B non-canonical fabrication vs microhom_mismatch:
  mh=0.10 -> FDR 5.8% (22 FP) ; mh=0.20 -> FDR 4.3% (from sweep) ; mh=0.35 -> FDR 5.5% (21 FP)
Fabrication persists ~4-6% across the microhomology range (not monotone, but robust — not a knife-edge
of one parameter value). mh=0.0 is REJECTED by the genuine-FP guard (perfect repeat = ambiguity slide).

## CHECKPOINT 9 (2026-07-09) — FP DECOMPOSITION (advisor-flagged conflation fixed)

Re-scored fab_sweep + fab_noseed BAMs, splitting non-canonical FP into k-bp-exact-drift vs 1bp-nearshift:
  SEEDED  : arm-B k-bp-drift=24 (1.28%) | 1bp-nearshift=50 (2.66%)
  NO-SEED : arm-B k-bp-drift= 0 (0.00%) | 1bp-nearshift=54 (2.86%)
The ~4% total was CONFLATING the SMA k-bp drift (~1.2%, POOL-DEPENDENT -> 0% no-seed) with generic 1bp
boundary noise (~2.5%, intrinsic, NOT the 6-28bp phenomenon). Verdict corrected: SMA-relevant number = ~1.2%.

---

# ★★★ FINAL VERDICT — spike-in precision benchmark (task #20)

## The consolidated arm ladder (ground truth, unbiased by the STAR/motif circularity)

### RECALL (non-canonical junction recovery) — recall_panel + decoy-distance sweep
| arm | R3 non-canonical recall | vs raw-mm2 | vs incumbent arm-A |
|---|---|---|---|
| raw-mm2 (minimap2) | 0.82-0.86 (flat across decoy 3-28bp) | baseline | — |
| arm-A (motif-BIASED refiner, incumbent) | 0.23-0.26 | **-0.58** (net HARMFUL) | baseline |
| arm-B (motif-BLIND refiner, shipped) | 0.81-0.85 | -0.003 to -0.010 (~= raw) | **+0.57 to +0.62** |
| arm-E (motif-blind + HP-drift guard) | == arm-B | == arm-B | == arm-B |

### PRECISION / FABRICATION (drift of a CANONICAL junction to a false non-canonical site) — fab_sweep + no-seed control
**★ The total non-canonical FP DECOMPOSES into two distinct phenomena (do NOT conflate):**
| component | seeded (5%) | no-seed | what it is |
|---|---|---|---|
| **k-bp DRIFT to exact non-canonical target (the SMA phenomenon, 6-28bp)** | arm-B **1.28%** / arm-A 1.06% | **0.00%** | the SNRPN/UBA1/PCBP2 mirror — POOL-DEPENDENT |
| 1bp boundary nearshift (generic refiner noise) | arm-B 2.66% / arm-A 2.45% | arm-B 2.86% / arm-A 2.44% | ordinary boundary jitter, NOT 6-28bp drift |
| raw-mm2 (either component) | 0 | 0 | minimap2 does not fabricate |
| arm-E (HP-drift guard) | == arm-B on BOTH | == arm-B | **guard reduces BOTH by ZERO** |

**The make-or-break SMA-relevant number = the k-bp drift = ~1.2% (arm-B), and it is POOL-DEPENDENT** (no-seed
= 0.00%): the re-placer only drifts a canonical read to a distant non-canonical site if that site ALREADY
EXISTS in the candidate pool. In real data the pool member comes from paralogs / other isoforms / minimap2
boundary variance (exactly the highly-expressed / ribosomal / SMN-paralog loci the SMA finding flagged). The
generic 1bp nearshift (~2.5%) IS intrinsic (present with or without seed) but is boundary noise, NOT the
6-28bp SMA drift. Robust ~4-6% TOTAL over microhom 0.1-0.35.

## ★ TOP LINE FOR THE PI (one sentence)
On this ground-truth sim, **arm-B (motif-blind, shipped) is strictly dominated by raw minimap2** — equal
non-canonical recall (~0.83), slightly worse precision — BUT this sim NEVER produces the minimap2-flattening
regime arm-B was built for (pbsim under-disperses ~13x vs real ONT), so the benchmark is STRUCTURALLY UNABLE
to test arm-B's intended benefit. What it DOES prove on truth: (a) the incumbent motif-BIASED refiner (arm-A)
is net-HARMFUL to non-canonical recall (flattens 0.83->0.23); (b) motif-blind fabricates a SMALL, pool-
dependent k-bp non-HP drift (~1.2%); (c) the HP-drift guard has ZERO coverage of that non-HP drift.

## The four answers the mission asked

1. **Does motif-blind (arm-B) RECOVER non-canonical junctions minimap2 flattens (recall win)?**
   NOT AS LITERALLY POSED. In this pbsim-error ground-truth sim, **minimap2 does NOT meaningfully flatten**
   (raw-mm2 recall ~0.83 flat across all decoy distances). arm-B matches minimap2, never beats it. The
   large motif-blind gain (+58-62pt) is entirely arm-B UNDOING the INCUMBENT REFINER's (arm-A) self-
   inflicted flattening — arm-A snaps non-canonical junctions to canonical (0% R3 recall with a decoy,
   48% without), which is WORSE than raw minimap2. So the ground-truth finding is: **the incumbent motif-
   biased re-placer is net-harmful to non-canonical recall; motif-blind re-placement (arm-B) is REQUIRED
   just to restore minimap2's baseline.** (Whether motif-blind beats minimap2 needs the real-flattening
   regime — error-free reads or real SIRV/dRNA — the pbsim over-dispersion is ~13x too low to reproduce
   minimap2 flattening; NONCANON_BENCHMARK_SURVEY's documented sim-realism gap.)

2. **At what FABRICATION cost?** DECOMPOSED (do not conflate): the SMA-relevant **k-bp non-HP drift
   (6-28bp) = ~1.2% (arm-B), and it is POOL-DEPENDENT** (0.00% without a pool member at the drift target);
   plus a generic, intrinsic **~2.5% 1bp boundary-nearshift** noise that is NOT the SMA phenomenon. raw
   minimap2 fabricates ~0 of either. The ~1.2% pool-dependent drift is the first UNBIASED, non-circular
   measurement of exactly the drift the Snaptron/recount3 oracle could not prove (it's motif-biased). The
   Sumner SMA real-data drift (SNRPN/UBA1/PCBP2, 6-28bp) is REPRODUCED in ground truth AND its mechanism
   pinned: it fires only when the false non-canonical site pre-exists in the pool (paralogs / other
   isoforms / mm2 boundary variance — the high-expression/ribosomal/SMN-paralog loci the SMA finding flagged).

3. **Does motif-blindness DRIVE the fabrication?** Only WEAKLY. arm-B fabricates slightly more than arm-A
   (74 vs 63 non-canonical FP), but both fabricate at ~similar rates. The fabrication is a property of the
   RE-PLACEMENT ENGINE (drift enabled by local microhomology + ONT error), NOT of motif-blindness — because
   canonical-origin reads keep the current junction canonical (tier<4), so arm-A's motif tiebreaker is
   mostly dormant. arm-A's motif bias buys a LITTLE precision but at a CATASTROPHIC recall cost (finding 1).

4. **Does the HP-drift guard (arm-E) help or MISS the non-HP drift?** It MISSES it COMPLETELY. arm-E ==
   arm-B EXACTLY on both fabrication components at every drift distance. This is TRUE BY CONSTRUCTION —
   the panel is HP-scrubbed (max run 2) so the guard's code path (veto a move INTO a homopolymer run) is
   never even entered; the empirical BAM-identity confirms but adds nothing beyond the construction. That
   IS the point: **the shipped HP-drift guard, by design, has ZERO coverage of non-HP drift — and the
   real-data drift (SNRPN/UBA1/PCBP2) is at non-HP positions.** A general drift/precision control (short-
   read-anchored or expression-aware, per the SMA doc) is needed before genome-wide novel discovery is
   trustworthy. (A stronger optional follow-up: an HP-ADJACENT drift contig where the guard DOES fire and
   help, alongside the non-HP contig where it can't — would show the coverage GAP rather than assume it.)

## The make-or-break number
**The SMA-relevant non-canonical k-bp DRIFT FDR of the motif-blind re-placer (arm-B) = ~1.2%, at non-HP
positions 6-28bp, POOL-DEPENDENT (0% without a pool member at the drift site), with the HP-drift guard
reducing it by ZERO.** (Plus a separate ~2.5% intrinsic 1bp boundary noise that is NOT the SMA drift.) raw
minimap2 = ~0. This is the precision cost of motif-blind discovery, measured on ground truth immune to the
motif circularity that made the Snaptron oracle one-way (confirm-only).

## Honest caveats (stated, not smoothed)
- Simulator realism: pbsim ERRHMM under-disperses vs real ONT (~13x); minimap2 does not flatten here, so
  the "motif-blind beats minimap2" regime is NOT demonstrated. The recall win is arm-B vs arm-A only.
  Real-error-on-non-canonical (the true minimap2-flattening + real-drift regime) needs wet-lab SIRVs /
  real dRNA (out of scope per SPIKEIN_DESIGN); this sim proves the MECHANISM + magnitude-at-a-parameter.
- The ~4% fabrication is at microhom_mismatch=0.2 and seed 5% (magnitude at a plausible-but-arbitrary
  local-similarity strength); it is robust 4-6% over microhom 0.1-0.35 and 2.4-2.9% with no seed, so it
  is an existence-proof + calibrated-range, not a single real-world rate. Real per-locus fabrication
  depends on the actual local sequence (ribosomal/high-expression loci carry more microhomology).
- All arms scored ambiguity-aware (normalize_junction); the genuine-FP guard rejects ambiguity-slides.

## Reproduce
- Fabrication: `spikein_fab.py --drift-dists 6,10,14,20,28 --microhom-mismatch 0.2 --seed-frac 0.05` ->
  `gen_reads.py --pbsim-accuracy 0.92` -> `run_arms.py --hp-drift-margin 2.0` -> `spikein_score.py`.
- Recall: `build_panel.py --panel full --decoy-offset <k>` -> gen_reads -> run_arms -> per-rung recall.
- No-seed control: `--seed-frac 0.0`. Data: fab_sweep/, fab_noseed/, recall_panel/, dsweep_k*/.

STATUS: COMPLETE.
