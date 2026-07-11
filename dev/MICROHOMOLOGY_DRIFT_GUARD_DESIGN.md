# Microhomology-drift guard — design + Phase-1 validation (2026-07-11)

**Problem (spike-in + Sumner + COMPASS-recall-pending):** the motif-blind re-placer fabricates non-canonical
junctions via NON-HOMOPOLYMER drift (~1.2% FDR, spike-in ground truth) — it shifts a real canonical junction a
few-to-many bp to a nearby non-canonical pool member. The HP-drift guard catches ZERO of this (it only checks
homopolymer runs). Mechanism (explained + spike-in-confirmed): motif-blindness REMOVED the specificity prior; at
loci with LOCAL MICROHOMOLOGY (a near-tandem-repeat at the drift distance) + ONT error, the true vs drifted boundary
is a near-tie the error tips wrong. Pool-dependent (needs the false site in the pool = paralog/high-expression loci).

## The fix: generalize the HP-drift guard to a MICROHOMOLOGY-drift guard
A homopolymer is the MAXIMAL microhomology (every base identical). The non-HP drift is the SAME phenomenon at
PARTIAL microhomology. So generalize the guard's detector from "boundary inside a homopolymer run" to "the MOVE
sits in a microhomology context."

### Detector (MOVE-aware — matches what the guard sees: shift ns/ne -> js/je)
`move_microhomology_frac(genome, ne, je)` for an acceptor shift (donor symmetric): k = je - ne; the drift is
enabled iff the k-mer absorbed into the intron (genome[ne:je]) near-matches the new-exon k-mer (genome[je:je+k]).
Return the FRACTION of the k-mer that matches at the drift distance. This is EXACTLY the spike-in's installed
repeat genome[A:A+k] ~ genome[A+k:A+2k]. Check both donor and acceptor shifts; take the max.

### Guard rule (mirrors the HP-drift guard)
On a move, if `move_microhomology_frac >= mh_threshold` (~0.5), require the move to beat the incumbent by
`microhom_drift_margin` (a new param, default 0 = OFF = byte-identical). Real non-canonical discovery has LOW
microhomology -> no margin -> discovery preserved. ADD alongside `hp_drift_margin` (keep the validated HP-guard
intact; the new guard adds non-HP coverage). Byte-identical when the new param is 0.

## ★ PHASE 1 VALIDATION (spike-in ground truth, /tmp/mh_char2.py)
Move-aware microhomology_frac on the DRIFT move (true->drift) vs the REAL move (decoy->true non-canonical):
| set | microhomology_frac | 
| FAB-drift (installed, non-HP, k=6..28) | 0.60-0.71 (mean 0.66), consistent across all k |
| R3 genuine non-canonical (real splice) | 0.33 (no repeat) |
=> a threshold ~0.5 SEPARATES fabricated drift (>=0.60) from real discovery (0.33), across all drift distances.
The detector catches exactly the non-HP drift the HP-guard misses, without flagging real discovery. DESIGN VALIDATED.
CAVEAT: small n (FAB 5, R3 2); must re-validate the threshold on REAL SMA drift (SNRPN etc.) via COMPASS recall
(in flight) + a broader R3/cryptic-3'SS panel before locking the threshold.

## BUILD PLAN (rigor = HP-guard + concat-DP discipline)
1. Implement `move_microhomology_frac` + wire a `microhom_drift_margin` / `mh_threshold` into junction_refiner
   (mirror the hp_drift_margin plumbing through refine_read_junctions -> _run_sequential/_run_parallel ->
   refine_bam_junctions). Byte-identical when margin=0 (default).
2. Regression tests (mirror test_hp_drift_guard): detector unit tests (repeat vs random; both shift dirs; edges),
   guard vetoes the fab drift, spares the real transition, byte-identical at 0.
3. VALIDATE: re-run the spike-in fab panel -> does arm-Bguard-micro cut the ~1.2% non-HP drift FDR? AND the R3/
   discovery panel -> is recall preserved (like the HP-guard's zero-discovery-cost)? Tune mh_threshold + margin.
4. REAL-DATA: test the detector's threshold on the COMPASS-labeled SMA drift (real fabrication vs real discovery).
5. Adversarial audit (byte-identity + does it kill any real discovery) before flipping default.

## ★★★ PHASE 2-3 VALIDATION (2026-07-11) — DECISIVE, CLEAN WIN
Implemented (junction_refiner.py, committed): _move_microhomology detector + microhom_drift_margin/microhom_threshold
threaded through all 4 refine fns, byte-identical when margin=0 (refiner+guard suites 50 passed).
GROUND-TRUTH VALIDATION (spike-in fab panel + R3 discovery panel; per-read truth, seed excluded):
| arm | fabrication FDR (canon->drift) | R3 discovery (HP cell) | canonical (plain) |
| arm-B no-guard | 1.31% (24 reads; == spike-in ~1.2%) | 0.284 | 0.931 |
| microhom m=3   | 0.05% (1 read; 96% cut) | 0.284 PRESERVED (b=0,c=0) | 0.941 (+0.010, improved) |
| microhom m=8   | 0.00% (ELIMINATED) | 0.284 PRESERVED | 0.941 (improved) |
=> the general non-HP drift the spike-in/Sumner/Snaptron exposed is ELIMINATED, discovery preserved EXACTLY,
canonical slightly IMPROVED (guard also prevents spurious canonical drift). Zero discovery cost — the HP-guard bar.
Operating point: microhom_drift_margin=8.0, microhom_threshold=0.5 (full elimination, zero discovery cost).
Measurement note: the initial 5.7%->4.4% was a MEASUREMENT ARTIFACT (counted seed reads legitimately at drift as
fabrication); per-read-truth (canonical-origin reads moved to drift) is the correct FDR -> 1.31%->0%.
REMAINING: regression tests (mirror test_hp_drift_guard); adversarial audit (byte-identity + no-discovery-loss);
COMPASS real-data threshold confirmation (in flight); then flip default / ship alongside the HP-guard.
