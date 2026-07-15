# Re-placer emits degenerate N-flanking compensating I/D CIGARs on multi-junction reads

**Date:** 2026-07-15 · Sumner SMA DRS (SMA_191), panel_deep re-placer BAMs · render_single.py / single-read IGV inspection
**Status:** CONFIRMED on 4 loci, ~85–95% of re-placer reads affected. Reframes the "fabrication/drift" estimate.

## What we found

Inspecting single junction-spanning reads (raw minimap2 vs RECTIFY re-placer) at the 4 microhomology-drift
leads, the re-placer's "acceptor shift" is **largely an artifact of the CIGAR it emits**, not a real relocation
of the read's exon sequence.

Example — read `a53fcd26…` at PCBP2 (chr12), the 2805 bp intron (donor≈53468808 region):

```
raw mm2 :  … 933N  56M  2805N        170M …     intron ref[53468832,53471637), exon2 M starts 53471637 (canonical)
re-placer: … 933N  17M  39I  2805N  39D  131M … intron ref[53468793,53471598), THEN 39D back to 53471637, exon2 there
```

The re-placer slid **both** intron boundaries 39 bp upstream and pinned the read's exon sequence back in place
with a **39 bp insertion on the donor side + 39 bp deletion on the acceptor side**. Net: the exon2 sequence
still lands at **53471637 — identical to raw** — but the reported N-op boundary (the junction *coordinate*) is
now 53471598, 39 bp off, and the alignment carries +78 bp of uncompensated indel penalty. It does this at
**multiple junctions in the same multi-exon read** (`68M 7I 2192N 7D 93M`, `17M 39I 2805N 39D 131M`,
`131M 39I 7598N 39D 262M`).

## Two metrics, two stories

Measuring the "acceptor" as the **N-op end** (what §4b / the drift panel did) vs. the **first query-aligned ref
base after the intron** (where the exon sequence actually sits, skipping leading D):

| Locus | N-op-end shift (WRONG) | True-acceptor shift (CORRECT) |
|---|---|---|
| PCBP2   | 54% moved ~40 bp | **4% moved (≤12 bp)** — exon stays at canonical 53471637 |
| UBA1    | 8% moved | **3% moved** — clean, stays at well-supported non-canonical |
| CACNA2D3| 27% moved | **15% moved (small)** |
| SNRPN   | 68% moved +6 bp | **64% moved +6 bp** — GENUINE, in-window |

The "PCBP2 drifts 40 bp beyond W=28" concern **evaporates** under the correct metric. The only genuine
relocation among the leads is SNRPN's +6 bp (in-window, well within W=28).

## Systematic prevalence (SMA_191 re-placer, reads with an N-flanking I+D pair)

| Locus | reads | with N-flanking I+D | total I/D pairs | modal I length |
|---|---|---|---|---|
| SNRPN    | 978 | 850 (87%) | 2049 | 26 |
| PCBP2    | 858 | 818 (95%) | 2231 | 36 / 48 / 42 |
| UBA1     | 378 | 344 (91%) | 2487 | 40 |
| CACNA2D3 |  27 |  23 (85%) |  194 | 19 |

The insertion lengths cluster exactly at the apparent "drift" distances — confirming the I/D pairs ARE the
mechanism generating the phantom drift signal.

## Why this matters

1. **The §4b "27% fabrication" number is likely dominated by this artifact** if it measured N-op coordinates
   (it classified drift by donor/acceptor position). The genuine cryptic-drift fraction is much smaller
   (SNRPN-style small in-window moves), which weakens the case that the microhomology guard is catching real
   fabrication. RE-RUN §4b classification on the TRUE acceptor (first query-aligned base), not the N boundary.
2. **The re-placer has a CIGAR-reconstruction defect**, separate from scoring/fabrication: it prefers/emits
   `NbpI … intron … NbpD` over a clean `N … M`, adding +2·Nbp indel penalty for zero benefit. No sane scorer
   should choose this — so either (a) the scorer never sees the global indel cost (scores the junction locally),
   or (b) the flank-stitching in bam_writer / CIGAR surgery leaves compensating I/D when it moves a boundary.
   This misreports every affected junction coordinate.
3. **The W-horizon worry is largely moot for these real loci** — genuine moves are ≤~6–12 bp (in-window). The
   scorer-level W-fix is lower priority than the CIGAR bug.

## Reproduce

```
cd /Users/kevinroy/sumner_igv_slices
python render_single.py          # single moved read per locus, both placements, N/D labelled
# per-read CIGAR dump: the a53fcd26 PCBP2 example above
```
Figures: `single_{PCBP2,SNRPN,UBA1,CACNA2D3}_SMA_191.png`. The PCBP2 panel shows the N2805 + black-D-bar
mechanism directly (re-placer exon2 lands at the same x as raw, with a 39 bp deletion bridging the slid N).

## FIX + §4b RE-RUN RESULT (2026-07-15)

Fixed in `_apply_junction_replacement` (commit e40ca00): refuse a BOTH-boundary move that RAISES the
read's I+D burden. Genome-wide over the 13-sample panel, **~94% of the refiner's junction moves were the
phantom** (reverted to raw); ~5.5% legit single-boundary, ~0.4% clean microhomology slide.

§4b re-run (recall_analyze_v2.py `--revealed-from-lrdir`, matched filter OLD lr_junc vs FIXED
lr_junc_fixed). **The fix nearly eliminates the fabrication while keeping — and PURIFYING — the genuine
discovery signal:**

| §4b (≥2 samples)          | OLD (buggy) | FIXED   | change |
|---------------------------|-------------|---------|--------|
| revealed non-canonical    | 170,066     | 23,277  | −86%   |
| drift (fabrication)       | 39,933      | 3,488   | **−91%** |
| recovery (genuine, SR-exact)| 8,101     | 1,462   | −82% (mostly phantom coincidences removed) |
| drift RATE                | 23.5%       | 15.0%   |        |

| §4b (≥5 samples, high-conf) | OLD    | FIXED  |        |
|-----------------------------|--------|--------|--------|
| revealed non-canonical      | 23,235 | 4,983  | −79%   |
| recovery RATE               | 5.7%   | **10.8%** | ~2× precision |
| drift RATE                  | 18.6%  | **9.1%**  | ~½     |

§3 recall: arm-B recall **17.5% → 4.1%** (the inflated 17.5% was phantom); still ~8× over raw-mm2 (0.5%).
SR-canonical control unchanged (0.946 both) — the fix does not touch canonical placement.

Interpretation: arm-B's apparent non-canonical "discovery" was ~80–90% CIGAR fabrication (phantom N-op
moves the read never supported). The fix removes it; the surviving reveals are ~2× more likely to be
real recoveries and ~½ as likely to be drift. A residual ~9–15% drift remains (single-boundary /
clean-slide moves to a canonical neighbor) — that is the genuine target for the microhomology guard /
positional close, now a small-set problem rather than the dominant signal.

NOTE: OLD_matched (170,066 revealed, 23.5% drift) differs from the published recall_result_fast.txt
(142,498, 27.3%) because the matched re-run rebuilds the revealed set from lr_junc *.refined.junc.tsv
(canonical==0) rather than the original panel_deep/*.revealed_noncanon.tsv generator. OLD_matched is the
correct apples-to-apples baseline for the fix; the published number used a different generator.

Result files (Sherlock): `/scratch/users/kevinroy/sma_recall/recall_result_{OLD,FIXED}_matched.txt`.

## Done / Next

- [x] Located the I/D-emitting code path: `_apply_junction_replacement` general path (junction_refiner.py).
- [x] Fixed (commit e40ca00) + regression tests; §4b re-run confirms −91% fabrication.
- [ ] (Optional) Full through-the-refiner re-run for publication (reconstruct the panel driver).
- [ ] Re-evaluate the microhomology guard / positional close against the residual 9–15% drift (small set).
