# FINDING — the v2 re-arbitration flattens genuine non-canonical junctions (measured; two mechanisms)

**Date:** 2026-08-09 (late evening). **Context:** the resolver's v2 junction re-arbitration
(`feat/overhang-resolver-641` @ `de98896`; Chanfreau planning/644b) added a splicing-grammar
tiebreak: a non-canonical-class incumbent junction loses to a canonical-class alternative at
match-or-beat ED. planning/644b's verdict criteria measure the accuracy upside (gold recovery,
junk budget) but NOT the discovery cost — whether genuine non-canonical junctions with a nearby
canonical decoy get flattened. The merged Re-aligner branch carries the by-construction ground
truth for exactly that (the noncanon_sim WT+cryptic mixture panel), so we ran the missing control.

## Setup

Scratch merge of `feat/realigner-triage` (7c0a8f6) + `feat/overhang-resolver-641` (de98896) —
merged CLEAN, combined suites 113 passed. Smoke panel (`build_panel --panel smoke`: R0/R3 ×
plain/HP × decoy-3bp, WT+cryptic mixture, INTRONFREE + canonical controls; 200 reads/cell),
fallback independent error injector (`--error-regime mid`), harness minimap2 (annotation-blind,
`-uf -k14 -G5000 --splice-flank=no`). Three arms on the SAME BAM:
`mm2` (raw) · `noarb` (resolver, `arb_enable=False`) · `v2` (resolver, full re-arbitration).
Driver: `scripts/benchmark/noncanon_sim/rearb_noncanon_control.py` (committed alongside;
requires the resolver branch merged). Outputs: scratchpad `rearb_control_out/`.

## Result

| cell | n | mm2 rec | noarb rec | v2 rec | Δ(v2−mm2) |
|---|---:|---:|---:|---:|---:|
| **R3 · HP · decoy** (cryptic non-canonical) | 200 | 0.735 | 0.735 | **0.675** | **−6.0 pp** |
| **R3 · plain · decoy** | 200 | 0.680 | 0.680 | **0.655** | −2.5 pp |
| R0 · plain · decoy | 200 | 0.880 | 0.880 | 0.860 | −2.0 pp |
| WT / INTRONFREE / canonical controls | 800 | — | — | — | 0 (unchanged; fjFDR ≈ 0) |

- `noarb` ≡ `mm2` on this panel (clip stage resolved 0 reads) ⇒ **the entire drop is the
  re-arbitration pass.** Counters: `arb_njunc_checked=1036`, `arb_shifted=33`, `clean_skip=29`.
- **XB-tag classification of the 33 shifts** (margin acceptance requires `ed_new ≤ ed_cur − 2.0`;
  anything else is the grammar tiebreak):
  - **15 grammar-tiebreak shifts — ALL on cryptic templates**, several at *equal* ED
    (e.g. `9.5→9.5`, `8.0→8.0`): the pure canonical-preference snap, arm-A behaviour.
  - **18 margin-beat shifts**: ~9% ONT-like error makes the canonical decoy genuinely score >2
    better on individual reads. Per-read arbitration CANNOT fix this class (the DP honestly
    favours the decoy on those reads); only POOL-level evidence can (the recurrent true site vs
    per-read-random drift) — the Station-C argument, measured.

## Caveats (honest scope)

Smoke scale (200/cell), fallback error model, and this panel measures ONLY the discovery cost —
the rearb's accuracy upside (canonical alt-SS recovery, D→N merges, the 644 gold gains) is real
and measured separately on the 617 gold set (H2 job `14271360`). Both numbers are needed; this
control supplies the missing one. minimap2's own raw recovery is 0.68–0.74 here — the decoy
drift is the same physics minimap2 suffers; rearb adds to it rather than causing it.

## Recommendations (sent to 641 via inbox)

1. **Split `arb_grammar` out of `arb_enable`** as its own config knob; grammar OFF (or tag-only)
   in discovery mode. Removes the 15-shift tiebreak class (~half the flattening) where the
   mission is non-canonical discovery (upf1Δ/prp18Δ).
2. **Tag grammar-tiebreak acceptances distinctly from margin acceptances** (the XB payload
   already suffices to reconstruct, but an explicit marker is cheaper downstream) — the triage
   layer will treat any grammar-moved read as TRIAGED, never high-confidence, so the flatten
   stays recoverable at Station B/C instead of invisible.
3. **Adopt this control into the resolver's verdict criteria** (it is one script on the merged
   tree) — the FP budget covers junk junctions; this covers junction IDENTITY loss, which junk
   counts cannot see.
4. The 18 margin-beat shifts are the standing argument for the pool-level gate: no per-read
   margin fixes them without also killing genuine rescues (the hold-margin no-sweet-spot result,
   re-observed at the arm level).
