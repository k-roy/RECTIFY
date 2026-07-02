# Real-human-data aligner behavior — Sumner/human RECTIFY evidence (Explore, 2026-07-02)

*Read-only synthesis of RECTIFY's documented real-human aligner behavior, commissioned
to reconcile a synthetic benchmark result (on toy reads mapPacBio recovered the
non-canonical novel junctions minimap2 flattens) against how the panel ACTUALLY
behaves on real human ONT data. Citations are to repo docs; full agent transcript in
the session tasks. The reconciliation (§Bottom line) is the load-bearing result.*

## Per-aligner real-human behavior (cited)

- **mapPacBio (BBMap)** — the synthetic "winner," but CATASTROPHIC on real human ONT-DRS:
  - Krizanović 2018: BBMap spans exon–exon junctions in only **~26.8%** of ONT reads
    (GMAP 87.1%). (`docs/aligners/mapPacBio.md:75-127`)
  - SMA_GSB2394 chr5 ONT-DRS (2026-05-25): **21.7%** of mapPacBio introns are
    mapPacBio-unique, and **97.7%** of those (23,037/23,590) are >5 bp from ANY
    annotated junction — i.e. novel-LOOKING but almost all SPURIOUS. It sole-wins 39%
    of chr5 reads but its junctions are only **77% annotated** vs 97–98% for the
    aligners it beat → **the consensus SELECTS mapPacBio's artifacts, not filters
    them.** (`docs/aligners/mapPacBio.md:86-99`)
  - Param-brittle; crashes on reads >6019 bp (RECTIFY splits); emits multiple
    all-PRIMARY records (unsafe for selection). (`mapPacBio.md:26-71,131-150`)
  - Mechanism: **no splice/GT-AG model** — introns are "any scored gap, D→N relabeled."
- **minimap2** — workhorse (~49% chr5 wins); GT-AG soft prior + `--junc-bonus 9` soft
  annotation bias (novel discoverable but disfavored); ±1–3 bp HP jitter.
  (`minimap2.md:43-44`, `dev/aligner_investigation/SPLICE_JUNCTION_PLACEMENT.md:39-53`)
- **uLTRA** — annotation-guided; excels on tiny exons; **structurally circular on
  novel** — snaps to exact GTF donor/acceptor (zero wobble), reads >10% off-annotation
  fall through to a minimap2 fallback → **cannot beat minimap2 on novel structure**
  ("snapped-to-GTF = correct" is assumed, not tested). (`ultra.md`, `SPLICE_JUNCTION_PLACEMENT.md:71-81`)
- **deSALT** — de-novo (no `-G`), cross-read consensus → homogeneous junctions, wins
  78.9% junction-quality; but chimera-prone + deterministic SIGSEGV/OOM on some
  structures; Linux-only. (`deSALT.md:75-91`, `SPLICE_JUNCTION_PLACEMENT.md:55-69`)
- **GMAP** — REJECTED (2026-06-14): 71% annotated (bar ≥95%), 5 bp anchor (bar ~13),
  STOLE 2655 ≫ ADDED 692 (net harm 3.8×); displaces better minimap2/deSALT alignments,
  spurious novels at mapPacBio-like rates. (`gmap.md`)
- **The 3 COMPASS A549 junctions** (SLC35A4/SQSTM1/TMED9): SLC35A4 = REAL non-canonical
  (168 short + mm2 101 + deSALT 102 + uLTRA 81 all at the NON-canonical coord — genuine
  novel, not a snap); TMED9 REAL; SQSTM1 inconclusive. (`dev/HANDOFF_SHORTREAD_P5.md:11-76`)

## BOTTOM LINE — the synthetic mapPacBio "win" is its real-data PATHOLOGY, not a virtue

mapPacBio recovers the synthetic non-canonical junction for **the exact reason it fails
on real human data**: it has **no splice-model gate** — it emits any scored gap as an
intron. On clean, error-free, single-true-intron synthetic reads that indiscriminate
behavior happens to land on the true site (looks like recovery). On real noisy human
ONT-DRS the same behavior yields a **97.7% spurious-novel rate** and consensus that
selects its artifacts. So:

- **The "pivot" branch (the panel already covers non-canonical discovery via mapPacBio)
  is REFUTED by real data.** mapPacBio's "coverage" is illusory — it emits true AND
  false junctions indiscriminately; its precision on real human data is catastrophic.
- **The panel has NO member that discovers novel non-canonical junctions at acceptable
  PRECISION:** the workhorses (minimap2/uLTRA/deSALT) are precise-but-annotation/GT-AG-
  biased (they flatten novel non-canonical); the de-novo members (mapPacBio/gmap) are
  sensitive-but-imprecise (spurious-novel-dominated, both effectively rejected on human).
  That precision–recall gap on novel non-canonical junctions is EXACTLY what a
  calibrated-empirical-−logP native member targets.
- **=> the native aligner IS justified**, and the synthetic benchmark MUST add an
  FDR/PRECISION dimension: recovery-alone is GAMED by indiscriminate emitters (mapPacBio
  "recovers" by emitting everything). The blindspot ladder needs a FALSE-non-canonical-
  junction control (spurious N-op rate on reads with NO true novel junction, and
  FP-junction rate alongside recovery) — the over-call lesson (C1 insertion-discount /
  the fp_canonical_snap track), applied to the panel-recovery measurement.

## Actions
1. Add an FDR/precision axis to the blindspot harness (per-aligner spurious-non-canonical
   N-op rate + FP-junction rate), so mapPacBio's synthetic "recovery" is scored against
   its false-emission rate. This is the missing control that reconciles synthetic vs real.
2. Interpret the in-flight rgen/rgerr/error runs THROUGH this lens: expect mapPacBio
   recovery to hold on clean but its FP rate to explode under error / real complexity.
3. The native-member case is strengthened, evidence-grounded; proceed to addressability
   formalization + the member design once the FDR-augmented measurement confirms the gap.
