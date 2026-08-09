# C4 — paralog / multi-copy loci (POA-pooled consensus): gate plan + design (2026-06-29)

The paralog facet. C4's premise: near-identical multi-copy loci (SMN1/SMN2-style)
defeat per-read **window/locus selection** — a read maps to the WRONG copy — and a
**POA-pooled** member (cluster reads at the locus → build a consensus → align the
consensus once → project back to the reads) recovers the true copy where a single
read cannot. The headline metric is **locus-selection accuracy**
(`scorer.locus_accuracy`: mapped contig == truth-origin contig).

Mirrors the C1/C2/C3 gate discipline (`dev/C3_DESIGN.md`): establish the incumbent
is BELOW CEILING on *recoverable* reads AND the gap is **pooling-addressable** (NOT
the paralog-zero-evidence trap) BEFORE building any pooling member, with a
pre-committed null, a zero-evidence guard, and a membership-preserving (minority-
collapse) fence. **Fitness = the truth set, NEVER any internal score.**

## The pinned incumbent (firsthand)

C4 is a **WHICH-CONTIG** (locus-selection) question, not a within-contig placement
question, so the single-contig DP arms used by C1/C3 (`align_exon_block_global`) do
**not** apply — there is no "flat vs law" arm for locus selection. The valid
incumbent is the **seed-chain-align aligner's contig choice**: `minimap2 -ax splice
-uf`, scored by `scorer._score_read`'s C4 readout (`locus_correct` /
`locus_incorrect` / `locus_mapq0`: did the read map to its TRUE origin contig?).

minimap2 is a **fair AND conservative** proxy for the full panel here:
- the decisive-SNP logic is **aligner-invariant** — with 2 copies differing only at
  the SNPs a fragment covers, an intact distinguishing base makes the true copy win
  by **exactly 1 mismatch**, so any seed-chain aligner picks it;
- minimap2/gapmm2/uLTRA **share the per-read blind spot** (OVERVIEW), so the panel
  cannot beat minimap2 on the reads where minimap2 is already at ceiling;
- the shipped scorer/smoke locus metric is itself minimap2-based.

Running the full 5-aligner panel would not change the verdict and costs cluster
compute — out of scope for a decoder-free Tier-1 gate.

## The C4-addressable slice (what the stratum was built to probe)

`controlled.gen_paralog_stratum` (3 families × 2 near-identical 600 bp contigs
differing at 3 **spread** distinguishing SNPs) emits two read kinds with
read-of-origin truth:

- **`span`** — full-length read covering ALL spread SNPs (redundant origin
  evidence) → the **at-ceiling control** (proves the locus metric is not trivially
  failing).
- **`weak`** — a 200 bp fragment centered on **exactly ONE** distinguishing SNP →
  per-read the lone SNP is sometimes noise-flipped, so minimap2 drops below ceiling;
  the SNP is recorded as a `VariantTruth` so the gap's recoverability can be tested.

The C4-addressable slice is the WEAK fragment — but **only** the sub-population that
still carries identifiable evidence. A fragment whose lone distinguishing base has
been **destroyed/inverted by noise** is informationally identical to (or points at)
the wrong copy → unrecoverable per-read by ANY method, incl. C4 → the **paralog-
zero-evidence trap** (the OVERVIEW's named null). The gate must separate these.

## Gate step 0 (RUN FIRST) — window-selection headroom (DECODER-FREE): `scripts/benchmark/c4_headroom.py`

The gate-deciding measurement, with NO pooling member (so zero construction-tuned-
win risk). Over the paralog WEAK fragments, scored vs TRUTH, classify each read on
two axes and compute `ceiling − incumbent`:

    evidence : intact      = read still carries THIS copy's SNP allele  (IDENTIFIABLE)
               looks_other = base flipped to the OTHER copy's allele    (wrong-evidence)
               third       = base flipped to a non-paralog base         (zero-evidence)
    placement: did minimap2 map it to its TRUE origin contig? correct / WRONG / unplaced

    ceiling   = the IDENTIFIABLE universe = intact-base reads
    incumbent = freq(minimap2 maps an intact-base read to the correct copy)
    HEADROOM  = freq(intact base BUT mis-placed) = the genuine C4-addressable gap

**Pre-committed null:** if HEADROOM ≈ 0 → the incumbent is at ceiling on every
identifiable read and the entire below-ceiling gap is corrupted-base (zero-evidence)
reads → **C4-as-window-selection REFUTED** on this stratum; document and stop. Do
NOT build a pooling member to chase a gap that is the zero-evidence null.

The probe ALSO audits the smoke **(F)** "pooling recovers 6/6 pools" addressability
claim two ways (see Verdict).

## RESULT — step 0 (2026-06-29, reps=100/copy/kind, families=3, noise=0.05, minimap2 -ax splice -uf)

```
SPAN control (all spread SNPs):                 locus_acc=0.9967  (598/600)   -> at ceiling

WEAK fragments (one covered SNP), n=600
  evidence       correct   WRONG  unplaced     n
  intact             574       0         0    574
  looks_other          0       9         0      9
  third                9       8         0     17

ceiling  (IDENTIFIABLE = intact-base reads)                 = 574
incumbent(minimap2 correct | intact base, placed)           = 1.0000 (574/574)
HEADROOM (intact base BUT mis-placed = genuine C4 gap)      = 0.0000 (mapq0=0, mapq>0=0)
zero-evidence null (looks_other 9 + third 17)               = 26
```

Robustness (probe sweep, reps=100): noise 0.05 → HEADROOM 0/574; noise 0.10 →
0/550 (3 intact reads UNPLACED — alignment sensitivity, not mis-cluster); noise
0.15 → 1/347 (and 184/531 intact reads UNPLACED — the failure mode has **migrated to
alignment sensitivity / C5**, not C4 locus selection). The genuine identifiable-but-
misplaced population is **~0 at every realistic noise level**. (Result file:
`scripts/benchmark/c4_headroom_result.txt`.)

## VERDICT — C4-as-window-selection DEFERRED (incumbent at ceiling on this stratum; the stratum is near-zero-power for C4's target regime)

Scoped precisely: on the shipped 2-copy/1-covered-SNP paralog stratum the incumbent
is at ceiling on every **identifiable** read and the entire below-ceiling gap is the
pre-committed zero-evidence null — so **no pooling member is justified now** — but
this stratum **cannot test** C4's actual target regime (≥3 near-identical copies with
multi-PSV linkage, where a single covered SNP leaves residual mapq0 ambiguity *with*
intact evidence). The honest facet-level label is therefore **DEFER, not a hard
REFUTE**: the measurement is valid but near-zero-power for the phenomenon C4 names —
on a 2-copy/1-SNP construct `identifiable ⟺ minimap2-correct`, so HEADROOM=0 is
near-structural, not an empirical refutation of C4's real regime. Do NOT tally C4 as
a "refuted" facet; it is **deferred-with-reason**.

Three independent nails, all decoder-free, all scored vs truth:

1. **Incumbent at ceiling on identifiable reads.** Every weak fragment that still
   carries its copy's distinguishing base (574/574 at noise 0.05) is placed on the
   correct copy by minimap2. HEADROOM = 0.000; `genuine_mapq0 = 0` (the incumbent
   isn't even *flagging* these ambiguous — it resolves them).

2. **The entire below-ceiling gap is zero-evidence.** Every mis-placed weak read is
   a corrupted-base read: `looks_other` (lone SNP flipped to the OTHER copy's base →
   the read genuinely IS the other copy across every position the fragment covers →
   wrong-evidence, unrecoverable per-read) or `third` (flipped to a non-paralog base
   → equidistant to both copies → mapq0 coin-flip). This is the pre-committed
   **paralog-zero-evidence NULL** and it fully reconciles smoke (F)'s own
   `weak_acc ≈ 0.94`: the `1 − weak_acc` deficit falls **100% on corrupted-base
   reads, 0% on intact-base reads**.

3. **The smoke (F) "pooling recovers 6/6" addressability proof does not hold up.**
   - **Truth-circular:** (F) forms its pools by `t.true_transcript` — the TRUE copy
     label — so it *presupposes* the very read→copy clustering a real C4 member must
     SOLVE. "Majority base within the true-copy pool recovers the allele" is
     necessary but not sufficient for "pooling recovers the copy assignment."
   - **Redundant:** minimap2's OWN per-read placement already recovers each pool
     (95–100/100 reads to the correct contig). Clustering by the incumbent's own
     placement — the realistic, truth-free move — yields the same pools (6/6).
     Pooling adds nothing the incumbent didn't already deliver per-read.
   - **Minority-collapse untested:** (F)'s pools are equal-reps-per-copy, so a
     "majority base" can never collapse a minority copy. A real pooling member faces
     both the clustering problem and abundance-skew minority-collapse (the
     membership-preserving trap) — neither exercised by (F).

C4 is **deferred-with-reason** (distinct from the C2/C3 hard refutes: those incumbents
were at ceiling on strata that *could* express the targeted gap; this stratum cannot).
The gate's two solid, transferable contributions are (1) **no pooling member is
justified now** and (2) the **(F) addressability-proof takedown** — most plausible
ideas do not survive contact with truth, and catching this before writing
pooling/clustering code is the prove-don't-assert discipline working as intended.

**Do NOT build the POA pooling member.** There is no identifiable-but-misplaced
population on this stratum for it to recover, and the only "recovery" smoke (F)
credits is truth-circular and already delivered by the incumbent per-read.

## Honest residue / scope (the named revisit locus)

- **The 2-copy / 1-covered-SNP construct *structurally cannot* produce the regime
  C4 targets.** With 2 copies and a single covered distinguishing SNP, an intact
  base is decisive by exactly 1 mismatch → mapq > 0 always → no identifiable-but-
  misplaced cell can exist. The genuine SMN1/SMN2 hard case is **≥3 near-identical
  copies** (or a paralog family) where a single covered SNP only *excludes one* copy
  → residual ambiguity (mapq0) **even with intact evidence**, exactly where pooling/
  read-linkage across multiple informative sites could earn its keep. This stratum
  does not exercise that regime, and `genuine_mapq0 = 0` confirms it.
- **Do NOT synthesize a ≥3-copy construct now to manufacture a gap** — that is the
  C2-`interrupted` construction-tuned-win trap (tuning the corpus until the facet
  "wins"). C4 becomes live again ONLY when a **measured real-data trigger** (e.g. a
  real SMN1/SMN2 ONT-DRS locus showing a mapq0-with-evidence read population the
  panel mis-clusters) demonstrates the regime exists outside a hand-built corpus.
  Until then: defer, do not build. This mirrors C3's "named revisit locus, do not
  build until built-and-shown."
- **Alignment-sensitivity ≠ window-selection.** At high noise the dominant weak
  failure is UNPLACED (alignment sensitivity), which is **C5 / FracMinHash-fallback**
  territory (reads with no acceptable window at all), not C4 locus selection. Even
  where the incumbent fails on paralogs, C4 is not the lever.

## Recommendation to the director (flagged separately — NOT applied here)

Smoke **(F)** currently PASSES and prints "provably C4-ADDRESSABLE." This gate shows
that addressability pass rests on a truth-circular + redundant pooling check.
Consider, as a follow-up (do not weaken a green fence unilaterally mid-stream):
- **downgrade (F)'s message** from "provably C4-ADDRESSABLE" to "MEASURED +
  DISCRIMINATING; addressability is the zero-evidence NULL on this 2-copy stratum
  (see `dev/C4_DESIGN.md`)"; and/or
- **add a truth-free clustering check**: cluster the weak reads by minimap2's own
  placement (NOT by `true_transcript`), then test recovery + minority-preservation —
  which would either surface a real gap or confirm the refute inside the fence.

## Files (this gate)
- `scripts/benchmark/c4_headroom.py` (+ `c4_headroom_result.txt`) — window-selection
  headroom: the decoder-free intact/corrupted × correct/WRONG 2×2, the (F) pooling
  audit, and the pre-committed verdict.
- REUSES `controlled.gen_paralog_stratum` unchanged (no edit to `controlled.py`).
- M1-light, decoder-free, fitness=truth. NO product/member code added (the C1-idiom
  refute outcome). Smoke gate GREEN (unchanged).
