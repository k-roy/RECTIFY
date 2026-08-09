# C3 — calibrated LLR arbitration: gate plan + design (DRAFT 2026-06-29)

The KEYSTONE facet. C1 (and any future soft CPA/discovery signal) only improves
the **pipeline** if the consensus can compare aligner paths by **likelihood**
instead of integer-max. C3 = the consensus picks the winning aligner per read by
a calibrated **−log P likelihood ratio (LLR)**, not by the shipped integer/edit
scalar. It is the named structural fix for the **0.09→1.07 artifact** (a pure
re-weighting of the internal score once flipped an aligner's apparent quality
with NO alignment change).

Mirrors the C1/C2 gate discipline (`dev/C1_DESIGN.md`, `dev/C2_DESIGN.md`):
establish the incumbent ARBITER is BELOW CEILING on recoverable reads AND the gap
is arbitration-addressable (a coherent LLR re-ranks toward truth, not generic
noise) BEFORE building any decoder, with a pre-committed null, a zero-evidence
guard, and an over-call/safety control. Fitness = the truth set, NEVER any
internal score.

## The two incumbent arbiters (BOTH pinned firsthand — the C2 "real incumbent" lesson)

There are **two** integer-max arbiters in the consensus path; C3 must beat the
one that actually decides each stratum (the C2 trap: raw-mm2 looked below ceiling
but the shipped *walkback* was at ceiling).

1. **Alignment-selection** — `rectify/core/consensus/select.py::select_best_alignment`
   picks `max(junction_score)` (from `scoring.py::score_alignment`) with
   tiebreakers (3'-agreement, `_n_annotated`, `canonical_count`). `junction_score`
   is a **flat, quality-blind** sum of integer penalties: −2 per unexplained 5'
   clip base, −2 per non-poly-A 3' terminal error, −1 per junction-proximity error
   (HP weighted 0.5). It is **blind to in-exon indel placement** (scores only clips
   + ±5 bp of junctions). **The 0.09→1.07 re-weighting bites HERE** — a change to
   those −2/−1 weights or the annotated/canonical tiebreakers flips which aligner
   sits at `max_score` with no alignment change.

2. **Corrected-consensus winner pick** — `corrected_consensus.py::merge_corrected_tsvs`
   (the `use_hp_ed` path, the one that actually selects the winning aligner for the
   corrected BAM) sorts `[_effective_chimera_ok ASC, hp_edit_distance ASC, _span
   DESC]` and takes `.first()` per read. `hp_edit_distance` (`_cigar_hp_edit_distance`)
   is **HP-aware on indels** (uses `penalty_table.del_cost/ins_cost`) but **flat on
   mismatch** (X = 1.0) and clips (S/H = 1.0/base), introns free. Two calibration
   defects vs a coherent LLR:
   - **flat mismatch**: an X always costs 1.0 regardless of context;
   - **`del_cost`/`ins_cost` return the `penalty_score` column, which C1 proved is
     NOT −log P** (a reciprocal-rate heuristic, `c/rate_mean`, **incoherent to sum**
     in an additive scalar). So summing it along a CIGAR mis-orders alignments with
     different event *counts/types* relative to a coherent −log P.

The C3-addressable gap, if any, is therefore **substitution-vs-indel / multi-event
arbitration on one coherent −log P scale**, where the flat-mismatch + mis-calibrated
(`penalty_score`) scalar ties or mis-orders members but a calibrated LLR matches
truth. The pure in-run HP-placement case is likely already at ceiling under
`hp_edit_distance` (HP-aware on indels) — the **C2-redux null**.

## Gate step 0 (RUN FIRST) — arbitration headroom (LLR-FREE): `scripts/benchmark/c3_headroom.py`

The gate-deciding measurement, with NO LLR (so zero construction-tuned-win risk —
the trap that sank C2's `interrupted`). Over the controlled corpus with a real
multi-member panel scored vs TRUTH (ambiguity-aware, the `_score_read`
position-exact rule), per stratum:

    arbitration headroom = ceiling − arbiter_accuracy
      ceiling          = freq( >=1 member is position-exact-correct )      (recoverable universe)
      arbiter_accuracy = freq( shipped argmin(hp_edit_distance)+span picks a correct member )
      headroom split:  mis-order  = strict argmin wrong, a STRICTLY-higher-ed member correct
                                    (the coherence defect — the genuine C3 target)
                       tie_arbitrary = hp_ed TIES a correct+wrong member; the span
                                    tiebreak landed wrong (recoverable by ANY tiebreak)

Panel members (controlled, truth-scored disagreement): **flat** =
`align_exon_block_global(penalty_table=None)` (the flat-affine family / C1 matched
baseline), **law** = `align_exon_block_global(penalty_table=LAW)` (the C1 native
member), **mm2** = real `minimap2 -ax splice` (an INDEPENDENT external view whose
placements are produced by no table — the hold-out that kills construction-tuned
wins). Scoped to the **indel strata (HP, HP_HARD, STR)** where the exon-block DP
arm is the valid tool; junction/paralog strata need splice aligners (future
`--with-mm2`-only extension).

**Pre-committed null:** if HEADROOM ≈ 0 → the shipped arbiter is already at ceiling
on recoverable reads → **C3-as-accuracy is REFUTED** on this corpus (the C2
outcome); document and stop — do NOT build an LLR to chase a phantom gap. The
artifact-replay (below) survives as a CI fence only. With three directions already
refuted/deferred, a refute here is a clean, expected outcome of prove-don't-assert.

## The coherent LLR (the C3 analog of C1's "use rate_mean, NOT penalty_score")

Only built if step 0 shows headroom. The single biggest fake-win/fake-null risk is
a mis-specified likelihood; two non-negotiables (advisor):

- **Full-query emission.** Every query base carries an emission cost, INCLUDING
  clipped/unexplained bases — subsume the −2/clip as a background-emission cost. A
  free soft-clip makes the clipping alignment always win → the LLR is garbage.
- **Coherent −log P over M/X/D/I** from the C1 `rate_mean` table (calibrated,
  sums to ~1 per (base_class, hp) context), NOT `penalty_score`. Score both
  hypotheses over comparable event sets (same query). **Unit-test against
  brute-force enumeration on a toy to 1e-6** (SPEC line 263). Note: the table
  currently loads `rate_mean` for D/I only; M/X emission must be sourced coherently
  (derive P(M) = 1 − P(X) − P(D) − P(I) per context, or load M/X rows) — a build
  prerequisite, not an afterthought.

## The matched-arm ablation (fitness — carries the entire gate)

Arm **INTEGER** (the shipped arbiter: `argmin hp_edit_distance` + span) vs arm
**LLR** (the coherent −log P arbiter), each picks a winner per read; score the
WINNER'S alignment vs TRUTH (position-exact concordance) on the **TEST split**.
Exit: LLR concordance ≥ INTEGER, strictly > on the addressable (mis-order) cells.

- **Zero-evidence guard (the C1/G analog):** a headroom read counts only if the
  truth-correct member is STRICTLY preferred by the calibrated likelihood. If both
  candidates are equiprobable under the calibrated model, no method can recover it
  (the iron-triangle / paralog-zero-evidence trap) — exclude/pre-commit as null.
  WATCH the sub-vs-indel case: if truth is a substitution but the context (long HP
  run) genuinely makes a deletion more likely, BOTH hp_ed AND a coherent LLR prefer
  the deletion — that is zero-evidence, not a C3 win.
- **Safety / inertness control:** on the clean/easy subset (members agree, or
  INTEGER already uniquely picks truth) the LLR must NOT flip to a worse pick (the
  C1 insertion-discount lesson — a permissive rule that rewrites correct calls).
- **Inertness fence:** LLR flag-OFF ⇒ byte-identical winner selection vs shipped
  (SPEC 276; the C3 backward-compat fence).

## The artifact-replay (REGRESSION FENCE, not fitness — advisor)

At FIXED placements (a frozen set of member alignments), re-weight the integer
penalties and show the INTEGER winner FLIPS while the LLR winner is INVARIANT —
replicating the named 0.09→1.07. This is **near-tautological** (the LLR doesn't
consume the −2/−1 weights, so invariance is by construction); it proves the
integer path is fragile and the LLR doesn't inherit *that* fragility, but says
NOTHING about whether the LLR is *correct*. File under CI regression-guards (SPEC
265/282); do NOT let its easy green inflate the verdict. The truth-accuracy
ablation carries the gate.

## Flat-Q scope caveat (SPEC:225)

The sim reads have FLAT quality (`'I'*len`). So the LLR's only lever over
integer-max on the current corpus is **error-TYPE/context −log P** (the C1 emission
table), NOT per-base Q. A C3 win here is scoped to **error-type arbitration**, not
"defeating the quality-blind panel" (the panel's headline orthogonality axis).
Injecting correlated Q (the hot-read model already exists) to test the
quality-blind axis is a separate, later validation.

## RESULT — step 0 (2026-06-29, reps=20, ALL splits, 1600 reads, panel=[flat,law,mm2])

| stratum | n | ceiling | arbiter | HEADROOM | mis-order | tie_arb |
| --- | --- | --- | --- | --- | --- | --- |
| HP | 960 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| HP_HARD | 400 | 0.998 | 0.998 | 0.000 | 0.000 | 0.000 |
| STR | 240 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| **TOTAL** | 1600 | 0.999 | 0.999 | **0.000** | 0.000 | 0.000 |

On the indel strata the shipped `hp_edit_distance` arbiter is **AT CEILING on
recoverable reads** (ceiling == arbiter == 0.999; headroom 0; zero mis-order, zero
tie-arbitrary). The C1 law member places truth-correct and the HP-aware `hp_ed`
already selects it — **C1's win flows through the existing arbiter with no LLR
needed.** This is the C2 saturation lesson one facet over: the indel-arbitration
gap C3 would target is already shipped. (Result file:
`scripts/benchmark/c3_headroom_result.txt`.)

## RESULT — step 0b, JUNCTION arbitration (2026-06-29, reps=30, 210 reads)

The 0.09→1.07 locus, tested on M1 (NOT a cluster panel — the SPEC line-127 "fixed
placements at fixed loci" case): members = {`mm2` real minimap2 (snaps), `truth`
the corpus true_cigar (orthogonal placer)}, governing arbiter = `hp_edit_distance`.
Governance confirmed from the wiring: the canonical **correct-first** pipeline
(`run/single_sample.py`, `split_command.py`) emits the `merge_corrected_tsvs`
(`hp_ed`) winner's full record — junctions included; `select_best_alignment`
(`junction_score`) is the older raw-consensus path.

| stratum | n | ceiling | arbiter | HEADROOM | mm2snap | snap&tru | pickSnap |
| --- | --- | --- | --- | --- | --- | --- | --- |
| JUNCTION_DISCOVERY | 120 | 1.000 | 1.000 | 0.000 | **0.467** | 0.467 | **0.000** |
| ANNOTATED / NIC / JUNCTION_AMB | 90 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| **TOTAL** | 210 | 1.000 | 1.000 | **0.000** | 0.267 | 0.267 | 0.000 |

minimap2 **snaps non-canonical→canonical on 46.7%** of JUNCTION_DISCOVERY reads
(the bias is real and live), but **given a truth-site member the `hp_ed` arbiter
picks the snap over truth 0.000 of the time**: introns are free (op 3 → 0 cost) and
the snap induces flanking-exon mismatches that raise the snap's `hp_ed`, so the
arbiter already prefers truth. (Result file:
`scripts/benchmark/c3_junction_headroom_result.txt`.)

## VERDICT — C3-as-accuracy FULLY REFUTED (both arbiters, both strata families)

- **Indel arbitration (`hp_ed`): REFUTED** — at ceiling even on `boundary_sub`
  where the 3 members disagree on 100% of reads (`hr|dis`=0).
- **Junction arbitration (`hp_ed`, the governing arbiter): REFUTED** — at ceiling
  even where minimap2 snaps 47% of reads; given a truth member it never picks the
  snap. The canonical-snap bias is NOT a calibrated-LLR-arbitration gap: where it
  bites the real panel (ALL members herd/snap → no truth member exists) it is
  **C5/discovery** territory (arbitration can't pick what no member produced), and
  the residual `canonical_count` TIEBREAKER in `select_best_alignment` is a
  **tiebreaker-reweight / Discovery-facet** concern — a one-line fix, not an LLR.
- **Falsifies part of C3's own premise:** the keystone justification was "C1 only
  improves the pipeline if C3 exists to consume probabilities." The probe shows
  C1's placement ALREADY flows through `hp_edit_distance` (which consumes the
  penalty table) with `ceiling == arbiter == 0.999` — **no LLR consumer needed for
  C1 to land.** The "C3 is the consumer that makes C1 useful" argument is
  empirically weakened.
- **Do NOT build the LLR.** The artifact-replay (re-weight → integer flips, LLR
  invariant) is near-tautological (the LLR doesn't consume the integer weights);
  it only earns its place if the LLR ships for an ACCURACY reason. No accuracy
  reason → no LLR → the replay fence has nothing to guard. C3 joins C2 as a
  gate-refuted facet (now 4 directions refuted/deferred, 1 confirmed).

## Two structural closures (prove-don't-assert, not new probes)

**Multi-event coherence — the one place a real gap could hide (RESOLVED).**
`hp_ed` scores indels with `penalty_score` (≈ `c/rate_mean`, reciprocal-rate), NOT
−logP — and a reciprocal-rate sum can rank MULTI-event paths differently from a
coherent −logP sum (the only mechanism by which an LLR could beat `hp_ed` on
accuracy). `scripts/benchmark/c3_multievent_check.py` confirms: **(a) the inversions
EXIST** — 89 of the 2-deletion decomposition pairs rank oppositely under
`penalty_score`-sum vs −logP-sum (confirming C1's "incoherent to sum"). **(b) but
they are UNREACHABLE as a truth-favoring tie on a realistic read:** for two
edit-distance-tied alignments of ONE read to assign deletions to DIFFERENT-length
runs, you'd need a cross-run-length-reassignable multi-deletion tie — which neither
same-base HP ambiguity (within-run shifts are ambiguity-equivalent → SAME run length
→ no inversion) nor different-base adjacent runs (each deletion sits unambiguously in
its own run) produces. The inversions are between hypothetical decompositions that
don't correspond to competing placements of a real read; margins are also small
(0.1–1.3 nat). Demonstrating any reachable gap would require a dedicated
**adjacent/interleaved-runs multi-event stratum** (not in the corpus). Absent that,
the refute stands — and this is the *precise, named* locus for any future C3 revisit.

**The `junction_score` path (the literal 0.09→1.07 home) is closed structurally.**
The named artifact was a GMAP *annotated* win/loss flip → `select_best_alignment`'s
`_n_annotated`/`canonical_count` TIEBREAKERS (which `hp_ed` doesn't have). But
`score_alignment` penalizes a snap on the **PRIMARY** score via
`_count_junction_proximity_errors` (−1 per mismatch within 5 bp of a junction); a
snap is ≥1 bp off → carries proximity errors → loses on score, so the
canonical/annotated tiebreakers only bite on exact *ties* a snap won't produce.
Spot-checked: `junction_proximity_errors` fires on 84/360 junction reads in the
corpus. So even in the `junction_score` path the snap is disfavored on the primary
score; any residual canonical-tiebreaker bias is a **one-line tiebreaker reweight**
(owned by the Discovery facet / WS-3), not a calibrated LLR.

### Honest residue / scope (not claimed)
- **Flat-Q (SPEC:225):** the per-base-Q arbitration axis (the panel's headline
  "quality-blind" orthogonality) is UNTESTED — sim reads have flat Q. The
  error-TYPE axis (the table's signal) shows no gap; whether correlated Q opens an
  arbitration gap the table doesn't is a separate, later check (the hot-read model
  exists). Unlikely to flip the verdict but honestly out of scope here.
- **The named revisit locus = the multi-event/adjacent-runs stratum** above. C3
  becomes live again ONLY if that stratum is built AND `hp_ed` is shown to mis-order
  a truth-vs-wrong tie there that a coherent −logP recovers. Until then: do not build.
- **Deferred architectural note (weakened, not dead):** a soft CPA/discovery
  posterior with a continuous-LR consumer remains a *possible* future need, but the
  measured indel finding removes its strongest motivation; gate it on a measured
  consumer, not on the C3 framing.

## Files (this gate)
- `scripts/benchmark/c3_headroom.py` (+ `c3_headroom_result.txt`) — indel-arbitration headroom.
- `scripts/benchmark/c3_junction_headroom.py` (+ `c3_junction_headroom_result.txt`) — junction-arbitration headroom (0.09→1.07 locus).
- `scripts/benchmark/c3_multievent_check.py` — the penalty_score-vs-−logP multi-event coherence/reachability check.
- All three are M1-light, LLR-free, fitness=truth. NO product/member code added (the C1-idiom refute outcome). Smoke gate GREEN.
