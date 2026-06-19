# Improving the Ensemble / Consensus / Winner-Selection Layer

> Verified against `origin/drs-validation-rebuild` @ 366c885 (2026-06-19).

**Scope.** Concrete, prioritized improvements to RECTIFY's per-read aligner-selection
(`merge_corrected_tsvs` in `rectify/core/consensus/corrected_consensus.py`) and the
deprecated raw-BAM consensus path (`scoring.py` / `select.py`). The selection layer in
production is the **HP-edit-distance path ("Path A")**: both merge call sites pass
per-aligner BAMs (corrected when available, raw otherwise) plus the genome, so
`use_hp_ed` is `True` and reads are ranked by a genome-anchored HP-aware edit distance
on the corrected CIGAR — not by a popularity vote. The legacy 5-level sort that scores
3'-end agreement via `_n_agree` is the **fallback only**, used when no BAMs reach the
merge. The remaining open levers are therefore *charging N-ops a calibrated cost*
(they are currently free) and *de-biasing / weighting the live HP-ED selector*, not
"turning Path A on." The win rates (78.9/18.2/2/0.8/0.1) remain single-dataset and
un-committed; they reflect the HP-ED selector, not the legacy sort.

**Convention.** **ESTABLISHED** = standard, published, or already-half-built in the
repo. **NOVEL** = a new design proposed here. Each item: *mechanism → bias it fixes →
feasibility (exact files) → expected impact → validation/ablation → risk.*

Code facts verified this pass:
- Path switch: `corrected_consensus.py:1262` `use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)`.
- Path A key: `[_effective_chimera_ok, hp_edit_distance, _span]`.
- Path B (fallback) key: `[_five_rescued, _chimera_ok, _conf_rank, _n_agree, _span, _n_junc]`.
- `_n_agree`: count of rows sharing `(read_id, corrected_3prime)` — popularity; used in
  the fallback path only.
- HP-ED cost model (`_cigar_hp_edit_distance`): D = HP-aware del_cost (fallback 1.0),
  I = HP-aware ins_cost (fallback 1.25), soft-clip = 1.0/base, **N (intron) = 0 — a free
  pass**. A junction-anchor gate (`_cigar_min_junction_anchor`) can veto under-anchored
  N-ops, but it defaults to `min_junction_anchor_bp=0` (off; the human profile uses 10 bp,
  yeast 0). So on yeast a spurious intron is currently under-penalized — the live open lever.
- Both call sites wire Path A: the single-sample path and `split_command.py` pass
  per-aligner BAMs + genome into `merge_corrected_tsvs`. When only raw BAMs are available,
  a lazy raw-BAM HP path computes HP-ED in memory rather than from pre-written corrected BAMs.
- `_run_correction_per_aligner` writes a corrected BAM `{stem}.rectified_corrected_3end.bam`
  into each `per_aligner_corrected/{aligner}/` dir; these feed the corrected-BAM HP path
  when present.

---

## TL;DR ranking (impact / effort)

| Rank | Item | Type | Effort | Impact |
|---|---|---|---|---|
| 1 | **§5 Win-rate harness + provenance** (do FIRST — unblocks every other claim) | ESTABLISHED | XS | Foundational |
| 2 | **§1 Calibrate the live HP-ED selector** (charge N-ops a cost; add a 3'-end term) | ESTABLISHED | S | High |
| 3 | **§2 Correlation-aware / weighted consensus** (de-herd cross-aligner agreement) | NOVEL | S–M | High |
| 4 | **§4 Calibrated per-read confidence + abstain** | NOVEL | M | High (downstream APA) |
| 5 | **§3 Ground-truth-free aligner weighting** (NET-seq / replicate concordance) | NOVEL | M | Medium–High |
| 6 | **§6 Learned selector** (gradient-boosted ranker over read features) | NOVEL | L | Medium (ceiling probe) |

Rationale: §5 is the cheapest and gates the credibility of everything else (you cannot
claim §1–§6 "helped" without a regenerable number). §1 fixes the two known biases in the
*already-running* HP-ED selector (free N-ops; no explicit 3'-end term). §2/§4 are small
principled changes with outsized correctness leverage. §3 supplies the weights §2/§6
consume. §6 is the research ceiling-probe, last because it needs §5's harness + §3's labels.

---

## §5 — Reproducible win-rate measurement harness + provenance  [ESTABLISHED, Rank 1]

**Mechanism.** A committed script `dev/aligner_investigation/04_discovery/measure_winrates.py`
that: (1) takes a directory of committed per-aligner `corrected_reads.tsv` (+ corrected or
raw BAMs), (2) calls `merge_corrected_tsvs` under both selectors — the production HP-ED
path (BAMs+genome+penalty table) and the legacy fallback sort (no BAMs) — (3) writes
`aligner_summary.tsv` with `winning_aligner` counts, %, n-present denominator, and per-read
tie statistics for each, and (4) emits a `PROVENANCE.json` (input dataset, git SHA, both
v3.3.0 bug fixes confirmed-applied, command line, date). Add a pytest that asserts the
committed summary reproduces to ±0 reads from the committed TSVs.

**Bias it fixes.** The numbers in `CLAUDE.md` exist only as a bare assertion with no
committed artifact; they may predate the `index_col` and `_pt:i:N` fixes. Without
regeneration, no proposed change can be evaluated, and the live HP-ED selector's behaviour
cannot be separated from its calibration weaknesses (§1).

**Feasibility.** Pure orchestration of an existing function. No core changes. Inputs:
the committed `wt_by4742_rep1` validation per-aligner TSVs already in
`rectify/data/validation/`. ~80 lines.

**Expected impact.** Converts every downstream claim from "asserted" to "measured."
Likely outcome: the HP-ED selector and the legacy fallback produce **materially different**
win rates — quantifying how much of the published spread is selector-bound, and giving §1's
N-op-cost / 3'-end-term changes a measurable baseline.

**Validation/ablation.** The harness *is* the validation substrate for §1–§6. Add the
"disable `_n_agree`/`_conf_rank`, random tie-break" ablation (fallback path) and the
denominator-normalization (per-aligner n-present) as flags.

**Risk.** None to production (dev-only). Only risk is that regenerated numbers
contradict `CLAUDE.md` — in which case `CLAUDE.md` must be corrected (a feature).

---

## §1 — Calibrate the live HP-edit-distance selector  [ESTABLISHED, Rank 2]

**Mechanism.** The HP-edit-distance path already runs in production: both `merge_corrected_tsvs`
call sites pass per-aligner BAMs + genome, so `use_hp_ed=True` and reads are ranked on a
genome-anchored quality score over the corrected CIGAR. The work here is not to wire it on
but to fix the two known calibration weaknesses in the score so it ranks for CPA accuracy,
not merely CIGAR fit.

Two changes, both in `_cigar_hp_edit_distance` / the merge key in
`corrected_consensus.py`:

1. **Charge N-ops a calibrated cost.** Today an intron N-op costs **0** (a free pass);
   the only backstop is `_cigar_min_junction_anchor`, a flank-anchor gate that defaults to
   `min_junction_anchor_bp=0` on yeast (off; human uses 10 bp). So an aligner inventing a
   long false intron pays nothing and can *lower* its edit distance to win. **Fix:** add a
   per-N-op opening cost — a flat `+1.0`, or, better, a splice-strength-discounted cost that
   waives the charge for canonical / annotated / cross-read-supported junctions and applies
   it to unsupported ones. This is the single live open lever the build exposes. NOVEL sub-fix.
2. **Add an explicit 3'-end term to the key.** HP-edit-distance scores CIGAR fit, not CPA
   position. **Fix:** add a 3'-end tiebreaker = distance from the read's `corrected_3prime`
   to a **weighted** cross-aligner consensus position (see §2, not raw `_n_agree`), inserted
   *after* `hp_edit_distance`, *before* `_span`.

**Bias it fixes.** The free-N-op bias (a fabricated intron buys a lower edit distance with no
penalty on yeast) and the absence of a CPA-position term in the ranking key. Both are
properties of the *live* selector, not of a dead code path.

**Expected impact.** High but **must be measured (§5), not assumed**: a calibrated N-op cost
re-ranks exactly the contested spliced reads where a spurious intron currently wins, and the
3'-end term aligns the selector with the actual deliverable (CPA). Both changes can move the
win-rate spread.

**Validation/ablation.** §5 harness: ablate {N-cost on/off}, {3'-term on/off}, and sweep the
junction-anchor gate (`min_junction_anchor_bp`) from 0 → human-style 10. Validation set:
`tests/test_consensus_selection.py` (40 tests, real `wt_by4742_rep1`) must still pass; extend
it to assert the N-op cost re-ranks the known free-intron cases on the bundled data.

**Risk.** Medium. The HP path already runs on production BAMs, so the added cost is only the
extra per-N-op accounting, not new I/O. Behavioural risk: a too-aggressive N-op charge could
penalize genuine introns at noisy junctions → tie it to splice strength / the anchor gate, not
a flat penalty, and gate the change behind a flag defaulting **on** only after §5 confirms it
improves ground-truth concordance (§3).

---

## §2 — Correlation-aware / weighted consensus (de-bias cross-aligner agreement)  [NOVEL, Rank 3]

**Mechanism.** Any cross-aligner agreement signal — the §1 3'-end consensus term in the live
HP-ED path, and `_n_agree` in the fallback sort — counts aligners landing on the same
`corrected_3prime` with **equal weight 1.0**. But three of five aligners share *minimap2
lineage* (minimap2 native; gapmm2 wraps minimap2 for the body; uLTRA falls back to minimap2
on >10% out-of-annotation reads). deSALT (RdBG) and mapPacBio (BBMap, independent Java engine)
are the only structurally-independent callers. So a "3-vote majority" can be **one shared
minimap2 decision counted three times.**

Replace the raw count with an **effective vote** that down-weights correlated aligners:

```
effective_agree(pos) = Σ_a  w_a · I[corrected_3prime_a == pos]
```

where `w_a` discounts membership in a correlation cluster. Concretely, assign each
aligner a **lineage weight**: independent callers (deSALT, mapPacBio) get `w=1.0`;
minimap2-family members share a pooled budget (e.g. each gets `1/√k` or a fixed
`{minimap2:0.5, gapmm2:0.3, uLTRA:0.5}`), so the family cannot outvote a single
independent-but-correct aligner. A principled version estimates the correlation matrix
empirically (§3 / §5 harness): pairwise agreement rate `ρ_ab` over a held-out set;
weight = inverse of summed redundancy (a Vandermonde/Hoeffding-style decorrelation, or
simply `w_a ∝ 1 / Σ_b ρ_ab`). This is the standard *correlated-ensemble* correction.

**Bias it fixes.** Herd bias: an aligner whose *output style* makes the minimap2 family
agree wins an agreement tiebreaker independent of truth. deSALT's "78.9%" is partly
cross-read-pooling making others agree, not raw 3' accuracy (the deSALT dossier itself
calls deSALT's *raw* 3' ends "imprecise").

**Feasibility.** Localized to `corrected_consensus.py`. Replace the equal-weight
`groupby().size()` agreement count with a weighted aggregate (map `_aligner → w_a`, sum
weights per `(read_id, corrected_3prime)`). Add a `--aligner-weights` TSV (default = lineage
weights baked in). ~30 lines. Apply identically to the §1 3'-end consensus term and the
fallback `_n_agree` so the live HP-ED path and the fallback sort share one de-herded
consensus definition.

**Expected impact.** High on *correctness of contested reads* (reads where the
independent caller disagrees with the minimap2 herd are exactly the hard cases). Modest
on aggregate win-rate magnitude; large on *which* aligner wins contested reads.

**Validation/ablation.** §5 harness with {equal weights} vs {lineage weights} vs
{empirical-ρ weights}. Ground truth (§3): does de-herding move contested-read CPAs
*toward* NET-seq peaks? Leave-one-family-out test (redteam exp. 4): drop the minimap2
family; if independent-aligner accuracy is unchanged but its win rate jumps under equal
weights and is *stable* under lineage weights, the de-herding worked.

**Risk.** Low-medium. Weights are a prior; a bad weight vector could over-penalize a
family that is genuinely right on some loci. Mitigate: weights only break *ties* /
contribute one term, never gate; calibrate empirically (§3) rather than hand-set.

---

## §3 — Ground-truth-free per-aligner accuracy estimation → weights  [NOVEL, Rank 5]

**Mechanism.** Estimate each aligner's *reliability* without a labeled truth set, using
three orthogonal internal signals, then feed the result as `w_a` into §2 / §6:

1. **NET-seq concordance (3' end).** RECTIFY already refines/validates corrected CPAs
   against bundled NET-seq tables (`rectify/core/netseq/`, `data/*netseq*.tsv.gz`). For
   each aligner, on reads where it is *present*, compute the fraction whose corrected
   3' end lands within ±k bp of a NET-seq-supported position. This is an
   *aligner-independent* 3'-end oracle (an orthogonal assay, not another aligner).
2. **Replicate concordance.** On reads shared across biological replicates
   (`wt_*_rep1/rep2`), measure per-aligner reproducibility of the corrected CPA at a
   locus. A reliable aligner produces the same corrected position across replicates.
3. **Internal consistency / leave-one-out.** For each read, compute each aligner's
   agreement with the **weighted** consensus of the *others* (excluding itself, to avoid
   self-reinforcement). Aligners that systematically agree with the orthogonal majority
   on easy reads earn higher weight (a Dawid–Skene / EM-style annotator-reliability
   estimate, ESTABLISHED in crowd-labeling; NOVEL here as applied to aligners).

Combine into `w_a` = product (or logistic blend) of the three; recompute periodically,
commit the table.

**Bias it fixes.** Equal voting (every per-aligner accuracy claim is currently UNTESTED).
Replaces "all aligners agree-vote equally / self-report confidence" with data-driven
reliability, and supplies the empirical correlation/weight inputs §2 and §6 need.

**Feasibility.** Medium. NET-seq concordance reuses `netseq_refiner.py` machinery.
Dawid–Skene EM is ~50 lines (or `scikit-learn`-free hand roll). Lives in a new
`dev/.../estimate_aligner_reliability.py` + a committed `aligner_weights.tsv`.

**Expected impact.** Medium–high: turns §2's hand-set lineage weights into measured
ones, and gives the first *truth-anchored* per-aligner accuracy numbers RECTIFY has
ever had (directly answers redteam experiment 3). NET-seq is the field-standard 3'
oracle, so this is the most defensible weighting signal.

**Validation/ablation.** Cross-check the three signals against each other (do NET-seq,
replicate, and EM rankings agree?). Hold out one signal, predict it from the other two.
Multi-dataset replication (redteam exp. 7): weights should be stable across DRS samples;
if deSALT's NET-seq concordance ≠ its win rate, the win rate was a metric artifact.

**Risk.** Medium. NET-seq is yeast-WT-specific (`_CANONICAL_HP_PRIOR`, tables are
"R10.4.1-specific"); weights must be regenerated per organism/chemistry — document this
loudly. Replicate concordance needs replicates (not always available). EM can converge
to a degenerate "trust the herd" solution → must exclude self and seed with NET-seq.

---

## §4 — Calibrated per-read confidence + abstain  [NOVEL, Rank 4]

**Mechanism.** Emit a per-read `selection_confidence ∈ [0,1]` column in
`corrected_reads.tsv`, calibrated so it means "P(this corrected CPA is within ±k bp of
truth)." Features per read: winner's hp_edit_distance, margin to runner-up
(Δhp_edit_distance), **weighted** cross-aligner 3'-end agreement (§2), splice/canonical
strength of any N-op, winner's NET-seq support (§3), n-aligners-present, chimera flag.
Fit an isotonic / Platt calibration of a simple score against the NET-seq oracle (§3) on
held-out reads, so the number is *calibrated*, not a raw heuristic rank. Add an
**abstain** tier: when margin is tiny AND no orthogonal support (NET-seq/independent
aligner) confirms the position, set `selection_confidence` low and a `_abstain=1` flag
so downstream APA can filter rather than commit a coin-flip CPA.

**Bias it fixes.** Today the emitted `confidence` is `select.py`'s 3-bucket *aligner-agreement*
flag (`_CONFIDENCE_RANK high/med/low`), self-assessed and uncalibrated, and it carries no
notion of the HP-ED margin the live selector actually used. Downstream APA/DESeq2 treats every
selected CPA as equally trustworthy; contested coin-flip reads contaminate cluster boundaries
and differential 3'-UTR calls.

**Feasibility.** Medium. The features are already computed in the merge (or added by
§1–§3). Calibration is `sklearn.isotonic`/`IsotonicRegression` or hand-rolled. Add the
column in `merge_corrected_tsvs` result; thread a `--min-selection-confidence` filter
into `rectify analyze`.

**Expected impact.** High for the **downstream** science (the actual deliverable is APA
calls): filtering/weighting by calibrated confidence should sharpen cluster boundaries
and reduce false differential-APA hits, even if it does not change win rates at all.
Abstain prevents the ensemble from manufacturing precision it doesn't have.

**Validation/ablation.** Reliability diagram (predicted vs empirical accuracy vs
NET-seq). Show that high-confidence reads have higher NET-seq concordance than
low-confidence. Downstream: does filtering at confidence ≥ τ improve replicate
concordance of APA clusters? Sweep τ.

**Risk.** Low–medium. Mis-calibration on out-of-distribution data (novel organism)
could over/under-filter — gate abstain behind a flag, default conservative. Dropping
abstained reads reduces depth; expose as opt-in.

---

## §3-vs-§6 bridge / §6 — Could a learned selector beat the heuristic sort?  [NOVEL, Rank 6]

**Mechanism.** Frame selection as **learning-to-rank** over the per-read aligner
candidates. Features per (read, aligner): hp_edit_distance, Δ-to-runner-up,
N-op count / splice strength / canonical tier / junction-anchor length, **weighted**
cross-read support (§2), aligner identity (one-hot), aligner reliability weight (§3),
chimera flag, aligned_bases, span, 5'-rescue flag. Label = "is the corrected CPA within ±k bp
of the NET-seq / replicate-consensus truth" (§3). Train a gradient-boosted ranker (LightGBM
`lambdarank`) or a small pairwise logistic model; the live HP-ED selector is the baseline.

**Bias it fixes.** The live HP-ED selector ranks on a single scalar edit distance plus a
fixed lexicographic tail (`_effective_chimera_ok > hp_edit_distance > _span`), with N-ops
free and no learned feature interactions. A learned ranker discovers interactions (e.g.
"trust mapPacBio's CPA only when no splice nearby", "discount an N-op only when its anchor is
short") that the fixed HP-ED key plus the §1 hand-set N-op cost cannot express — it augments
or replaces the HP-ED scalar with a calibrated, multi-feature score.

**Feasibility.** Large — needs §5 (harness), §3 (labels), and a feature-extraction pass.
Keep it **offline/advisory first**: ship as `dev/.../learned_selector.py` that re-ranks
the §5 harness output; only promote into `merge_corrected_tsvs` (behind
`--selector model.txt`) if it beats the live HP-ED selector on held-out NET-seq concordance
by a margin that survives cross-dataset validation.

**Expected impact.** Medium and *uncertain* — this is a ceiling probe. If a learned
ranker barely beats a calibrated HP-ED selector (§1) + de-herded consensus (§2), that is
itself the valuable finding: it means the HP-ED score captures the signal and the win rates
are metric-bound, not model-bound. If it wins big, there is real interaction structure left
on the table.

**Validation/ablation.** Train/test split by **chromosome** (avoid locus leakage) and by
**dataset** (the real test: does it generalize off `wt_by4742_rep1`?). Feature-importance
+ ablation to ensure it is not just relearning "trust the annotation/herd." Compare to
§1+§2 heuristic head-to-head on NET-seq concordance, not on win rate.

**Risk.** High. Overfitting to one yeast DRS sample (redteam confounder 5); learning
the annotation-circularity bias (redteam claim 8) if uLTRA's GTF-snap leaks into labels;
opacity vs the current auditable sort. Mitigate: orthogonal labels (NET-seq, not
annotation), cross-dataset gate, keep heuristic as default fallback.

---

## Cross-cutting notes

- **Sequence-first, never gate.** Every weighting/prior above (lineage weights, canonical
  N-cost, reliability) must enter as a *tiebreaker/soft term*, never a hard filter — this
  is the established Module 2H / `junction_refiner.py` policy (`_CANONICAL_HP_PRIOR=0.5`,
  "annotation/canonical are tie-breakers, never gates"). Re-introducing a candidate gate
  is explicitly forbidden in CLAUDE.md.
- **The deprecated raw-BAM consensus path** (`scoring.py` junction-proximity penalty,
  `select.py` tiebreakers) is **not** on the correct-first path and should be either deleted
  or clearly marked deprecated so future readers stop attributing production behaviour to it.
  (Note this is distinct from the live raw-BAM HP-ED path inside `merge_corrected_tsvs`, which
  *is* production: when only raw BAMs are available, HP-ED is computed lazily in memory.)
- **Order of operations:** §5 → §1 → (§2 + §3 together, since §2's empirical weights come
  from §3) → §4 → §6. §5 must land first or none of the others can be evidenced.
