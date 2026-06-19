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

## §1 — Make the HP-aware / corrected-3'-end path actually run  [ESTABLISHED, Rank 2]

**Mechanism.** Collect the already-written per-aligner corrected BAMs and pass them
(plus genome dict + penalty table) into both `merge_corrected_tsvs` call sites, flipping
`use_hp_ed=True`.

Exact changes:
1. **`run/stages.py:_run_correction_per_aligner`** — additionally return a
   `Dict[aligner -> corrected_bam_path]` (the `{stem}.rectified_corrected_3end.bam`
   already produced at line ~194's `corrected_bam_path`). Glob it back per aligner
   exactly as the TSV is globbed at line 327.
2. **`run/single_sample.py:495`** — pass `per_aligner_corrected_bams=<that dict>`,
   `genome=<loaded genome dict>`, `penalty_table=<loaded HpPenaltyTable or None>` to
   `merge_corrected_tsvs`. The genome dict is already loaded upstream for correction;
   thread `--junction-penalty-table` through (it already exists as a CLI flag).
3. **`split_command.py:985`** — the generated script already builds `corrected_bams`
   (line 970). Add `per_aligner_corrected_bams=corrected_bams`, load the genome
   pickle/`--Scer` data, and the penalty table (`OVERHANG_TABLE_PATH` sibling). One
   added kwarg + two loads in the heredoc.

Decision forced by the redteam: the docstring (lines 9–17) says HP-edit-distance is
primary, production disables it. Either wire it (this item) **or** delete the docstring
claim. Wiring is correct — Path B's `_n_agree` is a popularity proxy, not accuracy.

**Bias it fixes.** Replaces majority-popularity (`_n_agree`, herd-amplified) and
self-reported `_conf_rank` with a genome-anchored quality score on the **corrected**
CIGAR. Removes the documented-vs-actual metric mismatch (redteam claims 2–3).

**Expected impact.** High but **must be measured (§5), not assumed**: Path A could
*change* the win-rate spread substantially. Two known Path-A weaknesses to fix in the
same PR (otherwise §1 trades one bias for another):
- **N-ops are free (cost 0)** (lines 124–125): an aligner inventing a long false intron
  pays nothing and can *lower* its edit distance. Currently the only guard is the binary
  chimera/overhang flag. **Fix:** add a small per-N-op opening cost (e.g. `+1.0`, or a
  splice-strength-discounted cost) so an N must "pay its way" with canonical signal /
  cross-read support. NOVEL sub-fix.
- **3'-end position is still absent from the Path-A key.** HP-edit-distance scores CIGAR
  fit, not CPA accuracy. **Fix:** add a 3'-end tiebreaker term = distance from the read's
  `corrected_3prime` to a **weighted** cross-aligner consensus position (see §2, not raw
  `_n_agree`), inserted *after* hp_edit_distance, *before* `_span`.

**Validation/ablation.** §5 harness: Path A vs Path B diff; then ablate
{N-cost on/off}, {3'-term on/off}. Validation set: `tests/test_consensus_selection.py`
(40 tests, real `wt_by4742_rep1`) must still pass; extend it to assert Path A runs (BAMs
present) on the bundled data.

**Risk.** Medium. Path A on 100M-read production means computing HP-edit-distance over
every aligner's full corrected BAM — added I/O + CPU at the merge. Mitigate: the corrected
BAMs are already on `$SCRATCH` at that point; compute per-chunk (split path) so it
parallelizes. Behavioural risk: win rates shift → downstream NET-seq/APA numbers shift;
gate behind a flag (`--hp-edit-selection`) defaulting **on** only after §5 confirms it
improves ground-truth concordance (§3).

---

## §2 — Correlation-aware / weighted consensus (de-bias `_n_agree`)  [NOVEL, Rank 3]

**Mechanism.** `_n_agree` counts aligners landing on the same `corrected_3prime` with
**equal weight 1.0**. But minimap2, gapmm2, and mapPacBio... — actually three of five
share *minimap2 lineage* (minimap2 native; gapmm2 wraps minimap2 for the body; uLTRA
falls back to minimap2 on >10% out-of-annotation reads). deSALT (RdBG) and mapPacBio
(BBMap, independent Java engine) are the only structurally-independent callers. So a
"3-vote majority" can be **one shared minimap2 decision counted three times.**

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

**Bias it fixes.** Herd bias (redteam confounder 1, claim 4): an aligner whose *output
style* makes the minimap2 family agree wins a popularity tiebreaker independent of truth.
deSALT's "78.9%" is partly cross-read-pooling making others agree, not raw 3' accuracy
(the deSALT dossier itself calls deSALT's *raw* 3' ends "imprecise").

**Feasibility.** Localized to `corrected_consensus.py`. Replace the `groupby().size()`
at lines 693–698 with a weighted aggregate (map `_aligner → w_a`, sum weights per
`(read_id, corrected_3prime)`). Add a `--aligner-weights` TSV (default = lineage
weights baked in). ~30 lines. Apply identically as the §1 3'-end tiebreaker term so
Path A and Path B share one de-herded consensus definition.

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

**Bias it fixes.** Equal voting (redteam: every accuracy claim is UNTESTED). Replaces
"all aligners vote equally / self-report confidence" with data-driven reliability,
and supplies the empirical correlation/weight inputs §2 and §6 need.

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

**Bias it fixes.** Today `confidence` is `select.py`'s 3-bucket *aligner-agreement* flag
(`_CONFIDENCE_RANK high/med/low`), self-assessed and uncalibrated (redteam confounder 7).
Downstream APA/DESeq2 treats every selected CPA as equally trustworthy; contested
coin-flip reads contaminate cluster boundaries and differential 3'-UTR calls.

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
N-op count / splice strength / canonical tier, **weighted** cross-read support (§2),
aligner identity (one-hot), aligner reliability weight (§3), chimera flag, aligned_bases,
span, 5'-rescue flag. Label = "is the corrected CPA within ±k bp of the NET-seq /
replicate-consensus truth" (§3). Train a gradient-boosted ranker (LightGBM `lambdarank`)
or a small pairwise logistic model; the heuristic sort is the baseline.

**Bias it fixes.** The hand-tuned lexicographic sort imposes a fixed feature *priority*
(`_five_rescued > _chimera_ok > … `) and hard tier boundaries that the redteam shows are
fragile (e.g. the `int(score)` floor bug history in `junction_refiner.py`). A learned
ranker discovers feature *interactions* (e.g. "trust mapPacBio's CPA only when no splice
nearby") that a fixed sort cannot express.

**Feasibility.** Large — needs §5 (harness), §3 (labels), and a feature-extraction pass.
Keep it **offline/advisory first**: ship as `dev/.../learned_selector.py` that re-ranks
the §5 harness output; only promote into `merge_corrected_tsvs` (behind
`--selector model.txt`) if it beats the heuristic on held-out NET-seq concordance by a
margin that survives cross-dataset validation.

**Expected impact.** Medium and *uncertain* — this is a ceiling probe. If a learned
ranker barely beats a well-wired Path A (§1) + de-herded consensus (§2), that is itself
the valuable finding: it means the heuristic captures the signal and the win rates are
metric-bound, not model-bound. If it wins big, there is real interaction structure left
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
- **The deprecated raw-BAM path** (`scoring.py` junction-proximity penalty, `select.py`
  tiebreakers) is **not** on the correct-first path and should be either deleted or
  clearly marked deprecated so future readers stop attributing production behaviour to it
  (redteam claim 6 shows the dossiers already made this error).
- **Order of operations:** §5 → §1 → (§2 + §3 together, since §2's empirical weights come
  from §3) → §4 → §6. §5 must land first or none of the others can be evidenced.
