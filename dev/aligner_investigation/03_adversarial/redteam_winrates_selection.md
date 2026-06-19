# Red-Team: Win Rates vs. Selection Mechanism

**Scope:** Adversarial review of the causal chain linking aligner algorithms to RECTIFY's
reported correct-first win rates (deSALT 78.9%, mapPacBio 18.2%, uLTRA 2%, gapmm2 0.8%,
minimap2 0.1%), reading the dossiers in `01_investigation/`, the synthesis in
`02_synthesis/COMPARISON.md`, and the actual selection code.

**Bottom line:** The win-rate numbers are an artifact of a *single dataset* run through a
*selection metric the dossiers describe incorrectly*. The dossiers (and COMPARISON.md) assert
selection happens "on corrected 3'-end position" and the module docstring says HP-aware edit
distance is the primary criterion — **but both production call sites of `merge_corrected_tsvs`
pass no BAMs, which disables the HP-edit-distance path entirely and runs a legacy 5-level sort
in which the 3'-end position is not even a sort key.** Almost every "this aligner wins because
it is more accurate" claim is therefore an UNTESTED HYPOTHESIS, and several are LIKELY BIASED
toward an aligner's output *style*, not its accuracy.

---

## (a) The actual winner-selection mechanism (from the code)

Production winner selection lives in
`rectify/core/consensus/corrected_consensus.py::merge_corrected_tsvs`. There are **two distinct
sort paths**, switched by `use_hp_ed = bool(per_aligner_corrected_bams)` (line 662):

### Path A — HP-edit-distance sort (the one the docstring + dossiers describe)
Active only when `per_aligner_corrected_bams` is supplied. Sort key (lines 740–747):

```
sort by: read_id, _effective_chimera_ok (asc), hp_edit_distance (asc), _span (desc)
```

- `hp_edit_distance` (`_cigar_hp_edit_distance`, lines 57–132): sum over the corrected CIGAR of
  `X`/`M`-mismatch = 1.0, `D` = penalty-table del_cost (fallback 1.0), `I` = ins_cost (fallback
  1.25), `S` = 1.0/base, **`N` (intron) = 0 (free)**, `=` = 0.
- `_effective_chimera_ok` = `_chimera_ok & (_five_rescued==0)` — overhang/anchor filter
  (`_add_chimera_flag`, lines 331–482), gated by an `OverhangTable`.
- 3'-end position is **not** in this key. (It only appears as a *tiebreaker* inside the
  separate raw-BAM `select.py` path; see below.)

### Path B — Legacy 5-level sort (what actually runs in production)
Active when no BAMs are passed → `use_hp_ed = False`. Sort key (lines 749–752):

```
sort by: read_id, _five_rescued (desc), _chimera_ok (asc), _conf_rank (desc),
         _n_agree (desc), _span (desc), _n_junc (desc)
```

- `_five_rescued` = the `five_prime_rescued` flag (Cat3 rescue).
- `_conf_rank` = `confidence` ∈ {high:2, medium:1, low:0} — a flag set *during correction by the
  aligner's own read*, not a cross-aligner accuracy measure.
- `_n_agree` = count of rows sharing the same `(read_id, corrected_3prime)` — i.e. how many
  aligners landed on the *same corrected 3' position* (lines 693–698). This is the only place
  3'-end agreement enters, and it is a **majority-vote popularity** term, not an accuracy term.
- `_span`, `_n_junc` = wider alignment / more junctions preferred.

### Which path runs in production?

Both call sites pass **no** `per_aligner_corrected_bams`:

- `rectify/core/commands/run/single_sample.py:495` — `merge_corrected_tsvs(per_aligner_tsvs=…,
  output_tsv=…, summary_tsv=…, overhang_table=_overhang_table)`. No BAMs, no genome, no penalty
  table.
- `rectify/core/commands/split_command.py:985` (the chunked SLURM array path — the mandated
  production path per CLAUDE.md) — `merge_corrected_tsvs(valid_tsvs, …, overhang_table=…)`. No
  BAMs. (It *builds* `corrected_bams` at line 970 and uses them to assemble the consensus BAM,
  but never feeds them to the merge.)

**Therefore the reported win rates were almost certainly produced by Path B (legacy sort),
whose primary discriminators are `_five_rescued`, `_chimera_ok`, `_conf_rank`, and 3'-position
*popularity* (`_n_agree`) — not HP-edit-distance and not 3'-end accuracy.** The
`junction-proximity penalty`, `score_alignment`, and `select_best_alignment` in `scoring.py` /
`select.py` belong to the **raw-BAM `consensus.py` path** (the deprecated consensus-first
ordering), and are **not invoked at all** by the correct-first merge. The dossiers repeatedly
cite `scoring.py:745` (junction-proximity penalty) as a cause of the win rates — that penalty
does not run on the correct-first path that produced the numbers.

> Note: even in Path A, `N` ops are free (cost 0). An aligner that invents a long false intron
> pays nothing for the skipped reference and can *lower* its edit distance by replacing
> mismatched bases with a free `N`. The chimera/overhang filter exists precisely to paper over
> this — but in production (Path B) the overhang filter is the *only* defense against false
> introns, and it is a binary flag, not a graded score.

---

## (b) Causal-claim ratings

Rating key: **SUPPORTED** = backed by code structure or a cited measurement on this data;
**UNTESTED-HYPOTHESIS** = plausible, self-labeled `[INFERENCE]`, no ablation/ground-truth;
**LIKELY-BIASED** = the claimed cause is better explained by the selection metric rewarding an
output *style* than by accuracy.

| # | Claim (source) | Rating | Why |
|---|---|---|---|
| 1 | The win-rate numbers themselves (78.9/18.2/2/0.8/0.1) (CLAUDE.md L81; every dossier) | **UNTESTED-HYPOTHESIS (as "biology")** | Numbers exist only in CLAUDE.md L81. No committed script, count table, notebook, or `aligner_summary.tsv` in the repo derives them. `split_command.py:994` *can* print them per-chunk to a log, but no artifact is checked in. COMPARISON.md §6.5 concedes they "trace to one validated lineage… one experiment." |
| 2 | "RECTIFY selects per-read winners on **corrected 3'-end position**" (COMPARISON.md L100/L115/L121; deSALT §3 L144) | **LIKELY-BIASED / FACTUALLY WRONG** | The production sort (Path B) never sorts on 3'-end position. 3'-end enters only as `_n_agree` (popularity vote) and only as a *tiebreaker*. The entire "two properties dominate that metric" narrative (COMPARISON §4) is built on a metric that isn't the one running. |
| 3 | "HP-aware edit distance is the primary winner criterion" (corrected_consensus.py docstring L9–17; ont_drs L121) | **UNTESTED / NOT-IN-PRODUCTION** | True of Path A only. Path A is not wired in either production call site. So the docstring describes a code path the pipeline doesn't take. |
| 4 | deSALT 78.9% because cross-read exon pooling yields homogeneous/consistent junctions "exactly what RECTIFY's correct-first selection scores" (deSALT §3 L142–157; COMPARISON L115) | **LIKELY-BIASED** | This is the clearest style↔metric coupling. In Path B, homogeneous junctions → many aligners land on the same `corrected_3prime` → high `_n_agree`. An aligner whose *style* makes other aligners agree with it wins a *popularity* tiebreaker, independent of whether the shared position is biologically correct. The dossier itself says deSALT's *raw* 3' ends are "imprecise" (§7 L246–251) — so it cannot be winning on raw 3' accuracy. Self-labeled `[INFERENCE]`. |
| 5 | mapPacBio 18.2% because full affine DP gives "column-optimal… best raw 3'-end (CPA) precision," "closer to the true CPA," "hard to beat" (mapPacBio §"Why 2nd-best" L91–101, L231–252; COMPARISON L116) | **UNTESTED-HYPOTHESIS (circular)** | No measurement of CPA accuracy per aligner exists. "Closer to true CPA" is inferred *from* the win rate, then used to *explain* the win rate. Self-labeled `[INFERENCE from … empirical win rate]`. Also: Path B does not score CPA precision at all. |
| 6 | minimap2 0.1% / gapmm2 0.8% because RECTIFY's junction-proximity penalty "docks" them (gapmm2 §5.6 L146, §7 L166c; minimap2 §6 L303–311; COMPARISON L118) | **LIKELY-BIASED (wrong mechanism) / NOT-IN-PRODUCTION** | The junction-proximity penalty lives in `scoring.py` (raw-BAM consensus path) and is **not called** by the correct-first merge. So whatever suppresses gapmm2/minimap2 in production, it is not this penalty. The dossiers attribute the loss to a metric that doesn't run. |
| 7 | minimap2 loses because `--end-bonus 0` → soft-clips poly-A, drifting `reference_end`; HP indels distribute arbitrarily (minimap2 §6 L296–306) | **UNTESTED-HYPOTHESIS (partly moot)** | Plausible as raw-alignment behavior, but RECTIFY's `find_polya_boundary` reassigns the CPA regardless of where minimap2 stopped (scoring.py L724–732 explicitly says raw 3' endpoint is NOT scored). So raw soft-clipping cannot directly explain a *corrected*-output loss. Mechanism likely misattributed. |
| 8 | uLTRA 2% because GTF-coordinate snapping gives "zero wobble" that "beats the noisier de-novo aligners" (uLTRA §"Strengths" L116–117; §"Novel vs Annotated" L92–94) | **LIKELY-BIASED (annotation-circularity)** | If the metric rewards canonical/annotated junctions (Path A tiebreakers do; Path B's `_n_agree`/`_conf_rank` can correlate), an aligner that *hard-snaps to the annotation by construction* will win on annotated reads *whether or not the read's true junction matches the GTF*. "Snapped-to-GTF = correct" is assumed, not tested. uLTRA admits **no** 3'-end precision (§"Weaknesses" L125), yet supposedly wins on a "3'-end metric" — internal contradiction. |
| 9 | "Win rates reflect biology" (CLAUDE.md L81) | **LIKELY-BIASED** | They reflect a selection-metric output on one dataset. The same reads run through Path A vs Path B would likely give *different* win rates (the docstring vs production divergence proves the metric is not unique). |
| 10 | deSALT's `aln` runs the `null` (~13% error) model, not an ONT preset (COMPARISON §6.1) | **SUPPORTED (as a fact) / UNTESTED (as "fine")** | Source-verified that no `-x` is passed. But "still wins 78.9% so it's fine" is post-hoc; no A/B vs `-x ont2d`. |
| 11 | mapPacBio's `=`/`X` CIGARs + clean `N` reclassification "give RECTIFY clean ops to operate on" (mapPacBio L108–113, L244–246) | **SUPPORTED (fact) → LIKELY-BIASED (as win-rate cause)** | The CIGAR style is real. But Path A's `_cigar_hp_edit_distance` treats `=` as 0 and `M`-with-genome-match as 0 — so an aligner that emits `=`/`X` gets the *same* score as one emitting `M` for identical bases. The style helps downstream *surgery*, not the *score*. Citing it as a win-rate cause conflates consumability with winning. |

---

## (c) Confounders that could explain the win rates other than raw accuracy

1. **`_n_agree` is a popularity vote, not accuracy.** In Path B the dominant non-flag
   discriminator is how many aligners share a read's `corrected_3prime`. Whichever aligner sits
   in the *majority cluster* wins ties — so a *correlated-error* family (e.g. deSALT + gapmm2 +
   uLTRA all built on minimap2-like chaining, or all snapping to the same annotation) inflates
   the winner that matches the herd, regardless of truth. This is herd bias, not accuracy.

2. **Per-aligner read-count denominator / coverage gaps.** Win rate = wins / reads-where-this-
   aligner-was-present-and-others-too. uLTRA produced empty SAMs on chunks 0–6 of
   `wt_by4742_rep1_chunked_20260412` (BUGS_TO_FIX.md L179) and CLAUDE.md notes uLTRA/deSALT
   "stay at 34/36" on validation. deSALT also failed on chunks 2/11/15 (penalty-table notes).
   If aligners are differentially *present*, the win-rate denominators differ and the
   percentages are not a like-for-like comparison.

3. **`_pt:i:N` read-id normalization (CLAUDE.md v3.3.0).** mapPacBio reads are `UUID_pt:i:N`.
   `_normalize_read_id` (lines 507–528) and `_bare_uuid` strip both space and underscore forms.
   If *any* historical merge ran before this fix, mapPacBio rows didn't join to other aligners
   → mapPacBio won 100% of its *own isolated singleton groups* with no competition (CLAUDE.md
   says exactly this). A win rate measured pre-fix is inflated for mapPacBio; post-fix it should
   *drop*. **Whether the cited 18.2% is pre- or post-fix is unconfirmed.**

4. **`index_col` bug (CLAUDE.md v3.3.0).** Before the fix, minimap2/gapmm2/uLTRA TSVs had every
   column shifted one right (read_id → chromosome), so their `read_id` was a chromosome name
   (~17 unique values) and they **never properly competed**. The 0.1%/0.8%/2% figures are
   *suspiciously consistent* with "these three barely competed." If the cited rates predate the
   fix, the low-three numbers are a parsing artifact, not biology. **Provenance (pre/post-fix)
   is unverified** — COMPARISON.md §6.5 raises exactly this and marks it "⚠ Confirm."

5. **Single dataset, single organism, single chemistry.** All numbers derive from one
   `wt_by4742_rep1` DRS run (R10.4.1, S. cerevisiae, compact well-annotated genome, short
   introns). The penalty tables and `_CANONICAL_HP_PRIOR=0.5` are "R10.4.1-specific" (pacbio
   §6.3). uLTRA's niche ("complex multi-isoform genomes," uLTRA §"Weaknesses") is *suppressed*
   in yeast — so its 2% is a yeast-specific floor, not a property of uLTRA. Generalizing the
   spread to "DRS data" or "biology" is over-generalization.

6. **Chimera/overhang filter as a deSALT/gapmm2 *gate*.** `_add_chimera_flag` exists specifically
   to stop "chimeric deSALT/gapmm2 long-range false junctions from winning" (docstring L584).
   This is a *correction* to a known deSALT/gapmm2 failure mode — i.e., deSALT's win rate is
   already the product of a filter tuned against deSALT's own artifacts. The 78.9% is
   filter-conditioned, not raw.

7. **Confidence flag is self-assessed.** `_conf_rank` derives from each aligner's *own*
   `confidence` column, set during its *own* correction (bam_processor.py). An aligner that
   self-reports "high" more readily wins the `_conf_rank` tiebreaker without being more correct.

---

## (d) Experiments to disentangle accuracy from selection bias

1. **Provenance audit (do first, cheap).** Confirm whether the cited 78.9/18.2/2/0.8/0.1 were
   measured *after* both v3.3.0 fixes (`index_col=False` and `_pt:i:N` normalization). Re-run
   `merge_corrected_tsvs` on the committed `wt_by4742_rep1` per-aligner TSVs and diff the
   `winning_aligner` counts against CLAUDE.md L81. If they don't match, L81 is stale and must be
   corrected. **Commit the resulting `aligner_summary.tsv` so the numbers stop being a bare
   assertion.**

2. **Path A vs Path B ablation.** Run the *same* TSVs through `merge_corrected_tsvs` with and
   without `per_aligner_corrected_bams` (HP-edit-distance vs legacy sort). If the win-rate spread
   changes materially, the numbers are a *metric artifact*, definitively. This also forces a
   decision: the docstring claims HP-edit-distance is primary but production disables it — one of
   the two must be wrong. **This is the single most important experiment.**

3. **Ground-truth concordance, per aligner.** The win rate measures "won the selection," not
   "was correct." Build a held-out truth set (e.g. Quant-seq / 3'-seq CPA calls, or annotated
   intron coordinates for spliced reads) and measure each aligner's *corrected* 3'-end / junction
   accuracy **independent of selection**. Only then can "deSALT is more accurate" be separated
   from "deSALT wins the metric." Without this, every accuracy claim stays UNTESTED.

4. **Leave-one-aligner-out / herd-bias test.** Remove deSALT from the ensemble and re-measure. If
   mapPacBio's accuracy-vs-truth (experiment 3) is unchanged but its *win rate* jumps, the win
   rate was a function of the competitor field (popularity), not accuracy. Likewise drop the
   minimap2-family (minimap2+gapmm2+uLTRA share chaining) to test `_n_agree` herd inflation.

5. **Disable `_n_agree`/`_conf_rank` tiebreakers.** Re-rank using only quality terms (HP-edit-
   distance, span) with random tie-breaking. Quantify how many wins were decided *purely* by the
   popularity/self-confidence tiebreakers. A large fraction would show the spread is a tiebreaker
   artifact.

6. **Annotation-circularity test for uLTRA.** Restrict to reads whose true junction is *novel*
   (not in the GTF) using an orthogonal junction truth. uLTRA's GTF-snapping advantage should
   vanish; if uLTRA still "wins" novel reads, the metric is rewarding canonicality over truth.

7. **Multi-dataset / multi-chemistry replication.** Re-measure on a second DRS sample, a
   different organism, and (if possible) a HiFi sample. Per pacbio §4, HiFi has no HP-deletion
   bias, so the spread *should* collapse if it reflects ONT-specific correction; if deSALT still
   wins 79% on HiFi, the metric is selecting on something other than the artifact it claims to
   correct.

8. **False-intron stress test (N-op free-pass).** Inject reads with a known-false long intron and
   verify whether the overhang filter (Path B) or hp_edit_distance (Path A, where `N`=0) lets
   deSALT/gapmm2 win. Quantifies how much of deSALT's lead is filter-conditioned vs. genuine.

---

## Summary of what is and isn't established

- **Established (SUPPORTED):** the *code* mechanism of both sort paths; that the junction-
  proximity penalty / `select.py` belong to the deprecated raw-BAM path and are not on the
  correct-first path; that production runs Path B; that v3.3.0 fixed two bugs that previously
  prevented minimap2/gapmm2/uLTRA and mapPacBio from competing fairly.
- **Not established (UNTESTED / LIKELY-BIASED):** that any specific aligner is *more accurate*;
  that the win rates reflect "biology"; that selection is on "corrected 3'-end position"; that
  the junction-proximity penalty drives gapmm2/minimap2 losses; that the rates generalize beyond
  one yeast DRS sample; and that the cited numbers are even post-v3.3.0-fix.

The dossier authors deserve credit: each self-labels its causal chain `[INFERENCE]` and
COMPARISON.md §6.5 already flags the single-dataset and `index_col` concerns. The unflagged
problems are (i) the *factual* mismatch between the described metric (3'-end / HP-edit-distance)
and the *production* metric (legacy 5-level sort), and (ii) the unexamined style↔metric coupling
in the deSALT (`_n_agree` herd) and uLTRA (annotation-circularity) stories.
