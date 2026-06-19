# Red-Team: Win Rates vs. Selection Mechanism

> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** This file was written against the
> `master`-derived working tree. Its headline finding ("HP-edit-distance never runs in
> production; production runs the legacy popularity vote / Path B") is **retracted on the real
> build**. See `../CORRECTIONS_vs_DRS_BUILD.md` (sections A, D) for the build evidence. Inline
> callouts below mark every retracted / superseded claim.

**Scope:** Adversarial review of the causal chain linking aligner algorithms to RECTIFY's
reported correct-first win rates (deSALT 78.9%, mapPacBio 18.2%, uLTRA 2%, gapmm2 0.8%,
minimap2 0.1%), reading the dossiers in `01_investigation/`, the synthesis in
`02_synthesis/COMPARISON.md`, and the actual selection code.

**Bottom line:** The win-rate numbers are from a *single dataset* and an *un-committed
artifact*, so "this aligner wins because it is more accurate" remains an UNTESTED HYPOTHESIS
and several claims are LIKELY BIASED toward an aligner's output *style*, not its accuracy —
**that critique stands.** What does NOT stand (on the build) is the original mechanism story.

> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** The original bottom line asserted that
> "both production call sites of `merge_corrected_tsvs` pass no BAMs, which disables the
> HP-edit-distance path entirely and runs a legacy 5-level sort." **This is now FALSE on the
> build.** Both production call sites pass `per_aligner_raw_bams` + `genome`, and the gate is
> `use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)`
> (`corrected_consensus.py:1262`); a lazy raw-BAM HP-ED path
> (`corrected_consensus.py:1264-1364`) computes HP-edit-distance in memory without
> materializing corrected BAMs. Wiring: `single_sample.py:238-250` and
> `split_command.py:1085-1094`. **So Path A (HP-edit-distance) runs in production on the build.**
> Path B is the *fallback only* (no BAMs/genome supplied, or staging fails → `ed_rows` empty →
> `use_hp_ed=False`, `:1377`). The "single-dataset, un-committed artifact" caution survives;
> the "produced by Path B / legacy sort" lineage argument does not (see §(a)).

---

## (a) The actual winner-selection mechanism (from the code)

Production winner selection lives in
`rectify/core/consensus/corrected_consensus.py::merge_corrected_tsvs`. There are **two distinct
sort paths**, switched by `use_hp_ed = bool(per_aligner_corrected_bams)` (line 662 *on master*).

> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** On the build the gate is widened to
> `use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)`
> (`corrected_consensus.py:1262`). Passing raw BAMs alone is now enough to activate Path A via a
> lazy raw-BAM HP path (`:1264-1364`). The two-path structure below is still accurate; **which
> path is the default is the part that changed** — see "Which path runs in production?".

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

### Path B — Legacy 5-level sort (the FALLBACK; runs only when no BAMs/genome are supplied or staging fails)
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

> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** The original conclusion below
> ("production runs Path B; the win rates were produced by the legacy sort") is **RETRACTED**.
> On the build, **Path A (HP-edit-distance) runs in production.** Both production call sites pass
> per-aligner **raw** BAMs + genome:
> - **Single-sample:** `rectify/core/commands/run/single_sample.py:238-250` —
>   `merge_corrected_tsvs(…, per_aligner_raw_bams=_staged_s1 if not per_aligner_corrected_bams
>   else None, genome=_merge_genome, …)`.
> - **Chunked/split (mandated production path):** generated script body
>   `rectify/core/commands/split_command.py:1085-1094` —
>   `per_aligner_corrected_bams=… if corrected_bams else None,
>   per_aligner_raw_bams=_staged if (raw_bams and not corrected_bams) else None, genome=genome`.
>
> With raw BAMs + genome supplied, `use_hp_ed = bool(per_aligner_corrected_bams or
> per_aligner_raw_bams)` is **True** (`corrected_consensus.py:1262`); the lazy raw-BAM path
> (`:1264-1364`) computes HP-ED in memory. So the reported win rates, **on a fresh build run**,
> are produced by Path A (HP-edit-distance + overhang + anchor gate), **not** Path B. Path B is
> the documented FALLBACK only (no BAMs/genome supplied, e.g. a manual `merge_corrected_tsvs`
> call, or BAM staging fails → `ed_rows` empty → `use_hp_ed=False`, `:1377`).
>
> The text that follows describes the *master* behavior and is retained for the historical
> record; treat it as a description of the **fallback path**, not production.

*(Historical, master-tree description — superseded above.)* Both call sites *on master* pass
**no** `per_aligner_corrected_bams`:

- `rectify/core/commands/run/single_sample.py:495` — `merge_corrected_tsvs(per_aligner_tsvs=…,
  output_tsv=…, summary_tsv=…, overhang_table=_overhang_table)`. No BAMs, no genome, no penalty
  table.
- `rectify/core/commands/split_command.py:985` (the chunked SLURM array path — the mandated
  production path per CLAUDE.md) — `merge_corrected_tsvs(valid_tsvs, …, overhang_table=…)`. No
  BAMs. (It *builds* `corrected_bams` at line 970 and uses them to assemble the consensus BAM,
  but never feeds them to the merge.)

On master, therefore, the win rates were produced by Path B (legacy sort), whose primary
discriminators are `_five_rescued`, `_chimera_ok`, `_conf_rank`, and 3'-position *popularity*
(`_n_agree`). The `junction-proximity penalty`, `score_alignment`, and `select_best_alignment`
in `scoring.py` / `select.py` belong to the **raw-BAM `consensus.py` path** (the deprecated
consensus-first ordering), and are **not invoked at all** by the correct-first merge — this part
still holds on the build. The dossiers repeatedly cite `scoring.py:745` (junction-proximity
penalty) as a cause of the win rates — that *specific* penalty does not run on the correct-first
path on either tree.

> Note: even in Path A, `N` ops are free (cost 0; `corrected_consensus.py:142-143`). An aligner
> that invents a long false intron pays nothing for the skipped reference and can *lower* its
> edit distance by replacing mismatched bases with a free `N`.
>
> **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** On the build the free-`N` exploit is
> backstopped by a NEW **junction-anchor gate** beyond the binary overhang flag:
> `_cigar_min_junction_anchor` (`:166-220`) measures the perfect-match anchor on both flanks of
> every N-op, and `_add_chimera_flag` (`:868-882`) flags an aligner row chimeric when
> `anchor < min_junction_anchor_bp` (unless every junction is cross-read-supported). **The gate
> is organism-gated** (`rectify/data/__init__.py:1053-1062`; `_MIN_JUNCTION_ANCHOR=0` at
> `corrected_consensus.py:650`): human default **10 bp**, **yeast 0 (gate OFF)**. So "the
> overhang filter is the *only* defense" holds **for yeast only**; on human the anchor gate adds
> a graded structural defense. And because Path A runs in production, hp_edit_distance (with the
> free-`N` caveat) *is* a live selection term, not a dead one.

---

## (b) Causal-claim ratings

Rating key: **SUPPORTED** = backed by code structure or a cited measurement on this data;
**UNTESTED-HYPOTHESIS** = plausible, self-labeled `[INFERENCE]`, no ablation/ground-truth;
**LIKELY-BIASED** = the claimed cause is better explained by the selection metric rewarding an
output *style* than by accuracy.

| # | Claim (source) | Rating | Why |
|---|---|---|---|
| 1 | The win-rate numbers themselves (78.9/18.2/2/0.8/0.1) (CLAUDE.md L81; every dossier) | **UNTESTED-HYPOTHESIS (as "biology")** | Numbers exist only in CLAUDE.md L81. No committed script, count table, notebook, or `aligner_summary.tsv` in the repo derives them. `split_command.py:994` *can* print them per-chunk to a log, but no artifact is checked in. COMPARISON.md §6.5 concedes they "trace to one validated lineage… one experiment." |
| 2 | "RECTIFY selects per-read winners on **corrected 3'-end position**" (COMPARISON.md L100/L115/L121; deSALT §3 L144) | **SUPERSEDED ON BUILD — HP-ED now runs** ~~LIKELY-BIASED / FACTUALLY WRONG~~ | *(master)* The legacy fallback sort (Path B) never sorts on 3'-end position; 3'-end enters only as `_n_agree` popularity. **⚠️ BUILD: on `drs-validation-rebuild` Path A (HP-edit-distance) runs in production (`corrected_consensus.py:1262`; `single_sample.py:238-250`; `split_command.py:1085-1094`), so the production sort key is `(read_id, _effective_chimera_ok, hp_edit_distance, _span)`.** It still is not literally "corrected 3'-end position" — the criterion is HP-edit-distance of the corrected CIGAR, not the 3'-end coordinate — so the dossier phrasing remains imprecise, but the "FACTUALLY WRONG / metric doesn't run" basis is retracted. |
| 3 | "HP-aware edit distance is the primary winner criterion" (corrected_consensus.py docstring L9–17; ont_drs L121) | **SUPERSEDED ON BUILD — HP-ED now runs** ~~UNTESTED / NOT-IN-PRODUCTION~~ | *(master)* True of Path A only, which master did not wire. **⚠️ BUILD: Path A IS wired in both production call sites via the lazy raw-BAM path (`corrected_consensus.py:1264-1364`); the docstring now describes the code path the build pipeline actually takes.** Whether HP-ED selects the *biologically* correct aligner remains UNTESTED (no ground-truth concordance committed). |
| 4 | deSALT 78.9% because cross-read exon pooling yields homogeneous/consistent junctions "exactly what RECTIFY's correct-first selection scores" (deSALT §3 L142–157; COMPARISON L115) | **LIKELY-BIASED** | This is the clearest style↔metric coupling. In Path B, homogeneous junctions → many aligners land on the same `corrected_3prime` → high `_n_agree`. An aligner whose *style* makes other aligners agree with it wins a *popularity* tiebreaker, independent of whether the shared position is biologically correct. The dossier itself says deSALT's *raw* 3' ends are "imprecise" (§7 L246–251) — so it cannot be winning on raw 3' accuracy. Self-labeled `[INFERENCE]`. |
| 5 | mapPacBio 18.2% because full affine DP gives "column-optimal… best raw 3'-end (CPA) precision," "closer to the true CPA," "hard to beat" (mapPacBio §"Why 2nd-best" L91–101, L231–252; COMPARISON L116) | **UNTESTED-HYPOTHESIS (circular)** | No measurement of CPA accuracy per aligner exists. "Closer to true CPA" is inferred *from* the win rate, then used to *explain* the win rate. Self-labeled `[INFERENCE from … empirical win rate]`. Also: Path B does not score CPA precision at all. |
| 6 | minimap2 0.1% / gapmm2 0.8% because RECTIFY's junction-proximity penalty "docks" them (gapmm2 §5.6 L146, §7 L166c; minimap2 §6 L303–311; COMPARISON L118) | **PARTIALLY HOLDS — wrong mechanism named** | The specific `scoring.py:745` junction-proximity penalty (raw-BAM `consensus.py` path) is **still not called** by the correct-first merge — that attribution is wrong on both trees. **⚠️ BUILD: but the correct-first merge DOES now run HP-edit-distance + the overhang filter + (on human) the anchor gate, so a real production mechanism that can suppress gapmm2/minimap2 exists — it just isn't the junction-proximity penalty the dossiers cited.** |
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
   without per-aligner BAMs (HP-edit-distance vs legacy sort). If the win-rate spread changes
   materially, the numbers are sensitive to the selection metric — quantifying *how much* of the
   spread is metric-dependent. **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** on the build
   this is an **ablation between two live paths** (Path A is the production default; Path B is the
   fallback), **not** "wire a dead path on vs. the running one." The original framing — "the
   docstring claims HP-ED but production disables it, one must be wrong" — is moot on the build:
   both are real, Path A is the default, and the docstring matches production. Still a high-value
   experiment, but no longer "the single most important" tie-breaker between docstring and code.

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
  correct-first path; that v3.3.0 fixed two bugs that previously prevented
  minimap2/gapmm2/uLTRA and mapPacBio from competing fairly.
  - **⚠️ BUILD CORRECTION (vs `drs-validation-rebuild`):** "production runs Path B" is
    **RETRACTED**. On the build, production runs **Path A (HP-edit-distance)** — both call sites
    pass raw BAMs + genome (`corrected_consensus.py:1262`; `single_sample.py:238-250`;
    `split_command.py:1085-1094`). Path B is the documented fallback only.
- **Not established (UNTESTED / LIKELY-BIASED):** that any specific aligner is *more accurate*;
  that the win rates reflect "biology"; that the rates generalize beyond one yeast DRS sample;
  and that the cited numbers are even post-v3.3.0-fix. **These cautions are unaffected by the
  build correction** — they concern the *biology*/provenance, not which sort path runs.

The dossier authors deserve credit: each self-labels its causal chain `[INFERENCE]` and
COMPARISON.md §6.5 already flags the single-dataset and `index_col` concerns. The unflagged
problems are (i) ~~the *factual* mismatch between the described metric (3'-end / HP-edit-distance)
and the *production* metric (legacy 5-level sort)~~ — **⚠️ RETRACTED on the build: HP-edit-distance
(Path A) IS the production metric on `drs-validation-rebuild`, so there is no metric mismatch;
the residual imprecision is only that the dossiers say "3'-end position" when the key is
HP-edit-distance of the corrected CIGAR** — and (ii) the unexamined style↔metric coupling
in the deSALT (`_n_agree` herd, which applies only on the Path B *fallback*) and uLTRA
(annotation-circularity) stories. The style↔metric coupling concern (ii) is reframed but not
eliminated: under Path A, an output style that lowers HP-edit-distance (clean `=`/`X`, free `N`)
can still couple to wins independent of biological accuracy — see claims 10/11 and §(c).
