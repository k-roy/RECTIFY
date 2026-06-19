# Improvement Roadmap — Long-Read Aligner / 3′-End Program for RECTIFY

**What this is.** The forward-looking deliverable: the six discovery dossiers in `04_discovery/`,
synthesized **through the lens of `04_discovery/REDTEAM_proposals.md`** (which fact-checked every code
claim against RECTIFY source, found the redundancies, and ordered the work by dependency). The six
files proposed ~30 sub-ideas; after de-duplication they reduce to **~9 distinct work items**. This
roadmap lists each retained item once, with its red-team verdict, then a dependency-ordered phased plan.

**Tagging.** **FACT** = code/source-verified. **INFERENCE** = reasoned. **HYPOTHESIS** = plausible,
untested on committed ground truth. Effort: **XS / S / M / L**. Verdict (red-team):
**KEEP / KEEP-WITH-CAVEATS / DESCOPE / REJECT**.

---

## ⚠️ The measurement problem comes first (Phase 0) — and be honest about it

Every discovery file says "build the NET-seq concordance harness first." That is **necessary but not
sufficient**, and the red-team (`REDTEAM_proposals.md §3.3, §4`) exposes why:

- **The cheapest real unblocker needs NO oracle at all.** The win-rate **provenance audit + Path A vs
  Path B ablation** is pure code: re-run `merge_corrected_tsvs` on the committed per-aligner TSVs with
  and without per-aligner BAMs, confirm whether the cited 79/18/2/0.8/0.1 are post-v3.3.0-fix, and
  commit an `aligner_summary.tsv` + `PROVENANCE.json`. This **settles whether the win-rate spread is a
  metric artifact with zero accuracy assumptions** (FACT — both call sites verified to omit BAMs;
  the corrected BAMs Path A needs *already exist on disk*, `stages.py:194`, "simply never collected").
  **Do this first. It is days of work and unblocks every other claim.**

- **NET-seq is circular and is NOT a sufficient oracle.** NET-seq is used as *both* the tuning signal
  (for every proposed model/weight/prior) *and* the validation oracle, *and* it is already partially
  baked into the corrected positions (`netseq_refiner.py` refines against it). Holding out `_wt` vs
  `_pan` does **not** break this — they are not independent experiments (same assay biology, likely
  shared peaks). **Any proposal that trains on NET-seq and reports NET-seq concordance as its win is
  unfalsifiable** until a genuinely orthogonal truth set exists.

- **The real gating cost is a truly independent truth set.** Quant-seq / 3′-seq CPA calls (or, at
  minimum, a held-out NET-seq sample from a *different run*) plus annotated-intron junction truth for
  spliced reads. **None of this is committed in the repo today.** The discovery files treat this as
  "curate the truth set" (low effort); it is in fact the **weeks-long critical path of the whole
  program**. Until it exists, everything below Tier 0 is a hypothesis.

---

## The ~9 retained items (de-duplicated)

Each appears multiple times across the six files; consolidated here once.

| # | Item | Mechanism (1 line) | Failure mode it fixes | Effort | Risk | Verdict | Key caveat |
|---|---|---|---|---|---|---|---|
| **1** | **Win-rate provenance + Path A/B ablation harness** | Re-run the merge both ways on committed TSVs; commit `aligner_summary.tsv` + `PROVENANCE.json` | Win rates are unprovenanced and possibly a metric artifact | XS | Low | **KEEP (do first)** | Needs no oracle; the genuine cheapest unblocker |
| **2** | **Wire Path A** (HP-edit-distance selection) | Collect the per-aligner corrected BAMs that already exist on disk; flip `use_hp_ed=True` | Production sorts on `_n_agree` popularity, not corrected-3′ quality | S | Med | **KEEP-WITH-CAVEATS** | MUST ship with item 3 (N-op cost) or it trades herd bias for the free-intron exploit; flag-gated, not default-on by assumption |
| **3** | **Calibrated N-op cost** | Charge each intron −logP(novel intron here) instead of 0 | `N` ops are free → an aligner invents a false intron to lower its edit distance | S | Low–Med | **KEEP** | Highest correctness-per-effort; must be a *calibrated* cost, not a fixed gate |
| **4** | **Cheap config A/Bs** | deSALT `-x ont2d`; mapPacBio `maxindel`; minimap2 `-k`/`-w` & `--end-bonus` sweep; `pt:i` prior in 2E | Possible free sensitivity on the top aligner; HP 3′ drift | XS–S | Low | **KEEP** | One-line each via `extra_args`; `--end-bonus` is **diagnostic only** (a blunt bonus extends into A-runs — fights item 6) |
| **5** | **Mispriming terminal veto (cDNA)** | Reuse the genomic-downstream-A detector as a bounded terminal-placement penalty under `--dT-primed-cDNA` | Oligo-dT false CPAs on genomic A-runs (QuantSeq) | S | Low | **KEEP** | DRS-exempt by construction; penalty/flag, never a hard gate; self-contained |
| **6** | **ONE unified poly-A-aware terminal CPA solver** | A two-state (genomic→tail) terminal DP OR a window re-scoring pass — the CPA is the body→tail transition | Every aligner soft-clips / slips the 3′ end into a genomic A-run | M | **High** | **KEEP-WITH-CAVEATS** | The only *principled* replacement for the 2B/2E/2G heuristic stack — but it **rewrites the most heavily-patched code path in RECTIFY**; pick ONE design, validate cat1/cat2 shifts exactly; gate on numba availability |
| **7** | **Calibrated junction scorer** | Replace 2H's brittle lexicographic tuple with a log-odds: HP-seq support + yeast MaxEnt/PWM splice strength + cross-read support + small annotation prior | Hand-tuned tier ordering is brittle and not comparable across reads | M | Med | **KEEP-WITH-CAVEATS** | MaxEnt is metazoan-trained → retrain a yeast PWM; only helps *if wired into selection* (depends on item 2) |
| **8** | **Cross-read junction/CPA consensus + de-herded weights** | Pool all 5 aligners' junctions/CPAs into a per-locus consensus; replace raw `_n_agree` with lineage-/reliability-weighted agreement | Per-read idiosyncrasy; `_n_agree` herd bias (3 of 5 aligners share minimap2 lineage) | M | **High** | **KEEP-WITH-CAVEATS** | **Over-homogenization destroys the APA signal that is the product's reason to exist** — snap only within the HP-noise radius; weights must be empirical (no Dawid–Skene EM "trust the herd"), tiebreaker-only |
| **9** | **Calibrated per-read confidence + abstain** | Emit a calibrated `selection_confidence`; flag coin-flip reads so downstream APA can filter | `confidence` is a self-assessed 3-bucket flag; contested CPAs pollute APA/DESeq2 | M | Low–Med | **KEEP-WITH-CAVEATS** | Calibration target is NET-seq → inherits the circularity; `sklearn.isotonic` not importable in this checkout |

### REJECTED (with one-line reasons)

- **KSW2 fork** (in-aligner per-base gap costs) — forks 5 aligners RECTIFY doesn't control; `mappy`
  not even importable; maintenance/reproducibility liability. A windowed RECTIFY-side DP reaches the
  same target.
- **Transformer / seq2seq alignment post-correction** — labeled training set is **36 reads**; replaces
  auditable, validated DP with an opaque black box full of rare CIGAR edge cases.
- **DNA/RNA foundation-model embeddings as a core method** — label-rich, motif-simple problem; GPU cost
  violates the chunked-SLURM design; a GBDT over 15 hand-features wins on accuracy, speed, and
  interpretability. (One speculative chimera-scoring aside, not a pipeline.)
- **Monolithic unified graph index** — a multi-month deSALT-fragility rebuild; the author concedes the
  shippable form is "ensemble glue, not a new index" (= items 5/8). Build the layered approximation only.
- **Standalone syncmer / strobemer aligner** — a from-scratch chainer carries the same correctness tail
  the dossiers catalogue for gapmm2; yeast micro-exon upside is small (real floor ~6 nt, rare in yeast).
- **Remora / POD5 dwell fusion (now)** — requires squiggle retention RECTIFY does not keep + signal↔ref
  anchoring + per-read I/O against the SLURM budget. Defer to a funded R&D track; do the `pt:i` prior
  (item 4) first.

Also deferred-not-rejected: **SpliceAI/Pangolin junction scorer** (NIL for yeast's ~300 canonical
introns; build only the swappable interface, defer the CNN to a metazoan port) and **de-novo micro-exon
recovery** (low yeast value, no validation read exercises it).

---

## Dependency-ordered phased plan

```
PHASE 0 — measurement, no oracle needed (DO FIRST, days)
  [1] Win-rate provenance + Path A/B ablation harness
      → answers "is the spread a metric artifact?" with ZERO accuracy assumptions.
  [2]+[3] Wire Path A behind a flag, SHIPPED TOGETHER WITH the calibrated N-op cost
      → do not trade _n_agree herd bias for the N=0 free-intron exploit.

PHASE 1 — build the REAL oracle (gating cost, WEEKS; the true critical path)
  [O] Procure a genuinely orthogonal CPA truth set (Quant-seq/3'-seq, or a cross-run
      held-out NET-seq sample) + annotated-intron junction truth. Regenerate the absent
      empirical penalty table as part of this.
  *** Until [O] exists, NET-seq concordance is a sanity check, not an oracle. ***

PHASE 2 — cheap, orthogonal, low-risk wins (config + already-ingested data)
  [4] deSALT -x ont2d A/B; mapPacBio maxindel; minimap2 k/w sweep; pt:i prior in 2E walk-back
  [5] mispriming terminal veto for the cDNA path (reuses the existing detector)

PHASE 3 — the ONE terminal CPA solver (HIGH regression risk; depends on [O])
  [6] Unified poly-A-aware terminal boundary solver — choose EITHER the two-state DP
      OR the window re-scoring pass, NOT both. Validate cat1/cat2 shifts EXACTLY before merge.

PHASE 4 — calibrated selection / junction quality (depends on [O] + [2])
  [7] Calibrated junction scorer (yeast PWM/MaxEnt; sklearn-free fit)
  [8] Cross-read junction/CPA consensus + de-herded weights — WITH a hard APA-preservation radius
  [9] Calibrated per-read confidence + abstain

PHASE 5 — ceiling probes / deferred
  GBDT selector (advisory until it beats a wired Path A on [O]); SpliceAI interface stub
  (metazoan); scoped HPC re-seeding (research).
```

**Why this order.** Phase 0 is pure code and settles the metric question that every accuracy claim
depends on. Phase 1 is the unglamorous, expensive truth-set work that the discovery files under-cost;
nothing in Phases 2–5 is *creditable* without it. Phases 2 is safe config wins that need no rewrite.
Phase 3 is the one principled-but-risky structural change. Phase 4 makes selection defensible. Phase 5
is research, gated on everything above.

---

## The single biggest risk, and the environment caveats

**Biggest risk — NET-seq circularity.** The entire program is validated against an oracle (NET-seq)
that is simultaneously (a) the tuning signal for nearly every proposed model/weight/prior, (b) already
partially baked into the corrected positions being scored, and (c) the *only* CPA truth set in the repo.
**Every accuracy claim — including the "deSALT 78.9%" the whole investigation rests on — can be made to
"improve" by construction.** The two-state mitigation: do the zero-oracle Path A/B ablation first (it
needs no truth set), and treat procuring a genuinely orthogonal truth set as the real, weeks-long gating
investment — not the cheap "curate the truth set" line every file assumes.

**Secondary risk — two "low-effort" items are actually high-regression-risk:**
- **Item 6** rewrites the 2B/2E/2G terminal stack, whose 2G-over-2E priority, terminal-D/N stripping,
  and `find_polya_boundary` N-op guards were fixed across many CLAUDE.md releases (v2.9.x–v3.0.x) on
  real data. "Decide the boundary once" is cleaner but is a rewrite, not a cheap add.
- **Item 8** can erase the alternative-polyadenylation heterogeneity that is the product's reason to
  exist (and is self-reinforcing with `_n_agree` herd bias). The HP-noise-radius bound is load-bearing,
  not optional.

**Environment caveats (this checkout, verified by the red-team).** Several proposals understate plumbing
because the tools they assume are missing here:
- `edlib` imports; **`mappy`, `parasail`, `lightgbm`, `scikit-learn`, and `numba` do NOT import**, and
  **`minimap2` is not on `PATH`**. So: Module 2H's Numba DP falls back to pure Python ("only the call
  site is new" hides a perf dependency); the GBDT/sklearn-calibration items (9, Phase 5) are not
  runnable as-is; the minimap2 `--print-seeds` / flag-sweep items need the binary installed first.
- The empirical **`penalty_scores.tsv` / `str_penalty_scores.tsv` and the entire
  `common/scripts/nanopore/` tree are ABSENT.** Any DP proposal that consumes the empirical table must
  **regenerate and commit it first** — and regeneration itself depends on a working aligner panel + a
  labelled run, i.e. it is gated on Phase 1.

---

## If you do only three things

1. **Run the Phase 0 provenance + Path A/B ablation** (item 1). It is days of pure code, needs no
   oracle, and tells you whether the win rates are real or a metric artifact — the question every other
   decision hinges on. Commit the `aligner_summary.tsv` so the numbers stop being a bare assertion.

2. **Wire Path A with the calibrated N-op cost** (items 2 + 3, together). This makes the documented
   metric actually run and closes the code-verified free-intron exploit in the same change — turning a
   popularity vote into a defensible quality score. Keep it flag-gated until item 1's measurement
   confirms it helps.

3. **Procure a genuinely orthogonal CPA truth set** (Phase 1, item O). This is the unglamorous,
   weeks-long investment the discovery files all under-cost, and it is the precondition for *any*
   accuracy claim being more than a hypothesis. Without it, NET-seq is circular and every "we improved
   accuracy" result is unfalsifiable.

*(The cheap config A/Bs in item 4 — especially deSALT `-x ont2d` on the top aligner — are nearly free
and worth doing alongside, but the three above are the load-bearing moves.)*
