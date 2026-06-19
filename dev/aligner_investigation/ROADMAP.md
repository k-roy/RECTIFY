# Improvement Roadmap — Long-Read Aligner SPLICE-JUNCTION Program for RECTIFY

> **⚠️ READ FIRST.** This roadmap was re-led **junction-first** (the project's primary goal) and
> reconciled against the build the team runs (`origin/drs-validation-rebuild`). Two structural
> corrections from [`CORRECTIONS_vs_DRS_BUILD.md`](./CORRECTIONS_vs_DRS_BUILD.md) change the work:
> (1) **Path A (HP-edit-distance) already runs in production** — "wire Path A on" is DONE, only the
> calibrated N-op cost remains open; (2) **the empirical penalty/STR/overhang tables are bundled and
> auto-resolved by `--Scer`** — the "regenerate the absent tables" item becomes "validate / version
> the bundled tables." Both originals survive uncorrected in `04_discovery/`; this roadmap is authoritative.

**What this is.** The forward-looking deliverable. **The primary goal is validating and improving
SPLICE-JUNCTION PLACEMENT**; 3′-end / CPA resolution is **already well-handled** by RECTIFY and is
demoted to a clearly secondary section. This synthesizes the six discovery dossiers in `04_discovery/`
**through the lens of `04_discovery/REDTEAM_proposals.md`** (which fact-checked every code claim against
RECTIFY source, found the redundancies, and ordered the work by dependency) and
**`SPLICE_JUNCTION_PLACEMENT.md`** (the headline junction deliverable). The six files proposed ~30
sub-ideas; after de-duplication they reduce to **~9 distinct work items**. This roadmap lists each
retained item once, with its red-team verdict, then a dependency-ordered phased plan that **leads with
junction-placement work**.

**Tagging.** **FACT** = code/source-verified. **INFERENCE** = reasoned. **HYPOTHESIS** = plausible,
untested on committed ground truth. Effort: **XS / S / M / L**. Verdict (red-team):
**KEEP / KEEP-WITH-CAVEATS / DESCOPE / REJECT**.

---

## ⚠️ The measurement problem comes first (Phase 0) — and be honest about it

Even though the primary goal is junction placement, **the cheapest unblocker is still a
zero-oracle measurement step**, and the red-team (`REDTEAM_proposals.md §3.3, §4`) exposes why
NET-seq alone is not enough:

- **The cheapest real unblocker needs NO oracle at all.** The win-rate **provenance audit + Path A vs
  Path B ablation** is pure code: re-run `merge_corrected_tsvs` on the committed per-aligner TSVs with
  and without per-aligner raw BAMs, confirm whether the cited 79/18/2/0.8/0.1 are post-v3.3.0-fix, and
  commit an `aligner_summary.tsv` + `PROVENANCE.json`. This **settles whether the win-rate spread moves
  between the two live paths, with zero accuracy assumptions.**

  > **⚠️ BUILD CORRECTION (`CORRECTIONS_vs_DRS_BUILD.md` Claim 2 / §A).** This is now an **ablation of
  > two LIVE paths**, not "wire a dead path on." On the build, both production call sites pass
  > per-aligner **raw** BAMs + genome (`single_sample.py:238-250`; generated split body
  > `split_command.py:1085-1094`), and `use_hp_ed = bool(per_aligner_corrected_bams or
  > per_aligner_raw_bams)` (`corrected_consensus.py:1262`) → **Path A (HP-edit-distance) runs in
  > production**; the `_n_agree` popularity sort is the documented **fallback only**. The ablation still
  > has value (does the metric change between Path A and the Path-B fallback?), but the premise "the
  > corrected BAMs Path A needs exist on disk but are never collected" is **superseded** — Path A is
  > already wired via the lazy raw-BAM path. **Do this first; it is days of work.**

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

## The retained items (de-duplicated) — JUNCTION-PLACEMENT FIRST

Each appears multiple times across the six files; consolidated here once. **Junction-placement
items (the primary goal) lead; 3′/CPA and infra items follow in a clearly-secondary section.**
Item numbers are kept stable to the source dossiers/`REDTEAM_proposals.md` (hence non-sequential
here) so cross-references still resolve.

### A. Primary — splice-junction placement

| # | Item | Mechanism (1 line) | Junction failure it fixes | Effort | Risk | Verdict | Key caveat |
|---|---|---|---|---|---|---|---|
| **3** | **Calibrated N-op cost** (do first among junction work) | Charge each intron −logP(novel intron here) instead of 0 in `_cigar_hp_edit_distance` | `N` ops are **cost 0** → an aligner invents a false intron to lower its HP-edit distance and win | S | Low–Med | **KEEP (highest junction priority)** | Must be a *calibrated* cost, not a fixed gate (PERMANENT no-gate policy). ⚠ Path A (which consumes this score) **already runs** — see build correction below; and a junction-anchor gate already backstops free-N on human only |
| **8j** | **Cross-read junction consensus + de-herded weights** | Pool every `N`-op across all 5 aligner BAMs + annotation into a per-locus consensus; re-snap with `_score_junction` as a per-read veto; weight agreement by reliability, not raw count | Per-read junction inhomogeneity (the trait deSALT wins on) — now available to all five, de-biased | M | **High** | **KEEP-WITH-CAVEATS** | **Over-snapping collapses genuine alt-5′/3′SS isoforms <50 bp apart (real in yeast RP genes)** — snap only within `max_boundary_shift=50`; per-read sequence veto mandatory; weights empirical, tiebreaker-only |
| **7** | **Calibrated junction scorer + splice-strength (MaxEnt/PWM)** | Replace 2H's brittle lexicographic tuple with a log-odds: HP-seq support + yeast MaxEnt/PWM **splice strength** + cross-read support + small annotation prior | Hand-tuned tier ordering is brittle/non-comparable; RECTIFY has **NO positive splice-strength term** (only canonical tiers) | M | Med | **KEEP-WITH-CAVEATS** | MaxEnt is metazoan-trained → retrain a yeast PWM; `w_ann` small (annotation circularity); build behind a swappable `splice_strength()` interface |
| **2t** | **Validate / version the bundled penalty + overhang tables** | Confirm the auto-resolved `penalty_scores.tsv` / `str_penalty_scores.tsv` / `junction_overhang_table.tsv` are correct for the run; pin versions + provenance | Mis-scored HP-deletion junctions (dominant ONT boundary error) under a stale or wrong-protocol table | S | Low | **KEEP (re-scoped)** | ⚠ Tables are **bundled and auto-resolved by `--Scer`**, NOT absent — see build correction. Validate, do not "regenerate"; R10.4.1 + *S. cerevisiae*-specific, must not transfer to HiFi/other organisms |
| **5j** | **De-novo micro-exon recovery (≤~6–10 nt)** | GMAP-style: probe a ≤10 nt exon inside an `N`-op when flanking consensus exons leave unexplained boundary error; accept only if both splice sites clear a MaxEnt gate (item 7) AND cross-read support (item 8j) ≥ N | Micro-exon misses (de-novo, works where uLTRA's annotation-injection can't) | M | Med–High | **DESCOPE → metazoan** | Yeast has few micro-exons; **no validation read exercises this** (needs a cat10 + SIRV control); double gate essential; defer to metazoan port |

### B. Secondary — 3′/CPA (already-handled) + selection infrastructure

3′-end/CPA placement is **already well-handled** by RECTIFY's terminal stack and is **not** the
primary goal. These items are retained for completeness but are **lower priority** than Section A.

| # | Item | Mechanism (1 line) | Failure mode it fixes | Effort | Risk | Verdict | Key caveat |
|---|---|---|---|---|---|---|---|
| **1** | **Win-rate provenance + Path A/B ablation harness** | Re-run the merge both ways on committed TSVs; commit `aligner_summary.tsv` + `PROVENANCE.json` | Win rates are un-committed and possibly a metric artifact | XS | Low | **KEEP (do first, zero-oracle)** | Needs no oracle. ⚠ Now an ablation of **two live paths** — Path A already runs (build correction) |
| **2** | **Wire Path A** (HP-edit-distance selection) | — | Production allegedly sorts on `_n_agree` popularity, not corrected quality | — | — | **DONE ON BUILD** | ⚠ **Already wired** via the lazy raw-BAM path (`corrected_consensus.py:1262`; `single_sample.py:238-250`; `split_command.py:1085-1094`). Only the calibrated N-op cost (item 3) remains open. See build correction |
| **6** | **ONE unified poly-A-aware terminal CPA solver** | A two-state (genomic→tail) terminal DP OR a window re-scoring pass — the CPA is the body→tail transition | Every aligner soft-clips / slips the 3′ end into a genomic A-run | M | **High** | **DEMOTED — 3′/CPA already handled** | A 3′-end solver, **not a junction improvement**; rewrites the most heavily-patched code path; validate cat1/cat2 shifts exactly if ever attempted |
| **4** | **Cheap config A/Bs** | deSALT `-x ont2d`; minimap2 `-k`/`-w` & `--end-bonus` sweep; `pt:i` prior in 2E | Possible free sensitivity on the top aligner; HP 3′ drift | XS–S | Low | **KEEP (mixed)** | ⚠ mapPacBio `maxindel` is **already set** on the build (item 4e correction) — drop it. `pt:i` prior is a 3′ signal (demoted). deSALT `-x ont2d` (junction-relevant) is the keeper |
| **5** | **Mispriming terminal veto (cDNA)** | Reuse the genomic-downstream-A detector as a bounded terminal-placement penalty under `--dT-primed-cDNA` | Oligo-dT false CPAs on genomic A-runs (QuantSeq) | S | Low | **DEMOTED — 3′/CPA** | A cDNA 3′-end artifact, not a junction mechanism; DRS-exempt; penalty/flag, never a hard gate |
| **9** | **Calibrated per-read confidence + abstain** | Emit a calibrated `selection_confidence`; flag coin-flip reads so downstream APA can filter | `confidence` is a self-assessed 3-bucket flag; contested calls pollute APA/DESeq2 | M | Low–Med | **KEEP-WITH-CAVEATS** | Calibration target is NET-seq → inherits the circularity; `sklearn.isotonic` not importable in this checkout |

> **⚠️ BUILD CORRECTION (`CORRECTIONS_vs_DRS_BUILD.md`).** Three items above changed materially on
> `origin/drs-validation-rebuild`: **(Claim 2/§A)** item 2 "wire Path A" is **DONE** — both call sites
> pass raw BAMs + genome and `use_hp_ed=True` runs in production; the legacy `_n_agree` popularity sort
> is the **fallback only**. **(Claim 1/§B)** the empirical penalty/STR/overhang tables are **bundled**
> (`rectify/data/genomes/*/penalty_tables/`) and **auto-resolved by `--Scer`** (protocol-routed,
> `data/__init__.py:1188-1208`, hooked at `cli.py:198-199`) — item 2t is "validate," not "regenerate."
> **(Claim 3/§D)** the **N-op-cost-0** exploit (item 3) now has a **junction-anchor gate backstop**
> (`_cigar_min_junction_anchor` + `_add_chimera_flag`), default **10 bp on human, 0/off on yeast**
> (`data/__init__.py:1053-1062`) — so on yeast the overhang filter is still the only structural defense
> and item 3 retains full priority; on human the gate already adds a graded defense.

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
recovery** (item 5j; low yeast value, no validation read exercises it — defer to the metazoan port).

---

## Dependency-ordered phased plan

```
PHASE 0 — measurement, no oracle needed (DO FIRST, days)
  [1] Win-rate provenance + Path A/B ablation harness
      → answers "does the spread move between the two LIVE paths?" with ZERO accuracy assumptions.
      ⚠ Path A already runs in production (build correction) — this is an ablation, not a wiring task.

PHASE 1 — build the REAL junction oracle (gating cost, WEEKS; the true critical path)
  [O] Procure a genuinely orthogonal JUNCTION truth set — short-read RNA-seq split-read
      junctions (STAR SJ.out.tab) on the SAME sample, and/or a SIRV/sequin control with
      documented multi-exon transcripts + micro-exons — PLUS curated SGD intron coords.
      (A 3′/CPA truth set is secondary — CPA is already handled.)
  [2t] Validate / version the BUNDLED empirical penalty + overhang tables (NOT regenerate —
       they are bundled and auto-resolved by --Scer; build correction).
  *** Until [O] exists, ≥4-aligner concordance + cat3/cat7/cat9 are sanity checks, not an oracle. ***

PHASE 2 — junction-placement core (PRIMARY; depends on [1], partly on [O])
  [3]  Calibrated N-op cost in _cigar_hp_edit_distance — charge −logP(novel intron here) instead of 0,
       feeding the Path-A HP-edit-distance score that ALREADY runs. Highest junction priority.
       Validate cat7/cat9 do not regress before default-on.
  [8j] Cross-read junction consensus across all 5 aligners + annotation, de-biased by a calibrated
       support score, protected by a per-read _score_junction veto + a hard alt-SS-preservation radius.
  [7]  Calibrated junction scorer with a yeast MaxEnt/PWM splice-strength term behind a swappable
       interface; replaces 2H's lexicographic tuple with a comparable log-odds. (sklearn-free fit.)

PHASE 3 — cheap, orthogonal, low-risk wins (config + already-ingested data)
  [4]  deSALT -x ont2d A/B (junction-relevant on the top aligner); minimap2 k/w sweep.
       (Drop mapPacBio maxindel — already set on build; pt:i / --end-bonus are 3′ signals, demoted.)
  [9]  Calibrated per-read selection confidence + abstain.

PHASE 4 — SECONDARY: 3′/CPA (already handled) + research ceilings
  [6]  Unified poly-A-aware terminal CPA solver — DEMOTED (3′/CPA already well-handled; HIGH
       regression risk; only if a concrete CPA gap is shown). [5] mispriming terminal veto (cDNA 3′).
  [5j] De-novo micro-exon recovery — DESCOPE to metazoan port (no yeast validation read).
       GBDT selector (advisory); SpliceAI interface stub (metazoan); scoped HPC re-seeding (research).
```

**Why this order.** Phase 0 is pure code and settles whether the win-rate spread is metric-sensitive
(now an ablation of two live paths, not a wiring task). Phase 1 builds the **junction** oracle — the
unglamorous, expensive truth-set work the discovery files under-cost (and validates the already-bundled
tables rather than regenerating absent ones). Phase 2 is the **primary junction-placement program**:
the calibrated N-op cost feeding the live Path-A score, cross-read junction consensus, and a calibrated
splice-strength scorer. Phase 3 is safe config + confidence wins. Phase 4 is explicitly **secondary** —
3′/CPA (already handled), micro-exon (descoped to metazoan), and research ceilings.

---

## The single biggest risk, and the environment caveats

**Biggest risk — oracle circularity for the JUNCTION metric.** The junction program must be validated
against a signal that is independent of the aligners being judged. The two available internal signals —
≥4-aligner concordance and the cat3/cat7/cat9 reads — are **sanity checks, not oracles** (concordance is
partly herd-defined; 3 of 5 aligners share minimap2 lineage). For the **3′/CPA** side (secondary here),
NET-seq is simultaneously (a) the tuning signal, (b) partially baked into the corrected positions, and
(c) the only CPA truth set — so CPA accuracy claims can "improve" by construction. Mitigation: do the
zero-oracle Path A/B ablation first, then procure a genuinely orthogonal **junction** truth set
(short-read split-read junctions / SIRV) as the real, weeks-long gating investment — not the cheap
"curate the truth set" line every file assumes.

**Secondary risk — high-regression-risk items:**
- **Item 8j** (cross-read junction consensus) can erase genuine alt-5′/3′SS isoform heterogeneity
  (real in yeast RP genes, <50 bp apart) and is self-reinforcing with any agreement-sensitive
  tiebreaker. The `max_boundary_shift=50` radius + per-read `_score_junction` veto are load-bearing,
  not optional.
- **Item 6** (the demoted terminal CPA solver) rewrites the 2B/2E/2G terminal stack, whose
  2G-over-2E priority, terminal-D/N stripping, and `find_polya_boundary` N-op guards were fixed across
  many CLAUDE.md releases (v2.9.x–v3.0.x) on real data. It is a rewrite, not a cheap add — and 3′/CPA is
  already handled, so it stays demoted.

**Environment caveats (this checkout, verified by the red-team).** Several proposals understate plumbing
because the tools they assume are missing here:
- `edlib` imports; **`mappy`, `parasail`, `lightgbm`, `scikit-learn`, and `numba` do NOT import**, and
  **`minimap2` is not on `PATH`**. So: Module 2H's Numba DP falls back to pure Python ("only the call
  site is new" hides a perf dependency); the GBDT/sklearn-calibration items (9) are not runnable as-is;
  the minimap2 `--print-seeds` / flag-sweep items need the binary installed first.
- > **⚠️ BUILD CORRECTION (`CORRECTIONS_vs_DRS_BUILD.md` Claim 1/§B).** The old caveat that
  > `penalty_scores.tsv` / `str_penalty_scores.tsv` and the `common/scripts/nanopore/` tree are
  > **ABSENT** describes only the *master* checkout. **On the build they are bundled** at
  > `rectify/data/genomes/{saccharomyces_cerevisiae,homo_sapiens}/penalty_tables/` (with protocol
  > variants) and **auto-resolved by `--Scer`** (`data/__init__.py:1188-1208`; `cli.py:198-199`). Any
  > DP proposal consuming the empirical table should **validate/version the bundled table** (item 2t),
  > not regenerate an absent one. (Numba availability remains environment-dependent — unchanged.)

---

## If you do only three things (junction-placement first)

1. **Fix the N-op-cost-0 free-intron exploit + cross-read junction consensus** (items 3 + 8j). This is
   the primary-goal core: charge each `N`-op a *calibrated* −logP(novel intron here) so it feeds the
   **Path-A HP-edit-distance score that already runs** (no wiring needed — see build correction), and
   pool junctions across all 5 aligners into a de-biased consensus re-snapped under a per-read sequence
   veto. Together these make **junction quality**, not output style, decide placement. Validate
   cat7/cat9 do not regress before default-on. *(Note: "wire Path A on" is DONE on the build; only the
   N-op cost remains open.)*

2. **Procure a genuinely orthogonal JUNCTION truth set** (Phase 1, item O). Short-read RNA-seq
   split-read junctions (STAR `SJ.out.tab`) on the same sample, and/or a SIRV/sequin control with
   documented multi-exon transcripts + micro-exons. This is the unglamorous, weeks-long investment the
   discovery files under-cost, and it is the precondition for *any* junction-accuracy claim being more
   than a hypothesis. ≥4-aligner concordance and cat3/cat7/cat9 are sanity checks, not an oracle —
   never validate against internal agreement alone.

3. **Run the zero-oracle Path A/B ablation, then add the splice-strength signal** (items 1, then 7).
   Commit `aligner_summary.tsv` + `PROVENANCE.json` from re-running the merge through both live paths
   (settles whether the win-rate spread is metric-sensitive). Then build a yeast MaxEnt/PWM
   donor/acceptor model behind a swappable interface and fold it — plus the **already-bundled** empirical
   penalty table (validate, don't regenerate) — into a single calibrated junction scorer that replaces
   2H's brittle lexicographic tuple.

*(The cheap config A/Bs in item 4 — especially deSALT `-x ont2d` on the top aligner — are nearly free
and junction-relevant, worth doing alongside; mapPacBio `maxindel` is already set on the build, so drop
it. 3′/CPA items 5/6 are secondary — CPA is already well-handled.)*
