# Long-Read Aligner Investigation — Index & Executive Overview

> **⚠️ READ FIRST — build corrections.** Several load-bearing findings in this
> investigation were verified against a `master`-derived working tree and are **now wrong on
> the build the team actually runs** (`origin/drs-validation-rebuild`). The authoritative
> reconciliation — with build `file:line` evidence — is
> [`CORRECTIONS_vs_DRS_BUILD.md`](./CORRECTIONS_vs_DRS_BUILD.md). The most affected claims
> ("HP-edit-distance never runs in production"; "empirical penalty tables absent") have been
> **retracted/softened inline** below and in `ROADMAP.md`; the originals survive only in the
> unedited `03_adversarial/` and `04_discovery/` docs (cross-reference that file when reading them).

**What this is.** A multi-agent investigation into the long-read RNA-seq aligners RECTIFY
uses (and one it doesn't) for ONT direct-RNA / cDNA *S. cerevisiae* alignment. **The PRIMARY
objective is validating and improving SPLICE-JUNCTION PLACEMENT** — where each aligner lands
the `N` (intron-skip) op, and whether RECTIFY's refinement/selection picks the
sequence-supported boundary. **3′-end (cleavage-and-polyadenylation site, CPA) resolution is
treated as already well-handled by RECTIFY** (the 2B/2E/2G terminal stack + NET-seq
refinement) and is **secondary** throughout — it appears mainly where it interacts with
junction placement. The headline deliverable is **[`SPLICE_JUNCTION_PLACEMENT.md`](./SPLICE_JUNCTION_PLACEMENT.md)**.
The work proceeds in five layers: per-aligner source-level dossiers → synthesis → adversarial
red-team → a forward-looking improvement program (consolidated in `ROADMAP.md`). Every
non-trivial claim is tagged **FACT** (source/paper/code-verified), **INFERENCE** (reasoned),
or **HYPOTHESIS** (plausible, untested).

---

## Directory map

```
dev/aligner_investigation/
├── README.md                         ← you are here (index + executive overview)
├── SPLICE_JUNCTION_PLACEMENT.md      ← ★ HEADLINE DELIVERABLE: junction-placement validation & improvement (PRIMARY goal)
├── CORRECTIONS_vs_DRS_BUILD.md       ← ⚠ build re-verification (origin/drs-validation-rebuild); retracts/softens several findings
├── ROADMAP.md                        ← prioritized, dependency-ordered improvement plan (junction-first)
├── 01_investigation/                 ← per-aligner source-level dossiers
│   ├── minimap2.md                   minimizer seed-chain-align; the ONT-DRS baseline
│   ├── deSALT.md                     RdBG de-Bruijn index + two-pass cross-read junction pooling (top winner)
│   ├── uLTRA.md                      annotation-segment colinear chaining; micro-exon specialist
│   ├── gapmm2.md                     minimap2 wrapper + edlib terminal-exon refinement (PAF-only)
│   ├── mapPacBio_bbmap.md            BBMap full affine-SW seed-and-extend; best raw 3′ boundary
│   ├── gmap.md                       sandwich-DP splice-model-first aligner (CONTRAST — not in panel)
│   ├── pacbio_ecosystem.md           pbmm2/Iso-Seq landscape; poly-A trimmed BEFORE alignment
│   └── ont_drs_ecosystem.md          ONT error model, 3′-end ambiguity, the RECTIFY niche
├── 02_synthesis/
│   ├── COMPARISON.md                 master comparison table; what they share vs differ; win-rate mapping
│   └── DEEP_DIVE.md                  function/equation-level walkthroughs; 3′-end failure-mode crux
├── 03_adversarial/                   ← skeptical red-team of the above
│   ├── redteam_denovo.md             fact-checks minimap2/deSALT/gapmm2/mapPacBio dossiers
│   ├── redteam_annotation_ecosystem.md  fact-checks uLTRA/pacbio/ont_drs dossiers
│   └── redteam_winrates_selection.md ⚠ the load-bearing audit: win rates ≠ accuracy ranking
└── 04_discovery/                     ← improvement proposals (synthesized in ROADMAP.md)
    ├── improve_seeding_indexing.md   syncmer/HPC/strobemer seeding; CPA-anchor priors
    ├── improve_chaining_dp.md        poly-A-aware terminal DP; HP-aware realign; error-model gaps
    ├── improve_splice_junctions.md   cross-read junction consensus; calibrated scorer; N-op cost  (⚠ Path-A "wiring" already DONE on build — see CORRECTIONS)
    ├── improve_3prime_cpa.md         unified terminal CPA solver; pt:i fusion; mispriming veto  (secondary — 3′/CPA already handled)
    ├── improve_ensemble_selection.md de-herd _n_agree; reliability weights; abstain  (⚠ "wire Path A" item already DONE on build — see CORRECTIONS)
    ├── improve_ml_learned.md         GBDT selector; signal/dwell fusion; honest "where ML is overkill"
    └── REDTEAM_proposals.md          adversarial review of the 6 dossiers → 9 de-duped items + roadmap
```

## Reading order for a newcomer

1. **`SPLICE_JUNCTION_PLACEMENT.md`** — the headline deliverable and the project's primary
   goal: how each aligner places junctions, how RECTIFY refines/selects them, and the
   junction-placement validation strategy. Start here.
2. **`CORRECTIONS_vs_DRS_BUILD.md`** — the build re-verification. Read before trusting any
   "HP-edit-distance / penalty-table" claim in the older red-team/discovery docs.
3. **`02_synthesis/COMPARISON.md`** — the master table and the share/differ analysis.
4. **`03_adversarial/redteam_winrates_selection.md`** — read before you believe any win-rate
   claim (but apply the corrections in #2 — its "Path B is what production runs" conclusion is
   retracted on the build).
5. **`02_synthesis/DEEP_DIVE.md` §7** — the 3′-end / CPA failure-mode crux (secondary; explains
   why an ensemble exists).
6. Per-aligner dossiers in `01_investigation/` as needed (start with `deSALT.md`, `mapPacBio_bbmap.md`,
   `minimap2.md`).
7. **`04_discovery/REDTEAM_proposals.md`** then **`ROADMAP.md`** for what to do next.

---

# EXECUTIVE OVERVIEW

## Primary goal: splice-junction placement

The deliverable this investigation exists to support is **base-accurate, sequence-supported
splice-junction placement** — getting the `N`-op donor/acceptor boundary onto the true intron
edge, and getting RECTIFY's ensemble selection to pick the aligner that did so. **3′-end / CPA
resolution is already well-handled** by RECTIFY's terminal stack (Modules 2B/2E/2G +
`find_polya_boundary` + NET-seq refinement, hardened across v2.9.x–v3.0.x on real data) and is
**secondary** here. The aligner-mechanism discussion below is retained because the panel's
junction *and* 3′-end behaviour are entangled in the same CIGAR — but the forward-looking
program (`ROADMAP.md`) now **leads with junction-placement work** and demotes CPA-specific
items. See **`SPLICE_JUNCTION_PLACEMENT.md`** for the full junction analysis.

## The panel

RECTIFY runs a **5-aligner correct-first ensemble** on ONT DRS/cDNA yeast reads:
**minimap2, deSALT, uLTRA, gapmm2, mapPacBio (BBMap)**. A sixth, **GMAP**, was investigated as an
**algorithmic contrast / candidate only** — it is *not* in the panel and has no win-rate.

All six share the **seed → chain → align** paradigm: find short exact/approximate anchors (seeds),
order collinear anchors into a chain, then fill gaps with base-level dynamic programming. They differ
in *which stage they specialize* and *how they place the intron and the 3′ end* (FACT, verified across
the dossiers + `multi_aligner.py`):

- **minimap2** — (w,k)-minimizer seeds + sparse-DP chaining + KSW2 splice DP with a dedicated intron
  state (free intron extension) and a GT-AG donor/acceptor prior. The de-facto ONT-DRS baseline.
  *Weakness:* `--end-bonus 0` → soft-clips poly-A, `reference_end` drifts; HP indel noise placed
  arbitrarily in A-runs.
- **deSALT** — a **Reduced de-Bruijn-graph (RdBG)** index (one lookup returns every genomic copy of a
  unipath → breadth of candidates) + a **two-pass** strategy whose Pass-2 **integrates the alignment
  skeletons of *all* reads in the batch** to infer consensus exon boundaries (FACT-PAPER, verified:
  "integrates all the alignment skeletons of all the reads"). This **cross-read junction homogeneity**
  is its unique structural property. *Weakness:* no poly-A awareness; engineering-fragile in RECTIFY.
- **mapPacBio (BBMap)** — an all-k-mer index + multi-k-mer-voting seed-and-extend + a **full
  affine-transform global DP** (match +100 / mismatch −127) of every candidate window → the
  **column-optimal terminal boundary** = best *raw* 3′-end precision in the panel. *Weakness:* no
  GT-AG/splice-graph model (introns found as scored gaps); 6 kb read ceiling; JVM-slow.
- **uLTRA** — an **annotation-segment** index (GTF → "parts" → "segments") + two-pass colinear chaining
  that **injects all annotated exons of the candidate gene** (recovers micro-exons). *Weakness:*
  annotation-bound; strips poly-A (no CPA); falls back to minimap2 on novel structure.
- **gapmm2** — a thin **minimap2 wrapper** + an **edlib** terminal-exon refinement pass; PAF-only
  output (RECTIFY reconstructs the BAM). *Weakness:* its body alignment *is* minimap2, so it rarely
  adds anything distinctive.
- **GMAP (contrast)** — sampled-oligomer seeds + **sandwich DP** that brackets the intron between two
  anchored exon ends and chooses the boundary jointly maximizing base + splice score — the cleanest
  external example of **splice-model-first** placement, and the closest analogue to RECTIFY's own
  sequence-first Module 2H. *Not panel-fit:* ~30× slower than minimap2, no 3′-end benefit.

## The reported win rates — and the caveat that must come first

CLAUDE.md reports correct-first win rates of **deSALT 78.9% · mapPacBio 18.2% · uLTRA 2% ·
gapmm2 0.8% · minimap2 0.1%**, and frames them as biology. **The adversarial review
(`redteam_winrates_selection.md`) shows this framing is not safe** — but note one of its
central planks has been **retracted on the build** (see callout):

> **⚠️ BUILD CORRECTION (see `CORRECTIONS_vs_DRS_BUILD.md` Claim 2 / §A).** The original
> red-team conclusion — *"selection is a legacy popularity vote; the HP-edit-distance path never
> runs; the win rates were produced by Path B"* — is **OUTDATED-FALSE on
> `origin/drs-validation-rebuild`**. Both production merge call sites now pass per-aligner **raw**
> BAMs + genome (`single_sample.py:238-250`; generated split body `split_command.py:1085-1094`),
> and `merge_corrected_tsvs` gained a lazy raw-BAM path that sets
> `use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)`
> (`corrected_consensus.py:1262`) → **Path A (HP-edit-distance) runs in production.** The legacy
> 5-level `_n_agree` popularity sort (Path B) still exists but is the **documented FALLBACK only**
> (no BAMs/genome, or BAM staging fails). The "lineage" argument below is therefore retracted; the
> single-dataset / un-committed caution still stands.

- **These are NOT a clean accuracy ranking.** They are the output of a **single dataset**
  (`wt_by4742_rep1`, one R10.4.1 yeast DRS run) with **no committed provenance** — no script,
  count table, or verified `aligner_summary.tsv` artifact in the repo derives these exact figures.
  (FACT, code-verified.)
- **The selection metric on the build is HP-edit-distance (Path A), not the popularity vote** —
  so the older "produced by the legacy sort" lineage argument **no longer applies to a fresh build
  run** (retracted; see callout). The figures remain a **single-dataset, un-committed metric
  output**, and whether they are pre- or post-v3.3.0-fix is still unconfirmed (an `index_col`
  parsing bug corrupted minimap2/gapmm2/uLTRA read_ids and a `_pt:i:N` normalization bug let
  mapPacBio win isolated groups — both fixed in v3.3.0, but the lineage of these specific rates is
  unverified).
- **deSALT's lead may still reflect output style, not raw 3′ accuracy** (HYPOTHESIS). Its
  cross-read homogeneity makes other aligners' junctions agree, which can advantage it under any
  agreement-sensitive tiebreaker; the deSALT dossier calls its *raw* 3′ ends "imprecise," so the
  win is unlikely to be raw-3′-accuracy-driven. This is now a hypothesis about **junction
  homogeneity**, not an attribution to `_n_agree` specifically.

**Bottom line:** the *mechanisms* in the dossiers (RdBG breadth, cross-read pooling, full-SW boundary,
annotation snapping) are real and well-sourced; the *attribution of a specific win-rate percentage to
a specific mechanism* remains **untested hypothesis** on a single dataset. Treat 79/18/2/0.8/0.1 as
"won an HP-edit-distance sort on one yeast sample, un-committed," never as a validated accuracy ranking.

## Master comparison table (carried from COMPARISON.md)

| Axis | minimap2 | deSALT | uLTRA | gapmm2 | mapPacBio (BBMap) | GMAP *(contrast — not in panel)* |
|---|---|---|---|---|---|---|
| **Seeding** | (w,k)-minimizers; w=5, RECTIFY k=14 | MEMs from RdBG unipath hits | NAMs via strobemers (namfinder) | inherits minimap2 minimizers | all-k-mer index (k=13); multi-k-mer voting | sampled oligomer hash, k=15, stride-3 |
| **Index** | minimizer-hash → (pos,strand) | **Reduced de-Bruijn graph** of unipaths (k=22) | annotation segmentation (GTF→parts→segments) | minimap2 `.mmi` | k-mer → site list; disk-cached, shareable | gamma-pointer hash + IIT splice/SNP files |
| **Chaining** | sparse DP, h=50; 2021 RMQ | SDP over Match Blocks → skeleton | two-pass colinear (RMaxQST trees, overlaps+gaps) | none (minimap2's) | multi-k-mer hit voting | oligomer chaining + diagonalization |
| **Base DP** | KSW2 SIMD splice DP | ksw2-style SIMD banded SW | parasail + edlib | edlib (terminal only) | **full affine-transform global DP** (+100/−127) | **sandwich DP** (both-ends) + nucleotide DP |
| **Splice model** | intron state, two-piece affine, GT-AG prior | two-tier gap; GT-AG via `-R 9`; **cross-read exon pooling** | exact annotated junctions; novel = gap>10 bp | synthesizes `~gt..ag` tokens | **none** (scored gaps → D→N at intronlen≥10; explicit `maxindel=max(200000,max_intron)`) | **sequence-based splice model** (joint base+splice optimum) |
| **3′-end / poly-A** | `--end-bonus 0` → soft-clips; drifts | no poly-A awareness; slips into A-run | **strips poly-A** → no CPA | inherits minimap2 | **column-optimal boundary** (best raw CPA) | no poly-A awareness; slips into A-run |
| **Annotation** | optional (`--junc-bed`) | GTF-free de novo | **REQUIRED** | none | none | not required (de novo) |
| **Output** | SAM | SAM (→ dedup) | ultra.sam | **PAF + cs only** | SAM with `=`/`X` | SAM/GFF/etc. |
| **Speed / memory** | **fastest** | slower; RdBG ~35 GB on human, tiny on yeast | heaviest/brittle | cheap marginal | slow (JVM); RAM-hungry; 6 kb ceiling | **~30× slower than minimap2** (disqualifying) |
| **Determinism** | deterministic | needs post-hoc dedup | deterministic | name-sort; ~0.02% dropped | `ambiguous=best`; emits unmapped dups | deterministic per read |
| **Win rate** | 0.1% | 78.9% | 2% | 0.8% | 18.2% | n/a (contrast only) |

**What they share:** the seed-chain-align skeleton; a sparse-DP chaining lineage (minimap2/deSALT/uLTRA);
the **two-piece affine gap cost** that makes a long intron cheap (the shared "intron trick");
a **GT-AG canonical prior used as a soft signal, never a hard gate** (4 of 5; mapPacBio is the
exception); and RECTIFY-side `=`/`X` CIGAR normalization via `calmd`. **What genuinely differs:**
deSALT's *graph* index and its *cross-read* junction pooling (its exclusive mechanism); uLTRA's
*annotation-required* segmentation; mapPacBio's *full-SW* (not heuristic) boundary; gapmm2's
*wrapper* status; and the splice-placement *philosophy* (minimap2 gap-cost-first, uLTRA
annotation-first, deSALT cross-read-consensus-first, mapPacBio score-first, GMAP splice-model-first).

## What the adversarial review changed

Four corrections materially change how the dossiers should be read:

1. **Win rates are single-dataset and un-committed — treat the specific percentages as hypothesis,
   not biology.** Every "aligner X wins because it is more accurate" claim is downgraded to **untested
   hypothesis** (one yeast sample, no committed provenance), and several (deSALT junction-homogeneity
   advantage, uLTRA annotation-circularity) are **plausibly biased toward output style, not accuracy**.
   *(redteam_winrates)*
   > **⚠️ BUILD CORRECTION (`CORRECTIONS_vs_DRS_BUILD.md` Claim 2 / §A).** The red-team's *"selection
   > runs a legacy popularity sort; the HP-edit-distance / corrected-3′-end path does not even run"*
   > conclusion is **retracted on the build**: both call sites pass per-aligner raw BAMs + genome →
   > `use_hp_ed=True` (`corrected_consensus.py:1262`) → **Path A (HP-edit-distance) runs in
   > production**; the `_n_agree` popularity sort is the **fallback only**.
2. **minimap2 `--splice-flank=no` is "compatibility," not "3′ accuracy."** The live code comment is
   `# Disable for compatibility` (`multi_aligner.py:252`); CLAUDE.md's "important for 3′ end accuracy"
   rationale is asserted, not the author's stated reason. The elaborate mechanistic story built on it
   is hypothesis. *(redteam_denovo M7)*
3. **The "8-nt micro-exon" headline is unsupported.** The uLTRA paper bins small-exon accuracy
   (~60% at ≤10 nt); the canonical tiny-exon failure case (minimap2 issue #253) is **6 nt**. No source
   supports a specific "8 nt" example. *(redteam_annotation U9)*
4. **The AT>CG deletion-rate claim is RECTIFY-internal and circular.** External literature confirms A/T
   homopolymers are *basecalled more accurately overall* than C/G; it does **not** corroborate the
   AT>CG *deletion-specific* direction in the penalty table. The "ratcheting-speed" mechanism is an
   unproven hypothesis, and the dossier uses RECTIFY's own table to validate RECTIFY's table.
   *(redteam_annotation O6)*

One lower-tier flag worth carrying: **deSALT runs no `-x` ONT preset** (the `null` ~13% model, never
A/B-tested against `-x ont2d`).

> **⚠️ BUILD CORRECTION (`CORRECTIONS_vs_DRS_BUILD.md` Claim 4e / §C).** The old "mapPacBio sets no
> `maxindel`, relies on BBMap's soft ~16000 default" flag is **stale on the build**: mapPacBio now
> runs `intronlen=10` (not 50) **and** an explicit `maxindel=max(200000, max_intron)`
> (`multi_aligner.py:749,754`) — the previously-missing cap is set. The `redteam_denovo` B5 concern
> is RESOLVED.

---

*FACT claims trace to the dossiers, the cited papers/repos, or `rectify/core/align/multi_aligner.py`
and `corrected_consensus.py` (read directly). INFERENCE/HYPOTHESIS claims — above all the win-rate
causal chains — are synthesis, not assertions by the tool authors, and are explicitly flagged as such.
Code claims marked with a **⚠ BUILD CORRECTION** callout were re-verified against
`origin/drs-validation-rebuild` in `CORRECTIONS_vs_DRS_BUILD.md`; where the build differs from the
older `master`-derived dossiers, the build is authoritative.*
