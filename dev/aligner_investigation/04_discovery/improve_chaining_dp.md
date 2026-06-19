# Improving Chaining & Base-Level DP for bp-Accurate 3′-End Alignment in RECTIFY

> Verified against `origin/drs-validation-rebuild` @ 366c885 (2026-06-19).

**Goal.** Propose concrete, rankable improvements to the *state of the art* in chaining and
base-level dynamic programming for long-read RNA-seq alignment, specifically to win the metric
RECTIFY actually cares about: **bp-accurate 3′ end / cleavage-and-polyadenylation site (CPA)** on
ONT DRS/cDNA (and PacBio), where every panel aligner currently fails in a characteristic way.

**Grounding (verified, not assumed).** Read against
`02_synthesis/{COMPARISON.md,DEEP_DIVE.md}`, the five dossiers, and the three `03_adversarial/`
reports. The proposals below respect the verified caveats:

- minimap2's `--splice-flank=no` carries the live source comment `# Disable for compatibility`
  (`multi_aligner.py:252`); it is **not** a proven 3′-accuracy mechanism (redteam_denovo M7).
- Win rates (deSALT 78.9 / mapPacBio 18.2 / uLTRA 2 / gapmm2 0.8 / minimap2 0.1) are
  **single-dataset and un-committed**. They were produced by the live **HP-edit-distance selector**
  (`use_hp_ed=True`; both merge call sites pass per-aligner BAMs + genome), not by the legacy
  popularity sort — but the "aligner X is more accurate because it wins" story is still an
  **untested hypothesis** on one yeast DRS sample. Proposals here therefore (a) do **not** rely on
  the causal story, and (b) include a ground-truth concordance harness so improvements are measured
  against truth, not against the selection metric.
- deSALT runs with **no `-x` ONT preset** (`null` ≈13% divergence model) — possibly mistuned, never
  A/B-tested (COMPARISON §6.1).
- mapPacBio sets `intronlen=10` (the minimum deletion length relabeled D→N) and an explicit
  `maxindel=max(200000, max_intron)`, so it does **not** rely on BBMap's soft 16000 default — the
  search bound is sized from `max_intron` (multi_aligner.py:749,754).

**What RECTIFY already has** (verified on disk — this drives feasibility):
- `rectify/core/splice/hp_penalty.py::HpPenaltyTable` — parses the AT/CG `penalty_scores.tsv`
  (`del_cost(hp, ref_base)`, `ins_cost`, STR variants), with legacy + base-class formats.
- `hp_penalty._score_hp_dp_numba` — a **Numba-compiled rolling NW semi-global DP** with
  per-reference-position `del_costs` and per-query `ins_costs` arrays (left end fixed, right suffix
  free). This is *already* an error-model-aware DP kernel — it is the engine of Module 2H.
- `rectify/core/align/local_aligner.py` — affine-gap (Gotoh 1982) local aligner used for Cat3 exon
  CIGAR.
- `rectify/core/splice/junction_refiner.py` (Module 2H) — post-consensus N-op refinement.
- `edlib` is importable; **`mappy` and `parasail` are NOT currently importable** in this env
  (matters for "patch the aligner" vs "shell out" feasibility).
- The empirical penalty / STR / overhang tables **are bundled** at
  `rectify/data/genomes/{saccharomyces_cerevisiae,homo_sapiens}/penalty_tables/` (with protocol
  variants — `penalty_scores.tsv`, `penalty_scores_cdna.tsv`, `str_penalty_scores.tsv`,
  `junction_overhang_table.tsv`) and are **auto-resolved by `--Scer`**. Any proposal consuming
  them does not need to regenerate them; the open work is validating/versioning them (see P0).

---

## P0 (prerequisite, not a proposal): make the metric measurable

One cheap blocker gates everything below; do it first or the impact numbers are unmeasurable.

1. **Commit a ground-truth 3′-end harness.** Build a held-out CPA truth set (Quant-seq / 3′-seq
   peaks, or NET-seq-refined CPAs already bundled in `rectify/data/*netseq*`) and a per-aligner
   *corrected*-3′-end concordance script. Report `|corrected_3prime − truth|` distribution **per
   aligner, independent of selection**. Without this, "improved accuracy" cannot be distinguished
   from "won the selection metric". **ESTABLISHED method.**

The penalty / STR tables are already bundled and `--Scer`-resolved (`penalty_scores.tsv`,
`str_penalty_scores.tsv` under `rectify/data/genomes/.../penalty_tables/`), so no regeneration is
needed before P1–P8; the only table work is validating/versioning them for the target chemistry.

This harness is effort-light and impact-multiplying; it is assumed done before P1–P8 are validated.

---

## The five proposals (then four more), each: mechanism · failure fixed · feasibility · impact · validation · risk · ESTABLISHED/NOVEL

### 1. Error-model-aware gap costs injected into the aligner DP (R10.4.1 AT/CG tables)

**Mechanism.** minimap2's KSW2 uses a *generic two-piece affine* gap cost `γ(l)=min(q+l·e, q̃+l·ẽ)`
with `q=2,e=1` (short) — it has **no idea** that a deletion inside an A-run is cheap/expected while
a deletion inside a CG run is rarer. The R10.4.1 reality (from RECTIFY's own tables) is the
opposite of generic: `D/AT` penalty 0.37→0.28 across HP=1→8 vs `D/CG` 0.58→0.37; insertions all
~10.0; substitutions rise with HP length. **Inject these as position-specific gap costs into the
DP** so the alignment *expects* HP deletions in A-runs and stops distributing them arbitrarily.

Two injection points, in increasing invasiveness:
- **(1a) Aligner-internal — minimap2 `--junc-bed` *score mode* is the only sanctioned hook, and it
  is for junctions, not HP gaps.** There is **no public minimap2 flag** to supply per-base gap
  costs; you would have to fork KSW2 (`ksw2_exts2_sse.c`) to make `q/e` a function of the
  reference HP-context at each column. High effort, fragile, and `mappy` isn't even importable here.
  **Not recommended as the first move.**
- **(1b) RECTIFY-side realignment DP (recommended).** RECTIFY *already* runs an HP-aware,
  table-driven semi-global DP (`_score_hp_dp_numba`, fed by `HpPenaltyTable`). Today it scores
  *junctions* (Module 2H). Generalize the **same kernel** to re-place the **3′-terminal window's
  indels** (P2 below) and the N-vs-D boundary — i.e. apply the empirical tables to the part of the
  alignment that decides the CPA, where it matters, instead of post-hoc heuristics (2C/2E).

**Failure fixed.** "HP indel noise distributed arbitrarily within the A-run; affine gaps don't pin
a unique boundary" (minimap2 §6; DEEP_DIVE §7) — the staircase artifact CLAUDE.md already fights
with `_CANONICAL_HP_PRIOR`. Making the *placement* (not just the tie-break) table-aware removes the
ambiguity at its source.

**Feasibility.** **(1b) HIGH** — kernel, loader, and Numba JIT exist; only the *call site* (extend
from junction window to 3′ window) is new. **(1a) LOW** — KSW2 fork + no mappy.

**Expected impact.** High on the exact metric (CPA), because it directly governs where the terminal
indel — and thus `reference_end` — lands. Medium-effort for (1b).

**Validation.** P0 concordance harness, ablation: tables-on vs tables-off on the 3′ window;
expect the corrected-3′ error distribution to tighten and cross-aligner agreement to grow *for the
right reason* (reads converge on the table-supported CPA).

**Risk.** Tables are **R10.4.1/yeast-specific** (COMPARISON §5; must not transfer to HiFi/other
chemistries — gate on chemistry). Over-cheap A-run deletions could *over*-walk-back the 3′ end —
bound by the same `max_boundary_shift`/`scan_limit` guards already in 2E/2H.

**ESTABLISHED** (empirical penalty matrices are standard; Guppy/medaka/`f5c` use per-context error
models). **NOVEL** *in this codebase only* as a DP-placement input rather than a post-hoc edit-
distance tie-breaker.

---

### 2. HP-aware targeted full-SW re-alignment of the 3′-terminal window

**Mechanism.** The single robust empirical signal in the whole investigation is that **mapPacBio's
full affine DP beats minimap2's heuristic chain-anchored extension on terminal-boundary precision**
(COMPARISON §3c; DEEP_DIVE §3) — *and* that this is exactly the property RECTIFY needs. But
mapPacBio is slow, JVM-heavy, 6-kb-capped, and splice-blind. **Capture its one advantage cheaply:**
take the chosen alignment (from any panel aligner), extract a small reference window around the
reported 3′ end (± a few hundred bp, bounded by the last junction), and **re-align just that window
with a full HP-aware DP** (the `_score_hp_dp_numba` / `local_aligner` Gotoh engine, or `edlib` for
a fast first pass), reporting the column-optimal terminal-exon boundary. This is the
gapmm2/edlib idea (HW infix realignment of an orphan end) but **HP-aware and applied to the CPA**,
not to recovering a missing 5′ exon.

**Failure fixed.** minimap2/deSALT/gapmm2 "stop at first good-enough chain → indel a few bases off
in a homopolymer" (DEEP_DIVE §3, §7). Gives every aligner mapPacBio-class terminal precision
*without* running mapPacBio.

**Feasibility.** **HIGH.** `edlib` is importable now; the affine + HP-DP kernels exist. The window
is O(hundreds of bp) so cost is negligible vs a whole-read DP — runs inside `rectify correct`'s
existing per-read loop. No new aligner dependency.

**Expected impact.** **Highest of all** — it directly manufactures the "column-optimal boundary"
that the synthesis identifies as the *only* mechanistic edge actually tied to CPA accuracy, and
applies it to all five aligners' outputs (so it also reduces dependence on mapPacBio being in the
panel — a hedge against its 6-kb/OOM fragility).

**Validation.** P0 harness; ablate window-realign on/off; specifically test against mapPacBio's own
3′ ends (does the cheap re-align match mapPacBio's column-optimal boundary on reads where mapPacBio
won?). Use the bundled Cat1/Cat2 (poly-A walk-back / soft-clip rescue) validation reads.

**Risk.** Window must be **bounded by the nearest upstream N-op** so the re-align can't slide the
3′ end across an intron (the exact bug `find_polya_boundary` v2.9.4 already guards). Poly-A bases in
the window must be masked (see P3) or the DP will happily extend into a genomic A-run.

**ESTABLISHED** (targeted local realignment = ABRA2/GATK indel-realignment, gapmm2's edlib pass).
**NOVEL** in making it HP-aware and CPA-targeted.

---

### 3. Poly-A-aware terminal DP that stops at the true CPA (the core failure mode)

**Mechanism.** Every aligner surveyed — minimap2 (`--end-bonus 0` → soft-clip), deSALT/GMAP (slip
into genomic A-run), uLTRA (strips poly-A entirely) — mishandles the poly-A→genome boundary because
**none of them know the A-run is a poly-A tail, not genomic sequence.** Propose a **two-state
terminal DP** for the 3′ window: state `G` (genomic alignment) and state `T` (poly-A tail, matches
A on `+` / T on `−` at *zero cost regardless of reference base*). The DP is forced to transition
`G→T` exactly once and never back; the **CPA is the transition column**. The tail no longer competes
with the genome for alignment score, so the boundary cannot slide into a downstream genomic A-run,
and there is no incentive to soft-clip. This is a principled replacement for minimap2's
`--end-bonus` hack and for RECTIFY's current heuristic chain (2B trim → 2E walk-back → 2G soft-clip
rescue), unifying them into one optimal boundary call.

**Failure fixed.** *The* core RECTIFY failure mode: 3′ end soft-clipped or slid into a genomic
A-run (COMPARISON §3'-end row; DEEP_DIVE §7, all aligners). Directly replaces the
`find_polya_boundary` false-stop guards (v3.0.3) and the 2E/2G priority dance with a single optimal
decision.

**Feasibility.** **MEDIUM-HIGH.** It is a ~10-line extension of the existing rolling-DP kernel: add
one column-state and a one-way `G→T` transition; the poly-A length is already in the parquet trim
metadata (DRS Step 0) to bound the tail. Pure RECTIFY-side; no aligner change. (Doing it *inside*
minimap2 would require a KSW2 fork + a tail-flag — not worth it given mappy isn't even available.)

**Expected impact.** **High** — this is the named root cause; a model-based CPA boundary should
dominate the stack of heuristics it replaces and remove their mutual-cancellation bugs (2G-vs-2E
priority, terminal-D stripping, etc. — a whole class of CLAUDE.md fixes).

**Validation.** Cat1 (atract_ambiguity/polya_walkback) and Cat2 (softclip_rescue) bundled reads
have exact expected shifts; assert the two-state DP reproduces them. P0 harness on real CPAs.
A/B against the current 2B/2E/2G chain.

**Risk.** dT-primed cDNA (`--dT-primed-cDNA`) has the A-context *off-read* (downstream of the 3′
end) — the `T`-state must be **disabled** for that protocol and AG-mispriming handled separately, as
today. DRS internal-priming false tails: bound `T`-state length by the parquet `pt` length.

**ESTABLISHED in spirit** (PacBio `isoseq refine --require-polya` is a poly-A-aware boundary;
nanopolish/`f5c` segment poly-A). **NOVEL** as a two-state terminal *alignment* DP that emits the
CPA as a transition column.

---

### 4. Chaining improvements: 3′-end-anchored chaining + sub-optimal terminal-chain retention

**Mechanism.** minimap2's SDP keeps the single best chain (`--secondary=no` in RECTIFY) and its
`h=50` predecessor cap + Z-drop can **drop a short, noisy 3′-terminal exon** whose seeds are sparse
— the terminal exon is precisely where seeds are weakest (poly-A, HP errors) and where RECTIFY's
value lives. Three concrete, separable changes:

- **(4a) Retain sub-optimal terminal chains.** Run minimap2 with `--secondary=yes -N` (or lower
  `--secondary` suppression) **only to harvest alternative 3′-terminal-exon placements**, then let
  the P2 window-realign + P3 poly-A DP pick the best-supported terminus. Fixes "`--secondary=no`
  discards alternatives that might have placed the 3′ end better" (minimap2 §6). **ESTABLISHED**
  (flag-only). **HIGH feasibility** — pure invocation change in `run_minimap2`.
- **(4b) 3′-end-anchored chaining bias.** When the read is stranded DRS (`-uf`), the 3′ end is the
  most informative anchor (it's the CPA). Add an **end-anchor bonus** at the 3′-most anchor so the
  chain is built to *reach* the read's 3′ end rather than truncate it — conceptually
  minimap2's `--end-bonus` but applied to chaining, not base DP, and *only* at the transcript 3′.
  minimap2 has `--end-bonus` (base-DP) but **no chaining-end-anchor flag** → would need a small
  patch or a RECTIFY-side re-chain of harvested anchors. **NOVEL. MEDIUM feasibility.**
- **(4c) RMQ vs range chaining for short reads.** For QuantSeq/short-read dT-primed cDNA (the
  `--short-read` path: bbmap+bwa), long-gap RMQ chaining (minimap2 2021) is irrelevant — reads
  rarely span introns — and range chaining with a tight band is faster and less prone to spurious
  long-join introns. Tradeoff: **use range/short-gap chaining for short reads, RMQ only for
  long-read intron-spanning DRS.** **ESTABLISHED tradeoff; LOW-effort** (already implicitly true via
  separate short-read aligners; document + verify band settings).

**Failure fixed.** Dropped/truncated terminal exons and intron-boundary jitter from committing to
one chain on noisy ends.

**Expected impact.** Medium. (4a) is cheap insurance; (4b) is the highest-value but highest-risk;
(4c) is hygiene.

**Validation.** Count reads whose 3′-terminal exon changes when secondaries are harvested; measure
net CPA-concordance gain (P0). Guard against (4a) inflating false termini via the existing
overhang/chimera filter.

**Risk.** (4a) more candidates → more chances for a wrong terminus to win selection (now via a
lower HP-edit-distance, since N-ops are free) — **must** be paired with P0 truth measurement and
the chimera filter. (4b) end-anchor bonus can pull the 3′ end into a genomic A-run — must compose
with P3's poly-A state, not fight it.

---

### 5. deSALT-style cross-read chain/junction consensus, made the explicit selection signal

**Mechanism.** The synthesis's strongest *structural* differentiator is deSALT's **cross-read exon
pooling**: it integrates all reads' skeletons, agrees on consensus exon boundaries, and snaps each
noisy read to them — producing *homogeneous* junctions (DEEP_DIVE §2-3; verified FACT-PAPER,
redteam_denovo D4). RECTIFY's live HP-edit-distance selector ranks each read's aligners on a
genome-anchored edit distance, so cross-read homogeneity enters only *implicitly* (and, in the
fallback sort, via the `_n_agree` popularity tiebreaker — herd bias, not accuracy). **Make it
explicit and principled:** after correction, build a per-locus **cross-read CPA/junction histogram**
(pool all reads' corrected 3′ ends and N-op boundaries), then **snap each read's corrected 3′ end to
the nearest well-supported peak** within a small radius — the same operation RECTIFY's `analyze`
clustering already does for count matrices, but fed *back* into per-read correction. This gives
deSALT's homogeneity to *every* aligner's reads, and supplies a real, defensible cross-read consensus
signal that the HP-ED selector's 3′-end term (and the fallback `_n_agree`) can consume.

**Failure fixed.** Per-read junction/CPA idiosyncrasy (minimap2/gapmm2/mapPacBio "commit per-read in
isolation", COMPARISON §3d); and the absence of an explicit, de-herded cross-read consensus the
selector can lean on.

**Feasibility.** **MEDIUM.** The clustering primitives exist (`analyze` IntervalTree). The work is
plumbing a two-pass correction (pass 1 = correct + collect peaks; pass 2 = snap), which fits the
existing streaming/chunked architecture but needs a cross-chunk merge of peaks (non-trivial under
SLURM arrays — peaks must be globally pooled before pass 2, like deSALT's `-B 655350` batch).

**Expected impact.** Medium-high on homogeneity and on *defensibility* of the win-rate metric; it
turns the most-cited (and most-criticized) selection term into a sound one.

**Validation.** Measure junction-placement variance across reads at a locus before/after; P0 CPA
concordance. Critically: **leave-one-aligner-out** (redteam (d)#4) to confirm the consensus is
read-driven, not deSALT-driven herd.

**Risk.** Cross-read snapping can **erase genuine biological CPA heterogeneity** (alternative
polyadenylation is real and is often the *point* of a 3′ study). Must snap only within tight HP-noise
radius (≤ the empirical HP-deletion length), never across distinct APA peaks. This is the dominant
risk of the whole proposal set — over-homogenization destroys the signal.

**ESTABLISHED** (deSALT proves the mechanism; APA peak-calling is standard). **NOVEL** as a feedback
loop into per-read correction + as a principled replacement for `_n_agree`.

---

## Additional proposals (briefer)

### 6. Fix deSALT's missing `-x ont2d` preset (the cheapest experiment in the whole investigation)
deSALT runs `null` (~13%) on ONT DRS (COMPARISON §6.1; verified). **Mechanism/feasibility:** add
`-x ont2d` (or `ont1d`) in `run_desalt` — one line. **Impact:** unknown but possibly large on the
78.9%-winning aligner; it is currently the single most concrete "is this intentional?" item.
**Validation:** A/B `null` vs `ont2d` on P0 harness. **Risk:** could *worsen* sensitivity (null
evidently works); purely empirical. **ESTABLISHED** (vendor preset). **Effort: trivial.**

### 7. Audit mapPacBio's `maxindel` sizing for yeast introns
mapPacBio already sets `intronlen=10` and an explicit `maxindel=max(200000, max_intron)`
(multi_aligner.py:749,754), so it does **not** rely on BBMap's soft 16000 default — the search
bound is derived from the same `max_intron` that feeds minimap2 `-G` / gapmm2 `-i`. The residual
work is a **soundness audit, not a fix:** confirm `max_intron` is set appropriately for the target
organism (yeast introns <1 kb sit far inside 200k) and that an oversized `maxindel` does not inflate
runaway-gap search cost. **Mechanism:** verify/tune the `max_intron` input. **Feasibility:** config.
**Impact:** low/safety. **Risk:** too-small `max_intron` → missed long introns. **ESTABLISHED.**

### 8. `--end-bonus` sweep for minimap2 (controlled, not the splice-flank red herring)
`--end-bonus 0` is the *documented* reason minimap2 soft-clips poly-A (minimap2 §6, M9 — properly
hedged). **Do not** lean on `--splice-flank` (caveated as "compatibility", not 3′-accuracy).
Instead **sweep `--end-bonus`** (e.g. 0/5/10) and measure CPA concordance. **Mechanism:** raises the
reward for extending to the read 3′ end. **Risk:** a blunt end-bonus extends into genomic A-runs
(the opposite failure) — which is *why* P3 (poly-A-aware DP) is the principled fix and this is only a
diagnostic. **ESTABLISHED (flag); diagnostic value only.**

### 9. Replace heuristic Z-drop terminal break with explicit terminal-exon model
minimap2's Z-drop decides where the terminal exon ends by score-collapse, not by sequence model
(minimap2 §2.3). Subsumed by P2+P3 (window realign + poly-A DP make the terminal-exon end an
*explicit* optimum). Listed for completeness; **no separate work** beyond P2/P3.

---

## In-aligner DP vs RECTIFY post-correction — the tradeoff (explicit argument)

**Argue for post-correction (RECTIFY-side) as the default, with two narrow exceptions.**

*Why post-correction wins for almost everything here:*
- **The aligners can't be cheaply changed.** `mappy`/`parasail` aren't even importable in this env;
  the only sanctioned minimap2 hooks are `--junc-bed` (junctions) and `--end-bonus` (blunt). Real
  gap-cost/poly-A injection means **forking KSW2** — a maintenance and reproducibility liability
  across five aligners, four of which RECTIFY doesn't control.
- **RECTIFY already owns the right machinery** at the right altitude: `HpPenaltyTable`, a Numba
  HP-aware DP, the affine `local_aligner`, Module 2H, and the DRS poly-A metadata. P1–P3 are
  *extensions of existing kernels*, not new engines.
- **Correct-FIRST + ensemble is the whole architecture.** Fixing one aligner's DP helps that
  aligner's reads; a RECTIFY-side terminal DP helps **all five aligners' reads uniformly** and
  composes with cross-read consensus (P5). It also **decouples** the fix from any single aligner's
  fragility (mapPacBio's 6-kb crash, deSALT's dedup/SIGSEGV, gapmm2's PAF reconstruction).
- **Speed.** GMAP proves in-aligner per-read splice DP is ~30× slower and disqualifying
  (DEEP_DIVE §5b). A *windowed* RECTIFY-side DP (P2/P3, hundreds of bp) is negligible.

*The two cases where in-aligner is genuinely better:*
1. **Seeding/chaining of short terminal exons (P4).** No amount of post-correction recovers a
   terminal exon the aligner never seeded/chained. Chaining-side fixes (4a harvest secondaries, 4b
   end-anchor) must happen *in or around the aligner* — post-correction can only refine what exists.
   This is the strongest argument for touching the aligner at all.
2. **Cross-read consensus (P5).** deSALT shows this is most powerful *during* alignment (it re-aligns
   against pooled exons). RECTIFY's post-hoc version (snap-to-peak) is a strict approximation; it's
   the right tradeoff only because forking five aligners to share a skeleton pool is infeasible.

**Net:** put the *boundary precision* fixes (P1–P3) in RECTIFY's correction DP; put the *candidate
generation* fixes (P4) at the aligner invocation; make consensus (P5) explicit post-hoc. Never fork
KSW2 unless a windowed RECTIFY DP provably can't reach a measured target.

---

## Ranking by impact / effort

| Rank | Proposal | Impact | Effort | Risk | Tag |
|---|---|---|---|---|---|
| **1** | **P0** ground-truth harness (penalty tables already bundled) | (gates all) | Low | Low | ESTABLISHED |
| **2** | **P3** poly-A-aware two-state terminal DP (CPA = transition) | **Highest** (core failure) | Med | Med (protocol gating) | NOVEL kernel, established spirit |
| **3** | **P2** HP-aware full-SW realign of 3′ window | **Highest** (mapPacBio edge, cheap, all aligners) | Low-Med | Low-Med | ESTABLISHED+NOVEL |
| **4** | **P6** deSALT `-x ont2d` A/B | Unknown, possibly large on top aligner | **Trivial** | Low | ESTABLISHED |
| **5** | **P1b** empirical gap costs into the realign DP | High (CPA placement) | Med | Med (chemistry-specific) | ESTABLISHED+NOVEL |
| **6** | **P4a** harvest sub-optimal terminal chains | Medium | Low | Med (false termini) | ESTABLISHED |
| **7** | **P5** explicit cross-read CPA consensus (feed the selector's 3′ term) | Med-High (defensibility) | Med-High | **High (kills APA signal)** | ESTABLISHED+NOVEL |
| **8** | **P8** minimap2 `--end-bonus` sweep (diagnostic only) | Low (diagnostic) | Trivial | Med if shipped | ESTABLISHED |
| **9** | **P7** mapPacBio `maxindel` safety | Low/safety | Trivial | Low | ESTABLISHED |
| **10** | **P4b** 3′-end-anchored chaining bonus | Med (if it works) | Med-High | Med-High | NOVEL |
| —  | **P1a** KSW2 fork for in-DP gap costs | High in theory | **Very High** | High | NOVEL — **defer** |

**Recommended sprint:** P0 → P2 → P3 (the three together turn RECTIFY's stack of mutually-cancelling
3′ heuristics into one measured, poly-A-aware, HP-aware optimal terminal DP applied to every
aligner's output) → P6 (one-line, high-leverage on the top aligner) → P1b → P4a. Defer P1a (KSW2
fork) and approach P5 only with hard APA-preservation guards.
