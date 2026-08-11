# DIRECTOR ALGO EVAL — novel-isoform discovery on the RECTIFY long-read panel (Fable)

*Independent DIRECTOR-level algorithmic evaluation. Mandate: (1) tear down the
5-aligner long-read panel for NOVEL-ISOFORM discovery; (2) rank orthogonal
algorithmic approaches, each gate-testable against constructed truth. Neutral +
adversarial. Report back; do NOT commit. Render date 2026-07-01.*

Sources read: `RECTIFY_STRATEGIC_FRAME.md`, `NATIVE_ALIGNER_OVERVIEW.md`,
`SIMULATION_BENCHMARK_SPEC.md` + `scorer.py`/`novel_junction_blindspot.py` (gate),
`HANDOFF_ALIGNER_BENCHMARK.md` WS-1 + `error_injector.py` (empirical error),
`docs/aligners/{minimap2,gapmm2,deSALT}.md`, `dev/DELIVERABLE_B_FINDINGS.md` +
`dev/ALIGNER_INVESTIGATION_SYNTHESIS.md` (the gmap-111 work), `docs/ALIGNER_RECOMMENDATIONS.md`
(minisplice GT-AG bias), `multi_aligner.py` (actual invocations), `penalty_tables_quickref.md`.

*Provenance note (Fable OUTWARD pass, 2026-07-01): independently re-verified the load-bearing
claims. CONFIRMED: minisplice is GT-AG-canonical-only (`docs/ALIGNER_RECOMMENDATIONS.md:145`, 96%
GT-AG); cDNA error tables + UMI dedup-before-align exist (`penalty_scores_cdna.tsv`,
`cdna_correct_command.py --per-cluster-cap`); gapmm2 = minimap2+PAF-cs-tag wrapper
(`docs/aligners/gapmm2.md`); deSALT SIGSEGV/`-G` (`docs/aligners/deSALT.md:75`); the gmap-111 →
~0/111 validate as artifacts (`dev/DELIVERABLE_B_FINDINGS.md:60`). The gapmm2 "≥85% forced-canonical
/ 1.3–1.8% unique" and minimap2 "~49% win-share" figures below are QUALITATIVE (directionally
supported by the panel-redundancy work; exact percentages not re-traced to a committed source).*

---

## 0. The two load-bearing measurements (everything hangs on these)

**A. Isoform-flattening is real and now SIZED** (`novel_junction_blindspot.py`,
minimap2 `-ax splice -uf`):

| donor–acceptor | deviation from GT-AG | recovery | **blindspot (snap)** |
|---|---|---|---|
| GT-AG (control) | 0 | 0.983 | 0.017 |
| AT-AC | 2 | 0.533 | **0.467** |
| GA-AG (1-off) | 1 | 0.300 | **0.700** |
| CT-AC (2-off) | 2 | 0.217 | **0.783** |
| CA-TC (deep) | 4 | 0.100 | **0.900** |

Blindspot rises monotonically with motif deviation — on **error-free** reads. This
is the headline weakness, empirically quantified, not asserted.

**B. The empirical error model says WHERE independence can come from** (WS-1):
RNA004 = clean bulk ~1% + a ~1–2%-of-reads **hot tail at ~12%** (13× bulk, all
MAPQ 60, genuine chemistry). Deletions ~2× insertions (sub 0.55 / del 0.30 / ins
0.15); ≥2 bp indels = 0.39. **Flat-Q null (commit d412750): per-*base* quality adds
nothing over the error-type table.** Over-dispersion is LOCAL (autocorr r≈0.30) and
is a per-*read* property. → The orthogonality lever is **per-read reliability +
calibrated −logP emission**, NOT per-base Q.

---

## 1. Panel teardown (per-aligner, for discovery)

All five share one engine DNA: **seed→chain→flat-affine DP with a hard/である near-hard
GT-AG canonical bonus, quality-blind.** They co-fail on exactly the axis in §0A.

**minimap2** (`-ax splice -uf -k14 -G5000 --splice-flank=no --junc-bonus 9`).
Seed-chain + 2-piece affine, `splice` model rewards GT-AG and (with `--junc-bed`)
annotated sites. **Discovery failure:** snaps non-canonical/unannotated donors to the
nearest canonical motif (§0A). `-G5000` yeast default silently clips >5 kb human
introns (alt-TSS/distal-exon loss) unless `--max-intron` raised. No base-Q use.
The consensus workhorse (the plurality of per-read wins) → its bias dominates the panel.

**gapmm2** (minimap2 core + edlib terminal-exon rescue, PAF→BAM). Same optimum as
minimap2 **by construction**, plus a **forced-GT-AG terminal-rescue artifact**:
the large majority of its unique introns are clustered forced-canonical terminal rescues; only
a small fraction of introns are gapmm2-unique and it is largely redundant with uLTRA/deSALT (exact
percentages qualitative — see provenance note). Already
**dropped from the human-DRS panel**. For discovery it is *anti-orthogonal* — it
manufactures spurious canonical ends. Near-zero novelty value.

**uLTRA** (collinear chaining over an **annotation-derived** exon graph; needs GFF).
Genuinely strong on tiny exons (11–20 nt) seed-chainers miss — a real orthogonal
axis for **short internal/5′-terminal exons**. But it is **annotation-snapping**: a
novel exon/junction absent from the GFF is off-graph and cannot be reported. Discovers
alt-splicing *within known exon inventories*, not truly novel sites.

**deSALT** (De Bruijn RdBG index, two-pass, high-sensitivity spliced). Best raw
sensitivity for spliced placement and small-exon spanning; annotation `-G` disabled
(SIGSEGV) so it is less annotation-bound than uLTRA. Still a flat-affine optimum with
a canonical splice preference → snaps non-canonical sites like minimap2; emits each
aln N× (deduped). Orthogonal in *seeding* (kmer-DBG vs minimizer) but NOT in the
scoring axis that flattens novelty.

**gmap** (`-f samse --npaths 1`, canonical-preferring sandwich DP, pinned 2021-05-27).
**High-sensitivity/high-noise:** ~97% of its unique junctions are non-canonical and
mostly spurious; kept only because the integer `score_segment` **fence** suppresses
that mass. Validated unique real-novel yield is ~0–1 canonical junction/chromosome
(3/111 on chr5). It is the panel's only member that *emits* deep non-canonical
placements — but without calibrated scoring that signal is 96% artifact. Value =
corroborator, not discoverer.

**The herd:** minimap2 ≈ gapmm2 (same core) ≫ shared canonical-affine optimum with
deSALT/gmap; uLTRA/deSALT add *seeding/short-exon* orthogonality but re-converge at
the scoring step. **All five agree for the wrong reason on non-canonical novelty and
none consumes the empirical error model.** The panel-failure surface is precisely §0A
× the hot-read tail: reads whose true site is non-canonical AND whose junction-region
evidence is bursty.

**Genuine orthogonality that already exists:** uLTRA's short-exon graph chaining;
deSALT's DBG seeding; gmap's willingness to emit non-canonical (wasted without
calibration). Everything else is correlated. (**mapPacBio** — a 6th live tier-1 member
the task's "5" omits — is another correlated BBTools affine placer, badly intron-
misconfigured on human; it adds no orthogonal axis.)

### 1b. Per-aligner × per-target failure grid (● herd-snap · ◐ partial-orthogonal · ○ n/a→correct-step)

| target | minimap2 | gapmm2 | uLTRA | deSALT | gmap |
|---|---|---|---|---|---|
| alt-TSS / 5′ end | ● clips distal exon at `-G`; snaps | ● + forced canonical terminal rescue | ◐ short 5′-exon recovery if annotated | ◐ good sensitivity, still affine | ● canonical-preferring |
| unannotated intron | ● snaps to nearest canonical | ● anti-orthogonal (fabricates) | ● off-graph → cannot report | ◐ emits, but canonical-biased | ◐ emits (96% artifact w/o calib) |
| non-canonical intron | ● 0.47–0.90 blindspot (§0A) | ● | ● motif+annotation snap | ● snaps | ◐ emits deep non-canon, uncalibrated |
| cryptic pA extend/truncate | ○ CPA not placed by ANY panel member — walkback correct-step | ○ | ○ | ○ | ○ |
| variant-adjacent splice | ● del near SS → spurious intron | ● | ● | ● | ● none is variant-aware |
| paralog / allele (SMN) | ● lone distinguishing base ambiguous | ● | ● | ● | ● no per-molecule pooling |

The grid makes the herd visible: five near-identical ● columns on the splice targets,
a unanimous ○ on CPA (no member even attempts it — it is purely a correct-step), and a
unanimous ● on variant/paralog (no member models either). The only ◐ orthogonality is
seeding/short-exon (uLTRA/deSALT) and raw non-canonical emission (gmap) — and all of it
re-converges at the flat-affine, motif-biased *scoring* step.

## 1.5 STEELMAN — "the panel + shipped correct-steps are already near-optimal"

The strongest case AGAINST building anything, stated fairly: (a) **C3 was refuted
because the arbiter is AT CEILING on recoverable reads** — it already picks the
truth-correct member even at 100% member disagreement and with minimap2 snapping 47%,
so no calibrated arbitration adds accuracy. (b) **Validated unique novel yield is tiny**
— gmap's best 111 chr5 candidates reduce to ~3-4 real; genome-wide extrapolation is
tens, not thousands. (c) **Much of the deep-deviation blindspot is ZERO-EVIDENCE** — at
dev≥2 with a canonical motif nearby, the spurious placement is a 0-mismatch (NM==0) tie
with truth, so NO function of the read can recover it; the honest ceiling, not a
fixable gap. (d) **Flat-Q is null** — the obvious "add quality" lever provably buys
nothing. (e) C1 (the one confirmed facet) **already flows through the existing arbiter**
with no new machinery.

**Why it does not fully hold — the one hinge:** every one of those refutations tested
whether the arbiter *picks* the truth **from a menu the panel produced**. None tested
whether the panel **produces the non-canonical truth at all**. §0A measures exactly
that gap — at dev≥2 the true site is ABSENT from every member's output on 47–90% of
*error-free* reads. That is not an arbitration problem (C3) and not zero-evidence (these
rungs are constructed addressable: truth is 0-mismatch, snap is NM≥1) — it is a
*generation* problem no panel member solves. So the near-optimal case is correct about
arbitration and about the zero-evidence floor, and wrong about generation. The only
proposals worth building are the ones that put addressable non-canonical truth ON the
menu — which is precisely the gate-move Rank 1 requires.

---

## 2. Ranked candidate approaches

Bar: independent error axis vs the flat-affine/quality-blind/snapping family +
DISCOVERS what it flattens + dependency-light + **gate-testable against constructed
truth** (metric = ambiguity-aware position concordance; guards = ADDRESSABILITY
[truth strictly wins on a motif-blind −logP, i.e. snap carries NM≥1] vs ZERO-EVIDENCE
[NM==0 tie, honest ceiling] + FDR no-over-call on clean/annotated).

*(literature grounding appended below in §2-lit)*

### Rank 1 — Calibrated-likelihood local REALIGNMENT in the panel window (the native member, C-facet spine)
**Mechanism:** take the panel's placement cluster as a localization window, re-score
every candidate junction/end inside it on the empirical **−logP** scale (HP/STR
deletion-law tables already shipped) with **NO motif prior** (or a soft, separable
one), emit a posterior + runner-up. **Orthogonality:** its error axis IS the empirical
error model, independent of flat-affine; it converts gmap/deSALT's raw non-canonical
emissions into *calibrated* calls and de-snaps minimap2 within the window.
**Exploits our data:** directly consumes the −logP tables + per-read hotness as a soft
down-weight; the hot-tail becomes an FDR guard (novel-call support must not be
tail-enriched). **Dependency:** none new (pure DP + existing tables). **Gate:** the
blindspot ladder §0A — recovery must rise on *addressable* rungs (snap NM≥1) with FDR
flat on canonical/annotated controls; zero-evidence rungs must NOT rise (honesty).
*This is the single most promising thing and it is already the C1/C3 spine — but the
gate must move from "arbiter picks truth from the menu" to "member PUTS the
non-canonical truth ON the menu" (the ladder), which C3's refutation did not test.*

### Rank 2 — Learned SEQUENCE-GRAMMAR splice prior as a separable emission term (SpliceAI/Pangolin/**Splam**, NOT minisplice)
**Mechanism:** a learned donor/acceptor scorer supplies a per-position splice log-odds
that REPLACES the hard GT-AG bonus with a graded, motif-agnostic prior inside Rank 1's
DP. **Splam (Chao 2024) is the closest published precedent — a DL model built expressly
to RESCORE spliced alignments**; SpliceAI/Pangolin predict non-canonical sites de novo.
**Critical correction:** the already-wrapped **minisplice models GT-AG ONLY**, so as the
prior it would *re-snap* against the exact non-canonical target of §0A — it is the wrong
learned model here; prefer a non-canonical-capable grammar model. **Orthogonality:**
learned sequence-grammar axis, independent of the affine score; but only if kept
*separable* so ADDRESSABILITY is scored motif-blind (a hard prior just re-herds toward
canonical). **Exploits our data:** calibrate on yeast R64 truth junctions; the prior
must be *soft* so non-canonical truth still wins on read evidence + −logP. **Dependency:**
medium (a small CNN; Splam/Pangolin are pip-installable, no graph engine). **Gate:**
ablate prior-on vs prior-off on the ladder; win ONLY if recovery rises WITHOUT raising
non-canonical FDR (the re-snapping failure mode is the thing to catch).

### Rank 3 — UMI / duplicate-molecule PARTIAL-ORDER consensus before placement (abPOA/SPOA)
**Mechanism:** for cDNA-UMI (and DRS PCR duplicates), pool reads sharing a UMI, build
an abPOA consensus per molecule, align the consensus once. **Orthogonality:** attacks
error at the SOURCE — averages away the per-read hot tail so the junction/end evidence
is clean, independent of any placer. **Exploits our data:** the hot-tail structure is
exactly what per-molecule consensus collapses (1–2% hot reads drowned by clean copies),
and **cDNA error characterization ALREADY exists** — `penalty_scores_cdna.tsv` /
`error_rates_cdna.tsv` + the `run_profiler_cdna_*` calibration scripts, and the cDNA
pipeline already does UMI dedup-before-align + cross-orient merge (`--per-cluster-cap`,
`cdna_umi_architecture`). So this is "**exploit the existing cDNA-UMI characterization**,"
NOT an unmodeled gap (an earlier draft overstated this). The genuine gaps are two: (i)
that dedup is *collapse*, not a calibrated POA *consensus* that beats down the hot tail
before placement; (ii) **the benchmark/gate has NO molecule-level UMI stratum** (per the
sim digest) — so this needs a NEW molecule-truth harness, not the existing gate.
**Dependency:** abPOA is light (C lib / pip). **Gate (new harness required):** simulate
N reads/molecule at the empirical cDNA hot-tail mixture, truth = molecule sequence;
measure junction/CPA recovery of calibrated-consensus-then-align vs the current
dedup+align across coverage. **This is a CORRECT-STEP / pre-aligner improvement, not a
native-aligner facet** — and the highest-leverage cDNA-specific win.

### Rank 4 — Variant-aware emission (C6) for cis-variant & SMN1/SMN2 paralog isoforms
**Mechanism:** given a matched-DNA VCF/haplotype, bake known SNVs/indels into the
reference emission so a splice-site-adjacent variant is scored as *expected*, not as a
spurious intron; for paralogs, POA-pool per-paralog then project. **Orthogonality:**
high — no panel member is variant-aware; a deletion near a splice site is re-expressed
as a fake intron by all five. **Exploits our data:** the −logP emission already
supports per-position cost; add a variant-conditioned term. **Dependency:** light
(VCF, no graph engine — reject GraphAligner/vg as dependency-violating heavyweights
that also import the affine optimum). **Gate:** C6 stratum is BUILT + discriminating
(`gen_variant_stratum`, `fp_variant_adjacent`); arms = variant-BLIND (fabricates
intron) vs variant-AWARE (recovers), + specificity control (must not suppress a REAL
near-variant novel junction). **SMN1/SMN2 specifically needs C4 POA-pooling** (per-
paralog molecule consensus at the lone distinguishing base), not a graph aligner —
it is a native-member feature.

### Rank 5 — FracMinHash containment fallback for the panel-failure TAIL (C5)
**Mechanism:** for reads with NO acceptable panel window (~12% unmapped), a
k-mer-containment localizer (sourmash-style) proposes a locus, then Rank 1 realigns.
**Orthogonality:** total — it is the ONLY mechanism for reads the whole panel misses;
alignment-free, independent of affine chaining. **Exploits our data:** none directly
(localization only). **Dependency:** medium (FracMinHash sketch). **Gate:** GATED on a
MEASURED depletion trigger — Tier-2 must first SIZE the no-window tail (spec flags this
as unmeasured). Do not build until the tail is sized.

### Rejected (do not re-propose)
- **A 6th flat-affine/WFA placer** — shares minimap2's optimum; adds correlated votes,
  not information (WFA is fine as *infrastructure* inside Rank 1, not a member).
- **GraphAligner / vg / pangenome / r-index** — dependency-violating heavyweights that
  ALSO inherit an affine optimum; the variant payoff is reachable via Rank 4's emission
  term at a fraction of the cost.
- **StringTie2-lr / Bambu / FLAIR / IsoQuant / ESPRESSO as the discovery engine** —
  these consume a minimap2 BAM (see §2-lit), so they INHERIT the §0A snapping bias;
  useful as downstream isoform *assemblers*, useless as the orthogonal placer. Bambu's
  learned novel-junction classifier is a good FDR-model reference for Rank 1's guard,
  not an aligner.

---

## 3. Cis-variant / SMN1-SMN2 / trans-splicing specifics
- **cis-variant intron/exon isoforms:** Rank 4 emission term (CORRECT-step + member).
  A splice-adjacent SNV that creates/destroys a site is the exact case all five snap
  wrong; variant-aware emission is the fix, gate-tested on the built C6 stratum.
- **SMN1/SMN2 exon-7 skipping is BOTH a paralog AND a cis-variant problem:** (i) the
  paralog axis — SMN1 vs SMN2 differ by a lone distinguishing base, so per-read placement
  is ambiguous → needs **C4 POA-pooled per-paralog consensus** (cluster → consensus →
  align once → project), NOT a graph aligner; (ii) the cis-variant axis — the SMN2
  c.840C>T substitution *weakens the exon-7 splice enhancer* and drives the skip, so
  recovering *why* the isoform switches is a **Rank 4 variant-aware emission** case. Both
  are native-member/correct-step features; graph aligners are not needed for either.
  Allele/paralog-aware POA pooling + UMI per-molecule consensus (Rank 3) compound here.
- **Chimeric / trans-splicing:** NO stratum exists (only `chimeric_consensus.py` in the
  aligner). This is a genuine gate GAP — add a chimera/fusion truth stratum before
  claiming trans-splicing discovery. Flag, don't build yet.

## 4. Top-3 BUILD-AND-GATE-FIRST
1. **Rank 1 native realigner, gated on the §0A blindspot ladder** (member puts
   non-canonical truth on the menu). Highest leverage, dependency-free, gate exists.
2. **Rank 3 UMI/molecule POA consensus** — exploits the EXISTING cDNA-UMI
   characterization; a correct-step win. Needs a NEW molecule-level truth harness (the
   current gate has no UMI stratum) — build that harness first, then gate.
3. **Rank 4 variant-aware emission** — the only route to cis-variant + SMN isoforms;
   C6 stratum already built and discriminating.

---

## §2-lit — literature grounding (2023–2026)

**The decisive finding:** every surveyed isoform tool CONSUMES a minimap2 `-ax splice`
BAM, so all INHERIT the §0A snapping bias — none is an independent base-level placer,
and none uses ONT basecall quality. RECTIFY's multi-aligner + empirical-−logP consensus
is not duplicated by any of them.

- **Correlated placers (share minimap2 affine/snapping) — reject as the discovery
  engine:** IsoQuant (Prjibelski, Nat Biotechnol 2023; intron-graph on minimap2 junctions),
  Bambu (Chen, Nat Methods 2023; XGBoost novel-discovery-rate on minimap2 BAM — and
  *explicitly cannot separate TSS/TES isoforms*, a direct hit to 3'-end discovery),
  FLAIR2/FLAMES (Tang, Genome Biol 2024; *actively snap unsupported novel junctions to
  the nearest annotated site* — worst case for novel fidelity), StringTie2 `-L` (least
  annotation-biased consumer but ends/junctions are still exactly minimap2's),
  Winnowmap2 (weighted minimizers, no spliced mode, no Q model), WFA2/BiWFA (same affine
  objective — infrastructure only).
- **Partial independent layer (useful as reference, not engine):** ESPRESSO (Gao, Sci Adv
  2023; realigns to high-confidence SJs + empirical per-read error profile — can override
  *some* minimap2 junctions but requires canonical motif for novel SJ) — a good model for
  Rank 1's realign-in-window; Sicelore (Lebrigand, Nat Commun 2020; POA+racon per-UMI
  consensus) — precedent for Rank 3.
- **Genuinely orthogonal axes:** **Splam** (Chao, Genome Biol 2024) + SpliceAI (2019) /
  Pangolin (Zeng 2022) — learned sequence grammar, predict non-canonical de novo → Rank 2;
  abPOA (Gao 2021) / SPOA — per-molecule redundancy → Rank 3; GraphAligner / vg Giraffe —
  graph topology BUT *no de novo junction discovery* (introns must be pre-encoded as
  edges) → confirms rejection for discovery, and its variant payoff is cheaper via Rank 4.
- **Paralog/SMN + UMI-POA (Sumner-relevant):** paralog resolution is PSV-based bespoke —
  **Paraphase**-style masked-reference remap + phasing, PacBio-HiFi SMN profiling (Chen,
  AJHG 2023), a dual-mode targeted nanopore SMN1/SMN2 assay (J Mol Diagn 2025); transcript-
  level SMN paralog disambiguation remains bespoke → supports the C4-POA-pooling route over
  a graph aligner. UMI-POA beyond Sicelore: wf-single-cell, scNanoGPS, Longcell, Scywalker
  (2025 benchmark) → Rank 3 references.
