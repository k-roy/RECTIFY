# Splice-Junction Placement in RECTIFY — Validation & Improvement

> Verified against `origin/drs-validation-rebuild` @ 366c885 (2026-06-19).

**Headline deliverable.** Primary goal: **validating and improving splice-junction
placement** in RECTIFY's correct-first long-read RNA-seq pipeline (ONT DRS/cDNA,
*S. cerevisiae*). 5'/3'-end (CPA/poly-A) resolution is treated as **already solved**
by the project owner and is **demoted** throughout — it appears only where it directly
interacts with junction placement.

**Convention.** **FACT** = verified against source/paper this pass. **INFERENCE** =
reasoned synthesis. **HYPOTHESIS** = plausible but untested (no committed measurement).
Code claims were re-verified against RECTIFY source (`junction_refiner.py`,
`junction_scoring.py`, `hp_penalty.py`, `junction_validator.py`,
`false_junction_filter.py`, `calibrate_junction_overhang.py`, `utils/junction_bed.py`,
`corrected_consensus.py`, and the two `merge_corrected_tsvs` call sites).

**Sources.** `01_investigation/{minimap2,deSALT,uLTRA,gmap,mapPacBio_bbmap,gapmm2}.md`;
`02_synthesis/{COMPARISON,DEEP_DIVE}.md`; `03_adversarial/{redteam_winrates_selection,
redteam_denovo}.md`; `04_discovery/{improve_splice_junctions,REDTEAM_proposals}.md`;
RECTIFY source + CLAUDE.md (Module 2H v3.1.0–v3.3.0, empirical penalty tables,
minimap2 `--junc-bed`/`--splice-flank`).

---

## 1. How each aligner places splice junctions

This is a **junction-placement** comparison, *not* a 3'-end comparison. The central
distinction is **what primary criterion decides where the `N` (intron-skip) op lands**:
gap-cost, cross-read consensus, annotation, sequence-splice-model, or "no model at all."

### minimap2 — GT-AG DP with soft annotation prior (gap-cost-first)
**FACT.** Splice DP is `ksw_exts2_sse` with a dedicated intron state (`x2`, open
`q2=32`, extend `e2=0` → intron extension is free once opened). Donor/acceptor are
scored by position-specific tables: donor `0` for `GTA/GTG`, `p/2` for `GTC/GTT`, `p`
otherwise; acceptor symmetric (`CAG/TAG`=0); `p = -C = 9` in the `splice` preset.
RECTIFY runs `-ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD
--junc-bed annotation.junc.bed --junc-bonus 9` (`run_minimap2`).
- **`--junc-bed`/`--junc-bonus 9`** (FACT): annotated junctions get an additive `+9`
  bonus in the DP — a **soft tie-breaker bias toward annotated boundaries, never a hard
  constraint**, so novel junctions remain discoverable. The BED is built from the GFF by
  `utils/junction_bed.py::generate_junction_bed` (cached as `annotation.junc.bed`).
- **`--splice-flank=no`** (FACT/CONTRADICTION): disables the `GTA/GTG..C/T AG` flanking
  bonus that nudges boundary placement ±1 bp toward the statistically-favored flank.
  **The live code comment says `# Disable for compatibility`** (`multi_aligner.py:252`),
  whereas CLAUDE.md and the minimap2 dossier frame it as "important for 3' end accuracy."
  Per `redteam_denovo.md` M7 / `REDTEAM_proposals.md` C7 this is a genuine internal
  contradiction: the elaborate "hand junction placement to RECTIFY" rationale is
  **HYPOTHESIS**, not the author's stated reason. Either way the *effect* is real:
  minimap2 scores only the raw `GT..AG` signal, leaving fine boundary placement to RECTIFY.
- **Failure mode** (FACT/INFERENCE): intron-boundary **jitter ±1–few bp** when
  exon-end+intron-start resembles intron-end+exon2-start (common at HP-rich junctions);
  low `A=1/B=2` scoring does not pin a unique column. Wins **0.1%** (but see §2 caveat).

### deSALT — cross-read junction consensus (the accuracy lead's mechanism)
**FACT-PAPER.** Genome Biology 2019: deSALT's Pass-2 exon inference *"integrates all the
alignment skeletons of all the reads"* (independently verified by `redteam_denovo.md` D4 —
this is NOT an overstatement). Mechanism: Pass-1 builds coarse skeletons (Match Blocks +
SDP) for the whole batch (`-B 655350`); their pooled projections define **consensus exon
boundaries** that are applied back to each individual noisy read in Pass-2. GT-AG enters as
the `-R 9` non-canonical penalty (canonical = free). RECTIFY runs `deSALT aln` with **no
`-x` preset** (`null` ~13% model) and **no `-G`** (yeast GTF → SIGSEGV).
- **Why this is the junction-placement lead** (INFERENCE): cross-read pooling produces
  **homogeneous** junctions — reads from one isoform get the *same* CIGAR junction
  placement rather than landing a few bp apart. This is a mini cross-read junction
  correction performed *internally, before* RECTIFY's own correction runs. **This is the
  single mechanism RECTIFY does not itself implement across aligners** (see §5/P-CRC).
- **Failure mode:** can place long-range / chimeric false junctions; specifically guarded
  against by RECTIFY's overhang filter (§2). Wins **78.9%** (caveat §2).

### uLTRA — annotation exon-segment colinear chaining (exact annotated junctions)
**FACT.** uLTRA tiles the GTF into "parts"→"segments"; a gene becomes a linear tiling and
an isoform a subset path. Pass-2 MAM chaining **injects all annotated exons of the
candidate gene — including seedless ones** — letting the chain tile a micro-exon no seed
hit. `find_exons()` classifies each MAM-gap: `(m1.y,m2.x) in valid_introns_sites` →
annotated `N` **snapped to the exact GTF donor/acceptor** (zero wobble); `>10 bp` → novel
`N`; else `D`. RECTIFY runs `uLTRA pipeline --ont --disable_infer`.
- **Junction behavior** (FACT): exact canonical junctions on FSM/known reads (its
  strength). **Structurally circular on novel** (INFERENCE): reads >10% out-of-annotation
  fall through to a minimap2 fallback, so uLTRA *cannot beat minimap2* on novel structure —
  "snapped-to-GTF = correct" is assumed, not tested (`redteam_winrates` claim 8). Wins **2%**.

### GMAP — sandwich-DP sequence-based splice model (closest external precedent to Module 2H)
**FACT-PAPER (NOT in the panel; no win-rate).** GMAP brackets each intron between two
anchored exon ends and runs **two DPs simultaneously from both sides**, choosing the single
boundary that **jointly maximizes (base-alignment score + sequence splice strength)**.
Canonical is a **reward, never a gate**: `GT-AG > GC-AG/AT-AC > non-canonical`, with
`--canonical-mode 2` *down-weighting* the canonical bonus on clean reads (won't bend a good
alignment to fake GT-AG) and up-weighting it on noisy reads. Modern builds add **MaxEnt**
novel-splice scoring. **INFERENCE:** this "sequence-first, canonical-as-tiebreaker, slide-to-
best-supported-site" policy is the **closest external analogue to RECTIFY's Module 2H** —
GMAP does it inside the aligner per read; RECTIFY decouples it as post-hoc refinement (and
avoids GMAP's ≥30× speed penalty, which disqualifies GMAP from the production panel).

### gapmm2 — minimap2 body + edlib terminal only (no junction benefit)
**FACT.** Body alignment **is** minimap2 (`-x splice`); gapmm2 only re-places **orphan
terminal exons** via edlib (`mode=HW`), synthesizing `~gtNNNNag`/`~ctNNNNac` intron tokens
biased toward canonical dinucleotides. It touches only ~7% of termini (5'-biased) and adds
**no internal-junction capability** — its junctions are minimap2's junctions. Wins **0.8%**.

### mapPacBio (BBMap) — NO splice model; N-ops only as large indels
**FACT.** `align2.BBMapPacBio` has **no GT-AG / splice-graph model**. Introns are discovered
as **scored long deletions** by its full affine DP (match +100 / mismatch −127 / 2nd-consec
−51), then any `D ≥ intronlen=10` is **reclassified to `N`** afterward. RECTIFY sets
`intronlen=10` and an explicit `maxindel=max(200000, max_intron)` (`multi_aligner.py:749,754`),
so very long yeast-cluster introns stay within the indel cap. Junctions are sequence-optimal
*within the DP window* but
canonically blind. Wins **18.2%** — primarily on unspliced terminal exons (a 3'-end
property, demoted here), not on splice quality.

### Junction-specific comparison table

| Axis (junction) | minimap2 | deSALT | uLTRA | GMAP *(contrast)* | gapmm2 | mapPacBio |
|---|---|---|---|---|---|---|
| **Primary placement criterion** | gap-cost + GT-AG DP | **cross-read consensus** | **annotation snap** | **sequence splice model (sandwich DP)** | inherits minimap2 | DP score then D→N |
| **GT-AG model** | yes (donor/acceptor tables) | yes (`-R 9` noncanon) | snaps to GTF coords | yes (reward, not gate; +MaxEnt) | synth. `gt..ag` token | **none** |
| **Annotation use** | soft `--junc-bonus 9` | de-novo (no `-G`) | **required** | optional IIT | none | none |
| **Cross-read consensus** | no | **yes (batch B=655350)** | pools *annotation*, not reads | no (per-read) | no | no |
| **Novel-junction capable** | yes (jittery) | yes | falls back to mm2 | yes (de-novo + MaxEnt) | yes (terminal only) | yes (as scored gap) |
| **Micro-exon (≤~6–10 nt)** | drops <~13 nt | k=8 local-hash rescues some | **inject seedless annotated exon** | de-novo splice-prob gate | recovers short 5' exon | drops |
| **Signature junction failure** | ±1–few bp jitter | long-range chimeric N | circular on novel | per-read inhomogeneity; slow | no internal-junction value | splice-blind misplacement |
| **Junction win contribution** | 0.1% | **78.9%** | 2% | n/a | 0.8% | 18.2% |

*(All win rates are RECTIFY-internal and single-dataset. They are produced by Path A
(HP-edit-distance, a junction-aware quality metric — see §2). The single-dataset /
un-committed-artifact caution stands.)*

---

## 2. How RECTIFY currently validates / refines / selects junctions

RECTIFY has a genuinely sophisticated **junction-refinement** layer (Module 2H) and a
**validation** layer. Junction quality drives production winner selection via HP-edit-distance
(Path A) — both production call sites pass per-aligner raw BAMs + genome, so `use_hp_ed=True`
and Path A runs (detailed in the CRITICAL FINDING block at the end of this section).

### Module 2H — `refine_read_junctions` (FACT, verified in `junction_refiner.py`)
For every `N`-op in every read, 2H scores candidate junctions and re-snaps the boundary:
- **Sequence-first HP-anchored scoring** (`_score_junction`, `junction_scoring.py`):
  bilateral `score(k) = t1(k) + t2(k)` where `t1` is `rescue[k:]` aligned left-anchored to
  `g[intron_end:intron_end+buf]` (exon2 body) and `t2` is `rescue[:k]` reversed against
  `g[intron_end-buf:intron_end]` reversed (intron-end-proximal). t2 prevents degenerate
  `k=L-1` coincidental single-base matches from scoring 0 (the v3.1.1 `range(L)` fix). The
  DP is `_score_hp_anchored` (Numba-JIT when available, pure-Python fallback otherwise).
- **Candidate discovery / boundary constraint** (FACT): `search_radius=5000` discovers
  candidates (incl. adjacent-isoform annotated junctions); `max_boundary_shift=50` constrains
  each endpoint shift; `_apply_junction_replacement` re-asserts a `_MAX_BOUNDARY_SHIFT=50`
  guard and encodes shifts as I/D ops (preserving ref+query spans), stripping MD/cs.
- **Canonical/annotation as tie-breakers ONLY** (FACT): `_CANONICAL_HP_PRIOR = 0.5` gives
  canonical-tier (`tier<4`) alternatives a 0.5-edit discount **only when the current N-op is
  non-canonical** (`tier_beats_alt`). The tuple order flips on `current_tier>=4`:
  `(score, tier, is_alt, is_novel, |delta|)` vs `(score, is_alt, tier, is_novel, |delta|)`.
  Per the PERMANENT v3.1.7 policy: **no candidate guards** — non-canonical/novel junctions
  are scored and can win on clear sequence evidence (≥1.0 edit difference); annotation/
  canonical never gate.
- **Boundary-error pre-filter** (FACT): `_has_boundary_error` (window=10) skips clean M-op
  junctions unless they are non-canonical (`tier>=4`, always scored). A 255× fast path.
- **Wiring** (FACT): 2H runs only when `rectify correct --aligner-bams` is supplied (plus
  annotation); the candidate pool (`build_junction_pool`) is the **union of annotated +
  per-aligner-BAM `N`-ops** — so the raw material for cross-read recovery exists, but **no
  consensus site or cross-read support count is ever computed** (each read is refined in
  isolation against a static list).

### Empirical `--junction-penalty-table` (FACT — bundled and auto-resolved)
`HpPenaltyTable.from_tsv` (`hp_penalty.py`) loads base-class-split (AT/CG) HP del/ins
penalties; `--str-penalty-table` adds STR slippage costs. When present, these replace the
heuristic step-function (`del(HP=1)=1.0, del(HP≥4)=0.5, ins=1.25`).
- **FACT:** the empirical tables are **bundled** for both organisms with protocol variants:
  `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/{penalty_scores.tsv,
  penalty_scores_cdna.tsv, penalty_scores_cdna_umi1/2/3plus.tsv, penalty_scores_qsrev.tsv,
  str_penalty_scores.tsv, junction_overhang_table.tsv}` and a parallel `homo_sapiens/` set; the
  generator is `scripts/calibration/empirical_cigar_error_profiler.py`. They are **auto-resolved
  by `--Scer` by default** (protocol-routed): `resolve_reference_paths()` fills
  `args.junction_penalty_table`/`str_penalty_table`/`junction_overhang_table` from bundled data
  when empty (`rectify/data/__init__.py:1188-1208`, routed via `get_bundled_junction_penalty_table`
  `:835`), and the global CLI hook calls it (`rectify/cli.py:198-199`). So the empirical tables —
  not the heuristic fallback — are what runs by default.
- **Numba (sub-claim, environment-dependent):** when Numba is not importable, 2H runs the
  pure-Python DP fallback (`REDTEAM_proposals.md` C11). Numba availability is environment-specific
  and independent of the table bundling.

### Chimera-overhang filter (FACT — `calibrate_junction_overhang.py` + `_add_chimera_flag`)
Empirically calibrates `min_overhang(intron_size)` from **concordant multi-aligner reads
(≥4 aligners agree on a junction = silver-standard ground truth)**, isotonic-smoothed
(PAV). At merge time, `_add_chimera_flag` sets `_chimera_ok=1` for aligner rows whose
junction has insufficient flanking overhang for its intron size; chimeric aligners sort
*after* non-chimeric ones (a soft penalty, not a hard drop). Short-intron relaxation
(<500 bp, ≥5 support, max_overhang≥10) allows 1 bp. When no table is supplied,
`OverhangTable.default()` is used (never silently skipped).

### `junction_validator.py` + `false_junction_filter.py` (FACT — separate, supporting)
- **`junction_validator.py`**: a 3-pass COMPASS-style **cross-sample** validator
  (extract → `filter_cross_sample_junctions` requiring `min_samples≥2`, `min_total_reads≥3`,
  canonical motif unless unknown → `apply_junction_filter` **downgrades the `Xc` confidence
  tag**, never discards). This is a separate CLI (`extract/validate/apply-junction-filter`),
  not wired into the core correct-first merge.
- **`false_junction_filter.py`**: flags **poly(A)-artifact** `N`-ops near the 3' end
  (A-rich downstream + genomic A-tract target + non-GT-AG) and absorbs them as large
  deletions — a 3'-end-protection module, demoted here.

### CRITICAL FINDING — junction quality drives production selection
`merge_corrected_tsvs` (`corrected_consensus.py`) has two sort paths, switched by
`use_hp_ed = bool(per_aligner_corrected_bams or per_aligner_raw_bams)` (`corrected_consensus.py:1262`):
- **Path A (HP-edit-distance, production default)** — sort `(read_id, _effective_chimera_ok,
  hp_edit_distance, _span)`. The junction-quality selection criterion.
- **Path B (legacy popularity vote, FALLBACK)** — sort `(read_id, _five_rescued, _chimera_ok,
  _conf_rank, _n_agree, _span, _n_junc)`. `_n_agree` = count of rows sharing
  `(read_id, corrected_3prime)` — a **majority-vote popularity** term.

Both production call sites pass `per_aligner_raw_bams` + `genome`
(`run/single_sample.py:238-250` and `split_command.py:1085-1094`), so `use_hp_ed=True` and a lazy
raw-BAM HP path (`:1264-1364`) computes HP-edit-distance in memory from raw BAM + TSV. **Path A
runs in production.** Path B is the documented fallback only (no BAMs/genome supplied, or staging
fails → `ed_rows` empty → `use_hp_ed=False`, `:1377`).

Therefore:
1. **Junction quality IS a selection key in production** — `hp_edit_distance` (the
   junction-aware quality term) is computed and is the primary non-flag discriminator. (In the
   Path B fallback it is not computed and selection falls back to rescue flag / chimera flag /
   self-reported confidence / 3'-position popularity.)
2. **N-ops cost 0** in `_cigar_hp_edit_distance` (verified: `op == 3: ref_pos += length
   # free pass`; `corrected_consensus.py:142-143`). An aligner that invents a long false intron
   pays **nothing** for the skipped reference and can *lower* its edit distance by replacing
   mismatches with a free `N` — a "free false-intron" exploit. Beyond the binary overhang flag, a
   **junction-anchor gate** backstops this: `_cigar_min_junction_anchor` (`:166-220`) + the
   `_add_chimera_flag` anchor branch (`:868-882`) flag an N-op winner chimeric when its
   junction-local perfect-match anchor is below `min_junction_anchor_bp`. The gate is
   **organism-gated** (`rectify/data/__init__.py:1053-1062`; `_MIN_JUNCTION_ANCHOR=0` at
   `corrected_consensus.py:650`): **human default 10 bp, yeast 0 (gate OFF).** So on **yeast** the
   overhang filter is the only defense and the free-false-intron exploit is live; on **human** the
   anchor gate adds a graded structural defense. The calibrated N-op cost remains the open item
   (J1).
3. **Aligner-style coupling** (HYPOTHESIS): under Path A, an output style minimizing
   HP-edit-distance (clean `=`/`X`, free `N`) can couple to wins independent of biological
   accuracy. On the Path B *fallback*, deSALT's homogeneous junctions → many aligners share
   `corrected_3prime` → high `_n_agree` → herd-bias win, and uLTRA's GTF-snap is
   annotation-circular; `_n_agree` is not a Path A sort key, so the herd-bias mechanism does not
   apply to production runs. The 78.9/18.2/2/0.8/0.1 spread is a **metric output on one dataset,
   possibly pre-v3.3.0-fix** (`index_col`/`_pt:i:N` bugs that stopped 3 aligners competing) —
   **HYPOTHESIS as biology.**

---

## 3. Where junction placement goes wrong

Per-aligner junction failure modes (FACT mechanism / INFERENCE attribution), mapped to the
validation reads that exercise them:

| Failure mode | Aligner(s) | Mechanism | Validation read(s) |
|---|---|---|---|
| **Boundary jitter ±1–few bp** | minimap2, gapmm2 | low `A=1/B=2`; HP-rich junction where exon-end+intron-start ≈ intron-end+exon2-start; no column is pinned | **cat9** (00a1c9b3, 00a1e01e, 0b3b593b, d3357db5 — 2H corrects 4) |
| **Non-canonical drift** | minimap2, mapPacBio | aligner places a locally-optimal split at a non-GT-AG site when a nearby canonical one exists (e.g. `54=` into a GG donor) | **cat9_plus_2** (00a1e01e: GA-GG tier-6 → annotated GT-AG) |
| **Micro-exon miss (≤~6 nt)** | minimap2, mapPacBio, gapmm2 | seed too short to anchor a chain across the intron; exon dropped/soft-clipped. *Real floor ~6 nt, not the "8 nt" headline* (`improve_splice_junctions` F2) | **none** (no cat exercises micro-exon insertion — a validation gap) |
| **Novel-junction noise** | minimap2 | emits jittery novel `N` + free-N-cost exploit (§2) | (covered by overhang filter, not a single read) |
| **Annotation circularity** | uLTRA | hard-snaps to GTF → "correct" on annotated reads by construction; blind to true novel sites | (no novel-junction truth set committed) |
| **Long-range chimeric `N`** | deSALT, gapmm2 | spurious far junction to skip mismatches; the overhang filter exists *specifically* to gate this | **cat7** (XU=1 chimeric: f8050895, 7d5e8dc2, 72557a9a) |
| **5' truncation across a junction** | mapPacBio intrusion, all | aligner ends inside the upstream intron, dropping exon-1; needs 3'SS rescue + local exon CIGAR | **cat3** (5' rescue: 79f61403, 28ea9379; `rescue_3ss_truncation` + `local_aligner.py`) |

**Net:** the failure modes are **complementary** — minimap2/mapPacBio jitter and miss
micro-exons; deSALT/uLTRA recover structure but deSALT chimeras and uLTRA is circular;
gapmm2 adds only terminal exons. This complementarity is the justification for the ensemble
+ 2H union pool — and it is realized **because junction quality drives selection**: junction
quality (HP-edit-distance, Path A) drives production winner selection (§2), so the precondition
for the ensemble paying off is met.

---

## 4. A junction-placement VALIDATION strategy

Since validation is the user's stated primary goal, this section is concrete. The aim is a
**reproducible, committed harness** that measures junction-placement accuracy *independent
of selection*, so "deSALT places junctions better" can be separated from "deSALT wins the
popularity vote."

### Truth signals (ranked by orthogonality)
1. **Annotated GT-AG concordance (FACT, available now).** For every refined `N`-op, compute
   distance to the nearest annotated SGD intron and the donor/acceptor dinucleotide. Metrics:
   (a) % of reads whose junction is *exactly* an annotated boundary; (b) median |boundary
   shift| from annotation; (c) canonical-dinucleotide rate. Caveat: rewards annotation →
   must be paired with a **novel-junction arm** (below) to avoid the uLTRA circularity.
2. **Cross-aligner junction concordance (FACT, available now — the silver standard).**
   Reuse the **≥4-aligner-agreement** definition already in `calibrate_junction_overhang.py`
   (`N_CONCORDANT_DEFAULT=4`). A junction reported by ≥4 of the 5 aligners for the same read
   is treated as ground-truth; measure each aligner's (and 2H's) deviation from it. This is
   the strongest *internally available* signal and already trusted by the overhang calibrator.
   **Caveat (INFERENCE):** correlated-error families (3 of 5 aligners share minimap2 lineage)
   can make "concordance" partly herd-defined — report concordance both with and without the
   minimap2 family.
3. **The cat3/cat7/cat9 validation reads (FACT).** The committed `validation_reads.bam`
   already encodes the three junction scenarios: **cat3** (5' rescue, `five_prime_rescued=1`),
   **cat7** (chimeric, `XU=1` retained), **cat9** (2H junction-refine, 4 reads corrected to
   annotated GT-AG). The junction harness must assert these do not regress and ideally tighten.
   `evaluate_hp_pen_values` (`junction_refiner.py`) is a ready scaffold: it takes a
   `ground_truth_junctions` list and reports TP/FN/true-corrected/false-corrected per `hp_pen`.
4. **NET-seq — note the limits for junctions.** NET-seq is a **3'-end (nascent-CPA) assay**;
   it validates corrected CPAs, *not* internal junction boundaries. It is also partly **baked
   into** the corrected positions (`netseq_refiner.py`), so it cannot serve as an independent
   junction oracle. Use it only for the 3'-terminal-exon boundary, not for splice junctions.
5. **Orthogonal junction truth set (the real gating investment).** What is genuinely missing:
   an **independent intron-coordinate truth set** — e.g. short-read RNA-seq split-read
   junctions (STAR SJ.out.tab) on the *same* biological sample, or curated SGD introns plus a
   spiked **SIRV/sequin** control with known multi-exon transcripts and documented micro-exons.
   Short-read split-read junctions are the field-standard orthogonal call for long-read
   junction accuracy and are not subject to the long-read HP/jitter error mode. Until this
   exists, items 1–3 are sanity checks, not an oracle (`REDTEAM_proposals.md` §4).

### How to measure reproducibly
- **Commit a junction harness** (`dev/junction_eval/`) that, given per-aligner BAMs + a
  truth set, emits a `junction_accuracy.tsv` (per aligner and per pipeline stage: raw → 2H):
  exact-match rate, median |shift|, canonical rate, micro-exon recall, false-novel rate.
- **Commit provenance**: the exact dataset, commit hash, aligner versions, penalty-table
  presence (bundled + auto-resolved by `--Scer`), Numba availability, and which sort path ran
  (Path A by default). This directly answers the `redteam_winrates` "is the spread a metric
  artifact?" question.
- **Path-A/B ablation** (zero-oracle, do first): run the same TSVs through
  `merge_corrected_tsvs` with and without `per_aligner_corrected_bams`/`per_aligner_raw_bams`
  and diff `winning_aligner` counts. Both are live paths (Path A default, Path B fallback). If
  the spread moves, the win rates are a selection-metric artifact — settled with no truth set
  needed.

---

## 5. Prioritized junction-placement improvements

Re-ranked **junction-first** across `improve_splice_junctions.md` → `REDTEAM_proposals.md`.
3'/CPA items are **explicitly demoted** at the end. Each: mechanism · failure fixed · effort ·
risk · red-team verdict · caveat.

### J0 (do first, ~free, zero-oracle) — Selection provenance + Path-A/B ablation
- **Mechanism:** confirm whether the win rates are post-v3.3.0-fix and which sort path
  produced them; commit `aligner_summary.tsv`. **Fixes:** the §2 measurement blindness.
  **Effort:** days. **Risk:** none. **Red-team verdict:** KEEP — "the single correct do-first"
  (`REDTEAM_proposals` E0). **Caveat:** needs no truth set; gates every value claim below.

### J1 (HIGHEST) — Calibrated N-op open cost (the BAM wiring already exists)
- **Mechanism:** charge each `N`-op a **small calibrated open cost** = −logP(novel intron
  of this length at this site) instead of 0, in `_cigar_hp_edit_distance`. The BAM-wiring half
  is already done: both production call sites pass `per_aligner_raw_bams` + `genome` and a lazy
  raw-BAM HP path (`corrected_consensus.py:1264-1364`) runs Path A without materializing corrected
  BAMs (`single_sample.py:238-250`, `split_command.py:1085-1094`), so HP-edit-distance already
  ranks aligners. **Only the calibrated N-op open cost remains open.** The human-only
  junction-anchor gate already backstops the free-`N` exploit structurally (see §2); a calibrated
  N-cost would generalize that to yeast.
- **Failure fixed:** the free-false-intron exploit (§2). **Effort:** low (the wiring is done;
  only the N-cost calibration remains). **Risk:** medium — must validate cat7/cat9 before
  default-on. **Red-team verdict:** KEEP-WITH-CAVEATS (the most-emphasized correctness gap).
  **Caveat:** the N-cost must be *calibrated*, not a fixed penalty — a fixed penalty re-introduces
  a gate, violating the PERMANENT no-gate policy.

### J2 — Validate / version the bundled empirical junction penalty tables
The tables are bundled (per-organism, protocol-routed) and **auto-resolved by `--Scer`** (see §2:
`rectify/data/genomes/*/penalty_tables/`, `data/__init__.py:1188-1208`, `cli.py:198-199`). J2 is
therefore to **validate / version the already-bundled tables**, not to regenerate absent ones.
- **Mechanism:** confirm the bundled `penalty_scores.tsv` / `str_penalty_scores.tsv`
  (and protocol variants `penalty_scores_cdna*.tsv`, `penalty_scores_qsrev.tsv`) match the
  current chemistry/sample; record their provenance + generation command
  (`scripts/calibration/empirical_cigar_error_profiler.py`); pin a version so a recalibration is
  auditable. **Failure fixed:** silent drift between the bundled tables and the run's
  chemistry/protocol (the tables are R10.4.1 + *S. cerevisiae*-specific and must not transfer to
  HiFi/other organisms). **Effort:** low–medium (no regeneration needed unless the bundled
  tables are found stale; if they are, regeneration is gated on a working aligner panel + a
  labelled run — deSALT fails on some chunks). **Risk:** low. **Red-team verdict:** KEEP.
  **Caveat:** tables are R10.4.1 + *S. cerevisiae*-specific; verify protocol routing
  (drs/cdna/qsrev) picks the intended variant.

### J3 — Cross-read junction consensus as a first-class step (generalize deSALT)
- **Mechanism:** insert a `JunctionConsensus` build between alignment and per-read 2H: pool
  every `N`-op across **all 5 aligner BAMs + annotation**, aggregate `(chrom,strand,donor,
  acceptor) → n_reads, n_aligners, canonical_tier, mean_overhang`; cluster within
  `max_boundary_shift=50`; elect a **consensus site by a calibrated support score** (not raw
  count — the herd trap); re-snap each `N`-op to its cluster consensus **with the existing
  `_score_junction` as a per-read veto** (move only if the read's own sequence supports the
  consensus ≥ its current placement). **Failure fixed:** minimap2/gapmm2/mapPacBio per-read
  inhomogeneity (the exact trait deSALT wins on, now available to all five and de-biased).
- **Effort:** medium (primitives exist: `build_junction_pool`, `_candidates_near`, fork-shared
  worker pool). **Risk:** medium-HIGH — **CONFLICT-3** in `REDTEAM_proposals`: over-snapping
  can collapse genuine alt-5'/3'SS isoforms <50 bp apart (real in yeast RP genes) and is
  self-reinforcing with `_n_agree` herd bias. **Red-team verdict:** KEEP-WITH-CAVEATS. **Caveat:**
  the per-read sequence veto + shared-boundary-anchoring (cluster only junctions sharing one
  exact boundary) are mandatory, not optional; validate against truth, never internal agreement.

### J4 — Unified calibrated junction scorer with positive splice-strength (MaxEnt/GeneSplicer)
- **Mechanism:** replace 2H's brittle lexicographic tuple with a calibrated log-odds
  `S(r,J) = w_seq·logP_seq + w_don·logP_donor + w_acc·logP_acceptor + w_xr·log(1+n_aligners)
  + w_ann·1[annotated]`, weights fit by logistic regression on ≥4-aligner-concordant
  positives. **`logP_donor/acceptor` = MaxEntScan/GeneSplicer** — RECTIFY has **NO positive
  splice-strength term today** (only canonical/non-canonical *tiers*); minimap2's
  `--splice-flank=no` deliberately discards the flanking-strength prior and hands it to
  RECTIFY, which then doesn't model it. MaxEnt fills exactly that gap, and is the model GMAP's
  modern builds adopted. **Failure fixed:** the un-comparable lexicographic score (can't feed
  selection) + the missing splice-strength signal. **Effort:** medium (yeast PWM from SGD is
  small; the logistic fit needs sklearn, **not importable here**). **Risk:** medium — MaxEnt is
  metazoan-trained → **must retrain a yeast PWM**; `w_ann` must stay small (annotation
  circularity). **Red-team verdict:** KEEP-WITH-CAVEATS, gated on J1 wiring. **Caveat:** only
  helps "if wired into selection" — same dependency as J1; build the PWM/MaxEnt backend behind
  a `splice_strength(chrom,pos,strand)` interface so a metazoan SpliceAI/Pangolin CNN is a
  drop-in later (defer the CNN — NIL value on yeast's ~300 canonical introns).

### J5 (LOW for yeast) — De-novo micro-exon recovery (≤~6–10 nt)
- **Mechanism:** GMAP-style — when two consensus exons flank an `N`-op with unexplained
  boundary error, probe a ≤10 nt exon inside the intron; accept only if both flanking splice
  sites clear a MaxEnt splice-prob gate (J4) AND cross-read support (J3) ≥ N. **Failure fixed:**
  micro-exon misses, **de-novo** (works on novel micro-exons uLTRA's annotation-injection
  can't). **Effort:** medium-low (needs J4). **Risk:** medium-high (spurious micro-exons →
  false isoforms; double gate essential). **Red-team verdict:** DESCOPE → metazoan. **Caveat:**
  yeast has few micro-exons; **no validation read exercises this** (a new cat10 + SIRV control
  would be needed) — lowest-validated, defer to the metazoan port.

### Explicitly DEMOTED (3'/CPA — already solved; out of scope for junction placement)
- **Poly-A-aware terminal CPA solver / two-state G→T DP / terminal re-scoring** (chaining P3,
  3prime P1a, chaining P2): a 3'-end boundary solver — **NOT a junction improvement**, and a
  high-regression rewrite of the most-patched code path (`REDTEAM` CONFLICT-1). Demoted.
- **`pt:i` poly-A prior in 2E walk-back** (ml §2a): a 3'-end signal. Demoted.
- **NET-seq selection prior / CPA-anchor seeding** (3prime P3, seeding 5a): 3'-end + circular
  oracle risk. Demoted. **mispriming terminal veto** (3prime P5): a cDNA 3'-end artifact. Demoted.
- **GBDT/transformer/foundation-model selectors**: REJECT/advisory (overfit on 36 reads,
  libs absent); not a junction-placement mechanism. Demoted.

---

## If you do only three things for junction placement

1. **J0 + J1 — make junction quality decide selection well.** Run the zero-oracle Path-A/B
   ablation and commit `aligner_summary.tsv`; then **fix the N-op-cost-0 free-intron exploit**
   so `hp_edit_distance` ranks aligners robustly. The BAM wiring is already done (lazy raw-BAM
   path; Path A is the production default), so junction quality already ranks aligners; the
   remaining J1 work is the **calibrated N-op cost** — relevant for yeast, where the anchor gate
   is off by default.
2. **J3 — cross-read junction consensus across all five aligners**, de-biased by a calibrated
   support score and protected by a per-read sequence veto + a hard APA-preservation radius.
   This generalizes the single mechanism that makes deSALT the leader, with the over-
   homogenization risk explicitly bounded.
3. **J4 (+ J2) — add the splice-strength signal RECTIFY lacks.** Build a yeast MaxEnt/PWM
   donor/acceptor model behind a swappable interface and fold it (plus the already-bundled,
   `--Scer`-auto-resolved empirical penalty tables — see §2/J2) into a single calibrated junction
   scorer that feeds selection — replacing the brittle lexicographic tuple with a comparable,
   calibrated quantity. **Validate
   everything against ≥4-aligner concordance + cat3/cat7/cat9 + an orthogonal short-read
   split-read junction set — never against internal agreement alone.**
