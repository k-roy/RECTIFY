# Making Alignment Itself CPA-Accurate — Discovery Proposals

> Verified against `origin/drs-validation-rebuild` @ 366c885 (2026-06-19).

**Goal.** bp-accurate 3'-end / cleavage-and-polyadenylation (CPA) site determination from ONT
direct-RNA (DRS) / dT-primed cDNA long reads, by making the **aligner / DP terminal behavior**
report a CPA-correct `reference_end` *at alignment time*, reducing reliance on post-hoc correction
(Modules 2B/2C/2E/2G/2H, DRS Steps 0/4).

**Author role.** Algorithmic Discovery — ambitious but rigorous. Each proposal is marked
**ESTABLISHED** (the mechanism or a close analogue exists in published tools/literature) or
**NOVEL** (no surveyed tool does this; it is a new design).

**Sources.** `02_synthesis/{COMPARISON.md,DEEP_DIVE.md}`; `01_investigation/{ont_drs_ecosystem,
pacbio_ecosystem}.md`; `03_adversarial/{redteam_annotation_ecosystem,redteam_denovo,
redteam_winrates_selection}.md`; `CLAUDE.md` (Modules 2B/2C/2E/2G/2H, DRS Step 0/4,
`--dT-primed-cDNA` vs DRS, the AT>CG penalty table); RECTIFY source (`multi_aligner.py`,
`indel_corrector.py`, `walkback.py`, `junction_refiner.py`, `netseq/netseq_refiner.py`,
`polya/polya_trimmer.py`, `drs_trim_command.py`).

---

## 0. Framing, caveats, and the honest baseline

Three caveats from the adversarial pass constrain every claim below and must travel with this doc:

1. **The AT>CG deletion-rate claim is RECTIFY-internal and unproven.** External literature confirms
   A/T homopolymers are *basecalled more accurately overall* than C/G; it does **not** independently
   corroborate the AT>CG *deletion-specific* direction in `penalty_scores.tsv`. The "ratcheting-speed"
   mechanism is a hypothesis (redteam O6, "circular"). Any proposal that *reuses* that table inherits
   its uncertainty; proposals that re-derive penalties from signal/orthogonal data are *more* robust.

2. **"No mainstream tool corrects per-read 3' ends" must be framed as triad-uniqueness, not an
   absolute** (redteam O5). The defensible claim: *no surveyed tool combines (a) multi-aligner
   ensemble + (b) per-read 3'-end correction against an HP error model + (c) correct-first selection.*
   Signal tools (nanopolish polya/tailfindr/dorado) call **tail length**, not genomic CPA coordinate;
   SQANTI3/LAPA flag/cluster, they do not place per-read CPA. The gap below is **CPA-aware terminal
   alignment**, which genuinely no aligner does.

3. **The win-rate spread (deSALT 78.9 / mapPacBio 18.2 / …) is a single-dataset number, not a
   measured per-aligner CPA-accuracy ranking.** Production runs the live HP-edit-distance selector
   (`use_hp_ed=True`; both call sites pass per-aligner BAMs + genome), which ranks on a
   genome-anchored edit distance over the corrected CIGAR — not the legacy popularity sort. But that
   score still does not *directly* sort on 3'-end position (it scores CIGAR fit, with N-ops free), so
   a win does not establish CPA accuracy. **Consequence for THIS document:** we must not assume
   "deSALT is the most CPA-accurate backbone." Every proposal that touches selection or claims an
   accuracy improvement must be validated against an **orthogonal CPA truth set** (NET-seq,
   QuantSeq/3'-seq, validation cat1/cat2), independent of the selection metric. This is not optional
   rigor — it is the precondition for any of these proposals being more than a hypothesis.

**Why move correction into alignment at all?** Post-hoc correction is provably good at *length-
independent* fixes (2E walk-back uses the **reference** genome, so HP miscount doesn't matter). But it
is fundamentally limited where the aligner has *already discarded the evidence*: minimap2 soft-clips
the poly-A and `reference_end` drifts; deSALT/GMAP slip the terminal exon into a genomic A-run; uLTRA
*strips* poly-A before aligning (no CPA at all). Once the A-run is soft-clipped or slid, the only
recoverable signal is the **genomic** A-context (what 2E/2G use) — the *read's own* tail-boundary
information is gone. A CPA-aware terminal DP keeps that information in the alignment, where it is
strongest, and gives correction a better starting point (or makes it a near-no-op). That is the thesis.

---

## 1. Poly-A-aware terminal alignment mode (`--end-bonus`-style CPA DP)  — **NOVEL** (core)

### Mechanism
Replace the aligner's terminal-exon DP behavior with a **CPA-anchored terminal mode** that explicitly
models the read as `[exon body][poly-A tail][adapter/soft-clip]` and places the reference_end at the
**body→tail transition**, not at the first base where the A-run stops scoring. Concretely, a custom DP
terminal recurrence (or a post-chain terminal re-scoring pass over the existing aligner BAM) that:

- knows the read carries an in-read poly-A tail of estimated length `pt:i:N` (dorado), and
- scores the terminal columns with an **asymmetric A-run model**: extending the alignment through a
  genomic A-run is penalized (it is *tail*, not *body*), while the last *non-A genomic base before the
  A-run* is rewarded as the CPA candidate. This is the algorithmic inverse of minimap2's `--end-bonus`:
  instead of "reward extension to the read end," it is "reward stopping at the body/tail boundary."

This is the **alignment-time generalization** of Module 2E (`find_atract_boundaries` /
`find_polya_boundary` in `indel_corrector.py`): 2E already walks back over a genomic A-run to the last
non-A position using the reference. Proposal 1 performs that walk-back *as a terminal DP state* so the
emitted `reference_end` is the CPA directly, and a `tail_softclip` carries the A-run + adapter.

**Two feasible implementation forms (ranked by effort):**
- **(1a) Terminal re-scoring pass** (low effort, recommended first): a post-alignment module that takes
  each aligner's primary BAM record and re-places the 3' boundary using the genome A-context + the
  `pt:i` length + the soft-clipped read tail. This is ~90% of what `find_polya_boundary` +
  `rescue_softclip_at_homopolymer` (2G) already do — but run as a **unified terminal solver inside the
  aligner wrapper** (`multi_aligner.py`), before TSV correction, so *every* aligner's BAM enters
  correction already CPA-snapped. Net effect: collapses 2B/2E/2G into one boundary decision and removes
  the "opposite-direction corrections cancel" problem CLAUDE.md notes between 2E and 2G.
- **(1b) True DP terminal state** (high effort): patch a single aligner's DP. minimap2's KSW2
  (`ksw_exts2_sse`) and BBMap's affine DP are the candidates; minimap2 already exposes `--end-bonus`,
  `-z`/`zdrop`, and a flank model — a CPA terminal state is a natural extension. Higher risk (touches
  C/Java internals, must be re-verified per release).

### Specific failure mode fixed
- **minimap2** `--end-bonus 0` soft-clips poly-A → `reference_end` drifts; HP indel noise distributes
  arbitrarily in the A-run (low A=1/B=2 doesn't pin a column). *(DEEP_DIVE §7, minimap2)*
- **deSALT / GMAP / mapPacBio** slip the terminal-exon end into a downstream **genomic** A-run because
  they align the in-read A-tail as ordinary sequence. *(DEEP_DIVE §7)*
- **uLTRA** strips poly-A → no CPA candidate at all. *(DEEP_DIVE §4)*

### Feasibility
**High for (1a)** — it reuses `indel_corrector.find_polya_boundary`, `atract_detector`, and the genome
loader, and the `pt:i` tag is already parsed in `drs_trim_command.py:411`. The novelty is *where* it
runs (per-aligner, pre-correction, unified) not *what* it computes. **Medium-low for (1b).**

### Expected impact on CPA bp-accuracy
**Highest of all proposals.** This directly attacks the dominant residual error class (HP deletion at
the A-tract, DEEP_DIVE §7). Because the genomic-context walk-back is **length-independent**, it is
robust to the exact HP-miscount magnitude. Plausible to move the *raw* aligner 3' end from "noisy by
several bp" to "within 0–1 bp of the 2E-corrected position" *before* correction — i.e. it makes the
aligner output what correction currently has to manufacture.

### Validation
- **NET-seq concordance:** the in-pipeline gold standard (`netseq_refiner.py`, bundled
  `saccharomyces_cerevisiae_netseq_wt.tsv.gz`). Compare distance-to-nearest-NET-seq-peak distribution
  for raw aligner ends vs (1a)-snapped ends vs current 2E-corrected ends. Success = (1a) matches
  2E-corrected to within the NET-seq peak width on a large random sample.
- **Validation cat1/cat2:** cat1 reads exercise `atract_ambiguity,indel_correction,polya_walkback`
  (shifts −16/−4/+3/+3 bp); cat2 exercise `softclip_rescue` (+11/+4/−8 bp). (1a) must reproduce these
  exact corrected positions **from the raw alignment alone**, with no separate 2E/2G pass.

### Risk
- **Over-correction into real genomic A-runs** at true internal-priming sites (the `--dT-primed-cDNA`
  problem) — mitigated by Proposal 5's mispriming gate.
- **Cancellation with 2G** (CLAUDE.md already documents 2G-takes-priority-over-2E); a *unified* solver
  removes this by deciding the boundary once. **This is a feature, not a risk, of (1a).**
- **DRS Step 0/4 interaction:** under `--drs` the tail is already trimmed pre-alignment and restored
  post-correction. (1a) must run on the *trimmed* representation and hand off to Step 4 cleanly — i.e.
  it refines the trim boundary, it does not re-introduce the tail. Low risk, but must be wired between
  Step 0 and correction.

---

## 2. Signal-level / basecaller-integrated CPA fusion (`pt:i` + dwell → genomic CPA)  — **NOVEL** (high value)

### Mechanism
Fuse the **raw-signal / basecaller poly-A boundary** with the alignment so the genomic CPA is pinned by
two independent estimators. dorado `--estimate-poly-a` already emits `pt:i:N` (tail length in nt;
`pt:i:0` = anchor found/length unknown; `pt:i:-1` = no anchor) — RECTIFY stores these in parquet
(`drs_trim_command.py`). Signal tools (nanopolish polya, tailfindr, NanoTimer) locate the
**body→tail dwell-time transition** in the squiggle, which is *independent of basecalled sequence* and
therefore not corrupted by HP-deletion basecalling.

The fusion: the in-read poly-A tail occupies read coordinates `[q_body_end, q_body_end + N)`. The
signal/`pt` estimate fixes `q_body_end` (where the A-tail *starts* in the read) more reliably than the
basecalled A-run, which loses bases. Projecting `q_body_end` through the CIGAR to the genome gives a
**signal-anchored CPA** that can be used as:
- a **prior/anchor** in Proposal 1's terminal DP (constrain the body/tail boundary to agree with `pt`),
  and/or
- a **per-read confidence weight**: reads where signal-CPA, basecalled-CPA, and aligner-end agree get
  high confidence; disagreements flag ambiguous CPAs (a calibrated quantity, unlike the emitted
  self-assessed `confidence`/`_conf_rank` flag, which carries no notion of the HP-ED margin the live
  selector used).

### Specific failure mode fixed
The *root cause* in ONT DRS: a long poly-A produces a low-variance current segment whose **length is in
the dwell time, not the basecall** (`ont_drs §1.3`). The basecaller drops A's; the signal does not. No
panel aligner uses signal at all — they all consume the lossy basecalled sequence. This is the one
proposal that adds **genuinely new information** rather than re-using the genome.

### Feasibility
**Medium.** `pt:i` fusion is **easy and high-value** (the tag is already in the pipeline; just project
it). Full raw-signal (`.pod5`/dwell) integration is **heavier** — requires the squiggle, a move-table
or `tailfindr`/`nanopolish` dependency, and signal↔reference projection. **Recommended phasing:**
Phase A = `pt:i`-only fusion (cheap, mostly plumbing); Phase B = dwell/move-table fusion (research-grade,
gated on Phase A showing the basecalled boundary is the bottleneck). The adversarial caveat: independent
benchmarking finds **dorado the preferred poly(A) estimator** (`ont_drs §1.4`), so Phase A may capture
most of the available signal benefit at a fraction of Phase B's cost — do Phase A, measure, then decide.

### Expected impact
**High, and orthogonal to Proposal 1.** Proposal 1 uses *genomic* A-context; Proposal 2 uses *read/signal*
tail-length. They are independent error channels — combining them is exactly the "two independent
estimators agree → trust" logic that makes correct-first robust. Largest marginal gain on reads where the
genomic context is itself ambiguous (CPA *inside* or adjacent to a genomic A-run), which is precisely
where Proposal 1 alone is weakest.

### Validation
- NET-seq concordance, stratified by `pt` agreement: reads where signal-CPA agrees with aligner end
  should show tightest NET-seq concordance; the proposal's value = lifting the *disagreement* reads.
- cat1 reads carry `polya_walkback` + known shifts; check that `pt`-projected `q_body_end` predicts the
  same boundary the genomic walk-back finds (cross-validation of the two estimators on ground truth).
- Synthetic-RNA poly-A benchmark (cited in `ont_drs`, PMC12406214) as an external length-truth check.

### Risk
- `pt:i:0` / `pt:i:-1` reads (no usable estimate) — must fall back to Proposal 1 (genomic) gracefully.
- **cDNA caveat:** under `--dT-primed-cDNA` the tail is **not in the read**, so `pt`-tail fusion does
  not apply; signal fusion is a **DRS-only** enhancement. Must be gated on protocol (parallels the
  existing module-activation table).
- Phase B adds a heavy dependency and `.pod5` retention cost — do not commit to it without Phase A data.

---

## 3. Orthogonal-data-guided alignment / selection (NET-seq / 3'-seq CPA priors)  — **ESTABLISHED analogue, NOVEL placement**

### Mechanism
Use bundled **NET-seq** (and, for cDNA, QuantSeq/3'-seq) signal as a **soft CPA prior** that influences
*alignment-stage terminal placement and aligner selection*, not only post-hoc refinement. RECTIFY
already does NET-seq refinement *after* correction (`netseq_refiner.py`, NNLS deconvolution). The
proposal moves a **soft** version of this signal *upstream*:
- **Terminal-boundary prior:** when Proposal 1's terminal solver has two near-equal CPA candidates
  (the HP-ambiguous case), break the tie toward the candidate nearest a NET-seq peak — exactly the
  "annotation/canonical as tie-breaker, never a gate" policy CLAUDE.md mandates for Module 2H, applied
  to CPA instead of splice junctions.
- **Selection prior:** add **distance-to-orthogonal-CPA** as a tiebreaker term in the live HP-edit-distance
  selector (and the fallback `_n_agree`/`_conf_rank` tiebreakers). The HP-ED score ranks CIGAR fit and
  leaves N-ops free, so it does not directly reward landing on the true 3'-end position; an external
  3'-end anchor supplies the missing CPA-correctness signal. This is the field's documented gold standard
  (`ont_drs §5.1`: 3'-Seq/Iso-Seq quantify PAS more reliably than computational tools).

### Specific failure mode fixed
- HP-ambiguous CPA ties that Proposals 1/2 cannot resolve from the read alone.
- **The selection-metric gap:** the HP-ED score does not sort on 3'-end position (and the fallback's
  `_n_agree`/`_conf_rank` are popularity/self-report). An orthogonal anchor turns the tiebreaker into a
  truth-correlated quantity.

### Feasibility
**Medium.** The signal loader, deconvolution, and per-position query already exist
(`netseq_refiner.NetSeqSignal.get_signal`). The work is wiring it as a **soft prior** into the terminal
solver and the merge tiebreaker rather than a post-hoc snap. **Critical guardrail:** it must be a
*tie-breaker within the noise floor*, never a *gate* — otherwise it manufactures the same
annotation-circularity the red team flags for uLTRA (rewarding "matches the prior" over truth). Use the
same `_CANONICAL_HP_PRIOR = 0.5`-style bounded discount that Module 2H uses.

### Expected impact
**Medium on bp-accuracy directly** (most reads are already resolved by Proposals 1/2); **high on
selection correctness and on the *validity of the win-rate claim itself*.** This is the proposal that
gives the HP-ED selector a CPA-position term and turns the win-rate numbers from "single-dataset and
selector-bound" into "orthogonally-validated CPA accuracy."

### Validation
- **Leave-one-out NET-seq:** hold out a NET-seq sample (use `_pan` to guide, validate against `_wt`, or
  vice versa) so the prior and the truth are not the same file — avoids the circularity the red team
  warns about. Concordance must improve *on the held-out signal*, not the guiding one.
- cat reads: NET-seq tie-break must not move any cat1/cat2 read off its known-correct corrected position.

### Risk
- **Circularity** (using NET-seq to both guide and validate) — mitigated by held-out validation.
- **Organism/condition specificity:** bundled NET-seq is yeast WT; for other organisms/conditions the
  prior is unavailable or wrong. Must be strictly optional and protocol/organism-gated, with a clean
  no-prior fallback. **Never let the prior override read evidence** (gate vs tie-break is the whole risk).

---

## 4. Explicit DRS 3'→5' directionality & 5' motor fall-off handling  — **ESTABLISHED facts, NOVEL alignment treatment**

### Mechanism
DRS sequences **3'→5'**, starting at the poly-A end; the motor releases **~10–12 nt before the 5'
terminus**, truncating the 5' end (`ont_drs §2.1–2.2`). This has two alignment-stage consequences RECTIFY
currently handles only *partially*:
- **3' end is the high-confidence, well-anchored end; 5' end is `rescue`, not `trust`.** Make this
  explicit in the terminal solver: apply **strong** CPA anchoring at the 3' end (Proposals 1/2) and
  **asymmetric, lenient** end-gap handling at the 5' end (do not penalize the missing ~11 nt; do not let
  a truncated 5' exon drag the whole alignment). Cat3 rescue (`splice_aware_5prime.py`, `local_aligner.py`)
  already recovers truncated 5' splice boundaries — this proposal makes the **asymmetry a first-class DP
  property** (3' tight / 5' loose) rather than two unrelated modules.
- **Strand correctness from protocol, not inference.** minimap2 already gets `-uf` (forward-only, DRS is
  stranded). Extend this discipline: feed the known 3'→5' polarity to *every* aligner's terminal model so
  the CPA end is always the *sequencing-start* (high-quality) end and the 5' is always the
  *motor-fall-off* (low-quality) end. This prevents a terminal solver from spending precision budget on
  the end where the data is structurally absent.

### Specific failure mode fixed
- Treating both ends symmetrically wastes correction effort on the 5' end (where ~11 nt are simply
  missing) and can let 5' truncation noise perturb 3' placement. *(ont_drs §2.2)*
- 5' truncation generating spurious short-exon drops that the current Cat3 machinery rescues read-by-read
  rather than as a systematic DP prior.

### Feasibility
**Medium-low.** Much of the 5' machinery exists (Cat3, `local_aligner` affine-gap exon CIGAR). The new
work is **encoding the 3'-tight/5'-loose asymmetry explicitly** in the terminal solver and ensuring the
protocol polarity reaches every aligner. Mostly a consolidation + correctness proposal, not new algorithm.

### Expected impact
**Low-to-medium on 3' bp-accuracy directly** (the 3' end is already the anchored end), **medium on
robustness** — it prevents 5' noise from contaminating 3' placement and removes a class of redundant
correction. Its real value is *defensive*: it stops the terminal solver from being misled by the
structurally-missing 5' end.

### Validation
- cat3 reads (`five_prime_rescued=1`, `five_prime_exon_cigar`) must remain correctly rescued.
- Check that 3' CPA placement is **invariant** to 5' truncation length (inject synthetic 5' truncation,
  confirm CPA does not move) — the core asymmetry guarantee.

### Risk
- Low. Over-leniency at 5' could mask real 5' biology (TSS) — but DRS 5' is acknowledged unreliable, so
  this is consistent with the platform. Keep 5' coordinates flagged as low-confidence (already implied by
  Cat3 being *rescue*).

---

## 5. Principled internal-priming / AG-mispriming treatment at the alignment stage  — **ESTABLISHED concept, NOVEL alignment-time placement**

### Mechanism
For `--dT-primed-cDNA` (QuantSeq / oligo-dT), the oligo-dT primer misprimes onto **genomic A-rich runs**,
producing **false CPA sites** with no poly-A tail in the read (the A-content is *downstream of the
reported 3' end, off-read*). RECTIFY's AG module already detects this **post-correction** by querying the
**genome** (`get_downstream_sequence()`), not the read (`CLAUDE.md`, `ont_drs §2.3`). SQANTI3 does the
analogous downstream-A check but *post-hoc, per-isoform* (`pacbio §5.3`).

The proposal: make the downstream-genomic-A check a **terminal-alignment gate / penalty**, so the
aligner (or the terminal solver of Proposal 1) **refuses to anchor a CPA** at a position whose downstream
genomic window is A-rich (≥6 consecutive A, or downstream A-fraction >50% — the documented internal-
priming thresholds, `ont_drs §2.3`) **unless** there is independent support (a splice junction nearby, or
a `pt`/signal tail under DRS). This is the mispriming-aware inverse of Proposal 1: Proposal 1 *snaps to*
the body/tail boundary; Proposal 5 *vetoes* boundaries that are mispriming artifacts.

Crucially this must be a **penalty/flag, not a hard gate** (same policy as canonical junctions): a real
CPA can legitimately sit upstream of a genomic A-run. The gate fires only to **break ties away from**
A-rich-downstream candidates and to **flag** suspect CPAs, mirroring SQANTI3's intra-priming flag but at
per-read alignment time instead of per-isoform post hoc.

### Specific failure mode fixed
False 3' ends from oligo-dT mispriming on genomic A-runs — a "substantial fraction" of long-read
transcriptome data (`ont_drs §2.3`). Currently caught only after correction; catching it at terminal
placement prevents the false CPA from ever being selected as a winner and from polluting the count matrix.

### Feasibility
**High.** The detector exists (genomic downstream-A query, AG module); the work is invoking it as a
terminal-placement penalty under `--dT-primed-cDNA`. Reuses the existing protocol gate (module is already
off for DRS, on for cDNA — `CLAUDE.md` activation table).

### Expected impact
**High for cDNA / short-read QuantSeq specifically** (the primary `--dT-primed-cDNA` use case, per
`CLAUDE.md`); **zero for DRS** (no priming step — correctly disabled). Targeted but decisive on its
protocol.

### Validation
- QuantSeq REV reads at known internal-priming loci (downstream A >50%): the terminal solver must *not*
  place a CPA there without splice/tail support.
- A-tract walk-back (2E) interaction: 2E is enabled for short reads (CSP internal priming, `CLAUDE.md`);
  ensure the mispriming veto and 2E walk-back do not double-correct or cancel (same class of interaction
  as 2E/2G — decide the boundary once, per Proposal 1).

### Risk
- **False positives at genuine A-rich-context CPAs** (real poly-A signals do occur upstream of genomic
  A). This is exactly why it must be a bounded penalty + flag, never a hard reject — and why DRS reads
  (which carry the real tail) must be exempt. Mirrors the canonical-tie-breaker discipline already
  validated in Module 2H.

---

## Ranking by impact / effort

| Rank | Proposal | Impact on CPA bp-accuracy | Effort | Est/Novel | One-line rationale |
|---|---|---|---|---|---|
| **1** | **P1a — unified terminal re-scoring pass** | **Highest** | **Low** | NOVEL (placement) | Reuses 2E/2G logic as one per-aligner pre-correction CPA solver; attacks the dominant HP-A-tract error; makes raw ends ≈ corrected ends. |
| **2** | **P2 Phase A — `pt:i` signal fusion** | High (orthogonal) | Low-Med | NOVEL | Adds genuinely new (signal) information the basecall loses; cheap because `pt` is already in the pipeline; gates ambiguous CPAs better than the emitted self-confidence flag. |
| **3** | **P5 — mispriming terminal veto (cDNA)** | High (cDNA only) | Low-Med | ESTABLISHED concept / NOVEL placement | Stops oligo-dT false CPAs before selection; reuses the existing genomic-downstream-A detector. |
| **4** | **P3 — NET-seq/3'-seq soft prior in selection** | Med (bp); High (selection validity) | Medium | ESTABLISHED analogue | Adds a truth-correlated CPA-position term to the live HP-ED selector (which leaves N-ops free and does not sort on 3'-end position); converts win-rate claim from selector-bound to validated. |
| **5** | **P4 — 3'-tight/5'-loose DP asymmetry** | Low-Med (defensive) | Med-Low | NOVEL treatment | Stops 5' motor-fall-off noise from perturbing 3' CPA; consolidates Cat3 into a DP prior. |
| (later) | **P1b / P2 Phase B** — true DP terminal state / dwell fusion | Highest ceiling | High | NOVEL | Research-grade; do only after the low-effort versions quantify the remaining gap. |

**Recommended program.** Do **P1a first** (it is mostly a consolidation of code that already exists,
moved upstream and unified). Then **P2 Phase A** (`pt` projection) as an orthogonal anchor and a better
confidence signal. **P5** in parallel for the cDNA path. **P3** next — it is the proposal that makes the
accuracy claims *provable* (held-out NET-seq), which the red team shows is currently the project's
weakest link. **P4** as a correctness/robustness consolidation. Defer the high-effort DP-internals
versions (P1b, P2 Phase B) until P1a/P2A measurements show they are worth it.

**The unifying design principle** (consistent with Module 2H's validated policy): every CPA decision is
**sequence/signal-first**; genomic A-context, `pt`/dwell, NET-seq, and canonical-tail priors are
**bounded tie-breakers, never gates**. This is the same stance that keeps Module 2H from manufacturing
false junctions, applied to the 3' terminus — and it is the stance that survives the adversarial
critique, because it never lets a prior override read evidence.

---

## Validation harness common to all proposals

1. **Orthogonal CPA truth set** (the precondition the red team demands): held-out NET-seq (`_wt` vs
   `_pan`), plus QuantSeq/3'-seq for cDNA. Measure **per-aligner corrected-CPA accuracy independent of
   selection** — this is the experiment that separates "more accurate" from "won the metric"
   (redteam_winrates experiment 3).
2. **Validation cat1/cat2** (and cat3 for P4): each proposal must reproduce the known corrected positions
   and shifts (cat1 −16/−4/+3/+3; cat2 +11/+4/−8) **from the raw alignment**, ideally collapsing the
   separate 2B/2E/2G/2H passes into the terminal solver.
3. **No-regression gate:** the full 708-test suite (incl. `test_validation_reads.py`,
   `test_junction_refiner.py`, `test_bam_writer.py`) must stay green; CPA shifts that change cat
   expectations require explicit, justified test updates (per the v3.2.x precedent in `CLAUDE.md`).
4. **Selector calibration baseline** (redteam_winrates experiment 2): any selection-touching proposal
   (P3) must measure against the live HP-edit-distance selector that actually runs — including its known
   biases (N-ops cost 0; no explicit 3'-end term) — so the metric being improved is the production metric,
   not the legacy fallback sort.
