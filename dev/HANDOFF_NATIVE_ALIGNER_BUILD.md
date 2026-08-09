# HANDOFF — build the native −logP junction re-aligner (arm C)

**Branch:** `worktree-agent-a25a2c1e784ad37dc`. **Written:** 2026-07-02 (polished, durable).
**Supersedes** the working notes in `dev/HUMAN_SIM_RESCOPE_AND_FINDINGS.md` (kept for the
blow-by-blow; this doc is the clean pickup). Report to the PI (Kevin); commit surgically;
NEVER commit to `drs-validation-rebuild`.

---

## TL;DR — read this screen before anything

**What the member is.** A **motif-blind, −logP-scored junction RE-ALIGNER over the panel's
localization window.** The panel (minimap2/uLTRA/deSALT/…) localizes each read correctly but can
mis-call the junction *inside* that window; the member re-places the junction there, scoring
candidate placements by calibrated ONT error probability (−logP) instead of flat penalties, with
the GT-AG/annotation preference removed. It is a downstream re-aligner, **not** a 6th standalone
placer.

**The decision.** BUILD it (PI decided, twice). Rationale = "build to understand" + mapPacBio was
never updated for ONT. The build is cheap: it is arm B + a cost-table swap in an existing DP engine.

**THE MAKE-OR-BREAK — and it splits into two halves with very different testability** (advisor,
folding COMPASS): the member's job at a snap-vs-hold decision has a **do-no-harm half** (make the
*evidence-based* choice, don't inherit the panel's) and an **addressability half** (recover a flatten
the panel MISSES).
- **Do-no-harm IS real-data-testable NOW** — the 3 cross-platform anchors, snap-or-hold (Reconciliation
  §). Run it FIRST, before any sim; it's the first thing the built member should do.
- **Addressability (does arm C beat arm B) is the hard, UNPROVEN half** — currently un-validatable by
  any instrument we have:
- Real flattened junctions are *unobservable by construction* (the whole thesis — an aligner that
  snaps a true non-canonical to a nearby canonical hides it). So addressability needs *constructed*
  truth.
- The constructed sim can't test it either: `error_injector.py` is **not homopolymer-context-aware**,
  so it under-represents exactly the ONT error signature −logP exploits.
- **If you make the sim realistic using RECTIFY's own `error_rates`/`penalty_scores` tables — the same
  tables the −logP scorer uses — C wins BY CONSTRUCTION** (the 0.09→1.07 hill-climb, one level up). Do
  not do this.

**Two honest ways forward (name which one you take):**
1. **Make the injector HP-context-aware from an INDEPENDENT error model** — SIRV/spike-in-derived or a
   published ONT profile, NOT RECTIFY's tables. Then arm C vs arm B on that sim is a real test. This is
   the key enabling work; it is not free.
2. **Otherwise arm C's addressability is BUILD-TO-UNDERSTAND — explicitly unvalidated.** Legitimate
   given the PI's decision, but it MUST be labeled, never blurred into "validated."

**The honest split (do not conflate — this drifted once already):**
- **PROVEN:** on real A549, **96% of read-genome differences are ONT errors, 4% variants** → treating a
  difference as an error (not a variant to preserve) is correct → justifies using an error model *at all*.
- **UNPROVEN:** that −logP's probability weighting beats flat cost *at junction placement* (the C>B test).
  The 96% number is NOT evidence arm C works. Never quote it as such.

---

## RECONCILIATION — the native aligner's PURPOSE vs "C1" (folds `dev/HANDOFF_FOR_NATIVE_ALIGNER.md`)

The COMPASS/short-read session's handoff and this one describe ONE layered build with DIFFERENT emphasis, and
the difference is load-bearing — the next agent must get this, not deconflict it themselves (advisor + PI).

**The native aligner's DEFINING purpose is the junction RE-ALIGNER, deliberately agnostic to BOTH (i) prior
intron/exon annotation AND (ii) what counts as a "canonical" splice site** (PI). It ENUMERATES and MOVES
junction placements on read evidence — motif term zeroed (canonical-agnostic), `is_novel`/`tier` priors
neutralized (annotation-agnostic). **Moving junctions is the whole point.**

**"C1" (the COMPASS/design-doc framing) is the −logP EMISSION SCORER — a component, NOT that purpose.** It is
the per-position empirical gap cost inside an exon block (`hp_penalty` −logP replacing flat affine). As an
*indel corrector* it is **boundary-invariant** — re-attributes indels within FIXED junction edges, so it
*cannot move a junction* (FDR axis = false-INDELS). That boundary-invariance is a real safety/validatability
property, but it is **at odds with the native aligner's purpose, exactly as the PI flagged.** C1 is the
ENGINE the junction re-aligner calls to score candidate placements; it is **not** the native aligner.

**The build is LAYERED, and the emphasis is the junction-mover:**
- **Inner = the −logP scorer ("C1"):** boundary-safe, validatable on the HP benchmark. Build FIRST *because
  it is the safe, testable piece* — NOT because it is the goal.
- **Outer = the junction re-aligner (the native aligner):** enumerate candidate junctions, zero the motif
  term, neutralize the annotation priors, score each candidate's flanking blocks with the inner scorer, pick
  the best on evidence. **This is what fulfills the defining purpose.**

**CRITICAL — the layering is additive in CODE, NOT in SAFETY (advisor).** The moment the outer loop *chooses*
junction positions it DISCARDS C1's boundary-invariance: it can move/fabricate junctions → FDR axis becomes
false-JUNCTIONS (the harder problem, and the whole reason moving junctions is powerful AND dangerous). **Do
NOT let "C1 passed its indel-FDR test" imply the junction re-aligner is safe** — the outer layer needs its
OWN stricter false-junction control. Same hole as the **C6 variant-adjacent trap**: at a splice-site-proximal
genomic variant (the 4% real variants) a junction-MOVER can re-express the variant as a spurious junction; a
boundary-invariant indel-corrector cannot. **The do-no-harm test MUST include variant-adjacent junctions.**

**OPEN IMPLEMENTATION DECISION — resolve by READING CODE, not hypothesis (first build step).** Which DP is
the inner scorer? The COMPASS redteam warns `align_exon_block_global` (indel-within-block, `local_aligner.py:~522`)
is a DIFFERENT DP from `align_clip_to_exon` (movable junction boundary, the 5' rescue). The source docs CONFLICT
(member-design digest: "repoint align_exon_block_global to junction placement"; COMPASS: that engine can't move
junctions). READ BOTH before writing — picking wrong wastes days.

**−logP scorer wiring (Tier A MVP, from COMPASS §3):** reuse `_precompute_del_costs` (`hp_penalty.py:315`) +
add symmetric `_precompute_ins_costs`; inject per-ref-position cost into the D/I recurrence; add
`penalty_table: Optional[HpPenaltyTable]=None` where **None ⇒ byte-identical flat-affine** (backward-compat
fence). Tier B = run-unit length-transition emission.

**First real-data test = snap-or-hold on the 3 anchors** (corrects this doc's earlier "no real-data
substrate"; the junction layer is HALF real-data-testable). The 3 cross-platform anchors (SLC35A4/TMED9/SQSTM1;
`dev/COMPASS_2corroborated_CROSSPLATFORM.md`) test the DO-NO-HARM half — does the member SNAP to the nearby
canonical or HOLD the non-canonical? **Either is a PASS if defensible from the reads; the ONLY failure is
silently inheriting the panel's choice.** (J2 = a "should-abstain" case.) ADDRESSABILITY (recover a flatten the
panel MISSES) has NO real-data substrate. **Run snap-or-hold FIRST.**

**GMAP drop-test (downstream deliverable):** GMAP is low-unique-value (~1 well-supported unique junction of
201; `dev/GMAP_PANEL_APPROPRIATENESS.md`). Retire it by showing C1 + minimap2 + (deSALT|mapPacBio) recovers
its genuine unique finds (chr1:19219782-class) WITHOUT its noise; until then GMAP stays fenced
(`test_gmap_fence_regression.py` must pass). Mechanics/data: `dev/HANDOFF_FOR_NATIVE_ALIGNER.md §2,§5`.

**Headline caveat (short-read bias):** "recall-gap is 96% artifact" INHERITS the COMPASS gate's requirement of
independent SHORT-READ corroboration → a genuine LONG-READ-ONLY junction is undercounted as artifact. The gap
size is a **ceiling-on-real / floor-on-artifact**, not final; the benchmark can adjudicate the long-read-only
tail without short-read dependence (a capability to exploit).

---

## DONE (established this session, with evidence)

1. **Rescoped the task + corrected two committed-record errors** (both code-verified):
   - **Panel-wiring bug:** `run_panel_on_corpus.py` → `run_multi_aligner()` whose dispatcher
     (`rectify/core/align/multi_aligner.py:3052–3132`) has **no branch for deSALT/uLTRA/gmap** → they
     `RuntimeError("All requested aligners failed")` (line 3144). Every prior panel/blindspot run was
     **minimap2 + mapPacBio only**; the committed "the *panel* flattens 40–100%" is really "*minimap2*
     flattens." `run_desalt`/`run_ultra` are wired only into `align_command._run_one_aligner`.
   - **Real-data target is small/mostly-artifact:** v3 COMPASS adjudication of the 111 GMAP-only novels =
     **107 artifact / 3 real / 1 inconclusive** (`compass_a549/adjudication_111_v3.json`); recall-gap
     998/1070 = pure {GMAP,mapPacBio}. The 3 real anchors (SLC35A4/TMED9/SQSTM1) are placed *non-canonically*
     by minimap2/deSALT/uLTRA (agree with 168/323/2959 short reads) — the workhorses recover them.
2. **Human ground-truth sim** (`human_noncanon_sim.py`, GRCh38_chr5): minimap2/uLTRA **do flatten** real
   non-canonical junctions by ~2 bp (SLC35A4 verified: minimap2 `497M597N4I499M` vs truth `500M593N500M`;
   deSALT/mapPacBio place correctly). Reconciles with real data — the ~2 bp snap is *within* the real-data
   eval's ±3 bp tolerance, so real data can't see it but the exact-truth sim can. **This is a confirmed
   mechanism, not a sized gap.**
3. **Addressability gate** (`addressability_precheck.py`, edit-distance, no engine): of reads
   minimap2/uLTRA flatten AND deSALT fails, under realistic error **motif-blindness (arm B) recovers
   79–90%** (ED_true < ED_snap); the **−logP-specific tie sliver is 4–5%** (grows with error); zero-evidence
   6–15%. → arm B is the validated core; arm C is a small increment.
4. **PI's bet confirmed on real data** (`real_diff_composition.py` on real A549 minimap2 BAM, 3×1 Mb
   expressed windows): **ONT error 96.0% / genomic variant 4.0%.** (Justifies the error-model approach —
   see the "honest split" above for what it does NOT prove.)

## VERIFIED (code-verified facts + gotchas — expensive to rediscover; keep)

- **Wiring gap** (above): the sim panel silently drops deSALT/uLTRA. To run the full panel, use
  `rectify align … --aligners all --junction-aligners uLTRA deSALT --organism homo_sapiens` (the
  `align_command` path that produced the real A549 BAMs) — NOT `run_multi_aligner`/`run_panel_on_corpus`.
  (Fixing `run_multi_aligner`'s dispatch is a real, separate bug-fix worth committing.)
- **BAM `=`-convention:** the A549 minimap2 BAM stores SEQ with SAM `=` ("matches reference"). Any per-base
  read-vs-ref comparison must treat `qb=='='` as a match (else 100% phantom mismatch). Bit `real_diff_composition.py`.
- **Soft-masking:** GRCh38 FASTA is soft-masked (lowercase repeats); the sim uppercases its genome →
  UPPERCASE both sides in any edit-distance/pileup comparison. Bit `addressability_precheck.py` (TMED9 ED=67).
- **Read orientation:** get sim read sequences from `reads.fastq` (forward), NOT the BAM `query_sequence` —
  minimap2 flips non-canonical reads to the minus strand (rev-comp'd SEQ). Bit `addressability_precheck.py`.
- **true_cigar reconstruction:** build a read's true spliced reference by walking `ReadTruth.true_cigar`
  (M/=/X emit ref, N/D advance), not the normalized junction coords (a normalization slide inflates ED).
- **R1 deSALT-substrate lesson:** ALWAYS verify an aligner *placed* the read (right locus, sane CIGAR)
  before reading its recall — a mis-placement reads as 0-recall + high `.fpnc` and *looks* like flattening
  but is a seeding failure. deSALT is unusable on yeast-scale short-intron constructs (mis-mapped to wrong
  loci); it partially mis-seeds even at 1 kb human reads → contaminates the "deSALT-failed" candidate set.
- **`error_injector.py` is NOT homopolymer-context-aware** (fires indels at uniform hazard, per its own
  docs) → it cannot fairly test HP-correction / −logP. Magnitudes are PLACEHOLDER-PENDING-SIRV, not RNA004.
- **The tie-probe's `−logP` weights (0.5/3.0/4.0) are REPRESENTATIVE, not the empirical tables** — the sign
  of a reorder can flip with real magnitudes, so nothing from `tie_logp_probe.py` transfers quantitatively.
  The real build MUST use the actual `hp_penalty`/`error_rates` tables.
- **Model:** `model:fable` subagents silently fall back to **Opus 4.8** (Fable content safeguard on benign
  biology); JSONL `model` field is authoritative; use **Sonnet 4.6** for genuine model diversity.
- **SSH:** serialize (Sherlock `MaxSessions` throttle self-heals; `Permission denied (gssapi-with-mic)` may
  mean the ControlMaster dropped → **ask Kevin to re-open**, never re-open it yourself = 2FA cost). rsync
  chokes on the login banner → `tar czf - … | ssh 'tar xzf -'` or `ssh 'cat > dest' < file`.

## OPEN (the make-or-break + the two other traps)

- **C>B is unvalidated** and needs one of the two honest paths (TL;DR). Decide + label which.
- **The 4% variant sites → the C6 failure mode.** Arm C treats *every* difference as a correctable ONT
  error, so at a splice-site-proximal genomic variant it can re-express the variant as a spurious junction.
  **The do-no-harm test MUST include variant-adjacent junctions**, or the member ships a NEW error mode. A
  variant term is the eventual fix (deferred).
- **Gap SIZE is an upper bound** until a **multi-exon realistic sim** (short ~120 nt exons, ~2 kb reads,
  many exons) de-confounds the deSALT-seeding contamination in the candidate set. (Human internal exons
  median ~120 nt; reads span ~8 exons; introns median ~1.5 kb.) `human_noncanon_sim.py` currently builds
  single-junction ~1 kb reads — extend to multi-exon for the definitive run.

## RESUME — start-here for the build (runs without this transcript)

**Build (see RECONCILIATION § for the layering + the OPEN DP-choice; one engine, pluggable cost; A/B/C are
cost CONFIGURATIONS, built together — C = B + table swap):**
- **Inner −logP scorer ("C1"):** `hp_penalty` −logP costs in a boundary-safe block DP. **FIRST decide the DP**
  (Reconciliation §: `align_exon_block_global` `local_aligner.py:~544` indel-within-block vs `align_clip_to_exon`
  movable-boundary — READ BOTH). Tier A wiring: `_precompute_del_costs`+ symmetric `_precompute_ins_costs`,
  `penalty_table=None ⇒ byte-identical flat-affine` fence.
- **Outer junction re-aligner (the native aligner):** enumerate candidate junction placements in the panel
  window (the aligners' placement cluster ± local search), zero the motif term, neutralize the annotation
  priors, score each candidate's flanking blocks with the inner scorer, pick the best. This layer MOVES
  junctions → owns a stricter **false-junction** FDR (not C1's false-indel one).
- **Arm A** = flat affine + GT-AG/motif preference (incumbent). **Arm B** = flat affine, motif term = 0.
  **Arm C** = arm B with the empirical `HpPenaltyTable` (`rectify/core/splice/hp_penalty.py`) +
  `error_rates` costs (the real `penalty_scores*.tsv`, NOT representative weights).
- Neutralize the canonical priors (member config): `_CANONICAL_HP_PRIOR=0.5`
  (`junction_scoring.py:293`, consumed `junction_refiner.py:647`), the `is_novel` tiebreaker
  (`junction_refiner.py:638`, sort priority 4 at 660/663), AND the residual `tier` tiebreaker (priority 2/3,
  same tuple — the adversarial review flagged an ablation of only the first two as incomplete).
- BAM-ize each arm via `scorer.cigar_records_to_bam` (`rectify/core/benchmark/scorer.py:~266`); adjudicate
  with `score_bam` against constructed truth. **The metric is C-minus-B on the SAME substrate** — that delta
  is the −logP contribution, in isolation.

**Validation substrate — decomposed (NO single substrate covers arm C; state this plainly):**
- **SIRV/spike-in** → emission-model calibration + canonical do-no-harm (real ONT errors on known sequence;
  the instrument that caught the C1 bug). **First verify SIRV actually contains non-canonical junctions**
  before planning any flattening-recovery test on it (do not assume the Lexogen sets do).
- **3 anchors (SLC35A4/TMED9/SQSTM1)** → do-no-harm / specificity ONLY (workhorses already recover them;
  arm C must not make them worse). NOT an addressability substrate.
- **A549 vs COMPASS** → confirm COMPASS can actually emit corroborated non-canonical junctions first
  (it was mid-build as of `dev/COMPASS_HANDOFF.md`).
- **Addressability (C>B) has NO real-data substrate.** Either build path 1 (independent-error sim) or label
  arm C build-to-understand (TL;DR). Do NOT imply "validate C>B on SIRV/real" closes it.

**First concrete command (re-run the gate that sized arm B/C; template for the build harness):**
```
ssh sherlock 'source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh; conda activate rectify;
  export PYTHONPATH=/home/groups/larsms/users/kevinroy/nj_panel_code;
  C=/scratch/users/kevinroy/nj_panel/human_sim/corpus_err;
  python /home/groups/larsms/users/kevinroy/nj_panel_code/scripts/benchmark/addressability_precheck.py --corpus $C \
    --bam minimap2=$(ls $C/panel/*minimap2*.bam|grep -v namesort|head -1) \
    --bam uLTRA=$(ls $C/panel/*uLTRA*.bam|grep -v namesort|head -1) \
    --bam deSALT=$(ls $C/panel/*deSALT*.bam|grep -v namesort|head -1)'
```

## FILES

- **Scripts (repo `scripts/benchmark/`, staged on cluster `nj_panel_code/scripts/benchmark/`):**
  `human_noncanon_sim.py` (human sim generator), `addressability_precheck.py` (the arm-B/C gate),
  `tie_logp_probe.py` (tie −logP preview — representative weights, non-transferable), `real_diff_composition.py`
  (real-data ONT-error-vs-variant test), plus existing `panel_blindspot_score.py`, `run_panel_on_corpus.py`,
  `novel_junction_realgenome.py`.
- **Docs:** this file (pickup); `dev/HUMAN_SIM_RESCOPE_AND_FINDINGS.md` (full blow-by-blow);
  `dev/NATIVE_ALIGNER_OVERVIEW.md`, `dev/ALIGNER_MEMBER_DESIGN.md`, `dev/DIRECTOR_ALGO_EVAL_*`,
  `dev/REVIEW_SYNTHESIS_*` (member design consensus); `dev/HANDOFF_ALIGNER_BENCHMARK.md` (running log).
- **COMPASS session (folded above; consult for depth, don't re-transcribe):**
  `dev/HANDOFF_FOR_NATIVE_ALIGNER.md` (its briefing — 3-anchor acceptance test, GMAP drop-test, C1 wiring),
  `dev/GMAP_PANEL_APPROPRIATENESS.md` (GMAP low-unique-value + the drop-test), `dev/COMPASS_2corroborated_CROSSPLATFORM.md`
  (the 3 anchors' DNA-confirmed derivation), `rectify/data/validation/CASE_STUDIES.md` (acceptance criteria).
  Its data: `$W=/scratch/users/kevinroy/compass_a549`, `$OAK=/oak/…/compass_a549_out` (adjudication_111_v3.json,
  reval_chr{5,1,11,17,19}.json, confirm_13.json, wgs/a549_wgs_deep.bam); morechrom LR panels under
  `…/sgnex_a549/morechrom/alignments/`.
- **Cluster:** corpora `/scratch/users/kevinroy/nj_panel/human_sim/corpus_{ef,err,hot}/` (reports:
  `panel_report.txt`, `precheck.txt`, `tie_logp.txt`); `real_diff_composition.txt`; real A549 BAMs (read-only)
  `/scratch/users/kevinroy/rectify_human_validation/sgnex_a549/alignments/`; refs
  `/scratch/users/kevinroy/sumner_lab/references/GRCh38_chr5.fa` (+ `desalt_index/`, `gencode.v44.basic.chr5.gtf`);
  COMPASS `/scratch/users/kevinroy/compass_a549/`; code `/home/groups/larsms/users/kevinroy/nj_panel_code`.
- **M1:** base env `/Users/kevinroy/miniconda3/bin/python` for tests (pysam env lacks pandas/Bio).
