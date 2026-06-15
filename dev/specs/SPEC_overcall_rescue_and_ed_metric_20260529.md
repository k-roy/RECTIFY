# SPEC — 3′-end over-call rescue + raw-ED metric weakness (mapPacBio)

**Status:** DRAFT / for a dedicated session (do NOT interleave with the upf1Δ
validation-set expansion — this touches `rectify/core/correct/` production code
and must land on its own commit boundary, clear of the concurrent human track).
**Date:** 2026-05-29
**Author trigger:** review of `cat1_minus_1` (mpb 5′ junk extension) and
`cat1_plus_1` (terminal non-A 3′ base) validation-read PNGs while evaluating
dropping mapPacBio from the yeast panel.

---

## 1. Finding: raw-ED under-penalizes forced-mismatch-over-softclip (mapPacBio)

### Observation
`cat1_minus_1.png` shows mapPacBio with a visually junk alignment — a 9-bp
forced mismatch run at the 5′ where the other aligners soft-clip (`19S`), plus
scattered forced deletions/mismatches across the body. Yet the displayed **raw
ED is only +4** vs the minimap2/deSALT/gapmm2/uLTRA winner group (174 vs 170).

### Reconciliation (computed, `=`-decoded, ref-position thirds)
| aligner | total mm | 5′ soft-clip | raw ED |
|---|---|---|---|
| minimap2 | 62 (5p13 / mid19 / 3p30) | **19 S** | 167 |
| mapPacBio | **91** (5p34 / mid19 / 3p38) | 0 S | 172 |

mapPacBio has **~29 more mismatches** (the "20–30 more in the 3′ region"
observation is correct, and it's ~29 read-wide). But total ED is only +5
because **`_cigar_raw_edit_distance` charges a soft-clip and a forced mismatch
the same flat 1.0/base** (see `scripts/validation_data/render_read_alignment.py`
`_cigar_raw_edit_distance`: `X`, `I`, `D`, `S`, `H` all cost 1/base; `=`, `N`
cost 0). minimap2 "spends" 19 edits on an honest `19S` ("I don't know"), while
mapPacBio force-aligns those bases as ~29 mismatches ("confident wrong call") —
near-equal cost on a flat metric.

### Interpretation
- The renderer is **faithful** — `walk_cigar` decode == pysam
  `get_aligned_pairs(with_seq=True)` byte-for-byte over the 5′ region. The
  "makes no sense" picture reflects mpb's real alignment, not a parse bug.
- mapPacBio's alignment **is** junk (it manufactures 29 mismatches to avoid
  soft-clipping), and raw ED **fails to penalize it**. Two consequences:
  1. **Argument to drop mapPacBio** from the yeast panel (corroborates its
     prior removal from the human-DRS panel for HP-ED gaming, memory
     `feedback_hped_overcall_gaming`).
  2. **The ED metric should distinguish honest soft-clips from forced
     mismatches.** A soft-clip is a principled abstention; a forced mismatch is
     a wrong base call. Charging them equally lets an over-aggressive aligner
     look competitive. Candidate fix: soft-clip penalty < mismatch penalty
     (e.g. S=0.5, X=1.0), or an explicit "abstention is cheaper than wrong"
     asymmetry. **Caveat:** HP-ED already tried to address aligner gaming and
     mpb games it too via free N-ops — so model the gaming surface before
     changing penalties.

### Action items
- [ ] Decide soft-clip vs mismatch penalty asymmetry for the ED/HP-ED metric;
      re-score the 5-aligner panel; confirm mpb no longer wins on force-aligned
      reads.
- [ ] Quantify how many panel reads mpb "wins" via forced-mismatch-over-softclip
      (the cat1_minus_1 pattern) before/after the penalty change.

---

## 2. Feature request: homopolymer over-call rescue at a non-A 3′ end

### Goal
When rescuing a soft-clipped 3′ end, also consider a **homopolymer over-call**:
the read's true 3′ terminal base is a **non-A** (e.g. a genomic `T`) that should
be MATCHED, even though the basecaller over-called the upstream poly-A run and
pushed extra A's between the body and that terminal base.

### Example (`cat1_plus_1.png`, chrXIV, strand +)
orig 3′ 10617 → corr 10611. The genomic context at the 3′ junction is an A-run
followed by a terminal non-A base (a `T`). Today the rescue/walkback resolves
into the A-run; the request is to recognize the over-called A's and **match the
terminal `T`** as the true 3′ end. This is the mirror image of the existing
`softclip_rescue` (which handles UNDER-called homopolymers,
`rescue_softclip_at_homopolymer`) and of `rescue_overcall_terminal_match`
(`rectify/core/correct/indel_corrector.py:702`, emits `oc_terminal_base` /
`oc_overcall_count`, surgery in `read_edits.extend_read_3prime_for_overcall_rescue`).
Scope: extend the over-call path to fire when the terminal base to match is
**non-A** (currently the over-call logic appears A-centric — verify the gate).

### Empirical grounding — REQUIRED before building, and CURRENTLY UNMET
The feature should only be built if the empirical error model shows Nanopore
**over-calls** homopolymers across the board. Check of the bundled yeast table
`rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv`:

| op | what | rows | counts | rate |
|---|---|---|---|---|
| `D` (under-call) | HP deletion | 33 | **millions** (8.2M at AT hp=1) | 0.73 → 0.99 |
| `I` (over-call) | HP insertion | **5** | **1–3** | ~0.000 |
| `X` | mismatch | 32 | — | — |

**The table shows under-calls (deletions) dominate; over-calls (insertions) are
essentially absent.** So the over-call premise is **not supported by the current
error table.** Two possibilities, both must be resolved first:
1. Nanopore genuinely doesn't over-call HPs much on this chemistry → the
   feature is low-value; the `cat1_plus_1` case may be better explained by the
   poly-A tail boundary, not a generic over-call.
2. The profiler **under-captures insertions** (5 I-rows with single-digit counts
   is implausibly sparse next to millions of D events) → the table must be
   rebuilt with proper insertion/over-call profiling before the feature is
   justified.

### Action items
- [ ] Audit the error-profiler's insertion handling — why only 5 `I` rows? Is
      HP over-call profiling implemented at all? (start: the profiler that emits
      `penalty_scores.tsv`).
- [ ] If over-calls are real but unprofiled: rebuild the table with I-rows, then
      quantify HP over-call rate vs length. Only then design the rescue.
- [ ] Verify the `rescue_overcall_terminal_match` gate's base-identity
      assumption; extend to non-A terminal if grounded.
- [ ] Add a validation read for the non-A-terminal over-call case (the
      `cat1_plus_1` exemplar) once the behavior is defined.

---

## Notes
- Keep this OFF the upf1Δ validation-data commits. Separate branch/commit.
- Both items are correctness-sensitive → same discipline as the rescue DP cap
  (pre/post equivalence diff on real data, focused tests, advisor review).

---

## INVESTIGATION FINDINGS — 2026-05-29 (read-only audit, no code changed)

A grounding pass against the codebase at HEAD `726ebe1` corrected several of the
spec's premises. Net: **most of the assumed implementation work is already done
or aimed at the wrong layer.** What remains is one genuine arbitration/biology
question (Item 2) and an optional cosmetic display fix (Item 1). Evidence below.

### Item 1 — the gaming mechanism is already neutralized in PRODUCTION; flat S=X is a DISPLAY-only artifact

- `_cigar_raw_edit_distance` (`scripts/validation_data/render_read_alignment.py:908`)
  is **display-only** — it feeds the "raw ED" review columns (`PLOTTING.md:407`).
  It does charge `S = X = D = I = 1/base` and `= , N = 0`, exactly as the spec
  says. Editing it changes **no winner** — it is cosmetic review signal.
- The **production winner selection** is `score_alignment`
  (`rectify/core/consensus/scoring.py:665`), called by `select_best_alignment`
  (`rectify/core/consensus/select.py:62`). It already defeats the cat1_minus_1
  gaming pattern by construction:
  - 5′ penalty uses `effective_five_prime_clip = max(explicit_clip,
    terminal_error_length)` (scoring.py:709) — so mapPacBio's *forced mismatches
    in lieu of a soft-clip* are penalized **equivalently** to an honest soft-clip
    (`−2/base` either way). The "force-align instead of clip to look cheaper"
    move does not work on the production scorer.
  - Body mismatch count and 3′ raw-endpoint depth are **deliberately NOT scored**
    (scoring.py:684–732, with documented rationale: aligners maximise their own
    score; `find_polya_boundary()` assigns the CPA regardless of where the raw
    alignment stopped; penalising deeper 3′ extension would disfavour the *more*
    informative aligner). The spec's "~29 more mismatches go unpenalized" is true
    of the display ED but is an intentional non-criterion in production.
- **Consequence:** the spec's proposed fix (`S < X` asymmetry to stop mpb
  "winning") is (a) unnecessary for production selection — already handled via
  `effective_*_clip`; (b) in tension with the documented production design if
  pushed into `score_alignment`; (c) only meaningful as a *review-readability*
  change in the renderer. Lowering the soft-clip penalty would also shift
  intron-finder rankings (memory `feedback_hp_edit_distance_semantics`: the
  per-base soft-clip cost is precisely what biases against intron-finders), so a
  production-side change is a gated, cluster-re-scored, advisor-reviewed effort —
  not a quick edit, and not justified by the cat1_minus_1 observation alone.

### Item 2 — the non-A over-call rescue ALREADY EXISTS, is wired, fires at HEAD, and `cat1_plus_1` is ALREADY CORRECT. The spec's "resolves into the A-run" is a coordinate misread.

- `rescue_overcall_terminal_match` (`rectify/core/correct/indel_corrector.py:702`)
  is exactly the requested feature, and it is **already non-A-centric**: it
  returns `None` when the terminal base IS the stop-base (`indel_corrector.py:798`),
  i.e. it fires *only* for a non-A (+) / non-T (−) terminal. **cat1_plus_1 is its
  literal docstring example** (`indel_corrector.py:712`). The spec's "verify the
  gate / appears A-centric" is resolved: the gate already requires non-A.
- It is **fully wired**: called at `rectify/core/bam/bam_processor.py:765`, emits
  `oc_terminal_base` / `oc_overcall_count`, surgery applied via
  `extend_read_3prime_for_overcall_rescue` at `bam_writer.py:321/505/668/708`.
  Landed in commit `a1728eb` (2026-05-18). It is sequence/geometry-based and does
  NOT read the penalty table at runtime — so the "empirical grounding UNMET"
  concern does not gate whether the feature *works*, only whether it is *justified*.
- **The arbitration is ALREADY overcall-first.** Module 2G.5 (overcall rescue,
  `bam_processor.py:758`) runs **before** Module 2E (poly-A walkback,
  `bam_processor.py:802`); walkback is explicitly gated on
  `not overcall_rescue_applied` (`bam_processor.py:807`). There is **no
  arbitration bug** and no ordering change to make.
- **The coordinate, stated exactly** (0-based, genome `chrXIV`): the A-run is
  **10600–10610**; **`10611 = T`** (non-A); then `A G C T C…`. `corrected_3prime
  = 10611` is therefore the **terminal T**, i.e. the spec's *desired* result —
  NOT "into the A-run." (The spec author read 10611 as an A-run position; it is
  the T one base past the run.)
- **HEAD behaviour, reproduced 2026-05-29** (`rectify correct` re-run on each
  bundled per-aligner BAM at HEAD `726ebe1`, all 5 aligners → `corrected_3prime`
  on the terminal T):

  | aligner   | orig 3′ | corr 3′ | correction_applied            |
  |-----------|---------|---------|-------------------------------|
  | minimap2  | 10610   | **10611** | `atract_ambiguity,overcall_rescue` |
  | gapmm2    | 10610   | **10611** | `atract_ambiguity,overcall_rescue` |
  | uLTRA     | 10610   | **10611** | `atract_ambiguity,overcall_rescue` |
  | deSALT    | 8703    | **8704**  | `atract_ambiguity,overcall_rescue` *(chrVI map; 8704 = the T)* |
  | mapPacBio | 10617   | **10611** | `polya_walkback`              |

  **`overcall_rescue` fires at HEAD for the soft-clipping aligners** and lands on
  the terminal T. (minimap2/gapmm2/uLTRA map this read to chrXIV:10611=T; deSALT
  maps it cross-chrom to chrVI:8704, the homologous T — also via overcall_rescue;
  so the chrXIV-mapping aligners agree on 10611=T.) mapPacBio force-aligns through
  to 10617 (no 3′ soft-clip → overcall GATE1 not met) and reaches the **same** T
  via walkback
  (which stops at 10611=T, the first non-A read=ref agreement past the run — it
  does NOT enter the A-run). The bundled `corrected_3ends.tsv` (consensus winner,
  commit `e39089e`) likewise shows `orig=corr=10611`. Every path lands on the T.
- **Conclusion: `cat1_plus_1` needs no production change.** The feature works.
  The stale review markdown (`cat1_indel_review.md`, showing the mapPacBio winner
  via `polya_walkback`, "10617→10611") is accurate *for mapPacBio* but was
  misinterpreted: 10611 is the terminal T, the correct answer.
- **Biology (user-confirmed 2026-05-29): the terminal T IS the true 3′ end** for
  this geometry — which is exactly what HEAD already produces. No regression.

### Item 2 sub-finding — the "5 I-rows" scarcity is a PROFILER-VERSION artifact, not evidence over-calls are absent (and NOT the isolation filter)

- The bundled DRS `penalty_scores.tsv` has 5 trivial I-rows (counts 1–3) and
  README:139 records "no insertion events at HP=1 anywhere in 6.99M × 2 reads";
  insertions consequently fall back to a flat `--default-ins 1.25` baseline
  (README:143–144).
- **Initial hypothesis (the `--isolation-flank 10` filter hid them) was TESTED
  and REFUTED.** A controlled re-run on the same upf1Δ 3-aligner data at flank=10
  vs flank=0 (only that variable changed) shows flank=10 still yields **56 I-rows
  / 1.23M I-events** (vs 60 / 2.13M at flank=0) — the filter reduces insertions
  ~42% but does NOT hide them. Strict (≥2-of-3) ≈ union at both flanks. (Caveat:
  the control panel is minimap2 + gapmm2 (mm2-family) + deSALT at `min_aligners=2`;
  the *bundled* table is **5-aligner STRICT** — a config not separately tested, and
  5-way exact-position agreement could itself suppress I-rows.)
- **Cause — profiler version (strongly indicated by git, not definitively proven).**
  The insertion path itself
  (`empirical_cigar_error_profiler.py:1290` "Phase 5" + the `any_ins_before`
  union counting) was added in commit `39458d3` (2026-05-17 **13:43**), **7 minutes
  AFTER** `penalty_scores.tsv` was committed (`ba602c3`, 13:36), and the DRS table
  was never regenerated afterward. So the bundled table was produced by a profiler
  build with no working insertion counting — hence the near-zero I-rows and the
  README "no insertion events" note. The current profiler, at the *same* documented
  flags (flank=10, union), finds 55,673 `I A` hp1 events in 500k reads where the
  bundling run reported zero. The 50,000× gap points to a code-version difference
  rather than data or the isolation filter (both of which the control rules out);
  the one alternative not separately excluded is 5-aligner strict agreement.
  Caveat: the README documents a `--union` flag that didn't exist at the 13:36
  commit, so it isn't a reliable record of how the bundled table was actually made.
- **Consequence:** the bundled production penalty tables simply lack empirical
  insertion penalties (built before insertion profiling existed). Re-profiling
  with the current profiler — at *any* isolation_flank — populates them. To ground
  Item 2 (done above): re-profile insertions with the current profiler (e.g. an
  HP-aware insertion tally with `isolation_flank=0`, or a dedicated 3′-boundary
  over-call survey), then measure over-call rate vs HP length. The spec's
  possibility #2 ("profiler under-captures insertions") is the correct branch.

### Revised action items (supersede §1 and §2 above)

**Item 1 (display ED) — safe, in-scope, optional:**
- [x] **DONE 2026-05-30.** Renderer's `_cigar_raw_edit_distance` now charges
      `S/H = _SOFTCLIP_ED_COST (0.5)` per base vs `X = 1.0` — honest abstention <
      forced wrong call. Display-only (the two secondary annotation sites bumped
      `.0f`→`.1f` to show the half-integers). Verified: `19S`→9.5 vs `19X`→19.0;
      renderer parses; no test asserts this metric so nothing broke.
- [x] Confirmed NOT changing `score_alignment` — the 5′ gaming path is already
      equalized via `effective_five_prime_clip`; body-mismatch non-scoring is
      intentional. Any production-scorer change remains a separate, gated effort.

**Item 2 (over-call rescue) — feature works at HEAD; only the GROUNDING survives:**
- [x] Biology decided (user, 2026-05-29): terminal-T is the true 3′ end — which is
      what HEAD already produces for `cat1_plus_1` (all 5 aligners → 10611 = T).
- [x] ~~Implement walkback ↔ overcall arbitration~~ **MOOT** — arbitration is
      already overcall-first (`bam_processor.py:758` before `:802`) and overcall
      rescue already fires at HEAD. No production change.
- [x] **Re-profiled HP insertions with `isolation_flank=0`** (Sherlock job
      `26829646`, upf1Δ rep1 DRS, minimap2+gapmm2+deSALT, 500k reads,
      `--isolation-flank 0 --union --str-repeat`; output
      `$SCRATCH/oc_grounding_isoflank0_20260529/`).
      **Result — HP over-calls are abundant and HP-length-dependent at flank=0:**
      - I-rows: **60** (vs 5 in the bundled flank=10 table), ~2.12M I-events.
      - Over-call rate (`I AT`, base-class aggregated) rises with HP length:
        hp1 0.26% → hp3 1.08% → **hp6 1.35% (peak)** → declines to ~0.4–0.7% by
        hp10–12. Per-base curves match (`I A` hp1 0.24%→hp6 1.36%; `I T` similar).
        Support is large (`I A` hp2 = 261k, hp3 = 191k; `I T` hp2 = 275k).
      - **Not single-aligner noise:** strict (≥2-of-3 aligners agree on the
        insertion position) gives **2,118,184** I-events vs union's **2,125,903**
        — 0.4% apart, both 60 rows. The over-calls are multi-aligner-concordant,
        which argues against the "flank=0 admits messy-alignment artifacts" reading.
      - Magnitude vs under-call: deletions still dominate and rise monotonically
        (D `A`: hp1 1.1% → hp5 3.0% → hp10 10.4%); HP over-calls run ~⅓–½ the
        deletion rate at short-to-mid HPs.
      - **Causal attribution — control RAN (job `26834583`, COMPLETED).** Same
        upf1Δ data at `--isolation-flank 10` (only that var changed): **56 I-rows /
        1.23M events** (strict and union both). So the filter is **NOT** why the
        bundled table is empty — it only trims insertions ~42%. The real cause is
        **profiler version**: insertion counting (commit `39458d3`) was added 7 min
        *after* the bundled table was committed (`ba602c3`) and the DRS table was
        never regenerated. Over-call abundance is robust to flank choice and to
        strict-vs-union → the feature is solidly justified.
- [x] **ALREADY COMMITTED** — `tests/test_validation_reads.py::test_3prime_exact_position`
      already asserts `cat1_plus_1 → corrected_3prime == 10611` ("via overcall_rescue,
      terminal T past genomic A-run"). Re-verified passing 2026-05-30 (8/8 cat1 tests).
      No new test needed; the working behaviour is already locked.

**Downstream finding (separate track, NOT in this commit) — production penalty
tables under-represent insertions:**
- The bundled DRS `penalty_scores.tsv` was built by a profiler version that
  predated working insertion counting (see Item 2 sub-finding: insertion code
  `39458d3` landed 7 min after the table commit `ba602c3`; never regenerated). So
  its HP-insertion penalties are absent and fall back to the flat `--default-ins
  1.25` baseline. The current profiler (grounding run, any isolation_flank) shows
  the empirical insertion penalty varies strongly with HP length (e.g. union
  penalty_score: hp1 1.25 → hp6 0.24 → hp15 1.65). Regenerating the production
  tables with the current profiler would make HP-ED scoring insertion-accurate —
  but this **shifts production HP-ED rankings**, so it is a gated change (same
  discipline as Item 1: pre/post on real data + advisor review), to be raised with
  the penalty-table maintainers, not bundled here.
- **CHECKED 2026-05-30:** the cDNA (`penalty_scores_cdna.tsv`, 33 I-rows) and qsrev
  (`penalty_scores_qsrev.tsv`, 25 I-rows) tables — regenerated 05-20 in `ed331e7`,
  after the insertion code — **already carry proper I-rows**. Only the **DRS**
  `penalty_scores.tsv` (5 trivial I-rows) is stale. So the fix is narrow: regenerate
  the DRS table alone.
- **NOT done in this session (deliberately).** Doing it right needs (a) the original
  production provenance — WT by4742 full-read **5-aligner** panel, not the upf1Δ
  3-aligner 500k subset used for grounding (that subset is wrong genotype + wrong
  panel for a production WT table); (b) those input BAMs re-staged (the
  `wt_by4742_rep1_chunked_20260412` dev-run was cleaned); and (c) a pre/post HP-ED
  winner-selection diff on real DRS data + advisor review before the bundled table
  is swapped. This is its own task for the penalty-table maintainers, not a
  wind-down item.
