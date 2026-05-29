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

  **`overcall_rescue` fires at HEAD for all four soft-clipping aligners** and
  lands on the terminal T. mapPacBio force-aligns through to 10617 (no 3′
  soft-clip → overcall GATE1 not met) and reaches the **same** T via walkback
  (which stops at 10611=T, the first non-A read=ref agreement past the run — it
  does NOT enter the A-run). The bundled `corrected_3ends.tsv` (consensus winner,
  commit `e39089e`) likewise shows `orig=corr=10611`. Every path lands on the T.
- **Conclusion: `cat1_plus_1` needs no production change.** The feature works.
  The stale review markdown (`cat1_indel_review.md`, showing the mapPacBio winner
  via `polya_walkback`, "10617→10611") is accurate *for mapPacBio* but was
  misinterpreted: 10611 is the terminal T, the correct answer.
- **Biology (user-confirmed 2026-05-29): the terminal T IS the true 3′ end** for
  this geometry — which is exactly what HEAD already produces. No regression.

### Item 2 sub-finding — the "5 I-rows" scarcity is a PROFILER METHODOLOGY artifact, not evidence over-calls are absent

- The bundled DRS `penalty_scores.tsv` was generated with
  `--isolation-flank 10 --union` (README:156–160). The profiler's insertion path
  (`empirical_cigar_error_profiler.py:1290` "Phase 5") requires `isolation_flank`
  clean flanking matches on **both** sides of the inserted position
  (`:1313–1318`). **HP over-calls occur precisely where there is no clean 10 bp
  flank** (inside/adjacent to the homopolymer run), so they are structurally
  excluded — even on the permissive `--union` path. README:139 even records "no
  insertion events at HP=1 anywhere in 6.99M × 2 reads," and insertions are
  consequently handled by a flat `--default-ins 1.25` baseline (README:143–144).
- **Consequence:** the current error table **cannot** confirm or deny the HP
  over-call premise — it is blind to it by design. To ground Item 2 you must
  re-profile insertions WITHOUT the isolation requirement in HP context (e.g. an
  HP-aware insertion tally with `isolation_flank=0`, or a dedicated 3′-boundary
  over-call survey), then measure over-call rate vs HP length. The spec's
  possibility #2 ("profiler under-captures insertions") is the correct branch.

### Revised action items (supersede §1 and §2 above)

**Item 1 (display ED) — safe, in-scope, optional:**
- [ ] (cosmetic) Give the *renderer's* `_cigar_raw_edit_distance` an S<X (e.g.
      S=0.5) asymmetry so reviewers see honest-abstention vs forced-wrong-call
      reflected in the review columns. Does NOT touch production selection.
- [ ] Do NOT change `score_alignment` for this reason — the 5′ gaming path is
      already equalized via `effective_five_prime_clip`; body-mismatch non-scoring
      is intentional. Any production-scorer change is a separate, gated effort.

**Item 2 (over-call rescue) — feature works at HEAD; only the GROUNDING survives:**
- [x] Biology decided (user, 2026-05-29): terminal-T is the true 3′ end — which is
      what HEAD already produces for `cat1_plus_1` (all 5 aligners → 10611 = T).
- [x] ~~Implement walkback ↔ overcall arbitration~~ **MOOT** — arbitration is
      already overcall-first (`bam_processor.py:758` before `:802`) and overcall
      rescue already fires at HEAD. No production change.
- [ ] **Re-profile HP insertions with `isolation_flank=0`** (HP-aware) to quantify
      the real over-call rate vs HP length — the one piece of authorized work that
      survives. Independent of `cat1_plus_1`; it answers "how common are HP
      over-calls / is dedicated over-call handling worth keeping" and fixes the
      `isolation_flank=10` blindness documented above. *(In progress.)*
- [ ] (optional) Promote `cat1_plus_1` to a committed regression read that asserts
      `corrected_3prime == terminal-T` so the working behaviour can't silently
      regress (the in-tree review markdown was already stale).
