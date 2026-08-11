# INSCOST REFERENCE-COLUMN AUDIT — refcol model correctness lens

## LENS
Is reference-column ins_cost the RIGHT MODEL? An insertion has NO reference base.
"Genome HP context at the aligned position" is a CHOICE (base before? after?
run being inserted into?). Scrutinize the builder's exact indexing:
- Does it charge the over-call against the correct genome run, or an off-by-one/
  adjacent context?
- Construct the case where ref-column is WRONG: read poly-A insertion adjacent to
  genome NON-homopolymer (ref hp=1 -> expensive) vs inside genome poly-A
  (ref hp=long -> cheap) — which is biologically right for a PORE over-call?
- Does ref-column just RELOCATE gaming rather than remove it?
Certify ref-column sound OR expose where it mis-charges.

## PLAN
1. Read the durable records: INSCOST_REFCOL_BUILD.md (builder), INSCOST_SYNTHESIS.md,
   INSCOST_AUDIT_model-correctness.md, INSCOST_INVESTIGATION.md.
2. Read the actual code: junction_scoring.py _score_hp_anchored (693),
   _score_junction (807), _precompute_del_costs (the deletion mirror the refcol
   claims to copy). Find whether refcol is BUILT here (flag) or only proposed.
3. Determine the EXACT indexing refcol uses for an insertion position:
   which genome base/run does it read the HP length from?
4. Construct the biological-correctness argument for a PORE over-call insertion:
   is the physically-correct charge based on the READ's run, the GENOME's run,
   or neither?
5. Adversarial cases: (a) poly-A insertion at genome non-HP boundary,
   (b) poly-A insertion inside genome poly-A, (c) relocation-of-gaming check.
6. Verdict: sound / mis-charges / relocates gaming.

## CHECKPOINTS
(append one line per sub-step)
- [read] BUILD.md + model-correctness audit read. CRITICAL: builder ran the A/B and
  FALSIFIED the audit hypothesis. percut->refcol on fair no-guard: 358 wins (4->0) BUT
  196 losses (0->4) -- demotion class BALLOONS (196 vs full-run's 6), does NOT vanish.
  12.39% churn (16x full-run). Net recovery +0.011 (vs full-run +0.005) but via massive
  bidirectional churn. Ref-column does NOT cleanly dominate. So the audit's central
  prediction is already empirically wrong. My job: explain WHY via the indexing model.
- [design-confirmed-per-builder] Refcol indexing (BUILD.md lines 30-48): ins cost is a
  PER-GAP quantity, length R+1 vector. gap 0 = right neighbor only (ref[0]); gap j (1..R-1)
  = max(run(gp_{j-1}), run(gp_j)); gap R = left neighbor only (ref[R-1]). hp_base per
  empirical_cigar_error_profiler Phase 5 (1319-1324): hp_left=hp_arr[pos-1], hp_right=hp_arr[pos],
  hp_len=max(hp_left,hp_right). This is the KEY object to scrutinize for my lens.
- [code-read] refcol impl (hp_penalty.py 432-519) FAITHFULLY reproduces the calibrator
  axis: _run_and_base uses genome HP run via _hp_run_length(genome_seq, gp); per-gap
  max(run_l, run_r); hp_base = base of the LONGER-run side. Matches profiler Phase 5
  (1319-1324) EXACTLY. So "indexes the calibration axis" is TRUE -- no off-by-one on the
  axis definition itself. DP wiring (junction_scoring 834-852): above=prev[j]+ins_col[j],
  the insertion at DP COLUMN j (between ref[j-1] and ref[j]) charges ins_col[j].
- [KEY MECHANISM] The DP MINIMIZES over placement. It will place the inserted poly-A run
  at whichever COLUMN j gives the cheapest ins_col[j]. ins_col is max(run_l,run_r) of the
  GENOME. So the DP is free to slide the insertion to the column adjacent to the genome's
  LONGEST HP run within the ref window -> cheapest ins_cost. THIS is the relocation-of-gaming
  candidate: per-cut split the RUN; refcol lets the DP CHOOSE THE COLUMN. Need to test whether
  this re-introduces gaming via column choice.

- [ROOT DEFECT identified w/ advisor] ins_col[j] is BASE-IDENTITY-INDEPENDENT. DP transition
  (junction_scoring:848) `above = prev[j] + ins_col[j]` NEVER consults query[i-1]. Insertion
  cost depends ONLY on genome context, not WHICH base is inserted. Consequence:
    * inserted A inside/adjacent genome A-run -> cheap -> CORRECT (that IS a pore over-call).
    * inserted NON-A adjacent genome A-run -> ALSO charged cheap ins_cost(8,'A') -> WRONG.
      A pore over-call = extra copies of the run's OWN base. A C parked next to an A-run is a
      real insertion/fabrication and must be expensive.
  max(run_l,run_r) compounds: extends cheap zone one column PAST each HP boundary.
  base-agnostic + max + DP-minimization => gaming RELOCATED, not removed. HOW a non-canonical
  alt beats the true canonical -> the 196 demotions. Calibrator measured DESCRIPTIVE marginal
  rate UNCONDITIONED on inserted base; a max-based cost in a MINIMIZING DP inverts that
  descriptive axis into an exploitable discount-finder. Axis-match is NECESSARY not SUFFICIENT.
- PLAN NOW: (1) confirm code fact [DONE by read: no query dep in above=]. (2) witness genome
  AAAAAAAA + non-A tag; read=8A + 4 inserted C + tag; expect refcol charges 4*ins_cost(8,'A')
  ~0.79 NOT 4*ins_cost(1)~5.0 -> mis-charge witness. (3) tie to 196 losses via existing fair
  arms (percut_ng vs refcol_ng bams): cheap insertions base-MATCHED (aggressive-sound) or
  base-MISMATCHED/at max-boundary (mis-charge)?

## WITNESS RESULT (in-process, worktree flag, DRS penalty_scores.tsv) — DECISIVE
ins_cost table (base-class split IS present): ins_cost(8,'A')=0.1972; ins_cost(8,'C')=0.5729;
ins_cost(1,anything)=1.25.
Genome: [C*50][A*8][GCTAGCTAGCTT][G*50]; ref window starts at the 8-A run.
refcol_ins[0:9]=0.1972 (A-run cols), then 1.25 from col 9 (non-HP). Faithful axis.

  over-call 4xA (REAL pore error): refcol_OFF=2.7528  refcol_ON(refcol)=0.7888
  insert   4xC (FABRICATION)     : refcol_OFF=1.8492  refcol_ON(refcol)=0.7888   <-- SAME cheap cost

MIS-CHARGE CONFIRMED. Refcol charges 4 inserted C's the SAME 0.7888 (=4*ins_cost(8,'A'))
as 4 genuine A over-calls, because ins_col[j] never consults the inserted base. Biologically
correct cost for 4 non-A insertions in this context ~ 4*ins_cost(1,'C')=5.0 (or at minimum the
C-run cost, not the A-run cost). The table HAS the base-class info (ins_cost(8,'C')=0.5729 exists)
but refcol looks up the GENOME base (A, the longer-run side via max()), discarding the QUERY base.
NOTE the DIRECTION vs legacy: legacy(OFF) charges the fabrication (1.8492) LESS than the real
over-call (2.7528) -- also imperfect -- but refcol collapses BOTH to the A-run floor, so a
fabricated non-A bridge into a genome A-run is as cheap as a true poly-A pore stutter. This is
exactly the discount a non-canonical alternative needs to beat the true canonical -> mechanism
for the 196 (0->4) demotions the builder measured. RELOCATION of gaming, not removal.

## TIE-TO-196-LOSSES (fair arms percut_ng vs refcol_ng, already built) — IN PROGRESS
Arms confirmed: mix_fair_out/arm_mbF_percut_ng.bam, arm_mbF_refcol_ng.bam.
Plan: pull reads that went tier 0 (canonical, percut) -> tier 4 (non-canonical, refcol) on the
has_true_junction set; for each, examine the refcol-chosen junction and whether the discount that
let the alternative win comes from insertions parked at a genome A/T HP boundary (base-MISMATCHED
= mis-charge/relocation) vs genuine within-HP ambiguity (base-MATCHED = aggressive-but-sound).
Adjudicator _ins_compare.py buckets tier_shift (toff->ton) using _canonical_tier; I extend it to
dump the per-read refcol junction + local genome context for the 0->4 class.

## TIE-TO-196 RESULT (fair percut_ng vs refcol_ng) — MECHANISM CONFIRMED
196 reads went tier 0 (canonical, percut) -> tier 4 (non-canonical, refcol) [matches builder].
Dominant pattern: TRUE acceptor = 270; refcol slides acceptor to 271/272/273/275 -- INTO a
genome 9-A run (A(9,'A') flanking the refcol acceptor). 132/196 (67%) have a poly-A/T run>=4
flanking the refcol junction end; 65/196 sit at a genome HP>=4 EDGE.
The junction is being DRAGGED deeper into a genome poly-A run because ins_col there is the
cheap A-run cost (ins_cost(9,'A') ~0.15-0.20) regardless of what the read bases are. The DP
converts would-be matched/true-acceptor alignment into a non-canonical junction + cheap
in-A-run insertions. THIS is the mechanism of the 196 demotions, and it is precisely the
base-identity-independent + max()-boundary discount from the witness, now firing on real
simulated reads. NEXT: confirm whether the shifted-in read bases are A (matched, arguably a
real slide) or non-A (mismatched, hard mis-charge) to sharpen aggressive-vs-mischarge.

## PER-READ BASE-IDENTITY EVIDENCE (5 demotion reads) — BOTH FLAVORS PRESENT
Genome accep region: ...CCCAGAAAAAAAAACTAAGT  (9-A run at ~262-270; true acceptor 270 = AG, tier 0).
refcol slides acceptor 270 -> 271/272/273/275 (all give acc dinu 'AA'/'GA' = tier 4 non-canonical).
  tid_0_r001: refcol (178,272), ins@272='A'   -> A into A-run (base-MATCHED slide) BUT flips
              canonical AG(270) -> non-canonical AA(272). True call DEMOTED.
  tid_0_r025: refcol (180,273), ins@273='CA'  -> a C parked at ref273 charged at A-run cost
              (base-MISMATCHED = HARD MIS-CHARGE). true 270 demoted.
  tid_0_r323: refcol (180,275), ins@275='AAAA'-> 4 A's parked in A-run (matched) but acceptor
              275 = tier4. true 270 demoted.
  tid_0_r176: refcol (180,271) tier4, no ins right at jn (slide-only into run edge).
  tid_1_r027: refcol (178,272), ins@272='T'   -> a T parked at A-run col charged cheap
              (base-MISMATCHED mis-charge).
CONCLUSION: the 196 demotions include BOTH (a) base-matched poly-A slides that still convert a
TRUE canonical GT/AG call into a fabricated non-canonical AA junction (the max()-boundary +
DP-minimization pulling the acceptor into the run), AND (b) base-MISMATCHED insertions (C/T)
charged the cheap A-run cost (the base-identity-independence hard mis-charge). Either way the
TRUE canonical acceptor is lost. Ref-column does NOT eliminate the demotion class; it CREATES a
far larger one (196 vs full-run's 6) precisely because it made in-run insertions cheap without
conditioning on the inserted base or protecting the canonical dinucleotide.

## FULL 196-READ CLASSIFICATION (fair percut_ng vs refcol_ng) — QUANTIFIED
total tier 0->4 demotions: 196
  base-MISMATCHED insertion parked at genome HP>=4 column, charged cheap HP cost (HARD MIS-CHARGE): 61/196 (31%)
     mismatched inserted bases charged at HP cost: T=46, G=36, C=2  (i.e. non-A/T bases getting the poly-A/T discount)
  only base-MATCHED insertions near jn (poly-A slide, still DEMOTES a true canonical call): 133/196 (68%)
  no ins within +-6 of jn (pure slide): 2
  acceptor slide delta (refcol - true): +1:51, +2:64, +3:11, +5:5, +10:65  (junction dragged INTO the run)
=> BOTH failure modes are real and quantified: 31% are the hard base-identity-independence
   mis-charge (a genome-context discount applied to a base the pore did NOT stutter), and the
   remaining 68% are matched slides that still fabricate a non-canonical junction because refcol
   makes in-run insertions cheap and does not protect the canonical dinucleotide. NET: refcol
   trades run-splitting gaming (per-cut) for column-choice gaming (drag-into-run) + a base-blind
   discount. It RELOCATES the gaming, and enlarges the canonical-demotion class 33x (196 vs 6).

## ADJUDICATOR RECONFIRM (fair no-guard percut->refcol)
CALL-CHANGE 694/5600 = 12.39%. tier shift: 0->4=196, 4->0=358, 0->2=20, 2->0=2, 4->4=118.
Direction: recov_LOST x became_LESS_canonical = 196 (the demotions ARE true recovery losses,
not neutral reshuffles). recov_GAINED x MORE_canonical = 255. Matches builder exactly.
=> The 196 losses are genuine: a truly-canonical read's correct junction is LOST. R0flank
   recovery net +0.011 but bought with 196 real canonical losses -- the audit's predicted
   "0->4 losses vanish" is FALSE by 33x.

## TASK-CASE WITNESS (poly-A adjacent NON-HP vs inside poly-A) — CONCEDE WHERE REFCOL IS RIGHT
(a) 4 A's inserted adjacent to a genome NON-HP column: refcol = 4.5832 (~ 4*ins_cost(1,'A')=5.0)
    -> EXPENSIVE. Biologically CORRECT: extra A's with NO genome A-run to explain them are not a
    plausible pore stutter, so charging them fully is right.
(b) 4 A's inserted inside a genome poly-A(8) run: refcol = 0.7888 (= 4*ins_cost(8,'A'))
    -> CHEAP. Biologically CORRECT: a pore over-call of a run's own base IS cheap.
=> For an ISOLATED pore over-call, refcol's per-column charge is biologically CORRECT in BOTH
   sub-cases. The fault is NOT the isolated charge magnitude.

## VERDICT (refcol-model-correctness lens) — FAULT FOUND; recommendation does NOT survive
fault_found = TRUE. The audit's recommendation ("refcol keeps the ~35 anti-fab wins AND
ELIMINATES the 4-6/panel canonical-demotion losses") is FALSIFIED by the builder's own A/B and
re-confirmed here: refcol produces 196 tier 0->4 demotions on fair no-guard (vs full-run's 6) --
a 33x LARGER canonical-demotion class, not its elimination. These are genuine recovery losses
(recov_LOST x LESS_canonical = 196), i.e. a truly-canonical read's correct GT/AG junction is
lost.

ROOT CAUSE (STRUCTURAL, the headline): refcol charges insertion by GENOME-COLUMN context ALONE.
Inside/adjacent a genome poly-A run the per-column ins_cost is cheap, so placing the acceptor
DEEPER in the run (with cheap in-run insertions) is ALWAYS <= placing it at the run EDGE (the
true canonical acceptor). The junction search has NO term asserting "the canonical dinucleotide
is here" and no gradient penalizing drift off it, so the minimizing DP slides the acceptor into
the poly-A run (slide-delta +1..+10) and fabricates a non-canonical AA junction. This is the
DOMINANT mechanism: 133/196 (68%) of demotions are base-MATCHED A-into-A-run slides -- a fully
base-AWARE ins_cost gives the IDENTICAL cheap A-into-A cost and would STILL demote them (~133,
still 22x full-run's 6). => refcol is STRUCTURALLY wrong for junction comparison, not fixably
wrong; it RELOCATES per-cut's run-splitting gaming into column-choice (drag-into-run) gaming.

SECONDARY defect (proven, minority): ins_col[j] is BASE-IDENTITY-INDEPENDENT (junction_scoring:848
`above=prev[j]+ins_col[j]`, never reads query[i-1]). 61/196 (31%) of demotions park a base-
MISMATCHED insertion (T=46,G=36,C=2) at a genome HP>=4 column and charge it the cheap poly-A/T
cost -- a base the pore did NOT stutter getting the stutter discount. Witness: 4 inserted C's next
to a genome 8-A run scored 0.7888 (= 4*ins_cost(8,'A')), IDENTICAL to 4 genuine A over-calls,
though ins_cost(8,'C')=0.5729 and ins_cost(1,'C')=1.25 exist. This is real but NOT the driver.

WHERE REFCOL IS RIGHT (concede): for an ISOLATED pore over-call the per-column charge is
biologically correct (expensive adjacent to non-HP: 4.58; cheap inside poly-A: 0.79). The task's
constructed dichotomy resolves in refcol's favor IN ISOLATION. The fault is that this
locally-correct per-column cost becomes a DRIFT-INDUCER in the junction comparison.

STRONGEST CHALLENGE TO MY OWN VERDICT (and why it fails): "just make ins_cost base-aware."
Killed by the 133 base-MATCHED slides: base-awareness leaves ~133 demotions (22x full-run).
The demotion class is structural (genome-column-only cost + no canonical-site protection), not a
base-identity bug.

STATUS: COMPLETE.
