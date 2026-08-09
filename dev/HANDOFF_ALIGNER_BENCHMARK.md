# SCOPED HANDOFF — Deliverable A: the simulation ground-truth benchmark (the GATE)

**For:** a dedicated agent building the native-aligner program's benchmark, in PARALLEL with the
main session (which is re-deriving the COMPASS 111 verdict).
**Created:** 2026-06-26 by the COMPASS/short-read session. **Branch policy: ISOLATED WORKTREE only**
(see Constraints) — the shared `drs-validation-rebuild` branch has a concurrent DRS-arm agent; do NOT
commit there.

---

## OBJECTIVE (one sentence)
Build the **validated simulation ground-truth benchmark** that is the documented GATE for any RECTIFY
native-aligner member code (facet C1 = "de-novo realign within the panel-mapped locus, scored by the
empirical HpPenaltyTable"). No member code is built until this benchmark exists and its ablations are
runnable. **Fitness = the simulation truth set, NEVER the internal score** (a pure-scoring re-weight
once flipped GMAP annotated win/loss 0.09→1.07 with no aligner change — that artifact is the entire
reason this gate exists).

## READ FIRST (canonical design — do not re-derive)
- `dev/ALIGNER_PROGRAM_SCOPING.md` §"Deliverable A" — the benchmark's reuse/build split + sequencing.
- `dev/ALIGNER_MEMBER_DESIGN.md` §9 "Benchmark coupling" — the per-facet `benchmark_requirement`
  blocks (C1–C5 truth/stratum requirements) and §7 the shared train/test split discipline.
- `dev/ALIGNER_IDEATION_SYNTHESIS.md` §"Benchmark requirements" — headline metric spec.
- `docs/EMPIRICAL_HP_PENALTY_SCORING.md` — the empirical tables the member will score with (context).
- `dev/SIMULATION_BENCHMARK_SPEC.md` if present (C6 variant stratum) — fold in.

## THE DELIVERABLE — three build components (all ABSENT today; rest is reuse)
1. **Read-simulator wrapper** — drive pbsim3 / badread / nanosim (DRS + cDNA error models) over an
   LRGASP-style transcript GFF → FASTQ. No simulator is referenced anywhere in the repo yet; pick one,
   justify, wrap it. Must emit reads with a realistic ONT DRS length-error model (for the C1 HP-length law).
2. **Per-read junction/indel TRUTH propagation** — carry each simulated read's KNOWN donor/acceptor
   coords + NIC/NNC class + exact indel positions as SAM aux tags / a sidecar truth table. (The existing
   `XV` tag is a correction-CATEGORY label, NOT junction truth — do not overload it.)
3. **Per-junction + per-indel accuracy SCORER** — TP/FP/FN per aligner, stratified by junction class
   (canonical/annotated/NIC/NNC) AND read-class, against truth. **MUST** use the ambiguity-aware match
   (`normalize_junction` / `_canonical_within_window` / `junction_ambiguity_window`,
   `rectify/core/consensus/chimeric_consensus.py:59-155`) so a call 1bp into a donor/acceptor repeat is
   NOT charged FP (the exact trap that produced the GMAP 0.09 artifact). **Framing metric = EXACT
   INDEL-POSITION CONCORDANCE WITH TRUTH, ambiguity-aware — NEVER edit distance** (edit distance ties by
   construction at the contested HP/repeat positions).

## MUST-HAVES (from the design's §7/§9 — non-negotiable for the gate to be valid)
- **NIC/NNC junction truth labels** present from the start (else §8 discovery-FDR is unmeasurable).
- **One genomic-region-disjoint TRAIN/TEST split**, tagged per read, SHARED across all facets.
- **C1 cell coverage:** every `(HP_run_length 1-12, base_class AT/CG)` and STR `(unit, n_copies)` cell
  sized to clear `min_count=100` (`hp_penalty.py:184,213`) — else the length-law silently nullifies to
  flat and produces a FALSE REFUTE. Include clean (no-error) runs to measure the false-indel rate.
- **Strata to contain** (from scoping doc): HP runs A/C/G/T×1-12; genomic-A-abutting CPA; NIC/NNC
  canonical+noncanonical novels; edit-distance-tied STR positions; paralog loci (SMN1/SMN2-style);
  constructed panel-failure reads (size the discovery-ceiling tail); coverage strata + calibrated/
  miscalibrated phred deciles. (C6 variant stratum if SIMULATION_BENCHMARK_SPEC.md calls for it.)

## REUSE (already in-tree — wire, don't rebuild)
- `scripts/benchmark/aligner_contribution_analysis.py` — per-aligner win-rate cross-tab (report layer scaffold).
- `rectify/core/commands/validate_command.py` — 3'-end (CPA) truth scoring pattern (`CorrectedPosition`,
  `ValidationResult`); mirror for junctions (it does NOT do junctions today).
- `rectify/core/splice/junction_validator.py` (COMPASS 3-pass) + `junction_scoring.py` motif/anchor
  classification — for the scorer's junction-call side.
- `rectify/core/consensus/chimeric_consensus.py:59-155` — the ambiguity-aware match primitives (above).

## CONSTRAINTS
- **Run heavy work on Sherlock** (data locality, AVX-512 for the rectify env), chunked, owners partition,
  CPU-constraint per `sherlock-sbatch` skill traps. Do NOT relay BAMs through the M1.
- **Isolated worktree** (this task was launched `isolation: worktree`). Stage only NEW files under
  `scripts/benchmark/`, `dev/`, or a NEW `rectify/core/benchmark/` — never edit shared hot files
  (`rectify/core/bam/*`, `correct/*`, `consensus/*`, `splice/*`) that the DRS agent may be touching.
- Coordinate via `.claude/inbox/` if you need anything from the main session.
- **No member (C1) code** — that is the NEXT cycle, gated on THIS benchmark passing. Phase 0 of C1 (the
  `penalty_table=None` equivalence harness) may be PROTOTYPED as a consumer of the benchmark but not shipped.

## RESUME (concrete first step)
1. Read the three design docs above (esp. `ALIGNER_MEMBER_DESIGN.md` §9 + `ALIGNER_PROGRAM_SCOPING.md`).
2. Survey what simulators are installable on Sherlock (`module avail`, conda) — pbsim3 vs badread vs
   nanosim — and pick one; record the choice + DRS/cDNA error-model justification in
   `dev/SIMULATION_BENCHMARK_SPEC.md` (create/extend).
3. Build component (2) the truth-propagation schema FIRST (it constrains 1 and 3), then the simulator
   wrapper (1), then the ambiguity-aware scorer (3) reusing the chimeric_consensus primitives.
4. Validate end-to-end on a tiny transcript set (e.g. a few chr5 genes) before scaling: confirm a known
   junction round-trips truth→sim-read→aligner→scorer as TP, and a deliberately shifted call as the
   correct TP-not-FP under the ambiguity-aware match.
5. Drop a sentinel + refresh THIS handoff when the smoke round-trips.

## DONE / VERIFIED / OPEN (seed — update as you go)
- DONE: nothing yet (fresh scope).
- VERIFIED: the gate philosophy + reuse/build split are from the user-approved scoping doc.
- OPEN: all three components; simulator choice; Sherlock simulator availability.
