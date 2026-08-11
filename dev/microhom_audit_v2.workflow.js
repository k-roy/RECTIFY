export const meta = {
  name: 'microhom-guard-fix-audit-v2',
  description: 'Triple Opus-Max adversarial audit of the microhomology-drift guard AND its near-tie-cap fix',
  phases: [
    { title: 'Audit', detail: '3 orthogonal Opus-Max auditors, heavy incremental checkpointing' },
    { title: 'Synthesize', detail: 'all-clear verdict from the 3 durable auditor records' },
  ],
}

// ---------------------------------------------------------------------------
// Shared context handed to every auditor (verbatim), so each is self-contained.
// ---------------------------------------------------------------------------
const ROOT = '/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc'
const COMMON = `
You are an ADVERSARIAL AUDITOR (Opus-Max) of RECTIFY's microhomology-drift guard AND its new
near-tie read-evidence cap. Per project policy this byte-identity-critical core-scoring change is
triple-audited before ANY default is flipped. YOUR JOB IS TO FIND THE FAULT, not to confirm.

Work READ-ONLY in ${ROOT} (branch worktree-agent-a25a2c1e784ad37dc). Do NOT edit code, do NOT commit,
do NOT touch defaults. Python: /Users/kevinroy/miniconda3/bin/python (or the repo's env). You MAY write
ONLY your own durable audit record file (below).

WHAT WAS BUILT (read the code yourself; do not trust this summary):
- rectify/core/splice/junction_refiner.py: motif-blind re-placer + two read-BLIND drift guards
  (hp_drift via _hp_run_across; microhom_drift via _move_microhomology) that ADD to eff_margin when a
  move's boundary sits in a homopolymer / microhomology genomic context. Veto rule (~line 800):
  a move is held at the incumbent when best_score_cmp > incumbent_score - veto_margin
  (score_cmp: LOWER = better; delta_improve := incumbent_score - best_score_cmp >= 0).
- NEW FIX (Phase 4, commit on this branch): _effective_veto_margin(hold_margin, eff_margin,
  drift_near_tie_cap) = max(hold_margin, min(eff_margin, drift_near_tie_cap)) when the cap is active,
  else eff_margin. drift_near_tie_cap default 0.0 = disabled = byte-identical. Threaded through
  refine_read_junctions / _run_sequential / _run_parallel / refine_bam_junctions.

THE FAULT THE FIX RESPONDS TO (from dev/MICROHOM_AUDIT_SYNTHESIS.md): _move_microhomology is READ-BLIND
(genome-only), so the mh>=0.5 veto TRIGGER (genomic) and the delta_improve that must clear the margin
(READ evidence) are INDEPENDENT. A real cryptic the read distinguishes (delta_improve>0, discovered with
guard OFF) can still trip mh>=0.5 and be vetoed whenever delta_improve < margin. The near-tie cap is
claimed to BOUND this (a move with delta_improve>=cap is never drift-vetoed) but NOT to close it (the cap
adds no discriminating signal inside the (0,cap) band — real cryptics and fab-drift overlap there).

DURABILITY MANDATE (NON-NEGOTIABLE — the last audit lost 2 of 3 agents to API stalls with unrecoverable
work): WRITE your durable record file FIRST (plan), then APPEND a checkpoint after EACH sub-step
(orientation done -> harness built -> first numbers -> verdict). The checkpoint cadence IS the recovery
granularity. Persist every number the moment it lands. If you stall, your file is the resume point.

Existing prior partial records to READ and build on (do NOT restart from scratch):
dev/MICROHOM_AUDIT_discovery-loss.md, dev/MICROHOM_AUDIT_detector-correctness.md,
dev/MICROHOM_AUDIT_byte-identity.md, dev/MICROHOM_AUDIT_SYNTHESIS.md,
dev/MICROHOMOLOGY_DRIFT_GUARD_DESIGN.md. Benchmark harnesses:
scripts/benchmark/noncanon_sim/ (build_panel.py, gen_reads.py, spikein_fab.py, spikein_score.py,
mh_characterize.py, q1_*.py). Regression tests: tests/test_microhom_drift_guard.py.

INDEPENDENCE: build your OWN adversarial cases. Do NOT accept the design doc's "zero discovery cost"
or "byte-identical" claims without independently reproducing them. If a claim is an artifact of an
under-powered panel, SAY SO with evidence.
`

const SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['lens', 'fault_found', 'strongest_challenge', 'challenge_holds', 'evidence',
             'verdict_survives', 'recommendation', 'durable_record'],
  properties: {
    lens: { type: 'string' },
    fault_found: { type: 'boolean', description: 'true if a genuine holding fault was found' },
    strongest_challenge: { type: 'string', description: 'the hardest attack you mounted' },
    challenge_holds: { type: 'boolean', description: 'does that challenge survive your own rebuttal attempt' },
    evidence: { type: 'string', description: 'concrete numbers / cases / commands that establish the finding' },
    verdict_survives: { type: 'boolean', description: 'does the guard+fix survive THIS lens (true=clear)' },
    recommendation: { type: 'string', description: 'ship / hold / specific change; incl. any (margin,cap) finding' },
    durable_record: { type: 'string', description: 'path to your appended durable audit .md' },
  },
}

phase('Audit')

const LENSES = [
  {
    key: 'discovery-loss',
    file: 'dev/MICROHOM_AUDIT_V2_discovery-loss.md',
    task: `LENS = DISCOVERY-LOSS (the LOAD-BEARING leg; the prior auditor STALLED at checkpoint 3 before
building the decisive case — RESUME it). The prior auditor proved: a PURE microhomology tie (delta_improve=0)
is guard-IRRELEVANT (scorer soft-clips the drift-distance prefix -> refiner tie-break holds incumbent
regardless of guard). Genuine guard-CAUSED discovery loss therefore requires a cryptic that scores STRICTLY
BETTER (delta_improve>0, discovered with guard OFF) yet still trips mh>=0.5 and is vetoed because
delta_improve < veto_margin. BUILD THAT PANEL (the prior auditor named it and stalled):

  1. Construct REAL non-canonical (cryptic 3'SS / 5'SS) reads whose true site sits next to a SHORT direct
     repeat (incidental microhomology mh>=0.5) at VARIED delta_improve (the read genuinely distinguishes the
     cryptic from the canonical/incumbent by a controlled number of in-window bases -> sweep delta_improve
     across ~0.5 .. 8+). Inject the CRYPTIC exon as the read query while the incumbent N-op is at the wrong
     (canonical) site, so the read WANTS to move (mirror _acceptor_after_refine in tests, but read-cryptic).
  2. Measure the DISCOVERY-LOSS RATE (fraction of real cryptics wrongly vetoed) as a function of
     (microhom_drift_margin in {3,8}) x (drift_near_tie_cap in {off,1,2,3}) x (microhomology_frac).
  3. In the SAME sweep, re-measure the fab suppression (spike-in fab panel) so you see the precision/recall
     tradeoff at each (margin,cap). Use per-read truth, seed excluded.
  4. VERDICT: does the near-tie cap at margin=3 BOUND discovery loss to an acceptable level? Is there a
     (margin,cap) with fab~0 AND ~zero discovery loss, or is the (0,cap) overlap IRREDUCIBLE with an
     aggregate-delta signal (=> the cap is insufficient and a DIFFERENT signal — per-base positional
     distinctiveness the incumbent soft-clip cannot absorb — is required)? Give the numbers, not a hunch.
  Recommend a concrete (margin, cap) operating point IF one is defensible, else state the cap is insufficient.`,
  },
  {
    key: 'detector-correctness',
    file: 'dev/MICROHOM_AUDIT_V2_detector-correctness.md',
    task: `LENS = DETECTOR / HELPER CORRECTNESS (the prior auditor STALLED at PLAN — never ran; do it).
Two pure functions decide the guard: _move_microhomology(seq, ns, ne, js, je) and
_effective_veto_margin(hold_margin, eff_margin, drift_near_tie_cap). ADVERSARIALLY find inputs where each
scores WRONG:
  A. _move_microhomology: false-neg (misses real drift) / false-pos (flags a real transition as drift).
     Edge matrix: donor-side vs acceptor-side shifts; genome edges (lo-k<0, hi+k>n); k=0, k=1, huge k;
     non-ACGT / N / lowercase bases; overlapping & nested repeats; BOTH-boundary shift; the max-over-both
     masking case (a benign shift on one boundary masked by an unrelated repeat on the other — does taking
     the MAX over the two moved boundaries cause a spurious veto?). Verify the acceptor(downstream-repeat)
     vs donor(upstream-repeat) geometry is correct for BOTH strands.
  B. _effective_veto_margin: verify max(hold, min(eff, cap)) is correct across ALL regimes, especially the
     hold_margin>cap interaction (a drift-flagged move with delta_improve in (cap, hold_margin): is it still
     vetoed by hold, as intended? confirm the veto never drops below hold_margin). Confirm cap<=0 and
     eff==hold are exact no-ops (byte-identical). Try to find a parameter combo where the helper produces a
     veto_margin that is nonsensical (negative, > eff, or that changes behavior at default).
  Reproduce the unit tests in tests/test_microhom_drift_guard.py and try to BREAK them with a case they miss.`,
  },
  {
    key: 'byte-identity-and-fix',
    file: 'dev/MICROHOM_AUDIT_V2_byte-identity.md',
    task: `LENS = BYTE-IDENTITY AT DEFAULT + FIX STRUCTURAL HONESTY. Two duties:
  A. Re-verify the guard+cap is TRULY INERT at default (drift_near_tie_cap=0.0, all margins 0.0): output
     byte-identical to the pre-fix code. Independent BAM diff vs the pre-fix commit on a sim panel
     (mix_r3b / mix_fair / fab_sweep), SEQUENTIAL path with raw record order preserved (SHA256 over full
     SAM record text) AND the parallel path position-sorted (raw parallel order is known guard-independent
     scheduling noise — do NOT falsely flag it). Confirm the new _effective_veto_margin call does not run /
     has no effect at default. pytest -m "not slow" must be green. Isolate confounders via the DIRECT parent
     commit (only junction_refiner.py + tests should differ).
  B. FIX STRUCTURAL HONESTY: the code + design doc CLAIM the cap "bounds but does not close" the read-blind
     fault and adds "no discriminating signal inside (0,cap)". Independently verify this is HONEST, not an
     over- or under-claim: is delta_improve/eff_margin/cap really one scalar axis (=> cap == reducing the
     effective margin, no new signal)? Does the cap EVER protect a fab-drift move it shouldn't (i.e. does
     raising the cap re-admit fabrication)? Confirm the cap only ever LOWERS the veto margin (never raises
     it), so it can only ADD moves, never remove a real discovery relative to no-cap. State plainly whether
     the honest-framing claim in _effective_veto_margin's docstring is correct.`,
  },
]

const findings = await parallel(LENSES.map((L) => () =>
  agent(`${COMMON}\n\nYOUR DURABLE RECORD FILE: ${L.file} (write plan first, append a checkpoint per sub-step).\n\n${L.task}`,
    { label: `audit-v2:${L.key}`, phase: 'Audit', model: 'opus', effort: 'max', schema: SCHEMA })
))

phase('Synthesize')

const SYNTH_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['all_clear', 'final_call', 'per_auditor', 'fault_to_fix', 'recommended_operating_point'],
  properties: {
    all_clear: { type: 'boolean', description: 'true ONLY if all 3 auditors COMPLETED and none found a holding fault' },
    final_call: { type: 'string', description: 'SHIP-ENABLED (flip default) vs HOLD, with the decisive reason' },
    per_auditor: { type: 'array', items: { type: 'string' } },
    fault_to_fix: { type: 'string', description: 'the genuine holding fault if any, else "none"' },
    recommended_operating_point: { type: 'string', description: 'concrete (margin, threshold, cap) IF defensible, else the blocking gap' },
  },
}

const synth = await agent(
  `You are the SYNTHESIZER for the triple Opus-Max audit of the microhomology-drift guard + near-tie-cap fix.
Read ALL THREE durable auditor records on disk (they are the source of truth, more complete than any
returned summary — some agents may have STALLED and only their .md persists):
  ${ROOT}/dev/MICROHOM_AUDIT_V2_discovery-loss.md
  ${ROOT}/dev/MICROHOM_AUDIT_V2_detector-correctness.md
  ${ROOT}/dev/MICROHOM_AUDIT_V2_byte-identity.md
Also read dev/MICROHOM_AUDIT_SYNTHESIS.md (the prior HOLD) for continuity.

Returned auditor findings: ${JSON.stringify(findings)}

POLICY: all_clear=true ONLY if all three auditors COMPLETED (reached a verdict) AND none found a holding
fault. An incomplete/stalled audit is NOT all-clear (it failed to run) — say so explicitly and name which
leg is missing so the orchestrator can relaunch it from its partial .md. The default flip requires
all_clear=true. If discovery-loss shows the near-tie cap is INSUFFICIENT (the (0,cap) overlap is
irreducible), the call is HOLD + specify that a positional-distinctiveness signal is needed. If a
defensible (margin, threshold, cap) exists with fab~0 and bounded discovery loss, state it precisely.
Write your synthesis to ${ROOT}/dev/MICROHOM_AUDIT_V2_SYNTHESIS.md before returning.`,
  { label: 'audit-v2:synthesis', phase: 'Synthesize', model: 'opus', effort: 'max', schema: SYNTH_SCHEMA }
)

return { findings, synth }
