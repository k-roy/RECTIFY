export const meta = {
  name: 'microhom-guard-audit-v3-retry',
  description: 'Retry the 2 stalled load-bearing audit legs (discovery-loss + detector) with stall-retry + resume-from-partial',
  phases: [
    { title: 'Audit', detail: 'discovery-loss + detector-correctness, each retried on stall, resuming from partial .md' },
    { title: 'Synthesize', detail: 'combine V2 byte-identity (cleared) + V3 completed legs → all-clear verdict' },
  ],
}

const ROOT = '/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc'

const COMMON = `
You are an ADVERSARIAL AUDITOR (Opus-Max) of RECTIFY's microhomology-drift guard AND its near-tie
read-evidence cap (commit 05664bc). Per project policy this byte-identity-critical core-scoring change is
triple-audited before ANY default is flipped. YOUR JOB IS TO FIND THE FAULT with real numbers, not confirm.
This is a RETRY: two prior rounds STALLED this leg mid-stream (API errors) before it produced a verdict.

Work READ-ONLY in ${ROOT} (branch worktree-agent-a25a2c1e784ad37dc). Do NOT edit repo code, do NOT commit,
do NOT touch defaults. Python: /Users/kevinroy/miniconda3/bin/python (repo importable from ROOT). You MAY
write your own durable audit record + your own scratch harness under
/private/tmp/claude-501/-Users-kevinroy-work-rectify/798e5203-2185-4c27-87b4-f4e50c56c3ba/scratchpad/ .

GUARD MECHANICS (verify in rectify/core/splice/junction_refiner.py yourself):
- motif-blind re-placer + two READ-BLIND drift guards: hp_drift (_hp_run_across) and microhom_drift
  (_move_microhomology, lines ~499-532) ADD to eff_margin when a move's boundary sits in a homopolymer /
  microhomology GENOMIC context. Veto (~line 800): move held at incumbent when
  best_score_cmp > incumbent_score - veto_margin. score_cmp LOWER=better;
  delta_improve := incumbent_score - best_score_cmp >= 0.
- NEW CAP (05664bc): _effective_veto_margin(hold, eff, cap) = max(hold, min(eff, cap)) when active, else eff.
  drift_near_tie_cap default 0.0 = disabled = byte-identical. In refine_read_junctions /
  _run_sequential / _run_parallel / refine_bam_junctions.
- ESTABLISHED (V2 byte-identity leg, CLEARED — do not re-litigate): the cap is byte-identical off, and
  structurally BOUNDS but cannot CLOSE the read-blind fault (delta_improve/eff/cap = one score axis; a real
  cryptic and a fab drift with the SAME delta_improve are treated identically inside (0,cap)). YOUR job is
  the EMPIRICAL magnitude, not the structural claim.

DURABILITY MANDATE (NON-NEGOTIABLE — this is why prior rounds lost the verdict): WRITE your durable record
FIRST, then APPEND a checkpoint after EACH sub-step (harness written -> panel built + persisted to disk ->
first numbers -> full grid -> verdict). Persist your harness SCRIPT + panel DATA to disk BEFORE running any
long sweep, so a stall resumes from built artifacts, never from scratch. Persist every number the moment it
lands. Be ECONOMICAL: build the MINIMAL decisive panel (a few hundred reads per cell is plenty), not an
exhaustive one — a shorter run is far less likely to stall.

RESUME: your prior partial record already exists on disk (path in your task). READ it, find the last
checkpoint, and CONTINUE from there. If you already persisted a harness/panel file, REUSE it — do NOT rebuild.
`

const SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['lens', 'fault_found', 'strongest_challenge', 'challenge_holds', 'evidence',
             'verdict_survives', 'recommendation', 'durable_record'],
  properties: {
    lens: { type: 'string' },
    fault_found: { type: 'boolean' },
    strongest_challenge: { type: 'string' },
    challenge_holds: { type: 'boolean' },
    evidence: { type: 'string', description: 'concrete measured numbers / cases / commands' },
    verdict_survives: { type: 'boolean', description: 'does the guard+cap survive THIS lens (true=clear)' },
    recommendation: { type: 'string' },
    durable_record: { type: 'string' },
  },
}

const LEGS = [
  {
    key: 'discovery-loss',
    file: 'dev/MICROHOM_AUDIT_V3_discovery-loss.md',
    prior: 'dev/MICROHOM_AUDIT_V2_discovery-loss.md',
    task: `LENS = DISCOVERY-LOSS (LOAD-BEARING — the empirical rate the flip's safety depends on; stalled at
CP1/orientation twice). The structural fact is settled (cap bounds != closes); your job is the MAGNITUDE.
DO THE EMPIRICAL WORK (PLAN steps 3-5 the prior rounds never ran):
  1. Build a READ-CRYPTIC discovery panel: REAL non-canonical (cryptic 3'SS/5'SS) reads whose TRUE site is
     genuinely discoverable (the refiner with guard OFF moves the incumbent to it -> delta_improve>0), but
     whose boundary sits next to a SHORT direct repeat so it incidentally trips mh>=0.5. The distinguishing
     evidence MUST live OUTSIDE the soft-clipped drift prefix (else delta_improve=0 and the guard is moot —
     the prior auditor's key mechanistic finding). Sweep delta_improve via a controlled # of in-window bases.
     Mirror spikein_fab.py's build (its microhom_mismatch knob is the same microhomology axis) but with truth
     at the CRYPTIC site. Own read builder; incumbent N-op at the canonical/wrong site. Persist builder+panel.
  2. Measure the DISCOVERY-LOSS RATE (fraction of real cryptics wrongly vetoed = held at incumbent) across
     microhom_drift_margin in {3,8} x drift_near_tie_cap in {off,1,2,3} x microhomology_frac. Per-read truth.
  3. On the SAME grid re-measure FAB suppression (spikein_fab panel) so the precision/recall tradeoff is
     visible per (margin,cap).
  4. VERDICT: at each (margin,cap), what is the discovery-loss rate and the fab FDR? Does the cap at margin=3
     bound discovery loss to an acceptable level, or do real cryptics and fab drift OVERLAP in the (0,cap)
     delta_improve band (=> cap insufficient, positional signal required)? Give the grid of NUMBERS + a
     concrete recommendation (a defensible (margin,cap) if one exists, else "cap insufficient, need positional").
RESUME from ${ROOT}/dev/MICROHOM_AUDIT_V2_discovery-loss.md (@CP1). Append into ${ROOT}/dev/MICROHOM_AUDIT_V3_discovery-loss.md.`,
  },
  {
    key: 'detector-correctness',
    file: 'dev/MICROHOM_AUDIT_V3_detector-correctness.md',
    prior: 'dev/MICROHOM_AUDIT_V2_detector-correctness.md',
    task: `LENS = DETECTOR / HELPER CORRECTNESS (stalled at PLAN twice — never ran; RUN it). Adversarially test
_move_microhomology(seq, ns, ne, js, je) and _effective_veto_margin(hold, eff, cap) for inputs that score WRONG.
Run the full A1-A12 matrix, PRIORITIZING the two named-but-untested risks first (they may be real faults):
  * A5: _frac_match scores N==N (and any non-ACGT==non-ACGT) as a MATCH -> does a genomic ambiguity/repeat-N
    run produce a FALSE-POSITIVE microhomology >= 0.5 that vetoes a real move? Build the case; measure.
  * A8: max-over-both-boundaries masking — a benign transition on ONE boundary masked by an unrelated repeat
    on the OTHER, so the MAX trips >=0.5 and vetoes a move that should be spared. Build the case; measure.
  Then the rest: donor vs acceptor geometry both strands; genome edges (lo-k<0, hi+k>n); k=0/1/huge;
  overlapping/nested repeats; both-boundary shift. AND _effective_veto_margin across all regimes incl.
  hold>cap (veto must never drop below hold_margin) and cap<=0 / eff==hold exact no-ops.
  VERDICT: is either helper INCORRECT (a false-pos veto that kills real discovery, or a false-neg that lets
  fab through)? If A5/A8 are real, they are HOLDING faults — say so with the reproducing case.
RESUME from ${ROOT}/dev/MICROHOM_AUDIT_V2_detector-correctness.md (@PLAN). Append into ${ROOT}/dev/MICROHOM_AUDIT_V3_detector-correctness.md.`,
  },
]

async function auditWithRetry(leg) {
  const MAX = 3
  for (let attempt = 1; attempt <= MAX; attempt++) {
    const prompt = `${COMMON}\n\nYOUR DURABLE RECORD: ${leg.file}. YOUR PRIOR PARTIAL (resume from it): ${leg.prior}.\n[ATTEMPT ${attempt}/${MAX} — if a prior attempt stalled, read your ${leg.file} / ${leg.prior} + any persisted harness on disk and CONTINUE from the last checkpoint; do NOT restart.]\n\n${leg.task}`
    const r = await agent(prompt, {
      label: `audit-v3:${leg.key} #${attempt}`, phase: 'Audit',
      model: 'opus', effort: 'max', schema: SCHEMA,
    })
    if (r) return { ...r, attempts: attempt }
    log(`audit-v3:${leg.key} attempt ${attempt}/${MAX} stalled (null) — retrying from partial`)
  }
  log(`audit-v3:${leg.key} exhausted ${MAX} attempts — leaving partial for orchestrator relaunch`)
  return null
}

phase('Audit')
const legResults = await parallel(LEGS.map((leg) => () => auditWithRetry(leg)))

phase('Synthesize')
const SYNTH_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['all_clear', 'final_call', 'per_auditor', 'fault_to_fix', 'recommended_operating_point'],
  properties: {
    all_clear: { type: 'boolean' },
    final_call: { type: 'string' },
    per_auditor: { type: 'array', items: { type: 'string' } },
    fault_to_fix: { type: 'string' },
    recommended_operating_point: { type: 'string' },
  },
}

const synth = await agent(
  `You are the SYNTHESIZER for the microhomology-guard + near-tie-cap audit. The BYTE-IDENTITY +
structural-honesty leg already COMPLETED and CLEARED in V2 (read dev/MICROHOM_AUDIT_V2_byte-identity.md +
dev/MICROHOM_AUDIT_V2_SYNTHESIS.md) — cap is byte-identical off and honestly framed (bounds != closes). This
V3 round RETRIED the two previously-stalled legs. Read ALL relevant durable records on disk (source of truth,
more complete than any returned summary):
  ${ROOT}/dev/MICROHOM_AUDIT_V3_discovery-loss.md   (+ V2 predecessor for continuity)
  ${ROOT}/dev/MICROHOM_AUDIT_V3_detector-correctness.md   (+ V2 predecessor)
  ${ROOT}/dev/MICROHOM_AUDIT_V2_byte-identity.md
Returned V3 leg findings: ${JSON.stringify(legResults)}

POLICY: all_clear=true ONLY if ALL THREE lenses (byte-identity [V2, cleared] + discovery-loss [V3] +
detector-correctness [V3]) COMPLETED with no holding fault. A still-stalled V3 leg is NOT all-clear — name it
so the orchestrator can relaunch from its partial. Decide, with the measured numbers:
  (1) Is the guard+cap safe to enable as a NON-default (opt-in) with documented (margin,cap)?
  (2) Is any DEFAULT flip defensible, or does discovery-loss / A5 / A8 block it?
  (3) If discovery-loss shows the (0,cap) overlap is irreducible, state that a positional-distinctiveness
      signal is REQUIRED to close (cap only bounds). If A5 or A8 is a real detector fault, that is a HOLDING
      fault — state the fix.
State the concrete (margin, threshold, cap) IF one is defensible from the numbers, else the blocking gap.
Note COMPASS real-data confirmation remains an independent ship prerequisite regardless. Write your synthesis
to ${ROOT}/dev/MICROHOM_AUDIT_V3_SYNTHESIS.md before returning.`,
  { label: 'audit-v3:synthesis', phase: 'Synthesize', model: 'opus', effort: 'max', schema: SYNTH_SCHEMA }
)

return { legResults, synth }
