export const meta = {
  name: 'microhom-guard-audit-v4-redundant',
  description: '4 Opus-Max auditors: 2 independent per task (discovery-loss, detector) → robust consensus',
  phases: [
    { title: 'Audit', detail: '4 independent Opus-Max auditors (2/task), own harness + own verdict, retry-on-stall' },
    { title: 'Synthesize', detail: 'consensus per task (A vs B agree?) + V2 byte-identity → all-clear verdict' },
  ],
}

const ROOT = '/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc'
const SCRATCH = '/private/tmp/claude-501/-Users-kevinroy-work-rectify/798e5203-2185-4c27-87b4-f4e50c56c3ba/scratchpad/audit_v4'

const COMMON = `
You are an ADVERSARIAL AUDITOR (Opus-Max) of RECTIFY's microhomology-drift guard AND its near-tie read-
evidence cap (commit 05664bc). YOU ARE ONE OF TWO INDEPENDENT AUDITORS assigned the SAME task for redundancy
and cross-check. Build your OWN harness, reach your OWN verdict — do NOT coordinate with or assume the other
auditor's result; independent agreement is the whole point. YOUR JOB IS TO FIND THE FAULT with real numbers.

Work READ-ONLY in ${ROOT} (branch worktree-agent-a25a2c1e784ad37dc). Do NOT edit repo code, do NOT commit,
do NOT touch defaults. Python: /Users/kevinroy/miniconda3/bin/python (repo importable from ROOT). Write your
own durable record + your own scratch harness ONLY under your assigned paths (given in your task) — do not
touch the other auditor's files (avoids write races).

GUARD MECHANICS (verify in rectify/core/splice/junction_refiner.py yourself):
- motif-blind re-placer + two READ-BLIND drift guards: hp_drift (_hp_run_across) + microhom_drift
  (_move_microhomology, ~lines 499-532) ADD to eff_margin when a move's boundary sits in a homopolymer /
  microhomology GENOMIC context. Veto (~line 800): move held at incumbent when
  best_score_cmp > incumbent_score - veto_margin. score_cmp LOWER=better;
  delta_improve := incumbent_score - best_score_cmp >= 0.
- CAP (05664bc): _effective_veto_margin(hold, eff, cap) = max(hold, min(eff, cap)) when active, else eff.
  drift_near_tie_cap default 0.0 = disabled = byte-identical. Threaded through refine_read_junctions /
  _run_sequential / _run_parallel / refine_bam_junctions.
- ESTABLISHED (V2 byte-identity leg, CLEARED — do not re-litigate): cap byte-identical off; structurally
  BOUNDS but cannot CLOSE the read-blind fault (delta_improve/eff/cap = one score axis; a real cryptic and a
  fab drift with the SAME delta_improve are treated identically inside (0,cap)). YOUR job is the EMPIRICAL
  magnitude / detector correctness, NOT the structural claim.

PRIOR PARTIALS (read for ORIENTATION only — the mechanistic insights — then build your OWN harness fresh;
do NOT resume another agent's harness, that would couple the two independent verdicts):
dev/MICROHOM_AUDIT_V2_discovery-loss.md, dev/MICROHOM_AUDIT_V2_detector-correctness.md,
dev/MICROHOM_AUDIT_V2_byte-identity.md, dev/MICROHOM_AUDIT_SYNTHESIS.md. Harnesses to LEARN FROM (reuse
patterns, not results): scripts/benchmark/noncanon_sim/{spikein_fab.py,spikein_score.py,mh_characterize.py,
build_panel.py,gen_reads.py}; tests/test_microhom_drift_guard.py.

DURABILITY MANDATE (NON-NEGOTIABLE — prior rounds lost verdicts to API stalls): WRITE your durable record
FIRST, APPEND a checkpoint after EACH sub-step (harness written -> panel persisted to disk -> first numbers ->
full grid -> verdict). Persist your harness SCRIPT + panel DATA to disk BEFORE any long sweep. Be ECONOMICAL:
minimal decisive panel (a few hundred reads/cell), NOT exhaustive — short runs stall less.
`

const SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['task', 'auditor_id', 'fault_found', 'strongest_challenge', 'evidence',
             'verdict', 'key_numbers', 'recommendation', 'durable_record'],
  properties: {
    task: { type: 'string', description: 'discovery-loss | detector-correctness' },
    auditor_id: { type: 'string', description: 'A or B' },
    fault_found: { type: 'boolean', description: 'true if a genuine HOLDING fault was found' },
    strongest_challenge: { type: 'string' },
    evidence: { type: 'string', description: 'concrete MEASURED numbers / reproducing cases / commands' },
    verdict: { type: 'string', enum: ['CLEAR', 'HOLD'], description: 'this auditor\'s independent call' },
    key_numbers: { type: 'string', description: 'the decisive quantities (discovery-loss rate + fab FDR per (margin,cap); or A5/A8 pass/fail) so the twin can be compared' },
    recommendation: { type: 'string', description: 'defensible (margin,cap) if any, or the blocking gap / detector fix' },
    durable_record: { type: 'string' },
  },
}

const TASKS = {
  'discovery-loss': `LENS = DISCOVERY-LOSS (LOAD-BEARING — the empirical rate the flip's safety depends on).
The structural fact is settled (cap bounds != closes); measure the MAGNITUDE independently:
  1. Build YOUR OWN read-cryptic discovery panel: REAL non-canonical (cryptic 3'SS/5'SS) reads whose TRUE
     site is genuinely discoverable (refiner with guard OFF moves the incumbent to it -> delta_improve>0) but
     whose boundary incidentally trips mh>=0.5 (short direct repeat). Distinguishing evidence MUST live
     OUTSIDE the soft-clipped drift prefix (else delta_improve=0 and the guard is moot). Sweep delta_improve
     via a controlled # of in-window bases. Own read builder; incumbent N-op at the wrong (canonical) site.
  2. Measure DISCOVERY-LOSS RATE (fraction of real cryptics wrongly vetoed) across
     microhom_drift_margin in {3,8} x drift_near_tie_cap in {off,1,2,3} x microhomology_frac. Per-read truth.
  3. Same grid: FAB suppression (mirror spikein_fab.py) so precision/recall is visible per (margin,cap).
  4. VERDICT (CLEAR/HOLD): at each (margin,cap), discovery-loss rate + fab FDR. Does cap@margin=3 bound
     discovery loss acceptably, or do real cryptics & fab drift OVERLAP in the (0,cap) delta_improve band
     (=> cap insufficient, positional signal required)? Give the GRID OF NUMBERS in key_numbers.`,
  'detector-correctness': `LENS = DETECTOR / HELPER CORRECTNESS. Adversarially test
_move_microhomology(seq,ns,ne,js,je) + _effective_veto_margin(hold,eff,cap) for inputs that score WRONG.
Run the full matrix, PRIORITIZING the two named-but-untested risks (may be real HOLDING faults):
  * A5: _frac_match scores N==N (any non-ACGT==non-ACGT) as a MATCH -> can a genomic ambiguity/repeat-N run
    yield FALSE-POSITIVE microhomology>=0.5 that vetoes a REAL move? Build the reproducing case; measure.
  * A8: max-over-both-boundaries masking -> a benign transition on one boundary masked by an unrelated repeat
    on the other trips >=0.5 and vetoes a move that should be spared. Build the case; measure.
  Then: donor vs acceptor geometry both strands; genome edges (lo-k<0, hi+k>n); k=0/1/huge; overlapping/
  nested repeats; both-boundary shift. AND _effective_veto_margin all regimes incl. hold>cap (veto never
  below hold) + cap<=0 / eff==hold exact no-ops. VERDICT (CLEAR/HOLD): is either helper INCORRECT? If A5/A8
  are real, they are HOLDING faults — reproducing case in evidence + key_numbers.`,
}

async function auditWithRetry(task, variant) {
  const MAX = 2
  const rec = `dev/MICROHOM_AUDIT_V4_${task}_${variant}.md`
  const scr = `${SCRATCH}/${task}_${variant}`
  for (let attempt = 1; attempt <= MAX; attempt++) {
    const prompt = `${COMMON}\n\nYOU ARE AUDITOR ${variant} ON TASK "${task}".\nYOUR DURABLE RECORD (yours alone): ${rec}\nYOUR SCRATCH DIR (yours alone, mkdir -p it): ${scr}\n[ATTEMPT ${attempt}/${MAX} — if your prior attempt stalled, read ${rec} + any harness you persisted in ${scr} and CONTINUE from the last checkpoint; do NOT restart from zero.]\n\n${TASKS[task]}\n\nReturn your INDEPENDENT structured verdict (task="${task}", auditor_id="${variant}").`
    const r = await agent(prompt, {
      label: `audit-v4:${task}/${variant} #${attempt}`, phase: 'Audit',
      model: 'opus', effort: 'max', schema: SCHEMA,
    })
    if (r) return { ...r, attempts: attempt }
    log(`audit-v4:${task}/${variant} attempt ${attempt}/${MAX} stalled — retrying from partial`)
  }
  return null
}

phase('Audit')
const jobs = []
for (const task of ['discovery-loss', 'detector-correctness'])
  for (const variant of ['A', 'B'])
    jobs.push(() => auditWithRetry(task, variant))
const results = (await parallel(jobs)).filter(Boolean)

phase('Synthesize')
const SYNTH_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['all_clear', 'final_call', 'discovery_loss_consensus', 'detector_consensus',
             'per_auditor', 'fault_to_fix', 'recommended_operating_point'],
  properties: {
    all_clear: { type: 'boolean' },
    final_call: { type: 'string' },
    discovery_loss_consensus: { type: 'string', description: 'do A & B agree? agreement level + reconciled finding' },
    detector_consensus: { type: 'string', description: 'do A & B agree? agreement level + reconciled finding' },
    per_auditor: { type: 'array', items: { type: 'string' } },
    fault_to_fix: { type: 'string' },
    recommended_operating_point: { type: 'string' },
  },
}
const synth = await agent(
  `You are the SYNTHESIZER for a REDUNDANT audit of the microhomology-guard + near-tie-cap (05664bc). FOUR
independent Opus-Max auditors ran: TWO on discovery-loss (A,B) and TWO on detector-correctness (A,B). The
BYTE-IDENTITY + structural-honesty lens already COMPLETED & CLEARED in V2 (read
dev/MICROHOM_AUDIT_V2_byte-identity.md + dev/MICROHOM_AUDIT_SYNTHESIS.md).

Read ALL FOUR durable records on disk (source of truth — some agents may have stalled; only .md persists):
  ${ROOT}/dev/MICROHOM_AUDIT_V4_discovery-loss_A.md
  ${ROOT}/dev/MICROHOM_AUDIT_V4_discovery-loss_B.md
  ${ROOT}/dev/MICROHOM_AUDIT_V4_detector-correctness_A.md
  ${ROOT}/dev/MICROHOM_AUDIT_V4_detector-correctness_B.md
Returned findings: ${JSON.stringify(results)}

YOUR JOB = ROBUST CONSENSUS. For EACH task, compare the two independent auditors:
  - Do A and B AGREE on the verdict (CLEAR/HOLD) and are their key_numbers consistent? Report the agreement
    level. If they DISAGREE, that is itself important — reconcile it (whose evidence is stronger and why),
    and treat unresolved disagreement as NOT-clear (a split verdict cannot license a flip).
  - A HOLD or a real fault from EITHER auditor on a task blocks that task (adversarial: one solid fault is
    enough). all_clear=true ONLY if: byte-identity cleared (V2) AND both discovery-loss auditors reach CLEAR
    with consistent numbers AND both detector auditors reach CLEAR — i.e. consensus CLEAR on every lens.
  - If discovery-loss shows the (0,cap) overlap is irreducible, state a positional-distinctiveness signal is
    REQUIRED to close (cap only bounds). If A5/A8 is a confirmed detector fault, that is a HOLDING fault.
State the concrete (margin, threshold, cap) IF the numbers defensibly support one, else the blocking gap.
COMPASS real-data confirmation remains an independent ship prerequisite regardless. Write your synthesis to
${ROOT}/dev/MICROHOM_AUDIT_V4_SYNTHESIS.md before returning.`,
  { label: 'audit-v4:synthesis', phase: 'Synthesize', model: 'opus', effort: 'max', schema: SYNTH_SCHEMA }
)

return { results, synth }
