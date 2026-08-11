export const meta = {
  name: 'microhom-positional-close-audit-v5',
  description: '8 Opus-Max auditors (2/task ×4) rigorously vet the positional-distinctiveness close',
  phases: [
    { title: 'Audit', detail: '8 independent Opus-Max auditors (2 per task), own harness, retry-on-stall' },
    { title: 'Synthesize', detail: 'consensus per task → is the close CORRECT, COMPLETE, and worth enabling' },
  ],
}

const ROOT = '/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc'
const SCRATCH = '/private/tmp/claude-501/-Users-kevinroy-work-rectify/798e5203-2185-4c27-87b4-f4e50c56c3ba/scratchpad/audit_v5'

const COMMON = `
You are an ADVERSARIAL AUDITOR (Opus-Max) vetting RECTIFY's POSITIONAL-DISTINCTIVENESS CLOSE for the
motif-blind junction re-placer's microhomology-drift guard. YOU ARE ONE OF TWO INDEPENDENT AUDITORS on the
SAME task — build your OWN harness, reach your OWN verdict, do NOT coordinate. FIND THE FAULT.

Work READ-ONLY in ${ROOT} (branch worktree-agent-a25a2c1e784ad37dc). Do NOT edit repo code / commit / touch
defaults. Python /Users/kevinroy/miniconda3/bin/python (repo importable from ROOT). Write ONLY your own
durable record + scratch (paths in your task).

WHAT IS UNDER AUDIT (read the code + docs yourself; verify, do not trust this summary):
- rectify/core/splice/junction_refiner.py: the motif-blind re-placer + drift guards. Recent commits:
  05664bc near-tie cap (drift_near_tie_cap); d1fd08d detector fixes (A5 ACGT-only frac/HP, A8 min-over-moved);
  cd3de46 the POSITIONAL CLOSE; a97ff6d donor investigation (acceptor-only by design).
- The CLOSE: _positional_signal(genome, q, q_split, ne, new_je) = _semiglobal_ed(rescue, genome[ne:]) -
  _semiglobal_ed(rescue, genome[new_je:]), rescue = q[q_split:]. Hard-anchored edit distance (free ref
  suffix, NO free-prefix split) — the missing free-k split removes the scorer's soft-clip escape that hid
  the discriminating exon2 bases from delta_improve. Wired as drift_positional_gate (default 0.0 = OFF =
  byte-identical): a drift-flagged would-be-veto is SPARED when _positional_signal >= gate.
- CLAIMED result (dev/DISCOVERY_LOSS_PANEL_RESULT.md, dev/discovery_loss_panel.py): margin=3 alone loses
  ~24% discovery; cap=2 ~10%; the positional signal (ed) separates the delta overlap band at 98-99%
  balanced accuracy; WIRED m3/cap2/gate1 → ~0.4% discovery loss / ~4.3% fab-residual = "the fault CLOSED".
- CONTEXT: the guard is DEFAULT-OFF and UNWIRED (correct_command.py passes no drift kwargs); validated on
  SYNTHETIC panels only; COMPASS real-data is the untested enablement gate. The refiner is ACCEPTOR-CENTRIC
  (_score_junction ignores intron_start/donor) → the close is acceptor-only by design.

DURABILITY MANDATE (prior audit rounds lost verdicts to API stalls): WRITE your durable record FIRST, APPEND
a checkpoint after EACH sub-step (harness written -> data persisted -> numbers -> verdict). Persist harness +
data to disk BEFORE any long run. Be ECONOMICAL (minimal decisive panel). Independent agreement is the point.
`

const SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['task', 'auditor_id', 'fault_found', 'strongest_challenge', 'evidence',
             'verdict', 'key_numbers', 'recommendation', 'durable_record'],
  properties: {
    task: { type: 'string' }, auditor_id: { type: 'string' },
    fault_found: { type: 'boolean' },
    strongest_challenge: { type: 'string' },
    evidence: { type: 'string', description: 'concrete MEASURED numbers / reproducing cases / commands' },
    verdict: { type: 'string', enum: ['CLEAR', 'HOLD'] },
    key_numbers: { type: 'string', description: 'decisive quantities so the twin can be compared' },
    recommendation: { type: 'string' },
    durable_record: { type: 'string' },
  },
}

const TASKS = {
  'signal-correctness': `TASK = POSITIONAL-SIGNAL + _semiglobal_ed CORRECTNESS. Adversarially test that
_semiglobal_ed and _positional_signal are CORRECT and do what the docstring claims. Probe: (a) _semiglobal_ed
edge cases — empty query/ref, query longer than ref, all-mismatch, free-suffix behaviour, indels vs
substitutions, is it truly hard-anchored at ref[0]? (b) _positional_signal — window W truncation, rescue
running off genome end, ties (signal=0), the None cases; does >0 REALLY mean "read matches moved acceptor
better"? Build reads where it gives the WRONG sign (false-positive spare of a fab drift, or false-negative
veto of a real cryptic). (c) MINUS-STRAND reads: is the acceptor/exon2 geometry correct for is_reverse reads
(genome coords vs transcript orientation)? (d) the acceptor-centric assumption — construct a both-boundary
move where ignoring the donor causes a wrong spare. VERDICT CLEAR/HOLD with reproducing cases.`,
  'independent-remeasure': `TASK = INDEPENDENT DISCOVERY-LOSS RE-MEASUREMENT. Do NOT reuse
dev/discovery_loss_panel.py — build your OWN panel with a DIFFERENT construction/error model and independently
re-measure the close. Questions: (1) Does the ~0.4% discovery-loss / ~4.3% fab-residual at m3/cap2/gate1 hold
on YOUR panel? (2) Is the original panel's construction BIASED — e.g. does U+U'+TAIL over-represent
easy-to-separate cases; what about longer repeats (k>10), imperfect tails, higher ONT error (10-15%),
homopolymer-adjacent boundaries, multi-candidate pools (not just 2)? (3) At what microhomology / error regime
does the ed signal's separation DEGRADE below usefulness? Report the frontier. VERDICT: is "CLOSED" robust or
construction-dependent?`,
  'byte-identity-architecture': `TASK = BYTE-IDENTITY + ARCHITECTURE CRITIQUE. (a) Re-verify the guard+cap+
positional stack is byte-identical at the default (all of microhom_drift_margin/hp_drift_margin/
drift_near_tie_cap/drift_positional_gate = 0.0): independent BAM diff vs the pre-close parent, sequential
(raw-order SHA256) + parallel (position-sorted). pytest -m "not slow" green. Confirm _positional_signal /
_semiglobal_ed are never called at default. (b) ARCHITECTURE: is a 4th stacked guard PARAMETER
(margin/threshold/cap/positional_gate) the right design, or should the discrimination live IN the scorer
(fix _score_junction's free-k soft-clip so delta_improve itself discriminates)? Is computing a SECOND
alignment (edit distance) of the same rescue redundant with what the scorer already computes? Does the gate
interact cleanly with the cap (spare-then-veto ordering)? Assess complexity/maintainability vs the benefit on
a default-off feature.`,
  'strategic-realdata': `TASK = STRATEGIC / REAL-DATA VALIDITY (adversarial "should this ship" review). The
close is validated on SYNTHETIC panels only; the guard is default-OFF and unwired; it suppresses a ~1-2%
fabrication rate on a tool whose value is DISCOVERY (the COMPASS 32x recall win). Adversarially assess: (1)
Is the synthetic panel representative of REAL ONT cryptic-in-microhomology loci (paralogs, ribosomal, the SMA
fabrication-enriched set)? What real-data property could break the ed signal that no synthetic panel tested?
(2) What EXACTLY must COMPASS real-data show to justify enabling the guard, and is that measurable with the
existing COMPASS harness (dev/COMPASS_RECALL_RESULT.md)? (3) Is the guard NET-POSITIVE for a discovery tool
at all — does suppressing ~2% fabrication justify ANY discovery risk + 4 params of complexity? (4) Steelman
the case that the honest stop is "acceptor close done, default-OFF, gate on COMPASS" vs "keep building".
VERDICT: is the close's VALIDITY established, or is it synthetic-only insurance whose worth is unproven?`,
}

async function auditWithRetry(task, variant) {
  const MAX = 2
  const rec = `dev/MICROHOM_AUDIT_V5_${task}_${variant}.md`
  const scr = `${SCRATCH}/${task}_${variant}`
  for (let attempt = 1; attempt <= MAX; attempt++) {
    const prompt = `${COMMON}\n\nYOU ARE AUDITOR ${variant} ON TASK "${task}".\nYOUR DURABLE RECORD (yours alone): ${rec}\nYOUR SCRATCH DIR (yours alone, mkdir -p it): ${scr}\n[ATTEMPT ${attempt}/${MAX} — if your prior attempt stalled, read ${rec} + any harness in ${scr} and CONTINUE from the last checkpoint.]\n\n${TASKS[task]}\n\nReturn your INDEPENDENT structured verdict (task="${task}", auditor_id="${variant}").`
    const r = await agent(prompt, {
      label: `audit-v5:${task}/${variant} #${attempt}`, phase: 'Audit',
      model: 'opus', effort: 'max', schema: SCHEMA,
    })
    if (r) return { ...r, attempts: attempt }
    log(`audit-v5:${task}/${variant} attempt ${attempt}/${MAX} stalled — retrying from partial`)
  }
  return null
}

phase('Audit')
const TASK_KEYS = ['signal-correctness', 'independent-remeasure', 'byte-identity-architecture', 'strategic-realdata']
const jobs = []
for (const t of TASK_KEYS) for (const v of ['A', 'B']) jobs.push(() => auditWithRetry(t, v))
const results = (await parallel(jobs)).filter(Boolean)

phase('Synthesize')
const SYNTH_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['close_is_correct', 'close_is_complete', 'worth_enabling', 'final_call',
             'per_task_consensus', 'faults_to_fix', 'real_data_gate'],
  properties: {
    close_is_correct: { type: 'boolean', description: 'is the positional signal implemented correctly (no sign faults / edge bugs)' },
    close_is_complete: { type: 'boolean', description: 'does it close the discovery-loss fault (robust across independent panels)' },
    worth_enabling: { type: 'string', description: 'YES / NO / GATED-ON-COMPASS with reasoning' },
    final_call: { type: 'string' },
    per_task_consensus: { type: 'array', items: { type: 'string' }, description: 'per task: do A&B agree? reconciled finding' },
    faults_to_fix: { type: 'string' },
    real_data_gate: { type: 'string', description: 'exactly what COMPASS/real-data must show before enabling' },
  },
}
const synth = await agent(
  `You are the SYNTHESIZER for a REDUNDANT 8-auditor vet of RECTIFY's positional-distinctiveness CLOSE. FOUR
tasks, TWO independent Opus-Max auditors each (signal-correctness, independent-remeasure,
byte-identity-architecture, strategic-realdata). Read ALL EIGHT durable records on disk (source of truth —
some may have stalled; only .md persists):
${TASK_KEYS.map(t => `  ${ROOT}/dev/MICROHOM_AUDIT_V5_${t}_A.md , _B.md`).join('\n')}
Returned findings: ${JSON.stringify(results)}

ROBUST CONSENSUS per task: do the two independent auditors AGREE (verdict + key_numbers)? A HOLD or real
fault from EITHER auditor blocks that task (adversarial: one solid fault is enough); unresolved A/B
disagreement = NOT-clear. Decide: (a) close_is_correct — is the signal implemented correctly (signal-
correctness consensus, no sign/edge faults); (b) close_is_complete — does it robustly close discovery-loss on
INDEPENDENT panels (independent-remeasure consensus), or is "0.4%" construction-dependent; (c) worth_enabling
— synthesize the architecture + strategic auditors: is the gate the right design, is it net-positive on a
default-off discovery tool, is it validated on anything real; (d) real_data_gate — exactly what COMPASS must
show. Be honest: "correct + complete on synthetic, but worth/real-data unproven" is a legitimate and likely
outcome. Write your synthesis to ${ROOT}/dev/MICROHOM_AUDIT_V5_SYNTHESIS.md before returning.`,
  { label: 'audit-v5:synthesis', phase: 'Synthesize', model: 'opus', effort: 'max', schema: SYNTH_SCHEMA }
)

return { results, synth }
