export const meta = {
  name: 'rectify-aligner-design-doc-crafter',
  description: 'Elite-crafter fan-out: turn the 5 ideation concepts into a phased, file:line-grounded design doc for a deliberately-orthogonal native RECTIFY aligner member. 5 concept-architects (each grounded + position-exact ablation) -> architect assembles the unified doc -> adversarial redteam pass for confabulation/overclaim/paradigm-rename.',
  phases: [
    { title: 'Draft', detail: '5 concept-architects, parallel, grounded at file:line' },
    { title: 'Architect', detail: 'assemble the unified phased design doc + risk register' },
    { title: 'Redteam', detail: 'adversarial pass: confabulation/overclaim/constraint-violation' },
  ],
}

const CONTEXT = `
RECTIFY native-aligner DESIGN-DOC crafting. Cycle scope (user-approved): produce a VALIDATED SIMULATION
BENCHMARK + this DESIGN DOC only. The native member is BUILT LATER, gated on the benchmark. So your
output is DESIGN, not code: a phased plan with file:line integration points and a position-exact
ablation per concept. EVERY primitive is a HYPOTHESIS PENDING ABLATION on the benchmark, never a finding.

=== ARCHITECTURE ALREADY DECIDED (do not relitigate) ===
- LOCALIZATION = reuse the 5 aligners' placement cluster as the realignment window (verdict:
  ADOPT_WITH_FALLBACK). Discovery ceiling is WINDOW-level, not junction-level: a read in the correct
  window but with a junction NO aligner called IS discoverable by native realignment inside that window;
  only reads with NO acceptable window (all-5-misplace / GMAP-unmapped ~12% / empty union) are out of
  reach -> that residual is the ONLY place the independent fallback localizer earns its dep, GATED behind
  a measured depletion trigger.
- REFINEMENT = a native member that does local realignment in RECTIFY's correct/consensus framework,
  NOT a 6th independent placer (it is a downstream arbitration/refinement layer, friendly to
  consensus-collapse). The 5-aligner panel is RETAINED as the default consensus fan-out; goal is to ADD
  an orthogonal member and cut external deps 5 -> 2-3 (keep minimap2), NOT replace the panel, NOT 5->0.
- META-LEVER (governs the whole design): convert RECTIFY's hard/flat/quality-blind costs to CALIBRATED
  LIKELIHOODS and emit POSTERIORS, so consensus arbitrates on a validated scale. This is BOTH the
  orthogonality source (the panel shares a flat-affine/quality-blind error family) AND the defense
  against the scoring-artifact family that already bit the project (a pure-scoring fix flipped GMAP's
  annotated win/loss 0.09 -> 1.07 with NO aligner change).

=== VERIFIED GROUNDING (you may cite; re-verify anything load-bearing) ===
- rectify/utils/local_aligner.py: live 791-line semi-global affine-gap DP. Main recurrence is FLAT
  (_MATCH=2,_MISMATCH=-4,_GAP_OPEN=-4,_GAP_EXTEND=-1 at :73-76; homo_mismatch=-2 with a BINARY homo_mask
  at :604). Already imported + called in the live 5' rescue path (read_edits.py:811,
  junction_refiner.py:717, splice_aware_5prime.py:1930/2052). THIS is the DP your refinement concepts
  extend.
- rectify/core/consensus/scoring.py: consumes NO per-base phred/quality in its objective.
- query_qualities are PROPAGATED to output but NOT consumed in scoring (chimeric_consensus.py:952,
  consensus.py:532) -> soft-decision concepts must WIRE Q in; ONT Q is often miscalibrated (recalibrate
  + ablate the recalibration separately).
- HpPenaltyTable (rectify/.../hp_penalty.py:143,179,261,270) HAS per-length del/ins costs (from_tsv, fed
  by empirical_cigar_error_profiler.py) but the gap RECURRENCE does not consult them -> real internal
  headroom. NOTE: it IS imported widely (NOT "walled in the splice refiner" - that ideation claim was
  wrong; the SPIRIT - flat recurrence - holds).
- Junction tooling to reuse: collect_junction_counts_from_bam (junction_scoring.py:473-485) ->
  per-junction {raw_count, anchored_count, read_ids}; normalize_junction / junction_ambiguity_window /
  _canonical_within_window / _normalized_annotation_set (chimeric_consensus.py:59-155); junction_validator.py
  COMPASS 3-pass; load_annotated_junctions (consensus.py:1222+).

=== HARD CONSTRAINTS (a section that violates these is rejected) ===
1. THE 87-90% REALITY: ~87-90% of contested per-segment decisions are INDEL/RESIDUE micro-arbitration
   (edit-distance bookkeeping at homopolymer/repeat/poly-A boundaries), NOT splice biology. Anchor every
   concept to which slice it targets.
2. DEPENDENCY-LIGHT IS HARD: no torch/basecaller/large weights (violates 5->2-3 + fights the AVX-512/SIGILL
   build trap; any native engine path must be AVX2/scalar, NOT AVX-512). Learned components = TINY priors
   only (KB tables, CPU). Cost every dep honestly.
3. NO HILL-CLIMB-INTO-ARTIFACT: fitness = the simulation truth set (being built), NEVER the current
   internal score. Any calibrated table needs a HELD-OUT train/test split (else the win is overfitting)
   and re-estimation per chemistry/basecaller version (a silent maintenance dep).
4. CO-PRIMARY DISCOVERY (user decision): de-novo NON-canonical junction discovery is a CO-PRIMARY
   deliverable alongside indel/CPA refinement. This RAISES FDR risk and is ANTI-CORRELATED with
   annotation-reward, so the canonical/non-canonical FDR guard must be a FIRST-CLASS design element, NOT
   an afterthought. The anti-minisplice tension (a learned GT-AG prior fights de-novo discovery) must be
   addressed head-on.
5. FRAMING METRIC for every ablation: EXACT INDEL-POSITION CONCORDANCE WITH TRUTH, not edit distance
   (edit distance is tied by construction at contested positions). Use the ambiguity-aware match so a
   call one bp into a donor/acceptor repeat is not charged FP.
6. NO PARADIGM-RENAME: localization=inherited window + gated fallback is settled; do not re-propose a
   from-scratch seeder/ chain/ pangenome / WFA-as-a-member (WFA is enabling INFRASTRUCTURE, not a concept).

Source is at /Users/kevinroy/work/rectify (branch drs-validation-rebuild). You MAY read to ground
integration points. Your final text IS the return value (structured JSON via the tool) - no preamble.`;

// The 5 concepts from the ideation synthesis (dev/ALIGNER_IDEATION_SYNTHESIS.md)
const CONCEPTS = [
  { key: 'hp-length-law', title: 'Calibrated HP-length-law emission cost in the live affine-gap recurrence',
    brief: `Replace the flat gap cost INSIDE homopolymer runs with -log P(obs_run|true_run,base) from the
already-in-tree empirical table (HpPenaltyTable). Targets the dominant HP/poly-A/STR micro-indel slice of
the 87-90%. Delta from seed-chain-extend: this is the base-DP half with the scoring function changed from
flat affine to a calibrated length-error log-likelihood (run-length-aware, where flat affine is agnostic).
Ablation: position-exact indel concordance at HP boundaries; must beat BOTH flat affine AND current
homo_mismatch=-2; held-out region split; refuted if it only ties the existing knob. Dep: light (KB table,
numba path exists). The 5-panel-convergent concept.` },
  { key: 'cpa-decoder', title: 'Channel-aware 3-prime / poly-A CPA boundary decoder (JOINT localization+refinement)',
    brief: `2-state templated-vs-tail change-point under the A-run length law to recover the true cleavage
site the panel's soft-clip heuristics drift past. JOINT (not sequential): the inherited panel END is itself
the drifted quantity, so you cannot inherit-then-refine the terminus - re-decide it jointly. RECTIFY's core
mission; fires on ~every poly-A read. Reuse the existing walkback/poly-A extractors + HP length table.
Ablation: |est-true CPA| stratified by genomic-A-tract abutment vs heuristic walkback and vs each panel end;
refuted if it does not reduce drift into genomic A-tracts. Dep: light (sequence + optional Q, no signal).` },
  { key: 'posterior-lr', title: 'Posterior-emitting refinement + likelihood-ratio consensus arbitration',
    brief: `The refiner emits, per contested decision, a POSTERIOR + a runner-up; consensus arbitrates by
calibrated log-likelihood RATIO between competing aligner placements instead of integer-score max. This is
the structural defense against the scoring-artifact family (the 0.09->1.07 integer-score flip). Not a 6th
placer - a downstream arbitration layer (aligns with the reframe), so it cannot collapse panel diversity.
Reuses local_aligner + the chimeric sync-point machinery; forward-backward over a bounded window is cheap
CPU, no weights. Ablation: on NIC/NNC-truth simulation, compare consensus accuracy under integer-score vs
calibrated-LR at the SAME panel placements; REPLICATE the 0.09->1.07 scenario and show LR does not flip;
inert if integer-score and LR pick the same winners on truth.` },
  { key: 'poa-pooled', title: 'POA-pooled per-locus consensus refinement (project-back)',
    brief: `Change the UNIT from per-read to per-locus read POOL: cluster reads at a panel-union locus bin ->
POA consensus -> align the consensus ONCE to the genome -> project placements back to each read. Denoises
stochastic per-read HP/indel error by majority BEFORE the genome decision. STRONGEST orthogonality: new
error family = mis-clustering/paralog collapse, which NO panel member exhibits (all align reads
independently). Dep: light (abPOA/spoa KB-scale C lib; reuses local_aligner for the single genome
alignment). Ablation: per-read HP-boundary placement accuracy independent vs POA-pooled-then-projected at
>=5-read loci, holding the genome DP identical, stratified by coverage; refuted if pooling does not reduce
HP-boundary error OR if mis-clustering-injected errors (test on paralog loci) exceed indel errors removed.
KEY RISK: panel-union bin IMPURITY (paralogs/distinct isoforms merged) -> needs a dispersion gate.` },
  { key: 'fracminhash-fallback', title: 'FracMinHash containment fallback localizer (panel-failure tail, gated)',
    brief: `The discovery-ceiling escape hatch: the ONLY mechanism that recovers reads with NO acceptable
panel window (all-5-misplace / GMAP-unmapped ~12% / empty union). FracMinHash containment estimates set
overlap, robust where exact-seed chains vanish at 5-10% ONT error, then feeds RECTIFY's existing DP.
Orthogonal by construction (fires only where no panel placement exists). GATE it behind a MEASURED
depletion trigger (e.g. empty/sub-threshold union OR uniformly-low post-refinement likelihood); its size
is unknown pre-benchmark, so DESIGN it but make committing the index dep conditional on the benchmark
sizing the tail. Dep: light (integer inverted sketch index, single-digit MB yeast, evictable for human;
reimplementable dependency-free). Ablation: on reads where all 5 return no/sub-threshold placement, measure
correct-localization rate vs random-window null AND measure the SIZE of this tail; refuted if panel-unplaced
reads are unplaceable by containment too (truly junk).` },
];

const INTEGRATION = { type: 'object', additionalProperties: false, required: ['file_line', 'what_changes'],
  properties: { file_line: { type: 'string' }, what_changes: { type: 'string' } } };
const PHASE = { type: 'object', additionalProperties: false, required: ['phase', 'deliverable', 'exit_criterion'],
  properties: { phase: { type: 'string' }, deliverable: { type: 'string' }, exit_criterion: { type: 'string' } } };

const SECTION_SCHEMA = { type: 'object', additionalProperties: false,
  required: ['concept', 'mechanism', 'integration_points', 'build_phases', 'ablation',
    'benchmark_requirement_consumed', 'dependency_cost', 'orthogonality', 'fdr_or_regression_guards', 'open_risks'],
  properties: {
    concept: { type: 'string' },
    mechanism: { type: 'string', description: 'Precise, buildable how-it-works grounded in the live code.' },
    integration_points: { type: 'array', minItems: 1, items: INTEGRATION, description: 'Where in the live code this attaches; re-verify file:line.' },
    build_phases: { type: 'array', minItems: 2, items: PHASE, description: 'MVP -> full, each with an exit criterion.' },
    ablation: { type: 'string', description: 'The POSITION-EXACT held-out test that proves/refutes causal contribution.' },
    benchmark_requirement_consumed: { type: 'string', description: 'Which simulation-benchmark CONTAIN/MEASURE requirement this ablation needs (so the benchmark is built to support it).' },
    dependency_cost: { type: 'string' },
    orthogonality: { type: 'string', description: 'Why its error family is independent of the 5 panel members.' },
    fdr_or_regression_guards: { type: 'string', description: 'For discovery-relevant concepts: the canonical/non-canonical FDR guard + anti-minisplice handling. For all: the regression-guard (fence-style invariant) that pins the behavior.' },
    open_risks: { type: 'array', items: { type: 'string' } },
  } };

phase('Draft');
const sections = await parallel(CONCEPTS.map(c => () =>
  agent(`${CONTEXT}\n\n=== YOUR CONCEPT TO ARCHITECT ===\n${c.title}\n${c.brief}\n\nProduce the design-doc SECTION for THIS concept only. Re-verify your integration_points at file:line in the live code. Make build_phases genuinely phased (MVP that proves the mechanism cheaply -> full). The ablation MUST be position-exact and name the benchmark requirement it consumes. If this concept touches discovery (non-canonical junctions), the FDR guard is mandatory and first-class.`,
    { label: `draft:${c.key}`, phase: 'Draft', schema: SECTION_SCHEMA })
));
const drafted = sections.filter(Boolean);
log(`Draft: ${drafted.length}/5 concept sections grounded`);

phase('Architect');
const ARCHITECT_SCHEMA = { type: 'object', additionalProperties: false,
  required: ['design_doc_markdown', 'architecture_summary', 'phased_roadmap', 'benchmark_coupling',
    'risk_register', 'regression_guard_plan', 'discovery_fdr_strategy', 'open_decisions_for_user'],
  properties: {
    design_doc_markdown: { type: 'string', description: 'The FULL design doc as markdown: architecture, the 5 concept sections woven into a coherent member design, the phased roadmap, benchmark coupling, regression-guard plan, risk register, discovery/FDR strategy. This is the deliverable.' },
    architecture_summary: { type: 'string', description: 'The member in 1 paragraph: inherited-window localization + gated fallback; calibrated/posterior refinement DP; LR consensus arbitration; deps 5->2-3.' },
    phased_roadmap: { type: 'array', items: PHASE, description: 'Cross-concept ordered roadmap (what to build first: benchmark -> which concept MVP -> ...). Sequence by leverage x risk.' },
    benchmark_coupling: { type: 'string', description: 'How each concept ablation maps onto the simulation benchmark CONTAIN/MEASURE requirements; what the benchmark must include for the doc to be testable.' },
    risk_register: { type: 'array', items: { type: 'object', additionalProperties: false, required: ['risk', 'severity', 'mitigation'], properties: { risk: { type: 'string' }, severity: { type: 'string', enum: ['high', 'medium', 'low'] }, mitigation: { type: 'string' } } } },
    regression_guard_plan: { type: 'string', description: 'How to extend the test_gmap_fence_regression.py fence-invariant discipline to pin each new behavior so improvements cannot silently regress (and so a calibration table cannot be overfit into the guards).' },
    discovery_fdr_strategy: { type: 'string', description: 'The first-class canonical/non-canonical FDR strategy for CO-PRIMARY de-novo discovery; how the anti-minisplice tension is resolved (de-novo discovery vs learned GT-AG prior).' },
    open_decisions_for_user: { type: 'array', items: { type: 'string' } },
  } };

const architect = await agent(
  `${CONTEXT}\n\n=== YOUR ROLE: ARCHITECT ===\nYou received ${drafted.length} grounded concept sections (below). Assemble them into ONE coherent, phased design doc for the native member. Resolve cross-concept interactions (e.g. the HP-length-law DP feeds the posterior; the CPA decoder and the DP share the length table; POA-pooling sits UPSTREAM of the per-read DP; the fallback localizer feeds the same DP). Sequence the roadmap by leverage x risk, with the simulation benchmark as the gating first step. Make the CO-PRIMARY discovery FDR strategy first-class (canonical/non-canonical guard + anti-minisplice resolution). Specify the regression-guard plan that extends the existing fence-test discipline. Write design_doc_markdown as the full deliverable doc. Keep every claim a hypothesis-pending-ablation; cite file:line for integration points.\n\n=== GROUNDED SECTIONS ===\n${JSON.stringify(drafted, null, 1)}`,
  { label: 'architect', phase: 'Architect', schema: ARCHITECT_SCHEMA, effort: 'high' });

phase('Redteam');
const REDTEAM_SCHEMA = { type: 'object', additionalProperties: false,
  required: ['overall_verdict', 'punch_list', 'unsupported_file_lines', 'must_fix_before_build'],
  properties: {
    overall_verdict: { type: 'string', enum: ['sound', 'sound_with_fixes', 'needs_rework'] },
    punch_list: { type: 'array', items: { type: 'object', additionalProperties: false, required: ['severity', 'category', 'location', 'issue', 'fix'],
      properties: {
        severity: { type: 'string', enum: ['blocker', 'major', 'minor'] },
        category: { type: 'string', enum: ['confabulation', 'overclaim_vs_prior_art', 'constraint_violation', 'paradigm_rename', 'unsupported_file_line', 'untestable_ablation', 'gap', 'fdr_risk'] },
        location: { type: 'string' }, issue: { type: 'string' }, fix: { type: 'string' } } } },
    unsupported_file_lines: { type: 'array', items: { type: 'string' }, description: 'Any file:line cited in the doc that you checked and found wrong/nonexistent/mischaracterized.' },
    must_fix_before_build: { type: 'array', items: { type: 'string' }, description: 'The blocker subset, distilled.' },
  } };

const redteam = await agent(
  `${CONTEXT}\n\n=== YOUR ROLE: ADVERSARIAL REDTEAM ===\nThe design doc below was assembled by the architect. Attack it. SPOT-CHECK its cited file:line claims against the live source (sample the load-bearing ones; report any that are wrong/mischaracterized in unsupported_file_lines). Hunt for: confabulation (authoritative claims with no grounding), novelty overclaim vs prior art (pair-HMM/Dindel/GATK for the DP; sourmash/mash for containment; FLAIR/StringTie/POA for pooling), constraint violations (deps, AVX-512, hill-climb-into-artifact, the 87-90% reality), paradigm-renames (seed-chain-extend reskinned as novel), untestable/edit-distance ablations (must be position-exact), and FDR risk in the co-primary discovery path. Be specific and cite locations in the doc. Default to skeptical: if a claim cannot be tied to file:line + a provable read-class win, flag it.\n\n=== DESIGN DOC ===\n${architect.design_doc_markdown}`,
  { label: 'redteam', phase: 'Redteam', schema: REDTEAM_SCHEMA, effort: 'high' });

return { sections: drafted, architect, redteam };
