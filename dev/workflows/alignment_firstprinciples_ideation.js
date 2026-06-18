export const meta = {
  name: 'rectify-alignment-firstprinciples-ideation',
  description: '6-lens first-principles creativity panel reimagining ONT DRS/cDNA read alignment as localization + refinement, then a synthesis judge that clusters/filters against the hard-won constraints and emits a ranked concept shortlist + benchmark_requirements for the crafter design doc.',
  phases: [
    { title: 'Ideate', detail: '6 divergent creative lenses, parallel, strict per-concept schema' },
    { title: 'Synthesize', detail: 'cluster, reject paradigm-rename/dep-violations, rank, emit benchmark_requirements' },
  ],
}

// ── Shared brief: the hard-won ground truth every lens must internalize ────────
const CONTEXT = `
YOU ARE ON A FIRST-PRINCIPLES IDEATION PANEL for RECTIFY, a long-read RNA aligner
CONSENSUS + 3'-end-correction tool (ONT direct-RNA "DRS" and PCR-cDNA reads, yeast
+ human). The PI wants to reimagine how reads get aligned, from scratch, possibly as
a two-stage pipeline: a read "LOCALIZATION" stage (find where on the genome a read
belongs) then an alignment "REFINEMENT" stage (settle the exact base-level / junction
placement). You may ADOPT, EXTEND, or CHALLENGE that decomposition.

THIS IS DIVERGENT IDEATION. Novelty is the goal. But this panel runs BEFORE a
validated benchmark exists, so EVERY concept you propose is a HYPOTHESIS PENDING
ABLATION, never a finding. Label it as such. Mechanism-reading (what a technique DOES)
is NOT evidence of causal contribution (what it WINS).

=== THE EXISTING 5-ALIGNER PANEL (what you must be ORTHOGONAL to) ===
RECTIFY today runs 5 external aligners and forms a consensus:
  - minimap2   : minimizer seed -> chain -> base-level affine DP w/ splice kernel. 1.1% non-canonical. The incumbent.
  - uLTRA      : MEM + annotation-guided 2-pass. Strong WITH a GTF; snaps to catalog.
  - deSALT     : de-Bruijn-graph 2-pass. (Operationally fragile: SIGSEGV/OOM, Linux-only.)
  - mapPacBio  : k-mer + affine DP (BBMap family). 59.8% non-canonical (noisy).
  - GMAP       : sampled-oligomer hash -> diagonalization -> own genomic splice DP.
                 55.6% non-canonical BUT uniquely finds ~127 recurrent GT-AG novel junctions.
NOTE: minimap2 / winnowmap2 / graphmap2 / gapmm2 ALL share minimap2's splice DP kernel
(noncan=9). Real splice-model diversity exists only in uLTRA / deSALT / GMAP / mapPacBio.

=== HARD-WON CONSTRAINTS (violating these makes a concept dead-on-arrival) ===
1. THE 87-90% REALITY (anchor everything to this): when all aligners are made to
   disagree, ~87-90% of contested per-segment decisions are INDEL / RESIDUE
   MICRO-ARBITRATION (edit-distance bookkeeping near homopolymers/repeats), NOT splice
   biology. Most "where do reads actually get decided" lives in micro-indel placement,
   poly-A/homopolymer boundaries, and 3'-end (CPA) walkback — not exotic splice models.
   A concept that only improves splice-junction calling addresses <15% of decisions.
   State explicitly which slice of this reality your concept targets.
2. CONSENSUS IS LOAD-BEARING ON ~100% OF READS. Cross-aligner selection decides every
   read (67-71% winner-take-all fallback + 22-25% chimeric stitch). So a single
   "best-of-all" aligner is one CORRELATED error mode that cannot consensus-correct
   against itself. Any native member must be DELIBERATELY ORTHOGONAL by design (independent
   error modes), NOT a distillation of the 5 (a distillation is maximally correlated = the
   single worst possible new consensus member). The goal is to ADD an orthogonal member and
   cut external deps 5 -> 2-3 (keep minimap2), NOT to replace the panel and NOT 5 -> 0.
3. DEPENDENCY-LIGHT IS A HARD TARGET. Anything that drags torch / a basecaller / large
   model weights VIOLATES the 5->2-3 dep goal and the stack already fights an AVX-512/SIGILL
   build trap. Learned components are allowed ONLY as TINY priors (kilobytes, CPU, no
   framework). If your idea needs a GPU or a 100MB model, say so and cost it honestly.
4. NO HILL-CLIMB-INTO-ARTIFACT. The scoring criteria were provably artifacts 6 commits ago
   (a pure-scoring fix with NO aligner change flipped GMAP's annotated win/loss 0.09 -> 1.07).
   Fitness MUST be a validated benchmark (being built now: simulation with known NIC/NNC
   junction truth), not the current internal score. The annotation-reward and de-novo novel-
   junction-DISCOVERY missions are ANTI-CORRELATED; do not optimize one away.

=== THE SEED-CHAIN-EXTEND TRAP (read carefully) ===
The "localization then refinement" decomposition the PI seeded IS LITERALLY minimap2's
paradigm: localization = seed + chain, refinement = base-level DP. If you propose a
localization or refinement concept, you MUST state its EXPLICIT DELTA from minimap2
seed-chain-extend. "Seed differently then DP" reskinned is NOT a new concept — it is the
incumbent renamed. Be honest: if your idea IS seed-chain-extend with a tweak, say so and
justify the tweak; do not dress it as a new paradigm.

=== THE PI'S LEADING HYPOTHESIS FOR THE LOCALIZATION STAGE (engage with this directly) ===
The PI proposes that the LOCALIZATION stage may not need a new seeder/indexer AT ALL:
instead, REUSE the loci where the existing 5 aligners have ALREADY placed each read (the
union/cluster of their placements), and have the new RECTIFY aligner do a LOCAL REALIGNMENT
within that bounded window, grounded in RECTIFY's existing correct/consensus framework
(local_aligner.py affine-gap DP + chimeric scoring + poly-A walkback). Implications you must
weigh:
  - This SIDESTEPS the seed-chain-extend trap entirely: there is no seeding to reinvent, so
    the ENTIRE novelty budget moves to REFINEMENT. localization becomes "free" (inherited).
  - It is maximally dependency-light (no new index) and tractable now (local_aligner.py exists).
  - It reframes the native component as a REFINEMENT/ARBITRATION layer over the diverse panel's
    placements, NOT a 6th independent placer - which is friendlier to the consensus-collapse
    concern (it operates DOWNSTREAM of placement diversity, not competing with it).
  - BUT interrogate the failure modes it inherits: reads ALL 5 misplace or fail to place (GMAP
    maps 12% FEWER reads); loci where the panel's union is empty or wrong; whether local
    realignment in a panel-seeded window can DISCOVER novel junctions the panel never proposed,
    or only ARBITRATE among the ones it did. Does inheriting localization cap the discovery
    ceiling? When is a genuinely independent localization still needed as a fallback?
For your lens: treat this as the DEFAULT localization hypothesis. Either strengthen/extend it,
or make the explicit case for where a novel localization stage still earns its keep despite it.

=== GROUNDING ===
Source is at /Users/kevinroy/work/rectify. You MAY read to ground ORTHOGONALITY claims
(rectify/core/consensus/chimeric_consensus.py, rectify/utils/local_aligner.py [already a
native 791-line affine-gap DP, live in the 5' rescue path], docs/aligners/*.md). BUDGET
your reads: ideation is the goal, not an audit. Cite file:line only where it backs a claim.
Your final text IS the return value (structured JSON via the tool) - no preamble.`;

// ── 6 deliberately non-overlapping creative lenses ────────────────────────────
const LENSES = [
  {
    key: 'info-theory',
    title: 'Information-theory / noisy-channel decoding',
    prompt: `LENS: treat read alignment as DECODING a message sent through the ONT noisy
channel. Localization = coarse channel/position estimation; refinement = MAP/ML decoding
under an explicit ONT error model (homopolymer length errors, basecaller-correlated
substitutions, signal-level ambiguity). Think: posterior over placements, soft-decision
decoding, sequential probability ratio tests for early localization, error-model-aware
edit costs that beat flat affine-gap on the 87-90% indel reality. What does an explicit
generative ONT error model buy that minimap2's flat scoring leaves on the table?`,
  },
  {
    key: 'biology-first',
    title: 'Biology-first (splice / CPA / poly-A / RNA structure)',
    prompt: `LENS: exploit what the MOLECULE tells us that a generic DNA-string aligner
ignores. Co-transcriptional splicing, branch-point/poly-pyrimidine context, CPA + poly-A
walkback (RECTIFY's core mission), secondary-structure-induced basecaller errors, the
donor/acceptor ambiguity window. CAUTION: the 87-90% reality says splice-only ideas hit
<15% of decisions - so you must connect biology to the MICRO-ARBITRATION majority (e.g.
poly-A/homopolymer boundary placement, internal-priming, 3'-end truth) to earn weight.
Where does biological prior beat string-generic DP on the decisions that actually dominate?`,
  },
  {
    key: 'learned-tiny',
    title: 'Learned-representation / ML (dependency-light ONLY)',
    prompt: `LENS: learned representations for localization and/or refinement, under a HARD
constraint - NO torch, NO basecaller, NO large weights (violates the 5->2-3 dep target).
You are limited to TINY learned priors: kilobyte lookup tables, learned k-mer weights /
hashing, small linear/logistic models, learned edit-cost matrices, distilled motif priors
(flag the minisplice-style GT-AG motif-bias conflict: a learned splice prior fights the
de-novo discovery mission). Be brutally honest about dependency weight - if it needs a
framework, mark it hard-constraint-violating. What is the MOST a kilobyte-scale learned
prior can do for localization or micro-indel arbitration?`,
  },
  {
    key: 'classical-index',
    title: 'Classical algorithms & novel index/data structures',
    prompt: `LENS: reinvent the seeding/indexing/localization machinery with modern
sequence data structures - strobemers, syncmers, open/closed syncmers, minimizer variants,
MinHash/FracMinHash sketching, FM-index/r-index, wavelet trees, sparse suffix arrays,
spaced/gapped seeds. Target sublinear or more-robust localization, especially for the
error profile of ONT (where exact minimizers break). For refinement, consider banded /
wavefront (WFA) / bit-parallel Myers edit distance tuned to the homopolymer reality.
REQUIRED: state the explicit delta from minimap2 seed-chain-extend for every concept.`,
  },
  {
    key: 'cross-domain',
    title: 'Cross-domain analogy (SLAM / DTW / vision / ECC)',
    prompt: `LENS: borrow primitives from OTHER fields that solve "align a noisy observation
to a reference." Robotics SLAM / GPS (coarse localize then refine pose; loop closure;
particle filters); audio/speech Dynamic Time Warping & subsequence-DTW; computer-vision
image registration (coarse-to-fine pyramids, RANSAC geometric verification, feature
matching); error-correcting codes (belief propagation, Viterbi/trellis decoding);
compression (BWT, context modeling). For each borrowed primitive, map it concretely onto
DRS/cDNA localization or refinement, and name what it does that current aligners don't.
Cross-domain != hand-wave: give the concrete mechanism and the read class it targets.`,
  },
  {
    key: 'contrarian',
    title: 'Adversarial-contrarian (challenge the premise itself)',
    prompt: `LENS: attack the framing. Maybe "align each read independently to a linear
genome" is the wrong primitive. Explore: per-locus read-vs-read assembly/consensus FIRST
then align the consensus; pangenome/variation-graph or personalized-reference alignment;
alignment-FREE localization (k-mer profile -> locus) with alignment only where needed;
clustering reads into isoform groups before placement; co-estimating the transcriptome and
the alignments jointly (EM). Also challenge the localization->refinement split itself: is a
two-stage pipeline even right, or should localization and refinement be JOINT? Your job is
to find the idea that makes the other 5 lenses look conventional - but it must still respect
the 4 hard constraints (esp. dependency-light and orthogonality).`,
  },
];

const CONCEPT = {
  type: 'object',
  additionalProperties: false,
  required: ['name', 'stage', 'one_line', 'mechanism', 'target_failure_mode',
    'delta_from_seed_chain_extend', 'orthogonality_to_panel', 'dependency_cost',
    'tractability', 'ablation', 'novelty_honest'],
  properties: {
    name: { type: 'string' },
    stage: { type: 'string', enum: ['localization', 'refinement', 'both', 'reframe'] },
    one_line: { type: 'string' },
    mechanism: { type: 'string', description: 'Concrete how-it-works, not a slogan.' },
    target_failure_mode: { type: 'string', description: 'Which slice of the 87-90% indel/residue reality OR which named read class this targets, and roughly what fraction of decisions it touches.' },
    delta_from_seed_chain_extend: { type: 'string', description: 'Explicit difference from minimap2 seed+chain+base-DP. If it IS seed-chain-extend with a tweak, say so and justify the tweak.' },
    orthogonality_to_panel: { type: 'string', description: 'Why its error modes are INDEPENDENT of minimap2/uLTRA/deSALT/mapPacBio/GMAP. Name the families.' },
    dependency_cost: { type: 'string', enum: ['light', 'medium', 'heavy'], },
    dependency_detail: { type: 'string', description: 'What it drags in. For learned ideas, the honest framework/weight cost.' },
    tractability: { type: 'string', enum: ['buildable_now', 'research', 'multi_year'] },
    ablation: { type: 'string', description: 'The held-out test on validated truth that would prove or refute the causal contribution.' },
    novelty_honest: { type: 'string', enum: ['novel', 'repackaged_prior_art', 'uncertain'] },
    prior_art: { type: 'string', description: 'Named prior work if known; "none-known" otherwise.' },
  },
};

const CONCEPTS_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['lens', 'concepts', 'wildcard', 'challenge_to_decomposition'],
  properties: {
    lens: { type: 'string' },
    concepts: { type: 'array', minItems: 3, maxItems: 6, items: CONCEPT },
    wildcard: { ...CONCEPT, description: 'One deliberately risky/high-variance idea.' },
    challenge_to_decomposition: { type: 'string', description: 'Does localization->refinement hold for your lens? If not, what instead (joint? more stages? different primitive)?' },
  },
};

phase('Ideate');
const panels = await parallel(LENSES.map(l => () =>
  agent(`${CONTEXT}\n\n=== YOUR ASSIGNED LENS: ${l.title} ===\n${l.prompt}\n\nProduce 3-6 grounded concepts + 1 wildcard + your challenge to the decomposition. Every field required. Anchor target_failure_mode to the 87-90% reality. Fill delta_from_seed_chain_extend honestly.`,
    { label: `ideate:${l.key}`, phase: 'Ideate', schema: CONCEPTS_SCHEMA })
));

const good = panels.filter(Boolean);
log(`Ideate done: ${good.length}/6 lenses returned; ${good.reduce((n, p) => n + (p.concepts ? p.concepts.length : 0) + 1, 0)} concepts total`);

// ── Synthesis: cluster, reject, rank, and WIRE INTO THE BENCHMARK ─────────────
const SYNTH_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['pi_localization_hypothesis_verdict', 'clusters', 'localization_shortlist',
    'refinement_shortlist', 'reframes', 'cross_cutting_principles', 'top_concepts_for_design_doc',
    'rejected', 'benchmark_requirements', 'confabulation_flags', 'open_questions_for_user'],
  properties: {
    pi_localization_hypothesis_verdict: { type: 'object', additionalProperties: false,
      required: ['stance', 'rationale', 'discovery_ceiling', 'fallback_needed'],
      properties: {
        stance: { type: 'string', enum: ['adopt', 'adopt_with_fallback', 'reject'], description: 'On reusing panel loci for localization + native local realignment.' },
        rationale: { type: 'string' },
        discovery_ceiling: { type: 'string', description: 'Whether panel-seeded windows can DISCOVER novel junctions the panel never proposed, or only ARBITRATE - and what that caps.' },
        fallback_needed: { type: 'string', description: 'When (if ever) a genuinely independent localization is still needed (reads all 5 misplace / GMAP 12% unmapped / empty union).' },
      } },
    clusters: { type: 'array', items: { type: 'object', additionalProperties: false,
      required: ['theme', 'stage', 'member_concepts', 'verdict'],
      properties: {
        theme: { type: 'string' },
        stage: { type: 'string', enum: ['localization', 'refinement', 'both', 'reframe'] },
        member_concepts: { type: 'array', items: { type: 'string' } },
        verdict: { type: 'string' },
      } } },
    localization_shortlist: { type: 'array', description: 'Ranked best localization concepts.',
      items: { type: 'object', additionalProperties: false, required: ['concept', 'why', 'delta_from_sce', 'tractability'],
        properties: { concept: { type: 'string' }, why: { type: 'string' }, delta_from_sce: { type: 'string' }, tractability: { type: 'string' } } } },
    refinement_shortlist: { type: 'array', description: 'Ranked best refinement concepts.',
      items: { type: 'object', additionalProperties: false, required: ['concept', 'why', 'delta_from_sce', 'tractability'],
        properties: { concept: { type: 'string' }, why: { type: 'string' }, delta_from_sce: { type: 'string' }, tractability: { type: 'string' } } } },
    reframes: { type: 'array', description: 'Ideas that credibly challenge the 2-stage split (joint, assembly-first, alignment-free, etc.).', items: { type: 'string' } },
    cross_cutting_principles: { type: 'array', description: 'Design principles that recur across lenses and should govern the member regardless of which concepts win.', items: { type: 'string' } },
    top_concepts_for_design_doc: { type: 'array', description: 'The ranked set to carry into the crafter design doc. Each must clear orthogonality + dependency-light + an ablation.',
      items: { type: 'object', additionalProperties: false, required: ['concept', 'stage', 'rationale', 'orthogonality', 'dependency_cost', 'ablation', 'novelty_honest'],
        properties: { concept: { type: 'string' }, stage: { type: 'string' }, rationale: { type: 'string' }, orthogonality: { type: 'string' }, dependency_cost: { type: 'string' }, ablation: { type: 'string' }, novelty_honest: { type: 'string' } } } },
    rejected: { type: 'array', description: 'Concepts cut, with the reject category.',
      items: { type: 'object', additionalProperties: false, required: ['concept', 'reason'],
        properties: { concept: { type: 'string' }, reason: { type: 'string', enum: ['standard_paradigm_renamed', 'hard_constraint_violation_deps', 'not_orthogonal', 'untractable', 'confabulated', 'addresses_wrong_slice'] } } } },
    benchmark_requirements: { type: 'array', description: 'WHAT the simulation truth set must CONTAIN and MEASURE to discriminate these concepts (failure modes to simulate, metrics to score). This is how the panel feeds the gating deliverable.', items: { type: 'string' } },
    confabulation_flags: { type: 'array', description: 'Concepts that SOUND authoritative but lack file:line / read-class grounding or claim novelty for a known technique.', items: { type: 'string' } },
    open_questions_for_user: { type: 'array', items: { type: 'string' } },
  },
};

phase('Synthesize');
const synthesis = await agent(
  `${CONTEXT}\n\n=== YOUR ROLE: SYNTHESIS JUDGE ===\nYou received ${good.length} lens panels of concepts (below). Cluster them, then apply the hard-won filters RUTHLESSLY:\n` +
  `- REJECT "standard_paradigm_renamed": any localization/refinement idea that is minimap2 seed-chain-extend with a cosmetic tweak (check delta_from_seed_chain_extend honestly - this is a DISTINCT check from orthogonality; a concept can be orthogonal in error-mode yet still be the same paradigm).\n` +
  `- REJECT "hard_constraint_violation_deps": anything dragging torch / basecaller / large weights (violates 5->2-3).\n` +
  `- REJECT "not_orthogonal": distillations of the 5, or ideas correlated with an existing family's error mode.\n` +
  `- REJECT "addresses_wrong_slice": splice-only ideas that miss the 87-90% indel/residue majority without justification.\n` +
  `- FLAG confabulation: authoritative-sounding claims with no file:line / read-class grounding, or novelty claimed for known prior art.\n` +
  `FIRST deliver pi_localization_hypothesis_verdict: adopt / adopt_with_fallback / reject the PI's proposal to REUSE panel loci for localization + native local realignment - and be specific about the DISCOVERY CEILING (can panel-seeded windows find novel junctions the panel never proposed, or only arbitrate?) and when an independent localization fallback is still needed.\n` +
  `Then RANK a localization shortlist + a refinement shortlist + surface credible reframes. Produce top_concepts_for_design_doc (each must clear orthogonality + dependency-light + a concrete ablation).\n` +
  `CRITICAL DELIVERABLE - benchmark_requirements: state what the simulation-with-known-truth set must CONTAIN (which failure modes to simulate: homopolymer-boundary indels, internal-priming poly-A, NIC/NNC novel junctions, etc.) and MEASURE (which metrics discriminate these concepts) so that this idea-pile becomes EVALUABLE the moment the benchmark lands. Every concept stays a hypothesis-pending-ablation.\n\n` +
  `=== LENS PANELS ===\n${JSON.stringify(good, null, 1)}`,
  { label: 'synthesis', phase: 'Synthesize', schema: SYNTH_SCHEMA, effort: 'high' });

return { panels: good, synthesis };
