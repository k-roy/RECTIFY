export const meta = {
  name: 'rectify-compass-shortread-junction-scoping',
  description: 'Scope wrapping COMPASS-style multi-aligner SHORT-READ novel-junction discovery into RECTIFY (bbmap + minimap2 + a third, NOT STAR-alone), to robustly validate/discover novel junctions as the short-read complement to the long-read consensus. 4 grounded scoping agents -> synthesis.',
  phases: [
    { title: 'Scope', detail: '4 parallel: COMPASS method, short-read aligner panel, RECTIFY integration, application/re-validation' },
    { title: 'Synthesize', detail: 'integration plan + immediate experiment + open questions' },
  ],
}

const CONTEXT = `
RECTIFY is a long-read (ONT DRS/cDNA) RNA aligner-CONSENSUS + 3'-end-correction tool (Chanfreau lab).
This effort scopes adding SHORT-READ (Illumina) novel-junction discovery to RECTIFY via a COMPASS-style
MULTI-ALIGNER consensus, as the orthogonal complement to the long-read panel.

=== WHY (the trigger) ===
We have ~111 GMAP-only recurrent GT-AG "novel" junctions on A549 chr5 from the long-read panel. We tried to
validate them with SG-NEx A549 Illumina short reads, but used the SG-NEx-provided STAR 2.6.0c SINGLE-PASS
BAMs. The PI's point (well-documented): **STAR is poor for NOVEL junctions (high FDR; 2-pass can make it
worse)** — so a STAR-alone validation is not trustworthy. A STAR-INDEPENDENT coverage test (inside/flank
coverage ratio >0.5 for all 111) suggests the 111 are likely artifacts, BUT the split-read validation must
be REDONE with a robust MULTI-ALIGNER short-read consensus. The PI wants this built INTO RECTIFY.

=== WHAT EXISTS IN RECTIFY (verify at file:line; /Users/kevinroy/work/rectify, branch drs-validation-rebuild) ===
- **COMPASS is already partially present**: rectify/core/splice/junction_validator.py implements a "COMPASS
  3-pass" cross-sample junction pipeline: extract_sample_junctions() -> filter_cross_sample_junctions()
  (min-samples/min-reads/max-intron/canonicality) -> apply_junction_filter(). READ IT — understand what
  COMPASS does today and what's missing for short-read MULTI-ALIGNER discovery.
- Junction extraction + concordance: collect_junction_counts_from_bam (junction_scoring.py:473) ->
  per-junction {raw_count, anchored_count}; build_junction_pool (junction_scoring.py:489) with cross-FAMILY
  concordance relaxation (ALIGNER_FAMILY map); the ambiguity-aware helpers normalize_junction /
  junction_ambiguity_window / _canonical_within_window / _normalized_annotation_set
  (chimeric_consensus.py:59-155); load_annotated_junctions (consensus.py:1222).
- The long-read multi-aligner consensus (multi_aligner.py, chimeric_consensus.py) is the architectural
  precedent — the short-read COMPASS should mirror its cross-aligner-concordance philosophy.

=== DATA AVAILABLE (Sherlock) ===
- SG-NEx A549 Illumina: 3 PE reps, genome BAMs already pulled to
  /scratch/users/kevinroy/sgnex_a549_illumina/replicate{1,3,5}/ (STAR 1-pass). RAW FASTQ on public S3
  s3://sg-nex-data/data/sequencing_data_illumina/fastq/SGNex_A549_Illumina_replicate{1,3,5}_run1/ (~25GB)
  — re-alignable with other aligners. Reference: GRCh38 (Ensembl '5' naming; our long-read junctions are
  'chr5' — harmonize). PacBio A549 (1 SMRTcell) also available.
- The 111 candidates: rectify/dev/gmap_only_recurrent_novels_chr5.tsv; the 609 corroborated + the
  non-canonical set are in /scratch/users/kevinroy/deliverable_b/gmap_corroboration.json.
- Sherlock: conda envs include rectify, aligner_bench, base. aws CLI available.

=== HARD REQUIREMENTS ===
- GROUND every claim: file:line for RECTIFY code; verify aligner availability with which/conda on Sherlock
  (do NOT assume a tool is installed); cite sources for COMPASS / aligner novel-junction behavior. Flag
  unknowns as unknowns. NO confabulation — this exercise has already been burned by ungrounded claims.
- Dependency-light + reproducible (the lab runs chunked SLURM on Sherlock owners partition w/ an AVX-512
  constraint; short-read alignment of 25GB FASTQ is a real compute job — scope it as chunked).
- Short-read SPLICED alignment caveats: minimap2 '-ax sr' is NOT spliced-aware; assess what each candidate
  aligner actually does for novel spliced junctions. STAR may still be ONE panel member but NOT the sole one.
Your final text IS the return value (structured JSON via the tool). Read before you write.`;

const TASKS = [
  { key: 'compass-method', title: 'COMPASS method — what it is + what RECTIFY already implements',
    prompt: `Define precisely what a COMPASS-style multi-aligner junction-consensus does. (1) Research
"COMPASS" splice-junction / novel-junction discovery (web: is there a published COMPASS method? the
Chanfreau/Pleiss/Guttman lineage? what is its core algorithm?). (2) READ rectify/core/splice/junction_validator.py
and document the existing "COMPASS 3-pass" at file:line — exactly what it does (extract -> cross-sample
filter -> apply), its thresholds, and its concordance logic. (3) Define the DELTA needed for SHORT-READ
MULTI-ALIGNER novel-junction discovery: cross-ALIGNER (not just cross-sample) concordance, how to combine
per-aligner junction sets robustly, how novel (non-annotated) junctions are admitted with controlled FDR.
Output the algorithm RECTIFY should implement, grounded in what already exists.` },
  { key: 'shortread-panel', title: 'Short-read spliced-aligner panel (bbmap + minimap2 + a third)',
    prompt: `Scope the MULTI-ALIGNER short-read SPLICED panel for orthogonal novel-junction discovery. For
EACH candidate, state precisely how it handles NOVEL spliced junctions and its known novel-junction FDR
behavior, with sources: bbmap (which mode/flags — maxindel/intronlen; does it emit N-ops?), minimap2 (note
'-ax sr' is NOT spliced — is there any spliced short-read mode? if not, say so), and candidate THIRDS
(HISAT2, subread/subjunc, MapSplice2, segemehl, gsnap, olego, etc.). The PI says STAR is poor for novel
junctions (2-pass can worsen) — assess and decide whether STAR stays as ONE member or is dropped. Recommend
the orthogonal panel (which give genuinely independent novel-junction error modes) + EXACT invocations.
VERIFY availability on Sherlock (which/conda across rectify, aligner_bench, base; or bioconda installability).
Output the recommended panel with per-aligner rationale + commands + an availability table.` },
  { key: 'rectify-integration', title: 'RECTIFY integration architecture',
    prompt: `Scope how to wrap short-read multi-aligner COMPASS junction discovery INTO RECTIFY. Ground at
file:line. Cover: where it attaches (extend junction_validator.py? a new module rectify/core/splice/
shortread_compass.py? a new CLI subcommand e.g. 'rectify junctions discover'?); reuse of
collect_junction_counts_from_bam (junction_scoring.py:473), normalize_junction + the ambiguity helpers
(chimeric_consensus.py:59-155), build_junction_pool's cross-family concordance (junction_scoring.py:489),
load_annotated_junctions; the data flow (short-read FASTQ -> N aligners -> per-aligner junction sets ->
COMPASS concordance -> novel-junction set with support); and crucially HOW the short-read novel-junction set
CROSS-VALIDATES the long-read consensus novels (the 111/87/609) — the join logic incl. chrom-name (chr5<->5)
+ ambiguity normalization. Output the module/CLI design + the integration points + the long-read<->short-read
join design.` },
  { key: 'application', title: 'Immediate re-validation experiment + benchmark feed',
    prompt: `Design the IMMEDIATE experiment: properly re-validate the 111 (and the 87 "validated" + 609
corroborated) on SG-NEx A549 using the multi-aligner short-read COMPASS approach instead of STAR-alone.
Concretely: re-align the 3 A549 Illumina reps (FASTQ on S3, ~25GB) with the recommended panel on Sherlock
(chunked, owners, AVX-512 where rectify-env is used); run COMPASS; recompute validation with proper controls
(positive = annotated junctions in the same loci; negative = the non-canonical set). Specify thresholds,
the chrom-name harmonization, the anchor/recurrence gates, and what result would CONFIRM vs REFUTE the 111
as artifacts (and re-assess the keep-GMAP decision). Also: how this short-read junction truth feeds the P0
simulation benchmark (the NIC/NNC truth + the FDR controls) and the design-doc discovery-FDR track. Output
the experiment plan (steps, jobs, thresholds, decision criteria) + benchmark/design-doc implications.` },
];

const SCOPE_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['area', 'findings', 'grounding', 'recommendation', 'open_questions', 'unknowns'],
  properties: {
    area: { type: 'string' },
    findings: { type: 'array', items: { type: 'string' }, description: 'Key grounded findings.' },
    grounding: { type: 'array', items: { type: 'string' }, description: 'file:line refs, verified tool availability, cited sources backing the findings.' },
    recommendation: { type: 'string', description: 'The concrete recommended design/plan for this area.' },
    open_questions: { type: 'array', items: { type: 'string' } },
    unknowns: { type: 'array', items: { type: 'string' }, description: 'Things NOT verified / need checking — flagged honestly, not glossed.' },
  },
};

phase('Scope');
const scopes = await parallel(TASKS.map(t => () =>
  agent(`${CONTEXT}\n\n=== YOUR SCOPING AREA: ${t.title} ===\n${t.prompt}`,
    { label: `scope:${t.key}`, phase: 'Scope', schema: SCOPE_SCHEMA })
));
const good = scopes.filter(Boolean);
log(`Scope: ${good.length}/4 areas returned`);

phase('Synthesize');
const SYNTH_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['compass_shortread_design', 'recommended_aligner_panel', 'rectify_integration_plan',
    'immediate_experiment', 'benchmark_and_designdoc_impact', 'phased_plan', 'open_questions_for_user', 'risks_and_unknowns'],
  properties: {
    compass_shortread_design: { type: 'string', description: 'The COMPASS multi-aligner short-read junction-discovery algorithm RECTIFY should implement (grounded in the existing junction_validator.py).' },
    recommended_aligner_panel: { type: 'string', description: 'The chosen short-read panel (bbmap + minimap2 + the third), per-aligner rationale, STAR keep/drop, exact invocations, Sherlock availability.' },
    rectify_integration_plan: { type: 'string', description: 'Module/CLI design + integration points (file:line) + the long-read<->short-read cross-validation join.' },
    immediate_experiment: { type: 'string', description: 'The re-validation experiment (jobs, thresholds, controls, decision criteria for confirming/refuting the 111 + the keep-GMAP call).' },
    benchmark_and_designdoc_impact: { type: 'string', description: 'How short-read junction truth feeds the P0 benchmark + the discovery-FDR design track.' },
    phased_plan: { type: 'array', items: { type: 'object', additionalProperties: false, required: ['phase', 'deliverable', 'exit_criterion'], properties: { phase: { type: 'string' }, deliverable: { type: 'string' }, exit_criterion: { type: 'string' } } } },
    open_questions_for_user: { type: 'array', items: { type: 'string' } },
    risks_and_unknowns: { type: 'array', items: { type: 'string' }, description: 'Aggregated unverified items + risks; do not gloss.' },
  },
};

const synthesis = await agent(
  `${CONTEXT}\n\n=== YOUR ROLE: SYNTHESIS ===\nYou received ${good.length} grounded scoping areas (below). Produce a coherent scoping plan for wrapping COMPASS multi-aligner short-read novel-junction discovery into RECTIFY. Be skeptical: carry forward each area's UNKNOWNS honestly (do not present unverified tool availability or COMPASS claims as fact), and make the immediate re-validation experiment concrete enough to launch. Lead with the recommended aligner panel + the COMPASS design, then the RECTIFY integration, then the experiment, then the benchmark/design-doc feed. Keep every claim grounded; flag what still needs verification.\n\n=== SCOPING AREAS ===\n${JSON.stringify(good, null, 1)}`,
  { label: 'synthesis', phase: 'Synthesize', schema: SYNTH_SCHEMA, effort: 'high' });

return { scopes: good, synthesis };
