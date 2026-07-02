# RECTIFY Native-Aligner Program — Director Algorithmic Evaluation
**Reviewer:** claude-sonnet-4-6 (independent; primary-material only; no cross-reading of prior Director docs)
**Date:** 2026-07-01
**Status:** DRAFT — EARLY WRITE; appendix sections will be added

---

## 1. Per-Aligner Teardown (Novel-Isoform Discovery Lens)

### 1.1 Minimap2 (canonical + `--junc-bed`)
**Engine family:** Seed-and-chain (k=14 k-mers → anchors → dynamic programming). Flat affine gap model; quality-blind.

**Splice prior (mechanism corrected after manual verification):** Two distinct biases are stacked in RECTIFY's minimap2 config:
- **(a) Motif snapping** — the `-x splice` preset enables `--splice-flank=yes`, which models the tendency of the +1/+2 positions of GT donors (A/G) and the −1/−2 positions of AG acceptors (C/T) as alignment signals. This is ALWAYS active in splice mode; it biases the DP toward canonical flanking nucleotides independently of any annotation file.
- **(b) Annotation snapping** — `--junc-bed <junc.bed> --junc-bonus 9` applies +9 to junctions matching the supplied annotation BED. This is NOT a general canonical-motif bonus; it rewards annotation-matching junctions specifically. [Verified: Ubuntu man page: "Score bonus for a splice donor or acceptor found in annotation (effective with --junc-bed)"; the `-x splice` motif preference is separate.]

Together these create a two-layer snapping force: motif snapping (structural, always active) plus annotation snapping (additive, annotation-file-conditional). For novel non-annotated AND non-canonical junctions, both layers penalize the true site.

**Error handling:** No per-base quality consumption in the affine DP. Homopolymer indels are penalized uniformly regardless of run length. The empirical penalty tables (RECTIFY's HP-law) show that at HP=8, a deletion costs ~0.033× a substitution — but minimap2's affine model charges a flat gap-open + gap-extend regardless. This means minimap2 over-penalizes deletions in long HP runs and under-penalizes them in short ones.

**Discovery failure modes:**
- *Novel non-canonical junctions:* The measured blind-spot (error-free reads) is 47% for real U12 AT-AC junctions, 70% for 1-off GA-AG, 78% for 2-off CT-AC, 90% for deep CA-TC. This is pure motif-snapping, not error-induced.
- *Novel NIC (novel intronic combination):* Annotation snapping via `--junc-bed` may bias toward known donor/acceptor combinations even when true splice is shifted by 1-3 bp.
- *Alt TSS / short 5' exons:* k=14 seed misses exons shorter than 14 nt; short terminal exons soft-clip. The Cat3 rescue in RECTIFY corrects many of these downstream, but the fundamental placement comes from minimap2 and is already compromised.
- *Cryptic pA extension/truncation:* Quality-blind model doesn't weight the poly-A tail region's high error correctly; the walk-back module handles post-hoc but alignment is imprecise at 3'.
- *Cis-variant reads:* A 40+ bp GT-AG-flanked deletion near a splice site is re-expressed as a fabricated intron (verified: C6 stratum, `fp_variant_adjacent=60` at 40bp driver). The flat cost model cannot distinguish a variant-containing deletion from an intron.

**Panel role:** Canonical snapper and anchor. The most sensitive and fastest. Functions as the "herd leader" — other panel members that share seed-chain + flat-affine converge with it.

### 1.2 gapmm2 (dropped as ~identical to minimap2)
**Status:** Explicitly noted in the codebase as "~identical to minimap2" for discovery purposes; dropped from the panel. Was designed for terminal-exon HP refinement via edlib post-processing — a correct-step add-on, not an independent alignment engine. Correlated placer confirmed, correctly excluded.

### 1.3 mapPacBio (BBTools/BBMap splice-aware)
**Engine family:** BBMap's hash-based local aligner with a splice-aware extension. Uses "perfect" hash seeds in an inverted index over the reference, then does local DP assembly for gapped alignments. Important distinction from minimap2: BBMap uses a different seeding strategy (exact-match hash, not minimizer) and potentially a different chaining model.

**Splice prior:** BBMap's splice-aware mode uses `intronlen` and `maxindel` parameters. It scores junctions via a canonical donor/acceptor scan (GT-AG/GC-AG/AT-AC detection) in post-processing. The `intronlen=10` + `maxindel=200000` wiring bug noted in the codebase (historical misconfiguration) was a real source of correlated failures with minimap2 at large-intron loci.

**Error handling:** BBMap's DP is also flat-affine; quality-blind by default. Does not use per-base Q for alignment decisions.

**Discovery failure modes:**
- Shares the flat-affine quality-blind family with minimap2 → likely correlated snapping bias on non-canonical junctions. Unconfirmed by the ladder (not yet run on the 5-aligner panel), but mechanistic reasoning is strong.
- Historical configuration bugs (intronlen) have been a source of non-orthogonal failure vs minimap2.
- On the chr5 A549 panel run, mapPacBio showed a "~1k introns total, a disjoint read set" from the misconfiguration — suggesting it can find DIFFERENT placements than minimap2, but those differences were partly configuration artifacts.

**Panel value for discovery:** Uncertain. If mapPacBio's hash-seeding + local DP genuinely differs from minimap2 on the motif-snapping axis, it could be partially orthogonal. But mechanistic grounds for believing this are weak. The pending panel blind-spot run will decide.

### 1.4 uLTRA (annotation-guided)
**Engine family:** Collinear chaining over an annotation-derived genome graph. Builds a splice graph from the GTF annotation; alignment is graph-constrained (must traverse exon nodes). Uses MEM (maximal exact match) seeds with collinear chaining on the graph.

**Splice prior:** The annotation graph IS the structural prior. uLTRA excels at exons too short for seed-chain aligners (11-20 nt) because the graph provides connectivity. But for truly novel junctions (nodes not in the annotation), uLTRA must either (a) extend along the "unguided" sequence path or (b) create novel edges — its behavior on unannotated junctions is annotation-constrained and likely conservative.

**Error handling:** Still uses standard DP scoring for the alignment within graph paths; quality-blind.

**Discovery failure modes:**
- *Annotation snapping:* The genome graph is the annotation — novel junctions that don't match annotated exon boundaries require uLTRA to score "off-graph" paths, for which it likely has no special capability. It may actually INCREASE annotation snapping vs minimap2 for truly novel NNC junctions.
- *Novel TSS / 5' ends:* Annotation-constrained; unannotated TSS exons have no graph node to chain to.
- *Non-canonical motifs:* Graph edges are derived from annotation, which is predominantly GT-AG. Novel non-canonical junctions are doubly disfavored: no graph node AND unusual motif.
- **For the blind-spot ladder: uLTRA is the WORST expected aligner on truly novel non-canonical junctions.** It provides annotation coverage diversity, not discovery diversity.

**Panel value for discovery:** High value for annotated isoforms with short exons; negative value for novel junction discovery (it adds annotation snapping on top of motif snapping).

### 1.5 deSALT (de Bruijn graph / deBGA-based)
**Engine family:** De Bruijn graph (RdBG) index built on the reference; seed-and-extend using de Bruijn graph walks. This is a fundamentally different seeding paradigm from k-mer minimizer indexing — deBG traversal can handle repeat regions and high-error reads differently.

**Splice prior:** deSALT uses a GT-AG/AT-AC/GC-AG scan in its splice-site detection. The DP scoring is still flat-affine. High-sensitivity mode: it reports more candidate splice sites than minimap2, but the final scored call uses similar canonical-motif preference.

**Error handling:** Quality-blind flat-affine DP (similar to minimap2). The deBG seeding may handle high-error regions better in terms of seed coverage (fewer seed deserts), but the scoring model has the same biases.

**Discovery failure modes:**
- Flat-affine motif-snapping: shares the core bias with minimap2.
- Fragility (SIGSEGV, duplicate alignments requiring dedup) means it drops reads in specific inputs — non-orthogonal failures where deSALT simply produces no output.
- The deBG seeding is a genuine difference: it may seed differently at repeat-adjacent junctions or in high-error reads, giving it some independence from minimap2's minimizer-based seed failures.

**Panel value for discovery:** Potentially the most orthogonal of the existing panel members on the seed-placement axis (deBG vs minimizer), but NOT on the scoring axis (both flat-affine). Its fragility makes it an unreliable discovery source. The pending panel blind-spot run will show whether deSALT independently recovers non-canonical junctions better than minimap2 (high-value measurement).

### 1.6 GMAP (currently uninstalled)
**Engine family:** GMAP uses a branch-point analysis approach: it identifies candidate donor/acceptor pairs using a polynomial-time dynamic programming over a sliding window, then computes a "gene model score" that includes canonical-splice-site penalties. Importantly, GMAP was designed for EST/cDNA spliced alignment and has a more explicit probabilistic model of splice signals than minimap2.

**Splice prior:** GMAP computes a genomic IIT (interval index tree) and uses donor/acceptor frequency tables. Its intron scoring includes both sequence features (the +1/+2 donor and −1/−2 acceptor positions) and intron-length distributions. This is MORE informative than minimap2's flat bonus, but still canonical-biased.

**Discovery value:** GMAP may find junctions minimap2 misses due to its explicit intron-model, but its canonical-bias is baked into the scoring. For novel non-canonical junctions, GMAP's explicit canonical bonus may be less manipulable than minimap2's flat affine — unclear direction.

**Assessment:** Not in the panel (uninstalled); the codebase has a wrapper ready. Worth evaluating before including — specifically, does it independently recover more non-canonical junctions than minimap2 on the deviation ladder?

---

## 2. Panel-Level Herding Analysis

### 2.1 Engine family convergence
The panel has two layers of herding:

**Layer 1 — Flat-affine scoring (ALL five aligners share this):**
Every panel member uses a flat gap-open + gap-extend model for local DP. None uses per-base quality. This means they all systematically over-penalize HP-region indels and apply motif bonuses that dominate quality-evidence on non-canonical junctions. They herd on the SCORING axis.

**Layer 2 — GT-AG motif preference (minimap2/gapmm2/deSALT/mapPacBio via splice model; uLTRA via annotation graph; minimap2 also has a third annotation-snapping layer via `--junc-bed`/`--junc-bonus`):**
Each panel aligner has its own mechanism but the outcome is shared: junctions near a canonical GT-AG motif score better than those at non-canonical sites, even on error-free reads. In minimap2 this is the `--splice-flank` model (motif snapping) stacked with `--junc-bonus` annotation snapping (verified: these are mechanistically distinct — see Section 1.1 correction). For truly novel non-canonical junctions, the entire panel biases toward canonical motifs. The measured 47-90% single-aligner blind-spot on minimap2 is likely replicated panel-wide (exact numbers pending from the cluster run), giving panel-native recovery well below 100% even for error-free reads of genuinely non-canonical biology.

**Where diversity exists (genuine):**
- **Seeding:** deBG (deSALT) vs minimizer (minimap2) vs hash (mapPacBio) vs MEM-graph (uLTRA). Independent seed failures → different unmapped-read populations.
- **Short-exon coverage:** uLTRA (annotation-graph) vs seed-chain (others). uLTRA genuinely recovers exons too short to seed minimap2.
- **Indel sensitivity at high error:** deSALT's deBG traversal may handle high-error regions differently (not measured).

**Where herding is structural (not fixable by parameter tuning):**
- Junction scoring: all use flat-affine + canonical motif bonus. No parameter change makes minimap2 or deSALT score non-canonical junctions without the GT-AG penalty.
- Quality blindness: no panel member can weight read evidence by empirical HP-law. This gap is structural.

### 2.2 The blind-spot confirms structural herding
The 47-90% minimap2 blind-spot on error-free reads is NOT correctable by running more panel members that share the same scoring family. The panel-native recovery for AT-AC (real U12 junctions) is bounded above by the fraction that deSALT/mapPacBio independently recover WITHOUT motif snapping — which, given their shared flat-affine motif-preference, is probably not dramatically better than minimap2.

**Key unresolved measurement:** the 5-aligner panel blind-spot ladder run on the cluster. Until this is done, the claim of "47-90% panel-native blind-spot" is extrapolated from minimap2 alone. The actual panel-native number could be anything from 5% (deSALT independently recovers most) to 85%+ (all five aligners snap). The native-aligner case is strongest if panel-native is still >30% for AT-AC.

---

## 3. Novel / Superior Algorithm Candidates — Ranked

### Rank 1: Calibrated-Likelihood DP with No Motif Prior (the proposed native member)
**Mechanism:** Replace flat affine with the empirical HP-penalty table (already exists: AT/CG × HP-length deletion, insertion, substitution costs). No GT-AG bonus. Score each putative junction by summing per-position empirical -log P over the read-vs-reference alignment through the junction window. The junction with the highest calibrated likelihood wins regardless of motif.

**Why orthogonal:** The scoring error axis is entirely different from the panel's shared flat-affine family. Where minimap2 charges a flat gap-open regardless of HP context, the calibrated DP charges near-zero cost for HP deletions (at HP=8, empirical del cost ≈ 0.033 vs substitution). More importantly: it has NO GT-AG motif bonus, so a non-canonical junction supported by clean read evidence wins over a canonical junction with evidence against it.

**How it exploits the empirical error model:** Directly: the HP-penalty table IS the deletion-dominant DRS error model. The `_score_hp_anchored` DP in `hp_penalty.py` already implements this; the native aligner extends it to junction scoring rather than just local exon scoring. This is not a new mechanism — it's the same DP used in Module 2H junction refinement, extended to INITIAL junction placement rather than just correction.

**cDNA-UMI exploitation:** For cDNA with UMI, reads of the same UMI family share a molecule; their per-read HP-penalty scores can be pooled (soft-union) to get a higher-confidence junction call from lower-quality reads. The calibrated DP makes this pooling principled: the joint log-likelihood over UMI family reads.

**Gate test design:**
- Gate: the deviation ladder, error-free (already measured for minimap2). Extend to calibrated-DP member.
- Metric: junction recall of TRUE site (as in `novel_junction_blindspot.py`), scored ambiguity-aware.
- Guard 1 — FDR on canonical reads: canonical GT-AG junctions must NOT be destabilized (recovery ≥ 0.98 on the control rung).
- Guard 2 — Addressability: error-free reads with clean flanks → the true site strictly wins on calibrated score even for non-canonical motifs (verify by scoring both true and snap site with the DP; the true site must win).
- Guard 3 — Error overlay: re-run ladder with the 3-layer injector at SIRV-calibrated parameters; confirm recovery degrades gracefully (not catastrophically) vs panel.
- Dependency cost: LIGHT — the HP-penalty table already ships; the DP already exists in `hp_penalty.py`. The new code is wiring the DP as a junction placer (not just a local corrector) over the panel's localization window.

**Assessment:** STRONGEST CASE. Directly addresses the measured blind-spot, exploits existing infrastructure, low dependency cost, clear gate design. The discovery ceiling is bounded by the localization window (panel must place the read in the right genomic region); junction recovery within that window is the native aligner's domain.

### Rank 2: Read-Quality-Stratified Novel Junction Discovery (the WS-3 guardrail, correctly positioned)
**Mechanism:** Compute per-read "hotness" (error density in exonic regions away from junctions) as a soft prior. Down-weight reads in the hot-read tail when calling novel junction support. A novel junction supported 90% by low-Q reads is suspicious; one supported 90% by high-Q reads is credible.

**Why orthogonal:** This is a meta-layer over the panel output, not an aligner per se. It exploits the RECTIFY error model (per-read quality is inferable from exonic regions) and targets a different failure mode: FDR inflation from noisy reads calling false novel junctions, not from motif snapping. It's complementary to the calibrated DP (which addresses snapping) by guarding against a different FDR source.

**How it exploits empirical error model:** The per-read hotness covariate has confirmed signal: real head-vs-tail autocorrelation r≈0.30 (measured on real BY4742 DRS). A hot read is globally unreliable for novel junction calls. The prior should be SOFT (per the measurement — r≈0.30 means per-transcript alignability also contributes; hard filtering would suppress reads from hard-to-align transcripts).

**Gate test:** Stratum G (novel junction discovery, WS-3 prototype) with the error overlay. Show that soft hot-read down-weighting reduces false novel junction calls without reducing true novel junction recovery proportionally. FDR lift must be >0 in the error overlay regime.

**Dependency cost:** Very light. The error injector and per-read hotness probe already exist. The guardrail is a post-processing weight on junction support counts.

**Assessment:** CORRECT POSITIONING as a guardrail, not a primary aligner feature. Build and gate alongside the calibrated DP as a FDR complement. It should NOT be the primary mechanism for discovery (it filters reads, not re-aligns junctions).

### Rank 3: POA-Pooled Paralog Consensus (C4 — SMN1/SMN2 + other multi-copy loci)
**Mechanism:** Cluster reads at a paralog locus by their primary placement (which copy they map to), then run partial-order alignment (POA) consensus within each cluster to denoise the distinguishing-SNP evidence. The pooled consensus call at the paralog-distinguishing position has higher confidence than any single read.

**Why orthogonal:** This addresses a completely different failure mode from motif snapping: per-read placement ambiguity at a locus where the distinguishing information is weak on any single read. It's not about junction scoring at all — it's about locus disambiguation.

**How it exploits the empirical error model:** Per-read Q and HP-penalty scores feed into the POA cost matrix. In a cDNA-UMI system, UMI families provide the clustering signal directly. For DRS (no UMI), clustering must rely on the aligner's primary placement confidence (MAPQ) + POA consistency.

**Gate test:** The PARALOG/C4 stratum is already built and VERIFIED DISCRIMINATING: weak 1-SNP fragments drop to locus_acc≈0.94 with minimap2, and pooling-majority call recovers the true copy in 6/6 pools from truth. The C4 ablation is to show this holds with the actual POA pipeline, not just with truth labels. FDR control: must not destabilize the strong-read arm (window-spanning reads at ceiling must stay at ceiling).

**Dependency cost:** Medium. Requires a POA library (abPOA or spoa, both Python-accessible). This IS a new dependency. The gate is already built; the implementation is the bottleneck.

**SMN1/SMN2 specific application:** SMN1/SMN2 differ at ~6 positions across the ~50 kb locus with the key exon-7 C6T change. Single reads that span only one distinguishing position have MAPQ 0-20 (minimap2 is uncertain). POA pooling over reads spanning the same distinguishing region denoises the signal. Combined with C6 (variant-aware alignment using a catalog VCF of SMN1 vs SMN2 distinguishing positions), this is the complete solution for SMN allele-specific isoform quantification.

**Assessment:** CRITICAL for the cis-variant / paralog-specific discovery use case, but NOT the highest-priority for general novel-isoform discovery. Build after Rank 1.

### Rank 4: FracMinHash Containment Fallback Localizer (C5)
**Mechanism:** For reads that NO panel member maps (unmapped or confidence-rejected), use MinHash sketch containment to find approximate genomic loci — similar to MASH Screen or sourmash. Then run the calibrated DP within those windows.

**Why partially orthogonal:** Truly orthogonal on the localization axis (no chaining, no splice model — pure k-mer containment). But it only adds value where ALL panel members fail to localize, which is the true panel-failure tail. The size of this tail is the gating measurement.

**Gate test:** Panel-failure tail measurement (Tier-2, cluster-scale). Need the full 5-aligner panel run on real reads to measure what fraction are placed by NO aligner. If this tail is <1%, the C5 case is weak. Tier-1 construction of "defeat-minimap2 reads" is not sufficient for C5 sizing.

**Dependency cost:** Medium. `datasketch` or `sourmash` for MinHash; otherwise uses existing infrastructure.

**Assessment:** Gate on tail size FIRST. If panel-failure tail < 2% of reads, C5 is a niche improvement. Deprioritize unless the cluster run reveals a larger tail.

### Rank 5: Variant-Aware Emission (C6)
**Mechanism:** Given a VCF of known variants (gnomAD/ClinVar/SMN paralog positions), modify the reference emitted-sequence at het/hom variant sites before running the DP. A 40bp deletion variant is scored as a known deletion (D) rather than triggering the intron-length threshold (fabricated N). Reduces the `fp_variant_adjacent` rate measured at 60 FPs in the C6 stratum.

**Why orthogonal (targeted):** Addresses a specific FDR source (variant-induced false introns) that is independent of motif snapping. The snapping problem is about the scorer preferring GT-AG; the variant problem is about the length threshold converting long deletions to introns.

**Gate test:** C6 stratum already built and discriminating (fp_variant_adjacent = 60 driver FPs vs 0 control FPs). The ablation: does variant-aware alignment reduce this to near 0? With correct VCF-informed DP, a 40bp deletion with GT-AG flanks is scored as a KNOWN deletion, not an intron candidate → the intron-length threshold is never triggered.

**Dependency cost:** Light to medium. Requires a VCF reader and reference modification; the DP already handles D ops correctly. The VCF itself is the external dependency — for SMN1/SMN2, the distinguishing positions are known and small.

**Assessment:** HIGH VALUE for the SMN1/SMN2 and clinical-variant isoform discovery use cases. Medium priority. Blocks the cis-variant research program but is not needed for the general novel-junction discovery case.

### Rejected approaches (with reasons)

**WFA-banded aligner as standalone member:** Shares minimap2's affine optimum by construction. Not orthogonal on the scoring axis. Correctly rejected already.

**Pangenome/variation-graph alignment (vg, Minigraph-cactus):** Requires a pre-built variation graph — heavy dependency (workflow-changing). No demonstrated discovery gain for novel junction sites vs the calibrated-DP approach. Paradigm rename with high dependency cost.

**Strobemer re-seeding:** Addresses seed coverage in repetitive regions; a seeding improvement, not a scoring improvement. Doesn't address the motif-snapping bias. Could be layered on top of the calibrated DP as a localization enhancement, but is not itself the orthogonal concept.

**Full-transcript generative model (HMM over transcript set):** High dependency (needs a transcript model for every possible transcript), computationally expensive, and for NOVEL transcripts (the very target) the model is underspecified. Conceptually elegant but not dependency-light and not gate-testable without a predefined model.

---

## 4. Cis-Variant, SMN1/SMN2, and Trans-Splicing Recommendations

### 4.1 Cis-variant isoforms (gnomAD/ClinVar catalog variants near splice sites)
**Correct-step vs native-aligner:** BOTH.
- **Correct-step (C6):** Variant-aware emission to prevent variant-induced false introns. This is a correction to the existing aligner output, not a new aligner feature. Build C6 into the correct-step pipeline: when a VCF is provided, re-score reads at variant positions with the known-variant allele before applying the intron-vs-deletion length threshold.
- **Native-aligner (calibrated DP):** The calibrated DP naturally handles variant-containing reads better than flat affine, because it scores the evidence (deletion cost at empirical rate) rather than applying a length threshold. A 40bp deletion in a non-HP region costs ~40 × 0.44 = 17.6 HP units, while a 40bp intron costs 0. The calibrated DP will still prefer the intron interpretation if the flanks are GT-AG — but this is a scoring preference, not a hard threshold flip. The GT-AG-free calibrated DP removes even this preference.
- **Recommendation:** Build C6 (variant-aware) as a correct-step feature in the correct pipeline; the native calibrated-DP member handles the residual.

### 4.2 SMN1/SMN2 and allele-specific/paralog isoforms
**Key challenge:** SMN1 and SMN2 differ at ~6 positions, with exon 7 containing the C6T change that drives exon-7 skipping in SMN2. Single-read placement is ambiguous (MAPQ 0-20 for reads spanning only one distinguishing position). Allele-specific isoform quantification requires:
1. Locus disambiguation (which copy) — C4 (POA pooling) or C6 (variant-aware alignment)
2. Junction calling within the correct copy — native calibrated DP
3. Isoform-level integration — downstream of alignment, correct-step

**Recommendation sequence:**
1. C6 (variant-aware emission with SMN distinguishing-position VCF) as the primary locus-disambiguation step: align to a two-copy reference (SMN1 + SMN2 loci concatenated) with known SNP positions soft-anchoring reads to the correct copy.
2. C4 (POA pooling) as the denoising step: pool reads per copy, then call the exon-7 junction at the copy level.
3. The calibrated DP provides junction-level calling within the correct copy, replacing flat-affine.

**Phase order:** C6 before C4 (C6 provides the locus anchor; C4 denoises; calibrated DP refines).

### 4.3 Chimeric/trans-splicing
**Classification:** Chimeric reads are a SEPARATE problem class from novel isoforms. A trans-spliced read spans two genomic loci; all aligners will either (a) place it at one locus with a large soft-clip, or (b) produce a split/supplementary alignment. The RECTIFY chimeric_consensus.py already handles chimeric stitching (uses D not N for gaps ≤10bp).

**Key points:**
- Trans-spliced reads appear in the PANEL-FAILURE tail (C5 domain) or as split-read candidates
- Detection requires: identify reads where no single-locus placement is satisfactory; find the two loci via MinHash or split-read alignment; stitch
- This is a CORRECT-STEP feature (split-read handler), not a native-aligner feature
- For ONT DRS, chimeric reads also arise from template switching and RT artifacts — these are noise, not biology. Need FDR control.

**Recommendation:** Handle trans-splicing as a post-consensus step: reads with evidence for two loci (high chimeric score) → split-read analysis → chimeric annotation. Do NOT conflate with the novel-junction discovery problem.

---

## 5. Top-3 Build-and-Gate-First Shortlist

### Priority 1: Full 5-Aligner Panel Blind-Spot Ladder (the DECISIVE measurement)
**Build nothing yet.** Run the `novel_junction_blindspot.py` script extended to all 5 panel members on the deviation ladder. This is the gating measurement that decides the native aligner's marginal gain.

**What it answers:** What is panel-native recovery for AT-AC (real U12) and GT-AG 1-off/2-off motifs? If panel-native ≥ 80% for AT-AC, the Rank-1 calibrated DP's marginal gain is small and the program should pivot to correct-step improvements. If panel-native < 50% for AT-AC, the Rank-1 calibrated DP is strongly justified.

**Cost:** Low (cluster run, existing script, minor extension for multi-aligner scoring). This should be the VERY NEXT step.

### Priority 2: Calibrated-DP Junction Placer (Rank 1, if panel blind-spot > 30%)
**Gate test:** Deviation ladder (error-free first; error overlay second). Gated on Priority 1 showing a real gap.

**Implementation path:**
- The calibrated DP kernel (`_score_hp_anchored` in `hp_penalty.py`) already exists as a numba-compiled function
- Extension needed: run it in JUNCTION DISCOVERY mode (scan candidate donor/acceptor positions within the panel's localization window; score each; return the best empirical-likelihood site)
- No new dependency: uses the shipped HP-penalty table
- Gate: deviation ladder recovery ≥ 0.85 on AT-AC (vs minimap2's 0.533), with FDR ≤ 0.02 on canonical reads

### Priority 3: Read-Quality-Stratified FDR Guard (Rank 2, concurrent with Priority 2)
**Gate test:** Stratum G (WS-3 prototype) with error overlay. Shows the guardrail reduces false novel-junction calls in high-error reads without suppressing true calls.

**Implementation path:** Lightweight — compute per-read exonic error density (already done in the reliability probe); use as a soft prior (NOT a hard filter) on junction-support counting.

**Gate:** In the error-overlay regime (injected errors), soft down-weighting hot reads reduces FDR by >10% with <5% recall cost on true novel junctions.

---

## 6. Synthesis and Short Assessment

**Most orthogonal approach:** Calibrated-likelihood DP with no GT-AG motif prior. It is the only candidate whose error axis is structurally independent of the panel's flat-affine/motif-snapping family. Low dependency cost (uses existing HP-penalty tables + existing DP kernel). The measured blind-spot (47% for real U12 AT-AC on error-free reads) is its exact addressable regime, and the gate design is clear.

**Highest-value panel weakness:** The shared flat-affine + GT-AG motif bonus across all panel members is a structural, non-tunable bias. The measured 47% blind-spot on real U12 AT-AC junctions is the headline number — these are biologically real minor-spliceosome junctions being silently flattened to canonical GT-AG calls on ERROR-FREE reads. For mutant/disease isoform discovery (SMN2 exon skipping, cryptic splice sites in cancer), this bias directly suppresses the biology of interest.

**What to gate-test first:** The 5-aligner PANEL blind-spot ladder (not just minimap2). Until the panel-native recovery at each deviation rung is measured, the calibrated DP's marginal gain is unquantified. This run is low-cost (extends existing script) and completely decisive.

**Whether the blind-spot justifies the native aligner:**
- If panel-native AT-AC recovery is < 50% (likely, given shared flat-affine family): YES, clearly justified. The calibrated DP fills a structural gap.
- If panel-native AT-AC recovery is > 80% (possible if deSALT/mapPacBio are genuinely orthogonal on this axis): the case WEAKENS for the general novel-junction facet; focus shifts to correct-step (dirty-pA walkback, C6 variant-aware) and domain-specific features (C4 paralog pooling for SMN1/SMN2).
- For U12 minor-spliceosome junctions specifically (AT-AC at 47% blind-spot on minimap2 alone): the case for a calibrated motif-free member is already established even before the panel run, because AT-AC is a REAL biological motif that should be recoverable.

**Key risk to the program:** The novel-junction blind-spot benchmark is currently minimap2-only. The "build signal" conclusion from `NOVEL_JUNCTION_BLINDSPOT.md` correctly notes this as caveat #1. If the cluster panel run shows deSALT independently recovers 70%+ of AT-AC junctions (plausible given its different seeding via deBG), the net panel blind-spot shrinks. Gate the calibrated-DP implementation on that run to avoid building into a gap that no longer exists.

**One additional risk:** The calibrated DP as a junction PLACER (not just corrector) requires the panel's localization window to contain the true locus. If the true locus is soft-clipped entirely (true panel failure, C5 domain), the calibrated DP cannot recover it. The proportion of reads in this regime (panel-failure tail) should be measured concurrently with the panel blind-spot run.

---

---

## 7. Web Survey — Published Prior Art (cited, 2021–2025)

Surveyed to: (a) confirm orthogonality claims against prior art; (b) identify competing approaches not in the codebase; (c) verify mechanism claims. Searches conducted 2026-07-01.

### 7.1 The AT-AC / U12 biology correction

The measured "47% blind-spot on real U12 AT-AC junctions" requires an important clarification. The minor (U12-type) spliceosome recognizes **both AT-AC and GT-AG introns**; AT-AC represents only a minority of U12-type introns (~25%), with the majority being GT-AG in sequence but recognized by U12 machinery. U12-type introns total < 1% of all introns. So "AT-AC = U12" is not wrong but is partial: most U12 biology operates through GT-AG introns that are indistinguishable in sequence from U1/U2 junctions. The AT-AC rung tests an extreme case; the broader U12 discovery problem is actually harder to gate (GT-AG U12 introns are invisible to a motif-snapping detection approach). This STRENGTHENS the blind-spot case: the portion of minor-spliceosome biology that IS AT-AC is already at 47% blind-spot; the GT-AG-sequence U12 fraction is invisible to the gate measurement entirely.

### 7.2 minimap2 splice model (mechanism correction confirmed)

Ubuntu man page (`manpages.ubuntu.com/manpages/focal/man1/minimap2.1.html`) confirms:
- `--junc-bonus`: "Score bonus for a splice donor or acceptor **found in annotation** (effective with `--junc-bed`)." — annotation-snapping, NOT motif-snapping.
- `--splice-flank=yes` (set by `-x splice`): uses the tendency of GT+1/+2 to be A/G and AG-1/-2 to be C/T as alignment signals — SEPARATE motif-snapping mechanism.
- These are mechanistically distinct; prior description in Section 1.1 conflated them. Corrected in place.

### 7.3 ESPRESSO (Yan et al., Science Advances 2023)

ESPRESSO (Error Statistics PRomoted Evaluator of Splice Site Options; doi:10.1126/sciadv.abq5072) considers per-read error profiles and aggregates evidence across reads to refine splice junction calls. Key algorithmic feature: "considers alignments of all long reads aligned to a gene and uses error profiles of individual reads to improve the identification of splice junctions."

Critically, ESPRESSO still gates high-confidence novel junctions on canonical motifs: "a putative splice junction is defined as high-confidence if it has the canonical splice site dinucleotide motif (GT/AG, GC/AG, or AT/AC) and is supported by at least two reads with perfect alignments around splice sites." This means ESPRESSO's error-profile aggregation is prior art for the multi-read evidence concept, but it RETAINS canonical-motif gating — it does not address truly non-canonical junctions. The calibrated-DP proposal is distinguished from ESPRESSO precisely by removing the motif prior entirely.

Performance: on synthetic spike-in and human RNA-seq, ESPRESSO outperforms minimap2/FLAMES/StringTie for isoform discovery/quantification; multi-read evidence aggregation has demonstrated value.

### 7.4 2passtools (Parker, Knop, Barton; Genome Biology 2021)

2passtools (doi:10.1186/s13059-021-02296-0; PMC7919322) applies alignment metrics and ML-derived sequence features to filter spurious splice junctions from long-read alignments, then uses the filtered set to guide a second-pass realignment. Improves intron detection accuracy for both annotated and unannotated junctions.

Key distinction from calibrated-DP: 2passtools FILTERS spurious canonical junctions using ML features (no quality-weighted scoring of alternative sites) and re-aligns using the filtered annotation as a guide. It does not independently score non-canonical motifs. The calibrated DP operates in a different regime: it re-scores junction ALTERNATIVES using HP-aware empirical costs, addressing motif-snapping bias at the scoring level rather than the filtering level.

### 7.5 IsoQuant (Prjibelski et al., Nature Biotechnology 2022)

IsoQuant (doi:10.1038/s41587-022-01565-y) is the leading reference-guided isoform discovery/quantification tool: 90% precision, 98% sensitivity in reference-guided mode (benchmark from Oxford Bioinformatics Advances 2025). It outperforms FLAIR and FLAMES for de novo discovery. IsoQuant is the transcript-level assembly/quantification layer; it runs downstream of alignment and does not address the junction-placement bias problem. It would benefit from a motif-free calibrated aligner producing its input.

### 7.6 Long-read aligner landscape (2021-2025 survey)

No published aligner (2021-2025) was found that uses HP-aware, quality-calibrated per-position DP scoring for RNA splice junction placement. Published tools in this period addressing long-read splice accuracy:
- **2passtools** (2021): ML junction filter + two-pass annotation guide — accuracy improvement via filtering, not scoring
- **uLTRA** (2021, Bioinformatics): annotation-graph collinear chaining — diversity on short-exon axis, worse on non-canonical novel junctions
- **ESPRESSO** (2023): multi-read error-profile aggregation with canonical-motif gating for novel junctions
- **IsoQuant** (2022): transcript-level isoform assembly, downstream of alignment

None of these use an empirical HP-law deletion cost model for junction placement. The calibrated-DP approach (using RECTIFY's empirical AT/CG × HP-length penalty tables in a no-motif-prior DP) has no direct published analogue. The closest conceptual neighbor is ESPRESSO's "error profiles for junction refinement," but ESPRESSO retains canonical gating and operates post-alignment, while the calibrated DP is proposed as an initial placer or refinement layer without motif gating.

### 7.7 No quality-aware splice aligner found

A targeted search for "base quality aware splice alignment long read ONT quality-weighted junction scoring 2023 2024 2025" returned no results for ONT splice-specific quality-weighted scoring. The benchmarking literature evaluates junction accuracy but does not use per-base quality as a scoring input for junction placement in any published aligner. This confirms the novelty of RECTIFY's per-position HP-law DP as a junction scoring mechanism.

---

*Deliverable written to `dev/DIRECTOR_ALGO_EVAL_SONNET.md`. Primary materials read: NOVEL_JUNCTION_BLINDSPOT.md, RECTIFY_STRATEGIC_FRAME.md, NATIVE_ALIGNER_OVERVIEW.md, SIMULATION_BENCHMARK_SPEC.md (full), docs/EMPIRICAL_HP_PENALTY_SCORING.md, docs/aligners/minimap2.md, docs/aligners/deSALT.md, docs/aligners/ultra.md, rectify/core/align/multi_aligner.py (header), rectify/core/correct/walkback.py (header), rectify/core/splice/hp_penalty.py (kernel), scripts/benchmark/novel_junction_blindspot.py (design). Independent assessment; Fable Director docs not read. Web survey conducted 2026-07-01 (see Section 7).*
