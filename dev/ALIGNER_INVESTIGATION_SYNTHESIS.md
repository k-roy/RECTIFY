# Aligner-investigation synthesis — integration with our native-aligner design + COMPASS rework (2026-06-18)

Cross-referencing the cloud-agent investigation (`dev/aligner_investigation/`, branch
`claude/rna-seq-aligners-investigation-bi79w9`) against this session's native-aligner design
(`ALIGNER_MEMBER_DESIGN.md`) + GMAP/111 work + the COMPASS rework. CORRECTIONS_vs_DRS_BUILD is the
authoritative override (see memory `reference-aligner-investigation`). Retracted claims NOT relayed as fact.

## The headline flip (changes our framing)
**Junction quality ALREADY drives production selection** — HP-edit-distance (Path A) is LIVE on the build
(`use_hp_ed=True`, both merge call sites pass raw BAMs). So our native design + COMPASS rework are
**sharpening a live metric, not wiring a dead one.** Two genuinely-open correctness gaps carry forward:
the **calibrated N-op cost** and the **missing positive splice-strength term** (below).

## Critical ADDs to the native-aligner design (from the investigation's converging analysis)
1. **Concept #1 (HP-length-law DP) MUST charge each N-op −logP(novel intron here), CALIBRATED — not 0, not a
   fixed gate.** N-ops cost 0 today (`corrected_consensus.py:142-143`, "free pass") = the live
   free-false-intron exploit (J1a, the #1 open correctness gap). A FIXED penalty re-introduces a gate, which
   violates the PERMANENT v3.1.7 no-candidate-gate policy (Module 2H is sequence-first, canonical = tiebreaker
   never gate — the design MUST honor this). So the N-cost must be a calibrated likelihood term in the DP.
2. **Concept #3 (posterior+LLR arbitration) MUST de-herd correlated aligners.** 3 of 5 panel aligners share
   minimap2 lineage (minimap2 native; gapmm2 wraps it; uLTRA falls back to it >10% out-of-annotation). Only
   **deSALT (RdBG) + mapPacBio (BBMap)** are structurally independent. Any independence assumption (LLR
   posterior, `_n_agree`, POA vote, COMPASS family-gate) is miscalibrated ~3× toward the minimap2 family.
   → use lineage / empirical-ρ weights; this is the single highest-value add to the arbitration concept.
3. **MISSED — a positive splice-strength term.** RECTIFY has only canonical/non-canonical TIERS, no positive
   donor/acceptor strength model. Add a yeast-retrained PWM/MaxEnt behind a swappable
   `splice_strength(chrom,pos,strand)` interface, as a log-odds term INSIDE concept #1's DP and concept #3's
   arbitration. minimap2's `--splice-flank=no` deliberately discards the flanking-strength prior and hands it
   to RECTIFY, which then doesn't model it. (Defer SpliceAI/Pangolin CNN to a metazoan port — ~nil value on
   yeast's ~300 canonical introns.)
4. **Concept #4 (POA refinement) — over-homogenization guards are MANDATORY, not optional.** POA inherits
   CONFLICT-3 severely: it erases genuine alt-5′/3′SS isoforms <50bp apart (real in yeast RP genes) and is
   self-reinforcing with any agreement-sensitive tiebreaker. REQUIRED: per-read sequence veto + tight
   `max_boundary_shift` + shared-boundary-anchoring (cluster only junctions sharing one exact boundary) +
   validate against ORTHOGONAL truth, never internal agreement.
5. **Concept #2 (CPA decoder) — scope DOWN to the realizable near-term version.** The investigation rates a
   unified terminal CPA solver HIGH-regression-risk (rewrites the hard-won 2B/2E/2G stack) and notes full
   dwell/Remora fusion is BLOCKED on POD5 squiggle retention RECTIFY doesn't keep. Realizable now: the
   **`pt:i` dorado tail-length Bayesian-shrink prior** in 2E walkback (tag already parsed,
   `drs_trim_command.py:411`) + **3′-tight/5′-loose DP asymmetry** + **mispriming terminal veto for cDNA**
   (reuses the genomic-downstream-A detector, DRS-exempt). Validate cat1/cat2 shifts exactly (regression).
6. **Concept #5 (FracMinHash fallback) — no investigation analogue; clean divergence, genuinely novel** (their
   seeding menu is syncmer/HPC/strobemer for short-exon recall, not containment localization).
7. **Two more MISSED additive items:** (a) scoped-HPC re-seeding of poly-A-proximal soft-clips (the one
   seeding idea in native territory; feeds concept #2); (b) a calibrated per-read `selection_confidence` +
   `_abstain` flag so contested coin-flip reads are filtered downstream (APA/DESeq2) not committed (pairs
   with concept #3's posterior).

## COMPASS rework refinements (from the same insights)
- **The ALIGNER_FAMILY dedup in our family-concordance gate is VINDICATED** — de-herding is exactly why
  modes/wrappers must not double-vote. For the human COMPASS panel, family-map STAR×2→star, HISAT2×2→hisat2,
  and weight by independence (minimap2-lineage members share a budget).
- **Gate on a CALIBRATED support score, not a raw count** (the herd trap): raw agreement counts shared
  systematic error as truth.
- **The orthogonal junction truth set is the program's true critical path** (the investigation flags it as
  the weeks-long real investment) — the COMPASS multi-aligner short-read detector IS that truth set. This is
  why the rework matters beyond the 111.

## GMAP / the 111 — refined framing (tempers my earlier read further)
- **The STAR-1-pass short-read test (14.5% sensitivity on known positives) is NEARLY UNINFORMATIVE for
  NON-validation:** if all 111 were real you'd expect only ~14.5% to validate. Report the 111's rate RELATIVE
  to the 14.5% floor, not the raw count. (My STAR-independent COVERAGE test — all 111 introns NOT spliced out,
  real splicing at different coords with thousands of reads — remains the STRONGER evidence and is NOT subject
  to this sensitivity caveat.)
- **LEAD with intrinsic evidence**, not the short-read test: the canonical-motif rate, GMAP's documented
  ONT-HP / poly-A blindness + per-read-isolation (no cross-read consensus), and the coverage test.
- **"Recurrent" does NOT rescue the 111** — systematic per-read GMAP misplacement reproduces the SAME false
  junction across reads, so recurrence is EXPECTED for a systematic artifact. Pre-empt this counter-argument.
- **Verdict unchanged: likely artifacts, INCONCLUSIVE.** The COMPASS sensitive multi-aligner detector (the
  rework) is the resolution. GMAP-not-panel-fit is independently corroborated (gmap dossier: ~30× slower, no
  ONT-HP tuning, no poly-A awareness, per-read isolation).

## ROADMAP §A junction-placement (re-scoped by CORRECTIONS)
- J0 win-rate provenance + Path-A/B ablation harness (now an ablation of two LIVE paths). DO FIRST; commit
  `aligner_summary.tsv` + `PROVENANCE.json`.
- **J1a calibrated N-op cost = HIGHEST junction priority** (J1b wiring already DONE on the build).
- J2 re-scoped: VALIDATE/VERSION the already-bundled `--Scer`-auto-resolved penalty tables (NOT "regenerate
  absent tables").
- J3 cross-read junction consensus de-biased by calibrated support (= concept #4 territory; CONFLICT-3 guards).
- J4 calibrated junction scorer with a yeast MaxEnt/PWM splice-strength term (= the MISSED positive term).
- Phase-1 O: procure an orthogonal junction truth set (= the COMPASS short-read detector). The real critical path.
- Dropped: mapPacBio maxindel (already set on build); demoted: unified terminal CPA solver (regression risk).
