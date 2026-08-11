# Spike-in track — DESIGN (2026-07-07; PI: "both, spike-in first")

**Goal.** The make-or-break Q2 (does motif-blind re-placement RECOVER real novel non-canonical junctions that
minimap2 flattens, PRECISELY?) with GROUND TRUTH immune to the motif-circularity that rigs the short-read Q2.
Spike SYNTHETIC reads carrying KNOWN unannotated non-canonical junctions into a REAL ONT-dRNA background →
truth is the spike design, not any aligner's call. (COMPASS revival later for real-biology corroboration.)

**Survey basis (dev/NONCANON_BENCHMARK_SURVEY.md):** SQANTI-SIM is the published precedent (novel junctions by
transcript-removal → NanoSim/pbsim3 reads); SIRV-in-real-dRNA (SG-NEx RNA002 / LongBench RNA004) is the
real-background precedent; the SPECIFIC combo — UNANNOTATED NON-CANONICAL in a real dRNA background — is
genuinely novel design space (the defensible contribution). Metric convention = SQANTI3 NIC/NNC + motif-stratified.

## Design (SIRV-model, synthetic contigs — clean truth separation)
1. **Synthetic transcripts + junctions:** reuse `noncanon_sim/build_panel.py` — non-canonical (non-YAG, GC-AG,
   AT-AC, cryptic 3'SS) + unannotated junctions across the graded axes we already build (motif-deviation rung ×
   HP context × canonical-DECOY presence × flanking exon size). These become synthetic TRANSCRIPT sequences.
2. **Reference:** add the synthetic transcript source sequences as EXTRA CONTIGS to the background genome
   (SIRV-model). Spiked reads align to the synthetic contigs (unambiguous truth); real background reads align to
   the real genome. No collision.
3. **Spiked reads:** `gen_reads.py` (our VALIDATED ONT-dRNA error model — pbsim3 ERRHMM or the fallback with
   HP undercalls) from the synthetic transcripts, at a controlled fraction (~1-5%) of the background run.
4. **Pipeline:** align the combined FASTQ with the PANEL (minimap2 -ax splice + others) to the combined
   reference; refine to the FOUR arms — **raw-minimap2 → arm-A (motif-biased) → arm-B (motif-blind) →
   arm-Bguard**. (arm-A is the vet's mandatory control: motif-blind contribution = arm-B − arm-A, NEVER
   arm-B − raw-minimap2, which confounds motif-blindness with the re-placement engine.)
5. **Metrics (on the synthetic contigs, known truth), SQANTI3-style NIC/NNC + MOTIF-STRATIFIED:**
   - **RECALL** = fraction of spiked non-canonical junctions each arm places at the TRUE site (vs snapping to a
     nearby canonical decoy). The flattening curve: raw-mm2 (worst) → arm-A → arm-B (best). Stratify by
     motif-deviation rung + decoy presence (no-decoy = clean; with-decoy = the hard flattening case).
   - **PRECISION / fabrication** = junctions each arm calls on the synthetic contigs that are NOT true (the
     guard's job: arm-Bguard should cut fabrication vs arm-B without cutting recall).
   - Position-exact AND windowed (±ambiguity), ambiguity-aware normalization (our normalize_junction).

## Background-data choice (a decision — see the AskUserQuestion)
- **LongBench RNA004** (`s3://longbench-data/`) — RECTIFY's TARGET chemistry (RNA004); best realism; needs a
  download. The strongest external-validity background.
- **SG-NEx A549 RNA002** — already staged (`sgnex_a549_illumina` + the A549 DRS); older chemistry (RNA002).
- **Reuse the BY4742 DRS** we already have (RNA004, native organism) — cheapest; the background genome is
  irrelevant to the spiked-contig truth, so this is a valid fast first pass.
NOTE: the spiked reads carry SIMULATED error (our validated model), so the background's role is run-context
realism (read-length/quality mix, panel behavior in a real run), not real error on the spiked junctions —
that distinction is honest and worth stating. Real-error-on-non-canonical would need wet-lab custom SIRVs (out of scope).

## What this adds beyond the pure sim (which already showed minimap2 flattens 40-100% + arm-B recovers)
- A FABRICATION/precision measurement in a realistic context (pure sim measured recall; this adds precision +
  the arm-A control + the SQANTI3 convention).
- External validity: the panel + refiner behaving in a real-run mix, not isolation.
- A publishable, SQANTI-SIM-lineage benchmark (the novel unannotated-non-canonical-in-real-dRNA contribution).

## First concrete build (after the background decision)
1. Build a spike panel (build_panel: non-canonical × decoy × HP, ~synthetic transcripts) → synthetic contigs FASTA.
2. gen_reads → spiked FASTQ; cat into a subset of the chosen background run at ~2-5%.
3. Align (panel) to (background genome + synthetic contigs); refine arm-A/B/Bguard.
4. Score recall/precision on the synthetic contigs, motif-stratified, arm ladder. Adopt SQANTI3 category labels.
Reuses: build_panel.py, gen_reads.py, run_arms.py (add arm-A + the combined-reference + spike-mixing wiring).
