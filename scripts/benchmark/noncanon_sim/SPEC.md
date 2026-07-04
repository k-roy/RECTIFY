# SPEC — yeast non-canonical LONG-READ ground-truth sim + arm-A/B/C vetting (coordination contract)

Goal: a **ground-truth simulator of non-canonical splice junctions in ONT-like LONG reads**, in yeast,
plus a harness that vets three re-aligner arms and produces a **recovery-vs-false-junction-FDR trade curve**.
Motivation: prp18Δ (Roy et al. 2023 NAR) activates alternative 3′SS including **non-YAG acceptors**; a
motif-biased aligner would SNAP a real non-YAG acceptor to a nearby canonical YAG, flattening the signal.

Work in this worktree: `/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc`.
Put ALL new code under `scripts/benchmark/noncanon_sim/`. Base python (has pysam, no pandas/numba needed):
`/Users/kevinroy/miniconda3/bin/python`. minimap2 + samtools are on PATH (M1). Yeast genome bundled:
`rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa{,.gz}`;
annotation `saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz`.

## BIOLOGY — get this right (PI-corrected)
- The non-canonical event is at the **3′ splice site (acceptor)**. Donor stays canonical **GT**.
- The **poly(U) tract and first-intronic-nt complementarity are INTRONIC** → spliced out → **NOT in the read**.
  Do NOT put poly(U) in the read. They only motivate WHICH cryptic 3′SS is used.
- The read is a spliced **exon1 ++ exon2** (no intron). The aligner must place the intron (N-op); the
  **acceptor position determines the genomic 3′SS dinucleotide**. Snap = aligner shifts the acceptor to a
  nearby canonical **YAG** in the genome. Hold = keeps the true non-YAG acceptor.
- Ambiguity source = **exonic sequence near the junction + realistic long-read errors** + a canonical-YAG
  decoy a few bp from the true acceptor.

## PANEL DESIGN (by-construction truth)
Axes (Cartesian panel; each cell ≥ N reads, N≥100):
- **acceptor motif rung** (increasing non-canonicity of the TRUE acceptor):
  R0 = `YAG` (canonical CONTROL — must never be moved by any arm),
  R1 = `RAG` (AAG/GAG, tier-1 semi-canonical),
  R2 = `non-RAG NAG` (e.g. ...AG but with a non-canonical −1, still ends AG but low tier),
  R3 = `non-AG` (a genuinely non-canonical acceptor dinucleotide, e.g. AC/AA/AT — the deep prp18Δ case).
- **decoy offset**: a canonical YAG placed at ±k bp from the true acceptor, k ∈ {2, 3, 4} (the snap target).
- **exon sizes**: short vs long flanking exons, e.g. exon5 ∈ {40, 150} nt, exon3 ∈ {40, 150} nt (anchor length).
- **junction-flank context**: (a) plain unique sequence, (b) an EXONIC homopolymer run (A/T, len 5–8) within
  ~3 bp of the junction (the realistic HP-ambiguity case). NOTE: HP is EXONIC here, not the intronic poly(U).
- **intron length**: yeast-scale, ∈ {~90, ~250, ~450} nt.

FDR CONTROLS (must be in the panel):
- **INTRONFREE**: single-exon reads spanning the locus, NO true junction. Any arm emitting a junction = FP.
- **YAG-canonical, no decoy**: true acceptor is canonical YAG with NO nearby non-canonical site. Any arm
  MOVING it to a non-canonical placement = false-junction FP.

## INTERFACE CONTRACT (the files that flow A → B → C → D)
1. `templates.fa` — spliced transcript sequences (one per truth row), header `>tid_<id>`.
2. `sim_ref.fa` — the GENOMIC reference the reads align to (yeast chroms or synthetic contigs); reads' truth
   junction coords are IN THIS FRAME. If synthetic contigs, name them `chrSIM_<n>`.
3. `panel_truth.tsv` (tab-sep, header) — one row per template:
   `tid  chrom  true_donor  true_acceptor  strand  intron_len  motif_rung  acceptor_motif  decoy_offset
    decoy_acceptor  exon5_len  exon3_len  context  has_true_junction`
   (coords 0-based half-open intron [true_donor, true_acceptor); has_true_junction=0 for INTRONFREE.)
4. `reads.fastq` — ONT-like long reads; read id = `<tid>_r<NNN>`.
5. `read_truth.tsv` (tab-sep, header) — one row per READ, joining read→template truth:
   `read_id  tid  chrom  true_donor  true_acceptor  strand  motif_rung  context  has_true_junction`
6. `arm_{A,B,C}.bam` — sorted+indexed, reads aligned then junction-refined per arm.
7. `trade_curve.json` — per (motif_rung × context [× decoy_offset]) cell: recovery + false_junction_FDR for
   each arm, plus overall + the control FP rates.

## APIS (do not reinvent)
- Re-aligner: `from rectify.core.splice.junction_refiner import refine_bam_junctions`. Arms:
  A = `refine_bam_junctions(..., motif_blind=False)` (incumbent),
  B = `refine_bam_junctions(..., motif_blind=True)` (motif-blind — JUST ADDED, byte-identical when False),
  C = B + `penalty_table_path=<yeast penalty_scores.tsv>` (the −logP law).
  Call with `n_workers=1` (small panels). It needs `aligner_bams` (the per-read alignment BAM(s)) +
  `annotated_junctions` (build from the yeast GFF via `load_annotated_junctions`) + `genome` dict.
- Scorer: `rectify/core/benchmark/scorer.py` — `extract_junctions`, `net_indel_in_span`, ambiguity-aware
  matching via `chimeric_consensus.normalize_junction` / `_canonical_within_window`. Read it; reuse it.
- Error model (INDEPENDENCE — load-bearing): reads' errors MUST come from a source INDEPENDENT of the
  scorer's `penalty_scores.tsv`. Prefer **pbsim3** (`scripts/benchmark/sim/pbsim3_wrapper.py`, ERRHMM-ONT —
  check if pbsim3 is on PATH locally; if not, note it and use the fallback). Fallback = clean reads from
  templates + a documented independent per-base ONT error profile (marginal del≫ins~sub, with an EXONIC
  HP-context deletion elevation) — NOT derived from penalty_scores.tsv. Never feed penalty_scores into read
  generation.

## METRICS (the deliverable — advisor-mandated framing)
- **Recovery** (per non-canonical cell): fraction of reads whose refined junction matches the TRUE non-YAG
  acceptor (ambiguity-aware, within the normalize window). Expected: arm-A LOW (snaps to decoy YAG), arm-B/C
  HIGH. But recovery alone is near-tautological (arm-B removed the prior) — so:
- **False-junction FDR**: on the CONTROLS (INTRONFREE + YAG-canonical-no-decoy) and on R0, fraction of reads
  where an arm emits/moves-to a non-canonical junction that ISN'T true. This is what arm-B/C must NOT inflate.
- **LEAD WITH THE TRADE CURVE**: recovery vs false-junction-FDR, stratified by motif_rung × context (× decoy
  offset). The aligner earns its keep only if recovery rises without proportional FDR. Report per-arm curves.

## STATUS OF THE ARMS
- Arm-B toggle (`motif_blind`) is wired into junction_refiner.py THIS session, byte-identical when off,
  59 tests green (uncommitted). Arm-C uses the shipped yeast penalty table + C1 length-law (already built).
