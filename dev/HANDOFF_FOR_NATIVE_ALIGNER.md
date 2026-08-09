# HANDOFF → Native Aligner Build agent (C1) — from the COMPASS/GMAP-validation session

**Date:** 2026-07-02. **Author:** compass-111-rna session (short-read COMPASS mode + GMAP appropriateness +
human-junction validation). **Audience:** the in-flight agent building the RECTIFY native aligner member (C1).

This is a briefing, not a task list — it gives you (1) the acceptance targets this session produced, (2) the
GMAP context that bears on the 5→2-3 panel cut, (3) where all the data/tools live, (4) the design + gate refs.

---

## 0. TL;DR — what changed under you this session

- **The benchmark GATE (Deliverable A) is built** (separate fork agent) on branch
  `worktree-agent-a25a2c1e784ad37dc` (off `drs-validation-rebuild`; benchmark-only paths). Smoke green:
  truth→sim→aligner→scorer round-trips; a 1bp-shifted call scores TP-not-FP; the live flat-affine DP scores
  0.980 on HP_HARD-noisy with 24/24 failures = the indel-vs-substitution misplacement **C1 targets**. One open
  item: the pbsim3 LIVE run wasn't executed (bioconda solve was slow). Resume in `dev/HANDOFF_ALIGNER_BENCHMARK.md`.
- **You now have 3 real human non-canonical junctions as a concrete C1 acceptance test** (see §1).
- **GMAP is confirmed low-unique-value** → it is the DROP target your C1 must justify retiring (see §2).

---

## 1. C1 ACCEPTANCE TEST — 3 real non-canonical junctions (the sharp one)

Validated real novel intragenic RNA splice junctions in A549 (DNA-confirmed genomic — NOT deletions/SV/circRNA/
fusion; near-triploid but flat-copy at these loci). Full derivation: `dev/COMPASS_2corroborated_CROSSPLATFORM.md`.
Each has a canonical GT-AG/GC-AG sitting **1–4 bp off** the placement every aligner (incl. accurate short reads)
chose — the exact ambiguity C1's motif-aware empirical-penalty realignment exists to resolve:

| junction (host gene, +) | dominant placement motif | canonical nearby | short-rd | long-rd corrob. |
|---|---|---|---|---|
| chr5:179824400-179832205 (SQSTM1, 7805bp) | `GT..GA` non-canon | GT-AG 4bp | 2959 | deSALT/uLTRA/GMAP |
| chr5:177592500-177593474 (TMED9, 974bp)   | `CG..CA` non-canon | GC-AG 1bp | 323  | mm2/deSALT/uLTRA/GMAP |
| chr5:140564954-140565547 (SLC35A4, 593bp) | `AG..CA` non-canon | GT-AG 3bp | 168  | mm2/deSALT/uLTRA |

**THE TEST:** given the supporting reads at each locus, does C1's realignment (a) **SNAP** the junction to the
nearby canonical motif → truth was canonical, seed-and-chain aligners mis-placed by breakpoint read-errors; or
(b) **HOLD** the non-canonical placement → genuinely non-canonical splicing? **Either answer is a PASS if
defensible from the reads.** The FAILURE mode is silently inheriting whichever the panel preferred. RT-PCR/Sanger
is the wet-lab arbiter. (4th junction chr5:179823051 = J2, a +2bp alt-acceptor of an annotated SQSTM1 intron,
3 short reads in ONT acceptor smear → INCONCLUSIVE; use as a "should abstain" case.)

Also in `rectify/data/validation/CASE_STUDIES.md` §"human non-canonical junction discovery" as acceptance criteria.

---

## 2. GMAP context for the 5→2-3 panel cut (C1's drop-justification)

Full analysis: `dev/GMAP_PANEL_APPROPRIATENESS.md`. Headline for you:

- GMAP is **low-unique-value**: on the corrected (ambiguity-tolerant) 5-chromosome validation (chr5+1/11/17/19),
  201 gmap-only candidates → 5 short-read-SUPPORTED → **~1 well-supported** (chr1:19219782-19221557, sr99,
  non-canonical CA..AC). ~97% of its output is non-canonical noise, fully suppressed by the `score_segment` fences.
- **The 3 real junctions above are NOT gmap-unique** — deSALT/uLTRA place them too (they were mislabeled
  "gmap-only" by exact-coordinate matching = the normalization gap). So GMAP earned ~no unique credit for them.
- **YOUR DROP-TEST (the clean way to retire GMAP from production):** show **C1 + minimap2 + (deSALT|mapPacBio)
  recovers GMAP's genuine unique finds (chr1:19219782-class non-canonical junctions) WITHOUT GMAP's noise.**
  Until you prove that, GMAP stays (fenced) — do NOT relax the `score_segment` non-canonical fences
  (`test_gmap_fence_regression.py` must keep passing; it's also your fence adversary — keep GMAP in the TEST panel).
- **Stated bias to design around:** our "real" gate required INDEPENDENT SHORT-READ corroboration, so a genuine
  LONG-READ-ONLY junction is undercounted as "artifact." C1 + the benchmark can adjudicate that long-read-only
  tail without short-read dependence — a capability the panel-vs-truth benchmark gives you that this session couldn't.
- De-herding note (from `ALIGNER_MEMBER_DESIGN.md` #3, benchmark-earned not settled): minimap2/gapmm2/uLTRA
  share lineage (uLTRA uses minimap2; gapmm2 wraps it); deSALT + mapPacBio are the independent algorithms.
  The A549 validation panel is minimap2/deSALT/uLTRA/mapPacBio/GMAP (deSALT, NOT gapmm2).

---

## 3. Where C1 lives + what to wire (from the design)

- **C1 = facet C1 of the native aligner member** = per-position empirical gap cost in the exon-block realignment
  DP `align_exon_block_global` (`rectify/core/align/local_aligner.py:522`), replacing the flat
  `_MATCH=2/_MISMATCH=-4/_GAP_OPEN=-4/_GAP_EXTEND=-1` with `−log P` from the empirical `HpPenaltyTable`.
  - **Tier A (MVP):** reuse `_precompute_del_costs` (`hp_penalty.py:315`) + add symmetric `_precompute_ins_costs`;
    inject per-ref-position cost into the D/I recurrence. Add `penalty_table: Optional[HpPenaltyTable]=None` —
    **None ⇒ byte-identical flat-affine (backward-compat fence).**
  - **Tier B (full):** run-unit length-transition emission (A12→A10 = one length error).
  - **Boundary invariance:** C1 only re-attributes indels WITHIN a block bounded by fixed junction edges → it
    structurally CANNOT create/move a junction (FDR axis = false-indels, not false-junctions).
  - REDTEAM must-fix already noted in the doc: `align_exon_block_global` (read_edits.py:811) is a DIFFERENT DP
    from the 5' rescue `align_clip_to_exon` (junction_refiner/splice_aware_5prime); C1 extends the FORMER only.
- Design: `dev/ALIGNER_MEMBER_DESIGN.md` (§2 = C1; §7 shared calibration; §10 phased roadmap; §8 FDR; appendix redteam).
- Gate/scope: `dev/ALIGNER_PROGRAM_SCOPING.md` (member build gated on the benchmark; fitness = truth set NEVER internal score).
- Empirical tables: `docs/EMPIRICAL_HP_PENALTY_SCORING.md`; bundled `rectify/data/genomes/*/penalty_tables/`.

---

## 4. Recommended sequence for you (C1)

1. **Finish the benchmark's pbsim3 live run first** (the gate) — `dev/HANDOFF_ALIGNER_BENCHMARK.md` resume. No
   C1 fitness claim is valid until the truth-set benchmark is live. Then C1 Phase 0 (the `penalty_table=None`
   byte-identical equivalence harness) → Phase 1 (per-position gap cost; must STRICTLY beat flat-affine AND the
   existing `homo_mismatch` knob on held-out exact-indel-position concordance).
2. **Add a benchmark stratum** for the §1 failure mode: "canonical site flanked by an aligner-preferred
   NON-canonical placement N bp away." This session's 3 junctions are the real-data instances; the benchmark
   needs synthetic ones (with truth) to score C1's snap-or-hold behavior at controlled FDR.
3. **Run C1 on the §1 junctions** once it can realign — report snap-vs-hold per junction (the acceptance test).
4. **Run the GMAP drop-test** (§2) once C1 is validated — this is the deliverable that justifies 5→2-3.

---

## 5. FILES / DATA (all on Sherlock unless noted)

- `$W = /scratch/users/kevinroy/compass_a549`; `$OAK = /oak/stanford/groups/larsms/Users/kevinroy/compass_a549_out`.
- Verdicts: `$OAK/adjudication_111_v3.json` (chmod 444; note: its "3 real" are multi-aligner not gmap-unique).
- 5-chrom validation: `$W/reval_chr{5,1,11,17,19}.json`; confirm: `$W/confirm_13.json`.
- Long-read panels: chr5 = `…/sgnex_a549/alignments/a549_chr5_trimmed.{minimap2,deSALT,uLTRA,mapPacBio,GMAP}.bam`;
  chr1/11/17/19 = `…/sgnex_a549/morechrom/alignments/morechrom_trimmed.{minimap2,deSALT,uLTRA}.bam` + `$W/morechrom_trimmed.GMAP.bam`.
- WGS (DNA confirm): `$W/wgs/a549_wgs_deep.bam` (~10×, ENCODE ENCSR521ELB).
- Tools: `dev/gmap_validate_harness.py` (chrom-agnostic; STEP-1 fixed to ambiguity-tolerant), `dev/confirm13*.py`,
  `dev/lr_probe_4loci.py`, `dev/dna_split.py`, `dev/rederive_111.py`.
- Docs (this session): `dev/GMAP_PANEL_APPROPRIATENESS.md`, `dev/COMPASS_2corroborated_CROSSPLATFORM.md`,
  `dev/HANDOFF_SHORTREAD_P5.md` (full chronology).
- **SSH note:** the M1's Sherlock ControlMaster is flaky (no Stanford Kerberos ticket → every reconnect Duos;
  master dies on laptop sleep). Durable fix pending (SSH key or Stanford TGT). Cluster jobs are unaffected.

---

## 6. Open (non-blocking) follow-ups this session left

- Optional genome-wide GMAP confirm: harness ready (`--chrom all`) but needs a genome-wide GMAP+indep LR panel
  (only chr5 + morechrom chr1/11/17/19 aligned). 5-chrom is already decision-grade → LOW marginal value.
- One-line note in `adjudication_111_v3.json` that the 3 real junctions are multi-aligner (reality unchanged).
- Harness follow-up: the confirm helper OOMs if it caches many full-chrom junction dicts — free per-chrom
  (`dev/confirm_chr19.py` shows the lean pattern).
