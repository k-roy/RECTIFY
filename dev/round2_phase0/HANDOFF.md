# HANDOFF — Round-2 cDNA Phase-0 empirical kill-gate

**Session:** 2026-06-14 · **Branch:** `walkback-guard-refactor-todo` · **Repo:** `~/work/rectify`
**One-line state:** Phase-0 apparatus is built + validated (15/15 tests), the human DRS substrate
is staged & deep, and **Step B (data-driven locus selection) is DONE** — 28 covered, short-exon /
multi-junction loci in `selected_loci.tsv`. But **the empirical GO/NO-GO verdict has NOT been
produced yet** (no real read aligned→lifted→scored). Next agent runs RUNBOOK Steps C→G.

Read first: `RUNBOOK.md` (this dir) and the MASTER spec
`dev/specs/SPEC_round2_cdna_discovery_assignment_20260614.md`. **Phase 0 gates everything: if it
returns NO-GO, STOP — do not build the production library/lift-over.**

---

## 1. What was done

- Read all 3 specs (master + 2 designer drafts); mapped the integration surface
  (`merge_corrected_tsvs`, the anchor gate, the `_cigar_*` scorers) via two Explore agents.
- **Substrate decision (Kevin): fresh HUMAN DRS run.** No existing human Round-1 5-aligner run
  with uLTRA exists on either cluster → build a fresh, locus-scoped one. Dataset chosen via
  research agent: **SG-NEx A549 directRNA** (rep1+4+5+6).
- Built the Phase-0 harness (`liftover.py`, `score_phase0.py`, `test_phase0.py`) — **all
  uncommitted** (see §5). Reuses the **production** `_cigar_*` scorers so the verdict is trustworthy.
- Fixed an advisor-caught lift-over bug (N dropped when a CIGAR op breaks exactly at an exon
  boundary; worst on micro-exon reads) + added a regression test.
- Fetched SG-NEx A549 rep1 genome BAM, subset chr5+6, ran a micro-exon coverage check → the 6
  a-priori micro-exon genes are **uncovered in A549**; pivoted to data-driven locus selection.
- Deepened the pool (rep1+4+5+6) → `a549_pooled_chr5_6.bam` (**COMPLETE**, see §2).
- **Step B DONE** (`select_loci.py`): from 1838 chr5/6 MANE transcripts, 205 have a short internal
  exon (≤30 bp); **28 pass coverage + inclusion** → `selected_loci.tsv`. Top: TGFBI (1895 spliced
  reads, 25 bp exon, 16 junc), IK, HSD17B4, UBE2B, ANXA6; dense-junction stressors TRIO (56 junc),
  XPO5 (31), SNX14 (28). (True ≤15 bp micro-exons that are *included* are scarce in A549 — SRRM4-low
  — so the set leans on 18–30 bp short exons + dense junctions; honest A549 reality, still a real
  seed-aligner stressor.)
- Wrote `RUNBOOK.md` (Steps A–G) and memory `project_round2_cdna_phase0.md`.

## 2. What's verified

- Tests: `cd ~/work/rectify && /Users/kevinroy/miniconda3/bin/python dev/round2_phase0/test_phase0.py`
  → **15/15 passed** (lift both strands, boundary-break N regression, ≥3 N-ops, E5/E6 I/D at
  junction, minus-strand multi-junction E11, orientation-invariant catch, BLOCKER-1 reject,
  gate-off raise, trivial-win defeat, micro-exon win, IR-lose, no-shrink veto, verdict GO).
- Substrate staged on Sherlock `/oak/.../Users/kevinroy/projects/rectify_round2_phase0/data/`:
  `a549_pooled_chr5_6.bam` — `samtools idxstats` → **chr5 = 664,021 reads, chr6 = 485,468**
  (rep1 alone had 41k+30k; pool is ~16×). rep1 chr5+6 had **27,357 spliced reads**.
- Refs confirmed on oak (Ensembl chr-naming `5`, matches BAM):
  `…/Projects/split_tap/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` +
  `Homo_sapiens.GRCh38.109.gtf.gz`. Aligners present in rectify env: samtools/minimap2/uLTRA/deSALT.
- **NOT VERIFIED (the whole point of Phase 0):** no cDNA library built, no Round-1 5-aligner run,
  no Round-2 alignment, no lift on a real read, **no GO/NO-GO verdict**. The 15/15 tests validate
  the machinery on synthetic data only — NOT the hypothesis.
- **NOT VERIFIED:** lift-over on real multi-error reads. RUNBOOK **Step F.0** (real-read
  round-trip) must pass before any verdict is trusted.

## 3. Open items

- **Locus selection DONE** (was the riskiest discovery step) → `selected_loci.tsv` (28 loci).
  Recommended Step-C starting subset (mix of depth + junction-density): **TGFBI, IK, HSD17B4,
  ANXA6, TRIO, XPO5** (+ DIAPH1 as a sanity cross-check — its a-priori 9 bp micro-exon was
  uncovered, but its MANE 27 bp internal exon IS included at 305 spliced reads). Don't re-run
  selection; `select_loci.py` is rerunnable if you want different thresholds
  (`MAX_SHORT/MIN_SPLICED/MIN_INCL_*`).
- **Harness is uncommitted** (`?? dev/round2_phase0/`). Why deferred: I don't commit without
  Kevin's say-so, and the branch HEAD moved during the session (concurrent-session caution).
  Stage explicitly (`git add dev/round2_phase0/`), never `git add -A`.
- **RNA002 not RNA004:** only public A549 DRS option. Fine for intron/micro-exon placement;
  flag if the verdict is marginal.
- **gmap missing** from the rectify env — it's only a Round-2 tie-breaker; minimap2-map-ont +
  mapPacBio are the Round-2 panel. Don't block on gmap.
- **Round-1 must be CHUNKED + gate ON** (`min_junction_anchor_bp=10`, `--partition=owners` +
  AVX-512 constraint) — CLAUDE.md mandatory rule; use the `sherlock-sbatch` skill.

## 4. Resume command

**Resume:** pool is built and Step B is done — go straight to RUNBOOK **Step C** (build the hand cDNA library for the selected loci).
1. Verify the harness still passes: `cd ~/work/rectify && /Users/kevinroy/miniconda3/bin/python dev/round2_phase0/test_phase0.py` (expect 15/15).
2. Read `selected_loci.tsv` (28 loci, committed in this dir AND at `<P0>/coverage/`). Start with **TGFBI, IK, HSD17B4, ANXA6, TRIO, XPO5** (+DIAPH1 cross-check).
3. Step C — for each selected gene, build the padded-cDNA + `liftover.Block` block-map from its MANE exon chain (Ensembl GTF.109) ∪ any de-novo gate-passed chain seen in Round-1; `CdnaModel.validate(genome)` the 1:1 identity; spike-test a non-canonical gate-passed junction survives.
4. Step D — extract the selected-locus reads from `<P0>/data/a549_pooled_chr5_6.bam` → fastq → CHUNKED 5-aligner Round-1 (`min_junction_anchor_bp=10`, owners + AVX-512, `sherlock-sbatch` skill), keeping uLTRA's per-read record.
5. Steps E→G per RUNBOOK. **Before any verdict, run Step F.0** (real-read lift round-trip) — if it fails, fix `liftover.py` before scoring. NO-GO ⇒ STOP and record the falsifier.
- No background jobs are running (the deepen `nohup` finished). Nothing to poll.
- ControlMaster discipline: only Kevin re-establishes the Sherlock master; batch SSH calls; if "Permission denied" appears, use `ssh -O check sherlock` (don't hammer).

## 5. Files touched

All under `~/work/rectify/dev/round2_phase0/` — **ALL [uncommitted]** (`?? dev/round2_phase0/`):
- `liftover.py` [uncommitted] — transcript→genome lift-over + N-op insertion + minus-strand + §2.5 invariants; `LiftedRead`, `Block`, `CdnaModel`.
- `score_phase0.py` [uncommitted] — win-guard + BLOCKER-1 fix + verdict; imports production `_cigar_*`.
- `test_phase0.py` [uncommitted] — 15 unit tests (run with base miniconda python).
- `RUNBOOK.md` [uncommitted] — Steps A–G, substrate, findings, GO criteria.
- `apriori_stress_loci_NOTES.md` — the (uncovered) a-priori loci, relocated from the Sumner scratch dir (wrong workspace).
- `select_loci.py` — Step-B selection script (rerunnable on Sherlock in rectify env).
- `selected_loci.tsv` — **the 28 selected loci** (Step B output; also at `<P0>/coverage/`).
- `HANDOFF.md` — this file.

(`liftover.py`/`score_phase0.py`/`test_phase0.py`/`RUNBOOK.md`/`apriori_stress_loci_NOTES.md`
were committed `c607524`; `select_loci.py`/`selected_loci.tsv`/this updated HANDOFF are the new
uncommitted delta — stage with `git add dev/round2_phase0/`.)

Cluster artifacts (Sherlock, `/oak/.../Users/kevinroy/projects/rectify_round2_phase0/`):
`data/a549_pooled_chr5_6.bam` (+rep1/4/5/6 sources), `select_loci.py`,
`coverage/selected_loci.tsv`, `coverage/microexon_inclusion.tsv`, `microexons.tsv`,
`01_fetch_subset_coverage.sh`, `02_deepen.sh`, `logs/`.

Memory written: `project_round2_cdna_phase0.md` (+ index line in `MEMORY.md`).
