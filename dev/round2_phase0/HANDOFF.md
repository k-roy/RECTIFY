# HANDOFF — Round-2 cDNA Phase-0 empirical kill-gate

**Session:** 2026-06-14 · **Branch:** `walkback-guard-refactor-todo` · **Repo:** `~/work/rectify`
**One-line state:** Phase-0 apparatus is built + validated (15/15 tests) and the human DRS
substrate is staged & deep — but **the empirical GO/NO-GO verdict has NOT been produced yet**
(no real read has been aligned→lifted→scored). Next agent runs RUNBOOK Steps B→G.

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

- **Data-driven locus selection still owed** (RUNBOOK Step B). Why this way: the 6 a-priori
  micro-exon genes (DIAPH1/MYO10/ABLIM3/KIF3A/TRIO/MPC1 — `apriori_stress_loci_NOTES.md`) have
  ~0 A549 coverage, and A549 is SRRM4-low so neuronal microexons are skipped. Must select loci by
  **coverage ∩ short-internal-exon structure** in the pooled BAM, not an a-priori list. Top
  covered spliced loci already seen: EEF1A1, ATG10, RPS18, TRIM41.
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

**Resume:** the depth pool is DONE — go straight to RUNBOOK **Step B**.
1. Confirm pool: `ssh sherlock 'bash --norc --noprofile -c "source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh; conda activate rectify; samtools idxstats /oak/stanford/groups/larsms/Users/kevinroy/projects/rectify_round2_phase0/data/a549_pooled_chr5_6.bam | head"'` → expect chr5≈664k, chr6≈486k reads.
2. Run Step B (data-driven locus select): scan the Ensembl GTF.109 for internal exons ≤30 bp on chr5/6, intersect with per-gene spliced-read depth in the pooled BAM, keep genes with ≥~50 spliced reads AND a short internal exon, **confirm short-exon inclusion** via `samtools depth` (PSI proxy). Write `coverage/selected_loci.tsv`.
3. Then Steps C→G per RUNBOOK. Before any verdict, run **Step F.0** (real-read lift round-trip) — if it fails, fix `liftover.py` before scoring.
4. Verify the harness still passes first: `cd ~/work/rectify && /Users/kevinroy/miniconda3/bin/python dev/round2_phase0/test_phase0.py` (expect 15/15).
- No background jobs are still running (the deepen `nohup` finished). Nothing to poll.
- ControlMaster discipline: only Kevin re-establishes the Sherlock master; batch SSH calls; if "Permission denied" appears, use `ssh -O check sherlock` (don't hammer).

## 5. Files touched

All under `~/work/rectify/dev/round2_phase0/` — **ALL [uncommitted]** (`?? dev/round2_phase0/`):
- `liftover.py` [uncommitted] — transcript→genome lift-over + N-op insertion + minus-strand + §2.5 invariants; `LiftedRead`, `Block`, `CdnaModel`.
- `score_phase0.py` [uncommitted] — win-guard + BLOCKER-1 fix + verdict; imports production `_cigar_*`.
- `test_phase0.py` [uncommitted] — 15 unit tests (run with base miniconda python).
- `RUNBOOK.md` [uncommitted] — Steps A–G, substrate, findings, GO criteria.
- `apriori_stress_loci_NOTES.md` [uncommitted] — the (uncovered) a-priori loci, relocated from the Sumner scratch dir (wrong workspace).
- `HANDOFF.md` [uncommitted] — this file.

Cluster artifacts (Sherlock, `/oak/.../Users/kevinroy/projects/rectify_round2_phase0/`):
`data/a549_pooled_chr5_6.bam` (+rep1/4/5/6 sources), `microexons.tsv`,
`coverage/microexon_inclusion.tsv`, `01_fetch_subset_coverage.sh`, `02_deepen.sh`, `logs/`.

Memory written: `project_round2_cdna_phase0.md` (+ index line in `MEMORY.md`).
