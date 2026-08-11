# HANDOFF — human-transcript simulation → confirm workhorse flattening → scope the native aligner build

**For the next (fresh) agent.** Branch `worktree-agent-a25a2c1e784ad37dc` (worktree
`/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc`). Written
2026-07-02 by the outgoing Opus director as its context filled. Report to the PI
(Kevin); commit surgically (explicit paths); NEVER commit to `drs-validation-rebuild`.

> **Model note:** the PI asked for a "Fable" agent, but `model:fable` SILENTLY FALLS
> BACK to Opus 4.8 here (a Fable-SPECIFIC content SAFEGUARD flags benign biology — see
> the CORRECTION block in `HANDOFF_ALIGNER_BENCHMARK.md`; verify a subagent's REAL model
> from its transcript `"model"` field, NOT its self-report). Sonnet 4.6 IS served
> (verified, no fallback) if you want genuine model diversity. Write everything in
> plain, benign framing; keep intent clear (routine transcriptomics methods work).

---

## 1. Where we are (one screen)

RECTIFY corrects ONT direct-RNA/cDNA reads via a 5-aligner PANEL → per-read consensus.
We are gating a proposed **native de-novo aligner** (a 6th, orthogonal member) against
constructed truth. **Read first for the full story (plain English):**
`dev/NATIVE_ALIGNER_OVERVIEW.md` §"In plain English". Deeper: `dev/RECTIFY_STRATEGIC_FRAME.md`,
`dev/NOVEL_JUNCTION_BLINDSPOT.md`, `dev/REALDATA_NOVEL_JUNCTION_EVAL.md`, `dev/SUMNER_HUMAN_ALIGNER_EVAL.md`,
and the running log in `dev/HANDOFF_ALIGNER_BENCHMARK.md`.

**Settled findings (evidence, committed):**
- The gate refuted/deferred essentially every facet (C2/C3/C4/C5/C6, quality-axis) EXCEPT
  C1 (shipped) — and it surfaced+fixed a real consensus bug (the Discovery junction-proximity
  fix in `rectify/core/consensus/scoring.py`). The panel+consensus already handle most PLACEMENT well.
- **The real problem is DISCOVERY of novel non-canonical junctions.** minimap2 flattens
  **40–100%** of them (snaps to the nearest GT-AG) even on ERROR-FREE reads — confirmed on
  a REAL YEAST genome + multi-intron (`novel_junction_blindspot.py` / `novel_junction_realgenome.py`;
  see `dev/NOVEL_JUNCTION_BLINDSPOT.md`). This is the native aligner's reason to exist.
- **mapPacBio is NOT the answer.** On synthetic it "recovers" non-canonical (looked like the
  panel covers it), but on REAL A549 human data it FABRICATES: ~89–95% of its supported
  junctions are non-canonical + aligner-unique (real-human ~97% spurious; drops half its reads).
  Its synthetic "win" IS its real-data pathology (no splice-site sanity gate). GMAP likewise
  (rejected). See `dev/REALDATA_NOVEL_JUNCTION_EVAL.md` (real A549 per-aligner table).
- **Real-data verdict (SG-NEx A549 chr5, support≥2):** WORKHORSES (minimap2 1.6% / uLTRA 2.2%
  / winnowmap2 3.8% non-canonical) FLATTEN; deSALT best precise retainer (14%); mapPacBio/GMAP
  FABRICATE. So NO tool is both precise AND sensitive to novel non-canonical junctions — the gap.
- **HONEST CAVEAT (do not overstate):** the corrected #2 run showed the workhorse-MISSED
  non-canonical set is DOMINATED by GMAP+mapPacBio SHARED ARTIFACTS (48/50 top recall-gap =
  those two only), AND the workhorses DO recover ~686 real novel non-canonical junctions. So
  the ADDRESSABLE gap SIZE is unknown, pending cross-validation, and the CURE (a member recovers
  precisely) is UNPROVEN.

**PI decision:** BUILD the native aligner regardless (Feynman: "what I cannot create, I do
not understand"; mapPacBio was never updated for ONT). GUARDRAIL: fitness is ALWAYS the truth
set — build-to-learn must never drift into self-deception (the 0.09→1.07 trap). The
addressability ablation is a test we run ON the member as we build, not a blocker to starting.

---

## 2. TASK A — human-transcript simulation (do this first; it CONFIRMS the target with GROUND TRUTH)

**Why:** on real human data we can't get per-read truth for novel junctions (the tools hide
them). But if we SIMULATE reads from human-realistic transcript structures that INCLUDE the
non-canonical junctions we observed, we KNOW the truth and can measure exactly how many *real-type*
non-canonical events the workhorses miss. The yeast + real-A549 evidence predicts they flatten a
large fraction; this confirms it on human-realistic structure with ground truth.

**Design (concrete):**
1. Genome = the REAL human reference the A549 BAMs used (find it: the a549_chr5 BAMs were aligned
   to `sumner_lab/references/GRCh38_chr5.fa` — confirm via `@SQ`/prov.json in
   `/scratch/users/kevinroy/rectify_human_validation/sgnex_a549/alignments/`). Use chr5 (fast) or
   full GRCh38.
2. Transcript structures = REAL observed human junctions, not random-scanned. Source them from the
   real-data eval outputs `/scratch/users/kevinroy/nativealigner_realeval/jxn_eval_results.json`:
   - the `strong_noncanon_novel_top50` + `strong_noncanon_annotated_top30` (real-ish non-canonical
     junctions, workhorse-corroborated) → the "true non-canonical events" to simulate;
   - realistic intron-length + exon-length distributions from the observed junctions;
   - include the 3 COMPASS anchors (SLC35A4 non-canonical, TMED9, SQSTM1) as named ground-truth cases.
3. Build multi-exon transcripts on the human genome USING those real junction coordinates/motifs
   (mirror `scripts/benchmark/novel_junction_realgenome.py`, but feed OBSERVED junctions instead of
   the motif-scan). Read = spliced transcript from the real genome; truth = the (known) junctions.
4. Overlay REALISTIC ONT error (`scripts/benchmark/sim/error_injector.py` — but at a HUMAN-realistic
   RNA004 regime incl. the hot-read TAIL, NOT the too-benign 2% uniform that made the synthetic gate
   inconclusive; the tail is why the earlier rgfdr under-stressed mapPacBio). Also vary read length +
   exon size (short 5' exons) — the stressors the yeast run showed matter.
5. Align with the panel (minimap2 + deSALT + uLTRA + mapPacBio) to the human genome; score with
   `scripts/benchmark/panel_blindspot_score.py` (already has the RECALL + FDR axes: per-aligner
   `.rec` recall of the true junction + `.fpnc` spurious-non-canonical/read + the INTRONFREE control).
**Expected result (the confirmation):** workhorses recover the canonical/annotated junctions but
FLATTEN a large fraction of the true non-canonical ones (low `.rec` on non-canonical rungs);
mapPacBio recovers more but with high `.fpnc` (fabrication) — reproducing, WITH GROUND TRUTH on
human-realistic structure, what yeast + real-A549 showed. If workhorses miss a large real-non-canonical
fraction AND mapPacBio only "covers" it by fabricating → the native member's target is CONFIRMED with truth.

**Harness reuse (all built + committed):** `novel_junction_realgenome.py` (real-genome ladder, `--genome`
/`--n-introns`/`--exon`/`--error-rate`/`--intronfree`), `panel_blindspot_score.py` (recall+FDR),
`run_panel_on_corpus.py` (drives run_multi_aligner per aligner), `run_panel_blindspot.sbatch` (NJ_CORPUS-
parameterized, sentinel-rc). Cluster code synced at `/home/groups/larsms/users/kevinroy/nj_panel_code`;
corpora under `/scratch/users/kevinroy/nj_panel/<tag>/`. FIX NEEDED: deSALT+uLTRA failed on the tiny
synthetic contigs ("All requested aligners failed") — on a REAL human genome index they should work;
if not, debug run_multi_aligner index/annotation (uLTRA needs a GTF).

**Also (cheap, sharpens the target):** cross-validate the real-data recall-gap candidates
(`recall_gap_top50`, mostly GMAP+mapPacBio) against COMPASS SHORT reads (independent, non-ONT) to
split real-novel from shared-artifact → the honest addressable-gap SIZE. COMPASS data:
`/scratch/users/kevinroy/compass_a549/`.

---

## 3. TASK B — scope + start the native aligner build (after Task A confirms)

**Design in one phrase:** a Nanopore-calibrated, MOTIF-BLIND, −logP-scored RE-ALIGNER over the panel's
localization window — judge junctions by empirical EVIDENCE, not by "does it look GT-AG."

**Reuse (do NOT build a new DP from scratch):** the C1 engine already exists —
`rectify/core/align/local_aligner.py::align_exon_block_global(..., penalty_table, lam)` +
`rectify/core/splice/hp_penalty.py` (the empirical −logP HP/indel tables incl. per-cDNA-UMI-bin:
`rectify/data/genomes/*/penalty_tables/penalty_scores*.tsv`). The member REPOINTS this from HP-indel
placement to JUNCTION placement: take the panel's window (the 5 aligners' placement cluster), enumerate
candidate junction placements, score each by summed −logP emission over the read WITH the GT-AG/annotation
bonus set to ZERO (that's the orthogonality — it does NOT snap). The 2 canonical priors to NEUTRALIZE:
`_CANONICAL_HP_PRIOR` (junction_scoring.py:293, ±50bp knob, junction_refiner.py:647) and the `is_novel`
tiebreaker (junction_refiner.py:638). See the 3 Sonnet reviews `dev/REVIEW_SYNTHESIS_*_SONNET.md` +
`dev/DIRECTOR_ALGO_EVAL_{FABLE,FABLE_2,SONNET,SYNTHESIS}.md` for the member design consensus.

**THE make-or-break test (run it EARLY, on the member):** the ADDRESSABILITY ABLATION — on the flattened
non-canonical cases (from Task A, ground truth), does the −logP re-aligner recover the true junction at
HIGH recall AND LOW fabrication (`.fpnc`), vs a matched baseline? This is the C1-idiom ablation (arm
motif-blind-−logP vs arm flat/GT-AG-gated). PASS => the member earns its place. FAIL (can't beat panel
precision / flattened junctions are zero-evidence / a correct-step fix like the `_CANONICAL_HP_PRIOR`
ablation suffices) => document + it joins the refuted facets. Fitness = truth, never an internal score.

---

## 4. Loose ends / in-flight
- Real-data #2: a further "corrected" jxn_eval re-run may still be queued (job ~32531616); the KEY numbers
  are already in `dev/REALDATA_NOVEL_JUNCTION_EVAL.md` — don't block on it.
- Prior synthetic panel jobs (rgen/rgerr/err01/err05/short30/rgfdr) are DONE; rgfdr is the decisive one
  (recall+FDR; result `scripts/benchmark/rgfdr_result.txt`).
- The Discovery `_count_junction_proximity_errors` fix is SHIPPED (commit history) but owes a real-DRS recall spot-check.

## 5. Environment / discipline
- M1 gates: `PATH="/Users/kevinroy/miniconda3/bin:/opt/homebrew/bin:$PATH" PYTHONPATH=. /Users/kevinroy/miniconda3/envs/pysam/bin/python`
  (pysam env is THIN — no pandas/Bio; base env `/Users/kevinroy/miniconda3/bin/python` has full deps for tests).
  Smoke gate: `scripts/benchmark/smoke_roundtrip.py --out /tmp/x --reps 20` (exit 0). Tests via base env pytest.
- Sherlock: `ssh sherlock` (larsms; conda `/home/groups/larsms/users/kevinroy/anaconda3`, envs rectify/pbsim3).
  SERIALIZE ssh (MaxSessions throttles concurrent channels → transient "Permission denied" that SELF-HEALS —
  the master is 7-day ControlPersist, do NOT re-open it). rsync CHOKES on the login banner → use
  `ssh 'cat > dest' < file` or `tar czf - … | ssh 'tar xzf -'`. Sentinel-rc + robust watchers per CLAUDE.md.
- Budget/agents: ≤3 concurrent heavy agents, NO nested panel fan-out (over-spawning caused API stalls + a weekly
  limit). Prefer director-run INLINE gates. Advisor at every fork + before any verdict. Report-back-don't-commit.
- Consult the `advisor` before committing to the Task-A error model + the member scoring design (the two
  circularity-critical choices). RESUME `git -C <worktree> log --oneline -20` for the full trail.
