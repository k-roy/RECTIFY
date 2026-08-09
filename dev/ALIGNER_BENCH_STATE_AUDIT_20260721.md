# Aligner Benchmark — State Audit (2026-07-21)

> ⚠ PROVENANCE NOTE (2026-08-09): this document was written in the benchmark agent worktree
> (`.claude/worktrees/agent-a25a2c1e784ad37dc`, branch `worktree-agent-a25a2c1e784ad37dc`) and was
> never committed there. The physical worktree directory was deleted on 2026-08-09 (all
> `.claude/worktrees/` dirs were removed by an external cleanup while the branch + commits
> survive in the shared object store). This copy was restored verbatim from a same-day read of
> the file and committed into the MAIN tree so the definitive program state survives.

Reconstruction of the "gate-before-member" native-aligner program's true state, for
the orchestrator to pick the next increment without redoing finished work. Read+
reconstruct only; no code modified. Citations are `file:line-range` / commit hashes.

## TL;DR — the spine (read this first)

The orchestrator's brief is anchored to the **STALE 2026-06-30 SESSION-12 handoff**.
That handoff described a 6-facet "gate-before-member" program (C1–C6) with three facets
"lost in an 8-agent fan-out, to re-run inline." **That phase is now CLOSED, and the
program was REOPENED on 2026-07-01 (SESSION-13) onto a completely different axis.** The
live program is no longer "add a 6th consensus panel member scoring on the −logP axis";
it is a **motif-blind junction-DISCOVERY re-aligner implemented as a junction refiner
(Module 2H, `junction_refiner.py`)**, validated on REAL Nanopore data (Sumner SMA A549
chr5 DRS + short-read corroboration). As of **2026-07-17 (Phase 6)** the guard track
built on top of it was **overturned by a real-data gate**.

Three things the brief assumes that are now WRONG (contradiction flags, detailed below):
1. C6 / Discovery / Flat-Q were **NOT lost** — all three were committed same-day
   2026-06-30 (`28f6aac`, `fc3043e`, `d412750`). Verdicts recovered below.
2. The "3 COMPASS human junctions" (SQSTM1/TMED9/SLC35A4) listed as real-data grounding
   are **SPURIOUS** — reads place all three at canonical annotated GT-AG sites, not novel
   non-canonical (`dev/SNAP_OR_HOLD_ANCHOR_COORD_FINDING.md`, 2026-07-03). Deliverable B's
   premise collapsed.
3. The celebrated Phase-5 "32× recall win" was **inflated by a bookkeeping bug**;
   corrected it is **~4% vs minimap2 0.5% ≈ 8×** (Phase 6, `fc58950`).

**Member-vs-refiner (settles the brief's central question):** the native work is wired
as a **junction refiner**, NOT a 6th consensus panel member. `junction_refiner.py:648`
`refine_junctions(..., motif_blind=False, microhom_drift_margin=0.0,
drift_positional_gate=0.0)`. `multi_aligner.py` has **no native/6th member**. So the honest
answer to "does the native aligner earn *panel membership*" is: **it did not manifest as
a panel member at all — it manifested as a motif-blind refiner**, and its real-data win is
a *discovery* (recall) win, not a *placement-arbitration* member win.

---

## 1. True per-facet scoreboard

The original C1–C6 "facet" program (the −logP native aligner as a panel member). Fitness =
the ground-truth set, never the internal score. All verdicts below are from committed design
docs + committed `_result.txt` files.

| Facet | Verdict | Evidence (file:lines / commit) | Real-data-vetted? |
|---|---|---|---|
| **C1-del** (HP-length-law deletion gap cost) | **CONFIRMED — ships** | Claim-A ablation: boundary_sub placement flat 0.00→law 0.55/B0 0.78, clean false-indel ≈0 (`dev/C1_DESIGN.md:95-117`, `scripts/benchmark/c1_ablation.py`). Ships as `align_exon_block_global(penalty_table, lam, ins_lengthlaw=…)`; byte-identical when `penalty_table=None`. | **YES (safety)** — real-SIRV over-call ablation on LRGASP RNA002 + SG-NEx sub-only stratum = 0.0000 over-call (`c1_real_sirv_ablation.py`; `C1_DESIGN.md:157-178`). |
| **C1-ins** (insertion length-law) | **REFUTED → GATED OFF** | Real SIRV caught a bug: FULL law hallucinates indels 3–7%, growing with run length; `--zero-ins` isolated cause entirely to `ins_open_delta` (`C1_DESIGN.md:157-178`). Now `ins_lengthlaw=False` default + regression test `test_ins_discount_gated_off_by_default_no_hallucination`. | YES — LRGASP RNA002 + SG-NEx (the source of the refute). |
| **C1 Claim-B** (per-length SHAPE transfer) | **DEFERRED / UNVALIDATED** | Real-SIRV placement ablation is underpowered (~9 distinct HP runs ≥9; saturates at L8) AND truth-confounded ("iron triangle") → pre-committed NULL, do not build (`C1_DESIGN.md:120-155`; `aff9660` scaffold). Genuine test = held-out injection sim at SIRV-measured rates; multi-night, not done. | NO (deferred; no valid substrate). |
| **C2** (3'/poly-A CPA placement) | **REFUTED** | Shipped guarded walkback AT CEILING (err 0.00 all g) on identifiable genomic-A drift; apparent decoder "win" is a truth-definition artifact; decoder over-calls A-ending bodies (2.72) (`dev/C2_DESIGN.md:131-160`, `c2_gate.py`). No product facet. | N/A — DRS has no per-read CPA truth (`C2_DESIGN.md:170-174`); moot given refute. |
| **C3** (calibrated −logP/LLR arbiter) — the KEYSTONE | **REFUTED (both arbiters)** | Indel: ceiling==arbiter==0.999, HEADROOM 0 even on `HP_HARD-bsub` where members disagree 100% (`c3_headroom_result.txt`; `C3_DESIGN.md:141-153`). Junction: mm2 snaps 47% on DISCOVERY yet `hp_ed` picks truth 100% given a truth member, HEADROOM 0 (`c3_junction_headroom_result.txt`). C1's panel win already flows through the shipped `hp_ed` arbiter with NO LLR. | Sim-only; no real-data arm (refute stands on sim). |
| **C4** (paralog / POA-pooled window selection) | **DEFERRED-with-reason** (NOT a hard refute) | HEADROOM=0.0000; minimap2 at ceiling on identifiable (intact-base) paralog reads; entire below-ceiling gap is corrupted-/zero-evidence; smoke-(F) pooling proof truth-circular + redundant (`c4_headroom_result.txt`; `C4_DESIGN.md:102-147`). Stratum is near-zero-power for C4's real regime (2-copy/1-SNP → mapq>0 always). | NO — deferred to a measured real-data trigger (true SMN1/SMN2 ≥3-copy/mapq0). |
| **C5** (panel-failure tail / FracMinHash fallback) | **DEFERRED** (no trigger at realistic error) | Recoverable C5 tail ~0 at realistic error (0–0.39% up to 20% err; only balloons at 35–50% injected error); bulk below-ceiling is `rwwi` (C3/refiner) not C5 (`c5_tail_measure_result.txt`; `C5_DESIGN.md:39-56`). Dep-commit gate = no build. | NO — real trigger would be a measured panel-failure tail on real data; not built. |
| **C6** (variant-aware junctions) | **DEFERRED as VCF-integration / REFUTED as read-evidence member** | HEADROOM=1.0 (mm2 fabricates 180/180 variant-adjacent FP introns) BUT every driver FP intron is NM==0 == the truth D (zero-evidence trap); the only read-level recovery rule FNs length-matched real introns ~1:1; a blunt near-variant rule suppresses real junctions (recall 0.98) (`c6_headroom_result.txt`; `C6_DESIGN.md:93-155`). Read-evidence member REFUTED; catalog-gated (matched germline VCF) N→D relabel deferred to a real-data trigger. | NO — deferred to matched-VCF real-data (independent evidence). |
| **Discovery** (canonical-snap tiebreaker in `select.py`) | **OPEN — real residual found, surgical fix PROPOSED (deferred, uncommitted-to-shipped)** | `tb->snap=win=snap=0.009` residual > 0 (`discovery_tiebreak_probe_result.txt`; commit `fc3043e`). Primary recommendation: a surgical `scoring.py::_count_junction_proximity_errors` post-N `prev_rp = intron_end-1 → intron_end` fix that makes 21 snaps LOSE on primary score (zero truth-recall cost) — **DOMINATES** the −ED tiebreaker, no canonical_count interaction. NOT yet applied to shipped code. | NO — sim probe; touches shipped `scoring.py`, deferred pending real-data need. |
| **Flat-Q** (per-base quality arbitration axis) | **REFUTED (decisive NULL)** | Even ORACLE Q adds Qgain=0.0000 arbitration headroom over the error-type table across HP/HP_HARD/STR (`flatq_headroom_result.txt`; commit `d412750`). The one untested orthogonality lever is closed: the orthogonal axis is per-*read* reliability + −logP emission, NOT per-base Q. | N/A (oracle-Q upper bound already fails). |

**Facet tally:** 1 confirmed-and-ships (C1-del) · 4 refuted (C1-ins, C2, C3, Flat-Q) ·
3 deferred-with-reason (C4, C5, C6) · 1 open-proposed (Discovery surgical fix). **No C1–C6
facet produced a shipped panel MEMBER** — C1-del ships as an inner-scorer knob, not a member.

### The REOPENED program (2026-07-01 →) — the motif-blind junction-discovery re-aligner

This is the live work and is NOT captured by the C1–C6 table. Chronology (commit dates):

| Milestone | Verdict | Evidence |
|---|---|---|
| Novel-junction blindspot gate | **BUILD SIGNAL** — panel flattens 47–90% of non-canonical novels (blindspot rises monotonically GT-AG 0.017 → deep 0.90) | `novel_junction_blindspot_result.txt`; `6162419` (2026-07-01); `DIRECTOR_ALGO_EVAL_FABLE.md §0A` |
| Synthetic FDR gate (`rgfdr`) | **INCONCLUSIVE** — on synthetic yeast mapPacBio recovers non-canonical precisely (0 fabrication on intronfree); gate too benign vs real ONT | `rgfdr_result.txt`; `fc2193e` (2026-07-02). Native-aligner case rests on **real-data**, not the gate. |
| **arm-B** (motif-blind local re-placement inside panel window) | **the DISCOVERY ENGINE — real-data recall win** | `junction_refiner.py:648` `motif_blind` toggle (`69a230f`). Real recall ~4% vs mm2 0.5% ≈ **8×** (corrected; Phase 6). |
| **arm-C** (−logP inner scorer applied to junction placement) | **REFUTED — HURTS placement near HP runs; NOT salvageable in the re-placer** | `7ceed77`, `0984693` (2026-07-05). ⚠ This means "C1-del confirmed" does NOT license extending the −logP scorer to junctions. |
| Compensating-indel invariant fix | **SHIPS — always-on** | `e40ca00` (2026-07-15), `junction_refiner.py` +31 lines, `tests/test_junction_refiner.py`. Refuses a junction move that adds a matched compensating indel (read never supported the move). Reverts ~94% of moves, removes ~91% of apparent fabrication at zero real-discovery cost. |
| Homopolymer + microhomology drift GUARD | **SHELVED (real-data gate, 2026-07-17)** | Built + synthetic-validated (fab 1.31→0%; `8f50740`, `3a716aa`); but `_positional_signal` cannot separate genuine recovery (48% >0 = coin flip) from drift on real ONT data (`dev/GUARD_REEVAL_RESULT.md`; `d5b25d3`). Params default 0.0 (dormant/byte-identical). Leave OFF; candidate for future cleanup removal. |

---

## 2. What is GENUINELY still OPEN

1. **Real-data membership of the discovery re-aligner on MODERN chemistry / independent
   corpora.** The Phase-6 discovery win (~4% vs 0.5%) was measured on the Sumner SMA A549
   chr5 DRS panel with short-read corroboration. It is not yet reproduced on (a) a second
   independent DRS sample, or (b) with a clean, non-circular truth set. **Probe:** re-run
   the corrected re-placer (compensating-indel fix ON) on a second real DRS corpus with an
   orthogonal SR truth; measure fix/harm at high-confidence junctions. (Cluster; data
   already on Sherlock.)

2. **Discovery canonical-snap tiebreaker surgical fix (cheap M1-inline).** The proposed
   `scoring.py::_count_junction_proximity_errors` post-N `prev_rp` fix
   (`discovery_tiebreak_probe_result.txt`) is validated in-probe to dominate the −ED
   tiebreaker with zero truth-recall cost, but is NOT applied to shipped `select.py`/`scoring.py`.
   **Probe:** apply the one-line change behind the existing test harness, re-run the probe +
   smoke, confirm byte-identical on non-DISCOVERY strata. Distinct from the compensating-indel
   fix (different file/date).

3. **C1 Claim-B (length-SHAPE) validation** — still deferred; needs the held-out injection
   simulator (natural long-HP templates, inject length-correlated del+boundary-sub at a
   SIRV-measured *table-independent* rate). Multi-night; only earns effort if a real consumer
   needs the shape (currently none).

4. **Deferred-to-real-trigger facets (do NOT build until a measured trigger):** C4 (≥3-copy /
   mapq0-with-evidence paralog stratum — true SMN1/SMN2 on Sumner), C5 (measured panel-failure
   tail on real data), C6 (matched germline VCF as independent evidence for catalog-gated N→D).

5. **Guard-machinery cleanup (optional):** the shelved `microhom_drift_margin` /
   `drift_positional_gate` / `drift_near_tie_cap` are dormant (default 0.0). Candidate for
   removal in a simplification pass; no behavior change.

---

## 3. Real-data membership status (for the confirmed / live work)

**What real-data grounding EXISTS:**

- **C1-del inner scorer (safety):** real-SIRV over-call ablation on **LRGASP RNA002** (SIRV)
  and **SG-NEx HEYA8** (RNA00x DRS, 5,329 reads via the remote-BAM trick) — deletion-only law
  = 0.0000 over-call on the sub-only stratum, matching flat (`C1_DESIGN.md:157-178`). The
  error-model grounding (indel≥2 ≈ 0.37–0.39, low clustering) reproduces across both SIRV
  sources.
- **RNA004 modern-chemistry (error model):** **DONE** — the LongBench RNA004 H69 ~11.3 GB
  download WAS approved and processed (job 31838452 COMPLETED rc=0, 17,461 spike-in reads,
  session-9). Finding: RNA004 = clean bulk ~1% error + a ~1–2%-of-reads hot tail at ~12%
  (MAPQ 60, genuine chemistry, not mis-mapping) → the injector is NOT recalibrated toward high
  clustering; the hot tail is disregarded for facet-building per the read-quality principle
  (`NATIVE_ALIGNER_OVERVIEW.md` WS-1; `754a8ec`). **⚠ The ~11.3 GB approval gate is NOT
  pending — it was granted and consumed.**
- **Discovery re-aligner (recall win):** **Sumner SMA A549 chr5 direct-RNA (Nanopore)** vs
  GENCODE, with independent **short-read** corroboration. Corrected real recall ~4% vs mm2
  0.5% ≈ 8×; compensating-indel fix transfers with ~91% fabrication removed at zero discovery
  cost (Phase 6, `fc58950`).

**⚠ CONTRADICTION with the brief — the 3 COMPASS anchors are SPURIOUS.** The brief lists "the
3 COMPASS human junctions / Deliverable B" (SQSTM1/TMED9/SLC35A4) as real-data grounding.
`dev/SNAP_OR_HOLD_ANCHOR_COORD_FINDING.md` (2026-07-03) shows both COMPASS short reads and
SG-NEx long reads place all three at **canonical, annotated GT-AG sites 175–1017 bp away** from
the harness "non-canonical anchor" coordinate, which has ~0 read support. The anchors are not
novel and not non-canonical — Deliverable B's premise collapsed. Do NOT treat these as a
surviving real-data membership test.

**What real-data question REMAINS:**
- **(a) Modern-chemistry RNA004 for DISCOVERY (not just the error model):** the error-model
  characterization is done, but the discovery re-aligner's recall/precision has not been
  measured on a clean RNA004 discovery corpus with independent truth. (Sumner DRS is the
  modern real substrate actually used; a second independent RNA004 sample is the open check.)
- **(b) Human-data transfer / SMN biology:** the corrected engine is now cleared to be pointed
  at the actual SMA splicing science on Sumner (SMN1/SMN2 chr5) — this is the stated next prize
  (`fc58950`, `2b8d2ed`), and doubles as the C4 (paralog) real-data trigger. No clean-truth
  human discovery-precision number exists yet.

---

## 4. What DIRECTOR_ALGO_EVAL_FABLE* is

A set of three **independent DIRECTOR-level algorithmic evaluations** (2026-07-01) commissioned
when the program reopened, mandated to (1) tear down the 5-aligner long-read panel for
NOVEL-ISOFORM discovery and (2) rank orthogonal, gate-testable algorithmic approaches. Despite
the "FABLE" filenames, a provenance check proved all three actually ran on **Opus 4.8** (the
harness declared Fable 5 but the transcripts show the `claude-fable-5→claude-opus-4-8` fallback;
`DIRECTOR_ALGO_EVAL_SYNTHESIS.md §0`) — so they are three independent *same-model* runs, not a
cross-model check. **Conclusion (converged across all three, `SYNTHESIS.md §1`):** the panel's
blind spot is a *scoring* convergence (all five members share a flat-affine, quality-blind,
GT-AG/annotation-scored emission model) that flattens non-canonical novel junctions (blindspot
rises to 0.90); and the single most promising orthogonal fix is **motif-blind, empirical-−logP
local re-alignment inside the panel-defined window** (the C1/C3 engine repointed from HP-indel/
arbitration to junction discovery, GT-AG/annotation bonus zeroed). That recommendation is exactly
what became the reopened program's arm-B re-aligner.

---

## 5. Recommended next increments (ranked, falsifiable)

Framed to the real question: **"does the native (motif-blind) re-aligner earn its place on REAL
data as a discovery engine?"** (Panel-*membership* is the wrong frame — it manifested as a refiner.)

1. **[CHEAP · M1-inline] Apply + gate the Discovery surgical `scoring.py` fix.** Apply the
   `_count_junction_proximity_errors` post-N `prev_rp` correction, re-run
   `discovery_tiebreak_probe.py` + `smoke_roundtrip.py`, and assert byte-identical on
   non-DISCOVERY strata. **Falsifiable:** snap-FDR must drop to ~0 with zero truth-recall loss;
   if any non-DISCOVERY stratum changes, revert. Low-risk, already validated in-probe, closes
   the one OPEN sim item. (Distinct from the shipped compensating-indel fix.)

2. **[REAL-DATA · cluster] Reproduce the discovery recall/precision win on a SECOND independent
   real DRS corpus.** Run the corrected re-placer (compensating-indel fix ON, guard OFF) on a
   second real DRS sample (a non-Sumner DRS, or a held-out Sumner sample) with orthogonal
   short-read truth; measure high-confidence fix/harm and corrected recall. **Falsifiable:** the
   ~8× recall win and ~2× high-confidence precision must reproduce; if recall collapses to ≈mm2
   or fabrication returns, the Sumner win was corpus-specific. This is THE membership question.
   Data largely staged on Sherlock; no new download.

3. **[REAL-DATA · science, cluster] Point the corrected engine at SMN1/SMN2 splicing on Sumner
   chr5** (the stated prize, `2b8d2ed`). Doubles as the C4 paralog real-data trigger (true
   ≥2-copy/mapq0 locus). **Falsifiable:** does motif-blind re-placement recover SMN2 exon-7-skip
   / cryptic junctions the panel flattens, corroborated by short reads? A null here says the
   discovery win doesn't transfer to the paralog-hard regime it was ultimately built for.

Do NOT: extend the −logP inner scorer to junction placement (arm-C REFUTED — hurts near HP
runs); re-enable the shelved drift guard; rebuild C2/C3/Flat-Q (decisively refuted); or build
C4/C5/C6 members absent a measured real-data trigger.

---

## 6. Smoke / regression state

**GREEN.** Ran the read-only smoke from the worktree:
`PATH=… python scripts/benchmark/smoke_roundtrip.py --out /tmp/smoke_audit --reps 20` →
**exit code 0**, "SMOKE PASSED — all GATE assertions green" (all A/B/B2 junction-scoring gate
assertions pass). Note: the run took **>2 min** wall (not the "quick" the brief assumed) — it
completed cleanly in the background. No FAIL. The broader `pytest -m "not slow"` suite is
reported green on this branch per `CLAUDE.md` (not re-run here — read-only audit).
