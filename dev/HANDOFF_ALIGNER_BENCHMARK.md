# HANDOFF — Deliverable A: the simulation ground-truth benchmark (the GATE)

**Agent:** dedicated benchmark-builder (isolated worktree, branch
`worktree-agent-a25a2c1e784ad37dc`, based on `drs-validation-rebuild` so the reuse
primitives + design docs are present). **NEVER commit to `drs-validation-rebuild`.**
**Updated:** 2026-06-30 (session 12 — fan-out recovery).

---

## SESSION-12 (2026-06-30) — 8-AGENT FAN-OUT hit API instability; RECOVERED to reliable inline

The end-of-SESSION-11 8-agent parallel fan-out (C4/C5/C6/Discovery + Deliverable-B/
Flat-Q/Claim-B/WS-1) ran into a **systemic API/streaming instability**: 8 heavy
concurrent agents, each directed to spawn its own triple adversarial panel,
over-subscribed concurrency → most agents did real work (17–38 tool-rounds) but lost
their FINAL report to "connection closed mid-response"; two (C5, Flat-Q) failed early
at the 600s stall-watchdog. **LESSON: do not run >2–3 heavy agents concurrently, and
do not have them spawn nested panels. Prefer director-run INLINE gates** (the
c3/c4/c5 method never failed) + at most 1–2 controlled agents.

**DONE + COMMITTED this session (recovered from worktrees + re-run inline):**
- **C4 paralog/POA = DEFERRED** (commit c4e8207). HEADROOM=0.0000; minimap2 at
  ceiling on identifiable paralog reads; below-ceiling gap is all zero-evidence; smoke
  (F) pooling proof is truth-circular + redundant. Cross-checked by inline re-run.
  `dev/C4_DESIGN.md` + `scripts/benchmark/c4_headroom.py` (+result).
- **C5 panel-failure tail = DEFERRED** (commit 658b72c). Recovered the agent's
  `c5_tail_measure.py` (fixed 1-line bug candidates→cand_meta), ran Stage-1: C5 tail
  ~0–0.4% at realistic error (rate≤0.20), only 17–24% at EXTREME injected rates
  (0.35–0.50). Below-ceiling reads are overwhelmingly C3-refiner (rwwi) + C4 (paralog
  wrong-window), NOT C5. `dev/C5_DESIGN.md` (+result). Stage-2 containment-vs-null
  positive control NOT yet captured (the --stage2 run is heavy/slow).

**LOST (agents stalled, worktrees auto-cleaned / empty) → RE-RUN INLINE:**
- **C6 variant-aware** (worktree had nothing). **Discovery tiebreaker** (worktree gone).
  **Flat-Q** (failed early, gone). All three are M1 gates — re-run inline.

**BLOCKED on Sherlock auth (DOWN: `ssh sherlock` = Permission denied gssapi):**
- **Deliverable B real-data** (3 COMPASS junctions) — agent ran on cluster; M1
  worktree auto-cleaned; outputs (if any) at `/scratch/users/kevinroy/c1_realdata_dB/`.
- **WS-1 cleanup** — outputs at `/scratch/users/kevinroy/combined_ref/ws1_cleanup/`.
- **Claim-B** — needs the SIRV-measured rate from Sherlock; nothing recovered on M1.
- ACTION NEEDED FROM KEVIN: re-auth via `! ssh sherlock true` (do NOT re-open the
  master yourself). THEN check those two scratch dirs.

### RESUME (concrete, branch logic) — from worktree `worktree-agent-a25a2c1e784ad37dc`
Env prefix for all: `PATH="/Users/kevinroy/miniconda3/bin:/opt/homebrew/bin:$PATH" PYTHONPATH=. /Users/kevinroy/miniconda3/envs/pysam/bin/python`
- **C6 / Discovery / Flat-Q (do inline now):** write probes mirroring
  `scripts/benchmark/c3_headroom.py` (the exemplar) + `c4_headroom.py`. C6 → VARIANT
  stratum + scorer `junction.fp_variant_adjacent` (smoke E); strong prior it's
  variant-adjacent-FP-addressable OR zero-evidence. Discovery → drive the real
  `select.py::select_best_alignment` on JUNCTION_DISCOVERY {mm2-snap, truth-site};
  strong C3 prior = tiebreaker never decides (snap loses on primary score). Flat-Q →
  extend `error_injector.py` to emit correlated Q, headroom of Q-aware vs error-type.
  Commit each as `feat(aligner-bench): C{n} ... ` + re-run smoke (must stay exit 0).
- **C5 Stage-2 (optional control):** `… scripts/benchmark/c5_tail_measure.py --out
  /tmp/c5 --reps 60 --stage2` (BACKGROUND — heavy 42930-window sketch; ~>2min). Append
  result to `c5_tail_measure_result.txt` + the C5_DESIGN residue. Does not change DEFER.
- **After `ssh sherlock true` re-auth:** `ssh sherlock 'ls -la /scratch/users/kevinroy/c1_realdata_dB/ /scratch/users/kevinroy/combined_ref/ws1_cleanup/'`
  → if Deliverable-B outputs exist, recover the per-junction snap-vs-hold verdicts;
  else re-run the realign (single controlled agent or inline). WS-1 similar.
- **Gate scoreboard so far:** C1 del-law CONFIRMED; C1-ins / SIRV-ClaimB / C2 / C3 /
  C4 / C5 refuted-or-deferred. Remaining live: C6, Discovery, Flat-Q (inline);
  Deliverable-B, WS-1, Claim-B (cluster, blocked).

---

## SESSION-11 (2026-06-29) — C3 (calibrated LLR arbitration, the KEYSTONE) → GATE-REFUTED

C3's premise: replace the integer-max consensus ARBITER with a calibrated −logP LLR
arbiter (refiner emits posterior + runner-up; consensus compares paths by likelihood
ratio, not integer-max). The named justification is the **0.09→1.07 artifact**. Ran
the C1/C2 gate-first discipline (advisor-shaped, two advisor passes). **VERDICT:
C3-as-accuracy REFUTED** — no member code written (the C1-idiom refute). Now **4
directions refuted/deferred (C1-ins, SIRV Claim-B, C2, C3), 1 confirmed (C1-del)**.

**What was pinned (both incumbent arbiters, firsthand):**
- `select.py::select_best_alignment` = `max(junction_score)` (flat −2/−1 penalties,
  quality-BLIND, blind to in-exon indel placement) + annotated/canonical tiebreakers
  → **the literal 0.09→1.07 home** (a penalty/tiebreak re-weight flips the winner).
- `corrected_consensus.py::merge_corrected_tsvs` (`use_hp_ed`) = `argmin
  hp_edit_distance` + span tiebreak → **the GOVERNING arbiter** in the canonical
  correct-first pipeline (`write_corrected_consensus_bam` emits its winner's full
  record, junctions included). `hp_ed` is HP-aware on indels (via the penalty table)
  but flat on mismatch (X=1.0), clips, introns FREE.

**The gate (3 M1-light, LLR-free probes; fitness=truth; ALL show NO addressable gap):**
- **Indel arbitration** (`scripts/benchmark/c3_headroom.py`, panel=[flat,law,mm2],
  1600 reads): ceiling==arbiter==0.999, HEADROOM=0.000. STRONG null — on
  `HP_HARD-bsub` members disagree on **100%** of reads (`disagr`=1.0) yet `hp_ed`
  picks truth every time (`hr|dis`=0). **C1's panel-level win already flows through
  the shipped arbiter with no LLR** → C3's keystone premise ("C1 needs C3 to land") is
  empirically WEAKENED.
- **Junction arbitration** (`c3_junction_headroom.py`, members=[mm2 real snap, truth
  site], the governing `hp_ed`, 210 reads): on JUNCTION_DISCOVERY minimap2 **snaps 47%**
  yet given a truth member `hp_ed` picks the snap over truth **0.000** (introns free +
  snap induces flanking mismatches). HEADROOM=0.000.
- **Multi-event coherence** (`c3_multievent_check.py`): `penalty_score`-sum DOES rank
  multi-event paths oppositely to −logP-sum (89/55-choose-2 inversions — confirms C1's
  "incoherent to sum"), BUT they're UNREACHABLE as a truth-favoring tie (cross-run-length
  reassignable ties don't arise from same-base HP ambiguity or different-base adjacent
  runs). The **named locus** for any future revisit = a dedicated adjacent/interleaved-runs
  multi-event stratum; until built, do not build the LLR.

**Structural closures (not assertions):** the `junction_score` path is closed —
`_count_junction_proximity_errors` penalizes a snap on the PRIMARY score (fires on
84/360 corpus junction reads), so the canonical tiebreaker only bites on ties a snap
won't produce → a one-line tiebreaker reweight (Discovery/WS-3), not an LLR. The
artifact-replay (re-weight → integer flips, LLR invariant) is near-tautological → a CI
fence only, earns its place ONLY if the LLR ships for an accuracy reason (it won't).

**Files:** `dev/C3_DESIGN.md` (full gate + verdict + residue), `scripts/benchmark/c3_{headroom,junction_headroom,multievent_check}.py` (+ `*_result.txt`). Smoke gate GREEN.
RESUME: C3 is CLOSED (refuted). Next facet candidates per OVERVIEW: C4 (POA paralog —
stratum exists + discriminates) or C5 (panel-failure-tail FracMinHash — gated on a
MEASURED tail; the JUNCTION_DISCOVERY all-herd/snap case is C5 territory, not C3). The
flat-Q quality-arbitration axis (SPEC:225) remains the one untested orthogonality lever.

**IN FLIGHT (launched end of SESSION-11): an 8-agent parallel fan-out** — all
worktree-ISOLATED background agents, each running the C3 gate-first playbook,
"report back, don't commit." Each was also sent a director directive (via
SendMessage) to use the `advisor` LIBERALLY/FREQUENTLY and run a TRIPLE adversarial
subagent panel before finalizing any verdict. Roster:
- **C4 paralog/POA** (a317eb683ef499cb1): identifiable-but-misplaced population pooling recovers? (excl. zero-evidence fragments + minority-collapse trap).
- **C5 panel-failure tail** (ab1fc5de61e2c39e3): MEASURE the tail first (empty-union + all-misplaced, recoverable-vs-zero-evidence) at baseline & elevated error via `score_panel`; tail size = dep-commit trigger; FracMinHash-vs-random-window-null is the Stage-2 kill-gate.
- **C6 variant-aware junctions** (a2c9de0d12c209649): does variant-adjacent FP fabrication (smoke E) survive as an addressable, specific gap?
- **Discovery canonical-snap tiebreaker** (ab95f6701cc1b7d7e): does `select.py`'s `canonical_count` tiebreaker EVER decide a snap-vs-truth tie (C3 predicts ~never; snaps lose on the primary score)? reweight only if it harms — touches shipped code, propose diff.
- **Deliverable B real-data** (a9e9b7a45009c3a3f, CLUSTER): C1 realign on the 3 COMPASS human junctions (SQSTM1/TMED9/SLC35A4) — snap-to-canonical vs hold-non-canonical. Own scratch `/scratch/users/kevinroy/c1_realdata_dB/`, read-only on `compass_a549`.
- **Flat-Q quality axis** (adfc06a1d85fc50b5): does correlated per-base Q open any arbitration/discovery gap the error-type table doesn't? (the one untested orthogonality lever).
- **C1 Claim-B injection simulator** (a481726d77ac06d9a, MULTI-NIGHT): held-out injection at SIRV-measured (table-INDEPENDENT) rates → does the Scer length-SHAPE transfer?
- **WS-1 cleanup** (a18b9947050c10219, CLUSTER): debug the empty `measure-bam` on `h69_endo_1.bam` (contig 1) + the qscore SIGPIPE; endogenous-vs-SIRV contrast. Low value/cosmetic.

RESUME (concrete, branch logic):
- **Same session** (the agents report into this session via task-notifications):
  read each agent's final message → REVIEW its VERDICT + the full-text scripts it
  returned → in THIS worktree (`worktree-agent-a25a2c1e784ad37dc`) write
  `dev/C{4,5}_DESIGN.md` + `scripts/benchmark/c{4,5}_*.py` from that text → run the
  regression gate `PATH="/Users/kevinroy/miniconda3/bin:/opt/homebrew/bin:$PATH"
  PYTHONPATH=. /Users/kevinroy/miniconda3/envs/pysam/bin/python
  scripts/benchmark/smoke_roundtrip.py --out /tmp/x --reps 20` (exit 0) → surgical
  `git add` those explicit paths → commit. If a verdict is PROCEED, the next step is
  to SKETCH (not build) the member design; if REFUTE/DEFER, document and stop.
- **Fresh session** (agent IDs a317eb683ef499cb1 / ab1fc5de61e2c39e3 are NOT
  resumable outside this session): run `git worktree list` to find the two agents'
  isolated worktrees and inspect their `dev/C{4,5}_DESIGN.md` + `scripts/benchmark/`
  in-progress files; integrate the durable ones as above. If their worktrees were
  auto-cleaned (unchanged), the gates are fully re-runnable from the durable specs —
  the C3 exemplar (`dev/C3_DESIGN.md` + `scripts/benchmark/c3_headroom.py`), the
  SPEC §benchmark-coupling C4/C5 blocks, and the OVERVIEW facet table — re-launch the
  two gate investigations from those. Neither agent has authority to commit to this
  branch; the director integrates + re-runs smoke before any commit.

---

## SESSION-10 (2026-06-29) — PARALLEL FAN-OUT (advisor-shaped: aim, not breadth)

Three workstreams launched as background agents (directing agent integrates; all "report back, don't
commit"). Advisor reshaped the plan: WS-1 GATES WS-2 (pre-competitive BAMs' error tail is partly
contamination); bimodality must use reference-FREE QSCORE not error-rate (a 2-Gaussian on skewed count
data fakes bimodality); and the real aim is WS-1 (settle RNA004) + WS-3 (guardrail) + C2 (the actual next
FACET — C1 is consolidated so C2 is back on).
- **WS-1 — competitive alignment + qscore + endogenous contrast. (Agent ad7cf26d06de18d55 DIED on an API
  stall before doing cluster work → re-encoded as a SELF-CONTAINED SBATCH: job 31922545,
  `scripts/benchmark/run_competitive_fourway.sbatch`.)** Downloads the COMBINED hg38+SIRV ref, competitively
  re-aligns RNA004(H69)+LRGASP-ENCODE, measures the SIRVomeERCCome spike-in subset (FAIR four-way) AND a
  chr1 endogenous-human subset (biological-excess contrast), + qscore distributions. Branch on the sentinel:
  `ssh sherlock "sacct -j 31922545 -X -n -o State; cat /scratch/users/kevinroy/combined_ref/.comp_rc 2>/dev/null"`
  → if 0, read `/scratch/users/kevinroy/combined_ref/{competitive_fourway,qscore_dist}.report.txt`.
  **RESULT (2026-06-29): job ended rc=141 = SIGPIPE on the cosmetic qscore step ONLY; the alignment +
  four-way measure COMPLETED.** Headline: **RNA004 divergence is NOT mis-mapping contamination** — the
  competitive SIRVomeERCCome count is 17,449 vs the old SIRV-only 17,461 (−0.07%) and the burstiness
  PERSISTS identically (overdisp_v 2.03 / gap5x 9.34); RNA002 (encff) stays low (0.094 / 1.66). So the
  advisor's contamination-asymmetry hypothesis is REFUTED → the live explanation is a genuine HOT-READ
  SUBPOPULATION (dispersion 17.5 + p90/median 2.58). **Two recoverable gaps (BAMs persist on scratch):**
  (a) the endogenous contrast got 0 reads — the combined ref names human contigs `1`..`22`/`MT` (Ensembl),
  NOT `chr1`; re-run on contig `1`. (b) qscore step SIGPIPE'd; re-run with `head`/pipefail-off. **The
  settling diagnostic = the per-read hot-tail decomposition on `h69_spk.bam`** (pull top 1-2% by error
  rate; if removing them collapses dispersion toward RNA002 AND they are low-mapq → hot/low-Q subpop, not
  broad chemistry). **DECOMPOSITION DONE → RESOLVED:** RNA004 = clean bulk ~1% error + a ~1-2% tail at ~12%
  (13× bulk) at MAPQ 60 (genuine hot reads, not mis-mapped); removing the tail → ~1% bulk. So the divergence
  is a bimodal HOT-READ TAIL, NOT chemistry. **DECISION: do NOT recalibrate the injector toward high
  clustering**; the decent bulk is low-clustering (matches RNA002). The ~1-2% tail is disregarded for facet
  building (per the read-quality principle) — it's QC/flow-cell characterization + WS-3 guardrail input.
  Cosmetic-recoverable: (a) endogenous biological-excess contrast — contig `1` extraction CONFIRMED
  (h69_endo_1.bam, 517,487 primary reads; the `chr1`→`1` fix is done), but `measure-bam` on it produced
  NO output → debug WITHOUT `2>/dev/null` (likely a per-read error on a big spliced BAM, or the [600,1000]
  exonic window excludes most; the BAM persists at /scratch/users/kevinroy/combined_ref/h69_endo_1.bam).
  (b) qscore split (SIGPIPE → use `head`). **WS-1 effectively COMPLETE; clean fan-out boundary → reset
  director for C3.**
- **WS-3 — novel-feature support guardrail prototype: DONE + COMMITTED** (`scripts/benchmark/novel_support_probe.py`).
  TE = −log10 P(Binom(n,p0)≥k) per novel call (support-count-aware; truth-blind metric). AUC 0.998 (gt) /
  0.93-0.96 (reference-free exonic-density proxy); FDR 0.50→0.20 (gt) / 0.50→0.33 (proxy); recall ~0.99;
  seed-stable; advisor-vetted; reproduced locally. Honest caveats baked in (soft prior; magnitude
  SIRV-gated; separation partly by construction = mechanism demo; real-data risk = per-molecule-vs-
  per-transcript hotness confound). Feeds C3 / discovery-FDR. NOT a shipped feature — a vetted prototype.
- **C2 — DONE: REFUTED as a placement facet + COMMITTED** (`scripts/benchmark/c2_gate.py`,
  `dev/C2_DESIGN.md`, GENOMIC_A_CPA stratum in `controlled.py`). Shipped guarded walkback already at ceiling
  (0.00 |est−true_cpa|) on the identifiable genomic-A drift; read-through NULL holds; the apparent win is a
  truth-definition artifact (gate COMPUTES the truth-flip → arms swap); decoder over-calls A-ending bodies
  (2.72). NO product facet. **Regression caught+fixed:** agent wired `include_cpa=True` BEFORE variant →
  consumed shared RNG → broke smoke (E); I made it default-OFF + append-only LAST (smoke GREEN; c2_gate
  reproduces). **C2 promotes C3:** the CPA uncertainty (`ambiguity_range`, bam_processor.py:851-859) ships
  with NO consumer → **C3 (LLR arbitration) is now THE active next facet**; soft CPA posterior deferred behind it.
RESUME: each agent reports back; check its output, review, then I commit. If an agent died, its task is
re-runnable from its prompt (the specs are durable: WS-3→OVERVIEW §novel-feature; C2→C1_DESIGN playbook +
SPEC GENOMIC_A_CPA; WS-1→SPEC §FOUR-WAY UPDATE diagnostics). Combined ref:
`https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/annotations/genome_fasta/hg38_sequins_SIRV_ERCCs_longSIRVs.fa`.

---

## NIGHT PUSH (2026-06-29, session-9) — N-in-span FIXED; C1 de-novo-aligner facet BUILT+TESTED+PROVEN; SG-NEx grounding

Autonomous overnight push (user: "build + test the de novo aligner by morning; spare no expense;
adversarial panels on ambiguity"). All committed; smoke GREEN throughout.

- **N-in-span measurement bug FIXED** (RESUME #1 blocker). `events_from_alignment(..., exonic_coords=True)`
  excludes introns from span + exon-local positions; re-measured EXONIC both sides (the bug inflated
  yeast clustering 5.28→4.28). SPEC §"CORRECTED (session-7b)".
- **C1 facet (HP-length-law gap cost) IMPLEMENTED + unit-tested + ABLATION-PROVEN (Claim A).** Design
  LOCKED by a triple adversarial panel + multiple advisor passes (`dev/C1_DESIGN.md`). Cost = `rate_mean`
  baseline-anchored LOG-ODDS delta `λ·ln(rate(hp)/rate(1))` on gap-OPEN (NOT the `penalty_score` column —
  reviewer-1 proved it is reciprocal-rate, not −logP, incoherent to sum in a DP); `penalty_table=None` is
  byte-identical (Cat3 guard); gated on `homo_mask`. **Ablation (reps=120, TEST split, 3 matched arms):
  HP_HARD-noisy concordance flat 0.962 → B0 0.990 / law 0.985; boundary_sub flat 0.000 → B0 0.78 / law
  0.55; clean false-indel-rate 0.000 all arms (1854 clean reads). PASS Claim A.** HEADLINE HAND-VERIFIED
  (advisor-flagged the flat=0.000 could be a scorer artifact): flat places the deletion OUT-of-run (run
  [80,90)→D@90), law IN-run (D@80), e.g. `hph_A_10`: flat `[90M1D79M]`→law `[80M1D89M]` — real placement
  fix, not a metric quirk. **B0 ≥ law as pre-committed** (length-SHAPE mildly anti-helpful on flat-in-L —
  expected, NOT a failure; do NOT raise λ to "fix" it = the sim-overfit trap). SHIP the law on
  principle+real-data-deferred grounds, NOT sim perf (B0's constant is itself the mean of law's deltas =
  a shape-ablation control, not an independent alternative). Claim B (length-SHAPE) next test = the
  real-SIRV PLACEMENT ablation (NOT per-HP-rate-on-reads; see dev/C1_DESIGN.md). Files:
  `hp_penalty.py`, `local_aligner.py`, `tests/test_c1_lengthlaw.py` (4 pass), `scripts/benchmark/c1_ablation.py`.
  (reps=400 confirmatory run in flight.)
- **THREE-WAY real spike-in grounding (exonic, thin 0.0153, [600,1000]):** two INDEPENDENT SIRV sources
  AGREE — LRGASP-SIRV (RNA002) overdisp_v 0.10 / gap5x 1.71 / indel≥2 0.37; **SG-NEx HEYA8 (RNA00x,
  measured via the remote-BAM trick, 5,329 reads) overdisp_v 0.17 / gap5x 1.76 / indel≥2 0.38 /
  autocorr 0.62** — vs yeast read-vs-genome 0.54 / 4.28 / 0.39 / 0.47. The mod-free SIRV clean clustering
  is genuinely LOW and reproduces across labs/chemistries; `indel≥2 ≈ 0.37-0.39` is rock-solid. Caveats:
  both SIRV share IVT composition; autocorr likely per-transcript-confounded (few SIRV transcripts); RNA004
  (LongBench, GB download, needs approval) still the modern-chemistry check.

### FOUR-WAY DONE (2026-06-29, session-9 cont.) — RNA004 DIVERGES (key result)
- **LongBench RNA004 arm DONE — job 31838452 COMPLETED rc=0** (17,461 spike-in reads). It did NOT confirm
  the "clean clustering is LOW" story — RNA004 is lower-rate (0.0106) but shows a much higher dispersion
  number (overdisp_v 2.03, gap5x 9.34). **Likely NOT chemistry — a HOT-READ-SUBPOPULATION signature**
  (advisor): dispersion_index 17.5 with p90/median only 2.58 ⇒ a small ~1-2% extreme tail, not homogeneous
  burstiness. **Do NOT pivot the injector calibration** (RNA002-low vs RNA004-high stays OPEN; rate-matching
  won't resolve a mixture). **#1 OPEN QUESTION — next-session diagnostics (fresh agent):** (a) per-read rate
  histogram on the RNA004 [600,1000] set, pull top 1-2%, check MAPQ/identity/transcript (low → artifact);
  (b) UNIFY the alignment pipeline — SG-NEx used competitive hg38+spike-in; the SIRV-only sources can
  mis-map human reads → re-align ALL to a combined human+SIRV ref, take SIRV-contig reads, before any fair
  four-way magnitude comparison. Full detail: SPEC §"FOUR-WAY UPDATE".
- **Claim B coordination DONE → DECISION recorded (dev/C1_DESIGN.md §"Claim B — VETTED CONCLUSION").**
  Partner agent + 3-way adversarial fan-out + 2 advisor passes converged: the real-SIRV PLACEMENT/shape
  ablation is **underpowered** (only ~9 distinct HP runs ≥9 in SIRV; law delta saturates at L8) AND
  **truth-confounded** (the "iron triangle": discriminating ⟺ boundary mismatch ⟺ per-read truth
  unknowable) → a shape "win" would be CONFOUNDED, so it is a PRE-COMMITTED NULL; do NOT build it.
  **Genuine Claim B = a held-out INJECTION simulator** (natural long-HP templates; inject length-correlated
  del + boundary-sub at a SIRV-MEASURED, table-independent rate; injected=known truth) — DEFERRED, multi-night.
- **Over-call deliverable BUILT + RUN → caught a real bug (`scripts/benchmark/c1_real_sirv_ablation.py`).**
  On the honest sub-only stratum (length-preserving windows w/ mismatches → any indel = true hallucination),
  replicated on BOTH real SIRV BAMs: flat 0% / B0 0% / **FULL law 3-7%** (grows with run length); `--zero-ins`
  isolated the cause ENTIRELY to `ins_open_delta`. Reproduced locally (`flat [28M]` → full-law `[10M1D4M1I13M]`
  on a 1-sub A8 window). **FIX: `ins_open_delta` GATED OFF by default** (`align_exon_block_global(...,
  ins_lengthlaw=False)`); the deletion-only law is SAFE (0% over-call = flat) AND Claim-A-proven → that ships.
  Re-run with the gated default CONFIRMS PASS (SG-NEx sub-only 0.0000; LRGASP + sim-ablation re-run finalizing).
  Regression test `test_ins_discount_gated_off_by_default_no_hallucination`. Real data caught what the sim
  could not — the grounding working as intended. (C1_DESIGN.md §"OVER-CALL RESULT".)

### RESUME (night-push boundary) — concrete next steps
- **C1 ablation CONFIRMED at reps=400** (cell-size floor): HP_HARD-noisy flat 0.965→B0 0.994/law 0.985;
  boundary_sub flat 0.000→B0 0.812/law 0.572; clean false_indel 0.000 (6207 reads). PASS Claim A — same
  pattern as reps=120, hand-verified. Re-run any time: `PATH="/Users/kevinroy/miniconda3/bin:$PATH"
  PYTHONPATH=. /Users/kevinroy/miniconda3/envs/pysam/bin/python scripts/benchmark/c1_ablation.py --out
  /tmp/c1_abl_400 --reps 400 --penalty-table rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores.tsv`.
- **NEXT increment (deepen C1, do NOT widen to C2 — advisor): the real-SIRV PLACEMENT ablation (Claim B).**
  Spec in `dev/C1_DESIGN.md` §"Claim B". Inputs staged on Sherlock scratch
  `/scratch/users/kevinroy/sirv_work/{sirv.sorted.bam, sgnex_heya8_spikein.bam, sirv4.fasta}`. Enumerate
  HP-run spans from `sirv4.fasta`, re-align each real SIRV read's ref segment with flat/B0/law, score
  in-run vs out-of-run placement per HP run (`net_indel_in_span`). Multi-night; scope deliberately (no
  boundary_sub analog in real reads; co-occurring errors). This is where the length-SHAPE is earned/refuted.
- **LongBench RNA004 arm** (modern-chemistry confirmation): NEEDS USER APPROVAL for the ~11.3 GB H69
  download (SESSION-8 PLAN B). SG-NEx PLAN A already done (the HEYA8 measurement above).
- **Regression gate (any session):** `pysam-python scripts/benchmark/smoke_roundtrip.py --out /tmp/x
  --reps 20` → exit 0; unit runners `python tests/test_{c1_lengthlaw,error_injector,gff_panel_gtf,lrgasp_truth}.py`.

---

## SESSION-8 (2026-06-29) — REAL DRS SPIKE-IN RECON (SG-NEx RNA00x + LongBench RNA004); refs STAGED, two plans READY

Reconnoitered the two external-validity DRS spike-in corpora (RESUME #2 of session-7:
"confirm on RNA004 SIRV"; session-4b §EXTERNAL-VALIDITY two-chemistry plan). **No code
touched, no job submitted, no multi-GB read pulled.** SG-NEx spike-in reference + GTF are
STAGED on scratch (tiny); LongBench spike-in presence CONFIRMED for free (no read download).
Both plans are ready to submit; the LongBench one needs explicit GB approval (flagged below).

### Access facts (both buckets)
- **No `aws`/`aws2` CLI on Sherlock; no module.** `rclone v1.55.1` present but unneeded.
  Both buckets are **public and list/GET-able anonymously over the S3 REST https endpoint
  with `wget`** (Sherlock curl 7.29 still errors on https — use wget, as the CLAUDE.md says):
  - SG-NEx: `https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/` (ListObjectsV2:
    `?list-type=2&prefix=…&delimiter=/&max-keys=…`; paginate via `NextContinuationToken`).
  - LongBench: `https://longbench-data.s3.ap-southeast-2.amazonaws.com/`.
- **`samtools 1.23.1` in the `pbsim3` env has `libcurl=yes S3=yes`** → it reads **remote
  indexed BAMs directly** (`samtools idxstats/view/faidx <https-url>` does range-GETs against
  the object + its `.bai`/`.fai`). This is the key lever for SG-NEx (below).
- Parse gotcha: `tr '<' '\n'` splits `<Size>123</Size>` into `Size>123` AND `/Size>` — grep
  `'^(Key|Size)>'` (anchored) to drop the closing tags, then `paste - -` (Key precedes Size
  in each `<Contents>`). wget `--header=Range:` gave a 0-byte file here — use `samtools faidx`
  for ref subsetting instead.

### SG-NEx (RNA00x DRS) — spike-in ref STAGED; recommended first sample chosen
- **Combined genome+spike-in FASTA** (the "ready-made" ref): key
  `data/annotations/genome_fasta/hg38_sequins_SIRV_ERCCs_longSIRVs.fa` = **3,110,869,710 B
  (2.9 GB)** + `.fa.fai` (7,728 B). It bundles full hg38; **only 2 contigs are spike-in**,
  both at the TOP of the file: **`chrIS`** (10,567,884 bp = Sequins decoy chr) and
  **`SIRVomeERCCome`** (538,365 bp = SIRV + ERCC). GTF options:
  `…/gtf_file/hg38_sequins_SIRV_ERCCs_longSIRVs_corrected.gtf` (159 MB) and
  `…_v5_reformatted.gtf` (264 MB).
- **STAGED (durable, on scratch) — no full 2.9 GB download needed:**
  - `/scratch/users/kevinroy/sirv_work/sgnex_spikein.fa` (**11.1 MB**, chrIS+SIRVomeERCCome)
    + `.fai` — pulled via `samtools faidx <https-url> chrIS SIRVomeERCCome` in **3.8 s**.
  - `/scratch/users/kevinroy/sirv_work/sgnex_spikein.gtf` (**198 KB**; 165 Sequin + 176
    SIRV/ERCC transcripts) — `wget -qO- <corrected.gtf> | awk '$1=="chrIS"||$1=="SIRVomeERCCome"'`
    (streamed the 159 MB, stored only the 198 KB subset).
- **SG-NEx's OWN alignment = `minimap2 -ax splice -k14 -uf`** (from the BAM `@PG`) — IDENTICAL
  to our sirv/yeast pipeline flags → re-alignment is apples-to-apples (their BAMs are even
  directly reusable). minimap2 VN 2.17.
- **Which samples are spiked** (swept all 44 directRNA genome BAMs via remote `idxstats`,
  ~2 s each — reads only the `.bai`+header): the spiked cell lines are **HEYA8 and H9**
  (plus trace 1–583-read contamination in some Hct116/K562/MCF7/Hek293T). Top by SIRV+ERCC
  read count: **HEYA8 rep3_run1 = 15,855 `SIRVomeERCCome` + 197 `chrIS`**; rep2_run1 (14,407
  total spike-in); rep1_run1 (12,205 extracted). H9 rep2–4: 810–2,947 each.
- **CHOSEN FIRST SAMPLE: `SGNex_HEYA8_directRNA_replicate3_run1`** — most SIRV+ERCC reads
  (15,855; comparable to the yeast `sirv.sorted.bam`'s 4,424). Genome BAM key:
  `data/sequencing_data_ont/bam/genome/SGNex_HEYA8_directRNA_replicate3_run1/SGNex_HEYA8_directRNA_replicate3_run1.bam`
  (the matching fastq `…/fastq/…/…fastq.gz` is only ~few-hundred-MB, but is NOT needed — see plan).
- ⚠ **CHEMISTRY UNRESOLVED (one open item):** `samples.tsv` is collapsed/incomplete (44 dRNA
  rows, 26×SQK-RNA001 / 18×SQK-RNA002; **HEYA8 not listed, H9 only rep1_run1=RNA001**), and the
  BAM header carries no kit. So HEYA8's RNA001-vs-RNA002 is not pinned from the bucket — confirm
  from GoekeLab/sg-nex-data docs before labelling. Either way it is an **older chemistry** than
  LongBench RNA004, satisfying the two-chemistry contrast.

### LongBench (RNA004 DRS) — spike-ins CONFIRMED; ref must be local; reads are GBs
- **No reference / GTF / BAM anywhere in the bucket** (only `raw/{fastq,pod5}` + `processed/
  {cell_line_annotation,count_matrices,rds}`) → the cheap remote-BAM trick is UNAVAILABLE here;
  spike-in reads can only come from processing the full DRS fastq.
- **RNA004 DRS fastq keys + sizes** (`raw/fastq/{Sample}_dRNA_ONT.fastq.gz`, 8 cell lines;
  extra `raw/fastq/fastqs_topup/{Sample}_dRNA_ONT_topup.fastq.gz`): **H69 = 11.3 GB (SMALLEST)**,
  H1975 = 12.0 GB, HCC827 = 14.8 GB, **H146 = 15.3 GB**, H211 = 18.3 GB, SHP77 = 21.5 GB,
  H2228 = 22.8 GB, H526 = 26.8 GB (topups +15–28 GB each).
- **Spike-ins CONFIRMED present WITHOUT any read download** — grepped the DRS bulk count matrix
  `processed/count_matrices/dRNA_bulk_H69_matrix.tsv.gz` (6.3 MB, staged at
  `/scratch/users/kevinroy/sirv_work/dRNA_bulk_H69_matrix.tsv.gz`): **84 `SIRV*` transcript rows
  + Sequin `R1_*/R2_*` rows**; spike-in reads = **81,373 / 27,150,336 = 0.30 %** of H69 DRS.
  (Caveat: 0 `ERCC-00xxx` rows in the matrix — its quant annotation likely omits ERCC names; the
  reads still align to ERCC contigs in `sirv4.fasta`. Don't read "0 ERCC" as "no ERCC reads".)
- **Spike-in reference = REUSE LOCAL** `/scratch/users/kevinroy/sirv_work/sirv4.fasta` (430 KB,
  **114 contigs: ERCC-000xx + SIRV** = Lexogen SIRV-Set4) + `sirv4.gtf` (54 KB). This IS the same
  SIRV-Set4 LongBench used; it is also the SAME ref the yeast `sirv.sorted.bam` was built against
  → directly comparable. (Sequins are NOT in `sirv4.fasta`; for the SIRV↔SIRV cross-check that is
  fine — keep Sequins out, per comparability note below. Add the Garvan sequin ref from
  sequinstandards.com later only if a Sequin arm is wanted.)

### ★ COMPARABILITY (the whole point — do NOT skip)
The yeast SIRV reference (`/scratch/users/kevinroy/sirv_work/sirv.sorted.bam`, 4,424 reads,
SIRV+ERCC, length-window **[600,1000]**) and the two new measurements are only cross-checkable
if measured IDENTICALLY:
1. **Same minimap2 invocation** (`-ax splice -uf -k14 --MD`, DRS U→T fixup) and **same
   `measure-bam`/`autocorr-bam` length window [600,1000]** as `run_sirv_errstruct.sbatch`.
2. **Restrict the cross-check to SIRV+ERCC.** SG-NEx pools Sequins (`chrIS`) which the yeast/
   LongBench refs lack → **align/measure only `SIRVomeERCCome` (SG-NEx) and `sirv4.fasta`
   (LongBench); report Sequins separately, never pooled.**
3. ⚠ **The session-7 N-in-span `measure_error_structure` bug (RESUME #1) is STILL OPEN.** Fix it
   FIRST and re-measure the yeast SIRV BAM the same way, or none of these three numbers are
   magnitude-comparable. These plans produce the BAMs; the magnitude read awaits that fix.

### PLAN A — SG-NEx RNA00x (CHEAP; no GB download; ready to run tonight)
Because the genome BAM is remotely indexed, extract ONLY the spike-in reads as fastq via a remote
region query (downloads just the spike-in BGZF blocks — **verified: 12,205 reads in 3.8 s** for
rep1_run1), then realign with OUR params to the staged 11 MB ref. No fastq download.
```
ssh sherlock 'source <conda>/etc/profile.d/conda.sh; conda activate pbsim3
 W=/scratch/users/kevinroy/sirv_work; cd $W
 U=https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/sequencing_data_ont/bam/genome/SGNex_HEYA8_directRNA_replicate3_run1/SGNex_HEYA8_directRNA_replicate3_run1.bam
 # STEP 1 ALREADY DONE THIS SESSION — spike-in (SIRVomeERCCome, primary-only) fastqs are STAGED:
 #   sgnex_HEYA8_directRNA_replicate3_run1.sirvE.fq  (15,604 reads, 30M)   <- primary
 #   sgnex_H9_directRNA_replicate2_run1.sirvE.fq     ( 2,767 reads, 5.8M)
 #   sgnex_H9_directRNA_replicate4_run1.sirvE.fq     ( 3,739 reads, 8.1M)
 # (BAM SEQ is 4-bit ACGT so these are already T — no U->T fixup needed. samtools fastq restored
 #  original orientation. To re-extract / add a Sequin arm: samtools view -b -F 0x900 "$U" chrIS | samtools fastq -.)
 # STEP 2 — realign with the SAME alignment path as run_sirv_errstruct.sbatch (do NOT hand-roll;
 #  invoke that script parameterized for ref=sgnex_spikein.fa + fq, so flags -ax splice -uf -k14 --MD
 #  + processing are byte-identical = the comparability that is the whole point). Direct form:
 minimap2 -ax splice -uf -k14 --MD -t 6 sgnex_spikein.fa sgnex_HEYA8_directRNA_replicate3_run1.sirvE.fq | samtools sort -o sgnex_heya8.sorted.bam
 samtools index sgnex_heya8.sorted.bam'
 # STEP 3 — measure-bam + autocorr-bam on SIRVomeERCCome, length window [600,1000] (AFTER the N-in-span fix).
```
Download cost: a few MB of BGZF blocks (already paid — fastqs staged). **Only STEP 2+3 remain;
trivially runnable now, no approval needed.** Chemistry de-risk: H9 rep2/rep4 fastqs are staged
too, so RESUME #2 chemistry-pinning just picks among already-extracted samples (no re-run).

### PLAN B — LongBench RNA004 (⚠ NEEDS GB APPROVAL: ~11.3 GB download for H69)
No in-bucket BAM → must fetch the whole DRS fastq. **First sample = H69 (smallest, 11.3 GB),
spike-ins confirmed (0.30 % ≈ tens of thousands of SIRV+Sequin reads).** Prefer download-to-scratch
over `wget|minimap2` streaming (a transcontinental ap-southeast-2→Stanford stream that drops
mid-align wastes the compute and isn't resumable; 11 GB on scratch is fine and re-usable).
```
ssh sherlock 'source <conda>/etc/profile.d/conda.sh; conda activate pbsim3
 W=/scratch/users/kevinroy/sirv_work; cd $W
 # 1) download (FLAG: 11.3 GB) — idempotent .part->mv like run_sirv_errstruct.sbatch:
 wget -O H69_dRNA_ONT.fastq.gz.part "https://longbench-data.s3.ap-southeast-2.amazonaws.com/raw/fastq/H69_dRNA_ONT.fastq.gz" && mv H69_dRNA_ONT.fastq.gz.part H69_dRNA_ONT.fastq.gz
 # 2) align via the SAME alignment path as run_sirv_errstruct.sbatch (parameterize ref=sirv4.fasta
 #    + fq) so the U->T-fixup/--MD/flag handling is identical. NOTE (verified this session by peeking
 #    the first reads): the LongBench RNA004 fastq SEQ is ALREADY T (no U) -> the fixup is a no-op
 #    here, but route through the script anyway for byte-identical comparability. Direct form:
 minimap2 -ax splice -uf -k14 --MD -t 8 sirv4.fasta H69_dRNA_ONT.fastq.gz | samtools view -b -F4 | samtools sort -o longbench_h69_sirv.sorted.bam
 samtools index longbench_h69_sirv.sorted.bam'
 # 3) measure-bam + autocorr-bam, SIRV+ERCC, window [600,1000] (AFTER the N-in-span fix).
```
Estimated wall: download (network-bound, transcontinental — possibly 20–60+ min) + minimap2 over
27 M reads (the bulk cost; multi-thread, owners/larsms node). Submit via `sherlock-sbatch`.
**Cheaper de-risk option if 11 GB is unwanted now:** spike-in presence is ALREADY proven from the
count matrix, so a head-sample is optional; if desired, `wget` only the first ~150 MB (Range or
`--max-filesize`-style truncated gz) and align that to estimate the per-run SIRV read yield before
committing the full 11 GB. (A null head result is inconclusive at 0.30 % — don't infer absence.)

### STAGED THIS SESSION (durable on scratch `/scratch/users/kevinroy/sirv_work/`)
- `sgnex_spikein.fa` (11.1 MB) + `.fai` — SG-NEx chrIS + SIRVomeERCCome ref
- `sgnex_spikein.gtf` (198 KB) — 165 Sequin + 176 SIRV/ERCC transcripts
- `sgnex_HEYA8_directRNA_replicate3_run1.sirvE.fq` (30M, 15,604 SIRV+ERCC reads) — PLAN A first sample
- `sgnex_H9_directRNA_replicate2_run1.sirvE.fq` (5.8M, 2,767 reads) + `sgnex_H9_directRNA_replicate4_run1.sirvE.fq` (8.1M, 3,739) — chemistry-pinning alternates
- `dRNA_bulk_H69_matrix.tsv.gz` (6.3 MB) — LongBench spike-in presence proof
- `sgnex_samples.tsv` (88 lines) — SG-NEx metadata (incomplete; chemistry caveat above)
- (pre-existing, reuse) `sirv4.fasta`/`sirv4.gtf` (SIRV-Set4) and `sirv.sorted.bam` (yeast SIRV ref)

### RESUME (session-8 boundary) — nothing in flight; no marker set
1. **Decide tonight:** PLAN A is free → run it now (no approval). **PLAN B needs the go-ahead for
   the ~11.3 GB H69 download** — if approved, submit the PLAN B block via `sherlock-sbatch` (owners/
   larsms, scratch staging). 2. Pin SG-NEx HEYA8 chemistry from GoekeLab/sg-nex-data docs. 3. Both
   measurements MUST wait on / use the session-7 N-in-span `measure_error_structure` fix and the
   [600,1000] window + SIRV+ERCC-only restriction (★ COMPARABILITY above), then feed
   `error_injector.calibrate_params` per RESUME #3 of session-7.

---

## SESSION-7 (2026-06-28) — EXON-GTF loader BUILT + validated on real SIRV-Set 4 (the SIRV-fit gate dep)

The one M1-local dependency the SIRV/LRGASP absolute-truth fit was blocked on (the
exon-feature GTF loader variant of `gff_panel`, "still unbuilt" across sessions 4→6)
is **BUILT, tested, and validated on the real SIRV-Set 4 reference**. Everything else in
RESUME #1 is Sherlock-gated/heavy. Nothing in flight; smoke gate GREEN.

**DONE (benchmark-only paths; gate stays green):**
- `scripts/benchmark/sim/gff_panel.py` — added `parse_gtf_exons` (GTF: exon rows grouped by
  `transcript_id`, introns DERIVED from adjacent-exon gaps using the SAME gap test as
  `TranscriptModel.introns()`), `build_panel_from_gtf` (drop-in for `build_panel`, returns the
  identical `(models, pairs, donors, acceptors)`), a GTF attribute parser `_gtf_attrs`
  (`key "value";` — NOT GFF3 `key=value`), and `_build_models` (shared model construction the GFF
  and GTF paths both call — `build_panel` refactored onto it, behavior-preserving). `main()` gained
  `--gtf`. Transcript key = VERBATIM `transcript_id` (version included) so the LRGASP
  `read_to_isoform` join matches exactly.
- `tests/test_gff_panel_gtf.py` (6 tests, hermetic — synthetic inline exon-GTF, no network). The
  load-bearing invariant (advisor): every model junction classifies ANNOTATED against the catalogue
  from its OWN GTF (non-ANNOTATED == 0) — validates coords + intron-from-gaps simultaneously. Plus
  attr-parse, span/intron derivation, ±strand canonicity, monoexonic drop, missing-contig skip.
- **Real-SIRV validation (one-time, manual; tiny 4.6 KB GTF + 132 KB FASTA fetched to scratch — refs,
  not reads):** 176 transcripts (matches SPEC), 61 multi-exon, **288 junctions all ANNOTATED**
  (invariant holds on real data), **99.3% canonical** (GT-AG 276 / GC-AG 7 / AT-AC 3) at the derived
  boundaries → coordinates exact. (Note: `JunctionTruth.motif` STRING is a left-normalized
  genome-strand display field that reads oddly for minus-strand/slid junctions — a pre-existing
  representation quirk, NOT the loader; `.canonical` is the strand-aware truth.) **Count
  reconciliation:** `.canonical` flags 5/288 non-canonical = 2 genuinely non-canonical (CT-AC,
  CT-AG) + 3 `AT-AC` (U12 minor-spliceosome) that `_canonical_within_window` excludes by design
  (it tests GT/GC..AG only). So biological-canonical = 286/288 (99.3%); the 3 AT-AC SIRV introns
  will land in the **non-canonical FDR track** in the SIRV cycle — expected existing-scorer
  behavior (same for the GFF path), NOT a loader artifact. Don't misread it later.
- Regression: standalone tests pass (`python tests/test_gff_panel_gtf.py`); the unchanged GFF path
  re-verified through the refactored `_build_models`; smoke `--reps 20` GREEN.

**ALSO DONE this session (parallel tracks):**
- **SIRV measurement job SUBMITTED** — `scripts/benchmark/run_sirv_errstruct.sbatch` (larsms
  partition/account). Aligns real WTC11 DRS reads (ENCODE ENCSR392BGY rep `ENCFF155CFF`, ~910 MiB)
  to SIRV-Set 4 → keeps mapped spike-ins (absolute truth) → `measure-bam` + `autocorr-bam` (length
  window [600,1000]). Bakes in minimap2 `--MD` (events_from_bam_read needs MD or raises) + DRS U→T
  fixup. **Job 31822729 FAILED in 13s** (Sherlock curl 7.29.0 → error 43 on the CGL/ENCODE URLs);
  **fixed to `wget` (`.part`→`mv` idempotent)**; **re-submitted as job 31823230** (PENDING — larsms
  node busy with a549_wgs_deep). Outputs land in `/scratch/users/kevinroy/sirv_work/` (BAM,
  `sirv_errstruct.report.txt`, `.sirv_rc` sentinel rc-inside).
- **LRGASP/NanoSim sim-truth join BUILT (parallel track A)** — `scripts/benchmark/sim/lrgasp_truth.py`
  + `tests/test_lrgasp_truth.py` (6 tests, hermetic). Parses `read_to_isoform.tsv` (mouse-ENSMUST
  decoy visibility + ENST filter), builds the per-transcript truth catalogue via
  `build_panel_from_gtf` (spliced_only=False — sim reads come from all transcripts), and restricts
  each read's truth to its **spanned+anchored** junctions using the read's ALIGNED span (the
  session-3 fragment-FN bug: the scorer counts every unmatched truth junction as FN with no span
  intersection). Mirrors `pbsim3_wrapper.MIN_JUNCTION_ANCHOR=10` so the NanoSim and pbsim3 arms score
  identically. Unblocks the NanoSim arm of the pbsim3-vs-NanoSim-vs-real-SIRV three-way.
- **Plain-language overview doc** `dev/NATIVE_ALIGNER_OVERVIEW.md` (the de-novo-aligner-member
  framing: orthogonal algorithm complementing the panel for novel-isoform recovery + DRS error modes;
  C1-C6 facet table; end game). Living doc — revisit/polish.

**SIRV job — DONE (job 31823846 COMPLETED rc=0, 2m13s; 31822729 & 31823230 died on curl/gzip env traps,
both fixed → wget + `gunzip -c`).** 4424 spike-in reads mapped; report at
`/scratch/users/kevinroy/sirv_work/sirv_errstruct.report.txt`; numbers + reading in SPEC §"SIRV
ABSOLUTE-TRUTH RESULT". **Headline: the MAGNITUDE is NOT yet usable** — advisor + a check found an
**N-in-span bug** in `measure_error_structure` (`events_from_alignment` adds the intron N length to
`ref_span` and jumps `rpos` across it → spliced reads get a deflated rate + intron-crossing error gaps,
corrupting gap5x/overdisp). Defensible results: pipeline works; the window IS junction-spanning-dominated
(369 N-op vs 139 ERCC); `indel_run≥2` 0.36≈yeast 0.39 (transferable); autocorr r 0.20 positive (soft
down-weight confirmed on abs-truth). Refit deferred.

**RESUME (session-7 boundary) — SIRV measured; the magnitude needs a measurement fix first:**
1. **FIX `measure_error_structure` for spliced reads (the blocker)** — add a per-exonic-base mode: exclude
   N from `ref_span` and compute error gaps WITHIN exons (never straddling an intron). Re-measure the SIRV
   spliced reads (`/scratch/users/kevinroy/sirv_work/sirv.sorted.bam`) AND re-measure the real-yeast BAM the
   SAME way (the yeast 5.28 has the same bug, smaller because yeast introns are short) → only then are
   SIRV↔yeast comparable. (Unit-test the per-exon mode on a synthetic spliced read first.)
2. **Confirm on RNA004 SIRV** (LongBench `s3://longbench-data/`) — the ENCODE WTC11 DRS used here is
   RNA002-era (~4× the modern rate), shape-comparable only; RNA004 is the real target chemistry.
3. **THEN** re-fit `error_injector.calibrate_params` (or confirm the placeholder) + the error TYPE split
   (`empirical_cigar_error_profiler.py --isolation-flank 0` on the SIRV BAM). Do NOT drop the
   3-state/Hawkes upgrade on current evidence — the magnitude question is still open.
- **Re-run the SIRV job any time:** `ssh sherlock "cd $D && sbatch scripts/benchmark/run_sirv_errstruct.sbatch"`
  (idempotent: skips existing downloads/BAM). More reads: add reps ENCFF600LIU/ENCFF771DIX.
- **Independent of the SIRV result:** the LRGASP NanoSim arm — fetch the LRGASP sim reads +
  `read_to_isoform.tsv` + the **`with_mm_M27`** GENCODE GTF (mouse decoys!) on scratch, align, run
  `lrgasp_truth.py --bam ... --out-truth ...`, then `score_bam` for the NanoSim junction recall/FDR
  to compare against pbsim3 Tier-2 (DRS 0.816 / cDNA 0.843) and the SIRV arm.
Smoke regression gate as before: `pysam-python scripts/benchmark/smoke_roundtrip.py --out /tmp/x --reps 20` → exit 0.

---

## SESSION-6 (2026-06-27) — SHERLOCK VALIDATION RAN; verdict empirically nailed; placeholder LOCKED

Sherlock auth restored. Ran the injector validation with the NEW `measure_error_structure`
(self-consistent code across real/pbsim/injector). Full numbers in SPEC §"SHERLOCK VALIDATION
RESULT". This is a CLEAN CYCLE BOUNDARY (advisor): the deliverable is validated; the remaining
magnitude fit is the SIRV gate working as designed, not incompleteness.

**RESULTS (all gates green: smoke, injector self_check, 8 injector tests):**
- **Guard #1 PASSED** — new code reproduces the lost-err_corr verdict (disp 9.64 vs 9.95, run≥2 0.39
  vs 0.39, p90/med 2.03 vs 1.98). `sub5_gap_excess` = 5.28 (vs prior 2.83 — definitional). Real rate
  0.0153, length-invariant `overdisp_v` 0.70.
- **pbsim DECISIVELY shape-deficient** (thinned to real rate, self-consistent code): overdisp_v 0.054
  vs 0.70 (~13×), gap5x 1.16 vs 5.28, run≥2 0.25 vs 0.39. SPEC verdict empirically confirmed.
- **Real within-read autocorrelation r ≈ 0.30** (PI-#2(i) on REAL data) — moderate global hotness →
  covariate has REAL but PARTIAL signal → SOFT down-weight not hard filter (confirms PI refinement (a)
  on real reads). This r is the third constraint pinning the global/local split.
- **Model-expressiveness finding:** the 2-state HMM CANNOT jointly match real (overdisp_v 0.70, gap5x ~5,
  autocorr 0.30) — clustering↔autocorr are coupled in-model, real decouples (mostly LOCAL over-dispersion).
  Locked a HAND-PICKED `placeholder_params` matching the autocorr split (0.30); achieved overdisp_v ~0.24
  is the documented gap. SIRV fit likely needs a 3-state/self-exciting burst. Auto-fit NOT used (advisor).
- **Guard #3 (alignment-only inflation) — GAP-FREE case ~0; junction case UNRESOLVED.** Inject clean
  errors (truth gap5x 4.44) → realign → 3.99 (LOWER). Establishes only: minimap2 doesn't fabricate
  clustering from scattered clean errors in gap-free alignment. Does NOT speak to real's 5.28 (that is
  `-ax splice` read-vs-genome crossing introns; this test was `-ax map-ont` to gap-free templates — the
  junction pile-up the SPEC worried about is untested). → SIRV (junction-spanning) decides whether the
  true target is near 3 or 5; do NOT pre-bias the SIRV fit toward 5.

**CODE CHANGES this session** (committed): calibration metric → length-invariant `overdisp_v` (not
length-dependent `dispersion_index`); `REAL_TARGETS` re-measured + `head_tail_autocorr` 0.30 added;
hand-picked `placeholder_params` (PLACEHOLDER-PENDING-SIRV); `self_check` (3) decoupled from the
contaminated magnitude (asserts reachable composition bands + autocorr split); `error_realism_validate.py`
gained `autocorr-bam`, length-window filter, truth-track report on `inject-fastq`; fixed the burst
monotonicity bug (lower p_cold_to_hot → more clustering). New stat `overdisp_v`/`head_tail_autocorr`.

**RESUME (session-6 boundary) — SIRV-gated work only; nothing in flight. When resuming:**
1. **SIRV/LRGASP absolute-truth fit** (SPEC §EXTERNAL-VALIDITY has URL-verified accessions): build the
   exon-GTF `gff_panel` loader variant (still unbuilt), align SIRV reads to the SIRV reference,
   `measure-bam` + `autocorr-bam` → the CLEAN magnitude targets (junction-spanning, so it also re-checks
   guard #3's alignment-inflation on junctions). Re-fit; derive the error TYPE split (DRS deletion-dominant).
2. **3-state / self-exciting burst upgrade** — to decouple clustering (gap5x) from autocorr, which the
   2-state HMM cannot (this session's finding). Only if the SIRV targets confirm the real decoupling.
3. **PI-#2 FDR-LIFT on stratum (G)** — inject hot reads (calibrated, post-SIRV) into JUNCTION_DISCOVERY,
   measure novel-junction FDR with vs without hot-read down-weighting. Prove the lift before shipping.
   **CONTROL (advisor):** the real r≈0.30 conflates per-MOLECULE hotness with per-TRANSCRIPT alignability;
   down-weighting hot reads could suppress novel-junction support from hard-to-align transcripts (a DISCOVERY
   BIAS). Stratify by transcript / compare within- vs across-transcript autocorr to isolate per-molecule hotness
   before claiming the lift. (SPEC §SHERLOCK VALIDATION "ATTRIBUTION CAVEAT".)
4. **ORGANISM-SPECIFIC error model (design note recorded, decisions DEFERRED)** — before any human panel:
   yeast/human differ in RNA mods (the one organism-specific ERROR mechanism: yeast m6A-poor, human DRACH-pervasive),
   pA length (out-of-calibration HP regime), exon/intron architecture (panel property, not error model). Plan: ONE
   engine + per-organism PARAM SETS; "one model vs two" stays UNDECIDED. CONFOUNDS to honor: `human−SIRV` does NOT
   isolate the mod term (SIRV is IVT, own composition) — use a WITHIN-MOTIF modified-vs-unmodified contrast
   (miCLIP/GLORI); yeast-real is lower-mod, NOT mod-free; transfer requires a CONTEXT-CONDITIONED error map applied
   to human sequence (not ported yeast aggregates). Full note: SPEC §"ORGANISM-SPECIFIC ERROR MODEL".
5. C1 member code remains the NEXT GATED CYCLE.
Sherlock paths: real BAM `/scratch/users/kevinroy/rectify_wt_by4742_rep1_26167419_0/...minimap2.namesorted.bam`;
pbsim reads + templates `/scratch/users/kevinroy/err_corr_work/`; bench code
`/home/groups/larsms/users/kevinroy/aligner_bench_live/` (rsync the worktree `scripts/benchmark/` first).
Gates: `pysam-python scripts/benchmark/smoke_roundtrip.py --out /tmp/x --reps 20` + `pysam-python
scripts/benchmark/sim/error_injector.py` → both exit 0.

---

## SESSION-5 (2026-06-27) — ERROR INJECTOR BUILT (M1-local DONE; fitting/validation SHERLOCK-gated)

The FIRST TASK (the error injector the measurement proved is needed) is BUILT and
composition-verified on the M1; the Sherlock fitting/validation is turnkey. **No member
code started** (correctly — that is the next gated cycle). Regression gate GREEN.

**DONE this session (all benchmark-only paths; gate stays green):**
- `scripts/benchmark/sim/error_injector.py` — the 3-layer generative injector: (1) per-read
  multiplier (Gamma/mixture) → over-dispersion, (2) 2-state cold/hot burst HMM (marginal-preserving)
  → within-read clustering, (3) fat-tailed geometric indel-run pmf. `InjectorParams.null()` =
  backward-compat marginal Poisson. Source-DECOUPLED `inject(seq,...)` works on clean OR pbsim reads.
  `measure_error_structure` + a pysam BAM adapter (`events_from_bam_read` / pure
  `events_from_alignment`) compute the SAME 3 stats on real/pbsim/pbsim+injector.
- **TRUTH RULE enforced:** injected errors = a SEPARATE `List[ErrorEvent]` track, NEVER `IndelTruth`
  (the SESSION-2 `has_unexplained`-zeroing lesson). The track IS PI-#2's hotness label.
- **Calibration = coordinate descent on the COMBINED output** (NOT independent MoM — layers
  cross-contaminate; burst inflates dispersion so Gamma shape converges ~1.26). Robust across seeds:
  disp 9.45–10.24 / gap5x 2.71–2.95 / run≥2 0.40 (targets 9.95/2.83/0.39).
- `tests/test_error_injector.py` (8 fast tests, all PASS run directly — no local pytest env) +
  `error_injector.self_check` (9s; the local composition gate, NOT realism — labeled as such).
- **PI-#2 MECHANISM probe** `scripts/benchmark/read_reliability_probe.py`: over-dispersion drives
  exonic→junction predictability (r 0.955 vs null −0.02), burst-only stays local (r 0.011) → the
  reliability covariate has signal but is imperfect → SOFT down-weight not hard filter (PI refinement
  (a) CONFIRMED). FDR-lift NUMBER deliberately NOT computed (SIRV-gated/magnitude-sensitive).
- `scripts/benchmark/error_realism_validate.py` — turnkey Sherlock `measure-bam`/`inject-fastq`
  (full protocol in its docstring). `inject-fastq` verified on a synthetic FASTQ; **`measure-bam`
  NOT VERIFIED on a real BAM** (Sherlock down) — only the pure `events_from_alignment` is unit-tested.

**KEY DISCIPLINE (advisor-reinforced, don't regress):** the M1 self-check proves COMPOSITION, not
realism. Realism = SIRV absolute-truth magnitude + DISTRIBUTIONAL comparison, both Sherlock-gated.
Magnitude is UNKNOWN (2.83× is an upper bound — read-vs-ref conflates RNA-mod + minimap2 pile-up).

**RESUME (session-5 boundary) — Sherlock auth was DOWN (Permission denied gssapi); ask Kevin to
re-auth (`ssh sherlock true`), do NOT re-open the master. Nothing in flight. Then, in order:**
1. **Re-sync code to Sherlock** (worktree changed): `rsync -az scripts/benchmark/ sherlock:$D/scripts/benchmark/`
   where `D=/home/groups/larsms/users/kevinroy/aligner_bench_live`.
2. **Re-profile real BAM with the burst gate OFF:** `empirical_cigar_error_profiler.py --isolation-flank 0`
   (default is already 0) on `/scratch/users/kevinroy/rectify_wt_by4742_rep1_*/...minimap2.namesorted.bam`.
3. **Validate pbsim+injector vs real** (the FIRST-TASK validation): `error_realism_validate.py measure-bam`
   on {real, pbsim round-trip, pbsim+injector (via `inject-fastq` then re-align)} all with `--thin 0.019`;
   success = +injector closes the pbsim→real gap on disp/gap5x/run≥2.
   **4 GUARDS (advisor, the self-checks are circular — these make validation MEANINGFUL; full text in
   SPEC §INJECTOR BUILD "GUARDS"):** (i) RE-DERIVE targets with THIS module's `measure_error_structure`
   (the hardcoded 9.95/2.83/0.39 are from the LOST err_corr.py — don't trust); (ii) LENGTH-CONTROL the
   comparison (dispersion ≈ linear in read len → calibrate to real reads' length dist, not flat 600bp);
   (iii) verify the two measurement paths AGREE (inject→align→`measure-bam` reproduces injected-truth
   stats — `events_from_alignment` can emit ins+sub at the same ref pos → gap=0, which `inject` never
   does); (iv) derive `frac_sub/ins/del`+`base_rate` from the real profiler (DRS is deletion-dominant,
   not the PLACEHOLDER 0.55/0.15/0.30).
4. **Re-fit `calibrate_params` to SIRV ABSOLUTE-truth targets** (the real magnitude; not the read-vs-ref
   upper bound) — needs the SIRV/LRGASP reference + the exon-GTF `gff_panel` loader variant (still unbuilt).
5. **THEN PI-#2 FDR-lift on stratum (G):** inject hot reads into JUNCTION_DISCOVERY, measure novel-junction
   FDR with vs without hot-read down-weighting. Prove the lift BEFORE shipping the member feature.
6. C1 member code remains the NEXT GATED CYCLE — not a benchmark session.
Regression gate (any session): `PATH=... pysam-python scripts/benchmark/smoke_roundtrip.py --out /tmp/x
--reps 20` → exit 0; plus `pysam-python scripts/benchmark/sim/error_injector.py` → exit 0 (composition).

---

## SESSION-4b (2026-06-27) — REDIRECT: in-silico error-REALISM first (real-data = calibrator, deferred)

User redirected: prioritize **in-silico** synthetic reads (errors injected into perfect
reads = absolute truth) and get the ERROR MODEL right, BEFORE real-data. Real/spike-in
data is reframed as the **error-process CALIBRATOR** (the only contamination-free way to
measure real error correlation structure), not the primary truth set. Two user priorities:
- **(4, most important) junction-discovery bias stratum** — canonical/non-canonical ×
  annotated/unannotated, to measure which aligners over-SNAP to GT-AG/annotated vs
  reliably discover de-novo introns. PROBED + CONFIRMED: minimap2 -ax splice -uf snaps
  non-canonical junctions to nearby canonical motifs (canonical GT-AG=EXACT; CA-TC→shift
  (197,397); AT-AC→(205,400); even with -C0), and snaps even ERROR-FREE → member-addressable
  motif-bias. NOT yet built (design in this session's transcript: gen_junction_discovery_stratum
  + scorer `fp_canonical_snap` + smoke (G)).
- **error-distribution realism** — we have MARGINAL rates only; missing (2) spatial
  clustering/burst + (3) per-read over-dispersion. CONFIRMED IN CODE: the empirical profiler's
  `--isolation-flank 10` recipe (its documented usage) counts only errors with ≥10 exact
  matches both sides → by construction EXCLUDES clusters/bursts. So an injector driven by it
  under-produces the bursty/hot-read regime where the convergent-error problem lives.

**ERROR-STRUCTURE MEASUREMENT: DONE** (was in-flight; result in SPEC §ERROR-REALISM).
VERDICT: pbsim3 is MORE iid/uniform than real ONT on every axis (rate-matched to the real
1.9%). Injector needs THREE layers: (1) per-read over-dispersion (real disp-index 9.95 vs
pbsim 3.55), (2) within-read burst/clustering (real 2.83x excess sub-5bp gaps vs pbsim ~1.2x),
(3) longer multi-base indel runs (pbsim 19% vs real 39% ≥2bp). Magnitude must be calibrated
vs SIRV ABSOLUTE truth (real-yeast 2.83x is an upper bound — read-vs-ref conflates RNA-mod +
alignment artifacts a pre-alignment injector must not reproduce). Injector BUILD = next cycle.

**POINT 4 (JUNCTION_DISCOVERY) — DONE** (commits 6c40f60 + 478ad88). gen_junction_discovery_stratum
(canonical/non-canonical × annotated/unannotated) + scorer `fp_canonical_snap` + smoke (G).
VERIFIED: minimap2 recovers canonical at ceiling (1.000), SNAPS non-canonical (recall 0.000,
25 snaps), and addressability is PROVEN not asserted (all 40/40 snaps cost NM>=1, so the true
non-canonical site strictly wins on a motif-blind score → evidence-weighing member recovers it;
the smoke fails if >10% snaps are NM==0 = the zero-evidence trap). Annotation axis flat for
minimap2 (annotation-agnostic); moves uLTRA/gmap in the panel run.

**DONE this session-4b:** junction-snapping probe (above); empirical-table clean-flank bias
confirmed in code; LRGASP URL-verified acquisition plan (committed to SPEC §EXTERNAL-VALIDITY,
the pbsim3-vs-NanoSim-vs-real-SIRV three-way; RNA002 so it answers sim-realism, LongBench=RNA004
for transfer); the 4-bug red-team fixes + error-realism framework.

**RESUME (session-4b — at the HANDOFF BOUNDARY; steps 1-2 DONE this session):**
The gate phase is complete + documented. Remaining work is fresh-session/next-cycle:
1. **Error injector build** (highest priority — the error-realism verdict above). Add to the
   in-silico path: (a) per-read rate multiplier (over-dispersed mixture), (b) self-exciting /
   Markov BURST process, (c) longer multi-base indel-run model — on top of the marginal table;
   HP-only marginal stays the backward-compat default. CALIBRATE the 3 params against SIRV
   ABSOLUTE-truth error structure (read-vs-known-sequence; NOT read-vs-reference). Measure-first
   (the err_corr.py approach on Sherlock) before/after to confirm pbsim+injector matches real.
2. **EXON-feature GTF loader** (`gff_panel` variant: group exon rows by transcript_id) — the
   shared dep unlocking ALL real-data integration (SIRV + LRGASP-GENCODE are exon-GTF). Small.
3. **Real-data acquisition** (cluster, never M1): LRGASP three-way (SPEC §EXTERNAL-VALIDITY has
   curl-verified URLs, pilot ~38GB; pbsim3-vs-NanoSim-vs-real-SIRV, DISTRIBUTIONAL comparison,
   use the with_mm_M27 GTF) for sim-realism; LongBench RNA004 (s3://longbench-data/) for the
   current-chemistry transfer check. Both gated on Sherlock auth being up.
4. **C1 member code** = the NEXT GATED CYCLE (do NOT start in a benchmark session): arm-LAW vs
   arm-flat ablation via `align_exon_block_global(..., penalty_table=...)`; plug-in point is
   `run_flat_affine_arm`/`cigar_records_to_bam` in the smoke.
5. **Open red-team hardening** (from the RED-TEAM section): non-canonical/minus-strand smoke
   coverage, (E) addressability, discriminating-scale smoke, Tier-2 indel-projection coords,
   verify the reverse-strand MAF risk on Sherlock.
Smoke regression gate (any session): `PATH=... pysam-python scripts/benchmark/smoke_roundtrip.py
--out /tmp/x --reps 20` → exit 0 = (A)(A2)(B)(B2)(C)(C')(D)(E)(F)(G) green.

---

## PI DESIGN NOTES (2026-06-27) — assessed + AGREED; fold into the error-injector + C3/§8 cycles

Two PI ideas off the error-distribution measurement. Both agreed (with refinements); both are
benchmark-TESTABLE (keep the prove-don't-assert discipline).

**(1) Model the error CORRELATION structure, not just marginals (the injector model).**
AGREE — but it is NOT more columns on the marginal table (clustering + per-read hotness are
sequential/joint, not per-position). Build a lightweight generative ERROR-PROCESS model:
  - per-read rate MIXTURE → over-dispersion (the hot-read tail; real max/median ~13x vs pbsim ~4x);
  - a 2-3-state BURST HMM / self-exciting process → within-read clustering (real 2.83x excess
    sub-5bp gaps vs pbsim ~1.2x);
  - a fat-tailed INDEL-RUN-LENGTH distribution (real 39% vs pbsim 19% indels ≥2bp).
  Fit by EM/HMM/method-of-moments — **dependency-light, NOT a neural net** (project no-heavy-deps
  rule; escalate only if the simple model fails). **CRITICAL DATA CAVEAT:** the current empirical
  table CANNOT see bursts — it was built with `--isolation-flank 10`, which excludes clustered
  errors by construction. Re-profile with the gate OFF, ideally on **SIRV absolute-truth** reads
  (so real variants / alignment artifacts don't masquerade as bursts). This model feeds the
  injector (benchmark realism) AND enables (2).

**(2) Per-read "hotness" as a NOVEL-JUNCTION-DISCOVERY FDR signal (a member/§8 facet).**
STRONGLY AGREE — likely the higher-payoff idea. Mechanism: a globally-hot (error-riddled) read's
junction region is also likely error-riddled → it can map "best" to a spurious UNANNOTATED junction
by chance → counting it as support for a novel junction INFLATES discovery FDR. The enabler: a
read's hotness is ESTIMABLE FROM ITS EXONIC PORTION (error density vs reference in well-anchored
exon blocks) WITHOUT junction truth → a self-calibrating per-read RELIABILITY covariate. Works
because of the measured per-read over-dispersion. Refinements: (a) SOFT down-weight / posterior
input, NOT a hard filter — the localized-burst component means a read can be clean in the exon yet
bursty at the junction, so exonic density predicts junction reliability only probabilistically
(fits §8.2 "soft prior, never hard gate"); (b) slots into §8 discovery-FDR + the C3 calibrated
posterior as a NEW ORTHOGONAL read-reliability signal, complementary to the abstain band + anchor
gate + canonical/non-canonical tracks. **BENCHMARK TEST (build it):** once the injector makes hot
reads, measure (i) does exonic error density predict junction-region error (within-read error
autocorrelation), and (ii) does down-weighting hot reads improve novel-junction FDR on the
JUNCTION_DISCOVERY stratum (G)? The over-dispersion injector layer + stratum (G) are exactly the
pieces to test it — do NOT ship the signal as a member feature until the benchmark shows the lift.

**Thread connection:** (1) makes hot reads exist in the benchmark → which lets us validate (2).
So sequence: injector over-dispersion layer first, then the read-reliability FDR test, then (if it
shows lift) the member facet.

---

## SESSION-4 (2026-06-27) — chosen next: real-data transfer check (LongBench); BLOCKED on Sherlock re-auth

User chose to run the real-ONT external-validity transfer check (does the pbsim3
Tier-2 junction recall/FDR transfer to real reads?). Started recon; then Sherlock's
ControlMaster died (daily Kerberos/2FA lapse at the date boundary) — **needs Kevin
to re-auth** (`ssh sherlock true`); do NOT re-open the master yourself.

**Recon DONE (from M1 via public-S3 HTTP listing — no aws cli on M1):**
- LongBench bucket `s3://longbench-data/` (ap-southeast-2, --no-sign-request).
  Reads in `raw/fastq/`: per cell line (H146,H1975,H211,H2228,H526,H69,HCC827,SHP77)
  `{S}_dRNA_ONT.fastq.gz` = **RNA004 DRS** (~12–28 GB ea), `{S}_bulk_ONT.fastq.gz` =
  cDNA (~25 GB), plus PacBio/Illumina. A `fastqs_topup/` dir has more dRNA.
- ⚠ **OPEN/UNVERIFIED:** the bucket README does NOT mention SIRV/Sequin spike-ins and
  no reference/annotation/SIRV/sequin keys exist in the bucket. MUST verify spike-ins
  are actually present in these reads BEFORE downloading 16+GB (align a small sample
  to a SIRV/Sequin reference and check for mapped reads). The SIRV-Set 4 + Sequin
  reference (FASTA+GTF) must come from Lexogen (SIRVsuite / lexogen.com) + Garvan
  (sequinstandards.com), or LRGASP Synapse syn25683367/syn25683630.

**RESUME (real-data transfer check):**
1. Confirm Sherlock is back: `ssh sherlock true` (Kevin re-auths if it 2FA-prompts).
2. Build the SIRV/Sequin loader FIRST (M1-local, no cluster): extend `gff_panel.py`
   with an exon-feature GTF loader (`build_panel_from_gtf`: group exon features by
   transcript_id, sort, derive introns) — SIRV/Sequin are standard exon-GTF, unlike
   the yeast mRNA+intron GFF. Test on a synthetic exon-GTF.
3. On Sherlock (NOT M1 for the reads): fetch the SIRV-Set 4 + Sequin reference;
   `aws s3 cp --no-sign-request s3://longbench-data/raw/fastq/H146_dRNA_ONT.fastq.gz`
   (+ a cDNA `H146_bulk_ONT.fastq.gz`) to scratch; align to the SIRV/Sequin reference
   with minimap2 -ax splice (endogenous human reads won't map → mapped = spike-in
   subset); build truth via the GTF loader; `score_bam` → junction recall/FDR.
4. COMPARE to pbsim3 Tier-2 (DRS 0.816 / cDNA 0.843). Small gap = sim transfers;
   large gap = error model needs re-tuning/caveat. Use `sherlock-sbatch`; heavy →
   cluster only, never relay reads through M1.
Full dataset plan + ranked alternatives (SG-NEx, LRGASP): `SIMULATION_BENCHMARK_SPEC.md`
§"EXTERNAL-VALIDITY DATA PLAN". Nothing is mid-flight (no running job); marker cleared.

---

## SESSION-3 ADDENDUM (2026-06-26) — VARIANT/C6 stratum BUILT + verified discriminating

Continued the next-increment strata. Per advisor triage, led with **VARIANT/C6**
(the only remaining stratum that is fully wired to a scorer metric the incumbent is
below ceiling on); the other three are deferred/scaffold with reasons recorded.

- **`gen_variant_stratum` in `scripts/benchmark/sim/controlled.py`** (wired into
  `generate_corpus`, `include_variant=True`) + smoke assertion **(E)**.
- **Discriminating construct, verified EMPIRICALLY vs minimap2 -ax splice -uf:** a
  GT..AG-flanked **deletion variant ≥ ~40bp** is re-expressed by minimap2 as a
  spurious **intron (N-op)** → a FABRICATED **variant-adjacent** junction FP
  (`fp_variant_adjacent`) + the deletion mis-scored. Probed thresholds: 20/30bp
  stay `D` (correct), 40/60/100bp all become `N`. Smoke reps=20: driver
  `fp_variant_adjacent=60`, driver junction FDR=1.0.
- **SPECIFICITY proven (not just sensitivity):** SNP 3bp from donor, SNP ~100bp
  away, and a 25bp deletion are ALL handled correctly by minimap2 (0 FP). Smoke (E)
  asserts controls = 0 variant-adjacent FP, so the stratum shows the FDR is
  SPECIFIC to the splice-mimic-deletion context — guarding the future C6 member
  from a blunt abstain-near-every-variant rule. het/hom + non-Mendelian VAF (0.33)
  recorded per `VariantTruth`.
- Whole-corpus (C) indel concordance moves 0.999→0.940 — EXPECTED + correct: the 60
  splice-mimic deletions minimap2 turned into introns are now counted as genuine
  incumbent indel failures the benchmark measures (not a regression; (C) only
  requires ≥1 correct). Schema variant round-trip is exercised end-to-end (the
  `fp_variant_adjacent` count depends on variants surviving the TSV reload).
- **Anti-overcount (advisor-caught):** each VARIANT read now gets its OWN
  freshly-randomized contig carrying the same variant KIND (HP/STR already vary
  per-read; VARIANT had been replicating ONE construct `reps` times → effective
  n=#sub-cases). Now `3×reps` INDEPENDENT driver constructs; region-disjoint split
  preserved (per-contig). Verified: 60 distinct drivers all fabricate, 60 distinct
  controls give 0 variant-adjacent FP. Caveat recorded: 25bp→D / 40bp→N threshold
  is pinned to minimap2 2.28 `-k 14` — re-probe on a version/flag change.
- **PARALOG/C4 also landed this session** — see DONE/OPEN below + SPEC §Triage.
- **Smoke GREEN** with assertions (A)(A2)(E)(F)(D)(C)(C')(B)(B2) (`--reps 20`).

### SESSION-3 TIER-2 ADDENDUM (2026-06-26) — transcriptome run drafted, verified, SUBMITTED
The Tier-2 external-validity tier (Branch A) is built, verified on REAL yeast
coords, and submitted to Sherlock.
- **Aligner inventory decided the branch** (rectify env, verified): minimap2 +
  uLTRA + deSALT (binaries) + gapmm2 + mappy present; **gmap absent** (lives in
  `aligner_bench`/`apanel`/`compass` envs), **mapPacBio absent**, spoa/abpoa absent.
  None wired into the benchmark driver. So a true 5-aligner panel TAIL is NOT
  runnable today → **Branch A = minimap2 baseline** (honest labels in the JSON);
  **Branch B (panel + hard-read injection) enumerated, gated.**
- **NEW code:** `scripts/benchmark/sim/gff_panel.py` (GFF→TranscriptModel loader;
  yeast has no `exon` feature → exons = mRNA span minus `intron` features),
  `scripts/benchmark/tier2_run.py` (driver), `pbsim3_wrapper` annotation-aware
  classification + the projection fix below. `run_tier2.sbatch` (DRS+cDNA).
- **Projection bug FOUND+FIXED during real-coord verification** (the reason to
  verify before cluster spend): `project_maf_record` assigned EVERY transcript
  intron to EVERY read; pbsim reads are fragments → false FN, recall deflated
  0.79. Fixed to per-read spanned+anchored (≥10bp flank) junctions only → 0
  out-of-span junctions; real minimap2 DRS ANNOTATED recall ≈0.80–0.85 (short
  yeast introns are genuinely hard — minimap2 calls some as deletions), FDR
  ≈0.03–0.07. **This is the saturation control (harness validation), NOT
  discrimination.**
- **Job 31628436 COMPLETED (2m17s, rc=0).** DRS: 9810/10000 placed, ANNOTATED
  recall=0.816 FDR=0.051. cDNA: 9976/10000 placed, recall=0.843 FDR=0.077.
  Summaries in `$D/tier2_results/{drs,cdna}_summary.json` (with `_scope` labels);
  numbers + reading recorded in `SIMULATION_BENCHMARK_SPEC.md` §"TIER-2 BRANCH-A
  RESULT". Saturation control PASSED (harness validated on real coords at scale);
  ~0.82–0.84 recall is real minimap2/short-yeast-intron behavior, not artifact.
  (Infra: landed on AVX-512 node sh03-07n10 despite --exclude but the rectify env
  did NOT SIGILL — exclude was harmless-but-not-honored; use owners for Branch B.)
- **Remaining strata deferred with reasons** (PARALOG needs a scorer locus-
  concordance metric first; COVERAGE×Q can't discriminate until C3 consumes phred;
  PANEL_FAILURE's real validity is the Tier-2 5-aligner panel) — see
  `SIMULATION_BENCHMARK_SPEC.md` §"VARIANT/C6 stratum … Triage".

---

## SESSION-2 ADDENDUM (2026-06-26) — gate was NOT green on a clean checkout; now is

The prior VERIFIED claim "GATE smoke is GREEN" **overclaimed**: it was green only
because that session's `.fai` happened to be fresh. On a clean re-run the smoke
**crashed** (`KeyError 'hph_A_04'` in `run_flat_affine_arm`). Two fixes landed +
the cell-size MUST-HAVE was audited:

- **`021ef97` stale-.fai guard** — `pysam.FastaFile` reuses a `.fai` older than its
  FASTA (no mtime check). The benchmark regenerates `ref.fa` in place; when the
  HP_HARD (`hph_*`) contigs were added, a leftover 61-contig index from a prior run
  exposed only the OLD contig set. `score_bam` tolerated it (`genome.get(chrom,'')`
  masked the miss) but `run_flat_affine_arm` did `genome[t.chrom]` → crash, so the
  gate's two-arm DP path **never actually ran clean**. New `scorer.open_fasta()`
  drops a stale `.fai` so pysam rebuilds; routed `load_genome` + `cigar_records_to_bam`
  (header `get_tid` would also silently truncate) + the smoke's shifted-junction
  opener through it. Verified by a targeted guard test (rewrite FASTA + new contig
  with a backdated `.fai` → guard rebuilds, plain pysam does not).
- **`c680906` missing-contig now LOUD** — `genome.get(chrom,'')` silently disabled
  `normalize_junction` (junctions compared RAW = the GMAP-0.09 FP the gate prevents).
  Now counted (`reads_missing_contig`), warned once, surfaced in `summary()`.
- **Cell-size MUST-HAVE AUDITED at `--reps 400`** (the documented gate-validity bar):
  HP indel cells span **lengths 2–12 × {A,C,G,T}** (min n=183), HP_HARD {4,6,8,10,12}
  (n=400), STR ≥155 — **all clear `min_count=100`**. **Length 1 has no deletion cell
  by construction** (`_draw_k(L=1)`→0; a length-1 run can't undercall its length) but
  is fully present as **clean reads** (400/base — the FP/false-indel control). So
  `reps=400` produces a VALID C1 corpus. (Default `--reps 120` is fine for the smoke
  but the C1 ablation corpus MUST use ≥400.)
- **Smoke re-verified GREEN** after both fixes (`--reps 20`, all GATE assertions incl.
  the (D) two-arm HP_HARD-noisy DP path that was crashing).
- **`9f71770` pbsim3 LIVE round-trip — now VERIFIED on Sherlock** (was the top OPEN
  item). The fresh `pbsim3` conda env built clean (rc=0; `ERRHMM-ONT.model` DRS +
  `ERRHMM-ONT-HQ.model` cDNA present). The live run exposed **3 real mechanical bugs
  in `pbsim3_wrapper.py`** that the synthetic-MAF unit test could not — all fixed:
  - this pbsim3 build emits **gzipped single-file** outputs (`sim.maf.gz`/`sim.fq.gz`),
    not `<prefix>*.maf`/`*.fastq` → MAF discovery + FASTQ concat found nothing → 0
    truth. `parse_maf` now gz-aware; `run_pbsim3` matches `.maf(.gz)` and **RAISES**
    (not silent 0) if none; concat handles `.fq/.fastq(.gz)`.
  - templ mode emits **ONE read per template RECORD** (`--depth` does NOT multiply
    count; `--pass-num` does). The new live driver replicates each gene into N
    records → N reads/gene.
  - New `scripts/benchmark/sim/live_roundtrip.py` (self-contained ±strand spliced
    genes → pbsim → minimap2 → scorer). **VERIFIED end-to-end:** pbsim ran, gzipped
    MAF parsed, **166/180 reads placed, junction recall=0.667 FDR=0.277 on real
    pbsim reads, 0 missing contigs.**
  - **Indel exact-concordance is REPORTED but NOT gated** — it is a TIER-1 metric.
    pbsim's MAF encodes every stochastic error as an indel, so the projected per-read
    "indel truth" is thousands of scattered edits; minimap2 redistributes them and the
    scorer's per-read `has_unexplained` gate (built for Tier-1's single known indel)
    zeroes the read. This is the SPEC's two-tier split working as designed (Tier-2 =
    junction recall/FDR; Tier-1 = exact-indel), NOT a defect.

---

## RED-TEAM (2026-06-26, 4 independent auditor agents) — gate is TRUSTWORTHY, no blocker
4 read-only auditors swept truth-construction (Tier-1 + Tier-2), the scorer, and
gate-validity/claims. Verdict: fitness genuinely scores vs truth (never the internal
score), ambiguity-aware match works, the flat-affine ablation path runs. **4 verified
bugs FIXED** (commit `2f80ee5`): VARIANT eq-span too narrow (193/480, ambiguity-span
class); scorer insertion right-boundary FN; `fp_variant_adjacent` donor-only (missed
acceptor); HP `true_cigar` missing `-k` (latent).

**OPEN red-team findings (NOT yet fixed) — for the next hardening pass:**
- **Coverage: non-canonical + minus-strand FDR tracks UNEXERCISED.** Whole corpus is
  `+`-strand/canonical, so the §8 non-canonical track and the minus-strand canonicity
  fix are never hit by the smoke (classifier verified correct, but a regression would
  pass green). Add a non-canonical-NNC + a minus-strand spliced locus to the smoke.
- **(E) VARIANT proves discriminating+specific but NOT addressable** (the bar (D)/(F)
  meet). A ≥40bp-deletion-as-intron may leave no pileup signal a C6-MVP emission could
  consume — add an addressability proof or downgrade the claim.
- **reps=20 smoke is not discriminating-scale:** (D)'s C1-addressability verdict is
  decided on ~3 reads (fragile); `min_cell_reads`=9 vs the 100 floor; "cell-audited"
  has no CI guard. Make (C') a hard assert at C1 scale; run (D) at a reps with a
  non-trivial failure denominator.
- **Tier-2 indel-projection coords (REPORT-ONLY, not gated):** insertion off-by-one +
  not coalesced; minus-strand multi-base UNIQUE deletion `eq` anchored at the wrong
  end (`pbsim3_wrapper.py:142-170`). Fix for truth.tsv honesty.
- **UNVERIFIED risk (verify on Sherlock):** Tier-2 reverse-strand MAF projection
  assumes pbsim keeps the template 's' line forward; if pbsim revcomps it, `-`-read
  JUNCTIONS (gated) would be wrong. Check a real `.maf` block: template strand for a
  `-` read line.
- MINOR: `reads_missing_contig` still folds raw-coord TP/FP/FN into totals (warns but
  doesn't exclude); `locus_accuracy` diluted by unique-contig strata; STR AT-only
  cells <100 at reps=120; SPEC strata table still calls STR discriminating.

## DONE
Three ABSENT components built (everything else is reuse), all under the
benchmark-only paths the brief allows (`rectify/core/benchmark/`, `scripts/benchmark/`):

1. **Truth-propagation schema (built FIRST — it constrains the other two)** —
   `rectify/core/benchmark/truth_schema.py`. Per-read `ReadTruth` with:
   - **exact-indel truth** as a per-position `IndelTruth` carrying an
     ambiguity-equivalence span `[eq_start, eq_end)` (HP run / STR period / unique)
     — the net-indel-in-run model, so the framing metric is ambiguity-aware;
   - **junction truth** `JunctionTruth` left-normalized via
     `chimeric_consensus.normalize_junction`, with canonicity + **NIC/NNC class**
     (ANNOTATED / NIC / NNC) and the ambiguity window;
   - **CPA truth** + downstream-A context (C2); **standing-variant truth** (C6);
   - the **shared genomic-region-disjoint TRAIN/TEST split tag** (region-decided);
   - C1 HP/STR cell metadata (`base_class`, `run_copies`) for the `min_count=100` audit.
   - Lossless TSV round-trip (`write_truth_table`/`read_truth_table`). NOT overloaded
     onto the `XV` tag (a sidecar table, per the brief).

2. **Read-simulator wrapper (pbsim3)** — `scripts/benchmark/sim/`:
   - `transcript_model.py` — the truth-propagation backbone: `TranscriptModel`
     (exon structure -> spliced FASTA, `transcript_pos_to_genome` map, introns,
     NIC/NNC `junction_truths`). Junction/CPA/NIC-NNC truth is derived from THIS
     construction, NOT from the simulator.
   - `pbsim3_wrapper.py` — Tier-2 realistic reads. **pbsim3 chosen because it is
     the ONLY one of pbsim3/badread/nanosim that emits a per-read MAF
     (read<->template alignment)**, which the EXACT-INDEL framing metric requires;
     badread/nanosim give origin+identity but no per-read edit script. Parses the
     MAF and composes (read->transcript) ∘ (transcript->genome) into exact
     `ReadTruth` (genome CIGAR + indels + junctions). `--method errhmm`;
     ERRHMM-ONT model = DRS, ERRHMM-ONT-HQ = PCR-cDNA.
   - `controlled.py` — **Tier-1 controlled-error generators (the DISCRIMINATING
     tier)**: HP (A/C/G/T x 1-12, del-dominant), STR (ED-tied placements),
     JUNCTION_AMB (a GT-AG intron whose donor sits in a repeat → ambiguity test),
     CLEAN (false-indel control). Emits ref.fa + reads.fastq + truth.tsv.

3. **Ambiguity-aware per-junction + per-indel scorer** —
   `rectify/core/benchmark/scorer.py`. Reuses
   `chimeric_consensus.normalize_junction / junction_ambiguity_window /
   _canonical_within_window` so a call 1bp into a donor/acceptor repeat is NOT
   charged FP. Junction TP/FP/FN stratified by class (ANNOTATED/NIC/NNC) AND
   canonicity (separate FDR tracks, §8); indel **position-exact concordance** +
   false-indel rate stratified by C1 cell; CPA `|est-true|`; variant-adjacent FP
   tagging (C6); `score_panel` fills which-aligners-placed + panel-unplaced fraction.

4. **End-to-end smoke** — `scripts/benchmark/smoke_roundtrip.py`.

## VERIFIED — the gate is BUILT, CONTAINS a C1-addressable stratum (validated vs
## TRUTH), and the internal-DP ablation harness is RUNNABLE. (Proving the length-law
## CLOSES the gap — arm-LAW > arm-flat — is the next-cycle ablation; see OPEN.)
- **GATE smoke is GREEN** (`scripts/benchmark/smoke_roundtrip.py`, minimap2 -ax
  splice on the Tier-1 corpus):
  - (A) known junction round-trips truth->reads->minimap2->scorer as **NNC TP=30**;
  - (A2) **NIC TP=30 and ANNOTATED TP=30** — the discovery-class classifier is
    verified end-to-end (not just present); separate FDR tracks exercised;
  - (B) a deliberately 1bp-shifted junction call (into the repeat) scores
    **TP=1 FP=0** (ambiguity-aware match, the GMAP-0.09 trap defended);
  - (B2) NEGATIVE control: a junction shifted 40bp (beyond the window) is correctly **FP**;
  - (C) indel **position-exact concordance = 0.997**, false-indel-rate=0.0;
  - (D) **the gate CONTAINS a C1-addressable stratum, validated against TRUTH:**
    the live flat-affine DP (`align_exon_block_global`, the arm C1 upgrades) scores
    **0.980 on HP_HARD-noisy (BELOW ceiling vs truth)**, and **24/24 (100%) of its
    failures are systematic out-of-run deletion misplacement OR a substitution
    "repaired" with a spurious indel** — exactly the indel-vs-substitution tradeoff
    C1 targets (a run/length-aware cost keeps the deletion in-run and stops the
    mismatch-repair). The internal-DP ablation path RUNS — the DP is BAM-ized
    (`scorer.cigar_records_to_bam`) and scored on 2400 reads — so the arm-LAW vs
    arm-flat comparison is runnable; isolated-HP control = 1.000 (the SPEC
    VERTICAL-SLICE FINDING), STR = 1.000.
- **KEY method note (advisor-driven):** minimap2 and the flat-affine DP are the
  SAME error family, so they AGREE by construction (boundary_sub = 1.000 == 1.000
  after the truth-corruption fix below). C1-discrimination is therefore decided vs
  **TRUTH** (incumbent below ceiling + failures C1-addressable, both shown above),
  NOT by a flat-vs-flat separation. Whether the length-law CLOSES the measured gap
  is the next cycle's ablation, not claimable here.
- **Two scorer correctness fixes:**
  1. a per-indel TP now requires net-in-span == truth net AND no unexplained indel
     OUTSIDE every truth span (the vertical-slice ``out_run == 0`` rule); insertion
     off-by-one made consistent with the half-open span test;
  2. **indel ambiguity span now uses the base-equality SLIDE** (`truth_schema.
     deletion_ambiguity_span`, the same primitive style as the junction slide), not
     just full-period run detection — this captures HALF-PERIOD bleeds (a 'T' before
     an (AT)n run) that were charging ambiguity-EQUIVALENT STR placements as false
     errors (STR 0.91 -> 1.000). This separated the real C1 signal (HP_HARD out-of-
     run misplacement) from the scorer artifact (STR bleed).
- **Truth-corruption fix:** the boundary substitution in HP_HARD is never allowed to
  become the run base (it would silently change the effective run length).
- Schema lossless TSV round-trip verified (junction normalization, canonicity,
  HP-cell accounting, variant/split round-trip).
- pbsim3 MAF→genome **projection** verified locally on a synthetic MAF AND now on a
  **LIVE pbsim3 run on Sherlock** (session 2, `9f71770`) — see the addendum above.

## OPEN
- **pbsim3 live run: DONE (session 2, `9f71770`).** Fresh `pbsim3` conda env on
  Sherlock (rc=0). Live round-trip VERIFIED via
  `scripts/benchmark/sim/live_roundtrip.py` (166/180 reads placed, junction
  recall=0.667 FDR=0.277, 0 missing contigs); 3 wrapper I/O bugs found + fixed. See
  the addendum for the full result + the Tier-1-vs-Tier-2 indel-metric caveat.
  Code is synced to Sherlock at
  `/home/groups/larsms/users/kevinroy/aligner_bench_live/` (run with the `rectify`
  env python + `PYTHONPATH` there + pbsim3-env binaries; see RESUME).
- **Tier-1 cell sizing: RESOLVED (audited session 2).** `--reps 400` clears
  `min_count=100` for every indel cell (HP lengths 2–12, HP_HARD {4,6,8,10,12}, STR)
  and every clean cell incl. L=1 (the FP control). Use `--reps 400` for the C1 corpus.
- **NEXT-CYCLE ABLATION (gated on this gate) — does the length-law CLOSE the gap?**
  The gate now CONTAINS a validated C1-addressable stratum (HP_HARD-noisy: flat-
  affine 0.980 vs truth, 24/24 failures C1-addressable). The remaining proof is the
  C1 cycle's job: build the length-law arm and show **arm-LAW > arm-flat** on
  HP_HARD-noisy (recovering those 24/24 out-of-run misplacements). The harness is
  ready: `run_flat_affine_arm` (smoke) + `scorer.cigar_records_to_bam` are the
  plug-in point — add a `run_lengthlaw_arm` that calls
  `align_exon_block_global(..., penalty_table=...)` once C1 Phase-1 exists.
- **STR is non-discriminating after the scorer fix** (1.000 — the earlier 0.57/0.91
  was the half-period-bleed scorer artifact, now corrected). If an
  edit-distance-TIED STR discriminator is wanted, construct a case where the true
  unit-deletion placement is NOT the leftmost slide (so left-alignment genuinely
  diverges from truth) — the current whole-unit deletion is slide-equivalent and
  correctly scores 1.0.
- **Tier-2 (Branch A, yeast minimap2 baseline): DONE — job 31628436 COMPLETED rc=0.**
  Results recorded in SPEC §"TIER-2 BRANCH-A RESULT" (DRS recall 0.816/FDR 0.051;
  cDNA 0.843/0.077). Re-run any time: `ssh sherlock 'cd <D> && sbatch run_tier2.sbatch'`.
  Branch B (the real tail + cross-aligner FDR) is the gated follow-on: NOT this
  job. It needs (1) the multi-aligner panel WIRED + index-prepped (minimap2 +
  gapmm2 + uLTRA + deSALT are in the env; gmap is in `aligner_bench`, mapPacBio
  absent) so `score_panel`'s `_panel_unplaced_fraction` is "placed by NO aligner",
  and (2) an INJECTED hard sub-population (elevated error / repeat) — a clean run
  reports tail≈0, which is a false negative, not "no C5 needed". Novel-junction
  (NIC/NNC) recall also needs isoform injection (exon-skip→NIC, novel-site→NNC);
  Branch A measures ANNOTATED-recall + spurious-FDR only.
- **standing-variant/C6: DONE (session 3)** — `gen_variant_stratum`, smoke (E),
  verified discriminating + specific (see session-3 addendum).
- **paralog/C4: DONE (session 3)** — added the scorer locus-concordance readout
  (`AlignerScore.locus_accuracy`/`locus_correct`/`locus_incorrect`/`locus_mapq0`)
  FIRST, then `gen_paralog_stratum` + smoke (F). **Advisor-corrected:** the first
  cut used a window-EXCLUDING (zero-evidence) fragment — below ceiling at 0.5 but
  UNRECOVERABLE by any method incl. C4 (the vertical-slice trap; "C4 could recover"
  was false). Fixed to a WEAK-evidence fragment (covers exactly ONE distinguishing
  SNP): per-read minimap2 below ceiling (locus_acc≈0.94) AND pooling-majority at the
  lone SNP recovers the true copy 6/6 pools — the gap is provably C4-addressable
  from truth (the (D)-equivalent). Spanning reads at ceiling 1.000 (control).
  Effective diversity = n_families (3); reps = per-locus depth (a C4-pooling
  feature). See SPEC §Triage PARALOG bullet.
- Remaining strata (panel-failure/C5, coverage×Q) DEFERRED with recorded reasons
  (Tier-2-only / next-cycle-consumer) — see `SIMULATION_BENCHMARK_SPEC.md`
  §"Triage of the remaining SPEC strata". Best next: GENOMIC_A_CPA/C2 (scorer
  already emits `|est−true CPA|`; needs a walkback comparator arm to fully
  discriminate) OR the real Tier-2 transcriptome run on Sherlock.
- No C1/member code built (correctly — that is the next, gated cycle).

## RESUME
Both gating items (pbsim3 live run + Tier-1 cell sizing) are DONE. The gate is
green, valid, cell-audited, and Tier-2-mechanically verified — a clean stopping
point. The next increment is the remaining strata + the real Tier-2 run.

0. **Local smoke = regression gate first (≈30s at low reps):**
   `cd <this worktree> && PATH="/Users/kevinroy/miniconda3/bin:/opt/homebrew/bin:$PATH" \
    /Users/kevinroy/miniconda3/envs/pysam/bin/python scripts/benchmark/smoke_roundtrip.py \
    --out /tmp/bench_smoke_chk --reps 20` (exit 0 = green). Full-scale is `--reps 120`
   (slow: ~2+ min — the two-arm DP runs the Python aligner on ~2400 reads).
1. **Re-confirm the live pbsim3 round-trip (regression, any time):**
   `ssh sherlock 'D=/home/groups/larsms/users/kevinroy/aligner_bench_live; \
    CB=/home/groups/larsms/users/kevinroy/anaconda3; P3=$CB/envs/pbsim3; cd $D; \
    PYTHONPATH=$D $CB/envs/rectify/bin/python scripts/benchmark/sim/live_roundtrip.py \
    --errhmm-model $P3/data/ERRHMM-ONT.model --pbsim-bin $P3/bin/pbsim \
    --minimap2 $P3/bin/minimap2 --samtools $P3/bin/samtools --out /tmp/pbsim_live_chk \
    --seed 7 --copies 60'` → exit 0 = green (junction recall reported; indel concordance
   is Tier-1-only, see addendum). NOTE: run with the **`rectify`** env python (full deps
   incl. numpy/pandas for `import rectify`), pbsim/minimap2/samtools from the **`pbsim3`**
   env. If the worktree code changed, re-`rsync` it first:
   `rsync -az --exclude data/ --exclude __pycache__/ --exclude '*.pyc' --exclude bin/ \
    rectify/ sherlock:$D/rectify/` and `rsync -az scripts/benchmark/ sherlock:$D/scripts/benchmark/`
   (mkdir the remote `$D/scripts` parent first — rsync won't create missing parents).
2. **NEXT INCREMENT — remaining SPEC strata** (M1-light Tier-1 controlled generators;
   schema already supports all): paralog/SMN, panel-failure/C5, coverage×Q, and the
   standing-variant **C6** generator (the one with first-class variant-adjacent junction
   FDR). Build in `scripts/benchmark/sim/controlled.py` (mirror the HP/STR/JUNCTION
   generators), add a smoke assertion per stratum, size cells ≥100 at `--reps 400`.
3. **NEXT INCREMENT — real Tier-2 transcriptome run** (Sherlock, `sherlock-sbatch`
   skill, owners partition, AVX-512 exclude; do NOT relay BAMs through the M1): drive
   `live_roundtrip.py`-style propagation over a real transcript panel (yeast saturation
   control + human SMN1/SMN2 + a NIC/NNC-rich set mirroring A549 chr5) for global
   junction recall/FDR + panel-failure tail sizing. Build the genome+transcript panel
   from `TranscriptModel`s off the real GFF; the projection path is verified.
4. **C1/member code is the NEXT GATED CYCLE — not here.** The arm-LAW vs arm-flat
   plug-in point is `run_flat_affine_arm` + `scorer.cigar_records_to_bam` (smoke);
   add `run_lengthlaw_arm` calling `align_exon_block_global(..., penalty_table=...)`
   once C1 Phase-1 exists.

## FILES (all NEW, benchmark-only paths — none touch shared hot files)
- `rectify/core/benchmark/__init__.py`
- `rectify/core/benchmark/truth_schema.py`   (component 2 — schema)
- `rectify/core/benchmark/scorer.py`          (component 3 — ambiguity-aware scorer)
- `scripts/benchmark/sim/transcript_model.py` (component 1 — truth backbone)
- `scripts/benchmark/sim/pbsim3_wrapper.py`   (component 1 — Tier-2 pbsim3 wrapper; gz/fq-robust)
- `scripts/benchmark/sim/live_roundtrip.py`   (Tier-2 LIVE mechanical-integration check — session 2)
- `scripts/benchmark/sim/controlled.py`       (component 1 — Tier-1 generators)
- `scripts/benchmark/smoke_roundtrip.py`      (end-to-end GATE smoke)
- `dev/SIMULATION_BENCHMARK_SPEC.md`          (UPDATED: simulator decision recorded)

### SESSION-12 cont. — controlled re-run after Sherlock re-auth (batches of ≤3 leads, full rigor)
Infra lesson applied: ≤3 concurrent leads, each with advisor + triple adversarial panel + WRITE-TO-FILES reliability (recover from worktree if final report stalls). BATCH 1 (in flight):
- Deliverable B real-data RE-RUN (agent adf2053e3caf83776) — prior run left c1_realdata_dB/ EMPTY; cluster outputs to /scratch/users/kevinroy/c1_realdata_dB/.
- C6 variant-aware (agent ae9b0abea706d856e) — clean re-run, VARIANT stratum.
- Discovery tiebreaker (agent ac5b9ae2566a87f08) — drive real select_best_alignment on JUNCTION_DISCOVERY.
BATCH 2 (next): Flat-Q (M1, injector Q) + C1 Claim-B (M1 sim + SIRV rate read) + WS-1 (cluster, debug empty measure-bam on encff_endo_1.bam — its endo_full_measure.out is 0 bytes; encff_extract.err nonzero). On each agent report: recover from its worktree (`git worktree list`), cross-check inline, integrate (surgical add of dev/C{n}_DESIGN.md + scripts), re-run smoke, commit.

### SESSION-12 status update (live)
DONE+COMMITTED: C3 refuted (2bbe559); C4 deferred (c4e8207); C5 deferred (658b72c); **C6 deferred — read-evidence member REFUTED, VCF-residue deferred (28f6aac)**. Each recovered-from-worktree + reproduced inline (independent cross-check) + agent triple-panel converged; smoke GREEN throughout.
IN FLIGHT (3 leads, full advisor+triple-panel rigor, write-to-files reliability):
- Deliverable B real-data (agent adf2053e3caf83776, CLUSTER → /scratch/users/kevinroy/c1_realdata_dB/).
- Discovery tiebreaker (agent ac5b9ae2566a87f08, M1).
- Flat-Q quality axis (agent a789a05ca8a2aa51c, M1).
BATCH-2 QUEUE (launch as slots free, keep ≤3 leads): C1 Claim-B (M1 sim + SIRV-rate read from Sherlock) + WS-1 cleanup (CLUSTER — debug empty measure-bam: /scratch/users/kevinroy/combined_ref/ws1_cleanup/encff_endo_1.bam has endo_full_measure.out=0 bytes; check encff_extract.err).
INTEGRATION RECIPE (per agent report): `git worktree list` → its worktree; diff its controlled.py/smoke vs this worktree (expect IDENTICAL); `cp` dev/C{n}_DESIGN.md + scripts/benchmark/{probe}.py + result.txt here; re-run the probe inline as cross-check; `smoke_roundtrip.py --reps 20` (exit 0); surgical `git add` + commit `feat(aligner-bench): C{n} ...`. Worktrees may sit on stale base 255a06d (lacks scripts/benchmark/) — the agents `git reset --soft` to the branch tip; their FILES are valid regardless.
EMERGING PATTERN: every facet gate hits zero-evidence (C4,C6) or incumbent-at-ceiling (C3,C5); the genuine opportunities are REAL-DATA triggers (SMN1/SMN2 for C4, matched VCF for C6, measured panel-failure tail for C5).

### SESSION-12 CLOSE-OUT (2026-06-30, budget-paused at ~14% of 5hr limit)
COMPLETE + COMMITTED this session: C3 refuted (2bbe559); C4 deferred (c4e8207); C5 deferred (658b72c); C6 deferred (28f6aac); Discovery — REAL residual snap-FDR found + surgical scoring.py fix PROPOSED/deferred (fc3043e); Deliverable B — DEFERRED on chr5 reference-frame mismatch (64c0c28); Flat-Q — NULL, quality caveat moot (this commit). Every facet gate now has a definitive verdict. Pattern: all refute/defer via zero-evidence or incumbent-at-ceiling; the live opportunities are REAL-DATA triggers + the one Discovery scoring.py fix.
DEFERRED TO POST-RESET (3:59 PM; do NOT spawn agents/panels — INLINE only, budget-conscious):
1. **Discovery scoring.py fix** (HIGHEST value, low risk): in `_count_junction_proximity_errors` (rectify/core/consensus/scoring.py) attribute a post-`N` (intron) insertion to `intron_end` (first exon base) not `intron_end-1`. Closes 2/221 tiebreaker snap-wins; 0/570 truth penalized; no canonical_count/fence touch. GATE: full `pytest -m "not slow"` + the gmap-fence test + (ideally) a real-DRS recall spot-check BEFORE committing. Diff/rationale in `dev/DISCOVERY_TIEBREAK.md` + `scripts/benchmark/discovery_tiebreak_probe.py`.
2. **Deliverable B re-alignment** (CLUSTER): re-align `sumner_lab/.../a549_chr5_trimmed.fastq.gz` → gencode_v44 (minimap2 -ax splice + GMAP/deSALT for SQSTM1) into `/scratch/users/kevinroy/c1_realdata_dB/align/`, then run the agent's `c1_realign_3junctions.py`. Honest prior = HOLD (short reads support non-canonical at all 3 loci).
3. **C1 Claim-B** (M1 sim + SIRV-rate read): held-out injection simulator on natural long-HP templates at SIRV-measured (table-INDEPENDENT) rates; law>B0 on known edits = shape transfers. Multi-night.
4. **WS-1 cleanup** (CLUSTER, LOW value): debug empty `endo_full_measure.out` on `/scratch/users/kevinroy/combined_ref/ws1_cleanup/encff_endo_1.bam` (check encff_extract.err; likely the [600,1000] window excludes most endo reads or a per-read error).
RESUME: `git -C <worktree> log --oneline -12` shows all verdicts. For each deferred item, the named files + commands above are self-contained. Smoke GREEN at every commit.

### SESSION-12 RESUMED (budget OK: 14% USED not remaining) — Discovery fix + deferred queue in flight
⚠️ UNCOMMITTED CHANGE: the Discovery fix IS APPLIED to `rectify/core/consensus/scoring.py`
(the N-branch: `prev_rp = ref_pos + length` instead of `+ length - 1`, attributing a post-N
insertion to the first exon base). Validation so far GREEN: smoke + 49 consensus/fence/tiebreaker
tests (test_gmap_fence_regression, test_consensus_selection, test_corrected_consensus_tiebreaker).
FULL fast pytest RUNNING (base env — /Users/kevinroy/miniconda3/bin/python, has pandas+Bio+pysam;
pysam env does NOT). RESUME: check `pytest -q -m "not slow"` result; if all pass →
`git add rectify/core/consensus/scoring.py && git commit` (feat: close the Discovery snap-FDR
blind spot); if ANY fail → investigate before commit (do NOT commit a red suite). Real-DRS recall
spot-check still advisable per dev/DISCOVERY_TIEBREAK.md.
IN FLIGHT (3 agents, reliability=write-to-files): Deliverable-B re-align→gencode_v44 (adf2053e3caf83776,
RESUMED, cluster /scratch/users/kevinroy/c1_realdata_dB/align/); C1 Claim-B injection sim
(a49b24df86163d868, M1); WS-1 cleanup (a001f9fb2d9efd216, cluster combined_ref/ws1_cleanup/).
On each report: recover from worktree (`git worktree list`), cross-check, integrate, smoke, commit.

### UPDATE: Discovery fix COMMITTED (5e18960) — supersedes the "uncommitted" flag above.
Validated: smoke + 49 fence/consensus + 182 consensus/scoring/splice/junction tests (base env
/Users/kevinroy/miniconda3/bin/python). The full `-m "not slow"` suite HANGS on M1 on a heavy
validation-correction integration test (rectify.cli correct on upf1d, ~unrelated to this change) —
use the targeted module set instead, or run the full suite on the cluster. Real-DRS recall spot-check
still advisable before production human reliance. STILL IN FLIGHT: adf2053e (Deliverable-B re-align),
a49b24df (Claim-B), a001f9fb (WS-1).

### Deliverable-B: WAITING on Sherlock re-align job 32180508 (PENDING, normal partition)
The agent diagnosed the chr5 reference-frame drift, submitted job 32180508 (minimap2+deSALT de-novo
re-align of A549 chr5 reads → gencode_v44 chr5) → /scratch/users/kevinroy/c1_realdata_dB/align/, and
staged the scorer. RECOVERED to THIS worktree: scripts/benchmark/c1_realign_3junctions.py +
c1_realign_align_chr5.sbatch + dev/DELIVERABLE_B_REALDATA.md.
RESUME (concrete): `ssh sherlock "sacct -j 32180508 -X -n -o State; ls /scratch/users/kevinroy/c1_realdata_dB/align/"`
→ if COMPLETED + BAMs present: run the agent's c1_realign_3junctions.py on those BAMs at the gencode_v44
candidate coords (exact cmd in dev/DELIVERABLE_B_REALDATA.md) → SNAP/HOLD per junction (honest prior=HOLD:
short reads support non-canonical at all 3). If FAILED/CANCELLED: check the job .err, resubmit
c1_realign_align_chr5.sbatch. If still PENDING/RUNNING: keep waiting.
STILL IN FLIGHT: a49b24df (Claim-B), a001f9fb (WS-1) + this cluster job.

### IN-FLIGHT UPDATE (latest): Claim-B a49b24df DIED (stalled, zero durable output) → re-launched as a00cde628f2d37a20 (attempt #3, anti-stall: scaffold-to-files-first, minimal panels; SCAFFOLD-ONLY is an acceptable result — Claim-B is multi-night by design & C1 ships regardless on coherence). WS-1 a001f9fb still running. Deliverable-B job 32180508 still PENDING (Sherlock normal queue). If Claim-B attempt #3 also yields nothing, DEFER it (spec durable in dev/C1_DESIGN.md §"GENUINE non-circular Claim B path"); do not keep re-spawning the heaviest/lowest-value item.

### WS-1 cleanup — DONE (results on Sherlock; "bug" was OPERATIONAL not code)
The empty `endo_full_measure.out` was NOT a code bug: `measure-bam` runs clean on the completed full
BAM. Demonstrated empty-output cause = window [600,1000] on the chr1-subset endo BAM (reads=0); the
0-byte full-BAM artifact = process terminated before flush (stderr /dev/null'd, mechanism unrecoverable).
Recovered RNA002-endo error structure (disp 1.79, gap5x 1.52 at thin=0.0153/[600,1000]) shows NO
over-dispersion excess over its SIRV baseline (2.10/1.66) → the massive over-dispersion is RNA004
CHEMISTRY (hot-read tail: h69 SIRV disp 17.5/gap5x 9.34), NOT biology — confirms the RNA004 finding.
qscore fixed via qscore_from_bam.py (prior rc=141 = zcat|head SIGPIPE): RNA004 high-mean(~31)-but-bimodal
(low-q hot tail + q40 pileup), RNA002 low-mean(~17)-but-tight. Outputs: /scratch/users/kevinroy/combined_ref/ws1_cleanup/
(WS1_CLEANUP_FINDINGS.txt + report files). STILL RUNNING (confirmatory, low-pri): h69 RNA004-endo cell
jobs 32185601+32187048 → `sacct -j 32185601,32187048 -X -o State`; when COMPLETED cat
error_realism_h69endo_{measure,autocorr}.report.txt. OPTIONAL code hardening flagged (not applied): in
cmd_measure_bam move the len_min/len_max window check BEFORE events_from_bam_read (~6x less work/mem for
windowed runs); RNA004-scale measure-bam must run on a compute node ≥32G (login OOMs).

### SESSION-13 (2026-06-30 eve, after weekly-limit account switch) — Deliverable-B cluster job the ONLY open item
NOTHING LOST to the account switch (all conclusions committed). 7 of 8 gate investigations DONE+COMMITTED
(C3 refuted 2bbe559; C4 c4e8207; C5 658b72c; C6 28f6aac; Discovery FIX SHIPPED 5e18960; Flat-Q d412750;
Claim-B aff9660; WS-1 recorded 60ef525). ONLY OPEN: Deliverable-B real-data snap/hold verdict, gated on
Sherlock re-align job 32180508 (RUNNING ~1h; minimap2 done=mm2.bam, deSALT in sort/merge — the slow SQSTM1 arm).
DURABLE RESUME (no watcher needed):
  1) `ssh sherlock "sacct -j 32180508 -X -n -o State; ls -la /scratch/users/kevinroy/c1_realdata_dB/align/*.bam"`
  2) if COMPLETED + mm2.bam & deSALT bam present: run the recovered scorer
     `PATH="/Users/kevinroy/miniconda3/bin:$PATH" PYTHONPATH=. /Users/kevinroy/miniconda3/envs/pysam/bin/python \
      scripts/benchmark/c1_realign_3junctions.py` on those BAMs at the gencode_v44 candidate coords
     (exact cmd + coords in dev/DELIVERABLE_B_REALDATA.md) → SNAP/HOLD per junction (SQSTM1/TMED9/SLC35A4),
     WITH and WITHOUT the _CANONICAL_HP_PRIOR=0.5 discount. Honest prior=HOLD (short reads support non-canonical
     at all 3). Commit `feat(deliverable-B): real-data snap/hold verdict`.
  3) if FAILED/CANCELLED: check align.32180508.out + the job .err; resubmit c1_realign_align_chr5.sbatch.
  4) if still RUNNING: keep waiting (deSALT on human chr5 is slow, allow 1-2h).
A background watcher (/tmp/db_watch_32180508.sh, sentinel /tmp/db_32180508.sentinel) is armed but is CONVENIENCE
ONLY — this RESUME block is the durable signal (survives laptop sleep / session death).
FYI also-running (WS-1 confirmatory, low-pri, results land in files on scratch, NO integration needed):
h69 RNA004-endo cells jobs 32185601+32187048 → cat error_realism_h69endo_{measure,autocorr}.report.txt when done.

### WS-1 FINAL — confirmatory h69 RNA004-endo cell COMPLETED (jobs 32185601+32187048); full 2×2 measured
The 2×2 error-structure contrast is now complete (thin=0.0153, window [600,1000]):
| cell | disp | overdisp_v | gap5x | autocorr_r |
| --- | --- | --- | --- | --- |
| RNA002 endo | 1.79 | 0.065 | 1.52 | 0.43 |
| RNA002 SIRV | 2.10 | 0.094 | 1.66 | 0.48 |
| RNA004 SIRV (h69) | 17.5 | 2.03 | 9.34 | 0.57 |
| RNA004 endo (h69) | 33.9 | 2.75 | 7.43 | 0.62 |
CONCLUSION: the extreme over-dispersion / sub-5-gap-excess / high global-hotness is an RNA004 CHEMISTRY
property — present in BOTH RNA004 arms (~10-20× RNA002), driven by the hot-read tail; RNA002 endo shows NO
over-dispersion excess over its SIRV baseline → systematic-error/chemistry, NOT biology. BENCHMARK IMPLICATION
(actionable): the error injector must reproduce a per-read GLOBAL-HOTNESS tail (layer-1 multiplier, autocorr_r
~0.4-0.6), not a uniform per-base rate bump — this validates the existing error_injector.py 3-layer design
(over-dispersion + burst HMM) and says the magnitude knob targets the RNA004 hot-tail, NOT the decent bulk
(which stays RNA002-low, per the read-quality principle). qscore: RNA004 high-mean(~31)-bimodal, RNA002
low-mean(~17)-tight. Cluster CLEAR (no jobs/procs left). Results on Sherlock ws1_cleanup/ (WS1_CLEANUP_FINDINGS.txt).
Optional code hardening (flagged, not applied): cmd_measure_bam should apply the len window BEFORE
events_from_bam_read (~6× less work on windowed runs). WS-1 fully DONE.

### Deliverable-B — re-align job DONE-but-FAILED-arm; realign NOT yet valid; advisor caught a snapping-confound error
Job 32180508 FAILED at deSALT (`FAIL desalt_aln`); minimap2 arm OK (mm2.bam). Ran the scorer → n_pooled=0 at
all 3 candidate coords: the reads splice at coords +175 (TMED9) / different-intron (SLC35A4) / >400bp
(SQSTM1) off the rederive_111 candidates. ADVISOR CAUGHT MY ERROR: the "reads splice at canonical GT..AG"
motif read was measured on minimap2's SNAPPED placements — the exact confound Deliverable B exists to see
past — so it CANNOT conclude the junctions are canonical (H1 unsupported); the C1 realign on RAW SEQUENCE
(the real test) never ran (n=0). Correction committed to dev/DELIVERABLE_B_REALDATA.md. NO valid Deliverable-B
verdict yet.
⚠️ BLOCKER: Sherlock auth DROPPED after the account switch (`ssh sherlock` = Permission denied gssapi).
ACTION NEEDED FROM KEVIN: `! ssh sherlock true` to re-auth (do NOT re-open the master yourself).
RESUME (concrete, the DECISIVE path — after re-auth):
  1) Find the exact COMPASS reference: `ssh sherlock "ls /scratch/users/kevinroy/compass_a549/genome_references_latest/GRCh38*gencode*v44*.fa*"`
  2) Cheap frame check: compare its chr5 length to align/chr5.fa (181,538,259) + `samtools faidx` the TMED9
     locus from BOTH — identical seq ⇒ same build ⇒ +175 is a candidate-derivation problem (rederive coords
     off), NOT a reference shift; differing ⇒ genuine build mismatch.
  3) Re-align the reads to THAT gencode_v44 reference (minimap2 -ax splice -uf + a splice-aware aligner for
     SQSTM1's 7805bp intron; the deSALT arm FAILED — debug align.32180508.out or use GMAP) into
     /scratch/users/kevinroy/c1_realdata_dB/align_v44/ → frame matches COMPASS by construction.
  4) Run `scripts/benchmark/c1_realign_3junctions.py` (staged at /home/groups/larsms/users/kevinroy/c1_realdata_dB_code/,
     rectify importable via oak checkout — NB pre-C1 but the scorer uses _hp_edit_distance+HpPenaltyTable only,
     which works) on the v44 BAMs → per-read canon-vs-nonc likelihood = THE deliverable. Honest prior=HOLD
     (COMPASS accurate SHORT reads support non-canonical at all 3; those are untouched here).
  NB SLC35A4 candidate span 593bp vs reads 1014bp = possibly a DIFFERENT intron — verify the candidate
  identity in-frame before scoring. Result so far: /scratch/users/kevinroy/c1_realdata_dB/c1_realign_mm2only.json (n=0, documented).

### CORRECTION (ssh): Sherlock auth did NOT drop — no re-auth needed. `ssh -O check` = Master running (pid alive);
ControlPersist 168h working; both Kerberos tickets valid to Jul 1 03:23; a serial `ssh sherlock` returns
HEALED_OK with no kinit. The earlier "Permission denied (gssapi-with-mic)" was sshd MaxSessions channel
THROTTLING from too many rapid/concurrent ssh channels (watcher polling + probes + rsync retries) — the
documented self-healing transient; it cleared once the burst settled. rsync "unexpected end of file" is a
SEPARATE issue: the login banner corrupts rsync's protocol → use `ssh 'cat > dest' < file` instead.
LEVER = SERIALIZE ssh (one at a time, no fan-out), do NOT re-open the master. The Deliverable-B in-frame
realign RESUME above can proceed NOW (no re-auth required); genuine re-auth only needed if the ticket ages
out (~03:23) or the laptop sleeps/VPN drops.

### SESSION-13 STRATEGIC REDIRECT (2026-07-01) — PI reopened the program; Fable Director assessing
PI (Kevin) pushed back on "gate complete/exhausted" — correctly. Two under-tests acknowledged:
(a) NOVEL-JUNCTION DISCOVERY was tested too narrow — the C3 junction probe INJECTED a truth member and
measured arbitration, NOT the real 5-aligner panel's NATIVE recovery across a graded novelty ladder (the
"isoform flattening" headline was never sized); C5's crude uniform-error injection doesn't cover it.
(b) C2/CPA was refuted only on CLEAN genomic-A drift, never on DIRTY tails (empirical Dorado: adapter stub +
intra-tail errors, e.g. true A20 → AAAAATAAAAAGAAAAAAAAATC) — a CORRECT-step question, re-opened.
(c) C6/variant deferred too fast — injecting catalogue variants into constructed reads IS a legitimate
non-circular mechanism test.
FRAME COMMITTED: dev/RECTIFY_STRATEGIC_FRAME.md (2420a80) — whole RECTIFY package, PI's TWO ANGLES
[(1) improve trim/correct/consensus; (2) de-novo-aligner features for novel-isoform discovery: alt TSS,
unannotated+noncanonical introns, cryptic extending/truncating pA, mutant-specific], simulation-central
rationale, the plan (novel-junction ladder rungs 1-4 + error overlay; dirty-pA injector; variant injection;
cryptic-pA/alt-TSS), empirical anchors (WS-1 hot-tail, dirty-pA, variant spectrum).
IN FLIGHT: fresh FABLE Director agent ad1da5a89e4470b03 (model=fable, run_in_background) — neutral/adversarial
strategic assessment across both angles + cross-cutting "should any banked verdict re-open"; deliverable
dev/DIRECTOR_ASSESSMENT_FABLE.md; report-back-don't-commit; write-to-files reliability; ≤3 subagents no swarm.
RESUME: when ad1da5a89e4470b03 reports → recover dev/DIRECTOR_ASSESSMENT_FABLE.md from its worktree
(`git worktree list`), review its verdict + recommended first increment per angle, integrate (commit the
assessment doc), THEN — pending its go/no-go — build the FIRST increment INLINE or via ≤1 controlled agent:
likely either (Angle 1) the dirty-pA walkback stratum+injector in scripts/benchmark/sim/ + a walkback-recovery
probe, or (Angle 2) the novel-junction ladder generator (rung 1: multi-intron gene, move one SS to a nearby
unannotated canonical GT-AG) + a real-5-aligner-panel native-recovery probe (cluster run). Lock the design with
one advisor pass before building (the standing gate discipline). Env: PATH="/Users/kevinroy/miniconda3/bin:/opt/homebrew/bin:$PATH"
PYTHONPATH=. /Users/kevinroy/miniconda3/envs/pysam/bin/python; smoke gate scripts/benchmark/smoke_roundtrip.py --reps 20 (exit 0).

### Blind-spot gate DONE (6162419) + outward Fable Director launched (2026-07-01)
NOVEL-JUNCTION BLIND-SPOT gate (the Fable Director's #1, injector-free) BUILT+RUN+COMMITTED (6162419):
minimap2 flattens 47-90% of non-canonical novel junctions on ERROR-FREE reads (AT-AC/U12 0.467 blindspot,
1off 0.70, 2off 0.78, deep 0.90; canonical control 0.017) → BUILD SIGNAL, reverses "exhausted", injector-
independent + real-data-transferable. Files: scripts/benchmark/novel_junction_blindspot.py(+result), dev/NOVEL_JUNCTION_BLINDSPOT.md.
IN FLIGHT: OUTWARD Fable Director agent a40e1b0f0b3e30561 (model=fable, run_in_background) — independent
critical eval of the 5-aligner panel + survey/propose SUPERIOR orthogonal algorithms for novel-isoform
discovery (incl. cis-variant / SMN1-SMN2 paralog / trans-splicing), grounded in the measured blind-spot +
our empirical error model; deliverable dev/DIRECTOR_ALGO_EVAL_FABLE.md; reports to director; ≤3 subagents,
write-to-files, WebSearch grounded. (The INWARD Director ad1da5a89e4470b03 already reported: dev/DIRECTOR_ASSESSMENT_FABLE.md.)
RESUME on its report: recover dev/DIRECTOR_ALGO_EVAL_FABLE.md, review, commit; then the confirmed NEXT
increments are (a) full 5-aligner panel confirmation of the blind-spot on the CLUSTER (honest panel-native
recovery = ≥1 aligner), (b) exon-size/multi-intron/5'-terminal ladder rungs + isoform-level truth, (c)
addressability formalization on the hp_penalty −logP scale (the native-member design spec), (d) dirty-pA
clustered-stub gate (walkback.py:632 _max_non_stop_in_tail=1) as the parallel cheap correct-step.

### FULL-PANEL blind-spot confirmation SUBMITTED (2026-07-01) — the deciding measurement
Fable-availability VERIFIED: model="fable" FALLS BACK to claude-opus-4-8 in this env (explicit fallback
event in every transcript; the probe agent even self-reported "I am Fable 5" while actually being Opus —
self-report is unreliable, the JSONL model field is authoritative). So both Director assessments +
outward eval were fresh-context OPUS, not model-diverse. Fable needs enabling on the account.
Outward Director eval COMMITTED (dev/DIRECTOR_ALGO_EVAL_FABLE.md): panel shares a GT-AG/annotation motif
prior at SCORING (gapmm2==minimap2; uLTRA/deSALT add seeding not scoring orthogonality; gmap deep-noncanon
but ~0/111 real) → co-fails on non-canonical; most ORTHOGONAL approach = UMI/duplicate-molecule POA consensus
BEFORE placement (Sumner cDNA-UMI/SMN on-target); most PROMISING = calibrated-LLR realignment; panel-native
recovery = THE deciding number.
BUILT+SUBMITTED: full 5-aligner panel native-recovery harness (commit 9114a1c) — novel_junction_blindspot.py
--emit-corpus, panel_blindspot_score.py (per-rung panel-UNION recovery, ambiguity-aware), run_panel_on_corpus.py
(robust per-aligner via run_multi_aligner), run_panel_blindspot.sbatch. Corpus (300 reads, 5 rungs) + CURRENT
code synced to Sherlock (nj_panel_code, C1 present); job 32419195 running minimap2+gapmm2+deSALT+uLTRA (gmap
not installed on cluster). Sentinel /scratch/users/kevinroy/nj_panel/.nj_panel_rc; report
/scratch/users/kevinroy/nj_panel/corpus/panel_blindspot.report.txt; watcher armed (/tmp/nj_panel.sentinel).
RESUME: `ssh sherlock "sacct -j 32419195 -X -n -o State; cat /scratch/users/kevinroy/nj_panel/.nj_panel_rc 2>/dev/null; cat /scratch/users/kevinroy/nj_panel/corpus/panel_blindspot.report.txt 2>/dev/null"`
→ if COMPLETED rc=0: read the per-rung PANEL column. If BLINDSP stays HIGH on non-canonical across ALL 4
aligners => panel HERDS => native member JUSTIFIED (proceed to addressability formalization + member design).
If deSALT/uLTRA RECOVER what minimap2 flattens => gain is arbitration not a new placer => re-scope. If FAILED:
check nj_panel.<job>.out (likely uLTRA needs annotation/index on tiny contigs — the driver try/excepts per
aligner so others still score; re-read the report for whichever succeeded).

### CORRECTION (Fable cause) — it IS a content SAFEGUARD, not a provisioning/entitlement issue
Earlier this session I wrongly concluded (from transcript forensics) that model="fable" falls back due to
account provisioning. That is WRONG. The CLI displayed an explicit banner the director did not see/incorporate:
"Fable 5's safeguards flagged this message. The safeguards are intentionally broad right now and may flag safe
and routine coding, cybersecurity, or biology work ... Switched to Opus 4.8." So the Fable→Opus fallback is a
Fable-SPECIFIC content SAFEGUARD over-triggering on our benign yeast/human RNA-seq work (the PI's original
hypothesis was correct). My "content-free probe also fell back" and "Sonnet works so not content" arguments were
both FLAWED: the safeguard evaluates broader session context (not the literal sub-prompt), and Sonnet lacks
Fable's safeguards (so Sonnet-resolving / Fable-flagging is exactly what a Fable content safeguard looks like).
PRACTICAL: Sonnet 4.6 is VERIFIED to resolve genuinely (no fallback) and carries no such safeguard → use Sonnet
for real model-diversity cross-checks now. Fable safeguards are being refined by Anthropic; a benign-framing
prompt MAY help but is not reliable (broad classifier, evaluates context). Do NOT repeat the "provisioning"
claim.

### CROSS-MODEL synthesis (2026-07-01) — 2xOpus + 2xSonnet Directors; panel job the pivot
Verified-Sonnet Directors (genuine, unanchored from the Opus/fake-Fable ones) BOTH landed with the Opus pair:
- CONVERGENT (model-independent): native aligner CONDITIONALLY justified — blind-spot proves the MECHANISM
  (mm2 flattens 47-90% non-canonical novels, error-free, strictly recoverable) but mm2 is an UPPER BOUND;
  FULL-PANEL native recovery is THE deciding number; #1 next step = the panel blind-spot run (in flight);
  leading member = calibrated-likelihood DP on the empirical -logP scale with NO GT-AG motif prior; C1/C3/
  Discovery verdicts hold; C2 re-opens for dirty-tail (MEASURE, walkback guard indeterminate under static).
- Decision rule (both models): panel-native AT-AC recovery <~60% => member JUSTIFIED; >~80% => gain is
  arbitration/union not a new placer => pivot to correct-step + C4/C6.
- NEW integrated: (Sonnet-out) the 47% AT-AC blindspot is a FLOOR (most U12 are GT-AG/invisible); no published
  tool removes the motif prior/uses HP-law del costs (member concept is novel); mm2 has 2 snap mechanisms
  (--splice-flank motif vs --junc-bonus annotation). (Sonnet-in) add a UNION-of-aligners floor control (==
  the PANEL column in panel_blindspot_score.py); the blind-spot ("herd to WRONG junction") is a DISTINCT
  failure mode from C5 ("no placement", 0-0.4%). Docs: dev/DIRECTOR_ASSESSMENT_SONNET.md, dev/DIRECTOR_ALGO_EVAL_SONNET.md,
  updated dev/NOVEL_JUNCTION_BLINDSPOT.md.
FABLE: unavailable via a Fable-SPECIFIC content SAFEGUARD (per CLI banner) flagging biology; Sonnet 4.6 is the
verified model-diversity option (no fallback, no such safeguard). Do NOT re-run "Fable" agents expecting Fable.

### PANEL BLIND-SPOT job — RESUBMITTED as 32422876 (two prior fails fixed)
Fails fixed: 32419195 = set -u vs conda java_home.sh unbound JAVA_HOME (→ removed set -u); 32420606 =
run_multi_aligner missing output_dir+sample_name positionals (→ fixed run_panel_on_corpus.py). Panel =
minimap2+mapPacBio+deSALT+uLTRA (gapmm2 DROPPED per PI == minimap2; gmap not installed on cluster).
RESUME: `ssh sherlock "sacct -j 32422876 -X -n -o State,Elapsed; cat /scratch/users/kevinroy/nj_panel/.nj_panel_rc 2>/dev/null; cat /scratch/users/kevinroy/nj_panel/corpus/panel_blindspot.report.txt 2>/dev/null"`
→ COMPLETED rc=0: read the PANEL column per rung. AT-AC PANEL recovery <~60% => member JUSTIFIED (proceed:
addressability formalization on hp_penalty -logP + member design + exon-size/multi-intron/5'-terminal rungs +
isoform-level truth + dirty-pA correct-step gate). >~80% => pivot to correct-step + C4/C6. FAILED again: read
/scratch/users/kevinroy/nj_panel/nj_panel.<job>.out; the driver try/excepts per aligner so partial BAMs still
score — re-read the report for whichever aligners succeeded (uLTRA may need annotation/index on tiny contigs).
Code synced at /home/groups/larsms/users/kevinroy/nj_panel_code (CURRENT, C1 present); corpus at
/scratch/users/kevinroy/nj_panel/corpus. Resubmit: `sbatch scripts/benchmark/run_panel_blindspot.sbatch`.

### Independent multi-agent review of DIRECTOR_ALGO_EVAL_SYNTHESIS.md (2026-07-01)
User-supplied dev/DIRECTOR_ALGO_EVAL_SYNTHESIS.md (Opus-authored synthesis of 3 Opus algo-evals; NOT
cross-model) surfaces a NEW verified finding: a SECOND canonical preference in the CORRECTOR —
`_CANONICAL_HP_PRIOR=0.5` (junction_scoring.py:293) consumed at junction_refiner.py:647 (tier-gated OFF at
tier>=4) + an `is_novel` tiebreaker (junction_refiner.py:638,660/663) — whose discovery cost is UNMEASURED.
Also verified firsthand: cDNA-UMI penalty tables exist (penalty_scores_cdna_umi{1,2,3plus}.tsv). DISPUTED
(my quick grep found only minimap2/mapPacBio/gapmm2 via name=): the synthesis's "8 aligners wrapped incl.
winnowmap2/minisplice already present" claim.
LAUNCHED 3 SONNET reviewers (genuine cross-model; write-to-files):
- afcaa2b58ba4f3b3f = CODE-VERIFY every code claim (double-prior tier-gating + WHEN it fires; is_novel; UMI
  tables; the disputed 8-aligner-wrap; gapmm2 "85% forced-canonical"; minisplice GT/AG-only). → dev/REVIEW_SYNTHESIS_CODE_SONNET.md
- a6dc21d9456f7b705 = ADVERSARIAL steelman (is the double-prior ablation LOW-leverage since tier-gated off at
  deep deviation + lives in corrector not aligner => distinct population from the blindspot ladder?; is
  measurement->member deferring forever; reject-list too hasty; same-model non-independence discount). → dev/REVIEW_SYNTHESIS_ADVERSARIAL_SONNET.md
- a23805945182aa994 = INTEGRATION (what's already-in-plan vs genuinely NEW; updated top-3). → dev/REVIEW_SYNTHESIS_INTEGRATION_SONNET.md
KEY OPEN QUESTION for integration: is the `_CANONICAL_HP_PRIOR` ablation the SAME experiment as the panel
blindspot (aligner generation) or a DISTINCT one (corrector flattening of correctly-placed reads)? The
adversarial + integration reviewers address this. RESUME: when the 3 report, recover their
dev/REVIEW_SYNTHESIS_*_SONNET.md, reconcile (esp. whether the double-prior ablation is co-first or a distinct
corrector-step gate), update the plan/top-3, commit. Panel job 32422876 still the pivot (check
`ssh sherlock "sacct -j 32422876 -X -n -o State; cat /scratch/users/kevinroy/nj_panel/.nj_panel_rc 2>/dev/null; cat /scratch/users/kevinroy/nj_panel/corpus/panel_blindspot.report.txt 2>/dev/null"`).

### PANEL RESULT + advisor verdict (2026-07-01) — mapPacBio breaks the herd; DIRECTIONAL DECISION pending Kevin
Job 32422876 COMPLETED: panel-native recovery of non-canonical novels = ~1.000 because mapPacBio recovers what
minimap2 flattens (AT-AC/1off/2off/deep all ~1.0; VERIFIED genuine N-op at true coords). deSALT+uLTRA FAILED
(setup on synthetic per-contig corpus). Recorded + committed to dev/NOVEL_JUNCTION_BLINDSPOT.md (PANEL RESULT
section). Also RETRACTED the "47%-is-a-floor" claim (adversarial Sonnet reviewer: logical error conflating
spliceosome-type-ID with coordinate placement).
ADVISOR VERDICT: not build/pivot yet — the result RELOCATES the question from PLACEMENT to ARBITRATION+CORRECTION.
Bounded deciding set (all cheap, artifacts exist at /scratch/users/kevinroy/nj_panel/corpus/panel_bams/):
- TEST A (M1): does the shipped consensus hp_edit_distance PICK mapPacBio(true) over minimap2(snap)? (expect yes
  — snap carries flanking mismatches). Compute _cigar_hp_edit_distance per read for both BAMs; arbiter=argmin.
- TEST B: does the corrector _CANONICAL_HP_PRIOR RE-SNAP mapPacBio's recovered non-canonical junction back to
  canonical? (the double-prior ablation's natural home — on mapPacBio-recovered reads, NOT standalone; adversarial
  demoted it from co-first). ~50-line probe reusing the ladder corpus.
- ERROR OVERLAY (load-bearing, gates the pivot): re-run the ladder with RNA004-bulk error; mapPacBio's advantage
  may hold (pivot) or become fragile (calibrated-logP member revives).
Both A+B branches => "correct-step, not a new placer." Coverage rests on ONE heavy Java member (mapPacBio) => Test A matters more.
THROUGH-LINE surfaced to Kevin (his decision, not drift): C2/C3/C4/C5/C6 + (pending A+B) non-canonical-intron
placement all refuted-or-covered => evidence leans toward "panel+consensus already cover most placement; value =
(i) harden arbitration/correction to surface what mapPacBio produces, (ii) unmeasured targets: alt-TSS, cryptic
pA, SMN paralog, variant-aware." Asked Kevin: run the bounded A+B+error set now (recommended) vs discuss direction first.
STILL IN FLIGHT: code-verify Sonnet reviewer afcaa2b58ba4f3b3f (dev/REVIEW_SYNTHESIS_CODE_SONNET.md) — orthogonal to
the panel interpretation (advisor: don't block on it); resolves the disputed "8-aligner-wrap / winnowmap2+minisplice
already present" claim.
RESUME: await Kevin's direction; if "go", run TEST A (extend panel_blindspot_score.py or a small probe: per-read
_cigar_hp_edit_distance mm2 vs mapPacBio on non-canonical rungs → arbiter-picks-recoverer rate), TEST B (corrector
ablation on mapPacBio BAM), then the error-overlay ladder. Fix deSALT/uLTRA (run_multi_aligner index/annotation on
tiny contigs) + re-run panel for the full 4. Recover the code-verify review when it lands.

### Substrate critique (Kevin) + stressor runs (2026-07-01) — the prior panel result was the EASY slice
Kevin flagged: the blindspot ladder is SYNTHETIC RANDOM sequence (NOT yeast/human), SINGLE-intron, 200bp exons,
ERROR-FREE, one-contig-per-read (no decoy genome). So "mapPacBio breaks the herd" is on the easiest substrate
and hits NONE of the troublepoints. Honest state: only rungs 1-2 (motif deviation) were tested.
LAUNCHED 3 stressor panel jobs (minimap2+mapPacBio; deSALT/uLTRA still failing on tiny contigs, separate fix):
- err01 (job 32430052): RNA004-bulk 1% error overlay.  - err05 (32430054): 5% error.  - short30 (32430061): 30bp exons (short-anchor).
Corpora at /scratch/users/kevinroy/nj_panel/{err01,err05,short30}/corpus; sentinels .corpus_rc in each parent;
reports corpus/panel_blindspot.report.txt. RESUME: `ssh sherlock "for t in err01 err05 short30; do echo $t; sacct ... ; cat /scratch/users/kevinroy/nj_panel/$t/.corpus_rc; cat /scratch/users/kevinroy/nj_panel/$t/corpus/panel_blindspot.report.txt; done"`.
READ: does mapPacBio's non-canonical recovery SURVIVE error/short-anchors? If it COLLAPSES => native member
justified (matches PI's mapPacBio-real-weakness point). If it HOLDS => the pivot strengthens (but real-genome
still untested). Local minimap2-only sanity at 5% error: minimap2 still snaps (AT-AC 0.475, deep 0.05).
REMAINING BUILDS (the real troublepoints, not yet done): (1) MULTI-INTRON transcripts (compounding); (2) REAL
GENOME context — port the ladder onto bundled yeast S288C (rectify/data/S288C_*.fsa) then human, embedding loci
in real sequence so aligners face real composition + DECOY GT-AG sites (the per-contig setup removes decoy
competition where snapping/mis-mapping actually happens); (3) 5'-terminal-exon-near-TSS rung (ties Cat3). Plus the
bounded set from the advisor: TEST A (consensus hp_ed picks mapPacBio over minimap2?) + TEST B (corrector
_CANONICAL_HP_PRIOR ±50bp re-snap on mapPacBio-recovered reads). deSALT/uLTRA run_multi_aligner fix (index/annotation on tiny contigs).

### REAL-GENOME multi-intron ladder BUILT + panel runs submitted (2026-07-02)
The honest substrate (Kevin: "clearly necessary"): novel_junction_realgenome.py builds novel introns at REAL
yeast S288C loci (scanned for the target motif rung), reads spliced from real sequence + aligned to the FULL
12MB genome (real composition + DECOY GT-AG sites), MULTI-INTRON (--n-introns 2). COMMITTED. KEY CONFIRMATION
(local minimap2 on real yeast + 2-intron): minimap2 STILL flattens 40-100% of non-canonical novels (GT-AG
control 1.000; AT-AC 0.60; 1off 0.25; 2off 0.05; deep 0.00) — the snapping bias is ROBUST to real genome +
multi-intron + decoys, even slightly WORSE at deep deviation. So the blindspot is real-genome-confirmed.
IN-FLIGHT PANEL JOBS (all minimap2+mapPacBio[+deSALT/uLTRA attempt]; sentinels .corpus_rc, reports
corpus/panel_blindspot.report.txt under /scratch/users/kevinroy/nj_panel/<tag>/):
- rgen (32434900): REAL-GENOME 2-intron, error-free — THE decisive test: does mapPacBio's synthetic non-canonical
  recovery TRANSFER to real genome + multi-intron + decoys?
- rgerr (32434904): REAL-GENOME 2-intron + 2% error.
- err01/err05/short30 (32430052/54/61): synthetic error-overlay + short-anchor stressors.
ALSO IN FLIGHT: Explore agent ae22bd7c69acd3902 evaluating the Sumner-lab HUMAN rectify work for documented
per-aligner real-data strengths/weaknesses (mapPacBio + gmap documented; mm2/uLTRA/deSALT annotation-biased per
Kevin) — grounds whether mapPacBio's synthetic win transfers to human.
RESUME: `ssh sherlock 'for t in rgen rgerr err01 err05 short30; do echo == $t ==; cat /scratch/users/kevinroy/nj_panel/$t/.corpus_rc 2>/dev/null; cat /scratch/users/kevinroy/nj_panel/$t/corpus/panel_blindspot.report.txt 2>/dev/null; done'`
READ (rgen is the pivotal one): if mapPacBio recovery COLLAPSES on real genome / multi-intron / error =>
native member JUSTIFIED (mapPacBio's synthetic win was an artifact, matching PI's real-human-weakness point +
the Explore findings). If it HOLDS across all => the pivot strengthens (but human port still owed). Integrate
the Explore agent's per-aligner real-data report when it lands. REMAINING: human-genome port; Test A (consensus
picks mapPacBio?) + Test B (corrector _CANONICAL_HP_PRIOR re-snap); deSALT/uLTRA run_multi_aligner fix.

### FDR/PRECISION axis BUILT + decisive run (2026-07-02) — reconciles synthetic vs real, scores mapPacBio properly
Sumner-human eval showed mapPacBio's synthetic recall "win" IS its real-data pathology (no splice gate → emits
any gap as an intron → 97.7% spurious-novel on real human). So RECALL-ALONE is gamed. BUILT the FDR axis
(committed): panel_blindspot_score.py now reports per-aligner `.rec` (recall) + `.fpnc` (spurious NON-canonical
junctions/read); novel_junction_realgenome.py `--intronfree` control (single-exon reads, NO true junction => any
N-op is FABRICATED — reproduces the mapPacBio pathology, most discriminating under error). Local minimap2 (real
yeast + 2% err + intronfree): recall collapses on non-canonical (AT-AC 0.27, deep 0.03) BUT fabricates 0.000 on
INTRONFREE (GT-AG gate) = the precise-but-flattening workhorse profile. mapPacBio (no gate) should be the
opposite (recovers by fabricating; high .fpnc on INTRONFREE). If so => NO precise novel-non-canonical discoverer
exists in the panel => native member JUSTIFIED (on real-grounded evidence).
DECISIVE RUN: rgfdr (job 32435780) — real yeast, 2-intron, 2% error, +INTRONFREE, recall+FDR scored.
(Superseded prior recall-only runs rgen/rgerr/err01/err05/short30 for the mapPacBio question — rgfdr is the one.)
RESUME: `ssh sherlock 'cat /scratch/users/kevinroy/nj_panel/rgfdr/.corpus_rc; cat /scratch/users/kevinroy/nj_panel/rgfdr/corpus/panel_blindspot.report.txt'`
READ: mapPacBio .rec high + .fpnc HIGH (esp INTRONFREE) => fabricates, not discovers => native member justified.
mm2/uLTRA/deSALT .rec LOW on non-canonical + .fpnc ~0 => precise but flatten. Neither = precise discoverer => the
member's target. If some member has HIGH .rec + LOW .fpnc on non-canonical => panel covers it precisely (pivot).
REMAINING: human-genome port; Test A (consensus picks?) + Test B (corrector re-snap); deSALT/uLTRA run_multi_aligner fix.
