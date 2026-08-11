# HANDOFF — Phase 0 (WIP landing + master sync) → Phase 1 (Rectify Re-aligner) — CURRENT 2026-08-09

**Branch:** `chore/vendor-desalt-chanfreau1` (main checkout). Kevin approved: land all dirty WIP,
fast-forward master, then build the Re-aligner consensus-triage layer from the validated
native-aligner work. Decision record: `dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md`.
Program state of record: `dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md` (restored copy — see Files).

## Done (this session, all commits on `chore/vendor-desalt-chanfreau1`)

- **`70e9664` fix(run-all)** — the three Chanfreau-630 manifest-era bugs (sample column into
  region TSVs; `_readable_corrected_tsv()`; empty-cluster schema; `_is_manifest` superset check).
- **`5792b22` feat(qc,browser-pack)** — 630's Rbrowse unit as run-all's fail-soft final stage.
- **`3d893ef` feat(ont-cdna)** — Path A (UMI-collapse to molecules) default, both entry points;
  chunked-alignment protocol-flag refusal; UMI neighbour-sets/star-split.
- **`b9ae0e1` perf(splice)** — `_hp_edit_distance` numba fast path (memory-gated OFF) + cutoff.
- **`7648725` feat(analyze)** — 377's "group 3": containment-first gene attribution (NEW DEFAULT),
  region_class, transcript model, ncRNA atlases. Pending since 2026-07-21.
- **`999ceb5` docs** — accumulated dev record (native-aligner docs, COMPASS/GMAP, reassessment).
- Inbox: all 6 messages actioned + archived to `.claude/inbox/.read/`. Coordination note sent to
  the Chanfreau workspace inbox (630/653/641/377) with commit hashes.

## Verified

- Every commit's staged .py content parse-checked; the three mixed files (single_sample /
  multi_sample / stages) were split hunk-level across units A/B/C and reassemble **byte-identical**
  to the pre-commit working tree (zero unstaged residue after the splits).
- 630 had verified the identical tree content end-to-end (full suite `1986 passed`, run-all
  single+multi green) before handing off; my own clean-tree suite run is the in-flight gate below.

## ⚠ Incidents / concurrent-agent facts (2026-08-09 afternoon)

1. **`.claude/worktrees/` (all 11 agent worktrees) was DELETED mid-session by an external cleanup**
   — not this session. All branches/commits survive (shared object store). Known uncommitted
   casualty: the benchmark worktree's state audit, restored verbatim as
   `dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md` (provenance note at top). If other worktrees held
   uncommitted work, it is gone — flag to Kevin.
2. **A live Chanfreau session (653) is editing this checkout** (`cdna/cluster.py`,
   `cdna/consensus.py`, `tests/test_cdna_consensus_653.py`, `tests/test_umi_clustering_cdist.py`
   — left uncommitted, theirs to land). Do not touch; do not `git checkout` branches in this tree.
3. **Shell gotcha (cost an hour):** this environment's `grep`→ugrep aliasing + git pipes returned
   FALSE ZEROS (`git grep`/`git diff | wc -l` claimed the worktree branch lacked motif_blind — it
   does not). For load-bearing checks, dump blobs to files (`git show ref:path > f`) and use
   `/usr/bin/grep` / `cmp` / `diff` on the files.

## Open / in flight

- ✅ **Clean-tree full suite PASSED** (2026-08-09 evening): `1986 passed, 41 skipped, 4 deselected,
  1 xfailed` in 57:48 at `999ceb5` — identical counts to 630's baseline.
- ✅ **master FAST-FORWARDED** `255a06d → a6bdd69` (ref-update only; first sync in ~377 commits).
- Still-open bugs (filed, NOT this session's scope): AGENT_FIXES.md 2026-07-21 CRITICAL
  single-aligner path race (`align_command.py`); 2026-08-07 junction-pool-density cost scaling;
  Oak editable-install drift (Sherlock env runs uncommitted vintage).
- `feat/overhang-resolver-641`: 641/643's live branch (tip `b028e35`, 2 commits behind our line) —
  they land it themselves after cluster acceptance T3–T8; the triage layer adopts it as its
  overhang leg.

## Resume (concrete) — Phase 0 CLOSED; Phase 1 state below

**Phase 1 all landed on `feat/realigner-triage`** (worktree
`/Users/kevinroy/work/rectify_worktrees/realigner-triage`, based at `999ceb5`):

- `bd5872c` — MERGE of the native re-aligner branch (`worktree-agent-a25a2c1e784ad37dc` @
  `2b8d2ed`): motif_blind refiner, compensating-indel invariant (e40ca00), concat-DP (default ON),
  C1-del knob, dormant guard machinery, noncanon_sim harness + fence suites. Only 3 doc conflicts.
- `f2eec4a` — reimplemented `_precompute_refcol_ins_costs` (56addde landed only its call sites;
  the definition died uncommitted with the deleted worktree). Fence suites: 84 passed / 17 skipped.
  `smoke_roundtrip.py` gate: PASSED. The scoring.py post-N `prev_rp` fix arrived via the merge.
- `fc3f0c9` — the TRIAGE LAYER MVP (`rectify/core/consensus/triage.py`, `rectify triage` CLI,
  `tests/test_triage.py` = 11 passed incl. end-to-end bundled smoke).
- `eb573f9` + `34be42e` — `dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md`: the three-station
  architecture, probe-not-arbiter, scout mode, and the sequenced path.

**Next (per the path doc §4):**
1. ✅ **Full suite GREEN on `feat/realigner-triage`** (2026-08-09 evening): `2067 passed, 41
   skipped, 4 deselected, 1 xfailed, 0 failed` in 10:01 (pipefail-verified). +81 tests vs the main
   line (merged fence suites + triage). One merge collision was caught and fixed en route
   (`7c0a8f6` — the simple-slide fall-through test now pins the e40ca00 refusal outcome).
   **The branch is MERGE-READY (now at `7b931dc`).** No auto-merge — Kevin decides when it
   lands, ideally AFTER `feat/overhang-resolver-641` (merge-order note sent to 641/643).
1a. **[2026-08-10 night] 644b reviewed at v5.1 FINAL; control re-run CONFIRMS collapse; 644f
   Station-C census IN FLIGHT.** 641's v5.1 ladder (H2, @`a7bdd7e`): gold 769→852 (79% of
   mapPacBio's beyond-mm2 gold at 334 junk vs 5,456; D-merge did ~all the work — 1,421 accepts,
   `arb_shifted`=2, grammar fired 0, B1/B2/B3 accepted 0; `arb_frame_unsafe_skip`=1,446,929 =
   the stage-2 sizing). Our control re-run at `a7bdd7e`: flattening GONE (v2 ≡ mm2 every cell;
   931/1036 visits frame-refused) — addendum `ed4a805`: safety-by-IMPOTENCE; the control is the
   pre-registered STAGE-2 gate. Reply + stage-2-use-the-refiner suggestion sent to 641
   (their 3 messages archived). ✅ **644g phase-1 contest DONE + reported** (2026-08-10 ~23:59): free-N gaming MEASURED —
   junk read-slots win-or-tie the flat contest 94% vs gold 96% (no discrimination; the free
   N-op absorbs clip/mismatch cost). Gold recall through consensus CONFIRMED (35/37 residual
   gold keep winning reads). 5'-clip routing REFUTED (9% gold-read recall at 14% compute; the
   unique-gold reads are UNCLIPPED confidently-misplaced — SRC1 mechanism); distress-routing
   (644d M-D-N census) is the economical trigger. Finding `dev/PHASE1_CONTEST_644G_20260810.md`
   + tool (triage branch `4bbf870`); 641 notified. NEXT MEASUREMENTS (open): (i) the
   PRODUCTION-fenced phase-1 contest — correct+consensus over mm2+mapPacBio arms on a
   gold-window subset (do the score_segment fences refute what the flat score cannot?);
   (ii) prototype the PI's phase-2 short-exon-overhang likelihood against the 644g census
   (35 gold vs 1,737 recurrent junk separation).
✅ **644f census DONE + reported** (2026-08-10 ~21:30): true mapPacBio-only
   residual = **37 gold** (23 canonical-class incl. 20 alt introns; 14 Gould, 0/14 canonical)
   vs 5,414 junk. support≥2 + canonical-in-class → **12 gold / 13 junk (~48% precision)**;
   non-canonical track ~0.3% at ANY support threshold because **33% of mapPacBio junk RECURS**
   (systematic, not per-read random — scout-mode statistic corrected in the landscape doc §3b).
   Station-C spec = two-track admission. Finding `dev/STATIONC_MAPPACBIO_HARVEST_20260810.md` +
   tool `scripts/benchmark/644f_stationc_gate_census.py` (both on the triage branch, `687d1b1`);
   641 notified (incl. the 22→37 residual correction for their 644b doc). Next measured probe:
   mm2-side distress co-localization at the recurrent non-canonical junk vs Gould loci
   (644d machinery).
1b. **[2026-08-09 late] Resolver check + the missing discovery control RUN.** 641 is ACTIVE
   (v2 junction re-arbitration, commits through `de98896` 20:06; planning/644b/644c — Kevin's
   junction-proximal-ED hunch adjudicated SRC1: two REAL alt donors, 67/121 read split, minimap2
   forces one). Our control (`dev/REARB_NONCANON_CONTROL_20260809.md` + committed driver):
   **v2 rearb flattens cryptic non-canonical junctions −6.0 pp (HP) / −2.5 pp (plain)** — 15
   grammar-tiebreak shifts (all cryptic, some at EQUAL ED) + 18 margin-beats (pool-level-only
   fixable). Three asks sent to 641's inbox (split arb_grammar knob; distinct tag on grammar
   moves so triage never bypasses them; adopt the control into their verdict criteria).
   RESUME: check 641's inbox reply + their H2 job `14271360` (their sentinel `.p644b_v2_rc`);
   the control re-runs via
   `python scratchpad/rearb_noncanon_control.py` on any merged scratch of the two branches.
2. Triage-policy tuning on real corpora (upf1Δ 617 gold windows = ready truth set); wire the
   triage clip legs to Cat3/rescue-gate + resolver machinery.
3. Station C (pool-level discovery gate) + the second-corpus 8×-recall reproduction (Sherlock).
4. 641/643 own landing `feat/overhang-resolver-641` (rebase onto master `a6bdd69`, T6 decoy
   mitigation before human use); the `[[617]]` leave-one-out decides mapPacBio's panel status.

## Files

- Decision + design: `dev/DENOVO_ALIGNER_REASSESSMENT_20260809.md`
- Program state of record (restored): `dev/ALIGNER_BENCH_STATE_AUDIT_20260721.md`
- Native-aligner sources: branch `worktree-agent-a25a2c1e784ad37dc` @ `2b8d2ed`; key commits
  `69a230f` (motif_blind), `e40ca00` (compensating-indel invariant), `dd257b8`/`e1ed90c`
  (concat-DP), guard-shelve verdict `d5b25d3`, Phase-6 overview `fc58950`.
- Prior CMA-session handoff: superseded; see `dev/CMA_PROGRESS.md` + git history of this file.
