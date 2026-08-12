# HANDOFF — 679 cDNA trim bug (F1 root cause) · DRS handed to unit 682 · 313 GB move ~1/3 done — CURRENT 2026-08-11 ~23:00 PDT

## Delta ~23:55 — candidate-blowup ceiling IMPLEMENTED (Kevin: "implement it now")

- **`ResolverConfig.max_candidates_per_clip` (default 2000)** + **`ResolverStats.
  refused_candidate_blowup`** (flows into `stats.json` via `as_dict()` automatically).
  In-loop breaker in BOTH enumeration branches: past the budget the clip is abandoned, the
  read passes through untouched, the counter increments.
- **Chose the in-loop breaker over a pre-count** because the DP's pruning cutoff tightens as
  `best_ed` improves, so the loops must keep their order; a pre-count would have to duplicate
  the range arithmetic and could drift from it. Bail cost is bounded by the ceiling itself.
- **LOUD by construction** (the guard must not become the false-green class it prevents):
  a per-contig-throttled WARNING naming the contig + suggesting RECTIFY_SKIP_REGIONS, the
  count in the INFO summary, and a separate escalated WARNING so an acceptance gate skimming
  the summary cannot miss a non-zero field buried mid-row.
- **Default 2000 is PROVISIONAL and deliberately loose** — ~300x the post-trim mean
  (6.7 cand/clip) and >5x the worst legitimate figure ever seen (361/clip, 681 CP6).
  A test pins `> 361` so nobody tightens it into a silent-loss setting. Tighten from real
  percentiles when 682's DRS stats.json lands.
- **Measured on a GT/AG-rich tandem-repeat contig: 680 candidates unbounded vs 6 at
  ceiling=5** — the stall converted to a bound.
- 🔴 **PROCESS NOTE: my first version of the test passed VACUOUSLY.** The synthetic genome had
  no splice sites, so `candidates_evaluated=0` and the breaker never executed — a test
  advertising coverage it did not have, i.e. the exact class of this whole sweep. Caught by
  printing the stats instead of trusting the green. The fixture now asserts
  `candidates_evaluated > 100` as an anti-vacuity guard, and the module docstring warns
  anyone editing the fixture. **Tests: 7/7.**

## Delta ~23:35 — 🟢 DRS SIZING ANSWERED (~15x): unit 682's sampler says ~1.4 h/sample, not ~22 h

- **682 instrumented a throughput sampler** at
  `/u/project/guillom/kevinroy/682_sentinels/wt_rep1_progress.tsv`
  (epoch/elapsed_s/resolver_bam_bytes/minimap2_bam_bytes/frac). At 1,202 s it is at
  **frac=0.23211**, i.e. **~5,180 s ~= 1 h 26 m projected per sample.**
- **Versus the pre-numba 668b run: 0.1875 of the arm in 14,917 s => ~22.1 h/sample.**
  ⇒ **~15.4x from `RECTIFY_HP_ED_NUMBA=1` (277c708) on DRS.** The 9-sample cohort goes from
  untenable (22 h/sample against a 24 h h_rt — it would have failed on walltime regardless of
  what killed it) to roughly 3 waves at `-tc 4`.
  ⚠️ Caveats: BAM-size is a rough progress proxy (compression varies) and the rate may drift on
  later contigs; treat as +-. **This 15.4x is the NUMBA kernel on DRS — do not confuse it with
  the 15.5x the TRIM fix gave on cDNA (681 CP6). Different mechanisms, independent, and they
  happen to land on nearly the same factor.**
- Move `14293343`: **4/12 files freed, scratch 1,977.59 -> 1,845.47 GB (132 GB reclaimed)**,
  files 5-8 in flight. Working exactly as intended (per-file verify-then-rm).

## Delta ~23:20 — false-green sweep committed; 681 was self-committed and VERIFIED

- **681 needed no commit from me — the agent did it: `b3a8c35` (23:04:27),
  6 files / 574 insertions, and it ran the full suite: `2228 passed, 0 failed`
  (baseline 2216 at e7499ef).** Better than their own estimate of 2225. That agent is
  STILL ACTIVE (`rectify/core/splice/splice_aware_5prime.py` +72 lines at 23:07) working the
  1-bp G-overhang / 5'-rescue-arbiter thread — **leave their working tree alone.**
- **MY COMMIT `4533de5` — the false-green class:**
  - `--require-aligners` on `rectify align` + an always-on greppable `DROPPED-ALIGNER`
    summary naming requested-vs-produced. Gate sits AFTER the aligner loop and BEFORE the
    resolver post-pass. Default stays permissive. `tests/test_align_require_aligners.py`
    (5 tests) pins the contract AND asserts the real call site still carries the gate and
    still precedes the resolver.
  - `run_overhang_resolver(threads=)` now WARNS instead of silently ignoring — the body is
    single-threaded; `align_command` forwards `args.threads` in good faith and the 668b DRS
    job requested 8 slots for a 1-core stage.
  - Verified: 9 passed (align CLI set) + 119 passed/1 skipped (overhang/resolver/junction-scoring).
- **644k sentinels — BOTH shapes found and fixed on H2** (not repo code):
  - `644k_blast_radius.sh` ended with a hardcoded `echo 0 > .644k_rc` while **5 of 8 legs
    failed** (drs100k cliplegs SIGTERM 143; cdna_wt1 all four: 1/139-SIGSEGV/1/1). Sentinel
    corrected 0 → 5, original preserved at `.644k_rc.orig_false_green`; script now derives
    `$FAILS` via a `_rec` helper.
  - `644k_round3.sh` had the SAME hardcoded zero but **never reached it** — it died at
    18:58:54 mid-`cliplegs`, so `.644k_round3_rc` is ABSENT. Now derives from STATUS3.
  - ⚠️ **Correction to my earlier "cdna_wt1 total loss":** true only of the blast-radius run.
    Round 3 re-ran cdna_wt1 `classify` (rc=0, 580 s) and `realign` (rc=0, 7834 s) — those ARE
    usable. Full write-up on H2: `644k_blast_radius/README_644k_SENTINEL_CORRECTED.md`.
  - ⚠️ While patching I broke `status()`'s own `>> $WORK/STATUS` redirect with a too-broad
    string replace, caught it on the diff, restored it. `bash -n` would NOT have caught it —
    diff every generated patch.
- **STILL OPEN (not done):** the resolver **candidate-blowup ceiling** + a
  `refused_candidate_blowup` counter, so an unanticipated repetitive contig degrades to
  passthrough instead of a silent 555x stall. Today's mitigation is opt-in
  `RECTIFY_SKIP_REGIONS`, which requires the operator to know in advance. The 681 trim fix
  cut candidates 6,137 -> 60 on normal data but does not bound the pathological case.
- Also still open: **674 needs `h_data` raised** before it is unfrozen (OOM'd at 22.469 G
  against a 24 G ceiling).

## Delta ~23:15 — 🔴 CORRECTION: 681 WAS claimed (Kevin was right) — and the cDNA fix is DONE

- **I said "`HANDOFF_679_cdna_trim.md` is UNCLAIMED". WRONG.** NN **681** was claimed by session
  `cdna-trim-fix` at 04:52:09Z, ~1 h before I wrote that. I inferred it from an empty `qstat` —
  but a code fix does not appear in the queue. **The NN registry is the authoritative test:**
  `~/work/bin/nn_list.sh ~/work/UCLA/Chanfreau_Lab/planning`.
- **THE cDNA BLOCKER IS FIXED.** Unit 681, all four phases landed, acceptance measured:
  - mean max soft-clip **65.8 → 3.8 nt** on the bulk population (148.0 → 7.6 nt on the
    previously-untrimmed stratum)
  - `XQ==0` **52.5 % → 13.5 %**; the genuine-defect stratum **39.7 % → 0.65 %**
  - resolver **15.5× faster from the trim alone; candidate DP work −99.0 %** (6,137 → 60);
    ~361 candidates per clip quantified
  - **`resolved = 0` on BOTH arms** — independently corroborates 679 CP7 on a different slice:
    waste bug, not integrity bug, and the fix changes no call
  - over-trim guard clean: mapped 99.9 % both arms, mean aligned span **−1 nt** while ~58 nt of
    clip disappears
  - records: `HANDOFF_681_cdna_trim_fix.md` + Chanfreau `planning/681_*` (CP0–CP8)
- **🔴 THREE OF MY CLAIMS WERE WRONG — corrected in place; do not act on the originals:**
  1. **"pileup branch broken at 94.8 %" — RETRACTED (679 CP9(b), 681 CP2).** `pileup_consensus`
     uses `get_aligned_pairs(matches_only=True)` and excludes soft-clips by construction, so it
     has no adapter to strip; `XQ==0` there is CORRECT (13.7 nt clips vs 145.2 nt for real
     failures). **My acceptance criterion 2 would have BROKEN a correct branch.** I read a rate
     as a defect without crossing it against clip length.
  2. **My mechanism candidates (a)/(b)/(c) were incomplete** — (a) and (c) refuted. Real cause:
     **`orient` is computed in BAM-SEQ frame while the emitted consensus is flipped to basecalled
     frame, and nothing flips the label alongside the sequence.**
  3. **Quoted suite baseline 2,259 was stale** — real branch-point baseline is **2,216**.
- **⚠️ THE FIX IS UNCOMMITTED WORKING-TREE STATE** on `chore/vendor-desalt-chanfreau1` @ e7499ef:
  `M rectify/core/cdna/{_constants,consensus,io}.py`, `M
  rectify/core/commands/cdna_correct_command.py`, `M tests/test_cdna_correct.py`, `??
  tests/test_cdna_pretrim_frame.py`. **Do not `git add -A`** (Kevin has unrelated WIP). The cDNA
  freeze lifts only from a COMMITTED, suite-green tree.
- **New finding from 681 worth Kevin's eye:** the fixed 53-nt strip **deletes 1–3 bases of real
  mRNA at the TSS on 10.47 %** of molecules, and a 1-bp G overhang reaches the 5′ rescue arbiter
  (`DEFAULT_JUNCTION_PROXIMITY_BP = 10`, no minimum clip length) on ~34 % — supporting Kevin's
  "1 bp G overhangs deserve extra suspicion". Owned by 681, not me.

## Delta ~23:00 — SESSION WRAP. What is still in flight, and what is NOT mine.

**In flight from THIS session — exactly one thing:**
- **H2 job `14293343` move676par** — 313 GB scratch→project move. **Working as designed:
  4 of 12 files DONE+FREED, scratch 1,977.59 → 1,882.82 GB (94.8 GB reclaimed), 8 files left,
  174 G landed on project.** Sentinel `/u/project/guillom/kevinroy/.676b_move_ubam_rc`.
  Per-file verify-then-rm, so relief keeps arriving without me. Full branch logic in §Resume.

**In flight but NOT MINE — do not duplicate or interfere:**
- **H2 job `14294705` `p682_drs_res`** (qw, 8 slots, cwd `/u/scratch/k/kevinroy/668_drs`,
  script `682_drs_resolver.sh`). **A dedicated DRS agent has picked up
  `HANDOFF_668_drs_arm.md`, claimed NN 682, and queued the resolver rerun** — it also created
  `/u/scratch/k/kevinroy/682_numba_cache`, so it is running with the numba kernel as briefed.
  The DRS arm is THEIRS now. Do not resubmit 668b, do not touch `668_drs/` outputs.
- The cDNA trim fix belongs to whoever takes `HANDOFF_679_cdna_trim.md` — not yet claimed as of
  this writing.

**Confirmed quiet:** SCG has nothing running (my `scancel 52353209_17` took effect).
Sherlock has only two unrelated `BeginTime`-pending housekeeping jobs (`compass_`, `oak_usag`).
No cDNA arm is queued anywhere — the freeze is actually in force, not just declared.

**Nothing else from this session is running.** 679 census (14293031) and 679b frame diagnostic
(14293585) both finished rc=0 and their results are folded into `planning/679` (CP8/CP10).



## Delta ~22:35 — TWO dedicated-agent briefs written; 668b failure re-diagnosed

- **`HANDOFF_668_drs_arm.md`** (this repo) — brief for a dedicated DRS agent. Owns the 9-sample
  668 cohort (wt/rna15/ysh1 x3 = Kevin's P1 multi-strain set) end to end.
- **`HANDOFF_679_cdna_trim.md`** (this repo) — brief for a dedicated cDNA-trim agent.
- **668b exit-1 RE-DIAGNOSED from the script + log (better than the earlier "undiagnosed"):**
  the stdout dies inside the align call *before* `ALIGN_SECONDS` prints, AND the script's
  `trap _finish EXIT` never fired (no `.rc668b`, no `=== EXIT rc= ===` line), AND there is no
  `FATAL: align+resolver` line. ⇒ **the shell was killed hard, not exited** — not walltime
  (h_rt 24h), not memory (828M/24G), not the script's own error branch. Node-level event or the
  align child taking the shell with it. Still needs confirming, but the search space is now small.
- **🔴 SIZING: the DRS cohort was planned against a wrong figure.** At job end the resolver BAM
  was 462,880,768 B vs a 2,468,849,530 B mm2 arm = **~19% of the arm in 4h08m** (BAM-size proxy,
  +-) ⇒ **~22 h/sample single-threaded against a 24 h h_rt** — it would have failed on walltime
  anyway. The "DRS = 1,600 r/s class" assumption is ~2 orders off on this substrate.
- **The lever is landed and verified:** `origin/master` = **`277c708`** (numba twin, 289x,
  bit-identical, suite **2262 passed/41 skipped/1 xfailed, SUITE_RC=0**); **numba 0.66.0 +
  edlib both importable in H2 `~/.conda/envs/rectify`**. The smoke did NOT have it —
  `668b_resolver_smoke.sh` pins `670_tree` @ **`e7499ef`** (pre-11ed647 gate, pre-numba), so a
  master tree is both faster AND a different config; re-validate, don't assume carry-over.
  `tree_11ed647` exists but is one commit behind master.
- ⚠️ **Rerun trap flagged in the brief:** delete the TRUNCATED
  `668_drs/wt_rep1/align/wt_rep1.overhang_resolver.bam` (462 MB, no BGZF EOF) + any
  `resolver_delta.bam` by **explicit named path** before rerunning, or `--trust-existing-bams`
  may consume a truncated file. Aligner BAMs + `multialigned.bam` are COMPLETE — keep them; the
  rerun is resolver+consensus+prescan only.
- **662-frozen side effect:** the live-tree swap is unblocked (the script's "662 still imports
  from it" constraint has lapsed).

## Delta ~22:15 — fresh-agent brief written; dataset arms triaged (DRS GO / cDNA FROZEN)

- **`HANDOFF_679_cdna_trim.md` (this repo) = the standalone brief for a dedicated successor**
  on the cDNA trim. Self-contained: proven facts, the two DISCARDED mechanisms, the unpinned
  line + how to pin it, 4 phases, 6 acceptance criteria, instruments, traps, durability rule.
- **Dataset triage: Chanfreau `planning/680_dataset_status_drs_go_cdna_hold.md`.** Kevin:
  DRS proceeds unblocked, cDNA waits. Correct because 679 lives in `rectify/core/cdna/` — DRS
  never enters that path (no SSP/UMI/consensus).
  - 🔴 **668b (DRS gate) FAILED: exit 1 @ 4h08m — NOT walltime (h_rt 24h), NOT memory
    (maxvmem 828M). Cause UNDIAGNOSED.** All 3 aligner BAMs + multialigned.bam COMPLETE;
    only the resolver post-pass is partial (462 MB). Rerun is cheap via `--trust-existing-bams`.
    Resolver managed <=260 r/s on 3.88M reads — far under the 1,600 r/s "DRS class" sizing;
    turn on `RECTIFY_HP_ED_NUMBA=1` before re-timing.
  - 674 cold_shock (cDNA): exit **137 OOM** — maxvmem 22.469G vs `h_data=3G x 8` = 24G ceiling.
  - 662 (cDNA): 10/15 done, 5 unrun. 678 (cDNA): qdel'd smoke. Both frozen.
  - 673 (cDNA, reporter): Sherlock 38646035_11 COMPLETED; **SCG 52353209_17 CANCELLED BY ME**
    at 6h17m — it was producing consensus molecules the fix invalidates. Confirmed affected
    (planning/673 runs `correct-cdna` -> stage1_consensus.fastq.gz).
  - ⚠️ **644k `.644k_rc`=0 is a FALSE GREEN** — the log carries `Segmentation fault` on
    `rectify triage` (cdna_wt1 realign leg) yet the sentinel says 0. Do not quote 644 results
    until rc propagation is fixed.
- **All cDNA output already produced through `correct-cdna` is suspect** (incl. the 10 finished
  662 samples + the Sherlock 673 smoke): ~half of records carry untrimmed adapter. Not
  fabricated (679 CP7 = 0 junctions placed off adapter clips) but 5' ends/clips are wrong.
  **Do NOT flip alignment_inventory rows to YES and do NOT build browser packs on them.**

## Done

- **679** — answered Kevin's question: the "genuinely high-complexity" cDNA clips were
  untrimmed adapter + 27-nt UMI + poly(A), not genomic sequence. Correctness cleared (0
  junctions were placed off them). Record: Chanfreau `planning/679` (CP0–CP10).
- **681** — the fix, by a dedicated agent: `b3a8c35`, suite **2228 passed / 0 failed**
  (baseline 2216). Soft-clip 65.8→3.8 nt; `XQ==0` 52.5→13.5 % (real defect 39.7→**0.65 %**);
  resolver **15.5×** faster from the trim alone, candidate DP **−99.0 %**.
- **4533de5** — `--require-aligners` + always-on `DROPPED-ALIGNER` summary; `threads` no-op on
  `run_overhang_resolver` now warns instead of silently lying.
- **24c7805** — pathological-contig breaker: `max_candidates_per_clip` (2000) +
  `refused_candidate_blowup`, loud by construction. 680 candidates unbounded vs 6 at ceiling=5.
- **644k sentinels** — both shapes fixed on H2; `.644k_rc` corrected 0→5.
- **676 storage — COMPLETE.** 12/12 uBAMs moved scratch→project, "ALL 12 MOVED + VERIFIED",
  **scratch 1,977.59 → 1,657.13 GB (320 GB reclaimed; 343 GB headroom)**, 0 left on scratch.
  The orphaned rsync temp from the qdel'd 14293025 (`.prp28_rep1.ubam.bam.lAyJ9o`, 11.98 GB)
  removed by explicit named path; project dir now exactly 313 G / 12 files.
- **680** — dataset triage. **The cDNA freeze has LIFTED** (its blocker, 681, is committed+green).

## Verified

| claim | evidence |
|---|---|
| 681 fix real + green | `b3a8c35`; suite 2228/0; focused sets re-run by me 37 passed/2 skipped |
| my two commits safe | align CLI 9 passed; resolver/junction set **126 passed / 1 skipped** (was 119 pre-change) |
| **DRS ~15× faster with numba** | 682 smoke finished **frac 0.99966 at 4,927 s ≈ 1 h 22 m/sample**, vs pre-numba ~22.1 h. Projection made at frac 0.232 (~1 h 26 m) held to within ~5 % |
| ceiling guard fires and bounds | 680 candidates unbounded vs 6 at ceiling=5; 7/7 tests, incl. an anti-vacuity assert |
| 644k was a false green | STATUS: 8 legs, 5 non-zero (SIGSEGV 139, SIGTERM 143, three rc=1) under a hardcoded `echo 0` |
| `pt` constraint satisfied | all 12 retired uBAMs have dorado uBAM deposits carrying `pt:i:`; 9 differ by exactly 27 bytes = one `@PG` line |
| move completed cleanly | sentinel `.676b_move_ubam_rc`=0; log "ALL 12 MOVED + VERIFIED"; 12 files, 0 dotfile temps |

## Open

1. **676 (a) deletion — NOT started, needs Kevin.** Data is safe on project; nothing deleted.
   Run the count proof (§Resume 1) and **show Kevin the per-file list before any delete.**
2. **Candidate ceiling constant** — 2000 is provisional. Tighten from real percentiles once a
   682 `stats.json` exists (§Resume 3). Must stay **> 361** (test-pinned).
3. **674 needs `h_data` raised** before it reruns (OOM'd at 22.469 G against a 24 G ceiling).
4. **668b exit-1 root cause** still unconfirmed ("shell killed hard"); unit 682 owns it — its
   off-scratch sentinels will show whether it recurs.
5. **This branch is BEHIND `origin/master`.** `chore/vendor-desalt-chanfreau1` = e7499ef + three
   new commits; `origin/master` = `277c708` (numba twin). All three must be reconciled onto
   master — **nothing here is on master yet.**
6. **Other sessions own live work — do not touch their trees or jobs:** unit **682** (DRS,
   array `14295895` running, `HANDOFF_682_drs_arm.md`), unit **684** (cDNA P1 wave, job
   `14295902`, `HANDOFF_684_p1cdna_1M.md`), and the 681 agent (still editing
   `splice_aware_5prime.py` on the G-overhang thread).

## Resume

```bash
# 1. 676 (a) — the deletion gate. Data is on project; scratch copies are GONE.
ssh hoffman2 'ls -l /u/project/guillom/kevinroy/676_ubam_retire/'
```
Count proof — **md5 is the WRONG test** (9 files differ from their deposit twin by exactly 27
bytes of `@PG`; the 3 wtaa are genuine shard merges):
- 9 single-file cases: `samtools view -c` retired vs
  `shared/raw/cDNA/intronic_pa_snp1_prp28_prp5_cdna_pcb114_2026/<sample>.bam`
- wtaa_rep1/2/3: summed `samtools view -c` over
  `shared/raw/cDNA/cdna_wt_rna15_ysh1_20260711/wtaa_rep*/PBM45482_*.bam` vs the merged count
- deposit-side integrity from the deposit's own `MANIFEST.sha256`
**Then show Kevin the per-file list. Delete only on his go-ahead, explicit named paths, no
globs.** `pt` is already confirmed preserved in every deposit (planning/676 §9).

```bash
# 2. Watch the two live waves (NOT yours — recover only if their sessions died)
ssh hoffman2 'qstat -u kevinroy; ls -l /u/project/guillom/kevinroy/682_sentinels/'
```
- **682 (DRS)**: `.rc682_<sample>`=0 + `.done682_<sample>` → read the acceptance block
  (records==primary==distinct QNAMEs; resolver BAM + `stats.json`; XB delta) before trusting a
  wave. Job gone with NO sentinel = the 668b failure recurring → `qacct -j <id>` for the node,
  and delete any truncated `*.overhang_resolver.bam` by explicit path before rerunning.
- **684 (cDNA P1)**: must be running a tree at **`b3a8c35` or later** — anything earlier
  reintroduces the untrimmed-adapter defect. Verify before trusting its output.
- Re-runs of frozen cDNA arms (662's 5 remaining + the 10 already-done, 673, 674, 678) are
  resolver+consensus+prescan only via `--trust-existing-bams`; add `--require-aligners`
  (4533de5) while touching those scripts.

```bash
# 3. Tighten the candidate ceiling (guard is LANDED; only the constant is open)
ssh hoffman2 'cat /u/scratch/k/kevinroy/668_drs/*/align/*overhang_resolver.stats.json'
```
Set `max_candidates_per_clip` from `candidates_evaluated`/`clips_assessed` percentiles across
DRS (that JSON), cDNA post-trim (~6.7/clip; pre-trim 361/clip), and the reporter contigs (673).
**Keep it > 361** — test-pinned, because too low silently refuses legitimate clips. A non-zero
`refused_candidate_blowup` on healthy data means the constant is too tight, not that the data
is bad.

```bash
# 4. Reconcile onto master (nothing from tonight is on master)
cd ~/work/rectify && git log --oneline origin/master -1   # 277c708
```
Merge/rebase `chore/vendor-desalt-chanfreau1` (b3a8c35, 4533de5, 24c7805) onto `277c708`, run
`pytest -m "not slow"` on the result, and push. **Surgical staging only** — Kevin has unrelated
WIP (`Rectify_readme_changes_KR_proposed.txt`, `TODO_KR_NOTES.txt`,
`scripts/validation_data/*`), and the 681 agent has live edits in `splice_aware_5prime.py`.

## Files

- Briefs: `HANDOFF_668_drs_arm.md` (→ picked up as 682) · `HANDOFF_679_cdna_trim.md`
  (superseded banner listing its 3 errors) · `HANDOFF_681_cdna_trim_fix.md` ·
  `HANDOFF_682_drs_arm.md` · `HANDOFF_684_p1cdna_1M.md` (other sessions' records)
- Chanfreau `planning/`: `679_*` (evidence) · `680_*` (triage) · `681_*` (the fix, CP0–CP8)
- Commits on `chore/vendor-desalt-chanfreau1`: `b3a8c35` · `4533de5` · `24c7805`
- H2: `/u/project/guillom/kevinroy/{676_ubam_retire/,676b_move_ubam_par.sh,679_*,682_sentinels/,logs676b/}`
  · `/u/scratch/k/kevinroy/644k_blast_radius/{README_644k_SENTINEL_CORRECTED.md,*.pre644kfix}`
- Sentinels: `.676b_move_ubam_rc`=**0 (done)** · `.679_rc`=0 · `.679b_rc`=0 ·
  `644k_blast_radius/.644k_rc`=5 (corrected) · `682_sentinels/.rc682_*` (pending)

# HANDOFF — cDNA resolver throughput ROOT-CAUSED + numba kernel landing; 673 taken over — CURRENT 2026-08-11 ~19:50 PDT (NB: earlier "08-12 00:xx" stamps in this file were clock-skewed; H2 log times are authoritative)

## Delta ~19:50 — the cDNA-class slowdown: root-caused, kernel built, landing

- **673 TAKEOVER EXECUTED** (staging agent halted cleanly): SCG + Sherlock ControlMasters
  re-established per policy (one attempt each, Kevin approved both Duo pushes). Smokes
  sacct'd: SCG 52353209.17 RUNNING 3h16m, Sherlock 38646035.11 RUNNING 2h50m — both
  crawling at the PRE-gate e7499ef trees. Artifacts: Chanfreau planning/673* @ 2b82ebe.
- **The 1–4 r/s cDNA pathology is NOT the pA/period hypothesis** — measured at the
  period-gated tree: 2k reads didn't finish in 9.5 min; arb-off still 1 r/s; cProfile (50
  production reads): **122 s of 125 s inside pure-Python hp_edit_distance_bounded — 35,275
  calls = ~570 candidates/clip × 3.5 ms.** Long high-complexity cDNA clips legitimately
  pass the gate; W_max clamps at max_intron → hundreds of GT/AG candidates × 200-len DP.
- **FIX: numba twin kernel, commit `277c708`** (branch perf/hp-ed-bounded-numba):
  bit-identical (4,000 randomized pairs exact; every-row pruning preserved), **289× on
  resolver-shaped inputs** (→ ~1 r/s becomes ~200+ r/s). Same RECTIFY_HP_ED_NUMBA knob/
  guard as splice_aware (default OFF, spawn-RSS rationale; job scripts opt in). Targeted
  74/74; **full gate RUNNING** (log scratchpad/suite_numba_kernel.log, pid .gate_pid).
- numba installs: H2 user env IN PROGRESS (/tmp/numba_install.log, "Executing
  transaction"); SCG env TBD (envs at /oak/.../Users/kevinroy/envs — check rectify_scg);
  Sherlock TBD.
- **RESUME on gate green:** ff master + push 277c708 → update H2 tree_11ed647 (git fetch
  --depth 1 + reset per git-1.8 recipe) → add `export RECTIFY_HP_ED_NUMBA=1` to
  678/668b wave scripts → re-measure the 2k slice (target ≥100 r/s) → resubmit 678 smoke
  (14289428 was qdel'd — futile at 1 r/s) → THEN the 673 decision: kill both crawling
  smokes, update SCG/Sherlock trees + numba envs, rerun smokes at the fixed tree
  (planning/673 has refs/commands), arrays only on both-green. 668b (DRS) unaffected by
  the cDNA pathology (DRS = 1,600 r/s class) — its acceptance read still pending.

# (prev) beta staging waves launched — (clock-skewed stamp)

## Delta ~00:50 — Kevin's dataset priorities mapped + first wave submitted

- **Kevin's order**: P1 = WT/Ysh1-AA/Rna15-AA in DRS+cDNA (rbrowse multi-strain +
  multi-modality test); P2 = Rrp6-AA/Sen1-AA/Sen1×Rrp6/Dis3×Rrp6/Ski7/Upf1/Ski2/4NQO both
  modalities where available. **Full plan + inventory mapping: Chanfreau
  `planning/678_p1cdna_resolver_prescan.md`** (nn 678 claimed atomically).
- **P1 maps PERFECTLY to staged data**: 668 DRS manifest = exactly wt/rna15/ysh1 ×3;
  587 cDNA deposit holds WT_BY4742/wtaa/rna15aa/ysh1aa ×3 (P1) + ski7d/upf1d ×3 (P2 head
  start), all stage-1+panel-complete.
- **LAUNCHED: 678 cDNA wave smoke = H2 `14289428.1`** (WT_BY4742_rep1; resolver+gated
  prescan @ fresh clone `tree_11ed647`; script `587_ms2_cdna/678_p1cdna_resolver_prescan.sh`,
  sentinels `.rc678/.done678`). On acceptance READ → `qsub -t 2-12 -tc 4`, du-gated.
- **DRS P1 gate = 668b acceptance** (14287090.1 still RUNNING ~01:00; started 15:07 PT —
  nearing the 6–9 h window). On rc=0 + acceptance → `qsub -t 2-9 -tc 4 668b_resolver_smoke.sh`.
- Storage truth: freed space is filesystem-wide; MY user quota = 1.8/2.0 TB — still the
  gate. Sherlock scratch = authorized overflow; SCG available (Kevin's reminder).
- NEXT UNITS: (a) correct→consensus→analyze→browser-pack chain script (stage 2; gated pool;
  the expensive stage) — write fresh, smoke-first; (b) P2 deposit inventory pass
  (Rrp6-AA/Sen1-AA/doubles/Ski2 locations + modality availability, Hoffman2_datasets tab);
  (c) 644k round-3 still running (cdna_wt1 realign); (d) clip-leg phases + 669 items queued.
- RESUME: poll `.rc678` + `.rc668b` + `.644k_round3_rc` (+ 674/662/673 per task #13); read
  acceptance blocks BEFORE fanning; planning/678 carries the wave commands.

# (prev) 11ed647 LANDED; processing control active — ~23:55 PDT

## Delta ~23:55 — the landing is complete; overnight state

- **origin/master = `11ed647`** (suite 2,259/0): 649 complexity gate DEFAULT-ON
  (--complexity-alpha 0.01) + 305daff threaded determinism + §4c gate reconciliation +
  --trust-existing-bams + 669 cost docs. Lines notified (inbox 0650Z). **670_tree pin
  LIFTS after 668b/674 complete** — converge cluster trees on 11ed647 at the next safe
  window. rpl20b_myc_reporter launchable UNMASKED on trees ≥ b1b033f.
- Sentinels at 23:50: 668b RUNNING (acceptance read due — expect completion overnight;
  read per HANDOFF_cdna_umi_and_runall.md before any fan-out), 674 RUNNING, 662 tasks
  10–15, 644k round-3 in cdna_wt1 realign (long — big consensus BAM), 673 self-running.
- **RESUME (next session):** (1) poll all sentinels (commands in task #13 + earlier
  deltas); 668b done → READ acceptance block, record wall time → DRS fan-out decision
  (threading/sharding vs serial). (2) 644k round-3 done → compile cDNA blast radius +
  clip-leg timing → report Kevin + inbox. (3) Start clip-leg Phase 1/2 (task #11) on a
  branch from 11ed647. (4) 669 Items 2/3 (task #12). (5) Beta packs sequencing with
  630-analysis once first datasets clear correct.

# (prev) PROCESSING CONTROL HANDED; 649 landing in flight — ~22:50 PDT

## Delta ~22:50 — Chanfreau dataset-processing control is now MINE (Kevin via umi-perf-accuracy)

- **Handover ACK'd** (their inbox note 2026-08-12T0400Z + live SendMessage). Authoritative
  state: `~/work/UCLA/Chanfreau_Lab/HANDOFF_cdna_umi_and_runall.md` — READ IT before any
  acceptance action. I own: 668b resolver smoke (H2 14287090.1, sentinel
  `668_drs/wt_rep1/.rc668b`, acceptance read due on completion, wall time = DRS fan-out
  sizing datum), 674 cold_shock (14287595.8, `.rc674`, then ≤4-waves -tc 2 du-gated), 662
  tasks 10–15 (`.done662`; Gsheet flips ONLY after reading acceptance), 673 SCG/Sherlock
  cohort (self-running by their subagent; takeover artifacts = Chanfreau planning/673* if
  their session dies). Storage HARD STOP 1.85 TB (measured 1.6 TB at 22:45). Five open
  decisions on task #13; first = DRS fan-out threading vs serial (after 668b's wall).
- Their note pre-dates my pushes: master ALREADY has RN fix/resolver CLI/tag-test fix.
  Genuinely new on `chore@e7499ef` (origin): `--trust-existing-bams` — diff vs master
  before next landing.
- Both handed sentinels still RUNNING at 22:45 (668b, 674); 10 jobs in queue; scratch 1.6 TB.

# (prev header) 649 complexity-gate landing gating; clip-leg cost plan active — ~22:20 PDT

## Delta ~22:20 — 669-prodgap handoff actioned (Kevin-directed); Item 1 landing

- **669-prodgap handoff received + ACK'd** (their 2135Z note, my ACK 0510Z): three
  correct-cost items now mine. correct = 99% of the per-sample bill; production pool 648k
  junctions (support-1, no size cap) → hp_ed DP monoculture (85% cum).
- **Item 1 IN FLIGHT: branch `feat/land-complexity-gate`** (worktree realigner-triage) =
  master + cherry-picks `5278610` (305daff: pool-admission complexity gate + deterministic
  rescue candidates + threaded-output determinism) + `70e009c` (e34500d: gate ON default,
  --complexity-alpha 0.01) + §4c gate-reconciliation doc. Cherry-picked NOT merged (their
  branch base = stale pre-rebase b028e35). Targeted 96/96. Interaction verified: complexity
  gate calls effective_information_bits directly — period cap does NOT affect it.
  **Full gate running** (nohup pid `scratchpad/.gate_pid`, log
  `scratchpad/suite_complexity_gate.log`). RESUME: SUITE_RC=0 → ff master + push (same
  pattern), notify 669+umi (runbook: prescan adds --junction-max-size 10000, gate
  default-on). Non-zero → only NEW failures vs 2,244/0 matter.
- **Clip-leg cost plan (Kevin: commence)**: Phase 0 profile superseded by 669's
  production-shape profile + code confirmation (5' leg = linear scan over
  candidate_junctions per read, triage.py/_terminal_peel_rescue). 5k-read local profile
  timed out at 10 min (M1) — itself the datapoint. NEXT: Phase 1 (leg-entry bisect
  preconditions + pA/pt skip + threshold tuning) + Phase 2 (bisect-indexed candidate
  lookup in rescue — also speeds Cat3-in-correct) + Phase 3 memoization, on a branch after
  the complexity gate lands. Acceptance: ≤2× junction leg on drs100k, gold-window recall
  unchanged.
- **644k round-3 (H2 `14287945`, n1814): cdna_wt1 classify DONE 580s rc=0** (compute node
  fixed the OOM); realign running; drs clip-legs LAST (h_rt-bounded). Sentinel
  `.644k_round3_rc`, STATUS3.
- 669 Items 2 (numba default-on) + 3 (7–20× residual: uLTRA arm first, SCG calibration
  second) sequenced behind clip-leg work — task list #12.

# (prev) HANDOFF — period-gate fix landed — 2026-08-11 ~20:30 PDT

## Delta ~21:00 — b1b033f LANDED+PUSHED (suite 2,244/0); 644k round-1 numbers in; round-2 queued

- **origin/master = b1b033f** (period-gate fix; ff + push verified). H2 644k clone updated
  to b1b033f (git 1.8 on H2: no `-C`, use cd + `git fetch --depth 1 && git reset --hard
  FETCH_HEAD`).
- **644k round-1 (drs100k = 100k-read mm2-only run-all output, 98,991 primaries):**
  Station B classify: 5,032 bypass (5.1%) / 93,959 triaged (94.9% — placeholder thresholds
  confirm the "tune on real corpora" note); junction leg 12,944; refiner moved 4,662;
  **hp_ed accepted 1,478 = 1.49% of primaries modified**. Station C on input consensus:
  1,555 junctions (230 annotated); non-annotated verdicts 38 admit / 341 review / 946
  demote. C modifies nothing.
- **644k round-2 QUEUED: H2 job `14287597`** (compute node; login-node MemoryError killed
  the cDNA legs in round 1): period-FIXED clip-legs timing on drs100k (before = killed at
  >14 min), full cdna_wt1 leg set, mod-censuses, pool-gates. Sentinel
  `/u/scratch/k/kevinroy/644k_blast_radius/.644k_round2_rc`; progress `STATUS2`; log
  `round2_qsub.log`. NB cdna_wt1 uses --Scer: reporter-construct contigs absent from the
  bundled genome — construct reads skipped in census; if triage KeyErrors on the missing
  contig, rerun with the 587 run's own genome+construct fasta (find in 587_ms2_cdna).
- **⚠ ROUND-2 FINDING (~21:40): the period fix does NOT rescue the clip legs on plain
  DRS.** Fixed clip-legs ran 48+ min at 100% CPU on 100k reads without finishing (junction
  leg: 2 min). Two distinct diseases: (1) reporter-construct/resolver = period-driven —
  FIXED by b1b033f (the UMI agent's actual case; resolver at align has no clip legs);
  (2) clip legs = VOLUME-driven — placeholder classifier triages 94.9% of reads, so the
  legs run per-clip rescue/DP on ~90k reads. Clip legs stay OFF for beta (already the
  guidance). Engineering next: classifier threshold tuning + Fix 2 memoization + cheap
  leg pre-filters.
- Round 2 (`14287597`) QDEL'd (clip-legs blocked the serial script); **ROUND 3 = H2 job
  `14287945`**, reordered: cdna_wt1 legs FIRST, drs clip-legs last (h_rt-bounded).
  Sentinel `.644k_round3_rc`, progress `STATUS3`, log `round3_qsub.log`.
- RESUME: poll `ssh hoffman2 'tail -15 /u/scratch/k/kevinroy/644k_blast_radius/STATUS3;
  cat /u/scratch/k/kevinroy/644k_blast_radius/.644k_round3_rc 2>/dev/null'` → done:
  compile cDNA blast radius + clip-leg timing (or its h_rt-truncation, itself a number),
  report to Kevin + UMI inbox. Then Fix 2 (depth memoization) + classifier threshold
  tuning as the clip-leg cost work.

## (prev) Delta ~20:30 — resolver slowdown ROOT-CAUSED + FIXED; gate in flight

- **d7982cc LANDED + PUSHED** (~19:05): origin/master = d7982cc; suite was 2,235/0.
  Stations B+C in run-all; runbook addendum sent (Chanfreau inbox 0210Z).
- **>500× reporter-construct resolver slowdown (UMI agent's report, Kevin escalated):**
  ROOT CAUSE = I_eff estimators blind to tandem periodicity unit >~3 nt (MS2 19-mer array
  scores 76–81 bits, W_max 5000, true in-array discrimination ~0 → full candidate DP ×
  depth → rejected_ambiguous the expensive way). Same disease as rDNA; leaks into clip
  legs on plain DRS (measured: >14 min vs 2 min junction leg on 100k reads — that step was
  killed, PID verified first). Masking = band-aid AND silences the MS2 splicing reporters'
  own signal.
- **FIX committed `b1b033f`** (branch fix/overhang-period-gate, worktree realigner-triage):
  `min_self_match_period` (shift≤32, ≥80% identity, ≥12 bp overlap) caps W_max at
  period−1. MS2×3 5000→18; 8-mer→7; polyA/(AG)n refuse via the general mechanism; random/
  exonic untouched; 10%-error arrays still caught. Targeted suites 101/101 (one 641 test
  updated: alpha-monotonicity exempts period-capped seqs). UMI agent told to HOLD masking
  (inbox 0320Z).
- **IN FLIGHT: full gate on b1b033f** — nohup pid in `scratchpad/.gate_pid` (78757), log
  `scratchpad/suite_period_gate.log` ends SUITE_RC=. RESUME: SUITE_RC=0 → `cd
  ~/work/rectify_worktrees/land-master && git merge --ff-only fix/overhang-period-gate &&
  git push origin master` (beta-critical, Kevin: "proceed with haste"), then on H2: `git -C
  /u/scratch/k/kevinroy/644k_blast_radius/rectify_src pull` and rerun the clip-legs step
  (the line in 644k_blast_radius.sh) for the before/after timing; report numbers to Kevin +
  UMI inbox. Non-zero → only NEW failures matter (tree was 2,235/0 at d7982cc).
- **644k blast-radius**: drs100k classify+realign done rc=0; clip-legs killed (pathology
  measured); mod-census failed harmlessly (no triaged BAM); cdna_wt1 legs running. Results
  land in `/u/scratch/k/kevinroy/644k_blast_radius/` (STATUS + per-substrate dirs).
- **Fix 2 (depth memoization)**: designed, NOT started — next after Fix 1 lands.

# (prev) HANDOFF — Stations B+C built (d7982cc), gate running — 2026-08-11 ~18:15 PDT

## Delta ~18:15 — Kevin: entire B and C in beta run-all; BUILT, gate in flight

- **origin/master = 4a304f0 PUSHED** (~17:40): clip legs + UMI branch (RN fix, resolver CLI
  + substitution, tag-test adjudication). Suite was FULLY green: 2,225/0. UMI + 641 + rbrowse
  lines notified (Chanfreau inbox 0030Z).
- **Station C v0 BUILT + Station B wired into run-all** — branch `feat/station-bc-beta` @
  `d7982cc` (worktree `rectify_worktrees/realigner-triage`), one commit on top of 4a304f0:
  - `rectify/core/consensus/station_c.py` — pool-gate: junction census (min-anchor 8,
    ambiguity-canonicalised), 644h short-side overhang q ('='-decode), 644i flags (GFF
    repeat classes + bundled `R64_selfhomology.bed`, force-added past the genomes/ ignore),
    two-track verdicts (q_canon 40 / q_noncanon 80 / support≥2), REPORT-ONLY TSV+JSON.
  - `rectify pool-gate` CLI (registered in cli.py); `rectify triage --clip-legs` flag.
  - run-all: `_run_station_bc` fail-soft post-stage — C default ON (`--no-pool-gate`),
    B opt-in (`--triage`, `--triage-clip-legs`); v0 limitation documented: B emits review
    artifacts (triage.tsv + triaged BAM), analyze tables not rewritten.
  - Targeted tests 30/30 (10 new Station-C incl. synthetic mini-genome verdict classes,
    bundled-track rDNA coverage, CLI smoke; 9 clip-leg; 11 triage).
- **IN FLIGHT: full gate** on the d7982cc tree — NB two earlier harness-tracked attempts
  were externally stopped (second at 85%, zero failures visible); now running as plain
  nohup `nice -10`, **PID in `scratchpad/.gate_pid` (62099)**, log
  `scratchpad/suite_station_bc.log`, ends SUITE_RC=. Poll:
  `tail -2 …/suite_station_bc.log; kill -0 $(cat …/.gate_pid) && echo alive`.
  RESUME: SUITE_RC=0 → `cd ~/work/rectify_worktrees/land-master && git merge --ff-only
  feat/station-bc-beta && git push origin master` (Kevin's beta directive; same push
  pattern as 4a304f0), then notify Chanfreau inbox (beta runbook addendum: pool_gate.tsv
  artifact + --triage/--triage-clip-legs flags) and update memory. SUITE_RC!=0 → only NEW
  failures matter (tree was 2,225/0 at 4a304f0); fix before landing. Gate dead with no
  SUITE_RC → relaunch the same nohup line from the worktree.
- 644j smoke on H2 + its watcher: still pending (see earlier delta).

# (prev) HANDOFF — beta landing: clip legs + UMI branch at 4a304f0 — 2026-08-11 ~17:00 PDT

## Delta ~17:00 — Kevin prioritized clip-leg testing for the beta; landing tree assembled

- **Clip-leg implementation recovered from the killed agent** (worktree
  `.claude/worktrees/agent-ab99769edf5ec64f8`, real branch `feat/triage-clip-legs` @
  `5be89dd` — the auto worktree ref stayed at base, don't be fooled). Reviewed: backward-
  compatible `resolve_read(sides=)`, `TriagePolicy.clip_legs_enable` OFF-default, both legs
  on read COPIES, strict hp_ed adjudication, refusals first-class, two-sided seam stubbed.
  Targeted suites **55/55** (9 clip-leg + 11 triage + 35 resolver).
- **UMI agent's inbox notes actioned + archived**: f20a8e6 adjudicated (test updated —
  overwrite is scoped to `_CDNA_COMMENT_TAGS`, old test premise stale); RN-order consensus
  group-split FIXED (f652331 + monotonicity guard; explains 657's 1.5 primaries/read);
  resolver CLI `rectify align --junction-aligners overhang_resolver` (idempotent mtime-
  gated reuse) + family map + ARM SUBSTITUTION (965394b — no duplicate arm).
- **Landing tree: local master = `4a304f0`** = 2701f3b + merge chore@965394b (55298fb) +
  merge feat/triage-clip-legs (4a304f0). Both merges conflict-free. **Full gate in flight**
  (log `scratchpad/suite_beta_landing.log`, ends SUITE_RC=; expect FULLY green — the old
  tag-restoration red is adjudicated on this tree).
- **run-all fill-in-the-blanks (Kevin's ask): mostly ALREADY EXISTS** — sample-level
  resume-on-restart skip-if-output-exists gates at every stage (alignment multialigned/
  per-aligner checkpoints, cDNA prep, correction manifest/TSV, provenance reuse gate); the
  resolver CLI post-pass is mtime-idempotent. Gap: run-all doesn't yet invoke the resolver
  post-pass itself (production scripts use the align-command path) — follow-up wiring, not
  a beta blocker.
- **RESUME**: `tail -3 scratchpad/suite_beta_landing.log` → SUITE_RC=0 (or only NEW
  failures absent) → `git push origin 4a304f0:master` (beta pulls origin; Kevin's beta
  directive + prior push pattern authorize), notify umi-perf-accuracy (their inbox) +
  641-impl, add clip_legs_enable note to the beta runbook thread, clear marker. Suite red →
  diagnose; the tree is 3 independent lines, bisect by leaving out the last merge.
- 644j smoke on H2 still in flight (job 14285400, sentinel `.644j_rc_smoke` — see below).

# (prev) HANDOFF — Re-aligner Phase 1: 644h measured; landing gates green so far — 2026-08-11 ~12:50 PDT

## Delta ~14:00 — all fronts launched (supersedes "nothing in flight" below)

- **644i DONE + committed** (`5df199e` local master, NOT pushed): Station-C repeat flag
  measured — annotation flag free (−14% junk, 0 gold cost), genome self-homology track
  adopted (−37% junk at 2–3 demoted singleton gold; `minimap2 -DP -k19 -w19 -m200 G G`,
  seconds), variant-multiplicity REJECTED (conflates repeat spray with contested loci).
  `dev/STATIONC_REPEAT_FLAG_644I_20260811.md`.
- **Docs-Agent handoff written + committed** (`dev/DOCS_AGENT_HANDOFF_20260811.md`):
  two-pass plan, pass-1 scope, Kevin's Rectify_readme_changes_KR_proposed.txt reconciliation.
- **Resolver production-go relayed** (Kevin's directive): durable note in Chanfreau inbox
  (2050Z, to the cDNA-UMI/DRS line) — yeast DRS GO with v5.1 defaults +
  RECTIFY_SKIP_REGIONS=yeast-rdna; cDNA shadow-first (one library, XB census back to us).
  Live SendMessage failed (their Remote Control session gone); inbox note is the record.
- **IN FLIGHT — two subagents** (will notify on completion):
  1. Clip-leg wiring: branch `feat/triage-clip-legs` from local master `5df199e`, own
     worktree; plan doc dev/CLIP_LEG_WIRING_PLAN_20260811.md; default-off gate; full-suite
     gate (only tag-restoration red acceptable). Deliverable = branch, no merge.
  2. 644j production-fenced phase-1 contest on H2 [status ~14:50]: SMOKE JOB `14285400`
     RUNNING on `pod_smp.q@n1812` since 13:33 PT (8 slots; enumerate → subset → pool →
     smoke correct×2 --streaming → merge → census). Durable: `644_accept/644j/PLAN.md`
     (job appends checkpoints), smoke sentinel `644_accept/644j/.644j_rc_smoke`, local
     skeleton dev/PHASE1_FENCED_CONTEST_644J_20260811.md (uncommitted). A local watcher is
     armed on the sentinel; if it dies, poll: `ssh hoffman2 'cat
     /u/scratch/k/kevinroy/644_accept/644j/.644j_rc_smoke 2>/dev/null || qstat -j 14285400
     | head -3'` — rc=0 → verify smoke output then proceed per PLAN.md to the full run;
     nonzero → read the job log named in PLAN.md. Comparison row = 644g flat contest
     (junk 94% win-or-tie).
- Local master (`5df199e`) is 2 commits ahead of origin (644i + docs handoff) — push policy:
  Kevin approved push-after-suite-green for the LANDING; these docs commits ride the next
  approved push.
- RESUME if session dies: check both agents' plan docs + the H2 sentinel; read
  `dev/DOCS_AGENT_HANDOFF_20260811.md` before any docs work; 644i/644h JSONs in scratchpad +
  H2 t3/full/.

## FINAL ~13:30 — LANDED AND PUSHED. `origin/master` = `03bc44f`. Phase-1 landing sequence COMPLETE.

- Suite on `03bc44f`: **2,210 passed / 1 failed** — the failure is exactly and only the
  pre-existing tag-restoration one (ONT-cDNA line's intent call, filed 1940Z).
- Kevin approved the push in-session ("push after suite green"); 641-impl pushed the
  resolver leg (255a06d→345f6ca) with Kevin's word in THEIR session; I pushed
  345f6ca→03bc44f (ff, verified post-push). Public repo carries the full Re-aligner line.
- 641-impl notified (SendMessage); [[617]] LOO unblocked from origin refs.
- **Next phase (nothing in flight):** (1) wire triage clip legs against the landed resolver
  API (landscape §4b division of labor); (2) prototype the Station-C repeat-context flag +
  phase-2 scorer (644h design); (3) production-fenced phase-1 contest (644g follow-up (i));
  (4) triage-policy tuning on the upf1Δ gold windows; (5) fix/adjudicate the
  tag-restoration test (theirs). The worktree `rectify_worktrees/land-master` sits on master
  = `03bc44f`; `feat/realigner-triage` and `feat/overhang-resolver-641` are both fully
  merged (branch refs can be deleted at leisure).

## (superseded) Delta ~13:15 — BOTH LANDINGS DONE LOCALLY

641-impl answered by live SendMessage: they rebased `feat/overhang-resolver-641` onto the
chore tip `81ea36b` → `345f6ca` (pushed to origin as the branch), handed me the ref-update.
Verified everything (ancestry, ff-ability from origin/master@255a06d, 13 commits, content ≡
my abandoned `e605791` merge modulo exactly the 3 chore-tip commits). Then:

- **local `master` → `345f6ca`** (resolver landed; my `e605791` merge discarded, theirs
  carries the same content + the chore tip).
- **`feat/realigner-triage`@`425c2b7` merged on top → `master` = `03bc44f`** (clean; zero
  both-sides files; preview branch deleted).
- Noncanon control already PASSED on this merged content (v2 ≡ noarb ≡ mm2, fjFDR 0/600).
- **IN FLIGHT: full suite on `03bc44f`** (log `scratchpad/suite_final_master.log`, ends
  `SUITE_RC=`). Only acceptable failure = the pre-existing tag-restoration one (filed to
  ONT-cDNA line; attributions differ: me f20a8e6, 641-impl aa9fd59; nobody bisected).
- **ORIGIN PUSH NOT DONE — Kevin's call** (publishes ~400 commits to public GitHub;
  question asked in chat). If Kevin approves: `git push origin 03bc44f:master` (ff from
  255a06d). 641-impl notified of all of this by SendMessage; their [[617]] LOO hold lifts
  per their needs (local/cluster refs are landed; origin waits on Kevin).
- RESUME if session dies: check suite log; green (or only tag-restoration red) → report,
  ask Kevin re origin push, then clip-leg wiring (landscape §4b) + repeat-flag prototype +
  production-fenced phase-1 contest.

## Delta since 11:30 (superseded by the above — coordination freeze is OVER)

- **Full suite at `e605791` (resolver merge): REGRESSION-FREE.** 2,128 passed, 1 failed —
  and the failure (`test_consensus_tag_restoration.py::…without_overwriting`) fails
  IDENTICALLY at pre-merge `1f21a3f`: bisected, attributed to `f20a8e6` (uLTRA tag-overwrite
  fix vs the 6813bd1-era no-overwrite assertion). Filed to the ONT-cDNA line
  (Chanfreau inbox note 1940Z); suite-verdict follow-up to 641-impl sent (1935Z). Their
  keep-or-unwind reply still pending (my 1825Z ping unread as of ~12:35); **master frozen at
  `e605791`, do not move it.**
- **644h amended per Kevin (commit `425c2b7` on triage branch):** Station C gets a
  direct-repeat flag (rDNA, CDS tandem arrays, CUP1/ENA, Ty/LTR, subtelomeres) instead of a
  full uniqueness leg — measured: the ≥80-bit non-canonical junk collapses to rDNA + ~251
  loci; rDNA+top-100 = 85%. Flag DEMOTES (orthogonal evidence), never discards (one Gould
  gold sits inside a junk cluster). Architecture answer to Kevin: correct stays untouched —
  resolver = align arm; triage/Station C = consensus module; only two default-dark gates in
  shared splice/ utilities (`RECTIFY_OVERHANG_INFO_GATE`, `RECTIFY_SKIP_REGIONS`);
  correct_command's refiner call keeps `motif_blind=False`.
- **Triage-merge PREVIEW built:** temp branch `tmp/triage-merge-preview` @ `a7014ed`
  (= triage `425c2b7` merged onto `e605791`; clean, zero both-sides files, master untouched).
- **Noncanon control RERUN on the preview tree: PASSED — flattening GONE** (v2 ≡ noarb ≡ mm2
  every cell; fjFDR 0/600; recovery 0.708 all arms). Integration §4 step-3 control satisfied
  on the actual landed content. Output: scratchpad `rearb_control_out/trade_curve.json`.
- **IN FLIGHT: full suite on `a7014ed`** (background task; log
  `scratchpad/suite_preview_merge.log`, ends with `SUITE_RC=`). Expect ~25 min; the ONLY
  acceptable failure is the same pre-existing tag-restoration one.

## Resume (updated)

1. `tail -3 …/scratchpad/suite_preview_merge.log` → `SUITE_RC=1` with ONLY the
   tag-restoration failure = preview merge is regression-free (2,209+ passed expected:
   2,128 + 81 triage-branch tests). Any OTHER failure → diagnose before any landing.
2. Check `.claude/inbox/` for 641-impl's reply:
   - **keep `e605791`** → `cd ~/work/rectify_worktrees/land-master && git branch -f master
     a7014ed && git checkout master && git branch -d tmp/triage-merge-preview` — BUT only
     re-point master to `a7014ed` if the preview suite was clean; otherwise land `e605791`
     only and fix first.
   - **unwind** → `git branch -f master 1f21a3f`; wait for their landing; then re-merge
     triage onto THEIR master (preview results transfer if content identical; re-run suite).
3. Then: wire clip legs (landscape §4b), phase-2 repeat-flag prototype, production-fenced
   phase-1 contest (644g follow-up (i)).

# (previous header) HANDOFF — Re-aligner Phase 1: 644h measured; resolver-landing COORDINATION PENDING — 2026-08-11 ~11:30 PDT

**Owner:** rectify-realigner agent. Program docs: `dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md`
(now with §4b two-sided-enumeration design) + the 641→Re-aligner integration package
`~/work/UCLA/Chanfreau_Lab/planning/644h_realigner_integration.md` (the governing spec).
Prior session state: see git history of this file (`1f21a3f`).

## Done (this session, 2026-08-11)

- **644h MEASURED + committed** (`11dc0c3` on `feat/realigner-triage`): the PI's phase-2
  short-exon-side overhang-quality likelihood, scored against the 644f/644g census
  (tool `scripts/benchmark/644h_phase2_overhang_likelihood.py`; data H2
  `644_accept/t3/full/644h_phase2_overhang.json`; finding `dev/PHASE2_OVERHANG_644H_20260811.md`).
  **Headline: canonical track upgrades to ~79% precision at 11/12 gold (support≥2 + q≥40,
  AUC 0.91); the non-canonical (Gould) track CANNOT be purified by overhang quality alone**
  (AUC 0.77; 1,357 junk carry ≥80-bit clean overhangs, rDNA only ~12% of them — systematic
  copy-misplacement junk is locally perfect). Phase-2 scorer must be overhang-quality ×
  placement-UNIQUENESS × recurrence; the T6-style uniqueness leg is the next build.
  Implementation trap caught by the smoke: the mapPacBio BAM's SEQ is `=`-compressed —
  decode against the reference or all scores are 0.
- **641 inbox messages actioned + archived** (one-sided→two-sided correction folded into
  landscape doc §4b; integration handoff read; both in `.claude/inbox/.read/`).
- **Sherlock cross-check `38619572` verified COMPLETED 0:0** (integration §4 step-1 precondition).
- **Resolver merged onto local master as `e605791`** (merge of
  `feat/overhang-resolver-641`@`124c84a` onto `1f21a3f`; conflict-free; both-sides check:
  numba hp_ed fast path + region_skip coexist in `splice_aware_5prime.py`). Done under
  Kevin's "everything is in your hands" — THEN Kevin flagged 641-impl is in-flight doing the
  landing themselves → **I stopped, detached my worktree (`rectify_worktrees/land-master`),
  freed the master ref, and pinged their inbox** (Chanfreau inbox,
  `2026-08-11T1825Z__from-rectify-realigner__master-already-carries-resolver-merge...`):
  their call = keep `e605791` or unwind to `1f21a3f`. **Do NOT move master until they reply.**

## Verified

- 88 resolver-suite tests pass at `e605791` (overhang_resolver / informativeness /
  region_skip / splice_site_index / calibrate_overhang_header).
- Triage-merge preview vs resolver-bearing master: **ZERO both-sides files** since merge-base
  `999ceb5` (verified via dumped file lists + `/usr/bin/comm`, not git pipes — the false-zero
  shell trap).
- 644h frames reproduce 644f/644g exactly (37 gold / 5,414 junk; surv1 35/5,068; surv2 17/1,737).
- NOT VERIFIED yet: full `-m "not slow"` suite at `e605791` — in flight (84% at last check,
  0 failures so far; log `…/scratchpad/suite_land_master.log`, ends with `SUITE_RC=`).

## Open

1. **Full-suite verdict at `e605791`** (in flight, this session's gate).
2. **641-impl's keep-or-unwind reply** (their inbox ping; master frozen until then).
3. **Triage merge** (`feat/realigner-triage`@`11dc0c3` → master) AFTER resolver landing
   confirmed; then full suite again; then wire the clip legs against the landed resolver API.
4. **Noncanon control rerun on the landed tree** (driver committed:
   `scripts/benchmark/noncanon_sim/rearb_noncanon_control.py`; NB its `TREE`/`WORK` consts
   point at a dead scratchpad — repoint or parametrize; clear stale WORK outputs first).
5. **Phase-2 uniqueness leg** (the 644h consequence): prototype placement-uniqueness for the
   non-canonical track (overhang+context genome-wide uniqueness; T6 decoy check generalized).
6. Production-fenced phase-1 contest (644g follow-up (i)) — still open, unstarted.
7. `[[617]]` leave-one-out — fast-arm-design owns; waits on landing.

## Resume (concrete)

1. `tail -5 …/scratchpad/suite_land_master.log` (this session's scratchpad:
   `/private/tmp/claude-501/-Users-kevinroy-work-rectify/a5790346-5092-4428-9774-4f6d83e50315/scratchpad/`).
   `SUITE_RC=0` → suite green at `e605791`; append that to the 641-impl ping thread.
   Non-zero → read failures; if they implicate the merge, say so in the ping and recommend unwind.
2. Check `.claude/inbox/` for 641-impl's reply.
   - **Keep `e605791`** → proceed: `cd ~/work/rectify_worktrees/land-master && git checkout master`
     (should still be `e605791`; if they moved it, take THEIR ref), `git merge --no-ff
     feat/realigner-triage`, full suite, then open items 4–5.
   - **Unwind** → `git branch -f master 1f21a3f`, stand down until they announce their landing,
     then merge triage onto theirs.
3. H2 644h data is durable (`644_accept/t3/full/644h_phase2_overhang.json`, sentinel
   `.644h_rc` = 0). Local copy: scratchpad `644h.json`.
4. Tasks #4/#5/#6 in the session task list mirror this section.

## Files

- 644h: `dev/PHASE2_OVERHANG_644H_20260811.md` + `scripts/benchmark/644h_phase2_overhang_likelihood.py` (branch `feat/realigner-triage`@`11dc0c3`)
- Design record: `dev/REALIGNER_LANDSCAPE_AND_PATH_20260809.md` §4b (two-sided enumeration, corrected)
- Integration spec: `~/work/UCLA/Chanfreau_Lab/planning/644h_realigner_integration.md`
- Landing worktree: `~/work/rectify_worktrees/land-master` (DETACHED at `e605791`)
- Coordination ping: Chanfreau `.claude/inbox/2026-08-11T1825Z__from-rectify-realigner__…`
