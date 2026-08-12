# HANDOFF — 679 cDNA trim bug (F1 root cause) · DRS handed to unit 682 · 313 GB move ~1/3 done — CURRENT 2026-08-11 ~23:00 PDT

## Delta ~03:30 — overnight autonomous run: impossible-junction guard LANDED; three peer-driven fixes

Kevin signed off ~02:15 with "proceed autonomously through the night and coordinate with rbrowse
and DRS". **The 676 uBAM DELETION HAS NOT BEEN DONE and will not be — it needs him to approve the
file list.** Everything below is additive or test-covered.

- **`origin/master` = `ea4401e`** (full suite **2301 passed / 0 failed**), plus `d14cd22` queued
  behind a running suite. Landed tonight:
  - **`d0e3a0f` impossible-intron guard.** 10 kb max reportable intron for yeast; a winning
    alignment with a larger N-op is TRUNCATED at that junction (rest of query → one soft clip,
    query length preserved exactly), tagged **`Xn:i:<bp>`**, counted as
    `stats['impossible_intron_truncated']`. Env `RECTIFY_MAX_REPORTABLE_INTRON`; 0 disables.
    Runs inside consensus selection, so a plain `rectify consensus ALIGNER:BAM ...` re-run picks
    it up **with no realignment**. Ships WITH the `ambiguity_window` bounds fix deliberately —
    alone that would hide the defect.
  - **`ae69e79` write-time invariant** (adopted from 668-drs-arm's argument): no alignment may
    run past its contig end, in `_validate_bam_sample`, naming up to 5 offenders.
  - **`d14cd22` pool-gate `-o` silently overwrote versioned outputs.** `Path.with_suffix`
    REPLACES a suffix, so `sample.v1` and `sample.v2` both wrote `sample.pool_gate.tsv` — two
    runs, one file, no warning. A "one-line help fix" that was a real defect underneath.
- **`dev/683_G_ARBITER_EXPLORED_NOT_NEEDED.md`** filed per Kevin — patch stays parked.

### The evidence that decided the design (both peers, both independently measured)

- **cDNA cohort survival (27 samples, 19.28 M reads): they survive COMPLETELY.** 247 past-contig
  reads in → **247 rows in `corrected_reads.tsv` → 179 CLUSTER ANCHORS** (called 3' ends).
  `cdna-analyze` filters **zero**. That is what made invariant-only disproportionate.
- **DRS per-arm census (400k reads/arm):** minimap2 **0/0**, overhang_resolver **0/0** (Station A
  is NOT a contributor), uLTRA 0 past-contig but **394 N-ops >10 kb** (max 346 kb), deSALT **36
  past-contig / 2,067 >10 kb**, multialigned 28 / 860. **DRS runs 2–3× worse than cDNA**
  (0.215 % vs 0.067 %). Critically: **a past-contig check alone would catch NONE of uLTRA's 394**
  — which is why the guard keys on N-op size OR past-contig, two independent conditions.

### 🔴 OPEN RESIDUAL — do not read "guard landed" as "class closed"

The 10 kb bound catches >10 kb; the **2–10 kb band passes untouched** (DRS: 1,013 vs 860 per 400k
⇒ ~15 % of the >2 kb population survives). It **cannot** be closed by lowering the threshold —
minimap2's max N-op on a clean arm is **5,064 bp**, so below ~5.1 kb you clip legitimate output.
Needs a different discriminator (per-read disagreement with the minimap2 arm). Unfixed.

### Also open / flagged

- **`--junction-pool-max-intron-len` is exposed ONLY on `rectify consensus`**, not on `align` or
  `run-all` — i.e. unreachable in the production path. Plumbing gap, not closed tonight.
- **`pre-681 data may be ~24× worse`** (1.593 % vs 0.067 % N-ops >5 kb). Mechanism unproven;
  cDNA session's controlled test `684e` is running. **Do NOT build the 24× into any rationale
  until it lands.** If it holds, prioritise re-inspection by whether a run predates `b3a8c35`.
- **Both peers plan a consensus-only re-run** of completed samples (arm BAMs are on disk) before
  drawing Station B/C conclusions. Coordinate so they don't duplicate the mechanics.
- My earlier claim that `pool-gate` "does not exist in this tree" was **wrong about master** —
  true only of my branch, which forked before `d7982cc`. Branch now merged up to master.

## Delta ~02:10 — 🔴 NEW CORRECTNESS FINDING (peer-reported, INDEPENDENTLY VERIFIED): consensus selects physically impossible alignments

**Reported by the cDNA session (`planning/684c_impossible_junctions_desalt_ultra.md`); I
re-measured it myself on a different slice rather than take it on trust.**

My verification — `684_p1cdna_1M/WT_BY4742_rep1/align/*.multialigned.bam`, 400,001 primary reads:

| check | count |
|---|---:|
| alignments ending **PAST the contig end** | **3** |
| N-op **> 10 kb** | **268** (0.067 %) |
| **max N-op** | **261,350 bp** |

The longest annotated intron in *S. cerevisiae* is **~1 kb**. Example CIGARs:
- `r030_7056` chrIV: `…190957N…` and ends at **1,593,636 vs contig length 1,531,933** — 61,703 bp
  past the end of the chromosome.
- `r036_7091` chrV: `…150801N…99483N…102507N…` — three separate ~100–150 kb "introns" in one read.

Peer's per-arm census (3 samples, ~2.7 M reads/arm) attributes it: **minimap2 = 0 past-end, 0
>10 kb**; uLTRA 2,372 >10 kb; deSALT 47 past-end / 6,022 >10 kb; and `multialigned` inherits
36 / 1,886. Read-level proof in their doc shows minimap2 and uLTRA agreeing in-bounds while
deSALT appends a 192,339 bp N-op running off the chromosome — **and the consensus picks deSALT.**
Whatever the scorer optimises, "spans more query" is beating "is physically possible".

**Two distinct defects:**
1. **Consensus selection** admits a reference span past the contig end and N-ops far beyond the
   organism's max intron. `-G` constrains the minimap2 arm; nothing equivalent constrains
   uLTRA/deSALT or the selection step.
2. **`overhang_informativeness.ambiguity_window:402` CRASHES** on such input
   (`IndexError: string index out of range` at `g[end - 1 - l_amb]`) — that is how the peer
   found it, via `rectify pool-gate`.

⚠️ **Fix BOTH or neither** — a bounds check alone would silence the crash and HIDE defect 1.

**NOT [[681]]:** the trim acts on the consensus sequence pre-alignment, and minimap2 on the same
trimmed input is clean. Contigs are not mismatched (17/17 identical name+length).

**NOT YET MEASURED — the number that decides nuisance vs. correctness:** how many of these
survive `cdna-analyze`'s own filters into a junction table. The peer offered to run it on 684
output; I have accepted.

**Design fork worth Kevin's call before implementing:** rejecting an arm outright vs. flagging +
counting it. A hard reject is the correct-looking move but it changes aligner-selection semantics
for every dataset, and mis-set it would discard legitimate alignments — the same shape of trap as
the candidate ceiling. Detection (an invariant + counter at consensus write time, same shape as
`--require-aligners`) is safe to land unilaterally; policy is not.

## Delta ~01:50 — peer review from the cDNA session; two of my statements corrected

- **`684` runs NO resolver.** Its align line is `--aligners minimap2 --junction-aligners uLTRA
  deSALT` — the resolver is not in the panel, so **684 will never produce an
  `overhang_resolver.stats.json`.** My Resume step 3 asked for one from it; that ask is void.
- **🔴 CORRECTION — the 683 G-overhang arbiter has ZERO MEASURED IMPACT.** I earlier relayed to
  Kevin that the 1-bp bridge-G exposure was "frequent, not rare" (~34% of molecules reach the
  5' rescue arbiter), citing the 681 tail. The cDNA session's **683 CP2 supersedes that**:
  **0 of 23,663 suspect bridge-G clips fall within 10 bp of a 3'SS**, because the bridge G sits
  at the TSS while the rescue arbitrates 3' splice sites.
  ⚠️ **Both numbers are true and measure different things — do not conflate them:** ~34% is the
  1-bp G *PREVALENCE* across molecules; 0/23,663 is how often one *REACHES THE ARBITER*. The
  first was stated before the second was measured. The guard is correct and well-gated
  but fires on nothing. Parked for Kevin's decision; patch + design at Chanfreau
  `planning/683_arbiter.patch` · `683_test_bridge_g_arbiter.py` · `683_g_aware_5p_clip_arbiter.md`
  (reapply with `git apply` — deliberately extracted, NOT lost).
- **Ceiling calibration corrected (`b153303`, docs/test only).** The peer — who produced the 361
  measurement — pointed out that **361 is the PRE-FIX PATHOLOGY, not a healthy maximum**, and
  that my rationale invited exactly the wrong future edit. Both CP6 figures are now recorded in
  `ResolverConfig` and the test: pre-fix 17 clips/6,137 candidates = **361/clip**; post-fix
  9 clips/60 = **6.7/clip**. The `> 361` test is a **deliberately conservative floor, NOT a
  calibration target** — tightening toward the post-fix mean would refuse legitimate clips on
  repetitive/custom references, the case the ceiling exists for.
  **Histogram ask WITHDRAWN**: pre-fix percentiles describe a distribution that no longer
  exists, and healthy data cannot calibrate a bound on the pathological tail. The ceiling stays
  a safety BOUND, not a tuned parameter; 2000 stands.
  Also recorded: the guard bounds **COUNT, not TIME** — 2000 candidates is ~24 ms/clip with
  numba but **~7 s/clip pure-Python**, so a numba-less run on a pathological contig is bounded
  yet still slow. **Prefer enabling numba over lowering the ceiling.**
- **Gap (a) checked by the peer, not bitten:** every completed 684 sample has all three arm BAMs
  + `.multialigned.bam`. The 36 "not found at" log lines are OLD pre-PATH-fix runs (an SGE job
  does not inherit the conda env). That failure was loud only because ALL THREE dropped; with
  one missing it would have silently become a 2-aligner panel at exit 0. They will add
  `--require-aligners` on the next wave.
- **My process error, recorded:** I staged their `HANDOFF_681_*.md` while they were mid-thread,
  which pushed them into defensively extracting 683 to a patch file. In a shared working tree,
  stage only your OWN paths — being explicit is not sufficient if the path is someone else's.

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

**Overnight autonomous run (Kevin signed off ~02:15: "proceed autonomously... coordinate with
rbrowse and DRS"). 🔴 THE 676 uBAM DELETION WAS NOT DONE and must not be until Kevin approves the
file list.** Everything landed is additive or test-covered.

- **`origin/master` = `cdf4bd9`**, full suite **2311 passed / 0 failed**. Tonight's landings:
  | commit | what |
  |---|---|
  | `b3a8c35` | cDNA consensus trim fix (from the cDNA session) — `XQ==0` 52.5→13.5 %, clip 65.8→3.8 nt |
  | `4533de5` | `--require-aligners` + `DROPPED-ALIGNER` summary; resolver `threads` no-op now warns |
  | `24c7805` | per-clip candidate ceiling + `refused_candidate_blowup` |
  | `b153303` | ceiling calibration: 361 is the PRE-FIX pathology, a floor, not a tuning target |
  | `d0e3a0f` | **impossible-intron guard** — 10 kb bound, soft-clip at the junction, `Xn:i:<bp>`, + `ambiguity_window` bounds fix in the same commit |
  | `ae69e79` | write-time invariant: no alignment may run past its contig end |
  | `d14cd22` | pool-gate `-o` silently overwrote versioned outputs (`with_suffix` REPLACES) |
  | `cdf4bd9` | `analyze --manifest` no longer requires the unused positional |
- `dev/683_G_ARBITER_EXPLORED_NOT_NEEDED.md` filed per Kevin; patch stays parked.
- **Storage: 676 move COMPLETE** — 320 GB scratch reclaimed, orphan temp removed, nothing deleted.
- rbrowse notified of the `Xn` tag; DRS and cDNA sessions coordinated throughout.

## Verified

| claim | evidence |
|---|---|
| master green | 2311 passed / 0 failed at `cdf4bd9`; Station C files confirmed present after every push |
| impossible junctions are real | MY re-measurement: 400,001 reads → 3 past-contig, 268 N-op >10 kb, max 261,350 bp |
| **they survive to results** | cDNA cohort, 27 samples / 19.28 M reads: **247 in → 247 rows in `corrected_reads.tsv` → 179 CLUSTER ANCHORS**. `cdna-analyze` filters ZERO |
| which arms | minimap2 **0/0**, overhang_resolver **0/0** (Station A not a contributor), uLTRA 0 past-contig but **394 N-ops >10 kb**, deSALT **36 / 2,067**. DRS 2–3× worse than cDNA |
| the guard works in situ | DRS ran Station C clean at `ea4401e`: 34,122 junctions in 34 s, byte-identical verdicts to their own patch |
| cDNA trim fix live | 684 stage-1 `XQ==0` = 12.4 % (vs 52.5 % pre-fix) |
| DRS ~15× faster w/ numba | 682 smoke frac 0.99966 at 4,927 s ≈ 1 h 22 m/sample vs ~22.1 h |

## Open

1. **🔴 RESIDUAL — "guard landed" ≠ "class closed".** The 10 kb bound misses the **2–10 kb band**
   (DRS: 1,013 vs 860 per 400k ⇒ ~15 % of the >2 kb population). **Cannot** be closed by lowering
   the bound: minimap2's own clean max is **5,064 bp**. Needs a different discriminator —
   per-read disagreement with the minimap2 arm. Tracked in DRS `planning/689a`.
2. **`--junction-pool-max-intron-len` is unreachable in production** — exposed only on
   `rectify consensus`, absent from `align`/`run-all`. "The guard did not fail, it was never
   reachable." Plumbing gap, not closed.
3. **676 (a) deletion — NEEDS KEVIN.** Nothing deleted; data safe on project. Count proof in
   §Resume 1; show him the per-file list first.
4. ✅ **`684e` LANDED — quote 13.6×, RETIRE the 24×.** The 24× was an uncontrolled field
   comparison (587 vs 684, different trees and depths). The controlled test — same 60,000 reads,
   same aligners, same flags, **only the trim differing** — gives:
   | arm | untrimmed N>5 kb | trimmed | reduction |
   |---|---:|---:|---:|
   | minimap2 | **0** | **0** | clean either way (`-G` refuses) |
   | uLTRA | 185 | 17 | 10.9× |
   | deSALT | 252 | 46 | 5.5× |
   | **multialigned** | **341** | **25** | **13.6×** |
   Mechanism confirmed: an untrimmed molecule carries ~150 nt of non-genomic adapter+UMI that an
   unbounded aligner must place somewhere; minimap2's `-G` refuses, uLTRA/deSALT have no such
   bound. **Scope: cDNA only** — `b3a8c35` touches only `rectify/core/cdna/*` and the DRS correct
   path imports nothing from `core.cdna`, so DRS carries no pre/post era (which is what allowed
   the confound to be demonstrated rather than assumed). Prioritise re-inspection of cDNA runs by
   whether they predate `b3a8c35`.
5. **Both peers are doing consensus-only re-runs** of completed samples (arm BAMs on disk, no
   realignment) before drawing Station B/C conclusions. Baselines banked.
6. **✅ SOLVED overnight — and it needs a DESIGN DECISION FROM KEVIN: the `intronlen=40` D→N
   threshold manufactures novel junctions.** Full write-up: **`dev/DN_THRESHOLD_ARTEFACT_20260812.md`**
   (`5ade34a`).
   **What it is:** the panel converts a reference gap to an N-op purely on LENGTH
   (`multi_aligner.py:1016-1018`, `intronlen=40`). An N-op consumes reference without consuming
   query, so it cannot absorb errors *in* the read — and DRS's dominant error mode is deletion
   relative to reference. **Confirmed by histogram: 99.5 % of singleton novel junctions sit at or
   above 40 bp, ZERO between 20 and 29, a step change at exactly 40, peak at 50–54.** A
   splicing-length distribution has no reason to know about 40. These are **deletions relabelled
   as introns** — not splicing, not spurious splicing, not a junction-gate phenomenon.
   **🔴 THE OBVIOUS FIX IS WRONG — I checked it against the bundled annotation.** Raising
   `intronlen` far enough to remove the false population (~100 bp, since it peaks at 50–54 with a
   long tail) would suppress **~107 real spliceosomal introns, ~a third of them**. Length cannot
   separate the populations; they overlap almost completely. (58 of the 60 sub-40 bp annotated
   introns are **tRNA** introns, so `intronlen=40` itself costs only **2 of 318** spliceosomal
   introns and is well chosen AS a threshold.)
   **✅ THE RIGHT LEVER IS THE MOTIF.** Real spliceosomal introns are essentially all GT-AG; the
   false population is **81.5 % non-canonical**. Require a canonical donor/acceptor before a SHORT
   gap may be called an intron. Two placements: **(1) reporting-side** — filter in the junction
   table using Station C's existing `canonical_in_class`, nearly free, reversible, changes no
   alignment; **(2) alignment-side** — at the D/N decision, correct at source but changes CIGARs
   for every dataset and needs its own validation pass.
   **NEITHER IS IMPLEMENTED. Nothing landed overnight touches this — it is Kevin's call.**
   Not done unilaterally because the blast radius is dataset-wide.
   **✅ CLOSED — the driver is PER-LIBRARY BASECALL QUALITY, with a quantitative link:**
   - **Strain background REFUTED:** `wtaa_rep1/2/3` and `wt_rep1/2/3` share the HHY168-AA
     background across different deposits and sit **22× apart** (1,318 vs 29,299 mean). Genotype
     refuted too — `wt_rep1` is the WT-AA *control* at 34k while `wtaa`, also WT-AA, is at 1.3k.
   - **The extremes are the `pod5_skip` batch:** `ysh1_rep1` 178,469 and `ysh1_rep2` 196,643 are
     the corpus's only two `pod5_skip` libraries (two barcodes of one flow cell), already flagged
     INDEPENDENTLY for untailed-read excess; their own third replicate is 23,915.
   - **🔑 The two QC signals scale ~1:1** — junction-fold **7.46× / 8.22×** against untailed-fold
     **7.74× / 7.32×**, agreeing within **12 %** across an 8× range. Two independently-derived
     metrics tracking that closely indicate a **common upstream driver**, not coincidence.
   - **`annotated` stays 262–271 across all nine samples spanning a 155× range of totals** — real
     signal flat, false population moving 155×.
   - **🔴 ACTIONABLE: quarantine `ysh1_rep1` / `ysh1_rep2` from junction-level analysis.** Any
     cross-genotype junction comparison including them compares LIBRARY QUALITY, not biology.
   - **Useful inversion: junction count is itself a cheap library-quality QC metric** — Station C
     already emits it, so the artefact doubles as a diagnostic screen.
   **Arm signatures are OPPOSITE, and this is the dangerous misreading:** the long impossible class
   is uLTRA/deSALT-specific with **minimap2 at zero**; the short D→N class is **panel-wide with
   minimap2 participating equally (~34,000)**. ⚠️ **Switching to minimap2-only — the natural
   mitigation for the long class — carries the short population entirely untouched.**
   My per-arm-threshold prediction was **refuted**: every arm cliffs at ~40 including minimap2,
   which is not on the BBMap path, and we pass minimap2 no min-intron parameter. Revised (untested)
   hypothesis: a **scoring crossover** — affine deletion cost grows with length while intron cost is
   roughly fixed, so gaps flip above the crossover with no threshold required, which also explains
   the steep rise rather than a hard zero. **Test: vary minimap2 `-O`/`-E` and re-histogram.**
   Still open: **uLTRA's sub-20 population** (132 N-ops < 20 bp, the only arm with any) (arms carry different thresholds — 40 BBMap-family, 20
   STAR/HISAT2 — so a per-arm cliff at each arm's own threshold would confirm outright; the
   sub-40 residue of 19+156 is consistent with a second path). And `wtaa`/`ysh1` remain queued,
   though deposit-vs-background is largely reframed: the driver is **per-library deletion rate**,
   and "deposit" was the label for it.

7. **(the trail that got there — kept because the refutations are the useful part)** ORIGINALLY:
   "THE BIGGEST UNRESOLVED THING IN THIS WORKSPACE: ~70 % of spliced reads
   in some DRS samples each invent a UNIQUE, HIGH-QUALITY junction.**
   *(Started as "a 46× spread in Station C junction counts"; two rounds of discriminating tests
   have moved it a long way — the noise explanation is REFUTED.)*
   Across 16 samples at matched ~1M depth: total censused 741 → 34,122 (**46×**), novel **70.7×**,
   while **`annotated` is flat at 222–270 (1.21×)** — so real biology is detected consistently and
   the entire spread is in NOVEL junctions.
   **My analysis of their table:** the gate absorbs most of it — `admit_candidate` spread is only
   **10.6×**, and the admit RATE runs **1.12 % in the four high-count samples vs 4.84 % in the
   five low-count ones (4.3× lower)**. Shape is consistent with the excess being low-evidence junk.
   **✅ DISCRIMINATOR RUN (25 samples) — NOISE IS REFUTED, it is the dilution branch:**
   | metric | HIGH-count | LOW-count |
   |---|---:|---:|
   | median `q_max` | **93.3** | 28.5 |
   | % q < 50 | **12.0** | 62.4 |
   | % repeat_flag | **12.6** | 42.6 |
   | **% support == 1** | **96.2** | 81.9 |
   The excess junctions are BETTER informed and LESS repeat-associated than the low-count
   samples' — they simply have one read each. This also kills `pool_gate.py`'s 3-nt-anchor prior
   as the driver (that predicts low q, high repeat_flag; both go the other way). **Read length is
   also refuted** by the peer: `wt_rep1` and `dst1d_rep1` have near-identical length distributions
   and differ 46×. N-op-carrying read fraction differs only 1.3× (23.4 % vs 17.9 %) — so it is
   not a splicing-rate difference either.

   **MY DECOMPOSITION of their figures — the anomaly is the singleton tail:**
   | sample | novel | singletons | multi-read core | reads/core junction |
   |---|---:|---:|---:|---:|
   | wt_rep1 | 33,853 | **32,567** | 1,286 | **11.1** |
   | dst1d_rep1 | 479 | 392 | 87 | **407.0** |
   Singleton tail **83×** apart; multi-read core only **14.8×**. **69.6 % of wt_rep1's
   junction-carrying reads produce a junction seen exactly once, vs 1.1 % in dst1d_rep1.** At
   matched depth that is reads being SPREAD, not more junctions existing.

   **🔴 UPDATE — H1 and H2 BOTH REFUTED, and H2 came out INVERTED. There are TWO phenomena in
   different samples, not one:**
   - **H1 coordinate jitter: REFUTED.** Proximity-merge collapses `wt_rep1` only **1.7×** at
     ±50 bp (33,853 → 20,093) and 1.2×10⁴ even at ±200 bp; low-count samples collapse by a
     similar 1.3×, so the behaviour does not distinguish them. **These are genuinely distinct
     junctions.**
   - **H2 intron length: REFUTED, and BACKWARDS.** The HIGH-count excess is **SHORT — 81 bp
     median, only 1.7 % above 1 kb**. The **LOW**-count samples carry the implausible lengths
     (`ski7d_rep1` core median **24,313 bp**, 72.8 % > 1 kb). **The impossible-junction class the
     10 kb guard addresses is concentrated in the LOW-count samples.**
   ⇒ **Two problems:** HIGH-count = tens of thousands of distinct, short (~81 bp), high-q,
   **81.5 % NON-canonical**, singleton junctions. LOW-count = few junctions, large fraction from
   the long spurious class (Open #1's residual).

   **MY CURRENT HYPOTHESIS for the 81 bp class — a threshold artefact, testable in one histogram.**
   The panel relabels DELETIONS as introns at a fixed length, explicitly in source:
   `multi_aligner.py:1016` *"BBMap `intronlen` is the MIN deletion length relabeled D->N"*,
   `:1018 'intronlen=40'`, `:1177` *"any reference gap >= 40 bp is encoded as an N-op"*;
   `overhang_resolver.py:85 min_intron = 40`; `--min-intronlen` default 20 on the STAR/HISAT2 path.
   **An 81 bp median sits just above a 40 bp cutoff.**
   **TEST: histogram singleton intron lengths across 20–120 bp. A CLIFF at exactly 40 (nothing
   below, pile-up above, decaying tail) ⇒ these are DELETIONS RECLASSIFIED AS INTRONS, not
   spurious splicing. Smooth through 40 ⇒ the peer's error-cluster reading stands.**
   Mechanistic note: an N-op consumes REFERENCE without consuming query, so it does not absorb
   errors *in* the read — it absorbs a reference gap, and DRS's dominant error mode is deletion
   relative to reference. **If this holds the fix is NOT in Station C** — it is the D/N boundary
   (raise `intronlen`, or require a canonical motif before a SHORT gap may be called an intron) —
   different code, different owner, and untouched by anything landed tonight.
   ⚠️ It also predicts H3 for a DIFFERENT reason: arms carry different D→N thresholds (40 / 20 /
   minimap2's own), so an arm-dependent singleton rate may be threshold difference, not error
   handling. **Per-arm intron-length histograms, each showing a cliff at its own arm's threshold,
   would settle it outright.**
   **`review` counts are NOT comparable across samples in either direction** — same column, two
   different phenomena.

   **(superseded) the three candidate mechanisms + their tests (all groupbys on existing data, all orthogonal to
   the deposit-vs-background test, so they can run in parallel):**
   - **H1 coordinate jitter beyond the ambiguity window** — the same true junction called at
     slightly different coordinates per read never coalesces. `PoolGateConfig.max_ambiguity_shift
     = 30`; `ambiguity_window` defaults `max_shift=60`. Fits everything: high q, support 1, low
     repeat_flag. **Test: proximity-cluster novel junctions, merging when BOTH boundaries fall
     within ±50 bp then ±200 bp. If 33,853 collapses to ~10³ clusters it is jitter.**
   - **H2 the impossible-junction class below threshold** — long clips placed idiosyncratically;
     long clip ⇒ high q. The ~860/sample figure counts only the >10 kb tail; the 2–10 kb residual
     (Open #1) and everything below is uncounted. **Test: intron-LENGTH distribution of singleton
     junctions; excess above the ~1 kb yeast maximum ⇒ spurious placement.**
   - **H3 arm attribution** — **test: winning aligner for reads contributing singleton vs core
     junctions.** Disproportionate uLTRA/deSALT ⇒ same aligner story in a new guise.
   H1 and H2 are not exclusive; jitter could be the downstream signature of long-clip placement.

   ⚠️ **IF H1 HOLDS, THE 46× IS SUBSTANTIALLY AN ARTEFACT OF COUNTING** — true junction diversity
   would be far closer between samples than the totals suggest, and Station C's `support >= 2`
   gate would be defeated by jitter rather than by genuine singleton biology. **Worth settling
   before any beta reviewer reads a junction count.**

   **Deposit-vs-background is still open and is the peer's test:** `wtaa_rep{1,2,3}` (same
   HHY168-AA background, DIFFERENT deposit) is the clean discriminator — high ⇒ background, low ⇒
   deposit/run effect. `ysh1_rep{2,3}` is the reciprocal. `rrp6aa_rep2` at 2,136 in the same
   background already argues against a pure background story. Note that even if "deposit" wins,
   it may be the LABEL for whatever drives the jitter (basecaller version, chemistry, read
   quality) rather than the explanation.
8. **Station B costs ~30× Station C** (751–1,062 s vs 34 s, ~1 row per read). Relevant to run-all
   capacity: if a reviewer needs only junction verdicts, C alone is ~3 % of the cost.
9. 674 still needs `h_data` raised; 668b exit-1 root cause still unconfirmed (682 owns it).

## Lesson worth carrying (three instances tonight)

**A guard that references a name not in scope fails only when the guard fires — so it is the
NEGATIVE test that catches it, not the happy path.** `os` in `station_c.py` and `logger` in
`analyze_command.py` both raised `NameError` instead of doing their job, and both were caught only
because a test exercised the failure case. Related: my resolver ceiling test passed **vacuously**
until an anti-vacuity assert was added, and again when a tandem-repeat fixture met master's period
gate. **Write the negative case, and assert the guard's precondition was actually reached.**

## Resume

```bash
# 1. 676 (a) — the deletion gate. NEEDS KEVIN'S APPROVAL OF THE FILE LIST.
ssh hoffman2 'ls -l /u/project/guillom/kevinroy/676_ubam_retire/'
```
**md5 is the WRONG test** (9 files differ from their deposit twin by exactly 27 bytes of `@PG`; the
3 wtaa are genuine shard merges). Use `samtools view -c` retired vs
`shared/raw/cDNA/intronic_pa_snp1_prp28_prp5_cdna_pcb114_2026/<sample>.bam`, and for wtaa sum over
`cdna_wt_rna15_ysh1_20260711/wtaa_rep*/PBM45482_*.bam`. `pt` already confirmed preserved
(planning/676 §9). **Show Kevin the list; delete only on his go-ahead, explicit paths, no globs.**

```bash
# 2. Reproduce the impossible-junction census (expect "3 268")
ssh hoffman2 'export PATH=~/.conda/envs/rectify/bin:$PATH; python -c "
import pysam
b=pysam.AlignmentFile(\"/u/scratch/k/kevinroy/684_p1cdna_1M/WT_BY4742_rep1/align/WT_BY4742_rep1.multialigned.bam\")
L=dict(zip(b.references,b.lengths)); past=big=0
for i,r in enumerate(b):
    if i>400000: break
    if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
    if r.reference_end and r.reference_end>L[r.reference_name]: past+=1
    if max([l for o,l in r.cigartuples or [] if o==3] or [0])>10000: big+=1
print(past,big)"'
```
Post-fix reruns should trend toward 0 past-contig; `stats['impossible_intron_truncated']` and
`Xn:i:` tags are the positive signal that the guard fired.

```bash
# 3. Peers' waves — NOT yours; recover only if their sessions died
ssh hoffman2 'qstat -u kevinroy; ls -l /u/project/guillom/kevinroy/682_sentinels/'
```
682 (DRS): 32/51 chains at last report, zero failure sentinels. 684 (cDNA): 27 samples through
consensus. Both have their own handoffs (`HANDOFF_682_drs_arm.md`, `HANDOFF_684_p1cdna_1M.md`).

**Awaiting from peers:** the remaining Station C samples (16 of ~32 in), the `q_max`/`min_anchor`
stratification that decides the 46× question, and the post-re-run survival number vs the banked
3,789/3,789 + 84/84 baseline. `684e` has LANDED (see Open #4).

## Files

- `origin/master` = `cdf4bd9`; branch `chore/vendor-desalt-chanfreau1` merged up to master
- `rectify/core/consensus/intron_sanity.py` (new) · `tests/test_intron_sanity.py` ·
  `tests/test_station_c_output_prefix.py` · `tests/test_analyze_manifest_input.py` ·
  `tests/test_align_require_aligners.py` · `tests/test_resolver_candidate_ceiling.py`
- `dev/683_G_ARBITER_EXPLORED_NOT_NEEDED.md`
- Chanfreau `planning/`: `679_*` `680_*` `681_*` `683_*` `684c_*` `689a_*` (DRS's shared caveat)
- H2: `/u/project/guillom/kevinroy/{676_ubam_retire/,682_sentinels/,684_chain.sh}` ·
  `/u/scratch/k/kevinroy/{684_p1cdna_1M/,644k_blast_radius/}`
- Sentinels: `.676b_move_ubam_rc`=**0 (done)** · `644k_blast_radius/.644k_rc`=5 (corrected)

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
