# HANDOFF — DRS arm (unit 668). Brief for a dedicated agent.

**Written 2026-08-11 ~22:30 PDT by session `rectify-realigner`.**
You own the **DRS cohort end to end**: get the gate green, then fan out all 9 samples.

**Kevin's directive (2026-08-11): "The DRS arm should proceed unhinged, while the cDNA
datasets must await this fix."** You are the arm that moves. Do not touch cDNA
(662 / 673 / 674 / 678) — those are frozen behind `HANDOFF_679_cdna_trim.md`. Do not touch the
676 storage move (job 14293343) — that stays with `rectify-realigner`.

Dataset triage context: `~/work/UCLA/Chanfreau_Lab/planning/680_dataset_status_drs_go_cdna_hold.md`.

---

## 1. Why DRS is not blocked

The 679 defect (adapter/UMI never stripped from cDNA consensus molecules) lives entirely in
`rectify/core/cdna/`. **DRS never enters that code** — no PCB114 SSP, no UMI, no clustering, no
consensus molecule, no `pretrim_consensus`. DRS goes raw → trim → align → resolver → consensus →
prescan. Nothing in the 679 chain touches it, so there is nothing for you to wait on.

What DRS *does* share with cDNA is the **resolver throughput** problem — but for DRS that is a
cost issue with a landed lever, not a correctness issue. See §3.

## 2. The cohort

`/u/scratch/k/kevinroy/668_drs/manifest_668.tsv` — 9 samples, all from the locked deposit
`shared/raw/DRS/PRJNA1229592_cpa_depletion_DRS_2025/`:

`wt_rep1..3` · `rna15_rep1..3` · `ysh1_rep1..3`

This is **Kevin's P1 cohort** (WT / Rna15-AA / Ysh1-AA) — the multi-strain half of the rbrowse
multi-strain × multi-modality test. Its cDNA counterpart is frozen, so DRS is currently the only
modality that can deliver P1.

## 3. 🔴 The blocking issue: the 668b smoke FAILED, and the sizing assumption was wrong

**Job 14287090.1, `wt_rep1`: `exit_status 1`, `failed 0`, wallclock 4 h 08 m, maxvmem 828 M.**

Neither cause people assume:
- **NOT walltime** — the script requests `h_rt=24:00:00`.
- **NOT memory** — 828 M peak against `h_data=3G × 8 slots` = 24 G.

**What the evidence actually says.** The stdout ends at `=== [S2b] align w/ overhang_resolver
substitution ===`. The next line the script would print is `ALIGN_SECONDS=… rc=…`, and it never
appears — so the shell died *inside* the align call. Critically, the script has
`trap _finish EXIT` which writes `.rc668b` and echoes `=== [wt_rep1] EXIT rc=… ===`; **neither
happened**, so the shell did not exit normally — it was killed hard. There is also no
`FATAL: align+resolver` line, which rules out the script's own error branch. Candidates: a
node-level event, or the align child dying in a way that took the shell with it. **Diagnose
before fanning out — there is no green DRS gate until this is understood.**

**The sizing datum nobody has stated.** At job end the resolver BAM was 462,880,768 B and still
growing, against a 2,468,849,530 B minimap2 arm — i.e. roughly **19 % of the arm in 4 h 08 m**
(BAM-size is a rough progress proxy; compression varies, treat as ±). That projects to
**~22 h per sample single-threaded**, against a 24 h `h_rt`. So this job was going to fail on
walltime even if nothing had killed it. The cohort was planned against a "DRS = 1,600 reads/s
class" figure; the real rate here is ~2 orders below that. **Do not fan out until you have
re-timed.**

## 4. The lever, and it is already landed

**`origin/master` = `277c708`** — the numba twin for `hp_edit_distance_bounded`, **289×** on
resolver-shaped inputs, bit-identical (4,000 randomized pairs exact, every-row pruning
preserved). Full suite **2,262 passed / 41 skipped / 1 xfailed, SUITE_RC=0**. Verified tonight.

**numba 0.66.0 IS installed** in the H2 env `~/.conda/envs/rectify` (also `edlib` — both
confirmed importable). The kernel is **opt-in**: `export RECTIFY_HP_ED_NUMBA=1`.

⚠️ **The smoke did NOT have it.** `668b_resolver_smoke.sh` pins
`CLONE=/u/scratch/k/kevinroy/670_tree/rectify`, which is at **`e7499ef`** — before `11ed647`
(649 complexity gate default-ON, threaded-output determinism) *and* before `277c708` (numba).
A tree at `origin/master` is therefore both **much faster and a different configuration**; the
partial e7499ef results do not carry over. Re-validate acceptance, don't assume.

`/u/scratch/k/kevinroy/tree_11ed647` exists but is one commit behind master.

## 5. Your first moves

1. **Update a tree to `origin/master` (277c708)** and repoint `CLONE` in the script, or update
   `670_tree` in place. The script's import assert (`assert p.startswith('$CLONE')`) protects
   you from the 596 editable-install trap — keep it, and keep the deSALT md5 assert
   (`1f255a9b2ca783669f1520cfd5ffa01a`).
   The comment in the script says the live tree is left alone because "the 662 array still
   imports from it" — **662 is now frozen (cDNA)**, so that constraint has lapsed; the live-tree
   swap to master is unblocked if you want it.
2. **`export RECTIFY_HP_ED_NUMBA=1`** in the job script. Keep `RECTIFY_SKIP_REGIONS=yeast-rdna`
   (already there, and correct — it is the same pathological-contig class as the reporter
   constructs).
3. 🔴 **Delete the partial resolver outputs before rerunning, by explicit named path:**
   `/u/scratch/k/kevinroy/668_drs/wt_rep1/align/wt_rep1.overhang_resolver.bam` (462 MB,
   TRUNCATED — no BGZF EOF) and any `wt_rep1.resolver_delta.bam`. The align step reuses
   existing BAMs via `--trust-existing-bams`; if the substitution path also honours an existing
   resolver BAM you would silently consume a truncated file. **No globs — name the paths.**
   Do NOT delete the three aligner BAMs or `multialigned.bam`.
4. **Rerun is cheap.** All of S1 (trim) and the three aligner arms are complete on disk:
   `wt_rep1.minimap2.bam` (2.47 G), `.uLTRA.bam` (2.32 G), `.deSALT.bam` (2.44 G), plus the
   `.rnsorted.bam` trio and `multialigned.bam` (2.37 G). Only resolver + consensus + prescan
   need to run.
5. **Re-time on a bounded slice first** before committing 9 samples. Report reads/s.

## 6. Acceptance gate (the script already enforces most of it — read the block, don't trust exit 0)

The smoke's own checks, which you must see *pass in the log*:

- `overhang_resolver.bam` exists **and** `overhang_resolver.stats.json` exists (the script only
  `WARN`s on a missing stats.json — treat a missing one as a **failure**, it is the beta-ledger
  record).
- **The 670 invariant:** `consensus_records == primary == distinct_qnames`. This is the guard
  for the DRS group-split bug (`f652331`, RN-keyed merges need RN-sorted inputs) that was
  silently emitting ~2.2 primary records per read on **all DRS**. `e7499ef` already contains
  that fix, so this should pass — but it is the single most important line in the log.
- `resolver_delta.bam` (XB-tagged records only) extracted, `NDELTA` and `NXB` reported.
- `junction_pool.pkl` written by the S3b prescan.
- Sentinels `.s2b` / `.s3b` / `.rc668b` / `.done668b`.

**Read the acceptance block. Do not flip anything to green on exit code alone** — tonight's
`644k` unit wrote `rc=0` while its log contained a `Segmentation fault`, which is exactly how
this workspace gets fooled.

## 7. Fan-out

Only after `wt_rep1` passes acceptance **and** you have a measured per-sample wall:

```bash
cd /u/scratch/k/kevinroy/668_drs && qsub -t 2-9 -tc 4 668b_resolver_smoke.sh
```

Task IDs map to `manifest_668.tsv` line numbers. Keep `-tc 4`. **du-gate every wave.**

⚠️ **H2 scratch is at 1,977.6 / 2,000 GB — 22 GB of headroom.** A 313 GB move
(`move676par`, job 14293343) is freeing space incrementally as each file lands; check
`myquota | grep scratch` before each wave. **Do not stage new data on scratch until that
sentinel (`/u/project/guillom/kevinroy/.676b_move_ubam_rc`) reads 0.** Never `rm` a BAM to make
room — retire it to another path. Sherlock scratch is Kevin-authorized overflow; SCG is
available and shares the Oak mount.

## 8. Durability (non-negotiable)

Claim your own `NN` atomically:
`~/work/bin/nn_claim.sh ~/work/UCLA/Chanfreau_Lab/planning <owner> <slug>` — do **not** write
into 668/679/680, which are other sessions' records. Write your plan to `NN_*.md` **before**
any compute, then append a checkpoint after each sub-step (tree updated → partial deleted →
re-timed → acceptance read → wave N submitted). Drop `.claude/.handoff-needed` in the rectify
workspace the moment you `qsub`, and refresh `HANDOFF.md` at every stopping point. Sentinels the
job writes itself are the durable signal — a background watcher is convenience only.

## 9. Traps

- **`XB` means two different things in this cohort.** In the resolver arm it marks rewritten
  records (`resolver_delta.bam`); on cDNA consensus molecules `XB:Z:n/m` is the strand split.
  Don't cross-read them.
- **`rectify align` silently DROPS a requested junction aligner and exits 0 if its binary is
  absent from PATH.** The script guards this with the deSALT md5 assert — keep it, and assert
  all three arm BAMs exist. A "3-aligner run" with 2 arms passes an exit-code-only gate.
- **`run_overhang_resolver(threads=N)` accepts `threads` but never uses it** — the body is a
  single-threaded stream. The script passes `-t 8` in good faith; that buys nothing for the
  resolver leg. Plan capacity single-threaded (this is a known backlog item, not your fix).
- Surgical staging: `git add <explicit paths>`, never `git add -A` — Kevin routinely has WIP in
  `rectify/core/bam/*`, `rectify/core/correct/*`, `scripts/validation_data/`.

## 10. Report back

Drop results in `~/work/rectify/.claude/inbox/` addressed to `rectify-realigner`. Kevin wants
the finding in **text, in the same message you first mention it** — the number, the direction,
the n, and the one caveat that could overturn it. Never let figure/report polish gate a result.

**The two numbers he is waiting on: the diagnosed cause of the 668b exit-1, and the re-timed
per-sample resolver wall with numba on.**
