# HANDOFF — rectify: fixing fa7f315's two defects (cDNA tails XK=0; records past a contig end); nothing pushed

Session 6cfbe12b (rectify workspace), opened 2026-10-12 14:06 PDT. Thread folder: `~/work/rectify/dev/fixes_20261012/` (gitignored).
Builds on `../prod_20261009/HANDOFF.md` (§3 items 1–3, §4 Resume) and `../prod_20261009/CDNA_EQ_BUG.md`.

**State in one line:** origin/master = h2/master = `f15b523`. Push candidate `land/fixes-20261012b` = `2d0f00b` is VERIFIED and
Kevin was asked (PushNotification 19:2x PDT); no answer yet, nothing pushed. NEW (Kevin, 2026-10-12 evening): investigate the poly(A)
alignment artifact on atlas review card a1 (YNL246W, chrXIV:187,654) → `SPLICED_TAIL.md`: mechanism found at read level (tail + ONT
adapter `CT` spliced by a long N onto a genomic A₁₃`CT`; walkback blocked); scope scan DONE (H2 15143622): 63 atlas
clusters dominated (16 in the 90 % set, 2 rank-1), 6 of 23 review cards; the current pretrim pipeline makes 200–600× fewer.

## 1. Done
- §4 steps 0–1: login is gmail (queue unavailable, not needed); fetched origin + h2 (f15b523); read the four 1090 notes and
  1090z32 (notes NOT archived yet: archive after the replies go out).
- Tail rule measured on un_60S_bc19, RPL20B 20b_1 and the route smoke against Dorado XP (`TAIL_RULE.md`, H2 jobs 15139807.3,
  15139846). Decision: carried `XA` when the record carries `XK` (any value) and `XA > 0`, else dfbda14's rule.
- Past-contig origin (H2 job 15139832, `oob_origin.py`, outputs on H2 `route_smoke/tailrule/oob_origin.json`): the 875 records in
  July Path-A wtaa_rep1 `align2_fullpool` are spread over every chromosome (chrXV 360, chrIV 174, chrV 73, …, chrI 28; 1090's
  "all on chrI" came from the first 20 of a coordinate-sorted BAM), all primary, all `Xa=deSALT` (`Xn=1`, `Xc=low`). The
  deSALT arm holds 1,094 (incl. secondary); gapmm2 and minimap2 hold 0. Bases past each contig's end match the NEXT contig in @SQ
  order at (pos − contig length): identity 0.990 over 280,190 bases vs 0.261 at a random shift. deSALT aligns across contig
  boundaries of its concatenated index; the July consensus (pre-d0e3a0f) accepted the placement. Master: d0e3a0f (2026-08-12)
  `intron_sanity` soft-clips at the first N that runs past the contig end (but leaves the spurious upstream fragment on the WRONG
  contig) and only tags (`Xi`) an M/D overrun without truncating it.
- Found + proven by test: `correct`'s BAM writers decode strictly, so `--corrected-bam` aborted on any `=` record past its contig
  (fixed in b1ef47e) and STILL stops on a pre-#6 `=`-in-clip record. The latter is the landing's pinned design (cd200a9:
  `test_historical_unplaced_equals_*`, "refuse without inventing sequence"); my first attempt (N-masking in the writers) broke
  those 2 tests and was withdrawn. fa7f315's "no longer aborts" holds for TSV-only `correct` runs only. DECISION for Kevin
  (not blocking): should the writers mask such bases N too?

## 2. Verified
- See TAIL_RULE.md log (numbers, jobs, read-level examples).
- Fail-first (H2 dev trees, login node): the 5 new tests fail at fa7f315 for their own reasons (tail 4 ≠ 33; decode abort on the
  off-contig record ×3; writer refusal) and all 12 fail at f15b523; the writer test fails at 99eddd9 on `off_eq` only.
- At b1ef47e (dev tree `route_smoke/dev_fixcdna`): eq-seq 12, cdna_analyze, xo_frame, bam_writer_parallel_smoke (incl. the
  refusal controls), writer audit, hardclip, 2F routing: 59 passed / 4 skipped; earlier at the C-v1 state the broader writer +
  validation-read set ran 196 passed, the only failures being the 2 refusal controls (since restored).
- 4 tests in test_ultracode_writer_alignment_tags.py / test_ultracode_consensus_equals.py fail on the H2 LOGIN node at fa7f315
  and at my commits alike (htslib "Failed to open file ''"): environment, not code; the qsub harness ran them green before.

## 3. Open (with why)
0. **Spliced poly(A)-tail artifact (card a1)**: UPDATE 2026-10-13 00:10 PDT: Kevin asked for a full audit first
   (`../audit_20261012/`); its §6 measured the broader A17 class: current rectify cDNA 0.002–0.009 %, DRS 0.0004–0.007 %, May cDNA
   1.10 % (RPS26B 44 % vs 1 of 4,423 now), so re-deriving the atlas inputs is the fix and a rectify guard is optional.
   Reported to Kevin (SPLICED_TAIL.md). Waiting on his choice: (1) send the scan to
   the Chanfreau 1105/1070g atlas session for an atlas-level fix (re-anchor at the donor side of the N, or re-derive the cDNA inputs with
   the current pretrim pipeline); (2) the same rule in rectify's 3′-end logic, with a real review-queue card set (read-review skill). Peer note (chanfreau-lab-01, 19:4x): route-smoke inputs on H2 scratch are past the 14-day purge; re-cut with
   `dev/prod_20260928/route_smoke_prep.sh` lines 16–22 if a re-verification needs them.
1. **Kevin: push?** Options: (a) all = `2d0f00b`; (b) the CI fix alone first = `7dfc59b` (tree identical to the suite-verified
   4a21955, ancestor of 2d0f00b, so 2d0f00b stays a fast-forward afterwards); (c) wait. Why gated: D1 approved only the landing; the
   relayed "top priority" covers the remediation, not a master push.
2. **Kevin: should the BAM writers N-mask a pre-83af3d1 '=' in a clip** like correct's intake does? Today they refuse by design
   (cd200a9's pinned refusal controls), so `correct` with a BAM output on such an old BAM stops (verified, task 6b). Not blocking.
3. **Kevin: open ISSUE-094 on GitHub** (sister to #8)? Outward-facing, so it waits for his word.
4. Lab-wide cDNA re-derivation scope (beyond 1090's IPA set) and ISSUE-087 before any DRS re-pin: unchanged, Kevin's.
5. Evidence retention: H2 `$SCRATCH` (`route_smoke/verify_2d0f00b/`, `tailrule/`, `trees/2d0f00b`) is purged ~10-26; Sherlock
   scratch after 90 days untouched. Everything is reproducible from this folder's scripts.

## 4. Resume
```bash
cd ~/work/rectify && git fetch origin && git fetch h2 && git log --oneline -1 origin/master        # f15b523 unless moved
git rev-list --count 2d0f00b..origin/master                                                        # must be 0 before any push
```
- **Kevin says push all:** `git push origin 2d0f00b:refs/heads/master && git push h2 2d0f00b:refs/heads/master`; then
  `gh run list --branch master --limit 1` and compare with run 38030311161 (its 4 failures should be gone on 3.10–3.12); drop a
  one-line sha note in `~/work/UCLA/Chanfreau_Lab/.claude/inbox/` (1090 already has the full reply); append to the shared HANDOFF.
- **Kevin says CI fix alone first:** same, with `7dfc59b` in place of `2d0f00b`; CI check; then wait for his word on 2d0f00b
  (still a fast-forward from 7dfc59b; no re-verification needed unless master moves).
- **origin/master moved:** rebase `land/fixes-20261012b` onto it as a NEW branch, re-run `verify_job.sh` (qsub -v CAND=<sha>) and
  both suite jobs (`suite_job_fixes.sh` with SUITE_NAMES, `ci_fixes_job.sh` with CAND/BUNDLE/BUNDLE_REF/TAG) before asking again.
- **Kevin says wait / no:** nothing to do; the branches stay local.

## 5. Files
- Branches / worktrees (all local; nothing pushed):
  - `fix/cdna-analyze-eq-decode` (`.worktrees/fix-cdna-eq`) = fa7f315 → e67da88 → 99eddd9 → b1ef47e → 9bb1fde.
  - `land/fixes-20261012` = a86fd0f (`.worktrees/land-fixes-20261012`): FAILED suites; kept as the record. Do not push.
  - `land/fixes-20261012b` = 2d0f00b (`.worktrees/land-fixes-20261012b`): the VERIFIED push candidate.
  - `land/fixes-20261010` = e12d22f: untouched (previous session's record).
- This folder: `HANDOFF.md`, `TAIL_RULE.md` (measurement + decision), `VERIFY.md` (pre-registration + results), `tail_rule_measure.py`,
  `tail_rule_job.sh`, `summarize_tailrule.py`, `*.json` / `*.examples.tsv` (tail measurement), `oob_origin.py`, `oob_origin_job.sh`,
  `oob_origin.json`, `arm_oob_scan.py`, `arm_oob_scan_job.sh`, `verify_job.sh`, `verify_compare.py`, `xa_of_changed_rows.py`,
  `suite_job_fixes.sh` (Sherlock), `ci_fixes_job.sh` (H2 CI stack), `fixes_20261012{,b}.bundle`.
- `dev/BUGS_TO_FIX.md`: ISSUE-094 added. Shared `~/work/rectify/HANDOFF.md`: two appended entries (16:3x and 19:2x PDT).
- Chanfreau inbox: `2026-10-13T021525Z__from-rectify-6cfbe12b__fa7f315-defects-fixed-2d0f00b-verified-not-pushed.md`.
  rectify inbox: the four 1090 notes archived to `.claude/inbox/.read/`.
- Clusters: H2 `/u/scratch/k/kevinroy/prod_20260928/route_smoke/` (`tailrule/`, `verify_a86fd0f/`, `verify_2d0f00b/`, `trees/{a86fd0f,
  2d0f00b}`, `dev_fixcdna*/`, `fixtests_20261012/`), H2 `/u/scratch/k/kevinroy/ci_newdeps_20261009/tree_{a86fd0f,2d0f00b}`, Sherlock
  `/scratch/users/kevinroy/prod_20260928/trees/fixes_{a86fd0f,2d0f00b}` + logs.
- Project-space copy of my H2 scratch evidence (made by chanfreau-lab-01 at Kevin's request, 2026-10-12 20:50–20:55 PDT; originals
  untouched): `/u/project/guillom/kevinroy/scratch_keep/prod_20260928/route_smoke/` (everything, 9,953 files, 8.5 GB) and
  `scratch_keep/ci_newdeps_20261009/` (logs, job scripts, bundles; no venv/trees/tmp). Path rule: `/u/scratch/k/kevinroy/<X>` →
  `/u/project/guillom/kevinroy/scratch_keep/<X>`; see its README.md. Files written after 20:50 (the spliced-tail scan `tailrule/st_*`)
  are NOT in it: the small results are copied here (`st_atlas.summary.json`, `st_atlas.atlas_clusters.tsv`, `st_atlas.attractors.tsv`).
  H2 scratch is at 1,958 of 2,000 GB: check `myquota` before any large write.
- Spliced-tail investigation: `SPLICED_TAIL.md`, `spliced_tail_scan.py`, `spliced_tail_job.sh`, `review_pins.tsv`, `st_atlas.*`.
- Markers: cleared (the spliced-tail scan 15143622 finished; nothing running).
