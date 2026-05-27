# Handoff — 3'SS-rescue performance fixes + A549 human jobs

**Date:** 2026-05-27
**Branch:** `drs-validation-rebuild`
**HEAD at handoff:** `07df0ac` — `docs(figures): fix hardcoded #ffffff fills in dark-mode figures`
**My commits:** `fe1bb69`, `25d7a30`, `5ceb243` (see §2)
**Scope:** Fix `_rescue_3ss_truncation_body` / `_terminal_peel_rescue` stalls on
human chr5 (SMN1/SMN2 locus); monitor and iterate on Sherlock A549 jobs.

---

## 1. Session goal and outcome

**Goal:** Get four Sherlock RECTIFY jobs running without stalling on the dense
SMN1/SMN2 locus of human chr5 (SG-NEx A549):

| Job | Aligner | Status at end of session |
|-----|---------|--------------------------|
| WM2 | winnowmap2 | cancelled (not resubmitted — see §4) |
| MSP | minisplice + minimap2 | submitted as 26153569; no longer in queue |
| A549 correct | full panel | cancelled (not resubmitted) |
| nopool deSALT | deSALT, no junction pool | finished |

**Outcome:** Three performance fixes were committed and deployed to Sherlock;
smoke instrumentation revealed the root cause (read-side repeat expansions), and
a separate agent session (`91d8336`) fixed it definitively before this session
ended. The Sherlock rsync is stale by 11 commits — **the most important pending
action is to rsync HEAD to Sherlock and resubmit all three unfinished jobs.**

---

## 2. Fixes committed this session

### Fix 1 — `fe1bb69`: K=20 cap in `_terminal_peel_rescue`

**File:** `rectify/core/splice/splice_aware_5prime.py`

On junction-dense human loci, `_terminal_peel_rescue` collects up to 200+
candidate junctions. Each peel depth called `_rescue_3ss_truncation_body`
O(N × shifts × offsets × HP-ED) times, stalling 8 workers indefinitely.

Fix: sort candidates by edge distance from `align_5prime`, keep the K=20
closest. Recompute `max_edge_dist` from the capped set.

### Fix 2 — `25d7a30`: Partitioned K=25 cap in `_rescue_3ss_truncation_body`

**File:** `rectify/core/splice/splice_aware_5prime.py`

The **direct rescue path** (separate from the peel path) also admitted 50–200+
junctions. Fix: partitioned K=25 cap immediately after `_nearby_junctions`
collection.

Partition logic:
- When `_n_intervals` populated (mapPacBio leading-D path): preserve all
  N-op-matched junctions first (they can sit far from `align_5prime` and would
  be wrongly dropped by edge-distance sort); apply the K budget to the
  proximity-only remainder.
- When `_n_intervals` empty (MSP, WM2, most non-mapPacBio): simple
  sort-by-edge-distance + slice.

### Fix 3 — `5ceb243`: Offset-loop ED=0 early exit

**File:** `rectify/core/splice/splice_aware_5prime.py` (~lines 1524, 1658)

Added `if _best_local_ed == 0: break` inside the `_off` loop (both plus and
minus strand). The tiebreaker tuple `(not _in_amb, not _donor_ok, _shift_abs)`
depends only on the shift — no later offset can change the winner once ED=0.

**Safety:** Shift-loop and junction-loop breaks on ED=0 are UNSAFE — different
shifts/junctions can tie at ED=0 but differ on tiebreaker components. Only the
inner offset loop is safe to break early. A test regression (`test_plus_softclip_rescued`)
caught this during development.

This commit also added empirical mapPacBio ONT-DRS quality data to
`docs/aligners/mapPacBio.md` (SMA_GSB2394 chr5, 2026-05-25): 21.7% of mapPacBio
introns are unique-to-aligner, 97.7% of those are novel vs GENCODE → mapPacBio
dropped from human DRS panel.

---

## 3. Smoke instrumentation and key finding

To decide between option b (`_off_limit` cap) and option d (junction-density
gate), temporary instrumentation was added to `_rescue_3ss_truncation_body`
logging `3SS_SMOKE K=... hped=... max_off=... best_ed=...` per rescue call.

**The smoke ran on MSP job 26153569 (resumed from region_0357, 30,427 calls).**

| Metric | Value |
|--------|-------|
| Calls with hped > 0 | 14,780 (48%) |
| `max_off` ≥ 110 | **0** — never |
| `max_off` unique values | 0–10 only |
| `best_ed = 0` (perfect match) | 1,440 (4.7%) |
| **`best_ed = -1` (no rescue found)** | **22,678 (74.5%)** |
| hped median (nonzero) | 121 (= 1 junction × 11 shifts × 11 offsets) |
| hped max | 3,182 |

**Critical finding: `max_off` is never ≥ 110.** My earlier calculation of
30,250 HP-ED calls/read (K=25 × 11 × 110) was based on a wrong assumption.
`_off_limit = min(junction_proximity_bp, dist)` = `min(110, dist)` where `dist`
is always ≤ 10 in practice (the 5' alignment end is never deeply intronic for
these reads). The actual per-rescue max is ~3,182 calls (≈ 26 junctions × 11 × 11).

**Option b (`_off_limit` cap) was irrelevant** — the loop is already capped at
10 by `dist`.

**Root cause of the 74.5% failure rate:** reads entering rescue with long
low-period repeat-expansion soft-clips (Nanopore motor-slippage artifacts: AAG/CTT
multimers, hundreds of bp). These match nothing in the junction pool and produce
`best_ed = -1` after exhausting all K × 11 × 11 calls. Commit `91d8336` (by a
parallel Claude Opus 4.7 session) detected and short-circuits these before the
candidate loop — eliminating ~86–92% of futile human RNA004 rescue calls.

**The smoke instrumentation was NOT committed and is NOT in the current working
tree.** It was superseded by HEAD changes from the parallel session.

---

## 4. Sherlock state — action required

**Sherlock rsync is stale by 11 commits.** The Sherlock copy at
`/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/` has only
up to commit `5ceb243`. It is missing:

- `91d8336` feat(splice): detect read-side repeat expansions — **the key fix**
- `8fd0bcd` docs(handoff): human empirical error-model task
- `c1d0e83` docs(handoff): rRNA-exclusion fix
- `befbc3a` + `c0702c4` perf profiling / task #12
- `ae3ec41` calibration: variant masking
- `c52966c` docs: mapPacBio param fix + HP-ED gaming
- `1190129`, `fbde573`, `08f2696`, `d89b10b`, `07df0ac` — branding, figures, winnowmap2 fix

**To deploy and resubmit:**

```bash
# 1. Rsync the full rectify/core to Sherlock (safest: sync the whole package)
rsync -avz /Users/kevinroy/work/rectify/rectify/ \
  sherlock:/oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/

# 2. Verify repeat_expansion landed
ssh sherlock 'grep -c "is_repeat_expansion" \
  /oak/stanford/groups/larsms/Users/kevinroy/software/rectify/rectify/core/splice/splice_aware_5prime.py'
# Should be > 0

# 3. Resubmit MSP (resume from region_0357 checkpoint)
ssh sherlock 'cd /scratch/users/kevinroy/rectify_human_validation/sgnex_a549/correct && sbatch run_minisplice.sh'

# 4. Resubmit WM2
ssh sherlock 'cd /scratch/users/kevinroy/rectify_human_validation/sgnex_a549/correct && sbatch run_winnowmap2.sh'

# 5. Resubmit A549 correct (full multi-aligner) — check if the script exists
ssh sherlock 'ls /scratch/users/kevinroy/rectify_human_validation/sgnex_a549/correct/*.sh'
```

After submitting, watch for MSP checkpoint progress past region_0357. With the
repeat-expansion gate, `91d8336` should eliminate most of the futile 74.5% calls
and let region_0357 complete.

**nopool deSALT** finished (job 26131887 is no longer in the queue). Its output
is at `/scratch/users/kevinroy/rectify_human_validation/sgnex_a549/correct/`
(check `nopool_desalt/corrected_reads.tsv` or equivalent).

---

## 5. Uncommitted working tree on M1

`git status` shows four modified files from **other sessions** (not this one):

- `AGENT_FIXES.md` — new entries from parallel perf-profiling session
- `dev/PERF_AUDIT.md` — 178-line expansion (task #12 findings)
- `dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md` — §8 addendum
  (periodicity test results, DROP verdict)
- `dev/profile_correct_reads.py` — profiler improvements

These are not this session's work. Stage and commit them as a separate docs/perf
commit if you confirm their content is accurate, or ask Kevin whether to commit
them first.

---

## 6. Broader context (where this fits)

This session's work is part of the human DRS validation pipeline documented in
`dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`. The remaining goals
from that doc are:

1. **Re-align A549 with `--max-intron 500000`** — the current BAMs used `-G 5000`
   (yeast ceiling). Affects Cat3/5/6/7/9 junction categories.
2. **`correct --annotation` tractability** — validate that `91d8336` + pool
   anchor floor (`ed3df74`) makes `--annotation` usable on human (was 475×
   slower / 2 h timeout; without `--annotation` running correctly now).
3. **Build the committed validation set** — 72 reads × 9 categories, sourced from
   SG-NEx A549 (public), no patient reads.
4. **Human empirical error model** — high priority; see §5 of the prior handoff.

---

## 7. AGENT_FIXES.md entries from this session

Added entries (committed in `25d7a30`/`5ceb243`):

- `[2026-05-26] PERF: _rescue_3ss_truncation_body O(N) junction iteration — FIXED`
- `[2026-05-26] PERF: offset-loop ED=0 early exit in _rescue_3ss_truncation_body — FIXED`
- `[2026-05-26] PERF: junction_pool.pkl inflated by deSALT nosec artifacts — DIAGNOSED` (committed `001963c`)

Update the deSALT pool entry status: the nopool job (26131887) completed, so the
workaround (run without `--junction-pool-cache`) is confirmed working. The pool
itself still needs rebuild from non-deSALT aligners.
