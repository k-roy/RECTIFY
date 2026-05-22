# HANDOFF 2026-05-20: set2 minimap2-only salvage restart

**For:** the next agent picking up set2 production after this session.
**Mission:** get COMPLETED `corrected_reads.tsv` + `corrected.bam` for all 12 set2 samples by restarting with **single-aligner (minimap2-only) `rectify correct`** that bypasses the in-flight axis 2 htslib state corruption bug. Then run analyze and produce final per-condition bedgraphs for Kevin's downstream analysis.

**Don't re-investigate** the bugs we discovered today. They are documented; the goal is execution.

---

## 1. Current state (start of your session)

### 1a. Cascade in flight, mostly failing

A 99-job set2 cascade with 5-aligner consensus correct (`--aligner-bams` for each of mapPacBio/minimap2/gapmm2/uLTRA/deSALT) is RUNNING on Sherlock. It's hitting **axis 2: scale-induced htslib state corruption** documented in `AGENT_FIXES.md` entries 2 and 6 and root-cause-briefed in `dev/specs/briefings/axis2_htslib_state_corruption_briefing.md`. Of ~530 expected per-chunk correct tasks, **roughly 0 are producing real output** (the few "COMPLETED" tasks are empty-BAM-guard fallbacks for deSALT SIGSEGV chunks — they exit 0 with a header-only TSV in seconds, no real correction). Most live tasks are silent-hung in corrupted memory state at ~2-3 hours elapsed.

Job IDs at `/oak/stanford/groups/larsms/Users/kevinroy/set2_postfilter_ids` on Sherlock.

A separate **wt_rep2 5-aligner smoke test** (chain IDs at `/oak/stanford/groups/larsms/Users/kevinroy/wt_rep2_smoke_ids`) is also in flight; it WILL hit the same axis 2 bug in its correct stage. Treat it the same as the 11-sample cascade — cancel and restart minimap2-only.

### 1b. What IS already done and trustable

- **mapPacBio QNAME bug fixes** are in `rectify` HEAD (commits `382fcc7` + `e8c8070`). No re-alignment needed for set2's existing per-chunk BAMs (we have minimap2/gapmm2/uLTRA/deSALT chunks from earlier aligned-against-clean-FASTQ runs; the mapPacBio chunks are tainted but we won't use mapPacBio in the salvage).
- **Empty-SEQ filter** has been applied to all 11 samples' merged mapPacBio BAMs (366k empty-SEQ reads removed per sample). This filter is NOT needed for the minimap2-only salvage — the merged BAMs only matter for `--aligner-bams` cross-aligner lookups, which the salvage drops.
- **Prescan PKLs** (junction_pool.pkl + rescue_scan.pkl) have been built for all 11 samples + wt_rep2 from filtered merged BAMs. **You CAN reuse these for minimap2-only correct** — the prescan PKLs are universal across aligners; they're just junction-pool caches.
- **Per-chunk minimap2 BAMs** at `/oak/.../<sample>/chunks/aligner_chunks/minimap2/chunk_*/` are clean (minimap2 doesn't go through the broken QNAME sanitizer; its BAMs were never tainted).
- **Per-chunk minimap2 chunk indexes** (`*.minimap2.bam.bai`) exist alongside the BAMs.

### 1c. RECTIFY repo state

- HEAD (M1, GitHub, H2, Sherlock all): `e0c9dc3 docs(briefing): axis 2 htslib state corruption — bug debrief for dedicated agent`
- **M1 has one uncommitted WT change** in `rectify/core/analyze/bedgraph.py` (the 1-bp left-shift fix). YOU NEED TO COMMIT + PUSH this before running `rectify analyze` at the end. See §5 below.

---

## 2. Pre-flight checks (5 min)

Don't skip these. Cluster state can shift between sessions; verify before acting.

```bash
# 1. ControlMaster up?
ssh -O check sherlock

# 2. Repo in sync across clusters? (per CLAUDE.md §RECTIFY essentials)
cd /Users/kevinroy/work/rectify && git log --oneline -1
git ls-remote origin drs-validation-rebuild | awk '{print substr($1,1,7)}'
ssh hoffman2 'bash -lc "cd /u/home/k/kevinroy/software/rectify && git log --oneline -1"'
ssh sherlock 'bash --norc --noprofile -c "cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify && git log --oneline -1"'
# All four should be e0c9dc3 (or the same later HEAD if newer work has landed).

# 3. Read AGENT_FIXES.md entries 1, 2, 6 before touching anything pipeline-y.
cat /Users/kevinroy/work/rectify/AGENT_FIXES.md | less

# 4. Confirm the in-flight cascade is still running (so you can cancel it cleanly).
ssh sherlock 'bash --norc --noprofile -c "
export PS1=; source ~/.bashrc
squeue -u kevinroy -h -o \"%T\" | sort | uniq -c
"'
```

---

## 3. The salvage plan in five steps

### Step 1: Cancel the failing cascade

```bash
ssh sherlock 'bash --norc --noprofile -c "
export PS1=; source ~/.bashrc
# Cancel everything _trimmed_correct_* still in flight from the 11-sample post-filter cascade
squeue -u kevinroy -h -o %i,%j | grep _trimmed_correct_ | cut -d, -f1 | xargs -r scancel
# Also cancel wt_rep2 chain
source /oak/stanford/groups/larsms/Users/kevinroy/wt_rep2_smoke_ids
scancel \$CORRECT_mapPacBio \$CORRECT_minimap2 \$CORRECT_gapmm2 \$CORRECT_uLTRA \$CORRECT_deSALT \$CM \$FM 2>&1
# Confirm
sleep 5
squeue -u kevinroy -h -o \"%T\" | sort | uniq -c
"'
```

You're aiming for the only set2 jobs left in queue to be the prescan/merge/align that already COMPLETED. Anything `_trimmed_correct_*` should be gone.

### Step 2: Patch the minimap2 correct scripts to drop --aligner-bams

The salvage path is single-aligner: `rectify correct` reads ONLY its own input BAM, no cross-aligner consensus lookups. Removes both axis 1 and axis 2 trigger paths.

For each of 12 samples (including wt_rep2), edit `<sample>/chunks/run_array_correct_minimap2.sh` to:
- Remove all 5 `--aligner-bams` flag lines
- Keep everything else (constraint, time, mem, threads, --streaming, etc. — they're fine for single-aligner work)

Concrete patch via Python (more reliable than sed across line continuations):

```bash
cat <<'PYPATCH' | ssh sherlock 'cat > /tmp/minimap2_only_patch.py && python3 /tmp/minimap2_only_patch.py'
import glob, re
base = "/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery"
scripts = glob.glob(f"{base}/*/chunks/run_array_correct_minimap2.sh")
patched = 0
for s in scripts:
    with open(s) as f:
        content = f.read()
    # Drop all --aligner-bams flag lines (they have a line continuation '\' at end)
    new = re.sub(r'\s*--aligner-bams\s+"[^"]+"\s*\\\n', '\n', content)
    # Make sure no stray continuation lines remain
    if new != content:
        with open(s, "w") as f:
            f.write(new)
        patched += 1
        print(f"  patched: {s}")
print(f"\nTotal patched: {patched}/{len(scripts)}")
PYPATCH
```

Verify on one script:

```bash
ssh sherlock 'grep -E "aligner-bams|rectify correct|--threads|--streaming|-o " /oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_rep2/chunks/run_array_correct_minimap2.sh | head -20'
```

You want NO `--aligner-bams` lines and the rest of the rectify correct invocation intact.

### Step 3: Submit minimap2-only correct for all 12 samples

The merged BAMs and prescan PKLs are already in place. We don't need to re-run filter/prescan/merge — the salvage uses just the minimap2 chunk BAMs as input + the prescan PKLs.

```bash
cat <<'ORCH' | ssh sherlock 'cat > /oak/stanford/groups/larsms/Users/kevinroy/launch_set2_minimap2_only.sh && chmod +x /oak/stanford/groups/larsms/Users/kevinroy/launch_set2_minimap2_only.sh && echo wrote'
#!/bin/bash
export PS1=""
set -e
source ~/.bashrc
set -eu

SAMPLES="wt_rep1 wt_rep2 wt_rep3 rna15_rep1 rna15_rep2 rna15_rep3 ysh1_rep1 ysh1_rep2 ysh1_rep3 wt_tfiiib_rep1 wt_tfiiib_rep2 wt_tfiiib_rep3"
BASE_ROOT=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery
IDS_FILE=/oak/stanford/groups/larsms/Users/kevinroy/set2_minimap2_only_ids
> $IDS_FILE

for SAMPLE in $SAMPLES; do
    BASE=$BASE_ROOT/$SAMPLE/chunks
    [ ! -d "$BASE" ] && continue

    # NOTE: chunk_merge for minimap2-only doesn't need the multi-aligner consensus path.
    # run_array_chunk_merge.sh as generated expects 5 aligners' TSVs. Either:
    #   (a) Skip chunk_merge entirely and concatenate per-chunk minimap2 TSVs directly
    #   (b) Run chunk_merge anyway — it'll work but produce trivial single-aligner output
    # Recommendation: (a) — simpler, faster.
    CORRECT=$(sbatch --parsable $BASE/run_array_correct_minimap2.sh)
    # Skip CM and FM as generated; concat after CORRECT completes.
    echo "$SAMPLE correct_minimap2=$CORRECT" | tee -a $IDS_FILE
done
ORCH
```

Then run it:
```bash
ssh sherlock 'bash --norc --noprofile -c "
export PS1=; source ~/.bashrc
bash /oak/stanford/groups/larsms/Users/kevinroy/launch_set2_minimap2_only.sh
"'
```

12 array jobs submitted (1 array per sample × N chunks per array = ~120 total chunk tasks).

### Step 4: Monitor + concatenate per-chunk TSVs (after all correct arrays COMPLETED)

Poll every 5 min. Single-aligner correct should be FAST — no cross-aligner fetch, no axis 2 bug. Expected ~30-45 min per chunk (vs ~80-90 min for 5-aligner). With 8 workers, ~5-10 min wall per chunk. Total cascade wall: ~15-30 min for the small samples, up to ~60 min for wt_tfiiib_rep1 (24 chunks).

Poll command:
```bash
ssh sherlock 'bash --norc --noprofile -c "
export PS1=; source ~/.bashrc
sacct -u kevinroy -S 2026-05-20T<RESTART_TIME> -E now --format=State -X -P -n | grep -v ^\$ | sort | uniq -c
"'
```

For each sample, when all chunks COMPLETED, concatenate the TSVs:

```bash
cat <<'CONCAT' | ssh sherlock 'cat > /oak/stanford/groups/larsms/Users/kevinroy/concat_minimap2_tsvs.sh && chmod +x /oak/stanford/groups/larsms/Users/kevinroy/concat_minimap2_tsvs.sh && echo wrote'
#!/bin/bash
export PS1=""
source ~/.bashrc
set -eu

BASE_ROOT=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery
SAMPLES="wt_rep1 wt_rep2 wt_rep3 rna15_rep1 rna15_rep2 rna15_rep3 ysh1_rep1 ysh1_rep2 ysh1_rep3 wt_tfiiib_rep1 wt_tfiiib_rep2 wt_tfiiib_rep3"

for SAMPLE in $SAMPLES; do
    BASE=$BASE_ROOT/$SAMPLE/chunks
    OUT_DIR=$BASE_ROOT/$SAMPLE
    [ ! -d "$BASE/aligner_chunks/minimap2" ] && continue

    OUT_TSV=$OUT_DIR/${SAMPLE}_trimmed.corrected_reads.tsv
    # Header from chunk_000, then bodies (no header) from all chunks in order
    FIRST=1
    for CHUNK_DIR in $(ls -d $BASE/aligner_chunks/minimap2/chunk_*/ | sort); do
        TSV=$CHUNK_DIR/corrected_reads.tsv
        [ ! -f "$TSV" ] && { echo "MISSING $TSV — abort $SAMPLE"; continue 2; }
        if [ $FIRST -eq 1 ]; then
            cp $TSV $OUT_TSV
            FIRST=0
        else
            tail -n +2 $TSV >> $OUT_TSV
        fi
    done
    echo "$SAMPLE -> $OUT_TSV ($(wc -l < $OUT_TSV) lines)"
done
CONCAT
ssh sherlock 'bash /oak/stanford/groups/larsms/Users/kevinroy/concat_minimap2_tsvs.sh'
```

This produces `<sample>_trimmed.corrected_reads.tsv` per sample.

### Step 5: Commit the bedgraph fix, then run analyze

**5a. Commit the pending bedgraph 1-bp shift fix on M1.** This is critical because `rectify analyze` writes per-condition bedgraphs and the bug is in committed HEAD. M1 has the fix in WT.

```bash
cd /Users/kevinroy/work/rectify
git diff rectify/core/analyze/bedgraph.py   # confirm it's the start = int(pos), end = int(pos) + 1 change
git add rectify/core/analyze/bedgraph.py
git commit -m "fix(analyze/bedgraph): 0-based half-open BED convention (start=pos, end=pos+1)

The committed code shifted bedgraph positions 1 bp left of the
corrected_3prime value (start = int(pos) - 1, end = int(pos)),
producing systematic 1-bp shift in every per-condition bedgraph
emitted by rectify analyze across DRS, QuantSeq REV, and PCR-cDNA
workflows. corrected_3prime is 0-based-inclusive (from pysam
reference_end - 1 or reference_start, both 0-based); the correct BED
half-open representation of a single base at 0-based pos is [pos, pos+1).

Discovered in M1 WT triage 2026-05-20 (already-fixed in WT, never
committed).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
git push origin drs-validation-rebuild
```

**5b. Pull on Sherlock + H2** so they have the fix before analyze:

```bash
ssh sherlock 'bash --norc --noprofile -c "
export PS1=; source ~/.bashrc
cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify
git fetch origin drs-validation-rebuild
git pull --ff-only origin drs-validation-rebuild
git log --oneline -2
"'
ssh hoffman2 'bash -lc "
cd /u/home/k/kevinroy/software/rectify
git fetch origin drs-validation-rebuild
git pull --ff-only origin drs-validation-rebuild
git log --oneline -2
"'
```

**5c. Run analyze.** The Sherlock `run_analyze.sh` script expects per-sample directories with corrected output; you may need to update its manifest to point at the new `<sample>_trimmed.corrected_reads.tsv` files (single-aligner output), not the original per-chunk-per-aligner structure.

Inspect the manifest first:
```bash
ssh sherlock 'cat /oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/dry_run_20260518/manifest.tsv'
```

The manifest format is `sample_id\tpath\tcondition`. Update `path` for each sample to point at the new single-aligner TSV directory, or write a fresh manifest. Then submit run_analyze.sh:

```bash
ssh sherlock 'bash --norc --noprofile -c "
export PS1=; source ~/.bashrc
BASE=/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery
sbatch \$BASE/run_analyze.sh
"'
```

(You may need to edit run_analyze.sh's `--manifest` argument to point at your updated manifest, or generate a fresh one.)

---

## 4. Monitoring during step 3 (correct stage)

5-min poll cadence is plenty — single-aligner correct is much faster than 5-aligner. Watch for:

- **Failed-132 (SIGILL)** — should not happen with the existing constraint (`CPU_GEN:SIE|GEN|BGM|SKX`), but if it does, the node landed on AMD Rome and the constraint needs the `7543` SKU refinement (see `feedback_mkl_amd_sigill.md` memory, or `methods/aligner_environments.md`).
- **Timeout cascade** — if any task TIMES OUT at the 3:30 wall, the cluster nodes are unusually slow. Resubmit the array (rectify correct is idempotent — skips already-completed chunks).
- **free() errors / silent hang** — should NOT happen with minimap2-only since we dropped --aligner-bams. If you see them, axis 2 isn't fully isolated to the merged-BAM lookups and the dedicated agent's investigation (`dev/specs/briefings/axis2_htslib_state_corruption_briefing.md`) needs to extend.

Poll command:
```bash
ssh sherlock 'bash --norc --noprofile -c "
export PS1=; source ~/.bashrc
sacct -u kevinroy -S 2026-05-20T<RESTART_TIME> -E now --format=State -X -P -n \
  | grep -v ^\$ | sort | uniq -c
"'
```

Check for free() in logs:
```bash
ssh sherlock 'grep -l "free()" /oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/*/chunks/logs/*.err 2>/dev/null | wc -l'
```

Should be 0 after the minimap2-only restart. If non-zero, surface to Kevin immediately.

---

## 5. Definition of done

You're done when:

1. All 12 samples have `<sample>_trimmed.corrected_reads.tsv` files at the sample root, with line counts roughly matching the per-chunk read counts summed across that sample's chunks (verify a few — wt_rep2 should have ~4.2M lines for ~4.2M reads).
2. `rectify analyze` has completed and per-condition bedgraphs exist (the `<condition>_plus.bedgraph` / `<condition>_minus.bedgraph` files; conditions per the doc's §1 sample table).
3. Per-condition bedgraphs have correct 0-based half-open BED ranges (start = pos, end = pos+1) — verified by spot-check against a known 3' end position in a corrected_reads.tsv row.
4. No `free()` errors in any log.
5. The Session 5 handoff section in `project_status_markdowns/DRS_CPA_PROJECT_STATUS.md` is updated with the final-output paths and a "completed" status note.

Optional but nice:
- Delete the `*.pre_empty_seq_filter` backup BAMs (after the analyze output is validated; that's ~20 GB freed). They're at `<sample>/chunks/{merged_bams,aligner_chunks/mapPacBio/chunk_*}/*.pre_empty_seq_filter`.

---

## 6. If things go badly

### 6a. minimap2-only correct ALSO hits axis 2

Highly unlikely (each correct task only reads its own chunk BAM, no merged-BAM lookups). But if it happens:
- Reduce `--threads` to 1 or 2 in the wrapper scripts and retry. May expose a worker-fan-out version of the bug independent of --aligner-bams.
- Surface this to Kevin as new evidence for the axis 2 briefing (extend the scope from "merged-BAM lookups" to "any pysam fan-out at scale").

### 6b. minimap2 chunk BAMs themselves are tainted

Pysam-scan one of them for empty-SEQ:
```bash
ssh sherlock 'python -c "
import pysam
with pysam.AlignmentFile(\"/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/wt_rep1/chunks/aligner_chunks/minimap2/chunk_000/wt_rep1_trimmed_chunk_000_of_009.minimap2.bam\", \"rb\") as b:
    e = t = 0
    for r in b:
        t += 1
        if r.query_sequence is None: e += 1
    print(f\"empty/total: {e}/{t}\")"'
```

Expected: 0 empty. minimap2 doesn't go through the broken mapPacBio QNAME sanitizer. If for some reason it shows empty-SEQ reads, filter them with the recipe from AGENT_FIXES entry 1's "Faster recovery" section (adapted to minimap2 chunk BAMs — though those aren't coordinate-sorted; you'd need to sort first, or filter without re-indexing if rectify correct iterates them).

### 6c. Sherlock owners-partition preemption is too aggressive

Symptom: tasks bouncing PENDING ↔ RUNNING ↔ PREEMPTED, never making progress. If preemption is killing tasks faster than they can finish:
- Switch the partition to `larsms` only (drop `,owners`). Slower throughput but no preemption.
- Or rerun on a different cluster (H2 also has the data layout if needed; mex67 work is already there).

### 6d. The pending bedgraph fix doesn't cleanly commit

The diff is documented in §1c; if M1's WT has additional unrelated changes in that file, surgical-staging the relevant hunks with `git add -p` is the right move. The two-line change is around line 99-105 of `rectify/core/analyze/bedgraph.py`.

---

## 7. File index

### On Sherlock

- `/oak/stanford/groups/larsms/Users/kevinroy/projects/TRT/processed_data/rectify/v3_20260429/set2_cpa_machinery/<SAMPLE>/chunks/` — per-sample data root
  - `aligner_chunks/minimap2/chunk_NNN/*.minimap2.bam` — clean per-chunk minimap2 BAMs (the salvage input)
  - `merged_bams/<sample>_trimmed.minimap2.bam` — NOT used in minimap2-only salvage (only mattered for `--aligner-bams`)
  - `junction_pool.pkl` + `rescue_scan.pkl` — prescan outputs, reusable
  - `run_array_correct_minimap2.sh` — script to patch and resubmit
- `/oak/stanford/groups/larsms/Users/kevinroy/set2_postfilter_ids` — IDs from the failing 11-sample cascade (cancel these)
- `/oak/stanford/groups/larsms/Users/kevinroy/wt_rep2_smoke_ids` — IDs from wt_rep2 smoke (cancel these)
- `/oak/stanford/groups/larsms/Users/kevinroy/set2_minimap2_only_ids` — IDs from your new salvage runs (you create this)

### In the rectify repo (this checkout)

- `AGENT_FIXES.md` — entries 1, 2, 6 are required reading
- `dev/specs/briefings/axis2_htslib_state_corruption_briefing.md` — the root-cause briefing for the bug you're avoiding
- `docs/handoffs/HANDOFF_2026-05-20_set2_minimap2_only_salvage.md` — this file
- `rectify/core/analyze/bedgraph.py` — has uncommitted M1 WT fix; commit it before running analyze
- `project_status_markdowns/DRS_CPA_PROJECT_STATUS.md` — Session 5 handoff at top has the operational context

---

## 8. Why minimap2-only is acceptable for this dataset

Per Kevin's earlier analysis in the session conversation, the multi-aligner consensus in `rectify correct` mainly improves:
- Splice junction boundary accuracy (matters less in yeast which has few introns)
- Chimeric/multi-mapper disambiguation (handled OK by minimap2 alone for ONT direct RNA)
- Cases where one aligner truncates a read (rare in clean DRS data)

The 3' end calling itself — poly-A walkback, A-tract ambiguity, indel correction at the tail — operates on the primary alignment regardless of aligner. minimap2's primary alignments are reliable for these operations.

Kevin's verbatim assessment from the session: "3' end accuracy will be ~95-99% as good with minimap2 alone... If your downstream is CPA-site clustering + count-based comparison across conditions, the loss is negligible — the rare misaligned reads get clustered out as noise anyway."

This is acceptable for set2's goal (CPA site comparison across CPA-machinery depletion conditions).

If multi-aligner consensus is needed later (e.g., for the splice-y wt_tfiiib samples), the Stage E.5 per-chunk-consensus pattern or Commit B's region-pre-partitioning fix can be applied retroactively when axis 2 is properly resolved.
