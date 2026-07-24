# Scope: reclaiming ALL previous rectify outputs (H2 scratch) — 2026-07-23

**Ask:** "compress all previous rectify outputs … open up free space soon."
**Bottom line:** the fast CMA/aligner_chunks dedup lever is **already banked (~98 GB, decay run)**. There is **no
large, clean, mechanical reclaim left.** The remaining big space is **primary data (off-limits)**, **data-project
intermediates (your retention call)**, and **stale scratch (low-risk, needs your ack)** — not a CMA target.
Evidence: H2 `~/m0_scratch/rectify_scope_20260723_1448.txt` + `rectify_scope2_20260723_1452.txt` (read-only).

Scope boundary: H2 scratch only. `/u/home` + `/u/project` have **no chunked-rectify layout**. Oak/Sherlock/SCG =
SGTC/Sumner workspaces (separation rule) — **not scanned/crossed**.

## Ledger (H2 `/u/scratch/k/kevinroy/`), sorted by reclaim tractability

| Target | GB | Kind | Redundant vs durable product? | Verdict |
|---|--:|---|---|---|
| `408d_plants/raw` | 69 | primary sequencing reads | **NO — this IS the durable data** | 🔴 OFF-LIMITS |
| `invitro_rnt1/fastq` | 19 | primary fastq | NO | 🔴 OFF-LIMITS |
| `compare_pare/fastq`+`pare_invitro_381/fastq` | 18 | primary fastq | NO | 🔴 OFF-LIMITS |
| `sy14_p5b_retired` (SRR3706786x ×3) | 47 | **retired** SRA downloads | re-downloadable from SRA | 🟡 YOUR CALL (likely safe) |
| `408d_plants/trim` | 42 | trimmed fastq (from raw) | regenerable if trim is deterministic | 🟡 YOUR CALL (data project) |
| `398l_quantseq_epap/bam` | 29 | short-read BAMs | needs an inside look | 🟡 YOUR CALL |
| `408d_plants/bam` | 16 | derived alignments | derived from trim | 🟡 YOUR CALL |
| **orphaned `rectify_*` array scratch** | **31** | leftover per-chunk BAMs, DECAY samples | **PARTIAL** — Dcp2_repA & Xrn1_repB/minimap2 `==merged`; **Xrn1_repB/uLTRA = 3.2× merged (le=NO)** | 🟡 STALE scratch — reclaim needs your ack OR a per-job dedup pass |
| decay chunk `*.fastq.gz` (~30 GB across samples) | ~30 | regenerable via `rectify split` | regenerable (deterministic) | 🟡 low-risk, your call |
| decay `merged_bams` (2-aligner) | ~31 | retained pre-correct store | CMA lever A only ~1.5–2× → ~10–15 GB | ⚪ small, optional |

## Why the orphaned `rectify_*` (31 GB) is NOT a clean sweep
`rectify_<jobid>_<task>/` are SLURM array-task SCRATCH dirs (pre-gather per-chunk BAMs). By design they're
ephemeral; the durable store is each run's `merged_bams` + `consensus`. Gate result (scattered reads ≤ merged):
```
Dcp2_repA/minimap2  3158533 == 3158533  YES     Xrn1_repB/minimap2  2030361 == 2030361  YES
Dcp2_repA/uLTRA      514857 ==  514857  YES     Xrn1_repB/uLTRA     3334264 >  1037369  NO ⚠
```
The read-count `NO` meant those dirs accumulated chunks from **multiple/superseded Xrn1_repB attempts** (173 dirs,
`_of_039`) — DUPLICATE copies, not extra reads.

**RESOLVED by QNAME-subset verification (H2 `~/m0_scratch/orphan_qname_verify2.txt`, rc=0, 2026-07-24):** all
**15 distinct (sample,aligner) pairs across every `rectify_*` BAM are `SAFE=YES` with `missing=0`** — every
scattered UNIQUE read QNAME is already present in that sample's verified `merged_bams`. The earlier read-count
`le=NO` (Xrn1_repB/uLTRA 3.2×) was a FALSE ALARM: `scattered_uniq=1037369 == merged_uniq=1037369, missing=0`.
Several samples are partial subsets (e.g. total_mRNA_digested/minimap2 scattered_uniq=1.38M ⊆ merged 3.85M) —
still `missing=0`, still safe. ⇒ **the full ~31 GB of orphaned `rectify_*` scratch is verified-redundant and
reclaimable** (delete list build + completeness cross-check in progress). Kevin OK'd "verify-then-delete".

## ✅ EXECUTED: orphaned `rectify_*` reclaim (2026-07-24)
Kevin chose "verify-then-delete". QNAME-subset verify passed (all 15 pairs `missing=0`). Completeness
cross-check + tightening produced an explicit list (H2 `~/m0_scratch/orphan_safe_delete_dirs.txt`):
- **397 dirs = 28 GB** staged (mv→`_TRASH_orphan_rectify_20260724`) then purged. `moved=397=expected`.
- **10 dirs EXCLUDED** (NOT deleted): 5 `rectify_single_141133*` (ysh1_mdj1/ysh1_antisense outputs, non-decay,
  unverified), 3 empty `rectify_single_13761452_*`, `rectify_bench_20260703` (benchmark), **`rectify_deploy`**
  (deploy dir — a bad-delete narrowly avoided by requiring affirmative decay-chunk content).
- Decay `merged_bams`+`consensus` (separate path) confirmed UNTOUCHED. Evidence: `~/m0_scratch/orphan_purge_*.txt`.

**Running total reclaimed this effort: ~98 GB (decay aligner_chunks) + ~28 GB (orphaned scratch) = ~126 GB.**

## Recommended path to "free space soon" (biggest, lowest-risk first)
1. **`sy14_p5b_retired` (47 GB)** — if you confirm it's retired and re-downloadable, this is the single biggest
   clean win. One `mv → .TRASH → purge`. **Needs your OK.**
2. **Orphaned `rectify_*` scratch (31 GB)** — either (a) you ack "these are stale array scratch, delete them," or
   (b) I run a per-dir dedup pass (verify each scattered read is a duplicate of a merged read) to reach the
   verified-redundant standard, then reclaim. **Your call which.**
3. **Data-project intermediates** (`plants/trim` 42 GB, `398l/bam` 29 GB, `plants/bam` 16 GB) — regenerable/derived
   but they belong to analyses; **retention decision is yours**, per project.
4. **Primary fastq (~106 GB)** — do NOT delete; back up/offload if space-critical.

**Nothing here is deleted.** Awaiting your decisions on rows 1–3.
