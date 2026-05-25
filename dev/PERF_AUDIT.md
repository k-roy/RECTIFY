# PERF_AUDIT.md — hunting per-read over-computation in RECTIFY

**Created:** 2026-05-24 (Kevin's request). **Status:** living doc.

**Why this exists:** the 3'SS-rescue bottleneck below (run-all manifest inline-correct
stall → OOM) is almost certainly *not* the only place RECTIFY does far more work than
necessary. This doc (a) records that case as a worked example, and (b) gives the next
agent a concrete methodology + grep targets + suspect list to find and fix similar
bottlenecks. **Read this before any perf work; add findings as you go.**

RECTIFY is v0.9.x; manifest mode is experimental. These are perf bugs, not correctness
bugs — but at production scale (full junction pools, high-coverage loci, eventually
human introns >100 kb) they make pipelines hang. Fix them, but verify they don't change
output (see Verification below).

---

## Case study: 3'SS rescue = O(reads × full-pool × ambiguity-window) HP-edit-distance

**Where:** `rectify/core/splice/splice_aware_5prime.py` (`_rescue_3ss_truncation_body`,
`_hp_edit_distance`) fed by `rectify/core/bam/bam_processor.py:correct_read_3prime`.

**Symptom:** `rectify run-all --manifest` (DRS, 5% wt+upf1 subset) hung at
`Processing 126 regions across 16 workers / 0%` for 6.5 h then OOM-killed at 64 GB
(Sherlock job 25846844). Align + trim + Module-2H refine all succeeded; the stall is the
inline `correct` parallel region step.

**Cost decomposition (py-spy, 29,999 samples, single-thread full-pool subset = job DIAGE):**
`_hp_edit_distance` (an O(n×m) DP, splice_aware_5prime.py:693-723) = **~86% of CPU**.
- Baseline rescue `_rescue_3ss_truncation_body` = **87%** of that. Two loops scan the
  candidate set and call the DP per candidate: the sequence-rescue loop (line ~1409) and
  the Case-4 intronic-snap loop (~1877/1878).
- The sequence-rescue loop is worse than "per candidate": for EACH candidate it runs
  `for _shift in [±~15]` × `for _off in [0..10]` ≈ **~341 DP calls per candidate**.
- Terminal peel `_terminal_peel_rescue` = 12.9% (its own per-depth scan, fixed first).

**Three compounding root causes — the generalizable lessons:**
1. **Per-read scan of a global collection.** `candidate_junctions` (the junction pool)
   can hold ~17k entries in run-all mode; the body looped over all of them per read.
2. **Vestigial oversized bound.** `_POOL_SEARCH_RADIUS` was **10000** bp, with a comment
   claiming it must cover `junction_proximity_bp` "default 5000" — but the real default
   is **10**. The 10k radius was sized for a parameter that no longer exists, so it pulled
   every intron in a ±10 kb window per read.
3. **Bound keyed on the wrong axis.** The pool index is sorted by `intron_start`; a radius
   on `intron_start` bounds the *intron length*. That happens to work for yeast (introns
   ≤~1.2 kb) but **breaks for human** (introns >100 kb): the donor sits 100 kb from the
   acceptor near the read, so it falls outside any sane radius. The correct axis is
   **splice-site proximity to the read's 5' boundary**, leaving the *other* site free.

**Fixes (status as of 2026-05-25):**
- DONE (commit 8e8dc8c): peel candidate-narrowing + depth cap. **Eliminated the OOM**
  (full-pool RSS 64 GB → ~14 GB, verified by live /proc probe).
- DONE (uncommitted): body-level `_nearby_junctions` narrowing (union of the loops' own
  cheap gates, computed once) + `_POOL_SEARCH_RADIUS` 10k→2k. Correctness-clean
  (246 focused tests). **Did NOT fix wall-time**: DIAGF (full-pool, fixed code) still
  stalled at 0/126 regions after 9.5 min — because each surviving nearby candidate still
  triggers the ~341-call ambiguity block, and the full pool leaves more nearby candidates
  than GFF-only.
- DONE (commit `961c844`): **dual-site index + either-site fetch.** Each junction is
  indexed by BOTH splice sites (per-chrom `IntervalTree`); the fetch surfaces any junction
  with *either* site within the peeled-back genomic span + buffer, NOT a radius on one
  site. Intron-length-independent (human-ready) AND tighter. The restructure is localized
  to `_build_pool_chrom_index` + the `correct_read_3prime` fetch.
- DONE (commit `ed3df74`): **pool-build anchor floor + cross-family concordance
  relaxation** — cuts the candidate *count* at its source (the complement to the fetch
  narrowing above). A novel observed junction enters the shared rescue/refinement pool only
  if (a) some read crosses it with a clean ≥10 bp exon anchor on both flanks (D/N-collapsed
  first), OR (b) it is reported by ≥2 *independent algorithm families* (`ALIGNER_FAMILY`:
  minimap2/gapmm2 = one family since gapmm2 wraps minimap2, vs uLTRA / BBMap / deSALT).
  Drops single-family tiny-anchor splits (gapmm2 `4M4250N223M`-style) without losing real
  short-exon-1 junctions. **Verified on the 10k yeast 5-aligner set (committed code):**
  pool 648 (no floor) → 245 (floor only) → 287 (floor+relax); the 361 net-dropped junctions
  are **all single-family** (zero ≥2-family losses) and the 42 relaxation-recovered are
  **all ≥2-family** — i.e. the floor removes only junctions lacking independent
  corroboration. 1315 fast tests green; 800-read Sherlock end-to-end verify earlier showed
  all dropped rescues spurious.
- DEFERRED to human DRS (handoff `dev/handoffs/HANDOFF_2026-05-25_human_drs_validation.md`
  §7): a graded **periodicity/complexity** dimension for the consensus soft gate
  (`_add_chimera_flag`). Helpers built (`_seq_periodicity`,
  `_junction_worst_flank_periodicity`) but **not wired** — it is a no-op on yeast (≥0.9
  periodicity fired on 0/210 gapmm2-only and 0/139 real junctions; the gene-dense genome
  has no long low-complexity anchors). It must prove its worth on human (repeats, long
  introns) before wiring. **No splice-motif gating** anywhere — pruning is anchor-quality /
  family-support only, to preserve unbiased non-canonical junction discovery.

**Note:** the advisor's standing guidance is *don't* touch the shift×offset ambiguity
windows themselves (intentional splice-site-slide search). Reduce the candidate *count*
that enters them.

---

## Anti-patterns to hunt for elsewhere

1. **Per-item scan of a global collection** — `for x in <big_set>:` inside a per-read /
   per-record function, where a cheap gate (`if dist > ...: continue`, containment) rejects
   most iterations. Fix: pre-narrow ONCE via a sorted/bisect index, interval tree, or a
   single-pass filter, then loop the small subset.
2. **Nested windows multiplying an expensive op.** An O(n×m) DP / local alignment / edit
   distance inside ≥2 nested loops (shift × offset × candidate). The multiplier hides the
   real cost; cut the outer count.
3. **Vestigial / oversized constants.** A magic number whose comment cites a default that
   no longer exists (the 10k-for-5000 case). Grep magic constants; re-derive against the
   *current* defaults.
4. **Bound keyed on the wrong axis.** A radius/window that scales with the wrong quantity
   (intron length instead of splice-site proximity). Ask: does this break for human / long
   features? If the bound is on coordinate A but the biology cares about coordinate B,
   it's wrong.
5. **Recompute-per-item what's constant-per-read/region.** Hoist read-level work out of
   per-candidate loops (RECTIFY already did this once — the `_leading_del` flag hoist,
   "Bug 3 patch perf fix" — look for more).

---

## How to hunt (methodology — profile, don't guess)

We theorized the cause twice and were wrong both times; only the sampling profile settled
it. Always measure.

- **py-spy (installed in the Sherlock rectify env).** `gdb py-bt` is UNUSABLE on this
  conda python3.9 (DWARF v5 vs gdb v2-4). Recipe:
  - Launch the target single-threaded (`--legacy-single-threaded -j 1`) in the background,
    `sleep 25` to pass setup, then
    `py-spy record -f raw -o stacks.folded --pid $PID --duration 300 --rate 100 --nonblocking`.
  - Aggregate by leaf frame:
    `awk '{c=$NF; s=$0; sub(/ [0-9]+$/,"",s); n=split(s,a,";"); sum[a[n]]+=c} END{for(k in sum) print sum[k]"\t"k}' stacks.folded | sort -rn | head`.
  - Split by caller: grep the folded lines for the suspect frame, sum counts with/without
    a calling-frame substring.
- **faulthandler for hangs:** `timeout -s ABRT <sec> python -X faulthandler -m rectify.cli ...`
  prints the current stack on the SIGABRT. Good for "where is it stuck right now."
- **Reproduce at scale on a dense input.** Synthetic unit tests miss this entirely (they
  use tiny pools). Use a real high-coverage region. Diag artifacts (subset BAM, 5 aligner
  BAMs, py-spy folded stacks, sbatch templates) live in
  `$SCRATCH/rectify_runall_diag_20260524/` and the trimmed/aligned BAMs in
  `$SCRATCH/rectify_wt_by4742_rep1_25846844_0/` on Sherlock.
- **grep targets to start from:**
  - `for .* in candidate_`, `for .* in .*junctions`, `for .* in .*pool`
  - nested `for _shift` / `for _off`; call sites of `_hp_edit_distance`,
    `align_clip_to_exon`, anything named `*edit_distance*` / `*align*`
  - `_SEARCH_RADIUS`, large integer literals as window/radius bounds
  - `.get_aligned_pairs()` or genome slicing inside per-candidate loops
  - intervaltree `.search` / `.overlap` / `search_point` called per read

---

## Suspect hotspots to audit next (from the DIAGE profile + code structure)

These appeared in the profile (minor here, but per-read and worth checking at other scales
/ organisms) or are structurally similar to the 3'SS case:

- **`walkback_3prime_guarded` / `walkback_drs_full`** (`core/correct/walkback.py`) — showed
  up in the DIAGE profile; per-read poly-A walkback.
- **gene-attribution interval-tree** — `search_point (intervaltree/node.py)` appeared in
  the profile; per-read interval-tree lookups for gene assignment (`gene_interval_trees`).
- **`indel_corrector`** (`detect_indel_artifacts`, `extract_insertions/deletions`,
  `utils/alignment.py`) — per-read, in the profile tail.
- **`junction_refiner`** (Module 2H) — builds + scans the junction pool; already partly
  bisect-optimized, but re-audit the per-N-op candidate search for the same radius-vs-axis
  issue, especially for human.
- **variant-aware rescue scan** and **per-aligner consensus scoring** — any path that
  receives `annotated_junctions` / the pool in full and loops it per read.

---

## Verification discipline (REQUIRE before declaring any perf fix "done")

A passing unit test is NOT evidence a perf fix works — the tests use tiny inputs. Require:
1. Focused blast-radius tests green, then `pytest -m "not slow"`.
2. **Scale re-run COMPLETES** on the dense representative input AND peak memory stays
   bounded (target <32 GB for the 5% subset).
3. **Measured** wall-time / call-count drop (e.g. count `_hp_edit_distance` calls before vs
   after; must drop ~10× in a dense region) — not "it feels faster."
4. **Output unchanged** — rescue counts / corrected-TSV within tolerance of the pre-fix
   run (no over-narrowing that silently drops real rescues).

---

## Pointers

- `AGENT_FIXES.md` → `[2026-05-24] ... 3'SS rescue` entry (full diagnosis + commands).
- Commit `8e8dc8c` (peel fix, OOM eliminated). Task #10 (baseline narrowing, in progress).
- Diagnostic jobs A–F + py-spy folded stacks: `$SCRATCH/rectify_runall_diag_20260524/`.
