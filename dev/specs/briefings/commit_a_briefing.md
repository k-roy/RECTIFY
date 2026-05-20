# Sonnet Briefing — Commit A (Shared parallel infrastructure)

**You are a Sonnet 4.6 subagent.** Your goal: land Commit A on M1 working tree (do **not** commit or push — Opus reviews + commits). All work happens at `/Users/kevinroy/work/rectify/`, branch `drs-validation-rebuild` (currently at `cb2fe6c`).

Read these files before starting:
1. `/Users/kevinroy/work/rectify/dev/specs/parallel_refactor_plan.md` §1 (Architectural overview), §4 Commit A subsection, and §6.5 (you do NOT wire any sidecars in this commit, but you must know what's coming so your APIs don't preclude it).
2. `/Users/kevinroy/work/rectify/dev/specs/profile_results.md` — the Commit Zero decision. Particularly: realign_exon_blocks is 88% of WRITE phase for QSrev, write_corrected_bam parallelism is the QSrev win, BGZF compresslevel changes are REFUTED and must NOT appear in this commit.
3. `/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/CLAUDE.md` — particularly the "Verify-don't-rewrite for agents" rule and "Surgical staging" rule. **Do not run `git add -A` or `git add .`** — but you also will not commit at all, so this is just a sanity rule.

---

## 1. Critical rules

1. **No commits, no pushes.** Modify working tree only. Opus reviews the diff via `git diff` and commits manually.
2. **No behavior change to existing code paths.** New modules exist but nothing routes through them yet. `correct_command.py:run` still calls the existing `write_corrected_bam`, not `write_corrected_bam_parallel`. That rewire happens in Commit B.
3. **Verify, don't recite.** Every claim about an existing function's location/signature/behavior must be backed by a tool call (Read, grep). The lab has a documented failure mode: agents that confabulate produce phantom flags 30% of the time. Don't be that.
4. **No `git add -A` / `git add .`.** Not committing anyway, but if you stash or unstash, stage surgically.
5. **Do NOT add provenance sidecars in this commit.** That's Commit A.5. Your APIs should not preclude it but should not require it either.
6. **Do NOT change `tracker.register_staged`.** Same reason.
7. **Do NOT add a `compresslevel=1` BGZF tweak.** Refuted by Commit Zero profiling (BGZF is 1.9% of WRITE phase, not the ≥30% threshold). The plan §3.3 explicitly drops this escape hatch.
8. **Tests must pass.** `pytest -m "not slow"` must still pass after your changes (no regression). The new smoke test you add must also pass.

---

## 2. Current state of the tree (verified, do not re-derive)

- M1 HEAD: `cb2fe6c fix(multi-aligner+analyze): normalize qnames in consensus; keep AG_RICH reads + per-cluster flag`
- M1 = GitHub `drs-validation-rebuild` (in sync).
- `rectify/core/bam/regions.py` already exists and contains `find_coverage_gaps()` + `get_processing_regions()`. **The spec's "extract get_processing_regions" step is already done.** You do NOT need to touch `parallel.py:35` (the existing import line).
- `rectify/core/bam/parallel.py:244` is `def process_bam_file_parallel(...)`. This is the existing region-parallel correction pipeline. Read it to understand the worker pattern you're mirroring.
- `rectify/core/bam/parallel.py:114` is `def _process_region_worker(...)`. This is the existing per-region worker. Read it to understand pysam-fetch + region-dedup conventions.
- `rectify/core/bam/bam_writer.py:200` is `def write_corrected_bam(...)`. This is the single-threaded function you're parallelizing in concept (but NOT wiring through yet).
- `rectify/core/splice/junction_scoring.py:200` is `def build_junction_pool(...)`. Sequential across aligner BAMs — you parallelize this.

---

## 3. Files to add or modify

### 3.1 ADD `rectify/core/bam/regions.py` (~30 LOC additions — append to existing file)

Add a `RegionPlan` dataclass + a helper that wraps `get_processing_regions()`:

```python
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

@dataclass(frozen=True)
class RegionPlan:
    """A single region's work plan for parallel BAM processing.

    Carries an id (stable across runs given the same input + min_gap_size)
    and a per-region temp dir so workers can write outputs without colliding.
    """
    region_id: str          # e.g. "r000", "r001"; left-padded so lex-sort == coord-sort
    chrom: str
    start: int              # 0-based, inclusive
    end: int                # 0-based, exclusive (pysam fetch convention)
    tmp_dir: Path           # per-region scratch dir (workers write region_NNN.bam etc. here)

    @property
    def region_bam(self) -> Path:
        return self.tmp_dir / f"{self.region_id}.bam"

    @property
    def ok_sentinel(self) -> Path:
        """Marker file written ONLY after region_bam is fully sorted + fsync'd.
        Existence -> this region's work is durably complete."""
        return self.tmp_dir / f"{self.region_id}.bam.ok"


def plan_regions(
    bam_path: str,
    tmp_dir: Path,
    min_gap_size: int = 10000,
    max_region_size: int = 100000,
) -> List[RegionPlan]:
    """Wrap get_processing_regions() into RegionPlan objects.

    Region ids are assigned in coordinate order using fixed-width zero-padding
    so that lex-sort matches coord-sort (matters for samtools merge of
    'region_*.bam' globbed in shell order).
    """
    raw = get_processing_regions(bam_path, min_gap_size=min_gap_size, max_region_size=max_region_size)
    width = max(3, len(str(len(raw) - 1)))  # at least 3 digits, more if needed
    plans = []
    for i, (chrom, start, end) in enumerate(raw):
        plans.append(RegionPlan(
            region_id=f"r{i:0{width}d}",
            chrom=chrom,
            start=start,
            end=end,
            tmp_dir=tmp_dir,
        ))
    return plans
```

Notes:
- `@dataclass(frozen=True)` so RegionPlan is hashable (needed if anyone uses it as a dict key).
- Fixed-width id (`r000`, `r001`, ..., `r099`) — keeps shell-glob lex order matching coord order.
- `region_bam` and `ok_sentinel` paths use `tmp_dir` directly; multiple parallel BAM-write calls would need different tmp_dirs (caller's responsibility).
- Existing `get_processing_regions()` returns tuples, not RegionPlans. Don't change its signature — `plan_regions()` wraps it.

### 3.2 NEW `rectify/core/bam/bam_writer_parallel.py` (~300-400 LOC)

This is the heart of Commit A. Three components:

**A) `_process_region_for_bam_write(plan, input_bam_path, corrections, genome, hp_penalty_table, ...)` — module-level function (NOT a method; must be picklable for `mp.Pool`).**

For each region:
1. Open `input_bam_path` with `pysam.AlignmentFile(input_bam_path, 'rb')`. Each worker MUST open its own — pysam handles do not survive fork.
2. Skip work if `plan.ok_sentinel.exists()`. Log `"[region {plan.region_id}] skip: sentinel exists"` and return early. This is the resume mechanism.
3. Open a coord-sorted output writer at `plan.tmp_dir / f"{plan.region_id}.unsorted.bam"` using `pysam.AlignmentFile(... 'wb', header=bam_in.header)`. **Use the input BAM's header verbatim** — no per-region @PG additions; the main thread will add one combined @PG line after merge.
4. Fetch reads scoped to the region: `for read in bam_in.fetch(plan.chrom, plan.start, plan.end):`
5. **Region-boundary dedup:** skip the read if `read.reference_start < plan.start` or `read.reference_start >= plan.end`. (A read whose alignment overlaps the region but whose start is in the adjoining region is NOT this region's responsibility — the adjoining worker will pick it up.) Reads with `read.is_unmapped` are also skipped here (unmapped reads will be handled by a separate "tail" pass; see (C) below).
6. Apply the same per-read mutations as `bam_writer.py:200-368`'s loop body. **Read that function first** to see the exact sequence of calls:
   - `realign_exon_blocks(...)` if applicable
   - `extend_read_5prime_for_junction_rescue(...)`
   - `softclip_intronic_tail_5prime(...)` or `reroute_intronic_tail_5prime_via_junction(...)`
   - `extend_read_3prime_for_softclip_rescue(...)`
   - `clip_read_to_corrected_3prime(...)`
   - `_hardclip_trailing_a_run(...)`
   - `set_tag('cp', ...)`
   - `_decode_eq_seq_inplace(...)` (lines around bam_writer.py:60-100)
   - `bam_out.write(read)`
   Mirror this sequence exactly. Do NOT refactor it — copy semantics. Any deviation is a correctness regression.
7. After the fetch loop, close `bam_in` and `bam_out`. The unsorted region BAM is now complete.
8. **Sort the unsorted region BAM in place** to `plan.region_bam`:
   ```python
   pysam.sort('-@', '1', '-o', str(plan.region_bam), str(unsorted_path))
   unsorted_path.unlink()
   ```
   Use 1 thread per worker — total CPU is already saturated by N workers, no benefit to multi-threading each sort.
9. **fsync the region BAM** to ensure durability before writing the sentinel:
   ```python
   with open(plan.region_bam, 'rb') as fh:
       os.fsync(fh.fileno())
   ```
10. **Touch the sentinel last:** `plan.ok_sentinel.touch()`. Sentinel existence MUST imply region_bam is durably complete (sentinel goes after fsync).
11. Return a small dict of stats: `{"region_id": plan.region_id, "n_reads_in": ..., "n_reads_out": ..., "n_reads_skipped_dedup": ..., "wall_seconds": ...}`.

**B) `write_corrected_bam_parallel(input_bam_path, corrections_table, output_bam_path, n_threads, genome, hp_penalty_table=None, tmp_dir=None, allow_resume=True, **kwargs)` — public entry point.**

Behavior:
1. Resolve `tmp_dir` (default: a new tempdir under `$TMPDIR/rectify_regions/<pid>_<random>/`). Create the dir if it doesn't exist.
2. Pre-flight: `shutil.disk_usage(tmp_dir).free >= 1.5 * os.path.getsize(input_bam_path)`. If not, raise `RuntimeError(...)` with the actual numbers.
3. Plan regions: `plans = plan_regions(input_bam_path, tmp_dir)`. Log the region count.
4. Load corrections table ONCE in the parent process. Pass it to each worker as a dict (copy-on-write via fork on Linux; explicit pickle on macOS where fork isn't default).
5. Spawn `mp.Pool(n_threads)`. Use `pool.imap` (NOT `imap_unordered`) so regions complete in coord order — this is what makes the final `samtools merge` heap-merge produce a sorted output without a re-sort pass.
6. Collect per-region stats into a list as workers complete.
7. **Sort-then-merge.** After all workers join:
   - Verify each `plan.region_bam` exists and `plan.ok_sentinel` exists. If any missing, raise — that worker silently failed.
   - Run: `pysam.merge('-@', str(n_threads), '-p', '-c', '-f', str(output_bam_path), *[str(p.region_bam) for p in plans])`
     - `-p` preserves the @PG chain.
     - `-c` consolidates colliding @PG IDs (which there won't be since all region BAMs share the same input header, but defensive).
     - `-f` overwrites if `output_bam_path` exists.
     - `-@ N` parallel decompression of inputs (cheap; the merge itself is already heap-merge).
   - The merge output is **fully coord-sorted** because every input region BAM is coord-sorted and regions were processed in coord order. **No final sort pass needed.**
8. Index: `pysam.index(str(output_bam_path))`.
9. Add a single combined @PG line to the merged output (samtools merge already does this).
10. If `not kwargs.get("keep_tmp"):` clean up `tmp_dir` recursively. Else log the path.
11. Return aggregated stats: `{"n_regions": ..., "n_reads_in_total": ..., "n_reads_out_total": ..., "wall_seconds_total": ..., "stats_per_region": [...]}`.

**C) Unmapped-reads tail pass.** Reads with `is_unmapped=True` are not fetched by `bam.fetch(chrom, start, end)`. After all region workers finish, you need ONE more pass:
- Fetch `bam_in.fetch(until_eof=True)` and filter to `is_unmapped` reads.
- Apply the same mutations (most don't apply to unmapped reads — they'll be near-pass-through).
- Write to a `tmp_dir / "unmapped.bam"`.
- Add it to the merge input list.

Look at `bam_writer.py:200-368` to see how the existing single-threaded function handles unmapped reads. Mirror that exactly. If unmapped reads aren't given special handling there (they may just flow through the fetch loop because the existing function iterates the BAM without a region restriction), document that and either: (a) handle them in a final separate worker, or (b) include them in the LAST region's worker. Either is fine; just be explicit.

**Imports:** carefully replicate from `bam_writer.py`. Don't import `bam_writer.write_corrected_bam` (no need; we're parallelizing in a separate module).

### 3.3 NEW `rectify/core/bam/tsv_partition.py` (~150 LOC)

```python
"""Per-region TSV emission + manifest.

Each region worker writes its own corrected_reads.region_NNN.tsv.
The manifest is the canonical artifact; a single concatenated TSV is
optionally emitted for back-compat via --emit-merged-tsv (which lives in
Commit B's correct_command wiring, NOT here).
"""
```

Components:

**A) `RegionTsvWriter` — context manager.**
```python
class RegionTsvWriter:
    def __init__(self, region_plan: RegionPlan, header: List[str]):
        self.plan = region_plan
        self.header = header
        self.path = region_plan.tmp_dir / f"{region_plan.region_id}.tsv"
        self._fh = None
        self._n_rows = 0
        self._sha256 = hashlib.sha256()  # hash as we write, no separate pass

    def __enter__(self):
        self._fh = open(self.path, 'w')
        line = '\t'.join(self.header) + '\n'
        self._fh.write(line)
        self._sha256.update(line.encode())
        return self

    def write_row(self, row: List[str]) -> None:
        line = '\t'.join(row) + '\n'
        self._fh.write(line)
        self._sha256.update(line.encode())
        self._n_rows += 1

    def __exit__(self, *args):
        self._fh.close()
        # Caller reads .n_rows and .sha256.hexdigest() for manifest entry

    @property
    def n_rows(self) -> int:
        return self._n_rows

    @property
    def sha256_hex(self) -> str:
        return self._sha256.hexdigest()
```

**B) `write_manifest(manifest_path, region_writers)` — main-thread function.**

After all workers complete, the parent collects each worker's (region_id, path, n_rows, sha256) and writes a TSV manifest:
```
region_id	chrom	start	end	tsv_path	n_rows	sha256
r000	chrI	0	100000	r000.tsv	12345	abc123...
r001	chrI	100000	200000	r001.tsv	23456	def456...
...
```
- Header on first line.
- Rows sorted by region_id (which is already coord-ordered).
- Paths are stored RELATIVE to the manifest file's directory — so moving the directory doesn't break the manifest. Use `os.path.relpath(tsv_path, manifest_path.parent)`.

**C) `load_manifest(manifest_path) -> List[Dict]` — reader.**

Reads a manifest, returns list of dicts (one per region). Verifies header schema. Used by Commit B's `analyze` partial-streaming + Commit A.5's PortablePath validation. Keep the API minimal — Commit D extends it.

### 3.4 MODIFY `rectify/core/splice/junction_scoring.py:200-236` (`build_junction_pool`) (~30 LOC change)

Current code (lines 226-230):
```python
for bam_path in aligner_bams:
    novel = collect_junctions_from_bam(bam_path, chrom_filter=chrom_filter)
    all_j.update(novel)
```

Replace with:
```python
if not aligner_bams:
    # No aligner BAMs — pool is just annotated set.
    pass
elif len(aligner_bams) == 1:
    # Single aligner — no fork overhead.
    novel = collect_junctions_from_bam(aligner_bams[0], chrom_filter=chrom_filter)
    all_j.update(novel)
else:
    # Multiple aligners — process in parallel (one process per BAM).
    from concurrent.futures import ProcessPoolExecutor
    n_workers = min(len(aligner_bams), os.cpu_count() or 4)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = [ex.submit(collect_junctions_from_bam, bp, chrom_filter) for bp in aligner_bams]
        for fut in futures:
            all_j.update(fut.result())
```

Add `import os` at the top if not already present.

**Add test** at `tests/test_junction_scoring_parallel.py`:
- Mock or use a tiny fixture: two BAM files with a known set of junctions each.
- Call `build_junction_pool([bam1, bam2], annotated_set={})`.
- Assert the returned set is the union of the two BAMs' junctions plus the annotated set.
- Verify the result is identical whether the function uses the sequential or parallel path (parametrize over n_workers=1 vs n_workers=2, or just confirm correctness once at len(aligner_bams)=2).

### 3.5 NEW `tests/utils/__init__.py` and `tests/utils/bam_compare.py` (~120 LOC)

```python
"""QNAME-sort + record-by-record BAM equivalence helper.

Used by all parallel-vs-legacy equivalence tests across Commits A-F.
"""
def assert_bams_equivalent(bam_a: str, bam_b: str, *, ignore_pg: bool = True) -> None:
    """Assert two BAMs are equivalent modulo:
      - sort tie order (we QNAME-sort both before comparison)
      - the @PG chain (if ignore_pg=True; the parallel path adds extra @PG lines)

    Compares: CIGAR, FLAG, POS, MAPQ, all tags, query_name, reference_name.

    Raises AssertionError with a precise location on first divergence.
    """
    # Implementation:
    # 1. With pysam, fetch all reads from bam_a + bam_b
    # 2. Sort each list by (query_name, reference_name, reference_start, is_reverse)
    # 3. Assert lengths match (n_reads). If not, raise with both counts.
    # 4. For each pair: compare CIGAR string, FLAG, POS, MAPQ, all tags by name.
    #    If tag set differs, list which tags are present in A but not B.
    # 5. If everything matches, return None.
```

API design:
- One public function: `assert_bams_equivalent(a, b)`.
- A private helper to sort reads.
- A precise `_compare_one_read(read_a, read_b)` that returns a list of difference strings (empty if equal).
- The assert message should include the offending QNAME and a summary of the first difference.

### 3.6 NEW `tests/test_bam_writer_parallel_smoke.py` (~100 LOC)

```python
"""Smoke test: write_corrected_bam_parallel produces output equivalent to
the existing single-threaded write_corrected_bam on a small fixture."""
```

Test strategy:
1. Find an existing small BAM fixture in `rectify/data/validation/` — there should be one with ~100-1000 reads. Use whichever exists; don't create a new one.
2. Find or create a tiny `corrected_reads.tsv` for that fixture (the existing `bam_writer.write_corrected_bam` requires this input). If one doesn't already exist for the fixture, generate one by running `rectify correct` on the fixture once at test setup (cache the result).
3. Run the legacy path: `write_corrected_bam(fixture_bam, corrected_tsv, "/tmp/legacy.bam", ...)`.
4. Run the new path: `write_corrected_bam_parallel(fixture_bam, corrected_tsv, "/tmp/parallel.bam", n_threads=2, ...)`.
5. Assert equivalence: `assert_bams_equivalent("/tmp/legacy.bam", "/tmp/parallel.bam")`.

Skip the test (with a clear `pytest.skip` reason) if the fixture doesn't exist. Don't ship a broken test.

Mark with `@pytest.mark.parametrize("n_threads", [1, 2, 4])` so we exercise multiple worker counts.

---

## 4. Acceptance criteria (every box must be checked before reporting done)

- [ ] `rectify/core/bam/regions.py` has `RegionPlan` dataclass + `plan_regions()` function added.
- [ ] `rectify/core/bam/bam_writer_parallel.py` exists with `write_corrected_bam_parallel()` + `_process_region_for_bam_write()` (module-level, picklable).
- [ ] `rectify/core/bam/tsv_partition.py` exists with `RegionTsvWriter` + `write_manifest()` + `load_manifest()`.
- [ ] `rectify/core/splice/junction_scoring.py:build_junction_pool` uses `ProcessPoolExecutor` when `len(aligner_bams) > 1`, falls back to sequential otherwise.
- [ ] `tests/utils/__init__.py` and `tests/utils/bam_compare.py` exist with `assert_bams_equivalent()`.
- [ ] `tests/test_bam_writer_parallel_smoke.py` exists and PASSES on M1.
- [ ] `tests/test_junction_scoring_parallel.py` exists and PASSES.
- [ ] `pytest -m "not slow"` passes overall — NO regressions on the existing ~934-test suite. Run from `/Users/kevinroy/work/rectify/` with `source $HOME/.pyenv/versions/rectify/bin/activate` or whichever env Kevin uses. If you can't find the env, ask via TaskUpdate; do not guess.
- [ ] `correct_command.py:run` still calls the existing `write_corrected_bam` at line ~848. Do NOT change that line. (`git diff rectify/core/commands/correct_command.py` should show no changes from this commit.)
- [ ] All imports use absolute paths (`from rectify.core.bam.regions import ...`), not relative (`from .regions import ...`), unless the file you're editing already uses relative — then match local convention.
- [ ] No `git add` or `git commit` runs. Working tree changes only.
- [ ] Final `git status` output is captured in your handoff note so Opus can review what changed.

---

## 5. When to escalate back to Opus (stop and ask)

Hand back immediately if:
- You can't find the test fixture BAM in `rectify/data/validation/`. Don't fabricate one; ask.
- The legacy `write_corrected_bam` signature doesn't match the spec's description. Don't bend the spec to match a different reality — surface the discrepancy.
- `assert_bams_equivalent` fails on the smoke test and you can't determine if it's a real divergence vs. a tag ordering artifact. Don't paper over.
- `pytest -m "not slow"` regresses. The regression IS the bug — do not work around it by skipping tests.
- The pysam version on M1 lacks `pysam.merge`. Check `import pysam; pysam.merge.__doc__` first. If absent, hand back.
- Any cross-cluster sync state seems off. M1 should be at `cb2fe6c`; if `git log --oneline -1` shows something else, hand back.
- A subtle correctness concern emerges (e.g., region-boundary read crossover that violates `read.reference_start in [plan.start, plan.end)`). Read `parallel.py:_process_region_worker` to see how the existing function handles this and mirror it.

---

## 6. Reporting back

When done (all acceptance boxes checked), report via this conversation:
1. Final `git status -s` output (untracked + modified files).
2. `git diff --stat` summary (lines added/removed per file).
3. The pytest summary line for `pytest -m "not slow"` ("N passed, M skipped" etc.).
4. The smoke test's per-parameter output (n_threads=1 / 2 / 4 each passed).
5. Any rough edges you papered over with a TODO that Opus should address before Commit B.

Do NOT report success if any acceptance box is unchecked.

---

## 7. Time budget

- Reading the relevant source files: ~30 min.
- Implementing `bam_writer_parallel.py`: ~60-90 min.
- Implementing `tsv_partition.py`: ~30 min.
- Implementing `bam_compare.py` + smoke test: ~30-45 min.
- Implementing `junction_scoring.py` change + test: ~15 min.
- Running pytest, fixing regressions, iterating: ~30-60 min.
- Total: ~3-4 hours.

**End of briefing.**
