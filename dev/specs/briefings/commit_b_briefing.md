# Sonnet Briefing — Commit B (DRS default flip + first stage sidecar wiring)

**You are a Sonnet 4.6 subagent.** Goal: land Commit B on M1. This is where the parallel infrastructure from Commit A becomes the DEFAULT code path, the first stage-level provenance sidecars get emitted, and resume becomes user-visible for the `correct` and `restore_polya` stages.

**Hard prerequisites** — verify all THREE before doing anything else:
1. Commit A (`feat(parallel): Commit A - shared parallel BAM-write infrastructure`) is in `git log --oneline -5`. Sha was `39fbb63` at draft time but may have moved.
2. Commit A.5 (`feat(provenance): ...`) is in `git log --oneline -5`. Sha TBD when A.5 lands.
3. Cross-cluster heads match (M1 = GitHub = H2 = Sherlock). Probe Sherlock ControlMaster first: `ssh sherlock 'echo ok'` returns within 2s without a Duo prompt — if it hangs, STOP and tell Opus.

If any prereq fails, halt and report to Opus.

Read these before starting:
- `/Users/kevinroy/work/rectify/dev/specs/parallel_refactor_plan.md` §0 (success metric — note the PROTOCOL-CONDITIONAL split: QSrev ~6× speedup, DRS ~1.1-1.3×), §4 Commit B, §6.1 (`corrected_reads.tsv` back-compat policy), §6.5 (sidecar schema), §11.1 (Opus owns first sidecar wiring).
- `/Users/kevinroy/work/rectify/dev/specs/profile_results.md` — confirms region parallelism is the right knob for QSrev/cDNA-short; refutes BGZF compresslevel as a cheap-win.
- `/Users/kevinroy/work/rectify/dev/specs/briefings/commit_a_briefing.md` and `commit_a5_briefing.md` for context on what landed.
- `/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/CLAUDE.md` — "Verify-don't-rewrite" rule, "Surgical staging" rule, "Cross-cluster commit-status check" rule.

---

## 1. Pre-recorded Opus decisions

### 1.1 Back-compat: manifest-only by default, with `rectify export-merged-tsv` shim

**Decision (Kevin, 2026-05-20, Option B):** flip the default to manifest-only. `--emit-merged-tsv` opts back in for users who want both files. Ship a new `rectify export-merged-tsv <manifest>` shim subcommand so back-compat consumers can recreate the merged TSV on demand. Migrate internal loaders so they accept EITHER format transparently.

Reasoning behind the deviation from §6.1's "if >5 consumers keep emit ON" rule:
- The grep found 20+ consumers but most are TRANSITIVE: internal rectify code paths that ultimately load through `_load_corrections_from_tsv` (`bam_writer.py`) or `load_corrected_positions` (`loaders.py`). Update the loader → all internal consumers work.
- External consumers (the 6 scripts and the H2 lab shell scripts) can opt in via `--emit-merged-tsv` OR use the export shim. Cost = one flag or one extra command per workflow.
- Keeping the merged TSV emitted by default permanently anchors a back-compat tax on every future commit. Better to take the cutover hit now while the parallelism refactor is the headline change.

Architectural shape:
- Default: `corrected_reads.manifest.tsv` + per-region `corrected_reads.region_NNN.tsv` files.
- `--emit-merged-tsv` (off by default): additionally emit `corrected_reads.tsv` (concatenated from the manifest, ~5% of runtime).
- `rectify export-merged-tsv <manifest> [-o output.tsv]`: subcommand that consumes a manifest and emits the concatenated TSV. ~30 LOC. Idempotent.
- Internal loaders auto-detect: if input path ends with `.manifest.tsv`, dispatch to manifest-aware load; else single-TSV load.
- Deprecation warning when a consumer reads a single `corrected_reads.tsv` produced by the pre-cutover pipeline (just log; don't error).

### 1.2 First stage-level sidecar — the template for Commits D, E, F

This commit emits the FIRST two stage-level sidecars per the §6.5.2 schema:
- `<sample_output>/<sample_id>.correct.provenance.json`
- `<sample_output>/<sample_id>.restore_polya.provenance.json`

The exact pattern you establish here is what Commits D (analyze), E (drs_trim + align + cdna_stage1), and F (NETSEQ) will copy. Get it right. Specifically:
- Use `ProvenanceRecord.from_components()` from `rectify.core.provenance.sidecar` to build the record from raw inputs/outputs/argv. Do NOT manually populate fields.
- Use `write_stage_sidecar()` to write. It's atomic (`.tmp` + rename); errors that crash mid-write leave no partial sidecar.
- Skip-check at stage entry via `can_skip_stage()` (lazy-imports sidecar internally, no circular import). Returns `SkipDecision.SKIP` or `SkipDecision.RUN(reason=...)`.
- Add a `[RESUME] sample={sid} stage={stage} decision={SKIP|RUN} reason={reason}` log line on every stage entry.
- Sidecars MUST NOT be written until the stage's FINAL output (merged BAM + manifest + merged TSV if emitting) has been successfully created. Mid-flight crashes leave region-level `.ok` sentinels but no stage-level sidecar → next run reruns the stage but workers skip done regions (cheap rerun).

### 1.3 Small-input gate threshold

Hard-coded constants in `rectify/core/bam/bam_writer_parallel.py:write_corrected_bam_parallel` (or in `correct_command.py:run` if you prefer the dispatch live there):
- If `input_bam_size_bytes < 100 * 1024**2` (100 MB) OR `estimated_read_count < 500_000` → route to sequential `write_corrected_bam`.
- Estimated read count: `pysam.idxstats(input_bam)` returns per-chrom counts; sum them. Fast.

This gate prevents the existing test suite from regressing on ProcessPool spawn cost. Document it in the function docstring with a sentinel log line: `"[bam_writer] small-input gate fired (size=%d MB, est_reads=%d) — using sequential path"`.

### 1.4 `--legacy-single-threaded` escape hatch

Add a flag to `rectify correct` (and `rectify run-all`) that forces the pre-Commit-A code path. Routes to existing `write_corrected_bam` regardless of input size. For:
- Debugging: bisecting correctness regressions.
- CI environments with insufficient cores.
- Memory-constrained environments (e.g., M1 with other heavy work).

Document with a deprecation note: "This flag will be removed in a future release once the parallel path is the default for at least one minor version."

---

## 2. Critical rules

1. **NO git commits, NO pushes.** Working tree only. Opus reviews + commits.
2. **Verify, don't recite.** Per CLAUDE.md.
3. **No `git add -A`.** Surgical only (you won't commit anyway).
4. **No changes to `rectify/utils/provenance.py`.** That's the legacy run-log accumulator; A.5 left it untouched and Commit B does too.
5. **Watchdog mitigation:** the Commit A subagent stalled at 600s due to pytest transcript silence. Use `pytest -v --tb=line` and print ~60s progress updates while long commands execute. Don't go silent.
6. **`pytest -m "not slow"` MUST pass.** Current baseline: 902 + whatever A.5 added. Plus your new resume tests.
7. **Han wt_R1 smoke** is the headline acceptance run (see §3.4). Budget ~45 min wall on a 16-core H2 node.
8. **Sidecar emission is correctness-critical.** If a stage emits a sidecar BEFORE its outputs are durably written, resume can return SKIP on a stage that's actually incomplete. Verify the ordering via a SIGTERM injection test (see §3.6).

---

## 3. File-by-file work

### 3.1 MODIFY `rectify/core/commands/correct_command.py`

**Verify current state first:**
```bash
grep -n "def run\|write_corrected_bam\|restore_polya\|process_bam_file_parallel" \
  rectify/core/commands/correct_command.py | head -30
```
Spec cites line numbers (run:339, write_corrected_bam call:848, restore_polya:789) from `cb2fe6c`; HEAD may have moved them. Use the grep output as truth.

**Changes:**

A. **Add resume flags via `rectify.core.provenance.cli.add_resume_args(parser)`** in `create_correct_parser()`. Single line; argparse helper already does the work.

B. **Add `--emit-merged-tsv` and `--legacy-single-threaded` flags** to the parser. Defaults: `emit_merged_tsv=False` (manifest-only by default per §1.1), `legacy_single_threaded=False`.

C. **In `run(args)`:** at stage entry, before any work, call:
```python
from rectify.core.provenance import can_skip_stage, detect_current_cluster
import sys, subprocess

decision = can_skip_stage(
    stage='correct',
    sample_output=args.output_dir,  # or wherever the stage's outputs land
    sample_id=_derive_sample_id(args),
    current_inputs={
        'bam':       args.input_bam,
        'genome':    args.genome,
        'annotation': args.annotation,
        # plus prescan caches if applicable
    },
    current_argv=sys.argv,
    rectify_git_sha=_get_rectify_git_sha(),  # see §3.5 below
    force=args.force_all,
    force_stages=set((args.force_stage or '').split(',')) if args.force_stage else None,
    accept_prior_provenance=args.accept_prior_provenance,
)
logger.info("[RESUME] sample=%s stage=correct decision=%s reason=%s",
            sample_id, 'SKIP' if decision.skip else 'RUN', decision.reason)
if decision.skip:
    # Reuse prior outputs; return early.
    return 0
```

D. **Dispatch to parallel vs sequential write:**
```python
input_size = os.path.getsize(args.input_bam)
n_reads_est = _estimate_n_reads(args.input_bam)  # pysam.idxstats sum
use_parallel = (
    not args.legacy_single_threaded
    and input_size >= 100 * 1024**2
    and n_reads_est >= 500_000
)
if use_parallel:
    stats = write_corrected_bam_parallel(args.input_bam, corrected_tsv, output_bam,
                                          n_threads=args.n_threads, ...)
else:
    if args.legacy_single_threaded:
        logger.info("[bam_writer] using legacy single-threaded path (--legacy-single-threaded)")
    else:
        logger.info("[bam_writer] small-input gate fired (size=%d MB, est_reads=%d) — sequential",
                    input_size // 1024**2, n_reads_est)
    stats = write_corrected_bam(args.input_bam, corrected_tsv, output_bam, ...)
```

E. **TSV emission switch:** the parallel correction loop (`process_bam_file_parallel` / streaming variants) currently writes a single concatenated `corrected_reads.tsv`. Switch to per-region writes via `RegionTsvWriter` (from `rectify.core.bam.tsv_partition`), then at the end:
- ALWAYS write `corrected_reads.manifest.tsv` (the canonical artifact going forward).
- If `args.emit_merged_tsv` (default FALSE per §1.1): also write the concatenated `corrected_reads.tsv` for back-compat. Use `cat`-style concat from the manifest's region TSVs — cheap (~5% of total runtime).
- Log a one-line `[CHANGED] <stage> now writes corrected_reads.manifest.tsv by default; pass --emit-merged-tsv for the legacy concatenated form or use 'rectify export-merged-tsv <manifest>' on demand.` on the first run after upgrade, so users notice the default flip.

F. **`restore_polya` parallel path:** the same `bam_writer_parallel` infra applies (per-read CIGAR rewrite, same shape as `write_corrected_bam`). Route through it conditionally on the same gate as (D). Emit its own sidecar (`<sample_id>.restore_polya.provenance.json`).

G. **At stage exit**, before returning:
```python
from rectify.core.provenance import ProvenanceRecord, write_stage_sidecar

record = ProvenanceRecord.from_components(
    stage='correct',
    stage_subtype=_derive_subtype(args),  # 'drs' / 'cdna' / 'short' / etc.
    sample_id=sample_id,
    rectify_git_sha=rectify_git_sha,
    inputs={
        'bam': args.input_bam,
        'genome': args.genome,
        'annotation': args.annotation,
        # plus prescan caches if applicable
    },
    outputs={
        'corrected_bam': output_bam,
        'corrected_bam_index': output_bam + '.bai',
        'corrected_tsv_manifest': manifest_path,
        ('corrected_tsv_merged', {'optional': True}): merged_tsv_path if args.emit_merged_tsv else None,
    },
    sub_outputs=[
        # one entry per region BAM (from bam_writer_parallel's stats), transient=True
    ],
    argv=sys.argv,
    stats={
        'n_reads_input': n_reads_in_total,
        'n_reads_corrected': n_reads_out_total,
        'wall_seconds': wall_total,
        # …
    },
    ignore_argv={'--n-threads', '--tmp-dir', '--verbose', '--keep-tmp',
                 '--aligner-concurrency', '--legacy-single-threaded',
                 '--scratch-dir', '--temp-dir', '--output-dir'},
)
write_stage_sidecar(record, sample_output=Path(args.output_dir))
```

The exact `from_components()` signature comes from A.5's `sidecar.py`; verify it matches what's in `rectify/core/provenance/sidecar.py` after pulling A.5.

### 3.1bis NEW `rectify/core/commands/export_merged_tsv_command.py` + parser registration

The back-compat shim subcommand (`rectify export-merged-tsv`). ~50 LOC including parser.

```python
"""rectify export-merged-tsv — concatenate a manifest into a single merged TSV.

For users with downstream scripts that expect a single concatenated
corrected_reads.tsv. Reads a manifest, validates its sha256s, and emits
the merged form. Idempotent: re-running on the same manifest produces
the same bytes.
"""
import argparse
from pathlib import Path

from rectify.core.bam.tsv_partition import load_manifest


def run(args) -> int:
    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"error: manifest not found: {manifest_path}", file=sys.stderr)
        return 1

    entries = load_manifest(manifest_path)
    out_path = Path(args.output or manifest_path.parent / "corrected_reads.tsv")

    # Validate sha256s upfront; fail fast on any mismatch.
    for entry in entries:
        region_tsv = (manifest_path.parent / entry['tsv_path']).resolve()
        if not region_tsv.exists():
            print(f"error: region TSV missing: {region_tsv}", file=sys.stderr)
            return 2
        # sha256 check inline (cheap on TSV files; protects against partial pulls).
        from rectify.core.provenance.hashing import sha256_of_file
        observed = sha256_of_file(region_tsv)
        expected = entry['sha256']
        if observed != expected:
            print(f"error: sha256 mismatch on {region_tsv}: expected {expected[:16]}.., got {observed[:16]}..",
                  file=sys.stderr)
            return 3

    # Concatenate. Header from the first region TSV; skip headers on the rest.
    with out_path.open('w') as out_fh:
        for i, entry in enumerate(entries):
            region_tsv = (manifest_path.parent / entry['tsv_path']).resolve()
            with region_tsv.open() as in_fh:
                first_line = next(in_fh)
                if i == 0:
                    out_fh.write(first_line)  # keep header on first file only
                for line in in_fh:
                    out_fh.write(line)

    print(f"wrote merged TSV: {out_path} ({sum(e['n_rows'] for e in entries)} rows)")
    return 0


def create_export_merged_tsv_parser(subparsers):
    parser = subparsers.add_parser(
        'export-merged-tsv',
        help='Concatenate a corrected_reads.manifest.tsv into a single merged TSV (back-compat shim)',
        description=(
            'Reads a corrected_reads.manifest.tsv (produced by rectify correct in default '
            'manifest-only mode) and emits the legacy concatenated corrected_reads.tsv form '
            'expected by older downstream scripts. Idempotent. Equivalent to `cat` over the '
            'per-region TSVs but with sha256 validation against the manifest.'
        ),
    )
    parser.add_argument('manifest', type=str, help='Path to corrected_reads.manifest.tsv')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output TSV path (default: <manifest_dir>/corrected_reads.tsv)')
    return parser
```

Register in `rectify/cli.py` alongside the other subcommands. Idempotent + sha256-validating per §1.1.

### 3.1ter MODIFY `rectify/core/bam/bam_writer.py:_load_corrections_from_tsv`

This is the load-bearing internal loader. Make it manifest-aware:

```python
def _load_corrections_from_tsv(path: str | Path) -> dict:
    """Load corrections from a corrected_reads.tsv OR a corrected_reads.manifest.tsv.

    The merged-TSV format is the legacy single-file form; the manifest format is
    the new default (Commit B, 2026-05-20). This function auto-detects which
    form the input is by reading the first line: a manifest has the header
        region_id\tchrom\tstart\tend\ttsv_path\tn_rows\tsha256
    while a corrected_reads.tsv has the rectify-correct row header
        read_id\tchrom\t<...>

    Returns: the same {read_id: correction_dict} mapping regardless of input form.
    """
    path = Path(path)
    with path.open() as fh:
        header = fh.readline().rstrip('\n').split('\t')

    if header == ['region_id', 'chrom', 'start', 'end', 'tsv_path', 'n_rows', 'sha256']:
        # Manifest — load per-region TSVs and merge.
        from rectify.core.bam.tsv_partition import load_manifest
        entries = load_manifest(path)
        merged = {}
        for entry in entries:
            region_tsv = (path.parent / entry['tsv_path']).resolve()
            merged.update(_load_corrections_from_single_tsv(region_tsv))
        return merged
    else:
        # Legacy single TSV.
        return _load_corrections_from_single_tsv(path)
```

Refactor: extract the existing parsing into `_load_corrections_from_single_tsv` (same logic, takes a single TSV). The dispatch above wraps it. ALL existing callers of `_load_corrections_from_tsv` continue to work without changes; the manifest case "just works" transparently. No other internal-consumer migration needed.

For consumers that pandas-read the TSV directly (e.g., `scripts/validation_data/check_regression.py`), add a parallel `--manifest` flag if needed, OR document the `rectify export-merged-tsv` shim as the supported path. Don't migrate the scripts in this commit; track separately if needed.

### 3.2 MODIFY `rectify/core/commands/analyze_command.py:95`

Currently `load_corrected_positions(args.input, args.sample_column)` assumes a single concatenated TSV.

**Change:** detect whether `args.input` is a manifest (`*.manifest.tsv` with the header `region_id\tchrom\tstart\tend\ttsv_path\tn_rows\tsha256`) or a single merged TSV.

- If manifest: use the new `load_manifest()` from `rectify.core.bam.tsv_partition` to iterate region TSVs. Concatenate them in-memory for now (Commit D restructures analyze to stream partials; that's not this commit).
- If single TSV: existing behavior, unchanged.

Add `args.input` auto-detection — don't require the user to specify which.

### 3.3 MODIFY `rectify/core/bam/loaders.py`

Add manifest-aware path to `load_corrected_positions()`. Same heuristic as 3.2: detect manifest header, dispatch.

### 3.4 NEW `tests/test_correct_command_parallel_default.py`

Smoke test (NOT in `slow` marker; should run in <30s):

1. Use the bundled validation BAM (~36 reads) — small-input gate WILL fire (sequential path). Verify with a sentinel log capture.
2. Run `correct_command.run(args_simulating_cli)` with `--emit-merged-tsv` (default).
3. Assert outputs: `corrected.bam`, `corrected.bam.bai`, `corrected_reads.manifest.tsv`, `corrected_reads.tsv` ALL exist.
4. Assert `<sample_id>.correct.provenance.json` exists. Load it. Verify schema_version=1.0, exit_status=0, all inputs/outputs have non-empty sha256s, paths are PortablePath envelopes (not bare strings).
5. **Resume test:** rerun immediately. `can_skip_stage` should return SKIP. Total wall < 3 sec.
6. **Resume after input change:** `touch fixture.bam` (mtime change but content unchanged), rerun → should still SKIP (we hash content, not mtime).
7. **Resume after argv change (ignored flag):** rerun with `--n-threads 1` → SKIP.
8. **Resume after argv change (load-bearing flag):** rerun with `--no-emit-merged-tsv` → RUN (the output set changed).

Add a second test that exercises `--legacy-single-threaded`:
- Run the smoke twice — once with default, once with `--legacy-single-threaded`.
- Assert the corrected.bam files are equivalent via `tests.utils.bam_compare.assert_bams_equivalent`.
- Confirms the escape hatch produces identical bytes.

### 3.5 Helper: get the rectify git sha at runtime

Add to `rectify/utils/version.py` (or create if absent):

```python
import os
import subprocess
from pathlib import Path

def get_rectify_git_sha() -> str:
    """Return current HEAD sha of the rectify install's git repo, or 'unknown'.

    Used in PROVENANCE sidecars for skip-check git_sha matching.
    Falls back to 'unknown' if rectify was pip-installed without a .git dir
    (e.g., wheel install) — in that case, sidecar's git_sha mismatch will
    always fire, which is the safe behavior.
    """
    import rectify
    pkg_dir = Path(rectify.__file__).parent.parent
    git_dir = pkg_dir / '.git'
    if not git_dir.exists():
        return 'unknown'
    try:
        return subprocess.check_output(
            ['git', '-C', str(pkg_dir), 'rev-parse', 'HEAD'],
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return 'unknown'
```

Import + use in 3.1 step C and G.

### 3.6 Resume correctness test (CRITICAL)

Add `tests/test_resume_correctness.py`:

Tests that **sidecar emission must follow output durability**. Pseudo-code:

```python
def test_sidecar_not_written_on_crash_after_outputs():
    """If correct_command crashes AFTER outputs are written but BEFORE sidecar,
    next run must RERUN the stage (no SKIP on a half-finished run)."""
    # Setup: run correct once successfully.
    # Then: monkey-patch write_stage_sidecar to raise RuntimeError.
    # Run correct again — it should crash AT sidecar write.
    # Verify: corrected.bam EXISTS but sidecar DOES NOT.
    # Run a third time (real sidecar write restored).
    # Verify: this third run reruns (no SKIP because no sidecar from second).
```

This is the test that protects the §1.2 invariant.

---

## 4. Acceptance criteria

- [ ] Three prerequisites verified (Commits A + A.5 in log, cross-cluster heads match).
- [ ] `rectify/core/commands/correct_command.py` has the dispatch logic, sidecar emission, resume flags.
- [ ] `rectify/core/commands/analyze_command.py` auto-detects manifest vs single-TSV input.
- [ ] `rectify/core/bam/loaders.py` manifest-aware.
- [ ] `rectify/utils/version.py` has `get_rectify_git_sha()`.
- [ ] `--legacy-single-threaded` flag works; output matches default path via `assert_bams_equivalent`.
- [ ] **Default invocation produces `corrected_reads.manifest.tsv` ONLY** (no `corrected_reads.tsv` unless `--emit-merged-tsv` is passed).
- [ ] **`rectify export-merged-tsv <manifest>` shim** produces a TSV byte-identical to the `--emit-merged-tsv`-default output.
- [ ] `_load_corrections_from_tsv` auto-detects manifest vs single-TSV format; both code paths produce identical correction dicts on the same data.
- [ ] **Internal-loader migration verified:** with manifest-only correct output, downstream `restore_polya` and `bedgraph_writers` (which call `_load_corrections_from_tsv`) work unchanged.
- [ ] Small-input gate fires on the 36-read validation BAM; verified via sentinel log.
- [ ] `<sample_id>.correct.provenance.json` and `<sample_id>.restore_polya.provenance.json` are emitted.
- [ ] All 4 path fields in the sidecar are PortablePath envelopes (`kind` is one of sample_relative/env_relative/absolute, NOT bare strings).
- [ ] `can_skip_stage` returns SKIP on immediate rerun (verify the 5 mutation scenarios in §3.4 test 5-8).
- [ ] **Resume correctness test in §3.6 passes** — sidecar emission ordering invariant verified.
- [ ] `pytest -m "not slow"` passes — no regressions.
- [ ] `pytest -m slow` cDNA smoke passes (this commit changes `correct_command` which the cDNA path uses).
- [ ] **Han wt_R1 smoke on H2 (16-core) — this is the STRUCTURAL RESOLUTION TEST for the heap-corruption bug in `rectify/AGENT_FIXES.md` (2026-05-20 "STILL OPEN" section).**

  A second agent on a parallel investigation (note received 2026-05-20) flagged that this commit's region-parallel + per-region-sentinel design IS the structural fix the lab needs. Specifically:
   - The crash happens in `rectify.core.bam.parallel` at the 16-region × 16-worker scale on a monolithic 6.7M-read BAM. Commit A's `write_corrected_bam_parallel` partitions work into per-region BAMs (≤100k ref bp each) with idempotent `.ok` sentinels. This is, at the architectural level, exactly the AGENT_FIXES.md mitigation (`--cpus-per-task=8 --checkpoint-dir`).
   - Without this commit, every future ONT/cDNA dataset producing a >4M-read merged BAM will hit the same wall.

  Run the smoke at FULL Han wt_R1 6.7M-read scale on H2 16-core. Record the outcome in three buckets:

  | Outcome | Interpretation | AGENT_FIXES.md action |
  |---|---|---|
  | A. No crash, completes < 75 min | Bug RESOLVED by region partitioning. | Move "STILL OPEN" → resolved entry, cite commit sha. |
  | B. Crash mid-flight, clean resume from sentinels on second invocation | Bug STRUCTURALLY MITIGATED. Resume is the canonical fix. | Move → resolved via mitigation, document the resume protocol. |
  | C. Crash + resume fails (sentinels insufficient / data corrupted) | Deeper pysam-C bug independent of partitioning. | Keep "STILL OPEN" with sharpened scope; flag the residual failure mode. |

  For outcomes A and B, append a closing entry to `rectify/AGENT_FIXES.md` with the commit sha + scale numbers + wall time. For outcome C, document what specifically still fails (which region's worker, what malformed-read pattern if any) so the next agent has a starting point.

  Default run config (works for either A/B/C outcome):
  ```
  rectify correct han_wt_R1.rectified.bam \
    --short-read --dT-primed-cDNA \
    --genome ... --annotation ... \
    -o results/han_wt_R1/ \
    --write-corrected-bam results/han_wt_R1/wt_R1.corrected.bam \
    -j 16 \
    --tmp-dir $L_SCRATCH/rectify_regions
  # if it crashes mid-flight, RE-RUN the same command; per-region .ok sentinels skip done regions.
  ```

  Wall expectation per Commit Zero: ~6× speedup over baseline → ≤ 75 min if outcome A; ≤ 90 min total (first run + resume) if outcome B.

  Document the chosen scale + outcome in the commit message AND in `rectify/AGENT_FIXES.md`.

  **CLUSTER AVAILABILITY UPDATE (added 2026-05-20 after dispatch — Opus probed both clusters):**

  ```
  H2 pod_smp.q   : 55,760/86,400 slots used + 25,472 jobs WAITING.  DO NOT submit here today.
  H2 campus_prj.q: 2,609 jobs waiting.                              Same — don't submit.
  Sherlock larsms: 2 nodes; both running Kevin's wt_R1_co (25515911), tasks 8-23 PENDING.
  Sherlock owners: 86 nodes (mix of alloc / mix / drain), 7,555 queued.
                   Preemptible; my Commit A `.ok` sentinel-based resume handles preemption cleanly.
  ```

  **Implication:** the queue-submitted smoke is BLOCKED, BUT a workaround exists. Four options for the smoke step:

  1. **DIRECT-SSH-TO-NODE BYPASS (RECOMMENDED).** See the new "Sherlock: direct SSH to active allocation (queue bypass)" section in the lab CLAUDE.md (added 2026-05-20). Kevin already has a wt_R1_co job allocation on `sh03-07n10` and `sh03-08n22`; from the Sherlock login node we can SSH directly to either node and run our chrI-V smoke there, sharing the active allocation. **Use -j 2** (etiquette: don't blow up Kevin's job; nodes are 8-core but load avg 30 means CPU-contested) and **chrI-V subset (~2M reads, below the heap-corruption threshold)**.
     ```bash
     ssh sherlock 'bash --norc --noprofile' << 'EOSSH'
       ssh sh03-07n10 'bash -lc "
         source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh
         conda activate rectify
         cd /oak/stanford/groups/larsms/Users/kevinroy/scratch/commit_b_smoke/  # create this dir first
         # slice the Han wt_R1 BAM to chrI-V (samtools view ... > chrI-V.bam)
         time rectify correct chrI-V.bam \
           --short-read --dT-primed-cDNA \
           --genome ... --annotation ... \
           -o results/ \
           --write-corrected-bam results/wt_R1.chrI-V.corrected.bam \
           -j 2 \
           --tmp-dir /lscratch/rectify_b_smoke
       "'
     EOSSH
     ```
     Expected wall: ~30 min (vs ~15 min on dedicated 8 cores; 2× slower due to CPU contention). All Commit B acceptance items EXCEPT the full-scale heap-corruption structural test are validated.

  2. **DEFER the full-scale heap-corruption resolution test.** With Option 1 covering the architecture-level acceptance, the full 6.7M-read structural test is deferred to whenever the queues clear (or when Kevin's wt_R1_co finishes, freeing larsms). Update AGENT_FIXES.md's "STILL OPEN" entry to note that Commit B's design is hypothesized to fix/mitigate the bug AND that the test is deferred behind cluster availability — explicitly do not claim resolution yet.

  3. **Sherlock owners preemptible submit (fallback only).** Submit the full smoke with `-p owners` and `--requeue`; Commit A's `.ok` sentinels make preemption safe. Risk: 7,555 queued jobs in owners. Don't pursue unless Option 1 fails.

  4. **Piggyback on Kevin's wt_R1_co outputs (advanced).** Wait until 25515911 finishes, inspect outputs, see if Commit B's resume can SKIP-validate them. Requires understanding Kevin's chunked-pipeline output layout. ASK OPUS before pursuing.

  **Recommended sequence:** Option 1 (direct-SSH, chrI-V, -j 2) to validate the architecture immediately; document Option 2 deferral in commit message and AGENT_FIXES.md. The full-scale structural test happens in a follow-up coordinated session.

  If Sonnet hits this section: **proceed with Option 1 directly** — Opus has already cleared the SSH-to-node bypass with Kevin and added it to CLAUDE.md. No further stop-and-confirm needed for the chrI-V subset run.
- [ ] No `git add` / `git commit` runs.

---

## 5. Stop-and-ask triggers

- Either prereq (A, A.5) is not at HEAD.
- `rectify.core.provenance.sidecar.ProvenanceRecord.from_components()` signature doesn't match what your wiring expects — read it and surface the mismatch; don't bend your wiring around an unexpected API.
- The Han wt_R1 smoke crashes or exceeds 90 min wall.
- `assert_bams_equivalent` reports a divergence between parallel and `--legacy-single-threaded` paths.
- Resume returns SKIP after a real input change (sha256 false-positive) — STOP, this is a correctness bug.
- A pre-existing test regresses for a reason you can't pin to your changes — STOP, investigate, don't paper over.

---

## 6. Reporting back

1. `git status -s` output.
2. `git diff --stat` summary.
3. `pytest -m "not slow"` summary.
4. `pytest -m slow` cDNA smoke result.
5. Han wt_R1 smoke result (wall time, output BAM equivalence vs HEAD's reference output, sidecar emitted with correct fields).
6. Resume test results: all 8 acceptance scenarios from §3.4 + §3.6 pass.
7. Any TODOs you left.

---

## 7. Time budget

- Read briefings + verify prereqs + grep current file states: ~30 min.
- `correct_command.py` wiring + dispatch + sidecar emission: ~90-120 min (this is the heavy work; first sidecar wiring sets the template for D/E/F so be careful).
- `analyze_command.py` + `loaders.py` manifest auto-detect: ~30 min.
- `version.py` helper: ~10 min.
- Tests (test_correct_command_parallel_default + test_resume_correctness): ~60-90 min.
- Han wt_R1 smoke on H2 (run + observe): ~75 min wall, ~10 min agent time.
- pytest + iterate: ~30 min.
- **Total: ~5-7 hours.** Plan accordingly.

**End of briefing.**
