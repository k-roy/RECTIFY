# Handoff: cDNA UMI-binning — Phase D ran, XC-tag bug, Option 1 chosen

**Created:** 2026-05-18 (mid-session)
**Last updated:** 2026-05-18 (post-Phase-D, polished for new session)
**Status:** Phase A + B landed in the working tree. Phase D Sherlock job 25308818
completed at 12:27 PT but per-bin tables came out empty (XC tag dropped at FASTQ
extraction — see §5). **User has chosen Option 1: sidecar BAM lookup.** The
next session should pick up at §3 "Start here".

Longer-form project plan is in
`docs/handoffs/cdna_umi_stratified_calibration.md`; this doc is the live-state
handoff.

---

## 1. Start here (cold-session checklist)

```bash
cd "/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/rectify"

# 1a. Verify Phase A+B working-tree state is intact.
python -m pytest tests/test_data_bundling.py tests/test_profiler_umi_binning.py --quiet
# Expect: 29 passed

# 1b. Verify the 8-chunk pooled cDNA table is bundled.
wc -l rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv
# Expect: 126 lines (125 data rows + header)
head -3 rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv

# 1c. Confirm the per-bin TSVs are NOT yet present (bug from Phase D).
ls rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna*.tsv
# Expect: only penalty_scores_cdna.tsv

# 1d. Confirm Sherlock ControlMaster is alive; if not, Kevin re-establishes it.
ssh sherlock 'hostname'
```

Then jump to **§4 — Option 1 work plan**.

---

## 2. What's done (Phase A + B in working tree)

### Phase A — profiler XC-aware binning

`scripts/calibration/empirical_cigar_error_profiler.py`:
- `--umi-tag TAG` (default `XC`), `--umi-bins SPEC` (default `1,2,3+`)
- `parse_umi_bins()`, `assign_bin()`, `_sanitize_bin_label()`
- `profile_chunk()` extended with `umi_tag`, `umi_bins`, `counts_by_bin`,
  `union_counts_by_bin`, `str_counts_by_bin`, `reads_per_bin` (all kwarg,
  default `None` → back-compat). Per-read bin lookup done once after
  chromosome agreement check; route increments at every pooled site.
- `_worker()` builds + returns per-bin defaultdicts; `_merge_counts_by_bin()`.
- `main()` writes one `error_counts_umi<bin>.tsv` per bin (suffix form),
  empty bins emit header-only TSVs, one-line summary per bin (label,
  n_reads, n_events).

**Untagged reads** route to pooled-only (no per-bin output). **First-aligner-hit**
is used for the XC lookup (not all-aligners-agree) since aux tags typically
propagate from one aligner's BAM only — this assumption breaks for cDNA,
see §5.

Test: `tests/test_profiler_umi_binning.py` — 10 tests, all passing.

### Phase B — bundler + resolver bin-aware

`scripts/calibration/bundle_protocol_tables.py`:
- `--umi-bin LABEL` CLI flag. Default behavior unchanged.
- With `--umi-bin`: reads `error_counts_<LABEL>.tsv` (NOT pooled), writes
  `penalty_scores_<protocol>_<LABEL>.tsv`, skips STR (Phase A keeps STR
  pooled-only).
- `_read_counts(hp_filename=...)` parametrized so the per-bin variant flows
  through without touching the STR path.

`rectify/data/__init__.py`:
- `BUNDLED_GENOMES['saccharomyces_cerevisiae']['protocols']['cdna']`
  restructured from flat to nested `bin → key → path`:
  ```python
  'cdna': {
      'pooled':    {'junction_penalty_table': '...penalty_scores_cdna.tsv',
                    'str_penalty_table':      '...str_penalty_scores_cdna.tsv',
                    'junction_overhang_table':'...junction_overhang_table_cdna.tsv'},
      'umi1':      {'junction_penalty_table': '...penalty_scores_cdna_umi1.tsv'},
      'umi2':      {'junction_penalty_table': '...penalty_scores_cdna_umi2.tsv'},
      'umi3plus':  {'junction_penalty_table': '...penalty_scores_cdna_umi3plus.tsv'},
  }
  ```
  QSrev's protocols entry stays FLAT (resolver shape-detects).
- `_resolve_bundled_table(organism, key, protocol=None, bin=None)`.
  Resolution: (1) per-bin nested entry on disk → (2) `pooled` sub-key →
  (3) flat shape → (4) org-level flat (DRS) → (5) None.
- `get_bundled_junction_penalty_table` / `_str_` / `_junction_overhang_`
  all accept `bin=None`. Default preserves back-compat.
- `resolve_reference_paths()` NOT changed (Phase C territory).

Test: `tests/test_data_bundling.py` extended with 4 new cases (one
parametrized × 3 bins = 6 items). **29 tests in test_data_bundling +
test_profiler_umi_binning, all passing.**

### Other uncommitted working-tree state (from yesterday)

The empirical-error tables session
(`docs/handoffs/empirical_error_tables_cdna_qsrev.md`) left these:
- `scripts/calibration/empirical_cigar_error_profiler.py` — `=` SEQ +
  QNAME-whitespace fixes
- `scripts/calibration/run_profiler_qsrev_h2.sh` — `module load bwa` fix
- `scripts/calibration/run_profiler_{qsrev,cdna}_sherlock.sh` — Sherlock
  SLURM wrappers
- `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv`
  — now overwritten by Phase D's 8-chunk pool (126 lines)
- `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_qsrev.tsv`
  — 579 M-event QSrev table
- `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/README.md`
  — provenance for both

---

## 3. What's verified

- `pytest tests/test_data_bundling.py tests/test_profiler_umi_binning.py -v`
  → **29 passed, 0 failed** as of 2026-05-18 ~08:25 (re-verified at
  ~13:00 after Phase D — still green).
- `python scripts/calibration/bundle_protocol_tables.py --help` shows the
  new `--umi-bin` flag.
- Resolver verified to return:
  - DRS table for `get_bundled_*(ORG)` (no protocol)
  - cDNA pooled for `protocol='cdna'` (no bin)
  - cDNA pooled FALLBACK for `protocol='cdna', bin='umiN'` (since per-bin
    TSVs don't exist yet — once Option 1 lands, this fallback turns into
    a direct hit)
  - QSrev for `protocol='qsrev'` (flat-shape detection)
- **Pooled 8-chunk cDNA table is bundled** at
  `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv`
  (125 rows, 477 M strict-agreed events, overwrote yesterday's 117-row
  2-chunk version). `[uncommitted]`

---

## 4. Option 1 work plan — sidecar BAM lookup (~60 min code + ~1 h Sherlock)

The Phase D failure is upstream of the profiler: `samtools fastq` in
`run_profiler_cdna_sherlock.sh` doesn't pass `-T XC`, so XC drops at FASTQ
extraction and never reaches the aligned BAMs. The profiler's
`read.get_tag('XC')` then raises `KeyError` for every read → all reads
route to the untagged-only path → per-bin counts stay at 0.

The aligned BAMs are still valid CIGAR-wise — only the aux tag is missing.
So we patch the profiler to consult a **sidecar BAM** (the
`stage1_consensus.bam` that DOES carry XC), and re-run the profile step
only. Alignment cache stays warm.

### 4.1 Profiler change

Edit `scripts/calibration/empirical_cigar_error_profiler.py`:

1. Add CLI: `--umi-source-bam PATH` (default `None`).
2. At startup (in `main()`, before worker dispatch), if
   `--umi-source-bam` is set:
   - Open it once with `pysam.AlignmentFile(path, 'rb',
     check_sq=False)`.
   - Build `umi_lookup: dict[str, int|str] = {}` — `read.query_name →
     read.get_tag(umi_tag)`. Skip reads lacking the tag.
   - Pass `umi_lookup` into `_worker()` → `profile_chunk()`.
3. In `profile_chunk()`, when resolving the per-read bin:
   - If `umi_lookup is not None` and `read.query_name in umi_lookup`,
     use that value.
   - Else fall back to the existing `read.get_tag(umi_tag)` path
     (which still works for QSrev where the tag IS on the aligned BAM).
   - Untagged path unchanged.

### 4.2 Profiler test

Extend `tests/test_profiler_umi_binning.py` with one new test
`test_sidecar_bam_lookup_routes_per_bin`:
- Build a tiny synthetic sidecar BAM with three reads (XC=1, XC=2, XC=5)
  using `pysam.AlignmentFile` in-memory or a tmp_path.
- Build a tiny per-aligner BAM with the SAME three read names but NO XC
  tag.
- Invoke `profile_chunk(..., umi_tag='XC', umi_bins=parse_umi_bins('1,2,3+'),
  counts_by_bin=..., umi_lookup={...})` and assert each read's events
  land in the right bin TSV.

Run: `pytest tests/test_profiler_umi_binning.py -v` — expect 11 passed.

### 4.3 Wrapper change

Edit `scripts/calibration/run_profiler_cdna_sherlock.sh`:

- Add `--umi-source-bam "$WORK/raw/stage1_consensus.bam"` to every
  per-chunk profile invocation AND the pooled profile invocation.
- Do NOT touch the `samtools fastq` line or the alignment loop —
  skip-if-exists must hold.

### 4.4 Re-run on Sherlock

```bash
# Stage updated files to Sherlock.
rsync -a scripts/calibration/empirical_cigar_error_profiler.py \
        scripts/calibration/run_profiler_cdna_sherlock.sh \
    sherlock:/scratch/users/kevinroy/calibration_2026_05/cdna/scripts/

# Sanity-check the staged copies before sbatching.
ssh sherlock 'grep -c -- "--umi-source-bam" /scratch/users/kevinroy/calibration_2026_05/cdna/scripts/empirical_cigar_error_profiler.py'
# Expect: >= 1
ssh sherlock 'grep -c -- "--umi-source-bam" /scratch/users/kevinroy/calibration_2026_05/cdna/scripts/run_profiler_cdna_sherlock.sh'
# Expect: 9 (8 per-chunk + 1 pooled)  -- adjust if the wrapper batches differently

# Submit. Skip-if-exists will jump straight to profiling (~1 h wall).
ssh sherlock 'cd /scratch/users/kevinroy/calibration_2026_05/cdna && \
    export PS1=""; source ~/.bashrc; \
    sbatch scripts/run_profiler_cdna_sherlock.sh'
# Returns: Submitted batch job <JID>
```

Background-monitor with a Bash `run_in_background` polling
`squeue -j <JID>` until empty. (One job → one notification on completion.
Do NOT poll-loop in the foreground.)

### 4.5 Post-completion: rsync, bundle, flip test, verify

```bash
ssh sherlock 'sacct -j <JID> --format=JobID,State,Elapsed,ExitCode | head -5'
# Expect: COMPLETED 0:0
```

If COMPLETED, follow the original §4.1 → 4.5 path (preserved below as
§4.6 → §4.10), but verify FIRST that per-bin TSVs are now NON-empty:

```bash
ssh sherlock 'head -3 /scratch/users/kevinroy/calibration_2026_05/cdna/error_profile_cdna_20260518/pooled/error_counts_umi1.tsv'
# Expect: header + at least one data row
```

If the per-bin TSVs are STILL empty, the sidecar lookup didn't fire —
re-check that `--umi-source-bam` reached the profiler invocation
(see §4.4 sanity check) and that the sidecar BAM has XC tags
(`ssh sherlock 'samtools view /scratch/users/kevinroy/calibration_2026_05/cdna/raw/stage1_consensus.bam | head -3 | grep -c XC:'`).

### 4.6 Rsync results to M1

```bash
cd "/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/My Drive/Work/Chanfreau Lab/rectify"
rsync -a sherlock:/scratch/users/kevinroy/calibration_2026_05/cdna/error_profile_cdna_20260518 \
    scripts/calibration/runs/
ls scripts/calibration/runs/error_profile_cdna_20260518/pooled/
# expect: error_counts.tsv  error_counts_umi1.tsv  error_counts_umi2.tsv  error_counts_umi3plus.tsv  (all non-empty)
```

### 4.7 Bundle per-bin + re-bundle pooled

```bash
RUN=scripts/calibration/runs/error_profile_cdna_20260518/pooled

for bin in umi1 umi2 umi3plus; do
    python scripts/calibration/bundle_protocol_tables.py \
        --protocol cdna --umi-bin $bin \
        --inputs "$RUN"
done

# Re-bundle pooled (no-op if file already current, but cheap to redo).
python scripts/calibration/bundle_protocol_tables.py \
    --protocol cdna \
    --inputs "$RUN"

ls rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna*.tsv
# expect: penalty_scores_cdna.tsv  penalty_scores_cdna_umi1.tsv  penalty_scores_cdna_umi2.tsv  penalty_scores_cdna_umi3plus.tsv
```

### 4.8 Flip the fallback test

`tests/test_data_bundling.py` has a
`test_cdna_per_bin_falls_back_to_pooled_when_file_absent` test that
asserts the per-bin TSVs DON'T exist. Now they do. Replace with:

```python
@pytest.mark.parametrize('umi_bin', ['umi1', 'umi2', 'umi3plus'])
def test_cdna_per_bin_resolves_to_specific_file(umi_bin):
    bin_specific = (d.BUNDLED_DATA_DIR /
                    f'genomes/{ORG}/penalty_tables/penalty_scores_cdna_{umi_bin}.tsv')
    assert bin_specific.exists(), f'{bin_specific} missing'
    p = d.get_bundled_junction_penalty_table(ORG, protocol='cdna', bin=umi_bin)
    assert p is not None
    assert p.name == f'penalty_scores_cdna_{umi_bin}.tsv'
```

Then: `pytest tests/test_data_bundling.py -v` should still be all green.

### 4.9 Convergence check across chunks per bin

Drop this into `/tmp/convergence.py` and run from the repo root:

```python
import sys, statistics
from pathlib import Path
sys.path.insert(0, 'scripts/calibration')
from empirical_cigar_error_profiler import compute_rates
import pandas as pd

RUN = Path('scripts/calibration/runs/error_profile_cdna_20260518')
rows = []
for bin_label in ['umi1', 'umi2', 'umi3plus']:
    per_chunk_rates = {}
    per_chunk_counts = {}
    for c in range(8):
        tsv = RUN / f'chunk_{c}' / f'error_counts_{bin_label}.tsv'
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep='\t')
        counts = {(r.op_type, r.ref_base, int(r.hp_length)): int(r['count']) for _, r in df.iterrows()}
        rate_df = compute_rates(counts)
        for _, r in rate_df.iterrows():
            cell = (r.op_type, r.ref_base, int(r.hp_length))
            per_chunk_rates.setdefault(cell, []).append(r.rate)
            per_chunk_counts.setdefault(cell, []).append(int(r.count))
    per_cell_maxabs = []
    for cell, rates in per_chunk_rates.items():
        if max(per_chunk_counts[cell]) < 100:
            continue
        med = statistics.median(rates)
        if med == 0:
            continue
        max_abs_pct = max(abs(100 * (r - med) / med) for r in rates)
        per_cell_maxabs.append(max_abs_pct)
    if per_cell_maxabs:
        rows.append((bin_label, len(per_cell_maxabs),
                     statistics.median(per_cell_maxabs),
                     max(per_cell_maxabs)))
    else:
        rows.append((bin_label, 0, float('nan'), float('nan')))

print(f"| Bin       | n_cells | median |Δ%| | max |Δ%| | Verdict |")
print(f"|-----------|--------:|------------:|---------:|---------|")
for label, n, med, mx in rows:
    verdict = "8 chunks sufficient" if med < 10 else "extend to 16 (>10% med)"
    print(f"| {label:9s} | {n:7d} | {med:10.2f}% | {mx:7.2f}% | {verdict} |")
```

Save the printed table as
`scripts/calibration/runs/error_profile_cdna_20260518/CONVERGENCE.md`.

### 4.10 Final sanity check

```bash
pytest tests/test_data_bundling.py tests/test_profiler_umi_binning.py -v
```

Should be 30 passing now (29 + the new sidecar test from §4.2), with the
test edit from §4.8 in place.

---

## 5. Phase D diagnostic record (preserved for context)

### Sherlock job 25308818 — what ran

- **State:** COMPLETED 0:0
- **Wall:** 3 h 59 min (08:28 → 12:27 PT 2026-05-18, node sh03-08n22)
- **Stages:** alignment of chunks 2-7 ~3 h, profile (per-chunk + pooled)
  ~1 h. Chunks 0+1 cached from yesterday — skip-if-exists fired.

### Pooled 8-chunk table — GOOD

- 477,123,752 strict agreed events pooled across all 8 chunks (vs 119 M
  from yesterday's 2-chunk pool — strict 4× improvement on event counts).
- 7,427,483 STR events.
- Reference rates: `sub(HP=1) AT=0.0019 CG=0.0041`,
  `ins(HP=1) AT=0.0010 CG=0.0012` (consistent with yesterday —
  convergence on the pooled estimator was already strong).
- Pooled table bundled at
  `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv`
  (125 data rows). `[uncommitted]`

### Per-bin tables — empty (root cause)

```
UMI-bin summary  (--umi-tag XC):
  umi1        n_reads=       0  n_events=         0
  umi2        n_reads=       0  n_events=         0
  umi3plus    n_reads=       0  n_events=         0
```

`samtools fastq` in `run_profiler_cdna_sherlock.sh` line 54 doesn't pass
`-T XC`, so the XC tag from `stage1_consensus.bam` is dropped at FASTQ
extraction. Aligners then have no source for the tag, so per-aligner BAMs
carry only native aux (NM, MD, AS, ts, tp, ...). Profiler's
`read.get_tag('XC')` raises `KeyError` for every read → all reads routed
to the untagged-only path.

Diagnosis verified on `chunks/wt_rep2_consensus_chunk_000_of_016.fastq.gz`:
- FASTQ header: `@3988451a-7bee-49e7-8584-264971d6ae06` (no comment).
- Per-aligner BAM aux: `NM:i:15 AS:i:1200` only.

Source `stage1_consensus.bam` first record DOES carry `XC:i:1` — confirming
the tag exists upstream and is lost only at the FASTQ-extraction step.

### Salvage value of the failed run

- All per-aligner BAMs at
  `/scratch/users/kevinroy/calibration_2026_05/cdna/chunks/align/{minimap2,mapPacBio,gapmm2,uLTRA,deSALT}/`
  are valid CIGAR-wise. Option 1 reuses 100% of the ~3 h alignment compute.
- The M1 rsync at `scripts/calibration/runs/error_profile_cdna_20260518/`
  has all per-chunk + pooled `error_counts.tsv` + the empty
  `error_counts_umi{1,2,3plus}.tsv` files. Pooled has been bundled.

### Fallback recovery options (if Option 1 fails)

**Option 2 — Re-extract FASTQ with `-T XC` and re-align (~4 h)**
- Modify `run_profiler_cdna_sherlock.sh` line 54:
  `samtools fastq -@ 4 -T XC "$CONS_BAM" | gzip -1 > "$FQ"`.
- Delete `chunks/*.fastq.gz` and `chunks/align/*` to bust skip-if-exists.
- Verify each aligner propagates FASTQ comments to BAM aux:
  - minimap2 needs `-y` and TAB-separated comment tags. `samtools fastq -T`
    writes SPACE-separated by default — see lab memory
    `feedback_minimap2_y_tab_sep.md`. Will need a sed/awk step to convert
    separator.
  - mapPacBio (bbmap variant), gapmm2, deSALT, uLTRA: unknown — needs
    per-aligner verification.
- Re-run full pipeline. ~4 h Sherlock wall. More architecturally pure but
  expensive and has unknown aligner-by-aligner propagation issues.

**Option 3 — Phase C without per-bin support**
- Accept the consensus-weighted-average cDNA table as canonical. `rectify
  correct` loads one table and uses it for every read regardless of UMI
  cluster size. Weakest semantics (over-confident for raw-cDNA users),
  least work.

---

## 6. Out of scope for the next session

- **Phase C** (rectify correct integration — per-read table selection in
  `junction_refiner.py`) is deferred. Surgical change to load-bearing
  Module 2H; deserves interactive review, not agent dispatch. See
  `docs/handoffs/cdna_umi_stratified_calibration.md` §3.3 Phase C.
- **Phase E** (README updates + integration smoke test) waits until
  Phase C lands.

---

## 7. Files (current state at handoff write)

**Phase A + B working tree (uncommitted):**
- `[uncommitted]` `scripts/calibration/empirical_cigar_error_profiler.py`
- `[uncommitted]` `scripts/calibration/bundle_protocol_tables.py`
- `[uncommitted]` `scripts/calibration/run_profiler_cdna_sherlock.sh`
  (CHUNKS 0-7 + `--umi-tag XC --umi-bins 1,2,3+`)
- `[uncommitted]` `rectify/data/__init__.py`
- `[uncommitted]` `tests/test_data_bundling.py` (4 new tests; one to
  flip in §4.8)
- `[uncommitted, new]` `tests/test_profiler_umi_binning.py`

**Phase D-produced bundled tables:**
- `[uncommitted]` `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv`
  — 8-chunk pool, 125 rows. Already written; Option 1 will overwrite with
  re-bundle in §4.7 (no-op-equivalent unless inputs change).
- `[pending]` `penalty_scores_cdna_umi1.tsv` — Option 1 §4.7
- `[pending]` `penalty_scores_cdna_umi2.tsv` — Option 1 §4.7
- `[pending]` `penalty_scores_cdna_umi3plus.tsv` — Option 1 §4.7

**Sherlock staging (durable across sessions):**
- `/scratch/users/kevinroy/calibration_2026_05/cdna/raw/stage1_consensus.bam`
  — sidecar BAM for Option 1 (920 MB, XC-tagged).
- `/scratch/users/kevinroy/calibration_2026_05/cdna/chunks/align/{minimap2,mapPacBio,gapmm2,uLTRA,deSALT}/`
  — cached per-aligner BAMs, valid CIGAR, reused by Option 1.
- `/scratch/users/kevinroy/calibration_2026_05/cdna/scripts/` — staging
  area for the patched profiler + wrapper.
- `/scratch/users/kevinroy/calibration_2026_05/cdna/slurm_logs/` — prior
  job logs (25308818.out / .err) for forensics.
