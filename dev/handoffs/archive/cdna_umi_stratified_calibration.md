# Handoff: UMI-bin-stratified cDNA penalty calibration

**Created:** 2026-05-17 (post `empirical_error_tables_cdna_qsrev.md` session)
**Driver:** Kevin Roy
**Status:** Design + open work. The current single bundled cDNA penalty
table is a useful first cut but is **biased optimistic for UMI singletons**
(the largest cohort). Stratifying by UMI cluster size is the next step.

---

## 1. What was done (foundation this builds on)

| Commit | Summary |
|---|---|
| `ba602c3` | `feat(data): bundle yeast penalty tables; protocol-aware --Scer auto-resolve` — `BUNDLED_GENOMES['saccharomyces_cerevisiae']['protocols']` map, `get_bundled_*_table(organism, protocol=...)` accessors, `_protocol_for_args()` helper, `resolve_reference_paths()` routes via `args.dT_primed_cDNA` / `args.short_read`. 13 tests in `test_data_bundling.py`. |
| `39458d3` | `chore(scripts+docs): validation review tooling, calibration plumbing, TODO` — `scripts/calibration/empirical_cigar_error_profiler.py` with `--protocol {drs,cdna,qsrev}`, `--exclude-5prime-bp N`, `--exclude-3prime-bp N`; `get_exonic_ranges_inner()`; SGE/SLURM submit scripts; `bundle_protocol_tables.py`. |
| uncommitted | `=` SEQ + QNAME-whitespace bugfixes (profiler) + `module load bwa` (H2 script) + README updates + bundled cDNA/QSrev TSVs. See §5. |

**First-pass calibration outputs (live on disk):**
- `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv` — 119 M HP events pooled across 2 of 16 chunks of `stage1_consensus.bam` (Sherlock JID 25269616). 2-chunk-vs-full convergence checked: median |Δ rate%| = 2.17% on count_total≥100 cells.
- `penalty_scores_qsrev.tsv` — 579 M HP events from Han 2023 WT (wt_R1 + wt_R2), H2 JID 13421248.

---

## 2. What's verified

- `pytest tests/test_data_bundling.py tests/test_validation_reads_cdna.py tests/test_quantseq_rev_walkback.py tests/test_junction_refiner.py` → 67 passed / 17 skipped. Bundling resolvers correctly route `--Scer --dT-primed-cDNA` → `penalty_scores_cdna.tsv` and `--Scer --short-read --dT-primed-cDNA` → `penalty_scores_qsrev.tsv`. STR + overhang fall back to DRS as designed.
- DRS bundled BAMs (`/oak/.../dev_runs/wt_by4742_rep1_drs_trim_20260417/`) verified clean of the `=` SEQ encoding bug — DRS table is **not** affected and does not need regeneration.
- **NOT VERIFIED:** the cDNA table's behavior on raw (non-UMI-consensus) cDNA — it is calibrated only on consensus reads, so singleton reads will see optimistic penalties. This is the gap this handoff is for.

---

## 3. Open items (the UMI stratification work)

The user observation (2026-05-17):

> "In our initial datasets, roughly 75% of reads are UMI singletons, and
> there is an exponential decay pattern for the remaining 25%. Might be
> worth incorporating the # of UMIs for each read into the error profiler."

### 3.1 Why this matters

The current `penalty_scores_cdna.tsv` averages reads of different consensus
quality:

| UMI cluster size | Expected accuracy | % of stage1 reads (Kevin's estimate) |
|---|---|---|
| 1 (singleton)    | raw cDNA, ~95–96% | ~75 % |
| 2                | doubleton consensus, ~98 % | ~15 % |
| ≥3               | high-confidence consensus, ≥99 % | ~10 % |

A single penalty table over a 75/15/10 mix is dominated by the singleton
error rate but with the doubleton/triplet tail pulling rates DOWN. A
RECTIFY run on raw cDNA against this table will trust junction
refinements more than the data supports.

### 3.2 Design questions to answer before coding

These need an explicit decision before starting Phase A below:

1. **Bin granularity.** Three bins (1 / 2 / ≥3) is the obvious split given
   the 75/15/10 distribution. A finer split (1 / 2 / 3 / 4 / ≥5) buys
   resolution at the cost of low-count cells. Recommend **3 bins** unless
   a histogram of cluster sizes in `stage1_consensus.bam` shows a wider
   tail than expected.
2. **Tag carrying UMI count.** **`XC:i:N`** (verified 2026-05-17 against
   `/u/scratch/k/kevinroy/ont_cdna_v4/wt_rep2/out/stage1_consensus.bam`).
   First read has `XC:i:1` (singleton), consistent with the 75% singleton
   estimate. Use `--umi-tag XC` in the profiler.
3. **Fallback strategy.** If a read has no UMI tag (e.g. someone runs
   `rectify correct` on raw cDNA without going through the UMI pipeline),
   should it default to the singleton table (conservative) or the
   pooled table (matching today's behavior)? Recommend **singleton table
   with a one-time warning log**.
4. **Sparse-cell handling.** The ≥3 bin will have only ~10 % of events
   (~12 M of the 119 M pooled today). High-HP cells will be low-count.
   Re-use `min_count_flag` semantics from `derive_penalties()` — flag rows
   with `count_total < 50` and let the runtime fall through to a wider
   bin.

### 3.3 Work plan

**Phase A — profiler stratification.** Estimate: 2–3 h.

- In `scripts/calibration/empirical_cigar_error_profiler.py`:
  - Add `--umi-tag TAG` (default: confirmed tag from §3.2 Q2).
  - Add `--umi-bins SPEC` (default: `1,2,3+`). Parses to a list of bin edges.
  - Extend the per-read `counts` dict to a dict-of-dicts keyed by bin label:
    `counts[bin_label][(op, base, hp)] += 1`. The existing single dict
    becomes `counts['pooled']` (kept for back-compat).
  - In `profile_chunk()`, look up `read.get_tag(umi_tag)` for each read,
    map to a bin label, then route counts.
  - Output `error_counts_<bin>.tsv` per bin, plus the pooled file.
- Add a sentinel test: profile against the bundled cDNA validation BAM
  with a synthetic `XK` tag injected — confirm per-bin tables roll up
  to the pooled table when summed.

**Phase B — bundler stratification.** Estimate: 1 h.

- `scripts/calibration/bundle_protocol_tables.py`: accept `--umi-bin BIN`
  argument. Writes the per-bin file as
  `penalty_scores_cdna_umi<bin>.tsv` (e.g. `_umi1`, `_umi2`, `_umi3plus`).
- Extend `BUNDLED_GENOMES[...]['protocols']['cdna']` from a single map to a
  nested structure:
  ```python
  'cdna': {
      'umi1':    {'junction_penalty_table': 'penalty_tables/penalty_scores_cdna_umi1.tsv', ...},
      'umi2':    {'junction_penalty_table': 'penalty_tables/penalty_scores_cdna_umi2.tsv', ...},
      'umi3plus':{'junction_penalty_table': 'penalty_tables/penalty_scores_cdna_umi3plus.tsv', ...},
      'pooled':  {'junction_penalty_table': 'penalty_tables/penalty_scores_cdna.tsv', ...},
  }
  ```
- Today's `_resolve_bundled_table(organism, key, protocol=None)` API
  needs another arg: `bin=None`. Default `bin=None` returns the `pooled`
  entry — preserves today's behavior for any caller that doesn't know
  about bins yet.

**Phase C — `rectify correct` integration.** Estimate: 4–6 h. **This is the
load-bearing work.**

- `rectify correct` currently loads ONE junction penalty table at startup.
  Needs to load all 3 (or N) bin-specific tables.
- Per-read at correction time: read the UMI-count tag, pick the matching
  table. `junction_refiner.py::_score_junction()` is the call site —
  threading a per-read penalty table into it without breaking the
  existing API will require a small refactor (likely a
  `PenaltyTableSet` class holding all bins).
- Smoke test: synthesise a BAM with reads carrying `XC:i:1`, `XC:i:2`,
  `XC:i:5`. Confirm `rectify correct` picks the right table per read
  (use `--verbose` and grep the log for the "loaded penalty table" line
  to verify).
  (Note: `XC:i:N` is the UMI cluster size; `XK:i:N` is the 3′ trim
  length — the original spec had a typo here.)

**Phase D — recalibrate on H2 / Sherlock.** Estimate: ~1 h compute + ~30
min orchestration. The Sherlock 2-chunk cDNA panel (`/scratch/users/kevinroy/calibration_2026_05/cdna/chunks/`)
already has all 5 aligners' BAMs cached for chunks 0 and 1. Re-run the
profiler with `--umi-tag XK --umi-bins 1,2,3+` against that panel — the
alignment phase is skipped entirely. Pool both chunks, bundle each bin,
re-verify convergence with chunk_0-vs-chunk_1 |Δ%| (HP=1..4 high-count
cells should agree within 5% per bin).

**Phase E — docs + tests.** Estimate: 1 h.

- Update `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/README.md`
  with the three new files + the routing logic.
- Update `tests/test_data_bundling.py` with bin-routing tests
  (synthetic args, expected resolved path per bin).
- Update `docs/CLAUDE.md` "Empirical Penalty Tables" if Kevin grants
  permission (today's session was auto-blocked from editing it).

### 3.4 Why this was deferred (not done in the originating session)

- Context window of the session was already long after a debugging
  marathon (8 unrelated bugs surfaced and fixed in one stretch).
- The Phase C refactor of `junction_refiner.py` is a surgical change to
  Module 2H (the load-bearing correction step) and deserves a fresh
  head, not the tail end of a session.
- The current single-table bundle ships a useful average — not strictly
  wrong, just over-confident on singletons. There's no production
  emergency.

---

## 4. Resume

**Check:** open `stage1_consensus.bam` on H2 at
`/u/scratch/k/kevinroy/ont_cdna_v4/wt_rep2/out/stage1_consensus.bam` and
inspect the aux tags on the first read with
`samtools view <bam> | head -1 | tr '\t' '\n' | grep -E '^X[A-Z]:'` —
record which `X*:i:N` tag carries the UMI cluster size. Then read
`memory/project_cdna_pipeline.md` and `memory/project_rectify_tag_namespace.md`
to confirm the tag's semantics.

Once the tag is identified, start at **Phase A §3.3** above. The Sherlock
cached panel at `/scratch/users/kevinroy/calibration_2026_05/cdna/chunks/`
is ready — you can land a per-bin calibration in ~3–4 h without burning
new alignment compute.

If the tag is NOT what's documented in memory (or the count distribution
is very different from 75/15/10), STOP and re-check §3.2 design questions
1 and 2 with Kevin before implementing.

---

## 5. Files touched

**Committed (today's session, `ba602c3` + `39458d3`):**
- `rectify/data/__init__.py` — `BUNDLED_GENOMES['protocols']` map, accessors, `_protocol_for_args`
- `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/{penalty_scores,str_penalty_scores,junction_overhang_table}.tsv` — bundled DRS (unchanged content, just bundled location)
- `scripts/calibration/empirical_cigar_error_profiler.py` — `--protocol` / `--exclude-{5,3}prime-bp` flags, `get_exonic_ranges_inner`
- `scripts/calibration/bundle_protocol_tables.py` — first version
- `scripts/calibration/run_profiler_{qsrev,cdna}_h2.sh` — H2 SGE wrappers
- `tests/test_data_bundling.py` — 13 tests
- `docs/handoffs/empirical_error_tables_cdna_qsrev.md` — the prior handoff (§9 progress addendum)

**Uncommitted [working tree]:**
- `[uncommitted]` `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_cdna.tsv` — first-pass cDNA table (pooled, no UMI binning yet)
- `[uncommitted]` `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/penalty_scores_qsrev.tsv` — QSrev table
- `[uncommitted]` `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/README.md` — provenance for both new tables, convergence-check rationale, `=` SEQ note for DRS clean verification
- `[uncommitted]` `scripts/calibration/empirical_cigar_error_profiler.py` — `=` SEQ + QNAME-whitespace bugfixes (these are MANDATORY for the next session; don't revert)
- `[uncommitted]` `scripts/calibration/run_profiler_qsrev_h2.sh` — `module load bwa` fix
- `[uncommitted]` `scripts/calibration/run_profiler_{qsrev,cdna}_sherlock.sh` — Sherlock SLURM wrappers (work as-is)

**H2 + Sherlock staging (will outlive the session):**
- H2: `/u/scratch/k/kevinroy/calibration_2026_05/{cdna,qsrev}/chunks/align/` — 16-chunk cDNA partial (only ~4 chunks aligned before kill) + complete 8-chunk QSrev panel.
- Sherlock: `/scratch/users/kevinroy/calibration_2026_05/cdna/chunks/align/` — 2-chunk cDNA panel, all 5 aligners. **This is the panel to re-profile for the UMI binning work — no realignment needed.**
- Sherlock: `/scratch/users/kevinroy/calibration_2026_05/qsrev/chunks/align/` — 8-chunk QSrev panel (bbmap+bwa). Useful if QSrev needs a re-profile but not required for the cDNA work.
