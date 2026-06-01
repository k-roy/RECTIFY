# Empirical Nanopore Error Penalty Table

**Current production run:** `error_profile_20260422/`  
**Dataset:** `wt_by4742_rep1` — S. cerevisiae wild-type, R10.4.1 chemistry, Dorado basecaller  
**Method:** Multi-aligner strict positional consensus — all five aligners (minimap2, mapPacBio,
deSALT, uLTRA, gapmm2) must agree on the exonic region AND the CIGAR operation type at each
reference position. Positions flanked by fewer than 10 exact-match bases on each side
(`--isolation-flank 10`) are excluded to avoid alignment-artifact inflation near poorly-mapped
regions.

---

## Output Files

| File | Description |
|------|-------------|
| `error_profile_20260422/penalty_scores.tsv` | **Production table** — HP penalties (AT/CG split, isotonic-smoothed) |
| `error_profile_20260422/str_penalty_scores.tsv` | STR slippage penalties (dinucleotide/trinucleotide repeats) |
| `error_profile_20260422/error_rates_by_base.tsv` | Per-base rates (conditional D/(D+X) — M not tracked in this run) |
| `error_profile_20260422/error_counts.tsv` | Raw agreed-op counts |
| `error_profile_20260422/union_*` | "Any aligner" validation counts |
| `error_profile_20260422/diagnostics/` | 12 diagnostic plots |
| `error_profile_strand_20260501/strand_error_rates_by_base.tsv` | **Strand-split absolute rates** (D/(M+D+X), +/− strand separately) |
| `error_profile_strand_20260501/penalty_scores.tsv` | Penalty table from 3-aligner run (validation/reference) |
| `error_profile_strand_20260501/plots/` | Strand-split and absolute rate diagnostic plots |

---

## Penalty Table Format

`penalty_scores.tsv` uses a **base-class split** (AT vs CG) — not individual bases — because
AT and CG have meaningfully different substitution rates relative to their deletion rates
(see "Why AT/CG Penalties Differ" below).

```
op_type  base_class  hp_length  rate_mean  count_total  penalty_score  low_count
D        AT          1          0.0054     9101210      0.4442         False
D        CG          1          0.0052     10000942     0.8478         False
X        AT          1          0.0024     4043850      1.0000         False
X        CG          1          0.0044     7527911      1.0000         False
```

*(Values from `error_profile_strand_20260501` run — see note on rate definitions below.)*

**Columns:**
- `op_type`: D=deletion, I=insertion, X=substitution (mismatch), M=match (not emitted)
- `base_class`: AT or CG (reference base grouped by purine/pyrimidine content)
- `hp_length`: homopolymer run length (1 = isolated base)
- `rate_mean`: **absolute** error rate = op_count / (M+D+X total positions)
- `count_total`: number of agreed positions used to estimate the rate
- `penalty_score`: normalized penalty (AT sub at HP=1 → 1.0 reference per base_class)
- `low_count`: True if count_total < 50 (treat penalties as unreliable)

### ⚠️ Rate Definition Note

The `error_profile_20260422` production run did **not** accumulate M (match) events in its
counts dict, so `rate_mean` there is the **conditional** rate D/(D+X), not the absolute rate.
The `error_profile_strand_20260501` run properly tracks M, giving **absolute** rates D/(M+D+X).

| Metric | Formula | Example (AT, HP=1) | Meaning |
|--------|---------|-------------------|---------|
| Conditional (old) | D/(D+X) | 0.727 | 72.7% of errors are deletions |
| Absolute (new) | D/(M+D+X) | 0.0054 | 0.54% of positions are deleted |

> **⚠️ Correction (2026-06-01):** an earlier version of this note claimed "the normalization
> cancels the difference … both tables produce valid penalties." **That is wrong across HP
> lengths.** The penalty `sub_rate(HP=1)/rate(op,HP)` cancels the denominator only at *fixed* HP;
> across HP, `penalty_cond/penalty_abs = error_frac(HP)/error_frac(HP=1)`, which is 1.0 at HP=1 but
> grows with run length (the HP error fraction rises). Demonstrated on identical pooled counts
> (5.2 B DRS events): the D/AT penalty ratio conditional÷absolute is 1.01 at HP=1, 9.7× at HP=8,
> 16× at HP=12. So a **conditional**-derived table (e.g. the bundled DRS `penalty_scores.tsv` from
> `error_profile_20260422`) **under-weights long-homopolymer deletions by ~8–16× at HP≥8** and does
> NOT reproduce this repo's documented empirical penalties (0.44/0.17/0.033 at HP 1/4/8); an
> absolute-metric regen does (0.438/0.173/0.032). **Always derive penalties from M-tracking
> (absolute) counts.** The conditional table is mis-calibrated at long HP, not merely a different
> presentation of the same penalties.

### Penalty Formula

```
penalty(op, base_class, HP) = min(sub_rate(base_class, HP=1) / rate(op, base_class, HP), 10.0)
```

Normalized per base_class: AT substitution at HP=1 = 1.0, CG substitution at HP=1 = 1.0.
A higher error rate → lower penalty (operation is cheap/expected).
Capped at 10.0 for robustness at sparse HP lengths.

Deletion rates are further isotonic-smoothed (monotone non-decreasing with HP length)
to remove noise artifacts from small sample sizes at long HP lengths.

---

## Key Empirical Findings

### Why AT/CG Penalties Differ (2026-05-01 analysis)

The penalty difference between AT and CG is **not** primarily driven by different deletion
rates — at HP=1, absolute deletion rates are nearly identical (A=0.55%, T=0.53%, C=0.58%,
G=0.58%). The difference comes from **substitution rates**:

| Base | Abs del rate (HP=1) | Abs sub rate (HP=1) | D/X ratio | Penalty (HP=1) |
|------|---------------------|---------------------|-----------|----------------|
| A    | 0.549%              | 0.238%              | 2.3×      | —              |
| T    | 0.528%              | 0.240%              | 2.2×      | —              |
| **AT pooled** | **0.54%** | **0.24%**      | **2.3×**  | **0.44**       |
| C    | 0.584%              | 0.442%              | 1.3×      | —              |
| G    | 0.580%              | 0.434%              | 1.3×      | —              |
| **CG pooled** | **0.58%** | **0.44%**      | **1.3×**  | **0.85**       |

CG bases have ~2× higher substitution rates than AT bases. The penalty formula normalizes
deletions against the HP=1 substitution baseline, so higher substitution rate → higher
penalty for deletions (because deletions are relatively less expected compared to the
alternative error type).

Concisely: AT deletions are expected (2.3× more likely than AT substitutions); CG deletions
are less expected (only 1.3× more likely than CG substitutions). This drives the ~2× penalty
difference.

### AT/CG Grouping Validity (strand-split analysis, 2026-05-01)

To verify the AT/CG grouping, the profiler was run with `--strand-split` to separate + and −
strand reads (`error_profile_strand_20260501`). Key results:

- **A ≈ T** and **C ≈ G** on both strands for both absolute deletion and substitution rates
- No significant + vs − strand asymmetry detected — the yeast genome's roughly balanced +/−
  strand gene representation means neither strand has a systematic error bias
- The AT/CG grouping is **empirically supported** at HP=1–5 where counts are well-sampled

Strand-split outputs: `error_profile_strand_20260501/strand_error_rates_by_base.tsv`  
Strand-split plots: `error_profile_strand_20260501/plots/`

### AT vs CG Absolute Deletion Rates by HP Length

| HP | D_A (%) | D_T (%) | D_C (%) | D_G (%) | CG count note |
|----|---------|---------|---------|---------|---------------|
|  1 | 0.549   | 0.528   | 0.584   | 0.580   | well sampled  |
|  2 | 0.555   | 0.569   | 0.445   | 0.457   | well sampled  |
|  4 | 1.364   | 1.380   | 1.253   | 1.249   | well sampled  |
|  5 | 1.961   | 2.015   | 2.441   | 2.544   | adequate      |
|  6 | 2.807   | 3.035   | 5.425   | 5.341   | C:37k, G:18k  |
|  7 | 4.109   | 4.416   | 8.798   | 8.553   | C:4.7k, G:5.3k |
|  8 | 6.412   | 7.710   | 16.457  | 7.901   | C:6.8k, **G:519** |
|  9 | 8.410   | 8.103   | 21.591  | 8.411   | C:4k, **G:108** |

**C/G sparsity:** The yeast genome has very few long C or G homopolymers. Deletion event
counts collapse rapidly at HP≥7 for C/G (G has only 519 events at HP=8, 108 at HP=9).
The elevated C/G rates at HP≥6 are **noise from sparse data**, not biology. The penalty
table correctly handles this via isotonic smoothing + `low_count` flagging — CG penalties
at HP≥7 are extrapolated from the HP=5–6 trend and should be treated as approximations.

### Deletion/Substitution Rate Ratio by HP Length

The D/X ratio shows how deletion-dominated errors become at longer HP lengths:

| HP | D/X for AT | D/X for CG |
|----|-----------|-----------|
|  1 | 2.3×      | 1.3×      |
|  3 | 2.4×      | 2.9×      |
|  5 | 6.1×      | 7.8×      |
|  7 | 16.6×     | 24.4×     |

Both classes converge toward purely deletion-dominated errors at HP≥6. Substitutions at
long HP runs are strong evidence of real biological variants, not sequencing noise.

### Insertions Are Extremely Sparse

Insertion counts are below the `low_count` threshold for nearly all HP lengths.
All insertion penalties are effectively capped at 10.0, which is the correct behavior:
insertions are ~100× rarer than deletions in Nanopore reads.

### Insertions Are Extremely Sparse

Insertion counts are below the `low_count` threshold for nearly all HP lengths.
All insertion penalties are effectively capped at 10.0, which is the correct behavior:
insertions are ~100× rarer than deletions in Nanopore reads.

### Comparison to Previous Heuristic

| Operation | Context | Empirical penalty (AT) | Old heuristic | Note |
|-----------|---------|------------------------|---------------|------|
| Deletion  | HP=1    | **0.44**               | 1.00          | Heuristic ~2.3× too expensive |
| Deletion  | HP=4    | **0.17**               | 0.50          | Heuristic ~2.9× too expensive |
| Deletion  | HP=8    | **0.033**              | 0.25          | Heuristic ~7.6× too expensive |
| Substitution | HP=1 | **1.00**               | 1.00          | Identical by design |
| Insertion | any     | **10.0** (cap)         | 1.25          | Heuristic far too cheap |

**Key finding:** The old `del_normal=1.0` was substantially too expensive across all HP
lengths. The empirical table reflects that even a single isolated-base deletion is ~2×
more likely than a substitution in Nanopore reads, and deletions in long HP runs are
overwhelmingly the expected error type.

---

## How to Use with Rectify

**CLI flags:**
```bash
rectify correct \
    --junction-penalty-table common/scripts/nanopore/error_profile_20260422/penalty_scores.tsv \
    --str-penalty-table      common/scripts/nanopore/error_profile_20260422/str_penalty_scores.tsv \
    ...
```

**Batch profile YAML:**
```yaml
junction_penalty_table: /oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422/penalty_scores.tsv
str_penalty_table:      /oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422/str_penalty_scores.tsv
```

---

## How to Regenerate

Full 16-chunk run (~60 min interactive, or see `slurm/run_profiler_20260422.sh` for batch):

```bash
cd /oak/stanford/groups/larsms/Users/kevinroy

export OMP_NUM_THREADS=8; export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8; export LOKY_MAX_CPU_COUNT=8

OUT=common/scripts/nanopore/error_profile_$(date +%Y%m%d)
mkdir -p $OUT

/home/groups/larsms/users/kevinroy/anaconda3/bin/python \
    common/scripts/nanopore/empirical_cigar_error_profiler.py \
    --run-dir software/rectify/dev_runs/wt_by4742_rep1_chunked_20260412 \
    --reference software/rectify/rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa \
    --output-dir $OUT \
    --isolation-flank 10 \
    --union \
    --str-repeat \
    --workers 8 \
    2>&1 | tee $OUT/profiler.log
```

Key flags:
- `--isolation-flank 10`: only count errors flanked by ≥10 exact matches on both sides
- `--union`: also compute "any aligner" counts for validation
- `--str-repeat`: separate dinucleotide/trinucleotide repeat errors from the HP=1 baseline
- `--workers 8`: parallel chunk processing (16 chunks processed in 2 rounds of 8)

Then generate diagnostic plots:
```bash
/home/groups/larsms/users/kevinroy/anaconda3/bin/python \
    common/scripts/nanopore/plot_error_profiles.py \
    --profile-dir $OUT
```

### To Add Strand-Split Output to a New Run

Add `--strand-split` to the profiler command. Also note that `deSALT` BAMs are empty in
all but one chunk in `wt_by4742_rep1_drs_trim_20260417`, so use `--aligners minimap2 mapPacBio uLTRA`:

```bash
/home/groups/larsms/users/kevinroy/anaconda3/bin/python \
    common/scripts/nanopore/empirical_cigar_error_profiler.py \
    --run-dir software/rectify/dev_runs/wt_by4742_rep1_drs_trim_20260417 \
    --reference software/rectify/rectify/data/S288C_reference_sequence_R64-5-1_20240529.fsa \
    --output-dir $OUT \
    --isolation-flank 10 \
    --union \
    --str-repeat \
    --strand-split \
    --aligners minimap2 mapPacBio uLTRA \
    --workers 8 \
    2>&1 | tee $OUT/profiler.log
```

Output adds `strand_error_rates_by_base.tsv` with columns:
`op_type  ref_base  hp_length  strand(+/-)  count  total_positions  rate`

### Known Issues / Caveats

**deSALT failures:** Chunks 002, 011, and 015 had deSALT BAM failures (0-byte output or
SAM corruption). These chunks were skipped (min_aligners=5 not met). 13/16 chunks
contributed data. The deSALT failures are an interaction effect in its De Bruijn graph
construction — not fixable by read filtering.

**Insertion scarcity:** Insertion penalties (`I` rows) all have `low_count=True`.
Treat as 10.0 (cap) in all cases.

**Chemistry/species specificity:** These tables were derived from S. cerevisiae R10.4.1
data. HP error rates are chemistry- and species-dependent. Do not use for non-yeast or
non-R10.4.1 data without regenerating from your own dataset.

**plot_error_profiles.py column mapping:** The plot script expects `base_class` and
`rate_mean` columns but `error_rates.tsv` uses `ref_base` (A/C/G/T) and `rate`. The
script now handles this mapping automatically (added 2026-04-22).

---

*Author: Kevin R. Roy — kevinrjroy@gmail.com*  
*Last updated: 2026-05-01 (strand-split analysis; AT/CG grouping validation; absolute vs conditional rate clarification)*
