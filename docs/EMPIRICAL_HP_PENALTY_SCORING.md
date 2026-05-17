# Empirical Homopolymer-Context Penalty Scoring

**Status:** Production — implemented in `junction_refiner.py` (Module 2H) and `splice_aware_5prime.py`  
**Author:** Kevin R. Roy  
**Date:** 2026-05-01  
**Dataset:** `wt_by4742_rep1` — *S. cerevisiae* wild-type, R10.4.1 chemistry, Dorado basecaller (9.7M reads)

---

## What It Is

Rectify's junction refinement (Module 2H) and 3'SS rescue scoring use a
**homopolymer-context-aware edit distance** when comparing candidate splice junctions
against read sequence. By default these modules use fixed heuristics (`del=1.0` for
all contexts, `del=0.5` for HP≥4, `ins=1.25`). Empirical penalty tables derived from
multi-aligner consensus on real Nanopore data replace these heuristics with
data-driven values.

### What the Penalty Score Means

A penalty score quantifies how *surprising* a given edit operation is in a given
sequence context. The normalization is:

```
penalty(op, base_class, HP) = min( sub_rate(base_class, HP=1) / rate(op, base_class, HP), 10.0 )
```

Where:
- `rate(op, base_class, HP)` is the absolute error rate D/(M+D+X) for operation `op` at
  homopolymer length `HP` for bases in `base_class` (AT or CG)
- `sub_rate(base_class, HP=1)` is the substitution rate at HP=1 for the same base class
  (the reference baseline — substitutions at isolated bases are the most "generic" error)
- Higher deletion rate → lower penalty (deletions are cheap/expected in that context)
- Capped at 10.0 for sparse HP lengths where rates are unreliable

**Reference point:** `penalty(X, base_class, HP=1) = 1.0` for both AT and CG by definition.
All other operations are expressed relative to this substitution baseline.

---

## Why It Matters: Old Heuristics vs Empirical Values

| Operation | Context | Empirical penalty (AT) | Old heuristic | Factor off |
|-----------|---------|:----------------------:|:-------------:|:----------:|
| Deletion  | HP=1    | **0.44**               | 1.00          | 2.3× too expensive |
| Deletion  | HP=4    | **0.17**               | 0.50          | 2.9× too expensive |
| Deletion  | HP=8    | **0.033**              | 0.25          | 7.6× too expensive |
| Substitution | HP=1 | **1.00** (ref)         | 1.00          | Identical by design |
| Insertion | any     | **10.0** (cap)         | 1.25          | Heuristic far too cheap |

The old `del_normal=1.0` treated a deletion at HP=1 as equally costly as a substitution.
Empirically, even a single isolated-base deletion is **2.3× more likely** than a substitution
in Nanopore R10.4.1 reads. At HP=8 the disparity reaches 7.6×. Using over-priced deletion
costs causes the scorer to incorrectly prefer junction candidates that avoid deletions even
when those deletions are the expected Nanopore error mode.

Insertions are ~100× rarer than deletions in Nanopore reads. All insertion penalties are
effectively at the 10.0 cap — the old heuristic of 1.25 was dramatically too cheap.

---

## Base Class Split: AT vs CG

The penalty table splits reference bases into two classes — **AT** and **CG** — rather than
treating all four bases identically. This split is empirically justified and driven by
**substitution rates**, not deletion rates.

### The Numbers (HP=1, absolute rates)

| Base class | Abs del rate (HP=1) | Abs sub rate (HP=1) | D/X ratio | Penalty (HP=1) |
|------------|:-------------------:|:-------------------:|:---------:|:--------------:|
| AT pooled  | 0.54%               | 0.24%               | 2.3×      | **0.44**       |
| CG pooled  | 0.58%               | 0.44%               | 1.3×      | **0.85**       |

Deletion rates at HP=1 are nearly identical across all four bases (A=0.55%, T=0.53%,
C=0.58%, G=0.58%). The penalty difference arises because **CG bases have ~2× higher
substitution rates** than AT bases. Since the penalty formula normalizes deletion cost
against the substitution baseline, a higher substitution rate means deletions are
relatively less expected compared to the alternative error type — hence a higher penalty.

Concisely: AT deletions are the dominant error mode (2.3× more likely than AT subs);
CG deletions are less dominant (only 1.3× more likely than CG subs). This drives the
~2× penalty difference between AT and CG at HP=1.

**Figure — Absolute deletion rate vs HP length, pooled strands:**  
`common/scripts/nanopore/error_profile_strand_20260501/plots/abs_del_rate_by_base.png`

---

## Absolute Deletion Rates by Base and HP Length

| HP | D_A (%) | D_T (%) | D_C (%) | D_G (%) | CG count note |
|----|:-------:|:-------:|:-------:|:-------:|:-------------:|
|  1 | 0.549   | 0.528   | 0.584   | 0.580   | well sampled  |
|  2 | 0.555   | 0.569   | 0.445   | 0.457   | well sampled  |
|  4 | 1.364   | 1.380   | 1.253   | 1.249   | well sampled  |
|  5 | 1.961   | 2.015   | 2.441   | 2.544   | adequate      |
|  6 | 2.807   | 3.035   | 5.425   | 5.341   | C:37k, G:18k  |
|  7 | 4.109   | 4.416   | 8.798   | 8.553   | C:4.7k, G:5.3k |
|  8 | 6.412   | 7.710   | 16.457  | 7.901   | C:6.8k, **G:519** |
|  9 | 8.410   | 8.103   | 21.591  | 8.411   | C:4k, **G:108** |

Rate = D/(M+D+X). Absolute rates: fraction of all covered reference positions deleted.

Note: A≈T and C≈G across all HP lengths in the well-sampled range (HP=1–6). The
AT/CG grouping is empirically valid — there is no justification for a 4-way A/T/C/G
split based on these data.

---

## D/X Ratio by HP Length

Both base classes converge toward purely deletion-dominated errors at longer HP runs:

| HP | D/X for AT | D/X for CG |
|:--:|:----------:|:----------:|
|  1 |    2.3×    |    1.3×    |
|  3 |    2.4×    |    2.9×    |
|  5 |    6.1×    |    7.8×    |
|  7 |   16.6×    |   24.4×    |

At HP≥6, substitutions at homopolymers are extremely rare (>90% of errors are
deletions). A substitution at a long HP run is strong evidence of a real biological
variant, not sequencing noise.

---

## C/G Sparsity at High HP Lengths

The *S. cerevisiae* genome has very few long C or G homopolymers. Deletion event
counts collapse rapidly:

| HP | AT events | C events | G events |
|:--:|:---------:|:--------:|:--------:|
|  6 | >2M       | 37k      | 18k      |
|  7 | >1M       | 4.7k     | 5.3k     |
|  8 | >500k     | 6.8k     | **519**  |
|  9 | >200k     | 4k       | **108**  |

The elevated C/G deletion rates at HP≥7 (e.g., C at HP=8: 16.5%) are a **data sparsity
artifact**, not biology. With only 519 G deletion events at HP=8, the rate estimate is
dominated by sampling noise.

The penalty table handles this correctly:
1. **Isotonic smoothing**: deletion rates are forced monotone non-decreasing with HP length.
   The CG rate at HP=8 (noisy) is smoothed down to be consistent with HP=5–6.
2. **`low_count` flagging**: rows with fewer than 50 events are flagged as unreliable.
3. **Cap at 10.0**: sparse, unreliable high-HP penalties are capped at 10.0.

CG penalties at HP≥7 should be treated as extrapolations from the HP=5–6 empirical trend.

**Figure — Deletion event counts by HP length and base (log scale):**  
`common/scripts/nanopore/error_profile_strand_20260501/plots/hp_deletion_counts_sparsity.png`

---

## Strand-Split Validation

To verify that the AT/CG grouping is not confounded by strand biases, the profiler was
re-run with `--strand-split` on `wt_by4742_rep1_drs_trim_20260417` using three aligners
(minimap2, mapPacBio, uLTRA; deSALT excluded due to near-zero coverage outside chunk_003).

Results (`error_profile_strand_20260501/strand_error_rates_by_base.tsv`):

- **A ≈ T** and **C ≈ G** on both + and − strands for both absolute deletion and
  substitution rates at HP=1–5 (well-sampled range)
- No significant + vs − strand asymmetry in deletion or substitution rates
- The yeast genome's roughly balanced +/− strand gene representation means neither
  strand has a systematic error bias at the read level

These findings empirically support the AT/CG grouping and confirm that the penalty
table is not strand-biased.

**Figure — Absolute deletion rate by base and strand (8 lines):**  
`common/scripts/nanopore/error_profile_strand_20260501/plots/abs_del_rate_by_base_strand.png`

**Figure — Conditional deletion rate by base and strand (fraction D/(D+X)):**  
`common/scripts/nanopore/error_profile_strand_20260501/plots/del_rate_by_base_strand.png`

---

## Penalty Curves

**Figure — Empirical HP deletion penalty scores (AT and CG):**  
`common/scripts/nanopore/error_profile_strand_20260501/plots/hp_penalty_summary.png`

The curves show:
- AT penalties fall from 0.44 at HP=1 to ~0.03 at HP=8 (deletions become nearly free)
- CG penalties fall from 0.85 at HP=1 to ~0.18 at HP=8 (slower decay due to higher sub rate)
- X marks indicate extrapolated (low-count) points for CG at HP≥7
- Both curves are isotonic-smoothed (monotone non-decreasing HP→penalty-score is
  inverted: monotone non-increasing from HP=1 to longer runs)

---

## How the Profiler Works

Error rates are derived from multi-aligner strict positional consensus using
`empirical_cigar_error_profiler.py`. The method:

1. Load BAMs from 3–5 independent aligners (minimap2, mapPacBio, uLTRA, deSALT, gapmm2)
2. For each read present in all aligners, find positions where ALL aligners agree on:
   - The exonic region assignment (read is in a transcribed region)
   - The CIGAR operation type (M/D/X at that reference position)
3. Require ≥10 exact-match bases flanking each position (`--isolation-flank 10`) to
   exclude alignment artifacts near poorly-mapped regions
4. Accumulate per-position (op_type, ref_base, hp_length) counts
5. Compute absolute rates D/(M+D+X) per (op, base, HP) cell
6. Normalize to penalty scores; apply isotonic smoothing to deletion rates

The multi-aligner consensus approach filters out aligner-specific artifacts: a deletion
that only one aligner calls is not counted. This makes the rates conservative (true
deletion rate is slightly higher) but artifact-free.

**Production run:** `error_profile_20260422/` — 5 aligners, 16 chunks, ~13/16 chunks
contributed (deSALT failures on chunks 002, 011, 015).  
**Strand-split validation run:** `error_profile_strand_20260501/` — 3 aligners, 8 chunks,
absolute rate definition (M properly tracked).

---

## Penalty Table Format

`penalty_scores.tsv` columns:

| Column | Description |
|--------|-------------|
| `op_type` | D=deletion, I=insertion, X=substitution |
| `base_class` | AT or CG (reference base grouped) |
| `hp_length` | Homopolymer run length (1 = isolated base) |
| `rate_mean` | Absolute error rate = op_count / (M+D+X) |
| `count_total` | Number of agreed positions used |
| `penalty_score` | Normalized penalty (sub HP=1 = 1.0 per base_class) |
| `low_count` | True if count_total < 50 (treat as unreliable) |

Example rows:

```
op_type  base_class  hp_length  rate_mean  count_total  penalty_score  low_count
D        AT          1          0.0054     9101210      0.4442         False
D        CG          1          0.0052     10000942     0.8478         False
X        AT          1          0.0024     4043850      1.0000         False
X        CG          1          0.0044     7527911      1.0000         False
D        AT          8          0.0641     ...          0.0374         False
D        CG          8          ~est       ...          ~0.18          True
I        AT          1          ...        ...          10.0           True
```

---

## How to Use with Rectify

### Command line

```bash
rectify correct \
    --junction-penalty-table /oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422/penalty_scores.tsv \
    --str-penalty-table      /oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422/str_penalty_scores.tsv \
    [other flags ...]
```

### Batch profile YAML

```yaml
junction_penalty_table: /oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422/penalty_scores.tsv
str_penalty_table:      /oak/stanford/groups/larsms/Users/kevinroy/common/scripts/nanopore/error_profile_20260422/str_penalty_scores.tsv
```

Both flags affect Module 2H (`junction_refiner.py`) and the 3'SS rescue scorer
(`splice_aware_5prime.py`). If neither flag is set, the built-in heuristics are used.

### The canonical HP prior

When the empirical table is active, a non-canonical junction must beat the best canonical
alternative by more than `_CANONICAL_HP_PRIOR = 0.5` edit-distance units to win. This
equals approximately one HP deletion penalty equivalent (del_cost ≈ 0.44 at HP=1),
ensuring canonical junctions win within the noise floor. This prior was calibrated for
the empirical table and works correctly with the default heuristics as well.

---

## Caveats

- **Chemistry and species specificity**: These tables are derived from *S. cerevisiae*
  R10.4.1 data. HP error rates are chemistry- and species-dependent. Do not use for
  non-yeast or non-R10.4.1 data without regenerating from your own dataset.
- **Insertion scarcity**: All insertion penalties (`I` rows) have `low_count=True`.
  Treat as 10.0 in all cases.
- **CG at HP≥7**: Extrapolated via isotonic smoothing from HP=5–6 empirical trend.
  The actual penalty for CG deletions at HP≥7 is uncertain due to data sparsity in
  the *S. cerevisiae* genome.

---

## See Also

- [Junction-Shift FDR](JUNCTION_SHIFT_FDR.md) — theoretical risk model for
  HP-context junction-shift false positives; uses these rates as inputs.
- The penalty profiler and the production tables are maintained outside this
  repository. To regenerate the tables, run
  `rectify/calibration/calibrate_shift_corrections.py` on a CIGAR-error
  profile produced from a high-coverage WT BAM and a reference of the same
  organism; the resulting JSON is consumed by
  `--junction-penalty-table` / `--str-penalty-table` on `rectify correct`
  and `rectify run-all`.

---

*Author: Kevin R. Roy — kevinrjroy@gmail.com*  
*Last updated: 2026-05-01*
