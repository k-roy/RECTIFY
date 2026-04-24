# Junction-Shift FDR: Future Implementation Proposal

**Status:** Design proposal — not yet implemented  
**Author:** Kevin R. Roy  
**Date:** 2026-04-22  
**Depends on:** `empirical_cigar_error_profiler.py`, base-specific penalty tables (`penalty_scores.tsv`)

---

## The Problem

A nanopore HP-context deletion near a donor (GT) or acceptor (AG) dinucleotide can shift
the apparent junction boundary by δ bases. If an alternative canonical motif exists at
position d+δ or a−δ in the reference, rectify may select the wrong site — producing a
junction that passes all canonical-motif and annotation filters but is offset from the
true splice site.

This failure mode is most likely when:
- The exonic flank immediately adjacent to the splice boundary contains a long A/T
  homopolymer (HP≥5), where empirical deletion rates reach 2–8%
- An alternative GT or AG dinucleotide exists within the shift distance

---

## Analytical Risk Model (First-Pass Approximation)

For junction J with donor at d and acceptor at a:

```
P_shift(b)    = del_rate(base(b), h(b)) + ins_rate(base(b), h(b))
              [from empirical penalty table, base-class-specific]

C(b, δ)       = 1 if canonical motif exists at b ± δ in reference R, else 0

P_false(J, b) = P_shift(b) × (1/h(b)) × Σ_{δ=1}^{h(b)} C(b, δ)

P_false(J)    = 1 - (1 - P_false(J, donor)) × (1 - P_false(J, acceptor))

FDR_est       = (1/N) × Σ_j P_false(J_j)      [across all N reported junctions]
```

With n independent reads supporting the same junction:
```
P_false_adjusted(J) ≈ P_false(J)^n
```
(approximation — see Caveats below)

### Rough Envelope Values

| Exonic HP context | del_rate | Max shift | P_false per boundary |
|-------------------|----------|-----------|----------------------|
| HP=1 (any base)   | 0.005    | 1 bp      | ~0.0006              |
| HP=3 (A/T)        | 0.007    | 3 bp      | ~0.002               |
| HP=5 (A/T)        | 0.018    | 5 bp      | ~0.009               |
| HP=8 (A/T)        | 0.081    | 8 bp      | ~0.05                |
| HP=10 (A/T)       | 0.13     | 10 bp     | ~0.10                |
| HP=5 (C/G)        | 0.003    | 5 bp      | ~0.001               |

These use a genome-wide GT/AG background frequency (~1/8) as a placeholder for C(b,δ).
**These numbers are upper-bound estimates** — see Caveats.

---

## Critical Caveats Before Implementation

### 1. The model ignores rectify's existing defenses

Rectify's HP-anchored DP scoring with empirical penalties is specifically designed to
prefer the site with the best sequence match. The canonical HP prior (`_CANONICAL_HP_PRIOR
= 0.5` in v3.1.7) requires alternative sites to score >0.5 better than canonical
junctions to displace them. The analytical model above assumes P_rectify_selects ≈ 1
(i.e., rectify selects any nearby canonical site blindly), which is a gross overestimate.

**The table values above may be overestimates by 3–10× for this reason.**

### 2. Wrong input rates for junction-shift

The CIGAR profiler rates measure per-position errors where all 5 aligners agree. Junction
shift is an alignment ambiguity problem: a read with an HP deletion may produce a sequence
that aligners consistently map to the shifted site without any disagreement. P_shift from
the profiler is not directly the per-read junction-shift probability.

### 3. GT/AG background frequency should be local, not genome-wide

The ~1/8 background rate for GT or AG at a random position ignores:
- Evolutionary depletion of spurious splice sites in exonic sequence near real junctions
- Yeast's actual local sequence composition (~62% AT intronic content vs ~38% exonic)

Compute C(b, δ) empirically from the reference at annotated junction positions, not from
a genome-wide average.

### 4. P_false^n independence assumption is violated for alignment bias

If a long HP run systematically biases the aligner toward a shifted junction for every
read at that locus, then P_false^n is extremely optimistic — the error is correlated
across reads, not independent. P_false^n captures independent sequencing errors well but
not systematic aligner preferences.

### 5. STR context uses str_del_cost, not HP=1 rate

Dinucleotide/trinucleotide repeats near splice boundaries (e.g., TATATA before a GT)
should use `HpPenaltyTable.str_del_cost(unit, n_copies)` from `str_penalty_scores.tsv`,
not the HP=1 deletion rate. The STR profiler (`--str-repeat` flag in
`empirical_cigar_error_profiler.py`) generates the required table.

### 6. Vulnerability radius and HP context are orthogonal

The maximum HP run length that matters is determined by the error rate (longer runs →
higher error rate), not by rectify's `max_boundary_shift`. A 10-base HP run near a donor
is relevant even if the search radius is only ±3 bp, because a 1-base deletion (the most
common case) still shifts the boundary by 1 bp.

---

## Recommended Implementation Roadmap

### Step 1 (FIRST): Simulation calibration on truth set (Option C)

Before implementing any annotation or correction logic, validate the analytical model:

1. Select 200–500 high-confidence junctions with ≥10 independent reads
2. Introduce synthetic HP-context deletions at empirical rates from `penalty_scores.tsv`
3. Run `rectify correct` on the perturbed reads
4. Measure actual junction misassignment rate by HP context
5. Compare observed rates to analytical model predictions

If observed rates match the analytical model within 2×, the model is useful.
If the model overestimates by >5×, recalibrate before building Option A.

**If simulation shows <1% actual misassignment even at HP=8 A/T contexts (plausible given
rectify's scoring), revise the flagging thresholds accordingly before any output annotation.**

### Step 2: Per-junction annotation (Option A)

Add `junction_shift_risk` column to `corrected_3ends.tsv`:

```python
def compute_junction_shift_risk(
    donor_pos: int,
    acceptor_pos: int,
    ref_seq: str,
    penalty_table: HpPenaltyTable,
    str_penalty_table: Optional[dict] = None,
    n_reads: int = 1,
    flank_window: int = 20,
) -> float:
    """Return P_false(J) adjusted for read depth."""
    ...
```

Threshold tiers (to be calibrated against Step 1 simulation):

| Risk tier   | P_false(J) | Recommended action                                       |
|-------------|------------|----------------------------------------------------------|
| Low         | < 0.01     | No concern                                               |
| Moderate    | 0.01–0.05  | Flag; consider multi-read validation                     |
| High        | 0.05–0.15  | Require ≥2 independent reads supporting same junction    |
| Very high   | > 0.15     | Do not report as single-read junction                    |

### Step 3: Scoring prior integration (Option B, if Step 1 warrants it)

Only implement if simulation (Step 1) shows meaningful misassignment rates (>2% at HP≥5).
Modify `refine_bam_junctions` to penalize candidate junctions reachable by a single
HP-deletion from a known long-HP exonic flank. This is a Bayesian prior on the candidate
selection, not a hard filter.

---

## Files to Modify

| File | Change |
|------|--------|
| `core/junction_refiner.py` | Add `_junction_shift_risk()` helper; call from `refine_bam_junctions` |
| `core/bam_processor.py` | Pass `junction_shift_risk` into unified record |
| `core/unified_record.py` | Add `junction_shift_risk` field |
| New: `tests/test_junction_shift_fdr.py` | Unit tests; simulation fixture |

---

## Related Code

- `common/scripts/nanopore/empirical_cigar_error_profiler.py` — generates `penalty_scores.tsv` and `str_penalty_scores.tsv`
- `rectify/core/junction_refiner.py` — `HpPenaltyTable`, `_hp_run_length`, `_str_repeat_info`, `_score_hp_anchored`
- `rectify/core/junction_refiner.py:_CANONICAL_HP_PRIOR` — existing defense mechanism that this model must account for
