# DRS Poly(A) and Adapter Analysis — Findings & Algorithm

**Author:** Kevin R. Roy  
**Date:** 2026-04-17  
**Status:** Production implementation complete  
**Implementation:** `projects/roadblocks/scripts/nanopore_analysis/analyze_drs_adapter.py`

---

## 1. DRS Orientation and Dorado Output Convention

Oxford Nanopore DRS (SQK-RNA004) sequences RNA **3′→5′** through the nanopore. The physical molecule structure during sequencing is:

```
[DNA adapter: GGCTTCTTCTT] — [poly(A)] — [gene body] — 5' cap
                                          → sequencing direction
```

Dorado basecalls in the pore direction and then **reverse-complements** (for RNA) to output reads in the conventional **5′→3′ RNA orientation**. In BAM files:

| Strand (BAM) | Query sequence orientation | Adapter location |
|---|---|---|
| `is_reverse=False` (plus) | 5′→3′ RNA | 3′ end of query_sequence |
| `is_reverse=True` (minus) | stored rev-comp | After `reverse_complement(query)`, adapter is at 3′ end |

**Rule for analysis:** Always orient the read to RNA 5′→3′ (rev-comp minus-strand reads), then the **poly(A) tail + adapter stub** are always at the **3′ end** of the resulting sequence.

---

## 2. The GGCTTCTTCTT Adapter and Dorado Pre-Trimming

The ONT SQK-RNA004 **DNA adapter** sequence is:

```
5'-GGCTTCTTCTT-3'
```

This hybridizes to the poly(A) tail to initiate capture. After sequencing, this adapter signal immediately follows the poly(A) current signature in the raw electrical trace.

**Critical behavior: Dorado pre-trims high-confidence adapter matches.** Dorado's internal adapter detection removes the adapter from the basecalled read when the signal match is unambiguous. What appears in BAM files is therefore a *residual stub* — the portion of the adapter that Dorado was not confident enough to trim.

### The Residual Stub Pattern

Empirically, what is observed after the poly(A) tail in basecalled reads is:

```
...AAAAAAA T [CT]{0,10}
```

- **Position 1 always T:** The first decoded base of the residual stub is consistently T. The adapter begins with `G`; how the basecaller renders that G at the poly(A)/adapter signal transition is not fully understood — it may reflect a boundary-current artifact or a systematic model behavior. The mechanism has not been empirically characterized in the literature (as of 2026).
- **Positions 2+: C or T only** (no G, no A). Importantly, the adapter sequence after the first two G's is `CTTCTTCTT` — which is entirely C and T. So the C/T observed at stub positions 2+ likely reflects **correctly-called adapter bases**, not model confusion. Reads that retain more than 1–2 bases of adapter beyond Dorado's trim point simply show the natural C/T content of the adapter.
- **Absence of G at positions 2+:** The second adapter base is G, which would appear at stub position 2 only in reads where Dorado left ≥2 bases untrimmed. These appear to be a small minority, and when G is emitted at that position it may be miscalled or absent below detection threshold. The data do not resolve whether G→T/C substitution occurs or whether Dorado simply trims more aggressively when G is confidently detected.

---

## 3. Why We See T (Not G) at Position 1

The adapter begins with `GG`, but the first decoded stub base is consistently `T` across all datasets. Two plausible explanations exist, and they are not mutually exclusive:

1. **Boundary-current artifact:** The electrical transition from the poly(A) homopolymer current into the adapter current occurs over 1–3 bases. The model, trained on RNA, maps this transitional signal to the nearest learned state, which appears to be T (a common RNA base with a distinctive current). This would be a transient decoding artifact at the poly(A)/adapter junction.

2. **Structural feature of the adapter ligation:** The SQK-RNA004 adapter ligation involves a short RNA/DNA junction. The precise nucleotide at position 1 of the residual stub may reflect the ligation chemistry rather than the adapter G itself.

**The mechanism is unverified** — this T is an empirical observation, not a published finding. What is established: (a) it is extremely consistent across 400k reads and 4 independent datasets, and (b) it is not part of the poly(A) tail. Its presence as a near-universal anchor makes `T[CT]{0,10}$` a practical and specific detection regex regardless of mechanistic understanding.

Also note: the C/T pattern at positions 2+ does **not** require invoking model confusion — the adapter sequence `GGCTTCTTCTT` is itself predominantly C and T after the first two G's. See Section 2 for details.

---

## 4. Adapter Detection Algorithm

### 4.1 Overview: Two-Pass Approach

```
Input: 3′-oriented read sequence (last 150 bp window)

Pass 0 (no regex needed): Dorado already removed adapter → poly(A) runs to end of read
Pass 1: regex T[CT]{0,10}$ detects residual stub
Pass 2: iterative peel rescues stubs with A-basecalling errors
```

### 4.2 Constants

```python
_ADAPTER_RE   = re.compile(r'T[CT]{0,10}$')
_PASS2_MAX_STUB = 15    # maximum stub length to try peeling
MIN_POLYA_PASS2 = 5     # minimum poly(A) required to accept a pass-2 call
```

### 4.3 Strict 3′-Anchored Poly(A) Scan

The poly(A) scan **must** be anchored at the actual 3′ end of the (adapter-trimmed) sequence. An earlier implementation found the "rightmost A anywhere in the sequence" and extended left — this produced false positives when isolated A's appeared in the gene body (e.g., a read ending `...GCAACGT` would report poly(A)=2 from the AA upstream of CGT).

**Correct implementation:**

```python
def _scan_polya(seq: str, max_error_rate: float = 0.15) -> int:
    """Strict right-to-left poly(A) scan anchored at the actual 3′ end.
    
    Walks from seq[-1] leftward; stops as soon as the cumulative non-A
    rate exceeds max_error_rate. A sequence ending in C/G/T immediately
    returns 0 (no tolerable extension).
    """
    n = len(seq)
    errors = 0
    total = 0
    polya_start = n
    for i in range(n - 1, -1, -1):
        total += 1
        if seq[i] != 'A':
            errors += 1
        if errors / total > max_error_rate:
            break
        polya_start = i
    return n - polya_start
```

**Key invariant:** `seq[-1]` must be `A` for `_scan_polya` to return > 0. Any read whose last base is non-A gets poly(A)=0 unless pass-2 rescues it.

### 4.4 Full `find_polya_and_adapter` Function

```python
def find_polya_and_adapter(seq, max_error_rate=0.15, min_polya=1, adapter_window=50):
    """Detect poly(A) tail and adapter stub from a 3′-oriented read window.
    
    Returns
    -------
    (polya_len, adapter_seq, last_base, adapter_pass)
    
    adapter_pass values:
        0 = no adapter stub found (Dorado pre-trimmed, or genuine 3′ end)
        1 = Pass 1 regex T[CT]{0,10}$ matched
        2 = Pass 2 iterative peel recovered hidden poly(A)
    """
    n = len(seq)
    if n == 0:
        return 0, "", "", 0

    # --- Pass 1: regex adapter detection ---
    m = _ADAPTER_RE.search(seq)
    if m:
        adapter_seq = seq[m.start():m.start() + adapter_window]
        seq = seq[:m.start()]
        adapter_pass = 1
    else:
        adapter_seq = ""
        adapter_pass = 0

    n = len(seq)
    if n == 0:
        return 0, adapter_seq, "", adapter_pass

    last_base = seq[-1]
    polya_len = _scan_polya(seq, max_error_rate)

    if polya_len >= min_polya:
        return polya_len, adapter_seq, last_base, adapter_pass

    # --- Pass 2: iterative peel for stubs with A-basecalling errors ---
    # Example: seq = "AAAAAATACC"
    #   k=1: candidate="AAAAAAATAC", end='C' → not A, skip
    #   k=2: candidate="AAAAAAATA",  end='A' → check; stub="CC" has no T → skip
    #   k=3: candidate="AAAAAAATA",  stub="ACC" has no T → skip
    #   k=4: candidate="AAAAAA",     end='A' → stub="TACC" has T → polya=6 ✓
    if last_base != 'A':
        for k in range(1, _PASS2_MAX_STUB + 1):
            if len(seq) <= k:
                break
            candidate = seq[:-k]
            if not candidate or candidate[-1] != 'A':
                continue
            stub = seq[-k:]
            if 'T' not in stub:
                # No adapter boundary T → likely genuine gene body termination
                continue
            polya_len2 = _scan_polya(candidate, max_error_rate)
            if polya_len2 >= MIN_POLYA_PASS2:
                combined = (stub + adapter_seq)[:adapter_window]
                return polya_len2, combined, candidate[-1], 2

    return 0, adapter_seq, last_base, adapter_pass
```

### 4.5 Why the T-Requirement in Pass 2

Pass-2 requires that the peeled stub contains `T`. This prevents false positives:

- **Genuine gene body terminations** end in non-A bases (e.g., `AAAAACGT`). Peeling `CGT` leaves `AAAAA`, which has no T in the stub → correctly rejected.
- **Real adapter stubs** always contain T (the adapter boundary base), even if surrounded by other errors (e.g., `TACC`, `TCC`, `TTCT`).

Without the T-requirement, pass-2 would accept any sequence ending with a short non-A run after sufficient poly(A), producing high false-positive rates on gene body reads.

---

## 5. Empirical Findings

### 5.1 Pass Breakdown (400k reads, 4 datasets)

| Pass | Description | % of reads |
|---|---|---|
| 0 | Dorado pre-trimmed adapter; poly(A) runs to end | 53.2% |
| 1 | `T[CT]{0,10}$` regex matched residual stub | 45.7% |
| 2 | Iterative peel recovered hidden poly(A) | 1.1% |
| — | Genuine poly(A)=0 (gene body 3′ ends) | 2.14% |

### 5.2 The 2.14% Poly(A)=0 Population

Reads with poly(A)=0 are **not** an artifact — they represent transcripts that lack a poly(A) tail (e.g., stable 3′-processed RNA intermediates, degradation products, or genuine non-polyadenylated isoforms). These reads:
- Terminate in non-A bases (C, G, T) — the sequence ends before any A-tract
- Have low or zero `pt:i:` Dorado poly(A) estimates
- Are distributed across all samples, not concentrated in any single experimental condition

### 5.3 Adapter Position Composition

Stacked base-frequency analysis across positions 1–50 after the poly(A) tail:

| Position | Predominant base | Notes |
|---|---|---|
| +1 | T (~95%) | Adapter boundary base; always T |
| +2–+10 | C or T (~90%) | RNA-model decoding of GCTTCTTCTT |
| +11+ | Sparse / mixed | Fewer reads extend this far (deeper Dorado trimming) |

No G at position 2+. This falsified the rationale for the broader `T[^A]{0,10}$` regex.

### 5.4 Regex FDR Analysis

Evaluated on 1M random 150-mers (uniform ACGT):

| Regex | FDR (% random sequences matched) |
|---|---|
| `T[CT]{0,10}$` | 33.3% |
| `T[^A]{0,10}$` | 50.0% |

The `T[^A]` expansion was rejected because:
1. Empirical data shows **no G** at positions 2+ of real adapter stubs
2. FDR increases by ~17 percentage points with no sensitivity gain
3. G residuals in real data are captured by pass-2 (they arise from A-error stubs)

### 5.5 Sequence-Based vs. Signal-Based Poly(A) Length

Two independent measurements:

| Metric | Source | Behavior |
|---|---|---|
| `polya_len` | Sequence-based (this algorithm) | Sensitive to ≥1 bp; can have false positives from isolated A's if scan is not strictly 3′-anchored |
| `pt:i:` tag | Dorado signal-based | Effective minimum threshold ~5–6 bp; systematically ~60% longer than sequence-based (median ~35 vs ~20 bp) |

**Validation of the strict-anchor fix:** Under the old "rightmost A" approach, reads with seq poly(A) 1–5 bp showed 31.6% pt:i:=0. Under strict 3′-anchored scan, this drops dramatically, as short poly(A) calls from isolated upstream A's no longer occur.

**Why pt:i: is systematically longer:** The signal-level current captures the poly(A) homopolymer more fully, whereas sequence-based counting is bounded by the 3′ end of the aligned/basecalled read. Some poly(A) bases may be consumed by the adapter trimming step or soft-clipped by the aligner.

---

## 6. Implementation Notes

### 6.1 File Locations

```
scripts/nanopore_analysis/analyze_drs_adapter.py    # Main analysis script
scripts/nanopore_analysis/build_drs_html_report.py  # HTML report builder
scripts/nanopore_analysis/plot_drs_supplementary.py # Supplementary plots

processed_data/drs_adapter_analysis/
    adapter_analysis_{dataset}.tsv     # Per-read: polya_len, pt_tag, adapter_seq,
                                       #           last_base, adapter_pass, sample, ...
    position_frequency_{dataset}.tsv   # Position × base frequency matrix
    position_frequency_combined.tsv    # Combined across all datasets

plots/drs_adapter_analysis/
    executive_summary.html             # Full interactive HTML report (13 plots)
    adapter_position_composition.png   # Stacked bar: base fraction at each position
    polya_length_per_sample.png        # Two-panel violin: seq poly(A) + pt:i: per sample
```

### 6.2 Output Columns in `adapter_analysis_{dataset}.tsv`

| Column | Description |
|---|---|
| `read_name` | BAM query name |
| `sample` | Sample identifier |
| `strand` | `+` or `-` (original BAM strand) |
| `polya_len` | Empirical poly(A) length (bp), 0 if not detected |
| `pt_tag` | Dorado `pt:i:` signal-based poly(A) length (-1 if absent) |
| `adapter_seq` | First 50 bases of detected adapter stub (empty if pass-0) |
| `last_base` | Last base of poly(A)-trimmed sequence (quality check) |
| `adapter_pass` | 0/1/2 (see pass breakdown above) |
| `dataset` | Dataset identifier |

### 6.3 Running the Analysis

```bash
# Full run: 100k reads per dataset (4 datasets × 100k = 400k total)
PYTHON="python"  # use your environment's Python
$PYTHON analyze_drs_adapter.py

# Regenerate HTML report after re-running analysis
$PYTHON build_drs_html_report.py

# Regenerate supplementary plots
$PYTHON plot_drs_supplementary.py
```

---

## 7. Key Lessons and Pitfalls

1. **Always anchor poly(A) scan strictly at 3′ end.** Any implementation that finds the "rightmost A" and extends left will produce false positives on reads terminating in gene bodies. The scan must start at `seq[-1]` and stop immediately when error rate is exceeded.

2. **Dorado pre-trimming changes the problem.** The analysis is not finding the full `GGCTTCTTCTT` adapter — it is finding a short residual stub that Dorado left because it was below its confidence threshold. Longer stubs are rarer; most reads show 0–6 adapter bases.

3. **The T at position 1 is a boundary artifact, not a poly(A) base.** Do not count it as part of the poly(A) tail. The regex `T[CT]{0,10}$` correctly removes it as part of the adapter stub.

4. **Pass-2 catches ~1.1% additional reads** — mostly reads with A-basecalling errors in the adapter stub (e.g., `TACC` instead of `TCC`). These are genuine poly(A) tails; the stub is just misread.

5. **`T[^A]{0,10}$` is too permissive.** No G at adapter positions 2+ in empirical data. Expanding to `[^A]` raises FDR by 17 pp. Do not use it.

6. **MinKNOW-basecalled samples (dst1d replicates)** have lower poly(A) detection rates and different `pt:i:` behavior than Dorado-basecalled samples. These are marked in reports/plots with `*MinKNOW`.
