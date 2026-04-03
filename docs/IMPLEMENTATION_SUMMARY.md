# RECTIFY Implementation Summary

**Version:** 2.7.0b1
**Date:** 2026-03-31
**Status:** Beta Validation — 5-Aligner Pipeline, Fair Scoring, SLURM Batch Fixes

---

## Current State

RECTIFY is a fully-implemented RNA 5'/3' end correction framework for direct RNA nanopore sequencing. The core pipeline (`run-all`) runs a 5-aligner consensus alignment, corrects 3' end artifacts, filters spike-ins, and (with ≥2 samples) performs CPA clustering, DESeq2 differential expression, and APA shift analysis.

**Active development focus:** Beta validation on 12 chrI samples (4 conditions × 3 replicates) on the Sherlock HPC cluster.

---

## Architecture Overview

```
Input: FASTQ (direct RNA nanopore)
   │
   ▼
[Step 1] 5-Aligner Consensus (multi_aligner.py + consensus.py)
   ├── minimap2     (splice-aware, fast)
   ├── mapPacBio    (PacBio RNA mode — forces mismatches at splice junctions)
   ├── gapmm2       (gap-aware minimap2 variant)
   ├── uLTRA        (junction-aware long-read aligner)
   └── deSALT       (de novo splice detection, de novo mode only — annotation-
                     guided mode disabled due to SIGSEGV in yeast GTF parser)
        │
        ▼ Per-read scoring:
        ├── 5' effective clip penalty (-2/bp): max(explicit soft-clip,
        │   terminal mismatch/indel region for mapPacBio fairness)
        ├── 5' junction rescue: sequence-match to upstream exon → no penalty
        ├── 3' A-tract depth penalty (-1/downstream A, capped at 10)
        ├── 3' non-polyA terminal error penalty (-2/bp, capped at 10):
        │   aligners that stop before true 3' end with non-A content penalized
        └── Tiebreakers: canonical GT-AG motifs, annotated junction support,
            majority 3' position vote across aligners
   │
   ▼
[Step 2] 3' End Correction (bam_processor.py)
   ├── A-tract ambiguity walk-back (universal)
   ├── Indel artifact correction (deletion artifacts in poly(A) alignment)
   ├── False junction removal (N operations in poly(A) tail → discarded)
   ├── Poly(A) tail length measurement
   └── Spike-in filtering (ENO2 k-mer matching)
   │
   ▼
[Step 3] Downstream Analysis (analyze_command.py + analyze/)
   ├── CPA cluster detection
   ├── DESeq2 gene-level and cluster-level differential expression
   ├── APA shift analysis
   ├── Genomic distribution, motif discovery, GO enrichment
   └── HTML report generation
```

---

## Critical: Two `run-all` Implementations

**Do not confuse these two files.**

| File | Wired to CLI? | Purpose |
|------|--------------|---------|
| `core/run_command.py` | **Yes** — `rectify run-all` dispatches here | Production pipeline |
| `core/run_all_command.py` | **No** — experimental redesign, not active | Refactor-in-progress |

`run_command.py` dispatches to `_run_single_sample()` or `_run_multi_sample()` depending on whether `--manifest` is provided. **Always edit `run_command.py` for pipeline changes.**

DESeq2, motif discovery, and GO enrichment only run in multi-sample (`--manifest`) mode. Single-sample runs skip them by design (`run_deseq2 = n_samples > 1`).

---

## Module Map

| File | Role |
|------|------|
| `core/run_command.py` | **Active** `run-all` dispatcher (single + multi-sample) |
| `core/run_all_command.py` | Experimental redesign — NOT wired to CLI |
| `core/multi_aligner.py` | Runs all 5 aligners, builds per-aligner BAMs |
| `core/consensus.py` | Scores alignments, selects best per read |
| `core/bam_processor.py` | 3' end correction orchestration (parallel + streaming) |
| `core/atract_detector.py` | A-tract ambiguity detection |
| `core/indel_corrector.py` | Deletion artifact correction in poly(A) |
| `core/netseq_refiner.py` | NET-seq-guided ambiguity resolution |
| `core/splice_aware_5prime.py` | 5' junction rescue logic |
| `core/align_command.py` | `align` subcommand |
| `core/correct_command.py` | `correct` subcommand |
| `core/analyze_command.py` | `analyze` subcommand |
| `core/analyze/` | Analysis modules (DESeq2, APA shift, distribution, etc.) |
| `slurm.py` | CPU detection, thread limits, `$SCRATCH` utilities |
| `slurm_profiles/sherlock_larsms.yaml` | Standard CPU partition profile (use_scratch, streaming on) |

---

## Minimap2 Alignment: Junction Annotation

`rectify align` generates a junction BED from the GFF annotation and passes it to minimap2 via `--junc-bed`. The junction BED is cached as `annotation.junc.bed` in the sample output directory.

```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD \
    -t <threads> \
    --junc-bed <sample_dir>/annotation.junc.bed \
    --junc-bonus 9 \
    <genome.fsa.gz> <reads.fastq.gz>
```

Key flags:
- `-uf`: forward-strand-only (correct for direct RNA / cDNA sense reads)
- `-k14`: smaller k-mer for sensitivity on noisy nanopore reads
- `-G 5000`: max intron size (yeast introns < 1 kb)
- `--splice-flank=no`: disables GT-AG bonus (important for 3' end accuracy)
- `--MD`: required for indel artifact correction downstream

---

## Memory: Use `--streaming` for Large BAMs

**`rectify correct` without `--streaming` accumulates all results before writing.**

The default parallel path (`process_bam_file_parallel`) holds all corrected reads in RAM until the final TSV write. For a 7 GB BAM with 40M reads this requires ~30–40 GB RAM.

With `--streaming` (`process_bam_streaming`), reads are processed and written one chunk at a time (default 10,000 reads/chunk). Peak RSS drops to ~4–5 GB regardless of BAM size.

**Always use `--streaming` in SLURM jobs.** The `sherlock_larsms.yaml` profile sets `streaming: true` by default.

> **Known bug (fixed 2026-03-30):** `process_bam_streaming()` called `correct_read_3prime()` which returns `List[Dict]`, then passed the list directly to `stats.update_from_result()` which expected a single `Dict`. Fixed: use `.extend()` and iterate for stats.

---

## Multi-Sample Analysis: Streaming Pipeline

### Current state (rectify-beta `_run_multi_sample`)

The current manifest-mode pipeline (Stage 2) concatenates all per-sample corrected TSVs into a single combined file before running analysis. This works for small experiments but OOMs for large ones (e.g. 21 samples / 150M reads requires loading the entire dataset into RAM).

### Target state (ported from stable `rectify/`)

The stable repo implements a **two-pass streaming pipeline** that should be ported to rectify-beta. It completely avoids loading all samples simultaneously.

#### Pass 1 — position aggregation for clustering

Each per-sample TSV (or its position index, see below) is read sequentially. Only `chrom`, `strand`, `corrected_3prime` are loaded; all other columns are dropped immediately (**column pruning**). Positions are aggregated to unique (chrom, strand, pos) counts per sample and combined for clustering. Never loads more than one sample at a time.

#### Pass 2 — count matrix accumulation

For each sample, positions are streamed in 100k-row chunks and looked up against the cluster IntervalTree. Counts accumulate in a `defaultdict[cluster_id][sample_id]` — a ~1 MB structure regardless of dataset size. The final count matrix is built from this dict.

DESeq2, PCA, GO, and motifs then run on the count matrix (~10k × n_samples), which fits in a few MB.

**Peak RAM: O(clusters × samples) ≈ a few MB, regardless of read depth or sample count.**
Validated in stable repo on 21 samples / 150M reads on a 16 GB node.

#### Column pruning (Tier 1)

`load_corrected_positions()` drops all columns not needed downstream immediately after loading. Columns retained: `chrom`, `strand`, `corrected_position`, `sample`, `fraction` (if present), `alignment_start`, `alignment_end`. Everything else (`read_id`, `confidence`, `polya_length`, QC flags) is dropped on load.

The chunked-loading threshold should be lowered from 5 GB → 500 MB so aggregation kicks in for real datasets.

#### Position index (Tier 3) — even faster

`rectify correct` should write `corrected_3ends_index.bed.gz` alongside `corrected_3ends.tsv`. This is a pre-aggregated position count file:

```
chrom   corrected_3prime  strand  count
chrI    12836             +       1.0
chrI    41981             +       3.0
...
```

~300× smaller than the full per-read TSV. When an index is present, both passes use it instead of the full TSV — Pass 1 and Pass 2 become near-instant.

#### Back-filling indices on existing corrected TSVs

For samples corrected before the index feature is added:

```python
from rectify.core.bam_processor import write_position_index
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

samples = ['wt_rep1', 'wt_rep2', ...]
base = Path('results/')

def gen(s):
    tsv = str(base / s / 'corrected_3ends.tsv')
    write_position_index(tsv, tsv)

with ThreadPoolExecutor(max_workers=4) as ex:
    list(ex.map(gen, samples))
```

#### Reference condition matching

`--reference` is matched **case-insensitively** against the manifest `condition` column. `--reference wt` will match `WT`, `Wt`, etc. If no match is found, a warning is printed and the value is used as-is.

#### Manual post-correction analysis (streaming mode)

```bash
rectify analyze /dev/null \
    --manifest manifest.tsv \
    --genome genome.fsa.gz \
    --annotation genes.gff.gz \
    --reference wt \
    --run-deseq2 \
    --go-annotations go.tsv.gz \
    --run-motif \
    --threads 8 \
    -o results/combined/
```

> **Note:** Bedgraph and genomic distribution steps are skipped in manifest streaming mode (they require per-read alignment coordinates). Generate bedgraphs separately with `rectify export` if needed.

---

## Key Algorithms

### 5' Junction Rescue (`consensus.py: _rescue_5prime_softclip`)

Reads soft-clipped at the 5' end may have their clipped bases matching the upstream exon end (the aligner found the junction but couldn't extend fully into the upstream exon). Rescue checks whether the soft-clipped sequence matches the genomic sequence upstream of the nearest donor site within a candidate junction pool (annotated junctions ∪ all aligners' observed junctions for this read). Edit distance threshold: ≤20% mismatches over the clip length.

If rescued → **no penalty**. The aligner correctly identified the junction.
If not rescued → **penalty of -2 per clipped base**.

### Fair mapPacBio Scoring (`consensus.py: _get_effective_5prime_clip`)

mapPacBio forces mismatches/indels at splice junction boundaries instead of soft-clipping (it doesn't support soft-clipping at splice sites). This is functionally identical to a soft-clip but encoded differently in the CIGAR string. To avoid unfairly penalizing aligners that DO soft-clip (minimap2, gapmm2), the `effective_five_prime_clip` is computed as:

```
effective_clip = max(explicit_soft_clip, explicit_soft_clip + terminal_error_length)
```

where `terminal_error_length` is detected via a greedy sliding-window scan (window=8, threshold=40% error density) over the first 25 aligned bases from the 5' end.

Rescue is applied only if `explicit_soft_clip > 0` (we need the clip sequence to compare against). mapPacBio's forced-mismatch region is penalized but not rescuable via sequence comparison.

### Fair 3' End Scoring (`consensus.py: _get_effective_3prime_clip`)

Aligners should only soft-clip the poly(A) tail at the 3' end. Clipping or force-mismatching real exon sequence indicates the aligner stopped before the true 3' end. The `effective_three_prime_clip` detects terminal non-poly(A) errors using the same greedy sliding-window scan as the 5' version, but scanning inward from the 3' end and counting only non-A (+ strand) or non-T (− strand) mismatches as errors. Penalty: **-2/bp, capped at 10**.

This is complementary to the A-tract depth penalty (which catches going too far *into* downstream poly(A)).

### 3' End Correction Modules (`bam_processor.py: correct_read_3prime`)

**Module 1 — A-tract ambiguity (`atract_detector`):**
Counts downstream genomic A's (or T's for − strand) in a window after the apparent 3' end. If the poly(A) tail has extended INTO a downstream genomic A-tract, the apparent 3' end can be shifted too far downstream. Module 1 computes an *ambiguity range* used later by NET-seq refinement. It does NOT itself reposition the 3' end — it labels the uncertainty.

**Module 2C+2D — Comprehensive poly-A / homopolymer / indel correction (`indel_corrector.correct_3prime_position`):**
Applies three sub-strategies in order, taking the position that moves the 3' end furthest from the poly-A region:

1. **Poly-A boundary walk-back** (`find_polya_boundary`): walks backwards from the raw 3' end through the aligned region. Includes D-op positions (marked as `None` in the aligned pairs) and skips I-op bases; stops at the first non-A position where read and genome agree. Handles deletion artifacts (BDH1 case) and insertion artifacts — D ops shift the reference position downstream; I ops inflate poly-A length but don't shift reference coordinates. Both are detected and counted for QC.

2. **Soft-clip rescue at homopolymer boundaries** (`rescue_softclip_at_homopolymer`): handles the T-tract under-call case (HSP31). When the basecaller under-calls a T-tract, the aligner can't fit the true CPA base and soft-clips it along with the poly(A) tail. The apparent 3' end falls *inside* the T-tract. The fix: walk forward through the soft-clipped poly(A) bases, find the first non-A soft-clipped base, match it to the reference downstream of the T-tract, extend the 3' end there.

3. **Mismatch rescue inside homopolymers** (`rescue_mismatch_inside_homopolymer`): for global aligners that force-align rather than soft-clip; finds non-homopolymer bases that were force-mismatched and treats them as the true 3' end. Optional variant filter (`VariantAwareHomopolymerRescue`) skips positions where high read frequency suggests a true SNP.

**Module 2E — False junction filter (`false_junction_filter`):**
Detects poly(A)-artifact N ops (spurious splice junctions created when poly(A) tail aligns to downstream genomic A-tract) and removes them. Corrects 3' position to upstream exon boundary.

---

## Validation Suite

The following reads from the published S. cerevisiae direct RNA dataset (SRR32518284, Churchman lab) serve as gold-standard examples for each correction category. These should be verified in IGV after each pipeline change.

### 1. 3' End Correction: Indel Artifacts in Poly(A) Regions

Aligners introduce deletion artifacts when poly(A) tails align to downstream genomic A-tracts, shifting the apparent 3' end downstream. RECTIFY walks back through the A-tract, absorbing deletions to find the true CPA.

| Gene | Read Name | Artifact | Expected correction |
|------|-----------|----------|---------------------|
| BDH1 | `SRR32518284.2303904` | Deletion(s) in poly(A) region shift 3' end downstream | Walk-back absorbs D operations, stops at first non-A genome/read agreement |
| ECM1 | `SRR32518284.1550388` | Deletion artifact near CPA site | Walk-back corrects to annotated CPA |
| FLC2 | `SRR32518284.383494` | Deletion artifact in poly(A)-adjacent A-tract | Walk-back corrects 3' end |

**Verification:** In IGV, soft-clip boundary should be at the CIGAR D artifact. After correction, 3' end should be at the true annotated CPA. The `correction_applied` column in output TSV should contain `atract` or `indel`.

### 2. 5' End Correction: Splice Junction Soft-Clips

Long reads spanning the first intron of RPL20B are frequently soft-clipped at the 5' end because aligners cannot fully extend through the splice junction into the upstream exon. RECTIFY rescues these by matching the soft-clipped bases against the upstream exon sequence.

RPL20B is on the **minus strand** (`chrXV:900,112−901,419` in R64, − strand; YOR312C). The 5' end of the RNA is at the right (high coordinates). Soft-clips appear at the 5' (right) end. (Note: RPL20A paralog is YMR242C on chrXIII.)

Intron: chrXV:900,768−901,193 (1-based). Exon 1: 901,194−901,419. Exon 2: 900,250−900,767.

> **Note:** These reads are NOT present in the chrI-only test BAMs. RPL20B validation requires reads from chrXV.

| Read Name | Clip length | Expected outcome |
|-----------|-------------|-----------------|
| `SRR32518284.1771702` | ~10-15 bp | Rescued — matches upstream exon end, no penalty |
| `SRR32518284.722302` | ~10-15 bp | Rescued — matches upstream exon end, no penalty |
| `SRR32518284.2796136` | ~10-15 bp | Rescued — matches upstream exon end, no penalty |
| `SRR32518284.811112` | ~10-15 bp | Rescued — matches upstream exon end, no penalty |

**Verification:** In IGV, before correction the reads show a soft-clipped 5' end that does NOT span the 3' splice site. After consensus selection with rescue, the winning aligner should have the smallest `five_prime_softclip` (ideally 0) or `effective_five_prime_clip = 0`. The `consensus_stats.tsv` should show reduced 5' soft-clip rates for RPL20B reads.

**Plus-strand complement (EFB1):** EFB1 (`chrI`) shows the same artifact on the + strand, confirming strand-symmetric handling. In the 5-aligner validation run (1,110 reads: EFB1 + strand + RPL20B − strand), 5' clip rate dropped from 33% → 2.2% at EFB1 after rescue (v4 consensus scoring).

#### Module 2F: Exon 2 / 3'SS Truncation Rescue (post-consensus, second-pass)

**Status: Implemented** (`core/splice_aware_5prime.py: rescue_3ss_truncation`, wired into `bam_processor.py: correct_read_3prime()` as Module 2F).

Handles reads truncated or soft-clipped at the exon 2 / 3' splice site (3'SS) boundary that were not rescued during consensus scoring. Three cases in priority order:

1. **Soft-clip** — explicit S operation at 5' end: clip sequence compared against upstream exon end for each candidate junction (≤20% mismatch threshold).
2. **mapPacBio forced-mismatch** — no soft-clip but terminal alignment errors at 5' end: errors extracted via greedy sliding-window scan (same as `_get_effective_5prime_clip`) and compared against upstream exon.
3. **Proximity-only (zero-clip truncation)** — alignment starts at or within `junction_proximity_bp=10 bp` of a known 3'SS with no clipped or mismatched bases: recorded as `rescue_type='proximity'` without changing the attributed 5' position (no sequence evidence available).

**Junction universe:**
```
annotated 3'SS from GFF  ∪  novel 3'SS from consensus BAM CIGAR
```

Sets `five_prime_rescued=True` and updates `five_prime_position` for cases 1 & 2. Case 3 sets `rescue_type='proximity'` without position change.

**Expected impact for RPL20B:** Nearly all reads landing at the exon 2/3'SS boundary should be rescuable via case 1 or 2, since exon 1 / exon 2 sequences are unambiguous. Requires full-genome run to verify (RPL20B is chrXIII).

### 3. 3' False Junction Handling

Poly(A) tails can create spurious N (skip) operations in CIGAR strings when the aligner bridges from the poly(A) tail into a downstream genomic A-tract. These phantom junctions must be detected and removed.

| Gene | Read Name | Artifact | Expected correction |
|------|-----------|----------|---------------------|
| MCH5 | `SRR32518284.711345` | False 3' splice junction (N in CIGAR) from poly(A) tail aligning to downstream A-tract | Walk-back discards N operation; true 3' end upstream of false junction |

**Verification:** Raw CIGAR should contain an `N` operation near the 3' end. After correction, `corrected_3prime` should be upstream of the false junction start. The `correction_applied` column should indicate the N was discarded. No entry for MCH5 reads should appear in the corrected junction calls.

### 4. Soft-Clip Rescue at Homopolymer Boundaries

Nanopore basecallers systematically under-call homopolymer runs (e.g., calling 8 U's instead of 10). When this happens at CPA sites with upstream T-tracts, the aligner soft-clips non-T bases instead of placing them correctly, causing a shortened apparent 3' end.

| Gene | Read Name | Artifact | Expected correction |
|------|-----------|----------|---------------------|
| HSP31 | `SRR32518284.2969590` | Poly-U under-call at CPA site; aligner soft-clips non-T base at homopolymer boundary | 3' end extended past the under-called T-tract to match true CPA |

**Mechanism:** The basecaller under-calls the T-tract (e.g., calls 9 T's when the reference has 11). The aligner reaches the first non-T base (the CPA base, e.g. C) but cannot fit it because the T-tract isn't fully traversed, so it soft-clips that base along with the poly(A) tail. The apparent 3' end falls inside the T-tract. `rescue_softclip_at_homopolymer()` in Module 2C walks through the soft-clipped poly(A) bases (A, A, A, …) and finds the first non-A soft-clipped base, matches it to the reference downstream of the T-tract, and extends the 3' end there. `correction_applied` will contain `softclip_homopolymer_rescue`.

**Note:** Module 1 (A-tract ambiguity) detects a *different* problem — the poly(A) tail aligning INTO a downstream genomic A-tract, pushing the end too far downstream. The T-tract case pushes the end too far upstream. These are complementary failures.

**Verification:** Raw `reference_end` should be inside the annotated T-tract, upstream of the true HSP31 CPA. Soft-clip sequence should contain the CPA base followed by poly(A). After correction, `corrected_3prime` should match the annotated CPA. `correction_applied` should contain `softclip_homopolymer_rescue`.

---

## Known Issues & Active Fixes

### Memory: 5-Aligner Mode Requires >48GB

**Status:** Under investigation (2026-03-31)

The 5-aligner consensus step (minimap2 + mapPacBio + gapmm2 + uLTRA + deSALT) consumes approximately 48GB of RAM for typical yeast chrI datasets (13K–250K reads). After alignment completes, residual memory is not fully released before the correction step starts, causing OOM kills at `--mem=48G`.

**Workaround:** Increase `--mem=96G` in SLURM scripts.

**Fix needed:** Explicit `del` / garbage collection of large alignment data structures in `run_all_command.py` between alignment and correction steps, and/or BAM streaming instead of loading all 5 aligner outputs simultaneously.

### deSALT: Annotation-Guided Mode Disabled

**Status:** Permanent workaround (2026-03-30)

deSALT's `-G` flag (GTF/GFF annotation) causes a SIGSEGV in its GTF parser when used with yeast annotation files (both gzipped GFF3 and unzipped GTF trigger the crash). The crash occurs in deSALT binary regardless of annotation format.

**Workaround:** `-G` flag removed from `multi_aligner.py`. deSALT runs in de novo splice detection mode only. This is sufficient — deSALT's de novo mode correctly aligned and won 35/1,110 reads (3.15%) in the 5-aligner junction rescue validation.

### ProcessPoolExecutor (forkserver) + Spawn Subprocess Isolation

**Status:** Fixed (2026-04-01)

`multiprocessing.Pool` (fork-based) deadlocked when called after alignment subprocesses ran. `ThreadPoolExecutor` was used as an intermediate fix, but htslib calls in the main process (read counting, `get_processing_regions`) left residual malloc metadata corruption that caused `free(): invalid next size` SIGABRT in worker result-handler threads ~50–130 regions later (exit code 134).

**Final fix (both issues):**
1. `_count_reads_and_get_regions()` now runs in a `multiprocessing.get_context('spawn').Pool(1)` subprocess — main process heap is never touched by htslib, eliminating the corruption source.
2. `process_bam_file_parallel()` uses `ProcessPoolExecutor(forkserver)` with `_worker_initializer` that loads genome, netseq, and opens the BAM once per worker process and stores them in module-level globals (`_WORKER_GENOME`, `_WORKER_BAM`, etc.).
3. `_process_region_worker` reads exclusively from module-level globals — no large data structures are passed through the executor's IPC path.

### `run_3prime_distribution_analysis` / `run_transcript_body_distribution_analysis`

**Status:** Fixed (2026-03-31)

Both functions were removed or renamed in `analyze/genomic_distribution.py` but stale references remained in `analyze_command.py` and `analyze/__init__.py`, causing `ImportError` at CLI startup.

**Fix:** Updated call sites to use `run_genomic_distribution_analysis`; removed `run_transcript_body_distribution_analysis` (function no longer exists).

---

## Test Coverage

```
284 passing, 0 failing (as of 2026-04-01)
```

Tests are run from `/oak/stanford/groups/larsms/Users/kevinroy/software/rectify-beta`:

```bash
PYTHON="/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify/bin/python3"
$PYTHON -m pytest tests/ -q
```

**Previously failing tests now fixed:**
- `test_atract.py` — 3 off-by-one assertions (position convention: last-aligned-base, not first-unaligned)
- `test_parallel_processing.py` — 2 tests updated for `_process_region_worker` module-globals interface
- `test_splice_junction.py` — import path was stale; confirmed passing (51/51)
- `test_analyze.py` PCA tests — all passing; scikit-learn present in env

---

## SLURM Batch Configuration

### Current Batch Script
`/oak/stanford/groups/larsms/Users/kevinroy/projects/roadblocks/rectify_beta_output/chrI_test_new/run_rectify_beta_chrI_new.sh`

**Key features:**
- Stages input FASTQ from OAK → `$SCRATCH` before running (75 GB/s local SSD vs. ~1-5 GB/s OAK NFS)
- Runs all I/O (5-aligner BAMs, consensus BAM, correction outputs) on `$SCRATCH`
- Copies final outputs back to OAK via `rsync` after completion
- Cleans up `$SCRATCH` regardless of exit code

**Current resource allocation:**
- `--mem=48G` ← **needs to be increased to 96G** (5-aligner mode hits 48G limit)
- `--cpus-per-task=16`
- `--time=4:00:00`
- `--partition=larsms,owners`

### Conda Environment
```bash
PYTHON="/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify/bin/python3"
export PATH="/home/groups/larsms/users/kevinroy/anaconda3/envs/rectify/bin:$PATH"
```

Includes: pysam, numpy, tqdm, uLTRA (`/home/users/kevinroy/.local/bin/uLTRA`), namfinder (required by uLTRA)

---

## CLI Reference

```bash
# Single sample (bundled yeast data)
rectify run-all reads.fastq.gz --Scer -o results/

# Full 5-aligner mode with junction aligners
rectify run-all reads.fastq.gz --Scer \
    --parallel-aligners \
    --junction-aligners uLTRA deSALT \
    --ultra-path /home/users/kevinroy/.local/bin/uLTRA \
    --filter-spikein ENO2 \
    --threads 16 \
    -o results/

# Multi-sample via manifest
rectify run-all --manifest samples.tsv --Scer \
    --parallel-aligners --junction-aligners uLTRA deSALT \
    --ultra-path /home/users/kevinroy/.local/bin/uLTRA \
    --filter-spikein ENO2 --threads 16 \
    -o results/
```

---

## Scoring Summary (consensus.py)

| Penalty | Value | Notes |
|---------|-------|-------|
| 5' unrescued clip | −2/bp | Applied to `effective_five_prime_clip` |
| 5' rescued clip | 0 | Soft-clipped bases match upstream exon |
| 3' A-tract depth | −1/downstream A | Capped at −10 |
| 3' non-polyA terminal | −2/bp | Non-A/T terminal errors, capped at −10 |
| Canonical GT-AG | tiebreaker only | Not counted in primary score |
| Annotated junction | tiebreaker only | Not counted in primary score |

**Tiebreaker order:**
1. Fewer unrescued 5' clips (uses `effective_five_prime_clip`)
2. Majority 3' position vote across aligners
3. Canonical splice site count
4. Annotated junction count

---

## Next Steps

### Immediate (blocking batch validation)
- [x] Increase `--mem=96G` in SLURM batch script and resubmit *(done 2026-03-31)*
- [x] Fix htslib heap corruption (SIGABRT exit 134) — spawn subprocess isolation *(done 2026-04-01)*
- [x] Memory optimization: `gc.collect()` between alignment and correction in `run_command.py` *(done 2026-04-01)*

### Short-term (validation)
- [ ] Run full 12-sample chrI validation to completion (job 20209831 running)
- [ ] Validate all 4 correction categories using the reads listed in the Validation Suite above
- [x] Fix 8 failing unit tests *(done 2026-04-01 — 284/284 passing)*
- [x] Fix `tests/test_splice_junction.py` import error *(was never broken — 51/51 passing)*

### Medium-term (port from stable `rectify/`)
- [x] **Two-pass streaming manifest pipeline** — `_run_analysis_manifest()` added; `_run_multi_sample()` now writes manifest and calls manifest mode *(done 2026-04-01)*
- [x] **Position index (`corrected_3ends_index.bed.gz`)** — `write_position_index()` ported and wired in all three `bam_processor.py` output paths *(done 2026-04-01)*
- [x] **Column pruning in `load_corrected_positions()`** — drops non-essential columns after load; chunked-loading threshold lowered to 500 MB *(done 2026-04-01)*
- [x] **Case-insensitive `--reference` matching** — case-insensitive lookup against available conditions with fallback warning *(done 2026-04-01)*
- [x] **`--streaming` as default in SLURM** — `sherlock_larsms.yaml` sets `streaming: true` *(done 2026-04-01)*
- [x] **False junction filter (`false_junction_filter.py`)** — poly(A)-artifact N operations detected and removed; `n_false_junctions` + `five_prime_rescued` columns in TSV output *(done 2026-04-01)*
- [ ] **Exon 2 / 3'SS truncation rescue** — post-consensus module in `splice_aware_5prime.py`; rescues reads truncated at the 3' splice site boundary even when no soft-clip sequence is present; uses annotated junctions ∪ novel junctions from all aligners' first pass; see Validation Suite §2 for design detail. Priority: **HIGH** — large fraction of RPL20B reads affected.
- [ ] Extend validation from chrI → full genome
- [ ] Profile consensus scoring performance (`_get_effective_5prime_clip` and `_get_effective_3prime_clip` both call `read.get_aligned_pairs()` which is O(read length); may dominate at high coverage)
- [ ] Add validation suite reads to unit tests as fixtures

### Publication prep
- [ ] Benchmark 5-aligner consensus vs. 3-aligner and single-aligner on held-out dataset
- [ ] Characterize false-positive junction rescue rate
- [ ] Document scoring parameter sensitivity (window size, error threshold)

---

## Validation Run History

| Job | Date | Config | Result | Notes |
|-----|------|--------|--------|-------|
| Interactive | 2026-03-29 | 1,110-read rescue test (EFB1 + RPL20B) | ✅ v5 scoring: mapPacBio win rate 36.3%→32.3%; EFB1 5'-clip 33%→2.2% | First 5-aligner test |
| 19947271 | 2026-03-30 | 12×chrI, 5-aligner, OAK I/O | ❌ Deadlock in correction step | `multiprocessing.Pool` fork deadlock |
| 20009640 | 2026-03-30 | 12×chrI, 5-aligner, $SCRATCH | ❌ ImportError at startup | `run_3prime_distribution_analysis` missing |
| 20066871 | 2026-03-31 | 12×chrI, 5-aligner, $SCRATCH, ThreadPoolExecutor | ❌ OOM (48GB) | 5-aligner mode peaks at ~50GB |
| 20069496 | 2026-03-31 | Same + `--mem=96G` | In progress | — |
