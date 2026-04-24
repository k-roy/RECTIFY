# RECTIFY — Developer & Agent Context

This file documents non-obvious design decisions, known pitfalls, and
patterns for future developers and AI agents working on this codebase.

---

## Architecture: `run-all` dispatcher

`rectify run-all` dispatches to `rectify/core/run_command.py`, which calls
`_run_single_sample()` or `_run_multi_sample()` based on whether `--manifest`
is provided.

**DESeq2, motif discovery, and GO enrichment only run in multi-sample mode.**
Single-sample runs skip them by design (`run_deseq2 = n_samples > 1`).

### DRS workflow: `--drs` flag

For Dorado direct RNA-seq BAMs, pass `--drs` to `run-all`:

```bash
rectify run-all sample_dorado.bam --drs --Scer -o results/sample/
```

This automatically inserts two DRS-specific steps around the standard pipeline:

| Step | Action | Implementation |
|------|--------|----------------|
| **0** | Poly(A)+adapter pre-trimming | `trim_drs_bam_polya()` → `samtools fastq -T pt` |
| 1 | Multi-aligner alignment | (normal, on trimmed FASTQ) |
| 2 | Correction | (normal) |
| 3 | Analysis | (normal) |
| **4** | Restore poly(A) as soft-clips | `restore_polya_softclips()` |

**Key implementation details:**
- DRS Steps 0 and 4 run in **both** single-sample (`_run_single_sample`) and multi-sample manifest (`_process_one_sample`) paths
- Single-sample: trimmed FASTQ written to `$SCRATCH/drs_trim/`; trim metadata parquet to Oak output dir (persistent, needed for Step 4)
- Multi-sample: trimmed FASTQ and temp BAM written to `<sample_output>/drs_trim/`; metadata parquet to `<sample_output>/`
- `--drs` has no effect on FASTQ inputs (assumed already trimmed)
- `softclipped_bam` path for Step 4 is inferred from the trimmed FASTQ stem (stripping `_trimmed`/`.fastq`/`.gz`): `{stem}.rectified_pA_tail_trimmed.bam`
- Single-sample Steps 0 and 4 are recorded in the provenance tracker as `drs_trim` and `drs_restore`
- Multi-sample Step 4 restore is non-fatal: failure is logged but the sample still returns success

---

## HPC I/O: always stage through $SCRATCH

**Never run correction I/O directly on shared NFS storage in a SLURM array job.**

NFS-backed shared storage under concurrent load causes severe I/O contention and
inflates wall time 2–3×. On HPC clusters, `$SCRATCH` provides high-bandwidth
local storage and is node-local to each task.

### Pattern (Python pipeline)

`run_command._run_single_sample()` handles this automatically:

```python
from rectify.slurm import make_job_scratch_dir, sync_outputs

scratch = make_job_scratch_dir('rectify_single')  # None if no scratch
work_dir = scratch if scratch else output_dir

# ... all alignment + correction + analysis writes to work_dir ...

if scratch:
    sync_outputs(scratch, output_dir)  # rsync everything back
    shutil.rmtree(scratch)
```

### Pattern (generated SLURM scripts)

`batch_command.generate_slurm_scripts()` inserts staging blocks when
`use_scratch=True` in the profile (default in `hpc_cpu.yaml`):

- BAM inputs are `cp`'d to `$SCRATCH` before `rectify correct`
- FASTQ inputs are left on the source filesystem (read sequentially once, no benefit to staging)
- All outputs are `rsync`'d back to the output directory after correction

### $SCRATCH caveats

- Auto-purged after **90 days** on most HPC clusters — final outputs must land on persistent storage
- Not shared across nodes — do not use `$L_SCRATCH` for multi-node jobs
- `make_job_scratch_dir()` uses `$SCRATCH` → `$SLURM_TMPDIR` → `$TMPDIR` priority

---

## Memory: use `--streaming` for large BAMs

**`rectify correct` without `--streaming` loads all reads into RAM.**

The default parallel path (`process_bam_file_parallel`) accumulates all
results before writing. For a 7 GB BAM with 40M reads this requires
~30–40 GB RAM.

With `--streaming` (`process_bam_streaming`), reads are processed and
written one chunk at a time (default 10,000 reads/chunk). Peak RSS drops
to ~4–5 GB regardless of BAM size.

**Always use `--streaming` in SLURM jobs.** The `hpc_cpu.yaml`
profile sets `streaming: true` by default.

---

## Thread limits: set before importing numpy

**SLURM will ban accounts that spawn more processes than allocated CPUs.**

`set_thread_limits()` in `slurm.py` must be called before any import of
numpy, sklearn, or pydeseq2. These libraries auto-spawn threads via OpenMP,
OpenBLAS, MKL, and joblib's loky backend.

The loky backend (used by sklearn/pydeseq2) ignores `JOBLIB_WORKERS` —
you **must** set `LOKY_MAX_CPU_COUNT`.

```bash
# In SLURM scripts — set BEFORE running Python
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LOKY_MAX_CPU_COUNT=$SLURM_CPUS_PER_TASK   # critical for pydeseq2
```

---

## Protocol flag: `--dT-primed-cDNA`

**Nanopore direct RNA-seq (DRS):** sequences the entire molecule from the
adapter through the poly-A tail into the RNA body.  The poly-A tail IS in the
read.  This is the default mode — no protocol flag needed.

**dT-primed cDNA (QuantSeq, oligo-dT cDNA):** the oligo-dT primer hybridizes
at the start of the poly-A tract; sequencing starts at the first non-A base.
The poly-A tail is NOT in the read.  Pass `--dT-primed-cDNA` to enable
Module AG (AG mispriming detection), which is a cDNA-synthesis artifact caused
by reverse-transcriptase slippage at internal AG runs.  This artifact does not
occur in DRS and is disabled by default.

`--polya-sequenced` is a deprecated alias for `--dT-primed-cDNA` and will be
removed in a future version.

### Module activation by protocol

| Module | DRS (default) | dT-primed cDNA (`--dT-primed-cDNA`) |
|--------|:---:|:---:|
| A-tract ambiguity (2E) | ✓ | ✓ |
| Poly(A) trimming (2B) | ✓ | ✓ |
| Indel correction (2C) | ✓ (if MD tags present) | ✓ (if MD tags present) |
| Soft-clip rescue (2G) | ✓ | ✓ |
| AG mispriming | ✗ | ✓ |
| Variant-aware rescue | ✓ | ✓ |

Poly(A) trimming and indel correction are always enabled because the poly-A
tail is present in the read sequence for both protocols (DRS: native RNA;
dT-primed cDNA: the primer hybridizes *at* the poly-A, so poly-A bases appear
at the read start after alignment).  AG mispriming is the only module that is
protocol-specific.

Module 2G (`rescue_softclip_at_homopolymer`) runs for all protocols.  It
detects 3' soft-clips ≥ 3 bp adjacent to a genomic homopolymer and extends
the 3' end outward.  It takes priority over poly-A walk-back (Module 2E) to
prevent opposite-direction corrections from cancelling.

---

## DRS-trimmed FASTQ: strip auxiliary tags before alignment

`trim_drs_bam_polya()` + `samtools fastq -T pt` produces FASTQ with tab-separated
auxiliary tags embedded in the read header, e.g.:

```
@6c606d1b-2310-4292-a285-d519fbd52502	pt:i:25
```

**Problem:** gapmm2 runs minimap2 to produce a PAF, then looks up each query
sequence by name using `query_idx.seq(paf[0])`. minimap2 strips everything after
the first whitespace when building the read name in the PAF, so `paf[0]` = just
the UUID. But mappy's sequence index retains the full header including the tab
and `pt:i:N` suffix, so `query_idx.seq(UUID)` returns `None` →
`TypeError: object of type 'NoneType' has no len()`.

**Fix:** Strip auxiliary tags from FASTQ headers before passing to gapmm2 (or any
aligner that round-trips through PAF):

```python
clean_header = "@" + read_id + "\n"  # keep UUID only
```

**mapPacBio** has a different problem: it embeds `pt:i:N` into the READ NAME of
the aligned BAM.  The exact separator depends on processing stage:

- **Direct mapPacBio output** (before `samtools sort`): space-separated,
  `UUID pt:i:25`.  This is the format in live production data handled by
  `consensus.py` and `corrected_consensus.py`.
- **After `samtools sort/merge`**: samtools converts the space to an
  underscore (BAM spec forbids spaces in QNAME), yielding `UUID_pt:i:25`.
  This is the format found in merged dev-run BAMs and the sorted validation
  aligner BAMs.

Use whichever form is present:

```python
# Direct mapPacBio output (pre-sort):
if " pt:i:" in name:
    clean = name.split(" pt:i:")[0]

# Merged / sorted BAMs:
if "_pt:i:" in name:
    clean = name.split("_pt:i:")[0]
```

This has been encountered multiple times when extracting validation reads from
DRS-trimmed dev-run outputs.

---

## Minimap2 alignment: junction annotation

The `rectify align` command generates a junction BED from the GFF annotation
and passes it to minimap2 via `--junc-bed`. This improves splice junction
accuracy for long-read RNA-seq.

The junction BED is cached as `annotation.junc.bed` in the sample output
directory. The exact minimap2 command:

```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD
    -t <threads>
    --junc-bed <sample_dir>/annotation.junc.bed
    --junc-bonus 9
    <genome.fsa.gz> <reads.fastq.gz>
```

Key flags:
- `-uf`: forward-strand-only (correct for direct RNA / cDNA sense reads)
- `-k14`: smaller k-mer for sensitivity on noisy nanopore reads
- `-G 5000`: max intron size (yeast introns < 1 kb)
- `--splice-flank=no`: disables GT-AG bonus (important for 3' end accuracy)
- `--MD`: required for indel artifact correction downstream

---

## Multi-sample analysis: memory-efficient streaming pipeline

For DESeq2 to run, you must use `--manifest` with a TSV containing
`sample_id`, `path`, and `condition` columns.

**Do NOT pre-combine corrected TSVs into a single file.** The old approach
(concat all per-sample TSVs → one file → `rectify analyze`) requires loading
the entire dataset into RAM and OOMs for large experiments. Use manifest mode
instead.

### The correct pipeline (streaming, low RAM)

`rectify run-all --manifest` automatically uses manifest mode end-to-end.
For manual post-correction analysis:

```bash
rectify analyze /dev/null \
    --manifest manifest.tsv \       # columns: sample_id, path, condition
    --genome genome.fsa.gz \
    --annotation genes.gff.gz \
    --reference wt \                # case-insensitive match against condition column
    --run-deseq2 \
    --go-annotations go.tsv.gz \
    --run-motif \
    --threads 8 \
    -o results/combined/
```

### How manifest mode works (two-pass streaming)

**Pass 1 — position aggregation for clustering:**
Each per-sample TSV (or its index, see below) is read sequentially.
Only `chrom`, `strand`, `corrected_3prime` are loaded; all other columns
are dropped. Positions are aggregated to unique (chrom, strand, pos) counts
per sample and combined for clustering. Never loads >1 sample at a time.

**Pass 2 — count matrix accumulation:**
For each sample, positions are streamed in 100k-row chunks and looked up
against the cluster IntervalTree. Counts accumulate in a
`defaultdict[cluster_id][sample_id]` — a ~1 MB dict regardless of dataset
size. The final count matrix is built from this dict.

**DESeq2, PCA, GO, motifs** then run on the count matrix (~10k × n_samples),
which fits in a few MB.

**Note:** Bedgraph and genomic distribution steps are skipped in manifest
mode (they require per-read alignment coordinates). Generate bedgraphs
separately with `rectify export` if needed.

Peak RAM: O(clusters × samples) ≈ a few MB, regardless of read depth or
sample count. Validated on 21 samples / 150M reads: runs on a 16 GB node.

### Position index (Tier 3): even faster

`rectify correct` now writes `corrected_3ends_index.bed.gz` alongside
`corrected_3ends.tsv`. This is a pre-aggregated position count file:

```
chrom   corrected_3prime  strand  count
chrI    12836             +       1.0
chrI    41981             +       3.0
...
```

~300× smaller than the full per-read TSV. When an index is present, both
passes of manifest mode use it instead of the full TSV — Pass 1 and Pass 2
become near-instant (index files are tiny and already aggregated).

**Generating indices for existing corrected TSVs** (those corrected before
this feature was added):

```python
from rectify.core.bam_processor import write_position_index
from concurrent.futures import ThreadPoolExecutor

samples = ['wt_rep1', 'wt_rep2', ...]
base = Path('results/')

def gen(s):
    tsv = str(base / s / 'corrected_3ends.tsv')
    write_position_index(tsv, tsv)

with ThreadPoolExecutor(max_workers=4) as ex:
    list(ex.map(gen, samples))
```

### Column pruning (Tier 1)

`load_corrected_positions()` now drops all columns not needed downstream
immediately after loading. Columns retained: `chrom`, `strand`,
`corrected_position`, `sample`, `fraction` (if present), `alignment_start`,
`alignment_end`. Everything else (`read_id`, `confidence`, `polya_length`,
QC flags, etc.) is dropped on load.

The chunked-loading threshold was also lowered from 5 GB → 500 MB, so
aggregation now kicks in for most real datasets.

### Reference condition matching

`--reference` is matched **case-insensitively** against the manifest's
`condition` column. `--reference wt` will match `WT`, `Wt`, etc. If no
match is found, a warning is printed and the value is used as-is.

---

## Coordinate conventions

All coordinates are **0-based, half-open** (consistent with pysam / BED):

| Strand | 5' end (TSS) | 3' end (CPA) |
|---|---|---|
| `+` | `reference_start` | `reference_end - 1` |
| `-` | `reference_end - 1` | `reference_start` |

GFF files use 1-based coordinates — subtract 1 when loading.

---

## Cat3 5' rescue: local alignment for exon CIGAR (v2.8.0)

`rescue_3ss_truncation()` in `splice_aware_5prime.py` now calls `align_clip_to_exon()`
from `local_aligner.py` to compute a proper M/I/D CIGAR for the exon segment instead
of emitting a flat `nM` block.

**Result dict**: `five_prime_exon_cigar` (SAM CIGAR string, e.g. `"8M1D3M"`).
Stored in `corrected_3ends.tsv` as the `five_prime_exon_cigar` column (v2.8.0).
`bam_writer.extend_read_5prime_for_junction_rescue()` uses it when writing Cat3 reads.

**Aligner design** (`local_aligner.py`):
- `_align_right_anchored(query, ref)` — free prefix in ref, right end fixed.
  Plus strand: clip must end at `intron_start`.
- `_align_left_anchored(query, ref)` — left end fixed, free suffix.
  Minus strand: clip must start at `intron_end`.
- **Affine gap (Gotoh 1982)**: match=+2, mismatch=-4, gap_open=-4, gap_extend=-1.
  Total cost for a gap of length k: -4 + k×(-1). Three-matrix DP (H/D/I states).
  Affine gap prevents the "staircase" artifact where many isolated 1-base deletions
  outscore a single consolidated deletion (e.g., `3D` scores -7; three separate `1D`
  ops separated by matching bases score ≤ -9 total).
- Buffer = `max_indel=5` bp added to each side of the expected exon window.

**Chimeric stitch fix (v2.8.0)**: `chimeric_consensus.build_chimeric_cigar()` now
uses `D` (deletion, op=2) instead of `N` (intron skip, op=3) for reference gaps
≤ 10 bp at segment boundaries. Larger gaps are still classified as `N`.

---

## File layout

```
rectify/
├── rectify/
│   ├── slurm.py              # CPU detection, thread limits, scratch utilities
│   ├── slurm_profiles/
│   │   ├── hpc_cpu.yaml          # Standard CPU partition (use_scratch, streaming on)
│   │   └── hpc_gpu.yaml          # GPU partition
│   ├── core/
│   │   ├── run_command.py    # Active run-all dispatcher (single + multi-sample)
│   │   ├── batch_command.py  # Generates SLURM array scripts; reads profile YAMLs
│   │   ├── correct_command.py
│   │   ├── bam_processor.py  # process_bam_file_parallel + process_bam_streaming
│   │   ├── bam_writer.py     # CIGAR surgery + BAM output writers (Cat3/6)
│   │   ├── position_index.py # Writes corrected_3ends_index.bed.gz
│   │   ├── local_aligner.py          # Semi-global NW for Cat3 exon CIGAR (v2.8.0)
│   │   ├── multi_aligner.py          # minimap2/mapPacBio/gapmm2/uLTRA/deSALT wrappers
│   │   ├── split_command.py          # rectify split — FASTQ chunker for SLURM arrays
│   │   ├── consensus_command.py      # rectify consensus — aligner selection on pre-built BAMs
│   │   └── install_aligners_command.py  # rectify install-aligners — download/compile aligners
│   ├── data/                             # Bundled reference data — use --Scer to activate
│   │   ├── S288C_reference_sequence_R64-5-1_20240529.fsa     # S. cerevisiae genome FASTA
│   │   ├── S288C_reference_sequence_R64-5-1_20240529.fsa.fai # samtools FASTA index
│   │   ├── S288C_reference_sequence_R64-5-1_20240529.pkl     # pre-loaded genome dict (fast)
│   │   ├── saccharomyces_cerevisiae_R64-5-1_20240529.gff     # S. cerevisiae gene annotation
│   │   ├── saccharomyces_cerevisiae_R64-5-1_20240529.gtf     # GTF version of same annotation
│   │   ├── saccharomyces_cerevisiae_R64-5-1_20240529.junc.bed  # pre-built junction BED for minimap2 --junc-bed
│   │   ├── saccharomyces_cerevisiae_atract_netseq.tsv.gz     # A-tract NET-seq signal (yeast)
│   │   ├── saccharomyces_cerevisiae_netseq_pan.tsv.gz        # Pan NET-seq table
│   │   ├── saccharomyces_cerevisiae_netseq_wt.tsv.gz         # WT-only NET-seq table
│   │   ├── models/                                            # Poly-A model weights
│   │   ├── motif_databases/                                   # JASPAR-format motif databases
│   │   ├── validation/                                        # Bundled validation reads + TSVs
│   │   ├── genomes/saccharomyces_cerevisiae/                  # GFF.GZ copy (for rectify prescan)
│   │   │   └── saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
│   │   └── bin/
│   │       └── linux_x86_64/
│   │           └── deSALT            # Vendored deSALT v1.5.6 binary (773 KB)
│   └── ...
└── CLAUDE.md                         # This file
```

### Using bundled data directly

When `--Scer` (or `--organism saccharomyces_cerevisiae`) is passed, these paths
are resolved automatically.  To reference them explicitly (e.g., in `rectify prescan`
or test scripts):

```python
import rectify
from pathlib import Path
DATA = Path(rectify.__file__).parent / 'data'

GENOME     = DATA / 'S288C_reference_sequence_R64-5-1_20240529.fsa'
ANNOTATION = DATA / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff'
JUNC_BED   = DATA / 'saccharomyces_cerevisiae_R64-5-1_20240529.junc.bed'
```

Or from a shell script:
```bash
RECTIFY_DATA=$(python -c "import rectify; from pathlib import Path; print(Path(rectify.__file__).parent / 'data')")
GENOME=$RECTIFY_DATA/S288C_reference_sequence_R64-5-1_20240529.fsa
ANNOTATION=$RECTIFY_DATA/saccharomyces_cerevisiae_R64-5-1_20240529.gff
```

---

## Empirical Penalty Tables (`--junction-penalty-table`)

Production tables are in
`common/scripts/nanopore/error_profile_20260422/` (generated 2026-04-22).

### Format: AT/CG base-class split

`penalty_scores.tsv` groups reference bases into two classes: **AT** and **CG**.
This matters because AT runs have ~10–20% higher deletion rates than CG runs at the
same HP length (Nanopore pore-ratcheting asymmetry).

Columns: `op_type`, `base_class`, `hp_length`, `rate_mean`, `count_total`,
`penalty_score`, `low_count`.

Key numbers:

| op_type | base_class | HP=1 penalty | HP=4 penalty | HP=8 penalty |
|---------|------------|:------------:|:------------:|:------------:|
| D       | AT         | 0.37         | 0.33         | 0.28         |
| D       | CG         | 0.58         | 0.42         | 0.37         |
| X       | AT         | 1.00 (ref)   | 1.55         | 10.0 (cap)   |
| X       | CG         | 1.00 (ref)   | 2.56         | 10.0 (cap)   |

Insertions: all `low_count=True`, treat all as 10.0.

### Usage

```bash
rectify correct \
    --junction-penalty-table .../error_profile_20260422/penalty_scores.tsv \
    --str-penalty-table      .../error_profile_20260422/str_penalty_scores.tsv \
    ...
```

### Regeneration

See `common/scripts/nanopore/PENALTY_TABLE.md` for the full regeneration command,
known caveats (deSALT failures on chunks 2/11/15), and diagnostic plot instructions.

### Design notes

- Isotonic smoothing ensures deletion penalties are monotone non-decreasing with HP length.
- `low_count=True` rows (count < 100) should not be used as reliable estimates.
- Tables are S. cerevisiae R10.4.1-specific — regenerate for other organisms/chemistries.
- The `_CANONICAL_HP_PRIOR = 0.5` in `junction_refiner.py` was calibrated so that one
  HP deletion (del_cost ≈ 0.37–0.58 at HP=1) gives canonical junctions the tie-breaking
  advantage. This works correctly with both default heuristics and empirical tables.

---

## Audit history & remaining open issues

Two rounds of systematic codebase audits were completed 2026-04-08.
All findings are tracked in `docs/BUGS_TO_FIX.md`.

**v2.7.8 (2026-04-09):** Fixed 35 Round 2 findings including all confirmed
CRITICALs and HIGHs. Most impactful fixes:
- 6 minus-strand correctness bugs (junction validation donor/acceptor swap,
  junction aggregate strand handling, NET-seq double-subtraction, chimeric
  consensus CIGAR anchoring, 5' splice rescue window, full-length classifier)
- Shell injection in generated SLURM scripts (`shlex.quote()` added)
- `streaming: true` / `use_scratch: true` silently ignored by profile apply
- `NCBI_TO_CHR` reversal overwrote `chrI` with `I`
- `mapPacBio` timeout/deadlock/return-code issues

**v2.7.9 (2026-04-09):** Fixed 4 additional findings from validation audit
(N-op absorption in CIGAR surgery, `five_prime_rescued` TSV gap, chimeric
walkback exemption, NET-seq flag guard).

**v2.8.0 (2026-04-11):** Cat3 5' junction rescue now uses semi-global NW local
alignment (`local_aligner.py`) for the exon CIGAR instead of a flat `nM` block.
New `five_prime_exon_cigar` column in `corrected_3ends.tsv`. Chimeric stitch
gaps ≤ 10 bp use CIGAR op `D` instead of `N`. Local aligner upgraded to affine
gap scoring (Gotoh 1982, gap_open=-4, gap_extend=-1) to prevent staircase artifacts
where isolated 1-base deletions outscore a single consolidated deletion.

**v2.8.1 (2026-04-11):** New Module 2G: standalone soft-clip rescue at
homopolymer boundaries (`rescue_softclip_at_homopolymer`) runs for all protocols
(not just `--dT-primed-cDNA`). Renamed `--polya-sequenced` → `--dT-primed-cDNA`
(deprecated alias preserved). Cat2 validation reads replaced with T-tract
homopolymer examples. Direct RNA / dT-primed cDNA protocol distinction clarified.

**v0.9.0-dev (2026-04-11):** All 4 previously open items resolved:
- Bug 37 (HIGH): `tests/test_terminal_exon_refiner.py` added (51 tests, real wt_by4742_rep1 data)
- Bug 38 (HIGH): `tests/test_consensus_selection.py` added (40 tests, real wt_by4742_rep1 data)
- Bug 41 (MEDIUM): `--polya-model` wired through correction pipeline; `pt_tag`/`polya_score`/`polya_source` columns in TSV; `rectify tag-polya` subcommand; unaligned dorado BAM auto-detection
- Bug 55 (MEDIUM): `--max-cluster-radius`, `--min-peak-sep`, `--min-cluster-samples` added to `analyze` subparser

**v3.1.7 (2026-04-21):** `refine_read_junctions` — bilateral t2, no-candidate-guard policy, canonical HP prior:

- **No candidate guards**: All junctions in the candidate pool are scored; non-canonical, non-annotated (novel) alternatives are no longer filtered out before scoring. The previous guard (`if is_alt==1 and tier>=4 and is_novel==1: continue`) was removed because it would silently discard reads that genuinely belong at non-canonical junctions — e.g. when many reads from the same splice isoform all score perfectly at a novel non-canonical site. Annotation and canonical tier remain as TIE-BREAKERS only, never as gates. This policy is PERMANENT and must not be re-introduced.

- **Bilateral t2 scoring**: `_score_junction` now uses `score(k) = t1(k) + t2(k)` where:
  - `t1(k) = _score_hp_anchored(rescue[k:],       g[je      : je+buf ])` — intron_end-anchored
  - `t2(k) = _score_hp_anchored(rescue[:k][::-1], g[je-buf  : je     ][::-1])` — intron_end-proximal
  Both anchored to `intron_end` (je). t2 prevents degenerate k=L-1 coincidental single-base matches from scoring 0 at the wrong junction without any pre-filtering.

- **Canonical HP prior** (replaces `int(score)` floor binning): when the current N-op is non-canonical (`tier_beats_alt=True`), canonical-tier alternatives (tier < 4) receive a 0.5 edit-distance discount (`_CANONICAL_HP_PRIOR = 0.5`). This equals one Nanopore HP deletion equivalent — the expected noise floor — and ensures canonical junctions win within the noise floor regardless of which penalty table (default or empirical) is in use. Non-canonical alternatives must exceed the canonical score by >0.5 to win. The old `int(score)` floor was fragile: it worked only when both scores happened to fall in the same [n,n+1) integer bucket, which fails when the empirical del_cost(1)≈0.43 and scores straddle an integer boundary.

- **Impact**: all 146 tests pass with both default and empirical penalty table (`penalty_scores.tsv`); all 4 cat9 reads correctly refined with `--junction-penalty-table`.

**v3.2.3 (2026-04-24):** Validation read sequences synced to DRS-trimmed run for Cat5-9:
- **Root cause (v3.1.8 partial update)**: `update_validation_drs.py`'s "update in place" path updated CIGAR/N-op boundaries from the DRS-trimmed mapPacBio run but left read *sequences* from the old chunked-consensus BAM. This caused visual discrepancies in IGV between `rectified_corrected_3end.bam` (corrected from chunked sequences) and `validation_reads.mapPacBio.bam` (aligner BAM with DRS sequences). cat5_minus_1 retained `SEQ='*'` from the gapmm2 PAF issue (fix was never committed).
- **Fix**: targeted 9-read update preserving all XV/XG/XU tags and MD/cs/NM stripped:
  - cat5_minus_1: filled 561-base chimeric sequence from chunked consensus (chimeric CIGAR qlen=561; cannot use DRS mapPacBio seqlen=538 without breaking the chimeric structure)
  - cat6_plus_1, cat6_minus_2: seq+CIGAR+start from DRS mapPacBio merged BAM
  - cat7_minus_1, cat7_plus_1, cat7_plus_2: same
  - cat9_plus_1, cat9_plus_2, cat9_minus_1: from DRS per-chunk mapPacBio BAMs
  - **Cat1–4, Cat8 intentionally kept on chunked sequences**: a full swap to DRS-trimmed caused Cat4 reads to map to wrong genomic regions (minimap2 DRS alignment chose a different locus for cat4_plus_1: 19589 vs expected 22072), and Cat1–2 reads no longer exhibited the indel/soft-clip artifacts required for their validation scenarios.
- **Test updates**: Cat7 EXPECTED_JUNCTIONS shifted 3–8 bp (cat7_plus_1: 138856→138864, cat7_plus_2: 595736→595739); Cat9 RAW_JUNCTIONS shifted 1–3 bp (cat9_plus_1 start: 555824→555825; cat9_plus_2 end: 439324→439321) — same GT-AG annotated junctions, same correction outcomes.
- All 110 tests pass (3 skipped).

**v3.2.2 (2026-04-22):** Validation dataset fully certified for production — XV/XG tag fix + aligner BAM N-op correction:
- **Root cause**: `rebuild_aligner_bams_cat679.py` (v3.2.0) sourced replacement read alignments from `validation_reads.bam` before that BAM's replacement reads had proper N-ops. Additionally, the 3 replacement reads (f8050895, 7d5e8dc2, 72557a9a) were missing XV and XG auxiliary tags because they were inserted with full UUIDs while the source (stale BAM) had 8-char prefix names — the tag-copy logic couldn't match them.
- **XV/XG tag fix**: both tags are now present on all 36 reads in `validation_reads.bam`. XV (read label, e.g. `cat6_plus_2`) is required by `load_reads()` in the test fixture; XG (category, e.g. `cat6_chimeric`) drives category-level test assertions.
- **N-op fix for aligner BAMs**: `dev_runs/wt_by4742_rep1_drs_trim_20260417/fix_replacement_read_cigars.py` re-sourced correct alignments directly from dev run BAMs — `wt_by4742_rep1_drs_trim_20260417` (mapPacBio) and `wt_by4742_rep1_chunked_20260412` (deSALT/gapmm2/minimap2/uLTRA) — preserving full UUID QNAMEs from aligner BAMs and stripping `_pt:i:N` suffixes when matching. Correct N-op intervals: f8050895 all-aligners (45644–45977); 7d5e8dc2 mapPacBio-only (60193–60697); 72557a9a mapPacBio-only (104435–104495).
- **Full pipeline re-run**: `corrected_3ends.tsv` (36 rows: 32 high + 4 chimeric), `rectified_corrected_3end.bam`, `rectified_pA_tail_trimmed.bam`, `rectified_pA_tail_soft_clipped.bam` all regenerated from the fixed validation_reads.bam.
- **Independent certification**: 4 independent agent reviews (sequence-level genomic context, structural checks, category-specific assertions, cross-agent consensus) all confirmed 36/36 PASS. Cat1/Cat2 shifts match VALIDATION_READS.md exactly; Cat3 five_prime_rescued=1; Cat4 n_junctions=1; Cat5 correction=none+chimeric; Cat6/7 XU=1; Cat9 junction refinement corrects all 4 reads.
- All 708 tests pass.

**v3.2.1 (2026-04-22):** `rescue_3ss_truncation` — soft-clip exon CIGAR body-borrowing fix:
- `_get_intronic_query_bases` includes one extra body base when the alignment ends exactly at `intron_start` (boundary straddle: last M op's final ref base == intron_start). For soft-clip rescues this produced a CIGAR with one more query-consuming op than `five_prime_soft_clip_length` (e.g. `14M1D9M` for 79f61403, query_span=23 vs soft_clip=22). The `bam_writer` guard then fell back to a flat `22M` block, applying wrong exon geometry in IGV.
- Fix: for `rescue_type_candidate == 'softclip'`, use `rescue_seq` (exactly the soft-clip, already truncated to `five_clip` bases) as `_align_seq` instead of `_intronic_seq`. For `mpb_mismatch` rescues `_intronic_seq` is still used (it correctly excludes exon-2 body bases beyond the intron boundary).
- Impact: 79f61403 (`cat3_plus_2`, YAL003W) now correctly produces `14M1D7M1D1M` (query_span=22) instead of flat `22M`; exon geometry in IGV matches mapPacBio alignment. 28ea9379 (`cat3_minus_2`, YBR062C) was unaffected (its `_intronic_seq` was already empty in the current validation data).
- All 708 tests pass.

**v3.2.0 (2026-04-22):** Validation aligner BAMs rebuilt to match `validation_reads.bam` after DRS rebuild (v3.1.8):
- The DRS rebuild (v3.1.8) replaced 3 old Cat6/Cat7 reads in `validation_reads.bam` but did NOT update the individual aligner BAMs (`aligners/validation_reads.*.bam`). All 5 aligner BAMs still contained the old reads (ba761413, 64f4da08, 5c59f0bc) and lacked the replacements (f8050895, 7d5e8dc2, 72557a9a).
- `dev_runs/wt_by4742_rep1_drs_trim_20260417/rebuild_aligner_bams_cat679.py` removes the 3 old reads from each aligner BAM and adds the 3 new reads. NOTE: this script sourced alignments from `validation_reads.bam` which at the time also lacked N-ops for the replacement reads — the N-ops were fixed in a subsequent pass by `fix_replacement_read_cigars.py` (see v3.2.2).
- This inconsistency did not break any tests (Cat6 tests use only the basic `corrected` fixture; Cat9 tests use `--aligner-bams` but Cat9 reads were not changed). It was confusing because ba761413 appeared in IGV from aligner BAMs but not from the main rectify list.
- All 708 tests pass.

**v3.1.9 (2026-04-22):** `rectified_pA_tail_soft_clipped.bam` added as default Step 4 output; sort+index added to both Step 4 code paths:
- `run_command.py` Step 4 in both `_process_one_sample` (multi-sample) and `_run_single_sample` (single-sample): after `restore_polya_softclips()` completes, the output BAM is now sorted (via `pysam.sort`) and indexed (via `pysam.index`). Without sort+index, the BAM was written unsorted and had no `.bai` — unusable in IGV/samtools.
- Multi-sample path also gained the missing `threads=getattr(args,'threads',4)` argument to `restore_polya_softclips()`.
- Bundled validation data: `rectified_pA_tail_soft_clipped.bam` generated from `rectified_pA_tail_trimmed.bam` + `wt_by4742_rep1_polya_trim_metadata.parquet`. 32/36 reads have poly-A/adapter soft-clips restored; 3 replacement reads (f8050895, 7d5e8dc2, 72557a9a) have no trim metadata (new UUIDs) and are left unchanged.
- `TestPolyASoftClippedBam` added to `tests/test_validation_reads.py` (4 tests): checks 36 primary reads, coordinate-sorted+indexed, reads-with-metadata have ≥1 soft-clip base, replacement reads present.
- All 708 tests pass.

**v3.1.8 (2026-04-22):** Validation reads updated to DRS pre-trim mapPacBio alignments for Cat6/Cat7/Cat9:
- `validation_reads.bam`: all 12 Cat6/Cat7/Cat9 reads now use DRS-trimmed mapPacBio alignments. Three reads replaced entirely because their pre-DRS mapPacBio alignments lost the N op after DRS trimming: `ba761413→f8050895` (cat6_plus_2, chrII:45644–45977), `64f4da08→7d5e8dc2` (cat6_minus_1, same intron 60193–60697), `5c59f0bc→72557a9a` (cat7_minus_2, chrIII:104435–104495). Nine reads updated in place (same UUID, DRS-trimmed alignment with shifted junction boundaries).
- BAM header reduced to 16 nuclear chromosomes only (chrI–chrXVI); `chrmt`/`chrMito` removed.
- `corrected_3ends.tsv` and rectified BAMs regenerated from DRS-trimmed validation set.
- Test updates: Cat7 `EXPECTED_JUNCTIONS` updated (cat7_plus_1: 138856→138864 start; cat7_plus_2: 595736→595739; cat7_minus_2: 882352-882702→104435-104495 new read); Cat9 `RAW_JUNCTIONS` updated for plus-strand reads (cat9_plus_1: 555824→555825; cat9_plus_2: 439324→439321). Minus-strand raw junctions unchanged (test formula coincidentally gives same values via insertion-length counting artifact).
- Root cause of mapPacBio BAM unmapped-duplicate issue documented: mapPacBio emits unmapped (flag=4) copies of reads alongside the primary alignment; these pass `not is_secondary and not is_supplementary` checks and must be explicitly skipped with `read.is_unmapped` guard.
- All 704 tests pass.

**v3.1.6 (2026-04-21):** DRS pre-trimming integrated as default Step 0+4 in `run-all`:
- `--drs` flag added to `rectify run-all` (cli.py). When set with a BAM input, Step 0 (`trim_drs_bam_polya`) runs before alignment and Step 4 (`restore_polya_softclips`) runs after correction. FASTQ inputs are unaffected.
- `_run_single_sample()` docstring updated to list all 5 steps (0–4).
- `ARCHITECTURE.md` updated: mermaid diagram corrected (Step 3=analysis, Step 4=restore-softclip), `run-all` description expanded with DRS wiring details, `junction_refiner.py` added to directory listing and Layer 4 narrative (Module 2H, ⑦).
- `README.md` accuracy fixes: `rectify run` → `rectify run-all` (command didn't exist); Module 2H scoring description updated (W/max_slide unused, HP-anchored semi-global DP replacing "split-alignment"); tiebreaker priority corrected (score → canonical → annotation → shift, not score → shift → canonical → annotation).
- `junction_refiner.py` `_score_junction` docstring Returns line corrected to match tier-1-only implementation (`min_k t1(k)` not `t1+t2`).

**v3.1.5 (2026-04-21):** `process_bam_streaming_parallel` — two-level checkpoint/resume for `rectify correct`:
- New `--checkpoint-dir DIR` flag (performance group in `cli.py`). Has no effect without `--streaming`.
- **Scan-phase checkpoint**: after `run_variant_aware_scan` (~30–60 min), the `VariantAwareHomopolymerRescue` object is pickled to `<checkpoint_dir>/rescue_scan.pkl`. On resume the pickle is loaded and the scan is skipped entirely.
- **Region-phase checkpointing**: after each region batch completes, a `region_NNNN.done` sentinel file is written to the checkpoint dir. On resume, done regions are filtered from `_regions_to_run` and `_orig_idxs`; the partial output TSV is opened in append mode; the header line is not re-written; `_rebuild_pos_counts_from_partial` re-reads the partial TSV to restore `_pos_counts` for the final position index.
- **Map function**: `pool.imap` (ordered) is used when `checkpoint_dir` is set so that `_batch_num → _orig_idxs[_batch_num]` is stable; `pool.imap_unordered` is used otherwise (unchanged).
- **Failure behaviour**: when `--checkpoint-dir` is set, partial output and sentinels are preserved on failure (warning logged) instead of being deleted. Re-running the identical command resumes from the last completed region.
- Wired through `correct_command.py` via `getattr(args, 'checkpoint_dir', None)`.

**v3.1.4 (2026-04-21):** `refine_read_junctions` — candidate guard + adaptive tie-break scoring:
- **Candidate guard**: non-canonical (tier ≥ 4) AND non-annotated (novel) junctions are now excluded as replacement candidates during scoring. They were previously considered and could "win" by scoring 0.0 (aligner placed a non-canonical junction with a perfect split-alignment score), then get blocked by a post-selection guard that fired `continue` on the N-op, skipping all other candidates including the correct annotated GT-AG junction. The guard is now applied during candidate selection instead: `if is_alt==1 and tier>=4 and is_novel==1: continue`.
- **Adaptive tie-break ordering**: when the current N-op junction is non-canonical (tier ≥ 4), the scoring tuple uses `(score, tier, is_alt, ...)` so a canonical alternative wins at equal edit-distance. When the current junction is acceptably canonical (tier < 4), the tuple uses `(score, is_alt, tier, ...)` so the current junction is preferred at equal score (prevents displacing correctly-placed reads at annotated non-GT-AG junctions like TFC3's RAG site). Controlled by `tier_beats_alt = current_tier >= 4`.
- **Impact**: `cat9_plus_2` (read 00a1e01e, chrVII RPL30, current non-canonical GA-GG tier=6) now correctly refined to annotated GT-AG (439093,439323). TFC3 reads at annotated RAG 3'SS (151006,151096, tier=1) remain stable against the nearby YAG alternative (150989,151096, tier=0). All 704 tests pass.

**v3.1.3 (2026-04-21):** `--aligner-bams` `aligner:path` prefix stripping bug fix:
- `correct_command.py` `_strip_aligner_prefix()`: new private helper that strips the `aligner:` prefix (e.g. `"minimap2:/path/to/file.bam"` → `"/path/to/file.bam"`) before storing BAM paths in `config['aligner_bams']`. Without the fix, pysam received `"minimap2:/path/..."` as a URL scheme → "[Errno 93] Protocol not supported" → 0 novel junctions loaded from aligner BAMs.
- Root cause: `--aligner-bams` accepts the same `aligner:path` format as `rectify consensus`, but `collect_junctions_from_bam` needs a plain file path. The prefix was never stripped before pysam.
- Impact: With the fix, `rectify correct --aligner-bams minimap2:/path gapmm2:/path` loads all aligner-specific novel junctions into the pool (377 vs 362 for the validation dataset). Previously annotation-only junctions were used even when aligner BAMs were specified.
- `_strip_aligner_prefix` checks `'/' not in prefix` to distinguish aligner names from Windows-style drive letters (future-proofing); plain paths without `:` pass through unchanged.

**v3.1.2 (2026-04-21):** Cat9 validation reads added (Module 2H junction refinement):
- `validation_reads.bam` expanded from 32 → 36 reads; aligner BAMs updated to 36 reads each (uLTRA/deSALT stay at 34/36 due to dev-run coverage).
- Four Cat9 reads: `cat9_plus_1` (00a1c9b3, chrVII:555824→555830, +), `cat9_plus_2` (00a1e01e, chrVII:439089→439093, +), `cat9_minus_1` (0b3b593b, chrXV:900760→900767, −), `cat9_minus_2` (d3357db5, chrXV:900760→900767, −). Each has `XG=cat9_junction_refine`.
- `TestCategory9JunctionRefinement` class in `tests/test_validation_reads.py`: `corrected_with_aligner_bams` fixture runs `rectify correct --aligner-bams --annotation`; verifies corrected junction ≠ wrong junction for all 4 reads.
- `TestBamIntegrity` updated: 32→36 reads, cat9 labels added to `test_all_labels_present`, `cat9→cat9_junction_refine` in `test_category_tags`, cat9 added to `test_strand_balance`.
- Bundled validation outputs regenerated: `corrected_3ends.tsv`, rectified BAMs, PROVENANCE.json all updated to reflect 36 reads.

**v3.1.1 (2026-04-21):** `_score_junction` + `refine_read_junctions` — two correctness fixes:
- **`range(L)` fix**: the k-loop in `_score_junction` previously ran `range(L+1)`. At k=L, `q1 = rescue[L:]` is empty → `_score_hp_anchored` returns 0.0 for every candidate regardless of junction quality. All candidates tied at 0.0, making selection arbitrary. Changed to `range(L)`: q1 always has ≥1 base; scoring is always discriminating.
- **`is_alt` tiebreaker in `refine_read_junctions`**: added `is_alt` (0 if candidate = current N-op, 1 otherwise) as the second tuple element after `score`. When multiple candidates tie (e.g. two junctions share the same `intron_end` and both score 0.0), the existing junction is preferred over alternatives, preventing spurious displacement of already-correct reads.
- **Tests**: 9 previously failing tests now pass; all 41 tests in `test_junction_refiner.py` pass. Total suite: 698 tests pass, 4 skipped. Key reads: RPL20B `0b3b593b` corrected to `[900767,901193)`; TFC3 annotated junction stable against alternative `[150989,151096)`.

**v3.1.0 (2026-04-20):** Module 2H — post-consensus N-op junction refinement (`junction_refiner.py`):
- New module `rectify/core/junction_refiner.py`: for every N-op in every consensus read, tests all candidate junctions within a search radius and replaces imprecise N-op boundaries with the best sequence-supported junction.
- Scoring is **sequence-first**: hp_score (split-alignment with homopolymer-aware edit distance) is the primary criterion; canonical GT-AG and annotated status are tie-breakers only. Annotation NEVER overrides a better-scoring junction.
- Split-alignment: query window of ±W bp around the current N-op split point; query split can slide ±max_slide bp; scores both exon2 and exon1 genomic context simultaneously.
- Fast path: reads already at an annotated canonical-tier-0 junction skip scoring entirely (255× speedup over naive approach).
- `max_boundary_shift=50`: prevents false-positive matches from junctions in neighbouring genes. The `search_radius=5000` discovers candidates; `max_boundary_shift` constrains individual endpoint shifts.
- CIGAR surgery encodes boundary changes as I/D ops (not M adjustments) to preserve both ref and query spans. MD/cs tags stripped from modified reads to avoid pysam arithmetic errors on CentOS 7.
- Wired into `rectify correct` via `--aligner-bams`; also exposes `--junction-hp-pen`, `--junction-search-radius`, `--junction-window`, `--junction-max-slide`, `--junction-max-boundary-shift`.
- Key fixed reads: b2c4d195 (27S → 426N+exon1), f6590560 (422N1M1D1M → 426N2M), 0b3b593b (431N wrong boundaries → 430N correct exon1 end).
- 41 tests in `tests/test_junction_refiner.py` covering RPL20B, GCR1, TFC3, RPL22B, SRC1. (9 were failing at initial merge; fixed in v3.1.1.)

**v3.0.4 (2026-04-20):** `rescue_3ss_truncation` — minus-strand soft-clip truncation bug fix:
- `_extract_5prime_rescue_seq` returns `query_seq[n - last_imp_q:]` for minus-strand reads (the rightmost portion). The truncation to `five_clip` bases was using `[:five_clip]` (leftmost), taking aligned exon2 body instead of the actual soft-clip bases.
- Fix: `rescue_seq[-five_clip:]` for minus strand, `rescue_seq[:five_clip]` for plus strand (unchanged).
- Impact: b2c4d195 (27S), 7bf94550 (1S), f6590560 (2S), 08f6ddf7 (5S) all now correctly rescued to `five_prime_corrected=901193`. Previously all fell through to `rescue_type='proximity'` with `rescued=False`.

**No open bugs remaining** — see `docs/BUGS_TO_FIX.md`.

**v3.0.3 (2026-04-16):** `find_polya_boundary` — poly-A tail trailing-base false-stop guard:
- `find_polya_boundary` in `indel_corrector.py` now detects false stops where the trailing base of a poly-A tail (e.g., a T at the very end of `...AAAAAAAAAAAAAAAAT`) coincidentally matches a genomic base (T=T) at the alignment boundary, causing the backward scan to stop prematurely at the poly-A/exon junction rather than continuing to the true exon body.
- Fix (plus strand): Before accepting a candidate stop (`rb==gb, gb!='A'`), inspect the K=4 positions to the left. If all K have `rb='A'` AND at least one has `gb≠'A'` (unmistakably poly-A tail context), the candidate is skipped and scanning continues leftward.
- Fix (minus strand): Symmetric check — inspect K=4 positions to the right; if all have `rb='T'` AND at least one has `gb≠'T'`, continue scanning rightward.
- Root cause in cat6_plus_2 (RPL19B, ba761413): alignment ends `...AAAAAAAAAAAAAAAAT` at 169492 (T=T exact match). The 14 A-positions before it (169478–169491) all have `rb='A'` with several mismatches (`gb≠'A'`) — pure poly-A tail. With the guard, the scan continues to 169476 (G=G, true exon end). `corrected_3prime` moves 16 bp upstream from 169492 to 169476; NET-seq refinement stays at 169476 (signal=75.0). Bedgraph now shows signal at 169476 rather than 169491.
- All 659 tests pass.

**v3.0.2 (2026-04-16):** `clip_read_to_corrected_3prime` / `softclip_read_to_corrected_3prime` — terminal D/N stripping:
- Both functions now strip trailing D/N ops (plus strand) or leading D/N ops (minus strand) that are left dangling after the CIGAR walk loop. The bug arose when `corrected_3prime` fell within a deletion span: the loop removed all query-consuming ops to the right of the D (satisfying `n_ref_removed >= n_ref_clip`) and exited before reaching the D, then appended `H`/`S` directly after it, producing invalid CIGAR strings like `4D6H`.
- Example: read 299e1402 (chrII plus strand) CIGAR previously ended `5X3=4D6H` — a terminal deletion before a hard-clip that violated SAM spec. Now the 4D is stripped before the H is appended.
- Root cause: `corrected_3prime` from `find_polya_boundary` can land inside a deletion span (a reference position with no corresponding query base), which is a valid CPA but requires CIGAR surgery to snap to the last real query base.
- Fix applied in four locations: `clip_read_to_corrected_3prime` (plus + minus strand paths) and `softclip_read_to_corrected_3prime` (plus + minus strand paths).
- New test file: `tests/test_bam_writer.py` (11 tests covering normal clipping, terminal-D stripping, in-deletion-span corrected positions, and N-op stripping for spliced reads).

**v3.0.1 (2026-04-15):** `clip_intronic_tail_5prime` — off-by-one fix + trailing-I/S stripping + existing-H handling + `_MIN_SC_FOR_JUNCTION_EXTENSION` guard:
- **Off-by-one fix (minus strand)**: Exit condition changed from `<= clip_boundary + 1` → `<= clip_boundary`. Previously reads ending at `reference_end = intron_start + 1` (last mapped base = `intron_start` = first intron base) were left unclipped, appearing as red T→C mismatches in IGV. Now they are trimmed so the last mapped base is `intron_start - 1` (last exon base).
- **Trailing I/S stripping**: Before the main ref-consuming trim loop, `clip_intronic_tail_5prime` now explicitly strips trailing I (insertion) and S (soft-clip) ops from the 5' end. These bases lie at/past the intron boundary and were silently left in the CIGAR under the old code when `reference_end` was already at the boundary.
- **Existing H handling**: Any pre-existing trailing H op is extracted before processing and merged back at the end. This prevents double-counting query bases that were already removed from `query_sequence` in a prior call.
- **`_MIN_SC_FOR_JUNCTION_EXTENSION = 3`**: New module-level constant in `bam_writer.py`. All three write functions (`write_corrected_bam`, `write_softclipped_bam`, `write_dual_bam`) now only call `extend_read_5prime_for_junction_rescue` when the 5' soft clip is ≥ 3 bp. Reads with 1-2 bp soft clips fall through to `clip_intronic_tail_5prime` instead, preventing spurious Cat3 exon-N-exon CIGAR surgery on single-base alignment artefacts. An `_extended` flag prevents `clip_intronic_tail_5prime` from undoing a successful Cat3 extension.
- **`bam_processor.py` `_sc_at_5p == 0` restriction removed**: `five_prime_intron_clip_pos` is now set for ALL `_in_intron` reads regardless of whether a soft clip is present at the 5' end. Reads with small soft clips now get the clip applied rather than being directed to Cat3 extension.
- Impact: zero reads with last base at intron_start (X mismatch), zero trailing I ops at boundary, 292 reads properly clipped to `ref_end = intron_start` (last base = last exon base); all 2,073 reads pass CIGAR/sequence-length validation.

**v3.0.0 (2026-04-15):** `clip_intronic_tail_5prime` — BAM sequence trimming fix + generalised intron-clip trigger:
- `clip_intronic_tail_5prime` in `bam_writer.py` now trims `query_sequence` and `query_qualities` to match the new CIGAR after adding H ops. Previously only the CIGAR was updated, leaving `query_sequence` too long → pysam wrote malformed BAM records → `samtools sort: truncated file` error and silent loss of all intronic-snap clip corrections.
  - Minus strand (clip from right): `read.query_sequence = seq[:-clipped_query_bases]`
  - Plus strand (clip from left): `read.query_sequence = seq[clipped_query_bases:]`
- `bam_processor.py` `five_prime_intron_clip_pos` assignment generalised: previously only set for `rescue_type == 'intronic_snap'` (Case 4). Now set for ALL rescues (Cases 1/2/4) where the alignment's 5' end falls inside the rescued intron AND there is no 5' soft-clip. Case 1/2 reads rescued with no soft clip (e.g. d36f9748: `16M` exon match, alignment ending 16 bp inside intron) also need the intronic tail hard-clipped.
- Impact: clipped BAM count rose from 1,287 → 1,355; all 2,073 reads pass CIGAR/sequence-length validation.

**v2.9.1 (2026-04-12):** Cat2 soft-clip rescue CIGAR surgery:
- `rescue_softclip_at_homopolymer` now stops at `A` (plus strand) or `T` (minus strand) in the soft-clip, preventing poly-A tail bases from matching genomic A-runs. Fixes shifted corrected positions for cat2_plus_2 (+10→+9), cat2_minus_1 (-17→-10), cat2_minus_2 (-12→-11).
- New `extend_read_3prime_for_softclip_rescue()` in `bam_writer.py`: converts the 3' soft-clip to `{D}D{M}M{poly-A}H|S`, making the true RNA 3' end visible in IGV.
- Cat2 rescue metadata (`sc_homopolymer_extension`, `sc_rescued_seq`, `sc_original_softclip_len`) now written to corrected_3ends.tsv and read back by bam_writer.
- Bundled BAMs renamed: `rectified.bam` → `rectified_corrected_3end.bam`; `rectified_softclip.bam` → `rectified_pA_tail_trimmed.bam`.

**v2.9.7 (2026-04-14):** `_iter_name_grouped_bams` — natural sort K-way merge fix:
- The K-way merge in `consensus.py` used Python `min()` (lexicographic) to pick the next read ID, but the name-sorted BAMs use `SS:queryname:natural` (samtools natural sort, digit runs compared as integers). For UUID-format read names, these orderings diverge: `98297e97` (key=98297) sorts BEFORE `0633141e` (key=633141) in natural sort, but AFTER it lexicographically ('9' > '0').
- When a read present only in aligner A sits between two reads in natural sort, the lexicographic merge incorrectly pulled the shared read from aligner B first (as a single-aligner group), then later processed aligner A's copy as a separate single-aligner group. Both single-aligner groups "won" by default → duplicate records with different aligners/strands in the output BAM.
- Fix: added `_natural_sort_key()` (splits on `(\d+)`, converts digit tokens to `int`) and passed it as `key=` to the `min()` call in `_iter_name_grouped_bams`. The merge now uses the same natural sort as the BAMs.
- Impact on RPL19B+RPL20B test: mixed-strand duplicate QNAMEs in consensus BAM dropped from 51 → 1 (the 1 remaining is an inherent FASTQ-duplicate edge case where copies map to opposite strands). Plus-strand soft-clipped reads on chrXV dropped from 53 → 9 (residual are other-region reads + FASTQ-dup pairs).

**v2.9.9 (2026-04-15):** `rescue_3ss_truncation` — homopolymer-aware edit distance + 3'SS acceptor tiebreaker:
- New `_hp_edit_distance(s1, s2) -> float`: indels within a homopolymer run (indel base equals its immediate neighbour in the same string) cost 0.5 instead of 1.0; substitutions always cost 1.0. Used for all four edit-distance calls in `rescue_3ss_truncation` (inner-loop candidate scoring, outer `ed_exon`, and intron comparison `ed_intron`). Replaces `_edit_distance` (integer Levenshtein) for all splice-rescue scoring; `_edit_distance` is retained for other callers.
- New `_ACCEPTOR_PRIORITY_PLUS` and `_ACCEPTOR_PRIORITY_MINUS` module-level dicts encode 3'SS acceptor dinucleotide priority: AG=0, CG=1, TG=2, AT=3, other=4 (plus strand, last 2 bases of intron); minus strand uses RC equivalents (CT=0, CG=1, CA=2, AT=3).
- `_acceptor_priority` is computed once per junction (outside the shift loop) and added as the 5th element of the outer tiebreaking tuple: `(ed_exon, not_canonical_donor, not_in_amb, shift_abs, acceptor_priority)`.
- Outer tracking now also records `best_in_amb`, `best_shift_abs`, `best_acceptor_priority` alongside the existing `best_ed` and `best_is_canonical`.

**v2.9.8 (2026-04-15):** `rescue_3ss_truncation` — two-phase discovery + canonical refinement, ±5 bp baseline for both strands:
- Plus and minus strand branches now both use `max(5, amb)` as the shift range floor (previously `max(1, amb)`). The ±5 bp baseline catches imprecise aligner annotations (e.g. off by 3-4 bp) even when local sequence ambiguity is only 1 bp.
- Minus strand canonical check expanded: `in ('AC', 'GC')` (RC of GT/GC) instead of `== 'AC'` only. Mirrors plus strand which uses `in ('GT', 'GC')`.
- Minus strand now uses the same tuple scoring as plus strand: `(not_canonical, not_in_amb, shift_abs)` — canonical wins first, then within-ambiguity-window shift, then smallest |shift|. Previously used a simpler chain lacking `_in_amb` refinement.
- Three production-data failed rescues motivating this fix (not yet in validation set): chrI:142,611-142,653 (EFB1/YAL003W, plus strand), chrII:366,454-366,654 (YBR062C, minus strand), chrXI:20,510-20,584 (plus strand, poly-A run).

**v2.9.6 (2026-04-14):** `rescue_3ss_truncation` — splice-site ambiguity resolution via data-driven shift range and canonical GT-AG preference:
- For each candidate junction, the exon-end scan tries a range of donor/acceptor position shifts derived from the actual run-length of matching bases at the boundary: `r_amb` = consecutive intron bases (going rightward from the annotated donor) that equal the last exon base; `l_amb` = consecutive exon bases (going leftward) that equal the first intron base. The shift range is `[-max(5, l_amb), max(5, r_amb)]`, capped at ±15. The minimum ±5 ensures off-by-one errors (e.g. 1-based GFF vs 0-based) and imprecise aligner annotations are always caught.
- Scoring priority: lower edit distance wins first; among ties, canonical GT/GC donor (plus strand) or AC/GC (minus strand genomic) wins over non-canonical; among remaining ties, within-ambiguity-window shifts win over outside; then smallest |shift| wins.
- The data-driven range is critical for correctness: a fixed wide range (e.g. ±15) would spuriously try shifts deep into the intron, finding false perfect matches when exon1-end + intron-start looks identical to intron-end + exon2-start (e.g. `TTTCAG-GTA` at both the 5'SS and 3'SS). The run-length bound prevents this.
- For long homopolymer runs at the junction, the range automatically extends to cover the full ambiguous window (up to 15 bp). For typical 1-bp repeats (like G|G at YAL003W), the range is ±5 (dominated by the baseline).
- Example: YAL003W (chrI) read 79f61403. ref[142252]=G (last exon base), ref[142253]=G (first intron base, canonical GT donor at 142253-142254). With the wrong junction (intron_start=142254), shift=-1 gives intron_start=142253 (GT canonical) with the same edit distance → canonical wins → five_prime_corrected=142252, exon_cigar=14M1D7M1D1M. Deterministic across all set iteration orders.

**v2.9.5 (2026-04-14):** `rescue_3ss_truncation` — Case 4 intronic-snap rescue:
- New `_get_n_op_intervals()` helper extracts (start, end) genomic intervals for every N-op in a CIGAR.
- New Case 4 in `rescue_3ss_truncation()`: fires when `align_5prime` is strictly inside `[intron_start, intron_end)` AND no existing N-op already covers that intron (checked with ±`junction_proximity_bp` tolerance). Cases 1–2 already handle this via terminal-error sequence matching; Case 4 handles the remainder (clean alignments with D-ops or plain truncation inside the intron, no detectable terminal errors).
- Action: snaps `five_prime_corrected` to `intron_end` (minus strand) or `intron_start` (plus strand) — the exon-1-side splice donor boundary. Records `rescue_type = 'intronic_snap'`, `rescued = True`.
- Works for all aligner representations: does not depend on CIGAR op type, only on final `align_5prime` position.
- Impact on RPL20B test set: reduced reads with `five_prime_position` inside annotated intron from 170 → 0; added 281 `intronic_snap` rescues at `intron_end = 901193`.

**v2.9.4 (2026-04-14):** `find_polya_boundary` — N-op boundary guard for spliced minus-strand reads:
- `find_polya_boundary` in `indel_corrector.py` now records the first N-op (intron skip) boundary encountered while parsing the CIGAR for minus-strand reads and limits the forward scan to positions before that boundary (`scan_limit`). Previously the scan silently skipped N ops (they are absent from `aligned_positions`) and could find a spurious non-T match in a downstream exon, producing a huge erroneous correction that spanned the intron.
- If no non-T match is found before the N boundary and the last pre-N position is T (poly-A zone), `first_n_start` is returned as the CPA — the intron boundary is the natural exon-end for the terminal exon.
- Plus-strand backward scan is unchanged: for plus-strand reads the 3' terminal exon is to the right and crossing N ops leftward into earlier exons is legitimate (the entire last exon may be in the poly-A zone, requiring the scan to continue into exon N-1 to find the true non-A boundary).
- Example: YIL145C (chrIX) minus-strand read a9706bbe with CIGAR `1S11M223N1D34M...`. All 11M positions (76016-76026) are T's; scan previously crossed the 223N and found seq[12]='G'==genome[76251]='G', giving a 235bp correction. With the fix, scan stops at scan_limit (exon 1 positions only), finds no non-T match, and returns first_n_start=76027 as the CPA (11bp correction). The 12-row ambiguity problem is resolved.

**v2.9.3 (2026-04-14):** `find_polya_boundary` — large-deletion pre-scan for poly-A over-calling artifacts:
- `find_polya_boundary` in `indel_corrector.py` now detects large deletions (≥5bp) within 50bp of the RNA 3' end that are poly-A over-calling artifacts. The artifact: minimap2 over-extends the poly-A tail alignment into a mismatching region, then uses a large deletion to bridge back to the next exon match. The false-positive non-T (minus strand) or non-A (plus strand) matches before the deletion caused the walkback to stop too early (under-correcting by tens of bp).
- Pre-scan is gated on whether the alignment starts in a poly-T/poly-A context (first real aligned position has `gb='T'` for minus strand, `gb='A'` for plus strand), preventing false activation for reads whose 3' end is NOT in a poly-A zone (e.g., exon body deletions).
- Example: RPL20B (YOR312C) minus strand reads with CIGAR `3=1X1=2X3=18D...` at chrXV:900,120 — previous correction was 2bp (stopped at coincidental G match before the 18D), now correctly 28bp (stops at first exon match after the 18D, within 2bp of the annotated CPA at 900,150).

**v2.9.2 (2026-04-14):** Cat3 5' rescue — mapPacBio intrusion fix:
- `rescue_3ss_truncation()` in `splice_aware_5prime.py` now handles the case where mapPacBio's alignment extends INTO the upstream intron (`align_5prime` lands between `intron_start` and `intron_end`). Previously the `dist < 0` guard skipped these junctions entirely before any sequence comparison; now they are treated as `dist = 0` (touching the boundary) and proceed to exon-sequence alignment. Applies symmetrically for both strands.

**v0.9.0-dev (2026-04-12):** Parallel alignment infrastructure:
- `rectify split` — splits FASTQ.GZ into N equal chunks (round-robin interleave) for SLURM array jobs; `--generate-slurm` emits `run_array_align.sh` (array 0-79) and `run_merge_and_consensus.sh`
- `rectify consensus` — aligner selection on pre-built BAMs; accepts `aligner:path` pairs; used after merging per-chunk per-aligner BAMs from array jobs
- `rectify install-aligners` — downloads/installs external aligners; `--check` shows status; `--all` installs everything; deSALT gets vendored binary (Linux/x86_64) or source build fallback
- deSALT vendored binary bundled at `rectify/data/bin/linux_x86_64/deSALT` (v1.5.6, 773 KB); `_get_vendored_desalt()` in `multi_aligner.py` resolves automatically when deSALT not on PATH
