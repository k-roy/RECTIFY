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

---

## HPC I/O: always stage through $SCRATCH

**Never run correction I/O directly on Oak in a SLURM array job.**

Oak is NFS-backed shared storage. Running 20 concurrent array tasks each
streaming a 7 GB BAM through Oak causes severe I/O contention and inflates
wall time 2–3×. On Sherlock, `$SCRATCH` provides ~75 GB/s aggregate
bandwidth and is node-local to each task.

### Pattern (Python pipeline)

`run_command._run_single_sample()` handles this automatically:

```python
from rectify.slurm import make_job_scratch_dir, sync_to_oak

scratch = make_job_scratch_dir('rectify_single')  # None if no scratch
work_dir = scratch if scratch else output_dir

# ... all alignment + correction + analysis writes to work_dir ...

if scratch:
    sync_to_oak(scratch, output_dir)   # rsync everything back
    shutil.rmtree(scratch)
```

### Pattern (generated SLURM scripts)

`batch_command.generate_slurm_scripts()` inserts staging blocks when
`use_scratch=True` in the profile (default in `sherlock_larsms.yaml`):

- BAM inputs are `cp`'d to `$SCRATCH` before `rectify correct`
- FASTQ inputs are left on Oak (read sequentially once, no benefit to staging)
- All outputs are `rsync`'d back to Oak after correction

### $SCRATCH caveats

- Auto-purged after **90 days** on Sherlock — final outputs must land on Oak
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

**Always use `--streaming` in SLURM jobs.** The `sherlock_larsms.yaml`
profile sets `streaming: true` by default.

### Known bug (fixed 2026-03-30)

`process_bam_streaming()` called `correct_read_3prime()` which returns
`List[Dict]`, then passed that list directly to `stats.update_from_result()`
which expected a single `Dict`. Fix: use `.extend()` and iterate for stats.

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

## File layout

```
rectify/
├── rectify/
│   ├── slurm.py              # CPU detection, thread limits, scratch utilities
│   ├── slurm_profiles/
│   │   ├── sherlock_larsms.yaml  # Standard CPU partition (use_scratch, streaming on)
│   │   └── sherlock_gpu.yaml     # GPU partition
│   ├── core/
│   │   ├── run_command.py    # Active run-all dispatcher (single + multi-sample)
│   │   ├── batch_command.py  # Generates SLURM array scripts; reads profile YAMLs
│   │   ├── correct_command.py
│   │   ├── bam_processor.py  # process_bam_file_parallel + process_bam_streaming
│   │   └── multi_aligner.py  # minimap2/mapPacBio/gapmm2 wrappers
│   └── ...
└── CLAUDE.md                 # This file
```
