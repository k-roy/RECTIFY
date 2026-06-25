# Multi-sample analysis: memory-efficient streaming pipeline

For DESeq2 to run, you must use `--manifest` with a TSV containing
`sample_id`, `path`, and `condition` columns.

**Do NOT pre-combine corrected TSVs into a single file.** The old
approach (concat all per-sample TSVs → one file → `rectify analyze`)
requires loading the entire dataset into RAM and OOMs for large
experiments. Use manifest mode instead.

## The correct pipeline (streaming, low RAM)

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

## How manifest mode works (two-pass streaming)

**Pass 1 — position aggregation for clustering:**
Each per-sample TSV (or its index, see below) is read sequentially.
Only `chrom`, `strand`, `corrected_3prime` are loaded; all other columns
are dropped. Positions are aggregated to unique (chrom, strand, pos)
counts per sample and combined for clustering. Never loads >1 sample at
a time.

**Pass 2 — count matrix accumulation:**
For each sample, positions are streamed in 100k-row chunks and looked up
against the cluster IntervalTree. Counts accumulate in a
`defaultdict[cluster_id][sample_id]` — a ~1 MB dict regardless of dataset
size. The final count matrix is built from this dict.

**DESeq2, PCA, GO, motifs** then run on the count matrix
(~10k × n_samples), which fits in a few MB.

**Note:** Bedgraph and genomic distribution steps are skipped in
manifest mode (they require per-read alignment coordinates). Generate
bedgraphs separately with `rectify export` if needed.

Peak RAM: O(clusters × samples) ≈ a few MB, regardless of read depth or
sample count. Validated on 21 samples / 150M reads: runs on a 16 GB node.

## Position index (Tier 3): even faster

`rectify correct` now writes `corrected_reads_index.bed.gz` alongside
`corrected_reads.tsv`. Pre-aggregated position count file:

```
chrom   corrected_3prime  strand  count   count_ag_rich
chrI    12836             +       1.0     0.0
chrI    41981             +       3.0     1.0
...
```

(`count_ag_rich` is the AG-rich subset of `count`; older indices without the
column load with `count_ag_rich = 0`.)

~300× smaller than the full per-read TSV. When an index is present,
both passes of manifest mode use it instead of the full TSV — Pass 1 and
Pass 2 become near-instant.

**Generating indices for existing corrected TSVs:**

```python
from rectify.core.bam.bam_processor import write_position_index
from concurrent.futures import ThreadPoolExecutor

samples = ['wt_rep1', 'wt_rep2', ...]
base = Path('results/')

def gen(s):
    tsv = str(base / s / 'corrected_reads.tsv')
    write_position_index(tsv, tsv)

with ThreadPoolExecutor(max_workers=4) as ex:
    list(ex.map(gen, samples))
```

## Column pruning (Tier 1)

`load_corrected_positions()` drops all columns not needed downstream
immediately after loading. Retained: `chrom`, `strand`,
`corrected_position`, `sample`, `fraction` (if present),
`alignment_start`, `alignment_end`. Everything else (`read_id`,
`confidence`, `polya_length`, QC flags, etc.) is dropped on load.

Chunked-loading threshold lowered from 5 GB → 500 MB; aggregation now
kicks in for most real datasets.

## Reference condition matching

`--reference` is matched **case-insensitively** against the manifest's
`condition` column. `--reference wt` will match `WT`, `Wt`, etc. If no
match is found, a warning is printed and the value is used as-is.
