# rectify split

Split a FASTQ/FASTQ.GZ file into N equal chunks for parallel SLURM array alignment.

Reads are assigned to chunks round-robin (interleaved), so each chunk receives an even
distribution of read lengths even when reads are coordinate-sorted.

---

## Usage

```bash
rectify split <reads> -n <N> -o <output_dir> [options]
```

## Examples

```bash
# Split into 16 chunks (dry run — counts reads without writing)
rectify split reads.fastq.gz -n 16 -o chunks/ --dry-run

# Split and write chunks
rectify split reads.fastq.gz -n 16 -o chunks/

# Split and generate SLURM array scripts for 5 aligners
rectify split reads.fastq.gz \
    -n 16 \
    -o /scratch/chunks/ \
    --generate-slurm \
    --aligners minimap2 mapPacBio gapmm2 uLTRA deSALT \
    --genome /ref/genome.fa.gz \
    --annotation /ref/genes.gff.gz \
    --slurm-partition my-partition \
    --slurm-cpus 16 \
    --slurm-mem 64G \
    --slurm-time 12:00:00
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `reads` | — | Input FASTQ or FASTQ.GZ file |
| `-n, --n-chunks` | 16 | Number of chunks to create |
| `-o, --output-dir` | — | Output directory |
| `--prefix` | derived | Output file prefix (default: input filename stem) |
| `--dry-run` | off | Count reads and print chunk sizes; do not write files |
| `--generate-slurm` | off | Write SLURM array script and merge script alongside chunks |
| `--verbose` | off | Verbose logging |

### SLURM script options (with `--generate-slurm`)

| Argument | Default | Description |
|----------|---------|-------------|
| `--genome` | — | Reference genome path (written into generated scripts) |
| `--annotation` | — | Annotation GFF/GTF path |
| `--aligners` | minimap2 mapPacBio gapmm2 | Aligners to include in array job |
| `--slurm-partition` | — | SLURM partition(s) |
| `--slurm-account` | — | SLURM account |
| `--slurm-cpus` | 16 | CPUs per array task |
| `--slurm-mem` | 64G | Memory per array task |
| `--slurm-time` | 12:00:00 | Time limit per array task |
| `--python-path` | conda base | Explicit path to Python interpreter |

---

## Output files

After running with `-n 16`:

```
output_dir/
├── sample_chunk_000_of_016.fastq.gz   # chunk 0 (~1/16 of reads)
├── sample_chunk_001_of_016.fastq.gz
│   ...
├── sample_chunk_015_of_016.fastq.gz   # chunk 15
└── chunks_manifest.json               # JSON listing all chunk paths
```

With `--generate-slurm`:

```
output_dir/
├── ...chunk files...
├── chunks_manifest.json
├── run_array_align.sh          # submit this first
├── run_merge_and_consensus.sh  # run after array completes
└── slurm_logs/                 # log directory (created at submit time)
```

---

## Full chunked-alignment workflow

### Step 1 — Split

```bash
rectify split reads.fastq.gz \
    -n 16 -o /scratch/chunks/ \
    --generate-slurm \
    --aligners minimap2 mapPacBio gapmm2 uLTRA deSALT \
    --genome genome.fa.gz --annotation genes.gff.gz
```

### Step 2 — Submit array (80 tasks: 16 chunks × 5 aligners)

```bash
sbatch /scratch/chunks/run_array_align.sh
```

Each task ID decodes as:

```bash
CHUNK_IDX=$(( SLURM_ARRAY_TASK_ID % N_CHUNKS ))   # 0–15
ALIGNER_IDX=$(( SLURM_ARRAY_TASK_ID / N_CHUNKS ))  # 0–4
```

Tasks run `rectify align --no-consensus`, writing one sorted BAM per (chunk, aligner) pair.

### Step 3 — Merge and consensus (after array completes)

```bash
bash /scratch/chunks/run_merge_and_consensus.sh
```

This script:
1. `samtools merge` per aligner — combines all chunk BAMs into one sorted BAM per aligner
2. `rectify consensus` — selects the best aligner per read and writes the final rectified BAM

---

## Scheduler compatibility

The generated `run_array_align.sh` uses `SLURM_ARRAY_TASK_ID` and `SLURM_CPUS_PER_TASK`.

For non-SLURM schedulers, add a shim at the top of the array script body:

| Scheduler | Task ID | CPUs |
|-----------|---------|------|
| SLURM | `SLURM_ARRAY_TASK_ID` (0-based) | `SLURM_CPUS_PER_TASK` |
| UGE/SGE | `SGE_TASK_ID` (1-based) | `NSLOTS` |
| PBS/Torque | `PBS_ARRAY_INDEX` | `PBS_NUM_PPN` |

**UGE/SGE shim** (add before the `CHUNK_IDX=` line):

```bash
# UGE compatibility shim
SLURM_ARRAY_TASK_ID=$(( SGE_TASK_ID - 1 ))   # SGE is 1-based
SLURM_CPUS_PER_TASK=${NSLOTS:-8}
```

**PBS/Torque shim**:

```bash
# PBS compatibility shim
SLURM_ARRAY_TASK_ID=$(( PBS_ARRAY_INDEX - 1 ))   # PBS is 1-based
SLURM_CPUS_PER_TASK=${PBS_NUM_PPN:-8}
```

---

## Notes

- Chunks use round-robin assignment so read-length distributions are equal across chunks even when the input is coordinate-sorted. This prevents some chunks being fast (short reads) and others slow (long reads).
- The `--dry-run` flag reads the entire file to count reads but writes nothing. Use it to preview chunk sizes before committing to a long split.
- Chunk files are gzip-compressed regardless of whether the input was gzipped.
