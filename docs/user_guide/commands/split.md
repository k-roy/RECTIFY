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
# Auto-size chunks (default 50000 reads/chunk; dry run — counts reads without writing)
rectify split reads.fastq.gz -o chunks/ --dry-run

# Fixed 16 chunks
rectify split reads.fastq.gz -n 16 -o chunks/

# Auto-size + generate SLURM array scripts for the non-mapPacBio aligners
rectify split reads.fastq.gz \
    -o /scratch/chunks/ \
    --generate-slurm \
    --other-aligners minimap2 gapmm2 uLTRA deSALT \
    --genome /ref/genome.fa.gz \
    --annotation /ref/genes.gff.gz \
    --slurm-partition my-partition

# Bundled yeast reference (auto-fills --genome and --annotation)
rectify split reads.fastq.gz -o chunks/ --generate-slurm --Scer
```

---

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `reads` | — | Input FASTQ or FASTQ.GZ file |
| `-n, --n-chunks` | auto-sized from `--target-reads-per-chunk` (default 50000) | Number of chunks (overrides auto-sizing). Mutually exclusive with `--target-reads-per-chunk`. |
| `--target-reads-per-chunk` | 50000 | Auto-size chunks to ~this many reads each (min 4, max 500) |
| `-o, --output-dir` | — | Output directory |
| `--prefix` | derived | Output file prefix (default: input filename stem) |
| `--dry-run` | off | Count reads and print chunk sizes; do not write files |
| `--generate-slurm` | off | Generate scheduler array scripts alongside chunks |
| `--verbose` | off | Verbose logging |

### Script generation options (with `--generate-slurm`)

| Argument | Default | Description |
|----------|---------|-------------|
| `--genome` | — | Reference genome path (written into generated scripts) |
| `--annotation` | — | Annotation GFF/GTF path |
| `--Scer` / `--organism` | — | Use bundled *S. cerevisiae* genome + annotation (autofills `--genome`/`--annotation`) |
| `--other-aligners` | minimap2 gapmm2 uLTRA deSALT | Non-mapPacBio aligners to include in the "others" array job (mapPacBio runs in its own array) |
| `--skip-map-pacbio` | off | Omit the mapPacBio array script (e.g. not installed) |
| `--scheduler` | slurm | Target scheduler for generated headers (`slurm`, `uge`, `pbs`) |
| `--slurm-partition` | — | SLURM partition(s) |
| `--slurm-account` | — | SLURM account |
| `--uge-queue` | `long.q` | UGE/SGE queue name (`-q`) |
| `--uge-pe` | `smp` | UGE/SGE parallel environment (`-pe`) |
| `--pbs-queue` | `workq` | PBS queue name (`-q`) |
| `--python-path` | `python` | Explicit path to Python interpreter |
| `--rectify-src` | `.` | Path to RECTIFY source checkout (working directory for generated scripts) |
| `--junction-penalty-table` | — | Empirical HP junction penalty table forwarded to correction scripts |
| `--str-penalty-table` | — | STR slippage penalty table forwarded to correction scripts |
| `--junction-overhang-table` | — | Empirical junction overhang table forwarded to chunk-merge scripts |
| `--short-read` | off | Generate the short-read (bbmap + bwa) pipeline instead of the long-read panel; implies `--skip-map-pacbio` |
| `--dT-primed-cDNA` | off | Pass `--dT-primed-cDNA` through to generated `rectify correct` calls (auto-enabled by `--short-read`) |
| `--oak-output-dir` | `{output-dir}/final/` | Destination on persistent storage for the final merged consensus BAM |

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
    --other-aligners minimap2 gapmm2 uLTRA deSALT \
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
