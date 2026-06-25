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
| `--uge-queue` | — | UGE/SGE queue name (`-q`); omit on Hoffman2 campus jobs |
| `--uge-pe` | `shared` | UGE/SGE parallel environment (`-pe`; Hoffman2 = `shared`) |
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

With `--generate-slurm` (correct-first long-read pipeline):

```text
output_dir/
├── ...chunk files...
├── chunks_manifest.json
├── run_array_mapPacBio.sh          # mapPacBio alignment array (omitted with --skip-map-pacbio)
├── run_array_others.sh             # minimap2/gapmm2/uLTRA/deSALT alignment array
├── run_merge_aligners.sh           # merge per-aligner BAMs
├── run_prescan.sh                  # variant + junction prescan
├── run_array_correct_<aligner>.sh  # per-aligner correction arrays
├── run_array_chunk_merge.sh        # per-chunk merge → consensus
├── run_final_merge.sh              # final merge
├── submit_pipeline.sh              # submits the whole DAG with dependencies
└── logs/                           # log directory
```

The simplest way to run everything is `bash output_dir/submit_pipeline.sh`,
which submits each stage with the correct inter-job dependencies.

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

### Step 2 — Submit the pipeline

```bash
bash /scratch/chunks/submit_pipeline.sh
```

`submit_pipeline.sh` submits each generated stage (alignment arrays → merge →
prescan → per-aligner correction arrays → chunk-merge → final merge →
consensus) with the correct inter-job dependencies, following the correct-first
ordering. The alignment array tasks run `rectify align --no-consensus`, writing
one sorted BAM per (chunk, aligner) pair; the correction and merge stages then
produce the final consensus BAM.

You can also submit the individual array scripts by hand (e.g.
`sbatch run_array_others.sh`) if you need finer control.

---

## Scheduler compatibility

The generated array scripts use `SLURM_ARRAY_TASK_ID` and `SLURM_CPUS_PER_TASK`.
Generate non-SLURM headers directly with `--scheduler {uge,pbs}`, or add a shim
at the top of an array script body:

| Scheduler | Task ID | CPUs |
|-----------|---------|------|
| SLURM | `SLURM_ARRAY_TASK_ID` (0-based) | `SLURM_CPUS_PER_TASK` |
| UGE/SGE | `SGE_TASK_ID` (1-based) | `NSLOTS` |
| PBS/Torque | `PBS_ARRAY_INDEX` | `PBS_NUM_PPN` |

**UGE/SGE shim**:

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

## Read-Number Sidecar

`rectify split` adds an `RN:i:<read_num>` FASTQ comment tag to every derived
chunk record and writes `<sample>.read_num_sidecar.parquet` beside the chunks.
The sidecar maps RN back to the original FASTQ QNAME and full FASTQ comment,
including cDNA metadata tags. Existing comments are preserved after the RN tag.

Old chunk FASTQs and BAMs without RN still work; consensus falls back to the
legacy normalized-QNAME merge path whenever any input BAM lacks RN.
