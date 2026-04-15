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
│   │   ├── sherlock_larsms.yaml  # Standard CPU partition (use_scratch, streaming on)
│   │   └── sherlock_gpu.yaml     # GPU partition
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
│   ├── data/
│   │   └── bin/
│   │       └── linux_x86_64/
│   │           └── deSALT            # Vendored deSALT v1.5.6 binary (773 KB)
│   └── ...
└── CLAUDE.md                         # This file
```

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

**No open bugs remaining** — see `docs/BUGS_TO_FIX.md`.

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
- Bundled BAMs renamed: `rectified.bam` → `rectified_pA_hardclip.bam`; `rectified_softclip.bam` → `rectified_pA_softclip.bam`.

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
