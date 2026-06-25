# RECTIFY — Developer & Agent Context

This file is the **index and hot-path rules** for the RECTIFY codebase.
Per-topic deep-dives live in `docs/` — fetch on demand, don't preload.

---

## READ FIRST — canonical references

| Topic | File |
| --- | --- |
| Chunking + pipeline order + scratch staging + run-all + HPC I/O + thread limits | `docs/architecture/pipeline_and_io.md` |
| Multi-sample streaming pipeline (manifest mode, position index, DESeq2) | `docs/architecture/multi_sample_pipeline.md` |
| Protocol flag `--dT-primed-cDNA` + module activation by protocol | `docs/protocols/dt_primed_cdna.md` |
| gapmm2 DRS aux-tag / dup-UUID gotchas + PAF→BAM sequence injection | `docs/aligners/gapmm2.md` |
| Minimap2 long-read RNA flags | `docs/aligners/minimap2.md` |
| Empirical penalty tables (`--junction-penalty-table`) quick reference | `docs/penalty_tables_quickref.md` |
| Empirical penalty tables — algorithm + design rationale | `docs/EMPIRICAL_HP_PENALTY_SCORING.md` |
| Per-read validation plotting (junction zoom, pA-rest, HTML reports) | `scripts/validation_data/PLOTTING.md` |
| Pre-public version history (v2.x / v3.x — internal) | `dev/audits/audit_history.md` |
| Current session handoff (open bugs, in-flight design) | `HANDOFF.md` |
| Algorithm walk-throughs (per-module) | `docs/algorithms/*.md` |
| Validation read bundle layout + per-read expectations | `rectify/data/validation/VALIDATION_READS.md` |

---

## Critical rules (summary — fetch the linked doc for detail)

### Chunking is mandatory for FASTQ on HPC

NEVER `rectify run-all <fastq>` without `--chunked-alignment`. Use
`rectify split` for short-read or `--chunked-alignment` for long-read.
BGZF BAM writes to Oak NFS are ~75,000× slower than `$SCRATCH`. See
**`docs/architecture/pipeline_and_io.md` § Chunking**.

### Pipeline order: correct-first, never consensus-first

```
align → correct (per aligner, --aligner-bams [all 5]) → merge → consensus
```

Consensus-first ordering is WRONG and produces worse aligner selection.
`split_command.py` was fixed to correct-first 2026-04-30. See
**`docs/architecture/pipeline_and_io.md` § Pipeline order**.

### Scratch staging: long-read array scripts need manual patching

Generated long-read array scripts write to Oak NFS by default — that
melts the cluster. Patch every `run_array_*.sh` to write to
`$SCRATCH_DIR` then `rsync` back. See
**`docs/architecture/pipeline_and_io.md` § Scratch staging**.

### Always use `--streaming` for `rectify correct` in SLURM

Without `--streaming`, a 40M-read BAM peaks at 30–40 GB RAM. With it,
~4–5 GB regardless of size. `hpc_cpu.yaml` profile sets `streaming:
true` by default.

### Set thread limits BEFORE importing numpy

SLURM bans accounts that spawn more processes than allocated CPUs.
`set_thread_limits()` in `slurm.py` must run before any
numpy/sklearn/pydeseq2 import. **MUST** set `LOKY_MAX_CPU_COUNT` (loky
ignores `JOBLIB_WORKERS`).

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

`rescue_3ss_truncation()` in `splice_aware_5prime.py` calls
`align_clip_to_exon()` from `local_aligner.py` to compute a proper M/I/D
CIGAR for the exon segment instead of emitting a flat `nM` block.

**Result dict**: `five_prime_exon_cigar` (SAM CIGAR string, e.g.
`"8M1D3M"`). Stored in `corrected_reads.tsv` as the
`five_prime_exon_cigar` column.
`bam_writer.extend_read_5prime_for_junction_rescue()` uses it when
writing Cat3 reads.

**Aligner design** (`local_aligner.py`):
- `_align_right_anchored(query, ref)` — free prefix in ref, right end
  fixed. Plus strand: clip must end at `intron_start`.
- `_align_left_anchored(query, ref)` — left end fixed, free suffix.
  Minus strand: clip must start at `intron_end`.
- **Affine gap (Gotoh 1982)**: match=+2, mismatch=-4, gap_open=-4,
  gap_extend=-1. Affine gap prevents the "staircase" artifact where
  many isolated 1-base deletions outscore a single consolidated
  deletion.
- Buffer = `max_indel=5` bp added to each side of the expected exon
  window.

**Chimeric stitch fix**: `chimeric_consensus.build_chimeric_cigar()`
uses `D` (deletion, op=2) instead of `N` (intron skip, op=3) for
reference gaps ≤ 10 bp at segment boundaries. Larger gaps stay `N`.

**Resolved (2026-06): ** the downstream exon-CIGAR mangling
(`14=1D9=` → `22D1M21I1M`) and the materializer `intron_end`/acceptor-short
issue were fixed by the Cat3-rescue commits (`a11ef8e`, `00297d7`,
`91d8336`). The cat3 fixtures in `test_validation_reads.py` pass; the
suite is green. See "Open work (live)" below.

---

## File layout (top-level)

```
rectify/
├── rectify/
│   ├── slurm.py                       # CPU detection, thread limits, scratch utilities
│   ├── slurm_profiles/                # hpc_cpu.yaml (streaming on), hpc_gpu.yaml
│   ├── core/
│   │   ├── commands/                  # CLI command implementations
│   │   │   ├── correct_command.py
│   │   │   ├── run/ (single_sample.py, multi_sample.py, stages.py, chunked_batch.py)
│   │   │   ├── batch_command.py       # generates SLURM array scripts
│   │   │   ├── split_command.py       # rectify split — FASTQ chunker
│   │   │   ├── consensus_command.py
│   │   │   ├── restore_polya_command.py
│   │   │   └── install_aligners_command.py
│   │   ├── bam/                       # bam_processor, bam_writer (CIGAR surgery), read_edits
│   │   ├── consensus/                 # consensus.py, corrected_consensus.py, scoring, select, extract
│   │   ├── correct/                   # indel_corrector, walkback, etc.
│   │   ├── splice/                    # junction_refiner (Module 2H), junction_scoring, splice_aware_5prime, calibrate_junction_overhang
│   │   ├── align/                     # multi_aligner, local_aligner
│   │   ├── aggregate/                 # five_prime, three_prime
│   │   ├── analyze/                   # manifest streaming, clustering, DESeq2
│   │   └── ...
│   ├── data/                          # Bundled reference data — use --Scer
│   │   ├── S288C_reference_sequence_R64-5-1_20240529.fsa{,.fai,.gz}
│   │   ├── saccharomyces_cerevisiae_R64-5-1_20240529.{gff,gff.gz,gtf,junc.bed}
│   │   ├── genomes/saccharomyces_cerevisiae/penalty_tables/   # bundled DRS penalty tables
│   │   ├── validation/                # bundled validation reads + TSVs + per-aligner BAMs
│   │   └── bin/linux_x86_64/deSALT    # vendored deSALT binary
│   └── visualize/                     # Plotting (dev/PLOT_SKILLS.md)
├── scripts/validation_data/           # Per-read renderer + report generator (PLOTTING.md)
├── tests/
├── docs/                              # Algorithm + architecture documents
└── CLAUDE.md                          # This file
```

### Using bundled data directly

When `--Scer` (or `--organism saccharomyces_cerevisiae`) is passed,
data paths are resolved automatically. To reference explicitly:

```python
import rectify
from pathlib import Path
DATA = Path(rectify.__file__).parent / 'data'

GENOME     = DATA / 'S288C_reference_sequence_R64-5-1_20240529.fsa'
ANNOTATION = DATA / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff'
JUNC_BED   = DATA / 'saccharomyces_cerevisiae_R64-5-1_20240529.junc.bed'
```

---

## Project basics

- **Default branch**: `master` (not `main`).
- **Pre-1.0**: `0.9.0` is the first public release. Don't reference 2.x
  or 3.x in user-facing docs — those were internal pre-public; see
  `dev/audits/audit_history.md` for the history.
- **Test suite**: `pytest -m "not slow"` runs in ~1 minute, covers
  ~934 tests + 28 skipped + 4 deselected. The `slow` marker covers
  cDNA-pipeline smoke and chain canary (~5 min each).
- **Surgical staging**: users routinely have WIP in the working tree
  (especially `rectify/core/bam/*`, `rectify/core/correct/*`,
  `rectify/data/validation/*`, `scripts/validation_data/`). Always
  `git add <explicit paths>` — never `git add -A` or `git add .`.

---

## Open work (live)

See `HANDOFF.md` for the current session's state. The DRS validation
build on `drs-validation-rebuild` is **green** (`pytest -m "not slow"` =
1603 passed; 69 Cat3/Cat6/refiner tests + the cat1–cat9
`test_validation_reads.py` set all pass). The Cat3/Cat6 items that were
formerly in flight are **RESOLVED**:

- ✅ Materializer `intron_end`/acceptor-short, `reroute_intronic_tail_5prime_via_junction`,
  and the downstream exon-CIGAR mangling were fixed by the Cat3-rescue
  commits (`a11ef8e`, `00297d7`, `91d8336` on `splice_aware_5prime.py`).
- ✅ Cat6 cascade (1→9) closed; the cat1–cat9 fixtures pass.
- ✅ `merge_corrected_tsvs` settled on **hp_edit_distance** as the primary
  sort key (`corrected_consensus.py`), with the chrom-aware `_n_agree`
  tie-breaker as fallback. No HP-vs-legacy decision remains pending.

**Remaining (genuine, non-blocking):**

- **Future enhancement**: pre-compute up_amb / down_amb fields on
  annotated junctions so soft-clip rescue can flex the match length
  within the ambiguity window.
- **GLASS aligner eval** (the one open `dev/TODO.md` item): not installed
  on the cluster; evaluate before adding a wrapper. Winnowmap2 / Minisplice
  / GMAP are already wrapped.
