# RECTIFY short-read (COMPASS) mode — IMPLEMENTATION PLAN (2026-06-19)

**DECISION (user, 2026-06-19): Option B.** Build the short-read multi-aligner splice-junction caller as a
**mode inside RECTIFY**, skipping 5'/3' end-correction. The COMPASS algorithm survives as RECTIFY's existing
consensus module (`core/consensus/select.py`) — stays citable; we shed the standalone yeast-era COMPASS
scaffolding (renumbering, sequential alignment, no chunking, conda/samtools/introns/bbmap grief).
**Tiebreak order = configurable, DEFAULT to COMPASS published order** (`ungapped > gapped > annotated >
shorter-intron`); allow override. This REVERSES the earlier "RECTIFY = end-correction only" scope call —
RECTIFY is now deliberately a two-modality (long-read end-correction + short-read junction) tool.

**End goal (unchanged):** a sensitive multi-aligner short-read chr5 junction set for A549 to adjudicate the
**111 GMAP-only recurrent novels** (`dev/gmap_only_recurrent_novels_chr5.tsv`; currently "likely artifacts,
INCONCLUSIVE" because STAR-1-pass is only 14.5% sensitive). Validate: positives = annotated chr5 junctions in
expressed loci (must score HIGH), negatives ~0, then 111 ∩ COMPASS ≈ 0 → confirmed artifacts.

## Why B (audit, cross-validated by 2 independent agents + direct code read)
- **Arbitration already exists** = `core/consensus/select.py::select_best_alignment` (line 19): annotation-blind
  `score_alignment` → max-score → tiebreakers (`_tiebreak_key`, line 81-91). Junction pool = annotated ∪ all
  aligners' observed (line 45-58). This IS COMPASS's per-read best-alignment design.
- **Short-read path half-built** = `--short_read` flag real end-to-end (`run/stages.py:48,107-120`): switches
  base aligners to `['bbmap','bwa']`, junction aligners to `[]`, calls `run_align`; correction is a SEPARATE
  skippable stage (`stages.py` docstring; `_run_correction` independent of `_run_alignment`).
- **RN sidecar** = `core/chunking/sidecar.py`: additive `RN:i:<int>` tag in FASTQ comment (original QNAME
  PRESERVED) + Parquet sidecar (`read_num,original_qname,fastq_comment,chunk_id,seq_md5,qual_md5`; zstd-3,
  100k row groups). Propagated into BAMs (`multi_aligner.py:199-227 _inject_rn_into_bam`), consumed by an
  **RN-keyed K-way merge** (`consensus.py:300-312`) that natively tolerates aligners having DIFFERENT read
  subsets → moots the chr5-filter desync entirely.
- **Chunking/arrays** = `split_command.py`: round-robin FASTQ split (`DEFAULT_TARGET_READS=50_000`,
  `MIN_CHUNKS=4`, `MAX_CHUNKS=500`, clamp; 42M reads→500 chunks). **chunk×aligner grid array**
  (`_others_array_body:636-743`, `n_tasks=n_chunks*n_aligners`, `ALIGNER_IDX=TASK_ID/N_CHUNKS`,
  `CHUNK_IDX=TASK_ID%N_CHUNKS`). Dependency-chained sbatch generation, idempotent `.done` skips, scheduler-
  agnostic task-ID shim, scratch staging + rsync-back. CPU/thread portability = `slurm.py:49-142`.

## ⚠️ TWO REAL WORK ITEMS (neither a blocker)
1. **Paired-end adaptation (RECTIFY infra is SINGLE-END).** Verified: no R2/mate/pair handling in
   `core/chunking`, `split_command.py`, `multi_aligner.py`. Round-robin `read_idx % n_chunks` would scatter
   R1/R2 into different chunks with different RNs → breaks proper-pairs. **Fix:** chunk by PAIR index; share
   ONE `RN:i:` across both mates; emit paired chunk FASTQs (`_chunk_NNN_R1/_R2`); paired aligner invocations
   (`bbmap in1=/in2=`, `hisat2 -1/-2`, `STAR --readFilesIn R1 R2`, etc.); shared-RN keys consensus. This is the
   bulk of the new code.
2. **Extend the short-read panel** from `bbmap+bwa` to the full COMPASS set for sensitivity: add STAR_default,
   STAR_noncanonical, HISAT2_default, HISAT2_noncanonical, Magic-BLAST, GSNAP wrappers to `multi_aligner.py`
   (paired invocations + RN injection). Reuse the exact human params from the COMPASS work: end-to-end /
   `--no-softclip` / `alignEndsType EndToEnd`, `MAX_INTRON_LENGTH=200000` (bbmap `maxindel=200000 pairlen=200000
   intronlen=20`; STAR/HISAT2 `--max-intronlen 200000`), `--rna-strandness RF`. **CHECK:** whether RECTIFY's
   `scoring.score_alignment` needs CIGAR `=/X` (COMPASS needed samfixcigar) or uses NM/MD — if NM/MD, the
   samfixcigar step is UNNEEDED (one less thing to port). Read `core/consensus/scoring.py` + `extract.py` first.

## TIEBREAK CHANGE (small, localized — user's explicit decision)
`core/consensus/select.py` `_tiebreak_key` (line 81-91) currently returns
`(_count_3prime_agreement, n_annotated, canonical_count)`. The 3'-agreement term NO-OPS in short-read mode
(`corrected_3prime` is None without end-correction). **Make tiebreak order configurable; default to COMPASS
published order** `ungapped > gapped > annotated > shorter-intron`. Need from `AlignmentInfo` (read `extract.py`):
a gapped/ungapped flag (has-intron) and total intron length. Add a `tiebreak='compass'|'rectify'` param threaded
from the short-read entrypoint. Keep RECTIFY's order available for long-read.

## PHASED PLAN
- **P0 — recon (read before coding, ~30 min):** `core/consensus/{select,scoring,extract,consensus}.py`,
  `core/chunking/sidecar.py`, `core/align/multi_aligner.py` (aligner wrappers + RN inject points),
  `core/commands/{align_command,split_command,run/stages}.py`. Confirm: (a) does scoring need CIGAR =/X?
  (b) exact AlignmentInfo fields for the COMPASS tiebreak; (c) where base_aligners list is consumed so new
  aligners slot in.
- **P1 — tiebreak config** (smallest, self-contained): make `_tiebreak_key` order a param, default COMPASS;
  unit-test on a synthetic tie. Branch off `drs-validation-rebuild`.
- **P2 — paired-end RN/chunking:** pair-aware `split_fastq` (shared RN, same-chunk mates, paired chunk FASTQs)
  + sidecar schema note (mate in `fastq_comment` or a `mate` col). Test on a 10k-pair subset: RN identical
  across mates, both in same chunk, round-trip reconstruct.
- **P3 — aligner panel:** add the 6 COMPASS aligners (paired) to `multi_aligner.py` short-read path with the
  human params above + RN injection; reuse the already-built Sherlock indices (see DATA below).
- **P4 — wire + run:** short-read array pipeline (paired chunks × N aligners → RN-merge → consensus → junction
  aggregation, SKIP correction) via `split_command` generator; run on A549 rep1.
- **P5 — validate the 111:** intersect chr5 consensus junctions with the 111 + calibrated controls.

## DATA / ENV STATE (Sherlock, all PRESERVED — `$W=/scratch/users/kevinroy/compass_a549/`)
- **Raw paired FASTQ (USE THIS, not the renumbered copy):** `$W/COMPASS/fastq/A549_rep1_R{1,2}.fastq.gz`
  (~42M pairs, 3.5/3.7 GB). reps 3 & 5 NOT staged yet.
- **DISCARD/ignore** `$W/COMPASS/numbered_reads_fastq/A549_rep1_numbered_R*.fastq` — COMPASS's destructive
  `@N_R1` renumber; superseded by RN sidecar. (Safe to delete to reclaim ~26 GB.)
- **Prebuilt human indices reusable** at `$W/COMPASS/genome_references/`: bbmap/ (clean-FASTA, headers fixed),
  BLAST/, GSNAP/, HISAT2_annotated_index/, STAR_annotated_150_bp_SJDB_index/. Clean genome
  `GRCh38_gencode_v44.fasta` (real file, `>chr5`-style headers, faidx'd), `.gtf`, chr5 introns
  `GRCh38_gencode_v44_chr5_introns.tsv` (18,451) + chr5 FASTA `GRCh38_gencode_v44_chr5.fasta`.
- **conda env `compass`** (`/home/groups/larsms/users/kevinroy/anaconda3/envs/compass`) has all 7 aligners +
  samtools 1.17 + pysam. RECTIFY's own env may differ — confirm which env the RECTIFY short-read mode runs in.
- **Cluster:** Sherlock `larsms` partition (256GB non-preempt AMD Milan, no AVX-512 needed). conda solves on
  compute nodes only. NEVER tear down the SSH ControlMaster.

## NON-NEGOTIABLES carried over
BBMap stays in the panel (2023 NAR workhorse). Chunked idempotent arrays (this whole pivot). RN sidecar (not
renumbering). Reuse RECTIFY's proven infra, don't re-solve. Tiebreak defaults to COMPASS for the 111.
