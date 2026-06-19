# HANDOFF — RECTIFY short-read (COMPASS) mode (2026-06-19)

Implements `dev/RECTIFY_SHORTREAD_MODE_PLAN.md` P1–P4 on branch `drs-validation-rebuild`.
End goal: a sensitive paired multi-aligner short-read chr5 junction set for A549 to
adjudicate the **111 GMAP-only recurrent novels** (`dev/gmap_only_recurrent_novels_chr5.tsv`).

## Done
- **P1 — configurable tiebreak** (`core/consensus/select.py`): `select_best_alignment(tiebreak=...)`,
  COMPASS order `ungapped > gapped > annotated > shorter-intron`. Function default stays `rectify`
  (long-read unchanged); short-read entrypoints pass `compass`. `rectify consensus --tiebreak`,
  and `align --short-read` auto-selects compass. Threaded through `run_consensus_selection`.
- **P2 — paired chunking** (`core/commands/split_command.py::split_fastq_paired`): round-robin over
  PAIRS → shared RN, same-chunk mates, paired chunk FASTQs (`_chunk_NNN_of_MMM_R1/_R2`), identical
  bare QNAME per mate. `rectify split --read2` wired. Mate-desync / count-mismatch raise.
- **P3 — COMPASS aligner panel** (`core/align/multi_aligner.py`): `run_star`, `run_hisat2`,
  `run_magicblast`, `run_gsnap` (paired) + paired `run_bbmap`; invocations ported verbatim from
  Sherlock `…/compass_a549/COMPASS/process_reads_and_align.sh`. Index paths derived to match the
  prebuilt `genome_references/` layout. Threaded `reads2_path`/`read_length` through
  `run_multi_aligner` + `align_command`. Panel = **COMPASS-7** (user-confirmed: keep non-canonical).
- **P4 — paired array pipeline** (`split_command.py`): per chunk, one `rectify align --short-read
  --read2` = panel + compass consensus, **NO correction**; per-chunk rectified BAM → `*.consensus.bam`
  → shared final-merge. `--generate-slurm` wired; idempotent skip; atomic chunk copy.
- **Adversarial review fixes** (commit `2e0388c`):
  - **HIGH** paired-consensus mate-drop fixed (merge now emits one group per mate).
  - bbmap `maxindel=200000 pairlen=200000` for paired; atomic `cp`; Magic-BLAST gz decompress.

## Verified (locally)
- `pytest` green for: tiebreak (discriminating compass-vs-rectify), paired chunking + round-trip,
  paired-consensus both-mates-survive, COMPASS command builders (exact flags vs source), and the
  touched-module regressions. Latest focused run: 62 passed, 4 skipped.
- `rectify split --read2 … --generate-slurm --short-read` produces R1/R2 chunks + sidecar + manifest
  + 3 scheduler scripts; all pass `bash -n`.
- Sourced exact COMPASS params from Sherlock (not memory) → added gz adaptations
  (STAR `--readFilesCommand zcat`, GSNAP `--gunzip`, Magic-BLAST decompress).

## NOT VERIFIED (the RUN — this is your job, P5)
- **P3/P4 RUN**: no STAR/HISAT2/magicblast/gsnap binaries or human indices exist on the M1, so the
  aligners have never executed — only their assembled command lines are unit-tested. Running them on
  A549 is P5 and must happen on Sherlock. Everything below de-risks that run.

## Pre-flight — COMPLETE (2026-06-19, read-only SSH; you do NOT need to re-check)
All confirmed on Sherlock under `$REF=$W/COMPASS/genome_references`, env
`/home/groups/larsms/users/kevinroy/anaconda3/envs/compass`:
- ✅ `GRCh38_gencode_v44.fasta` (+ `.fai`) and `.gtf` (symlink) present. The code's
  `_genome_version`/`_compass_index_paths` derivation matches the on-disk layout **exactly**:
  HISAT2 `HISAT2_annotated_index/GRCh38_gencode_v44`, BLAST `BLAST/GRCh38_gencode_v44`,
  GSNAP `-D GSNAP -d GRCh38_gencode_v44`, STAR `STAR_annotated_150_bp_SJDB_index`, bbmap `bbmap/`.
- ✅ **HISAT2 `_splice_sites.txt` sidecar present** (10 MB) alongside 8 `.ht2` files + `_exons.txt`.
- ✅ All 7 binaries in the `compass` env: STAR, hisat2, magicblast, gsnap, bbmap.sh, samtools.
- ✅ BLAST (`.nhr/.nin/.nsq/…`), GSNAP (`GRCh38_gencode_v44/` + `.iit` + `.splicesites`) indices present.
- ⚠️ **STAR version: index `versionGenome 2.7.4a`, binary `2.7.10a`.** Same 2.7.x genome format, so
  2.7.10a should read it. If STAR aborts with `versionGenome … is INCOMPATIBLE`, regenerate the index
  (`STAR --runMode genomeGenerate … --sjdbOverhang 149`, ~30 min) — otherwise no action.
- ⚠️ **magicblast 1.5.0** (Dec 2019) — gz support at this version is uncertain; this is exactly why
  `run_magicblast` decompresses the gz chunk to scratch first (commit `2e0388c`). No action needed.

**Net: pre-flight is GREEN. Use the `compass` env. The only thing that could bite is the STAR
genome-version warning above — and that's a 30-min regenerate if it even triggers.**

## Open (not blocking P5, noted by review)
- If R1 and R2 pick different best aligners, the output BAM's pair flags (RNEXT/PNEXT/proper-pair)
  are inconsistent. Harmless for junction extraction (the end goal); the BAM is not a clean paired
  BAM for other tools.
- `_strip_mate_suffix` strips only `/1`,`/2` (NOT `.1`/`.2` — that would collide with SRA read
  numbering `SRR.1`/`SRR.2`, a worse bug). A549 raw is `_R{1,2}` + Casava headers → fine.
- Junction aggregation/extraction from the merged consensus BAM is a P5 step (no change made here).

## Resume — TO THE FOLLOW-UP AGENT picking up P5
You are running the short-read COMPASS pipeline on A549 and adjudicating the 111. The code is done,
committed, and locally tested; **pre-flight is already done (see above) — do not repeat it.** Context:
`$W=/scratch/users/kevinroy/compass_a549`, Sherlock `larsms`, env `compass`. Never tear down the SSH
ControlMaster (re-2FA cost); retry transient sshd errors serially. The pipeline aligns to the FULL
genome (`GRCh38_gencode_v44.fasta`); the chr5 restriction is a step-5 junction filter, NOT an
alignment input. Your steps:

1. Get this branch (`drs-validation-rebuild`, through commit `2e0388c`) onto wherever RECTIFY runs on
   Sherlock, and make the array body's `$PYTHON` resolve to the `compass` env (so all 7 binaries +
   samtools 1.17 + pysam are on PATH inside each task). Confirm with `STAR --version` → `2.7.10a`.

2. SMOKE FIRST (do NOT launch the 500-chunk run blind). Generate + submit a 4-chunk job on a 50k-pair
   subset and confirm all 7 aligners run, the RN-keyed merge keeps BOTH mates, and consensus emits:
   ```
   rectify split $W/COMPASS/fastq/A549_rep1_R1.fastq.gz -2 $W/COMPASS/fastq/A549_rep1_R2.fastq.gz \
     --short-read --read-length 150 -n 4 -o $W/rectify_sr_smoke \
     --genome   $W/COMPASS/genome_references/GRCh38_gencode_v44.fasta \
     --annotation $W/COMPASS/genome_references/GRCh38_gencode_v44.gtf --generate-slurm
   bash $W/rectify_sr_smoke/submit_pipeline.sh        # array (7 aligners/chunk) → final merge
   ```
   ACCEPTANCE for the smoke: per-chunk log shows all 7 aligners succeeding (watch the STAR
   `versionGenome` line — regenerate the index only if it errors INCOMPATIBLE); the merged BAM has
   **two records per RN** (both mates survived the consensus fix — this is the thing the adversarial
   review caught, so verify it explicitly with e.g. `samtools view … | grep RN:i: | …`); and N-op
   junctions are present. If any aligner silently produces 0 reads, check its index/env before scaling.

3. FULL A549 rep1 (~42M pairs → ~500 chunks): same command without `-n 4` (auto-size). Submit on
   `larsms`; tasks are 20–60 min and idempotent (`.consensus.bam` skip + atomic copy make requeue
   safe). Drop a sentinel + refresh THIS handoff while it runs (per `~/work` discipline).

4. VALIDATE THE 111 (the actual deliverable; `dev/RECTIFY_SHORTREAD_MODE_PLAN.md` §P5): from the merged
   consensus BAM, extract chr5 junctions. (a) POSITIVE control: annotated chr5 junctions in expressed
   loci must score HIGH — if they don't, the panel/consensus is broken, stop and debug. (b) NEGATIVE
   control ≈ 0. (c) Then intersect the 111 (`dev/gmap_only_recurrent_novels_chr5.tsv`) with the COMPASS
   junction set: **`111 ∩ COMPASS ≈ 0` ⇒ the 111 are confirmed artifacts** (the whole point). Report
   the three numbers together — a near-zero intersection only means "artifact" if the positive control
   passed (otherwise it just means the panel found nothing).

If you change the aligner panel or params, the source of truth is the COMPASS shell at
`$W/COMPASS/process_reads_and_align.sh`; the wrappers in `core/align/multi_aligner.py` were ported from
it verbatim (with gz adaptations) — keep them in sync.

## Files
- Code: `core/consensus/{select,consensus}.py`, `core/align/multi_aligner.py`,
  `core/commands/{split_command,align_command,consensus_command}.py`.
- Tests: `tests/test_{consensus_selection,chunker_paired,compass_aligner_cmds,paired_consensus_merge}.py`.
- Commits: `ba06e4f` (P1) `7132cb3` (P2) `c65ed2c` (P3) `646ddc1` (P4) `2e0388c` (review fixes).
- Source of truth for aligner params: Sherlock `…/compass_a549/COMPASS/process_reads_and_align.sh`.
