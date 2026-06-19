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

## NOT VERIFIED (cluster-only — needs Sherlock conda env + human indices)
- **P3/P4 RUN**: no STAR/HISAT2/magicblast/gsnap binaries or human indices locally. The aligners
  have never executed; only their assembled commands are tested. **This is P5.**
- Cluster pre-flight checks before the A549 run (from the adversarial review):
  1. `magicblast --version` — if recent, the decompress is belt-and-suspenders; either way OK now.
  2. Confirm `genome_references/HISAT2_annotated_index/GRCh38_gencode_v44_splice_sites.txt` exists
     (HISAT2 `--known-splicesite-infile`; the `.ht2` index existing doesn't guarantee this sidecar).
  3. Confirm the run uses the **`compass` conda env** (STAR refuses version-incompatible indices;
     `STAR_annotated_150_bp_SJDB_index` was built by COMPASS's STAR 2.7.x). Plan left "which env" open.
  4. `--read-length` matches the prebuilt STAR index (150 → `STAR_annotated_150_bp_SJDB_index`).

## Open (not blocking P5, noted by review)
- If R1 and R2 pick different best aligners, the output BAM's pair flags (RNEXT/PNEXT/proper-pair)
  are inconsistent. Harmless for junction extraction (the end goal); the BAM is not a clean paired
  BAM for other tools.
- `_strip_mate_suffix` strips only `/1`,`/2` (NOT `.1`/`.2` — that would collide with SRA read
  numbering `SRR.1`/`SRR.2`, a worse bug). A549 raw is `_R{1,2}` + Casava headers → fine.
- Junction aggregation/extraction from the merged consensus BAM is a P5 step (no change made here).

## Resume (P5 — run on Sherlock; `larsms`, env `compass`, `$W=/scratch/users/kevinroy/compass_a549`)
1. Sync this branch to where RECTIFY runs on Sherlock; activate the `compass` conda env.
2. Run pre-flight checks 1–4 above. If the HISAT2 splice-sites sidecar is missing, regenerate:
   `hisat2_extract_splice_sites.py $GTF > …_splice_sites.txt`.
3. SMOKE on a tiny subset first (e.g. 50k pairs) to confirm all 7 aligners run + RN-merge + consensus:
   ```
   rectify split $W/COMPASS/fastq/A549_rep1_R1.fastq.gz -2 $W/COMPASS/fastq/A549_rep1_R2.fastq.gz \
     --short-read --read-length 150 -n 4 -o $W/rectify_sr_smoke \
     --genome $W/COMPASS/genome_references/GRCh38_gencode_v44.fasta \
     --annotation $W/COMPASS/genome_references/GRCh38_gencode_v44.gtf --generate-slurm
   bash $W/rectify_sr_smoke/submit_pipeline.sh   # array (7 aligners/chunk) → final merge
   ```
   Verify the merged BAM has BOTH mates per read (RN appears on 2 records) and junctions present.
4. FULL A549 rep1 (~42M pairs → ~500 chunks): same command without `-n 4` (auto-size), on `larsms`.
   Chunk the array per `~/work` cluster discipline (20–60 min idempotent tasks; `.consensus.bam`
   skip makes requeue safe).
5. **Validate the 111** (`dev/RECTIFY_SHORTREAD_MODE_PLAN.md` P5): extract chr5 junctions from the
   merged consensus BAM; positives = annotated chr5 junctions in expressed loci must score HIGH;
   negatives ≈ 0; then `111 ∩ COMPASS` — if ≈ 0, the 111 are confirmed artifacts.

## Files
- Code: `core/consensus/{select,consensus}.py`, `core/align/multi_aligner.py`,
  `core/commands/{split_command,align_command,consensus_command}.py`.
- Tests: `tests/test_{consensus_selection,chunker_paired,compass_aligner_cmds,paired_consensus_merge}.py`.
- Commits: `ba06e4f` (P1) `7132cb3` (P2) `c65ed2c` (P3) `646ddc1` (P4) `2e0388c` (review fixes).
- Source of truth for aligner params: Sherlock `…/compass_a549/COMPASS/process_reads_and_align.sh`.
