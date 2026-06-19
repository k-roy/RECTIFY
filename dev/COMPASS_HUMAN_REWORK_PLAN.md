# COMPASS human-A549 rework plan (2026-06-18)

Goal: run full per-read COMPASS arbitration on SG-NEx A549 Illumina short reads, back-propagate RECTIFY's
family-concordance gate, and re-validate the 111/609 with a SENSITIVE multi-aligner detector (STAR-alone
is only 14.5% sensitive on known positives — STEP 0). Work happens INSIDE the COMPASS repo
(`~/work/COMPASS`, cloned from k-roy/COMPASS; NOT pushed to the public repo without PI OK). Compute on Sherlock.

## What the codebase map established (grounding for this plan)
- COMPASS core (`compare_splice_junctions_from_multiple_aligners.py`) IS a full per-read arbitration:
  edit-distance score per aligner (end-to-end, soft-clip disabled → cross-aligner-comparable), keep best
  `< MAX_EDIT_DIST_TO_CONSIDER=8`, select by ungapped>gapped → annotated>unannotated → shortest-intron →
  no-junction-proximal-mismatch → seeded RNG (lines 549-571).
- **Ambiguity normalization ALREADY exists** (`adjust_ambiguous_junctions`, COMPASS_functions.py:112-184)
  → snaps to annotated-or-signal-optimal within the slide window. RECTIFY's canonical-within-window is a
  *refinement of the tie-resolution*, not a missing capability.
- **Cross-family concordance is the genuine MISSING piece** — COMPASS records per-aligner agreement
  (`junction_counts`, lines 660-679) but NEVER gates on it (pure per-read best-pick). This is the primary
  back-propagation payload.
- Cross-sample combine module named in the README **does not exist**; the de-facto one is
  `compare_individual_aligners.R` (QC filter + plots).

## Phase 0 — prerequisites
- ✅ **Strandedness = RF/dUTP** (0.999 on existing BAM) → matches COMPASS default (COMPASS_functions.py:216). No flip.
- Aligner availability (from scoping, VERIFIED): STAR 2.7.x, GSNAP/gmap, BBMap, minimap2 splice:sr present
  on Sherlock; HISAT2 via Lmod module; **Magic-BLAST ABSENT** (the weakest + heaviest member → DEFER it).
- Build the COMPASS conda env on Sherlock (env yml pins bbmap=38.84/hisat2=2.2.1/star=2.7.10a/magicblast/
  gmap=2021.05.27/cutadapt/samtools=1.7/picard/pysam) + **samfixcigar (jvarkit, separate install)**. Or
  adapt to the existing aligner_bench env where versions allow.
- Pull A549 Illumina FASTQ (~25GB, 3 reps) from S3 (pattern proven in pull.sbatch).
- Build a GRCh38 introns TSV in COMPASS's format (`chrom start end strand type Name intron_type`) from
  gencode v44 — chrom names MUST match the FASTA (use chr-named GRCh38 to avoid the chr5↔chrV / 5↔chr5 trap).

## Phase 1 — yeast→human adaptation (code edits in ~/work/COMPASS)
- `COMPASS.sh`: GENOME_VERSION/FASTA/GFF/GTF/INTRONS_FILE → GRCh38; READ_LENGTH → A549 (~150);
  fix the `analyze_exonic_and_intronic_sequence.py` → `_elements.py` filename mismatch (line 104).
- `process_reads_and_align.sh`: **MAX_INTRON_LENGTH 2000 → 500000** (BBMap maxindel/pairlen, STAR/HISAT2
  --max-intronlen — this is where the cap is EFFECTIVE; the Python arg is inert); cutadapt — DROP the
  Roy-prep `-g "T{100}"` / `-A "A{100}"` polyT/polyA arms (A549 is standard TruSeq), keep AGATCGGAAGAGC.
- Human SS consensus/penalties (compare-script 140-144) — retune for human; the human-oriented SS-pair
  list already in `analyze_exonic_and_intronic_elements.py` (lines 98-127) is a template.
- `compare_individual_aligners.R`: remove yeast-locus filters (rRNA chrXII, LSR1 chrII, HIS3/etc.,
  `chrom!='chrmt'`, `intron_size<2000`) → human equivalents (rRNA repeats, snRNA, chrM).
- Branchpoint module (`analyze_exonic_and_intronic_elements.py`): yeast CUNAN → human degenerate yUnAy,
  widen MAX_DISTANCE_BETWEEN_3SS_AND_BP. (Lower priority — not needed for the 111 re-validation.)

## Phase 2 — multi-aligner alignment (heavy, chunked SLURM on Sherlock owners + AVX-512 constraint)
- Build human indices once: STAR (~30GB RAM at build; sjdbOverhang=READ_LENGTH-1), HISAT2 (--ss/--exon),
  GSNAP (gmap_build + gtf_splicesites + iit_store), Magic-BLAST (makeblastdb), BBMap (auto). Idempotent skip-checks.
- Align 3 reps × panel, end-to-end/no-softclip, numbered reads, samfixcigar→SAM1.4, name-sort. Chunk per
  (rep × aligner). **FULL 6-aligner COMPASS panel (PI decision 2026-06-18 — match COMPASS spec):** BBMap,
  STAR default+noncanonical, HISAT2 default+noncanonical, Magic-BLAST, GSNAP — + minimap2 splice:sr
  (RECTIFY addition). GSNAP/HISAT2/Magic-BLAST SIMD → sse42/nosimd to avoid SIGILL on AVX-512-excluded
  nodes; smoke-test one chunk before the array. Magic-BLAST on the 3GB human genome is the heaviest member
  (cap `-num_threads 12` for memory) — budget for it; it is the slow pole, not deferred.
- **Push authorized (PI 2026-06-18):** the human-adapted COMPASS may be pushed to github.com/k-roy/COMPASS
  (no active watchers) — use a branch (e.g. `human-a549`), not a force-push to master.

## Phase 3 — per-read arbitration
- Run `compare_splice_junctions_from_multiple_aligners.py` (human-adapted) → per-rep COMPASS junction TSV.

## Phase 4 — back-propagate RECTIFY's family-concordance gate (the genuine payload)
- At `junction_counts` (lines 660-679): map aligners → FAMILIES (star=STAR×2, hisat2=HISAT2×2, gsnap,
  bbmap, minimap2) so modes don't double-vote; emit a family-vote count alongside COMPASS_counts. A NOVEL
  (non-annotated) junction requires ≥N independent families. Mirror an R-side threshold for cross-sample.
- (Optional) refine `adjust_ambiguous_junctions` tie-resolution with canonical-within-window.

## Phase 5 — re-validate the 111/609 (settle the inconclusive verdict)
- Intersect the human-COMPASS junction set (family-gated, ambiguity-normalized) against the 111 / 609.
- CALIBRATED controls (the FDR fix — concordance alone isn't FDR-proof; aligners can share the short-anchor-
  over-GT-AG failure mode): positives = annotated gencode junctions in the same expressed loci (the panel
  must validate these HIGH — fixing the 14.5% STAR-alone failure); negatives = gmap_noncanonical sample +
  a shuffled-junction null of matched intron length. Tune the family gate so negatives stay ~0.
- DECISION: if positives(expressed) HIGH and negatives ~0 → trustworthy; then 111 ~0 → confirmed artifacts
  (finalize keep-GMAP-behind-fences); nontrivial fraction validate → real GMAP-unique novels (revise
  Deliverable B). If positives stay low even with the full panel → still inconclusive → escalate to the
  P0 sim benchmark + the C6 variant-aware check on the clustered loci.

## Critical path & risks
- Critical path = Phase 1 (config) → Phase 2 (index builds + alignment, the heavy/long compute) → Phase 5.
- Magic-BLAST deferred (weakest + heaviest on a 3GB genome). minimap2 splice:sr added (cheap, available).
- Biggest risk = compute scale (human STAR/Magic-BLAST memory; ~25GB FASTQ → ~10× intermediate disk).
- Concordance ≠ FDR control → the calibrated controls in Phase 5 are mandatory, not optional.
