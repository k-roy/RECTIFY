# HANDOFF — Deliverable A: the simulation ground-truth benchmark (the GATE)

**Agent:** dedicated benchmark-builder (isolated worktree, branch
`worktree-agent-a25a2c1e784ad37dc`, based on `drs-validation-rebuild` so the reuse
primitives + design docs are present). **NEVER commit to `drs-validation-rebuild`.**
**Updated:** 2026-06-26.

---

## DONE
Three ABSENT components built (everything else is reuse), all under the
benchmark-only paths the brief allows (`rectify/core/benchmark/`, `scripts/benchmark/`):

1. **Truth-propagation schema (built FIRST — it constrains the other two)** —
   `rectify/core/benchmark/truth_schema.py`. Per-read `ReadTruth` with:
   - **exact-indel truth** as a per-position `IndelTruth` carrying an
     ambiguity-equivalence span `[eq_start, eq_end)` (HP run / STR period / unique)
     — the net-indel-in-run model, so the framing metric is ambiguity-aware;
   - **junction truth** `JunctionTruth` left-normalized via
     `chimeric_consensus.normalize_junction`, with canonicity + **NIC/NNC class**
     (ANNOTATED / NIC / NNC) and the ambiguity window;
   - **CPA truth** + downstream-A context (C2); **standing-variant truth** (C6);
   - the **shared genomic-region-disjoint TRAIN/TEST split tag** (region-decided);
   - C1 HP/STR cell metadata (`base_class`, `run_copies`) for the `min_count=100` audit.
   - Lossless TSV round-trip (`write_truth_table`/`read_truth_table`). NOT overloaded
     onto the `XV` tag (a sidecar table, per the brief).

2. **Read-simulator wrapper (pbsim3)** — `scripts/benchmark/sim/`:
   - `transcript_model.py` — the truth-propagation backbone: `TranscriptModel`
     (exon structure -> spliced FASTA, `transcript_pos_to_genome` map, introns,
     NIC/NNC `junction_truths`). Junction/CPA/NIC-NNC truth is derived from THIS
     construction, NOT from the simulator.
   - `pbsim3_wrapper.py` — Tier-2 realistic reads. **pbsim3 chosen because it is
     the ONLY one of pbsim3/badread/nanosim that emits a per-read MAF
     (read<->template alignment)**, which the EXACT-INDEL framing metric requires;
     badread/nanosim give origin+identity but no per-read edit script. Parses the
     MAF and composes (read->transcript) ∘ (transcript->genome) into exact
     `ReadTruth` (genome CIGAR + indels + junctions). `--method errhmm`;
     ERRHMM-ONT model = DRS, ERRHMM-ONT-HQ = PCR-cDNA.
   - `controlled.py` — **Tier-1 controlled-error generators (the DISCRIMINATING
     tier)**: HP (A/C/G/T x 1-12, del-dominant), STR (ED-tied placements),
     JUNCTION_AMB (a GT-AG intron whose donor sits in a repeat → ambiguity test),
     CLEAN (false-indel control). Emits ref.fa + reads.fastq + truth.tsv.

3. **Ambiguity-aware per-junction + per-indel scorer** —
   `rectify/core/benchmark/scorer.py`. Reuses
   `chimeric_consensus.normalize_junction / junction_ambiguity_window /
   _canonical_within_window` so a call 1bp into a donor/acceptor repeat is NOT
   charged FP. Junction TP/FP/FN stratified by class (ANNOTATED/NIC/NNC) AND
   canonicity (separate FDR tracks, §8); indel **position-exact concordance** +
   false-indel rate stratified by C1 cell; CPA `|est-true|`; variant-adjacent FP
   tagging (C6); `score_panel` fills which-aligners-placed + panel-unplaced fraction.

4. **End-to-end smoke** — `scripts/benchmark/smoke_roundtrip.py`.

## VERIFIED — the gate is BUILT for the external-aligner path + the internal-DP
## ablation path is RUNNABLE. (C1-arm SEPARATION is the named remaining proof — see OPEN.)
- **GATE smoke is GREEN** (`scripts/benchmark/smoke_roundtrip.py`, minimap2 -ax
  splice on the Tier-1 corpus):
  - (A) known junction round-trips truth->reads->minimap2->scorer as **NNC TP=30**;
  - (A2) **NIC TP=30 and ANNOTATED TP=30** — the discovery-class classifier is
    verified end-to-end (not just present); separate FDR tracks exercised;
  - (B) a deliberately 1bp-shifted junction call (into the repeat) scores
    **TP=1 FP=0** (ambiguity-aware match, the GMAP-0.09 trap defended);
  - (B2) NEGATIVE control: a junction shifted 40bp (beyond the window) is correctly **FP**;
  - (C) indel **position-exact concordance = 0.987**, false-indel-rate=0.0;
  - (D) **internal-DP ablation path RUNS** — the live flat-affine DP
    (`align_exon_block_global`, the arm C1 upgrades) is BAM-ized
    (`scorer.cigar_records_to_bam`) and scored on 2400 HP_HARD reads
    (concordance 0.990). This is the 'ablations runnable' exit criterion AND the
    exact harness the future arm-LAW vs arm-flat comparison plugs into.
- **HONEST validity finding (advisor-driven):** minimap2 and the flat-affine DP are
  the SAME error family, so they AGREE by construction (boundary_sub = 1.000 ==
  1.000 after a truth-corruption fix where the boundary substitution was never
  allowed to become the run base). A flat-vs-flat smoke therefore CANNOT
  demonstrate C1 discrimination — a genuine C1 win is **arm-LAW vs arm-flat**, and
  the length-law arm is the NEXT cycle. The smoke reports this rather than falsely
  asserting flat-affine 'headroom' (which would be a false gate). isolated-HP =
  1.000 (the SPEC VERTICAL-SLICE FINDING), STR = 0.907.
- **Scorer correctness fix:** a per-indel TP now requires net-in-span == truth net
  AND no unexplained indel OUTSIDE every truth span (the vertical-slice
  ``out_run == 0`` rule); insertion-boundary off-by-one made consistent with the
  half-open span test.
- Schema lossless TSV round-trip verified (junction normalization, canonicity,
  HP-cell accounting, variant/split round-trip).
- pbsim3 MAF→genome **projection** verified locally on a synthetic MAF (2bp deletion
  + 1 intron correctly propagated to genome coords) — no live pbsim needed for this.

## OPEN
- **pbsim3 live run NOT yet executed.** Install on Sherlock env `aligner_bench`
  was still solving (slow conda solve on the large base env) at handoff time —
  the wrapper is code-complete + projection-validated but a real pbsim3 round-trip
  is unverified. (badread also requested in the same install for cross-check.)
- **Tier-1 cell sizing:** at `--reps 120` the min indel-bearing cell is ~43 reads
  (clean reads are ~50% by the del-dominant K_DIST; HP_HARD alternates two modes).
  The Sherlock scale-up must set `--reps` so every `(base_class, run_copies)` cell
  clears `min_count=100` (recommend `--reps ~400`).
- **NAMED REMAINING PROOF — C1-arm separation (the gate is not yet C1-VALID):**
  HP_HARD's boundary-substitution + combined-noise cases do NOT create a
  flat-affine error that a length-law could fix (flat-affine scores ~0.99-1.00 on
  them) — so they cannot, by a flat-vs-flat comparison, prove C1 discrimination.
  To make the gate C1-VALID, EITHER (a) build the C1 length-law arm (next cycle)
  and show arm-LAW > arm-flat on HP_HARD, OR (b) construct a stratum where
  flat-affine demonstrably MIS-PLACES the indel out of the run (e.g. a competing
  alignment that ties flat-affine but the length-law breaks toward truth — the
  run-bleeds-into-flank / adjacent-different-base-run cases). Until one exists,
  treat HP_HARD as exercising the scorer + the ablation harness, NOT as a proven
  C1 discriminator. The BAM-ization + two-arm harness (`run_flat_affine_arm` in
  the smoke, `cigar_records_to_bam` in the scorer) is the plug-in point.
- **STR stratum** is ~0.907 (mildly discriminating at this seed; an earlier seed
  showed 0.57). Harden STR (force the ED-tied whole-unit-deletion placement to
  differ from left-alignment) before relying on it for a C1 ablation.
- Tier-2 realistic transcriptome run (yeast saturation control + human SMN1/SMN2 +
  NIC/NNC panel) NOT yet run — that is the recall/FDR + tail-sizing tier (Sherlock).
- Other strata from the SPEC not yet generated (paralog, panel-failure/C5,
  coverage×Q, standing-variant/C6) — schema supports all; generators are the next
  build increment once the pbsim3 live run is confirmed.
- No C1/member code built (correctly — that is the next, gated cycle).

## RESUME
1. **Check the pbsim3 install:**
   `ssh sherlock 'ls /home/groups/larsms/users/kevinroy/anaconda3/envs/aligner_bench/bin/ | grep -iE "pbsim|badread"; tail -3 /tmp/pbsim3_install.log'`
   - If `pbsim` present → run a live mini round-trip (below).
   - If still solving / failed → install into a FRESH minimal env (faster solve):
     `ssh sherlock '/home/groups/larsms/users/kevinroy/anaconda3/bin/conda create -y -n pbsim3 -c bioconda -c conda-forge pbsim3 badread minimap2 samtools pysam'`
2. **Live pbsim3 round-trip** (find the packaged model first):
   `MODEL=$(ssh sherlock 'find /home/groups/larsms/users/kevinroy/anaconda3/envs/*/ -name "ERRHMM-ONT.model" | head -1')`
   then drive `scripts/benchmark/sim/pbsim3_wrapper.py:simulate_and_propagate`
   over a few `TranscriptModel`s, align the reads with minimap2 -ax splice, and
   run `scorer.score_bam` — confirm junction TP + indel concordance on REAL pbsim
   reads (the projection-truth path end-to-end). Use the `sherlock-sbatch` skill
   (owners partition, AVX-512 exclude); do NOT relay BAMs through the M1.
3. **Scale Tier-1** to `--reps 400` and audit every C1 cell ≥100 before any C1
   ablation (else the length-law nullifies to flat → FALSE REFUTE).
4. **Local smoke (regression gate), any time:**
   `PATH="/Users/kevinroy/miniconda3/bin:/opt/homebrew/bin:$PATH" \
    /Users/kevinroy/miniconda3/envs/pysam/bin/python \
    scripts/benchmark/smoke_roundtrip.py --out /tmp/bench_smoke` (exit 0 = green).
   (Local `pysam` env = pysam+numpy+minimap2+samtools; the rectify env is on Sherlock.)

## FILES (all NEW, benchmark-only paths — none touch shared hot files)
- `rectify/core/benchmark/__init__.py`
- `rectify/core/benchmark/truth_schema.py`   (component 2 — schema)
- `rectify/core/benchmark/scorer.py`          (component 3 — ambiguity-aware scorer)
- `scripts/benchmark/sim/transcript_model.py` (component 1 — truth backbone)
- `scripts/benchmark/sim/pbsim3_wrapper.py`   (component 1 — Tier-2 pbsim3 wrapper)
- `scripts/benchmark/sim/controlled.py`       (component 1 — Tier-1 generators)
- `scripts/benchmark/smoke_roundtrip.py`      (end-to-end GATE smoke)
- `dev/SIMULATION_BENCHMARK_SPEC.md`          (UPDATED: simulator decision recorded)
