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

## VERIFIED
- **GATE smoke is GREEN** (`scripts/benchmark/smoke_roundtrip.py`, minimap2 -ax
  splice on the Tier-1 corpus, 7230 reads):
  - (A) known junction round-trips truth->reads->minimap2->scorer as **TP=30**;
  - (B) a deliberately 1bp-shifted junction call (into the repeat) scores
    **TP=1 FP=0** (ambiguity-aware match, the GMAP-0.09 trap defended);
  - (B2) NEGATIVE control: a junction shifted 40bp (beyond the window) is correctly **FP**;
  - (C) indel **position-exact concordance = 0.920** (correct=2975, incorrect=258),
    false-indel-rate=0.0 over 3997 clean reads.
- Schema lossless TSV round-trip verified (junction normalization, canonicity,
  HP-cell accounting, variant/split round-trip).
- pbsim3 MAF→genome **projection** verified locally on a synthetic MAF (2bp deletion
  + 1 intron correctly propagated to genome coords) — no live pbsim needed for this.

## OPEN
- **pbsim3 live run NOT yet executed.** Install on Sherlock env `aligner_bench`
  was still solving (slow conda solve on the large base env) at handoff time —
  the wrapper is code-complete + projection-validated but a real pbsim3 round-trip
  is unverified. (badread also requested in the same install for cross-check.)
- **Tier-1 cell sizing:** at `--reps 120` the min indel-bearing cell is ~50 reads
  (clean reads are ~50% by the del-dominant K_DIST). The Sherlock scale-up must set
  `--reps` so every `(base_class, run_copies)` cell clears `min_count=100`
  (recommend `--reps ~400`).
- Tier-2 realistic transcriptome run (yeast saturation control + human SMN1/SMN2 +
  NIC/NNC panel) NOT yet run — that is the recall/FDR + tail-sizing tier (Sherlock).
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
