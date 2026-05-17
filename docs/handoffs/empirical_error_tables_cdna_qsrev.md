# Handoff: Empirical Error Tables for cDNA and QSrev Protocols

**Status:** Plumbing done (2026-05-17 PM). cDNA + QSrev calibration jobs
are staged on H2 but not yet submitted — see §9 below.

**Goal:** Generate `penalty_scores.tsv`, `str_penalty_scores.tsv`, and
`junction_overhang_table.tsv` for **cDNA (PCR/ONT cDNA)** and **QSrev
(QuantSeq REV short read)** protocols, bundle them into rectify, and wire
`--Scer` auto-resolution to pick the right table per protocol.

DRS tables now live at
`rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/` —
mirror that layout per protocol.

---

## 1. What already works (DRS reference)

**Reference run (DRS):**
- Script: `/u/project/guillom/kevinroy/common/scripts/nanopore/empirical_cigar_error_profiler.py` on Hoffman2 (1,845 lines)
- Calibration data: `/u/project/guillom/kevinroy/common/scripts/nanopore/error_profile_20260422/`
- Overhang table calibration: `rectify/core/splice/calibrate_junction_overhang.py` (function `calibrate()` at line 554), produces `junction_overhang_table.tsv`
- Production SLURM submit script: `run_profiler_20260422.sh` (same dir on H2)
- Docs: `PENALTY_TABLE.md` in same dir, mirrored at CLAUDE.md → "Empirical Penalty Tables" section
- Bundling commit (this handoff's prior step):
  - `rectify/data/__init__.py` — `get_bundled_junction_penalty_table()`, `get_bundled_str_penalty_table()`, `get_bundled_junction_overhang_table()`, plus auto-resolution in `resolve_reference_paths()`
  - `BUNDLED_GENOMES['saccharomyces_cerevisiae']` now has `junction_penalty_table`, `str_penalty_table`, `junction_overhang_table` entries

**Method (DRS) — the algorithm to replicate per protocol:**

1. Align reads with all 5 aligners (minimap2, mapPacBio, gapmm2, uLTRA, deSALT)
2. For each read present in ≥ N aligners (4 of 5 for production; 3 for strand-split where deSALT was unavailable):
   - Intersect exonic intervals across aligners (regions every aligner judges as non-intronic)
   - At each position where ALL aligners agree on the CIGAR op AND the reference base, score one error event by op type (D/X/I) × base × HP length
3. `--isolation-flank 10` excludes positions within 10 bp of any non-agreed position — prevents alignment-artifact inflation near poorly-mapped regions
4. Isotonic-smooth penalties so deletion penalty is monotone non-decreasing with HP length (or non-increasing for the rate)
5. Calibrate junction-overhang minimum thresholds from junctions that ≥4 aligners agree on, binned by intron size, isotonic-smoothed

The DRS calibration was on wt_by4742_rep1 (S. cerevisiae R10.4.1 Dorado). For
cDNA and QSrev, use the analogous lab-curated wild-type S. cerevisiae datasets
in `dev_runs/` (paths below).

---

## 2. Per-protocol differences that need handling

The DRS profiler reads a 5-aligner BAM panel from `dev_runs/<sample>_chunked_*/`.
For each new protocol you'll point at a different panel AND likely change which
read regions count toward the error pool.

| Aspect | DRS | cDNA | QSrev |
|---|---|---|---|
| **Aligner panel** | minimap2 + mapPacBio + gapmm2 + uLTRA + deSALT (5) | Same 5 (long-read cDNA) | bbmap + bwa (2, short-read) |
| **Calibration sample** | `wt_by4742_rep1` DRS Dorado | `wt_by4742_rep1` (or equivalent) ONT cDNA | `wt_by4742_rep1` (or equivalent) QSrev |
| **Read length** | ~700–2500 bp long-read | ~500–2500 bp long-read (cDNA template-switch) | ~75–150 bp short-read |
| **Poly-A region status** | trimmed pre-alignment (Step 0) | trimmed pre-alignment (cdna_trim) | NOT in read at all (oligo-dT primer hybridizes at A-tract) |
| **5'/3' region quirks** | template-switch adapter at 3', truncations at 5' | template-switch on both ends; PCR chimera ends; SSP/RT primer at one end | reverse-complement orientation; CSP internal priming; A-tract walkback at 3' (after RC, looks like 3' A-run) |
| **Strand orientation** | Native RNA, `-uf` minimap2 | cDNA both strands; aligners produce paired plus/minus | RC of mRNA; A-tract on the right after orientation |
| **Region to exclude from error pool** | last 50 bp of aligned 3' end (poly-A boundary noise) | last 50 bp of aligned 3' end + first 50 bp of aligned 5' end (template-switch noise) | last 50 bp of aligned 3' end (CSP A-tract walkback noise); ALSO first 30 bp at 5' (CSP adapter remnant if present) |

**Why region exclusion matters:** the error profiler intentionally avoids
positions where the protocol introduces systematic, non-Dorado-noise errors.
For DRS the pA boundary is the main offender. For cDNA, template-switch
oligo bleed at BOTH ends. For QSrev, the CSP A-tract bleed near the 3' end.

You will need to:

1. **Add `--protocol {drs,cdna,qsrev}` flag** to `empirical_cigar_error_profiler.py` that toggles these region masks. Default = `drs` to preserve current behavior.
2. **Add `--exclude-3prime-bp N` and `--exclude-5prime-bp N`** integer flags so the protocol mask is honest about what's excluded. Defaults driven by `--protocol`:
   - `drs`: `--exclude-3prime-bp 50 --exclude-5prime-bp 0`
   - `cdna`: `--exclude-3prime-bp 50 --exclude-5prime-bp 50`
   - `qsrev`: `--exclude-3prime-bp 50 --exclude-5prime-bp 30`
3. **Inside the per-read scoring loop**, compute `aligned_start_on_read`, `aligned_end_on_read` (where on the read body the aligned region starts/ends after soft-clip stripping), and skip CIGAR positions whose read-coordinate falls within `--exclude-5prime-bp` of the read start or `--exclude-3prime-bp` of the read end. The DRS profiler does NOT do this today — it counts every aligned position. Reusing it for cDNA/QSrev without this guard would inflate end-of-read error rates.

---

## 3. Calibration datasets — where to find them

These are the lab wild-type S. cerevisiae datasets the lab has been using.
Confirm paths with Kevin if you can't find them at the listed locations.

**cDNA (long-read ONT cDNA):**
- On Hoffman2: ask Kevin for the chunked dev run path (analogue of
  `dev_runs/wt_by4742_rep1_chunked_20260412/` but cDNA). If it doesn't exist,
  the smallest viable workflow is:
  1. `rectify split <fastq.gz> -n 16` on the cDNA wild-type FASTQ
  2. Run all 5 aligners per chunk via the generated SLURM scripts
  3. Pass the chunked output dir as `--run-dir`
- See `tests/test_validation_reads_cdna.py` for the bundled validation reads
  format — those 12 reads can be a sanity check for protocol detection.

**QSrev (QuantSeq REV short-read):**
- Aligner panel changes: bbmap + bwa instead of the 5-aligner long-read panel
- Use `rectify split --short-read --dT-primed-cDNA -n <N>` to chunk the input
  FASTQ, then run-per-chunk + collect aligner BAMs
- Calibration sample: lab QSrev wild-type. Kevin can point at the right
  `/oak/.../short_read_data/wt_*` path
- The profiler currently hard-codes the 5-aligner DRS panel via
  `_DEFAULT_ALIGNERS`. Make `--aligners` mandatory when `--protocol qsrev`,
  or set `_DEFAULT_ALIGNERS_QSREV = ['bbmap', 'bwa']` and gate it on protocol.

---

## 4. Concrete deliverables

For EACH of `cdna` and `qsrev`:

1. **Patch `empirical_cigar_error_profiler.py`** with the protocol flag and
   region-exclusion knobs described in §2. Push the patched script back to
   `common/scripts/nanopore/` on H2 AND add the same file to the rectify
   repo at `scripts/calibration/empirical_cigar_error_profiler.py` so the
   tool isn't trapped on H2 (the current arrangement is a known smell).
2. **Run the profiler** on the calibration sample. Production scale: full
   wild-type sample, 8-worker parallel, ~3 h on a Sherlock larsms node
   (DRS reference timing). Output goes to
   `common/scripts/nanopore/error_profile_<protocol>_<YYYYMMDD>/`.
3. **Run the overhang calibrator** via
   `python -m rectify.core.splice.calibrate_junction_overhang` on the same
   per-chunk corrected BAMs from the calibration sample. (For QSrev short
   reads, junctions are short and rare — confirm with Kevin whether an
   overhang table is meaningful here; we may bundle an empty/stub for QSrev
   and rely on the DRS table as fallback.)
4. **Bundle the three tables** into rectify:
   - Copy to `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/`
     with protocol-suffixed names:
     `penalty_scores_<protocol>.tsv`,
     `str_penalty_scores_<protocol>.tsv`,
     `junction_overhang_table_<protocol>.tsv`
   - Extend `BUNDLED_GENOMES['saccharomyces_cerevisiae']` in
     `rectify/data/__init__.py` with a nested `protocols` map:

     ```python
     'protocols': {
         'drs':   {'junction_penalty_table': 'penalty_tables/penalty_scores.tsv',
                   'str_penalty_table':      'penalty_tables/str_penalty_scores.tsv',
                   'junction_overhang_table':'penalty_tables/junction_overhang_table.tsv'},
         'cdna':  {...},
         'qsrev': {...},
     },
     ```

     Update `get_bundled_junction_penalty_table()` etc. to take an optional
     `protocol` kwarg (default `'drs'` for back-compat), and update
     `resolve_reference_paths()` to pick the protocol from `args.dT_primed_cDNA`
     (`'cdna'` if set on a long-read run, `'qsrev'` if `args.short_read` is
     also set, else `'drs'`).
5. **Add a per-protocol README** at
   `rectify/data/genomes/saccharomyces_cerevisiae/penalty_tables/README.md`
   summarising:
   - Which sample each table was calibrated on
   - Date + commit of the profiler used
   - Region-exclusion settings
   - Top-10 deletion-rate rows for sanity inspection
6. **Update `docs/CLAUDE.md`** "Empirical Penalty Tables" section with the
   new tables and the protocol-aware resolution.

---

## 5. Testing the new tables

For each protocol-table set, run a small smoke test:

```bash
# Build a tiny test BAM input (~500 reads from the calibration sample)
samtools view -bh -s 0.001 wt_<protocol>.bam > /tmp/wt_test.bam

# Run correct with auto-resolved tables
rectify correct /tmp/wt_test.bam --Scer --dT-primed-cDNA \
    -o /tmp/wt_test.tsv
# (or --short-read --dT-primed-cDNA for QSrev)

# Inspect the corrected TSV — Module 2H should fire with the new tables
# (junction-refinement events appear in correction_applied when --aligner-bams
#  is also set; for the smoke test, just verify no crash and reasonable
#  high-confidence rate vs. heuristic baseline)
```

Then verify:

1. **Bundled-data unit test** — `tests/test_data_bundling.py` (create if absent):
   verify `get_bundled_junction_penalty_table(organism, protocol='cdna')`
   resolves to the right file path under all three protocols.
2. **Validation smoke tests** — `tests/test_validation_reads_cdna.py` and
   `tests/test_quantseq_rev_walkback.py` should still pass with default
   `--Scer` (auto-resolution now picks the cDNA/qsrev table). No expected-
   value updates needed — the bundled validation reads were chosen to be
   robust to penalty-table changes within reasonable bounds.

---

## 6. Open questions to ask Kevin before you start

1. **Calibration dataset paths** for cDNA and QSrev wild-type runs. The DRS
   profiler used `wt_by4742_rep1`; what are the cDNA and QSrev analogues?
2. **Is short-read overhang calibration meaningful?** Short reads have few
   introns; the overhang table may be a no-op or harmful if undercalibrated.
   Decide: skip, bundle the DRS table as fallback, or generate a calibrated
   QSrev-specific table.
3. **AT/CG base-class split** — DRS uses AT/CG. cDNA may have different
   substitution-rate asymmetry due to PCR error. QSrev (Illumina-class
   chemistry) is well-known for substitution-dominated errors and may not
   need an HP penalty table at all. Confirm whether to derive HP penalties
   for QSrev or just emit a stub indicating "use unit-cost edit distance".

---

## 7. What NOT to do

- Do not regenerate the DRS penalty tables. They're calibrated, bundled,
  and tested. Touch only the cDNA + QSrev pipelines.
- Do not change `_hp_edit_distance` defaults or the heuristic-fallback
  penalty values in `rectify/core/splice/junction_refiner.py`. The bundled
  tables override these at runtime; the heuristics are the fallback path
  for non-S. cerevisiae users.
- Do not introduce a new `rectify calibrate-penalty` subcommand for this
  pass — keep the profiler as a standalone script under `scripts/calibration/`.
  Promoting it to a first-class subcommand is a separate PR after Kevin
  reviews the protocol-aware version.
- Do not run the profiler on the bundled validation reads (only 36 reads,
  far too few). Use the full lab wild-type samples.

---

## 8. Estimated effort

- Patch profiler + add protocol flag: ~2 h (the script is well-structured,
  the changes are localized to the per-read scoring loop)
- Run cDNA calibration: ~3 h wall on Sherlock larsms (8 workers)
- Run QSrev calibration: ~1 h wall (short reads are faster)
- Bundle + wire `--Scer` resolution: ~1 h
- Tests + README: ~1 h

Total: ~half a day of agent time, plus calibration job wall time.

---

## 9. Progress update — 2026-05-17 PM

**Plumbing landed (verified locally):**

- `scripts/calibration/empirical_cigar_error_profiler.py` — patched copy of the
  H2 v1 profiler with three new flags. `--protocol {drs,cdna,qsrev}` resolves
  defaults for `--aligners` and the two new read-coord exclusion knobs
  `--exclude-5prime-bp` / `--exclude-3prime-bp`. The DRS default is `(0, 0)` to
  preserve byte-identical reproducibility of the bundled DRS tables; pass
  `--exclude-3prime-bp 50` to regenerate DRS with a cDNA-style mask. Per-read
  exclusion happens via a new `get_exonic_ranges_inner()` that trims each
  aligner's exonic intervals to the read-coord interior before the consensus
  intersection. Unit-tested against 8 hand-computed CIGAR cases (M/=/X/I/D/N/S
  combinations), and end-to-end smoke-tested on the bundled
  `validation_reads_cdna.bam` (`--protocol cdna --max-reads 5`).

- `rectify/data/__init__.py` — `BUNDLED_GENOMES['saccharomyces_cerevisiae']`
  gained a nested `protocols` map. The three `get_bundled_*_table()` accessors
  now accept an optional `protocol` kwarg (default `None` keeps back-compat).
  Resolution falls back to the DRS flat keys when a protocol-specific TSV is
  absent — so the cdna/qsrev entries are safe to wire **before** the actual
  TSVs land. `resolve_reference_paths()` derives the protocol from
  `args.dT_primed_cDNA` and `args.short_read` via a new `_protocol_for_args()`
  helper.

- `tests/test_data_bundling.py` — 13 new tests covering back-compat, explicit
  protocol kwargs, `_protocol_for_args()` mapping, and end-to-end resolution.
  Full suite: 945 passed / 30 skipped / 4 deselected (was 934 / 28 / 4).

- `scripts/calibration/{run_profiler_qsrev_h2.sh,run_profiler_cdna_h2.sh}` —
  end-to-end SGE submit scripts for H2. Each one: env activation, BAM->FASTQ
  (cDNA only), `rectify split`, per-chunk `rectify align`, post-align move
  (see note below), then the patched profiler. Both staged at
  `/u/scratch/k/kevinroy/calibration_2026_05/scripts/` on H2 alongside the
  patched profiler. Bash-syntax-validated; an end-to-end layout test on H2
  with the bundled `validation_reads.fastq.gz` confirmed
  `split → align → move → profile` flows correctly.

  **Empirical layout note:** `rectify align -o $D` writes BAMs as
  `$D/<chunk_id>.<aligner>.bam`, not in per-aligner subdirs. The profiler's
  `find_chunk_bams` expects `$D/<aligner>/<chunk_id>*.<aligner>.bam`, so
  after each chunk completes the scripts `mv` the BAMs into
  `align/<aligner>/`. Don't drop this step.

- `scripts/calibration/bundle_protocol_tables.py` — pools per-sample
  `error_counts.tsv` files, re-derives penalties via the profiler's own
  `compute_rates` + `derive_penalties`, and writes
  `penalty_scores_<protocol>.tsv` + `str_penalty_scores_<protocol>.tsv` into
  the bundled `penalty_tables/` directory. Optional `--overhang PATH` bundles
  an overhang TSV; omit for QSrev so the resolver falls back to the DRS
  overhang.

**Decisions confirmed with Kevin during execution:**

| Question | Decision |
|---|---|
| cDNA calibration sample | `/u/scratch/k/kevinroy/ont_cdna_v4/wt_rep2/out/stage1_consensus.bam` — UMI-consensus reads (median ~1.9 kb). Resulting table is calibrated for **consensus-grade** cDNA inputs; raw (non-UMI) cDNA will exceed these penalties. |
| QSrev calibration sample | Han 2023 `wt_R1` + `wt_R2` FASTQs at `/u/project/guillom/shared/raw/SRA/external/han_2023_quantseq_rev/fastq/` |
| QSrev scope | Calibrated HP only. Overhang falls back to DRS (short-read junctions too sparse, per §6 Q2). |
| Compute site | Hoffman2 (H2 conda env `rectify`, SGE/UGE scheduler). |
| QSrev exclusion budget | **(5'=10, 3'=16)** — the handoff's (30, 50) was infeasible for 76-bp Han 2023 reads (would drop every read). (10, 16) keeps ~50 bp interior. |
| QSrev panel source | Re-align WT FASTQs on H2 — the Han 2023 H2 deposit only has post-RECTIFY `.cpfixed.consensus.bam`, no per-aligner outputs. |

**Outstanding gap not yet addressed:** the patched profiler lives at
`/u/scratch/k/kevinroy/calibration_2026_05/scripts/` on H2, not at the
handoff §4 deliverable location of `/u/project/guillom/kevinroy/common/scripts/nanopore/`.
The auto-mode classifier blocked the overwrite of the shared script during
this session — Kevin can either approve the overwrite explicitly or accept
that the patched copy lives in scratch until the calibration runs complete
and the bundled tables land.

**Open / pending (require H2 wall time):**

1. Submit `qsub /u/scratch/k/kevinroy/calibration_2026_05/scripts/run_profiler_qsrev_h2.sh`
   (~1-2 h wall: align ~30 min + profile ~30 min per sample × 2 samples,
   h_rt=4h).
2. Submit `qsub /u/scratch/k/kevinroy/calibration_2026_05/scripts/run_profiler_cdna_h2.sh`
   (~6-12 h wall depending on aligner parallelism: BAM->FASTQ ~10 min, split
   <5 min, 5-aligner panel ~4-6 h serial-per-chunk even with
   `--parallel-aligners`, profile ~2-3 h, overhang ~30 min. h_rt=24h for
   safety; convert to a `qsub -t 1-16` chunk-array if alignment time
   becomes the gating factor).
3. rsync results back to M1, then per protocol:
   ```bash
   python scripts/calibration/bundle_protocol_tables.py \
       --protocol qsrev --inputs runs/error_profile_qsrev_*/wt_R{1,2}
   python scripts/calibration/bundle_protocol_tables.py \
       --protocol cdna --inputs runs/error_profile_cdna_*/ \
       --overhang runs/overhang_table_cdna_*/junction_overhang_table.tsv
   ```
4. Smoke-test via `pytest tests/test_data_bundling.py
   tests/test_validation_reads_cdna.py tests/test_quantseq_rev_walkback.py`.
5. Write `penalty_tables/README.md` with per-protocol calibration provenance
   and update `docs/CLAUDE.md` "Empirical Penalty Tables" section.

**Known caveat to document in step 5's README:** The cDNA table is calibrated
from UMI-consensus reads, so RECTIFY users running on raw (non-UMI) ONT cDNA
will see penalties that are optimistic relative to their data. Flag this
prominently — anyone hitting it will get suspiciously-confident junction calls.
