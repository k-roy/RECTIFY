# S. cerevisiae bundled penalty tables

Empirical per-base error penalties for `rectify correct`'s junction-refinement
module (Module 2H). Resolution is **protocol-aware** via
`rectify/data/__init__.py::resolve_reference_paths()` — passing `--Scer`
together with the right combination of `--dT-primed-cDNA` / `--short-read`
picks the matching table:

| Caller flags | Resolved tables | Underlying calibration |
|---|---|---|
| `--Scer` (no protocol flags) | `*.tsv` (DRS — back-compat default) | BY4742 wild-type direct RNA-seq, R10.4.1 / Dorado |
| `--Scer --dT-primed-cDNA` | `*_cdna.tsv` pooled (falls back to DRS if absent) | ONT PCR-cDNA, UMI-consensus, BY4742 wt_rep2 |
| `--Scer --ONT-cDNA` | per-UMI-bin `*_cdna_umi{1,2,3plus}.tsv` + pooled fallback, routed per-read by `XC` tag | ONT PCR-cDNA, UMI-consensus, BY4742 wt_rep2 |
| `--Scer --short-read --dT-primed-cDNA` | `*_qsrev.tsv` (overhang falls back to DRS) | Han 2023 QuantSeq REV WT (76 bp Illumina) |

The `--ONT-cDNA` path uses a `PenaltyTableSet` (Module 2H, built at startup)
rather than a single static table: each read's `XC:i:N` tag (UMI cluster size
written by the cDNA pipeline) selects the matching calibrated table at scoring
time.  Reads without an `XC` tag fall back to the pooled table.

The fallback to DRS for missing protocol-specific tables is intentional — it
means a partially-calibrated bundle still resolves to a working table set,
rather than crashing on `None`.

---

## Files

### DRS (default — calibrated 2026-04-22; verified clean of `=` SEQ encoding)

| File | Purpose |
|---|---|
| `penalty_scores.tsv`         | HP-length-aware D/X/I penalties, AT vs CG split |
| `str_penalty_scores.tsv`     | Dinucleotide / trinucleotide STR slippage penalties |
| `junction_overhang_table.tsv`| Minimum overhang threshold per intron-size bin |

**Calibration sample:** `wt_by4742_rep1` (BY4742 wild-type DRS, R10.4.1
flowcell, Dorado basecaller). 5-aligner concordance panel:
minimap2 + mapPacBio + deSALT + uLTRA + gapmm2.
**Read-coord exclusion:** none (5'=0, 3'=0).
**Calibrator script:** `scripts/calibration/empirical_cigar_error_profiler.py`
(commit at calibration time predates the protocol-aware refactor).

**`=` SEQ encoding check (2026-05-17):** the cDNA/QSrev calibration cycle
surfaced a profiler bug where `=` characters in BAM SEQ (SAM/BAM spec
"matches reference at this position" compression) were being miscounted as
substitutions, silently zeroing M counts. Verified the DRS dev_run BAMs at
`/oak/.../dev_runs/wt_by4742_rep1_drs_trim_20260417/` contain real
nucleotide bases (no `=`), so the bundled DRS table is unaffected and does
not need regeneration.

### cDNA — calibrated 2026-05-17/18 (UMI-bin-stratified)

| File | Status | UMI cluster size | Read count | HP events |
|---|---|---|---|---|
| `penalty_scores_cdna_umi1.tsv`    | ✓ bundled (125 rows) | 1 (singleton) | 399,644 | 368 M |
| `penalty_scores_cdna_umi2.tsv`    | ✓ bundled (117 rows) | 2 (doubleton) | 58,166 | 74 M |
| `penalty_scores_cdna_umi3plus.tsv`| ✓ bundled (106 rows) | ≥3 (consensus) | 25,605 | 35 M |
| `penalty_scores_cdna.tsv`         | ✓ bundled (125 rows, pooled) | all | 483,415 | 477 M |
| `str_penalty_scores_cdna.tsv`     | not bundled — STR profiling yielded 0 events; resolver falls back to DRS |
| `junction_overhang_table_cdna.tsv`| not bundled — falls back to DRS table |

#### UMI-bin routing (Phase C)

`rectify correct --ONT-cDNA` builds a `PenaltyTableSet` at Module 2H startup
that holds all four tables. For each read, the `XC:i:N` tag (UMI cluster size,
written by `rectify consensus`) selects the matching table:

| `XC` value | Table used |
|---|---|
| 1 | `penalty_scores_cdna_umi1.tsv` — calibrated on raw-accuracy singletons |
| 2 | `penalty_scores_cdna_umi2.tsv` — calibrated on doubleton consensus |
| ≥3 | `penalty_scores_cdna_umi3plus.tsv` — calibrated on high-confidence consensus |
| absent / untagged | `penalty_scores_cdna.tsv` — pooled fallback |

Tag note: `XC:i:N` = UMI cluster size (number of reads merged). `XK:i:N` is a
different tag (3′ trim length) — do not confuse the two.

#### Why stratify by UMI cluster size?

The pooled cDNA table averages over a ~75/15/10 singleton/doubleton/triplet mix.
Singletons have raw per-base error rates ~0.2–0.5 %, while ≥3-read consensus
reads are nearly error-free (<0.05 %). A single pooled table is biased optimistic
for singletons and conservative for high-consensus reads. Stratification lets
RECTIFY apply appropriate junction-scoring confidence per read.

**Calibration sample:** `stage1_consensus.bam` from
`/u/scratch/k/kevinroy/ont_cdna_v4/wt_rep2/out/` — **UMI-consensus** ONT
PCR-cDNA reads (median ~1.9 kb). 5-aligner concordance panel:
minimap2 + mapPacBio + gapmm2 + uLTRA + deSALT (same as DRS).
**Read-coord exclusion:** (5'=50 bp, 3'=50 bp) to mask template-switch
adapter bleed at both ends.
**Compute:** Sherlock larsms partition (JID 25334819), 8 chunks, all 5 aligners.

> **Important caveat:** the per-bin tables are calibrated on UMI-consensus
> reads from a single replicate (`wt_rep2`). Running RECTIFY on raw
> (non-UMI-collapsed) PCR-cDNA will apply the `umi1` (singleton) table, which
> is calibrated on consensus-quality singletons — still more accurate than the
> pooled table for single-pass reads, but not specifically calibrated for
> fully-raw ONT cDNA. If you don't have a UMI consensus step, prefer the DRS
> table or pass `--junction-penalty-table` explicitly.

**Convergence summary (8-chunk run, median |Δ%| across consecutive chunk pairs):**

| Bin | Median \|Δ%\| | Verdict |
|---|---|---|
| umi1 | 3.12% | 8 chunks sufficient |
| umi2 | 3.03% | 8 chunks sufficient |
| umi3+ | 1.98% | 8 chunks sufficient |

### QSrev — calibrated 2026-05-17

| File | Status |
|---|---|
| `penalty_scores_qsrev.tsv`     | ✓ bundled (579 M HP events pooled across wt_R1 + wt_R2) |
| `str_penalty_scores_qsrev.tsv` | not bundled — STR profiling yielded 0 events; resolver falls back to DRS |
| `junction_overhang_table_qsrev.tsv` | not bundled — falls back to DRS table |

**Calibration sample:** Han et al. 2023 QuantSeq REV WT (`wt_R1` + `wt_R2`,
76 bp single-end Illumina). 2-aligner concordance panel: bbmap + bwa.
**Read-coord exclusion:** (5'=10 bp, 3'=16 bp) — keeps ~50 bp interior; the
handoff's original (30, 50) would have dropped every read at 76 bp.
**Compute:** H2 SGE (JID 13421248), ~50 min total wall (split + bbmap + bwa
+ profile per replicate). Bundled via `bundle_protocol_tables.py --protocol qsrev`.

**Notes:**
- Illumina-class chemistry is substitution-dominated. Empirical D and I rates
  are ~3 × 10⁻⁴ — far below sub rates — so penalty_score for D and I clamps
  at 10.0 (the cap) for most HP lengths. The table is therefore mostly
  *substitution-driven*; using it is roughly equivalent to a unit-cost-edit
  model that down-weights insertions and deletions strongly. See handoff
  §6 Q3 — Kevin chose "calibrated HP, skip overhang" precisely because the
  short-read chemistry doesn't need overhang nuance.
- STR profiling emitted 0 events because the (10, 16) read-coord exclusion
  + 5+ copy threshold left no qualifying contexts in 76 bp reads. The
  bundled DRS STR table is used as fallback — fine for QSrev since STR
  slippage is rare in Illumina reads.
- The bundled `penalty_scores_qsrev.tsv` has `ins(HP=1) = 0` for both AT and
  CG (no insertion events at HP=1 anywhere in 6.99M × 2 reads). This is
  partly real (Illumina rarely inserts) and partly a calibration artifact:
  the (10, 16) read-coord mask excises both ends of the 76 bp reads, and
  the residual ~50 bp middle has too few HP=1 contexts to sample a rare
  event. The bundled DRS HP=1 insertion rate is used implicitly via the
  `--default-ins 1.25` baseline.

**Why no QSrev overhang table:** short-read junctions are too rare in
QuantSeq REV to calibrate a meaningful intron-size-binned threshold.

---

## Regeneration

DRS — uses the legacy invocation (preserved for byte-identical reproducibility):

```bash
python scripts/calibration/empirical_cigar_error_profiler.py \
    --run-dir <DRS dev_run> \
    --reference .../S288C_reference_sequence_R64-5-1_20240529.fsa.gz \
    --output-dir ./error_profile_drs_<date> \
    --isolation-flank 10 --union --str-repeat --workers 8
# DRS defaults to --protocol drs which keeps --exclude-{5,3}prime-bp at 0.
```

cDNA / QSrev — see the SGE/SLURM submit scripts in `scripts/calibration/`:
- `run_profiler_cdna_sherlock.sh` (Sherlock, SLURM, ~1-2 h wall for 8 chunks)
- `run_profiler_cdna_h2.sh`       (H2, SGE/UGE, ~6-12 h wall)
- `run_profiler_qsrev_h2.sh`      (H2, SGE/UGE, ~1-2 h wall)

The profiler writes per-bin error count TSVs when `--umi-tag XC` and
`--umi-bins 1,2,3+` are passed (the Sherlock script already sets these).
After the per-sample profiler outputs land, bundle **each bin separately**:

```bash
# Per-bin tables (run once per bin):
for BIN in umi1 umi2 umi3plus; do
    python scripts/calibration/bundle_protocol_tables.py \
        --protocol cdna \
        --umi-bin "$BIN" \
        --inputs path/to/error_profile_cdna_*/
done

# Pooled table (all reads together):
python scripts/calibration/bundle_protocol_tables.py \
    --protocol cdna \
    --inputs path/to/error_profile_cdna_*/
```

For QSrev (no UMI binning):
```bash
python scripts/calibration/bundle_protocol_tables.py \
    --protocol qsrev \
    --inputs path/to/error_profile_qsrev_*/  \
    [--overhang path/to/junction_overhang_table.tsv]
```

Each `bundle_protocol_tables.py` invocation writes the protocol-suffixed TSV
directly into this directory. The `BUNDLED_GENOMES[...]['protocols']` map
already points at the expected filenames, so resolution starts working
immediately after the files are written.
