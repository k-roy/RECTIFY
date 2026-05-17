# S. cerevisiae bundled penalty tables

Empirical per-base error penalties for `rectify correct`'s junction-refinement
module (Module 2H). Resolution is **protocol-aware** via
`rectify/data/__init__.py::resolve_reference_paths()` — passing `--Scer`
together with the right combination of `--dT-primed-cDNA` / `--short-read`
picks the matching table:

| Caller flags | Resolved tables | Underlying calibration |
|---|---|---|
| `--Scer` (no protocol flags) | `*.tsv` (DRS — back-compat default) | BY4742 wild-type direct RNA-seq, R10.4.1 / Dorado |
| `--Scer --dT-primed-cDNA` | `*_cdna.tsv` (falls back to DRS if absent) | ONT PCR-cDNA, UMI-consensus, BY4742 wt_rep2 |
| `--Scer --short-read --dT-primed-cDNA` | `*_qsrev.tsv` (overhang falls back to DRS) | Han 2023 QuantSeq REV WT (76 bp Illumina) |

The fallback to DRS for missing protocol-specific tables is intentional — it
means a partially-calibrated bundle still resolves to a working table set,
rather than crashing on `None`.

---

## Files

### DRS (default — calibrated 2026-04-22)

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

### cDNA — calibration in flight as of 2026-05-17

Output filenames (resolved automatically by `--Scer --dT-primed-cDNA`):

| File | Status |
|---|---|
| `penalty_scores_cdna.tsv`         | PENDING — calibration running on H2 (JID 13417082) |
| `str_penalty_scores_cdna.tsv`     | PENDING |
| `junction_overhang_table_cdna.tsv`| PENDING |

**Calibration sample:** `stage1_consensus.bam` from
`/u/scratch/k/kevinroy/ont_cdna_v4/wt_rep2/out/` — **UMI-consensus** ONT
PCR-cDNA reads (median ~1.9 kb). 5-aligner concordance panel: same as DRS.
**Read-coord exclusion:** (5'=50 bp, 3'=50 bp) to mask template-switch
adapter bleed at both ends.

> **Important caveat:** because the calibration is on UMI-consensus reads,
> the cDNA penalties are tuned for **consensus-grade** inputs. Running
> RECTIFY on raw (non-UMI-collapsed) PCR-cDNA against this table will yield
> optimistic penalties — junction calls will look more confident than they
> should. If you don't have a UMI consensus step in your pipeline, prefer
> the DRS table or pass an explicit custom `--junction-penalty-table`.

### QSrev — calibration in flight as of 2026-05-17

| File | Status |
|---|---|
| `penalty_scores_qsrev.tsv`     | PENDING — calibration running on H2 (JID 13417081) |
| `str_penalty_scores_qsrev.tsv` | PENDING |
| `junction_overhang_table_qsrev.tsv` | **NOT BUNDLED** — falls back to DRS table (see below) |

**Calibration sample:** Han et al. 2023 QuantSeq REV WT (`wt_R1` + `wt_R2`,
76 bp single-end). 2-aligner concordance panel: bbmap + bwa.
**Read-coord exclusion:** (5'=10 bp, 3'=16 bp) — keeps ~50 bp interior; the
handoff's original (30, 50) would have dropped every read at 76 bp.

**Why no QSrev overhang table:** short-read junctions are too rare in
QuantSeq REV to calibrate a meaningful intron-size-binned threshold. The
QSrev resolver falls back to the DRS overhang via the protocols-map
mechanism — the DRS minimums are conservative enough to not produce
false-positive splice rescues in short-read context.

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

cDNA / QSrev — see the SGE submit scripts in `scripts/calibration/`:
- `run_profiler_cdna_h2.sh`  (H2, SGE/UGE, ~6-12 h wall)
- `run_profiler_qsrev_h2.sh` (H2, SGE/UGE, ~1-2 h wall)

After the per-sample profiler outputs land:

```bash
python scripts/calibration/bundle_protocol_tables.py \
    --protocol {cdna,qsrev} \
    --inputs path/to/error_profile_*/  \
    [--overhang path/to/junction_overhang_table.tsv]
```

This writes the protocol-suffixed TSVs directly into this directory; the
`BUNDLED_GENOMES[...]['protocols']` map already points at the expected
filenames, so resolution starts working immediately.
