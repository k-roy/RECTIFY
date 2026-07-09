# COMPASS short-read re-processing of SMA data — motif-agnostic junction validation (2026-07-09)

**PI insight:** Snaptron/recount3 (STAR-built, motif-filtered) is a CIRCULAR oracle for non-canonical junctions —
it flattens exactly what RECTIFY targets. The non-circular test is COMPASS short-read mode (per-read best-alignment
across BBMap/STAR×2/HISAT2×2/Magic-BLAST/GSNAP, NO canonical-motif post-filter) on REAL SMA short-read FASTQ, then
compare its motif-agnostic non-canonical junctions to the long-read native re-placer's calls (motif-agnostic vs
motif-agnostic = fair). Snaptron stays as a one-way CONFIRM-only oracle (positive match = real).

## Inputs (from dev/SMA_SHORTREAD_SURVEY.md)
- GSE108094 / SRP126773 (Rizzo 2019 Brain): human iPSC-MN, SMA-I vs control, 8 samples, HiSeq PE100, open FASTQ.
  ★ BEST match to the Sumner iPSC-MN panel. PULL FIRST.
- SRP334251 / PRJNA758038 (Ottesen 2023 NAR): SMA-I fibroblast (SMN1-null), 15 samples, NovaSeq deep. Depth for
  high-sensitivity SMN-locus / shared-gene confirmation (cell-type mismatch caveat).

## Plan
1. Pull GSE108094 FASTQ (SRA prefetch/fasterq-dump on Sherlock scratch).
2. Run COMPASS short-read mode (human GRCh38) — leverage the partially-built human COMPASS infra (compass_a549 on
   Sherlock; human scripts at oak .../collaborations/COMPASS/main/scripts/; smoke test had failed rc=1 -> revive).
   OR rectify --short-read mode if that's the current path. Family-gated, ambiguity-normalized, NO motif filter.
3. Extract COMPASS short-read junctions (canonical + non-canonical) per sample.
4. COMPARE to the long-read native re-placer's non-canonical calls (the sweep-B leads): for each re-placer novel
   non-canonical junction, is there a COMPASS short-read split read at the SAME coordinate (motif-agnostic)?
   - COMPASS-supported -> REAL novel non-canonical (the tool found genuine biology STAR/Snaptron missed) = the WIN.
   - COMPASS-unsupported AND high short-read coverage at the locus -> genuine fabrication/drift (now non-circular).
   - Low coverage -> inconclusive.
5. Calibrated controls (COMPASS_HANDOFF): annotated junctions in expressed loci (must validate high); shuffled null.

## Caveats
- COMPASS reduces but may not fully eliminate short-read placement bias; Magic-BLAST/BBMap are the least motif-biased.
- Short reads have less power at non-canonical junctions (short anchors) -> a real dRNA-only junction may still be
  short-read-invisible. The SYNTHETIC SPIKE-IN remains the ultimate unbiased ground truth (task #20).
- iPSC-MN (GSE108094) is a different cohort than Sumner -> shared/SMN-locus junctions comparable; sample-specific not.
