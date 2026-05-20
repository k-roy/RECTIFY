# Nanopore Aligner Recommendations - 2026 Update

**Date**: 2026-03-10  
**Updated**: 2026-05-19
**Context**: Five-aligner production panel active on Sherlock; deSALT restored via conda symlink (see deSALT section)

---

## Production Aligner Panel (Current)

Rectify uses a five-aligner consensus panel in production. deSALT has a
platform-specific binary issue on Sherlock (see
[deSALT platform notes](#desalt--platform-specific-binary-issue) below); the conda
binary on H2 works correctly.

| Aligner | Version | Notes |
|---------|---------|-------|
| **minimap2** | 2.28 | System PATH; splice-aware, general-purpose |
| **mapPacBio** | — | System PATH; PacBio/ONT long-read aligner |
| **gapmm2** | — | System PATH; gap-aware minimap2 wrapper |
| **uLTRA** | 0.1 | `/path/to/uLTRA`; annotation-guided, collinear chaining; requires sibling GTF + namfinder (see [ultra.md](aligners/ultra.md)) |
| **deSALT** | — | Conda binary on both H2 and Sherlock; Sherlock uses `~/bin/deSALT` symlink (see below) |

uLTRA requires a separate installation; pass the path via `rectify install-aligners`
or `--ultra-path`.

---

## uLTRA — GTF Annotation Requirements

uLTRA requires every `exon` feature in the GTF to carry a `gene_id` attribute.
GTFs produced by `gffread -T` from a yeast GFF fail this requirement: ~10% of
exon lines (tRNAs, ncRNAs, some multi-exon mRNAs) have `transcript_id` and
`gene_name` but no `gene_id`, causing a crash at `create_augmented_gene.py:323`.

**Symptom**: `database.db` is created (prep_splicing succeeds) but no
`.uLTRA.bam` is emitted.  RECTIFY logs `uLTRA failed: ... KeyError: 'gene_id'`
at ERROR level; the pipeline continues with the remaining aligners.

**Fix (built into RECTIFY since 2026-05-19)**: `run_ultra()` always runs
`_normalize_gtf_for_ultra()` for GFF-sourced annotations, which collects
existing exon intervals and injects `gene_id` from the parent transcript on
any line missing it.  No manual intervention required with current RECTIFY.

For full details, additional failure modes, and verification steps, see
[docs/aligners/ultra.md](aligners/ultra.md).

---

## deSALT — Sherlock Binary (RESOLVED 2026-03-10)

deSALT works correctly on both H2 and Sherlock using the conda-installed binary.

**H2 (conda binary) — working:**
`~/.conda/envs/rectify/bin/deSALT` (856 KB, built 2025-09-29). Confirmed producing
real alignments: 5065 records, 99.5% mapped, correct CIGAR strings with splice N-ops
on a wt_rep1 subsample.

**Sherlock (conda symlink) — working with rare residual crashes:**
`~/bin/deSALT` → `anaconda3/pkgs/desalt-1.5.6-he4a0461_5/bin/deSALT` (symlink created
2026-03-10). `run_desalt()` in `multi_aligner.py` resolves via `shutil.which('deSALT')`
before the vendored fallback; `run_array_others.sh` exports `$HOME/bin:$PATH`. Production
job 25072789 (2026-05-15) produced valid deSALT BAMs for the great majority of chunks;
merged per-sample BAMs are valid for all 12 set2 and 12 set3 samples.

Two chunks hit SIGSEGV with the conda binary during set2 re-alignment (May 2026):
- `wt_tfiiib_rep3/chunk_005` — 420,447 reads → exit -11
- `rna15_rep3/chunk_003` — 203,832 reads → exit -11 (original SAM also truncated from prior SIGTERM)

Both triggered the empty-BAM fallback; those chunks use 4-aligner consensus in the
correction step. The crash is non-monotonic and read-composition dependent (not purely a
read-count threshold). See `docs/desalt_crash_investigation_handoff.md` for investigation
notes and reproduction approach.

**Background (pre-2026-03-10):**
The vendored binary at `rectify/data/bin/linux_x86_64/deSALT` (v1.5.6) produced SIGSEGV
(exit -11 / shell 139) on Sherlock for certain read batches. Root cause: platform ABI
incompatibility; consistent with upstream deSALT#49 (SIGSEGV in `Loop-ProcessReads`).
The conda bioconda build (`he4a0461_5`) avoids this crash. The correction hang observed
in smoke test 25410215 was caused by the empty-BAM fallbacks from that era; it does not
occur with real deSALT BAMs.

**The vendored binary is still present** but is never used on Sherlock as long as
`~/bin/deSALT` is on PATH. For portability, any future Sherlock environment setup should
include: `ln -sf $(conda run -n rectify which deSALT) ~/bin/deSALT`

---

## Historical Installation Notes

### ❌ Removed

| Aligner | Issue | Recommendation |
|---------|-------|----------------|
| **GraphMap2** | Not updated since 2018, compilation issues | **Replaced** by current panel |

---

## Notes on Candidate Aligners (2025–2026)

The following are candidates to evaluate as additions to the production panel.
deSALT is no longer a stability concern on Sherlock; the rationale for evaluating
these is to improve consensus accuracy with diverse algorithmic approaches, not
to replace a broken aligner.

| Aligner | Status | Notes |
|---------|--------|-------|
| **Minisplice** | Released ~June 2025; minimap2 + deep-learning splice signals | High priority — DL splice signals complement the existing rule-based panel |
| **GLASS** | Published April 2025; graph-learning splice-aware alignment | High priority — graph approach closest in spirit to deSALT's De Bruijn strategy |
| **Winnowmap2** | Actively updated; optimized for repetitive/human genomics | Lower priority for yeast but worth benchmarking for multi-organism generality |
| **GraphMap2** | Last updated 2018, compilation issues | Replaced |

None have been integrated into the production panel. See `rectify/core/align/multi_aligner.py`
for the canonical list of supported aligners. See TODO for the evaluation work item.

---

## Adding Aligners to the Production Panel

To integrate a new aligner:

1. Add a wrapper function in `rectify/core/align/multi_aligner.py` following the existing
   `run_minimap2`, `run_mappacbio`, `run_gapmm2`, `run_ultra`, `run_desalt` patterns.
2. Register the aligner name in `SUPPORTED_ALIGNERS` and the consensus scoring logic
   in `rectify/core/commands/consensus_command.py`.
3. Run `rectify install-aligners --check` to verify the new binary is on PATH.
4. Validate against the bundled test dataset: `pytest tests/test_consensus_selection.py`

---

## Benchmarking Criteria

When comparing aligner panels, evaluate:

### Accuracy Metrics
- Junction detection sensitivity
- Canonical dinucleotide percentage
- Annotation match rate
- False positive rate

### Performance Metrics
- Alignment time per chunk
- Memory usage
- Parallelization efficiency
- Overall throughput

### Consensus Metrics
- Agreement rate between aligners
- Consensus confidence distribution
- Novel junction validation

---

## References

**Production Aligners**:
- minimap2: Li, H. (2018). Bioinformatics.
- uLTRA: Sahlin et al. (2021). Genome Biology.
- deSALT: Liu et al. (2019). Genome Biology.

**Candidate Aligners (not yet integrated)**:
- Minisplice: Released ~June 2025, enhances minimap2 with DL splice signals
- Winnowmap2: Actively updated, optimized for human genomics
- GLASS: Published April 2025, graph learning for splice-aware alignment

**Obsolete**:
- GraphMap2: Last updated 2018, no longer maintained

---

**Last Updated**: 2026-05-19
**Status**: Five-aligner production panel (minimap2, mapPacBio, gapmm2, uLTRA, deSALT). deSALT restored on Sherlock via `~/bin/deSALT` conda symlink; merged per-sample BAMs valid for all set2/set3 samples. 2/~200 set2 chunks hit residual SIGSEGV (empty-BAM fallback handles them; see `docs/desalt_crash_investigation_handoff.md`). uLTRA GTF normalization bug (KeyError: gene_id) fixed 2026-05-19.
