# Nanopore Aligner Recommendations - 2026 Update

**Date**: 2026-03-10  
**Updated**: 2026-04-27  
**Context**: Updating aligner panel based on latest developments

---

## Production Aligner Panel (Current)

Rectify uses a five-aligner consensus panel in production. All five are required for
`rectify consensus` aligner selection:

| Aligner | Version | Notes |
|---------|---------|-------|
| **minimap2** | 2.28 | System PATH; splice-aware, general-purpose |
| **mapPacBio** | — | System PATH; PacBio/ONT long-read aligner |
| **gapmm2** | — | System PATH; gap-aware minimap2 wrapper |
| **uLTRA** | 0.1 | `/path/to/uLTRA`; annotation-guided, collinear chaining |
| **deSALT** | 1.5.6 | Vendored binary at `rectify/data/bin/linux_x86_64/deSALT`; De Bruijn graph |

deSALT is vendored as a binary within the rectify package. uLTRA requires a separate
installation; pass the path via `rectify install-aligners` or `--ultra-path`.

---

## Historical Installation Notes

### ❌ Removed

| Aligner | Issue | Recommendation |
|---------|-------|----------------|
| **GraphMap2** | Not updated since 2018, compilation issues | **Replaced** by current panel |

---

## Notes on Candidate Aligners (2025–2026)

The following were evaluated as potential additions to the panel:

| Aligner | Status | Outcome |
|---------|--------|---------|
| **Minisplice** | Released ~June 2025; minimap2 + deep-learning splice signals | Under evaluation |
| **GLASS** | Published April 2025; graph-learning splice-aware alignment | Under evaluation |
| **Winnowmap2** | Actively updated; optimized for human genomics | Low priority for yeast |
| **GraphMap2** | Last updated 2018, compilation issues | Replaced |

These have not been integrated into the production panel. See `rectify/core/align/multi_aligner.py`
for the canonical list of supported aligners.

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

**Last Updated**: 2026-04-27  
**Status**: Five-aligner production panel (minimap2, mapPacBio, gapmm2, uLTRA, deSALT) in use.
