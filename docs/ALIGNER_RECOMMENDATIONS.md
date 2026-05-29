# Nanopore Aligner Recommendations

**Date**: 2026-03-10  
**Updated**: 2026-05-25

---

## Production Aligner Panel

Rectify uses a five-aligner consensus panel for yeast. For human/large-genome DRS, gapmm2 is dropped (see [gapmm2.md](aligners/gapmm2.md) — its per-read terminal-refinement loop costs ~28–50 h single-threaded on human scale).

| Aligner | Notes | Doc |
|---------|-------|-----|
| **minimap2** | Splice-aware, general-purpose. Key pitfall: `-G` defaults to 5000 (yeast); pass `--max-intron 500000` for human data. | [minimap2.md](aligners/minimap2.md) |
| **mapPacBio** | BBMap long-read mode. Needs `intronlen=10 maxindel=200000` for mammalian introns. `pt:i` embedded in QNAME. | [mapPacBio.md](aligners/mapPacBio.md) |
| **gapmm2** | Gap-aware minimap2 wrapper. **Yeast only in practice** — see gapmm2.md for timing, duplicate-UUID fix, PAF→BAM sequence injection, and minimum-overhang gate. | [gapmm2.md](aligners/gapmm2.md) |
| **uLTRA** | Annotation-guided collinear chaining; catches small exons seed aligners miss. Requires sibling GTF + namfinder; GTF normalization (gene_id) auto-applied since 2026-05-19. | [ultra.md](aligners/ultra.md) |
| **deSALT** | De Bruijn graph-based, high sensitivity. Linux only; on Sherlock use `~/bin/deSALT` conda symlink. Quirks: duplicate-alignment dedup, SIGSEGV fallback, must run outside thread pool. | [deSALT.md](aligners/deSALT.md) |

To add an aligner, see [Adding Aligners](#adding-aligners-to-the-production-panel) below.

> **Secondary / supplementary / duplicate-primary handling** is documented once,
> canonically, in [minimap2.md → Duplicate primary alignments](aligners/minimap2.md#-duplicate-primary-alignments--2-double-counted-3-ends-external-bam-hazard)
> (cross-aligner table + the by4742 2× external-BAM incident). Each aligner doc
> carries a short "Primary-alignment & duplicate handling" note pointing there.
> Summary: `rectify correct` filters `is_secondary`/`is_supplementary` for all
> aligners, but does **not** dedup identical **primary** records or honor `0x400`
> — so a doubled external BAM 2×-inflates 3′-end counts. Re-align through rectify
> (clean `--secondary=no` + DRS pre-trim) rather than trusting external BAMs.

---

## Benchmark: Human chr5 10k DRS Reads (SG-NEx A549, 2026-05-25)

**Setup:** 10k trimmed ONT DRS reads from SG-NEx A549 aligned to GRCh38 chr5.
Reference: `GRCh38_chr5.fa`. Reads: `a549_chr5.primary.random10k.seed20260525_trimmed.fastq.gz`.
Cluster: Sherlock `owners` partition. gapmm2 on sh02-03n47 Intel; others on AMD Milan sh04 nodes.

> **Strand-correction note**: an initial run of `compute_junction_stats.py` omitted strand handling —
> `get_splice_dinucs` fetched reference dinucleotides in reference orientation without accounting for
> `aln.is_reverse`. For `-` strand reads the left N-op boundary is the acceptor and the right is the
> donor; both need `rev_comp`. The ~50/50 GT-AG/non-canonical split in the first run was entirely this
> bug. Numbers below are from the corrected script. gapmm2's numbers (from `gapmm2_paf_junction_stats.py`,
> which uses PAF strand column + cs:Z: dinucleotides) were unaffected by the bug.

### Timing (wall seconds, single node)

| Aligner | Wall (s) | Notes |
|---------|----------|-------|
| minimap2 | 9 | baseline |
| winnowmap2 | 9 | |
| minisplice_mm2 | 16 | minimap2 + DL splice scores |
| gmap | 169 | |
| gapmm2 | **4,644** | sh02 Intel (0.46 s/read); AMD Milan ≈ 0.19–0.21 s/read |
| STARlong | incompatible | ONT reads >1000 bp not supported |

### Junction quality (strand-corrected)

| Aligner | Mapped | % Mapped | % GT-AG | % GC-AG | % Non-canonical | Notes |
|---------|--------|----------|---------|---------|----------------|-------|
| minimap2 | 8,937/10,512 | 85.0% | 91.6% | 3.0% | 4.7% | baseline |
| minisplice_mm2 | 8,967/10,512 | 85.3% | 96.0% | 1.7% | 2.1% | DL splice signals clearly improve over minimap2 |
| gapmm2 | 9,915/10,000 | 99.2%† | 97.2% | 2.4% | 0.3% | forces GT-AG via edlib; timing impractical |
| winnowmap2 | 9,072/10,526 | 86.2% | **99.0%** | 0.6% | **0.4%** | highest GT-AG of all; same speed as minimap2 |
| gmap | 8,727/8,727 | 100.0%‡ | 38.2% | 1.3% | 59.6% | `--nofails` inflates mapping; disqualified |
| STARlong | — | — | — | — | — | incompatible with ONT reads >1000 bp |

† gapmm2 deduplicates the FASTQ first (10,000 unique UUIDs after stripping DRS headers + empty-seq placeholders); PAF-based stats were already strand-correct via `cs:Z:` dinucleotides.  
‡ GMAP `--nofails` forces every read to produce an alignment; disqualified.

**Winnowmap2** (99.0% GT-AG, 0.4% non-canonical) is the strongest aligner — matches/exceeds gapmm2's GT-AG at 1/516th the wall time. Weighted minimizers suppress false positives in chr5's repetitive regions (SMN1/SMN2 locus). **Recommended addition for human panel.**

**minisplice_mm2** (96.0% GT-AG) clearly outperforms minimap2 baseline (91.6%). 1.8× slower than minimap2 (16 s vs 9 s on 10k reads) — negligible at scale. **Recommended addition.**

**gapmm2**: timing impractical at scale (~28–50 h single-threaded for a full sample; see [gapmm2.md](aligners/gapmm2.md)). Its 97.2% GT-AG is now matched/exceeded by winnowmap2. Dropped from human-DRS panel.

**Decision (2026-05-25):** gapmm2 dropped from human-DRS panel. Remains in yeast panel unchanged.

---

## Benchmark: Yeast 10k DRS Reads (2026-05-25, strand-corrected)

Same benchmark setup applied to a yeast 10k DRS read subset aligned to S288C R64-5-1.
gapmm2 was not re-benchmarked here (it is already in the yeast production panel).

| Aligner | Mapped | % Mapped | % GT-AG | % Non-canonical | Notes |
|---------|--------|----------|---------|----------------|-------|
| minimap2 | 8,726/14,287 | 61.1% | 71.1% | 22.1% | baseline |
| minisplice_mm2 | 8,726/14,288 | 61.1% | 76.6% | 18.7% | modest improvement over minimap2 on yeast |
| winnowmap2 | 8,620/13,769 | 62.6% | **90.8%** | 9.2% | strong improvement; weighted minimizers help on yeast repeats |
| gmap | 9,735/9,735 | 100.0%† | 31.6% | 67.2% | `--nofails`; disqualified |

† GMAP `--nofails` forces alignment of every read; disqualified.

**Winnowmap2** shows the strongest GT-AG on yeast (90.8% vs 71.1% minimap2 baseline) — worthwhile addition to the yeast panel as well. **minisplice_mm2** gives a modest boost over minimap2 at negligible speed cost.

---

## Candidate Aligners

| Aligner | Status | Notes |
|---------|--------|-------|
| **Winnowmap2** | **Integrate** | Highest GT-AG on both organisms; same speed as minimap2. Implemented in RECTIFY as `--aligners winnowmap2` (opt-in). Full-genome GRCh38 validation recommended before adding to `"all"`. See [winnowmap2.md](aligners/winnowmap2.md). |
| **Minisplice** (minisplice_mm2) | **Integrate** | minimap2 + DL splice signals; clear improvement over minimap2 baseline on DRS at ~1.8× cost. Implemented in RECTIFY as `--aligners minisplice_mm2` (opt-in; requires `--minisplice-model` or `--minisplice-scores`). See [minisplice.md](aligners/minisplice.md). |
| **GLASS** | **Not available** | Post-alignment BAM filter, NOT a standalone aligner (bioRxiv 2025.04.07.647681). Code not publicly released; cannot be benchmarked or installed. Revisit when code is released. |
| **2passtools** | Untested | Scripted in benchmark but produced no output; not investigated. |

---

## Adding Aligners to the Production Panel

1. Add a wrapper function in `rectify/core/align/multi_aligner.py` following the existing `run_minimap2`, `run_mappacbio`, `run_gapmm2`, `run_ultra`, `run_desalt` patterns.
2. Register the aligner name in `SUPPORTED_ALIGNERS` and the consensus scoring logic in `rectify/core/commands/consensus_command.py`.
3. Run `rectify install-aligners --check` to verify the new binary is on PATH.
4. Validate against the bundled test dataset: `pytest tests/test_consensus_selection.py`

---

## References

- minimap2: Li, H. (2018). Bioinformatics.
- uLTRA: Sahlin et al. (2021). Genome Biology.
- deSALT: Liu et al. (2019). Genome Biology.
- Minisplice: Released ~June 2025; minimap2 + deep-learning splice signals.
- Winnowmap2: Jain et al. (2022). Nature Methods; actively updated for human genomics.
- GLASS: Preprint 2025.04.07.647681; post-alignment filter, code not released.

---

**Last Updated**: 2026-05-25  
**Status**: Five-aligner yeast panel (minimap2, mapPacBio, gapmm2, uLTRA, deSALT). Human/large-genome: 4-aligner (gapmm2 dropped). Winnowmap2 + minisplice_mm2 benchmarked and recommended for integration. Chrom-naming fix (`chr5`→`chrV`) required for non-yeast runs (commit `99ff250`).
