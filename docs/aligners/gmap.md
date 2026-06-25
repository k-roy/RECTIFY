# GMAP — candidate aligner (oligomer-hash + splice-scoring DP)

**Status (2026-06-15): REJECTED (final).** GMAP's first "disqualification" was a
**configuration artifact (`--nofails`), not a real failure** (see "The `--nofails`
dead end" below), so the 2026-06-14 algorithmic-orthogonality panel re-opened it as
a **top orthogonal candidate** and re-evaluated it correctly (`-n 1`). The correct
run then **failed all three pre-registered bars** on the SG-NEx A549 chr5 full-panel
test (2026-06-14), and a most-charitable re-test (2026-06-15, Sherlock job 29696736:
K-filtered self-rescue + shared junction pool + junction-guided scoring OFF)
**confirmed the fail** — annotated_rate 71.0% (bar ≥95%), anchor median 5 bp
(bar ~13 bp), STOLE-correct 2655 ≫ ADDED-correct 692. GMAP wins displace better
minimap2/deSALT alignments and introduce spurious novel junctions; the K-filter
cannot lift the intrinsic anchor placement. The wrapper remains shipped and opt-in
(`--junction-aligners gmap`) for further study, but GMAP is **not** a default panel
member. See ALIGNER_RECOMMENDATIONS.md → "Candidate Aligners" for the full decision.

## Why it was interesting (algorithmic orthogonality)

GMAP is one of the few long-read RNA aligners whose core algorithm is
**genuinely orthogonal to minimap2's minimizer seed-chain-align** — and, unlike
Graphmap2, it is orthogonal in its **splice model too**, not just its seeding:

- **Index/seeding:** a genome **sampled-oligomer hash table** (not minimizers,
  not FM-index, not suffix array) → **diagonalization** to localize candidate
  exon regions → **oligomer chaining** of short k-mers.
- **Splice DP:** nucleotide-level **genomic "sandwich" dynamic programming with
  its own explicit splice-site scoring** (canonical donor/acceptor models +
  intron-length maximization). This is a distinct DP from minimap2's
  `ksw_exts2_sse` kernel (which minimap2, winnowmap2, minisplice, gapmm2, and
  Graphmap2 all share). So GMAP adds **splice-model diversity**, the scarcest
  axis in the panel.

Family map: minimizer-chain (minimap2) · MEM+annotation (uLTRA) · de-Bruijn-graph
2-pass (deSALT) · k-mer-affine (BBMap/mapPacBio) · **oligomer-hash + own splice-DP
(GMAP)**.

## ONT suitability (evidence)

- Krizanović et al. 2018 (*Bioinformatics* 34(5):748) ranked **GMAP #1 on ONT**
  for exon-spanning. Synthetic ONT: 87.1% hit-all-exons, 98.8% aligned. **Real
  ONT (R9 2D, obsolete chemistry): 68.8% exon-hit, only ~23–30% positioned within
  ±5 nt.** So "87%" is the synthetic ceiling; real-ONT exon-spanning ≈ 69%.
- Caveat: that benchmark predates minimap2 and used R9 2D. Re-benchmark on
  current RNA004/PCB114 chemistry before production weighting.
- RECTIFY does its own 3′ poly-A walkback, so GMAP's exon-level positional slop
  matters far less here than for a naive pipeline — what we want from GMAP is the
  orthogonal **junction/exon-structure** call, which the consensus + anchor gate
  then quality-screen.

## The `--nofails` dead end (do not repeat)

The 2026-05-25 benchmark ran GMAP with **`--nofails`**, which forces *every* read
to emit an alignment. Result: 100% "mapping" but **38.2% GT-AG (human) / 31.6%
(yeast)** — i.e. garbage junctions from forced alignment of unmappable reads.
GMAP was disqualified on those numbers. **This was the flag, not the aligner.**
Correct usage drops `--nofails` and uses `-n 1` (report the single best path):
unmappable reads are simply not emitted, and the junction quality of the reads
GMAP *does* place is what the consensus should see.

## Correct invocation

```bash
# one-time index build (writable -D dir; the reference dir may be read-only)
gmap_build -d GRCh38_chr5 -D <index_dir> GRCh38_chr5.fa

# align (NO --nofails; -n 1 = best path only; SAM output emits N-op introns)
gmap -d GRCh38_chr5 -D <index_dir> -f samse -n 1 -t <threads> reads.fastq > out.sam
```

- Output: SAM with **N-op introns** (e.g. `9M32N8M`); M-CIGAR (not =/X) — fine,
  RECTIFY's HP-ED handles M by base-comparison against the genome.
- Splice is intrinsic (GMAP is a cDNA→genome spliced aligner by design); no
  special preset needed.

## Install / license / maintenance

- bioconda: `conda install -c bioconda gmap` (latest **2025.07.31** as of writing).
- License: **Apache-2.0**. Actively maintained (multiple 2024–2025 releases).
- Lowest-risk of the three 2026-06-14 candidates.

## Integration into RECTIFY

A `run_gmap()` wrapper now ships in `rectify/core/align/multi_aligner.py` and is
exposed as an opt-in junction aligner via `rectify align --junction-aligners gmap`.
It requires a pre-built db: pass `--gmap-db <dir>/<name>` or place a
`gmap_db/<genome_stem>` directory adjacent to the genome (the wrapper does **not**
auto-build it). On any GMAP/samtools failure it emits an empty name-sorted BAM (like
deSALT) so the consensus proceeds. For a direct orthogonality *measurement* you can
still run GMAP by hand (above) and feed its BAM to
`rectify correct --aligner-bams gmap:<bam>` + `merge_corrected_tsvs` (both are
aligner-name-agnostic). Per the rejection above, GMAP is not enabled in the default
panel.

## Sources

- Wu & Watanabe 2005, *Bioinformatics* (GMAP). bioconda `gmap`.
- Krizanović et al. 2018, *Bioinformatics* 34(5):748 ([PMC6192213](https://pmc.ncbi.nlm.nih.gov/articles/PMC6192213/)).
- 2026-06-14 orthogonality panel: memory `project-aligner-orthogonality-panel`.
