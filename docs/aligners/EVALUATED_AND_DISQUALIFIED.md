# Evaluated aligners — disqualified / dead ends (do not re-investigate)

This registry records long-read aligners we have evaluated and **rejected**, with
the specific disqualifying reason and source, so the investigation is not repeated.
Primary survey: the 2026-06-14 algorithmic-orthogonality Opus panel (memory
`project-aligner-orthogonality-panel`) + the 2026-05-25 in-house benchmark
(ALIGNER_RECOMMENDATIONS.md).

**Two hard gates** every RNA-panel candidate must pass: (1) **splice-aware** —
emits N-op introns in linear SAM; (2) **noisy-ONT-suitable** — works on ONT
direct-RNA error rates, not just PacBio HiFi / accurate transcripts. Most rejects
fail gate (1).

## Disqualified — NOT splice-aware (DNA/SV mappers; no N-op introns)

| Tool | Family / what it is | Why rejected |
|---|---|---|
| **lra** | minimizer + concave-gap sparse-DP chaining | SV/DNA only; no splice mode, no N-ops. Novel chaining, wrong domain. |
| **NGMLR** | convex-gap banded SW | SV breakpoint mapper; no intron model. |
| **VACmap** (2024) | variant-aware non-linear chaining (DAG) | DNA structural-variation only; no splice mode. Its non-linear chaining would mis-model introns as rearrangements. |
| **BLEND** | fuzzy-seed (SimHash over strobemers), inside minimap2 | No splice preset; overlap/read-mapping focused. Most novel *seeding* in the survey but no RNA path; stale (Nov 2022). |
| **mapquik** | k-min-mer (minimizer-space) seeding | No splice; AND HiFi-only (collapses on noisy reads). Double-disqualified. |
| **strobealign** | randstrobe/syncmer seeding | Short-read DNA mapper; no spliced long-read tool exists. (The strobemer *idea* is promising but unpackaged for spliced ONT.) |
| **Sublong** (Subread family) | seed-and-vote 16-mer hash | Long-read DNA mapper; does NOT detect junctions (Subjunc, the splice variant, is short-read only). |

## Disqualified — graph-path output (not linear SAM with N-ops)

| Tool | Why rejected |
|---|---|
| **minigraph** | Sequence-to-graph; outputs GAF graph paths, no spliced/RNA mode, no N-op SAM. |
| **GraphAligner** | Sequence-to-graph; GAF/.gam output, requires pre-built GFA, no de-novo splice. (Revisit only if RECTIFY ever moves to a spliced-pangenome-graph representation.) |

## Disqualified — catastrophic ONT failure (splice-aware but unusable on long noisy reads)

| Tool | Why rejected (source) |
|---|---|
| **STAR / STARlong** | Uncompressed suffix array; splice-aware, but aligned only **5.5%** of ONT reads (Krizanović 2018); STARlong needs heavy hand-tuning, user reports ~0.46% uniquely-mapped on ONT DRS. In-house: "incompatible — ONT reads >1000 bp not supported." Suffix-array exact seeds are structurally hostile to ONT error density. |
| **HISAT2 / HISAT-3N** | Hierarchical graph FM-index (the most orthogonal *index* surveyed) — but "unable to align almost any read," **<7%** on long reads (Krizanović 2018); indels capped <20 bp; built for short reads. |

## Redundant with an existing / chosen panel member

| Tool | Redundant with | Note |
|---|---|---|
| **GSNAP** | GMAP | Same oligomer-hash lineage, short-read-tuned. If GMAP is in, GSNAP adds ~zero orthogonality. |
| **deBGA** | deSALT | deSALT's RdBG-index ancestor; deBGA itself is not splice-aware. |
| **mm2-fast** | minimap2 | SIMD/ISA reimplementation, **identical output** → zero consensus diversity (faster binary only). |
| **mm2-plus** | minimap2 | Parallel-chaining perf fork, whole-genome-alignment oriented; no diversity. |

## Low-priority / unverified on ONT

| Tool | Note |
|---|---|
| **Spaln / Spaln3** | Block-index + banded spliced DP; cDNA/protein focused, **no ONT/noisy-read mode**. Spaln3 (2024) is protein-query-focused. Would need empirical test before any inclusion; not pursued. |
| **BLAT** | Clean-EST-era spliced aligner; no ONT mode, no active long-noisy-read support. |
| **2passtools** | Scripted in the 2026-05-25 benchmark but produced no output; not investigated further. |

## Out of scope — different problem class (not basecalled-read aligners)

| Tool | Why out of scope |
|---|---|
| **RawHash / RawHash2, Sigmap, UNCALLED** | Raw-signal (squiggle-level) mappers — real-time, pre-basecall, **DNA-only, no splice model** (RNA is explicit future work in all three). Orthogonal *input* modality, but they don't produce the per-read genomic SAM a basecalled-read consensus consumes. |
| **NanoSplicer, Uncalled4** | Signal-level junction/analysis tools, not read-to-genome SAM aligners. (NanoSplicer additionally underperforms on direct-RNA vs cDNA.) |

## Not an aligner (adjacent ideas worth remembering)

| Tool | What it is | Disposition |
|---|---|---|
| **GLASS** (bioRxiv 2025.04) | Post-alignment BAM **filter** (graph-ML "Read-AS Map" removes falsely-spliced reads from minimap2 output) — no seeding of its own. Code not publicly released. | Not a panel member. Possibly interesting as a *post-consensus* false-splice filter (conceptually adjacent to RECTIFY's own anchor gate) if code is released. |
| **minisplice** | A 1D-CNN that scores **canonical GT/AG** sites genome-wide, fed to minimap2/miniprot as a splice prior. | Already an opt-in panel member (`--aligners minisplice_mm2`), but **note the motif-bias caveat**: it is a *learned GT-AG* prior → at odds with the lab's motif-agnostic unbiased-discovery rule. Keep as an *additional* member, never a replacement; weigh against non-canonical-junction recall. |

## The seeding-frontier negative result (important)

The novel-seeding frontier (strobemers, syncmers, k-min-mers, fuzzy/SimHash seeds)
yielded **zero** usable orthogonal ONT-RNA panel members: every novel-seeding tool
is DNA/SV-only (no splice), HiFi-only, or a minimap2 derivative. Syncmer seeding is
reachable only as a DNA-only minimap2/winnowmap2 fork (Shaw–Yu 2022), not as a
standalone spliced ONT aligner. **Do not re-scour the seeding literature for an
RNA panel member** — the orthogonal gains are in the classical (GMAP, Magic-BLAST)
and graph/ONT-specialized (Graphmap2) families, not in new seeds.

See ALIGNER_RECOMMENDATIONS.md → "Algorithmic Orthogonality Survey" for the
accepted candidates and the headline benchmark-gap finding.
