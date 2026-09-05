<p align="center">
  <img src="docs/branding/rectify_banner_light.png#gh-light-mode-only" alt="RECTIFY — Nanopore RNA Alignment Correction" width="700">
  <img src="docs/branding/rectify_banner_dark.png#gh-dark-mode-only"  alt="RECTIFY — Nanopore RNA Alignment Correction" width="700">
  <br>
  <em><b>R</b>NA 5' and 3' <b>E</b>nd <b>C</b>orrection <b>T</b>ool with <b>I</b>ntron re<b>F</b>inement and ambiguit<b>Y</b> resolution</em>
</p>
<p align="center">
  <a href="https://pypi.org/project/rectify-rna/"><img src="https://img.shields.io/pypi/v/rectify-rna?color=blue&label=PyPI" alt="PyPI"></a>
  <a href="https://rectify-rna.readthedocs.io"><img src="https://img.shields.io/readthedocs/rectify-rna?label=docs" alt="Docs"></a>
  <a href="https://github.com/k-roy/RECTIFY/actions"><img src="https://img.shields.io/github/actions/workflow/status/k-roy/RECTIFY/tests.yml?branch=master&label=tests" alt="Tests"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.9%2B-blue.svg" alt="Python 3.9+"></a>
</p>

Off-the-shelf aligners often misplace the ends of poly(A)-RNA reads: 5' splice-junction overhangs that are soft-clipped, junctions forced to annotated sites when an alternative is a better match, and poly(A) tails that align over genomic A-tracts so the apparent 3' end overshoots the true cleavage site — sometimes by thousands of bases with an artifactual intron added. **RECTIFY corrects each of these in one pass**: it runs multiple aligners in parallel, corrects each independently using **chemistry-specific empirical error models** (indel penalties calibrated by homopolymer length and base class, with bundled tables for *S. cerevisiae* and *H. sapiens*), then selects the single best-corrected alignment per read. The bases at the read ends then mean what biology says they should — 5' end → transcription start site, N-cigar op → splice junction, 3' end → cleavage and polyadenylation (CPA) site.

One pipeline, two sequencing families — **ONT long reads** (DRS, PCR-cDNA) and **NGS short reads** (TruSeq/CORALL v2, QuantSeq REV 3'-end).

<p align="center">
  <img src="docs/figures/pipeline_overview.png#gh-light-mode-only" alt="RECTIFY pipeline overview" width="820">
  <img src="docs/figures/pipeline_overview_dark.png#gh-dark-mode-only" alt="RECTIFY pipeline overview" width="820">
</p>

---

## Inputs

|                       | **ONT DRS**            | **ONT PCR-cDNA**                            | **TruSeq / CORALL v2**                        | **QuantSeq REV**           |
|---                    |---                     |---                                          |---                                            |---                         |
| Family                | ONT long read          | ONT long read                               | NGS short read                                | NGS short read             |
| Chemistry             | Direct RNA             | PCB114.24 (SSP + 27-nt UMI)                 | Whole-transcriptome (CORALL v2: 12-nt UMI at R1 start) | dT-primed **3'-end**       |
| Poly(A) in read?      | Yes                    | Yes (poly-A *or* poly-T, depending on orient) | No                                            | No                         |
| 3' end side of read   | right                  | right (`fwd`) / left (`rev`)                | n/a — reads are internal fragments            | left                       |
| Gene strand vs BAM    | matches                | matches (`fwd`) / inverted (`rev`)          | kit-dependent (stranded kits)                 | **inverted (antisense)**   |

Both strands of each PCR-cDNA amplicon are sequenced; Dorado does **not** reverse-complement antisense reads. About half come out as `orient=fwd` (SSP+UMI at basecalled 5′, poly-A at 3′ — same layout as DRS) and half as `orient=rev` (SSP_RC+UMI_RC at basecalled 3′, poly-T at 5′ — mirror layout). `rectify correct-cdna` normalises both so all downstream stages see a single canonical orientation per read. All four chemistries converge at the multi-aligner stage.

---

## Pipeline

### (a) Per-protocol pre-process

Reads carry non-genomic sequence — basecalled poly(A), template-switching primers, UMIs, adapter stubs — that causes aligners to misplace the exact 5' and 3' ends. Pre-processing strips this material so every downstream step sees a clean mRNA body.

| Protocol     | Step                    | What it does |
|---           |---                      |---           |
| **DRS**      | `rectify trim-polya`    | Strip basecalled poly(A) + adapter from the 3' end; cache as soft-clip metadata for later restoration |
| **ONT cDNA** | `rectify correct-cdna`  | Extract 27-nt UMI, directional cluster reads by molecule, build abPOA consensus per cluster, strip SSP+UMI+GGG (5') and poly(A) (3') |
| **TruSeq / CORALL** | `rectify split --umi` | Slice the CORALL 12-nt UMI off the R1 start (`--umi-location r1-start --umi-length 12`) into an `RX` tag; the COMPASS aligners carry it through, and `rectify umi-dedup` collapses reads to one record per molecule after alignment |
| **QSrev**    | *(none)*                | Reads carry no poly(A); strand inversion is handled at the correction step |

<p align="center">
  <img src="docs/figures/polya_pretrim.png#gh-light-mode-only" alt="DRS poly(A) pre-trim" width="660">
  <img src="docs/figures/polya_pretrim_dark.png#gh-dark-mode-only" alt="DRS poly(A) pre-trim" width="660">
  <img src="docs/figures/cdna_umi_architecture.png#gh-light-mode-only" alt="ONT PCR-cDNA UMI architecture" width="660">
  <img src="docs/figures/cdna_umi_architecture_dark.png#gh-dark-mode-only" alt="ONT PCR-cDNA UMI architecture" width="660">
</p>

### (b) Multi-aligner (parallel) + overhang resolver + junction pool

`rectify align` runs **minimap2 + uLTRA + deSALT** in parallel for long reads. Short reads (`--short-read`) pick their panel by protocol: TruSeq-style whole-transcriptome RNA-seq uses the **COMPASS splice-aware panel** (bbmap + STAR ×2 + HISAT2 ×2, joined by magicblast + gsnap when paired via `--read2`), while dT-primed 3'-end libraries (`--dT-primed-cDNA`, e.g. QuantSeq REV) use **bbmap + bwa** — there, exact 3'-end placement matters more than splice recall. For UMI kits (Lexogen CORALL v2), the R1 UMI extracted at `rectify split --umi` rides through the COMPASS panel as an `RX` tag; `rectify umi-dedup` then collapses PCR duplicates to molecules. Aligners disagree on junctions within reads for several reasons — inherent ambiguity at homopolymer donor/acceptor sites, differing gap-open penalties, different handling of soft-clips — so running them in parallel gives RECTIFY a panel of candidate alignments per read to reconcile. **Each aligner produces its own BAM** with its own junction set; junctions observed across all BAMs are unioned with the annotated splice database to form a shared **junction pool** that the correction step uses to rescue partial alignments.

The panel's fourth arm is RECTIFY's own **overhang resolver** — a native re-alignment stage (Station A of the re-aligner; no external binary) that consumes the finished minimap2 BAM, re-places terminal soft-clipped overhangs across splice junctions, and substitutes the resolved BAM for the raw minimap2 arm. It was inspired by [gapmm2](https://github.com/nextgenusfs/gapmm2)'s insight that a minimap2 alignment can be post-processed to recover terminal splice overhangs, with several deliberate upgrades:

- the search window for each overhang is **bounded by its information content** — an explicit false-discovery budget instead of a fixed genomic window;
- low-information overhangs (poly(A) tails, repeat expansions, short or low-complexity clips) are **refused rather than force-placed** — refusal is a counted, first-class outcome;
- candidate sites come from a **precomputed splice-site index** (binary-search range query, never a scan), with configurable acceptor classes for alternative-3'SS studies (`--resolver-acceptor-classes prp18`) and an opt-in paired AT-AC intron class (`--resolver-atac`; yeast splices AT-AC through its major spliceosome, and in human it is the U12-type class);
- placements are scored with RECTIFY's **chemistry-calibrated homopolymer-aware penalties**, and accepted only when the winner is unambiguous after junction-ambiguity canonicalization — competitors within the score margin must be the *same* junction, not merely a nearby coordinate.

Additional arms (mapPacBio, gapmm2, GMAP, winnowmap2, minisplice) remain available via `--aligners` / `--junction-aligners`. mapPacBio in particular is retained as a discovery scout for non-canonical splicing studies: its motif-free placement can find junctions the index-driven resolver cannot enumerate, and Station C (below) is the gate that adjudicates such candidates.

<p align="center">
  <img src="docs/figures/multi_aligner_consensus.png#gh-light-mode-only" alt="Multi-aligner pipeline" width="720">
  <img src="docs/figures/multi_aligner_consensus_dark.png#gh-dark-mode-only" alt="Multi-aligner pipeline" width="720">
</p>

### (c) `rectify correct` — repair each aligner's output

The three corrections (5' / introns / 3') run **per aligner** so every alignment gets the same rescues on the same standardized scoring before any winner is picked. This way the consensus step in (e) compares apples to apples — an aligner whose raw output happens to soft-clip the 5' end isn't unfairly penalised against one that didn't, because the soft-clip is rescued (or scored consistently) on every lane.

#### 5' ends — junction rescue
Long reads spanning splice junctions often carry soft-clipped 5' bases that match the upstream exon — the aligner couldn't extend through the junction. RECTIFY checks each soft-clipped overhang against the **junction pool** from (b) and, if it matches a nearby donor, extends the alignment across the junction to that donor.

<p align="center">
  <img src="docs/figures/5prime_junction_rescue.png#gh-light-mode-only" alt="5' junction rescue" width="660">
  <img src="docs/figures/5prime_junction_rescue_dark.png#gh-dark-mode-only" alt="5' junction rescue" width="660">
</p>

> *Note:* in direct RNA, the read's 5' end is affected by two factors: (1) 5'→3' degradation before capture, and (2) incomplete sequencing as the RNA translocates through the nanopore — so 5' junction rescue corrects alignment artifacts but does not guarantee full TSS recovery. For ONT cDNA, the presence of a terminal template-switching primer distinguishes full-length from truncated reads.

#### Introns — splice junction refinement and classification
Each N-cigar op is refined against the **junction pool** (annotated + observed-across-aligners) with homopolymer-aware scoring, then classified per read:

| Class (output label)                | Condition                                                       |
|---                                  |---                                                              |
| `unspliced`                         | Read spans an intron but has no N-cigar op                      |
| `annotated`                         | Every N-op matches an annotated donor + acceptor exactly        |
| one-side novel (output `alternative`) | One end of an N-op matches an annotated site, the other does not (a novel donor *or* acceptor) |
| both-side novel (output `novel`)    | Neither end of any N-op matches any annotated site              |

Here **"novel" means "not in the annotation"** — *not* "never observed before." A
junction absent from a given GENCODE *basic* build is often catalogued elsewhere
(recount3/intropolis) or present in the *comprehensive* set. The output TSV keeps the
short labels `alternative`/`novel`; the figure and prose use the more standard
*one-side novel* / *both-side novel* (see `docs/ARCHITECTURE.md` for the mapping).

<p align="center">
  <img src="docs/figures/splice_classification.png#gh-light-mode-only" alt="Splice classification" width="620">
  <img src="docs/figures/splice_classification_dark.png#gh-dark-mode-only" alt="Splice classification" width="620">
</p>

#### 3' ends — two complementary corrections
The aligner can land the 3' end in the wrong place in either direction. RECTIFY runs **read-vs-reference walkback** (upstream — for poly(A) over-extension into genomic A-tracts) and **HP-aware soft-clip rescue** (downstream — for homopolymer under-calls), so the rectified 3' end converges on the true CPA regardless of which artifact caused the original miss.

**Walkback (upstream).** Basecalling errors in the pA tail mean that the trimmed sequences can still contain non-genomic, pA tail fragments which align to genomic regions. Scan the alignment from the 3' end inward. Walk past every A — whether the genome at that position is also A or not — and stop at the first non-A base where read matches reference. This catches **internal priming** at genomic A-tracts, where the read's own poly(A) tail aligns over genomic A's so the apparent 3' end overshoots the true CPA. The pre-trim in step (a) anchors the scan at the boundary of basecalled poly(A), letting walkback focus on genomic A-tract ambiguity rather than the poly(A) tail itself.

<p align="center">
  <img src="docs/figures/walkback_readvsref.png#gh-light-mode-only" alt="Read-vs-reference walkback" width="680">
  <img src="docs/figures/walkback_readvsref_dark.png#gh-dark-mode-only" alt="Read-vs-reference walkback" width="680">
</p>

**HP-aware soft-clip rescue (downstream).** Nanopore basecallers systematically under-call long homopolymers — a genomic 9-T tract may be called as only 6 T's, leaving the next 3 reference bases (often the CPA site itself, plus its immediate context) **soft-clipped** by the aligner because they no longer have a contiguous home in the read. RECTIFY identifies reads whose alignment terminates inside a reference homopolymer with a soft-clipped overhang, and checks whether the overhang matches the genome immediately downstream of the HP. If it does, the alignment is **extended through the HP-gap** (deletion ops for the under-called positions) so the soft-clip is recovered as a proper match — moving the 3' end **further downstream** to the true CPA. This is the opposite directional fix from walkback and runs on the same per-aligner BAMs (`indel_corrector.rescue_softclip_at_homopolymer`).

<p align="center">
  <img src="docs/figures/softclip_rescue.png#gh-light-mode-only" alt="HP-aware soft-clip rescue" width="700">
  <img src="docs/figures/softclip_rescue_dark.png#gh-dark-mode-only" alt="HP-aware soft-clip rescue" width="700">
</p>

### (d) Empirical error-rate scoring

When refinement chooses between candidate junctions or end positions, edit operations (mismatch, insertion, deletion) are scored by their **empirical per-chemistry error rates** as a function of homopolymer context and base class (AT vs CG). For Nanopore R10.4.1 reads, a deletion at HP=8 carries a penalty of ~0.032 — nearly free — versus 0.44 at HP=1, because Nanopore basecallers routinely under-call long homopolymer runs; penalties reach a floor around HP ≥ 9 (AT) and HP ≥ 7 (CG). Each chemistry uses its own calibrated table (DRS, ONT cDNA, QuantSeq REV), derived from *S. cerevisiae* (R10.4.1) reads. **Human *H. sapiens* DRS tables are now bundled** (GM12878 IVT RNA004): AT penalties closely match yeast, and both species show the same floor behavior.

<p align="center">
  <img src="docs/figures/hp_scoring.png#gh-light-mode-only" alt="Empirical HP deletion penalty curves" width="620">
  <img src="docs/figures/hp_scoring_dark.png#gh-dark-mode-only" alt="Empirical HP deletion penalty curves" width="620">
</p>

See [docs/EMPIRICAL_HP_PENALTY_SCORING.md](docs/EMPIRICAL_HP_PENALTY_SCORING.md).

### (e) Consensus → rectified BAM

Now that every aligner's output has been independently rescued and scored on the same scale, the consensus step picks the **best corrected alignment per read** — primary key: lowest homopolymer-aware edit distance over the corrected CIGAR (the same empirical error model as (d)), with alignment span as the tie-breaker. Optional chimeric reconstruction stitches complementary junction sets across aligners to recover junctions no single aligner produced on its own. Output: one BAM where every record has been independently rescued, scored, and reconciled across the full aligner panel.

### (f) Re-aligner stations B and C

The rectified BAM then passes through the remaining two stations of RECTIFY's **re-aligner architecture**. All three stations share one contract: a *constrained proposer* plus an *evidence arbiter*, with refusal as a first-class outcome — no station ever force-places an alignment the read evidence cannot support.

- **Station B — consensus triage** (`rectify triage`; opt-in via `run-all --triage`). Classifies every corrected alignment from read evidence alone (junction-proximal errors, unexplained soft clips, annotation status). Clean alignments are bypassed untouched; only the distressed minority is re-aligned (motif-blind junction refinement over the full evidence pool, plus optional clip legs that reuse the overhang resolver and the 5' rescue machinery). A re-aligned candidate replaces the original **only if it strictly improves the homopolymer-aware edit distance** — the re-aligner proposes, the arbiter disposes.
- **Station C — pool-level junction gate** (`rectify pool-gate`; on by default in `run-all`, report-only). Judges each junction over its whole pool of supporting reads: overhang quality, repeat-context flags (rRNA/Ty/LTR annotation + genome self-homology), an annotation-derived intron-length pre-gate, and separate admission bars for canonical vs non-canonical junctions. It annotates — writing a per-junction admission table (`admit_candidate` / `review` / `demote_orthogonal_evidence`) — and never deletes a read or a junction.

---

## `rectify analyze` — downstream

Run on the rectified BAM.

### 5' and 3' end clustering
A valley-based adaptive clustering algorithm groups nearby corrected 5' ends (TSS) and 3' ends (CPA) into discrete usage clusters. The same algorithm runs separately on each end.

<p align="center">
  <img src="docs/figures/adaptive_clustering.png#gh-light-mode-only" alt="Adaptive clustering" width="600">
  <img src="docs/figures/adaptive_clustering_dark.png#gh-dark-mode-only" alt="Adaptive clustering" width="600">
</p>

### Splice junction analysis
Per-read splice classifications from step (c) are aggregated into per-junction tables: annotated / one-side-novel / both-side-novel counts (output labels `annotated`/`alternative`/`novel`) and junction-shift frequencies. The novel (one-side + both-side) fraction is the readout for the NMD-AS isoform audit (surveillance-mutant studies where unproductive transcripts accumulate).

### Isoform characterization (ONT cDNA only)
`rectify cdna-analyze` uses the matched 5' + 3' coordinates available from full-length cDNA reads to assemble per-read isoforms. Type-1 reads (full-length, UMI captured) cluster by both TSS and CPA; Type-2 reads (truncated, no UMI) cluster by CPA only. 5' and 3' end variation within each cluster is bounded by a configurable window (default ±5 bp). Same-molecule Type-1 ↔ Type-2 cluster pairs are linked by gene + 3' end proximity, recovering deep coverage that would otherwise be discarded as random truncation noise.

<p align="center">
  <img src="docs/figures/cdna_isoform_clustering.png#gh-light-mode-only" alt="cDNA isoform clustering" width="720">
  <img src="docs/figures/cdna_isoform_clustering_dark.png#gh-dark-mode-only" alt="cDNA isoform clustering" width="720">
</p>

---

## Differential analysis (DESeq2)

The corrected BAMs feed four independent DESeq2 analyses. A gene that looks flat at the gene level can show a clear 3' shift or a switch in intron retention — the four views are complementary, not redundant.

| Level                  | What changes                          | Status     |
|---                     |---                                    |---         |
| **Gene**               | Total expression (classic DE)         | shipping   |
| **3' end cluster**     | Differential CPA usage (APA shifts)   | shipping   |
| **5' end cluster**     | Differential TSS usage                | planned (aggregation infra present) |
| **Splice junction**    | Differential intron usage             | planned (aggregation infra present) |

---

## Quick start

> Full quickstarts: [ONT DRS](docs/quickstart.md) · [ONT cDNA](docs/quickstart_cdna.md) · [QuantSeq REV](docs/quickstart_quantseq_rev.md)

```bash
pip install rectify-rna
```

### ONT DRS

```bash
# One command — bundled yeast genome, full pipeline + DESeq2
rectify run-all reads.bam --drs --organism yeast --output-dir results/

# Step by step (any organism)
rectify trim-polya reads.bam                  --genome genome.fa -o trimmed.bam
rectify align      trimmed.bam                 --genome genome.fa --annotation genes.gff -o aligned/
rectify correct    aligned/trimmed.rectified.bam --genome genome.fa --annotation genes.gff -o corrected.tsv
rectify analyze    corrected.tsv              --annotation genes.gff -o results/
```

### ONT PCR-cDNA (PCB114.24)

The cDNA pipeline runs as three explicit stages (Stage 1 deduplicates by UMI
*before* alignment, so there is no single `run-all` shortcut):

```bash
rectify correct-cdna  pcb114.bam                      --reference genome.fa --gff genes.gff -o out/
rectify align         out/stage1_consensus.fastq.gz   --genome genome.fa --prefix stage1 -o out/
rectify cdna-analyze  out/stage1.rectified.bam        --reference genome.fa --gff genes.gff -o out/
```

### QuantSeq REV

```bash
# One command
rectify run-all reads.fastq.gz --short-read --dT-primed-cDNA --genome genome.fa \
    --annotation genes.gff --output-dir results/

# Step by step — pass FASTQ to `rectify align --short-read` first;
# piping straight to `rectify correct` bypasses the proper short-read aligner panel.
rectify align   reads.fastq.gz --short-read --dT-primed-cDNA --genome genome.fa --prefix sample -o aligned/
rectify correct aligned/sample.rectified.bam --short-read --dT-primed-cDNA --genome genome.fa -o corrected.tsv
rectify analyze corrected.tsv  --annotation genes.gff -o results/
```

### Multi-sample and HPC

```bash
# Manifest mode — run all samples in one invocation
rectify run-all --manifest samples.tsv --genome genome.fa --annotation genes.gtf -o results/

# Large datasets: split BAM into chunks and run as a SLURM/SGE job array, then merge.
# See docs/quickstart.md for an HPC job-array template.
```

**Bundled for *S. cerevisiae*:** genome, SGD annotations, GO terms, WT NET-seq, 64K pre-computed A-tract CPA sites. No external files needed for yeast.

---

## Installation

```bash
pip install rectify-rna                                            # core
pip install rectify-rna[visualize]                                 # +plots (matplotlib, seaborn)
conda install -c kevinrjroy -c conda-forge -c bioconda rectify-rna # +MEME for motif discovery
```

The `[visualize]` extra adds `matplotlib` and `seaborn` for metagene plots, genome-browser figures, and heatmaps. The cDNA UMI-consensus extras (`pip install rectify-rna[cdna-correct]`) add `edlib` + `pyabpoa`.

---

## Output

`rectify correct` writes a per-read `corrected_3ends.tsv` (one row per read). Selected columns:

```
read_id │ chrom │ strand │ original_3prime │ corrected_3prime │ polya_length │ confidence │ qc_flags
read001 │ chrI  │   +    │     147592      │      147585      │      42      │    high    │   PASS
read002 │ chrI  │   +    │     147594      │      147591      │      38      │   medium   │   PASS
read003 │ chrII │   +    │     283109      │      283104      │      31      │    low     │ AG_RICH
```

`rectify analyze` produces (per condition): `cpa_clusters.tsv`, `deseq2_genes_<condition>.tsv`, `deseq2_clusters_<condition>.tsv`, `shift_analysis_<condition>.tsv` (APA shifts), `go_enrichment_up_<condition>.tsv` / `go_enrichment_down_<condition>.tsv`, and `motif_summary_<condition>.tsv`.

---

## Commands

| Command                | Description                                                                    |
|---                     |---                                                                             |
| `rectify run-all`      | End-to-end pipeline (manifest support, provenance, step-skip)                  |
| `rectify trim-polya`   | DRS pre-process: strip poly(A) + adapter                                       |
| `rectify correct-cdna` | ONT cDNA Stage 1: UMI cluster + abPOA consensus + pre-trim                     |
| `rectify align`        | Multi-aligner panel + overhang resolver (long-read or `--short-read`)          |
| `rectify correct`      | 5'/intron/3' end correction per read                                           |
| `rectify cdna-analyze` | ONT cDNA Stage 3: post-align isoform assembly + T1↔T2 linkage                 |
| `rectify triage`       | Station B: classify corrected alignments, re-align the distressed minority     |
| `rectify pool-gate`    | Station C: per-junction admission report over the whole read pool              |
| `rectify analyze`      | Clustering, DESeq2, GO, motifs                                                 |
| `rectify attribute-reads` | Per-read gene-attribution sidecar (readthrough-aware)                       |
| `rectify qc`           | Per-library sequencing QC (N50, base quality, error rate, clips, tails)        |
| `rectify umi-dedup`    | Collapse a UMI-tagged (`RX`) short-read BAM to one record per molecule         |
| `rectify export`       | Corrected positions → bigWig/bedGraph                                          |
| `rectify aggregate`    | Reads → 3'-end / 5'-end / junction datasets                                    |
| `rectify netseq`       | NET-seq A-tract refinement (see below)                                         |

Full reference at [readthedocs](https://rectify-rna.readthedocs.io).

---

## Species-specific optimizations

### NET-seq A-tract refinement (yeast)

Nascent Elongating Transcript Sequencing data captures RNA polymerase II–associated transcripts at single-nucleotide resolution and provides direct ground-truth evidence of CPA intermediates within genomic A-tracts. RECTIFY uses NNLS deconvolution with a point-spread function derived from 5000+ zero-A sites to recover true cleavage positions where DRS alone is ambiguous.

<p align="center">
  <img src="docs/figures/netseq_deconvolution.png#gh-light-mode-only" alt="NET-seq oligo(A) deconvolution" width="760">
  <img src="docs/figures/netseq_deconvolution_dark.png#gh-dark-mode-only" alt="NET-seq oligo(A) deconvolution" width="760">
</p>

```bash
# Bundled WT yeast NET-seq is auto-detected
rectify correct reads.bam --organism yeast -o corrected.tsv

# Other organisms or mutant conditions
rectify correct reads.bam --genome genome.fa --netseq-dir my_netseq/ -o corrected.tsv
```

For organism-specific poly(A) models and custom A-tract priors, see [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and the [user guide](https://rectify-rna.readthedocs.io).

---

## Documentation

- **Algorithms** — [docs/algorithms/](docs/algorithms/) · [empirical HP scoring](docs/EMPIRICAL_HP_PENALTY_SCORING.md) · [3' indel correction](docs/algorithms/3prime_indel_correction.md) · [multi-aligner consensus](docs/algorithms/multi_aligner_consensus.md)
- **Architecture / internals** — [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) (SAM tag namespace, walkback wrappers, module call graph)
- **Quickstarts** — [DRS](docs/quickstart.md) · [ONT cDNA](docs/quickstart_cdna.md) · [QuantSeq REV](docs/quickstart_quantseq_rev.md)

---

## Citation

> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020;176:4-13. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 0.9.0** (first public release): manuscript in preparation. · [CITATION.cff](CITATION.cff)

> RECTIFY is pre-1.0; the CLI and output schema may change before the 1.0 release.

---

## License

MIT — see [LICENSE](LICENSE).

## Contact

Kevin R. Roy — [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com) · GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)
