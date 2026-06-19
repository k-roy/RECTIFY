# ONT Direct-RNA-seq & cDNA Alignment Ecosystem — Error Model, 3'-End Precision, and the RECTIFY Problem

**Investigation date:** 2026-06-18
**Scope:** The Oxford Nanopore (ONT) direct-RNA-seq (DRS) and cDNA alignment ecosystem, with emphasis on the error model and 3'-end / cleavage-and-polyadenylation (CPA) precision challenges that motivate downstream correction in RECTIFY.
**Anchoring:** RECTIFY aligns ONT yeast RNA with an ensemble (minimap2 / mapPacBio[BBMap] / gapmm2 / deSALT / uLTRA), then corrects 3' ends (poly-A walkback Module 2E, indel correction 2C, soft-clip rescue 2G, junction refinement 2H). The DRS path pre-trims poly-A+adapter (Step 0, `trim_drs_bam_polya`) and restores it as soft-clips post-correction (Step 4). DRS is the default protocol; `--dT-primed-cDNA` (QuantSeq / oligo-dT) enables the AG-mispriming module. Coordinates are 0-based half-open (pysam/BED convention).

**FACT / INFERENCE convention:** Claims sourced to a cited paper, tool doc, or repo are marked **[FACT]**. Claims that are reasoned syntheses across sources or extrapolations to RECTIFY's situation are marked **[INFERENCE]**.

---

## 1. ONT RNA Error Model (homopolymer / poly-A / junction)

### 1.1 Signal → basecall path

ONT sequencing measures ionic-current disruption as a nucleic-acid strand is **ratcheted** base-by-base through a protein nanopore by a motor enzyme. The motor controls translocation speed; published phi29-controlled ratcheting rates are on the order of 2.5–40 nt/s, and the current signal reflects a *k-mer* (multiple bases simultaneously occupying the pore constriction), not a single base **[FACT]**. Basecalling (guppy historically, now **dorado**, ONT's production basecaller) converts the raw current squiggle into a base sequence using a neural network **[FACT]**.

Because the readout is a k-mer-level current, the central algorithmic difficulty is **segmentation**: deciding how many identical bases produced a flat current segment. This is the root cause of homopolymer error.

### 1.2 R9.4.1 vs R10.4.1 pore error profiles

| Property | R9.4.1 | R10.4.1 |
|---|---|---|
| Reader head | single | **dual reader head**; elongated pore [FACT] |
| Typical raw error (matched kits) | ~5–6% (V10/Kit10) [FACT] | **<2%** (V14/Kit12/Kit14), model-dependent [FACT] |
| Homopolymer 4–9 bp accuracy | lower | **higher** on 4–9 bp runs [FACT] |
| Base-class asymmetry | — | A/T homopolymers basecalled **more accurately** than C/G homopolymers [FACT] |

The R10.4.1 dual-reader-head design specifically targets homopolymer resolution by sampling the strand at two points, improving the segmentation of long runs **[FACT]**. Even so, homopolymer deletion is the dominant residual error class on both chemistries **[FACT]**.

> **Note on the CLAUDE.md AT-vs-CG claim.** CLAUDE.md's empirical penalty table states **AT runs have ~10–20% higher *deletion* rates than CG runs** at equal HP length, and the D-penalty table (D/AT HP=1 = 0.37 vs D/CG HP=1 = 0.58) encodes *lower* penalty (i.e., *more expected* deletions) for AT. The literature I found reports the **basecalling-accuracy** direction the opposite way for *substitution/overall* accuracy — A/T homopolymers are *basecalled more accurately* than C/G overall [FACT]. These are **not contradictory** **[INFERENCE]**: the literature's "A/T more accurate" is an aggregate accuracy statement (driven heavily by substitutions and short runs), whereas RECTIFY's table isolates the **deletion operation** specifically. AT homopolymers having a *higher deletion-specific rate* (pore ratcheting through low-variance A/T current passes through identical k-mers faster, dropping bases) while still having *higher overall accuracy* (fewer substitutions) is internally consistent. The two statements measure different error operations. RECTIFY's table is the operative one for 3'-end correction because poly-A walkback and HP-anchored junction scoring care about **deletions in A/T context**, exactly where the table says risk is highest.

### 1.3 Why poly-A tails and homopolymers are systematically mis-resolved

- **Dwell-time / segmentation ambiguity.** A long poly(A) tract produces a "characteristic, relatively low-variance current segment" whose *duration* (dwell time) encodes length, but the basecaller cannot reliably count identical k-mers from a flat signal [FACT]. Signal-level tools (nanopolish polya, tailfindr, NanoTimer) estimate HP/tail length from **dwell time normalized to a calibration strand** precisely because the basecall sequence loses this information [FACT].
- **Ratcheting asymmetry.** Translocation speed is base-/context-dependent; A/T-rich k-mers pass faster and are under-segmented, biasing toward **deletion** in A/T runs **[INFERENCE, consistent with FACT on ratcheting + RECTIFY table]**.
- **GC homopolymers.** GC-rich tracts >3 bp are independently prone to mis-calling and false deletion [FACT]; this is the substitution/mis-call axis rather than the count-loss axis.

For RECTIFY this means the **3' boundary of an aligned read is intrinsically fuzzy**: the poly-A tail (DRS) or the A-run immediately upstream of the CPA is where deletion error concentrates, so the aligner's reported `reference_end` is noisy by several bp **[INFERENCE]**. This is the precise error that Module 2E (A-tract walk-back), 2C (indel correction, requires MD tags), and the empirical penalty tables address.

### 1.4 dorado poly-A estimation and the `pt:i` tag

`dorado basecaller --estimate-poly-a` (disabled by default) estimates poly(A)/poly(T) tail length for cDNA (PCS/PCB kits) and RNA, configurable for custom primers, interrupted tails, and plasmids [FACT]. The estimate is written to the **`pt:i` tag** of each BAM record [FACT]:

- `pt:i:N` (N>0) — estimated tail length in nt
- `pt:i:0` — primer anchor found but length not estimable
- `pt:i:-1` — primer anchor for the tail not found [FACT]

Independent benchmarking with synthetic RNA finds dorado the preferred poly(A) estimator (fast, low mean error, integrated with basecalling) [FACT]. **RECTIFY consequence:** CLAUDE.md correctly stores `pt` lengths in parquet trim metadata rather than `samtools fastq -T pt` (which embeds `@UUID\tpt:i:N` and corrupts QNAMEs); the `pt:i` value is the dorado tail estimate, and RECTIFY's Step 0/Step 4 trim-then-restore is built around it.

---

## 2. DRS vs cDNA Alignment Differences

### 2.1 DRS is sequenced 3'→5' with the poly-A in-read

DRS proceeds **3'→5'**, starting from the adapter ligated to the **poly(A) tail end**; sequencing runs from the poly(T) adapter through the poly(A) tail into the RNA body toward the 5' cap [FACT]. Therefore:

- **The poly(A) tail is physically present in the read** (DRS). This is exactly CLAUDE.md's "DRS — the poly-A tail IS in the read" default.
- **No priming step exists** — the native RNA is sequenced directly, so oligo-dT internal mispriming **cannot occur** in DRS. This is why RECTIFY disables the AG-mispriming module for DRS and enables it only with `--dT-primed-cDNA` [FACT alignment with CLAUDE.md].

### 2.2 5' truncation from motor fall-off

The motor enzyme releases the strand **~10–12 nt (≈11 nt) before the 5' terminus**; the final 5' bases are not threaded through the pore and are lost [FACT]. Researchers report being unable to resolve the final ~6 nt of the 5' end, and there is no native marker for the true 5' end [FACT]. **RECTIFY consequence:** 5' ends are *systematically truncated*, which is why RECTIFY's 5' work (Cat3 `rescue_3ss_truncation`, `local_aligner.py` exon-CIGAR) is about *recovering* a truncated/soft-clipped 5' splice boundary, not trusting the aligned 5' coordinate. The 3' end (CPA) is the biologically meaningful, well-anchored end in DRS — which is why RECTIFY's primary product is 3'-end/CPA precision.

### 2.3 dT-primed cDNA (QuantSeq / oligo-dT): poly-A removed, mispriming added

In dT-primed cDNA the oligo-dT primer hybridizes at the **start of the poly(A) tract**; sequencing starts at the first non-A base, so **the poly-A tail is NOT in the read** [FACT alignment with CLAUDE.md]. Two consequences:

1. The A-content distinguishing real CPA from genomic context is **downstream of the reported 3' end** (off-read), so mispriming detection must query the **genome** (`get_downstream_sequence()`), not read sequence — RECTIFY's design.
2. **Internal/intra-priming** becomes a major artifact: the oligo-dT (or RT) misprimes onto genomic A-rich runs (e.g., ≥6 consecutive A's, or downstream A-content >50%), producing **false CPA sites** that align without splice junctions to A-rich loci [FACT]. This is the artifact RECTIFY's AG-mispriming module targets, and benchmarks confirm internal-priming reads are a "substantial fraction" of long-read transcriptome data [FACT].

### 2.4 Summary table

| Axis | DRS (default) | dT-primed cDNA (`--dT-primed-cDNA`) |
|---|---|---|
| Direction | 3'→5' [FACT] | cDNA, library-dependent |
| Poly-A in read? | **Yes** [FACT] | **No** [FACT] |
| Internal mispriming? | No (no priming step) [FACT] | **Yes** (oligo-dT/RT on A-runs) [FACT] |
| 5' end | truncated ~11 nt by motor [FACT] | dependent on full-length capture (pychopper) |
| RECTIFY AG module | off | **on** |
| RECTIFY poly-A trim (2B) | on | on |

---

## 3. Community Pipelines & Which Aligners They Use

Nearly the entire ONT RNA ecosystem standardizes on **minimap2** as the genome/transcriptome aligner; tools differ in *what they do with the alignment*, not in the aligner itself **[FACT, multiple sources]**.

| Tool | Role | Aligner used | 3'-end handling |
|---|---|---|---|
| **dorado aligner** | basecaller-integrated alignment | **wraps minimap2 internally** [FACT] | inherits minimap2 soft-clip behavior |
| **minimap2** | de facto genome/transcriptome aligner | itself; `-ax splice -uf -k14` for DRS [FACT] | soft-clips poly-A; 3' end = `reference_end` |
| **pychopper** | cDNA full-length detection/orientation | preprocessor (no align) → feeds minimap2/FLAIR [FACT] | identifies/orients full-length cDNA before alignment |
| **FLAIR** | isoform definition/quant | **minimap2** → BED12; collapse groups by splice junctions and defines common TSS/TES [FACT] | "common transcription end sites" — **clusters** 3' ends, does not refine per-read |
| **FLAMES** | full-length isoform + sc | minimap2 | isoform-level |
| **NanoCount** | DRS expression quant | **minimap2 to transcriptome**, primary+secondary, EM | transcriptome-relative; 3'-end precision not its goal [FACT] |
| **Bambu** | context-aware quant/discovery | **minimap2** `-ax splice` → Bambu [FACT] | isoform-level counts |
| **Sicelore** | single-cell long-read | minimap2; scans fastq for ≥15 nt polyA/T, ≥75% A within 100 nt of ends [FACT] | poly-A used for barcode/UMI, not CPA precision |
| **Nano3P-seq** | 3'-end + tail dynamics (end-capture cDNA) | minimap2 | **3'-end-aware by protocol design** [FACT] |

**Key gap (the RECTIFY niche) [INFERENCE]:** every mainstream pipeline either (a) treats the aligner-reported 3' end as ground truth, (b) *clusters* 3' ends into TES/PAS bins (FLAIR, LAPA), or (c) works at the isoform/transcriptome level where exact CPA position is averaged out. **None performs per-read 3'-end *correction*** against the homopolymer/indel/soft-clip error model. RECTIFY's ensemble-then-correct design fills exactly this gap.

---

## 4. Why 3'-Ends Are Ambiguous (the RECTIFY problem)

### 4.1 minimap2 `-ax splice -uf` is the de facto standard

For DRS the canonical command is `minimap2 -ax splice -uf -k14 ref reads` [FACT]:
- `-ax splice`: splice-aware (introns as `N` CIGAR ops)
- `-uf`: force forward-transcript-strand mapping (DRS is stranded) [FACT]
- `-k14`: smaller k-mer for sensitivity to first/last exons on noisy reads [FACT]

minimap2 aligns 10M human DRS reads in <1 wall-hour on 16 cores with 94.2% of junctions annotation-consistent [FACT] — fast and good on *junctions*, but junction-consistency ≠ 3'-end-precision. (RECTIFY's production command adds `-G 5000 --splice-flank=no --secondary=no --MD --junc-bed --junc-bonus 9` for yeast; `--splice-flank=no` disabling the GT-AG bonus is explicitly chosen for 3'-end accuracy per CLAUDE.md.)

### 4.2 The four sources of 3'-end ambiguity

1. **Homopolymer deletion at the A-tract.** As in §1, the A-run immediately 5' of the CPA loses bases to deletion, so `reference_end` lands a few bp short/long. RECTIFY Module 2E walks back over the A-tract using the *reference* genome (length-independent), and 2C corrects indels where MD tags exist.
2. **Poly-A soft-clipping ambiguity.** minimap2 community guidance is to **soft-clip poly-A**, but users repeatedly question whether the soft-clip boundary (and thus the 3' coordinate) "can be trusted" [FACT]. The soft-clip start is itself error-prone — RECTIFY Module 2G (`rescue_softclip_at_homopolymer`) detects ≥3 bp 3' soft-clips adjacent to a genomic homopolymer and *extends* the 3' end outward, then 2G takes priority over 2E to prevent cancelling corrections.
3. **Internal priming (cDNA).** False 3' ends on genomic A-runs (downstream A >50%, or ≥6 consecutive A) [FACT] — handled by the AG module under `--dT-primed-cDNA`.
4. **Junction imprecision near the 3' terminal exon.** Aligners place intron `N`-op boundaries imprecisely; in noisy A/T context the donor/acceptor can slide. RECTIFY Module 2H (`junction_refiner.py`) rescores every `N`-op with HP-anchored semi-global DP, using the empirical penalty tables and a 0.5 canonical-HP prior.

### 4.3 Why an ensemble (deSALT / BBMap) recovers what minimap2 misses

- **deSALT** uses a **de Bruijn-graph-based two-pass strategy**: it builds graph-based alignment skeletons, sensitively infers exons (including small exons), generates a spliced reference, and re-aligns — explicitly designed for "small exons, serious sequencing errors, and consensus spliced alignment," producing **homogeneous full-length alignments** [FACT]. Where minimap2 truncates or mis-places a noisy terminal exon, deSALT's exon-recovery step can place it correctly. (RECTIFY reports deSALT winning **78.9%** of reads in correct-first consensus.)
- **uLTRA** is a two-pass collinear-chaining aligner (standalone or minimap2 wrapper) with **substantially higher accuracy on small exons** on simulated/synthetic data, and recovers novel exon structures other aligners miss [FACT].
- **mapPacBio (BBMap)** and **gapmm2** contribute orthogonal error/gap models; the consensus selects per-read winners on **corrected** 3' positions (correct-first), not raw scores.

**[INFERENCE]** No single aligner is uniformly best at the noisy 3' terminal exon; the failure modes are *complementary* (minimap2 over-extends/soft-clips poly-A; deSALT/uLTRA recover small/terminal exons). RECTIFY's ensemble-then-correct-then-select architecture is the correct response to a problem where the *aligner choice itself* is read-dependent. The empirical win rates (deSALT 78.9 / mapPacBio 18.2 / uLTRA 2 / gapmm2 0.8 / minimap2 0.1) quantify this complementarity.

### 4.4 Chimera artifacts

A 2024 genomic-language-model study finds DRS produces **chimera artifacts** (two molecules basecalled as one), an additional 3'/junction confounder [FACT]. RECTIFY's chimeric-consensus handling (`build_chimeric_cigar`, `D` vs `N` for ≤10 bp gaps) and the v3.3.0 overhang filter address this class.

---

## 5. State-of-the-Art 3'-End / CPA / APA Tools

| Tool | Input | 3'-end strategy | Precision claim |
|---|---|---|---|
| **LAPA** | long- & short-read | clusters 3' ends → PAS; APA quant [FACT] | cluster/peak-level, not per-read |
| **APALORD** | 3'-primed long-read | leverages "precise 3' end information," quantifies PAS usage per site [FACT] | site-level on 3'-primed LR |
| **scPAISO** | single-cell | identifies mRNA 3' ends from Read1; novel PAS; **~½ of PASs span <20 bp** [FACT] | <20 bp PAS width |
| **Nano3P-seq** | end-capture cDNA | protocol-level 3'-end + tail dynamics at single-molecule resolution [FACT] | molecule-level tails |
| **Sicelore** | single-cell ONT | poly-A/T scan for barcode/UMI [FACT] | not CPA-precision |
| **nanopolish polya / tailfindr / NanoTimer** | signal-level | dwell-time tail length vs calibration strand [FACT] | tail *length*, not CPA coordinate |
| **RECTIFY** | DRS/cDNA BAM | **per-read 3'-end correction** (2E/2C/2G/2H) + ensemble consensus | bp-level per-read CPA |

### 5.1 Orthogonal validation: NET-seq / QuantSeq / 3'-Seq

Benchmarking of APA methods concludes that **3'-Seq and Iso-Seq quantify PAS usage more reliably than computational tools on short-read data** [FACT] — i.e., a *dedicated 3'-end assay* is the gold standard. RECTIFY's use of **NET-seq tables** (`saccharomyces_cerevisiae_netseq_*.tsv.gz`) as an orthogonal refinement/validation signal, and **QuantSeq REV (dT-primed)** as the short-read protocol, is methodologically aligned with the field's best practice **[INFERENCE]**: correct the aligner-reported CPA, then refine/validate against an orthogonal 3'-end signal. CLAUDE.md's NET-seq refinement (e.g., signal=75.0 at a corrected boundary) is the in-pipeline instance of this.

### 5.2 Where RECTIFY sits relative to SOTA

**[INFERENCE]** The published SOTA operates one level *above* RECTIFY: LAPA/APALORD/scPAISO consume 3' ends and *cluster* them into PAS; they assume the per-read 3' coordinate is approximately correct. RECTIFY operates one level *below*: it makes the per-read 3' coordinate correct *before* clustering, by modeling the homopolymer/indel/soft-clip/junction error explicitly and by selecting the best aligner per read on corrected positions. The two are complementary — RECTIFY's corrected per-read CPA is the ideal *input* to a LAPA-style PAS clustering step. No surveyed tool combines (a) multi-aligner ensemble, (b) explicit per-read 3'-end error correction against an empirical Nanopore HP penalty model, and (c) correct-first aligner selection. That combination is RECTIFY-specific.

---

## 6. Key Takeaways for RECTIFY

1. **The error model justifies the architecture.** Homopolymer deletion concentrated in A/T context (RECTIFY's penalty tables) is the documented dominant residual error of even R10.4.1; it lands exactly at the poly-A/CPA boundary, making the aligner-reported 3' end intrinsically noisy. Per-read correction (2E/2C/2G) is the correct response, not better aligner flags alone. **[INFERENCE from FACTs in §1.]**
2. **DRS 5' truncation (~11 nt motor fall-off) is why 5' work is *rescue*, not trust** — and why 3' (CPA) is the meaningful, well-anchored end in DRS. **[FACT §2.2.]**
3. **The `--dT-primed-cDNA` vs DRS split is biologically grounded:** internal/intra-priming on genomic A-runs is real and substantial in primed protocols but impossible in DRS [FACT], exactly matching RECTIFY's module activation table.
4. **No mainstream pipeline corrects per-read 3' ends.** They cluster (FLAIR/LAPA), quantify (NanoCount/Bambu), or estimate tail length (signal tools). RECTIFY's niche is real and unfilled. **[INFERENCE from §3, §5.]**
5. **Ensemble is justified by complementary aligner failure modes** at the noisy terminal exon (deSALT/uLTRA small-exon recovery vs minimap2 poly-A over-extension/soft-clip). **[FACT on tool capabilities; INFERENCE on complementarity.]**
6. **NET-seq/QuantSeq orthogonal validation matches field best practice** (dedicated 3'-end assays are the APA gold standard). **[FACT §5.1.]**

---

## References

ONT error model & chemistry:
- Benchmarking of Nanopore R10.4 and R9.4.1 flow cells — https://pmc.ncbi.nlm.nih.gov/articles/PMC10070092/ (homopolymer 4–9 bp; A/T vs C/G accuracy)
- Comparison of R9.4.1/Kit10 and R10/Kit12 ONT flowcells — https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000910 (5–6% → <2% error)
- Evaluation of ONT R10.4.1 bacterial genome reconstruction — https://pmc.ncbi.nlm.nih.gov/articles/PMC11170131/
- Sequencing DNA with nanopores: Troubles and biases (PLOS ONE) — https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0257521 (GC-rich homopolymer deletion; low- vs high-GC error)
- Systematic and stochastic influences on MinION across nucleotide bias — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5816649/
- Calling Homopolymer Stretches from Raw Nanopore Reads by Analyzing k-mer Dwell Times — https://link.springer.com/chapter/10.1007/978-981-10-5122-7_61

dorado poly-A / basecalling:
- Dorado Poly(A) Estimation docs — https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/ (pt:i tag; -1/0/N semantics; --estimate-poly-a)
- nanoporetech/dorado (aligner wraps minimap2; --estimate-poly-a) — https://github.com/nanoporetech/dorado
- Using synthetic RNA to benchmark poly(A) length inference from DRS — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12406214/
- TAILcaller (Bioinformatics Advances) — https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf235/8266338

DRS direction / 5' truncation / poly-A:
- Advances in nanopore direct RNA sequencing — https://pmc.ncbi.nlm.nih.gov/articles/PMC11388133/
- Identification of high-confidence human poly(A) RNA isoform scaffolds (motor fall-off ~11 nt) — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8906549/
- Nanopore DRS maps Arabidopsis mRNA processing & m6A (eLife) — https://elifesciences.org/articles/49658

Aligners:
- minimap2 (Li 2018; -ax splice -uf -k14; 94.2% junctions) — https://arxiv.org/pdf/1708.01492 , https://github.com/lh3/minimap2
- deSALT (de Bruijn two-pass; small-exon/full-length) — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6913027/ , https://github.com/ydLiu-HIT/deSALT
- uLTRA (two-pass collinear chaining; small exons) — https://academic.oup.com/bioinformatics/article/37/24/4643/6327681
- 2passtools (ML-filtered two-pass junctions) — https://link.springer.com/article/10.1186/s13059-021-02296-0
- minimap2 iso-seq soft-clip & polyA discussion — https://github.com/lh3/minimap2/issues/459

Community pipelines:
- FLAIR docs (minimap2 → BED12; TSS/TES) — https://flair.readthedocs.io/en/latest/modules.html
- NanoCount (DRS quant, minimap2, EM) — https://academic.oup.com/nar/article/50/4/e19/6439677
- Bambu (context-aware quant; minimap2 -ax splice) — https://genomebiology... (ResearchGate PDF) https://www.researchgate.net/publication/371505702
- Sicelore (poly-A/T scan ≥15 nt, ≥75% A) — https://github.com/ucagenomix/sicelore
- Nano3P-seq (end-capture cDNA; 3'-end + tails) — https://www.nature.com/articles/s41592-022-01714-w

3'-end / CPA / APA tools & chimera:
- LAPA (long/short-read APA) — https://www.biorxiv.org/content/10.1101/2022.11.08.515683
- APALORD (3'-primed long-read APA) — https://www.biorxiv.org/content/10.1101/2025.06.11.658931v1.full
- scPAISO (single-cell PAS; <20 bp width) — https://www.biorxiv.org/content/10.1101/2025.08.20.669565
- Benchmarking APA methods (3'-Seq/Iso-Seq gold standard) — https://link.springer.com/article/10.1186/s13059-021-02502-z
- Internal priming / intra-priming artifact (downstream A >50%; ≥6 consecutive A) — https://elifesciences.org/articles/49658
- Genomic language model for DRS chimera artifacts — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12923543/
