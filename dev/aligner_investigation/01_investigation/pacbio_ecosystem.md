# The PacBio Long-Read RNA-seq Alignment Ecosystem (Iso-Seq / HiFi)

**Audience:** RECTIFY developers and agents.
**Purpose:** Map the PacBio-side alignment landscape for Iso-Seq (HiFi/CCS) transcriptomics
and contrast it algorithmically with the ONT (DRS/cDNA) ecosystem RECTIFY currently targets.
**Date:** 2026-06-18.

**Convention:** Each claim is tagged **[FACT]** (sourced from a primary repo/paper/man-page,
cited in References) or **[INFERENCE]** (analytical synthesis, including comparisons that draw
on RECTIFY's own `CLAUDE.md`). Numeric preset values are version-sensitive — always re-verify
against the installed tool's `--help` / `options.c`.

---

## 0. TL;DR

- PacBio's reference aligner is **pbmm2**, a C++ wrapper around the **minimap2 C API** plus
  PacBio BAM (pbbam) I/O. For Iso-Seq it runs `pbmm2 align --preset ISOSEQ --sort`, which is
  minimap2's spliced mode (the `splice` / `splice:hq` family) with HiFi-tuned penalties.
- The canonical Iso-Seq pipeline is `ccs → lima → isoseq refine → isoseq cluster2 →
  pbmm2 align (ISOSEQ) → isoseq collapse → SQANTI3/pigeon`. **Alignment sits between
  clustering and collapse.**
- **The single most important algorithmic difference for RECTIFY:** PacBio **trims the
  poly(A) tail BEFORE alignment** (at `isoseq refine`), whereas ONT DRS sequences the poly(A)
  tail *in-read* and aligns with it present. The "true 3' end" of a PacBio aligned read is a
  software trim boundary; the true 3' end of an ONT DRS read must be recovered from a noisy
  in-read A-run.
- **HiFi (Q20–Q33) reads have no systematic homopolymer-deletion bias**; ONT DRS does. The
  exact artifacts RECTIFY corrects (A-tract walk-back, poly-A indel correction, homopolymer
  soft-clip rescue, junction refinement) are ONT-specific and largely vanish on HiFi.

---

## 1. pbmm2 & minimap2 PacBio Presets

### 1.1 pbmm2 vs raw minimap2

**[FACT]** pbmm2 is "a minimap2 SMRT wrapper for PacBio data: native PacBio data in ⇨ native
PacBio BAM out." It links **minimap2's C API directly** (it does not shell out to the
`minimap2` binary) and adds pbbam for native PacBio BAM I/O.

Differences from raw minimap2:
- **[FACT] Native PacBio BAM I/O**: consumes `.subreads.bam`, `.ccs.bam`, `.hifi_reads.bam`,
  `.fasta/.fastq`, dataset `.xml`, and `.fofn`; preserves PacBio per-read BAM tags; emits a
  `.pbi` for dataset XML output.
- **[FACT] PacBio-tuned output defaults**: soft-clipping on (`-Y`), `X/=` CIGAR ops instead of
  `M` (`--eqx`), long-CIGAR `CG` tag handling.
- **[FACT] On-the-fly coordinate sorting** via `--sort` (integrated samtools-sort code; no
  external pipe). Sort threads default to 25% of alignment threads, capped at 8 (`-J`);
  per-thread sort memory via `-m`. Produces `.bai` by default. (Contrast RECTIFY's
  `minimap2 | samtools sort` pattern in `multi_aligner.py:run_minimap2`.)
- **[FACT] Reusable `.mmi` index** (`pbmm2 index ref.fa ref.mmi --preset ISOSEQ`). Constraint:
  once an index is built, **`-k`, `-w`, and `-u` cannot be overridden at align time** — they
  are baked into the index, so the index preset must match the align preset.

**[FACT] minimap2 version bundled:** pbmm2 ≥ v1.13.0 → minimap2 v2.26; pbmm2 < v1.13.0 →
minimap2 v2.15. Version history highlights: v1.0.0 (first stable, SMRT Link v7.0); v1.12.0
(CCS became default preset); v26.1.0 (calendar versioning; **ISOSEQ preset parameters updated**;
`AS` tag added).

### 1.2 pbmm2 presets

**[FACT]** Presets: `SUBREAD`, `CCS` (alias `HIFI`; default since v1.12.0), `ISOSEQ`,
`UNROLLED`. Homopolymer-compressed (HPC) minimizers are **always on for SUBREAD and UNROLLED**
and **off for CCS/HIFI and ISOSEQ** (pbmm2 passes `-u` to disable HPC for the accurate presets).

**[FACT, version-dependent]** Numeric strings from pbmm2 v1.x `pbmm2 align --help`
(pbmm2 splits minimap2's comma `-O q,q2` / `-E e,e2` into separate single-letter flags):

| Preset | Parameters (pbmm2 v1.x) |
|---|---|
| **SUBREAD** | `-k19 -w10 -o5 -O56 -e4 -E1 -A2 -B5 -z400 -Z50 -r2000 -L0.5 -g5000` (HPC on) |
| **CCS / HIFI** (default) | `-k19 -w10 -u -o5 -O56 -e4 -E1 -A2 -B5 -z400 -Z50 -r2000 -L0.5 -g5000` |
| **ISOSEQ** | `-k15 -w5 -u -o2 -O32 -e1 -E0 -A1 -B2 -z200 -Z100 -r200000 -L0.5 -g2000 -C5 -G200000` |
| **UNROLLED** | `-k15 -w15 -o2 -O32 -e1 -E0 -A1 -B2 -z200 -Z100 -r2000 -L0.5 -g10000` (HPC on) |

- `-o/-O` = two-piece gap-open 1/2; `-e/-E` = gap-extend 1/2; `-A/-B` = match/mismatch;
  `-z/-Z` = Z-drop / inversion Z-drop; `-r` = chaining bandwidth; `-L` = `--lj-min-ratio`;
  `-g` = max chain gap; `-G` = max intron (ISOSEQ only); `-C` = non-canonical GT-AG cost
  (ISOSEQ only).
- **[FACT]** `HIFI` is an exact **alias of `CCS`** — same parameters.
- **[FACT]** ISOSEQ is the only spliced preset; it exposes `-G`, `-C`, and `--no-splice-flank`
  as user overrides. **[FACT]** v26.1.0 changed ISOSEQ's exact numbers — re-verify with
  `pbmm2 align --help` on the installed version.
- **[INFERENCE]** pbmm2 ISOSEQ ≈ minimap2 `splice:hq` adapted for Iso-Seq (k15/w5, `-C5`, large
  `-G`/`-r`); i.e. the PacBio analogue of `minimap2 -ax splice -uf -C5`.

### 1.3 minimap2 PacBio-relevant presets (exact, from `options.c` / man page)

**[FACT]**
- **`map-pb`** (PacBio CLR genomic): `-Hk19` — **HPC minimizers ON**, k=19, w=10, global default
  match/gap. HPC compresses homopolymer runs to tolerate CLR's ~15% homopolymer indel error.
- **`map-hifi` / `map-ccs`** (HiFi/CCS genomic): expands `lr:hq` (k19, w19, `-g10000`,
  `-U50,500`) then overrides `-A1 -B4 -O6,26 -E2,1 -s200`. **No HPC** — HiFi resolves
  homopolymers accurately, so compression would only lose information.
- **`map-ont`** (Nanopore genomic): default init, non-HPC minimizers, k=15, w=10.
- **`splice`** (long-read spliced, strand unknown): `-k15 -w5 --splice -g2000 -G200000 -A1 -B2
  -O2,32 -E1,0 -C9 -z200 -ub --junc-bonus=9 --splice-flank=yes`.
- **`splice:hq`** (PacBio CCS/HiFi spliced — the Iso-Seq aligner preset): start from `splice`,
  then override `-C5 -O6,24 -B4`.

**[FACT]** `splice` vs `splice:hq` side-by-side:

| Param | `splice` (ONT-style) | `splice:hq` (HiFi/CCS) |
|---|---|---|
| `-B` (mismatch) | 2 | **4** |
| `-O` (gap open) | 2,32 | **6,24** |
| `-C` (non-canonical) | 9 | **5** |
| `-k -w -A -E -G -z --junc-bonus --splice-flank` | 15 / 5 / 1 / 1,0 / 200k / 200 / 9 / yes | inherited |

**[INFERENCE]** `splice:hq` *raises* mismatch and gap-open penalties (trust the accurate read;
a mismatch is more likely a real variant than an error) and *lowers* the non-canonical splice
cost (a non-canonical junction in a clean read is more believable). This is the **opposite
tuning direction** from RECTIFY's ONT settings.

### 1.4 `--splice-flank`

**[FACT]** `--splice-flank=yes` (default with `--splice`, inherited by `splice:hq`) models one
*extra* conserved base around a canonical junction: the base after a GT donor tends to be A/G
(~91–92% in human/mouse) and the base before an AG acceptor tends to be C/T — i.e. it prefers
`GT[A/G]..[C/T]AG`. minimap2's cookbook recommends `--splice-flank=no` for SIRV spike-ins.

**[INFERENCE]** The extra-base flank model can nudge a junction boundary to satisfy the A/G·C/T
preference. For 3'-end / terminal-exon precision this is undesirable, which is why RECTIFY sets
`--splice-flank=no`. PacBio Iso-Seq leaves it on (default), accepting the canonical bias because
HiFi accuracy makes the nudges reliable.

### 1.5 What RECTIFY actually runs (codebase-grounded)

**[FACT]** `rectify/core/align/multi_aligner.py::run_minimap2` builds:
```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD -y -t N \
    [--junc-bed <annotation.junc.bed> --junc-bonus 9] <genome> <reads>
```
This is the **base `splice` preset (ONT-tuned)** — B=2, O=2,32, C=9 — with three deliberate
departures: `-k14` (below splice's default k15, for noisy nanopore), `-G 5000` (yeast introns
< 1 kb, vs splice's 200k), and `--splice-flank=no`. RECTIFY pipes to `samtools sort -n`
(name-sorted) for cross-aligner consensus streaming. **It is the opposite tuning direction from
`splice:hq`** and has no PacBio/HiFi code path today.

---

## 2. Iso-Seq3 Pipeline & Where Alignment Fits

**[FACT]** The modern PacBio Iso-Seq stack (binary renamed `isoseq3 → isoseq` at v4.0;
`cluster2` added in v4.0; current ~v4.3.0):

| # | Stage | Command (canonical) | What it does |
|---|---|---|---|
| 1 | **CCS / HiFi** | `ccs subreads.bam ccs.bam --min-rq 0.9` | Collapse each ZMW's subreads into one HiFi consensus read (predicted ≥ Q20). Revio/Sequel II emit `hifi_reads.bam` directly. `ccs` ≥ v4.2.0 recommended. |
| 2 | **lima** | `lima hifi_reads.bam primers.fasta fl.bam --isoseq --peek-guess` | Primer/barcode demux; removes 5'/3' primers; orients reads 5'→3'. Output = **FL (full-length)** reads. Does NOT trim poly(A) or detect concatemers. |
| 3 | **isoseq refine** | `isoseq refine fl.bam primers.fasta flnc.bam --require-polya` | **Trims poly(A) tail** AND **removes concatemers** (5'-3'-5' artifacts). `--require-polya` keeps reads with ≥ `--min-polya-length` (default **20 bp**) poly(A), then **strips the tail**. Output = **FLNC**. |
| 4 | **isoseq cluster2** | `isoseq cluster2 flnc.fofn clustered.bam` | Hierarchical n·log(n) isoform clustering + QV-guided POA consensus. (Legacy `isoseq cluster` → `clustered.hq/lq.fasta.gz`.) Default outputs isoforms with ≥ 2 FLNC reads; `--singletons` to include 1-read. |
| 5 | **pbmm2 align** | `pbmm2 align --preset ISOSEQ --sort ref.mmi <in.bam> mapped.bam` | **GENOME ALIGNMENT** (minimap2 spliced, HiFi-tuned). Produces the sorted BAM collapse needs. |
| 6 | **isoseq collapse** | `isoseq collapse [--do-not-collapse-extra-5exons] mapped.bam [flnc.bam] collapsed.gff` | Collapse redundant isoforms by intron-chain (exonic) structure. Requires the **sorted, genome-aligned BAM** from step 5. |
| 7 | **pigeon / SQANTI3** | `pigeon prepare → classify → filter` | SQANTI-style structural classification + artifact filtering (downstream QC). |

**[FACT] Where alignment sits — the exact order is:**
```
refine → cluster2 → pbmm2 align (--preset ISOSEQ --sort) → collapse
```
Alignment is a **distinct middle step between clustering and collapse**. `isoseq collapse`
operates on genomic exon/intron structure and therefore **cannot run on unaligned reads** —
it requires the coordinate-sorted aligned BAM produced by `pbmm2 align --sort`.

**[FACT] FLNC = Full-Length Non-Concatemer** (also rendered "Non-Chimeric"). An FLNC read has:
(1) a valid 5' primer, (2) a valid 3' primer (both found+removed by lima), (3) a poly(A) tail
(detected+removed by refine; required under `--require-polya`), and (4) no concatemer (removed
by refine). lima establishes full-length status; refine establishes non-concatemer + poly-A.

**[FACT] Historical context:** Pre-isoseq3 clustering used **ICE (Iterative Clustering and
Error correction)** in the `pbtranscript` suite, which relied on NP-hard clique-finding;
isoseq3 replaced it with the hierarchical n·log(n) approach. cDNA_Cupcake's
`collapse_isoforms_by_sam.py` is the direct ancestor of `isoseq collapse`.

---

## 3. Poly-A Handling: PacBio vs ONT

This is the **central architectural divergence** for a 3'-end pipeline.

### 3.1 PacBio: poly-A trimmed BEFORE alignment

**[FACT]** On the Iso-Seq path the poly(A) tail is **physically removed from the read at
`isoseq refine`** (`--require-polya`, `--min-polya-length 20`), which runs **before** any genome
alignment (`pbmm2`). By the time reads reach the aligner, the poly(A) tail is gone. Refine
detects a ≥20 bp poly(A) homopolymer downstream of the validated 3' primer and strips it.
(The exact internal detector is not published; the historical standalone `trim_isoseq_polyA`
used an HMM — **[INFERENCE]** for refine's exact algorithm.)

**[INFERENCE]** Consequence: in an aligned Iso-Seq read, the **reported 3' end = the poly-A trim
boundary**, defined by software, not by in-read homopolymer basecalling. CPA-position confidence
depends on **trim-tool boundary accuracy**, and HiFi's high base accuracy makes the body/tail
junction sharp and easy to locate. The **poly-A length is discarded** unless separately measured.

### 3.2 ONT DRS: poly-A in-read at alignment time

**[FACT, RECTIFY CLAUDE.md + DRS literature]** ONT direct RNA-seq sequences the entire molecule
from adapter through the poly(A) tail into the RNA body — the tail **is in the read** when
aligned. The 3' end must be inferred from where the in-read poly-A run ends and genomic sequence
begins, and that boundary is exactly where ONT's systematic homopolymer deletion bias corrupts
the A-run length, sliding the apparent 3' end.

**[FACT]** RECTIFY's `--drs` workflow mirrors the PacBio representation deliberately: **Step 0**
(`trim_drs_bam_polya`) trims poly-A pre-alignment, and **Step 4** (`restore_polya_softclips`)
restores it as soft-clips post-correction, keeping tail-length metadata in a parquet sidecar.
In effect RECTIFY *manufactures* the PacBio-style trim-boundary representation for ONT, while
preserving the tail-length information PacBio throws away.

### 3.3 The trade-off

**[INFERENCE]** ONT DRS retains poly-A length (valuable for APA / poly-A-length biology) at the
cost of 3'-end positional precision; HiFi gives precise 3'-end position at the cost of tail
length. RECTIFY's protocol flags encode this: `--dT-primed-cDNA` (QuantSeq, tail NOT in read,
enables AG-mispriming) vs default DRS (tail in read). PacBio Iso-Seq is structurally closest to
the "tail already removed" regime.

---

## 4. HiFi Accuracy Implications for Splice & 3'-End

### 4.1 Base accuracy

**[FACT]**
- **PacBio HiFi/CCS**: defined as predicted ≥ **Q20 (<1% error)**; ~4 passes ≈ Q20, ≥9 passes ≈
  **Q30 (99.9%)**. HG002 benchmark Q27 (99.8%); Revio markets median **Q30–Q33 (up to 99.95%)**.
- **ONT**: R10.4.1 DNA > 99% single-read (one 16S study: raw Q18.8 → Q22.4 post-primer); duplex
  ~Q30. **DRS** historically far worse: **RNA001/002 median 87–92% identity** with deletions
  dominating; **RNA004 (2024)** median alignment identity **98.67%**, with median indel
  proportion falling **7.19% (RNA002) → 0.88% (RNA004)** (~8× reduction).

**[INFERENCE]** RECTIFY's design assumptions (severe poly-A/homopolymer indels) match RNA002-era
DRS. On RNA004 the indel *magnitude* shrinks ~8× but the *character* (systematic homopolymer
deletions) persists. HiFi sidesteps it entirely.

### 4.2 Homopolymers (the load-bearing difference)

**[FACT]** Mechanism: PacBio errors are **stochastic** (random per-pass), so CCS averages them
out → near-zero *systematic* homopolymer bias. ONT errors are **systematic**: once a homopolymer
exceeds the pore's sensing span, additional identical bases don't change the ionic current, so
length is inferred from dwell time (std dev > median) → undercalling.

**[FACT]** Quantitative: homopolymers ≥5 bp correctly called ≈ **PacBio HiFi 99%, raw ONT ~91%,
PacBio CLR ~88%**. HiFi: >99.5% of homopolymers up to 5 bp accurate. ONT DRS: poly-A/poly-T
deletion errors exceed 2% on R9 (poly-A is the *best* basecalled homopolymer yet still the
dominant error source by abundance).

**[INFERENCE]** The exact artifact RECTIFY corrects — systematic deletion undercalling inside
A-tracts/poly-A that slides the reported 3' end — is **the signature ONT homopolymer error and
is largely absent in HiFi**. This is the core platform asymmetry justifying RECTIFY's
homopolymer-aware modules (2C indel correction, 2E A-tract walk-back, the AT/CG-split empirical
penalty tables). Those tables are explicitly "S. cerevisiae R10.4.1-specific" and **would not
transfer to HiFi**.

### 4.3 Splice junctions

**[FACT]** LRGASP (Nature Methods 2024, >427M long reads, PacBio + ONT): "Libraries with longer,
more accurate sequences (cDNA-PacBio and R2C2-ONT) produced more precise transcripts," while
higher-depth ONT libraries improved quantification. Junction support by Illumina + canonical
splice sites was very high except for some ONT pipelines; **LyRic on PacBio showed the highest
% of Illumina-supported junctions**. ONT reads' ~5–10% error rate makes mapping near splice
sites hard and generates false novel junctions.

**[INFERENCE]** HiFi's near-perfect base accuracy yields cleaner junction detection and fewer
spurious novel junctions, reducing the need for junction-refinement machinery. ONT noise is
precisely why RECTIFY implements **Module 2H** (`junction_refiner.py`), the junction-overhang
chimera filter, and the canonical-HP-prior tie-break — ONT-specific countermeasures a HiFi
pipeline would largely not need.

### 4.4 3'-end / CPA concordance

**[FACT]** PolyA_DB v4 integrates 135 PacBio Iso-Seq + 92 ONT samples (~1.4M PAS); ~20% of the
full PAS set matches long-read transcript-end sites. **[GAP]** No clean nucleotide-level
PacBio-vs-ONT 3'-end positional-concordance percentage appears to be published in citable form.
The strongest defensible statement: long-read 3'-ends recover ~20% of curated PAS; ONT DRS
3'-end error is dominated by the >2% poly-A deletion rate (R9) / 0.88% indel rate (RNA004),
whereas HiFi 3'-end error is bounded by its overall Q27–Q33 accuracy and poly-A trim-tool
precision.

---

## 5. Other PacBio Aligners & Downstream Tools

| Tool | Role | Algorithm | PacBio status | 3'/5' relevance |
|---|---|---|---|---|
| **pbmm2 / minimap2** | Spliced aligner | Minimizer chaining + DP | **Current standard** (`--preset ISOSEQ` / `splice:hq`) | Junction + 3'-end via HiFi-tuned penalties |
| **GMAP** | Spliced aligner | Oligomer chaining + sandwich DP | Legacy standard (pre-minimap2); still offered | De novo splice precision, microexons |
| **STARlong** | Spliced aligner | Seed-extend (MMP) + stitching | ENCODE legacy | Dense seeding, `--scoreGapNoncan` |
| **deSALT** | Spliced aligner | RdBG-index + 2-pass SDP | Supported (`-x ccs`/`clr`) | Error-model presets per platform |
| **uLTRA** | Spliced aligner | Annotation + MEM collinear chaining | Supported (`--isoseq`) | Small-exon (<30 nt) precision |
| **cDNA_Cupcake** | Collapse | SJ-identity collapse | **Deprecated** → `isoseq collapse` | `--max_3_diff` end handling |
| **TAMA** | Collapse/merge | Threshold-based end/SJ collapse | Active alternative | Tunable 5'/3'/SJ wobble (`-a`/`-z`/`-m`) |
| **SQANTI3** | QC/classify | Structural categories vs ref | **Current standard** | Intra-priming flag, CAGE/polyA support |

### 5.1 GMAP

**[FACT]** Wu & Watanabe, *Bioinformatics* 2005. Algorithm: minimal-sampling oligomer lookup →
**oligomer chaining** (diagonalization) → **sandwich DP** (bidirectional DP that solves for an
intron position bounded between two anchored exon segments → precise canonical GT-AG placement)
→ microexon identification. Iso-Seq usage:
`gmap -d <db> -f samse -n 0 --no-chimeras -t N input.fasta > aln.sam`. Key options: `-f samse`,
`-n` (npaths, default 5), `--no-chimeras` (≡ `--chimera-margin=0`), `--cross-species`,
`-z` (cDNA direction), `-K` (max intron). PacBio's old ToFU/Cupcake used GMAP before
`collapse_isoforms_by_sam.py`. **[INFERENCE]** Slower than minimap2 (nucleotide-level sandwich
DP vs minimizer chaining + SIMD); minimap2 (2018) displaced it on speed.

### 5.2 STARlong

**[FACT]** Long-read build of STAR (Dobin 2013). ENCODE PacBio Iso-Seq protocol used STAR
v2.5.2a, GRCh38, GENCODE v24. Recommended params: `--seedSearchStartLmax 50` (seed start every
50 bases — denser seeding for long reads), `--seedPerReadNmax 100000`,
`--alignTranscriptsPerReadNmax 100000`, `--outFilterMismatchNmax 2000`,
`--scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8`, `--alignEndsType Local`. Known
issues with very long reads (max-length ceiling; BAM indexing). **[INFERENCE]** Largely ENCODE
legacy now; pbmm2/minimap2 dominate.

### 5.3 SQANTI3

**[FACT]** Conesa Lab; Pardo-Palacios et al., *Nature Methods* 2024. Runs **downstream of
alignment + collapse** on isoform GTFs. Structural categories: **FSM** (Full Splice Match — all
internal junctions match a reference; 5'/3' ends may differ), **ISM** (Incomplete SM — subset of
consecutive ref junctions; typically 5'/3' truncation), **NIC** (Novel In Catalog — known
donors/acceptors in a new combination), **NNC** (Novel Not in Catalog — ≥1 novel splice site),
plus Antisense / Genic-intron / Genic-genomic / Intergenic / Fusion.
3'/5' orthogonal support: CAGE peaks (5' TSS, `dist_peak`/`within_peak`), TSS ratio,
**polyA_motif** (top hexamer within 50 bp upstream of 3' end) + PolyASite/Quant-seq support.
**Intra-priming flag**: A-stretch downstream of the TTS (≥60% A's downstream) → likely false 3'
end. SQANTI3 filter (rules + ML) removes intra-priming / RT-switching / low-coverage /
non-canonical artifacts.
**[INFERENCE]** SQANTI3's intra-priming detection is conceptually the same genomic-downstream-A
check as RECTIFY's AG-mispriming / A-tract modules — but SQANTI3 *flags/filters whole isoforms
post hoc*, whereas RECTIFY *corrects per-read 3' coordinates*.

### 5.4 TAMA

**[FACT]** GenomeRIK; Kuo et al., *BMC Genomics* 2020. TAMA Collapse thresholds: `-a` (5',
default 10 bp), `-z` (3', default 10 bp), `-m` (splice junction, default 10 bp); `-x capped`
(cap-selected, strict 5') vs `no_cap` (standard Iso-Seq with 5'-degraded reads). Modes: TSSC
(merge sharing 3' ends + junctions, allow 5' variation) and ECC. **[FACT]** TAMA explicitly
**preserves more 5'/3' end diversity** than Cupcake (treats end differences > threshold as
distinct isoforms) → more transcript models; its benchmark reports lower 5'/3' wobble than
Cupcake. TAMA Merge unifies annotations across samples/tools.

### 5.5 cDNA_Cupcake

**[FACT]** Magdoll/cDNA_Cupcake (E. Tseng). Core: `collapse_isoforms_by_sam.py` (collapse by
SJ identity; `--max_5_diff` 1000, `--max_3_diff` 100, `--max_fuzzy_junction` 5),
`get_abundance_post_collapse.py`, `chain_samples.py`, `fusion_finder.py`. **[FACT] Deprecated
(~2023); collapse functionality folded into official `isoseq collapse`.** New analyses should
use `isoseq collapse` + pigeon.

### 5.6 deSALT & uLTRA on PacBio

**[FACT] deSALT** (Liu et al., *Genome Biology* 2019): RdBG-index (de Bruijn graph, default
k=22) + 2-pass skeleton alignment. **`-x` preset accepts `ccs` (HiFi/Iso-Seq), `clr`, `ont1d`,
`ont2d`**; `-T` restricts splice detection to the forward transcript strand.
**[FACT] uLTRA** (Sahlin & Mäkinen, *Bioinformatics* 2021): annotation-guided MEM-seed +
collinear chaining; includes all annotated exons of candidate genes (excels at exons <30 nt);
wraps minimap2 for out-of-annotation reads. **Presets `--isoseq` (PacBio) and `--ont`.**
**[FACT]** So both Tier-2 aligners in RECTIFY's panel **already expose PacBio presets** — but
RECTIFY does not pass them (see §6).

---

## 6. Porting RECTIFY to PacBio HiFi Iso-Seq

RECTIFY is ONT-first throughout. Supporting PacBio HiFi Iso-Seq is **not a parameter tweak** —
it changes the alignment presets, removes most of the correction modules' reason to exist, and
shifts where the 3' end comes from. Codebase-grounded change list:

### 6.1 Alignment presets (multi_aligner.py)

- **[FACT/INFERENCE] minimap2**: switch `run_minimap2` from base `-ax splice ... -k14
  --splice-flank=no` to **`-ax splice:hq`** (adds `-C5 -O6,24 -B4`), drop the `-k14` override
  (HiFi wants k15), and raise `-G` if not yeast. Keep `--secondary=no --MD`. Reconsider `-uf`:
  Iso-Seq cDNA is stranded after lima orientation, so `-uf` may still be valid, but verify.
  `--splice-flank` can revert to default (yes) since HiFi makes the canonical nudge reliable.
- **[INFERENCE] mapPacBio (BBMap)**: BBMap's `mapPacBio.sh` is already PacBio-oriented; the
  RECTIFY-specific `_pt:i:N` QNAME stripping and DRS dedup logic are ONT/Dorado artifacts that
  do not apply to Iso-Seq FLNC reads (bare read names). The `>6 kb` read splitting may still be
  needed for full-length transcripts.
- **[FACT] deSALT**: `run_desalt` currently calls `aln` with **no `-x` preset**. For HiFi, add
  **`-x ccs`** (and consider `-T` for stranded). The `-G` annotation flag is already disabled
  (SIGSEGV on yeast GTF) — unchanged.
- **[FACT] uLTRA**: `run_ultra` passes no platform flag. Add **`--isoseq`** instead of the ONT
  default.
- **[INFERENCE] gapmm2**: a minimap2 wrapper; its terminal-exon homopolymer refinement targets
  ONT homopolymer error and is largely unnecessary for HiFi. Either pass an `:hq`-equivalent
  preset or drop gapmm2 from the HiFi panel.

### 6.2 Poly-A handling (the big one)

- **[FACT/INFERENCE]** PacBio poly-A is trimmed at `isoseq refine` **pre-alignment**. So for
  Iso-Seq input, RECTIFY's poly-A-tail-dependent modules should be **disabled by default**,
  analogous to how `--short-read` already disables poly(A) trimming (2B) and indel correction
  (2C). The aligned 3' end *is* the CPA candidate directly — no in-read walk-through needed.
- A new protocol mode (e.g. `--isoseq` / `--pacbio-hifi`) is the clean abstraction, parallel to
  the existing `--drs` / `--dT-primed-cDNA` / `--short-read` switches. It would:
  - skip Step 0 DRS trim and Step 4 soft-clip restore (no in-read tail);
  - disable Module 2B (poly-A trim) and 2C (indel correction);
  - disable Module 2E (A-tract walk-back) and 2G (homopolymer soft-clip rescue) — these correct
    ONT homopolymer slippage that HiFi does not exhibit;
  - **keep** AG-mispriming *only if* the upstream library is oligo-dT-primed cDNA with internal
    priming risk (Iso-Seq uses a 3' primer + poly-A selection, so intra-priming is the SQANTI3
    concern; evaluate whether RECTIFY's genomic-downstream-A check adds value here);
  - **keep** junction handling but expect far fewer corrections.

### 6.3 Error model / penalty tables

- **[FACT]** The empirical penalty tables (`penalty_scores.tsv`, `str_penalty_scores.tsv`) and
  the `_CANONICAL_HP_PRIOR = 0.5` calibration are **"S. cerevisiae R10.4.1-specific"** (CLAUDE.md)
  — calibrated to ONT homopolymer deletion rates. **They must NOT be reused for HiFi.** HiFi
  would need either no penalty table (defaults) or a HiFi-recalibrated table (deletion penalties
  much higher because HiFi indels are rare and meaningful — mirroring `splice:hq`'s `-O6,24`).
- **[INFERENCE]** Module 2H junction refinement and the overhang chimera filter remain
  conceptually valid but will fire on far fewer reads; their HP-aware edit-distance scoring is
  tuned to ONT noise and may need re-tuning or simple disabling for HiFi.

### 6.4 Input plumbing

- **[INFERENCE]** PacBio inputs are FLNC/clustered BAMs (already CCS-consensus, primers + poly-A
  removed). RECTIFY's FASTQ-chunking mandate (`rectify split`) still applies for I/O scaling on
  HPC, but the DRS-specific Dorado handling (`samtools fastq -T pt`, duplicate-UUID dedup in
  `run_gapmm2`) is ONT-only and should be bypassed for Iso-Seq inputs.

### 6.5 Net assessment

**[INFERENCE]** Most of RECTIFY's value is *correcting ONT-specific homopolymer/poly-A
artifacts*. On HiFi those artifacts are largely absent and the poly-A is already trimmed
upstream, so a HiFi port is mostly about (a) flipping aligner presets to the `:hq`/`ccs`/`isoseq`
family, (b) adding an `--isoseq` protocol mode that disables the homopolymer/poly-A correction
modules, and (c) refusing to apply the ONT-calibrated penalty tables. The pipeline's distinctive
correction logic would run, but mostly as a near-no-op — which is itself a useful validation
signal (a HiFi run should show very few corrections). The clearer long-term value on the PacBio
side is **3'-end clustering / CPA quantification on the already-clean aligned ends**, not
artifact correction.

---

## 7. References

**pbmm2 / minimap2**
- PacificBiosciences/pbmm2 — README & changelog. https://github.com/PacificBiosciences/pbmm2
- lh3/minimap2 — README, `options.c`, `minimap2.1` man page, `cookbook.md`.
  https://github.com/lh3/minimap2
- minimap2 Issue #769 — mapping pbmm2 presets to minimap2 (HPC k-mer). https://github.com/lh3/minimap2/issues/769
- Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 2018;34:3094–3100. https://doi.org/10.1093/bioinformatics/bty191
- Snakemake-wrappers pbmm2 (preset/override corroboration). https://snakemake-wrappers.readthedocs.io/

**Iso-Seq pipeline**
- PacificBiosciences/IsoSeq — clustering/classification/collapse docs. https://github.com/PacificBiosciences/IsoSeq
- Iso-Seq documentation. https://isoseq.how/
- Magdoll/cDNA_Cupcake (Cupcake → isoseq collapse; ICE history). https://github.com/Magdoll/cDNA_Cupcake
- bowhan/trim_isoseq_polyA (historical HMM poly-A trimmer). https://github.com/bowhan/trim_isoseq_polyA
- nf-core/isoseq (FLNC definition, aligner menu). https://github.com/nf-core/isoseq / PMC10199315

**HiFi vs ONT accuracy / homopolymers / junctions / 3'-end**
- LRGASP Consortium. Systematic assessment of long-read RNA-seq methods. *Nature Methods* 2024. https://www.nature.com/articles/s41592-024-02298-3 / bioRxiv 2023.07.25.550582
- Hwang et al. Sequencing accuracy and systematic errors of nanopore direct RNA sequencing.
  *BMC Genomics* 2024. https://doi.org/10.1186/s12864-024-10440-w / bioRxiv 2023.03.29.534691
- RNA004 direct RNA improvement. bioRxiv 2025.05.01.651717 / PMC12663090.
- R10.4.1 full-length 16S. *Appl Environ Microbiol* 2023. https://journals.asm.org/doi/10.1128/aem.00605-23
- PacBio HiFi "how it works" + Revio brochure. https://www.pacb.com/technology/hifi-sequencing/how-it-works/
- de novo mutation discovery with HiFi (HG002 Q27). *Genome Medicine* 2023. https://doi.org/10.1186/s13073-023-01183-6
- PacBio vs Nanopore error-rate review. https://www.cd-genomics.com/blog/pacbio-nanopore-error-rate-correction-strategies/
- PolyA_DB v4 (135 Iso-Seq + 92 ONT). PubMed 41316728 / PMC12807684.
- Benchmarking APA sequencing methods. *Genome Biology* 2021. https://doi.org/10.1186/s13059-021-02502-z

**Other aligners / QC**
- Wu TD, Watanabe CK. GMAP. *Bioinformatics* 2005;21:1859–1875. https://doi.org/10.1093/bioinformatics/bti310
- Dobin A et al. STAR. *Bioinformatics* 2013;29:15–21. https://doi.org/10.1093/bioinformatics/bts635 ; STARlong: https://github.com/alexdobin/STAR ; ENCODE Iso-Seq protocol (encodeproject.org)
- Tardaguila M et al. SQANTI. *Genome Res* 2018;28:396–411. https://doi.org/10.1101/gr.222976.117 ; SQANTI3 (Pardo-Palacios et al., *Nat Methods* 2024, s41592-024-02229-2): https://github.com/ConesaLab/SQANTI3
- Kuo RI et al. TAMA. *BMC Genomics* 2020;21:751. https://doi.org/10.1186/s12864-020-07123-7 ; https://github.com/GenomeRIK/tama
- Liu B et al. deSALT. *Genome Biology* 2019;20:274. https://doi.org/10.1186/s13059-019-1850-9 ; https://github.com/ydLiu-HIT/deSALT
- Sahlin K, Mäkinen V. uLTRA. *Bioinformatics* 2021;37:4643–4651. https://doi.org/10.1093/bioinformatics/btab540 ; https://github.com/ksahlin/uLTRA

**RECTIFY codebase (grounding for §1.5, §6)**
- `rectify/core/align/multi_aligner.py` — `run_minimap2`, `run_desalt`, `run_ultra`,
  `run_gapmm2`, `run_mapPacBio` (ONT presets, no PacBio path).
- `CLAUDE.md` — protocol flags, DRS Step 0/4, penalty-table calibration notes.

---

## 8. Confidence & Caveats

- **High confidence (FACT):** pipeline order (refine→cluster2→pbmm2→collapse), FLNC definition,
  `--require-polya` / `--min-polya-length 20`, minimap2 `splice` vs `splice:hq` exact deltas
  (`-C5 -O6,24 -B4`), deSALT `-x ccs` and uLTRA `--isoseq` presets, RECTIFY's own minimap2 args.
- **Version-sensitive:** pbmm2 numeric preset strings (ISOSEQ changed in v26.1.0 — verify with
  `pbmm2 align --help`); minimap2 preset numbers settled ~v2.27.
- **From search-snippet extraction (re-verify against source PDFs):** decimal accuracy figures
  (98.67% RNA004 identity, 7.19%→0.88% indel, Q27/Q33, >2% poly-A deletion).
- **Genuine gap:** no published nucleotide-level PacBio-vs-ONT 3'-end positional-concordance
  percentage.
- **Inference framing:** all "implications for RECTIFY" and the §6 porting plan are analytical
  synthesis grounded in the cited facts plus RECTIFY's `CLAUDE.md`, not direct quotes.
