# How RECTIFY compares to related tools

RECTIFY is one of several open tools that touch RNA 3'-end correction,
poly(A) handling, or UMI consensus. The aim of this page is to help you pick
the right tool for the right task — not to argue that any one tool is
"better" in the abstract.

**Reading these tables.** A ✓ marks a feature whose presence we have
verified in the tool's published documentation (cited below). An empty cell
means we could not find that feature in the public docs as of 2026-05; it
does not necessarily mean the tool cannot do it via scripting. Where a tool
takes an opinionated stance that differs from RECTIFY's, the cell explains
the difference rather than scoring a win/loss.

If you maintain one of the tools below and a cell is wrong, please open an
[issue](https://github.com/k-roy/RECTIFY/issues) — we will update with a
direct citation to your README or paper.

---

## Bucket A — Direct RNA-seq (Nanopore DRS) 3'-end & poly(A) tools

These tools all operate on Nanopore DRS data (poly(A) is sequenced as part
of the read) and have something to say about 3' end accuracy or poly(A)
tail length.

| Capability | **RECTIFY** | [nanopolish polya][nanopolish] | [tailfindr][tailfindr] | [dorado polya][dorado] | [Nano3P-seq][nano3p] |
|:-----------|:-----------:|:------------------------------:|:----------------------:|:----------------------:|:--------------------:|
| Poly(A) tail length estimation | ✓ (model + dorado `pt:i`) | ✓ | ✓ | ✓ | ✓ |
| 3'-end **position** correction (walk-back) | ✓ | | | | |
| Multi-aligner consensus | ✓ (minimap2 + mapPacBio + gapmm2 ± uLTRA/deSALT) | n/a (works on minimap2 BAMs) | n/a (signal-level) | n/a (basecall-time) | n/a |
| False junction handling for poly(A) on genomic A-tracts | ✓ | | | | |
| 5'-splice-junction rescue (terminal-exon truncation) | ✓ | | | | |
| Adaptive valley-based CPA clustering | ✓ | | | | |
| DESeq2 at gene + cluster resolution | ✓ | | | | |
| Bundled S. cerevisiae reference | ✓ | | | | |
| Open source license | MIT | MIT | MIT | BSD-3 | MIT |
| Active maintenance as of 2026-05 | ✓ | (paper 2019) | ✓ | ✓ (ONT-supported) | ✓ |

**When to pick something other than RECTIFY:**

- For poly(A) tail length only — no 3'-end-position concern — `dorado
  polya` is built into the basecaller and is the lowest-friction option;
  its `pt:i` tag is consumed by RECTIFY transparently when present.
- For raw-signal tail length on non-Dorado basecallers, `tailfindr` or
  `nanopolish polya` are appropriate.
- For full Nano3P-seq protocol-specific analyses (including 3'-end
  adenylation profiling), the Nano3P-seq pipeline ships protocol-aware
  pre-processing not present in RECTIFY.

[nanopolish]: https://github.com/jts/nanopolish "Loman lab, jts/nanopolish"
[tailfindr]: https://github.com/adnaniazi/tailfindr "Krause lab, adnaniazi/tailfindr"
[dorado]: https://github.com/nanoporetech/dorado "Oxford Nanopore Technologies, nanoporetech/dorado"
[nano3p]: https://github.com/novoalab/Nano3P_seq "Novoa lab, novoalab/Nano3P_seq"

---

## Bucket B — Short-read 3'-end (oligo-dT primed) tools

These tools target short-read 3'-end protocols (QuantSeq REV, 3' Tag-Seq,
PAS-seq, A-seq, 3'READS). The shared challenge is internal-priming: the
oligo-dT primer can anneal to genomic A-tracts and produce reads that look
like 3'-end fragments but represent mid-gene sequence.

| Capability | **RECTIFY** | [polyAsite/APAtrap][apatrap] | [QuantSeq Lexogen analyzer][lexogen] | [DaPars2][dapars] | [QAPA][qapa] |
|:-----------|:-----------:|:----------------------------:|:------------------------------------:|:-----------------:|:------------:|
| Internal-priming filtering (genomic A-tract downstream) | ✓ (walk-back + NET-seq deconv) | ✓ (A-tract motif filter) | ✓ (Lexogen-proprietary) | | n/a (uses annotated polyA sites) |
| 3'-end **position** correction | ✓ (walk-back) | | | | n/a |
| Differential 3' UTR usage (DESeq2 at cluster level) | ✓ | | | ✓ (DaPars own model) | ✓ (per-isoform PAU) |
| NET-seq-informed A-tract refinement | ✓ | | | | |
| Bundled S. cerevisiae reference | ✓ | | | | |
| Works on QuantSeq REV (antisense) | ✓ (`rectify correct --dT-primed-cDNA`) | partial (motif filter is strand-agnostic) | ✓ | partial | n/a |
| Open source license | MIT | GPL-2 | proprietary | GPL-3 | GPL-3 |

**When to pick something other than RECTIFY:**

- For an annotated-poly(A)-site framework (no de novo CPA calling, focus on
  pre-annotated APA), QAPA is the canonical tool.
- For human-focused APA analysis with established statistical pipelines,
  DaPars2 has a large installed base.
- The Lexogen analyzer integrates UMI counting natively for QuantSeq
  libraries; RECTIFY does not consume Lexogen UMIs (use `umi_tools dedup`
  before `rectify correct`).

[apatrap]: https://github.com/SCBaylor/APAtrap "Wang lab, SCBaylor/APAtrap"
[lexogen]: https://www.lexogen.com/quantseq-3mrna-sequencing/ "Lexogen QuantSeq"
[dapars]: https://github.com/3UTR/DaPars2 "Xia lab, 3UTR/DaPars2"
[qapa]: https://github.com/morrislab/qapa "Morris lab, morrislab/qapa"

---

## Bucket C — Poly(A)-aware / UMI consensus pipelines

These tools focus on per-read poly(A) handling or UMI-based consensus for
long-read sequencing. RECTIFY's `rectify correct-cdna` subcommand operates
in this space when given ONT PCR-cDNA data with PCB114 or similar UMI
adapters.

| Capability | **RECTIFY correct-cdna** | [pychopper][pychopper] | [isONcorrect][isoncorrect] | [scNanoGPS][scnanogps] | [oxbow / VSL][oxbow] |
|:-----------|:------------------------:|:----------------------:|:--------------------------:|:----------------------:|:--------------------:|
| ONT cDNA adapter classification | ✓ | ✓ | | ✓ | ✓ |
| UMI extraction + clustering | ✓ (PCB114) | | | ✓ | partial |
| POA consensus across UMI cluster | ✓ (when `[cdna-correct]` extras installed) | | ✓ (isONcorrect for splice-aware) | ✓ | |
| Strand orientation rescue | ✓ | ✓ | | | ✓ |
| Re-alignment after consensus | ✓ (mappy + edlib fuzzy anchors) | | | | |
| 3'-end correction integrated | ✓ | | | | |
| Open source license | MIT | MPL-2 | GPL-3 | MIT | Apache-2 |

**When to pick something other than RECTIFY:**

- For high-error-rate cDNA *isoform* consensus only (no UMI), `isONcorrect`
  + `isONclust` is the canonical option.
- For single-cell Nanopore (per-cell UMIs alongside reads), `scNanoGPS`
  handles cell barcode + UMI clustering in one step.
- For adapter trimming only (no correction), `pychopper` is lighter-weight
  and runs in seconds per file.

[pychopper]: https://github.com/epi2me-labs/pychopper "EPI2ME Labs, epi2me-labs/pychopper"
[isoncorrect]: https://github.com/ksahlin/isONcorrect "Sahlin lab, ksahlin/isONcorrect"
[scnanogps]: https://github.com/gaolabtools/scNanoGPS "Gao lab, gaolabtools/scNanoGPS"
[oxbow]: https://github.com/timplab/oxbow "Timp lab, timplab/oxbow"

---

## Verification policy

Every ✓ cell on this page should be backed by either:

1. The tool's published README or documentation page (cited via the
   footnote link in each table header), or
2. The tool's first journal paper if no online docs exist for the feature.

If a cell is wrong as of your version, please file an issue with a link to
the README / paper section that contradicts it. We will update with the
citation in the next docs release.
