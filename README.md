# RECTIFY: Unified RNA 3' End Correction Framework

[![PyPI version](https://badge.fury.io/py/rectify-rna.svg)](https://pypi.org/project/rectify-rna/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-147%20passing-brightgreen.svg)](tests/)

**RECTIFY** (**R**NA 3' **E**nd **C**orrection **T**ool **I**ntegrating **F**alse-priming and pol**y**(A) ambiguity) is a unified framework for correcting 3' end mapping artifacts in poly(A)-tailed RNA sequencing data.

## Overview

RECTIFY addresses two fundamental problems affecting RNA 3' end mapping:

1. **A-tract Ambiguity (Universal)**: Genomic A-tracts near true 3' ends create positional uncertainty affecting ALL poly(A)-tailed RNA-seq technologies
2. **Technology-Specific Artifacts**:
   - **AG mispriming**: Internal priming on A/G-rich regions (oligo-dT methods)
   - **Poly(A) tail alignment**: Tail bases align to genomic A-tracts creating systematic shifts (when poly(A) is sequenced)

### Key Features

- **Modular correction strategies** that apply based on sequencing technology
- **Universal A-tract ambiguity detection** for all poly(A)-tailed RNA-seq
- **AG mispriming screening** (from original RECTIFY, Roy & Chanfreau 2019)
- **Poly(A) tail trimming and indel artifact correction** (for direct RNA-seq: nanopore, Helicos, QuantSeq)
- **NET-seq refinement** (optional, technology-independent)
- **Unified output format** with confidence scores and QC flags

## How It Works

RECTIFY corrects 3' end mapping artifacts by walking a read through multiple correction steps.
Here's a single Nanopore read traversing the full pipeline:

```
═══════════════════════════════════════════════════════════════════════════════
 STEP 1: RAW ALIGNMENT
═══════════════════════════════════════════════════════════════════════════════

 Genome:        ...CTAGTGACAGTCAAAAAAAA-AAACAAAAGT|CTAGCGATC...
                                        └───────┘
                                     genomic A-tract (11 bp)

 Nanopore read: ...CTAGTGACAGTCAAAAAAAATAAA-AAAAA--AAAAAAAAAAAAAAAAAAAAAA
                                       ↑       ↑↑  └──────────────────┘
                                   seq error  indels    soft-clipped
                                   (T in As)            poly(A) tail

 Aligner output:
   • Mapped 3' end: position 42 (at the |)
   • Soft-clip starts after last aligned base
   • Indels in A-tract are alignment artifacts (homopolymer confusion)
   • The 'T' in the read's A-tract is a sequencing error

═══════════════════════════════════════════════════════════════════════════════
 STEP 2: INDEL CORRECTION — Walk Backwards to True 3' End
═══════════════════════════════════════════════════════════════════════════════

 Starting from mapped 3' end, walk upstream until genome and read AGREE on
 the first non-A base:

 Genome:  ...CTAGTGACAGT C A A A A A A A A - A A A C A A A A G T | ...
 Read:    ...CTAGTGACAGT C A A A A A A A A T A A A - A A A A - - | ...
                         ↑                 ↑       ↑         ↑↑
                      AGREE             seq err   del       del
                      (C=C)             (ignore)  (+1bp)   (+2bp)

 RECTIFY logic:
   • Walk back from position 42
   • Skip A's (ambiguous with poly(A) tail)
   • Find first non-A agreement: 'C' at position 31
   • Count deletions in A-tract: 3 bp total
   • Corrected range: [31, 42] → ambiguity window = 11 bp

═══════════════════════════════════════════════════════════════════════════════
 STEP 3: NET-seq REFINEMENT (Optional)
═══════════════════════════════════════════════════════════════════════════════

 For species with NET-seq data (captures Pol II at nascent 3' ends), we can
 resolve the ambiguity window. Note: oligo(dT) priming spreads signal upstream
 due to internal A-tract priming.

 Position:           31      34      37      40      42
                     |       |       |       |       |
 Raw NET-seq:     ▂▄▆█▆▄▃▃▄▅▇█▇▅▄▃▃▄▅▆▇█▇▅▄▂▁▁▁▁▁▁▁▁
                  └──────┘ └──────┘ └──────┘
                   Peak 1   Peak 2   Peak 3
                  (spread)  (spread) (spread)

 After deconvolution (removing oligo(dT) spreading artifact):

 True CPA signal:    █           █           █
                     ↑           ↑           ↑
                   pos 31      pos 35      pos 39
                   (25%)       (50%)       (25%)

═══════════════════════════════════════════════════════════════════════════════
 STEP 4: PROPORTIONAL ASSIGNMENT
═══════════════════════════════════════════════════════════════════════════════

 The Nanopore read has ambiguity window [31, 42]. Using NET-seq deconvolved
 peaks, we assign probability to each CPA site:

 Nanopore read:    ════════════════════════════════
                   [          ambiguity           ]
                   31                             42

 NET-seq peaks:    █ (25%)      █ (50%)      █ (25%)
                   31           35           39

 Final output:
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  read_id    position   confidence   ambiguity_range   netseq_support     │
 │  ────────   ────────   ──────────   ───────────────   ──────────────     │
 │  read001    35         HIGH         11                0.50               │
 │  read001    31         MEDIUM       11                0.25               │
 │  read001    39         MEDIUM       11                0.25               │
 └──────────────────────────────────────────────────────────────────────────┘

 Without NET-seq: Reports ambiguity window [31, 42] with uniform distribution
```

## Installation

### From source (development)

```bash
git clone https://github.com/k-roy/RECTIFY.git
cd RECTIFY
pip install -e .
```

### From PyPI

```bash
pip install rectify-rna
```

## Quick Start

### QuantSeq (oligo-dT short-read)

```bash
rectify correct quantseq.bam \
  --genome sacCer3.fa \
  --annotation genes.gtf \
  --polya-sequenced \
  --output corrected_3ends.tsv
```

### Nanopore direct RNA-seq with NET-seq refinement

```bash
rectify correct nanopore.bam \
  --genome sacCer3.fa \
  --annotation genes.gtf \
  --polya-sequenced \
  --aligner minimap2 \
  --netseq-dir churchman_bigwigs/ \
  --output corrected_3ends.tsv
```

## Output Format

RECTIFY produces a TSV file with corrected 3' end positions and QC metrics:

```
read_id  chrom  strand  raw_position  corrected_position  ambiguity_min  ambiguity_max  ambiguity_range  correction_type  confidence  qc_flags
read001  chrI   +       147588        147585              147583         147588         5                polya_trim       high        PASS
read002  chrI   +       147593        147591              147591         147593         2                ag_mispriming    medium      AG_RICH
```

## Module Architecture

RECTIFY applies corrections modularly based on your data:

1. **Module 1: A-tract Ambiguity** (always applied)
   - Identifies genomic A-tracts near 3' ends
   - Calculates ambiguity windows

2. **Module 2A: AG Mispriming** (when oligo-dT priming used)
   - Screens for downstream AG-richness
   - Flags likely misprimed reads

3. **Module 2B+2C: Poly(A) Corrections** (when poly(A) IS sequenced)
   - Models and trims poly(A) tails
   - Detects and removes indel artifacts

4. **Module 3: NET-seq Refinement** (optional)
   - Resolves ambiguity using NET-seq data
   - Assigns confidence scores

## Citation

If you use RECTIFY, please cite:

**Original RECTIFY (AG mispriming correction):**
> Roy KR, Chanfreau GF. Robust mapping of polyadenylated and non-polyadenylated RNA 3' ends at nucleotide resolution by 3'-end sequencing. *Methods*. 2020 Apr 1;176:4-13. doi: 10.1016/j.ymeth.2019.05.016. [PMID: 31128237](https://pubmed.ncbi.nlm.nih.gov/31128237/)

**RECTIFY 2.0 (unified framework):**
> Manuscript in preparation

## License

MIT License - See [LICENSE](LICENSE) for details

## Contact

- Kevin R. Roy - [kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com)
- GitHub: [k-roy/RECTIFY](https://github.com/k-roy/RECTIFY)

## Acknowledgments

- Original RECTIFY development supported by Chanfreau Lab, UCLA
- NET-seq data from Churchman Lab, Harvard Medical School
