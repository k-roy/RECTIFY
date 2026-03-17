# RECTIFY: Unified RNA 3' End Correction Framework

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Implementation Status](https://img.shields.io/badge/status-complete-brightgreen.svg)](https://github.com/k-roy/RECTIFY)
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

## Installation

### From source (development)

```bash
git clone https://github.com/k-roy/RECTIFY.git
cd RECTIFY
pip install -e .
```

### From PyPI (future release)

```bash
pip install rectify
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
