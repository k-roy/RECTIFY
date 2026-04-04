# Installation

## Requirements

- Python 3.8 or later
- A Unix-like system (Linux, macOS)

## PyPI (recommended)

```bash
pip install rectify-rna
```

### With visualization support

```bash
pip install rectify-rna[visualize]
```

This adds `matplotlib` and `seaborn` for genome browser figures, metagene plots, and heatmaps.

## Conda

The conda install also brings in the [MEME Suite](https://meme-suite.org/) for motif discovery:

```bash
conda install -c conda-forge -c bioconda rectify-rna
```

## From source

```bash
git clone https://github.com/k-roy/RECTIFY
cd RECTIFY
pip install -e .            # Core install
pip install -e ".[dev]"     # Development install (adds pytest, black, mypy)
```

Run the test suite to verify the install:

```bash
pytest
```

---

## External aligners

`rectify align` and `rectify run` require one or more of the following aligners to be on your `PATH`:

| Aligner | Required? | Notes |
|---------|-----------|-------|
| [minimap2](https://github.com/lh3/minimap2) | Yes (default) | Long-read splice-aware aligner |
| [pbmm2](https://github.com/PacificBiosciences/pbmm2) | Recommended | PacBio RNA mode; improves junction accuracy |
| [gapmm2](https://github.com/vpc-ccg/gapmm2) | Recommended | Gap-aware minimap2 variant |
| [uLTRA](https://github.com/ksahlin/ultra) | Optional | Annotation-guided collinear chaining; best for small exons (11–20 nt) |
| [deSALT](https://github.com/ydLiu-HIT/deSALT) | Optional | De Bruijn graph mapper; adds a 5th aligner to consensus |

Install minimap2 via conda:

```bash
conda install -c bioconda minimap2 pbmm2
```

!!! note
    `rectify correct` (correction only, no alignment) does **not** require any external aligner — it works directly on pre-aligned BAM files.

---

## Core dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `pysam` | ≥ 0.19 | BAM file I/O |
| `numpy` | ≥ 1.20 | Numerical computing |
| `pandas` | ≥ 1.3 | Data manipulation |
| `pyBigWig` | ≥ 0.3.18 | BigWig export |
| `biopython` | ≥ 1.79 | Sequence utilities |
| `tqdm` | ≥ 4.60 | Progress bars |

---

## Verifying the install

```bash
rectify --version
# RECTIFY 2.7.6

rectify --help
# Usage: rectify [correct|run|align|analyze|export|extract|aggregate|netseq|validate|train-polya|batch] ...
```

For a quick end-to-end test using bundled yeast data, see the [Quick Start](quickstart.md) guide.
