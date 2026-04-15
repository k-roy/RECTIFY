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

| Aligner | Required? | Install | Notes |
|---------|-----------|---------|-------|
| [minimap2](https://github.com/lh3/minimap2) | Yes (default) | conda/bioconda | Long-read splice-aware aligner |
| [mapPacBio](https://sourceforge.net/projects/bbmap/) (BBMap) | Recommended | conda/bioconda (`bbmap`) | PacBio RNA mode; improves junction accuracy |
| [gapmm2](https://github.com/vpc-ccg/gapmm2) | Recommended | pip or conda | Gap-aware minimap2 variant |
| [uLTRA](https://github.com/ksahlin/ultra) | Optional | pip or conda | Annotation-guided collinear chaining; best for small exons (11–20 nt) |
| [deSALT](https://github.com/ydLiu-HIT/deSALT) | Optional | **vendored** (Linux/x86_64) | De Bruijn graph mapper; bundled binary used automatically |

### Quickest install

```bash
# minimap2, mapPacBio (BBMap), gapmm2 via conda
conda install -c bioconda minimap2 bbmap
pip install gapmm2 ultra-bioinformatics    # pip-installable aligners

# deSALT: vendored Linux/x86_64 binary is bundled — no install needed.
# For other platforms, or to install system-wide:
rectify install-aligners --desalt
```

### rectify install-aligners

RECTIFY ships a built-in installer for all supported aligners:

```bash
rectify install-aligners --check          # show what's already available
rectify install-aligners --all            # install everything possible
rectify install-aligners --desalt         # deSALT only
rectify install-aligners --desalt --force # reinstall even if already on PATH
```

Binaries are installed to `~/.rectify/bin/` by default. Add that directory to
`PATH` once (e.g. in `~/.bashrc`):

```bash
export PATH="$HOME/.rectify/bin:$PATH"
```

### deSALT bundled binary

A pre-compiled deSALT v1.5.6 binary for Linux/x86_64 (glibc 2.6.32+) is bundled
with the package at `rectify/data/bin/linux_x86_64/deSALT`. RECTIFY uses it
automatically when `deSALT` is not found on `PATH` — no action required on
standard Linux HPC systems.

For macOS or Linux/aarch64, `rectify install-aligners --desalt` downloads the
source from GitHub and compiles it (requires `gcc` and `zlib`).

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
# RECTIFY 0.9.0

rectify --help
# Usage: rectify [correct|run|align|split|consensus|analyze|export|extract|aggregate|netseq|validate|train-polya|batch|install-aligners] ...

# Check aligner availability
rectify install-aligners --check
```

For a quick end-to-end test using bundled yeast data, see the [Quick Start](quickstart.md) guide.
