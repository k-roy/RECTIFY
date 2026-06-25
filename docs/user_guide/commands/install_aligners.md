# rectify install-aligners

Download, install, or check external aligners used by `rectify align`. The
default long-read panel is minimap2 + mapPacBio + gapmm2; uLTRA, deSALT, and
GMAP are opt-in (`--junction-aligners`).

---

## Usage

```bash
rectify install-aligners [--all] [--desalt] [--minimap2] [--gapmm2] [--ultra] \
    [--gmap] [--mapp-pacbio] [--check] [--install-dir DIR] [--force]
```

## Examples

```bash
# Check what's already available
rectify install-aligners --check

# Install everything possible
rectify install-aligners --all

# Install deSALT only (vendored binary or source build)
rectify install-aligners --desalt

# Reinstall to a custom directory
rectify install-aligners --desalt --minimap2 --install-dir ~/.local/bin --force
```

---

## Options

| Flag | Description |
|------|-------------|
| `--check` | Print availability of all aligners (PATH + vendored); do not install |
| `--all` | Install all supported aligners |
| `--desalt` | Install deSALT (see below for platform logic) |
| `--minimap2` | Download pre-built Linux/x86_64 minimap2 binary |
| `--gapmm2` | Install gapmm2 via `pip install gapmm2==25.4.5` (pinned — see note below) |
| `--ultra` | Install uLTRA via `pip install ultra-bioinformatics` (its `namfinder` helper is vendored with RECTIFY) |
| `--gmap` | Check for GMAP and print the conda install + one-time database-build reminder (not pip-installable) |
| `--mapp-pacbio` | Print conda install instructions for BBMap/mapPacBio (not pip-installable) |
| `--install-dir DIR` | Target directory (default: `~/.rectify/bin/`) |
| `--force` | Reinstall even if the aligner is already on `PATH` |

---

## Aligner install matrix

| Aligner | Vendored binary | pip | conda | Source build |
|---------|-----------------|-----|-------|--------------|
| **minimap2** | No | No | bioconda | Yes (GitHub release binary) |
| **mapPacBio** | No | No | bioconda (`bbmap`) | No |
| **gapmm2** | No | Yes (`==25.4.5`) | bioconda | — |
| **uLTRA** | No | Yes | bioconda | — |
| **GMAP** | No | No | bioconda (`gmap`) | No |
| **deSALT** | **Yes (Linux/x86_64)** | No | bioconda | Yes (GitHub source) |

!!! warning "gapmm2 is pinned to `==25.4.5`"
    `--gapmm2` installs exactly `gapmm2==25.4.5` to match `pyproject.toml`.
    Versions `25.4.13`+ silently drop single-exon (unspliced) reads via a
    `ts:` transcript-strand tag regression — do not upgrade past the pin.

---

## deSALT install logic

1. **Already on PATH?** — use it, skip install (unless `--force`).
2. **Vendored binary available?** (`rectify/data/bin/linux_x86_64/deSALT`) — used automatically by `rectify align` even without install; `--desalt` copies it to `--install-dir` so it is on `PATH`.
3. **Other platform** (macOS, Linux/aarch64) — downloads source from [GitHub](https://github.com/ydLiu-HIT/deSALT) and compiles with `make` (requires `gcc` and `zlib`).

!!! note "Vendored binary is used automatically"
    When `deSALT` is not on `PATH`, RECTIFY falls back to the bundled Linux/x86_64
    binary silently. You do **not** need to run `install-aligners` on standard Linux
    HPC systems — the vendored binary will be used for alignment.

    Running `rectify install-aligners --desalt` copies the binary to `~/.rectify/bin/`
    so that `deSALT` is also available as a standalone command (useful for manual
    index building or inspection).

---

## Adding the install directory to PATH

After install, add the directory to your shell's `PATH` once:

```bash
echo 'export PATH="$HOME/.rectify/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

Then verify:

```bash
rectify install-aligners --check
```

---

## mapPacBio / BBMap

mapPacBio is part of the [BBMap](https://sourceforge.net/projects/bbmap/) Java toolkit.
It cannot be installed via pip. Install via conda:

```bash
conda install -c bioconda bbmap
```

The `--mapp-pacbio` flag prints this instruction and exits.

---

## GMAP

GMAP is a C toolkit (`gmap` / `gsnap` / `gmap_build`) and is not pip-installable.
Install via conda and build a database once before use:

```bash
conda install -c bioconda gmap
gmap_build -D <dir> -d <name> genome.fa
```

The `--gmap` flag checks whether `gmap` is on `PATH` and prints these
reminders; it does not install GMAP itself.

---

## Notes

- `--check` works without any network access — it only inspects `PATH` and the vendored binary directory.
- The vendored deSALT binary links against `libz` and `libgomp` from system paths (conda's `LD_LIBRARY_PATH` is stripped at runtime to avoid conflicts).
- Source builds require `gcc`, `make`, and `zlib-devel` (or equivalent OS package).
