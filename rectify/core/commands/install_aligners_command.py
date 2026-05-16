"""
Install-Aligners Command for RECTIFY.

Downloads and/or compiles external aligners required by rectify align.

Aligner install matrix:
  aligner     pip-installable  conda (bioconda)  vendored binary  source build
  ----------  ---------------  ----------------  ---------------  ------------
  minimap2    no (C binary)    yes               no               yes (GitHub)
  mapPacBio   no (Java/BBMap)  yes (bbmap)       no               no
  gapmm2      yes              yes               no               no
  uLTRA       yes              yes               no               no
  deSALT      no               yes               linux/x86_64     yes (GitHub)

The vendored deSALT binary (linux/x86_64, glibc 2.6+) is bundled at
  rectify/data/bin/linux_x86_64/deSALT
and is used automatically when deSALT is not on PATH.

For other platforms (macOS, linux/aarch64) this command downloads the deSALT
source from GitHub and compiles it (requires gcc and zlib-dev).

Author: Kevin R. Roy
"""

import argparse
import logging
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ─── aligner metadata ────────────────────────────────────────────────────────

DESALT_VERSION = '1.5.6'
DESALT_GITHUB = 'https://github.com/ydLiu-HIT/deSALT'
DESALT_RELEASE_URL = (
    f'https://github.com/ydLiu-HIT/deSALT/archive/refs/tags/v{DESALT_VERSION}.tar.gz'
)

MINIMAP2_VERSION = '2.28'
MINIMAP2_RELEASE_URL = (
    f'https://github.com/lh3/minimap2/releases/download/v{MINIMAP2_VERSION}/'
    f'minimap2-{MINIMAP2_VERSION}_x64-linux.tar.bz2'
)


# ─── helpers ─────────────────────────────────────────────────────────────────

def _install_dir() -> Path:
    """Return the default install directory: ~/.rectify/bin/"""
    d = Path.home() / '.rectify' / 'bin'
    d.mkdir(parents=True, exist_ok=True)
    return d


def _run(cmd, **kwargs):
    """Run a command, logging it first, raising on non-zero exit."""
    logger.info(f"  $ {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, check=True, **kwargs)


def _check_available(name: str) -> Optional[str]:
    path = shutil.which(name)
    if path:
        logger.info(f"  {name}: already on PATH → {path}")
    return path


def _add_to_path_hint(install_dir: Path):
    logger.info("")
    logger.info(f"Add to PATH to use installed aligners:")
    logger.info(f"  export PATH=\"{install_dir}:$PATH\"")
    logger.info("")
    logger.info("To make permanent, add the above line to ~/.bashrc or ~/.zshrc")


# ─── deSALT ──────────────────────────────────────────────────────────────────

def _install_desalt_from_source(install_dir: Path) -> Path:
    """Download deSALT source from GitHub, compile, install binary."""
    logger.info(f"Building deSALT {DESALT_VERSION} from source...")

    # Check build dependencies
    for dep in ('gcc', 'make'):
        if not shutil.which(dep):
            raise RuntimeError(
                f"Build tool '{dep}' not found. "
                "Install with your OS package manager (apt/yum/brew)."
            )

    with tempfile.TemporaryDirectory(prefix='rectify_desalt_build_') as tmpdir:
        src_archive = os.path.join(tmpdir, f'desalt-{DESALT_VERSION}.tar.gz')

        # Download
        logger.info(f"  Downloading: {DESALT_RELEASE_URL}")
        try:
            import urllib.request
            urllib.request.urlretrieve(DESALT_RELEASE_URL, src_archive)
        except Exception as e:
            raise RuntimeError(
                f"Download failed: {e}\n"
                f"Try manually: wget {DESALT_RELEASE_URL} -O {src_archive}"
            ) from e

        # Extract
        import tarfile
        with tarfile.open(src_archive, 'r:gz') as tf:
            tf.extractall(tmpdir)

        # Find extracted directory
        extracted = [
            d for d in os.listdir(tmpdir)
            if os.path.isdir(os.path.join(tmpdir, d)) and d.startswith('deSALT')
        ]
        if not extracted:
            raise RuntimeError("Could not find extracted deSALT source directory")
        src_dir = os.path.join(tmpdir, extracted[0], 'src')

        # Compile
        logger.info("  Compiling deSALT (make)...")
        _run(['make', '-C', src_dir, '-j', '4'])

        # Install
        built_binary = os.path.join(src_dir, 'deSALT')
        if not os.path.isfile(built_binary):
            raise RuntimeError(f"Compiled binary not found at {built_binary}")
        dest = install_dir / 'deSALT'
        shutil.copy2(built_binary, str(dest))
        dest.chmod(0o755)

    logger.info(f"  deSALT installed: {dest}")
    return dest


def install_desalt(install_dir: Path, force: bool = False) -> bool:
    """Install deSALT. Returns True if newly installed or already present."""
    logger.info("=== deSALT ===")

    # Already on PATH?
    if _check_available('deSALT') and not force:
        return True

    # Vendored binary available?
    from rectify.core.multi_aligner import _get_vendored_desalt
    vendored = _get_vendored_desalt()
    if vendored and not force:
        logger.info(f"  Vendored binary available: {vendored}")
        logger.info("  (will be used automatically — no install needed)")
        return True

    # Platform-specific action
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == 'linux' and machine in ('x86_64', 'amd64'):
        # Use vendored binary: copy to install_dir so it's on PATH
        logger.info("  Copying vendored Linux/x86_64 binary to install dir...")
        here = Path(__file__).parent.parent
        vendored_src = here / 'data' / 'bin' / 'linux_x86_64' / 'deSALT'
        if vendored_src.exists():
            dest = install_dir / 'deSALT'
            shutil.copy2(str(vendored_src), str(dest))
            dest.chmod(0o755)
            logger.info(f"  deSALT installed: {dest}")
            return True
        # Vendored file missing (shouldn't happen in a proper install)

    # Fall back to source build
    logger.info(f"  No pre-built binary for {system}/{machine} — building from source...")
    try:
        _install_desalt_from_source(install_dir)
        return True
    except Exception as e:
        logger.error(f"  Source build failed: {e}")
        logger.error(
            f"  Manual install: conda install -c bioconda desalt={DESALT_VERSION}"
        )
        return False


# ─── minimap2 ────────────────────────────────────────────────────────────────

def install_minimap2(install_dir: Path, force: bool = False) -> bool:
    """Install minimap2 pre-built binary (Linux/x86_64) or guide user."""
    logger.info("=== minimap2 ===")

    if _check_available('minimap2') and not force:
        return True

    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == 'linux' and machine in ('x86_64', 'amd64'):
        logger.info(f"  Downloading minimap2 {MINIMAP2_VERSION} pre-built binary...")
        try:
            import urllib.request, tarfile as _tf
            with tempfile.TemporaryDirectory() as tmpdir:
                archive = os.path.join(tmpdir, 'minimap2.tar.bz2')
                urllib.request.urlretrieve(MINIMAP2_RELEASE_URL, archive)
                with _tf.open(archive, 'r:bz2') as tf:
                    tf.extractall(tmpdir)
                # Find the binary
                for root, dirs, files in os.walk(tmpdir):
                    if 'minimap2' in files:
                        src = os.path.join(root, 'minimap2')
                        dest = install_dir / 'minimap2'
                        shutil.copy2(src, str(dest))
                        dest.chmod(0o755)
                        logger.info(f"  minimap2 installed: {dest}")
                        return True
        except Exception as e:
            logger.warning(f"  Pre-built download failed: {e}")

    logger.info("  Install minimap2 via conda:")
    logger.info("    conda install -c bioconda minimap2")
    return False


# ─── pip-installable tools ───────────────────────────────────────────────────

def install_gapmm2(force: bool = False) -> bool:
    logger.info("=== gapmm2 ===")
    if _check_available('gapmm2') and not force:
        return True
    logger.info("  Installing gapmm2 via pip...")
    try:
        _run([sys.executable, '-m', 'pip', 'install', 'gapmm2>=0.2'])
        return True
    except subprocess.CalledProcessError:
        logger.error("  pip install gapmm2 failed")
        return False


def install_ultra(force: bool = False) -> bool:
    logger.info("=== uLTRA ===")
    ultra_ok = _check_available('uLTRA')
    if ultra_ok and not force:
        logger.info("  uLTRA already installed.")
        # namfinder is vendored — always available via rectify's bundled binary
        logger.info("  namfinder: bundled with rectify (no separate install needed).")
        return True
    logger.info("  Installing uLTRA via pip...")
    try:
        _run([sys.executable, '-m', 'pip', 'install', 'ultra-bioinformatics'])
        logger.info("  uLTRA installed.")
        logger.info("  namfinder: bundled with rectify (no separate install needed).")
        return True
    except subprocess.CalledProcessError:
        logger.error("  pip install ultra-bioinformatics failed")
        logger.error("  Try: conda install -c bioconda ultra")
        return False


def check_mapp_pacbio() -> bool:
    logger.info("=== mapPacBio (BBMap) ===")
    if _check_available('mapPacBio.sh'):
        return True
    logger.info("  mapPacBio is part of BBMap (large Java toolkit).")
    logger.info("  Cannot be pip-installed; use conda:")
    logger.info("    conda install -c bioconda bbmap")
    return False


# ─── CLI ─────────────────────────────────────────────────────────────────────

def create_install_aligners_parser(subparsers) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        'install-aligners',
        help='Download and install external aligners required by rectify align',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
The rectify align pipeline uses up to 5 aligners:
  - minimap2    (C binary, bioconda or downloaded)
  - mapPacBio   (Java/BBMap, bioconda only)
  - gapmm2      (pip-installable)
  - uLTRA       (pip-installable)
  - deSALT      (vendored Linux binary OR built from source)

By default --all installs all pip-installable tools and deSALT, and
prints conda instructions for minimap2/mapPacBio.

Examples:
  rectify install-aligners --all
  rectify install-aligners --desalt --minimap2
  rectify install-aligners --desalt --force --install-dir ~/.local/bin
        """
    )

    target = parser.add_argument_group('Aligners to install')
    target.add_argument('--all', action='store_true',
                        help='Install all supported aligners')
    target.add_argument('--desalt', action='store_true',
                        help='Install deSALT (vendored binary or source build)')
    target.add_argument('--minimap2', action='store_true',
                        help='Install minimap2 (pre-built binary, Linux/x86_64)')
    target.add_argument('--gapmm2', action='store_true',
                        help='Install gapmm2 via pip')
    target.add_argument('--ultra', action='store_true',
                        help='Install uLTRA via pip')
    target.add_argument('--mapp-pacbio', action='store_true',
                        help='Print conda install instructions for mapPacBio/BBMap')
    target.add_argument('--check', action='store_true',
                        help='Check which aligners are already available (do not install)')

    opts = parser.add_argument_group('Options')
    opts.add_argument(
        '--install-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Directory to install binaries into (default: ~/.rectify/bin/)'
    )
    opts.add_argument(
        '--force',
        action='store_true',
        help='Reinstall even if aligner is already on PATH'
    )
    opts.add_argument('--verbose', action='store_true', help='Verbose logging')
    return parser


def _check_all_aligners():
    """Print availability status of all aligners."""
    aligners = [
        ('minimap2',    'minimap2'),
        ('mapPacBio',   'mapPacBio.sh'),
        ('gapmm2',      'gapmm2'),
        ('uLTRA',       'uLTRA'),
        ('deSALT',      'deSALT'),
    ]
    logger.info("Aligner availability:")
    for name, cmd in aligners:
        path = shutil.which(cmd)
        if path:
            logger.info(f"  {name:12s}  OK  → {path}")
        else:
            # Check for vendored deSALT specifically
            if name == 'deSALT':
                from rectify.core.multi_aligner import _get_vendored_desalt
                vendored = _get_vendored_desalt()
                if vendored:
                    logger.info(f"  {name:12s}  VENDORED → {vendored}")
                    continue
            logger.info(f"  {name:12s}  MISSING")


def run_install_aligners(args: argparse.Namespace) -> int:
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s %(levelname)s %(message)s'
    )

    if args.check:
        _check_all_aligners()
        return 0

    install_dir = args.install_dir or _install_dir()
    install_dir.mkdir(parents=True, exist_ok=True)
    force = args.force

    do_all = args.all
    results = {}

    if do_all or args.desalt:
        results['deSALT'] = install_desalt(install_dir, force=force)

    if do_all or args.minimap2:
        results['minimap2'] = install_minimap2(install_dir, force=force)

    if do_all or args.gapmm2:
        results['gapmm2'] = install_gapmm2(force=force)

    if do_all or args.ultra:
        results['uLTRA'] = install_ultra(force=force)

    if do_all or args.mapp_pacbio:
        results['mapPacBio'] = check_mapp_pacbio()

    if not results:
        logger.info("No aligners specified. Use --all or --desalt/--minimap2/etc.")
        logger.info("Use --check to see current aligner availability.")
        return 1

    _add_to_path_hint(install_dir)

    logger.info("Summary:")
    ok = True
    for name, success in results.items():
        status = "OK" if success else "FAILED"
        logger.info(f"  {name:12s}  {status}")
        if not success:
            ok = False

    return 0 if ok else 1


def run(args: argparse.Namespace):
    sys.exit(run_install_aligners(args))
