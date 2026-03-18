"""
RECTIFY bundled data management.

Handles downloading and caching of organism-specific data (NET-seq, etc.)
for automatic refinement without requiring user-provided files.

Author: Kevin R. Roy
Date: 2026-03-18
"""

import os
import shutil
from pathlib import Path
from typing import Optional, Dict, List
import urllib.request
import hashlib

# Default cache directory
DEFAULT_CACHE_DIR = Path.home() / '.rectify' / 'data'

# Available organism data
# Format: organism -> {data_type -> {url, md5, files}}
AVAILABLE_DATA = {
    'saccharomyces_cerevisiae': {
        'netseq_wt': {
            'description': 'Wild-type NET-seq data (Churchman lab)',
            'citation': 'Churchman & Weissman, Nature 2011; Harlen et al, Mol Cell 2016',
            'files': [
                'netseq_wt.plus.bw',
                'netseq_wt.minus.bw',
            ],
            # URL will point to GitHub releases or Zenodo
            'url_base': 'https://github.com/k-roy/RECTIFY/releases/download/data-v1/',
            'total_size_mb': 20,
        }
    },
    # Future organisms can be added here
    # 'homo_sapiens': {...},
    # 'mus_musculus': {...},
}

# Organism aliases
ORGANISM_ALIASES = {
    'yeast': 'saccharomyces_cerevisiae',
    's_cerevisiae': 'saccharomyces_cerevisiae',
    'saccer': 'saccharomyces_cerevisiae',
    'saccer3': 'saccharomyces_cerevisiae',
    'sc': 'saccharomyces_cerevisiae',
}


def detect_organism_from_genome(genome_path: Path) -> Optional[str]:
    """
    Auto-detect organism from genome FASTA file.

    Detection is based on:
    - Genome size
    - Chromosome naming conventions
    - Known sequence signatures

    Args:
        genome_path: Path to genome FASTA (or .fai index)

    Returns:
        Detected organism name, or None if unknown
    """
    import re

    # Try to read .fai index first (faster)
    fai_path = Path(str(genome_path) + '.fai')
    if fai_path.exists():
        try:
            chroms = {}
            total_size = 0
            with open(fai_path) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        chroms[parts[0]] = int(parts[1])
                        total_size += int(parts[1])

            # S. cerevisiae detection:
            # - ~12 Mb genome
            # - Chromosomes named chrI-chrXVI or I-XVI
            yeast_chroms = {'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI',
                           'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII',
                           'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrM', 'chrMito',
                           'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX',
                           'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'Mito'}

            yeast_matches = sum(1 for c in chroms if c in yeast_chroms)

            # Yeast: 12-13 Mb, mostly yeast chromosome names
            if 10_000_000 < total_size < 15_000_000 and yeast_matches >= 10:
                return 'saccharomyces_cerevisiae'

            # Human: ~3 Gb, chr1-chr22 + chrX/Y
            if 2_500_000_000 < total_size < 3_500_000_000:
                human_chroms = {'chr1', 'chr2', 'chr3', 'chrX', 'chrY'}
                if any(c in chroms for c in human_chroms):
                    return 'homo_sapiens'

            # Mouse: ~2.7 Gb
            if 2_000_000_000 < total_size < 3_000_000_000:
                mouse_chroms = {'chr1', 'chr2', 'chrX', 'chrY'}
                if any(c in chroms for c in mouse_chroms):
                    return 'mus_musculus'

        except Exception:
            pass

    # Fallback: read first few lines of FASTA
    try:
        with open(genome_path) as f:
            header_lines = []
            for i, line in enumerate(f):
                if line.startswith('>'):
                    header_lines.append(line.strip())
                if i > 100:  # Only check first 100 lines
                    break

            # Check for yeast chromosome patterns in headers
            yeast_pattern = re.compile(r'>chr(I{1,3}|IV|V|VI{1,3}|IX|X{1,3}|XI{1,3}|XIV|XV|XVI|M)', re.IGNORECASE)
            yeast_matches = sum(1 for h in header_lines if yeast_pattern.match(h))
            if yeast_matches >= 5:
                return 'saccharomyces_cerevisiae'

    except Exception:
        pass

    return None


def detect_organism_from_annotation(annotation_path: Path) -> Optional[str]:
    """
    Auto-detect organism from annotation GTF/GFF file.

    Args:
        annotation_path: Path to annotation file

    Returns:
        Detected organism name, or None if unknown
    """
    try:
        with open(annotation_path) as f:
            for i, line in enumerate(f):
                if i > 50:  # Only check first 50 lines
                    break

                line_lower = line.lower()

                # Check for organism indicators in comments/metadata
                if 'saccharomyces' in line_lower or 'cerevisiae' in line_lower:
                    return 'saccharomyces_cerevisiae'
                if 'homo sapiens' in line_lower or 'human' in line_lower:
                    return 'homo_sapiens'
                if 'mus musculus' in line_lower or 'mouse' in line_lower:
                    return 'mus_musculus'

                # Check for SGD (yeast) gene naming
                if line.startswith('#'):
                    continue
                if 'SGD' in line or 'gene_id "Y' in line:
                    return 'saccharomyces_cerevisiae'

    except Exception:
        pass

    return None


def detect_organism(
    genome_path: Optional[Path] = None,
    annotation_path: Optional[Path] = None,
) -> Optional[str]:
    """
    Auto-detect organism from genome and/or annotation files.

    Args:
        genome_path: Optional path to genome FASTA
        annotation_path: Optional path to annotation GTF/GFF

    Returns:
        Detected organism name, or None if unknown
    """
    # Try annotation first (often has explicit organism info)
    if annotation_path and annotation_path.exists():
        org = detect_organism_from_annotation(annotation_path)
        if org:
            return org

    # Try genome
    if genome_path and genome_path.exists():
        org = detect_organism_from_genome(genome_path)
        if org:
            return org

    return None


def get_cache_dir() -> Path:
    """Get the cache directory, creating if needed."""
    cache_dir = Path(os.environ.get('RECTIFY_CACHE_DIR', DEFAULT_CACHE_DIR))
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def normalize_organism(organism: str) -> str:
    """Normalize organism name to canonical form."""
    org_lower = organism.lower().replace(' ', '_').replace('-', '_')
    return ORGANISM_ALIASES.get(org_lower, org_lower)


def get_available_organisms() -> List[str]:
    """Get list of organisms with available data."""
    return list(AVAILABLE_DATA.keys())


def is_data_available(organism: str, data_type: str = 'netseq_wt') -> bool:
    """Check if data is available for an organism."""
    org = normalize_organism(organism)
    return org in AVAILABLE_DATA and data_type in AVAILABLE_DATA[org]


def is_data_cached(organism: str, data_type: str = 'netseq_wt') -> bool:
    """Check if data is already downloaded and cached."""
    org = normalize_organism(organism)
    if not is_data_available(org, data_type):
        return False

    cache_dir = get_cache_dir() / org / data_type
    if not cache_dir.exists():
        return False

    # Check all expected files exist
    data_info = AVAILABLE_DATA[org][data_type]
    for filename in data_info['files']:
        if not (cache_dir / filename).exists():
            return False

    return True


def get_netseq_dir(organism: str, custom_dir: Optional[Path] = None) -> Optional[Path]:
    """
    Get NET-seq directory for an organism.

    Priority:
    1. Custom directory if provided
    2. Cached bundled data if available
    3. None if no data available

    Args:
        organism: Organism name (e.g., 'yeast', 'saccharomyces_cerevisiae')
        custom_dir: Optional custom NET-seq directory (overrides bundled data)

    Returns:
        Path to NET-seq directory, or None if not available
    """
    # Custom dir takes priority
    if custom_dir is not None:
        return custom_dir

    org = normalize_organism(organism)

    # Check if bundled data is cached
    if is_data_cached(org, 'netseq_wt'):
        return get_cache_dir() / org / 'netseq_wt'

    # Check if data is available for download
    if is_data_available(org, 'netseq_wt'):
        # Data available but not downloaded - will need to download first
        return None

    return None


def download_data(
    organism: str,
    data_type: str = 'netseq_wt',
    force: bool = False,
    verbose: bool = True,
) -> Path:
    """
    Download organism data to cache.

    Args:
        organism: Organism name
        data_type: Type of data to download (default: 'netseq_wt')
        force: Force re-download even if cached
        verbose: Print progress

    Returns:
        Path to downloaded data directory

    Raises:
        ValueError: If organism/data_type not available
        RuntimeError: If download fails
    """
    org = normalize_organism(organism)

    if not is_data_available(org, data_type):
        available = get_available_organisms()
        raise ValueError(
            f"No {data_type} data available for '{organism}'. "
            f"Available organisms: {', '.join(available)}"
        )

    cache_dir = get_cache_dir() / org / data_type

    # Check if already cached
    if is_data_cached(org, data_type) and not force:
        if verbose:
            print(f"Using cached {data_type} data for {org}: {cache_dir}")
        return cache_dir

    # Download
    data_info = AVAILABLE_DATA[org][data_type]
    cache_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"Downloading {data_type} data for {org}...")
        print(f"  Source: {data_info.get('citation', 'See documentation')}")
        print(f"  Size: ~{data_info['total_size_mb']} MB")

    for filename in data_info['files']:
        url = data_info['url_base'] + filename
        dest = cache_dir / filename

        if verbose:
            print(f"  Downloading {filename}...", end=' ', flush=True)

        try:
            urllib.request.urlretrieve(url, dest)
            if verbose:
                print("OK")
        except Exception as e:
            if verbose:
                print(f"FAILED: {e}")
            raise RuntimeError(f"Failed to download {filename}: {e}")

    if verbose:
        print(f"Data cached to: {cache_dir}")

    return cache_dir


def ensure_netseq_data(
    organism: str,
    custom_dir: Optional[Path] = None,
    auto_download: bool = True,
    verbose: bool = True,
) -> Optional[Path]:
    """
    Ensure NET-seq data is available, downloading if needed.

    Args:
        organism: Organism name
        custom_dir: Optional custom NET-seq directory
        auto_download: Whether to auto-download if not cached
        verbose: Print progress

    Returns:
        Path to NET-seq directory, or None if not available
    """
    # Custom dir takes priority
    if custom_dir is not None:
        if not custom_dir.exists():
            raise ValueError(f"Custom NET-seq directory not found: {custom_dir}")
        return custom_dir

    org = normalize_organism(organism)

    # Check cache
    if is_data_cached(org, 'netseq_wt'):
        cache_dir = get_cache_dir() / org / 'netseq_wt'
        if verbose:
            print(f"Using cached NET-seq data: {cache_dir}")
        return cache_dir

    # Check if downloadable
    if is_data_available(org, 'netseq_wt'):
        if auto_download:
            return download_data(org, 'netseq_wt', verbose=verbose)
        else:
            if verbose:
                print(f"NET-seq data available for {org} but not downloaded.")
                print(f"Run: rectify download-data --organism {org}")
            return None

    # Not available
    if verbose:
        print(f"No bundled NET-seq data for '{organism}'.")
        print("Provide custom data with --netseq-dir, or run without NET-seq refinement.")

    return None


def clear_cache(organism: Optional[str] = None, data_type: Optional[str] = None):
    """
    Clear cached data.

    Args:
        organism: Optional organism to clear (None = all)
        data_type: Optional data type to clear (None = all)
    """
    cache_dir = get_cache_dir()

    if organism is None:
        # Clear all
        if cache_dir.exists():
            shutil.rmtree(cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
    else:
        org = normalize_organism(organism)
        org_dir = cache_dir / org

        if data_type is None:
            # Clear all data for organism
            if org_dir.exists():
                shutil.rmtree(org_dir)
        else:
            # Clear specific data type
            type_dir = org_dir / data_type
            if type_dir.exists():
                shutil.rmtree(type_dir)


def get_data_info() -> Dict:
    """Get information about available and cached data."""
    info = {
        'cache_dir': str(get_cache_dir()),
        'organisms': {}
    }

    for org, data_types in AVAILABLE_DATA.items():
        info['organisms'][org] = {
            'data_types': {}
        }
        for data_type, data_info in data_types.items():
            info['organisms'][org]['data_types'][data_type] = {
                'description': data_info.get('description', ''),
                'size_mb': data_info.get('total_size_mb', 0),
                'cached': is_data_cached(org, data_type),
            }

    return info
