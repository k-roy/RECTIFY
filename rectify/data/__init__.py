"""
RECTIFY bundled data management.

Provides pre-processed NET-seq data bundled with the package for
automatic A-tract ambiguity resolution without requiring external downloads.

Author: Kevin R. Roy
Date: 2026-03-18
"""

import gzip
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import numpy as np

# Path to bundled data files (relative to this module)
BUNDLED_DATA_DIR = Path(__file__).parent

# Bundled organism data
BUNDLED_DATA = {
    'saccharomyces_cerevisiae': {
        'netseq_pan': {
            'file': 'saccharomyces_cerevisiae_netseq_pan.tsv.gz',
            'description': 'Pan-mutant NET-seq consensus (WT + DST1D, 6 replicates)',
            'citation': 'Churchman LS, Weissman JS. Nature 2011; Harlen KM et al. Cell 2016',
            'has_wt_signal': True,  # Contains separate WT signal column for weighting
        },
        'netseq_wt': {
            'file': 'saccharomyces_cerevisiae_netseq_wt.tsv.gz',
            'description': 'Wild-type NET-seq only (3 replicates)',
            'citation': 'Churchman LS, Weissman JS. Nature 2011',
            'has_wt_signal': False,
        }
    },
}

# Organism aliases
ORGANISM_ALIASES = {
    'yeast': 'saccharomyces_cerevisiae',
    's_cerevisiae': 'saccharomyces_cerevisiae',
    'saccer': 'saccharomyces_cerevisiae',
    'saccer3': 'saccharomyces_cerevisiae',
    'sc': 'saccharomyces_cerevisiae',
}

# Cache for loaded NET-seq data
# Maps organism -> (signal_dict, wt_signal_dict)
# signal_dict: (chrom, strand, pos) -> total signal
# wt_signal_dict: (chrom, strand, pos) -> WT signal (for preferential weighting)
_NETSEQ_CACHE: Dict[str, Tuple[Dict[Tuple[str, str, int], float], Dict[Tuple[str, str, int], float]]] = {}


def normalize_organism(organism: str) -> str:
    """Normalize organism name to canonical form."""
    org_lower = organism.lower().replace(' ', '_').replace('-', '_')
    return ORGANISM_ALIASES.get(org_lower, org_lower)


def get_bundled_netseq_path(organism: str, dataset: str = 'netseq_pan') -> Optional[Path]:
    """
    Get path to bundled NET-seq data file for an organism.

    Args:
        organism: Organism name (e.g., 'yeast', 'saccharomyces_cerevisiae')
        dataset: Which dataset to use ('netseq_pan' or 'netseq_wt')

    Returns:
        Path to bundled data file, or None if not available
    """
    org = normalize_organism(organism)
    if org not in BUNDLED_DATA:
        return None

    # Try requested dataset first, fall back to netseq_pan, then netseq_wt
    for ds in [dataset, 'netseq_pan', 'netseq_wt']:
        data_info = BUNDLED_DATA[org].get(ds)
        if data_info:
            data_path = BUNDLED_DATA_DIR / data_info['file']
            if data_path.exists():
                return data_path

    return None


def load_bundled_netseq(
    organism: str,
    return_wt_signal: bool = False
) -> Dict[Tuple[str, str, int], float]:
    """
    Load bundled NET-seq signal data for an organism.

    Returns a dictionary mapping (chrom, strand, position) -> signal value.
    This is used for A-tract ambiguity resolution via NNLS deconvolution.

    Args:
        organism: Organism name
        return_wt_signal: If True, return WT-specific signal dict (for preferential weighting)

    Returns:
        Dict mapping (chrom, strand, position) -> signal
        If return_wt_signal=True, returns WT signal dict instead of total

    Raises:
        ValueError: If no bundled data available for organism
    """
    org = normalize_organism(organism)

    # Check cache
    if org in _NETSEQ_CACHE:
        signal_dict, wt_signal_dict = _NETSEQ_CACHE[org]
        return wt_signal_dict if return_wt_signal else signal_dict

    data_path = get_bundled_netseq_path(org)
    if data_path is None:
        raise ValueError(
            f"No bundled NET-seq data for '{organism}'. "
            f"Available: {list(BUNDLED_DATA.keys())}"
        )

    # Load the TSV
    signal_dict: Dict[Tuple[str, str, int], float] = {}
    wt_signal_dict: Dict[Tuple[str, str, int], float] = {}
    open_func = gzip.open if str(data_path).endswith('.gz') else open

    with open_func(data_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        has_wt_col = len(header) >= 5 and header[4] == 'wt_signal'

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom, pos, strand, signal = parts[:4]
                key = (chrom, strand, int(pos))
                signal_dict[key] = float(signal)

                # Load WT signal if available
                if has_wt_col and len(parts) >= 5:
                    wt_sig = float(parts[4])
                    if wt_sig > 0:
                        wt_signal_dict[key] = wt_sig
                else:
                    # No WT column = all signal is WT
                    wt_signal_dict[key] = float(signal)

    # Cache both dicts
    _NETSEQ_CACHE[org] = (signal_dict, wt_signal_dict)
    return wt_signal_dict if return_wt_signal else signal_dict


def get_netseq_signal(
    organism: str,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    wt_only: bool = False,
) -> np.ndarray:
    """
    Get NET-seq signal for a genomic region from bundled data.

    Args:
        organism: Organism name
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (exclusive)
        strand: Strand ('+' or '-')
        wt_only: If True, return only WT signal (for preferential weighting)

    Returns:
        Array of signal values (one per base position)
    """
    signal_dict = load_bundled_netseq(organism, return_wt_signal=wt_only)

    result = np.zeros(end - start)
    for pos in range(start, end):
        key = (chrom, strand, pos)
        if key in signal_dict:
            result[pos - start] = signal_dict[key]

    return result


def get_netseq_signal_with_wt(
    organism: str,
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get both total and WT-specific NET-seq signal for a region.

    This is useful for apportionment where WT ratios should be preferred.

    Args:
        organism: Organism name
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (exclusive)
        strand: Strand ('+' or '-')

    Returns:
        Tuple of (total_signal, wt_signal) arrays
    """
    total_signal = get_netseq_signal(organism, chrom, start, end, strand, wt_only=False)
    wt_signal = get_netseq_signal(organism, chrom, start, end, strand, wt_only=True)
    return total_signal, wt_signal


def is_bundled_data_available(organism: str) -> bool:
    """Check if bundled NET-seq data is available for an organism."""
    org = normalize_organism(organism)
    return get_bundled_netseq_path(org) is not None


def get_available_organisms() -> List[str]:
    """Get list of organisms with bundled data available."""
    return [org for org in BUNDLED_DATA.keys()
            if get_bundled_netseq_path(org) is not None]


def ensure_netseq_data(
    organism: str,
    custom_dir: Optional[Path] = None,
    auto_download: bool = True,  # Kept for API compatibility, but ignored
    verbose: bool = True,
) -> Optional[str]:
    """
    Ensure NET-seq data is available for an organism.

    For organisms with bundled data (e.g., yeast), returns the organism name
    to signal that bundled data should be used.

    For custom data, validates the directory exists and returns 'custom'.

    Args:
        organism: Organism name
        custom_dir: Optional custom NET-seq directory (BigWig files)
        auto_download: Ignored (kept for API compatibility)
        verbose: Print progress messages

    Returns:
        'bundled' if using bundled data
        'custom' if using custom directory
        None if no data available
    """
    # Custom dir takes priority
    if custom_dir is not None:
        custom_path = Path(custom_dir)
        if not custom_path.exists():
            raise ValueError(f"Custom NET-seq directory not found: {custom_dir}")
        if verbose:
            print(f"Using custom NET-seq data: {custom_dir}")
        return 'custom'

    org = normalize_organism(organism)

    # Check for bundled data
    if is_bundled_data_available(org):
        if verbose:
            # Prefer pan-mutant data if available
            data_info = BUNDLED_DATA[org].get('netseq_pan') or BUNDLED_DATA[org].get('netseq_wt')
            print(f"Using bundled NET-seq data for {org}")
            print(f"  Dataset: {data_info['description']}")
            print(f"  Source: {data_info['citation']}")
        return 'bundled'

    # Not available
    if verbose:
        print(f"No bundled NET-seq data for '{organism}'.")
        print("Provide custom data with --netseq-dir, or run without NET-seq refinement.")

    return None


# =============================================================================
# Organism detection functions
# =============================================================================

def detect_organism_from_genome(genome_path: Path) -> Optional[str]:
    """
    Auto-detect organism from genome FASTA file.

    Detection is based on:
    - Genome size
    - Chromosome naming conventions

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
    if annotation_path and Path(annotation_path).exists():
        org = detect_organism_from_annotation(Path(annotation_path))
        if org:
            return org

    # Try genome
    if genome_path and Path(genome_path).exists():
        org = detect_organism_from_genome(Path(genome_path))
        if org:
            return org

    return None


# =============================================================================
# Motif Database Functions
# =============================================================================

# Bundled motif databases
BUNDLED_MOTIF_DATABASES = {
    'saccharomyces_cerevisiae': {
        'cpa_factors': {
            'file': 'motif_databases/scerevisiae_tf_motifs.meme',
            'description': 'Yeast CPA factors, NNS pathway, and general TFs',
        },
    },
}


def get_motif_database_path(
    organism: str,
    database: str = 'cpa_factors',
) -> Optional[Path]:
    """
    Get path to bundled motif database for an organism.

    Args:
        organism: Organism name (e.g., 'yeast', 'saccharomyces_cerevisiae')
        database: Database name (default: 'cpa_factors')

    Returns:
        Path to MEME-format motif database, or None if not available
    """
    org = normalize_organism(organism)

    if org not in BUNDLED_MOTIF_DATABASES:
        return None

    db_info = BUNDLED_MOTIF_DATABASES[org].get(database)
    if db_info is None:
        return None

    db_path = BUNDLED_DATA_DIR / db_info['file']
    if db_path.exists():
        return db_path

    return None


def get_available_motif_databases(organism: str = None) -> Dict[str, List[str]]:
    """
    Get available motif databases.

    Args:
        organism: Optional organism to filter by

    Returns:
        Dict mapping organism -> list of database names
    """
    result = {}

    for org, databases in BUNDLED_MOTIF_DATABASES.items():
        if organism and normalize_organism(organism) != org:
            continue

        available = []
        for db_name, db_info in databases.items():
            db_path = BUNDLED_DATA_DIR / db_info['file']
            if db_path.exists():
                available.append(db_name)

        if available:
            result[org] = available

    return result


# =============================================================================
# Bundled Reference Genome and Annotation Functions
# =============================================================================

# Bundled genomes and annotations
BUNDLED_GENOMES = {
    'saccharomyces_cerevisiae': {
        'genome': {
            'file': 'genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz',
            'version': 'R64-5-1',
            'source': 'SGD (Saccharomyces Genome Database)',
            'size_mb': 12,  # Uncompressed size
        },
        'annotation': {
            'file': 'genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz',
            'version': 'R64-5-1',
            'source': 'SGD',
            'format': 'GFF3',
        },
        'go_annotations': {
            'file': 'genomes/saccharomyces_cerevisiae/go_annotations.tsv.gz',
            'version': 'R64-5-1',
            'source': 'SGD',
        },
    },
}


def get_bundled_genome_path(organism: str) -> Optional[Path]:
    """
    Get path to bundled reference genome for an organism.

    Args:
        organism: Organism name (e.g., 'yeast', 'saccharomyces_cerevisiae')

    Returns:
        Path to bundled genome FASTA (gzipped), or None if not available
    """
    org = normalize_organism(organism)
    if org not in BUNDLED_GENOMES:
        return None

    genome_info = BUNDLED_GENOMES[org].get('genome')
    if genome_info is None:
        return None

    genome_path = BUNDLED_DATA_DIR / genome_info['file']
    if genome_path.exists():
        return genome_path

    return None


def get_bundled_annotation_path(organism: str) -> Optional[Path]:
    """
    Get path to bundled annotation file for an organism.

    Args:
        organism: Organism name (e.g., 'yeast', 'saccharomyces_cerevisiae')

    Returns:
        Path to bundled GFF annotation file (gzipped), or None if not available
    """
    org = normalize_organism(organism)
    if org not in BUNDLED_GENOMES:
        return None

    ann_info = BUNDLED_GENOMES[org].get('annotation')
    if ann_info is None:
        return None

    ann_path = BUNDLED_DATA_DIR / ann_info['file']
    if ann_path.exists():
        return ann_path

    return None


def get_bundled_go_annotations_path(organism: str) -> Optional[Path]:
    """
    Get path to bundled GO annotations for an organism.

    Args:
        organism: Organism name

    Returns:
        Path to GO annotations file (gzipped TSV), or None if not available
    """
    org = normalize_organism(organism)
    if org not in BUNDLED_GENOMES:
        return None

    go_info = BUNDLED_GENOMES[org].get('go_annotations')
    if go_info is None:
        return None

    go_path = BUNDLED_DATA_DIR / go_info['file']
    if go_path.exists():
        return go_path

    return None


def is_bundled_genome_available(organism: str) -> bool:
    """Check if bundled genome is available for an organism."""
    return get_bundled_genome_path(organism) is not None


def is_bundled_annotation_available(organism: str) -> bool:
    """Check if bundled annotation is available for an organism."""
    return get_bundled_annotation_path(organism) is not None


def get_bundled_genome_info(organism: str) -> Optional[Dict]:
    """Get metadata about bundled genome."""
    org = normalize_organism(organism)
    if org in BUNDLED_GENOMES:
        return BUNDLED_GENOMES[org].get('genome')
    return None


def get_available_bundled_genomes() -> List[str]:
    """Get list of organisms with bundled genomes available."""
    available = []
    for org in BUNDLED_GENOMES.keys():
        if get_bundled_genome_path(org) is not None:
            available.append(org)
    return available


def ensure_reference_data(
    organism: str,
    custom_genome: Optional[Path] = None,
    custom_annotation: Optional[Path] = None,
    verbose: bool = True,
) -> Tuple[Optional[Path], Optional[Path], str]:
    """
    Ensure reference genome and annotation are available.

    Uses bundled data for supported organisms, or validates custom paths.

    Args:
        organism: Organism name
        custom_genome: Optional path to custom genome FASTA
        custom_annotation: Optional path to custom annotation GFF/GTF
        verbose: Print progress messages

    Returns:
        Tuple of (genome_path, annotation_path, data_source)
        data_source is 'bundled', 'custom', or 'none'
    """
    org = normalize_organism(organism)

    # Custom paths take priority
    if custom_genome is not None or custom_annotation is not None:
        genome_path = Path(custom_genome) if custom_genome else None
        ann_path = Path(custom_annotation) if custom_annotation else None

        if genome_path and not genome_path.exists():
            raise ValueError(f"Custom genome not found: {custom_genome}")
        if ann_path and not ann_path.exists():
            raise ValueError(f"Custom annotation not found: {custom_annotation}")

        if verbose:
            if genome_path:
                print(f"Using custom genome: {genome_path}")
            if ann_path:
                print(f"Using custom annotation: {ann_path}")

        return genome_path, ann_path, 'custom'

    # Check for bundled data
    genome_path = get_bundled_genome_path(org)
    ann_path = get_bundled_annotation_path(org)

    if genome_path or ann_path:
        if verbose:
            genome_info = get_bundled_genome_info(org)
            print(f"Using bundled data for {org}")
            if genome_info:
                print(f"  Genome version: {genome_info.get('version', 'unknown')}")
                print(f"  Source: {genome_info.get('source', 'unknown')}")

        return genome_path, ann_path, 'bundled'

    # No data available
    if verbose:
        print(f"No bundled reference data for '{organism}'.")
        print("Provide custom paths with --genome and/or --annotation.")

    return None, None, 'none'


def load_bundled_go_annotations(organism: str) -> Optional['pd.DataFrame']:
    """
    Load bundled GO annotations for an organism.

    Args:
        organism: Organism name

    Returns:
        DataFrame with columns: gene_name, go_id, go_category, description
        Returns None if not available or pandas not installed
    """
    try:
        import pandas as pd
    except ImportError:
        return None

    go_path = get_bundled_go_annotations_path(organism)
    if go_path is None:
        return None

    return pd.read_csv(go_path, sep='\t', compression='gzip')
