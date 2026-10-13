#!/usr/bin/env python3
"""
De Novo Motif Discovery Module

Discovers enriched sequence motifs around CPA clusters using STREME/MEME.

Author: Kevin R. Roy
Date: 2026-03-17
"""

from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
import subprocess
import tempfile
import os
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Default window sizes (configurable)
DEFAULT_UPSTREAM_WINDOW = 100  # bp upstream of CPA site
DEFAULT_DOWNSTREAM_WINDOW = 50  # bp downstream of CPA site

# STREME/MEME parameters
DEFAULT_MIN_MOTIF_WIDTH = 6
DEFAULT_MAX_MOTIF_WIDTH = 15
DEFAULT_N_MOTIFS = 10

# Complement mapping for reverse complement
COMPLEMENT = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')


def extract_sequences_around_clusters(
    clusters_df: pd.DataFrame,
    genome_fasta: str,
    upstream_window: int = DEFAULT_UPSTREAM_WINDOW,
    downstream_window: int = DEFAULT_DOWNSTREAM_WINDOW,
    region: str = 'both',
) -> Dict[str, List[Dict]]:
    """
    Extract sequences around CPA cluster positions.

    Args:
        clusters_df: DataFrame with cluster definitions (must have
            chrom, strand, modal_position columns)
        genome_fasta: Path to genome FASTA file
        upstream_window: bp upstream of CPA site (default: 100)
        downstream_window: bp downstream of CPA site (default: 50)
        region: Which region to extract ('upstream', 'downstream', 'both')

    Returns:
        Dict with keys 'upstream' and/or 'downstream', each containing
        list of dicts with 'cluster_id', 'sequence', 'chrom', 'start', 'end', 'strand'
    """
    # Try to import pysam for FASTA access
    try:
        import pysam
    except ImportError:
        raise ImportError(
            "pysam is required for sequence extraction. "
            "Install with: pip install pysam"
        )

    # Use context manager for proper resource cleanup
    with pysam.FastaFile(genome_fasta) as fasta:
        chrom_lengths = dict(zip(fasta.references, fasta.lengths))

        results = {}

        if region in ('upstream', 'both'):
            results['upstream'] = []

        if region in ('downstream', 'both'):
            results['downstream'] = []

        skipped_chroms = set()
        extraction_errors = 0

        for _, cluster in clusters_df.iterrows():
            chrom = cluster['chrom']
            strand = cluster['strand']
            pos = cluster['modal_position']
            cluster_id = cluster.get('cluster_id', f"{chrom}_{pos}_{strand}")

            # Skip if chromosome not in FASTA
            if chrom not in chrom_lengths:
                skipped_chroms.add(chrom)
                continue

            chrom_len = chrom_lengths[chrom]

            # Calculate coordinates (strand-aware)
            if strand == '+':
                upstream_start = max(0, pos - upstream_window)
                upstream_end = pos
                downstream_start = pos
                downstream_end = min(chrom_len, pos + downstream_window)
            else:
                # For minus strand, upstream is higher coordinates
                upstream_start = pos
                upstream_end = min(chrom_len, pos + upstream_window)
                downstream_start = max(0, pos - downstream_window)
                downstream_end = pos

            # Extract upstream
            if 'upstream' in results:
                try:
                    seq = fasta.fetch(chrom, upstream_start, upstream_end)
                    if strand == '-':
                        seq = _reverse_complement(seq)
                    if len(seq) >= upstream_window * 0.8:  # At least 80% of requested
                        results['upstream'].append({
                            'cluster_id': cluster_id,
                            'sequence': seq.upper(),
                            'chrom': chrom,
                            'start': upstream_start,
                            'end': upstream_end,
                            'strand': strand,
                        })
                except Exception as e:
                    logger.warning(f"Skipping locus due to error: {e}")
                    extraction_errors += 1
                    continue

            # Extract downstream
            if 'downstream' in results:
                try:
                    seq = fasta.fetch(chrom, downstream_start, downstream_end)
                    if strand == '-':
                        seq = _reverse_complement(seq)
                    if len(seq) >= downstream_window * 0.8:
                        results['downstream'].append({
                            'cluster_id': cluster_id,
                            'sequence': seq.upper(),
                            'chrom': chrom,
                            'start': downstream_start,
                            'end': downstream_end,
                            'strand': strand,
                        })
                except Exception as e:
                    logger.warning(f"Skipping locus due to error: {e}")
                    extraction_errors += 1
                    continue

        # Log statistics
        if skipped_chroms:
            print(f"  Skipped {len(skipped_chroms)} chromosomes not in FASTA: {list(skipped_chroms)[:5]}...")
        if extraction_errors > 0:
            print(f"  {extraction_errors} extraction errors (e.g., invalid coordinates)")

    return results


def _reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def write_fasta(
    sequences: List[Dict],
    output_path: str,
    id_format: str = '{cluster_id}',
) -> int:
    """
    Write sequences to FASTA file.

    Args:
        sequences: List of dicts with 'cluster_id' and 'sequence'
        output_path: Path to output FASTA file
        id_format: Format string for sequence IDs

    Returns:
        Number of sequences written
    """
    with open(output_path, 'w') as f:
        for seq_dict in sequences:
            seq_id = id_format.format(**seq_dict)
            seq = seq_dict['sequence']
            f.write(f">{seq_id}\n{seq}\n")

    return len(sequences)


def run_motif_discovery(
    foreground_sequences: List[Dict],
    background_sequences: Optional[List[Dict]] = None,
    output_dir: str = '.',
    tool: str = 'streme',
    min_width: int = DEFAULT_MIN_MOTIF_WIDTH,
    max_width: int = DEFAULT_MAX_MOTIF_WIDTH,
    n_motifs: int = DEFAULT_N_MOTIFS,
    region_name: str = 'region',
) -> Dict:
    """
    Run de novo motif discovery using STREME or MEME.

    Args:
        foreground_sequences: List of sequence dicts (foreground/test set)
        background_sequences: Optional background sequences (for STREME)
        output_dir: Directory for output files
        tool: 'streme' or 'meme'
        min_width: Minimum motif width
        max_width: Maximum motif width
        n_motifs: Number of motifs to find
        region_name: Name for this analysis (used in output files)

    Returns:
        Dict with:
            - success: bool
            - output_dir: Path to results
            - n_foreground: Number of foreground sequences
            - n_background: Number of background sequences (if applicable)
            - motifs_found: Number of motifs found (if successful)
            - error: Error message (if failed)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write foreground FASTA
    fg_fasta = output_dir / f'{region_name}_foreground.fa'
    n_fg = write_fasta(foreground_sequences, str(fg_fasta))

    if n_fg < 10:
        return {
            'success': False,
            'error': f'Insufficient foreground sequences ({n_fg} < 10)',
            'n_foreground': n_fg,
        }

    # Write background FASTA if provided
    bg_fasta = None
    n_bg = 0
    if background_sequences:
        bg_fasta = output_dir / f'{region_name}_background.fa'
        n_bg = write_fasta(background_sequences, str(bg_fasta))

    # Check if tool is available
    tool_available = _check_tool_available(tool)
    if not tool_available:
        return {
            'success': False,
            'error': f'{tool} not found in PATH. Install MEME Suite.',
            'n_foreground': n_fg,
            'n_background': n_bg,
        }

    # Run motif discovery
    results_dir = output_dir / f'{region_name}_{tool}_results'

    try:
        if tool == 'streme':
            success = _run_streme(
                str(fg_fasta),
                str(results_dir),
                background_fasta=str(bg_fasta) if bg_fasta else None,
                min_width=min_width,
                max_width=max_width,
            )
        else:  # meme
            success = _run_meme(
                str(fg_fasta),
                str(results_dir),
                min_width=min_width,
                max_width=max_width,
                n_motifs=n_motifs,
            )

        # Count motifs found
        n_motifs_found = _count_motifs_in_results(str(results_dir), tool)

        return {
            'success': success,
            'output_dir': str(results_dir),
            'n_foreground': n_fg,
            'n_background': n_bg,
            'motifs_found': n_motifs_found,
        }

    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'n_foreground': n_fg,
            'n_background': n_bg,
        }


def _check_tool_available(tool: str) -> bool:
    """Check if MEME Suite tool is available."""
    try:
        result = subprocess.run(
            [tool, '--version'],
            capture_output=True,
            timeout=10,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def _run_streme(
    foreground_fasta: str,
    output_dir: str,
    background_fasta: Optional[str] = None,
    min_width: int = 4,
    max_width: int = 12,
    n_threads: int = 4,
) -> bool:
    """Run STREME for short motif discovery."""
    cmd = [
        'streme',
        '--p', foreground_fasta,
        '--oc', output_dir,
        '--dna',
        '--minw', str(min_width),
        '--maxw', str(max_width),
    ]

    if background_fasta:
        cmd.extend(['--n', background_fasta])

    result = subprocess.run(
        cmd,
        capture_output=True,
        timeout=600,  # 10 minute timeout
    )

    return result.returncode == 0


def _run_meme(
    fasta_file: str,
    output_dir: str,
    min_width: int = 6,
    max_width: int = 20,
    n_motifs: int = 15,
    n_threads: int = 4,
) -> bool:
    """Run MEME for de novo motif discovery with parallelization."""
    cmd = [
        'meme',
        fasta_file,
        '-dna',
        '-oc', output_dir,
        '-nmotifs', str(n_motifs),
        '-minw', str(min_width),
        '-maxw', str(max_width),
        '-mod', 'zoops',  # Zero or one per sequence
        '-revcomp',       # Search both strands
        '-p', str(n_threads),  # Parallel threads
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        timeout=1800,  # 30 minute timeout
    )

    return result.returncode == 0


def _count_motifs_in_results(results_dir: str, tool: str) -> int:
    """Count number of motifs found in results."""
    results_dir = Path(results_dir)

    if tool == 'streme':
        # STREME outputs streme.txt
        txt_file = results_dir / 'streme.txt'
        if txt_file.exists():
            with open(txt_file, 'r') as f:
                content = f.read()
                # Count "MOTIF" lines
                return content.count('\nMOTIF ')
    else:
        # MEME outputs meme.txt
        txt_file = results_dir / 'meme.txt'
        if txt_file.exists():
            with open(txt_file, 'r') as f:
                content = f.read()
                return content.count('\nMOTIF ')

    return 0


def _run_single_motif_analysis(args):
    """Worker function for parallel motif discovery."""
    name, foreground, background, output_dir = args
    result = run_motif_discovery(
        foreground_sequences=foreground,
        background_sequences=background if background else None,
        output_dir=output_dir,
        region_name=name,
    )
    return name, result


def run_differential_motif_analysis(
    enriched_clusters: pd.DataFrame,
    depleted_clusters: pd.DataFrame,
    genome_fasta: str,
    output_dir: str,
    upstream_window: int = DEFAULT_UPSTREAM_WINDOW,
    downstream_window: int = DEFAULT_DOWNSTREAM_WINDOW,
    n_parallel: int = 4,
    background_clusters: Optional[pd.DataFrame] = None,
    min_foreground: int = 0,
    max_background_ratio: int = 5,
) -> Dict[str, Dict]:
    """
    Run motif discovery for enriched and depleted clusters, upstream and downstream (up to 4 analyses).

    Each side is judged on its own: a side with fewer than *min_foreground* clusters is skipped (its two
    entries carry the reason), and the other side still runs. Before 2026-10-13 the caller required both
    sides to reach 20, so one thin side silenced all four scans (the deliverable spec in TODO_KR_NOTES asks
    for enriched vs background and depleted vs background).

    Background for a side: *background_clusters* (e.g. the unchanged clusters), subsampled with a fixed
    seed to at most *max_background_ratio* times the side's size, when it holds at least *min_foreground*
    clusters; otherwise the other side when that side qualifies (the earlier enriched-vs-depleted
    contrast); otherwise none (STREME then uses its shuffled control).

    Args:
        enriched_clusters: Clusters with increased signal
        depleted_clusters: Clusters with decreased signal
        genome_fasta: Path to genome FASTA
        output_dir: Output directory
        upstream_window: bp upstream (default: 100)
        downstream_window: bp downstream (default: 50)
        n_parallel: Number of parallel workers (default: 4); 1 runs in this process
        background_clusters: Optional background set (same columns as the foreground tables)
        min_foreground: Smallest side that is scanned (0 = scan whatever is given)
        max_background_ratio: Background cap as a multiple of the side's size

    Returns:
        Dict mapping analysis_name -> motif discovery results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    def _seqs(df):
        if df is None or len(df) == 0:
            return {}
        return extract_sequences_around_clusters(
            df, genome_fasta, upstream_window=upstream_window, downstream_window=downstream_window,
        )

    print("Extracting sequences around clusters...")
    sides = {'enriched': enriched_clusters, 'depleted': depleted_clusters}
    ok = {name: df is not None and len(df) >= max(min_foreground, 1) for name, df in sides.items()}
    seqs = {name: _seqs(df) if ok[name] else {} for name, df in sides.items()}
    bg_ok = background_clusters is not None and len(background_clusters) >= max(min_foreground, 1)

    analyses = []
    for name, other in (('enriched', 'depleted'), ('depleted', 'enriched')):
        if not ok[name]:
            n = 0 if sides[name] is None else len(sides[name])
            for region in ('upstream', 'downstream'):
                results[f'{name}_{region}'] = {'success': False, 'n_foreground': n, 'n_background': 0,
                                               'error': f'skipped: {n} {name} clusters < {min_foreground}'}
            continue
        if bg_ok:
            cap = max_background_ratio * len(sides[name])
            bg_df = background_clusters
            if len(bg_df) > cap:
                bg_df = bg_df.sample(n=cap, random_state=0)
            bg = _seqs(bg_df)
        elif ok[other]:
            bg = seqs[other]
        else:
            bg = {}
        for region in ('upstream', 'downstream'):
            analyses.append((f'{name}_{region}', seqs[name].get(region, []), bg.get(region, []), str(output_dir)))

    if n_parallel <= 1:
        for args in analyses:
            result_name, result = _run_single_motif_analysis(args)
            results[result_name] = result
            print(f"  {result_name}: {result.get('motifs_found', 0)} motifs found")
        return results

    # Run analyses in parallel
    print(f"Running {len(analyses)} motif analyses in parallel ({n_parallel} workers)...")
    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = {executor.submit(_run_single_motif_analysis, args): args[0] for args in analyses}
        for future in as_completed(futures):
            name = futures[future]
            try:
                result_name, result = future.result()
                results[result_name] = result
                print(f"  {result_name}: {result.get('motifs_found', 0)} motifs found")
            except Exception as e:
                results[name] = {'success': False, 'error': str(e)}
                print(f"  {name}: FAILED - {e}")

    return results


def summarize_motif_results(results: Dict[str, Dict]) -> pd.DataFrame:
    """
    Create summary table of motif discovery results.

    Args:
        results: Results from run_differential_motif_analysis()

    Returns:
        Summary DataFrame
    """
    rows = []
    for name, result in results.items():
        rows.append({
            'analysis': name,
            'success': result.get('success', False),
            'n_foreground': result.get('n_foreground', 0),
            'n_background': result.get('n_background', 0),
            'motifs_found': result.get('motifs_found', 0),
            'output_dir': result.get('output_dir', ''),
            'error': result.get('error', ''),
        })

    return pd.DataFrame(rows)


# =============================================================================
# Motif Database Matching
# =============================================================================

# Default motif database paths (MEME format)
YEAST_MOTIF_DB = {
    'jaspar_fungi': 'JASPAR2022_CORE_fungi_non-redundant_pfms_meme.txt',
    'yeastract': 'yeastract_motifs.meme',
    'scertf': 'ScerTF_motifs.meme',
}

# Known yeast transcription factor binding sites
YEAST_TF_MOTIFS = {
    # CPA-related factors
    'Hrp1': 'UAUAUA',  # Hrp1/Nab4 - poly(A) signal recognition
    'Rna15': 'AAUAAA',  # Rna15 - canonical poly(A) signal
    'Fip1': 'UUUUUU',  # Fip1 - U-rich element
    'Yth1': 'AAGAA',   # Yth1 - near CPA site
    'Pab1': 'AAAAAAA', # Pab1 - poly(A) binding

    # Termination-related
    'Nrd1': 'GUAA|UGUA',  # Nrd1-Nab3-Sen1 pathway
    'Nab3': 'UCUU',
    'Sen1': 'G-rich',

    # General TFs with A/T-rich motifs
    'Reb1': 'TTACCCG',
    'Abf1': 'RTCRYNNNNNACG',
    'Rap1': 'ACACCCATACATC',
    'Gcr1': 'CTTCC',
}


def get_bundled_motif_database(organism: str = 'scerevisiae') -> Optional[str]:
    """
    Get path to bundled motif database for an organism.

    Args:
        organism: Organism name ('scerevisiae', 'human', etc.)

    Returns:
        Path to MEME-format motif database, or None if not available
    """
    # Check if motif database is bundled with rectify
    pkg_dir = Path(__file__).parent.parent.parent / 'data' / 'motif_databases'

    if organism == 'scerevisiae':
        db_path = pkg_dir / 'scerevisiae_tf_motifs.meme'
        if db_path.exists():
            return str(db_path)

        # Try to find JASPAR fungi database
        for db_name, db_file in YEAST_MOTIF_DB.items():
            test_path = pkg_dir / db_file
            if test_path.exists():
                return str(test_path)

    return None


def create_yeast_cpa_motif_database(output_path: str) -> str:
    """
    Create a MEME-format motif database for yeast CPA-related factors.

    This includes known binding sites for:
    - Poly(A) signal recognition factors (Hrp1, Rna15)
    - U-rich element binding (Fip1)
    - NNS termination pathway (Nrd1, Nab3)
    - Other relevant RNA-binding proteins

    Args:
        output_path: Path to write the MEME database

    Returns:
        Path to created database
    """
    meme_content = """MEME version 5.0

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.31 C 0.19 G 0.19 T 0.31

MOTIF Hrp1_polyA_signal
letter-probability matrix: alength= 4 w= 6
 0.1  0.1  0.1  0.7
 0.8  0.1  0.0  0.1
 0.1  0.1  0.1  0.7
 0.8  0.1  0.0  0.1
 0.1  0.1  0.1  0.7
 0.8  0.1  0.0  0.1

MOTIF Rna15_AAUAAA
letter-probability matrix: alength= 4 w= 6
 0.9  0.0  0.05 0.05
 0.9  0.0  0.05 0.05
 0.05 0.05 0.0  0.9
 0.9  0.0  0.05 0.05
 0.9  0.0  0.05 0.05
 0.9  0.0  0.05 0.05

MOTIF Fip1_Urich
letter-probability matrix: alength= 4 w= 8
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7

MOTIF Nrd1_GUAA
letter-probability matrix: alength= 4 w= 4
 0.1  0.1  0.7  0.1
 0.1  0.1  0.1  0.7
 0.8  0.1  0.0  0.1
 0.8  0.1  0.0  0.1

MOTIF Nab3_UCUU
letter-probability matrix: alength= 4 w= 4
 0.1  0.1  0.1  0.7
 0.1  0.7  0.1  0.1
 0.1  0.1  0.1  0.7
 0.1  0.1  0.1  0.7

MOTIF Pab1_polyA
letter-probability matrix: alength= 4 w= 10
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05

MOTIF T_tract_terminator
letter-probability matrix: alength= 4 w= 8
 0.05 0.05 0.05 0.85
 0.05 0.05 0.05 0.85
 0.05 0.05 0.05 0.85
 0.05 0.05 0.05 0.85
 0.05 0.05 0.05 0.85
 0.05 0.05 0.05 0.85
 0.05 0.05 0.05 0.85
 0.05 0.05 0.05 0.85

MOTIF Efficiency_element_EE
letter-probability matrix: alength= 4 w= 6
 0.05 0.05 0.05 0.85
 0.85 0.05 0.05 0.05
 0.05 0.05 0.05 0.85
 0.85 0.05 0.05 0.05
 0.05 0.05 0.05 0.85
 0.85 0.05 0.05 0.05

MOTIF Positioning_element_PE
letter-probability matrix: alength= 4 w= 6
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.05 0.05 0.05 0.85
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05
 0.85 0.05 0.05 0.05

"""
    with open(output_path, 'w') as f:
        f.write(meme_content)

    return output_path


def run_tomtom_motif_comparison(
    query_motifs: str,
    target_database: str,
    output_dir: str,
    evalue_threshold: float = 10.0,
    min_overlap: int = 5,
) -> Dict:
    """
    Run Tomtom to compare discovered motifs against a database.

    Args:
        query_motifs: Path to MEME-format file with query motifs
        target_database: Path to MEME-format motif database
        output_dir: Output directory for results
        evalue_threshold: E-value threshold for matches
        min_overlap: Minimum overlap between motifs

    Returns:
        Dict with:
            - success: bool
            - matches: List of matches (if successful)
            - output_dir: Path to results
            - error: Error message (if failed)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if tomtom is available
    if not _check_tool_available('tomtom'):
        return {
            'success': False,
            'error': 'tomtom not found in PATH. Install MEME Suite.',
        }

    # Run Tomtom
    cmd = [
        'tomtom',
        '-oc', str(output_dir),
        '-evalue',
        '-thresh', str(evalue_threshold),
        '-min-overlap', str(min_overlap),
        '-dist', 'pearson',
        query_motifs,
        target_database,
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            timeout=300,  # 5 minute timeout
        )

        if result.returncode != 0:
            return {
                'success': False,
                'error': f'tomtom failed: {result.stderr.decode()[:500]}',
                'output_dir': str(output_dir),
            }

        # Parse results
        matches = _parse_tomtom_results(output_dir / 'tomtom.tsv')

        return {
            'success': True,
            'matches': matches,
            'output_dir': str(output_dir),
            'n_matches': len(matches),
        }

    except subprocess.TimeoutExpired:
        return {
            'success': False,
            'error': 'tomtom timed out',
            'output_dir': str(output_dir),
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e),
            'output_dir': str(output_dir),
        }


def _parse_tomtom_results(tsv_path: Path) -> List[Dict]:
    """Parse Tomtom TSV results file."""
    matches = []

    if not tsv_path.exists():
        return matches

    try:
        df = pd.read_csv(tsv_path, sep='\t', comment='#')

        # Standardize column names
        col_map = {
            'Query_ID': 'query_motif',
            'Target_ID': 'target_motif',
            'Optimal_offset': 'offset',
            'p-value': 'pvalue',
            'E-value': 'evalue',
            'q-value': 'qvalue',
            'Overlap': 'overlap',
            'Query_consensus': 'query_consensus',
            'Target_consensus': 'target_consensus',
        }

        df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

        for _, row in df.iterrows():
            matches.append(row.to_dict())

    except Exception:
        pass

    return matches


def annotate_motifs_with_database(
    motif_results_dir: str,
    motif_database: Optional[str] = None,
    organism: str = 'scerevisiae',
    evalue_threshold: float = 0.1,
) -> pd.DataFrame:
    """
    Annotate discovered motifs by comparing to known motif databases.

    Args:
        motif_results_dir: Directory containing MEME/STREME results
        motif_database: Path to custom motif database (MEME format).
            If None, uses bundled database for organism.
        organism: Organism name for bundled database
        evalue_threshold: E-value threshold for matches

    Returns:
        DataFrame with motif annotations:
            - motif_id: Discovered motif ID
            - consensus: Discovered motif consensus
            - match_tf: Matched transcription factor
            - match_evalue: E-value of match
            - match_consensus: Consensus of matched motif
    """
    results_dir = Path(motif_results_dir)

    # Find query motifs file
    meme_file = results_dir / 'meme.txt'
    streme_file = results_dir / 'streme.txt'

    if meme_file.exists():
        query_file = meme_file
    elif streme_file.exists():
        query_file = streme_file
    else:
        print(f"Warning: No motif file found in {results_dir}")
        return pd.DataFrame()

    # Get motif database
    if motif_database is None:
        motif_database = get_bundled_motif_database(organism)

        # If no bundled database, create CPA-specific one
        if motif_database is None and organism == 'scerevisiae':
            motif_database = str(results_dir / 'yeast_cpa_motifs.meme')
            create_yeast_cpa_motif_database(motif_database)

    if motif_database is None:
        print(f"Warning: No motif database available for {organism}")
        return pd.DataFrame()

    # Run Tomtom
    tomtom_dir = results_dir / 'tomtom_results'
    result = run_tomtom_motif_comparison(
        str(query_file),
        motif_database,
        str(tomtom_dir),
        evalue_threshold=evalue_threshold,
    )

    if not result['success']:
        print(f"Warning: Tomtom failed: {result.get('error', 'Unknown error')}")
        return pd.DataFrame()

    # Format results
    matches = result.get('matches', [])
    if not matches:
        return pd.DataFrame()

    df = pd.DataFrame(matches)

    # Select relevant columns
    output_cols = ['query_motif', 'query_consensus', 'target_motif',
                   'target_consensus', 'evalue', 'pvalue', 'overlap']
    output_cols = [c for c in output_cols if c in df.columns]

    return df[output_cols]


def run_motif_discovery_with_annotation(
    foreground_sequences: List[Dict],
    background_sequences: Optional[List[Dict]] = None,
    output_dir: str = '.',
    organism: str = 'scerevisiae',
    motif_database: Optional[str] = None,
    **kwargs,
) -> Dict:
    """
    Run motif discovery and automatically annotate with known motif database.

    This is the recommended entry point for motif analysis in rectify.

    Args:
        foreground_sequences: Foreground sequences
        background_sequences: Background sequences (optional)
        output_dir: Output directory
        organism: Organism for motif database
        motif_database: Custom motif database path (optional)
        **kwargs: Additional arguments passed to run_motif_discovery

    Returns:
        Dict with motif discovery results plus annotations
    """
    # Run motif discovery
    result = run_motif_discovery(
        foreground_sequences,
        background_sequences,
        output_dir,
        **kwargs,
    )

    if not result.get('success', False):
        return result

    # Annotate motifs
    annotations = annotate_motifs_with_database(
        result['output_dir'],
        motif_database=motif_database,
        organism=organism,
    )

    result['annotations'] = annotations
    result['n_annotated'] = len(annotations)

    # Summarize top matches
    if not annotations.empty and 'target_motif' in annotations.columns:
        top_matches = annotations.groupby('query_motif')['target_motif'].first().to_dict()
        result['top_matches'] = top_matches

    return result


# =============================================================================
# Motif logos (2026-10-13: the deliverable spec asks for plots of the discovered motifs; until now only the
# STREME/MEME HTML reports existed)
# =============================================================================

LOGO_COLORS = {'A': '#2a9d3f', 'C': '#1f5fa8', 'G': '#e08a00', 'T': '#c1272d'}


def parse_meme_minimal(path: str) -> List[Dict]:
    """Motifs from a MEME-format text file (streme.txt / meme.txt), in file order.

    Each: ``{'name', 'alt', 'evalue', 'nsites', 'pwm'}`` with ``pwm`` a list of ``[pA, pC, pG, pT]`` rows.
    """
    motifs: List[Dict] = []
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('MOTIF '):
            parts = line.split()
            cur = {'name': parts[1], 'alt': parts[2] if len(parts) > 2 else '', 'evalue': None,
                   'nsites': None, 'pwm': []}
            j = i + 1
            while j < len(lines) and not lines[j].strip().startswith('letter-probability matrix'):
                if lines[j].strip().startswith('MOTIF '):
                    break
                j += 1
            if j < len(lines) and lines[j].strip().startswith('letter-probability matrix'):
                header = lines[j]
                w = None
                for key, val in (kv.split('=') for kv in
                                 header.split(':', 1)[1].replace('= ', '=').split() if '=' in kv):
                    if key == 'w':
                        w = int(float(val))
                    elif key == 'nsites':
                        cur['nsites'] = int(float(val))
                    elif key in ('E', 'S') and cur['evalue'] is None:   # MEME writes E=, STREME S= or E=
                        cur['evalue'] = float(val)
                k = j + 1
                while k < len(lines) and (w is None or len(cur['pwm']) < w):
                    row = lines[k].split()
                    if len(row) != 4:
                        break
                    try:
                        cur['pwm'].append([float(x) for x in row])
                    except ValueError:
                        break
                    k += 1
                i = k
            else:
                i = j
            if cur['pwm']:
                motifs.append(cur)
            continue
        i += 1
    return motifs


def _draw_logo(ax, pwm: List[List[float]]) -> None:
    """Information-content logo of one PWM (bits; letters stacked smallest first)."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D

    font = FontProperties(family='DejaVu Sans', weight='bold')
    for x, row in enumerate(pwm):
        total = sum(row) or 1.0
        probs = [p / total for p in row]
        ic = 2.0 + sum(p * np.log2(p) for p in probs if p > 0)
        y = 0.0
        for letter, p in sorted(zip('ACGT', probs), key=lambda t: t[1]):
            h = p * ic
            if h <= 0:
                continue
            tp = TextPath((0, 0), letter, size=1, prop=font)
            ext = tp.get_extents()
            sx = 0.9 / max(ext.width, 1e-6)
            sy = h / max(ext.height, 1e-6)
            tf = (Affine2D().translate(-ext.x0, -ext.y0).scale(sx, sy).translate(x + 0.05, y))
            ax.add_patch(PathPatch(tf.transform_path(tp), facecolor=LOGO_COLORS[letter], edgecolor='none'))
            y += h
    ax.set_xlim(0, max(len(pwm), 1))
    ax.set_ylim(0, 2)
    ax.set_xticks([i + 0.5 for i in range(len(pwm))])
    ax.set_xticklabels([str(i + 1) for i in range(len(pwm))], fontsize=7)
    ax.set_ylabel('bits')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def plot_motif_logos(results: Dict[str, Dict], output_path: str, condition: str = '') -> Optional[str]:
    """One logo per analysis: the top motif STREME (or MEME) reported, in a row of panels.

    Returns the PNG path, or None when no analysis produced a motif (nothing is written).
    """
    tops = []
    for name in ('enriched_upstream', 'enriched_downstream', 'depleted_upstream', 'depleted_downstream'):
        res = results.get(name) or {}
        out = res.get('output_dir')
        if not res.get('success') or not out:
            continue
        for fname in ('streme.txt', 'meme.txt'):
            p = Path(out) / fname
            if p.exists():
                motifs = parse_meme_minimal(str(p))
                if motifs:
                    tops.append((name, motifs[0], p))
                break
    if not tops:
        return None
    import datetime
    import matplotlib
    matplotlib.use('Agg', force=False)
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(tops), figsize=(3.2 * len(tops), 2.4), squeeze=False)
    for ax, (name, motif, src) in zip(axes[0], tops):
        _draw_logo(ax, motif['pwm'])
        e = motif.get('evalue')
        ax.set_xlabel(f"{name.replace('_', ' ')}: {motif.get('alt') or motif['name']}"
                      + (f" (E = {e:.1e})" if e is not None else ''), fontsize=8)
    fig.text(0.01, 0.01, f"rectify analyze motif discovery {condition} · top motif per scan from STREME/MEME output "
             f"· {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}", fontsize=5.5, color='0.45')
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return output_path
