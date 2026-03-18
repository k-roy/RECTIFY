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
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

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
                except Exception:
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
                except Exception:
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
) -> Dict[str, Dict]:
    """
    Run motif discovery for enriched vs depleted clusters.

    Performs 4 analyses IN PARALLEL:
    1. Enriched upstream vs background
    2. Enriched downstream vs background
    3. Depleted upstream vs background
    4. Depleted downstream vs background

    Args:
        enriched_clusters: Clusters with increased signal
        depleted_clusters: Clusters with decreased signal
        genome_fasta: Path to genome FASTA
        output_dir: Output directory
        upstream_window: bp upstream (default: 100)
        downstream_window: bp downstream (default: 50)
        n_parallel: Number of parallel workers (default: 4)

    Returns:
        Dict mapping analysis_name -> motif discovery results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    # Extract sequences
    print("Extracting sequences around clusters...")
    enriched_seqs = extract_sequences_around_clusters(
        enriched_clusters, genome_fasta,
        upstream_window=upstream_window,
        downstream_window=downstream_window,
    )

    depleted_seqs = extract_sequences_around_clusters(
        depleted_clusters, genome_fasta,
        upstream_window=upstream_window,
        downstream_window=downstream_window,
    )

    # Prepare analyses
    analyses = [
        ('enriched_upstream', enriched_seqs.get('upstream', []), depleted_seqs.get('upstream', []), str(output_dir)),
        ('enriched_downstream', enriched_seqs.get('downstream', []), depleted_seqs.get('downstream', []), str(output_dir)),
        ('depleted_upstream', depleted_seqs.get('upstream', []), enriched_seqs.get('upstream', []), str(output_dir)),
        ('depleted_downstream', depleted_seqs.get('downstream', []), enriched_seqs.get('downstream', []), str(output_dir)),
    ]

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
