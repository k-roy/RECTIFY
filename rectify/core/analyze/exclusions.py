"""
Lightweight exclusion-region detection for the analyze pipeline.

Returns (mito_chroms, rdna_regions) tuples consumed by the boolean-mask
filter inline in `run_analyze` / `_run_analyze_manifest`. A richer
class-based detector lives in `rectify/core/exclusion_regions.py`
(`ExclusionRegionDetector`) — the two could potentially be consolidated
behind a single API in a future pass, but the call sites differ in shape
enough that they have not been merged yet.
"""

import pandas as pd


def _detect_exclusion_regions(
    annotation_df: pd.DataFrame,
    data_chroms: list,
) -> tuple:
    """
    Auto-detect mitochondrial chromosomes and rDNA regions from annotation.

    Uses annotation to identify:
    1. Mitochondrial chromosomes (by name or gene content)
    2. rDNA loci (by gene names like RDN*, ribosomal RNA genes)

    Args:
        annotation_df: Gene annotation DataFrame (may be None)
        data_chroms: List of chromosome names in the data

    Returns:
        Tuple of (mito_chroms: set, rdna_regions: list of (chrom, start, end))
    """
    # Default mitochondrial chromosome patterns
    mito_patterns = {'chrM', 'chrMT', 'chrmt', 'MT', 'Mt', 'Mito', 'mitochondrion'}
    mito_ncbi = {'ref|NC_001224|'}  # Yeast mito

    # Find mitochondrial chromosomes in data
    mito_chroms = set()
    for chrom in data_chroms:
        chrom_lower = str(chrom).lower()
        if any(pat.lower() in chrom_lower for pat in mito_patterns):
            mito_chroms.add(chrom)
        if chrom in mito_ncbi:
            mito_chroms.add(chrom)

    # rDNA regions - detect from annotation if available
    rdna_regions = []

    if annotation_df is not None and not annotation_df.empty:
        # Look for mitochondrial genes in annotation
        if 'chrom' in annotation_df.columns:
            for chrom in annotation_df['chrom'].unique():
                chrom_lower = str(chrom).lower()
                if any(pat.lower() in chrom_lower for pat in mito_patterns):
                    mito_chroms.add(chrom)

        # Look for rDNA genes (yeast: RDN5, RDN18, RDN25, RDN37, RDN58, ETS, ITS, NTS)
        # These are the actual ribosomal RNA genes in the rDNA repeat unit
        # Pattern matches gene names starting with RDN, ETS, ITS, NTS (rDNA-specific)
        rdna_gene_patterns = ['^RDN', '^ETS', '^ITS', '^NTS']

        gene_col = 'gene_name' if 'gene_name' in annotation_df.columns else 'gene_id'
        if gene_col in annotation_df.columns:
            combined_pattern = '|'.join(rdna_gene_patterns)
            rdna_genes = annotation_df[
                annotation_df[gene_col].str.match(combined_pattern, case=False, na=False)
            ]

            if not rdna_genes.empty and 'chrom' in rdna_genes.columns:
                # Group by chromosome and find bounds
                for chrom in rdna_genes['chrom'].unique():
                    chrom_genes = rdna_genes[rdna_genes['chrom'] == chrom]
                    if 'start' in chrom_genes.columns and 'end' in chrom_genes.columns:
                        # Get region bounds with padding
                        start = int(chrom_genes['start'].min()) - 1000
                        end = int(chrom_genes['end'].max()) + 1000
                        start = max(0, start)

                        # Add region if substantial (>5kb suggests rDNA locus)
                        if end - start > 5000:
                            rdna_regions.append((chrom, start, end))

    # Deduplicate and sort rDNA regions
    rdna_regions = list(set(rdna_regions))
    rdna_regions.sort()

    # If no rDNA detected from annotation, use default yeast coordinates
    if not rdna_regions:
        # Yeast rDNA locus on chrXII: ~451,000-469,000 bp
        yeast_rdna_chroms = {'chrXII', 'chr12', 'ref|NC_001144|'}
        for chrom in data_chroms:
            if chrom in yeast_rdna_chroms:
                rdna_regions.append((chrom, 450_000, 490_000))
                break

    return mito_chroms, rdna_regions
