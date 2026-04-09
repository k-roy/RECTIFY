"""
Exclusion region detection for RECTIFY.

Identifies and provides coordinates for regions to exclude from CPA analysis:
- rDNA locus (chrXII ~450,000-490,000 in yeast)
- Pol III genes (tRNAs, SNR6, RDN5, RPR1, SCR1) + flanking
- Mitochondrial genome (optional)

Author: Kevin R. Roy
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Dict, Set, Optional
import re

import pandas as pd


@dataclass
class ExclusionRegion:
    """Single genomic region to exclude."""
    chrom: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    reason: str  # 'rDNA', 'tRNA', 'snRNA', 'ncRNA', 'mito'
    gene_name: Optional[str] = None

    def __hash__(self):
        return hash((self.chrom, self.start, self.end, self.reason))

    def contains(self, chrom: str, position: int) -> bool:
        """Check if position falls within this region."""
        return self.chrom == chrom and self.start <= position < self.end


class ExclusionRegionDetector:
    """
    Detect and manage genomic exclusion regions.

    Usage:
        detector = ExclusionRegionDetector(flanking_bp=100)
        detector.load_from_gff(gff_path)
        detector.add_rdna_region()

        # Check if position should be excluded
        if detector.is_excluded('chrXII', 455000):
            continue

        # Or filter a DataFrame
        kept_df, excluded_df = detector.filter_dataframe(df)
    """

    # Default yeast rDNA locus coordinates
    YEAST_RDNA_CHROMS = {'chrXII', 'chr12', 'ref|NC_001144|'}
    YEAST_RDNA_START = 450_000
    YEAST_RDNA_END = 490_000

    # Pol III gene patterns
    POL3_FEATURE_TYPES = {'tRNA_gene', 'tRNA'}
    POL3_GENE_PATTERNS = [
        r'^SNR6$',   # U6 snRNA
        r'^RDN5',    # 5S rRNA genes
        r'^RPR1$',   # RNase P RNA
        r'^SCR1$',   # Signal recognition particle RNA
    ]

    # rDNA gene patterns (for auto-detection)
    RDNA_GENE_PATTERNS = [r'^RDN', r'^ETS', r'^ITS', r'^NTS']

    # Mitochondrial patterns
    MITO_PATTERNS = {'chrM', 'chrMT', 'chrmt', 'MT', 'Mt', 'Mito', 'mitochondrion'}
    MITO_NCBI = {'ref|NC_001224|'}

    def __init__(self, flanking_bp: int = 100, chrom_sizes: Optional[Dict[str, int]] = None):
        """
        Args:
            flanking_bp: Bases to add on each side of Pol III genes
            chrom_sizes: Optional dict of {chrom: size} used to clamp downstream
                         flanking to valid chromosome boundaries
        """
        self.flanking_bp = flanking_bp
        self.chrom_sizes: Optional[Dict[str, int]] = chrom_sizes
        self.regions: List[ExclusionRegion] = []
        self._region_index: Dict[str, List[Tuple[int, int, str]]] = {}

    def add_region(
        self,
        chrom: str,
        start: int,
        end: int,
        reason: str,
        gene_name: Optional[str] = None,
    ) -> None:
        """Add a single exclusion region."""
        # Clamp upstream boundary to 0 (chromosome start).
        clamped_start = max(0, start)
        # Clamp downstream boundary to chromosome size when available.
        if self.chrom_sizes and chrom in self.chrom_sizes:
            clamped_end = min(end, self.chrom_sizes[chrom])
        else:
            clamped_end = end
        region = ExclusionRegion(
            chrom=chrom,
            start=clamped_start,
            end=clamped_end,
            reason=reason,
            gene_name=gene_name,
        )
        self.regions.append(region)

        # Update index for fast lookup
        if chrom not in self._region_index:
            self._region_index[chrom] = []
        self._region_index[chrom].append((region.start, region.end, reason))

    def add_rdna_region(
        self,
        chrom: str = 'chrXII',
        start: int = None,
        end: int = None,
    ) -> None:
        """Add rDNA locus (with default yeast coordinates)."""
        if start is None:
            start = self.YEAST_RDNA_START
        if end is None:
            end = self.YEAST_RDNA_END
        self.add_region(chrom, start, end, 'rDNA', 'rDNA_locus')

    def add_mito_chrom(self, chrom: str) -> None:
        """Mark entire chromosome as mitochondrial (exclude all)."""
        # Use a very large end coordinate to cover the whole chromosome
        self.add_region(chrom, 0, 100_000_000, 'mito', 'mitochondrion')

    def load_from_gff(
        self,
        gff_path: Path,
        exclude_tRNA: bool = True,
        exclude_snRNA: bool = True,
        exclude_rDNA: bool = True,
        exclude_mito: bool = False,
        data_chroms: Optional[Set[str]] = None,
        chrom_sizes: Optional[Dict[str, int]] = None,
    ) -> int:
        """
        Auto-detect exclusion regions from GFF annotation.

        Detects:
        - tRNA_gene features (Pol III)
        - snRNA_gene features matching SNR6 (U6)
        - ncRNA_gene with RPR1, SCR1
        - rRNA_gene features with RDN5
        - rDNA locus from RDN*, ETS*, ITS*, NTS* genes

        Args:
            gff_path: Path to GFF/GFF3 annotation file
            exclude_tRNA: Include tRNA genes in exclusion regions
            exclude_snRNA: Include SNR6 (U6 snRNA) in exclusion regions
            exclude_rDNA: Include rDNA locus in exclusion regions
            exclude_mito: Include mitochondrial chromosome in exclusion regions
            data_chroms: Set of chromosome names in data (for rDNA fallback)
            chrom_sizes: Optional dict of {chrom: size} to clamp downstream
                         flanking regions to valid chromosome boundaries.
                         Stored on self.chrom_sizes for all subsequent add_region calls.

        Returns:
            Number of exclusion regions detected
        """
        import gzip

        # Update stored chrom_sizes if provided (used by add_region for clamping).
        if chrom_sizes is not None:
            self.chrom_sizes = chrom_sizes

        gff_path = Path(gff_path)
        open_fn = gzip.open if gff_path.suffix == '.gz' else open

        initial_count = len(self.regions)
        rdna_genes = []  # Collect for locus detection
        mito_chroms = set()

        with open_fn(gff_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chrom, source, feature_type, start_str, end_str, _, strand, _, attrs = fields

                try:
                    start = int(start_str) - 1  # Convert to 0-based
                    end = int(end_str)          # Already exclusive in 1-based GFF
                except ValueError:
                    continue

                # Parse GFF3 attributes
                attr_dict = {}
                for attr in attrs.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                gene_name = attr_dict.get('gene', attr_dict.get('Name', attr_dict.get('ID', '')))

                # Check for mitochondrial chromosomes
                chrom_lower = chrom.lower()
                if any(pat.lower() in chrom_lower for pat in self.MITO_PATTERNS):
                    mito_chroms.add(chrom)
                if chrom in self.MITO_NCBI:
                    mito_chroms.add(chrom)

                # Check for tRNA genes
                if exclude_tRNA and feature_type in self.POL3_FEATURE_TYPES:
                    self.add_region(
                        chrom,
                        start - self.flanking_bp,
                        end + self.flanking_bp,
                        'tRNA',
                        gene_name,
                    )
                    continue

                # Check for Pol III gene patterns (SNR6, RDN5, RPR1, SCR1)
                for pattern in self.POL3_GENE_PATTERNS:
                    if re.match(pattern, gene_name, re.IGNORECASE):
                        if 'SNR6' in gene_name.upper() and exclude_snRNA:
                            self.add_region(
                                chrom,
                                start - self.flanking_bp,
                                end + self.flanking_bp,
                                'snRNA',
                                gene_name,
                            )
                        elif 'RDN5' in gene_name.upper() and exclude_rDNA:
                            self.add_region(
                                chrom,
                                start - self.flanking_bp,
                                end + self.flanking_bp,
                                'rRNA',
                                gene_name,
                            )
                        elif gene_name.upper() in ('RPR1', 'SCR1') and exclude_snRNA:
                            self.add_region(
                                chrom,
                                start - self.flanking_bp,
                                end + self.flanking_bp,
                                'ncRNA',
                                gene_name,
                            )
                        break

                # Collect rDNA genes for locus detection
                if exclude_rDNA:
                    for pattern in self.RDNA_GENE_PATTERNS:
                        if re.match(pattern, gene_name, re.IGNORECASE):
                            rdna_genes.append((chrom, start, end, gene_name))
                            break

        # Detect rDNA locus from collected genes
        if exclude_rDNA and rdna_genes:
            self._detect_rdna_locus(rdna_genes)
        elif exclude_rDNA and data_chroms:
            # Fall back to default yeast coordinates
            for chrom in data_chroms:
                if chrom in self.YEAST_RDNA_CHROMS:
                    self.add_rdna_region(chrom)
                    break

        # Add mitochondrial exclusion if requested
        if exclude_mito:
            for chrom in mito_chroms:
                self.add_mito_chrom(chrom)

        return len(self.regions) - initial_count

    def _detect_rdna_locus(self, rdna_genes: List[Tuple[str, int, int, str]]) -> None:
        """Detect rDNA locus boundaries from collected rDNA genes."""
        if not rdna_genes:
            return

        # Group by chromosome
        by_chrom: Dict[str, List[Tuple[int, int]]] = {}
        for chrom, start, end, _ in rdna_genes:
            if chrom not in by_chrom:
                by_chrom[chrom] = []
            by_chrom[chrom].append((start, end))

        # Find locus on each chromosome
        for chrom, coords in by_chrom.items():
            if not coords:
                continue
            starts, ends = zip(*coords)
            locus_start = min(starts) - 1000  # Padding
            locus_end = max(ends) + 1000
            locus_start = max(0, locus_start)

            # Only add if substantial (>5kb suggests rDNA locus)
            if locus_end - locus_start > 5000:
                self.add_region(chrom, locus_start, locus_end, 'rDNA', 'rDNA_locus')

    def is_excluded(self, chrom: str, position: int, reason_filter: Optional[str] = None) -> bool:
        """
        Check if genomic position falls within any exclusion region.

        Args:
            chrom: Chromosome name
            position: 0-based position
            reason_filter: Optional filter for specific reason (e.g., 'rDNA', 'tRNA')

        Returns:
            True if position should be excluded
        """
        if chrom not in self._region_index:
            return False

        for start, end, reason in self._region_index[chrom]:
            if reason_filter and reason != reason_filter:
                continue
            if start <= position < end:
                return True
        return False

    def filter_dataframe(
        self,
        df: pd.DataFrame,
        chrom_col: str = 'chrom',
        position_col: str = 'three_prime_raw',
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Filter DataFrame to remove excluded positions.

        Args:
            df: DataFrame with chromosome and position columns
            chrom_col: Name of chromosome column
            position_col: Name of position column

        Returns:
            Tuple of (kept_df, excluded_df)
        """
        if df.empty or not self.regions:
            return df, pd.DataFrame(columns=df.columns)

        # Coordinate validation: if the minimum position in the column is exactly 1,
        # this is a strong indicator that 1-based (GFF/GTF) coordinates were passed
        # without the required -1 conversion. Raise immediately so the caller fixes
        # the source rather than silently filtering the wrong genomic positions.
        if position_col in df.columns and len(df) > 0:
            col_min = int(df[position_col].min())
            if col_min == 1 and (df[position_col] == 1).sum() / len(df) > 0.01:
                import logging
                logging.getLogger(__name__).error(
                    "filter_dataframe(): column '%s' has a minimum value of 1 with "
                    "%.1f%% of rows at position 1 — this strongly indicates 1-based "
                    "(GFF/GTF) coordinates were passed to a 0-based function. "
                    "Subtract 1 from '%s' before calling filter_dataframe() "
                    "(e.g. via load_annotation() which performs this conversion). "
                    "Proceeding, but exclusion filtering may target incorrect positions.",
                    position_col,
                    100.0 * (df[position_col] == 1).sum() / len(df),
                    position_col,
                )

        # Build exclusion mask
        excluded_mask = pd.Series(False, index=df.index)

        for region in self.regions:
            region_mask = (
                (df[chrom_col] == region.chrom) &
                (df[position_col] >= region.start) &
                (df[position_col] < region.end)
            )
            excluded_mask = excluded_mask | region_mask

        kept_df = df[~excluded_mask].copy()
        excluded_df = df[excluded_mask].copy()

        return kept_df, excluded_df

    def get_summary(self) -> pd.DataFrame:
        """Return DataFrame summarizing all exclusion regions."""
        if not self.regions:
            return pd.DataFrame(columns=['chrom', 'start', 'end', 'reason', 'gene_name', 'size'])

        data = []
        for r in self.regions:
            data.append({
                'chrom': r.chrom,
                'start': r.start,
                'end': r.end,
                'reason': r.reason,
                'gene_name': r.gene_name or '',
                'size': r.end - r.start,
            })

        df = pd.DataFrame(data)
        return df.sort_values(['chrom', 'start']).reset_index(drop=True)

    def get_stats_by_reason(self) -> Dict[str, int]:
        """Return count of regions by reason."""
        stats: Dict[str, int] = {}
        for region in self.regions:
            stats[region.reason] = stats.get(region.reason, 0) + 1
        return stats

    def __len__(self) -> int:
        return len(self.regions)

    def __repr__(self) -> str:
        stats = self.get_stats_by_reason()
        parts = [f"{reason}={count}" for reason, count in sorted(stats.items())]
        return f"ExclusionRegionDetector({', '.join(parts) if parts else 'empty'})"
