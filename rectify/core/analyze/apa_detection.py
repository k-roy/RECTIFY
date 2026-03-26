#!/usr/bin/env python3
"""
Alternative Polyadenylation (APA) Isoform Detection for RECTIFY

Implements Isosceles-style APA isoform detection by grouping reads based on
their 3' end (TES - Transcript End Site) and junction pattern. This enables
quantification of alternative polyadenylation usage per gene.

Key Features:
- Groups reads by (gene, junction_signature, 3' cluster)
- Quantifies APA usage as fraction per gene
- Detects condition-specific 3' UTR changes (proximal vs distal)

Biology Context:
- Alternative polyadenylation is a major gene regulation mechanism
- Stress responses often involve 3' UTR shortening
- Different APA isoforms can have different mRNA stability/localization

Reference:
    Isosceles: Accurate long-read transcript discovery and quantification
    (Nature Communications, 2024)

Author: Kevin R. Roy
Email: kevinrjroy@gmail.com
Date: 2026-03-25
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import numpy as np

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

# Import UnifiedReadRecord for type hints
from rectify.core.unified_record import UnifiedReadRecord


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class APAIsoform:
    """
    A transcript isoform defined by its 3' end (TES) and junction pattern.

    Attributes:
        isoform_id: Unique identifier for this isoform
        gene_id: Gene identifier (systematic name)
        gene_name: Gene common name
        chrom: Chromosome
        strand: Strand (+/-)
        tes_position: Transcript End Site (3' end) position
        tes_cluster_id: ID of the CPA cluster this TES belongs to
        junction_signature: Sorted tuple of (donor, acceptor) junctions
        supporting_reads: List of read IDs supporting this isoform
        is_canonical_tes: True if TES matches annotated gene end
    """
    isoform_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    tes_position: int
    tes_cluster_id: Optional[str] = None
    junction_signature: Tuple[Tuple[int, int], ...] = field(default_factory=tuple)
    supporting_reads: List[str] = field(default_factory=list)
    is_canonical_tes: bool = False

    @property
    def n_reads(self) -> int:
        """Number of reads supporting this isoform."""
        return len(self.supporting_reads)

    @property
    def n_junctions(self) -> int:
        """Number of splice junctions."""
        return len(self.junction_signature)

    @property
    def is_intronless(self) -> bool:
        """True if isoform has no splice junctions."""
        return self.n_junctions == 0

    @property
    def junction_str(self) -> str:
        """String representation of junction signature."""
        if not self.junction_signature:
            return "intronless"
        return ";".join(f"{d}-{a}" for d, a in self.junction_signature)

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame construction."""
        return {
            'isoform_id': self.isoform_id,
            'gene_id': self.gene_id,
            'gene_name': self.gene_name,
            'chrom': self.chrom,
            'strand': self.strand,
            'tes_position': self.tes_position,
            'tes_cluster_id': self.tes_cluster_id,
            'junction_signature': self.junction_str,
            'n_junctions': self.n_junctions,
            'n_reads': self.n_reads,
            'is_canonical_tes': self.is_canonical_tes,
        }


@dataclass
class GeneAPAProfile:
    """
    APA profile for a single gene.

    Captures all APA isoforms and their relative usage.

    Attributes:
        gene_id: Gene identifier
        gene_name: Gene common name
        chrom: Chromosome
        strand: Strand
        isoforms: List of APAIsoform objects
        total_reads: Total reads attributed to this gene
    """
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    isoforms: List[APAIsoform] = field(default_factory=list)
    total_reads: int = 0

    @property
    def n_isoforms(self) -> int:
        """Number of APA isoforms."""
        return len(self.isoforms)

    @property
    def has_apa(self) -> bool:
        """True if gene has multiple TES positions."""
        tes_positions = {iso.tes_position for iso in self.isoforms}
        return len(tes_positions) > 1

    def get_isoform_fractions(self) -> Dict[str, float]:
        """Get fraction of reads per isoform."""
        if self.total_reads == 0:
            return {}
        return {
            iso.isoform_id: iso.n_reads / self.total_reads
            for iso in self.isoforms
        }

    def get_tes_usage(self) -> Dict[int, float]:
        """Get fraction of reads per TES position."""
        if self.total_reads == 0:
            return {}

        tes_counts = defaultdict(int)
        for iso in self.isoforms:
            tes_counts[iso.tes_position] += iso.n_reads

        return {
            tes: count / self.total_reads
            for tes, count in tes_counts.items()
        }


# =============================================================================
# CPA Cluster Assignment
# =============================================================================

def assign_reads_to_clusters(
    records: List[UnifiedReadRecord],
    cluster_df: 'pd.DataFrame',
    position_col: str = 'three_prime_corrected',
) -> Dict[str, Tuple[str, int]]:
    """
    Assign reads to CPA clusters based on their 3' end position.

    Args:
        records: List of UnifiedReadRecord objects
        cluster_df: DataFrame with cluster definitions (from clustering module)
        position_col: Column to use for read 3' position

    Returns:
        Dictionary mapping read_id to (cluster_id, cluster_modal_position)
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for cluster assignment")

    # Build cluster lookup by chrom/strand
    cluster_lookup = defaultdict(list)
    for _, cluster in cluster_df.iterrows():
        key = (cluster['chrom'], cluster['strand'])
        cluster_lookup[key].append({
            'cluster_id': cluster['cluster_id'],
            'start': cluster['start'],
            'end': cluster['end'],
            'modal': cluster['modal_position'],
        })

    # Sort clusters by position for efficient lookup
    for key in cluster_lookup:
        cluster_lookup[key].sort(key=lambda x: x['start'])

    # Assign reads to clusters
    assignments = {}

    for record in records:
        key = (record.chrom, record.strand)
        if key not in cluster_lookup:
            continue

        # Get 3' position
        pos_3prime = record.three_prime_corrected

        # Find overlapping cluster
        for cluster in cluster_lookup[key]:
            if cluster['start'] <= pos_3prime <= cluster['end']:
                assignments[record.read_id] = (
                    cluster['cluster_id'],
                    cluster['modal']
                )
                break

    return assignments


# =============================================================================
# APA Isoform Detection
# =============================================================================

def detect_apa_isoforms(
    records: List[UnifiedReadRecord],
    gene_attributions: Dict[str, Tuple[str, str]],
    cluster_assignments: Optional[Dict[str, Tuple[str, int]]] = None,
    min_reads_per_isoform: int = 3,
    annotated_tes: Optional[Dict[str, List[int]]] = None,
    tes_tolerance: int = 50,
) -> List[APAIsoform]:
    """
    Detect APA isoforms by grouping reads by gene, junction pattern, and 3' cluster.

    Isosceles insight: Group reads by (gene, junction_signature, TES_cluster)
    to define distinct transcript isoforms. This enables:
    - Quantification of APA usage per gene
    - Detection of condition-specific 3' UTR shortening/lengthening
    - Linking junction patterns to specific 3' ends

    Args:
        records: List of UnifiedReadRecord objects
        gene_attributions: Dict mapping read_id to (gene_id, gene_name)
        cluster_assignments: Dict mapping read_id to (cluster_id, cluster_modal)
        min_reads_per_isoform: Minimum reads to call an isoform
        annotated_tes: Dict mapping gene_id to list of annotated TES positions
        tes_tolerance: Distance tolerance for matching to annotated TES

    Returns:
        List of APAIsoform objects
    """
    if annotated_tes is None:
        annotated_tes = {}

    # Group reads by (gene, junction_signature, tes_cluster/position)
    isoform_groups = defaultdict(lambda: {
        'reads': [],
        'chrom': None,
        'strand': None,
        'tes_positions': [],
    })

    for record in records:
        read_id = record.read_id

        # Get gene attribution
        if read_id not in gene_attributions:
            continue
        gene_id, gene_name = gene_attributions[read_id]

        # Get junction signature
        junction_sig = tuple(sorted(record.junctions))

        # Get TES (cluster or raw position)
        if cluster_assignments and read_id in cluster_assignments:
            cluster_id, tes_modal = cluster_assignments[read_id]
            tes_key = cluster_id
        else:
            tes_modal = record.three_prime_corrected
            tes_key = f"pos_{tes_modal}"

        # Group key
        key = (gene_id, gene_name, junction_sig, tes_key)

        isoform_groups[key]['reads'].append(read_id)
        isoform_groups[key]['chrom'] = record.chrom
        isoform_groups[key]['strand'] = record.strand
        isoform_groups[key]['tes_positions'].append(tes_modal)

    # Create APAIsoform objects
    isoforms = []
    isoform_counter = 0

    for (gene_id, gene_name, junction_sig, tes_key), data in isoform_groups.items():
        n_reads = len(data['reads'])
        if n_reads < min_reads_per_isoform:
            continue

        # Calculate modal TES position
        tes_position = int(np.median(data['tes_positions']))

        # Check if canonical TES
        is_canonical = False
        if gene_id in annotated_tes:
            for ann_tes in annotated_tes[gene_id]:
                if abs(tes_position - ann_tes) <= tes_tolerance:
                    is_canonical = True
                    break

        # Create isoform
        isoform_id = f"APA_{gene_id}_{isoform_counter:03d}"
        isoform_counter += 1

        isoforms.append(APAIsoform(
            isoform_id=isoform_id,
            gene_id=gene_id,
            gene_name=gene_name,
            chrom=data['chrom'],
            strand=data['strand'],
            tes_position=tes_position,
            tes_cluster_id=tes_key if tes_key.startswith('cluster') else None,
            junction_signature=junction_sig,
            supporting_reads=data['reads'],
            is_canonical_tes=is_canonical,
        ))

    return isoforms


def build_gene_apa_profiles(
    isoforms: List[APAIsoform],
) -> Dict[str, GeneAPAProfile]:
    """
    Build APA profiles for each gene.

    Groups isoforms by gene and calculates total reads.

    Args:
        isoforms: List of APAIsoform objects

    Returns:
        Dictionary mapping gene_id to GeneAPAProfile
    """
    profiles = {}

    for iso in isoforms:
        if iso.gene_id not in profiles:
            profiles[iso.gene_id] = GeneAPAProfile(
                gene_id=iso.gene_id,
                gene_name=iso.gene_name,
                chrom=iso.chrom,
                strand=iso.strand,
                isoforms=[],
                total_reads=0,
            )

        profiles[iso.gene_id].isoforms.append(iso)
        profiles[iso.gene_id].total_reads += iso.n_reads

    return profiles


# =============================================================================
# APA Quantification
# =============================================================================

def quantify_apa_usage(
    isoforms: List[APAIsoform],
    by_gene: bool = True,
) -> 'pd.DataFrame':
    """
    Quantify APA usage as fraction of reads per TES per gene.

    Args:
        isoforms: List of APAIsoform objects
        by_gene: If True, calculate fractions within each gene

    Returns:
        DataFrame with columns:
        - gene_id, gene_name
        - tes_position, junction_signature
        - n_reads, fraction
        - is_canonical
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for APA quantification")

    if not isoforms:
        return pd.DataFrame(columns=[
            'gene_id', 'gene_name', 'chrom', 'strand',
            'tes_position', 'junction_signature', 'n_reads',
            'fraction', 'is_canonical'
        ])

    # Build profiles
    profiles = build_gene_apa_profiles(isoforms)

    # Create output rows
    rows = []
    for gene_id, profile in profiles.items():
        for iso in profile.isoforms:
            fraction = iso.n_reads / profile.total_reads if profile.total_reads > 0 else 0

            rows.append({
                'gene_id': iso.gene_id,
                'gene_name': iso.gene_name,
                'chrom': iso.chrom,
                'strand': iso.strand,
                'tes_position': iso.tes_position,
                'junction_signature': iso.junction_str,
                'n_junctions': iso.n_junctions,
                'n_reads': iso.n_reads,
                'fraction': fraction,
                'is_canonical': iso.is_canonical_tes,
            })

    df = pd.DataFrame(rows)

    # Sort by gene and fraction
    df = df.sort_values(
        ['gene_id', 'fraction'],
        ascending=[True, False]
    ).reset_index(drop=True)

    return df


def identify_proximal_distal_tes(
    profiles: Dict[str, GeneAPAProfile],
    min_tes_difference: int = 100,
) -> 'pd.DataFrame':
    """
    Identify proximal (shorter 3' UTR) vs distal (longer 3' UTR) TES usage.

    For genes with multiple TES, categorize as proximal/distal based on
    genomic position and strand.

    Args:
        profiles: Dictionary of GeneAPAProfile objects
        min_tes_difference: Minimum distance between TES to classify

    Returns:
        DataFrame with proximal/distal TES fractions per gene
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required for TES analysis")

    rows = []

    for gene_id, profile in profiles.items():
        if not profile.has_apa:
            continue

        tes_usage = profile.get_tes_usage()
        if len(tes_usage) < 2:
            continue

        # Sort TES positions
        sorted_tes = sorted(tes_usage.keys())

        # Check if difference is significant
        if sorted_tes[-1] - sorted_tes[0] < min_tes_difference:
            continue

        # Determine proximal vs distal based on strand
        # Plus strand: proximal = lower position (closer to CDS)
        # Minus strand: proximal = higher position (closer to CDS)
        if profile.strand == '+':
            proximal_tes = sorted_tes[0]
            distal_tes = sorted_tes[-1]
        else:
            proximal_tes = sorted_tes[-1]
            distal_tes = sorted_tes[0]

        proximal_fraction = tes_usage[proximal_tes]
        distal_fraction = tes_usage[distal_tes]

        rows.append({
            'gene_id': gene_id,
            'gene_name': profile.gene_name,
            'chrom': profile.chrom,
            'strand': profile.strand,
            'n_tes': len(tes_usage),
            'proximal_tes': proximal_tes,
            'distal_tes': distal_tes,
            'proximal_fraction': proximal_fraction,
            'distal_fraction': distal_fraction,
            'tes_distance': abs(distal_tes - proximal_tes),
            'total_reads': profile.total_reads,
        })

    return pd.DataFrame(rows)


# =============================================================================
# Annotation Loading
# =============================================================================

def load_annotated_tes_from_gff(
    gff_path: str,
    feature_type: str = 'gene',
) -> Dict[str, List[int]]:
    """
    Load annotated transcript end sites (TES) from GFF.

    TES is the 3' end of the gene:
    - Plus strand: gene end (rightmost coordinate)
    - Minus strand: gene start (leftmost coordinate)

    Args:
        gff_path: Path to GFF/GTF file
        feature_type: Feature type to extract TES from

    Returns:
        Dictionary mapping gene_id to list of TES positions
    """
    tes_positions = defaultdict(list)

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == feature_type:
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]

                # Extract gene ID
                gene_id = None
                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1]
                        break
                    elif attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1] if '"' in attr else attr.split('=')[1]
                        break

                if gene_id:
                    # TES is 3' end
                    if strand == '+':
                        tes = end - 1
                    else:
                        tes = start

                    tes_positions[gene_id].append(tes)

    return dict(tes_positions)


# =============================================================================
# Summary Statistics
# =============================================================================

def summarize_apa_detection(
    isoforms: List[APAIsoform],
    profiles: Dict[str, GeneAPAProfile],
) -> Dict:
    """
    Generate summary statistics for APA detection.

    Args:
        isoforms: List of detected APA isoforms
        profiles: Dictionary of gene APA profiles

    Returns:
        Dictionary with summary statistics
    """
    n_isoforms = len(isoforms)
    n_genes = len(profiles)
    n_genes_with_apa = sum(1 for p in profiles.values() if p.has_apa)

    # Count by junction type
    n_intronless = sum(1 for iso in isoforms if iso.is_intronless)
    n_spliced = n_isoforms - n_intronless

    # Count canonical vs novel
    n_canonical = sum(1 for iso in isoforms if iso.is_canonical_tes)
    n_novel = n_isoforms - n_canonical

    # TES statistics
    tes_per_gene = [p.n_isoforms for p in profiles.values()]

    return {
        'total_isoforms': n_isoforms,
        'total_genes': n_genes,
        'genes_with_apa': n_genes_with_apa,
        'genes_single_tes': n_genes - n_genes_with_apa,
        'intronless_isoforms': n_intronless,
        'spliced_isoforms': n_spliced,
        'canonical_tes': n_canonical,
        'novel_tes': n_novel,
        'mean_isoforms_per_gene': np.mean(tes_per_gene) if tes_per_gene else 0,
        'max_isoforms_per_gene': max(tes_per_gene) if tes_per_gene else 0,
    }


# =============================================================================
# DataFrame Export
# =============================================================================

def isoforms_to_dataframe(
    isoforms: List[APAIsoform],
) -> 'pd.DataFrame':
    """
    Convert APA isoforms to pandas DataFrame.

    Args:
        isoforms: List of APAIsoform objects

    Returns:
        DataFrame with isoform information
    """
    if not HAS_PANDAS:
        raise ImportError("pandas is required")

    if not isoforms:
        return pd.DataFrame(columns=[
            'isoform_id', 'gene_id', 'gene_name', 'chrom', 'strand',
            'tes_position', 'tes_cluster_id', 'junction_signature',
            'n_junctions', 'n_reads', 'is_canonical_tes'
        ])

    rows = [iso.to_dict() for iso in isoforms]
    return pd.DataFrame(rows)


if __name__ == '__main__':
    # Quick test
    print("Testing apa_detection module...")

    # Create test records with different junction patterns and 3' ends
    test_records = [
        # Gene A, junction pattern 1, TES position 1
        UnifiedReadRecord(
            read_id="read_001", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=3000, three_prime_corrected=2995,
            alignment_start=1000, alignment_end=3000,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_002", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=3000, three_prime_corrected=2998,
            alignment_start=1000, alignment_end=3000,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_003", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=3000, three_prime_corrected=3002,
            alignment_start=1000, alignment_end=3000,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        # Gene A, junction pattern 1, TES position 2 (APA)
        UnifiedReadRecord(
            read_id="read_004", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=2500, three_prime_corrected=2495,
            alignment_start=1000, alignment_end=2500,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_005", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=2500, three_prime_corrected=2498,
            alignment_start=1000, alignment_end=2500,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
        UnifiedReadRecord(
            read_id="read_006", chrom="chrI", strand="+",
            five_prime_raw=1000, five_prime_corrected=1000,
            first_exon_start=1000, starts_in_intron=False,
            three_prime_raw=2500, three_prime_corrected=2502,
            alignment_start=1000, alignment_end=2500,
            junctions=[(1500, 1700)], n_junctions=1,
        ),
    ]

    # Create gene attributions
    gene_attributions = {
        'read_001': ('YAL001C', 'TFC3'),
        'read_002': ('YAL001C', 'TFC3'),
        'read_003': ('YAL001C', 'TFC3'),
        'read_004': ('YAL001C', 'TFC3'),
        'read_005': ('YAL001C', 'TFC3'),
        'read_006': ('YAL001C', 'TFC3'),
    }

    # Detect isoforms
    isoforms = detect_apa_isoforms(
        test_records,
        gene_attributions,
        min_reads_per_isoform=3,
    )

    print(f"\nDetected {len(isoforms)} APA isoforms:")
    for iso in isoforms:
        print(f"  {iso.isoform_id}: TES={iso.tes_position}, "
              f"junctions={iso.junction_str}, reads={iso.n_reads}")

    # Build profiles
    profiles = build_gene_apa_profiles(isoforms)

    print(f"\nGene APA profiles:")
    for gene_id, profile in profiles.items():
        print(f"  {gene_id} ({profile.gene_name}): "
              f"{profile.n_isoforms} isoforms, "
              f"has_apa={profile.has_apa}")
        for tes, fraction in profile.get_tes_usage().items():
            print(f"    TES {tes}: {fraction:.1%}")

    # Quantify
    if HAS_PANDAS:
        df = quantify_apa_usage(isoforms)
        print(f"\nAPA quantification DataFrame:")
        print(df.to_string())

    # Summary
    summary = summarize_apa_detection(isoforms, profiles)
    print(f"\nSummary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")

    print("\napa_detection module ready!")
