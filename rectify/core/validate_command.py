#!/usr/bin/env python3
"""
RECTIFY 'validate' command implementation.

Validates corrections against ground truth data. Supports multiple ground
truth sources (all optional, but at least one required):
- NET-seq BigWig files
- Gene annotations (GTF/GFF) with known 3' ends
- Known positions TSV file

Author: Kevin R. Roy
Date: 2026-03-16
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import logging
import sys

import numpy as np
from tqdm import tqdm

from ..slurm import set_thread_limits

logger = logging.getLogger(__name__)


@dataclass
class CorrectedPosition:
    """A corrected 3' end position from RECTIFY output."""
    read_id: str
    chrom: str
    strand: str
    original_3prime: int
    corrected_3prime: int
    ambiguity_min: Optional[int] = None
    ambiguity_max: Optional[int] = None
    ambiguity_range: int = 0
    correction_applied: str = ""
    confidence: str = ""
    qc_flags: str = ""


@dataclass
class ValidationResult:
    """Validation result for a single position."""
    read_id: str
    chrom: str
    strand: str
    original_3prime: int
    corrected_3prime: int
    ground_truth_position: Optional[int] = None
    ground_truth_source: str = "none"
    distance_from_truth: Optional[int] = None
    original_distance: Optional[int] = None
    is_correct_exact: bool = False
    is_correct_1bp: bool = False
    is_correct_2bp: bool = False
    is_correct_5bp: bool = False
    improvement_bp: int = 0
    validation_status: str = "NO_TRUTH"
    correction_applied: str = ""
    confidence: str = ""


@dataclass
class ValidationMetrics:
    """Aggregate validation metrics."""
    total_positions: int = 0
    with_ground_truth: int = 0

    # Corrected accuracy
    exact_match: int = 0
    within_1bp: int = 0
    within_2bp: int = 0
    within_5bp: int = 0

    # Original accuracy (for comparison)
    original_exact_match: int = 0
    original_within_1bp: int = 0
    original_within_2bp: int = 0

    # Improvement
    improved: int = 0
    unchanged: int = 0
    worsened: int = 0

    # By correction type
    by_correction: Dict[str, Dict] = field(default_factory=dict)

    # By confidence level
    by_confidence: Dict[str, Dict] = field(default_factory=dict)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for command execution."""
    level = logging.DEBUG if verbose else logging.INFO

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def load_corrected_tsv(filepath: Path) -> List[CorrectedPosition]:
    """
    Load corrected positions from RECTIFY output TSV.

    Args:
        filepath: Path to corrected TSV file

    Returns:
        List of CorrectedPosition objects
    """
    positions = []
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Corrected TSV not found: {filepath}")

    with open(filepath, 'r') as f:
        # Read header
        header = f.readline().strip().split('\t')
        header_lower = [h.lower().replace('_', '') for h in header]

        # Map columns
        col_map = {}
        for i, col in enumerate(header_lower):
            col_map[col] = i

        # Read data
        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')

            try:
                # Required fields
                read_id = fields[col_map.get('readid', col_map.get('readname', 0))]
                chrom = fields[col_map.get('chrom', col_map.get('chromosome', 1))]
                strand = fields[col_map.get('strand', 2)]

                # Position fields - try various column names
                orig_idx = col_map.get('original3prime', col_map.get('originalposition', col_map.get('original', -1)))
                corr_idx = col_map.get('corrected3prime', col_map.get('correctedposition', col_map.get('corrected', -1)))

                if orig_idx >= 0 and orig_idx < len(fields):
                    original_3prime = int(fields[orig_idx])
                else:
                    continue

                if corr_idx >= 0 and corr_idx < len(fields):
                    corrected_3prime = int(fields[corr_idx])
                else:
                    corrected_3prime = original_3prime

                # Optional fields
                ambiguity_min = None
                ambiguity_max = None
                ambiguity_range = 0
                correction_applied = ""
                confidence = ""
                qc_flags = ""

                if 'ambiguitymin' in col_map and col_map['ambiguitymin'] < len(fields):
                    try:
                        ambiguity_min = int(fields[col_map['ambiguitymin']])
                    except ValueError:
                        pass

                if 'ambiguitymax' in col_map and col_map['ambiguitymax'] < len(fields):
                    try:
                        ambiguity_max = int(fields[col_map['ambiguitymax']])
                    except ValueError:
                        pass

                if 'ambiguityrange' in col_map and col_map['ambiguityrange'] < len(fields):
                    try:
                        ambiguity_range = int(fields[col_map['ambiguityrange']])
                    except ValueError:
                        pass

                if 'correctionapplied' in col_map and col_map['correctionapplied'] < len(fields):
                    correction_applied = fields[col_map['correctionapplied']]

                if 'confidence' in col_map and col_map['confidence'] < len(fields):
                    confidence = fields[col_map['confidence']]

                if 'qcflags' in col_map and col_map['qcflags'] < len(fields):
                    qc_flags = fields[col_map['qcflags']]

                positions.append(CorrectedPosition(
                    read_id=read_id,
                    chrom=chrom,
                    strand=strand,
                    original_3prime=original_3prime,
                    corrected_3prime=corrected_3prime,
                    ambiguity_min=ambiguity_min,
                    ambiguity_max=ambiguity_max,
                    ambiguity_range=ambiguity_range,
                    correction_applied=correction_applied,
                    confidence=confidence,
                    qc_flags=qc_flags,
                ))

            except (ValueError, IndexError, KeyError) as e:
                logger.debug(f"Line {line_num}: parse error ({e}), skipping")
                continue

    logger.info(f"Loaded {len(positions)} corrected positions from {filepath}")
    return positions


def load_ground_truth_tsv(filepath: Path) -> Dict[Tuple[str, int, str], Dict]:
    """
    Load known ground truth positions from TSV.

    Expected format:
        chrom    position    strand    [gene_name]    [source]

    Args:
        filepath: Path to ground truth TSV

    Returns:
        Dict keyed by (chrom, position, strand) -> {gene_name, source}
    """
    ground_truth = {}
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Ground truth file not found: {filepath}")

    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')
        header_lower = [h.lower() for h in header]

        chrom_idx = 0
        pos_idx = 1
        strand_idx = 2

        gene_idx = None
        if 'gene_name' in header_lower:
            gene_idx = header_lower.index('gene_name')
        elif 'gene' in header_lower:
            gene_idx = header_lower.index('gene')

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 3:
                continue

            try:
                chrom = fields[chrom_idx]
                position = int(fields[pos_idx])
                strand = fields[strand_idx]

                gene_name = None
                if gene_idx is not None and gene_idx < len(fields):
                    gene_name = fields[gene_idx]

                key = (chrom, position, strand)
                ground_truth[key] = {
                    'gene_name': gene_name,
                    'source': 'known_positions',
                }
            except (ValueError, IndexError):
                continue

    logger.info(f"Loaded {len(ground_truth)} ground truth positions from {filepath}")
    return ground_truth


def load_annotation_3prime_ends(filepath: Path) -> Dict[Tuple[str, int, str], Dict]:
    """
    Load annotated 3' ends from GTF/GFF file.

    Extracts gene/transcript 3' ends based on strand:
    - + strand: end coordinate
    - - strand: start coordinate

    Args:
        filepath: Path to GTF/GFF file

    Returns:
        Dict keyed by (chrom, position, strand) -> {gene_name, source}
    """
    ground_truth = {}
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Annotation file not found: {filepath}")

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]

            # Only process gene or transcript features
            if feature_type not in ['gene', 'transcript', 'mRNA']:
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])  # Already 0-based exclusive
            strand = fields[6]

            # Extract gene name from attributes
            attrs = fields[8]
            gene_name = None

            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split('=')[1]
                    break
                elif attr.startswith('Name='):
                    gene_name = attr.split('=')[1]
                    break
                elif attr.startswith('gene='):
                    gene_name = attr.split('=')[1]
                    break

            # Determine 3' end position
            if strand == '+':
                three_prime = end - 1  # 0-based inclusive
            else:
                three_prime = start

            key = (chrom, three_prime, strand)
            ground_truth[key] = {
                'gene_name': gene_name,
                'source': 'annotation',
            }

    logger.info(f"Loaded {len(ground_truth)} annotated 3' ends from {filepath}")
    return ground_truth


def find_nearest_netseq_peak(
    netseq_loader,
    chrom: str,
    position: int,
    strand: str,
    search_window: int = 10,
    min_signal: float = 0.5
) -> Optional[Dict]:
    """
    Find nearest NET-seq peak to a position.

    Args:
        netseq_loader: NetseqLoader instance
        chrom: Chromosome
        position: Position to search around
        strand: Strand
        search_window: Window size for search
        min_signal: Minimum signal threshold

    Returns:
        Dict with ground_truth_position, signal, distance or None
    """
    start = max(0, position - search_window)
    end = position + search_window + 1

    try:
        signal = netseq_loader.get_signal(chrom, start, end, strand)
    except Exception:
        return None

    if len(signal) == 0:
        return None

    max_signal = np.max(signal)
    if max_signal < min_signal:
        return None

    # Find position of maximum signal
    max_idx = np.argmax(signal)
    peak_position = start + max_idx
    peak_signal = signal[max_idx]

    return {
        'ground_truth_position': peak_position,
        'signal': float(peak_signal),
        'distance': peak_position - position,
        'source': 'netseq',
    }


def find_nearest_annotation(
    annotations: Dict[Tuple[str, int, str], Dict],
    chrom: str,
    position: int,
    strand: str,
    search_window: int = 50
) -> Optional[Dict]:
    """
    Find nearest annotated 3' end to a position.

    Args:
        annotations: Dict of annotated positions
        chrom: Chromosome
        position: Position to search around
        strand: Strand
        search_window: Maximum distance to search

    Returns:
        Dict with ground_truth_position, gene_name, distance or None
    """
    best_match = None
    best_distance = search_window + 1

    for (ann_chrom, ann_pos, ann_strand), info in annotations.items():
        if ann_chrom != chrom or ann_strand != strand:
            continue

        distance = abs(ann_pos - position)
        if distance < best_distance:
            best_distance = distance
            best_match = {
                'ground_truth_position': ann_pos,
                'gene_name': info.get('gene_name'),
                'distance': ann_pos - position,
                'source': 'annotation',
            }

    return best_match


def classify_correction_outcome(
    original: int,
    corrected: int,
    ground_truth: int
) -> str:
    """
    Classify outcome of correction.

    Returns:
        'IMPROVED': Correction moved position closer to truth
        'UNCHANGED': No change (or same distance)
        'WORSENED': Correction moved position farther from truth
    """
    orig_dist = abs(original - ground_truth)
    corr_dist = abs(corrected - ground_truth)

    if corr_dist < orig_dist:
        return 'IMPROVED'
    elif corr_dist > orig_dist:
        return 'WORSENED'
    else:
        return 'UNCHANGED'


def validate_position(
    position: CorrectedPosition,
    netseq_loader=None,
    annotations: Optional[Dict] = None,
    known_positions: Optional[Dict] = None,
    tolerance: int = 1,
    search_window: int = 10,
    min_signal: float = 0.5
) -> ValidationResult:
    """
    Validate a single corrected position against ground truth.

    Tries ground truth sources in order of specificity:
    1. Known positions (exact match)
    2. NET-seq peaks
    3. Annotation 3' ends

    Args:
        position: CorrectedPosition to validate
        netseq_loader: Optional NetseqLoader instance
        annotations: Optional annotation positions dict
        known_positions: Optional known positions dict
        tolerance: Distance tolerance for "correct" classification
        search_window: Window for searching nearby ground truth
        min_signal: Minimum NET-seq signal

    Returns:
        ValidationResult
    """
    result = ValidationResult(
        read_id=position.read_id,
        chrom=position.chrom,
        strand=position.strand,
        original_3prime=position.original_3prime,
        corrected_3prime=position.corrected_3prime,
        correction_applied=position.correction_applied,
        confidence=position.confidence,
    )

    ground_truth = None
    source = "none"

    # Try known positions first (highest priority)
    if known_positions:
        key = (position.chrom, position.corrected_3prime, position.strand)
        if key in known_positions:
            ground_truth = position.corrected_3prime
            source = "known_exact"
        else:
            # Search nearby
            for pos_offset in range(-search_window, search_window + 1):
                test_key = (position.chrom, position.corrected_3prime + pos_offset, position.strand)
                if test_key in known_positions:
                    ground_truth = position.corrected_3prime + pos_offset
                    source = "known_positions"
                    break

    # Try NET-seq if available and no known position found
    if ground_truth is None and netseq_loader is not None:
        netseq_result = find_nearest_netseq_peak(
            netseq_loader,
            position.chrom,
            position.corrected_3prime,
            position.strand,
            search_window=search_window,
            min_signal=min_signal
        )
        if netseq_result:
            ground_truth = netseq_result['ground_truth_position']
            source = "netseq"

    # Try annotations if still no ground truth
    if ground_truth is None and annotations:
        ann_result = find_nearest_annotation(
            annotations,
            position.chrom,
            position.corrected_3prime,
            position.strand,
            search_window=search_window * 5  # Larger window for annotations
        )
        if ann_result:
            ground_truth = ann_result['ground_truth_position']
            source = "annotation"

    # Populate result
    result.ground_truth_source = source

    if ground_truth is not None:
        result.ground_truth_position = ground_truth
        result.distance_from_truth = position.corrected_3prime - ground_truth
        result.original_distance = position.original_3prime - ground_truth

        abs_dist = abs(result.distance_from_truth)
        abs_orig = abs(result.original_distance)

        result.is_correct_exact = abs_dist == 0
        result.is_correct_1bp = abs_dist <= 1
        result.is_correct_2bp = abs_dist <= 2
        result.is_correct_5bp = abs_dist <= 5

        result.improvement_bp = abs_orig - abs_dist
        result.validation_status = classify_correction_outcome(
            position.original_3prime,
            position.corrected_3prime,
            ground_truth
        )
    else:
        result.validation_status = "NO_TRUTH"

    return result


def validate_batch(
    positions: List[CorrectedPosition],
    netseq_loader=None,
    annotations: Optional[Dict] = None,
    known_positions: Optional[Dict] = None,
    tolerance: int = 1,
    min_signal: float = 0.5
) -> List[ValidationResult]:
    """
    Validate batch of corrected positions.

    Args:
        positions: List of CorrectedPosition objects
        netseq_loader: Optional NetseqLoader
        annotations: Optional annotation positions
        known_positions: Optional known positions
        tolerance: Distance tolerance
        min_signal: Minimum NET-seq signal

    Returns:
        List of ValidationResult objects
    """
    results = []

    for pos in tqdm(positions, desc="Validating positions"):
        result = validate_position(
            pos,
            netseq_loader=netseq_loader,
            annotations=annotations,
            known_positions=known_positions,
            tolerance=tolerance,
            min_signal=min_signal
        )
        results.append(result)

    return results


def calculate_metrics(results: List[ValidationResult]) -> ValidationMetrics:
    """
    Calculate aggregate validation metrics.

    Args:
        results: List of ValidationResult objects

    Returns:
        ValidationMetrics object
    """
    metrics = ValidationMetrics()
    metrics.total_positions = len(results)

    for r in results:
        if r.ground_truth_position is not None:
            metrics.with_ground_truth += 1

            if r.is_correct_exact:
                metrics.exact_match += 1
            if r.is_correct_1bp:
                metrics.within_1bp += 1
            if r.is_correct_2bp:
                metrics.within_2bp += 1
            if r.is_correct_5bp:
                metrics.within_5bp += 1

            # Original accuracy
            abs_orig = abs(r.original_distance) if r.original_distance else 0
            if abs_orig == 0:
                metrics.original_exact_match += 1
            if abs_orig <= 1:
                metrics.original_within_1bp += 1
            if abs_orig <= 2:
                metrics.original_within_2bp += 1

            # Improvement
            if r.validation_status == 'IMPROVED':
                metrics.improved += 1
            elif r.validation_status == 'WORSENED':
                metrics.worsened += 1
            else:
                metrics.unchanged += 1

            # By correction type
            corr_type = r.correction_applied or 'none'
            if corr_type not in metrics.by_correction:
                metrics.by_correction[corr_type] = {
                    'total': 0, 'exact': 0, 'within_1bp': 0,
                    'improved': 0, 'worsened': 0
                }
            metrics.by_correction[corr_type]['total'] += 1
            if r.is_correct_exact:
                metrics.by_correction[corr_type]['exact'] += 1
            if r.is_correct_1bp:
                metrics.by_correction[corr_type]['within_1bp'] += 1
            if r.validation_status == 'IMPROVED':
                metrics.by_correction[corr_type]['improved'] += 1
            elif r.validation_status == 'WORSENED':
                metrics.by_correction[corr_type]['worsened'] += 1

            # By confidence
            conf = r.confidence or 'unknown'
            if conf not in metrics.by_confidence:
                metrics.by_confidence[conf] = {
                    'total': 0, 'exact': 0, 'within_1bp': 0
                }
            metrics.by_confidence[conf]['total'] += 1
            if r.is_correct_exact:
                metrics.by_confidence[conf]['exact'] += 1
            if r.is_correct_1bp:
                metrics.by_confidence[conf]['within_1bp'] += 1

    return metrics


def write_validation_tsv(results: List[ValidationResult], output_path: Path) -> None:
    """Write validation results to TSV file."""
    output_path = Path(output_path)

    with open(output_path, 'w') as f:
        # Header
        header = [
            'read_id', 'chrom', 'strand',
            'original_3prime', 'corrected_3prime', 'ground_truth_position',
            'ground_truth_source', 'distance_from_truth', 'original_distance',
            'is_correct_exact', 'is_correct_1bp', 'is_correct_2bp',
            'improvement_bp', 'validation_status',
            'correction_applied', 'confidence'
        ]
        f.write('\t'.join(header) + '\n')

        for r in results:
            row = [
                r.read_id, r.chrom, r.strand,
                str(r.original_3prime), str(r.corrected_3prime),
                str(r.ground_truth_position) if r.ground_truth_position is not None else '',
                r.ground_truth_source,
                str(r.distance_from_truth) if r.distance_from_truth is not None else '',
                str(r.original_distance) if r.original_distance is not None else '',
                str(r.is_correct_exact), str(r.is_correct_1bp), str(r.is_correct_2bp),
                str(r.improvement_bp), r.validation_status,
                r.correction_applied, r.confidence
            ]
            f.write('\t'.join(row) + '\n')

    logger.info(f"Wrote validation results to {output_path}")


def generate_summary_report(metrics: ValidationMetrics) -> str:
    """Generate human-readable summary report."""
    lines = []
    lines.append("=" * 70)
    lines.append("RECTIFY Validation Report")
    lines.append("=" * 70)
    lines.append("")

    # Overall
    lines.append("Overall:")
    lines.append(f"  Total positions:        {metrics.total_positions:,}")
    lines.append(f"  With ground truth:      {metrics.with_ground_truth:,} "
                f"({100*metrics.with_ground_truth/max(1, metrics.total_positions):.1f}%)")
    lines.append("")

    if metrics.with_ground_truth > 0:
        n = metrics.with_ground_truth

        # Corrected accuracy
        lines.append("Corrected Position Accuracy:")
        lines.append(f"  Exact match:            {metrics.exact_match:,} ({100*metrics.exact_match/n:.1f}%)")
        lines.append(f"  Within 1 bp:            {metrics.within_1bp:,} ({100*metrics.within_1bp/n:.1f}%)")
        lines.append(f"  Within 2 bp:            {metrics.within_2bp:,} ({100*metrics.within_2bp/n:.1f}%)")
        lines.append(f"  Within 5 bp:            {metrics.within_5bp:,} ({100*metrics.within_5bp/n:.1f}%)")
        lines.append("")

        # Original accuracy (comparison)
        lines.append("Original Position Accuracy (Before Correction):")
        lines.append(f"  Exact match:            {metrics.original_exact_match:,} ({100*metrics.original_exact_match/n:.1f}%)")
        lines.append(f"  Within 1 bp:            {metrics.original_within_1bp:,} ({100*metrics.original_within_1bp/n:.1f}%)")
        lines.append(f"  Within 2 bp:            {metrics.original_within_2bp:,} ({100*metrics.original_within_2bp/n:.1f}%)")
        lines.append("")

        # Improvement
        improvement = metrics.within_1bp - metrics.original_within_1bp
        improvement_pct = 100 * improvement / n if n > 0 else 0

        lines.append("Improvement Over Original:")
        lines.append(f"  Improved:               {metrics.improved:,} ({100*metrics.improved/n:.1f}%)")
        lines.append(f"  Unchanged:              {metrics.unchanged:,} ({100*metrics.unchanged/n:.1f}%)")
        lines.append(f"  Worsened:               {metrics.worsened:,} ({100*metrics.worsened/n:.1f}%)")
        lines.append(f"  Net ±1bp improvement:   {improvement:+,} ({improvement_pct:+.1f} percentage points)")
        lines.append("")

        # By correction type
        if metrics.by_correction:
            lines.append("By Correction Type:")
            for corr_type, stats in sorted(metrics.by_correction.items()):
                total = stats['total']
                exact_pct = 100 * stats['exact'] / total if total > 0 else 0
                bp1_pct = 100 * stats['within_1bp'] / total if total > 0 else 0
                lines.append(f"  {corr_type:25s} n={total:6,}  exact={exact_pct:5.1f}%  ±1bp={bp1_pct:5.1f}%")
            lines.append("")

        # By confidence
        if metrics.by_confidence:
            lines.append("By Confidence Level:")
            for conf, stats in sorted(metrics.by_confidence.items()):
                total = stats['total']
                exact_pct = 100 * stats['exact'] / total if total > 0 else 0
                bp1_pct = 100 * stats['within_1bp'] / total if total > 0 else 0
                lines.append(f"  {conf:15s}  n={total:6,}  exact={exact_pct:5.1f}%  ±1bp={bp1_pct:5.1f}%")
            lines.append("")

    lines.append("=" * 70)
    return "\n".join(lines)


def validate_inputs(args) -> Dict:
    """Validate command arguments and return config."""
    errors = []

    # Check corrected TSV
    if not args.corrected.exists():
        errors.append(f"Corrected TSV not found: {args.corrected}")

    # Check at least one ground truth source
    has_ground_truth = False

    if hasattr(args, 'netseq_dir') and args.netseq_dir:
        if args.netseq_dir.exists():
            has_ground_truth = True
        else:
            errors.append(f"NET-seq directory not found: {args.netseq_dir}")

    if hasattr(args, 'annotation') and args.annotation:
        if args.annotation.exists():
            has_ground_truth = True
        else:
            errors.append(f"Annotation file not found: {args.annotation}")

    if hasattr(args, 'ground_truth') and args.ground_truth:
        if args.ground_truth.exists():
            has_ground_truth = True
        else:
            errors.append(f"Ground truth file not found: {args.ground_truth}")

    if not has_ground_truth:
        errors.append("At least one ground truth source required: "
                     "--netseq-dir, --annotation, or --ground-truth")

    # Check output directory
    output_dir = args.output.parent
    if output_dir and not output_dir.exists():
        try:
            output_dir.mkdir(parents=True)
        except Exception as e:
            errors.append(f"Cannot create output directory {output_dir}: {e}")

    if errors:
        for error in errors:
            logger.error(error)
        sys.exit(1)

    return {
        'tolerance': getattr(args, 'tolerance', 1),
        'min_signal': getattr(args, 'min_signal', 0.5),
        'search_window': getattr(args, 'search_window', 10),
    }


def run(args) -> None:
    """
    Execute the 'validate' command.

    Args:
        args: Parsed command-line arguments from argparse
    """
    # Setup
    set_thread_limits(getattr(args, 'threads', None))
    setup_logging(getattr(args, 'verbose', False))

    logger.info("=" * 70)
    logger.info("RECTIFY - Validation Against Ground Truth")
    logger.info("=" * 70)
    logger.info("")

    # Validate inputs
    config = validate_inputs(args)

    logger.info("Inputs:")
    logger.info(f"  Corrected TSV:    {args.corrected}")
    if hasattr(args, 'netseq_dir') and args.netseq_dir:
        logger.info(f"  NET-seq dir:      {args.netseq_dir}")
    if hasattr(args, 'annotation') and args.annotation:
        logger.info(f"  Annotation:       {args.annotation}")
    if hasattr(args, 'ground_truth') and args.ground_truth:
        logger.info(f"  Ground truth:     {args.ground_truth}")
    logger.info(f"  Output:           {args.output}")
    logger.info(f"  Tolerance:        {config['tolerance']} bp")
    logger.info("")

    # Load corrected positions
    logger.info("Loading corrected positions...")
    positions = load_corrected_tsv(args.corrected)

    if len(positions) == 0:
        logger.error("No positions loaded from corrected TSV")
        sys.exit(1)

    # Load ground truth sources
    netseq_loader = None
    annotations = None
    known_positions = None

    # NET-seq (optional)
    if hasattr(args, 'netseq_dir') and args.netseq_dir:
        try:
            from .netseq_refiner import NetseqLoader
            logger.info("Loading NET-seq data...")
            netseq_loader = NetseqLoader()
            netseq_loader.load_directory(str(args.netseq_dir))
            logger.info(f"  Loaded {len(netseq_loader.bigwigs)} BigWig files")
        except ImportError:
            logger.warning("pyBigWig not available, skipping NET-seq validation")

    # Annotations (optional)
    if hasattr(args, 'annotation') and args.annotation:
        logger.info("Loading gene annotations...")
        annotations = load_annotation_3prime_ends(args.annotation)

    # Known positions (optional)
    if hasattr(args, 'ground_truth') and args.ground_truth:
        logger.info("Loading known positions...")
        known_positions = load_ground_truth_tsv(args.ground_truth)

    # Validate positions
    logger.info("Validating positions...")
    results = validate_batch(
        positions,
        netseq_loader=netseq_loader,
        annotations=annotations,
        known_positions=known_positions,
        tolerance=config['tolerance'],
        min_signal=config['min_signal']
    )

    # Calculate metrics
    logger.info("Calculating metrics...")
    metrics = calculate_metrics(results)

    # Write results
    write_validation_tsv(results, args.output)

    # Generate and print summary
    summary = generate_summary_report(metrics)
    logger.info("")
    logger.info(summary)

    # Cleanup
    if netseq_loader:
        netseq_loader.close()

    logger.info("")
    logger.info("=" * 70)
    logger.info("Validation complete!")
    logger.info("=" * 70)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('corrected', type=Path)
    parser.add_argument('--netseq-dir', type=Path)
    parser.add_argument('--annotation', type=Path)
    parser.add_argument('--ground-truth', type=Path)
    parser.add_argument('--netseq-samples', nargs='+')
    parser.add_argument('-o', '--output', type=Path, required=True)
    parser.add_argument('--tolerance', type=int, default=1)
    parser.add_argument('--min-signal', type=float, default=0.5)
    parser.add_argument('--search-window', type=int, default=10)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()
    run(args)
