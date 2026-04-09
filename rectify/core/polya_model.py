#!/usr/bin/env python3
"""
Poly(A) tail model for RECTIFY.

This module defines the PolyAModel class for representing learned poly(A)
tail characteristics. Models can be trained from control sites (0 downstream
A's) where any A-rich soft-clip must be genuine poly(A) tail.

Author: Kevin R. Roy
Date: 2026-03-16
"""

from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
import json
import logging
import numpy as np

from ..config import (
    POLYA_RICHNESS_THRESHOLD,
    POLYA_WINDOW_SIZE,
    MIN_POLYA_LENGTH,
)

logger = logging.getLogger(__name__)


@dataclass
class PolyAModelParameters:
    """Core parameters of a poly(A) tail model."""

    # A-richness parameters
    a_richness_threshold: float = POLYA_RICHNESS_THRESHOLD
    a_richness_mean: float = 0.91
    a_richness_std: float = 0.08

    # Length parameters
    min_tail_length: int = MIN_POLYA_LENGTH
    mean_tail_length: float = 28.5
    median_tail_length: float = 25.0
    std_tail_length: float = 15.2
    max_tail_length: int = 200

    # Window size for A-richness calculation
    window_size: int = POLYA_WINDOW_SIZE


@dataclass
class PolyAModelStats:
    """Training statistics for a poly(A) model."""

    n_control_sites: int = 0
    n_reads_used: int = 0
    n_soft_clips_analyzed: int = 0
    n_sites_filtered: int = 0


@dataclass
class PolyAModel:
    """
    Poly(A) tail model learned from control data.

    Control sites are CPA positions with 0 downstream A's, where any
    A-rich soft-clipped sequence must genuinely be a poly(A) tail
    (not genomic sequence misaligned).

    Attributes:
        version: Model format version
        created: Timestamp when model was created
        technology: Sequencing technology (e.g., 'nanopore', 'quantseq')
        genome: Reference genome used for training
        parameters: Core model parameters
        training_stats: Statistics from training
        position_profile: Position-dependent A-frequency (0 = start of tail)
        non_a_frequencies: Frequency of non-A bases in poly(A) tails
        adapter_patterns: Detected adapter characteristics
    """

    version: str = "1.0"
    created: str = field(default_factory=lambda: datetime.now().isoformat())
    technology: str = "unknown"
    genome: str = "unknown"
    parameters: PolyAModelParameters = field(default_factory=PolyAModelParameters)
    training_stats: PolyAModelStats = field(default_factory=PolyAModelStats)
    position_profile: Dict[int, float] = field(default_factory=dict)
    non_a_frequencies: Dict[str, float] = field(default_factory=lambda: {
        'G': 0.04, 'C': 0.02, 'T': 0.01, 'N': 0.0
    })
    adapter_patterns: Dict[str, Any] = field(default_factory=lambda: {
        'poly_t_min_detected': 6,
        'has_tc_motifs': True,
    })

    def to_dict(self) -> Dict[str, Any]:
        """Convert model to dictionary for JSON serialization."""
        return {
            'version': self.version,
            'created': self.created,
            'technology': self.technology,
            'genome': self.genome,
            'parameters': asdict(self.parameters),
            'training_stats': asdict(self.training_stats),
            'position_profile': {str(k): v for k, v in self.position_profile.items()},
            'non_a_frequencies': self.non_a_frequencies,
            'adapter_patterns': self.adapter_patterns,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PolyAModel':
        """Create model from dictionary (JSON deserialization)."""
        model = cls()
        model.version = data.get('version', '1.0')
        model.created = data.get('created', datetime.now().isoformat())
        model.technology = data.get('technology', 'unknown')
        model.genome = data.get('genome', 'unknown')

        # Load parameters
        if 'parameters' in data:
            model.parameters = PolyAModelParameters(**data['parameters'])

        # Load training stats
        if 'training_stats' in data:
            model.training_stats = PolyAModelStats(**data['training_stats'])

        # Load position profile (convert string keys to int)
        if 'position_profile' in data:
            position_profile = {}
            for k, v in data['position_profile'].items():
                try:
                    position_profile[int(k)] = v
                except (ValueError, TypeError):
                    logger.debug(f"Skipping non-integer position_profile key: {k!r}")
            model.position_profile = position_profile

        # Load other fields
        if 'non_a_frequencies' in data:
            model.non_a_frequencies = data['non_a_frequencies']

        if 'adapter_patterns' in data:
            model.adapter_patterns = data['adapter_patterns']

        return model

    def to_json(self, path: Path) -> None:
        """Save model to JSON file."""
        path = Path(path)
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
        logger.info(f"Saved poly(A) model to {path}")

    @classmethod
    def from_json(cls, path: Path) -> 'PolyAModel':
        """Load model from JSON file."""
        path = Path(path)
        with open(path, 'r') as f:
            data = json.load(f)
        logger.info(f"Loaded poly(A) model from {path}")
        return cls.from_dict(data)

    def score_sequence(self, sequence: str) -> Dict[str, Any]:
        """
        Score a sequence using learned parameters.

        Args:
            sequence: DNA sequence (typically soft-clip)

        Returns:
            Dict with:
                - a_richness: Maximum A-richness across windows
                - is_polya: True if classified as poly(A)
                - confidence: Confidence score (0-1)
                - length_score: Score based on length distribution
                - profile_score: Score based on position profile
        """
        if not sequence:
            return {
                'a_richness': 0.0,
                'is_polya': False,
                'confidence': 0.0,
                'length_score': 0.0,
                'profile_score': 0.0,
            }

        sequence = sequence.upper()
        params = self.parameters

        # Calculate A-richness using sliding window
        a_richness = self._calculate_a_richness(sequence, params.window_size)

        # Basic classification
        is_polya = a_richness >= params.a_richness_threshold

        # Length score (how well length matches training distribution)
        length_score = self._score_length(len(sequence))

        # Position profile score
        profile_score = self._score_position_profile(sequence)

        # Combined confidence
        if is_polya:
            # Weight factors
            richness_weight = 0.5
            length_weight = 0.3
            profile_weight = 0.2

            confidence = (
                richness_weight * min(1.0, a_richness / params.a_richness_threshold) +
                length_weight * length_score +
                profile_weight * profile_score
            )
        else:
            confidence = 0.0

        return {
            'a_richness': a_richness,
            'is_polya': is_polya,
            'confidence': confidence,
            'length_score': length_score,
            'profile_score': profile_score,
        }

    def _calculate_a_richness(self, sequence: str, window: int) -> float:
        """Calculate maximum A-richness using sliding window."""
        if len(sequence) == 0:
            return 0.0

        if len(sequence) < window:
            return sequence.count('A') / len(sequence)

        max_richness = 0.0
        a_count = sequence[:window].count('A')
        max_richness = a_count / window

        # Sliding window with incremental update
        for i in range(1, len(sequence) - window + 1):
            # Remove leftmost base, add rightmost base
            if sequence[i - 1] == 'A':
                a_count -= 1
            if sequence[i + window - 1] == 'A':
                a_count += 1

            richness = a_count / window
            max_richness = max(max_richness, richness)

        return max_richness

    def _score_length(self, length: int) -> float:
        """
        Score sequence length against learned distribution.

        Returns score between 0 and 1, where 1 means perfect match
        to the expected length distribution.
        """
        params = self.parameters

        if length < params.min_tail_length:
            # Too short - penalize proportionally
            return length / params.min_tail_length

        if length > params.max_tail_length:
            # Too long - slight penalty
            return 0.8

        # Score based on distance from mean
        if params.std_tail_length > 0:
            z_score = abs(length - params.mean_tail_length) / params.std_tail_length
            # Convert to 0-1 score (higher is better)
            score = np.exp(-0.5 * z_score ** 2)
            return float(score)

        return 1.0 if length >= params.min_tail_length else 0.0

    def _score_position_profile(self, sequence: str) -> float:
        """
        Score sequence against learned position profile.

        The position profile represents expected A-frequency at each
        position from the start of the poly(A) tail.
        """
        if not self.position_profile:
            # No profile learned - return neutral score
            return 0.5

        if len(sequence) == 0:
            return 0.0

        # Calculate observed A-frequency at each position
        matches = 0
        total = 0

        for i, base in enumerate(sequence):
            expected_a_freq = self.position_profile.get(i, 0.9)  # Default high A-freq

            if base == 'A':
                # Match if we expected high A-frequency
                if expected_a_freq >= 0.5:
                    matches += 1
            else:
                # Match if we expected low A-frequency (near adapter)
                if expected_a_freq < 0.5:
                    matches += 1

            total += 1

        return matches / total if total > 0 else 0.0

    def is_polya(self, sequence: str) -> bool:
        """Classify sequence as poly(A) or not."""
        return self.score_sequence(sequence)['is_polya']

    def get_threshold(self) -> float:
        """Get the A-richness threshold."""
        return self.parameters.a_richness_threshold

    def summary(self) -> str:
        """Return human-readable model summary."""
        lines = [
            "=" * 60,
            "Poly(A) Tail Model Summary",
            "=" * 60,
            "",
            f"Version:              {self.version}",
            f"Created:              {self.created}",
            f"Technology:           {self.technology}",
            f"Genome:               {self.genome}",
            "",
            "Training Statistics:",
            f"  Control sites:      {self.training_stats.n_control_sites:,}",
            f"  Reads used:         {self.training_stats.n_reads_used:,}",
            f"  Soft-clips analyzed:{self.training_stats.n_soft_clips_analyzed:,}",
            "",
            "Learned Parameters:",
            f"  A-richness threshold: {self.parameters.a_richness_threshold:.2f}",
            f"  A-richness mean:      {self.parameters.a_richness_mean:.2f}",
            f"  Min tail length:      {self.parameters.min_tail_length} bp",
            f"  Mean tail length:     {self.parameters.mean_tail_length:.1f} bp",
            f"  Median tail length:   {self.parameters.median_tail_length:.1f} bp",
            "",
            "Non-A Base Frequencies:",
        ]

        for base, freq in sorted(self.non_a_frequencies.items()):
            lines.append(f"  {base}: {freq:.3f}")

        lines.extend(["", "=" * 60])
        return "\n".join(lines)


def load_model(path: Optional[Path]) -> Optional[PolyAModel]:
    """
    Load poly(A) model from file, or return None for built-in defaults.

    Args:
        path: Path to model JSON file, or None for defaults

    Returns:
        PolyAModel instance, or None to use built-in defaults
    """
    if path is None:
        return None

    path = Path(path)
    if not path.exists():
        logger.warning(f"Model file not found: {path}, using defaults")
        return None

    try:
        return PolyAModel.from_json(path)
    except Exception as e:
        logger.error(f"Failed to load model from {path}: {e}")
        return None


def get_default_model() -> PolyAModel:
    """Get default poly(A) model with standard parameters."""
    return PolyAModel(
        technology="default",
        genome="generic",
    )
