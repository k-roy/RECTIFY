"""RECTIFY native-aligner simulation ground-truth benchmark (Deliverable A).

This package is the **fitness function** (the GATE) for any future native-aligner
member code: every primitive is a hypothesis-pending-ablation against THIS truth
set, NEVER the internal score (which was provably artifact-prone — the GMAP
0.09->1.07 re-weighting flip with no aligner change is the entire reason this gate
exists).

Modules
-------
truth_schema
    The per-read ground-truth table + sidecar format (component 2 — built FIRST,
    it constrains the simulator wrapper and the scorer). Carries exact-indel
    truth (per-position true edit + ambiguity-equivalence window), junction truth
    (donor/acceptor + NIC/NNC class), CPA truth, the shared genomic-region-disjoint
    train/test split tag, and the C1 HP/STR cell metadata.
scorer
    The ambiguity-aware per-junction + per-indel accuracy scorer (component 3).
    Framing metric = EXACT INDEL-POSITION CONCORDANCE WITH TRUTH, ambiguity-aware
    (reusing chimeric_consensus.normalize_junction / _canonical_within_window /
    junction_ambiguity_window), NEVER edit distance.

The simulator wrapper (component 1) lives under ``scripts/benchmark/sim/`` because
it shells out to pbsim3 / badread and is not importable library code.
"""

from .truth_schema import (  # noqa: F401
    JunctionTruth,
    IndelTruth,
    VariantTruth,
    ReadTruth,
    JunctionClass,
    IndelKind,
    SplitTag,
    write_truth_table,
    read_truth_table,
)

__all__ = [
    "JunctionTruth",
    "IndelTruth",
    "VariantTruth",
    "ReadTruth",
    "JunctionClass",
    "IndelKind",
    "SplitTag",
    "write_truth_table",
    "read_truth_table",
]
