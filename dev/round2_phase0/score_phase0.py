"""Round-2 cDNA Phase-0 harness — scoring, win-guard, and GO/NO-GO verdict.

Implements the Designer-2 win-guard (§1.4) and the master-spec BLOCKER-1 fix (§7):
a lifted cDNA candidate may replace the genome winner for a read IFF, after lift-over:
  (1) it passes the anchor gate on its TRANSFERRED N-ops, gated on the RAW
      ``min_junction_anchor >= min_junction_anchor_bp`` -- NOT on ``_effective_chimera_ok``
      (the Cat3 / five_prime_rescued exemption would otherwise switch the gate off
      vacuously: BLOCKER-1); AND
  (2) strict improvement: hp_ed_cdna <= hp_ed_genome - epsilon (epsilon = 1.0); AND
  (3) no-shrink: aligned_bases_cdna >= aligned_bases_genome.

All three metrics are computed with the PRODUCTION functions imported from
``rectify.core.consensus.corrected_consensus`` -- so a Phase-0 GO/NO-GO verdict is
produced by the same code production winner-selection uses (advisor requirement: if the
verdict uses the same metric production will use, the verdict is trustworthy).

The verdict also enforces the two GO criteria the spec flags as easy-to-drop:
  - **uLTRA-specific win** (§9): a cDNA win must beat uLTRA's OWN per-read record, not
    just the consensus winner -- else the feature is redundant with annotation-guided uLTRA.
  - **per-win attribution** (§9): every cDNA win must carry an attributable cause
    (recovered-clip / placed-micro-exon / spanned-multi-junction). >5% no-cause leak is a
    NO-GO falsifier (reverse-gaming).
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

# Import the PRODUCTION scorers. Adjust sys.path to the repo root if run standalone.
_REPO = Path(__file__).resolve().parents[2]   # .../rectify
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from rectify.core.consensus.corrected_consensus import (  # noqa: E402
    _cigar_hp_edit_distance,
    _cigar_aligned_bases,
    _cigar_min_junction_anchor,
    _NO_JUNCTION_ANCHOR,
)

EPSILON = 1.0            # win-guard (2): must beat by >= one mismatch-equivalent
WIN_CAUSES = ("recovered_clip", "placed_micro_exon", "spanned_multi_junction")


@dataclass
class Candidate:
    """A scored alignment candidate for one read (genome aligner OR lifted cDNA)."""
    read_id: str
    aligner: str                 # e.g. 'minimap2', 'uLTRA', 'cdna_round2'
    hp_ed: float
    aligned_bases: int
    min_junction_anchor: int     # raw per-read min junction anchor (sentinel if no N-ops)
    five_prime_rescued: int = 0  # Cat3 flag; BLOCKER-1: must NOT switch off the gate
    # diagnostics for attribution (optional)
    genome_softclip_bp: int = 0  # the genome winner's effective soft-clip this read had
    n_junctions: int = 0

    @classmethod
    def from_read(cls, read, genome, penalty_table, *, aligner, read_id,
                  five_prime_rescued=0, genome_softclip_bp=0):
        return cls(
            read_id=read_id,
            aligner=aligner,
            hp_ed=_cigar_hp_edit_distance(read, genome, penalty_table),
            aligned_bases=_cigar_aligned_bases(read),
            min_junction_anchor=_cigar_min_junction_anchor(read, genome),
            five_prime_rescued=five_prime_rescued,
            genome_softclip_bp=genome_softclip_bp,
            n_junctions=sum(1 for op, _ in (read.cigartuples or []) if op == 3),
        )


def anchor_gate_passes(cand: Candidate, min_junction_anchor_bp: int) -> bool:
    """BLOCKER-1 FIX: gate on the RAW anchor, never on _effective_chimera_ok.

    A candidate with no N-ops returns the sentinel and is never gated (passes).
    A candidate carrying five_prime_rescued=1 does NOT get a free pass here -- that is
    precisely the bypass the master spec calls BLOCKER-1.
    """
    if min_junction_anchor_bp <= 0:
        # Safety precondition (§3.2 / §7): Round 2 must not run with the gate off.
        raise ValueError(
            "anchor gate is OFF (min_junction_anchor_bp<=0); Round-2 cDNA must refuse to run"
        )
    if cand.min_junction_anchor >= _NO_JUNCTION_ANCHOR:
        return True   # no junctions -> not gated
    return cand.min_junction_anchor >= min_junction_anchor_bp


@dataclass
class WinDecision:
    read_id: str
    cdna_wins: bool
    reason: str                       # why it won, or why it was kept on genome
    cause: Optional[str] = None       # attribution (one of WIN_CAUSES) when cdna_wins
    beats_ultra: Optional[bool] = None
    hp_ed_genome: float = 0.0
    hp_ed_cdna: float = 0.0


def decide(
    cdna: Candidate,
    genome_winner: Candidate,
    *,
    ultra: Optional[Candidate] = None,
    min_junction_anchor_bp: int = 10,
    epsilon: float = EPSILON,
) -> WinDecision:
    """Apply the win-guard to one read. ``genome_winner`` is the best genome candidate;
    ``ultra`` is uLTRA's own record for the same read (for the redundancy metric)."""
    # (1) anchor gate on the cDNA's transferred N-ops (raw anchor)
    if not anchor_gate_passes(cdna, min_junction_anchor_bp):
        return WinDecision(cdna.read_id, False, "cdna_failed_anchor_gate",
                           hp_ed_genome=genome_winner.hp_ed, hp_ed_cdna=cdna.hp_ed)
    # (2) strict improvement
    if not (cdna.hp_ed <= genome_winner.hp_ed - epsilon):
        return WinDecision(cdna.read_id, False, "not_strictly_better",
                           hp_ed_genome=genome_winner.hp_ed, hp_ed_cdna=cdna.hp_ed)
    # (3) no-shrink guard
    if not (cdna.aligned_bases >= genome_winner.aligned_bases):
        return WinDecision(cdna.read_id, False, "no_shrink_veto",
                           hp_ed_genome=genome_winner.hp_ed, hp_ed_cdna=cdna.hp_ed)

    cause = attribute_win(cdna, genome_winner)
    beats_ultra = None
    if ultra is not None:
        # did the cDNA beat uLTRA's OWN record specifically (§9 must-fix)?
        beats_ultra = cdna.hp_ed <= ultra.hp_ed - epsilon
    return WinDecision(cdna.read_id, True, "cdna_wins", cause=cause,
                       beats_ultra=beats_ultra,
                       hp_ed_genome=genome_winner.hp_ed, hp_ed_cdna=cdna.hp_ed)


def attribute_win(cdna: Candidate, genome_winner: Candidate) -> Optional[str]:
    """Assign an attributable cause to a cDNA win (§9 / §4.1). Returns None if no cause
    can be attributed -- those count toward the trivial-win-leak NO-GO falsifier.

    Most-specific structural cause first: spanning an extra junction is unambiguous
    evidence the genome aligner missed a micro-exon/junction, so it outranks the generic
    'recovered more aligned bases over a soft-clip' signal."""
    # placed-micro-exon: the cDNA spans MORE junctions than the genome winner
    # (genome aligner skipped/fused a micro-exon it couldn't seed)
    if cdna.n_junctions > genome_winner.n_junctions:
        return "placed_micro_exon"
    # spanned-multi-junction: same junction count (>=2) but materially more aligned bases
    if cdna.n_junctions >= 2 and cdna.aligned_bases > genome_winner.aligned_bases:
        return "spanned_multi_junction"
    # recovered-clip: the genome winner carried a real soft-clip the cDNA placed as matches
    if genome_winner.genome_softclip_bp >= 20 and cdna.aligned_bases > genome_winner.aligned_bases:
        return "recovered_clip"
    return None


@dataclass
class Verdict:
    n_eligible: int = 0
    n_cdna_wins: int = 0
    n_wins_beating_ultra: int = 0
    n_wins_no_cause: int = 0
    net_hp_ed_reduction: float = 0.0
    decisions: List[WinDecision] = field(default_factory=list)

    @property
    def trivial_win_leak(self) -> float:
        return (self.n_wins_no_cause / self.n_cdna_wins) if self.n_cdna_wins else 0.0

    @property
    def ultra_specific_win_frac(self) -> float:
        return (self.n_wins_beating_ultra / self.n_cdna_wins) if self.n_cdna_wins else 0.0

    def go(self) -> bool:
        """GO criteria (Phase 0, §10 + falsifiers §4.2):
          - a non-trivial subset wins (>0),
          - those wins beat uLTRA specifically (majority),
          - >95% of wins are attributable (<=5% trivial-win leak),
          - measurable net HP-ED reduction (>0).
        """
        return (
            self.n_cdna_wins > 0
            and self.net_hp_ed_reduction > 0
            and self.trivial_win_leak <= 0.05
            and self.ultra_specific_win_frac > 0.5
        )

    def summary(self) -> str:
        v = "GO" if self.go() else "NO-GO"
        return (
            f"[Phase-0 verdict: {v}]\n"
            f"  eligible (weak) reads      : {self.n_eligible}\n"
            f"  cDNA wins                  : {self.n_cdna_wins}\n"
            f"  wins beating uLTRA-specific: {self.n_wins_beating_ultra} "
            f"({self.ultra_specific_win_frac:.0%})\n"
            f"  wins with NO attributable cause (leak): {self.n_wins_no_cause} "
            f"({self.trivial_win_leak:.1%})  [NO-GO if >5%]\n"
            f"  net HP-ED reduction        : {self.net_hp_ed_reduction:.1f}\n"
        )


def run_verdict(triples, *, min_junction_anchor_bp: int = 10, epsilon: float = EPSILON) -> Verdict:
    """triples: iterable of (cdna_cand, genome_winner_cand, ultra_cand_or_None)."""
    v = Verdict()
    for cdna, gwin, ultra in triples:
        v.n_eligible += 1
        d = decide(cdna, gwin, ultra=ultra,
                   min_junction_anchor_bp=min_junction_anchor_bp, epsilon=epsilon)
        v.decisions.append(d)
        if d.cdna_wins:
            v.n_cdna_wins += 1
            v.net_hp_ed_reduction += (d.hp_ed_genome - d.hp_ed_cdna)
            if d.cause is None:
                v.n_wins_no_cause += 1
            if d.beats_ultra:
                v.n_wins_beating_ultra += 1
    return v
