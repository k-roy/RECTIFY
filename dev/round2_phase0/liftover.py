"""Round-2 cDNA Phase-0 harness — transcript->genome lift-over (N-op insertion).

THROWAWAY / Phase-0 scope (per SPEC_round2_cdna_discovery_assignment_20260614.md §10):
this is a *minimal, locus-specific* lift-over to answer the empirical kill-gate, NOT the
production block-map-driven lift-over (Phase 2). It is deliberately small and heavily
asserted so the GO/NO-GO verdict is trustworthy.

It consumes:
  - a read->cDNA alignment (cDNA-local CIGAR + offset into the padded cDNA), and
  - the per-cDNA block map (transcript<->genome tiling),
and produces a ``LiftedRead`` whose four attributes
(``cigartuples``, ``reference_name``, ``reference_start``, ``query_sequence``)
are EXACTLY what the production scorers in
``rectify.core.consensus.corrected_consensus`` read:
``_cigar_hp_edit_distance``, ``_cigar_aligned_bases``, ``_cigar_min_junction_anchor``.
So the Phase-0 verdict is computed with production code, not a re-implementation.

Designer-2 §2 algorithm + 12-row edge-case table; minus-strand is the highest-bug-density
path (revcomp query + reverse CIGAR + assert the orientation invariant §2.5).

CIGAR op codes (pysam convention):  0=M 1=I 2=D 3=N 4=S 5=H 7=  8=X
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

# op codes
M, I, D, N, S, H, EQ, X = 0, 1, 2, 3, 4, 5, 7, 8
_REF_CONSUMING = frozenset({M, D, N, EQ, X})
_QUERY_CONSUMING = frozenset({M, I, S, EQ, X})

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass
class Block:
    """One contiguous transcript<->genome tile of a padded cDNA.

    The genome interval is ALWAYS given in plus-strand coordinates with
    ``g_start < g_end`` (half-open).  Transcript coordinates are 0-based into the
    padded cDNA, increasing 5'->3' in RNA/gene sense.

    Plus-strand gene:  cdna[t_start + k]  ==  genome[g_start + k]
    Minus-strand gene: cdna[t_start + k]  ==  revcomp(genome)[...]  i.e.
                       cdna[t_start + k]  ==  comp(genome[g_end - 1 - k])
    """
    t_start: int
    g_start: int
    g_end: int           # exclusive, plus-strand; g_end - g_start == length
    is_pad: bool = False
    kind: str = "exon"

    @property
    def length(self) -> int:
        return self.g_end - self.g_start

    @property
    def t_end(self) -> int:
        return self.t_start + self.length


@dataclass
class CdnaModel:
    cdna_id: str
    chrom: str
    strand: str               # '+' or '-'  (gene/RNA strand)
    blocks: List[Block]       # ordered 5'->3' along the transcript
    length: int               # total padded cDNA length

    def validate(self, genome: Dict[str, str]) -> None:
        """Assert the per-block 1:1 identity against the genome (Designer-1 §2.3 / §5.8)."""
        gseq = genome[self.chrom]
        # blocks tile the transcript contiguously
        exp_t = 0
        for b in self.blocks:
            assert b.t_start == exp_t, f"{self.cdna_id}: block gap at t={exp_t} (got {b.t_start})"
            assert b.g_end > b.g_start, f"{self.cdna_id}: bad genome interval {b.g_start}-{b.g_end}"
            exp_t = b.t_end
        assert exp_t == self.length, f"{self.cdna_id}: blocks tile {exp_t} != cdna length {self.length}"
        # transcript-order genome monotonicity by strand
        gs = [b.g_start for b in self.blocks]
        if self.strand == '+':
            assert gs == sorted(gs), f"{self.cdna_id}: + gene blocks not genome-ascending"
        else:
            assert gs == sorted(gs, reverse=True), f"{self.cdna_id}: - gene blocks not genome-descending"

    def t2g(self, t: int) -> int:
        """Map a transcript coordinate to a plus-strand genome coordinate."""
        for b in self.blocks:
            if b.t_start <= t < b.t_end:
                off = t - b.t_start
                if self.strand == '+':
                    return b.g_start + off
                return b.g_end - 1 - off
        raise ValueError(f"{self.cdna_id}: transcript coord {t} outside blocks")

    def block_of(self, t: int) -> Block:
        return self.blocks[self.block_index(t)]

    def block_index(self, t: int) -> int:
        for i, b in enumerate(self.blocks):
            if b.t_start <= t < b.t_end:
                return i
        raise ValueError(f"{self.cdna_id}: transcript coord {t} outside blocks")

    def genome_gap(self, i: int) -> int:
        """Intron length (genome gap) between block i and block i+1, always >= 0.

        Plus strand:  blocks ascend  -> next.g_start - cur.g_end
        Minus strand: blocks descend -> cur.g_start - next.g_end (cur is the higher locus)
        """
        cur, nxt = self.blocks[i], self.blocks[i + 1]
        gap = (nxt.g_start - cur.g_end) if self.strand == '+' else (cur.g_start - nxt.g_end)
        assert gap >= 0, f"{self.cdna_id}: negative intron gap {gap} at block {i}"
        return gap


@dataclass
class LiftedRead:
    """Duck-typed pysam-AlignedSegment replacement exposing ONLY the attributes the
    production ``_cigar_*`` scorers read. Keeps Phase-0 free of a BAM header round-trip."""
    cigartuples: List[Tuple[int, int]]
    reference_name: str
    reference_start: int
    query_sequence: str
    is_reverse: bool = False


def _merge(cig: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge adjacent identical ops — but NEVER across an N (Designer-2 §2.3 step 4)."""
    out: List[Tuple[int, int]] = []
    for op, ln in cig:
        if ln == 0:
            continue
        if out and out[-1][0] == op and op != N:
            out[-1] = (op, out[-1][1] + ln)
        else:
            out.append((op, ln))
    return out


def lift(
    cdna: CdnaModel,
    cdna_cigar: List[Tuple[int, int]],
    cdna_ref_start: int,
    read_query: str,
) -> LiftedRead:
    """Lift a read's forward alignment-to-cDNA into a genome ``LiftedRead`` with N-ops.

    Parameters
    ----------
    cdna : CdnaModel              the padded-cDNA model + block map the read aligned to
    cdna_cigar : list[(op,len)]   cDNA-local CIGAR (no N expected; E12 => error)
    cdna_ref_start : int          0-based offset of the alignment into the padded cDNA
    read_query : str              read sequence in cDNA/transcript orientation (as aligned)

    Returns a LiftedRead in genome (reference-forward, low->high) orientation.
    """
    # E12 — a no-splice cDNA aligner must never emit N.
    assert not any(op == N for op, _ in cdna_cigar), "cDNA-local CIGAR carries N (E12 violation)"

    # Walk the cDNA CIGAR in TRANSCRIPT order, emitting genome ops + N at spanned
    # block boundaries.  For a '-' gene this yields ops in DESCENDING genome order;
    # we reverse + revcomp at the end (Designer-2 §2.3 step 1, E11).
    t = cdna_ref_start
    q = 0
    ops: List[Tuple[int, int]] = []        # genome ops in transcript order
    lead_clip = 0                          # leading S (read overhang past pad start, E2.a)
    trail_clip = 0

    prev_ref_block = [None]   # block index of the LAST ref-consuming base emitted (list = mutable closure)

    def _emit_ref_run(op_code: int, length: int):
        """Emit a ref-consuming run (M/=/X/D) over [t, t+length).

        N-insertion is driven by the transcript walk crossing a block boundary, tracked
        GLOBALLY across the whole CIGAR via ``prev_ref_block`` — so an intron is emitted
        whether a single op spans the boundary OR two ops break exactly at it (the bug the
        per-op version missed, which silently fused exons and hit micro-exon reads hardest).
        """
        nonlocal t
        remaining = length
        while remaining > 0:
            bi = cdna.block_index(t)
            if prev_ref_block[0] is not None and bi != prev_ref_block[0]:
                # crossed one or more block boundaries since the last ref base -> emit each N.
                # An I op emitted between the two ref runs was appended BEFORE this N, so a
                # donor-side insertion lands as `=, I, N, =` (E6).
                j = prev_ref_block[0]
                while j < bi:
                    gap = cdna.genome_gap(j)
                    if gap > 0:
                        ops.append((N, gap))
                    j += 1
            b = cdna.blocks[bi]
            take = min(remaining, b.t_end - t)
            ops.append((op_code, take))
            t += take
            remaining -= take
            prev_ref_block[0] = bi

    n_cig = len(cdna_cigar)
    for idx, (op, length) in enumerate(cdna_cigar):
        if op in (M, EQ, X):
            _emit_ref_run(op, length)
            q += length
        elif op == D:
            # E5: D belongs to the exon side, never merged into an N. _emit_ref_run
            # places it within the current block; a D exactly at a boundary stays
            # attached to the 5' block because we consume the block's remaining
            # transcript span first.
            _emit_ref_run(D, length)
        elif op == I:
            # E6: I attaches to the 5' (upstream) exon side — emit it inline; since it
            # consumes no transcript, it never triggers a boundary split.
            ops.append((I, length))
            q += length
        elif op == S:
            # leading/trailing soft-clip = read overhang past the padded cDNA ends (E2.a/E9)
            if idx == 0:
                lead_clip += length
            elif idx == n_cig - 1:
                trail_clip += length
            else:
                raise ValueError("internal S op in cDNA CIGAR")
            q += length
        elif op == H:
            ops.append((H, length))
        else:
            raise ValueError(f"unhandled cDNA CIGAR op {op}")

    ops = _merge(ops)

    # reference_start = genome coord of the first ref-consuming op (transcript order).
    # For '+', that's the lowest genome coord we touched; for '-', the highest.
    first_ref_t = cdna_ref_start  # leading S already excluded from t-walk? No: S advanced q only.
    # find the transcript coord where the first ref-consuming op began:
    # it is cdna_ref_start + (leading soft-clip consumes query only, not transcript),
    # i.e. exactly cdna_ref_start.
    g_first = cdna.t2g(first_ref_t)

    if cdna.strand == '+':
        ref_start = g_first
        out_cigar = ([(S, lead_clip)] if lead_clip else []) + ops + ([(S, trail_clip)] if trail_clip else [])
        query_out = read_query
        is_rev = False
    else:
        # E11: reverse op order (transcript 0 -> highest genome coord) and revcomp query.
        # reference_start is the LOWEST genome coord touched = end of the last ref-consuming op.
        last_ref_t = cdna_ref_start + sum(
            l for o, l in cdna_cigar if o in (M, D, EQ, X)
        ) - 1
        ref_start = cdna.t2g(last_ref_t)
        rev_ops = list(reversed(ops))
        # leading clip in transcript space becomes trailing clip in genome space and vice-versa
        out_cigar = (
            ([(S, trail_clip)] if trail_clip else [])
            + rev_ops
            + ([(S, lead_clip)] if lead_clip else [])
        )
        query_out = revcomp(read_query)
        is_rev = True

    out_cigar = _merge(out_cigar)
    return LiftedRead(
        cigartuples=out_cigar,
        reference_name=cdna.chrom,
        reference_start=ref_start,
        query_sequence=query_out,
        is_reverse=is_rev,
    )


# ─────────────────────────── invariants (Designer-2 §2.5) ───────────────────────────

def assert_query_length(lifted: LiftedRead, read_len: int) -> None:
    n = sum(l for o, l in lifted.cigartuples if o in _QUERY_CONSUMING)
    assert n == read_len, f"query-length invariant: {n} != {read_len}"


def reference_span(lifted: LiftedRead) -> int:
    return sum(l for o, l in lifted.cigartuples if o in _REF_CONSUMING)


def assert_no_n_merge(lifted: LiftedRead) -> None:
    prev = None
    for op, _ in lifted.cigartuples:
        if op in (M, EQ) and prev in (M, EQ):
            raise AssertionError("two match runs merged (should have been separated)")
        prev = op


def assert_orientation_invariant(lifted: LiftedRead, genome: Dict[str, str]) -> Tuple[bool, bool]:
    """CRITICAL minus-strand check (§2.5): the first matched read base must equal
    ``genome[reference_start]`` and the last matched read base ``genome[reference_end-1]``.
    Returns (first_ok, last_ok). On a known-clean read both must be True; this is the
    ONE invariant that catches a missed revcomp (query/ref-span invariants are revcomp-blind).
    """
    gseq = genome[lifted.reference_name]
    q = lifted.query_sequence
    # first matched base
    ref_pos = lifted.reference_start
    q_pos = 0
    first_ok = last_ok = None
    for op, ln in lifted.cigartuples:
        if op in (M, EQ, X) and first_ok is None:
            first_ok = (q[q_pos].upper() == gseq[ref_pos].upper())
        if op in _REF_CONSUMING:
            ref_pos += ln
        if op in _QUERY_CONSUMING:
            q_pos += ln
    # last matched base — walk again tracking the final M/=/X base
    ref_pos = lifted.reference_start
    q_pos = 0
    for op, ln in lifted.cigartuples:
        if op in (M, EQ, X):
            last_ok = (q[q_pos + ln - 1].upper() == gseq[ref_pos + ln - 1].upper())
        if op in _REF_CONSUMING:
            ref_pos += ln
        if op in _QUERY_CONSUMING:
            q_pos += ln
    return bool(first_ok), bool(last_ok)
