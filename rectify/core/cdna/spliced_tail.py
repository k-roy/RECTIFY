"""Report-only flag for an RNA 3' end placed on a genomic A-run after an intron (the A17 / spliced poly(A)-tail class).

When a molecule is aligned with its poly(A) tail still in SEQ, a splice-aware aligner can place the tail (or the
molecule's last bases plus the tail) as a short terminal "exon" on a downstream genomic A-run, joined by an N
(RPS26B: 744 of 1,675 molecules in the May cDNA tables; the joint atlas's clusters on review cards a1, a2, c1, c2,
c3 and h2). An adapter base that happens to match the genome after the A-run then blocks the walkback, so the false
3' end survives correction. Rectify's current pipeline cuts the tail before alignment and makes almost none
(0.002-0.009 % of cDNA records, 0.0004-0.007 % of direct RNA; dev/audit_20261012/AUDIT.md section 6); BAMs aligned
with their tails (legacy tables, other pipelines) carry 1-2 %.

The rule is the audit's pre-registered one: the RNA-3' terminal block (after the LAST N in the 'fwd' frame, before the
FIRST N in 'rev') holds at most 30 aligned read bases, sits next to an N of at least 30 nt, and the genome under its 8
RNA-3'-most reference positions holds at least 6 of the tail base (A, or T for 'rev'). It almost never matches an
annotated intron, which is how a tail placed by chance behaves. This module only reports; re-anchoring the end at
the donor side of the N waits for Kevin's review of cards (queue r1 D4).
"""

from typing import Optional, Tuple

MAX_BLOCK_ALIGNED = 30
MIN_INTRON = 30
GENOME_WINDOW = 8
MIN_GENOME_TAIL = 6

_REF_OPS = (0, 2, 3, 7, 8)


def terminal_arun_junction(cigartuples, reference_start: int, chrom_seq: str,
                           rna_reverse: bool) -> Optional[Tuple[int, int, int]]:
    """``(donor_end, block_aligned, intron_len)`` when the RNA 3' end is a short block on a genomic A-run after an N.

    ``donor_end`` is the alternative 3' end on the donor side of the N: the last base before it ('fwd'), or the first
    base after it ('rev'). Returns None when the rule does not hold.
    """
    if not cigartuples or not chrom_seq:
        return None
    starts = []
    pos = reference_start
    for op, n in cigartuples:
        starts.append(pos)
        if op in _REF_OPS:
            pos += n
    n_idx = [i for i, (op, _) in enumerate(cigartuples) if op == 3]
    if not n_idx:
        return None
    i_n = n_idx[0] if rna_reverse else n_idx[-1]
    intron_len = cigartuples[i_n][1]
    if intron_len < MIN_INTRON:
        return None
    block = range(0, i_n) if rna_reverse else range(i_n + 1, len(cigartuples))
    aligned = 0
    ref_positions = []
    for i in block:
        op, n = cigartuples[i]
        if op in (0, 7, 8):
            aligned += n
            ref_positions.extend(range(starts[i], starts[i] + n))
        elif op == 2:
            ref_positions.extend(range(starts[i], starts[i] + n))
    if aligned == 0 or aligned > MAX_BLOCK_ALIGNED:
        return None
    ref_positions.sort()
    window = ref_positions[:GENOME_WINDOW] if rna_reverse else ref_positions[-GENOME_WINDOW:]
    tail = "T" if rna_reverse else "A"
    if sum(1 for p in window if 0 <= p < len(chrom_seq) and chrom_seq[p].upper() == tail) < MIN_GENOME_TAIL:
        return None
    n_start = starts[i_n]
    donor_end = n_start + intron_len if rna_reverse else n_start - 1
    return donor_end, aligned, intron_len
