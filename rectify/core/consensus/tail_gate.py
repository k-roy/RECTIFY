"""ISSUE-087: a forced 3' tail beats a soft clip only when it fits the genome better than chance.

Consensus scoring charges an RNA-3' soft clip (non-poly(A) bases, -2 per base, cap 10) but charges nothing when
another aligner places the same bases with mismatches or deletions, so a forced tail always beats the clip. On the
2026-09-28 route smoke (wt_by4742_rep1 DRS, chrI, 30 %) the landing extended 925 reads; tails of 16 nt or more that
fit no better than chance were the one class it lost on balance against the library's clean reads (46 pin-better,
21 landing-better), and one 16-nt tail was written as a 5,875-nt deletion that moved the 3' end 5.9 kb.

The gate, applied per read among arms on the same contig and strand. Arm A is the one with the longest RNA-3' soft
clip (``clip``, 16 to 150 nt); every arm E with a shorter clip aligns some of those bases. The clip is tested
against the genome just past A's aligned end with the audit's chance model (dev/prod_20260928/chance_3prime_ext.py):
the best score in bits (match +2, mismatch -2, gap open -2, gap base -0.5; mismatches and gaps inside a genome
homopolymer run of 5 or more charged half) of a placement that consumes the whole clip, left-anchored at the aligned
end (S1), and of a placement that ends with the clip's last base (S2), each against composition-preserving shuffles
of the clip over the same window (clip + 12 nt). Two modes (``RECTIFY_TAIL_GATE``):

* ``filed``: the clip fails the chance test (p >= 0.01 on both S1 and S2). Costed on the route smoke: reverts 69,
  fixes 46, breaks 21 (net +25 of 925); keeps the positive control P3-20, reverts P3-3.
* ``written``: what E drew over those bases (flat scoring, no homopolymer forgiveness) scores below the 99th
  percentile of the S1 null. Reverts 65, fixes 45, breaks 18 (net +27); reverts P3-20.

In both modes a deletion of 100 nt or more inside E's tail always gates (only one read on the route smoke; fixes 1,
breaks 0). A gated arm is charged the clip arm's 3' clip, so its forced tail confers no advantage and the usual
tie-breaks decide. Unset (the default): no change. The random draws are seeded by the read name, so a decision
does not depend on read order or thread count.
"""

from __future__ import annotations

import hashlib
import os
from typing import Dict, Optional, Tuple

import numpy as np

BM, BX, BGO, BGE = 2.0, -2.0, -2.0, -0.5
HP_RUN, EXTRA, NEG = 5, 12, -1e9
MIN_TAIL, MAX_TAIL = 16, 150
P_CHANCE = 0.01
BIG_DELETION = 100
_ENC = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate('ACGT'):
    _ENC[ord(_b)] = _i
_COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def gate_mode() -> Optional[str]:
    """'filed', 'written', or None (off: the default) from ``RECTIFY_TAIL_GATE``."""
    mode = os.environ.get('RECTIFY_TAIL_GATE', '').strip().lower()
    return mode if mode in ('filed', 'written') else None


def n_shuffles() -> int:
    try:
        return max(50, int(os.environ.get('RECTIFY_TAIL_GATE_SHUFFLES', '200')))
    except ValueError:
        return 200


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _enc(s: str) -> np.ndarray:
    return _ENC[np.frombuffer(s.encode('ascii'), dtype=np.uint8)]


def clip3(read) -> int:
    """The RNA-3' soft clip of a DRS/RNA-sense record: right end on +, left end on -."""
    ct = read.cigartuples or []
    if not ct:
        return 0
    op, n = ct[0] if read.is_reverse else ct[-1]
    return n if op == 4 else 0


def decoded_seq(read, chrom_seq: str) -> str:
    s = list(read.query_sequence or '')
    if '=' in s:
        for q, p in read.get_aligned_pairs(matches_only=True):
            if s[q] == '=' and 0 <= p < len(chrom_seq):
                s[q] = chrom_seq[p]
    return ''.join(s).upper()


def window(read, chrom_seq: str, n: int) -> Tuple[str, np.ndarray]:
    """The RNA-orientation genome window of n bases just past the aligned 3' end, and its homopolymer mask (runs
    are measured with 8 extra bases on each side, so a run straddling the window edge counts)."""
    L = len(chrom_seq)
    if not read.is_reverse:
        a = read.reference_end
        lo, hi = max(0, a - 8), min(L, a + n + 8)
        s = chrom_seq[lo:hi].upper()
        off = a - lo
    else:
        a = read.reference_start
        lo, hi = max(0, a - n - 8), min(L, a + 8)
        s = _rc(chrom_seq[lo:hi].upper())
        off = hi - a
    mask = np.zeros(len(s), bool)
    i = 0
    while i < len(s):
        j = i
        while j < len(s) and s[j] == s[i]:
            j += 1
        if j - i >= HP_RUN:
            mask[i:j] = True
        i = j
    return s[off:off + n], mask[off:off + n]


def best_full(Q: np.ndarray, ref: str, hp: np.ndarray, local_end: bool = False) -> np.ndarray:
    """Gotoh DP over the rows of Q (N x L codes), as in the audit's model. S1 (default): the whole query consumed,
    left end fixed at ref[0], free reference suffix. S2 (local_end): the alignment ends with the query's last base and
    may start anywhere. Returns the best bits per row."""
    N, L = Q.shape
    R = len(ref)
    r = _enc(ref)
    half = np.where(hp, 0.5, 1.0)
    H = np.full((N, R + 1), NEG)
    D = np.full((N, R + 1), NEG)
    I = np.full((N, R + 1), NEG)
    H[:, 0] = 0.0
    if local_end:
        H[:, :] = 0.0
    else:
        for j in range(1, R + 1):
            c = half[j - 1]
            D[:, j] = np.maximum(H[:, j - 1] + (BGO + BGE) * c, D[:, j - 1] + BGE * c)
    ins_c = np.concatenate([[half[0] if R else 1.0], half])
    for i in range(1, L + 1):
        q = Q[:, i - 1]
        Hn = np.full((N, R + 1), NEG)
        Dn = np.full((N, R + 1), NEG)
        prev = np.maximum(np.maximum(H[:, :-1], D[:, :-1]), I[:, :-1])
        s = np.where(q[:, None] == r[None, :], BM, BX * half[None, :])
        Hn[:, 1:] = prev + s
        In = np.maximum(H + (BGO + BGE) * ins_c[None, :], I + BGE * ins_c[None, :])
        for j in range(1, R + 1):
            c = half[j - 1]
            Dn[:, j] = np.maximum(Hn[:, j - 1] + (BGO + BGE) * c, Dn[:, j - 1] + BGE * c)
        if local_end and i < L:
            Hn = np.maximum(Hn, 0.0)
        H, D, I = Hn, Dn, In
    if local_end:
        return H.max(axis=1)
    return np.maximum(np.maximum(H, D), I).max(axis=1)


def chance_test(clip: str, ref: str, hp: np.ndarray, seed: int, n: int) -> Dict[str, float]:
    """p values of the clip's whole-tail (S1) and end (S2) fits against n composition-preserving shuffles."""
    rng = np.random.default_rng(seed)
    q = _enc(clip)
    Q = np.vstack([q] + [rng.permutation(q) for _ in range(n)])
    s1 = best_full(Q, ref, hp)
    s2 = best_full(Q, ref, hp, local_end=True)
    p1 = (1 + int((s1[1:] >= s1[0] - 1e-9).sum())) / (1 + n)
    p2 = (1 + int((s2[1:] >= s2[0] - 1e-9).sum())) / (1 + n)
    return {'p_whole': p1, 'p_end': p2, 'bits': float(s1[0]), 'null_q99': float(np.quantile(s1[1:], 0.99))}


def drawn_tail(read, chrom_seq: str, seq: str, tail_len: int) -> Tuple[float, int]:
    """Flat bits of what *read* drew over its last *tail_len* RNA-3' query bases, and its longest deletion there."""
    Lq = len(seq)
    pairs = read.get_aligned_pairs()
    if not read.is_reverse:
        first_q = Lq - tail_len
        cols = []
        for q, p in pairs:
            if q is not None and q < first_q:
                cols = []
                continue
            cols.append((q, p))
    else:
        cols = []
        for q, p in pairs:
            if q is not None and q >= tail_len:
                break
            cols.append((q, p))
        cols = cols[::-1]
    c3 = clip3(read)
    bits, prev, run, longest_d = 0.0, None, 0, 0
    for q, p in cols:
        if q is not None and p is not None:
            g = chrom_seq[p] if 0 <= p < len(chrom_seq) else 'N'
            bits += BM if seq[q] == g else BX
            prev, run = 'M', 0
        elif p is None:
            soft = (not read.is_reverse and q >= Lq - c3) or (read.is_reverse and q < c3)
            if soft:
                prev, run = 'S', 0
            else:
                bits += BGE + (BGO if prev != 'I' else 0.0)
                prev, run = 'I', 0
        else:
            bits += BGE + (BGO if prev != 'D' else 0.0)
            run = run + 1 if prev == 'D' else 1
            longest_d = max(longest_d, run)
            prev = 'D'
    return bits, longest_d


def _seed(read_id: str) -> int:
    return int(hashlib.md5(str(read_id).encode()).hexdigest()[:12], 16)


def apply_tail_gate(alignments, reads, genome, mode: Optional[str] = None) -> Dict[str, str]:
    """Charge every arm whose forced 3' tail the gate rejects the clip arm's 3' clip (see the module docstring).

    *alignments*: aligner -> AlignmentInfo (scored); *reads*: aligner -> the decoded pysam record. Returns
    aligner -> reason for each gated arm (empty when nothing is gated or the gate is off).
    """
    mode = mode or gate_mode()
    gated: Dict[str, str] = {}
    if mode is None or not reads or len(reads) < 2:
        return gated
    groups: Dict[Tuple[str, bool], list] = {}
    for name, r in reads.items():
        if r is None or r.is_unmapped or not r.cigartuples or name not in alignments:
            continue
        groups.setdefault((r.reference_name, r.is_reverse), []).append(name)
    for (chrom, _rev), names in groups.items():
        if len(names) < 2:
            continue
        chrom_seq = genome.get(chrom) or ''
        if not chrom_seq:
            continue
        clips = {n: clip3(reads[n]) for n in names}
        a_name = max(names, key=lambda n: (clips[n], n))
        clip_len = clips[a_name]
        if clip_len < MIN_TAIL or clip_len > MAX_TAIL:
            continue
        a_read = reads[a_name]
        a_seq = decoded_seq(a_read, chrom_seq)
        others = [n for n in names if clips[n] < clip_len and len(decoded_seq(reads[n], chrom_seq)) == len(a_seq)]
        if not others:
            continue
        clip = a_seq[-clip_len:] if not a_read.is_reverse else _rc(a_seq[:clip_len])
        ref, hp = window(a_read, chrom_seq, clip_len + EXTRA)
        if len(ref) < clip_len:
            continue
        test = chance_test(clip, ref, hp, _seed(a_read.query_name), n_shuffles())
        a_info = alignments[a_name]
        for e in others:
            e_read = reads[e]
            e_seq = decoded_seq(e_read, chrom_seq)
            bits, longest_d = drawn_tail(e_read, chrom_seq, e_seq, clip_len)
            if longest_d >= BIG_DELETION:
                reason = f'tail_deletion_{longest_d}'
            elif mode == 'filed' and test['p_whole'] >= P_CHANCE and test['p_end'] >= P_CHANCE:
                reason = f"chance_fail_p{test['p_whole']:.3f}_{test['p_end']:.3f}"
            elif mode == 'written' and bits < test['null_q99']:
                reason = f"written_tail_{bits:.1f}_below_{test['null_q99']:.1f}"
            else:
                continue
            e_info = alignments[e]
            old_eff = e_info.effective_three_prime_clip
            new_eff = max(old_eff, a_info.effective_three_prime_clip)
            if new_eff != old_eff:
                e_info.effective_three_prime_clip = new_eff
                # score_alignment charges min(2 x effective 3' clip, 10): charge the difference
                e_info.junction_score -= min(new_eff * 2, 10) - min(old_eff * 2, 10)
            gated[e] = reason
    return gated
