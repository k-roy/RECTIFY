"""Phase-1 junction-proximal indel realignment (planning/776; Kevin's two-phase
architecture, 2026-08-24).

Phase 1 assigns the read AS OBSERVED to its ed-minimal split alignment —
objective (plain-Levenshtein) minimization, motif- and annotation-BLIND, up to
``MAXD`` net indel bases and boundary/split shifts of ``SHIFT`` around each
existing N-op. Junction ADMISSION is deliberately not this module's job:
Station C (core/consensus/station_c.py) arbitrates junctions over pooled
evidence, gating on junction-proximal overhang cleanliness. Run WITHOUT that
gate, this module will happily move reads onto systematic-error
pseudo-junctions — planning/776b measured 29 recurrent ones, every single one
canonical->non-canonical jitter of an annotated RPG intron. That is BY DESIGN:
per-read edit distance says nothing about isoforms. Never wire this into a
pipeline whose junction consumers do not pass through Station C.

DEV STATUS (2026-08-24): standalone command, NOT registered in the resolver,
run-all, or any production path. Production use awaits Kevin's sign-off.

    python -m rectify.core.align.indel_realigner IN.bam OUT.bam \
        --genome GENOME.fa [--x 50 --shift 6 --maxd 4 --min-gain 1]

Input contract: a BAM straight from ``rectify align`` (``samtools calmd -e``
applied, CIGARs unmodified since calmd) — '=' bytes in SEQ are decoded under
the record's OWN placement. 🔴 Do NOT feed overhang-resolver output: rescued
blocks' '=' bytes were written under the PRE-rescue placement and decode
wrongly at the rescued coordinates ([[772]] defect 4). After a rewrite this
module RE-ENCODES '=' under the new placement and recomputes NM (MD dropped),
so the calmd contract survives the rewrite.

Scoring model (validated in planning/776b's scan, gate-tested in 776c):
    s      query-split shift        a, b   intron start / end shifts
    dl, dr net query-indel window compensation per side, |dl|+|dr| <= MAXD
Per-side global Levenshtein over flank-clamped windows (Xl/Xr <= X), combined
additively. A candidate wins only if it beats the emitted alignment's window
cost by >= min_gain AND its winning value is achieved by exactly ONE distinct
junction after ``same_junction`` collapse (planning/770's ambiguity
discipline). The rewrite realigns the ENTIRE flanking M-blocks anchored at
their outer edges (reference_start and downstream ops never move) and is
accepted only if the full-block edit total also improves by >= min_gain.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import logging
import sys
import time
from typing import Dict, List, Optional, Tuple

import pysam

from ..splice.overhang_informativeness import canonicalize_junction, same_junction
from ...utils.genome import load_genome, standardize_chrom_name
from ...config import CHROM_TO_GENOME

try:
    from rapidfuzz.distance import Levenshtein as _Lev
except ImportError as _e:  # pragma: no cover
    _Lev = None

logger = logging.getLogger("rectify.indel_realigner")

# The phase-1 metadata tag. Audited free on feat/resolver-v3-773 @ db975a7:
# not set anywhere in the tree and absent from cma_schema's reserved list
# (XI = cdna isoform, XU = UMI, XV = existing aux, XQ/XW = 773's collision
# lesson — every one of those was "obviously free" until grepped).
TAG = "XP"

MIN_INTRON = 20


def decode_query(cigartuples, ref_start: int, seq: str, chrom: str) -> str:
    """Literal query bases, '=' resolved under the record's OWN placement."""
    out: List[str] = []
    qp = 0
    rp = ref_start
    for op, ln in cigartuples:
        if op in (0, 7, 8):
            for t in range(ln):
                c = seq[qp + t]
                out.append(chrom[rp + t] if c == "=" else c)
            qp += ln
            rp += ln
        elif op in (1, 4):
            out.append(seq[qp:qp + ln])
            qp += ln
        elif op in (2, 3):
            rp += ln
    return "".join(out)


def _iter_junctions(ct, ref_start: int):
    """(i, d, e, qp, m_l, m_r) for every N-op flanked by plain-M ops."""
    qp = 0
    rp = ref_start
    for i, (op, ln) in enumerate(ct):
        if op == 3 and 0 < i < len(ct) - 1 \
                and ct[i - 1][0] == 0 and ct[i + 1][0] == 0:
            yield i, rp, rp + ln, qp, ct[i - 1][1], ct[i + 1][1]
        if op in (0, 1, 4, 7, 8):
            qp += ln
        if op in (0, 2, 3, 7, 8):
            rp += ln


def scan(dq: str, qp: int, chrom: str, d: int, e: int, Xl: int, Xr: int,
         shift: int, maxd: int, min_gain: int, stats) -> Optional[dict]:
    """776b scoring: return the unique winning candidate, or None."""
    L = len(chrom)
    ql = len(dq)
    cur = sum(1 for t in range(Xl) if dq[qp - Xl + t] != chrom[d - Xl + t]) \
        + sum(1 for t in range(Xr) if dq[qp + t] != chrom[e + t])
    if cur < min_gain:  # nothing to gain (cur==0 for min_gain 1)
        stats["clean"] += 1
        return None

    s_lo = max(-shift, Xl - qp)
    s_hi = min(shift, ql - qp - Xr)
    if s_lo > 0 or s_hi < 0:
        s_lo = s_hi = 0

    cutoff = cur - min_gain

    def side_tab(left: bool, dlt: int):
        tab: Dict[int, Dict[int, int]] = {}
        for s in range(s_lo, s_hi + 1):
            q = dq[qp + s - Xl:qp + s] if left else dq[qp + s:qp + s + Xr]
            row: Dict[int, int] = {}
            for sh in range(-shift, shift + 1):
                if left:
                    lo, hi = d + sh - Xl + dlt, d + sh
                else:
                    lo, hi = e + sh, e + sh + Xr - dlt
                if lo < 0 or hi > L or hi <= lo:
                    continue
                v = _Lev.distance(q, chrom[lo:hi], score_cutoff=cutoff)
                if v <= cutoff:
                    row[sh] = v
            if row:
                tab[s] = row
        return tab

    Lt = {dl: side_tab(True, dl) for dl in range(-maxd, maxd + 1)}
    Rt = {dr: side_tab(False, dr) for dr in range(-maxd, maxd + 1)}
    pairs = [(dl, dr) for dl in range(-maxd, maxd + 1)
             for dr in range(-maxd, maxd + 1) if abs(dl) + abs(dr) <= maxd]

    best = None
    cands: List[Tuple[int, int, int, int, int, int]] = []  # (v, a, b, dl, dr, s)
    for dl, dr in pairs:
        for s, Ls in Lt[dl].items():
            Rs = Rt[dr].get(s)
            if not Rs:
                continue
            for a, vl in Ls.items():
                if vl > cutoff:
                    continue
                for b, vr in Rs.items():
                    v = vl + vr
                    if v > cutoff:
                        continue
                    if (e + b) - (d + a) < MIN_INTRON:
                        continue
                    cands.append((v, a, b, dl, dr, s))
                    if best is None or v < best:
                        best = v
    if best is None:
        stats["no_gain"] += 1
        return None

    ties = [c for c in cands if c[0] == best]
    distinct: List[Tuple[int, int]] = []
    for _v, a, b, _dl, _dr, _s in ties:
        j = (d + a, e + b)
        if not any(same_junction(chrom, j, g) for g in distinct):
            distinct.append(j)
    if len(distinct) > 1:
        stats["ambiguous"] += 1
        return None
    # deterministic, conservative representative: closest to the current placement
    ties.sort(key=lambda c: (abs(c[1]) + abs(c[2]), abs(c[5]),
                             abs(c[3]) + abs(c[4]), c[1:]))
    v, a, b, dl, dr, s = ties[0]
    return dict(cur=cur, win=v, gain=cur - v, a=a, b=b, dl=dl, dr=dr, s=s,
                moved=not same_junction(chrom, (d, e), (d + a, e + b)))


def _opcodes_to_cigar(q: str, r: str) -> Tuple[List[Tuple[int, int]], int]:
    """Levenshtein alignment of q->r as cigar ops [(op, len)] + edit count.
    opcodes tags: equal/replace -> M; delete (extra query) -> I; insert -> D."""
    ops: List[Tuple[int, int]] = []
    edits = 0
    for oc in _Lev.opcodes(q, r):
        if oc.tag in ("equal", "replace"):
            ops.append((0, oc.src_end - oc.src_start))
            if oc.tag == "replace":
                edits += oc.src_end - oc.src_start
        elif oc.tag == "delete":
            ops.append((1, oc.src_end - oc.src_start))
            edits += oc.src_end - oc.src_start
        elif oc.tag == "insert":
            ops.append((2, oc.dest_end - oc.dest_start))
            edits += oc.dest_end - oc.dest_start
    merged: List[Tuple[int, int]] = []
    for op, ln in ops:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + ln)
        else:
            merged.append((op, ln))
    return merged, edits


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def rebuild(ct, i: int, dq: str, chrom: str, d: int, e: int, qp: int,
            m_l: int, m_r: int, winner: dict, min_gain: int, stats):
    """Rewrite flank blocks + N for one junction.

    Each flank block realigns over its NATURAL ref span (query length minus its
    net indel), so the alignment carries exactly the net indels and no spurious
    edge ops. The OUTER ref edge of a side may therefore move by
    ``dl_shift = a - s + dl`` (left) / ``dr_shift = b - s - dr`` (right) —
    allowed only when nothing beyond that side consumes reference (the same
    frame-safety constraint arm 1 enforces); a moving left edge shifts
    reference_start. Returns (new_ct, ref_start_delta, old_edits, new_edits,
    nlen) or None."""
    a, b, s = winner["a"], winner["b"], winner["s"]
    dl, dr = winner["dl"], winner["dr"]
    qL_len = m_l + s
    qR_len = m_r - s
    if qL_len < 1 or qR_len < 1:
        stats["guard_block_empty"] += 1
        return None
    dl_shift = a - s + dl               # movement of the left OUTER ref edge
    dr_shift = b - s - dr               # movement of the right OUTER ref edge
    prefix_ref_free = all(op in (1, 4, 5, 6) for op, _ in ct[:i - 1])
    suffix_ref_free = all(op in (1, 4, 5, 6) for op, _ in ct[i + 2:])
    if (dl_shift != 0 and not prefix_ref_free) or \
            (dr_shift != 0 and not suffix_ref_free):
        stats["guard_frame"] += 1
        return None
    lstart = (d + a) - (qL_len - dl)
    rend = (e + b) + (qR_len - dr)
    if lstart < 0 or rend > len(chrom):
        stats["guard_bounds"] += 1
        return None
    qL, rL = dq[qp - m_l:qp + s], chrom[lstart:d + a]
    qR, rR = dq[qp + s:qp + m_r], chrom[e + b:rend]
    if not rL or not rR:
        stats["guard_block_empty"] += 1
        return None
    old_edits = _hamming(dq[qp - m_l:qp], chrom[d - m_l:d]) \
        + _hamming(dq[qp:qp + m_r], chrom[e:e + m_r])
    fragL, edL = _opcodes_to_cigar(qL, rL)
    fragR, edR = _opcodes_to_cigar(qR, rR)
    if old_edits - (edL + edR) < min_gain:
        stats["guard_full_block"] += 1
        return None
    # ops abutting the N must be M so the junction stays base-anchored
    if fragL[-1][0] != 0 or fragR[0][0] != 0:
        stats["guard_edge_op"] += 1
        return None
    nlen = (e + b) - (d + a)
    if nlen < MIN_INTRON:
        stats["guard_intron"] += 1
        return None
    new_ct = list(ct[:i - 1]) + fragL + [(3, nlen)] + fragR + list(ct[i + 2:])
    merged: List[Tuple[int, int]] = []
    for op, ln in new_ct:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + ln)
        else:
            merged.append((op, ln))
    return merged, dl_shift, old_edits, edL + edR, nlen


def _reencode_and_nm(read, dq: str, chrom: str) -> int:
    """Re-encode '=' under the CURRENT (new) placement and recompute NM."""
    out: List[str] = []
    nm = 0
    qp = 0
    rp = read.reference_start
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            for t in range(ln):
                if dq[qp + t] == chrom[rp + t]:
                    out.append("=")
                else:
                    out.append(dq[qp + t])
                    nm += 1
            qp += ln
            rp += ln
        elif op == 1:
            out.append(dq[qp:qp + ln])
            nm += ln
            qp += ln
        elif op == 4:
            out.append(dq[qp:qp + ln])
            qp += ln
        elif op == 2:
            nm += ln
            rp += ln
        elif op == 3:
            rp += ln
    qual = read.query_qualities
    read.query_sequence = "".join(out)
    if qual is not None:
        read.query_qualities = qual
    return nm


def realign_read(read, chrom: str, cfg: dict, stats, log_rows) -> bool:
    if read.query_sequence is None or not read.cigartuples:
        return False
    dq = decode_query(read.cigartuples, read.reference_start, read.query_sequence, chrom)
    changed = False
    done: List[Tuple[int, int]] = []
    for _round in range(6):
        hit = None
        for i, d, e, qp, m_l, m_r in _iter_junctions(read.cigartuples, read.reference_start):
            if any(abs(d - dd) <= cfg["shift"] + cfg["maxd"]
                   and abs(e - ee) <= cfg["shift"] + cfg["maxd"] for dd, ee in done):
                continue
            hit = (i, d, e, qp, m_l, m_r)
            break
        if hit is None:
            break
        i, d, e, qp, m_l, m_r = hit
        done.append((d, e))
        stats["junctions"] += 1
        Xl = min(cfg["x"], m_l, qp, d)
        Xr = min(cfg["x"], m_r, len(dq) - qp, len(chrom) - e)
        if Xl < 1 or Xr < 1:
            stats["zero_flank"] += 1
            continue
        w = scan(dq, qp, chrom, d, e, Xl, Xr,
                 cfg["shift"], cfg["maxd"], cfg["min_gain"], stats)
        if w is None:
            continue
        rb = rebuild(read.cigartuples, i, dq, chrom, d, e, qp, m_l, m_r,
                     w, cfg["min_gain"], stats)
        if rb is None:
            continue
        new_ct, ref_start_delta, old_ed, new_ed, nlen = rb
        read.cigartuples = new_ct
        if ref_start_delta:
            read.reference_start += ref_start_delta
        d2, e2 = d + w["a"], e + w["b"]
        done.append((d2, e2))
        entry = f"{d}-{e}>{d2}-{e2}:{old_ed}>{new_ed}"
        prev = str(read.get_tag(TAG)) if read.has_tag(TAG) else ""
        read.set_tag(TAG, f"{prev},{entry}" if prev else entry)
        stats["rewritten"] += 1
        stats[f"gain_{min(w['gain'], 5)}"] += 1
        if w["moved"]:
            stats["moved"] += 1
        log_rows.append(dict(read=read.query_name, chrom=read.reference_name,
                             src=f"{d}-{e}", win=f"{d2}-{e2}", nlen=nlen,
                             cur=w["cur"], winv=w["win"], gain=w["gain"],
                             s=w["s"], a=w["a"], b=w["b"], dl=w["dl"], dr=w["dr"],
                             old_block_edits=old_ed, new_block_edits=new_ed,
                             moved=w["moved"]))
        changed = True
    if changed:
        nm = _reencode_and_nm(read, dq, chrom)
        read.set_tag("NM", nm)
        if read.has_tag("MD"):
            read.set_tag("MD", None)
    return changed


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("in_bam")
    ap.add_argument("out_bam")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--x", type=int, default=50)
    ap.add_argument("--shift", type=int, default=6)
    ap.add_argument("--maxd", type=int, default=4)
    ap.add_argument("--min-gain", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after N reads with junctions (smoke)")
    args = ap.parse_args(argv)
    if _Lev is None:
        sys.exit("FATAL: rapidfuzz is required for the phase-1 indel realigner")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    cfg = dict(x=args.x, shift=args.shift, maxd=args.maxd, min_gain=args.min_gain)
    genome = load_genome(args.genome)
    t0 = time.time()
    stats = collections.Counter()
    log_rows: List[dict] = []
    with pysam.AlignmentFile(args.in_bam, check_sq=False) as inb, \
            pysam.AlignmentFile(args.out_bam, "wb", template=inb) as outb:
        n_junc_reads = 0
        for read in inb:
            stats["reads"] += 1
            if (not read.is_unmapped and not read.is_secondary
                    and not read.is_supplementary
                    and read.cigartuples
                    and any(o == 3 for o, _ in read.cigartuples)):
                chrom_name = standardize_chrom_name(read.reference_name)
                key = chrom_name if chrom_name in genome \
                    else CHROM_TO_GENOME.get(chrom_name or "", "")
                if key in genome:
                    n_junc_reads += 1
                    if realign_read(read, genome[key], cfg, stats, log_rows):
                        stats["reads_rewritten"] += 1
            outb.write(read)
            if args.limit and n_junc_reads >= args.limit:
                stats["limited"] = 1
                for rest in inb:
                    outb.write(rest)
                break
    with gzip.open(args.out_bam + ".xp.jsonl.gz", "wt") as fh:
        for row in log_rows:
            fh.write(json.dumps(row) + "\n")
    with open(args.out_bam + ".xp_stats.json", "w") as fh:
        json.dump(dict(stats), fh, indent=1, sort_keys=True)
    logger.info("done in %.0fs: %s", time.time() - t0,
                {k: stats[k] for k in ("reads", "junctions", "rewritten",
                                       "moved", "ambiguous") if k in stats})
    return 0


if __name__ == "__main__":
    sys.exit(main())
