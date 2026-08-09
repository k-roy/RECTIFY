#!/usr/bin/env python
"""Discovery-loss panel: does the microhomology guard veto REAL cryptics that sit in
microhomology?  Quantifies discovery-loss rate vs (microhom_drift_margin, drift_near_tie_cap),
binned by delta_improve.  Go/no-go on whether margin=3 suffices or a positional signal is needed.

Real cryptic construction: EXON2 = U (k-mer) + U' (U with `mm` mismatches) + TAIL (random).
Incumbent acceptor ne = start of U (canonical GT..AG decoy); CRYPTIC acceptor je = ne+k = start
of U' (the real non-canonical site).  Microhomology U~U' (frac = 1 - mm/k).  A read genuinely from
the cryptic carries exon2 bases = genome[je:], but is PLACED at the incumbent → the refiner should
move ne->je (discovery); the guard may veto it.

Author: Kevin R. Roy (inline, per advisor: panel-first is step 1 of building the close).
"""
from __future__ import annotations
import random
import sys
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from rectify.core.splice.junction_refiner import (
    refine_read_junctions, _move_microhomology,
)
from rectify.core.splice.junction_scoring import _score_junction, _EQ, _N

CHROM = "chrTEST"
BASES = "ACGT"


def _rand_seq(n, rng):
    return "".join(rng.choice(BASES) for _ in range(n))


def _mutate(s, nmm, rng):
    """Return s with exactly nmm substituted positions (distinct base each)."""
    s = list(s)
    if nmm <= 0:
        return "".join(s)
    pos = rng.sample(range(len(s)), min(nmm, len(s)))
    for p in pos:
        s[p] = rng.choice([b for b in BASES if b != s[p]])
    return "".join(s)


def _ont_errors(seq, sub_rate, indel_rate, rng):
    """Apply a simple ONT-like error model (substitutions + short indels)."""
    out = []
    for b in seq:
        r = rng.random()
        if r < indel_rate * 0.5:
            continue                                  # deletion
        if r < indel_rate:
            out.append(b); out.append(rng.choice(BASES))   # insertion
            continue
        if rng.random() < sub_rate:
            out.append(rng.choice([x for x in BASES if x != b]))
        else:
            out.append(b)
    return "".join(out)


def build_case(k, mm, tail_len, e1, e2_tail_take, rng):
    """Build a genome + a cryptic-origin read placed at the incumbent.  Returns
    (genome, ns, ne, je, read, q_split) or None if geometry invalid."""
    LPAD = 40
    EXON1_LEN = 60
    INTRON = "GT" + _rand_seq(86, rng) + "AG"          # 90 bp canonical intron
    U = _rand_seq(k, rng)
    Uprime = _mutate(U, mm, rng)
    TAIL = _rand_seq(tail_len, rng)
    exon2 = U + Uprime + TAIL + _rand_seq(40, rng)     # extra so windows stay in-bounds
    genome = "T" * LPAD + _rand_seq(EXON1_LEN, rng) + INTRON + exon2 + "T" * 40
    ns = LPAD + EXON1_LEN                               # donor (fixed)
    ne = ns + len(INTRON)                               # incumbent acceptor = start of U
    je = ne + k                                         # cryptic acceptor = start of U'
    # read: exon1 (last e1 bases of exon1) + N-op[ns,ne) + exon2 = genome[je:je+E2]
    #       (the read carries CRYPTIC exon2 bases but is aligned at the incumbent)
    E2 = k + e2_tail_take                               # U' (k) + part of TAIL
    exon1_bases = genome[ns - e1:ns]
    cryptic_exon2 = genome[je:je + E2]                  # read genuinely FROM the cryptic je
    # fab (canonical-origin) read: exon2 from the incumbent ne (= U + U' + tail...). Its TRUE
    # site is ne; a drift to je is a FABRICATION.  Same length window past je for fairness.
    canon_exon2 = genome[ne:ne + k + E2]
    read_start = ns - e1
    return genome, ns, ne, je, exon1_bases, cryptic_exon2, canon_exon2, read_start


def make_read(genome, exon1_bases, cryptic_exon2, read_start, ns, ne):
    q = exon1_bases + cryptic_exon2
    header = pysam.AlignmentHeader.from_references([CHROM], [len(genome)])
    r = pysam.AlignedSegment(header)
    r.query_name = "cryptic_read"
    r.reference_id = 0
    r.reference_start = read_start
    r.mapping_quality = 60
    r.is_unmapped = False
    r.is_reverse = False
    r.query_sequence = q
    # CIGAR: exon1 = M(len exon1), N(intron), exon2 = M(len exon2). M (op 0) allows mismatch.
    r.cigartuples = [(0, len(exon1_bases)), (_N, ne - ns), (0, len(cryptic_exon2))]
    return r, len(exon1_bases)   # q_split = len(exon1_bases)


def _semiglobal_ed(query, ref):
    """Edit distance aligning ALL of `query` to a PREFIX of `ref` (free ref suffix; no free
    query prefix — hard-anchored at ref[0]).  Indel-tolerant.  Lower = query matches ref better."""
    m, n = len(query), len(ref)
    if m == 0:
        return 0
    prev = list(range(n + 1))          # aligning empty query: 0 cost to consume any ref prefix? No:
    # row 0 = deleting query chars only; but we FORCE query fully consumed and ref free-suffix.
    # Standard semi-global: dp[i][j] = min edits to align query[:i] to ref[:j].
    # Answer = min over j of dp[m][j] (query consumed, ref up to any j = free suffix).
    prev = [0] * (n + 1)               # dp[0][j] = 0 (free ref prefix? NO — hard anchor at ref[0])
    # Hard anchor at ref[0]: query[0] must align at ref[0], so ref has NO free prefix.
    # dp[0][j] = j? no. Let me do: dp[i][0] = i (query prefix vs empty ref = all insertions).
    #            dp[0][j] = j (empty query vs ref prefix = all deletions) BUT we take free suffix
    #            at the END only, so intermediate dp[0][j]=j is correct; free suffix = min last row.
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        cur = [i] + [0] * n
        qi = query[i - 1]
        for j in range(1, n + 1):
            cost = 0 if qi == ref[j - 1] else 1
            cur[j] = min(prev[j - 1] + cost, prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return min(prev)                   # free ref suffix: query fully consumed, ref up to any j


def ed_signal(genome, ns, ne, je, q, q_split, W=28):
    """Indel-robust positional-distinctiveness signal: hard-anchored edit distance of the read's
    rescue to the incumbent exon2 vs the cryptic exon2.  signal = ed(rescue, genome[ne:]) -
    ed(rescue, genome[je:]).  >0 ⇒ read matches the CRYPTIC better (real discovery); <0 ⇒ matches
    the INCUMBENT better (error-driven fab drift).  No free-prefix split ⇒ removes the scorer's
    soft-clip escape that hides the discriminating mismatches."""
    rescue = q[q_split:q_split + W]
    n = len(genome)
    ref_inc = genome[ne:min(n, ne + W + 6)]
    ref_cry = genome[je:min(n, je + W + 6)]
    return _semiglobal_ed(rescue, ref_inc) - _semiglobal_ed(rescue, ref_cry)


def disc_signal(genome, ns, ne, je, q, q_split, W=24):
    """Positional-distinctiveness signal (naive per-index; indel-sensitive by design — this
    is the orthogonality CHECK, not the production impl).  At read offsets r where the two
    candidate exon2s DIVERGE (genome[je+r] != genome[ne+r]), tally whether the read base votes
    CRYPTIC (matches genome[je+r]) or INCUMBENT (matches genome[ne+r]).  A real cryptic should
    vote cryptic (>0); an error-driven fab drift should vote incumbent (<=0)."""
    rescue = q[q_split:q_split + W]
    cry = inc = ndisc = 0
    n = len(genome)
    for r in range(len(rescue)):
        if je + r >= n or ne + r >= n:
            break
        gj, gn = genome[je + r], genome[ne + r]
        if gj == gn:
            continue
        ndisc += 1
        if rescue[r] == gj:
            cry += 1
        elif rescue[r] == gn:
            inc += 1
    return cry - inc, ndisc


def acceptor_after_refine(read, genome, ns, ne, je, **guard_kw):
    pool = {CHROM: [(ns, ne), (ns, je)]}
    reps = refine_read_junctions(
        read, pool, set(), genome, "+", motif_blind=True,
        boundary_error_window=0, **guard_kw,
    )
    acc = ne
    for _idx, ons, one, njs, nje in reps:
        if ons == ns and one == ne:
            acc = nje
    return acc


def run_panel(n_per_cell=40, with_errors=True, seed=1):
    rng = random.Random(seed)
    MARGINS = [3.0, 8.0]
    CAPS = [0.0, 1.0, 2.0, 3.0]
    grid = []
    # sweep microhomology strength (k, mm) and tail-distinctiveness (tail_len, e2_tail_take)
    for k in (4, 6, 8, 10):
        for mm in range(0, k // 2 + 1):                # 0..k/2 mismatches → mh 1.0..0.5
            for e2_tail_take in (2, 4, 8, 16):
                grid.append((k, mm, e2_tail_take))
    rows = []
    for (k, mm, e2_tail_take) in grid:
        for _rep in range(n_per_cell):
            g, ns, ne, je, e1b, cx, canx, rs = build_case(
                k, mm, tail_len=30, e1=40, e2_tail_take=e2_tail_take, rng=rng)
            mh = _move_microhomology(g, ns, ne, ns, je)
            for origin, exon2 in (("cry", cx), ("fab", canx)):
                e2b = _ont_errors(exon2, 0.06, 0.03, rng) if with_errors else exon2
                read, q_split = make_read(g, e1b, e2b, rs, ns, ne)
                q = read.query_sequence
                s_inc, _ = _score_junction(q, q_split, ns, ne, g, 0.25, 15, 10, current_ns=ns)
                s_cry, _ = _score_junction(q, q_split, ns, je, g, 0.25, 15, 10, current_ns=ns)
                delta = s_inc - s_cry
                sig, ndisc = disc_signal(g, ns, ne, je, q, q_split)
                esig = ed_signal(g, ns, ne, je, q, q_split)
                # guard OFF: does the refiner move ne->je?  cry: discovery (good); fab: drift (bad)
                moved_off = acceptor_after_refine(read, g, ns, ne, je) == je
                row = dict(origin=origin, k=k, mm=mm, e2=e2_tail_take, mh=mh, delta=delta,
                           sig=sig, esig=esig, ndisc=ndisc, moved_off=moved_off)
                for m in MARGINS:
                    for c in CAPS:
                        acc = acceptor_after_refine(read, g, ns, ne, je,
                                                    microhom_drift_margin=m, drift_near_tie_cap=c)
                        row[f"moved_m{int(m)}_c{int(c)}"] = (acc == je)
                # TRUE end-to-end SHIPPED operating point through the real refiner: margin=3,
                # cap=2, positional_gate=1 (the full close)
                row["moved_close"] = (acceptor_after_refine(
                    read, g, ns, ne, je, microhom_drift_margin=3.0,
                    drift_near_tie_cap=2.0, drift_positional_gate=1.0) == je)
                rows.append(row)
    return rows, MARGINS, CAPS


def _dist(ds):
    ds = sorted(ds)
    if not ds:
        return "n=0"
    def pct(p): return ds[min(len(ds) - 1, int(p * len(ds)))]
    return (f"n={len(ds)} min={ds[0]:.2f} p10={pct(.1):.2f} median={pct(.5):.2f} "
            f"p90={pct(.9):.2f} max={ds[-1]:.2f} | <1:{sum(d<1 for d in ds)/len(ds):.0%} "
            f"<2:{sum(d<2 for d in ds)/len(ds):.0%} <3:{sum(d<3 for d in ds)/len(ds):.0%}")


def summarize(rows, MARGINS, CAPS):
    out = []
    out.append(f"total reads: {len(rows)}")
    # DISCOVERY side: real cryptics that are discoverable guard-OFF AND trip mh>=0.5 (at-risk)
    disc_at_risk = [r for r in rows if r["origin"] == "cry" and r["moved_off"] and r["mh"] >= 0.5]
    # FAB side: canonical reads that DRIFT guard-OFF (error-driven FP) AND trip mh>=0.5
    fab_drift = [r for r in rows if r["origin"] == "fab" and r["moved_off"] and r["mh"] >= 0.5]
    fab_total = [r for r in rows if r["origin"] == "fab" and r["mh"] >= 0.5]
    out.append(f"DISCOVERY at-risk (cry, moved-OFF, mh>=0.5): {len(disc_at_risk)}")
    out.append(f"  delta_improve dist: {_dist([r['delta'] for r in disc_at_risk])}")
    out.append(f"FAB drift (canon, moved-OFF=drifted, mh>=0.5): {len(fab_drift)} "
               f"of {len(fab_total)} canon-in-mh (guard-OFF fab rate {len(fab_drift)/max(1,len(fab_total)):.1%})")
    out.append(f"  delta_improve dist: {_dist([r['delta'] for r in fab_drift])}")
    out.append("")
    out.append("TRADEOFF per (margin,cap):  discovery-LOSS = cry at-risk now HELD (bad);")
    out.append("                            fab-RESIDUAL   = fab drift still MOVED (bad, want low)")
    out.append(f"  {'m':>3} {'cap':>3}   {'disc-loss':>9}   {'fab-residual':>12}")
    for m in MARGINS:
        for c in CAPS:
            key = f"moved_m{int(m)}_c{int(c)}"
            loss = sum(not r[key] for r in disc_at_risk) / max(1, len(disc_at_risk))   # cry no longer moved
            fabr = sum(r[key] for r in fab_drift) / max(1, len(fab_drift))             # fab still moved
            out.append(f"  {int(m):>3} {int(c):>3}   {loss:>9.1%}   {fabr:>12.1%}")
    # ORTHOGONALITY CHECK: within the delta overlap band, does the discriminating-column signal
    # separate real cryptics (should vote >0) from fab drifts (should vote <=0)?
    out.append("")
    out.append("ORTHOGONALITY of the positional signal (disc-column vote) within the delta OVERLAP band:")
    for lo, hi in ((0.0, 2.0), (0.5, 1.5)):
        cry_band = [r for r in disc_at_risk if lo <= r["delta"] <= hi]
        fab_band = [r for r in fab_drift if lo <= r["delta"] <= hi]
        out.append(f"  delta in [{lo},{hi}]:  cry n={len(cry_band)}  fab n={len(fab_band)}")
        out.append(f"    cry signal: {_dist([r['sig'] for r in cry_band])}")
        out.append(f"    fab signal: {_dist([r['sig'] for r in fab_band])}")
        if cry_band and fab_band:
            # separability at signal>0: what frac of each is on the correct side
            cry_pos = sum(r["sig"] > 0 for r in cry_band) / len(cry_band)
            fab_le0 = sum(r["sig"] <= 0 for r in fab_band) / len(fab_band)
            out.append(f"    sep @sig>0:  cry votes cryptic {cry_pos:.0%} | fab votes incumbent(<=0) {fab_le0:.0%}")
    # INDEL-ROBUST signal separation (align-first ed_signal) vs naive, in the overlap band
    out.append("")
    out.append("INDEL-ROBUST signal (ed_signal, align-first) separation in delta band [0.5,1.5]:")
    cryb = [r for r in disc_at_risk if 0.5 <= r["delta"] <= 1.5]
    fabb = [r for r in fab_drift if 0.5 <= r["delta"] <= 1.5]
    for name, fld in (("naive sig", "sig"), ("ed_signal", "esig")):
        if cryb and fabb:
            cp = sum(r[fld] > 0 for r in cryb) / len(cryb)
            fl = sum(r[fld] <= 0 for r in fabb) / len(fabb)
            out.append(f"    {name:10s}: cry votes cryptic {cp:.0%} | fab votes incumbent {fl:.0%} "
                       f"| balanced-acc {(cp+fl)/2:.0%}")
    # DECISIVE: cap-alone vs cap+positional (naive vs ed) at m=3, cap=2
    out.append("")
    out.append("CAP-ALONE vs CAP+POSITIONAL (m=3, cap=2; positional spares a veto when signal>thr):")
    key = "moved_m3_c2"
    for label, fld, thr in (("cap-alone", None, None),
                            ("cap+naive(sig>0)", "sig", 0),
                            ("cap+ed(esig>0)", "esig", 0),
                            ("cap+ed(esig>=2)", "esig", 1)):
        def moved(r):
            if r[key]:
                return True
            if fld is not None and r["moved_off"] and r[fld] > thr:
                return True
            return False
        loss = sum(not moved(r) for r in disc_at_risk) / max(1, len(disc_at_risk))
        fabr = sum(moved(r) for r in fab_drift) / max(1, len(fab_drift))
        out.append(f"    {label:18s}  disc-loss {loss:6.1%}   fab-residual {fabr:6.1%}")
    # TRUE end-to-end shipped operating point (margin=3, cap=2, positional_gate=1) via real refiner
    cl_loss = sum(not r["moved_close"] for r in disc_at_risk) / max(1, len(disc_at_risk))
    cl_fabr = sum(r["moved_close"] for r in fab_drift) / max(1, len(fab_drift))
    out.append(f"    {'WIRED m3/c2/gate1':18s}  disc-loss {cl_loss:6.1%}   fab-residual {cl_fabr:6.1%}  <- shipped close")
    return "\n".join(out)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--no-errors", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rows, M, C = run_panel(n_per_cell=a.n, with_errors=not a.no_errors, seed=a.seed)
    print(summarize(rows, M, C))
