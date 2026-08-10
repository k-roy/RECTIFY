"""Per-library UMI error-profile calibration from the fixed-T spacers (658).

The PCB114 UMI is ``TT VVVV TT VVVV TT VVVV TT VVVV TTT``: 11 positions are fixed ``T`` by
construction, so they are ground truth in every read and the library's UMI error profile is
measurable directly, with no clustering assumptions (Chanfreau planning/650). Error is
stratified by the read's sequencing-time REGIME — ``orient ⊕ is_reverse`` = UMI at pore entry
vs exit — which sets a ~1.7x rate difference (entry dirtier). Runs as a capped pre-pass
(~1-2 min per library) and is recorded as a per-sample sidecar for provenance/QC and for
downstream calibrated merge rules.
"""
from __future__ import annotations

from collections import Counter
from typing import Dict, Optional

import pysam

from ._constants import SSP_FWD, SSP_RC, UMI_LEN
from .read_info import revcomp

FIXED_T = (0, 1, 6, 7, 12, 13, 18, 19, 24, 25, 26)
MAX_SHIFT = 3
EXT = MAX_SHIFT


def _spacer_mismatches(u: str) -> int:
    return sum(1 for j in FIXED_T if u[j] != "T")


def _classify(ext: str):
    """clean / frameshift / sub_k / frameshift_sub / corrupt for one UMI window."""
    m = {}
    for s in range(-MAX_SHIFT, MAX_SHIFT + 1):
        w = ext[EXT + s: EXT + s + UMI_LEN]
        if len(w) == UMI_LEN:
            m[s] = _spacer_mismatches(w)
    m0 = m.get(0, 99)
    if m0 == 0:
        return "clean", 0
    best_s = min((s for s in m if s != 0), key=lambda s: (m[s], abs(s)), default=None)
    if best_s is not None and m[best_s] == 0:
        return "frameshift", best_s
    if m0 <= 2 and (best_s is None or m0 <= m[best_s]):
        return f"sub_{m0}", 0
    if best_s is not None and m[best_s] <= 2 and m[best_s] < m0:
        return "frameshift_sub", best_s
    return "corrupt", 0


def derive_umi_error_profile(bam_path: str, limit: int = 500_000) -> Dict:
    """Scan up to `limit` UMI-bearing primary reads; return the per-regime profile dict.

    Coordinate-sorted scan bias is acceptable: the error regime is a read-level basecall
    property, not a locus property (planning/650).
    """
    cells: Dict[tuple, dict] = {}

    def cell(orient: str, is_rev: bool) -> dict:
        key = (orient, is_rev)
        if key not in cells:
            cells[key] = dict(n=0, cat=Counter(),
                              inframe_pos_nonT=[0] * UMI_LEN, inframe_n=0,
                              qsum=0.0, qn=0)
        return cells[key]

    n_primary = n_umi = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            n_primary += 1
            orient = None
            p = seq.find(SSP_FWD)
            if p >= 0 and p + len(SSP_FWD) + UMI_LEN + EXT <= len(seq):
                u_start = p + len(SSP_FWD)
                if u_start - EXT >= 0:
                    ext = seq[u_start - EXT: u_start + UMI_LEN + EXT]
                    orient = "fwd"
                    umi_span = (u_start, u_start + UMI_LEN)
            if orient is None:
                p = seq.find(SSP_RC)
                if p >= UMI_LEN + EXT and p + EXT <= len(seq):
                    ext = revcomp(seq[p - UMI_LEN - EXT: p + EXT])
                    orient = "rev"
                    umi_span = (p - UMI_LEN, p)
            if orient is None:
                continue
            n_umi += 1
            c = cell(orient, read.is_reverse)
            c["n"] += 1
            cat, _shift = _classify(ext)
            c["cat"][cat] += 1
            if cat == "clean" or cat.startswith("sub_"):
                c["inframe_n"] += 1
                u = ext[EXT: EXT + UMI_LEN]
                for j in FIXED_T:
                    if u[j] != "T":
                        c["inframe_pos_nonT"][j] += 1
            quals = read.query_qualities
            if quals is not None:
                qs = quals[umi_span[0]: umi_span[1]]
                if len(qs):
                    c["qsum"] += sum(qs) / len(qs)
                    c["qn"] += 1
            if n_umi >= limit:
                break

    def regime_of(orient: str, is_rev: bool) -> str:
        return "entry" if (orient == "fwd") == (not is_rev) else "exit"

    out: Dict = dict(bam=str(bam_path), n_primary_scanned=n_primary, n_umi=n_umi,
                     cells={}, regimes={})
    for (orient, is_rev), c in sorted(cells.items()):
        n = max(1, c["n"])
        fs = c["cat"]["frameshift"] + c["cat"]["frameshift_sub"]
        out["cells"][f"{orient}_rev{int(is_rev)}"] = dict(
            orient=orient, is_reverse=is_rev, regime=regime_of(orient, is_rev),
            n=c["n"], clean=c["cat"]["clean"] / n, frameshift=fs / n,
            corrupt=c["cat"]["corrupt"] / n,
            mean_umi_Q=(c["qsum"] / c["qn"]) if c["qn"] else None,
            inframe_n=c["inframe_n"],
            inframe_pos_subrate={str(j): c["inframe_pos_nonT"][j] / max(1, c["inframe_n"])
                                 for j in FIXED_T},
        )
    for reg in ("entry", "exit"):
        ns = [c for k, c in cells.items() if regime_of(*k) == reg]
        if not ns:
            continue
        n = sum(c["n"] for c in ns)
        infn = sum(c["inframe_n"] for c in ns)
        sub = sum(sum(c["inframe_pos_nonT"][j] for j in FIXED_T) for c in ns)
        out["regimes"][reg] = dict(
            n=n,
            clean=sum(c["cat"]["clean"] for c in ns) / max(1, n),
            frameshift=sum(c["cat"]["frameshift"] + c["cat"]["frameshift_sub"]
                           for c in ns) / max(1, n),
            corrupt=sum(c["cat"]["corrupt"] for c in ns) / max(1, n),
            inframe_sub_per_base=sub / max(1, infn * len(FIXED_T)),
            mean_umi_Q=(sum(c["qsum"] for c in ns) / max(1, sum(c["qn"] for c in ns))),
        )
    return out
