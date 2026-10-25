"""Run one abPOA call in a throwaway child process.

WHY THIS EXISTS. abPOA sizes its DP allocation from the partial-order graph, and for some inputs it
asks for a single enormous block — measured at exactly 2^38 bytes (256 GiB) — and calls
``posix_memalign``. On failure its C code prints

    [SIMDMalloc] posix_memalign fail!
    Size: 274877906944, Error: ENOMEM

and **aborts the process**. Inside ``correct-cdna``'s ``ProcessPoolExecutor`` that kills a pool
worker, the parent sees only ``concurrent.futures.process.BrokenProcessPool``, and the whole library
fails. Two libraries of a 35-library campaign were lost this way (Chanfreau ``[[1190]]`` tasks 8
and 28); the real cause is invisible in the Python traceback because the child printed it before the
parent's traceback began.

WHY NOT A CAP. The obvious fix is to refuse oversized clusters, and it does not work:

  * A **cluster-size** cap is the wrong axis — 16,384 sequences of 2,000 nt complete fine, while
    FOUR sequences of 100,000 nt reproduce the failure byte for byte.
  * A **length** cap has no viable value. Across 31 completed libraries (101,100,207 molecules) the
    longest consensus abPOA produced successfully is **286,560 nt**, yet a synthetic 4 x 100,000 nt
    probe dies. So no threshold is simultaneously above what real data needs and below where abPOA
    breaks: 20 kb would silently change 19 of 31 libraries, 50 kb would change 7, 100 kb would
    change 3. The request tracks graph size, i.e. sequence DIVERGENCE, not length alone.

So this module catches the abort instead of predicting it. Evidence: Chanfreau
``planning/1195_abpoa_probe.md`` (jobs 15266506, 15266567).

PROTOCOL. stdin: UTF-8 JSON list of sequence strings. stdout: the consensus sequence, raw bytes, or
nothing at all when abPOA declined to produce one. Exit 0 means the stdout payload is authoritative
(empty payload = no consensus). Any other exit, or death by signal, means the caller must fall back.
Diagnostics go to stderr so the abPOA message survives in the run log.

Run as:  python -m rectify.core.cdna._poa_child
"""
from __future__ import annotations

import json
import sys


def main() -> int:
    try:
        seqs = json.loads(sys.stdin.buffer.read().decode("utf-8"))
    except Exception as exc:                                  # noqa: BLE001
        sys.stderr.write(f"_poa_child: unreadable stdin payload: {exc}\n")
        return 64
    if not isinstance(seqs, list) or not all(isinstance(s, str) for s in seqs):
        sys.stderr.write("_poa_child: payload is not a list of strings\n")
        return 65
    if len(seqs) < 2:
        return 0                                              # nothing to align; empty payload

    try:
        import pyabpoa
    except ImportError:
        sys.stderr.write("_poa_child: pyabpoa is not importable in this interpreter\n")
        return 66

    # 🔴 No try/except around the msa() call can help: abPOA aborts the PROCESS rather than raising,
    # which is exactly why this runs in a child. The guard that matters is the process boundary.
    aligner = pyabpoa.msa_aligner()
    res = aligner.msa(seqs, out_cons=True, out_msa=False)
    if getattr(res, "cons_seq", None):
        sys.stdout.buffer.write(res.cons_seq[0].encode("ascii"))
    sys.stdout.buffer.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
