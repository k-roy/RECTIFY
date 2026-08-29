#!/usr/bin/env python3
"""Cross-sample QC over a cDNA cohort: PCA + correlation heatmap + a QC-metric table.

WHY THIS IS SEPARATE FROM ``correct-cdna``
------------------------------------------
Per-sample QC (read types, sub-types, UMI duplication) is emitted by every Stage-1 run —
see ``rectify/core/cdna/qc.py``. PCA and clustering are **inherently cross-sample**: they
need N libraries to compare, so they cannot live inside a single-sample command. This
script is the cohort pass; it consumes the ``stage1_qc.json`` sidecars Stage 1 writes.

🔴 **WHAT THIS PCA IS, AND WHAT IT IS NOT.** Stage 1 runs *before* gene assignment (that
happens in Stage 3 ``cdna-analyze``), so there are no expression counts here. This PCA is
over **library QC metrics** — duplication, type composition, XF tiers, tail-length shape.
It answers *"is one library technically unlike its replicates?"* — batch effects, a failed
prep, a mis-assigned barcode. It is **NOT** the expression PCA you would show in a paper;
for that, run ``rectify cdna-analyze`` and PCA the CPA-site/gene count matrix. Both are
worth doing and they answer different questions. Do not substitute one for the other.

Usage
-----
    python scripts/cdna_cohort_qc.py OUTDIR STAGE1_DIR [STAGE1_DIR ...]
    python scripts/cdna_cohort_qc.py OUTDIR --glob '/path/to/*/stage1'

Outputs ``cohort_qc.tsv`` (one row per sample), and — if matplotlib/sklearn are present —
``cohort_qc_pca.png`` and ``cohort_qc_heatmap.png``.
"""
from __future__ import annotations

import argparse
import glob as _glob
import json
import math
import sys
from pathlib import Path
from typing import Dict, List

# Metrics that describe library technical quality. Deliberately excludes raw depth:
# depth differences are expected between samples and would otherwise dominate PC1.
FEATURES = [
    ("type1_read_frac", lambda q: q["read_type"].get("type1_read_frac")),
    ("type1_cluster_frac", lambda q: q["read_type"].get("type1_cluster_frac")),
    ("umi_duplication_rate", lambda q: q.get("umi_duplication_rate")),
    ("mean_reads_per_molecule", lambda q: q.get("mean_reads_per_molecule")),
    ("xf_ge1_frac", lambda q: q["xf_tier_clusters"].get("xf_ge1_frac")),
    ("frame_flip_frac", lambda q: q["pretrim_health"]["frame_flipped"]
     / max(1, q["output_molecules"])),
    ("noop_3p_frac", lambda q: q["pretrim_health"]["noop_3p"]
     / max(1, q["output_molecules"])),
    ("rdna_masked_frac", lambda q: q["dropped"]["rdna_masked_reads"]
     / max(1, q["input_reads"])),
]


def load(paths: List[Path]) -> List[Dict]:
    out = []
    for p in paths:
        f = p / "stage1_qc.json" if p.is_dir() else p
        if not f.exists():
            print(f"  skip (no stage1_qc.json): {p}", file=sys.stderr)
            continue
        q = json.loads(f.read_text())
        q["_sample"] = q.get("sample") or f.parent.parent.name
        q["_path"] = str(f)
        out.append(q)
    return out


def tail_vector(q: Dict) -> List[float]:
    h = q.get("polya_tail_hist_xf2", {})
    tot = sum(h.values()) or 1
    return [v / tot for v in h.values()]


def pearson(a: List[float], b: List[float]) -> float:
    n = len(a)
    if n < 2:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    va = math.sqrt(sum((x - ma) ** 2 for x in a))
    vb = math.sqrt(sum((x - mb) ** 2 for x in b))
    if va == 0 or vb == 0:
        return float("nan")
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (va * vb)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("stage1_dirs", nargs="*", type=Path)
    ap.add_argument("--glob", default=None, help="glob for stage1 dirs")
    a = ap.parse_args()

    paths = list(a.stage1_dirs)
    if a.glob:
        paths += [Path(p) for p in _glob.glob(a.glob)]
    qcs = load(paths)
    if not qcs:
        print("no stage1_qc.json found", file=sys.stderr)
        return 2
    a.outdir.mkdir(parents=True, exist_ok=True)

    names = [q["_sample"] for q in qcs]
    rows = []
    for q in qcs:
        r = {"sample": q["_sample"], "input_reads": q["input_reads"],
             "output_molecules": q["output_molecules"],
             "shallow": q.get("shallow_library"),
             "flags": ";".join(q.get("flags", [])) or "-"}
        for k, fn in FEATURES:
            try:
                r[k] = fn(q)
            except Exception:
                r[k] = None
        rows.append(r)

    cols = (["sample", "input_reads", "output_molecules", "shallow"]
            + [k for k, _ in FEATURES] + ["flags"])
    tsv = a.outdir / "cohort_qc.tsv"
    with tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    print(f"wrote {tsv}  ({len(rows)} samples)")

    for r in rows:
        if r["flags"] != "-":
            print(f"  FLAG {r['sample']}: {r['flags']}")

    # correlation over the poly(A) tail-shape profile — a compositional fingerprint that
    # separates libraries by prep behaviour rather than by depth
    vecs = {q["_sample"]: tail_vector(q) for q in qcs}
    corr = a.outdir / "cohort_tailshape_corr.tsv"
    with corr.open("w") as fh:
        fh.write("sample\t" + "\t".join(names) + "\n")
        for n1 in names:
            fh.write(n1 + "\t" + "\t".join(f"{pearson(vecs[n1], vecs[n2]):.4f}"
                                           for n2 in names) + "\n")
    print(f"wrote {corr}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib/numpy unavailable — TSVs written, plots skipped", file=sys.stderr)
        return 0

    feat_names = [k for k, _ in FEATURES]
    M = np.array([[float(r[k]) if r.get(k) is not None else np.nan for k in feat_names]
                  for r in rows], dtype=float)
    colmean = np.nanmean(M, axis=0)
    idx = np.where(np.isnan(M))
    M[idx] = np.take(colmean, idx[1])
    sd = M.std(axis=0)
    sd[sd == 0] = 1.0
    Z = (M - M.mean(axis=0)) / sd

    if Z.shape[0] >= 3:
        U, S, Vt = np.linalg.svd(Z - Z.mean(axis=0), full_matrices=False)
        PC = U[:, :2] * S[:2]
        var = (S ** 2) / max(1e-12, (S ** 2).sum())
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(PC[:, 0], PC[:, 1], s=60)
        for i, n in enumerate(names):
            ax.annotate(n, (PC[i, 0], PC[i, 1]), fontsize=7,
                        xytext=(4, 3), textcoords="offset points")
        ax.set_xlabel(f"PC1 ({100 * var[0]:.0f}% var)")
        ax.set_ylabel(f"PC2 ({100 * var[1]:.0f}% var)")
        ax.set_title("cDNA Stage-1 QC-metric PCA\n(technical QC, NOT expression)", fontsize=10)
        fig.tight_layout()
        fig.savefig(a.outdir / "cohort_qc_pca.png", dpi=150)
        print(f"wrote {a.outdir / 'cohort_qc_pca.png'}")

    C = np.corrcoef(np.array([vecs[n] for n in names]))
    fig, ax = plt.subplots(figsize=(1 + 0.5 * len(names), 1 + 0.45 * len(names)))
    im = ax.imshow(C, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7)
    ax.set_title("poly(A) tail-shape correlation", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(a.outdir / "cohort_qc_heatmap.png", dpi=150)
    print(f"wrote {a.outdir / 'cohort_qc_heatmap.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
