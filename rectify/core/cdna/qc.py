"""Stage-1 (pre-trim) library QC for ONT PCR-cDNA.

Every ``rectify correct-cdna`` run emits this block, regardless of ``--workers``.
That invariance is the point: before this module the serial path printed a read-type
breakdown, XF tiers and a tail-length histogram while the parallel path printed none of
it, because ``_cdna_region_task`` computed those metrics per region and discarded them.
Since parallel is the path any real run takes, production runs silently shipped without
the QC that the docs describe.

The block answers four questions a first pass over a cDNA library must answer:

1. **Read-type composition** — Type-1 (SSP+UMI captured) vs Type-2 (SSP-less,
   5'-truncated), plus the ``XY`` sub-types (``1a`` fwd / ``1b`` rev / ``2``). A Type-1
   read share well below the ~82 % documented in ``algorithms/cdna_correct.md`` means
   heavy truncation (degraded RNA / failed prep) or a chemistry mismatch.
2. **UMI duplication** — how much amplification the library actually carries.
3. **Full-length confidence** (``XF`` tiers) and the poly(A) tail-length distribution.
4. **What was dropped** — rDNA masking, polyA-pileup buckets over ``--per-cluster-cap``,
   region-boundary reads.

🔴 **Read-level vs molecule-level is the trap.** ``type1_reads`` is a READ count; the
number of Type-1 *clusters* is a MOLECULE count. Type-1 reads collapse by UMI, Type-2
reads do not, so the two fractions are NOT interchangeable and a threshold written
against one does not transfer to the other. Both are reported, always labelled.

Depth-awareness matters too: ``umi_duplication_rate`` is meaningless on a shallow
library — at ~50 k reads you rarely sample the same molecule twice, so a near-zero
duplication rate is expected and is not evidence of a dedup failure. ``qc_depth_caveat``
carries that judgement so downstream consumers do not have to re-derive it.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

# Below this many input reads, duplication/saturation metrics are not interpretable.
SHALLOW_READ_THRESHOLD = 250_000
# Documented Type-1 read share for clean PCB114 yeast data (algorithms/cdna_correct.md).
EXPECTED_TYPE1_READ_FRAC = 0.817
# Warn below this; the docs call a much lower fraction a chemistry-mismatch signal.
TYPE1_READ_FRAC_WARN = 0.70

TAIL_BINS = [(0, 1), (1, 10), (10, 20), (20, 30), (30, 50), (50, 75),
             (75, 100), (100, 150), (150, 250), (250, 10_000)]


def _median(xs: List[int]) -> int:
    if not xs:
        return 0
    s = sorted(xs)
    return s[len(s) // 2]


def collect_qc(
    *,
    fastq_stats: Dict,
    stats: Dict,
    n_clusters: int,
    cluster_xf_tier: Optional[Dict[int, int]] = None,
    cluster_tail_len: Optional[Dict[int, int]] = None,
    tier_counts: Optional[Dict[int, int]] = None,
    tail_hist: Optional[List[int]] = None,
    n_input_reads: Optional[int] = None,
    workers: int = 1,
    sample: Optional[str] = None,
) -> Dict:
    """Build the Stage-1 QC record.

    Accepts either the raw per-cluster maps (serial path, which holds them in memory) or
    pre-aggregated ``tier_counts`` / ``tail_hist`` (parallel path, which sums them across
    regions). Either way the emitted record has the same shape, so a consumer never has
    to know which path produced it.
    """
    t1_reads = int(stats.get("type1_reads", 0))
    t2_reads = int(stats.get("type2_reads", 0))
    t1_clusters = int(stats.get("type1_clusters", 0))
    t2_clusters = int(stats.get("type2_clusters", 0))
    typed_reads = t1_reads + t2_reads

    written = int(fastq_stats.get("written", 0))
    in_reads = int(n_input_reads if n_input_reads is not None
                   else fastq_stats.get("input_reads", 0) or typed_reads)

    if tier_counts is None:
        tier_counts = defaultdict(int)
        for v in (cluster_xf_tier or {}).values():
            tier_counts[int(v)] += 1
    tier_counts = {int(k): int(v) for k, v in dict(tier_counts).items()}

    if tail_hist is None:
        tail_hist = [0] * len(TAIL_BINS)
        for cid in range(n_clusters):
            if (cluster_xf_tier or {}).get(cid, 0) < 2:
                continue
            tl = (cluster_tail_len or {}).get(cid, 0)
            for i, (lo, hi) in enumerate(TAIL_BINS):
                if lo <= tl < hi:
                    tail_hist[i] += 1
                    break

    n_clust = max(1, n_clusters)
    # A molecule seen more than once is a PCR duplicate; on a shallow library this is
    # depth-limited, not amplification-limited (see qc_depth_caveat).
    dup_rate = (1.0 - written / in_reads) if in_reads > 0 else 0.0
    mean_cluster_size = (in_reads / written) if written > 0 else 0.0

    sub = {k: int(v) for k, v in (stats.get("subtype_reads") or {}).items()}

    qc: Dict = {
        "sample": sample,
        "workers": workers,
        "input_reads": in_reads,
        "output_molecules": written,
        "consensus_source": {
            "singletons": int(fastq_stats.get("from_singletons", 0)),
            "multi_read_poa_or_pileup": int(fastq_stats.get("from_multi_pileup", 0)),
            "multi_read_rep_fallback": int(fastq_stats.get("from_multi_fallback", 0)),
        },
        "read_type": {
            "type1_reads": t1_reads,
            "type2_reads": t2_reads,
            "type1_read_frac": (t1_reads / typed_reads) if typed_reads else None,
            "type1_clusters": t1_clusters,
            "type2_clusters": t2_clusters,
            "type1_cluster_frac": (t1_clusters / (t1_clusters + t2_clusters))
                                  if (t1_clusters + t2_clusters) else None,
        },
        "subtype_reads": sub,
        "umi_duplication_rate": dup_rate,
        "mean_reads_per_molecule": mean_cluster_size,
        "xf_tier_clusters": {
            "xf2_anchored_high": tier_counts.get(2, 0),
            "xf1_unanchored_medium": tier_counts.get(1, 0),
            "xf0_not_detected": tier_counts.get(0, 0),
            "xf_ge1_frac": (tier_counts.get(1, 0) + tier_counts.get(2, 0)) / n_clust,
        },
        "polya_tail_hist_xf2": {f"{lo}-{hi}": n for (lo, hi), n in zip(TAIL_BINS, tail_hist)},
        "pretrim_health": {
            "frame_flipped": int(fastq_stats.get("trim_frame_flipped", 0)),
            "frame_vs_xo_mismatch": int(fastq_stats.get("trim_frame_mismatch", 0)),
            "noop_5p_type1_only": int(fastq_stats.get("trim_noop_5p", 0)),
            "noop_3p": int(fastq_stats.get("trim_noop_3p", 0)),
        },
        "dropped": {
            "rdna_masked_reads": int(fastq_stats.get("n_rdna_masked",
                                                     stats.get("n_rdna_masked", 0)) or 0),
            "polya_pileup_buckets": int(stats.get("buckets_dropped_polyA_pileup", 0)),
            "reads_in_dropped_buckets": int(stats.get("reads_in_dropped_buckets", 0)),
            "region_boundary_reads": int(fastq_stats.get("n_boundary_dropped", 0)),
        },
    }

    flags: List[str] = []
    shallow = in_reads < SHALLOW_READ_THRESHOLD
    qc["shallow_library"] = shallow
    qc["qc_depth_caveat"] = (
        f"input_reads={in_reads} < {SHALLOW_READ_THRESHOLD}: umi_duplication_rate and "
        "mean_reads_per_molecule are DEPTH-LIMITED, not amplification-limited. A near-zero "
        "duplication rate here is expected and is NOT evidence that dedup failed."
    ) if shallow else None

    t1f = qc["read_type"]["type1_read_frac"]
    if t1f is not None and t1f < TYPE1_READ_FRAC_WARN:
        flags.append(
            f"LOW_TYPE1_READ_FRAC: {t1f:.1%} of typed READS are Type-1 (expected "
            f"~{EXPECTED_TYPE1_READ_FRAC:.0%}). Suggests heavy 5' truncation (degraded RNA "
            "or failed prep) or a chemistry mismatch. NOTE: compare against the READ-level "
            "expectation only — the cluster/molecule fraction is a different quantity."
        )
    if qc["pretrim_health"]["frame_vs_xo_mismatch"] > 0.01 * max(1, written):
        flags.append("FRAME_VS_XO_MISMATCH >1% of molecules — orientation handling suspect.")
    if not shallow and dup_rate < 0.05:
        flags.append(
            f"LOW_UMI_DUPLICATION: {dup_rate:.1%} on a non-shallow library "
            f"({in_reads} reads). Expected amplification duplication; near-zero may mean "
            "UMI extraction or clustering is not working."
        )
    qc["flags"] = flags
    return qc


def render_qc(qc: Dict) -> str:
    """Human-readable Stage-1 QC block. Identical for serial and parallel runs."""
    L: List[str] = []
    a = L.append
    w = qc.get("workers", 1)
    a("=" * 70)
    a(f"cdna_correct v1 — Stage-1 dedup complete"
      f"{f' (parallel {w} workers)' if w and w > 1 else ''}")
    a("=" * 70)
    inr, out = qc["input_reads"], qc["output_molecules"]
    a(f"input reads:                       {inr:>10d}")
    a(f"output records (one per molecule): {out:>10d}"
      + (f"  ({100 * out / inr:.1f}% of input)" if inr else ""))
    cs = qc["consensus_source"]
    a(f"  singletons (passed through):     {cs['singletons']:>10d}")
    a(f"  multi-read (pileup/POA):         {cs['multi_read_poa_or_pileup']:>10d}")
    a(f"  multi-read (rep fallback):       {cs['multi_read_rep_fallback']:>10d}")
    a("")
    a(f"UMI duplication rate:              {qc['umi_duplication_rate']:>10.1%}"
      f"   (mean {qc['mean_reads_per_molecule']:.2f} reads/molecule)")
    if qc.get("shallow_library"):
        a("  ⚠ DEPTH-LIMITED — see qc_depth_caveat; a low rate here is expected, not a failure.")
    a("")
    rt = qc["read_type"]
    a("Read-type breakdown (READS — not interchangeable with clusters below):")
    if rt["type1_read_frac"] is not None:
        a(f"  Type 1 (SSP+UMI captured):       {rt['type1_reads']:>10d}  ({rt['type1_read_frac']:.1%})")
        a(f"  Type 2 (SSP-less, 5'-truncated): {rt['type2_reads']:>10d}  ({1 - rt['type1_read_frac']:.1%})")
    else:
        a("  (no typed reads)")
    if rt["type1_cluster_frac"] is not None:
        a("Read-type breakdown (MOLECULES / clusters):")
        a(f"  Type 1 clusters:                 {rt['type1_clusters']:>10d}  ({rt['type1_cluster_frac']:.1%})")
        a(f"  Type 2 clusters:                 {rt['type2_clusters']:>10d}  ({1 - rt['type1_cluster_frac']:.1%})")
    if qc.get("subtype_reads"):
        a("Sub-type (XY) breakdown (READS):")
        tot = max(1, sum(qc["subtype_reads"].values()))
        for k in sorted(qc["subtype_reads"]):
            v = qc["subtype_reads"][k]
            a(f"  XY={k:<26s} {v:>10d}  ({100 * v / tot:.1f}%)")
    a("")
    x = qc["xf_tier_clusters"]
    a("XF tier breakdown (full-length confidence, clusters):")
    a(f"  XF=2 (anchored, HIGH):           {x['xf2_anchored_high']:>10d}")
    a(f"  XF=1 (unanchored, MEDIUM):       {x['xf1_unanchored_medium']:>10d}")
    a(f"  XF=0 (not detected):             {x['xf0_not_detected']:>10d}")
    a(f"  XF≥1 (any full-length):          "
      f"{x['xf1_unanchored_medium'] + x['xf2_anchored_high']:>10d}  ({x['xf_ge1_frac']:.1%})")
    a("")
    a("PolyA tail-length distribution (per-cluster median, XF=2 only):")
    for k, v in qc["polya_tail_hist_xf2"].items():
        a(f"  {k:>12s} nt: {v:>10d}")
    a("")
    ph = qc["pretrim_health"]
    a("Adapter pretrim health:")
    for label, key in [("frame flipped (basecalled seq)", "frame_flipped"),
                       ("frame vs XO label mismatch", "frame_vs_xo_mismatch"),
                       ("5' trim NO-OP (Type-1 only)", "noop_5p_type1_only"),
                       ("3' trim NO-OP", "noop_3p")]:
        v = ph[key]
        a(f"  {label:<32s} {v:>10d}  ({100 * v / max(1, out):.1f}%)")
    d = qc["dropped"]
    a("")
    a("Dropped / masked:")
    a(f"  rDNA-masked reads:               {d['rdna_masked_reads']:>10d}")
    a(f"  polyA-pileup buckets:            {d['polya_pileup_buckets']:>10d}"
      f"  ({d['reads_in_dropped_buckets']} reads)")
    a(f"  region-boundary reads:           {d['region_boundary_reads']:>10d}")
    if qc.get("flags"):
        a("")
        a("🔴 QC FLAGS:")
        for f in qc["flags"]:
            a(f"  - {f}")
    return "\n".join(L)


def write_qc_json(qc: Dict, out_dir: Path) -> Path:
    """Persist the QC record next to the Stage-1 FASTQ.

    Machine-readable on purpose: the cross-sample pass (PCA / correlation heatmap over a
    cohort) consumes these sidecars, and console text that scrolls away cannot be
    re-analysed later.
    """
    p = Path(out_dir) / "stage1_qc.json"
    p.write_text(json.dumps(qc, indent=2, sort_keys=True))
    return p
