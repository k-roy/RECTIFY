"""Mapped-read CPA concordance — the headline generalization metric.

What fraction of nascent genome-aware poly-A 3' ends (``oaNT>=2``) from a Tier-1
pileup land on known DRS CPA sites (+- ``reg``), versus the genome null (the
fraction of the 2-strand genome covered by CPA windows)? The enrichment over
null (the orthogonal DRS map being the ground truth) is the per-dataset
validation that these are real cleaved + polyadenylated intermediates.

Default ``nuclear_bp`` is the S. cerevisiae R64 nuclear genome (excl. mito);
override for other organisms.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Dict

from . import SCER_NUCLEAR_BP

REG_DEFAULT = 5


def load_cpa_set(drs_clusters: str | Path, *, reg: int = REG_DEFAULT) -> set:
    """Load DRS CPA modal positions into a ``{(chrom, strand, pos)}`` set.

    The single source of truth for "is this coordinate a known CPA site" — used
    both by :func:`mapped_cpa_concordance` (genome-null enrichment) and by the
    per-read ``at_cpa`` flag in the reads.parquet sidecar, so the two agree
    exactly (each modal position is expanded to a +-``reg`` window).
    """
    cpa: set = set()
    with open(drs_clusters) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c, s, m = r["chrom"], r["strand"], int(float(r["modal_position"]))
            for g in range(m - reg, m + reg + 1):
                cpa.add((c, s, g))
    return cpa


def mapped_cpa_concordance(
    pileup_path: str | Path,
    drs_clusters: str | Path,
    *,
    label: str = "",
    reg: int = REG_DEFAULT,
    nuclear_bp: int = SCER_NUCLEAR_BP,
) -> Dict[str, object]:
    """Concordance of a Tier-1 pileup's ``oaNT>=2`` reads with the DRS CPA map.

    Returns ``{label, mapped_reads, oaNT_reads, at_cpa, frac_at_cpa, null,
    enrichment, mean_polya_len}``.
    """
    cpa = load_cpa_set(drs_clusters, reg=reg)

    tot_reads = tot_oa = at_oa = sum_len = n_len = 0
    with gzip.open(pileup_path, "rt") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for row in rd:
            nr = int(row[3]); oa1 = int(row[9]); oa2 = int(row[10]); soa = int(row[12])
            tot_reads += nr
            tot_oa += oa2
            sum_len += soa
            n_len += oa1
            if (row[0], row[2], int(row[1])) in cpa:
                at_oa += oa2

    null = len(cpa) / (nuclear_bp * 2)
    frac = at_oa / max(tot_oa, 1)
    return {
        "label": label,
        "mapped_reads": tot_reads,
        "oaNT_reads": tot_oa,
        "at_cpa": at_oa,
        "frac_at_cpa": frac,
        "null": null,
        "enrichment": frac / max(null, 1e-9),
        "mean_polya_len": sum_len / max(n_len, 1),
    }
