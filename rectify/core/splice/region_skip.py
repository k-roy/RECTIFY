"""Reference-region skip list for splice-junction rescue machinery.

Motivation (planning/644b, 2026-08-10): on the v5.1 full-genome resolver run,
chrXII alone was 47% of ALL CPU (16,436 s of 34,910 s aggregate; 38.3 reads/s
vs 100-500 elsewhere) — and the cost concentrates in the rDNA repeat
(RDN37-1/RDN37-2 tandem array). rRNA molecules are not spliceosomal
substrates: every junction-rescue attempt there is by construction chasing
alignment noise through two near-identical 9.1 kb repeats, the worst possible
substrate for candidate enumeration (deep coverage x repeat homology). The
same pathology hits ``correct``'s 5' rescue loop (planning/649: pool
composition drives the per-read DP bill).

A read whose alignment overlaps a skip region bypasses junction rescue
entirely (resolver: clip resolution + re-arbitration; correct: 3'SS-truncation
rescue). The read itself is untouched and written through — this is a compute
filter, not a data filter.

Region spec format (``RECTIFY_SKIP_REGIONS`` env var or explicit config):
``chrom:start-end[,chrom:start-end...]`` with 0-based half-open reference
coordinates, or the shorthand ``yeast-rdna`` for the S. cerevisiae R64 rDNA
envelope (chrXII 451,574-468,812: RDN37-1 through NTS2-2, from the lab
R64-5-1 GFF). Shorthands and explicit specs may be mixed.
"""

import os
from typing import Dict, List, Tuple

# S. cerevisiae R64 rDNA repeat envelope on chrXII (0-based half-open).
# GFF (1-based): RDN37-1 starts 451,575; NTS2-2 ends 468,812.
YEAST_RDNA_SPEC = 'chrXII:451574-468812'

_SHORTHANDS = {
    'yeast-rdna': YEAST_RDNA_SPEC,
}

SkipRegions = Dict[str, List[Tuple[int, int]]]

ENV_VAR = 'RECTIFY_SKIP_REGIONS'


def parse_skip_regions(spec: str) -> SkipRegions:
    """Parse a region spec string into {chrom: [(start, end), ...]}.

    Raises ValueError on malformed entries — a silently-ignored typo here
    would quietly disable the filter (or worse, skip the wrong locus).
    """
    regions: SkipRegions = {}
    for raw in spec.split(','):
        entry = raw.strip()
        if not entry:
            continue
        entry = _SHORTHANDS.get(entry.lower(), entry)
        try:
            chrom, span = entry.rsplit(':', 1)
            start_s, end_s = span.split('-')
            start, end = int(start_s), int(end_s)
        except ValueError:
            raise ValueError(
                f'malformed skip-region entry {entry!r} '
                f'(expected chrom:start-end or a shorthand '
                f'{sorted(_SHORTHANDS)})') from None
        if not chrom or start < 0 or end <= start:
            raise ValueError(f'invalid skip-region interval {entry!r}')
        regions.setdefault(chrom, []).append((start, end))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def skip_regions_from_env() -> SkipRegions:
    """Skip regions from the RECTIFY_SKIP_REGIONS env var (empty if unset)."""
    spec = os.environ.get(ENV_VAR, '')
    return parse_skip_regions(spec) if spec.strip() else {}


def overlaps_skip_region(regions: SkipRegions, chrom: str,
                         start: int, end: int) -> bool:
    """True when [start, end) overlaps any skip interval on chrom.

    Interval lists are short (typically one), so a linear scan is fine.
    """
    if not regions:
        return False
    for r_start, r_end in regions.get(chrom, ()):
        if start < r_end and end > r_start:
            return True
    return False
