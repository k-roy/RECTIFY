"""Per-clip resolver verdict tag — the seam between Station A and Module 2F.

Station A (``align.overhang_resolver``) and Module 2F (``splice_aware_5prime``)
ask the same question of a terminal soft clip — "does this overhang belong
across a junction?" — with different evidence: the resolver searches a
genome-wide motif index on both ends with no annotation; 2F searches the
annotation + junction pool near the RNA 5' end. Until 2026-09-05 the resolver's
answer died with the resolver: a clip it had refused as uninformative was
assessed again in ``correct``, and a clip it had rejected on edit distance or
ambiguity could be rescued by 2F with nobody counting the disagreement.

The resolver now records what it decided about EVERY assessed terminal clip on
the record it leaves behind — refusals included — in the ``XW`` tag. ``XJ``
still marks placements; ``XW`` is the refusal side of the same ledger. 2F reads
it at its clip-assessment entry and:

* SKIPS its sequence search on ``low_info`` — that verdict depends only on the
  clip sequence. Both sides assess the same junction-proximal slice
  (``ResolverConfig.max_clip_match`` = ``_CLIP_ASSESS_BP`` = 200 bp) with the
  same information bound (``overhang_informativeness``), so recomputing it
  would give the same answer. S2's 10-nt informative-clip floor (ISSUE-006)
  was derived from that same bound.
* still TRIES annotation-guided placement after ``rejected_edit`` /
  ``ambiguous`` / ``no_candidates`` / ``blowup`` / ``repeat`` — annotation and
  the pool are evidence the resolver deliberately lacks — and COUNTS a rescue
  there as ``rescued_over_refusal`` (``ProcessingStats``
  ``ends_2f_rescued_over_resolver_refusal``).
* honors a verdict only while the clip is still the length it was formed on:
  the tag carries ``clip_len`` and a clip that has since been rewritten (a
  resolved placement leaves a shorter remainder, a reanchor grows it) is
  assessed afresh.

Tag format — one entry per assessed genomic side, joined by ``;``::

    XW:Z:<side>:<token>:<bits>:<window>:<clip_len>[;<side>:...]

``side`` L/R (genomic clip side, as XJ); ``token`` one of
:data:`VERDICT_TOKENS`; ``bits`` = effective information bits of the assessed
slice with one decimal, ``-`` when refused before assessment (``repeat``);
``window`` = the information-bounded W_max in bp (0 = refused);
``clip_len`` = the full soft-clip length at verdict time.

The tag is advisory: a malformed entry is ignored, never an error.
"""

from typing import Dict, Optional

import pysam

VERDICT_TAG = 'XW'

#: Resolver outcomes for one clip, in the order resolve_clip can reach them.
VERDICT_TOKENS = (
    'repeat',          # is_repeat_expansion() refused it before assessment
    'low_info',        # assess_overhang(): W_max < 1 bp — cannot localize in any window
    'blowup',          # every enabled class overflowed its candidate budget
    'no_candidates',   # nothing within max_edit_frac in the window
    'rejected_edit',   # best candidate over max_edit_frac
    'ambiguous',       # two different junctions within min_margin
    'resolved',        # placed (XJ carries the placement)
)
VERDICT_REFUSALS = tuple(t for t in VERDICT_TOKENS if t != 'resolved')

#: The one verdict 2F may reuse without re-deriving it (sequence-only).
SEQUENCE_ONLY_VERDICTS = ('low_info',)


def format_verdict_entry(side: str, token: str, bits: Optional[float],
                         window: int, clip_len: int) -> str:
    b = '-' if bits is None else f'{float(bits):.1f}'
    return f'{side}:{token}:{b}:{int(window)}:{int(clip_len)}'


def parse_verdict_tag(value: Optional[str]) -> Dict[str, Dict[str, object]]:
    """``{side: {'token', 'bits', 'window', 'clip_len'}}`` from an XW value.

    Malformed entries are skipped, and a later entry for the same side wins.
    """
    out: Dict[str, Dict[str, object]] = {}
    for entry in (value or '').split(';'):
        parts = entry.split(':')
        if len(parts) != 5 or parts[0] not in ('L', 'R') or parts[1] not in VERDICT_TOKENS:
            continue
        try:
            out[parts[0]] = {
                'token': parts[1],
                'bits': None if parts[2] == '-' else float(parts[2]),
                'window': int(parts[3]),
                'clip_len': int(parts[4]),
            }
        except ValueError:
            continue
    return out


def verdict_for_clip(read: pysam.AlignedSegment, side: str,
                     clip_len: int) -> Optional[Dict[str, object]]:
    """The resolver's verdict for the clip on genomic ``side`` of ``read``,
    or None when there is no tag, no entry for that side, or the clip is no
    longer the length the verdict was formed on."""
    if clip_len <= 0:
        return None
    try:
        raw = read.get_tag(VERDICT_TAG)
    except KeyError:
        return None
    entry = parse_verdict_tag(raw if isinstance(raw, str) else None).get(side)
    if entry is None or entry['clip_len'] != int(clip_len):
        return None
    return entry
