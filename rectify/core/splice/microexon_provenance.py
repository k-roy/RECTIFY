"""Placement-specific provenance for successfully drawn micro-exons.

Xb:Z stores compact JSON with version, chrom, strand and calls. Each call names
the original intron, the exons actually drawn, and that call's equally good
alternatives. Coordinates are genomic, zero-based and half-open on both strands.
XB belongs to cDNA strand-split counts; XV belongs to validation labels. Neither
legacy tag is written here or treated as proof of a successful micro-exon draw.
"""

import json
import re

TAG_MICROEXON_PROVENANCE = 'Xb'
MICROEXON_PROVENANCE_VERSION = 1

_LEGACY_SEGMENTS = re.compile(r'[^,|]+:[0-9]+-[0-9]+(?:[,|][^,|]+:[0-9]+-[0-9]+)*')


def _span(value):
    return (isinstance(value, list) and len(value) == 2
            and all(type(x) is int for x in value) and 0 <= value[0] < value[1])


def _valid_call(call):
    if not isinstance(call, dict) or not _span(call.get('intron')):
        return False
    exons = call.get('exons')
    if not isinstance(exons, list) or not exons:
        return False
    cursor, end = call['intron']
    for exon in exons:
        if not _span(exon) or not cursor < exon[0] < exon[1] < end:
            return False
        cursor = exon[1]
    if 'junctions' in call:
        # Chimeric selection can retain only part of a previously drawn call.
        selected = call['junctions']
        if (not isinstance(selected, list) or not selected
                or not all(_span(j) for j in selected)
                or not {tuple(j) for j in selected} <= _call_junctions(call, full=True)):
            return False
    return isinstance(call.get('alternatives'), str)


def _call_junctions(call, full=False):
    if not full and 'junctions' in call:
        return {tuple(j) for j in call['junctions']}
    cursor, end = call['intron']
    junctions = set()
    for start, stop in call['exons']:
        junctions.add((cursor, start))
        cursor = stop
    junctions.add((cursor, end))
    return junctions


def load_microexon_calls(read):
    """Validated calls for this placement; [] when absent, None when untrusted.

    Version/shape/contig/strand mismatches are explicit unknown provenance. A
    caller must still intersect the recorded introns with the current CIGAR,
    because later clipping or junction refinement can remove a recorded edge.
    """
    if not read.has_tag(TAG_MICROEXON_PROVENANCE):
        return []
    raw = read.get_tag(TAG_MICROEXON_PROVENANCE)
    if not isinstance(raw, str):
        return None
    try:
        payload = json.loads(raw)
    except (ValueError, RecursionError):
        return None
    if (not isinstance(payload, dict)
            or type(payload.get('v')) is not int
            or payload['v'] != MICROEXON_PROVENANCE_VERSION
            or payload.get('chrom') != read.reference_name
            or payload.get('strand') != ('-' if read.is_reverse else '+')):
        return None
    calls = payload.get('calls')
    if not isinstance(calls, list) or not all(_valid_call(call) for call in calls):
        return None
    return calls


def record_microexon_draws(read, calls):
    """Append successful calls without changing XB/XV or duplicating a replay."""
    if not calls:
        return
    previous = load_microexon_calls(read) or []
    for call in calls:
        if not _valid_call(call):
            raise ValueError('Invalid successful micro-exon call')
        if call not in previous:
            previous.append(call)
    payload = {'v': MICROEXON_PROVENANCE_VERSION, 'chrom': read.reference_name,
               'strand': '-' if read.is_reverse else '+', 'calls': previous}
    read.set_tag(TAG_MICROEXON_PROVENANCE,
                 json.dumps(payload, separators=(',', ':'), sort_keys=True), value_type='Z')


def microexon_junction_provenance(read):
    """Return (recorded junction -> alternatives, has_unverified_provenance).

    Legacy XB coordinate strings cannot distinguish a planned call from a
    successful rewrite. Flag them for recount/reprocessing, never credit them
    as verified draws. cDNA n_top/n_bottom strings do not match that grammar.
    The mapping is in genomic coordinates; the caller checks its live N-ops.
    """
    calls = load_microexon_calls(read)
    unverified = calls is None
    if read.has_tag('XB'):
        legacy = read.get_tag('XB')
        unverified |= isinstance(legacy, str) and _LEGACY_SEGMENTS.fullmatch(legacy) is not None
    junctions = {}
    for call in calls or []:
        alternatives = {alt for alt in call['alternatives'].split(';') if alt}
        for junction in _call_junctions(call):
            junctions.setdefault(junction, set()).update(alternatives)
    return junctions, unverified


def copy_selected_microexon_provenance(out, result, aligner_reads, anchor):
    """Carry provenance from the source of each selected N, never the SEQ donor.

    The optional per-call junctions list restricts credit when only part of a
    micro-exon configuration survives segment selection. Newly constructed
    bridge introns without a matching source N receive no micro-exon credit.
    """
    # Most reads have no micro-exon call. Avoid walking their segment events or
    # output CIGAR just to rediscover the absent tag.
    if not any(r.has_tag(TAG_MICROEXON_PROVENANCE)
               for r in (aligner_reads or {'anchor': anchor}).values()):
        return
    selected = {}
    if result.all_segment_scores:
        for segment in result.all_segment_scores:
            selected.setdefault(segment.winning_aligner, set()).update(
                (event.r_start, event.r_end)
                for event in segment.cigar_events.get(segment.winning_aligner, ())
                if event.op == 3)
        sources = aligner_reads or {}
    elif not result.is_chimeric:
        # Fallback and single-arm results use the entire winning placement.
        sources = {result.anchor_aligner: anchor}
        selected[result.anchor_aligner] = None
    else:
        return  # No source-event attribution available; do not guess.

    live = set()
    pos = out.reference_start
    for op, length in out.cigartuples or ():
        if op == 3:
            live.add((pos, pos + length))
        if op in (0, 2, 3, 7, 8):
            pos += length

    retained = []
    for name, junctions in selected.items():
        source = sources.get(name)
        if (source is None or source.reference_name != out.reference_name
                or source.is_reverse != out.is_reverse):
            continue
        allowed = live if junctions is None else live & junctions
        for call in load_microexon_calls(source) or []:
            surviving = _call_junctions(call) & allowed
            if surviving:
                retained.append(dict(call, junctions=[list(j) for j in sorted(surviving)]))
    record_microexon_draws(out, retained)
