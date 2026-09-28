"""Keep alignment-derived tags consistent with the final explicit placement."""

import logging

from ...utils.genome import get_chrom_sequence

logger = logging.getLogger(__name__)
_STALE_TAGS = ('NM', 'MD', 'AS', 'ms', 'cs', 'de', 'dv', 'UQ')
_NT16 = {base: code for code, base in enumerate('=ACMGRSVTWYHKDBN')}
_NT16.update({'0': 1, '1': 2, '2': 4, '3': 8})


def placement_state(read):
    """Immutable state relevant to reference-dependent and quality-based tags."""
    return (read.reference_id, read.reference_start, read.flag,
            tuple(read.cigartuples or ()), read.query_sequence,
            bytes(read.query_qualities) if read.query_qualities is not None else None)


def decode_original_sequence(read, genome):
    """Decode SEQ '=' before surgery, or refuse without changing the record.

    Compressed bases only have meaning in aligned M/=/X cells of the original
    placement. H and N never allocate query/reference-sized pair arrays.
    """
    sequence = read.query_sequence
    if not sequence or '=' not in sequence:
        return False

    def refuse(reason):
        raise ValueError(
            f"Cannot decode original SEQ '=' for read {read.query_name!r}: {reason}"
        )

    reference, _ = get_chrom_sequence(genome, read.reference_name)
    if reference is None:
        refuse('matching reference sequence is required before correction')
    if read.is_unmapped or read.reference_start < 0 or not read.cigartuples:
        refuse('mapped CIGAR with a nonnegative start is required')
    chars = list(sequence)
    query_pos, ref_pos = 0, read.reference_start
    for op, length in read.cigartuples:
        if length <= 0:
            refuse('CIGAR operations must have positive length')
        if op in (0, 7, 8):
            if query_pos + length > len(chars) or ref_pos + length > len(reference):
                refuse('CIGAR exceeds stored query or available reference')
            for offset in range(length):
                if chars[query_pos + offset] == '=':
                    base = reference[ref_pos + offset].upper()
                    if base not in _NT16 or base == '=':
                        refuse('reference has no encodable base at a compressed cell')
                    chars[query_pos + offset] = base
            query_pos += length
            ref_pos += length
        elif op in (1, 4):
            if query_pos + length > len(chars):
                refuse('CIGAR exceeds stored query')
            if '=' in chars[query_pos:query_pos + length]:
                refuse('equals in insertion or soft clip has no reference placement')
            query_pos += length
        elif op in (2, 3):
            ref_pos += length
            if ref_pos > len(reference):
                refuse('CIGAR exceeds available reference')
        elif op not in (5, 6):
            refuse('unsupported CIGAR operation')
    if query_pos != len(chars):
        refuse('CIGAR does not consume stored query')
    quality = read.query_qualities
    read.query_sequence = ''.join(chars)
    if quality is not None:
        read.query_qualities = quality
    return True


def calculate_nm_md(read, genome):
    """Return calmd-compatible NM/MD from complete final reference and SEQ.

    Match semantics follow samtools bam_fillmd1_core: identical nt16 symbols
    except N match, while N/N is a mismatch. Inspect M, = and X bases alike.
    A skipped N advances only the reference cursor and contributes no edits or
    MD length. Refuse incomplete data instead of returning a partial MD string.
    """
    sequence = read.query_sequence
    reference, _ = get_chrom_sequence(genome, read.reference_name)
    if not sequence or reference is None:
        raise ValueError('final literal SEQ and matching reference are required')
    if '=' in sequence:
        raise ValueError('final SEQ still contains unresolved equals')
    if read.is_unmapped or read.reference_start < 0 or not read.cigartuples:
        raise ValueError('final mapped CIGAR and nonnegative start are required')
    query_pos, ref_pos = 0, read.reference_start
    matched, nm, md = 0, 0, []
    for op, length in read.cigartuples:
        if length <= 0:
            raise ValueError('final CIGAR operations must have positive length')
        if op in (0, 7, 8):
            if query_pos + length > len(sequence) or ref_pos + length > len(reference):
                raise ValueError('final CIGAR exceeds stored query or available reference')
            for offset in range(length):
                query_base = sequence[query_pos + offset].upper()
                ref_base = reference[ref_pos + offset].upper()
                if ref_base == '=' or ref_base == '\0':
                    raise ValueError('final reference contains an unresolved base')
                query_code, ref_code = _NT16.get(query_base, 15), _NT16.get(ref_base, 15)
                if query_code == ref_code and query_code != 15:
                    matched += 1
                else:
                    md.extend((str(matched), ref_base))
                    matched = 0
                    nm += 1
            query_pos += length
            ref_pos += length
        elif op == 2:
            if ref_pos + length > len(reference):
                raise ValueError('final deletion exceeds available reference')
            deleted = reference[ref_pos:ref_pos + length].upper()
            if '=' in deleted or '\0' in deleted:
                raise ValueError('final deletion reference contains an unresolved base')
            md.extend((str(matched), '^', deleted))
            matched = 0
            nm += length
            ref_pos += length
        elif op in (1, 4):
            query_pos += length
            if query_pos > len(sequence):
                raise ValueError('final CIGAR exceeds stored query')
            if op == 1:
                nm += length
        elif op == 3:
            ref_pos += length
            if ref_pos > len(reference):
                raise ValueError('final skipped region exceeds available reference')
        elif op not in (5, 6):
            raise ValueError('unsupported final CIGAR operation')
    if query_pos != len(sequence):
        raise ValueError('final CIGAR does not consume stored query')
    md.append(str(matched))
    return nm, ''.join(md)


def finalize_alignment_tags(read, before, genome):
    """Refresh final tags after an actual edit, preserving unchanged records.

    With missing/truncated reference or missing SEQ, clear obsolete fields and
    report why NM/MD are unavailable. Other metadata and typed tags are untouched.
    """
    if placement_state(read) == before:
        return False
    try:
        nm_md = calculate_nm_md(read, genome)
    except ValueError as exc:
        nm_md = None
        logger.warning('Edited read %s has unavailable NM/MD: %s', read.query_name, exc)
    for tag in _STALE_TAGS:
        while read.has_tag(tag):
            read.set_tag(tag, None)
    if nm_md is not None:
        read.set_tag('NM', nm_md[0], value_type='i')
        read.set_tag('MD', nm_md[1], value_type='Z')
    return True
