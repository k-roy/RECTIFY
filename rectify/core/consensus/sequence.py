"""Resolve placement-relative SEQ before consensus reads nucleotide evidence."""

import copy
from typing import Dict

import pysam


def decoded_alignment_copy(
    read: pysam.AlignedSegment,
    genome: Dict[str, str],
) -> pysam.AlignedSegment:
    """Return explicit SEQ at the original placement without changing ``read``.

    Literal and SEQ-absent records take the no-copy fast path. SAM ``=`` bases
    have meaning only in aligned M/=/X operations with an available reference.
    Refuse unresolved input instead of treating it as mismatch evidence or
    transplanting its reference-relative spelling to another placement. The
    CIGAR walk skips N spans in constant time; it never materializes intronic
    reference positions. QUAL and all auxiliary tags survive the private copy.
    """
    sequence = read.query_sequence
    if not sequence or '=' not in sequence:
        return read

    def refuse(reason):
        raise ValueError(
            f"Cannot decode consensus SEQ '=' for read {read.query_name!r} "
            f"at {read.reference_name}:{read.reference_start}: {reason}"
        )

    reference = genome.get(read.reference_name) if genome else None
    if reference is None:
        refuse('matching reference sequence is required')
    if read.is_unmapped or not read.cigartuples or read.reference_start < 0:
        refuse('a mapped CIGAR and nonnegative reference start are required')

    chars = list(sequence)
    query_pos, ref_pos = 0, read.reference_start
    for op, length in read.cigartuples:
        if length <= 0:
            refuse('CIGAR operations must have positive length')
        if op in (0, 7, 8):
            if query_pos + length > len(chars) or ref_pos + length > len(reference):
                refuse('CIGAR exceeds stored query or reference sequence')
            for offset in range(length):
                if chars[query_pos + offset] == '=':
                    base = reference[ref_pos + offset].upper()
                    if base == '=':
                        refuse('reference contains an unresolved equals base')
                    chars[query_pos + offset] = base
            query_pos += length
            ref_pos += length
        elif op in (1, 4):
            if query_pos + length > len(chars):
                refuse('CIGAR exceeds stored query sequence')
            if '=' in chars[query_pos:query_pos + length]:
                refuse('equals in insertion or soft clip has no reference base')
            query_pos += length
        elif op in (2, 3):
            ref_pos += length
            if ref_pos > len(reference):
                refuse('CIGAR exceeds reference sequence')
        elif op not in (5, 6):
            refuse('unsupported CIGAR operation')
    if query_pos != len(chars):
        refuse('CIGAR does not consume the stored query sequence')

    decoded = copy.deepcopy(read)
    quality = decoded.query_qualities
    decoded.query_sequence = ''.join(chars)
    decoded.query_qualities = quality
    return decoded
