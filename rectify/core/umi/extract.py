"""Standard (fixed-length, random N-mer) UMI extraction for short-read protocols.

This is the SHORT-READ counterpart to the ONT PCB114 UMI path in
``core/cdna/read_info.py``, and it is deliberately the INVERSE of it:

  * cdna (ONT): the UMI is a STRUCTURED motif (``TT-VVVV-TTT`` + ``GGG`` bridge,
    27 nt) recovered POST-alignment by an exact ``str.find`` of a 23-nt SSP anchor
    in the aligned BAM SEQ. It needs an anchor because its position is unknown.
  * short read (here): the UMI is a fixed-length RANDOM N-mer at a KNOWN offset,
    sliced PRE-alignment straight off the FASTQ record. There is no anchor to
    find, and none is needed.

None of the cdna extraction machinery ports; the clustering (``core/umi/clustering``)
is shared.

Chemistry reference -- Lexogen CORALL Total RNA-Seq V2 (User Guide 171UG394V0111,
App. D p.34 / App. E p.36; V1 UG 095UG190V0110 is word-for-word identical, i.e. V2
did NOT change the UMI):

    5'...ACACTCTTTCCCTACACGACGCTCTTCCGATCT- N[12] -(Insert...

  * UMI = **N12**: 12 nt, fully random, at the START of Read 1.
  * **No linker/spacer** between UMI and insert -- the N12 abuts the insert
    directly, so a positional slice is exact and no anchoring regex is required.
  * **Read 2 carries no UMI** (single-UMI design, not duplex).
  * NOTE: the "12 nt" in Lexogen's *UDI 12 nt* refers to the INDEX length -- an
    independent 12-mer. Do not conflate it with the UMI.

MODES (``--umi-location``):
  ``read-id``  -- the UMI is already in the read name (``umi_tools extract``
                  convention: ``@READID_UMI``). Required for datasets that arrive
                  pre-extracted, which includes Angel's RES libraries (their raw
                  FASTQs no longer exist).
  ``r1-start`` -- slice the first N bases off R1 (raw CORALL).
  ``r2-start`` -- slice the first N bases off R2 (for chemistries that put it there).
"""
from __future__ import annotations

from typing import NamedTuple, Optional, Tuple

# Valid UMI alphabet. N is permitted (an uncalled base) but such UMIs are
# suspect: an N matches nothing exactly and inflates apparent diversity.
_UMI_ALPHABET = frozenset("ACGTN")

UMI_LOCATIONS = ("read-id", "r1-start", "r2-start")

#: Lexogen CORALL V2 (and V1) UMI length. See module docstring.
CORALL_UMI_LENGTH = 12

#: ``umi_tools extract`` default separator between read name and UMI.
DEFAULT_UMI_SEPARATOR = "_"


class ExtractedUmi(NamedTuple):
    """Result of pulling a UMI off a FASTQ record.

    ``seq``/``qual`` are the record MINUS the UMI (for ``read-id`` mode they are
    returned unchanged -- the UMI was never in the sequence).
    """
    umi: str
    umi_qual: Optional[str]
    seq: str
    qual: str


def is_valid_umi(umi: str, expected_length: Optional[int] = None,
                 allow_n: bool = True) -> bool:
    """True if ``umi`` looks like a real UMI.

    Guards against the classic silent failure: a mis-parsed read ID yields a
    plausible-looking string that then clusters into garbage families. Cheap to
    check, and the alternative is a wrong answer that looks fine.
    """
    if not umi:
        return False
    if expected_length is not None and len(umi) != expected_length:
        return False
    alphabet = _UMI_ALPHABET if allow_n else (_UMI_ALPHABET - {"N"})
    return all(c in alphabet for c in umi)


def read_name_of(read_id: str) -> str:
    """The read NAME: everything before the first whitespace.

    A FASTQ header is ``@<name> <comment>``; aligners key on the name only, and
    Illumina's comment field (``1:N:0:INDEX``) must never reach the UMI parser.
    """
    return read_id.split(None, 1)[0].lstrip("@")


def extract_umi_from_read_id(
    read_id: str,
    expected_length: Optional[int] = None,
    separator: str = DEFAULT_UMI_SEPARATOR,
) -> Optional[str]:
    """Parse a UMI appended to the read name by ``umi_tools extract``.

    ``umi_tools`` writes ``@<name><sep><UMI> <comment>`` -- the UMI is appended to
    the NAME, so the comment must be stripped first and the split must be on the
    LAST separator (Illumina names contain ':' but instrument/run names can and do
    contain '_', e.g. ``@RUN_1:...``; splitting on the first '_' would silently
    return a fragment of the run name).

    Returns None if no plausible UMI is present -- callers decide whether that is
    fatal. Never raises on a malformed header.

    >>> extract_umi_from_read_id("@A00197:591:H:3:1101:23475:1031_ACGTACGTACGT 1:N:0:AT", 12)
    'ACGTACGTACGT'
    """
    if not read_id:
        return None
    name = read_name_of(read_id)
    if separator not in name:
        return None
    candidate = name.rsplit(separator, 1)[1]
    if not is_valid_umi(candidate, expected_length=expected_length):
        return None
    return candidate


def strip_umi_from_read_id(read_id: str, separator: str = DEFAULT_UMI_SEPARATOR) -> str:
    """The read name with a trailing ``<sep><UMI>`` removed.

    Needed because mates must share an identical QNAME for pair-aware grouping,
    and because RECTIFY keys its QNAME->tag injection map on the bare name.
    """
    name = read_name_of(read_id)
    if separator not in name:
        return name
    head, tail = name.rsplit(separator, 1)
    return head if is_valid_umi(tail) else name


def extract_umi_from_sequence(
    seq: str,
    qual: str,
    length: int,
) -> Optional[ExtractedUmi]:
    """Slice a fixed-length UMI off the 5' end of a read (the CORALL case).

    Returns None when the read is too short to carry the UMI -- such a record has
    no usable insert anyway and must not be silently emitted with a truncated UMI.

    The slice is exact: CORALL has no linker between the N12 and the insert, so
    bases ``[0:length]`` are UMI and ``[length:]`` are biological sequence.
    """
    if length <= 0:
        raise ValueError(f"UMI length must be positive, got {length}")
    if len(seq) < length:
        return None
    umi = seq[:length]
    return ExtractedUmi(
        umi=umi,
        umi_qual=qual[:length] if qual else None,
        seq=seq[length:],
        qual=qual[length:] if qual else "",
    )


def extract_umi(
    location: str,
    read_id: str,
    seq: str,
    qual: str,
    length: Optional[int] = None,
    separator: str = DEFAULT_UMI_SEPARATOR,
    is_read2: bool = False,
) -> Optional[ExtractedUmi]:
    """Dispatch UMI extraction by ``--umi-location``.

    For ``r1-start``/``r2-start`` the UMI is sliced only from the mate that carries
    it; the other mate is returned untouched but still receives the UMI string, as
    the UMI identifies the FRAGMENT, not the read. That is why both mates of a pair
    end up with the same ``RX``.
    """
    if location not in UMI_LOCATIONS:
        raise ValueError(f"unknown umi location {location!r}; expected one of {UMI_LOCATIONS}")

    if location == "read-id":
        umi = extract_umi_from_read_id(read_id, expected_length=length, separator=separator)
        if umi is None:
            return None
        # The UMI was never in the sequence -- nothing to strip.
        return ExtractedUmi(umi=umi, umi_qual=None, seq=seq, qual=qual)

    if length is None:
        raise ValueError(f"--umi-length is required for --umi-location={location}")

    carries_umi = (location == "r2-start") == is_read2
    if not carries_umi:
        # This mate does not carry the UMI; the caller supplies it from the mate.
        return None
    return extract_umi_from_sequence(seq, qual, length)
