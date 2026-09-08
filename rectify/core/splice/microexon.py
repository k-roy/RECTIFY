"""STATION B, first pass — recover ORPHANED MICROEXONS from an insertion beside a junction (ISSUE-040).

A splice aligner will not split an intron for a 6- or 9-nt internal exon: chaining needs a seed, and a
15-mer minimizer cannot live inside a 6-nt exon. Tested directly on card 3's read (minimap2 2.28,
`-x splice:hq -u f`): supplying that transcript's OWN junctions with `--junc-bed` at `--junc-bonus` 9 and
at 30 gives a BYTE-IDENTICAL CIGAR. The aligner does the reasonable thing and leaves the microexon's
bases as an insertion beside the intron it did find. This is a structural limit, not a tuning failure,
and no splice-aligner parameter reaches it.

Measured in one library (GSB_191, i020n_5d70a92 T1): of 153 insertions of >= 3 nt adjacent to an N-op,
**32 (21 %) exactly contain an annotated exon <= 30 nt from inside that same intron**, and 59 % are a
multiple of 3 — so losing them is usually a silent in-frame deletion rather than a frameshift, invisible
to an ORF check. GENCODE basic holds only 5,870 exons <= 30 nt, so the candidate set is enumerable.

**What carries the evidence is the CONFIGURATION, not the motif.** Measured inside card 3's own 10,182-bp
intron: 52 positions could host a 6-nt segment with canonical flanks and 48 could host a 9-nt one — about
one every 200 bp, so canonical flanks alone are nearly free (and admitting more motifs LOWERS the ante).
What makes the answer unique is requiring an EXACT match to the read's own bases in ORDER, consuming the
insertion exactly: that leaves exactly one home for each segment. A 6-mer alone would still be borderline
de novo (52/4096 ~ 1.3 % per intron); the 9-mer is 48/262144 ~ 0.02 %. Hence the two rules this module
enforces and `station_b` (a de-novo search) will need to earn:

1. annotated microexons are candidates FIRST — this pass proposes nothing that is not already annotated;
2. the split must consume the inserted bases EXACTLY and in transcript order. Partial consumption is a
   refusal, not a partial credit: leftover bases would have to be glued to an N as an indel, which is the
   shape ISSUE-031/038 ban.

Two microexons in one intron is not exotic — card 3 (a0d80d8a) has exactly that, `CAGCTC` (6 nt) and
`TTGTGCAAA` (9 nt), both in ENST00000434715.7 — so nothing here is built for exactly one.
"""
from __future__ import annotations

import gzip
import hashlib
import logging
import os
import random
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class MicroexonCall:
    """What station B found for one junction-adjacent insertion.

    ``segments`` is the configuration to DRAW; ``alternatives`` are the other configurations that
    explain the same bases — the tied ones the draw did not win (Kevin: keep them as noted
    equally-good alternatives) followed by any lower-scoring ones. ``n_tied`` is how many were tied
    for best, so a reader can tell an arbitrary pick (n_tied > 1) from a unique answer.
    """
    chrom: str
    segments: List[Tuple[int, int]]
    alternatives: List[List[Tuple[int, int]]]
    n_tied: int
    bits: float
    new_cigar: List[Tuple[int, int]]
    insertion_len: int
    intron: Tuple[int, int]

    @property
    def ambiguous(self) -> bool:
        return self.n_tied > 1

#: An exon this short cannot hold an aligner seed, so it is the class this pass exists for.
MAX_MICROEXON_LEN = 30
#: A recovered flanking intron shorter than this is not a splicing event worth proposing.
MIN_FLANKING_INTRON = 30
#: An insertion shorter than this is ordinary alignment noise, not an orphaned exon.
MIN_INSERTION_LEN = 3
#: Bits added when every segment of a split shares one annotated transcript (see split_score).
TX_COHERENCE_BITS = 4.0
#: Two splits within this many bits of each other are TIED, and the draw picks one at random.
TIE_EPSILON = 1e-9
#: Stop enumerating configurations past this many: beyond a handful, the split is not evidence.
MAX_SPLITS = 16

_M, _I, _D, _N, _S, _H, _P, _EQ, _X = range(9)

# Intron ends are judged as a PAIR and in TRANSCRIPT orientation, never as two independent genomic
# dinucleotides: checking sets separately would admit a GT donor with an AC acceptor, which is one U2
# end and one U12 end, not a splice class.
#
# 🔴 THE 3'SS SIDE IS SPECIES-DEPENDENT AND THIS MODULE USED TO PRETEND IT WAS NOT (Kevin,
# 2026-09-08: "both human and yeast can utilize GT-AG/GC-AG/AT-AC, but yeast also have the
# non-canonical 3' SS arms BG and AT"). The hierarchy he names is already in the codebase —
# `junction_scoring._3ss_tier_from_rna_trinucleotide`, derived from yeast splicing observations:
#
#   tier 0  YAG  (C/T)AG   most common, highest efficiency
#   tier 1  RAG  (A/G)AG
#   tier 2  NBG  B = C/G/T — "non-canonical but observed in yeast"   <- Kevin's "BG"
#   tier 3  NAT                  very rare non-canonical             <- Kevin's "AT"
#   tier 4  other
#
# So the acceptor is judged by TIER against a per-organism ceiling rather than by a hard-coded
# dinucleotide set, and there is one model in the tree instead of two that can drift apart.
#
# The DONOR side is species-neutral: GT and GC for U2, AT for U12. AT-AC is a real class in every
# organism rectify runs on — yeast splices it through its major spliceosome (Talkish 2019, Kevin's
# call 2026-09-05) — so it is admitted everywhere, but only PAIRED with an AC acceptor.
_U2_DONORS = ('GT', 'GC')
_U12_DONOR = 'AT'
_U12_ACCEPTOR = 'AC'
#: 3'SS tier ceiling per organism. Human keeps the U2 hierarchy at RAG; yeast reaches NBG/NAT.
MICROEXON_MAX_3SS_TIER = {
    'saccharomyces_cerevisiae': 3,
}
MICROEXON_MAX_3SS_TIER_DEFAULT = 1
_SPECIES_MAX_3SS_TIER = MICROEXON_MAX_3SS_TIER_DEFAULT


def set_species(organism) -> None:
    """Set the 3'SS tier ceiling from the run's organism; unknown organisms keep the default."""
    global _SPECIES_MAX_3SS_TIER
    key = (organism or '').strip().lower().replace(' ', '_')
    _SPECIES_MAX_3SS_TIER = MICROEXON_MAX_3SS_TIER.get(key, MICROEXON_MAX_3SS_TIER_DEFAULT)


def max_3ss_tier() -> int:
    raw = os.environ.get('RECTIFY_MICROEXON_MAX_3SS_TIER', '').strip()
    try:
        return int(raw) if raw else _SPECIES_MAX_3SS_TIER
    except ValueError:
        return _SPECIES_MAX_3SS_TIER


def _rc(x: str) -> str:
    return x[::-1].translate(str.maketrans('ACGTacgtN', 'TGCAtgcaN'))


def intron_ends_rna(genome_seq: str, start: int, end: int, strand: str):
    """``(donor_dinucleotide, acceptor_trinucleotide)`` of ``[start, end)`` in TRANSCRIPT
    orientation — the orientation every motif rule in this file and in junction_scoring is written
    in. Reading a minus-strand intron as if it were plus is the bug ISSUE-038 was written for."""
    if strand == '-':
        return _rc(genome_seq[end - 2:end]).upper(), _rc(genome_seq[start:start + 3]).upper()
    return genome_seq[start:start + 2].upper(), genome_seq[end - 3:end].upper()

#: Installed once per run (and per spawned worker) from the annotation; empty = station B inert.
_MICROEXON_INDEX: Dict[str, List[Tuple[int, int]]] = {}


def set_microexon_index(index) -> None:
    """Install the annotated micro-exon index (``load_microexons``); ``None`` clears it."""
    global _MICROEXON_INDEX
    _MICROEXON_INDEX = dict(index or {})


def microexon_index() -> Dict[str, List[Tuple[int, int]]]:
    return _MICROEXON_INDEX


def station_b_mode() -> str:
    """``'apply'`` (DEFAULT since 2026-09-08 — draw the micro-exons) or ``'report'`` (record what the
    search found, rewrite nothing). Env RECTIFY_STATION_B=report opts out.

    Kevin flipped the default on the integration wave: "Let's have all stations be default ON."
    The evidence for this one specifically: 304 draws over the SMA panel with ZERO of 304 applied
    rows disagreeing between the TSV and the BAM, against a pre-existing ISSUE-024 rate of 0.8 %,
    and it is the one class a splice aligner structurally cannot reach.
    """
    return 'report' if os.environ.get('RECTIFY_STATION_B', '').strip().lower() == 'report' else 'apply'


def _transcript_id(attrs: str) -> str:
    for key in ('transcript_id', 'Parent'):
        i = attrs.find(key)
        if i == -1:
            continue
        rest = attrs[i + len(key):].lstrip(' =')
        if rest.startswith('"'):
            j = rest.find('"', 1)
            return rest[1:j] if j != -1 else rest[1:]
        for sep in (';', '\t', ' '):
            k = rest.find(sep)
            if k != -1:
                rest = rest[:k]
        return rest.strip()
    return ''


def load_microexons(annotation_path: str,
                    max_len: int = MAX_MICROEXON_LEN) -> Dict[str, List[Tuple[int, int]]]:
    """``{chrom: sorted [(start, end)]}`` of annotated exons no longer than *max_len* (0-based,
    half-open). Deduplicated across transcripts — the same short exon appears in many.

    The transcripts each exon belongs to are kept in a side table (``transcripts_of``) rather than in
    the tuples, so the index stays a plain sorted list a bisect can search and a worker can pickle
    cheaply. Transcript coherence is what breaks ties between equally exact splits (ISSUE-035:
    card 3's two micro-exons are both in ENST00000434715.7 — that is the evidence, not the motif)."""
    from ...utils.genome import standardize_chrom_name

    per_chrom: Dict[str, set] = defaultdict(set)
    tx: Dict[Tuple[str, int, int], set] = defaultdict(set)
    _open = gzip.open if str(annotation_path).endswith('.gz') else open
    with _open(annotation_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2].lower() != 'exon':
                continue
            start = int(parts[3]) - 1                 # GTF/GFF are 1-based inclusive
            end = int(parts[4])
            if 0 < end - start <= max_len:
                chrom = standardize_chrom_name(parts[0])
                per_chrom[chrom].add((start, end))
                t = _transcript_id(parts[8])
                if t:
                    tx[(chrom, start, end)].add(t)
    out = {c: sorted(v) for c, v in per_chrom.items()}
    out['__transcripts__'] = {k: frozenset(v) for k, v in tx.items()}   # side table, not a contig
    logger.debug("load_microexons: %d exons <= %d nt over %d contigs",
                 sum(len(v) for k, v in out.items() if k != '__transcripts__'),
                 max_len, len(out) - 1)
    return out


def transcripts_of(index, chrom: str, start: int, end: int) -> frozenset:
    """Transcripts carrying this annotated exon; empty when the index has no side table."""
    return (index.get('__transcripts__') or {}).get((chrom, start, end), frozenset())


def exons_inside(index: Dict[str, List[Tuple[int, int]]], chrom: str,
                 intron_start: int, intron_end: int) -> List[Tuple[int, int]]:
    """Annotated micro-exons lying wholly inside ``[intron_start, intron_end)``, in genomic order."""
    ex = index.get(chrom) or ()
    if not ex or chrom == '__transcripts__':
        return []
    lo = bisect_left(ex, (intron_start, -1))
    out = []
    for s, e in ex[lo:]:
        if s >= intron_end:
            break
        if e <= intron_end:
            out.append((s, e))
    return out


def _intron_pair_ok(genome_seq: str, start: int, end: int, strand: str) -> bool:
    """Whether ``[start, end)`` is a legal intron — judged as a PAIR, in transcript orientation, with
    the 3'SS graded by the shared yeast-derived tier model at this organism's ceiling."""
    from .junction_scoring import _3ss_tier_from_rna_trinucleotide

    if end - start < MIN_FLANKING_INTRON:
        return False
    donor, acc3 = intron_ends_rna(genome_seq, start, end, strand)
    if len(donor) != 2 or len(acc3) != 3:
        return False
    if donor == _U12_DONOR:                       # U12: AT pairs only with AC
        return acc3[1:] == _U12_ACCEPTOR
    if donor not in _U2_DONORS:
        return False
    return _3ss_tier_from_rna_trinucleotide(acc3) <= max_3ss_tier()


def split_introns(split: List[Tuple[int, int]], intron_start: int,
                  intron_end: int) -> List[Tuple[int, int]]:
    """The introns a split creates, in genomic order: the drawn intron partitioned by the segments."""
    out = []
    cursor = intron_start
    for (s_, e_) in split:
        out.append((cursor, s_))
        cursor = e_
    out.append((cursor, intron_end))
    return out


def split_is_canonical(split: List[Tuple[int, int]], intron_start: int, intron_end: int,
                       strand: str, genome_seq: str) -> bool:
    """Every intron the split creates is a legal PAIR — including the two that reuse the original
    donor and acceptor, which is what makes this a test of the whole configuration rather than of
    the segments in isolation."""
    return all(_intron_pair_ok(genome_seq, s_, e_, strand)
               for s_, e_ in split_introns(split, intron_start, intron_end))


def find_microexon_splits(inserted_seq: str,
                          chrom: str,
                          intron_start: int,
                          intron_end: int,
                          strand: str,
                          genome_seq: str,
                          index: Dict[str, List[Tuple[int, int]]],
                          limit: int = MAX_SPLITS,
                          ) -> List[List[Tuple[int, int]]]:
    """EVERY ordered set of annotated micro-exons that consumes *inserted_seq* exactly.

    *inserted_seq* is the insertion's bases in GENOMIC orientation (as they appear in SEQ), and the
    consumption test runs in genomic order for both strands — on the minus strand the transcript order
    is simply the reverse, which the caller never has to redo.

    Every returned split satisfies: exact sequence identity per segment, exact and complete
    consumption of the insertion, strictly increasing genomic order, all segments inside the intron,
    each resulting intron at least ``MIN_FLANKING_INTRON``, and canonical flanks on both sides of
    every segment. The search stops at *limit* splits — beyond a handful the configuration is not
    evidence for anything and the read is better left alone.
    """
    n = len(inserted_seq)
    if n < MIN_INSERTION_LEN:
        return []
    cands = [(s_, e_) for (s_, e_) in exons_inside(index, chrom, intron_start, intron_end)
             if e_ - s_ <= n]
    if not cands:
        return []
    up = inserted_seq.upper()
    found: List[List[Tuple[int, int]]] = []

    def _search(consumed: int, min_start: int, prev_end: int, chosen: List[Tuple[int, int]]) -> None:
        if len(found) >= limit:
            return
        if consumed == n:
            if chosen and split_is_canonical(chosen, intron_start, intron_end, strand, genome_seq):
                found.append(list(chosen))
            return
        for (s_, e_) in cands:
            if s_ < min_start:
                continue
            length = e_ - s_
            if consumed + length > n:
                continue
            if not chosen:
                if s_ - intron_start < MIN_FLANKING_INTRON:
                    continue
            elif s_ - prev_end < MIN_FLANKING_INTRON:
                continue
            if genome_seq[s_:e_].upper() != up[consumed:consumed + length]:
                continue
            # Cheap pre-filter only: the intron this segment opens must at least begin with a
            # dinucleotide that can open one, read in transcript orientation. The full PAIR test
            # runs on the complete split.
            _open = (_rc(genome_seq[s_ - 2:s_]) if strand == '-' else genome_seq[e_:e_ + 2]).upper()
            if _open not in (_U2_DONORS + (_U12_DONOR,)):
                continue
            _search(consumed + length, e_ + MIN_FLANKING_INTRON, e_, chosen + [(s_, e_)])
            if len(found) >= limit:
                return

    _search(0, intron_start, intron_start, [])
    return found


def split_score(split: List[Tuple[int, int]], chrom: str, index) -> float:
    """How much evidence a configuration carries, in bits.

    Two terms, both derived rather than tuned:

    * **sequence** — each segment's own exact match is ``2 * length`` bits of agreement, and the
      chance of finding a canonically-flanked home for a k-mer inside one intron falls as 4^-k
      (measured in card 3's intron: a 6-mer lands once in 79, a 9-mer once in 5,600). This is the
      term that makes a 9-nt segment worth far more than a 6-nt one.
    * **transcript coherence** — ``TX_COHERENCE_BITS`` when every segment of the split is carried by
      one common transcript. Two micro-exons that co-occur in an annotated transcript are a spliced
      structure someone has already observed; two that never co-occur are a coincidence of the
      annotation. Set to 4 bits: it separates coherent from incoherent splits without ever
      outweighing three matched bases.

    A split into FEWER, LONGER segments therefore outscores a split of the same bases into more,
    shorter ones only through the coherence term — the sequence term is identical, since both consume
    the same insertion. That is deliberate: the bases are the same evidence either way.
    """
    seq_bits = 2.0 * sum(e - s_ for s_, e in split)
    tx = None
    for (s_, e) in split:
        t = transcripts_of(index, chrom, s_, e)
        tx = set(t) if tx is None else (tx & set(t))
    coherent = bool(tx)
    return seq_bits + (TX_COHERENCE_BITS if coherent else 0.0)


def choose_split(splits: List[List[Tuple[int, int]]], chrom: str, index, read_name: str):
    """``(chosen, alternatives, n_tied)`` — the configuration to draw, and the equally good ones.

    Kevin, 2026-09-07: *where several candidates are equally plausible, pick one at random for the
    actual alignment and keep the others as noted equally-good alternatives.* "At random" here is
    seeded by the READ NAME, so the choice is arbitrary with respect to the biology (no positional
    or ordering bias creeps in) and yet reproducible: the same BAM re-run gives the same answer, and
    a reviewer can re-derive the pick. Only splits within ``TIE_EPSILON`` bits of the best are tied;
    a lower-scoring split is an alternative that lost, and is reported as such but never drawn.
    """
    if not splits:
        return None, [], 0
    scored = sorted(((split_score(sp, chrom, index), sp) for sp in splits),
                    key=lambda t: (-t[0], t[1]))
    best = scored[0][0]
    tied = [sp for sc, sp in scored if best - sc <= TIE_EPSILON]
    rest = [sp for sc, sp in scored if best - sc > TIE_EPSILON]
    rng = random.Random(hashlib.sha1(read_name.encode('utf-8')).hexdigest())
    chosen = rng.choice(tied)
    # EQUALLY GOOD means equally good. `rest` lost on evidence — a shorter total match, or segments
    # that never co-occur in a transcript — and recording a loser beside a tie would blur exactly the
    # distinction Kevin asked to preserve. `n_tied` says how many were in the draw.
    alternatives = [sp for sp in tied if sp != chosen]
    return chosen, alternatives, len(tied)


def format_segments(chrom: str, split) -> str:
    """``chrom:start-end`` per segment, comma-joined in genomic order (0-based, half-open)."""
    return ','.join(f'{chrom}:{s_}-{e}' for s_, e in (split or []))


def format_alternatives(chrom: str, splits) -> str:
    """Alternative configurations, ';'-separated; each formatted as ``format_segments``."""
    return ';'.join(format_segments(chrom, sp) for sp in (splits or []))


def find_microexon_split(inserted_seq: str,
                         chrom: str,
                         intron_start: int,
                         intron_end: int,
                         strand: str,
                         genome_seq: str,
                         index: Dict[str, List[Tuple[int, int]]],
                         ) -> Optional[List[Tuple[int, int]]]:
    """The single best split, or ``None`` — the thin wrapper the tests and simple callers use."""
    splits = find_microexon_splits(inserted_seq, chrom, intron_start, intron_end,
                                   strand, genome_seq, index)
    if not splits:
        return None
    return max(splits, key=lambda sp: (split_score(sp, chrom, index), [-x for x in sp[0]]))


def _cigar_insertions_beside_n(cigar: Sequence[Tuple[int, int]]) -> List[Tuple[int, int, int]]:
    """``[(insertion_index, n_index, insertion_len)]`` for every I of >= MIN_INSERTION_LEN that touches
    an N. Both orders are reported (``I N`` and ``N I``); which one is the transcript's 5' side depends
    on the strand and is the caller's business."""
    out = []
    for i, (op, ln) in enumerate(cigar):
        if op != _I or ln < MIN_INSERTION_LEN:
            continue
        if i + 1 < len(cigar) and cigar[i + 1][0] == _N:
            out.append((i, i + 1, ln))
        elif i > 0 and cigar[i - 1][0] == _N:
            out.append((i, i - 1, ln))
    return out


def rewrite_with_microexons(cigar: Sequence[Tuple[int, int]],
                            insertion_idx: int,
                            n_idx: int,
                            segments: List[Tuple[int, int]],
                            n_start: int) -> List[Tuple[int, int]]:
    """Replace the ``I`` at *insertion_idx* and the ``N`` at *n_idx* with
    ``N seg1 N seg2 … N``, in genomic order.

    *n_start* is the reference position where the N op begins. The rewrite conserves both the query
    length (the inserted bases become M) and the reference span (the N is partitioned by the segments),
    so downstream geometry is unchanged apart from the recovered exons. It never leaves an I or a D
    beside an N — the whole insertion is consumed as M — which is the ISSUE-031/038 invariant.
    """
    n_len = cigar[n_idx][1]
    n_end = n_start + n_len
    total_seg = sum(e - s for s, e in segments)
    assert total_seg == cigar[insertion_idx][1], "segments must consume the insertion exactly"
    assert segments and segments[0][0] >= n_start and segments[-1][1] <= n_end

    middle: List[Tuple[int, int]] = []
    cursor = n_start
    for (s, e) in segments:
        middle.append((_N, s - cursor))
        middle.append((_M, e - s))
        cursor = e
    middle.append((_N, n_end - cursor))
    assert all(ln > 0 for _op, ln in middle), "a zero-length op means the flanking-intron floor leaked"

    lo, hi = min(insertion_idx, n_idx), max(insertion_idx, n_idx)
    return list(cigar[:lo]) + middle + list(cigar[hi + 1:])


#: A read carrying more than this many separately-explainable junction-adjacent insertions is not a
#: micro-exon story any more; stop and leave it alone.
MAX_CALLS_PER_READ = 4


def recover_all_microexons(read, genome_seq: str, strand: str,
                           index: Dict[str, List[Tuple[int, int]]],
                           max_calls: int = MAX_CALLS_PER_READ):
    """Every junction-adjacent insertion on this read that annotated micro-exons explain, not just
    the first.

    ⚠️ THE FIRST-ONLY BEHAVIOUR WAS A REAL LOSS, and it is the reason this exists. Measured on the
    SMA panel (25,999 reads, T1 at 936ee76): 20 reads carry a SECOND explainable insertion at a
    different intron, all 20 of them on reads where the first call had already been drawn — 6.6 % of
    the 304 station-B reads. `recover_read_microexons` returned after the first hit and the rest were
    silently left as insertions.

    Each call names its own intron by COORDINATES and the calls are independent: drawing one only
    subdivides its own intron, so it cannot move another's. That is what lets the writer apply them
    one at a time against the live record, and it is why the search may be run against the original
    CIGAR and the results applied later.
    """
    import copy as _copy

    out = []
    work = read
    for _ in range(max(1, max_calls)):
        call = recover_read_microexons(work, genome_seq, strand, index)
        if call is None:
            break
        if any(c.intron == call.intron for c in out):
            break                                     # no progress; refuse to spin
        out.append(call)
        work = _copy.copy(work)
        try:
            work.cigartuples = call.new_cigar
        except Exception:
            break
    return out


def recover_read_microexons(read, genome_seq: str, strand: str,
                            index: Dict[str, List[Tuple[int, int]]],
                            ):
    """``MicroexonCall`` for the first junction-adjacent insertion this read can explain as annotated
    micro-exons, or ``None``.

    Reads the RECORD — never a TSV, never a coordinate computed from one (CLAUDE.md). One insertion
    per call: a read carrying two of them is rare enough that the caller can call again on the
    rewritten record.
    """
    cigar = list(read.cigartuples or [])
    if not cigar:
        return None
    # The index is keyed by the STANDARDIZED contig name (load_microexons), which equals the raw
    # name in any run that registered its genome's contigs; try both so a caller that did not
    # register still gets the right list instead of silently empty candidates.
    chrom = read.reference_name
    if chrom not in index:
        try:
            from ...utils.genome import standardize_chrom_name
            if standardize_chrom_name(chrom) in index:
                chrom = standardize_chrom_name(chrom)
        except Exception:
            pass
    for (i_idx, n_idx, i_len) in _cigar_insertions_beside_n(cigar):
        ref = read.reference_start
        qpos = 0
        n_start = None
        i_qstart = None
        for k, (op, ln) in enumerate(cigar):
            if k == n_idx:
                n_start = ref
            if k == i_idx:
                i_qstart = qpos
            if op in (_M, _EQ, _X, _D, _N):
                ref += ln
            if op in (_M, _EQ, _X, _I, _S):
                qpos += ln
        if n_start is None or i_qstart is None:
            continue
        seq = read.query_sequence or ''
        inserted = seq[i_qstart:i_qstart + i_len]
        if len(inserted) != i_len:
            continue
        splits = find_microexon_splits(inserted, chrom, n_start, n_start + cigar[n_idx][1],
                                       strand, genome_seq, index)
        if not splits:
            continue
        chosen, alternatives, n_tied = choose_split(splits, chrom, index, read.query_name or '')
        if not chosen:
            continue
        return MicroexonCall(
            chrom=chrom,
            segments=chosen,
            alternatives=alternatives,
            n_tied=n_tied,
            bits=split_score(chosen, chrom, index),
            new_cigar=rewrite_with_microexons(cigar, i_idx, n_idx, chosen, n_start),
            insertion_len=i_len,
            intron=(n_start, n_start + cigar[n_idx][1]),
        )
    return None

def format_calls(chrom: str, calls) -> Tuple[str, str, str, str, int]:
    """The row encoding for a list of MicroexonCall: ``'|'`` separates CALLS, ``','`` separates the
    segments within one call, ``';'`` separates a call's equally-good alternatives. Returns
    ``(segments, alternatives, intron_starts, intron_ends, max_n_tied)``; a single call round-trips
    to exactly the strings the one-call format used, so every existing consumer keeps working."""
    if not calls:
        return '', '', '', '', 0
    segs = '|'.join(format_segments(chrom, c.segments) for c in calls)
    alts = '|'.join(format_alternatives(chrom, c.alternatives) for c in calls)
    starts = ','.join(str(c.intron[0]) for c in calls)
    ends = ','.join(str(c.intron[1]) for c in calls)
    return segs, alts, starts, ends, max(c.n_tied for c in calls)


def parse_calls(segments: str, intron_starts: str, intron_ends: str):
    """Inverse of :func:`format_calls` for the writer: ``[(intron, [(s, e), ...]), ...]``."""
    if not segments:
        return []
    groups = segments.split('|')
    starts = [t for t in (intron_starts or '').split(',') if t != '']
    ends = [t for t in (intron_ends or '').split(',') if t != '']
    if not (len(groups) == len(starts) == len(ends)):
        return []
    out = []
    for g, a, b in zip(groups, starts, ends):
        segs = []
        for tok in g.split(','):
            tok = tok.strip()
            if ':' not in tok or '-' not in tok:
                continue
            try:
                span = tok.rsplit(':', 1)[1]
                lo, hi = span.split('-', 1)
                segs.append((int(lo), int(hi)))
            except ValueError:
                return []
        try:
            out.append(((int(a), int(b)), segs))
        except ValueError:
            return []
    return out
