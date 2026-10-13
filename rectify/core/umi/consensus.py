"""Consensus collapse of UMI families, using the random-primer bases as a molecular signature.

``rectify umi-dedup --umi-collapse consensus`` (the extension point left open in ``dedup.py``).

WHAT IT DOES
------------
Every read of one UMI family is a PCR copy of one original molecule. Instead of keeping one
representative (dedup mode), this mode writes one CONSENSUS fragment per family: per query
position, a quality-weighted posterior over A/C/G/T across the members, so a sequencing or
late-PCR error present in one member is outvoted by the others.

WHY THE PRIMER BASES MATTER (Lexogen CORALL V2; Kevin 2026-10-13)
-----------------------------------------------------------------
CORALL reverse transcription starts from a random primer; read 2 begins with the primer's
bases (9 nt in CORALL V2, user guide App. E: "first 9 nt of R2 have elevated error"). They
are the primer's sequence, not the RNA's, wherever the primer annealed with mismatches. Every
copy of a molecule comes from one RT event, so the copies SHARE those bases. That makes them
useful three ways:

1. **Masking.** The primer positions of a fragment are masked in the consensus
   (``--primer-action mask``: N, Q2), on WHICHEVER mate covers them. On short inserts read 1's
   3′ end reads the same bases as read 2's first 9, so read-2-only masking would leave the
   artifact in read 1 (unit 1124 found COMPASS counting it twice for exactly this reason). The
   primer span is defined in reference coordinates from the primer mate's soft-clip-corrected
   5′ end.
2. **An independent check on UMI clustering.** Directional clustering merges a rare UMI into
   a common neighbour when it looks like a sequencing error of it. A sequencing error in the
   UMI cannot change the primer, which was fixed at RT before PCR. So when a cluster holds
   several RAW UMIs and the members of two raw UMIs carry different primer sequences (each
   supported by >= ``primer_min_support`` members, differing at a base where both have
   consensus Q >= 20), the cluster merged two molecules and is split. Primer differences
   WITHIN one raw UMI are not split on (an early PCR error in the primer looks like that); they
   are counted.
3. **Mismatch classes.** At each reference position of a family with >= 2 members (each
   fragment's base = its two mates' agreeing base where they overlap):
     * unanimous non-reference base inside the primer span   -> ``primer_template``
     * unanimous non-reference base elsewhere                 -> ``molecule`` (a real variant,
       an RT error or non-templated addition, pre-PCR damage; the optional recurrence check
       moves sites most molecules share into ``recurrent``: a strain variant or a systematic
       sequencing error)
     * a minority base in a family of >= 3 that the consensus outvotes -> ``minority``
       (sequencing or late PCR error, corrected). A family of 2 cannot have a minority base.
   Each class is binned by distance from both fragment ends: the primer end (primer mate's 5′)
   and the ligation end (the other mate's 5′), which is the direct test of whether the
   mismatches crowding COMPASS's low-quality pairs are end artifacts.

The molecule key can be the dedup key (``span_key='tlen'``: R1 5′ + |TLEN|) or both
soft-clip-corrected 5′ ends (``span_key='unclipped'``, default here): primer bases left in the
reads are where aligners soft-clip differently copy to copy, which shifts TLEN but not the
unclipped ends. The run reports molecule counts under both keys.

PASS STRUCTURE (memory-bounded)
--------------------------------
1. R1 records -> fragment keys (the mate's 5′ end from the ``MC`` tag; run ``samtools fixmate``
   first) -> directional UMI clusters -> molecule ids.
2. Every read of a kept fragment is written with ``MI`` and sorted by ``MI`` (``samtools sort
   -t MI``), so each family is contiguous and chimeric or orphan mates are handled explicitly.
3. Families are streamed, consensus called, records written, then coordinate-sorted + indexed.
"""
from __future__ import annotations

import json
import logging
import math
import os
import re
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

from .clustering import canonical_umi, umi_components, umi_components_directional
from .grouping import is_spliced, leading_soft_clip, trailing_soft_clip

logger = logging.getLogger(__name__)

_MOLECULE_TAG = "MI"
_FAMILY_SIZE_TAG = "cD"
_BASES = "ACGT"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_REF_OPS = frozenset("MDN=X")
_QUERY_OPS = frozenset("MIS=X")
MAX_CONSENSUS_Q = 40
MASK_Q = 2
DIST_BINS = 50          # end-distance histogram: 0..49, then ">=50"


# ------------------------------------------------------------------------------------------ primitives
def call_column(bases: Sequence[str], quals: Sequence[int], min_posterior: float = 0.9,
                max_q: int = MAX_CONSENSUS_Q) -> Tuple[str, int, float]:
    """Quality-weighted consensus of one column.

    Each observation of base b with Phred q has error e = 10^(-q/10) (floored at Q40, so no
    single base is treated as certain). With a uniform prior over A/C/G/T, the log-likelihood
    of true base B is sum(log(1-e) if b == B else log(e/3)). Returns (base, phred, posterior);
    the base is 'N' (Phred 2) when the best posterior is below ``min_posterior`` or no
    observation is informative. Pure function.
    """
    ll = {B: 0.0 for B in _BASES}
    n = 0
    for b, q in zip(bases, quals):
        b = b.upper()
        if b not in ll:
            continue
        e = min(0.75, max(10 ** (-max(q, 0) / 10.0), 1e-4))
        good, bad = math.log(1.0 - e), math.log(e / 3.0)
        for B in _BASES:
            ll[B] += good if b == B else bad
        n += 1
    if n == 0:
        return "N", MASK_Q, 0.0
    mx = max(ll.values())
    z = sum(math.exp(v - mx) for v in ll.values())
    best = max(_BASES, key=lambda B: (ll[B], -_BASES.index(B)))
    post = math.exp(ll[best] - mx) / z
    if post < min_posterior:
        return "N", MASK_Q, post
    q = MAX_CONSENSUS_Q if post >= 1.0 else int(-10.0 * math.log10(max(1.0 - post, 1e-12)))
    return best, max(MASK_Q, min(max_q, q)), post


def parse_cigar_string(cigar: str) -> List[Tuple[str, int]]:
    return [(op, int(n)) for n, op in _CIGAR_RE.findall(cigar or "")]


def mate_five_prime_unclipped(read) -> Optional[int]:
    """Inclusive 0-based 5′ end of the MATE, soft-clip corrected, from ``next_reference_start``
    and the ``MC`` (mate CIGAR) tag. None when the mate is unmapped or MC is absent."""
    if not read.is_paired or read.mate_is_unmapped or read.next_reference_start is None or read.next_reference_start < 0:
        return None
    try:
        mc = read.get_tag("MC")
    except KeyError:
        return None
    ops = parse_cigar_string(mc)
    if not ops:
        return None
    lead = 0
    for op, n in ops:
        if op == "S":
            lead += n
        elif op == "H":
            continue
        else:
            break
    trail = 0
    for op, n in reversed(ops):
        if op == "S":
            trail += n
        elif op == "H":
            continue
        else:
            break
    reflen = sum(n for op, n in ops if op in _REF_OPS)
    start = read.next_reference_start
    if read.mate_is_reverse:
        return start + reflen - 1 + trail
    return start - lead


def five_prime_inclusive(read) -> Optional[int]:
    """Inclusive 0-based 5′ end of a mapped read, soft-clip corrected."""
    if read.is_unmapped or read.reference_start is None:
        return None
    if read.is_reverse:
        if read.reference_end is None:
            return None
        return read.reference_end - 1 + trailing_soft_clip(read)
    return read.reference_start - leading_soft_clip(read)


def unclipped_ref_coords(read) -> List[Optional[int]]:
    """Reference coordinate of every query base (BAM SEQ order). Aligned bases map to their
    position, insertions to None, soft-clipped bases to the position they would occupy
    unclipped (start - k ... / end + k), so primer bases the aligner clipped still have a place."""
    n = read.query_length or 0
    out: List[Optional[int]] = [None] * n
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if qpos is not None and rpos is not None:
            out[qpos] = rpos
    lc, tc = leading_soft_clip(read), trailing_soft_clip(read)
    for i in range(lc):
        out[i] = read.reference_start - (lc - i)
    for j in range(tc):
        out[n - tc + j] = read.reference_end + j
    return out


def primer_span(primer_mate_five: int, primer_mate_reverse: bool, length: int, source: str) -> Tuple[int, int]:
    """Inclusive reference span the primer occupied, from the primer mate's soft-clip-corrected
    5′ end. ``source='seq'``: the primer bases are still the read's first ``length`` bases.
    ``source='tag'``: they were trimmed before alignment, so the primer sat just upstream."""
    p = primer_mate_five
    if source == "seq":
        return (p - length + 1, p) if primer_mate_reverse else (p, p + length - 1)
    if source == "tag":
        return (p + 1, p + length) if primer_mate_reverse else (p - length, p - 1)
    raise ValueError(f"unknown primer source {source!r}")


def primer_bases_read_orientation(read, length: int, source: str, ref_seq: Optional[str] = None) -> Tuple[str, List[int]]:
    """The primer in READ orientation (5′->3′ as sequenced) with qualities: the read's first
    ``length`` bases (``seq``; '=' resolved against ``ref_seq``) or its ``pR``/``pQ`` tags
    (``tag``). Empty when unavailable."""
    if source == "tag":
        try:
            s = read.get_tag("pR")
        except KeyError:
            return "", []
        try:
            q = [ord(c) - 33 for c in read.get_tag("pQ")]
        except KeyError:
            q = [30] * len(s)
        return s[:length].upper(), q[:length]
    seq = resolve_seq(read, ref_seq) or ""
    qual = list(read.query_qualities or [30] * len(seq))
    if len(seq) < length:
        return "", []
    if read.is_reverse:
        s = seq[-length:][::-1].translate(_COMP)
        return s.upper(), qual[-length:][::-1]
    return seq[:length].upper(), qual[:length]


def resolve_seq(read, ref_seq: Optional[str]) -> Optional[str]:
    """Query sequence with '=' (samtools calmd -e) replaced by the reference base."""
    s = read.query_sequence
    if s is None:
        return None
    if "=" not in s:
        return s.upper()
    if ref_seq is None:
        return None
    coords = unclipped_ref_coords(read)
    out = []
    for c, x in zip(s, coords):
        if c == "=":
            out.append(ref_seq[x].upper() if x is not None and 0 <= x < len(ref_seq) else "N")
        else:
            out.append(c.upper())
    return "".join(out)


def split_by_primer(raw_umis: Sequence[str], primers: Sequence[Tuple[str, List[int]]], min_support: int = 2,
                    min_q: int = 20) -> Tuple[List[List[int]], bool]:
    """Partition a UMI cluster's members (indices) when its RAW UMIs carry different primers.

    Members are grouped by raw UMI; each raw-UMI group with >= ``min_support`` members gets a
    consensus primer. Two groups whose consensus primers differ at a position where both have
    consensus Q >= ``min_q`` are different molecules. Returns (groups, heterogeneous_within_umi)
    where ``groups`` is the partition (a single group when nothing splits) and the flag records
    a primer difference WITHIN one raw UMI (not split on; an early PCR error looks like that).
    Members of small raw-UMI groups join the supported group whose primer they match best.
    """
    idx_by_umi: Dict[str, List[int]] = defaultdict(list)
    for i, u in enumerate(raw_umis):
        idx_by_umi[u].append(i)
    cons: Dict[str, Tuple[str, List[int]]] = {}
    hetero = False
    for u, idx in idx_by_umi.items():
        ps = [primers[i] for i in idx if primers[i][0]]
        if len(ps) < min_support:
            continue
        L = min(len(p[0]) for p in ps)
        seq, qs = [], []
        for k in range(L):
            b, q, _ = call_column([p[0][k] for p in ps], [p[1][k] for p in ps])
            seq.append(b)
            qs.append(q)
            obs = Counter(p[0][k] for p in ps if p[1][k] >= min_q)
            if len([c for c in obs.values() if c >= min_support]) >= 2:
                hetero = True
        cons[u] = ("".join(seq), qs)
    if len(cons) < 2:
        return [list(range(len(raw_umis)))], hetero

    def differ(a: Tuple[str, List[int]], b: Tuple[str, List[int]]) -> bool:
        return any(x != y and x != "N" and y != "N" and qa >= min_q and qb >= min_q
                   for x, y, qa, qb in zip(a[0], b[0], a[1], b[1]))

    # Connected components over "same primer" among the supported raw UMIs.
    us = sorted(cons, key=lambda u: (-len(idx_by_umi[u]), u))
    comp = {u: u for u in us}

    def find(u):
        while comp[u] != u:
            comp[u] = comp[comp[u]]
            u = comp[u]
        return u
    for i, a in enumerate(us):
        for b in us[i + 1:]:
            if not differ(cons[a], cons[b]):
                comp[find(b)] = find(a)
    roots = sorted({find(u) for u in us}, key=lambda r: (-sum(len(idx_by_umi[u]) for u in us if find(u) == r), r))
    if len(roots) < 2:
        return [list(range(len(raw_umis)))], hetero
    groups: Dict[str, List[int]] = {r: [] for r in roots}
    for u, idx in idx_by_umi.items():
        if u in cons:
            groups[find(u)].extend(idx)
            continue
        for i in idx:   # unsupported raw UMI: join the root whose primer it matches best
            p = primers[i]

            def mism(r):
                rc = cons[r]
                return sum(1 for x, y in zip(p[0], rc[0]) if x != y and y != "N") if p[0] else 0
            best = min(roots, key=lambda r: (mism(r), roots.index(r)))
            groups[best].append(i)
    return [sorted(g) for g in groups.values() if g], hetero


# ------------------------------------------------------------------------------------------ statistics
def _dist_bin(d: Optional[int]) -> str:
    if d is None:
        return "na"
    return str(d) if d < DIST_BINS else f">={DIST_BINS}"


def _size_bin(n: int) -> str:
    return "1" if n == 1 else "2" if n == 2 else "3" if n == 3 else "4-5" if n <= 5 else "6+"


@dataclass
class ConsensusStats:
    n_input_fragments: int = 0
    n_no_umi: int = 0
    n_no_mate_info: int = 0
    n_molecules: int = 0
    n_molecules_other_key: int = 0
    span_key: str = "unclipped"
    family_size_hist: Dict[int, int] = field(default_factory=lambda: defaultdict(int))
    clusters_checked_for_primer: int = 0
    clusters_split_by_primer: int = 0
    molecules_from_splits: int = 0
    primer_hetero_within_umi: int = 0
    members_not_voting: int = 0
    reads_written: int = 0
    # class -> size bin -> qc ("pass"/"qcfail") -> count of mismatch events
    events: Dict[str, Dict[str, Dict[str, int]]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    # class -> "primer_end"/"ligation_end" -> distance bin -> count
    end_hist: Dict[str, Dict[str, Dict[str, int]]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    # positions examined (denominators): size bin -> qc -> count
    positions: Dict[str, Dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    # position denominators by end distance (for rates): end -> bin -> count, families >= 2 only
    end_positions: Dict[str, Dict[str, int]] = field(default_factory=lambda: defaultdict(lambda: defaultdict(int)))
    mismatches_before_after: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    recurrence: Dict[str, int] = field(default_factory=lambda: defaultdict(int))

    def as_dict(self) -> dict:
        def plain(d):
            if isinstance(d, dict):
                return {str(k): plain(v) for k, v in sorted(d.items(), key=lambda kv: str(kv[0]))}
            return d
        return {
            "span_key": self.span_key,
            "n_input_fragments": self.n_input_fragments,
            "n_no_umi": self.n_no_umi,
            "n_no_mate_info": self.n_no_mate_info,
            "n_molecules": self.n_molecules,
            "n_molecules_other_key": self.n_molecules_other_key,
            "duplication_rate": round(1 - self.n_molecules / self.n_input_fragments, 6) if self.n_input_fragments else 0.0,
            "family_size_hist": plain(dict(self.family_size_hist)),
            "primer": {"clusters_checked": self.clusters_checked_for_primer, "clusters_split": self.clusters_split_by_primer,
                       "molecules_from_splits": self.molecules_from_splits, "hetero_within_raw_umi": self.primer_hetero_within_umi},
            "members_not_voting": self.members_not_voting,
            "reads_written": self.reads_written,
            "mismatch_events": plain(self.events),
            "end_distance_hist": plain(self.end_hist),
            "positions_examined": plain(self.positions),
            "end_positions_examined": plain(self.end_positions),
            "mismatches_representative_vs_consensus": plain(dict(self.mismatches_before_after)),
            "recurrence": plain(dict(self.recurrence)),
        }


# ------------------------------------------------------------------------------------------ pass 1
class Frag(NamedTuple):
    qname: str
    contig: str
    strand: str
    r1_five: int
    mate_contig: Optional[str]
    mate_five: Optional[int]
    tlen: int
    tlen_five: int          # dedup.py's R1 5′ (exclusive end for reverse reads), for the tlen key
    spliced: bool
    umi: str
    score: int


def umi_of(read, umi_tag: str, from_qname: bool) -> Optional[str]:
    try:
        u = read.get_tag(umi_tag)
        if u:
            return str(u)
    except KeyError:
        pass
    if from_qname and "_" in read.query_name:
        u = read.query_name.rsplit("_", 1)[1]
        if u and set(u) <= set("ACGTN"):
            return u
    return None


def _score(read) -> int:
    q = read.query_qualities
    return int(read.mapping_quality or 0) + (int(sum(q) / len(q)) if q else 0)


def extract(bam, umi_tag: str, from_qname: bool, r1_sense: bool, stats: ConsensusStats) -> List[Frag]:
    from .grouping import five_prime_unclipped
    frags: List[Frag] = []
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_paired and not read.is_read1:
            continue
        umi = umi_of(read, umi_tag, from_qname)
        if not umi:
            stats.n_no_umi += 1
            continue
        f5 = five_prime_inclusive(read)
        if f5 is None:
            continue
        m5 = mate_five_prime_unclipped(read) if read.is_paired else None
        if read.is_paired and not read.mate_is_unmapped and m5 is None:
            stats.n_no_mate_info += 1
        sense = (not read.is_reverse) if r1_sense else read.is_reverse
        frags.append(Frag(read.query_name, read.reference_name, "+" if sense else "-", f5,
                          read.next_reference_name if (read.is_paired and not read.mate_is_unmapped) else None, m5,
                          abs(read.template_length or 0), five_prime_unclipped(read), is_spliced(read), umi, _score(read)))
    stats.n_input_fragments = len(frags)
    return frags


def cluster(frags: List[Frag], key: str, edit_distance: int, clustering: str) -> List[List[int]]:
    fn = umi_components_directional if clustering == "directional" else umi_components
    buckets: Dict[tuple, List[int]] = defaultdict(list)
    for i, f in enumerate(frags):
        if key == "unclipped":
            k = (f.contig, f.strand, f.r1_five, f.mate_contig, f.mate_five, f.spliced)
        else:
            k = (f.contig, f.strand, f.tlen_five, f.tlen, f.spliced)
        buckets[k].append(i)
    out: List[List[int]] = []
    for k in sorted(buckets, key=lambda t: tuple("" if x is None else str(x) for x in t)):
        idx = buckets[k]
        for members in fn([frags[i].umi for i in idx], edit_distance):
            out.append([idx[m] for m in members])
    return out


# ------------------------------------------------------------------------------------------ pass 3
class _Ref:
    def __init__(self, fasta: Optional[str]):
        self.fa = None
        self.cache: Dict[str, str] = {}
        if fasta:
            import pysam
            self.fa = pysam.FastaFile(fasta)

    def seq(self, contig: str) -> Optional[str]:
        if self.fa is None or contig is None:
            return None
        if contig not in self.cache:
            if contig not in self.fa.references:
                return None
            self.cache = {contig: self.fa.fetch(contig).upper()}   # keep one contig (input is MI-sorted, so this can thrash; contigs are small in yeast)
        return self.cache.get(contig)


def _mate_consensus(reads: List, ref: Optional[str], min_posterior: float) -> Tuple[object, List[str], List[int], int]:
    """Consensus of one mate's members: the modal-CIGAR group votes per query position.
    Returns (template_read, seq, quals, n_not_voting)."""
    by_cigar: Dict[str, List] = defaultdict(list)
    for r in reads:
        by_cigar[r.cigarstring].append(r)
    rep = max(reads, key=lambda r: (_score(r), r.query_name))
    modal = max(by_cigar.values(), key=lambda g: (len(g), rep in g))
    template = max(modal, key=lambda r: (_score(r), r.query_name))
    seqs = [resolve_seq(r, ref) for r in modal]
    ok = [(s, list(r.query_qualities or [30] * len(s))) for s, r in zip(seqs, modal) if s is not None]
    if not ok:
        s = resolve_seq(template, ref) or (template.query_sequence or "")
        return template, list(s), list(template.query_qualities or [30] * len(s)), len(reads) - 1
    L = min(len(s) for s, _ in ok)
    cseq, cq = [], []
    for i in range(L):
        b, q, _ = call_column([s[i] for s, _ in ok], [qq[i] for _, qq in ok], min_posterior=min_posterior)
        cseq.append(b)
        cq.append(q)
    return template, cseq, cq, len(reads) - len(ok)


def _fragment_bases(frag_reads: Dict[str, object], ref: Optional[str], contig: str) -> Dict[int, Tuple[str, int]]:
    """Per reference coordinate on ``contig``, the fragment's base (its mates agreeing where they
    overlap; 'N' where they disagree) and the better quality. A mate on another contig (a chimeric
    pair) is left out."""
    out: Dict[int, Tuple[str, int]] = {}
    for r in frag_reads.values():
        if r is None or r.is_unmapped or r.reference_name != contig:
            continue
        s = resolve_seq(r, ref)
        if s is None:
            continue
        q = r.query_qualities or [30] * len(s)
        for b, qq, x in zip(s, q, unclipped_ref_coords(r)):
            if x is None:
                continue
            if x in out:
                b0, q0 = out[x]
                out[x] = (b0, max(q0, qq)) if b0 == b else ("N", 0)
            else:
                out[x] = (b, qq)
    return out


def consensus_bam(
    input_bam: str,
    output_bam: str,
    reference: Optional[str] = None,
    umi_tag: str = "RX",
    umi_from_qname: bool = False,
    edit_distance: int = 1,
    clustering: str = "directional",
    r1_sense: bool = True,
    span_key: str = "unclipped",
    primer_mate: str = "r2",
    primer_length: int = 0,
    primer_source: str = "seq",
    primer_action: str = "mask",
    primer_min_support: int = 2,
    min_posterior: float = 0.9,
    recurrence_check: bool = False,
    threads: int = 1,
    tmp_dir: Optional[str] = None,
) -> ConsensusStats:
    """Write one consensus fragment per UMI family (see the module docstring) and return stats."""
    import pysam

    if span_key not in ("unclipped", "tlen"):
        raise ValueError(f"span_key must be 'unclipped' or 'tlen', not {span_key!r}")
    if primer_mate not in ("r1", "r2"):
        raise ValueError("primer_mate must be r1 or r2")
    if primer_source not in ("seq", "tag"):
        raise ValueError("primer_source must be seq or tag")
    if primer_action not in ("mask", "keep"):
        raise ValueError("primer_action must be mask or keep")
    stats = ConsensusStats(span_key=span_key)
    refobj = _Ref(reference)

    # ---- pass 1: fragments, clusters under both keys
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        frags = extract(bam_in, umi_tag, umi_from_qname, r1_sense, stats)
    if stats.n_no_mate_info and span_key == "unclipped":
        logger.warning("%d paired R1 records lack an MC tag; run samtools fixmate first (their mate end is unknown "
                       "and they are keyed as mate-less)", stats.n_no_mate_info)
    clusters = cluster(frags, span_key, edit_distance, clustering)
    stats.n_molecules_other_key = len(cluster(frags, "tlen" if span_key == "unclipped" else "unclipped", edit_distance, clustering))
    mol_of: Dict[str, str] = {}
    raw_umi: Dict[str, str] = {}
    for c_no, members in enumerate(clusters):
        f0 = frags[members[0]]
        mid = f"{f0.contig}:{f0.strand}:{f0.r1_five}:{c_no}"
        for i in members:
            mol_of[frags[i].qname] = mid
            raw_umi[frags[i].qname] = frags[i].umi
    del frags

    tmp_root = tempfile.mkdtemp(prefix="umicons_", dir=tmp_dir)
    tagged = os.path.join(tmp_root, "tagged.bam")
    by_mi = os.path.join(tmp_root, "by_mi.bam")
    unsorted_out = os.path.join(tmp_root, "consensus.unsorted.bam")

    # ---- pass 2: tag MI on every read of a kept fragment, sort by MI
    with pysam.AlignmentFile(input_bam, "rb") as bam_in, pysam.AlignmentFile(tagged, "wb", header=bam_in.header) as out:
        header = bam_in.header.to_dict()
        for read in bam_in:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            mid = mol_of.get(read.query_name)
            if mid is None:
                continue
            read.set_tag(_MOLECULE_TAG, mid, value_type="Z")
            read.set_tag("RX", raw_umi[read.query_name], value_type="Z")
            out.write(read)
    del mol_of
    pysam.sort("-t", _MOLECULE_TAG, "-@", str(max(threads, 1)), "-o", by_mi, tagged)
    os.remove(tagged)

    header.setdefault("PG", []).append({"ID": "rectify-umi-consensus", "PN": "rectify", "CL":
                                        f"umi-dedup --umi-collapse consensus span_key={span_key} primer={primer_mate}:"
                                        f"{primer_length}:{primer_source}:{primer_action}"})
    sites: Dict[Tuple[str, int, str], int] = defaultdict(int)

    def flush(mid: str, reads: List, out):
        frs: Dict[str, Dict[str, object]] = defaultdict(lambda: {"r1": None, "r2": None})
        for r in reads:
            frs[r.query_name]["r1" if (r.is_read1 or not r.is_paired) else "r2"] = r
        names = sorted(frs)
        groups = [names]
        if primer_length > 0 and len({raw_umi_of(frs[n]) for n in names}) >= 2:
            stats.clusters_checked_for_primer += 1
            prim = [_primer_of(frs[n]) for n in names]
            parts, hetero = split_by_primer([raw_umi_of(frs[n]) for n in names], prim, primer_min_support)
            if hetero:
                stats.primer_hetero_within_umi += 1
            if len(parts) > 1:
                stats.clusters_split_by_primer += 1
                stats.molecules_from_splits += len(parts)
                groups = [[names[i] for i in p] for p in parts]
        for gi, g in enumerate(groups):
            _emit(mid if len(groups) == 1 else f"{mid}/p{gi + 1}", [frs[n] for n in g], out)

    def raw_umi_of(fr) -> str:
        r = fr["r1"] or fr["r2"]
        try:
            return str(r.get_tag("RX"))
        except KeyError:
            return ""

    def _primer_of(fr) -> Tuple[str, List[int]]:
        r = fr[primer_mate]
        if r is None or r.is_unmapped:
            return "", []
        return primer_bases_read_orientation(r, primer_length, primer_source, refobj.seq(r.reference_name))

    def _emit(mid: str, members: List[Dict[str, object]], out):
        n = len(members)
        stats.n_molecules += 1
        stats.family_size_hist[n] += 1
        anchor = next((m["r1"] or m["r2"]) for m in members)
        contig = anchor.reference_name
        ref = refobj.seq(contig)
        qc = "qcfail" if any((m["r1"] is not None and m["r1"].is_qcfail) or (m["r2"] is not None and m["r2"].is_qcfail)
                             for m in members) else "pass"
        # fragment ends (reference, inclusive): primer end = primer mate's 5′, ligation end = the other mate's 5′
        pm = next((m[primer_mate] for m in members if m[primer_mate] is not None and not m[primer_mate].is_unmapped), None)
        other = "r1" if primer_mate == "r2" else "r2"
        om = next((m[other] for m in members if m[other] is not None and not m[other].is_unmapped), None)
        p_end = five_prime_inclusive(pm) if pm is not None and pm.reference_name == contig else None
        l_end = five_prime_inclusive(om) if om is not None and om.reference_name == contig else None
        span = primer_span(p_end, pm.is_reverse, primer_length, primer_source) if (p_end is not None and primer_length > 0) else None

        # ---- classification over fragment-level bases
        fb = [_fragment_bases(m, ref, contig) for m in members]
        sb = _size_bin(n)
        if ref is not None:
            xs = set()
            for d in fb:
                xs.update(d)
            for x in xs:
                if not (0 <= x < len(ref)):
                    continue
                rb = ref[x]
                obs = [d[x][0] for d in fb if x in d and d[x][0] in _BASES]
                if not obs:
                    continue
                stats.positions[sb][qc] += 1
                dp = abs(x - p_end) if p_end is not None else None
                dl = abs(x - l_end) if l_end is not None else None
                if n >= 2:
                    stats.end_positions["primer_end"][_dist_bin(dp)] += 1
                    stats.end_positions["ligation_end"][_dist_bin(dl)] += 1
                alts = [b for b in obs if b != rb]
                if not alts:
                    continue
                in_primer = span is not None and span[0] <= x <= span[1]
                if len(obs) >= 2 and len(alts) == len(obs) and len(set(alts)) == 1:
                    cls = "primer_template" if in_primer else "molecule"
                    if cls == "molecule":
                        sites[(contig, x, alts[0])] += 1
                elif n >= 3 and len(alts) < len(obs) and Counter(obs).most_common(1)[0][0] == rb:
                    cls = "minority"
                elif n == 1:
                    cls = "singleton_in_primer" if in_primer else "singleton"
                else:
                    cls = "mixed"
                k = len(alts) if cls == "minority" else 1
                stats.events[cls][sb][qc] += k
                stats.end_hist[cls]["primer_end"][_dist_bin(dp)] += k
                stats.end_hist[cls]["ligation_end"][_dist_bin(dl)] += k

        # ---- consensus records, one per mate
        for mate in ("r1", "r2"):
            reads = [m[mate] for m in members if m[mate] is not None]
            if not reads:
                continue
            if n == 1:
                tpl = reads[0]
                s = resolve_seq(tpl, ref)
                seq = list(s if s is not None else (tpl.query_sequence or ""))
                qual = list(tpl.query_qualities or [30] * len(seq))
                not_voting = 0
            else:
                tpl, seq, qual, not_voting = _mate_consensus(reads, ref, min_posterior)
                stats.members_not_voting += not_voting
            if len(seq) != (tpl.query_length or 0):
                seq = list(resolve_seq(tpl, ref) or tpl.query_sequence or "")
                qual = list(tpl.query_qualities or [30] * len(seq))
            masked = 0
            if span is not None and primer_action == "mask" and tpl.reference_name == contig:
                for i, x in enumerate(unclipped_ref_coords(tpl)):
                    if x is not None and span[0] <= x <= span[1]:
                        seq[i], qual[i] = "N", MASK_Q
                        masked += 1
            if ref is not None and n >= 2 and not tpl.is_unmapped:
                coords = unclipped_ref_coords(tpl)
                rs = resolve_seq(tpl, ref) or ""
                before = sum(1 for b, x in zip(rs, coords) if x is not None and 0 <= x < len(ref) and b in _BASES and b != ref[x])
                after = sum(1 for b, x in zip(seq, coords) if x is not None and 0 <= x < len(ref) and b in _BASES and b != ref[x])
                stats.mismatches_before_after[f"{sb}|representative"] += before
                stats.mismatches_before_after[f"{sb}|consensus"] += after
            new = pysam.AlignedSegment(out.header)
            new.query_name = tpl.query_name
            new.flag = tpl.flag
            new.reference_id = tpl.reference_id
            new.reference_start = tpl.reference_start
            new.mapping_quality = tpl.mapping_quality
            new.cigartuples = tpl.cigartuples
            new.next_reference_id = tpl.next_reference_id
            new.next_reference_start = tpl.next_reference_start
            new.template_length = tpl.template_length
            new.query_sequence = "".join(seq)
            new.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in qual))
            for tag, val, typ in tpl.get_tags(with_value_type=True):
                if tag in ("MD", "NM", _MOLECULE_TAG, _FAMILY_SIZE_TAG):
                    continue
                if typ == "B":          # arrays: pysam infers the subtype from the array itself
                    new.set_tag(tag, val)
                else:
                    new.set_tag(tag, val, value_type=typ)
            new.set_tag(_MOLECULE_TAG, mid, value_type="Z")
            new.set_tag(_FAMILY_SIZE_TAG, n, value_type="i")
            new.set_tag("cX", int(not_voting), value_type="i")
            new.set_tag("pK", int(masked), value_type="i")
            out.write(new)
            stats.reads_written += 1

    # ---- pass 3: stream families
    with pysam.AlignmentFile(by_mi, "rb") as bam_in, pysam.AlignmentFile(unsorted_out, "wb", header=header) as out:
        cur, buf = None, []
        for read in bam_in:
            mid = read.get_tag(_MOLECULE_TAG)
            if mid != cur:
                if buf:
                    flush(cur, buf, out)
                cur, buf = mid, []
            buf.append(read)
        if buf:
            flush(cur, buf, out)
    os.remove(by_mi)
    pysam.sort("-@", str(max(threads, 1)), "-o", output_bam, unsorted_out)
    os.remove(unsorted_out)
    pysam.index(output_bam)
    try:
        os.rmdir(tmp_root)
    except OSError:
        pass

    # ---- optional recurrence check: how many consensus molecules covering a 'molecule' site carry the same base
    if recurrence_check and sites and reference:
        with pysam.AlignmentFile(output_bam, "rb") as bam:
            for (contig, x, alt), nfam in sites.items():
                carriers, total = set(), set()
                for col in bam.pileup(contig, x, x + 1, truncate=True, stepper="nofilter", min_base_quality=0,
                                      ignore_orphans=False, max_depth=100000):
                    for pr in col.pileups:
                        if pr.is_del or pr.is_refskip or pr.query_position is None:
                            continue
                        mi = pr.alignment.get_tag(_MOLECULE_TAG)
                        total.add(mi)
                        if pr.alignment.query_sequence[pr.query_position].upper() == alt:
                            carriers.add(mi)
                frac = len(carriers) / len(total) if total else 0.0
                k = "recurrent" if (len(total) >= 3 and frac >= 0.5) else "molecule_specific"
                stats.recurrence[k] += 1
                stats.recurrence[f"{k}_families"] += nfam
    logger.info("UMI consensus %s: %d fragments -> %d molecules (other key %d); %d reads written",
                input_bam, stats.n_input_fragments, stats.n_molecules, stats.n_molecules_other_key, stats.reads_written)
    return stats


def write_stats(stats: ConsensusStats, path: str) -> None:
    with open(path, "w") as fh:
        json.dump(stats.as_dict(), fh, indent=1)
