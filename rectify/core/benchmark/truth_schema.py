#!/usr/bin/env python3
"""Per-read ground-truth schema for the RECTIFY simulation benchmark (component 2).

Built FIRST because it constrains the simulator wrapper (component 1) and the
scorer (component 3). The schema carries, per simulated read:

* **Exact-indel truth** as a per-position true edit with an **ambiguity-equivalence
  span** (the net-indel-in-run model). The framing metric is EXACT INDEL-POSITION
  CONCORDANCE WITH TRUTH, ambiguity-aware (a homopolymer/STR deletion may sit
  anywhere in its run), NEVER edit distance (tied by construction at the contested
  positions). Each ``IndelTruth`` therefore records both the canonical leftmost
  position AND the genomic span within which the net indel is positionally
  equivalent, so the scorer can credit any placement inside the span as a TP.
* **Junction truth** — donor/acceptor (intron ``[start, end)`` half-open genome
  coords, left-normalized via ``chimeric_consensus.normalize_junction`` so the
  scorer's membership test is ambiguity-aware), the canonical/non-canonical motif,
  and the **NIC/NNC class label** (required from P0 — without it the co-primary
  discovery FDR is unmeasurable). NIC = novel combination of *known* sites; NNC =
  at least one novel (uncatalogued) splice site.
* **CPA truth** — the 3' cleavage genome coord + downstream genomic-A context.
* **Standing-variant truth (C6)** — known SNP/indel set per read (het/hom, VAF,
  distance to nearest junction) for the variant-induced-FDR stratum.
* The **single shared genomic-region-disjoint TRAIN/TEST split tag** (§7) — a
  first-class column so a calibration win cannot leak across concepts.
* C1 **HP/STR cell metadata** (``run_unit``/``run_copies`` per indel) so cell
  coverage against ``min_count=100`` (``hp_penalty.py:184,213``) is auditable.

Truth is NEVER overloaded onto the existing ``XV`` aux tag (that is a
correction-CATEGORY label). Truth is a sidecar table keyed by ``read_id``.

Serialization: a TSV with JSON-encoded list columns (junctions / indels /
variants). Stdlib only, plus the in-tree ambiguity primitives. Round-trips
losslessly via ``write_truth_table`` / ``read_truth_table``.

Author: Kevin R. Roy
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable, Union

# In-tree ambiguity-aware primitives (REUSE — do not reimplement). Importing the
# schema must not pull heavy deps, and chimeric_consensus only needs pysam/numpy
# at module import, which are already RECTIFY deps.
from ..consensus.chimeric_consensus import (
    normalize_junction,
    junction_ambiguity_window,
    _canonical_within_window,
)

SCHEMA_VERSION = "1.0"


# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------
class JunctionClass(str, Enum):
    """Discovery class of a truth junction (§8 — separate FDR tracks).

    ANNOTATED  both sites + the pairing are catalogued (in the GTF/GFF).
    NIC        Novel In Catalog: both donor AND acceptor are individually
               catalogued, but THIS pairing is novel (a new isoform from known
               sites). Recall here measures isoform discovery without new sites.
    NNC        Novel Not in Catalog: at least one splice site is uncatalogued.
               This is the hardest, FDR-sensitive track (esp. non-canonical NNC).
    """

    ANNOTATED = "ANNOTATED"
    NIC = "NIC"
    NNC = "NNC"


class IndelKind(str, Enum):
    DEL = "DEL"  # bases present in reference, absent in read (consumes ref)
    INS = "INS"  # bases present in read, absent in reference (consumes query)


class SplitTag(str, Enum):
    """The ONE shared genomic-region-disjoint train/test split (§7).

    A locus/region used to ESTIMATE any calibration table is TRAIN; ALL ablation
    scoring is on TEST. The split is the SAME across every facet so a win cannot
    leak across concepts.
    """

    TRAIN = "TRAIN"
    TEST = "TEST"


# ---------------------------------------------------------------------------
# Per-feature truth records
# ---------------------------------------------------------------------------
@dataclass
class JunctionTruth:
    """One true intron crossed by a read.

    Coordinates are 0-based half-open genome coords for the intron span
    ``[intron_start, intron_end)`` (the N-op convention used throughout RECTIFY:
    ``Event.r_start``/``r_end``). Stored **left-normalized** (leftmost
    ambiguity-equivalent placement) so the scorer's set-membership test is
    ambiguity-aware by construction.
    """

    intron_start: int
    intron_end: int
    strand: str
    canonical: bool          # any ambiguity-equivalent placement is GT/GC..AG
    klass: JunctionClass
    motif: str = ""          # e.g. "GT-AG", "AT-AC", "GT-TG" (genome-strand donor..acceptor)
    l_amb: int = 0           # how far the junction may slide LEFT (bp)
    r_amb: int = 0           # how far it may slide RIGHT (bp)

    @classmethod
    def from_intron(cls, intron_start: int, intron_end: int, strand: str,
                    klass: Union[JunctionClass, str], genome_seq: str,
                    motif: str = "") -> "JunctionTruth":
        """Build from raw intron coords + the contig sequence: left-normalize,
        compute the ambiguity window, and decide canonicity via the in-tree
        ambiguity-aware primitives."""
        ns, ne = normalize_junction(intron_start, intron_end, genome_seq)
        l_amb, r_amb = junction_ambiguity_window(ns, ne, genome_seq)
        canonical = _canonical_within_window(ns, ne, genome_seq, l_amb, r_amb)
        if not motif:
            donor = genome_seq[ns:ns + 2].upper() if ns + 2 <= len(genome_seq) else ""
            acceptor = genome_seq[ne - 2:ne].upper() if ne - 2 >= 0 else ""
            motif = f"{donor}-{acceptor}"
        if isinstance(klass, str):
            klass = JunctionClass(klass)
        return cls(ns, ne, strand, canonical, klass, motif, l_amb, r_amb)


@dataclass
class IndelTruth:
    """One true indel carried by a read, with its ambiguity-equivalence span.

    ``pos`` is the canonical **leftmost** genome position of the edit (0-based).
    ``length`` is the NET base count (always positive; sign carried by ``kind``).
    ``eq_start``/``eq_end`` bound the genomic span within which the net indel is
    positionally equivalent (the HP run, the STR period-tiled region, or just
    ``[pos, pos+length)`` for a unique-context indel). The scorer credits any
    aligner placement whose net (D-I) within ``[eq_start, eq_end)`` equals the
    truth net AND which introduces no indel outside the span (mirrors the
    vertical-slice ``net_indel_in_run`` ambiguity-aware match).
    """

    pos: int
    length: int
    kind: IndelKind
    eq_start: int
    eq_end: int
    context: str = "UNIQUE"      # HP | STR | UNIQUE | BOUNDARY
    run_unit: str = ""           # "A" for an HP run; "AG" for an STR dinuc unit
    run_copies: int = 0          # HP run length, or STR copy number (for C1 cells)

    @property
    def base_class(self) -> str:
        """AT vs CG class of the run unit (the C1 ``base_class`` axis). Empty for
        non-HP/STR context."""
        if not self.run_unit:
            return ""
        u = self.run_unit.upper()
        if all(c in "AT" for c in u):
            return "AT"
        if all(c in "CG" for c in u):
            return "CG"
        return "MIXED"


@dataclass
class VariantTruth:
    """A standing germline/somatic variant carried by a read (C6 stratum).

    A variant near a splice site can fabricate a 'novel' non-canonical junction
    or be mis-charged as a mismatch the DP 'repairs' with a spurious micro-indel,
    inflating the de-novo discovery FDR. This record lets the scorer measure
    variant-induced FDR stratified het/hom and variant-adjacent vs -distant.
    """

    pos: int                       # 0-based genome position
    ref_allele: str
    alt_allele: str
    zygosity: str = "HET"          # HET | HOM
    vaf: float = 0.5               # non-Mendelian VAFs allowed (aneuploid A549)
    dist_to_junction: Optional[int] = None  # bp to nearest truth junction (None=far)


# ---------------------------------------------------------------------------
# The per-read truth row
# ---------------------------------------------------------------------------
@dataclass
class ReadTruth:
    """Ground truth for one simulated read.

    The exact read<->genome alignment is captured in ``true_cigar`` (a SAM CIGAR
    string over ``[genome_start, genome_end)``). For Tier-1 constructed reads it
    is known by construction; for Tier-2 pbsim3 reads it is composed from the
    per-read MAF (read<->transcript) and the known transcript<->genome exon
    structure. ``junctions`` / ``indels`` are the *salient* truth features the
    scorer evaluates (derived from / consistent with ``true_cigar``).
    """

    read_id: str
    true_locus: str                       # locus / gene identifier
    true_transcript: str                  # transcript-of-origin (simulator reports this)
    chrom: str
    strand: str
    genome_start: int                     # 0-based, half-open read span on genome
    genome_end: int
    true_cigar: str = ""                  # SAM CIGAR, read<->genome (exact-edit source)
    junctions: List[JunctionTruth] = field(default_factory=list)
    indels: List[IndelTruth] = field(default_factory=list)
    variants: List[VariantTruth] = field(default_factory=list)
    true_cpa: Optional[int] = None        # 3' cleavage genome coord (0-based)
    downstream_a_count: Optional[int] = None  # genomic-A run length 3' of true CPA (C2)
    stratum: str = ""                     # HP | GENOMIC_A_CPA | NIC_NNC | STR |
                                          # PARALOG | PANEL_FAILURE | COVERAGE_Q | VARIANT
    split: SplitTag = SplitTag.TEST
    coverage: Optional[int] = None        # depth of the locus pool the read belongs to
    q_calibration: str = "calibrated"     # calibrated | miscalibrated (C3 soft-decision)
    panel_unplaced_intended: bool = False  # constructed C5 panel-failure read
    schema_version: str = SCHEMA_VERSION

    # ---- HP/STR cell accounting (C1 min_count=100 audit) ------------------
    def hp_cells(self) -> List[Tuple[str, int, str]]:
        """The (base_class, run_copies, context) C1 cells this read exercises —
        used to verify every cell clears ``min_count=100`` across the corpus."""
        cells = []
        for ind in self.indels:
            if ind.context in ("HP", "STR") and ind.run_copies:
                cells.append((ind.base_class, ind.run_copies, ind.context))
        return cells


# ---------------------------------------------------------------------------
# (De)serialization — TSV with JSON-encoded list columns
# ---------------------------------------------------------------------------
_COLUMNS = [
    "read_id", "true_locus", "true_transcript", "chrom", "strand",
    "genome_start", "genome_end", "true_cigar", "junctions", "indels",
    "variants", "true_cpa", "downstream_a_count", "stratum", "split",
    "coverage", "q_calibration", "panel_unplaced_intended", "schema_version",
]


def _enc_junctions(js: List[JunctionTruth]) -> str:
    return json.dumps([
        {**asdict(j), "klass": j.klass.value} for j in js
    ], separators=(",", ":"))


def _enc_indels(ins: List[IndelTruth]) -> str:
    return json.dumps([
        {**asdict(i), "kind": i.kind.value} for i in ins
    ], separators=(",", ":"))


def _enc_variants(vs: List[VariantTruth]) -> str:
    return json.dumps([asdict(v) for v in vs], separators=(",", ":"))


def _row_to_fields(r: ReadTruth) -> List[str]:
    def s(x):
        return "" if x is None else str(x)
    return [
        r.read_id, r.true_locus, r.true_transcript, r.chrom, r.strand,
        s(r.genome_start), s(r.genome_end), r.true_cigar,
        _enc_junctions(r.junctions), _enc_indels(r.indels),
        _enc_variants(r.variants), s(r.true_cpa), s(r.downstream_a_count),
        r.stratum, r.split.value, s(r.coverage), r.q_calibration,
        "1" if r.panel_unplaced_intended else "0", r.schema_version,
    ]


def write_truth_table(rows: Iterable[ReadTruth], path: Union[str, Path]) -> int:
    """Write truth rows as a TSV (JSON-encoded list columns). Returns row count."""
    path = Path(path)
    n = 0
    with path.open("w") as fh:
        fh.write("\t".join(_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(_row_to_fields(r)) + "\n")
            n += 1
    return n


def _dec_junctions(blob: str) -> List[JunctionTruth]:
    if not blob:
        return []
    out = []
    for d in json.loads(blob):
        d = dict(d)
        d["klass"] = JunctionClass(d["klass"])
        out.append(JunctionTruth(**d))
    return out


def _dec_indels(blob: str) -> List[IndelTruth]:
    if not blob:
        return []
    out = []
    for d in json.loads(blob):
        d = dict(d)
        d["kind"] = IndelKind(d["kind"])
        out.append(IndelTruth(**d))
    return out


def _dec_variants(blob: str) -> List[VariantTruth]:
    if not blob:
        return []
    return [VariantTruth(**d) for d in json.loads(blob)]


def _opt_int(x: str) -> Optional[int]:
    return int(x) if x not in ("", "None") else None


def read_truth_table(path: Union[str, Path]) -> List[ReadTruth]:
    """Load a truth TSV back into ``ReadTruth`` objects (lossless round-trip)."""
    path = Path(path)
    rows: List[ReadTruth] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(ReadTruth(
                read_id=f[idx["read_id"]],
                true_locus=f[idx["true_locus"]],
                true_transcript=f[idx["true_transcript"]],
                chrom=f[idx["chrom"]],
                strand=f[idx["strand"]],
                genome_start=int(f[idx["genome_start"]]),
                genome_end=int(f[idx["genome_end"]]),
                true_cigar=f[idx["true_cigar"]],
                junctions=_dec_junctions(f[idx["junctions"]]),
                indels=_dec_indels(f[idx["indels"]]),
                variants=_dec_variants(f[idx["variants"]]),
                true_cpa=_opt_int(f[idx["true_cpa"]]),
                downstream_a_count=_opt_int(f[idx["downstream_a_count"]]),
                stratum=f[idx["stratum"]],
                split=SplitTag(f[idx["split"]]),
                coverage=_opt_int(f[idx["coverage"]]),
                q_calibration=f[idx["q_calibration"]],
                panel_unplaced_intended=f[idx["panel_unplaced_intended"]] == "1",
                schema_version=f[idx["schema_version"]] if "schema_version" in idx else SCHEMA_VERSION,
            ))
    return rows


# ---------------------------------------------------------------------------
# Helpers for builders (Tier-1 constructors + Tier-2 MAF composer)
# ---------------------------------------------------------------------------
def homopolymer_run(genome_seq: str, pos: int) -> Tuple[int, int, str]:
    """Return ``(run_start, run_end, base)`` of the maximal homopolymer run that
    contains ``pos``. ``[run_start, run_end)`` half-open. Used to compute the
    ambiguity-equivalence span of an in-run indel."""
    if pos < 0 or pos >= len(genome_seq):
        return pos, pos, ""
    base = genome_seq[pos].upper()
    s = pos
    while s > 0 and genome_seq[s - 1].upper() == base:
        s -= 1
    e = pos + 1
    while e < len(genome_seq) and genome_seq[e].upper() == base:
        e += 1
    return s, e, base


def str_run(genome_seq: str, pos: int, unit_len: int) -> Tuple[int, int, str, int]:
    """Return ``(run_start, run_end, unit, copies)`` for a short-tandem-repeat of
    period ``unit_len`` containing ``pos``. Falls back to a 1-copy window when no
    repeat is detected."""
    if unit_len < 1 or pos < 0 or pos + unit_len > len(genome_seq):
        return pos, pos, "", 0
    unit = genome_seq[pos:pos + unit_len].upper()
    s = pos
    while s - unit_len >= 0 and genome_seq[s - unit_len:s].upper() == unit:
        s -= unit_len
    e = pos + unit_len
    while e + unit_len <= len(genome_seq) and genome_seq[e:e + unit_len].upper() == unit:
        e += unit_len
    copies = (e - s) // unit_len
    return s, e, unit, copies


def deletion_ambiguity_span(genome_seq: str, pos: int, length: int) -> Tuple[int, int]:
    """Genomic span ``[eq_start, eq_end)`` covering every ambiguity-EQUIVALENT
    placement of a ``length``-bp deletion that starts at ``pos`` — the SAME
    base-equality slide used for junctions (``chimeric_consensus``): the deletion
    slides left while ``seq[s-1] == seq[s-1+length]`` and right while
    ``seq[p] == seq[p+length]``, preserving the read.

    This generalizes homopolymer/STR run detection: it also captures HALF-period
    bleeds (e.g. a 'T' immediately before an (AT)n run extends an 'AT' deletion's
    equivalence one base left), which full-period run tiling misses and which
    otherwise charge an equivalent placement as a false error.
    """
    s = pos
    while (s - 1 >= 0 and s - 1 + length < len(genome_seq)
           and genome_seq[s - 1].upper() == genome_seq[s - 1 + length].upper()):
        s -= 1
    p = pos
    while (p + length < len(genome_seq)
           and genome_seq[p].upper() == genome_seq[p + length].upper()):
        p += 1
    return s, p + length


def make_hp_indel(genome_seq: str, run_start: int, length: int,
                  kind: IndelKind = IndelKind.DEL) -> IndelTruth:
    """Construct an ``IndelTruth`` for a homopolymer-run indel, computing the
    ambiguity-equivalence span from the base-equality slide (== the run for a pure
    HP, but also captures flank bleeds)."""
    rs, re, base = homopolymer_run(genome_seq, run_start)
    es, ee = deletion_ambiguity_span(genome_seq, rs, length) if kind == IndelKind.DEL else (rs, re)
    return IndelTruth(
        pos=es, length=length, kind=kind, eq_start=min(es, rs), eq_end=max(ee, re),
        context="HP", run_unit=base, run_copies=re - rs,
    )


def make_unique_indel(pos: int, length: int,
                      kind: IndelKind = IndelKind.DEL) -> IndelTruth:
    """Construct an ``IndelTruth`` in unique context (ambiguity span = the edit
    itself, no in-run slide)."""
    return IndelTruth(pos=pos, length=length, kind=kind,
                      eq_start=pos, eq_end=pos + length, context="UNIQUE")
