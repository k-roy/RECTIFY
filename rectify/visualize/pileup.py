"""
Compressed read views, the union-exon (spliced) axis, isoform chains and
junction arcs -- rbrowse's read-alignment vocabulary on the PUBLICATION path.

rbrowse (``~/work/rbrowse``, the ONT long-read browser) draws reads in RNA
space with three representations this module reproduces in matplotlib, from
the same read fields the slicer emits, in the house tokens:

* **merged row** (``merged_reads``): a group of reads compressed to one raster
  strip. Per column: the fraction of reads whose ALIGNED BODY covers it sets the
  alpha of the body colour; the A-rich poly(A) PREFIX of the 3' clip paints the
  tail colour where tails are at least a quarter of the bodies present (or the
  bodies have ended); the non-A clip remainder paints ``mute`` -- "unaligned,
  not poly(A)". A per-column MAJORITY is what a compressed raster can honestly
  claim (rbrowse, 2026-08-21, the rDNA lesson: a minority tail must not repaint
  the truth).
* **ribbon** (``ribbon``): one 3'-end cluster as a WEDGE -- column height = the
  fraction of the cluster's reads present at that x, so the free 5' ends render
  as a real ECDF ramp; the poly(A) prefix stacks on top; the shared anchor gets
  a sharp line; a 1-px ink rule sits at the body/tail junction because that
  boundary is the most-read edge in the tool and a rule survives any palette.
* **spliced axis** (``SplicedAxis``): union-exon space. Cut points = every
  observed junction end + every annotated exon/intron boundary; a segment is
  KEPT if any read has bases in it or it is annotated-exonic, else collapsed to
  a fixed gap; kept non-exonic segments longer than ``seg_cap_bp`` are
  compressed and the ruler says so. ``to_t`` / ``from_t`` are the piecewise-
  linear maps; pass ``xform=axis.to_t`` to the drawing functions.
* **isoform chains** (``chain_clusters``): reads grouped by their EXACT
  junction chain plus retention tokens, named relative to the canonical exon
  structure (FL, dE<k>, ret-i<k>, alt5'SS/alt3'SS, novel, unspliced, shared),
  DRS-suffix compatibility, the observed-dominant chain starred.
* **junction arcs** (``junction_counts`` + ``junction_arcs``): a sashimi-style
  view -- counts per junction as arcs whose width scales with support.

THE COLOUR RULE HERE (tokens.json ``layers``): a read body is ``subtle`` grey,
or the SAMPLE's layer-A role when samples are compared (signal rule); the
poly(A) tail and the 3'-end anchor are the ``polya`` identity -- ONE hue for
the 3' end, its mark and its tail (rbrowse paints tails base-A green on screen;
the publication path keeps one identity per hue); the soft-clip remainder is
``mute``; isoform IDENTITY is a LETTER, never a hue (Kevin, rbrowse
2026-08-12: "too many colors on the figure" -- the cluster palette left the
figure and the letter is the sole cross-panel key); strand is position.

Author: Kevin R. Roy
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    from matplotlib.patches import Rectangle
    MATPLOTLIB_AVAILABLE = True
except ImportError:  # pragma: no cover
    MATPLOTLIB_AVAILABLE = False

from . import tokens as TOK
from .tracks import Transcript, _apply_xlim, _as_region, arc

__all__ = [
    "Read", "merged_reads", "ribbon", "SplicedAxis", "chain_clusters", "Chain", "isoform_rows",
    "junction_counts", "junction_arcs",
]

Interval = Tuple[int, int]
XForm = Optional[Callable[[float], float]]


# ============================================================================
# The read carrier -- rbrowse's slicer fields, plain
# ============================================================================
@dataclass
class Read:
    """One aligned read in genomic coordinates (0-based, half-open).

    ``blocks`` are the aligned blocks as ``(start, length)`` -- rbrowse's ``bl``;
    ``junctions`` the introns as ``(start, length)`` -- rbrowse's ``n`` (derived from
    the block gaps when omitted); ``tail`` the A-rich prefix length of the 3' soft
    clip -- rbrowse's ``tp``, the tail gate; ``clip3`` the full 3' soft-clip length
    (so ``clip3 - tail`` is the non-A remainder). ``sample`` is free text.
    """
    id: str
    strand: str
    start: int
    end: int
    blocks: List[Interval]
    junctions: Optional[List[Interval]] = None
    tail: int = 0
    clip3: int = 0
    clip5: int = 0
    sample: Optional[str] = None

    def __post_init__(self):
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}")
        self.blocks = [(int(s), int(l)) for s, l in self.blocks]
        if not self.blocks:
            self.blocks = [(int(self.start), int(self.end - self.start))]
        if self.junctions is None:
            self.junctions = [(self.blocks[i][0] + self.blocks[i][1],
                               self.blocks[i + 1][0] - (self.blocks[i][0] + self.blocks[i][1]))
                              for i in range(len(self.blocks) - 1)
                              if self.blocks[i + 1][0] > self.blocks[i][0] + self.blocks[i][1]]
        self.tail = int(max(0, self.tail))
        self.clip3 = int(max(self.clip3, self.tail))

    @property
    def three_prime(self) -> int:
        """The aligned 3' end coordinate (the anchor)."""
        return self.end if self.strand == "+" else self.start

    @property
    def five_prime(self) -> int:
        return self.start if self.strand == "+" else self.end

    def has_bases_in(self, s: int, e: int) -> bool:
        return any(bs < e and bs + bl > s for bs, bl in self.blocks)

    def spans_point(self, p: int) -> bool:
        """Does an aligned block run CONTIGUOUSLY through position p?"""
        return any(bs < p < bs + bl for bs, bl in self.blocks)

    @classmethod
    def from_rbrowse(cls, rec: dict) -> "Read":
        """From a slicer record (``id st s e bl [n] tp cr cl``)."""
        clip3 = len(rec.get("cr", "") if rec.get("st", "+") == "+" else rec.get("cl", ""))
        clip5 = len(rec.get("cl", "") if rec.get("st", "+") == "+" else rec.get("cr", ""))
        return cls(rec["id"], rec["st"], rec["s"], rec["e"], [tuple(b) for b in rec["bl"]],
                   junctions=[tuple(j) for j in rec.get("n", [])] or None,
                   tail=rec.get("tp", 0) or 0, clip3=clip3, clip5=clip5, sample=rec.get("sample"))

    @classmethod
    def from_aligned_segment(cls, seg, *, tail: int = 0, sample: Optional[str] = None,
                             strand: Optional[str] = None) -> "Read":
        """From a ``pysam.AlignedSegment`` (blocks from ``get_blocks``, clips from the
        CIGAR). ``tail`` is the caller's poly(A) prefix length -- calling it is the
        slicer's job, not this module's."""
        st = strand or ("-" if seg.is_reverse else "+")
        blocks = [(s, e - s) for s, e in seg.get_blocks()]
        cig = seg.cigartuples or []
        left = cig[0][1] if cig and cig[0][0] == 4 else 0
        right = cig[-1][1] if cig and cig[-1][0] == 4 else 0
        return cls(seg.query_name, st, seg.reference_start, seg.reference_end, blocks,
                   tail=tail, clip3=(right if st == "+" else left), clip5=(left if st == "+" else right),
                   sample=sample)


# ============================================================================
# Column accumulation (difference arrays, O(blocks))
# ============================================================================
def _columns(reads: Sequence[Read], x0: float, x1: float, nbins: int, xform: XForm = None):
    """Per-bin counts of body / tail / non-A clip across [x0, x1) (in display space
    when ``xform`` is given). Returns three arrays of length nbins."""
    X = xform or (lambda p: p)
    lo, hi = (X(x0), X(x1)) if xform else (x0, x1)
    lo, hi = min(lo, hi), max(lo, hi)
    span = max(hi - lo, 1e-9)
    body = np.zeros(nbins + 1); tail = np.zeros(nbins + 1); clip = np.zeros(nbins + 1)

    def add(arr, a, b):
        ia = int(np.floor((min(a, b) - lo) / span * nbins)); ib = int(np.floor((max(a, b) - lo) / span * nbins))
        if ib < 0 or ia > nbins:
            return
        ia = max(0, min(nbins, ia)); ib = max(0, min(nbins, ib))
        arr[ia] += 1; arr[ib] -= 1

    for r in reads:
        for bs, bl in r.blocks:
            add(body, X(bs), X(bs + bl))
        if r.tail > 0:
            t0 = r.end if r.strand == "+" else r.start - r.tail
            add(tail, X(t0), X(t0 + r.tail))
        rest = r.clip3 - r.tail
        if rest > 0:
            c0 = r.end + r.tail if r.strand == "+" else r.start - r.clip3
            add(clip, X(c0), X(c0 + rest))
    return np.cumsum(body)[:nbins], np.cumsum(tail)[:nbins], np.cumsum(clip)[:nbins]


def _nbins_for(ax, x0, x1) -> int:
    fig = ax.figure
    width_px = ax.get_position().width * fig.get_figwidth() * fig.dpi
    return int(max(8, min(4000, round(width_px))))


def _body_color(role: Optional[str]) -> str:
    return TOK.color("subtle") if role is None else TOK.role(role)


# ============================================================================
# merged row
# ============================================================================
def merged_reads(ax, reads: Sequence[Read], *, y: float = 0.0, h: float = 1.0, role: Optional[str] = None,
                 region=None, xform: XForm = None, nbins: Optional[int] = None, zorder: int = 3):
    """A group of reads as one compressed raster strip from ``y`` to ``y + h``.

    Per column: body alpha = 0.2 + 0.8 * (bodies / n) in the body colour (``subtle`` or the
    sample's layer-A ``role``); the poly(A) prefix wins the column in ``polya`` when tails are
    >= a quarter of the bodies present (or bodies have ended); the non-A clip remainder in
    ``mute`` by the same rule. Returns the ``AxesImage``.
    """
    from matplotlib.colors import to_rgb
    r = _apply_xlim(ax, region)
    x0, x1 = (r.start, r.end) if r is not None else sorted(ax.get_xlim())
    n = len(reads)
    if n == 0:
        return None
    nb = nbins or _nbins_for(ax, x0, x1)
    cb, ct, cc = _columns(reads, x0, x1, nb, xform)
    G = TOK.track_geometry()
    body_rgb = np.array(to_rgb(_body_color(role)))
    tail_rgb = np.array(to_rgb(TOK.color("polya")))
    clip_rgb = np.array(to_rgb(TOK.color("mute")))
    img = np.zeros((1, nb, 4))
    tail_wins = (ct > 0) & (ct * 4 >= cb)
    clip_wins = (cc > 0) & (cc * 4 >= cb) & ~tail_wins
    body = (cb > 0) & ~tail_wins & ~clip_wins
    a0 = G.get("merged_alpha_floor", 0.2)
    img[0, tail_wins, :3] = tail_rgb; img[0, tail_wins, 3] = np.minimum(1, 0.35 + 0.65 * ct[tail_wins] / n)
    img[0, clip_wins, :3] = clip_rgb; img[0, clip_wins, 3] = np.minimum(1, 0.25 + 0.55 * cc[clip_wins] / n)
    img[0, body, :3] = body_rgb; img[0, body, 3] = np.minimum(1, a0 + (1 - a0) * cb[body] / n)
    X = xform or (lambda p: p)
    lo, hi = sorted((X(x0), X(x1)))
    im = ax.imshow(img, extent=(lo, hi, y, y + h), aspect="auto", interpolation="nearest",
                   origin="lower", zorder=zorder)
    if r is not None:
        ax.set_xlim(*r.xlim())
    return im


# ============================================================================
# ribbon (wedge)
# ============================================================================
def ribbon(ax, reads: Sequence[Read], *, y: float = 0.0, h: float = 1.0, anchor: Optional[int] = None,
           role: Optional[str] = None, minor: bool = False, capped: bool = False, region=None,
           xform: XForm = None, nbins: Optional[int] = None, letter: Optional[str] = None,
           zorder: int = 3) -> list:
    """One cluster of reads as a WEDGE: at each x the stacked heights are the fractions of
    the cluster's reads whose body / poly(A) prefix / non-A clip cover it, scaled to ``h``.
    ``anchor`` (the shared 3' end) draws a sharp ink line through the band; ``letter`` is the
    isoform identity badge at the anchor; ``minor`` dims the body; ``capped`` marks a
    proportionally shaved band with a dashed rule over its extent.
    """
    r = _apply_xlim(ax, region)
    x0, x1 = (r.start, r.end) if r is not None else sorted(ax.get_xlim())
    n = len(reads)
    out: list = []
    if n == 0:
        return out
    nb = nbins or _nbins_for(ax, x0, x1)
    cb, ct, cc = _columns(reads, x0, x1, nb, xform)
    X = xform or (lambda p: p)
    lo, hi = sorted((X(x0), X(x1)))
    edges = np.linspace(lo, hi, nb + 1)
    S = TOK.stroke()
    body_col = _body_color(role)
    if minor:
        from matplotlib.colors import to_rgb
        rgb = np.array(to_rgb(body_col)); body_col = tuple(rgb + (1 - rgb) * 0.42)
    hb = cb / n * h; ht = ct / n * h; hc = cc / n * h
    base = np.full(nb, y)
    out.append(ax.fill_between(edges[:-1], base, base + hb, step="post", color=body_col, lw=0, zorder=zorder))
    top_b = base + hb
    out.append(ax.fill_between(edges[:-1], top_b, top_b + ht, step="post", color=TOK.color("polya"), lw=0,
                               zorder=zorder))
    # the body/tail junction rule: a 1-px ink line wherever a tail sits on a body
    has = ht > 0
    if has.any():
        yy = np.where(has, top_b, np.nan)
        out.append(ax.step(edges[:-1], yy, where="post", color=TOK.color("ink"), lw=S["hairline"], alpha=0.5,
                           zorder=zorder + 1)[0])
    top_t = top_b + ht
    out.append(ax.fill_between(edges[:-1], top_t, top_t + hc, step="post", color=TOK.color("mute"), lw=0,
                               alpha=0.7, zorder=zorder))
    if anchor is not None:
        ax_x = X(anchor)
        out.append(ax.plot([ax_x, ax_x], [y - 0.02 * h, y + 1.02 * h], color=TOK.color("ink"), lw=S["hairline"],
                           alpha=0.65, zorder=zorder + 2)[0])
        if letter:
            out.append(ax.annotate(letter, xy=(ax_x, y + h), xytext=(0, 2), textcoords="offset points",
                                   ha="center", va="bottom", fontsize=TOK.typography()["annotation"],
                                   fontweight="bold", color=TOK.color("ink"), zorder=zorder + 3))
    if capped:
        a = min(min(X(rd.start), X(rd.end)) for rd in reads); b = max(max(X(rd.start), X(rd.end)) for rd in reads)
        out.append(ax.plot([a, b], [y + h, y + h], color=TOK.color("ink"), lw=S["hairline"], alpha=0.4,
                           ls=(0, (3, 3)), zorder=zorder + 2)[0])
    if r is not None:
        ax.set_xlim(*r.xlim())
    return out


# ============================================================================
# the union-exon (spliced) axis
# ============================================================================
@dataclass
class Segment:
    s: int
    e: int
    cov: int = 0
    exonic: bool = True
    intronic: bool = False
    keep: bool = False
    compressed: bool = False
    wbp: float = 0.0
    t0: float = 0.0
    t1: float = 0.0


class SplicedAxis:
    """Union-exon space for a window (rbrowse ``buildAxis``, the compressed level).

    ``scope`` = (start, end); ``transcripts`` = the annotated models in the window (their
    exon boundaries are cut points; a segment is exonic when it lies inside an exon of
    the FIRST transcript listed per gene -- pass the canonical one first). Kept segments
    draw at 1 bp = 1 unit; collapsed ones take ``gap_bp`` units; kept non-exonic segments
    longer than ``seg_cap_bp`` are compressed to that width (``ruler`` marks both).
    ``strict=True`` keeps only covered segments with no gaps (the inspector's transcript
    space); ``exon_only=True`` keeps annotated exons only and removes introns entirely.
    """

    def __init__(self, scope: Interval, reads: Sequence[Read], transcripts: Sequence[Transcript] = (),
                 *, strict: bool = False, exon_only: bool = False, gap_bp: Optional[float] = None,
                 seg_cap_bp: Optional[int] = None):
        G = TOK.track_geometry()
        self.scope = (int(scope[0]), int(scope[1]))
        self.gap_bp = G.get("spliced_gap_bp", 60) if gap_bp is None else gap_bp
        self.seg_cap_bp = G.get("spliced_seg_cap_bp", 500) if seg_cap_bp is None else seg_cap_bp
        cuts = {self.scope[0], self.scope[1]}
        exons: List[Interval] = []
        introns: List[Interval] = []
        for tx in transcripts:
            for s, e in tx.exons:
                cuts.add(s); cuts.add(e); exons.append((s, e))
            for s, e in tx.introns:
                introns.append((s, e))
        for r in reads:
            for s, l in r.junctions or []:
                cuts.add(s); cuts.add(s + l)
        pts = sorted(p for p in cuts if self.scope[0] <= p <= self.scope[1])
        segs = [Segment(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
        for sg in segs:
            sg.exonic = any(sg.s >= s and sg.e <= e for s, e in exons) if exons else \
                not any(sg.s >= s and sg.e <= e for s, e in introns)
            sg.intronic = (not sg.exonic) and any(sg.s >= tx.start and sg.e <= tx.end for tx in transcripts)
        for r in reads:
            seen = set()
            for bs, bl in r.blocks:
                for i, sg in enumerate(segs):
                    if sg.s < bs + bl and sg.e > bs:
                        seen.add(i)
            for i in seen:
                segs[i].cov += 1
        t = 0.0
        gap = 0.0 if strict else self.gap_bp
        for sg in segs:
            sg.keep = (sg.exonic if exon_only else (sg.cov > 0 if strict else (sg.cov > 0 or sg.exonic)))
            ln = sg.e - sg.s
            if not sg.keep:
                gu = 0.0 if (exon_only and sg.intronic) else gap
                sg.wbp = 0.0; sg.t0 = t; t += gu; sg.t1 = t; sg.compressed = False
                continue
            sg.compressed = (not sg.exonic) and ln > self.seg_cap_bp
            sg.wbp = float(self.seg_cap_bp if sg.compressed else ln)
            sg.t0 = t; t += sg.wbp; sg.t1 = t
        self.segs = segs
        self.t_end = t
        self._starts = [sg.s for sg in segs]

    def _seg_at(self, pos: float) -> int:
        import bisect
        return max(0, bisect.bisect_right(self._starts, pos) - 1)

    def to_t(self, pos: float) -> float:
        """Genome -> display coordinate (piecewise linear, monotone)."""
        s0, s1 = self.scope
        if not self.segs:
            return pos - s0
        if pos <= s0:
            return float(pos - s0)
        if pos >= s1:
            return self.t_end + (pos - s1)
        sg = self.segs[self._seg_at(pos)]
        f = (pos - sg.s) / max(1, sg.e - sg.s)
        return sg.t0 + f * (sg.t1 - sg.t0)

    def from_t(self, tv: float) -> float:
        s0, s1 = self.scope
        if not self.segs:
            return s0 + tv
        if tv <= 0:
            return s0 + tv
        if tv >= self.t_end:
            return s1 + (tv - self.t_end)
        import bisect
        t0s = [sg.t0 for sg in self.segs]
        i = max(0, bisect.bisect_right(t0s, tv) - 1)
        sg = self.segs[i]
        f = (tv - sg.t0) / max(1e-9, sg.t1 - sg.t0)
        return sg.s + f * (sg.e - sg.s)

    @property
    def xlim(self) -> Tuple[float, float]:
        return (0.0, self.t_end)

    def ruler(self, ax, *, y: float = 0.0, height: float = 0.25, label_kb: bool = True):
        """A segment strip under the axis: exonic segments as ink bars, kept non-exonic as
        hairline bars, a '~' at every collapse, a break mark on every compressed segment,
        and the genome coordinate at each kept segment's start. Returns the artists."""
        S = TOK.stroke(); T = TOK.typography()
        out = []
        for sg in self.segs:
            if not sg.keep:
                if sg.t1 > sg.t0:
                    out.append(ax.text((sg.t0 + sg.t1) / 2, y + height / 2, "∼", ha="center", va="center",
                                       fontsize=T["annotation"], color=TOK.color("subtle"),
                                       fontfamily=[TOK.load()["typography"]["family"][0], "DejaVu Sans"]))
                continue
            col = TOK.color("ink") if sg.exonic else TOK.color("hairline")
            hh = height if sg.exonic else height * 0.5
            out.append(ax.add_patch(Rectangle((sg.t0, y + (height - hh) / 2), sg.t1 - sg.t0, hh,
                                              facecolor=col, edgecolor="none")))
            if sg.compressed:
                out.append(ax.plot([sg.t0, sg.t1], [y - 0.02, y - 0.02], color=TOK.color("polya"), lw=S["hairline"],
                                   ls=(0, (2, 1.5)))[0])
            if label_kb:
                out.append(ax.annotate(f"{sg.s / 1000:,.1f}", xy=(sg.t0, y), xytext=(0, -2), textcoords="offset points",
                                       ha="left", va="top", fontsize=T["tick_label"], color=TOK.color("subtle")))
        ax.set_xlim(*self.xlim)
        return out


# ============================================================================
# isoform chains
# ============================================================================
MIN_PARENT_N = 3
MIN_PARENT_FRAC = 0.5
NEAR_DUP_NT = 3


@dataclass
class Chain:
    key: str
    junctions: List[Interval]
    ret: List[int]
    reads: List[Read] = field(default_factory=list)
    n: int = 0
    parents: List[str] = field(default_factory=list)
    maximal: bool = False
    shared: bool = False
    assigned_to: Optional[str] = None
    cls: str = "fl"
    glyph: str = ""
    label: str = ""
    tx_name: Optional[str] = None
    dominant: bool = False
    support: int = 0
    near_dup_of: Optional[str] = None

    @property
    def jset(self) -> set:
        return {f"{s}-{s + l}" for s, l in self.junctions}


def _scoped_junctions(r: Read, scope: Interval) -> List[Interval]:
    return sorted((s, l) for s, l in (r.junctions or []) if s >= scope[0] and s + l <= scope[1])


def chain_clusters(reads: Sequence[Read], scope: Interval, model: Optional[Transcript] = None,
                   *, exon_names: Optional[Sequence[str]] = None,
                   annotated: Sequence[Transcript] = ()) -> Tuple[List[Chain], Dict[str, int]]:
    """Cluster reads by exact junction chain within ``scope`` and name each cluster relative
    to ``model`` (the canonical transcript). Returns ``(chains, summary)``; chains are ordered
    dominant first, then maximal, assigned, retention, shared, unspliced.

    Port of rbrowse ``clusterByChain`` (SPLICED_VIEW_DESIGN.md section 2): retention = an
    aligned block running CONTIGUOUSLY through a canonical intron boundary the read does not
    splice at; DRS reads are 3'-anchored suffixes, so a chain whose junctions are a subset of
    exactly one maximal chain is ASSIGNED to it (compatible), a chain extending >= 2 maximal
    chains stays its own SHARED cluster; unspliced reads are uninformative, never shared.
    """
    ref_introns = list(model.introns) if model else []
    exons = list(model.exons) if model else []
    minus = bool(model and model.strand == "-")
    n_ex = len(exons)

    def exon_name(i):
        k = (n_ex - i) if minus else (i + 1)
        return (exon_names[k - 1] if exon_names and k - 1 < len(exon_names) else str(k))

    def intron_name(i):
        k = (n_ex - 1 - i) if minus else (i + 1)
        return (exon_names[k - 1] if exon_names and k - 1 < len(exon_names) else str(k))

    annot_j = {f"{s}-{e}" for s, e in ref_introns}
    tx_by_chain: Dict[str, str] = {}
    for tx in annotated:
        chain = "|".join(f"{s}-{e}" for s, e in tx.introns)
        annot_j |= {f"{s}-{e}" for s, e in tx.introns}
        if chain:
            tx_by_chain[chain] = tx.label or tx.name

    clusters: Dict[str, Chain] = {}
    for r in reads:
        js = _scoped_junctions(r, scope)
        jk = [f"{s}-{s + l}" for s, l in js]
        ret = [i for i, (s, e) in enumerate(ref_introns)
               if f"{s}-{e}" not in jk and (r.spans_point(s) or r.spans_point(e))]
        key = ("|".join(jk) if jk else "unspliced") + ("||ret:" + ",".join(map(str, ret)) if ret else "")
        c = clusters.setdefault(key, Chain(key, js, ret))
        c.reads.append(r)
    lst = list(clusters.values())
    for c in lst:
        c.n = len(c.reads)

    def compatible(A: Chain, B: Chain) -> bool:
        if A is B or A.ret or B.ret:
            return False
        if len(A.jset) >= len(B.jset):
            return False
        if B.n < max(MIN_PARENT_N, MIN_PARENT_FRAC * A.n):
            return False
        if not A.jset <= B.jset:
            return False
        for s, l in B.junctions:
            if f"{s}-{s + l}" in A.jset:
                continue
            if any(r.spans_point(s) or r.spans_point(s + l) for r in A.reads):
                return False
        return True

    for c in lst:
        c.parents = [o.key for o in lst if compatible(c, o)]
    for c in lst:
        c.maximal = bool(c.junctions) and not c.ret and not c.parents
    max_keys = {c.key for c in lst if c.maximal}
    for c in lst:
        mp = [k for k in c.parents if k in max_keys]
        informative = bool(c.junctions)
        c.shared = informative and len(mp) >= 2
        c.assigned_to = mp[0] if (informative and len(mp) == 1) else None

    canon_j = {f"{s}-{e}" for s, e in ref_introns}
    for c in lst:
        parts: List[str] = []; cls = "fl"; glyph = ""
        if not c.junctions and not c.ret:
            cls, glyph = "unspliced", "~"; parts.append("unspliced")
        else:
            skips, novel, alt = [], [], []
            n_alt = 0
            for s, l in c.junctions:
                js, je, k = s, s + l, f"{s}-{s + l}"
                for i, (a, b) in enumerate(exons):
                    if a >= js and b <= je:
                        skips.append(i)
                if k in canon_j:
                    continue
                skip_only = any(a >= js and b <= je for a, b in exons) and \
                    any(x == js for x, _ in ref_introns) and any(y == je for _, y in ref_introns)
                if not skip_only:
                    d_shared = any(x == js for x, _ in ref_introns); a_shared = any(y == je for _, y in ref_introns)
                    plus = not minus
                    if d_shared and not a_shared:
                        alt.append("alt3′SS" if plus else "alt5′SS")
                    elif a_shared and not d_shared:
                        alt.append("alt5′SS" if plus else "alt3′SS")
                    elif not d_shared and not a_shared:
                        n_alt += 1
                if model is not None and k not in annot_j:
                    novel.append(k)
            if n_alt:
                alt.append("alt junction" + (f"×{n_alt}" if n_alt > 1 else ""))
            if c.ret:
                cls, glyph = "ret", "▮"; parts.append("ret-i" + ",".join(intron_name(i) for i in c.ret))
            if skips:
                cls = "skip" if cls == "fl" else cls; glyph = glyph or "⊘"
                parts.append("Δe" + ",".join(exon_name(i) for i in sorted(set(skips))))
            if alt:
                cls = "alt" if cls == "fl" else cls; glyph = glyph or "⇢"
                parts.append("+".join(dict.fromkeys(alt)))
            if novel:
                cls = "novel" if cls == "fl" else cls; glyph = "✱"
                parts.append("novel" + (f"×{len(novel)}" if len(novel) > 1 else ""))
            if not parts:
                parts.append("FL" if c.maximal else "FL-compatible")
        if c.shared:
            glyph = "≈"; parts.insert(0, "shared")
        ck = "|".join(sorted(c.jset, key=lambda k: int(k.split("-")[0])))
        if ck in tx_by_chain:
            c.tx_name = tx_by_chain[ck]
        c.cls, c.glyph, c.label = cls, glyph, " ".join(parts)

    for c in lst:
        if not c.junctions:
            continue
        for o in lst:
            if o is c or o.n <= c.n or len(o.junctions) != len(c.junctions):
                continue
            diff, near = 0, True
            for (s1, l1), (s2, l2) in zip(c.junctions, o.junctions):
                if (s1, l1) == (s2, l2):
                    continue
                diff += 1
                if abs(s1 - s2) > NEAR_DUP_NT or abs((s1 + l1) - (s2 + l2)) > NEAR_DUP_NT:
                    near = False
            if diff == 1 and near:
                c.near_dup_of = o.key; c.label += " ≈dup"
                break

    for c in lst:
        c.support = c.n + sum(o.n for o in lst if o.assigned_to == c.key)
    dominant = None
    for c in lst:
        if c.maximal and (dominant is None or c.support > dominant.support):
            dominant = c
    for c in lst:
        c.dominant = c is dominant

    def rank(c: Chain) -> int:
        return 0 if c.dominant else 1 if c.maximal else 2 if c.assigned_to else 3 if c.ret else 4 if c.shared else 5

    lst.sort(key=lambda c: (rank(c), -c.support))
    summary = {
        "n_reads": len(reads),
        "unambiguous": sum(c.n for c in lst if c.maximal),
        "compatible": sum(c.n for c in lst if c.assigned_to),
        "shared": sum(c.n for c in lst if c.shared),
        "unspliced": sum(c.n for c in lst if not c.junctions and not c.ret),
        "retention": sum(c.n for c in lst if c.ret),
    }
    return lst, summary


def isoform_rows(ax, chains: Sequence[Chain], axis: SplicedAxis, *, top: int = 6, y0: float = 0.0,
                 pitch: float = 1.0, height: float = 0.5, role: Optional[str] = None,
                 letters: str = "ABCDEFGHIJ", denominator: Optional[int] = None,
                 chip_space: float = 0.42) -> Dict[str, float]:
    """The top-k chains as rows on the spliced axis: covered kept segments as blocks (the
    body colour; hatched when shared), a hairline connector across the chain's extent, a
    LETTER badge (the identity -- never a hue), and a monochrome chip: glyph, label,
    transcript name, n and the share of the denominator. ``chip_space`` reserves that
    fraction of the axis width to the right of the rows for the chips (the x-limits become
    ``(0, t_end * (1 + chip_space))``; draw a companion axis with the same limits). Returns
    ``{chain key: y}``.
    """
    T = TOK.typography(); S = TOK.stroke()
    kept = [sg for sg in axis.segs if sg.keep]
    ys: Dict[str, float] = {}
    rows = [c for c in chains if c.maximal or c.ret or c.shared or c.assigned_to][:top]
    ink = TOK.color("ink")
    body = _body_color(role)
    for k, c in enumerate(rows):
        y = y0 - k * pitch
        cov = set()
        for r in c.reads:
            for i, sg in enumerate(kept):
                if r.has_bases_in(sg.s, sg.e):
                    cov.add(i)
        lo = min(min(r.start, r.end) for r in c.reads); hi = max(max(r.start, r.end) for r in c.reads)
        ax.plot([axis.to_t(lo), axis.to_t(hi)], [y, y], color=TOK.color("hairline"), lw=S["hairline"], zorder=2)
        for i in cov:
            sg = kept[i]
            ax.add_patch(Rectangle((sg.t0, y - height / 2), sg.t1 - sg.t0, height, facecolor=body,
                                   edgecolor=ink if c.shared else "none", lw=S["hairline"] if c.shared else 0,
                                   hatch="////" if c.shared else None, zorder=3))
        letter = letters[k] if k < len(letters) else str(k)
        share = f" · {100 * c.n / denominator:.0f}%" if denominator else ""
        # the class glyphs (★ ⊘ ▮ ⇢ ✱ ≈) are not in Arial/Helvetica: they go in their own Text
        # with an explicit DejaVu fallback so they never render as tofu, whatever the rc family
        glyph = ("★" if c.dominant else "") + (c.glyph if c.glyph else "")
        chip = f"{c.label}" + (f" ({c.tx_name})" if c.tx_name else "") + f" · n={c.n}{share}"
        # the letter sits left of the first COVERED segment (a block can start before the reads)
        x_first = min((kept[i].t0 for i in cov), default=axis.to_t(lo))
        ax.annotate(letter, xy=(x_first, y), xytext=(-6, 0), textcoords="offset points", ha="right",
                    va="center", fontsize=T["in_figure"], fontweight="bold", color=ink)
        gl = ax.annotate(glyph, xy=(axis.to_t(hi), y), xytext=(4, 0), textcoords="offset points", ha="left",
                         va="center", fontsize=T["annotation"], color=ink,
                         fontfamily=[TOK.load()["typography"]["family"][0], "DejaVu Sans"])
        ax.annotate(chip, xy=(axis.to_t(hi), y), xytext=(4 + (len(glyph) * T["annotation"] * 0.95), 0),
                    textcoords="offset points", ha="left", va="center", fontsize=T["annotation"], color=ink)
        ys[c.key] = y
    ax.set_xlim(0.0, axis.t_end * (1.0 + chip_space))
    ax.set_ylim(y0 - (len(rows) - 1) * pitch - 0.8, y0 + 0.8)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    return ys


# ============================================================================
# junction arcs (sashimi)
# ============================================================================
def junction_counts(reads: Sequence[Read], scope: Optional[Interval] = None) -> Dict[Interval, int]:
    """``{(intron_start, intron_end): n_reads}`` within ``scope``."""
    out: Dict[Interval, int] = {}
    for r in reads:
        for s, l in r.junctions or []:
            if scope is None or (s >= scope[0] and s + l <= scope[1]):
                out[(s, s + l)] = out.get((s, s + l), 0) + 1
    return out


def junction_arcs(ax, counts: Dict[Interval, int], *, y: float = 0.0, height: Optional[float] = None,
                  role: Optional[str] = None, min_support: int = 1, lw_range: Tuple[float, float] = (0.6, 2.6),
                  label: bool = True, annotated: Iterable[Interval] = (), xform: XForm = None,
                  direction: str = "up") -> list:
    """Sashimi arcs: one per junction with >= ``min_support`` reads, width scaled between
    ``lw_range`` by support, in the sample's layer-A ``role`` (signal) or ``subtle``; junctions
    in ``annotated`` that no read supports are drawn as thin ``splice`` annotation arcs.
    ``height`` is the apex above ``y`` in data units (default: the token ratio of the y-range).
    """
    X = xform or (lambda p: p)
    out = []
    if not counts and not annotated:
        return out
    mx = max(counts.values()) if counts else 1
    for (s, e), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        if n < min_support:
            continue
        lw = lw_range[0] + (lw_range[1] - lw_range[0]) * (n / mx)
        out += arc(ax, X(s), X(e), y=y, height=height, role=role or "neutral", lw=lw,
                   label=str(n) if label else None, direction=direction)
    for s, e in annotated:
        if (s, e) not in counts:
            out += arc(ax, X(s), X(e), y=y, height=height, lw=TOK.stroke()["hairline"], direction=direction,
                       alpha=0.8)
    return out
