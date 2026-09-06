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
3'-end anchor (the site mark, the ribbon's shared edge, the chain badge) is the
``polya`` identity, and the poly(A) TAIL on a read is the ``tail`` identity -- the
green rbrowse paints on screen, searched for print (tokens 1.4.0; Kevin's ruling R9,
2026-09-06: a green tail at the 3' end of every read is the reader's marker of
orientation and of a full RNA head to tail); the soft-clip remainder is
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
    # space economy (Chanfreau planning/880)
    "Cluster", "BandPlan", "PanelPlan", "BandScale", "cluster_by_three_prime", "union_clusters",
    "solve_band_scale", "select_bands", "effective_counts", "plan_read_panel",
    "read_panel", "read_stack", "stack_height_mm",
    # the dense pile, quantification windows, truncation survival (Chanfreau planning/885)
    "sort_reads", "dense_pile", "Window", "windows_from_clusters", "window_bands", "window_counts",
    "rep_values", "Survival", "spliced_distance", "corrected_coverage", "five_prime_profile",
    # the compact block (planning/885; the 887 landing)
    "pile_stack", "strip_rows", "halflife_ladder",
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
        span = int(self.end) - int(self.start)
        if any(l > span or l <= 0 for _, l in self.blocks):
            raise ValueError(f"Read {self.id!r}: blocks are (start, LENGTH) -- rbrowse's wire format -- "
                             f"not (start, end); a block of length {max(l for _, l in self.blocks)} cannot "
                             f"fit a read spanning {span} bp (audit 879, F3.8)")
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


def coverage_columns(reads: Sequence[Read], x0: float, x1: float, nbins: int, xform: XForm = None):
    """Per-bin (body, tail, clip) read counts across [x0, x1): the aligned bodies, the poly(A)
    tails that hang PAST each read's 3' end, and the non-A clip remainder. The public form of
    the raster's column accumulator, for a coverage panel that wants to draw the tails as a
    ``tail`` fill beyond the 3' edge (R9; Chanfreau planning/889a)."""
    return _columns(reads, x0, x1, nbins, xform)


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
                 region=None, xform: XForm = None, nbins: Optional[int] = None, zorder: int = 3,
                 survival: Optional["Survival"] = None):
    """A group of reads as one compressed raster strip from ``y`` to ``y + h``.

    Per column: body alpha = 0.2 + 0.8 * (bodies / n) in the body colour (``subtle`` or the
    sample's layer-A ``role``); the poly(A) prefix wins the column in ``polya`` when tails are
    >= a quarter of the bodies present (or bodies have ended); the non-A clip remainder in
    ``mute`` by the same rule. With ``survival`` (the library's own truncation half-life) the
    body fraction is the survival-CORRECTED coverage over n, capped at 1: the strip then shows
    the molecules the library would have covered, not the reads it did -- a model, which the
    caller must declare (chip / legend). Returns the ``AxesImage``.
    """
    from matplotlib.colors import to_rgb
    r = _apply_xlim(ax, region)
    x0, x1 = (r.start, r.end) if r is not None else sorted(ax.get_xlim())
    n = len(reads)
    if n == 0:
        return None
    nb = nbins or _nbins_for(ax, x0, x1)
    cb, ct, cc = _columns(reads, x0, x1, nb, xform)
    if survival is not None:
        cb = np.minimum(float(n), _weighted_columns(reads, x0, x1, nb, survival, xform)[0])
    G = TOK.track_geometry()
    body_rgb = np.array(to_rgb(_body_color(role)))
    tail_rgb = np.array(to_rgb(TOK.color("tail")))
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
    out.append(ax.fill_between(edges[:-1], top_b, top_b + ht, step="post", color=TOK.color("tail"), lw=0,
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
    hidden_n: int = 0        # exon_only: reads with bases in this REMOVED segment
    across: int = 0          # exon_only: reads spanning it with no bases
    hidden_frac: float = 0.0


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
                # only junctions INSIDE the scope cut the axis -- the same rule chain_clusters
                # applies (`_scoped_junctions`). A cDNA alignment can carry a spurious intron
                # from the 3' end to tens of kb away; its start cut the 3' flank into 14
                # segments inside 100 bp at RPL24B while the chains ignored it (planning/881d)
                if s >= self.scope[0] and s + l <= self.scope[1]:
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
        n_reads = len(reads) or 1
        # 🔴 THE 3' FLANK STAYS 1:1 WHERE READS END THERE (rbrowse spliced.js:155-170, Kevin
        # 2026-08-16: "the read alignments look cut off"). Segments on a gene's 3' side that >= 1 %
        # of the reads cover keep their real width up to `downstream_keep_bp` in total; beyond that
        # the flank compresses like any non-exonic stretch. Without this every 3' end past the
        # annotation -- a longer-than-annotated UTR, the very landscape this lab works on -- piles
        # onto one compression boundary and reads as truncated. Missing from the first port
        # (audit 879, F7.3).
        self.downstream_keep_bp = G.get("spliced_downstream_keep_bp", 3000)
        keep3 = set()
        for tx in transcripts:
            plus = tx.strand != "-"
            end3 = tx.end if plus else tx.start
            budget = self.downstream_keep_bp
            for sg in (segs if plus else reversed(segs)):
                past = sg.s >= end3 if plus else sg.e <= end3
                if not past or sg.exonic:
                    continue
                if any(o is not tx and sg.s >= o.start and sg.e <= o.end for o in transcripts):
                    break                     # reached the next gene
                if sg.cov < 0.01 * n_reads:
                    break
                ln = sg.e - sg.s
                if ln > budget:
                    break
                keep3.add(id(sg)); budget -= ln
        t = 0.0
        gap = 0.0 if strict else self.gap_bp
        for sg in segs:
            k3 = id(sg) in keep3
            sg.keep = ((sg.exonic or k3) if exon_only else (sg.cov > 0 if strict else (sg.cov > 0 or sg.exonic)))
            ln = sg.e - sg.s
            if not sg.keep:
                gu = 0.0 if (exon_only and sg.intronic) else gap
                sg.wbp = 0.0; sg.t0 = t; t += gu; sg.t1 = t; sg.compressed = False
                continue
            sg.compressed = (not sg.exonic) and ln > self.seg_cap_bp and not k3
            sg.wbp = float(self.seg_cap_bp if sg.compressed else ln)
            sg.t0 = t; t += sg.wbp; sg.t1 = t
        # 🔴 NOTHING HIDDEN SILENTLY (spliced.js:205-225). For every REMOVED segment in the
        # exon-only level: how many reads have bases in it (retention / alt site / novel exon --
        # invisible once the axis drops it) against how many span it without bases. Exposed as
        # `hidden` and drawn by `ruler` as a faint tick. Missing from the first port (879, F7.4).
        self.hidden: List[Segment] = []
        if exon_only:
            for sg in segs:
                if sg.keep:
                    continue
                across = sum(1 for r in reads if r.start < sg.s and r.end > sg.e and not r.has_bases_in(sg.s, sg.e))
                sg.hidden_n = sg.cov; sg.across = across
                sg.hidden_frac = (sg.cov / (sg.cov + across)) if across else 0.0
                if sg.cov:
                    self.hidden.append(sg)
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

    def ruler(self, ax, *, y: float = 0.0, height: float = 0.25, label_kb: bool = True,
              min_label_gap: Optional[float] = None):
        """A segment strip under the axis: exonic segments as ink bars, kept non-exonic as
        hairline bars, a '~' at every collapse, a break mark on every compressed segment,
        and the genome coordinate at each kept segment's start. Returns the artists.

        ``min_label_gap`` (display units; default 7 % of the axis) thins the coordinate
        labels: a segment whose start lies closer than that to the last labelled start gets
        no label. Deep real data cuts the 3' flank at every read end ≥ 1 % of the reads
        (RPL24B: 14 segments inside 100 bp), and a label per segment overprinted itself
        (planning/881d, render 2)."""
        S = TOK.stroke(); T = TOK.typography()
        out = []
        gap = (0.07 * self.t_end) if min_label_gap is None else min_label_gap
        last_label_t = None
        for sg in self.segs:
            if not sg.keep:
                if sg.t1 > sg.t0:
                    out.append(ax.text((sg.t0 + sg.t1) / 2, y + height / 2, "∼", ha="center", va="center",
                                       fontsize=T["annotation"], color=TOK.color("subtle"),
                                       fontfamily=[TOK.load()["typography"]["family"][0], "DejaVu Sans"]))
                if sg.hidden_n:
                    # a removed segment some reads had bases in: a faint tick and the rate, never silence
                    out.append(ax.plot([sg.t0, sg.t0], [y - height * 0.3, y + height * 1.3], color=TOK.color("polya"),
                                       lw=S["hairline"], alpha=0.6)[0])
                    out.append(ax.annotate(f"{100 * sg.hidden_frac:.0f}% ({sg.hidden_n})" if sg.across else f"({sg.hidden_n})",
                                           xy=(sg.t0, y + height * 1.3), xytext=(1, 1), textcoords="offset points",
                                           ha="left", va="bottom", fontsize=T["annotation"], color=TOK.color("polya")))
                continue
            col = TOK.color("ink") if sg.exonic else TOK.color("hairline")
            hh = height if sg.exonic else height * 0.5
            out.append(ax.add_patch(Rectangle((sg.t0, y + (height - hh) / 2), sg.t1 - sg.t0, hh,
                                              facecolor=col, edgecolor="none")))
            if sg.compressed:
                out.append(ax.plot([sg.t0, sg.t1], [y - 0.02, y - 0.02], color=TOK.color("polya"), lw=S["hairline"],
                                   ls=(0, (2, 1.5)))[0])
            if label_kb and (last_label_t is None or sg.t0 - last_label_t >= gap):
                out.append(ax.annotate(f"{sg.s / 1000:,.1f}", xy=(sg.t0, y), xytext=(0, -2), textcoords="offset points",
                                       ha="left", va="top", fontsize=T["tick_label"], color=TOK.color("subtle")))
                last_label_t = sg.t0
        ax.set_xlim(*self.xlim)
        return out


# ============================================================================
# isoform chains
# ============================================================================
# 🔴 THE SAME NUMBERS AS rbrowse spliced.js:24-25, AND THE INCIDENT THAT SET THEM (spliced.js:18-23):
# "at a long gene (SDHA: 15 exons, 1,823 reads) the 989-read all-canonical suffix chain was declared
# 'shared' among a handful of n<=7 novel chains that happen to extend it -- the rare chains hijacked
# the main isoform." A parent chain must carry >= max(MIN_PARENT_N, MIN_PARENT_FRAC * n) reads to
# claim a suffix. The first port shipped 3 / 0.5 with no comment (audit 879, F7.1): 2.5x stricter,
# so the browser and the figure called the same reads differently. Never change one side alone.
MIN_PARENT_N = 2
MIN_PARENT_FRAC = 0.2
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
                # an all-canonical chain: FL, saying where its reads START -- DRS reads are
                # 5'-truncated, so "FL" means canonical junctions over its span, not that exon 1
                # was reached (spliced.js; missing from the first port, audit 879 F7.2)
                frm = ""
                if exons and c.junctions:
                    first, last = c.junctions[0], c.junctions[-1]
                    if minus:
                        idx = next((i for i, (a, b) in enumerate(exons) if a == last[0] + last[1]), -1)
                    else:
                        idx = next((i for i, (a, b) in enumerate(exons) if b == first[0]), -1)
                    if idx >= 0:
                        ord_ = exon_name(idx)
                        if ord_ != exon_name(n_ex - 1 if minus else 0):
                            frm = f" from e{ord_}"
                parts.append(("FL" if c.maximal else "FL-compatible") + frm)
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
        # clip to the scope: a cDNA alignment can run tens of kb past it through a spurious
        # intron, and `to_t` extends linearly beyond the scope, which put the chip off the
        # axis (planning/881d, render 1: the FL row had no visible chip)
        lo, hi = max(lo, axis.scope[0]), min(hi, axis.scope[1])
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
        chip = f"{c.label}" + (f" ({c.tx_name})" if c.tx_name else "") + f" · n={c.n:,}{share}"
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
    if role is not None:
        TOK.role(role)                      # refuse a layer-B token even when nothing is drawn
    out = []
    if not counts and not annotated:
        return out
    mx = max(counts.values()) if counts else 1
    for (s, e), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        if n < min_support:
            continue
        lw = lw_range[0] + (lw_range[1] - lw_range[0]) * (n / mx)
        out += arc(ax, X(s), X(e), y=y, height=height, role=role or "neutral", lw=lw,
                   label=f"{n:,}" if label else None, direction=direction)
    for s, e in annotated:
        if (s, e) not in counts:
            out += arc(ax, X(s), X(e), y=y, height=height, lw=TOK.stroke()["hairline"], direction=direction,
                       alpha=0.8)
    return out


# ============================================================================
# SPACE ECONOMY -- a stack of sample read panels on an mm budget
#
# Kevin, 2026-09-05: "I like showing the raw reads as stacked alignments in
# merged/ribbon views but I think we need to consider space economy."
#
# The rules and their measurements are Chanfreau ``planning/880_read_view_space
# _economy.md``; the mechanisms are ported from rbrowse's figure path (cited by
# file:line below, no code copied):
#
#  * a panel has NO content height -- it gets a BUDGET and the drawing fits into
#    it (rbrowse.js:10298-10317, the adaptive `need` sum);
#  * band height is proportional to the cluster's read count on ONE scale shared
#    across every panel, solved by bisection (rbrowse.js:5169, :13805-13867,
#    `bandScale:'shared'` :305), with a hard floor so a real minor cluster is
#    never a hairline (`RIBBON_MIN_PX` :4588, here `band_floor_mm`);
#  * keep-top k with everything else pooled into one VISIBLE `other` band
#    (:205, :8771 -- "an orphan beyond every capped reach falls to 'other',
#    which is visible; a window silently spanning read-free space is not");
#  * a shaved band declares its cut with a dashed rule over its own extent
#    (:5229, :17577);
#  * ONE letter strip for the whole stack instead of a label per band
#    (:17255-17282, superseding per-band badges :17461) -- letters, never hues;
#  * merged vs squished chosen from what the budget can show at a legible pitch
#    (`DENSITY` :5706, `AUTO_DENSITY` :5742).
#
# EVERY height in this section is in MILLIMETRES ON THE PAGE. A panel's axes is
# given ``ylim = (0, budget_mm)``, so one data unit is one millimetre and the
# arithmetic below is the geometry, not a proxy for it.
# ============================================================================

MM_PER_IN = 25.4


@dataclass
class Cluster:
    """One 3'-end cluster inside one sample."""
    anchor: int
    reads: List[Read] = field(default_factory=list)
    letter: Optional[str] = None
    pooled: bool = False          # the `other` band

    @property
    def n(self) -> int:
        return len(self.reads)

    @property
    def span(self) -> Optional[Interval]:
        """(min, max) 3' end over the members -- the cluster's real extent. Clusters are
        single-linkage chains, so a member can sit further than ``win`` from the modal
        anchor (RPL24B, planning/881: two CPA sub-sites 60 nt apart chain into one cluster)."""
        if not self.reads:
            return None
        ends = [r.three_prime for r in self.reads]
        return (min(ends), max(ends))


@dataclass
class BandPlan:
    """One drawn band, in mm on the page."""
    letter: Optional[str]
    anchor: Optional[int]
    n: int
    h_mm: float
    reads: List[Read] = field(default_factory=list)
    pooled: bool = False
    capped: bool = False          # shaved to fit; draws the dashed cut mark
    minor: bool = False


@dataclass
class PanelPlan:
    """The arithmetic for ONE sample panel. Pure geometry -- no matplotlib."""
    name: str
    budget_mm: float
    merged_mm: float
    merged_gap_mm: float
    bands_mm: float               # the region the bands may use
    band_gap_mm: float
    bands: List[BandPlan] = field(default_factory=list)
    mode: str = "merged"          # "merged" | "squished"
    pitch_mm: float = 0.0         # squished only
    lw_pt: float = 0.0            # squished only
    n_rows: int = 0               # squished only
    n_reads: int = 0
    n_pooled: int = 0             # reads that fell into `other`
    k_used: int = 0
    shaved_mm: float = 0.0

    @property
    def used_mm(self) -> float:
        """Height the panel actually draws into: the merged row plus either the bands and
        their gaps, or (squished) the per-read rows."""
        if self.mode == "squished":
            content = self.n_rows * self.pitch_mm
        else:
            nb = len(self.bands)
            content = sum(b.h_mm for b in self.bands) + max(0, nb - 1) * self.band_gap_mm
        return content + (self.merged_mm + self.merged_gap_mm if self.merged_mm > 0 else 0.0)

    @property
    def slack_mm(self) -> float:
        return self.budget_mm - self.used_mm


@dataclass
class BandScale:
    """mm of band height per read -- ONE scale for every panel in a stack, so a
    band's height means the same thing in every sample (rbrowse `bandScale:
    'shared'`). ``floor_mm`` is the hard minimum a kept band may have."""
    mm_per_read: float
    floor_mm: float

    def height_mm(self, n: int) -> float:
        return 0.0 if n <= 0 else max(self.floor_mm, self.mm_per_read * n)


# ---------------------------------------------------------------------------
# clustering by 3' end
# ---------------------------------------------------------------------------
def cluster_by_three_prime(reads: Sequence[Read], *, win: int = 32) -> List[Cluster]:
    """Group reads whose 3' ends fall within ``win`` bp, largest cluster first.

    ``win`` defaults to rbrowse's shipped ``clusterWin`` (32 bp, rbrowse.js:505-515):
    wide enough that a quantification window is a rectangle rather than a hairline.
    The anchor is the MODAL 3' end of the cluster -- the one coordinate the members
    share -- not the mean, which would sit where no read ends.
    """
    if not reads:
        return []
    ends = sorted(((r.three_prime, r) for r in reads), key=lambda t: t[0])
    groups: List[List[Read]] = []
    cur = [ends[0][1]]
    last = ends[0][0]
    for pos, r in ends[1:]:
        if pos - last <= win:
            cur.append(r)
        else:
            groups.append(cur); cur = [r]
        last = pos
    groups.append(cur)
    out = []
    for g in groups:
        counts: Dict[int, int] = {}
        for r in g:
            counts[r.three_prime] = counts.get(r.three_prime, 0) + 1
        anchor = max(counts.items(), key=lambda kv: (kv[1], -kv[0]))[0]
        out.append(Cluster(anchor=anchor, reads=list(g)))
    out.sort(key=lambda c: (-c.n, c.anchor))
    return out


def union_clusters(samples: Sequence[dict], *, win: int = 32, k: int = 4,
                   letters: str = "ABCDEFGHIJ") -> List[Cluster]:
    """The cross-panel cluster ranking: pool every sample's reads, cluster once, keep
    the top ``k`` by POOLED count and letter them A, B, C... in that order.

    One ranking for the whole stack is what makes a letter mean the same 3' end in
    every panel -- rbrowse's ``assignCrossPanelRanks`` (rbrowse.js:13803) and the
    reason its letter strip can be a single row (:17255). Deriving the letters twice
    is how two panels drift apart.
    """
    pooled: List[Read] = []
    for s in samples:
        pooled.extend(s.get("reads", ()))
    cl = cluster_by_three_prime(pooled, win=win)[:max(0, k)]
    for i, c in enumerate(cl):
        c.letter = letters[i] if i < len(letters) else str(i)
    return cl


def _panel_bands(reads: Sequence[Read], keep: Sequence[Cluster], *, win: int) -> Tuple[List[Cluster], Cluster]:
    """Split one sample's reads over the stack's kept clusters plus one `other`."""
    assigned: List[Cluster] = [Cluster(anchor=c.anchor, letter=c.letter) for c in keep]
    other = Cluster(anchor=None, letter=None, pooled=True)
    spans = [c.span for c in keep]
    for r in reads:
        tp = r.three_prime
        best, bestd = None, None
        for slot, c, sp in zip(assigned, keep, spans):
            d = abs(tp - c.anchor)
            # a read belongs to a kept cluster when its 3' end lies inside the cluster's
            # EXTENT (the chain it was pooled into) -- not only within `win` of the modal
            # anchor. Membership by anchor +- win disagreed with single-linkage clustering
            # and pooled 50 % of RPL24B's reads into `other` (planning/881d, render 1).
            inside = sp is not None and sp[0] <= tp <= sp[1]
            if (inside or d <= win) and (bestd is None or d < bestd):
                best, bestd = slot, d
        (best or other).reads.append(r)
    return assigned, other


# ---------------------------------------------------------------------------
# the budget arithmetic -- this is what the tests pin
# ---------------------------------------------------------------------------
def _band_room(nb: int, bands_mm: float, band_gap_mm: float) -> float:
    """mm available to the band FILLS once the inter-band gaps are paid."""
    return bands_mm - max(0, nb - 1) * band_gap_mm


def _fits(counts: Sequence[int], scale: BandScale, bands_mm: float, band_gap_mm: float) -> bool:
    live = [n for n in counts if n > 0]
    if not live:
        return True
    return sum(scale.height_mm(n) for n in live) <= _band_room(len(live), bands_mm, band_gap_mm) + 1e-9


def _reduce_k(counts: Sequence[int], *, floor_mm: float, bands_mm: float, band_gap_mm: float,
              pooled_n: int) -> int:
    """The largest number of bands whose FLOORS fit the budget. Everything dropped rolls
    into `other`, which stays visible -- rbrowse's rule that a pooled remainder is drawn
    and a silently-absent one is not (rbrowse.js:8771)."""
    live = sum(1 for n in counts if n > 0) + (1 if pooled_n > 0 else 0)
    while live > 0 and live * floor_mm > _band_room(live, bands_mm, band_gap_mm) + 1e-9:
        live -= 1
    return live


def solve_band_scale(panel_counts: Sequence[Sequence[int]], bands_mm: float, *,
                     floor_mm: float, band_gap_mm: float, cap_mm_per_read: float = 1.0,
                     iters: int = 40) -> BandScale:
    """The largest mm-per-read at which EVERY panel's bands fit ``bands_mm``.

    Bisection, like rbrowse's ~22-probe search (rbrowse.js:13840-13866), but without
    the warm start: a figure is rendered once, not at 60 fps. The upper bound is seeded
    from the largest single cluster so the search resolves the real answer instead of
    collapsing every panel to its floor -- the failure rbrowse's at-scale smoke found.
    """
    biggest = max((max(c) for c in panel_counts if c), default=1)
    biggest = max(1, biggest)
    hi = min(cap_mm_per_read, max(1e-6, 1.5 * bands_mm / biggest))
    lo = 0.0

    def ok(s: float) -> bool:
        sc = BandScale(s, floor_mm)
        return all(_fits(c, sc, bands_mm, band_gap_mm) for c in panel_counts)

    if ok(hi):
        return BandScale(hi, floor_mm)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return BandScale(lo, floor_mm)


def select_bands(reads: Sequence[Read], keep: Sequence[Cluster], *, bands_mm: float,
                 floor_mm: float, band_gap_mm: float, win: int = 32
                 ) -> Tuple[List[Cluster], Optional[Cluster]]:
    """Which bands a panel will actually draw, and what falls into `other`.

    Depends only on the FLOOR and the budget -- never on the band scale -- so the scale
    solver and the planner agree by construction. (They did not in the first version: the
    solver ignored the pooled band, the planner drew it, and the panel came out shaved.)
    """
    assigned, other = _panel_bands(reads, keep, win=win)
    counts = [c.n for c in assigned]
    live = _reduce_k(counts, floor_mm=floor_mm, bands_mm=bands_mm, band_gap_mm=band_gap_mm,
                     pooled_n=other.n)
    order = [i for i in sorted(range(len(assigned)), key=lambda j: (-counts[j], j)) if counts[i] > 0]
    room_for_other = 1 if other.n > 0 else 0
    n_keep = max(0, live - room_for_other)
    keep_idx, drop_idx = order[:n_keep], order[n_keep:]
    for i in drop_idx:
        other.reads.extend(assigned[i].reads)
    kept = [assigned[i] for i in keep_idx]
    pooled = other if (other.n > 0 and live > len(kept)) else None
    return kept, pooled


def effective_counts(reads: Sequence[Read], keep: Sequence[Cluster], *, bands_mm: float,
                     floor_mm: float, band_gap_mm: float, win: int = 32) -> List[int]:
    """The band counts a panel will draw -- what :func:`solve_band_scale` must be fed."""
    kept, pooled = select_bands(reads, keep, bands_mm=bands_mm, floor_mm=floor_mm,
                                band_gap_mm=band_gap_mm, win=win)
    return [c.n for c in kept] + ([pooled.n] if pooled is not None else [])


def plan_read_panel(name: str, reads: Sequence[Read], keep: Sequence[Cluster], *,
                    budget_mm: float, scale: BandScale, merged_mm: float = 2.5,
                    merged_gap_mm: float = 0.8, band_gap_mm: float = 0.4,
                    win: int = 32, mode: str = "auto",
                    pitch_floor_mm: Optional[float] = None,
                    max_lw_pt: float = 4.0) -> PanelPlan:
    """Lay ONE sample out inside ``budget_mm``. Pure arithmetic; draws nothing.

    Invariants (pinned by ``tests/test_visualize_pileup.py``):

    1. ``plan.used_mm <= plan.budget_mm`` -- a panel never exceeds its budget;
    2. every band in ``plan.bands`` has ``h_mm >= scale.floor_mm`` -- a kept band is
       never below the floor (bands that cannot clear it are POOLED, not shrunk);
    3. in squished mode ``pitch_mm >= pitch_floor_mm`` and
       ``n_rows * pitch_mm <= bands_mm`` -- the drawn ink never leaves the strip.
    """
    G = TOK.track_geometry()
    if pitch_floor_mm is None:
        pitch_floor_mm = G.get("pitch_floor_mm", 0.423)
    merged_block = (merged_mm + merged_gap_mm) if merged_mm > 0 else 0.0
    bands_mm = max(0.0, budget_mm - merged_block)
    plan = PanelPlan(name=name, budget_mm=budget_mm, merged_mm=merged_mm,
                     merged_gap_mm=merged_gap_mm, bands_mm=bands_mm, band_gap_mm=band_gap_mm,
                     n_reads=len(reads))

    # -- mode: squished only when the budget can show every read at a legible pitch
    want_squished = mode == "squished"
    if mode == "auto" and reads:
        want_squished = len(reads) * pitch_floor_mm <= bands_mm
    if want_squished and reads and len(reads) * pitch_floor_mm <= bands_mm:
        plan.mode = "squished"
        plan.n_rows = len(reads)
        plan.pitch_mm = bands_mm / plan.n_rows
        pitch_pt = plan.pitch_mm / MM_PER_IN * 72.0
        floor_pt = TOK.load()["geometry"].get("stroke_min_pt", 0.6)
        plan.lw_pt = max(floor_pt, min(max_lw_pt, pitch_pt - floor_pt))
        return plan
    if want_squished and mode == "squished":
        # asked for per-read rows the budget cannot honour: say so by falling back
        plan.mode = "merged"

    # -- merged + bands
    kept_cl, pooled_band = select_bands(reads, keep, bands_mm=bands_mm, floor_mm=scale.floor_mm,
                                        band_gap_mm=band_gap_mm, win=win)
    plan.k_used = len(kept_cl)
    plan.n_pooled = pooled_band.n if pooled_band is not None else 0

    bands: List[BandPlan] = []
    for c in kept_cl:
        bands.append(BandPlan(letter=c.letter, anchor=c.anchor, n=c.n, reads=c.reads,
                              h_mm=scale.height_mm(c.n)))
    if pooled_band is not None:
        bands.append(BandPlan(letter=None, anchor=None, n=pooled_band.n, reads=pooled_band.reads,
                              h_mm=scale.height_mm(pooled_band.n), pooled=True, minor=True))

    # -- shave, only when a caller PINS a scale that does not fit (rbrowse `capped`)
    room = _band_room(len(bands), bands_mm, band_gap_mm)
    total = sum(b.h_mm for b in bands)
    if bands and total > room + 1e-9:
        excess = total - room
        plan.shaved_mm = excess
        headroom = [max(0.0, b.h_mm - scale.floor_mm) for b in bands]
        tot_head = sum(headroom)
        for b, hd in zip(bands, headroom):
            if tot_head <= 0:
                break
            cut = excess * hd / tot_head
            if cut > 1e-9:
                b.h_mm = max(scale.floor_mm, b.h_mm - cut)
                b.capped = True
    plan.bands = bands
    return plan


# ---------------------------------------------------------------------------
# drawing: one panel
# ---------------------------------------------------------------------------
def _squished_rows(ax, reads: Sequence[Read], plan: PanelPlan, *, role, region, xform, zorder):
    """Per-read rows at a pitch the budget can honour, with the LINE WIDTH DERIVED FROM
    THE PITCH. ``tracks.reads``'s fixed 4.0 pt overdraws 2.2x at 40 reads and 22x at 400
    in a 1 in strip (Chanfreau planning/880 checkpoint 1b) -- the same class of defect
    rbrowse fixed on its raster (HANDOFF.md:1961)."""
    S = TOK.stroke()
    X = xform or (lambda p: p)
    col = _body_color(role)
    tail_col = TOK.color("tail")
    order = sorted(reads, key=lambda r: (r.three_prime, r.five_prime))
    for i, r in enumerate(order):
        y = plan.bands_mm - (i + 0.5) * plan.pitch_mm
        for bs, bl in r.blocks:
            ax.plot([X(bs), X(bs + bl)], [y, y], color=col, lw=plan.lw_pt, solid_capstyle="butt",
                    zorder=zorder)
        prev = None
        for bs, bl in r.blocks:
            if prev is not None:
                ax.plot([X(prev), X(bs)], [y, y], color=TOK.color("hairline"), lw=S["hairline"],
                        zorder=zorder - 1)
            prev = bs + bl
        if r.tail > 0:
            t0 = r.end if r.strand == "+" else r.start - r.tail
            ax.plot([X(t0), X(t0 + r.tail)], [y, y], color=tail_col, lw=plan.lw_pt,
                    solid_capstyle="butt", zorder=zorder + 1)


def read_panel(ax, reads: Sequence[Read], keep: Sequence[Cluster] = (), *,
               budget_mm: float, scale: Optional[BandScale] = None, name: str = "",
               role: Optional[str] = None, region=None, xform: XForm = None,
               merged_mm: float = 2.5, merged_gap_mm: float = 0.8, band_gap_mm: float = 0.4,
               win: int = 32, mode: str = "auto", plan: Optional[PanelPlan] = None,
               guides: bool = True, zorder: int = 3) -> PanelPlan:
    """ONE sample inside ``budget_mm`` of page: a merged raster row over up to k ribbon
    bands (or per-read rows when the budget can show them at a legible pitch).

    The axes is given ``ylim = (0, budget_mm)``, so **one data unit is one millimetre**.
    Pass ``plan`` to reuse a plan computed by :func:`read_stack` (the shared band scale);
    otherwise the scale is solved for this panel alone and heights are NOT comparable
    across panels -- rbrowse says so out loud in that case (rbrowse.js:16648) and so
    should any caller.

    Returns the :class:`PanelPlan` actually drawn.
    """
    G = TOK.track_geometry()
    if plan is None:
        if scale is None:
            bm = max(0.0, budget_mm - merged_mm - merged_gap_mm)
            fl = G.get("band_floor_mm", 0.8)
            scale = solve_band_scale([effective_counts(reads, keep, bands_mm=bm, floor_mm=fl,
                                                       band_gap_mm=band_gap_mm, win=win)],
                                     bm, floor_mm=fl, band_gap_mm=band_gap_mm)
        plan = plan_read_panel(name, reads, keep, budget_mm=budget_mm, scale=scale,
                               merged_mm=merged_mm, merged_gap_mm=merged_gap_mm,
                               band_gap_mm=band_gap_mm, win=win, mode=mode)
    _apply_xlim(ax, region)
    ax.set_ylim(0.0, budget_mm)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)

    if plan.mode == "squished":
        _squished_rows(ax, reads, plan, role=role, region=region, xform=xform, zorder=zorder)
    else:
        y = plan.bands_mm
        for b in plan.bands:
            y -= b.h_mm
            ribbon(ax, b.reads, y=y, h=b.h_mm, anchor=b.anchor, role=role, minor=b.minor,
                   capped=b.capped, region=region, xform=xform, letter=None, zorder=zorder)
            y -= plan.band_gap_mm
    if plan.merged_mm > 0 and reads:
        merged_reads(ax, reads, y=plan.bands_mm + plan.merged_gap_mm, h=plan.merged_mm,
                     role=role, region=region, xform=xform, zorder=zorder)
    if guides:
        _panel_guides(ax, keep, xform=xform, top=budget_mm)
    return plan


def _panel_guides(ax, keep: Sequence[Cluster], *, xform: XForm, top: float):
    """The faint vertical registration line at each kept cluster's anchor, run down through
    the panel. rbrowse.js:17262: 'one line at cluster D's anchor lets the reader see D in
    all four genotypes at once, which IS the comparison the figure exists for.' Grey and
    faint by design -- it is a registration mark, never a data channel (:17293)."""
    X = xform or (lambda p: p)
    S = TOK.stroke()
    for c in keep:
        if c.anchor is None:
            continue
        ax.plot([X(c.anchor), X(c.anchor)], [0, top], color=TOK.color("mute"), lw=S["hairline"],
                alpha=0.55, zorder=1, clip_on=True)


# ---------------------------------------------------------------------------
# drawing: N samples under ONE model
# ---------------------------------------------------------------------------
def stack_height_mm(n_samples: int, *, panel_mm: float, panel_gap_mm: float,
                    model_mm: float = 0.0, model_gap_mm: float = 0.0,
                    axis_mm: float = 0.0, top_mm: float = 0.0,
                    reserved_mm: float = 0.0) -> float:
    """The page height an N-sample read stack needs -- computed BEFORE anything is drawn.

    rbrowse.js:10298: a banded panel has no content height; it takes a budget. So the
    figure height is the sum of the budgets plus the declared chrome, never a function of
    the read count. ``reserved_mm`` is the caller's legend + footer block, which is NOT
    part of the stack (`figstyle.concise_legend` only knows its height after layout).
    """
    n = max(0, int(n_samples))
    body = n * panel_mm + max(0, n - 1) * panel_gap_mm
    return top_mm + model_mm + model_gap_mm + body + axis_mm + reserved_mm


def read_stack(fig, samples: Sequence[dict], model=None, *, region=None, xform: XForm = None,
               panel_mm: float = 12.0, panel_gap_mm: float = 2.5, model_mm: float = 7.0,
               model_gap_mm: float = 2.0, axis_mm: float = 8.0, top_mm: float = 3.0,
               reserved_mm: float = 0.0, left_mm: float = 4.0, label_mm: float = 20.0,
               right_mm: float = 3.0, k: int = 4, win: int = 32, mode: str = "auto",
               letters: str = "ABCDEFGHIJ", band_gap_mm: float = 0.4,
               merged_mm: float = 2.5, merged_gap_mm: float = 0.8,
               size_figure: bool = True, axis_label: Optional[str] = None,
               origin_mm: Optional[float] = None) -> dict:
    """N sample read panels under ONE gene model, ONE letter key and ONE axis.

    ``samples`` is ``[{"name": str, "reads": [Read, ...], "role": Optional[str]}, ...]``.

    The economy rules (Chanfreau planning/880):

    (a) the model is drawn ONCE above the stack -- it is annotation, identical in every
        sample, and N copies cost N x 19.1 mm for one model's worth of content;
    (b) each panel's height is ``panel_mm``, a budget, not a function of its read count;
    (c) band heights come from ONE ``BandScale`` solved across every panel, floored at
        ``tokens.geometry.track.band_floor_mm``, minor clusters pooled into `other`;
    (d) the isoform letters ride the MODEL's 3' marks and a faint guide runs down the
        stack -- no per-band labels;
    (e) one ``name / n=`` chip per panel in the left gutter, nothing else in words inside;
    (f) ``mode="auto"`` gives per-read rows only when
        ``n_reads * pitch_floor_mm <= panel band budget``, else the merged raster.

    With ``size_figure`` the figure is RESIZED to ``stack_height_mm(...)`` at its current
    width; pass ``False`` to lay the stack into a figure you sized yourself, and
    ``origin_mm`` (the stack's TOP edge, in mm from the figure bottom) to place it in a band
    of a taller page. A sample dict may carry ``"note"``, appended to its chip -- use it when
    the panel draws fewer reads than the sample has (rbrowse's rule: the legend reports the
    number actually DRAWN, HANDOFF.md:1961).

    Returns ``{"height_mm", "panels": [...], "model_ax", "axis_ax", "scale", "keep",
    "plans", "used_mm"}``.
    """
    G = TOK.track_geometry()
    T = TOK.typography()
    n = len(samples)
    has_model = model is not None
    m_mm = model_mm if has_model else 0.0
    m_gap = model_gap_mm if has_model else 0.0
    H_mm = stack_height_mm(n, panel_mm=panel_mm, panel_gap_mm=panel_gap_mm, model_mm=m_mm,
                           model_gap_mm=m_gap, axis_mm=axis_mm, top_mm=top_mm,
                           reserved_mm=reserved_mm)
    W_in = fig.get_figwidth()
    if size_figure:
        fig.set_size_inches(W_in, H_mm / MM_PER_IN)
    W_mm = W_in * MM_PER_IN
    Hf_mm = fig.get_figheight() * MM_PER_IN

    x0 = (left_mm + label_mm) / W_mm
    x1 = 1.0 - right_mm / W_mm
    wf = x1 - x0

    def rect(bottom_mm, h_mm):
        return [x0, bottom_mm / Hf_mm, wf, h_mm / Hf_mm]

    keep = union_clusters(samples, win=win, k=k, letters=letters)

    # -- one shared band scale across every panel (rbrowse bandScale:'shared')
    bands_mm = max(0.0, panel_mm - (merged_mm + merged_gap_mm if merged_mm > 0 else 0.0))
    floor_mm = G.get("band_floor_mm", 0.8)
    panel_counts = [effective_counts(s.get("reads", ()), keep, bands_mm=bands_mm,
                                     floor_mm=floor_mm, band_gap_mm=band_gap_mm, win=win)
                    for s in samples]
    scale = solve_band_scale(panel_counts, bands_mm, floor_mm=floor_mm, band_gap_mm=band_gap_mm)

    # the stack occupies the band [axis_mm + reserved_mm, origin]
    y = (Hf_mm - top_mm) if origin_mm is None else origin_mm
    out: dict = {"height_mm": H_mm, "panels": [], "plans": [], "scale": scale, "keep": keep,
                 "model_ax": None, "axis_ax": None}

    if has_model:
        y -= m_mm
        axm = fig.add_axes(rect(y, m_mm))
        _draw_stack_model(axm, model, keep, region=region, xform=xform, height_mm=m_mm)
        out["model_ax"] = axm
        y -= m_gap

    used = 0.0
    for i, s in enumerate(samples):
        y -= panel_mm
        ax = fig.add_axes(rect(y, panel_mm))
        plan = plan_read_panel(s.get("name", f"s{i}"), s.get("reads", ()), keep,
                               budget_mm=panel_mm, scale=scale, merged_mm=merged_mm,
                               merged_gap_mm=merged_gap_mm, band_gap_mm=band_gap_mm,
                               win=win, mode=mode)
        read_panel(ax, s.get("reads", ()), keep, budget_mm=panel_mm, plan=plan,
                   role=s.get("role"), region=region, xform=xform, guides=True)
        ax.tick_params(labelbottom=False, bottom=False)
        ax.spines["bottom"].set_visible(False)
        # (e) ONE chip per panel, in the left gutter -- reserved width, never an overhang
        chip = f"{s.get('name', '')}\nn={plan.n_reads:,}"
        if plan.n_pooled:
            chip += f" · other {plan.n_pooled:,}"
        if plan.mode == "squished":
            chip += f"\n{plan.pitch_mm:.2f} mm/read"
        if s.get("note"):
            chip += f"\n{s['note']}"
        fig.text(left_mm / W_mm, (y + panel_mm / 2) / Hf_mm, chip, ha="left", va="center",
                 fontsize=T["annotation"], color=TOK.color("ink"), linespacing=1.35)
        out["panels"].append(ax); out["plans"].append(plan)
        used = max(used, plan.used_mm)
        if i < len(samples) - 1:
            y -= panel_gap_mm

    if axis_mm > 0:
        # the axis axes is a HAIRLINE strip at the TOP of its band, so its tick labels and
        # xlabel hang DOWN into the band that was budgeted for them -- placing a full-height
        # axes at the bottom of the band puts the labels off the page (found in the smoke).
        rule_mm, rule_gap_mm = 0.4, 2.0
        axa = fig.add_axes(rect(y - rule_mm - rule_gap_mm, rule_mm))
        axa.set_ylim(0, 1); axa.set_yticks([])
        for side in ("left", "right", "top"):
            axa.spines[side].set_visible(False)
        _apply_xlim(axa, region)
        if xform is None and region is not None:
            from .tracks import region_axis
            region_axis(axa, region, label=False)
        if axis_label:
            axa.set_xlabel(axis_label, fontsize=T["axis_label"])
        out["axis_ax"] = axa
    out["used_mm"] = used
    return out


def _draw_stack_model(ax, model, keep: Sequence[Cluster], *, region, xform, height_mm: float):
    """The one gene model above the stack, with the isoform LETTERS riding its ``polya``
    marks -- rule (d): the letter key and the model's 3' marks are the same object, so the
    strip costs no row of its own."""
    from .tracks import gene_model, mark
    T = TOK.typography()
    ax.set_ylim(-1.0, 1.6)
    ax.set_yticks([])
    for side in ("left", "right", "top", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelbottom=False, bottom=False)
    _apply_xlim(ax, region)
    gene_model(ax, model, y=0.0, height=0.55, region=region, xform=xform, tss=True,
               label_pos="left")
    for c in keep:
        if c.anchor is None:
            continue
        mark(ax, c.anchor, "polya", y=0.0, height=0.55, label=c.letter, xform=xform)


# ============================================================================
# the DENSE PILE -- rbrowse's merged density (drawBandRaster) on the publication path
# ============================================================================
def sort_reads(reads: Sequence[Read], order: str = "3p") -> List[Read]:
    """rbrowse's row order. ``"5p"`` is the browser's flat merged default (``sortBy:
    'five'``: rows by 5' end, 3' ends unsorted within); ``"3p"`` groups by 3' end first,
    5' end within -- the order that turns the pile's 3' edge into the isoform staircase.
    Both run in TRANSCRIPT direction so the most 5'-reaching read is on top either way."""
    if order not in ("3p", "5p"):
        raise ValueError("order must be '3p' or '5p'")
    def k(r):
        sgn = 1 if r.strand == "+" else -1
        a, b = sgn * r.three_prime, sgn * r.five_prime
        return (a, b) if order == "3p" else (b, a)
    return sorted(reads, key=k)


def dense_pile(ax, reads: Sequence[Read], *, y: float = 0.0, h: float = 1.0, role: Optional[str] = None,
               region=None, xform: XForm = None, order: str = "3p", row_mm: Optional[float] = None,
               nbins: Optional[int] = None, zorder: int = 3):
    """rbrowse's **merged density** as a publication glyph: the reads in ``order`` are cut into
    as many chunks as the strip has pixel rows, and each chunk is one row of the per-column
    MAJORITY raster (:func:`merged_reads`). Every read contributes, nothing is sampled, and
    a row is never thinner than a device pixel -- so 5,000 reads in a 12 mm strip draw as a
    ~140-row density picture whose edges are the 5' and 3' end distributions
    (rbrowse ``drawBandRaster``, HANDOFF.md:1961's fix for the overflowing per-read pile).

    ``row_mm`` is the minimum row height on the page (default one 300-dpi pixel, 0.0847 mm);
    :func:`merged_reads` is the ``n_rows == 1`` case. The pile carries no legible single read,
    by design: what it shows is the shape of the population. Returns the ``AxesImage``.
    """
    from matplotlib.colors import to_rgb
    r = _apply_xlim(ax, region)
    x0, x1 = (r.start, r.end) if r is not None else sorted(ax.get_xlim())
    n = len(reads)
    if n == 0:
        return None
    nb = nbins or _nbins_for(ax, x0, x1)
    fig = ax.figure
    strip_mm = h / abs(np.diff(ax.get_ylim())[0]) * ax.get_position().height * fig.get_figheight() * MM_PER_IN \
        if np.diff(ax.get_ylim())[0] else 0.0
    rm = row_mm if row_mm is not None else 25.4 / 300.0
    n_rows = int(max(1, min(n, np.floor(strip_mm / rm)))) if strip_mm > 0 else min(n, 200)
    per = int(np.ceil(n / n_rows))
    n_rows = int(np.ceil(n / per))
    G = TOK.track_geometry()
    body_rgb = np.array(to_rgb(_body_color(role)))
    tail_rgb = np.array(to_rgb(TOK.color("tail")))
    clip_rgb = np.array(to_rgb(TOK.color("mute")))
    a0 = G.get("merged_alpha_floor", 0.2)
    ordered = sort_reads(reads, order)
    img = np.zeros((n_rows, nb, 4))
    for k in range(n_rows):
        chunk = ordered[k * per:(k + 1) * per]
        m = len(chunk)
        cb, ct, cc = _columns(chunk, x0, x1, nb, xform)
        tail_wins = (ct > 0) & (ct * 4 >= cb)
        clip_wins = (cc > 0) & (cc * 4 >= cb) & ~tail_wins
        body = (cb > 0) & ~tail_wins & ~clip_wins
        row = img[k]
        row[tail_wins, :3] = tail_rgb; row[tail_wins, 3] = np.minimum(1, 0.35 + 0.65 * ct[tail_wins] / m)
        row[clip_wins, :3] = clip_rgb; row[clip_wins, 3] = np.minimum(1, 0.25 + 0.55 * cc[clip_wins] / m)
        row[body, :3] = body_rgb; row[body, 3] = np.minimum(1, a0 + (1 - a0) * cb[body] / m)
    X = xform or (lambda p: p)
    lo, hi = sorted((X(x0), X(x1)))
    # row 0 (the most 5'-reaching reads) at the TOP of the strip, like the browser
    im = ax.imshow(img, extent=(lo, hi, y, y + h), aspect="auto", interpolation="nearest",
                   origin="upper", zorder=zorder)
    if r is not None:
        ax.set_xlim(*r.xlim())
    return im


# ============================================================================
# quantification WINDOWS -- the studio's shaded bands, and what they count
# ============================================================================
@dataclass
class Window:
    """One quantification window: every 3' end in ``[lo, hi]`` (inclusive) counts for
    ``letter``. rbrowse: "a letter is an INTERVAL, not a point" (rbrowse.js:17397)."""
    letter: str
    lo: int
    hi: int
    anchor: Optional[int] = None

    def holds(self, r: Read) -> bool:
        return self.lo <= r.three_prime <= self.hi


def windows_from_clusters(keep: Sequence[Cluster], *, strand: str = "+", grow: int = 0,
                          letters: str = "ABCDEFGHIJ", order: str = "position") -> List[Window]:
    """Windows from the stack's kept clusters: each window is the cluster's own 3'-end extent
    (``Cluster.span``), grown outward by ``grow`` bp but never past the midpoint to its
    neighbour (rbrowse ``figUsage`` orphan absorption, without the annotation-class rule).
    ``order="position"`` letters them in TRANSCRIPT order 5'->3' (the studio's
    ``clusterOrder: 'position'``); ``"rank"`` keeps the pooled-count order the stack uses."""
    cl = [c for c in keep if c.span is not None]
    if not cl:
        return []
    if order == "position":
        cl = sorted(cl, key=lambda c: c.anchor if strand == "+" else -c.anchor)
    elif order != "rank":
        raise ValueError("order must be 'position' or 'rank'")
    spans = [list(c.span) for c in cl]
    by_pos = sorted(range(len(cl)), key=lambda i: spans[i][0])
    for k, i in enumerate(by_pos):
        lo, hi = spans[i]
        lo2 = lo - grow; hi2 = hi + grow
        if k > 0:
            j = by_pos[k - 1]; lo2 = max(lo2, (spans[j][1] + lo) // 2 + 1)
        if k < len(by_pos) - 1:
            j = by_pos[k + 1]; hi2 = min(hi2, (hi + spans[j][0]) // 2)
        spans[i] = [min(lo, lo2), max(hi, hi2)]
    return [Window(letters[i] if i < len(letters) else str(i), int(s[0]), int(s[1]), anchor=c.anchor)
            for i, (c, s) in enumerate(zip(cl, spans))]


def window_bands(ax, windows: Sequence[Window], *, y0: float, y1: float, xform: XForm = None,
                 zorder: int = 1) -> list:
    """The studio's shaded quantification windows over a read panel: a ``wash``-grey band
    across ``[lo, hi]`` with a slightly stronger edge, a registration mark and never a data
    colour (rbrowse.js:17397, Kevin 2026-08-13: "faint grey rectangles ... rather than a
    single vertical line"). Draw BEHIND the reads (low ``zorder``)."""
    X = xform or (lambda p: p)
    S = TOK.stroke()
    out = []
    for w in windows:
        a, b = sorted((X(w.lo), X(w.hi + 1)))
        out.append(ax.add_patch(Rectangle((a, y0), b - a, y1 - y0, facecolor=TOK.color("neutral"), alpha=0.14,
                                          edgecolor="none", zorder=zorder)))
        for x in (a, b):
            out.append(ax.plot([x, x], [y0, y1], color=TOK.color("neutral"), lw=S["hairline"], alpha=0.5,
                               zorder=zorder)[0])
    return out


def window_counts(samples: Sequence[dict], windows: Sequence[Window]) -> Dict[str, Dict[str, int]]:
    """``{sample name: {letter: n, "other": n, "total": n}}`` -- what each window counts in
    each sample. Every read is counted exactly once (windows are disjoint by construction)."""
    out: Dict[str, Dict[str, int]] = {}
    for s in samples:
        d = {w.letter: 0 for w in windows}
        d["other"] = 0
        for r in s.get("reads", ()):
            hit = next((w for w in windows if w.holds(r)), None)
            d[hit.letter if hit else "other"] += 1
        d["total"] = len(s.get("reads", ()))
        out[s.get("name", "")] = d
    return out


def rep_values(samples: Sequence[dict], windows: Sequence[Window], *, depth: Optional[Dict[str, float]] = None,
               scale: Optional[Dict[str, float]] = None) -> List[dict]:
    """Per replicate, per window (plus ``"other"``): the count, the SHARE of that replicate's
    reads at the locus, and -- when ``depth[name]`` (the library's denominator, e.g. its
    protein-coding read count) is given -- reads per million (``ppm``). ``scale[name]``
    multiplies the count first: the slicer's uniform subsample factor ``n_before_cap /
    n_kept``, so a capped shard is restored before it is normalised (Chanfreau planning/885b).

    Rows are ``{"rep", "group", "letter", "n", "share", "total"[, "ppm"]}`` -- the input
    :func:`panels.reps_as_points` takes. A replicate with no reads has share 0 in every
    window, not NaN, so it still counts as one point.
    """
    counts = window_counts(samples, windows)
    rows: List[dict] = []
    for s in samples:
        name = s.get("name", "")
        c = counts[name]
        tot = c["total"] or 1
        for w in list(windows) + [None]:
            key = w.letter if w is not None else "other"
            n = c[key]
            row = {"rep": name, "group": s.get("group"), "letter": key, "n": n, "share": n / tot,
                   "total": c["total"]}
            if depth and name in depth:
                f = (scale or {}).get(name, 1.0)
                row["ppm"] = 1e6 * n * f / depth[name]
            rows.append(row)
    return rows


# ============================================================================
# TRUNCATION SURVIVAL -- correcting a library's signal for its own 5' truncation hazard
# ============================================================================
@dataclass
class Survival:
    """A library's 5'-truncation survival, S(d) = 2^(-d / t_half): the probability that a
    molecule's aligned body still covers a position ``d`` SPLICED nucleotides upstream of its
    3' end.

    The model is Sumner planning/136: an exponential per-nucleotide hazard ``lambda`` that the
    observed 5' end is missing (t_half = ln 2 / lambda), fitted per library on the spliced span
    of reads whose 3' end sits at an annotated transcript end. Direct RNA is 3'-anchored, so
    every 5'-end signal a library shows is the true 5' ends TIMES this survival, plus the
    truncation events themselves; the correction below undoes both.

    ``max_weight`` caps 1/S: beyond ~3 half-lives an observed read stands for eight molecules
    and the estimate is noise, so the cap is declared, not hidden. The default is
    ``tokens.geometry.track.survival_max_weight`` (8); print the cap wherever the correction
    is applied (R6, Chanfreau planning/885d).
    """
    t_half_nt: float
    max_weight: Optional[float] = None

    def __post_init__(self):
        if self.t_half_nt <= 0:
            raise ValueError(f"Survival: t_half_nt must be > 0, got {self.t_half_nt}")
        if self.max_weight is None:
            # the declared cap: tokens.geometry.track.survival_max_weight (8 = three half-lives)
            self.max_weight = float(TOK.track_geometry().get("survival_max_weight", 8.0))

    @property
    def hazard(self) -> float:
        """lambda, per spliced nucleotide."""
        return float(np.log(2.0) / self.t_half_nt)

    def s(self, d):
        return np.power(2.0, -np.asarray(d, dtype=float) / self.t_half_nt)

    def w(self, d):
        """Inverse-probability weight 1/S(d), capped."""
        return np.minimum(self.max_weight, 1.0 / self.s(d))


def spliced_distance(r: Read, x) -> np.ndarray:
    """Spliced distance (aligned nucleotides, introns excluded) from the read's 3' end to each
    position ``x`` -- the exposure the hazard runs on. Positions past the 5' end still get the
    distance to the 5' end (the read's full spliced span)."""
    x = np.asarray(x, dtype=float)
    d = np.zeros_like(x)
    if r.strand == "+":
        e = r.end
        for bs, bl in r.blocks:
            lo, hi = bs, bs + bl
            d += np.clip(np.minimum(hi, e) - np.maximum(x, lo), 0, None)
    else:
        s0 = r.start
        for bs, bl in r.blocks:
            lo, hi = bs, bs + bl
            d += np.clip(np.minimum(x + 1, hi) - np.maximum(lo, s0), 0, None)
    return d


def _weighted_columns(reads: Sequence[Read], x0: float, x1: float, nbins: int, survival: Survival,
                      xform: XForm = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-bin (weighted body coverage, weighted 5'-end count, raw 5'-end count): each read
    covering a bin centre counts 1/S(d) molecules, d its spliced distance from the read's 3'
    end to that centre. Bins are in display space when ``xform`` is given."""
    X = xform or (lambda p: p)
    lo, hi = sorted((X(x0), X(x1)))
    edges = np.linspace(lo, hi, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    cov_w = np.zeros(nbins); ends_w = np.zeros(nbins); ends_raw = np.zeros(nbins)
    for r in reads:
        xs = [(X(bs), X(bs + bl), bs, bl) for bs, bl in r.blocks]
        for a, b, bs, bl in xs:
            a, b = min(a, b), max(a, b)
            i0 = int(np.searchsorted(edges, a, side="right") - 1)
            i1 = int(np.searchsorted(edges, b, side="left"))
            i0 = max(0, i0); i1 = min(nbins, i1)
            if i1 <= i0:
                continue
            # genome positions of the covered bin centres (inverse of a linear xform is not
            # available in general, so the distance is taken along the block in genome space)
            frac = (centers[i0:i1] - a) / max(b - a, 1e-9)
            gx = bs + frac * bl if xform is None else bs + frac * bl
            cov_w[i0:i1] += survival.w(spliced_distance(r, gx))
        p5 = r.five_prime
        j = int(np.searchsorted(edges, X(p5), side="right") - 1)
        if 0 <= j < nbins:
            ends_raw[j] += 1
            ends_w[j] += float(survival.w(spliced_distance(r, np.array([p5])))[0])
    return cov_w, ends_w, ends_raw


def corrected_coverage(reads: Sequence[Read], x0: float, x1: float, nbins: int, survival: Survival,
                       xform: XForm = None) -> Tuple[np.ndarray, np.ndarray]:
    """(raw coverage, survival-corrected coverage) per bin -- the second is the number of
    molecules that would cover each position had the library not truncated them
    (Horvitz-Thompson: each surviving read stands for 1/S(d))."""
    cb, _, _ = _columns(reads, x0, x1, nbins, xform)
    cov_w, _, _ = _weighted_columns(reads, x0, x1, nbins, survival, xform)
    return cb, cov_w


def five_prime_profile(reads: Sequence[Read], x0: float, x1: float, nbins: int,
                       survival: Optional[Survival] = None, xform: XForm = None) -> Dict[str, np.ndarray]:
    """The 5'-end density per bin, raw and, with a ``survival``, deconvolved:

        true_ends(x) = [ends_w(x) - lambda * binwidth * cov_w(x)]_+

    where ``ends_w`` are the observed 5' ends each weighted by 1/S (de-attenuated) and the
    subtracted term is the hazard's expected TRUNCATION ends at x (lambda per nt times the
    de-attenuated molecules still covering x). Negative bins are clipped to zero and counted in
    ``n_clipped`` -- a place where the model overshoots the data, to be stated, not hidden.
    """
    X = xform or (lambda p: p)
    lo, hi = sorted((X(x0), X(x1)))
    edges = np.linspace(lo, hi, nbins + 1)
    out: Dict[str, np.ndarray] = {"edges": edges}
    if survival is None:
        raw = np.zeros(nbins)
        for r in reads:
            j = int(np.searchsorted(edges, X(r.five_prime), side="right") - 1)
            if 0 <= j < nbins:
                raw[j] += 1
        out["raw"] = raw
        return out
    cov_w, ends_w, ends_raw = _weighted_columns(reads, x0, x1, nbins, survival, xform)
    bw = (hi - lo) / nbins                      # in display units == genome nt when xform is None
    trunc = survival.hazard * bw * cov_w
    corr = ends_w - trunc
    out.update({"raw": ends_raw, "deattenuated": ends_w, "expected_truncation": trunc,
                "corrected": np.clip(corr, 0, None), "n_clipped": np.array([int((corr < 0).sum())])})
    return out


# ============================================================================
# THE COMPACT BLOCK -- piles under one model, replicate strips, the t1/2 ladder
# (Chanfreau planning/885; the 887 landing; Kevin's rulings R1-R8, 2026-09-05/06)
# ============================================================================
def _bare_axes(ax, *, bottom: bool = False):
    """A read / strip axes: no y axis, no spines except (optionally) the bottom rule."""
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    if not bottom:
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(labelbottom=False, bottom=False)


def _column(fig, *, left_mm: float, label_mm: float, right_mm: float) -> Tuple[float, float]:
    """(x0, width) of the drawing column as figure fractions, from the mm gutters: a label
    gutter on the left (the chips), a right margin, the same for every block of a page so
    the blocks register column for column."""
    W_mm = fig.get_figwidth() * MM_PER_IN
    x0 = (left_mm + label_mm) / W_mm
    return x0, 1.0 - right_mm / W_mm - x0


def _check_role(role: Optional[str]) -> None:
    """A sample's role is a LAYER-A role or None -- refused up front, even for an empty
    sample (the junction_arcs precedent: an empty panel must not silently accept a
    molecular identity as a sample colour)."""
    if role is not None:
        TOK.role(role)


def pile_stack(fig, samples: Sequence[dict], model=None, *, region=None, xform: XForm = None,
               origin_mm: Optional[float] = None, panel_mm: float = 12.0, panel_gap_mm: float = 2.5,
               model_mm: float = 7.0, model_gap_mm: float = 2.0, axis_mm: float = 8.0,
               top_mm: float = 2.0, left_mm: float = 4.0, label_mm: float = 22.0,
               right_mm: float = 3.0, order: str = "3p", strip_mm: float = 0.0,
               strip_gap_mm: float = 0.6, windows: Sequence[Window] = (),
               axis_label: Optional[str] = None, chip: bool = True, model_marks: bool = True,
               neighbours: Sequence = ()) -> dict:
    """N sample panels under ONE model, each panel a DENSE PILE (:func:`dense_pile`,
    rbrowse's merged density) on a mm budget, the quantification ``windows`` shaded down
    the stack and lettered on the model's ``polya`` marks, one chip per panel, one axis.

    This is the read block of the COMPACT BLOCK -- the default showcase for DRS / cDNA in
    a figure (Chanfreau planning/885b/885c, Kevin's ruling R5, 2026-09-06): one model with
    the window letters, one strip per library (:func:`strip_rows`), one 3'-sorted pile per
    group in the two sample roles, the windows shaded, one axis, then per-window bars with
    the libraries as points (``panels.reps_as_points``).

    ``samples`` is ``[{"name", "reads", "role"[, "note"]}, ...]``; the geometry mirrors
    :func:`read_stack` (same gutters, the same ``panel_mm`` budget, ``top_mm`` above the
    model, ``axis_mm`` below the last panel) so the two forms compare mm for mm, and
    :func:`stack_height_mm` gives the block's height BEFORE drawing. ``origin_mm`` is the
    block's TOP edge in mm from the figure bottom (default: the page top); ``strip_mm`` > 0
    puts the one-row majority strip of every read above the pile INSIDE the budget;
    ``order`` is ``"3p"`` (the isoform staircase; the figure default) or ``"5p"`` (the
    browser's truncation ramp). A layer-B token as a role is refused. Only the axis band
    may carry a label, and only when ``axis_label`` is given -- one axis label per page, a
    stack's axis band holds tick labels only (R8; 881a render 1).

    Returns ``{"panels", "model_ax", "axis_ax", "windows", "rows", "top_mm", "bottom_mm",
    "height_mm"}``; ``rows`` is the raster row count each pile drew.
    """
    from .tracks import gene_model, mark, region_axis
    T = TOK.typography()
    W_mm = fig.get_figwidth() * MM_PER_IN
    H_mm = fig.get_figheight() * MM_PER_IN
    x0, wf = _column(fig, left_mm=left_mm, label_mm=label_mm, right_mm=right_mm)
    for s in samples:
        _check_role(s.get("role"))

    def rect(bottom_mm, h_mm):
        return [x0, bottom_mm / H_mm, wf, h_mm / H_mm]

    y_top = H_mm if origin_mm is None else origin_mm
    y = y_top - top_mm
    out: dict = {"panels": [], "model_ax": None, "axis_ax": None, "windows": list(windows), "rows": [],
                 "top_mm": y_top}
    if model is not None:
        y -= model_mm
        axm = fig.add_axes(rect(y, model_mm))
        axm.set_ylim(-1.0, 1.6)
        _bare_axes(axm)
        _apply_xlim(axm, region)
        gene_model(axm, model, y=0.0, height=0.55, region=region, xform=xform, tss=True, label_pos="left")
        for nb in neighbours:
            gene_model(axm, nb, y=0.0, height=0.55, region=region, xform=xform, tss=True, label_pos="above")
        if model_marks:
            for w in windows:
                at = w.anchor if w.anchor is not None else (w.lo + w.hi) // 2
                mark(axm, at, "polya", y=0.0, height=0.55, label=w.letter, xform=xform)
        out["model_ax"] = axm
        y -= model_gap_mm
    for i, s in enumerate(samples):
        y -= panel_mm
        ax = fig.add_axes(rect(y, panel_mm))
        ax.set_ylim(0.0, panel_mm)              # one data unit = one millimetre (read_panel's convention)
        _bare_axes(ax)
        _apply_xlim(ax, region)
        pile_h = panel_mm - (strip_mm + strip_gap_mm if strip_mm > 0 else 0.0)
        window_bands(ax, windows, y0=0.0, y1=panel_mm, xform=xform, zorder=1)
        reads = s.get("reads", ())
        n_rows = 0
        if reads:
            im = dense_pile(ax, reads, y=0.0, h=pile_h, role=s.get("role"), region=region, xform=xform,
                            order=order, zorder=3)
            n_rows = int(im.get_array().shape[0]) if im is not None else 0
            if strip_mm > 0:
                merged_reads(ax, reads, y=pile_h + strip_gap_mm, h=strip_mm, role=s.get("role"),
                             region=region, xform=xform, zorder=3)
        out["rows"].append(n_rows)
        if chip:
            txt = f"{s.get('name', '')}\nn={len(reads):,}"
            if s.get("note"):
                txt += f"\n{s['note']}"
            fig.text(left_mm / W_mm, (y + panel_mm / 2) / H_mm, txt, ha="left", va="center",
                     fontsize=T["annotation"], color=TOK.color("ink"), linespacing=1.35)
        out["panels"].append(ax)
        if i < len(samples) - 1:
            y -= panel_gap_mm
    if axis_mm > 0:
        # a hairline strip at the TOP of the axis band so the tick labels hang DOWN into it
        rule_mm, rule_gap_mm = 0.4, 2.0
        axa = fig.add_axes(rect(y - rule_mm - rule_gap_mm, rule_mm))
        axa.set_ylim(0, 1)
        axa.set_yticks([])
        for side in ("left", "right", "top"):
            axa.spines[side].set_visible(False)
        _apply_xlim(axa, region)
        if xform is None and region is not None:
            region_axis(axa, region, label=False)
        if axis_label:
            axa.set_xlabel(axis_label, fontsize=T["axis_label"])
        out["axis_ax"] = axa
        y -= axis_mm
    out["bottom_mm"] = y
    out["height_mm"] = y_top - y
    return out


def strip_rows(fig, samples: Sequence[dict], *, region=None, xform: XForm = None, origin_mm: float,
               row_mm: Optional[float] = None, gap_mm: float = 0.5, left_mm: float = 4.0,
               label_mm: float = 22.0, right_mm: float = 3.0, windows: Sequence[Window] = (),
               label_every: bool = True, group_labels: bool = True):
    """One majority strip per sample, one row each, every read -- the replicate heatmap
    (Chanfreau planning/885: "do the replicates agree" is what a strip block answers; the
    tail blocks and extents line up row for row, or they do not).

    ``row_mm`` defaults to ``tokens.geometry.track.strip_row_mm`` (1.6 mm): the smallest
    row at which a tail block and a body fraction still read in print. A 7.5 pt label is
    taller than that pitch, so rows are labelled by GROUP -- a sample dict carrying
    ``"group"`` gets one label per run of consecutive rows (name, library count, reads) with
    a hairline bracket -- and never one label per row; samples without a group get a
    per-row label when ``label_every``. A sample may carry ``"survival"`` (a
    :class:`Survival`): its strip is then the survival-CORRECTED coverage, a MODEL, which
    the page must declare (R6: the raw strip stays on the page, the parameter is printed
    where it applies -- :func:`halflife_ladder`). Returns ``(axes, bottom_mm)``.
    """
    T = TOK.typography()
    S = TOK.stroke()
    G = TOK.track_geometry()
    rm = float(G.get("strip_row_mm", 1.6)) if row_mm is None else float(row_mm)
    W_mm = fig.get_figwidth() * MM_PER_IN
    H_mm = fig.get_figheight() * MM_PER_IN
    x0, wf = _column(fig, left_mm=left_mm, label_mm=label_mm, right_mm=right_mm)
    for s in samples:
        _check_role(s.get("role"))
    n = len(samples)
    h = n * rm + max(0, n - 1) * gap_mm
    y = origin_mm - h
    ax = fig.add_axes([x0, y / H_mm, wf, h / H_mm])
    ax.set_ylim(0.0, h)
    _bare_axes(ax)
    _apply_xlim(ax, region)
    window_bands(ax, windows, y0=0.0, y1=h, xform=xform, zorder=1)
    tops: Dict[int, float] = {}
    for i, s in enumerate(samples):
        yy = h - (i + 1) * rm - i * gap_mm
        tops[i] = yy
        merged_reads(ax, s.get("reads", ()), y=yy, h=rm, role=s.get("role"), region=region, xform=xform,
                     zorder=3, survival=s.get("survival"))
        if label_every and not s.get("group"):
            fig.text(left_mm / W_mm, (y + yy + rm / 2) / H_mm,
                     f"{s.get('name', '')} · n={len(s.get('reads', ())):,}", ha="left", va="center",
                     fontsize=T["annotation"], color=TOK.color("ink"))
    if group_labels and any(s.get("group") for s in samples):
        xl = ax.get_xlim()[0]
        i = 0
        while i < n:
            g = samples[i].get("group")
            j = i
            while j + 1 < n and samples[j + 1].get("group") == g:
                j += 1
            if g:
                y_top, y_bot = tops[i] + rm, tops[j]
                ntot = sum(len(samples[k].get("reads", ())) for k in range(i, j + 1))
                fig.text(left_mm / W_mm, (y + (y_top + y_bot) / 2) / H_mm,
                         f"{g} · {j - i + 1} libs\nn={ntot:,}", ha="left", va="center",
                         fontsize=T["annotation"], color=TOK.color("ink"), linespacing=1.15)
                ax.plot([xl, xl], [y_bot, y_top], color=TOK.color("hairline"), lw=S["hairline"],
                        clip_on=False, zorder=5)
            i = j + 1
    return ax, y


def halflife_ladder(fig, samples: Sequence[dict], thalf: Dict[str, float], *, y_bottom_mm: float,
                    row_mm: Optional[float] = None, gap_mm: float = 0.5, x_left_mm: float, width_mm: float,
                    lo: float = 250.0, hi: float = 5000.0, label: bool = True):
    """The truncation half-life ladder beside a strip block: one dot per row on a log axis,
    ALIGNED to the rows, in the sample's role -- the Sumner planning/136 idiom placed next
    to the heatmap it explains, so a corrected strip carries its parameter where it is
    applied (R6; planning/885d). ``thalf`` maps sample name -> t1/2 in nt. Every sample
    needs a value (``KeyError`` otherwise), and a value outside ``[lo, hi]`` is REFUSED
    rather than drawn off the axis -- widen the range. Returns the axes.
    """
    T = TOK.typography()
    S = TOK.stroke()
    G = TOK.track_geometry()
    rm = float(G.get("strip_row_mm", 1.6)) if row_mm is None else float(row_mm)
    missing = [s.get("name", "") for s in samples if s.get("name", "") not in thalf]
    if missing:
        raise KeyError(f"halflife_ladder: no t1/2 for {missing}")
    outside = [(s["name"], thalf[s["name"]]) for s in samples if not (lo <= thalf[s["name"]] <= hi)]
    if outside:
        raise ValueError(f"halflife_ladder: t1/2 outside [{lo}, {hi}] nt would draw off the axis: {outside}")
    W_mm = fig.get_figwidth() * MM_PER_IN
    H_mm = fig.get_figheight() * MM_PER_IN
    n = len(samples)
    h = n * rm + max(0, n - 1) * gap_mm
    ax = fig.add_axes([x_left_mm / W_mm, y_bottom_mm / H_mm, width_mm / W_mm, h / H_mm])
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(0.0, h)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(S["hairline"])
    ticks = [t for t in (300, 1000, 3000) if lo <= t <= hi]
    ax.set_xticks(ticks)
    ax.set_xticks([], minor=True)
    if label:
        ax.set_xticklabels([f"{t / 1000:g}" for t in ticks], fontsize=T["tick_label"])
        ax.tick_params(axis="x", length=2, width=S["hairline"], pad=1)
        ax.set_xlabel("t½ (kb)", fontsize=T["annotation"], labelpad=1)
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis="x", length=2, width=S["hairline"])
    for x in ticks:
        ax.plot([x, x], [0, h], color=TOK.color("wash"), lw=S["hairline"], zorder=1)
    for i, s in enumerate(samples):
        _check_role(s.get("role"))
        yy = h - (i + 1) * rm - i * gap_mm + rm / 2
        col = TOK.role(s["role"]) if s.get("role") else TOK.color("subtle")
        ax.scatter([thalf[s["name"]]], [yy], s=9, color=col, edgecolor=TOK.color("paper"),
                   linewidth=S["hairline"], zorder=4)
    return ax
