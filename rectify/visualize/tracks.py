"""
Layer-B track glyphs for RECTIFY visualization -- gene models, site marks,
coverage and junction arcs, drawn from PLAIN COORDINATES with colours from the
house token file.

This is the promoted, token-driven form of ``gene_track.py`` / ``coverage.py``
(Chanfreau planning/871, 2026-09-05). The older modules keep their APIs and
now delegate their colours here; this module is what a ``planning/`` script
imports so it stops hand-rolling an exon model out of ``Rectangle``::

    from rectify.visualize import tracks as T
    tx = T.Transcript("PDH1", "chrXVI", "+", exons=[(558219, 558423), (559811, 560137)])
    reg = T.Region("chrXVI", 558000, 560400)
    T.gene_model(ax, tx, region=reg, tss=True)
    T.mark(ax, 559811, "ss3")             # a 3'SS tick, layer-B `splice`
    T.mark(ax, 560137, "polya", label="pA")
    T.coverage(ax2, pos, depth, role="focal", region=reg)   # SIGNAL: a layer-A role
    T.region_axis(ax2, reg)

THE COLOUR RULE THIS MODULE ENFORCES
------------------------------------
* ``gene_model`` and ``mark`` take NO colour argument. A gene model is drawn in
  the grey ramp (exon = ``ink``, UTR = ``ink`` at half height, intron =
  ``hairline``) and a mark in its own layer-B identity colour (``splice``,
  ``polya``, ``mod``). Nothing about a molecule can be made to carry an argument.
* ``coverage``, ``strand_coverage``, ``reads`` and ``arc`` are SIGNAL glyphs --
  data about a sample -- so they take a ``role=`` (a layer-A argument role, or
  ``"neutral"``). With no role a signal is ``subtle`` grey. ``arc`` with no role
  is an ANNOTATION arc and is drawn in ``splice``.
* Strand is carried by POSITION (a minus-strand locus reverses the x-axis so it
  reads 5' -> 3'; ``strand_coverage`` mirrors below the axis) and by arrow
  direction -- never by hue.

COORDINATES are 0-based half-open, as everywhere in rectify. ``Region`` is
``(chrom, start, end[, strand])``; passing ``strand="-"`` makes ``region_axis``
reverse the x-axis so the locus reads 5' -> 3' left to right.

Author: Kevin R. Roy
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

try:
    from matplotlib.patches import PathPatch, Rectangle
    from matplotlib.path import Path as MPath
    from matplotlib.transforms import offset_copy
    MATPLOTLIB_AVAILABLE = True
except ImportError:  # pragma: no cover
    MATPLOTLIB_AVAILABLE = False

from . import tokens as TOK

__all__ = [
    "Transcript", "Region", "gene_model", "gene_track", "coverage", "strand_coverage",
    "arc", "mark", "reads", "region_axis", "MARK_KINDS",
]

Interval = Tuple[int, int]
XForm = Optional[Callable[[float], float]]


# ============================================================================
# Data carriers -- plain coordinates, no loader
# ============================================================================
@dataclass
class Transcript:
    """A transcript / gene model in genomic coordinates (0-based, half-open).

    ``exons`` are genomic intervals in ascending order. ``cds`` is the genomic
    (min, max) of the coding region; when given, the parts of each exon outside
    it are drawn as UTR at half height. ``introns`` is derived.
    """
    name: str
    chrom: str
    strand: str
    exons: List[Interval]
    cds: Optional[Interval] = None
    label: Optional[str] = None

    def __post_init__(self):
        if self.strand not in ("+", "-"):
            raise ValueError(f"strand must be '+' or '-', got {self.strand!r}")
        self.exons = sorted((int(s), int(e)) for s, e in self.exons)
        if not self.exons:
            raise ValueError("a Transcript needs at least one exon")
        for s, e in self.exons:
            if e <= s:
                raise ValueError(f"exon {s}-{e} is empty or reversed")
        if self.label is None:
            self.label = self.name

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def introns(self) -> List[Interval]:
        return [(self.exons[i][1], self.exons[i + 1][0]) for i in range(len(self.exons) - 1)]

    @property
    def tss(self) -> int:
        """The transcription start as a coordinate on the axis (5' end)."""
        return self.start if self.strand == "+" else self.end

    @property
    def tes(self) -> int:
        return self.end if self.strand == "+" else self.start

    @classmethod
    def from_introns(cls, name: str, chrom: str, strand: str, start: int, end: int,
                     introns: Iterable[Interval] = (), cds: Optional[Interval] = None,
                     label: Optional[str] = None) -> "Transcript":
        """Build from a gene span plus intron intervals -- the shape most Chanfreau
        scripts already hold (``gstart, gend, intron0=(s, e)``)."""
        introns = sorted((int(s), int(e)) for s, e in introns)
        exons, cur = [], int(start)
        for s, e in introns:
            if s < cur or e > end:
                raise ValueError(f"intron {s}-{e} lies outside {start}-{end} or overlaps another")
            exons.append((cur, s))
            cur = e
        exons.append((cur, int(end)))
        return cls(name, chrom, strand, exons, cds=cds, label=label)

    def utr_cds_parts(self) -> List[Tuple[Interval, bool]]:
        """Each exon split into (interval, is_cds) pieces."""
        if self.cds is None:
            return [((s, e), True) for s, e in self.exons]
        c0, c1 = self.cds
        parts = []
        for s, e in self.exons:
            if e <= c0 or s >= c1:
                parts.append(((s, e), False))
                continue
            if s < c0:
                parts.append(((s, c0), False))
            parts.append(((max(s, c0), min(e, c1)), True))
            if e > c1:
                parts.append(((c1, e), False))
        return parts


@dataclass
class Region:
    """A genomic window. ``strand="-"`` asks ``region_axis`` to reverse the x-axis."""
    chrom: str
    start: int
    end: int
    strand: Optional[str] = None

    def __post_init__(self):
        if self.end <= self.start:
            raise ValueError(f"empty region {self.start}-{self.end}")

    @property
    def span(self) -> int:
        return self.end - self.start

    @property
    def reversed(self) -> bool:
        return self.strand == "-"

    def xlim(self) -> Tuple[int, int]:
        return (self.end, self.start) if self.reversed else (self.start, self.end)


def _as_region(region) -> Optional[Region]:
    if region is None or isinstance(region, Region):
        return region
    if len(region) == 3:
        return Region(*region)
    if len(region) == 4:
        return Region(*region)
    raise TypeError("region must be a Region or (chrom, start, end[, strand])")


# ============================================================================
# Geometry helpers
# ============================================================================
def _apply_xlim(ax, region):
    r = _as_region(region)
    if r is not None:
        ax.set_xlim(*r.xlim())
    return r


def _bp_per_inch(ax) -> float:
    """Data units per inch along x, from the axes' current size and limits."""
    fig = ax.figure
    bbox = ax.get_position()
    width_in = bbox.width * fig.get_figwidth()
    x0, x1 = ax.get_xlim()
    return abs(x1 - x0) / max(width_in, 1e-9)


def _x_direction(ax) -> int:
    """+1 when data x increases to the right, -1 on a reversed axis."""
    x0, x1 = ax.get_xlim()
    return 1 if x1 >= x0 else -1


def _transcription_sign(ax, strand: str) -> int:
    """Screen direction of transcription: +1 = rightwards."""
    return _x_direction(ax) * (1 if strand == "+" else -1)


# ============================================================================
# Layer B: the gene model
# ============================================================================
def gene_model(
    ax,
    tx: Transcript,
    *,
    y: float = 0.0,
    height: float = 1.0,
    region=None,
    label: Union[bool, str] = True,
    label_pos: str = "left",
    intron_style: Optional[str] = None,
    tss: bool = False,
    show_utr: bool = True,
    zorder: int = 3,
    alpha: float = 1.0,
    xform: XForm = None,
) -> Dict[str, list]:
    """Draw one transcript as exon boxes joined by an intron line, in the grey ramp.

    Takes NO colour argument (layer-B rule). ``height`` is the full exon height in
    data units; UTRs draw at ``geometry.track.utr_height_ratio`` of it. ``intron_style``
    is ``"line"`` or ``"chevron"`` (chevrons point in the direction of transcription and
    are spaced in inches, so set the axes limits -- pass ``region`` -- before calling).
    ``tss=True`` adds a bent TSS arrow at the 5' end. ``label`` is the gene name in
    italics (a string overrides ``tx.label``); ``label_pos`` is ``"left"``, ``"right"``,
    ``"above"`` or ``"inside"``. ``xform`` maps genome coordinates to display coordinates
    (e.g. ``pileup.SplicedAxis(...).to_t``) so the model draws in union-exon space.

    Returns the artists by kind, for callers that need to restyle or exempt them.
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for gene_model")
    G = TOK.track_geometry()
    S = TOK.stroke()
    T = TOK.typography()
    ink, hairline = TOK.color("ink"), TOK.color("hairline")
    if xform is None:
        _apply_xlim(ax, region)
    X = xform or (lambda p: p)
    style = intron_style or G.get("intron_style", "line")
    out: Dict[str, list] = {"exons": [], "utrs": [], "introns": [], "chevrons": [], "tss": [], "label": []}

    # exons (CDS at full height, UTR at the ratio)
    for (s, e), is_cds in tx.utr_cds_parts():
        h = height if (is_cds or not show_utr) else height * G["utr_height_ratio"]
        xa, xb = sorted((X(s), X(e)))
        patch = Rectangle((xa, y - h / 2), xb - xa, h, facecolor=ink, edgecolor="none",
                          alpha=alpha, zorder=zorder)
        ax.add_patch(patch)
        (out["exons"] if is_cds else out["utrs"]).append(patch)

    # introns
    for s, e in tx.introns:
        xa, xb = X(s), X(e)
        (line,) = ax.plot([xa, xb], [y, y], color=hairline, lw=S["hairline"], solid_capstyle="butt",
                          zorder=zorder - 1, alpha=alpha)
        out["introns"].append(line)
        if style == "chevron":
            sign = _transcription_sign(ax, tx.strand)
            spacing = G["chevron_spacing_in"] * _bp_per_inch(ax)
            n = int(abs(xb - xa) // spacing)
            if n >= 1:
                xs = xa + (xb - xa) / (n + 1) * np.arange(1, n + 1)
                (chev,) = ax.plot(xs, np.full(n, y), linestyle="none",
                                  marker=">" if sign > 0 else "<", ms=G["mark_size_pt"] * 0.8,
                                  color=hairline, mec=hairline, zorder=zorder - 1, alpha=alpha)
                out["chevrons"].append(chev)

    # TSS: a bent arrow at the 5' end, rising from the exon and pointing with transcription
    if tss:
        sign = _transcription_sign(ax, tx.strand)
        x0, ytop = X(tx.tss), y + height * 0.5
        rise_pt = G["tss_arrow_in"] * 72 * 0.9
        head_pt = G["tss_arrow_in"] * 72
        stem = ax.annotate("", xy=(x0, ytop), xytext=(0, rise_pt), textcoords="offset points",
                           arrowprops=dict(arrowstyle="-", color=ink, lw=S["secondary"],
                                           shrinkA=0, shrinkB=0), zorder=zorder + 1)
        head = ax.annotate("", xy=(x0, ytop), xytext=(sign * head_pt, rise_pt),
                           textcoords="offset points",
                           arrowprops=dict(arrowstyle="<|-", color=ink, lw=S["secondary"],
                                           mutation_scale=8, shrinkA=0, shrinkB=0,
                                           connectionstyle="angle,angleA=0,angleB=90,rad=0"),
                           zorder=zorder + 1)
        out["tss"] += [stem, head]

    # label: italic gene name
    if label:
        text = tx.label if label is True else str(label)
        xdir = _x_direction(ax)
        xs_, xe_ = X(tx.start), X(tx.end)
        if label_pos == "inside":
            t = ax.text((xs_ + xe_) / 2, y, text, ha="center", va="center",
                        fontsize=T["in_figure"], fontstyle="italic",
                        color=TOK.color("paper"), zorder=zorder + 2)
        elif label_pos == "above":
            lift = 2.5 + (G["tss_arrow_in"] * 72 * 0.9 + 1.0 if tss else 0.0)
            t = ax.annotate(text, xy=((xs_ + xe_) / 2, y + height / 2), xytext=(0, lift),
                            textcoords="offset points", ha="center", va="bottom",
                            fontsize=T["in_figure"], fontstyle="italic", color=ink, zorder=zorder + 2)
        else:
            left = (label_pos == "left")
            # "left"/"right" are SCREEN sides; on a reversed axis the data end flips
            xdata = (min if (left == (xdir > 0)) else max)(xs_, xe_)
            t = ax.annotate(text, xy=(xdata, y), xytext=(-4 if left else 4, 0),
                            textcoords="offset points", ha="right" if left else "left",
                            va="center", fontsize=T["in_figure"], fontstyle="italic",
                            color=ink, zorder=zorder + 2)
        out["label"].append(t)
    return out


def gene_track(
    ax,
    transcripts: Sequence[Transcript],
    region=None,
    *,
    height: Optional[float] = None,
    min_gap_bp: Optional[int] = None,
    label: bool = True,
    label_pos: str = "above",
    intron_style: Optional[str] = None,
    tss: bool = False,
    show_utr: bool = True,
    priority: Sequence[str] = (),
    row_pitch: Optional[float] = None,
    xform: XForm = None,
) -> Dict[str, float]:
    """Stack several transcripts in one axes without overlaps (greedy levels), draw each
    with :func:`gene_model`, and hide the y-axis.

    Returns ``{name: y}`` -- the y at which each transcript was drawn (level 0 is y=0, the
    next row y=-row_pitch, ...), which is what :func:`mark` needs as its ``y=``. Transcripts
    named in ``priority`` are placed first, so the locus the panel is about takes the top row.
    ``row_pitch`` defaults to 1.0, or 1.75 when ``label_pos="above"`` -- a label above a lower
    row otherwise rises into the row above it (found on the 877 reference sheet).
    Transcripts outside ``region`` are skipped; those straddling it are clipped by the axes.
    """
    G = TOK.track_geometry()
    r = _as_region(region) if xform is not None else _apply_xlim(ax, region)
    height = G["exon_height"] if height is None else height
    min_gap = G["min_gap_bp"] if min_gap_bp is None else min_gap_bp
    keep = [t for t in transcripts
            if r is None or (t.end > r.start and t.start < r.end)]
    pri = {n: i for i, n in enumerate(priority)}
    # priority transcripts claim rows first (assign_feature_levels is greedy in list order
    # AFTER its own sort by start, so seed their rows explicitly)
    levels: Dict[str, int] = {}
    ends: List[int] = []
    ordered = sorted(keep, key=lambda t: (pri.get(t.name, len(pri)), t.start))
    for t in ordered:
        placed = False
        for lv, e in enumerate(ends):
            if t.start > e + min_gap:
                levels[t.name] = lv
                ends[lv] = t.end
                placed = True
                break
        if not placed:
            levels[t.name] = len(ends)
            ends.append(t.end)
    # a label above a lower row rises into the row above it; a TSS arrow lifts the label further
    # (877 sheet, audit 879 F4.8) -- both are paid for in the default pitch
    if row_pitch is not None:
        pitch = row_pitch
    elif label and label_pos == "above":
        pitch = 1.75 + (0.9 if tss else 0.0)
    else:
        pitch = 1.0 + (0.5 if tss else 0.0)
    ys: Dict[str, float] = {}
    for t in keep:
        y = -pitch * levels[t.name] if levels[t.name] else 0.0
        gene_model(ax, t, y=y, height=height, label=label, label_pos=label_pos,
                   intron_style=intron_style, tss=tss, show_utr=show_utr, xform=xform)
        ys[t.name] = y
    n_levels = (max(levels.values()) + 1) if levels else 1
    ax.set_ylim(-pitch * (n_levels - 1) - 0.9, 0.9 if label_pos != "above" else 1.6)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    return ys


# ============================================================================
# Layer B marks -- three shapes, three identities
# ============================================================================
#: kind -> (token, shape). ``pas`` is the hexamer position (open circle) as opposed to the
#: cleavage site itself (filled circle); both are the poly(A) identity.
MARK_KINDS: Dict[str, Tuple[str, str]] = {
    "tss": ("ink", "arrow"),
    "ss5": ("splice", "tick"),
    "ss3": ("splice", "tick"),
    "bp": ("splice", "dot"),
    "polya": ("polya", "lollipop_o"),
    "cpa": ("polya", "lollipop_o"),
    "pas": ("polya", "lollipop_open"),
    "mod": ("mod", "lollipop_d"),
}


def mark(
    ax,
    x: float,
    kind: str,
    *,
    y: float = 0.0,
    height: float = 1.0,
    label: Optional[str] = None,
    side: str = "above",
    strand: str = "+",
    zorder: int = 6,
    xform: XForm = None,
    ha: str = "center",
) -> list:
    """A site mark on a track at data ``x``: a molecular identity, so NO colour argument.

    ``kind`` is one of :data:`MARK_KINDS`. ``y``/``height`` locate the track the mark sits on
    (the same values passed to :func:`gene_model`). Lollipop stems are in points, so the
    mark keeps its size when the axes is resized. ``side`` is ``"above"`` or ``"below"``.
    ``strand`` orients the splice-site triangle so it points INTO the intron. ``ha`` places the
    label centred on the mark, or hanging ``"left"`` / ``"right"`` of it (two long labels 500 bp
    apart cannot both be centred -- the 176 rebuild had to hand-place them; audit 879 F4.10).
    """
    if kind not in MARK_KINDS:
        raise ValueError(f"unknown mark kind {kind!r}; kinds are {sorted(MARK_KINDS)}")
    token, shape = MARK_KINDS[kind]
    col = TOK.color(token)
    if xform is not None:
        x = xform(x)
    G = TOK.track_geometry()
    S = TOK.stroke()
    T = TOK.typography()
    sgn = 1 if side == "above" else -1
    edge = y + sgn * height / 2
    ms = G["mark_size_pt"]
    stem_pt = G["lollipop_stem_in"] * 72
    out: list = []

    if shape == "tick":
        (ln,) = ax.plot([x, x], [y - height * 0.45, y + height * 0.45], color=col,
                        lw=S["secondary"], solid_capstyle="butt", zorder=zorder)
        out.append(ln)
        # triangle at the tick's outer end pointing into the intron
        into = _x_direction(ax) * (1 if strand == "+" else -1)
        if kind == "ss3":
            into = -into
        (tri,) = ax.plot([x], [edge], linestyle="none", marker=">" if into > 0 else "<",
                         ms=ms * 0.9, color=col, mec=col, zorder=zorder,
                         transform=offset_copy(ax.transData, ax.figure, 0, sgn * 1.5, units="points"))
        out.append(tri)
        text_anchor, text_off = (x, edge), sgn * (ms * 0.9 + 2.5)
    elif shape == "dot":
        (dot,) = ax.plot([x], [y], linestyle="none", marker="o", ms=ms, color=col, mec=col,
                         zorder=zorder)
        out.append(dot)
        text_anchor, text_off = (x, edge), sgn * 2.5
    elif shape.startswith("lollipop"):
        stem = ax.annotate("", xy=(x, edge), xytext=(0, sgn * stem_pt), textcoords="offset points",
                           arrowprops=dict(arrowstyle="-", color=col, lw=S["secondary"],
                                           shrinkA=0, shrinkB=0), zorder=zorder)
        marker = {"lollipop_o": "o", "lollipop_open": "o", "lollipop_d": "D"}[shape]
        face = TOK.color("paper") if shape == "lollipop_open" else col
        (head,) = ax.plot([x], [edge], linestyle="none", marker=marker, ms=ms, mfc=face, mec=col,
                          mew=S["secondary"], zorder=zorder + 1,
                          transform=offset_copy(ax.transData, ax.figure, 0, sgn * stem_pt, units="points"))
        out += [stem, head]
        text_anchor, text_off = (x, edge), sgn * (stem_pt + ms * 0.6 + 2.0)
    elif shape == "arrow":
        into = _x_direction(ax) * (1 if strand == "+" else -1)
        rise = G["tss_arrow_in"] * 72 * 0.9
        head_pt = G["tss_arrow_in"] * 72
        stem = ax.annotate("", xy=(x, edge), xytext=(0, sgn * rise), textcoords="offset points",
                           arrowprops=dict(arrowstyle="-", color=col, lw=S["secondary"],
                                           shrinkA=0, shrinkB=0), zorder=zorder)
        head = ax.annotate("", xy=(x, edge), xytext=(into * head_pt, sgn * rise),
                           textcoords="offset points",
                           arrowprops=dict(arrowstyle="<|-", color=col, lw=S["secondary"],
                                           mutation_scale=8, shrinkA=0, shrinkB=0,
                                           connectionstyle="angle,angleA=0,angleB=90,rad=0"),
                           zorder=zorder)
        out += [stem, head]
        text_anchor, text_off = (x, edge), sgn * (rise + 3.0)
    else:  # pragma: no cover
        raise AssertionError(shape)

    if label:
        if ha not in ("center", "left", "right"):
            raise ValueError("ha must be 'center', 'left' or 'right'")
        dx = 0 if ha == "center" else (ms * 0.6 + 2) * (1 if ha == "left" else -1)
        t = ax.annotate(label, xy=text_anchor, xytext=(dx, text_off), textcoords="offset points",
                        ha=ha, va="bottom" if sgn > 0 else "top", fontsize=T["annotation"],
                        color=col, zorder=zorder + 2)
        out.append(t)
    return out


# ============================================================================
# Signal glyphs -- data about a sample, so they take a ROLE
# ============================================================================
def _signal_color(role: Optional[str]) -> str:
    if role is None:
        return TOK.color("subtle")
    return TOK.role(role)


def coverage(
    ax,
    positions,
    depths,
    *,
    role: Optional[str] = None,
    region=None,
    strand: str = "+",
    mirror: bool = False,
    fill_alpha: Optional[float] = None,
    outline: bool = False,
    label: Optional[str] = None,
    normalize: bool = False,
    log: bool = False,
    smooth: Optional[int] = None,
    zorder: int = 2,
    xform: XForm = None,
):
    """A coverage / signal fill for ONE sample.

    ``role`` is a layer-A argument role (``"focal"``, ``"reference"``, ``"stratum_a"``, ...) or
    ``"neutral"``; with no role the fill is ``subtle`` grey. A layer-B token is refused.
    Draw one sample per axes row and label the row: two samples overlaid on one axes must
    take ``stratum_a`` / ``stratum_b`` (hatched by ``panels.bars``), never focal / reference,
    because those two are dL* 8 apart and separate in mono print only by position.
    ``strand="-"`` with ``mirror=True`` draws the signal below the axis.
    """
    G = TOK.track_geometry()
    col = _signal_color(role)
    if xform is None:
        _apply_xlim(ax, region)
    pos = np.asarray(positions)
    if xform is not None:
        pos = np.array([xform(float(p)) for p in pos])
    dep = np.asarray(depths, dtype=float)
    if smooth is not None and smooth > 1:
        dep = np.convolve(dep, np.ones(smooth) / smooth, mode="same")
    if normalize and dep.max() > 0:
        dep = dep / dep.max()
    if log:
        dep = np.log10(dep + 1)
    if strand == "-" and mirror:
        dep = -dep
    alpha = G["coverage_fill_alpha"] if fill_alpha is None else fill_alpha
    poly = ax.fill_between(pos, 0, dep, color=col, alpha=alpha, lw=0, label=label, zorder=zorder)
    out = [poly]
    if outline:
        (ln,) = ax.plot(pos, dep, color=col, lw=TOK.stroke()["hairline"], zorder=zorder + 1)
        out.append(ln)
    return out


def strand_coverage(
    ax,
    positions_plus, depths_plus, positions_minus, depths_minus,
    *,
    role: Optional[str] = None,
    region=None,
    fill_alpha: Optional[float] = None,
    labels: Tuple[str, str] = ("+ strand", "− strand"),
    xform: XForm = None,
):
    """Plus above, minus below the axis, in ONE colour: strand is position, not hue."""
    out = coverage(ax, positions_plus, depths_plus, role=role, region=region, strand="+",
                   fill_alpha=fill_alpha, label=labels[0], xform=xform)
    out += coverage(ax, positions_minus, depths_minus, role=role, region=region, strand="-",
                    mirror=True, fill_alpha=fill_alpha, label=labels[1], xform=xform)
    ax.axhline(0, color=TOK.color("hairline"), lw=TOK.stroke()["hairline"], zorder=3)
    return out


def arc(
    ax,
    x0: float,
    x1: float,
    *,
    y: float = 0.0,
    height: Optional[float] = None,
    role: Optional[str] = None,
    lw: Optional[float] = None,
    label: Optional[str] = None,
    direction: str = "up",
    alpha: float = 1.0,
    zorder: int = 5,
    xform: XForm = None,
    on_signal: bool = False,
) -> list:
    """A junction arc from ``x0`` to ``x1`` with its apex ``height`` data units above ``y``.

    With no ``role`` it is an ANNOTATION arc in layer-B ``splice``. With a role it is a
    SIGNAL arc (junction reads counted per sample) in that layer-A colour; scale ``lw`` with
    the count. ``on_signal=True`` is for an arc drawn OVER its own sample's coverage: it is
    ``ink`` at hairline weight with an ``ink`` count, whatever the role -- a signal glyph never
    overprints the signal it summarizes in the same colour (audit 879 F4.2; SKILL.md section 7).
    ``height`` defaults to ``geometry.track.arc_height_ratio`` of the y-range.
    """
    G = TOK.track_geometry()
    S = TOK.stroke()
    col = TOK.color("splice") if role is None else TOK.role(role)
    if on_signal:
        col = TOK.color("ink")
        lw = S["hairline"] if lw is None else lw
    if xform is not None:
        x0, x1 = xform(x0), xform(x1)
    if height is None:
        lo, hi = ax.get_ylim()
        height = G["arc_height_ratio"] * abs(hi - lo)
    sgn = 1 if direction == "up" else -1
    xm = 0.5 * (x0 + x1)
    verts = [(x0, y), (xm, y + sgn * 2 * height), (x1, y)]   # quadratic: apex = y + height
    codes = [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]
    patch = PathPatch(MPath(verts, codes), fill=False, ec=col, lw=S["secondary"] if lw is None else lw,
                      joinstyle="round", capstyle="round", alpha=alpha, zorder=zorder)
    ax.add_patch(patch)
    out = [patch]
    if label is not None:
        lw_pt = S["secondary"] if lw is None else lw
        t = ax.annotate(str(label), xy=(xm, y + sgn * height), xytext=(0, sgn * (lw_pt / 2 + 2.5)),
                        textcoords="offset points", ha="center",
                        va="bottom" if sgn > 0 else "top", fontsize=TOK.typography()["annotation"],
                        color=col, zorder=zorder + 1)
        out.append(t)
    return out


def reads(
    ax,
    starts,
    ends,
    *,
    junctions: Optional[List[List[Interval]]] = None,
    roles: Optional[Sequence[Optional[str]]] = None,
    region=None,
    read_lw: float = 4.0,
    alpha: float = 0.9,
    xform: XForm = None,
) -> int:
    """Stacked individual reads (a signal glyph). ``roles`` is one layer-A role (or ``None``
    for ``subtle``) per read; wraps :mod:`read_browser`. ``xform`` draws in spliced space
    (starts, ends and junctions are all mapped). Returns the number of rows."""
    from .read_browser import assign_rows, draw_stacked_reads
    if xform is None:
        _apply_xlim(ax, region)
    starts = np.asarray(starts)
    ends = np.asarray(ends)
    if xform is not None:
        starts = np.array([xform(float(v)) for v in starts]); ends = np.array([xform(float(v)) for v in ends])
        starts, ends = np.minimum(starts, ends), np.maximum(starts, ends)
        if junctions is not None:
            junctions = [[(xform(float(a)), xform(float(b))) for a, b in js] for js in junctions]
    if roles is None:
        colors = [TOK.color("subtle")] * len(starts)
    else:
        colors = [_signal_color(r) for r in roles]
    rows = assign_rows(starts, ends)
    draw_stacked_reads(ax, starts, ends, rows, colors, junction_lists=junctions,
                       read_lw=read_lw, intron_lw=TOK.stroke()["hairline"], alpha=alpha)
    n = int(rows.max()) + 1 if len(rows) else 0
    ax.set_ylim(-0.5, max(n, 1) - 0.5)
    ax.set_yticks([])
    return n


# ============================================================================
# The genomic axis
# ============================================================================
def region_axis(ax, region, *, unit: str = "auto", label: bool = True, n_ticks: int = 5):
    """Ticks in kb (or bp / Mb) with the chromosome in the label; reversed for a
    minus-strand region so the locus reads 5' -> 3' left to right."""
    from .figure_utils import format_genomic_axis
    r = _as_region(region)
    format_genomic_axis(ax, r.start, r.end, show_unit=label, n_ticks=n_ticks, unit=unit)
    if r.reversed:
        ax.set_xlim(r.end, r.start)
    if label:
        u = ax.get_xlabel()
        u = u[u.find("(") + 1:u.find(")")] if "(" in u else unit
        ax.set_xlabel(f"{r.chrom} position ({u})" + ("  5′→3′ (− strand)" if r.reversed else ""))
    return ax
