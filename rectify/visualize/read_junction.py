"""
One read, one junction, every arm, at BASE resolution -- the per-read review PNG.

Kevin reviews individual reads on his phone (the CLAUDE.md rule: vet alignments at the
individual-read level at every stage). For a read in a review bundle (``stock.bam`` plus one
``<sha>.bam`` per arm, ``manifest.tsv``) this module draws, per junction the arms disagree
about, ONE panel with the same x-window for every arm:

* the reference letters (RNA sense; a minus-strand read is drawn 5' -> 3' left to right with
  complemented letters), the coordinates of the junction ends;
* the annotation model -- union exons in ``ink``, the intron as a ``hairline`` with chevrons
  in the direction of transcription, every annotated exon end inside the window as a
  ``splice`` tick;
* one row per arm: the read as LETTERS where it is aligned (``subtle`` = match, an ``ink``
  letter with a one-column tick = mismatch), a NOTCHED body for a deletion, a RAISED
  half-height block for an insertion, the soft clip as hatched letters continuing ungapped
  from the aligned end, the N-op as a thin ``hairline`` connector across the intron gap;
* a numbers-only verdict strip per arm (Kevin adjudicates; no verdict words): that arm's
  junction nearest the panel and its motif, the 5' block's matched/aligned and identity, its
  I/D, the junction-proximal clean run on each side, the 5' clip's best ungapped fit to an
  annotated exon end within 5 kb, MAPQ.

A junction is drawn as two base-resolution FLANKS (the donor end and the acceptor end, each
``window`` letters wide, centered on the intron boundary) joined by a short gap, so a 100-kb
intron and a 90-bp intron get the same picture; an intron shorter than the two flanks is
drawn continuously.

The data model is ``dev/todo_run_20260905/replay_scripts/inspect_read.py`` (blocks, introns
with motif + annotation status, the 5' block base by base, the clip's ungapped placement
against annotated exon ends); its logic is ported here as functions over a genome accessor
because the script hard-codes the Sherlock FASTA path and runs at import.

COLOURS, TYPE AND STROKES come from the rna-figure tokens. ``rectify.visualize.tokens`` is
used when this tree has it; a tree that predates the token layer reads ``$RNA_FIGURE_TOKENS``
or ``~/.claude/skills/rna-figure/tokens.json`` (the same file). The exon / intron glyph is
drawn here with the tokens' track geometry (``tracks.gene_model``'s conventions: ``ink`` exon,
``hairline`` intron, chevrons; no colour argument anywhere) because ``tracks`` is not on every
branch this module has to run on. Per rna-figure section 7: bases are letters in ``ink``, a
mismatch is one ``ink`` tick, strand is position, never a hue.

WINDOW. The house type floor (7.5 pt = 1.89 mm cap) caps a 7.2-in figure at ~75 letter
columns, so ``window`` is the number of letters per junction END (default 32: 16 each side of
the donor boundary and 16 each side of the acceptor boundary, 64 columns + the gap).

Usage::

    python -m rectify.visualize.read_junction --bundle DIR --read 26f8fb45 \\
        --genome GRCh38_primary.fa --gtf gencode.v48.basic.annotation.gtf --out 26f8fb45.png
    # --all-junctions: one panel per junction present in ANY arm (default: the junction(s)
    # at which the arms disagree, else the 5'-most)

``--genome`` is any ``pysam.FastaFile``; a slice FASTA whose contigs are named
``chrX:12969411-12983142`` (the ``offset0based=`` header convention of the d4 replay slices)
is mapped back to genome coordinates automatically.

Author: Kevin R. Roy
"""

from __future__ import annotations

import argparse
import collections
import datetime
import glob
import json
import os
import re
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = ["render_read", "Genome", "Annotation", "ArmView", "arm_view", "select_junctions",
           "check_floors", "apply_style", "main"]

# ----------------------------------------------------------------------------- tokens
_TOK: Optional[dict] = None
_SKILL_DIR = os.path.expanduser("~/.claude/skills/rna-figure")


def tokens() -> dict:
    """The rna-figure token file: ``rectify.visualize.tokens`` when this tree has it, else
    ``$RNA_FIGURE_TOKENS`` / the skill's master copy. Raises with the remedy when neither exists."""
    global _TOK
    if _TOK is not None:
        return _TOK
    try:
        from . import tokens as _t  # the token layer (feat/rna-figure-871 and later)
        _TOK = _t.load()
        return _TOK
    except ImportError:
        pass
    path = os.environ.get("RNA_FIGURE_TOKENS") or os.path.join(_SKILL_DIR, "tokens.json")
    if not os.path.exists(path):
        raise ImportError("read_junction needs the rna-figure tokens: this rectify tree has no "
                          "rectify.visualize.tokens and %s does not exist; set RNA_FIGURE_TOKENS to a copy "
                          "of ~/.claude/skills/rna-figure/tokens.json" % path)
    with open(path) as fh:
        _TOK = json.load(fh)
    return _TOK


def color(name: str) -> str:
    t = tokens()
    pal = {**t["palette"]["semantic"], **t["palette"]["structural"]}
    return pal[name]


def _type() -> Dict[str, float]:
    return dict(tokens()["typography"]["manuscript_pt"])


def _stroke() -> Dict[str, float]:
    return dict(tokens()["geometry"]["stroke_pt"])


def _skill_panels():
    """The skill's ``panels`` module (``apply`` / ``check_floors``) when the skill directory exists."""
    d = os.environ.get("RNA_FIGURE_SKILL") or _SKILL_DIR
    if not os.path.exists(os.path.join(d, "panels.py")):
        return None
    if d not in sys.path:
        sys.path.insert(0, d)
    try:
        import panels  # type: ignore
        return panels
    except Exception:
        return None


def apply_style(matplotlib_module) -> None:
    """The house rcParams (``panels.apply`` when the skill is present, else the same values from the tokens)."""
    P = _skill_panels()
    if P is not None:
        P.apply(matplotlib_module)
        return
    t = tokens()
    T, S = _type(), _stroke()
    fam = t["typography"]["family"]
    matplotlib_module.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": fam, "font.size": T["in_figure"],
        "axes.labelsize": T["axis_label"], "xtick.labelsize": T["tick_label"], "ytick.labelsize": T["tick_label"],
        "legend.fontsize": T["annotation"], "axes.linewidth": S["hairline"], "lines.linewidth": S["secondary"],
        "patch.linewidth": S["hairline"], "text.color": color("ink"), "axes.edgecolor": color("hairline"),
        "savefig.bbox": None, "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def check_floors(fig, **kw) -> bool:
    """Type / stroke / on-page floors. Delegates to the skill's ``panels.check_floors``; without the
    skill, a local check of the same three rules against the tokens."""
    P = _skill_panels()
    if P is not None:
        return P.check_floors(fig, **kw)
    import numpy as np
    t = tokens()
    cap_min, ratio = t["typography"]["cap_min_mm"], t["typography"]["cap_height_ratio"]
    stroke_min = t["geometry"]["stroke_min_pt"]
    bad = []
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    texts = list(fig.texts)
    for ax in fig.axes:
        texts += list(ax.texts)
        for a in list(ax.patches) + list(ax.lines):
            lw = float(np.atleast_1d(a.get_linewidth())[0])
            if 0 < lw < stroke_min - 1e-9:
                bad.append(("stroke", type(a).__name__, lw))
    for tx in texts:
        if not tx.get_text().strip() or not tx.get_visible():
            continue
        if tx.get_fontsize() * ratio * 25.4 / 72.0 < cap_min - 1e-3:
            bad.append(("type", tx.get_text()[:30], tx.get_fontsize()))
        bb = tx.get_window_extent(r)
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > W + 0.5 or bb.y1 > H + 0.5:
            bad.append(("off-page", tx.get_text()[:30], (bb.x0, bb.y0, bb.x1, bb.y1)))
    if kw.get("verbose", True):
        print("read_junction: floors -- %s" % ("OK" if not bad else "%d BELOW" % len(bad)))
        for b in bad[:12]:
            print("   ", b)
    return not bad


# ----------------------------------------------------------------------------- genome + annotation
_SLICE_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<s>\d+)-(?P<e>\d+)$")


class Genome:
    """A ``pysam.FastaFile`` that also understands slice contigs named ``chrom:start-end``
    (0-based half-open, the ``offset0based=`` header convention). ``fetch`` returns ``None``
    when the interval is not covered (a slice file), so callers can skip instead of failing."""

    def __init__(self, path):
        import pysam
        self.path = str(path)
        self.fa = pysam.FastaFile(self.path)
        self.slices: Dict[str, List[Tuple[int, int, str]]] = collections.defaultdict(list)
        self.direct = set()
        for ref, ln in zip(self.fa.references, self.fa.lengths):
            m = _SLICE_RE.match(ref)
            if m and int(m["e"]) - int(m["s"]) == ln:
                self.slices[m["chrom"]].append((int(m["s"]), int(m["e"]), ref))
            else:
                self.direct.add(ref)

    def fetch(self, chrom: str, start: int, end: int) -> Optional[str]:
        if start < 0 or end <= start:
            return None
        if chrom in self.direct:
            if end > self.fa.get_reference_length(chrom):
                return None
            return self.fa.fetch(chrom, start, end).upper()
        for s, e, ref in self.slices.get(chrom, ()):
            if s <= start and end <= e:
                return self.fa.fetch(ref, start - s, end - s).upper()
        return None


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _comp(seq: str) -> str:
    return seq.translate(_RC)


class Annotation:
    """Exons of a GTF (``exon`` lines, ``transcript_id``), per chromosome: annotated introns,
    exon starts / ends, and the union exon intervals (for the model row)."""

    def __init__(self, gtf, chroms: Optional[Sequence[str]] = None):
        self.path = str(gtf)
        ex = collections.defaultdict(list)
        self.gene_of: Dict[str, str] = {}
        want = set(chroms) if chroms else None
        with open(self.path) as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "exon":
                    continue
                if want is not None and f[0] not in want:
                    continue
                m = re.search(r'transcript_id "([^"]+)"', f[8])
                tid = m.group(1) if m else f[8]
                g = re.search(r'gene_name "([^"]+)"', f[8])
                self.gene_of[tid] = g.group(1) if g else tid
                ex[(f[0], tid)].append((int(f[3]) - 1, int(f[4])))
        self.introns: Dict[str, set] = collections.defaultdict(set)
        self.exon_starts: Dict[str, set] = collections.defaultdict(set)
        self.exon_ends: Dict[str, set] = collections.defaultdict(set)
        self._exons: Dict[str, List[Tuple[int, int]]] = collections.defaultdict(list)
        #: chrom -> [(transcript_id, sorted exons)]
        self.transcripts: Dict[str, List[Tuple[str, List[Tuple[int, int]]]]] = collections.defaultdict(list)
        for (c, tid), e in ex.items():
            e.sort()
            self.transcripts[c].append((tid, e))
            for (s1, e1), (s2, e2) in zip(e, e[1:]):
                self.introns[c].add((e1, s2))
            for s, t in e:
                self.exon_starts[c].add(s)
                self.exon_ends[c].add(t)
                self._exons[c].append((s, t))
        self.union: Dict[str, List[Tuple[int, int]]] = {}
        for c, lst in self._exons.items():
            lst.sort()
            merged: List[List[int]] = []
            for s, t in lst:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], t)
                else:
                    merged.append([s, t])
            self.union[c] = [(s, t) for s, t in merged]

    def models_in(self, chrom: str, segments: Sequence[Tuple[int, int]], junction: Tuple[int, int],
                  max_rows: int = 3) -> List[Tuple[str, List[Tuple[int, int]]]]:
        """Up to ``max_rows`` DISTINCT exon structures inside the frame: transcripts carrying the frame's
        junction first, then those with an exon end inside a segment, then any overlapping one; two
        transcripts identical inside the frame are one row. ``[(label, exons)]``."""
        lo, hi = min(a for a, _ in segments), max(b for _, b in segments)
        ranked = []
        for tid, e in self.transcripts.get(chrom, ()):
            if e[-1][1] <= lo or e[0][0] >= hi:
                continue
            inside = tuple((max(s, a), min(t, b)) for s, t in e for a, b in segments if s < b and t > a)
            if not inside:
                continue
            has_j = any(e1 == junction[0] and s2 == junction[1] for (s1, e1), (s2, e2) in zip(e, e[1:]))
            ends = any(a < x < b for s, t in e for x in (s, t) for a, b in segments)
            ranked.append(((not has_j, not ends, tid), inside, tid))
        ranked.sort()
        out, seen = [], set()
        for _, inside, tid in ranked:
            if inside in seen:
                continue
            seen.add(inside)
            out.append((self.gene_of.get(tid, tid), [(s, t) for s, t in inside]))
            if len(out) >= max_rows:
                break
        return out

    def exonic(self, chrom: str, pos: int) -> bool:
        for s, t in self.union.get(chrom, ()):
            if s <= pos < t:
                return True
            if s > pos:
                break
        return False


# ----------------------------------------------------------------------------- the per-arm data model (ported from inspect_read.py)
@dataclass
class Block:
    rs: int
    re: int
    qs: int
    qe: int
    M: int = 0
    X: int = 0
    I: int = 0
    D: int = 0

    @property
    def aligned(self) -> int:
        return self.M + self.X

    @property
    def identity(self) -> float:
        return 100.0 * self.M / self.aligned if self.aligned else 0.0


@dataclass
class Intron:
    start: int
    end: int
    motif: str
    annotated: bool


@dataclass
class ArmView:
    name: str
    read_id: str
    chrom: str
    strand: str
    mapq: int
    cigar: str
    ref_start: int
    ref_end: int
    lead_clip: int
    trail_clip: int
    blocks: List[Block]
    introns: List[Intron]
    #: genome position -> read base (aligned columns; '-' = deletion)
    aligned: Dict[int, str]
    #: boundary position (the base after the insertion point) -> inserted length
    insertions: Dict[int, int]
    #: clip letters continued ungapped from the aligned ends: genome position -> base
    clip_letters: Dict[int, str]
    #: best ungapped placement of the 5' clip against an annotated exon end within 5 kb
    clip_fit: Optional[dict] = None

    @property
    def five_clip(self) -> int:
        return self.lead_clip if self.strand == "+" else self.trail_clip

    @property
    def five_block(self) -> Optional[Block]:
        if not self.blocks:
            return None
        return self.blocks[0] if self.strand == "+" else self.blocks[-1]

    def junctions(self) -> set:
        return {(i.start, i.end) for i in self.introns}


def _walk(a, genome: Genome):
    """One pass over the CIGAR: blocks with = / X / I / D counts (inspect_read.blocks_of), the
    aligned column map, insertions and the clip letters."""
    q = a.query_sequence or ""
    blocks: List[Block] = []
    aligned: Dict[int, str] = {}
    ins: Dict[int, int] = {}
    p, qi, cur = a.reference_start, 0, None
    ct = a.cigartuples or []
    for op, ln in ct:
        if op in (4, 5):  # soft / hard clip
            if op == 4:
                qi += ln
            continue
        if op == 3:
            if cur:
                blocks.append(cur)
            cur = None
            p += ln
            continue
        if cur is None:
            cur = Block(rs=p, re=p, qs=qi, qe=qi)
        if op in (0, 7, 8):
            ref = genome.fetch(a.reference_name, p, p + ln)
            rd = q[qi:qi + ln].upper()
            for k in range(ln):
                aligned[p + k] = rd[k] if k < len(rd) else "N"
                if ref is not None and k < len(rd):
                    if ref[k] == rd[k]:
                        cur.M += 1
                    else:
                        cur.X += 1
            p += ln
            qi += ln
        elif op == 1:
            cur.I += ln
            ins[p] = ins.get(p, 0) + ln
            qi += ln
        elif op == 2:
            cur.D += ln
            for k in range(ln):
                aligned[p + k] = "-"
            p += ln
        cur.re, cur.qe = p, qi
    if cur:
        blocks.append(cur)
    lead = ct[0][1] if ct and ct[0][0] == 4 else 0
    trail = ct[-1][1] if ct and ct[-1][0] == 4 else 0
    clip_letters: Dict[int, str] = {}
    for k in range(lead):
        clip_letters[a.reference_start - lead + k] = q[k].upper()
    for k in range(trail):
        clip_letters[a.reference_end + k] = q[len(q) - trail + k].upper()
    return blocks, aligned, ins, lead, trail


def _motif(genome: Genome, chrom: str, s: int, e: int, strand: str) -> str:
    dn, ac = genome.fetch(chrom, s, s + 2), genome.fetch(chrom, e - 2, e)
    if dn is None or ac is None:
        return "??-??"
    return f"{dn}-{ac}" if strand == "+" else f"{_rc(ac)}-{_rc(dn)}"


def _ungapped_best(genome: Genome, seq: str, chrom: str, lo: int, hi: int, anchor: str):
    """inspect_read.ungapped_best: the best ungapped placement of ``seq`` with its LAST base ending at
    any position in [lo, hi) (``anchor='end'``) or its FIRST base starting there (``'start'``)."""
    best = (-1, None)
    for x in range(lo, hi):
        ref = genome.fetch(chrom, x - len(seq), x) if anchor == "end" else genome.fetch(chrom, x, x + len(seq))
        if ref is None or len(ref) != len(seq):
            continue
        m = sum(1 for u, v in zip(ref, seq) if u == v)
        if m > best[0]:
            best = (m, x)
    return best


def _clip_fit(a, genome: Genome, ann: Annotation, strand: str, lead: int, trail: int) -> Optional[dict]:
    """The 5' soft clip's best ungapped placement against annotated exon ends within 5 kb upstream
    (plus) / downstream (minus) -- inspect_read's last block, over the candidates it prints."""
    q = a.query_sequence or ""
    five = lead if strand == "+" else trail
    if not five:
        return None
    chrom = a.reference_name
    clip = q[:lead].upper() if strand == "+" else q[len(q) - trail:].upper()
    edge = a.reference_start if strand == "+" else a.reference_end
    if strand == "+":
        cands = sorted(c for c in ann.exon_ends.get(chrom, ()) if 0 < edge - c <= 5000)
        seg = clip[-min(len(clip), 40):]
    else:
        cands = sorted(c for c in ann.exon_starts.get(chrom, ()) if 0 < c - edge <= 5000)
        seg = clip[:min(len(clip), 40)]
    best = None
    for c in cands:
        m, pos = _ungapped_best(genome, seg, chrom, c - 3, c + 4, "end" if strand == "+" else "start")
        if pos is None:
            continue
        cand = dict(exon_end=c, matches=m, n=len(seg), pos=pos, offset=pos - c, clip_len=five)
        if best is None or (m, -abs(pos - c)) > (best["matches"], -abs(best["offset"])):
            best = cand
    return best


def arm_view(name: str, a, genome: Genome, ann: Annotation) -> ArmView:
    """The inspector's numbers for one alignment record."""
    strand = "-" if a.is_reverse else "+"
    chrom = a.reference_name
    blocks, aligned, ins, lead, trail = _walk(a, genome)
    introns = []
    for b1, b2 in zip(blocks, blocks[1:]):
        s, e = b1.re, b2.rs
        introns.append(Intron(s, e, _motif(genome, chrom, s, e, strand), (s, e) in ann.introns.get(chrom, ())))
    q = a.query_sequence or ""
    clip_letters = {}
    for k in range(lead):
        clip_letters[a.reference_start - lead + k] = q[k].upper()
    for k in range(trail):
        clip_letters[a.reference_end + k] = q[len(q) - trail + k].upper()
    return ArmView(name=name, read_id=a.query_name, chrom=chrom, strand=strand, mapq=a.mapping_quality,
                   cigar=a.cigarstring or "", ref_start=a.reference_start, ref_end=a.reference_end,
                   lead_clip=lead, trail_clip=trail, blocks=blocks, introns=introns, aligned=aligned,
                   insertions=ins, clip_letters=clip_letters,
                   clip_fit=_clip_fit(a, genome, ann, strand, lead, trail))


def _manifest_arm_order(bundle_dir) -> List[str]:
    """The arms in the manifest header's ``corrected_junction_<arm>`` column order (stock, baseline, candidate)."""
    path = os.path.join(str(bundle_dir), "manifest.tsv")
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                return [c[len("corrected_junction_"):] for c in line.rstrip("\n").split("\t") if c.startswith("corrected_junction_")]
    return []


def load_arms(bundle_dir, read_id: str, arms: Optional[Sequence[str]] = None):
    """``{arm: AlignedSegment}`` for the first record whose name starts with ``read_id``. Order: ``stock``
    first, then the manifest header's column order (baseline before candidate), then name order."""
    import pysam
    d = str(bundle_dir)
    names = list(arms) if arms else sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(d, "*.bam")))
    order = {n: k for k, n in enumerate(_manifest_arm_order(d))}
    if not arms:
        names = sorted(names, key=lambda x: (x != "stock", order.get(x, len(order)), x))
    out = {}
    for arm in names:
        path = os.path.join(d, arm + ".bam")
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        with pysam.AlignmentFile(path, check_sq=False) as fh:
            for a in fh:
                if a.is_unmapped or a.is_secondary or a.is_supplementary:
                    continue
                if a.query_name.startswith(read_id):
                    out[arm] = a
                    break
    if not out:
        raise KeyError(f"read {read_id!r} not found in any arm of {d}")
    return out


def read_manifest(bundle_dir, read_id: str) -> dict:
    path = os.path.join(str(bundle_dir), "manifest.tsv")
    if not os.path.exists(path):
        return {}
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            if f and f[0].startswith(read_id):
                return dict(zip(header, f))
    return {}


# ----------------------------------------------------------------------------- junction selection
def _five_most(j: Tuple[int, int], strand: str) -> int:
    return j[0] if strand == "+" else -j[1]


def select_junctions(views: Dict[str, ArmView], mode: str = "auto", window: int = 32) -> List[List[Tuple[int, int]]]:
    """The junctions to draw, grouped into panels. ``auto``: the junctions at which the arms DISAGREE
    (present in some arms, absent from others), else the 5'-most junction; ``all``: every junction in
    any arm. Junctions whose donor and acceptor ends both fall inside one another's flanks share a
    panel (a 4-nt slide is one picture, not two). Each group is sorted with the annotated junction
    first, so it frames the panel."""
    if not views:
        return []
    strand = next(iter(views.values())).strand
    sets = {n: v.junctions() for n, v in views.items()}
    union = set().union(*sets.values())
    if mode == "all":
        chosen = set(union)
    else:
        common = set.intersection(*sets.values()) if sets else set()
        chosen = union - common
        if not chosen and union:
            chosen = {min(union, key=lambda j: _five_most(j, strand))}
    if not chosen:
        return []
    half = max(1, window // 2)
    annotated = {(i.start, i.end) for v in views.values() for i in v.introns if i.annotated}
    ordered = sorted(chosen, key=lambda j: _five_most(j, strand))
    groups: List[List[Tuple[int, int]]] = []
    for j in ordered:
        for g in groups:
            if any(abs(j[0] - k[0]) < half and abs(j[1] - k[1]) < half for k in g):
                g.append(j)
                break
        else:
            groups.append([j])
    for g in groups:
        g.sort(key=lambda j: (j not in annotated, _five_most(j, strand)))
    return groups


# ----------------------------------------------------------------------------- the frame (genome position -> letter column)
@dataclass
class Frame:
    chrom: str
    strand: str
    junction: Tuple[int, int]
    half: int
    gap: int = 3
    segments: List[Tuple[int, int]] = field(default_factory=list)
    offsets: List[int] = field(default_factory=list)

    def __post_init__(self):
        s, e = self.junction
        h = self.half
        if e - s <= 2 * h:  # a short intron: one continuous window
            segs = [(s - h, e + h)]
        else:
            segs = [(s - h, s + h), (e - h, e + h)]
        if self.strand == "-":
            segs = segs[::-1]
        self.segments = segs
        off, self.offsets = 0, []
        for a, b in segs:
            self.offsets.append(off)
            off += (b - a) + self.gap
        self.ncols = off - self.gap

    def col(self, pos: int) -> Optional[int]:
        """The column whose LEFT edge is at x=col for genome position ``pos``; None outside the frame."""
        for (a, b), off in zip(self.segments, self.offsets):
            if a <= pos < b:
                return off + (pos - a if self.strand == "+" else b - 1 - pos)
        return None

    def boundary(self, pos: int) -> Optional[float]:
        """x of the boundary between genome positions pos-1 and pos, in RNA order."""
        c = self.col(pos)
        if c is not None:
            return float(c) if self.strand == "+" else float(c + 1)
        c = self.col(pos - 1)
        if c is not None:
            return float(c + 1) if self.strand == "+" else float(c)
        return None

    def boundary_clipped(self, pos: int) -> float:
        """As ``boundary`` but clipped to the frame's edges (for a connector that leaves the frame)."""
        b = self.boundary(pos)
        if b is not None:
            return b
        # before the first segment in RNA order?
        first_a, first_b = self.segments[0]
        before = pos < first_a if self.strand == "+" else pos >= first_b
        return 0.0 if before else float(self.ncols)

    def positions(self):
        for (a, b), off in zip(self.segments, self.offsets):
            for pos in range(a, b):
                yield pos, self.col(pos)

    def gap_ranges(self) -> List[Tuple[int, int]]:
        out = []
        for (a, b), off in list(zip(self.segments, self.offsets))[:-1]:
            out.append((off + (b - a), off + (b - a) + self.gap))
        return out


# ----------------------------------------------------------------------------- numbers for the verdict strip
def _clean_run(view: ArmView, genome: Genome, j: Tuple[int, int], cap: int = 200) -> Tuple[int, int]:
    """Consecutive matching bases walking away from the junction on its 5' side and its 3' side."""
    s, e = j
    chrom = view.chrom

    def run(start: int, step: int) -> int:
        n = 0
        p = start
        while n < cap:
            b = view.aligned.get(p)
            if b is None or b == "-":
                break
            ref = genome.fetch(chrom, p, p + 1)
            if ref is None or ref != b:
                break
            n += 1
            p += step
        return n

    left = run(s - 1, -1)   # the block ending at the donor (genome left)
    right = run(e, +1)      # the block starting at the acceptor (genome right)
    return (left, right) if view.strand == "+" else (right, left)


def _nearest_junction(view: ArmView, j: Tuple[int, int], within: int) -> Optional[Intron]:
    best, dist = None, None
    for i in view.introns:
        d = min(abs(i.start - j[0]), abs(i.end - j[1]))
        if d <= within and (dist is None or d < dist):
            best, dist = i, d
    return best


def _fmt(n: int) -> str:
    return f"{n:,}"


def verdict_line(view: ArmView, genome: Genome, frame: Frame) -> str:
    """Numbers only. ``J start-end motif · 5′ M/aligned id% · I D · run 5′|3′ · clip n fit m/n @pos±off · MQ``."""
    near = _nearest_junction(view, frame.junction, frame.half)
    if near is not None:
        jtxt = f"J {_fmt(near.start)}–{_fmt(near.end)} {near.motif}" + ("" if near.annotated else " novel")
        l5, l3 = _clean_run(view, genome, (near.start, near.end))
        run = f"run {l5}|{l3}"
    else:
        jtxt = "J none"
        run = "run –"
    b = view.five_block
    blk = f"5′ {b.M}/{b.aligned} {b.identity:.0f} % · I{b.I} D{b.D}" if b else "5′ –"
    if view.five_clip:
        f = view.clip_fit
        if f:
            clip = f"clip {view.five_clip} fit {f['matches']}/{f['n']} @{_fmt(f['exon_end'])} ({f['offset']:+d})"
        else:
            clip = f"clip {view.five_clip} fit –"
    else:
        clip = "clip 0"
    return f"{jtxt} · {blk} · {run} · {clip} · MQ {view.mapq}"


# ----------------------------------------------------------------------------- drawing
_ROW = 3.4          # mm per row unit
_GUTTER_COLS = 6.5  # label gutter, in columns
_MARGIN_MM = 4.0


def _draw_panel(ax, frame: Frame, views: Dict[str, ArmView], genome: Genome, ann: Annotation, *, letter_pt: float):
    """One panel: reference, coordinates, model, one read row + verdict per arm. y grows DOWNWARD
    (row units; the axes' y-axis is inverted by the caller)."""
    from matplotlib.patches import Rectangle
    S = _stroke()
    ink, hair, subtle, mute, wash, splice = (color(n) for n in ("ink", "hairline", "subtle", "mute", "wash", "splice"))
    chrom, strand = frame.chrom, frame.strand
    letter_kw = dict(ha="center", va="center", fontsize=letter_pt)
    fam = tokens()["typography"]["family"]
    sym_fam = [fam[0], "DejaVu Sans"]  # for glyphs Arial lacks (the prime, the en dash are fine; ∼ is not used)

    def L(pos_or_base: str, strand_: str) -> str:
        return _comp(pos_or_base) if strand_ == "-" else pos_or_base

    # ---- row 0: coordinates of the junction ends (+ one more label per flank for scale)
    y_coord, y_ref = 0.55, 1.5
    labels = []
    s, e = frame.junction
    for seg_i, (a, b) in enumerate(frame.segments):
        mid = a + (b - a) // 2  # the boundary each flank is centered on (s or e), or the window center
        for pos in (mid - 10, mid, mid + 10) if (b - a) >= 24 else (mid,):
            bx = frame.boundary(pos)
            if bx is None:
                continue
            labels.append((bx, pos))
    for bx, pos in labels:
        ax.plot([bx, bx], [y_coord + 0.35, y_ref - 0.45], color=hair, lw=S["hairline"], solid_capstyle="butt")
        ax.text(bx, y_coord, _fmt(pos), ha="center", va="center", fontsize=_type()["tick_label"], color=hair)

    # ---- row 1: the reference, RNA sense
    for pos, c in frame.positions():
        base = genome.fetch(chrom, pos, pos + 1)
        ax.text(c + 0.5, y_ref, L(base or "N", strand), color=ink, **letter_kw)
    for g0, g1 in frame.gap_ranges():
        ax.text((g0 + g1) / 2, y_ref, "//", color=subtle, ha="center", va="center", fontsize=_type()["annotation"])

    # ---- annotated exon ends inside the frame, as splice ticks on the ruler line under the reference
    for (a, b), off in zip(frame.segments, frame.offsets):
        for pos in range(a, b + 1):
            if pos in ann.exon_ends.get(chrom, ()) or pos in ann.exon_starts.get(chrom, ()):
                bx = frame.boundary(pos)
                if bx is not None:
                    ax.plot([bx, bx], [y_ref + 0.42, y_ref + 0.62], color=splice, lw=S["secondary"],
                            solid_capstyle="butt", zorder=4)

    # ---- the annotation model rows: up to three distinct exon structures inside the frame
    # (the transcript carrying the frame's junction first); exon = ink block, intron = hairline + chevrons
    h_ex = 0.62
    models = ann.models_in(chrom, frame.segments, frame.junction) or [("", [])]
    y_model = 2.1 + _MODEL_PITCH / 2
    for label, exons in models:
        for (a, b), off in zip(frame.segments, frame.offsets):
            x0, x1 = off, off + (b - a)
            ax.plot([x0, x1], [y_model, y_model], color=hair, lw=S["hairline"], solid_capstyle="butt", zorder=1)
            exonic = set()
            for s_, t_ in exons:
                for pos in range(max(s_, a), min(t_, b)):
                    exonic.add(pos)
            runs: List[List[int]] = []
            for pos in range(a, b):
                c = frame.col(pos)
                if pos in exonic:
                    if runs and runs[-1][1] == c:
                        runs[-1][1] = c + 1
                    elif runs and runs[-1][0] == c + 1:
                        runs[-1][0] = c
                    else:
                        runs.append([c, c + 1])
            for c0, c1 in runs:
                ax.add_patch(Rectangle((c0, y_model - h_ex / 2), c1 - c0, h_ex, facecolor=ink, edgecolor="none", zorder=2))
            intronic = sorted(frame.col(p) for p in range(a, b) if p not in exonic)
            for c in intronic[2::5]:
                ax.plot([c + 0.2, c + 0.7, c + 0.2], [y_model - 0.22, y_model, y_model + 0.22], color=hair,
                        lw=S["hairline"], solid_capstyle="butt", zorder=3)
        short = label if len(label) <= 9 else label[:8] + "…"
        ax.text(-_GUTTER_COLS + 0.2, y_model, short, ha="left", va="center", fontsize=_type()["annotation"],
                color=subtle, fontstyle="italic", fontfamily=sym_fam)
        y_model += _MODEL_PITCH
    ax.text(-_GUTTER_COLS + 0.2, y_ref, "ref", ha="left", va="center", fontsize=_type()["annotation"], color=subtle)

    # ---- one row per arm
    y = 2.1 + _MODEL_PITCH * len(models) + 0.55
    body_h = 0.66
    for name, v in views.items():
        y_body = y + 0.45
        top, bot = y_body - body_h / 2, y_body + body_h / 2
        ax.text(-_GUTTER_COLS + 0.2, y_body, name, ha="left", va="center", fontsize=_type()["annotation"], color=ink)
        # body band over contiguous aligned (non-deletion) columns; a deletion = notch (hairline through the gap)
        cols = {}
        for pos, c in frame.positions():
            b = v.aligned.get(pos)
            if b is not None:
                cols[c] = (pos, b)
        run_start = None
        for c in range(frame.ncols + 1):
            here = c in cols and cols[c][1] != "-"
            if here and run_start is None:
                run_start = c
            if (not here) and run_start is not None:
                ax.add_patch(Rectangle((run_start, top), c - run_start, body_h, facecolor=wash, edgecolor="none", zorder=1))
                run_start = None
        for c, (pos, b) in cols.items():
            ref = genome.fetch(chrom, pos, pos + 1) or "N"
            if b == "-":
                ax.plot([c, c + 1], [y_body, y_body], color=ink, lw=S["secondary"], solid_capstyle="butt", zorder=3)
                continue
            if b == ref:
                ax.text(c + 0.5, y_body, L(b, strand), color=subtle, zorder=4, **letter_kw)
            else:
                ax.text(c + 0.5, y_body, L(b, strand), color=ink, fontweight="bold", zorder=4, **letter_kw)
                ax.plot([c + 0.5, c + 0.5], [top - 0.26, top - 0.04], color=ink, lw=S["secondary"], solid_capstyle="butt", zorder=5)
        # insertions: a raised half-height block at the boundary
        for pos, ln in v.insertions.items():
            bx = frame.boundary(pos)
            if bx is None:
                continue
            ax.add_patch(Rectangle((bx - 0.22, top - 0.3), 0.44, 0.3, facecolor=ink, edgecolor="none", zorder=5))
        # soft clip letters, hatched, continued ungapped from the aligned end
        clip_cols = []
        for pos, b in v.clip_letters.items():
            c = frame.col(pos)
            if c is None:
                continue
            clip_cols.append(c)
            ref = genome.fetch(chrom, pos, pos + 1) or "N"
            if b == ref:
                ax.text(c + 0.5, y_body, L(b, strand), color=subtle, zorder=4, **letter_kw)
            else:
                ax.text(c + 0.5, y_body, L(b, strand), color=ink, fontweight="bold", zorder=4, **letter_kw)
                ax.plot([c + 0.5, c + 0.5], [top - 0.26, top - 0.04], color=ink, lw=S["secondary"], solid_capstyle="butt", zorder=5)
        if clip_cols:
            c0, c1 = min(clip_cols), max(clip_cols) + 1
            ax.add_patch(Rectangle((c0, top), c1 - c0, body_h, facecolor="none", edgecolor=mute, hatch="////",
                                   linewidth=0, zorder=1))
        # N-op connectors (any intron of this arm touching the frame)
        for i in v.introns:
            x0 = frame.boundary_clipped(i.start) if strand == "+" else frame.boundary_clipped(i.end)
            x1 = frame.boundary_clipped(i.end) if strand == "+" else frame.boundary_clipped(i.start)
            if x0 == x1:
                continue
            ax.plot([x0, x1], [y_body, y_body], color=hair, lw=S["hairline"], solid_capstyle="butt", zorder=2)
        # the verdict strip
        ax.text(-_GUTTER_COLS + 0.2, y + 1.38, verdict_line(v, genome, frame), ha="left", va="center",
                fontsize=_type()["annotation"], color=ink, fontfamily=sym_fam)
        y += _ARM_PITCH
    return y


_ARM_PITCH = 2.15
_MODEL_PITCH = 0.85


def _rows_per_panel(n_arms: int, n_models: int = 1) -> float:
    return 2.1 + _MODEL_PITCH * max(1, n_models) + 0.55 + _ARM_PITCH * n_arms + 0.2


def render_read(bundle_dir, read_id: str, genome, gtf, out_png, *, arms: Optional[Sequence[str]] = None,
                window: int = 32, junction_index="auto", sidecar: bool = True) -> dict:
    """Render the per-read junction PNG. Returns ``{"png", "panels", "junctions", "arms", "floors_ok"}``.

    ``junction_index``: ``'auto'`` (the junction(s) the arms disagree about, else the 5'-most),
    ``'all'`` (one panel per junction in any arm), or an int index into the ``'all'`` list.
    ``window``: letters per junction end (centered on the boundary); the type floor caps a 7.2-in
    figure at ~75 columns, so > 36 is refused."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if window > 36:
        raise ValueError("window > 36 cannot hold floor-size letters at COL_DOUBLE (7.2 in); use <= 36")
    apply_style(matplotlib)
    T = _type()
    G = genome if isinstance(genome, Genome) else Genome(genome)
    recs = load_arms(bundle_dir, read_id, arms)
    chroms = {a.reference_name for a in recs.values()}
    A = gtf if isinstance(gtf, Annotation) else Annotation(gtf, chroms=chroms)
    views = {n: arm_view(n, a, G, A) for n, a in recs.items()}
    if junction_index == "all":
        groups = select_junctions(views, "all", window)
    elif junction_index == "auto":
        groups = select_junctions(views, "auto", window)
    else:
        allg = select_junctions(views, "all", window)
        groups = [allg[int(junction_index)]]
    if not groups:
        raise ValueError(f"read {read_id}: no junction in any arm to draw")
    man = read_manifest(bundle_dir, read_id)
    v0 = next(iter(views.values()))
    strand = v0.strand
    frames = [Frame(v0.chrom, strand, g[0], max(1, window // 2)) for g in groups]

    # ---- geometry, in mm
    width_in = tokens()["geometry"]["column_in"]["double"]
    width_mm = width_in * 25.4
    n_arms = len(views)
    n_models = max(1, max(len(A.models_in(f.chrom, f.segments, f.junction)) for f in frames))
    panel_rows = _rows_per_panel(n_arms, n_models)
    title_mm = 9.0
    panel_gap_mm = 5.0
    ncols = max(f.ncols for f in frames)
    # legend + footer band (measured on a probe when figstyle is available)
    legend = _legend_text(read_id, man, frames, list(views))
    fs = _figstyle()
    if fs is not None:
        probe = plt.figure(figsize=(width_in, 4.0))
        top = fs.concise_legend(probe, legend, y0=0.0)
        plt.close(probe)
        legend_mm = top * 4.0 * 25.4 + 2.0
    else:
        legend_mm = 4.2 * len(textwrap.wrap(legend, 120)) + 2.0
    footer_mm = 9.0
    panels_mm = len(frames) * panel_rows * _ROW + (len(frames) - 1) * panel_gap_mm
    height_mm = title_mm + panels_mm + legend_mm + footer_mm + 2 * _MARGIN_MM
    fig = plt.figure(figsize=(width_in, height_mm / 25.4))

    # ---- title (two lines: identity; per-arm 5'-most junction + motif from the manifest)
    ident = " · ".join(x for x in (read_id[:8], man.get("library", ""), man.get("class", ""), f"{v0.chrom} {strand}") if x)
    fig.text(_MARGIN_MM / width_mm, 1 - _MARGIN_MM / height_mm, ident, ha="left", va="top", fontsize=T["in_figure"],
             fontweight="bold", color=color("ink"))
    arm_line = " · ".join(_arm_summary(n, v) for n, v in views.items())
    fig.text(_MARGIN_MM / width_mm, 1 - (_MARGIN_MM + 4.2) / height_mm, arm_line, ha="left", va="top",
             fontsize=T["annotation"], color=color("hairline"))

    # ---- panels
    usable_mm = width_mm - 2 * _MARGIN_MM
    letter_pt = T["annotation"]
    axes = []
    y_top_mm = _MARGIN_MM + title_mm
    for k, (frame, g) in enumerate(zip(frames, groups)):
        h_mm = panel_rows * _ROW
        ax = fig.add_axes([_MARGIN_MM / width_mm, 1 - (y_top_mm + h_mm) / height_mm, usable_mm / width_mm, h_mm / height_mm])
        ax.set_xlim(-_GUTTER_COLS, ncols)
        ax.set_ylim(panel_rows, 0)
        ax.axis("off")
        _draw_panel(ax, frame, views, G, A, letter_pt=letter_pt)
        if len(frames) > 1:
            ax.text(-_GUTTER_COLS, 0.55, "abcdefgh"[k], ha="left", va="center", fontsize=T["panel_letter"],
                    fontweight="bold", color=color("ink"))
        axes.append(ax)
        y_top_mm += h_mm + panel_gap_mm

    # ---- legend + provenance footer + stamp
    script = "rectify.visualize.read_junction"
    if fs is not None:
        fs.concise_legend(fig, legend, script=script, y0=(footer_mm + _MARGIN_MM) / height_mm)
    else:
        fig.text(_MARGIN_MM / width_mm, (footer_mm + _MARGIN_MM) / height_mm, "\n".join(textwrap.wrap(legend, 120)),
                 ha="left", va="bottom", fontsize=T["legend"], color="0.15")
        fig.text(0.995, 0.004, f"{script} · {datetime.datetime.now():%Y-%m-%d %H:%M}", ha="right", va="bottom",
                 fontsize=T["provenance"] + 1.5, color="0.5")
    prov = [f"bundle {_short(bundle_dir)} · arms {', '.join(views)} · window {window} letters per junction end · "
            f"{script} · {datetime.date.today():%Y-%m-%d}",
            f"genome {_short(G.path)} · annotation {_short(A.path)}"]
    for k, line in enumerate(prov):
        fig.text(_MARGIN_MM / width_mm, (_MARGIN_MM + 0.5 + 2.6 * (len(prov) - 1 - k)) / height_mm, line, ha="left",
                 va="bottom", fontsize=T["provenance"], color="0.45")

    ok = check_floors(fig)
    out_png = str(out_png)
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300)
    if sidecar:
        _write_sidecar(out_png, read_id, man, views, frames, G, A, bundle_dir, window)
    result = dict(png=out_png, panels=len(frames), junctions=[f.junction for f in frames], arms=list(views),
                  floors_ok=ok, fig=fig)
    return result


def _short(p) -> str:
    p = str(p)
    parts = p.rstrip("/").split("/")
    return "…/" + "/".join(parts[-2:]) if len(parts) > 2 else p


def _arm_summary(name: str, v: ArmView) -> str:
    if not v.introns:
        return f"{name} no junction (clip {v.five_clip})"
    i = v.introns[0] if v.strand == "+" else v.introns[-1]
    return f"{name} {_fmt(i.start)}–{_fmt(i.end)} {i.motif}"


def _legend_text(read_id: str, man: dict, frames: List[Frame], arms: List[str]) -> str:
    n = len(frames)
    parts = [f"Figure | read {read_id[:8]}, {len(arms)} arms, one junction per panel at base resolution."]
    for k, f in enumerate(frames):
        call = f"({'abcdefgh'[k]}) " if n > 1 else ""
        s, e = f.junction
        parts.append(f"{call}{f.chrom}:{_fmt(s)}–{_fmt(e)} ({f.strand}), {f.half} letters each side of each end.")
    parts.append("Grey letters match, dark letters with a tick mismatch; notch = deletion, raised block = insertion, "
                 "hatched = soft clip continued ungapped, thin line = N-op; blue ticks = annotated exon ends. "
                 "Per arm: junction and motif · 5′ block matched/aligned · I, D · clean run 5′|3′ of the junction · "
                 "5′ clip fit to the nearest annotated exon end within 5 kb · MAPQ.")
    return " ".join(parts)


def _figstyle():
    d = os.path.expanduser("~/.claude/skills/figure-qa")
    if not os.path.exists(os.path.join(d, "figstyle.py")):
        return None
    if d not in sys.path:
        sys.path.insert(0, d)
    try:
        import figstyle  # type: ignore
        return figstyle
    except Exception:
        return None


def _write_sidecar(out_png, read_id, man, views, frames, G, A, bundle_dir, window):
    """The figure-standard .md sidecar next to the PNG. Refuses to overwrite a file it did not write."""
    path = os.path.splitext(out_png)[0] + ".md"
    marker = "<!-- rectify.visualize.read_junction sidecar -->"
    if os.path.exists(path):
        with open(path) as fh:
            if fh.readline().strip() != marker:
                raise FileExistsError(f"{path} exists and is not a read_junction sidecar -- refusing to overwrite")
    lines = [marker, f"# {read_id[:8]} — per-read junction PNG", "",
             f"- read: `{read_id}`", f"- bundle: `{bundle_dir}`", f"- genome: `{G.path}`", f"- annotation: `{A.path}`",
             f"- window: {window} letters per junction end", f"- rendered: {datetime.datetime.now():%Y-%m-%d %H:%M}", ""]
    if man:
        lines += ["## manifest row", ""] + [f"- {k}: {v}" for k, v in man.items() if v] + [""]
    lines += ["## panels", ""]
    for k, f in enumerate(frames):
        s, e = f.junction
        lines.append(f"- panel {'abcdefgh'[k]}: {f.chrom}:{s}-{e} ({f.strand}); segments {f.segments}")
    lines += ["", "## arms (the inspector's numbers)", ""]
    for n, v in views.items():
        lines.append(f"### {n}")
        lines.append(f"- {v.chrom}:{v.ref_start}-{v.ref_end} {v.strand} MAPQ {v.mapq} CIGAR `{v.cigar}`")
        lines.append(f"- soft clips lead {v.lead_clip} trail {v.trail_clip} (5′ clip {v.five_clip})")
        for i, b in enumerate(v.blocks):
            lines.append(f"- block {i + 1}: {b.rs}-{b.re} ({b.re - b.rs} bp) ={b.M} X={b.X} I={b.I} D={b.D} identity {b.identity:.0f} %")
        for i in v.introns:
            lines.append(f"- intron {i.start}-{i.end} ({i.end - i.start} bp) motif {i.motif} {'ANNOTATED' if i.annotated else 'novel'}")
        if v.clip_fit:
            f = v.clip_fit
            lines.append(f"- 5′ clip fit: {f['matches']}/{f['n']} vs annotated exon end {f['exon_end']} at {f['pos']} (offset {f['offset']:+d})")
        for f in frames:
            lines.append(f"- verdict ({f.chrom}:{f.junction[0]}-{f.junction[1]}): {verdict_line(v, G, f)}")
        lines.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# ----------------------------------------------------------------------------- CLI
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m rectify.visualize.read_junction",
                                 description="One read, one junction, every arm, at base resolution (phone-review PNG).")
    ap.add_argument("--bundle", required=True, help="review bundle dir: stock.bam + <sha>.bam per arm (+ manifest.tsv)")
    ap.add_argument("--read", required=True, help="read id or prefix")
    ap.add_argument("--genome", required=True, help="FASTA (pysam.FastaFile); slice contigs named chrom:start-end are mapped back")
    ap.add_argument("--gtf", required=True, help="GTF with exon lines (transcript_id)")
    ap.add_argument("--out", required=True, help="output PNG (a .md sidecar is written beside it)")
    ap.add_argument("--arms", nargs="*", default=None, help="arm names (default: every *.bam, stock first)")
    ap.add_argument("--window", type=int, default=32, help="letters per junction end (default 32; max 36)")
    ap.add_argument("--all-junctions", action="store_true", help="one panel per junction present in any arm")
    ap.add_argument("--junction", default=None, help="index into the --all-junctions list (0 = 5'-most)")
    ap.add_argument("--no-sidecar", action="store_true")
    a = ap.parse_args(argv)
    mode = "all" if a.all_junctions else ("auto" if a.junction is None else int(a.junction))
    res = render_read(a.bundle, a.read, a.genome, a.gtf, a.out, arms=a.arms, window=a.window, junction_index=mode,
                      sidecar=not a.no_sidecar)
    print(f"{res['png']}  panels {res['panels']}  junctions {res['junctions']}  arms {res['arms']}  floors {'OK' if res['floors_ok'] else 'BELOW'}")
    return 0 if res["floors_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
