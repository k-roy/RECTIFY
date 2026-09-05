"""
Design tokens for RECTIFY visualization -- the RNA figure house standard.

The token file is the SINGLE SOURCE OF TRUTH for colour, type and geometry in
``rectify.visualize``. Two colour layers live in it, and the rule that keeps them
apart is the point of the standard (Chanfreau planning/871, ruled 2026-09-05):

* **Layer A -- argument roles** (``focal``, ``reference``, ``match``, ``stratum_a``,
  ``stratum_b``, ``control``): applied to DATA panels and to SIGNAL glyphs coloured
  by sample (coverage fills, junction arcs drawn as counts, pileups).
* **Layer B -- molecular identities** (``splice``, ``polya``, ``mod``): applied to
  RNA glyphs -- gene models, site marks, annotation arcs. Neutral-first: exons,
  introns, reads and single-sample coverage are drawn in the grey ramp.

A layer-A colour never touches a molecular glyph; a layer-B colour never encodes a
sample, a condition or a claim. The drawing API enforces it: ``tracks.gene_model``
and ``tracks.mark`` take no colour argument, ``tracks.coverage`` takes a ROLE.

Resolution order for the token file:

1. ``$RNA_FIGURE_TOKENS`` (an explicit path -- how the ``rna-figure`` skill points
   rectify at its own master copy on a workstation);
2. the bundled ``rna_tokens.json`` next to this module, which is a byte-identical
   copy of ``~/.claude/skills/rna-figure/tokens.json`` (checked by that skill's
   ``tests/sync_tokens.sh --check``).

Author: Kevin R. Roy
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List, Optional

BUNDLED = Path(__file__).with_name("rna_tokens.json")
ENV_VAR = "RNA_FIGURE_TOKENS"

#: the API level of this token-driven layer; the rna-figure skill shim checks it
API_VERSION = 1

_TOK: Optional[dict] = None
_PATH: Optional[Path] = None


def load(path: Optional[os.PathLike] = None, *, force: bool = False) -> dict:
    """Load (and cache) the token file. ``path`` overrides ``$RNA_FIGURE_TOKENS`` and the bundle."""
    global _TOK, _PATH
    if path is None:
        # no explicit path: keep whatever is loaded (a `use()` override must survive a
        # later bare `load()`), else resolve env -> bundle
        if _TOK is not None and not force:
            return _TOK
        env = os.environ.get(ENV_VAR)
        path = Path(env).expanduser() if env else BUNDLED
    path = Path(path)
    if _TOK is not None and not force and path == _PATH:
        return _TOK
    with open(path) as fh:
        _TOK = json.load(fh)
    _PATH = path
    return _TOK


def use(path: os.PathLike) -> dict:
    """Point this process at another token file (the skill's master copy, a project overlay)."""
    return load(path, force=True)


def source() -> Path:
    """Which file the current tokens came from."""
    load()
    return _PATH  # type: ignore[return-value]


# ------------------------------------------------------------------ accessors
def palette() -> Dict[str, str]:
    """Every colour token -> hex (semantic + structural)."""
    t = load()
    return {**t["palette"]["semantic"], **t["palette"]["structural"]}


def color(name: str) -> str:
    """Hex for a token name. Raises ``KeyError`` naming the valid tokens -- a typo must not
    silently fall through to a default."""
    pal = palette()
    try:
        return pal[name]
    except KeyError:
        raise KeyError(f"unknown token {name!r}; valid tokens: {sorted(pal)}") from None


def layer_a() -> List[str]:
    return list(load()["layers"]["A"]["tokens"])


def layer_b() -> List[str]:
    return list(load()["layers"]["B"]["tokens"])


def role(name: str) -> str:
    """Hex for a LAYER-A argument role. Refuses a layer-B token: a molecular identity must not be
    handed to a data series (the boundary rule)."""
    if name in layer_b():
        raise ValueError(f"{name!r} is a layer-B molecular identity, not an argument role; "
                         f"data series take one of {layer_a()} or 'neutral'")
    if name not in layer_a() and name != "neutral":
        raise ValueError(f"unknown role {name!r}; roles are {layer_a()} or 'neutral'")
    return color(name)


def role_cycle() -> List[str]:
    """The default order in which multi-series data panels take roles."""
    t = load()
    order = t["layers"]["A"].get("cycle") or ["reference", "focal", "stratum_a", "stratum_b", "control", "match"]
    return [color(n) for n in order]


def typography(scale: str = "manuscript") -> Dict[str, float]:
    t = load()["typography"]
    key = {"manuscript": "manuscript_pt", "slide": "slide_pt"}[scale]
    return dict(t[key])


def stroke() -> Dict[str, float]:
    return dict(load()["geometry"]["stroke_pt"])


def track_geometry() -> dict:
    return dict(load()["geometry"]["track"])


def apply_rc(matplotlib_module, *, scale: str = "manuscript") -> dict:
    """Apply the house family, the floored type sizes and the stroke tiers to rcParams.

    Deliberately does NOT set figure size, dpi or ``savefig.bbox`` -- those belong to the
    calling figure. Mirrors ``magestic-figure.glyphs.house_rc`` so the two house styles share one
    type scale on a page.
    """
    t = load()
    T = typography(scale)
    S = stroke()
    fam = t["typography"]["family"]
    rc = {
        "font.family": "sans-serif",
        "font.sans-serif": fam,
        "font.size": T["in_figure"],
        "axes.labelsize": T["axis_label"],
        "axes.titlesize": T["axis_label"],
        "xtick.labelsize": T["tick_label"],
        "ytick.labelsize": T["tick_label"],
        "legend.fontsize": T["annotation"],
        "axes.linewidth": S["hairline"],
        "xtick.major.width": S["hairline"],
        "ytick.major.width": S["hairline"],
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": color("hairline"),
        "xtick.color": color("hairline"),
        "ytick.color": color("hairline"),
        "axes.labelcolor": color("ink"),
        "text.color": color("ink"),
        "legend.frameon": False,
        "lines.linewidth": S["secondary"],
        "lines.solid_capstyle": "round",
        "patch.linewidth": S["hairline"],
        # mathtext otherwise renders in DejaVu Sans and ships as a second embedded font
        "mathtext.fontset": "custom",
        "mathtext.rm": fam[0],
        "mathtext.it": fam[0] + ":italic",
        "mathtext.bf": fam[0] + ":bold",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
    matplotlib_module.rcParams.update(rc)
    return rc


__all__ = [
    "BUNDLED", "ENV_VAR", "API_VERSION", "load", "use", "source", "palette", "color",
    "layer_a", "layer_b", "role", "role_cycle", "typography", "stroke", "track_geometry", "apply_rc",
]
