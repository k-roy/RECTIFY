#!/usr/bin/env python3
"""
Generate dark-mode versions of all v3 figures.

Replaces each script's PAL dict with a dark-theme palette and exec()s it
with a `_dark` output suffix.  No changes to the individual generation scripts.
"""
import os
import re

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ── Unified dark palette (superset of all keys used across all 12 scripts) ──
DARK_PAL = dict(
    # Structural / text
    title   = "#f1f5f9",   # slate-100
    heading = "#cbd5e1",   # slate-300
    label   = "#94a3b8",   # slate-400
    muted   = "#64748b",   # slate-500
    border  = "#334155",   # slate-700
    divider = "#1e293b",   # slate-800
    bg      = "#0f172a",   # slate-950
    # Blues
    blue    = "#60a5fa",   # blue-400
    blue_l  = "#172554",   # blue-950
    # Teals
    teal    = "#2dd4bf",   # teal-400
    teal_l  = "#042f2e",   # teal-950
    # Oranges
    orange  = "#fb923c",   # orange-400
    orange_l= "#431407",   # orange-950
    # Greens
    green   = "#4ade80",   # green-400
    green_l = "#052e16",   # green-950
    # Reds
    red     = "#f87171",   # red-400
    red_l   = "#450a0a",   # red-950
    # Indigos
    indigo  = "#818cf8",   # indigo-400
    indigo_l= "#1e1b4b",   # indigo-950
    # Slates
    slate   = "#cbd5e1",   # slate-300 (label text on dark fills)
    slate_l = "#1e293b",   # slate-800 (subdued box fill)
    # Purples
    purple  = "#c084fc",   # purple-400
    purple_l= "#3b0764",   # purple-950
    # Ambers
    amber   = "#fbbf24",   # amber-400
    amber_l = "#451a03",   # amber-950
    # Base nucleotide colors (A/T/G/C boxes)
    A_f = "#064e3b",  A_t = "#6ee7b7",   # green on dark
    T_f = "#431407",  T_t = "#fdba74",   # orange on dark
    G_f = "#172554",  G_t = "#93c5fd",   # blue on dark
    C_f = "#3b0764",  C_t = "#d8b4fe",   # purple on dark
    # Gap / deletion fill (softclip_rescue, walkback)
    gap_f   = "#1e293b",
    gap_t   = "#64748b",
)


def dark_pal_src():
    lines = ["PAL = dict("]
    for k, v in DARK_PAL.items():
        lines.append(f"    {k:<8}= \"{v}\",")
    lines.append(")")
    return "\n".join(lines)


PAL_RE     = re.compile(r"PAL\s*=\s*dict\(.*?\)", re.DOTALL)
RENDER_RE  = re.compile(r'render\(build\(\),\s*"([^"]+)"\)')

SCRIPTS = [
    ("generate_5prime_junction_v3.py",      "5prime_junction_rescue"),
    ("generate_adaptive_clustering_v3.py",   "adaptive_clustering"),
    ("generate_cdna_isoform_v3.py",          "cdna_isoform_clustering"),
    ("generate_cdna_umi_v3.py",              "cdna_umi_architecture"),
    ("generate_hp_scoring_v3.py",            "hp_scoring"),
    ("generate_multi_aligner_v3.py",         "multi_aligner_consensus"),
    ("generate_netseq_deconvolution_v3.py",  "netseq_deconvolution"),
    ("generate_pipeline_overview_v3.py",     "pipeline_overview"),
    ("generate_polya_pretrim_v3.py",         "polya_pretrim"),
    ("generate_softclip_rescue_v3.py",       "softclip_rescue"),
    ("generate_splice_classification_v3.py", "splice_classification"),
    ("generate_walkback_v3.py",              "walkback_readvsref"),
]


for script_file, fig_name in SCRIPTS:
    path = os.path.join(OUTDIR, script_file)
    src  = open(path).read()

    # 1. Swap in the dark palette
    src = PAL_RE.sub(dark_pal_src(), src, count=1)

    # 2. Redirect output to _dark names
    src = RENDER_RE.sub(f'render(build(), "{fig_name}_dark")', src)

    # 3. Force the __main__ block to always run
    src = src.replace('if __name__ == "__main__":', "if True:")

    exec(compile(src, path, "exec"), {"__name__": "__main__", "__file__": path})
    print(f"  ✓ {fig_name}_dark")
