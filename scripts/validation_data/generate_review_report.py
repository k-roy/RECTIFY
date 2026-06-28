#!/usr/bin/env python3
"""Generate an HTML report for visual review of all validation reads.

For each primary read in the bundled validation BAM, drives IGV via HTTP
to produce a screenshot, then assembles a single self-contained HTML file
with per-read metadata, derived interpretation of what RECTIFY did for
that read, and the embedded screenshot.

The HTML report is mobile-friendly: single-column layout, large fonts,
images sized to phone screen width. It auto-detects light/dark mode via
``prefers-color-scheme`` and also offers a manual toggle button.

Output lands in the user's Google Drive folder so it auto-syncs to phone:
    .../Chanfreau Lab/validation_read_review/drs_validation_review.html

Usage::

    python scripts/validation_data/generate_review_report.py [--arm drs] [--limit N] [--html-only]

``--html-only`` skips the screenshot pass and rebuilds the HTML from
existing PNGs in ``/tmp/igv_snapshots/<xv>.png``.
"""
from __future__ import annotations

import argparse
import base64
import csv
import io
import re
import shutil
import socket
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path

import pysam
from PIL import Image

# Add script dir to path so we can import render_read_alignment
sys.path.insert(0, str(Path(__file__).resolve().parent))


REPO = Path("/Users/kevinroy/work/rectify")
VAL_DIR = REPO / "rectify" / "data" / "validation"
DRIVE_OUT = Path(
    "/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinrjroy@gmail.com/"
    "My Drive/Work/Chanfreau Lab/validation_read_review"
)
SNAPSHOT_DIR = Path("/tmp/igv_snapshots")
# Render DPI for the per-read PNGs (overridable via --dpi). Higher = more pixels
# / sharper when opened at native size; layout proportions are unchanged.
RENDER_DPI = 150

GENOME = REPO / "rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
GFF    = REPO / "rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
CORR_BAM    = VAL_DIR / "rectified/rectified_pA_tail_trimmed.bam"
PARESTORE_BAM = VAL_DIR / "rectified/rectified_pA_tail_soft_clipped.bam"
BG_PLUS     = VAL_DIR / "rectified/corrected_3ends.plus.bedgraph"
BG_MINUS    = VAL_DIR / "rectified/corrected_3ends.minus.bedgraph"
SUMMARY_TSV = VAL_DIR / "rectified/per_aligner_summary.tsv"
ALIGNER_DIR = VAL_DIR / "aligners"

CAT_TO_ALIGNER_DRS = {
    "cat1_indel":             "minimap2",
    "cat2_softclip":          "minimap2",
    "cat3_junction":          "minimap2",
    "cat4_false_junc":        "minimap2",
    "cat5_chimeric":          "minimap2",
    "cat6_chimeric":          "mapPacBio",
    "cat7_alt_splice":        "mapPacBio",
    "cat8_netseq_refine":     "minimap2",
    "cat9_junction_refine":   "mapPacBio",
}

MODULE_INTERP = {
    "atract_ambiguity": (
        "<strong>A-tract ambiguity window detected</strong> — the genome at the 3' end "
        "contains a poly-A homopolymer; the walkback's leftmost-possible-CPA principle applies."
    ),
    "indel_correction": (
        "<strong>Module 2C (indel correction)</strong> fired — the aligner had introduced "
        "indels in the terminal alignment over the poly-A region; the walkback recognized the "
        "homopolymer context and resolved the read=ref anchor."
    ),
    "polya_walkback": (
        "<strong>Module 2E (poly-A walkback)</strong> fired — walked back through poly-A "
        "stop-base matches until finding the first non-stop read=ref agreement = "
        "leftmost-possible-CPA."
    ),
    "polya_walkback_readgenome": (
        "<strong>Module 2E (read-vs-reference walkback)</strong> fired through the unified "
        "<code>walkback_3prime</code> core."
    ),
    "softclip_rescue": (
        "<strong>Module 2G (soft-clip rescue at homopolymer)</strong> fired — a 3' soft-clip "
        "adjacent to a genomic homopolymer was extended outward, re-anchoring the 3' end at "
        "the leftmost-possible-CPA."
    ),
    "five_prime_rescued": (
        "<strong>5' junction rescue</strong> fired — a leading soft-clip in the mapped track "
        "was extended through an annotated splice junction into the upstream exon."
    ),
    "homopolymer_rescue": (
        "<strong>Module 2D (variant-aware homopolymer rescue)</strong> fired — a mismatch "
        "inside a homopolymer was rescued after the two-pass variant scan confirmed it's not a SNP."
    ),
    "none": (
        "No correction modules fired — the read's terminal alignment was already at a clean "
        "read=ref non-A match."
    ),
}

CAT_HINTS = {
    "cat1_indel": (
        "Indel correction: terminal poly-A indel artifact in the aligner output. The walkback "
        "should anchor at the leftmost-possible-CPA = first non-A read=ref match. Verify the "
        "genome at the corrected_3prime position is a non-A and there's an A-tract starting "
        "one base to its right."
    ),
    "cat2_softclip": (
        "Soft-clip rescue: 3' soft-clip at a homopolymer boundary. RECTIFY extends the alignment "
        "outward through the soft-clip, anchoring at the leftmost-possible-CPA. Check that the "
        "soft-clip bases in the mapped track are M bases (or restored as 3' soft-clip) in the "
        "corrected track."
    ),
    "cat3_junction": (
        "5' junction rescue: leading soft-clip in the mapped track gets extended back through an "
        "annotated splice junction into the upstream exon. The corrected track should show M ops "
        "where the mapped track had soft-clip."
    ),
    "cat4_false_junc": (
        "False junction filter: aligner introduces a spurious N-op (false intron) near the 3' end. "
        "Walkback should ignore the false junction and land at the leftmost-possible-CPA in the "
        "terminal exon."
    ),
    "cat5_chimeric": (
        "Chimeric consensus: different aligners disagree on the junction set; consensus stitches "
        "the best segments. Look for SA-tagged segments in the mapped track and a clean stitched "
        "read in the corrected track."
    ),
    "cat6_chimeric": (
        "Simple chimeric: shorter chimeric structure, usually one stitch. Same logic as Cat5 but "
        "simpler topology."
    ),
    "cat7_alt_splice": (
        "Alternative splice: junction boundaries that are off by 1–5 bp from the canonical "
        "annotated site. Junction refinement (Module 2H) should adjust to the annotated GT-AG."
    ),
    "cat8_netseq_refine": (
        "NET-seq refinement: walkback alone can't resolve the CPA inside a long A-tract; NET-seq "
        "deconvolution narrows it to a specific position. The corrected position should be inside "
        "the ambiguity_min..ambiguity_max window."
    ),
    "cat9_junction_refine": (
        "Module 2H junction refinement: aligner's N-op boundaries off by 1–5 bp from annotated "
        "splice site. After refinement, the corrected N-op should match an annotated junction "
        "exactly."
    ),
}


def interpret_correction(corr_applied: str, shift_bp: int) -> str:
    modules = [m.strip() for m in corr_applied.split(",") if m.strip()] or ["none"]
    out = []
    if shift_bp == 0:
        out.append(
            '<p><em>Zero-shift correction.</em> The terminal alignment position was already at '
            'the leftmost-possible-CPA; the listed module(s) verified this but did not move '
            'the 3\' end.</p>'
        )
    else:
        direction = "upstream" if shift_bp < 0 else "downstream"
        out.append(
            f'<p><em>{abs(shift_bp)} bp {direction} shift.</em> The walkback moved the 3\' end '
            f'from the aligner\'s original position to the leftmost-possible-CPA by absorbing '
            f'{abs(shift_bp)} bp of poly-A artifact.</p>'
        )
    out.append("<ul>")
    for m in modules:
        out.append(f"<li>{MODULE_INTERP.get(m, f'<code>{m}</code> module fired.')}</li>")
    out.append("</ul>")
    return "\n".join(out)


def url_encode(s: str) -> str:
    return urllib.parse.quote(s)


def igv_request(path: str, timeout: float = 5.0) -> int:
    """Send a GET to IGV's command port. Returns HTTP status code, or 0 on connection failure."""
    try:
        with urllib.request.urlopen(f"http://localhost:60151{path}", timeout=timeout) as resp:
            return resp.status
    except urllib.error.HTTPError as e:
        return e.code
    except (urllib.error.URLError, socket.timeout, ConnectionResetError):
        return 0


def is_port_up() -> bool:
    s = socket.socket()
    s.settimeout(0.5)
    try:
        return s.connect_ex(("127.0.0.1", 60151)) == 0
    finally:
        s.close()


def render_plot(xv: str, qname: str, chrom: str, left: int, right: int,
                orig_3p: int, corr_3p: int, mapped_bam: Path, title: str,
                minimap2_bam: Path | None = None,
                winner_label: str = "",
                gene_id: str = "") -> Path:
    """Render a per-base alignment plot using render_read_alignment.

    When ``minimap2_bam`` is given and differs from ``mapped_bam``, the
    renderer adds a separate minimap2 row above the winner row so the
    reviewer can compare both upstream alignments side-by-side.

    When ``SUMMARY_TSV`` exists, the End-correction + Winner-selection
    panel is automatically added at the top of the figure.
    """
    from render_read_alignment import render as render_alignment
    out = SNAPSHOT_DIR / f"{xv}.png"
    render_alignment(
        qname=qname, chrom=chrom,
        win_start=left, win_end=right,
        orig_3p=orig_3p, corr_3p=corr_3p,
        genome_fa=GENOME,
        mapped_bam=mapped_bam,
        minimap2_bam=minimap2_bam,
        winner_label=winner_label,
        corrected_bam=CORR_BAM,
        parestore_bam=PARESTORE_BAM,
        bg_plus=BG_PLUS, bg_minus=BG_MINUS,
        out_path=out, title=title,
        summary_tsv=SUMMARY_TSV if SUMMARY_TSV.exists() else None,
        aligner_bam_dir=ALIGNER_DIR if ALIGNER_DIR.exists() else None,
        gene_id=gene_id,
        dpi=RENDER_DPI,
    )
    return out


# Keep the IGV-based path around as a fallback (unused by default)
def take_screenshot(xv: str, chrom: str, left: int, right: int, mapped_bam: Path) -> Path:
    """LEGACY: drive IGV via HTTP + screencap. Slow (~1 min/read). Kept for fallback."""
    files = [GFF, mapped_bam, CORR_BAM, PARESTORE_BAM, BG_PLUS, BG_MINUS]
    files_csv = ",".join(url_encode(str(p)) for p in files)
    igv_request(f"/load?file={files_csv}&genome={url_encode(str(GENOME))}")
    igv_request(f"/goto?locus={url_encode(f'{chrom}:{left}-{right}')}")
    subprocess.run(
        ["osascript", "-e",
         'tell application "System Events" to set frontmost of (first process whose name is "java") to true'],
        check=False, capture_output=True,
    )
    time.sleep(1.2)
    out = SNAPSHOT_DIR / f"{xv}.png"
    subprocess.run(["screencapture", "-x", str(out)], check=True)
    return out


def png_to_base64(png_path: Path, max_width: int = 1100) -> str:
    img = Image.open(png_path)
    if img.size[0] > max_width:
        ratio = max_width / img.size[0]
        img = img.resize((max_width, int(img.size[1] * ratio)), Image.Resampling.LANCZOS)
    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return base64.b64encode(buf.getvalue()).decode("ascii")


HTML_HEAD = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>{title}</title>
<style>
:root {{
  --bg: #fafafa;
  --bg-card: #ffffff;
  --bg-toc: #f0f4f9;
  --text: #1a1a1a;
  --text-muted: #555;
  --border: #d0d7de;
  --border-strong: #888;
  --link: #0a4ea0;
  --hint-bg: #fff8dc;
  --hint-border: #d9b738;
  --hint-text: #5a4500;
  --interp-bg: #e8f4ff;
  --interp-border: #2a6fb8;
  --interp-text: #0e2e58;
  --note-bg: #fffce8;
  --note-text: #555;
  --strand-plus: #2e7d32;
  --strand-minus: #c62828;
  --shift-zero: #888;
  --shift-nonzero: #1976d2;
  --code-bg: #f3f3f3;
  --shadow: 0 1px 3px rgba(0,0,0,0.06);
}}
[data-theme="dark"] {{
  --bg: #0d1117;
  --bg-card: #161b22;
  --bg-toc: #1a2333;
  --text: #e6edf3;
  --text-muted: #8b949e;
  --border: #30363d;
  --border-strong: #586069;
  --link: #58a6ff;
  --hint-bg: #2a2418;
  --hint-border: #d9b738;
  --hint-text: #e3c878;
  --interp-bg: #0f2440;
  --interp-border: #58a6ff;
  --interp-text: #a5cfff;
  --note-bg: #1f1d10;
  --note-text: #b8b8a0;
  --strand-plus: #56d364;
  --strand-minus: #ff7b72;
  --shift-zero: #6e7681;
  --shift-nonzero: #58a6ff;
  --code-bg: #2d333b;
  --shadow: 0 1px 3px rgba(0,0,0,0.4);
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg: #0d1117;
    --bg-card: #161b22;
    --bg-toc: #1a2333;
    --text: #e6edf3;
    --text-muted: #8b949e;
    --border: #30363d;
    --border-strong: #586069;
    --link: #58a6ff;
    --hint-bg: #2a2418;
    --hint-border: #d9b738;
    --hint-text: #e3c878;
    --interp-bg: #0f2440;
    --interp-border: #58a6ff;
    --interp-text: #a5cfff;
    --note-bg: #1f1d10;
    --note-text: #b8b8a0;
    --strand-plus: #56d364;
    --strand-minus: #ff7b72;
    --shift-zero: #6e7681;
    --shift-nonzero: #58a6ff;
    --code-bg: #2d333b;
    --shadow: 0 1px 3px rgba(0,0,0,0.4);
  }}
}}
body {{
  font-family: -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
  background: var(--bg);
  color: var(--text);
  max-width: 1200px;
  margin: 0 auto;
  padding: 12px 14px 80px;
  line-height: 1.5;
  font-size: 16px;
  transition: background-color 0.2s ease, color 0.2s ease;
}}
.header-bar {{
  display: flex; align-items: center; justify-content: space-between;
  margin-bottom: 16px; padding-bottom: 12px;
  border-bottom: 1px solid var(--border);
}}
.header-bar h1 {{ margin: 0; font-size: 22px; }}
.theme-toggle {{
  background: var(--bg-card); border: 1px solid var(--border);
  border-radius: 20px; padding: 6px 14px; font-size: 14px;
  cursor: pointer; color: var(--text);
}}
.theme-toggle:hover {{ border-color: var(--link); }}
.summary {{ color: var(--text-muted); font-size: 14px; margin: 10px 0 18px; }}
.principle {{
  background: var(--interp-bg); color: var(--interp-text);
  border-left: 4px solid var(--interp-border);
  padding: 10px 14px; margin: 16px 0 24px;
  border-radius: 6px; font-size: 14px;
}}
h2 {{
  font-size: 19px; margin: 36px 0 8px; padding: 14px 16px;
  background: var(--bg-card); border: 1px solid var(--border);
  border-radius: 8px; box-shadow: var(--shadow);
}}
h2 small {{ font-size: 14px; color: var(--text-muted); font-weight: normal; }}
h3 {{ font-size: 15px; margin: 18px 0 6px; color: var(--text-muted);
      text-transform: uppercase; letter-spacing: 0.5px; }}
table {{ border-collapse: collapse; margin: 6px 0; font-size: 14px;
        background: var(--bg-card); box-shadow: var(--shadow);
        border-radius: 6px; overflow: hidden; }}
td, th {{ border: 1px solid var(--border); padding: 6px 10px; vertical-align: top; }}
th {{ background: var(--bg-toc); text-align: left; font-weight: 500;
      color: var(--text-muted); }}
code {{ font-family: "SF Mono", "Monaco", monospace; font-size: 13px;
       background: var(--code-bg); color: var(--text);
       padding: 1px 6px; border-radius: 3px; }}
.hint {{
  background: var(--hint-bg); color: var(--hint-text);
  border-left: 4px solid var(--hint-border);
  padding: 10px 14px; margin: 8px 0 14px;
  border-radius: 6px; font-size: 14px;
}}
.interp {{
  background: var(--interp-bg); color: var(--interp-text);
  border-left: 4px solid var(--interp-border);
  padding: 10px 14px; margin: 8px 0 14px;
  border-radius: 6px; font-size: 14px;
}}
.interp ul {{ margin: 6px 0; padding-left: 22px; }}
.interp li {{ margin: 5px 0; }}
.interp p {{ margin: 6px 0; }}
img.screenshot {{
  width: 100%; max-width: 1100px;
  border: 1px solid var(--border-strong);
  border-radius: 8px; display: block; margin: 12px 0;
  box-shadow: 0 2px 12px rgba(0,0,0,0.18);
}}
.toc {{
  background: var(--bg-toc);
  border: 1px solid var(--border);
  padding: 14px 18px; border-radius: 8px;
  margin: 14px 0 28px;
}}
.toc h3 {{ margin: 0 0 10px; color: var(--text); text-transform: none; letter-spacing: 0; }}
.toc-group {{ margin: 8px 0; }}
.toc-group .cat {{ display: inline-block; font-weight: 600; min-width: 220px;
                  color: var(--text-muted); font-size: 13px; }}
.toc a {{
  display: inline-block; margin: 2px 4px; padding: 4px 10px;
  background: var(--bg-card); border: 1px solid var(--border);
  border-radius: 14px; text-decoration: none;
  font-size: 13px; font-family: "SF Mono", monospace;
}}
.toc a.plus  {{ color: var(--strand-plus); }}
.toc a.minus {{ color: var(--strand-minus); }}
.toc a:hover {{ border-color: var(--link); }}
.meta-strand-plus  {{ color: var(--strand-plus); font-weight: 600; }}
.meta-strand-minus {{ color: var(--strand-minus); font-weight: 600; }}
.shift-zero    {{ color: var(--shift-zero); }}
.shift-nonzero {{ color: var(--shift-nonzero); font-weight: 600; }}
hr {{ border: none; border-top: 1px solid var(--border); margin: 40px 0; }}
.review-status {{
  margin: 8px 0; padding: 10px 14px;
  border: 2px dashed var(--border-strong);
  border-radius: 6px;
  background: var(--note-bg);
  color: var(--note-text);
  font-style: italic; font-size: 14px;
}}
.back-to-top {{
  position: fixed; bottom: 18px; right: 18px;
  background: var(--bg-card); border: 1px solid var(--border);
  border-radius: 50%; width: 44px; height: 44px;
  display: flex; align-items: center; justify-content: center;
  text-decoration: none; color: var(--text);
  box-shadow: var(--shadow); font-size: 18px;
}}
</style></head><body>
<div class="header-bar">
  <h1>{title}</h1>
  <button class="theme-toggle" id="theme-toggle">🌓 theme</button>
</div>
<div class="summary">Generated {now} · {n_reads} reads across {n_cats} categories</div>
<div class="principle">
<p><strong>Leftmost-possible-CPA principle</strong> — in mature mRNA, the boundary between
genomic A's and post-transcriptionally added poly(A) is irrecoverably lost. RECTIFY's walkback
returns the principled answer to that ambiguity: anchor at the first non-stop read=ref agreement
walking inward from the alignment 3' edge.</p>
<p>When NET-seq data is available (a separate nascent-RNA dataset that captures CPA cleavage
intermediates as a mixture of ~50% non-adenylated and ~50% oligo-adenylated (1–12 nt) species),
<code>rectify correct --netseq-dir</code> deconvolves the NET-seq peaks within A-tracts and
apportions the corrected 3' end within the CPA ambiguity window to NET-seq–informed CPA sites in
the observed proportions (see RECTIFY README). The walkback alone gives you the principled left
edge of that window; NET-seq refinement pins the position inside it.</p>
<p>For each read below, verify visually that the corrected position sits at the
leftmost-possible-CPA anchor (or, for Cat8 reads, at the NET-seq peak inside the ambiguity window).</p>
</div>
"""

HTML_TOOLBAR_JS = """<a class="back-to-top" href="#top" title="Back to top">↑</a>
<script>
(function() {
  const btn = document.getElementById('theme-toggle');
  const root = document.documentElement;
  const stored = localStorage.getItem('theme');
  if (stored) root.setAttribute('data-theme', stored);
  btn.addEventListener('click', () => {
    const isDark = root.getAttribute('data-theme') === 'dark' ||
                   (!root.getAttribute('data-theme') &&
                    window.matchMedia('(prefers-color-scheme: dark)').matches);
    const next = isDark ? 'light' : 'dark';
    root.setAttribute('data-theme', next);
    localStorage.setItem('theme', next);
  });
})();
</script>
"""


def build_toc(reads_meta: list[dict]) -> str:
    by_cat: dict[str, list[dict]] = {}
    for m in reads_meta:
        by_cat.setdefault(m["category"], []).append(m)
    lines = ['<div class="toc"><h3>Jump to category</h3>']
    for cat in sorted(by_cat):
        lines.append('<div class="toc-group">')
        lines.append(f'<span class="cat">{cat}</span>')
        for m in by_cat[cat]:
            cls = "plus" if m["strand"] == "+" else "minus"
            lines.append(f'<a href="#{m["xv"]}" class="{cls}">{m["xv"]}</a>')
        lines.append("</div>")
    lines.append("</div>")
    return "\n".join(lines)


def _html_to_md(s: str) -> str:
    """Convert the limited HTML subset used in MODULE_INTERP / CAT_HINTS /
    interpret_correction output into plain Markdown."""
    s = re.sub(r"</?(strong|b)>", "**", s)
    s = re.sub(r"</?(em|i)>", "*", s)
    s = re.sub(r"<code>([^<]*)</code>", r"`\1`", s)
    s = re.sub(r"<br\s*/?>", "  \n", s)
    s = re.sub(r"<li>\s*", "- ", s)
    s = re.sub(r"\s*</li>", "\n", s)
    s = re.sub(r"</?(ul|ol)>\s*", "", s)
    s = re.sub(r"<p>\s*", "", s)
    s = re.sub(r"\s*</p>", "\n\n", s)
    return s.strip()


def build_read_section_md(m: dict, png_rel: str) -> str:
    shift = int(m["corrected_3prime"]) - int(m["original_3prime"])
    interp = _html_to_md(interpret_correction(m["correction_applied"], shift))
    hint = _html_to_md(CAT_HINTS.get(m["category"], ""))
    rows = [
        ("read_id", f"`{m['read_id']}`"),
        ("category", m["category"]),
        ("strand", m["strand"]),
        ("aligner shown", m["aligner"]),
        ("chrom", m["chrom"]),
        ("locus shown", m["locus"]),
        ("alignment span", f"[{m['alignment_start']}, {m['alignment_end']})"),
        ("5′ position", m["five_prime_position"] or "—"),
        ("original 3′", m["original_3prime"]),
        ("corrected 3′", m["corrected_3prime"]),
        ("shift (corr − orig)", f"**{shift:+d} bp**" if shift else "0 bp (zero-shift)"),
        ("correction_applied", f"`{m['correction_applied']}`"),
    ]
    table = "| field | value |\n| --- | --- |\n" + "\n".join(
        f"| {k} | {v} |" for k, v in rows
    )
    return (
        f"## {m['xv']} <a id=\"{m['xv']}\"></a> — {m['category']}\n\n"
        f"{table}\n\n"
        f"**What RECTIFY did for this read**\n\n{interp}\n\n"
        f"**What to verify visually**\n\n> {hint}\n\n"
        f"![{m['xv']}]({png_rel})\n\n"
        f"_Reviewer notes:_ [ ] OK   [ ] flag — note here\n"
    )


MD_PRINCIPLE = (
    "> **Leftmost-possible-CPA principle** — in mature mRNA, the boundary between "
    "genomic A's and post-transcriptionally added poly(A) is irrecoverably lost. "
    "RECTIFY's walkback anchors at the first non-stop read=ref agreement walking "
    "inward from the alignment 3′ edge. When NET-seq is available, "
    "`rectify correct --netseq-dir` pins the corrected 3′ within that ambiguity window.\n"
)


def _copy_png(src: Path, dest_dir: Path) -> Path:
    """Copy a PNG into dest_dir/<src.name> (idempotent by mtime+size)."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / src.name
    if dest.exists() and dest.stat().st_size == src.stat().st_size \
            and dest.stat().st_mtime >= src.stat().st_mtime:
        return dest
    shutil.copy2(src, dest)
    return dest


def write_md_report(reads: list[dict], title: str, out_path: Path) -> None:
    """Write a Markdown report. PNGs are copied next to the .md in a
    sibling folder named ``<out_path stem>_pngs/`` and referenced relatively
    so the file works when opened directly in VS Code's markdown preview."""
    png_dir = out_path.parent / f"{out_path.stem}_pngs"
    sections: list[str] = []
    for m in reads:
        copied = _copy_png(Path(m["png_path"]), png_dir)
        rel = f"{png_dir.name}/{copied.name}"
        sections.append(build_read_section_md(m, rel))

    by_cat: dict[str, list[dict]] = {}
    for m in reads:
        by_cat.setdefault(m["category"], []).append(m)
    toc_lines = ["## Jump to read\n"]
    for cat in sorted(by_cat):
        links = " · ".join(f"[{m['xv']}](#{m['xv']})" for m in by_cat[cat])
        toc_lines.append(f"- **{cat}** — {links}")
    toc = "\n".join(toc_lines)

    header = (
        f"# {title}\n\n"
        f"_Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} · "
        f"{len(reads)} reads across {len(by_cat)} categor"
        f"{'y' if len(by_cat) == 1 else 'ies'}_\n\n"
        f"{MD_PRINCIPLE}\n"
        f"{toc}\n\n---\n\n"
    )
    out_path.write_text(header + "\n---\n\n".join(sections), encoding="utf-8")
    size_kb = out_path.stat().st_size / 1024
    print(f"Wrote {out_path} ({size_kb:.1f} KB markdown + PNGs in {png_dir.name}/)")


def build_read_section(m: dict, png_b64: str) -> str:
    shift = int(m["corrected_3prime"]) - int(m["original_3prime"])
    shift_class = "shift-zero" if shift == 0 else "shift-nonzero"
    strand_class = "meta-strand-plus" if m["strand"] == "+" else "meta-strand-minus"
    interp = interpret_correction(m["correction_applied"], shift)
    return f"""<h2 id="{m['xv']}">{m['xv']} <small>({m['category']})</small></h2>
<table>
  <tr><th>read_id</th><td><code>{m['read_id']}</code></td></tr>
  <tr><th>category</th><td>{m['category']}</td></tr>
  <tr><th>strand</th><td class="{strand_class}">{m['strand']}</td></tr>
  <tr><th>aligner shown</th><td>{m['aligner']}</td></tr>
  <tr><th>chrom</th><td>{m['chrom']}</td></tr>
  <tr><th>locus shown</th><td>{m['locus']}</td></tr>
  <tr><th>alignment span</th><td>[{m['alignment_start']}, {m['alignment_end']})</td></tr>
  <tr><th>5' position</th><td>{m['five_prime_position']}</td></tr>
  <tr><th>original 3'</th><td>{m['original_3prime']}</td></tr>
  <tr><th>corrected 3'</th><td>{m['corrected_3prime']}</td></tr>
  <tr><th>shift (corr − orig)</th><td class="{shift_class}">{shift:+d} bp</td></tr>
  <tr><th>correction_applied</th><td><code>{m['correction_applied']}</code></td></tr>
</table>
<h3>What RECTIFY did for this read</h3>
<div class="interp">{interp}</div>
<h3>What to verify visually</h3>
<div class="hint">{CAT_HINTS.get(m['category'], '')}</div>
<img class="screenshot" src="data:image/png;base64,{png_b64}" alt="{m['xv']} IGV view">
<div class="review-status">Reviewer notes: [&nbsp;&nbsp;] OK &nbsp; [&nbsp;&nbsp;] flag — note here</div>
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=["drs"], default="drs")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--pad", type=int, default=50)
    parser.add_argument("--dpi", type=int, default=150,
                        help="Render DPI for the per-read PNGs (default 150). "
                             "Higher = sharper at native size; layout unchanged.")
    parser.add_argument("--out", type=Path, default=None,
                        help="Output HTML path (single-file mode only). Ignored with --per-category.")
    parser.add_argument("--html-only", action="store_true",
                        help="Skip screenshot pass; rebuild HTML from existing PNGs")
    parser.add_argument("--category", action="append", default=None,
                        help="Filter to this XG category (repeatable). Examples: "
                             "cat1_indel, cat2_softclip, cat3_junction, cat4_false_junc, "
                             "cat5_chimeric, cat6_chimeric, cat7_alt_splice, "
                             "cat8_netseq_refine, cat9_junction_refine.")
    parser.add_argument("--per-category", action="store_true",
                        help="Write one HTML file per category instead of a single combined report. "
                             "Useful for iterative phone review of one 4-read category at a time. "
                             "Files are named <category>_review.html in the same output dir.")
    parser.add_argument("--format", choices=["html", "md", "both"], default="html",
                        help="Output format. md is VS-Code-friendly (markdown preview + "
                             "sibling PNG folder). both writes html and md side-by-side.")
    parser.add_argument("--fresh", action="store_true",
                        help="Re-correct the validation reads with the CURRENT code "
                             "(regen_pa_rest_bundle.py) before rendering, so the figure "
                             "reflects HEAD. NOTE: this currently regenerates the committed "
                             "rectified/ fixtures IN PLACE (mutates tracked files + takes a "
                             "few minutes). A scratch-dir version is a TODO. Without --fresh, "
                             "the renderer uses the committed fixtures and prints a loud "
                             "staleness banner so you are never silently fooled by old data.")
    args = parser.parse_args()

    global RENDER_DPI
    RENDER_DPI = args.dpi

    if args.arm != "drs":
        print(f"--arm {args.arm} not yet supported", file=sys.stderr)
        return 1

    if args.fresh and not args.html_only:
        print("[--fresh] re-correcting validation reads with current code "
              "(regen_pa_rest_bundle.py) — this mutates the committed rectified/ "
              "fixtures in place...", file=sys.stderr, flush=True)
        _rc = subprocess.run(
            [sys.executable, str(Path(__file__).parent / "regen_pa_rest_bundle.py")]
        )
        if _rc.returncode != 0:
            print("[--fresh] re-correction FAILED — refusing to render a possibly-stale "
                  "figure. Fix the regen, or omit --fresh to render committed fixtures.",
                  file=sys.stderr)
            return 1
    elif not args.html_only:
        # Loud staleness banner: the committed fixtures can drift from the code.
        import subprocess as _sp
        _date = _sp.run(["git", "log", "-1", "--format=%ci",
                         "--", str(VAL_DIR / "rectified" / "per_aligner_summary.tsv")],
                        capture_output=True, text=True).stdout.strip() or "unknown"
        print("=" * 78, file=sys.stderr)
        print(f"[STALENESS WARNING] Rendering from COMMITTED fixtures last regenerated "
              f"{_date}.\n  These can drift from current code (winner/ED/3'-end may be "
              f"obsolete).\n  Pass --fresh to re-correct with HEAD first.", file=sys.stderr)
        print("=" * 78, file=sys.stderr)

    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    DRIVE_OUT.mkdir(parents=True, exist_ok=True)
    bam_path = VAL_DIR / "validation_reads.bam"
    tsv_path = VAL_DIR / "corrected_reads.tsv"

    # Load XV → (uuid, category, strand) once
    xv_map: dict[str, tuple[str, str, str]] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            try:
                xv = r.get_tag("XV")
                xv_map[xv] = (r.query_name, r.get_tag("XG"), "-" if r.is_reverse else "+")
            except KeyError:
                pass

    # Load corrected TSV once
    tsv_map: dict[str, dict] = {}
    with tsv_path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tsv_map[row["read_id"]] = row

    xvs = sorted(xv_map.keys())
    # Filter by category if --category was given (single-file mode); for
    # --per-category mode we keep all reads and group at output time.
    if args.category:
        wanted_cats = set(args.category)
        xvs = [xv for xv in xvs if xv_map[xv][1] in wanted_cats]
        missing = wanted_cats - {xv_map[xv][1] for xv in xvs}
        if missing:
            print(f"[warn] no reads found for categories: {sorted(missing)}",
                  file=sys.stderr)
        if not xvs:
            print(f"[error] no reads after --category filter; exiting", file=sys.stderr)
            return 1
        print(f"[--category] filtered to {len(xvs)} reads in categories {sorted(wanted_cats)}",
              file=sys.stderr)
    if args.limit:
        xvs = xvs[: args.limit]

    # IGV-port check disabled: this generator now uses the matplotlib-based
    # render_plot path which doesn't need IGV at all. Only the legacy
    # take_screenshot() path requires IGV. Keep is_port_up() callable for
    # future use but don't gate on it.
    _ = is_port_up  # silence unused warning

    reads_meta = []
    t0 = time.time()
    for i, xv in enumerate(xvs, 1):
        uuid, category, strand = xv_map[xv]
        aligner = CAT_TO_ALIGNER_DRS.get(category, "minimap2")
        row = tsv_map.get(uuid)
        if not row:
            print(f"[{i}/{len(xvs)}] WARN: {uuid} not in corrected_reads.tsv", file=sys.stderr)
            continue
        corr = int(row["corrected_3prime"])
        left = max(1, corr - args.pad)
        right = corr + args.pad

        if args.html_only:
            png_path = SNAPSHOT_DIR / f"{xv}.png"
            if not png_path.exists():
                print(f"[{i}/{len(xvs)}] --html-only: no PNG for {xv}", file=sys.stderr)
                continue
            print(f"[{i}/{len(xvs)}] {xv} ({category}, {strand}, {aligner}) — using existing PNG")
        else:
            mapped_bam = VAL_DIR / "aligners" / f"validation_reads.{aligner}.bam"
            # Always pass the minimap2 BAM separately. The renderer compares
            # paths and only adds a second upstream row when winner != minimap2
            # (i.e. cat6_chimeric / cat7_alt_splice / cat9_junction_refine).
            minimap2_bam = VAL_DIR / "aligners" / "validation_reads.minimap2.bam"
            title = f"{xv} ({category}, strand={strand}) — {row['chrom']}:{left}-{right} (orig {row['original_3prime']} -> corr {corr})"
            print(f"[{i}/{len(xvs)}] {xv} ({category}, {strand}, {aligner}) → {row['chrom']}:{left}-{right}", flush=True)
            png_path = render_plot(
                xv, uuid, row["chrom"], left, right,
                int(row["original_3prime"]), corr,
                mapped_bam, title,
                minimap2_bam=minimap2_bam,
                winner_label=aligner,
                gene_id=row.get("gene_id", ""),
            )

        reads_meta.append({
            "xv": xv,
            "category": category,
            "strand": strand,
            "aligner": aligner,
            "read_id": uuid,
            "chrom": row["chrom"],
            "locus": f"{row['chrom']}:{left}-{right}",
            "alignment_start": row["alignment_start"],
            "alignment_end": row["alignment_end"],
            "five_prime_position": row.get("five_prime_position", ""),
            "original_3prime": row["original_3prime"],
            "corrected_3prime": row["corrected_3prime"],
            "correction_applied": row.get("correction_applied", "none"),
            "png_path": str(png_path),
        })
    print(f"Screenshot pass done in {time.time()-t0:.1f}s. Encoding PNGs...", flush=True)

    def _write_html(reads: list[dict], title: str, out_path: Path):
        sections = [build_read_section(m, png_to_base64(Path(m["png_path"])))
                    for m in reads]
        n_cats = len({m["category"] for m in reads})
        head = HTML_HEAD.format(
            title=title,
            now=datetime.now().strftime("%Y-%m-%d %H:%M"),
            n_reads=len(reads),
            n_cats=n_cats,
        )
        body = "<a id='top'></a>\n" + build_toc(reads) + "\n<hr>\n".join(sections)
        html = head + body + HTML_TOOLBAR_JS + "\n</body></html>"
        out_path.write_text(html, encoding="utf-8")
        size_mb = out_path.stat().st_size / 1024 / 1024
        print(f"Wrote {out_path} ({size_mb:.1f} MB)")

    want_html = args.format in ("html", "both")
    want_md = args.format in ("md", "both")

    if args.per_category:
        # Group by category and write one report per category. Categories are
        # rendered alphabetically; within a category, reads stay in the same
        # XV sort order as the single-file mode (cat1_minus_1 → cat1_plus_2).
        from collections import defaultdict
        by_cat: dict[str, list[dict]] = defaultdict(list)
        for m in reads_meta:
            by_cat[m["category"]].append(m)
        out_dir = (args.out.parent if args.out else DRIVE_OUT)
        out_dir.mkdir(parents=True, exist_ok=True)
        for cat in sorted(by_cat):
            cat_reads = by_cat[cat]
            cat_title = f"RECTIFY DRS validation — {cat} ({len(cat_reads)} reads)"
            if want_html:
                _write_html(cat_reads, cat_title, out_dir / f"{cat}_review.html")
            if want_md:
                write_md_report(cat_reads, cat_title, out_dir / f"{cat}_review.md")
    else:
        title_suffix = f" — {','.join(sorted(args.category))}" if args.category else ""
        combined_title = f"RECTIFY DRS validation read review{title_suffix}"
        if want_html:
            out_path = args.out or (DRIVE_OUT / "drs_validation_review.html")
            if args.out and args.out.suffix == ".md":
                out_path = args.out.with_suffix(".html")
            _write_html(reads_meta, combined_title, out_path)
        if want_md:
            md_path = (args.out.with_suffix(".md") if args.out
                       else DRIVE_OUT / "drs_validation_review.md")
            write_md_report(reads_meta, combined_title, md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
