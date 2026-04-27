#!/usr/bin/env python3
"""
generate_igv_html.py — regenerate validation_reads_igv.html from live data.

Reads validation_reads.bam (read spans, CIGAR N-ops, tags) and
corrected_3ends.tsv (original/corrected positions, junction strings)
to compute all IGV loci automatically.

Run:
    cd rectify/data/validation/
    python generate_igv_html.py

Author: Kevin R. Roy
"""

import csv
import pysam
from pathlib import Path

VAL_DIR   = Path(__file__).parent
BAM_PATH  = VAL_DIR / "validation_reads.bam"
TSV_PATH  = VAL_DIR / "corrected_3ends.tsv"
HTML_PATH = VAL_DIR / "validation_reads_igv.html"
IGV_PORT  = 60151

# ---------------------------------------------------------------------------
# Per-read metadata that requires annotation / aligner-BAM lookup
# Update this dict if Cat3-9 reads are ever replaced.
# ---------------------------------------------------------------------------
READ_META = {
    # Cat3: soft-clip column text, raw 5'-end detail
    "cat3_plus_1":  {"clip": "10S",  "fp_detail": "168808 (reference_start); intron 168424–168808 upstream"},
    "cat3_plus_2":  {"clip": "22S",  "fp_detail": "142618 (reference_start)"},
    "cat3_minus_1": {"clip": "— (mapPacBio forced-mismatch; mpb_mismatch path)",
                     "fp_detail": "900766 (reference_end−1); intron 900767–901193 upstream"},
    "cat3_minus_2": {"clip": "25S",  "fp_detail": "366502 (reference_end−1)"},

    # Cat4: description of spurious N-op, whether it was absorbed, corrected 3' display
    "cat4_plus_1":  {"n_desc": "chrXI:20527–22047 (1520 bp)",  "absorbed": "Yes — snapped to 20526",                    "corr_disp": "20526"},
    "cat4_plus_2":  {"n_desc": "chrX:393725–393825 (100 bp)",  "absorbed": "Yes — window clipped to [393721,393724]",   "corr_disp": "393721"},
    "cat4_minus_1": {"n_desc": "chrI:128521–129021 (500 bp)",  "absorbed": "No — N is 424 bp from 3′",                  "corr_disp": "128097"},
    "cat4_minus_2": {"n_desc": "chrIX:76027–76250 (223 bp)",   "absorbed": "No — N remains at CIGAR start (12H 223N …)", "corr_disp": "76027"},

    # Cat5: aligner names, source intron text, complementary intron text + coords
    "cat5_plus_1":  {"aligners": "mapPacBio, minimap2", "src": "423590–423951 (mapPacBio)", "cmp": "424421–425030 (minimap2)", "cmp_s": 424421, "cmp_e": 425030},
    "cat5_plus_2":  {"aligners": "minimap2, mapPacBio", "src": "332875–333386 (minimap2)",  "cmp": "334050–334122 (mapPacBio)", "cmp_s": 334050, "cmp_e": 334122},
    "cat5_minus_1": {"aligners": "minimap2, gapmm2",   "src": "437941–438397 (gapmm2)",    "cmp": "436480–437396 (minimap2)",  "cmp_s": 436480, "cmp_e": 437396},
    "cat5_minus_2": {"aligners": "mapPacBio, gapmm2",  "src": "177906–178213 (mapPacBio)", "cmp": "176709–177362 (gapmm2)",    "cmp_s": 176709, "cmp_e": 177362},

    # Cat6: minimap2 soft-clip length (the clip that mapPacBio resolves)
    "cat6_plus_1":  {"clip": "14S"},
    "cat6_plus_2":  {"clip": "36S"},
    "cat6_minus_1": {"clip": "47S"},
    "cat6_minus_2": {"clip": "31S"},

    # Cat7: gene name, junction description, splice motif (RNA strand)
    "cat7_plus_1":  {"gene": "YCR012W (chrIII)", "junc_desc": "138864–138952 (88 bp)", "motif": "AC-AG"},
    "cat7_plus_2":  {"gene": "— (chrXII)",        "junc_desc": "595739–595853 (114 bp)","motif": "CA-TT"},
    "cat7_minus_1": {"gene": "YBR101C (chrII)",  "junc_desc": "443720–443833 (113 bp)","motif": "GT-CG"},
    "cat7_minus_2": {"gene": "YCL009C (chrIII)", "junc_desc": "104435–104495 (60 bp)", "motif": "GT-CG"},

    # Cat8: expected output description for each read
    "cat8_plus_single":  {"expected": "1 row, fraction=1.0"},
    "cat8_plus_multi":   {"expected": "primary 3′ at chrIV:234059 (G); fractions sum=1.0"},
    "cat8_minus_single": {"expected": "1 row, fraction=1.0"},
    "cat8_minus_multi":  {"expected": "primary 3′ at chrVIII:100520 (A ref); fractions sum=1.0"},

    # Cat9: gene name, splice motif (wrong N-op and corrected N-op derived from BAM/TSV)
    "cat9_plus_1":  {"gene": "RPL26B (chrVII)", "motif": "GT-AG"},
    "cat9_plus_2":  {"gene": "RPL30 (chrVII)",  "motif": "GT-AG"},
    "cat9_minus_1": {"gene": "RPL20B (chrXV)",  "motif": "GT-AG"},
    "cat9_minus_2": {"gene": "RPL20B (chrXV)",  "motif": "GT-AG"},
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def parse_n_ops(read):
    """Return list of (intron_start, intron_end) from N CIGAR ops (0-based)."""
    ops, ref = [], read.reference_start
    for op, length in (read.cigartuples or []):
        if op == 3:
            ops.append((ref, ref + length))
        if op in (0, 2, 3, 7, 8):
            ref += length
    return ops


def load_bam(path):
    """Return dict label → {chrom, strand, rs, re, n_ops, tags}."""
    reads = {}
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            tags = dict(read.tags)
            label = tags.get("XV", "")
            if not label:
                continue
            reads[label] = {
                "read_id": read.query_name,
                "chrom":   read.reference_name,
                "strand":  "-" if read.is_reverse else "+",
                "rs":      read.reference_start,
                "re":      read.reference_end,
                "n_ops":   parse_n_ops(read),
                "tags":    tags,
            }
    return reads


def load_tsv(path):
    """Return dict read_id → row dict."""
    with open(str(path)) as f:
        return {r["read_id"]: r for r in csv.DictReader(f, delimiter="\t")}


def match_tsv(bam_reads, tsv_rows):
    """Attach TSV row to each BAM read by UUID prefix match."""
    for rd in bam_reads.values():
        rid = rd["read_id"]
        row = tsv_rows.get(rid)
        if not row:
            for k, v in tsv_rows.items():
                if k.startswith(rid[:8]):
                    row = v
                    break
        rd["tsv"] = row or {}


# ---------------------------------------------------------------------------
# Locus helpers
# ---------------------------------------------------------------------------

def loc(chrom, start, end):
    return f"{chrom}:{start}-{end}"


def igv_link(locus, text=None):
    url = f"http://localhost:{IGV_PORT}/goto?locus={locus}"
    label = text or locus
    return f'<a class="igv" href="{url}">{label}</a>'


def td(content, cls=""):
    attr = f' class="{cls}"' if cls else ""
    return f"  <td{attr}>{content}</td>"


def strand_tag(strand):
    cls = "plus" if strand == "+" else "minus"
    sign = "+" if strand == "+" else "−"
    return f'<span class="tag {cls}">{sign}</span>'


# ---------------------------------------------------------------------------
# Per-category row renderers
# ---------------------------------------------------------------------------

def cat1_row(label, rd):
    t = rd["tsv"]
    orig  = int(t["original_3prime"])
    corr  = int(t["corrected_3prime"])
    shift = corr - orig
    sign  = "+" if shift > 0 else "−" if shift < 0 else "0"
    focused = loc(rd["chrom"], orig - 100, orig + 100)
    full    = loc(rd["chrom"], rd["rs"],   rd["re"])
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(f"{sign}{abs(shift)} bp (3′ at {orig})", "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


def cat2_row(label, rd):
    t = rd["tsv"]
    orig  = int(t["original_3prime"])
    corr  = int(t["corrected_3prime"])
    shift = corr - orig
    sign  = "+" if shift > 0 else "−" if shift < 0 else "0"
    full_end = max(rd["re"], corr + 1)
    focused  = loc(rd["chrom"], corr - 100, corr + 100)
    full     = loc(rd["chrom"], rd["rs"], full_end)
    hp       = "T"  # T-tract homopolymer (RNA poly-A = genomic T for minus; T for plus)
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(f"{sign}{abs(shift)} bp → {corr}", "detail")}\n'
        f'    {td(hp, "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


def cat3_row(label, rd):
    m   = READ_META[label]
    t   = rd["tsv"]
    raw_5p  = rd["rs"] if rd["strand"] == "+" else rd["re"] - 1
    corr_5p = int(t.get("five_prime_position", raw_5p))
    lo  = min(raw_5p, corr_5p) - 50
    hi  = max(raw_5p, corr_5p) + 50
    focused = loc(rd["chrom"], lo, hi)
    full    = loc(rd["chrom"], rd["rs"], rd["re"])
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(m["clip"], "detail")}\n'
        f'    {td(m["fp_detail"], "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


def cat4_row(label, rd):
    m = READ_META[label]
    n_ops = rd["n_ops"]
    n_s, n_e = n_ops[0] if n_ops else (rd["rs"], rd["rs"])
    focused = loc(rd["chrom"], n_s - 50, n_e + 50)
    full    = loc(rd["chrom"], rd["rs"], rd["re"])
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(m["n_desc"], "detail")}\n'
        f'    {td(m["absorbed"], "detail")}\n'
        f'    {td(m["corr_disp"], "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


def cat5_row(label, rd):
    m = READ_META[label]
    n_ops = rd["n_ops"]
    src_s, src_e = n_ops[0] if n_ops else (rd["rs"], rd["re"])
    lo  = min(rd["rs"], m["cmp_s"]) - 50
    hi  = max(rd["re"], m["cmp_e"]) + 50
    full = loc(rd["chrom"], lo, hi)
    tags = rd["tags"]
    xa   = tags.get("XA", "?")
    locus_str = rd["chrom"]
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(locus_str, "detail")}\n'
        f'    {td(xa, "detail")}\n'
        f'    {td(m["src"], "detail")}\n'
        f'    {td(m["cmp"], "detail")}\n'
        f'    {td(igv_link(full, full))}\n'
        f'  </tr>'
    )


def cat6_row(label, rd):
    m = READ_META[label]
    n_ops = rd["n_ops"]
    n_s, n_e = n_ops[0] if n_ops else (rd["rs"], rd["rs"])
    intron_desc = f'{rd["chrom"]}:{n_s}–{n_e} ({n_e - n_s} bp)'
    focused = loc(rd["chrom"], n_s - 50, n_e + 50)
    full    = loc(rd["chrom"], rd["rs"], rd["re"])
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(intron_desc, "detail")}\n'
        f'    {td(m["clip"], "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


def cat7_row(label, rd):
    m = READ_META[label]
    n_ops = rd["n_ops"]
    n_s, n_e = n_ops[0] if n_ops else (rd["rs"], rd["rs"])
    focused = loc(rd["chrom"], n_s - 200, n_e + 200)
    full    = loc(rd["chrom"], rd["rs"], rd["re"])
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(m["gene"])}\n'
        f'    {td(f'{m["junc_desc"]} · <strong>{m["motif"]}</strong>', "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


def cat8_row(label, rd):
    m = READ_META[label]
    t = rd["tsv"]
    corr    = int(t["corrected_3prime"])
    focused = loc(rd["chrom"], corr - 100, corr + 100)
    full    = loc(rd["chrom"], rd["rs"], rd["re"])
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(m["expected"], "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


def cat9_row(label, rd):
    m = READ_META[label]
    t = rd["tsv"]
    n_ops = rd["n_ops"]
    wrong_s, wrong_e = n_ops[0] if n_ops else (0, 0)
    # corrected junction from TSV junctions column: "start-end"
    junc_str = t.get("junctions", "")
    if junc_str and "-" in junc_str:
        parts = junc_str.split("-")
        corr_s, corr_e = int(parts[0]), int(parts[1])
    else:
        corr_s, corr_e = wrong_s, wrong_e
    wrong_desc = f"{wrong_s}–{wrong_e}"
    corr_desc  = f"{corr_s}–{corr_e}"
    focused = loc(rd["chrom"], corr_s - 200, corr_e + 200)
    full    = loc(rd["chrom"], rd["rs"], rd["re"])
    return (
        f'  <tr>\n'
        f'    {td(f"<code>{label}</code>")}\n'
        f'    {td(strand_tag(rd["strand"]))}\n'
        f'    {td(m["gene"])}\n'
        f'    {td(f"{wrong_desc} → {corr_desc}", "detail")}\n'
        f'    {td(m["motif"], "detail")}\n'
        f'    {td(igv_link(focused, focused))}\n'
        f'    {td(igv_link(full, "full"))}\n'
        f'  </tr>'
    )


ROW_FN = {
    "cat1": cat1_row, "cat2": cat2_row, "cat3": cat3_row,
    "cat4": cat4_row, "cat5": cat5_row, "cat6": cat6_row,
    "cat7": cat7_row, "cat8": cat8_row, "cat9": cat9_row,
}


def render_table(labels, reads, header_cols):
    th = "".join(f"<th>{c}</th>" for c in header_cols)
    rows = [f"  <tr>{th}</tr>"]
    for label in labels:
        rd  = reads[label]
        cat = label.split("_")[0]
        rows.append(ROW_FN[cat](label, rd))
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

HTML = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>IGV Links — RECTIFY Validation Reads</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; margin: 2rem; color: #1a1a1a; background: #fafafa; }}
  h1 {{ font-size: 1.4rem; margin-bottom: 0.3rem; }}
  p.subtitle {{ color: #666; margin-top: 0; font-size: 0.9rem; }}
  h2 {{ font-size: 1.1rem; margin-top: 2rem; border-bottom: 2px solid #e0e0e0; padding-bottom: 0.3rem; }}
  table {{ border-collapse: collapse; width: 100%; margin-bottom: 1.5rem; }}
  th, td {{ text-align: left; padding: 6px 12px; border: 1px solid #ddd; font-size: 0.88rem; }}
  th {{ background: #f0f0f0; font-weight: 600; }}
  tr:hover {{ background: #f5f8ff; }}
  a.igv {{ color: #0969da; text-decoration: none; font-weight: 500; }}
  a.igv:hover {{ text-decoration: underline; }}
  .tag {{ display: inline-block; font-size: 0.75rem; padding: 1px 6px; border-radius: 3px; font-weight: 600; }}
  .plus  {{ background: #dff0d8; color: #3c763d; }}
  .minus {{ background: #f2dede; color: #a94442; }}
  code {{ background: #eee; padding: 1px 4px; border-radius: 3px; font-size: 0.85rem; }}
  .note {{ background: #fff8e1; border-left: 3px solid #f9a825; padding: 0.6rem 1rem; margin-bottom: 1.5rem; font-size: 0.88rem; }}
  td.detail {{ color: #555; font-size: 0.82rem; }}
</style>
</head>
<body>

<h1>IGV Links — RECTIFY Validation Reads</h1>
<p class="subtitle">36 reads from DRS minimap2 / mapPacBio runs (S. cerevisiae, direct RNA nanopore) — nine correction categories</p>

<div class="note">
  <strong>Prerequisite:</strong> In IGV, enable the HTTP server via
  <em>View &rarr; Preferences &rarr; Advanced &rarr; Enable port</em> (default 60151).
  Load <code>validation_reads.bam</code> (and optionally <code>rectified/rectified_corrected_3end.bam</code>,
  <code>rectified/rectified_pA_tail_trimmed.bam</code>, or <code>rectified/rectified_pA_tail_soft_clipped.bam</code>)
  and select the <em>S. cerevisiae (sacCer3)</em> genome.
  The soft-clipped BAM restores the full poly(A) tail and adapter stub as soft-clips so the complete Dorado read is visible in IGV.
  <br><br>
  <strong>Corrected 3′-end bedgraphs (all categories):</strong> Load
  <code>rectified/corrected_3ends.plus.bedgraph</code> and
  <code>rectified/corrected_3ends.minus.bedgraph</code> to see
  corrected 3′-end signal for all 36 reads (Cat1–9). These tracks show a single signal peak
  at the corrected position for each read.
  Cat9 reads have N-ops refined by Module 2H — load <code>rectified/rectified_corrected_3end.bam</code>
  alongside <code>validation_reads.bam</code> to see the junction boundary shift.
</div>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 1 — poly-A walkback (<code>cat1_indel</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Reads ending in genomic A-tracts; corrected 3′ end shifts <strong>inward</strong> (+ strand: leftward; − strand: rightward).
Focused view is ±100 bp around the original 3′ end.</p>
<table>
{cat1_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 2 — soft-clip rescue at homopolymer (<code>cat2_softclip</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Aligner stopped inside a T-tract; soft-clipped bases match downstream reference. Corrected 3′ end shifts <strong>outward</strong>.
All four reads have T-tract boundaries. Focused view is ±100 bp around the corrected 3′ end.</p>
<table>
{cat2_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 3 — 5′ junction rescue (<code>cat3_junction</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Reads where the aligner's 5′ end is truncated at or inside an annotated 3′ splice site.
Four aligners (minimap2, gapmm2, deSALT, uLTRA) produce a 5′ soft-clip; mapPacBio may extend into
the intron without soft-clipping. The 5′ end is rescued to the exon boundary using the affine-gap
semi-global aligner. Focused view spans both the raw 5′ end and the rescued exon boundary.</p>
<table>
{cat3_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 4 — false N op near 3′ end (<code>cat4_false_junc</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Spurious intron (N CIGAR op) inserted by the aligner in the poly-T/A region near the 3′ end. Corrected 3′ end walks back past the N; for cat4_plus_1 and cat4_minus_2 the N is fully absorbed.
Focused view spans the spurious N op with ±50 bp padding.</p>
<table>
{cat4_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 5 — chimeric reconstruction (<code>cat5_chimeric</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Reads spanning a two-intron gene where each aligner only finds one GT-AG intron (tag <code>XK=1</code>).
Source alignment has one N op; chimeric consensus would add the complementary aligner's intron.
Tests verify <code>XA</code>, <code>XK</code>, <code>XS</code> tags and N op in source CIGAR.
Views span the full read with both introns visible (±50 bp).</p>
<table>
{cat5_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 6 — simple chimeric / single-aligner intron rescue (<code>cat6_chimeric</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Reads where mapPacBio (annotation-agnostic) correctly spans the 5′ intron while minimap2 and
gapmm2 soft-clip the same region. Stored in the BAM using the mapPacBio alignment (<code>XU=1</code>),
making the junction directly visible. Focused view is ±50 bp around the intron.</p>
<table>
{cat6_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 7 — Non-canonical, unannotated splice junctions (<code>cat7_alt_splice</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Reads with a single biologically-plausible but non-canonical splice junction aligned by mapPacBio (annotation-agnostic).
Focused view is ±200 bp around the junction. Motif = &#123;5′SS&#125;-&#123;3′SS&#125; on RNA strand.</p>
<table>
{cat7_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 8 — NET-seq A-tract refinement (<code>cat8_netseq_refine</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Reads ending in A-tracts at positions with NET-seq signal; fractional weights assigned across nearby peaks.
Focused view is ±100 bp around the corrected 3′ end.</p>
<table>
{cat8_table}
</table>

<!-- ═══════════════════════════════════════════════════════════════ -->
<h2>Category 9 — junction N-op boundary refinement (<code>cat9_junction_refine</code>)</h2>
<p style="font-size:0.88rem; color:#555;">Reads where the consensus aligner placed the intron boundary a few bp off from the true GT-AG splice site.
Module 2H (<code>--aligner-bams</code> + <code>--annotation</code>) rescores all candidate junctions within 5 kb
using HP-aware split alignment and corrects the N-op boundaries. Load <strong>both</strong>
<code>validation_reads.bam</code> (wrong boundary) and <code>rectified/rectified_corrected_3end.bam</code>
(corrected) to see the shift. Focused view is ±200 bp around the corrected junction.</p>
<table>
{cat9_table}
</table>

</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    reads = load_bam(BAM_PATH)
    tsv   = load_tsv(TSV_PATH)
    match_tsv(reads, tsv)

    def table(labels, headers):
        return render_table(labels, reads, headers)

    html = HTML.format(
        cat1_table=table(
            ["cat1_plus_1", "cat1_plus_2", "cat1_minus_1", "cat1_minus_2"],
            ["Label", "Strand", "Shift", "Focused (3′ ±100 bp)", "Full Span"],
        ),
        cat2_table=table(
            ["cat2_plus_1", "cat2_plus_2", "cat2_minus_1", "cat2_minus_2"],
            ["Label", "Strand", "Shift / corrected 3′", "Homopolymer", "Focused (corrected ±100 bp)", "Full Span"],
        ),
        cat3_table=table(
            ["cat3_plus_1", "cat3_plus_2", "cat3_minus_1", "cat3_minus_2"],
            ["Label", "Strand", "5′ soft-clip", "Raw 5′ end", "Focused", "Full Span"],
        ),
        cat4_table=table(
            ["cat4_plus_1", "cat4_plus_2", "cat4_minus_1", "cat4_minus_2"],
            ["Label", "Strand", "Spurious N op", "N absorbed?", "Corrected 3′", "Focused (N ±50 bp)", "Full Span"],
        ),
        cat5_table=table(
            ["cat5_plus_1", "cat5_plus_2", "cat5_minus_1", "cat5_minus_2"],
            ["Label", "Strand", "Locus", "Aligners (XA)", "Source intron", "Complementary intron", "Full Span ±50 bp"],
        ),
        cat6_table=table(
            ["cat6_plus_1", "cat6_plus_2", "cat6_minus_1", "cat6_minus_2"],
            ["Label", "Strand", "Intron (mapPacBio)", "mm2 soft-clip", "Focused (intron ±50 bp)", "Full Span"],
        ),
        cat7_table=table(
            ["cat7_plus_1", "cat7_plus_2", "cat7_minus_1", "cat7_minus_2"],
            ["Label", "Strand", "Gene", "Junction / Motif", "Focused (junction ±200 bp)", "Full Span"],
        ),
        cat8_table=table(
            ["cat8_plus_single", "cat8_plus_multi", "cat8_minus_single", "cat8_minus_multi"],
            ["Label", "Strand", "Expected output", "Focused (3′ ±100 bp)", "Full Span"],
        ),
        cat9_table=table(
            ["cat9_plus_1", "cat9_plus_2", "cat9_minus_1", "cat9_minus_2"],
            ["Label", "Strand", "Gene", "Wrong N-op → Corrected", "Motif", "Focused (junction ±200 bp)", "Full Span"],
        ),
    )

    HTML_PATH.write_text(html)
    print(f"Written: {HTML_PATH}")


if __name__ == "__main__":
    main()
