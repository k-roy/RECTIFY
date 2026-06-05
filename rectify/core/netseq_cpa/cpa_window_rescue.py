"""Tier 2 — targeted CPA-window mini-reference rescue (residual unmapped anchors).

For anchors too short to place uniquely genome-wide (below the k-mer uniqueness
floor; yeast: 12 bp -> 18%, 16 bp -> 90% unique), align ONLY against a small
candidate set: a mini-reference of CPA windows. An anchor ambiguous genome-wide
becomes unique within the set -> place it -> reveal the full (possibly long)
poly-A tail the genome-wide step cannot see.

The mini-reference = one contig per CPA, the genomic window spanning
``gene_bodyward`` bp into the gene body and ``downstream`` bp past the cleavage
site (default 50/10):

  + gene: genome [modal - gene_bodyward, modal + downstream]
  - gene: genome [modal - downstream,    modal + gene_bodyward]

so the contig always carries the gene-body-ward anchor sequence (where an
intermediate's genomic portion lies) plus a short downstream margin. Each contig
header encodes ``chrom|strand|modal|win_start`` so a within-contig placement maps
back to genome coordinates and the CPA.

⚠ CIRCULARITY GUARD (non-negotiable): Tier-2 places reads AT known CPAs, so
Tier-2 reads are CIRCULAR for site discovery, DRS-concordance, and the novelty
claim. **Tier 2 is for poly-A LENGTH characterization at already-established
CPAs ONLY.** Placements are tagged ``tier=2``; callers must never pool them into
discovery/concordance numbers. :func:`build_shuffled_window_reference` builds a
coordinate-shuffled null to bound chance placement (report the FP rate).

bowtie2 lives in the dedicated ``netseq_align`` env (the rectify env is
untouched); pass its binary paths via ``bowtie2_bin`` / ``bowtie2_build_bin``.
This module is import-safe without bowtie2; only :func:`align_anchors_to_windows`
invokes it.
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple

import pysam


def _windows_from_clusters(
    drs_clusters: str | Path,
    gene_bodyward: int,
    downstream: int,
) -> List[Tuple[str, str, int, int, int]]:
    """Return ``[(chrom, strand, modal, win_start, win_end), ...]`` (0-based,
    end-exclusive) from a DRS cluster TSV."""
    out = []
    with open(drs_clusters) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            chrom, strand = r["chrom"], r["strand"]
            modal = int(float(r["modal_position"]))
            if strand == "+":
                ws, we = modal - gene_bodyward, modal + downstream + 1
            else:
                ws, we = modal - downstream, modal + gene_bodyward + 1
            out.append((chrom, strand, modal, max(0, ws), we))
    return out


def build_cpa_window_reference(
    drs_clusters: str | Path,
    genome_path: str | Path,
    out_fa: str | Path,
    *,
    gene_bodyward: int = 50,
    downstream: int = 10,
    chrom_map: Optional[Mapping[str, str]] = None,
    extra_sites: Optional[List[Tuple[str, str, int]]] = None,
) -> int:
    """Write the CPA-window mini-reference FASTA. Returns the number of contigs.

    Each contig is the genomic ``+``-strand sequence of the window; the header
    ``{chrom}|{strand}|{modal}|{win_start}`` lets a placement be mapped back to
    genome coordinates. ``chrom_map`` translates the cluster ``chrom`` to the
    genome FASTA's reference name if they differ. ``extra_sites`` (e.g. cDNA /
    Tier-1-supported CPAs) are appended as ``(chrom, strand, modal)`` tuples.
    """
    rev_map = {v: k for k, v in dict(chrom_map or {}).items()}
    fa = pysam.FastaFile(str(genome_path))
    windows = _windows_from_clusters(drs_clusters, gene_bodyward, downstream)
    if extra_sites:
        for chrom, strand, modal in extra_sites:
            if strand == "+":
                ws, we = modal - gene_bodyward, modal + downstream + 1
            else:
                ws, we = modal - downstream, modal + gene_bodyward + 1
            windows.append((chrom, strand, modal, max(0, ws), we))

    n = 0
    with open(out_fa, "w") as out:
        for chrom, strand, modal, ws, we in windows:
            ref_name = rev_map.get(chrom, chrom)
            try:
                seq = fa.fetch(ref_name, ws, we)
            except (KeyError, ValueError):
                continue
            if not seq:
                continue
            out.write(f">{chrom}|{strand}|{modal}|{ws}\n{seq}\n")
            n += 1
    fa.close()
    return n


def build_shuffled_window_reference(
    drs_clusters: str | Path,
    genome_path: str | Path,
    out_fa: str | Path,
    *,
    gene_bodyward: int = 50,
    downstream: int = 10,
    shift: int = 5000,
    chrom_map: Optional[Mapping[str, str]] = None,
) -> int:
    """Coordinate-shuffled null: the same number/size of windows, each shifted
    ``shift`` bp from a real CPA, to bound chance placement. Returns contig count.
    """
    rev_map = {v: k for k, v in dict(chrom_map or {}).items()}
    fa = pysam.FastaFile(str(genome_path))
    windows = _windows_from_clusters(drs_clusters, gene_bodyward, downstream)
    n = 0
    with open(out_fa, "w") as out:
        for chrom, strand, modal, ws, we in windows:
            ref_name = rev_map.get(chrom, chrom)
            ws2, we2 = ws + shift, we + shift
            try:
                seq = fa.fetch(ref_name, ws2, we2)
            except (KeyError, ValueError):
                continue
            if not seq:
                continue
            out.write(f">{chrom}|{strand}|{modal}|{ws2}|shuffled\n{seq}\n")
            n += 1
    fa.close()
    return n


def align_anchors_to_windows(
    anchors_fq: str | Path,
    window_fa: str | Path,
    out_bam: str | Path,
    *,
    bowtie2_bin: str = "bowtie2",
    bowtie2_build_bin: str = "bowtie2-build",
    samtools_bin: str = "samtools",
    threads: int = 4,
    min_anchor: int = 12,
) -> str:
    """Build a bowtie2 index over ``window_fa`` and align ``anchors_fq`` to it.

    Uses ``--norc``? No — anchors may be either strand, so default (both) is
    kept; ``-L`` seed is dropped to ``min_anchor`` and ``-k 2`` reports up to two
    hits so within-set multiplicity can be checked. Returns the sorted BAM path.
    Requires bowtie2 (``netseq_align`` env); not invoked at import.
    """
    window_fa = str(window_fa)
    idx = window_fa + ".idx"
    subprocess.run(
        [bowtie2_build_bin, "--quiet", window_fa, idx], check=True
    )
    seed = max(6, min(min_anchor, 20))
    bt2 = subprocess.Popen(
        [
            bowtie2_bin, "-x", idx, "-U", str(anchors_fq),
            "-L", str(seed), "-N", "0", "-k", "2",
            "--score-min", "L,-0.6,-0.6", "--no-unal",
            "-p", str(threads),
        ],
        stdout=subprocess.PIPE,
    )
    with open(out_bam, "wb") as fh:
        sort = subprocess.Popen(
            [samtools_bin, "sort", "-@", str(threads), "-o", "-", "-"],
            stdin=bt2.stdout, stdout=fh,
        )
        bt2.stdout.close()
        sort.communicate()
    if bt2.wait() != 0:
        raise RuntimeError("bowtie2 alignment failed")
    subprocess.run([samtools_bin, "index", str(out_bam)], check=True)
    return str(out_bam)


def quantify_tier2(
    bam_path: str | Path,
    *,
    reg: int = 5,
    min_mapq: int = 1,
) -> Dict[str, object]:
    """Quantify Tier-2 placements: poly-A length distribution at the windows'
    CPAs. CIRCULAR — length characterization only, never discovery.

    A placement's CPA is the window's encoded ``modal``; a placement counts if
    the anchor's poly-T-side end falls within ``reg`` of ``modal`` (i.e. the
    anchor really ends at the CPA, not elsewhere in the window). Returns
    ``{placed, at_cpa, palen_atcpa_median/mean, tier: 2, circular: True}``.
    """
    import re
    import statistics as st

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    placed = at_cpa = 0
    palen_atcpa: List[int] = []
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        mm = re.search(r"_pa(\d+)$", read.query_name)
        if not mm:
            continue
        L = int(mm.group(1))
        parts = read.reference_name.split("|")
        if len(parts) < 4:
            continue
        _chrom, strand, modal_s, ws_s = parts[0], parts[1], parts[2], parts[3]
        modal, ws = int(modal_s), int(ws_s)
        placed += 1
        # The anchor's RNA-3' end within the window, mapped back to genome coord.
        if read.is_reverse:
            genome_pos = ws + (read.reference_end - 1)
        else:
            genome_pos = ws + read.reference_start
        if abs(genome_pos - modal) <= reg:
            at_cpa += 1
            palen_atcpa.append(min(L, 200))
    bam.close()
    return {
        "placed": placed,
        "at_cpa": at_cpa,
        "palen_atcpa_median": st.median(palen_atcpa) if palen_atcpa else 0,
        "palen_atcpa_mean": st.mean(palen_atcpa) if palen_atcpa else 0.0,
        "tier": 2,
        "circular": True,
    }
