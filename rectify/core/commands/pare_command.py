"""``rectify pare`` — 5'P PARE / degradome processing orchestrator.

Processes 5'-monophosphate PARE (degradome) libraries from raw FASTQ (or a
local-mode BAM that retains unmapped reads): the strand-aware 5'P cut-site
pileup + bedgraphs, the mapped read-length census, and — the discovery arm —
soft-clipped poly(A) tails at the RNA 3' side as evidence of endonucleolytic
cleavage + polyadenylation (CPA). PARE is SENSE (read base 1 = the 5'P cut
site); see :mod:`rectify.core.pare`.

Pipeline (mirrors the validated ``netseq-cpa`` arm):

  [1] (FASTQ in)   3' adapter trim                 bbduk
  [2] (FASTQ in)   align, keep unmapped            bbmap local=t
  [3] Tier-1 single pass (mapped reads)            pare.pileup
        5'P pileup (clip-gated) + bedgraphs, read-length census,
        genome-aware poly-A CPA pileup, per-read parquet incl. 5'P
  [4] mapped CPA concordance vs DRS map            netseq_cpa.concordance
  [5] rescue: trim 3' poly-A -> remap -> CPA       pare.rescue      (unmapped)
  [6] summary.json + PROVENANCE.json

BAM input is assumed to already retain unmapped reads (bbmap local=t output);
steps [1]-[2] are skipped. ``--drs-clusters`` enables [4] and the at-CPA
rescue fraction.

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

from .. import pare as pare_mod
from ..netseq_cpa import concordance as _conc
from ..pare import pileup as _pileup
from ..pare import rescue as _rescue
from ..pare.reads_parquet import PareReadWriter


def _log(msg: str) -> None:
    print(f"[pare] {msg}", flush=True)


def _run(cmd: list, **kw) -> None:
    _log("$ " + " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kw)


def _bbtool(bin_name: str, mem: Optional[str]) -> list:
    """Start a bb-tool command, capping the JVM heap when ``mem`` is set (see
    the netseq-cpa arm: auto-detected heaps overshoot on cgroup-limited nodes)."""
    cmd = [bin_name]
    if mem:
        cmd.append(f"-Xmx{mem}")
    return cmd


def _resolve_chrom_map(name: Optional[str]):
    if name in (None, "none", "identity"):
        return None
    if name in ("ncbi_yeast", "scer", "yeast"):
        return pare_mod.NCBI_TO_CHROM
    raise SystemExit(f"Unknown --chrom-map '{name}' (use none | ncbi_yeast)")


def _is_bam(path: Path) -> bool:
    return path.suffix.lower() in (".bam", ".cram", ".sam")


def run(args: argparse.Namespace) -> int:
    t0 = time.time()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    genome = Path(args.genome)
    chrom_map = _resolve_chrom_map(getattr(args, "chrom_map", None))
    drs = Path(args.drs_clusters) if getattr(args, "drs_clusters", None) else None
    threads = int(getattr(args, "threads", 4) or 4)
    label = getattr(args, "label", None) or out.name
    summary: dict = {"label": label, "input": str(args.input)}

    inp = Path(args.input)

    # ----------------------------------------------------------- [1][2] align
    if _is_bam(inp):
        aligned = inp
        _log(f"input is a BAM; using as aligned (unmapped reads must be retained): {aligned}")
    else:
        adapter = getattr(args, "adapter", None)
        if not adapter:
            raise SystemExit("--adapter is required for FASTQ input")
        trimmed = out / "trimmed.fq.gz"
        k = min(len(adapter), 23)
        _log(f"[1] bbduk 3' adapter trim (k={k})")
        bbduk = _bbtool(getattr(args, "bbduk_bin", "bbduk.sh"), getattr(args, "bbmap_mem", None)) + [
            f"in={inp}", f"out={trimmed}", f"literal={adapter}",
            "ktrim=r", f"k={k}", f"mink={min(11, k - 1)}", "hdist=1",
            f"minlen={int(getattr(args, 'minlen', 18))}",
            f"threads={threads}", "overwrite=t",
        ]
        _run(bbduk)
        aligned = out / "aligned.bam"
        _log("[2] bbmap local=t (keep unmapped)")
        bbmap = _bbtool(getattr(args, "bbmap_bin", "bbmap.sh"), getattr(args, "bbmap_mem", None)) + [
            f"in={trimmed}", f"out={aligned}", f"ref={genome}",
            "local=t", f"threads={threads}", "overwrite=t",
        ]
        if getattr(args, "bbmap_index", None):
            bbmap.append(f"path={args.bbmap_index}")
        _run(bbmap)

    # ------------------------------------------- per-read parquet sidecar
    reads_parquet_path = out / "reads.parquet"
    cpa_set = _conc.load_cpa_set(drs) if drs else None
    rpw = None
    if getattr(args, "reads_parquet", True):
        rpw = PareReadWriter(reads_parquet_path, log=_log)

    # ------------------------------------------------------- [3] Tier-1 pass
    _log("[3] Tier-1: 5'P pileup + read-length census + poly-A CPA pileup")
    pstats = _pileup.pare_pileup(
        aligned, genome, out, label, chrom_map=chrom_map,
        limit=int(getattr(args, "limit", 0) or 0),
        max_fivep_clip=int(getattr(args, "max_fivep_clip",
                                   _pileup.MAX_FIVEP_CLIP_DEFAULT)),
        reads_parquet=rpw, cpa_set=cpa_set,
    )
    _log(f"    {pstats}")
    _log(f"    5'P clip-excluded fraction: {100 * pstats['fivep_clip_fraction']:.2f}% "
         "(reads whose RNA-5' terminus was clipped; invariant: report always)")
    summary["tier1"] = pstats

    # ------------------------------------------------------- [4] concordance
    cpa_pileup_path = out / "cpa_pileup.tsv.gz"
    if drs:
        _log("[4] mapped CPA concordance vs DRS map")
        cstats = _conc.mapped_cpa_concordance(
            cpa_pileup_path, drs, label=label,
            nuclear_bp=int(getattr(args, "nuclear_bp", pare_mod.SCER_NUCLEAR_BP)),
        )
        _log(f"    at CPA {100*cstats['frac_at_cpa']:.1f}%  ({cstats['enrichment']:.1f}x null)"
             f"  mean poly-A {cstats['mean_polya_len']:.1f} nt")
        summary["concordance"] = cstats

    # ------------------------------------------------------------- [5] rescue
    anchors_fq = out / "anchors.fq.gz"
    _log("[5] rescue: trim 3' poly-A from unmapped reads")
    rstats = _rescue.trim_unmapped_polya(
        aligned, anchors_fq,
        min_pa=int(getattr(args, "min_pa", _rescue.MIN_PA_DEFAULT)),
        min_anchor=int(getattr(args, "min_anchor", _rescue.MIN_ANCHOR_DEFAULT)),
    )
    _log(f"    unmapped={rstats['unmapped']} rescued_anchors={rstats['rescued']}")
    summary["rescue_trim"] = rstats
    anchors_bam = out / "anchors.bam"
    if rstats["rescued"] > 0:
        _log("    re-mapping rescue anchors (bbmap)")
        bbmap2 = _bbtool(getattr(args, "bbmap_bin", "bbmap.sh"), getattr(args, "bbmap_mem", None)) + [
            f"in={anchors_fq}", f"out={anchors_bam}", f"ref={genome}",
            "local=t", "maxindel=5", "ambiguous=random",
            f"threads={threads}", "overwrite=t",
        ]
        if getattr(args, "bbmap_index", None):
            bbmap2.append(f"path={args.bbmap_index}")
        _run(bbmap2)
        qstats = _rescue.quantify_rescue(
            anchors_bam, drs, min_mapq=int(getattr(args, "min_mapq", 3)),
            chrom_map=chrom_map, reads_parquet=rpw,
        )
        _log(f"    rescued mapped={qstats['mapped']} at_CPA={qstats['at_cpa']} "
             f"({100*qstats['frac_at_cpa']:.1f}%) poly-A median={qstats['palen_median']}")
        summary["rescue_quantify"] = qstats

    # --------------------------------------------- close per-read parquet
    if rpw is not None:
        n_rows = rpw.close()
        summary["reads_parquet"] = {
            "path": str(reads_parquet_path), "rows": n_rows,
            "disabled": rpw.disabled,
        }
        if rpw.disabled:
            _log("    reads.parquet: DISABLED (pyarrow missing or a write failed) "
                 "— pileups are unaffected")
        else:
            _log(f"    reads.parquet: {n_rows} per-read rows -> {reads_parquet_path.name}")

    # --------------------------------------------------------- [6] summary
    summary["elapsed_sec"] = round(time.time() - t0, 1)
    (out / "pare_summary.json").write_text(json.dumps(summary, indent=2))
    try:
        from rectify.utils.provenance import ProvenanceTracker
        tr = ProvenanceTracker(out, description="RECTIFY pare 5'P + poly-A CPA capture")
        tr.set_command(sys.argv)
        tr.set_config(vars(args))
        for f in (cpa_pileup_path, out / "fivep_pileup.tsv.gz",
                  out / "read_lengths.tsv.gz", reads_parquet_path, anchors_fq,
                  out / "pare_summary.json"):
            if Path(f).exists():
                tr.add_output_file(Path(f), source_files=[Path(args.input)])
        tr.save()
    except Exception as e:  # provenance is best-effort, never fail the run
        _log(f"provenance sidecar skipped: {e}")
    _log(f"DONE {label} in {summary['elapsed_sec']}s -> {out}")
    return 0


def create_pare_parser(subparsers) -> None:
    """Wire the ``pare`` subcommand."""
    p = subparsers.add_parser(
        "pare",
        help="Process 5'P PARE/degradome reads: 5'P pileup + soft-clipped poly(A) CPA capture",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "PARE (5'-monophosphate degradome) arm: strand-aware 5'P cut-site "
            "pileup + bedgraphs (clip-gated), mapped read-length census, and "
            "genome-aware capture of soft-clipped poly(A) tails at the RNA 3' "
            "side as CPA evidence — Tier-1 on mapped reads plus rescue of "
            "poly-A-dominated unmapped reads (trim 3' poly-A -> remap -> CPA). "
            "Optionally validated against the DRS CPA map."
        ),
    )
    p.add_argument("input", type=Path,
                   help="PARE FASTQ.gz (single-end, 3' adapter still on) OR a "
                        "bbmap local=t BAM that retains unmapped reads.")
    p.add_argument("--genome", type=Path, required=True, help="Reference genome FASTA.")
    p.add_argument("-o", "--output-dir", type=Path, required=True, dest="output_dir")
    p.add_argument("--adapter", default=None,
                   help="3' adapter literal for bbduk (FASTQ input). comPARE "
                        "PRJNA663967: TGGAATTCTCGGGTGCCAAGG; PRJNA1330880: AGATCGGAAGAGC.")
    p.add_argument("--label", default=None, help="Sample label (default: output dir name).")
    p.add_argument("--drs-clusters", type=Path, default=None, dest="drs_clusters",
                   help="DRS CPA cluster TSV (chrom, strand, modal_position) for "
                        "concordance + at-CPA rescue. Optional but recommended.")
    p.add_argument("--chrom-map", default="none", dest="chrom_map",
                   choices=["none", "ncbi_yeast"],
                   help="Translate BAM reference names to chrI-style (ncbi_yeast) or identity (none).")
    p.add_argument("--max-fivep-clip", type=int, default=_pileup.MAX_FIVEP_CLIP_DEFAULT,
                   dest="max_fivep_clip",
                   help="Max RNA-5'-side clip for a read to enter the 5'P pileup "
                        "(workspace invariant: 0; the excluded fraction is always reported).")
    p.add_argument("--minlen", type=int, default=18,
                   help="bbduk minlen= (min read length after adapter trim).")
    p.add_argument("--min-pa", type=int, default=_rescue.MIN_PA_DEFAULT, dest="min_pa",
                   help="Min 3' poly-A length to attempt rescue of an unmapped read.")
    p.add_argument("--min-anchor", type=int, default=_rescue.MIN_ANCHOR_DEFAULT, dest="min_anchor",
                   help="Min genomic anchor length to keep a rescued read.")
    p.add_argument("--min-mapq", type=int, default=3, dest="min_mapq",
                   help="Min MAPQ for re-mapped rescue anchors.")
    p.add_argument("--nuclear-bp", type=int, default=pare_mod.SCER_NUCLEAR_BP, dest="nuclear_bp",
                   help="Nuclear genome bp for the concordance genome null (default: S. cerevisiae R64).")
    p.add_argument("--limit", type=int, default=0, help="Cap reads processed in Tier-1 (0 = all).")
    p.add_argument("--no-reads-parquet", action="store_false", dest="reads_parquet",
                   default=True,
                   help="Disable the per-read reads.parquet sidecar (one row per "
                        "Tier-1/rescued read: read_id, chrom, cpa_pos, gene_strand, "
                        "oaNT_tail_len, at_cpa, tier, mapq, five_p_pos, five_p_clip). "
                        "On by default; needs pyarrow.")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--bbduk-bin", default="bbduk.sh", dest="bbduk_bin")
    p.add_argument("--bbmap-bin", default="bbmap.sh", dest="bbmap_bin")
    p.add_argument("--bbmap-index", default=None, dest="bbmap_index",
                   help="Pre-built bbmap index dir (path=...). Built on the fly if omitted.")
    p.add_argument("--bbmap-mem", default=None, dest="bbmap_mem",
                   help="Cap the bbmap/bbduk JVM heap (e.g. 8g). Forwards -Xmx. "
                        "Set on login / cgroup-limited nodes; omit on dedicated compute nodes.")
