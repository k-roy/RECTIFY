"""``rectify netseq-cpa`` — NET-seq CPA-intermediate capture orchestrator.

Recovers cleavage-and-polyadenylation (CPA) intermediates from NET-seq /
mNET-seq libraries: Pol II-associated nascent RNA 3' ends that map to poly-A
sites and carry a non-templated poly(A) tail. Standard NET-seq pipelines discard
these (poly-A-dominated reads fail the aligner min-ratio filter); this arm
recovers them via a two-tier strategy. See :mod:`rectify.core.netseq_cpa`.

Pipeline (mirrors the validated ``process_netseq.sbatch``, productionized):

  [0] (paired-end) isolate read 2          pe_reduce
  [1] (FASTQ in)   adapter-trim            bbduk
  [2] (FASTQ in)   align, keep unmapped    bbmap local=t
  [3] Tier-1 genome-aware poly-A pileup    netseq_cpa.pileup        (mapped)
  [4] mapped CPA concordance vs DRS map    netseq_cpa.concordance
  [5] rescue: trim poly-T -> remap -> CPA  netseq_cpa.rescue        (unmapped)
  [6] (--tier2) CPA-window mini-ref rescue netseq_cpa.cpa_window_rescue  (CIRCULAR)
  [7] summary.json + PROVENANCE.json

BAM input is assumed to already retain unmapped reads (bbmap local=t output);
the alignment steps [1]-[2] are skipped. ``--drs-clusters`` enables steps [4]/[6]
and the at-CPA rescue fraction.

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

from .. import netseq_cpa
from ..netseq_cpa import concordance as _conc
from ..netseq_cpa import pileup as _pileup
from ..netseq_cpa import rescue as _rescue
from ..netseq_cpa import reads_parquet as _rpq


def _log(msg: str) -> None:
    print(f"[netseq-cpa] {msg}", flush=True)


def _run(cmd: list, **kw) -> None:
    _log("$ " + " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kw)


def _bbtool(bin_name: str, mem: Optional[str]) -> list:
    """Start a bb-tool command, capping the JVM heap when ``mem`` is set.

    bbmap/bbduk auto-detect heap from physical RAM, which overshoots on login or
    cgroup-limited nodes (it requests a multi-GB heap the ulimit rejects -> "Unable
    to allocate ... heap"). ``--bbmap-mem 4g`` forwards ``-Xmx4g`` to cap it.
    """
    cmd = [bin_name]
    if mem:
        cmd.append(f"-Xmx{mem}")
    return cmd


def _resolve_chrom_map(name: Optional[str]):
    if name in (None, "none", "identity"):
        return None
    if name in ("ncbi_yeast", "scer", "yeast"):
        return netseq_cpa.NCBI_TO_CHROM
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

    # ----------------------------------------------------------------- [0] PE
    if getattr(args, "paired_end", False) and _is_bam(inp):
        from ..netseq_cpa import pe_reduce
        se_bam = out / "read2_single_end.bam"
        _log(f"[0] paired-end: isolating read 2 -> {se_bam.name}")
        stats = pe_reduce.reduce_bam_to_read2(inp, se_bam)
        _log(f"    read2 written: {stats['written']}/{stats['seen']}")
        _run(["samtools", "index", str(se_bam)])
        inp = se_bam
        summary["pe_reduce"] = stats

    # ----------------------------------------------------------- [1][2] align
    if _is_bam(inp):
        aligned = inp
        _log(f"input is a BAM; using as aligned (unmapped reads must be retained): {aligned}")
    else:
        if getattr(args, "paired_end", False):
            from ..netseq_cpa import pe_reduce
            inp = pe_reduce.select_read2_fastq(getattr(args, "read1", None), getattr(args, "read2", None) or inp)
            _log(f"paired-end FASTQ: using read 2 = {inp}")
        trimmed = out / "trimmed.fq.gz"
        adapter = getattr(args, "adapter", None)
        if not adapter:
            raise SystemExit("--adapter is required for FASTQ input")
        _log("[1] bbduk adapter trim")
        bbduk = _bbtool(getattr(args, "bbduk_bin", "bbduk.sh"), getattr(args, "bbmap_mem", None)) + [
            f"in={inp}", f"out={trimmed}", f"literal={adapter}",
            "ktrim=r", "k=22", "mink=11", "hdist=1", "minlen=15", "tbo",
            f"threads={threads}", "overwrite=t",
        ]
        utrim = int(getattr(args, "utrim", 0) or 0)
        if utrim > 0:
            bbduk.append(f"ftl={utrim}")
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

    # ----------------------------------------- per-read parquet sidecar (unified)
    # one row per CPA-calling read (tier mapped|rescued); spans [3] + [5]. The CPA
    # set is the single source of truth shared with concordance so at_cpa agrees.
    reads_parquet_path = out / "reads.parquet"
    cpa_set = _conc.load_cpa_set(drs) if drs else None
    rpw = None
    if getattr(args, "reads_parquet", True):
        rpw = _rpq.ParquetReadWriter(reads_parquet_path, log=_log)

    # --------------------------------------------------------- [3] Tier-1 pileup
    pileup_path = out / "pileup.tsv.gz"
    _log("[3] Tier-1 genome-aware poly-A pileup (mapped reads)")
    pstats = _pileup.walkback_pileup(
        aligned, genome, pileup_path, label, chrom_map=chrom_map,
        limit=int(getattr(args, "limit", 0) or 0),
        reads_parquet=rpw, cpa_set=cpa_set,
    )
    _log(f"    {pstats}")
    summary["tier1_pileup"] = pstats

    # ------------------------------------------------------- [4] concordance
    if drs:
        _log("[4] mapped CPA concordance vs DRS map")
        cstats = _conc.mapped_cpa_concordance(
            pileup_path, drs, label=label,
            nuclear_bp=int(getattr(args, "nuclear_bp", netseq_cpa.SCER_NUCLEAR_BP)),
        )
        _log(f"    at CPA {100*cstats['frac_at_cpa']:.1f}%  ({cstats['enrichment']:.1f}x null)"
             f"  mean poly-A {cstats['mean_polya_len']:.1f} nt")
        summary["concordance"] = cstats

    # ------------------------------------------------------------- [5] rescue
    anchors_fq = out / "anchors.fq.gz"
    _log("[5] rescue: trim 5' poly-T from unmapped reads")
    rstats = _rescue.trim_unmapped_polyt(
        aligned, anchors_fq,
        min_pt=int(getattr(args, "min_pt", _rescue.MIN_PT_DEFAULT)),
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
        if drs:
            qstats = _rescue.quantify_rescue(
                anchors_bam, drs, min_mapq=int(getattr(args, "min_mapq", 3)),
                chrom_map=chrom_map, reads_parquet=rpw,
            )
            _log(f"    rescued mapped={qstats['mapped']} at_CPA={qstats['at_cpa']} "
                 f"({100*qstats['frac_at_cpa']:.1f}%) poly-A median={qstats['palen_median']}")
            summary["rescue_quantify"] = qstats

    # --------------------------------------------- close per-read parquet sidecar
    if rpw is not None:
        n_rows = rpw.close()
        summary["reads_parquet"] = {
            "path": str(reads_parquet_path), "rows": n_rows,
            "disabled": rpw.disabled,
        }
        if rpw.disabled:
            _log("    reads.parquet: DISABLED (pyarrow missing or a write failed) "
                 "— pileup.tsv.gz is unaffected")
        else:
            _log(f"    reads.parquet: {n_rows} per-read rows -> {reads_parquet_path.name}")

    # --------------------------------------------------------- [6] Tier-2 (opt)
    if getattr(args, "tier2", False):
        if not drs:
            _log("[6] Tier-2 skipped: --tier2 requires --drs-clusters")
        else:
            from ..netseq_cpa import cpa_window_rescue as _t2
            _log("[6] Tier-2 CPA-window mini-ref rescue (CIRCULAR — length only)")
            win_fa = out / "cpa_windows.fa"
            n_win = _t2.build_cpa_window_reference(
                drs, genome, win_fa, chrom_map=chrom_map,
                gene_bodyward=int(getattr(args, "tier2_bodyward", 50)),
                downstream=int(getattr(args, "tier2_downstream", 10)),
            )
            _log(f"    built {n_win} CPA-window contigs")
            t2_bam = out / "tier2_anchors.bam"
            _t2.align_anchors_to_windows(
                anchors_fq, win_fa, t2_bam,
                bowtie2_bin=getattr(args, "bowtie2_bin", "bowtie2"),
                bowtie2_build_bin=getattr(args, "bowtie2_build_bin", "bowtie2-build"),
                threads=threads, min_anchor=int(getattr(args, "tier2_min_anchor", 12)),
            )
            t2stats = _t2.quantify_tier2(t2_bam)
            _log(f"    Tier-2 placed={t2stats['placed']} at_CPA={t2stats['at_cpa']} "
                 f"poly-A median={t2stats['palen_atcpa_median']} (CIRCULAR)")
            summary["tier2"] = t2stats

    # ------------------------------------------------------- [7] summary + prov
    summary["elapsed_sec"] = round(time.time() - t0, 1)
    (out / "netseq_cpa_summary.json").write_text(json.dumps(summary, indent=2))
    try:
        from rectify.utils.provenance import ProvenanceTracker
        tr = ProvenanceTracker(out, description="RECTIFY netseq-cpa CPA-intermediate capture")
        tr.set_command(sys.argv)
        tr.set_config(vars(args))
        for f in (pileup_path, reads_parquet_path, anchors_fq,
                  out / "netseq_cpa_summary.json"):
            if Path(f).exists():
                tr.add_output_file(Path(f), source_files=[Path(args.input)])
        tr.save()
    except Exception as e:  # provenance is best-effort, never fail the run
        _log(f"provenance sidecar skipped: {e}")
    _log(f"DONE {label} in {summary['elapsed_sec']}s -> {out}")
    return 0


def create_netseq_cpa_parser(subparsers) -> None:
    """Wire the ``netseq-cpa`` subcommand."""
    p = subparsers.add_parser(
        "netseq-cpa",
        help="Recover CPA cleavage intermediates from NET-seq via soft-clipped oligo(A)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Two-tier poly-A-aware NET-seq recovery: Tier-1 genome-aware poly-A "
            "pileup on mapped reads + rescue of poly-A-dominated unmapped reads "
            "(trim 5' poly-T -> remap -> CPA), validated against the DRS CPA map. "
            "Optional Tier-2 CPA-window rescue is CIRCULAR (poly-A length only)."
        ),
    )
    p.add_argument("input", type=Path,
                   help="NET-seq FASTQ.gz (single-end, or read 2 with --paired-end) "
                        "OR a bbmap local=t BAM that retains unmapped reads.")
    p.add_argument("--genome", type=Path, required=True, help="Reference genome FASTA.")
    p.add_argument("-o", "--output-dir", type=Path, required=True, dest="output_dir")
    p.add_argument("--drs-clusters", type=Path, default=None, dest="drs_clusters",
                   help="DRS CPA cluster TSV (chrom, strand, modal_position) for "
                        "concordance + at-CPA rescue + Tier-2. Optional but recommended.")
    p.add_argument("--label", default=None, help="Sample label (default: output dir name).")
    p.add_argument("--adapter", default=None, help="3' adapter literal for bbduk (FASTQ input).")
    p.add_argument("--bbmap-index", default=None, dest="bbmap_index",
                   help="Pre-built bbmap index dir (path=...). Built on the fly if omitted.")
    p.add_argument("--paired-end", action="store_true", dest="paired_end",
                   help="mNET-seq paired-end: isolate read 2 (nascent 3' end) first.")
    p.add_argument("--read1", type=Path, default=None, help="R1 FASTQ (paired-end FASTQ input).")
    p.add_argument("--read2", type=Path, default=None, help="R2 FASTQ (paired-end FASTQ input).")
    p.add_argument("--chrom-map", default="none", dest="chrom_map",
                   choices=["none", "ncbi_yeast"],
                   help="Translate BAM reference names to chrI-style (ncbi_yeast) or identity (none).")
    p.add_argument("--nuclear-bp", type=int, default=netseq_cpa.SCER_NUCLEAR_BP, dest="nuclear_bp",
                   help="Nuclear genome bp for the concordance genome null (default: S. cerevisiae R64).")
    p.add_argument("--utrim", type=int, default=0, help="bbduk ftl= (trim N leading bases; UMI).")
    p.add_argument("--min-pt", type=int, default=_rescue.MIN_PT_DEFAULT, dest="min_pt",
                   help="Min 5' poly-T length to attempt rescue.")
    p.add_argument("--min-anchor", type=int, default=_rescue.MIN_ANCHOR_DEFAULT, dest="min_anchor",
                   help="Min genomic anchor length to keep a rescued read.")
    p.add_argument("--min-mapq", type=int, default=3, dest="min_mapq",
                   help="Min MAPQ for re-mapped rescue anchors.")
    p.add_argument("--limit", type=int, default=0, help="Cap reads processed in Tier-1 (0 = all).")
    p.add_argument("--no-reads-parquet", action="store_false", dest="reads_parquet",
                   default=True,
                   help="Disable the per-read reads.parquet sidecar (one row per "
                        "CPA-calling read: read_id, chrom, cpa_pos, gene_strand, "
                        "oaNT_tail_len, at_cpa, tier, mapq). On by default; needs pyarrow.")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--bbduk-bin", default="bbduk.sh", dest="bbduk_bin")
    p.add_argument("--bbmap-bin", default="bbmap.sh", dest="bbmap_bin")
    p.add_argument("--bbmap-mem", default=None, dest="bbmap_mem",
                   help="Cap the bbmap/bbduk JVM heap (e.g. 8g). Forwards -Xmx. "
                        "Set this on login / cgroup-limited nodes where bbmap's "
                        "auto-detected heap overshoots the ulimit; omit on "
                        "dedicated compute nodes (auto-detect).")
    # Tier-2 (extension; bowtie2 from the netseq_align env)
    p.add_argument("--tier2", action="store_true",
                   help="Enable Tier-2 CPA-window mini-ref rescue (CIRCULAR — poly-A length only).")
    p.add_argument("--tier2-bodyward", type=int, default=50, dest="tier2_bodyward")
    p.add_argument("--tier2-downstream", type=int, default=10, dest="tier2_downstream")
    p.add_argument("--tier2-min-anchor", type=int, default=12, dest="tier2_min_anchor")
    p.add_argument("--bowtie2-bin", default="bowtie2", dest="bowtie2_bin")
    p.add_argument("--bowtie2-build-bin", default="bowtie2-build", dest="bowtie2_build_bin")
