#!/usr/bin/env python3
"""Component D — SCORER + RECOVERY-vs-FALSE-JUNCTION-FDR TRADE CURVE.

Consumes the per-arm refined BAMs (SPEC item 6: ``arm_{A,B,C}.bam``), the per-read
truth join (SPEC item 5: ``read_truth.tsv``), and the genomic reference the reads
were aligned to (SPEC item 2: ``sim_ref.fa``). Optionally consumes the template
panel (SPEC item 3: ``panel_truth.tsv``) to recover the ``decoy_offset`` axis that
``read_truth.tsv`` does not carry. Emits SPEC item 7: ``trade_curve.json`` — per
(motif_rung x context [x decoy_offset]) cell, RECOVERY + FALSE-JUNCTION-FDR per
arm, plus overall + control FP rates — and prints a compact human summary.

WHY THIS IS THE DELIVERABLE (SPEC METRICS / advisor framing): recovery ALONE is
near-tautological (arm-B simply removed the motif prior, so of course it "recovers"
more non-YAG acceptors). The scientific claim is the TRADE: does recovery rise
WITHOUT a proportional rise in false non-canonical junctions on the controls? So
the lead metric is the stratified recovery-vs-FDR curve, not recovery.

MATCHING LOGIC (reused primitives; NOT reinvented):
  * ``scorer.extract_junctions`` walks each read's CIGAR -> N-op introns
    ``[start, end)`` in genome coords.
  * ``chimeric_consensus.normalize_junction`` left-normalizes an intron to its
    leftmost ambiguity-equivalent coordinate. Two junctions with the SAME spliced
    product normalize to the same tuple; a decoy snap changes intron *length*
    (donor fixed, acceptor moved), and normalize preserves length, so a snap can
    never collapse into the truth. This is the ambiguity-aware match.
  * RECOVERY(read) = the read has >=1 called N-op whose normalized (start,end)
    equals the read's normalized truth intron (on the truth contig). Because the
    donor is a fixed canonical GT by construction, full-intron equality is
    identical to "matches the true acceptor" — and correctly excludes decoy snaps.
  * FALSE-JUNCTION(read) = the read emits >=1 called N-op that is (a) NOT the truth
    junction AND (b) NON-CANONICAL, i.e. no ambiguity-equivalent placement is a
    canonical GT/GC..AG (strand-aware, via truth_schema._canonical_strand_aware).
    This is the "emits / moves-to an untrue non-canonical junction" of SPEC METRICS
    — the failure arm-B/C must NOT inflate. For INTRONFREE we ALSO report the
    literal SPEC-control rate "any junction emitted = FP".

Author: Kevin R. Roy (component D)
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# --- make the in-tree rectify package importable regardless of install state ---
# scripts/benchmark/noncanon_sim/score_trade.py -> parents[3] is the worktree root
_WORKTREE_ROOT = Path(__file__).resolve().parents[3]
if str(_WORKTREE_ROOT) not in sys.path:
    sys.path.insert(0, str(_WORKTREE_ROOT))

import pysam  # noqa: E402

# REUSE the load-bearing benchmark + ambiguity primitives (do not reinvent).
from rectify.core.benchmark.scorer import (  # noqa: E402
    extract_junctions,
    load_genome,
    cigar_records_to_bam,
)
from rectify.core.consensus.chimeric_consensus import (  # noqa: E402
    normalize_junction,
    junction_ambiguity_window,
)
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    _canonical_strand_aware,
)

RUNGS_NONCANONICAL = ("R1", "R2", "R3")
_NODECOY_TOKENS = {"", "0", "na", "none", ".", "nan", "null"}


# ---------------------------------------------------------------------------
# Truth ingestion (SPEC item 5 read_truth.tsv + optional item 3 panel_truth.tsv)
# ---------------------------------------------------------------------------
class ReadTruthRow:
    """A single row of ``read_truth.tsv`` (SPEC item 5), plus a joined
    ``decoy_offset`` from ``panel_truth.tsv`` (SPEC item 3) when available."""

    __slots__ = ("read_id", "tid", "chrom", "true_donor", "true_acceptor",
                 "strand", "motif_rung", "context", "has_true_junction",
                 "decoy_offset")

    def __init__(self, read_id, tid, chrom, true_donor, true_acceptor, strand,
                 motif_rung, context, has_true_junction, decoy_offset):
        self.read_id = read_id
        self.tid = tid
        self.chrom = chrom
        self.true_donor = true_donor          # int or None
        self.true_acceptor = true_acceptor    # int or None
        self.strand = strand
        self.motif_rung = motif_rung
        self.context = context
        self.has_true_junction = has_true_junction  # bool
        self.decoy_offset = decoy_offset      # str token or None (no panel_truth)

    @property
    def has_decoy(self) -> bool:
        if self.decoy_offset is None:
            return False
        return str(self.decoy_offset).strip().lower() not in _NODECOY_TOKENS


def _read_tsv(path: Path) -> Tuple[List[str], List[List[str]]]:
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    return header, rows


def _opt_int(tok: str) -> Optional[int]:
    tok = (tok or "").strip()
    if tok.lower() in _NODECOY_TOKENS:
        return None
    try:
        return int(tok)
    except ValueError:
        return None


def load_panel_decoys(panel_path: Path) -> Dict[str, str]:
    """tid -> decoy_offset token, from ``panel_truth.tsv`` (SPEC item 3)."""
    header, rows = _read_tsv(panel_path)
    idx = {c: i for i, c in enumerate(header)}
    if "tid" not in idx:
        raise ValueError(f"{panel_path}: panel_truth.tsv missing 'tid' column")
    if "decoy_offset" not in idx:
        # Panel exists but has no decoy axis; treat as no-decoy everywhere.
        return {}
    out: Dict[str, str] = {}
    for f in rows:
        out[f[idx["tid"]]] = f[idx["decoy_offset"]] if idx["decoy_offset"] < len(f) else ""
    return out


def load_read_truth(read_truth_path: Path,
                    panel_decoys: Optional[Dict[str, str]]) -> Dict[str, ReadTruthRow]:
    """Load ``read_truth.tsv`` (SPEC item 5) keyed by read_id, joining the
    ``decoy_offset`` from ``panel_truth.tsv`` on ``tid`` when provided."""
    header, rows = _read_tsv(read_truth_path)
    idx = {c: i for i, c in enumerate(header)}
    required = ["read_id", "tid", "chrom", "true_donor", "true_acceptor",
                "strand", "motif_rung", "context", "has_true_junction"]
    missing = [c for c in required if c not in idx]
    if missing:
        raise ValueError(
            f"{read_truth_path}: read_truth.tsv missing required column(s) "
            f"{missing}; SPEC item 5 header is: {' '.join(required)}")

    def get(f, col):
        i = idx[col]
        return f[i] if i < len(f) else ""

    out: Dict[str, ReadTruthRow] = {}
    for f in rows:
        rid = get(f, "read_id")
        tid = get(f, "tid")
        htj_tok = get(f, "has_true_junction").strip().lower()
        has_true = htj_tok in ("1", "true", "yes", "t")
        decoy = None
        if panel_decoys is not None:
            decoy = panel_decoys.get(tid, "")
        out[rid] = ReadTruthRow(
            read_id=rid,
            tid=tid,
            chrom=get(f, "chrom"),
            true_donor=_opt_int(get(f, "true_donor")),
            true_acceptor=_opt_int(get(f, "true_acceptor")),
            strand=get(f, "strand") or "+",
            motif_rung=get(f, "motif_rung"),
            context=get(f, "context"),
            has_true_junction=has_true,
            decoy_offset=decoy,
        )
    return out


# ---------------------------------------------------------------------------
# Per-read evaluation (the reused matching logic)
# ---------------------------------------------------------------------------
class ReadEval:
    __slots__ = ("placed", "recovered", "emits_any_junction",
                 "emits_untrue_junction", "emits_untrue_noncanonical")

    def __init__(self, placed=False, recovered=False, emits_any_junction=False,
                 emits_untrue_junction=False, emits_untrue_noncanonical=False):
        self.placed = placed
        self.recovered = recovered
        self.emits_any_junction = emits_any_junction
        self.emits_untrue_junction = emits_untrue_junction
        self.emits_untrue_noncanonical = emits_untrue_noncanonical


def _norm_intron(start: int, end: int, seq: str,
                 max_shift: int) -> Tuple[int, int]:
    """min/max-order the intron span (advisor: free insurance against a sibling
    storing donor/acceptor in biological order on minus-strand genes), then
    left-normalize on the genome seq. A no-op when already ascending."""
    s, e = (start, end) if start <= end else (end, start)
    if seq:
        return normalize_junction(s, e, seq, max_shift=max_shift)
    return s, e


def evaluate_read(read, rt: ReadTruthRow, genome: Dict[str, str],
                  max_shift: int) -> ReadEval:
    """Evaluate one primary alignment against its truth row. Reuses
    ``extract_junctions`` + ``normalize_junction`` + strand-aware canonicity."""
    mapped_seq = genome.get(read.reference_name, "")
    called = extract_junctions(read.reference_start, read.cigartuples)
    norm_called = [_norm_intron(cs, ce, mapped_seq, max_shift) for (cs, ce) in called]

    # Normalized truth intron (on the truth contig), if this read has one.
    truth_norm: Optional[Tuple[int, int]] = None
    if rt.has_true_junction and rt.true_donor is not None and rt.true_acceptor is not None:
        tseq = genome.get(rt.chrom, "")
        truth_norm = _norm_intron(rt.true_donor, rt.true_acceptor, tseq, max_shift)

    # Recovery: only meaningful when the read is on its truth contig.
    on_truth_contig = (read.reference_name == rt.chrom)
    truth_set = {truth_norm} if (truth_norm is not None and on_truth_contig) else set()

    recovered = truth_norm is not None and truth_norm in set(norm_called) and on_truth_contig

    emits_any = len(norm_called) > 0
    emits_untrue = False
    emits_untrue_noncanon = False
    for (ns, ne) in norm_called:
        if (ns, ne) in truth_set:
            continue  # this is the true junction — never an FP
        emits_untrue = True
        if mapped_seq:
            l_amb, r_amb = junction_ambiguity_window(ns, ne, mapped_seq)
            canonical = _canonical_strand_aware(ns, ne, mapped_seq, l_amb, r_amb,
                                                rt.strand)
        else:
            canonical = False  # no seq -> cannot claim canonical; conservative
        if not canonical:
            emits_untrue_noncanon = True

    return ReadEval(
        placed=True,
        recovered=recovered,
        emits_any_junction=emits_any,
        emits_untrue_junction=emits_untrue,
        emits_untrue_noncanonical=emits_untrue_noncanon,
    )


def evaluate_arm(bam_path: Path, truth: Dict[str, ReadTruthRow],
                 genome: Dict[str, str], max_shift: int
                 ) -> Tuple[Dict[str, ReadEval], int]:
    """Return {read_id: ReadEval} for every truth read (unplaced reads default to
    placed=False), plus a count of reads mapped to a contig absent from genome."""
    evals: Dict[str, ReadEval] = {rid: ReadEval() for rid in truth}
    missing_contig = 0
    seen = set()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            rt = truth.get(read.query_name)
            if rt is None:
                continue
            if read.query_name in seen:
                continue  # only the first primary
            seen.add(read.query_name)
            if read.reference_name not in genome:
                missing_contig += 1
            evals[read.query_name] = evaluate_read(read, rt, genome, max_shift)
    return evals, missing_contig


# ---------------------------------------------------------------------------
# Aggregation into cells + controls + overall
# ---------------------------------------------------------------------------
def _rate(num: int, den: int) -> Optional[float]:
    return round(num / den, 4) if den else None


def build_trade_curve(truth: Dict[str, ReadTruthRow],
                      arm_evals: Dict[str, Dict[str, ReadEval]],
                      has_decoy_axis: bool) -> Tuple[dict, list]:
    """Assemble the SPEC item 7 trade-curve dict. Returns (json_obj, cell_keys)."""
    arms = sorted(arm_evals)

    # cell key -> {rung, context, decoy} ; membership lists of read_ids
    cell_reads: Dict[Tuple[str, str, Optional[str]], List[str]] = defaultdict(list)
    intronfree_reads: List[str] = []
    r0_reads: List[str] = []
    yag_nodecoy_reads: List[str] = []
    noncanon_reads: List[str] = []

    for rid, rt in truth.items():
        if not rt.has_true_junction:
            intronfree_reads.append(rid)
            continue
        decoy = (str(rt.decoy_offset) if (has_decoy_axis and rt.has_decoy) else
                 ("none" if has_decoy_axis else None))
        cell_reads[(rt.motif_rung, rt.context, decoy)].append(rid)
        if rt.motif_rung == "R0":
            r0_reads.append(rid)
            if not rt.has_decoy:
                yag_nodecoy_reads.append(rid)
        if rt.motif_rung in RUNGS_NONCANONICAL:
            noncanon_reads.append(rid)

    def arm_cell_stats(read_ids: List[str], arm: str,
                       want_recovery: bool) -> dict:
        ev = arm_evals[arm]
        n = len(read_ids)
        placed = sum(ev[r].placed for r in read_ids)
        rec = sum(ev[r].recovered for r in read_ids)
        fj = sum(ev[r].emits_untrue_noncanonical for r in read_ids)
        anyj = sum(ev[r].emits_any_junction for r in read_ids)
        d = {
            "n_reads": n,
            "n_placed": placed,
            "n_false_junction": fj,
            "false_junction_FDR": _rate(fj, n),
        }
        if want_recovery:
            d["n_recovered"] = rec
            d["recovery"] = _rate(rec, n)
        d["n_any_junction"] = anyj
        d["any_junction_rate"] = _rate(anyj, n)
        return d

    # ---- per-cell rows (SPEC item 7 core) --------------------------------
    cells = []
    for key in sorted(cell_reads,
                      key=lambda k: (k[0], k[1], str(k[2]))):
        rung, context, decoy = key
        rids = cell_reads[key]
        row = {
            "motif_rung": rung,
            "context": context,
            "n_reads": len(rids),
            "arms": {a: arm_cell_stats(rids, a, want_recovery=True) for a in arms},
        }
        if has_decoy_axis:
            row["decoy_offset"] = decoy
        cells.append(row)

    # ---- controls --------------------------------------------------------
    controls = {
        "INTRONFREE": {
            "definition": "has_true_junction=0; any emitted junction is a FP",
            "n_reads": len(intronfree_reads),
            "arms": {a: arm_cell_stats(intronfree_reads, a, want_recovery=False)
                     for a in arms},
        },
        "R0_CANONICAL": {
            "definition": "motif_rung=R0 (canonical YAG); FP=moved-to untrue non-canonical",
            "n_reads": len(r0_reads),
            "arms": {a: arm_cell_stats(r0_reads, a, want_recovery=True)
                     for a in arms},
        },
        "YAG_CANONICAL_NO_DECOY": {
            "definition": ("motif_rung=R0 with NO nearby decoy; FP=moved-to untrue "
                           "non-canonical" if has_decoy_axis else
                           "no panel_truth decoy axis -> identical to R0_CANONICAL"),
            "n_reads": len(yag_nodecoy_reads),
            "arms": {a: arm_cell_stats(yag_nodecoy_reads, a, want_recovery=True)
                     for a in arms},
        },
    }

    # ---- overall (the LEAD trade: recovery vs controls FDR) --------------
    control_pool = sorted(set(intronfree_reads) | set(r0_reads))
    overall = {"arms": {}}
    by_rung_reads = defaultdict(list)
    for rid, rt in truth.items():
        if rt.has_true_junction:
            by_rung_reads[rt.motif_rung].append(rid)
    for a in arms:
        ev = arm_evals[a]
        n_nc = len(noncanon_reads)
        rec_nc = sum(ev[r].recovered for r in noncanon_reads)
        n_ctrl = len(control_pool)
        fj_ctrl = sum(ev[r].emits_untrue_noncanonical for r in control_pool)
        overall["arms"][a] = {
            "recovery_noncanonical": _rate(rec_nc, n_nc),
            "n_recovered_noncanonical": rec_nc,
            "n_reads_noncanonical": n_nc,
            "false_junction_FDR_controls": _rate(fj_ctrl, n_ctrl),
            "n_false_junction_controls": fj_ctrl,
            "n_reads_controls": n_ctrl,
            "recovery_by_rung": {
                rung: _rate(sum(ev[r].recovered for r in rids), len(rids))
                for rung, rids in sorted(by_rung_reads.items())
            },
        }

    obj = {
        "cells": cells,
        "controls": controls,
        "overall": overall,
    }
    return obj, sorted(cell_reads)


# ---------------------------------------------------------------------------
# Human-readable summary
# ---------------------------------------------------------------------------
def print_summary(obj: dict, arms: List[str]) -> None:
    has_decoy = any("decoy_offset" in c for c in obj["cells"])
    print("\n=== TRADE CURVE: recovery vs false-junction FDR (per cell) ===")
    hdr = ["rung", "context"]
    if has_decoy:
        hdr.append("decoy")
    hdr.append("n")
    for a in arms:
        hdr += [f"{a}.rec", f"{a}.fjFDR"]
    widths = [6, 10] + ([6] if has_decoy else []) + [5] + [8, 8] * len(arms)
    line = "  ".join(h.ljust(w) for h, w in zip(hdr, widths))
    print(line)
    print("-" * len(line))
    for c in obj["cells"]:
        cols = [str(c["motif_rung"]), str(c["context"])]
        if has_decoy:
            cols.append(str(c.get("decoy_offset")))
        cols.append(str(c["n_reads"]))
        for a in arms:
            s = c["arms"][a]
            cols.append(_fmt(s.get("recovery")))
            cols.append(_fmt(s.get("false_junction_FDR")))
        print("  ".join(col.ljust(w) for col, w in zip(cols, widths)))

    print("\n=== CONTROLS (false-junction FP rates) ===")
    for name, blk in obj["controls"].items():
        cells = [name.ljust(24), f"n={blk['n_reads']}".ljust(8)]
        for a in arms:
            s = blk["arms"][a]
            cells.append(f"{a}: fjFDR={_fmt(s.get('false_junction_FDR'))} "
                         f"any={_fmt(s.get('any_junction_rate'))}")
        print("  ".join(cells))

    print("\n=== OVERALL (lead trade) ===")
    for a in arms:
        o = obj["overall"]["arms"][a]
        print(f"  arm {a}: recovery(R1-R3)={_fmt(o['recovery_noncanonical'])}  "
              f"false_junction_FDR(controls+R0)={_fmt(o['false_junction_FDR_controls'])}  "
              f"[rec {o['n_recovered_noncanonical']}/{o['n_reads_noncanonical']}, "
              f"fj {o['n_false_junction_controls']}/{o['n_reads_controls']}]")
    print()


def _fmt(x) -> str:
    return "  na" if x is None else f"{x:.3f}"


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def _resolve_inputs(args) -> dict:
    wd = Path(args.work_dir).resolve() if args.work_dir else Path.cwd()

    def pick(explicit, default_name):
        if explicit:
            return Path(explicit).resolve()
        return (wd / default_name)

    read_truth = pick(args.read_truth, "read_truth.tsv")
    sim_ref = pick(args.sim_ref, "sim_ref.fa")
    panel_truth = pick(args.panel_truth, "panel_truth.tsv") if args.panel_truth \
        else (wd / "panel_truth.tsv")

    arms: Dict[str, Path] = {}
    if args.arm:
        for spec in args.arm:
            if "=" not in spec:
                raise SystemExit(f"--arm expects NAME=PATH, got {spec!r}")
            name, path = spec.split("=", 1)
            arms[name] = Path(path).resolve()
    else:
        for name in ("A", "B", "C"):
            p = wd / f"arm_{name}.bam"
            if p.exists():
                arms[name] = p
    return {"read_truth": read_truth, "sim_ref": sim_ref,
            "panel_truth": panel_truth, "arms": arms,
            "out": Path(args.out).resolve() if args.out else (wd / "trade_curve.json")}


def run(args) -> int:
    inp = _resolve_inputs(args)

    def need(p: Path, what: str):
        if not p.exists():
            raise SystemExit(
                f"ERROR: required input {what} not found: {p}\n"
                f"       (run components A/B/C first, or pass an explicit path)")

    need(inp["read_truth"], "read_truth.tsv (SPEC item 5)")
    need(inp["sim_ref"], "sim_ref.fa (SPEC item 2)")
    if not inp["arms"]:
        raise SystemExit(
            "ERROR: no arm BAMs found. Pass --arm A=arm_A.bam [--arm B=...] "
            "or place arm_{A,B,C}.bam in --work-dir.")
    for name, bp in inp["arms"].items():
        need(bp, f"arm_{name}.bam (SPEC item 6)")

    panel_decoys = None
    has_decoy_axis = False
    if inp["panel_truth"].exists():
        panel_decoys = load_panel_decoys(inp["panel_truth"])
        has_decoy_axis = True

    truth = load_read_truth(inp["read_truth"], panel_decoys)
    if not truth:
        raise SystemExit(f"ERROR: read_truth.tsv had no rows: {inp['read_truth']}")

    genome = load_genome(inp["sim_ref"])

    arm_evals: Dict[str, Dict[str, ReadEval]] = {}
    missing_by_arm: Dict[str, int] = {}
    for name, bp in inp["arms"].items():
        ev, miss = evaluate_arm(bp, truth, genome, args.normalize_window)
        arm_evals[name] = ev
        missing_by_arm[name] = miss
        if miss:
            print(f"WARNING: arm {name}: {miss} read(s) mapped to a contig ABSENT "
                  f"from {inp['sim_ref'].name} — their junctions scored without the "
                  f"ambiguity-aware normalize. Check sim_ref covers all BAM contigs.",
                  file=sys.stderr)

    obj, _ = build_trade_curve(truth, arm_evals, has_decoy_axis)

    obj["meta"] = {
        "generated": datetime.date.today().isoformat(),
        "script": "scripts/benchmark/noncanon_sim/score_trade.py",
        "component": "D (scorer + trade curve)",
        "sim_ref": str(inp["sim_ref"]),
        "read_truth": str(inp["read_truth"]),
        "panel_truth": str(inp["panel_truth"]) if has_decoy_axis else None,
        "has_decoy_axis": has_decoy_axis,
        "arms": {n: str(p) for n, p in inp["arms"].items()},
        "reads_mapped_to_missing_contig": missing_by_arm,
        "normalize_window_bp": args.normalize_window,
        "n_reads_total": len(truth),
        "n_reads_intronfree": sum(1 for r in truth.values() if not r.has_true_junction),
        "recovery_match": ("read has >=1 called N-op whose min/max-ordered + "
                           "left-normalized (start,end) equals the read's normalized "
                           "truth intron on its truth contig; full-intron equality == "
                           "true-acceptor match because the donor is a fixed canonical "
                           "GT (a decoy snap changes intron length, which normalize "
                           "preserves, so a snap never matches truth)"),
        "false_junction_match": ("read emits >=1 called N-op that is NOT the truth "
                                 "junction AND non-canonical (no ambiguity-equivalent "
                                 "placement is GT/GC..AG, strand-aware via "
                                 "truth_schema._canonical_strand_aware); INTRONFREE "
                                 "also reports the literal 'any junction emitted = FP' "
                                 "rate as any_junction_rate"),
        "reused_primitives": [
            "scorer.extract_junctions", "scorer.load_genome",
            "chimeric_consensus.normalize_junction",
            "chimeric_consensus.junction_ambiguity_window",
            "truth_schema._canonical_strand_aware",
        ],
    }

    inp["out"].write_text(json.dumps(obj, indent=2) + "\n")
    print_summary(obj, sorted(inp["arms"]))
    print(f"wrote {inp['out']}")
    return 0


# ---------------------------------------------------------------------------
# Self-test (hand-computable fixture — no minimap2, no siblings needed)
# ---------------------------------------------------------------------------
def self_test() -> int:
    """Build a 3-read synthetic fixture with cigar_records_to_bam and assert the
    recovery / false-junction numbers by construction. Confirms in-frame with
    extract_junctions' intron_end convention (acceptor dinuc = genome[acc-2:acc])."""
    import tempfile

    # contig (0-based):
    #   0-9   exon1 "ACGTACGTAA"
    #   10-11 donor "GT"
    #   12-27 intron fill "C"*16
    #   28-29 true acceptor "AC"  (NON-canonical -> R3)
    #   30    exon2 "T"
    #   31-32 decoy acceptor "AG" (canonical; snap target = intron [10,33))
    #   33-59 exon2 fill (27 bp)
    genome_seq = ("ACGTACGTAA" + "GT" + "C" * 16 + "AC" + "T" + "AG"
                  + "GATCGATCGATCGATCGATCGATCGAT")
    assert len(genome_seq) == 60, len(genome_seq)
    assert genome_seq[10:12] == "GT"
    assert genome_seq[28:30] == "AC"      # true acceptor (non-canonical)
    assert genome_seq[31:33] == "AG"      # decoy acceptor (canonical)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        ref = td / "sim_ref.fa"
        ref.write_text(">chrSIM_0\n" + genome_seq + "\n")

        # read_truth.tsv (SPEC item 5)
        rt = td / "read_truth.tsv"
        rt.write_text(
            "read_id\ttid\tchrom\ttrue_donor\ttrue_acceptor\tstrand\tmotif_rung\tcontext\thas_true_junction\n"
            "tid_R3_r000\ttid_R3\tchrSIM_0\t10\t30\t+\tR3\tplain\t1\n"   # will be recovered
            "tid_R3_r001\ttid_R3\tchrSIM_0\t10\t30\t+\tR3\tplain\t1\n"   # snapped to decoy
            "tid_IF_r000\ttid_IF\tchrSIM_0\t\t\t+\tINTRONFREE\tplain\t0\n"  # INTRONFREE
        )
        # panel_truth.tsv (SPEC item 3) — exercise the decoy_offset join
        pt = td / "panel_truth.tsv"
        pt.write_text(
            "tid\tchrom\ttrue_donor\ttrue_acceptor\tstrand\tintron_len\tmotif_rung\t"
            "acceptor_motif\tdecoy_offset\tdecoy_acceptor\texon5_len\texon3_len\t"
            "context\thas_true_junction\n"
            "tid_R3\tchrSIM_0\t10\t30\t+\t20\tR3\tAC\t3\t33\t10\t30\tplain\t1\n"
            "tid_IF\tchrSIM_0\t\t\t+\t0\tINTRONFREE\t\t\t\t60\t0\tplain\t0\n"
        )

        # Fabricated arm BAM via reused cigar_records_to_bam:
        #  read a: true junction  -> 10M 20N 30M  (N-op [10,30))  => recovered
        #  read b: decoy snap      -> 10M 23N 27M  (N-op [10,33))  => NOT recovered,
        #          canonical (GT..AG) so NOT a non-canonical FP
        #  read c: INTRONFREE with fabricated non-canon N-op -> 12M 13N 35M (N-op [12,25))
        #          acceptor genome[23:25]="CC" -> non-canonical => any + non-canon FP
        records = [
            ("tid_R3_r000", "chrSIM_0", 0, [(0, 10), (3, 20), (0, 30)], "A" * 40),
            ("tid_R3_r001", "chrSIM_0", 0, [(0, 10), (3, 23), (0, 27)], "A" * 37),
            ("tid_IF_r000", "chrSIM_0", 0, [(0, 12), (3, 13), (0, 35)], "A" * 47),
        ]
        bam = td / "arm_A.bam"
        cigar_records_to_bam(records, str(ref), str(bam))

        class A:
            work_dir = str(td)
            read_truth = None
            sim_ref = None
            panel_truth = str(pt)
            arm = None
            out = str(td / "trade_curve.json")
            normalize_window = 15

        rc = run(A)
        obj = json.loads((td / "trade_curve.json").read_text())

    # ---- assertions -----------------------------------------------------
    r3 = [c for c in obj["cells"] if c["motif_rung"] == "R3"][0]
    assert r3["n_reads"] == 2, r3
    assert r3.get("decoy_offset") == "3", r3        # panel join worked
    a = r3["arms"]["A"]
    assert a["recovery"] == 0.5, a                  # 1 recovered, 1 snapped
    assert a["n_recovered"] == 1, a
    assert a["false_junction_FDR"] == 0.0, a        # decoy snap is CANONICAL -> no FP

    intronfree = obj["controls"]["INTRONFREE"]["arms"]["A"]
    assert intronfree["false_junction_FDR"] == 1.0, intronfree   # non-canon fabricated
    assert intronfree["any_junction_rate"] == 1.0, intronfree    # literal any-junction FP

    ov = obj["overall"]["arms"]["A"]
    assert ov["recovery_noncanonical"] == 0.5, ov
    # controls pool = INTRONFREE(1) + R0(0) = 1 read, 1 non-canon FP
    assert ov["false_junction_FDR_controls"] == 1.0, ov

    # ---- minus-strand canonicity path (feeds the lead FDR metric) --------
    # A canonical MINUS-strand intron has genome motif CT..AC (revcomp of GT..AG).
    # Raw _canonical_within_window on the forward genome would call that
    # non-canonical and wrongly charge a minus-strand canonical snap as a
    # false-junction FP — inflating arm-B/C FDR on half the panel. Verify the
    # strand-aware wrapper prevents that.
    #   0-9  "ACGTACGTAA" | 10-11 "GT" | 12-13 "CT" | 14-22 "G"*9 | 23-24 "AC" | 25-59 fill
    mg = ("ACGTACGTAA" + "GT" + "CT" + "G" * 9 + "AC"
          + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
    assert len(mg) == 60, len(mg)
    assert mg[12:14] == "CT" and mg[23:25] == "AC"   # CT..AC = minus-canonical
    assert mg[14:16] == "GG" and mg[21:23] == "GG"   # GG..GG = non-canonical both strands
    mgenome = {"chrSIM_0": mg}

    class _StubRead:
        is_secondary = is_supplementary = is_unmapped = is_reverse = False

        def __init__(self, cigartuples):
            self.reference_name = "chrSIM_0"
            self.reference_start = 0
            self.cigartuples = cigartuples
            self.mapping_quality = 60

    # minus-strand read, TRUE junction at [10,30), called an UNTRUE CT..AC N-op [12,25)
    rt_minus = ReadTruthRow("m0", "tid_m", "chrSIM_0", 10, 30, "-", "R3", "plain",
                            True, None)
    ev = evaluate_read(_StubRead([(0, 12), (3, 13), (0, 35)]), rt_minus, mgenome, 15)
    assert ev.emits_untrue_junction is True, ev.__dict__
    assert ev.emits_untrue_noncanonical is False, (
        "minus-strand CT..AC snap must NOT count as a non-canonical FP "
        f"(strand-aware canonicity broken): {ev.__dict__}")

    # minus-strand read, called an UNTRUE GG..GG N-op [14,23): non-canonical BOTH strands
    ev2 = evaluate_read(_StubRead([(0, 14), (3, 9), (0, 37)]), rt_minus, mgenome, 15)
    assert ev2.emits_untrue_noncanonical is True, (
        f"GG..GG must count as a non-canonical FP on either strand: {ev2.__dict__}")

    print("SELF-TEST PASSED: recovery/false-junction matching verified on the "
          "hand-computed fixture (in-frame with extract_junctions intron_end); "
          "minus-strand canonicity (CT..AC canonical, GG..GG not) confirmed.")
    return rc


# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Component D: score arm BAMs into a recovery-vs-false-junction-"
                    "FDR trade curve (SPEC item 7).")
    p.add_argument("--work-dir", help="dir holding read_truth.tsv / sim_ref.fa / "
                   "panel_truth.tsv / arm_{A,B,C}.bam (default: cwd)")
    p.add_argument("--read-truth", help="override path to read_truth.tsv (SPEC item 5)")
    p.add_argument("--sim-ref", help="override path to sim_ref.fa (SPEC item 2)")
    p.add_argument("--panel-truth", help="override path to panel_truth.tsv (SPEC item 3; "
                   "supplies the decoy_offset axis)")
    p.add_argument("--arm", action="append", metavar="NAME=PATH",
                   help="arm BAM as NAME=PATH (repeatable). Default: auto-discover "
                        "arm_{A,B,C}.bam in --work-dir.")
    p.add_argument("--out", help="output trade_curve.json (default: <work-dir>/trade_curve.json)")
    p.add_argument("--normalize-window", type=int, default=15,
                   help="max junction ambiguity slide in bp (default 15, matches "
                        "chimeric_consensus._JUNCTION_AMBIGUITY_MAX_SHIFT)")
    p.add_argument("--self-test", action="store_true",
                   help="run the built-in synthetic fixture + assertions and exit")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test:
        return self_test()
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
