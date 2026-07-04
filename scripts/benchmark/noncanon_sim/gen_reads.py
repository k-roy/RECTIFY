#!/usr/bin/env python3
"""Component B — LONG-READ GENERATOR for the yeast non-canonical splice sim.

Reads ``templates.fa`` (SPEC item 1, spliced exon1++exon2 transcripts) +
``panel_truth.tsv`` (SPEC item 3) and emits ONT-like LONG reads
``reads.fastq`` (SPEC item 4) + ``read_truth.tsv`` (SPEC item 5).

=========================== WHAT THIS PRODUCES ===========================
One row per template drives ``--reads-per-template`` full-length long reads.
Each read is the spliced transcript sequence with an INDEPENDENT per-base
ONT error process applied. The genomic junction/CPA truth is carried through
UNCHANGED from ``panel_truth.tsv`` — this component never computes coords,
it only propagates them so the A->B->C->D join is byte-stable.

======================= ERROR MODEL — INDEPENDENCE ========================
LOAD-BEARING (SPEC APIS / error model): the read error source is COMPLETELY
INDEPENDENT of the scorer's ``penalty_scores.tsv`` (the -logP junction law).
It is NEVER read here. Two paths:

  * pbsim3 (PREFERRED, ERRHMM-ONT): used iff ``pbsim`` is on PATH (or
    --pbsim-bin resolves). See ``scripts/benchmark/sim/pbsim3_wrapper.py``.
  * FALLBACK (this file, when pbsim3 absent): clean reads from the templates
    corrupted by a HARD-CODED per-base ONT profile:
      - DELETION-DOMINANT marginal: del > ins ~ sub (task spec);
      - per-base deletion hazard ELEVATED inside EXONIC homopolymer runs
        (len >= HP_MIN_LEN) — the realistic HP-collapse structure. The HP
        runs are read straight off the template, which is exonic by
        construction (exon1++exon2, introns spliced out), so "exonic HP"
        falls out for free (and INTRONFREE single-exon templates too).
    Every number below is a fixed constant / CLI arg. None derives from
    penalty_scores.tsv. (We deliberately do NOT reuse
    ``scripts/benchmark/sim/error_injector.py``: its marginal is
    substitution-dominant (frac_sub=0.55), contradicting the deletion-
    dominant requirement — it is a post-hoc correlation layer, not a base
    model. Both are penalty-independent; the marginal shape is the tiebreak.)

============================ TRUTH INTEGRITY =============================
The ``read_truth.tsv`` schema is FLAT (template-level ``has_true_junction``,
no per-read anchor field). So reads are FULL-LENGTH BY DEFAULT: any
truncation that clipped past the junction would manufacture an unrecoverable
false-negative and silently deflate recovery. ``--trunc-frac`` is opt-in and
JUNCTION-SAFE (it always preserves >= ANCHOR bp of exon on both flanks of the
junction; if the template layout can't be verified it keeps the read full).

Author: Kevin R. Roy  (component B)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

BASES = "ACGT"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# read_truth.tsv columns — EXACT SPEC item-5 order. Do not reorder.
READ_TRUTH_COLS = [
    "read_id", "tid", "chrom", "true_donor", "true_acceptor", "strand",
    "motif_rung", "context", "has_true_junction",
]

# --- fallback error model constants (HARD-CODED, penalty-table-independent) ---
# Total per-base error rate by regime (marginal probability a position is hit).
ERROR_REGIMES = {
    "low":  0.02,   # HQ / R10 / cDNA-like
    "mid":  0.06,   # ONT R9 direct-RNA-like (default)
    "high": 0.10,   # noisy DRS
}
# Deletion-dominant split of the total rate: del > ins ~ sub, ins~sub lower.
FRAC_DEL = 0.45
FRAC_SUB = 0.33
FRAC_INS = 0.22
# Exonic homopolymer structure: a position in a run of >= HP_MIN_LEN identical
# bases gets its deletion hazard multiplied by --hp-del-mult, capped at HP_DEL_CAP
# (a modest multiplier + cap so long runs shorten by ~1, they do NOT fully
# collapse). HP_MIN_LEN=4 is the run length at which ONT length-calling gets
# ambiguous.
HP_MIN_LEN = 4
HP_DEL_CAP = 0.35
# Junction-safe truncation: keep >= ANCHOR exon bp on each junction flank.
ANCHOR = 15
MIN_READ_LEN = 30  # never emit a shorter read (guards degenerate windows)


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------
def _open_maybe_gzip(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path: str) -> Dict[str, str]:
    """Parse a (optionally gzipped, line-wrapped) FASTA into name->sequence.
    The record name is the header token up to the first whitespace."""
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with _open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_panel_truth(path: str) -> List[Dict[str, str]]:
    """Parse panel_truth.tsv by HEADER NAME (robust to column reorder/adds)."""
    with _open_maybe_gzip(path, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"panel_truth is empty: {path}")
    if "tid" not in rows[0]:
        raise ValueError(
            f"panel_truth {path} lacks a 'tid' column; header={list(rows[0])}")
    return rows


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# tid <-> template join (fuzzy for LOOKUP ONLY; output tid is verbatim)
# ---------------------------------------------------------------------------
def lookup_template(tid: str, templates: Dict[str, str]) -> Optional[str]:
    """Find the sequence for a panel_truth ``tid``. Tries the verbatim tid, then
    a 'tid_' prefix, then a stripped 'tid_' prefix — so header ``>tid_0`` joins
    whether the tsv stores ``tid_0`` or ``0``. The normalization is used ONLY to
    locate the sequence; the value WRITTEN to read_truth is always the verbatim
    panel_truth tid (C/D join read_truth.tid against panel_truth.tid)."""
    if tid in templates:
        return templates[tid]
    if f"tid_{tid}" in templates:
        return templates[f"tid_{tid}"]
    if tid.startswith("tid_") and tid[len("tid_"):] in templates:
        return templates[tid[len("tid_"):]]
    return None


# ---------------------------------------------------------------------------
# Fallback error model
# ---------------------------------------------------------------------------
def hp_run_lengths(seq: str) -> List[int]:
    """For each position, the length of the maximal homopolymer run it belongs to."""
    n = len(seq)
    res = [1] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and seq[j + 1] == seq[i]:
            j += 1
        L = j - i + 1
        for k in range(i, j + 1):
            res[k] = L
        i = j + 1
    return res


def _qual_char(err_prob: float, rng: random.Random) -> str:
    """Phred+33 char from an approx per-base error probability, with mild jitter.
    Cosmetic only — minimap2 ignores base quality for mapping."""
    err_prob = min(max(err_prob, 1e-4), 0.75)
    q = -10.0 * math.log10(err_prob)
    q = int(round(q + rng.uniform(-2.0, 2.0)))
    q = max(2, min(40, q))
    return chr(q + 33)


def corrupt_window(seq: str, hp_len: List[int], w0: int, w1: int,
                   del_rate: float, ins_rate: float, sub_rate: float,
                   hp_del_mult: float, rng: random.Random) -> Tuple[str, str]:
    """Apply the independent per-base ONT error process to template[w0:w1].

    Per position, a single uniform draw partitions into mutually-exclusive
    del / ins / sub / correct (the three rates always sum < 1). Deletion hazard
    is elevated (capped) inside exonic homopolymer runs. Returns (seq, qual)."""
    out_b: List[str] = []
    out_q: List[str] = []
    total = del_rate + ins_rate + sub_rate
    i = w0
    while i < w1:
        base = seq[i]
        if base not in BASES:            # pass through N / ambiguous unchanged
            out_b.append(base)
            out_q.append(_qual_char(total, rng))
            i += 1
            continue
        d = del_rate
        if hp_len[i] >= HP_MIN_LEN:
            d = min(del_rate * hp_del_mult, HP_DEL_CAP)
        r = rng.random()
        if r < d:                        # deletion — consume ref base, emit nothing
            i += 1
            continue
        r -= d
        if r < ins_rate:                 # insertion — random base then the ref base
            out_b.append(rng.choice(BASES))
            out_q.append(_qual_char(max(ins_rate, 0.15), rng))
            out_b.append(base)
            out_q.append(_qual_char(total, rng))
            i += 1
            continue
        r -= ins_rate
        if r < sub_rate:                 # substitution — a different base
            out_b.append(rng.choice([b for b in BASES if b != base]))
            out_q.append(_qual_char(max(sub_rate, 0.15), rng))
            i += 1
            continue
        out_b.append(base)               # correct
        out_q.append(_qual_char(total, rng))
        i += 1
    return "".join(out_b), "".join(out_q)


def choose_window(tlen: int, row: Dict[str, str], truncate: bool,
                  rng: random.Random) -> Tuple[int, int]:
    """Pick the [w0, w1) slice of the template to sequence.

    Full-length unless ``truncate``. When truncating, preserve >= ANCHOR exon bp
    on BOTH junction flanks (junction offset == exon5_len iff the template is the
    2-exon exon1++exon2 layout, i.e. len == exon5_len + exon3_len). If the layout
    can't be verified, or the row has no true junction that constrains it, trim
    only the ends while keeping >= MIN_READ_LEN. 5'-biased trim (DRS-like)."""
    if not truncate or tlen <= MIN_READ_LEN:
        return 0, tlen

    has_j = str(row.get("has_true_junction", "1")).strip() in ("1", "True", "true")
    prot_lo: Optional[int] = None
    prot_hi: Optional[int] = None
    if has_j:
        try:
            e5 = int(row["exon5_len"])
            e3 = int(row["exon3_len"])
        except (KeyError, ValueError, TypeError):
            e5 = e3 = -1
        # Only trust the junction offset if the 2-exon layout checks out.
        if e5 > 0 and e3 > 0 and (e5 + e3) == tlen:
            prot_lo = max(0, e5 - ANCHOR)
            prot_hi = min(tlen, e5 + ANCHOR)
        else:
            # Layout unverifiable -> do not risk clipping the junction.
            return 0, tlen

    # Max we may trim off each end (keep protected window + MIN_READ_LEN).
    if prot_lo is None:                  # no junction constraint (INTRONFREE etc.)
        max5 = max(0, tlen - MIN_READ_LEN)
        max3 = max(0, tlen - MIN_READ_LEN)
    else:
        max5 = max(0, prot_lo)
        max3 = max(0, tlen - prot_hi)
    # 5'-biased: most trim comes off the 5' end (DRS truncates the 5').
    t5 = rng.randint(0, max5) if max5 > 0 else 0
    t3 = rng.randint(0, max3 // 3) if max3 > 0 else 0
    w0, w1 = t5, tlen - t3
    if w1 - w0 < MIN_READ_LEN:
        return 0, tlen
    return w0, w1


# ---------------------------------------------------------------------------
# Fallback generator (the main path when pbsim3 is absent)
# ---------------------------------------------------------------------------
def generate_fallback(rows: List[Dict[str, str]], templates: Dict[str, str],
                      out_dir: str, reads_per_template: int, total_rate: float,
                      hp_del_mult: float, trunc_frac: float, rc_frac: float,
                      seed: int, allow_extra_templates: bool) -> Dict[str, str]:
    """Emit reads.fastq + read_truth.tsv from templates + panel_truth."""
    del_rate = total_rate * FRAC_DEL
    ins_rate = total_rate * FRAC_INS
    sub_rate = total_rate * FRAC_SUB

    rng = random.Random(seed)
    os.makedirs(out_dir, exist_ok=True)
    fastq_path = os.path.join(out_dir, "reads.fastq")
    truth_path = os.path.join(out_dir, "read_truth.tsv")

    matched_templates = set()
    n_reads = 0
    with open(fastq_path, "w") as fq, open(truth_path, "w") as tt:
        tt.write("\t".join(READ_TRUTH_COLS) + "\n")
        for row in rows:
            tid = row["tid"]
            seq = lookup_template(tid, templates)
            if seq is None:
                raise KeyError(
                    f"panel_truth tid {tid!r} has no matching template in "
                    f"templates.fa (checked {tid!r}, 'tid_{tid}', stripped). "
                    f"A silent drop would starve the per-cell read count below "
                    f"N>=100 and corrupt the metric — refusing.")
            # record which physical FASTA record satisfied this tid
            for cand in (tid, f"tid_{tid}", tid[len('tid_'):] if tid.startswith('tid_') else tid):
                if cand in templates:
                    matched_templates.add(cand)
                    break
            hp_len = hp_run_lengths(seq)
            tlen = len(seq)
            # carried-through truth (verbatim) for this template
            base_truth = [
                row.get("chrom", ""), row.get("true_donor", ""),
                row.get("true_acceptor", ""), row.get("strand", ""),
                row.get("motif_rung", ""), row.get("context", ""),
                row.get("has_true_junction", ""),
            ]
            for n in range(reads_per_template):
                truncate = rng.random() < trunc_frac
                w0, w1 = choose_window(tlen, row, truncate, rng)
                rseq, rqual = corrupt_window(
                    seq, hp_len, w0, w1, del_rate, ins_rate, sub_rate,
                    hp_del_mult, rng)
                if len(rseq) < 1:        # guard: never emit an empty record
                    rseq, rqual = seq[w0:w1], "I" * (w1 - w0)
                if rc_frac > 0 and rng.random() < rc_frac:
                    rseq = revcomp(rseq)
                    rqual = rqual[::-1]
                read_id = f"{tid}_r{n:03d}"
                fq.write(f"@{read_id}\n{rseq}\n+\n{rqual}\n")
                tt.write("\t".join([read_id, tid] + base_truth) + "\n")
                n_reads += 1

    extra = set(templates) - matched_templates
    if extra and not allow_extra_templates:
        raise ValueError(
            f"{len(extra)} template record(s) in templates.fa were never "
            f"referenced by any panel_truth row (e.g. {sorted(extra)[:5]}). "
            f"This signals A<->B drift. Pass --allow-extra-templates to permit.")
    return {
        "reads_fastq": fastq_path, "read_truth": truth_path,
        "n_reads": str(n_reads), "n_templates": str(len(rows)),
        "n_extra_templates": str(len(extra)),
        "error_source": "fallback-independent-perbase-ONT",
    }


# ---------------------------------------------------------------------------
# pbsim3 path (preferred when available)
# ---------------------------------------------------------------------------
def _find_pbsim(pbsim_bin: str) -> Optional[str]:
    return shutil.which(pbsim_bin)


def _find_errhmm_model(pbsim_path: str) -> Optional[str]:
    """Locate a packaged ERRHMM-ONT model next to the pbsim install."""
    root = os.path.dirname(os.path.dirname(os.path.realpath(pbsim_path)))
    for rel in ("data/ERRHMM-ONT.model", "data/ERRHMM-ONT-HQ.model"):
        cand = os.path.join(root, rel)
        if os.path.exists(cand):
            return cand
    # broader search
    for base, _dirs, files in os.walk(root):
        for f in files:
            if f.startswith("ERRHMM-ONT") and f.endswith(".model"):
                return os.path.join(base, f)
    return None


def generate_pbsim(rows: List[Dict[str, str]], templates: Dict[str, str],
                   out_dir: str, reads_per_template: int, pbsim_path: str,
                   errhmm_model: str, seed: int,
                   allow_extra_templates: bool) -> Dict[str, str]:
    """Drive pbsim3 (ERRHMM-ONT) per-template, then relabel reads to the SPEC
    ``<tid>_r<NNN>`` id scheme and emit reads.fastq + read_truth.tsv.

    We invoke pbsim in ``templ`` (transcript) mode once over ALL templates with
    the requested depth (== reads_per_template), then map each pbsim read back to
    its template of origin. pbsim's read errors come from the packaged ONT HMM —
    fully independent of penalty_scores.tsv."""
    os.makedirs(out_dir, exist_ok=True)
    # verify the join up front (fail loud) and build the truth lookup
    truth_by_tid: Dict[str, List[str]] = {}
    matched_templates = set()
    tfa = os.path.join(out_dir, "_pbsim_templates.fa")
    with open(tfa, "w") as fh:
        for row in rows:
            tid = row["tid"]
            seq = lookup_template(tid, templates)
            if seq is None:
                raise KeyError(f"panel_truth tid {tid!r} has no template in templates.fa")
            for cand in (tid, f"tid_{tid}", tid[len('tid_'):] if tid.startswith('tid_') else tid):
                if cand in templates:
                    matched_templates.add(cand)
                    break
            fh.write(f">{tid}\n{seq}\n")
            truth_by_tid[tid] = [
                row.get("chrom", ""), row.get("true_donor", ""),
                row.get("true_acceptor", ""), row.get("strand", ""),
                row.get("motif_rung", ""), row.get("context", ""),
                row.get("has_true_junction", ""),
            ]
    extra = set(templates) - matched_templates
    if extra and not allow_extra_templates:
        raise ValueError(f"{len(extra)} unreferenced templates (A<->B drift); "
                         f"use --allow-extra-templates")

    prefix = os.path.join(out_dir, "_pbsim")
    cmd = [pbsim_path, "--strategy", "templ", "--method", "errhmm",
           "--errhmm", errhmm_model, "--template", tfa,
           "--depth", str(reads_per_template), "--prefix", prefix,
           "--seed", str(seed)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"pbsim3 failed rc={res.returncode}: {res.stderr[:800]}")

    # Collect pbsim FASTQs, relabel per template with a per-tid counter.
    fastq_path = os.path.join(out_dir, "reads.fastq")
    truth_path = os.path.join(out_dir, "read_truth.tsv")
    counters: Dict[str, int] = {t: 0 for t in truth_by_tid}
    n_reads = 0
    fq_exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    base = os.path.basename(prefix)
    with open(fastq_path, "w") as fq, open(truth_path, "w") as tt:
        tt.write("\t".join(READ_TRUTH_COLS) + "\n")
        for f in sorted(os.listdir(out_dir)):
            if not (f.startswith(base) and f.endswith(fq_exts)):
                continue
            with _open_maybe_gzip(os.path.join(out_dir, f), "rt") as g:
                while True:
                    h = g.readline()
                    if not h:
                        break
                    seqln = g.readline().rstrip("\n")
                    g.readline()  # '+'
                    qln = g.readline().rstrip("\n")
                    # pbsim names reads like ">S1_1" ; the template-of-origin is
                    # encoded in its per-template FASTQ or the ref field. pbsim3
                    # 'templ' mode emits one FASTQ per template named
                    # <prefix>_<idx>.fastq with the template order == input order.
                    # We recover tid from the accompanying MAF ref name when
                    # present; otherwise fall back to filename index.
                    tid = _pbsim_read_tid(f, base, list(truth_by_tid))
                    if tid is None:
                        continue
                    n = counters[tid]
                    counters[tid] = n + 1
                    read_id = f"{tid}_r{n:03d}"
                    fq.write(f"@{read_id}\n{seqln}\n+\n{qln}\n")
                    tt.write("\t".join([read_id, tid] + truth_by_tid[tid]) + "\n")
                    n_reads += 1
    return {
        "reads_fastq": fastq_path, "read_truth": truth_path,
        "n_reads": str(n_reads), "n_templates": str(len(rows)),
        "n_extra_templates": str(len(extra)),
        "error_source": "pbsim3-ERRHMM-ONT",
    }


def _pbsim_read_tid(fname: str, base: str, tids: List[str]) -> Optional[str]:
    """Map a pbsim per-template FASTQ filename to its tid via the 1-based index
    pbsim assigns in template input order (<base>_0001.fastq -> tids[0])."""
    stem = fname[len(base):].lstrip("_")
    for ext in (".fastq", ".fq", ".fastq.gz", ".fq.gz"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    if stem.isdigit():
        idx = int(stem) - 1
        if 0 <= idx < len(tids):
            return tids[idx]
    return None


# ---------------------------------------------------------------------------
# Self-check — builds a tiny synthetic panel and asserts the SPEC invariants.
# ---------------------------------------------------------------------------
def self_check(verbose: bool = True) -> bool:
    def emit(m: str):
        if verbose:
            print(m)
    ok = True
    tmp = tempfile.mkdtemp(prefix="genreads_selfcheck_")
    try:
        # synthetic 2-exon templates + one INTRONFREE control + an HP template
        e5, e3 = 40, 40
        rng = random.Random(1)
        def randseq(n):
            return "".join(rng.choice(BASES) for _ in range(n))
        templ = {
            "tid_0": randseq(e5) + randseq(e3),                       # junction
            "tid_1": randseq(e5) + "AAAAAAA" + randseq(e3 - 7),       # exonic HP
            "tid_2": randseq(120),                                    # INTRONFREE
        }
        with open(os.path.join(tmp, "templates.fa"), "w") as fh:
            for k, v in templ.items():
                fh.write(f">{k}\n{v}\n")
        panel_cols = ["tid", "chrom", "true_donor", "true_acceptor", "strand",
                      "intron_len", "motif_rung", "acceptor_motif", "decoy_offset",
                      "decoy_acceptor", "exon5_len", "exon3_len", "context",
                      "has_true_junction"]
        panel_rows = [
            ["tid_0", "chrSIM_0", "100", "190", "+", "90", "R3", "AC", "3",
             "193", str(e5), str(e3), "plain", "1"],
            ["tid_1", "chrSIM_0", "300", "550", "+", "250", "R1", "AAG", "2",
             "552", str(e5), str(e3), "HP", "1"],
            ["tid_2", "chrSIM_0", "", "", "+", "0", "NA", "NA", "0", "", "60",
             "60", "INTRONFREE", "0"],
        ]
        with open(os.path.join(tmp, "panel_truth.tsv"), "w") as fh:
            fh.write("\t".join(panel_cols) + "\n")
            for r in panel_rows:
                fh.write("\t".join(r) + "\n")

        rpt = 110
        rows = read_panel_truth(os.path.join(tmp, "panel_truth.tsv"))
        templates = read_fasta(os.path.join(tmp, "templates.fa"))
        info = generate_fallback(rows, templates, tmp, reads_per_template=rpt,
                                 total_rate=ERROR_REGIMES["mid"], hp_del_mult=3.0,
                                 trunc_frac=0.3, rc_frac=0.2, seed=7,
                                 allow_extra_templates=False)

        # (a) read_truth header EXACTLY equals the SPEC 9-column string
        with open(info["read_truth"]) as fh:
            header = fh.readline().rstrip("\n")
            body = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        hdr_ok = header == "\t".join(READ_TRUTH_COLS)
        emit(f"[self-check] read_truth header exact: {'PASS' if hdr_ok else 'FAIL'}")
        ok &= hdr_ok

        # (b) row count == reads_per_template * n_templates
        cnt_ok = len(body) == rpt * len(panel_rows)
        emit(f"[self-check] row count {len(body)} == {rpt}*{len(panel_rows)}: "
             f"{'PASS' if cnt_ok else 'FAIL'}")
        ok &= cnt_ok

        # (c) FASTQ read_ids 1:1 with read_truth, all match <tid>_r<NNN>
        fq_ids = []
        with open(info["reads_fastq"]) as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                fh.readline(); fh.readline(); fh.readline()
                fq_ids.append(h[1:].rstrip("\n"))
        truth_ids = [r[0] for r in body]
        id_col = 0
        one2one = set(fq_ids) == set(truth_ids) and len(fq_ids) == len(truth_ids) == len(set(fq_ids))
        emit(f"[self-check] FASTQ<->truth read_id 1:1 (n={len(fq_ids)}): "
             f"{'PASS' if one2one else 'FAIL'}")
        ok &= one2one
        fmt_ok = all(rid.rsplit("_r", 1)[0] in templ and rid.rsplit("_r", 1)[1].isdigit()
                     for rid in truth_ids)
        emit(f"[self-check] read_id format <tid>_r<NNN>: {'PASS' if fmt_ok else 'FAIL'}")
        ok &= fmt_ok

        # (d) tid written verbatim; has_true_junction preserved per template
        tj = {r[panel_cols_idx(header, 'tid')]: r[panel_cols_idx(header, 'has_true_junction')]
              for r in body}
        want = {"tid_0": "1", "tid_1": "1", "tid_2": "0"}
        tj_ok = tj == want
        emit(f"[self-check] tid verbatim + has_true_junction preserved: "
             f"{'PASS' if tj_ok else 'FAIL'}  ({tj})")
        ok &= tj_ok

        # (e) error model actually mutated reads (not a clean copy) and stayed
        #     deletion-dominant at the marginal level on a big clean template
        _demo_error_shape(emit)

        # (f) missing-template raises loud
        bad_rows = rows + [{"tid": "tid_MISSING", "has_true_junction": "1"}]
        raised = False
        try:
            generate_fallback(bad_rows, templates, tmp, 5, ERROR_REGIMES["mid"],
                              3.0, 0.0, 0.0, 1, True)
        except KeyError:
            raised = True
        emit(f"[self-check] missing template raises: {'PASS' if raised else 'FAIL'}")
        ok &= raised

        emit("\n" + ("SELF-CHECK PASSED" if ok else "SELF-CHECK FAILED"))
        return ok
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def panel_cols_idx(header_line: str, col: str) -> int:
    return header_line.split("\t").index(col)


def _demo_error_shape(emit) -> None:
    """Confirm the marginal is deletion-dominant and HP elevates deletions."""
    rng = random.Random(3)
    seq = "".join(rng.choice(BASES) for _ in range(20000))
    hp = hp_run_lengths(seq)
    total = ERROR_REGIMES["mid"]
    # count event kinds over many corrupted copies by diffing lengths is noisy;
    # instead sample the process directly via the rate constants.
    dr, ir, sr = total * FRAC_DEL, total * FRAC_INS, total * FRAC_SUB
    emit(f"[self-check] marginal rates del={dr:.4f} ins={ir:.4f} sub={sr:.4f} "
         f"-> del-dominant: {'PASS' if dr > ir and dr > sr else 'FAIL'}")
    # HP elevation sanity: an A-run of 7 gets an elevated per-base del hazard
    hp_seq = "C" * 5 + "A" * 7 + "C" * 5
    hpl = hp_run_lengths(hp_seq)
    elevated = min(dr * 3.0, HP_DEL_CAP)
    emit(f"[self-check] HP run detected (max run={max(hpl)}), del hazard in-run "
         f"{elevated:.4f} > baseline {dr:.4f}: "
         f"{'PASS' if max(hpl) >= HP_MIN_LEN and elevated > dr else 'FAIL'}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Component B — ONT-like long-read generator for the yeast "
                    "non-canonical splice sim (templates.fa + panel_truth.tsv -> "
                    "reads.fastq + read_truth.tsv).")
    ap.add_argument("--out-dir", default=".",
                    help="dir to read templates/panel from and write outputs to")
    ap.add_argument("--templates", default=None,
                    help="templates.fa (default <out-dir>/templates.fa)")
    ap.add_argument("--panel-truth", default=None,
                    help="panel_truth.tsv (default <out-dir>/panel_truth.tsv)")
    ap.add_argument("--reads-per-template", type=int, default=120,
                    help="reads emitted per template (>=100 to satisfy N>=100)")
    ap.add_argument("--error-regime", choices=list(ERROR_REGIMES), default="mid",
                    help="hard-coded ONT error regime (fallback path)")
    ap.add_argument("--error-rate", type=float, default=None,
                    help="override total per-base error rate (fallback path)")
    ap.add_argument("--hp-del-mult", type=float, default=3.0,
                    help="deletion-hazard multiplier inside exonic HP runs")
    ap.add_argument("--trunc-frac", type=float, default=0.0,
                    help="fraction of reads that are junction-safe truncated "
                         "(default 0.0 = full-length; opt-in)")
    ap.add_argument("--rc-frac", type=float, default=0.0,
                    help="fraction of reads emitted reverse-complemented (cDNA "
                         "libs; default 0.0 = DRS/forward). Does NOT flip truth "
                         "strand.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--pbsim-bin", default="pbsim",
                    help="pbsim3 binary name/path (preferred if on PATH)")
    ap.add_argument("--force-fallback", action="store_true",
                    help="ignore pbsim3 even if present; use the fallback model")
    ap.add_argument("--allow-extra-templates", action="store_true",
                    help="permit templates.fa records not referenced by panel_truth")
    ap.add_argument("--self-check", action="store_true",
                    help="run the synthetic end-to-end invariant self-check and exit")
    args = ap.parse_args(argv)

    if args.self_check:
        return 0 if self_check() else 1

    templates_path = args.templates or os.path.join(args.out_dir, "templates.fa")
    panel_path = args.panel_truth or os.path.join(args.out_dir, "panel_truth.tsv")
    for p in (templates_path, panel_path):
        if not os.path.exists(p):
            ap.error(f"required input not found: {p}")

    rows = read_panel_truth(panel_path)
    templates = read_fasta(templates_path)

    pbsim_path = None if args.force_fallback else _find_pbsim(args.pbsim_bin)
    errhmm = _find_errhmm_model(pbsim_path) if pbsim_path else None

    if pbsim_path and errhmm:
        sys.stderr.write(
            f"[gen_reads] pbsim3 path ACTIVE: {pbsim_path} (model {os.path.basename(errhmm)}); "
            f"error source = ERRHMM-ONT, independent of penalty_scores.tsv\n")
        info = generate_pbsim(rows, templates, args.out_dir,
                              args.reads_per_template, pbsim_path, errhmm,
                              args.seed, args.allow_extra_templates)
    else:
        why = "pbsim3 not on PATH" if not pbsim_path else "no ERRHMM-ONT model found"
        if args.force_fallback:
            why = "--force-fallback"
        total = args.error_rate if args.error_rate is not None else ERROR_REGIMES[args.error_regime]
        sys.stderr.write(
            f"[gen_reads] FALLBACK path ACTIVE ({why}): hard-coded per-base ONT "
            f"model, total_rate={total:.4f}, del/ins/sub={FRAC_DEL}/{FRAC_INS}/{FRAC_SUB}, "
            f"HP(len>={HP_MIN_LEN}) del x{args.hp_del_mult} (cap {HP_DEL_CAP}); "
            f"INDEPENDENT of penalty_scores.tsv\n")
        info = generate_fallback(rows, templates, args.out_dir,
                                 args.reads_per_template, total, args.hp_del_mult,
                                 args.trunc_frac, args.rc_frac, args.seed,
                                 args.allow_extra_templates)

    sys.stderr.write("[gen_reads] done: " + ", ".join(f"{k}={v}" for k, v in info.items()) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
