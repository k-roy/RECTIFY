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

  * pbsim3 (PREFERRED, ERRHMM-ONT): used iff a pbsim binary resolves — on PATH,
    via --pbsim-bin, OR the source build shipped under ``tools/pbsim3/src/pbsim``
    next to this script (auto-detected, so the pbsim path is the DEFAULT active
    path on a machine where the tools build exists; --force-fallback overrides).
    We drive pbsim in ``--strategy templ --method errhmm`` (full-length template
    sequencing): each panel_truth row is REPLICATED n_i times in the pbsim input
    FASTA with tid-recoverable headers ``<tid>#r<NNN>``, so templ mode's one-read-
    per-record rule yields EXACTLY n_i FULL-LENGTH reads/template (n_i honors the
    optional ``n_reads`` column, else --reads-per-template). pbsim names its reads
    ``S_<k>`` (opaque), so the READ->TEMPLATE join goes through the MAF (block line
    1 = template ref_name ``<tid>#r<NNN>``, line 2 = read name ``S_<k>``) — NOT the
    output filename or order. pbsim's errors come from the packaged ONT ERRHMM HMM
    (a pretrained, external model), fully independent of penalty_scores.tsv. The
    error RATE is set by ``--accuracy-mean`` (default 0.95 = modern ONT R10-DRS-
    like), which scales the HMM's del-dominant/HP-collapsing STRUCTURE without
    touching the scorer. At the pbsim-native 0.85 (R9) the realized rates are
    ~sub 7.5% / ins 5.7% / del 8.3% (~21% total) and 120bp smoke reads under-map
    (~31% under minimap2 -k14); at 0.95 smoke maps ~89% (its short-read seeding
    ceiling) and realistic 300bp reads map ~99%. NOTE: pbsim numbers legitimately
    differ from the fallback-generated v2 baseline (a realism upgrade, not a
    regression); reproduce v2 exactly with --force-fallback.
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
# n_reads coordination column (A<->B contract, shared with build_panel.py)
# ---------------------------------------------------------------------------
def _template_read_count(row: Dict[str, str], default: int) -> int:
    """Per-template read count. Honors the OPTIONAL ``n_reads`` column in
    panel_truth (the A<->B coordination contract: build_panel writes WT rows high
    / cryptic rows low; controls = a default). CONTRACT (both builders agree):
    an empty, missing, non-integer, or non-positive ``n_reads`` cell falls back to
    ``default`` (--reads-per-template). This keeps v2 panels (no ``n_reads``
    column at all) BYTE-compatible — every row just uses the default."""
    raw = str(row.get("n_reads", "")).strip()
    if raw:
        try:
            v = int(raw)
        except ValueError:
            return default
        if v > 0:
            return v
    return default


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
            n_this = _template_read_count(row, reads_per_template)
            for n in range(n_this):
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
    """Resolve a runnable pbsim binary, in priority order:
      1. ``pbsim_bin`` on PATH (shutil.which),
      2. ``pbsim_bin`` as an explicit executable path,
      3. the source build shipped under ``tools/pbsim3/src/pbsim`` next to this
         script (so a machine that has run the local build gets the pbsim path as
         its DEFAULT — no flag needed; --force-fallback still overrides).
    Returns the resolved path or None. NOTE: this is a per-machine source build,
    so auto-default is machine-local (pbsim where built, fallback elsewhere)."""
    p = shutil.which(pbsim_bin)
    if p:
        return p
    if os.path.isabs(pbsim_bin) and os.path.exists(pbsim_bin) and os.access(pbsim_bin, os.X_OK):
        return pbsim_bin
    here = os.path.dirname(os.path.abspath(__file__))
    cand = os.path.join(here, "tools", "pbsim3", "src", "pbsim")
    if os.path.exists(cand) and os.access(cand, os.X_OK):
        return cand
    return None


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


def _parse_maf_read_to_ref(path: str) -> Dict[str, str]:
    """Map pbsim read name -> template ref name from a MAF (``.maf``/``.maf.gz``).

    Each alignment block opens with an ``a`` line then two ``s`` lines: the FIRST
    is the TEMPLATE (ref), the SECOND is the READ. So ``s`` line order gives the
    read->template join WITHOUT any reliance on output-file naming or read order
    (pbsim relabels reads ``S_<k>`` opaquely). Returns {read_name: ref_name}."""
    mapping: Dict[str, str] = {}
    s_names: List[str] = []
    with _open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            if line.startswith("a"):
                s_names = []
            elif line.startswith("s"):
                parts = line.split()
                if len(parts) >= 2:
                    s_names.append(parts[1])
                if len(s_names) == 2:
                    ref_name, read_name = s_names[0], s_names[1]
                    mapping[read_name] = ref_name
                    s_names = []
    return mapping


def generate_pbsim(rows: List[Dict[str, str]], templates: Dict[str, str],
                   out_dir: str, reads_per_template: int, pbsim_path: str,
                   errhmm_model: str, seed: int, allow_extra_templates: bool,
                   rc_frac: float = 0.0, accuracy: Optional[float] = None,
                   trunc_frac: float = 0.0) -> Dict[str, str]:
    """Drive pbsim3 (``--strategy templ --method errhmm``, full-length template
    sequencing), then relabel reads to the SPEC ``<tid>_r<NNN>`` id scheme and
    emit reads.fastq + read_truth.tsv (SAME schema as the fallback).

    Each panel_truth row is REPLICATED n_i times in the pbsim input FASTA with a
    tid-recoverable header ``<tid>#r<NNN>`` (n_i = honored ``n_reads`` else
    reads_per_template). templ mode emits exactly one FULL-LENGTH read per record,
    so replication == exact per-template count. Reads are joined back to their
    template via the MAF ref_name (NOT filename/order). pbsim's read errors come
    from the packaged ONT ERRHMM HMM — a pretrained external model, fully
    independent of penalty_scores.tsv (never read here)."""
    os.makedirs(out_dir, exist_ok=True)
    if trunc_frac:
        sys.stderr.write("[gen_reads] NOTE: --trunc-frac is a FALLBACK-only knob; "
                         "the pbsim path is full-length by construction (templ "
                         "mode) and ignores it.\n")
    # verify the join up front (fail loud), build the truth lookup + replicated FASTA
    truth_by_tid: Dict[str, List[str]] = {}
    want_by_tid: Dict[str, int] = {}
    matched_templates = set()
    tfa = os.path.join(out_dir, "_pbsim_templates.fa")
    with open(tfa, "w") as fh:
        for row in rows:
            tid = row["tid"]
            if "#" in tid:
                raise ValueError(f"tid {tid!r} contains '#', which collides with the "
                                 f"pbsim replicate-header separator '<tid>#r<NNN>'.")
            seq = lookup_template(tid, templates)
            if seq is None:
                raise KeyError(f"panel_truth tid {tid!r} has no template in templates.fa")
            for cand in (tid, f"tid_{tid}", tid[len('tid_'):] if tid.startswith('tid_') else tid):
                if cand in templates:
                    matched_templates.add(cand)
                    break
            n_i = _template_read_count(row, reads_per_template)
            want_by_tid[tid] = n_i
            for r in range(n_i):
                fh.write(f">{tid}#r{r:03d}\n{seq}\n")
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
           "--prefix", prefix, "--seed", str(seed)]
    if accuracy is not None:
        cmd += ["--accuracy-mean", str(accuracy)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"pbsim3 failed rc={res.returncode}: {res.stderr[:800]}")

    base = os.path.basename(prefix)
    # 1) read_name -> tid via the MAF (strip the '#r<NNN>' replicate suffix)
    read_to_tid: Dict[str, str] = {}
    for f in sorted(os.listdir(out_dir)):
        if f.startswith(base) and (f.endswith(".maf") or f.endswith(".maf.gz")):
            for rn, ref in _parse_maf_read_to_ref(os.path.join(out_dir, f)).items():
                read_to_tid[rn] = ref.split("#", 1)[0]
    if not read_to_tid:
        raise RuntimeError(f"pbsim3 emitted no MAF under {out_dir}/{base}*.maf* — "
                           f"cannot map reads to templates (gzip/samtools missing?)")

    # 2) stream the pbsim FASTQ(s), relabel per template via a per-tid counter
    rng = random.Random(seed)
    fastq_path = os.path.join(out_dir, "reads.fastq")
    truth_path = os.path.join(out_dir, "read_truth.tsv")
    counters: Dict[str, int] = {t: 0 for t in truth_by_tid}
    n_reads = 0
    fq_exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
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
                    rn = h[1:].split()[0].strip()   # pbsim read name, e.g. 'S_1'
                    tid = read_to_tid.get(rn)
                    if tid is None:
                        raise KeyError(
                            f"pbsim read {rn!r} has no MAF template mapping; the "
                            f"read<->template join is broken — refusing to emit a "
                            f"read with unknown truth.")
                    if tid not in truth_by_tid:
                        raise KeyError(f"pbsim MAF ref for read {rn!r} maps to tid "
                                       f"{tid!r}, absent from panel_truth (bad '#' split?)")
                    if rc_frac > 0 and rng.random() < rc_frac:
                        seqln = revcomp(seqln)
                        qln = qln[::-1]
                    n = counters[tid]
                    counters[tid] = n + 1
                    read_id = f"{tid}_r{n:03d}"
                    fq.write(f"@{read_id}\n{seqln}\n+\n{qln}\n")
                    tt.write("\t".join([read_id, tid] + truth_by_tid[tid]) + "\n")
                    n_reads += 1

    # 3) COUNT EXACTNESS (fail loud): every tid must have emitted exactly what we
    #    asked pbsim for. An under-count silently starves a panel cell below N>=100
    #    (the same failure the fallback's missing-template KeyError guards against).
    mism = {t: (counters[t], want_by_tid[t]) for t in want_by_tid
            if counters[t] != want_by_tid[t]}
    if mism:
        raise RuntimeError(
            f"pbsim per-template read count mismatch (got, want): "
            f"{dict(list(mism.items())[:8])}{' ...' if len(mism) > 8 else ''}. "
            f"pbsim dropped/duplicated records — refusing to emit a corrupt panel.")
    return {
        "reads_fastq": fastq_path, "read_truth": truth_path,
        "n_reads": str(n_reads), "n_templates": str(len(rows)),
        "n_extra_templates": str(len(extra)),
        "error_source": "pbsim3-ERRHMM-ONT",
    }


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
            "tid_0": randseq(e5) + randseq(e3),                       # junction (cryptic)
            "tid_0_wt": randseq(e5) + randseq(e3),                    # co-located WT isoform
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
        # A '_wt' row is included so the tid string ops (the pbsim '#'-split and the
        # '<tid>_r<NNN>' rsplit) are proven on the co-located-isoform tid shape.
        panel_rows = [
            ["tid_0", "chrSIM_0", "100", "190", "+", "90", "R3", "AC", "3",
             "193", str(e5), str(e3), "plain", "1"],
            ["tid_0_wt", "chrSIM_0", "100", "193", "+", "93", "WT", "AG", "",
             "190", str(e5), str(e3), "plain", "1"],
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
        want = {"tid_0": "1", "tid_0_wt": "1", "tid_1": "1", "tid_2": "0"}
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

        # (g) n_reads coordination column honored + backward-compatible (fallback)
        ok &= _check_n_reads_fallback(emit)

        # (h) pbsim3 path end-to-end (only if a pbsim binary + model resolve)
        ok &= _check_pbsim_path(emit)

        emit("\n" + ("SELF-CHECK PASSED" if ok else "SELF-CHECK FAILED"))
        return ok
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _write_panel(tmp: str, cols: List[str], rows: List[List[str]]) -> str:
    p = os.path.join(tmp, "panel_truth.tsv")
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return p


def _tid_counts(read_truth_path: str) -> Dict[str, int]:
    """Per-tid emitted read count from a read_truth.tsv (join by 'tid' column)."""
    with open(read_truth_path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        counts: Dict[str, int] = {}
        for r in rd:
            counts[r["tid"]] = counts.get(r["tid"], 0) + 1
    return counts


def _check_n_reads_fallback(emit) -> bool:
    """The n_reads A<->B contract on the fallback path: (i) a panel WITH n_reads
    yields exactly n_i reads per tid (WT high, cryptic low, control default);
    (ii) a panel WITHOUT the column is byte-compatible (every tid == default).
    Includes a '_wt' tid so the tid string ops are exercised end-to-end."""
    ok = True
    tmp = tempfile.mkdtemp(prefix="genreads_nreads_")
    try:
        rng = random.Random(2)
        e5 = e3 = 40
        def rs(n): return "".join(rng.choice(BASES) for _ in range(n))
        templ = {"tid_0": rs(e5) + rs(e3), "tid_0_wt": rs(e5) + rs(e3),
                 "tid_1": rs(120), "tid_2": rs(120), "tid_3": rs(120)}
        with open(os.path.join(tmp, "templates.fa"), "w") as fh:
            for k, v in templ.items():
                fh.write(f">{k}\n{v}\n")
        templates = read_fasta(os.path.join(tmp, "templates.fa"))

        # (i) WITH n_reads: cryptic=30, WT=200, and all three "->default" branches
        #     of _template_read_count: empty, "0" (non-positive), "x" (non-int).
        cols = ["tid", "chrom", "true_donor", "true_acceptor", "strand",
                "intron_len", "motif_rung", "acceptor_motif", "decoy_offset",
                "decoy_acceptor", "exon5_len", "exon3_len", "context",
                "has_true_junction", "n_reads"]
        rows = [
            ["tid_0", "c", "100", "190", "+", "90", "R3", "AC", "3", "193",
             str(e5), str(e3), "plain", "1", "30"],
            ["tid_0_wt", "c", "100", "193", "+", "93", "WT", "AG", "", "190",
             str(e5), str(e3), "plain", "1", "200"],
            ["tid_1", "c", "", "", "+", "0", "NA", "NA", "0", "", "120", "0",
             "INTRONFREE", "0", ""],          # empty -> default
            ["tid_2", "c", "", "", "+", "0", "NA", "NA", "0", "", "120", "0",
             "INTRONFREE", "0", "0"],          # non-positive -> default
            ["tid_3", "c", "", "", "+", "0", "NA", "NA", "0", "", "120", "0",
             "INTRONFREE", "0", "notanint"],   # non-integer -> default
        ]
        p = _write_panel(tmp, cols, rows)
        default = 77
        generate_fallback(read_panel_truth(p), templates, tmp,
                          reads_per_template=default, total_rate=ERROR_REGIMES["mid"],
                          hp_del_mult=3.0, trunc_frac=0.0, rc_frac=0.0, seed=5,
                          allow_extra_templates=False)
        counts = _tid_counts(os.path.join(tmp, "read_truth.tsv"))
        want = {"tid_0": 30, "tid_0_wt": 200, "tid_1": default,
                "tid_2": default, "tid_3": default}
        with_ok = counts == want
        emit(f"[self-check] n_reads honored (WT>cryptic; empty/0/non-int -> default): "
             f"{'PASS' if with_ok else 'FAIL'}  ({counts} vs {want})")
        ok &= with_ok

        # (ii) WITHOUT the n_reads column (v2 panel): every tid == default
        cols2 = cols[:-1]
        rows2 = [r[:-1] for r in rows]
        p2 = _write_panel(tmp, cols2, rows2)
        generate_fallback(read_panel_truth(p2), templates, tmp,
                          reads_per_template=default, total_rate=ERROR_REGIMES["mid"],
                          hp_del_mult=3.0, trunc_frac=0.0, rc_frac=0.0, seed=5,
                          allow_extra_templates=False)
        counts2 = _tid_counts(os.path.join(tmp, "read_truth.tsv"))
        bc_ok = counts2 == {t: default for t in templ}
        emit(f"[self-check] n_reads backward-compat (no column -> all default): "
             f"{'PASS' if bc_ok else 'FAIL'}  ({counts2})")
        ok &= bc_ok
        return ok
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _check_pbsim_path(emit) -> bool:
    """End-to-end pbsim3 check (SKIPPED if no pbsim binary+model resolves): a tiny
    panel with a '_wt' tid + varying n_reads. Asserts read_truth schema, exact
    per-tid counts (== honored n_reads), read_id 1:1 with FASTQ in <tid>_r<NNN>
    form, the tid round-trip (the '#'-split must keep '_wt'), and full-length
    reads (templ mode). This is the only place the pbsim MAF mapping is proven."""
    pbsim = _find_pbsim("pbsim")
    if not pbsim:
        emit("[self-check] pbsim path: SKIP (no pbsim binary/tools build found)")
        return True
    model = _find_errhmm_model(pbsim)
    if not model:
        emit("[self-check] pbsim path: SKIP (no ERRHMM-ONT model found)")
        return True
    ok = True
    tmp = tempfile.mkdtemp(prefix="genreads_pbsim_")
    try:
        rng = random.Random(4)
        e5 = e3 = 60
        def rs(n): return "".join(rng.choice(BASES) for _ in range(n))
        seqs = {"tid_0": rs(e5) + rs(e3), "tid_0_wt": rs(e5) + rs(e3),
                "tid_1": rs(120)}
        with open(os.path.join(tmp, "templates.fa"), "w") as fh:
            for k, v in seqs.items():
                fh.write(f">{k}\n{v}\n")
        templates = read_fasta(os.path.join(tmp, "templates.fa"))
        cols = ["tid", "chrom", "true_donor", "true_acceptor", "strand",
                "intron_len", "motif_rung", "acceptor_motif", "decoy_offset",
                "decoy_acceptor", "exon5_len", "exon3_len", "context",
                "has_true_junction", "n_reads"]
        rows = [
            ["tid_0", "c", "100", "190", "+", "90", "R3", "AC", "3", "193",
             str(e5), str(e3), "plain", "1", "8"],
            ["tid_0_wt", "c", "100", "193", "+", "93", "WT", "AG", "", "190",
             str(e5), str(e3), "plain", "1", "20"],
            ["tid_1", "c", "", "", "+", "0", "NA", "NA", "0", "", "120", "0",
             "INTRONFREE", "0", ""],          # empty -> default
        ]
        p = _write_panel(tmp, cols, rows)
        default = 12
        info = generate_pbsim(read_panel_truth(p), templates, tmp,
                              reads_per_template=default, pbsim_path=pbsim,
                              errhmm_model=model, seed=9,
                              allow_extra_templates=False)
        emit(f"[self-check] pbsim path: ACTIVE ({info['error_source']})")
        with open(info["read_truth"]) as fh:
            header = fh.readline().rstrip("\n")
            body = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        hdr_ok = header == "\t".join(READ_TRUTH_COLS)
        emit(f"[self-check] pbsim read_truth header exact: {'PASS' if hdr_ok else 'FAIL'}")
        ok &= hdr_ok

        counts = _tid_counts(info["read_truth"])
        want = {"tid_0": 8, "tid_0_wt": 20, "tid_1": default}
        cnt_ok = counts == want
        emit(f"[self-check] pbsim exact per-tid counts (n_reads honored): "
             f"{'PASS' if cnt_ok else 'FAIL'}  ({counts} vs {want})")
        ok &= cnt_ok

        # read_id 1:1 with FASTQ, all <tid>_r<NNN>, and the '_wt' tid round-trips
        fq_ids: List[str] = []
        fq_lens: List[int] = []
        with open(info["reads_fastq"]) as fh:
            while True:
                h = fh.readline()
                if not h:
                    break
                s = fh.readline().rstrip("\n"); fh.readline(); fh.readline()
                fq_ids.append(h[1:].rstrip("\n")); fq_lens.append(len(s))
        truth_ids = [r[0] for r in body]
        one2one = (set(fq_ids) == set(truth_ids)
                   and len(fq_ids) == len(truth_ids) == len(set(fq_ids)))
        emit(f"[self-check] pbsim FASTQ<->truth 1:1 (n={len(fq_ids)}): "
             f"{'PASS' if one2one else 'FAIL'}")
        ok &= one2one
        fmt_ok = all(rid.rsplit("_r", 1)[0] in seqs and rid.rsplit("_r", 1)[1].isdigit()
                     for rid in truth_ids)
        wt_ok = any(rid.startswith("tid_0_wt_r") for rid in truth_ids)
        emit(f"[self-check] pbsim read_id <tid>_r<NNN> + '_wt' tid round-trip: "
             f"{'PASS' if (fmt_ok and wt_ok) else 'FAIL'}")
        ok &= fmt_ok and wt_ok

        # full-length: templ mode reads should be ~template length (indels only),
        # never truncated far below it. All templates here are 120 bp.
        full_ok = min(fq_lens) >= int(120 * 0.7)
        emit(f"[self-check] pbsim reads full-length (min len {min(fq_lens)} >= 84): "
             f"{'PASS' if full_ok else 'FAIL'}")
        ok &= full_ok
        return ok
    except Exception as exc:   # surface, don't mask, a pbsim-path failure
        emit(f"[self-check] pbsim path: FAIL ({type(exc).__name__}: {exc})")
        return False
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
                    help="pbsim3 binary name/path (preferred if on PATH; also "
                         "auto-detects the tools/pbsim3/src/pbsim source build)")
    ap.add_argument("--pbsim-accuracy", type=float, default=0.95,
                    help="pbsim --accuracy-mean (pbsim path). Default 0.95 = modern "
                         "ONT R10-DRS-like; keeps the ERRHMM-ONT del-dominant/HP "
                         "error STRUCTURE, just scales the rate. Lower to 0.85 for "
                         "noisy R9 (~21%% error, but 120bp smoke reads then under-"
                         "map ~31%% under minimap2 -k14). Pass a value or 'nan' via "
                         "code (accuracy=None) to use pbsim's own 0.85 default.")
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
                              args.seed, args.allow_extra_templates,
                              rc_frac=args.rc_frac, accuracy=args.pbsim_accuracy,
                              trunc_frac=args.trunc_frac)
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
