#!/usr/bin/env python3
"""SPIKE-IN FABRICATION panel builder + scorer (task #20, the precision half).

The recall harness (build_panel/gen_reads/run_arms) measures whether motif-blind
arm-B RECOVERS a non-canonical junction minimap2 flattens. That panel cannot
measure FABRICATION: its reads are spliced at the true non-canonical site, so
arm-B holds there and never drifts.

The real-data failure (dev/SUMNER_SMA_LEADS_B_splicing.md) is the MIRROR: a REAL
(canonical/dominant) junction gets DRIFTED to a false NON-CANONICAL site at a
NON-HP position 6-28bp away. Truth = canonical; false call = non-canonical.

This builds the mirror panel:
  * a CANONICAL GT-AG true junction, reads generated from it (spliced at truth);
  * a nearby NON-CANONICAL drift target at graded NON-HP distance (the sweep axis);
  * the drift target is NOT a true isoform — a drifted call is a NON-CANONICAL FP.
  * to seed the candidate pool (arm-B can only drift to a POOL member within
    max_boundary_shift), a controllable minor fraction of "contaminant" reads are
    spliced at the drift target (mimics minimap2 error-driven boundary variance
    that populates the real pool). This is a KNOB: --seed-frac. seed-frac=0 tests
    whether error alone (no pool seed) fabricates; seed-frac>0 tests the realistic
    pool-contaminated regime.

Truth is emitted in BOTH formats:
  * flat read_truth.tsv (gen_reads schema) — for reuse with paired_arm_test.
  * rich truth (truth_schema) — consumed by scorer.py:score_bam for the
    ambiguity-aware precision/FDR-by-canonicity fabrication metric.

Author: spike-in agent (task #20)
"""
from __future__ import annotations

import argparse
import os
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_here = os.path.dirname(os.path.abspath(__file__))
_repo = os.path.abspath(os.path.join(_here, "..", "..", ".."))
if _repo not in sys.path:
    sys.path.insert(0, _repo)

BASES = "ACGT"
LPAD = 40
RPAD = 40


def rand_seq(n: int, rng: random.Random) -> List[str]:
    return [rng.choice(BASES) for _ in range(n)]


def _no_run(seq: List[str], rng: random.Random, max_run: int = 3) -> List[str]:
    """Break any homopolymer run >= max_run so the panel is NON-HP by construction
    (the whole point: the drift is at NON-homopolymer positions, so the HP-drift
    guard cannot see it). Rewrites in place, returns seq."""
    i = 0
    n = len(seq)
    while i < n:
        j = i
        while j + 1 < n and seq[j + 1] == seq[i]:
            j += 1
        run = j - i + 1
        if run >= max_run:
            # break every (max_run-1)th base to a different base
            for k in range(i + max_run - 1, j + 1, max_run - 1):
                alt = rng.choice([b for b in BASES if b != seq[k]
                                  and (k == 0 or b != seq[k - 1])
                                  and (k == n - 1 or b != seq[k + 1])])
                seq[k] = alt
        i = j + 1
    return seq


class FabContig:
    def __init__(self, name: str, seq: str, true_row: Dict, drift_row: Dict,
                 true_template: str):
        self.name = name
        self.seq = seq
        self.true_row = true_row        # canonical truth (the real junction)
        self.drift_row = drift_row      # non-canonical drift target (FALSE site)
        self.true_template = true_template


def build_fab_contig(idx: int, drift_dist: int, e5: int, e3: int,
                     intron_len: int, rng: random.Random,
                     drift_side: str = "acceptor",
                     microhom_mismatch: float = 0.2) -> FabContig:
    """A CANONICAL GT-AG junction with a NON-CANONICAL drift target `drift_dist`
    bp away, at a NON-HP position, EMBEDDED IN LOCAL MICROHOMOLOGY.

    Layout (acceptor-side drift; the real acceptor is canonical AG, the drift
    target shifts the acceptor by drift_dist=k into a NON-AG dinucleotide):

        [lpad][exon1 : e5][GT..intron..AG][exon2 : e3][rpad]
                          ^D=donor        ^A=true acceptor (canonical AG)
                                          ^A+k = drift acceptor (non-canonical)

    THE MICROHOMOLOGY (the enabling structure — advisor-corrected): for a read
    truly spliced at canonical A to re-place to A+k, the aligner realigns the
    read's exon2 (== genome[A:A+e3]) starting at genome position A+k. The cost of
    that move = the mismatches between the read's first k exon2 bases (genome[A:A+k])
    and the reference bases they'd land on (genome[A+k:A+2k]). On RANDOM sequence
    that's ~k mismatches (huge) so the true site always wins — no drift. REAL
    drift-prone loci (ribosomal / highly-expressed, the SMA fabrication-enriched
    set) have LOCAL SEQUENCE SIMILARITY at the drift distance. We install it: make
    genome[A:A+k] a NEAR-COPY of genome[A+k:A+2k] with `microhom_mismatch` fraction
    of bases differing (0.0 = perfect repeat -> ambiguity-equivalent SLIDE, which
    normalize_junction collapses = NOT a genuine FP; ~0.15-0.3 = imperfect =
    evidence-competitive under ONT error yet a DISTINCT spliced product = a real
    FP). The mismatches keep it non-canonical and let pbsim error decide which
    reads drift — no base is hand-tuned to force a specific read over.

    Region is HP-scrubbed (max run 2) so this is a NON-HP drift the guard cannot catch.
    """
    name = f"chrFAB_{idx}"
    k = drift_dist
    lpad = rand_seq(LPAD, rng)
    exon1 = rand_seq(e5, rng)
    body_len = intron_len - 2 - 3
    body = rand_seq(body_len, rng)
    acc_tri = list("CAG")             # canonical YAG acceptor (tier 0)
    intron = ["G", "T"] + body + acc_tri
    assert len(intron) == intron_len
    exon2 = rand_seq(e3, rng)

    if drift_side == "acceptor":
        if not (2 <= k <= e3 // 2):
            raise ValueError(f"drift_dist {k} needs 2<=k<=e3//2 ({e3//2})")
        # --- microhomology: exon2[0:k] := near-copy of exon2[k:2k] ---
        # (exon2[j] is genome[A+j]; so exon2[0:k]==genome[A:A+k], exon2[k:2k]==
        #  genome[A+k:A+2k]). Copy exon2[k:2k] into exon2[0:k], then flip
        #  round(k*mismatch) positions so it's IMPERFECT (distinct, not a slide).
        for j in range(k):
            exon2[j] = exon2[k + j]
        n_mis = max(1, int(round(k * microhom_mismatch)))
        mis_pos = rng.sample(range(k), min(n_mis, k))
        for j in mis_pos:
            exon2[j] = rng.choice([b for b in BASES if b != exon2[j]])
        # force the drift-acceptor dinuc (exon2[k-2:k] == genome[A+k-2:A+k]) NON-canonical
        exon2[k - 2] = "T"
        exon2[k - 1] = "C"   # dinuc "TC" -> non-canonical acceptor at A+k

    seq_list = lpad + exon1 + intron + exon2 + rpad_list(rng)
    # scrub HP runs OUTSIDE the microhomology block (positions < A) and in padding;
    # inside the microhom block we tolerate runs<=2 already (random-derived).
    seq_list = _no_run(seq_list, rng, max_run=3)
    seq = "".join(seq_list)

    D = LPAD + e5
    A = D + intron_len
    seq = _restamp(seq, D, A, k, drift_side)

    assert seq[D:D + 2] == "GT", (seq[D:D + 2], name, "donor")
    assert seq[A - 2:A] == "AG", (seq[A - 2:A], name, "true acceptor")
    drift_A = A + k
    drift_dinuc = seq[drift_A - 2:drift_A]
    assert drift_dinuc != "AG", (drift_dinuc, name, "drift acceptor must be non-canonical")
    # GENUINE-FP GUARD (advisor): the drift must be a DISTINCT junction, not an
    # ambiguity-equivalent slide. Verify normalize_junction(D,A) != (D,A+k).
    import importlib
    _nj = importlib.import_module(
        "rectify.core.consensus.chimeric_consensus").normalize_junction
    if _nj(D, A, seq) == _nj(D, drift_A, seq):
        raise ValueError(
            f"{name}: drift @+{k} is ambiguity-EQUIVALENT to truth "
            f"(microhom_mismatch={microhom_mismatch} too low / perfect repeat) — "
            f"would score as a non-error slide, not a fabrication FP. Raise mismatch.")

    true_template = seq[D - e5:D] + seq[A:A + e3]
    true_row = {
        "chrom": name, "donor": D, "acceptor": A, "strand": "+",
        "dinuc": seq[A - 2:A], "canonical": True, "template": true_template,
        "e5": e5, "e3": e3, "intron_len": intron_len,
    }
    drift_row = {
        "chrom": name, "donor": D, "acceptor": drift_A, "strand": "+",
        "dinuc": drift_dinuc, "canonical": False, "drift_dist": drift_dist,
        "template": seq[D - e5:D] + seq[drift_A:drift_A + e3],
        "e5": e5, "e3": e3,
    }
    return FabContig(name, seq, true_row, drift_row, true_template)


_RPAD_RNG_CACHE: Dict[int, List[str]] = {}


def rpad_list(rng: random.Random) -> List[str]:
    return rand_seq(RPAD, rng)


def _restamp(seq: str, D: int, A: int, drift_dist: int, drift_side: str) -> str:
    """Re-stamp the load-bearing motifs after the HP scrub, keeping local context
    non-HP by choosing flanking bases that differ from the motif bases."""
    s = list(seq)
    s[D], s[D + 1] = "G", "T"
    s[A - 2], s[A - 1] = "A", "G"
    if drift_side == "acceptor":
        drift_A = A + drift_dist
        s[drift_A - 2], s[drift_A - 1] = "T", "C"   # non-canonical, non-HP
    # ensure no new HP created at the stamped positions: nudge immediate neighbours
    for p in (D - 1, D + 2, A - 3, A, (A + drift_dist) - 3, (A + drift_dist)):
        if 0 <= p < len(s):
            # break a 3-run centered near p
            if 0 < p < len(s) - 1 and s[p - 1] == s[p] == s[p + 1]:
                s[p] = next(b for b in BASES if b != s[p])
    return "".join(s)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Build the spike-in FABRICATION panel "
                                 "(canonical truth + non-HP non-canonical drift target).")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--drift-dists", default="6,10,14,20,28",
                    help="non-HP drift distances (bp) matching SNRPN/UBA1/PCBP2 6-28bp")
    ap.add_argument("--e5", type=int, default=80)
    ap.add_argument("--e3", type=int, default=80)
    ap.add_argument("--intron-len", type=int, default=120)
    ap.add_argument("--true-reads", type=int, default=400,
                    help="reads spliced at the TRUE canonical junction per contig")
    ap.add_argument("--seed-frac", type=float, default=0.05,
                    help="fraction of contaminant reads spliced at the DRIFT target "
                         "to SEED the candidate pool (mimics mm2 boundary variance). "
                         "0 = error-only (no pool seed).")
    ap.add_argument("--microhom-mismatch", type=float, default=0.2,
                    help="microhomology imperfection: fraction of the k drift-flank "
                         "bases that differ between the true and drift placements. "
                         "0.0 = perfect repeat (ambiguity slide, rejected); ~0.15-0.3 "
                         "= imperfect near-repeat (evidence-competitive, genuine FP). "
                         "The drift-ENABLING structure real drift-prone loci carry.")
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    dists = [int(x) for x in args.drift_dists.split(",") if x.strip()]
    contigs = [build_fab_contig(i, d, args.e5, args.e3, args.intron_len, rng,
                                microhom_mismatch=args.microhom_mismatch)
               for i, d in enumerate(dists)]

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # sim_ref.fa (genomic reference)
    with open(outdir / "sim_ref.fa", "w") as fh:
        for c in contigs:
            fh.write(f">{c.name}\n{c.seq}\n")

    # templates.fa + a flat panel_truth-like table for gen_reads. gen_reads needs
    # panel_truth with tid/chrom/true_donor/true_acceptor/motif_rung/context/
    # has_true_junction/exon5_len/exon3_len/n_reads. We emit ONE true-junction
    # template per contig (tid_<i>) at the canonical site, plus (if seed-frac>0) a
    # co-located drift template (tid_<i>_seed) at the drift target to seed the pool.
    n_seed = {}
    with open(outdir / "templates.fa", "w") as tfh, \
         open(outdir / "panel_truth.tsv", "w") as pfh:
        cols = ["tid", "chrom", "true_donor", "true_acceptor", "strand",
                "intron_len", "motif_rung", "acceptor_motif", "decoy_offset",
                "decoy_acceptor", "exon5_len", "exon3_len", "context",
                "has_true_junction", "n_reads"]
        pfh.write("\t".join(cols) + "\n")
        for i, c in enumerate(contigs):
            tr = c.true_row
            tid = f"tid_{i}"
            tfh.write(f">{tid}\n{tr['template']}\n")
            pfh.write("\t".join([
                tid, c.name, str(tr["donor"]), str(tr["acceptor"]), "+",
                str(tr["intron_len"]), "CANON", tr["dinuc"], "",
                str(c.drift_row["acceptor"]), str(tr["e5"]), str(tr["e3"]),
                f"FABD{c.drift_row['drift_dist']}", "1", str(args.true_reads)]) + "\n")
            if args.seed_frac > 0:
                seed_n = max(1, int(round(args.true_reads * args.seed_frac)))
                n_seed[c.name] = seed_n
                dr = c.drift_row
                stid = f"tid_{i}_seed"
                tfh.write(f">{stid}\n{dr['template']}\n")
                # seed reads carry the DRIFT junction as their "truth" so gen_reads
                # emits reads spliced there — but they are NOT scored (score only the
                # true-junction reads). We tag context SEED so the scorer/harness can
                # exclude them. Their has_true_junction=1 so gen_reads is happy.
                pfh.write("\t".join([
                    stid, c.name, str(dr["donor"]), str(dr["acceptor"]), "+",
                    str(dr["acceptor"] - dr["donor"]), "SEED", dr["dinuc"], "",
                    "", str(dr["e5"]), str(dr["e3"]), "SEED", "1", str(seed_n)]) + "\n")

    # rich truth (truth_schema) for scorer.py — ONLY the true canonical junctions
    # are registered as truth; the drift target is deliberately absent so a drifted
    # call scores as a NON-CANONICAL FP. Written AFTER gen_reads (needs read ids),
    # so we emit a sidecar mapping here that build_rich_truth.py consumes.
    with open(outdir / "fab_contig_truth.tsv", "w") as fh:
        fh.write("chrom\tdonor\tacceptor\tstrand\tcanonical\tdrift_acceptor\tdrift_dist\te5\te3\n")
        for c in contigs:
            tr, dr = c.true_row, c.drift_row
            fh.write("\t".join([c.name, str(tr["donor"]), str(tr["acceptor"]), "+",
                                "1", str(dr["acceptor"]), str(dr["drift_dist"]),
                                str(tr["e5"]), str(tr["e3"])]) + "\n")

    sys.stderr.write(
        f"[spikein_fab] {len(contigs)} contigs, drift dists {dists}, "
        f"true_reads={args.true_reads}, seed_frac={args.seed_frac} "
        f"(seed reads/contig {n_seed})\n"
        f"[spikein_fab] wrote sim_ref.fa, templates.fa, panel_truth.tsv, "
        f"fab_contig_truth.tsv -> {outdir}\n")
    for c in contigs:
        sys.stderr.write(
            f"    {c.name}: TRUE canonical (D={c.true_row['donor']},"
            f"A={c.true_row['acceptor']},{c.true_row['dinuc']}) | "
            f"DRIFT non-canon @+{c.drift_row['drift_dist']}bp "
            f"(A={c.drift_row['acceptor']},{c.drift_row['dinuc']})\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
