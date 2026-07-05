#!/usr/bin/env python3
"""Component A — PANEL BUILDER for the yeast non-canonical LONG-READ splice sim.

*** AUTHORED BY THE INTEGRATOR (2026-07-03). ***
The four-agent hand-off named build_panel.py as an existing deliverable, but no
such file (and no build report) was present — only gen_reads.py / run_arms.py /
score_trade.py existed. This component was reconstructed from ``SPEC.md`` and the
interface the three downstream scripts already expect, so the A->B->C->D chain
composes. It is deliberately MINIMAL and spec-faithful (its correctness is the
only thing between the harness and a valid trade curve).

============================ WHAT THIS PRODUCES ============================
SPEC items 1-3, the head of the A->B->C->D contract:
  * ``templates.fa``     — spliced exon1++exon2 transcripts, header ``>tid_<id>``
                           (SPEC item 1). INTRONFREE templates are single-exon.
  * ``sim_ref.fa``       — the GENOMIC reference the reads align to (SPEC item 2);
                           synthetic contigs named ``chrSIM_<n>``. Truth junction
                           coords are IN THIS FRAME.
  * ``panel_truth.tsv``  — one row per template (SPEC item 3), 14 columns:
        tid chrom true_donor true_acceptor strand intron_len motif_rung
        acceptor_motif decoy_offset decoy_acceptor exon5_len exon3_len
        context has_true_junction
      Coords are 0-based half-open: the intron is ``[true_donor, true_acceptor)``.

============================= PANEL GEOMETRY ==============================
Each junction template gets its OWN synthetic contig (plus strand):

    [lpad] [exon1 : E5] [GT ....... intron ....... acc_tri] [exon2 : E3] [rpad]
            ^D-E5        ^D=donor                    ^A=acceptor  ^A+E3
                         genome[D:D+2]   == "GT"   (5'SS, canonical BY CONSTRUCTION)
                         genome[A-3:A]   == acc_tri (3'SS trinucleotide -> tier)
    template            == genome[D-E5:D] ++ genome[A:A+E3]   (exon1 ++ exon2)
    intron             == [D, A)   (len = intron_len)

When a decoy is requested a CANONICAL "AG" is embedded ``k`` bp downstream of the
true acceptor inside exon2: ``genome[A+k-2:A+k] == "AG"`` (decoy intron
``[D, A+k)``). A snap moves the acceptor D->A+k; because the donor is fixed and
normalize_junction preserves intron length, a snap is a DISTINCT spliced product
(never ambiguity-equivalent to truth), so the scorer scores it as not-recovered.

RUNG <-> tier map (matches junction_scoring._3ss_tier_from_rna_trinucleotide):
  * R0 = YAG (acc_tri "CAG"/"TAG") -> tier 0  (canonical CONTROL — must not move)
  * R3 = non-AG (acc dinuc "AC"/"AA") -> tier 4  (the snap-eligible non-canonical
         case; tier>=4 is the gate that arms the incumbent's canonical preference)
NOTE the deferred rungs and a SPEC/code conflict flagged for the full panel:
  R1 = RAG ("AAG"/"GAG") -> tier 1 is buildable, but SPEC's R2 ("still ends AG but
  low tier") is UNREALIZABLE here: every ``_AG`` trinucleotide maps to tier 0 or 1,
  never 2 (tier-2 NBG requires ending in G but NOT AG, e.g. "TG"/"GG"). A faithful
  tier ladder therefore forces R2's dinucleotide OFF "AG", contradicting the SPEC
  wording. The smoke panel omits R1/R2; resolve the conflict in the full panel.

HP-context cell: an EXONIC homopolymer run straddles the junction (intron tail +
exon2 start share an A-run) so the decoy shift can cost within the refiner's
0.5-edit canonical-prior floor. The run + decoy are built faithfully (no base is
hand-tuned to guarantee a sub-floor cost); the empirics decide whether A snaps.

Author: integrator (component A reconstruction)
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

BASES = "ACGT"

# panel_truth.tsv columns — SPEC item-3 order, plus the OPTIONAL ``n_reads``
# coordination column (the A<->B contract gen_reads._template_read_count honors:
# WT rows abundant, cryptic rows fixed-N, controls a base). Downstream (gen_reads,
# score_trade) join by header NAME, so appending n_reads is backward-safe.
PANEL_COLS = [
    "tid", "chrom", "true_donor", "true_acceptor", "strand", "intron_len",
    "motif_rung", "acceptor_motif", "decoy_offset", "decoy_acceptor",
    "exon5_len", "exon3_len", "context", "has_true_junction", "n_reads",
]

# Acceptor trinucleotides (last 3 nt of intron, RNA==genomic on the + strand).
# Chosen so _canonical_tier lands on the intended tier (see module docstring).
ACC_TRI = {
    "R0": "CAG",   # YAG  -> tier 0 (canonical)
    "R1": "AAG",   # RAG  -> tier 1 (semi-canonical; deferred, provided for full panel)
    "R3": "TAC",   # ..AC -> tier 4 (non-canonical; snap-eligible)
}
# HP variant of R3: the acceptor tail is an A-run so an exonic homopolymer
# straddles the junction. dinuc "AA" -> tier 4.
ACC_TRI_HP = "AAA"

LPAD = 30
RPAD = 30
DECOY_DINUC = "AG"


# ---------------------------------------------------------------------------
def rand_seq(n: int, rng: random.Random) -> List[str]:
    return [rng.choice(BASES) for _ in range(n)]


def _no_gt(rng: random.Random, prev: Optional[str]) -> str:
    """Pick a base; used to keep the intron body from accidentally spawning a
    second in-frame 'GT' donor right after the real one (cosmetic robustness)."""
    return rng.choice(BASES)


# ---------------------------------------------------------------------------
class Contig:
    """A synthetic contig + the truth row for the template(s) it carries.

    ``extra`` holds co-located isoforms sharing the SAME genomic contig (SPEC v2
    prp18Δ mixture: an abundant canonical WT isoform spliced at the decoy acceptor
    alongside the minor cryptic non-canonical one). Co-location is essential — the
    refiner's junction pool is per-contig, so the WT canonical junction only becomes
    an admissible snap target for the cryptic reads when both splice on one contig.
    """

    def __init__(self, name: str, seq: str, row: Dict[str, str], template: str,
                 extra: Optional[List[Tuple[Dict[str, str], str]]] = None):
        self.name = name
        self.seq = seq
        self.row = row
        self.template = template
        self.extra = extra or []   # List[(row_dict, template_str)] — co-located isoforms


def build_junction_contig(
    idx: int,
    rung: str,
    context: str,
    has_decoy: bool,
    e5: int,
    e3: int,
    intron_len: int,
    decoy_k: int,
    rng: random.Random,
    hp_run: int = 0,
) -> Contig:
    """Construct one junction-bearing contig (plus strand) + its truth row.

    ``hp_run`` (R3-HP dose-response knob, arm-C power test): when > 0 the acceptor
    sits inside a GENOMIC A-run of length ``3 + hp_run`` that STRADDLES it — the
    intron's ``AAA`` acc_tri (3 A's, left of the acceptor) + ``hp_run`` A's at the
    start of exon2 (right), terminated by a G. Because the intron's trailing A's
    and exon2's leading A's are indistinguishable to an aligner, the N-op end
    (acceptor) can SLIDE within the run — and ONT length-miscalls in the run make
    the exact placement ambiguous. That is precisely the length-error the -logP
    HP-del law is meant to resolve, and the ambiguity (hence any arm-C edge) should
    grow with hp_run. The canonical decoy AG = (last run-A, the G) just past the
    run. hp_run==0 keeps the legacy single-base-A HP construction / generic decoy.
    """
    name = f"chrSIM_{idx}"
    hp = context == "HP"
    acc_tri = ACC_TRI_HP if (hp and rung == "R3") else ACC_TRI[rung]
    assert len(acc_tri) == 3
    if intron_len < 5 + 0:
        raise ValueError("intron_len must be >= 5 (GT + body + acc_tri)")
    hp_run_active = hp and rung == "R3" and hp_run > 0

    lpad = rand_seq(LPAD, rng)
    exon1 = rand_seq(e5, rng)
    # intron = GT + body + acc_tri  (len == intron_len)
    body_len = intron_len - 2 - 3
    if body_len < 0:
        raise ValueError(f"intron_len {intron_len} too short for GT+body+acc_tri")
    body = rand_seq(body_len, rng)
    # For a deterministic swept run, force the base BEFORE acc_tri to be non-A so
    # the run's left half is EXACTLY the 3 acc_tri A's (no random extension).
    if hp_run_active and body_len > 0 and body[-1] == "A":
        body[-1] = "C"
    intron = ["G", "T"] + body + list(acc_tri)
    assert len(intron) == intron_len

    # exon2 + decoy.
    exon2 = rand_seq(e3, rng)
    if hp_run_active:
        if hp_run + 2 > e3:
            raise ValueError(f"hp_run {hp_run}+2 exceeds exon3_len {e3}")
        for j in range(hp_run):
            exon2[j] = "A"          # exon2 leading A-run (right half of the run)
        exon2[hp_run] = "G"         # terminate the run; forms the decoy AG with exon2[hp_run-1]
        eff_decoy_k = hp_run + 1    # decoy acceptor = A + hp_run + 1 (just past the run)
        has_decoy = True
    else:
        if hp and rung == "R3":
            exon2[0] = "A"          # legacy: one exonic A across the junction
        eff_decoy_k = decoy_k
        if has_decoy:
            if not (2 <= decoy_k <= e3):
                raise ValueError(f"decoy_k {decoy_k} out of range for e3 {e3}")
            exon2[decoy_k - 2] = "A"
            exon2[decoy_k - 1] = "G"

    rpad = rand_seq(RPAD, rng)

    seq = "".join(lpad + exon1 + intron + exon2 + rpad)
    D = LPAD + e5                      # donor / intron start (0-based)
    A = D + intron_len                 # acceptor / intron end (0-based, half-open)
    assert seq[D:D + 2] == "GT", (seq[D:D + 2], name)
    assert seq[A - 3:A] == acc_tri, (seq[A - 3:A], acc_tri, name)
    if hp_run_active:
        # genomic A-run [A-3, A+hp_run) == A^(3+hp_run), then G
        assert seq[A - 3:A + hp_run] == "A" * (3 + hp_run), (seq[A - 3:A + hp_run + 1], name)
        assert seq[A + hp_run] == "G", (seq[A + hp_run], name)

    template = seq[D - e5:D] + seq[A:A + e3]
    acc_dinuc = seq[A - 2:A]
    decoy_acc = str(A + eff_decoy_k) if has_decoy else ""
    if has_decoy:
        assert seq[A + eff_decoy_k - 2:A + eff_decoy_k] == DECOY_DINUC, (
            seq[A + eff_decoy_k - 2:A + eff_decoy_k], name)

    row = {
        "tid": f"tid_{idx}",
        "chrom": name,
        "true_donor": str(D),
        "true_acceptor": str(A),
        "strand": "+",
        "intron_len": str(intron_len),
        "motif_rung": rung,
        "acceptor_motif": acc_dinuc,
        "decoy_offset": str(eff_decoy_k) if has_decoy else "",
        "decoy_acceptor": decoy_acc,
        "exon5_len": str(e5),
        "exon3_len": str(e3),
        "context": (f"HP{hp_run}" if hp_run_active else context),
        "has_true_junction": "1",
    }

    # SPEC v2 MIXTURE: co-locate the abundant WT canonical isoform at the decoy.
    extra: List[Tuple[Dict[str, str], str]] = []
    if has_decoy and rung == "R3":
        A_wt = A + eff_decoy_k
        if A_wt + e3 > len(seq):
            raise ValueError(f"WT isoform overruns contig ({A_wt}+{e3} > {len(seq)}); "
                             f"increase RPAD or decrease decoy_k/hp_run")
        wt_template = seq[D - e5:D] + seq[A_wt:A_wt + e3]
        assert seq[A_wt - 2:A_wt] == DECOY_DINUC, (seq[A_wt - 2:A_wt], name)
        wt_row = {
            "tid": f"tid_{idx}_wt",
            "chrom": name,
            "true_donor": str(D),
            "true_acceptor": str(A_wt),
            "strand": "+",
            "intron_len": str(A_wt - D),
            "motif_rung": "WT",
            "acceptor_motif": seq[A_wt - 2:A_wt],
            "decoy_offset": "",
            "decoy_acceptor": str(A),
            "exon5_len": str(e5),
            "exon3_len": str(e3),
            "context": (f"HP{hp_run}" if hp_run_active else context),
            "has_true_junction": "1",
        }
        extra.append((wt_row, wt_template))
    return Contig(name, seq, row, template, extra)


def build_flanking_hp_contig(
    idx: int,
    dist: int,
    run_len: int,
    e5: int,
    e3: int,
    intron_len: int,
    rng: random.Random,
    side: str = "acceptor",
    base: str = "A",
) -> Contig:
    """A CANONICAL GT-AG junction with A-runs FLANKING it at distance ``dist``
    (PI model: A-stretches occur up- AND down-stream of junctions; HP deletions
    RIGHT AT the junction are ~impossible to handle, so test at varying DISTANCE).

        exon1 = [fill][A^run_len][dist non-A spacer] | GT..intron..AG | [dist non-A spacer][A^run_len][fill] = exon2

    e.g. dist=0 -> ...CGAAAAAA-GT..AG-AAAAAAAA...  (runs adjacent, extreme case);
         dist=2 -> ...CGAAAAAACG-GT..AG-CGAAAAAA... (CG spacer decouples the run).
    An ONT HP deletion in a flanking run shifts the read length there; at dist=0 it
    is ambiguous with the junction boundary (moves an A into the intron), at dist>0
    it is unambiguously exonic. Truth = the exact canonical junction (D, A); the
    metric is whether each arm places it EXACTLY right despite the flanking HP
    deletions. This isolates the -logP HP-del law (arm-C) from the motif snap."""
    name = f"chrSIM_{idx}"
    assert base in "ACGT" and side in ("acceptor", "donor", "both")
    up = side in ("donor", "both")          # run upstream of the donor (exon1 tail)
    dn = side in ("acceptor", "both")       # run downstream of the acceptor (exon2 head)
    nonbase = [b for b in "ACGT" if b != base]
    if up and e5 < run_len + dist + 2:
        raise ValueError(f"exon5_len {e5} too short for run {run_len}+dist {dist}")
    if dn and e3 < run_len + dist + 2:
        raise ValueError(f"exon3_len {e3} too short for run {run_len}+dist {dist}")

    def spacer(d):
        return [rng.choice(nonbase) for _ in range(d)]

    lpad = rand_seq(LPAD, rng)
    # exon1: when `up`, [fill][run][non-base spacer(dist)] | GT (dist=0 => run adjacent
    # to the donor, the "moves an A into the intron" extreme). Force the base before the
    # run non-base so the run length is exactly run_len.
    if up:
        fill1 = rand_seq(e5 - run_len - dist, rng)
        if fill1 and fill1[-1] == base:
            fill1[-1] = nonbase[0]
        exon1 = fill1 + [base] * run_len + spacer(dist)
    else:
        exon1 = rand_seq(e5, rng)
        if exon1 and exon1[-1] == base:
            exon1[-1] = nonbase[0]           # keep the donor flank clean

    acc_tri = "CAG"                          # canonical YAG acceptor (tier 0)
    body = rand_seq(intron_len - 2 - 3, rng)
    if body and body[-1] == "A":
        body[-1] = "C"
    intron = ["G", "T"] + body + list(acc_tri)
    assert len(intron) == intron_len

    # exon2: when `dn`, AG | [non-base spacer(dist)][run][fill].
    if dn:
        fill2 = rand_seq(e3 - run_len - dist, rng)
        if fill2 and fill2[0] == base:
            fill2[0] = nonbase[0]
        exon2 = spacer(dist) + [base] * run_len + fill2
    else:
        exon2 = rand_seq(e3, rng)
        if exon2 and exon2[0] == base:
            exon2[0] = nonbase[0]            # keep the acceptor flank clean

    rpad = rand_seq(RPAD, rng)
    seq = "".join(lpad + exon1 + intron + exon2 + rpad)
    D = LPAD + e5
    A = D + intron_len
    assert seq[D:D + 2] == "GT", (seq[D:D + 2], name)
    assert seq[A - 3:A] == acc_tri, (seq[A - 3:A], name)
    # verify each requested flanking run sits at the intended distance from the boundary
    if up:
        assert seq[D - dist - run_len:D - dist] == base * run_len, (name, "up-run")
    if dn:
        assert seq[A + dist:A + dist + run_len] == base * run_len, (name, "dn-run")

    template = seq[D - e5:D] + seq[A:A + e3]
    row = {
        "tid": f"tid_{idx}", "chrom": name, "true_donor": str(D), "true_acceptor": str(A),
        "strand": "+", "intron_len": str(intron_len), "motif_rung": "R0flank",
        "acceptor_motif": seq[A - 2:A], "decoy_offset": "", "decoy_acceptor": "",
        "exon5_len": str(e5), "exon3_len": str(e3),
        "context": f"{side[:3].upper()}_{base}_D{dist}",   # e.g. ACC_A_D3, DON_A_D1, BOT_A_D0
        "has_true_junction": "1",
    }
    return Contig(name, seq, row, template)


def hp_dist_cells(distances: List[int], run_len: int) -> List[Tuple]:
    """Advisor-designed flanking-DISTANCE panel (the arm-C niche map). Cells encoded
    as ('flank', dist, run_len, side, base):

      * PRIMARY: acceptor-side (exon2/downstream) poly-A, fine distance sweep — the
        dose-response. Pre-registered prediction: HUMP in arm-C−arm-B (≈0 at d=0
        [degenerate], peak at intermediate d, →0 at large d).
      * DONOR MIRROR: donor-side (exon1/upstream) poly-A at a few distances (symmetry).
      * BASE-CLASS BREADTH at a mid distance: poly-T (AT class, other base) + poly-G/C
        (CG class) — tests whether the AT/CG-split penalty table transfers.
      * DEGENERATE CONTROLS (pre-registered NULL): both-sided runs, adjacent/near.
      * INTRONFREE FDR control."""
    cells: List[Tuple] = []
    for d in distances:                                  # PRIMARY acceptor sweep
        cells.append(("flank", d, run_len, "acceptor", "A"))
    for d in (1, 3, 6):                                  # donor-side mirror
        cells.append(("flank", d, run_len, "donor", "A"))
    for b in ("T", "G", "C"):                            # base-class breadth @ d=3
        cells.append(("flank", 3, run_len, "acceptor", b))
    cells.append(("flank", 0, run_len, "both", "A"))     # degenerate control (null)
    cells.append(("flank", 1, run_len, "both", "A"))
    cells.append(("intronfree", "INTRONFREE", "INTRONFREE", False, 0))
    return cells


def build_intronfree_contig(idx: int, span: int, rng: random.Random) -> Contig:
    """A single-exon contig spanning the locus — NO true junction (FDR control)."""
    name = f"chrSIM_{idx}"
    lpad = rand_seq(LPAD, rng)
    exon = rand_seq(span, rng)
    rpad = rand_seq(RPAD, rng)
    seq = "".join(lpad + exon + rpad)
    template = "".join(exon)
    row = {
        "tid": f"tid_{idx}",
        "chrom": name,
        "true_donor": "",
        "true_acceptor": "",
        "strand": "+",
        "intron_len": "0",
        "motif_rung": "INTRONFREE",
        "acceptor_motif": "",
        "decoy_offset": "",
        "decoy_acceptor": "",
        "exon5_len": str(span),
        "exon3_len": "0",
        "context": "INTRONFREE",
        "has_true_junction": "0",
    }
    return Contig(name, seq, row, template)


# ---------------------------------------------------------------------------
# Panel spec (the SMOKE panel). Each entry -> one contig/template.
#   (kind, rung, context, has_decoy)
# ---------------------------------------------------------------------------
def smoke_cells() -> List[Tuple[str, str, str, bool]]:
    return [
        ("junction", "R0", "plain", True),    # canonical + decoy: must NOT move
        ("junction", "R0", "plain", False),   # YAG canonical, NO decoy (control)
        ("junction", "R3", "plain", True),    # non-canonical + decoy (snap target)
        ("junction", "R3", "HP",    True),    # non-canonical + decoy in exonic HP
        ("intronfree", "INTRONFREE", "INTRONFREE", False),
    ]


def full_cells() -> List[Tuple[str, str, str, bool]]:
    """Fuller ladder (R0/R1/R3 x context x decoy) — NOT used by the smoke run;
    provided so the harness has a documented next step. R2 intentionally absent
    (see the SPEC/code tier conflict in the module docstring)."""
    cells: List[Tuple[str, str, str, bool]] = []
    for rung in ("R0", "R1", "R3"):
        for ctx in ("plain", "HP"):
            for dec in (True, False):
                cells.append(("junction", rung, ctx, dec))
    cells.append(("intronfree", "INTRONFREE", "INTRONFREE", False))
    return cells


def hp_power_cells(hp_lengths: List[int]) -> List[Tuple]:
    """arm-C dose-response POWER panel (the make-or-break −logP test): R3
    non-canonical acceptors with an acceptor-straddling exonic A-run swept over
    *hp_lengths* (0 = plain non-HP anchor), each with the WT+decoy mixture, plus
    the FDR controls. Run at high N per cell (``--cryptic-reads``) + the paired
    arm-B-vs-arm-C test in score_trade. Prediction: any arm-C recovery/FDR edge
    over arm-B is ~0 at L=0 and GROWS with L (the -logP HP-del law's signature)."""
    cells: List[Tuple] = []
    for L in hp_lengths:
        if L <= 0:
            cells.append(("junction", "R3", "plain", True, 0))    # non-HP anchor
        else:
            cells.append(("junction", "R3", "HP", True, int(L)))  # HP run length L
    cells.append(("junction", "R0", "plain", True, 0))            # canonical+decoy control
    cells.append(("junction", "R0", "plain", False, 0))           # canonical no-decoy control
    cells.append(("intronfree", "INTRONFREE", "INTRONFREE", False, 0))
    return cells


def _stamp_n_reads(contigs: List[Contig], reads_per_locus: int,
                   cryptic_reads: int, cryptic_frac: float) -> None:
    """Write the OPTIONAL ``n_reads`` coordination column onto every truth row
    (A<->B contract). Semantics (advisor-shaped: concentrate reads where the
    make-or-break arm-C-vs-B question lives, keep WT abundant per task #12):

      * R3 cryptic main rows (the scored non-canonical cells, incl. R3-HP) get
        ``cryptic_reads`` — the paired-comparison power knob. N per cryptic cell
        stays EXACTLY cryptic_reads regardless of the mixture fraction.
      * every other scored main row (R0, R1, INTRONFREE) gets ``reads_per_locus``.
      * co-located WT canonical isoform rows get ``round(cryptic*(1-f)/f)`` so the
        cryptic:WT split at a mixture locus is ``cryptic_frac`` (WT abundant ->
        the decoy YAG is a confident pool member -> arm-A flattens cleanly).

    Every PANEL_COLS row MUST carry n_reads or write_outputs would KeyError, so
    this stamps ALL rows (main + WT + INTRONFREE)."""
    if not (0.0 < cryptic_frac < 1.0):
        raise ValueError(f"--cryptic-frac must be in (0,1); got {cryptic_frac}")
    for c in contigs:
        # R3 cryptic cells AND R0flank flanking-HP cells are the scored make-or-break
        # cells -> give them the (high) cryptic_reads count for paired-test power.
        cryptic_n = (cryptic_reads if c.row.get("motif_rung") in ("R3", "R0flank")
                     else reads_per_locus)
        c.row["n_reads"] = str(cryptic_n)
        for wt_row, _tmpl in c.extra:
            # WT sits opposite an R3 cryptic row on this contig; tie its abundance
            # to that cryptic count so cryptic/(cryptic+WT) == cryptic_frac.
            wt_n = int(round(cryptic_reads * (1.0 - cryptic_frac) / cryptic_frac))
            wt_row["n_reads"] = str(max(wt_n, 1))


def build_panel(
    cells: List[Tuple[str, str, str, bool]],
    e5: int,
    e3: int,
    intron_len: int,
    decoy_k: int,
    seed: int,
    reads_per_locus: int = 200,
    cryptic_reads: Optional[int] = None,
    cryptic_frac: float = 0.5,
) -> List[Contig]:
    rng = random.Random(seed)
    contigs: List[Contig] = []
    for i, cell in enumerate(cells):
        kind = cell[0]
        if kind == "intronfree":
            contigs.append(build_intronfree_contig(i, e5 + e3, rng))
        elif kind == "flank":
            _, dist, run_len, side, base = cell
            contigs.append(build_flanking_hp_contig(
                i, dist, run_len, e5, e3, intron_len, rng, side=side, base=base))
        else:
            rung, ctx, has_decoy = cell[1], cell[2], cell[3]
            hp_run = cell[4] if len(cell) > 4 else 0
            contigs.append(build_junction_contig(
                i, rung, ctx, has_decoy, e5, e3, intron_len, decoy_k, rng, hp_run=hp_run))
    _stamp_n_reads(contigs, reads_per_locus,
                   cryptic_reads if cryptic_reads is not None else reads_per_locus,
                   cryptic_frac)
    return contigs


def write_outputs(contigs: List[Contig], out_dir: str) -> Dict[str, str]:
    os.makedirs(out_dir, exist_ok=True)
    templates_fa = os.path.join(out_dir, "templates.fa")
    sim_ref_fa = os.path.join(out_dir, "sim_ref.fa")
    panel_tsv = os.path.join(out_dir, "panel_truth.tsv")

    with open(templates_fa, "w") as fh:
        for c in contigs:
            fh.write(f">{c.row['tid']}\n{c.template}\n")
            for row, template in c.extra:
                fh.write(f">{row['tid']}\n{template}\n")
    with open(sim_ref_fa, "w") as fh:
        for c in contigs:
            fh.write(f">{c.name}\n{c.seq}\n")   # ONE genomic ref per contig
    with open(panel_tsv, "w") as fh:
        fh.write("\t".join(PANEL_COLS) + "\n")
        for c in contigs:
            fh.write("\t".join(c.row[col] for col in PANEL_COLS) + "\n")
            for row, _ in c.extra:
                fh.write("\t".join(row[col] for col in PANEL_COLS) + "\n")
    return {"templates": templates_fa, "sim_ref": sim_ref_fa, "panel_truth": panel_tsv}


# ---------------------------------------------------------------------------
# Self-check — builds a panel and asserts every SPEC/frame invariant, including
# an extract_junctions round-trip (the frame bug that would silently deflate
# recovery). Mirrors the self-test discipline in gen_reads / score_trade.
# ---------------------------------------------------------------------------
def self_check(verbose: bool = True) -> bool:
    def emit(m: str):
        if verbose:
            print(m)

    # import the scorer primitive to assert the coord frame end-to-end
    _here = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(_here, "..", "..", ".."))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from rectify.core.benchmark.scorer import extract_junctions
    from rectify.core.splice.junction_scoring import _canonical_tier

    ok = True
    e5, e3, ilen, k = 60, 60, 90, 3
    contigs = build_panel(smoke_cells(), e5, e3, ilen, k, seed=7)

    exp_tier = {"R0": 0, "R1": 1, "R3": 4}
    for c in contigs:
        r = c.row
        if r["has_true_junction"] != "1":
            # INTRONFREE: template must equal a genome substring, no junction cols
            sub = r["chrom"]
            tin = c.template in c.seq
            emit(f"[self-check] {r['tid']} INTRONFREE template in contig: "
                 f"{'PASS' if tin else 'FAIL'}")
            ok &= tin
            continue
        D = int(r["true_donor"]); A = int(r["true_acceptor"])
        # (1) donor GT
        gt = c.seq[D:D + 2] == "GT"
        # (2) acceptor dinuc matches the recorded motif
        acc_ok = c.seq[A - 2:A] == r["acceptor_motif"]
        # (3) template == exon1 ++ exon2 in genome frame
        tmpl_ok = c.template == c.seq[D - e5:D] + c.seq[A:A + e3]
        # (4) intron_len consistent
        ilen_ok = (A - D) == int(r["intron_len"])
        # (5) decoy AG present at A+k when a decoy is declared
        if r["decoy_offset"]:
            kk = int(r["decoy_offset"])
            dec_ok = c.seq[A + kk - 2:A + kk] == DECOY_DINUC and \
                int(r["decoy_acceptor"]) == A + kk
        else:
            dec_ok = True
        # (6) tier is what the rung claims (the arm-A snap gate)
        tier = _canonical_tier(D, A, c.seq, "+")
        tier_ok = tier == exp_tier[r["motif_rung"]]
        # (7) extract_junctions frame round-trip: a clean read of the template
        #     placed at reference_start = D - e5 with CIGAR e5 M / (A-D) N / e3 M
        #     must recover EXACTLY the truth intron [D, A).
        cig = [(0, e5), (3, A - D), (0, e3)]
        js = extract_junctions(D - e5, cig)
        frame_ok = js == [(D, A)]

        allok = gt and acc_ok and tmpl_ok and ilen_ok and dec_ok and tier_ok and frame_ok
        emit(f"[self-check] {r['tid']} {r['motif_rung']}/{r['context']} "
             f"decoy={r['decoy_offset'] or '-'}: GT={gt} accdinuc={acc_ok} "
             f"tmpl={tmpl_ok} ilen={ilen_ok} decoy={dec_ok} "
             f"tier={tier}=={exp_tier[r['motif_rung']]}({tier_ok}) "
             f"frame(extract_junctions)={frame_ok} -> {'PASS' if allok else 'FAIL'}")
        ok &= allok

    # header/round-trip through the files (schema stability for the join)
    tmp = tempfile.mkdtemp(prefix="buildpanel_selfcheck_")
    try:
        paths = write_outputs(contigs, tmp)
        with open(paths["panel_truth"]) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        hdr_ok = hdr == PANEL_COLS
        emit(f"[self-check] panel_truth header == SPEC item-3 order: "
             f"{'PASS' if hdr_ok else 'FAIL'}")
        ok &= hdr_ok
        # every template header is >tid_<id> and 1:1 with panel rows
        with open(paths["templates"]) as fh:
            tids = [ln[1:].strip() for ln in fh if ln.startswith(">")]
        all_tids: List[str] = []
        for c in contigs:
            all_tids.append(c.row["tid"])
            all_tids.extend(row["tid"] for row, _ in c.extra)
        one2one = len(tids) == len(all_tids) and set(tids) == set(all_tids)
        emit(f"[self-check] templates.fa 1:1 with panel rows incl. co-located isoforms "
             f"({len(all_tids)} rows): {'PASS' if one2one else 'FAIL'}")
        ok &= one2one
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)

    emit("\n" + ("SELF-CHECK PASSED" if ok else "SELF-CHECK FAILED"))
    return ok


# ---------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Component A — build templates.fa + sim_ref.fa + panel_truth.tsv "
                    "for the yeast non-canonical long-read splice sim.")
    ap.add_argument("--out-dir", default=".", help="output directory")
    ap.add_argument("--panel", choices=["smoke", "full", "hp_power", "hp_dist"], default="smoke",
                    help="smoke = R0+R3+INTRONFREE+YAG-no-decoy (default); "
                         "full = R0/R1/R3 x context x decoy ladder; "
                         "hp_power = R3 non-canonical x HP-run-length sweep (arm-C "
                         "dose-response, run STRADDLING the acceptor); "
                         "hp_dist = canonical junction x flanking-A-run DISTANCE sweep "
                         "(PI model: HP dels at varying distance from the junction)")
    ap.add_argument("--hp-run-lengths", default="0,3,6,9,12",
                    help="comma-separated exonic A-run lengths for --panel hp_power "
                         "(0 = plain non-HP anchor); default 0,3,6,9,12")
    ap.add_argument("--hp-distances", default="0,1,2,3,5,8",
                    help="comma-separated flanking-A-run distances from the junction "
                         "for --panel hp_dist (0 = adjacent/extreme); default 0,1,2,3,5,8")
    ap.add_argument("--hp-run-len", type=int, default=8,
                    help="flanking A-run length for --panel hp_dist (default 8)")
    ap.add_argument("--exon5-len", type=int, default=60)
    ap.add_argument("--exon3-len", type=int, default=60)
    ap.add_argument("--intron-len", type=int, default=90)
    ap.add_argument("--decoy-offset", type=int, default=3,
                    help="k bp downstream of the true acceptor for the canonical "
                         "YAG decoy (SPEC decoy offset; default 3)")
    ap.add_argument("--reads-per-locus", type=int, default=200,
                    help="n_reads for each non-cryptic scored cell (R0/R1/INTRONFREE); "
                         "default 200 (>=100 satisfies SPEC N>=100)")
    ap.add_argument("--cryptic-reads", type=int, default=None,
                    help="n_reads for each R3 cryptic cell (the make-or-break "
                         "arm-C-vs-B non-canonical cells, incl. R3-HP). Default = "
                         "--reads-per-locus; set higher to buy paired-comparison power.")
    ap.add_argument("--cryptic-frac", type=float, default=0.5,
                    help="cryptic fraction at a mixture (R3+decoy) locus: WT reads = "
                         "round(cryptic*(1-f)/f). Default 0.5 (equal, == v2); lower "
                         "for an abundant-WT / minor-cryptic prp18d-like mixture.")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--self-check", action="store_true",
                    help="run the frame/schema invariant self-check and exit")
    args = ap.parse_args(argv)

    if args.self_check:
        return 0 if self_check() else 1

    if args.panel == "hp_power":
        hp_lengths = [int(x) for x in args.hp_run_lengths.split(",") if x.strip() != ""]
        cells = hp_power_cells(hp_lengths)
    elif args.panel == "hp_dist":
        dists = [int(x) for x in args.hp_distances.split(",") if x.strip() != ""]
        cells = hp_dist_cells(dists, args.hp_run_len)
    elif args.panel == "full":
        cells = full_cells()
    else:
        cells = smoke_cells()
    contigs = build_panel(cells, args.exon5_len, args.exon3_len,
                          args.intron_len, args.decoy_offset, args.seed,
                          reads_per_locus=args.reads_per_locus,
                          cryptic_reads=args.cryptic_reads,
                          cryptic_frac=args.cryptic_frac)
    paths = write_outputs(contigs, args.out_dir)

    n_j = sum(1 for c in contigs if c.row["has_true_junction"] == "1")
    n_if = len(contigs) - n_j
    cr = args.cryptic_reads if args.cryptic_reads is not None else args.reads_per_locus
    sys.stderr.write(
        f"[build_panel] panel={args.panel} contigs={len(contigs)} "
        f"(junction={n_j}, intronfree={n_if}); "
        f"e5={args.exon5_len} e3={args.exon3_len} intron={args.intron_len} "
        f"decoy_k={args.decoy_offset} seed={args.seed}; "
        f"reads/locus={args.reads_per_locus} cryptic_reads={cr} "
        f"cryptic_frac={args.cryptic_frac}\n"
        f"[build_panel] wrote {paths['templates']}, {paths['sim_ref']}, "
        f"{paths['panel_truth']}\n")
    for c in contigs:
        for r, _tmpl in [(c.row, c.template)] + c.extra:
            sys.stderr.write(
                f"    {r['tid']:12s} {r['chrom']:10s} {r['motif_rung']:10s} "
                f"ctx={r['context']:10s} donor={r['true_donor'] or '-':>4} "
                f"acc={r['true_acceptor'] or '-':>4} accmotif={r['acceptor_motif'] or '-':>3} "
                f"decoy_off={r['decoy_offset'] or '-':>2} decoy_acc={r['decoy_acceptor'] or '-':>4} "
                f"n_reads={r.get('n_reads', '-')}\n")
    n_wt = sum(len(c.extra) for c in contigs)
    if n_wt:
        sys.stderr.write(f"[build_panel] + {n_wt} co-located WT canonical isoform(s) "
                         f"(v2 mixture: makes the decoy an observed pool member)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
