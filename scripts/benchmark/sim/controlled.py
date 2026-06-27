#!/usr/bin/env python3
"""Tier-1 CONTROLLED-error generators — the DISCRIMINATING benchmark tier.

The brief forbids hand-rolling a *read simulator*; controlled error injection on
constructed loci with truth recorded BY CONSTRUCTION is a different thing (it is
what the existing ``dev/bench/hp_vertical_slice.py`` already does) and is the ONLY
way to get exact-indel-position truth that does not depend on a black-box
simulator's unknown injected edits. This module generalizes the vertical slice
into ``ReadTruth`` rows over a synthetic mini-genome, emitting:

* a reference FASTA (the synthetic contigs),
* a reads FASTQ (constructed reads with controlled errors), and
* a per-read truth table (``truth_schema.ReadTruth``).

Strata implemented here (the discriminating ones; per the SPEC's VERTICAL-SLICE
FINDING isolated runs are non-discriminating, so the HP stratum engineers the
HARD boundary/STR cases):

* ``HP``            homopolymer runs A/C/G/T x len 1-12, del-dominant length error,
                    flanked; every (base_class, run_len) cell sized >= min_count.
* ``STR``           short tandem repeats where >=2 indel placements are ED-tied.
* ``JUNCTION_AMB``  a GT-AG intron whose donor sits in a 1bp repeat, so a junction
                    called one bp off is ambiguity-EQUIVALENT (the GMAP-0.09 trap).
* ``CLEAN``         no-error reads — the false-indel-rate (FP) control.

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import random
import sys
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    ReadTruth, IndelTruth, JunctionTruth, JunctionClass, IndelKind, SplitTag,
    VariantTruth, write_truth_table, homopolymer_run, make_hp_indel,
    make_unique_indel, deletion_ambiguity_span,
)
from scripts.benchmark.sim.transcript_model import TranscriptModel, Exon  # noqa: E402

BASES = "ACGT"
FLANK = 80
HP_LENS = list(range(1, 13))
K_DIST = [(0, 0.50), (1, 0.30), (2, 0.14), (3, 0.06)]  # del-dominant under-call


def _rand_unique(n: int, rng: random.Random, avoid: str = "") -> str:
    """Random non-homopolymeric flank that does not start/end with ``avoid``."""
    out = []
    for i in range(n):
        choices = [b for b in BASES if b != (out[-1] if out else "")
                   and not (i in (0, n - 1) and b == avoid)]
        out.append(rng.choice(choices))
    return "".join(out)


def _draw_k(L: int, rng: random.Random) -> int:
    r, acc = rng.random(), 0.0
    for k, p in K_DIST:
        acc += p
        if r <= acc:
            return min(k, L - 1)
    return 0


def _split_for(locus_idx: int) -> SplitTag:
    """Genomic-region-disjoint split: even loci TRAIN, odd loci TEST (the ONE
    shared split — region, never read, decides it so a calibration win can't
    leak across the train/test boundary)."""
    return SplitTag.TRAIN if locus_idx % 2 == 0 else SplitTag.TEST


# ---------------------------------------------------------------------------
# HP stratum
# ---------------------------------------------------------------------------
def gen_hp_stratum(reps: int, rng: random.Random, locus0: int = 0
                   ) -> Tuple[Dict[str, str], List[Tuple[str, str]], List[ReadTruth]]:
    refs: Dict[str, str] = {}
    reads: List[Tuple[str, str]] = []
    truth: List[ReadTruth] = []
    li = locus0
    for b in BASES:
        for L in HP_LENS:
            chrom = f"hp_{b}_{L:02d}"
            lf = _rand_unique(FLANK, rng, avoid=b)
            rf = _rand_unique(FLANK, rng, avoid=b)
            ref = lf + b * L + rf
            refs[chrom] = ref
            run_start, run_end = FLANK, FLANK + L
            split = _split_for(li)
            li += 1
            for i in range(reps):
                k = _draw_k(L, rng)
                rid = f"{chrom}_r{i:03d}_k{k}"
                reads.append((rid, lf + b * (L - k) + rf))
                indels = []
                if k > 0:
                    es, ee = deletion_ambiguity_span(ref, run_start, k)
                    indels.append(IndelTruth(
                        pos=es, length=k, kind=IndelKind.DEL,
                        eq_start=es, eq_end=ee,
                        context="HP", run_unit=b, run_copies=L))
                truth.append(ReadTruth(
                    read_id=rid, true_locus=chrom, true_transcript=chrom,
                    chrom=chrom, strand="+", genome_start=0, genome_end=len(ref),
                    true_cigar=f"{len(ref) - k}M" if k == 0 else f"{run_start}M{k}D{len(ref) - run_start}M",
                    indels=indels, stratum="HP", split=split, coverage=reps,
                ))
    return refs, reads, truth


# ---------------------------------------------------------------------------
# HP_HARD stratum — the DISCRIMINATING HP cases (per SPEC VERTICAL-SLICE FINDING)
# ---------------------------------------------------------------------------
# An isolated, cleanly-flanked HP run is NON-DISCRIMINATING: minimap2 AND the live
# flat-affine DP both score 1.000 (any in-run gap placement is ambiguity-
# equivalent). The action — and C1's claimed value — is at the HARD boundary cases
# where the flat-affine indel-vs-substitution tradeoff DIVERGES from a calibrated
# length-law. These generators construct those, so the incumbent is BELOW ceiling
# (a candidate ablation only counts if minimap2 is below ceiling on the stratum).
def _mutate(seq: str, pos: int, rng: random.Random, avoid: str = "") -> str:
    """Substitute the base at ``pos`` for a DIFFERENT base, never ``avoid`` (so a
    boundary substitution can never be turned INTO the run base, which would
    silently change the effective run length and corrupt the truth)."""
    forbidden = {seq[pos].upper()} | ({avoid.upper()} if avoid else set())
    alt = rng.choice([b for b in BASES if b not in forbidden])
    return seq[:pos] + alt + seq[pos + 1:]


def gen_hp_hard_stratum(reps: int, rng: random.Random, locus0: int = 50,
                        sub_rate: float = 0.04
                        ) -> Tuple[Dict[str, str], List[Tuple[str, str]], List[ReadTruth]]:
    """Two hard sub-cases, both with truth = the run deletion ONLY (substitutions
    are NOT indels):

    * ``boundary_sub`` — a forced mismatch at the base IMMEDIATELY flanking the run
      + a length error. Flat-affine may absorb the boundary mismatch into the gap
      (an indel OUTSIDE the true run span) instead of calling a mismatch + an
      in-run deletion — the indel-vs-substitution tradeoff.
    * ``noisy`` — a length error + ``sub_rate`` random substitutions over the whole
      read (combined background noise), so the aligner sometimes misplaces the
      indel relative to the true run span or absorbs a nearby mismatch.

    The scorer credits a placement only if net (D-I) in the run span == k AND no
    indel falls outside the truth span, so a mis-absorbed mismatch scores incorrect.
    """
    refs: Dict[str, str] = {}
    reads: List[Tuple[str, str]] = []
    truth: List[ReadTruth] = []
    li = locus0
    hard_lens = [4, 6, 8, 10, 12]
    for b in BASES:
        for L in hard_lens:
            chrom = f"hph_{b}_{L:02d}"
            lf = _rand_unique(FLANK, rng, avoid=b)
            rf = _rand_unique(FLANK, rng, avoid=b)
            ref = lf + b * L + rf
            refs[chrom] = ref
            run_start, run_end = FLANK, FLANK + L
            split = _split_for(li)
            li += 1
            for i in range(reps):
                k = max(1, _draw_k(L, rng))   # hard cases always carry an error
                core = lf + b * (L - k) + rf
                mode = "boundary_sub" if (i % 2 == 0) else "noisy"
                if mode == "boundary_sub":
                    # mutate the base immediately 3' of the (shortened) run; never
                    # to the run base (would change the effective run length).
                    bpos = FLANK + (L - k)
                    if bpos < len(core):
                        core = _mutate(core, bpos, rng, avoid=b)
                else:
                    n_sub = max(1, int(round(sub_rate * len(core))))
                    for _ in range(n_sub):
                        p = rng.randrange(len(core))
                        # don't perturb the run or its immediate boundary (that
                        # would change the true run length / ambiguity)
                        if run_start - 1 <= p <= run_start + (L - k):
                            continue
                        core = _mutate(core, p, rng, avoid=b)
                rid = f"{chrom}_{mode}_r{i:03d}_k{k}"
                reads.append((rid, core))
                es, ee = deletion_ambiguity_span(ref, run_start, k)
                indels = [IndelTruth(
                    pos=es, length=k, kind=IndelKind.DEL,
                    eq_start=es, eq_end=ee,
                    context="HP", run_unit=b, run_copies=L)]
                truth.append(ReadTruth(
                    read_id=rid, true_locus=chrom, true_transcript=chrom,
                    chrom=chrom, strand="+", genome_start=0, genome_end=len(ref),
                    true_cigar=f"{run_start}M{k}D{len(ref) - run_start - k}M",
                    indels=indels, stratum="HP_HARD", split=split, coverage=reps))
    return refs, reads, truth


# ---------------------------------------------------------------------------
# STR stratum (edit-distance-tied placements)
# ---------------------------------------------------------------------------
def gen_str_stratum(reps: int, rng: random.Random, locus0: int = 100
                    ) -> Tuple[Dict[str, str], List[Tuple[str, str]], List[ReadTruth]]:
    refs: Dict[str, str] = {}
    reads: List[Tuple[str, str]] = []
    truth: List[ReadTruth] = []
    li = locus0
    units = ["AT", "AG", "CAG", "GT"]
    copies_list = [4, 6, 8]
    for unit in units:
        for ncopy in copies_list:
            chrom = f"str_{unit}_{ncopy}"
            lf = _rand_unique(FLANK, rng)
            rf = _rand_unique(FLANK, rng)
            run = unit * ncopy
            ref = lf + run + rf
            refs[chrom] = ref
            run_start, run_end = FLANK, FLANK + len(run)
            split = _split_for(li)
            li += 1
            for i in range(reps):
                # delete a whole unit half the time (the ED-tied case)
                drop = rng.random() < 0.4
                if drop and ncopy > 1:
                    read = lf + unit * (ncopy - 1) + rf
                    es, ee = deletion_ambiguity_span(ref, run_start, len(unit))
                    indels = [IndelTruth(
                        pos=es, length=len(unit), kind=IndelKind.DEL,
                        eq_start=es, eq_end=ee,
                        context="STR", run_unit=unit, run_copies=ncopy)]
                    cig = f"{run_start}M{len(unit)}D{len(read) - run_start}M"
                else:
                    read = ref
                    indels = []
                    cig = f"{len(ref)}M"
                rid = f"{chrom}_r{i:03d}_{'d' if (drop and ncopy>1) else 'c'}"
                reads.append((rid, read))
                truth.append(ReadTruth(
                    read_id=rid, true_locus=chrom, true_transcript=chrom,
                    chrom=chrom, strand="+", genome_start=0, genome_end=len(ref),
                    true_cigar=cig, indels=indels, stratum="STR",
                    split=split, coverage=reps))
    return refs, reads, truth


# ---------------------------------------------------------------------------
# JUNCTION_AMB stratum — a spliced gene whose donor sits in a 1bp repeat
# ---------------------------------------------------------------------------
def gen_junction_ambiguity_locus(rng: random.Random, n_reads: int = 30,
                                 locus_idx: int = 200
                                 ) -> Tuple[str, str, TranscriptModel, List[Tuple[str, str]], List[ReadTruth]]:
    """Build a contig with ONE intron whose flanking bases repeat so the junction
    has a >=1bp ambiguity window, plus full-length spliced reads. Returns
    (chrom, contig_seq, model, reads, truth). Used by the smoke's
    shifted-call-is-TP-not-FP assertion."""
    chrom = "chrTjunc"
    # A spliced gene with ONE GT-AG intron, sized so minimap2 -ax splice reliably
    # places it (200bp exons, 200bp intron). The donor sits in a 1bp repeat so a
    # junction called one bp left is ambiguity-EQUIVALENT (the GMAP-0.09 trap):
    # we engineer seq[intron_start-1] == seq[intron_end-1] so normalize_junction
    # slides a left-shifted call back onto truth.
    e1 = _rand_unique(200, rng)
    # exon1 ends with 'G'; intron starts 'GT'; base just inside the acceptor is
    # 'G' so seq[start-1]=='G'==seq[end-1] -> a 1bp left slide preserves product.
    e1 = e1[:-1] + "G"
    intron_core = _rand_unique(196, rng)
    intron = "GT" + intron_core + "AG"
    intron = intron[:-3] + "G" + "AG"
    e2 = _rand_unique(200, rng)
    contig = e1 + intron + e2
    e1_start = 0
    e1_end = len(e1)
    intron_start = e1_end
    intron_end = e1_end + len(intron)
    e2_start = intron_end
    e2_end = intron_end + len(e2)
    model = TranscriptModel(
        name="tjunc", chrom=chrom, strand="+",
        exons=[Exon(e1_start, e1_end), Exon(e2_start, e2_end)],
        genome_seq=contig)
    spliced = model.spliced_transcript()
    reads = []
    truth = []
    jt = model.junction_truths()  # NNC (no annotation supplied)
    split = _split_for(locus_idx)
    for i in range(n_reads):
        rid = f"{chrom}_full_r{i:03d}"
        reads.append((rid, spliced))
        truth.append(ReadTruth(
            read_id=rid, true_locus="tjunc", true_transcript="tjunc",
            chrom=chrom, strand="+",
            genome_start=model.genome_start, genome_end=model.genome_end,
            true_cigar=model.fulllength_cigar(), junctions=jt,
            true_cpa=model.genome_end - 1, stratum="JUNCTION_AMB",
            split=split, coverage=n_reads))
    return chrom, contig, model, reads, truth


# ---------------------------------------------------------------------------
# NIC / ANNOTATED junction loci — verify the discovery-class classifier
# ---------------------------------------------------------------------------
def gen_junction_class_loci(rng: random.Random, n_reads: int = 30, locus0: int = 300
                            ) -> Tuple[Dict[str, str], List[Tuple[str, str]], List[ReadTruth]]:
    """One ANNOTATED locus (read intron is in the annotation) and one NIC locus
    (exon-skipping: a novel pairing of two KNOWN sites), so the NIC/ANNOTATED
    branches of ``TranscriptModel.junction_truths`` are exercised end-to-end (the
    JUNCTION_AMB locus only covers NNC)."""
    refs: Dict[str, str] = {}
    reads: List[Tuple[str, str]] = []
    truth: List[ReadTruth] = []

    def gt_intron(core_len: int) -> str:
        return "GT" + _rand_unique(core_len, rng) + "AG"

    # ---- ANNOTATED: 2-exon transcript, intron is catalogued ----
    chrom = "chrAnno"
    e1 = _rand_unique(200, rng)
    intron = gt_intron(200)
    e2 = _rand_unique(200, rng)
    contig = e1 + intron + e2
    refs[chrom] = contig
    m = TranscriptModel("anno", chrom, "+",
                        [Exon(0, 200), Exon(200 + len(intron), len(contig))], contig)
    introns = m.introns()
    annotated_pairs = {(chrom, introns[0][0], introns[0][1])}
    jt = m.junction_truths(annotated_pairs=annotated_pairs)
    assert jt[0].klass == JunctionClass.ANNOTATED, jt[0].klass
    spliced = m.spliced_transcript()
    split = _split_for(locus0)
    for i in range(n_reads):
        rid = f"{chrom}_r{i:03d}"
        reads.append((rid, spliced))
        truth.append(ReadTruth(
            read_id=rid, true_locus="anno", true_transcript="anno", chrom=chrom,
            strand="+", genome_start=m.genome_start, genome_end=m.genome_end,
            true_cigar=m.fulllength_cigar(), junctions=jt,
            true_cpa=m.genome_end - 1, stratum="ANNOTATED", split=split,
            coverage=n_reads))

    # ---- NIC: 3-exon gene, read SKIPS exon2 (novel pairing of known sites) ----
    chrom = "chrNic"
    E1 = _rand_unique(200, rng)
    I1 = gt_intron(200)
    E2 = _rand_unique(200, rng)
    I2 = gt_intron(200)
    E3 = _rand_unique(200, rng)
    contig = E1 + I1 + E2 + I2 + E3
    refs[chrom] = contig
    e1s, e1e = 0, 200
    e2s = e1e + len(I1)
    e2e = e2s + 200
    e3s = e2e + len(I2)
    e3e = e3s + 200
    # annotation: the two ADJACENT introns (known donors + acceptors)
    annotated_donors = {(chrom, e1e), (chrom, e2e)}
    annotated_acceptors = {(chrom, e2s), (chrom, e3s)}
    annotated_pairs = {(chrom, e1e, e2s), (chrom, e2e, e3s)}
    # the READ skips exon2: exon1 + exon3 (novel donor e1e <-> acceptor e3s pairing)
    skip = TranscriptModel("nic", chrom, "+", [Exon(e1s, e1e), Exon(e3s, e3e)], contig)
    jt = skip.junction_truths(annotated_pairs=annotated_pairs,
                              annotated_donors=annotated_donors,
                              annotated_acceptors=annotated_acceptors)
    assert jt[0].klass == JunctionClass.NIC, jt[0].klass
    spliced = skip.spliced_transcript()
    split = _split_for(locus0 + 1)
    for i in range(n_reads):
        rid = f"{chrom}_skip_r{i:03d}"
        reads.append((rid, spliced))
        truth.append(ReadTruth(
            read_id=rid, true_locus="nic", true_transcript="nic_skip", chrom=chrom,
            strand="+", genome_start=skip.genome_start, genome_end=skip.genome_end,
            true_cigar=skip.fulllength_cigar(), junctions=jt,
            true_cpa=skip.genome_end - 1, stratum="NIC", split=split,
            coverage=n_reads))
    return refs, reads, truth


# ---------------------------------------------------------------------------
# VARIANT stratum (C6) — standing-variant-induced discovery FDR
# ---------------------------------------------------------------------------
# Per the design (ALIGNER_MEMBER_DESIGN.md Addendum (b) / §8): a standing
# germline/somatic variant near a splice site can FABRICATE a "novel" junction
# or be mis-charged as a mismatch the DP "repairs" with a spurious micro-indel —
# both directly INFLATE the de-novo (esp. non-canonical) discovery FDR. C6 is a
# first-class discovery-FDR GUARD, so the benchmark must MEASURE the variant-
# induced FDR the flat-haploid-reference panel suffers.
#
# DISCRIMINATING (the incumbent is BELOW ceiling, verified empirically against
# minimap2 -ax splice -uf, 2026-06-26):
#   * a GT..AG-flanked DELETION variant >= ~40bp is called as a spurious INTRON
#     (N-op) instead of a deletion (D) -> a FABRICATED, variant-adjacent junction
#     FP (and the deletion is NOT scored position-exact). Shorter deletions
#     (<=30bp) stay a correct D (the control).
# NON-discriminating but kept as the SPECIFICITY control (minimap2 is robust):
#   * a SNP near or far from a real junction is correctly called as a mismatch
#     and does NOT shift/fabricate a junction. These rows let the scorer show the
#     variant-induced FP is SPECIFIC to the splice-mimicking-deletion context,
#     not "any variant near a junction" (so the C6 guard is targeted, not a blunt
#     abstain-near-every-variant rule).
#
# The metric the scorer already emits for this stratum: junction FP +
# ``fp_variant_adjacent`` (the false junction lands within variant_adjacency_bp
# of a truth variant), plus indel position-exact concordance on the deletion.
SPLICE_MIMIC_DEL_LENS = [25, 40, 60, 100]   # <=30 = control (stays D); >=40 = FDR driver


def gen_variant_stratum(reps: int, rng: random.Random, locus0: int = 400
                        ) -> Tuple[Dict[str, str], List[Tuple[str, str]], List[ReadTruth]]:
    """Build the C6 standing-variant stratum (see module note above).

    Three sub-cases, all with truth recorded BY CONSTRUCTION:

    * ``SPLICE_MIMIC_DEL`` — single-exon contig with a GT..AG-flanked block; the
      read DELETES the block (a deletion variant). Truth = the deletion (no
      junction). minimap2 fabricates a variant-adjacent junction FP for blocks
      >= ~40bp. zygosity/VAF varied across loci for the het/hom + non-Mendelian
      (aneuploid-A549) stratification the SPEC asks for.
    * ``SNP_NEAR_JUNC`` — spliced gene + a SNP 3bp into exon1 from the real donor.
      Truth = the real junction + a variant-adjacent SNP. Control: minimap2 keeps
      the junction and calls the SNP a mismatch.
    * ``SNP_DISTANT`` — spliced gene + a SNP mid-exon (~100bp from any junction).
      Truth = the real junction + a variant-DISTANT SNP. Control.
    """
    refs: Dict[str, str] = {}
    reads: List[Tuple[str, str]] = []
    truth: List[ReadTruth] = []
    li = locus0

    # ---- SPLICE_MIMIC_DEL: the variant-induced-FDR driver ----------------
    L = 200
    zyg_cycle = [("HET", 0.5), ("HOM", 1.0), ("HET", 0.33)]  # incl. non-Mendelian VAF
    for di, dlen in enumerate(SPLICE_MIMIC_DEL_LENS):
        chrom = f"var_del_{dlen:03d}"
        left = _rand_unique(L, rng)
        # the deleted block looks splice-like (GT..AG) so minimap2 is tempted to
        # call it an intron; the flanks are random so the deletion has no in-run
        # ambiguity slide (eq span == the edit).
        block = "GT" + _rand_unique(dlen - 4, rng) + "AG"
        right = _rand_unique(L, rng)
        ref = left + block + right
        refs[chrom] = ref
        var_pos = L                       # 0-based start of the deleted block
        read_seq = left + right           # the read lacks the block
        zyg, vaf = zyg_cycle[di % len(zyg_cycle)]
        split = _split_for(li)
        li += 1
        for i in range(reps):
            rid = f"{chrom}_r{i:03d}"
            reads.append((rid, read_seq))
            # truth: a unique-context deletion variant; NO junction.
            indel = make_unique_indel(var_pos, dlen, IndelKind.DEL)
            indel.context = "VARIANT"
            var = VariantTruth(
                pos=var_pos, ref_allele=block, alt_allele="",
                zygosity=zyg, vaf=vaf, dist_to_junction=None)
            truth.append(ReadTruth(
                read_id=rid, true_locus=chrom, true_transcript=chrom,
                chrom=chrom, strand="+", genome_start=0, genome_end=len(ref),
                true_cigar=f"{var_pos}M{dlen}D{len(ref) - var_pos - dlen}M",
                indels=[indel], variants=[var], stratum="VARIANT",
                split=split, coverage=reps))

    # ---- SNP controls (minimap2 robust) on a spliced gene ----------------
    def _spliced_snp(tag: str, snp_offset_from_donor: int, dist_label):
        nonlocal li
        chrom = f"var_snp_{tag}"
        e1 = _rand_unique(200, rng)
        intron = "GT" + _rand_unique(196, rng) + "AG"
        e2 = _rand_unique(200, rng)
        contig = e1 + intron + e2
        refs[chrom] = contig
        i_start, i_end = 200, 200 + len(intron)
        model = TranscriptModel(
            name=tag, chrom=chrom, strand="+",
            exons=[Exon(0, 200), Exon(i_end, len(contig))], genome_seq=contig)
        jt = model.junction_truths()      # NNC (no annotation supplied)
        spliced = model.spliced_transcript()
        # transcript coord of the SNP: in exon1, snp_offset_from_donor bp 5' of
        # the donor (so genome pos = 200 - off; transcript pos identical in exon1)
        tpos = 200 - snp_offset_from_donor
        gpos = tpos
        ref_base = spliced[tpos]
        alt_base = rng.choice([b for b in BASES if b != ref_base])
        read_seq = spliced[:tpos] + alt_base + spliced[tpos + 1:]
        dist = abs(gpos - i_start) if dist_label == "near" else None
        split = _split_for(li)
        li += 1
        for i in range(reps):
            rid = f"{chrom}_r{i:03d}"
            reads.append((rid, read_seq))
            var = VariantTruth(
                pos=gpos, ref_allele=ref_base, alt_allele=alt_base,
                zygosity="HET", vaf=0.5,
                dist_to_junction=dist)
            truth.append(ReadTruth(
                read_id=rid, true_locus=tag, true_transcript=tag, chrom=chrom,
                strand="+", genome_start=model.genome_start,
                genome_end=model.genome_end, true_cigar=model.fulllength_cigar(),
                junctions=jt, variants=[var], true_cpa=model.genome_end - 1,
                stratum="VARIANT", split=split, coverage=reps))

    _spliced_snp("near", snp_offset_from_donor=3, dist_label="near")
    _spliced_snp("distant", snp_offset_from_donor=100, dist_label="far")
    return refs, reads, truth


# ---------------------------------------------------------------------------
# Top-level corpus generator
# ---------------------------------------------------------------------------
def generate_corpus(out_dir: str, reps: int = 120, seed: int = 7,
                    include_junction: bool = True,
                    include_variant: bool = True) -> Dict[str, str]:
    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(seed)
    refs: Dict[str, str] = {}
    reads: List[Tuple[str, str]] = []
    truth: List[ReadTruth] = []

    r1, rd1, t1 = gen_hp_stratum(reps, rng)
    rh, rdh, th = gen_hp_hard_stratum(reps, rng)
    r2, rd2, t2 = gen_str_stratum(reps, rng)
    refs.update(r1); refs.update(rh); refs.update(r2)
    reads += rd1 + rdh + rd2
    truth += t1 + th + t2

    if include_junction:
        jc, jseq, _model, jrd, jt = gen_junction_ambiguity_locus(rng)
        refs[jc] = jseq
        reads += jrd
        truth += jt
        rj, rdj, tj = gen_junction_class_loci(rng)
        refs.update(rj)
        reads += rdj
        truth += tj

    if include_variant:
        rv, rdv, tv = gen_variant_stratum(reps, rng)
        refs.update(rv)
        reads += rdv
        truth += tv

    ref_fa = os.path.join(out_dir, "ref.fa")
    reads_fq = os.path.join(out_dir, "reads.fastq")
    truth_tsv = os.path.join(out_dir, "truth.tsv")
    with open(ref_fa, "w") as fh:
        for n, s in refs.items():
            fh.write(f">{n}\n{s}\n")
    with open(reads_fq, "w") as fh:
        for rid, s in reads:
            fh.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")
    n = write_truth_table(truth, truth_tsv)
    return {"ref_fa": ref_fa, "reads_fastq": reads_fq, "truth_tsv": truth_tsv,
            "n_reads": str(len(reads)), "n_truth": str(n), "n_refs": str(len(refs))}


def main():
    ap = argparse.ArgumentParser(description="Tier-1 controlled benchmark generator")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=120, help=">=100 to clear min_count floor")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    info = generate_corpus(args.out, args.reps, args.seed)
    for k, v in info.items():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
