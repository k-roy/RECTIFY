#!/usr/bin/env python
"""Round-2 cDNA Phase-0 — Step F.0: lift-over round-trip validation on REAL reads.

For reads that minimap2 spliced cleanly along a MANE multi-exon gene: project the genome
alignment into cDNA-local coords, lift it back, and assert the lifted CIGAR + reference_start
+ query reproduce the read's ORIGINAL genome alignment, and the orientation invariant holds.
This exercises N-insertion / merge / minus-strand revcomp+reversal on real multi-junction,
indel-bearing CIGARs on both strands — which the synthetic unit tests cannot.

MUST pass (both strands) before any GO/NO-GO verdict is trusted.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
from collections import defaultdict

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover import lift, assert_orientation_invariant, assert_query_length, N  # noqa: E402
from cdna_models import load_models, project_genome_to_cdna, norm_cigar  # noqa: E402

_MATCH = frozenset({0, 7, 8})       # M/=/X
_REF = frozenset({0, 2, 3, 7, 8})   # ref-consuming
_QRY = frozenset({0, 1, 4, 7, 8})   # query-consuming


def matched_pairs(cigartuples, ref_start):
    """Semantic (query_index, ref_pos) pairs over M/=/X — invariant to how I/N order at a
    junction is represented. Two CIGARs with equal matched-pair sets describe the SAME
    alignment of query bases to genome positions (the round-trip correctness criterion)."""
    pairs = []
    q = 0
    r = ref_start
    for op, ln in cigartuples:
        if op in _MATCH:
            for k in range(ln):
                pairs.append((q + k, r + k))
        if op in _REF:
            r += ln
        if op in _QRY:
            q += ln
    return pairs


def model_span(model):
    exon = [b for b in model.blocks if not b.is_pad]
    return model.chrom, min(b.g_start for b in exon), max(b.g_end for b in exon)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blockmaps", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bam", required=True, help="minimap2 genome BAM")
    ap.add_argument("--per-strand", type=int, default=10, help="min clean reads required per strand")
    ap.add_argument("--max-test", type=int, default=400, help="cap reads tested per gene")
    args = ap.parse_args()

    models = load_models(args.blockmaps)
    # only MANE models for the round-trip ground truth
    mane = {cid: m for cid, m in models.items() if "__MANE" in cid}
    spans = {cid: model_span(m) for cid, m in mane.items()}
    fa = pysam.FastaFile(args.genome)
    genome = {c: fa.fetch(c) for c in {m.chrom for m in mane.values()}}

    bam = pysam.AlignmentFile(args.bam, "rb")
    passes = defaultdict(int)   # strand -> n pass
    fails = []
    tested = 0
    examples = []

    for cid, m in mane.items():
        chrom, lo, hi = spans[cid]
        seen = 0
        for read in bam.fetch(chrom, lo, hi):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.cigartuples is None or not any(op == N for op, _ in read.cigartuples):
                continue
            if read.query_sequence is None:
                continue
            proj = project_genome_to_cdna(m, read.cigartuples, read.reference_start, read.query_sequence)
            if proj is None:
                continue
            cdna_cigar, cdna_ref_start, cdna_query = proj
            seen += 1
            if seen > args.max_test:
                break
            tested += 1
            try:
                lifted = lift(m, cdna_cigar, cdna_ref_start, cdna_query)
            except Exception as e:
                fails.append((read.query_name, cid, f"lift raised {e!r}"))
                continue
            # semantic comparison: same query-base -> genome-pos mapping (I/N order-invariant)
            ok_cigar = (matched_pairs(lifted.cigartuples, lifted.reference_start)
                        == matched_pairs(read.cigartuples, read.reference_start))
            ok_start = (lifted.reference_start == read.reference_start)
            ok_query = (lifted.query_sequence == read.query_sequence)
            f_ok, l_ok = assert_orientation_invariant(lifted, genome)
            ok_orient = f_ok and l_ok
            try:
                assert_query_length(lifted, len(read.query_sequence))
                ok_qlen = True
            except AssertionError:
                ok_qlen = False
            if ok_cigar and ok_start and ok_query and ok_orient and ok_qlen:
                passes[m.strand] += 1
            else:
                reason = []
                if not ok_cigar: reason.append(f"matched-pairs differ ({norm_cigar(lifted.cigartuples)[:5]} vs {norm_cigar(read.cigartuples)[:5]})")
                if not ok_start: reason.append(f"start {lifted.reference_start}!={read.reference_start}")
                if not ok_query: reason.append("query")
                if not ok_orient: reason.append(f"orient(f={f_ok},l={l_ok})")
                if not ok_qlen: reason.append("qlen")
                fails.append((read.query_name, cid, "; ".join(reason)))

    print(f"[F.0] tested {tested} clean spliced reads")
    print(f"[F.0] PASS  + strand : {passes['+']}")
    print(f"[F.0] PASS  - strand : {passes['-']}")
    print(f"[F.0] FAIL           : {len(fails)}")
    for fn in fails[:25]:
        print(f"   FAIL {fn[0]} [{fn[1]}] : {fn[2]}")

    # The HARD correctness gate is ZERO lift mismatches (a single fail = a real lift bug).
    # Both strands must be exercised (>=1 clean read each) so the minus path is actually tested.
    # Low per-strand counts with zero fails is a coverage note, not a correctness failure.
    no_lift_bug = (len(fails) == 0)
    both_strands = (passes['+'] >= 1 and passes['-'] >= 1)
    well_powered = (passes['+'] >= args.per_strand and passes['-'] >= args.per_strand)
    ok = no_lift_bug and both_strands
    if ok and not well_powered:
        print(f"[F.0] NOTE: under-powered ({passes['+']}+/{passes['-']}- < {args.per_strand}); "
              f"correctness OK (0 fails) but few exact-MANE reads.")
    if not both_strands:
        print(f"[F.0] WARN: a strand had 0 clean reads tested — minus-strand lift not exercised.")
    print(f"\n[F.0] {'PASS — lift-over reproduces real genome alignments (0 mismatches)' if ok else 'FAIL — lift-over disagrees with real genome alignments; fix liftover.py before any verdict'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
