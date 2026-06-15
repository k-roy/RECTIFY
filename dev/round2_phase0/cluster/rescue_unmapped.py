#!/usr/bin/env python
"""Selection-bias check scorer: among GENOME-UNMAPPED reads, how many does the cDNA round
place with a confident, gate-passing alignment that the fresh genome panel still can't?

A "rescue" = a read that (a) maps to a locus cDNA covering >=80% of the read at >=85%
identity, (b) whose lifted genome alignment passes the k>=min_anchor junction gate, and
(c) that the fresh genome minimap2 leaves UNMAPPED or maps with <50% of the read aligned.
Spliced rescues (lifted n_junctions>=1) are the ones that actually exercise the hypothesis.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
from collections import defaultdict
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover import lift, assert_orientation_invariant, N, S, M, I, D, EQ, X  # noqa
from cdna_models import load_models  # noqa
from rectify.core.consensus.corrected_consensus import (  # noqa
    _cigar_hp_edit_distance, _cigar_min_junction_anchor, _NO_JUNCTION_ANCHOR)
from rectify.core.splice.hp_penalty import HpPenaltyTable  # noqa

_QRY = frozenset({M, I, S, EQ, X})
_ALN = frozenset({M, I, D, EQ, X})


def aligned_query_frac(read):
    if read.cigartuples is None:
        return 0.0
    qlen = sum(l for o, l in read.cigartuples if o in _QRY)
    if qlen == 0:
        return 0.0
    aln = sum(l for o, l in read.cigartuples if o in (M, EQ, X, I))
    return aln / qlen


def identity(read):
    try:
        nm = read.get_tag("NM")
    except KeyError:
        return None
    block = sum(l for o, l in read.cigartuples if o in _ALN)
    return (1 - nm / block) if block else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blockmaps", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--cdna-bam", required=True)
    ap.add_argument("--genome-bam", required=True)
    ap.add_argument("--penalty-scores", required=True)
    ap.add_argument("--str-penalty", default=None)
    ap.add_argument("--min-anchor", type=int, default=10)
    ap.add_argument("--min-frac", type=float, default=0.80)
    ap.add_argument("--min-ident", type=float, default=0.85)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    models = load_models(args.blockmaps)
    fa = pysam.FastaFile(args.genome)
    genome = {c: fa.fetch(c) for c in fa.references}
    ptab = HpPenaltyTable.from_tsv(args.penalty_scores, str_path=args.str_penalty)

    # genome placement lookup for the cDNA-hit reads
    gstat = {}  # rid -> ('unmapped'|frac)
    gb = pysam.AlignmentFile(args.genome_bam, "rb")
    for r in gb.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:
            continue
        rid = r.query_name
        if r.is_unmapped:
            gstat.setdefault(rid, ("unmapped", 0.0))
        else:
            f = aligned_query_frac(r)
            # keep the best (highest aligned frac) genome placement
            prev = gstat.get(rid)
            if prev is None or (prev[0] != "mapped") or f > prev[1]:
                gstat[rid] = ("mapped", f)

    n_cdna = 0
    good = 0
    rescues = []
    cb = pysam.AlignmentFile(args.cdna_bam, "rb")
    for read in cb.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_reverse:
            continue
        cid = read.reference_name
        if cid not in models or read.query_sequence is None:
            continue
        n_cdna += 1
        frac = aligned_query_frac(read)
        idt = identity(read) or 0.0
        if frac < args.min_frac or idt < args.min_ident:
            continue
        m = models[cid]
        try:
            lifted = lift(m, read.cigartuples, read.reference_start, read.query_sequence)
        except Exception:
            continue
        f_ok, l_ok = assert_orientation_invariant(lifted, genome)
        if not (f_ok and l_ok):
            continue
        anchor = _cigar_min_junction_anchor(lifted, genome)
        njunc = sum(1 for o, l in lifted.cigartuples if o == N)
        # gate: junction-spanning lifts must pass the perfect-anchor gate
        if njunc >= 1 and anchor < args.min_anchor:
            continue
        good += 1
        gs = gstat.get(read.query_name, ("absent", 0.0))
        genome_fails = (gs[0] in ("unmapped", "absent")) or (gs[0] == "mapped" and gs[1] < 0.5)
        if genome_fails:
            rescues.append(dict(rid=read.query_name, cid=cid, gene=cid.split("__")[0],
                                njunc=njunc, anchor=(None if anchor >= _NO_JUNCTION_ANCHOR else anchor),
                                cdna_frac=round(frac, 3), cdna_ident=round(idt, 3),
                                genome=gs[0], genome_frac=round(gs[1], 3),
                                hp_ed=round(_cigar_hp_edit_distance(lifted, genome, ptab), 1)))

    spliced_rescues = [r for r in rescues if r["njunc"] >= 1]
    # write
    with open(out / "rescues.tsv", "w") as fh:
        cols = ["rid", "gene", "cid", "njunc", "anchor", "cdna_frac", "cdna_ident",
                "genome", "genome_frac", "hp_ed"]
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rescues, key=lambda x: -x["njunc"]):
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    from collections import Counter
    txt = (
        f"[bias-check: cDNA rescue of GENOME-UNMAPPED reads]\n"
        f"  unmapped reads mapped to a cDNA (primary,fwd): {n_cdna}\n"
        f"  confident gate-passing cDNA placements (frac>={args.min_frac}, ident>={args.min_ident}): {good}\n"
        f"  RESCUES (good cDNA placement AND genome unmapped/<50%%): {len(rescues)}\n"
        f"  of which SPLICED (lifted junction>=1): {len(spliced_rescues)}\n"
        f"  spliced-rescue genes: {dict(Counter(r['gene'] for r in spliced_rescues))}\n"
        f"  spliced-rescue junction counts: {dict(Counter(r['njunc'] for r in spliced_rescues))}\n"
    )
    print(txt)
    (out / "rescue_summary.txt").write_text(txt)
    with open(out / "rescue_summary.json", "w") as fh:
        json.dump(dict(n_cdna_mapped=n_cdna, n_good=good, n_rescue=len(rescues),
                       n_spliced_rescue=len(spliced_rescues)), fh, indent=1)


if __name__ == "__main__":
    main()
