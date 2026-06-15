#!/usr/bin/env python
"""Debug F.0 failures: for given read ids, dump original genome CIGAR, the projected
cDNA-local alignment, the lifted-back CIGAR, and the first op where they diverge."""
import argparse, sys
from pathlib import Path
import pysam
sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover import lift, assert_orientation_invariant
from cdna_models import load_models, project_genome_to_cdna, norm_cigar, g2t

OPN = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 7: "=", 8: "X"}
def fmt(cig): return "".join(f"{l}{OPN.get(o,o)}" for o, l in cig)

def model_span(m):
    ex = [b for b in m.blocks if not b.is_pad]
    return m.chrom, min(b.g_start for b in ex), max(b.g_end for b in ex)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blockmaps", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reads", nargs="+", required=True)
    args = ap.parse_args()
    models = load_models(args.blockmaps)
    mane = {c: m for c, m in models.items() if "__MANE" in c}
    spans = {c: model_span(m) for c, m in mane.items()}
    fa = pysam.FastaFile(args.genome)
    genome = {c: fa.fetch(c) for c in {m.chrom for m in mane.values()}}
    want = set(args.reads)
    bam = pysam.AlignmentFile(args.bam, "rb")
    found = set()
    for cid, m in mane.items():
        chrom, lo, hi = spans[cid]
        for read in bam.fetch(chrom, lo, hi):
            if read.query_name not in want or read.is_secondary or read.is_supplementary:
                continue
            if read.query_name in found:
                continue
            proj = project_genome_to_cdna(m, read.cigartuples, read.reference_start, read.query_sequence)
            print(f"\n=== {read.query_name} on {cid} strand={m.strand} is_rev={read.is_reverse} ===")
            print(f" orig start={read.reference_start} cigar={fmt(read.cigartuples)}")
            if proj is None:
                print("  projection=None (matched base off-cDNA)"); found.add(read.query_name); continue
            cc, crs, cq = proj
            print(f" proj  cdna_ref_start={crs} cdna_cigar={fmt(cc)} qlen={len(cq)}")
            try:
                lifted = lift(m, cc, crs, cq)
            except Exception as e:
                print(f"  lift raised {e!r}"); found.add(read.query_name); continue
            o = norm_cigar(read.cigartuples); g = norm_cigar(lifted.cigartuples)
            print(f" lift  start={lifted.reference_start} cigar={fmt(g)}")
            print(f" match: start={lifted.reference_start==read.reference_start} cigar={o==g} "
                  f"query={lifted.query_sequence==read.query_sequence}")
            if o != g:
                for i, (a, b) in enumerate(zip(o, g)):
                    if a != b:
                        print(f"   first cigar divergence at op {i}: orig {a} vs lift {b}")
                        print(f"     orig[{max(0,i-2)}:{i+3}]={o[max(0,i-2):i+3]}")
                        print(f"     lift[{max(0,i-2)}:{i+3}]={g[max(0,i-2):i+3]}")
                        break
                if len(o) != len(g):
                    print(f"   length differ: orig {len(o)} ops vs lift {len(g)} ops")
            f_ok, l_ok = assert_orientation_invariant(lifted, genome)
            print(f" orient first={f_ok} last={l_ok}")
            found.add(read.query_name)
    miss = want - found
    if miss:
        print(f"\nNOT FOUND in span: {miss}")

if __name__ == "__main__":
    main()
