"""Local test of cdna_models.project_genome_to_cdna (the F.0 round-trip path) on synthetic
multi-junction reads, both strands. project(lift(x)) must reproduce the genome alignment.
Run with base miniconda python on M1."""
import random, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # for liftover if needed
from liftover import Block, CdnaModel, lift, revcomp
from cdna_models import project_genome_to_cdna, norm_cigar, g2t

M, I, D, N, S, EQ, X = 0, 1, 2, 3, 4, 7, 8


def genome(seed=1, length=600):
    random.seed(seed)
    return {"T1": "".join(random.choice("ACGT") for _ in range(length))}


def plus4(g):
    blocks = [Block(0, 100, 120), Block(20, 200, 209), Block(29, 300, 320), Block(49, 400, 420)]
    cdna = CdnaModel("plus4", "T1", "+", blocks, 69)
    cdna.validate(g)
    seq = g["T1"]
    cseq = seq[100:120] + seq[200:209] + seq[300:320] + seq[400:420]
    return cdna, cseq


def minus4(g):
    blocks = [Block(0, 400, 420), Block(20, 300, 320), Block(40, 200, 209), Block(49, 100, 120)]
    cdna = CdnaModel("minus4", "T1", "-", blocks, 69)
    cdna.validate(g)
    seq = g["T1"]
    cseq = revcomp(seq[400:420]) + revcomp(seq[300:320]) + revcomp(seq[200:209]) + revcomp(seq[100:120])
    return cdna, cseq


def roundtrip(model, genome_cigar, ref_start, genome_query):
    """genome alignment -> project to cdna-local -> lift back -> must equal input."""
    proj = project_genome_to_cdna(model, genome_cigar, ref_start, genome_query)
    assert proj is not None, "projection returned None (read not clean on model?)"
    cdna_cigar, cdna_ref_start, cdna_query = proj
    lifted = lift(model, cdna_cigar, cdna_ref_start, cdna_query)
    assert norm_cigar(lifted.cigartuples) == norm_cigar(genome_cigar), \
        f"cigar {norm_cigar(lifted.cigartuples)} != {norm_cigar(genome_cigar)}"
    assert lifted.reference_start == ref_start, f"start {lifted.reference_start} != {ref_start}"
    assert lifted.query_sequence == genome_query, "query mismatch"
    return True


def main():
    g = genome()
    seq = g["T1"]
    fails = 0
    cases = []

    # --- plus, clean multi-junction (from lift of a clean read) ---
    cdna, cseq = plus4(g)
    lp = lift(cdna, [(M, 69)], 0, cseq)
    cases.append(("plus_clean", cdna, lp.cigartuples, lp.reference_start, lp.query_sequence))

    # --- minus, clean multi-junction ---
    cdnam, cseqm = minus4(g)
    lm = lift(cdnam, [(M, 69)], 0, cseqm)
    cases.append(("minus_clean", cdnam, lm.cigartuples, lm.reference_start, lm.query_sequence))

    # --- plus with an insertion at a junction (E6) ---
    read = cseq[:20] + "GGG" + cseq[20:]
    lp2 = lift(cdna, [(M, 20), (I, 3), (M, 49)], 0, read)
    cases.append(("plus_ins_junction", cdna, lp2.cigartuples, lp2.reference_start, lp2.query_sequence))

    # --- plus with deletion at a junction (E5) ---
    read = cseq[:18] + cseq[20:]
    lp3 = lift(cdna, [(M, 18), (D, 2), (M, 49)], 0, read)
    cases.append(("plus_del_junction", cdna, lp3.cigartuples, lp3.reference_start, lp3.query_sequence))

    # --- plus with leading + trailing soft-clip (read overhang) ---
    read = "TTTT" + cseq + "AAA"
    lp4 = lift(cdna, [(S, 4), (M, 69), (S, 3)], 0, read)
    cases.append(("plus_softclip_both", cdna, lp4.cigartuples, lp4.reference_start, lp4.query_sequence))

    # --- minus with insertion at a junction (E6 x E11) ---
    readm = cseqm[:20] + "CC" + cseqm[20:]
    lm2 = lift(cdnam, [(M, 20), (I, 2), (M, 49)], 0, readm)
    cases.append(("minus_ins_junction", cdnam, lm2.cigartuples, lm2.reference_start, lm2.query_sequence))

    # --- minus with soft-clips both ends ---
    readm = "GG" + cseqm + "TTTTT"
    lm3 = lift(cdnam, [(S, 2), (M, 69), (S, 5)], 0, readm)
    cases.append(("minus_softclip_both", cdnam, lm3.cigartuples, lm3.reference_start, lm3.query_sequence))

    for name, model, cig, rs, q in cases:
        try:
            roundtrip(model, cig, rs, q)
            print(f"PASS  {name}")
        except AssertionError as e:
            fails += 1
            print(f"FAIL  {name}: {e}")
    print(f"\n{len(cases) - fails}/{len(cases)} projection round-trips passed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
