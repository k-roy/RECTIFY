#!/usr/bin/env python3
"""Deliverable B — C1 motif-aware realign of 3 real A549 SUPPORTED_NONCANONICAL junctions.

Per junction there are TWO candidate placements (from $W/rederive_111.json):
  * canon  = the nearby canonical GT-AG / GC-AG motif (where GMAP, a minority, places it)
  * nonc   = the dominant supported NON-canonical placement (where the bulk of reads/aligners
             AND accurate short reads place it), 1-4 bp away.

Question: does the empirical del-law (HP-context edit distance, junction_refiner's scorer)
SNAP each supporting long read to the canonical motif, or HOLD the non-canonical placement?

Mechanism = the SHIPPED del-law scorer `_hp_edit_distance(read_query, spliced_ref, penalty_table)`
(NOT align_exon_block_global, whose real signature has no penalty_table arg). For each read we
build the read's contiguous junction-crossing query (the read has NO intron, so exon1-tail and
exon2-head bases are contiguous in query_sequence) anchored at two genome positions BOTH
hypotheses agree are exonic: [donor_min - W, donor_min) on exon1 and [acc_max, acc_max + W) on
exon2. We score that query against each hypothesis' spliced reference and compare cost.

  margin = dist(read, ref_nonc) - dist(read, ref_canon)
  margin > 0  -> canonical cheaper -> SNAP (this read)
  margin < 0  -> non-canonical cheaper -> HOLD (this read)
  margin == 0 -> TIE

CONTROLS (advisor):
  * Score under penalty_table=LAW (yeast del-law) AND penalty_table=None (plain unit-cost edit
    distance). If the verdict is identical under both, the yeast-table caveat is moot. Divergence
    = exactly where yeast-on-human bites; surfaced per junction.
  * Report the per-read margin distribution, not just a binary tally. Apply the shipped
    _CANONICAL_HP_PRIOR=0.5 (canonical discount) AS A SEPARATE PASS and report how many reads it
    flips — lead with the NO-prior (raw evidence) result.
  * Per-read sanity: classify each read by where GMAP put its OWN N-op (canon/nonc/other);
    reads GMAP placed canonically vs non-canonically should give the same verdict.

Author: Kevin R. Roy (Deliverable-B agent). Yeast table on human reads is FLAGGED; in-silico is
SUGGESTIVE, RT-PCR/Sanger is the gold standard.
"""
from __future__ import annotations
import argparse, json, os, sys
from collections import Counter, defaultdict

import pysam
from rectify.core.splice.hp_penalty import _hp_edit_distance, HpPenaltyTable

# 0-based, intron = [start, end); donor = fetch(chrom,start,start+2), acceptor = fetch(chrom,end-2,end)
# All 6 motifs verified to reproduce under this convention against GRCh38_gencode_v44.fasta.
JUNCS = {
    "SQSTM1":  dict(chrom="chr5", canon=(179824400, 179832205), nonc=(179824404, 179832209),
                    motif_canon="GT..AG", motif_nonc="GT..GA", strand="+"),
    "TMED9":   dict(chrom="chr5", canon=(177592500, 177593474), nonc=(177592499, 177593473),
                    motif_canon="GC..AG", motif_nonc="CG..CA", strand="+"),
    "SLC35A4": dict(chrom="chr5", canon=(140564954, 140565547), nonc=(140564957, 140565550),
                    motif_canon="GT..AG", motif_nonc="AG..CA", strand="+"),
}
PRIOR = 0.5  # _CANONICAL_HP_PRIOR


TOL = 8  # bp tolerance for matching a read's N-op to my junction (offset<=4 + aligner wobble)


def matching_Nop(read, J):
    """Return the read's N-op interval (dS,dE) that splices at MY junction (within TOL of
    canon OR nonc placement), else None. Excludes reads that splice at other introns or
    retain the intron."""
    cd, ca = J["canon"]
    nd, na = J["nonc"]
    ref = read.reference_start
    for op, ln in read.cigartuples:
        if op == 3:  # N
            dS, dE = ref, ref + ln
            donor_ok = min(abs(dS - cd), abs(dS - nd)) <= TOL
            acc_ok = min(abs(dE - ca), abs(dE - na)) <= TOL
            if donor_ok and acc_ok:
                return (dS, dE)
        if op in (0, 2, 3, 7, 8):  # M,D,N,=,X consume ref
            ref += ln
    return None


def classify_placement(nop, J):
    """Which placement does this read's N-op match? canon / nonc / other (by nearest donor)."""
    dS, dE = nop
    cd, ca = J["canon"]
    nd, na = J["nonc"]
    dc = abs(dS - cd) + abs(dE - ca)
    dn = abs(dS - nd) + abs(dE - na)
    if dc == dn:
        return "tie"
    return "canon" if dc < dn else "nonc"


def nearest_qpos(pairs, target):
    """pairs: sorted list of (qpos,rpos) M-positions. Return qpos of the pair whose rpos is
    nearest target, plus the distance. None if empty."""
    best = None
    bestd = None
    for qpos, rpos in pairs:
        d = abs(rpos - target)
        if bestd is None or d < bestd:
            bestd, best = d, qpos
    return best, bestd


def harvest_reads(bam, chrom, donor_min, acc_max, a_up, a_dn, W, J, source):
    """Realign-ready records for reads that splice at MY junction (N-op within TOL of canon
    or nonc). Read query = contiguous slice anchored at nearest-match to the two shared-exonic
    anchors (a_up, a_dn-1); the read has no intron so this slice crosses the splice contiguously."""
    out = []
    for read in bam.fetch(chrom, donor_min, acc_max):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None:
            continue
        nop = matching_Nop(read, J)
        if nop is None:
            continue
        pairs = [(q, r) for q, r in read.get_aligned_pairs(matches_only=True)]
        if not pairs:
            continue
        qs, ds = nearest_qpos(pairs, a_up)
        qe, de = nearest_qpos(pairs, a_dn - 1)
        # require the read to actually reach near both exonic anchors (else it doesn't flank)
        if qs is None or qe is None or ds > 4 or de > 4 or qe <= qs:
            continue
        rq = read.query_sequence[qs:qe + 1].upper()
        if not rq or len(rq) > 6 * W:  # guard against intron-retention slipping through
            continue
        out.append({"qname": read.query_name, "rq": rq, "rqlen": len(rq),
                    "placement": classify_placement(nop, J), "nop": list(nop),
                    "source": source})
    return out


def score_records(records, ref_canon, ref_nonc, LAW):
    for rec in records:
        rq = rec["rq"]
        for label, tab in (("law", LAW), ("plain", None)):
            d_can = _hp_edit_distance(rq, ref_canon, penalty_table=tab)
            d_non = _hp_edit_distance(rq, ref_nonc, penalty_table=tab)
            rec[label] = {"d_can": round(d_can, 4), "d_non": round(d_non, 4),
                          "margin": round(d_non - d_can, 4)}
    return records


def run(genome, bam_paths, penalty_tsv, W, out_json):
    fa = pysam.FastaFile(genome)
    LAW = HpPenaltyTable.from_tsv(penalty_tsv)
    results = {"window": W, "penalty_tsv": penalty_tsv, "bams": bam_paths, "junctions": {}}
    src_of = lambda p: os.path.basename(p).replace("a549_chr5_trimmed.", "").replace(".bam", "")

    for gene, J in JUNCS.items():
        chrom = J["chrom"]
        cd, ca = J["canon"]
        nd, na = J["nonc"]
        donor_min = min(cd, nd)
        acc_max = max(ca, na)
        a_up = donor_min - W           # left shared-exonic anchor (genome, inclusive)
        a_dn = acc_max + W             # right shared-exonic anchor (genome, exclusive)

        def ref_for(donor, acc):
            return (fa.fetch(chrom, a_up, donor) + fa.fetch(chrom, acc, a_dn)).upper()

        ref_canon = ref_for(cd, ca)
        ref_nonc = ref_for(nd, na)

        # harvest from every source BAM, keep per-source and a pooled-unique (by qname) set
        all_recs = []
        seen = {}
        for bp in bam_paths:
            bam = pysam.AlignmentFile(bp, "rb")
            src = src_of(bp)
            recs = harvest_reads(bam, chrom, donor_min, acc_max, a_up, a_dn, W, J, src)
            all_recs.extend(recs)
            for r in recs:
                seen.setdefault(r["qname"], r)  # first occurrence = pooled-unique
            bam.close()
        per_read = score_records(list(seen.values()), ref_canon, ref_nonc, LAW)
        # also score per-source for the cross-source consistency check
        by_source = defaultdict(list)
        for r in score_records(all_recs, ref_canon, ref_nonc, LAW):
            by_source[r["source"]].append(r)

        # aggregate
        agg = {}
        for label in ("law", "plain"):
            tally = Counter()
            tally_prior = Counter()
            place_split = defaultdict(Counter)
            margins = []
            for r in per_read:
                m = r[label]["margin"]
                margins.append(m)
                v = "SNAP" if m > 0 else ("HOLD" if m < 0 else "TIE")
                tally[v] += 1
                place_split[r["placement"]][v] += 1
                # canonical prior: subtract 0.5 from canonical cost => add 0.5 to margin
                mp = m + PRIOR
                vp = "SNAP" if mp > 0 else ("HOLD" if mp < 0 else "TIE")
                tally_prior[vp] += 1
            margins.sort()
            agg[label] = {
                "n": len(per_read),
                "tally_noprior": dict(tally),
                "tally_prior0.5": dict(tally_prior),
                "margin_min": round(margins[0], 3) if margins else None,
                "margin_med": round(margins[len(margins)//2], 3) if margins else None,
                "margin_max": round(margins[-1], 3) if margins else None,
                "margin_mean": round(sum(margins)/len(margins), 4) if margins else None,
                "placement_split": {k: dict(v) for k, v in place_split.items()},
            }
        # per-source verdict (cross-source consistency: harvesting-aligner bias check)
        per_source = {}
        for src, recs in sorted(by_source.items()):
            t = Counter()
            for r in recs:
                m = r["law"]["margin"]
                t["SNAP" if m > 0 else ("HOLD" if m < 0 else "TIE")] += 1
            per_source[src] = {"n": len(recs), "law_noprior": dict(t)}

        results["junctions"][gene] = {
            "chrom": chrom, "canon": J["canon"], "nonc": J["nonc"],
            "motif_canon": J["motif_canon"], "motif_nonc": J["motif_nonc"],
            "anchor": [a_up, a_dn], "ref_canon": ref_canon, "ref_nonc": ref_nonc,
            "n_reads_pooled_unique": len(per_read), "agg": agg,
            "per_source": per_source, "per_read": per_read,
        }
        # console summary
        a = results["junctions"][gene]["agg"]
        print(f"\n=== {gene} {chrom}:{J['canon']}(canon {J['motif_canon']}) vs "
              f"{J['nonc']}(nonc {J['motif_nonc']})  n_pooled={len(per_read)} ===")
        for label in ("law", "plain"):
            print(f"  [{label:5s}] no-prior {a[label]['tally_noprior']}  "
                  f"+0.5prior {a[label]['tally_prior0.5']}  "
                  f"margin(min/med/max)={a[label]['margin_min']}/{a[label]['margin_med']}/{a[label]['margin_max']}  "
                  f"place={a[label]['placement_split']}")
        print(f"  per-source(law,no-prior): " +
              "  ".join(f"{s}:n{v['n']}{v['law_noprior']}" for s, v in per_source.items()))

    with open(out_json, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"\n[written] {out_json}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bams", required=True, help="comma-separated per-aligner BAM paths")
    ap.add_argument("--penalty-tsv", required=True)
    ap.add_argument("--window", type=int, default=15)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    bam_paths = [b for b in args.bams.split(",") if b]
    run(args.genome, bam_paths, args.penalty_tsv, args.window, args.out)


if __name__ == "__main__":
    main()
