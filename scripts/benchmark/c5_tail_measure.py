#!/usr/bin/env python3
"""C5 panel-failure-tail MEASUREMENT — the dep-commit gate (decoder-free).

C5's facet (NATIVE_ALIGNER_OVERVIEW, row C5) is a **FracMinHash containment
fallback localizer** for the panel-failure tail — and it is EXPLICITLY GATED
behind a *measured depletion trigger*: build nothing until a measurement proves a
recoverable tail EXISTS at realistic error.  This script is that measurement.  It
writes NO FracMinHash production code; Stage 2 only PROTOTYPES containment to test
the pre-set kill-gate.

THE LOAD-BEARING DISTINCTION (reconciling the director's lead with the OVERVIEW)
-------------------------------------------------------------------------------
The director's lead calls the JUNCTION_DISCOVERY herd ("all members snap, no
member produces truth") C5 territory.  But OVERVIEW lines 88-92 are explicit: the
discovery ceiling is at the **WINDOW** level — a read in the RIGHT window with a
mis-called junction inside it is recoverable by **realignment** (C3/refiner), and
"only reads with no acceptable window at all ... are the only place an independent
localizer earns its keep."  The C5 facet is specifically the **localizer**.

So an "all-misplaced" read is only **C5/FracMinHash-addressable** if it is misplaced
at the WINDOW level (wrong contig/locus, or unmapped).  We therefore split the tail:

  * empty-union          : NO panel member placed the read at all          -> C5
  * wrong-window         : placed, but in the WRONG locus (reference_name !=
                           truth chrom, or no overlap of the truth span)
      - PARALOG sub-case : near-identical copies -> containment matches BOTH,
                           cannot disambiguate -> C4-pooling territory, NOT C5
      - non-paralog      : a genuine window-localization failure             -> C5
  * right-window-wrong-internal : placed at the CORRECT locus but junction/indel
                           mis-called (the JUNCTION_DISCOVERY herd) -> C3/refiner,
                           NOT C5 (a localizer adds nothing to a right window)
  * correct              : right window, position-exact

  C5-addressable tail = empty-union  UNION  (wrong-window AND non-paralog)

Lumping the junction-herd into C5 would OVERSTATE the tail and drift toward a false
PROCEED.  Classified correctly, the baseline tail on the tiny Tier-1 contigs is
expected ~0 (Tier-2 placed 98-99.8%).

WHY minimap2-ALONE is the conservative UPPER BOUND (advisor)
------------------------------------------------------------
A single seed-chain member's union is the LARGEST possible panel-unplaced set;
adding members can only SHRINK it.  So if minimap2-alone's *classified* recoverable
tail is ~0, DEFER is robust without ever wiring the multi-aligner panel.  (The DP
arms in c3_headroom are NOT valid here: align_exon_block_global is handed the TRUTH
contig, so it can never be locus_incorrect and cannot measure window selection.)

ELEVATED ERROR (the trigger needs a guaranteed tail; spec lines 122-124)
------------------------------------------------------------------------
The tail "will NOT appear at low simulated error."  We inject a tunable hot-read
subpopulation with **bursty/clustered** errors (the spec's own ERROR-REALISM
verdict, line 283: uniform noise is unrealistic and is also the easiest for
minimap2 to seed through).  Injection is POST-HOC in THIS script — controlled.py
and the smoke stay byte-identical by construction (no shared RNG touched).  We only
need the ORIGIN WINDOW truth, which survives noising; we report tail size as a
FUNCTION of escalating hot-read error, not one number.

RECOVERABLE vs ZERO-EVIDENCE (the discriminator that decides REFUTE vs PROCEED)
------------------------------------------------------------------------------
A tail read RETAINS recoverable origin truth iff its (noised) sequence still shares
k-mer signal with its TRUE origin window (a containment localizer COULD find it).
A read noised into oblivion shares no k-mers -> zero-evidence -> never fabricate a
placement.  The injection is NOT tuned to be containment-friendly: bursty errors on
a hot subpopulation naturally yield BOTH recoverable (partly-clean) and
zero-evidence (destroyed) reads, and we report the honest split.

  trigger = fraction of the classified C5-addressable tail that is RECOVERABLE.

STAGE 2 (only if a recoverable tail exists): does FracMinHash *containment*
localize those reads to the CORRECT origin window better than a genome-scale
random-window null?  Candidate universe = all corpus contig windows + S288C-tiled
DECOY windows (so the null is not trivially good on a handful of contigs).

Usage:
  python scripts/benchmark/c5_tail_measure.py --out /tmp/c5 --reps 120 \
      --hot-frac 0.3 --rates 0,0.05,0.10,0.20,0.35,0.50
Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from rectify.core.benchmark.scorer import load_genome  # noqa: E402
from rectify.core.benchmark.truth_schema import read_truth_table  # noqa: E402
from scripts.benchmark.sim import controlled  # noqa: E402
from scripts.benchmark.c3_headroom import _read_is_position_exact  # noqa: E402
import pysam  # noqa: E402

_BASES = "ACGT"


# --------------------------------------------------------------------------
# FracMinHash-style k-mer sketch (PROTOTYPE — NOT production C5 code).
# FNV-1a 64-bit hash (deterministic across processes, unlike Python str hash).
# --------------------------------------------------------------------------
_FNV_OFF = 0xCBF29CE484222325
_FNV_PRIME = 0x100000001B3
_MASK64 = (1 << 64) - 1


def _fnv1a(s: str) -> int:
    h = _FNV_OFF
    for ch in s:
        h ^= ord(ch)
        h = (h * _FNV_PRIME) & _MASK64
    return h


def fracminhash(seq: str, k: int = 15, scale: int = 1) -> set:
    """Return the FracMinHash sketch: hashes of all k-mers whose hash falls in the
    bottom 1/scale of hash space.  scale=1 keeps ALL k-mers (max sensitivity, the
    right choice for the SHORT Tier-1 reads/windows; scale>1 is the production
    sub-sample).  Canonical k-mer (min of fwd/rev-comp) so strand doesn't matter."""
    seq = seq.upper()
    if len(seq) < k:
        return set()
    thresh = _MASK64 // scale if scale > 1 else _MASK64
    out = set()
    rc = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if any(c not in rc for c in kmer):
            continue
        rev = "".join(rc[c] for c in reversed(kmer))
        canon = kmer if kmer <= rev else rev
        h = _fnv1a(canon)
        if h <= thresh:
            out.add(h)
    return out


def containment(query_sketch: set, target_sketch: set) -> float:
    if not query_sketch:
        return 0.0
    return len(query_sketch & target_sketch) / len(query_sketch)


# --------------------------------------------------------------------------
# Bursty hot-read error injection (post-hoc; truth window preserved).
# --------------------------------------------------------------------------
def inject_bursty(seq: str, rate: float, rng: random.Random,
                  burst_len: int = 6, n_indel_frac: float = 0.4) -> str:
    """Inject ~rate per-base errors CLUSTERED in bursts (spec line 283: real ONT
    has ~2.83x excess sub-5bp gaps over geometric).  A fraction are indels (multi-
    base runs), the rest substitutions.  Bursts of length ~burst_len are dropped at
    random anchors until the realized error budget is met.  Returns a NEW string;
    the origin window in truth is unchanged (we never claim per-base indel truth on
    noised reads — only WINDOW correctness is asked of minimap2 here)."""
    if rate <= 0:
        return seq
    chars = list(seq)
    n_target = int(round(rate * len(seq)))
    realized = 0
    guard = 0
    while realized < n_target and guard < 10 * len(seq) + 50:
        guard += 1
        anchor = rng.randrange(len(chars))
        blen = max(1, int(rng.expovariate(1.0 / burst_len)))
        for j in range(blen):
            p = anchor + j
            if p >= len(chars):
                break
            roll = rng.random()
            if roll < n_indel_frac * 0.5:          # deletion
                chars[p] = ""
            elif roll < n_indel_frac:               # insertion
                chars[p] = chars[p] + rng.choice(_BASES)
            else:                                    # substitution
                alt = rng.choice(_BASES)
                while alt == chars[p][:1].upper():
                    alt = rng.choice(_BASES)
                chars[p] = alt
            realized += 1
            if realized >= n_target:
                break
    return "".join(chars)


def write_fastq(reads, path):
    with open(path, "w") as fh:
        for rid, s in reads:
            if not s:
                s = "N"
            fh.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")


def run_minimap2(ref_fa, reads_fq, out_bam):
    p = subprocess.run(
        ["minimap2", "-ax", "splice", "-uf", "--eqx", "-k", "14", "-t", "2",
         ref_fa, reads_fq], capture_output=True)
    if p.returncode:
        raise RuntimeError("minimap2 failed: " + p.stderr.decode()[:400])
    s = subprocess.run(["samtools", "sort", "-o", out_bam], input=p.stdout,
                       capture_output=True)
    if s.returncode:
        raise RuntimeError("samtools sort failed: " + s.stderr.decode()[:300])
    subprocess.run(["samtools", "index", out_bam], check=True)


def classify(bam_path, truth_map, genome):
    """Per-read WINDOW-level classification against truth.
    Returns {rid: class} where class in
      {correct, right_window_wrong_internal, wrong_window, empty_union}."""
    seen = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            rid = read.query_name
            rt = truth_map.get(rid)
            if rt is None:
                continue
            if read.is_unmapped:
                seen.setdefault(rid, "empty_union")
                continue
            # locus / window correctness
            locus_ok = (read.reference_name == rt.chrom and
                        read.reference_start < rt.genome_end and
                        read.reference_end > rt.genome_start)
            if not locus_ok:
                seen[rid] = "wrong_window"
                continue
            # right window: is the internal placement position-exact vs truth?
            gseq = genome.get(rt.chrom, "")
            try:
                exact = _read_is_position_exact(read, rt, genome)
            except Exception:
                exact = False
            # a primary right-window read overrides any earlier non-primary note
            prev = seen.get(rid)
            cls = "correct" if exact else "right_window_wrong_internal"
            if prev in (None, "empty_union"):
                seen[rid] = cls
            elif prev == "wrong_window":
                seen[rid] = cls  # a right-window primary wins
            else:
                # keep the better of the two
                if prev == "right_window_wrong_internal" and cls == "correct":
                    seen[rid] = cls
    # truth reads never seen at all -> empty_union
    for rid in truth_map:
        seen.setdefault(rid, "empty_union")
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/c5")
    ap.add_argument("--reps", type=int, default=120)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--noise-seed", type=int, default=101)
    ap.add_argument("--hot-frac", type=float, default=0.30,
                    help="fraction of reads that are 'hot' (get the elevated rate)")
    ap.add_argument("--rates", default="0,0.05,0.10,0.20,0.35,0.50",
                    help="comma list of hot-read per-base error rates to sweep")
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--recover-thresh", type=float, default=0.05,
                    help="min read->true-window containment to call a tail read "
                         "RECOVERABLE (k-mer signal to its origin survives)")
    ap.add_argument("--stage2", action="store_true",
                    help="run the Stage-2 containment-vs-null kill-gate")
    ap.add_argument("--decoy-genome", default=None,
                    help="FASTA tiled into decoy windows for the genome-scale null "
                         "(default: bundled S288C)")
    args = ap.parse_args()
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    rates = [float(x) for x in args.rates.split(",") if x.strip() != ""]
    os.makedirs(args.out, exist_ok=True)

    print(f"[c5] generating corpus reps={args.reps} seed={args.seed} ...",
          file=sys.stderr)
    info = controlled.generate_corpus(args.out, reps=args.reps, seed=args.seed)
    genome = load_genome(info["ref_fa"])
    ref_fa = info["ref_fa"]
    truth = read_truth_table(info["truth_tsv"])
    truth_map = {t.read_id: t for t in truth}

    # original (un-noised) read sequences
    base_reads = {}
    with pysam.FastxFile(info["reads_fastq"]) as fq:
        for e in fq:
            base_reads[e.name] = e.sequence
    print(f"[c5] {len(truth_map)} truth reads; strata="
          f"{sorted(set(t.stratum for t in truth))}", file=sys.stderr)

    # which reads are 'hot' (deterministic from noise-seed)
    hot_rng = random.Random(args.noise_seed)
    rids = sorted(base_reads)
    is_hot = {rid: (hot_rng.random() < args.hot_frac) for rid in rids}

    # pre-sketch true-origin windows (read->window containment for recoverability)
    win_sketch = {}
    for rid, t in truth_map.items():
        gseq = genome.get(t.chrom, "")
        win = gseq[t.genome_start:t.genome_end]
        win_sketch[rid] = fracminhash(win, k=args.k)

    rows = []
    tail_reads_by_rate = {}
    for rate in rates:
        nrng = random.Random(args.noise_seed + int(rate * 1000) + 1)
        noised = []
        cur_seq = {}
        for rid in rids:
            s = base_reads[rid]
            if rate > 0 and is_hot[rid]:
                s = inject_bursty(s, rate, nrng)
            cur_seq[rid] = s
            noised.append((rid, s))
        fq = os.path.join(args.out, f"reads_r{rate:.2f}.fastq")
        bam = os.path.join(args.out, f"mm2_r{rate:.2f}.bam")
        write_fastq(noised, fq)
        run_minimap2(ref_fa, fq, bam)
        cls = classify(bam, truth_map, genome)

        # tally (overall + hot-only, since the tail lives in the hot subpop)
        cnt = defaultdict(int)
        cnt_hot = defaultdict(int)
        c5_tail = []          # rids in the C5-addressable tail
        c5_tail_recoverable = 0
        c5_tail_zeroev = 0
        c4_wrongwin = 0       # paralog wrong-window (NOT C5)
        c3_internal = 0       # right-window-wrong-internal (NOT C5)
        for rid in rids:
            c = cls[rid]
            t = truth_map[rid]
            cnt[c] += 1
            if is_hot[rid]:
                cnt_hot[c] += 1
            # C5 addressability
            is_para = (t.stratum == "PARALOG")
            if c == "right_window_wrong_internal":
                c3_internal += 1
                continue
            if c == "wrong_window" and is_para:
                c4_wrongwin += 1
                continue
            if c == "empty_union" or (c == "wrong_window" and not is_para):
                c5_tail.append(rid)
                # recoverability: does the (noised) read still contain k-mer signal
                # to its TRUE origin window?
                rsk = fracminhash(cur_seq[rid], k=args.k)
                cont = containment(rsk, win_sketch[rid])
                if cont >= args.recover_thresh:
                    c5_tail_recoverable += 1
                else:
                    c5_tail_zeroev += 1
        tail_reads_by_rate[rate] = (c5_tail, cur_seq)
        N = len(rids)
        rows.append(dict(
            rate=rate,
            n=N,
            empty=cnt["empty_union"], wrongwin=cnt["wrong_window"],
            rwwi=cnt["right_window_wrong_internal"], correct=cnt["correct"],
            c5_tail=len(c5_tail), c5_recov=c5_tail_recoverable,
            c5_zeroev=c5_tail_zeroev, c4=c4_wrongwin, c3=c3_internal,
        ))

    # ---- report Stage 1 ----
    print("\n================ C5 PANEL-FAILURE TAIL (minimap2-alone upper bound) ================")
    print(f"hot_frac={args.hot_frac}  k={args.k}  recover_thresh={args.recover_thresh}  reps={args.reps}")
    print(f"{'rate':>5s} {'N':>5s} {'empty':>6s} {'wrongW':>6s} {'rwwi':>5s} {'correct':>7s} "
          f"| {'C5tail':>6s} {'recov':>5s} {'zeroEv':>6s} {'C4wW':>5s} {'C3int':>5s} {'C5tail%':>8s}")
    for r in rows:
        print(f"{r['rate']:5.2f} {r['n']:5d} {r['empty']:6d} {r['wrongwin']:6d} "
              f"{r['rwwi']:5d} {r['correct']:7d} | {r['c5_tail']:6d} {r['c5_recov']:5d} "
              f"{r['c5_zeroev']:6d} {r['c4']:5d} {r['c3']:5d} {100*r['c5_tail']/r['n']:7.2f}%")
    print("\n---- READING ----")
    print("empty   = unmapped by minimap2 (single-member panel: empty-union)")
    print("wrongW  = mapped to WRONG locus/window (reference != truth chrom or no span overlap)")
    print("rwwi    = right window, wrong internal junction/indel (C3/refiner, NOT C5)")
    print("correct = right window + position-exact")
    print("C5tail  = empty UNION (wrongW AND non-paralog)  = the FracMinHash-addressable tail")
    print("recov   = C5tail reads whose noised seq still contains k-mer signal to TRUE window")
    print("zeroEv  = C5tail reads with NO origin signal (unplaceable; never fabricate)")
    print("C4wW    = paralog wrong-window (containment can't disambiguate copies -> C4)")
    print("C3int   = right-window-wrong-internal count (realignment territory)")

    # pre-committed Stage-1 verdict gate
    max_recov = max(r["c5_recov"] for r in rows)
    print("\n---- STAGE-1 TRIGGER (pre-committed) ----")
    if max_recov == 0:
        print(f"  max recoverable C5 tail across all rates = 0  => NO depletion trigger.")
        print("  => DEFER: do NOT commit any FracMinHash index dependency. The panel-failure")
        print("     tail is either empty, junction-internal (C3), paralog (C4), or zero-evidence")
        print("     (unplaceable). Record the measured fractions; revisit only if a future")
        print("     realistic-error corpus produces a recoverable window-level tail.")
    else:
        print(f"  recoverable C5 tail EXISTS (max={max_recov} reads). Trigger is LIVE.")
        print("  => run Stage 2 (--stage2): does FracMinHash containment localize these reads")
        print("     to the CORRECT origin window better than a genome-scale random-window null?")

    # ---- Stage 2 (kill-gate) ----
    if args.stage2 and max_recov > 0:
        run_stage2(rows, tail_reads_by_rate, truth_map, genome, ref_fa, args)
    elif args.stage2:
        print("\n[c5] --stage2 requested but recoverable tail = 0; nothing to localize "
              "(REFUTE branch is reached via the trigger above, not Stage 2).")


def run_stage2(rows, tail_reads_by_rate, truth_map, genome, ref_fa, args):
    """Containment-vs-null kill gate.  Candidate universe = all corpus contig
    windows + DECOY windows tiled from a large genome (S288C) so the random-window
    null is not trivially good.  For each recoverable tail read: argmax-containment
    window; correct iff it overlaps the truth origin span."""
    import rectify
    from pathlib import Path
    from collections import Counter
    W, STRIDE = 300, 150
    # candidate windows = (chrom, start, end). The true origin lives among the
    # corpus-contig windows; S288C tiles are genome-scale DECOYS so the null isn't
    # trivially good. We build an INVERTED INDEX kmer_hash -> [window_id] so
    # localization is O(read_kmers * postings), not O(reads * windows) (the latter
    # is ~45k*5k set-intersections = minutes; the index runs in seconds).
    cand_meta = []          # window_id -> (chrom, start, end, is_decoy)
    cand_size = []          # window_id -> sketch size (for containment denom, unused for argmax)
    inv = defaultdict(list)  # kmer_hash -> [window_id]

    def _add_window(chrom, s, e, seq, is_decoy):
        sk = fracminhash(seq, k=args.k)
        if not sk:
            return
        wid = len(cand_meta)
        cand_meta.append((chrom, s, e, is_decoy))
        cand_size.append(len(sk))
        for h in sk:
            inv[h].append(wid)

    for chrom, seq in genome.items():
        for s in range(0, max(1, len(seq) - 1), STRIDE):
            sub = seq[s:s + W]
            if len(sub) >= args.k:
                _add_window(chrom, s, s + len(sub), sub, False)
    n_corpus = len(cand_meta)
    decoy_path = args.decoy_genome or str(
        Path(rectify.__file__).parent / "data" / "genomes" /
        "saccharomyces_cerevisiae" /
        "S288C_reference_sequence_R64-5-1_20240529.fsa.gz")
    n_decoy = 0
    if os.path.exists(decoy_path):
        dg = load_genome(decoy_path)
        for chrom, seq in dg.items():
            for s in range(0, max(1, len(seq) - 1), W):  # non-overlap decoys, cheaper
                sub = seq[s:s + W]
                if len(sub) >= args.k:
                    _add_window(f"DECOY::{chrom}", s, s + len(sub), sub, True)
                    n_decoy += 1
    print(f"\n================ C5 STAGE-2 containment-vs-null (PROTOTYPE FracMinHash) ================",
          flush=True)
    print(f"candidate windows: {len(cand_meta)} total ({n_corpus} corpus + {n_decoy} S288C decoys), "
          f"W={W} stride={STRIDE}; genome-scale random null ~ origin_windows/{len(cand_meta)}", flush=True)

    print(f"{'rate':>5s} | {'recov_n':>7s} {'recov_acc':>9s} {'recov_decoy':>11s} "
          f"| {'zeroEv_n':>8s} {'zeroEv_acc':>10s} | {'null':>7s}", flush=True)
    rng = random.Random(args.noise_seed + 999)
    for r in rows:
        rate = r["rate"]
        tail_rids, cur_seq = tail_reads_by_rate[rate]
        if not tail_rids:
            continue
        recov, zeroev = [], []
        for rid in tail_rids:
            t = truth_map[rid]
            win = genome.get(t.chrom, "")[t.genome_start:t.genome_end]
            rsk = fracminhash(cur_seq[rid], k=args.k)
            if containment(rsk, fracminhash(win, k=args.k)) >= args.recover_thresh:
                recov.append(rid)
            else:
                zeroev.append(rid)

        def _localize(rid):
            """argmax shared-kmer count over candidate windows (== argmax containment
            since the read sketch size is a constant denominator per read). Returns
            (correct_origin?, best_is_decoy?)."""
            t = truth_map[rid]
            rsk = fracminhash(cur_seq[rid], k=args.k)
            if not rsk:
                return (False, False)
            hits = Counter()
            for h in rsk:
                for wid in inv.get(h, ()):
                    hits[wid] += 1
            if not hits:
                return (False, False)
            best_wid, _ = hits.most_common(1)[0]
            chrom, s, e, is_decoy = cand_meta[best_wid]
            ok = (chrom == t.chrom and s < t.genome_end and e > t.genome_start)
            return (ok, is_decoy)

        rhit = rdecoy = 0
        for rid in recov:
            ok, dec = _localize(rid)
            rhit += ok
            rdecoy += dec
        zhit = 0
        for rid in zeroev:
            ok, _ = _localize(rid)
            zhit += ok
        # genome-scale random-window null (over the whole tail)
        nh = 0
        allt = recov + zeroev
        for rid in allt:
            t = truth_map[rid]
            rc = rng.choice(cand_meta)
            if rc[0] == t.chrom and rc[1] < t.genome_end and rc[2] > t.genome_start:
                nh += 1
        nr, nz = max(1, len(recov)), max(1, len(zeroev))
        print(f"{rate:5.2f} | {len(recov):7d} {rhit/nr:9.3f} {rdecoy/nr:11.3f} "
              f"| {len(zeroev):8d} {zhit/nz:10.3f} | {nh/max(1,len(allt)):7.4f}")
    print("\n---- READING ----")
    print("recov_acc   = OF recoverable tail, freq containment localizes to TRUE origin window")
    print("recov_decoy = OF recoverable tail, freq best window is an S288C DECOY (false-localize)")
    print("zeroEv_acc  = OF zero-evidence tail, freq 'localized' to origin (SHOULD be ~0 = no signal;")
    print("              a high value would expose the recoverability filter as the real gate, not")
    print("              containment specificity)")
    print("null        = genome-scale random-window baseline over the whole tail")
    print("\n---- STAGE-2 VERDICT (pre-committed) ----")
    print("  recov_acc >> null AND zeroEv_acc ~ 0  => mechanism is sound + specific (PROCEED only")
    print("    if Stage-1 ALSO shows a recoverable tail at REALISTIC error; else this is a positive")
    print("    CONTROL of the localizer, and the dep-commit gate stays DEFER pending Tier-2)")
    print("  recov_acc <= null  => REFUTE the FracMinHash fallback (kill-gate failed)")


if __name__ == "__main__":
    main()
