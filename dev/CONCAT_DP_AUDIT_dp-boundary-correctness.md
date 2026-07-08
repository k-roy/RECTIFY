# Concat-DP Audit — Lens: DP BOUNDARY CORRECTNESS

## Task
Adversarial audit of concat-reference DP optimization of `_score_junction`
(gated behind `_USE_CONCAT_DP`, default False) in
`rectify/core/splice/junction_scoring.py`.

Find ONE input where OLD (`_USE_CONCAT_DP=False`) != NEW (`_USE_CONCAT_DP=True`)
on exact float. Byte-identity is sacred.

## Lens
DP boundary correctness. Attack junction handling at concat boundary
(intron_end-1 | intron_end): cross-boundary del/ins, HP-run del cost
computed across boundary vs split, k=0/k=L extremes, 1bp rescue,
rescue exactly at boundary.

## Plan
1. Read _score_junction, _score_hp_anchored, _score_junction_concat.
2. Understand del_costs / ins_costs construction (HP-aware).
3. Build in-process harness toggling _USE_CONCAT_DP.
4. Construct adversarial boundary cases; compare OLD vs NEW exact float.
5. Persist any divergence immediately.

## Checkpoints
- (start) file created, beginning orientation.

## STATUS OF BUILD
- concat DP was NEVER implemented in junction_scoring.py. No `_USE_CONCAT_DP`
  flag, no `_score_junction_concat`. Builder self-report = null; BUILD doc
  checkpoints log is EMPTY. Therefore I audited the DESIGN (from the FIX spec
  + dev/CONCAT_DP_BUILD.md), implementing NEW two independent ways and comparing
  to the REAL OLD `_score_junction`.
- numba OFF (`js._NUMBA_AVAILABLE == False`) — pure-Python DP path runs. Good.

## HARNESS FIDELITY
- new_decomp(kmax_extra=0) == OLD `_score_junction` on 2000/2000 random cases
  (0 mismatches). My replication of the split loop is byte-faithful, so any
  divergence at kmax_extra=1 is purely the added k=L term.

## DIVERGENCE FOUND — YES. checkpoint: DIVERGE

### Root cause (DP BOUNDARY): the k=L (all-intron-side) path
OLD loops `for k in range(L)` => k in [0, L-1]; it STRUCTURALLY EXCLUDES k=L.
The FIX's single concat DP has a FREE ref-suffix, so a path that consumes all
rescue bases within ref_intron_end (columns <= boundary R2) and then treats the
boundary + entire exon2 as free suffix IS a valid path — that path == k=L
(all rescue on the intron side, exon2 alignment empty=0). The design's own
decomposition ("last cell in boundary column R2 at row k") ranges k over [0, L],
i.e. INCLUDES k=L. Whenever rescue aligns to the intron-end window better than to
exon2, NEW < OLD.

### MINIMAL reproducing input (penalty_table=None; degenerate 1bp rescue)
  query="A", q_split=0, intron_end=100,
  genome_seq = "T"*50 + "GCGC...CGA"(50bp, genome[99]='A') + "C"*60
  genome[99]='A' (last intron base), exon2[0]=genome[100]='C'.
  Per-k: k=0 -> t1=1.0 t2=0 sum=1.0 (OLD's only term);
         k=1(=L) -> t1=0.0 t2=0.0 sum=0.0  (EXCLUDED by OLD range(L); INCLUDED by concat DP)
  OLD  `_score_junction` = 1.0
  NEW  (single concat DP / design min_k[t1+t2] over [0,L]) = 0.0
  ==> OLD != NEW. BYTE-IDENTITY BROKEN.

### Scale of the problem
- 512 / 5000 random cases (len 1..20, random genome) DIVERGE (OLD != design[0,L]).
- 8bp adversarial case (rescue == last 8 intron bases, exon2 all 'C'):
  OLD=2.0, NEW(decomp[0,L])=0.0, NEW(singleDP)=0.0.

### SECOND, independent flaw (bonus): design's equivalence claim is FALSE
My literal single concat DP (free intron prefix / penalized exon2 prefix / free
suffix, del_concat=reverse(del_costs_rev)++del_costs_fwd) DISAGREES with
min_k[t1+t2] on many cases (e.g. query='GG' intron_end=93: decomp=1.0 vs
singleDP=0.0). The naive single DP loses t2's near-junction right-anchor for k<L
(free prefix over the whole intron window lets rescue[:k] float to any intron
match), so it is EVEN MORE permissive than [0,L]. So BUILD doc line 51
("concat DP == min_k[t2(k)+t1(k)]") does not hold; whichever way NEW is built it
diverges from OLD, always NEW <= OLD, frequently NEW < OLD.

### VERDICT: byte-identity does NOT survive. Divergence is robust
(penalty_table=None, 1bp rescue). checkpoint: audit complete.
