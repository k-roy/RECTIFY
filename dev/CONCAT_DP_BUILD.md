# CONCAT-REFERENCE DP for `_score_junction` — build log (SYSTEM OF RECORD)

Agent worktree: agent-a25a2c1e784ad37dc. Target:
`rectify/core/splice/junction_scoring.py::_score_junction`.
Python: /Users/kevinroy/miniconda3/bin/python (numba NOT installed — keep off).

## Goal
Replace the bilateral k-sweep (~2L pure-Python `_score_hp_anchored` DP calls
per candidate) with a SINGLE concat-reference DP, gated behind module flag
`_USE_CONCAT_DP` (DEFAULT FALSE). Old loop stays production default.
Must be byte-identical (== float) to the old loop over >=20000 random cases.

## The split (flag=False, unchanged) — semantics
`_score_junction` returns `min_k [ t1(k) + t2(k) ]`, k in [0, L), L=len(rescue<=30).
- t1(k) = `_score_hp_anchored(rescue[k:], ref_exon2_start, del=del_costs_fwd)`
  left-anchored: query GLOBAL, ref-PREFIX penalized (leading exon2 dels paid),
  ref-SUFFIX FREE (`min(prev)`). ref_exon2_start = genome[intron_end : +BUF].
- t2(k>0) = `_score_hp_anchored(reverse(rescue[:k]), ref_intron_end_rev,
  del=del_costs_rev)`; t2(0)=0. Reversal trick ⇒ right-anchored: rescue[:k]'s
  last base anchored to end at intron_end-1; near-junction dels paid, far dels free.
- ins costs: penalty_table=None ⇒ constant 1.25. Else per-query-base
  `ins_cost(hp_run_length(query_slice, i), query_slice[i])` — SLICE-relative.
- Column-0 quirk: `curr[0] = i * ins_costs[i-1]` (NOT cumulative sum).

## The concat DP (flag=True) — design
Concat reference C = ref_intron_end (forward) ++ ref_exon2_start, i.e. a
CONTIGUOUS genome window genome[intron_end-BUF : intron_end+BUF]. Boundary
column = R2 = len(ref_intron_end). Single DP: query=rescue (global), FREE outer
ends (free ref prefix before intron_end window = t2 free far-left; free ref
suffix after exon2 = t1 free right). del_costs_concat =
reverse(del_costs_rev) ++ del_costs_fwd (each ref base's HP-aware cost stays at
its own genomic column — NO reordering across boundary).

### Boundary handling (deletions) — PROVEN clean
- A deletion consumes exactly ONE ref base; every ref base is in exactly one side
  (left C[<R2] or right C[>=R2]), so no single deletion spans the boundary.
- A *run* of deletions crossing R2 (delete trailing left bases then leading right
  bases at the same query row) IS allowed in BOTH: t2(k) pays trailing near-junction
  left dels, t1(k) pays leading exon2 dels — same row, no query base between.
  concat reproduces this via horizontal moves across column R2, each using its own
  column's del cost. del cost order preserved ⇒ byte-identical for deletions.
- del_costs_rev[j] is the cost of genome[intron_end-1-j]; reversing it gives the
  FORWARD left-ref del costs (cost depends only on genome pos via HP-run/STR, which
  is reversal-invariant). So reverse(del_costs_rev) == _precompute_del_costs(forward
  ref_intron_end). Use reverse(del_costs_rev) to guarantee identity + reuse.

### Decomposition argument (why min_k[t1+t2] == concat DP)
Every monotone path over C has a LAST cell in column R2, at row k. Prefix =
align rescue[:k] to left-ref ending at R2, free left prefix = t2(k) (fwd
right-anchored). Suffix = align rescue[k:] to right-ref starting at R2, free
right suffix = t1(k). concat DP = min over paths = min_k[t2(k)+t1(k)]. ✓ for del/sub.

### KNOWN RISK (insertions at the boundary) — must verify empirically
Two quirks make a *naive fixed-cost* concat DP potentially NOT byte-identical
when penalty_table != None:
1. ins costs are SLICE-relative: `hp_run_length(rescue[k:], ·)` truncates an HP
   run that straddles k; full-rescue context does not. k-dependent ⇒ not a fixed
   array.
2. Column-0 quirk `i*ins_costs[i-1]` is non-additive; only collapses to a cumulative
   sum when ins is CONSTANT (penalty_table=None).
⇒ Expectation: penalty_table=None should be byte-identical; penalty_table!=None MAY
diverge on straddling-HP + boundary-insertion cases. WILL TEST AND REPORT HONESTLY.

## Plan / checkpoints
- [ ] design done
- [ ] concat DP coded behind `_USE_CONCAT_DP` (default False)
- [ ] harness written (>=20000 cases: len 1..30, varied windows, HP runs,
      penalty_table None + real DRS table, fwd + reversed refs)
- [ ] harness run: N cases, mismatches=M, first=...
- [ ] pytest suite pass (flag False AND True)
- [ ] speedup measured (old vs new, batch of few thousand); DP-call/cell-op count
- [ ] NO commit, default stays False, numba NOT installed

## Checkpoints log
