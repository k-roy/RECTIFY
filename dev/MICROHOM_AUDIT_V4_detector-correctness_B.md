# MICROHOM AUDIT V4 — task=detector-correctness — AUDITOR B

Independent adversarial audit of `_move_microhomology` + `_effective_veto_margin`
(commit 05664bc) in `rectify/core/splice/junction_refiner.py`.

Auditor B, ATTEMPT 1. READ-ONLY. Python: /Users/kevinroy/miniconda3/bin/python
Scratch: /private/tmp/.../audit_v4/detector-correctness_B

## LENS
Helper/detector correctness. Do the two helpers ever score WRONG for some input?
Priority: A5 (N==N counts as match in `_frac_match`) and A8 (max-over-both-boundaries
masking). Then geometry both strands, genome edges, k=0/1/huge, nested repeats,
both-boundary shift. And `_effective_veto_margin` all regimes.

## CONFIRMED MECHANICS (read from source)
- `_frac_match(a,b)`: `sum(1 for x,y in zip(a,b) if x==y)/len(a)`; 0.0 if empty or
  len mismatch. **Raw `x==y`. NO `.upper()`, NO exclusion of 'N'.** So:
  - `N==N` scores as a MATCH (1.0 for an all-N pair).
  - lower-case vs upper-case ('a' vs 'A') scores as MISMATCH.
- `_move_microhomology(seq, ns, ne, js, je)`:
  - acceptor leg (je!=ne): lo,hi = sorted(ne,je); k=hi-lo; if hi+k<=n:
    best=max(best, frac(seq[lo:hi], seq[hi:hi+k]))  — DOWNSTREAM repeat check.
  - donor leg (js!=ns): lo,hi=sorted(ns,js); k=hi-lo; if lo-k>=0:
    best=max(best, frac(seq[lo-k:lo], seq[lo:hi]))  — UPSTREAM repeat check.
  - returns max over the two legs. **A8 substrate = the `max`.**
  - out-of-bounds leg contributes nothing (guarded by hi+k<=n / lo-k>=0), so returns
    0.0 if the only shifted leg is OOB. No crash.
- Caller (line 816-820): `if _move_microhomology(genome_seq, ns, ne, new_js, new_je)
  >= microhom_threshold: eff_margin += microhom_drift_margin`. **genome_seq passed
  RAW (whatever case/N the FASTA loader produced).** microhom_threshold default 0.5.
- `_effective_veto_margin(hold, eff, cap)`: `if cap>0 and eff>hold: return max(hold,
  min(eff,cap)); return eff`.
- Veto (line 838): `if best_score_cmp > incumbent_score - veto_margin: moves=False`.
  score_cmp LOWER=better; a LARGER veto_margin => easier to veto (holds incumbent).

## SUBSTRATE QUESTION for A5
Does the reference genome ever contain 'N'? junction_scoring.py line 232 comment:
"only A/C/G/T/N (real bases)". So the codebase EXPECTS N in genomes (assembly gaps,
telomeres, IUPAC). Also case: is genome_seq upper-cased before reaching the refiner?
-> MUST verify how genome dict is built (soft-masked lowercase repeats are common).

## SUBSTRATE RESOLUTION (decisive)
- `rectify/utils/genome.py::load_genome` (lines 152/154): `str(record.seq).upper()`.
  => genome_seq reaching the refiner is ALWAYS upper-cased.
  * CASE bug (lowercase soft-mask false-negative): NEUTRALIZED in production. Not reachable
    via load_genome. Record & drop (per advisor).
  * 'N' bug: `.upper()` does NOT strip N. All-N regions (assembly gaps / telomeres / IUPAC)
    reach the refiner as 'NNNN...'. `_frac_match(N..,N..)` counts N==N as MATCH. A5 REACHABLE.
- Upstream pre-filter check (BLOCKING per advisor): grepped lines 640-790. Candidates come
  from the observed+annotated junction POOL, gated ONLY by boundary_error_window /
  max_boundary_shift / max_junction_size / (canonical-tier for the boundary-CLEAN skip).
  **NO N-region or non-GT-AG exclusion before line 817.** In motif_blind (audited) mode,
  `_move_microhomology` at 817 sees the raw genome slice incl. N. => A5 NOT pre-filtered.

## CHECKPOINTS
- [x] source read, mechanics confirmed, durable record written
- [x] substrate resolved: load_genome uppercases (case neutralized), N survives, no N pre-filter
- [x] harness written to disk (scratch/harness.py)
- [x] unit numbers (unit_results.json)
- [x] A5 end-to-end DECISION FLIP proven (e2e2.py / e2e2_results.json)
- [x] A8 unit-level max-masking proven (harness.py)
- [x] _effective_veto_margin 10-row truth table: ALL PASS, invariant hold<=got<=eff holds
- [ ] A5 realism/severity (does a real junction move ever sit adjacent to N?)

## RESULTS — A5 (N==N phantom microhomology) — **HOLDING-CLASS DETECTOR FAULT**
`_frac_match` scores N==N as a match; `_move_microhomology` inherits it; genome_seq
reaches the refiner with N intact (load_genome uppercases but does not strip N); no
upstream pre-filter blocks N-region moves in motif_blind mode.

UNIT (harness.py):
- `_frac_match("NNNNNN","NNNNNN")` = **1.0** (correct N-excluding = 0.0).
- minimal tip "ACNTTT" vs "ACNGGG": actual **0.5** (TRIPS), correct 0.333 (does not).

END-TO-END DECISION FLIP (e2e2.py — honest-CIGAR read that genuinely supports a +6
acceptor move, delta_improve = 6.0 edit units):
| variant | seg[A:A+6] | seg[A+6:A+12] | actual mh | correct(no-N) mh | guard OFF acc | guard ON(m=8) acc | phantom veto |
|---|---|---|---|---|---|---|---|
| V1_real_repeat | ACGTAC | ACGTAC | 1.0 | 1.0 | A+6 (move) | A (held) | NO (real mh — intended) |
| V2_N_phantom | ACNTGA | ACNCTT | **0.5 TRIPS** | **0.333 spare** | A+6 (move) | A (held) | **YES** |
| V3_full_N | NNNNNN | NNNNNN | **1.0 TRIPS** | **0.0 spare** | A+6 (move) | A (held) | **YES** |

In V2/V3 a genuine, strongly read-supported move (delta_improve=6.0) is VETOED
purely because N==N fabricates microhomology >= 0.5. A correct detector (N excluded)
returns < 0.5 and would SPARE the move. FALSE-POSITIVE VETO — the helper computes a
WRONG value on a realizable input. Cap bounds it: with drift_near_tie_cap=2.0 the move
survives in V2/V3 (delta_improve 6.0 > cap), but the DEFAULT cap=0.0 leaves the phantom
veto UNBOUNDED across the whole (hold, hold+microhom_drift_margin) delta_improve band.

## RESULTS — A8 (max-over-both-boundaries) — over-broad BY DESIGN, direction = over-veto only
harness.py A8_move: dual-shift where acceptor leg = 0.0 (benign transition, should be
spared) but donor leg = 1.0 (genuine repeat). `max` => actual **1.0**, masking the benign
acceptor. Direction: max can only INFLATE -> only ever causes "real move vetoed", never
"fake move accepted". BUT: a dual-shift whose donor sits in genuine microhomology IS
partly the drift the guard targets, and the return value 1.0 is a CORRECT report of the
donor leg's real microhomology — the docstring states "returns the max over the moved
boundaries". So A8 is NOT a wrong-value computation; it is an over-broad SCOPE that
collapses into the ALREADY-CLEARED read-blind structural point (one score axis; a real
acceptor cryptic co-moving with a donor repeat is indistinguishable from a drift). NOT a
separate HOLDING fault. (With N in the donor flank, A8 stacks on top of A5 — the max can
be tripped by a PHANTOM donor repeat — but that is A5, not a distinct A8 fault.)

## RESULTS — _effective_veto_margin — CORRECT
10-row truth table all pass; invariant hold <= got <= eff holds every regime incl.
hold>cap (->hold, never below hold), cap<=0 (->eff), eff==hold (->eff), cap>=eff (->eff),
eff==hold==cap, eff==hold==0. No fault.
- [ ] A5 reproducing case (N-run false-positive veto of a real move)
- [ ] A8 reproducing case (unrelated repeat on other boundary masks benign transition)
- [ ] geometry/edges/k grid
- [ ] _effective_veto_margin regime table
- [ ] VERDICT
