# HANDOFF — guard cost/benefit re-eval on the corrected residual (2026-07-17)

## Context
The compensating-indel fix (commit e40ca00) removed ~91% of arm-B fabrication (see
`dev/REPLACER_COMPENSATING_INDEL_BUG.md`). Now re-evaluating whether the microhomology guard /
positional close (`_positional_signal`, `drift_positional_gate`) is still justified on the corrected
residual (§4b FIXED: ~3,488 drift + ~1,462 recovery at ≥2 samples), or should be shelved to simplify.
Advisor framing: the decision hinges on whether `_positional_signal` SEPARATES the residual drift
(fabrication, read matches incumbent → signal ≤0) from recovery (genuine, read matches moved → >0).
Bias toward simplification — make the close prove separation. Guard was originally tuned against a
signal that was 80–90% artifact.

## Done
- Deployed the FIXED `junction_refiner.py` to `/scratch/users/kevinroy/rectify_guard2` (has both the fix
  and `_positional_signal`; py3 imports OK).
- Patched `recall_analyze_v2.py` with `--dump-4b <tsv>` to dump (chrom,donor,acceptor,category,nsamp) for
  the FIXED revealed set.
- Wrote `pos_signal_eval.py` (on Sherlock `sma_recall/`): for every arm-B ACCEPTOR-ONLY move whose refined
  junction is in the dump, computes per-read `_positional_signal(genome, q, q_split, ne=raw_acc,
  new_je=refined_acc)`, aggregates per junction, and prints per-category separation (recovery vs drift).

## Open / in flight
- ✅ **Job 34292684** (`dump4b`) DONE. `fixed_4b_dump.tsv`: 3488 drift, 1462 recovery, 18327 inconclusive
  (≥2 samples). Drift is **DIFFUSE** — 3488 junctions across 2368 distinct donor-kb loci, all chromosomes
  (~1.5/locus). Diffuse = weaker case for the guard (not concentrated real biology).
- ⏳ **Job 34292864** (`poseval`) — the decisive `_positional_signal` separation test. Writes
  `pos_signal_eval.tsv` (per-junction) + `pos_signal_eval.summary.txt` (stderr: per-category separation).
  Sentinel `.poseval_rc`. IN FLIGHT (~20-40min, streams 13 BAM pairs).

## Resume
1. Check dump: `ssh sherlock 'cat /scratch/users/kevinroy/sma_recall/.dump4b_rc 2>/dev/null; wc -l < /scratch/users/kevinroy/sma_recall/fixed_4b_dump.tsv'`
   - If rc==0 and file non-empty → (a) check drift concentration:
     `ssh sherlock 'awk -F"\t" "\$4==\"drift\"" /scratch/users/kevinroy/sma_recall/fixed_4b_dump.tsv | cut -f1 | sort | uniq -c | sort -rn | head'`
     (concentrated in few loci = more likely real biology / tractable; diffuse = weaker case for the guard).
   - If PENDING/RUNNING → wait. If FAILED → check `logs/dump4b-*.log` + `dump4b.err`.
2. Run the positional-signal separation (submit as sbatch, 24-48G, conda `rectify`):
   `python pos_signal_eval.py --dump fixed_4b_dump.tsv --out pos_signal_eval.tsv` (streams 13 BAM pairs;
   stderr prints the per-category separation summary). Read the stderr summary + `pos_signal_eval.tsv`.
3. DECISION rule: if recovery median >0 and drift median ≤0 with clean separation (e.g. recovery frac>0
   high, drift frac<=0 high) → the close discriminates on the clean residual → KEEP + tune the gate.
   If distributions OVERLAP → the close does not earn its keep on real data → SHELVE the guard machinery,
   ship the fix alone. Report to Kevin either way. Update `dev/REPLACER_COMPENSATING_INDEL_BUG.md`.
   Then `rm -f .claude/.handoff-needed`.

## Files
- Sherlock `sma_recall/`: `recall_analyze_v2.py` (--dump-4b), `pos_signal_eval.py`, `dump_4b.sbatch`,
  `fixed_4b_dump.tsv` (pending), `lr_junc_fixed/` (corrected junctions), `recall_result_{OLD,FIXED}_matched.txt`.
- Sherlock `rectify_guard2/` — FIXED junction_refiner.py deployed.
- Worktree: this handoff; `dev/REPLACER_COMPENSATING_INDEL_BUG.md` (fix + §4b result).
