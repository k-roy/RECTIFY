# HANDOFF — compensating-indel refiner fix + §4b re-run (2026-07-15)

## Done
- **Root-caused the phantom junction "drift"** on Sumner SMA ONT DRS. Single-read IGV/composite
  inspection (render_single.py) + per-base CIGAR dump showed the re-placer emits DEGENERATE
  compensating-indel CIGARs: a both-boundary junction move rendered as `I(k)…N…D(k)` (pure slide) or
  `D(a)…N…D(b)` (asymmetric). The exon SEQUENCE stays put; only the N-op COORDINATE moves → a false
  junction position. ~85–95% of "moved" reads at the leads. Full writeup: `dev/REPLACER_COMPENSATING_INDEL_BUG.md`.
- **The §4b "27% fabrication" was largely this artifact** (it measured N-op coordinates via
  extract_sr_junctions.py). Corrected true-acceptor metric: PCBP2 54%→4%, UBA1 8%→3%, SNRPN's "+6
  genuine" is ALSO a degenerate pure slide (D6 before N + I6 after).
- **FIXED** `_apply_junction_replacement` (junction_refiner.py): refuse a BOTH-boundary move that
  RAISES the read's I+D burden (the compensating-indel signature). Single-boundary corrections and
  clean microhomology slides (fast path) untouched. Committed **e40ca00** on branch
  `worktree-agent-a25a2c1e784ad37dc` (fix + 4 regression tests + findings doc).
- **Verified**: `pytest tests/test_junction_refiner.py tests/test_microhom_drift_guard.py` (66 pass) +
  `tests/test_validation_reads.py` (108 pass). Confirmed the fix refuses the real PCBP2 `a53fcd26`
  −39/−39 slide locally (applied=False, CIGAR unchanged).

## Verified
- Corrected extractor (`extract_junctions_corrected.py`) faithfully applies the fix's exact criterion
  to EXISTING refined.bam+sub.bam (reverts both-boundary+adjacent-indel N-ops to raw). On the SMA_191
  lead slice: 10709/10993 moved N-ops (97%) reverted as phantom; PCBP2 phantom 53471598 gone, canonical
  53471637 restored (851 reads). NOT VERIFIED genome-wide yet (job running).

## Open / in flight
- ✅ **Job 34074461** (corrected LR-junction re-extraction, array 1-13) COMPLETED. Genome-wide result:
  **~94% of moved N-ops are phantom** (reverted to raw); ~5.5% single-boundary-kept, ~0.4% clean slide.
  Output in `/scratch/users/kevinroy/sma_recall/lr_junc_fixed/` (13 `.refined.junc.tsv` + `.sub.junc.tsv`).
- ✅ **Job 34076391** (`recall4bfixed`, --lr-dir lr_junc_fixed) DONE — but its §4b was UNCHANGED because
  recall_analyze_fast.py's §4b reverse check reads the HARDCODED `panel_deep/*.revealed_noncanon.tsv`
  (line 207), NOT --lr-dir (--lr-dir only feeds §3 recall). So §4b never saw the corrected junctions.
- Reverse-engineered: non-canonical arm-B junctions drop **76923 → 5986 at ≥2 samples (−92%)** under the
  fix (raw lr_junc vs lr_junc_fixed, canonical==0). Confirms the reveals were phantom.
- Patched script `recall_analyze_v2.py`: added `--revealed-from-lrdir --rev-min-split N` to build the §4b
  revealed set from `--lr-dir *.refined.junc.tsv` (canonical==0) for an apples-to-apples old-vs-fixed §4b.
- ⏳ **Job 34078097** (`recall4bcmp`) — runs v2 §4b on BOTH lr_junc (OLD matched) and lr_junc_fixed (FIXED
  matched) → `recall_result_{OLD,FIXED}_matched.txt` + sentinel `.recall4bcmp_rc`. IN FLIGHT.
- The original panel refine DRIVER that produced `panel_deep/*.refined.bam` was not found (ad-hoc/removed).
  The corrected-extractor path sidesteps it by applying the fix to the existing BAMs. If a full
  through-the-refiner re-run is wanted for publication, reconstruct the driver (inputs: sub.bam,
  aligner_bams=[sub.bam]?, GENCODE v44 annotated set, GRCh38) and deploy the fixed junction_refiner.py.

## Resume
1. Check the matched comparison job:
   `ssh sherlock 'cat /scratch/users/kevinroy/sma_recall/.recall4bcmp_rc 2>/dev/null; sacct -j 34078097 -X -o State'`
   - If `.recall4bcmp_rc` == 0 → read BOTH `recall_result_OLD_matched.txt` and
     `recall_result_FIXED_matched.txt`; compare §4b (recovery/drift/inconclusive counts + total revealed
     non-canonical at ≥2/≥3/≥5 samples) and §3 recall. The OLD_matched is the correct baseline (same
     revealed-from-lrdir filter), NOT the published recall_result_fast.txt (different generator).
   - If missing/RUNNING → wait. If FAILED → check `recall_result_*_matched.err` + `logs/recall4bcmp-*.log`.
   - RE-FIRE: `ssh sherlock 'cd /scratch/users/kevinroy/sma_recall && sbatch run_4b_compare.sbatch'`.
2. Report the OLD-vs-FIXED §4b (esp. total revealed non-canonical and drift count) to Kevin.
   Update `dev/COMPASS_RECALL_RESULT.md` §4b with the corrected numbers.
   Then `rm -f .claude/.handoff-needed` (once reported).

## Files
- `rectify/core/splice/junction_refiner.py` — the fix (commit e40ca00)
- `tests/test_junction_refiner.py` — regression tests
- `dev/REPLACER_COMPENSATING_INDEL_BUG.md` — mechanism writeup
- `/Users/kevinroy/sumner_igv_slices/{render_single,extract_junctions_corrected}.py` — local tools
- Sherlock: `/scratch/users/kevinroy/sma_recall/{extract_junctions_corrected.py,lr_extract_fixed.sbatch,run_4b_fixed.sbatch}`
- Sherlock outputs: `lr_junc_fixed/` (corrected junctions), `recall_result_FIXED.txt` (pending)
