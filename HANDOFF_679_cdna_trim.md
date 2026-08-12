> # ✅ ANSWERED AND SUPERSEDED by unit 681 (session `cdna-trim-fix`), 2026-08-11 ~22:30.
> Fix written, all four phases landed, acceptance measured. Read
> **`HANDOFF_681_cdna_trim_fix.md`** + `~/work/UCLA/Chanfreau_Lab/planning/681_cdna_consensus_frame_trim_fix.md`.
> **This brief contained three errors that 681 corrected — do not act on them:**
> 1. **§2 "second failing branch: `pileup` 94.8%" was WRONG.** `pileup_consensus` builds from
>    `get_aligned_pairs(matches_only=True)` and therefore excludes soft-clips by construction, so
>    a pileup consensus has no adapter to strip and `XQ==0` there is CORRECT. Its clips are
>    13.7 nt vs 145.2 nt for genuinely-untrimmed singletons. **Acceptance criterion 2 was not
>    merely unmet — acting on it would have BROKEN a correct branch.** Withdrawn.
> 2. **§3 candidates (a)/(b)/(c) were incomplete.** (a) and (c) were refuted; the real mechanism
>    is a fourth: **`orient` is computed in BAM-SEQ frame while the emitted consensus is flipped
>    to basecalled frame, and nothing flips the label alongside the sequence.**
> 3. **§5 criterion 1's "≲10%" floor is unreachable** (~13% is the real floor) and the quoted
>    suite baseline **2,259 was stale** — the real branch-point baseline is 2,216.

# HANDOFF — cDNA consensus adapter/UMI trim is broken (unit 679). Brief for a fresh agent.

**Written 2026-08-11 ~22:10 PDT by session `rectify-realigner` for a dedicated successor.**
You own ONLY this fix. Dataset processing, wave scheduling and the 676 storage work stay with
`rectify-realigner` — do not touch them. Kevin's standing decision: **the DRS arm proceeds
unblocked; every cDNA dataset waits for this fix.** That makes you the critical path.

Evidence record you inherit: **`~/work/UCLA/Chanfreau_Lab/planning/679_cdna_clip_adapter_census.md`**
(checkpoints CP0–CP10). Read it before touching code — it contains two mechanism hypotheses
that were *measured and discarded*, and repeating either wastes a day.

---

## 1. The defect, in one paragraph

`rectify/core/cdna/consensus.py::pretrim_consensus` is supposed to strip the PCB114 SSP + 27-nt
UMI + GGG bridge from the 5′ end and the poly(A)/adapter from the 3′ end of each cDNA consensus
molecule, so aligners see clean mRNA ("no adapter contamination at either end" — its docstring).
It is **silently no-opping on ~half of all molecules**. Those molecules go to the aligner with
~124 nt of adapter + random barcode + poly(A) still attached, which becomes a long terminal
soft-clip. The junction resolver then enumerates hundreds of GT/AG candidates against that clip
and runs a 200-length DP on each — the "1–4 reads/s cDNA pathology." Because a random 27-nt UMI
is *maximally high-complexity*, no complexity gate can reject it; the pipeline's own guard sees
legitimate-looking sequence.

## 2. What is PROVEN (do not re-derive)

All on `587_ms2_cdna/WT_BY4742_rep1`, 200k primary consensus records unless noted.

| fact | number |
|---|---|
| 5′ trim no-op (`XQ:i:0`) | **52.0 %** |
| 3′ trim no-op (`XK:i:0`) | 57.3 % |
| read_type-2 (legitimately has no SSP/UMI) | only **8.2 %** — does NOT explain it |
| clip length, `XQ==0` vs `XQ>0` | **124.8 nt** vs **4.1 nt** mean max soft-clip |
| clips that are adapter-ish | **95.9 %** of 240,229 clips ≥30 nt |
| clips containing a genomic 30-mer | **0.4 %** (shuffled control 0.08 %) |
| **root cause class** | **F1 = consensus FRAME is wrong. F2 (orient label wrong) is REFUTED.** |
| F1 evidence | 119,392 singleton `XO:Z:fwd` molecules: failures = 0 forward / 48,326 reverse (**100.0 % reverse**); successes = 71,058 forward / 8 reverse (**0.0 %**) |
| second failing branch | `XQ==0` by method: `rep` (singleton) 47.9 % vs **`pileup` (multi-read) 94.8 %** |
| `XQ`/`XK` consumers | **none** — only the write site (`cdna/io.py:168-186`) + `multialign/cma_schema.py:34` passthrough |
| junctions fabricated off adapter clips | **0** of 2,590 reads gained an N-op (1 resized). Waste bug, not integrity bug |

**Two mechanisms already tested and DISCARDED — do not re-propose:**
1. *"The exact 23-mer `find()` loses to ~5 % ONT error."* Only the **28.5 %** minority. The SSP
   is present *exactly* in 48.9 % of clips.
2. *"`orient` (XO) is mis-called."* Refuted at 100.0 %/0.0 %. **`XO` is sound, so `XA` poly-A
   length, 3′-end/CPA assignment and strand calls are NOT implicated.** Good news: the blast
   radius is the consensus SEQUENCE only.

## 3. What is NOT yet known — your first task

**The exact defective line is unpinned, and the obvious suspect reads as correct.**
`restore_eq_seq` (`cdna/consensus.py:48`) documents *"Returns the BAM-oriented sequence (same
orientation as `query_sequence`)"*, and `cdna/io.py:136-137` then reverse-complements it when
`seg.is_reverse` — which is exactly what it should do. Yet every failing molecule has a
minus-strand source read. Three live candidates:

- **(a)** that RC is not executing on the path these molecules take;
- **(b)** the segment's strand at consensus time differs from what we measured in `pre.bam`;
- **(c)** `restore_eq_seq`'s documented contract does not hold on the `=`-substitution path.

**Do NOT patch on a code read.** Hypothesis #1 above survived a careful code read and then died
on contact with measurement; that is the failure mode this repo's `CLAUDE.md` is scarred by
("vet alignments at the INDIVIDUAL-READ level at EVERY stage"). Pin it empirically:

```
Take ONE failing molecule: XQ:i:0, XO:Z:fwd, XC:i:1, SSP_RC present, from
  /u/scratch/k/kevinroy/587_ms2_cdna/WT_BY4742_rep1/stage1/stage1_consensus.fastq.gz
(its XR tag holds the source read id). Pull that read from the stage-1 input BAM and step
rectify/core/cdna/io.py:131-137 by hand, printing:
  seg.is_reverse | SSP orientation of restore_eq_seq()'s return | SSP orientation after the RC
That single run discriminates (a) / (b) / (c).
```

## 4. The fix, phased

**Phase 1 — defensive, safe under any verdict, land first.**
Make `pretrim_consensus` orientation-agnostic: search **both** `SSP_FWD` and `SSP_RC`, strip in
whichever frame hits, and treat `orient` as a hint rather than a selector. Recovers ~71 % of
failures by keying on the sequence actually present instead of a label.
**Widen the return contract** to carry left-trim and right-trim separately rather than
overloading `q_trim_5` — the rev branch's own comment (`consensus.py:365-367`) already admits
"right-side trim is not encoded in q_trim_5". This is safe *today* because nothing reads
`XQ`/`XK`, and it makes the tags truthful before anything starts to.

**Phase 2 — the causal fix. REQUIRED regardless of what §3 returns**, because the `pileup`
(multi-read) branch fails at **94.8 %** and is precisely the branch that never applies the frame
restore (`io.py:140-149`; only the singleton path `:136-137` and `rep_fallback` `:153-154` do).
Make frame handling uniform across all consensus branches and **assert** the assembled frame
agrees with `orient` before returning. Fixing only the singleton path leaves every multi-read
molecule untrimmed.

**Phase 3 — fuzzy match + fail-loud counter.** edlib `mode="HW"`, `k=3`, mirroring
`_find_adapter_anchor_pos`, which is *already* fuzzy on the 3′ anchor — the 5′ SSP path simply
never got the same treatment. Covers the ~29 % where neither pattern is exactly present.
`HAS_EDLIB=True` is confirmed in the H2 job env; keep the existing fallback. **Add a trim-no-op
counter to the stage-1 stats** so a future silent regression surfaces as a number rather than as
a 555×-class stall.

**Phase 4 — re-measure BOTH axes on the same slice.** The `XQ==0` rate (expect collapse toward
the ~8 % true type-2 share) **and** resolver reads/s. These are the acceptance numbers.

## 5. Acceptance criteria

1. `XQ==0` rate on `WT_BY4742_rep1` drops from 52.0 % to ≲10 %, and the residual is
   *demonstrably* read_type-2 (cross-tab against `XT`).
2. `XQ==0` rate for `XM=pileup` drops from 94.8 % to the same neighbourhood as `XM=rep`.
3. Mean max soft-clip on the realigned consensus falls from 124.8 nt toward the ~4 nt that
   correctly-trimmed molecules already show.
4. `pytest -m "not slow"` green (baseline: 2,259/0 at `11ed647`).
5. A regression test on the bundled validation reads asserting a known adapter-bearing consensus
   yields `XQ>0` **in both frames** (fwd-labelled and rev-labelled input).
6. Resolver reads/s re-measured on the 2k cDNA slice — reported separately from the numba number.

## 6. Instruments (all verified working tonight)

| what | where |
|---|---|
| trim-tag census (the cheap primary instrument) | `samtools view -F 0x900 <consensus BAM>` → parse `XQ:i:`/`XK:i:`/`XT:i:`/`XO:Z:`/`XM:Z:` |
| clip sequence census + genomic k-mer control | H2 `/u/project/guillom/kevinroy/679_cdna_clip_adapter_census.py`, driver `679_census.sh` → `679_out/` |
| frame diagnostic (F1 vs F2) | H2 `/u/project/guillom/kevinroy/679b_frame_diag.py`, driver `679b_jobs.sh` |
| test substrate — aligned consensus | `/u/scratch/k/kevinroy/587_ms2_cdna/WT_BY4742_rep1/align/WT_BY4742_rep1.minimap2.bam` |
| test substrate — the trim's own INPUT | `.../WT_BY4742_rep1/stage1/stage1_consensus.fastq.gz` (for `XQ==0` records the stored SEQ **is** the untrimmed input — that is why re-running the function on them works) |
| conda env | H2 `~/.conda/envs/rectify` (`HAS_EDLIB=True`) |
| clean tree | branch off `11ed647`; H2 isolated clone pattern in the 668b/674 scripts |

⚠️ **H2 scratch is at 1977.6/2000 GB.** Do not stage new data there. Write outputs to
`/u/project/guillom/kevinroy/`. A 313 GB move (`move676par`, job 14293343) is freeing space —
`rectify-realigner` owns it; leave it alone.

## 7. Durability requirement (non-negotiable)

Claim your own `NN` atomically (`~/work/bin/nn_claim.sh ~/work/UCLA/Chanfreau_Lab/planning
<owner> <slug>`) — **do not write into 679**, which is a closed evidence record. Write your
plan to `NN_*.md` **before** any compute, then append a checkpoint after each sub-step (line
pinned → Phase 1 written → tests green → re-measured). An agent that persists nothing is
unrecoverable when the API blips; the checkpoint cadence IS the recovery granularity.

## 8. Traps

- **Don't "fix" the raw-read side.** Raw reads are deliberately untrimmed — that is Path A by
  design ([[581]]), because `trim-cdna-polya` was destroying the XF signal ([[580]] CP1) and
  because `extract_read_info` must see the SSP in the raw read to slice the UMI. Only the
  *consensus* trim is broken.
- **The `679_census.py` auto-VERDICT line prints "FORK B" and is WRONG** — its threshold keyed
  on a fuzzy-but-not-exact signature that the frame-mismatch mechanism routes around. The tables
  under it are right.
- **`XM=poa` never appears in this data** — POA is not engaging; multi-read clusters land in
  `pileup`. Don't optimise a branch that isn't running, but do fix its frame handling too.
- Don't plan cDNA wave capacity off the numba kernel's rate (`277c708`, 289×, bit-identical,
  should still land). It makes the DP fast; **your fix removes most of the DP.** Phase 4's
  number is the one that gates wave sizing.
- Surgical staging only: `git add <explicit paths>`, never `git add -A` (Kevin routinely has WIP
  in `rectify/core/bam/*`, `rectify/core/correct/*`, `scripts/validation_data/`).

## 9. Report back

Drop your result in `~/work/rectify/.claude/inbox/` addressed to `rectify-realigner` — the cDNA
fan-out cannot be scheduled until Phase 4's two numbers exist. Kevin wants the finding in text
the moment you have it, not after the writeup.
