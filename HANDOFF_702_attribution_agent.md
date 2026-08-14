# HANDOFF — 702 Attribution Agent (rectify `attribute-reads` → rbrowse)

**Session:** 2026-08-13 · **Status:** scope complete, committed, not pushed
**Spec:** `~/work/UCLA/Chanfreau_Lab/planning/702_attribution_agent_spec.md` (from rbrowse)
**Design + findings:** `…/planning/703_attribution_sidecar_design.md`
**Peer:** rbrowse agent — coordinated live all session; their consumer is on PROD
**Commits:** `5bb2efd` (the sidecar) · `7f18e56` (a correction) on `feat/attribution-sidecar`

---

## What this is, in one paragraph

rbrowse decided which 3′-end clusters belong to a gene **from the 3′ end alone**. For a
readthrough molecule — one that initiates in gene A and terminates inside gene B — that
counts it *for* B and subtracts it *from* A. In a CPA-depletion experiment readthrough
rises by construction, so the error is **directional and correlated with the genotype being
compared**: it inflates the apparent decrease of exactly the genes the experiment is about.
This session built the per-read sidecar that carries both ends, so a read can be kept with
the gene it *initiated* in. Measured on the anchor-away DRS cohort: **in WT, 0/327
PLB3-initiated molecules reach RCL1; in rna15-AA, 30/58 do.**

---

## 🔴 TWO THINGS THAT WILL LAPSE IF NOBODY OWNS THEM

### 1. The correctness fix is available but NOT YET LIVE on the page

rbrowse has deliberately queued consuming the sidecars rather than starting now. Consuming
them needs a re-slice and bundle re-upload across all 15 libraries — a substantial deploy,
not a viewer change — and since the two implementations already agree on 1,400 of 1,401
reads, it changes almost nothing visible today. Their loader and the numeric registration
(`ATTR_INT`, shard key `ar3`) are in place, so it's a scheduling decision, not a blocker.
Worth knowing that the correctness fix is available but not yet live on the page.

*(Kevin, 2026-08-13, asked for this paragraph in BOTH handoffs so it does not go neglected.
Reproduced verbatim except one false clause — see below. Landed in rbrowse's `HANDOFF.md`
as `877c2e0`, routed via their `.claude/inbox/` rather than written directly, since that
session live-edits the file.)*

### 2. The molecules-vs-reads guard is **NOT BUILT**

An earlier version of the paragraph above claimed "a fail-closed guard against mixing
molecule and read denominators" was "already in place". **That was false, and it was my
error.** rbrowse wrote that such a guard *has to be* built when the page starts consuming
these; I read the future tense as done, asserted it back to them as fact, and it propagated
into Kevin's summary and into this file. rbrowse caught it before pasting.

**What actually exists:** a 🔴 comment beside the loader in `slice_genome.py`. No code.

Why this matters more than a wording slip: the one genuinely hazardous property of the cDNA
arm is that **a row is a MOLECULE, not a read**. "Loader in place, registration in place,
guard in place" reads as *the safety is handled, go ahead and consume*. A future session
trusting a guard that does not exist is worse off than one told plainly the guard is owed.
**Building it belongs inside the consumption job.**

---

## Done

- **Module** `rectify/core/analyze/read_attribution.py` — per-read readthrough-aware gene
  attribution. **No BAM, no re-run of `correct`.**
- **CLI** `rectify attribute-reads` (`core/commands/attribute_reads_command.py`, wired in
  `cli.py`). `--control-tsv` designates the CPA-reference libraries; warns loudly if omitted.
- **25 unit tests** in `tests/test_read_attribution.py`.
- **All 5 spec acceptance tests pass** (`planning/708`, `708b`).
- **Both arms emitted** on H2 scratch:

  | cohort | dir | libraries | rows | unit |
  |---|---|---:|---:|---|
  | anchor-away DRS | `/u/scratch/k/kevinroy/620_rbrowse/attribution/` | 9 | 17,627,566 | reads |
  | by4742 cDNA | `/u/scratch/k/kevinroy/620_rbrowse/attribution_cdna/` | 6 | 21,024,789 | **molecules** |

- **Per-read diff vs rbrowse's `figEscape()`** (`planning/710`): **1,400 of 1,401 agree**.
- **Committed** on `feat/attribution-sidecar`, branched off `chore/vendor-desalt-chanfreau1`
  HEAD. Staged by explicit path; Kevin's unrelated `align` WIP deliberately left out.

## Verified

- `pytest -m "not slow"`: **2343 passed**, 41 skipped, 4 deselected, 1 xfailed.
- `pytest tests/test_read_attribution.py`: **25 passed**.
- **AT1** PLB3→RCL1: rna15_rep1 = **30** (spec ≥30), wt_rep1 = **0**. Across all 9:
  wt 0/0/0 · rna15 30/34/30 · ysh1 36/51/56. ✅
- **AT2** nothing classed `initiating` without cDNA tier-2 evidence — structurally
  unreachable on DRS, so it cannot be reached by accident. ✅
- **AT3** cDNA strand: TSV `strand` attributes **98.4%** vs **51.2%** on its flip,
  reproducing the documented 99.1%-vs-53.4% signature. ✅ *(upper bound — see Open)*
- **AT4** split composition survives (3′ group at 309,250: `YOL010W` 84,
  `YOL011W|YOL010W` 17, `YOL010W|YOL011W` 15). ✅
- **AT5** escaping-read count **invariant** at 4/10/25/50/100/200 bp while cluster count
  moves 31,932 → 5,772. The test a `cluster_id`-keyed output cannot pass. ✅
- **Cross-check vs rbrowse's independent browser implementation** at PLB3:
  WT 0 vs 0 · rna15 94 vs 95 · ysh1 143 vs 143.
- Output inspected directly on both a WT and a mutant file: 12 columns, header comments
  present, `escapes_gene_cpa` populated (40,648 in rna15_rep1), cDNA attribution rate
  98.4% matching AT3 exactly.

## Key findings (all measured this session, none recalled)

1. **The sidecar needs no BAM and no re-correction.** `alignment_start`/`alignment_end` in
   `corrected_reads.tsv` are the raw input alignment. Read retention vs the raw BAM measured
   at **100%** on all 6 libraries (327/327, 58/58, 89/89 …). This *removes* the rbrowse
   ITEM-7 anchor-mismatch class rather than avoiding it.
2. **`five_prime_position` is NOT the observed 5′ end when `five_prime_rescued=1`** — the
   Cat3 rescue extends it upstream, making it model-derived. `origin5` is derived from the
   raw alignment; the rescue flag rides in `origin5_evidence`. The spec stated the
   underlying asymmetry claim unconditionally; it does not hold unconditionally.
3. **🔴 SGD R64 GFF `gene` features are UTR-INCLUSIVE, not CDS-ended** (RCL1: `gene`
   307,883–309,265 vs `CDS` 307,938–309,041 — the union of the per-condition `mRNA` models).
   I gave rbrowse the wrong cause first; they had propagated it into a code comment before I
   corrected it (`36dd3d3`).
4. **Containment-only INVENTS readthrough in WT** — PLB3 WT 0/845 → ~204/845. Confirmed
   independently by rbrowse at 201/845 from the browser side.
5. **No fixed window is defensible** (Kevin's call; the data agrees emphatically). Over
   5,607 genes: **71.6%** have their modal CPA *inside* the annotated interval; ~4% sit past
   +300, up to **+5,024**. `YGL007C-A`: n=2,136 with **46% of reads at one base +3,249**.
   → **Rule B**, chosen by measurement (`planning/707`):
   `boundary = max(annotated_end, observed_modal_CPA) + 281 bp`, where 281 is the *measured*
   genome-wide p95. Rules using each gene's own q95/q99 looked more principled but put
   25/845 and 6/845 **false** escapees into WT.
6. **RCL1's own CPA (309,142) is 133 bp upstream of where PLB3 readthrough terminates
   (309,275)** — the escaping molecules do not use RCL1's CPA at all. Independent positive
   evidence they are readthrough, stronger than the attribution argument the spec was built on.

## Open

- **The single diff disagreement is genuinely ambiguous, not a bug.** Read `dfe0ff08…`,
  3′ end 307,895 (+207 past PLB3): stopped 1,247 bp short of RCL1's own CPA *and* past
  PLB3's entire observed WT range (max +59). Neither rule is clearly right. The tolerance
  tradeoff is deliberate; do not add a third rule for one read in 1,401.
- **NOT VERIFIED: AT3's circularity check is impossible on the shipped cDNA cohort.** The
  11-column UMI schema has no `strand_evidence`, so molecules whose strand was resolved *by*
  gene overlap (tier 3, `ont_cdna.resolve_rna_strand`) cannot be excluded — for those,
  attributing by gene overlap is circular and cannot fail. **The 98.4% is an UPPER BOUND.**
  Needs `strand_evidence` in the cDNA TSV; `rectify correct` emits it (`bam/output.py:44`)
  but that cohort predates it.
- ~~NOT VERIFIED: upf1Δ/WT CPA equivalence~~ → **MEASURED** (`planning/711`). They do **not**
  fully agree, and the WT-only choice is now a justified decision rather than an assumption.
  - 81.4% of genes show **zero** modal shift (median 0, p95 +25 bp), but 7.3% of genes get a
    different boundary.
  - 🔴 **Denominator matters:** only 0.24% of *molecules* change `escapes_gene_cpa`, but that
    is **6.97% of ESCAPE CALLS** (21,088 of 302,497), with a systematic **5.33% net
    reduction** under a upf1Δ-derived reference.
  - **KEEP WT-only. Do NOT re-emit.** NMD degrades long-3′UTR and readthrough products, so
    in upf1Δ those are stabilised — a upf1Δ-derived reference would absorb that phenotype
    exactly as a rna15-derived one would absorb readthrough. The measured sign (fewer
    escapes under the upf1Δ reference) is what absorbing a real phenotype looks like.
  - The script's own auto-verdict says "re-emit with a per-condition reference" — **it is
    wrong**, fired off a threshold picked before seeing data, and is corrected in the
    `INTERPRETATION` section appended to `711_upf1d_cpa_equivalence.md`.
  - ⚠️ **Still not established:** that upf1Δ's shifted 3′-end usage is genuinely phenotypic
    rather than a depth/composition artifact. The argument is mechanistic, not measured.
    Settling it needs a designed test (are the shifted genes enriched for known NMD
    substrates?). **Do not upgrade this to "verified".**
- **`region_class` is empty in v1**, deliberately: the vocabulary is cluster-keyed and the
  spec forbids that join. rbrowse confirmed nothing breaks. Add only when derivable per-read
  from the rules `analyze` already uses — do **not** reimplement the vocabulary.
- **Not pushed.** Two commits sit on `feat/attribution-sidecar` locally.
- H2 `rectify` checkout is OLDER than the M1 tree — its `cli.py` lacks `qc_command`. I
  overwrote it once and **restored it with `git checkout --`**. Drive the module directly on
  H2 (`planning/709`); do **not** sync `cli.py` there.

## Resume

**Nothing is in flight.** Both arms are emitted, verified, committed and handed to rbrowse.

```bash
# Are both output sets still on scratch? (scratch is not backed up)
ssh hoffman2 'cd /u/scratch/k/kevinroy/620_rbrowse && \
  echo "DRS : $(ls attribution      2>/dev/null | wc -l) of 9"; \
  echo "cDNA: $(ls attribution_cdna 2>/dev/null | wc -l) of 6"'
```

| state | do |
|---|---|
| 9 and 6 | nothing to resume — work is complete |
| either short | re-emit; both drivers are idempotent (whole-file writes, no half-file state):<br>`ssh hoffman2 'module load conda/23.11.0; conda activate rectify; cd /u/scratch/k/kevinroy/620_rbrowse && python 709_emit_sidecars.py'` (DRS) · `… python 709b_emit_cdna.py` (cDNA) |
| scratch purged | both drivers rebuild from the corrected TSVs under `/u/project`, which are durable |

**If anything is re-emitted, verify the OUTPUT, not the exit code:**
```bash
ssh hoffman2 'cd /u/scratch/k/kevinroy/620_rbrowse/attribution_cdna && \
  head -1 WT_BY4742_rep1.attribution.tsv; \
  awk -F"\t" "!/^#/ {print NF; exit}" WT_BY4742_rep1.attribution.tsv'
```
MUST print `# unit: molecules` (**not** `reads`) and `12`. 🔴 If it says `reads`, every
downstream count means molecules while claiming reads — **re-emit; do not patch the header.**

**Remaining work is on rbrowse's side and is not blocked here:** the re-slice + bundle
re-upload (9 DRS, then 6 cDNA), and **building the molecules-vs-reads guard**, which does
not exist (§🔴 2). If the page disagrees with the sidecar after switching to `ec`, re-run
`planning/710` — it reimplements their rule and prints per-read *reasons*, not a count.

## Files

| path | what |
|---|---|
| `rectify/core/analyze/read_attribution.py` | the module |
| `rectify/core/commands/attribute_reads_command.py` | `rectify attribute-reads` |
| `rectify/cli.py` | wiring (2 hunks) |
| `tests/test_read_attribution.py` | 25 tests |
| `…/planning/703_attribution_sidecar_design.md` | design + findings F1–F6 |
| `…/planning/704, 705, 705b` | AT1 prototype + the 13-vs-30 reconciliation |
| `…/planning/706, 706b` | genome-wide 3′UTR survey + per-gene inspection |
| `…/planning/707` | tolerance-rule selection (why Rule B) |
| `…/planning/708, 708b` | acceptance tests (DRS, cDNA) |
| `…/planning/709` (+ H2 `709b_emit_cdna.py`) | emission drivers |
| `…/planning/710` | per-read diff vs `figEscape()` |
| H2 `…/620_rbrowse/attribution/`, `…/attribution_cdna/` | the 15 sidecars |
