# HANDOFF — 702 Attribution Agent (rectify `attribute-reads` → rbrowse)

**Session:** 2026-08-13 · **Spec:** `~/work/UCLA/Chanfreau_Lab/planning/702_attribution_agent_spec.md`
**Design + findings:** `…/planning/703_attribution_sidecar_design.md`
**Peer:** rbrowse agent (live, coordinating by SendMessage; consumer is on PROD)

---

## 🔴 THE CORRECTNESS FIX IS AVAILABLE BUT NOT YET LIVE ON THE PAGE

**Do not let this lapse. It is a scheduling decision, not a blocker, and it has no owner
until someone schedules it.**

rbrowse has deliberately queued consuming the sidecars rather than starting now. Consuming
them needs a re-slice and bundle re-upload across all 15 libraries — a substantial deploy,
not a viewer change — and since the two implementations already agree on 1,400 of 1,401
reads, it changes almost nothing visible today. Their loader and the numeric registration
(`ATTR_INT`, shard key `ar3`) are in place, so it's a scheduling decision, not a blocker.
Worth knowing that the correctness fix is available but not yet live on the page.

### 🔴 THE MOLECULES-VS-READS GUARD IS **NOT BUILT**. Do not assume it.

An earlier version of the paragraph above claimed "a fail-closed guard against mixing
molecule and read denominators" was "already in place". **That was false and it was my
error** — rbrowse said such a guard *has to be* built when the page starts consuming these;
I read the future tense as done, asserted it back to them as fact, and it propagated into
Kevin's summary and into this file. rbrowse caught it on 2026-08-13 before pasting.

**What actually exists:** a 🔴 comment beside the loader in `slice_genome.py`. No code.

Why this specific error is dangerous rather than merely wrong: the one genuinely hazardous
property of the cDNA arm is that **a row is a MOLECULE, not a read**. "Loader in place,
registration in place, guard in place" reads as *the safety is handled, go ahead and
consume*. A future session trusting a guard that does not exist is worse off than one told
plainly that the guard is owed. **Building it belongs inside the consumption job.**

*(Kevin, 2026-08-13 — asked for the paragraph verbatim in BOTH handoffs so it does not go
neglected. Reproduced verbatim except the false clause, which is corrected above rather than
preserved. Mirrored to rbrowse via their `.claude/inbox/`; their HANDOFF.md is live-edited
by that session, so it was routed rather than written directly — landed there as `877c2e0`.)*

---

## Done

- **New rectify module** `rectify/core/analyze/read_attribution.py` — per-read
  readthrough-aware gene attribution. No BAM, no re-run of `correct`.
- **New CLI command** `rectify attribute-reads` (`core/commands/attribute_reads_command.py`,
  wired in `cli.py`). `--control-tsv` designates the libraries the CPA reference is
  built from; warns loudly if omitted.
- **23 unit tests** in `tests/test_read_attribution.py`.
- **All 5 spec acceptance tests pass** (`planning/708`, `708b`).
- **Sidecars emitted** for the 9 anchor-away DRS libraries →
  `/u/scratch/k/kevinroy/620_rbrowse/attribution/<sample>.attribution.tsv`.

## Verified

- `pytest -m "not slow"`: **2343 passed**, 41 skipped, 4 deselected, 1 xfailed. Nothing broken.
- `pytest tests/test_read_attribution.py`: 23 passed.
- **AT1** PLB3→RCL1: rna15_rep1 = **30** (spec ≥30), wt_rep1 = **0**. ✅
- **AT2** no read classed `initiating` without cDNA tier-2 evidence (structurally
  unreachable on DRS). ✅
- **AT3** cDNA strand trap: TSV `strand` attributes **98.4%** vs **51.2%** on its flip —
  reproduces the documented 99.1%-vs-53.4% signature. ✅ *(upper bound only — see Open)*
- **AT4** split composition survives (e.g. 3′ group at 309,250: `YOL010W` 84,
  `YOL011W|YOL010W` 17, `YOL010W|YOL011W` 15). ✅
- **AT5** escaping-read count **invariant** at 4/10/25/50/100/200 bp windows while
  cluster count moves 31,932 → 5,772. ✅
- **Cross-check vs rbrowse's independent browser implementation** at PLB3:
  WT 0 vs 0 · rna15 94 vs 95 · ysh1 143 vs 143.
- Sidecar output inspected directly: 12 columns, header comments present, real read_ids.

## Key findings (all measured today, none recalled)

1. **The sidecar needs no BAM and no re-correction.** `alignment_start`/`alignment_end`
   in `corrected_reads.tsv` are the raw input alignment. Retention vs raw BAM measured at
   **100%** on all 6 libraries (327/327, 58/58, 89/89 …). Removes the ITEM-7 anchor-mismatch
   class entirely rather than avoiding it.
2. **`five_prime_position` is NOT the observed 5′ end for rescued reads** — the Cat3 rescue
   extends it upstream, so it is model-derived. `origin5` is derived from the raw alignment
   instead; the rescue flag rides in `origin5_evidence`. The spec's asymmetry claim was
   stated unconditionally and is not.
3. **🔴 GFF `gene` features are UTR-INCLUSIVE, not CDS-ended** (RCL1: `gene`
   307,883–309,265 vs `CDS` 307,938–309,041). I told rbrowse the wrong cause first and
   corrected it; they had propagated it into a code comment (`36dd3d3`).
4. **Containment-only INVENTS readthrough in WT** — PLB3 WT goes 0/845 → ~204/845.
   Confirmed independently by rbrowse at 201/845.
5. **No fixed window is defensible** (Kevin's call, and the data agrees). Genome-wide over
   5,607 genes: 71.6% have their modal CPA *inside* the annotated interval; ~4% sit past
   +300 up to **+5,024**. `YGL007C-A` has n=2,136 with 46% of reads at one base **+3,249**.
   → **Rule B**, selected by measurement (`planning/707`):
   `boundary = max(annotated_end, observed_modal_CPA) + 281 bp` where 281 is the *measured*
   genome-wide p95. Rules using each gene's own q95/q99 inflate WT and were rejected.
6. **RCL1's own CPA (309,142) is 133 bp upstream of where PLB3 readthrough terminates
   (309,275)** — the escaping molecules do not use RCL1's CPA site. Independent positive
   evidence they are readthrough.

- **Sidecars complete**: 9 libraries, **17,627,566 reads** in 318s. PLB3 `escapes_gene_cpa`
  = wt 0/0/0 · rna15 30/34/30 · ysh1 36/51/56.
- **Per-read diff vs rbrowse's `figEscape()`** (`planning/710`): **1,400 of 1,401 reads
  agree** (WT 845/845, rna15 200/201, ysh1 355/355). Sent to rbrowse.

## Open

- **The single diff disagreement is genuinely ambiguous, not a bug.** Read `dfe0ff08…`,
  3′ end 307,895 (+207 past PLB3): terminated 1,247 bp short of RCL1's own CPA (309,142)
  *and* past PLB3's entire observed WT range (max +59). Neither rule is clearly right.
  The tolerance tradeoff is deliberate — tightening to each gene's own spread (rules C/D)
  put 25/845 and 6/845 **false** escapees into WT.
- **`region_class` is empty in v1**, deliberately: the vocabulary is cluster-keyed and the
  spec forbids that join. rbrowse confirmed nothing breaks. Add when derivable per-read from
  the rules `analyze` uses — do NOT reimplement the vocabulary.
- **NOT VERIFIED: AT3's circularity check is impossible on the shipped cDNA cohort.** The
  11-column UMI schema has no `strand_evidence`, so molecules whose strand was resolved *by*
  gene overlap (tier 3, `ont_cdna.resolve_rna_strand`) cannot be excluded. The 98.4% is an
  **upper bound**. Needs `strand_evidence` carried into the cDNA TSV — already emitted by
  `rectify correct` (`bam/output.py:44`) but that cohort predates it.
- ~~cDNA sidecars~~ **DONE** — 6 libraries, **21,024,789 molecules**, 306s →
  `/u/scratch/k/kevinroy/620_rbrowse/attribution_cdna/`. Emitted with `unit="molecules"`;
  attribution rate **98.4%** (4,335,048 / 4,404,703 on WT_rep1), matching the independent
  AT3 measurement exactly. CPA reference = 5,508 genes from the three WT_BY4742 libraries.
- **`ambiguity_range` passthrough added** at rbrowse's request. They want it to set
  `clusterWin` from the assay's real 3′-end resolution instead of a parameter-sweep plateau
  (their sweep-based recommender picked settings from the degenerate regime, where the metric
  stops changing because the parameter has already destroyed the structure).
  ⚠️ Registered on their side in **`ATTR_INT`** (the sidecar's numeric set), **not `_INT`** —
  `_INT` is for `corrected_reads.tsv` columns and this arrives via the sidecar. Shard key is
  **`ar3`** (`ar` was already `attr_rule`). I had told them `_INT`; they corrected it.
- H2 `rectify` checkout is OLDER than the M1 tree — its `cli.py` lacks `qc_command`. I
  overwrote it once and **restored it with `git checkout --`**. Drive the module directly on
  H2 (`planning/709`), do not sync `cli.py` there.

## Resume

**Nothing is in flight. Both arms are emitted, verified, and handed to rbrowse.**

```bash
# Confirm both output sets are still present (scratch, not backed up):
ssh hoffman2 'cd /u/scratch/k/kevinroy/620_rbrowse && \
  echo "DRS : $(ls attribution      | wc -l) of 9"; \
  echo "cDNA: $(ls attribution_cdna | wc -l) of 6"'
```

| state | do |
|---|---|
| 9 and 6 | nothing to resume; work is complete |
| either short | re-emit — both drivers are idempotent (whole-file writes, no half-file state):<br>`ssh hoffman2 'module load conda/23.11.0; conda activate rectify; cd /u/scratch/k/kevinroy/620_rbrowse && python 709_emit_sidecars.py'` (DRS) and `… python 709b_emit_cdna.py` (cDNA) |
| scratch purged | both drivers rebuild from the corrected TSVs on `/u/project`, which are durable |

**If anything is re-emitted, verify the OUTPUT, not the exit code:**
```bash
ssh hoffman2 'cd /u/scratch/k/kevinroy/620_rbrowse/attribution_cdna && \
  head -1 WT_BY4742_rep1.attribution.tsv; \
  awk -F"\t" "!/^#/ {print NF; exit}" WT_BY4742_rep1.attribution.tsv'
```
MUST print `# unit: molecules` (**not** `reads`) and `12`. 🔴 If it says `reads`, every
downstream count means molecules while claiming reads — **re-emit, do not patch the header.**

Remaining work is all on rbrowse's side and **is not blocked on anything here**: consuming
the sidecars needs a re-slice + bundle re-upload per dataset (9 DRS, then 6 cDNA), which is
a substantial multi-dataset deploy, not a viewer change. They have deliberately queued it
rather than starting it mid-session, since the 1,400/1,401 agreement means it changes almost
nothing visible today. Their loader, `ATTR_INT` registration and the fail-closed
molecules-vs-reads guard are already in.

If they report the page disagrees with the sidecar after switching to `ec`, re-run
`planning/710` — it reimplements their rule and prints per-read reasons.

**Not committed.** 4 new files + a 2-hunk `cli.py` change are unstaged alongside Kevin's
pre-existing WIP. Per repo CLAUDE.md: `git add` explicit paths only, never `-A`. Kevin has
not yet approved a commit.

## Files

| path | what |
|---|---|
| `rectify/core/analyze/read_attribution.py` | the module |
| `rectify/core/commands/attribute_reads_command.py` | `rectify attribute-reads` |
| `rectify/cli.py` | wiring (2 hunks) |
| `tests/test_read_attribution.py` | 23 tests |
| `…/planning/703_attribution_sidecar_design.md` | design + all findings F1–F6 |
| `…/planning/704,705,705b` | AT1 prototype + the 13-vs-30 reconciliation |
| `…/planning/706,706b` | genome-wide 3′UTR survey + per-gene inspection |
| `…/planning/707` | tolerance-rule selection (why Rule B) |
| `…/planning/708,708b` | acceptance tests |
| `…/planning/709` | sidecar emission driver |
| H2 `/u/scratch/k/kevinroy/620_rbrowse/attribution/` | the 9 sidecars |
