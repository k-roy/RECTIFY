# ISSUE: BAM output names are misleading (`*.rectified.bam` vs `corrected_consensus.bam`)

## Problem

The two consensus-stage BAMs in a `run-all` output directory have names that invert their actual
meaning, which repeatedly confuses users (surfaced in the 2026-06-25 Chanfreau workshop):

- **`*.rectified.bam`** is the file that has **not** been rectified — it is the *pre-correction*
  merged multi-aligner alignment set. "rectified" is the single most misleading token in the output.
- **`corrected_consensus.bam`** is the actually-rectified product, but "consensus" implies agreement
  among aligners, when the selection is in fact **winner-take-all per read**.

## Verified data model (do not trust the names — trust this)

```
per-aligner BAMs (minimap2, mapPacBio, gapmm2, uLTRA, deSALT, gmap)
        │
        ├──► [merge / vote]  ───────────────►  *.rectified.bam          (PRE-correction)
        │                                       • one record per DISTINCT placement
        │                                       • Xa:Z:<aligner,aligner,…> = SUPPORTING aligners (a list)
        │                                       • a read appears multiple times if aligners disagree on placement
        │
        └──► [per-aligner end-correction] ──►  per_aligner_corrected/<aligner>/...
                        │
                        └──► [winner per read] ─►  corrected_consensus.bam   (the rectified output)
                                                   • one record per read (the WINNER's corrected alignment)
                                                   • winner = lowest HP-aware edit distance on the CORRECTED CIGAR
                                                   • tags: cp:i (corrected pos), RN:i (read number); no Xa
```

Key consequences, often misunderstood:

1. **`corrected_consensus.bam` does NOT derive from `*.rectified.bam`.** They are parallel branches
   off the per-aligner alignments. The final output is built from the per-aligner **corrected**
   outputs (`corrected_consensus.py`: *"select winning aligner per read across per-aligner corrected
   outputs … the entire read is taken from the winner"*), not from the merged `*.rectified.bam`.
2. **`*.rectified.bam` is effectively a QC / inspection artifact** — it carries the `Xa` aligner-vote
   tags (useful for "which aligners agreed") but is not on the correction path.
3. **The `Xa` tag is a list of supporting aligners, not "the winner."** The per-read winner is only
   resolved at the `corrected_consensus` stage, on corrected CIGARs.

### Evidence (wt_drs_rep1, chrXIII)

| file | primary records | unique read names | per-record tag |
| --- | --- | --- | --- |
| `*.rectified.bam` | 4334 | 2798 | `Xa:Z:mapPacBio,minimap2,gmap` (vote list) |
| `corrected_consensus.bam` | 2790 | 2781 | `cp:i`, `RN:i` (no `Xa`) |

Source: `*.rectified.bam` written in `rectify/core/commands/align_command.py` (~L815, L879);
`corrected_consensus.bam` written by `write_corrected_consensus_bam` in
`rectify/core/consensus/corrected_consensus.py` (docstring is authoritative on the winner-take-all
selection over per-aligner *corrected* outputs).

## Proposed rename

**Recommended (full fix, with one release of back-compat):**

| current | proposed | rationale |
| --- | --- | --- |
| `<sample>….rectified.bam` | **`<sample>.multialigned.bam`** (or `<sample>.aligner_votes.bam`) | says what it is: merged multi-aligner alignments with per-placement aligner votes; pre-correction |
| `corrected_consensus.bam` | **`<sample>.rectified.bam`** | frees "rectified" for the file that has actually been rectified (matches the tool name) **and adds the sample-name prefix the other BAMs already carry** |

### Sample-name prefix on the final BAM (do this regardless of the rename)

The final consensus BAM is the **only** output BAM without a sample prefix: the per-aligner BAMs are
`<sample>.minimap2.bam` etc., and the merged file is `<sample>…_trimmed.rectified.bam`, but the final
product is the bare `corrected_consensus.bam`. The moment it leaves its per-sample directory it becomes
unidentifiable — and it routinely does:

- the shared **`processed/alignments/`** reuse store collects BAMs from many samples in one place;
- **IGV** track names collapse to "corrected_consensus.bam" for every sample (the workshop's
  `stage_igv.sh` already works around this by copying it to `<sample>.corrected_consensus.bam`);
- any flat `find … -name corrected_consensus.bam` collection is ambiguous.

**Proposal:** prefix the final BAM (and its `.bai`) with the sample id — `<sample>.rectified.bam`
(or, if the rename is deferred, at minimum `<sample>.corrected_consensus.bam`). This is independent
of the rename and worth doing on its own.

**Minimal fix (lower churn, fixes the worst confusion only):** rename just the misleading file
`*.rectified.bam → *.multialigned.bam`, and leave `corrected_consensus.bam` as-is (or → `*.corrected.bam`
to drop the "consensus = agreement" implication). This removes "rectified" from the pre-correction file
without re-pointing the token, so there is no transient period where `*.rectified.bam` means two things.

> Avoid the swap-in-place hazard: if we move "rectified" onto the final output in the same release that
> a different file still ships as `*.rectified.bam` in older runs, downstream globs (`*.rectified.bam`)
> silently pick up the wrong file. Either do the minimal fix, or ship back-compat symlinks for one release.

## Migration / blast radius

Names are referenced in (non-exhaustive — grep before changing):

- **Code:** `align_command.py`, `consensus/corrected_consensus.py`, `commands/split_command.py`
  (the chunked pipeline writes `corrected_consensus.bam`), `commands/run_command.py`
  (`--trust-existing-bams` help text mentions `rectified.bam` / `consensus.bam`),
  `commands/restore_polya_command.py`, `commands/batch_command.py` / `run/multi_sample.py`.
- **Tests:** any fixture asserting on these filenames.
- **Docs:** `docs/ARCHITECTURE.md` and the output-layout sections; the Chanfreau workshop guide
  (`/u/project/guillom/shared/workshop/RECTIFY_WORKSHOP.md`) which currently (correctly) tells users
  the `Xa` query targets `*.rectified.bam` and the analysis-ready BAM is `corrected_consensus.bam`.
- **Downstream / user scripts:** IGV staging (`stage_igv.sh` globs `*.minimap2.bam` and
  `corrected_consensus.bam`), DESeq/splice scripts that consume `corrected_consensus.bam`.

Suggested rollout: add the new names as the written outputs, keep symlinks (or a `--legacy-bam-names`
flag) from the old names for one minor release, update docs + `ARCHITECTURE.md` output map, then drop
the symlinks. Provenance sidecars already record the invocation, so renamed outputs stay traceable.

---
*Filed 2026-06-25 from the Chanfreau RECTIFY workshop. Analysis verified against the v0.9.0 shared
env output (`wt_drs_rep1` chrXIII) and the consensus source.*
