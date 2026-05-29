# cat5_chimeric — plotter handoff notes

**Status**: deferred. The per-segment provenance overlay needs data
that isn't currently in the validation bundle.

---

## What was supposed to land

Per the user's earlier ask (task #5):

> "For Cat5, I'd like the read overview to somehow highlight which
> portions of the final winning alignment came from which aligner."

Plan was:
- Color-band the overview track by segment provenance (per-segment
  aligner that contributed).
- In the End-correction panel (now removed in Phase 1 redesign — but
  any successor display), break Cat5 rows into per-segment rows tagged
  with the contributing aligner.

---

## Blocker — segment provenance not in the bundle

Checked the bundle's `rectified_pA_tail_trimmed.bam` for the 4 Cat5
reads. Tags actually present:

| qname | XA | XS | Xz | XV |
| --- | --- | --- | --- | --- |
| 040195ff (cat5_plus_1) | — | + (strand!) | — | cat5_plus_1 |
| 4d1e5c19 (cat5_plus_2) | — | + (strand!) | — | cat5_plus_2 |
| 8f86cb34 (cat5_minus_1) | — | — | — | cat5_minus_1 |
| 02165816 (cat5_minus_2) | YCR031C_id002 (gene_id!) | — | — | cat5_minus_2 |

Per the tag-rename refactor mentioned in earlier handoff
(`9e1b9a1 refactor(align): rename consensus-selection BAM tags to
X[lower]`), the chimeric provenance tags appear to have been moved
to lowercase variants or lost in regen. The expected schema (per
VALIDATION_READS.md) was:
- `Xz`: chimeric flag (1 = chimeric)
- `XA`: comma-separated aligner list per segment
- `XS`: segment count

What's actually written conflicts with that schema (XS contains
strand char; XA contains gene_id for one read; Xz is absent
entirely).

This is a **data-side / pipeline-side issue**, not a renderer issue.

---

## Path to unblocking

Two options for the debugger to consider:

1. **Restore the chimeric provenance tags** to the rectified bundle
   BAMs. The chimeric consensus path (`chimeric_consensus.py`) computes
   per-segment provenance internally and persists it during
   `merge_corrected_tsvs`. Wire those tags onto the output BAM records
   (e.g. `XA=mm2,deSALT` listing aligners per segment; `XS=<count>`;
   `Xz=1` for chimeric reads).

2. **Side-channel TSV**: emit a per-read TSV in the bundle (e.g.
   `rectified/cat5_segment_provenance.tsv`) with columns `read_id`,
   `segment_idx`, `ref_start`, `ref_end`, `winning_aligner`. The plotter
   reads this file and overlays per-segment bands on the overview track.

Either works. Option 2 is decoupled from the existing BAM tag
namespace and easier to extend (could include per-segment HP-ED, etc.).

---

## Rendering plan once data is available

In `render_overview`, after the per-aligner-group rows are drawn,
overlay a thin strip BELOW the winner's row showing per-segment color
bands. One color per contributing aligner (use the same palette as the
existing aligner-grouping legend).

For the per-aligner panel rows: when a Cat5 read is detected (XS > 1
or segments tsv has > 1 row), expand the winner group's row into N
sub-rows, each tagged `winner seg N (<aligner>)`. The painted CIGAR of
each sub-row is restricted to that segment's ref range.

Both treatments would only fire when segment data is present;
non-Cat5 reads render unchanged.

---

## Tasks queued

- Task #5 (plotter): blocked on data restoration.
- Suggested new debugger task: "Restore Cat5 chimeric provenance to
  rectified bundle (BAM tags or sidecar TSV)" — file under
  `validation_read_review/cat5_chimeric_findings.md` once a debugger
  picks it up.
