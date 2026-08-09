# ncRNA atlases — pluggable CUT/SUT/XUT layer for the CPA classifier

This directory holds the **supplementary ncRNA annotation layer** (planning/169)
that the transcript-model CPA classifier (`rectify/core/analyze/transcript_model.py`)
composes on top of the SGD core GFF. The core GFF is never edited; this layer is a
separate, versioned, swappable set of study-defined tracks (CUTs, SUTs, XUTs, …).

## Files
- `atlases.yaml` — the registry. Maps a named atlas → its track files + per-track
  provenance (`source`, `ncrna_class`, `genome_build`, PMID, original build, liftover).
- `demo_synthetic.gff` — a tiny synthetic fixture (2 CUTs) used by the unit tests and
  `--ncrna-atlas demo_synthetic`. **Not real biology.**
- `xu2009_cuts.gff`, `xu2009_suts.gff`, `vandijk2011_xuts.gff` — the **production**
  `chanfreau_ncrna_v1` tracks. **Not committed** (they live on Sherlock Oak). Populate
  once before using `--ncrna-atlas chanfreau_ncrna_v1`.

## Usage
```
# named registry atlas
rectify analyze ... --annotation core.gff.gz --ncrna-atlas chanfreau_ncrna_v1
# ad-hoc, force-tag a single-class GFF/BED
rectify analyze ... --ncrna-annotations mycuts.gff:MyStudy2024:CUT
# default (no flag): SGD-core-only (snoRNA/tRNA/snRNA/rRNA come free)
```

## Populating the production tracks (once)
From Sherlock Oak (`/oak/stanford/groups/larsms/Users/kevinroy/common/annotation_files/non_coding_rna/`):
- `CUTs_Xu2009.gff`  → `xu2009_cuts.gff`  (925 CUT, Xu 2009, SGD-lifted to R64-1-1)
- `SUTs_Xu2009.gff`  → `xu2009_suts.gff`  (847 SUT; **strip the stray leading `CUTs` line**, 170b §2c)
- `XUTs_VanDijk2011.gff` → `vandijk2011_xuts.gff` (1658 XUT, van Dijk 2011, sacCer3)

All three are R64/sacCer3 nuclear coordinates — coordinate-compatible with the bundled
R64-5-1 genome, **no runtime liftover needed** (planning/170b §3). Each `.gff`'s own
`source` column already carries the study tag; the classifier force-tags `ncrna_class`
from the registry so a single-class file needs no special col-3 typing.

## Build compatibility
The loader asserts an atlas is compatible with the working genome on the shared **R64 /
sacCer3 nuclear** frame — tolerant of the SGD `#!assembly R64-4-1` label lag (a hard
exact `-5-1` match would false-fail; planning/170b §5e). A genuine build change (e.g. a
future non-R64 assembly) errors loudly → re-lift the atlas + bump the registry.

Provenance: planning/169 (design), 170b (atlas recon), 170d (oracle).
