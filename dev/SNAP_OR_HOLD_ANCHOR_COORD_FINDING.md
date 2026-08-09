# FINDING — the "3 non-canonical anchors" don't reproduce: reads place all 3 at CANONICAL GT-AG

**Written:** 2026-07-03. **Context:** asked to run the snap-or-hold do-no-harm test (the real-data half of
the native-aligner make-or-break) on the 3 anchors SQSTM1/TMED9/SLC35A4 via
`scripts/benchmark/c1_realign_3junctions.py`. It returned n=0 (no reads harvested). Investigating that n=0
uncovered a coordinate-derivation error that undermines the "non-canonical anchor" premise itself.

## What the harness expects vs what the reads actually do

The harness hardcodes each anchor as a non-canonical placement with a canonical motif 1–4 bp away
(`JUNCS` dict, from the COMPASS `rederive_111` derivation). But **both the COMPASS short reads and the
SG-NEx long reads place all three junctions at a CANONICAL GT-AG site 175–1017 bp away** from the harness
coordinate — and the harness coordinate has ~0 read support in either platform.

| anchor | harness coord (donor,acc) | harness motif | TRUE coord (donor,acc) | offset | motif @ true | short-rd n | long-rd n |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SQSTM1  | (179824400, 179832205) | GT..GA nonc | (179825226, 179833031) | +826  | **GT..AG** | 8161 | 17010 |
| TMED9   | (177592500, 177593474) | CG..CA nonc | (177592675, 177593649) | +175  | **GT..AG** | 2351 |  4006 |
| SLC35A4 | (140564954, 140565547) | AG..CA nonc | (140565971, 140566564) | +1017 | **GT..AG** |  521 |   906 |

- short-rd = `A549_rep1.consensus.bam` (COMPASS Illumina consensus — the exact BAM `rederive_111` used).
- long-rd  = sum over `sgnex_a549/alignments/a549_chr5_trimmed.{minimap2,deSALT,uLTRA,GMAP}.bam`.
- All aligners agree on the true placement; the harness coord gets essentially 0 reads at the right intron length.

## The true junctions are ANNOTATED (known introns), not novel

Checked the true coordinates against gencode v44 (`gencode.v44.basic.chr5.gtf`, 13,547 chr5 junctions) with
the same `load_annotated_junctions`→`normalize_junction`→membership test `rederive_111` uses:

| anchor | true coord (normalized) | annotated_exact | harness coord annotated? |
| --- | --- | --- | --- |
| SQSTM1  | (179825225, 179833030) | **True** | False (nearest annot start 77 bp, end 826 bp away) |
| TMED9   | (177592673, 177593647) | **True** | False (nearest 102 / 175 bp) |
| SLC35A4 | (140565968, 140566561) | **True** | False (nearest 96 / 200 bp) |

So the reads' true junctions are **known, catalogued, canonical GT-AG introns** — not novel and not
non-canonical. The harness coordinates land on unannotated positions with ~0 read support. This upgrades the
verdict from "mis-coordinated" to "the anchors are spurious": there is no real novel non-canonical junction
here to hold or snap.

## It is NOT a reference-frame difference (ruled out)

The COMPASS short reads were aligned to `GRCh38_gencode_v44.fasta`; the sgnex long reads to the Sumner
`GRCh38_chr5.fa`. Both chr5 are length 181,538,259 and are **byte-identical at all 3 loci** — a 60 bp window
around each true donor is found at the *exact same position* in both references (offset +0). So the 175/826/
1017 bp offsets are NOT a reference shift; they are the harness/derivation coordinates being wrong in *both*
references. Motifs verified independently by `samtools faidx` (1-based) AND pysam.

## Conclusion

- The **snap-or-hold do-no-harm test cannot run as specified**: the harness coordinates match neither
  platform's reads (n=0 harvest), and the reads' actual junctions are already **canonical GT-AG** — there is
  no nearby non-canonical placement to "snap" from.
- The **"3 non-canonical anchors" premise does not reproduce.** The COMPASS doc's open question
  (`dev/COMPASS_2corroborated_CROSSPLATFORM.md`: "genuinely non-canonical vs canonical-placed-a-few-bp-off")
  is answered by the data: both short + long reads place these at canonical GT-AG. The "non-canonical coord"
  was a mis-derived coordinate 175–1017 bp from the real (canonical) junction, not a 1–4 bp ambiguity.
- The true junctions are **annotated known introns** (above), so the anchors are spurious, not merely
  "mischaracterized novel" ones.
- Root cause is a coordinate-derivation error **upstream** of the harness — NOT the reference and NOT a
  short-vs-long discrepancy (both platforms agree on the canonical annotated site). Loose thread for whoever
  fixes it: the COMPASS doc claims short-read support at 177592500 via a `|shift|≤5` relaxed window, yet the
  short reads are 175 bp away at 177592675 — `rederive`'s ≤5 window could NOT have matched them. So the error
  sits in the coordinates that SEEDED rederive (the GMAP-only-novels TSV — consistent with GMAP's known
  non-canonical mis-placement, the very GMAP-drop-test context) or a transcription step into
  `c1_realign_3junctions.py` JUNCS. Not fully traced (would need `rederive_111.json` +
  `gmap_only_recurrent_novels_chr5.tsv`); not needed to establish the blocker.

## Why this matters (for the PI)

These 3 anchors are cited as the real-data do-no-harm substrate for the native junction re-aligner (handoff:
"the workhorses recover them" non-canonically). If they are in fact canonical GT-AG (as both platforms show),
then: (a) the do-no-harm snap-or-hold test needs new/valid substrate, and (b) this specific "real non-canonical
junctions in A549" evidence line for the native aligner should be re-examined. It does not by itself
invalidate the broader native-aligner case (which has other lines), but it removes/So-changes this one.

## Repro (all on Sherlock; data read-only)

- Scratch: `/scratch/users/kevinroy/c1_realdata_dB/` (this session): `c1_realign_3junctions.py` (n=0 run →
  `snap_hold_3junctions.json`), `motif_probe.py` (true placements + motifs), `ref_compare.py` (reference
  identity check).
- env: `source …/anaconda3/etc/profile.d/conda.sh; conda activate rectify; PYTHONPATH=…/nj_panel_code`.
- Data: sgnex BAMs `…/rectify_human_validation/sgnex_a549/alignments/`; short-read
  `/oak/…/compass_a549_out/A549_rep1.consensus.bam`; refs
  `/scratch/users/kevinroy/sumner_lab/references/GRCh38_chr5.fa` and
  `/scratch/users/kevinroy/compass_a549/COMPASS/genome_references/GRCh38_gencode_v44.fasta`.
