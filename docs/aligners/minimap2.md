# Minimap2: junction annotation + the `-G` / `--max-intron` human pitfall

The `rectify align` command generates a junction BED from the GFF
annotation and passes it to minimap2 via `--junc-bed`. This improves
splice junction accuracy for long-read RNA-seq.

The junction BED is cached as `annotation.junc.bed` in the sample output
directory. Exact minimap2 command:

```
minimap2 -ax splice -uf -k14 -G 5000 --splice-flank=no --secondary=no --MD
    -t <threads>
    --junc-bed <sample_dir>/annotation.junc.bed
    --junc-bonus 9
    <genome.fsa.gz> <reads.fastq.gz>
```

Key flags:
- `-uf`: forward-strand-only (correct for direct RNA / cDNA sense reads).
- `-k14`: smaller k-mer for sensitivity on noisy nanopore reads.
- `-G 5000`: max intron (chainable gap) size. **This is the yeast default and
  the #1 misconfiguration on mammalian data — see below.**
- `--splice-flank=no`: disables GT-AG bonus (important for 3' end
  accuracy).
- `--MD`: required for indel artifact correction downstream.

---

## ⚠️ `-G` comes from `--max-intron` (default 5000 = yeast) and `--organism` does NOT raise it

`-G` is minimap2's maximum intron length: minimap2 will not chain a spliced
alignment across a gap larger than `-G`. RECTIFY sets `-G` from the
`rectify align --max-intron` flag, which **defaults to 5000**
(`align_command.py`, help text: *"appropriate for S. cerevisiae. Use 500000 or
larger for human data."*). Critically, **`--organism homo_sapiens` selects the
bundled genome/annotation but does NOT raise `--max-intron`** — so unless you
pass it explicitly, minimap2 runs with `-G 5000` on human data.

**Consequence on mammalian RNA:** human introns are routinely 10 kb–1 Mb. With
`-G 5000`, a read spanning a larger intron cannot be chained as one spliced
alignment — minimap2 soft-clips the distal exon or maps a fragment. minimap2 is
the consensus workhorse (~49 % of per-read wins on A549 chr5), so a wrong `-G`
silently degrades the *primary* alignment for large-intron transcripts.

**How bad, measured (A549 chr5 ONT DRS, 2026-05-25, panel run WITHOUT
`--max-intron` → `-G 5000`):** the impact on minimap2 is real but **bounded**,
not catastrophic. The consensus combo `deSALT+minimap2+uLTRA` still agreed on
312,796 reads, i.e. minimap2 tracked the splice-aware aligners on most reads —
the chr5 read set is dominated by transcripts whose introns are all <5 kb.
Damage is concentrated in the minority of reads spanning a >5 kb intron. (Contrast
mapPacBio, which was badly broken by its own intron misconfiguration — ~1k introns
total, a disjoint read set — see [mapPacBio.md](mapPacBio.md).)

**Caveat for junction analysis:** at large-intron loci a `-G 5000` minimap2
clipping at a 5 kb boundary while uLTRA/deSALT place the real 50 kb intron looks
like a genuine algorithmic *disagreement* (Cat6/Cat9-style) when it is only a
parameter artifact. Do not trust minimap2-vs-{uLTRA,deSALT} junction
disagreements at large-intron genes on human data until `-G` is set correctly.

**Fix — always pass `--max-intron 500000` (or larger) for human/mammalian:**

```
rectify align <reads> --genome <GRCh38> --organism homo_sapiens \
    --max-intron 500000  ...
```

**`--max-intron` is overloaded across the panel** — know what it does to each:
- **minimap2** `-G` ← this is the one it controls correctly; set it to 500000 for human.
- **gapmm2** `-i` (edlib terminal-rescue window) — no wall-time effect, larger is
  free; see [gapmm2.md](gapmm2.md).
- **mapPacBio** `intronlen` — historically (and incorrectly) wired here; mapPacBio
  needs `intronlen=10` + `maxindel=200000` instead; see [mapPacBio.md](mapPacBio.md).

---

## Quick check

```bash
rectify align <reads> --genome <GRCh38> --organism homo_sapiens \
    --max-intron 500000 --aligners minimap2 -o /tmp/mm2_smoke -t 8
# spliced reads should show N ops well over 5 kb on a mammalian sample:
samtools view /tmp/mm2_smoke/*.minimap2.bam \
  | awk '$6 ~ /[0-9]{4,}N/' | head
```

If the largest `N` op tops out near 5,000 on human data, `-G` is still at the
yeast default — `--max-intron` was not raised.
