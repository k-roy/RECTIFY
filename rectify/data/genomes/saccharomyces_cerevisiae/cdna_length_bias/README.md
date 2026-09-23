# S. cerevisiae ONT PCR-cDNA length-bias calibration

Bundled calibration for `rectify cdna-length-correct` with `--Scer`. Machine-readable provenance, file checksums and
the certified scope are in `PROVENANCE.json`.

| File | Rows | What it is |
| --- | --- | --- |
| `gene_lengths.tsv` | 5,899 | CDS-defined genes (stop codon coordinates, spliced CDS length) with spliced transcript lengths |
| `panel.tsv` | 405 | Genes whose direct RNA abundance is stable across 11 anchor-away genotype sets (SD ≤ 0.5 log2), with ≥ 100 molecules per library in both assays |
| `drs_reference_cpm.tsv` | 5,875 | Direct RNA cohort-mean CPM (33 RNA004 libraries), counted with the same rule as the cDNA |

## Where it is certified

The panel's stability was measured only in **W303 anchor-away strains grown in YPD and depleted with rapamycin**.
Declare each library's `background`, `medium` and `condition` in `fit --library-meta`; a library is certified only
when they are `W303-AA`, `YPD` and `rapamycin`. Every other strain (for example BY4742), medium (for example SD),
condition (for example a stress) and every undeclared library is reported as **not certified**. The correction
still runs, but it is an extrapolation.

## Gene universe

ORFs with a CDS in the bundled GFF and a transcript model in the ONT cDNA atlas used for transcript lengths, minus
mitochondrial genes, genes overlapping the rDNA locus (chrXII:450,000–470,000), and eight genes whose loci collect
reporter-plasmid transcripts in the calibration cohort: PSP2 (YML017W), MDJ1 (YFL016C), ADH1 (YOL086C), PDC1
(YLR044C), URA3 (YEL021W), CYC1 (YJR048W), GAL1 (YBR020W) and GAL10 (YBR019C). Those eight are therefore absent
from counts made with this gene table. To count them, build your own table with `cdna-length-correct genes` and
build a matching panel and reference.

## Rebuilding

```bash
rectify cdna-length-correct genes --gff GENOME.gff3 --transcript-models MODELS.tsv --require-transcript-model \
    --exclude-genes YML017W,YFL016C,YOL086C,YLR044C,YEL021W,YJR048W,YBR020W,YBR019C \
    --exclude-region chrXII:450000-470000 -o gene_lengths.tsv
rectify cdna-length-correct panel --cdna-counts CDNA_COUNTS.tsv --reference-counts DRS_COUNTS.tsv \
    --reference-groups DRS_GROUPS.tsv --gene-table gene_lengths.tsv -o panel_out/
```

`panel_out/panel.tsv` lists every considered gene with `panel` and `strict` flags; the bundled `panel.tsv` keeps
the 405 panel genes. `panel_out/reference_cpm.tsv` is the reference.
