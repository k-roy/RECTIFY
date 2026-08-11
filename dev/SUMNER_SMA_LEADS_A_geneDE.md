# Sumner SMA leads — SWEEP A: gene abundance DE (SMA vs control, RIN-matched good subset)
Source: apa_3prime/RESULTS_gene_de_orient_20260613.tsv (June 2026; log2FC = mean log2CPM SMA - control).
CAVEAT: shallow — only 20 genes with >=30 reads (DRS 3'-biased, small n). SUGGESTIVE candidates, not powered hits.

## Positive control: SMN1 log2FC -4.18 (DOWN) — SMN1-null validated.

## Top DOWN in SMA (biology-flagged)
| gene | log2FC | reads | note |
|---|---|---|---|
| GTF2H2 | -2.60 | 81 | ★ 5q13 SMA LOCUS (adjacent SMN1/2; co-affected by 5q deletions) — strong anchor |
| FST (follistatin) | -2.31 | 147 | ★ myostatin antagonist, muscle-mass relevant in a MN/muscle disease |
| NMUR2 | -3.14 | 50 | neuromedin U receptor (neuronal) |
| ADAMTS6 | -2.57 | 41 | ECM protease |
| FBN2 | -2.57 | 71 | fibrillin-2 (ECM/neuromuscular) |
| ANXA2R-OT1 | -3.36 | 113 | annexin A2 receptor overlapping transcript |
| GUSBP3 | -2.84 | 47 | pseudogene |

## Top UP in SMA
| gene | log2FC | reads | note |
|---|---|---|---|
| TARS1-DT | +4.36 | 124 | ★ threonyl-tRNA synthetase divergent transcript (aaRS -> motor neuropathy link) |
| MIR4458HG | +4.20 | 180 | miRNA host gene |
| PCDHGA6 / PCDHGA3 / PCDHB14 | +2.97/+2.63/+1.80 | 139/144/53 | ★ protocadherins (neuronal adhesion, MN-relevant) |
| CEP72 | +3.12 | 98 | centrosomal |
| SLCO4C1 | +2.22 | 32 | solute carrier |

VALIDATION PATH: qPCR / targeted DRS depth on GTF2H2, FST, TARS1-DT, the protocadherin cluster.
