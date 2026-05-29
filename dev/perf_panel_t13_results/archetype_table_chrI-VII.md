# Task #13 expensive-read panel categorization
panel-dir: dev/perf_panel_t13_results

## Whole-set + panel sizing (per aligner)

| Aligner | reads scanned | stride | total correct time (s) | threshold (ms) | panel reads | panel time (s) | panel share of total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| deSALT | ? | ? | ? | ? | 200 | 883.7 | ? |
| gapmm2 | ? | ? | ? | ? | 200 | 1883.8 | ? |
| mapPacBio | ? | ? | ? | ? | 200 | 1039.3 | ? |
| minimap2 | ? | ? | ? | ? | 200 | 1705.9 | ? |
| uLTRA | ? | ? | ? | ? | 200 | 1175.3 | ? |
| **POOLED** | 0 | - | 0.0 | - | 1,000 | 6687.9 | 0.0% |

## Archetype breakdown per aligner (count | % of panel time)

| Archetype | deSALT | gapmm2 | mapPacBio | minimap2 | uLTRA | POOLED |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `excluded_BUG` | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) |
| `bigN_artifact_intron` | 15 (4.0%) | 7 (1.8%) | 0 (0.0%) | 0 (0.0%) | 3 (0.1%) | 25 (1.1%) |
| `large_5prime_clip` | 24 (24.8%) | 19 (35.2%) | 0 (0.0%) | 34 (58.8%) | 29 (43.7%) | 106 (35.9%) |
| `over_cap_candidates` | 2 (0.3%) | 9 (17.4%) | 12 (2.0%) | 4 (0.4%) | 1 (0.2%) | 28 (5.4%) |
| `high_cigar_op` | 20 (20.6%) | 38 (22.4%) | 74 (70.6%) | 19 (14.1%) | 46 (28.0%) | 197 (28.5%) |
| `dead_end_with_work` | 17 (1.3%) | 10 (0.4%) | 14 (1.0%) | 12 (0.5%) | 13 (0.9%) | 66 (0.8%) |
| `productive_candidate` | 122 (49.0%) | 116 (22.7%) | 99 (26.0%) | 131 (26.2%) | 108 (27.1%) | 576 (28.4%) |
| `other` | 0 (0.0%) | 1 (0.0%) | 1 (0.4%) | 0 (0.0%) | 0 (0.0%) | 2 (0.1%) |

## POOLED archetype time-share (the headline)

| Archetype | panel reads | panel time (s) | % of panel time |
| --- | ---: | ---: | ---: |
| `excluded_BUG` | 0 | 0.0 | 0.0% |
| `bigN_artifact_intron` | 25 | 71.1 | 1.1% |
| `large_5prime_clip` | 106 | 2399.1 | 35.9% |
| `over_cap_candidates` | 28 | 360.1 | 5.4% |
| `high_cigar_op` | 197 | 1906.3 | 28.5% |
| `dead_end_with_work` | 66 | 50.3 | 0.8% |
| `productive_candidate` | 576 | 1896.5 | 28.4% |
| `other` | 2 | 4.5 | 0.1% |
| **TOTAL** | 1,000 | 6687.9 | 100.0% |

## Pooled feature medians per archetype

| Archetype | n | med_elapsed_ms | med_five_clip | med_n_pool_nearby | med_n_cigar_ops | med_n_N | med_n_op_intervals | med_mapq |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `bigN_artifact_intron` | 25 | 1045.1 | 0.0 | 8.0 | 27.0 | 99791.0 | 1.0 | 47.0 |
| `large_5prime_clip` | 106 | 4635.7 | 244.0 | 5.0 | 18.0 | 0.0 | 0.0 | 60.0 |
| `over_cap_candidates` | 28 | 1303.8 | 0.0 | 33.0 | 32.0 | 435.0 | 1.0 | 60.0 |
| `high_cigar_op` | 197 | 1901.0 | 0.0 | 9.0 | 86.0 | 0.0 | 0.0 | 49.0 |
| `dead_end_with_work` | 66 | 696.4 | 0.0 | 8.0 | 20.0 | 0.0 | 0.0 | 60.0 |
| `productive_candidate` | 576 | 983.0 | 0.0 | 7.0 | 22.0 | 0.0 | 0.0 | 60.0 |
| `other` | 2 | 3695.2 | 12.0 | 0.0 | 43.0 | 272.0 | 1.0 | 31.0 |

