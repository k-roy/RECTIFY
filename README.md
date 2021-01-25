# RECTIFY

processed data can be found here: https://drive.google.com/drive/folders/1UETZGwg6iBHzm2H77nL9xPg4p_lkUQqk?usp=sharing

processed annotations can be found here: https://drive.google.com/drive/folders/1UHcrso8nN2RZ6hXgJ221ROnczVPCDdJG?usp=sharing
*annotations are integrated from multiple sources with priority according to  so that each coordinate in the genome is assigned to an annotation
*conflicts are re

Data processing steps:

1) BBMap_STAR_BWA.sh
2) collapse_quantseq_bam_to_three_prime_end.py
3) annotate_pA_sites_and_check_for_downstream_AG_richness.py
